!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PIC_ModParticles
  use PC_ModSize,ONLY: nPType, nElectronMax, x_, y_, z_
  use PC_ModSize,ONLY: nX, nY, nZ, nCell_D, nDim, MaxDim
  use PC_ModSize,ONLY: nHybrid
  use PIC_ModMain,ONLY: c, c2, Dt, Dx_D, DxInv_D, CellVolume, &
       SpeedOfLight_D, vInv, uTh_P
  use PC_ModParticleInField,ONLY: &
       State_VGBI,add_current, add_moments, add_density
  use PC_ModMpi, ONLY: &
       get_min_val_rho, get_max_val_rho, get_rho_avr, get_rho_int
  use PIC_ModFormFactor,ONLY: HighFF_ID, HighFFNew_ID,&
       Node_D, NodeNew_D, get_form_factors
  use PIC_ModProc,      ONLY:iProc,iError
  use ModNumConst,      ONLY: cHalf, cTwoPi
  use PC_BATL_particles
  use PC_BATL_lib, ONLY: CoordMin_DB, CoordMax_DB, &
       CoordMin_D, CoordMax_D
  implicit none

  integer,parameter :: W_ = nDim
  !Index for velocity(momentum)
  integer,parameter :: Wx_ = W_+x_, Wy_ = W_+y_, Wz_ = W_+z_

  integer,parameter::Electrons_ = 1 - nHybrid, Electron_=Electrons_

  !Structures
  real,dimension(Electrons_:max(nPType,Electrons_)) :: M_P, Q_P
  !Particle's mass and charge

  real,dimension(1:max(nPType,1))    :: Energy_P
  real,dimension(1:max(nPType,1))    :: OmegaPDtMax_P

  !Only at the root PE:
  real            :: nTotal_P(1:max(nPType,1)) = 0
  !\
  ! Add velocity
  !/
  logical :: DoAddVelocity_P(1:max(nPType,1)) = .false.
  real    :: VelocityToAdd_DP(Wx_:Wz_,1:max(nPType,1)) = 0.0

  !Methods
  public::set_particle_param       !Assigns M_P, Q_P and allocates coordinate
  ! arrays
  public::put_particle             !Add particle with known coordinates
  public::advance_particles   !Advance particles and collect moments info in 
                              !certain timestep
contains
  !======================================
  subroutine set_particle_param
    logical :: DoInit=.true.
    integer :: iSort
    !--------------------------
    if(nPType<1)RETURN
    if(DoInit)then
       DoInit=.false.
    else
       if(iProc==0)write(*,*)&
            'Particle arrays are deallocated, information is lost'
       do iSort=Electrons_,nPType
          deallocate(Particle_I(iSort)%State_VI,Particle_I(iSort)%iIndex_II)
       end do
    end if

    do iSort= 1, nPType
       Particle_I(iSort)%nVar      = Wz_
       Particle_I(iSort)%nIndex    = 0
       Particle_I(iSort)%nParticle = 0
       Particle_I(iSort)%nParticleMax=nElectronMax
    end do
    call allocate_particles
  end subroutine set_particle_param
  !================================
  subroutine put_particle(iSort,Coords_D,iBlock,W_D)
    integer,intent(in)::iSort, iBlock
    real,intent(in)::Coords_D(nDim)
    real,intent(in),optional::W_D(:)
    integer:: nParticle
    !---------------------------------------------
    nParticle = Particle_I(iSort)%nParticle + 1
    Particle_I(iSort)%nParticle = nParticle
    if( nParticle > Particle_I(iSort)%nParticleMax )then
       write(*,*)&
            'Cannot put particle of sort ', iSort,' at PE=',iProc
       call CON_stop('Simulation stops')
    end if
    Particle_I(iSort)%State_VI(x_:nDim,nParticle) = Coords_D
    Particle_I(iSort)%iIndex_II(0,     nParticle)      = iBlock
    if(present(W_D))then
       Particle_I(iSort)%State_VI(nDim+1:nDim+3,nParticle) = W_D
    else
       Particle_I(iSort)%State_VI(nDim+1:nDim+3,nParticle) = 0.0
    end if
  end subroutine put_particle
  !==========================
  !== This routine allows to get reproducible random distribution independent 
  !== to the number of processors involved. This is extremely expensive way,
  !== as long as the operation of the random number calculation actually 
  !== becomes serial 
  !===========================
  subroutine parallel_init_rand(nCall,iSeed,iIndex)
    use ModMpi
    use PIC_ModProc
    use PIC_ModRandom
    integer, intent(in)           :: nCall, iSeed
    integer, intent(in), optional :: iIndex

    integer :: iSend_I(1), iRecv_IP(1,0:nProc-1), iLoop
    real :: Aux
    !-----------------------------------------
    call init_rand(iSeed,iIndex)
    if(nProc==1)return
    iSend_I(1) = nCall
    iRecv_IP = 0
    call MPI_ALLGATHER(iSend_I, 1, MPI_INTEGER, &
         iRecv_IP, 1, MPI_INTEGER, iComm, iError)
    if(iProc==0)RETURN
    do iLoop =1,sum(iRecv_IP(1,0:iProc-1))
       Aux = RAND(iIndex)
    end do
  end subroutine parallel_init_rand
  !==========================================
  !=========reading command #UNIFORM=========
  !==========================================
  subroutine read_uniform
    use ModReadParam,ONLY: read_var
    use PIC_ModMain, ONLY: nPPerCellUniform_P 
    integer:: iSort
    character(len=12) :: NameVar
    !--------------
    do iSort = 1,nPType
       write(NameVar,'(a10,i1,a1)') 'nPPerCell(',iSort,')'
       call read_var(NameVar,nPPerCellUniform_P(iSort))
       if(nPPerCellUniform_P(iSort)<0.and.iProc==0)&
            write(*,*)'Particles will neutralize electrons'
    end do
  end subroutine read_uniform
  !================================
  subroutine uniform_shared_field
    use PIC_ModProc
    use PC_ModMpi
    use ModMpi
    use PIC_ModRandom
    use PIC_ModMain, ONLY: nPPerCellUniform_P 
    use PC_BATL_tree, ONLY: nRoot_D
    use PC_BATL_grid, ONLY: find_grid_block

    integer        :: nPPerCell_P(1:max(nPType,1))
    real :: n_P(1:max(nPType,1))
    integer        :: nPPerPE, nResidual, nPTotal
    real           :: Coord_D(MaxDim) = 0.0, W_D(MaxDim)= 0.0
    logical        :: UseQuasiNeutral
    integer        :: iBlockOut=1, iProcOut=0, iSort, iDim, iP
    !--------------------------
    do iSort = 1, nPType
       State_VGBI(:,:,:,:,:,iSort) = 0
       if(nPPerCellUniform_P(iSort)==0)CYCLE

       if(nPPerCellUniform_P(iSort) > 0)then
          nPPerCell_P(iSort) = nPPerCellUniform_P(iSort)
          UseQuasiNeutral = .false.
       else
          nPPerCell_P(iSort) = -nPPerCellUniform_P(iSort)
          UseQuasiNeutral = .true.
       end if

       nPTotal = product(nCell_D(1:nDim)*nRoot_D(1:nDim)) * NPPerCell_P(iSort)
       nPPerPE = nPTotal/nProc; nResidual = nPTotal - nProc*nPPerPE

       if(iProc+1<=nResidual)nPPerPE = nPPerPE + 1

       if(UseQuasiNeutral)then
          !Initialize the same random number sequence as 
          !for electrons
          call parallel_init_rand(nDim * nPPerPE, Electrons_)
       else
          call parallel_init_rand(nDim * nPPerPE, iSort)
       end if

       do iP = 1, nPPerPE
          do iDim = 1,nDim
             Coord_D(iDim) = (CoordMax_D(iDim)-CoordMin_D(iDim))&
                  * RAND() + CoordMin_D(iDim)
          end do
          !\
          ! This verion only works for UseSharedBlock
          !/
          call find_grid_block(Coord_D, iProcOut, iBlockOut)
          call put_particle(iSort, Coord_D(1:nDim), iBlockOut)
          Coord_D(1:nDim) = (Coord_D(x_:nDim) - &
               CoordMin_DB(1:nDim,iBlockOut))*DxInv_D
          call get_form_factors(Coord_D(1:nDim),Node_D,HighFF_ID)
          call add_density(Node_D,HighFF_ID,iBlockOut,iSort)
          
       end do
       call pass_density(iSort)
       if(iProc==0)write(*,*)'Total number of particles of sort ', iSort,&
            ' equals ',sum(State_VGBI(1,1:nX,1:nY,1:nZ,1:nBlock,iSort))
       n_P(iSort) = Particle_I(iSort)%nParticle
    end do
    if(nProc==1)then
       nTotal_P = n_P
    else
       call MPI_reduce(n_P, nTotal_P, nPType, MPI_REAL,&
            MPI_SUM, 0, iComm, iError)
    end if

    if(iProc==0)then
       write(*,*)'Particles are distributed'
       do iSort = 1, nPType
          write(*,*)'Totally ',nTotal_P(iSort),' particles of sort ',iSort
       end do
    end if
    do iSort = 1, nPType
       if(uTh_P(iSort)<=0.0)CYCLE
       call parallel_init_rand(6*Particle_I(iSort)%nParticle,iSort,2)
       do iP = 1, Particle_I(iSort)%nParticle
          call thermalize_particle(iSort, W_D)
          Particle_I(iSort)%State_VI(Wx_:Wz_,iP) = &
               Particle_I(iSort)%State_VI(Wx_:Wz_,iP) + W_D
       end do
    end do
  end subroutine uniform_shared_field
  !============================================================
  !Distribute particles by blocks.
  subroutine uniform
    use PIC_ModMain,  ONLY: nTotBlocks, UseSharedField, &
         nPPerCellUniform_P 
    use PC_ModSize,  ONLY: nBlock
    use PIC_ModProc
    use PIC_ModRandom
    use PC_ModMpi
    use PC_BATL_lib,  ONLY: iTree_IA, Block_, iNode_B

    integer             :: nPPerCell_P(1:max(nPType,1)), nBPerPE, nBResidual,&
         nPPerBlock, i, iNodeStart, iNodeEnd,  iNode, iP, &
         iPStart, iPEnd, iBlock, iDim, iSort
    real                :: n_P(1:max(nPType,1)), Coord_D(nDim), W_D(MaxDim)
    logical             :: UseQuasiNeutral
    !--------------------------
    DoAddVelocity_P = .false. !Added here
    SORTS:do iSort = 1, nPType
       State_VGBI(:,:,:,:,:,iSort) = 0
       if(uTh_P(iSort)==0.0)W_D = 0.0
       if(nPPerCellUniform_P(iSort)==0)CYCLE SORTS

       if(nPPerCellUniform_P(iSort) > 0)then
          nPPerCell_P(iSort) = nPPerCellUniform_P(iSort)
          UseQuasiNeutral = .false.
       else
          nPPerCell_P(iSort) = -nPPerCellUniform_P(iSort)
          UseQuasiNeutral = .true.
       end if
       nPPerBlock = product(nCell_D(1:nDim)) * NPPerCell_P(iSort)
       if(UseSharedField)then
          !\
          ! Blocks are distributed just to distribute particles
          !/
          nBPerPE    = nTotBlocks/nProc 
          nBResidual = nTotBlocks - nProc*nBPerPE
          if(iProc+1<=nBResidual)then
             nBPerPE = nBPerPE + 1
             iNodeStart = nBPerPE*iProc + 1
             iNodeEnd   = nBPerPE*iProc + nBPerPE
          else
             iNodeStart = nBPerPE*iProc + nBResidual + 1
             iNodeEnd   = nBPerPE*iProc + nBResidual + nBPerPE
          end if
       else
          !\
          ! Blocks are already distrbuted
          !/
          iNodeStart =1; iNodeEnd = nBlock
       end if
       BLOCKS:do i = iNodeStart, iNodeEnd
          if(UseSharedField)then
             !\
             ! Blocks are distributed just to distribute particles
             !/
             iNode = i
             iBlock = iTree_IA(Block_,iNode)
          else
             !\
             ! Blocks are already distrbuted
             !/
             iBlock = i
             iNode  = iNode_B(iBlock)
          end if
          if(UseQuasiNeutral)then 
             !Initialize the same random number sequence as 
             !for electrons
             call init_rand(nTotBlocks*Electrons_ + iNode)
          else
             call init_rand(nTotBlocks*iSort      + iNode)
          end if
          !\
          ! Init another sequence of random numbers to
          ! assign velocity
          !/ 
          call init_rand(nTotBlocks*(iSort + nPType) + iNode,2)
          iPStart = Particle_I(iSort)%nParticle + 1
          iPEnd   = Particle_I(iSort)%nParticle + nPPerBlock
          do iP = iPStart, iPEnd
             do iDim = 1,nDim
                Coord_D(iDim) = (CoordMax_DB(iDim,iBlock) -&
                     CoordMin_DB(iDim,iBlock))&
                     * RAND() + CoordMin_DB(iDim,iBlock)
             end do
             if(uTh_P(iSort)>0.0)&
                  call thermalize_particle(iSort, W_D)
             call put_particle(iSort, Coord_D(1:nDim), iBlock, &
                  W_D + VelocityToAdd_DP(:,iSort))
             Coord_D(1:nDim) = (Coord_D(x_:nDim) - &
                  CoordMin_DB(1:nDim,iBlock))*DxInv_D
          end do
       end do BLOCKS
       n_P(iSort) = Particle_I(iSort)%nParticle
    end do SORTS
    if(nProc==1)then
       nTotal_P = n_P
    else
       call MPI_reduce(n_P, nTotal_P, nPType, MPI_REAL,&
            MPI_SUM, 0, iComm, iError)
    end if
    
    if(iProc==0)then
       write(*,*)'Particles are distributed'
       do iSort = 1, nPType
          write(*,*)'Totally ',nTotal_P(iSort),' particles of sort ',iSort
       end do
    end if
  end subroutine uniform

  !==========================================
  !=========reading command #FOIL============
  !==========================================
  subroutine read_foil
    use ModReadParam,ONLY: read_var
    use PIC_ModProc
    use ModConst,    ONLY: cDegToRad
    use PC_ModMpi,  ONLY: pass_density
    use ModMpi
    use PIC_ModMain, ONLY: nPPerCellFoil_P, FoilCenter_D, &
         FoilWidth_D, AngleFoil
    integer:: iSort, iDim

    !--------------
    do iSort = 1,nPType
       call read_var('nPPerCellFoil',nPPerCellFoil_P(iSort))
       if(nPPerCellFoil_P(iSort)<0.and.iProc==0)&
            write(*,*)'Particles will neutralize electrons'
    end do
    do iDim = 1, nDim
       call read_var('FoilCenter_D',FoilCenter_D(iDim))
    end do
    do iDim = 1, nDim
       call read_var('FoilWidth_D',FoilWidth_D(iDim))
    end do
    call read_var('AngleFoil',AngleFoil)
  end subroutine read_foil
  !================================
  subroutine foil
    use PIC_ModProc
    use PIC_ModRandom
    use ModConst,    ONLY: cDegToRad
    use PC_ModMpi,  ONLY: pass_density
    use ModMpi
    use PIC_ModMain, ONLY: nPPerCellFoil_P, FoilCenter_D, &
         FoilWidth_D, AngleFoil
    integer :: nPPerPE, nResidual, nPTotal, iSort, iDim, iP
    integer :: nPPerCell_P(1:max(nPType,1))
    integer :: n_P(1:max(nPType,1))
    real    :: Coord_D(nDim), W_D(MaxDim) 
    real    :: PrimaryCoord_D(nDim)
    real    :: angleSin, angleCos
    logical :: UseQuasineutral
    integer :: iBlockOut = 1, iProcOut 
    !--------------------------
    angleSin=sin(angleFoil*cDegToRad)
    angleCos=cos(angleFoil*cDegToRad)
    do iSort = 1, nPType

       if(nPPerCellFoil_P(iSort)==0)CYCLE

       if(nPPerCellFoil_P(iSort) > 0)then
          nPPerCell_P(iSort) = nPPerCellFoil_P(iSort)
          UseQuasiNeutral = .false.
       else
          nPPerCell_P(iSort) = -nPPerCellFoil_P(iSort)
          UseQuasiNeutral = .true.
       end if

       NPTotal = product(FoilWidth_D)/CellVolume &
            *nPPerCell_P(iSort)
       nPPerPE = nPTotal/nProc; nResidual = nPTotal - nProc*nPPerPE

       if(iProc+1<=nResidual)nPPerPE = nPPerPE + 1

       if(UseQuasiNeutral)then
          !Initialize the same random number sequence as 
          !for electrons
          call parallel_init_rand(nDim * nPPerPE, Electrons_)
       else
          call parallel_init_rand(nDim * nPPerPE, iSort)
       end if
       !
       PART:        do iP = 1, nPPerPE
          !
          do iDim = 1,nDim
             PrimaryCoord_D(iDim) = FoilWidth_D(iDim)*(RAND()-0.5)
          end do
          !
          Coord_D(1) = PrimaryCoord_D(2)*angleSin+PrimaryCoord_D(1)*angleCos
          Coord_D(2) = PrimaryCoord_D(2)*angleCos-PrimaryCoord_D(1)*angleSin
          if(nDim==3)Coord_D(nDim) = PrimaryCoord_D(nDim)
          !
          do iDim = 1,nDim
             Coord_D(iDim) = (FoilCenter_D(iDim)+Coord_D(iDim))/Dx_D(iDim)
          end do
!          if(  Coord_D(1) .lt. 0.0 .or. Coord_D(1) .ge. real(nCell_D(1))&
!               .or.&
!               Coord_D(2) .lt. 0.0 .or. Coord_D(2) .ge. real(nCell_D(2))&
!               .or.&
!               Coord_D(3) .lt. 0.0 .or. Coord_D(3) .ge. real(nCell_D(3))&
!               ) cycle PART
          call put_particle(iSort, Coord_D, iBlockOut)
          !call get_form_factors(Coord_D,Node_D,HighFF_ID)
          !hyzhou: I removed add_density. Need to modify this part later
          !call add_density(Node_D,HighFF_ID,1)
       end do PART
       n_P(iSort) = Particle_I(iSort)%nParticle
    end do
       !
    !Rho_GB = 0.0
    if(nProc==1)then
       nTotal_P = n_P
    else
       call MPI_reduce(n_P, nTotal_P, nPType, MPI_REAL,&
            MPI_SUM, 0, iComm, iError)
    end if
    !call pass_density(iSort)
    if(iProc==0)then
       write(*,*)'Particles are distributed'
       do iSort = 1, nPType
          write(*,*)'Totally ',nTotal_P(iSort),' particles of sort ',iSort
       end do
    !   write(*,*)'Particle density max:',max_val_rho()
    !   write(*,*)'Particle density min:',min_val_rho()
    !   write(*,*)'Particle density average:',rho_avr()
    end if
    do iSort = 1, nPType
       if(uTh_P(iSort)<=0.0)CYCLE
       call parallel_init_rand(6*Particle_I(iSort)%nParticle,iSort)
       do iP = 1, Particle_I(iSort)%nParticle
          call thermalize_particle(iSort, W_D)
          Particle_I(iSort)%State_VI(Wx_:Wz_,iP) = &
               Particle_I(iSort)%State_VI(Wx_:Wz_,iP) + W_D
       end do
    end do
  end subroutine foil
  !====================
  subroutine show_density(iSort)
    use PC_ModSize,  ONLY: nBlock
    use PIC_ModProc
    use PC_ModMpi
    integer, intent(in):: iSort
    real :: RhoInt, RhoMin, RhoMax, RhoAvr
    Aux_CB(:,:,:,1:nBlock) = &
         State_VGBI(1,1:nX,1:nY,1:nZ,1:nBlock,iSort)
    call get_min_val_rho(RhoMin)
    call get_max_val_rho(RhoMax)
    call get_rho_avr(    RhoAvr)
    call get_rho_int(    RhoInt)
    if(iProc==0)then
       write(*,*)'Total number of particles of sort ', iSort,&
            ' equals ',RhoInt
       write(*,*)'Density min=',Q_P(iSort)*RhoMin,', density max=',&
            Q_P(iSort)*RhoMax,', average density=',Q_P(iSort)*RhoAvr
    end if
  end subroutine show_density
  !==================
  !\
  !Generate Gaussian distribution using Box-Muller method
  !/
  subroutine thermalize_particle(iSort,W_D)
    use PIC_ModRandom
    integer, intent(in ) :: iSort
    real   , intent(out) :: W_D(MaxDim)
    real                 :: Energy, MomentumAvr
    integer              :: iW
    !-----------------------------------
    do iW = x_, z_
       Energy      = -uTh_P(iSort)**2 *LOG(RAND(2))
       MomentumAvr = sqrt(2.0*Energy)
       W_D(iW)     = MomentumAvr*COS(cTwoPi*RAND(2))
    end do
  end subroutine thermalize_particle
  !===========Reading command #THERMALIZE============
  subroutine read_temperature
    use ModReadParam,    ONLY: read_var
    integer:: iSort
    character(len=8) :: NameVar
    !--------------
    do iSort = 1,nPType
       write(NameVar,'(a6,i1,a1)') 'uTh_P(',iSort,')'
       !uTh_P is the thermal velocity in normalized units.
       call read_var(NameVar,uTh_P(iSort))
    end do
  end subroutine read_temperature
  !=======================
  subroutine pass_energy
    use ModMpi
    use PIC_ModProc

    real :: E_P(1:max(nPType,1))
    !---------------------------
    if(nProc==1)RETURN
    if(nPType<1)RETURN
    E_P = Energy_P
    call MPI_Reduce(&
         E_P, Energy_P, nPType, MPI_REAL,MPI_SUM,  0, iComm, iError)
  end subroutine pass_energy
  !=========================
  subroutine read_initial_velocity
    use ModReadParam, ONLY: read_var
    real      :: W_D(Wx_:Wz_)
    integer   :: iSort
    real,dimension(:,:),pointer :: Coord_VI
    integer   :: iP, nParticle
    !-----------------------------
    call read_var('iSort',iSort)
    DoAddVelocity_P(iSort) = .true.
    call read_var('Wx'   ,VelocityToAdd_DP(Wx_,iSort))
    call read_var('Wy'   ,VelocityToAdd_DP(Wy_,iSort))
    call read_var('Wz'   ,VelocityToAdd_DP(Wz_,iSort))
  end subroutine read_initial_velocity
  !===========================
  subroutine add_initial_velocity(iSort)
    use ModReadParam, ONLY: read_var
    integer,intent(in) ::iSort
    real,dimension(:,:),pointer :: Coord_VI
    integer   :: iP, nParticle
    !-----------------------------
    call set_pointer_to_particles(iSort,Coord_VI,nParticle=nParticle)

    do iP = 1,nParticle
       Coord_VI(Wx_:Wz_,iP) = Coord_VI(Wx_:Wz_,iP) + &
            VelocityToAdd_DP(Wx_:Wz_,iSort)
    end do
  end subroutine add_initial_velocity
  !==============================================================
  !u = u_0*sin(kx)
  subroutine add_velocity_sine
    use ModReadParam, ONLY: read_var
    real    :: Ampl,WaveNumber
    character(len=3) :: Direction
    integer :: iSort
    integer :: iP, nParticle
    real,dimension(:,:),pointer :: Coord_VI
    character(len=17) :: NameSub='add_velocity_sine'
    !---------------------------------------
    call read_var('iSort',iSort)
    call read_var('Amplitude' ,Ampl)      ! u_0
    call read_var('WaveNumber',WaveNumber)! k
    call read_var('Direction' ,Direction)
    call set_pointer_to_particles(iSort,Coord_VI,nParticle=nParticle)

    select case(Direction)
    case('x')
       do iP=1,nParticle
          Coord_VI(Wx_,iP) = Coord_VI(Wx_,iP) + &
               Ampl*sin(WaveNumber*Coord_VI(x_,iP))
       end do

    case('y')
       do iP=1,nParticle
          Coord_VI(Wy_,iP) = Coord_VI(Wy_,iP) + &
               Ampl*sin(WaveNumber*Coord_VI(y_,iP))
       end do

    case('z')
       do iP=1,nParticle
          Coord_VI(Wz_,iP) = Coord_VI(Wz_,iP) + &
               Ampl*sin(WaveNumber*Coord_VI(z_,iP))
       end do

    case default
       if(iProc==0) then
           write(*,*) NameSub // ' WARNING: unknown direction ' // &
                trim(Direction),' !!!'
           call CON_stop('Correct #ADDSINEWAVEVELOCITY!')
        end if

    end select

  end subroutine add_velocity_sine
  
  !=======================PARTICLE MOVER=========================!
  !Advance the particles in one timestep; calculate cell-centered
  !number density and velocity moments if DoComputeMoments==.true.
  subroutine advance_particles(iSort, &
       DoComputeMoments, DoPredictorOnly)
    use ModCoordTransform, ONLY: cross_product
    use PC_ModParticleInField,ONLY: &
         b_interpolated_d,e_interpolated_d
    use PC_ModHybrid,      ONLY: UseHybrid, add_predictor_current
    integer,intent(in) :: iSort
    logical,intent(in) :: DoComputeMoments
    !\
    ! If present, then only the predictor step is done.
    ! This option may be used: (1) to make a final plot or
    ! to print the last iteration data into the log file
    ! (which both use the velocity data to be advanced by a half 
    ! time step); (2) in some schemes for hybrid simulations (to 
    ! advance the particle current to the beginning of time step)
    !/
    logical, optional, intent(in) :: DoPredictorOnly
    real:: QDtPerM
    real,dimension(MaxDim)::QPerVDx_D
    real:: W2
    integer::iParticle, iBlock, nParticle, iDim
    real,dimension(nDim)::X_D
    real,dimension(x_:z_)   ::W_D,W12_D,EForce_D,BForce_D
    real    :: Gamma, Gamma12Inv
    real,pointer::Coord_VI(:,:)
    integer,pointer::iIndex_II(:,:)
    !----------------------------------------------------
    !Initialize the simulation for this sort of particles

    !1/2 * Q * dt / M, 
    !used to calculate the force effect
    QDtPerM = cHalf * Dt * Q_P(iSort) / M_P(iSort)

    !Q / V * (/\Delta x, \Delta y, \Delta z/)
    !Used to calculate J*dt in the charge-conserving scheme
    !QPerVDx_D = (Q_P(iSort)/CellVolume) * Dx_D

    QPerVDx_D(1:nDim) = Q_P(iSort)*vInv*Dx_D
    !\
    ! For directions orthogonal to grid 
    ! (should be multiplied by 
    ! velocity-to-the speed of light ratio  
    !/
    do iDim = nDim+1, MaxDim
       QPerVDx_D(iDim)= Q_P(iSort)*vInv*c*Dt
    end do

    call set_pointer_to_particles(&
         iSort,Coord_VI,iIndex_II,nParticle=nParticle)

    !Looping over particles
    do iParticle=1,nParticle
       !\
       ! Get block number
       !/
       iBlock = iIndex_II(0,iParticle)
       !Get coordinates and momentum
       X_D = (Coord_VI(x_:nDim,iParticle) - CoordMin_DB(1:nDim,iBlock))*&
            DxInv_D
       W_D = Coord_VI(Wx_:Wz_,iParticle)
       
       W2 = sum(W_D**2)
       Gamma = sqrt(c2 + W2)
       
       !Now W_D is the initial momentum, divided by mass, W2=W_D^2
       !Now Gamma is the initial Gamma-factor multiplied by c

       call get_form_factors(X_D,Node_D,HighFF_ID)
       
       !Electric field force, divided by particle mass 
       !and multiplied by \Delta t/2
       EForce_D = QDtPerM * e_interpolated_d(iBlock)
       !\
       !Add kinetic energy, divided my mc^2 and 
       !advanced by half time step 
       !/
       Energy_P(iSort) = Energy_P(iSort) + &
            M_P(iSort)/c*(W2/(Gamma+c) + sum(W_D*EForce_D)/Gamma)
       !Acceleration from the electric field, for the
       !first half of the time step:
       W_D = W_D + EForce_D
       Gamma12Inv = 1/sqrt( c2+sum(W_D**2) )
       !Get magnetic force divided by particle mass and by c 
       !and multiplied by \Delta t/2
       BForce_D = QDtPerM*Gamma12Inv*b_interpolated_d(iBlock)       

       !Add a half of the magnetic rotation:

       W12_D = W_D + cross_product(W_D,BForce_D)
       !\
       !Contribute to number density and velocity moments
       !/
       if(DoComputeMoments)then 
          call add_moments(W12_D*Gamma12Inv,&
               Node_D,HighFF_ID,iBlock,iSort)
       else
          call add_density(Node_D,HighFF_ID,iBlock,iSort)
       end if       
       if(UseHybrid)&
            call add_predictor_current(&
            QPerVDx_D,W12_D*Gamma12Inv,iBlock)
       if(present(DoPredictorOnly))CYCLE
       !Multiply the magnetic force by 2 to take a whole
       !rotation and reduce its magnitude not to perturb energy
       BForce_D = (2.0/(1.0 + sum(BForce_D**2))) * BForce_D

       !Get a final momentum

       W_D = W_D + cross_product(W12_D,BForce_D) + EForce_D

       Gamma = sqrt(c2+sum(W_D**2))
       !Now Gamma is the final Gamma-factor multiplied by c

       !Save momentum
       Coord_VI(Wx_:Wz_,iParticle) = W_D
       W_D = (1.0/Gamma)*W_D
       !Now W_D is the velocity divided by c

       !Update coordinate
       X_D = X_D + SpeedOfLight_D(1:nDim)*W_D(1:nDim)

       !New form factor
       call get_form_factors(X_D,NodeNew_D,HighFFNew_ID)
          
       !Contribute to the current
       call add_current(QPerVDx_D,W_D,iBlock)
       
       Coord_VI(1:nDim,iParticle)= &
            X_D*Dx_D(1:nDIm) + CoordMin_DB(1:nDim,iBlock)
    end do
    !\
    ! If only the predictor step is done, the particle 
    ! coordinates did not change
    !/
    if(.not.present(DoPredictorOnly))&
         call message_pass_particles(iSort)
  end subroutine advance_particles
end module PIC_ModParticles

