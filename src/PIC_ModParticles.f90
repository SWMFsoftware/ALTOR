!^CFG COPYRIGHT UofM
module PIC_ModParticles
  use PIC_ModSize,ONLY: nPType, nElectronMax,nDim, x_, y_, z_
  use PIC_ModSize,ONLY: nX, nY, nZ, nCell_D
  use PIC_ModMain,ONLY: c, c2, Dt, Dx_D, CellVolume, SpeedOfLight_D
  use PIC_ModParticleInField,ONLY: Rho_GB,add_current
  use PIC_ModParticleInField,ONLY: add_DensityVelocity
  use PIC_ModParticleInField,ONLY: b_interpolated_d,e_interpolated_d
  use PIC_ModParticleInField,ONLY: min_val_rho, max_val_rho, rho_avr
  use PIC_ModFormFactor,ONLY: HighFF_ID, HighFFNew_ID,&
       Node_D, NodeNew_D, get_form_factors
  use PIC_ModProc,      ONLY:iProc,iError
  use ModNumConst,      ONLY: cHalf
  use PC_BATL_particles, ONLY:allocate_particles, Particle_I

  implicit none

  integer,parameter :: W_ = nDim
  !Index for velocity(momentum)
  integer,parameter :: Wx_ = W_+x_, Wy_ = W_+y_, Wz_ = W_+z_

  integer,parameter::Electrons_=1,Electron_=Electrons_

  !Structures
  real,dimension(nPType) :: M_P, Q_P
  !Particle's mass and charge

  integer,dimension(nPType) :: n_P,nMax_P
  real,dimension(nPType)    :: Energy_P
  real,dimension(nPType)    :: OmegaPDtMax_P

  !Only at the root PE:
  integer,dimension(nPType) :: nTotal_P

  !Methods
  public::set_pointer_to_particles !Set pointer to the coordinate array of 
                                   !electrons or ions
  public::set_particle_param       !Assigns M_P, Q_P and allocates coordinate
                                   ! arrays
  public::put_particle             !Add particle with known coordinates
  public::advance_particles   !Advance particles and collect moments info in 
                              !certain timestep
contains
  !=============================
  subroutine set_pointer_to_particles(iSort,PointerToSet)
    integer,intent(in)::iSort
    real,dimension(:,:),pointer::PointerToSet
    nullify(PointerToSet)
    PointerToSet=>Particle_I(iSort)%State_VI
  end subroutine set_pointer_to_particles

  !=============================

  subroutine set_particle_param(MIn_P,QIn_P)
    real,dimension(nPType),intent(in),optional::MIn_P,QIn_P
    logical::DoInit=.true.
    integer::iSort
    if(present(MIn_P))M_P=MIn_P
    if(present(QIn_P))Q_P=QIn_P
    n_P=0
    if(DoInit)then
       DoInit=.false.
    else
       if(iProc==0)write(*,*)&
            'Particle arrays are deallocated, information is lost'
       do iSort=Electrons_,nPType
          deallocate(Particle_I(iSort)%State_VI,Particle_I(iSort)%iIndex_II)
       end do
    end if
    !Max number of particles is max number of electrons per the
    !charge ratio
    nMax_P= nint(nElectronMax * abs(Q_P(Electron_)/Q_P) )

    do iSort=Electrons_,nPType
       Particle_I(iSort)%nVar=Wz_
       Particle_I(iSort)%nIndex=0
       Particle_I(iSort)%nParticleMax=nMax_P(iSort)
    end do
    call allocate_particles
  end subroutine set_particle_param
  !===============================

  subroutine put_particle(iSort,PhaseCoords_D)
    integer,intent(in)::iSort
    real,dimension(nDim),intent(in)::PhaseCoords_D
    !--------------
    n_P(iSort) = n_P(iSort)+1
    if( n_P(iSort) > nMax_P(iSort) )then
       write(*,*)&
            'Cannot put particle of sort ', iSort,' at PE=',iProc
       call CON_stop('Simulation stops')
    end if
    Particle_I(iSort)%State_VI(x_:nDim,n_P(iSort)) = PhaseCoords_D
    Particle_I(iSort)%State_VI(nDim+1:nDim+3,n_P(iSort)) = 0.0

  end subroutine put_particle
  !===========================
  !== This routine allows to get reproducible random distribution independent 
  !== to the number of processors involved. This is extremely expensive way,
  !== as long as the operation of the random number calculation actually 
  !== becomes serial 
  !===========================
  subroutine parallel_init_rand(nCall,iSeed)
    use ModMpi
    use PIC_ModProc
    use PIC_ModRandom
    integer, intent(in)           :: nCall
    integer, intent(in), optional :: iSeed

    integer :: iSend_I(1), iRecv_IP(1,0:nProc-1), iLoop
    real :: Aux
    !-----------------------------------------
    call init_rand(iSeed)
    if(nProc==1)return
    iSend_I(1) = nCall
    iRecv_IP = 0
    call MPI_ALLGATHER(iSend_I, 1, MPI_INTEGER, &
         iRecv_IP, 1, MPI_INTEGER, iComm, iError)
    if(iProc==0)return
    do iLoop =1,sum(iRecv_IP(1,0:iProc-1))
       Aux = RAND()
    end do
  end subroutine parallel_init_rand
  !==========================================
  !=========reading command #UNIFORM=========
  !==========================================
  subroutine read_uniform
    use ModReadParam,ONLY: read_var
    use PIC_ModProc
    use PIC_ModMpi,  ONLY: pass_density
    use ModMpi

    integer:: nPPerCell_P(nPType)=0
    integer:: iSort
    character(len=12) :: NameVar
    !--------------
    do iSort = 1,nPType
       write(NameVar,'(a10,i1,a1)') 'nPPerCell(',iSort,')'
       call read_var(NameVar,nPPerCell_P(iSort))
       if(nPPerCell_P(iSort)<0.and.iProc==0)&
            write(*,*)'Particles will neutralize electrons'
    end do

    call uniform(nPPerCell_P)

    if(nProc==1)then
       nTotal_P = n_P
    else
       call MPI_reduce(n_P, nTotal_P, nPType, MPI_INTEGER,&
            MPI_SUM, 0, iComm, iError)
    end if

    if(iProc==0)then
       write(*,*)'Particles are distributed'
       do iSort = 1, nPType
          write(*,*)'Totally ',nTotal_P(iSort),' particles of sort ',iSort
       end do
    end if
  end subroutine read_uniform
  !================================
  subroutine uniform(nPPerCellIn_P)
    use PIC_ModProc
    use PIC_ModRandom
    integer, intent(in) :: nPPerCellIn_P(nPType)
    integer             :: nPPerCell_P(nPType)
    integer :: nPPerPE, nResidual, nPTotal, iSort, iDim, iP
    real    :: Coord_D(nDim)
    logical :: UseQuasiNeutral

    integer :: i,j,k
    !--------------------------
    do iSort = 1, nPType
       if(nPPerCellIn_P(iSort)==0)CYCLE

       if(nPPerCellIn_P(iSort) > 0)then
          nPPerCell_P(iSort) = nPPerCellIn_P(iSort)
          UseQuasiNeutral = .false.
       else
          nPPerCell_P(iSort) = -nPPerCellIn_P(iSort)
          UseQuasiNeutral = .true.
       end if

       nPTotal = product(nCell_D) * NPPerCell_P(iSort)
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
             Coord_D(iDim) = nCell_D(iDim) * RAND()
          end do
          call put_particle(iSort, Coord_D)
       end do
    end do
  end subroutine uniform

  !==========================================
  !=========reading command #FOIL============
  !==========================================
  subroutine read_foil
    use ModReadParam,ONLY: read_var
    use PIC_ModProc
    use ModConst,    ONLY: cDegToRad
    use PIC_ModMpi,  ONLY: pass_density
    use ModMpi

    integer:: nPPerCell_P(nPType)=0
    integer:: iSort
    real :: xFoilCenter(nDim)
    real :: xFoilWidth(nDim)
    real :: angleFoil
    !--------------
    do iSort = 1,nPType
       call read_var('nPPerCell',nPPerCell_P(iSort))
       if(nPPerCell_P(iSort)<0.and.iProc==0)&
            write(*,*)'Particles will neutralize electrons'
    end do
    call read_var('xFoilCenter',xFoilCenter(1))
    call read_var('yFoilCenter',xFoilCenter(2))
    call read_var('zFoilCenter',xFoilCenter(3))
    call read_var('xFoilWidth',xFoilWidth(1))
    call read_var('yFoilWidth',xFoilWidth(2))
    call read_var('zFoilWidth',xFoilWidth(3))
    call read_var('angleFoil',angleFoil)
    !
    Rho_GB = 0.0
    call foil(nPPerCell_P)
    if(nProc==1)then
       nTotal_P = n_P
    else
       call MPI_reduce(n_P, nTotal_P, nPType, MPI_INTEGER,&
            MPI_SUM, 0, iComm, iError)
    end if
    call pass_density(0)
    if(iProc==0)then
       write(*,*)'Particles are distributed'
       do iSort = 1, nPType
          write(*,*)'Totally ',nTotal_P(iSort),' particles of sort ',iSort
       end do
       write(*,*)'Particle density max:',max_val_rho()
       write(*,*)'Particle density min:',min_val_rho()
       write(*,*)'Particle density average:',rho_avr()
    end if
  contains
    !
    !================================
    subroutine foil(nPPerCellIn_P)
      use PIC_ModProc
      use PIC_ModRandom
      integer, intent(in) :: nPPerCellIn_P(nPType)
      integer             :: nPPerCell_P(nPType)
      integer :: nPPerPE, nResidual, nPTotal, iSort, iDim, iP
      real    :: Coord_D(nDim) 
      real    :: xPrime(nDim)
      real    :: angleSin, angleCos
      logical :: UseQuasineutral
      !--------------------------
      angleSin=sin(angleFoil*cDegToRad)
      angleCos=cos(angleFoil*cDegToRad)
      do iSort = 1, nPType

         if(nPPerCellIn_P(iSort)==0)CYCLE

         if(nPPerCellIn_P(iSort) > 0)then
            nPPerCell_P(iSort) = nPPerCellIn_P(iSort)
            UseQuasiNeutral = .false.
         else
            nPPerCell_P(iSort) = -nPPerCellIn_P(iSort)
            UseQuasiNeutral = .true.
         end if

         NPTotal = product(xFoilWidth)/CellVolume &
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
               xPrime(iDim) = xFoilWidth(iDim)*(RAND()-0.5)
            end do
            !
            Coord_D(1) = xPrime(2)*angleSin+xPrime(1)*angleCos
            Coord_D(2) = xPrime(2)*angleCos-xPrime(1)*angleSin
            Coord_D(3) = xPrime(3)
            !
            do iDim = 1,nDim
               Coord_D(iDim) = (xFoilCenter(iDim)+Coord_D(iDim))/Dx_D(iDim)
            end do
            if(  Coord_D(1) .lt. 0.0 .or. Coord_D(1) .ge. real(nCell_D(1))&
                 .or.&
                 Coord_D(2) .lt. 0.0 .or. Coord_D(2) .ge. real(nCell_D(2))&
                 .or.&
                 Coord_D(3) .lt. 0.0 .or. Coord_D(3) .ge. real(nCell_D(3))&
                 ) cycle PART
            call put_particle(iSort, Coord_D)
            call get_form_factors(Coord_D,Node_D,HighFF_ID)
            !hyzhou: I removed add_density. Need to modify this part later
            !call add_density(Node_D,HighFF_ID,1)
         end do PART
      end do
    end subroutine foil
    !
  end subroutine read_foil
  !=======================
  subroutine get_energy
    use ModMpi
    use PIC_ModProc

    real :: E_P(nPType), P2
    integer:: iSort, iP
    real,dimension(:,:),pointer::Coord_VI
    !-----------------------------------
    E_P = 0.0
    do iSort = 1, nPType
       call set_pointer_to_particles(iSort,Coord_VI)
       do iP = 1, n_P(iSort)
          P2 = sum(Coord_VI(Wx_:Wz_,iP)**2)
          E_P(iSort) = E_P(iSort) + &
               M_P(iSort) * c*P2/(c + sqrt(c2 + P2))
       end do
    end do
    if(nProc==1)then
       Energy_P = E_P
    else
       call MPI_Reduce(&
            E_P, Energy_P, nPType, MPI_REAL,MPI_SUM,  0, iComm, iError)
    end if
  end subroutine get_energy
  !==============================
  subroutine pass_energy
    use ModMpi
    use PIC_ModProc

    real :: E_P(nPType)
    !---------------------------
    if(nProc==1)return
    E_P = Energy_P
    call MPI_Reduce(&
         E_P, Energy_P, nPType, MPI_REAL,MPI_SUM,  0, iComm, iError)
  end subroutine pass_energy
  !=========================
  subroutine add_velocity_init
    use ModReadParam, ONLY: read_var
    real      :: W_D(Wx_:Wz_)
    integer   :: iSort
    real,dimension(:,:),pointer :: Coord_VI
    integer   :: iP
    !-----------------------------
    call read_var('iSort',iSort)
    call read_var('Wx'   ,W_D(Wx_))
    call read_var('Wy'   ,W_D(Wy_))
    call read_var('Wz'   ,W_D(Wz_))
    call set_pointer_to_particles(iSort,Coord_VI)

    do iP = 1,n_P(iSort)
       Coord_VI(Wx_:Wz_,iP) = Coord_VI(Wx_:Wz_,iP) + W_D
    end do
  end subroutine add_velocity_init
  !==============================================================
  subroutine add_velocity_sine
    use ModReadParam, ONLY: read_var
    real    :: Ampl,WaveNumber
    character(len=3) :: Direction
    integer :: iSort
    integer :: iP
    real,dimension(:,:),pointer :: Coord_VI
    !---------------------------------------
    call read_var('iSort',iSort)
    call read_var('Amplitude' ,Ampl)      ! u_0
    call read_var('WaveNumber',WaveNumber)! k
    call read_var('Direction' ,Direction)
    call set_pointer_to_particles(iSort,Coord_VI)

    select case(Direction)
    case('x')
       do iP=1,n_P(iSort)
          Coord_VI(Wx_,iP) = Coord_VI(Wx_,iP) + &
               Ampl*sin(WaveNumber*Coord_VI(x_,iP))
       end do

    case('y')
       do iP=1,n_P(iSort)
          Coord_VI(Wy_,iP) = Coord_VI(Wy_,iP) + &
               Ampl*sin(WaveNumber*Coord_VI(y_,iP))
       end do

    case('z')
       do iP=1,n_P(iSort)
          Coord_VI(Wz_,iP) = Coord_VI(Wz_,iP) + &
               Ampl*sin(WaveNumber*Coord_VI(z_,iP))
       end do

    case default
       !Add something here

    end select

  end subroutine add_velocity_sine
  
  !=======================PARTICLE MOVER=========================!
  !Advance the particles in one timestep; calculate cell-centered
  !number density and velocity if DoComputeMoments==.true.
  subroutine advance_particles(iSort,DoComputeMoments)
    use ModCoordTransform, ONLY: cross_product

    integer,intent(in) :: iSort
    logical,intent(in) :: DoComputeMoments

    real:: QDtPerM, M
    real,dimension(nDim)::QPerVDx_D
    real:: W2
    integer::iParticle, iBlock
    real,dimension(nDim)::X_D
    real,dimension(x_:z_)   ::W_D,W12_D,EForce_D,BForce_D
    real    :: Gamma
    real,dimension(:,:),pointer::Coord_VI
    integer:: iShift_D(nDim)
    !-------------------------------
    !Initialize the simulation for this sort of particles

    !1/2 * Q * dt / M, 
    !used to calculate the force effect
    QDtPerM = cHalf * Dt * Q_P(iSort) / M_P(iSort)

    !Q / V * (/\Delta x, \Delta y, \Delta z/)
    !Used to calculate J*dt in the 
    !charge-conserving scheme
    QPerVDx_D = (Q_P(iSort)/CellVolume) * Dx_D
    M = M_P(iSort)

    call set_pointer_to_particles(iSort,Coord_VI)

    !Looping over particles
    do iParticle=1,n_P(iSort)
       !Get coordinates and momentum
       X_D = Coord_VI(x_:nDim,iParticle)
       W_D = Coord_VI(Wx_:Wz_,iParticle)
       
       W2 = sum(W_D**2)
       Gamma = sqrt(c2 + W2)
       
       !Now W_D is the initial momentum, W2=W_D^2
       !Mow Gamma is the initial Gamma-factor multiplied by c
       
       !\
       ! Get block number
       !/
       iBlock = 1

       !call timing_start('formfactor')
       call get_form_factors(X_D,Node_D,HighFF_ID)
       !call timing_stop('formfactor')
       
       !Electric field force
       !call timing_start('electric')
       EForce_D = QDtPerM * e_interpolated_d(iBlock)
       !call timing_stop('electric')

       !Add kinetic energy
       Energy_P(iSort) = Energy_P(iSort) + &
            M*c*(W2/(Gamma+c) + sum(W_D*EForce_D)/Gamma)
      
       !Acceleration from the electric field, for the
       !first half of the time step:

       W_D = W_D + EForce_D

       !Get magnetic force
       !call timing_start('magnetic')
       BForce_D = QDtPerM/sqrt( c2+sum(W_D**2) ) * b_interpolated_d(iBlock)
       !call timing_stop('magnetic')       

       !Add a half of the magnetic rotation:

       W12_D = W_D + cross_product(W_D,BForce_D)

       !Multiply the magnetic force by 2 to take a whole
       !rotation and reduce its magnitude not to perturb energy

       BForce_D = (2.0/(1.0 + sum(BForce_D**2))) * BForce_D

       !Get a final momentum

       W_D = W_D + cross_product(W12_D,BForce_D) + EForce_D

       Gamma = sqrt(c2+sum(W_D**2))
       !Now Gamma is the final Gamma-factor multiplied by c

       !Save momentum
       Coord_VI(1+nDim:3+nDim,iParticle) = W_D
       W_D = (1.0/Gamma)*W_D
       !Now W_D is the velocity divided by c

       !Update coordinate
       X_D = X_D + SpeedOfLight_D*W_D(1:nDim)

       !New form factor
       !call timing_start('formfactor')
       call get_form_factors(X_D,NodeNew_D,HighFFNew_ID)
       !call timing_stop('formfactor')
       
       !Contribute to number density and velocity
       if(DoComputeMoments)&
            call add_DensityVelocity(W_D*c,NodeNew_D,HighFFNew_ID,iBlock)

       !Contribute to the current
       !call timing_start('current')
       call add_current(QPerVDx_D,W_D,iBlock)
       !call timing_stop('current')
       
       iShift_D = floor(X_D/nCell_D)
       X_D = X_D - nCell_D*iShift_D
       Coord_VI(1:nDim,iParticle) = X_D
       
       !To be done: for non-zero iShift_D, depending on the choice of 
       !the whole scheme and/or boundary conditions, some more action
       !may be needed.
    end do
  end subroutine advance_particles

end module PIC_ModParticles

