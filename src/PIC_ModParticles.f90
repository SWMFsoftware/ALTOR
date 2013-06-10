!^CFG COPYRIGHT UofM
module PIC_ModParticles
  use PIC_ModSize,ONLY: nPType, nElectronMax,nDim, x_, y_, z_
  use PIC_ModSize,ONLY: nX, nY, nZ, nCell_D
  use PIC_ModMain,ONLY: c, c2, Dt, Dx_D, CellVolume, SpeedOfLight_D
  use PIC_ModParticleInField,ONLY: rho_G,add_density,add_current
  use PIC_ModParticleInField,ONLY: b_interpolated_d,e_interpolated_d
  use PIC_ModParticleInField,ONLY: min_val_rho, max_val_rho, rho_avr
  use PIC_ModFormFactor,ONLY: HighFF_ID, HighFFNew_ID,&
                              Node_D, NodeNew_D, get_form_factors
  use PIC_ModProc,      ONLY:iProc,iError
  use ModNumConst,      ONLY: cHalf
  implicit none
  integer,parameter :: W_ = nDim
  integer,parameter :: Wx_ = W_+x_, Wy_ = W_+y_, Wz_ = W_+z_

  integer,parameter::Electrons_=1,Electron_=Electrons_

  !Structures
  real,dimension(nPType) :: M_P, Q_P
  !Particle's mass and charge

  type particles
     real,dimension(:,:),pointer::Coords
     !Particle coordinates in the phase space. First nDim components 
     !are the cartezian coordinates, the other are velocities, or
     !the momentum components (if .not.SaveVelocity)
  end type particles
  type(particles),dimension(nPType) :: Of
 
  integer,dimension(nPType) :: n_P,nMax_P
  real,dimension(nPType)    :: Energy_P
  real,dimension(nPType)    :: OmegaPDtMax_P

  !Only at the root PE:
  integer,dimension(nPType) :: nTotal_P

  !Methods
  public::set_pointer_to_particles !Set pointer to the coordinate array of electrons or ions
  public::set_particle_param       !Assigns M_P, Q_P and allocates coordinate arrays
  public::put_particle             !Add particle with known coordinates
  public::advance_particles        
contains
  !=============================
  subroutine set_pointer_to_particles(iSort,PointerToSet)
    integer,intent(in)::iSort
    real,dimension(:,:),pointer::PointerToSet
    nullify(PointerToSet)
    PointerToSet=>Of(iSort)%Coords
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
          deallocate(Of(iSort)%Coords)
       end do
    end if
    !Max number of particles is max number of electrons per the
    !charge ratio
    nMax_P= nint(nElectronMax * abs(Q_P(Electron_)/Q_P) )
    
    do iSort=Electrons_,nPType
       nullify(Of(iSort)%Coords)
       allocate(Of(iSort)%Coords(&
            x_:Wz_,1:nMax_P(iSort)),&
            stat=iError)
       if(iError>0)write(*,*)'Cannot allocate of sort=',iSort,&
            'particles at PE=',iProc
       Of(iSort)%Coords = 0.0
    end do
  end subroutine set_particle_param
  
  !===============
  
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
    Of(iSort)%Coords(x_:nDim,n_P(iSort)) = PhaseCoords_D
    Of(iSort)%Coords(nDim+1:nDim+3,n_P(iSort)) = 0.0

  end subroutine put_particle
  !===========================
  !== This routine allows to get reporducible random disdtribution independent 
  !== on the number of processors invoved. This is extremely expensive way,
  !== as long as the operation of the random number calculation actually becomes
  !== serial 
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
    !--------------
    do iSort = 1,nPType
       call read_var('nPPerCell',nPPerCell_P(iSort))
       if(nPPerCell_P(iSort)<0.and.iProc==0)&
            write(*,*)'Particles will neutralize electrons'
    end do
    rho_G = 0.0
    call uniform(nPPerCell_P)
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
  end subroutine read_uniform
  !================================
  subroutine uniform(nPPerCellIn_P)
    use PIC_ModProc
    use PIC_ModRandom
    integer, intent(in) :: nPPerCellIn_P(nPType)
    integer             :: nPPerCell_P(nPType)
    integer :: nPPerPE, nResidual, NPTotal, iSort, iDim, iP
    real    :: Coord_D(nDim)
    logical :: UseQuasineutral
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

       NPTotal = product(nCell_D) * NPPerCell_P(iSort)
       nPPerPE = NPTotal/nProc; nResidual = nPTotal - nProc*nPPerPE

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
          call get_form_factors(Coord_D,Node_D,HighFF_ID)
          call add_density(Node_D,HighFF_ID)
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
    rho_G = 0.0
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
      integer :: nPPerPE, nResidual, NPTotal, iSort, iDim, iP
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
         nPPerPE = NPTotal/nProc; nResidual = nPTotal - nProc*nPPerPE
         
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
            if(  Coord_D(1) .lt. 0.0 .or. Coord_D(1) .ge. real(nCell_D(1)) .or.&
                 Coord_D(2) .lt. 0.0 .or. Coord_D(2) .ge. real(nCell_D(2)) .or.&
                 Coord_D(3) .lt. 0.0 .or. Coord_D(3) .ge. real(nCell_D(3)) &
                 ) cycle PART
            call put_particle(iSort, Coord_D)
            call get_form_factors(Coord_D,Node_D,HighFF_ID)
            call add_density(Node_D,HighFF_ID)
         end do PART
    end do
  end subroutine foil
  !
end subroutine read_foil
  !====================
  subroutine get_energy
    use ModMpi
    use PIC_ModProc

    real :: E_P(nPType), P2
    integer:: iSort, iP
    real,dimension(:,:),pointer::Coord_VI
    !--------------
    E_P = 0
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
    !-----------------

    if(nProc==1)return
    E_P = Energy_P
    call MPI_Reduce(&
         E_P, Energy_P, nPType, MPI_REAL,MPI_SUM,  0, iComm, iError)
  end subroutine pass_energy
  !=========================
  subroutine add_velocity
    use ModReadParam, ONLY: read_var
    real      :: W_D(Wx_:Wz_)
    integer   :: iSort

    real,dimension(:,:),pointer::Coord_VI
    integer :: iP
    !------------------
    call read_var('iSort',iSort)
    call read_var('Wx'   ,W_D(Wx_))
    call read_var('Wy'   ,W_D(Wy_))
    call read_var('Wz'   ,W_D(Wz_))
    call set_pointer_to_particles(iSort,Coord_VI)
    do iP = 1,n_P(iSort)
       Coord_VI(Wx_:Wz_,iP) = Coord_VI(Wx_:Wz_,iP) + W_D
    end do
  end subroutine add_velocity
   
  !================PARTICLE MOVER================!
  subroutine advance_particles(iSort)
    integer,intent(in)::iSort

    real:: QDtPerM, M
    real,dimension(nDim)::QPerVDx_D

    !real :: QcDtPerV

    !real::Energy,
    real::W2

    integer::iParticle

    real,dimension(nDim)::X_D

    real,dimension(x_:z_)   ::W_D,W12_D,EForce_D,BForce_D
    real    :: Gamma
    real,dimension(:,:),pointer::Coord_VI
    integer:: iShift_D(nDim)
    !-------------------------------
    !    Initialize the simulation for this sort of particles
    

    !Q * dt * c / V (dx*dy*dz)
    !QcDtPerV=QDt*c/CellVolume

    !1/2 * Q * dt / M, 
    !used to calculate the force effect

    QDtPerM = cHalf * Dt * Q_P(iSort) / M_P(iSort)

    !Q / V * (/\Delta x, \Delta y, \Delta z/)
    !Used to calculate J*dt in the 
    !charge-conserving scheme

    QPerVDx_D = (Q_P(iSort)/CellVolume) * Dx_D

    M = M_P(iSort)

    call set_pointer_to_particles(iSort,Coord_VI)

    !Start the loop over particles
    do iParticle=1,n_P(iSort)
       !Get coordinates and velocities!
       X_D=Coord_VI(x_:nDim,iParticle)
       W_D=Coord_VI(Wx_:Wz_,iParticle)

       
       W2 = sum(W_D**2)
       Gamma = sqrt(c2 + W2)
       
       !Now W_D is the initial momentum, W2=W_D^2
       !Mow Gamma is the initial Gamma-factor multiplied by c
       
       !call timing_start('formfactor')
       call get_form_factors(X_D,Node_D,HighFF_ID)
       !call timing_stop('formfactor')

       !Electric field force
       !call timing_start('electric')
       EForce_D = QDtPerM * e_interpolated_d()
       !call timing_stop('electric')

       !Add kinetic energy

       Energy_P(iSort) = Energy_P(iSort) + &
            M*c*(W2/(Gamma+c) + sum(W_D*EForce_D)/Gamma)
      
       !Acceleration from the electric field, for the
       !first half of the time step:

       W_D = W_D + EForce_D

       !Get magnetic force
       !call timing_start('magnetic')
       BForce_D = QDtPerM/sqrt( c2+sum(W_D**2) ) * b_interpolated_d()
       !call timing_stop('magnetic')       

       !Add a half of the magnetic rotation:

       W12_D = W_D + cross_product(W_D,BForce_D)

       !Multiply the magnetic force by 2 to take a whole
       !rotation and reduce its magnitude not to perturb energy

       BForce_D=(2.0/(1.0 + sum(BForce_D**2))) * BForce_D

       !Get a final momentum

       W_D=W_D + cross_product(W12_D,BForce_D)+EForce_D

       Gamma = sqrt(c2+sum(W_D**2))
       !Mow Gamma is the final Gamma-factor multiplied by c


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
      
       !Contribute to the charge densiry
       !call timing_start('density')
       call add_density(NodeNew_D,HighFFNew_ID)
       !call timing_stop('density')
       !      if(nDim<3)then
       !         W_D=QcDtPerV*W_D !For nDim=3 velocity is not used
       !         call add_current(QPerVDx_D,W_D)
       !      end if
       !Contribute to the current
       !call timing_start('current')
       call add_current(QPerVDx_D,W_D)
       !call timing_stop('current')

       iShift_D = floor(X_D/nCell_D)
       X_D = X_D - nCell_D*iShift_D
       Coord_VI(1:nDim,iParticle) = X_D
       
       !To be done: for non-zero iShift_D, depending on the choice of 
       !the whole scheme and/or boundary conditions, some more action
       !may be needed.

    end do
  contains
    !========================
    function cross_product(a,b)
      real,dimension(3)::cross_product
      real,dimension(3),intent(in)::a,b
      !----------------
      cross_product(1)=a(2)*b(3)-a(3)*b(2)
      cross_product(2)=a(3)*b(1)-a(1)*b(3)
      cross_product(3)=a(1)*b(2)-a(2)*b(1)
    end function cross_product
    !=======================
  end subroutine advance_particles
  !--------------------------------------------------------------!
end module PIC_ModParticles
