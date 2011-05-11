!^CFG COPYRIGHT UofM
module PIC_ModParticles
  use PIC_ModGrid,ONLY:nPType, nElectronMax
  !  use PIC_ModMain,ONLY:SaveVelocity
  use PIC_ModMain,ONLY:c,c2,Dt,Dx_D,CellVolume
  use PIC_ModParticleInField,ONLY: rho_G,add_density,add_current
  use PIC_ModParticleInField,ONLY: b_interpolated_d,e_interpolated_d
  use PIC_ModFormFactor,ONLY: HighFF_ID, HighFFNew_ID,&
                              Node_D,NodeNew_D 
  use PIC_ModProc,ONLY:iProc,iError
  use PIC_ModMpi
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
       Of(iSort)%Coords=cZero
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
  end subroutine put_particle
  
  !==========================
  
  subroutine add_velocity(W_D, iSort)
    real,    intent(in) :: W_D(Wx_:Wz_)
    integer, intent(in) :: iSort

    real,dimension(:,:),pointer::Coord_VI
    integer :: iP
    !------------------
    call set_pointer_to_particles(iSort,Coord_VI)
    do iP = 1,n_P(iSort)
       Coord_VI(Wx_:Wz_,iP) = Coord_VI(Wx_:Wz_,iP) + W_D
    end do
  end subroutine add_velocity
   
  !================PARTICLE MOVER================!
  subroutine advance_particles(iSort)

    integer,intent(in)::iSort

    real:: QDtPerM
    real,dimension(nDim)::QPerVDx_D

    !real :: QcDtPerV

    real::Energy,W2

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

    !Nullify counters: for energy:
    !NOT A GENERAL WAY - DISABLED
    !Should be nullified outside the routine.
    !Energy_P(iSort) = 0.0
    !rho_G = 0.0

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

       call get_form_factors(X_D,Node_D,HighFF_ID)

       !Electric field force
       EForce_D = QDtPerM * e_interpolated_d()

       !Add kinetic energy

       Energy_P(iSort) = Energy_P(iSort) + &
            M*c*(W2/(Gamma+c) + sum(W_D*EForce_D)/Gamma)

       !Acceleration from the electric field, for the
       !first half of the time step:

       W_D = W_D + EForce_D

       !Get magnetic force

       BForce_D = QDtPerM/sqrt( c2+sum(W_D**2) ) * b_interpolated_d()

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
       W_D = (cOne/Gamma)*W_D
       
       !end if

       !Now W_D is the velocity divided by c
       !Update coordinate
       X_D = X_D + SpeedOfLight_D*W_D(1:nDim)

       !New form factor
       call get_form_factors(X_D,NodeNew_D,HighFFNew_ID)
 
       !Contribute to the charge densiry
       call add_density(NodeNew_D,HighFFNew_ID)
   
       !      if(nDim<3)then
       !         W_D=QcDtPerV*W_D !For nDim=3 velocity is not used
       !         call add_current(QPerVDx_D,W_D)
       !      end if
       !Contribute to the current
       call add_current(QPerVDx_D)

       iShift_D = floor(X_D/nCell_D)
       X_D = X_D -nCell_D*iShift_D
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
