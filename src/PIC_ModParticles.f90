!^CFG COPYRIGHT UofM
module PIC_ModParticles
  use PIC_ModGrid,ONLY:nSortParticle,nElectronMax
  use PIC_ModMain,ONLY:SaveVelocity,UseBunemanScheme
  use PIC_ModMain,ONLY:c,c2,Dt,Dx_D,CellVolume
  use PIC_ModParticleInField
  use PIC_ModProc,ONLY:iProc,iError
  use PIC_ModMpi
  implicit none
  integer,parameter::W_=nDim,Wx_=W_+x_,Wy_=W_+y_,Wz_=W_+z_
  integer,parameter::Electrons_=1,Electron_=Electrons_
  !Structures
  real,dimension(nSortParticle)::M_P,Q_P
  !Particle's mass M_P, dt*Q_P, dt*Q_P/M_P
  type particles
     real,dimension(:,:),pointer::Coords
     !Particle coordinates in the phase space. First nDim components 
     !are the cartezian coordinates, the other are velocities, or
     !the momentum components (if .not.SaveVelocity)
  end type particles
  type(particles),dimension(nSortParticle)::Of
 
  integer,dimension(nSortParticle)::n_P,nMax_P
  real,dimension(nSortParticle)::Energy_P
  real,dimension(nSortParticle)::OmegaPDtMax_P
  !Methods
  public::set_pointer_to_particles
  public::set_particle_param
  public::put_particle
  public::advance_particles
contains
  !--------------------------------------------------------------!
  subroutine set_pointer_to_particles(iSort,PointerToSet)
    integer,intent(in)::iSort
    real,dimension(:,:),pointer::PointerToSet
    PointerToSet=>Of(iSort)%Coords
  end subroutine set_pointer_to_particles
  !--------------------------------------------------------------!
  subroutine set_particle_param(MIn_P,QIn_P)
    real,dimension(nSortParticle),intent(in)::MIn_P,QIn_P
    logical::DoInit=.true.
    integer::iSort
    M_P=MIn_P
    Q_P=QIn_P
    n_P=0
    if(DoInit)then
       DoInit=.false.
    else
       if(iProc==0)write(*,*)&
            'Particle arrays are deallocated, information is lost'
       do iSort=Electrons_,nSortParticle
          deallocate(Of(iSort)%Coords)
       end do
    end if
    !Max number of particles is max number of electrons per the
    !charge ratio
    nMax_P=nElectronMax*nint(abs(Q_P(Electron_)/Q_P))
    
    do iSort=Electrons_,nSortParticle
       nullify(Of(iSort)%Coords)
       allocate(Of(iSort)%Coords(&
            x_:Wz_,1:nMax_P(iSort)),&
            stat=iError)
       if(iError>0)write(*,*)'Cannot allocate of sort=',iSort,&
            'particles at PE=',iProc
       Of(iSort)%Coords=cZero
    end do
  end subroutine set_particle_param
  !--------------------------------------------------------------!
  subroutine put_particle(iSort,PhaseCoords_D)
    integer,intent(in)::iSort
    real,dimension(nDim),intent(in)::PhaseCoords_D
    n_P(iSort)=n_P(iSort)+1
    Of(iSort)%Coords(x_:nDim,n_P(iSort))=PhaseCoords_D
  end subroutine put_particle
  !--------------------------------------------------------------!
  subroutine advance_particles(iSort)
    integer,intent(in)::iSort
    real::M,QDt,QcDtPerV,QDtPerM,Energy,W2
    integer::iParticle
    real,dimension(nDim)::X_D,QPerVDx_D
    real,dimension(x_:z_)   ::W_D,W12_D,EForce_D,BForce_D
    real::Gamma
    real,dimension(:,:),pointer::CoordPhaseSp
    M=M_P(iSort)
    QDt=Dt*Q_P(iSort)
    QcDtPerV=QDt*c/CellVolume
    QDtPerM=cHalf*QDt/M
    QPerVDx_D=Q_P(iSort)/CellVolume*Dx_D
    Energy_P(iSort)=cZero
    rho_G=cZero
    call set_pointer_to_particles(iSort,CoordPhaseSp)
    do iParticle=1,n_P(iSort)
       X_D=CoordPhaseSp(x_:nDim,iParticle)
       W_D=CoordPhaseSp(Wx_:Wz_,iParticle)
       if(SaveVelocity)then
          Gamma=cOne/sqrt(c2-sum(W_D**2))
          W_D=Gamma*c*W_D
          W2=Gamma**2-c2
       else
          W2=sum(W_D**2)
          Gamma=sqrt(c2+W2)
       end if
       !Now W_D is the initial momentum, W2=W_D^2
       !Mow Gamma is the initial Gamma-factor multiplied by c

       call get_form_factors(X_D,Node_D,HighFF_ID)

       !Electric field force
       EForce_D=QDtPerM*e_interpolated_d()

       !Add kinetic energy

       Energy_P(iSort)=Energy_P(iSort)+&
            M*c*(W2/(Gamma+c)+sum(W_D*EForce_D)/Gamma)

       !Acceleration from the electric field, for the
       !first half of the time step:

       W_D=W_D+EForce_D

       !Get magnetic force

       BForce_D=(QDtPerM/sqrt(c2+sum(W_D**2)))*b_interpolated_d()

       !Add a half of the magnetic rotation:

       W12_D=W_D+cross_product(W_D,BForce_D)

       !Multiply the magnetic force by 2 to take a whole
       !rotation and reduce its magnitude not to perturb energy

       BForce_D=(2.0/(1.0 + sum(BForce_D**2)))*BForce_D

       !Get a final momentum

       W_D=W_D+cross_product(W12_D,BForce_D)+EForce_D

       Gamma=sqrt(c2+sum(W_D**2))
       !Mow Gamma is the final Gamma-factor multiplied by c

       if(SaveVelocity)then
          W_D=(cOne/Gamma)*W_D
          CoordPhaseSp(Wx_:Wz_,iParticle)=c*W_D
       else
          CoordPhaseSp(1+nDim:3+nDim,iParticle)=W_D
          W_D=(cOne/Gamma)*W_D
       end if
       !Now W_D is the velocity divided by c
       !Update coordinate
       X_D=X_D+SpeedOfLight_D*W_D(1:nDim)
       CoordPhaseSp(1:nDim,iParticle)=X_D

       !New form factor
       call get_form_factors(X_D,NodeNew_D,HighFFNew_ID)
       call add_density(NodeNew_D,HighFFNew_ID)
       if(UseBunemanScheme)CYCLE
       if(nDim<3)then
          W_D=QcDtPerV*W_D !For nDim=3 velocity is not used
          call add_current(QPerVDx_D,W_D)
       else
          call add_current(QPerVDx_D)
       end if
    end do
    call pass_density
    if(iProc==0)then
       call get_rho_max(OmegaPDtMax_P(iSort))
       OmegaPDtMax_P(iSort)=sqrt(QDt**2/M*OmegaPDtMax_P(iSort))
    end if
  contains
    function cross_product(a,b)
      real,dimension(3)::cross_product
      real,dimension(3),intent(in)::a,b
      cross_product(1)=a(2)*b(3)-a(3)*b(2)
      cross_product(2)=a(3)*b(1)-a(1)*b(3)
      cross_product(3)=a(1)*b(2)-a(2)*b(1)
    end function cross_product
  end subroutine advance_particles
  !--------------------------------------------------------------!
end module PIC_ModParticles
