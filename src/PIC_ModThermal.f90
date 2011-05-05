module PIC_ModThermal
  use PIC_ModMain,ONLY:nSortParticle,SaveVelocity
  use PIC_ModRandom
  implicit none
  private!Except
  real,dimension(nSortParticle)::vT2_P
  !Methods
  public::init_thermal
  public::thermalize
  interface thermalize
     module procedure thermalize_particle
     module procedure thermalize_all_sorts
  end interface
contains
  subroutine init_thermal(vT2In_P)
    real,dimension(nSortParticle),intent(in)::vT2In_P
    vT2_P=vT2In_P
  end subroutine init_thermal
  !----------------------------------------------------------!
  subroutine thermalize_particle(iP,iSort)
    use PIC_ModParticles
    integer,intent(in)::iP,iSort
    real::Energy,MomentumAvr,GammaInv,W_D(x_:z_)
    real,parameter::cTwoThird=cTwo*cThird
    Energy=-vT2_P(iSort)*LOG(ERAND())
    MomentumAvr=sqrt(cTwoThird*Energy)
    W_D(x_)=MomentumAvr*COS(cTwoPi*RAND())
    W_D(y_)=MomentumAvr*COS(cTwoPi*RAND())
    W_D(z_)=MomentumAvr*COS(cTwoPi*RAND())
    if(SaveVelocity)then
       GammaInv=c/sqrt(c2+sum(W_D**2))
       W_D=W_D*GammaInv
    end if
    Of(iSort)%Coords(Wx_:Wz_,iP)=W_D
  end subroutine thermalize_particle
  !----------------------------------------------------------!
  subroutine thermalize_all_sorts(iP)
    use PIC_ModParticles
    integer,intent(in)::iP
    integer::iSort
    do iSort=1,nSortParticle
       call thermalize_particle(iP,iSort)
    end do
  end subroutine thermalize_all_sorts
end module PIC_ModThermal
