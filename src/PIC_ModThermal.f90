Module PIC_ModThermal
  use PIC_ModMain,ONLY:nPType
  use PIC_ModRandom
  use PIC_ModParticles
  implicit none

  real,dimension(nPType) :: vT2_P = 0.0
contains
  subroutine thermalize_particle(iP,iSort)
    integer,intent(in)::iP,iSort
    real::Energy, MomentumAvr
    integer :: iW
    !==================================

    Energy=-vT2_P(iSort)*LOG(RAND())

    MomentumAvr=sqrt(2.0*Energy)
    do iW = Wx_, Wz_
       Of(iSort)%Coords(iW,iP) = Of(iSort)%Coords(iW,iP) + &
            MomentumAvr*COS(cTwoPi*RAND())
    end do
  end subroutine thermalize_particle
  !===============
  subroutine thermalize
    integer:: iSort, iP
    !------------------
    do iSort = 1, nPType
       do iP = 1, n_P(iSort)
          call thermalize_particle(iP, iSort)
       end do
    end do
  end subroutine thermalize
  !=============
end module PIC_ModThermal
