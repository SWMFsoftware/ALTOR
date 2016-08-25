Module PIC_ModThermal
  use PIC_ModSize,      ONLY: nPType
  use PIC_ModMain,      ONLY: c2
  use PIC_ModRandom
  use PIC_ModParticles, ONLY: Wx_, Wz_, n_P, parallel_init_rand
  use PIC_ModParticles, ONLY: Energy_P, nTotal_P, M_P
  use PC_BATL_particles, ONLY: Particle_I
  use ModNumConst,      ONLY: cTwoPi
  implicit none

  real,dimension(nPType) :: uTh_P = 0.0
contains
  !Generate Gaussian distribution using Box-Muller method
  subroutine thermalize_particle(iP,iSort)
    integer,intent(in)::iP,iSort
    real::Energy, MomentumAvr
    integer :: iW
    !-----------------------------------
    do iW = Wx_, Wz_
       Energy = -uTh_P(iSort)**2 *LOG(RAND())
       MomentumAvr = sqrt(2.0*Energy)
       Particle_I(iSort)%State_VI(iW,iP) = Particle_I(iSort)%State_VI(iW,iP)+&
            MomentumAvr*COS(cTwoPi*RAND())
    end do
  end subroutine thermalize_particle
  !==================================
  subroutine thermalize(iSort)
    integer,intent(in):: iSort
    integer:: iP
    !------------------
    do iP = 1, n_P(iSort)
       call thermalize_particle(iP, iSort)
    end do
  end subroutine thermalize
  !===========Reading command #THERMALIZE============
  subroutine read_temperature
    use ModReadParam,    ONLY: read_var
    use PIC_ModProc,     ONLY: iProc
    use PIC_ModParticles,ONLY: get_energy
    integer:: iSort
    character(len=8) :: NameVar
    !--------------
    do iSort = 1,nPType
       write(NameVar,'(a6,i1,a1)') 'uTh_P(',iSort,')'
       !uTh_P is the thermal velocity in normalized units.
       call read_var(NameVar,uTh_P(iSort))
       call parallel_init_rand(6*n_P(iSort),iSort)
       call thermalize(iSort)
    end do
    call get_energy
    if(iProc==0)then
       do iSort = 1, nPType
          write(*,*)'For particle of sort ',iSort,&
               ' averaged energy is ',Energy_P(iSort)/nTotal_P(iSort),&
               ' should be ', 1.50 * M_P(iSort) * uTh_P(iSort)**2 *&
               (1.0 -1.250*uTh_P(iSort)**2/c2)
       end do
    end if
  end subroutine read_temperature
  !==============================
end module PIC_ModThermal
