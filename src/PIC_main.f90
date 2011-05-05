program PIC
  use PIC_ModProc
  use PIC_ModParticles
  implicit none
  !---------------
  call init_mpi
  write(*,*)iProc,nProc
  call set_particle_param((/cOne,cOne/),(/cOne,cOne/))
  call init_field
  call clean_mpi
end program PIC
!********************
!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  use PIC_ModProc
  implicit none
  character (len=*), intent(in) :: StringError
  write(*,*)StringError
  call clean_mpi
  stop
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test
