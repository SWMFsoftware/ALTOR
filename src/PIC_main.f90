program PIC
  use PIC_ModProc
  use PIC_ModParticles
  use PIC_ModRandom
  use PIC_ModMain
  use ModReadParam
  implicit none

  integer:: iSession = 1
  !---------------
  !Initialize MPI
  call MPI_INIT(iError)
  iComm=MPI_COMM_WORLD
  call MPI_COMM_RANK(iComm,iProc,iError)
  call MPI_COMM_SIZE(iComm,nProc,iError)
  call read_file('PARAM.in',iComm)

  SESSIONLOOP: do
     call read_init('  ', iSessionIn=iSession)

     if(iProc==0)&
          write(*,*)'----- Starting Session ',iSession,' ------'
     !\
     ! Set and check input parameters for this session
     !/
     call PIC_set_param('READ')
     call PIC_set_param('CHECK')
  if(IsLastRead)EXIT SESSIONLOOP
  
  if(iProc==0) &
       write(*,*)'----- End of Session   ',iSession,' ------'
  iSession=iSession+1

  end do SESSIONLOOP
  stop
  call init_rand
 


  !Finalize MPI
  call MPI_Finalize(iError)

end program PIC
subroutine clean_mpi

  use PIC_ModProc
  implicit none
  !-----------------
  
  stop
end subroutine clean_mpi
!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  use ModMpi
  use PIC_ModProc
  implicit none
  !-------------------
  character (len=*), intent(in) :: StringError
  write(*,'(a)')StringError
  
  !Finalize MPI
  call MPI_Finalize(iError)
  
  stop
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test
