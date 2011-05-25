program PIC
  use PIC_ModProc
  use PIC_ModParticles
  use PIC_ModRandom
  use PIC_ModMain
  use ModReadParam
  use ModIoUnit,ONLY:UNITTMP_
  implicit none

  integer:: iSession = 1
  logical:: IsFound

  real(Real8_) :: CpuTimeStart

  !---------------
  !Initialize MPI
  call MPI_INIT(iError)
  iComm=MPI_COMM_WORLD
  call MPI_COMM_RANK(iComm,iProc,iError)
  call MPI_COMM_SIZE(iComm,nProc,iError)


  !\
  ! Initialize time which is used to check CPU time
  !/
  CpuTimeStart = MPI_WTIME()


  !\
  ! Delete BATSRUS.SUCCESS and BATSRUS.STOP files if found
  !/
  if(iProc==0)then

     inquire(file='BATSRUS.SUCCESS',EXIST=IsFound)
     if(IsFound)then
        open(UNITTMP_, file = 'ALTOR.SUCCESS')
        close(UNITTMP_,STATUS = 'DELETE')
     end if

     inquire(file='ALTOR.STOP',EXIST=IsFound)
     if (IsFound) then
        open(UNITTMP_, file = 'ALTOR.STOP')
        close(UNITTMP_, STATUS = 'DELETE')
     endif

  end if

  !\
  ! Read input parameter file. Provide the default restart file for #RESTART
  !/
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

  
     iSession=iSession+1

     TIMELOOP: do
        if(stop_condition_true())EXIT TIMELOOP
        if(tMax > 0.0)then
           call PIC_advance(tMax)
        else
           call PIC_advance(huge(0.0))
        end if
     end do TIMELOOP
     if(iProc==0) &
          write(*,*)'----- End of Session   ',iSession,' ------'   
     if(IsLastRead)EXIT SESSIONLOOP
  end do SESSIONLOOP
 
  !\
  ! Touch BATSRUS.SUCCESS
  !/
  if(iProc==0)then
     open(UNITTMP_, file = 'ALTOR.SUCCESS')
     close(UNITTMP_)
  end if


  !Finalize MPI
  call MPI_Finalize(iError)
contains
   function stop_condition_true() result(StopConditionTrue)

    logical :: StopConditionTrue

    StopConditionTrue = .false.

    
    if(nIter >= 0 .and. iIter >= nIter) StopConditionTrue = .true.
    if(tMax > 0.0 .and. tSimulation >= tMax) StopConditionTrue = .true.

  end function stop_condition_true
  !===============================

  function is_time_to_stop() result(IsTimeToStop)

    use PIC_ModMain, ONLY: CpuTimeMax, UseStopFile

    logical :: IsTimeToStop

    IsTimeToStop = .false.

    if(iProc==0)then
       if(CpuTimeMax > 0.0 .and. MPI_WTIME()-CpuTimeStart >= CpuTimeMax)then
          write(*,*)'CPU time exceeded:',CpuTimeMax,MPI_WTIME()-CpuTimeStart
          IsTimeToStop=.true.
       end if
       if(.not.IsTimeToStop .and. UseStopFile) then
          inquire(file='ALTOR.STOP',exist=IsTimeToStop)
          if (IsTimeToStop) &
               write(*,*)'ALTOR.STOP file exists: recieved stop signal'
       end if
    end if

    call MPI_BCAST(IsTimeToStop,1,MPI_LOGICAL,0,iComm,iError)

  end function is_time_to_stop

end program PIC

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
