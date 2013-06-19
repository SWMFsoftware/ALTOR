program PIC
  use PIC_ModProc
  use PIC_ModParticles
  use PIC_ModRandom
  use PIC_ModMain
  use ModReadParam
  use ModIoUnit,ONLY:UNITTMP_
  use ModMpi
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
  ! Delete PIC.SUCCESS and PIC.STOP files if found
  !/
  if(iProc==0)then

     inquire(file='PIC.SUCCESS',EXIST=IsFound)
     if(IsFound)then
        open(UNITTMP_, file = 'PIC.SUCCESS')
        close(UNITTMP_,STATUS = 'DELETE')
     end if

     inquire(file='PIC.STOP',EXIST=IsFound)
     if (IsFound) then
        open(UNITTMP_, file = 'PIC.STOP')
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

     !\
     ! Time execution (timing parameters were set by MH_set_parameters)
     !/
     if(iSession==1)then
        call timing_start('PIC')
        call timing_start('setup')
        call PIC_setup
        call PIC_init_session
        call timing_stop('setup')
        if(nTiming > -3)call timing_report_total
        if(iProc==0)write(*,*)'Resetting timing counters after setup.'
        call timing_reset('#all',3)
     else
        call PIC_init_session
     end if
  


     TIMELOOP: do
        if(stop_condition_true())EXIT TIMELOOP
        if(is_time_to_stop())exit SESSIONLOOP
        call timing_step(iStep + 1)
        if(tMax > 0.0)then
           call PIC_advance(tMax)
        else
           call PIC_advance(huge(0.0))
        end if

        call show_progress

     end do TIMELOOP
     if(IsLastRead)EXIT SESSIONLOOP
     if(iProc==0) &
          write(*,*)'----- End of Session   ',iSession,' ------'   
     iSession=iSession+1
     if (nTiming > -2) call timing_report
     call timing_reset_all
  end do SESSIONLOOP

  if(iProc==0)then
     write(*,*)
     write(*,'(a)')'    Finished Numerical Simulation'
     write(*,'(a)')'    -----------------------------'
  end if

  if (nTiming > -2) call timing_report

 
  call PIC_save_files('FINALWITHRESTART')
 

  if(iProc==0)then
     write(*,*)
     write(*,'(a)')'    Finished Saving Output Files'
     write(*,'(a)')'    ----------------------------'
  end if

  call timing_stop('PIC')

  if(nTiming > -3)call timing_report_total

  call PIC_finalize
  !\
  ! Touch PIC.SUCCESS
  !/
  if(iProc==0)then
     open(UNITTMP_, file = 'PIC.SUCCESS')
     close(UNITTMP_)
  end if


  !Finalize MPI
  call MPI_Finalize(iError)
  stop
contains
   function stop_condition_true() result(StopConditionTrue)

    logical :: StopConditionTrue

    StopConditionTrue = .false.

    
    if(nIter >= 0 .and. iStep >= nIter) StopConditionTrue = .true.
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
          inquire(file='PIC.STOP',exist=IsTimeToStop)
          if (IsTimeToStop) &
               write(*,*)'PIC.STOP file exists: recieved stop signal'
       end if
    end if
    if(nProc==1)return
    call MPI_BCAST(IsTimeToStop,1,MPI_LOGICAL,0,iComm,iError)

  end function is_time_to_stop
!===========================================================================

  subroutine show_progress

    use PIC_ModMain,      ONLY: nProgress1, nProgress2
    use PIC_ModMain,      ONLY: Dt

    real(Real8_), external :: timing_func_d
    real(Real8_) :: CpuTimePIC,CpuTimeAdvance

    !\
    ! Show timing results if required
    !/

    ! Show speed as cells/second/PE/step
    if( UseTiming .and. iProc==0 &
         .and. nProgress1>0 .and. mod(iStep,nProgress1) == 0 ) then

       CpuTimePIC = timing_func_d('sum',3,'PIC','PIC')
       CpuTimeAdvance=timing_func_d('sum',1,'advance','PIC')
  
       write(*,'(a,f9.1,a,f9.1,a,i8,a,1p,e10.4,a)') 'Speed is',&
            sum(n_P) &
            /max(1.D-10,CpuTimeAdvance),&
            ' p/s/pe after',&
            CpuTimePIC,&
            ' s at N =',iStep, ' (', tSimulation,' s)'
      
    endif

    ! Show timing tables
    if(nTiming>0.and.mod(iStep,nTiming)==0) then
       call timing_report
    else if(nProgress2>0.and.mod(iStep,nProgress2) == 0) then
       call timing_tree(2,2)
    end if

  end subroutine show_progress

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
