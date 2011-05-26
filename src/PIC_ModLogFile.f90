module PIC_ModLogFile
  use PIC_ModParticles, ONLY: nPType
  implicit none
  SAVE
  PRIVATE !Except

  !\
  !Unit for the log file
  !/
  integer :: iLogUnit = -1

  !\
  !Writing frequency for the logfile
  !/
  integer, public :: nLogFile = -1

  !\
  ! Starting energy (to be used for normalization)
  !/
  real :: Energy0 = -1.0


  public :: open_logfile
  public :: write_logfile
  public :: close_logfile
  
  integer, parameter ::         &
       Time_ = 1,               &
       PartFirst_ = 2,          &
       PartLast_  = 1 + nPType, &
       Bx_        = 2 + nPType, &
       Ez_        = 7 + nPType, &
       ETotal_    = 8 + nPType

  real  :: Value_V(Time_:ETotal_)
  character(LEN=30) ::NameFormat
contains
  !======================
  subroutine open_logfile
     use ModIoUnit,        ONLY: io_unit_new
     use PIC_ModParticles, ONLY: get_energy
     use PIC_ModProc,      ONLY: iProc
     use PIC_ModGrid,      ONLY: nPType
  
  character(LEN=100) :: Name
  character(LEN=1  ) :: Name1
  
  integer :: iP
  !---------------------
     !\
     !Open log file
     !/
     if(iProc==0)then
        write(Name,'(a,i4.4,a)')'log_n',nLogFile,'.out'
        iLogUnit = io_unit_new()
        open(iLogUnit,file=Name,status='replace')
        write(Name,'(a)')'iStep Time '
        do iP = 1, nPType
           write(Name1,'(i1.1)')iP
           Name = trim(Name)//'Energy'//Name1//' '
        end do
        Name = trim(Name)//'Bx2 By2 Bz2 Ex2 Ey2 Ez2 ETotal'
        write(iLogUnit,'(a)') Name
     end if

     !\
     !Calculate the particle energies
     !/
     call get_energy
     call write_logfile
  end subroutine open_logfile
  !===========================
  subroutine write_logfile
    use PIC_ModField,     ONLY: get_field_energy
    use PIC_ModMain,      ONLY: iStep, tSimulation
    use PIC_ModParticles, ONLY: Energy_P
    use ModUtilities, ONLY: flush_unit
    !-----------------------
    Value_V = 0.0
    Value_V(Time_) = tSimulation
    Value_V(PartFirst_:PartLast_) = Energy_P
    call get_field_energy(Value_V(Bx_:Ez_))
    Value_V(ETotal_) = sum(Value_V(PartFirst_:Ez_))
    if(Energy0 <=0)then
       !Initial stage, save ETotal and use it for normalization
       Energy0 =  Value_V(ETotal_)
       if(Energy0 <=0)Energy0 = 1.0
       !Make format statement
       if(nPType==1)then
          write(NameFormat,'(a)')'(i10,9es13.5)'
       else
          write(NameFormat,'(a,i2,a)')'(i10,',ETotal_,'es13.5)'
       end if
    end if
    Value_V(ETotal_) =  Value_V(ETotal_) / Energy0
    write(iLogUnit,trim(NameFormat))iStep, Value_V
    call flush_unit(iLogUnit)
  end subroutine write_logfile
  !===========================
  subroutine close_logfile
    use PIC_ModProc, ONLY: iProc
    !--------------------
    if(iProc==0)close(iLogUnit)
  end subroutine close_logfile
end module PIC_ModLogFile
