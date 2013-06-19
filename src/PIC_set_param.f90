subroutine PIC_set_param(TypeAction)
  use ModReadParam
  use PIC_ModProc,    ONLY: iProc
  use PIC_ModMain
  use PIC_ModSize,    ONLY: nDim
  use PIC_ModParticles
  use PIC_ModLogFile, ONLY: nLogFile, nToWrite, nToWrite_II
  use PIC_ModThermal, ONLY: read_temperature
  use PIC_ModField,   ONLY: add_e, add_b
  use ModConst
  implicit none

  character (len=13) :: NameSub='PIC_set_param'
  
  ! Arguments

  ! TypeAction determines if we read or check parameters
  character (len=*), intent(in) :: TypeAction

  ! The name of the command
  character (len=lStringLine) :: NameCommand, StringLine, NameDescription
  integer :: iSession
  integer :: iDim, iP, iVar, iSide
  
  character(LEN=10) :: NameNormalization
  integer:: nPPerCrit

  real   :: Value_V(nDim+3)

  integer            :: TimingDepth=-1
  character (len=10) :: TimingStyle='cumm'
  !------------------
  iSession = i_session_read()
  select case(TypeAction)
  case('CHECK')
     if(iProc==0)write(*,*) NameSub,': CHECK iSession =',iSession
     SpeedOfLight_D = Dt * c / Dx_D
     if(iProc==0)then
        call timing_active(UseTiming)
        if(iSession==1)then
           call timing_step(0)
        end if
        call timing_depth(TimingDepth)
        call timing_report_style(TimingStyle)
     end if

     RETURN
  case('read','Read','READ')
     if(iProc==0)then
        write(*,*) NameSub,': READ iSession =',iSession,&
             ' iLine=',i_line_read(),' nLine =',n_line_read()
        call read_echo_set(.true.)
     end if
  case default
     call CON_stop(NameSub//': TypeAction='//TypeAction// &
          ' must be "CHECK" or "READ"!')
  end select
 !\
  ! Read parameters from the text
  !/
  READPARAM: do
     if(.not.read_line(StringLine) )then
        IsLastRead = .true.
        EXIT READPARAM
     end if
     if(.not.read_command(NameCommand)) CYCLE READPARAM

     select case(NameCommand)
     case('#DESCRIPTION')
        call read_var('',NameDescription)
     case('#CHECKGRID')
        call read_var('nX',iP)
        if(iP /= nX)call CON_stop('nX differs, reconfigure ALTOR')
        call read_var('nY',iP)
        if(iP /= nY)call CON_stop('nY differs, reconfigure ALTOR')
        if(nDim==3)then
           call read_var('nZ',iP)
           if(iP /= nZ)call CON_stop('nZ differs, reconfigure ALTOR')
        end if
        call read_var('nPType',iP)
        if(iP /= nPType)call CON_stop('nPType differs, reconfigure ALTOR')
     case('#DXYZ')
        do iDim=1,nDim
           call read_var('Dx_D(iDim)', Dx_D(iDim))
        end do
        CellVolume = product(Dx_D)
        
     case('#NORMALIZATION')
        call read_var('Normalization',NameNormalization)
        select case(trim(NameNormalization))
        case('none','CGS')
           c = cLightSpeed * 0.010
           c2= c*c
           do iP = 1, nPType
              call read_var('Q_P(iP)',Q_P(iP))
              call read_var('M_P(iP)',M_P(iP))
           end do
        case('standard')
           !c=1, (e/m)=1 for electron and \omega=1
           !For critical density nPPerCell particles per
           !cell give the omega_Pe=\omega=1

           c = 1.0; c2 = 1.0
           call read_var('nPPerCrit',nPPerCrit)
           Q_P(Electron_) = - CellVolume/(4*cPi*nPPerCrit)
           M_P(Electron_) =   CellVolume/(4*cPi*nPPerCrit)
           do iP=2, nPType
              call read_var('Q_P/| Q_P(Electron_) |',Q_P(iP))
              Q_P(iP) = Q_P(iP) *  CellVolume/(4*cPi*nPPerCrit)
              call read_var('M_P/|M_P(Electron_)',M_P(iP))
              M_P(iP) = M_P(iP) *  CellVolume/(4*cPi*nPPerCrit)
           end do
        case default
           call CON_stop(NameSub//':Unknown normalizaton type='&
                //NameNormalization)
        end select
        call set_particle_param(M_P,Q_P)
     case('#UNIFORM')
        call read_uniform

     case('#FOIL')
        call read_foil

     case('#THERMALIZE')
        call read_temperature
        
     case('#ADDE')
        call add_e

     case('#ADDB')
        call add_b

     case('#ADDVELOCITY')
        call add_velocity

     case('#TIMESTEP')
        call read_var('Dt',Dt)
    
     case('#TESTPARTICLE')
        call read_var('iSort', iP)
        do iVar = 1, 3 + nDim
           call read_var('Coords(iVar)', Value_V(iVar))
        end do
        if (iProc==0)then
           nToWrite = nToWrite +1
           call put_particle(iP,Value_V(1:nDim))
           Of(iP)% Coords(Wx_:Wz_,n_P(iP)) = Value_V(Wx_:Wz_)
           nToWrite_II(2,nToWrite) = iP
           nToWrite_II(3,nToWrite) = n_P(iP)
        end if
  
     case('#FIELDBC')
        do iSide = 1, 2*nDim 
           call read_var('TypeFieldBc_S(iSide)',TypeFieldBc_S(iSide))
        end do
        !\
        !Assign array IsPeriodicField_D
        !/
        IsPeriodicField_D = .false.
        do iDim=1,nDim
           if(any(TypeFieldBC_S(2*iDim-1:2*iDim)=='periodic'))then
              IsPeriodicField_D(iDim)=.true.
              if(iProc==0)write(*,*)&
                   'Periodic Boundary Conditions for fields along the direction iDim=',iDim
              TypeFieldBC_S(2*iDim-1:2*iDim) = 'periodic'
           end if
        end do

     case("#END")
        IslastRead=.true.
        EXIT READPARAM

     case("#RUN")
        IslastRead=.false.
        EXIT READPARAM

     case("#STOP")
        call read_var('MaxIteration',nIter)
        call read_var('tSimulationMax',tMax)

     case("#CPUTIMEMAX")
        call read_var('CpuTimeMax',CpuTimeMax)

     case("#CHECKSTOPFILE")
        call read_var('UseStopFile',UseStopFile)

     case("#TIMING")
        if(iSession /= 1)CYCLE READPARAM
        call read_var('UseTiming',UseTiming)
        if(.not.UseTiming)CYCLE READPARAM
        call read_var('DnTiming',nTiming)
        call read_var('nDepthTiming',TimingDepth)
        call read_var('TypeTimingReport',TimingStyle)
     
     case('#LOGFILE')
        if(iSession /=1 )CYCLE READPARAM
        call read_var('nLogFile',nLogFile)

     case default
        if(iProc==0) then
           write(*,*) NameSub // ' WARNING: unknown #COMMAND ' // &
                trim(NameCommand),' !!!'
           call CON_stop('Correct PARAM.in!')
        end if

     end select
  end do READPARAM
end subroutine PIC_set_param
