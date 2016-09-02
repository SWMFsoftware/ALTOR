subroutine PIC_set_param(TypeAction)
  use ModReadParam
  use PIC_ModProc,    ONLY: iProc
  use PIC_ModMain
  use PIC_ModSize,    ONLY: nDim
  use PIC_ModParticles
  use PC_BATL_particles, ONLY: Particle_I
  use PIC_ModLogFile, ONLY: nLogFile, nToWrite, nToWrite_II
  use PIC_ModOutput, ONLY: nStepOut, nStepOutMin, TypeFile
  use PIC_ModThermal, ONLY: read_temperature
  use PIC_ModField,   ONLY: add_e, add_b 
  use PIC_ModLaserBeam, ONLY: read_laser_beam
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

  integer           :: TimingDepth=-1
  character(len=10) :: TimingStyle='cumm'
  character(len=30) :: NameVar
  !------------------
  iSession = i_session_read()
  select case(TypeAction)
  case('CHECK','Check','check')
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
        if(nDim==3)then
           call read_var('Dx_D(1)', Dx_D(1))
           call read_var('Dx_D(2)', Dx_D(2))
           call read_var('Dx_D(3)', Dx_D(3))
        else
           call read_var('Dx_D(1)', Dx_D(1))
           call read_var('Dx_D(2)', Dx_D(2))
        end if
        
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
           !First charge is normalized to electron charge, mass is normalized 
           !to electron mass. Furthermore they are both divided by 
           !4*cPi*nPPerCrit electron frequency is 1.
           !c=1, (e/m)=1 for electron and \omega=1
           !For critical density nPPerCell particles per
           !cell give the omega_Pe=\omega=1

           c = 1.0; c2 = 1.0
           call read_var('nPPerCrit',nPPerCrit)
           !hyzhou: I don`t understand why there is CellVolume here.
           !I checked the PARTICLE MOVER, and saw that there is actually
           !no difference if I remove it.
           !The reason is that this is inconsistent with the normalization.
           !I should ask Igor to confirm this.
           Q_P(Electron_) = - 1.0/(4*cPi*nPPerCrit)
           M_P(Electron_) =   1.0/(4*cPi*nPPerCrit)
           do iP=2, nPType
              write(NameVar,'(a23,i1,a1)') 'Q_P/| Q_P(Electron_) |(',iP,')'
              call read_var(NameVar,Q_P(iP))
              Q_P(iP) = Q_P(iP) /(4*cPi*nPPerCrit)
              write(NameVar,'(a20,i1,a1)') 'M_P/ M_P(Electron_)(',iP,')'
              call read_var(NameVar,M_P(iP))
              M_P(iP) = M_P(iP) /(4*cPi*nPPerCrit)
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

     case('#LASERBEAM')
        call read_laser_beam

     case('#THERMALIZE')
        call read_temperature
        
     case('#ADDE')
        call add_e

     case('#ADDB')
        call add_b

     case('#ADDVELOCITY')
        call add_velocity_init

     case('#ADDSINEWAVEVELOCITY')
        call add_velocity_sine

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
           Particle_I(iP)%State_VI(Wx_:Wz_,n_P(iP)) = Value_V(Wx_:Wz_)
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

     case('#SAVEFIELDMOMENTS')
        call read_var('TypeFile',TypeFile)
        call read_var('nStepOutMin', nStepOutMin)
        call read_var('nStepOut', nStepOut)

     case default
        if(iProc==0) then
           write(*,*) NameSub // ' WARNING: unknown #COMMAND ' // &
                trim(NameCommand),' !!!'
           call CON_stop('Correct PARAM.in!')
        end if

     end select
  end do READPARAM
end subroutine PIC_set_param
