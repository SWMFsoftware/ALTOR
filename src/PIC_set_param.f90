!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine PIC_set_param(TypeAction)
  use ModReadParam
  use PIC_ModProc,    ONLY: iProc,iComm
  use PIC_ModMain
  use PIC_ModSize,    ONLY: nDim, MaxDim
  use PIC_ModParticles
  use PC_BATL_particles, ONLY: Particle_I
  use PIC_ModLogFile, ONLY: nLogFile, nToWrite, nToWrite_II
  use PIC_ModOutput, ONLY: nStepOut, nStepOutMin, TypeFile
  use PIC_ModThermal, ONLY: read_temperature, thermalize 
  use PIC_ModField,   ONLY: add_e, add_b, iGCN 
  use PIC_ModLaserBeam, ONLY: read_laser_beam
  use ModConst
  use PC_BATL_lib,ONLY: init_mpi, init_batl, nG
  use PC_BATL_mpi,ONLY:BATL_iProc=>iProc, BATL_nProc=>nProc
  use PIC_ModBatlInterface
  implicit none

  character (len=13) :: NameSub='PIC_set_param'
  
  ! Arguments

  ! TypeAction determines if we read or check parameters
  character (len=*), intent(in) :: TypeAction

  logical :: IsUninitialized      = .true.

  ! The name of the command
  character (len=lStringLine) :: NameCommand, StringLine, NameDescription
  integer :: iSession
  integer :: iDim, iP, iVar, iSide, nRootRead_D(nDim)=1
  character(LEN=10) :: NameNormalization
  real   :: Value_V(nDim+MaxDim) = 0

  integer           :: TimingDepth=-1
  character(len=10) :: TimingStyle='cumm'
  character(len=30) :: NameVar
  !------------------
  iSession = i_session_read()
  if(iGCN/=nG.and.iProc==0)then
     write(*,'(a)')'Reconfigure ALTOR using the command'
     write(*,'(a,i1)')'./Config.pl -ng=',iGCN
     call CON_stop('Code stopped')
  end if
  select case(TypeAction)
  case('CHECK','Check','check')
     if(iProc==0)write(*,*) NameSub,': CHECK iSession =',iSession
     !\
     ! Initialize timing
     !/
     if(iProc==0)then
        call timing_active(UseTiming)
        if(iSession==1)then
           call timing_step(0)
        end if
        call timing_depth(TimingDepth)
        call timing_report_style(TimingStyle)
     end if
     if(IsUninitialized)then
        IsUninitialized = .false.
        if(UseSharedField)then
           BATL_nProc = 1
           BATL_iProc = 0
        else
           call init_mpi(iComm)
        end if
        if(i_line_command("#GRID", iSessionIn = 1) > 0) then
           call init_batl(XyzMin_D, XyzMax_D, MaxBlock, 'cartesian', &
                TypeFieldBC_S(1:2*nDim-1:2) == 'periodic', nRootRead_D)
           Dx_D = (XyzMax_D - XyzMin_D)/(nRootRead_D(1:nDim)*nCell_D(1:nDim))
           DxInv_D = 1/Dx_D
           CellVolume = product(Dx_D); vInv = 1/CellVolume
           Q_P(Electron_) = - CellVolume/nPPerCellCrit
           M_P(Electron_) =   CellVolume/nPPerCellCrit
           do iP = 2, nPType
              Q_P(iP) = Q_P(iP)*CellVolume/nPPerCellCrit
              M_P(iP) = M_P(iP)*CellVolume/nPPerCellCrit
           end do
           call set_altor_grid
        end if
        if(UseUniform)call uniform
        if(UseFoil)call foil
        if(UseThermalization)call thermalize
     end if
     !\
     ! Initialize parameters
     !/
     SpeedOfLight_D = Dt * c / Dx_D
     RETURN
  case('read','Read','READ')
     if(iProc==0)then
        write(*,*) NameSub,': READ iSession =',iSession,&
             ' iLine=',i_line_read(),' nLine =',n_line_read()
        call read_echo_set(.true.)
     end if
     if(iSession ==1)then
        call set_particle_param
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
        do iDim = 1, nDim
           call read_var('Dx_D', Dx_D(iDim))
        end do
        CellVolume = product(Dx_D); vInv = 1/CellVolume
        DxInv_D = 1/Dx_D
        
     case('#NORMALIZATION')
        call read_var('Normalization',NameNormalization)
        select case(trim(NameNormalization))
        case('standard')
           !First charge is normalized to electron charge, mass is normalized 
           !to electron mass. Furthermore they are both divided by 
           !nPPerCrit electron frequency is 1.
           !c=1, (e/m)=1 for electron and \omega=1
           !For critical density nPPerCrit particles per
           !cell give the omega_Pe=\omega=1
           c = 1.0; c2 = 1.0
           call read_var('nPPerCellCrit',nPPerCellCrit)
           do iP=2, nPType
              write(NameVar,'(a,i1,a)') 'Q_P(',iP,')/| Q_P(Electron_) |'
              call read_var(NameVar,Q_P(iP))
              write(NameVar,'(a,i1,a)') 'M_P(',iP,')/ M_P(Electron_)'
              call read_var(NameVar,M_P(iP))
           end do
        case default
           call CON_stop(NameSub//':Unknown normalizaton type='&
                //NameNormalization)
        end select
     case('#UNIFORM')
        UseUniform = .true.
        call read_uniform

     case('#FOIL')
        UseFoil = .true.
        call read_foil

     case('#LASERBEAM')
        call read_laser_beam

     case('#THERMALIZE')
        UseThermalization = .true.
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
           call put_particle(iP,Value_V(1:nDim),1,Value_V(Wx_:Wz_))
           nToWrite_II(2,nToWrite) = iP
           nToWrite_II(3,nToWrite) = Particle_I(iP)%nParticle
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
                   'Periodic BC along the direction iDim=',iDim
              TypeFieldBC_S(2*iDim-1:2*iDim) = 'periodic'
           end if
        end do

     case("#GRID")
        if(.not.is_first_session())CYCLE READPARAM
        do iDim = 1, nDim
           call read_var('nRoot_I', nRootRead_D(iDim)) 
        end do
        do iDim = 1, nDim
           call read_var('CoordMin', XyzMin_D(iDim))
           call read_var('CoordMax', XyzMax_D(iDim))
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
contains
!===========================================================================

  logical function is_first_session()

    is_first_session = iSession == 1

    if(iSession/=1 .and. iProc==0)then
       write(*,*)NameSub//' WARNING: command '//trim(NameCommand)// &
            ' can be used in the first session only !!!'
       call CON_stop('Correct PARAM.in')
    end if

  end function is_first_session
end subroutine PIC_set_param
