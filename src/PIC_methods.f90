!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! This file contains the top level methods for ALTOR
!==========================
subroutine PIC_setup
  use PIC_ModProc,      ONLY: iComm
  use PIC_ModMain,      ONLY: UseSharedField, UseUniform, UseFoil
  use PIC_ModLogFile,   ONLY: open_logfile, nLogFile
  use PIC_ModOutput,    ONLY: PIC_save_files
  use PIC_ModParticles, ONLY: uniform, foil, DoAddVelocity_P, &
       add_initial_velocity, nPType

  implicit none
  integer :: iSort
  character(len=*), parameter :: NameSub='PIC_setup'
  !------------------------------------------------
  if(UseUniform)call uniform
  if(UseFoil   )call foil
  if(any(DoAddVelocity_P))then
     do iSort = 1, nPType
        if(DoAddVelocity_P(iSort))&
             call add_initial_velocity(iSort)
     end do
  end if
  !Save the initial outputs
  call timing_start('output')
  if(nLogFile >=1) call open_logfile
  call PIC_save_files('INITIAL')
  call timing_stop('output')
end subroutine PIC_setup
!==========================
subroutine PIC_init_session(iSession)
  use PIC_ModMain, ONLY:IsInitialized
  integer, intent(in) :: iSession
  character(len=*), parameter :: NameSub='PIC_init_session'
  !-------------------------------------------------------------------
  if(.not.IsInitialized)then
     call PIC_setup
     IsInitialized = .false.
  end if
end subroutine PIC_init_session
!==============================
subroutine PIC_advance(tMax)
  use PIC_ModField,     ONLY: update_magnetic, Current_GDB, State_VGBI
  use PIC_ModParticles, ONLY: Energy_P, nPType, pass_energy
  use PIC_ModParticles, ONLY: advance_particles
  use PIC_ModField,     ONLY: update_e, field_bc, E_GDB, iGCN 
  use PIC_ModMain,      ONLY: tSimulation, iStep, Dt, IsPeriodicField_D
  use PIC_ModLogFile,   ONLY: write_logfile, nLogFile
  use PIC_ModOutput,    ONLY: nStepOutMin, nStepOut
  use PIC_ModMpi,       ONLY: pass_current, pass_density, pass_moments
  use PC_BATL_pass_face_field, ONLY: message_pass_field
  implicit none
  real,intent(in) :: tMax
  integer :: iSort, iBlock
  character(len=*), parameter :: NameSub='PIC_advance'
  !--------------------------------------------------------------------
  if(tSimulation > tMax)return
  call timing_start('advance')
  !\
  ! Electric field and the particle coordinates are at the beginning of 
  ! the time step, the magnetic field and particle velocities are half
  ! time step behind.
  !/
  !\
  ! Start update through the time step
  !/
  !1. Update the magnetic field through the half time step
  call timing_start('adv_b')
  call update_magnetic
  call timing_stop('adv_b')
  !\
  ! Electromagnetic fields and the particle coordinates are at the
  !beginning of the time step, the particle velocities are half 
  !time step behind.
  !/
  !2. Prepare to move particles
  Current_GDB = 0.0; Energy_P = 0.0
  !3. Move particles
  call timing_start('adv_particles')

  State_VGBI = 0.0
  if(nStepOut>=1.and.nStepOutMin<=iStep+1&
       .and.mod(iStep+1,nStepOut)==0)then
     do iSort=1, nPType
        !Calculate cell-centered number density and velocity while
        !advancing the particles
        call advance_particles(iSort,DoComputeMoments=.TRUE.)
        !Save the moments
        call pass_moments(iSort)
     end do
  else
     do iSort=1, nPType
        !Advance the particles without calculating cell-centered
        !number velocity 
        call advance_particles(iSort,DoComputeMoments=.FALSE.)
        call pass_density(iSort)
     end do
  end if

  call pass_energy
  !\
  ! Electromagnetic fields and the particle energies are at  
  ! the beginning of the time step. Density and the particle 
  ! coordinates are at the end of the timestep.
  ! The particle velocities are in the middle of the timestep.
  !/
  call timing_stop('adv_particles')
  if(nLogFile>=1)then
     !All energies in the logfile are at the beginning of the timestep.
     !For test particles velocities are in the middle of the timestep,
     !coordinates are at the end of it
     if(iStep/=0.and.mod(iStep,nLogFile)==0)&
          call write_logfile
  end if
  !\
  !4. Update Magnetic field through a half timestep
  !/
  call timing_start('adv_b')
  call update_magnetic
  call timing_stop('adv_b')
  !\
  ! Electric fields are at the end of the time step.
  ! The particle coordinates are at the end of the timestep.
  ! The particle velocities and magnetic fields are in  
  ! the middle of the timestep.
  !/
  !5. Collect currents
  call pass_current

  !6. Advance electric field
  !\
  ! 6.1 calculate time dependent boundary conditions and 
  !/
  iStep = iStep + 1
  tSimulation = tSimulation + Dt
  call field_bc

  call timing_start('adv_e')
  call update_e
  call timing_stop('adv_e')

  !\
  ! Electric fields are at the end of the time step.
  ! The particle coordinates are at the end of the timestep.
  ! The particle velocities and magnetic fields are in  
  ! the middle of the timestep.
  !/
  call message_pass_field(iGCN, E_GDB)

  call timing_stop('advance')

end subroutine PIC_advance
!==============================
subroutine PIC_finalize
  use PIC_ModLogFile, ONLY: nLogFile, close_logfile
  use PIC_ModProc, ONLY: iProc
  use PIC_ModOutput, ONLY: PIC_save_files
  implicit none
  character(len=*), parameter :: NameSub='PIC_finalize'
  !--------------------
  if(nLogFile >=1)call close_logfile
  call timing_start('save_final')
  call PIC_save_files('FINAL') 
  call timing_stop('save_final')

  if(iProc==0)then
     write(*,*)
     write(*,'(a)')'    Finished Saving Output Files'
     write(*,'(a)')'    ----------------------------'
  end if
end subroutine PIC_finalize
!=====================================================================
subroutine PC_user_specify_region(iArea, iBlock, nValue, &
     NameLocation, IsInside, IsInside_I, Value_I)
  implicit none

  integer,   intent(in):: iArea        ! area index in BATL_region
  integer,   intent(in):: iBlock       ! block index
  integer,   intent(in):: nValue       ! number of output values
  character, intent(in):: NameLocation ! c, g, x, y, z, or n

  logical, optional, intent(out) :: IsInside
  logical, optional, intent(out) :: IsInside_I(nValue)
  real,    optional, intent(out) :: Value_I(nValue)

  character(len=*), parameter :: NameSub = 'user_specify_region'
  !-------------------------------------------------------------------
end subroutine PC_user_specify_region

