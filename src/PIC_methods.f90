! This file contains the top level methods for ALTOR
!==========================
subroutine PIC_setup
  use PIC_ModMain,      ONLY: IsPeriodicField_D, IsInitialized
  use PIC_ModLogFile,   ONLY: open_logfile, nLogFile
  implicit none
  character(len=*), parameter :: NameSub='PIC_setup'
  !-------------------------
  if(nLogFile >=1) call open_logfile
  IsInitialized = .true.
end subroutine PIC_setup
!==========================
subroutine PIC_init_session(iSession)
  use PIC_ModMain, ONLY:IsInitialized
  integer, intent(in) :: iSession
  character(len=*), parameter :: NameSub='PIC_init_session'
  !--------------------------------------------------------------------------
  if(.not.IsInitialized)call PIC_setup
end subroutine PIC_init_session
!==========================
subroutine PIC_advance(tMax)
  use PIC_ModMain,      ONLY: tSimulation
  use PIC_ModField,     ONLY: update_magnetic, Rho_G, Counter_GD
  use PIC_ModParticles, ONLY: advance_particles, Energy_P, nPType,&
                              pass_energy
  use PIC_ModField,     ONLY: update_e, field_bc, &
                              current_bc_periodic, field_bc_periodic
  use PIC_ModMain,      ONLY: tSimulation, iStep, Dt, IsPeriodicField_D
  use PIC_ModLogFile,   ONLY: write_logfile, nLogFile
  use PIC_ModOutput,    ONLY: write_density,write_field, nStepOutMin,nStepOut
  use PIC_ModMpi,       ONLY: pass_current
  implicit none
  real,intent(in) :: tMax

  integer :: iSort,ii

  character(len=*), parameter :: NameSub='PIC_advance'
  !-----------------

  if(tSimulation > tMax)return
  call timing_start('advance')
  !\
  ! Electric field and the particle coordinates are at the beginning of 
  ! the time step, the magnetic field and particle velocities are half time step
  ! behind.
  !/
  !\
  ! Start update through the time step
  !/
  !1. Update the magnetic field through the half time step
  call timing_start('adv_b')
  call update_magnetic
  call timing_stop('adv_b')



  !\
  ! Electromagnetic fields and the particle coordinates are at the beginning of 
  ! the time step, the particle velocities are half time step
  ! behind.
  !/
  !2. Prepare to move particles
  Counter_GD = 0.0; Energy_P = 0.0

  !3. Move particles
  call timing_start('adv_particles')
  do iSort = 1, nPType
     Rho_G = 0.0
     call advance_particles(iSort)
     if(nStepOut>=1.and.nStepOutMin<=iStep+1)then
        if(mod(iStep+1,nStepOut)==0)&
             call write_density(iSort)
     end if
  end do
  call pass_energy
  !\
  ! Electromagnetic fields and the particle energies are at the beginning of 
  ! the time step. Density and the particle coordinates are at the end of the timestep.
  ! The particle velocities are in the middle of the timestep.
  !/
  call timing_stop('adv_particles')
  if(nLogFile>=1)then
     !All energies in the logfile are at the beginning of the timestep.
     !For test particles velocities are in the midlle of the timestep, coordinates are
     !at the end of it
     if(mod(iStep,nLogFile)==0)&
          call write_logfile
  end if
  
  !\
  !4. Update Magnetic field through a half timestep
  !/
  call timing_start('adv_b')
  call update_magnetic
  call timing_stop('adv_b')
  !\
  ! Electric fields are at the beginning of 
  ! the time step. The particle coordinates are at the end of the timestep.
  ! The particle velocities and magnetic fields are in the middle of the timestep.
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

  !call field_bc_absorption
  !if(laser_pulse > 0) call laser_pulse_init
  !\
  ! Electric fields are at the beginning of 
  ! the time step. The particle coordinates are at the end of the timestep.
  ! The particle velocities and magnetic fields are in the middle of the timestep.
  !/
  if(any(IsPeriodicField_D))call field_bc_periodic

  call timing_stop('advance')

  do ii = 1, 3
     if(nStepOut>=1.and.nStepOutMin<=iStep)then
        if(mod(iStep,nStepOut)==0)&
             call write_field(ii)
     end if
  end do
  ! 
end subroutine PIC_advance
!==============================

subroutine PIC_save_files(TypeSaveIn)
  use ModUtilities, ONLY : upper_case
  implicit none
  character(len=*), intent(in) :: TypeSaveIn
  character(len=len(TypeSaveIn)) :: TypeSave

  character(len=*), parameter :: NameSub='PIC_save_files'
  !-----------------------------------------

  TypeSave = TypeSaveIn
  call upper_case(TypeSave)
end subroutine PIC_save_files
!==============================
subroutine PIC_finalize
  use PIC_ModLogFile, ONLY: nLogFile, close_logfile
  implicit none
  character(len=*), parameter :: NameSub='PIC_finalize'
  !--------------------
  if(nLogFile >=1)call close_logfile
end subroutine PIC_finalize
!=============================
!=====================================================================
subroutine PC_user_specify_region(iArea, iBlock, nValue, NameLocation, &
     IsInside, IsInside_I, Value_I)
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

!=====================================================================
