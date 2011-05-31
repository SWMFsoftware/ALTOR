! This file contains the top level methods for ALTOR
!==========================
subroutine PIC_setup
  use PIC_ModLogFile,   ONLY: open_logfile, nLogFile
  implicit none
  character(len=*), parameter :: NameSub='PIC_setup'
  !-------------------------
  if(nLogFile >=1) call open_logfile
 
end subroutine PIC_setup
!==========================
subroutine PIC_init_session
  character(len=*), parameter :: NameSub='PIC_init_session'
  !-------------------------
end subroutine PIC_init_session
!==========================
subroutine PIC_advance(tMax)
  use PIC_ModMain,      ONLY: tSimulation
  use PIC_ModField,     ONLY: update_magnetic, Rho_G, Counter_GD
  use PIC_ModParticles, ONLY: advance_particles, Energy_P, nPType
  use PIC_ModField,     ONLY: update_e, field_bc_periodic, &
                              current_bc_periodic
  use PIC_ModMain,      ONLY: tSimulation, iStep, Dt
  use PIC_ModLogFile,   ONLY: write_logfile, nLogFile
  use PIC_ModMpi,       ONLY: pass_current
  implicit none
  real,intent(in) :: tMax

  integer :: iSort

  character(len=*), parameter :: NameSub='PIC_advance'
  !-----------------
  if(tSimulation > tMax)return
  call timing_start('advance')
  !\
  ! Start update through the time step
  !/
  !1. Update the magnetic field through the half time step
  call timing_start('adv_b')
  call update_magnetic
  call timing_stop('adv_b')
  !2. Prepare to move particles
  Rho_G = 0.0; Counter_GD = 0.0; Energy_P = 0.0

  !3. Move particles
  call timing_start('adv_particles')
  do iSort = 1, nPType
     call advance_particles(iSort)
  end do
  call timing_stop('adv_particles')
  if(nLogFile>=1)then
     if(mod(iStep,nLogFile)==0)&
          call write_logfile
  end if

  
  !\
  !4. Update Magnetic field through a half timestep
  !/
  call timing_start('adv_b')
  call update_magnetic
  call timing_stop('adv_b')
  !5. Collect currents
  call pass_current
  
  !6. Advance electric field
  call timing_start('adv_e')
  call update_e
  call timing_stop('adv_e')
  call field_bc_periodic
  call timing_stop('advance')
  iStep = iStep + 1
  tSimulation = tSimulation + Dt
  
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