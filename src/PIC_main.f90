program PIC
  use PIC_ModProc
  use PIC_ModParticles
  call init_mpi
  write(*,*)iProc,nProc
  call set_particle_param((/cOne,cOne/),(/cOne,cOne/))
  call init_field
  call clean_mpi
end program PIC
