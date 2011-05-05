module PIC_ModProc
  use ModMpi
  implicit none
  integer::iProc !PE rank
  integer::nProc !Number of PEs
  integer::iComm !Communicator
  integer::iError!Auxiliary integer to store an error number
  !Methods
  public::init_mpi !Intitialize MPI, iProc,nProc,iComm
  public::clean_mpi
contains
!-----------------------------------------------------------!
  subroutine init_mpi
    call MPI_INIT(iError)
    iComm=MPI_COMM_WORLD
    call MPI_COMM_RANK(iComm,iProc,iError)
    call MPI_COMM_SIZE(iComm,nProc,iError)
  end subroutine init_mpi
!-----------------------------------------------------------!
  subroutine clean_mpi
    call MPI_Finalize(iError)
    stop
  end subroutine clean_mpi
!-----------------------------------------------------------!
end module PIC_ModProc
    
    
