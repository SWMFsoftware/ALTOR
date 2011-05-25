module PIC_ModMain
  use ModNumConst
  use PIC_ModGrid
  implicit none
  !-----------------------
  !If UseSharedField==.true. at all processaros
  !The field is the same, only particles are distributed.
  logical, parameter:: UseSharedField = .true.
 

  logical::UseVectorPotential = .false.

  ! logical::DoAccelerateLight = .false.

  real :: c = cHalf, c2 = 0.25
  real,dimension(nDim)::SpeedOfLight_D = cHalf

  !\
  ! Time stepping parameters and values.
  !/
  real :: tSimulation = 0.0
  integer :: n_step, iIter =0
  real :: Dt= 0.0
  real ::CellVolume=cOne, Dx_D(nDim)=cOne

  !\
  ! Stopping conditions. These variables are only used in stand alone mode.
  !/
  real    :: tMax = -1.0, CpuTimeMax = -1.0
  integer :: nIter = -1
  logical :: UseStopFile = .true.
  logical :: IsLastRead=.false.
end module PIC_ModMain
