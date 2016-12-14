module PIC_ModMain
  use ModNumConst
  use PIC_ModSize, ONLY: nDim, MaxBlock
  implicit none
  SAVE
  !-----------------------
  integer:: nBlockMax
  !If UseSharedField==.true. at all processors
  !the field is the same, only particles are distributed.
  logical, parameter:: UseSharedField = .true.
 
  logical:: IsInitialized = .false.

  logical:: UseVectorPotential = .false.

  ! logical::DoAccelerateLight = .false.

  real :: c = cHalf, c2 = 0.25
  real,dimension(nDim):: SpeedOfLight_D = cHalf

  !\
  ! Time stepping parameters and values.
  !/
  real :: tSimulation = 0.0
  integer :: iStep = 0
  real :: Dt= 0.0

  !\
  ! Grid spacing
  !/
  real :: CellVolume=1.0, Dx_D(nDim)=1.0, DxInv_D(nDim) = 1.0
  !\
  ! Grid limits
  !/
  real :: XyzMin_D(nDim) = -0.50, XyzMax_D(nDim) = 0.50
  !\
  ! Stopping conditions. These variables are only used in stand alone mode.
  !/
  real    :: tMax = -1.0, CpuTimeMax = -1.0
  integer :: nIter = -1
  logical :: UseStopFile = .true.
  logical :: IsLastRead=.false.

  !\
  !Progress variables
  !/

  integer :: nProgress1=10, nProgress2=100

  
  !Timing variables
  logical:: UseTiming = .true.
  integer:: nTiming = -2

  !\
  ! Boundary condition for the electric field
  !/
  !Read by the command
  !#FIELDBC
  !Example:
  !
  !#FIELDBC
  !laserbeam            TypeFieldBC_S(x<0)
  !noreflect            TypeFieldBC_S(x>nX.dx)
  !periodic             TypeFieldBC_S(y<0)
  !periodic             TypeFieldBC_S(y>nY.dy)
  !periodic             TypeFieldBC_S(z<0)
  !periodic             TypeFieldBC_S(z>nZ.dz)
  !
  !
  !Implemented types: periodic, noreflect, laserbeam
  
  character(LEN=10) :: TypeFieldBC_S(1:2*nDim)='periodic'
  
  
  logical:: IsPeriodicField_D(nDim) = .true.
  
end module PIC_ModMain
