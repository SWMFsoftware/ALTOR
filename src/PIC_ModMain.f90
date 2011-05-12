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
  real :: Dt= 0.0
  real ::CellVolume=cOne, Dx_D(nDim)=cOne

  logical:: IsLastRead=.false.
end module PIC_ModMain
