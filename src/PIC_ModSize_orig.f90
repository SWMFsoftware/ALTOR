!--------------------------------------------------------------!
module PIC_ModSize
  use PC_BATL_size, ONLY: nX=>nI, nY=>nJ, nZ=>nK, jDim_, kDim_
  use PC_BATL_size, ONLY: nPType=>nKindParticle, nDim, MaxDim
  use PC_BATL_size, ONLY:nBlock
  integer, parameter :: x_=1,y_=2,z_=3
  integer, parameter :: MaxBlock=1
  integer, parameter :: nElectronMax=10000000
  integer, parameter :: nCell_D(MaxDim) = (/nX, nY, nZ/)

end module PIC_ModSize
