!--------------------------------------------------------------!
module PIC_ModGrid
  integer,parameter::nDim=3            !Dimensionality

  integer,parameter::MaxBlock = 1

  integer,parameter::nX=10,nY=10,nZ=10 !The numbers of the grid 
                                       !cells
  integer,parameter::x_=1,y_=2,z_=3
  integer,parameter::nSortParticle=2,nElectronMax=10000000

end module PIC_ModGrid
