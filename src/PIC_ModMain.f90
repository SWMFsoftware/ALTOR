module PIC_ModMain
  use ModNumConst
  use PIC_ModGrid
  implicit none
  logical::UseBunemanScheme =.false.
  logical::SaveVelocity = .false.
  logical::UseVectorPotential = .false.
  logical::DoAccelerateLight = .false.
  real,parameter::c=cHalf,c2=c*c
  real,parameter,dimension(nDim)::SpeedOfLight_D=c
  real::Dt
  real,parameter::CellVolume=cOne,Dx_D(nDim)=cOne
end module PIC_ModMain
