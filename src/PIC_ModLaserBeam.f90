module PIC_ModLaserBeam
  use PIC_ModMain, ONLY: tSimulation
  implicit none
contains
  subroutine laser_beam(iDir, x, y, z, EField)
    !\
    ! The electric field parameter EField has an intent inout
    ! The input value is found from 'noreflect' boundary condition 
    !/
    !\
    ! Calculate electric field iDir component in the point, x,y,z 
    !/
    integer,intent(in)        :: iDir
    real,   intent(in)        :: x,y
    real, optional, intent(in):: z !not used in 2D geometry
    real,   intent(inout)     :: EField !Before 
    !-----------
    
  end subroutine laser_beam
end module PIC_ModLaserBeam
