module PIC_ModLaserBeam
  implicit none
contains
  subroutine laser_beam(iDir, x, y, z, EField)
    !\
    ! Calculate electric field iDir component in the oint, x,y,z 
    !/
    integer,intent(in)        :: iDir
    real,   intent(in)        :: x,y
    real, optional, intent(in):: z !not used in 2D geometry
    real,   intent(inout)     :: EField !Before 
    !-----------
    
  end subroutine laser_beam
end module PIC_ModLaserBeam
