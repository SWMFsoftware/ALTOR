module PIC_ModFormFactor
  use PIC_ModGrid
  implicit none

  integer, parameter :: lOrderFF = 2 
  integer, parameter :: iDownFF = -1, iUpFF = 1

  !Structures
  !Discsrete particle coordinates
  integer,dimension(nDim) :: Node_D, NodeNew_D
 
  !Particle form-factors
  real,dimension(1:lOrderFF+1,nDim) :: HighFF_ID, HighFFNew_ID
               

  real,dimension(lOrderFF,nDim)::LowFF_ID
contains
  subroutine get_form_factors(X_D, NodeOut_D, HighFFOut_ID)
    use ModNumConst,ONLY: cHalf
    !Input parameters
    !X_D should be mormalized by (/\Delta x,\Delta y,\Delta z/)
    real,dimension(nDim),intent(in) :: X_D

    !Output parameters
    !The cell-center index    
    integer,dimension(nDim),intent(out)::NodeOut_D
    real,dimension(lOrderFF+1,nDim)::HighFFOut_ID
    !Normally both Node_D and HighFF_ID are members of this module
    real,dimension(nDim)::d_D
    integer::iDim
    !--------------------
    NodeOut_D = floor(X_D)

    d_D = X_D - real(NodeOut_D)

    Node_D=Node_D+1
    do iDim=1,nDim
       LowFF_ID(1,iDim) =1.0 - d_D(iDim)
       LowFF_ID(2,iDim) =      d_D(iDim)

       HighFFOut_ID(1, iDim) = LowFF_ID(1,iDim)**2
       HighFFOut_ID(3, iDim) = LowFF_ID(2,iDim)**2

       HighFFOut_ID(2,iDim)=2.0 - &
            (HighFFOut_ID(1,iDim) + HighFFOut_ID(3,iDim)) 
    end do

    HighFF_ID= HighFF_ID * cHalf
  end subroutine get_form_factors
end module PIC_ModFormFactor
