module PIC_ModFormFactor
  use PIC_ModSize, ONLY: nDim
  implicit none

  integer, parameter :: lOrderFF =  2 
  integer, parameter :: iDownFF  = -1, iUpFF = 1
  integer, parameter :: iExt     =  1

  !Structures
  !Discsrete particle coordinates
  integer,dimension(nDim) :: Node_D, NodeNew_D
 
  !Particle form-factors
  real,dimension(1:lOrderFF+1,nDim)   :: HighFF_ID, HighFFNew_ID
               

  real,dimension(0:lOrderFF+1,nDim)   :: LowFF_ID
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

    integer,dimension(nDim) :: NodeAux_D, iShift_D
    real                    :: Aux_I(lOrderFF+1)
    !--------------------
    NodeOut_D = floor(X_D)

    d_D = X_D - real(NodeOut_D)

    NodeOut_D = NodeOut_D + 1

    !Asterisk means non-zero FF value
    !Node-2!Node-1!Node  !Node+1!           
    !      !      !------!      ! '----' means the particle position
    !      0*     1*     2*     3* Low
    !      !      !      !      !
    !      !  1*  !  2*  !  3*  ! High
    !      !      !      !      !
    !   Node-2  Node-1 Node   Node+1
    do iDim=1,nDim
   
       HighFFOut_ID(1, iDim) = cHalf *( 1.0 - d_D(iDim))**2
       HighFFOut_ID(2, iDim) = 0.750 - (cHalf - d_D(iDim))**2
       HighFFOut_ID(3, iDim) = cHalf * d_D(iDim)**2

       LowFF_ID(0,iDim)      = cHalf * HighFFOut_ID(1, iDim)
       LowFF_ID(1,iDim)      = cHalf * (HighFFOut_ID(1, iDim) +&
                                        HighFFOut_ID(2, iDim)  )
       LowFF_ID(2,iDim)      = cHalf * (HighFFOut_ID(2, iDim) +&
                                        HighFFOut_ID(3, iDim)  )
       LowFF_ID(3,iDim)      = cHalf * HighFFOut_ID(3, iDim)
    end do

  end subroutine get_form_factors
end module PIC_ModFormFactor
