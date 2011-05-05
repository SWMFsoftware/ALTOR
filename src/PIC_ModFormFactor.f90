module PIC_ModFormFactor
  use PIC_ModGrid
  implicit none
  integer,parameter::lOrderFF=2,iDownFF=-1,iUpFF=1
  real,dimension(lOrderFF,nDim)::LowFF_ID
contains
  subroutine get_form_factors(X_D,Node_D,HighFF_ID)
    use ModNumConst
    real,dimension(nDim),intent(in)::X_D
    integer,dimension(nDim),intent(out)::Node_D
    real,dimension(lOrderFF+1,nDim)::HighFF_ID
    real,dimension(nDim)::d_D
    integer::iDim
    Node_D=floor(X_D)
    d_D=X_D-real(Node_D)
    Node_D=Node_D+1
    do iDim=1,nDim
       LowFF_ID(1,iDim)=cOne-d_D(iDim)
       LowFF_ID(2,iDim)=d_D(iDim)
       HighFF_ID(1:3:2,iDim)=LowFF_ID(:,iDim)**2
       HighFF_ID(2,iDim)=cTwo-sum(HighFF_ID(1:3:2,iDim))
    end do
    LowFF_ID = LowFF_ID *cHalf
    HighFF_ID= HighFF_ID*cHalf
  end subroutine get_form_factors
end module PIC_ModFormFactor
