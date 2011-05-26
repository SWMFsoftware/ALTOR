!This file, once added to 
!PIC_main.f90 
!PIC_solver.f90 and 
!PIC_form_factor.f90, 
!solves the Vlasov-Maxwell equations for 2D geometries. 
!Dependencies: uses ModFormFactor, location:PIC_form_factor.f90

module PIC_ModMpi
  use PIC_ModField
  use ModMpi
  use PIC_ModProc
  implicit none
  !structures
  integer,parameter::iBuffSizeMax=10 !In MByte
  integer,parameter,private::KByte=1024, M=KByte*KByte
  !Methods:
  public::pass_density
  public::pass_current
contains
  subroutine pass_density(iProcIn)
    integer,optional::iProcIn
    integer,parameter::iLength=max(&
       iBuffSizeMax*M/(4*(iRealPrec+1)*(nX+2*iGCN)),1)
    real,dimension(1-iGCN:nX+iGCN,iLength)::Buff_G
    integer::j,jNew
    !----------------

    !If iProcIn is given, the result is at PE=iProcIn only

    if(present(iProcIn))then
       j = -iGCN
       do while(j<nY+iGCN)
          jNew=min(nY+iGCN,j+iLength)
          call MPI_REDUCE(&
               rho_G(1-iGCN,j+1),&
               Buff_G(1-iGCN,1),&
               (nX+2*iGCN)*(jNew-j),&
               MPI_REAL,&
               MPI_SUM,&
               iProcIn,iComm,iError)
          if(iProc==iProcIn)rho_G(:,j+1:jNew)=Buff_G(:,1:jNew-j)
          j=jNew
       end do
    else
       j=-iGCN
       do while(j<nY+iGCN)
          jNew=min(nY+iGCN,j+iLength)
          call MPI_ALLREDUCE(&
               rho_G(1-iGCN,j+1),&
               Buff_G(1-iGCN,1),&
               (nX+2*iGCN)*(jNew-j),&
               MPI_REAL,&
               MPI_SUM,&
               iComm,iError)
          rho_G(:,j+1:jNew)=Buff_G(:,1:jNew-j)
          j=jNew
       end do
    end if
  end subroutine pass_density
  !---------------------------------------------------------------!
  subroutine pass_current
    integer,parameter::iLength=max(&
       iBuffSizeMax*M/(4*(iRealPrec+1)*(nX+2*iGCN)),1)
    real,dimension(1-iGCN:nX+iGCN,iLength)::Buff_G
    integer::j,jNew,iDim
    !---------------
    do iDim = 1,3
       j= - iGCN
       do while(j<nY+iGCN)
          jNew=min(nY+iGCN,j+iLength)
          call MPI_ALLREDUCE(&
               Counter_GD(1-iGCN,j+1,iDim),&
               Buff_G(1-iGCN,1),&
               (nX+2*iGCN)*(jNew-j),&
               MPI_REAL,&
               MPI_SUM,&
               iComm,iError)
          Counter_GD(:,j+1:jNew,iDim)=Buff_G(:,1:jNew-j)
          j=jNew
       end do
    end do
  end subroutine pass_current
  !--------------------------------------------------------------!
end module PIC_ModMpi
   
    
