
!
!This file, being added to 
!PIC_main.f90 
!PIC_solver.f90 and 
!PIC_form_factor.f90, 
!solves the Vlasov-Maxwell equations for 3D geometries. 
!Dependencies: uses ModFormFactor, location:PIC_form_factor.f90

module PIC_ModMpi
  use PIC_ModField
  use ModMpi
  use PIC_ModProc
  implicit none
  !structures
  integer,parameter::iBuffSizeMax=10 !In MByte
  integer,parameter::KByte=2**10, M=KByte*KByte
  !Methods:
  public::pass_density
  public::pass_current
contains
  subroutine pass_density(iProcIn)
    integer,optional::iProcIn
    integer,parameter::iLength=max(&
       iBuffSizeMax*M/(4*(iRealPrec+1)*(nX+2*iGCN)*(nY+2*iGCN)),1)
    real,dimension(1-iGCN:nX+iGCN,1-iGCN:nY+iGCN,iLength)::Buff_G
    integer::k,kNew
    !If iProcIn is given, the result is at PE=iProcIn only
    if(present(iProcIn))then
       k=-iGCN
       do while(k<nZ+iGCN)
          kNew=min(nZ+iGCN,k+iLength)
          call MPI_REDUCE(&
               rho_G(1-iGCN,1-iGCN,k+1),&
               Buff_G(1-iGCN,1-iGCN,1),&
               (nX+2*iGCN)*(nY+2*iGCN)*(kNew-k),&
               MPI_REAL,&
               MPI_SUM,&
               iProcIn,iComm,iError)
          if(iProc==iProcIn)rho_G(:,:,k+1:kNew)=Buff_G(:,:,1:kNew-k)
          k=kNew
       end do
    else
       k=-iGCN
       do while(k<nZ+iGCN)
          kNew=min(nZ+iGCN,k+iLength)
          call MPI_ALLREDUCE(&
               rho_G(1-iGCN,1-iGCN,k+1),&
               Buff_G(1-iGCN,1-iGCN,1),&
               (nX+2*iGCN)*(nY+2*iGCN)*(kNew-k),&
               MPI_REAL,&
               MPI_SUM,&
               iComm,iError)
          rho_G(:,:,k+1:kNew)=Buff_G(:,:,1:kNew-k)
          k=kNew
       end do
    end if
  end subroutine pass_density
  !---------------------------------------------------------------!
  subroutine pass_current
    integer,parameter::iLength=max(&
       iBuffSizeMax*M/(12*(iRealPrec+1)*(nX+2*iGCN+1)*(nY+2*iGCN+1)),1)
    real,dimension(3,0-iGCN:nX+iGCN,0-iGCN:nY+iGCN,iLength)::Buff_G
    integer::k,kNew
    !If iProcIn is given, the result is at PE=iProcIn only
    k=-iGCN-1
    do while(k<nZ+iGCN)
       kNew=min(nZ+iGCN,k+iLength)
       call MPI_ALLREDUCE(&
            Counter_DG(1,0-iGCN,0-iGCN,k+1),&
            Buff_G(1,0-iGCN,0-iGCN,1),&
            (nX+2*iGCN+1)*(nY+2*iGCN+1)*(kNew-k),&
            MPI_REAL,&
            MPI_SUM,&
            iComm,iError)
       Counter_DG(:,:,:,k+1:kNew)=Buff_G(:,:,:,1:kNew-k)
       k=kNew
    end do
  end subroutine pass_current
  !--------------------------------------------------------------!
end module PIC_ModMpi
   
    
