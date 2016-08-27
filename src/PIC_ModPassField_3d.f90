
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
  integer,parameter,private::KByte=1024, M=KByte*KByte
  !Methods:
  public::pass_density, pass_velocity
  public::pass_current
contains
  subroutine pass_density(iProcIn)
    integer,optional::iProcIn
    integer,parameter::iLength=max(&
       iBuffSizeMax*M/(4*(iRealPrec+1)*(nX+2*iGCN)*(nY+2*iGCN)),1)
    real,dimension(1-iGCN:nX+iGCN,1-iGCN:nY+iGCN,iLength)::Buff_G
    integer::k,kNew
    !-----------------------------------------------------------
    !Set up periodic BC in all three dimensions.
    call density_bc_periodic(1)
    if(nProc==1)return
    !If iProcIn is given, the result is at PE=iProcIn only
    if(present(iProcIn))then
       k=-iGCN
       do while(k<nZ+iGCN)
          kNew=min(nZ+iGCN,k+iLength)
          call MPI_REDUCE(&
               Rho_GB(1-iGCN,1-iGCN,k+1,1),&
               Buff_G(1-iGCN,1-iGCN,1),&
               (nX+2*iGCN)*(nY+2*iGCN)*(kNew-k),&
               MPI_REAL,&
               MPI_SUM,&
               iProcIn,iComm,iError)
          if(iProc==iProcIn)Rho_GB(:,:,k+1:kNew,1)=Buff_G(:,:,1:kNew-k)
          k=kNew
       end do
    else
       k=-iGCN
       do while(k<nZ+iGCN)
          kNew=min(nZ+iGCN,k+iLength)
          call MPI_ALLREDUCE(&
               Rho_GB(1-iGCN,1-iGCN,k+1,1),&
               Buff_G(1-iGCN,1-iGCN,1),&
               (nX+2*iGCN)*(nY+2*iGCN)*(kNew-k),&
               MPI_REAL,&
               MPI_SUM,&
               iComm,iError)
          Rho_GB(:,:,k+1:kNew,1)=Buff_G(:,:,1:kNew-k)
          k=kNew
       end do
    end if
  end subroutine pass_density
  !==============================================================
  !Test
  subroutine pass_velocity(iProcIn)
    integer,optional::iProcIn
    integer,parameter::iLength=max(&
         iBuffSizeMax*M/(4*(iRealPrec+1)*(nX+2*iGCN)*(nY+2*iGCN)),1)
    real,dimension(1:3,1-iGCN:nX+iGCN,1-iGCN:nY+iGCN,iLength)::Buff_G
    integer::k,kNew
    !-----------------------------------------------------------
    !Set up periodic BC in all three dimensions.
    call velocity_bc_periodic(1)
    if(nProc==1)return
    !If iProcIn is given, the result is at PE=iProcIn only
    if(present(iProcIn))then
       k=-iGCN
       do while(k<nZ+iGCN)
          kNew=min(nZ+iGCN,k+iLength)
          call MPI_REDUCE(&
               V_GDB(:,1-iGCN,1-iGCN,k+1,1),&
               Buff_G(:,1-iGCN,1-iGCN,1),&
               3*(nX+2*iGCN)*(nY+2*iGCN)*(kNew-k),&
               MPI_REAL,&
               MPI_SUM,&
               iProcIn,iComm,iError)
          if(iProc==iProcIn)&
             V_GDB(:,:,:,k+1:kNew,1) = Buff_G(:,:,:,1:kNew-k)
          k=kNew
       end do
    else
       k=-iGCN
       do while(k<nZ+iGCN)
          kNew=min(nZ+iGCN,k+iLength)
          call MPI_ALLREDUCE(&
               V_GDB(:,1-iGCN,1-iGCN,k+1,1),&
               Buff_G(:,1-iGCN,1-iGCN,1),&
               3*(nX+2*iGCN)*(nY+2*iGCN)*(kNew-k),&
               MPI_REAL,&
               MPI_SUM,&
               iComm,iError)
          V_GDB(:,:,:,k+1:kNew,1)=Buff_G(:,:,:,1:kNew-k)
          k=kNew
       end do
    end if
  end subroutine pass_velocity
  !=============================================================
  subroutine pass_current
    integer,parameter::iLength=max(&
       iBuffSizeMax*M/(4*(iRealPrec+1)*(nX+2*iGCN)*(nY+2*iGCN)),1)
!    real,dimension(1-iGCN:nX+iGCN,1-iGCN:nY+iGCN,iLength)::Buff_G
    real,dimension(0:nX,0:nY,0:nZ,3)::Buff_G
    integer::k,kNew,iDim
    !-------------------
    call current_bc_periodic(1)
    if(nProc==1)return
          call MPI_ALLREDUCE(&
               Counter_GDB(0:nX,0:nY,0:nZ,:,1),&
               Buff_G,&
               (nX+1)*(nY+1)*(nZ+1)*3,&
               MPI_REAL,&
               MPI_SUM,&
               iComm,iError)
          Counter_GDB(0:nX,0:nY,0:nZ,:,1) = Buff_G 
  end subroutine pass_current
  !--------------------------------------------------------------!
end module PIC_ModMpi
   
    
