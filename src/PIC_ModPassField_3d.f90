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
  use PIC_ModMain, ONLY: UseSharedField
  use PC_BATL_pass_face_field, ONLY: add_ghost_face_field, add_ghost_cell_field
  use PC_BATL_lib, ONLY: nBlock
  implicit none
  !structures
  integer,parameter::iBuffSizeMax=10 !In MByte
  integer,parameter,private::KByte=1024, M=KByte*KByte
  !Methods:
  public::pass_density
  public::pass_velocity
  public::pass_current
contains
  subroutine pass_density(iProcIn)
    integer,optional::iProcIn
    integer,parameter::iLength=max(&
       iBuffSizeMax*M/(4*(iRealPrec+1)*nX*nY),1)
    real,dimension(1:nX,1:nY,iLength)::Buff_G
    integer::k,kNew,iBlock
    !-----------------------------------------------------------
    !Set up periodic BC in all three dimensions.
    call add_ghost_cell_field(1,iGCN,Rho_GB)
    if(nProc==1)return
    !If iProcIn is given, the result is at PE=iProcIn only
    
    if(present(iProcIn))then
       do iBlock = 1, nBlock
          k=0
          do while(k<nZ)
             kNew=min(nZ,k+iLength)
             call MPI_REDUCE(&
                  Rho_GB(1,1,k+1,iBlock),&
                  Buff_G(1,1,1),&
                  nX*nY*(kNew-k),&
                  MPI_REAL,&
                  MPI_SUM,&
                  iProcIn,iComm,iError)
             if(iProc==iProcIn)Rho_GB(1:nX,1:nY,k+1:kNew,iBlock)=&
                  Buff_G(:,:,1:kNew-k)
             k=kNew
          end do
       end do
    else
       do iBlock = 1, nBlock
          k=0
          do while(k<nZ)
             kNew=min(nZ,k+iLength)
             call MPI_ALLREDUCE(&
                  Rho_GB(1,1,k+1,iBlock),&
                  Buff_G(1,1,1),&
                  nX*nY*(kNew-k),&
                  MPI_REAL,&
                  MPI_SUM,&
                  iComm,iError)
             Rho_GB(1:nX,1:nY,k+1:kNew,iBlock)=Buff_G(:,:,1:kNew-k)
             k=kNew
          end do
       end do
    end if
  end subroutine pass_density
  !==============================================================
  !
  subroutine pass_velocity(iProcIn)
    integer,optional::iProcIn
    integer,parameter::iLength=max(&
         iBuffSizeMax*M/(4*(iRealPrec+1)*nX*nY),1)
    real,dimension(1:MaxDim,1:nX,1:nY,iLength)::Buff_DG
    integer::k, kNew, iBlock
    !-----------------------------------------------------------
    !Set up periodic BC in all three dimensions.
    call add_ghost_cell_field(3,iGCN,V_DGB)
    if(nProc==1)return
    !If iProcIn is given, the result is at PE=iProcIn only
    if(present(iProcIn))then
       do iBlock = 1, nBlock
          k=0
          do while(k<nZ)
             kNew=min(nZ,k+iLength)
             call MPI_REDUCE(&
                  V_DGB(1,1,1,k+1,iBlock),&
                  Buff_DG(1,1,1,1),&
                  MaxDim*nX*nY*(kNew-k),&
                  MPI_REAL,&
                  MPI_SUM,&
                  iProcIn,iComm,iError)
             if(iProc==iProcIn)&
                  V_DGB(:,1:nX,1:nY,k+1:kNew,iBlock) = Buff_DG(:,:,:,1:kNew-k)
             k=kNew
          end do
       end do
    else
       do iBlock = 1, nBlock
          k=0
          do while(k<nZ)
             kNew=min(nZ,k+iLength)
             call MPI_ALLREDUCE(&
                  V_DGB(1,1,1,k+1,iBlock),&
                  Buff_DG(1,1,1,1),&
                  MaxDim*nX*nY*(kNew-k),&
                  MPI_REAL,&
                  MPI_SUM,&
                  iComm,iError)
             V_DGB(:,1:nX,1:nY,k+1:kNew,1)=Buff_DG(:,:,:,1:kNew-k)
             k=kNew
          end do
       end do
    end if
  end subroutine pass_velocity
  !=============================================================
  subroutine pass_current
    use PC_BATL_size
    real,dimension(0:nX,1-jDim_:nY,1-kDim_:nZ,MaxDim)::Buff_G
    integer::iBlock
    !-------------------
    call add_ghost_face_field(iGCN,Counter_GDB)
    if(nProc==1)RETURN
    !\
    ! if the blocks are distributed, the current through the internal
    ! faces is correctly calculated
    !/
    if(.not.UseSharedField)RETURN
    do iBlock = 1, nBlock
       call MPI_ALLREDUCE(&
            Counter_GDB(0:nX,1-jDim_:nY,1-kDim_:nZ,:,iBlock),&
            Buff_G,&
            (nX+1)*(nY+jDim_)*(nZ+kDim_)*MaxDim*nBlock,&
            MPI_REAL,&
            MPI_SUM,&
            iComm,iError)
       Counter_GDB(0:nX,1-jDim_:nY,1-kDim_:nZ,:,iBlock) = Buff_G 
    end do
  end subroutine pass_current
  !--------------------------------------------------------------!
end module PIC_ModMpi
   
    
