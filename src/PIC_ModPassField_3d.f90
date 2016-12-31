!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! This file contains the top level methods for ALTOR
!==========================
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
  public::pass_moments
  public::pass_current
contains
  subroutine pass_density(iSort)
    integer,intent(in) :: iSort
    integer,parameter::iProcIn = 0
    integer,parameter::iLength=max(&
       iBuffSizeMax*M/(4*(iRealPrec+1)*nX*nY),1)
    real,dimension(1,1:nX,1:nY,iLength)::Buff_G
    integer::k,kNew,iBlock
    !-----------------------------------------------------------
    !Set up periodic BC in all three dimensions.
    call add_ghost_cell_field(1,iGCN,State_VGBI(1:1,:,:,:,:,iSort))
    if(nProc==1)return
    !the result is at PE=iProcIn only
    
    do iBlock = 1, nBlock
       k=0
       do while(k<nZ)
          kNew=min(nZ,k+iLength)
          call MPI_REDUCE(&
               State_VGBI(1:1,1:nX,1:nY,k+1:kNew,iBlock,iSort),&
               Buff_G(1:1,1:nX,1:nY,1:kNew-k),&
               nX*nY*(kNew-k),&
               MPI_REAL,&
               MPI_SUM,&
               iProcIn,iComm,iError)
          if(iProc==iProcIn)State_VGBI(1,1:nX,1:nY,k+1:kNew,iBlock,iSort)=&
               Buff_G(1,:,:,1:kNew-k)
          k=kNew
       end do
    end do
  end subroutine pass_density
  !==============================================================
  !
  subroutine pass_moments(iSort)
    integer,intent(in)::iSort
    integer,parameter::iProcIn = 0
    integer,parameter::iLength=max(&
         iBuffSizeMax*M/(40*(iRealPrec+1)*nX*nY),1)
    real,dimension(1:10,1:nX,1:nY,iLength)::Buff_DG
    integer::k, kNew, iBlock
    !-----------------------------------------------------------
    !Set up periodic BC in all three dimensions.
    call add_ghost_cell_field(10,iGCN,State_VGBI(:,:,:,:,:,iSort))
    if(nProc==1)return
    if(.not.UseSharedField)RETURN
    !the result is at PE=iProcIn only
    do iBlock = 1, nBlock
       k=0
       do while(k<nZ)
          kNew=min(nZ,k+iLength)
          call MPI_REDUCE(&
               State_VGBI(1:10,1:nX,1:nY,k+1:kNew,iBlock,iSort),&
               Buff_DG(1:10,1:nX,1:nY,1:kNew-k),&
               10*nX*nY*(kNew-k),&
               MPI_REAL,&
               MPI_SUM,&
               iProcIn,iComm,iError)
          if(iProc==iProcIn)&
               State_VGBI(:,1:nX,1:nY,k+1:kNew,iBlock,iSort) = &
               Buff_DG(:,:,:,1:kNew-k)
          k=kNew
       end do
    end do
  end subroutine pass_moments
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
            (nX+1)*(nY+jDim_)*(nZ+kDim_)*MaxDim,&
            MPI_REAL,&
            MPI_SUM,&
            iComm,iError)
       Counter_GDB(0:nX,1-jDim_:nY,1-kDim_:nZ,:,iBlock) = Buff_G 
    end do
  end subroutine pass_current
  !--------------------------------------------------------------!
end module PIC_ModMpi
   
    
