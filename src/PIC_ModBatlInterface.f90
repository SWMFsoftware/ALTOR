!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PIC_ModBatlInterface
  use PIC_ModMain, ONLY: MaxBlock
  implicit none
  integer, dimension(1:6,MaxBlock):: neiLEV
  integer, dimension(4,1:6,MaxBlock) :: neiPE, neiBLK
contains
  !===========================================================================
  subroutine set_altor_grid

    use PC_BATL_lib, ONLY: nBlock, Unused_B, Unused_BP, iProc, iComm, &
         IsNewDecomposition, IsNewTree, nLevelMin, nLevelMax, &
         DomainSize_D, nRoot_D, nI, CellSizeRoot

    use PIC_ModMain, ONLY: nBlockMax
    use ModMpi

    integer:: iBlock, iError
    character(len=*), parameter:: NameSub = 'set_altor_grid'
    !-------------------------------------------------------------------------

    call MPI_allreduce(nBlock, nBlockMax, 1, MPI_INTEGER, MPI_MAX, &
         iComm, iError)
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       call set_altor_block(iBlock)
    end do
  end subroutine set_altor_grid
  !===========================================================================
  subroutine set_altor_block(iBlock)

    use PC_BATL_lib, ONLY: nDim, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         Xyz_DGB, CellSize_DB, CoordMin_DB, &
         iNode_B, iNodeNei_IIIB, DiLevelNei_IIIB, &
         iTree_IA, Block_, Proc_, Unset_

    integer, intent(in):: iBlock

    ! Convert from BATL to BATSRUS ordering of subfaces. 

    integer, parameter:: iOrder_I(4) = (/1,3,2,4/)
    integer:: iNodeNei, iNodeNei_I(4)
    integer:: i, j, k
    !-------------------------------------------------------------------------
    neiLEV(1,iBlock)  = DiLevelNei_IIIB(-1,0,0,iBlock)
    neiLEV(2,iBlock)  = DiLevelNei_IIIB(+1,0,0,iBlock)
    neiLEV(3,iBlock)  = DiLevelNei_IIIB(0,-1,0,iBlock)
    neiLEV(4,iBlock)  = DiLevelNei_IIIB(0,+1,0,iBlock)
    neiLEV(5,iBlock)  = DiLevelNei_IIIB(0,0,-1,iBlock)
    neiLEV(6,iBlock)  = DiLevelNei_IIIB(0,0,+1,iBlock)

   
    select case(DiLevelNei_IIIB(-1,0,0,iBlock))
    case(Unset_)
       neiBLK(:,1,iBlock)  = Unset_
       neiPE(:,1,iBlock)   = Unset_
    case(-1)
       iNodeNei_I = pack(iNodeNei_IIIB(0,1:2,1:2,iBlock),.true.)
       iNodeNei_I = iNodeNei_I(iOrder_I)
       if(nDim < 3) where(iNodeNei_I == Unset_) iNodeNei_I = iNode_B(iBlock)
       neiBLK(:,1,iBlock)  = iTree_IA(Block_,iNodeNei_I)
       neiPE(:,1,iBlock)   = iTree_IA(Proc_,iNodeNei_I)
    case default
       iNodeNei = iNodeNei_IIIB(0,1,1,iBlock)
       neiBLK(:,1,iBlock)  = iTree_IA(Block_,iNodeNei)
       neiPE(:,1,iBlock)   = iTree_IA(Proc_,iNodeNei)
    end select

    select case(DiLevelNei_IIIB(+1,0,0,iBlock))
    case(Unset_)
       neiBLK(:,2,iBlock)  = Unset_
       neiPE(:,2,iBlock)   = Unset_
    case(-1)
       iNodeNei_I = pack(iNodeNei_IIIB(3,1:2,1:2,iBlock),.true.)
       iNodeNei_I = iNodeNei_I(iOrder_I)
       if(nDim < 3) where(iNodeNei_I == Unset_) iNodeNei_I = iNode_B(iBlock)
       neiBLK(:,2,iBlock)  = iTree_IA(Block_,iNodeNei_I)
       neiPE(:,2,iBlock)   = iTree_IA(Proc_,iNodeNei_I)
    case default
       iNodeNei = iNodeNei_IIIB(3,1,1,iBlock)
       neiBLK(:,2,iBlock)  = iTree_IA(Block_,iNodeNei)
       neiPE(:,2,iBlock)   = iTree_IA(Proc_,iNodeNei)
    end select

    select case(DiLevelNei_IIIB(0,-1,0,iBlock))
    case(Unset_)
       neiBLK(:,3,iBlock)  = Unset_
       neiPE(:,3,iBlock)   = Unset_
    case(-1)
       iNodeNei_I = pack(iNodeNei_IIIB(1:2,0,1:2,iBlock),.true.)
       iNodeNei_I = iNodeNei_I(iOrder_I)
       if(nDim < 3) where(iNodeNei_I == Unset_) iNodeNei_I = iNode_B(iBlock)
       neiBLK(:,3,iBlock)  = iTree_IA(Block_,iNodeNei_I)
        neiPE(:,3,iBlock)  = iTree_IA(Proc_,iNodeNei_I)
    case default
       iNodeNei = iNodeNei_IIIB(1,0,1,iBlock)
       neiBLK(:,3,iBlock)  = iTree_IA(Block_,iNodeNei)
       neiPE(:,3,iBlock)   = iTree_IA(Proc_,iNodeNei)
    end select

    select case(DiLevelNei_IIIB(0,+1,0,iBlock))
    case(Unset_)
       neiBLK(:,4,iBlock)  = Unset_
       neiPE(:,4,iBlock)   = Unset_
    case(-1)
       iNodeNei_I = pack(iNodeNei_IIIB(1:2,3,1:2,iBlock),.true.)
       iNodeNei_I = iNodeNei_I(iOrder_I)
       if(nDim < 3) where(iNodeNei_I == Unset_) iNodeNei_I = iNode_B(iBlock)
       neiBLK(:,4,iBlock)  = iTree_IA(Block_,iNodeNei_I)
       neiPE(:,4,iBlock)   = iTree_IA(Proc_,iNodeNei_I)
    case default
       iNodeNei = iNodeNei_IIIB(1,3,1,iBlock)
       neiBLK(:,4,iBlock)  = iTree_IA(Block_,iNodeNei)
       neiPE(:,4,iBlock)   = iTree_IA(Proc_,iNodeNei)
    end select

    select case(DiLevelNei_IIIB(0,0,-1,iBlock))
    case(Unset_ )
       neiBLK(:,5,iBlock)  = Unset_
       neiPE(:,5,iBlock)   = Unset_
    case(-1)
       iNodeNei_I = pack(iNodeNei_IIIB(1:2,1:2,0,iBlock),.true.)
       neiBLK(:,5,iBlock)  = iTree_IA(Block_,iNodeNei_I)
       neiPE(:,5,iBlock)   = iTree_IA(Proc_,iNodeNei_I)
    case default
       iNodeNei = iNodeNei_IIIB(1,1,0,iBlock)
       neiBLK(:,5,iBlock)  = iTree_IA(Block_,iNodeNei)
       neiPE(:,5,iBlock)   = iTree_IA(Proc_,iNodeNei)
    end select

    select case(DiLevelNei_IIIB(0,0,+1,iBlock))
    case(Unset_)
       neiBLK(:,6,iBlock)  = Unset_
       neiPE(:,6,iBlock)   = Unset_
    case(-1)
       iNodeNei_I = pack(iNodeNei_IIIB(1:2,1:2,3,iBlock),.true.)
       neiBLK(:,6,iBlock)  = iTree_IA(Block_,iNodeNei_I)
       neiPE(:,6,iBlock)   = iTree_IA(Proc_,iNodeNei_I)
    case default
       iNodeNei = iNodeNei_IIIB(1,1,3,iBlock)
       neiBLK(:,6,iBlock)  = iTree_IA(Block_,iNodeNei)
       neiPE(:,6,iBlock)   = iTree_IA(Proc_,iNodeNei)
    end select
  end subroutine set_altor_block
  !=============================
end module PIC_ModBatlInterface
