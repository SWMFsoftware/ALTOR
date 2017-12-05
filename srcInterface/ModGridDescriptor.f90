!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!MHD grid in BATSRUS 
module PC_domain_decomposition

  use CON_grid_storage, ProcToolkit_ => Pe_

  implicit none

  SAVE

  interface PC_get_root_decomposition
     module procedure PC_get_roots_id
     module procedure PC_get_roots_dd
  end interface

  public:: PC_get_root_decomposition
  public:: PC_update_local_decomposition

  type(DomainType), public:: PC_Domain
    
  integer, parameter, private::    &
       PELast_      = 5, &
       LEV_         = 6, &
       LEVmin_      = 7, &
       LEVmax_      = 8

  ! Position of children relative to the parent block
  ! in the Morton ordering
  integer, parameter, private:: iShiftMorton_DI(3,8)= reshape( (/ &
       0,0,0, &
       1,0,0, &
       0,1,0, &
       1,1,0, &
       0,0,1, &
       1,0,1, &
       0,1,1, &
       1,1,1 /), (/3,8/))

  private:: show_domain_decomp, get_batl_tree, &
       PC_get_roots_dd, PC_get_roots_id

contains

  !===========================================================================
  subroutine show_domain_decomp(Dd)
    use PIC_ModProc, ONLY: iProc

    type(DomainType),intent(in):: Dd
    integer:: iNode, iChild
    !-------------------------------------------------------------------------
    if(iProc /= 0) RETURN

    write(*,*)'!!! Starting show_domain_decomp'
    write(*,*)'!!! CompID_            =', Dd%CompID_
    write(*,*)'!!! nDim               =', Dd%nDim
    write(*,*)'!!! XyzMin_D           =', Dd%XyzMin_D
    write(*,*)'!!! XyzMax_D           =', Dd%XyzMax_D
    write(*,*)'!!! iRootMapDim_D      =', Dd%iRootMapDim_D
    write(*,*)'!!! IsTreeDecomposition=', Dd%IsTreeDecomposition
    write(*,*)'!!! nDimTree           =', Dd%nDimTree
    write(*,*)'!!! nChildren          =', Dd%nChildren
    write(*,*)'!!! nDim               =', Dd%nDim
    write(*,*)'!!! IsTreeDecomposition=', Dd%IsTreeDecomposition
    write(*,*)'!!! nTreeNodes         =', Dd%nTreeNodes
    write(*,*)'!!! nAllocatedNodes    =', Dd%nAllocatedNodes
    write(*,*)'!!! IsPeriodic_D       =', Dd%IsPeriodic_D
    write(*,*)'!!! DoGlueMargins      =', Dd%DoGlueMargins
    write(*,*)'!!! iRealization       =', Dd%iRealization
    write(*,*)'!!! IsLocal            =', Dd%IsLocal
    write(*,*)'!!! MinBlock           =', Dd%MinBlock
    write(*,*)'!!! MaxBlock           =', Dd%MaxBlock
    write(*,*)'!!! nBlockAll          =', Dd%nBlockAll

    write(*,*)'!!! iChild, iShift_DI'
    do iChild = 1, Dd%nChildren
       write(*,*) iChild, Dd%iShift_DI(:,iChild)
    end do

    write(*,*)'!!! iNode, iDecomposition_II'
    do iNode = 1, Dd%nTreeNodes
       write(*,*) iNode, Dd%iDecomposition_II(:,iNode)
    end do

    write(*,*)'!!! iNode, XyzBlock_DI'
    do iNode = 1, Dd%nTreeNodes
       write(*,*) iNode, Dd%XyzBlock_DI(:,iNode), Dd%DXyzCell_DI(:,iNode)
    enddo

    write(*,*)'!!! Done with show_domain_decomp'

  end subroutine show_domain_decomp
  !===========================================================================
  subroutine get_batl_tree(Domain)

    ! Avoid name conflict with Parent_ in the SWMF coupling toolkit
    use PC_BATL_tree, ParentBatl_ => Parent_

    type(DomainType),intent(inout)::Domain

    integer:: iNode, iNodeParent, iChild
    !-------------------------------------------------------------------------

    ! Allocate arrays for nNode sized tree
    Domain%nTreeNodes = nNode
    call check_octree_grid_allocation(Domain)

    ! Here we assume that there are no holes in the BATL tree
    do iNode = 1, nNode

       iNodeParent = iTree_IA(ParentBatl_,iNode)
       if(iNodeParent == Unset_)then
          ! For root blocks coupling toolkit seems to set parent to itself
          Domain%iDecomposition_II(Parent_,iNode) = iNode
          ! For root blocks coupling toolkit seems to set child index to 0
          Domain%iDecomposition_II(MyNumberAsAChild_,iNode) = 0
       else
          Domain%iDecomposition_II(Parent_,iNode) = iNodeParent
          ! Find child index
          do iChild = 1, nChild
             if(iTree_IA(Child0_+iChild,iNodeParent) == iNode) then
                Domain%iDecomposition_II(MyNumberAsAChild_,iNode)&
                     =iChild
                EXIT
             end if
          end do
       end if

       if(iTree_IA(Status_,iNode) == Unused_)then 
          do iChild = 1, nChild
             ! iChildOrder_II may be required here !!!
             Domain%iDecomposition_II(iChild,iNode) = &
                  iTree_IA(Child0_+iChild,iNode) 
          end do
       else
          Domain%iDecomposition_II(FirstChild_,iNode) = &
               None_
          Domain%iDecomposition_II(GlobalBlock_,iNode) = &
               iMortonNode_A(iNode)
          Domain%iDecomposition_II(ProcToolkit_,iNode) = &
               iTree_IA(Proc_,iNode)
          Domain%iDecomposition_II(PELast_,iNode) = &
               iTree_IA(Proc_,iNode)
          Domain%iDecomposition_II(BLK_,iNode) = &
               iTree_IA(Block_,iNode)
          Domain%iDecomposition_II(LEV_,iNode) = &
               iTree_IA(Level_,iNode)
          Domain%iDecomposition_II(LEVmin_,iNode) = &
               iTree_IA(MinLevel_,iNode)
          Domain%iDecomposition_II(LEVmax_,iNode) = &
               iTree_IA(MaxLevel_,iNode)
       end if
    end do

    ! call show_domain_decomp(Domain)

  end subroutine get_batl_tree
  !===========================================================================
  subroutine PC_get_roots_dd(Domain)                         

    use PC_BATL_lib, ONLY: nIJK_D, IsPeriodic_D, nRoot_D, CoordMin_D, CoordMax_D

    type(DomainType),intent(inout)::Domain  
    !-------------------------------------------------------------------------

    call get_root_decomposition_dd(&
         Domain,       & ! Decomposition to be constructed
         nRoot_D,                   & ! As in DomainType
         CoordMin_D,                & ! As in DomainType
         CoordMax_D,                & ! As in DomainType
         nIJK_D,                    & ! As in DomainType
         IsPeriodic_D=IsPeriodic_D, &
         iShift_DI=iShiftMorton_DI)

  end subroutine PC_get_roots_dd
  !===========================================================================
  subroutine PC_get_roots_id(GridID_)                         

    use PC_BATL_lib, ONLY: nIJK_D, IsPeriodic_D, nRoot_D, CoordMin_D, CoordMax_D

    integer, intent(in):: GridID_  
    !-------------------------------------------------------------------------
    call get_root_decomposition_id(&
         GridID_,                   & ! Decomposition to be constructed
         nRoot_D,                   & ! As in DomainType
         CoordMin_D,                & ! As in DomainType
         CoordMax_D,                & ! As in DomainType
         nIJK_D,                    & ! As in DomainType
         IsPeriodic_D=IsPeriodic_D, &
         iShift_DI=iShiftMorton_DI)

  end subroutine PC_get_roots_id

  !==========================================================================
  subroutine PC_update_local_decomposition(Domain)

    type(DomainType), intent(inout):: Domain
    !-----------------------------------------------------------------------
    call get_batl_tree(Domain)

    Domain%iRealization = &
         mod(Domain%iRealization+1, 1000)
    call complete_grid(Domain)

  end subroutine PC_update_local_decomposition

end module PC_domain_decomposition
