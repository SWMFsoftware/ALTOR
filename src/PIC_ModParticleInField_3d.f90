module PIC_ModParticleInField
  use PIC_ModSize, ONLY:nBlock, nX, nY, nZ
  use PIC_ModField
  use PIC_ModFormFactor, ONLY: lOrderFF, iUpFF, iDownFF, iExt
  use PIC_ModFormfactor, ONLY: Node_D, NodeNew_D
  use PIC_ModFormFactor, ONLY: HighFF_ID, HighFFNew_ID, LowFF_ID
  implicit none
  SAVE
  
  !As long as Node_D and form-factors are defined, we
  !can get interpolated values from the structures of the
  !ModField class (electric and magnetic field) or calculate
  !the input from a given particle to the structures of the 
  !ModField class (density, electric current)
  !Methods
  public::get_b_from_a
  public::e_interpolated_d !Interpolated electric field
  public::b_interpolated_d !Interpolated magnetic field
  public::add_density  !Adds an input to cell-centered density
  public::add_moments  !Adds an input to cell-centered density
                           !and velocity moments from a given particle
  public::add_current !Adds an input to an electric current
contains
  !================================
  function e_interpolated_d(iBlock)
    integer, intent(in)::iBlock
    !Sets the pointer to the fragment of E_GDB array to interpolate E
    real,dimension(1:3)::e_interpolated_d
    integer::i, j, k, n_D(3), n1_D(3)
    !----------------

    e_interpolated_d = 0.0
    n_D = Node_D + iDownFF; n1_D = n_D -1
    
    do k=1, lOrderFF + 1
       do j=1, lOrderFF + 1
          e_interpolated_d(x_)= e_interpolated_d(x_)+ &
               sum(E_GDB(n_D(1)-iExt:lOrderFF+iExt+n1_D(1), &
               j+n1_D(2), k+n1_d(3), x_,iBlock)*LowFF_ID(:, x_))&
               *HighFF_ID(j, y_)*HighFF_ID(k, z_)
       end do
       do j=1-iExt,lOrderFF+iExt
          e_interpolated_d(y_)= e_interpolated_d(y_)+ &
               sum(E_GDB(n_D(1):lOrderFF+n_D(1), &
               j+n1_D(2), k+n1_D(3), y_,iBlock)*HighFF_ID(:,x_))&
               *LowFF_ID(j,y_)*HighFF_ID(k,z_)
       end do
    end do
    do k=1-iExt,lOrderFF+iExt
       do j=1,lOrderFF+1
          e_interpolated_d(z_)= e_interpolated_d(z_)+ &
               sum(E_GDB(n_D(1):lOrderFF+n_D(1),&
               j+n1_D(2),k+n1_D(3),z_,iBlock)*HighFF_ID(:,x_))&
               *HighFF_ID(j,y_)*LowFF_ID(k,z_)
       end do
    end do
  end function e_interpolated_d

  !===========================
  function b_interpolated_d(iBlock)
    integer,intent(in)::iBlock
    !Sets the pointer to the fragment of the magnetic field array to 
    !interpolate B
    real,dimension(1:3)::b_interpolated_d
 
    integer::i, j, k, n_D(3), n1_D(3)
    !--------------
    b_interpolated_d = 0.0

    n_D = Node_D + iDownFF; n1_D = n_D -1
    do k=1-iExt,lOrderFF+iExt
       do j=1-iExt,lOrderFF+iExt
          b_interpolated_d(x_)= b_interpolated_d(x_)+ &
               sum(B_GDB(n_D(1):lOrderFF+n_D(1),&
               j+n1_D(2),k+n1_D(3),x_,iBlock)*HighFF_ID(:,x_))&
               *LowFF_ID(j,y_)*LowFF_ID(k,z_)
       end do
       do j=1,lOrderFF+1
          b_interpolated_d(y_)= b_interpolated_d(y_)+ &
               sum(B_GDB(n_D(1)-iExt:lOrderFF+iExt+n1_D(1),&
               j+n1_D(2),k+n1_D(3),y_,iBlock)*LowFF_ID(:,x_))&
               *HighFF_ID(j,y_)*LowFF_ID(k,z_)
       end do
    end do
    do k=1,lOrderFF+1
       do j=1-iExt,lOrderFF+iExt
          b_interpolated_d(z_)= b_interpolated_d(z_)+ &
               sum(B_GDB(n_D(1)-iExt:n1_D(1)+lOrderFF+iExt,&
               j+n1_D(2),k+n1_D(3),z_,iBlock)*LowFF_ID(:,x_))&
               *LowFF_ID(j,y_)*HighFF_ID(k,z_)
       end do
    end do
 
  end function b_interpolated_d
  !==============================================================
  subroutine add_moments(W_D,NodeIn_D,HighFFIn_ID,iBlock,iSort)
    !Adds an input to the density and velocity, from a given particle
    !using the same form factors
    real,   dimension(MaxDim),intent(in) :: W_D  !Velocity
    integer,dimension(MaxDim),intent(in) :: NodeIn_D
    real,   dimension(1+lOrderFF,MaxDim),intent(in)::HighFFIn_ID
    integer,intent(in):: iBlock,iSort
    integer:: i,j,k,i1,j1,k1
    real:: FFProduct, Var_V(10)
    !-----------------------   
    Var_V = (/1.0, W_D, W_D**2, W_D(1)*W_D(2), W_D(1)*W_D(3), W_D(2)*W_D(3)/)
    k1=0
    do k=NodeIn_D(3)+iDownFF*kDim_,NodeIn_D(3)+iUpFF*kDim_
       k1=k1+1
       j1=0
       do j=NodeIn_D(2)+iDownFF*jDim_,NodeIn_D(2)+iUpFF*jDim_
          j1=j1+1
          i1=0
          FFProduct=HighFFIn_ID(k1,z_)*HighFFIn_ID(j1,y_)
          do i=NodeIn_D(1)+iDownFF,NodeIn_D(1)+iUpFF
             i1=i1+1
             !Add the product of formfactors
             State_VGBI(:,i,j,k,iBlock,iSort) = State_VGBI(:,i,j,k,iBlock,iSort) + &
                  HighFFIn_ID(i1,x_)*FFProduct*Var_V  
          end do
       end do
    end do

  end subroutine add_moments
  !=========================
  subroutine add_density(NodeIn_D,HighFFIn_ID,iBlock,iSort)
    !Adds an input to the density only.
    integer,dimension(MaxDim),intent(in) :: NodeIn_D
    real,   dimension(1+lOrderFF,MaxDim),intent(in)::HighFFIn_ID
    integer,intent(in):: iBlock,iSort
    integer:: i,j,k,i1,j1,k1
    real:: FFProduct
    !-----------------------   
    k1=0
    do k=NodeIn_D(3)+iDownFF*kDim_,NodeIn_D(3)+iUpFF*kDim_
       k1=k1+1
       j1=0
       do j=NodeIn_D(2)+iDownFF*jDim_,NodeIn_D(2)+iUpFF*jDim_
          j1=j1+1
          i1=0
          FFProduct=HighFFIn_ID(k1,z_)*HighFFIn_ID(j1,y_)
          do i=NodeIn_D(1)+iDownFF,NodeIn_D(1)+iUpFF
             i1=i1+1
             !Add the product of formfactors
             State_VGBI(1,i,j,k,iBlock,iSort) = State_VGBI(1,i,j,k,iBlock,iSort) + &
                  HighFFIn_ID(i1,x_)*FFProduct
          end do
       end do
    end do

  end subroutine add_density
  !============================================
  subroutine add_current(QPerVDx_D,W_D,iBlock) 
    real,intent(in)::QPerVDx_D(x_:z_),W_D(x_:z_)
    integer,intent(in)::iBlock
    optional::W_D
    logical:: IsExtended
    !\
    ! For extended stencil
    !/
    real :: CurrentFragmentExt_GD(0:iUpFF-iDownFF+2,&
         0:iUpFF-iDownFF+2,&
         0:iUpFF-iDownFF+2,1:3)
    integer :: iD_D(3), iU_D(3)
    real,parameter::sqrt13=1/1.732050808

    real,dimension(1:lOrderFF+1,nDim)::DeltaFPlus,DeltaFMinus
    real,dimension(1:lOrderFF  ,nDim)::FI
    integer::i, j, k, iDim, n_D(3), n1_D(3)
    !-------------------
    IsExtended = any(Node_D/=NodeNew_D)

    if(IsExtended)then
       call get_current_extended

       iD_D = min(NodeNew_D - Node_D,0)
       iU_D = max(NodeNew_D - Node_D,0)

       Current_GDB(&
            Node_D(1)+iDownFF+iD_D(1):Node_D(1)+iUpFF+iU_D(1),&
            Node_D(2)+iDownFF+iD_D(2):Node_D(2)+iUpFF+iU_D(2),&
            Node_D(3)+iDownFF+iD_D(3):Node_D(3)+iUpFF+iU_D(3),&
            1:3,iBlock) = &
            Current_GDB( Node_D(1)+iDownFF+iD_D(1):Node_D(1)+iUpFF+iU_D(1),&
            Node_D(2)+iDownFF+iD_D(2):Node_D(2)+iUpFF+iU_D(2),&
            Node_D(3)+iDownFF+iD_D(3):Node_D(3)+iUpFF+iU_D(3),&
            1:3,iBlock) + & 
            CurrentFragmentExt_GD(      1+iD_D(1):iUpFF-iDownFF+1+iU_D(1),&
            1+iD_D(2):iUpFF-iDownFF+1+iU_D(2),&
            1+iD_D(3):iUpFF-iDownFF+1+iU_D(3),&
            1:3)                
    else
       DeltaFPlus  = HighFFNew_ID + HighFF_ID
       DeltaFMinus = HighFFNew_ID - HighFF_ID

       FI(1,:)=DeltaFMinus(1,:)

       !Use recursive formula for FI
       do i=2,lOrderFF
          FI(i,:) = FI(i-1,:) + DeltaFMinus(i,:)
       end do

       !To optimize algebra, multiply FI by -q*dt/4
       !and divide \DeltaFMinus by sqrt(3):
       do iDim=x_,z_
          FI(:,iDim)=(-QPerVDx_D(iDim)*0.250)*FI(:,iDim)
       end do

       DeltaFMinus=sqrt13*DeltaFMinus

       n_D = Node_D + iDownFF; n1_d = n_D -1


       !add current density:
       do k=1,lOrderFF+1
          do j=1,lOrderFF+1
             !\
             !  x_ currents
             !/
             Current_GDB(n_D(1):lOrderFF+n1_D(1),&
                  j+n1_D(2),k+n1_D(3),x_,iBlock) = &
                  Current_GDB(n_D(1):lOrderFF+n1_D(1),&
                  j+n1_D(2),k+n1_D(3),x_,iBlock) +&
                  FI(:,x_)*(&
                  DeltaFPlus (j,y_)*DeltaFPlus (k,z_)+&
                  DeltaFMinus(j,y_)*DeltaFMinus(k,z_))
          end do
          do j=1,lOrderFF
             !\
             !  y_ currents
             !/
             Current_GDB(n_D(1):lOrderFF+n_D(1),&
                  j+n1_D(2),k+n1_D(3),y_,iBlock) = &
                  Current_GDB(n_D(1):lOrderFF+n_D(1),&
                  j+n1_D(2),k+n1_D(3),y_,iBlock) + &
                  FI(j,y_)*(&
                  DeltaFPlus (:,x_)*DeltaFPlus (k,z_)+&
                  DeltaFMinus(:,x_)*DeltaFMinus(k,z_))
          end do
       end do
       do k=1,lOrderFF
          do j=1,lOrderFF+1
             !\
             !  z_ currents
             !/
             Current_GDB(n_D(1):lOrderFF+n_D(1),&
                  j+n1_D(2),k+n1_D(3),z_,iBlock) = &
                  Current_GDB(n_D(1):lOrderFF+n_D(1),&
                  j+n1_D(2),k+n1_D(3),z_,iBlock) + &
                  FI(k,z_)*(&
                  DeltaFPlus (:,x_)*DeltaFPlus (j,y_)+&
                  DeltaFMinus(:,x_)*DeltaFMinus(j,y_))
          end do
       end do
    end if

  contains
    subroutine get_current_extended
      real,dimension(0:iUpFF-iDownFF+2,nDim)::DeltaFPlusExt,DeltaFMinusExt
      real,dimension(0:iUpFF-iDownFF+1,nDim)::FIExt
      !---------------
      CurrentFragmentExt_GD = 0.0
      DeltaFPlusExt            = 0.0

      !\     
      !First, add the properly shifted new formfactor
      !/
      DeltaFPlusExt(1+NodeNew_D(1)-Node_D(1):&
           1+NodeNew_D(1)-Node_D(1)+iUpFF-iDownFF,x_)=&
           HighFFNew_ID(:,x_)

      DeltaFPlusExt(1+NodeNew_D(2)-Node_D(2):&
           1+NodeNew_D(2)-Node_D(2)+iUpFF-iDownFF,y_) = &
           HighFFNew_ID(:,y_)

      DeltaFPlusExt(1+NodeNew_D(3)-Node_D(3):&
           1+NodeNew_D(3)-Node_D(3)+iUpFF-iDownFF,z_) = &
           HighFFNew_ID(:,z_)

      DeltaFMinusExt = DeltaFPlusExt

      !\
      !Then add and subtruct the unshifted old formfactor 
      !/          
      DeltaFPlusExt(1:iUpFF-iDownFF+1,:)  = &
           DeltaFPlusExt( 1:iUpFF-iDownFF+1,:) + HighFF_ID

      DeltaFMinusExt(1:iUpFF-iDownFF+1,:) = &
           DeltaFMinusExt(1:iUpFF-iDownFF+1,:) - HighFF_ID

      !\
      ! Using recurrent relationships, find F_I
      !/
      FIExt(0,:)=DeltaFMinusExt(0,:)

      do i=1,iUpFF-iDownFF+1
         FIExt(i,:) = FIExt(i-1,:) + DeltaFMinusExt(i,:)
      end do

      !To optimize algebra, multiply FI by -q*4*\pi*dt/4
      !and divide \DeltaFMinus by sqrt(3):
      do iDim=x_,z_
         FIExt(:,iDim) = (-QPerVDx_D(iDim)*0.250)*FIExt(:,iDim)
      end do
      DeltaFMinusExt = sqrt13*DeltaFMinusExt

      !add current density:
      do k=0,iUpFF-iDownFF+2
         do j=0,iUpFF-iDownFF+2
            !\
            !  x_ currents
            !/
            CurrentFragmentExt_GD(0:iUpFF-iDownFF+1,j,k,x_) = &
                 FIExt(:,x_)*(&
                 DeltaFPlusExt( j,y_)*DeltaFPlusExt( k,z_)+&
                 DeltaFMinusExt(j,y_)*DeltaFMinusExt(k,z_))
         end do
         do j=0,iUpFF-iDownFF+1
            !\
            !  y_ currents
            !/
            CurrentFragmentExt_GD(0:iUpFF-iDownFF+2,j,k,y_) = &
                 FIExt(j,y_)*(&
                 DeltaFPlusExt( :,x_)*DeltaFPlusExt( k,z_)+&
                 DeltaFMinusExt(:,x_)*DeltaFMinusExt(k,z_))
         end do
      end do
      do k=0,iUpFF-iDownFF+1
         do j=0,iUpFF-iDownFF+2
            !\
            !  z_ currents
            !/
            CurrentFragmentExt_GD(0:iUpFF-iDownFF+2,j,k,z_) = &
                 FIExt(k,z_)*(&
                 DeltaFPlusExt( :,x_)*DeltaFPlusExt( j,y_)+&
                 DeltaFMinusExt(:,x_)*DeltaFMinusExt(j,y_))
         end do
      end do
    end subroutine get_current_extended
  end subroutine add_current
  !=========================
end module PIC_ModParticleInField

