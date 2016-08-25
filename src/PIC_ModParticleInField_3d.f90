module PIC_ModParticleInField
  use PIC_ModField, get_b_from_a_global=>get_b_from_a
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
  interface get_b_from_a
     module procedure get_b_from_a_global
     module procedure get_b_from_a_local
  end interface
  public::e_interpolated_d !Interpolated electric field
  public::b_interpolated_d !Interpolated magnetic field
  public::add_density      !Adds an input to a density, 
                           !from a given particle
  public::add_DensityVelocity  !Adds an input to cell-centered velocity
  public::add_current !Adds an input to an electric current
contains
  !=============================================================
  subroutine get_b_from_a_local(iBlock,BOut_GD)  
    integer,intent(in)::iBlock
    real,dimension(1:lOrderFF+1,1:lOrderFF+1,1:lOrderFF+1,1:3),&
           intent(out)::BOut_GD
    integer::i,j,k,i1,j1,k1
    character(LEN=*),parameter:: NameSub = 'PIC_get_b_from_a_local'
    !--------------------
    !Applied only if UseVectorPotential
    !The vector potential is in Magnetic_GD
    !Calculates the magnetic field.
    !
    !The resulting magnetic field is saved to BOut_GD
    !Only those values of the field are calculated which are needed
    !to interpolated the magnetic field acting on the particle
    !around Node_D
    call CON_stop(NameSub//' is not tested')
    k1=0
    do k=Node_D(3)+iDownFF,Node_D(3)+iUpFF
       k1=k1+1
       j1=0
       do j=Node_D(2)+iDownFF,Node_D(2)+iUpFF
          j1=j1+1
          i1=0
          do i=Node_D(1)+iDownFF,Node_D(1)+iUpFF
             i1=i1+1
             BOut_GD(i1,j1,k1,x_) = &
                  SpeedOfLight_D(y_)*&
                  (B_GDB(i,j+1,k,z_,iBlock) - B_GDB(i,j,k,z_,iBlock)) &
                  - SpeedOfLight_D(z_)*&
                  (B_GDB(i,j,k+1,y_,iBlock) - B_GDB(i,j,k,y_,iBlock))
             
             BOut_GD(i1,j1,k1,y_) = &
                  SpeedOfLight_D(z_)*&
                  (B_GDB(i,j,k+1,x_,iBlock) - B_GDB(i,j,k,x_,iBlock)) &
                  - SpeedOfLight_D(x_)*&
                  (B_GDB(i+1,j,k,z_,iBlock) - B_GDB(i,j,k,z_,iBlock))
             
             BOut_GD(i1,j1,k1,z_) = &
                  SpeedOfLight_D(x_)*&
                  (B_GDB(i+1,j,k,y_,iBlock) - B_GDB(i,j,k,y_,iBlock)) &
                  - SpeedOfLight_D(y_)*& 
                  (B_GDB(i,j+1,k,x_,iBlock) - B_GDB(i,j,k,x_,iBlock))
          end do
       end do
    end do
  end subroutine get_b_from_a_local
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
  !============================
 
  subroutine add_density(NodeIn_D,HighFFIn_ID,iBlock)
    !Adds an input to the density, from a given particle
    integer,dimension(nDim),intent(in)::NodeIn_D
    real,dimension(1+lOrderFF,nDim),intent(in)::HighFFIn_ID
    integer,intent(in)::iBlock
    integer::i,j,k,i1,j1,k1
    real::FFProduct
    !------------------
    k1=0
    do k=NodeIn_D(3)+iDownFF,NodeIn_D(3)+iUpFF
       k1=k1+1
       j1=0
       do j=NodeIn_D(2)+iDownFF,NodeIn_D(2)+iUpFF
          j1=j1+1
          i1=0
          FFProduct=HighFFIn_ID(k1,z_)*HighFFIn_ID(j1,y_)
          do i=NodeIn_D(1)+iDownFF,NodeIn_D(1)+iUpFF
             i1=i1+1
             Rho_GB(i,j,k,iBlock)=Rho_GB(i,j,k,iBlock) + HighFFIn_ID(i1,x_)*&
                  FFProduct   !Add the product of formfactors
          end do
       end do
    end do
  end subroutine add_density
  !=========================
  subroutine add_DensityVelocity(V_D,NodeIn_D,HighFFIn_ID,iBlock)
    !Adds an input to the density and velocity, from a given particle
    !using the same form factors
    real, dimension(nDim),intent(in) :: V_D
    integer,dimension(nDim),intent(in)::NodeIn_D
    real,dimension(1+lOrderFF,nDim),intent(in)::HighFFIn_ID
    integer,intent(in):: iBlock
    integer:: i,j,k,i1,j1,k1
    real:: FFProduct
    !----------------------                                     
    k1=0
    do k=NodeIn_D(3)+iDownFF,NodeIn_D(3)+iUpFF
       k1=k1+1
       j1=0
       do j=NodeIn_D(2)+iDownFF,NodeIn_D(2)+iUpFF
          j1=j1+1
          i1=0
          FFProduct=HighFFIn_ID(k1,z_)*HighFFIn_ID(j1,y_)
          do i=NodeIn_D(1)+iDownFF,NodeIn_D(1)+iUpFF
             i1=i1+1
             !Add the product of formfactors
             Rho_GB(i,j,k,iBlock)=Rho_GB(i,j,k,iBlock) + HighFFIn_ID(i1,x_)*&
                  FFProduct
             V_GDB(:,i,j,k,iBlock)=V_GDB(:,i,j,k,iBlock) + HighFFIn_ID(i1,x_)*&
                  FFProduct*V_D  
          end do
       end do
    end do

  end subroutine add_DensityVelocity

  !=========================
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
    real,parameter::sqrt13=cOne/1.732050808

    real,dimension(1:lOrderFF+1,nDim)::DeltaFPlus,DeltaFMinus
    real,dimension(1:lOrderFF  ,nDim)::FI
    integer::i, j, k, iDim, n_D(3), n1_D(3)
    !-------------------
    IsExtended = any(Node_D/=NodeNew_D)



    if(IsExtended)then
       call get_current_extended

       iD_D = min(NodeNew_D - Node_D,0)
       iU_D = max(NodeNew_D - Node_D,0)

       Counter_GDB(&
            Node_D(1)+iDownFF+iD_D(1):Node_D(1)+iUpFF+iU_D(1),&
            Node_D(2)+iDownFF+iD_D(2):Node_D(2)+iUpFF+iU_D(2),&
            Node_D(3)+iDownFF+iD_D(3):Node_D(3)+iUpFF+iU_D(3),&
            1:3,iBlock) = &
            Counter_GDB( Node_D(1)+iDownFF+iD_D(1):Node_D(1)+iUpFF+iU_D(1),&
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

       !To optimize algebra, multiply FI by -q*4*\pi*dt/4
       !and divide \DeltaFMinus by sqrt(3):
       do iDim=x_,z_
          FI(:,iDim)=(-QPerVDx_D(iDim)*cPi)*FI(:,iDim)
       end do

       DeltaFMinus=sqrt13*DeltaFMinus

       n_D = Node_D + iDownFF; n1_d = n_D -1


       !add current density:
       do k=1,lOrderFF+1
          do j=1,lOrderFF+1
             !\
             !  x_ currents
             !/
             Counter_GDB(n_D(1):lOrderFF+n1_D(1),&
                  j+n1_D(2),k+n1_D(3),x_,iBlock) = &
                  Counter_GDB(n_D(1):lOrderFF+n1_D(1),&
                  j+n1_D(2),k+n1_D(3),x_,iBlock) +&
                  FI(:,x_)*(&
                  DeltaFPlus (j,y_)*DeltaFPlus (k,z_)+&
                  DeltaFMinus(j,y_)*DeltaFMinus(k,z_))
          end do
          do j=1,lOrderFF
             !\
             !  y_ currents
             !/
             Counter_GDB(n_D(1):lOrderFF+n_D(1),&
                  j+n1_D(2),k+n1_D(3),y_,iBlock) = &
                  Counter_GDB(n_D(1):lOrderFF+n_D(1),&
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
             Counter_GDB(n_D(1):lOrderFF+n_D(1),&
                  j+n1_D(2),k+n1_D(3),z_,iBlock) = &
                  Counter_GDB(n_D(1):lOrderFF+n_D(1),&
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
         FIExt(:,iDim) = (-QPerVDx_D(iDim)*cPi)*FIExt(:,iDim)
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
    !===================
    subroutine write_current(iBlock)
      use PIC_ModMain, ONLY: Dt, SpeedOfLight_D
      integer,intent(in)::iBlock
      integer::i,j,k,iDim
      !-------------------
      write(*,*)'cPi*QPerVDx_D, Dt=',cPi*QPerVDx_D, Dt
      write(*,*)'W_D=',W_D
      write(*,*)'Node_D=',Node_D,' FormfactorOld ='
      do iDim = 1,3
         write(*,*)HighFF_ID(:,iDim)
      end do
      write(*,*)'NodeNew_D=',NodeNew_D,' FormfactorNew ='
      do iDim = 1,3
         write(*,*)HighFFNew_ID(:,iDim)
      end do
      do  iDim = 1,3 
         write(*,*)'iDim=',iDim,' W_D(iDim)*Dt='!, W_D(iDim)*Dt
         if(IsExtended)then
            do k=Node_D(3)+iDownFF-1,Node_D(3)+iUpFF+1
               do j= Node_D(2)+iDownFF-1,Node_D(2)+iUpFF+1
                  write(*,*)'Current_G(',Node_D(1)+iDownFF-1,&
                       ':',Node_D(1)+iUpFF+1,',',j,',',k,')= ',&
                       Counter_GDB(Node_D(1)+iDownFF-1:Node_D(1)+iUpFF+1,&
                       j,k,iDim,iBlock)
               end do
            end do
            write(*,*)'Sum of currents=',sum(Counter_GDB(&
                 Node_D(1)+iDownFF-1:Node_D(1)+iUpFF+1,&
                 Node_D(2)+iDownFF-1:Node_D(2)+iUpFF+1,&
                 Node_D(3)+iDownFF-1:Node_D(3)+iUpFF+1,&
                 iDim,iBlock))

         else
            do k=Node_D(3)+iDownFF,Node_D(3)+iUpFF
               do j= Node_D(2)+iDownFF,Node_D(2)+iUpFF
                  write(*,*)'Current_G(',Node_D(1)+iDownFF,&
                       ':',Node_D(1)+iUpFF,',',j,',',k,')= ',&
                       Counter_GDB(Node_D(1)+iDownFF:Node_D(1)+iUpFF,&
                       j,k,iDim,iBlock)
               end do
            end do
            write(*,*)'Sum of currents=',sum(Counter_GDB(&
                 Node_D(1)+iDownFF  :Node_D(1)+iUpFF  ,&
                 Node_D(2)+iDownFF  :Node_D(2)+iUpFF  ,&
                 Node_D(3)+iDownFF  :Node_D(3)+iUpFF  ,&
                 iDim,iBlock))
         end if
         write(*,*)'Compare with 4*\pi*q*u*dt/Volume=',&
              4.0*cPi*QPerVDx_D(iDim)*W_D(iDim)*SpeedOfLight_D(iDim)
      end do
    end subroutine write_current
    !===========================
  end subroutine add_current
  !=========================
  real function min_val_rho()
    min_val_rho = &
         minval(Rho_GB(1:nX,1:nY,1:nZ,:))
  end function min_val_rho
  !=======================
  !=========================
  real function max_val_rho()
    max_val_rho = &
         maxval(Rho_GB(1:nX,1:nY,1:nZ,:))
  end function max_val_rho
  !=======================
  real function rho_avr()
    use PIC_ModSize, ONLY: nCell_D
    rho_avr = &
         sum(Rho_GB(1:nX,1:nY,1:nZ,:))/product(nCell_D)
  end function rho_avr
  !=======================
end module PIC_ModParticleInField
!=================================================================!
