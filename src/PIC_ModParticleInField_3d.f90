module PIC_ModParticleInField
  use PIC_ModField,get_b_from_a_global=>get_b_from_a
!-----------------------------------------------------------------------!
  !Structures
  integer,dimension(nDim)::Node_D,NodeNew_D
               !Discsrete particle coordinates

  real,dimension(1:lOrderFF+1,nDim)::HighFF_ID,HighFFNew_ID
               !Particle form-factor

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
  public::add_current !Adds an input to an electric current
contains
!----------------------------------------------------------------------!
  subroutine get_b_from_a_local(BOut_DG)    
    real,dimension(1:3,1:lOrderFF+1,1:lOrderFF+1,1:lOrderFF+1),&
           intent(out)::BOut_DG
    integer::i,j,k,i1,j1,k1
    !Applied only if UseVectorPotential
    !The vector potential is in Magnetic_DG
    !Calculates the magnetic field.
    !
    !The resulting magnetic field is saved to BOut_DG
    !Only those values of the field are calculated which are needed
    !to interpolated the magnetic field acting on the particle
    !around Node_D
    k1=0
    do k=Node_D(3)+iDownFF,Node_D(3)+iUpFF
       k1=k1+1
       j1=0
       do j=Node_D(2)+iDownFF,Node_D(2)+iUpFF
          j1=j1+1
          i1=0
          do i=Node_D(1)+iDownFF,Node_D(1)+iUpFF
             i1=i1+1
             BOut_DG(x_,i1,j1,k1)=&
                  SpeedOfLight_D(y_)*&
                  (Magnetic_DG(z_,i,  j+1,k  )-Magnetic_DG(z_,i,j,k))-&
                  SpeedOfLight_D(z_)*&
                  (Magnetic_DG(y_,i,  j,  k+1)-Magnetic_DG(y_,i,j,k))
             
             BOut_DG(y_,i1,j1,k1)=&
                  SpeedOfLight_D(z_)*&
                  (Magnetic_DG(x_,i,  j,  k+1)-Magnetic_DG(x_,i,j,k))-&
                  SpeedOfLight_D(x_)*&
                  (Magnetic_DG(z_,i+1,j,  k  )-Magnetic_DG(z_,i,j,k))
             
             BOut_DG(z_,i1,j1,k1)=&
                  SpeedOfLight_D(x_)*&
                  (Magnetic_DG(y_,i+1,j,  k  )-Magnetic_DG(y_,i,j,k))-&
                  SpeedOfLight_D(y_)*& 
                  (Magnetic_DG(x_,i,  j+1,k  )-Magnetic_DG(x_,i,j,k))
          end do
       end do
    end do
  end subroutine get_b_from_a_local
  !-------------------------------------------------------------------------!
  function e_interpolated_d()
    !Sets the pointer to the fragment of E_DG array to interpolate E
    real,dimension(1:3)::e_interpolated_d
    call interpolate_e( E_DG(1:3,Node_D(1)+iDownFF:Node_D(1)+iUpFF,&
                         Node_D(2)+iDownFF:Node_D(2)+iUpFF,&
                         Node_D(3)+iDownFF:Node_D(3)+iUpFF))
  contains
    subroutine interpolate_e(EIn_DG)
      !Calculate the interpolated value of the electric field
      !Optimized loops are used, the ranges of indexes being
      !declared as parameters
      use ModNumConst
      real,dimension(1:3,1:lOrderFF+1,1:lOrderFF+1,1:lOrderFF+1),&
           intent(in)::EIn_DG
      integer::i,j,k
      e_interpolated_d=cZero
      do k=1,lOrderFF+1
         do j=1,lOrderFF+1
            e_interpolated_d(x_)= e_interpolated_d(x_)+ &
            sum(EIn_DG(x_,1:lOrderFF,j,k)*LowFF_ID(:,x_))&
                 *HighFF_ID(j,y_)*HighFF_ID(k,z_)
         end do
         do j=1,lOrderFF
            e_interpolated_d(y_)= e_interpolated_d(y_)+ &
            sum(EIn_DG(y_,1:lOrderFF+1,j,k)*HighFF_ID(:,x_))&
                 *LowFF_ID(j,y_)*HighFF_ID(k,z_)
         end do
      end do
      do k=1,lOrderFF
         do j=1,lOrderFF+1
            e_interpolated_d(z_)= e_interpolated_d(z_)+ &
            sum(EIn_DG(z_,1:lOrderFF+1,j,k)*HighFF_ID(:,x_))&
                 *HighFF_ID(j,y_)*LowFF_ID(k,z_)
         end do
      end do
    end subroutine interpolate_e
  end function e_interpolated_d
!-------------------------------------------------------------------------!
  function b_interpolated_d()
    !Sets the pointer to the fragment of the magnetic field array to 
    !interpolate B
    real,dimension(1:3)::b_interpolated_d
    real,dimension(1:3,1:lOrderFF+1,1:lOrderFF+1,1:lOrderFF+1)&
         ::BOut_DG
    if(UseVectorPotential)then
       call get_b_from_a(BOut_DG)
       call interpolate_b(BOut_DG)
    else   
       call interpolate_b(Magnetic_DG(1:3,&
                         Node_D(1)+iDownFF:Node_D(1)+iUpFF,&
                         Node_D(2)+iDownFF:Node_D(2)+iUpFF,&
                         Node_D(3)+iDownFF:Node_D(3)+iUpFF))
    end if
  contains
    subroutine interpolate_b(BIn_DG)
      !calculate the interpolated value of the magnetic field
      !Optimized loops are used, the ranges of indexes being
      !declared as parameters
      real,dimension(1:3,1:lOrderFF+1,1:lOrderFF+1,1:lOrderFF+1),&
           intent(in)::BIn_DG
      integer::i,j,k
      b_interpolated_d=cZero
      do k=1,lOrderFF
         do j=1,lOrderFF
            b_interpolated_d(x_)= b_interpolated_d(x_)+ &
            sum(BIn_DG(x_,1:lOrderFF+1,j,k)*HighFF_ID(:,x_))&
                 *LowFF_ID(j,y_)*LowFF_ID(k,z_)
         end do
         do j=1,lOrderFF+1
            b_interpolated_d(y_)= b_interpolated_d(y_)+ &
            sum(BIn_DG(y_,1:lOrderFF,j,k)*LowFF_ID(:,x_))&
                 *HighFF_ID(j,y_)*LowFF_ID(k,z_)
         end do
      end do
      do k=1,lOrderFF+1
         do j=1,lOrderFF
            b_interpolated_d(z_)= b_interpolated_d(z_)+ &
            sum(BIn_DG(z_,1:lOrderFF,j,k)*LowFF_ID(:,x_))&
                 *LowFF_ID(j,y_)*HighFF_ID(k,z_)
         end do
      end do
    end subroutine interpolate_b
  end function b_interpolated_d
!-------------------------------------------------------------------------!
  subroutine add_density(NodeIn_D,HighFFIn_ID)
    !Adds an input to the density, from a given particle
    integer,dimension(nDim),intent(in)::NodeIn_D
    real,dimension(1+lOrderFF,nDim)::HighFFIn_ID
    integer::i,j,k,i1,j1,k1
    real::FFProduct
    k1=0
    do k=NodeIn_D(3)+iDownFF,NodeIn_D(3)+iUpFF
       k1=k1+1
       j1=0
       do j=NodeIn_D(2)+iDownFF,NodeIn_D(2)+iUpFF
          j1=j1+1
          i1=0
          FFProduct=HighFFIn_ID(k1,x_)*HighFFIn_ID(j1,x_)
          do i=NodeIn_D(1)+iDownFF,NodeIn_D(1)+iUpFF
             i1=i1+1
             rho_G(i,j,k)=rho_G(i,j,k)+HighFFIn_ID(i1,x_)*&
                  FFProduct   !Add the product of formfactors
          end do
       end do
    end do
  end subroutine add_density
!-----------------------------------------------------------------!
  subroutine add_current(QPerVDx_D,W_D) 
    real,intent(in)::QPerVDx_D(x_:z_),W_D(x_:z_)
    optional::W_D
    if(any(Node_D/=NodeNew_D))then
       call add_current_extended(Counter_DG(1:3,&
                         Node_D(1)+iDownFF-1:Node_D(1)+iUpFF+1,&
                         Node_D(2)+iDownFF-1:Node_D(2)+iUpFF+1,&
                         Node_D(3)+iDownFF-1:Node_D(3)+iUpFF+1))
    else
       call add_current_simple(Counter_DG(1:3,&
                         Node_D(1)+iDownFF  :Node_D(1)+iUpFF  ,&
                         Node_D(2)+iDownFF  :Node_D(2)+iUpFF  ,&
                         Node_D(3)+iDownFF  :Node_D(3)+iUpFF  ))
    end if

  contains
    subroutine add_current_simple(CurrentFragment_DG)
      real,dimension(1:3,1:lOrderFF+1,1:lOrderFF+1,1:lOrderFF+1),&
           intent(inout)::CurrentFragment_DG
      real,dimension(1:lOrderFF+1,nDim)::DeltaFPlus,DeltaFMinus
      real,dimension(1:lOrderFF,nDim)::FI
      real,parameter::sqrt13=cOne/1.732050808
      integer::i,j,k,iDim
      DeltaFPlus=HighFFNew_ID+HighFF_ID
      DeltaFMinus=HighFFNew_ID-HighFF_ID
      FI(1,:)=DeltaFMinus(1,:)
      !Use recursive formula for FI
      do i=2,lOrderFF
         FI(i,:)=FI(i-1,:)+DeltaFMinus(i,:)
      end do
      !To optimize the algebra, multiply FI by -q*dt/4
      !and divide \DeltaFMinus by sqrt(3):
      do iDim=x_,z_
         FI(:,iDim)=(-QPerVDx_D(iDim)*cQuarter)*FI(:,iDim)
      end do
      DeltaFMinus=sqrt13*DeltaFMinus
      !add current density:
      do k=1,lOrderFF+1
         do j=1,lOrderFF+1
            CurrentFragment_DG(x_,1:lOrderFF,j,k)=&
                 CurrentFragment_DG(x_,1:lOrderFF,j,k)+&
                 FI(:,x_)*(&
                 DeltaFPlus (j,y_)*DeltaFPlus (k,z_)+&
                 DeltaFMinus(j,y_)*DeltaFMinus(k,z_))
         end do
         do j=1,lOrderFF
            CurrentFragment_DG(y_,1:lOrderFF+1,j,k)=&
                 CurrentFragment_DG(x_,1:lOrderFF+1,j,k)+&
                 FI(j,y_)*(&
                 DeltaFPlus (:,x_)*DeltaFPlus (k,z_)+&
                 DeltaFMinus(:,x_)*DeltaFMinus(k,z_))
         end do
      end do
      do k=1,lOrderFF
         do j=1,lOrderFF+1
            CurrentFragment_DG(z_,1:lOrderFF+1,j,k)=&
                 CurrentFragment_DG(z_,1:lOrderFF+1,j,k)+&
                 FI(k,z_)*(&
                 DeltaFPlus (:,x_)*DeltaFPlus (j,y_)+&
                 DeltaFMinus(:,x_)*DeltaFMinus(j,y_))
         end do
      end do
    end subroutine add_current_simple
    !-------------------------------------------------------------!
    subroutine add_current_extended(CurrentFragment_DG)
      !The same, but with the extended stencil for the current
      real,dimension(1:3,0:lOrderFF+2,0:lOrderFF+2,0:lOrderFF+2),&
           intent(inout)::CurrentFragment_DG
      real,dimension(0:lOrderFF+2,nDim)::DeltaFPlus,DeltaFMinus
      real,dimension(0:lOrderFF+1,nDim)::FI
      real,parameter::sqrt13=cOne/1.732050808
      integer::i,j,k
      DeltaFPlus=cZero
      !First, add the properly shifted new formfactor
      DeltaFPlus(1+NodeNew_D(1)-Node_D(1):&
           1+NodeNew_D(1)-Node_D(1)+lOrderFF,x_)=&
           HighFFNew_ID(:,x_)
      DeltaFPlus(1+NodeNew_D(2)-Node_D(2):&
           1+NodeNew_D(2)-Node_D(2)+lOrderFF,y_)=&
           HighFFNew_ID(:,y_)
      DeltaFPlus(1+NodeNew_D(3)-Node_D(3):&
           1+NodeNew_D(3)-Node_D(3)+lOrderFF,z_)=&
           HighFFNew_ID(:,z_)
      DeltaFMinus=DeltaFPlus
      !Then add and subtruct the unshifted old formfactor           
      DeltaFPlus(1:lOrderFF+1,:) =&
           DeltaFPlus (1:lOrderFF+1,:)+HighFF_ID
      DeltaFMinus(1:lOrderFF+1,:)=&
           DeltaFMinus(1:lOrderFF+1,:)-HighFF_ID
      FI(0,:)=DeltaFMinus(0,:)
      !Use recursive formula for FI
      do i=1,lOrderFF+1
         FI(i,:)=FI(i-1,:)+DeltaFMinus(i,:)
      end do
      !To optimize the algebra, multiply FI by -q*dt/4
      !and divide \DeltaFMinus by sqrt(3):
      do iDim=x_,z_
         FI(:,iDim)=(-QPerVDx_D(iDim)*cQuarter)*FI(:,iDim)
      end do
      DeltaFMinus=sqrt13*DeltaFMinus
      !add current density:
      do k=0,lOrderFF+2
         do j=0,lOrderFF+2
            CurrentFragment_DG(x_,0:lOrderFF+1,j,k)=&
                 CurrentFragment_DG(x_,0:lOrderFF+1,j,k)+&
                 FI(:,x_)*(&
                 DeltaFPlus (j,y_)*DeltaFPlus (k,z_)+&
                 DeltaFMinus(j,y_)*DeltaFMinus(k,z_))
         end do
         do j=0,lOrderFF+1
            CurrentFragment_DG(y_,0:lOrderFF+2,j,k)=&
                 CurrentFragment_DG(x_,0:lOrderFF+2,j,k)+&
                 FI(j,y_)*(&
                 DeltaFPlus (:,x_)*DeltaFPlus (k,z_)+&
                 DeltaFMinus(:,x_)*DeltaFMinus(k,z_))
         end do
      end do
      do k=0,lOrderFF+1
         do j=0,lOrderFF+2
            CurrentFragment_DG(z_,0:lOrderFF+2,j,k)=&
                 CurrentFragment_DG(z_,0:lOrderFF+2,j,k)+&
                 FI(k,z_)*(&
                 DeltaFPlus (:,x_)*DeltaFPlus (j,y_)+&
                 DeltaFMinus(:,x_)*DeltaFMinus(j,y_))
         end do
      end do
    end subroutine add_current_extended
  end subroutine add_current
end module PIC_ModParticleInField
!=================================================================!
