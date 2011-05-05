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
    real,dimension(1:3,1:lOrderFF+1,1:lOrderFF+1),&
           intent(out)::BOut_DG
    integer::i,j,i1,j1
    !Applied only if UseVectorPotential
    !The vector potential is in Magnetic_DG
    !Calculates the magnetic field.
    !
    !The resulting magnetic field is saved to BOut_DG
    !Only those values of the field are calculated which are needed
    !to interpolated the magnetic field acting on the particle
    !around Node_D
    j1=0
    do j=Node_D(2)+iDownFF,Node_D(2)+iUpFF
       j1=j1+1
       i1=0
       do i=Node_D(1)+iDownFF,Node_D(1)+iUpFF
          i1=i1+1
          BOut_DG(x_,i1,j1)=&
               SpeedOfLight_D(y_)*&
                 (Magnetic_DG(z_,i,  j+1)-Magnetic_DG(z_,i,j))
             
          BOut_DG(y_,i1,j1)=-&
               SpeedOfLight_D(x_)*&
                 (Magnetic_DG(z_,i+1,j  )-Magnetic_DG(z_,i,j))
             
          BOut_DG(z_,i1,j1)=&
               SpeedOfLight_D(x_)*&
                 (Magnetic_DG(y_,i+1,j  )-Magnetic_DG(y_,i,j))-&
               SpeedOfLight_D(y_)*& 
                 (Magnetic_DG(x_,i,  j+1)-Magnetic_DG(x_,i,j))
       end do
    end do
  end subroutine get_b_from_a_local
  !-------------------------------------------------------------------------!
  function e_interpolated_d()
    !Sets the pointer to the fragment of E_DG array to interpolate E
    real,dimension(1:3)::e_interpolated_d
    call interpolate_e( E_DG(1:3,Node_D(1)+iDownFF:Node_D(1)+iUpFF,&
                         Node_D(2)+iDownFF:Node_D(2)+iUpFF))
  contains
    subroutine interpolate_e(EIn_DG)
      !Calculate the interpolated value of the electric field
      !Optimized loops are used, the ranges of indexes being
      !declared as parameters
      use ModNumConst
      real,dimension(1:3,1:lOrderFF+1,1:lOrderFF+1),&
           intent(in)::EIn_DG
      integer::i,j
      e_interpolated_d=cZero
      do j=1,lOrderFF+1
         e_interpolated_d(x_)= e_interpolated_d(x_)+ &
               sum(EIn_DG(x_,1:lOrderFF,j)*LowFF_ID(:,x_))&
                 *HighFF_ID(j,y_)
      end do
      do j=1,lOrderFF
         e_interpolated_d(y_)= e_interpolated_d(y_)+ &
               sum(EIn_DG(y_,1:lOrderFF+1,j)*HighFF_ID(:,x_))&
                 *LowFF_ID(j,y_)
      end do
      do j=1,lOrderFF+1
         e_interpolated_d(z_)= e_interpolated_d(z_)+ &
               sum(EIn_DG(z_,1:lOrderFF+1,j)*HighFF_ID(:,x_))&
                 *HighFF_ID(j,y_)
      end do
    end subroutine interpolate_e
  end function e_interpolated_d
!-------------------------------------------------------------------------!
  function b_interpolated_d()
    !Sets the pointer to the fragment of the magnetic field array to 
    !interpolate B
    real,dimension(1:3)::b_interpolated_d
    real,dimension(1:3,1:lOrderFF+1,1:lOrderFF+1)::BOut_DG
    if(UseVectorPotential)then
       call get_b_from_a(BOut_DG)
       call interpolate_b(BOut_DG)
    else   
       call interpolate_b(Magnetic_DG(1:3,&
                         Node_D(1)+iDownFF:Node_D(1)+iUpFF,&
                         Node_D(2)+iDownFF:Node_D(2)+iUpFF))
    end if
  contains
    subroutine interpolate_b(BIn_DG)
      !calculate the interpolated value of the magnetic field
      !Optimized loops are used, the ranges of indexes being
      !declared as parameters
      real,dimension(1:3,1:lOrderFF+1,1:lOrderFF+1),&
           intent(in)::BIn_DG
      integer::i,j
      b_interpolated_d=cZero
      do j=1,lOrderFF
         b_interpolated_d(x_)= b_interpolated_d(x_)+ &
               sum(BIn_DG(x_,1:lOrderFF+1,j)*HighFF_ID(:,x_))&
                 *LowFF_ID(j,y_)
      end do
      do j=1,lOrderFF+1
         b_interpolated_d(y_)= b_interpolated_d(y_)+ &
              sum(BIn_DG(y_,1:lOrderFF,j)*LowFF_ID(:,x_))&
                 *HighFF_ID(j,y_)
      end do
      do j=1,lOrderFF
         b_interpolated_d(z_)= b_interpolated_d(z_)+ &
              sum(BIn_DG(z_,1:lOrderFF,j)*LowFF_ID(:,x_))&
                 *LowFF_ID(j,y_)
      end do
    end subroutine interpolate_b
  end function b_interpolated_d
!-------------------------------------------------------------------------!
  subroutine add_density(NodeIn_D,HighFFIn_ID)
    !Adds an input to the density, from a given particle
    integer,dimension(nDim),intent(in)::NodeIn_D
    real,dimension(1+lOrderFF,nDim)::HighFFIn_ID
    integer::i,j,i1,j1
    real::FFProduct
    j1=0
    do j=NodeIn_D(2)+iDownFF,NodeIn_D(2)+iUpFF
       j1=j1+1
       i1=0
       FFProduct=HighFFIn_ID(j1,x_)
       do i=NodeIn_D(1)+iDownFF,NodeIn_D(1)+iUpFF
          i1=i1+1
          rho_G(i,j)=rho_G(i,j)+HighFFIn_ID(i1,x_)*&
                  FFProduct   !Add the product of formfactors
       end do
    end do
  end subroutine add_density
!-----------------------------------------------------------------!
  subroutine add_current(QPerVDx_D,W_D) 
    real,intent(in)::QPerVDx_D(nDim),W_D(x_:z_)
    optional::W_D
    if(any(Node_D/=NodeNew_D))then
       call add_current_extended(Counter_DG(x_:z_,&
                         Node_D(1)+iDownFF-1:Node_D(1)+iUpFF+1,&
                         Node_D(2)+iDownFF-1:Node_D(2)+iUpFF+1))
    else
       call add_current_simple(Counter_DG(x_:z_,&
                         Node_D(1)+iDownFF  :Node_D(1)+iUpFF  ,&
                         Node_D(2)+iDownFF  :Node_D(2)+iUpFF  ))
    end if

  contains
    subroutine add_current_simple(CurrentFragment_DG)
      real,dimension(1:3,1:lOrderFF+1,1:lOrderFF+1),&
           intent(inout)::CurrentFragment_DG
      real,dimension(1:lOrderFF+1,nDim)::DeltaFPlus,DeltaFMinus
      real,dimension(1:lOrderFF,nDim)::FI
      real,parameter::sqrt13=cOne/1.732050808
      integer::i,j,iDim
      DeltaFPlus=HighFFNew_ID+HighFF_ID
      DeltaFMinus=HighFFNew_ID-HighFF_ID
      FI(1,:)=DeltaFMinus(1,:)
      !Use recursive formula for FI
      do i=2,lOrderFF
         FI(i,:)=FI(i-1,:)+DeltaFMinus(i,:)
      end do
      !To optimize algebra, multiply FI by -q/V*dx/4
      !and divide \DeltaFMinus by sqrt(3):
      do iDim=x_,z_
         FI(:,iDim)=(-QPerVDx_D(iDim)*cQuarter)*FI(:,iDim)
      end do
      DeltaFMinus=sqrt13*DeltaFMinus
      !add current density:
      do j=1,lOrderFF+1
         CurrentFragment_DG(x_,1:lOrderFF,j)=&
              CurrentFragment_DG(x_,1:lOrderFF,j)+&
              FI(:,x_)*(&
              DeltaFPlus (j,y_))
      end do
      do j=1,lOrderFF
         CurrentFragment_DG(y_,1:lOrderFF+1,j)=&
              CurrentFragment_DG(x_,1:lOrderFF+1,j)+&
              FI(j,y_)*(&
              DeltaFPlus (:,x_))
      end do
      do j=1,lOrderFF+1
         CurrentFragment_DG(z_,1:lOrderFF+1,j)=&
              CurrentFragment_DG(z_,1:lOrderFF+1,j)+&
              W_D(3)*(&
              DeltaFPlus (:,x_)*DeltaFPlus (j,y_)+&
              DeltaFMinus(:,x_)*DeltaFMinus(j,y_))
      end do
    end subroutine add_current_simple
    !-------------------------------------------------------------!
    subroutine add_current_extended(CurrentFragment_DG)
      !The same, but with the extended stencil for the current
      real,dimension(1:3,0:lOrderFF+2,0:lOrderFF+2),&
           intent(inout)::CurrentFragment_DG
      real,dimension(0:lOrderFF+2,nDim)::DeltaFPlus,DeltaFMinus
      real,dimension(0:lOrderFF+1,nDim)::FI
      real,parameter::sqrt13=cOne/1.732050808
      integer::i,j
      DeltaFPlus=cZero
      !First, add the properly shifted new formfactor
      DeltaFPlus(1+NodeNew_D(1)-Node_D(1):&
           1+NodeNew_D(1)-Node_D(1)+lOrderFF,x_)=&
           HighFFNew_ID(:,x_)
      DeltaFPlus(1+NodeNew_D(2)-Node_D(2):&
           1+NodeNew_D(2)-Node_D(2)+lOrderFF,y_)=&
           HighFFNew_ID(:,y_)
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
      !To optimize algebra, multiply FI by -q/V*dx/4
      !and divide \DeltaFMinus by sqrt(3):
      do iDim=x_,z_
         FI(:,iDim)=(-QPerVDx_D(iDim)*cQuarter)*FI(:,iDim)
      end do
      DeltaFMinus=sqrt13*DeltaFMinus
      !add current density:
      do j=0,lOrderFF+2
         CurrentFragment_DG(x_,0:lOrderFF+1,j)=&
              CurrentFragment_DG(x_,0:lOrderFF+1,j)+&
              FI(:,x_)*(&
              DeltaFPlus (j,y_))
      end do
      do j=0,lOrderFF+1
         CurrentFragment_DG(y_,0:lOrderFF+2,j)=&
              CurrentFragment_DG(x_,0:lOrderFF+2,j)+&
              FI(j,y_)*(&
              DeltaFPlus (:,x_))
      end do
      do j=0,lOrderFF+2
         CurrentFragment_DG(z_,0:lOrderFF+2,j)=&
              CurrentFragment_DG(z_,0:lOrderFF+2,j)+&
              W_D(3)*(&
              DeltaFPlus (:,x_)*DeltaFPlus (j,y_)+&
              DeltaFMinus(:,x_)*DeltaFMinus(j,y_))
      end do
    end subroutine add_current_extended
  end subroutine add_current
end module PIC_ModParticleInField
!=================================================================!
