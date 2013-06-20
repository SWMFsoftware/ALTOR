!^CFG COPYRIGHT UofM
!--------------------------------------------------------------!

module PIC_ModField
  use PIC_ModSize, ONLY: nX, nY, nZ, x_, y_, z_, nDim, nCell_D
  use ModNumConst
  use PIC_ModMain,ONLY:&
       SpeedOfLight_D,&     !This is c\Delta t/\Delta x,...
       UseVectorPotential   !If we use vector-potential, or not
  use PIC_ModMain, ONLY: Dx_D, IsPeriodicField_D, TypeFieldBc_S
  use PIC_ModFormFactor, ONLY: lOrderFF
  implicit none

  SAVE


  !Introduce the number of ghostcells (iGCN)
  integer,parameter:: iGCN=&
       1 &    !because even with the first order formfactor the 
                                !particle can escape the domain within a time step
       +lOrderFF/2 
  !because for lOrderFF>1 the particle form-factor
  !is wider than a mesh size

  !Structures
  !Array index is the coordinate of the gridpoint with +1/2 being
  !approximated as 1
  !                                                Ez (:,:,nZ+2)    
  !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:)
  !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:)
  !       _!_!_!                x!x!_!_  Bz(:,nY+2,:),Bz(nX+2,:,:)       
  !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:),Bx(:,:,nZ+2)
  !       -2                             By(nX+2,:,:),By(:,:,nZ+2)
  real,dimension(&
       1-iGCN:nX+iGCN,&
       1-iGCN:nY+iGCN,&
       1-iGCN:nZ+iGCN,3)::&
       E_GD        = 0.0,& !This is the electric field
       Magnetic_GD = 0.0,& !vector-potential if used, magnetic field otherwise
       Counter_GD  = 0.0   !Counter for electric current

  real,dimension(&
       1-iGCN:nX+iGCN,&
       1-iGCN:nY+iGCN,&
       1-iGCN:nZ+iGCN)::rho_G
  real :: B0_D(3) = 0.0

  !Methods

  public::get_b_from_a   !Transforms the vector potential to magn. field
  public::update_magnetic!Updates the magnetic field, or vector potential
  public::update_e       !Updates the electric field
  public::get_rho_max     
contains
  !=======================
  subroutine add_e
    use ModReadParam, ONLY: read_var
    real:: E_D(3)
    integer:: i,j,k,iDim
    !------------
    call read_var('Ex',E_D(1))
    call read_var('Ey',E_D(2))
    call read_var('Ez',E_D(3))
    do iDim = 1,3
       do k=1-iGCN,nZ+iGCN
          do j=1-iGCN,nY+iGCN
             do i=1-iGCN,nX+iGCN
                E_GD(i,j,k,iDim) = &
                     E_GD(i,j,k,iDim) + E_D(iDim)
             end do
          end do
       end do
    end do
  end subroutine add_e
  !===================
  subroutine add_b
    use ModReadParam, ONLY: read_var
    integer:: i,j,k,iDim
    !------------
    call read_var('Bx',B0_D(1))
    call read_var('By',B0_D(2))
    call read_var('Bz',B0_D(3))
    do iDim = 1,3
       do k=1-iGCN,nZ+iGCN
          do j=1-iGCN,nY+iGCN
             do i=1-iGCN,nX+iGCN
                Magnetic_GD(i,j,k,iDim) = &
                     Magnetic_GD(i,j,k,iDim) + B0_D(iDim)
             end do
          end do
       end do
    end do
  end subroutine add_b
  !===================
  subroutine get_b_from_a
    integer::i,j,k
    !Applied only if UseVectorPotential
    !The vector potential is in Magnetic_GD
    !Calculates the magnetic field 

    !The whole magnetic field is calculated and 
    !the result is saved to Counter_GD
    !---------------------------------
    Counter_GD = 0.0

    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1
    !                                                Ez (:,:,nZ+2)    
    !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:)
    !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:)
    !       _!_!_!                x!x!_!_  Bz(:,nY+2,:),Bz(nX+2,:,:)       
    !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:),Bx(:,:,nZ+2)
    !       -2                             By(nX+2,:,:),By(:,:,nZ+2)
    do k=1-iGCN,nZ+iGCN-1; do j=1-iGCN,nY+iGCN-1; do i=1-iGCN,nX+iGCN
       Counter_GD(i,j,k,x_) = &
            SpeedOfLight_D(y_)*&
            (Magnetic_GD(i,  j+1,k  ,z_) - Magnetic_GD(i,j,k,z_))-&
            SpeedOfLight_D(z_)*&
            (Magnetic_GD(i,  j,  k+1,y_) - Magnetic_GD(i,j,k,y_))
    end do; end do; end do

    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1
    !                                                Ez (:,:,nZ+2)    
    !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:)
    !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:)
    !       _!_!_!                x!x!_!_  Bz(:,nY+2,:),Bz(nX+2,:,:)       
    !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:),Bx(:,:,nZ+2)
    !       -2                             By(nX+2,:,:),By(:,:,nZ+2)
    do k=1-iGCN,nZ+iGCN-1; do j=1-iGCN,nY+iGCN; do i=1-iGCN,nX+iGCN-1
       Counter_GD(i,j,k,y_) = &
            SpeedOfLight_D(z_)*&
            (Magnetic_GD(i,  j,  k+1,x_) - Magnetic_GD(i,j,k,x_))-&
            SpeedOfLight_D(x_)*&
            (Magnetic_GD(i+1,j,  k  ,z_) - Magnetic_GD(i,j,k,z_))
    end do; end do; end do

    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1
    !                                                Ez (:,:,nZ+2)    
    !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:)
    !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:)
    !       _!_!_!                x!x!_!_  Bz(:,nY+2,:),Bz(nX+2,:,:)       
    !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:),Bx(:,:,nZ+2)
    !       -2                             By(nX+2,:,:),By(:,:,nZ+2)
    do k=1-iGCN,nZ+iGCN; do j=1-iGCN,nY+iGCN-1; do i=1-iGCN,nX+iGCN-1
       Counter_GD(i,j,k,z_) = &
            SpeedOfLight_D(x_)*&
            (Magnetic_GD(i+1,j,  k  ,y_) - Magnetic_GD(i,j,k,y_))-&
            SpeedOfLight_D(y_)*& 
            (Magnetic_GD(i,  j+1,k  ,x_) - Magnetic_GD(i,j,k,x_))
    end do; end do; end do
  end subroutine get_b_from_a
  !==========================

  subroutine update_magnetic
    integer::i,j,k
    real,dimension(nDim):: SpeedOfLightHalf_D
    !--------------------------
    if(UseVectorPotential)then

       !Advance vector potential throught a half timestep
       !Array index is the coordinate of the gridpoint with +1/2 being
       !approximated as 1
       !                                                Ez (:,:,nZ+2)    
       !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:)
       !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:)
       !       _!_!_!                x!x!_!_  Bz(:,nY+2,:),Bz(nX+2,:,:)       
       !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:),Bx(:,:,nZ+2)
       !       -2                             By(nX+2,:,:),By(:,:,nZ+2)
       Magnetic_GD(     1-iGCN:nX+iGCN-1,:,:,x_) = &
            Magnetic_GD(1-iGCN:nX+iGCN-1,:,:,x_) - &
            cHalf*E_GD(1-iGCN:nX+iGCN-1,:,:,x_)

       Magnetic_GD(     :,1-iGCN:nY+iGCN-1,:,y_) = &
            Magnetic_GD(:,1-iGCN:nY+iGCN-1,:,y_) - &
            cHalf*E_GD(:,1-iGCN:nY+iGCN-1,:,y_)

       Magnetic_GD(     :,:,1-iGCN:nZ+iGCN-1,z_) = &
            Magnetic_GD(:,:,1-iGCN:nZ+iGCN-1,z_) - &
            cHalf*E_GD(:,:,1-iGCN:nZ+iGCN-1,z_)
       return
    end if

    SpeedOfLightHalf_D = &
         SpeedOfLight_D*cHalf

    !Advance vector potential throught a half timestep
    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1
    !                                                Ez (:,:,nZ+2)    
    !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:)
    !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:)
    !       _!_!_!                x!x!_!_  Bz(:,nY+2,:),Bz(nX+2,:,:)       
    !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:),Bx(:,:,nZ+2)
    !       -2                             By(nX+2,:,:),By(:,:,nZ+2)
    do k=1-iGCN,nZ+iGCN-1; do j=1-iGCN,nY+iGCN-1; do i=1-iGCN,nX+iGCN
       Magnetic_GD(i,j,k,x_) = Magnetic_GD(i,j,k,x_) - &
            SpeedOfLightHalf_D(y_)*&
            (E_GD(i,  j+1,k  ,z_) - E_GD(i,j,k,z_)) + &
            SpeedOfLightHalf_D(z_)*&
            (E_GD(i,  j,  k+1,y_) - E_GD(i,j,k,y_))
    end do; end do; end do

    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1
    !                                                Ez (:,:,nZ+2)    
    !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:)
    !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:)
    !       _!_!_!                x!x!_!_  Bz(:,nY+2,:),Bz(nX+2,:,:)       
    !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:),Bx(:,:,nZ+2)
    !       -2                             By(nX+2,:,:),By(:,:,nZ+2)
    do k=1-iGCN,nZ+iGCN-1; do j=1-iGCN,nY+iGCN; do i=1-iGCN,nX+iGCN-1
       Magnetic_GD(i,j,k,y_) = Magnetic_GD(i,j,k,y_) - &
            SpeedOfLightHalf_D(z_)*&
            (E_GD(i,  j,  k+1,x_) - E_GD(i,j,k,x_)) + &
            SpeedOfLightHalf_D(x_)*&
            (E_GD(i+1,j,  k  ,z_) - E_GD(i,j,k,z_))
    end do; end do; end do

    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1
    !                                                Ez (:,:,nZ+2)    
    !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:)
    !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:)
    !       _!_!_!                x!x!_!_  Bz(:,nY+2,:),Bz(nX+2,:,:)       
    !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:),Bx(:,:,nZ+2)
    !       -2                             By(nX+2,:,:),By(:,:,nZ+2)
    do k=1-iGCN,nZ+iGCN; do j=1-iGCN,nY+iGCN-1; do i=1-iGCN,nX+iGCN-1
       Magnetic_GD(i,j,k,z_) = Magnetic_GD(i,j,k,z_) - &
            SpeedOfLightHalf_D(x_)*&
            (E_GD(i+1,j,  k  ,y_) - E_GD(i,j,k,y_)) + &
            SpeedOfLightHalf_D(y_)*&
            (E_GD(i,  j+1,k  ,x_) - E_GD(i,j,k,x_))
    end do; end do; end do
  end subroutine update_magnetic
  !-----------------------------------------------------------------!
  subroutine update_e
    !Add current
    E_GD(0:nX,1:nY,1:nZ,x_) = E_GD(0:nX,1:nY,1:nZ,x_) - &
         Counter_GD(0:nX,1:nY,1:nZ,x_)
    E_GD(1:nX,0:nY,1:nZ,y_) = E_GD(1:nX,0:nY,1:nZ,y_) - &
         Counter_GD(1:nX,0:nY,1:nZ,y_)
    E_GD(1:nX,1:nY,0:nZ,z_) = E_GD(1:nX,1:nY,0:nZ,z_) - &
         Counter_GD(1:nX,1:nY,0:nZ,z_)

    if(UseVectorPotential)then
       call get_b_from_a()!Put the magnetic field to Counter_GD
       call use_magnetic_field_from(Counter_GD)
    else
       call use_magnetic_field_from(Magnetic_GD)
    end if
  contains
    subroutine use_magnetic_field_from(BIn_GD)
      real,dimension(&
           1-iGCN:nX+iGCN,&
           1-iGCN:nY+iGCN,&
           1-iGCN:nZ+iGCN,3),intent(in)::&
           BIn_GD
      integer::i,j,k
      integer :: i1_D(nDim)
      !---------------
      i1_D = 0
      !\
      ! Along the directions at which the periodic boundary condition is used
      ! only the electric field at the faces inside the computation domain
      ! are calculated here, all the outer face values are filled in with the
      ! periodic boundary conditions. Along other directions the computational
      ! domain is extended by iGCN-1 (usually 1 ) layer of the gostcells and all
      ! the faces inside the extended domain are filled in here
      !/
      where(.not.IsPeriodicField_D)i1_D=1
      do k=1-i1_D(z_),nZ+i1_D(z_)
         do j=1-i1_D(y_),nY+i1_D(y_)
            do i=0-i1_D(x_),nX+i1_D(x_)
               E_GD(i,j,k,x_) = E_GD(i,j,k,x_) + &
                    SpeedOfLight_D(y_)*&
                             (BIn_GD(i,j,k,z_) - BIn_GD(i,  j-1,k  ,z_)) - &
                    SpeedOfLight_D(z_)*&
                             (BIn_GD(i,j,k,y_) - BIn_GD(i,  j,  k-1,y_))
            end do
         end do
      end do
      do k=1-i1_D(z_),nZ+i1_D(z_)
         do j=0-i1_D(y_),nY+i1_D(y_);
            do i=1-i1_D(x_),nX+i1_D(x_)
               E_GD(i,j,k,y_) = E_GD(i,j,k,y_) + &
                    SpeedOfLight_D(z_)*&
                             (BIn_GD(i,j,k,x_) - BIn_GD(i,  j,  k-1,x_))-&
                    SpeedOfLight_D(x_)*&
                             (BIn_GD(i,j,k,z_) - BIn_GD(i-1,j,  k  ,z_))
            end do
         end do
      end do
      do k=0-i1_D(z_),nZ+i1_D(z_)
         do j=1-i1_D(y_),nY+i1_D(y_)
            do i=1-i1_D(x_),nX+i1_D(x_)
               E_GD(i,j,k,z_) = E_GD(i,j,k,z_)+&
                    SpeedOfLight_D(x_)*&
                             (BIn_GD(i,j,k,y_) - BIn_GD(i-1,j,  k  ,y_))-&
                    SpeedOfLight_D(y_)*& 
                             (BIn_GD(i,j,k,x_) - BIn_GD(i,  j-1,k  ,x_))
            end do
         end do
      end do
    end subroutine use_magnetic_field_from
  end subroutine update_e
  !=======================
  subroutine field_bc
    use PIC_ModLaserBeam, ONLY:laser_beam 
    integer :: i,j,k
    real    :: x, y, z
    integer :: i1_D(nDim), i0_D(nDim)
    !---------------
    i1_D = 0
    !\
    ! Along the directions at which the periodic boundary condition is used
    ! only the electric field at the faces inside the computation domain
    ! are calculated here, all the outer face values are filled in with the
    ! periodic boundary conditions. Along other directions the computational
    ! domain is extended by iGCN-1 (usually 1 ) layer of the gostcells and all
    ! the faces inside the extended domain are filled in here
    !/
    where(.not.IsPeriodicField_D)i1_D=1
      
    i0_D=0; where(IsPeriodicField_D)i0_D=0


    if(IsPeriodicField_D(x_))goto 100
    ! Face X < 0    
    i=1-iGCN
    do k=1-iGCN*i0_D(z_),nZ+iGCN*i0_D(z_)
       do j=0-i1_D(y_),nY+i1_D(y_);
          E_GD(i,j,k,y_) = E_GD(i,j,k,y_)*(1 - SpeedOfLight_D(x_)) + &
               SpeedOfLight_D(x_)*E_GD(i+1,j,k,y_)
          if(TypeFieldBC_S(1)=='laserbeam')&
               call laser_beam(iDir=y_,       &
               x=(i-0.50)*Dx_D(x_),&
               y=j       *Dx_D(y_),&
               z=(k-0.50)*Dx_D(z_),&
               EField=E_GD(i,j,k,y_) )


       end do
    end do

    do k=0-i1_D(z_),nZ+i1_D(z_)
       do j=1-iGCN*i0_D(y_),nY+iGCN*i0_D(y_)
          E_GD(i,j,k,z_) = E_GD(i,j,k,z_)*(1 - SpeedOfLight_D(x_))+&
               SpeedOfLight_D(x_)*E_GD(i+1,j,k,z_)
          if(TypeFieldBC_S(1)=='laserbeam')&
               call laser_beam(iDir=z_,       &
               x=(i-0.50)*Dx_D(x_),&
               y=(j-0.50)*Dx_D(y_),&
               z=k       *Dx_D(z_),&
               EField=E_GD(i,j,k,z_) )
       end do
    end do


    !face X > dx.nX    

    do k=1-iGCN*i0_D(z_),nZ+iGCN*i0_D(z_)
       do j=0-i1_D(y_),nY+i1_D(y_);
          i=nX+iGCN
          E_GD(i,j,k,y_) = E_GD(i,j,k,y_)*(1 - SpeedOfLight_D(x_)) + &
               SpeedOfLight_D(x_)* E_GD(i-1,j,k,y_)
          if(TypeFieldBC_S(2)=='laserbeam')&
               call laser_beam(iDir=y_,       &
               x=(i-0.50)*Dx_D(x_),&
               y=j       *Dx_D(y_),&
               z=(k-0.50)*Dx_D(z_),&
               EField=E_GD(i,j,k,y_) )
       end do
    end do
    do k=0-i1_D(z_),nZ+i1_D(z_)
       do j=1-iGCN*i0_D(y_),nY+iGCN*i0_D(y_)
          i=nX+iGCN
          E_GD(i,j,k,z_) = E_GD(i,j,k,z_)*(1 - SpeedOfLight_D(x_)) + &
               SpeedOfLight_D(x_)*E_GD(i-1,j,k,z_)
          if(TypeFieldBC_S(2)=='laserbeam')&
               call laser_beam(iDir=z_,       &
               x=(i-0.50)*Dx_D(x_),&
               y=(j-0.50)*Dx_D(y_),&
               z=k       *Dx_D(z_),&
               EField=E_GD(i,j,k,z_) )
       end do
    end do
100 continue

    if(IsPeriodicField_D(y_))goto 200
    !\
    ! Face Y<0
    !/ 
    j=1-iGCN
    do k=1-iGCN*i0_D(z_),nZ+iGCN*i0_D(z_)
       do i=0-i1_D(x_),nX+i1_D(x_)
          E_GD(i,j,k,x_) = E_GD(i,j,k,x_)*(1 - SpeedOfLight_D(y_)) + &
               SpeedOfLight_D(y_)*E_GD(i,j+1,k,x_)
          if(TypeFieldBC_S(3)=='laserbeam')&
               call laser_beam(iDir=x_,       &
               x= i      *Dx_D(x_),&
               y=(j-0.50)*Dx_D(y_),&
               z=(k-0.50)*Dx_D(z_),&
               EField=E_GD(i,j,k,x_) )
       end do
    end do

    do k=0-i1_D(z_),nZ+i1_D(z_)
       do i=1-iGCN*i0_D(x_),nX+iGCN*i0_D(x_)
          E_GD(i,j,k,z_) = E_GD(i,j,k,z_)*(1 - SpeedOfLight_D(y_)) +&
               SpeedOfLight_D(y_)*E_GD(i,j+1,k,z_)
          if(TypeFieldBC_S(3)=='laserbeam')&
               call laser_beam(iDir=z_,       &
               x=(i-0.50)*Dx_D(x_),&
               y=(j-0.50)*Dx_D(y_),&
               z=k       *Dx_D(z_),&
               EField=E_GD(i,j,k,z_) )
       end do
    end do
    !\
    ! Face Y>nY.dy
    !/    
    j=nY+iGCN
    do k=1-iGCN*i0_D(z_),nZ+iGCN*i0_D(z_)
       do i=0-i1_D(x_),nX+i1_D(x_)
          E_GD(i,j,k,x_) = E_GD(i,j,k,x_)*(1 - SpeedOfLight_D(y_)) + &
               SpeedOfLight_D(y_)*E_GD(i,j-1,k,x_)
          if(TypeFieldBC_S(4)=='laserbeam')&
               call laser_beam(iDir=x_,       &
               x= i      *Dx_D(x_),&
               y=(j-0.50)*Dx_D(y_),&
               z=(k-0.50)*Dx_D(z_),&
               EField=E_GD(i,j,k,x_) )
       end do
    end do

    do k=0-i1_D(z_),nZ+i1_D(z_)
       do i=1-iGCN*i0_D(x_),nX+iGCN*i0_D(x_)
          E_GD(i,j,k,z_) = E_GD(i,j,k,z_)*(1 - SpeedOfLight_D(y_)) + &
               SpeedOfLight_D(y_)*E_GD(i,j-1,k,z_)
          if(TypeFieldBC_S(4)=='laserbeam')&
               call laser_beam(iDir=z_,       &
               x=(i-0.50)*Dx_D(x_),&
               y=(j-0.50)*Dx_D(y_),&
               z=k       *Dx_D(z_),&
               EField=E_GD(i,j,k,z_) )
       end do
    end do

200 continue
    if(IsPeriodicField_D(z_))return
    !\
    ! Face Z<0
    !/
    k=1-iGCN
    do j=1-iGCN*i0_D(y_),nY+iGCN*i0_D(y_)
       do i=0-i1_D(x_),nX+i1_D(x_)
          E_GD(i,j,k,x_) = E_GD(i,j,k,x_)*(1 - SpeedOfLight_D(z_)) + &
               SpeedOfLight_D(z_)*E_GD(i,j,k+1,x_)
          if(TypeFieldBC_S(5)=='laserbeam')&
               call laser_beam(iDir=x_,       &
               x= i      *Dx_D(x_),&
               y=(j-0.50)*Dx_D(y_),&
               z=(k-0.50)*Dx_D(z_),&
               EField=E_GD(i,j,k,x_) )
       end do
    end do
    do j=0-i1_D(y_),nY+i1_D(y_)
       do i=1-iGCN*i0_D(x_),nX+iGCN*i0_D(x_)
          E_GD(i,j,k,y_) = E_GD(i,j,k,y_)*(1 - SpeedOfLight_D(z_)) + &
               SpeedOfLight_D(z_)*E_GD(i,j,k+1,y_)
          if(TypeFieldBC_S(5)=='laserbeam')&
               call laser_beam(iDir=y_,       &
               x=(i-0.50)*Dx_D(x_),&
               y=j       *Dx_D(y_),&
               z=(k-0.50)*Dx_D(z_),&
               EField=E_GD(i,j,k,y_) )
       end do
    end do

    !\
    ! Face Z>nZ.dz
    !/
    k=nZ+iGCN

    do j=1-iGCN*i0_D(y_),nY+iGCN*i0_D(y_)
       do i=0-i1_D(x_),nX+i1_D(x_)
          E_GD(i,j,k,x_) = E_GD(i,j,k,x_)*(1 - SpeedOfLight_D(z_)) + &
               SpeedOfLight_D(z_)*E_GD(i,j,k-1,x_)
          if(TypeFieldBC_S(6)=='laserbeam')&
               call laser_beam(iDir=x_,       &
               x= i      *Dx_D(x_),&
               y=(j-0.50)*Dx_D(y_),&
               z=(k-0.50)*Dx_D(z_),&
               EField=E_GD(i,j,k,x_) )
       end do
    end do
    do j=0-i1_D(y_),nY+i1_D(y_)
       do i=1-iGCN*i0_D(x_),nX+iGCN*i0_D(x_)
          E_GD(i,j,k,y_) = E_GD(i,j,k,y_)*(1 - SpeedOfLight_D(z_)) + &
               SpeedOfLight_D(z_)*E_GD(i,j,k-1,y_)
          if(TypeFieldBC_S(6)=='laserbeam')&
               call laser_beam(iDir=y_,       &
               x=(i-0.50)*Dx_D(x_),&
               y=j       *Dx_D(y_),&
               z=(k-0.50)*Dx_D(z_),&
               EField=E_GD(i,j,k,y_) )
       end do
    end do
  end subroutine field_bc
  !=======================
  subroutine get_rho_max(RhoMax,Coord_D)
    real,intent(out)::RhoMax
    real,dimension(nDim),optional,intent(out)::Coord_D
    !------------------
    RhoMax = maxval(rho_G(1:nX,1:nY,1:nZ))
    if(present(Coord_D))Coord_D=real(maxloc(rho_G(1:nX,1:nY,1:nZ))) - cHalf
  end subroutine get_rho_max
  !---------------------------------------------------------------------!
  subroutine get_max_intensity(EnergyMax,Coord_D)
    real,intent(out)::EnergyMax
    real,dimension(nDim),optional,intent(out)::Coord_D
    integer::i,j,k
    !-------------
    rho_G=cZero
    if(UseVectorPotential)then
       call get_b_from_a()
       do k=1,nZ; do j=1,nY; do i=1,nX
          rho_G(i,j,k)=0.1250*(&
               sum(Counter_GD(i,j-1:j,k-1:k,x_)**2)+&
               sum(Counter_GD(i-1:i,j,k-1:k,y_)**2)+&
               sum(Counter_GD(i-1:i,j-1:j,k,z_)**2))+&
               0.250 * (&
               sum(E_GD(i-1:i,j,k,x_)**2)+&
               sum(E_GD(i,j-1:j,k,y_)**2)+&
               sum(E_GD(i,j,k-1:k,z_)**2))
       end do; end do; end do
    else
       do k=1,nZ; do j=1,nY; do i=1,nX
          rho_G(i,j,k)=0.1250*(&
               sum(Magnetic_GD(i,j-1:j,k-1:k,x_)**2)+&
               sum(Magnetic_GD(i-1:i,j,k-1:k,y_)**2)+&
               sum(Magnetic_GD(i-1:i,j-1:j,k,z_)**2))+&
               0.250*(&
               sum(E_GD(i-1:i,j,k,x_)**2)+&
               sum(E_GD(i,j-1:j,k,y_)**2)+&
               sum(E_GD(i,j,k-1:k,z_)**2))
       end do; end do; end do
    end if
    call get_rho_max(EnergyMax,Coord_D)
  end subroutine get_max_intensity
  !====================
  subroutine get_field_energy(Energy_V)
    use PIC_ModMain, ONLY: CellVolume
    real,intent(out)::Energy_V(6)
    integer::i,j,k
    !-------------
    Energy_V = cZero
    if(UseVectorPotential)then
       call get_b_from_a()
       do k=1,nZ; do j=1,nY; do i=1,nX
          Energy_V(1) = Energy_V(1) + 0.1250*&
               sum((Counter_GD(i,j-1:j,k-1:k,x_) - B0_D(x_))**2)
          
          Energy_V(2) = Energy_V(2) + 0.1250*&
               sum((Counter_GD(i-1:i,j,k-1:k,y_) - B0_D(y_))**2)
          
          Energy_V(3) = Energy_V(3) + 0.1250*&
               sum((Counter_GD(i-1:i,j-1:j,k,z_) - B0_D(z_))**2)
          
          
          Energy_V(4) = Energy_V(4) + 0.250*&
               sum(E_GD(i-1:i,j,k,x_)**2)
          
          Energy_V(5) = Energy_V(5) + 0.250*&
               sum(E_GD(i,j-1:j,k,y_)**2)
          
          Energy_V(6) = Energy_V(6) + 0.250*&
               sum(E_GD(i,j,k-1:k,z_)**2)
       end do; end do; end do
    else
       do k=1,nZ; do j=1,nY; do i=1,nX
          Energy_V(1) = Energy_V(1) + 0.1250*&
               sum((Magnetic_GD(i,j-1:j,k-1:k,x_) - B0_D(x_))**2)
          
          Energy_V(2) = Energy_V(2) + 0.1250*&
               sum((Magnetic_GD(i-1:i,j,k-1:k,y_) - B0_D(y_))**2)
          
          Energy_V(3) = Energy_V(3) + 0.1250*&
               sum((Magnetic_GD(i-1:i,j-1:j,k,z_) - B0_D(z_))**2)
          
          
          Energy_V(4) = Energy_V(4) + 0.250*&
               sum(E_GD(i-1:i,j,k,x_)**2)
          
          Energy_V(5) = Energy_V(5) + 0.250*&
               sum(E_GD(i,j-1:j,k,y_)**2)
          
          Energy_V(6) = Energy_V(6) + 0.250*&
               sum(E_GD(i,j,k-1:k,z_)**2)
       end do; end do; end do
    end if
    Energy_V = (CellVolume/(4.0 * cPi)) * Energy_V
  end subroutine get_field_energy
  !===============================
  !ONLY FOR PERIODIC BC.
  !================
  subroutine density_bc_periodic
    integer :: i, j, k, iInner, jInner, kInner
    !-------------------------------!
    do k=1-iGCN,nZ+iGCN
       kInner = k - nZ*floor( (k - 0.50) /nZ )
       do j=1-iGCN,nY+iGCN
          jInner = j - nY*floor( (j - 0.50) /nY )
          do i=1-iGCN,nX+iGCN
             iInner = i - nX*floor( (i - 0.50) /nX )
             if(i==iInner.and.j==jInner.and.&
                  k==kInner)CYCLE
             Rho_G(iInner,jInner,kInner) = &
                  Rho_G(iInner,jInner,kInner) + &
                  Rho_G(i,j,k)
             Rho_G(i,j,k) = 0.0
          end do
       end do
    end do
  end subroutine density_bc_periodic
  !=================================
  subroutine current_bc_periodic
    real :: Current_F(0:nX,0:nY,0:nZ)
    integer :: i, j, k, iInner, jInner, kInner
    !-------------------------------!
    !++++++++++++x_ currents++++++++++++++++=
    Current_F = 0.0
    do k=1-iGCN,nZ+iGCN
       kInner = k - nZ*floor( (k - 0.50) /nZ )
       do j=1-iGCN,nY+iGCN
          jInner = j - nY*floor( (j - 0.50) /nY )
          !Select x_ faces.
          ! -1
          do i=1-iGCN, -1
             iInner = i + nX
             Current_F(iInner,jInner,kInner) = &
                  Current_F(iInner,jInner,kInner) + &
                  Counter_GD(i,j,k, x_)
          end do
          !  0
          Current_F( 0,jInner,kInner) = &
               Current_F( 0,jInner,kInner) + &
               Counter_GD( 0,j,k, x_)      + &
               Counter_GD(nX,j,k, x_)
          ! 1:nX-1
          do i=1,nX-1
             iInner = i
             Current_F(iInner,jInner,kInner) = &
                  Current_F(iInner,jInner,kInner) + &
                  Counter_GD(i,j,k, x_)
          end do
          !nX
          Current_F(nX,jInner,kInner) = &
               Current_F(nX,jInner,kInner) + &
               Counter_GD( 0,j,k, x_)      + &
               Counter_GD(nX,j,k, x_)
          !nX+1
          do i=nX+1, nX+iGCN-1
             iInner = i - nX
             Current_F(iInner,jInner,kInner) = &
                  Current_F(iInner,jInner,kInner) + &
                  Counter_GD(i,j,k, x_)
          end do
       end do
    end do
    Counter_GD(:,:,:,x_) = 0.0
    Counter_GD(0:nX,1:nY,1:nZ,x_) = Current_F(0:nX,1:nY,1:nZ)
    
    !++++++++++y_ currents+++++++++++++++
    
    Current_F = 0.0
    do k=1-iGCN,nZ+iGCN
       kInner = k - nZ*floor( (k - 0.50) /nZ )
       !Select y_ faces
       !-1
       do j=1-iGCN, -1
          jInner = j + nY
          do i=1-iGCN,nX+iGCN
             iInner = i - nX*floor( (i - 0.50) /nX )
             Current_F(iInner,jInner,kInner) = &
                  Current_F(iInner,jInner,kInner) + &
                  Counter_GD(i,j,k, y_)
          end do
       end do
       !0
       do i=1-iGCN,nX+iGCN
          iInner = i - nX*floor( (i - 0.50) /nX )
          Current_F( iInner,0,kInner) = &
               Current_F( iInner,0,kInner) + &
               Counter_GD( i, 0,k, y_)      + &
               Counter_GD( i,nY,k, y_)
       end do
       !1:nY-1
       do j=1,nY-1
          jInner = j 
          do i=1-iGCN,nX+iGCN
             iInner = i - nX*floor( (i - 0.50) /nX )
             Current_F(iInner,jInner,kInner) = &
                  Current_F(iInner,jInner,kInner) + &
                  Counter_GD(i,j,k, y_)
          end do
       end do
       !nY
       do i=1-iGCN,nX+iGCN
          iInner = i - nX*floor( (i - 0.50) /nX )
          Current_F( iInner,nY,kInner) = &
               Current_F( iInner,nY,kInner) + &
               Counter_GD( i, 0,k, y_)      + &
               Counter_GD( i,nY,k, y_)
       end do
       !-1
       do j=nY+1,nY+iGCN-1
          jInner = j - nY
          do i=1-iGCN,nX+iGCN
             iInner = i - nX*floor( (i - 0.50) /nX )
             Current_F(iInner,jInner,kInner) = &
                  Current_F(iInner,jInner,kInner) + &
                  Counter_GD(i,j,k, y_)
          end do
       end do
    end do
    Counter_GD(:,:,:,y_) = 0.0
    Counter_GD(1:nX,0:nY,1:nZ,y_) = Current_F(1:nX,0:nY,1:nZ)
    
    !++++++++++z_ currents+++++++++++++++
    
    Current_F = 0.0
    !Select z_ faces
    !-1
    do k=1-iGCN,-1
       kInner = k + nZ
       do j=1-iGCN,nY+iGCN
          jInner = j - nY*floor( (j - 0.50) /nY ) 
          do i=1-iGCN,nX+iGCN
             iInner = i - nX*floor( (i - 0.50) /nX )
             Current_F(iInner,jInner,kInner) = &
                  Current_F(iInner,jInner,kInner) + &
                  Counter_GD(i,j,k, z_)
          end do
       end do
    end do
    !0
    do j=1-iGCN,nY+iGCN
       jInner = j - nY*floor( (j - 0.50) /nY )
       do i=1-iGCN,nX+iGCN
          iInner = i - nX*floor( (i - 0.50) /nX )
          Current_F( iInner,jInner,0) = &
               Current_F( iInner,jInner,0) + &
               Counter_GD( i, j, 0, z_)      + &
               Counter_GD( i, j,nZ, z_)
       end do
    end do
    !1:nZ-1
    do k=1,nZ-1
       kInner = k
       do j=1-iGCN,nY+iGCN
          jInner = j - nY*floor( (j - 0.50) /nY ) 
          do i=1-iGCN,nX+iGCN
             iInner = i - nX*floor( (i - 0.50) /nX )
             Current_F(iInner,jInner,kInner) = &
                  Current_F(iInner,jInner,kInner) + &
                  Counter_GD(i,j,k, z_)
          end do
       end do
    end do
    !mZ
    do j=1-iGCN,nY+iGCN
       jInner = j - nY*floor( (j - 0.50) /nY )
       do i=1-iGCN,nX+iGCN
          iInner = i - nX*floor( (i - 0.50) /nX )
          Current_F( iInner,jInner,nZ) = &
               Current_F( iInner,jInner,nZ) + &
               Counter_GD( i, j, 0, z_)      + &
               Counter_GD( i, j,nZ, z_)
       end do
    end do
    !nZ+1
    do k=nZ+1,nZ+iGCN,-1
       kInner = k - nZ
       do j=1-iGCN,nY+iGCN
          jInner = j - nY*floor( (j - 0.50) /nY ) 
          do i=1-iGCN,nX+iGCN
             iInner = i - nX*floor( (i - 0.50) /nX )
             Current_F(iInner,jInner,kInner) = &
                  Current_F(iInner,jInner,kInner) + &
                  Counter_GD(i,j,k, z_)
          end do
       end do
    end do
    
    Counter_GD(:,:,:,z_) = 0.0
    Counter_GD(1:nX,1:nY,0:nZ,z_) = Current_F(1:nX,1:nY,0:nZ)
    
  end subroutine current_bc_periodic

    !======================
  subroutine field_bc_periodic
    use PIC_ModMain, ONLY: IsPeriodicField_D
    integer :: i, j, k, iInner, jInner, kInner
    integer :: nXyzIfPeriodic_D(nDim)
    !-------------------------------!
    nXyzIfPeriodic_D = 0
    where(IsPeriodicField_D)nXyzIfPeriodic_D=nCell_D
    !++++++++++++x_ field++++++++++++++++=
    do k=1-iGCN,nZ+iGCN
       kInner = k - nXyzIfPeriodic_D(z_)*floor( (k - 0.50) /nZ )
       do j=1-iGCN,nY+iGCN
          jInner = j -  nXyzIfPeriodic_D(y_)*floor( (j - 0.50) /nY )
          !\
          !Select x_ faces
          !/
          !-1
          do i=1-iGCN, -1
             iInner = i +  nXyzIfPeriodic_D(x_)
             E_GD(i,j,k,x_) = &
                  E_GD(iInner,jInner,kInner,x_)
          end do
          !nX+1
          do i=nX+1, nX+iGCN-1
             iInner = i -  nXyzIfPeriodic_D(x_)
             E_GD(i,j,k,x_) = &
                  E_GD(iInner,jInner,kInner,x_)
          end do
          if(k==kInner.and.j==jInner)CYCLE
          !0:nX
          E_GD(0:nX,j,k,x_) = &
               E_GD(0:nX,jInner,kInner,x_)            
       end do
    end do


    !++++++++++y_ field+++++++++++++++
    do k=1-iGCN,nZ+iGCN
       kInner = k -  nXyzIfPeriodic_D(z_)*floor( (k - 0.50) /nZ )
       !\
       !Select y_ face
       !/
       !-1
       do j=1-iGCN, -1
          jInner = j +  nXyzIfPeriodic_D(y_)
          do i=1-iGCN,nX+iGCN
             iInner = i -  nXyzIfPeriodic_D(x_)*floor( (i - 0.50) /nX )
             E_GD(i,j,k,y_) = &
                  E_GD(iInner,jInner,kInner,y_) 
          end do
       end do
       !nY+1
       do j=nY+1,nY+iGCN-1
          jInner = j -  nXyzIfPeriodic_D(y_)
          do i=1-iGCN,nX+iGCN
             iInner = i -  nXyzIfPeriodic_D(x_)*floor( (i - 0.50) /nX )
             E_GD(i,j,k,y_) = &
                  E_GD(iInner,jInner,kInner,y_)
          end do
       end do
       !0:nY
       do j=0,nY
          jInner = j 
          do i=1-iGCN,nX+iGCN
             iInner = i -  nXyzIfPeriodic_D(x_)*floor( (i - 0.50) /nX )
             if(iInner==i.and.k==kInner)CYCLE
             E_GD(i,j,k,y_) = &
                  E_GD(iInner,jInner,kInner,y_)
          end do
       end do
    end do


    !++++++++++z_ field+++++++++++++++
    !Select z_ faces
    !-1
    do k=1-iGCN,-1
       kInner = k +  nXyzIfPeriodic_D(z_)
       do j=1-iGCN,nY+iGCN
          jInner = j -  nXyzIfPeriodic_D(y_)*floor( (j - 0.50) /nY ) 
          do i=1-iGCN,nX+iGCN
             iInner = i -  nXyzIfPeriodic_D(x_)*floor( (i - 0.50) /nX )
             E_GD(i,j,k,z_) = &
                  E_GD(iInner,jInner,kInner,z_)
          end do
       end do
    end do
    !nZ+1
    do k=nZ+1,nZ+iGCN,-1
       kInner = k -  nXyzIfPeriodic_D(z_)
       do j=1-iGCN,nY+iGCN
          jInner = j -  nXyzIfPeriodic_D(y_)*floor( (j - 0.50) /nY ) 
          do i=1-iGCN,nX+iGCN
             iInner = i -  nXyzIfPeriodic_D(x_)*floor( (i - 0.50) /nX )
             E_GD(i,j,k,z_) = &
                  E_GD(iInner,jInner,kInner,z_)
          end do
       end do
    end do
    !0:nZ
    do k=0,nZ
       kInner = k
       do j=1-iGCN,nY+iGCN
          jInner = j -  nXyzIfPeriodic_D(y_)*floor( (j - 0.50) /nY ) 
          do i=1-iGCN,nX+iGCN
             iInner = i -  nXyzIfPeriodic_D(x_)*floor( (i - 0.50) /nX )
             if(i==iInner.and.j==jInner)CYCLE
             E_GD(i,j,k,z_) = &
                  E_GD(iInner,jInner,kInner,z_)
          end do
       end do
    end do
  end subroutine field_bc_periodic
end module PIC_ModField
