!^CFG COPYRIGHT UofM
!--------------------------------------------------------------!

module PIC_ModField
  use PIC_ModSize, ONLY: nX, nY, nZ, x_, y_, z_, &
       nDim, nCell_D, MaxBlock
  use PC_BATL_size,ONLY: MaxDim
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
  !                                                Ez (:,:,nZ+2,1)    
  !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:)
  !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:,1)
  !       _!_!_!                x!x!_!_  Bz(:,nY+2,:,1),Bz(nX+2,:,:,1)       
  !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:,1),Bx(:,:,nZ+2,1)
  !       -2                             By(nX+2,:,:,1),By(:,:,nZ+2,1)
  !real,dimension(&
  !     1-iGCN:nX+iGCN,&
  !     1-iGCN:nY+iGCN,&
  !     1-iGCN:nZ+iGCN,3,MaxBlock)::&
  !     E_GDB        = 0.0,& !This is the electric field
  !     B_GDB = 0.0,& !vector-potential if used, magnetic field otherwise
  !     Counter_GDB  = 0.0   !Counter for electric current
  !     V_GDB        = 0.0   !Cell-centered velocity
  real,allocatable:: E_GDB(:,:,:,:,:)
  real,allocatable:: B_GDB(:,:,:,:,:)
  real,allocatable:: Counter_GDB(:,:,:,:,:)

  real,allocatable:: Rho_GB(:,:,:,:)
  real,allocatable:: V_GDB(:,:,:,:,:)
  real :: B0_D(3) = 0.0

  !Methods

  public::get_b_from_a   !Transforms the vector potential to magn. field
  public::update_magnetic!Updates the magnetic field, or vector potential
  public::update_e       !Updates the electric field
  public::get_rho_max     
contains
  !=======================
  subroutine allocate_fields
    
    allocate(E_GDB(1-iGCN:nX+iGCN, 1-iGCN:nY+iGCN,&
         1-iGCN:nZ+iGCN, MaxDim, MaxBlock)); E_GDB = 0.0
    allocate(B_GDB(1-iGCN:nX+iGCN, 1-iGCN:nY+iGCN,&
         1-iGCN:nZ+iGCN, MaxDim, MaxBlock)); B_GDB = 0.0 
    allocate(Counter_GDB(1-iGCN:nX+iGCN, 1-iGCN:nY+iGCN,&
         1-iGCN:nZ+iGCN, MaxDim, MaxBlock)); Counter_GDB = 0.0 
    allocate(Rho_GB(1-iGCN:nX+iGCN, 1-iGCN:nY+iGCN,&
         1-iGCN:nZ+iGCN, MaxBlock)); Rho_GB = 0.0
    allocate(V_GDB(MaxDim, 1-iGCN:nX+iGCN,1-iGCN:nY+iGCN,&
         1-iGCN:nZ+iGCN, MaxBlock)); V_GDB = 0.0

  end subroutine allocate_fields
  !================
  subroutine add_e
    use ModReadParam, ONLY: read_var
    real:: E_D(3)
    integer:: i,j,k,iDim, iBlock
    !------------
    call read_var('Ex',E_D(1))
    call read_var('Ey',E_D(2))
    call read_var('Ez',E_D(3))
    do iBlock = 1, MaxBlock
    do iDim = 1,3
       do k=1-iGCN,nZ+iGCN
          do j=1-iGCN,nY+iGCN
             do i=1-iGCN,nX+iGCN
                E_GDB(i,j,k,iDim,iBlock) = &
                     E_GDB(i,j,k,iDim,iBlock) + E_D(iDim)
             end do
          end do
       end do
    end do
    end do
  end subroutine add_e
  !===================
  subroutine add_b
    use ModReadParam, ONLY: read_var
    integer:: i,j,k,iDim, iBlock
    !------------
    call read_var('Bx',B0_D(1))
    call read_var('By',B0_D(2))
    call read_var('Bz',B0_D(3))
    do iBlock = 1, MaxBlock
    do iDim = 1,3
       do k=1-iGCN,nZ+iGCN
          do j=1-iGCN,nY+iGCN
             do i=1-iGCN,nX+iGCN
                B_GDB(i,j,k,iDim,iBlock) = &
                     B_GDB(i,j,k,iDim,iBlock) + B0_D(iDim)
             end do
          end do
       end do
    end do
    end do
  end subroutine add_b
  !===================
  subroutine get_b_from_a(iBlock)
    integer, intent(in)::iBlock
    integer::i,j,k
    !Applied only if UseVectorPotential
    !The vector potential is in B_GDB
    !Calculates the magnetic field 

    !The whole magnetic field is calculated and 
    !the result is saved to Counter_GDB
    !---------------------------------
    Counter_GDB = 0.0

    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1
    !                                                Ez (:,:,nZ+2,1)    
    !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:,1)
    !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:,1)
    !       _!_!_!                x!x!_!_  Bz(:,nY+2,:,1),Bz(nX+2,:,:,1)       
    !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:,1),Bx(:,:,nZ+2,1)
    !       -2                             By(nX+2,:,:,1),By(:,:,nZ+2,1)
    do k=1-iGCN,nZ+iGCN-1; do j=1-iGCN,nY+iGCN-1; do i=1-iGCN,nX+iGCN
       Counter_GDB(i,j,k,x_,iBlock) = &
            SpeedOfLight_D(y_)*&
            (B_GDB(i,j+1,k,z_,iBlock) - B_GDB(i,j,k,z_,iBlock))-&
            SpeedOfLight_D(z_)*&
            (B_GDB(i,j,k+1,y_,iBlock) - B_GDB(i,j,k,y_,iBlock))
    end do; end do; end do

    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1
    !                                                Ez (:,:,nZ+2,1)    
    !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:,1)
    !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:,1)
    !       _!_!_!                x!x!_!_  Bz(:,nY+2,:,1),Bz(nX+2,:,:,1)       
    !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:,1),Bx(:,:,nZ+2,1)
    !       -2                             By(nX+2,:,:,1),By(:,:,nZ+2,1)
    do k=1-iGCN,nZ+iGCN-1; do j=1-iGCN,nY+iGCN; do i=1-iGCN,nX+iGCN-1
       Counter_GDB(i,j,k,y_,iBlock) = &
            SpeedOfLight_D(z_)*&
            (B_GDB(i,j,k+1,x_,iBlock) - B_GDB(i,j,k,x_,iBlock))-&
            SpeedOfLight_D(x_)*&
            (B_GDB(i+1,j,k,z_,iBlock) - B_GDB(i,j,k,z_,iBlock))
    end do; end do; end do

    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1
    !                                                Ez (:,:,nZ+2,1)    
    !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:,1)
    !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:,1)
    !       _!_!_!                x!x!_!_  Bz(:,nY+2,:,1),Bz(nX+2,:,:,1)       
    !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:,1),Bx(:,:,nZ+2,1)
    !       -2                             By(nX+2,:,:,1),By(:,:,nZ+2,1)
    do k=1-iGCN,nZ+iGCN; do j=1-iGCN,nY+iGCN-1; do i=1-iGCN,nX+iGCN-1
       Counter_GDB(i,j,k,z_,iBlock) = &
            SpeedOfLight_D(x_)*&
            (B_GDB(i+1,j,k,y_,iBlock) - B_GDB(i,j,k,y_,iBlock))-&
            SpeedOfLight_D(y_)*& 
            (B_GDB(i,j+1,k,x_,iBlock) - B_GDB(i,j,k,x_,iBlock))
    end do; end do; end do
  end subroutine get_b_from_a
  !==========================

  subroutine update_magnetic
    integer::i,j,k,iBlock
    real,dimension(nDim):: SpeedOfLightHalf_D
    !--------------------------
    if(UseVectorPotential)then

       !Advance vector potential throught a half timestep
       !Array index is the coordinate of the gridpoint with +1/2 being
       !approximated as 1
       !                                                Ez (:,:,nZ+2,1)    
       !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:,1)
       !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:,1)
       !       _!_!_!                x!x!_!_  Bz(:,nY+2,:,1),Bz(nX+2,:,:,1)       
       !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:,1),Bx(:,:,nZ+2,1)
       !       -2                             By(nX+2,:,:,1),By(:,:,nZ+2,1)
       do iBlock = 1,MaxBlock
       B_GDB(     1-iGCN:nX+iGCN-1,:,:,x_,iBlock) = &
            B_GDB(1-iGCN:nX+iGCN-1,:,:,x_,iBlock) - &
            cHalf*E_GDB(1-iGCN:nX+iGCN-1,:,:,x_,iBlock)

       B_GDB(     :,1-iGCN:nY+iGCN-1,:,y_,iBlock) = &
            B_GDB(:,1-iGCN:nY+iGCN-1,:,y_,iBlock) - &
            cHalf*E_GDB(:,1-iGCN:nY+iGCN-1,:,y_,iBlock)

       B_GDB(     :,:,1-iGCN:nZ+iGCN-1,z_,iBlock) = &
            B_GDB(:,:,1-iGCN:nZ+iGCN-1,z_,iBlock) - &
            cHalf*E_GDB(:,:,1-iGCN:nZ+iGCN-1,z_,iBlock)
       end do
       return
    end if

    SpeedOfLightHalf_D = &
         SpeedOfLight_D*cHalf

    !Advance vector potential throught a half timestep
    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1
    !                                                Ez (:,:,nZ+2,1)    
    !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:,1)
    !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:,1)
    !       _!_!_!                x!x!_!_  Bz(:,nY+2,:,1),Bz(nX+2,:,:,1)       
    !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:,1),Bx(:,:,nZ+2,1)
    !       -2                             By(nX+2,:,:,1),By(:,:,nZ+2,1)
    do iBlock = 1, MaxBlock
    do k=1-iGCN,nZ+iGCN-1; do j=1-iGCN,nY+iGCN-1; do i=1-iGCN,nX+iGCN
       B_GDB(i,j,k,x_,iBlock) = B_GDB(i,j,k,x_,iBlock) - &
            SpeedOfLightHalf_D(y_)*&
            (E_GDB(i,  j+1,k  ,z_,iBlock) - E_GDB(i,j,k,z_,iBlock)) + &
            SpeedOfLightHalf_D(z_)*&
            (E_GDB(i,  j,  k+1,y_,iBlock) - E_GDB(i,j,k,y_,iBlock))
    end do; end do; end do

    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1
    !                                                Ez (:,:,nZ+2)    
    !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:,1)
    !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:,1)
    !       _!_!_!                x!x!_!_  Bz(:,nY+2,:,1),Bz(nX+2,:,:,1)       
    !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:,1),Bx(:,:,nZ+2,1)
    !       -2                             By(nX+2,:,:,1),By(:,:,nZ+2,1)
    do k=1-iGCN,nZ+iGCN-1; do j=1-iGCN,nY+iGCN; do i=1-iGCN,nX+iGCN-1
       B_GDB(i,j,k,y_,iBlock) = B_GDB(i,j,k,y_,iBlock) - &
            SpeedOfLightHalf_D(z_)*&
            (E_GDB(i,  j,  k+1,x_,iBlock) - E_GDB(i,j,k,x_,iBlock)) + &
            SpeedOfLightHalf_D(x_)*&
            (E_GDB(i+1,j,  k  ,z_,iBlock) - E_GDB(i,j,k,z_,iBlock))
    end do; end do; end do

    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1
    !                                                Ez (:,:,nZ+2,1)    
    !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:,:,1)
    !       _!_!x!x               _!_!_!_            Ey (:,nY+2,:,1)
    !       _!_!_!                x!x!_!_  Bz(:,nY+2,:,1),Bz(nX+2,:,:,1)       
    !    -2  ! ! !                x!x!_!_  Bx(:,nY+2,:,1),Bx(:,:,nZ+2,1)
    !       -2                             By(nX+2,:,:,1),By(:,:,nZ+2,1)
    do k=1-iGCN,nZ+iGCN; do j=1-iGCN,nY+iGCN-1; do i=1-iGCN,nX+iGCN-1
       B_GDB(i,j,k,z_,iBlock) = B_GDB(i,j,k,z_,iBlock) - &
            SpeedOfLightHalf_D(x_)*&
            (E_GDB(i+1,j,  k  ,y_,iBlock) - E_GDB(i,j,k,y_,iBlock)) + &
            SpeedOfLightHalf_D(y_)*&
            (E_GDB(i,  j+1,k  ,x_,iBlock) - E_GDB(i,j,k,x_,iBlock))
    end do; end do; end do
    end do
  end subroutine update_magnetic
  !-----------------------------------------------------------------!
  subroutine update_e
    integer::iBlock
    !---------------
    do iBlock = 1, MaxBlock
    !Add current
    E_GDB(0:nX,1:nY,1:nZ,x_,iBlock) = E_GDB(0:nX,1:nY,1:nZ,x_,iBlock) - &
         Counter_GDB(0:nX,1:nY,1:nZ,x_,iBlock)
    E_GDB(1:nX,0:nY,1:nZ,y_,iBlock) = E_GDB(1:nX,0:nY,1:nZ,y_,iBlock) - &
         Counter_GDB(1:nX,0:nY,1:nZ,y_,iBlock)
    E_GDB(1:nX,1:nY,0:nZ,z_,iBlock) = E_GDB(1:nX,1:nY,0:nZ,z_,iBlock) - &
         Counter_GDB(1:nX,1:nY,0:nZ,z_,iBlock)

    if(UseVectorPotential)then
       call get_b_from_a(iBlock)!Put the magnetic field to Counter_GDB
       call use_magnetic_field_from(Counter_GDB(:,:,:,:,iBlock))
    else
       call use_magnetic_field_from(B_GDB(:,:,:,:,iBlock))
    end if
    end do
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
               E_GDB(i,j,k,x_,iBlock) = E_GDB(i,j,k,x_,iBlock) + &
                    SpeedOfLight_D(y_)*&
                             (BIn_GD(i,j,k,z_) - BIn_GD(i,j-1,k,z_)) - &
                    SpeedOfLight_D(z_)*&
                             (BIn_GD(i,j,k,y_) - BIn_GD(i,j,k-1,y_))
            end do
         end do
      end do
      do k=1-i1_D(z_),nZ+i1_D(z_)
         do j=0-i1_D(y_),nY+i1_D(y_);
            do i=1-i1_D(x_),nX+i1_D(x_)
               E_GDB(i,j,k,y_,iBlock) = E_GDB(i,j,k,y_,iBlock) + &
                    SpeedOfLight_D(z_)*&
                             (BIn_GD(i,j,k,x_) - BIn_GD(i,j,k-1,x_))-&
                    SpeedOfLight_D(x_)*&
                             (BIn_GD(i,j,k,z_) - BIn_GD(i-1,j,k,z_))
            end do
         end do
      end do
      do k=0-i1_D(z_),nZ+i1_D(z_)
         do j=1-i1_D(y_),nY+i1_D(y_)
            do i=1-i1_D(x_),nX+i1_D(x_)
               E_GDB(i,j,k,z_,iBlock) = E_GDB(i,j,k,z_,iBlock)+&
                    SpeedOfLight_D(x_)*&
                             (BIn_GD(i,j,k,y_) - BIn_GD(i-1,j,k,y_))-&
                    SpeedOfLight_D(y_)*& 
                             (BIn_GD(i,j,k,x_) - BIn_GD(i,j-1,k,x_))
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


    if(.not.IsPeriodicField_D(x_))then
       ! Face X < 0    
       i=1-iGCN
       do k=1-iGCN*i0_D(z_),nZ+iGCN*i0_D(z_)
          do j=0-i1_D(y_),nY+i1_D(y_);
             E_GDB(i,j,k,y_,1) = E_GDB(i,j,k,y_,1)*(1 - SpeedOfLight_D(x_)) + &
                  SpeedOfLight_D(x_)*E_GDB(i+1,j,k,y_,1)
             if(TypeFieldBC_S(1)=='laserbeam')&
                  call laser_beam(iDir=y_,       &
                  x=(i-0.50)*Dx_D(x_),&
                  y=j       *Dx_D(y_),&
                  z=(k-0.50)*Dx_D(z_),&
                  EField=E_GDB(i,j,k,y_,1) )


          end do
       end do

       do k=0-i1_D(z_),nZ+i1_D(z_)
          do j=1-iGCN*i0_D(y_),nY+iGCN*i0_D(y_)
             E_GDB(i,j,k,z_,1) = E_GDB(i,j,k,z_,1)*(1 - SpeedOfLight_D(x_))+&
                  SpeedOfLight_D(x_)*E_GDB(i+1,j,k,z_,1)
             if(TypeFieldBC_S(1)=='laserbeam')&
                  call laser_beam(iDir=z_,       &
                  x=(i-0.50)*Dx_D(x_),&
                  y=(j-0.50)*Dx_D(y_),&
                  z=k       *Dx_D(z_),&
                  EField=E_GDB(i,j,k,z_,1) )
          end do
       end do


       !face X > dx.nX    

       do k=1-iGCN*i0_D(z_),nZ+iGCN*i0_D(z_)
          do j=0-i1_D(y_),nY+i1_D(y_);
             i=nX+iGCN
             E_GDB(i,j,k,y_,1) = E_GDB(i,j,k,y_,1)*(1 - SpeedOfLight_D(x_)) + &
                  SpeedOfLight_D(x_)* E_GDB(i-1,j,k,y_,1)
             if(TypeFieldBC_S(2)=='laserbeam')&
                  call laser_beam(iDir=y_,       &
                  x=(i-0.50)*Dx_D(x_),&
                  y=j       *Dx_D(y_),&
                  z=(k-0.50)*Dx_D(z_),&
                  EField=E_GDB(i,j,k,y_,1) )
          end do
       end do
       do k=0-i1_D(z_),nZ+i1_D(z_)
          do j=1-iGCN*i0_D(y_),nY+iGCN*i0_D(y_)
             i=nX+iGCN
             E_GDB(i,j,k,z_,1) = E_GDB(i,j,k,z_,1)*(1 - SpeedOfLight_D(x_)) + &
                  SpeedOfLight_D(x_)*E_GDB(i-1,j,k,z_,1)
             if(TypeFieldBC_S(2)=='laserbeam')&
                  call laser_beam(iDir=z_,       &
                  x=(i-0.50)*Dx_D(x_),&
                  y=(j-0.50)*Dx_D(y_),&
                  z=k       *Dx_D(z_),&
                  EField=E_GDB(i,j,k,z_,1) )
          end do
       end do
    end if

    if(.not.IsPeriodicField_D(y_))then
       !\
       ! Face Y<0
       !/ 
       j=1-iGCN
       do k=1-iGCN*i0_D(z_),nZ+iGCN*i0_D(z_)
          do i=0-i1_D(x_),nX+i1_D(x_)
             E_GDB(i,j,k,x_,1) = E_GDB(i,j,k,x_,1)*(1 - SpeedOfLight_D(y_)) + &
                  SpeedOfLight_D(y_)*E_GDB(i,j+1,k,x_,1)
             if(TypeFieldBC_S(3)=='laserbeam')&
                  call laser_beam(iDir=x_,       &
                  x= i      *Dx_D(x_),&
                  y=(j-0.50)*Dx_D(y_),&
                  z=(k-0.50)*Dx_D(z_),&
                  EField=E_GDB(i,j,k,x_,1) )
          end do
       end do

       do k=0-i1_D(z_),nZ+i1_D(z_)
          do i=1-iGCN*i0_D(x_),nX+iGCN*i0_D(x_)
             E_GDB(i,j,k,z_,1) = E_GDB(i,j,k,z_,1)*(1 - SpeedOfLight_D(y_)) +&
                  SpeedOfLight_D(y_)*E_GDB(i,j+1,k,z_,1)
             if(TypeFieldBC_S(3)=='laserbeam')&
                  call laser_beam(iDir=z_,       &
                  x=(i-0.50)*Dx_D(x_),&
                  y=(j-0.50)*Dx_D(y_),&
                  z=k       *Dx_D(z_),&
                  EField=E_GDB(i,j,k,z_,1) )
          end do
       end do
       !\
       ! Face Y>nY.dy
       !/    
       j=nY+iGCN
       do k=1-iGCN*i0_D(z_),nZ+iGCN*i0_D(z_)
          do i=0-i1_D(x_),nX+i1_D(x_)
             E_GDB(i,j,k,x_,1) = E_GDB(i,j,k,x_,1)*(1 - SpeedOfLight_D(y_)) + &
                  SpeedOfLight_D(y_)*E_GDB(i,j-1,k,x_,1)
             if(TypeFieldBC_S(4)=='laserbeam')&
                  call laser_beam(iDir=x_,       &
                  x= i      *Dx_D(x_),&
                  y=(j-0.50)*Dx_D(y_),&
                  z=(k-0.50)*Dx_D(z_),&
                  EField=E_GDB(i,j,k,x_,1) )
          end do
       end do

       do k=0-i1_D(z_),nZ+i1_D(z_)
          do i=1-iGCN*i0_D(x_),nX+iGCN*i0_D(x_)
             E_GDB(i,j,k,z_,1) = E_GDB(i,j,k,z_,1)*(1 - SpeedOfLight_D(y_)) + &
                  SpeedOfLight_D(y_)*E_GDB(i,j-1,k,z_,1)
             if(TypeFieldBC_S(4)=='laserbeam')&
                  call laser_beam(iDir=z_,       &
                  x=(i-0.50)*Dx_D(x_),&
                  y=(j-0.50)*Dx_D(y_),&
                  z=k       *Dx_D(z_),&
                  EField=E_GDB(i,j,k,z_,1) )
          end do
       end do

    end if
    if(IsPeriodicField_D(z_)) RETURN
    !\
    ! Face Z<0
    !/
    k=1-iGCN
    do j=1-iGCN*i0_D(y_),nY+iGCN*i0_D(y_)
       do i=0-i1_D(x_),nX+i1_D(x_)
          E_GDB(i,j,k,x_,1) = E_GDB(i,j,k,x_,1)*(1 - SpeedOfLight_D(z_)) + &
               SpeedOfLight_D(z_)*E_GDB(i,j,k+1,x_,1)
          if(TypeFieldBC_S(5)=='laserbeam')&
               call laser_beam(iDir=x_,       &
               x= i      *Dx_D(x_),&
               y=(j-0.50)*Dx_D(y_),&
               z=(k-0.50)*Dx_D(z_),&
               EField=E_GDB(i,j,k,x_,1) )
       end do
    end do
    do j=0-i1_D(y_),nY+i1_D(y_)
       do i=1-iGCN*i0_D(x_),nX+iGCN*i0_D(x_)
          E_GDB(i,j,k,y_,1) = E_GDB(i,j,k,y_,1)*(1 - SpeedOfLight_D(z_)) + &
               SpeedOfLight_D(z_)*E_GDB(i,j,k+1,y_,1)
          if(TypeFieldBC_S(5)=='laserbeam')&
               call laser_beam(iDir=y_,       &
               x=(i-0.50)*Dx_D(x_),&
               y=j       *Dx_D(y_),&
               z=(k-0.50)*Dx_D(z_),&
               EField=E_GDB(i,j,k,y_,1) )
       end do
    end do

    !\
    ! Face Z>nZ.dz
    !/
    k=nZ+iGCN

    do j=1-iGCN*i0_D(y_),nY+iGCN*i0_D(y_)
       do i=0-i1_D(x_),nX+i1_D(x_)
          E_GDB(i,j,k,x_,1) = E_GDB(i,j,k,x_,1)*(1 - SpeedOfLight_D(z_)) + &
               SpeedOfLight_D(z_)*E_GDB(i,j,k-1,x_,1)
          if(TypeFieldBC_S(6)=='laserbeam')&
               call laser_beam(iDir=x_,       &
               x= i      *Dx_D(x_),&
               y=(j-0.50)*Dx_D(y_),&
               z=(k-0.50)*Dx_D(z_),&
               EField=E_GDB(i,j,k,x_,1) )
       end do
    end do
    do j=0-i1_D(y_),nY+i1_D(y_)
       do i=1-iGCN*i0_D(x_),nX+iGCN*i0_D(x_)
          E_GDB(i,j,k,y_,1) = E_GDB(i,j,k,y_,1)*(1 - SpeedOfLight_D(z_)) + &
               SpeedOfLight_D(z_)*E_GDB(i,j,k-1,y_,1)
          if(TypeFieldBC_S(6)=='laserbeam')&
               call laser_beam(iDir=y_,       &
               x=(i-0.50)*Dx_D(x_),&
               y=j       *Dx_D(y_),&
               z=(k-0.50)*Dx_D(z_),&
               EField=E_GDB(i,j,k,y_,1) )
       end do
    end do
  end subroutine field_bc
  !=======================
  subroutine get_rho_max(RhoMax,Coord_D)
    real,intent(out)::RhoMax
    real,dimension(nDim),optional,intent(out)::Coord_D
    real::Aux_D(4)
    !------------------
    RhoMax = maxval(Rho_GB(1:nX,1:nY,1:nZ,:))
    if(present(Coord_D))then
       Aux_D=real(maxloc(Rho_GB(1:nX,1:nY,1:nZ,:))) - cHalf
       Coord_D = Aux_D(1:nDim)
    end if
  end subroutine get_rho_max
  !---------------------------------------------------------------------!
  subroutine get_max_intensity(EnergyMax,Coord_D)
    real,intent(out)::EnergyMax
    real,dimension(nDim),optional,intent(out)::Coord_D
    integer::i,j,k,iBlock
    !-------------
    Rho_GB=cZero
    if(UseVectorPotential)then
       do iBlock = 1, MaxBlock
       call get_b_from_a(iBlock)
       do k=1,nZ; do j=1,nY; do i=1,nX
          Rho_GB(i,j,k,iBlock)=0.1250*(&
               sum(Counter_GDB(i,j-1:j,k-1:k,x_,iBlock)**2)+&
               sum(Counter_GDB(i-1:i,j,k-1:k,y_,iBlock)**2)+&
               sum(Counter_GDB(i-1:i,j-1:j,k,z_,iBlock)**2))+&
               0.250 * (&
               sum(E_GDB(i-1:i,j,k,x_,iBlock)**2)+&
               sum(E_GDB(i,j-1:j,k,y_,iBlock)**2)+&
               sum(E_GDB(i,j,k-1:k,z_,iBlock)**2))
       end do; end do; end do
       end do
    else
       do iBlock = 1, MaxBlock
       do k=1,nZ; do j=1,nY; do i=1,nX
          Rho_GB(i,j,k,iBlock)=0.1250*(&
               sum(B_GDB(i,j-1:j,k-1:k,x_,iBlock)**2)+&
               sum(B_GDB(i-1:i,j,k-1:k,y_,iBlock)**2)+&
               sum(B_GDB(i-1:i,j-1:j,k,z_,iBlock)**2))+&
               0.250*(&
               sum(E_GDB(i-1:i,j,k,x_,iBlock)**2)+&
               sum(E_GDB(i,j-1:j,k,y_,iBlock)**2)+&
               sum(E_GDB(i,j,k-1:k,z_,iBlock)**2))
       end do; end do; end do
       end do
    end if
    call get_rho_max(EnergyMax,Coord_D)
  end subroutine get_max_intensity
  !====================
  subroutine get_field_energy(Energy_V)
    use PIC_ModMain, ONLY: CellVolume
    real,intent(out)::Energy_V(6)
    integer::i,j,k,iBlock
    !-------------
    Energy_V = cZero
    if(UseVectorPotential)then
       do iBlock = 1, MaxBlock
       call get_b_from_a(iBlock)
       do k=1,nZ; do j=1,nY; do i=1,nX
          Energy_V(1) = Energy_V(1) + 0.1250*&
               sum((Counter_GDB(i,j-1:j,k-1:k,x_,iBlock) - B0_D(x_))**2)
          
          Energy_V(2) = Energy_V(2) + 0.1250*&
               sum((Counter_GDB(i-1:i,j,k-1:k,y_,iBlock) - B0_D(y_))**2)
          
          Energy_V(3) = Energy_V(3) + 0.1250*&
               sum((Counter_GDB(i-1:i,j-1:j,k,z_,iBlock) - B0_D(z_))**2)
          
          
          Energy_V(4) = Energy_V(4) + 0.250*&
               sum(E_GDB(i-1:i,j,k,x_,iBlock)**2)
          
          Energy_V(5) = Energy_V(5) + 0.250*&
               sum(E_GDB(i,j-1:j,k,y_,iBlock)**2)
          
          Energy_V(6) = Energy_V(6) + 0.250*&
               sum(E_GDB(i,j,k-1:k,z_,iBlock)**2)
       end do; end do; end do
       end do
    else
       do iBlock = 1, MaxBlock
       do k=1,nZ; do j=1,nY; do i=1,nX
          Energy_V(1) = Energy_V(1) + 0.1250*&
               sum((B_GDB(i,j-1:j,k-1:k,x_,iBlock) - B0_D(x_))**2)
          
          Energy_V(2) = Energy_V(2) + 0.1250*&
               sum((B_GDB(i-1:i,j,k-1:k,y_,iBlock) - B0_D(y_))**2)
          
          Energy_V(3) = Energy_V(3) + 0.1250*&
               sum((B_GDB(i-1:i,j-1:j,k,z_,iBlock) - B0_D(z_))**2)
          
          
          Energy_V(4) = Energy_V(4) + 0.250*&
               sum(E_GDB(i-1:i,j,k,x_,iBlock)**2)
          
          Energy_V(5) = Energy_V(5) + 0.250*&
               sum(E_GDB(i,j-1:j,k,y_,iBlock)**2)
          
          Energy_V(6) = Energy_V(6) + 0.250*&
               sum(E_GDB(i,j,k-1:k,z_,iBlock)**2)
       end do; end do; end do
       end do
    end if
    Energy_V = (CellVolume/(4.0 * cPi)) * Energy_V
  end subroutine get_field_energy
  !===============================
  !ONLY FOR PERIODIC BC.
  !================
  subroutine density_bc_periodic(iBlock)
    integer, intent(in)::iBlock
    integer :: i, j, k, iInner, jInner, kInner
    !----------------------------------------
    do k=1-iGCN,nZ+iGCN
       kInner = k - nZ*floor( (k - 0.50) /nZ )
       do j=1-iGCN,nY+iGCN
          jInner = j - nY*floor( (j - 0.50) /nY )
          do i=1-iGCN,nX+iGCN
             iInner = i - nX*floor( (i - 0.50) /nX )
             if(i==iInner.and.j==jInner.and.&
                  k==kInner)CYCLE
             Rho_GB(iInner,jInner,kInner,iBlock) = &
                  Rho_GB(iInner,jInner,kInner,iBlock) + &
                  Rho_GB(i,j,k,iBlock)
             Rho_GB(i,j,k,iBlock) = 0.0
          end do
       end do
    end do
  end subroutine density_bc_periodic
  !=================================
  subroutine current_bc_periodic(iBlock)
    integer,intent(in)::iBlock
    !\
    ! Local variables
    !/
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
                  Counter_GDB(i,j,k,x_,iBlock)
          end do
          !  0
          Current_F( 0,jInner,kInner) = &
               Current_F( 0,jInner,kInner) + &
               Counter_GDB( 0,j,k, x_,iBlock)      + &
               Counter_GDB(nX,j,k, x_,iBlock)
          ! 1:nX-1
          do i=1,nX-1
             iInner = i
             Current_F(iInner,jInner,kInner) = &
                  Current_F(iInner,jInner,kInner) + &
                  Counter_GDB(i,j,k, x_,iBlock)
          end do
          !nX
          Current_F(nX,jInner,kInner) = &
               Current_F(nX,jInner,kInner) + &
               Counter_GDB( 0,j,k, x_,iBlock)      + &
               Counter_GDB(nX,j,k, x_,iBlock)
          !nX+1
          do i=nX+1, nX+iGCN-1
             iInner = i - nX
             Current_F(iInner,jInner,kInner) = &
                  Current_F(iInner,jInner,kInner) + &
                  Counter_GDB(i,j,k, x_,iBlock)
          end do
       end do
    end do
    Counter_GDB(:,:,:,x_,iBlock) = 0.0
    Counter_GDB(0:nX,1:nY,1:nZ,x_,iBlock) = Current_F(0:nX,1:nY,1:nZ)
    
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
                  Counter_GDB(i,j,k, y_,iBlock)
          end do
       end do
       !0
       do i=1-iGCN,nX+iGCN
          iInner = i - nX*floor( (i - 0.50) /nX )
          Current_F( iInner,0,kInner) = &
               Current_F( iInner,0,kInner) + &
               Counter_GDB( i, 0,k, y_,iBlock)      + &
               Counter_GDB( i,nY,k, y_,iBlock)
       end do
       !1:nY-1
       do j=1,nY-1
          jInner = j 
          do i=1-iGCN,nX+iGCN
             iInner = i - nX*floor( (i - 0.50) /nX )
             Current_F(iInner,jInner,kInner) = &
                  Current_F(iInner,jInner,kInner) + &
                  Counter_GDB(i,j,k, y_,iBlock)
          end do
       end do
       !nY
       do i=1-iGCN,nX+iGCN
          iInner = i - nX*floor( (i - 0.50) /nX )
          Current_F( iInner,nY,kInner) = &
               Current_F( iInner,nY,kInner) + &
               Counter_GDB( i, 0,k, y_,iBlock)      + &
               Counter_GDB( i,nY,k, y_,iBlock)
       end do
       !-1
       do j=nY+1,nY+iGCN-1
          jInner = j - nY
          do i=1-iGCN,nX+iGCN
             iInner = i - nX*floor( (i - 0.50) /nX )
             Current_F(iInner,jInner,kInner) = &
                  Current_F(iInner,jInner,kInner) + &
                  Counter_GDB(i,j,k, y_,iBlock)
          end do
       end do
    end do
    Counter_GDB(:,:,:,y_,iBlock) = 0.0
    Counter_GDB(1:nX,0:nY,1:nZ,y_,iBlock) = Current_F(1:nX,0:nY,1:nZ)
    
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
                  Counter_GDB(i,j,k, z_,iBlock)
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
               Counter_GDB( i, j, 0, z_,iBlock)      + &
               Counter_GDB( i, j,nZ, z_,iBlock)
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
                  Counter_GDB(i,j,k, z_,iBlock)
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
               Counter_GDB( i, j, 0, z_,iBlock)      + &
               Counter_GDB( i, j,nZ, z_,iBlock)
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
                  Counter_GDB(i,j,k, z_,iBlock)
          end do
       end do
    end do
    
    Counter_GDB(:,:,:,z_,iBlock) = 0.0
    Counter_GDB(1:nX,1:nY,0:nZ,z_,iBlock) = Current_F(1:nX,1:nY,0:nZ)
    
  end subroutine current_bc_periodic

  !======================
  subroutine field_bc_periodic(iBlock)
    use PIC_ModMain, ONLY: IsPeriodicField_D
    integer,intent(in)::iBlock
    !\
    ! Local variables
    !/
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
             E_GDB(i,j,k,x_,iBlock) = &
                  E_GDB(iInner,jInner,kInner,x_,iBlock)
          end do
          !nX+1
          do i=nX+1, nX+iGCN-1
             iInner = i -  nXyzIfPeriodic_D(x_)
             E_GDB(i,j,k,x_,iBlock) = &
                  E_GDB(iInner,jInner,kInner,x_,iBlock)
          end do
          if(k==kInner.and.j==jInner)CYCLE
          !0:nX
          E_GDB(0:nX,j,k,x_,iBlock) = &
               E_GDB(0:nX,jInner,kInner,x_,iBlock)            
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
             E_GDB(i,j,k,y_,iBlock) = &
                  E_GDB(iInner,jInner,kInner,y_,iBlock) 
          end do
       end do
       !nY+1
       do j=nY+1,nY+iGCN-1
          jInner = j -  nXyzIfPeriodic_D(y_)
          do i=1-iGCN,nX+iGCN
             iInner = i -  nXyzIfPeriodic_D(x_)*floor( (i - 0.50) /nX )
             E_GDB(i,j,k,y_,iBlock) = &
                  E_GDB(iInner,jInner,kInner,y_,iBlock)
          end do
       end do
       !0:nY
       do j=0,nY
          jInner = j 
          do i=1-iGCN,nX+iGCN
             iInner = i -  nXyzIfPeriodic_D(x_)*floor( (i - 0.50) /nX )
             if(iInner==i.and.k==kInner)CYCLE
             E_GDB(i,j,k,y_,iBlock) = &
                  E_GDB(iInner,jInner,kInner,y_,iBlock)
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
             E_GDB(i,j,k,z_,iBlock) = &
                  E_GDB(iInner,jInner,kInner,z_,iBlock)
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
             E_GDB(i,j,k,z_,iBlock) = &
                  E_GDB(iInner,jInner,kInner,z_,iBlock)
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
             E_GDB(i,j,k,z_,iBlock) = &
                  E_GDB(iInner,jInner,kInner,z_,iBlock)
          end do
       end do
    end do
  end subroutine field_bc_periodic
end module PIC_ModField
