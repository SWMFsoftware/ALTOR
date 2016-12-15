!^CFG COPYRIGHT UofM
!--------------------------------------------------------------!

module PIC_ModField
  use PIC_ModSize, ONLY: nX, nY, x_, y_, z_, nDim
  use ModNumConst
  use PIC_ModMain,ONLY:&
       SpeedOfLight_D,&     !This is c\Delta t/\Delta x,...
       UseVectorPotential   !If we use vector-potential, or not

  use PIC_ModFormFactor, ONLY: lOrderFF
  implicit none



  !Introduce the number of ghostcells (iGCN)
  integer,parameter:: iGCN=&
       1 &    !because even with the first order formfactor the 
              !particle can escape the domain within a time step
  +lOrderFF/2 
              !because for lOrderFF>1 the particle form-factor
              !is wider than a mesh size
  !Structures
  !      
  !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:)
  !       _!_!x!x               _!_!_!_            Ey (:,nY+2)
  !       _!_!_!                x!x!_!_  Bz(:,nY+2),Bz(nX+2,:)       
  !    -2  ! ! !                x!x!_!_            Bx (:,nY+2)
  !       -2                                       By (nX+2,:)
  real,dimension(&
       1-iGCN:nX+iGCN,&
       1-iGCN:nY+iGCN,3)::&
       E_GD        = 0.0,& !This is the electric field
       Magnetic_GD = 0.0,& !vector-potential if used, magnetic field otherwise
       Counter_GD  = 0.0   !Counter for electric current
  real,dimension(&
       1-iGCN:nX+iGCN,&
       1-iGCN:nY+iGCN) :: rho_G = 0.0 

  !Methods

  public::get_b_from_a   !Transforms the vector potential to magn. field 
  public::update_magnetic!Updates the magnetic field, or vector potential
  public::update_e       !Updates the electric field
  public::get_rho_max    !
contains
  !=======================
  subroutine add_e
    use ModReadParam, ONLY: read_var
    real:: E_D(3)
    integer:: i,j,iDim
    !------------
    call read_var('Ex',E_D(1))
    call read_var('Ey',E_D(2))
    call read_var('Ez',E_D(3))
    do iDim = 1,3
       do j=1-iGCN,nY+iGCN
          do i=1-iGCN,nX+iGCN
             E_GD(i,j,iDim) = &
                  E_GD(i,j,iDim) + E_D(iDim)
          end do
       end do
    end do
  end subroutine add_e
  !===================
  subroutine add_b
    use ModReadParam, ONLY: read_var
    real:: B_D(3)
    integer:: i,j,iDim
    !------------
    call read_var('Bx',B_D(1))
    call read_var('By',B_D(2))
    call read_var('Bz',B_D(3))
    do iDim = 1,3
       do j=1-iGCN,nY+iGCN
          do i=1-iGCN,nX+iGCN
             Magnetic_GD(i,j,iDim) = &
                  Magnetic_GD(i,j,iDim) + B_D(iDim)
          end do
       end do
    end do
  end subroutine add_b
  !===================
  !=====================
  !Transforms the vector potential to magn. field
  subroutine get_b_from_a
    integer::i,j
    !Applied only if UseVectorPotential
    !The vector potential is in Magnetic_GD
    !Calculates the magnetic field 
    
    !The whole magnetic field is calculated and 
    !the result is saved to Counter_GD

    do j=1-iGCN,nY+iGCN-1; do i=1-iGCN,nX+iGCN
       !      
       !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:)
       !       _!_!x!x               _!_!_!_            Ey (:,nY+2)
       !       _!_!_!                x!x!_!_  Bz(:,nY+2),Bz(nX+2,:)       
       !    -2  ! ! !                x!x!_!_            Bx (:,nY+2)
       !       -2                                       By (nX+2,:)
       Counter_GD(i,j,x_) =   &
            SpeedOfLight_D(y_)*&
              (Magnetic_GD(i,  j+1,z_) - Magnetic_GD(i,j,z_))
    end do; end do 

    do j=1-iGCN,nY+iGCN; do i=1-iGCN,nX+iGCN-1
       !      
       !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:)
       !       _!_!x!x               _!_!_!_            Ey (:,nY+2)
       !       _!_!_!                x!x!_!_  Bz(:,nY+2),Bz(nX+2,:)       
       !    -2  ! ! !                x!x!_!_            Bx (:,nY+2)
       !       -2                                       By (nX+2,:)
       Counter_GD(i,j,y_) = - &
            SpeedOfLight_D(x_)*&
              (Magnetic_GD(i+1,j  ,z_) - Magnetic_GD(i,j,z_))
    end do; end do

    do j=1-iGCN,nY+iGCN-1; do i=1-iGCN,nX+iGCN-1  
       !      
       !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:)
       !       _!_!x!x               _!_!_!_            Ey (:,nY+2)
       !       _!_!_!                x!x!_!_  Bz(:,nY+2),Bz(nX+2,:)       
       !    -2  ! ! !                x!x!_!_            Bx (:,nY+2)
       !       -2                                       By (nX+2,:)   
       Counter_GD(z_,i,j) =   &
            SpeedOfLight_D(x_)*&
              (Magnetic_GD(i+1,j  ,y_) - Magnetic_GD(i,j,y_))-&
            SpeedOfLight_D(y_)*& 
              (Magnetic_GD(i,  j+1,x_) - Magnetic_GD(i,j,x_))
    end do; end do
  
  end subroutine get_b_from_a
  !==========================
  subroutine update_magnetic
    integer::i,j
    real,dimension(nDim) :: SpeedOfLightHalf_D
    if(UseVectorPotential)then
       !      
       !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:)
       !       _!_!x!x               _!_!_!_            Ey (:,nY+2)
       !       _!_!_!                x!x!_!_  Bz(:,nY+2),Bz(nX+2,:)       
       !    -2  ! ! !                x!x!_!_            Bx (:,nY+2)
       !       -2                                       By (nX+2,:)
       !Advance vector potential throught a half timestep

       Magnetic_GD(1-iGCN:nX+iGCN-1,:,x_) = &
            Magnetic_GD(1-iGCN:nX+iGCN-1,:,x_) - &
             cHalf*E_GD(1-iGCN:nX+iGCN-1,:,x_)

       Magnetic_GD(:,1-iGCN:nY+iGCN-1,y_) = &
            Magnetic_GD(:,1-iGCN:nY+iGCN-1,y_) - &
             cHalf*E_GD(:,1-iGCN:nY+iGCN-1,y_)

       Magnetic_GD(:,:,z_) = &
            Magnetic_GD(:,:,z_) - &
             cHalf*E_GD(:,:,z_)
       return
    end if
  
    SpeedOfLightHalf_D=&
         SpeedOfLight_D*cHalf
    do j=1-iGCN,nY+iGCN-1; do i=1-iGCN,nX+iGCN
       !      
       !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:)
       !       _!_!x!x               _!_!_!_            Ey (:,nY+2)
       !       _!_!_!                x!x!_!_  Bz(:,nY+2),Bz(nX+2,:)       
       !    -2  ! ! !                x!x!_!_            Bx (:,nY+2)
       !       -2                                       By (nX+2,:)
       !
       Magnetic_GD(i,j,x_) = Magnetic_GD(i,j,x_) - &
            SpeedOfLightHalf_D(y_)*&
                         (E_GD(i,  j+1,z_) - E_GD(i,j,z_))
    end do; end do

    do j=1-iGCN,nY+iGCN; do i=1-iGCN,nX+iGCN-1
       !      
       !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:)
       !       _!_!x!x               _!_!_!_            Ey (:,nY+2)
       !       _!_!_!                x!x!_!_  Bz(:,nY+2),Bz(nX+2,:)       
       !    -2  ! ! !                x!x!_!_            Bx (:,nY+2)
       !       -2                                       By (nX+2,:)
       !
       Magnetic_GD(i,j,y_) = Magnetic_GD(i,j,y_) + &
            SpeedOfLightHalf_D(x_)*&
                         (E_GD(i+1,j  ,z_) - E_GD(i,j,z_))

    end do; end do

    do j=1-iGCN,nY+iGCN-1; do i=1-iGCN,nX+iGCN-1
       !      
       !       _!_!x!x               _!_!_!_  Not used: Ex (nX+2,:)
       !       _!_!x!x               _!_!_!_            Ey (:,nY+2)
       !       _!_!_!                x!x!_!_  Bz(:,nY+2),Bz(nX+2,:)       
       !    -2  ! ! !                x!x!_!_            Bx (:,nY+2)
       !       -2                                       By (nX+2,:)
       Magnetic_GD(i,j,z_) = Magnetic_GD(i,j,z_) - &
            SpeedOfLightHalf_D(x_)*&
                         (E_GD(i+1,j  ,y_) - E_GD(i,j,y_))+&
            SpeedOfLightHalf_D(y_)*&
                         (E_GD(i,  j+1,x_) - E_GD(i,j,x_))
    end do; end do
    
  end subroutine update_magnetic
!-----------------------------------------------------------------!
  subroutine update_e
    !Add current
    E_GD(0:nX, 1:nY, x_) = E_GD(0:nX, 1:nY, x_) - &
                     Counter_GD(0:nX, 1:nY ,x_)
    E_GD(1:nX, 0:nY, y_) = E_GD(1:nX, 0:nY, y_) - &
                     Counter_GD(1:nX, 0:nY, y_)
    E_GD(1:nX, 1:nY, z_) = E_GD(1:nX, 1:nY, z_) - &
                     Counter_GD(1:nX, 1:nY, z_)
    if(UseVectorPotential)then
       call get_b_from_a()!Put the magnetic field to Counter_DB
       call use_magnetic_field_from(Counter_GD)
    else
       call use_magnetic_field_from(Magnetic_GD)
    end if
  contains
    subroutine use_magnetic_field_from(BIn_GD)
      real,dimension(&
           1-iGCN:nX+iGCN,&
           1-iGCN:nY+iGCN,3),intent(in)::&
           BIn_GD
      integer::i,j
      !-----------
      do j=1,nY; do i=0,nX
         E_GD(i,j,x_) = E_GD(i,j,x_) + &
              SpeedOfLight_D(y_)*&
                     (BIn_GD(i,j,z_) - BIn_GD(i,  j-1, z_))
      end do; end do

      do j=0,nY; do i=1,nX
         E_GD(i,j,y_) = E_GD(i,j,y_) - &
              SpeedOfLight_D(x_)*&
                     (BIn_GD(i,j,z_) - BIn_GD(i-1,j  , z_))
      end do; end do

      do j=1,nY; do i=1,nX
         E_GD(i,j,z_) = E_GD(z_,i,j) + &
              SpeedOfLight_D(x_)*&
                     (BIn_GD(i,j,y_) - BIn_GD(i-1,j  , y_)) - &
              SpeedOfLight_D(y_)*&
                     (BIn_GD(i,j,x_) - BIn_GD(i,  j-1, x_))
      end do; end do
    end subroutine use_magnetic_field_from
  end subroutine update_e
  !======================
  !====================
  subroutine get_field_energy(Energy_V)
    use PIC_ModMain, ONLY: CellVolume
    real,intent(out)::Energy_V(6)
    integer::i,j
    !-------------
    Energy_V = cZero
    if(UseVectorPotential)then
       call get_b_from_a()
       do j=1,nY; do i=1,nX
          Energy_V(1) = Energy_V(1) + 0.250*&
               sum(Counter_GD(i,j-1:j,x_)**2)

          Energy_V(2) = Energy_V(2) + 0.250*&
               sum(Counter_GD(i-1:i,j,y_)**2)

          Energy_V(3) = Energy_V(3) + 0.1250*&
               sum(Counter_GD(i-1:i,j-1:j,z_)**2)


          Energy_V(4) = Energy_V(4) + 0.250*&
               sum(E_GD(i-1:i,j,x_)**2)

          Energy_V(5) = Energy_V(5) + 0.250*&
               sum(E_GD(i,j-1:j,y_)**2)

          Energy_V(6) = Energy_V(6) + 0.50*E_GD(i,j,z_)**2
       end do; end do
    else
       do j=1,nY; do i=1,nX
          Energy_V(1) = Energy_V(1) + 0.250*&
               sum(Magnetic_GD(i,j-1:j,x_)**2)

          Energy_V(2) = Energy_V(2) + 0.250*&
               sum(Magnetic_GD(i-1:i,j,y_)**2)

          Energy_V(3) = Energy_V(3) + 0.1250*&
               sum(Magnetic_GD(i-1:i,j-1:j,z_)**2)


          Energy_V(4) = Energy_V(4) + 0.250*&
               sum(E_GD(i-1:i,j,x_)**2)

          Energy_V(5) = Energy_V(5) + 0.250*&
               sum(E_GD(i,j-1:j,y_)**2)

          Energy_V(6) = Energy_V(6) + 0.50*E_GD(i,j,z_)**2
       end do; end do
    end if
    Energy_V = (CellVolume/(4.0 * cPi)) * Energy_V
  end subroutine get_field_energy
  !===============================
  subroutine get_rho_max(RhoMax,Coord_D)
    real,intent(out) :: RhoMax
    real,dimension(nDim),optional,intent(out) :: Coord_D
    !-----------------
    RhoMax = maxval(rho_G(1:nX,1:nY))
    if(present(Coord_D))Coord_D=real(maxloc(rho_G(1:nX,1:nY))) - cHalf
  end subroutine get_rho_max
  !---------------------------------------------------------------------!
  subroutine get_intensity(EnergyMax,Coord_D)
    real,intent(out)::EnergyMax
    real,dimension(nDim),optional,intent(out)::Coord_D
    integer::i,j
    !------------
    rho_G=cZero
    if(UseVectorPotential)then

       call get_b_from_a()

       do j=1,nY; do i=1,nX

          rho_G(i,j)=cHalf*(cHalf*(&
               sum(E_GD(i-1:i,j,x_)**2)+&
               sum(E_GD(i,j-1:j,y_)**2)+&
               sum(Counter_GD(i,j-1:j,x_)**2)+&
               sum(Counter_GD(i-1:i,j,y_)**2))+&
               0.250*&
               sum(Counter_GD(i-1:i,j-1:j,z_)**2)+&
               E_GD(i,j,z_)**2)

       end do; end do
    else
       do j=1,nY; do i=1,nX
          rho_G(i,j)=cHalf*(cHalf*(&
               sum(E_GD(i-1:i,j,x_)**2)+&
               sum(E_GD(i,j-1:j,y_)**2)+&
               sum(Magnetic_GD(i,j-1:j,x_)**2)+&
               sum(Magnetic_GD(i-1:i,j,y_)**2))+&
               0.250*&
               sum(Magnetic_GD(i-1:i,j-1:j,z_)**2)+&
               E_GD(i,j,z_)**2)
       end do; end do
    end if
  end subroutine get_intensity
end module PIC_ModField
