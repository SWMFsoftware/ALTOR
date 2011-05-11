!^CFG COPYRIGHT UofM
!--------------------------------------------------------------!

module PIC_ModField
  use PIC_ModGrid
  use ModNumConst
  use PIC_ModMain,ONLY:&
       SpeedOfLight_D,&     !This is c\Delta t/\Delta x,...
       UseVectorPotential   !If we use vector-potential, or not

  use PIC_ModFormFactor
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

 !Methods
  ! public::evaluate_memory_for_field !Summarizes the memory requirements
  public::get_b_from_a   !Transforms the vector potential to magn. field
  public::update_magnetic!Updates the magnetic field, or vector potential
  public::update_e       !Updates the electric field
  public::get_rho_max     
contains
  !=================================
  !subroutine evaluate_memory_for_field
  !  use ModMpi,ONLY:iRealPrec
  !  use ModNumConst
  !  integer,parameter::K=2**10,M=K*K
  !  
  !  write(*,'(a,f6.1)')&
  !       'To store the fields, the memory requirement is: ',&
  !       real((nX+1+2*iGCN)*&
  !       (nY+1+2*iGCN)*&
  !       (nZ+1+2*iGCN)*9+&
  !       nX*nY*nZ)&          !For density
  !       *real(iRealPrec+1)& !*2, if compiled with a double prec-n
  !       *4.0/real(M),'+/-0.1 MB'
  !end subroutine evaluate_memory_for_field
  !=======================
  subroutine get_b_from_a
    integer::i,j,k
    !Applied only if UseVectorPotential
    !The vector potential is in Magnetic_GD
    !Calculates the magnetic field 
    
    !The whole magnetic field is calculated and 
    !the result is saved to Counter_GD

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
      !---------------
      do k=1,nZ; do j=1,nY; do i=0,nX
         E_GD(i,j,k,x_) = E_GD(i,j,k,x_) + &
              SpeedOfLight_D(y_)*&
                     (BIn_GD(i,j,k,z_) - BIn_GD(i,  j-1,k  ,z_)) - &
              SpeedOfLight_D(z_)*&
                     (BIn_GD(i,j,k,y_) - BIn_GD(i,  j,  k-1,y_))
      end do; end do; end do
      do k=1,nZ; do j=0,nY; do i=1,nX
         E_GD(i,j,k,y_) = E_GD(i,j,k,y_) + &
              SpeedOfLight_D(z_)*&
                     (BIn_GD(i,j,k,x_) - BIn_GD(i,  j,  k-1,x_))-&
              SpeedOfLight_D(x_)*&
                     (BIn_GD(i,j,k,z_) - BIn_GD(i-1,j,  k  ,z_))
      end do; end do; end do
      do k=0,nZ; do j=1,nY; do i=1,nX
         E_GD(i,j,k,z_) = E_GD(i,j,k,z_)+&
              SpeedOfLight_D(x_)*&
                     (BIn_GD(i,j,k,y_) - BIn_GD(i-1,j,  k  ,y_))-&
              SpeedOfLight_D(y_)*& 
                     (BIn_GD(i,j,k,x_) - BIn_GD(i,  j-1,k  ,x_))
      end do; end do; end do
    end subroutine use_magnetic_field_from
  end subroutine update_e
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
  !===============================
end module PIC_ModField
