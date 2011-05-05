!^CFG COPYRIGHT UofM
!--------------------------------------------------------------!

module PIC_ModField
  use PIC_ModGrid
  use ModNumConst
  use PIC_ModMain,ONLY:&
       SpeedOfLight_D,&     !This is c\Delta t/\Delta x,...
       UseVectorPotential,& !If we use vector-potential, or not
       DoAccelerateLight    !If we force the propagation speed
                            !of light in vacuum to be not
                            !less than c
  use PIC_ModFormFactor
  implicit none
  !Introduce the number of ghostcells (iGCN)
  integer,parameter:: iGCN=&
       1 &    !because even with the first order formfactor the 
              !particle can escape the domain within a time step
  +lOrderFF/2 
              !because for lOrderFF>1 the particle form-factor
              !is wider than a mesh size
  !Structures
  real,dimension(3,&
       0-iGCN:nX+iGCN,&
       0-iGCN:nY+iGCN,&
       0-iGCN:nZ+iGCN)::&
       E_DG,& !This is the electric field
       Magnetic_DG,& !vector-potential if used, magnetic field otherwise
       Counter_DG    !Counter for electric current
  real,dimension(&
       1-iGCN:nX+iGCN,&
       1-iGCN:nY+iGCN,&
       1-iGCN:nZ+iGCN)::rho_G
 !Methods
  public::evaluate_memory_for_field !Summarizes the memory requirements
  public::get_b_from_a   !Transforms the vector potential to magn. field
  public::update_magnetic!Updates the magnetic field, or vector potential
  public::update_e       !Updates the electric field
  public::get_rho_max     !
contains
  !----------------------------------------------------------------------!
  subroutine init_field
    use ModNumConst
    E_DG=cZero
    Magnetic_DG=cZero
    Counter_DG=cZero
    rho_G=cZero
  end subroutine init_field
  !----------------------------------------------------------------------!
  subroutine evaluate_memory_for_field
    use ModMpi,ONLY:iRealPrec
    use ModNumConst
    integer,parameter::K=2**10,M=K*K
    
    write(*,'(a,f6.1)')&
         'To store the fields, the memory requirement is: ',&
         real((nX+1+2*iGCN)*&
         (nY+1+2*iGCN)*&
         (nZ+1+2*iGCN)*9+&
         nX*nY*nZ)&          !For density
         *real(iRealPrec+1)& !*2, if compiled with a double prec-n
         *cFour/real(M),'+/-0.1 MB'
  end subroutine evaluate_memory_for_field
!----------------------------------------------------------------!
  subroutine get_b_from_a
    integer::i,j,k
    !Applied only if UseVectorPotential
    !The vector potential is in Magnetic_DG
    !Calculates the magnetic field 
    
    !The whole magnetic field is calculated and 
    !the result is saved to Counter_DG
    do k=0-iGCN,nZ+iGCN-1
       do j=0-iGCN,nY+iGCN-1
          do i=0-iGCN,nX+iGCN-1
             Counter_DG(x_,i,j,k)=&
                  SpeedOfLight_D(y_)*&
                  (Magnetic_DG(z_,i,  j+1,k  )-Magnetic_DG(z_,i,j,k))-&
                  SpeedOfLight_D(z_)*&
                  (Magnetic_DG(y_,i,  j,  k+1)-Magnetic_DG(y_,i,j,k))
             
             Counter_DG(y_,i,j,k)=&
                  SpeedOfLight_D(z_)*&
                  (Magnetic_DG(x_,i,  j,  k+1)-Magnetic_DG(x_,i,j,k))-&
                  SpeedOfLight_D(x_)*&
                  (Magnetic_DG(z_,i+1,j,  k  )-Magnetic_DG(z_,i,j,k))
             
             Counter_DG(z_,i,j,k)=&
                  SpeedOfLight_D(x_)*&
                  (Magnetic_DG(y_,i+1,j,  k  )-Magnetic_DG(y_,i,j,k))-&
                  SpeedOfLight_D(y_)*& 
                  (Magnetic_DG(x_,i,  j+1,k  )-Magnetic_DG(x_,i,j,k))
          end do
       end do
    end do
  end subroutine get_b_from_a
!----------------------------------------------------------------!

  subroutine update_magnetic
    integer::i,j,k
    real,dimension(nDim),parameter::SpeedOfLightHalf_D=&
         SpeedOfLight_D*cHalf
    if(UseVectorPotential)then
       !Advance vector potential throught a half timestep
       Magnetic_DG=Magnetic_DG-cHalf*E_DG
       return
    end if
    do k=0,nZ
       do j=0,nY
          do i=0,nX
             Magnetic_DG(x_,i,j,k)=Magnetic_DG(x_,i,j,k)-&
              SpeedOfLightHalf_D(y_)*&
                           (E_DG(z_,i,  j+1,k  )-E_DG(z_,i,j,k))+&
              SpeedOfLightHalf_D(z_)*&
                           (E_DG(y_,i,  j,  k+1)-E_DG(y_,i,j,k))

             Magnetic_DG(y_,i,j,k)=Magnetic_DG(y_,i,j,k)-&
              SpeedOfLightHalf_D(z_)*&
                           (E_DG(x_,i,  j,  k+1)-E_DG(x_,i,j,k))+&
              SpeedOfLightHalf_D(x_)*&
                           (E_DG(z_,i+1,j,  k  )-E_DG(z_,i,j,k))

             Magnetic_DG(z_,i,j,k)=Magnetic_DG(z_,i,j,k)-&
              SpeedOfLightHalf_D(x_)*&
                           (E_DG(y_,i+1,j,  k  )-E_DG(y_,i,j,k))+&
              SpeedOfLightHalf_D(y_)*&
                           (E_DG(x_,i,  j+1,k  )-E_DG(x_,i,j,k))
          end do
       end do
    end do
  end subroutine update_magnetic
!-----------------------------------------------------------------!
  subroutine update_e
    !Add current
    E_DG=E_DG-Counter_DG
    if(UseVectorPotential)then
       call get_b_from_a()!Put the magnetic field to Counter_DB
       call use_magnetic_field_from(Counter_DG)
    else
       call use_magnetic_field_from(Magnetic_DG)
    end if
  contains
    subroutine use_magnetic_field_from(BIn_DG)
      real,dimension(3,&
           0-iGCN:nX+iGCN,&
           0-iGCN:nY+iGCN,&
           0-iGCN:nZ+iGCN),intent(in)::&
           BIn_DG
      integer::i,j,k
      do k=0,nZ
         do j=0,nY
            do i=0,nX
               E_DG(x_,i,j,k)=E_DG(x_,i,j,k)+&
                    SpeedOfLight_D(y_)*&
                           (BIn_DG(z_,i,j,k)-BIn_DG(z_,i,  j-1,k  ))-&
                    SpeedOfLight_D(z_)*&
                           (BIn_DG(y_,i,j,k)-BIn_DG(y_,i,  j,  k-1))

               E_DG(y_,i,j,k)=E_DG(y_,i,j,k)+&
                    SpeedOfLight_D(z_)*&
                           (BIn_DG(x_,i,j,k)-BIn_DG(x_,i,  j,  k-1))-&
                    SpeedOfLight_D(x_)*&
                           (BIn_DG(z_,i,j,k)-BIn_DG(z_,i-1,j,  k  ))

               E_DG(z_,i,j,k)=E_DG(z_,i,j,k)+&
                    SpeedOfLight_D(x_)*&
                           (BIn_DG(y_,i,j,k)-BIn_DG(y_,i-1,j,  k  ))-&
                    SpeedOfLight_D(y_)*&
                           (BIn_DG(x_,i,j,k)-BIn_DG(x_,i,  j-1,k  ))
            end do
         end do
      end do
    end subroutine use_magnetic_field_from
  end subroutine update_e
  !---------------------------------------------------------------------!
  subroutine get_rho_max(RhoMax,Coord_D)
    real,intent(out)::RhoMax
    real,dimension(nDim),optional,intent(out)::Coord_D
    RhoMax = maxval(rho_G(1:nX,1:nY,1:nZ))
    if(present(Coord_D))Coord_D=real(maxloc(rho_G(1:nX,1:nY,1:nZ)))-cHalf
  end subroutine get_rho_max
  !---------------------------------------------------------------------!
  subroutine get_max_intensity(EnergyMax,Coord_D)
    real,intent(out)::EnergyMax
    real,dimension(nDim),optional,intent(out)::Coord_D
    integer::i,j,k
    rho_G=cZero
    if(UseVectorPotential)then
       call get_b_from_a()
       do k=1,nZ
          do j=1,nY
             do i=1,nX
                rho_G(i,j,k)=cEighth*(&
                     sum(Counter_DG(x_,i,j-1:j,k-1:k)**2)+&
                     sum(Counter_DG(y_,i-1:i,j,k-1:k)**2)+&
                     sum(Counter_DG(z_,i-1:i,j-1:j,k)**2))+&
                     cQuarter*(&
                     sum(E_DG(x_,i-1:i,j,k)**2)+&
                     sum(E_DG(y_,i,j-1:j,k)**2)+&
                     sum(E_DG(z_,i,j,k-1:k)**2))
             end do
          end do
       end do
    else
       do k=1,nZ
          do j=1,nY
             do i=1,nX
                rho_G(i,j,k)=cEighth*(&
                     sum(Magnetic_DG(x_,i,j-1:j,k-1:k)**2)+&
                     sum(Magnetic_DG(y_,i-1:i,j,k-1:k)**2)+&
                     sum(Magnetic_DG(z_,i-1:i,j-1:j,k)**2))+&
                     cQuarter*(&
                     sum(E_DG(x_,i-1:i,j,k)**2)+&
                     sum(E_DG(y_,i,j-1:j,k)**2)+&
                     sum(E_DG(z_,i,j,k-1:k)**2))
             end do
          end do
       end do
    end if
    call get_rho_max(EnergyMax,Coord_D)
  end subroutine get_max_intensity
  !---------------------------------------------------------------------!
end module PIC_ModField
!=======================================================================!
