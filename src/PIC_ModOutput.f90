module PIC_ModOutput
  use PIC_ModMain,ONLY: iStep,tSimulation,Dt,DX_D 
  use PIC_ModSize,ONLY: nDim, nCell_D, nX, nY, nZ, x_, y_, z_
  use PIC_ModSize,ONLY: MaxBlock, nPType

  use PIC_ModParticles,ONLY: M_P,Q_P
  use ModNumConst,Only: cPi

  use PIC_ModParticleInField,ONLY: Rho_GB, V_GDB
  use PIC_ModField, ONLY: iGCN, E_GDB,Magnetic_GDB
  use PIC_ModMpi,  ONLY: pass_density
  use PIC_ModProc,      ONLY: iProc
  use ModIoUnit,        ONLY: io_unit_new
  implicit none
  SAVE
  PRIVATE !Except  

  !\                                                                          
  !Unit for the log file                                                      
  !/                           

  integer :: iOutUnit = -1,iOutUnit_den=-1
  integer, public :: nStepOut=100000,nStepOutMin=100000
  !Array for saving coordinates, fields and moments in one timestep
  real   :: State_VCB(6+11*nPType,1:nX,1:nY,1:nZ,MaxBlock)=0.0

  public :: write_density,write_field

  public :: write_density_test, write_momentum

  public :: PIC_save_files

contains
  !======================================================================
  subroutine write_density(iSort)
    character(LEN=100) :: Name
    integer, intent(in) :: iSort
    real :: rmax
    call pass_density(0)
    rmax=maxval(Rho_GB)
    if(iProc==0.and.rmax>1.e-10)then
       write(Name,'(a,i1,a,i4.4,a)')'n',iSort,'_',iStep-1,'.dat'
       iOutUnit = io_unit_new()
       open(iOutUnit,file=trim(Name),status='replace',form='unformatted')
       write(iOutUnit) &
            iStep,nDim, iGCN, nCell_D,tSimulation-dt,dt,Dx_D &
            ,minval(Rho_GB),rmax

       write(iOutUnit) Rho_GB(:,:,nZ/2,1)
       write(iOutUnit) Rho_GB(:,nY/2,:,1)
       close(iOutUnit)
       write(*,*) ' Density iSort=',iSort,'  Max=',rmax
    end if
  end subroutine write_density
  !======================================================================
  subroutine write_field(iX)
    character(LEN=100) :: Name
    integer, intent(in) :: iX
    if(iProc==0)then
       write(Name,'(a,i1,a,i4.4,a)')'e',iX,'_',iStep,'.dat'
       iOutUnit = io_unit_new()
       open(iOutUnit,file=trim(Name),status='replace',form='unformatted')

       write(iOutUnit) &
            iStep,nDim, iGCN, nCell_D,tSimulation,dt,Dx_D &
            ,minval(E_GDB(1:nX,1:nY,1:nZ,iX,:)),maxval(E_GDB(:,:,:,iX,:))
       write(iOutUnit) E_GDB(:,:,nZ/2,iX,:)
       write(iOutUnit) E_GDB(:,nY/2,:,iX,:)                                 
       close(iOutUnit)
    end if
  end subroutine write_field

  !=====================================================================
  !Save the number density at each timestep 
  subroutine write_density_test(iSort)
    integer, intent(in) :: iSort
    !------------------------------------------------------------------
    call pass_density(0)            

    if(iProc==0)then
       State_VCB(iSort+6,:,:,:,1) = Rho_GB(1:nX,1:nY,1:nZ,1)
    end if

  end subroutine write_density_test
  !====================================
  !Save the velocity and pressure at each timestep
  subroutine write_momentum(iSort)
    integer, intent(in) :: iSort
    integer :: i,j,k
    !------------------------------------------------------------------
    if(iProc==0)then
       !Save velocity
       State_VCB(3*iSort+6,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,1,1)
       State_VCB(3*iSort+7,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,2,1)
       State_VCB(3*iSort+8,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,3,1)
       !Save pressure tensor: Pxx Pyy Pzz Pxy Pxz Pyz
       State_VCB(6*iSort+11,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,1,1)**2
       State_VCB(6*iSort+12,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,2,1)**2
       State_VCB(6*iSort+13,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,3,1)**2
       State_VCB(6*iSort+14,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,1,1)*&
            V_GDB(1:nX,1:nY,1:nZ,2,1)
       State_VCB(6*iSort+15,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,1,1)*&
            V_GDB(1:nX,1:nY,1:nZ,3,1)
       State_VCB(6*iSort+16,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,2,1)*&
            V_GDB(1:nX,1:nY,1:nZ,3,1)
       !Save scalar pressure
       do k=1,nX; do j=1,nY; do i=1,nX
          State_VCB(iSort+14,i,j,k,1) = sum(&
               State_VCB(6*iSort+11:6*iSort+13,i,j,k,1))/3.0
       end do; end do; end do
    end if

  end subroutine write_momentum


  !===========================================================================

  subroutine PIC_save_files(TypeSave)

    use ModPlotFile, ONLY: save_plot_file
    use PIC_ModParticles, ONLY: nPType,n_P
    use PIC_ModFormFactor,ONLY: HighFF_ID,Node_D, get_form_factors
    use PC_BATL_particles,ONLY: Particle_I
    use PIC_ModParticleInField,ONLY: add_density, add_velocity
    use PIC_ModMain, ONLY: c,c2   
    use ModMpi, ONLY: mpi_reduce_real_array, MPI_SUM
    use PIC_ModProc, ONLY: iProc, iComm, nProc

    character(len=*), intent(in):: TypeSave
    character(len=6) :: TypePosition
    real    :: Param_I(1:2*nPType)
    integer :: iX,i,j,k, iSort, iParticle
    character(len=15) :: NameFile='variable.outs'

    real,dimension(nDim):: X_D
    real,dimension(nDim):: V_D
    integer :: iBlock=1
    integer :: iError

    character(len=:), allocatable, save:: NameVar
    character(len=*), parameter:: NameSub = 'save_files'
    !----------------------------------------------------------------------

    !
    call pass_density(0)

    select case(TypeSave)
    case('INITIAL')
       !Calculate the initial number densities and velocities
       do iSort = 1, nPType
          Rho_GB = 0.0 ; V_GDB = 0.0
          do iParticle = 1, n_P(iSort)
             !Get the position of particle
             X_D = Particle_I(iSort)%State_VI(1:nDim,iParticle)
             !Get the velocity of particle from generalized momentum
             V_D = Particle_I(iSort)%State_VI(nDim+1:nDim+3,iParticle)
             V_D = c/(sqrt(c2+sum(V_D**2)))*V_D
             call get_form_factors(X_D,Node_D,HighFF_ID)
             call add_density(Node_D, HighFF_ID, iBlock)
             call add_velocity(V_D,Node_D, HighFF_ID,iBlock)
          end do
          !Save number densities
          State_VCB(iSort+6,:,:,:,1) = Rho_GB(1:nX,1:nY,1:nZ,1)
          !Save velocities
          State_VCB(3*iSort+6,:,:,:,iBlock) = V_GDB(1:nX,1:nY,1:nZ,1,iBlock)
          State_VCB(3*iSort+7,:,:,:,iBlock) = V_GDB(1:nX,1:nY,1:nZ,2,iBlock)
          State_VCB(3*iSort+8,:,:,:,iBlock) = V_GDB(1:nX,1:nY,1:nZ,3,iBlock)
          !Save pressure tensor: Pxx Pyy Pzz Pxy Pxz Pyz
          State_VCB(6*iSort+11,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,1,1)**2
          State_VCB(6*iSort+12,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,2,1)**2
          State_VCB(6*iSort+13,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,3,1)**2
          State_VCB(6*iSort+14,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,1,1)*&
               V_GDB(1:nX,1:nY,1:nZ,2,1)
          State_VCB(6*iSort+15,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,1,1)*&
               V_GDB(1:nX,1:nY,1:nZ,3,1)
          State_VCB(6*iSort+16,:,:,:,1) = V_GDB(1:nX,1:nY,1:nZ,2,1)*&
               V_GDB(1:nX,1:nY,1:nZ,3,1)
          !Save scalar pressure
          do k=1,nX; do j=1,nY; do i=1,nX
             State_VCB(iSort+14,i,j,k,1) = sum(&
                  State_VCB(6*iSort+11:6*iSort+13,i,j,k,1))/3.0
          end do; end do; end do
       end do
   
       !Collect information from all processors
       if(nProc>1) call mpi_reduce_real_array(State_VCB, &
            size(State_VCB), MPI_SUM, 0, iComm, iError)
       
       if(iProc/=0) RETURN

       ! Processor 0 has all information to be saved
       TypePosition = 'rewind'

       if(nDim==3) then
          allocate(character(len=24+73*nPType) :: NameVar)
          do iSort=1,nPType
             write(NameVar((6*iSort-5):(6*iSort)),'(a4,i1,1x)') 'rhoS',iSort-1
             write(NameVar((15*iSort-14+6*nPType):(15*iSort+6*nPType)),&
                  '(3(a3,i1,1x))') 'uxS',iSort-1,'uyS',iSort-1,'uzS',iSort-1
             write(NameVar((4*iSort-3+21*nPType):(4*iSort+21*nPType)),&
                  '(a2,i1,1x)') 'pS',iSort-1
             write(NameVar((36*iSort-35+25*nPType):(36*iSort+25*nPType)),&
                  '(6(a4,i1,1x))') 'pXXS',iSort-1,'pYYS',iSort-1,'pZZS',&
                  iSort-1,'pXYS',iSort-1,'pXZS',iSort-1,'pYZS',iSort-1
             write(NameVar((6*iSort-5+61*nPType):(6*iSort+61*nPType)),&
                  '(a4,i1,1x)') 'Q_PS',iSort-1
             write(NameVar((6*iSort-5+67*nPType):(6*iSort+67*nPType)),&
                  '(a4,i1,1x)') 'M_PS',iSort-1
          end do
          NameVar = 'x y z Ex Ey Ez Bx By Bz '//NameVar
       elseif(nDim==2) then
          !Understand which var to save
       
          NameVar = 'x y Ex Ey Ez Bx By Bz '//NameVar
       end if

    case('NORMAL')
       ! Check if this is time to save; if not, return
       if(nStepOut>=1.and.nStepOutMin<=iStep)then
          if(mod(iStep,nStepOut) /= 0)&
               RETURN
       end if

       TypePosition = 'append'

    case('FINAL')

       TypePosition = 'append'

    case default
       call CON_stop('ALTOR '//NameSub// &
            ': Incorrect value for TypeSave='//trim(TypeSave))
    end select

    do k=1,nZ; do j=1,nY; do i=1,nX
       ! Cell-centered averaged E field                                      
       State_VCB(1,i,j,k,1) = 0.5*(E_GDB(i,j,k,1,1)+E_GDB(i-1,j,k,1,1))
       State_VCB(2,i,j,k,1) = 0.5*(E_GDB(i,j,k,2,1)+E_GDB(i,j-1,k,2,1))
       State_VCB(3,i,j,k,1) = 0.5*(E_GDB(i,j,k,3,1)+E_GDB(i,j,k-1,3,1))

       ! Cell-centered averaged B field
       State_VCB(4,i,j,k,1) = 0.25*sum(magnetic_GDB(i,j-1:j,k-1:k,1,1))
       State_VCB(5,i,j,k,1) = 0.25*sum(magnetic_GDB(i-1:i,j,k-1:k,2,1))
       State_VCB(6,i,j,k,1) = 0.25*sum(magnetic_GDB(i-1:i,j-1:j,k,3,1))
    end do; end do; end do

    !Save the charge and mass for each species
    !The unit may be a problem
    Param_I(1:nPType)          = Q_P
    Param_I(nPType+1:2*nPType) = M_P

    call save_plot_file(NameFile,     &
         TypeFileIn='real4',          &
         TypePositionIn=TypePosition, &
         nStepIn = iStep, &
         TimeIn  = tSimulation, &
         ParamIn_I = Param_I, &
         NameVarIn = NameVar, &
         CoordMinIn_D = 0.5*Dx_D, &
         CoordMaxIn_D = (nCell_D - 0.5)*Dx_D, &
         VarIn_VIII = State_VCB(:,:,:,:,iBlock))


  end subroutine PIC_save_files

end module PIC_ModOutput
