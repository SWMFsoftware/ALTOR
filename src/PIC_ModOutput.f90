module PIC_ModOutput
  use PIC_ModMain,ONLY: iStep,tSimulation,Dt,DX_D 
  use PIC_ModSize,ONLY: nDim, nCell_D, nX, nY, nZ, x_, y_, z_
  use PIC_ModSize,ONLY: MaxBlock, nPType
  use PIC_ModParticles,ONLY: M_P,Q_P
  use ModNumConst,Only: cPi
  use PIC_ModParticleInField,ONLY: Rho_GB, V_GDB
  use PIC_ModField,   ONLY: iGCN, E_GDB,B_GDB
  use PIC_ModMpi,     ONLY: pass_density, pass_velocity
  use PIC_ModProc,    ONLY: iProc
  use ModIoUnit,      ONLY: io_unit_new
  use PIC_ModLogFile, ONLY: nToWrite     !Number of test particles

  implicit none
  SAVE
  PRIVATE !Except  

  !\                                                                          
  !Unit for the log file                                                      
  !/                           
  integer :: iOutUnit = -1,iOutUnit_den=-1
  integer,  public :: nStepOut=100000,nStepOutMin=100000
  character(len=5),public :: TypeFile='real4' 
  !Array for saving coordinates, fields and moments in one timestep
  real   :: State_VCB(6+11*nPType,1:nX,1:nY,1:nZ,MaxBlock)=0.0

  public :: write_moments

  public :: PIC_save_files
  
contains
  !======================================================================
  subroutine PIC_save_files(TypeSave)
    use ModPlotFile, ONLY: save_plot_file
    use PIC_ModParticles, ONLY: nPType,n_P
    use PC_BATL_particles,ONLY: Particle_I
    use ModMpi, ONLY: mpi_reduce_real_array, MPI_SUM

    character(len=*), intent(in):: TypeSave
    character(len=6) :: TypePosition
    real    :: Param_I(1:2*nPType)
    integer :: iSort
    character(len=15) :: NameFile='variable.outs'

    real,dimension(nDim):: X_D
    real,dimension(nDim):: V_D
    integer :: iBlock=1

    character(len=:), allocatable, save:: NameVar
    character(len=*), parameter:: NameSub = 'save_files'
    !----------------------------------------------------------------------
    select case(TypeSave)
    case('INITIAL')
       !Do not save in test particle runs
       if(nToWrite/=0) RETURN

       !Calculate the initial moments
       call compute_moments

       if(iProc==0.and.nDim==3) then
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
                  '(a4,i1,1x)') 'QS',iSort-1
             write(NameVar((6*iSort-5+67*nPType):(6*iSort+67*nPType)),&
                  '(a4,i1,1x)') 'MS',iSort-1
          end do
          NameVar = 'x y z Ex Ey Ez Bx By Bz '//NameVar
          TypePosition = 'rewind'
       elseif(nDim==2) then
          NameVar = 'x y Ex Ey Ez Bx By Bz '//NameVar
       end if

    case('NORMAL')
       !Check if this is time to save; if not, return
       !The calculation is done in advance_particles
       if(nStepOut<1.or.nStepOutMin>=iStep.or.mod(iStep,nStepOut)/=0)&
            RETURN

       !Do not save in test particle runs
       if(nToWrite/=0) RETURN

       TypePosition = 'append'

    case('FINAL')
       !Do not save in test particle runs
       if(nToWrite/=0) RETURN

       !Calculate the moments at final timestep
       call compute_moments

       TypePosition = 'append'

    case default
       call CON_stop('ALTOR '//NameSub// &
            ': Incorrect value for TypeSave='//trim(TypeSave))
    end select

    if(iProc==0) then
       !Get cell-center averaged E and B; fields are the same on each Proc
       call average_fields
       !Save the charge and mass for each species
       Param_I(1:nPType)          = Q_P
       Param_I(nPType+1:2*nPType) = M_P

       call save_plot_file(NameFile,     &
            TypeFileIn = TypeFile,         &
            TypePositionIn = TypePosition, &
            nStepIn = iStep, &
            TimeIn  = tSimulation, &
            ParamIn_I = Param_I, &
            NameVarIn = NameVar, &
            CoordMinIn_D = 0.5*Dx_D, &
            CoordMaxIn_D = (nCell_D - 0.5)*Dx_D, &
            VarIn_VIII = State_VCB(:,:,:,:,iBlock))
    end if
  end subroutine PIC_save_files
  !======================================================================
  !Calculate cell-centered number density and velocity for output
  subroutine compute_moments
    use PIC_ModParticles, ONLY: nPType,n_P
    use PIC_ModFormFactor,ONLY: HighFF_ID,Node_D, get_form_factors
    use PC_BATL_particles,ONLY: Particle_I
    use PIC_ModParticleInField,ONLY: add_DensityVelocity
    use PIC_ModMain, ONLY: c,c2
    use ModMpi, ONLY: mpi_reduce_real_array, MPI_SUM

    integer :: iSort,iParticle, i
    real,dimension(nDim):: X_D
    real,dimension(nDim):: V_D
    integer :: iBlock=1

    character(len=:), allocatable, save:: NameVar
    !----------------------------------------------------
    !Calculate the number densities and velocities                   
    do iSort = 1, nPType
       Rho_GB = 0.0 ; V_GDB = 0.0
       do iParticle = 1, n_P(iSort)
          !Get the position of particle                                      
          X_D = Particle_I(iSort)%State_VI(1:nDim,iParticle)
          !Get the velocity of particle from generalized momentum            
          V_D = Particle_I(iSort)%State_VI(nDim+1:nDim+3,iParticle)
          V_D = c/(sqrt(c2+sum(V_D**2)))*V_D

          call get_form_factors(X_D,Node_D,HighFF_ID)

          call add_DensityVelocity(V_D,Node_D,HighFF_ID,iBlock)
       end do
       
       !Save number density, velocity and pressure                           
       call write_moments(iSort)
    end do
  end subroutine compute_moments
  !=====================================================================
  !Save the number density, velocity and pressure at certain timestep
  subroutine write_moments(iSort)
    use ModMpi, ONLY: mpi_reduce_real_array, MPI_SUM
    use PIC_ModProc, ONLY: iComm, nProc

    integer, intent(in) :: iSort
    integer :: i,j,k
    integer :: iError
    real    :: U_GDB(1:nDim,1:nX,1:nY,1:nZ,1)  !cell-centered bulk velocity
    !------------------------------------------------------------------
    !Collect info from all processor to processor 0 when in parallel
    !assuming periodic BC.
    call pass_density(0)

    !Periodic velocity BC. Need to be merged into other function later.
    call pass_velocity(0)

    if(iProc==0)then
       !The cell-centered velocity should be normalized by number density
       !After applying periodic BC.
       do i=1,3
          V_GDB(i,1:nX,1:nY,1:nZ,1) = V_GDB(i,1:nX,1:nY,1:nZ,1)&
               / Rho_GB(1:nX,1:nY,1:nZ,1)
       end do

       State_VCB(iSort+6,:,:,:,1) = Rho_GB(1:nX,1:nY,1:nZ,1)

       !hyzhou: need to check if this is correct!
       do i=1,nDim
          U_GDB(i,:,:,:,1) = sum(V_GDB(i,:,:,:,1))/product(nCell_D)
       end do

       !Looping over all the cell centers
       do k=1,nZ; do j=1,nY; do i=1,nX
          !Save velocity
          State_VCB(nPType+3*iSort+4:nPType+3*iSort+6,i,j,k,1) = &
               V_GDB(:,i,j,k,1)
          !Save pressure tensor: Pxx Pyy Pzz Pxy Pxz Pyz
          State_VCB(4*nPType+6*iSort+1:4*nPType+6*iSort+3,i,j,k,1) = &
               M_P(iSort)*Rho_GB(i,j,k,1)*&
               (V_GDB(:,i,j,k,1)**2 - U_GDB(:,i,j,k,1)**2)

          State_VCB(4*nPType+6*iSort+4,:,:,:,1) = M_P(iSort)*Rho_GB(i,j,k,1)*&
               (V_GDB(1,i,j,k,1)*V_GDB(2,i,j,k,1) - U_GDB(1,i,j,k,1)*&
               U_GDB(2,i,j,k,1))
          State_VCB(4*nPType+6*iSort+5,:,:,:,1) = M_P(iSort)*Rho_GB(i,j,k,1)*&
               (V_GDB(1,i,j,k,1)*V_GDB(3,i,j,k,1) - U_GDB(1,i,j,k,1)*&
               U_GDB(3,i,j,k,1))
          State_VCB(4*nPType+6*iSort+6,:,:,:,1) = M_P(iSort)*Rho_GB(i,j,k,1)*&
               (V_GDB(2,i,j,k,1)*V_GDB(3,i,j,k,1) - U_GDB(2,i,j,k,1)*&
               U_GDB(3,i,j,k,1))
       end do; end do; end do

       !Save scalar pressure
       do k=1,nZ; do j=1,nY; do i=1,nX
          State_VCB(10*nPType+iSort+6,i,j,k,1) = sum(&
               State_VCB(4*nPType+6*iSort+4:4*nPType+6*iSort+6,i,j,k,1))/3.0
       end do; end do; end do
    end if
    !hyzhou: I think I need to remove this!!!
    !pass_density and pass_velocity have already done the mpi_reduce
    !so you don`t need to do it again here

    !Collect information from all processors after looping over all
    !particle species
    !if(nProc>1.and.iSort==nPType) call mpi_reduce_real_array(State_VCB, &
    !     size(State_VCB), MPI_SUM, 0, iComm, iError)
    !Processor 0 has all the information to be saved

  end subroutine write_moments
  !======================================================================
  subroutine average_fields
    integer :: i,j,k
    integer :: iBlock=1
    !--------------------
    do k=1,nZ; do j=1,nY; do i=1,nX
       ! Cell-centered averaged E field
       State_VCB(1,i,j,k,iBlock) = 0.5*(E_GDB(i,j,k,1,1)+E_GDB(i-1,j,k,1,1))
       State_VCB(2,i,j,k,iBlock) = 0.5*(E_GDB(i,j,k,2,1)+E_GDB(i,j-1,k,2,1))
       State_VCB(3,i,j,k,iBlock) = 0.5*(E_GDB(i,j,k,3,1)+E_GDB(i,j,k-1,3,1))
       ! Cell-centered averaged B field
       State_VCB(4,i,j,k,iBlock) = 0.25*sum(B_GDB(i,j-1:j,k-1:k,1,1))
       State_VCB(5,i,j,k,iBlock) = 0.25*sum(B_GDB(i-1:i,j,k-1:k,2,1))
       State_VCB(6,i,j,k,iBlock) = 0.25*sum(B_GDB(i-1:i,j-1:j,k,3,1))
    end do; end do; end do
  end subroutine average_fields

end module PIC_ModOutput
