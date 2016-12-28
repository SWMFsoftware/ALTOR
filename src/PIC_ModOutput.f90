module PIC_ModOutput
  use PIC_ModMain,ONLY: iStep,tSimulation,Dt,DX_D, vInv, DxInv_D 
  use PIC_ModSize,ONLY: nDim, nCell_D, nX, nY, nZ, x_, y_, z_
  use PIC_ModSize,ONLY: MaxBlock, nPType
  use PIC_ModParticles,ONLY: M_P,Q_P
  use PIC_ModParticleInField,ONLY: State_VGBI, add_moments
  use PIC_ModField,   ONLY: iGCN, E_GDB,B_GDB
  use PIC_ModMpi,     ONLY: pass_moments
  use PIC_ModProc,    ONLY: iProc
  use ModIoUnit,      ONLY: io_unit_new
  use PIC_ModLogFile, ONLY: nToWrite     !Number of test particles
  use PC_BATL_lib,    ONLY: CoordMin_D, CoordMax_D, CoordMin_DB

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
  real   :: PlotVar_VC(6+11*nPType,1:nX,1:nY,1:nZ)=0.0

  public :: write_moments

  public :: PIC_save_files
  
contains
  !======================================================================
  subroutine PIC_save_files(TypeSave)
    use ModPlotFile, ONLY: save_plot_file
    use PIC_ModParticles, ONLY: nPType!,n_P
    use PIC_ModSize, ONLY: nCell_D
    use PC_BATL_particles,ONLY: Particle_I
    use ModMpi, ONLY: mpi_reduce_real_array, MPI_SUM

    character(len=*), intent(in):: TypeSave
    character(len=6) :: TypePosition
    real    :: Param_I(1:2*nPType)
    integer :: iSort
    character(len=15) :: NameFile='variables.outs'
    character(len=*),parameter :: FileDir='PC/plots/'

    real,dimension(nDim):: X_D
    real,dimension(nDim):: V_D
    integer :: iBlock=1

    character(len=:), allocatable, save:: NameVar
    character(len=*), parameter:: NameSub = 'save_files'
    !----------------------------------------------------------------------
    select case(TypeSave)
    case('INITIAL')
       !Do not save in test particle runs alone (test particle number>1 and
       !density of electrons < 1 per cell)
       if(nToWrite/=0.and.Particle_I(1)%nParticle<product(nCell_D)) RETURN

       !Calculate the initial moments
       call compute_moments

       if(iProc==0.and.nDim==3) then
          allocate(character(len=69*nPType) :: NameVar)
          do iSort=1,nPType
             write(NameVar((6*iSort-5):(6*iSort)),'(a4,i1,1x)') 'rhoS',iSort-1
             write(NameVar((15*iSort-14+6*nPType):(15*iSort+6*nPType)),&
                  '(3(a3,i1,1x))') 'uxS',iSort-1,'uyS',iSort-1,'uzS',iSort-1
             write(NameVar((4*iSort-3+21*nPType):(4*iSort+21*nPType)),&
                  '(a2,i1,1x)') 'pS',iSort-1
             write(NameVar((36*iSort-35+25*nPType):(36*iSort+25*nPType)),&
                  '(6(a4,i1,1x))') 'pXXS',iSort-1,'pYYS',iSort-1,'pZZS',&
                  iSort-1,'pXYS',iSort-1,'pXZS',iSort-1,'pYZS',iSort-1
             write(NameVar((4*iSort-3+61*nPType):(4*iSort+61*nPType)),&
                  '(a2,i1,1x)') 'QS',iSort-1
             write(NameVar((4*iSort-3+65*nPType):(4*iSort+65*nPType)),&
                  '(a2,i1,1x)') 'MS',iSort-1
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
       if(nToWrite/=0.and.Particle_I(1)%nParticle<product(nCell_D)) RETURN

       TypePosition = 'append'

    case('FINAL')
       !Do not save in test particle runs alone
       if(nToWrite/=0.and.Particle_I(1)%nParticle<product(nCell_D)) RETURN

       !Calculate the moments at final timestep
       call compute_moments

       TypePosition = 'append'

    case default
       call CON_stop('ALTOR '//NameSub// &
            ': Incorrect value for TypeSave='//trim(TypeSave))
    end select

    if(iProc==0) then
       !Get cell-center averaged E and B; fields are the same on each Proc
       call average_fields(iBlock)
       !Save the charge and mass ratio for each species in normalized unit
       Param_I(1:nPType)          = Q_P
       Param_I(nPType+1:2*nPType) = M_P
              
       call save_plot_file(FileDir//NameFile,     &
            TypeFileIn = TypeFile,         &
            TypePositionIn = TypePosition, &
            nStepIn = iStep, &
            TimeIn  = tSimulation, &
            ParamIn_I = Param_I, &
            NameVarIn = NameVar, &
            CoordMinIn_D = 0.5*Dx_D, &
            CoordMaxIn_D = (nCell_D - 0.5)*Dx_D, &
            VarIn_VIII = PlotVar_VC(:,:,:,:))
    end if
  end subroutine PIC_save_files
  !======================================================================
  !Calculate cell-centered number density and velocity for output
  subroutine compute_moments
    use PIC_ModParticles, ONLY: nPType!,n_P
    use PIC_ModFormFactor,ONLY: HighFF_ID,Node_D, get_form_factors
    use PC_BATL_particles,ONLY: Particle_I
    !use PIC_ModParticleInField,ONLY: add_moments
    use PIC_ModMain, ONLY: c,c2
    use ModMpi, ONLY: mpi_reduce_real_array, MPI_SUM

    integer :: iSort,iParticle, i
    real,dimension(nDim):: X_D
    real,dimension(nDim):: V_D
    integer :: iBlock

    character(len=:), allocatable, save:: NameVar
    !----------------------------------------------------
    !Calculate the number densities and velocities                   
    do iSort = 1, nPType
       State_VGBI(:,:,:,:,:,iSort) = 0.0
       do iParticle = 1, Particle_I(1)%nParticle
          !\
          ! Get iBlock
          !/
          iBlock = Particle_I(iSort)%iIndex_II(0,iParticle)
          !Get the position of particle                                      
          X_D = (Particle_I(iSort)%State_VI(1:nDim,iParticle) - &
               CoordMin_DB(1:nDim,iBlock))*DxInv_D
          !Get the velocity of particle from generalized momentum            
          V_D = Particle_I(iSort)%State_VI(nDim+1:nDim+3,iParticle)
          V_D = c/(sqrt(c2+sum(V_D**2)))*V_D

          call get_form_factors(X_D,Node_D,HighFF_ID)

          call add_moments(V_D,Node_D,HighFF_ID,iBlock,iSort)
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
    !------------------------------------------------------------------
   
    !Periodic velocity BC. Need to be merged into other function later.
    call pass_moments(0,iSort)

    if(iProc==0)then
       !The cell-centered velocity should be normalized by number density
       !After applying periodic BC.
       do i=2,4
          State_VGBI(i,1:nX,1:nY,1:nZ,1,iSort) = &
              State_VGBI(i,1:nX,1:nY,1:nZ,1,iSort)/State_VGBI(1,1:nX,1:nY,1:nZ,1,iSort)
       end do

       PlotVar_VC(iSort+6,:,:,:) = abs(Q_P(iSort))*vInv*State_VGBI(1,1:nX,1:nY,1:nZ,1,iSort)

       !Looping over all the cell centers
       do k=1,nZ; do j=1,nY; do i=1,nX
          !Save velocity
          PlotVar_VC(nPType+3*iSort+4:nPType+3*iSort+6,i,j,k) = &
               State_VGBI(2:4,i,j,k,1,iSort)
          !Save pressure tensor: Pxx Pyy Pzz Pxy Pxz Pyz
          PlotVar_VC(5*nPType+6*iSort+1:5*nPType+6*iSort+3,i,j,k) = &
               M_P(iSort)*vInv*(State_VGBI(5:7,i,j,k,1,iSort) -&
               State_VGBI(1,i,j,k,1,iSort)*&
               State_VGBI(2:4,i,j,k,1,iSort)**2)
          
          PlotVar_VC(5*nPType+6*iSort+4,:,:,:) = M_P(iSort)*vInv*(State_VGBI(8,i,j,k,1,iSort)&
               - State_VGBI(1,i,j,k,1,iSort)*&
               State_VGBI(2,i,j,k,1,iSort)*State_VGBI(3,i,j,k,1,iSort))
          PlotVar_VC(5*nPType+6*iSort+5,:,:,:) = M_P(iSort)*vInv*(State_VGBI(9,i,j,k,1,iSort)&
               - State_VGBI(1,i,j,k,1,iSort)*&
               State_VGBI(2,i,j,k,1,iSort)*State_VGBI(4,i,j,k,1,iSort))
          PlotVar_VC(5*nPType+6*iSort+6,:,:,:) = M_P(iSort)*vInv*(State_VGBI(10,i,j,k,1,iSort)&
               - State_VGBI(1,i,j,k,1,iSort)*&
               State_VGBI(3,i,j,k,1,iSort)*State_VGBI(4,i,j,k,1,iSort))
       end do; end do; end do

       !Save scalar pressure
       do k=1,nZ; do j=1,nY; do i=1,nX
          PlotVar_VC(4*nPType+iSort+6,i,j,k) = sum(&
               PlotVar_VC(5*nPType+6*iSort+1:5*nPType+6*iSort+3,i,j,k))/3.0
       end do; end do; end do
    end if

  end subroutine write_moments
  !======================================================================
  subroutine average_fields(iBlock)
    integer :: i,j,k
    integer,intent(in) :: iBlock
    !--------------------
    do k=1,nZ; do j=1,nY; do i=1,nX
       ! Cell-centered averaged E field
       PlotVar_VC(1,i,j,k) = &
            0.5*(E_GDB(i,j,k,1,iBlock)+E_GDB(i-1,j,k,1,iBlock))
       PlotVar_VC(2,i,j,k) = &
            0.5*(E_GDB(i,j,k,2,iBlock)+E_GDB(i,j-1,k,2,iBlock))
       PlotVar_VC(3,i,j,k) = &
            0.5*(E_GDB(i,j,k,3,iBlock)+E_GDB(i,j,k-1,3,iBlock))
       ! Cell-centered averaged B field
       PlotVar_VC(4,i,j,k) = 0.25*sum(B_GDB(i,j-1:j,k-1:k,1,iBlock))
       PlotVar_VC(5,i,j,k) = 0.25*sum(B_GDB(i-1:i,j,k-1:k,2,iBlock))
       PlotVar_VC(6,i,j,k) = 0.25*sum(B_GDB(i-1:i,j-1:j,k,3,iBlock))
    end do; end do; end do
  end subroutine average_fields

end module PIC_ModOutput
