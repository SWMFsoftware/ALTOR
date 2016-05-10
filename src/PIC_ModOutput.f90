module PIC_ModOutput
  use PIC_ModMain,ONLY: iStep,tSimulation,dt,dX_D 
  use PIC_ModSize,ONLY: nDim, nCell_D, nX, nY, nZ
  use PIC_ModParticleInField,ONLY: Rho_GB
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
                                               
  integer :: iOutUnit = -1
  integer, public :: nStepOut=100000,nStepOutMin=100000
  public :: write_density,write_field
 
  character(LEN=30) ::NameFormat
contains
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
  !
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
  !
end module PIC_ModOutput
