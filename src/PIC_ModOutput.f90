module PIC_ModOutput
  use PIC_ModMain,ONLY: iStep,tSimulation,dt,dX_D 
  use PIC_ModSize,ONLY: nDim, nCell_D, nX, nY, nZ
  use PIC_ModParticleInField,ONLY: rho_G
  use PIC_ModField, ONLY: iGCN, E_GD,Magnetic_GD
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
    rmax=maxval(rho_G)
    if(iProc==0.and.rmax>1.e-10)then
        write(Name,'(a,i1,a,i4.4,a)')'n',iSort,'_',iStep-1,'.dat'
        iOutUnit = io_unit_new()
        open(iOutUnit,file=trim(Name),status='replace',form='unformatted')
        !open(iOutUnit,file=trim(Name),status='replace') !,form='formatted')
        !write(iOutUnit,'(6i10,7e13.5)') &
        write(iOutUnit) &
             iStep,nDim, iGCN, nCell_D,tSimulation-dt,dt,Dx_D &
             ,minval(rho_G),rmax
        !write(iOutUnit) rho_G(:,:,:) 
        !write(iOutUnit,'(8es13.5)')rho_G(:,:,:)
        write(iOutUnit) rho_G(:,:,nZ/2)
        write(iOutUnit) rho_G(:,nY/2,:)
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
        !open(iOutUnit,file=trim(Name),status='replace') !,form='formatted')    
        !write(iOutUnit,'(6i10,7e13.5)') &                                      
        write(iOutUnit) &
             iStep,nDim, iGCN, nCell_D,tSimulation,dt,Dx_D &
             !,minval(E_GD(1:nX,1:nY,1:nZ,iX)),maxval(E_GD(:,:,:,iX))
             ,minval(E_GD(1:nX,1:nY,1:nZ,iX)),maxval(E_GD(:,:,:,iX))
        write(iOutUnit) E_GD(:,:,nZ/2,iX)
        write(iOutUnit) E_GD(:,nY/2,:,iX) 
       !write(iOutUnit,'(8es13.5)')rho_G(:,:,:)                                
        close(iOutUnit)
     end if
  end subroutine write_field  
  !
end module PIC_ModOutput
