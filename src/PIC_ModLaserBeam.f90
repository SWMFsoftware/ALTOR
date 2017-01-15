module PIC_ModLaserBeam
  use PC_ModSize,  ONLY: nDim, x_, y_, z_, MaxDim
  use PC_BATL_lib, ONLY: CoordMin_D, CoordMax_D
  use PIC_ModMain, ONLY: tSimulation
  use ModNumConst, ONLY: cPi, cTwoPi, cDegToRad
  use PIC_ModProc,      ONLY:iProc
  implicit none
  SAVE
  PRIVATE
  !LaserAmplitude in focus
  !Polarization {Ey/=0 Ez=0} p-pol, {Ey=0 Ez/=0} s-pol 
  !Polarization {Ey/=0 Ez/=0} circular, elliptical
  real :: laserAmplitude(2:3)=0.0 !amplitude * polarization for Ey,Ez

  real :: phaseShift=0.0 !Ey(center)=cos(0), Ez(center)=sin(0)
  !\
  ! Upon reading, the transformation may be done as follows:
  ! if(PulseWidthFocus_D(2)<=0.0.or. PulseWidthFocus_D(3)<=0) then   
  !           PulseWidthFocus_D(2:3)=0.42*XyzFocus_D1(1)
  ! which corresponds to focusing in 1 wavelength        
  !/
  real :: PulseWidthFocus_D(nDim)= 1. !In cycles, wavelengths
  real :: XyzFocus_D(nDim)= 1.            
  real :: coeff=-1.!Gaussian distribution: Intensity/2 at pulse width          
  real :: timePulseBegin=0.0,timePulseEnd=0.0,xPulseCenter=0.0
  !
  !Pulse longitudinal envelope: Default-Steplike, 
  !              1-Gauss, 2-Cosine (FWHM - full width at half maximum)
  integer :: nEnvelope = 1
  real:: PulseWidthBoundary_D(1:nDim)
  public::read_laser_beam,laser_beam
  !                                 
contains
  subroutine laser_beam(iDir, Xyz_D, EField)
    !\
    ! The electric field parameter EField has an intent inout
    ! The input value is found from 'noreflect' boundary condition 
    !/
    !\
    ! Calculate electric field iDir component in the point, x,y,z 
    !/
    integer,intent(in)        :: iDir
    real,   intent(in)        :: Xyz_D(MaxDim)
    real,   intent(inout)     :: EField !Before 
    !-------------------------
    real :: AmplitudeFactor_D(nDim)
    real :: PulseWidthLocal_D(2:nDim)
    real :: focusDist,timeReal,tmp
    real :: laser_profile,laser_phase
    !==========================================
    !\
    ! Works for the left boundary along x only. 
    ! Beam is parallel to x-axis
    !/
    !
    if(tSimulation > timePulseEnd) RETURN
    focusDist=sqrt(sum((XyzFocus_D(1:nDim)-Xyz_D(1:nDim))**2))
    timeReal=tSimulation+focusDist-XyzFocus_D(1)
    AmplitudeFactor_D(:)=0.0
    !
    if(timeReal<=timePulseEnd.and.timeReal>=timePulseBegin) then
       !
       laser_phase=focusDist-XyzFocus_D(1)+(xPulseCenter+tSimulation)
       !
       !==========================along X
       if(nEnvelope==1) then 
          tmp=( laser_phase / (PulseWidthFocus_D(1)*coeff) )**2
          if(tmp < 30.) AmplitudeFactor_D(1)=exp(-tmp)
       else if(nEnvelope==2) then
          AmplitudeFactor_D(1)=cos(laser_phase/PulseWidthFocus_D(1)*(cPi/2.0))
       else
          AmplitudeFactor_D(1)=1.0
       end if
       !==========================along Y,Z
       !====de facto pulse width at the boundary
       PulseWidthLocal_D(2:nDim)=PulseWidthBoundary_D(2:nDim)&
            *(XyzFocus_D(1)-Xyz_D(x_))/XyzFocus_D(1)
       !==========================along Y
       !tmp=( (Xyz_D(2)-XyzFocus_D(2)) / (PulseWidthLocal_D(2)/2.0) )**2
       tmp=( (Xyz_D(2)-XyzFocus_D(2)) / (PulseWidthLocal_D(2)*coeff) )**2
       if(tmp < 30.) AmplitudeFactor_D(2)=exp(-tmp) &
            *PulseWidthFocus_D(2)/PulseWidthLocal_D(2) !due to focusing
       
       !==========================along Z for 3D 
       if(nDim==3) then
          tmp=((Xyz_D(nDim)-XyzFocus_D(nDim))/(PulseWidthLocal_D(nDim)*coeff))**2
          if(tmp < 30.) AmplitudeFactor_D(nDim)=exp(-tmp) &
               *PulseWidthFocus_D(nDim)/PulseWidthLocal_D(nDim) !due to focusing
       end if
       !tmp=(-laser_phase *2.0/PulseWidth_D(1))**2
       !if(tmp < 30.) AmplitudeFactor_D(1)=exp(-tmp)
       !tmp=(focusDist*asin((Xyz_D(2)-XyzFocus_D(2))/focusDist)&
       !     *2.0/PulseWidthFocus_D(2))**2
       !if(tmp < 30.) AmplitudeFactor_D(2)=exp(-tmp)
       !if(nDim==3)then
       !   tmp=(focusDist*asin((Xyz_D(nDim)-XyzFocus_D(nDim))/focusDist)&
       !        *2.0/PulseWidthFocus_D(nDim))**2
       !   if(tmp < 30.) AmplitudeFactor_D(nDim)=exp(-tmp)
       !end if
       !
       laser_profile=product(AmplitudeFactor_D(1:nDim))*laserAmplitude(iDir)
       laser_phase=laser_phase+phaseShift
       !
       if(iDir==y_) then
          EField=laser_profile*cos(laser_phase)
       end if
       if(iDir==z_) then
          EField=laser_profile*sin(laser_phase)
       end if
       !
    end if
  end subroutine laser_beam
  !                                
  !================================================                          
  !=========reading command #LASERBEAM=============                       
  !================================================                        
  subroutine read_laser_beam
    use ModReadParam,ONLY: read_var
    
    !======================================READING BEGIN          
    !LaserAmplitude in focus
    !Polarization {Ey/=0 Ez=0} p-pol, {Ey=0 Ez/=0} s-pol 
    !Polarization {Ey/=0 Ez/=0} circular, elliptical
    call read_var('laserAmplitude2',laserAmplitude(2))
    call read_var('laserAmplitude3',laserAmplitude(3))
    !
    !PhaseShift in respect of Ey(xCenter)=cos(0), Ez(xCenter)=sin(0); 
    !PhaseShift is given in fractions of +-1, to be * (2*cPi).
    call read_var('phaseShift',phaseShift)
    !
    !PulseWidth in focus is given in cycles/wavelengths  
    !PulseWidthFocus_D(2:3)==f/D
    !f/D==[focal radius...]/[diameter... of the focusing parabola] 
    !== N=0.94 (56^degrees), N=1., N=1.5, N=2,...

    call read_var('PulseWidthFocus_D1',PulseWidthFocus_D(1))
    call read_var('PulseWidthFocus_D2',PulseWidthFocus_D(2))
    if(nDim==3)call read_var('PulseWidthFocus_D3',PulseWidthFocus_D(nDim))
    !
    !Focal point: 
    !if XyzFocus_D(:) >=0 - no change
    !if -1<XyzFocus_D(:)<0 - fractions of [region size]
    !if XyzFocus_D(:)<=-1 - in wavelengths
    ! XyzFocus_D(2:3) should be given in the same range 
    call read_var('XyzFocus_D1',XyzFocus_D(1))
    call read_var('XyzFocus_D2',XyzFocus_D(2))
    if(nDim==3)call read_var('XyzFocus_D3',XyzFocus_D(nDim))
    !
    !Pulse longitudinal envelope: Default-Steplike, 
    !              1-Gauss, 2-Cosine (FWHM - full width at half maximum) 
    !Transverse envelope - Gauss, exp(-2) at +-Width/2 
    call read_var('nEnvelope',nEnvelope)
    !                                                                         
    timePulseBegin=0.0
    xPulseCenter=-PulseWidthFocus_D(1)
    timePulseEnd=PulseWidthFocus_D(1)*2.0
    !                                                                        
    coeff=1.0/sqrt(2.0*alog(2.0)) !coeff=0.849322                           
    !                                                                        
    if(iProc==0)then
       write(*,'(a,es13.5)')'Laser pulse of duration [cycles]=',PulseWidthFocus_D(1)
!       if(XyzFocus_D(1) <=0.0 ) then                                           
       if(PulseWidthFocus_D(2)<=0.0.or.PulseWidthFocus_D(nDim)<=0) then
          write(*,'(a,es13.5)')'is focused to a wavelength at x=',XyzFocus_D(1)
       else
          write(*,'(a,2es13.5)')'and of width [wavelengths]=',PulseWidthFocus_D(2:nDim)
       end if
    end if
    !                                                                         
    if(PulseWidthFocus_D(2)<=0.0.or.PulseWidthFocus_D(nDim)<=0) then
       PulseWidthFocus_D(2:nDim)=XyzFocus_D(1)*sin(45.0*cDegToRad)
    end if
    !                                                                    
  end subroutine read_laser_beam
  !                  
  !
end module PIC_ModLaserBeam
