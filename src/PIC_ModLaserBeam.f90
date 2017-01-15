!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==========================
module PIC_ModLaserBeam
  use PC_ModSize,  ONLY: nDim, x_, y_, z_, MaxDim
  use PC_BATL_lib, ONLY: CoordMin_D, CoordMax_D
  use PIC_ModMain, ONLY: tSimulation
  use ModNumConst, ONLY: cPi, cTwoPi, cDegToRad
  use PIC_ModProc,      ONLY:iProc
  implicit none
  SAVE
  PRIVATE
  !LaserAmplitude in focus (This comment is only in ModLaserBeam, while
  !       in ModFiledLaser.f90 the amplitude is rather at the boundary)
  !\
  !Polarization {Ey/=0 Ez=0} p-pol, {Ey=0 Ez/=0} s-pol 
  !Polarization {Ey/=0 Ez/=0} circular, elliptical
  !/
  real :: laserAmplitude(2:3)=0.0 !amplitude * polarization for Ey,Ez
  !\
  ! Phase shift, such that at the initial time instant and in the
  ! central point of the pulse the field is proportional to
  ! Ey_center ~ cos(PhaseShift), Ez_center ~ sin(PhaseShift)
  !/
  !\
  ! While reading: In cycles, wavelengths
  ! Upon reading is converted into dimensionless form k*width 
  ! by multiplying by cTwoPi
  !/
  real :: phaseShift=0.0 
  !\
  ! THIS COMMENT IS ONLY VALID FOR ModFieldLaser
  ! Upon reading, the transformation may be done as follows:
  ! if(PulseWidthFocus_D(2)<=0.0.or. PulseWidthFocus_D(3)<=0) then   
  !           PulseWidthFocus_D(2:3)=0.42*XyzFocus_D1(1)
  ! which corresponds to focusing in 1 wavelength   
  ! THE END OF COMMENT VALID FOR  ModFieldLaser    
  !/
  !\
  ! While reading: In cycles, wavelengths
  ! Upon reading is converted into dimensionless form k*width 
  ! by multiplying by cTwoPi
  !/
  real :: PulseWidthFocus_D(nDim)= 1. 
  !\
  ! Dimensionless coordinates of the focal point
  ! While reading, they are set with respect to the 
  ! left corner and may be expressed in terms of
  ! wavelengths or fractions of the computational 
  ! domain size.
  !/
  real :: XyzFocus_D(nDim)= 1.
  !\
  ! WidthHalfMax2Gaussian=1.0/sqrt(2.0*alog(2.0))
  ! If WidthGaussian = WidthHalfMax* WidthHalfMax2Gaussian,
  ! then the relative amplitude in the points 
  ! (+/-)0.5*WidthHalfMax from the beam axis equals 
  ! exp(-(0.5*WidthHalfMax/WidthGaussian)**2)=sqrt(0.5)
  ! so that the intensity in these points is 0.5 of the
  ! axial one, in a full agreement with the definition 
  ! of the Full Width at Half Maximal (Intensity) - FWHM
  !/
  real :: WidthHalfMax2Gaussian=-1.0
  !\
  ! Time moments when the beamed pulse starts and ends
  ! passing through the domain boundary at the beam axis.
  !/          
  real :: timePulseBegin=0.0,timePulseEnd=0.0
  !\
  ! X coordinate of the laser pulse center in the 
  ! initial time instant
  !/

  real :: xPulseCenter=0.0
  !
  !Pulse longitudinal envelope: Default-Steplike, 
  !              1-Gauss, 2-Cosine (FWHM - full width at half maximum)
  integer :: nEnvelope = 1
  !\
  ! Pulse FWHM at the domain boundary (well before focusing)
  !/
  real    :: PulseWidthBoundary_D(1:nDim)
  !\
  ! public members
  !/
  public  :: read_laser_beam !reads the laser pulsed beam parameters
  public  :: laser_beam      !modifies the boundary field for laser beam 
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
    real,   intent(inout)     :: EField !Before and after assignment (if any) 
    !-------------------------
    real :: AmplitudeFactor_D(nDim)
    real :: PulseWidthLocal_D(2:nDim)
    real :: focusDist,timeReal,tmp
    real :: laser_profile,laser_phase
    real:: focalrad(1:nDim)
    integer :: ii
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
          tmp=( laser_phase / (PulseWidthFocus_D(1)*WidthHalfMax2Gaussian) )**2
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
       tmp=( (Xyz_D(2)-XyzFocus_D(2)) / (PulseWidthLocal_D(2)*WidthHalfMax2Gaussian) )**2
       if(tmp < 30.) AmplitudeFactor_D(2)=exp(-tmp) &
            *PulseWidthFocus_D(2)/PulseWidthLocal_D(2) !due to focusing
       
       !==========================along Z for 3D 
       if(nDim==3) then
          tmp=((Xyz_D(nDim)-XyzFocus_D(nDim))/(PulseWidthLocal_D(nDim)*WidthHalfMax2Gaussian))**2
          if(tmp < 30.) AmplitudeFactor_D(nDim)=exp(-tmp) &
               *PulseWidthFocus_D(nDim)/PulseWidthLocal_D(nDim) !due to focusing
       end if
       !Version of PIC_ModFieldLaser.f90
       !tmp=(-laser_phase *2.0/PulseWidth)**2
       !if(tmp < 30.) AmplitudeFactor_D(1)=exp(-tmp)
       !tmp=(focusDist*asin((Xyz_D(2)-XyzFocus_D(2))/focusDist)&
       !     *2.0/PulseWidth(2))**2
       !if(tmp < 30.) AmplitudeFactor_D(2)=exp(-tmp)
       !if(nDim==3)then
       !   tmp=(focusDist*asin((Xyz_D(nDim)-XyzFocus_D(nDim))/focusDist)&
       !        *2.0/PulseWidth(nDim))**2
       !   if(tmp < 30.) AmplitudeFactor_D(nDim)=exp(-tmp)
       !end if
       !
       laser_profile=product(AmplitudeFactor_D(1:nDim))*laserAmplitude(iDir)
       laser_phase=laser_phase + phaseShift
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
    use PC_ModSize,  ONLY: nCell_D
    use PIC_ModMain, ONLY: Dt, Dx_D
    integer:: ii
    real :: focalrad(2:nDim)
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
    !========================================READING END
    !
    !=======================TRANSFORMATION OF DATA BEGIN
    ! 
    phaseShift=phaseShift*cTwoPi
    PulseWidthFocus_D(:)=PulseWidthFocus_D(:)*cTwoPi    !pulse width in c/\omega
    !
    do ii=1,nDim
       if(XyzFocus_D(ii) < 0.0) then 
          if(XyzFocus_D(ii) <= -1.0) then 
             XyzFocus_D(ii)=-XyzFocus_D(ii)*cTwoPi
          else
             XyzFocus_D(ii)=-XyzFocus_D(ii)*Dx_D(ii)*nCell_D(ii)
          end if
       end if
    end do
    !
    !=========================TRANSFORMATION OF DATA END
    !
    !pulse diameter at x=0 & FocalRadius in transverse direction
    !===================================IGOR====================
    !===========================================================
    !\
    ! I do not understand any of the formulae below
    !/
    PulseWidthBoundary_D(2:nDim)=XyzFocus_D(1)/sqrt((PulseWidthFocus_D(2:nDim)/cTwoPi)**2-0.25)
    !\
    ! Above there is a formula, relating a Gaussian width in the focus with that
    ! at the distance XyzFocus_D(1) from the focus. Where that particlar formula 
    ! is taking from
    !/

    !\
    ! And now:
    !/
    focalRad(2:nDim)=PulseWidthFocus_D(2:nDim)/cTwoPi*PulseWidthBoundary_D(2:nDim)
    !\
    !Why in the formula above focalRad is the focal distance - it does not include 
    !XyzFocus_D(1). And this is the focal distance for what?
    !/

    !\
    ! Now
    !/
    timePulseBegin=0.0 !start of the pulse
    !\
    ! What is the point to nullify the variable, which is reassigned just after?
    !/
    !\
    ! Now:
    !/
    timePulseBegin=maxval(focalRad)-XyzFocus_D(1)
    !\
    ! timePulseBegin is the time, when the pulse front reaches the domain boundary             
    ! 
    !Double pulse width for Gaussian and Cosine envelopes
    if(nEnvelope==1.or.nEnvelope==2) then 
       xPulseCenter=-timePulseBegin-PulseWidthFocus_D(1)
       timePulseEnd=timePulseBegin+PulseWidthFocus_D(1)*2.0
    else 
       xPulseCenter=-timePulseBegin-PulseWidthFocus_D(1)/2.0
       timePulseEnd=timePulseBegin+PulseWidthFocus_D(1)
    end if
    !                             
    !For Gaussian distribution: Intensity/2 at pulse width 
    WidthHalfMax2Gaussian=1.0/sqrt(2.0*alog(2.0)) !WidthHalfMax2Gaussian=0.849322 for Gauss          
    !                                                                        
    if(iProc==0)then
       write(*,'(a,es13.5)')'Laser pulse of duration [1/omega]=',PulseWidthFocus_D(1)
       write(*,'(a,2es13.5)')'is focused to a spot of [1/k] ',&
               PulseWidthFocus_D(2:nDim)
       write(*,'(a,3es13.5)')'  at the focal point Xyz=',XyzFocus_D(1:nDim)
       write(*,'(a,2es13.5)')'PulseWidthBoundary_D=',PulseWidthBoundary_D
       write(*,'(a,2es13.5)')'focalRad=',focalRad
       write(*,'(a,3es13.5)')'timePulseBegin,timePulseEnd,xPulseCenter=',&
            timePulseBegin,timePulseEnd,xPulseCenter
       write(*,'(a,2es13.5,i5)')'>>Pulse should be focused et time,it=',&
            XyzFocus_D(1)-xPulseCenter,floor((XyzFocus_D(1)-xPulseCenter)/dt)
    end if
    !\
    !  Version of PIC_ModFieldLaser.f90 (with some renamings)
    !/
    !timePulseBegin=0.0
    !xPulseCenter=-PulseWidth(1)
    !timePulseEnd=PulseWidth(1)*2.0
    !                                                                        
    !WidthHalfMax2Gaussian=1.0/sqrt(2.0*alog(2.0)) !=0.849322                           
    !                                                                        
    !if(iProc==0)then
    !   write(*,'(a,es13.5)')'Laser pulse of duration [cycles]=',PulseWidth(1)
    !                                           
    !   if(PulseWidthFocus_D(2)<=0.0.or.PulseWidthFocus_D(nDim)<=0) then
    !      write(*,'(a,es13.5)')'is focused to a wavelength at x=',XyzFocus_D(1)
    !   else
    !      write(*,'(a,2es13.5)')'and of width [wavelengths]=',PulseWidth(2:nDim)
    !   end if
    !end if
    !pulseWidth=pulseWidth*WidthHalfMax2Gaussian
    !if(PulseWidth(2)<=0.0.or.PulseWidth(nDim)<=0) then
    !   One version:     PulseWidth(2:nDim) = XyzFocus_D(1)*sin(45.0*cDegToRad)
    !   Another version: PulseWidth(2:nDim) = xFocus(1)*0.42/2.0
    !end if
    !                                                                    
  end subroutine read_laser_beam
  !                  
  !
end module PIC_ModLaserBeam
