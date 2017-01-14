module PIC_ModLaserBeam
  use PC_ModSize, ONLY: nDim, x_, y_, z_
  use PIC_ModMain, ONLY: tSimulation
  use ModNumConst, ONLY: cOne,cHalf,cDegToRad
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
  ! if(PulseWidth_D(2)<=0.0.or. PulseWidth_D(3)<=0) then   
  !           PulseWidth_D(2:3)=0.42*XyzFocus_D1(1)
  ! which corresponds to focusing in 1 wavelength        
  !/                                                            
  real :: PulseWidth_D(nDim)=-1. !In cycles, wavelengths
  real :: XyzFocus_D(nDim)=-1.            
  real :: coeff=-1.!Gaussian distribution: Intensity/2 at pulse width          
  real :: timePulseBegin=0.0,timePulseEnd=0.0,xPulseCenter=0.0
  !                                                                            
  public::read_laser_beam,laser_beam
  integer,public :: laser_pulse=-1
  !                                 
contains
  subroutine laser_beam(iDir, x, y, z, EField)
    !\
    ! The electric field parameter EField has an intent inout
    ! The input value is found from 'noreflect' boundary condition 
    !/
    !\
    ! Calculate electric field iDir component in the point, x,y,z 
    !/
    integer,intent(in)        :: iDir
    real,   intent(in)        :: x,y
    real, optional, intent(in):: z !not used in 2D geometry
    real,   intent(inout)     :: EField !Before 
    !-------------------------
    real :: Xyz_D(nDim)
    real :: pX(nDim)
    real :: focusDist,timeReal,tmp
    real :: laser_profile,laser_phase
    !==========================================
    !
    !Works for the left boundary only
    !
    Xyz_D(1:2)=(/x,y/)
    if(nDim==3)Xyz_D(nDim) = z
    EField=0.0
    !
    if(tSimulation>=timePulseBegin .and. tSimulation <= timePulseEnd) then
       if(x<0.0) then 
          focusDist=sqrt(sum((XyzFocus_D(1:nDim)-Xyz_D(1:nDim))**2))
          timeReal=tSimulation+focusDist-XyzFocus_D(1)
          pX(:)=0.0
          !
          if(timeReal<=timePulseEnd) then
             !
             laser_phase=focusDist-XyzFocus_D(1)+(xPulseCenter+tSimulation)
             !                                                          
             tmp=(-laser_phase *2.0/PulseWidth_D(1))**2
             if(tmp < 30.) pX(1)=exp(-tmp)
             tmp=(focusDist*asin((Xyz_D(2)-XyzFocus_D(2))/focusDist)&
                  *2.0/PulseWidth_D(2))**2
             if(tmp < 30.) pX(2)=exp(-tmp)
             if(nDim==3)then
                tmp=(focusDist*asin((Xyz_D(nDim)-XyzFocus_D(nDim))/focusDist)&
                     *2.0/PulseWidth_D(nDim))**2
                if(tmp < 30.) pX(nDim)=exp(-tmp)
             end if
             !
             laser_profile=product(pX(1:nDim))*laserAmplitude(iDir)
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
       end if
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
    !PulseWidth_D(2:3)==f/D
    !f/D==[focal radius...]/[diameter... of the focusing parabola] 
    !== N=0.94 (56^degrees), N=1., N=1.5, N=2,...

    call read_var('PulseWidth_D1',PulseWidth_D(1))
    call read_var('PulseWidth_D2',PulseWidth_D(2))
    if(nDim==3)call read_var('PulseWidth_D3',PulseWidth_D(nDim))
    !
    !Focal point: 
    !if XyzFocus_D(:) >=0 - no change
    !if -1<XyzFocus_D(:)<0 - fractions of [region size]
    !if XyzFocus_D(:)<=-1 - in wavelengths
    ! XyzFocus_D(2:3) should be given in the same range 
    call read_var('xFocus1',XyzFocus_D(1))
    call read_var('xFocus2',XyzFocus_D(2))
    if(nDim==3)call read_var('xFocus3',XyzFocus_D(nDim))
    !
    laser_pulse=2
    !                                                                         
    timePulseBegin=0.0
    xPulseCenter=-PulseWidth_D(1)
    timePulseEnd=PulseWidth_D(1)*2.0
    !                                                                        
    coeff=1.0/sqrt(2.0*alog(2.0)) !coeff=0.849322                           
    !                                                                        
    if(iProc==0)then
       write(*,'(a,es13.5)')'Laser pulse of duration [cycles]=',PulseWidth_D(1)
!       if(XyzFocus_D(1) <=0.0 ) then                                           
       if(PulseWidth_D(2)<=0.0.or.PulseWidth_D(nDim)<=0) then
          write(*,'(a,es13.5)')'is focused to a wavelength at x=',XyzFocus_D(1)
       else
          write(*,'(a,2es13.5)')'and of width [wavelengths]=',PulseWidth_D(2:nDim)
       end if
    end if
    !                                                                         
    PulseWidth_D=PulseWidth_D*coeff
    if(PulseWidth_D(2)<=0.0.or.PulseWidth_D(nDim)<=0) then
       PulseWidth_D(2:nDim)=XyzFocus_D(1)*sin(45.0*cDegToRad)
       laser_pulse=1
    end if
    !                                                                    
  end subroutine read_laser_beam
  !                  
  !
end module PIC_ModLaserBeam
