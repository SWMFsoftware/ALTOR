module PIC_ModLaserBeam
  use PIC_ModSize, ONLY: nDim, x_, y_, z_
  use PIC_ModMain, ONLY: tSimulation
  use ModNumConst, ONLY: cOne,cHalf,cDegToRad
  use PIC_ModProc,      ONLY:iProc
  implicit none
  SAVE
  PRIVATE
  real :: laserAmplitude(2:3)=0.0 !amplitude * polarization for Ey,Ez          
  real :: phaseShift=0.0 !Ey(center)=cos(0), Ez(center)=sin(0)                 
  real :: pulseWidth(nDim)=-1. !In cycles, wavelengths                         
  real :: xFocus(nDim)=-1. !if(pulseWidth(2)<=0.0.or. pulseWidth(3)<=0) then   
                           !           pulseWidth(2:3)=0.42*xFocus1(1)        
  !                                                                           
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
    real :: xx(nDim)
    real ::pX(nDim)
    real :: focusDist,timeReal,tmp
    real :: laser_profile,laser_phase
    !==========================================
    !
    !Works for the left boundary only
    !
    xx(:)=(/x,y,z/)
    !
    if(tSimulation>=timePulseBegin .and. tSimulation <= timePulseEnd) then
       if(x<0.0) then 
          focusDist=sqrt(sum((xFocus-xx)**2))
          timeReal=tSimulation+focusDist-xFocus(1)
          pX(:)=0.0
          !
          if(timeReal<=timePulseEnd) then
             !
             laser_phase=focusDist-xFocus(1)+(xPulseCenter+tSimulation)
             !                                                          
             tmp=(-laser_phase *2.0/pulseWidth(1))**2
             if(tmp < 30.) pX(1)=exp(-tmp)
             tmp=(focusDist*asin((xx(2)-xFocus(2))/focusDist)&
                  *2.0/pulseWidth(2))**2
             if(tmp < 30.) pX(2)=exp(-tmp)
             tmp=(focusDist*asin((xx(3)-xFocus(3))/focusDist)&
                  *2.0/pulseWidth(3))**2
             if(tmp < 30.) pX(3)=exp(-tmp)
             !
             laser_profile=product(pX)*laserAmplitude(iDir)
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
    !                                                                          
    call read_var('laserAmplitude2',laserAmplitude(2))
    call read_var('laserAmplitude3',laserAmplitude(3))
    call read_var('phaseShift',phaseShift)
    call read_var('pulseWidth1',pulseWidth(1))
    call read_var('pulseWidth2',pulseWidth(2))
    call read_var('pulseWidth3',pulseWidth(3))
    call read_var('xFocus1',xFocus(1))
    call read_var('xFocus2',xFocus(2))
    call read_var('xFocus3',xFocus(3))
    !
    laser_pulse=2
    !                                                                         
    timePulseBegin=0.0
    xPulseCenter=-pulseWidth(1)
    timePulseEnd=pulseWidth(1)*2.0
    !                                                                        
    coeff=1.0/sqrt(2.0*alog(2.0)) !coeff=0.849322                           
    !                                                                        
    if(iProc==0)then
       write(*,'(a,es13.5)')'Laser pulse of duration [cycles]=',pulseWidth(1)
!       if(xFocus(1) <=0.0 ) then                                           
       if(pulseWidth(2)<=0.0.or.pulseWidth(3)<=0) then
          write(*,'(a,es13.5)')'is focused to a wavelength at x=',xFocus(1)
          else
          write(*,'(a,2es13.5)')'and of width [wavelengths]=',pulseWidth(2:3)
       end if
    end if
    !                                                                         
    pulseWidth=pulseWidth*coeff
    if(pulseWidth(2)<=0.0.or.pulseWidth(3)<=0) then
       pulseWidth(2:3)=xFocus(1)*sin(45.0*cDegToRad)
       laser_pulse=1
    end if
    !                                                                    
  end subroutine read_laser_beam
  !                  
  !
end module PIC_ModLaserBeam
