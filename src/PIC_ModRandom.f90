!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!The random number generator.
!From the Buneman code.
module PIC_ModRandom
  implicit none
  private !Except
  integer::iSeed=0
  !Methods
  public::rand
  public::init_rand
contains
  subroutine init_rand(iSeedIn)
    integer,optional,intent(in)::iSeedIn
    if(present(iSeedIn))then
       iSeed=iSeedIn
    else
       iSeed=1
    end if
  end subroutine init_rand
  !-----------------------------------------------------------------!
  !===============================================================C
  REAL FUNCTION RAND()
    ISEED=ISEED*48828125
    IF(ISEED < 0) ISEED=(ISEED+2147483647)+1
    if(iSeed==0) iSeed=1
    RAND=FLOAT(ISEED)/2147483647
  END FUNCTION RAND
end module PIC_ModRandom
