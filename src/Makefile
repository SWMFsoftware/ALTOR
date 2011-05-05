#^CFG COPYRIGHT UM
SHELL = /bin/sh
OS    = Linux

#MF = Makefile.${OS}
MF = 

DEFAULT_TARGET=PIC.exe

default : ${DEFAULT_TARGET}
#
#	Makefile for parallel code BATSRUS
#
#       Linux version
#
FTN     = mpif90

#PRECISION = -r8
PRECISION =

MPILIB = 
##MPILIB = -L. -lNOMPI						#^CFG  IF NOT SIMPLE

MPIHEADERS = mpif90.h
MPIF90H = mpif90.h

# Set optimization / debug levels
#OPT = -O0 -C -gline #		Debugging
#OPT = -O1 -C -gline #		Development
OPT = -O3 #			Production

OPTSAFE = -O0

Cflag0  = -c ${OPTSAFE} ${PRECISION}
Cflag1  = -c ${OPT} ${PRECISION}

Lflag1  = ${PRECISION} ${MPILIB} ${OPT}
Lflag2  = ${PRECISION} ${OPT}

nDim = 3

# BLAS and LAPACK libraries

OBJECTS=ModMpi.o \
	ModNumConst.o \
	PIC_ModSize_${nDim}d.o \
	PIC_ModProc.o \
	PIC_ModMain.o \
	PIC_ModFormFactor.o \
	PIC_ModField_${nDim}d.o \
	PIC_ModParticleInField_${nDim}d.o \
	PIC_ModPassField_${nDim}d.o \
	PIC_ModParticles.o \
	PIC_ModRandom.o \
	PIC_ModThermal.o \
	PIC_main.o
#2dperiodicMPI.o

ANALYZE: Analyzer2D.o
	${FTN} -o Analyze.exe Analyzer2D.o	


PIC.exe: ${OBJECTS}
	${FTN} -o PIC.exe ${OBJECTS} 
clean:
	rm -f *.o *.mod *~ *.exe

.SUFFIXES:
.SUFFIXES: .f90 .f .o

.f90.o:
	${FTN} ${Cflag1} $<

.f.o:
	${FTN} ${Cflag1} $<

#
#	Specific Rules
#


