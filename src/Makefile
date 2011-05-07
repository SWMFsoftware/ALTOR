#^CFG COPYRIGHT UM
SHELL = /bin/sh

include ../Makefile.def

SEARCH_EXTRA =  

include ../Makefile.conf
include Makefile.DEPEND

#
#	Makefile for parallel code ALTOR
#


# BLAS and LAPACK libraries

OBJECTS= \
	PIC_ModSize.o \
	PIC_ModProc.o \
	PIC_ModMain.o \
	PIC_ModFormFactor.o \
	PIC_ModField_${nDim}d.o \
	PIC_ModParticleInField_${nDim}d.o \
	PIC_ModPassField_${nDim}d.o \
	PIC_ModParticles.o \
	PIC_ModRandom.o \
	PIC_set_param.o \
	PIC_ModThermal.o 


#2dperiodicMPI.o

ANALYZE: Analyzer2D.o
	${FTN} -o Analyze.exe Analyzer2D.o	

############################################################################

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCH} \
		${OBJECTS}


MY_LIB  = libALTOR.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

# Object files not included into the library
OBJECTS_EXE = PIC_main.o

# Other requireed libraries
LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a

# Libraries should be compiled first, because modules are used in main.
${OBJECTS_EXE}: ${LIBSHARE} ${LIBTIMING} ${MY_LIB}

EXE = ${BINDIR}/ALTOR.exe

ALTOR:
	make ${EXE}
	@echo ' '
	@echo Program ALTOR has been brought up to date.
	@echo ' '	

${EXE}: ${OBJECTS_EXE}
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ../${MY_LIB}; \
		ar -x ${LIBTIMING}; \
		ar -x ${LIBSHARE}
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} Tmp_/*.o \
		${Lflag1}
	rm -rf Tmp_ 

distclean: clean
	rm -f PIC_ModSize.f90




