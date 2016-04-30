#
# List the default target first for stand alone mode
#
DEFAULT_TARGET = ALTOR
DEFAULT_EXE    = ${DEFAULT_TARGET}.exe

default : ${DEFAULT_TARGET}

include Makefile.def

nDim = 3

FF = 1.5

#
# Menu of make options
#
help:
	@echo ' '
	@echo '  You can "make" the following:'
	@echo ' '
	@echo '    <default> ${DEFAULT_TARGET} in stand alone mode, help in SWMF'
	@echo ' '
	@echo '    help         (makefile option list)'
	@echo '    install      (install ALTOR)'
	@echo '    2d           (re-install the code with nDim=2)'
	@echo '    3d           (re-install the code with nDim=3)'
	@echo ' '
	@echo '    LIB     (Component library libPC for SWMF)'
	@echo '    ALTOR (ALTernating-ORder interpolation scheme for PIC)'
	@echo '    nDim=2 compiles the code with nDim=2'
	@echo '    NOMPI   (NOMPI library for compilation without MPI)'
	@echo ' '
	@echo '    rundir      (create run directory for standalone or SWMF)'
	@echo '    rundir RUNDIR=run_test (create run directory run_test)'
	@echo ' '
	@echo "    nompirun    (make and run ${DEFAULT_EXE} on 1 PE)"
	@echo "    mpirun      (make and mpirun ${DEFAULT_EXE} on 8 PEs)"
	@echo "    mpirun NP=7 RUNDIR=run_test (run on 7 PEs in run_test)"
	@echo "    mprun NP=5  (make and mprun ${DEFAULT_EXE} on 5 PEs)"
	@echo ' '	
	@echo '    clean     (remove temp files like: *~ *.o *.kmo *.mod *.T *.lst core)'
	@echo '    distclean (equivalent to ./Config.pl -uninstall)'
	@echo '    dist      (create source distribution tar file)'

INSTALLFILES =	src/Makefile.DEPEND srcBATL/Makefile.DEPEND


install: src/PIC_ModSize.f90 srcBATL/Makefile srcBATL/BATL_size.f90
	touch ${INSTALLFILES}

srcBATL/Makefile:
	rm -rf srcBATL; mkdir srcBATL 
	cd srcBATL_orig; cp BATL*.f90 Makefile* ../srcBATL; \
	cd ../srcBATL; ${SCRIPTDIR}/Methods.pl PC *.f90; \
	${SCRIPTDIR}/Rename.pl -w -r -common=PC *.f90; \
	rm -f *~

src/PIC_ModSize.f90: src/PIC_ModSize_${nDim}d.f90
	cp -f src/PIC_ModSize_${nDim}d.f90 src/PIC_ModSize.f90

srcBATL/BATL_size.f90: srcBATL/BATL_size_orig.f90
	cp -f srcBATL/BATL_size_orig.f90 srcBATL/BATL_size.f90

LIB:
	cd srcBATL; make LIB
	cd src; make LIB nDim=${nDim} FF=${FF}

ALTOR:
	cd ${SHAREDIR}; ${MAKE} LIB
	cd ${TIMINGDIR}; ${MAKE} LIB
	cd srcBATL; make LIB
	cd src; ${MAKE} LIB nDim=${nDim} FF=${FF}
	cd src; ${MAKE} ALTOR nDim=${nDim} FF=${FF}

QED:
	cd ${SHAREDIR}; ${MAKE} LIB
	cd ${TIMINGDIR}; ${MAKE} LIB
	cd srcQED; ${MAKE} LIB nDim=${nDim} FF=${FF}
	cd srcQED; ${MAKE} ALTOR

NOMPI:
	cd util/NOMPI/src; make LIB

2d:
	touch src/PIC_ModSize_2d.f90
	cd src;rm -f *3d.o
	make install nDim=2

3d:
	touch src/PIC_ModSize_3d.f90
	cd src;rm -f *2d.o
	make install

# The MACHINE variable holds the machine name for which scripts should
# be copied to the run directory when it is created.  This is used mostly
# when several different machines have the same operating system,
# but they require different batch queue scripts.
# If MACHINE is empty or not defined, all scripts for the current OS will
# be copied.
#
# The default is the short name of the current machine
MACHINE = `hostname | sed -e 's/\..*//;s/[0-9]*$$//'`
COMPONENT = PC

rundir:
	mkdir -p ${RUNDIR}/${COMPONENT}
	cd ${RUNDIR}/${COMPONENT}; \
		ln -s ${PCDIR}/Param .
	@(if [ "$(STANDALONE)" != "NO" ]; then \
		touch ${DIR}/share/JobScripts/job._TMP_${MACHINE}; \
		touch ${DIR}/share/JobScripts/_TMP_.${MACHINE}.pl; \
		cp ${DIR}/share/JobScripts/job.*${MACHINE}* ${RUNDIR}/; \
		cp ${DIR}/share/JobScripts/*.${MACHINE}.pl ${RUNDIR}/; \
		rm -f ${RUNDIR}/*_TMP_* ${DIR}/share/JobScripts/*_TMP_*; \
		cp -f Param/PARAM.DEFAULT ${RUNDIR}/PARAM.in; \
		touch ${RUNDIR}/core; chmod 444 ${RUNDIR}/core; \
		cd ${RUNDIR}; ln -s ${BINDIR}/${DEFAULT_EXE} .; \
		ln -s ${COMPONENT}/* .; \
	fi);
 
#
#       Run the default code on NP processors
#

NP=8

mpirun: ${DEFAULT_TARGET}
	cd ${RUNDIR}; mpirun -np ${NP} ./${DEFAULT_EXE}

mprun: ${DEFAULT_TARGET}
	cd ${RUNDIR}; mprun -np ${NP} ./${DEFAULT_EXE}

nompirun: ${DEFAULT_TARGET}
	cd ${RUNDIR}; ./${DEFAULT_EXE}


#
# Cleaning
#

clean:
	@touch ${INSTALLFILES}
	cd src; make clean
	@(if [ -d util  ]; then cd util;  make clean; fi);
	@(if [ -d share ]; then cd share; make clean; fi);

distclean: 
	./Config.pl -uninstall
	rm -rf srcBATL

allclean:
	@touch ${INSTALLFILES}
	cd src; make distclean

TESTDIR = run_test

test:
	-@(${MAKE} test_altor)

test_altor:
	@echo "test_altor_compile..." > test_altor.diff
	${MAKE} test_altor_compile
	@echo "test_altor_rundir..." >> test_altor.diff
	${MAKE} test_altor_rundir
	@echo "test_altor_run..." >> test_altor.diff
	${MAKE} test_altor_run
	@echo "test_altor_check..." >> test_altor.diff
	${MAKE} test_altor_check

test_altor_compile:
	./Config.pl -g=64,64,64,2,1
	${MAKE} 

test_altor_rundir: 
	rm -rf ${TESTDIR}
	${MAKE} rundir RUNDIR=${TESTDIR} STANDALONE=YES PCDIR=`pwd`
	cd ${TESTDIR}; cp -f Param/PARAM.TEST PARAM.in

test_altor_run:
	cd ${TESTDIR}; ${MPIRUN} ./ALTOR.exe > runlog

test_altor_check:
	${SCRIPTDIR}/DiffNum.pl -t -r=1e-5 -a=3e-8 \
	Param/TestOutput/log_noise.out \
	${TESTDIR}/log_n0001.out > test_altor.diff
	ls -l test_altor.diff
