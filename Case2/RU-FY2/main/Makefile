# Makefile for CLM initialization 

# name of executable

EXECUTABLE = clm.x

#compiler options
RM = rm
FF = ifort
CPP = ifort -EP

INCLUDE_DIR = -I. -I../.. -I../include
FOPTS = -i4 -r8 -free -CB -extend_source $(INCLUDE_DIR)

####################################################################

OBJS =  precision.o         \
	phycon_module.o     \
	CLM.o               \
	CLMDRIVER.o         \
	CLMMAIN.o           \
	GETMET.o            \
	LAKE.o              \
	SOCEAN.o            \
        THERMAL.o           \
	WATER.o             \
	albland.o           \
	albocean.o          \
	combo.o             \
	dewfraction.o       \
	eroot.o             \
	flxwrite.o          \
	groundfluxes.o      \
	groundtem.o         \
	hCapacity.o         \
	hConductivity.o     \
	lai_empirical.o     \
	leafinterception.o  \
	leaftemone.o        \
	leaftemtwo.o        \
        lpwrite.o           \
	meltf.o             \
	moninobuk.o         \
	netsolar.o          \
	newsnow.o           \
	orb_coszen.o        \
	qsadv.o             \
	rstFileMod.o        \
	snowage.o           \
	snowcompaction.o    \
	snowfraction.o      \
	snowlayerscombine.o \
	snowlayersdivide.o  \
	snowwater.o         \
	soilwater.o         \
	stomata.o           \
	subsurfacerunoff.o  \
	surfacerunoff.o     \
	ticktime.o          \
	tridia.o            \
	twostream.o         \
	vec2xy.o

####################################################################

${EXECUTABLE} : ${OBJS} 
	${FF} $(FOPTS) $(OBJS) -o $@

.SUFFIXES : .o .F90 .f90

.F90.o:
	${CPP} ${INCLUDE_DIR} $*.F90 >$*.i
	sed '/^ *$$/d' $*.i >$*.f90; rm $*.i
	${FF} -c ${FOPTS} $*.f90 

clean :
	- ${RM} -f *.i *.o *.f90 *.mod *.x core.*
	- ${RM} -rf rii_files
