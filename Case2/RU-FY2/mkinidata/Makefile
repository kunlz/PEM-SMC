# Makefile for CLM initialization 

# name of executable

EXECUTABLE = initial.x

#compiler options
RM = rm
FF = ifort
CPP = ifort -EP

INCLUDE_DIR = -I. -I../.. -I../include
FOPTS = -i4 -r8 -free -CB -extend_source $(INCLUDE_DIR)

####################################################################

.DEFAULT :
	ln -s ../main/precision.F90 .
	ln -s ../main/phycon_module.F90 .
	ln -s ../main/albland.F90 .
	ln -s ../main/albocean.F90 .
	ln -s ../main/twostream.F90 .
	ln -s ../main/snowfraction.F90 .
	ln -s ../main/lai_empirical.F90 .
	ln -s ../main/orb_coszen.F90 .
	ln -s ../main/rstFileMod.F90 .

OBJS =  precision.o phycon_module.o \
	initialize.o CLMINI.o \
	iniTimeConst.o iniTimeVar.o \
	orb_coszen.o lai_empirical.o rstFileMod.o \
	albland.o albocean.o twostream.o snowfraction.o

####################################################################

${EXECUTABLE} : ${OBJS} 
	${FF} ${FOPTS} ${OBJS} -o $@

.SUFFIXES : .o .F90 .f90

.F90.o:
	${CPP} ${INCLUDE_DIR} $*.F90 >$*.i
	sed '/^ *$$/d' $*.i >$*.f90; rm $*.i
	${FF} -c ${FOPTS} $*.f90 

clean :
	- ${RM} -f *.i *.o *.f90 *.mod *.x core.*
	- ${RM} -rf rii_files
