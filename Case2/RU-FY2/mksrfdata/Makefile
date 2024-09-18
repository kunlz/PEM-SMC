# Makefile for CLM surface data processing

# name of executable

EXECUTABLE = srf.x

#compiler options
RM = rm
FF = ifort
CPP = ifort -EP

INCLUDE_DIR = -I. -I../.. -I../include
FOPTS = -i4 -r8 -free -CB -extend_source $(INCLUDE_DIR)

####################################################################
.DEFAULT :
	ln -s ../main/precision.F90 .

OBJS =  precision.o  \
	mksrfdata.o  \
	rdgrid.o     \
	crgrid.o     \
	rdlanddata.o \
	celledge.o   \
	cellarea.o

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
