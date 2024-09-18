#! /bin/csh -f

# Set batch system options for SGI

#-------------------------------------------------------
# [1] Set necessary environment variables
#-------------------------------------------------------

# your HAVE TO change it with your directory root
setenv ROOTDIR /group_homes/lzu_public/home/u120220909911/Summer/Last/LE/RU-FY2

# 1) set clm include directory root
setenv CLM_INCDIR $ROOTDIR/include

# 2) set clm raw land data directory root
setenv CLM_RAWDIR /group_homes/lzu_public/home/u120220909911/CoLM/rawdata

# 3) set clm surface data rectory root
setenv CLM_SRFDIR $ROOTDIR/mksrfdata

# 4) set clm input data directory root
setenv CLM_DATADIR $ROOTDIR/data

# 5) set clm initial directory root
setenv CLM_INIDIR $ROOTDIR/mkinidata

# 6) set clm source directory root
setenv CLM_SRCDIR $ROOTDIR/main

# 7) set executable directory
setenv CLM_EXEDIR $ROOTDIR/run

# 8) set output directory
setenv CLM_OUTDIR $ROOTDIR/output
#mkdir -p $CLM_OUTDIR 2>/dev/null

# set debugging flag to TRUE or FALSE
setenv DEBUG  FALSE     

# set NTASKS to number of MPI tasks (if 1 than SPMD is off)   
# set NTHRDS to number of OPENMP threads

setenv ARCH `uname -s`

if ($ARCH == 'IRIX64') then
  # the following are default values that the script provides
  setenv NTASKS  1       
  setenv NTHRDS  28
endif

#------------------------------------------------------
# build define.h in ./include directory
#------------------------------------------------------

#\cat >! .tmp << EOF; cmp -s .tmp $CLM_INCDIR/define.h || mv -f .tmp $CLM_INCDIR/define.h
#undef coup_atmosmodel
#undef RDGRID
#define offline
#define USGS24
#define EcoDynamics
#define LANDONLY
#define WR_HOURLY
#EOF


#------------------------------------------------------
# [2] compling and executing clm surface data making
#-------------------------------------------------------

# Compile
cd $CLM_SRFDIR
# Executing clm initialization'
$CLM_SRFDIR/srf.x < $CLM_EXEDIR/srfdat.stdin 
#echo 'CLM Making Surface Data comleted'



#-------------------------------------------------------
# [3] compling and executing clm initialization
#-------------------------------------------------------

# Compile
cd $CLM_INIDIR

# Executing clm initialization'
$CLM_INIDIR/initial.x < $CLM_EXEDIR/inidat.stdin 

#echo 'CLM Initialization Completed'


#-------------------------------------------------------
# [4] compiling and executing clm model
#-------------------------------------------------------

# Compile
cd $CLM_SRCDIR
# Executing clm'
$CLM_SRCDIR/clm.x < $CLM_EXEDIR/timeloop.stdin 


#echo '-----------------------------------------------------------------'
#echo ' End of nqs shell script                                         '
# qsub -l ncpus=1 -l cput=00:10:00 jobclm.csh
#echo '-----------------------------------------------------------------'
