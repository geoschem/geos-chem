#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: build.sh
# 
# !DESCRIPTION: Cleans or compiles GEOS-Chem High Performance (GCHP) 
#  source code. Accepts keyword argument indicating clean or compile option.
#  Optionally accepts '--debug' argument to build with debug flags on. 
#\\
#\\
# !REMARKS:
#  (1) This script is used within the GCHP run directory Makefile but can
#      also be used as a standalone option to clean and compile source code.
#  (2) Type './build.sh' or './build.sh help' within the run directory for 
#      options.
#
# !REVISION HISTORY: 
#  19 Oct 2016 - E. Lundgren - Initial version
#  See git history for version changes.
#EOP
#------------------------------------------------------------------------------
#BOC

# Check usage
if [ $# == 0 ] || [ $1 == "help" ]; then
  echo "Script name: build.sh"
  echo "Arguments: "
  echo "   none: prints options to screen"
  echo "   1st:  clean and/or compile option"
  echo "   2nd:  Optional '--debug' to turn on debug flags"
  echo "Clean options:"
  echo "      clean_all        - clean GEOS-Chem, ESMF, MAPL, and FVdycore"
  echo "      clean_mapl       - clean GC, MAPL, and FVdycore (skip ESMF)"
  echo "      clean_core       - clean GC only (skip ESMF, MAPL, and FVdycore)"
  echo "Compile options:"
  echo "      build            - general build command"
  echo "Example usage:"

  echo "   ./build.sh build         # Build without debug flags"
  echo "   ./build.sh build --debug # Build with debug flags"
  exit 0
fi

# Set run directory path
rundir=$PWD
gcdir=$(readlink -f CodeDir)
gchpdir=${gcdir}/GCHP

# Check source your environment file. This requires first setting the gchp.env
# symbolic link using script setEnvironment in the run directory. 
# Be sure gchp.env points to the same file for both compilation and 
# running. You can copy or adapt sample environment files located in 
# ./environmentFileSamples subdirectory.
gchp_env=$(readlink -f gchp.env)
if [ ! -f ${gchp_env} ] 
then
   echo "ERROR: gchp.env symbolic link is not set!"
   echo "Copy or adapt an environment file from the ./environmentFileSamples "
   echo "subdirectory prior to running. Then set the gchp.env "
   echo "symbolic link to point to it using ./setEnvironment."
   echo "Exiting."
   exit 1
fi
echo "WARNING: You are using environment settings in ${gchp_env}"
source ${gchp_env}

# Check environment vars
if [[ "x${MPI_ROOT}" == x ]]; then
   echo "ERROR: MPI_ROOT is not set! See environmentFileSamples subdir for examples."
elif [[ "x${ESMF_COMM}" == x ]]; then
   echo "ERROR: ESMF_COMM is not set! See environmentFileSamples subdir for examples."
   exit 1
elif [[ "x${ESMF_COMPILER}" == x ]]; then
   echo "ERROR: ESMF_COMPILER is not set! See environmentFileSamples subdir for examples."
   exit 1
elif [[ ! -e CodeDir ]]; then
  echo "ERROR: CodeDir symbolic link to source code directory does not exist!"
  echo "You may use the setCodeDir function for this, e.g."
  echo "   ./setCodeDir /path/to/your/code"
  exit 1
fi

# Go to the source code directory

# Flag for when to exit
done=0

###############################
###     Clean Options       ###
###############################

#-----------------------------------------------------------------------
#   clean_all
#-----------------------------------------------------------------------
if [[ $1 == "clean_all"      ]]; then
   cd ${gcdir}
   make HPC=yes realclean
   cd ${gchpdir}
   make the_nuclear_option
   cd ${rundir}
   done=1

#-----------------------------------------------------------------------
#   clean_mapl
#-----------------------------------------------------------------------
elif [[ $1 == "clean_mapl"      ]]; then
   cd ${gcdir}
   make HPC=yes realclean
   cd ${gchpdir}
   make wipeout_fvdycore
   make wipeout_mapl
   cd ${rundir}
   done=1

#-----------------------------------------------------------------------
#   clean_core
#-----------------------------------------------------------------------
elif [[ $1 == "clean_core" ]]; then
   cd ${gcdir} 
   make HPC=yes realclean
   cd ${rundir}
   done=1

fi

if [[ ${done} == "1" ]]; then
   unset rundir
   exit 0
fi

###############################
###     Compile Options     ###
###############################

#-----------------------------------------------------------------------
#   build (compile) - with or without debug flags
#-----------------------------------------------------------------------
if [[ $1 == "build" ]]; then
   cd ${gcdir}
   if [[ $2 == "--debug" ]]; then
      make -j${NUM_JOB_SLOTS} CHEM=standard EXTERNAL_GRID=y TRACEBACK=y    \
                              DEBUG=y       BOUNDS=y        FPEX=y         \
                              BPCH_DIAG=n   MODEL_GCHP=y    hpc
   else 
      make -j${NUM_JOB_SLOTS} CHEM=standard EXTERNAL_GRID=y TRACEBACK=y    \
                              BPCH_DIAG=n   MODEL_GCHP=y    hpc
   fi
fi

cd ${rundir}

if [[ -e ${gcdir}/bin/geos ]]; then
   echo '###################################'
   echo '###    GCHP executable exists!  ###'
   echo '###################################'
   cp CodeDir/bin/geos .
else
   echo '###################################################'
   echo '###   WARNING: GCHP executable does not exist   ###'
   echo '###################################################'
fi

unset rundir
exit 0
