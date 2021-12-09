#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: build_mechanism.sh
#
# !DESCRIPTION: Runs KPP to create new chemical mechanism solver files
#  (in Fortran-90 format) while preserving the heterogeneous chemistry file.
#
# !CALLING SEQUENCE:
#  ./build_mechanism.sh fullchem   # builds mechanism in KPP/fullchem folder
#  ./build_mechanism.sh Hg         # builds mechanism in KPP/Hg folder
#  ./build_mechanism.sh custom     # builds mechanism in KPP/custom folder
#
# !AUTHOR:
#  Melissa Sulprizio (mpayer@seas.harvard.edu) -- Initial version
#  Bob Yantosca (yantosca@seas.harvard.edu) -- Updates for KPP 2.3.0_gc+
#
# !REMARKS:
#  (1) Requires KPP version 2.3.0_gc or later.
#
#  (2) KPP may have issues parsing the RHS of an equation when it needs
#      to be split up into more than one F90 line.  The quick solution
#      is to keep the length of the RHS at about ~100 characters or less.
#      https://github.com/geoschem/KPP/issues/1
#
# !REVISION HISTORY:
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

#============================================================================
# Check KPP version
#============================================================================
kppCodeDir=$(basename $KPP_HOME)
if [[ "x${kppCodeDir}" != "xkpp-code" ]]; then
    echo "ERROR: You must use KPP version 2.3.0_gc or later!"
    echo 
    exit 1
fi

#============================================================================
# Check that directory exists before proceeding
#============================================================================
if [ "x${1}" != "x" ]; then
    mechDir=${1}
    if [ ! -d "$mechDir" ]; then
        echo "ERROR: Mechanism directory '${mechDir}' does not exist."
        exit 1
    fi
else
    echo "ERROR: You must pass the mechanism directory name as argument."
    exit 1
fi

#============================================================================
# Remove prior files that have been built with KPP
# while leaving those files containing the chemistry mechanism specification
#============================================================================
cd ${mechDir}

# Find the mechanism name (which is at the top of the .eqn file)
mechName=$(grep ".eqn" *.eqn)
mechName="${mechName/\{ /}"
mechName="${mechName/\.eqn/}"

# Exit if mechanism name isn't found
if [[ "x${mechName}" == "x" ]]; then
    echo "Could not find the mechanism name in the ${mechanismDir}/*.eqn file"
    exit 1
fi

# Prefix for files
mechPrefix="gckpp"

# Remove these files, which will be will be regnerated by KPP
filesToRemove=(                     \
    ${mechPrefix}.map               \
    ${mechPrefix}_Function.F90      \
    ${mechPrefix}_Global.F90        \
    ${mechPrefix}_Initialize.F90    \
    ${mechPrefix}_Integrator.F90    \
    ${mechPrefix}_Jacobian.F90      \
    ${mechPrefix}_JacobianSP.F90    \
    ${mechPrefix}_LinearAlgebra.F90 \
    ${mechPrefix}_Model.F90         \
    ${mechPrefix}_Monitor.F90       \
    ${mechPrefix}_Parameters.F90    \
    ${mechPrefix}_Precision.F90     \
    ${mechPrefix}_Rates.F90         \
    ${mechPrefix}_Util.F90          \
)
for f in ${filesToRemove[@]}; do 
    rm -f $f
done

# Also remove any generated files
rm -f *.o *.mod *.a

#============================================================================
# Build the mechanism!
#============================================================================
if [[ -f ${mechPrefix}.kpp ]]; then
    kpp ${mechPrefix}.kpp
else
    echo "Could not find the ${mechPrefix}.kpp file... Aborting!"
    exit 1
fi

# Remove the GNU Makefile (not needed, since we use CMake)
[[ -f Makefile_${mechPrefix} ]] && rm -f Makefile_${mechPrefix}

# If the gckpp_Rates.F90 file is not found, there was an error
if [[ ! -f ${mechPrefix}_Rates.F90 ]]; then
  echo "KPP failed to build ${mechPrefix}_Rates.F90! Aborting."
  exit 1
fi

#============================================================================
# Strip unwanted characters in gckpp_Rates.F90
# These seem to be created by KPP due to issues in breaking long lines
# We might be able to get rid of this later on with thenew mechanism
#
# NOTE: This should be unnecessary if we apply the fix described at:
# https://github.com/geoschem/KPP/issues/1
#============================================================================
line1="         write(6,'(a)') 'GCJPLEQ: Missing parameters for P-dependent reaction.'I2O3"
line2="         write(6,'(a)') 'GCJPLEQ: Missing parameters for P-dependent reaction.'"
sed -i -e "s|${line1}|${line2}|" ${mechPrefix}_Rates.F90

#============================================================================
# Run python parser OHreactParser.py. This will create fortran code
# for subroutine Get_OHreactivity and insert it into gckpp_Util.F90
#
# TODO: Port this to C and include within KPP
# Also note: Don't do this for offline mechanisms such as Hg
#============================================================================
[[ "x${mechName}" != "xHg" ]] && python ../OHreact_parser.py

#============================================================================
# Change back to the prior directory and exit
#============================================================================
cd ..
exit 0
#EOC
