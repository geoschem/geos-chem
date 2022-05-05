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
# Check that directory exists before proceeding
#============================================================================
if [ "x${1}" != "x" ]; then
    mechDir=${1}
    if [ ! -d "$mechDir" ]; then
        echo "ERROR: Mechanism directory ${mechDir} does not exist."
        exit 1
    fi
else
    echo "ERROR: You must pass the mechanism directory name as argument."
    exit 1
fi

#============================================================================
# Remove prior files that have been built with KPP
# while leaving those files containing the chemistry mechanism specification
#
# NOTE: KPP-generated source code files for GEOS-Chem will have the
# prefix "gckpp".  This is necessary because a consistent naming scheme
# needs to used so that modules in GeosCore etc. can find KPP code.
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

# Remove these files, which will be will be regnerated by KPP
filesToRemove=(                     \
    gckpp.map               \
    gckpp_Function.F90      \
    gckpp_Global.F90        \
    gckpp_Initialize.F90    \
    gckpp_Integrator.F90    \
    gckpp_Jacobian.F90      \
    gckpp_JacobianSP.F90    \
    gckpp_LinearAlgebra.F90 \
    gckpp_Model.F90         \
    gckpp_Monitor.F90       \
    gckpp_Parameters.F90    \
    gckpp_Precision.F90     \
    gckpp_Rates.F90         \
    gckpp_Util.F90          \
)
for f in ${filesToRemove[@]}; do 
    rm -f $f
done

# Also remove any generated files
rm -f *.o *.mod *.a

#============================================================================
# Build the mechanism!
#============================================================================
if [[ -f gckpp.kpp ]]; then
    kpp gckpp.kpp
else
    echo "Could not find the 'gckpp.kpp' file... Aborting!"
    exit 1
fi

# Remove the GNU Makefile (not needed, since we use CMake)
[[ -f Makefile_gckpp ]] && rm -f Makefile_gckpp

# If the gckpp_Rates.F90 file is not found, there was an error
if [[ ! -f gckpp_Rates.F90 ]]; then
  echo "KPP failed to build 'gckpp_Rates.F90'! Aborting."
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
sed -i -e "s|${line1}|${line2}|" gckpp_Rates.F90

#============================================================================
# Run python parser OHreactParser.py. This will create fortran code
# for subroutine Get_OHreactivity and insert it into gckpp_Util.F90
#
# TODO: Port this to C and include within KPP
#============================================================================
python ../OHreact_parser.py ${mechName}

#============================================================================
# Change back to the prior directory and exit
#============================================================================
cd ..
exit 0
#EOC
