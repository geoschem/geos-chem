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
#  ./build_mechanism.sh custom     # builds mechanism in KPP/custom folder
#
# !AUTHOR:
#  Melissa Sulprizio (mpayer@seas.harvard.edu) -- Initial version
#  Bob Yantosca (yantosca@seas.harvard.edu) -- Updates for KPP 2.3.0_gc
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
    mechanismDir=${1}
    if [ ! -d "$mechanismDir" ]; then
        echo "ERROR: Mechanism directory '${1}' does not exist."
        exit 1
    fi
else
    echo "ERROR: You must pass the mechanism directory name as argument."
    exit 1
fi

#============================================================================
# Remove prior files, except gckpp_HetRates.F90
#============================================================================
cd ${mechanismDir}
mv gckpp_HetRates.F90 HETCODE
rm -f *.F90 *.f90 *.o

#============================================================================
# Build the mechanism and change extension from *.f90 to *.F90
# Halt if the gckpp_Rates file was not built
#============================================================================
kpp gckpp.kpp
if [[ ! -e gckpp_Rates.F90 ]]; then
  echo "KPP failed to build gckpp_Rates.F90! Aborting."
  exit 1
fi

#============================================================================
# Restore the preserved files to their original names
#============================================================================
mv HETCODE gckpp_HetRates.F90

#============================================================================
# Strip unwanted characters in gckpp_Rates.F90
# These seem to be created by KPP due to issues in breaking long lines
# We might be able to get rid of this later on with thenew mechanism
#============================================================================
line1="         write(6,'(a)') 'GCJPLEQ: Missing parameters for P-dependent reaction.'I2O3"
line2="         write(6,'(a)') 'GCJPLEQ: Missing parameters for P-dependent reaction.'"
sed -i -e "s|${line1}|${line2}|" gckpp_Rates.F90

#============================================================================
# Run python parser OHreactParser.py. This will create fortran code
# for subroutine Get_OHreactivity and insert it into gckpp_Util.F90
#============================================================================
python ../OHreact_parser.py

#============================================================================
# Change back to the prior directory and exit
#============================================================================
cd ..
exit 0
#EOC
