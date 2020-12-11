#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: build_mechanism.sh
#
# !DESCRIPTION:  Runs KPP while preserving the custom het-rate
#  calculation and KPP integrator
#
# !CALLING SEQUENCE:
#  ./build_mechanism.sh fullchem   # builds mechanism in KPP/fullchem folder
#  ./build_mechanism.sh custom     # builds mechanism in KPP/custom folder
#
# !AUTHOR:
#  Melissa Sulprizio (mpayer@seas.harvard.edu)
#
# !REVISION HISTORY:
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

# Check that directory exists
if [ "x${1}" != "x" ]; then
    mechanismDir=$1
    if [ ! -d "$mechanismDir" ]; then
        echo "ERROR: Mechanism directory '$1' does not exist."
        exit 1
    fi
else
    echo "ERROR: You must pass the mechanism directory name as argument."
    exit 1
fi

# Remove prior gckpp*.F90 files
# (but preserve a few files that do not need to be rebuilt)
cd $mechanismDir
mv gckpp_HetRates.F90 HETCODE
mv gckpp_Integrator.F90 INTEGRATOR
mv gckpp_Precision.F90 PRECISION
rm -f *.F90
rm -f *.o

# Build the mechanism and change extension from *.f90 to *.F90
# Halt if the gckpp_Rates file was not built
kpp gckpp.kpp
for a in $(ls *.f90); do mv -v $a ${a%.f90}.F90; done
if [[ ! -e gckpp_Rates.F90 ]]; then
  echo "KPP failed to build gckpp_Rates.F90! Aborting."
  exit 1
fi

# Restore the preserved files to their original names
mv HETCODE gckpp_HetRates.F90
mv INTEGRATOR gckpp_Integrator.F90
mv PRECISION gckpp_Precision.F90

# Insert code to enable reaction rate diagnostics in gckpp_Function.F90
line_new="SUBROUTINE Fun ( V, F, RCT, Vdot, Aout )"
line_orig="SUBROUTINE Fun ( V, F, RCT, Vdot )"
sed -i -e "s|${line_orig}|${line_new}|" gckpp_Function.F90

line1="  REAL(dp), optional :: Aout(NREACT)\n\n"
line2="! Computation of equation rates"
sed -i -e "s|${line2}|${line1}${line2}|" gckpp_Function.F90

line1="! Aout\n"
line2="  if(present(Aout)) Aout(:) = A(:)\n\n"
line3="! Aggregate function"
sed -i -e "s|${line3}|${line1}${line2}${line3}|" gckpp_Function.F90

# Run python parser OHreactParser.py. This will create fortran code
# for subroutine Get_OHreactivity and insert it into gckpp_Util.F90
python ../OHreact_parser.py

# Change back to the prior directory
cd ..

exit 0
#EOC
