#!/bin/bash
# Runs KPP while preserving the custom het-rate calculation and KPP integrator

if [ "$1" != "" ]; then
    # Check that directory exists
    mechanismDir=$1
    if [ ! -d "$mechanismDir" ]; then
        echo "ERROR: Mechanism directory '$1' does not exist."
        exit 1
    fi
else
    echo "ERROR: You must pass mechanism directory name as argument."
    exit 1
fi

cd $mechanismDir

mv gckpp_HetRates.F90 HETCODE
mv gckpp_Integrator.F90 INTEGRATOR
mv gckpp_Precision.F90 PRECISION
rm *F90
rm *o
kpp gckpp.kpp
if [[ ! -e gckpp_Rates.f90 ]]; then
  echo "KPP failure! Aborting."
  exit 1
fi
for a in $(ls *.f90); do mv -v $a ${a%.f90}.F90; done
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

cd ..

exit 0
