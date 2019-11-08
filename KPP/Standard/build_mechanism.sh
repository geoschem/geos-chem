#!/bin/bash
# Runs KPP while preserving the custom het-rate calculation and KPP integrator

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
./rename_f90.sh
mv HETCODE gckpp_HetRates.F90
mv INTEGRATOR gckpp_Integrator.F90
mv PRECISION gckpp_Precision.F90

# Insert code to enable reaction rate diagnostics in gckpp_Function.F90
# Currently this diagnostic is implemented for use in GEOS only
line1="#if defined( MODEL_GEOS )\n"
line2="SUBROUTINE Fun ( V, F, RCT, Vdot, Aout )\n"
line3="#else\n"
line4="SUBROUTINE Fun ( V, F, RCT, Vdot )"
line5="\n#endif\n"
sed -i -e "s|${line4}|${line1}${line2}${line3}${line4}${line5}|" gckpp_Function.F90

line1="#if defined( MODEL_GEOS )\n"
line2="  REAL(dp), optional :: Aout(NREACT)\n"
line3="#endif\n\n"
line4="! Computation of equation rates"
sed -i -e "s|${line4}|${line1}${line2}${line3}${line4}|" gckpp_Function.F90

line1="#if defined( MODEL_GEOS )\n"
line2="! Aout\n"
line3="  if(present(Aout)) Aout(:) = A(:)\n"
line4="#endif\n\n"
line5="! Aggregate function"
sed -i -e "s|${line5}|${line1}${line2}${line3}${line4}${line5}|" gckpp_Function.F90

# Run python parser OHreactParser.py. This will create fortran code
# for subroutine Get_OHreactivity and insert it into gckpp_Util.F90
python OHreact_parser.py

exit 0
