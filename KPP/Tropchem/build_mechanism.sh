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

exit 0