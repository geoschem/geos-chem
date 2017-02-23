#!/bin/bash
# Runs KPP while preserving the custom het-rate calculation and KPP integrator

function finish {
  if [[ -e HETCODE ]]; then
    mv HETCODE gckpp_HetRates.F90
  fi
  if [[ -e INTEGRATOR ]]; then
    mv INTEGRATOR gckpp_Integrator.F90
  fi
}
trap finish EXIT

mv gckpp_HetRates.F90 HETCODE
mv gckpp_Integrator.F90 INTEGRATOR
rm *F90
rm *o
kpp gckpp.kpp
if [[ ! -e gckpp_Rates.f90 ]]; then
  echo "KPP failure! Aborting."
  exit 1
fi
./rename_f90.sh
