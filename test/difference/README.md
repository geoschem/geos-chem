# README for Difference Tests

## Contents

`diffTest.sh`
- Script that looks for differences in two different integration test or parallel test directories.  
- Checks both diagnostic files and restart files in each integration test run directory for identicality.
- TODO (as of June 2023): Add capability to check outputs from two different GEOS-Chem Classic parallel test directories.  As of now `difftest.sh` can only be used with integration tests.

## Performing a difference test

Run the `difftest.sh` script as follows:

```console
$ cd test/difference
$ ./diffTest.sh /path/to/<REF> /path/to/<DEV>
```
where `<REF>` and `<DEV>` indicate the names of the integration tests that are being compared.

The script will compare diagnostic and restart files within each integration test run directory.  If no differences are found, you will see ouptut such as:

```console
Checking gc_4x5_merra2_fullchem_<SIM>
   -> No differences in OutputDir
   -> No differences in Restarts
```
for each type of simulation `<SIM>` (e.g. `fullchem`, `fullchem_benchmark`, `CH4`, etc.) that was included in the `<REF>` and `<DEV>` integration tests.

If differences are found, you will see output such as:

```console
Checking gc_4x5_merra2_<SIM>
   -> 2 differences found in OutputDir
      * <REF>/rundirs/gc_4x5_merra2_<SIM>/OutputDir/GEOSChem.Metrics.20190701_0000z.nc4 
        <DEV>/rundirs/gc_4x5_merra2_<SIM>/OutputDir/GEOSChem.Metrics.20190701_0000z.nc4 
      * <REF>/rundirs/gc_4x5_merra2_<SIM>/OutputDir/GEOSChem.SpeciesConc.20190701_0000z.nc4 
        <DEV>/rundirs/gc_4x5_merra2_<SIM>/OutputDir/GEOSChem.SpeciesConc.20190701_0000z.nc4 
   -> 1 difference found in Restarts
      * <REF>/rundirs/gc_4x5_merra2_<SIM>/Restarts/GEOSChem.Restart.20190701_0100z.nc4 
        <DEV>/rundirs/gc_4x5_merra2_<SIM>/Restarts/GEOSChem.Restart.20190701_0100z.nc4
```
