# README for GEOS-Chem Classic Parallelization Tests

## Directories

Components of GEOS-Chem Classic parallelization tests have been separated into these directories:

`bin`

  - Contains GEOS-Chem Classic executable files

`build`

  - Directories for building GEOS-Chem Classic executables

`CodeDir`

  - Symbolic link to the `GCClassic` superproject directory

`env`

  - Contains an environment file that loads the required libraries needed to run GEOS-Chem Classic.

`logs`

  - Contains log files from the parallelization tests.

    - `results.compile.log`: Results of GEOS-Chem Classic compilation tests.
    - `results.parallel.log`: Results of GEOS-Chem Classic parallel tests.
    - `compile.*.log`: Output of individual compilation tests
    - `parallel.*.log`: Output of individual parallel tests
    - `lsf-*.txt`: LSF scheduler job logs
    - `slurm*.out`: SLURM scheduler job logs

`rundirs`

  - Contains individual GEOS-Chem Classic run directories.

`scripts`

  - Contains scripts that are used for running the GEOS-Chem Classic parallelization tests.  These are copied from the `test/GCClassic/parallel` folder of the GEOS-Chem "science codebase" repository. 

## For more information

Please see the `scripts/README.md` file for detailed instructions on running GEOS-Chem Classic parallelization tests.
