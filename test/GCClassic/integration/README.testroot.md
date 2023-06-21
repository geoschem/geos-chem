# README for GEOS-Chem Classic Integration Tests

## Directories

Components of GEOS-Chem Classic integration tests have been separated into these directories:

`bin`

  - Contains GEOS-Chem Classic executable files.

`build`

  - Directories for building GEOS-Chem Classic executables.

`CodeDir`

  - Symbolic link to the `GCClassic` superproject directory.

`env`

  - Contains an environment file that loads the software libraries needed to run GEOS-Chem Classic.

`logs`

  - Contains log files from the integration tests.

    - `results.compile.log`: Results of GEOS-Chem Classic compilation tests.
    - `results.execute.log`: Results of GEOS-Chem Classic execution tests.
    - `compile.*.log`: Output of individual compilation tests
    - `execute.*.log`: Output of individual execution tests
    - `lsf-*.txt`: LSF scheduler job logs
    - `slurm*.out`: SLURM scheduler job logs

`rundirs`

  - Contains individual GEOS-Chem Classic run directories.

`scripts`

  - Contains the scripts that are used for running the GEOS-Chem Classic integration tests.  These are copied from the `test/GCClassic/integration` folder of the GEOS-Chem "science codebase" repository.

## For more information

Please see the `scripts/README.md` file for detailed instructions on running GEOS-Chem Classic integration tests.
