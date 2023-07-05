# README for GCHP Integration Tests

## Directories

Components of GCHP integration tests have been separated into these directories:

`bin`

  - Contains GCHP executable files.

`build`

  - Directories for building GCHP executables.

`CodeDir`

  - Symbolic link to the GCHP superproject directory.

`env`

  - Contains an environment file that loads the software libraries needed to run GCHP.

`logs`

  - Contains log files from the integration tests.

    - `results.compile.log`: Results of GCHP compilation tests.
    - `results.execute.log`: Results of GCHP Classic execution tests.
    - `compile.*.log`: Output of individual compilation tests
    - `execute.*.log`: Output of individual execution tests
    - `lsf-*.txt`: LSF scheduler job logs
    - `slurm*.out`: SLURM scheduler job logs

`rundirs`

  - Contains individual GCHP run directories.

`scripts`

  - Contains the scripts that are used for running the GCHP integration tests.  These are copied from the `test/GCHP/integration` folder of the GEOS-Chem "science codebase" repository.

## For more information

Please see the `scripts/README.md` file for detailed instructions on running GCHP integration tests.
