# README for GEOS-Chem Classic Integration Tests

## Overview:

This directory contains:

1. Scripts to submit automated GEOS-Chem Classic integration tests, which will test the following processes:

    - GEOS-Chem Classic run directory creation
    - GEOS-Chem Classic code configuration and compilation
    - GEOS-Chem Classic execution

Integration tests are short GEOS-Chem Classic simulations.  Executables are compiled with all debugging options.  This should reveal any coding errors or run-directory configuration errors.

## Files

### Integration Test Scripts

- `integrationTest.sh`
  - Driver script to run GEOS-Chem Classic integration tests
- `integrationTestCreate.sh`
  - Script to create GEOS-Chem Classic run directories for integration tests
- `integrationTestCompile.sh`
  - Script to compile GEOS-Chem Classic executables
- `integrationTestExecute.sh`
   - Script to run GEOS-Chem Classic integration test simulations

### Shared Scripts

- `commonFunctionsForTests.sh`
  - Link to `../../shared/commonFunctionsForTests.sh`, which contains global variables and functions for the integration and parallelization tests.

## Running GEOS-Chem Classic Integration and/or Parallelization Tests

Before you submit any GEOS-Chem Classic integration tests, please take a moment to:

1. Verify that the `GCClassic` superproject is checked out to the correct branch and commit.
2. Run `git submodule update --init --recursive` in order to update all submodules.
3. Verify that the `HEAD` commit of the `GEOS-Chem` submodule contains the code that you wish to test. (If not, then check out the proper branch.)
4. Verify that the `HEAD` commit of the `HEMCO` submodule contains the code that you wish to test. (If not, then check out the proper branch.)

### With the SLURM scheduler

```console
$ cd /path/to/GCClassic     # Path to GCClassic superproject directory
$ cd test/integration/GCClassic
$ ./integrationTest.sh -d /path/to/test/dir -e /path/to/env-file -s slurm -p partition
```

### With the LSF scheduler

```console
$ cd /path/to/GCClassic     # Path to GCClassic superproject directory
$ cd test/integration/GCClassic
$ ./integrationTest.sh -d /path/to/test/dir -e /path/to/env-file -s lsf -p partition
```

### Interactively at the command line

```console
$ cd /path/to/GCClassic     # Path to GCClassic superproject directory
$ cd test/integration/GCClassic
$ ./integrationTest.sh -d /path/to/test/dir -e /path/to/env-file
```

## Command-line Arguments

The integration test scripts accept the following command-line arguments

`-d /path/to/test/dir` specifies the root directory where integration test subdirectories and scripts will be placed.

`-e /path/to/env-file` specifies the file that is used to initialize the software environment.  This file will typically contain either `module load` or `spack load` commands.

`-s scheduler` specifies the choice of scheduler (case-insensitive). SLURM and LSF schedulers are currently supported.  You may omit this for running interactively.

`-p partition` specifies the partition for the scheduler.  You may omit this for running interactively.

`-s` specifies that the SLURM scheduler will be used to run the integration test scripts.

`-n` specifies that missing species in restart files will not be bootstrapped to a missing value.  This can be used to test if "out-of-the-box" simulations will fail before a version release.

`-h` displays a help screeen.


You can also use long names for the option switches:
- `--directory` instead of `-d`
- `--env-file` instead of `-e`
- `--help` instead of `h`
- `--no-bootstrap` instead of `n`
- `--partition partition` instead of `-p`
- `--scheduler scheduler` instead of `-s`

There is an additional option (`-q` or `--quick`) that will run only a couple of integration tests instead of the full suite.  This is intended for development and debugging.  You will normally not need to use this option.
