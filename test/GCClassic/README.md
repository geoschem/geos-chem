# README for Integration Tests and Parallelization Tests

## Overview:

This directory contains:

1. Scripts to submit automated GEOS-Chem Classic integration tests, which will test the following processes:

    - GEOS-Chem Classic run directory creation
    - GEOS-Chem Classic code configuration and compilation
    - GEOS-Chem Classic execution

2. Scripts to submit automated GEOS-Chem Classic parallelization tests, which check for improperly-parallelized loops.

Integration tests are short GEOS-Chem Classic simulations.  Executables are compiled with all debugging options.  This should reveal any coding errors or run-directory configuration errors.

Parallelization tests are similar to integration tests, except that each simulation is run twice with a different number of cores. The restart files from each run then are compared for bitwise-identicality.

## Files

### Integration Test Scripts

- `intTest.sh`
  - Driver script to run GEOS-Chem Classic integration tests
- `intTestCreate.sh`
  - Script to create GEOS-Chem Classic run directories for integration tests
- `intTestCompile.sh`
  - Script to compile GEOS-Chem Classic executables
- `intTestExecute.sh`
   - Script to run GEOS-Chem Classic integration test simulations

### Parallelization Test Scripts

- `parTest.sh`
  - Driver script to run GEOS-Chem Classic parallelization tests
- `parTestCreate.sh`
  - Script to create GEOS-Chem Classic run directories for parallelization tests
- `parTestCompile.sh`
  - Script to compile GEOS-Chem Classic executables
- `parTestExecute.sh`
   - Script to run GEOS-Chem Classic parallelization test simulations

### Shared Scripts

- `commonFunctionsForTests.sh`
  - Link to `../shared/commonFunctionsForTests.sh`, which contains global variables and functions for the integration and parallelization tests.

## Running GEOS-Chem Classic Integration and/or Parallelization Tests

Before you submit any GEOS-Chem Classic integration or parallelization tests, take a moment and make sure that:

1. The `GCClassic` superproject is checked out to the correct branch and commit.
2. The `GEOS-Chem` submodule is checked out to the correct branch and commit.
3. The `HEMCO` submodule is checked out to the correct branch and commit.

### With the SLURM scheduler

```console
$ cd test/GCClassic
$ ./intTest.sh -d /path/to/test/dir -e /path/to/env-file -s -p partition   # Integration tests
$ ./parTest.sh -d /path/to/test/dir -e /path/to/env-file -s -p partition   # Parallelization tests
```

### With the LSF scheduler

```console
$ cd test/GCClassic
$ ./intTest.sh -d /path/to/test/dir -e /path/to/env-file -l -p partition   # Integration tests
$ ./parTest.sh -d /path/to/test/dir -e /path/to/env-file -l -p partition   # Parallelization tests
```

### Interactively at the command line

```console
$ cd test/GCClassic
$ ./intTest.sh -d /path/to/test/dir -e /path/to/env-file   # Integration tests
$ ./parTest.sh -d /path/to/test/dir -e /path/to/env-file   # Parallelization tests
```

## Command-line Arguments

The integration test scripts and parallelization test scripts accept the following command-line arguments

`-d /path/to/test/dir` specifies the root directory where integration/parallelization test subdirectories and scripts will be placed.

`-e /path/to/env-file` specifies the file that is used to initialize the software environment.  This file will typically contain either `module load` or `spack load` commands.

`-l` specifies that the LSF scheduler will be used to run the integration/parallelization test scripts.

`-p partition` specifies the partition for the SLURM or LSF scheduler.  You may omit this for running interactively.

`-s` specifies that the SLURM scheduler will be used to run the integration/parallelization test scripts.

You can also use long names for the option switches:
- `--directory` instead of `-d`
- `--env-file` instead of `-e`
- `--lsf` instead of `-l`
- `--partition` instead of `-p`
- `--slurm` instead of `-s`

There is an additional option (`-q` or `--quick`) that will run only a couple of integration tests instead of the full suite.  This is intendend for development and debugging.  You will normally not need to use this option.
