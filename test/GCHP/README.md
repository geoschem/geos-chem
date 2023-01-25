# README for Integration Tests and Parallelization Tests

## Overview:

This directory contains scripts to submit automated GCHP integration tests, which will test the following processes:

- GCHP run directory creation
- GCHP code configuration and compilation
- GCHP execution

Integration tests are short GCHP simulations.  Executables are compiled with all debugging options.  This should reveal any coding errors or run-directory configuration errors.

## Files

### Integration Test Scripts

- `intTest.sh`
  - Driver script to run GCHP integration tests
- `intTestCreate.sh`
  - Script to create GCHP run directories for integration tests
- `intTestCompile.sh`
  - Script to compile GCHP executables
- `intTestExecute.sh`
   - Script to run GCHP integration test simulations

### Shared Scripts

- `commonFunctionsForTests.sh`
  - Link to `../shared/commonFunctionsForTests.sh`, which contains global variables and functions for the integration tests.

## Running GCHP Integration Tests

Before you submit any GCHP integration tests, please take a moment to:

1. Verify that the `GCHP` superproject is checked out to the correct branch and commit.
2. Run `git submodule update --init --recursive` in order to update all submodules.
3. Verify that the `HEAD` commit of the `GEOS-Chem` submodule contains the code that you wish to test. (If not, then check out the proper branch.)
4. Verify that the `HEAD` commit of the `HEMCO` submodule contains the code that you wish to test. (If not, then check out the proper branch.)

### With the SLURM scheduler

```console
$ cd /path/to/gchp   # Path to the GCHP superproject directory
$ cd test
$ ./intTest.sh -d /path/to/test/dir -e /path/to/env-file -s -p partition
```

### With the LSF scheduler

```console
$ cd /path/to/gchp   # Path to the GCHP superproject directory
$ cd test
$ ./intTest.sh -d /path/to/test/dir -e /path/to/env-file -l -p partition
```

### Interactively at the command line

```console
$ cd /path/to/gchp   # Path to the GCHP superproject directory
$ cd test
$ ./intTest.sh -d /path/to/test/dir -e /path/to/env-file
```

## Command-line Arguments

The integration test scripts and parallelization test scripts accept the following command-line arguments

`-d /path/to/test/dir` specifies the root directory where integration test subdirectories and scripts will be placed.

`-e /path/to/env-file` specifies the file that is used to initialize the software environment.  This file will typically contain either `module load` or `spack load` commands.

`-l` specifies that the LSF scheduler will be used to run the integration test scripts.

`-p partition` specifies the partition for the SLURM or LSF scheduler.  You may omit this for running interactively

`-s` specifies that the SLURM scheduler will be used to run the integration test scripts.

You can also use long names for the option switches:
- `--directory` instead of `-d`
- `--env-file` instead of `-e`
- `--lsf` instead of `-l`
- `--partition` instead of `-p`
- `--slurm` instead of `-s`

There is an additional option (`-q` or `--quick`) that will run only a couple of integration tests instead of the full suite.  This is intendend for development and debugging.  You will normally not need to use this option.
