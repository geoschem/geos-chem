# README for Integration Tests

## Overview:

This directory contains scripts to submit automated GCHP integration tests, which will test the following process:

1. GCClassic run directory creation
2. GCClassic code configuration and compilation
3. GCClassic execution

Integration tests are mostly intended to be run by GEOS-Chem developers, rather than by end users.

## Files

- `intTest.sh`
  - Driver script to run GCHP integration tests
- `intTestCreate.sh`
  - Script to create GCHP run directories for integration tests
- `intTestCompile.sh`
  - Script to compile GCHP Classic executables
- `intTestExecute.sh`
   - Script to run GCHP Classic simulations
- `commonFunctionsForTests.sh`
  - Link to `../shared/commonFunctionsForTests.sh`, which contains global variables and functions for the integration test scripts.

The integration tests are all short simulations (1 model hour).  Executables are compiled with all debugging options.  This should reveal any coding errors or run-directory configuration errors.

## Running GCHP Integration Tests

Before you submit any GEOS-Chem Classic integration tests, take a moment and make sure that:

1. The `GCHP` superproject is checked out to the correct branch and commit.
2. The `GEOS-Chem` submodule is checked out to the correct branch and commit.
3. The `HEMCO` submodule is checked out to the correct branch and commit.

### With the SLURM scheduler

```console
$ cd gchp/test
$ ./intTest.sh -d /path/to/test/dir -e /path/to/env-file -s -p partition
```

### With the LSF scheduler

```console
$ cd gchp/test
$ ./intTest.sh -d /path/to/test/dir -e /path/to/env-file -l -p partition
```

### Interactively at the command line

```console
$ cd gchp/test
$ ./intTest.sh -d /path/to/test/dir -e /path/to/env-file
```

### Notes

1. `-d /path/to/test/dir` specifies the root directory where integration test subdirectories and scripts will be placed.
2. `-e /path/to/env-file` specifies the file that is used to initialize the software environment.  This file will typically contain either `module load` or `spack load` commands.
3. `-p partition` specifies the partition for the SLURM or LSF scheduler.  You may omit this for running interactively
4. You can also use long names for the option switches:
   - `--directory` instead of `-d`
   - `--env-file` instead of `-e`
   - `--lsf` instead of `-l`
   - `--partition` instead of `-p`
   - `--slurm` instead of `-s`
5. There is an additional option (`-q` or `--quick`) that will run only a couple of integration tests instead of the full suite.  This is intendend for development and debugging.  You will normally not need to use this option.
