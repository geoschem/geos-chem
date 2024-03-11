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

## Before you begin

Please take a moment to:

1. Verify that the `GCClassic` superproject is checked out to the correct branch and commit.
2. Run `git submodule update --init --recursive` in order to update all submodules.
3. Verify that the `HEAD` commit of the `GEOS-Chem` submodule contains the code that you wish to test. (If not, then check out the proper branch.)
4. Verify that the `HEAD` commit of the `HEMCO` submodule contains the code that you wish to test. (If not, then check out the proper branch.)

## Command-line Arguments

The integration test scripts accept the following command-line arguments:

### Required arguments

`-d /path/to/test/dir` specifies the root directory where integration test subdirectories and scripts will be placed.

`-t compile|all` specifies the type of test to be run:
  - `compile` will run compilation-only tests.
  - `all` will run compilation and execution tests.

### Optional arguments

`-e /path/to/env-file` Specifies the file that is used to initialize the software environment on the Harvard Cannon cluster.  If omitted, a default file will be selected.

`-h` displays a help screeen.

`-q` will run only a couple of integration tests instead of the full suite.  This is intended for development and debugging.  You will normally not need to use this option.

You can also use long names for the option switches:
- `--directory` instead of `-d`
- `--env-file` instead of `-e`
- `--help` instead of `h`
- `--no-bootstrap` instead of `n`
- `--tests-to-run` instead of `-t`

## Examples

### Request compile-only tests on Cannon or Compute1
```console
$ cd /path/to/GCClassic     # Path to GCClassic superproject directory
$ cd test/integration/GCClassic
$ ./integrationTest.sh -d /path/to/test/dir -t compile
```

### Request compile and execution tests (Harvard Cannon)
```console
$ cd /path/to/GCClassic     # Path to GCClassic superproject directory
$ cd test/integration/GCClassic
$ ./integrationTest.sh -d /path/to/test/dir -t all -e /path/to/env-file
```
NOTE: If you omit the `-e /path/to/env/file` a default environment file will be used to load GNU Compiler Collection 10 and related libraries.

### Request compile & execution tests (WashU compute1)
```console
$ cd /path/to/GCClassic     # Path to GCClassic superproject directory
$ cd test/integration/GCClassic
$ ./integrationTest.sh -d /path/to/test/dir -t all
```
NOTE: No environment file is needed.  On Compute1 the tests will run inside a software container with all necessary libraries included.
