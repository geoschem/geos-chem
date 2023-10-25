# README: The `src/GEOS-Chem/test` directory

## Contents

`difference/`
- Directory containing scripts for performing a **difference test**.
- A difference test compares the output of two different integration tests for identicality.
- TODO (Jun 2023): Add capability to compare the output of two **parallel tests** for identicality.

`integration/`
- Directory containing scripts for performing **integration tests** for both GEOS-Chem Classic and GCHP.
- Integration tests compile and run several different "out-of-the-box" GEOS-Chem Classic and GCHP configurations in order to identify errors caused by source code updates.

`parallel/`
- Directory contianing scripts for performing **parallel tests** for GEOS-Chem Classic.
- Parallel tests compile and run several different "out-of-the-box" GEOS-Chem Classic configurations with different numbers of OpenMP threads in order to identify parallelization errors caused by source code updates.

`shared/`
- Directory containing scripts with common settings and functions that are used by the difference test, integration test, and parallel test scripts.
