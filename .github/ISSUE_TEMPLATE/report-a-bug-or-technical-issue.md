---
name: Report a bug or technical issue
about: This template allows users to report GEOS-Chem bugs and technical issues in the Github issue tracker.
title: "[BUG/ISSUE]"
labels: bug
assignees: ''

---

# Report a GEOS-Chem bug or technical issue

## First check to see if your issue has already been resolved
Before submitting a GEOS-Chem bug/issue report, please take a moment to check if a solution for your issue has already been reported at:

1. [The issue tracker attached to this Github repository](https://github.com/geoschem/geos-chem/issues);
2. Our [Guide to GEOS-Chem error messages](http://wiki.geos-chem.org/Guide_to_GEOS-Chem_error_messages) on the GEOS-Chem wiki; and/or
3. Our [list of GEOS-Chem bugs and fixes](http://wiki.geos-chem.org/Bugs_and_fixes) on the GEOS-Chem wiki.

## Describe the bug
Include a clear and concise description of the bug or issue that you have encountered.

## To Reproduce
Include the steps that must be done in order to reproduce the observed behavior:

**Compilation commands**
1. Step 1
2. Step 2
3. ... etc ...

**Run commands**
1. Step 1
2. Step 2
3. ... etc ...

## Expected behavior
Include a clear and concise description of what you expected to happen.

## Error messages
```
Cut and paste any error output here.
```

## Required information
Please include the following:
 - GEOS-Chem version you are using [e.g. 12.3.2]
 - Compiler version that you are using [e.g. gfortran 8.2.0, ifort 17.0.4] 
 - netCDF and netCDF-Fortran library version
 - Computational environment [e.g. a cluster, or the AWS cloud]
 - The Amazon Machine Image (AMI) ID that you used (if you ran on the AWS cloud)
 - Are you using "out of the box" code, or have you made modifications?

## Input and log files to attach
For more info, see: http://wiki.geos-chem.org/Submitting_GEOS-Chem_support_requests
 - The lastbuild file
 - The input.geos file
 - The HEMCO_Config.rc file
 - The GEOS-Chem "Classic" log file
 - The HEMCO.log file
 - Error output from your scheduler, if applicable [e.g. slurm*.out]
 - Any other error messages

## Additional context
Include any other context about the problem here.
