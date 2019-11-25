---
name: Report a bug or technical issue
about: Use this template to report bugs and technical issues encountered while using GEOS-Chem.
title: "[BUG/ISSUE]"
labels: bug
assignees: ''

---

# Report a GEOS-Chem bug or technical issue

## STOP!  Before you submit an issue, please follow these steps

### Please see our Github tutorial videos
If this is the first time you are submitting a GEOS-Chem issue via Github, we recommend that you first view our tutorial videos on [youtube.geos-chem.org](http://youtube.geos-chem.org).  These videos will demonstrate how to use the Github issue trackers and how to receive notifications about issues that have been posted:

  * [*Submitting GEOS-Chem issues on Github*](https://www.youtube.com/watch?v=dFBhdotYVf8&t=1103s)
  * [*Subscribing to Github notifications*](https://www.youtube.com/watch?v=RuH6zeYuzuY)
  * [*Searching for GEOS-Chem issues and pull requests*](https://www.youtube.com/watch?v=EiZC2vaXNnU&t=35s)

### Check to see if your issue has already been resolved

Before submitting a GEOS-Chem bug/issue report, please take a moment to check if a solution for your issue has already been reported at:

1. [The issue tracker attached to this Github repository](https://github.com/geoschem/geos-chem/issues);
2. Our [Guide to GEOS-Chem error messages](http://wiki.geos-chem.org/Guide_to_GEOS-Chem_error_messages) on the GEOS-Chem wiki; and/or
3. Our [list of GEOS-Chem bugs and fixes](http://wiki.geos-chem.org/Bugs_and_fixes) on the GEOS-Chem wiki.

### Try to resolve the issue yourself
Please see our [Debugging GEOS-Chem wiki page](http://wiki.geos-chem.org/Debugging_GEOS-Chem), which contains a checklist of things that you can do to try to solve your problem yourself.  In particular:

1. Make sure that your [Unix environment](http://wiki.seas.harvard.edu/geos-chem/index.php/Setting_Unix_environment_variables_for_GEOS-Chem) is set up properly;
2. Make sure that all of the [required libraries](http://wiki.geos-chem.org/Guide_to_netCDF_in_GEOS-Chem) have been installed;
3. Look at the error output (GEOS-Chem log, HEMCO.log, etc.) and error messages
4. Recompile GEOS-Chem with [debugging flags](http://wiki.seas.harvard.edu/geos-chem/index.php/Debugging_GEOS-Chem#Recompile_GEOS-Chem_with_debug_options_turned_on) and re-run your simulation;
5. [Use a debugger](http://wiki.geos-chem.org/Debugging_GEOS-Chem#Run_GEOS-Chem_in_a_debugger_to_find_the_source_of_error) to find more information about where the run is crashing;
6. Identify if the error happens at the same date & time on every run, or seems to happen at random;
7. Identify if the error happens in a clean, "out-of-the-box" version of GEOS-Chem, or if it happens in code that you added/modified yourself.
8. Try to isolate the issue to a particular operation (i.e. turn off chemistry, dry-depostion, transport one at a time)
9. Add print statements to your code to try to determine the source of the error.

### Provide all requested information
Please see our [Submitting GEOS-Chem Support Requests](http://wiki.geos-chem.org/Submitting_GEOS-Chem_Support_Requests) wiki page for a checklist of required configuration files and log files to attach to this issue.  The more information you provide, the better able we will be to help you diagnose and fix your issue.

^^^^^^ YOU MAY DELETE ALL TEXT ABOVE THIS LINE WHEN CREATING YOUR ISSUE ^^^^^^

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
