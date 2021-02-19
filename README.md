[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1343546.svg)](https://doi.org/10.5281/zenodo.1343546) [![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/geoschem/geos-chem/blob/master/LICENSE.txt)

# README for the GEOS-Chem science codebase repository

## Overview

This repository (https://github.com/geoschem/geos-chem) contains the __GEOS-Chem science codebase__ (aka the GEOS-Chem chemical module).  Included in this repository are:

  * The source code for GEOS-Chem science routines;

  * Scripts to create GEOS-Chem run directories;

  * Scripts to run GEOS-Chem tests;

  * Driver routines (e.g. "main.F90") that enable GEOS-Chem to be run in several different contexts, including:

    * __GEOS-Chem "Classic"__ (suitable for running with less than ~64 cores on a single node)
    * __GCHP__ (GEOS-Chem "High Performance, which can run across multiple nodes for high scalability)
    * __GEOS-GC__ (GEOS-Chem coupled to the NASA GEOS ESM)
    * __WRF-GC__ (GEOS-Chem coupled with WRF)
    * and GEOS-Chem coupled to other Earth System Models (e.g. CESM2, NCAR "next-generation" models, etc.)

## HEMCO is no longer part of the geoschem/GEOS-Chem repository

One important note: The __Harmonized Emissions Component (HEMCO)__ source code has been removed, and is maintained as a separate repository (geoschem/HEMCO). This is because HEMCO is now being developed independently of GEOS-Chem.  This will facilitate HEMCO being used in external models (such as at NOAA and NCAR).  

## Wrapper repositories

You should not clone geoschem/geos-chem directly.  Instead, you should clone one of the "wrapper" repositories:

  * __geoschem/GCClassic__

  * __geoschem/GCHP__

which will then treat geoschem/geos-chem (and geoschem/HEMCO for that matter) as a Git submodule.  For details, we refer you to our "Getting Started with GEOS-Chem" manual or our Youtube tutorial videos (discussed in the Documentation section below).

## GEOS-Chem Development

### Branches
This repository contains several branches.  Each branch contains code updates belonging to a particular line of development.

 * The __main__ branch always contains the __current stable version__.  You should never add new code directly into this branch.  Instead, open a new branch off of main and add your code there.

 * The __dev/X.Y.Z__ branches always contains in-development code for upcoming version X.Y.Z.  Code in dev/X.Y.Z is very much "work in progress" and should not be relied upon until it has been fully debugged, validated, and merged back into the master branch.

 * From time to time, you will see other branches pertaining to new lines of development being created.  These branches usually will start with __feature/__ or __bugfix/__.  Once the code in these branches has been sufficiently validated, these branches will be merged back into the master branch.  You should not use code in these branches.

 * You may also see branches beginning with e.g. __GEOS__, __CESM__, etc.  These branches are intended for ongoing development to connect GEOS_Chem within other Earth System Models.  You may ignore the code in these branches.

### Versions

GEOS-Chem versions are now denoted by 3 digits (X.Y.Z):

 * The __X__ digit is the __major version number__.  A change in X denotes that the current version contains a significant update that breaks backwards-compatibility with the prior series of GEOS-Chem versions.  Major structural updates typically will require an update to X.  In the past we have changed the X digit when replacing SMVGEAR with FlexChem (version 10 to version 11) and replacing legacy emissions with HEMCO (version 9 to version 10).

* The __Y__ digit is the __feature version number__.  A change in Y denotes that a 1-month benchmark has been performed to validate a new feature or set of features.  Some (but not all) Y versions will have 1-year benchmarks performed as well.  In general, the Y digit changes whenever a new feature  breaks backwards compatibility with one or more run directories from the prior version.

* The __Z__ digit is the __bug fix (or patch) version number__.   A change in Z denotes that a bug fix was made that does not break backwards compatibility with run directories from the prior verison.  Z will also be updated when bug fixes or minor updates are made to one or more of the GEOS-Chem "Specialty" simulations.  Updating specialty simulations should not affect the output of the GEOS-Chem 1-month or 1-year benchmark simulations.

For more information, please see this wiki page: http://wiki.geos-chem.org/GEOS-Chem_version_numbering_system

All benchmarked GEOS-Chem versions are tagged in the Git history. Use _git tag_ in your terminal to see a list of available tags. Tags will also be highlighted in the _gitk_ browser window.

### Citing GEOS-Chem versions with DOI's

You can now cite GEOS-Chem in publications with a Digital Object Identifier (DOI). Each GEOS-Chem release will be assigned its own individual DOI.  DOI's for each GEOS-Chem version will be posted on the [GEOS-Chem website](http://geos-chem.org) and [GEOS-Chem wiki](http://wiki.geos-chem.org).

We have also generated a concept DOI, which will always point to the current stable version of GEOS-Chem (i.e. corresponding to the __master__ branch): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1343546.svg)](https://doi.org/10.5281/zenodo.1343546)

# Documentation

### Web site
The __GEOS-Chem web site__ is a good place to get started.  It will point you to many important GEOS-Chem resources.

  * http://www.geos-chem.org

### Online user's manual
You can find the __The GEOS-Chem User's Guide__ online here:

  * http://manual.geos-chem.org

NOTE: The above link currently points to __Getting Started with GEOS-Chem__ on the GEOS-Chem wiki.  In the near future, we will migrate this documentation to readthedocs.io.

### Wiki
The most up-to-date information about GEOS-Chem is posted on the __GEOS-Chem wiki__.  Here you will find information about technical issues, bug fixes, and other pertinent topics.

  * http://wiki.geos-chem.org

### Youtube
We have created several tutorial videos at our GEOS-Chem Youtube channel (https://youtube.com/c/geos-chem).  Please take a few moments to view the following tutorials

  * [Getting started with GEOS-Chem "Classic" 13.0.0](https://www.youtube.com/watch?v=BV4BIj8WAxE)

  * [Getting started with GCHP 13.0.0 -- Building GCHP](https://www.youtube.com/watch?v=G_DMCv-mJ2k)
  * [Getting started with GCHP 13.0.0 -- Running GCHP](https://www.youtube.com/watch?v=K6frcfCjpds)

  * [Getting started with the HEMCO 3.0.0 standalone](https://www.youtube.com/watch?v=6Bup9V0ts6U&t=25s)

  * [Getting started with GCPy](https://www.youtube.com/watch?v=eC6203eF05g)

## Support
We encourage GEOS-Chem users to use [the Github issue tracker attached to this repository](https://github.com/geoschem/geos-chem/issues/new/choose) to report bugs or technical issues with the GEOS-Chem code.

## License

GEOS-Chem (and related software) is distributed under the MIT license. Please see the license documents LICENSE.txt and AUTHORS.txt in the root folder.


06 Jan 2021
GEOS-Chem Support Team
geos-chem-support@g.harvard.edu
