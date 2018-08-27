[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1343546.svg)](https://doi.org/10.5281/zenodo.1343546)
# README for the GEOS-Chem Source code repository

[![Build
Status](https://travis-ci.org/JiaweiZhuang/geos-chem.svg?branch=travis_ci)](https://travis-ci.org/JiaweiZhuang/geos-chem)

This repository (https://github.com/gcst/geos-chem) contains the source code for the GEOS-Chem model of atmospheric chemistry and composition. 

---
---
---
## Important notes

### We have migrated from Bitbucket to Github!
As of June 2018, we have migrated the GEOS-Chem source code repository to back Github.  Going forward, please make sure to clone or pull code updates ONLY from this repository.

### The GEOS-Chem version numbering scheme is changing!

We are migrating to a purely numeric versioning system in order to adhere more closely to software development best practices. In order to facilitate the switchover between versioning systems, the __GEOS-Chem v11-02-final__ version will also carry the designation __GEOS-Chem 12.0.0__.

For a complete description of the new versioning system, please see our GEOS-Chem version numbering system wiki page: http://wiki.geos-chem.org/GEOS-Chem_version_numbering_system

---
---
---

## GEOS-Chem Development

### Branches
This repository contains several branches, each of which contains several code updates belonging to a particular line of development.  In particular you will see:

 * The __master__ branch always contains the __current stable version__.  You should never add new code directly into this branch.  Instead, open a new branch off of master and add your code there.

 * The __Dev__ branch always contains in-development code for the next version to be benchmarked.  Code in Dev is very much "work in progress" and should not relied upon until it has been debugged and benchmarked.

 * The __GEOS__ branch contains updates that are specific to the interface between GEOS-Chem and the NASA GEOS-DAS Earth System Model.  Most GEOS-Chem users can simply ignore this branch.

 * From time to time, you will see other branches pertaining to new lines of development being created.  Once the code in these branches has been sufficiently validated, these branches will be merged back into the master branch.  You should not use code in these branches.

### Versions

GEOS-Chem versions correspond to those points in the code development history where benchmark simulations were performed in order to validate the scientific output of GEOS-Chem.  All versions are benchmarked with a 1-month simulation; some versions are benchmarked with both 1-month and 1-year simulations.  To learn more about GEOS-Chem versions and the features and/or bug fixes they contained, please visit the GEOS-Chem wiki pages listed below.

For __GEOS-Chem 12.\*.\*__ versions, see:
* http://wiki.geos-chem.org/GEOS-Chem_12
* http://wiki.geos-chem.org/GEOS-Chem_12_benchmark_history

All benchmarked GEOS-Chem versions are tagged in the Git history. Use _git tag_ in your terminal to see a list of available tags. Tags will also be highlighted in the _gitk_ browser window.

### Citing GEOS-Chem versions with DOI's

You can now cite GEOS-Chem in publications with a Digital Object Identifier (DOI). Each GEOS-Chem release will be assigned its own individual DOI.  DOI's for each GEOS-Chem version will be posted on the [GEOS-Chem website](http://geos-chem.org) and [GEOS-Chem wiki](http://wiki.geos-chem.org).

We also have generated a static DOI, which will always point to the current stable version of GEOS-Chem (i.e. corresponding to the __master__ branch): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1343546.svg)](https://doi.org/10.5281/zenodo.1343546)

# Documentation

### Web site
The __GEOS-Chem web site__ is a good place to get started.  It will point you to many important GEOS-Chem resources.
* http://www.geos-chem.org

### Online user's manual
You can find the __The GEOS-Chem User's Guide__ online here:
* http://manual.geos-chem.org

### Wiki
The most up-to-date information about GEOS-Chem is posted on the __GEOS-Chem wiki__.  Here you will find information about technical issues, bug fixes, and other pertinent topics.
* http://wiki-geos.chem.org

## GEOS-Chem run directories
To generate GEOS-Chem run directories, please clone the [__geos-chem-unittest__](https://github.com/geoschem/geos-chem-unittest) repository and follow the instructions as listed in the GEOS-Chem wiki pags listed below.

## Support 
We encourage GEOS-Chem users to use the Github issue tracker attached to this repository to report  bugs or technical issues with the GEOS-Chem code.

You are also invited to direct GEOS-Chem support requests to the GEOS-Chem Support Team at geos-chem-support@as.harvard.edu.

15 Aug 2018
GEOS-Chem Support Team
geos-chem-support@as.harvard.edu
