[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1343546.svg)](https://doi.org/10.5281/zenodo.1343546) [![Build
Status](https://travis-ci.org/JiaweiZhuang/geos-chem.svg?branch=travis_ci)](https://travis-ci.org/JiaweiZhuang/geos-chem) [![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/geoschem/geos-chem/blob/master/LICENSE.txt)

# README for the GEOS-Chem source code repository

This repository (https://github.com/geoschem/geos-chem) contains the source code for the GEOS-Chem model of atmospheric chemistry and composition. 

## GEOS-Chem Development

### Branches
This repository contains several branches.  Each branch contains code updates belonging to a particular line of development.

 * The __master__ branch always contains the __current stable version__.  You should never add new code directly into this branch.  Instead, open a new branch off of master and add your code there.

 * The __dev/X.Y.Z__ branches always contains in-development code for upcoming version X.Y.Z.  Code in dev/X.Y.Z is very much "work in progress" and should not be relied upon until it has been fully debugged, validated, and merged back into the master branch.

 * The __GEOS__ branch contains updates that are specific to the interface between GEOS-Chem and the NASA GEOS-DAS Earth System Model.  Most GEOS-Chem users can simply ignore this branch.

 * From time to time, you will see other branches pertaining to new lines of development being created.  These branches usually will start with __feature/__ or __bugfix/__.  Once the code in these branches has been sufficiently validated, these branches will be merged back into the master branch.  You should not use code in these branches.

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

### Wiki
The most up-to-date information about GEOS-Chem is posted on the __GEOS-Chem wiki__.  Here you will find information about technical issues, bug fixes, and other pertinent topics.

  * http://wiki.geos-chem.org

## GEOS-Chem run directories
To generate GEOS-Chem run directories, please clone the [__geos-chem-unittest__](https://github.com/geoschem/geos-chem-unittest) repository and follow the instructions as listed on the [Creating GEOS-Chem run directories wiki page](http://wiki.seas.harvard.edu/geos-chem/index.php/Creating_GEOS-Chem_run_directories).

## Support 
We encourage GEOS-Chem users to use the Github issue tracker attached to this repository to report  bugs or technical issues with the GEOS-Chem code.

You are also invited to direct GEOS-Chem support requests to the GEOS-Chem Support Team at geos-chem-support@g.harvard.edu.

## License

GEOS-Chem (and related software) is distributed under the MIT license. Please see the license documents LICENSE.txt and AUTHORS.txt in the root folder.


14 Nov 2018
GEOS-Chem Support Team
geos-chem-support@g.harvard.edu
