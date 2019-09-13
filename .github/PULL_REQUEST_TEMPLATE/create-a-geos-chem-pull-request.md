---
name: Submit updates to GEOS-Chem with a pull request
about: This template instructs users to submit a GEOS-Chem pull request
title: "[PULL REQUEST]"
labels: ''
assignees: ''

---

# Creating a GEOS-Chem Pull Request

## Overview

The [GEOS-Chem Steering Committee (GCSC)](http://geos-chem.org/geos_steering_cmte.html) encourages updates to the GEOS-Chem model. Updates to the GEOS-Chem model benefit both the developer and the entire GEOS-Chem community. The developer benefits through coauthorship and citations. Priority development needs are identified at   GEOS-Chem users' meetings with updates between meetings based on GCSC input through working groups.

When should you submit updates to the GEOS-Chem code? 

* Bug fixes should be submitted as soon as possible and be given high priority. 

* Code related to model developments should be submitted when it is mature (i.e. a paper has been submitted). Your working group chair can offer guidance on the timing of submitting code to the GCST.

The practical aspects of submitting updates to the GEOS-Chem Support Team are outlined below. 

## Instructions for submitting new code updates

Please make sure that all of the following points have been completed before submitting your pull request.

1. [ ] Contact your [GEOS-Chem Working Group](http://geos-chem.org/geos_working_groups.html) to request that your changes be included in the standard code. The GEOS-Chem Steering Committee will meet every three months and set the GEOS-Chem model development priorities (i.e. decide on the order in which updates will be added to GEOS-Chem).

2. [ ]  Create or log into your [Github account](https://github.com/GitHub).

3. [ ]  [Fork this repoitory (i.e. geoschem/geos-chem)](https://help.github.com/articles/fork-a-repo/). Navigate to each appropriate GEOS-Chem repository on GitHub and click Fork in the top-right corner of the page.

4. [ ] Clone your fork of the GEOS-Chem repository to your system, e.g.  
    ```   git clone https://github.com/YOUR-USERNAME/geos-chem Code.12.6.0```

5. [ ] Add your modifications into a new branch off of the master branch.

6. [ ] Test your updates thoroughly.

    * Use the [GEOS-Chem Unit Tester](http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_Unit_Tester) to make sure your code updates work for a given combination of met fields, resolutions, and simulation types.

     *  For structural updates, we recommend running [difference tests](http://wiki.geos-chem.org/Performing_Difference_Tests_with_GEOS-Chem) often to ensure your updates don't impact model output.

7. [ ] Make sure that your source code adheres to the following style guidelines:
    * [ ] Use F90 free format whenever possible.
    * [ ] Include thorough comments throughout code updates.
    * [ ] Remove extraneous code updates (e.g. testing options, other science)
    * [ ] Include full citations at the top of relevant source code modules
    * [ ] Include any required GCHP updates along with GEOS-Chem Classic updates

8. [ ]  Create a [Github pull request](https://help.github.com/articles/creating-a-pull-request/) (recommended) or a [Git patch file](https://www.devroom.io/2009/10/26/how-to-create-and-apply-a-patch-with-git/).

The GEOS-Chem Support Team will add your changes to the standard code. Your update may be bundled with other code updates and the GEOS-Chem version number (X.Y.Z) will be changed. The update will be validated following the GEOS-Chem benchmarking procedure