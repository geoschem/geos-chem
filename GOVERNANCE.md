# Governance

GEOS-Chem is a grass-roots, open-access model owned by its [users](https://geoschem.github.io/community.html). This document summarizes how development priorities are set and how code moves from a contributor's changes into an official release.

## Roles

### GEOS-Chem Steering Committee (GCSC)

The [GEOS-Chem Steering Committee](https://geoschem.github.io/steering-committee.html) provides scientific direction for the model. It meets quarterly to set [GEOS-Chem model development priorities](http://wiki.geos-chem.org/GEOS-Chem_model_development_priorities) and to slate submitted updates for inclusion into upcoming GEOS-Chem versions.

### Working Groups

[GEOS-Chem Working Groups](https://geoschem.github.io/working-groups.html) identify priority development needs within their science areas and advise the GCSC accordingly. Working Group chairs are the first point of contact for contributors who want to propose a model update, and can advise on timing for submitting code tied to a mature development.

### GEOS-Chem Support Team (GCST)

The [GEOS-Chem Support Team](https://geoschem.github.io/support-team.html), based at Harvard University and Washington University in St. Louis, manages this codebase day to day: it reviews and merges pull requests into the development branch for an upcoming version, and validates incoming updates with [benchmark simulations](http://wiki.geos-chem.org/GEOS-Chem_benchmarking) before release. If a benchmark reveals a problem with a submitted update, the GCST works with the contributor on corrective action.

## How a change becomes an official release

1. A contributor proposes an update to their relevant Working Group chair(s).

2. The Working Group forwards the request to the GCSC, which sets its priority and target version at a quarterly meeting.

3. The contributor submits the update as a pull request against this repository. Structural (non-science) changes should include a difference test (see `test/difference/`) confirming bit-for-bit identical results, and all changes should include a `CHANGELOG.md` entry.

4. The GCST reviews and merges the update into the development branch, then benchmarks and includes it in the next tagged release. Releases here are picked up as pinned submodule updates by the [GCClassic](https://github.com/geoschem/GCClassic) and [GCHP](https://github.com/geoschem/GCHP) wrapper repositories.

## Sponsors

GEOS-Chem development is supported by the US NASA Earth Science Division, the Canadian National Sciences and Engineering Research Council, and the Nanjing University of Information Science and Technology.
