---
name: GEOS-Chem version release checklist
about: Use this template to create a checklist for a new GEOS-Chem version release.
title: "[RELEASE-CHECKLIST]"
labels: 'checklist'
assignees: ''

---

# Checklist for releasing GEOS-Chem X.Y.0

## Benchmarks

The following benchmark simulations are needed to validate this version (edit as needed). For items marked at complete, benchmark simulations have finished, benchmark plots have been created and posted on the FTP site, and benchmark information has been posted on the GEOS-Chem wiki.

GEOS-Chem classic:
  - [ ] Complete unit test
  - [ ] 1-month full-chemistry benchmark
  - [ ] 1-year full-chemistry benchmark
  - [ ] 1-year transport tracer benchmark
  - [ ] Stratospheric benchmark

GCHP:
  - [ ] 1-month full-chemistry benchmark
  - [ ] 1-year full-chemistry benchmark
  - [ ] 1-year transport tracer benchmark
  - [ ] Stratospheric benchmark

[ ] Send benchmarks to developers and GCSC for approval
[ ] Update [GEOS-Chem benchmarks](http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_Benchmarks) 

## Github repositories
[ ] Push tags for this version to all repositories:
- [ ] GEOS-Chem classic
- [ ] GCHP
- [ ] Unit tester
- [ ] HEMCO

[ ] Publish new release on Github repository

## Data
[ ] Document updated data directories for current version on wiki
[ ] Copy updated data files to Compute Canada
[ ] Archive PCO and PCH4 fields from 1-year benchmark for tagged CO simulation
[ ] Archive restart files from 1-year full-chemistry benchmark

## Wiki and website updates
[ ] Update version number on the following pages:
- [ ] [GEOS-Chem website](http://acmg.seas.harvard.edu/geos/)
- [ ] [GEOS-Chem wiki](http://wiki.seas.harvard.edu/geos-chem/index.php/Main_Page)
- [ ] [GEOS-Chem versions](http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_versions)
- [ ] [GEOS-Chem 12](http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_12#Version_history_summary)

[ ] Send email to user list announcing release