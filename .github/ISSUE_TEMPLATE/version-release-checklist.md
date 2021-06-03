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

- GEOS-Chem classic benchmarks:
  - [ ] Integration tests
  - [ ] 1-month full-chemistry benchmark
  - [ ] 1-year full-chemistry benchmark
    - [ ] Archive PCO and PCH4 fields for tagged CO simulations in `ExtData/HEMCO/TAGGED_CO/`
    - [ ] Archive restart files in `ExtData/GEOSCHEM_RESTARTS/`
  - [ ] 1-year transport tracer benchmark
  - [ ] 10-year stratospheric benchmark
    - Archive monthly mean species concentrations in `ExtData/CHEM_INPUTS/`
    - Archive monthly P/L rates in `ExtData/HEMCO/UCX/`
- [ ] GCHP benchmarks:
  - [ ] 1-month full-chemistry benchmark
  - [ ] 1-year full-chemistry benchmark
  - [ ] 1-year transport tracer benchmark
  - [ ] Stratospheric benchmark
- [ ] Send benchmarks to developers and GCSC for approval

## Github repositories
- [ ] Push tags for this version to all repositories:
  - [ ] https://github.com/geoschem/geos-chem
  - [ ] https://github.com/geoschem/gchp
  - [ ] https://github.com/geoschem/geos-chem-unittest
  - [ ] https://github.com/geoschem/hemco
  - [ ] https://github.com/geoschem/gcpy
- [ ] Publish new release(s) on Github

## Data
- [ ] Document updated data directories for current version on wiki
- [ ] Copy updated data files to Compute Canada

## Wiki and website updates
- [ ] Update [Species in GEOS-Chem](http://wiki.seas.harvard.edu/geos-chem/index.php/Species_in_GEOS-Chem)
- [ ] Update version number on the following pages:
  - [ ] [GEOS-Chem website](http://acmg.seas.harvard.edu/geos/)
  - [ ] [GEOS-Chem wiki](http://wiki.seas.harvard.edu/geos-chem/index.php/Main_Page)
  - [ ] [GEOS-Chem versions](http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_versions)
  - [ ] [GEOS-Chem 12](http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_12#Version_history_summary)

## YouTube tutorials
- [ ] TBD

## Announcement
- [ ] Send email to user list announcing release