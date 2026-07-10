# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

This is the **GEOS-Chem science codebase** (`geoschem/geos-chem`) — the Fortran source for GEOS-Chem, a global 3-D model of atmospheric chemistry. It is almost never built or tested standalone; it is consumed as a git submodule by two superproject wrapper repos:

- **GCClassic** (`geoschem/GCClassic`) — builds this code as a standalone executable ("GEOS-Chem Classic")
- **GCHP** (`geoschem/GCHP`) — builds this code as an ESMF/MAPL gridded component inside the GEOS/NASA modeling framework ("GCHP")

Both superprojects vendor this repo at `src/GEOS-Chem`, and symlink `run/`, `test/`, and `spack/` from it up to their own top level. If you were pointed here from a GCClassic or GCHP checkout, you are actually editing *this* repo — commits/PRs belong here (geoschem/geos-chem), not in the wrapper repo.

Preprocessor macros `MODEL_CLASSIC` / `MODEL_GCHP` (set by the superproject's CMake) gate code paths that only apply to one implementation (e.g. `Interfaces/GCClassic/main.F90` is wrapped in `#ifdef MODEL_CLASSIC`). When reading or editing shared modules, check which macros surround a block before assuming it runs in both implementations.

## Repository layout

Source is organized by role, not by scientific topic — a given "feature" (e.g. dry deposition) typically spans several of these:

| Directory | Contents |
|---|---|
| `GeosCore/` | Core science drivers: chemistry, transport, convection, deposition, emissions coupling, per-simulation-type modules (`tagged_o3_mod.F90`, `mercury_mod.F90`, `tomas_mod.F90`, `carbon_mod.F90`, etc.) |
| `Headers/` | Shared derived types and utilities used everywhere: `state_chm_mod.F90`, `state_met_mod.F90`, `state_grid_mod.F90`, `state_diag_mod.F90`, `species_database_mod.F90`, `species_mod.F90`, `precision_mod.F90` |
| `History/` | The netCDF diagnostics framework (`history_mod.F90` and the `Hist*` container types) that HISTORY.rc/diagnostics are built on |
| `KPP/` | Chemical mechanisms. Each subdirectory (`fullchem/`, `carbon/`, `Hg/`, `custom/`, `aciduptake/`, `stubs/`) is a KPP-generated solver plus mechanism-specific Fortran (e.g. `fullchem_RateLawFuncs.F90`, `fullchem_SulfurChemFuncs.F90`). `KPP/standalone` is a separate submodule (KPP-Standalone box model, built when `-DKPPSA=y`). Never hand-edit KPP-generated solver files (`gckpp_*`) — regenerate with `build_mechanism.sh` |
| `Interfaces/GCClassic/` | GEOS-Chem Classic driver (`main.F90`), built only when `MODEL_CLASSIC` |
| `Interfaces/GCHP/` | GCHP gridded-component glue (`Chem_GridCompMod.F90`, `gchp_chunk_mod.F90`), built only when `MODEL_GCHP` |
| `GeosUtil/`, `NcdfUtil/` | Generic utilities (string parsing, netCDF I/O wrappers, YAML config parsing via `qfyaml_mod.F90`) with no science content |
| `GTMM/`, `APM/`, `GeosRad/`, `PKUCPL/`, `ObsPack/` | Optional/pluggable components (Global Terrestrial Mercury Model, aerosol microphysics, RRTMG radiative transfer, ObsPack diagnostics) enabled by their own CMake switches |
| `run/GCClassic/`, `run/GCHP/`, `run/shared/`, `run/CESM/`, `run/GEOS/`, `run/WRF/` | Run-directory creation scripts (`createRunDir.sh`) and config-file templates (`geoschem_config.yml`, `HEMCO_Config.rc`, `HISTORY.rc`, etc.) for each implementation |
| `test/` | Integration/parallel/difference test drivers (see below) |
| `CMakeScripts/` | `GC-Helpers.cmake` — shared CMake macros (`gc_pretty_print`, version detection) used by this repo's `CMakeLists.txt` |

`CMakeLists.txt` at the repo root is included by the superproject's build, not invoked standalone — it expects `GEOSChemBuildProperties` and options like `MECH`, `MODEL_CLASSIC`/`MODEL_GCHP` to already be defined by the caller.

## Building

There is no standalone build here — always build via a superproject run directory. From a GCClassic or GCHP checkout with this repo as its `src/GEOS-Chem` submodule:

```console
cd run/GCClassic && ./createRunDir.sh      # or run/GCHP/createRunDir.sh
cd /path/to/rundir/build
cmake ../CodeDir -DRUNDIR=..
make -j && make install
```

Relevant CMake options that live in *this* repo (`KPP/CMakeLists.txt`, root `CMakeLists.txt`):
- `MECH` — `fullchem` (default), `carbon`, `custom`, `Hg` — selects which `KPP/<mech>` subdirectory is built
- `KPPSA` — also builds `KPP/standalone` (KPP-Standalone box model) alongside `fullchem`/`custom`
- Supported Fortran compilers are Intel and GNU only (`GEOSChem_Fortran_FLAGS_{Intel,GNU}`); anything else is a hard CMake `FATAL_ERROR`

## Modifying a chemical mechanism

Do not hand-edit the KPP-generated solver files under `KPP/<mechanism>/` (files prefixed `gckpp_`). Instead edit the mechanism's `.eqn`/`.spc`/`.kpp` definition files, then regenerate:

```console
cd KPP
./build_mechanism.sh fullchem   # or Hg, custom
```

This requires KPP 3.4.0+ (see root `CHANGELOG.md` for the currently-required version) and preserves the heterogeneous-chemistry files while regenerating the solver. Hand-written mechanism support code (rate laws, heterogeneous chemistry hookups) lives alongside the generated files, e.g. `KPP/fullchem/fullchem_RateLawFuncs.F90`, `KPP/fullchem/fullchem_SulfurChemFuncs.F90`, `KPP/fullchem_HetStateFuncs.F90` / `KPP/stubs/stub_fullchem_HetStateFuncs.F90` (a stub is required for every mechanism that doesn't implement heterogeneous chemistry).

## Testing

Test infrastructure lives under `test/` and is symlinked to the same path in GCClassic/GCHP superprojects. All test scripts require a real (non-conda) netCDF — `conda deactivate` first, or the scripts abort immediately.

- `test/integration/GCClassic/` and `test/integration/GCHP/` — compile + optionally run several out-of-the-box run-directory configurations, to catch build/run regressions:
  ```console
  cd test/integration/GCClassic
  ./integrationTest.sh -d <root-dir> -t compile   # compile-only
  ./integrationTest.sh -d <root-dir> -t all        # compile and run short sims
  ./integrationTest.sh -d <root-dir> -t compile -q # quick subset, for local dev
  ```
- `test/parallel/GCClassic/` — same idea but sweeps OpenMP thread counts, to catch parallelization bugs.
- `test/difference/` — bit-for-bit diff of output between two integration/parallel test runs (e.g. before/after a structural, non-science change).
- `test/shared/commonFunctionsForTests.sh` — functions/settings shared by all of the above; source this rather than duplicating logic when adding new test scripts.

## Versioning and changes

- Root `CHANGELOG.md` follows Keep a Changelog / SemVer and documents changes to this repo specifically (GCClassic's own CHANGELOG.md separately tracks submodule-pointer bumps and wrapper-level changes).
  Individual KPP mechanisms may keep their own changelog, e.g. `KPP/fullchem/CHANGELOG_fullchem.md` — update that too when changing fullchem's chemistry.
- This is a community-governed, grass-roots model (see `README.md` and geos-chem.org). Substantive science/structural changes are expected to go through the GEOS-Chem Steering Committee / User Working Group process, not just a code review.
- Any structural (non-science) change should be accompanied by a difference test (`test/difference/`) against the prior version to confirm bit-for-bit identical results.
- Config/run-directory changes should be mirrored across `run/GCClassic/` and `run/GCHP/` (and `run/CESM`, `run/GEOS`, `run/WRF` where applicable) since they share the same underlying `geoschem_config.yml` / `HEMCO_Config.rc` / `HISTORY.rc` schema.
