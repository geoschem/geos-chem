# MAPL3 run directory configuration files (developers only)

**These files are for GCHP development with MAPL3, which is in progress and
not supported for general use.** Standard GCHP builds use MAPL2 and the
configuration file templates in the parent `run/GCHP` folder, not these.

Each `case.*` folder holds the MAPL3 configuration files (`*.yaml` and, where
needed, `input.nml` and `HISTORY.rc`) for one test scenario:

| Folder | Scenario |
|---|---|
| `case.no_children` | GCHP gridded component with no children |
| `case.one_child` | GCHP with GCHPctmEnv only |
| `case.two_children` | GCHP with GCHPctmEnv and AdvCore (advection) |
| `case.transport_tracers` | GCHP with GCHPctmEnv, AdvCore and GEOS-Chem (transport tracers simulation) |

To create a MAPL3 run directory, set `GCHP_MAPL3_DEV=y` before running
`createRunDir.sh` in `run/GCHP`, then choose option 5 and a test scenario:

```console
GCHP_MAPL3_DEV=y ./createRunDir.sh
```

The run directory must be used with a GCHP build configured with `-DMAPL3=y`
and with MAPL3-compatible versions of the MAPL and FVdycoreCubed_GridComp
submodules, which are not the versions pinned in standard GCHP releases.
