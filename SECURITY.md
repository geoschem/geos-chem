# Security Policy

## Supported Versions

GEOS-Chem does not maintain long-term-support branches. Security fixes are only provided for the most recently released version, listed in `CHANGELOG.md`.

## Reporting a Vulnerability

If you believe you have found a security vulnerability in this repository — for example, a supply-chain issue in a GitHub Actions workflow, or an issue in a run-directory/data-download script (`run/shared/download_data.py`, `createRunDir.sh`, etc.) that could lead to unintended code execution — please report it privately using GitHub's **[Report a vulnerability](https://github.com/geoschem/geos-chem/security/advisories/new)** feature (Security tab) rather than opening a public issue.

The threat class that matters most here is **arbitrary code execution when reading a data/config file**. GEOS-Chem and its helper scripts read many files that users download or share: YAML configuration files (`geoschem_config.yml` and others, parsed by `Headers/qfyaml_mod.F90`), `run/shared/download_data.yml` (read by `run/shared/download_data.py`), `HEMCO_Config.rc` and `HISTORY.rc`, and NetCDF input files. If any of these can be crafted to make the model or a script run code or shell commands, please report it.

If the issue is specific to the wrapper repositories that consume this codebase (GCClassic or GCHP) or to another submodule (HEMCO, Cloud-J, HETP), please report it in that repository instead.

This project is maintained by the **GEOS-Chem Support Team (GCST)** on a best-effort basis, so there is no guaranteed response SLA, but we will acknowledge reports as promptly as we can and work with you on a fix and coordinated disclosure.

## Out of Scope

Scientific-correctness bugs, numerical issues, and general "how do I..." questions are **not** security reports. Please use the normal channels — [GitHub issues](https://github.com/geoschem/geos-chem/issues/new/choose) or the [GEOS-Chem user manual](https://geos-chem.readthedocs.io/en/stable) — for those instead.
