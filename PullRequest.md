### Name and Institution (Required)

Name: Daisy Wang / Jordan Sun

Institution: Washington University in St. Louis

### Describe the update
<!-- Please provide a clear and concise overview of the update. -->
#### Describe
GCHP’s chemistry step exhibits **substantial per‑column workload variability** because the KPP solver takes a different number of sub‑steps depending on local atmospheric state and time. When columns are pinned to MPI ranks and ranks are sequentially mapped to cores in nodes, this variability creates **rank‑ and node‑level imbalance** (some cores and subsequenly hosts have mostly “heavy” columns, others mostly “light”), leaving resources idle while stragglers finish. We propose Hierarchical Balancing with MPI Shared Memory approch to addresses this by (i) **round‑robin rank to node mapping** to decorrelate node‑level loads, and (ii) **MPI RMA based intra‑node shared‑memory pooling** that dynamically redistributes column work across the ranks on a node during the chemistry phase. The redistribution follows a **Decouple - Compute - Remap** pipeline that packs per‑column KPP state into 1‑D arrays, executes chemistry from a node‑local work pool, then restores canonical ownership before transport, preserving ESMF/MAPL contracts and halo semantics. The approach yields large runtime reductions (often **~25–35%** in our experiments) without requiring a priori load information.  

#### Modified Scope
- `GeosCore/fullchem_mod.F90`: `Init_FullChem`, `Do_FullChem`, `Cleanup_FullChem`
- `GeosCore/gc_environment_mod.F90`: `GC_Init_Extra` (call update)

#### Key edits
- **`Init_FullChem`**
  - **Signature change:**  
    From: `SUBROUTINE Init_FullChem(Input_Opt, State_Chm, State_Diag, RC)`  
    To:   `SUBROUTINE Init_FullChem(Input_Opt, State_Chm, State_Diag, State_Grid, RC)`  
    (adds `State_Grid`).
  - **Intra-node MPI-RMA based Shared Memory setup:**
    - Create `shm_comm` via `MPI_Comm_split_type`.
    - Gather per-rank cell counts; compute **prefix offsets** and `NCELL_total`; broadcast offsets.
    - Allocate/attach MPI **shared-memory windows** for packed chemistry arrays (see list below).
    - Allocate `Idx_to_IJL(3, NCELL_MAX)` to record the mapping \( N \leftrightarrow (I,J,L) \).

- **`Do_FullChem`**
  - Replace the single triple loop over `(I,J,L)` with **Decouple - Compute - Remap**:
    - **Decouple (pack):** pack each rank’s per-column KPP state into its segment of node-wide 1-D arrays; 
    - **Compute:** enable **dynamic, node-local work sharing** via a shared atomic counter `next_cell_index` (`MPI_Fetch_and_op`); each rank claims the next cell and **integrates in place** in the shared arrays; optional `cell_status(N)` for tracing.
    - **Remap (unpack):** copy updated state for originally owned columns back to the rank’s grid; proceed to transport.

- **`Cleanup_FullChem`**
  - Free all `win_*` shared windows, deallocate bookkeeping arrays, and free `shm_comm`.

- **`GC_Init_Extra` (in `gc_environment_mod.F90`)**
  - Update the call to `Init_FullChem` to pass the new `State_Grid` argument.

#### Assumptions / constraints
- KPP is **cell-independent**; canonical ESMF/MAPL decomposition is restored before non-chemistry phases.
- Redistribution is **scoped to chemistry only**.
- Requires **MPI-3 shared memory / RMA** support.


### Expected changes
<!-- Please provide details on how this update will impact model output and include plots or tables as needed. -->
**Output: scientific results remain unchanged.** 
Example: C180, 2019-07-02:
```bash
Variable Xdim matches within the threshold of 1e-07.
Variable Ydim matches within the threshold of 1e-07.
Variable nf matches within the threshold of 1e-07.
Variable ncontact matches within the threshold of 1e-07.
Variable contacts matches within the threshold of 1e-07.
Variable anchor matches within the threshold of 1e-07.
Variable lons matches within the threshold of 1e-07.
Variable lats matches within the threshold of 1e-07.
Variable corner_lons matches within the threshold of 1e-07.
Variable corner_lats matches within the threshold of 1e-07.
Variable lev matches within the threshold of 1e-07.
Variable time matches within the threshold of 1e-07.
Variable Met_AIRDEN matches within the threshold of 1e-07.
Variable Met_BXHEIGHT matches within the threshold of 1e-07.
Variable Met_PBLH matches within the threshold of 1e-07.
Variable Met_PMID matches within the threshold of 1e-07.
Variable Met_PMIDDRY matches within the threshold of 1e-07.
Variable Met_PS1DRY matches within the threshold of 1e-07.
Variable Met_PS1WET matches within the threshold of 1e-07.
Variable Met_RH matches within the threshold of 1e-07.
Variable Met_T matches within the threshold of 1e-07.
Variable Met_TropHt matches within the threshold of 1e-07.
Variable Met_TropP matches within the threshold of 1e-07.
Variable PM25 matches within the threshold of 1e-07.
Variable SpeciesConcVV_CO matches within the threshold of 1e-07.
Variable SpeciesConcVV_NO matches within the threshold of 1e-07.
Variable SpeciesConcVV_NO2 matches within the threshold of 1e-07.
Variable SpeciesConcVV_O3 matches within the threshold of 1e-07.
Variable SpeciesConcVV_OH matches within the threshold of 1e-07.
All variables match within the given threshold.
```

**Performance: reduces wall time by ~24–41% across platforms and resolutions.**
> **Setup:** 576 ranks (cores) for C48–C180; 1152 ranks (cores) for C360; GCC + Intel MPI.

| Platform | Resolution | Original Time (s) | HB-SHM Time (s) | Time Saved (s) | Reduction | Speedup |
|---|---:|---:|---:|---:|---:|---:|
| Compute1 | C48  | 1825.2 | 1392.3 | 432.9  | 23.72% | 1.31× |
| Compute1 | C90  | 4840.6 | 3381.4 | 1459.2 | 30.15% | 1.43× |
| Compute1 | C180 | 16776.9 | 11540.4 | 5236.5 | 31.21% | 1.45× |
| AWS      | C48  | 1737.8 | 1275.7 | 462.0  | 26.59% | 1.36× |
| AWS      | C90  | 4925.0 | 3056.8 | 1868.2 | 37.93% | 1.61× |
| AWS      | C180 | 17352.5 | 10277.4 | 7075.2 | 40.77% | 1.69× |
| AWS      | C360 | 55499.7 | 34156.0 | 21343.7 | 38.46% | 1.62× |


**Memory overhead: additional memory is required for shared windows.** The table below >**Assumptions**: REAL/INTEGER = 8 B; Nspec=356, Nreact=1058, N(R)_state=20, N(I)_ctrl=20, N(R)_ctrl=20, N(I)_stat=20. Per-node values assume the noted node counts in our layouts.

- **C48:** 995,328 cells — **Total extra ≈ 11.12 GiB**; with 4 nodes ⇒ **≈ 2.78 GiB per node**  
- **C90:** 3,499,200 cells — **Total extra ≈ 39.08 GiB**; with 16 nodes ⇒ **≈ 2.44 GiB per node**  
- **C180:** 13,996,800 cells — **Total extra ≈ 156.32 GiB**; with 16 nodes ⇒ **≈ 9.77 GiB per node**  
- **C360:** 55,987,200 cells — **Total extra ≈ 625.29 GiB**; with 16 nodes ⇒ **≈ 39.08 GiB per node**  

Note: per-cell footprint ≈ 12 KB ≈ 1,500 doubles; constant across resolutions.

### Reference(s)
<!-- If this is a science update, please provide a literature citation. -->
Please refer to recently submitted work: [Exploring Load-Balancing Solutions for Earth-Scale Atmospheric Simulations](https://drive.google.com/file/d/1OGxPt1m9OVKSSxCQUO1YPXk-qUg5oZNs/view)


### Related Github Issue

<!-- Please link to the corresponding Github issue(s) here. If fixing a bug, there should be an issue describing it with steps to reproduce. -->
Our approach effectively addresses the load-imbalance problem discussed in [GCHP#52](https://github.com/geoschem/GCHP/issues/52).