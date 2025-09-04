# Hierarchical Load Balancing with MPI Shared Memory for GCHP v14.5.2

## Summary
The KPP chemistry solver takes a variable number of sub-steps per column, creating **rank and node imbalance** when columns are pinned to MPI ranks. This repo implemented Hierarchical Balancing with MPI Shared Memory mitigates GCHP load imbalance by:
1) **Round-robin rank-node mapping** to decorrelate node loads, and  
2) **Intra-node shared-memory pooling** so ranks on a node dynamically share column work **during chemistry**.

The pipeline is **Decouple - Compute - Remap**: pack per-column KPP state into node-wide 1-D arrays; integrate from a node-local work pool via an atomic counter; restore canonical ownership before transport (ESMF/MAPL semantics preserved). In our tests we see **~25–35%** walltime reductions for 7-day simulation with different resolutions.

---

## Prerequisites (same as stock GCHP)
- A working GCHP environment per the [official requirements](https://gchp.readthedocs.io/en/stable/getting-started/requirements.html).

---

## Integrate into your GCHP 

### Point the submodule to this branch
From your **GCHP repo root**:
```bash
# Make sure the submodule exists locally
git submodule update --init --recursive

# Point the geos-chem submodule at this repo
git submodule set-url geos-chem https://github.com/Jordan-Sun/geos-chem.git
git submodule sync
git submodule update --init --remote

# Check out the feature branch inside the submodule
cd src/GCHP_GridComp/GEOSChem_GridComp/geos-chem
git fetch --all
git switch shared_memory
cd -    # back to GCHP root

---

## Build (unchanged)
Build GCHP as usual:
```bash
mkdir -p build && cd build
cmake .. -DRUNDIR="/path/to/your/run"
make -j
```
Use standard compiler/MPI stack (e.g., Intel MPI or Open MPI).

---

## Run 
Map ranks to nodes in a round-robin fashion:
- **Open MPI:** `--map-by ppr:1:node`
- **Intel MPI:** `-ppn 1`
### Open MPI

```bash
time mpirun -np ${N_CORES} --map-by ppr:1:node ./gchp
```

### Intel MPI
```bash
time mpirun -np ${N_CORES} -ppn 1 ./gchp
```