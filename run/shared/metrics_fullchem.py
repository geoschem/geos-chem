#!/usr/bin/env python3
"""
metrics_fullchem.py: Prints key metrics (e.g. global mean OH, MCF
lifetime, and CH4 lifetime) for a GEOS-Chem full-chemistry simulation.
Requires Python3.

Calling sequence:
-----------------
./metrics.py
"""

# Imports
import os
import warnings
import numpy as np
import xarray as xr

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"] = "offscreen"

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def combine_dataset(file_list=None):
    """
    Wrapper for xarray.open_mfdataset, taking into account the
    extra arguments needed in xarray 0.15 and later.

    Args:
    -----
        file_list : list of str
    """

    # netCDF variables that we should skip reading
    # (These are from older versions of GCHP output)
    skip_these_vars = ["anchor",
                       "ncontact",
                       "orientation",
                       "contacts",
                       "cubed_sphere"]

    # Return a single Dataset containing data from all MeanOH files.
    # NOTE: Need to add combine="nested" and concat_dim="time"
    # for xarray 0.15 and higher!!!
    v = xr.__version__.split(".")
    if int(v[0]) == 0 and int(v[1]) >= 15:
        try:
            ds = xr.open_mfdataset(
                file_list,
                drop_variables=skip_these_vars,
                combine="nested",
                concat_dim="time"
            )
        except FileNotFoundError:
            msg = "Could not find one or more files in {}".format(file_list)
            raise FileNotFoundError(msg)
    else:
        try:
            ds = xr.open_mfdataset(
                file_list,
                drop_variables=skip_these_vars
            )
        except FileNotFoundError:
            msg = "Could not find one or more files in {}".format(file_list)
            raise FileNotFoundError(msg)

    return ds


def check_dataset(ds):
    """
    Makes sure the necessary variables are present in a dataset.

    Args:
    -----
        ds : xarray Dataset
    """

    # List of variables that we absolutely must have
    varlist = ["AirMassColumnFull", "MeanOHcolumnFull", "MCFlossInTrop"]

    # Error check inputs
    for v in varlist:
        if v not in ds.data_vars.keys():
            msg = "The '{}' diagnostic is not in this dataset!".format(v)
            raise ValueError(msg)

    return ds


def read_metrics_collection(data_dir="./OutputDir"):
    """
    Reads data from all "Metrics" collection netCDF files
    into a single xarray Dataset.

    Args:
    -----
        data_dir : str
            Directory containing data files.
            Default: "./OutputDir".

    Returns:
    --------
        ds : xarray Dataset
    """

    # Find a list of all MeanOH collection files in data_dir.
    # Walk through subdirectories of data_dir if they exist.
    file_list = []
    for root, dirs, files in os.walk(data_dir):
        if len(dirs) > 0:
            for d in dirs:
                for f in files:
                    if "GEOSChem.Metrics" in f:
                        file_list.append(os.path.join(root, d, f))
        else:
            for f in files:
                if "GEOSChem.Metrics" in f:
                    file_list.append(os.path.join(root, f))

    # Combine data into a single dataset
    # Exit if we do not have all necessary metrics variables
    ds = check_dataset(combine_dataset(file_list))

    return ds


def print_metrics(ds):
    """
    Prints the mass-weighted mean OH (full atmospheric column)
    from a GEOS-Chem simulation.

    Args:
    -----
        ds : xarray Dataset
    """

    # Conversion factors
    ten_to_minus_5 = np.float64(1.0e-5)
    s_per_yr = np.float64(86400.0) * np.float64(365.25)

    # Sum of air mass in the chemistry grid [molec]
    # NOTE: This is also the same as the total atmospheric burden of
    # methyl chloroform (MCF), which has a constant mixing ratio (=1).
    airmass = np.nansum(ds["AirMassColumnFull"].values)

    # Global Mean OH concentration [1e5 molec cm^-3]
    oh = (np.nansum(ds["MeanOHcolumnFull"].values) / airmass) * ten_to_minus_5

    # MCF lifetime [years]
    mcf = (airmass / np.nansum(ds["MCFlossInTrop"].values)) / s_per_yr

    # Get start and end time of run
    start = ds.attrs["simulation_start_date_and_time"]
    end = ds.attrs["simulation_end_date_and_time"]

    # Print results
    print("="*78)
    print("GEOS-Chem FULL-CHEMISTRY SIMULATION METRICS\n")
    print("Simulation start : {}".format(start))
    print("Simulation end   : {}".format(end))
    print("="*78, "\n")
    msg = "Mass-weighted Global Mean OH = {:.11f} x 10^5 molec cm-3\n".format(
        oh)
    print(msg)
    msg = "MCF lifetime in troposphere  = {:.4f} years\n ".format(mcf)
    print(msg)


def main():
    """
    Main program.  Calls all functions.
    """

    # Read
    ds = read_metrics_collection()

    # Compute global mean OHm across all
    print_metrics(ds)


if __name__ == "__main__":
    main()
