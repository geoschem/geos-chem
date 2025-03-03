#!/usr/bin/env python3
"""
Prints key metrics (e.g. global mean OH, MCF lifetime, and CH4 lifetimes)
for a GEOS-Chem full-chemistry simulation or methane simulation.
Requires Python3.

Calling sequence:
-----------------
./metrics.py
"""
# =====================================================================
# %%% IMPORTS ETC. %%%
# =====================================================================
import os
import sys
import warnings
import numpy as np
import xarray as xr
import yaml

# Tell matplotlib not to look for an X-window
os.environ["QT_QPA_PLATFORM"] = "offscreen"

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# =====================================================================
# %%% GLOBAL VARIABLES %%%
# =====================================================================

# Molecular weights (get from species database)
spec_db = yaml.load(open("./species_database.yml"), Loader=yaml.FullLoader)
mw_ch4_kg = spec_db["CH4"]["MW_g"] * 1.0e-3
mw_oh_kg = spec_db["OH"]["MW_g"] * 1.0e-3
mw_air_kg = 28.9644e-3                        # Same value as in GEOS-Chem

# Physical constants and conversion factors
avogadro = 6.022140857e+23                    # Same value as in GEOS-Chem
m3_to_cm3 = 1.0e6
s_per_yr = np.float64(86400.0) * np.float64(365.25)
ten_to_minus_5 = np.float64(1.0e-5)
kg_to_m_ch4 = mw_ch4_kg / avogadro

# =====================================================================
# %%% METHODS %%%
# =====================================================================

def combine_dataset(file_list=None):
    """
    Wrapper for xarray.open_mfdataset, taking into account the
    extra arguments needed in xarray 0.15 and later.

    Args:
    -----
        file_list : list of str

    Returns:
    --------
        ds : xarray Dataset
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


def ch4_or_fullchem(ds):
    """
    Determines if a Dataset contains variables for computing
    metrics from a CH4 simulation or a fullchem simulation.

    Args:
    -----
        ds : xarray Dataset

    Returns:
    --------
        is_ch4_sim : bool
    """

    # CH4 and fullchem simulations have these variables
    common_vars = [
        "AirMassColumnFull",
        "LossOHbyCH4columnTrop",
        "LossOHbyMCFcolumnTrop",
        "OHwgtByAirMassColumnFull",
    ]

    # CH4 simulations also have these variables
    extra_vars_ch4 = [
        "CH4emission",
        "CH4massColumnFull",
        "CH4massColumnTrop",
    ]

    # Keep a count
    count = 0

    # Look for the common variables in the dataset
    for v in common_vars:
        if v in ds.data_vars.keys():
            count += 1

    # Look for the CH4-only variables in the dataset
    for v in extra_vars_ch4:
        if v in ds.data_vars.keys():
            count += 1

    if count == len(common_vars):
        is_ch4_sim = False
    elif count == len(common_vars) + len(extra_vars_ch4):
        is_ch4_sim = True
    else:
        msg = "The '{}' diagnostic is not in this dataset!".format(v)
        raise ValueError(msg)

    return is_ch4_sim


def read_metrics_collection(data_dir):
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
                    if ".Metrics." in f:
                        file_list.append(os.path.join(root, d, f))
        else:
            for f in files:
                if ".Metrics." in f:
                    file_list.append(os.path.join(root, f))

    # Combine data into a single dataset
    # Exit if we do not have all necessary metrics variables
    ds = combine_dataset(file_list)

    return ds


def total_airmass(ds):
    """
    Computes the total airmass (in both kg and molec).

    Args:
    -----
        ds : xarray Dataset

    Returns:
    --------
        sum_airmass_kg, sum_airmass_m: numpy float64
    """
    sum_airmass_kg = np.nansum(ds["AirMassColumnFull"].values)
    sum_airmass_m = sum_airmass_kg * (avogadro / mw_air_kg)

    return sum_airmass_kg, sum_airmass_m


def global_mean_oh(sum_airmass_kg, ds):
    """
    Computes the global mean OH concentration (1e5 molec cm-3)

    Args:
    -----
        sum_airmass_kg : numpy float64
        ds : xarray Dataset

    Returns:
    --------
        sum_mean_oh : numpy float64
    """
    # Divide out total airmass to get total mean OH concentration [kg m-3]
    # Then convert mean OH from [kg m-3] to [1e5 molec cm-3]
    sum_mean_oh = np.nansum(ds["OHwgtByAirMassColumnFull"].values)
    sum_mean_oh = (sum_mean_oh / sum_airmass_kg)
    sum_mean_oh *= (avogadro / (mw_oh_kg * m3_to_cm3)) * ten_to_minus_5

    return sum_mean_oh


def lifetimes_wrt_oh(sum_airmass_m, ds):
    """
    Computes the lifetimes (in years) of CH4 and CH3CCl3 (aka MCF)
    against tropospheric OH.

    Args:
    -----
        sum_airmass_m : numpy float64
        ds : xarray Dataset

    Returns:
    --------
        ch4_life_wrt_oh, mcf_life_wrt_oh : numpy float64
    """
    # Loss of OH by CH4+OH and MCF+OH reactions [molec]
    oh_loss_by_ch4 = np.nansum(ds["LossOHbyCH4columnTrop"].values)
    oh_loss_by_mcf = np.nansum(ds["LossOHbyMCFcolumnTrop"].values)

    # CH4 and MCF lifetimes against OH [years]
    ch4_life_wrt_oh = (sum_airmass_m / oh_loss_by_ch4) / s_per_yr
    mcf_life_wrt_oh = (sum_airmass_m / oh_loss_by_mcf) / s_per_yr

    return ch4_life_wrt_oh, mcf_life_wrt_oh


def overall_ch4_lifetimes(ds):
    """
    Computes the overall lifetimes (in years) of CH4 (accounting for
    emissions) in the full-atmosphere and in the troposphere.

    Args:
    -----
        ds : xarray Dataset

    Returns:
    --------
        ch4_life_full, ch4_life_trop : numpy float64
    """
    # Sum of CH4 emissions [kg/s]
    ch4_emis = np.nansum(ds["CH4emission"])

    # CH4 mass in full-atm column and trop-only column [kg]
    ch4_mass_full = np.nansum(ds["CH4massColumnFull"].values)
    ch4_mass_trop = np.nansum(ds["CH4massColumnTrop"].values)

    # CH4 overall lifetimes [years] in full-atmosphere and troposphere
    ch4_life_full = (ch4_mass_full / ch4_emis) / s_per_yr
    ch4_life_trop = (ch4_mass_trop / ch4_emis) / s_per_yr

    return ch4_life_full, ch4_life_trop


def get_start_and_end_dates(ds):
    """
    Gets the start and end dates of a GEOS-Chem Classic or
    a GCHP simulation.

    Args:
    -----
        ds : xarray Dataset

    Returns:
    --------
        start, end : str
     """
    # Get start and end date of simulation for GCHP or GEOS-Chem Classic
    if "nf" in ds.dims:
        with open("./CAP.rc", "r") as cap_file:
            for line in cap_file:
                if "BEG_DATE:" in line:
                    substrs = (line.rstrip()).split()
                    start = substrs[1] + " " + substrs[2] + "z"
                elif "END_DATE:" in line:
                    substrs = (line.rstrip()).split()
                    end = substrs[1] + " " + substrs[2] + "z"
    else:
        start = ds.attrs["simulation_start_date_and_time"]   # GC-Classic
        end = ds.attrs["simulation_end_date_and_time"]

    return start, end


def print_metrics(ds, is_ch4_sim=False):
    """
    Prints the mass-weighted mean OH (full atmospheric column)
    from a GEOS-Chem simulation.

    Args:
    -----
        ds : xarray Dataset
        is_ch4_sim : bool
    """

    # Get total airmasses ([kg] and [molec])
    sum_airmass_kg, sum_airmass_m = total_airmass(ds)

    # Get mean OH [1e-5 molec cm-3]
    mean_oh = global_mean_oh(sum_airmass_kg, ds)

    # Get lifetimes of CH4 and MCF against tropospheric OH [years]
    ch4_life_wrt_oh, mcf_life_wrt_oh = lifetimes_wrt_oh(sum_airmass_m, ds)

    # Get overall lifetime of CH4 in full-atm and in trop [years[
    if is_ch4_sim:
        ch4_life_full, ch4_life_trop = overall_ch4_lifetimes(ds)

    # Get start and end dates of the simulation
    start, end = get_start_and_end_dates(ds)

    # Print results
    print("="*78)
    if is_ch4_sim:
        print("GEOS-Chem METHANE SIMULATION METRICS\n")
    else:
        print("GEOS-Chem FULL-CHEMISTRY SIMULATION METRICS\n")
    print("Simulation start : {}".format(start))
    print("Simulation end   : {}".format(end))
    print("="*78, "\n")
    msg = "Mass-weighted mean OH concentration    = {:.11f} ".format(mean_oh)
    msg += "x 10^5 molec cm-3\n"
    print(msg)
    msg = "CH3CCl3 lifetime w/r/t tropospheric OH = {:.4f} years\n".format(
        mcf_life_wrt_oh)
    print(msg)
    msg = "CH4 lifetime w/r/t tropospheric OH     = {:.4f} years\n".format(
        ch4_life_wrt_oh)
    print(msg)

    # These metrics are only valid for the CH4 simulation
    if is_ch4_sim:
        msg = "CH4 total lifetime (full atmosphere)   = {:.4f} years\n".format(
            ch4_life_full)
        print(msg)
        msg = "CH4 total lifetime (troposphere only)  = {:.4f} years\n".format(
            ch4_life_trop)
        print(msg)


def main():
    """
    Main program.  Calls all functions.
    """
    # Get the data directory from the passed arguments
    n_args = len(sys.argv)
    if n_args == 1:
        data_dir = "./OutputDir"
    elif n_args == 2:
        data_dir = sys.argv[1]
    else:
        raise ValueError("Usage: metrics.py DIRECTORY-NAME")

    # Combine all metrics output into a single Dataset
    ds = read_metrics_collection(data_dir)

    # Determine if this dataset is for a CH4 simulation
    # Otherwise, assume it is a full-chemistry simulation
    is_ch4_sim = ch4_or_fullchem(ds)

    # Print metrics information
    print_metrics(ds, is_ch4_sim)


if __name__ == "__main__":
    main()
