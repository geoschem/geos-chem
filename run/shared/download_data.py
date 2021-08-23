#!/usr/bin/env python3

"""
Description:
------------
This Python script (assumes Python3) reads a GEOS-Chem or
HEMCO-standalone log file containing dry-run output and does
the following:

    (1) Creates a list of unique files that are required for the
        GEOS-Chem or HEMCO-standalone simulation;

    (2) Creates a bash script to download missing files from either
        the ComputeCanada server (default) or the AWS s3://gcgrid
        bucket;

    (3) Executes the bash script to download the necessary data;

    (4) Removes the bash script upon successful download.


Remarks:
--------
    (1) This script only requires the "os", "sys", and "subprocess"
        packages, which are core Python.  Therefore, this script can
        be shipped with GEOS-Chem run directories.  It only requires
        Python 3 and not a full Anaconda/Miniconda environment (but
        you can run in an Anaconda environment if you have one).

    (2) Jiawei Zhuang found that it is much faster to issue aws s3 cp
        commands from a bash script than a Python script.  Therefore,
        in this routine we create a bash script with all of the
        download commands that will be executed by the main routine.
"""

# Imports
import os
import sys
import subprocess

# Exit with error if we are not using Python3
assert sys.version_info.major >= 3, \
"ERROR: Python 3 is required to run download_data.py!"

# Define global variables
INPUT_GEOS_FILE = "./input.geos"
DATA_DOWNLOAD_SCRIPT = "./auto_generated_download_script.sh"


def extract_pathnames_from_log(dryrun_log):
    """
    Returns a list of pathnames from a GEOS-Chem log file.

    Args:
    -----
        dryrun_log : str
            GEOS-Chem or HEMCO-standalone log file with dry-run output.

    Returns:
    --------
        paths : dict
            paths["comments"]: Dry-run comment lines.
            paths["found"] : List of file paths found on disk.
            paths["missing"]: List of file paths that are missing.
            paths["local_prefix"]: Local data directory root.

    Author:
    -------
        Jiawei Zhuang (jiaweizhuang@g.harvard.edu)
        Modified by Bob Yantosca (yantosca@seas.harvard.edu)
    """

    # Initialization
    comments = ["!"*79,
                "!!! LIST OF (UNIQUE) FILES REQUIRED FOR THE SIMULATION"]
    data_found = set()
    data_missing = set()

    # Open file (or die with error)
    try:
        f = open(dryrun_log, "r")
    except FileNotFoundError:
        raise FileNotFoundError("Could not find file {}".format(dryrun_log))

    # Read data from the file line by line.
    # Add file paths to the data_list set.
    line = f.readline()
    while line:

        # Convert line to uppercase for string match
        upcaseline = line.upper()

        # Search for data paths that have been found
        if (": OPENING" in upcaseline) or (": READING" in upcaseline):
            data_found.add(line.split()[-1])

        # Search for data paths that are missing
        elif "FILE NOT FOUND" in upcaseline:
            data_missing.add(line.split()[-1])

        # Search for certain dry-run comment strings
        # (and make sure to prevent duplicates)
        elif ("!!! STA" in upcaseline) or ("!!! END" in upcaseline) or \
             ("!!! SIM" in upcaseline) or ("!!! MET" in upcaseline) or \
             ("!!! GRI" in upcaseline):
            if line.rstrip() not in comments:
                comments.append(line.rstrip())

        else:
            pass

        # Read next line
        line = f.readline()

    # Add another line to the comment list
    comments.append("!"*79)

    # Convert sets to lists and sort in alphabetical order
    found = sorted(list(data_found))
    missing = sorted(list(data_missing))

    # Find the local data directory prefix (path to ExtData)
    local_prefix = ""
    for path in found + missing:
        if "ExtData" in path:
            index = path.find("ExtData")
            local_prefix = path[:index]
            break

    # Exit if the local path does not contain ExtData
    if len(local_prefix) == 0:
        msg = "Could not locate the ExtData folder in your local disk space!"
        raise ValueError(msg)

    # Close file and return
    # The "sorted" command will return unique values
    f.close()
    return {
        "comments": comments,
        "found": found,
        "missing": missing,
        "local_prefix": local_prefix
    }


def get_run_info():
    """
    Searches through the input.geos file for GEOS-Chem run parameters.

    Returns:
    -------
        run_info : dict of str
            Contains the GEOS-Chem run parameters: start_date,
            start_time, end_date, end_time, met, grid, and sim.
    """
    run_info = {}
    run_info["nest"] = ""
    run_info["tomas15"] = False
    run_info["tomas40"] = False

    try:
        with open(INPUT_GEOS_FILE, "r") as f:
            for line in f:
                if "Start YYYYMMDD" in line:
                    substr = line.split(":")[1]
                    run_info["start_date"] = (substr.split(" ")[1]).strip()
                    run_info["start_time"] = (substr.split(" ")[2]).strip()
                elif "End   YYYYMMDD" in line:
                    substr = line.split(":")[1]
                    run_info["end_date"] = (substr.split(" ")[1]).strip()
                    run_info["end_time"] = (substr.split(" ")[2]).strip()
                elif "Met field" in line:
                    run_info["met"] = (line.split(":")[1]).strip()
                elif "Simulation name" in line:
                    run_info["sim"] = (line.split(":")[1]).strip()
                elif "Grid resolution" in line:
                    grid = (line.split(":")[1]).strip()

                    # Adjust grid string to match file names
                    if "4.0x5.0" in grid:
                        run_info["grid"] = "4x5"
                    elif "2.0x2.5" in grid:
                        run_info["grid"] = "2x25"
                    elif "0.5x0.625" in grid:
                        run_info["grid"] = "05x0625"
                    elif "0.25x0.3125" in grid:
                        run_info["grid"] = "025x03125"
                elif "Longitude" in line:
                    if "-130.0" in line or "-140.0" in line:
                        run_info["nest"] = "na"
                        break
                    elif "60.0" in line or "70.0" in line:
                        run_info["nest"] = "as"
                elif "NK15" in line:
                    run_info["tomas15"] = True
                elif "NK40" in line:
                    run_info["tomas15"] = False
                    run_info["tomas40"] = True
            f.close()
    except FileNotFoundError:
        raise FileNotFoundError("Could not open {}".format(INPUT_GEOS_FILE))

    return run_info


def expand_restart_file_names(paths, run_info):
    """
    Tests if the GEOS-Chem restart file is a symbolic link to
    ExtData.  If so, will append the link to the remote file
    to the line in which the restart file name is found.

    Args:
    ----
        paths : dict
            Output of function extract_pathnames_from_log.

        run_info : dict
            Output of function get_run_info.
    """
    prefix = ""

    # Get the prefix to ExtData
    for path in paths["found"] + paths["missing"]:
        if "ExtData" in path:
            index = path.find("ExtData")+8

            # Fullchem simulations use restarts from benchmarks (GC_13.0.0)
            # but TOMAS simulations use restarts from v2020-02
            if "fullchem" in run_info["sim"]:
                if run_info["tomas15"] is True or run_info["tomas40"] is True:
                    prefix = path[0:index] + "GEOSCHEM_RESTARTS/v2020-02/"
                else:
                    prefix = path[0:index] + "GEOSCHEM_RESTARTS/GC_13.0.0/"
                break

            # TransportTracers sims use restarts from benchmarks (GC_13.0.0)
            elif "TransportTracers" in run_info["sim"]:
                prefix = path[0:index] + "GEOSCHEM_RESTARTS/GC_13.0.0/"
                break

            elif "aerosol" in run_info["sim"]:
                prefix = path[0:index] + "GEOSCHEM_RESTARTS/GC_13.0.0/"
                break

            # All other simulations use restarts from v2020-02
            else:
                prefix = path[0:index] + "GEOSCHEM_RESTARTS/v2020-02/"

    # Search for the restart file name in the found files
    new_list = []

    # Suffix string (takes into account nested grids)
    if run_info["nest"] == "":
        suffix = "{}.nc".format(run_info["sim"])
    else:
        suffix = "{}_{}.nc".format(run_info["sim"], run_info["nest"])

    for path in paths["found"]:
        if "GEOSChem.Restart" in path:
            if "fullchem" in run_info["sim"]:
                if run_info["tomas15"] is True:
                    realpath = prefix + "initial_GEOSChem_rst.4x5_TOMAS15.nc"
                elif run_info["tomas40"] is True:
                    realpath = prefix + "initial_GEOSChem_rst.4x5_TOMAS40.nc"
                else:
                    realpath = prefix + \
                        "GEOSChem.Restart.fullchem.20190701_0000z.nc4"
            elif "TransportTracers" in run_info["sim"]:
                realpath = prefix + \
                    "GEOSChem.Restart.TransportTracers.20190101_0000z.nc4"
            elif "aerosol" in run_info["sim"]:
                    realpath = prefix + \
                        "GEOSChem.Restart.fullchem.20190101_0000z.nc4"
            else:
                realpath = prefix + "initial_GEOSChem_rst.2x25_" + suffix
            path = path + " --> " + realpath
        new_list.append(path)
    paths["found"] = sorted(new_list)

    # Search for the restart file name in the missing files
    new_list = []
    for path in paths["missing"]:
        if "GEOSChem.Restart" in path:
            if "fullchem" in run_info["sim"]:
                if run_info["tomas15"] is True:
                    realpath = prefix + "initial_GEOSChem_rst.4x5_TOMAS15.nc"
                elif run_info["tomas40"] is True:
                    realpath = prefix + "initial_GEOSChem_rst.4x5_TOMAS40.nc"
                else:
                    realpath = prefix + \
                        "GEOSChem.Restart.fullchem.20190701_0000z.nc4"
            elif "TransportTracers" in run_info["sim"]:
                realpath = prefix + \
                    "GEOSChem.Restart.TransportTracers.20190101_0000z.nc4"
            elif "aerosol" in run_info["sim"]:
                    realpath = prefix + \
                        "GEOSChem.Restart.fullchem.20190101_0000z.nc4"
            else:
                realpath = prefix + "initial_GEOSChem_rst.2x25_" + suffix
            path = path + " --> " + realpath
        new_list.append(path)
    paths["missing"] = sorted(new_list)

    # Return the updated data paths
    return paths


def write_unique_paths(paths, unique_log):
    """
    Writes unique data paths from dry-run output to a file.

    Args:
    -----
        paths : dict
            Output of function extract_pathnames_from_log.

        unique_log : str
            Log file that will hold unique data paths.
    """
    combined_paths = paths["found"] + paths["missing"]
    combined_paths.sort()

    try:
        with open(unique_log, "w") as f:
            for comment in paths["comments"]:
                print(comment, file=f)
            for path in combined_paths:
                print(path, file=f)
            for comment in paths["comments"]:
                print(comment, file=f)
        f.close()
        print("Log with unique file paths written to: {}".format(unique_log))
    except FileNotFoundError:
        raise FileNotFoundError("Could not write {}".format(unique_log))


def create_download_script(paths, args):
    """
    Creates a data download script to obtain missing files
    from the ComputeCanada data archive (default), or the
    GEOS-Chem s3://gcgrid bucket on the AWS cloud,

    Args:
    -----
        paths : dict
            Output of function extract_pathnames_from_log.

        args : dict
            Arguments 
    """

    # Get the path to the data on the remote server
    remote_root = args["remote_root"]
    
    # Define variables to create data download commands
    # for either ComputeCanada or AWS  
    if "AWS" in args["remote"]:
        cmd_prefix = "aws s3 cp --request-payer=requester "
        quote = ""
    if ("WUSTL" in args["remote"]) or ("CC" in args["remote"]):
        cmd_prefix = 'wget -r -np -nH -R "*.html" -N -P ' + \
                     paths["local_prefix"] + " "
        quote = '"'
    if "UR" in args["remote"]:
        cmd_prefix = 'wget -r -np -nH -R "*.html" -N --cut-dirs=2 -P ' + \
                     paths["local_prefix"] + " "
        quote = '"'

    # Create the data download script
    with open(DATA_DOWNLOAD_SCRIPT, "w") as f:

        # Write shebang line to script
        print("#!/bin/bash\n", file=f)
        print("# This script was generated by download_data.py\n", file=f)

        # Write download commands for only the missing data files
        for path in paths["missing"]:

            if "-->" in path:

                # ------------------------------------------------------
                # Edge case: Linked restart files
                # ------------------------------------------------------

                # First copy the restart file to local ExtData
                remote_rst = (path.split("-->")[1]).strip()
                local_rst = (path.split("-->")[0]).strip()
                index1 = remote_rst.find("initial")
                index2 = remote_rst.find("ExtData") + 7
                prefix = local_rst
                extdata = remote_rst[:index2]
                remote_rst = remote_root + remote_rst[index2:]
                cmd = cmd_prefix + quote + remote_rst + quote
                if "AWS" in args["remote"]:
                    cmd += " " + prefix
                print(cmd, file=f)
                print(file=f)

                # If the file does not exist in the run directory,
                # then copy it from the restart folder.
                # This only has to be done from WUSTL, CC, or UR,
                # since AWS will download directly to the rundir.
                if "AWS" not in args["remote"]:
                    if not os.path.exists(local_rst):
                        index3 = remote_rst.find("GEOSCHEM_RESTARTS")
                        rst = os.path.join(extdata, remote_rst[index3:])
                        cmd = "cp -f " + rst + " " + local_rst
                        print(cmd, file=f)
                        print(file=f)

            elif "gmi.clim.IPMN.geos5.2x25.nc" in path:

                # ------------------------------------------------------
                # Edge case: GMI IPMN file is really the PMN file
                # ------------------------------------------------------

                # Download the PMN file
                index = path.find("ExtData") + 7
                local_dir = os.path.dirname(path)
                remote_path = remote_root + path[index:]
                remote_path = remote_path.replace("IPMN", "PMN")
                cmd = cmd_prefix + quote + remote_path + quote
                if "AWS" in args["remote"]:
                    cmd += " " + local_dir + "/"
                print(cmd, file=f)

                # Rename it to IPMN
                cmd = "mv " + local_dir + "/gmi.clim.PMN.geos5.2x25.nc " + \
                      local_dir + "/gmi.clim.IPMN.geos5.2x25.nc"
                print(cmd, file=f)

            elif "gmi.clim.NPMN.geos5.2x25.nc" in path:

                # ------------------------------------------------------
                # Edge case: GMI NPMN file is really the PMN file
                # ------------------------------------------------------

                # Download the PMN file
                index = path.find("ExtData") + 7
                local_dir = os.path.dirname(path)
                remote_path = remote_root + path[index:]
                remote_path = remote_path.replace("NPMN", "PMN")
                cmd = cmd_prefix + quote + remote_path + quote
                if "AWS" in args["remote"]:
                    cmd += " " + local_dir + "/"
                print(cmd, file=f)

                # Rename it to NPMN
                cmd = "mv " + local_dir + "/gmi.clim.PMN.geos5.2x25.nc " + \
                      local_dir + "/gmi.clim.NPMN.geos5.2x25.nc"
                print(cmd, file=f)
                print(file=f)

            elif "gmi.clim.RIPA.geos5.2x25.nc" in path:

                # ------------------------------------------------------
                # Edge case: GMI RIPA file is really the RIP file
                # ------------------------------------------------------

                # Download the RIP file
                index = path.find("ExtData")+7
                local_dir = os.path.dirname(path)
                remote_path = remote_root + path[index:]
                remote_path = remote_path.replace("RIPA", "RIP")
                cmd = cmd_prefix + quote + remote_path + quote
                if "AWS" in args["remote"]:
                    cmd += " " + local_dir + "/"
                print(cmd, file=f)

                # Rename it to NPMN
                cmd = "mv " + local_dir + "/gmi.clim.RIP.geos5.2x25.nc " + \
                      local_dir + "/gmi.clim.RIPA.geos5.2x25.nc"
                print(cmd, file=f)
                print(file=f)

            elif "gmi.clim.RIPB.geos5.2x25.nc" in path:

                # ------------------------------------------------------
                # Edge case: GMI RIPB file is really the RIP file
                # ------------------------------------------------------

                # Download the RIP file
                index = path.find("ExtData")+7
                local_dir = os.path.dirname(path)
                remote_path = remote_root + path[index:]
                remote_path = remote_path.replace("RIPB", "RIP")
                cmd = cmd_prefix + quote + remote_path + quote
                if "AWS" in args["remote"]:
                    cmd += " " + local_dir + "/"
                print(cmd, file=f)

                # Rename it to RIPB
                cmd = "mv " + local_dir + "/gmi.clim.RIP.geos5.2x25.nc " + \
                      local_dir + "/gmi.clim.RIPB.geos5.2x25.nc"
                print(cmd, file=f)
                print(file=f)

            elif "gmi.clim.RIPD.geos5.2x25.nc" in path:

                # ------------------------------------------------------
                # Edge case: GMI RIPD file is really the RIP file
                # ------------------------------------------------------

                # Download the RIP file
                index = path.find("ExtData")+7
                local_dir = os.path.dirname(path)
                remote_path = remote_root + path[index:]
                remote_path = remote_path.replace("RIPD", "RIP")
                cmd = cmd_prefix + quote + remote_path + quote
                if "AWS" in args["remote"]:
                    cmd += " " + local_dir + "/"
                print(cmd, file=f)

                # Rename it to RIPD
                cmd = "mv " + local_dir + "/gmi.clim.RIP.geos5.2x25.nc " + \
                      local_dir + "/gmi.clim.RIPD.geos5.2x25.nc"
                print(cmd, file=f)
                print(file=f)

            elif "ExtData" in path:

                # ------------------------------------------------------
                # All other files in ExtData
                # ------------------------------------------------------
                index = path.find("ExtData") + 7
                local_dir = os.path.dirname(path)
                remote_path = remote_root + path[index:]
                cmd = cmd_prefix + quote + remote_path + quote
                if "AWS" in args["remote"]:
                    cmd += " " + local_dir + "/"
                print(cmd, file=f)
                print(file=f)

        # Kludge: Create a ExtData/CHEM_INPUTS folder if it
        # does not exist. This will prevent abnormal exits.
        chem_inputs_dir = paths["local_prefix"] + 'ExtData/CHEM_INPUTS'
        cmd = "if [[ ! -d {} ]]; then mkdir {}; fi".format(
            chem_inputs_dir, chem_inputs_dir)
        print(cmd, file=f)
        print(file=f)

        # Close file and make it executable
        f.close()
        os.chmod(DATA_DOWNLOAD_SCRIPT, 0o755)


def download_the_data(args):
    """
    Downloads GEOS-Chem data files from the ComputeCanada server
    or the AWS s3://gcgrid bucket.

    Args:
    -----
        args : dict
            Output of runction parse_args.
    """

    # Get information about the run
    run_info = get_run_info()

    # Get a unique list of data paths, both found and missing:
    # Expand the data paths to include links to restart files
    paths = extract_pathnames_from_log(args["dryrun_log"])
    paths = expand_restart_file_names(paths, run_info)

    # Write a list of unique file paths
    write_unique_paths(paths,
                       args["dryrun_log"] + ".unique")

    # Exit without downloading if skip-download lag was specified
    if args["skip_download"]:
        return

    # Print a message
    if len(args["remote"]) > 0:
        print("Downloading data from " + args["remote"])

    # Create script to download missing files from AWS S3
    create_download_script(paths, args)

    # Run the data download script and return the status
    # Remove the file afterwards
    status = subprocess.call(DATA_DOWNLOAD_SCRIPT)
    os.remove(DATA_DOWNLOAD_SCRIPT)

    # Raise an exception if the data was not successfully downloaded
    if status != 0:
        err_msg = "Error downloading data from " + args["remote"]
        raise Exception(err_msg)


def parse_args():
    """
    Parses the arguments passed to the main program.

    Args:
    -----
        argv : list
            List of arguments passed from main().

    Returns:
    --------
        args : dict
            args["dryrun_log"]: Log file with dry-run output.
            args["remote"]: Name of the remote server.
            args["remote_root"]: Path to ther remote data
            args["skip-download"]: Skip downloading and only write
               out the log with unique file names.
    """
    dryrun_log = ""
    remote = ""
    remote_root = ""
    skip_download = False

    for i in range(1, len(sys.argv)):
        if "AWS" in sys.argv[i].upper():
            remote = "AWS"
            remote_root = "s3://gcgrid"
        elif "CC" in sys.argv[i].upper():
            remote = "CC"
            remote_root = "http://geoschemdata.computecanada.ca/ExtData"
        elif "WUSTL" in sys.argv[i].upper():
            remote = "WUSTL"
            remote_root = "http://geoschemdata.wustl.edu/ExtData"
        elif "UR" in sys.argv[i].upper():
            remote = "UR"
            remote_root = "http://atmos.earth.rochester.edu/input/gc/ExtData/"
        elif "SKIP" in sys.argv[i].upper():
            skip_download = True
        else:
            dryrun_log = sys.argv[i]

    if len(dryrun_log) == 0:
        raise ValueError("Need to specify the log file with dryrun output!")

    return {
        "dryrun_log": dryrun_log,
        "remote" : remote,
        "remote_root": remote_root,
        "skip_download" : skip_download
    }


def main():
    """
    Main program.  Gets command-line arguments and calls function
    download_the_data to initiate a data-downloading process.

    Calling sequence:
    -----------------
        ./download_data.py log -aws            # from AWS
        ./download_data.py log -wustl          # from Washington University in St. Louis
        ./download_data.py log -cc             # from ComputeCanada
        ./download_data.py log -ur             # from U. Rochester
        ./download_data.py log -skip-download  # Print unique log & exit
    """
    download_the_data(parse_args())


if __name__ == "__main__":
    main()
