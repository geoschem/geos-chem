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
import yaml

# Exit with error if we are not using Python3
assert sys.version_info.major >= 3, \
"ERROR: Python 3 is required to run download_data.py!"

# Define global variables
INPUT_GEOS_FILE = "./input.geos"
DATA_DOWNLOAD_SCRIPT = "./auto_generated_download_script.sh"


def extract_pathnames_from_log(args):
    """
    Returns a list of pathnames from a GEOS-Chem log file.

    Args:
    -----
    args : dict
        Contains output from function parse_args.

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
    dryrun_log = args["dryrun_log"]

    # Open file (or die with error)
    try:
        f = open(dryrun_log, "r")
    except FileNotFoundError:
        msg = "Could not find file " + dryrun_log
        raise FileNotFoundError(msg)

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

    paths = {
        "comments": comments,
        "found": found,
        "missing": missing,
        "local_prefix": local_prefix
    }
    return paths


def get_run_info():
    """
    Searches through the input.geos file for GEOS-Chem run parameters.

    Returns:
    -------
    run_info : dict
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
        msg = "Could not open " + INPUT_GEOS_FILE
        raise FileNotFoundError(msg)

    return run_info


def expand_restart_file_names(paths, args, run_info):
    """
    Tests if the GEOS-Chem restart file is a symbolic link to
    ExtData.  If so, will append the link to the remote file
    to the line in which the restart file name is found.

    Args:
    ----
    paths : dict
        Contains output from function extract_pathnames_from_log.
    args : dict
        Contains output from function parse_args.
    run_info : dict
        Contains output from function get_run_info.
    """
    remote_rst = ""
    rst = args["config"]["restarts"]

    # ------------------------------------------------------------------
    # Get the full name of the restart file in ExtData
    # ------------------------------------------------------------------
    for path in paths["found"] + paths["missing"]:
        if "ExtData" in path:
            index = path.find("ExtData")+8
            root = path[0:index] + rst["root"]

            if "aerosol" in run_info["sim"]:
                remote_rst = root + rst["aerosol"]["remote"]

            elif "fullchem" in run_info["sim"]:
                if run_info["tomas15"] is True:
                    remote_rst = root + rst["tomas15"]["remote"]
                elif run_info["tomas40"] is True:
                    remote_rst = root + rst["tomas40"]["remote"]
                else:
                    remote_rst = root + rst["fullchem"]["remote"]

            elif "TransportTracers" in run_info["sim"]:
                remote_rst = root + rst["transporttracers"]["remote"]

            else:
                remote_rst = root + rst["other"]["remote"]

    # Append a suffix string (e.g. for nested grids) if necessary
    if run_info["nest"] == "":
        suffix = "{}.nc".format(run_info["sim"])
    else:
        suffix = "{}_{}.nc".format(run_info["sim"], run_info["nest"])
    remote_rst = remote_rst.replace("@SUFFIX@", suffix)

    # ------------------------------------------------------------------
    # Search for the restart file name in the found files
    # ------------------------------------------------------------------
    new_list = []
    for path in paths["found"]:
        if "GEOSChem.Restart" in path:
            path = path + " --> " + remote_rst

        new_list.append(path)
    paths["found"] = sorted(new_list)

    # ------------------------------------------------------------------
    # Search for the restart file name in the missing files
    # ------------------------------------------------------------------
    new_list = []
    for path in paths["missing"]:
        if "GEOSChem.Restart" in path:
            path = path + " --> " + remote_rst
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
            Contains output from function extract_pathnames_from_log.

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
        Contains output from function extract_pathnames_from_log.
    args : dict
        Contains output from function parse_args.
    """

    # Extract mirror parameters
    mirror_name = args["mirror"]
    mirror = args["config"]["mirrors"][mirror_name]
    is_s3_bucket = mirror["s3_bucket"]
    remote_root = mirror["remote"]
    quote = mirror["quote"]
    cmd_prefix = mirror["command"]
    if "@PATH@" in cmd_prefix:
        cmd_prefix = cmd_prefix.replace("@PATH@", paths["local_prefix"])

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
                if is_s3_bucket:
                    cmd += " " + prefix
                print(cmd, file=f)
                print(file=f)

                # If the file does not exist in the run directory,
                # then copy it from the restart folder.
                # This only has to be done if not using the amazon mirror.
                if not is_s3_bucket:
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
                if is_s3_bucket:
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
                if is_s3_bucket:
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
                if is_s3_bucket:
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
                if is_s3_bucket:
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
                if is_s3_bucket:
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
                if is_s3_bucket:
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
    paths = extract_pathnames_from_log(args)
    paths = expand_restart_file_names(paths, args, run_info)

    # Write a list of unique file paths
    write_unique_paths(paths, args["dryrun_log"] + ".unique")

    # Exit without downloading if skip-download lag was specified
    if args["skip_download"]:
        return

    # Print a message
    if len(args["mirror"]) > 0:
        print("Downloading data from " + args["mirror"])

    # Create script to download missing files from AWS S3
    create_download_script(paths, args)

    #### DEBUG: Uncomment this if you want to see the download script
    #if args["skip_download"]:
    #    return

    # Run the data download script and return the status
    # Remove the file afterwards
    status = subprocess.call(DATA_DOWNLOAD_SCRIPT)
    os.remove(DATA_DOWNLOAD_SCRIPT)

    # Raise an exception if the data was not successfully downloaded
    if status != 0:
        err_msg = "Error downloading data from " + args["mirror"]
        raise Exception(err_msg)


def parse_args():
    """
    Reads global settings from the download_data.yml configuration file.
    Also parses command-line arguments and returns a dictionary
    containing all of these settings.

    Returns:
    --------
    args : dict
        args["config"] : Dict with global settings from download_data.yml
        args["dryrun_log"] Name of the GEOS-Chem dry-run log file
        args["mirror"]: Name of the remote mirror for download
        args["skip_download"]: Are we skipping the download? (T/F)
    """
    dryrun_log = None
    dryrun_found = False
    mirror_found = False
    mirror_remote = None
    skip_download = False
    skip_found = False

    # Read the YAML configuration file
    try:
        config = yaml.load(open("download_data.yml"), Loader=yaml.FullLoader)
    except FileNotFoundError:
        msg = "Could not find configuration file 'download_data.yml'!"
        raise FileNotFoundError(msg)

    # Get a list of mirror names + short names
    mirror_list = list(config["mirrors"].keys())
    short_name_list = []
    for m in mirror_list:
        short_name_list.append(config["mirrors"][m]["short_name"])

    # Parse command-line arguments (argument 0 is the program name)
    for i in range(1, len(sys.argv)):
        arg = sys.argv[i].lower()
        arg = arg.lstrip('-')

        if not dryrun_found:
            dryrun_log = arg
            dryrun_found = True
            continue

        if not mirror_found:
            for m in mirror_list:
                mirror = m.lower()
                short_name = config["mirrors"][m]["short_name"].lower()
                if arg in mirror or arg in short_name:
                    mirror_remote = mirror
                    mirror_found = True
                    continue

        if not skip_found:
            if "skip" in arg:
                skip_download = True
                skip_found = True
                continue


    if dryrun_log is None:
        msg = "The dryrun log file was not supplied!  Exiting ..."
        raise ValueError(msg)

    if mirror_remote is None and not skip_download:
        msg = "Mirror name missing or invalid!  Exiting ..."
        raise ValueError(msg)

    args = {
        "config": config,
        "dryrun_log": dryrun_log,
        "mirror": mirror_remote,
        "skip_download": skip_download
    }
    return args


def main():
    """
    Main program.  Gets command-line arguments and calls function
    download_the_data to initiate a data-downloading process.

    Calling sequence:
    -----------------
        ./download_data.py log MIRROR-NAME
        ./download_data.py log -skip-download  # Print unique log & exit
    """

    # Download the data files from the remote server
    download_the_data(parse_args())


if __name__ == "__main__":
    main()
