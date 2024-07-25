#!/usr/bin/env python3

"""
Description:
------------
This Python script (assumes Python3) reads a GEOS-Chem or
HEMCO-standalone log file containing dry-run output and does
the following:

    (1) Creates a list of unique files that are required for the
        GEOS-Chem or HEMCO-standalone simulation;

    (2) Creates a bash script to download missing files from the AWS
        s3://gcgrid bucket or from a specified server;

    (3) Executes the bash script to download the necessary data;

    (4) Removes the bash script upon successful download.


Remarks:
--------
    (1) This script only requires the "os", "sys", "subprocess", and
        PyYaml packages.

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
GEOSCHEM_INPUT_FILE = "./geoschem_config.yml"
DATA_DOWNLOAD_SCRIPT = "./auto_generated_download_script.sh"


def read_config_file(
        config_file,
        to_str=False
):
    """
    Reads configuration information from a YAML file.

    Args:
    -----
    config_file : str
        The configuration file in YAML format
    to_str : bool
        Set this to True if you wish to return the data in the YAML
        file as strings, or False otherwise.

    Returns:
    --------
    config : dict
        Dictionary with the contents of the YAML file
    """
    try:
        with open(config_file, encoding="UTF-8") as stream:
            if to_str:
                return yaml.load(stream, Loader=yaml.loader.BaseLoader)
            return yaml.load(stream, Loader=yaml.loader.SafeLoader)
    except FileNotFoundError as err:
        msg = f"Error reading configuration in {config_file}: {err}"
        raise FileNotFoundError(msg) from err


def extract_pathnames_from_log(
        args
):
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
    with open(dryrun_log, "r", encoding="UTF-8") as ifile:

        # Read data from the file line by line.
        # Add file paths to the data_list set.
        line = ifile.readline()
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
            line = ifile.readline()

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
            msg = \
                "Could not locate the ExtData folder in your local disk space!"
            raise ValueError(msg)

        # Close file and return
        # The "sorted" command will return unique values
        ifile.close()

        paths = {
            "comments": comments,
            "found": found,
            "missing": missing,
            "local_prefix": local_prefix
        }
        return paths


def get_run_info():
    """
    Searches through the geoschem_config.yml file for GEOS-Chem
    simulation parameters.

    Returns:
    -------
    run_info : dict
        Contains the GEOS-Chem run parameters: start_date,
        start_time, end_date, end_time, met, grid, and sim.
    """

    # Read GEOS-Chem configuration file
    config = read_config_file(
        GEOSCHEM_INPUT_FILE,
        to_str=True
    )

    # Create dictionary with GEOS-Chem simulation parameters
    # NOTE: Numbers are returned as strings, and need to be converted
    run_info = {}
    run_info["nest"] = ""
    run_info["tomas15"] = False
    run_info["tomas40"] = False
    run_info["start_date"] = int(
        config["simulation"]["start_date"][0]
    )
    run_info["start_time"] = int(
        config["simulation"]["start_date"][1]
    )
    run_info["end_date"] = int(
        config["simulation"]["end_date"][0]
    )
    run_info["end_time"] = int(
        config["simulation"]["end_date"][1]
    )
    run_info["met_field"] = config["simulation"]["met_field"]
    run_info["sim"] = config["simulation"]["name"]
    run_info["resolution"] = config["grid"]["resolution"]
    run_info["grid"] = get_grid_suffix(
        run_info["resolution"]
    )
    run_info["nest"] = get_nest_suffix(
        config["grid"]["longitude"]["range"]
    )
    run_info["tomas15"] = \
        "NK15" in config["operations"]["transport"]["transported_species"]
    run_info["tomas40"] = \
        "NK40" in config["operations"]["transport"]["transported_species"]

    return run_info


def get_grid_suffix(
        resolution
):
    """
    Given a model resolution, returns the grid filename suffix.

    Args:
    -----
    resolution : str
        The grid resolution (read from geoschem_config.yml)

    Returns:
    --------
    suffix : str
        The corresponding filename suffix
    """
    if "4.0x5.0" in resolution:
        return "4x5"
    if "2.0x2.5" in resolution:
        return "2x25"
    if "0.5x0.625" in resolution:
        return "05x0625"
    return "025x03125"


def get_nest_suffix(
        longitude
):
    """
    Given a model resolution, returns the nested-grid suffix.

    Args:
    -----
    resolution : str
        The grid resolution (read from geoschem_config.yml)

    Returns:
    --------
    suffix : str
        The corresponding nested-grid suffix
    """
    if "-130" in longitude or "-140" in longitude:
        return "na"
    if "60" in longitude or "70" in longitude:
        return "as"
    return ""


def get_remote_restart_filename(
        local_prefix,
        run_info,
        rst_info
):
    """
    Returns the remote restart file name for a given
    GEOS-Chem Classic simulation type.

    Args:
    -----
    local_prefix : str
        The root data folder.  ExtData is a subfolder of this folder.
    run_info : dict
        Information read from geoschem_config.yml
    rst_info : dict
        Restart file paths (local and remote), as read from the
        download_data.yml configuration file.

    Returns:
    --------
    remote_rst : str
        Path to the remote restart file
    """

    # Simulation type
    simulation = run_info["sim"].lower()

    # Remote restart file directory
    root = os.path.join(local_prefix, "ExtData", rst_info["root"])

    # Special handling for fullchem
    if "fullchem" in simulation:
        if run_info["tomas15"] is True:
            return os.path.join(root, rst_info["tomas15"]["remote"])
        if run_info["tomas40"] is True:
            return os.path.join(root, rst_info["tomas40"]["remote"])
        return os.path.join(root, rst_info["fullchem"]["remote"])

    # Special handling for mercury
    if "mercury" in simulation or "hg" in simulation:
        return os.path.join(root, rst_info["mercury"]["remote"])

    # All other simulations use the lowercase simulation name
    return os.path.join(root, rst_info[simulation]["remote"])


def replace_entry_in_list(
        the_list,
        old_entry,
        new_entry
):
    """
    Replaces a string entry in a list with a new entry.

    Args:
    -----
    the_list : list of str
       The list
    old_entry : (str
        Entry to replace
    new_entry : str
        Replacement text

    Returns:
    --------
    the_list : list of str
        The modified list
    """
    return list(map(lambda x: x.replace(old_entry, new_entry), the_list))


def expand_restart_file_names(
        paths,
        args,
        run_info
):
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

    # Get the name of the remote restart file for this simulation
    remote_rst = get_remote_restart_filename(
        paths["local_prefix"],
        run_info,
        args["config"]["restarts"]
    )

    # First, look for the restart file name in the found files
    do_exit = False
    for path in paths["found"]:
        if "GEOSChem.Restart" in path:
            new_path = path + " --> " + remote_rst
            do_exit = True
            paths["found"] = replace_entry_in_list(
                paths["found"],
                path,
                new_path
            )
            break
    paths["found"] = sorted(paths["found"])
    if do_exit:
        return paths

    # Then, look for the restart file name in the missing files
    for path in paths["missing"]:
        if "GEOSChem.Restart" in path:
            new_path = path + " --> " + remote_rst
            do_exit = True
            paths["missing"] = replace_entry_in_list(
                paths["missing"],
                path,
                new_path
            )
            break
    paths["missing"] = sorted(paths["missing"])

    # Return the updated data paths
    return paths


def write_unique_paths(
        paths,
        unique_log
):
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
        with open(unique_log, "w", encoding="UTF-8") as ofile:
            for comment in paths["comments"]:
                print(comment, file=ofile)
            for path in combined_paths:
                print(path, file=ofile)
            for comment in paths["comments"]:
                print(comment, file=ofile)
        ofile.close()
        print(f"Log with unique file paths written to: {unique_log}")
    except RuntimeError as exc:
        raise RuntimeError(f"Could not write {unique_log}") from exc


def create_download_script(
        paths,
        args
):
    """
    Creates a data download script to obtain missing files
    from the s3://gcgrid bucket on the AWS cloud or from a
    specified server.

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
    with open(DATA_DOWNLOAD_SCRIPT, "w", encoding="UTF-8") as ofile:

        # Write shebang line to script
        print("#!/bin/bash\n", file=ofile)
        print("# This script was generated by download_data.py\n", file=ofile)

        # Write download commands for only the missing data files
        for path in paths["missing"]:

            if "-->" in path:

                # ------------------------------------------------------
                # Edge case: Linked restart files
                # ------------------------------------------------------

                # First copy the restart file to local ExtData
                remote_rst = (path.split("-->")[1]).strip()
                local_rst = (path.split("-->")[0]).strip()
                index2 = remote_rst.find("ExtData") + 7
                prefix = local_rst
                extdata = remote_rst[:index2]
                remote_rst = remote_root + remote_rst[index2:]
                cmd = cmd_prefix + quote + remote_rst + quote
                if is_s3_bucket:
                    cmd += " " + prefix
                print(cmd, file=ofile)
                print(file=ofile)

                # If the file does not exist in the run directory,
                # then copy it from the restart folder.
                # This only has to be done if not using the amazon mirror.
                if not is_s3_bucket:
                    if not os.path.exists(local_rst):
                        index3 = remote_rst.find("GEOSCHEM_RESTARTS")
                        rst = os.path.join(extdata, remote_rst[index3:])
                        cmd = "cp -f " + rst + " " + local_rst
                        print(cmd, file=ofile)
                        print(file=ofile)

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
                print(cmd, file=ofile)

                # Rename it to IPMN
                cmd = "mv " + local_dir + "/gmi.clim.PMN.geos5.2x25.nc " + \
                      local_dir + "/gmi.clim.IPMN.geos5.2x25.nc"
                print(cmd, file=ofile)

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
                print(cmd, file=ofile)

                # Rename it to NPMN
                cmd = "mv " + local_dir + "/gmi.clim.PMN.geos5.2x25.nc " + \
                      local_dir + "/gmi.clim.NPMN.geos5.2x25.nc"
                print(cmd, file=ofile)
                print(file=ofile)

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
                print(cmd, file=ofile)

                # Rename it to NPMN
                cmd = "mv " + local_dir + "/gmi.clim.RIP.geos5.2x25.nc " + \
                      local_dir + "/gmi.clim.RIPA.geos5.2x25.nc"
                print(cmd, file=ofile)
                print(file=ofile)

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
                print(cmd, file=ofile)

                # Rename it to RIPB
                cmd = "mv " + local_dir + "/gmi.clim.RIP.geos5.2x25.nc " + \
                      local_dir + "/gmi.clim.RIPB.geos5.2x25.nc"
                print(cmd, file=ofile)
                print(file=ofile)

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
                print(cmd, file=ofile)

                # Rename it to RIPD
                cmd = "mv " + local_dir + "/gmi.clim.RIP.geos5.2x25.nc " + \
                      local_dir + "/gmi.clim.RIPD.geos5.2x25.nc"
                print(cmd, file=ofile)
                print(file=ofile)

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
                print(cmd, file=ofile)
                print(file=ofile)

        # Kludge: Create a ExtData/CHEM_INPUTS folder if it
        # does not exist. This will prevent abnormal exits.
        chem_dir = paths["local_prefix"] + 'ExtData/CHEM_INPUTS'
        cmd = f"if [[ ! -d {chem_dir} ]]; then mkdir {chem_dir}; fi"
        print(cmd, file=ofile)
        print(file=ofile)

        # Close file and make it executable
        ofile.close()
        os.chmod(DATA_DOWNLOAD_SCRIPT, 0o755)


def download_the_data(
        args
):
    """
    Downloads GEOS-Chem data files from the AWS s3://gcgrid bucket
    or from a specified server.

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

    # Exit without downloading if skip-download flag was specified
    if args["skip_download"]:
        return

    # Print a message
    if len(args["mirror"]) > 0:
        print(f"Downloading data from {args['mirror']}")

    # Create script to download missing files from AWS S3
    create_download_script(paths, args)

    #### DEBUG: Uncomment this if you wish to see the download script
    #if args["skip_download"]:
    #    return

    # Run the data download script and return the status
    # Remove the file afterwards
    status = subprocess.call(DATA_DOWNLOAD_SCRIPT)
    os.remove(DATA_DOWNLOAD_SCRIPT)

    # Raise an exception if the data was not successfully downloaded
    if status != 0:
        msg = f"Error downloading data from {args['mirror']}"
        raise RuntimeError(msg)


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
    config = read_config_file("download_data.yml")

    # Get a list of mirror names + short names
    mirror_list = list(config["mirrors"].keys())
    short_name_list = []
    for mir in mirror_list:
        short_name_list.append(config["mirrors"][mir]["short_name"])

    # Parse command-line arguments (argument 0 is the program name)
    for i in range(1, len(sys.argv)):
        arg = sys.argv[i].lower()
        arg = arg.lstrip('-')

        if not dryrun_found:
            dryrun_log = arg
            dryrun_found = True
            continue

        if not mirror_found:
            for mir in mirror_list:
                mirror = mir.lower()
                short_name = config["mirrors"][mir]["short_name"].lower()
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
