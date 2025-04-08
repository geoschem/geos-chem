#!/usr/bin/env python3
"""
Creates a ReadTheDocs-style list table containing a list of species
for a given type of GEOS-Chem simulation.

Calling sequence:
$ python -m rtd_species_by_simulation SIMULATION-NAME TABLE-FILE
"""
import os
from sys import argv
from yaml import safe_load
from gcpy.constants import ENCODING
from gcpy.util import read_config_file, verify_variable_type

# Constants
GC_CONFIG = f"../GCClassic/geoschem_config.yml.templates/geoschem_config.yml."
SPECIES_DB = f"./species_databaseXX.yml"
KPP_FILE = f"../../KPP/XX/XX.eqn"


def read_species_from_eqn_file(simulation_name):
    """
    Reads a geoschem_config.yml template file for a given GEOS-Chem
    simulation and returns the list of transported species

    Args
    simulation_name     : str : Name of the GEOS-Chem simulation

    returns
    transported_species : list : List of transported species
    """
    verify_variable_type(simulation_name, str)

    # Read YAML file, replacing variables with placeholder text
    filename = os.path.expanduser(f"{KPP_FILE}").replace("XX", simulation_name)

    # Species contain the word IGNORE; take the 1st substring
    species = []
    with open(filename, 'r', encoding=ENCODING) as ifile:
        for line in ifile:
            line = line.strip()
            if "IGNORE" in line:
                species.append(line.split("=")[0].strip())

    return species
                

def read_transported_species(simulation_name):
    """
    Reads a geoschem_config.yml template file for a given GEOS-Chem
    simulation and returns the list of transported species

    Args
    simulation_name     : str : Name of the GEOS-Chem simulation

    returns
    transported_species : list : List of transported species
    """
    verify_variable_type(simulation_name, str)

    # We need to define fake values for the various environment
    # variables in order to be able to read these with PyYaml
    os.environ['RUNDIR_SIM_NAME'] = 'X'
    os.environ['RUNDIR_SIM_START_DATE'] = 'X'
    os.environ['RUNDIR_SIM_START_TIME'] = 'X'
    os.environ['RUNDIR_SIM_END_DATE'] = 'X'
    os.environ['RUNDIR_SIM_END_TIME'] = 'X'
    os.environ['RUNDIR_DATA_ROOT'] = 'X'
    os.environ['RUNDIR_MET'] = 'X'
    os.environ['RUNDIR_USE_GCCLASSIC_TIMERS'] = 'X'
    os.environ['RUNDIR_GRID_RES_LONG'] = 'X'
    os.environ['RUNDIR_GRID_NLEV'] = 'X'
    os.environ['RUNDIR_GRID_LON_RANGE'] = 'X'
    os.environ['RUNDIR_CENTER_LON_180'] = 'X'
    os.environ['RUNDIR_GRID_LAT_RANGE'] = 'X'
    os.environ['RUNDIR_GRID_HALF_POLAR'] = 'X'
    os.environ['RUNDIR_GRID_NESTED_SIM'] = 'X'
    os.environ['RUNDIR_GRID_BUFFER_ZONE'] = 'X'
    os.environ['RUNDIR_TRANSPORT_TS'] = 'X'
    os.environ['RUNDIR_CHEMISTRY_TS'] = 'X'
    os.environ['RUNDIR_USE_NLPBL'] = 'X'

    # Read YAML file, replacing variables with placeholder text
    filename = os.path.expanduser(f"{GC_CONFIG}{simulation_name}")
    with open(filename, 'r', encoding=ENCODING) as ifile:
        config = safe_load(os.path.expandvars(ifile.read()))

    return config["operations"]["transport"]["transported_species"]


def read_species_database():
    """
    Reads the GEOS-Chem Species Database.

    Returns
    species_database : dict : GEOS-Chem Species Database object
    """
    count = 0
    for sim in ["", "_apm", "_hg", "_tomas"]:
        filename = os.path.expanduser(SPECIES_DB.replace("XX", sim))
        if count == 0:
            species_database = read_config_file(filename, quiet=True)
        else:
            species_database = species_database | read_config_file(
                filename, quiet=True
            )
        count += 1

    return species_database


def get_species_metadata(species, species_database):
    """
    Compares a species against the species database and returns
    selected metadata fields.

    Args
    species          : list : List of species names
    species_database : dict : GEOS-Chem Species Database

    Returns
    metadata         : dict : Selected species metadata fields
    """
    verify_variable_type(species, list)
    verify_variable_type(species_database, dict)

    metadata = {}
    for spc in species:
        if spc in species_database:
            formula = "not listed"
            if "Formula" in species_database[spc]:
                formula = species_database[spc]["Formula"]
            fullname = "not listed"
            if "FullName" in species_database[spc]:
                fullname =  species_database[spc]["FullName"]
            metadata[spc] = {
                "MW_g": species_database[spc]["MW_g"],
                "FullName": fullname,
                "Formula": formula,
            }

    return metadata


def create_rtd_list_table(metadata, table_file):
    """
    Creates a ReadTheDocs list table displaying the transported species
    for a given GEOS-Chem simulation.

    Args
    metadata   : dict : Selected species metadata fields
    table_file : str  : File where the list table will be written
    """
    verify_variable_type(metadata, dict)
    verify_variable_type(table_file, str)

    with open(table_file, "w", encoding=ENCODING) as ofile:

        # Table declaration
        print(".. list-table::", file=ofile)
        print("   :header-rows: 1", file=ofile)
        print("   :align: left\n", file=ofile)

        # Header
        print("   * - Species", file=ofile)
        print("     - Description", file=ofile)
        print("     - Formula", file=ofile)
        print("     - MW (g)", file=ofile)

        # Species and species metadata
        for (species, metadata_fields) in metadata.items():
            print(f"   * - {species}", file=ofile)
            print(f"     - {metadata_fields['FullName']}", file=ofile)
            print(f"     - {metadata_fields['Formula']}", file=ofile)
            print(f"     - {metadata_fields['MW_g']}", file=ofile)


def main(simulation_name, table_file):
    """
    Main program.  Calls subroutines to create the list table.

    Args
    simulation_name : str : Name of a GEOS-Chem simulation
    table_file      : str : File where the list table will be written
    """
    verify_variable_type(simulation_name, str)
    verify_variable_type(table_file, str)

    # Get a list of species for the simulation from either
    # the *.kpp file or from the geoschem_config.yml file
    if "fullchem" in simulation_name or \
       "Hg" in simulation_name or \
       "carbon" in simulation_name:
        species = read_species_from_eqn_file(simulation_name)
    else:
        species = read_transported_species(simulation_name)
        
    # Read the species database info
    species_database = read_species_database()

    # Compare the list of transported species to the
    # species database and return a dict of metadata
    metadata = get_species_metadata(species, species_database)

    # Create a list table in ReadTheDocs format
    create_rtd_list_table(metadata, table_file)


if __name__ == '__main__':

    if len(argv) != 3:
        raise ValueError("Usage: python -m sims SIMULATION-NAME TABLE-FILE")

    main(argv[1], argv[2])
