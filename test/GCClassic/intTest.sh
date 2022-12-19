#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: intTest.sh
#
# !DESCRIPTION: Runs integration tests on the various GEOS-Chem Classic
#  run directories (interactively, or with a scheduler).
#\\
#\\
# !CALLING SEQUENCE:
#  ./intTest.sh -d root-directory -e env-file [-l] [-p partition] [-q] [-s]
#
#  Where the command-line arguments are as follows:
#
#    -d root-directory : Specify the root folder for integration tests
#    -e env-file       : Specify the environment file (w/ module loads)
#    -l                : Select the LSF scheduler
#    -p partition      : Select partition for SLURM or LSF schedulers
#    -q                : Run a quick set of integration tests (for testing)
#    -s                : Use the SLURM scheduler
#
#  NOTE: you can also use the following long name options:
#
#    --directory root-directory (instead of -d root-directory)
#    --env-file  env-file       (instead of -e env-file      )
#    --lsf                      (instead of -l               )
#    --partition partition      (instead of -p partition     )
#    --quick                    (instead of -q               )
#    --slurm                    (instead of -s               )
#EOP
#------------------------------------------------------------------------------
#BOC

#=============================================================================
# Initialize
#=============================================================================
this="$(basename ${0})"
usage="Usage: ${this} -d root-directory -e env-file [-l] [-p partition] [-s] [-q]"
it_root="none"
env_file="none"
scheduler="none"
sed_cmd="none"
quick="no"

#=============================================================================
# Parse command-line arguments
# See https://www.baeldung.com/linux/bash-parse-command-line-arguments
#=============================================================================

# Call Linux getopt function to specify short & long input options
# (e.g. -d or --directory, etc).  Exit if not succesful
valid_args=$(getopt --options d:e:hlp:qs \
		    --long directory:,env-file:,help,lsf,partition:,quick,slurm -- "$@")
if [[ $? -ne 0 ]]; then
    exit 1;
fi

# Parse arguments and set variables accordingly
eval set -- "${valid_args}"
while [ : ]; do
    case "${1}" in

	# -d or --directory specifies the root folder for tests
	-d | --directory)
	    it_root="${2}"
            shift 2
            ;;

	# -e or --env-file specifies the environment file
	-e | --env-file)
	    env_file="${2}"
            shift 2
            ;;

	# -h or --help prints a help message
	-h | --help)
            echo "$usage"
            exit 1
            ;;

	# -l or --lsf selects the LSF scheduler
	-l | --lsf)
	    scheduler="LSF"
            shift
            ;;

	# -q or --quick runs a quick set of integration tests (for testing)
	-q | --quick)
	    quick="yes"
            shift
	    ;;

	# -p or --partition edits the SLURM and LSF tags
	-p | --partition)
	    sed_cmd="s/REQUESTED_PARTITION/${2}/"
	    shift 2
	    ;;

	# -s or --slurm selects the SLURM scheduler
	-s | --slurm)
	    scheduler="SLURM"
            shift
            ;;

	--) shift;
            break
            ;;
    esac
done

# Error check integration tests root path
if [[ "x${it_root}" == "xnone" ]]; then
    echo "ERROR: The integration test root directory has not been specified!"
    echo "${usage}"
    exit 1
fi

# Error check environment file
if [[ "x${env_file}" == "xnone" ]]; then
    echo "ERROR: The enviroment file (module loads) has not been specified!"
    echo "${usage}"
    exit 1
fi

# Exit if no partition has been selected for SLURM
if [[ "x${scheduler}" == "xSLURM" && "x${sed_cmd}" == "xnone" ]]; then
    echo "ERROR: You must specify a partition for SLURM."
    echo "${usage}"
    exit 1
fi

# Exit if no partition has been selected for SLURM
if [[ "x${scheduler}" == "xLSF" && "x${sed_cmd}" == "xnone" ]]; then
    echo "ERROR: You must specify a partition for LSF."
    echo "${usage}"
    exit 1
fi

#=============================================================================
# Load file with utility functions to setup configuration files
#=============================================================================

# Current directory
this_dir=$(pwd -P)

# Load common functions
. "${this_dir}/commonFunctionsForTests.sh"

#=============================================================================
# Create integration test directories in the root folder
#=============================================================================

# Convert integration test root folder to an absolute path
it_root=$(absolute_path "${it_root}")

# Create GEOS-Chem run directories in the integration test root folder
./intTestCreate.sh "${it_root}" "${env_file}" "${quick}"
if [[ $? -ne 0 ]]; then
   echo "ERROR: Could not create integration test run directories!"
   exit 1
fi

# Change to the integration test root folder
if [[ -d ${it_root} ]]; then
    cd "${it_root}"
else
    echo "ERROR: ${it_root} is not a valid directory!  Exiting..."
    exit 1
fi

#=============================================================================
# Compile the code and run the integration tests
#=============================================================================
if [[ "x${scheduler}" == "xSLURM" ]]; then

    #-------------------------------------------------------------------------
    # Integration tests will run via SLURM
    #-------------------------------------------------------------------------

    # Remove LSF #BSUB tags
    sed_ie '/#BSUB -q REQUESTED_PARTITION/d' "${it_root}/intTestCompile.sh"
    sed_ie '/#BSUB -n 8/d'                   "${it_root}/intTestCompile.sh"
    sed_ie '/#BSUB -W 00:30/d'               "${it_root}/intTestCompile.sh"
    sed_ie '/#BSUB -o lsf-%J.txt/d'          "${it_root}/intTestCompile.sh"
    sed_ie \
     '/#BSUB -R "rusage[mem=8GB] span[ptile=1] select[mem < 1TB]"/d' \
     "${it_root}/intTestCompile.sh"
    sed_ie \
    '/#BSUB -a ''docker(registry.gsc.wustl.edu/sleong/esm:intel-2021.1.2)''/d' \
     "${it_root}/intTestCompile.sh"
    sed_ie '/#BSUB -q REQUESTED_PARTITION/d' "${it_root}/intTestExecute.sh"
    sed_ie '/#BSUB -n 24/d'                  "${it_root}/intTestExecute.sh"
    sed_ie '/#BSUB -W 03:00/d'               "${it_root}/intTestExecute.sh"
    sed_ie '/#BSUB -o lsf-%J.txt/d'          "${it_root}/intTestExecute.sh"
    sed_ie \
    '/#BSUB -R "rusage[mem=90GB] span[ptile=1] select[mem < 2TB]"/d' \
     "${it_root}/intTestExecute.sh"
    sed_ie \
    '/#BSUB -a ''docker(registry.gsc.wustl.edu/sleong/esm:intel-2021.1.2)''/d' \
     "${it_root}/intTestExecute.sh"
    exit 1
    
    # Replace "REQUESTED_PARTITION" with the partition name
    sed_ie "${sed_cmd}" "${it_root}/intTestCompile.sh"
    sed_ie "${sed_cmd}" "${it_root}/intTestExecute.sh"
    
    # Submit compilation tests script
    output=$(sbatch intTestCompile.sh)
    output=($output)
    cmp_id=${output[3]}

    # Submit execution tests script as a job dependency
    output=$(sbatch --dependency=afterok:${cmp_id} intTestExecute.sh)
    output=($output)
    exe_id=${output[3]}

    echo ""
    echo "Compilation tests submitted as SLURM job ${cmp_id}"
    echo "Execution   tests submitted as SLURM job ${exe_id}"

elif [[ "x{$scheduler}" == "xLSF" ]]; then

    #-------------------------------------------------------------------------
    # Integration tests will run via LSF (TODO)
    #-------------------------------------------------------------------------

    # Remove SLURM #SBATCH tags
    sed_ie '/#SBATCH -c 8/d'                   "${it_root}/intTestCompile.sh"
    sed_ie '/#SBATCH -N 1/d'                   "${it_root}/intTestCompile.sh"
    sed_ie '/#SBATCH -t 0-00:30/d'             "${it_root}/intTestCompile.sh"
    sed_ie '/#SBATCH -p REQUESTED_PARTITION/d' "${it_root}/intTestCompile.sh"
    sed_ie '/#SBATCH --mem=8000/d'             "${it_root}/intTestCompile.sh"
    sed_ie '/#SBATCH -p REQUESTED_PARTITION/d' "${it_root}/intTestCompile.sh"
    sed_ie '/#SBATCH --mail-type=END/d'        "${it_root}/intTestCompile.sh"
    sed_ie '/#SBATCH -c 24/d'                  "${it_root}/intTestExecute.sh"
    sed_ie '/#SBATCH -N 1/d'                   "${it_root}/intTestExecute.sh"
    sed_ie '/#SBATCH -t 0-03:00/d'             "${it_root}/intTestExecute.sh" 
    sed_ie '/#SBATCH -p REQUESTED_PARTITION/d' "${it_root}/intTestExecute.sh" 
    sed_ie '/#SBATCH --mem=90000/d'            "${it_root}/intTestExecute.sh"
    sed_ie '/#SBATCH --mail-type=END/d'        "${it_root}/intTestExecute.sh" 

    # Replace "REQUESTED_PARTITION" with the partition name
    sed_ie "${sed_cmd}" "${it_root}/intTestCompile.sh"
    sed_ie "${sed_cmd}" "${it_root}/intTestExecute.sh"

else

    #----------------------------------------------
    # Integration tests will run interactively
    #----------------------------------------------

    # Run compilation tests
    echo ""
    echo "Compiliation tests are running..."
    ./intTestCompile.sh
    if [[ $? != 0 ]]; then
	echo ""
	echo "Compilation tests failed!  Exiting..."
	cd ${this_dir}
	exit 1
    fi
    echo ""
    echo "Compilation tests finished!"

    # Run execution tests
    echo ""
    echo "Execution tests are running..."
    ./intTestExecute.sh
    if [[ $? != 0 ]]; then
	echo ""
	echo "Execution tests failed!  Exiting..."
	cd ${this_dir}
	exit 1
    fi
    echo ""
    echo "Execution tests finished!"

    # Change back to this directory
    cd ${this_dir}

fi

#=============================================================================
# Cleanup and quit
#=============================================================================

# Free local variables
unset env_file
unset it_root
unset cmp_id
unset exe_id
unset quick
unset output
unset scheduler
unset this_dir

# Free imported variables
unset FILL
unset SEP_MAJOR
unset SEP_MINOR
unset CMP_PASS_STR
unset CMP_FAIL_STR
unset EXE_PASS_STR
unset EXE_FAIL_STR
#EOC
