#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: parTest.sh
#
# !DESCRIPTION: Runs parallelization tests on the various GEOS-Chem Classic
#  run directories (interactively, or with a scheduler).
#\\
#\\
# !CALLING SEQUENCE:
#  ./parTest.sh -d root-directory -e env-file [-l] [-p partition] [-q] [-s]
#
#  Where the command-line arguments are as follows:
#
#    -d root-directory : Specify the root folder for parallelization tests
#    -e env-file       : Specitify the environment file (w/ module loads)
#    -l                : Select the LSF scheduler
#    -p partition      : Select partition for SLURM or LSF schedulers
#    -q                : Run a quick set of parallelization tests (for testing)
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
ptRoot="none"
envFile="none"
scheduler="none"
sedCmd="none"
quick="no"

#=============================================================================
# Parse command-line arguments
# See https://www.baeldung.com/linux/bash-parse-command-line-arguments
#=============================================================================

# Call Linux getopt function to specify short & long input options
# (e.g. -d or --directory, etc).  Exit if not succesful
validArgs=$(getopt --options d:e:hlp:qs \
  --long directory:,env-file:,help,lsf,partition:,quick,slurm -- "$@")
if [[ $? -ne 0 ]]; then
    exit 1;
fi

# Parse arguments and set variables accordingly
eval set -- "${validArgs}"
while [ : ]; do
    case "${1}" in

	# -d or --directory specifies the root folder for tests
	-d | --directory)
	    ptRoot="${2}"
            shift 2
            ;;

	# -e or --env-file specifies the environment file
	-e | --env-file)
	    envFile="${2}"
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

	# -q or --quick runs a quick set of parallelization tests (for testing)
	-q | --quick)
	    quick="yes"
            shift
	    ;;

	# -p or --partition edits the SLURM and LSF tags
	-p | --partition)
	    sedCmd="s/REQUESTED_PARTITION/${2}/"
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

# Error check parallelization tests root path
if [[ "x${ptRoot}" == "xnone" ]]; then
    echo "ERROR: The parallelization test root directory has not been specified!"
    echo "${usage}"
    exit 1
fi

# Error check environment file
if [[ "x${envFile}" == "xnone" ]]; then
    echo "ERROR: The enviroment file (module loads) has not been specified!"
    echo "${usage}"
    exit 1
fi

# Exit if no partition has been selected for SLURM
if [[ "x${scheduler}" == "xSLURM" && "x${sedCmd}" == "xnone" ]]; then
    echo "ERROR: You must specify a partition for SLURM."
    echo "${usage}"
    exit 1
fi

# Exit if no partition has been selected for SLURM
if [[ "x${scheduler}" == "xLSF" && "x${sedCmd}" == "xnone" ]]; then
    echo "ERROR: You must specify a partition for LSF."
    echo "${usage}"
    exit 1
fi

#=============================================================================
# Load file with utility functions to setup configuration files
#=============================================================================

# Current directory
thisDir=$(pwd -P)

# Load common functions
. "${thisDir}/commonFunctionsForTests.sh"

#=============================================================================
# Create parallelization test directories in the root folder
#=============================================================================

# Convert parallelization test root folder to an absolute path
ptRoot=$(absolute_path "${ptRoot}")

# Create GEOS-Chem run directories in the parallelization test root folder
./parTestCreate.sh "${ptRoot}" "${envFile}" "${quick}"
if [[ $? -ne 0 ]]; then
   echo "ERROR: Could not create parallelization test run directories!"
   exit 1
fi

# Change to the parallelization test root folder
if [[ -d ${ptRoot} ]]; then
    cd "${ptRoot}"
else
    echo "ERROR: ${ptRoot} is not a valid directory!  Exiting..."
    exit 1
fi

#=============================================================================
# Compile the code and run the parallelization tests
#=============================================================================
if [[ "x${scheduler}" == "xSLURM" ]]; then

    #-------------------------------------------------------------------------
    # Parallelization tests will run via SLURM
    #-------------------------------------------------------------------------

    # Remove LSF #BSUB tags
    sed_ie '/#BSUB -q REQUESTED_PARTITION/d' "${ptRoot}/parTestCompile.sh"
    sed_ie '/#BSUB -n 8/d'                   "${ptRoot}/parTestCompile.sh"
    sed_ie '/#BSUB -W 00:30/d'               "${ptRoot}/parTestCompile.sh"
    sed_ie '/#BSUB -o lsf-%J.txt/d'          "${ptRoot}/parTestCompile.sh"
    sed_ie \
	'/#BSUB -R "rusage\[mem=8GB\] span\[ptile=1\] select\[mem < 1TB\]"/d' \
	"${ptRoot}/parTestCompile.sh"
    sed_ie \
	"/#BSUB -a 'docker(registry\.gsc\.wustl\.edu\/sleong\/esm\:intel\-2021\.1\.2)'/d" \
	"${ptRoot}/parTestCompile.sh"
    sed_ie '/#BSUB -q REQUESTED_PARTITION/d' "${ptRoot}/parTestExecute.sh"
    sed_ie '/#BSUB -n 24/d'                  "${ptRoot}/parTestExecute.sh"
    sed_ie '/#BSUB -W 3:30/d'                "${ptRoot}/parTestExecute.sh"
    sed_ie '/#BSUB -o lsf-%J.txt/d'          "${ptRoot}/parTestExecute.sh"
    sed_ie \
	'/#BSUB -R "rusage\[mem=90GB\] span\[ptile=1\] select\[mem < 2TB\]"/d' \
	"${ptRoot}/parTestExecute.sh"
    sed_ie \
	"/#BSUB -a 'docker(registry\.gsc\.wustl\.edu\/sleong\/esm\:intel\-2021\.1\.2)'/d" \
	"${ptRoot}/parTestExecute.sh"

    # Replace "REQUESTED_PARTITION" with the partition name
    sed_ie "${sedCmd}" "${ptRoot}/parTestCompile.sh"
    sed_ie "${sedCmd}" "${ptRoot}/parTestExecute.sh"

    # Submit compilation tests script
    output=$(sbatch parTestCompile.sh)
    output=($output)
    cmpId=${output[3]}

    # Submit execution tests script as a job dependency
    output=$(sbatch --dependency=afterok:${cmpId} parTestExecute.sh)
    output=($output)
    exeId=${output[3]}

    echo ""
    echo "Compilation tests submitted as SLURM job ${cmpId}"
    echo "Execution   tests submitted as SLURM job ${exeId}"

elif [[ "x${scheduler}" == "xLSF" ]]; then

    #-------------------------------------------------------------------------
    # Parallelization tests will run via LSF
    #-------------------------------------------------------------------------

    # Remove SLURM #SBATCH tags
    sed_ie '/#SBATCH -c 8/d'                   "${ptRoot}/parTestCompile.sh"
    sed_ie '/#SBATCH -N 1/d'                   "${ptRoot}/parTestCompile.sh"
    sed_ie '/#SBATCH -t 0-00:30/d'             "${ptRoot}/parTestCompile.sh"
    sed_ie '/#SBATCH -p REQUESTED_PARTITION/d' "${ptRoot}/parTestCompile.sh"
    sed_ie '/#SBATCH --mem=8000/d'             "${ptRoot}/parTestCompile.sh"
    sed_ie '/#SBATCH -p REQUESTED_PARTITION/d' "${ptRoot}/parTestCompile.sh"
    sed_ie '/#SBATCH --mail-type=END/d'        "${ptRoot}/parTestCompile.sh"
    sed_ie '/#SBATCH -c 24/d'                  "${ptRoot}/parTestExecute.sh"
    sed_ie '/#SBATCH -N 1/d'                   "${ptRoot}/parTestExecute.sh"
    sed_ie '/#SBATCH -t 0-03:30/d'             "${ptRoot}/parTestExecute.sh"
    sed_ie '/#SBATCH -p REQUESTED_PARTITION/d' "${ptRoot}/parTestExecute.sh"
    sed_ie '/#SBATCH --mem=90000/d'            "${ptRoot}/parTestExecute.sh"
    sed_ie '/#SBATCH --mail-type=END/d'        "${ptRoot}/parTestExecute.sh"

    # Replace "REQUESTED_PARTITION" with the partition name
    sed_ie "${sedCmd}" "${ptRoot}/parTestCompile.sh"
    sed_ie "${sedCmd}" "${ptRoot}/parTestExecute.sh"

    # Submit compilation tests script
    output=$(bsub parTestCompile.sh)
    output=($output)
    cmpId=${output[1]}
    cmpId=${cmpId/<}
    cmpId=${cmpId/>}

    # Submit execution tests script as a job dependency
    output=$(bsub -w "exit(${cmpId},0)" parTestExecute.sh)
    output=($output)
    exeId=${output[1]}
    exeId=${exeId/<}
    exeId=${exeId/>}

else

    #-------------------------------------------------------------------------
    # Parallelization tests will run interactively
    #-------------------------------------------------------------------------

    # Run compilation tests
    echo ""
    echo "Compiliation tests are running..."
    ./parTestCompile.sh &

    # Change back to this directory
    cd "${thisDir}"

fi

#=============================================================================
# Cleanup and quit
#=============================================================================

# Free local variables
unset cmpId
unset envFile
unset exeId
unset ptRoot
unset quick
unset output
unset scheduler
unset thisDir

# Free imported variables
unset FILL
unset SEP_MAJOR
unset SEP_MINOR
unset CMP_PASS_STR
unset CMP_FAIL_STR
unset EXE_PASS_STR
unset EXE_FAIL_STR
#EOC
