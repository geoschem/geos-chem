#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: parallelTest.sh
#
# !DESCRIPTION: Runs parallelization tests on the various GEOS-Chem Classic
#  run directories (interactively, or with a scheduler).
#\\
#\\
# !CALLING SEQUENCE:
#  ./parallelTest.sh -d root-dir -t compile|all [-e env-file] [-h] [-q]
#
#  Required arguments
#    -d root-dir     : Specify the root folder for integration tests
#    -t compile|all  : Specify the tests to run (compile-only or all)
#
#  Optional arguments
#    -e env-file     : Software environment file for Harvard Cannon
#    -h              : Display a help message
#    -n              : Do not bootstrap missing restart file variables
#    -q              : Run a quick set of integration tests (for testing)
#
#  NOTE: you can also use the following long name options:
#
#    --directory     (instead of -d)
#    --env-file      (instead of -e)
#    --help          (instead of -h)
#    --no-bootstrap  (instead of -n)
#    --quick         (instead of -q)
#    --tests-to-run  (instead of -t)
#EOP
#------------------------------------------------------------------------------
#BOC

#=============================================================================
# Initialize
#=============================================================================
this="$(basename ${0})"
usage="Usage: ${this} -d root-dir -t compile|all [-e env-file] [-h] [-q]"
quick="NO"

# Current directory
thisDir=$(pwd -P)
cd "${thisDir}"

# Load common functions
. "${thisDir}/../../shared/commonFunctionsForTests.sh"

#=============================================================================
# Parse command-line arguments
# See https://www.baeldung.com/linux/bash-parse-command-line-arguments
#=============================================================================

# Call Linux getopt function to specify short & long input options
# (e.g. -d or --directory, etc).  Exit if not succesful
validArgs=$(getopt --options d:e:hqt: \
  --long directory:,env-file:,help,quick,tests-to-run: -- "$@")
if [[ $? -ne 0 ]]; then
    exit 1;
fi

# Parse arguments and set variables accordingly
# NOTE: Convert some inputs to uppercase to facilitate comparisons
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

	# -q or --quick runs a quick set of parallelization tests (for testing)
	-q | --quick)
	    quick="YES"
            shift
	    ;;
	
	# -t or --tests-to-run specifies the type of tests to run
	-t | --tests-to-run)
            testsToRun="${2^^}"
            shift 2
            ;;

	--) shift;
            break
            ;;
    esac
done

#=============================================================================
# Sanity-check user input
#=============================================================================

# Get the site name from the node name
site=$(get_site_name)

# Error check integration tests root path
if [[ "X${ptRoot}" == "X" ]]; then
    echo "ERROR: The parallelization test root directory has not been specified!"
    echo "${usage}"
    exit 1
fi

# Error check the type of tests to run
if [[ "X${testsToRun}" == "X" ]]; then
    echo "ERROR: You must specify the test type: compile|all"
    echo "${usage}"
    exit 1
fi
if [[ "X${testsToRun}" != "XCOMPILE" && "X${testsToRun}" != "XALL" ]]; then
    echo "ERROR: Invalid selction for tests-to-run, must be: compile|all"
    echo "${usage}"
    exit 1
fi

# Error checks for tests that include compile & run phases
if [[ "X${testsToRun}" == "XALL" ]]; then
    
    # Use the default environment file for Cannon if not specified
    if [[ "X${site}" == "XCANNON" && "X${envFile}" == "X" ]]; then
	envFile=$(get_default_gcc_env_file)
    fi

    # Get the sed command that will replace the partition name
    sedPartitionCmd=$(get_sed_partition_cmd_from_site "${site}")
fi

#=============================================================================
# Create parallelization test directories in the root folder
#=============================================================================

# Convert parallelization test root folder to an absolute path
ptRoot=$(absolute_path "${ptRoot}")

# Prevent running integration tests in the source code directory tree
if [[ "$(absolute_path ${thisDir})" =~ "${ptRoot}" ]]; then
    echo "ERROR: You cannot run parallel tests in the source code directory!"
    exit 1
fi

# Create GEOS-Chem run directories in the parallelization test root folder
./parallelTestCreate.sh "${ptRoot}" "${envFile}" "${testsToRun}" "${quick}"
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

# Define local convenience variables
logsDir="${ptRoot}/${LOGS_DIR}"
scriptsDir="${ptRoot}/${SCRIPTS_DIR}"

# Navigate to the logs directory (so all output will be placed there)
cd "${logsDir}"

#=============================================================================
# Run the tests
#=============================================================================
if [[ "X${testsToRun}" == "XCOMPILE" ]]; then

    #-------------------------------------------------------------------------
    # Compilation-only tests (scheduler is not used)
    #-------------------------------------------------------------------------
    echo ""
    echo "Compiliation tests are running..."
    ${scriptsDir}/parallelTestCompile.sh &    

elif [[ "X${testsToRun}" == "XALL" && "X${site}" == "XCANNON" ]]; then

    #-------------------------------------------------------------------------
    # Compilation and execution tests on Harvard Cannon (via SLURM)
    #-------------------------------------------------------------------------

    # Remove LSF #BSUB tags
    sed_ie '/#BSUB -q REQUESTED_PARTITION/d' "${scriptsDir}/parallelTestCompile.sh"
    sed_ie '/#BSUB -n 8/d'                   "${scriptsDir}/parallelTestCompile.sh"
    sed_ie '/#BSUB -W 0:30/d'                "${scriptsDir}/parallelTestCompile.sh"
    sed_ie '/#BSUB -o lsf-%J.txt/d'          "${scriptsDir}/parallelTestCompile.sh"
    sed_ie \
	'/#BSUB -R "rusage\[mem=8GB\] span\[ptile=1\] select\[mem < 1TB\]"/d' \
	"${scriptsDir}/parallelTestCompile.sh"
    sed_ie \
	"/#BSUB -a 'docker(registry\.gsc\.wustl\.edu\/sleong\/esm\:intel\-2021\.1\.2)'/d" \
	"${scriptsDir}/parallelTestCompile.sh"
    sed_ie '/#BSUB -q REQUESTED_PARTITION/d' "${scriptsDir}/parallelTestExecute.sh"
    sed_ie '/#BSUB -n 24/d'                  "${scriptsDir}/parallelTestExecute.sh"
    sed_ie '/#BSUB -W 6:00/d'                "${scriptsDir}/parallelTestExecute.sh"
    sed_ie '/#BSUB -o lsf-%J.txt/d'          "${scriptsDir}/parallelTestExecute.sh"
    sed_ie \
	'/#BSUB -R "rusage\[mem=90GB\] span\[ptile=1\] select\[mem < 2TB\]"/d' \
	"${scriptsDir}/parallelTestExecute.sh"
    sed_ie \
	"/#BSUB -a 'docker(registry\.gsc\.wustl\.edu\/sleong\/esm\:intel\-2021\.1\.2)'/d" \
	"${scriptsDir}/parallelTestExecute.sh"

    # Replace "REQUESTED_PARTITION" with the partition name
    sed_ie "${sedPartitionCmd}" "${scriptsDir}/parallelTestCompile.sh"
    sed_ie "${sedPartitionCmd}" "${scriptsDir}/parallelTestExecute.sh"

    # Submit compilation tests script
    output=$(sbatch ${scriptsDir}/parallelTestCompile.sh)
    output=($output)
    cmpId=${output[3]}

    # Submit execution tests script as a job dependency
    output=$(sbatch --dependency=afterok:${cmpId} ${scriptsDir}/parallelTestExecute.sh)
    output=($output)
    exeId=${output[3]}

    echo ""
    echo "Compilation tests submitted as SLURM job ${cmpId}"
    echo "Execution   tests submitted as SLURM job ${exeId}"

elif [[ "X${testsToRun}" == "XALL" && "X${site}" == "XCOMPUTE1" ]]; then

    #-------------------------------------------------------------------------
    # Compilation and execution tests on WashU Compute1 (via LSF)
    #-------------------------------------------------------------------------

    # Remove SLURM #SBATCH tags
    sed_ie '/#SBATCH -c 8/d'                   "${scriptsDir}/parallelTestCompile.sh"
    sed_ie '/#SBATCH -N 1/d'                   "${scriptsDir}/parallelTestCompile.sh"
    sed_ie '/#SBATCH -t 0-0:30/d'              "${scriptsDir}/parallelTestCompile.sh"
    sed_ie '/#SBATCH -p REQUESTED_PARTITION/d' "${scriptsDir}/parallelTestCompile.sh"
    sed_ie '/#SBATCH --mem=8000/d'             "${scriptsDir}/parallelTestCompile.sh"
    sed_ie '/#SBATCH -p REQUESTED_PARTITION/d' "${scriptsDir}/parallelTestCompile.sh"
    sed_ie '/#SBATCH --mail-type=END/d'        "${scriptsDir}/parallelTestCompile.sh"
    sed_ie '/#SBATCH -c 24/d'                  "${scriptsDir}/parallelTestExecute.sh"
    sed_ie '/#SBATCH -N 1/d'                   "${scriptsDir}/parallelTestExecute.sh"
    sed_ie '/#SBATCH -t 0-6:00/d'              "${scriptsDir}/parallelTestExecute.sh"
    sed_ie '/#SBATCH -p REQUESTED_PARTITION/d' "${scriptsDir}/parallelTestExecute.sh"
    sed_ie '/#SBATCH --mem=90000/d'            "${scriptsDir}/parallelTestExecute.sh"
    sed_ie '/#SBATCH --mail-type=END/d'        "${scriptsDir}/parallelTestExecute.sh"

    # Replace "REQUESTED_PARTITION" with the partition name
    sed_ie "${sedPartitionCmd}" "${scriptsDir}/parallelTestCompile.sh"
    sed_ie "${sedPartitionCmd}" "${scriptsDir}/parallelTestExecute.sh"

    # Submit compilation tests script
    output=$(bsub $scriptsDir}/parallelTestCompile.sh)
    output=($output)
    cmpId=${output[1]}
    cmpId=${cmpId/<}
    cmpId=${cmpId/>}

    # Submit execution tests script as a job dependency
    output=$(bsub -w "exit(${cmpId},0)" ${scriptsDir}/parallelTestExecute.sh)
    output=($output)
    exeId=${output[1]}
    exeId=${exeId/<}
    exeId=${exeId/>}

else

    #-------------------------------------------------------------------------
    # Exit with error
    #-------------------------------------------------------------------------
    echo ""
    echo "ERROR! Invalid choice of arguments!"
    echo "${usage}"
    exit 1

fi

# Change back to this directory
cd "${thisDir}"
						      
#=============================================================================
# Cleanup and quit
#=============================================================================

# Free local variables
unset cmpId
unset envFile
unset exeId
unset logsDir
unset ptRoot
unset quick
unset output
unset scheduler
unset scriptsDir
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
