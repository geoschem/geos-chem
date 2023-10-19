#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: integrationTest.sh
#
# !DESCRIPTION: Runs integration tests on the various GCHP run directories
# (interactively, or with a scheduler).
#\\
#\\
# !CALLING SEQUENCE:
#  ./integrationTest.sh -d root-dir -e env-file [-h] [-p partition] [-q] [-s scheduler]
#
#  Where the command-line arguments are as follows:
#
#    -d root-dir  : Specify the root folder for integration tests
#    -e env-file  : Specify the environment file (w/ module loads)
#    -h           : Display a help message
#    -n           : Do not bootstrap missing restart file variables
#    -p partition : Select partition for SLURM or LSF schedulers
#    -q           : Run a quick set of integration tests (for testing)
#    -s scheduler : Specify the scheduler (SLURM or LSF)
#
#  NOTE: you can also use the following long name options:
#
#    --directory root-dir  (instead of -d root-dir )
#    --env-file  env-file  (instead of -e env-file )
#    --help                (instead of -h          )
#    --no-bootstrap        (instead of -n          )
#    --partition partition (instead of -p partition)
#    --quick               (instead of -q          )
#    --scheduler scheduler (instead of -s scheduler)
#EOP
#------------------------------------------------------------------------------
#BOC

#=============================================================================
# Initialize
#=============================================================================
this="$(basename ${0})"
usage="Usage: ${this} -d root-dir -e env-file [-h] [-p partition] [-q] [-s scheduler]"
itRoot="none"
envFile="none"
scheduler="none"
sedPartitionCmd="none"
quick="no"
bootStrap="yes"

#=============================================================================
# Parse command-line arguments
# See https://www.baeldung.com/linux/bash-parse-command-line-arguments
#=============================================================================

# Call Linux getopt function to specify short & long input options
# (e.g. -d or --directory, etc).  Exit if not succesful
validArgs=$(getopt --options d:e:hnp:qs: \
  --long directory:,env-file:,help,no-bootstrap,partition:,quick,scheduler: -- "$@")
if [[ $? -ne 0 ]]; then
    exit 1;
fi

# Parse arguments and set variables accordingly
# TODO: replace "scheduler" with "site" argument; set scheduler from site.
eval set -- "${validArgs}"
while [ : ]; do
    case "${1}" in

        # -d or --directory specifies the root folder for tests
        -d | --directory)
            itRoot="${2}"
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

	# -n or --no-bootstrap prevents bootstrapping missing variables in
        # restart files (i.e. do not change EFYO -> CYS in HEMCO_Config.rc)
	-n | --no-bootstrap)
            bootStrap="no"
	    shift
            ;;

	# -p or --partition replaces REQUESTED_PARTITION with the user's choice
	-p | --partition)
            sedPartitionCmd="s/REQUESTED_PARTITION/${2}/"
            shift 2
            ;;

	# -q or --quick runs a quick set of integration tests (for testing)
	-q | --quick)
            quick="yes"
            shift
	    ;;

	# -s or --scheduler selects the scheduler (case-insensitive)
	-s | --scheduler)
	    scheduler="${2^^}"
            shift 2
            ;;

	--) shift;
            break
            ;;
    esac
done

# Error check integration tests root path
if [[ "x${itRoot}" == "xnone" ]]; then
    echo "ERROR: The integration test root directory has not been specified!"
    echo "${usage}"
    exit 1
fi

# Error check environment file
# TODO: Add a test on site name rather than scheduler
if [[ "x${scheduler}" != "xLSF" ]]; then
    if [[ "x${envFile}" == "xnone" ]]; then
        echo "ERROR: The enviroment file (module loads) has not been specified!"
        echo "${usage}"
        exit 1
    fi
fi

# Exit if no partition has been selected for SLURM
if [[ "x${scheduler}" == "xSLURM" && "x${sedPartitionCmd}" == "xnone" ]]; then
    echo "ERROR: You must specify a partition for SLURM."
    echo "${usage}"
    exit 1
fi

# Exit if no partition has been selected for SLURM
if [[ "x${scheduler}" == "xLSF" && "x${sedPartitionCmd}" == "xnone" ]]; then
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
if [[ -f "../../shared/commonFunctionsForTests.sh" ]]; then
    . "${thisDir}/../../shared/commonFunctionsForTests.sh"
elif [[ -f "${thisDir}/commonFunctionsForTests.sh" ]]; then
    . "${thisDir}/commonFunctionsForTests.sh"
fi

#=============================================================================
# Create integration test directories in the root folder
#=============================================================================

# Convert integration test root folder to an absolute path
itRoot=$(absolute_path "${itRoot}")

# Prevent running integration tests in the source code directory tree
if [[ "$(absolute_path ${thisDir})" =~ "${itRoot}" ]]; then
    echo "ERROR: You cannot run integration tests in the source code directory!"
    exit 1
fi
# Create GEOS-Chem run directories in the integration test root folder
./integrationTestCreate.sh "${itRoot}" "${envFile}" "${quick}"
if [[ $? -ne 0 ]]; then
   echo "ERROR: Could not create integration test run directories!"
   exit 1
fi

# Navigate to the root test folder
if [[ -d "${itRoot}" ]]; then
    cd "${itRoot}"
else
    echo "ERROR: ${itRoot} is not a valid directory!  Exiting..."
    exit 1
fi

# Define local convenience variables
logsDir="${itRoot}/${LOGS_DIR}"
scriptsDir="${itRoot}/${SCRIPTS_DIR}"
rundirsDir="${itRoot}/${RUNDIRS_DIR}"

# Edit setCommonRunSettingss.sh scripts to enable or disable bootstrapping
# (i.e. to allow missing species in restart files or not)
gchp_enable_or_disable_bootstrap "${bootStrap}" "${rundirsDir}"

# Navigate to the logs directory (so all output will be placed there)
cd "${logsDir}"

#=============================================================================
# Compile the code and run the integration tests
#=============================================================================
if [[ "x${scheduler}" == "xSLURM" ]]; then

    #-------------------------------------------------------------------------
    # Integration tests will run via SLURM
    #-------------------------------------------------------------------------

    # Remove LSF #BSUB tags
    sed_ie '/#BSUB -q REQUESTED_PARTITION/d' "${scriptsDir}/integrationTestCompile.sh"
    sed_ie '/#BSUB -n 8/d'                   "${scriptsDir}/integrationTestCompile.sh"
    sed_ie '/#BSUB -W 1:30/d'                "${scriptsDir}/integrationTestCompile.sh"
    sed_ie '/#BSUB -o lsf-%J.txt/d'          "${scriptsDir}/integrationTestCompile.sh"
    sed_ie \
	'/#BSUB -R "rusage\[mem=8GB\] span\[ptile=1\] select\[mem < 1TB\]"/d' \
	"${scriptsDir}/integrationTestCompile.sh"
    sed_ie \
	"/#BSUB -a 'docker(registry\.gsc\.wustl\.edu\/sleong\/esm\:intel\-2021\.1\.2)'/d" \
	"${scriptsDir}/integrationTestCompile.sh"
    sed_ie '/#BSUB -q REQUESTED_PARTITION/d' "${scriptsDir}/integrationTestExecute.sh"
    sed_ie '/#BSUB -n 24/d'                  "${scriptsDir}/integrationTestExecute.sh"
    sed_ie '/#BSUB -W 3:30/d'                "${scriptsDir}/integrationTestExecute.sh"
    sed_ie '/#BSUB -o lsf-%J.txt/d'          "${scriptsDir}/integrationTestExecute.sh"
    sed_ie \
	'/#BSUB -R "rusage\[mem=90GB\] span\[ptile=1\] select\[mem < 2TB\]"/d' \
	"${scriptsDir}/integrationTestExecute.sh"
    sed_ie \
	"/#BSUB -a 'docker(registry\.gsc\.wustl\.edu\/sleong\/esm\:intel\-2021\.1\.2)'/d" \
	"${scriptsDir}/integrationTestExecute.sh"

    # Replace "REQUESTED_PARTITION" with the partition name
    sed_ie "${sedPartitionCmd}" "${scriptsDir}/integrationTestCompile.sh"
    sed_ie "${sedPartitionCmd}" "${scriptsDir}/integrationTestExecute.sh"

    # Submit compilation tests script
    output=$(sbatch ${scriptsDir}/integrationTestCompile.sh)
    output=($output)
    cmpId=${output[3]}

    # Submit execution tests script as a job dependency
    output=$(sbatch --dependency=afterok:${cmpId} ${scriptsDir}/integrationTestExecute.sh)
    output=($output)
    exeId=${output[3]}

    echo ""
    echo "Compilation tests submitted as SLURM job ${cmpId}"
    echo "Execution   tests submitted as SLURM job ${exeId}"

elif [[ "x${scheduler}" == "xLSF" ]]; then

    #-------------------------------------------------------------------------
    # Integration tests will run via LSF
    #-------------------------------------------------------------------------

    # Remove SLURM #SBATCH tags
    sed_ie '/#SBATCH -c 8/d'                   "${scriptsDir}/integrationTestCompile.sh"
    sed_ie '/#SBATCH -N 1/d'                   "${scriptsDir}/integrationTestCompile.sh"
    sed_ie '/#SBATCH -t 0-1:30/d'              "${scriptsDir}/integrationTestCompile.sh"
    sed_ie '/#SBATCH -p REQUESTED_PARTITION/d' "${scriptsDir}/integrationTestCompile.sh"
    sed_ie '/#SBATCH --mem=8000/d'             "${scriptsDir}/integrationTestCompile.sh"
    sed_ie '/#SBATCH -p REQUESTED_PARTITION/d' "${scriptsDir}/integrationTestCompile.sh"
    sed_ie '/#SBATCH --mail-type=END/d'        "${scriptsDir}/integrationTestCompile.sh"
    sed_ie '/#SBATCH -c 24/d'                  "${scriptsDir}/integrationTestExecute.sh"
    sed_ie '/#SBATCH -N 1/d'                   "${scriptsDir}/integrationTestExecute.sh"
    sed_ie '/#SBATCH -t 0-3:30/d'              "${scriptsDir}/integrationTestExecute.sh"
    sed_ie '/#SBATCH -p REQUESTED_PARTITION/d' "${scriptsDir}/integrationTestExecute.sh"
    sed_ie '/#SBATCH --mem=90000/d'            "${scriptsDir}/integrationTestExecute.sh"
    sed_ie '/#SBATCH --mail-type=END/d'        "${scriptsDir}/integrationTestExecute.sh"

    # Replace "REQUESTED_PARTITION" with the partition name
    sed_ie "${sedPartitionCmd}" "${scriptsDir}/integrationTestCompile.sh"
    sed_ie "${sedPartitionCmd}" "${scriptsDir}/integrationTestExecute.sh"

    # Submit compilation tests script
    output=$(bsub ${scriptsDir}/integrationTestCompile.sh)
    output=($output)
    cmpId=${output[1]}
    cmpId=${cmpId/<}
    cmpId=${cmpId/>}

    # Submit execution tests script as a job dependency
    output=$(bsub -w "exit(${cmpId},0)" ${scriptsDir}/integrationTestExecute.sh)
    output=($output)
    exeId=${output[1]}
    exeId=${exeId/<}
    exeId=${exeId/>}

    echo ""
    echo "Compilation tests submitted as LSF job ${cmpId}"
    echo "Execution   tests submitted as LSF job ${exeId}"

else

    #-------------------------------------------------------------------------
    # Integration tests will run interactively
    #-------------------------------------------------------------------------

    # Run compilation tests
    echo ""
    echo "Compiliation tests are running..."
    ${scriptsDir}/integrationTestCompile.sh &

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
unset itRoot
unset logsDir
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
