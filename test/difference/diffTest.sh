#!/bin/bash

#EOC
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: diffTest.sh
#
# !DESCRIPTION: Looks for differences in two integration test or parallel
#  test directories.  Checks both diagnostic files (in the OutputDir/ folder
#  of each run directory) and restart files (in the Restarts/ folder of
#  each run directory).
#\\
#\\
# !CALLING SEQUENCE
#  ./diffTest.sh <ref_it_dir> <dev_it_dir>
#
#  where <ref_it_dir> is the path to the "Ref" integration test root folder
#  and   <dev_it_dir> is the path to the "Dev" integration test root folder
#
# !REMARKS:
#  TODO: Add capability to check parallel tests for identicality
#
# !AUTHORS
#  Lizzie Lundgren (@lizziel)
#  Bob Yantosca (@yantosca)
#EOP
#------------------------------------------------------------------------------
#BOC


function check_for_diffs() {

    #========================================================================
    # Compares netCDF files (GEOS-Chem diagnostics or restart files)
    # in two different integration test folders.  Prints informational
    # messages denoting if differences are found, and in which files.
    #========================================================================

    # Directories to compare
    refRoot="${1}"
    devRoot="${2}"
    refRunDir="${3}"
    devRunDir="${4}"

    # Dir where files are found ("OutputDir" or "Restarts")
    refBase=$(basename "${refRunDir}")
    devBase=$(basename "${devRunDir}")

    # Compare files in data directories
    result=$(diff -r "${refRunDir}" "${devRunDir}")

    # No differences found!  Print and exit.
    if [[ "x${result}" == "x" ]]; then
	printf "   -> No differences in ${refBase}\n"
	return 0
    fi

    # Get number of differences
    nDiffs=$(printf "$result\n" | wc -l)
    if [[ "x${nDiffs}" == "x1" ]]; then
	printf "   -> ${nDiffs} difference found in ${refBase}\n"
    else
	printf "   -> ${nDiffs} differences found in ${refBase}\n"
    fi

    # Format the result to remove extraneous characters
    for (( i=1; i<=$nDiffs; i++ )) do
	result=${result/"${refRoot}/"}
	result=${result/"${refRoot}/"}
	result=${result/"${devRoot}/"}
	result=${result/"${devRoot}/"}
	if (( $i == 1 )); then
	   result=${result/'Binary files '/'* '}
	else
	   result=${result/'Binary files '/'      * '}
	fi
	result=${result/'and'/'\n       '}
	result=${result/'differ'}
    done
    printf "      $result\n"
    return 0
}


function main() {

    #========================================================================
    # Main program.  Loops over all run directories in the <ref_it_dir>
    # and <dev_it_dir> paths, and checks diagnostic & restart files
    # for differences.
    #========================================================================

    # Error check # of arguments
    if [[ $# -ne 2 ]]; then
	echo "Usage: ./diffTest.sh <ref_it_dir> <dev_it_dir>"
	exit 1
    fi

    # Paths to Ref & Def root directories
    refRoot=$(dirname "${1}")
    devRoot=$(dirname "${2}")

    # Paths to Ref/rundirs and Dev/rundirs directories
    refRunDirs="${1}/rundirs"
    devRunDirs="${2}/rundirs"

    # Get a list of all the run directories
    runDirs=$(ls "${refRunDirs}")

    # Loop over run directories
    for dir in ${runDirs[@]}; do
	printf "Checking ${dir}\n"

	# Check diagnostic files for differences
	ref="${refRunDirs}/${dir}/OutputDir"
	dev="${devRunDirs}/${dir}/OutputDir"
	check_for_diffs "${refRoot}" "${devRoot}" "${ref}" "${dev}"

	# Check restart files for differences
	ref="${refRunDirs}/${dir}/Restarts"
	dev="${devRunDirs}/${dir}/Restarts"
	check_for_diffs "${refRoot}" "${devRoot}" "${ref}" "${dev}"

	# Print a space between rundirs
	ref=""
	dev=""
	printf "\n"
    done

    # Return with success
    return 0
}

#========================================================================
# Call main with command-line arguments
#========================================================================
main $@

