#!/bin/bash

#EOC
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: integrationDiffTest.sh
#
# !DESCRIPTION: Looks for differences in two integration test or parallel
#  test directories.  Checks both diagnostic files (in the OutputDir/ folder
#  of each run directory) and restart files (in the Restarts/ folder of
#  each run directory).
#\\
#\\
# !CALLING SEQUENCE
#  ./integrationDiffTest.sh <ref_it_dir> <dev_it_dir>
#
#  where <ref_it_dir> is the path to the "Ref" integration test root folder
#  and   <dev_it_dir> is the path to the "Dev" integration test root folder
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
    rootDir="${1}"
    refDir="${2}"
    devDir="${3}"

    # Dir where files are found ("OutputDir" or "Restarts")
    refBase=$(basename "${refDir}")
    devBase=$(basename "${devDir}")

    # Compare files in data directories
    result=$(diff -r "${refDir}" "${devDir}")

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
	result=${result/"${rootDir}/"}
	result=${result/"${rootDir}/"}
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
	echo "Usage: ./intDiffTest <ref_it_dir> <dev_it_dir>"
	exit 1
    fi

    # Integration test root directory
    rootDir=$(pwd)

    # Ref & Dev root directories
    refRoot="${rootDir}/${1}/rundirs"
    devRoot="${rootDir}/${2}/rundirs"

    # List of run directories
    runDirs=$(ls "${refRoot}")

    # Loop over run directories
    for dir in ${runDirs[@]}; do
	printf "Checking ${dir}\n"

	# Check diagnostic files for differences
	ref="${refRoot}/${dir}/OutputDir"
	dev="${devRoot}/${dir}/OutputDir"
	check_for_diffs "${rootDir}" "${ref}" "${dev}"

	# Check restart files for differences
	ref="${refRoot}/${dir}/Restarts"
	dev="${devRoot}/${dir}/Restarts"
	check_for_diffs "${rootDir}" "${ref}" "${dev}"

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
