#!/bin/bash                                                                                                                   

# When running on AWS spot instances you may have your instance terminated and your job restarted.
# Rather than start from the beginning, you can output frequent checkpoint files and restart
# with the most recent checkpoint date. This utility script recursively loops over all Restarts
# directories in a base path, determines the earliest and latest restart file, and deletes
# all others. This can be run periodically to keep storage size down when doing long runs.
#
# Dryrun option can be turned on by passing --dryrun when executing the script.
# This will print what files are kept and what files will be deleted, without actually
# deleting any files.
#
# ~ Lizzie Lundgren, May 2026

# Set path within which to recursively search for */Restarts
base_path="/path/to/directory/containing/gchp/run/directories"

# Dryrun is off by default (files will be deleted!)
dryrun=false

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --dryrun)
            dryrun=true
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [--dryrun]"
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 [--dryrun]"
            exit 1
            ;;
    esac
done

# Loop over Restarts directories
for dir in "$base_path"/*/Restarts; do
    if [ -d "$dir" ]; then

        echo " "
        echo "Directory $dir"

        # Get all dated restart files, starting with either gcchem or GEOSChem.Restart
        files=($(ls $dir/GEOSChem.Restart.*.nc4 $dir/gcchem_internal_checkpoint.*.nc4 2>/dev/null))

        # Exit if not more than 2                                                                                             
        if [[ ${#files[@]} -lt 2 ]]; then
            echo "Fewer than 2 files. No need to delete."
            exit 0
        fi

        # Extract all unique YYYYMMDD dates across all files and sort them                                                    
        dates=($(for f in "${files[@]}"; do
            echo "$f" | grep -oE '[0-9]{8}'
        done | sort -u))

        first_date="${dates[0]}"
        last_date="${dates[-1]}"

        echo "First date: $first_date"
        echo "Last date:  $last_date"

        # Keep any file (either prefix) matching first or last date                                                           
        for file in "${files[@]}"; do
            if [[ "$file" == *"$first_date"* || "$file" == *"$last_date"* ]]; then
                echo "--> Keeping  $file"
            else
                echo "--> Deleting $file"
		if [ "$dryrun" = false ]; then
		    rm "$file"
		fi

            fi
        done
    fi
done

exit 0
