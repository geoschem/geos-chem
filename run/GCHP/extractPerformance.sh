#!/bin/bash

# extractPerformance.sh

# Prints summary of timing information based on timers in allPEs.log, including model throughput,
# total run time, total times for SetServices, Initialize, and Finalize, and Run time per
# gridded component (ExtData, GCHP and its child grid comps (GCHPctmEnv, DYNAMICS, GEOSchem),
# and History).
#
# Optional argument: path to log file allPEs.log; default is local directory
#
# Example calls:
#   $ path/to/extractPerformance.sh
#   $ path/to/extractPerformance.sh /path/to/allPEs.log

# If a directory is supplied, navigate to it
if [[ "x${1}" != "x" ]]; then
    cd ${1}
fi
curdir=$(pwd -P)
logFile="$curdir/allPEs.log"
if [ ! -f "$logFile" ]; then
    echo "Error: File 'allPEs.log' not found"
    exit 1
fi

# Find the line number right before the timers summary in the log
start_line=$(grep -n "Report on process: 0" "$logFile" | head -n 1 | cut -d: -f1)
if [ -z "$start_line" ]; then
    echo "Error: Could not locate 'Report on process: 0' in allPEs.log!"
    exit 1
fi

# Get the times from the entries starting at "ALL"
total_time_val=""
raw_all_line=$(tail -n +"$start_line" "$logFile" | grep -m 1 "INFO: All")
if [ -n "$raw_all_line" ]; then
    read -ra all_tokens <<< "$raw_all_line"
    for i in "${!all_tokens[@]}"; do
        if [[ "${all_tokens[i]}" == "INFO:" ]]; then
            total_time_val="${all_tokens[((i + 3))]}"
            break
        fi
    done
fi

# Get the Model Throughput value
throughput_val=""
raw_throughput=$(grep -m 1 "Model Throughput:" "$logFile")
if [ -n "$raw_throughput" ]; then
    throughput_val=$(echo "$raw_throughput" | awk '{print $(NF-3)}')
fi

# Print the total time ("ALL" in the file) and the model throughput
echo "=========================================================="
if [ -n "$total_time_val" ]; then
    printf "Total time:       %10s seconds\n" "$total_time_val"
else
    echo "Total time:              Not found"
fi
if [ -n "$throughput_val" ]; then
    printf "Model Throughput: %10s days per day\n" "$throughput_val"
else
    echo "Model Throughput:        Not found"
fi
echo "=========================================================="
echo ""

# Print the header for the timers summary
echo "GCHP Timing Summary By Lifecycle Phase"
echo "--------------------------------------------------------"
printf "%-15s | %-20s\n" "Component Name" "Inclusive Time T (sec)"
echo "--------------------------------------------------------"

# Clear any temporary files created from a previous run, just in case of early termination and no cleanup
rm -f .*.tmp

current_phase=""
skip_subsets=false

# Helper function to strip two dashes from a string, to clean up the look of the timers table
clean_name() {
    local name="$1"
    if [[ "$name" == --* ]]; then
        echo "${name#--}"
    else
        echo "$name"
    fi
}

# Read the timer lines.
tail -n +"$start_line" "$logFile" | while read -r line; do
    if [[ "$line" == *"Times for component"* ]]; then
        break
    fi
    
    if [[ "$line" =~ "INFO: All" || "$line" =~ "INFO: --" ]]; then
        read -ra tokens <<< "$line"
        
        info_idx=-1
        for i in "${!tokens[@]}"; do
            if [[ "${tokens[i]}" == "INFO:" ]]; then
                info_idx=$i
                break
            fi
        done

        # Print strings to temporary files, since will be assembled into table in order of GCHP execution
        if [ $info_idx -ne -1 ]; then
            comp_name="${tokens[((info_idx + 1))]}"
            time_val="${tokens[((info_idx + 3))]}"
            
            # Specialized handling for certain sections, e.g. only print subset timers for Run phase
            if [[ "$comp_name" == "--SetService" ]]; then
                current_phase="SetService"
                skip_subsets=true
                printf "%-15s | %-20s\n" "$(clean_name "$comp_name")" "$time_val" >> .tree_main.tmp
                continue
            elif [[ "$comp_name" == "--Initialize" ]]; then
                current_phase="Initialize"
                skip_subsets=true
                printf "%-15s | %-20s\n" "$(clean_name "$comp_name")" "$time_val" >> .tree_main.tmp
                continue
            elif [[ "$comp_name" == "--Run" ]]; then
                current_phase="Run"
                skip_subsets=false
                echo "INJECT_RUN_SUBSECTIONS" >> .tree_main.tmp
                printf "%-15s | %-20s\n" "$(clean_name "$comp_name")" "$time_val" >> .tree_main.tmp
                continue
            elif [[ "$comp_name" == "--Finalize" ]]; then
                current_phase="Finalize"
                skip_subsets=true  # Drop all Finalize subsections entirely
                printf "%-15s | %-20s\n" "$(clean_name "$comp_name")" "$time_val" >> .tree_main.tmp
                continue
            elif [[ "$comp_name" == "All" ]]; then
                continue
            fi

            # Skip deeper nested components in addition to the ones flagged above to skip

            if $skip_subsets && [[ "$comp_name" =~ ^---- ]]; then
                continue
            fi
            
            adjusted_name=$(clean_name "$comp_name")
            formatted_line=$(printf "%-15s | %-20s" "$adjusted_name" "$time_val")

            if [[ "$current_phase" == "Run" ]]; then
                if [[ "$comp_name" == *EXTDATA* ]]; then
                    echo "$formatted_line" >> .run_extdata.tmp
                elif [[ "$comp_name" == "----GCHP" ]]; then
                    echo "$formatted_line" >> .run_gchp_parent.tmp
                elif [[ "$comp_name" == *GCHPctmEnv* ]]; then
                    echo "$formatted_line" >> .run_gchp_env.tmp
                elif [[ "$comp_name" == *DYNAMICS* ]]; then
                    echo "$formatted_line" >> .run_dynamics.tmp
                elif [[ "$comp_name" == *GCHPchem* ]]; then
                    echo "$formatted_line" >> .run_gchpchem.tmp
                elif [[ "$comp_name" == *HIST* ]]; then
                    echo "$formatted_line" >> .run_hist.tmp
                fi
            else
                echo "$formatted_line" >> .tree_main.tmp
            fi
        fi
    fi
done

# Write the table in a specific order to match GCHP execution
while read -r line; do
    if [[ "$line" == "INJECT_RUN_SUBSECTIONS" ]]; then
        read -r header_line; echo "$header_line"
        [ -f .run_extdata.tmp ] && cat .run_extdata.tmp
        [ -f .run_gchp_parent.tmp ] && cat .run_gchp_parent.tmp
        [ -f .run_gchp_env.tmp ] && cat .run_gchp_env.tmp
        [ -f .run_dynamics.tmp ] && cat .run_dynamics.tmp
        [ -f .run_gchpchem.tmp ] && cat .run_gchpchem.tmp
        [ -f .run_hist.tmp ] && cat .run_hist.tmp
        continue
    fi

    echo "$line"
done < .tree_main.tmp

# Remove temporary files
rm -f .*.tmp

echo "--------------------------------------------------------"
echo "See log file allPEs.log for more detailed timing"
