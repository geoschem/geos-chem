#!/usr/bin/env bash
set -euo pipefail

# Helper script for checking input file availability on local machine based on 
# configured settings in extdata.yaml.
# By default, only warns for any missing input files in the simulation date range. 

# Output format:
# SUMMARY: [freq] collection_path: X / Y files missing (start_date to end)
#
# [freq] is the frequency files are expected in the provided folder, parsed from 
# the smallest date token (eg %d2) in the template file path. 
# note this is NOT the same as the sampling frequency as defined in yaml samplings. 
# 
# collection_path is the path where files are expected to be found, derived from 
# the collection template in the extdata.yaml.
#
# X is the number of missing files, Y is the total number of files expected based
# on the simulation date range and the file frequency. Minimum number is 1. 
#

#############################################################
# User config (you should not need to edit after this section)
#############################################################

# This setting determines whether to check for files on disk within a provided valid range
# for a given collection. This can be slow, especially for daily collections like MERRA2, 
# so it is disabled by default. 
# This is useful because valid_range currently fails if any files in valid range are missing,
# regardless if they are within simulation date range or not. 
CHECK_VALID_RANGE=false

# Lightweight check to determine if the provided simulation date is within the valid_range 
# for collections. Does not check disk for files, just compares simulation settings to yaml
# text. Independent of CHECK_VALID_RANGE.
CHECK_SIM_IN_VALID_RANGE=true


# When true, print file availability summary for all collections, 
# even those with 0 missing files 
DEBUG_ALL=false           

# Format: YYYYMMDD
SIM_START_DATE=20190701   
SIM_END_DATE=20190801

# Path to your YAML file relative to this script
YAML=extdata.yaml        

##########################################################
# Derived config
##########################################################

SIM_START_ISO="${SIM_START_DATE:0:4}-${SIM_START_DATE:4:2}-${SIM_START_DATE:6:2}"
SIM_END_ISO="${SIM_END_DATE:0:4}-${SIM_END_DATE:4:2}-${SIM_END_DATE:6:2}"

##########################################################
# Early failure checks
##########################################################
if [[ ! -f ${YAML} ]]; then
   echo "Could not find configuration file: ${YAML}!  Exiting ..."
   exit 1
fi

if (( $(awk '/^[[:space:]][[:space:]][^[:space:]]+:[[:space:]]*(#.*)?$/' "${YAML}" | wc -l) == 0 )); then
   echo "ERROR: no collections parsed from ${YAML} — check indentation/format" >&2
   exit 1
fi

#------------------------
# Config info
#------------------------

echo "INFO: Simulation date range: $SIM_START_ISO to $SIM_END_ISO"

if ! $DEBUG_ALL; then
    echo "INFO: DEBUG_ALL is false, only printing summaries for collections with missing files"
else
    echo "INFO: DEBUG_ALL is true, printing summaries for all collections (including those with 0 missing files)"
fi

if $CHECK_VALID_RANGE; then
    echo "INFO: Checking collections with valid_range (using YAML dates)"
else
    echo "INFO: Skipping collections with provided valid_range"
fi

########################################
# AWK: parse Exports -> collection_name<TAB>freq
########################################

AWK_COL_WITH_VR=$(cat << 'EOF'
function emit() {
    if (current_name != "" && template != "" && saw_valid) {
        print current_name "\t" start_date "\t" end_date "\t" template
    }
}

# 2-space-indented collection name -> flush previous, start new
/^[[:space:]][[:space:]][^[:space:]]+:[[:space:]]*(#.*)?$/ {
    emit()
    current_name = $1
    sub(/:$/, "", current_name)
    start_date = ""
    end_date   = ""
    template   = ""
    saw_valid  = 0
}

/^[[:space:]][[:space:]][[:space:]][[:space:]]valid_range:[[:space:]]*/ {
    line = $0
    sub(/^[[:space:]][[:space:]][[:space:]][[:space:]]valid_range:[[:space:]]*/, "", line)
    gsub(/^"[[:space:]]*|"[[:space:]]*$/, "", line)
    gsub(/^'[[:space:]]*|'[[:space:]]*$/, "", line)

    n = split(line, a, "/")
    if (n >= 1) {
        start_date = a[1]
        sub(/T.*/, "", start_date)
    }
    if (n >= 2) {
        end_date = a[2]
        sub(/T.*/, "", end_date)
    }
    saw_valid = 1
}

/^[[:space:]][[:space:]][[:space:]][[:space:]]template:[[:space:]]*/ {
    line = $0
    sub(/^[[:space:]][[:space:]][[:space:]][[:space:]]template:[[:space:]]*/, "", line)
    gsub(/^"[[:space:]]*|"[[:space:]]*$/, "", line)
    gsub(/^'[[:space:]]*|'[[:space:]]*$/, "", line)
    template = line
}

END { emit() }
EOF
)

AWK_COL_NO_VR=$(cat << 'EOF'
function emit() {
    if (current_name != "" && template != "" && !saw_valid) {
        print current_name "\t" template   # name<TAB>template
    }
}

# 2-space-indented collection name -> flush previous, start new
/^[[:space:]][[:space:]][^[:space:]]+:[[:space:]]*(#.*)?$/ {
    emit()
    current_name = $1
    sub(/:$/, "", current_name)
    template  = ""
    saw_valid = 0
}

/^[[:space:]][[:space:]][[:space:]][[:space:]]valid_range:[[:space:]]*/ {
    saw_valid = 1
}

/^[[:space:]][[:space:]][[:space:]][[:space:]]template:[[:space:]]*/ {
    line = $0
    sub(/^[[:space:]][[:space:]][[:space:]][[:space:]]template:[[:space:]]*/, "", line)
    gsub(/^"[[:space:]]*|"[[:space:]]*$/, "", line)
    gsub(/^'[[:space:]]*|'[[:space:]]*$/, "", line)
    template = line
}

END { emit() }
EOF
)


########################################
# Helper: infer file frequency from template tokens
########################################

infer_freq_from_template() {
    local template="$1"
    local freq

    if [[ "$template" == *%d2* ]]; then
        freq="daily"      # has day component
    elif [[ "$template" == *%m2* ]]; then
        freq="monthly"    # has month, but not day
    elif [[ "$template" == *%y4* ]]; then
        freq="annual"     # has only year
    else
        freq="constant"   # no date tokens at all
    fi

    printf '%s\n' "$freq"
}

########################################
# Helper: check sim range against a collection's valid_range
########################################

check_sim_in_valid_range() {
    local name="$1"       # collection name
    local vr_start="$2"   # valid_range start, YYYY-MM-DD (may be empty)
    local vr_end="$3"     # valid_range end,   YYYY-MM-DD (may be empty)

    local msg=""

    # Lower bound: sim starts before data begins
    if [[ -n "$vr_start" && "$SIM_START_ISO" < "$vr_start" ]]; then
        msg="sim start $SIM_START_ISO precedes valid_range start $vr_start"
    fi

    # Upper bound: sim ends after data ends (strict > so an exact match is OK)
    if [[ -n "$vr_end" && "$SIM_END_ISO" > "$vr_end" ]]; then
        [[ -n "$msg" ]] && msg="$msg; "
        msg="${msg}sim end $SIM_END_ISO exceeds valid_range end $vr_end"
    fi

    if [[ -n "$msg" ]]; then
        echo "WARNING: [valid_range] $name: $msg (valid_range ${vr_start:-open} to ${vr_end:-open})"
    elif [[ "$DEBUG_ALL" == "true" ]]; then
        echo "SUMMARY: [valid_range] $name: sim range within valid_range (${vr_start:-open} to ${vr_end:-open})"
    fi
}

########################################
# Helper: check range with frequency
########################################

check_range_with_freq() {
    # Avoid killing entire script if a single date -d fails
    set +e

    local name="$1"
    local d_start="$2"   # YYYY-MM-DD
    local d_end="$3"     # YYYY-MM-DD
    local template="$4"
    local freq="${5:-}"

    [[ -z "$freq" ]] && freq="daily"

    local missing_count=0
    local checked_count=0

    case "$freq" in
        constant)
            # For constant, we just test whether the template path exists once.
            # If it has %y4/%m2/%d2, substitute d_start; otherwise, use it literally.

            local path

            if [[ "$template" == *%y4* || "$template" == *%m2* || "$template" == *%d2* ]]; then
                local y4 m2 d2
                y4=$(date -d "$d_start" +%Y)
                m2=$(date -d "$d_start" +%m)
                d2=$(date -d "$d_start" +%d)
                path=${template//%y4/$y4}
                path=${path//%m2/$m2}
                path=${path//%d2/$d2}
            else
                path="$template"
            fi

            # if [[ -f "$path" ]]; then
            #     echo "SUMMARY [constant] $name: constant file present ($path)"
            # else
            #     echo "SUMMARY [constant] $name: constant file missing ($path)"
            # fi
            ;;

        monthly)
            local d
            d=$(date -d "$d_start" +%Y-%m-01)
            
            while [[ "$d" < "$d_end" || "$d" == "$d_end" ]]; do
                local y4 m2 d2 path
                y4=$(date -d "$d" +%Y)
                m2=$(date -d "$d" +%m)
                d2=$(date -d "$d" +%d)
                path=${template//%y4/$y4}
                path=${path//%m2/$m2}
                path=${path//%d2/$d2}

                ((checked_count++))
                if [[ ! -f "$path" ]]; then
                    ((missing_count++))
                #else
                    #echo "DEBUG: found monthly file for $name: $path" >&2
                fi

                d=$(date -d "$d +1 month" +%F)
            done
            ;;

        annual)
            # Check one annual file per year overlapping the simulation window.
            # For a range 2019-07-01..2019-08-01, we still check 2019-01-01.

            local year_start year_end y d
            year_start=$(date -d "$d_start" +%Y)
            year_end=$(date -d "$d_end" +%Y)

            for ((y = year_start; y <= year_end; y++)); do
                d="${y}-01-01"
                local y4 m2 d2 path
                y4=$(date -d "$d" +%Y)
                m2=$(date -d "$d" +%m)
                d2=$(date -d "$d" +%d)
                path=${template//%y4/$y4}
                path=${path//%m2/$m2}
                path=${path//%d2/$d2}

                ((checked_count++))
                if [[ ! -f "$path" ]]; then
                    ((missing_count++))
                fi
            done
            ;;

        daily|hourly|*)
            local d="$d_start"
            while [[ "$d" < "$d_end" || "$d" == "$d_end" ]]; do
                local y4 m2 d2 path
                y4=$(date -d "$d" +%Y)
                m2=$(date -d "$d" +%m)
                d2=$(date -d "$d" +%d)
                path=${template//%y4/$y4}
                path=${path//%m2/$m2}
                path=${path//%d2/$d2}

                ((checked_count++))
                if [[ ! -f "$path" ]]; then
                    ((missing_count++))
                #else
                    #echo "DEBUG: found $freq file for $name: $path" >&2
                fi

                d=$(date -d "$d +1 day" +%F)
            done
            ;;
    esac


    if [[ "$DEBUG_ALL" = "true" ]]; then
        if (( checked_count > 0 )); then
            echo "FILE COUNT: [$freq] $name $missing_count / $checked_count files missing ($d_start to $d_end)"
        fi
    fi
    
    if [[ "$DEBUG_ALL" = "false" ]]; then
        if (( missing_count > 0 )); then
            echo "FILE COUNT: [$freq] $template $missing_count / $checked_count files missing ($d_start to $d_end)"
        fi
    fi

    #if [[ "$freq" != "constant" && $checked_count -gt 0 ]]; then
    #    echo "SUMMARY [$freq] $name: $missing_count / $checked_count files missing ($d_start to $d_end)"
    #fi

    set -e
}

########################################
# Pass: sim range vs valid_range
########################################

if $CHECK_SIM_IN_VALID_RANGE; then
    echo "INFO: Checking simulation range against each collection's valid_range"
    awk "$AWK_COL_WITH_VR" "$YAML" | \
    while IFS=$'\t' read -r name start_date end_date template; do
        [[ -z "$name" ]] && continue
        check_sim_in_valid_range "$name" "$start_date" "$end_date"
    done
fi

########################################
# Run passes using a cached collections file
########################################

if $CHECK_VALID_RANGE; then
    awk "$AWK_COL_WITH_VR" "$YAML" | \
    while IFS=$'\t' read -r name start_date end_date template; do
        [[ -z "$name" || -z "$template" ]] && continue

        freq=$(infer_freq_from_template "$template")

        check_range_with_freq "$name" "$start_date" "$end_date" "$template" "$freq"
    done
else
    awk "$AWK_COL_NO_VR" "$YAML" | \
    while IFS=$'\t' read -r name template; do
        [[ -z "$name" || -z "$template" ]] && continue

        freq=$(infer_freq_from_template "$template")

        #echo "PASS2: name='$name' freq='$freq' template='$template'" >&2

        check_range_with_freq "$name" "$SIM_START_ISO" "$SIM_END_ISO" "$template" "$freq"
    done
fi


