#!/bin/bash
DESCRIPTION="""Initialize the current working directory as a run directory (minimal file set). 

usage: 
   init_rd.sh FILE...
      FILE...: Files defining variables used to initialize run directory. The
               rightmost file has the highest precedence.
               
               Each FILE is sourced and the loaded variables are used to
               fill the run directory templates using envsubst.

   init_rd.sh --rundir-vars
      Lists all the RUNDIR variables found in the templates.

   init_rd.sh [-h] [--help]
      Prints this help message.
"""
if [[ ( $* == --help ) ||  ( $* == -h ) || $# -eq 0 ]]; then 
   echo "$DESCRIPTION"
   exit 0
fi 

set -e

THIS_SCRIPTS_DIRECTORY=$(realpath $(dirname "$0"))

if [[ ( $* == --rundir-vars ) ]]; then 
   grep -roh 'RUNDIR_[A-Z_][A-Z_]*' $THIS_SCRIPTS_DIRECTORY | grep -v 'RUNDIR_VARS' |  sort | uniq
   exit 0
fi 

function get_rundir_vars_list() {
   sed -n 's#^\s*\([A-Za-z0-9_][A-Za-z0-9_]*\)=.*#\1#p'
}

# Source given files, and build variable list
variables=
for envfile in "$@"; do
   source $envfile
   variables+="$(cat $envfile | get_rundir_vars_list) "
   export $variables
done
# Source a second time to resolve dependent variables
for envfile in "$@"; do
   source $envfile
   export $variables
done
variables=$(echo $variables | sort | uniq)
envsubst_list="$(printf '${%s} ' $variables)"

COPY_LIST="""
HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.${RUNDIR_SIM_NAME}
"""

SUBST_LIST="""
geoschem_config.yml.templates/geoschem_config.yml.${RUNDIR_SIM_NAME}
HEMCO_Config.rc.templates/HEMCO_Config.rc.${RUNDIR_SIM_NAME}
HISTORY.rc.templates/HISTORY.rc.${RUNDIR_SIM_NAME}
HEMCO_Config.rc.templates/${RUNDIR_MET_FIELD_CONFIG}
"""

function filename_with_suffixes_removed() {
   basename $(basename $1 .${RUNDIR_SIM_NAME}) .template
}

# Copy files in COPY_LIST to cwd
for fpath in $COPY_LIST; do
   cp $THIS_SCRIPTS_DIRECTORY/$fpath $(filename_with_suffixes_removed $fpath)
done

# Make OutputDir
mkdir -p OutputDir

# Copy and make substitutions for each file in SUBST_LIST
for fpath in $SUBST_LIST; do
   envsubst "$envsubst_list" < $THIS_SCRIPTS_DIRECTORY/$fpath > $(filename_with_suffixes_removed $fpath)
done

# Make links
unlink runScriptSamples
ln -s $THIS_SCRIPTS_DIRECTORY/runScriptSamples runScriptSamples
