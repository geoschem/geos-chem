#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: newUserRegistration.sh
#
# !DESCRIPTION: Defines utility functions for first-time GEOS-Chem
#  user registration.  This code has been been abstracted out of
#  run/GCClassic/createRunDir.sh and run/GCHP/createRundir.sh.
#\\
#\\
# !REVISION HISTORY:
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

####=========================================================================
#### Common variables
####==========================================================================

# User prompt for reading in data
USER_PROMPT=">>> "

####=========================================================================
#### Send registration details to AWS database via curl
####==========================================================================
function postRegistration() {
    #### Sends registration details to API ####
    curl --location --request POST "https://gc-dashboard.org/registration" \
        --header "Content-Type: text/plain" \
        -d "{
            \"email\":             \"${1}\",
            \"name\":              \"${2}\",
            \"affiliation\":       \"${3}\",
            \"name_of_pi\":        \"${4}\",
            \"site\":              \"${5}\",
            \"git_username\":      \"${6/\@/}\",
            \"research_interest\": \"${7}\",
            \"model_type\":        \"${8}\",
            \"env_type\":          \"${9}\"
        }"
}

####=========================================================================
#### Read an entry from the command line
####==========================================================================
function userInput() {

    # Keep asking for the input until user enters a non-blank input
    while [[ -z "${val}" ]]; do
        IFS='\n' read -r -p "${USER_PROMPT}" val
    done

    # Return the answer
    printf "${val}"
}

####=========================================================================
#### Query user to provide registration information
####
#### NOTE: See https://www.baeldung.com/linux/ifs-shell-variable
#### for a description of what the IFS variable below does
####==========================================================================
function registerNewUser() {

    # Pass the model type from the GCClassic or GCHP createRunDir.sh script
    model_type="${1}"

    # Ask user several questions
    printf "\nInitiating User Registration:\n"
    printf "You will only need to fill this information out once.\n"
    printf "Please respond to all questions.\n"

    printf "${thinline}What is your name?${thinline}"
    name=$(userInput)

    # Ask for email
    printf "${thinline}What is your email address?${thinline}"
    email=$(userInput)

    printf "${thinline}What is the name of your research institution?${thinline}"
    institution=$(userInput)

    printf "${thinline}What is the name of your principal invesigator?\n"
    printf "(Enter 'self' if you are the principal investigator.)${thinline}"
    name_of_pi=$(userInput)

    printf "${thinline}Please provide the web site for your institution\n"
    printf "(group website, company website, etc.)?${thinline}"
    site=$(userInput)

    printf "${thinline}Please provide your github username (if any) so that we\n"
    printf "can recognize you in submitted issues and pull requests.${thinline}"
    git_username=$(userInput)

    printf "${thinline}Where do you plan to run GEOS-Chem?\n"
    printf "(e.g. local compute cluster, AWS, other supercomputer)?${thinline}"
    env_type=$(userInput)

    printf "${thinline}Please briefly describe how you plan on using GEOS-Chem\n"
    printf "so that we can add you to 'GEOS-Chem People and Projects'\n"
    printf "(https://geoschem.github.io/geos-chem-people-projects-map/)"
    printf "${thinline}"
    research_interest=$(userInput)

    # Send information to database on AWS
    postRegistration "${email}"             "${name}"       "${institution}"  \
                     "${name_of_pi}"        "${site}"       "${git_username}" \
                     "${research_interest}" "${model_type}" "${env_type}"

    # Update the .geoschem/config file and apply settings
    echo "export GC_USER_REGISTERED=true" >> "${HOME}/.geoschem/config"
    . ${HOME}/.geoschem/config
}
