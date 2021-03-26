#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile
#
# !DESCRIPTION: Makefile for GEOS-Chem run directories.
#\\
#\\
# !REVISION HISTORY:
#  Navigate to your unit tester directory and type 'gitk' at the prompt
#  to browse the revision history.
#EOP
#------------------------------------------------------------------------------
#BOC

# Unix shell (we'll assume Bash, which is on every Linux system)
SHELL        :=/bin/bash

###############################################################################
#####                                                                     #####
#####   CONFIGURABLE TOKENS: You can modify these for your environment.   #####
#####                                                                     #####
###############################################################################

# GEOS-Chem version number (X.Y.Z)
ifndef VERSION
 VERSION     :=12.9.3
endif

# If VERSION_TAG is not specified, then prefix "GC_" to the version number.
# This will prevent log files from beginning with a number.
ifndef VERSION_TAG
 VERSION_TAG :=GC_$(VERSION)
endif

# Source code location (you can modify as necessary)
ifndef CODE_DIR
 CODE_DIR    :=./CodeDir
endif

# Check if the code directory is not valid
NO_CODE_FOUND:=$(shell test -d $(CODE_DIR); echo $$?)
ifeq ($(NO_CODE_FOUND),1)
 $(error Could not find source code directory CODE_DIR=$(CODE_DIR))
endif

# Run directory path
ifndef RUN_DIR
 RUN_DIR     :=$(shell pwd)
endif

# getRunInfo perl script location (default is run directory)
ifndef PERL_DIR
 PERL_DIR    :=$(RUN_DIR)
endif

# Get start & end dates from the "getRunInfo" perl script
START        :=$(shell $(PERL_DIR)/getRunInfo $(RUN_DIR) 1)
END          :=$(shell $(PERL_DIR)/getRunInfo $(RUN_DIR) 2)
STARTDATE    :=$(shell $(PERL_DIR)/getRunInfo $(RUN_DIR) 3)
ENDDATE      :=$(shell $(PERL_DIR)/getRunInfo $(RUN_DIR) 4)
SIM          :=$(shell $(PERL_DIR)/getRunInfo $(RUN_DIR) 5)

###############################################################################
#####                                                                     #####
#####  DEFAULT COMPILATION OPTIONS: Set various switches if not defined   #####
#####                                                                     #####
###############################################################################

ifndef OMP
 OMP         :=y
endif

ifndef TRACEBACK
 TRACEBACK   :=y
endif

ifndef BOUNDS
 BOUNDS      :=n
endif

ifndef DEBUG
 DEBUG       :=n
endif

ifndef FPE
 FPE         :=n
endif

ifndef FPEX
 FPEX        :=n
endif

ifndef TAU_PROF
 TAU_PROF    :=n
endif

ifndef APM
 APM         :=n
endif

ifndef TOMAS12
 TOMAS12     :=n
endif

ifndef TOMAS15
 TOMAS15     :=n
endif

ifndef TOMAS30
 TOMAS30     :=n
endif

ifndef TOMAS40
 TOMAS40     :=n
endif

ifndef RRTMG
 RRTMG       :=n
endif

# Now turn OFF BPCH_DIAG by default
ifndef BPCH_DIAG
 BPCH_DIAG   :=n
endif

###############################################################################
#####                                                                     #####
#####  RUN OPTIONS: Get various settings from the run directory name,     #####
#####               or from the type of simulation that is being used.    #####
#####                                                                     #####
###############################################################################

# Define the proper setting for CHEM, which picks the correct KPP files
ifeq ($(SIM),TOMAS12)
 TOMAS12     :=y
 CHEM        :=Tropchem
 BPCH_DIAG   :=y
endif

ifeq ($(SIM),TOMAS15)
 TOMAS15     :=y
 CHEM        :=Tropchem
 BPCH_DIAG   :=y
endif

ifeq ($(SIM),TOMAS30)
 TOMAS30     :=y
 CHEM        :=Tropchem
 BPCH_DIAG   :=y
endif

ifeq ($(SIM),TOMAS40)
 TOMAS40     :=y
 CHEM        :=Tropchem
 BPCH_DIAG   :=y
endif

ifeq ($(SIM),RRTMG)
 RRTMG       :=y
 CHEM        :=Tropchem
 BPCH_DIAG   :=y
endif

ifeq ($(SIM),standard)
 CHEM        :=Standard
endif

ifeq ($(SIM),benchmark)
 CHEM        :=Standard
endif

ifeq ($(SIM),aciduptake)
 CHEM        :=Standard
endif

ifeq ($(SIM),marinePOA)
 CHEM        :=Standard
endif

ifeq ($(SIM),complexSOA)
 CHEM        :=Tropchem
endif

ifeq ($(SIM),complexSOA_SVPOA)
 CHEM        :=SOA_SVPOA
 NEST        :=n
endif

ifeq ($(SIM),tropchem)
 CHEM        :=Tropchem
endif

ifeq ($(SIM),custom)
 CHEM        :=Custom
endif

ifeq ($(SIM),APM)
 APM         :=y
 CHEM        :=Tropchem
endif

ifndef CHEM
 CHEM        :=Standard
endif

ifeq ($(SIM),Hg)
 BPCH_DIAG   :=y
endif

ifeq ($(SIM),tagHg)
 BPCH_DIAG   :=y
endif

ifeq ($(SIM),POPs)
 BPCH_DIAG   :=y
endif

# General run information
TIMESTAMP    :=$(shell date +"%Y/%m/%d %H:%M")

# Log files that will be written to the log directory
LOG_COMP     :=compile.log
LOG_GC       :=$(VERSION_TAG).log
LOG_HMCO     :=HEMCO.log
LOG_HMCO_SA  :=HEMCO.standalone.log
LOG_HMCO_DR  :=HEMCO.dryrun.log
BUILDINFO    :=lastbuild

# Executables
EXE          :=geos
EXE_HMCO_SA  :=hemco_standalone.x

# Get information about the Git version, because some features
# (like -C) are not available in older Git versions.  The git
# version command returns "git version X.Y.Z", so we will just take
# the 3rd word and strip all the dots. (bmy, 12/21/16)
GIT_VERSION    :=$(subst .,,$(word 3, $(shell git --version)))
NEWER_THAN_185 :=$(shell perl -e "print ($(GIT_VERSION) gt 185)")

ifeq ($(NEWER_THAN_185),1)

 # Git version 1.8.5 and higher uses the -C option to look in
 # directories other than the current directory.  Use this to
 # get info about the last committed state of the code.
 CODE_BRANCH  :=$(shell git -C $(CODE_DIR) rev-parse --abbrev-ref HEAD)
 LAST_COMMIT  :=$(shell git -C $(CODE_DIR) log -n 1 --pretty=format:"%s")
 COMMIT_DATE  :=$(shell git -C $(CODE_DIR) log -n 1 --pretty=format:"%cd")
 COMMIT_USER  :=$(shell git -C $(CODE_DIR) log -n 1 --pretty=format:"%cn")
 COMMIT_HASH  :=$(shell git -C $(CODE_DIR) log -n 1 --pretty=format:"%h")
 CODE_STATUS  :=$(shell git -C $(CODE_DIR) status -s)

else

 # Git versions older than 1.8.5 don't have the -C option,
 # so we have to manually change to the code dir to get information
 # about the last committed state of the code. (bmy, 12/21/16)
 CODE_BRANCH  :=$(shell cd $(CODE_DIR); git rev-parse --abbrev-ref HEAD; cd $(PWD))
 LAST_COMMIT  :=$(shell cd $(CODE_DIR); git log -n 1 --pretty=format:"%s"; cd $(PWD))
 COMMIT_DATE  :=$(shell cd $(CODE_DIR); git log -n 1 --pretty=format:"%cd"; cd $(PWD))
 COMMIT_USER  :=$(shell cd $(CODE_DIR); git log -n 1 --pretty=format:"%cn"; cd $(PWD))
 COMMIT_HASH  :=$(shell cd $(CODE_DIR); git log -n 1 --pretty=format:"%h"; cd $(PWD))
 CODE_STATUS  :=$(shell cd $(CODE_DIR); git status -s; cd $(PWD))

endif

# Error check, make sure to strip out quote characters from commit string
# or else GNU Make will think that the string is ending prematurely.
# See: https://stackoverflow.com/questions/10424645/how-to-convert-a-quoted-string-to-a-normal-one-in-makefile
LAST_COMMIT := $(subst $\",,$(LAST_COMMIT))
LAST_COMMIT := $(subst $\',,$(LAST_COMMIT))

# Get compiler version
COMPILER_VERSION_LONG :=$(shell $(FC) --version))
COMPILER_VERSION_LONG :=$(sort $(COMPILER_VERSION_LONG))

# For ifort, the 3rd substring of the sorted text is the version number.
# For pgfortran and gfortran, it's the 4th substring.
# NOTE: Future compiler updates may break this algorithm.
REGEXP      :=(^[Ii][Ff][Oo][Rr][Tt])
ifeq ($(shell [[ "$(FC)" =~ $(REGEXP) ]] && echo true),true)
 COMPILER_VERSION :=$(word 3, $(COMPILER_VERSION_LONG))
else
 COMPILER_VERSION :=$(word 4, $(COMPILER_VERSION_LONG))
endif

# Variable to hold the TAU_OPTIONS environment setting
TAU_OPT     :=""

# If we are compiling for use w/ TAU, then also set TAU_OPTIONS accordingly
# TAU_SF sets the -optTauSelectFile option for removing throttled files
REGEXP      :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(TAU_PROF)" =~ $(REGEXP) ]] && echo true),true)
  ifdef TAU_SF
    TAU_OPT :="$(TAU_OPTIONS) -optTauSelectFile=$(TAU_SF)"
  else
    TAU_OPT :="$(TAU_OPTIONS)"
  endif
endif

###############################################################################
#####                                                                     #####
#####                          MAKEFILE TARGETS                           #####
#####                                                                     #####
###############################################################################

# PHONY targets don't actually compile anything. They either are
# synonyms for other targets, they remove files, or they print output.
.PHONY: all             dataclean        logclean       execlean
.PHONY: cleanup_output  cleanup_data     cleanup_logs   cleanup_exe
.PHONY: printruninfo    printbuildinfo   printcodeinfo  printallinfo
.PHONY: hemco           hcobuild         hcorun         cleanup_hco
.PHONY: hcodataclean    printruninfohco  build

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Default make is to print help screen            %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default:
	@$(MAKE) help

#%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Compile Source Code    %
#%%%%%%%%%%%%%%%%%%%%%%%%%%

# Build and run
all: clean build run

#  Build only
build:
	@$(MAKE) cleanup_logs
	@$(MAKE) cleanup_exe
	@$(MAKE) -C $(CODE_DIR) APM=$(APM)                                    \
                                CHEM=$(CHEM)                                  \
                                RRTMG=$(RRTMG)                                \
                                TOMAS12=$(TOMAS12)                            \
                                TOMAS15=$(TOMAS15)                            \
                                TOMAS30=$(TOMAS30)                            \
                                TOMAS40=$(TOMAS40)                            \
                                OMP=$(OMP)                                    \
                                TAU_PROF=$(TAU_PROF)                          \
                                TAU_OPTIONS=$(TAU_OPT)                        \
                                BPCH_DIAG=$(BPCH_DIAG)                        \
                                2>&1 | tee -a $(LOG_COMP)
	cp -f $(CODE_DIR)/bin/geos $(EXE)
	@$(MAKE) printbuildinfo 2>&1 | tee -a $(LOG_COMP)
	@$(MAKE) printbuildinfo >> $(BUILDINFO)

# Run only
run:
	@$(MAKE) cleanup_output
	./geos > $(LOG_GC)
	@$(MAKE) printruninfo >> $(LOG_GC)
	@$(MAKE) printruninfo

# Run only
dryrun:
	@$(MAKE) cleanup_output
	./geos --dryrun > $(LOG_GC)
	@$(MAKE) printruninfo >> $(LOG_GC)
	@$(MAKE) printruninfo

#%%%%%%%%%%%%%%%%%%%%%%%%%%
#  HEMCO Standalone       %
#%%%%%%%%%%%%%%%%%%%%%%%%%%

# Synonyms
hemco:      hco
hemcobuild: hcobuild
hemcorun:   hcorun

# Build and run
hco:
	@$(MAKE) hcobuild
	@$(MAKE) hcorun

# Build only
hcobuild:
	@$(MAKE) -C $(CODE_DIR) APM=$(APM)                                    \
                                CHEM=$(CHEM)                                  \
                                RRTMG=$(RRTMG)                                \
                                TOMAS12=$(TOMAS12)                            \
                                TOMAS15=$(TOMAS15)                            \
                                TOMAS30=$(TOMAS30)                            \
                                TOMAS40=$(TOMAS40)                            \
                                OMP=yes                                       \
                                TAU_PROF=$(TAU_PROF)                          \
                                TAU_OPTIONS=$(TAU_OPT)                        \
                                BPCH_DIAG=$(BPCH_DIAG)                        \
                                libhemcosa 2>&1 | tee -a $(LOG_COMP)
	cp -f $(CODE_DIR)/bin/hemco_standalone.x $(EXE_HMCO_SA) 2>&1 | tee -a $(LOG_COMP)
	@$(MAKE) printbuildinfo 2>&1 | tee -a $(LOG_COMP)
	@$(MAKE) printbuildinfo > $(BUILDINFO)

# Run only
hcorun:
	@$(MAKE) cleanup_hco
	./$(EXE_HMCO_SA) -c HEMCO_sa_Config.rc  > $(LOG_HMCO_SA)
	@$(MAKE) printruninfohco >> $(LOG_HMCO_SA)
	@$(MAKE) printruninfohco

# Dry Run only
hcodryrun:
	@$(MAKE) cleanup_hco
	./$(EXE_HMCO_SA) -c HEMCO_sa_Config.rc --dryrun > $(LOG_HMCO_DR)
	@$(MAKE) printruninfohco >> $(LOG_HMCO_DR)
	@$(MAKE) printruninfohco

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Clean up run directory                          %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Synonyms
dataclean: cleanup_data
hcodataclean: cleanup_hco
logclean: cleanup_logs
execlean: cleanup_exe

cleanup_output: cleanup_data cleanup_logs

cleanup_data:
	rm -f HEMCO_diagnostics*
	if [ -f HEMCO_restart.$(START).nc ]; then mv -f HEMCO_restart.$(START).nc HC.nc; fi;
	rm -f HEMCO_restart.*;
	if [ -f HC.nc ]; then mv -f HC.nc HEMCO_restart.$(START).nc; fi;
	rm -f trac_avg*
	rm -f Ox.mass.*
	rm -f diaginfo.dat tracerinfo.dat
	rm -f tracer_mass*.dat
	rm -f GEOS-Chem_Species_Database.json
	rm -f core.*
	rm -f OutputDir/*
	@for f in ./*.nc4; do                               \
           if [[ "$$f" != "./GEOSChem.Restart.$(STARTDATE)_0000z.nc4" ]]; \
	      then rm -f $$f; fi;              \
        done

cleanup_hco:
	rm -f HEMCO_restart*
	rm -f diaginfo.dat tracerinfo.dat
	rm -f output/*.nc

cleanup_logs:
	rm -f *.log
	rm -f slurm-*

cleanup_exe:
	rm -f $(EXE)
	rm -f $(EXE_HMCO_SA)
	rm -f $(LOG_COMP)
	rm -f $(BUILDINFO)

#%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Clean Source           %
#%%%%%%%%%%%%%%%%%%%%%%%%%%

clean:
	@$(MAKE) -C $(CODE_DIR) clean

realclean:
	@$(MAKE) -C $(CODE_DIR) realclean

tauclean:
	@$(MAKE) -C $(CODE_DIR) tauclean


#%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Clean and Remove All   %
#%%%%%%%%%%%%%%%%%%%%%%%%%%

superclean:
	@$(MAKE) realclean
	@$(MAKE) tauclean
	@$(MAKE) cleanup_output
	@$(MAKE) cleanup_exe

#%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Print information      %
#%%%%%%%%%%%%%%%%%%%%%%%%%%

printruninfo:
	@echo "RUN SETTINGS:"
	@echo "  Run directory       : $(RUN_DIR)"
	@echo "  Run ID              : $(RUNID)"
	@echo "  Simulation Start    : $(START)"
	@echo "  Simulation End      : $(END)"
	@echo "  Version number      : $(VERSION)"
	@echo "  Version tag         : $(VERSION_TAG)"
	@echo "  GC Diagnostic File  : $(TRAC_AVG)"
	@echo "  GC Restart File     : $(GC_RST)"
	@echo "  HEMCO Restart File  : $(HMCO_RST)"
	@echo "  GC Run Log File     : $(LOG_SP)"
	@echo "  HEMCO Log File      : $(LOG_HMCO)"

printruninfohco:
	@echo "RUN SETTINGS:"
	@echo "  Run directory       : $(RUN_DIR)"
	@echo "  Run ID              : $(RUNID)"
	@echo "  Simulation Start    : $(START)"
	@echo "  Simulation End      : $(END)"
	@echo "  HEMCO End           : $(HMCO_END)"
	@echo "  Version             : $(VERSION)"
	@echo "  HEMCO Restart File  : $(HMCO_RST)"
	@echo "  HEMCO Run Log File  : $(LOG_HMCO_SA)"
	@echo "  HEMCO Log File      : $(LOG_HMCO)"

printbuildinfo:
	@echo "LAST BUILD INFORMATION:"
	@echo "  CODE_DIR     : $(CODE_DIR)"
	@echo "  CODE_BRANCH  : $(CODE_BRANCH)"
	@echo "  LAST_COMMIT  : $(LAST_COMMIT)"
	@echo "  COMMIT_DATE  : $(COMMIT_DATE)"
	@echo "  VERSION      : $(VERSION)"
	@echo "  VERSION_TAG  : $(VERSION_TAG)"
	@echo "  TRACEBACK    : $(TRACEBACK)"
	@echo "  BOUNDS       : $(BOUNDS)"
	@echo "  DEBUG        : $(DEBUG)"
	@echo "  FPE          : $(FPE)"
	@echo "  OMP          : $(OMP)"
	@echo "  CHEM         : $(CHEM)"
	@echo "  TOMAS12      : $(TOMAS12)"
	@echo "  TOMAS15      : $(TOMAS15)"
	@echo "  TOMAS30      : $(TOMAS30)"
	@echo "  TOMAS40      : $(TOMAS40)"
	@echo "  APM          : $(APM)"
	@echo "  RRTMG        : $(RRTMG)"
	@echo "  TAU_PROF     : $(TAU_PROF)"
	@echo "  BPCH_DIAG    : $(BPCH_DIAG)"
	@echo "  COMPILER     : $(FC) $(COMPILER_VERSION)"
	@echo "  Datetime     : $(TIMESTAMP)"

printallinfo:
	@$(MAKE) printbuildinfo
	@$(MAKE) printruninfo

printcodeinfo:
	@echo -e "Code directory:  $(CODE_DIR)"
	@echo -e "Git branch:      $(CODE_BRANCH)"
	@echo -e "Last commit:"
	@echo -e "   Message:      $(LAST_COMMIT)"
	@echo -e "   Date:         $(COMMIT_DATE)"
	@echo -e "   Committer:    $(COMMIT_USER)"
	@echo -e "   Hash abbrev:  $(COMMIT_HASH)"
	@echo -e "Uncommitted files (if any):\n$(CODE_STATUS)"

###############################################################################
#####                                                                     #####
#####                             HELP SCREEN                             #####
#####                                                                     #####
###############################################################################

help:
	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	@echo '%%%    GEOS-Chem Run Directory Makefile Options        %%%'
	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	@echo
	@echo 'Usage: make -jN TARGET [OPTIONAL-FLAGS]'
	@echo ''
	@echo '-jN               Compiles N files at a time (reduces compile time)'
	@echo ''
	@echo '----------------------------------------------------------'
	@echo 'TARGET may be one of the following:'
	@echo '----------------------------------------------------------'
	@echo ''
	@echo '%%% COMPILE AND RUN %%%'
	@echo 'all               Cleans, compiles, and runs GEOS-Chem'
	@echo 'hco               Cleans, compiles, and runs HEMCO'
	@echo ''
	@echo '%%% BUILD ONLY %%%'
	@echo 'build             Compiles GEOS-Chem'
	@echo 'hcobuild          Compiles HEMCO'
	@echo ''
	@echo '%%% RUN ONLY %%%'
	@echo 'run               Runs GEOS-Chem'
	@echo 'dryrun            Runs GEOS-Chem in dry-run mode'
	@echo 'hcorun            Runs HEMCO'
	@echo 'hcodryrun         Runs HEMCO in dry-run mode'
	@echo ''
	@echo '%%% CLEAN UP %%%'
	@echo 'cleanup_output    Synonym for: cleanup_data cleanup_logs'
	@echo 'cleanup_data      Removes all GEOS-Chem diagnostic and restart file'
	@echo 'cleanup_hco       Removes all HEMCO diagnostic and restart files'
	@echo 'cleanup_logs      Removes all GEOS-Chem and HEMCO output log files'
	@echo 'cleanup_exe       Removes all GEOS-Chem and HEMCO executable files'
	@echo ''
	@echo '%%% CLEAN SOURCE CODE %%%'
	@echo 'clean             Makes "clean" in source code directory $CODE_DIR'
	@echo 'realclean         Makes "realclean" in the source code directory $CODE_DIR'
	@echo 'tauclean          Removes all TAU files in the source code directory $CODE_DIR'
	@echo ''
	@echo '%%% CLEAN AND REMOVE ALL %%%'
	@echo 'superclean        Synonym for: realclean tauclean cleanup_output cleanup_exe'
	@echo ''
	@echo '%%% PRINT INFORMATION %%%'
	@echo 'printruninfo      Prints the run settings for GEOS-Chem simulations'
	@echo 'printruninfohco   Prints the run settings for HEMCO standalone simulations'
	@echo 'printbuildinfo    Prints the build settings for GEOS-Chem simulations'
	@echo 'printallinfo      Synonym for: printbuildinfo printruninfo'
	@echo 'printcodeinfo     Print code directory git information'
	@echo ''
	@echo '%%% UTILITIES %%%'
	@echo 'help              Prints this help screen'
	@echo ''
	@echo '----------------------------------------------------------'
	@echo 'OPTIONAL-FLAGS may be one of the following:'
	@echo '----------------------------------------------------------'
	@echo 'DEBUG=y           Compiles GEOS-Chem with various debugging options'
	@echo 'BOUNDS=y          Compiles GEOS-Chem with out-of-bounds error checks'
	@echo 'FPE=y             Compiles GEOS-Chem with floating-point math error checks'
	@echo 'TRACEBACK=y       Compiles GEOS-Chem with traceback error printout (this is the default)'
#EOC
