#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile
#
# !DESCRIPTION: Utility Makefile for the GCHP run directory.
#\\
#\\
# !REMARKS:
# 
# !REVISION HISTORY: 
#  See git history in the GCHP repository for version 12.1.0 and beyond.
#  For prior history, see the git history in geos-chem-unittest
#EOP
#------------------------------------------------------------------------------
#BOC

# Unix shell (we'll assume Bash, which is on every Linux system)
SHELL :=/bin/bash

###############################################################################
#####                                                                     #####
#####                              VARIABLES                              #####
#####                                                                     #####
###############################################################################

# Source code location
ifndef CODEDIR_GC
 CODEDIR_GC :=$(shell readlink -f ./CodeDir)
endif

# GCHP code directory path
CODEDIR_GCHP :=$(CODEDIR_GC)/GCHP

# Run directory path
RUN_DIR :=$(shell pwd)

# Log files that will be written to the log directory
BUILD_LOG  :="$(RUN_DIR)/compile.log"

# Executables
EXE :=geos

###############################################################################
#####                                                                     #####
#####                              TARGETS                                #####
#####                                                                     #####
###############################################################################

# PHONY targets don't actually compile anything. They either are
# synonyms for other targets, they remove files, or they print output.
.PHONY: cleanup_data     clean_all      build_all   build_all_debug
.PHONY: cleanup_logs     clean_core     build_core  build_core_debug
.PHONY: cleanup_output   clean_mapl     build_mapl  build_mapl_debug
.PHONY: cleanup_exe      clean_esmf     build_esmf  build_esmf_debug
.PHONY: superclean       help           printbuildinfo

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Default make is to print help screen            %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default:
	@$(MAKE) help

#%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Clean Source Code      %
#%%%%%%%%%%%%%%%%%%%%%%%%%%

# Clean all source code
clean_all:
	./build.sh clean_all

# Clean ESMF, MAPL, FVdycore, and GEOS-Chem core (same as clean_all)
clean_esmf:
	@$(MAKE) clean_all

# Clean only MAPL, FVdycore, and GEOS-Chem core
clean_mapl:
	./build.sh clean_mapl

# Clean only GEOS-Chem, including top-level GCHP directory files
clean_core:
	./build.sh clean_core

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Compile All Source Code  %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Build without cleaning
build_all:
	@$(MAKE) build_esmf

# Build with GEOS-Chem core debug flags and without cleaning
build_all_debug:
	@$(MAKE) build_esmf_debug

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Compile ESMF, MAPL, FVdycore, and GEOS-Chem  %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Build without cleaning
build_esmf:
	rm -f geos
	rm -f $(CODEDIR_GC)/bin/geos
	date > ${BUILD_LOG}
	rm -f $(CODEDIR_GCHP)/ESMF/esmf.install
	rm -f $(CODEDIR_GCHP)/Shared/mapl.install
	rm -f $(CODEDIR_GCHP)/FVdycoreCubed_GridComp/fvdycore.install
	./build.sh build 2>&1 | tee -a $(BUILD_LOG)
	@$(MAKE) printbuildinfo 2>&1 | tee -a $(BUILD_LOG)


# Build with GEOS-Chem core debug flags and without cleaning
build_esmf_debug:
	rm -f geos
	rm -f $(CODEDIR_GC)/bin/geos
	date > ${BUILD_LOG}
	rm -f $(CODEDIR_GCHP)/ESMF/esmf.install
	rm -f $(CODEDIR_GCHP)/Shared/mapl.install
	rm -f $(CODEDIR_GCHP)/FVdycoreCubed_GridComp/fvdycore.install
	./build.sh build --debug 2>&1 | tee -a $(BUILD_LOG)
	@$(MAKE) printbuildinfo 2>&1 | tee -a $(BUILD_LOG)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Compile MAPL, FVdycore, and GEOS-Chem but not ESMF %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Build without cleaning
build_mapl:
	rm -f geos
	rm -f $(CODEDIR_GC)/bin/geos
	date > ${BUILD_LOG}
	rm -f $(CODEDIR_GCHP)/Shared/mapl.install
	rm -f $(CODEDIR_GCHP)/FVdycoreCubed_GridComp/fvdycore.install
	./build.sh build 2>&1 | tee -a $(BUILD_LOG)
	@$(MAKE) printbuildinfo 2>&1 | tee -a $(BUILD_LOG)

# Build with GEOS-Chem core debug flags and without cleaning
build_mapl_debug:
	rm -f geos
	rm -f $(CODEDIR_GC)/bin/geos
	date > ${BUILD_LOG}
	rm -f $(CODEDIR_GCHP)/Shared/mapl.install
	rm -f $(CODEDIR_GCHP)/FVdycoreCubed_GridComp/fvdycore.install
	./build.sh build --debug 2>&1 | tee -a $(BUILD_LOG)
	@$(MAKE) printbuildinfo 2>&1 | tee -a $(BUILD_LOG)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Compile GEOS-Chem but not ESMF, MAPL, and FVdycore %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Build without cleaning
build_core:
	rm -f geos
	rm -f $(CODEDIR_GC)/bin/geos
	date > ${BUILD_LOG}
	./build.sh build 2>&1 | tee -a $(BUILD_LOG)
	@$(MAKE) printbuildinfo 2>&1 | tee -a $(BUILD_LOG)

# Build with GEOS-Chem core debug flags and without cleaning
build_core_debug:
	rm -f geos
	rm -f $(CODEDIR_GC)/bin/geos
	date > ${BUILD_LOG}
	./build.sh build --debug 2>&1 | tee -a $(BUILD_LOG)
	@$(MAKE) printbuildinfo 2>&1 | tee -a $(BUILD_LOG)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Clean up Run Directory  %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%

cleanup_data: 
	rm -f $(RUN_DIR)/OutputDir/*.nc4
	rm -f trac_avg.*
	rm -f tracerinfo.dat
	rm -f diaginfo.dat
	rm -f cap_restart
	rm -f gcchem*
	rm -f *.rcx
	rm -f *~

cleanup_logs: 
	rm -f gchp.log
	rm -f HEMCO.log
	rm -f PET*.log
	rm -f multirun.log
	rm -f logfile.000000.out
	rm -f slurm-*
	rm -f 1
	rm -f EGRESS

cleanup_output: cleanup_data cleanup_logs

cleanup_exe: 
	rm -f geos
	rm -f $(BUILD_LOG)

# Clean all source code and the run directory
superclean:
	@$(MAKE) clean_all
	@$(MAKE) cleanup_exe
	@$(MAKE) cleanup_output

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Print information                      %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printbuildinfo:
	$(eval CODE_BRANCH_GC :=$(shell git -C $(CODEDIR_GC) rev-parse --abbrev-ref HEAD))
	$(eval LAST_COMMIT_GC :=$(shell git -C $(CODEDIR_GC) log -n 1 --pretty=format:"%s")) 
	$(eval COMMIT_DATE_GC :=$(shell git -C $(CODEDIR_GC) log -n 1 --pretty=format:"%cd")) 
	$(eval COMMIT_USER_GC :=$(shell git -C $(CODEDIR_GC) log -n 1 --pretty=format:"%cn")) 
	$(eval COMMIT_HASH_GC :=$(shell git -C $(CODEDIR_GC) log -n 1 --pretty=format:"%h")) 
	$(eval CODE_STATUS_GC :=$(shell git -C $(CODEDIR_GC) status -s)) 
	$(eval CODE_BRANCH_GCHP:=$(shell git -C $(CODEDIR_GCHP) rev-parse --abbrev-ref HEAD))
	$(eval LAST_COMMIT_GCHP:=$(shell git -C $(CODEDIR_GCHP) log -n 1 --pretty=format:"%s"))
	$(eval COMMIT_DATE_GCHP:=$(shell git -C $(CODEDIR_GCHP) log -n 1 --pretty=format:"%cd"))
	$(eval COMMIT_USER_GCHP:=$(shell git -C $(CODEDIR_GCHP) log -n 1 --pretty=format:"%cn"))
	$(eval COMMIT_HASH_GCHP:=$(shell git -C $(CODEDIR_GCHP) log -n 1 --pretty=format:"%h")) 
	$(eval CODE_STATUS_GCHP:=$(shell git -C $(CODEDIR_GCHP) status -s)) 
	@echo "GEOS-Chem repository"
	@echo "   Path        : $(CODEDIR_GC)"
	@echo "   Branch      : $(CODE_BRANCH_GC)"
	@echo "   Last commit : $(LAST_COMMIT_GC)"
	@echo "   Date        : $(COMMIT_DATE_GC)"
	@echo "   User        : $(COMMIT_USER_GC)"
	@echo "   Hash        : $(COMMIT_HASH_GC)"
	@echo "   Git status  : $(CODE_STATUS_GC)"
	@echo "GCHP repository"
	@echo "   Path        : $(CODEDIR_GCHP)"
	@echo "   Branch      : $(CODE_BRANCH_GCHP)"
	@echo "   Last commit : $(LAST_COMMIT_GCHP)"
	@echo "   Date        : $(COMMIT_DATE_GCHP)"
	@echo "   User        : $(COMMIT_USER_GCHP)"
	@echo "   Hash        : $(COMMIT_HASH_GCHP)"
	@echo "   Git status  : $(CODE_STATUS_GCHP)"

help:
	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	@echo '%%%    GCHP Run Directory Makefile Options    %%%'
	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	@echo ' '
	@echo 'Usage: make TARGET [OPTIONAL-FLAGS]'
	@echo ' '
	@echo '----------------------------------------------------------'
	@echo 'TARGET may be one of the following:'
	@echo '----------------------------------------------------------'
	@echo ' '      
	@echo 'Print make options:' 
	@echo '  printbuildinfo     Print out git info for GC and GCHP repos'
	@echo '  help               Shows the display you are looking at now'
	@echo ' '
	@echo 'Clean source code:'
	@echo '  clean_all          Clean ESMF, MAPL, FVdycore, and GEOS-Chem (all)'
	@echo '  clean_esmf         Same as clean_esmf'
	@echo '  clean_mapl         Clean MAPL, FVdycore, and GEOS-Chem (not ESMF)'
	@echo '  clean_core         Clean GEOS-Chem only'
	@echo '  superclean         Clean all source code and the run directory'
	@echo ' '
	@echo 'Compile all (ESMF, MAPL, FVdycore, and GEOS-Chem):'
	@echo '  build_all          Compile without cleaning; GC core debug off'
	@echo '  build_esmf         Same as build_all'
	@echo '  build_esmf_debug   Compile without cleaning; GC core debug on'
	@echo ' '
	@echo 'Compile MAPL, FVdycore, and GEOS-Chem (not ESMF):'
	@echo '  build_mapl         Compile without cleaning; GC core debug off'
	@echo '  build_mapl_debug   Compile without cleaning; GC core debug on'
	@echo ' '
	@echo 'Compile only GEOS-Chem (not ESMF, MAPL, or FVdycore):'
	@echo '  build_core         Compile without cleaning; GC core debug off'
	@echo '  build_core_debug   Compile without cleaning; GC core debug on'
	@echo ' '
	@echo 'To remove run directory files:' 
	@echo '  cleanup_output     Remove output data, log files, and executable'
	@echo '  cleanup_data       Remove output data' 
	@echo '  cleanup_logs       Remove log files' 
	@echo '  cleanup_exe        Remove executable'
	@echo '  superclean         All of the above and clean all source code'
	@echo ' ' 
	@echo 'Print source code info:' 
	@echo '  printbuildinfo     Print out git info for GC and GCHP repos'
	@echo ' '
	@echo '----------------------------------------------------------'
	@echo 'Tips for when to use build targets:'
	@echo '----------------------------------------------------------'
	@echo ' '      
	@echo ' 1. For a fresh install:'
	@echo '      make clean_all; make build_all : cleans and compiles all source code'
	@echo ' 2. Subsequent compile depends on what you changed:'      
	@echo '  a. If you change GEOS-Chem files or files in GCHP source code directory (common):'      
	@echo '      make build_core                : skips recompiling MAPL, ESMF, and FV3 advection'      
	@echo '      make build_core_debug          : same as above but with debug compile flags on'      
	@echo '  b. If you change files in the Shared, FV3, or Registry directories (advanced):'      
	@echo '      make build_mapl                : skips recompiling ESMF'      
	@echo '      make build_mapl_debug          : same as above but with debug compile flags on'      
	@echo '  c. If you change files in the ESMF directory (should not happen; make github issue):'      
	@echo '      make build_esmf                : recompiles all source code'      
	@echo '      make build_esmf_debug          : same as above but with debug compile flags on'      
	@echo '  d. If you change libraries, e.g. switch compiler:'      
	@echo '      make clean_all; make build_all : cleans and recompiles all source code'      
	@echo ' '      
#EOC


