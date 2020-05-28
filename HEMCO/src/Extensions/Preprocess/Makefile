#------------------------------------------------------------------------------
#                  Harvard-NASA Emissions Component (HEMCO)                   !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the HEMCO/Extensions/Preprocess subdirectory)
#
# !DESCRIPTION: This Makefile builds header files that replace ASCII file
#  reads with hard-coded F90 assignment statements.  This facilitates
#  file I/O in the ESMF environment.
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make -jN TARGET [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#
# !REVISION HISTORY:
#  11 Aug 2014 - R. Yantosca - Initial version
#  26 Feb 2015 - R. Yantosca - Remove lightdist script, we now read CDF table
#                              data from a file to reduce compilation time
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Shell script
SHELL     :=/bin/bash

# Directories
HCOX      :=..

# Files for GFED
GFED_IN  ?=$(HEMCO_DATA_ROOT)/GFED3/v2014-07/GFED3_emission_factors.txt
GFED_OUT ?=$(HCOX)/hcox_gfed_include.H

# Files for FINN
FINN_CO2  ?=$(HEMCO_DATA_ROOT)/FINN/v2014-07/FINN_EFratios_CO2.csv
FINN_VOC  ?=$(HEMCO_DATA_ROOT)/FINN/v2014-07/FINN_VOC_speciation.csv
FINN_OUT  ?=$(HCOX)/hcox_finn_include.H

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: all gfed finn lightdist paranox

all: gfed finn lightdist paranox

gfed:
	@echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	@echo "%%% GFED.......creating '$(GFED_OUT)'"
	@echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	./gfed.pl $(GFED_IN) $(GFED_OUT)

finn:
	@echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	@echo "%%% FINN........creating '$(FINN_OUT)'"
	@echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	./finn.pl $(FINN_CO2) $(FINN_VOC) $(FINN_OUT)
#EOC
