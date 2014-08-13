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

# Files for GFED3
GFED3_IN  ?=$(HEMCO_DATA_ROOT)/GFED3/v2014-07/GFED3_emission_factors.txt
GFED3_OUT ?=$(HCOX)/hcox_gfed3_include.H

# Files for FINN
FINN_CO2  ?=$(HEMCO_DATA_ROOT)/FINN/v2014-07/FINN_EFratios_CO2.csv
FINN_VOC  ?=$(HEMCO_DATA_ROOT)/FINN/v2014-07/FINN_VOC_speciation.csv
FINN_OUT  ?=$(HCOX)/hcox_finn_include.H

# Files for Lightning CDF
LDIST_IN  ?=$(HEMCO_DATA_ROOT)/LIGHTNOX/v2014-07/light_dist.ott2010.dat
LDIST_OUT ?=$(HCOX)/lightning_cdf_include.H

# Files for PARANOX
PNOX_FN   ?=$(HEMCO_DATA_ROOT)/PARANOX/v2014-07/FracNOx_binary_5hrs_20gs.dat
PNOX_OPE  ?=$(HEMCO_DATA_ROOT)/PARANOX/v2014-07/IntOPE_binary_5hrs_20gs.dat
PNOX_OUT  ?=$(HCOX)/hcox_gfed3_include.H

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: all gfed finn lightdist paranox

all: gfed3 finn lightdist paranox

gfed3:
	@echo "%%% GFED3.......creating '$(GFED3_OUT)'"
	@gfed3.pl $(GFED3_IN) $(GFED3_OUT)

finn:
	@echo "%%% FINN........creating '$(FINN_OUT)'"
	@finn.pl $(FINN_CO2) $(FINN_VOC) $(FINN_OUT)

lightdist:
	@echo "%%% LIGHTDIST...creating '$(LDIST_OUT)'"
	@lightdist.pl $(LDIST_IN) $(LDIST_OUT)

paranox:
	@echo "%%% PARANOX.....creating '$(PNOX_OUT)'"
#	./paranox.pl $(PNOX_FN) $(PNOX_OPE) $(PNOX_OUT)

#EOC
