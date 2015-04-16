#------------------------------------------------------------------------------
#                  Harvard-NASA Emissions Component (HEMCO)                   !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the HEMCO/Core subdirectory)
#
# !DESCRIPTION: This Makefile builds the HEMCO core routines into library
#  file libHCOI.a.
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#
#   make -jN TARGET [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#
# !REVISION HISTORY: 
#  21 Jul 2014 - R. Yantosca - Initial version, based on HEMCO standalone
#                              and modified specifically for GEOS-Chem
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Directories
ROOT :=../..
BIN  :=$(ROOT)/bin
LIB  :=$(ROOT)/lib
MOD  :=$(ROOT)/mod

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# List of source code files to compile
# (i.e. everything ending in either .F or .F90)
SRC  := $(wildcard *.F) $(wildcard *.F90)

# Replace .F and .F90 extensions with *.o
TMP  := $(SRC:.F=.o)
OBJ  := $(TMP:.F90=.o)

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean debug

all: lib

lib: $(OBJ)
	$(AR) crs libHCO.a $(OBJ)
	mv libHCO.a $(LIB)

clean:
	rm -f *.o *.mod *.a *.x

debug:
	@echo "Targets : $(MAKECMDGOALS)"
	@echo "ROOT    : $(ROOT)"
	@echo "BIN     : $(BIN)"
	@echo "LIB     : $(LIB)"
	@echo "MOD     : $(MOD)"
	@echo "F90     : $(F90)"

###############################################################################
###                                                                         ###
###  Dependencies listing                                                   ###
###  (grep "USE " to get the list of module references!)                    ###
###                                                                         ###
###  From this list of dependencies, the "make" utility will figure out     ###
###  correct order of compilation (so we don't have to do that ourselves).  ###
###  This also allows us to compile on multiple processors with "make -j".  ###
###                                                                         ###
###  NOTES:                                                                 ###
###  (1) Only specify object-file dependencies that are within this         ###
###       directory.  Object files in other directories will be referenced  ### 
###       at link-time.                                                     ###
###  (2) For "make -jN" (i.e. compile N files simultaneously), all files    ###
###       in this directory must have a listed dependency.                  ###
###                                                                         ###
###############################################################################

hco_arr_mod.o          : hco_arr_mod.F90         \
                         hco_error_mod.o

hco_calc_mod.o         : hco_calc_mod.F90        \
                         hco_arr_mod.o           \
                         hco_datacont_mod.o      \
                         hco_diagn_mod.o         \
                         hco_extlist_mod.o       \
                         hco_emislist_mod.o      \
                         hco_error_mod.o         \
                         hco_filedata_mod.o      \
                         hco_state_mod.o         \
                         hco_tidx_mod.o

hco_chartools_mod.o    : hco_chartools_mod.F90   \
	                 hco_error_mod.o         \
                         hco_extlist_mod.o

hco_clock_mod.o        : hco_clock_mod.F90       \
                         hco_emislist_mod.o      \
                         hco_error_mod.o         \
                         hco_state_mod.o

hco_config_mod.o       : hco_config_mod.F90      \
                         hco_chartools_mod.o     \
                         hco_datacont_mod.o      \
                         hco_diagn_mod.o         \
                         hco_error_mod.o         \
                         hco_extlist_mod.o       \
                         hco_filedata_mod.o      \
                         hco_readlist_mod.o      \
                         hco_state_mod.o         \
                         hco_tidx_mod.o

hco_datacont_mod.o     : hco_datacont_mod.F90    \
                         hco_arr_mod.o           \
                         hco_error_mod.o         \
                         hco_filedata_mod.o      \
                         hco_state_mod.o

hco_diagn_mod.o        : hco_diagn_mod.F90       \
                         hco_arr_mod.o           \
                         hco_clock_mod.o         \
                         hco_error_mod.o         \
                         hco_extlist_mod.o       \
                         hco_state_mod.o         \
                         hco_unit_mod.o

hco_driver_mod.o       : hco_driver_mod.F90      \
                         hco_calc_mod.o          \
                         hco_clock_mod.o         \
                         hco_config_mod.o        \
                         hco_datacont_mod.o      \
                         hco_diagn_mod.o         \
                         hco_emislist_mod.o      \
                         hco_error_mod.o         \
                         hco_extlist_mod.o       \
                         hco_readlist_mod.o      \
                         hco_state_mod.o         \
                         hco_tidx_mod.o          \
                         hcoio_diagn_mod.o

hco_emislist_mod.o     : hco_emislist_mod.F90    \
                         hco_datacont_mod.o      \
                         hco_error_mod.o         \
                         hco_filedata_mod.o      \
                         hco_logfile_mod.o       \
                         hco_state_mod.o

hco_extlist_mod.o      : hco_extlist_mod.F90     \
                         hco_error_mod.o

hco_error_mod.o        : hco_error_mod.F90

hco_filedata_mod.o     : hco_filedata_mod.F90    \
                         hco_arr_mod.o           \
                         hco_chartools_mod.o     \
                         hco_error_mod.o         \
                         hco_unit_mod.o

hco_fluxarr_mod.o      : hco_fluxarr_mod.F90     \
                         hco_arr_mod.o           \
                         hco_error_mod.o         \
                         hco_state_mod.o

hco_geotools_mod.o     : hco_geotools_mod.F90    \
	                 hco_error_mod.o

hco_interp_mod.o       : hco_interp_mod.F90      \
                         hco_datacont_mod.o      \
                         hco_error_mod.o         \
                         hco_filedata_mod.o      \
                         hco_state_mod.o

hco_logfile_mod.o      : hco_logfile_mod.F90     \
                         hco_datacont_mod.o      \
                         hco_error_mod.o         \
                         hco_state_mod.o

hco_readlist_mod.o     : hco_readlist_mod.F90    \
                         hco_clock_mod.o         \
                         hco_datacont_mod.o      \
                         hco_emislist_mod.o      \
                         hco_error_mod.o         \
                         hco_filedata_mod.o      \
                         hco_logfile_mod.o       \
                         hco_state_mod.o         \
                         hco_tidx_mod.o          \
                         hcoio_dataread_mod.o

hco_restart_mod.o      : hco_restart_mod.F90     \
                         hco_diagn_mod.o         \
                         hco_emislist_mod.o      \
                         hco_error_mod.o         \
                         hco_state_mod.o

hco_state_mod.o        : hco_state_mod.F90       \
	                 hco_arr_mod.o           \
                         hco_chartools_mod.o     \
                         hco_error_mod.o         \
                         hco_extlist_mod.o       \
                         hco_unit_mod.o

hco_tidx_mod.o         : hco_tidx_mod.F90        \
                         hco_chartools_mod.o     \
                         hco_clock_mod.o         \
                         hco_datacont_mod.o      \
                         hco_error_mod.o         \
                         hco_filedata_mod.o      \
                         hco_state_mod.o

hco_unit_mod.o         : hco_unit_mod.F90        \
                         hco_chartools_mod.o     \
                         hco_error_mod.o         \
                         hco_extlist_mod.o

hcoio_diagn_mod.o      : hcoio_diagn_mod.F90     \
                         hco_clock_mod.o         \
                         hco_diagn_mod.o         \
                         hco_error_mod.o         \
                         hco_state_mod.o

hcoio_dataread_mod.o   : hcoio_dataread_mod.F90  \
                         hco_chartools_mod.o     \
                         hco_datacont_mod.o      \
                         hco_emislist_mod.o      \
                         hco_error_mod.o         \
                         hco_filedata_mod.o      \
                         hco_geotools_mod.o      \
                         hco_interp_mod.o        \
                         hco_state_mod.o         \
                         hco_tidx_mod.o          \
                         hco_unit_mod.o          \
                         hcoio_messy_mod.o

hcoio_messy_mod.o      : hcoio_messy_mod.F90     \
                         hco_datacont_mod.o      \
                         hco_error_mod.o         \
                         hco_filedata_mod.o      \
                         hco_interp_mod.o        \
                         hco_state_mod.o         \
                         hco_unit_mod.o          \
                         messy_ncregrid_base.o
                         
messy_ncregrid_base.o  : messy_ncregrid_base.F90 \
                         hco_error_mod.o         \
                         messy_ncregrid_mpi.o

messy_ncregrid_mpi.o   : messy_ncregrid_mpi.F90

#EOC
