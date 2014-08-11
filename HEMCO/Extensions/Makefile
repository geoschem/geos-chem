#------------------------------------------------------------------------------
#                  Harvard-NASA Emissions Component (HEMCO)                   !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the HEMCO/Extensions subdirectory)
#
# !DESCRIPTION: This Makefile builds the HEMCO extensions routines into 
#  library file libHCOX.a.
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
#  15 Jul 2014 - R. Yantosca - Initial version
#  22 Jul 2014 - R. Yantosca - Now compile lightning_cdf_mod.F90
#  22 Jul 2014 - R. Yantosca - Turn off optimization for lightning_cdf_mod.F90
#  08 Aug 2014 - R. Yantosca - Now reference include files w/ hardwired
#                              emissions data. This facilitates I/O for ESMF.
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Define variables
ROOT :=../..
BIN  :=$(ROOT)/bin
LIB  :=$(ROOT)/lib
MOD  :=$(ROOT)/mod

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# List of source files to compile (everything ending in .F and .F90)
SRC  :=$(wildcard *.F) $(wildcard *.F90)

# For now, skip commsoil_mod.F90 till we figure out what to do with it
SRC  :=$(filter-out commsoil_mod.F90,$(SRC))

# Replace .f and .f90 extensions with *.o
TMP  :=$(SRC:.F=.o)
OBJ  :=$(TMP:.F90=.o)

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean debug

all: lib

lib: $(OBJ)
	$(AR) crs libHCOX.a $(OBJ)
	mv libHCOX.a $(LIB)

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

drydep_toolbox_mod.o  : drydep_toolbox_mod.F90

EF_MGN20_mod.o        : EF_MGN20_mod.F

hcox_custom_mod.o     : hcox_custom_mod.F90       \
                        hcox_state_mod.o

hcox_driver_mod.o     : hcox_driver_mod.F90       \
                        hcox_custom_mod.o         \
                        hcox_dustdead_mod.o       \
                        hcox_dustginoux_mod.o     \
                        hcox_finn_mod.o           \
                        hcox_gc_RnPbBe_mod.o      \
                        hcox_gfed3_mod.o          \
                        hcox_paranox_mod.o        \
                        hcox_lightnox_mod.o       \
                        hcox_megan_mod.o          \
                        hcox_paranox_mod.o        \
                        hcox_seaflux_mod.o        \
                        hcox_seasalt_mod.o        \
                        hcox_soilnox_mod.o        \
                        hcox_state_mod.o

hcox_dustdead_mod.o   : hcox_dustdead_mod.F       \
                        hcox_state_mod.o

hcox_dustginoux_mod.o : hcox_dustginoux_mod.F90   \
                        hcox_state_mod.o

hcox_gc_RnPbBe_mod.o  : hcox_gc_RnPbBe_mod.F90    \
                        hcox_gc_RnPbBe_include.H  \
                        hcox_state_mod.o          

hcox_gfed3_mod.o      : hcox_gfed3_mod.F90        \
                        hcox_gfed3_include.H      \
                        hcox_state_mod.o

hcox_finn_mod.o       : hcox_finn_mod.F90         \
                        hcox_finn_include.H       \
                        hcox_state_mod.o

hcox_lightnox_mod.o   : hcox_lightnox_mod.F90     \
                        hcox_state_mod.o          \
                        lightning_cdf_mod.o

hcox_megan_mod.o      : hcox_megan_mod.F          \
                        hcox_state_mod.o          \
                        EF_MGN20_mod.o

hcox_paranox_mod.o    : hcox_paranox_mod.F90      \
                        hcox_state_mod.o          \
                        paranox_util_mod.o

hcox_seaflux_mod.o    : hcox_seaflux_mod.F90      \
                        hcox_state_mod.o          \
                        henry_mod.o               \
                        ocean_toolbox_mod.o

hcox_seasalt_mod.o    : hcox_seasalt_mod.F90      \
                        hcox_state_mod.o

hcox_soilnox_mod.o    : hcox_soilnox_mod.F90      \
                        hcox_state_mod.o          \
                        drydep_toolbox_mod.o

hcox_state_mod.o      : hcox_state_mod.F90


henry_mod.o           : henry_mod.F

lightning_cdf_mod.o   : lightning_cdf_mod.F90     \
                        lightning_cdf_include.H

ocean_toolbox_mod.o   : ocean_toolbox_mod.F90

paranox_util_mod.o    : paranox_util_mod.F

#EOC
