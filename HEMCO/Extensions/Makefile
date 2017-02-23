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
#  20 Aug 2014 - M. Sulprizio- Add hcox_gc_POPs_mod.F90
#  01 Oct 2014 - R. Yantosca - Add hcox_tomas_seasalt_mod.F90
#  03 Jun 2015 - R. Yantosca - Also remove *.mod, *.a files with "make clean"
#  04 Jun 2015 - R. Yantosca - Remove EF_MGN20.F, this isn't used anymore
#  07 Dec 2015 - R. Yantosca - Restore fast "clean" command; add "slowclean"
#  01 Nov 2016 - M. Sulprizio- Rename hcox_tomas_seasalt_mod.F90 to
#                              hcox_tomas_jeagle_mod.F90 (J. Kodros)
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Directories
ROOT    :=../..
LIB     :=$(ROOT)/lib
MOD     :=$(ROOT)/mod

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# List of source files to compile (everything ending in .F and .F90)
SOURCES :=$(wildcard *.F) $(wildcard *.F90)

# For now, skip commsoil_mod.F90 till we figure out what to do with it
SOURCES :=$(filter-out commsoil_mod.F90,$(SOURCES))

# List of object files (replace .f and .f90 extensions with *.o)
TMP     :=$(SOURCES:.F=.o)
OBJECTS :=$(TMP:.F90=.o)

# List of module files.  Convert to lowercase, then prefix directory name.
MODULES :=$(OBJECTS:.o=.mod)
MODULES :=$(shell echo $(MODULES) | tr A-Z a-z)
MODULES :=$(foreach I,$(MODULES),$(MOD)/$(I))

# Library file
LIBRARY :=libHCOX.a

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean debug slowclean

all: lib

lib: $(OBJECTS)
	$(AR) crs $(LIBRARY) $(OBJECTS)
	mv $(LIBRARY) $(LIB)

clean:
	@echo "===> Making clean in directory: HEMCO/Extensions <==="
	@rm -f *.o *.mod *.a *.x

slowclean:
	@echo "===> Making slowclean in directory: HEMCO/Extensions <==="
	@rm -f $(OBJECTS) $(MODULES) $(LIBRARY) $(LIB)/$(LIBRARY)

debug:
	@echo "Targets : $(MAKECMDGOALS)"
	@echo "ROOT    : $(ROOT)"
	@echo "LIB     : $(LIB)"
	@echo "MOD     : $(MOD)"
	@echo "F90     : $(F90)"
	@echo "OBJECTS : $(OBJECTS)"
	@echo "MODULES : $(MODULES)"
	@echo "LIBRARY : $(LIBRARY)"

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

drydep_toolbox_mod.o        : drydep_toolbox_mod.F90

hcox_aerocom_mod.o          : hcox_aerocom_mod.F90            \
                              hcox_state_mod.o		      \
                              hcox_tools_mod.o

hcox_ch4wetland_mod.o       : hcox_ch4wetland_mod.F90         \
                              hcox_state_mod.o		      \
                              hcox_tools_mod.o

hcox_custom_mod.o           : hcox_custom_mod.F90             \
                              hcox_state_mod.o		      

hcox_driver_mod.o           : hcox_driver_mod.F90             \
                              hcox_aerocom_mod.o              \
                              hcox_ch4wetland_mod.o           \
                              hcox_custom_mod.o               \
                              hcox_dustdead_mod.o             \
                              hcox_tomas_dustdead_mod.o       \
                              hcox_dustginoux_mod.o           \
                              hcox_finn_mod.o                 \
                              hcox_gc_POPs_mod.o              \
                              hcox_gc_RnPbBe_mod.o            \
                              hcox_gfed_mod.o                 \
                              hcox_iodine_mod.o               \
                              hcox_paranox_mod.o              \
                              hcox_lightnox_mod.o             \
                              hcox_megan_mod.o                \
                              hcox_paranox_mod.o              \
                              hcox_seaflux_mod.o              \
                              hcox_seasalt_mod.o              \
                              hcox_soilnox_mod.o              \
                              hcox_state_mod.o		      \
                              hcox_tomas_jeagle_mod.o        


hcox_dustdead_mod.o         : hcox_dustdead_mod.F             \
                              hcox_state_mod.o		      

hcox_dustginoux_mod.o       : hcox_dustginoux_mod.F90         \
                              hcox_state_mod.o		      

hcox_gc_POPs_mod.o          : hcox_gc_POPs_mod.F90            \
                              hcox_state_mod.o		      

hcox_gc_RnPbBe_mod.o        : hcox_gc_RnPbBe_mod.F90          \
                              hcox_state_mod.o                

hcox_gfed_mod.o             : hcox_gfed_mod.F90               \
                              hcox_gfed_include_gfed3.H       \
                              hcox_gfed_include_gfed4.H       \
                              hcox_state_mod.o		      \
                              hcox_tools_mod.o

hcox_iodine_mod.o           : hcox_iodine_mod.F90             \
                              hcox_state_mod.o		      

hcox_finn_mod.o             : hcox_finn_mod.F90               \
                              hcox_finn_include.H             \
                              hcox_state_mod.o		      \
                              hcox_tools_mod.o

hcox_lightnox_mod.o         : hcox_lightnox_mod.F90           \
                              hcox_state_mod.o                \
                              hcox_tools_mod.o

hcox_megan_mod.o            : hcox_megan_mod.F                \
                              hcox_state_mod.o

hcox_paranox_mod.o          : hcox_paranox_mod.F90            \
                              hcox_state_mod.o

hcox_seaflux_mod.o          : hcox_seaflux_mod.F90            \
                              hcox_state_mod.o                \
                              ocean_toolbox_mod.o	      

hcox_seasalt_mod.o          : hcox_seasalt_mod.F90            \
                              hcox_state_mod.o		      

hcox_soilnox_mod.o          : hcox_soilnox_mod.F90            \
                              hcox_state_mod.o                \
                              hcox_tools_mod.o                \
                              drydep_toolbox_mod.o

hcox_state_mod.o            : hcox_state_mod.F90

hcox_tomas_jeagle_mod.o     : hcox_tomas_jeagle_mod.F90      \
                              hcox_state_mod.o

hcox_tools_mod.o            : hcox_tools_mod.F90

hcox_tomas_dustdead_mod.o   : hcox_tomas_dustdead_mod.F       \
                              hcox_state_mod.o

ocean_toolbox_mod.o         : ocean_toolbox_mod.F90

#EOC
