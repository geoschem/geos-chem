#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the Headers subdirectory)
#
# !DESCRIPTION: This makefile compiles the various GEOS-Chem Header modules,
#  which contain many PARAMETERs and global arrays for GEOS-Chem routines.
#
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% NOTE: Normally you will not have to call this Makefile directly,     %%%
# %%% it will be called automatically from the Makefile in the directory   %%%
# %%% just above this one!                                                 %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
# ROOTDIR    Specifies the root-level directory of the GEOS-Chem code
# HDR        Specifies the directory where GEOS-Chem include files are found
# LIB        Specifies the directory where library files (*.a) are stored
# MOD        Specifies the directory where module files (*.mod) are stored
# AR         Sys var w/ name of library creator program (i.e., "ar", "ranlib")
# MAKE       Sys var w/ name of Make command (i.e, "make" or "gmake")
#
# !REVISION HISTORY:
#  23 Aug 2011 - M. Long     - Initial version
#  See https://github.com/geoschem/geos-chem for complete history
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Directories
ROOT    :=..
LIB     :=$(ROOT)/lib
MOD     :=$(ROOT)/mod

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# List of source files
SOURCES :=$(wildcard *.F) $(wildcard *.F90)

# List of object files (replace .F and .F90 extensions with *.o)
TMP     :=$(SOURCES:.F=.o)
OBJECTS :=$(TMP:.F90=.o)

# List of module files.  Convert to lowercase, then prefix directory name.
MODULES :=$(OBJECTS:.o=.mod)
MODULES :=$(shell echo $(MODULES) | tr A-Z a-z)
MODULES :=$(foreach I,$(MODULES),$(MOD)/$(I))

# Library file
LIBRARY :=libHeaders.a

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean debug slowclean

lib: $(OBJECTS)
	$(AR) crs $(LIBRARY) $(OBJECTS)
	mv $(LIBRARY) $(LIB)

clean:
	@echo "===> Making clean in directory: Headers <==="
	@rm -f *.o *.mod *.a *.x

slowclean:
	@echo "===> Making slowclean in directory: Headers <==="
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

charpak_mod.o          : charpak_mod.F90

CMN_DIAG_mod.o         : CMN_DIAG_mod.F90         \
                         CMN_SIZE_mod.o           \
                         CMN_FJX_MOD.o            \
                         errcode_mod.o            \
                         precision_mod.o          \
                         state_grid_mod.o

CMN_FJX_MOD.o          : CMN_FJX_MOD.F90          \
                         CMN_SIZE_mod.o           \
                         precision_mod.o          \
                         state_grid_mod.o

CMN_O3_mod.o           : CMN_O3_mod.F90           \
                         CMN_SIZE_mod.o           \
                         errcode_mod.o            \
                         precision_mod.o          \
                         state_grid_mod.o

CMN_SIZE_mod.o         : CMN_SIZE_mod.F90         \
                         errcode_mod.o            \
                         input_opt_mod.o          \
                         precision_mod.o          \
                         state_grid_mod.o

commsoil_mod.o         : commsoil_mod.F90         \
                         CMN_SIZE_mod.o           \
                         errcode_mod.o

diaglist_mod.o         : diaglist_mod.F90         \
                         charpak_mod.o            \
                         inquireMod.o             \
                         registry_mod.o

dictionary_m.o         : dictionary_m.F90

errcode_mod.o          : errcode_mod.F90          \
                         charpak_mod.o

input_opt_mod.o        : input_opt_mod.F90        \
                         errcode_mod.o            \
                         precision_mod.o

inquiremod.o           : inquireMod.F90

state_chm_mod.o        : state_chm_mod.F90        \
                         CMN_SIZE_mod.o           \
                         CMN_FJX_MOD.o            \
                         charpak_mod.o            \
                         dictionary_m.o           \
                         physconstants.o          \
                         errcode_mod.o            \
                         input_opt_mod.o          \
                         precision_mod.o          \
                         registry_mod.o           \
                         species_mod.o            \
                         species_database_mod.o   \
                         state_grid_mod.o

state_diag_mod.o       : state_diag_mod.F90       \
                         CMN_SIZE_mod.o           \
                         CMN_FJX_MOD.o            \
                         charpak_mod.o            \
                         diaglist_mod.o           \
                         dictionary_m.o           \
                         errcode_mod.o            \
                         precision_mod.o          \
                         input_opt_mod.o          \
                         registry_mod.o           \
                         state_chm_mod.o          \
                         state_grid_mod.o

state_grid_mod.o       : state_grid_mod.F90       \
                         errcode_mod.o            \
                         precision_mod.o          \
                         input_opt_mod.o

state_met_mod.o        : state_met_mod.F90        \
                         charpak_mod.o            \
                         errcode_mod.o            \
                         dictionary_m.o           \
                         CMN_SIZE_mod.o           \
                         precision_mod.o          \
                         registry_mod.o           \
                         state_grid_mod.o

taggeddiaglist_mod.o   : taggeddiaglist_mod.F90   \
                         charpak_mod.o            \
                         diaglist_mod.o           \

physconstants.o        : physconstants.F90        \
                         precision_mod.o

registry_mod.o         : registry_mod.F90         \
                         charpak_mod.o            \
                         dictionary_m.o           \
                         errcode_mod.o            \
                         input_opt_mod.o          \
                         precision_mod.o          \
                         registry_params_mod.o

registry_params_mod.o  : registry_params_mod.F90 \
                         precision_mod.o

species_mod.o          : species_mod.F90          \
                         charpak_mod.o            \
                         input_opt_mod.o          \
                         physconstants.o          \
                         precision_mod.o

species_database_mod.o : species_database_mod.F90 \
                         errcode_mod.o            \
                         input_opt_mod.o          \
                         precision_mod.o          \
                         species_mod.o

#EOC
