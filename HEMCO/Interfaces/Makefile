#------------------------------------------------------------------------------
#                  Harvard-NASA Emissions Component (HEMCO)                   !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the HEMCO/Interfaces subdirectory)
#
# !DESCRIPTION: This Makefile builds the HEMCO interface routines into 
#  library file libHCOI.a.
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
DOC  :=$(ROOT)/doc
LIB  :=$(ROOT)/lib
MOD  :=$(ROOT)/mod

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# List of source files to compile
SRC  := $(wildcard *.F) $(wildcard *.F90)

# Replace .f and .f90 extensions with *.o
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
	$(AR) crs libHCOI.a $(OBJ)
	mv libHCOI.a $(LIB)
clean:
	rm -f *.o *.mod *.a *.x

debug:
	@echo "Targets : $(MAKECMDGOALS)"
	@echo "ROOT    : $(ROOT)"
	@echo "BIN     : $(BIN)"
	@echo "LIB     : $(LIB)"
	@echo "MOD     : $(MOD)"
	@echo "F90     : $(F90)"
	@echo "LINK    : $(LINK)"

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

hcoi_esmf_mod.o : hcoi_esmf_mod.F90

hcoi_standalone_mod.o : hcoi_standalone_mod.F90

#EOC
