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
#  03 Jun 2015 - R. Yantosca - Also remove *.mod, *.a files with "make clean"
#  10 Jun 2015 - R. Yantosca - Also place a copy of the executable in bin
#  23 Jun 2015 - M. Sulprizio- Remove executable from bin with "make clean"
#  07 Dec 2015 - R. Yantosca - Restore fast "clean" command; add "slowclean"
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
BIN     :=$(ROOT)/bin

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# List of source files to compile
SOURCES :=$(wildcard *.F) $(wildcard *.F90)

# Replace .f and .f90 extensions with *.o
TMP     :=$(SOURCES:.F=.o)
OBJECTS :=$(TMP:.F90=.o)

# List of module files.  Convert to lowercase, then prefix directory name.
MODULES :=$(OBJECTS:.o=.mod)
MODULES :=$(shell echo $(MODULES) | tr A-Z a-z)
MODULES :=$(foreach I,$(MODULES),$(MOD)/$(I))

# Library file
LIBRARY :=libHCOI.a

# Executable file
EXE     :=hemco_standalone.x

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean debug check slowclean

all:
	@${MAKE} lib
	@${MAKE} exe

check:
	@${MAKE} exe

lib: $(OBJECTS)
	$(AR) crs $(LIBRARY) $(OBJECTS)
	mv $(LIBRARY) $(LIB)

clean:
	@echo "===> Making clean in directory: HEMCO/Interfaces <==="
	@rm -f *.o *.mod *.a *.x $(EXE) $(BIN)/$(EXE)

slowclean:
	@echo "===> Making slowclean in directory: HEMCO/Interfaces <==="
	@rm -f *.x *$(OBJECTS) $(MODULES) $(LIBRARY) $(LIB)/$(LIBRARY)
	@rm -f $(BIN)/$(EXE)

exe:
	$(LD) $(OBJECTS) $(LINK_HCO) -o $(EXE)
	cp -f $(EXE) $(BIN)

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

hcoi_esmf_mod.o : hcoi_esmf_mod.F90

hcoi_standalone_mod.o : hcoi_standalone_mod.F90

hemco_standalone.o : hemco_standalone.F90 \
                     hcoi_standalone_mod.o

#EOC
