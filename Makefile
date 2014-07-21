#------------------------------------------------------------------------------
#                  Harvard-NASA Emissions Component (HEMCO)                   !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the HEMCO/src directory)
#
# !DESCRIPTION: Calls makefiles in the subdirectories src\/Core, 
#  src\/Extensions, src\/Interfaces to compile the HEMCO source code into
#  library files and to create an executable.
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
#  14 Jul 2014 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Directories
ROOT :=..
BIN  :=$(ROOT)/bin
DOC  :=$(ROOT)/doc
HCO  :=$(ROOT)/Core
HCOI :=$(ROOT)/Interfaces
HCOX :=$(ROOT)/Extensions
LIB  :=$(ROOT)/lib
MOD  :=$(ROOT)/mod
RUN  :=$(ROOT)/run

# Include header file, which returns various Makefile variables, defines 
# the compilation rules, and sets the correct C-preprocessor switches.
include $(ROOT)/Makefile_header.mk

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete list!                ###
###                                                                         ###
###############################################################################

.PHONY: all check clean distclean debug test

#-----------------------------------------
# Targets for building code
#-----------------------------------------

all:
	@$(MAKE) lib
	@$(MAKE) exe

check:
	@$(MAKE) -C $(RUN) all

exe: check

lib:
	@$(MAKE) libHCO
#	@$(MAKE) libHCOX
#	@$(MAKE) libHCOI

libHCO:
	@$(MAKE) -C $(HCO) lib

libHCOI:
	@$(MAKE) -C $(HCOI) lib

libHCOX:
	@$(MAKE) -C $(HCOX) lib

test: check

#-----------------------------------------
# Targets for cleaning up code
#-----------------------------------------

clean:
	@$(MAKE) -C $(HCO)  clean
	@$(MAKE) -C $(HCOI) clean
	@$(MAKE) -C $(HCOX) clean
	@$(MAKE) -C $(RUN)  clean

distclean: realclean

realclean:
	@$(MAKE) clean
	@$(MAKE) docclean
	rm -f $(LIB)/*.a
	rm -f $(MOD)/*.mod
	rm -f $(BIN)/*.x
	rm -f $(RUN)/*.x

#-----------------------------------------
# Targets for building the documentation
#-----------------------------------------

doc:
	@$(MAKE) -C $(DOC) doc

docclean: 
	@$(MAKE) -C $(DOC) clean

#-----------------------------------------
# Targets for debugging and help screen
#-----------------------------------------

debug:
	@echo "%%% Makefile variable settings %%%"
	@echo "SHELL      : $(SHELL)"
	@echo "%%% Directories %%%"
	@echo "ROOT       : $(ROOT)"
	@echo "HCO        : $(HCO)"
	@echo "HCOI       : $(HCOI)"
	@echo "HCOX       : $(HCOX)"
	@echo "HELP       : $(HELP)"
	@echo "DOC        : $(DOC)"
	@echo "LIB        : $(LIB)"
	@echo "MOD        : $(MOD)"
	@echo "RUN        : $(RUN)"
	@echo "%%% For the NcdfUtilities %%%"
	@echo "NCU_BIN    : $(NCU_BIN)"
	@echo "NCU_MOD    : $(NCU_MOD)"
	@echo "NCU_LIB    : $(NCU_LIB)"
	@echo "%%% For netCDF library paths %%%"
	@echo "BIN_NETCDF : $(BIN_NETCDF)"
	@echo "INC_NETCDF : $(INC_NETCDF)"
	@echo "LIB_NETCDF : $(LIB_NETCDF)"
	@echo "%%%% Compilation commands %%%"
	@echo "CC         : $(CC)"
	@echo "R8         : $(R8)"
	@echo "FREEFORM   : $(FREEFORM)"
	@echo "F90        : $(F90)"
	@echo "%%% Linking commands %%%"
	@echo "NC_LINK    : $(NC_LINK)"
	@echo "LINK       : $(LINK)"
	@echo "LD         : $(LD)"

help:
	@$(MAKE) -C $(HELP) help
#EOC



