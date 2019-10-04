#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile
#
# !DESCRIPTION: 
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% NOTE: Normally you will not have to call this Makefile directly,     %%%
# %%% it will be called automatically from the Makefile in a directory     %%%
# %%% above this one!                                                      %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# !REVISION HISTORY: 
#   19 Sep 2019 - E. Lundgren - Initial version
#   See git history
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

# BASEDIR is a synonym for ROOTDIR
ifndef BASEDIR
 export BASEDIR=$(realpath $(ROOT))
endif 

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
LIBRARY :=libGCHPint.a

# Auto-generated header files
REGDIR    := Registry
ACGS      := GIGCchem_ExportSpec___.h GIGCchem_GetPointer___.h \
             GIGCchem_DeclarePointer___.h GIGCchem_History___.rc

# RCDIR is the directory where the registry files live
export RCDIR=$(REGDIR)

# Include MAPL files
include ../../GCHP/Shared/Config/ESMA_base.mk  # Generic stuff
include ../../GCHP/Shared/Config/ESMA_arch.mk  # System dependencies

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean debug slowclean

lib: $(ACGS) $(OBJECTS)
	$(AR) crs $(LIBRARY) $(OBJECTS)
	mv $(LIBRARY) $(LIB)

$(ACGS) : $(REGDIR)/Chem_Registry.rc $(REGDIR)/HEMCO_Registry.rc $(ACG)
	@echo "test"
	@echo $(ACG)
	@$(ACG) $(ACG_FLAGS) $(REGDIR)/Chem_Registry.rc
	@$(ACG) $(ACG_FLAGS) $(REGDIR)/HEMCO_Registry.rc

clean:
	@echo "===> Making clean in directory: Interfaces/GCHP <==="
	rm -f *.o *.mod *___.h *___.rc

slowclean:
	@echo "===> Making slowclean in directory: Interfaces/GCHP <==="
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

Chem_GridCompMod.o          : Chem_GridCompMod.F90                          \
			      gigc_chunk_mod.o                              \
                              gigc_historyexports_mod.o                     \
                              gigc_providerservices_mod.o

gigc_chunk_mod.o            : gigc_chunk_mod.F90                            \
                              gigc_historyexports_mod.o                     \

gigc_historyexports_mod.o   : gigc_historyexports_mod.F90

gigc_providerservices_mod.o : gigc_providerservices_mod.F90

#EOC

