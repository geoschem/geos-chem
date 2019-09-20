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

#
#
#
#
#
#
################################################################################
####                                                                         ###
####  Initialization section                                                 ###
####                                                                         ###
################################################################################
#
## Use the bash shell
#SHELL=/bin/bash
#
## ROOTDIR is the directory just above this one.
#ROOTDIR=..
#
## Directories where GEOS-Chem code live
#HDR=$(ROOTDIR)/Headers
#HELP=$(ROOTDIR)/help
#LIB=$(ROOTDIR)/lib
#MOD=$(ROOTDIR)/mod
#
## This code here should get the exact number of the intel version (MDY)
#ifeq ($(ESMF_COMPILER),intel)
#  INTELVERSIONTEXT :=$(shell ifort --version))
#  INTELVERSIONTEXT :=$(sort $(INTELVERSIONTEXT))
#  LONGVERSION      :=$(word 3, $(INTELVERSIONTEXT))
#  MAJORVERSION     :=$(subst ., ,$(LONGVERSION))
#  MAJORVERSION     :=$(firstword $(MAJORVERSION))
#  MAJORVERSION     :=$(strip $(MAJORVERSION))
#else ifeq ($(ESMF_COMPILER),gfortran)
#  GNUVERSIONTEXT   :=$(shell gfortran -dumpversion))
#  LONGVERSION      :=$(GNUVERSIONTEXT))
#  MAJORVERSION     :=$(subst ., ,$(LONGVERSION))
#  MAJORVERSION     :=$(firstword $(MAJORVERSION))
#  MAJORVERSION     :=$(strip $(MAJORVERSION))
#else
#  MAJORVERSION     :="0"
#endif
#
## Include header file.  This returns variables CC, F90, FREEFORM, LD, R8,
## as well as the default Makefile compilation rules for source code files.
#include $(ROOTDIR)/Makefile_header.mk
#
## BASEDIR is a synonym for ROOTDIR
#ifndef BASEDIR
# export BASEDIR=$(realpath $(ROOTDIR))
#endif 
#
##-----------------------------------------------------------------------------
## Settings for ESMF
##-----------------------------------------------------------------------------
#
## Directory where ESMF lives
#ifndef ESMF_DIR
#  export ESMF_DIR=$(CURDIR)/ESMF
#endif
#
## Compiler for ESMF - eg intel, intelgcc, gfortran
#ifndef ESMF_COMPILER
#   $(error ESMF_COMPILER is not defined)
#endif
#
## MPI type for ESMF
#ifndef ESMF_COMM
#   $(error ESMF_COMM is not defined)
#endif
#
## Operating system type for ESMF
#ifndef ESMF_OS
#  export ESMF_OS=$(ARCH)
#endif
#
## Optimization level for ESMF
#ifndef ESMF_BOPT
#  export ESMF_BOPT=O
#endif
#
## Other ESMF directory settings
#export ESMF_INSTALL_PREFIX=$(ESMF_DIR)/$(ARCH)
#export ESMF_INSTALL_LIBDIR=$(ESMF_DIR)/$(ARCH)/lib
#export ESMF_INSTALL_MODDIR=$(ESMF_DIR)/$(ARCH)/mod
#export ESMF_INSTALL_HEADERDIR=$(ESMF_DIR)/$(ARCH)/include
#
## Other ESMF compilation settings
#ifeq ($(ESMF_COMPILER),intel)
#  export ESMF_F90COMPILEOPTS=-align all -fPIC -traceback 
#else ifeq ($(ESMF_COMPILER),gfortran)
#  export ESMF_F90COMPILEOPTS=-falign-commons -fPIC -fbacktrace
#else
#  export ESMF_F90COMPILEOPTS=-fPIC
#endif
#export ESMF_CXXCOMPILEOPTS=-fPIC
#export ESMF_OPENMP=OFF
#
##-----------------------------------------------------------------------------
## Settings for MAPL
##-----------------------------------------------------------------------------
#
## ESMADIR is the directory where MAPL lives
#ifndef ESMADIR
#  ESMADIR:=$(CURDIR)/Shared
#endif
#
## Fortran compiler for MAPL
#export ESMA_FC=$(FC)
#
## RCDIR is the directory where the registry files live
#export RCDIR=$(CURDIR)/Registry
#
## Inline the proper include files
#ifeq ($(HPC),yes)
#  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
#  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
#else
#  include ./Shared/Config/ESMA_base.mk  # Generic stuff
#  include ./Shared/Config/ESMA_arch.mk  # System dependencies
#endif
#
##-----------------------------------------------------------------------------
## Settings for FVdycore
##-----------------------------------------------------------------------------
#
## FVDIR is where the FVDycore lives
#ifndef FVDIR
# export FVDIR=$(CURDIR)/FVdycoreCubed_GridComp
#endif
#
################################################################################
####                                                                         ###
####  List of files to compile                                               ###
####                                                                         ###
################################################################################
#
## List of source files
#SRC = $(wildcard *.F) $(wildcard *.F90)
#
## Replace .F and .F90 extensions with *.o
#TMP = $(SRC:.F=.o)
#OBJ = $(TMP:.F90=.o)
#
#REGDIR    := Registry
#ACGS      := GIGCchem_ExportSpec___.h GIGCchem_GetPointer___.h \
#             GIGCchem_DeclarePointer___.h GIGCchem_History___.rc
##LIB_ESMF  := $(ESMF_DIR)/$(ARCH)/lib/libesmf.so
##LIB_MAPL  := $(ESMADIR)/$(ARCH)/libMAPL_Base.a # At this point, we only check for MAPL_Base
#
################################################################################
####                                                                         ###
####  Makefile targets: type "make help" for a complete listing!             ###
####                                                                         ###
################################################################################
#
#.PHONY: clean help baselibs
#
#baselibs:
#	@$(MAKE) baselibs_esmf
#	@$(MAKE) baselibs_mapl
#	@$(MAKE) baselibs_fvdycore
#
#baselibs_esmf:
#ifeq ($(wildcard $(ESMF_DIR)/esmf.install),)
#	$(MAKE) -C $(ESMF_DIR)
#	$(MAKE) -C $(ESMF_DIR) install
#	cp -f  $(ESMF_DIR)/src/include/*.inc $(ESMF_INSTALL_HEADERDIR)
#	@touch $(ESMF_DIR)/esmf.install
#endif
#
#
#baselibs_mapl:
#ifeq ($(wildcard $(ESMADIR)/mapl.install),)
#	$(MAKE) -C $(ESMADIR) install
#	rm -f $(ESMADIR)/$(ARCH)/lib/libGMAO_gfio.a
#	ln -sf $(ESMADIR)/$(ARCH)/lib/libGMAO_gfio_r4.a $(ESMADIR)/$(ARCH)/lib/libGMAO_gfio.a
#	rm -f $(ESMADIR)/$(ARCH)/lib/libMAPL_cfio.a
#	ln -sf $(ESMADIR)/$(ARCH)/lib/libMAPL_cfio_r4.a $(ESMADIR)/$(ARCH)/lib/libMAPL_cfio.a
#	@touch $(ESMADIR)/mapl.install
#endif
#
#baselibs_fvdycore:
#ifeq ($(wildcard $(FVDIR)/fvdycore.install),)
#	@echo "<><>"
#	$(MAKE) -C $(FVDIR) ESMADIR=$(ESMADIR) install
#	@touch $(FVDIR)/fvdycore.install
#endif
#
#lib: $(ACGS) $(OBJ)
#	$(AR) crs libGIGC.a $(OBJ)
#	mv libGIGC.a $(LIB)
#
#$(ACGS) : $(REGDIR)/Chem_Registry.rc $(REGDIR)/HEMCO_Registry.rc $(ACG) #$(REGDIR)/Dyn_Registry.rc
#	@$(ACG) $(ACG_FLAGS) $(REGDIR)/Chem_Registry.rc
###	@$(ACG) $(ACG_FLAGS) $(REGDIR)/Dyn_Registry.rc
#	@$(ACG) $(ACG_FLAGS) $(REGDIR)/HEMCO_Registry.rc
#
#libesmf:
#	@$(MAKE) -C $(GCHP) esmf
#
#libmapl:
#	@$(MAKE) -C $(GCHP) mapl
#
#clean:
#	rm -f *.o *.mod *___.h *___.rc
#
#help:
#	@$(MAKE) -C $(HELP)
#
#gigc_debug gigc_help:
#	@echo "Directories:"
#	@echo "----------------------------------------------------------"
#	@echo "Current working dir    : $(CURDIR)"
#	@echo "ROOTDIR                : $(ROOTDIR)"
#	@echo "HDR                    : $(HDR)"
#	@echo "HELP                   : $(HELP)"
#	@echo "LIB                    : $(LIB)"
#	@echo "MOD                    : $(MOD)"
#	@echo "BASEDIR                : $(BASEDIR)"
#	@echo ""
#	@echo "Settings for ESMF"
#	@echo "----------------------------------------------------------"
#	@echo "ESMF_DIR               : $(ESMF_DIR)"
#	@echo "ESMF_COMPILER          : $(ESMF_COMPILER)"
#	@echo "ESMF_COMM              : $(ESMF_COMM)"
#	@echo "ESMF_OS                : $(ESMF_OS)"
#	@echo "ESMF_BOPT              : $(ESMF_BOPT)"
#	@echo "ESMF_INSTALL_PREFIX    : $(ESMF_INSTALL_PREFIX)"
#	@echo "ESMF_INSTALL_LIBDIR    : $(ESMF_INSTALL_LIBDIR)"
#	@echo "ESMF_INSTALL_MODDIR    : $(ESMF_INSTALL_MODDIR)"
#	@echo "ESMF_INSTALL_HEADERDIR : $(ESMF_INSTALL_HEADERDIR)"
#	@echo "ESMF_F90COMPILEOPTS    : $(ESMF_F90COMPILEOPTS)"
#	@echo "ESMF_CXXCOMPILEOPTS    : $(ESMF_CXXCOMPILEOPTS)"
#	@echo ""
#	@echo "Settings for MAPL"
#	@echo "----------------------------------------------------------"
#	@echo "ESMADIR                : $(ESMADIR)"
#	@echo "ESMA_FC                : $(ESMA_FC)"
#	@echo "RCDIR                  : $(RCDIR)"
#	@echo ""
#	@echo "Settings for FVdycore"
#	@echo "----------------------------------------------------------"
#	@echo "FVDIR                  : $(FVDIR)"
#	@$(MAKE) -C $(FVDIR) ESMADIR=$(ESMADIR) help
#
################################################################################
####                                                                         ###
####  Targets to remove ESMF, MAPL, and FVDYCORE!                            ###
####  USE WITH EXTREME CAUTION!!!                                            ###
####                                                                         ###
################################################################################
#
#.PHONY: the_nuclear_option
#
#the_nuclear_option:
#	@$(MAKE) wipeout_esmf
#	@$(MAKE) wipeout_mapl
#	@$(MAKE) wipeout_fvdycore
#
#wipeout_esmf:
#	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
#	@echo '%%%%%  Wiping out the ESMF installation    %%%%%'
#	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
#	rm -f $(ESMF_DIR)/esmf.install
#	@$(MAKE) -C $(ESMF_DIR) distclean
#	rm -rf $(ESMF_DIR)/$(ARCH)
#	rm -f $(ESMF_DIR)/esmf.install
#
#wipeout_mapl:
#	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
#	@echo '%%%%%  Wiping out the MAPL installation    %%%%%'
#	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
#	rm -f $(ESMADIR)/mapl.install
#	rm -rf $(ESMADIR)/$(ARCH)
#	rm -f $(ESMADIR)/Config/bin/*.x
#	rm -f ./*___.*
#	@$(MAKE) -C $(ESMADIR) distclean
#	rm -f $(ESMADIR)/mapl.install
#
#wipeout_fvdycore:
#	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
#	@echo '%%%%%  Wiping out the FVDYCORE installation  %%%%%'
#	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
#	rm -f $(FVDIR)/fvdycore.install
#	@$(MAKE) -C $(FVDIR) ESMADIR=$(ESMADIR) distclean
#	rm -f $(FVDIR)/fvdycore.install
#
################################################################################
####                                                                         ###
####  Dependencies listing                                                   ###
####  (grep "USE " to get the list of module references!)                    ###
####                                                                         ###
####  From this list of dependencies, the "make" utility will figure out     ###
####  correct order of compilation (so we don't have to do that ourselves).  ###
####  This also allows us to compile on multiple processors with "make -j".  ###
####                                                                         ###
####  NOTES:                                                                 ###
####  (1) Only specify object-file dependencies that are within this         ###
####       directory.  Object files in other directories will be referenced  ### 
####       at link-time.                                                     ###
####  (2) For "make -jN" (i.e. compile N files simultaneously), all files    ###
####       in this directory must have a listed dependency.                  ###
####                                                                         ###
################################################################################
#
#GEOSChem.o		    : GEOSChem.F90                                  \
#                              GIGC_GridCompMod.o
#
#Chem_GridCompMod.o          : Chem_GridCompMod.F90                          \
#			      gigc_chunk_mod.o                              \
#                              gigc_historyexports_mod.o                     \
#                              gigc_providerservices_mod.o
#
#GIGC_GridCompMod.o          : GIGC_GridCompMod.F90                          \
#                              Chem_GridCompMod.o                            \
#	                      GEOS_ctmEnvGridComp.o
#
#GEOS_ctmEnvGridComp.o	    : GEOS_ctmEnvGridComp.F90
#
#gigc_chunk_mod.o            : gigc_chunk_mod.F90                            \
#                              gigc_historyexports_mod.o                     \
#
#gigc_historyexports_mod.o   : gigc_historyexports_mod.F90
#
#gigc_providerservices_mod.o : gigc_providerservices_mod.F90

#EOC

