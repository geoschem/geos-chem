#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the GTMM subdirectory)
#
# !DESCRIPTION: This is main "router" makefile for the GTMM model.  It 
#  compiles the GTMM code for GEOS-Chem mercury simulations.
#
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make TARGET [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% NOTE: Normally you will not have to call this Makefile directly,     %%%
# %%% it will be called automatically from the main GEOS-Chem Makefile in  %%%
# %%% GeosCore directory!                                                  %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
# ROOTDIR    Specifies the root-level directory of the GEOS-Chem code
# DOC        Specifies the directory where GEOS-Chem documentation is found
# HDR        Specifies the directory where GEOS-Chem include files are found
# LIB        Specifies the directory where library files (*.a) are stored
# MOD        Specifies the directory where module files (*.mod) are stored
# AR         Sys var w/ name of library creator program (i.e., "ar", "ranlib")
# MAKE       Sys var w/ name of Make command (i.e, "make" or "gmake")
#
# !REVISION HISTORY: 
#  16 Sep 2009 - R. Yantosca - Initial version
#  18 Sep 2009 - P. Le Sager - Added kppintegrator target & commented
#                              "make -C int" calls
#  21 Sep 2009 - C. Carouge  - Adapted to use with GTMM model. 
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
SHELL   = /bin/sh
ROOTDIR = ..
HDR     = $(ROOTDIR)/Headers
HELP    = $(ROOTDIR)/help
LIB     = $(ROOTDIR)/lib
MOD     = $(ROOTDIR)/mod
LGTMM   = -L$(LIB) -lHg


# Include header file.  This returns variables CC, F90, FREEFORM, LD, R8,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

#=============================================================================
# List of files to compile.  Here the order is not important,
# as we will explicity define the dependencies listing below.
#=============================================================================

# Source code files
SRC  = $(wildcard *.f90)

# Object files
OBJ  = $(SRC:.f90=.o)

#=============================================================================
# Makefile targets: type "make help" for a complete listing!
#=============================================================================

.PHONY: clean help gtmm

lib: $(OBJ)
	$(AR) crs libHg.a $(OBJ)
	mv libHg.a $(LIB)

gtmm:
	@$(MAKE) lib
	@$(MAKE) exe

exe: 
	$(LD) $(OBJ) $(LGTMM) -o gtmm

clean:
	rm -f *.o *.mod gtmm

help:
	@$(MAKE) -C $(HELP)

#=============================================================================
# Dependencies listing (grep "USE " to get the list of module references!)
#
# From this list of dependencies, the "make" utility will figure out the
# correct order of compilation (so we don't have to do that ourselves!)
#=============================================================================

CleanupCASAarrays.o            : CleanupCASAarrays.f90 defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

GTMM.o                         : GTMM.f90              defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o        \
                                 dorestart_mod.o       input_gtmm_mod.o

GTMM_coupled.o                 : GTMM_coupled.f90      defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o        \
                                 HgdataforGC_mod.o

HgdataforGC_mod.o              : HgdataforGC_mod.f90   defineConstants.o     \
                                 defineArrays.o

HgOutForGEOS.o                 : HgOutForGEOS.f90      defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o        \
                                 CasaRegridModule.o

assignAgeClassToRunningPool.o  : assignAgeClassToRunningPool.f90             \
                                 defineConstants.o     loadCASAinput.o       \
                                 defineArrays.o          

assignRanPoolToAgeClass.o      : assignRanPoolToAgeClass.f90                 \
                                 defineConstants.o     loadCASAinput.o       \
                                 defineArrays.o          

defineArrays.o                 : defineArrays.f90      defineConstants.o

doFPARandLAI.o                 : doFPARandLAI.f90      defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doHerbCarbon.o                 : doHerbCarbon.f90      defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o
#	$(F90) -O0 -c $(FREEFORM) doHerbCarbon.f90

doHerbCarbonHg.o               : doHerbCarbonHg.f90    defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          
#	$(F90) -O0 -c $(FREEFORM) doHerbCarbonHg.f90

doHerbivory.o                  : doHerbivory.f90       defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doHgDeposition.o               : doHgDeposition.f90    defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doLatitude.o                   : doLatitude.f90        defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doLeafRootShedding.o           : doLeafRootShedding.f90 defineConstants.o    \
                                 loadCASAinput.o       defineArrays.o          

doMaxHg.o                      : doMaxHg.f90           defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doNPP.o                        : doNPP.f90             defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doOptimumTemperature.o         : doOptimumTemperature.f90 defineConstants.o  \
                                 loadCASAinput.o       defineArrays.o          

doPET.o                        : doPET.f90             defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doSoilMoisture.o               : doSoilMoisture.f90    defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doTreeCarbon.o                 : doTreeCarbon.f90      defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          
#	$(F90) -O0 -c $(FREEFORM) doTreeCarbon.f90

doTreeCarbonHg.o               : doTreeCarbonHg.f90    defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          
#	$(F90) -O0 -c $(FREEFORM) doTreeCarbonHg.f90

getAgeClassBF.o                : getAgeClassBF.f90     defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

getFireParams.o                : getFireParams.f90     defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

getFuelWood.o                  : getFuelWood.f90       defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

getSoilMoistParams.o           : getSoilMoistParams.f90 defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

getSoilParams.o                : getSoilParams.f90     defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

input_gtmm_mod.o               : input_gtmm_mod.f90    defineConstants.o     \
                                 defineArrays.o

loadCASAinput.o                : loadCASAinput.f90     defineConstants.o     \
                                 defineArrays.o        CasaRegridModule.o

loadHgDeposition.o             : loadHgDeposition.f90  defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o        \
                                 CasaRegridModule.o

organizeAgeClasses.o           : organizeAgeClasses.f90 defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

processData.o                  : processData.f90       defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          
#	$(F90) -O0 -c $(FREEFORM) processData.f90

sort_pick_veg.o                : sort_pick_veg.f90     defineConstants.o

dorestart_mod.o                : dorestart_mod.f90     defineConstants.o     \
                                 defineArrays.o

#EOC	
