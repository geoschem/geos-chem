#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
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
#   make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]
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
SRC  = $(wildcard *.F90)

# Object files
OBJ  = $(SRC:.F90=.o)

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

CleanupCASAarrays.o            : CleanupCASAarrays.F90 defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

GTMM.o                         : GTMM.F90              defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o        \
                                 dorestart_mod.o       input_gtmm_mod.o

GTMM_coupled.o                 : GTMM_coupled.F90      defineConstants.o     \
                                 defineArrays.o        dorestart_mod.o       \
	                         loadCASAinput.o       input_gtmm_mod.o      

HgOutForGEOS.o                 : HgOutForGEOS.F90      defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o        \
                                 CasaRegridModule.o

assignAgeClassToRunningPool.o  : assignAgeClassToRunningPool.F90             \
                                 defineConstants.o     loadCASAinput.o       \
                                 defineArrays.o          

assignRanPoolToAgeClass.o      : assignRanPoolToAgeClass.F90                 \
                                 defineConstants.o     loadCASAinput.o       \
                                 defineArrays.o          

defineArrays.o                 : defineArrays.F90      defineConstants.o

doFPARandLAI.o                 : doFPARandLAI.F90      defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doHerbCarbon.o                 : doHerbCarbon.F90      defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o
#	$(F90) -O0 -c $(FREEFORM) doHerbCarbon.F90

doHerbCarbonHg.o               : doHerbCarbonHg.F90    defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          
#	$(F90) -O0 -c $(FREEFORM) doHerbCarbonHg.F90

doHerbivory.o                  : doHerbivory.F90       defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doHgDeposition.o               : doHgDeposition.F90    defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doLatitude.o                   : doLatitude.F90        defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doLeafRootShedding.o           : doLeafRootShedding.F90 defineConstants.o    \
                                 loadCASAinput.o       defineArrays.o          

doMaxHg.o                      : doMaxHg.F90           defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doNPP.o                        : doNPP.F90             defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doOptimumTemperature.o         : doOptimumTemperature.F90 defineConstants.o  \
                                 loadCASAinput.o       defineArrays.o          

doPET.o                        : doPET.F90             defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doSoilMoisture.o               : doSoilMoisture.F90    defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

doTreeCarbon.o                 : doTreeCarbon.F90      defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          
#	$(F90) -O0 -c $(FREEFORM) doTreeCarbon.F90

doTreeCarbonHg.o               : doTreeCarbonHg.F90    defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          
#	$(F90) -O0 -c $(FREEFORM) doTreeCarbonHg.F90

getAgeClassBF.o                : getAgeClassBF.F90     defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

getFireParams.o                : getFireParams.F90     defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

getFuelWood.o                  : getFuelWood.F90       defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

getSoilMoistParams.o           : getSoilMoistParams.F90 defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

getSoilParams.o                : getSoilParams.F90     defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

input_gtmm_mod.o               : input_gtmm_mod.F90    defineConstants.o     \
                                 defineArrays.o

load_GC_data.o                 : load_GC_data.F90      defineConstants.o     \
                                 loadCASAinput.o       CasaRegridModule.o   

loadCASAinput.o                : loadCASAinput.F90     defineConstants.o     \
                                 defineArrays.o        CasaRegridModule.o

loadHgDeposition.o             : loadHgDeposition.F90  defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o        \
                                 CasaRegridModule.o

organizeAgeClasses.o           : organizeAgeClasses.F90 defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          

processData.o                  : processData.F90       defineConstants.o     \
                                 loadCASAinput.o       defineArrays.o          
#	$(F90) -O0 -c $(FREEFORM) processData.F90

sort_pick_veg.o                : sort_pick_veg.F90     defineConstants.o

dorestart_mod.o                : dorestart_mod.F90     defineConstants.o     \
                                 defineArrays.o

#EOC	
