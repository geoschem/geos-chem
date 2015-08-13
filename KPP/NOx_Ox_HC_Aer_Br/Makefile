#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the KPP/NOx_Ox_HC_Aer_Br subdirectory)
#
# !DESCRIPTION: This makefile compiles the KPP solver code for the GEOS-Chem
#  43 tracer simulation (i.e. without secondary organic aerosol tracers).  
#  Object files (*.o) are bundled into the libKpp.a library (located in the 
#  LIB directory).  Module files (*.mod) are copied to the MOD directory. 
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
#  16 Sep 2009 - R. Yantosca - Initial version
#  21 Sep 2009 - R. Yantosca - Now call Makefile in help directory to
#                              display the help screen options
#  23 Nov 2009 - R. Yantosca - Now don't copy module files; they will be
#                              automatically written to the mod directory
#  11 Dec 2009 - R. Yantosca - Now get SHELL from Makefile_header.mk
#  04 Jun 2015 - R. Yantosca - Also remove *.mod, *.a files with "make clean"
#  04 Jun 2015 - R. Yantosca - Add debug target, remove help
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Define variables
ROOT    := ../..
HDR     := $(ROOT)/Headers
HELP    := $(ROOT)/help
LIB     := $(ROOT)/lib
MOD     := $(ROOT)/mod

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# List of source code files
SOURCES := $(wildcard gckpp*.F90)

# List of object files (replace .F and .F90 extensions with *.o)
OBJECTS := $(SOURCES:.F90=.o)

# List of module files.  Convert to lowercase, then prefix directory name.
MODULES :=$(OBJECTS:.o=.mod)
MODULES :=$(shell echo $(MODULES) | tr A-Z a-z)
MODULES :=$(foreach I,$(MODULES),$(MOD)/$(I))

# Library file
LIBRARY :=libKpp.a

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean debug

lib: $(OBJECTS)
	$(AR) crs $(LIBRARY) $(OBJECTS)
	mv $(LIBRARY) $(LIB)

clean:
	@echo "===> Making clean in directory: KPP/NOx_Ox_HC_Aer_Br <==="
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

gckpp_Function.o      : gckpp_Parameters.o

gckpp_Global.o        : gckpp_Parameters.o

gckpp_Hessian.o       : gckpp_Parameters.o      \
                        gckpp_HessianSP.o

gckpp_Initialize.o    : gckpp_Parameters.o      \
                        gckpp_Global.o          \
                        gckpp_Util.o            \
                        gckpp_Monitor.o

gckpp_Integrator.o    : gckpp_Parameters.o      \
                        gckpp_Global.o          \
                        gckpp_Function.o        \
                        gckpp_Rates.o           \
                        gckpp_Jacobian.o        \
                        gckpp_LinearAlgebra.o

gckpp_Jacobian.o      : gckpp_Parameters.o      \
                        gckpp_JacobianSP.o

gckpp_LinearAlgebra.o : gckpp_Parameters.o      \
                        gckpp_JacobianSP.o

gckpp_Model.o         : gckpp_Precision.o       \
                        gckpp_Parameters.o      \
                        gckpp_Global.o          \
                        gckpp_Function.o        \
                        gckpp_Integrator.o      \
                        gckpp_Rates.o           \
                        gckpp_Jacobian.o        \
                        gckpp_Hessian.o         \
                        gckpp_Stoichiom.o       \
                        gckpp_Monitor.o         \
                        gckpp_Util.o            \
                        gckpp_LinearAlgebra.o

gckpp_Parameters.o    : gckpp_Precision.o

#gckpp_Rates.o         : gckpp_Parameters.o      \
#                        gckpp_Global.o          \
#                        gckpp_Monitor.o         \
#                        gckpp_comode_mod.o

gckpp_Rates.o         : gckpp_Parameters.o      \
                        gckpp_Global.o          \
                        gckpp_Monitor.o

gckpp_Stoichiom.o     : gckpp_Parameters.o      \
                        gckpp_StoichiomSP.o

gckpp_Util.o          : gckpp_Parameters.o      \
                        gckpp_Global.o          \
                        gckpp_Monitor.o

#EOC
