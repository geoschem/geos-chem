# $Id: Makefile,v 1.2 2009/12/11 19:27:42 bmy Exp $
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the KPP/SOA subdirectory)
#
# !DESCRIPTION: This makefile compiles the KPP solver code for the GEOS-Chem
#  SOA simulation (with aromatic formation of secondary organic aerosol
#  tracers).  
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
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
ROOTDIR = ../..
HDR     = $(ROOTDIR)/Headers
HELP    = $(ROOTDIR)/help
LIB     = $(ROOTDIR)/lib
MOD     = $(ROOTDIR)/mod

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

#=============================================================================
# List of files to compile.  Here the order is not important,
# as we will explicity define the dependencies listing below.
#=============================================================================

# Source code files
SRC  = $(wildcard gckpp*.F90)

# Object files
OBJ  = $(SRC:.F90=.o)

#=============================================================================
# Makefile targets: type "make help" for a complete listing!
#=============================================================================

.PHONY: clean help

lib: $(OBJ)
	$(AR) crs libKpp.a $(OBJ)
	mv libKpp.a $(LIB)

clean:
	rm -f *.o *.mod geos

help:
	@$(MAKE) -C $(HELP)

#=============================================================================
# Dependencies listing (grep "USE " to get the list of module references!)
#
# From this list of dependencies, the "make" utility will figure out the
# correct order of compilation (so we don't have to do that ourselves!)
#=============================================================================

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

gckpp_Rates.o         : gckpp_Parameters.o      \
                        gckpp_Global.o          \
                        gckpp_Monitor.o

gckpp_Stoichiom.o     : gckpp_Parameters.o      \
                        gckpp_StoichiomSP.o

gckpp_Util.o          : gckpp_Parameters.o      \
                        gckpp_Global.o          \
                        gckpp_Monitor.o

#EOC
