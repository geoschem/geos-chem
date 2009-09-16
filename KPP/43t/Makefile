# $Id: Makefile,v 1.1 2009/09/16 19:20:27 bmy Exp $
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the KPP/43t subdirectory)
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
#   make TARGET [ OPTIONAL-FLAGS ]
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
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
SHELL   = /bin/sh
ROOTDIR = ../..
HDR     = $(ROOTDIR)/Headers
LIB     = $(ROOTDIR)/lib
MOD     = $(ROOTDIR)/mod

# Include header file.  This returns variables CC, F90, FREEFORM, LD, R8,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

#=============================================================================
# List of files to compile.  Here the order is not important,
# as we will explicity define the dependencies listing below.
#=============================================================================

# Source code files
SRC  = $(wildcard gckpp*.f90)

# Object files
OBJ  = $(SRC:.f90=.o)

#=============================================================================
# Makefile targets: type "make help" for a complete listing!
#=============================================================================

.PHONY: clean help

lib: $(OBJ)
	cp *.mod $(MOD)
	$(AR) crs libKpp.a $(OBJ)
	mv libKpp.a $(LIB)

clean:
	rm -f *.o *.mod geos

help:
	@echo 'Usage: make TARGET [ OPTIONAL-FLAGS ]'
	@echo ''
	@echo 'TARGET may be one of the following:'
	@echo 'all    Default target.  Synonym for "lib".'
	@echo 'lib    Only build GEOS-Chem objects & libraries in KPP subdir (default)'
	@echo 'clean  Removes *.o, *.mod files in KPP subdir only'
	@echo 'help   Displays this help screen'
	@echo ''
	@echo 'OPTIONAL-FLAGS may be:'
	@echo 'COMPILER=___   Options: ifort pgi sun xlf (default is ifort)'
	@echo 'DEBUG=yes      Compiles GEOS-Chem for use w/ a debugger'
	@echo 'BOUNDS=yes     Turns on subscript-array checking (for debugging purposes)'
	@echo 'IPO=yes        Turns on optmization options -ipo -static (IFORT only)'
	@echo 'NTRAC=[43|54]  Specifies # of tracers for KPP chemical solver'
	@echo 'KPPSOLVER=___  Specifies the integrator used w/ KPP:'
	@echo '               Options: lsodes radau5 rosenbrock runge_kutta (default is rosenbrock)'

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
