# $Id: Makefile,v 1.1 2009/09/16 14:05:52 bmy Exp $
#==============================================================================
# ./KPP/Makefile: This sub-makefile compiles the KPP solver source code 
# and bundles all of the object files (*.o) into the /libKpp.a library (located
# in the LIB directory).  Module files (*.mod) are copied to the MOD directory. 
# (bmy, 6/18/09)
#
# To build the GEOS_Chem code, call "make" with the following syntax:
#
#   make TARGET [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% NOTE: Normally you will not have to call this Makefile directly,     %%%
# %%% it will be called automatically from the Makefile in the directory   %%%
# %%% just above this one!                                                 %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Makefile uses the following variables:
#
# Variable   Description
# --------   -----------
# CC         Contains the default C compilation commands (for PGI only)
# F90        Contains the Fortran compilation commands
# FREEFORM   Contains the command to force F90 "free format" compilation
# HDR        Specifies the directory where include files are found
# LD         Contains the command to link to libraries & make executable
# LIB        Specifies the directory where library files (*.a) are stored
# MOD        Specifies the directory where module files (*.mod) are stored
# OBJ        Specifies the list of object files (*.o) to be created.
# R8         Contains the command to force REAL -> REAL*8
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
#
# NOTE: CC, F90, FREEFORM, LD, R8 are included from "Makefile_header.mk".
#==============================================================================

# Define variables
SHELL  = /bin/sh
HDR    = ../Headers
LIB    = ../lib
MOD    = ../mod


# Include header file.  This returns variables CC, F90, FREEFORM, LD, R8,
# as well as the default Makefile compilation rules for source code files.
include ../Makefile_header.mk

# Make NTRAC=43 the default for KPP if it isn't passed!
ifndef NTRAC 
NTRAC = 43
endif

#=============================================================================
# List of files to compile.  Here the order is not important,
# as we will explicity define the dependencies listing below.
#=============================================================================

# Source code files
SRC  = $(wildcard gckpp_$(NTRAC)*.f90)

# # Object files
OBJ  = $(SRC:.f90=.o)

#=============================================================================
# Makefile targets: type "make help" for a complete listing!
#=============================================================================

.PHONY: clean help

lib: $(OBJ)
	cp *.mod $(MOD)
	$(AR) crs libKpp.a $(OBJ)
	mv libKpp.a $(LIB)
	@$(MAKE) -C int 

clean:
	@$(MAKE) -C int clean
	rm -f *.o *.mod geos

help:
	@echo 'Usage: make TARGET [ OPTIONAL-FLAGS ]'
	@echo ''
	@echo 'TARGET may be one of the following:'
	@echo 'lib    Only build GEOS-Chem objects & libraries in KPP subdir (default)'
	@echo 'clean  Removes *.o, *.mod files in KPP subdir only'
	@echo ''
	@echo 'OPTIONAL-FLAGS may be:'
	@echo 'COMPILER=___   Options: ifort pgi sun xlf (default is ifort)'
	@echo 'DEBUG=yes      Compiles GEOS-Chem for use w/ a debugger'
	@echo 'BOUNDS=yes     Turns on subscript-array checking (for debugging purposes)'
	@echo 'NTRAC=[43|54]  Specifies # of tracers for KPP chemical solver'
	@echo 'KPPSOLVER=___  Specifies the integrator used w/ KPP:'
	@echo '               Options: lsodes radau5 rosenbrock runge_kutta (default is lsodes)'	

#=============================================================================
# Dependencies listing (grep "USE " to get the list of module references!)
#
# From this list of dependencies, the "make" utility will figure out the
# correct order of compilation (so we don't have to do that ourselves!)
#=============================================================================

gckpp_$(NTRAC)t_Function.o      : gckpp_$(NTRAC)t_Parameters.o

gckpp_$(NTRAC)t_Global.o        : gckpp_$(NTRAC)t_Parameters.o

gckpp_$(NTRAC)t_Hessian.o       : gckpp_$(NTRAC)t_Parameters.o      \
                                  gckpp_$(NTRAC)t_HessianSP.o

gckpp_$(NTRAC)t_Initialize.o    : gckpp_$(NTRAC)t_Parameters.o      \
                                  gckpp_$(NTRAC)t_Global.o          \
                                  gckpp_$(NTRAC)t_Util.o            \
                                  gckpp_$(NTRAC)t_Monitor.o

gckpp_$(NTRAC)t_Integrator.o    : gckpp_$(NTRAC)t_Parameters.o      \
                                  gckpp_$(NTRAC)t_Global.o          \
                                  gckpp_$(NTRAC)t_Function.o        \
                                  gckpp_$(NTRAC)t_Rates.o           \
                                  gckpp_$(NTRAC)t_Jacobian.o        \
                                  gckpp_$(NTRAC)t_LinearAlgebra.o

gckpp_$(NTRAC)t_Jacobian.o      : gckpp_$(NTRAC)t_Parameters.o      \
                                  gckpp_$(NTRAC)t_JacobianSP.o

gckpp_$(NTRAC)t_LinearAlgebra.o : gckpp_$(NTRAC)t_Parameters.o      \
                                  gckpp_$(NTRAC)t_JacobianSP.o

gckpp_$(NTRAC)t_Model.o         : gckpp_$(NTRAC)t_Precision.o       \
                                  gckpp_$(NTRAC)t_Parameters.o      \
                                  gckpp_$(NTRAC)t_Global.o          \
                                  gckpp_$(NTRAC)t_Function.o        \
                                  gckpp_$(NTRAC)t_Integrator.o      \
                                  gckpp_$(NTRAC)t_Rates.o           \
                                  gckpp_$(NTRAC)t_Jacobian.o        \
                                  gckpp_$(NTRAC)t_Hessian.o         \
                                  gckpp_$(NTRAC)t_Stoichiom.o       \
                                  gckpp_$(NTRAC)t_Monitor.o         \
                                  gckpp_$(NTRAC)t_Util.o            \
                                  gckpp_$(NTRAC)t_LinearAlgebra.o

gckpp_$(NTRAC)t_Parameters.o    : gckpp_$(NTRAC)t_Precision.o

gckpp_$(NTRAC)t_Rates.o         : gckpp_$(NTRAC)t_Parameters.o      \
                                  gckpp_$(NTRAC)t_Global.o          \
                                  gckpp_$(NTRAC)t_Monitor.o

gckpp_$(NTRAC)t_Stoichiom.o     : gckpp_$(NTRAC)t_Parameters.o      \
                                  gckpp_$(NTRAC)t_StoichiomSP.o

gckpp_$(NTRAC)t_Util.o          : gckpp_$(NTRAC)t_Parameters.o      \
                                  gckpp_$(NTRAC)t_Global.o          \
                                  gckpp_$(NTRAC)t_Monitor.o
