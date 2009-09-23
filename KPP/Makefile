# $Id: Makefile,v 1.6 2009/09/23 14:36:36 bmy Exp $
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the KPP subdirectory)
#
# !DESCRIPTION: This is main "router" makefile for the KPP solver.  It 
#  compiles the KPP code for one of the following types of GEOS-Chem
#  simulations:
#
#  \begin{enumerate}
#  \item GEOS-Chem "standard" simulation (43 tracers)
#  \item GEOS-Chem "secondary organic aerosol" simulation (54 tracers)
#  \end{enumerate}
#
# The KPP code will be compiled using one of the following numerical solvers:
#
#  \begin{enumerate}
#  \item rosenbrock (This is the default option.)
#  \item lsodes
#  \item radau5
#  \item runge\_kutta
#  \end{enumerate}
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
# NTRAC      Cmd line argument; specifies either 43 or 54 tracer simulation
# KPPSOLVER  Cmd line argument; specifies the type of integrator to use
#
# !REVISION HISTORY: 
#  16 Sep 2009 - R. Yantosca - Initial version
#  18 Sep 2009 - P. Le Sager - Added kppintegrator target & commented
#                              "make -C int" calls 
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
SHELL   = /bin/sh
ROOTDIR = ..
DOC     = $(ROOTDIR)/doc
HDR     = $(ROOTDIR)/Headers
HELP    = $(ROOTDIR)/help
LIB     = $(ROOTDIR)/lib
MOD     = $(ROOTDIR)/mod

# Include header file.  This returns variables CC, F90, FREEFORM, LD, R8,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

# Make the 43-tracer simulation the default
ifndef NTRAC 
NTRAC = 43
endif

# Make rosenbrock the default solver
ifndef KPPSOLVER
KPPSOLVER = rosenbrock
endif

# solver (S=Source, T=Target)
SOLVER_SFILE=./int/gckpp_Integrator_$(KPPSOLVER).f90
SOLVER_TFILE=./$(NTRAC)t/gckpp_Integrator.f90


#=============================================================================
# Makefile targets: type "make help" for a complete listing!
#=============================================================================

.PHONY: all lib kppintegrator clean doc help

all: lib

lib: kppintegrator
	@$(MAKE) -C $(NTRAC)t

kppintegrator:
	@diff $(SOLVER_SFILE) $(SOLVER_TFILE) ;\
	if [ $$? == 1 ] ; then   \
		echo " copy $(SOLVER_SFILE) --> $(SOLVER_TFILE)";\
		cp $(SOLVER_SFILE) $(SOLVER_TFILE) ;    \
	fi

clean:
	@$(MAKE) -C $(NTRAC)t clean

help:
	@$(MAKE) -C $(HELP) 

#EOC	
