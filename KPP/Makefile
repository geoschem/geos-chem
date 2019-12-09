#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
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
# NTRAC      Cmd line argument; specifies either 43 or 54 tracer simulation
# KPPSOLVER  Cmd line argument; specifies the type of integrator to use
#
# !REVISION HISTORY:
#  16 Sep 2009 - R. Yantosca - Initial version
#  18 Sep 2009 - P. Le Sager - Added kppintegrator target & commented
#                              "make -C int" calls
#  20 Nov 2009 - P. Le Sager - Added CHEM option
#  23 Nov 2009 - R. Yantosca - Added realclean target
#  11 Dec 2009 - R. Yantosca - Now get SHELL from Makefile_header.mk
#  16 Oct 2013 - M. Sulprizio- Remove isoprene directory from realclean target.
#                              This directory has been removed because it is
#                              obsolete.
#  26 May 2015 - R. Yantosca - Bug fix, now point CHEM=benchmark by default
#  29 May 2015 - R. Yantosca - Add tropchem as a synonym for NOX_Ox_HC_Aer_Br
#  29 May 2015 - R. Yantosca - Now use KPP_CHEM instead of CHEM.  CHEM is an
#                              environment variable and can't be reset.
#  07 Dec 2015 - R. Yantosca - Restore fast "clean" command; add "slowclean"
#  07 Dec 2015 - R. Yantosca - Also add the slowrealclean target, which can
#                              be used for cleaning everything except RRTMG
#  12 May 2016 - R. Yantosca - Added "firstpass" target to compile Precision,
#                              Parameters, and Monitor modules first
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Define variables
ROOT := ..
DOC  := $(ROOT)/doc
HELP := $(ROOT)/help
LIB  := $(ROOT)/lib
MOD  := $(ROOT)/mod

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# Make the standard benchmark simulation if CHEM is not passed
# NOTE: In v10-01+, CHEM should always be passed from Makefile_header.mk!
ifndef KPP_CHEM
 KPP_CHEM :=Standard
endif

# Make rosenbrock the default solver
ifndef KPPSOLVER
 KPPSOLVER :=rosenbrock
endif

# solver (S=Source, T=Target)
SOLVER_SFILE :=./int/gckpp_Integrator_$(KPPSOLVER).F90
SOLVER_TFILE :=./$(KPP_CHEM)/gckpp_Integrator.F90

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: all       lib firstpass kppintegrator clean
.PHONY: realclean doc help      slowclean     slowrealclean

all: lib

lib: kppintegrator
	@$(MAKE) -C $(KPP_CHEM)

firstpass:
	@$(MAKE) -C $(KPP_CHEM) firstpass

kppintegrator:
#	@diff $(SOLVER_SFILE) $(SOLVER_TFILE) ;\
#	if [ $$? == 1 ] ; then   \
#		echo " copy $(SOLVER_SFILE) --> $(SOLVER_TFILE)";\
#		cp $(SOLVER_SFILE) $(SOLVER_TFILE) ;    \
#	fi

clean:
	@$(MAKE) -C $(KPP_CHEM) clean

slowclean:
	@$(MAKE) -C $(KPP_CHEM) slowclean

realclean:
	@$(MAKE) -C Standard      clean
	@$(MAKE) -C SOA_SVPOA     clean
	@$(MAKE) -C Tropchem      clean
	@$(MAKE) -C Custom        clean

slowrealclean:
	@$(MAKE) -C Standard      slowclean
	@$(MAKE) -C SOA_SVPOA     slowclean
	@$(MAKE) -C Tropchem      slowclean
	@$(MAKE) -C Custom        slowclean

help:
	@$(MAKE) -C $(HELP)
#EOC
