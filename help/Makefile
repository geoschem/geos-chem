# $Id: Makefile,v 1.4 2009/11/20 21:43:00 bmy Exp $
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the help subdirectory)
#
# !DESCRIPTION: Displays the makefile help screen for GEOS-Chem.
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make TARGET [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
#
# !REVISION HISTORY: 
#  21 Sep 2009 - R. Yantosca - Initial version
#  24 Sep 2009 - R. Yantosca - Added info about NONUMA option for PGI
#  24 Sep 2009 - R. Yantosca - Now list rosenbrock as default solver
#  19 Nov 2009 - R. Yantosca - Updated comments
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
SHELL   = /bin/sh

help:
	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	@echo '%%%      GEOS-Chem Help Screen      %%%'
	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	@echo
	@echo 'Usage: make TARGET [ OPTIONAL-FLAGS ]'
	@echo ''
	@echo 'TARGET may be one of the following:'
	@echo 'all       Default target.  Synonym for libkpp lib exe'
	@echo 'lib       Only build GEOS-Chem objects & libraries in GeosCore subdir'
	@echo 'libbpch   Only build GEOS-Chem objects & libraries in GeosBpch subdir'
	@echo 'libutil   Only build GEOS-Chem objects & libraries in GeosUtil subdir'
	@echo 'libkpp    Only build GEOS-Chem objects & libraries in KPP subdir'
	@echo 'exe       Build GEOS-Chem executable'
	@echo 'clean     Removes *.o, *.mod files in source code subdirs only'
	@echo 'realclean Removes *.o, *mod, *.lib *.a files source code, mod, lib, bin subdirs'
	@echo 'doc       Builds GEOS-Chem documentation (*.ps, *.pdf) with ProTeX in doc subdir'
	@echo 'docclean  Removes *.tex, *.pdf, *,ps from the doc subdir'
	@echo 'help      Displays this help screen'
	@echo ''
	@echo 'OPTIONAL-FLAGS may be:'
	@echo 'COMPILER=___   Options: ifort pgi sun xlf (default is ifort)'
	@echo 'DEBUG=yes      Compiles GEOS-Chem for use w/ a debugger'
	@echo 'BOUNDS=yes     Turns on subscript-array checking (for debugging purposes)'
	@echo 'OMP=[yes|no]   Turns OpenMP parallelization on/off (default is yes)'
	@echo 'IPO=yes        Turns on optmization options -ipo -static (ifort only)'
	@echo 'TRACEBACK=yes  Turns on -traceback option (ifort only)'
	@echo 'NONUMA=yes     Turns on -mp=nonuma option (pgi only)'
	@echo 'NTRAC=[43|54]  Specifies # of tracers for KPP chemical solver'
	@echo 'KPPSOLVER=___  Specifies the integrator used w/ KPP:'
	@echo '               Options: lsodes radau5 rosenbrock runge_kutta (default is rosenbrock)'