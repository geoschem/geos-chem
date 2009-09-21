# $Id: Makefile,v 1.1 2009/09/21 18:53:21 bmy Exp $
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
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
SHELL   = /bin/sh

help:
	@echo 'Usage: make TARGET [ OPTIONAL-FLAGS ]'
	@echo ''
	@echo 'TARGET may be one of the following:'
	@echo 'all       Default target.  Synonym for libkpp lib exe'
	@echo 'lib       Only build GEOS-Chem objects & libraries in GeosCore subdir'
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
	@echo 'IPO=yes        Turns on optmization options -ipo -static (IFORT only)'
	@echo 'TRACEBACK=yes  Turns on traceback option (IFORT only)'
	@echo 'NTRAC=[43|54]  Specifies # of tracers for KPP chemical solver'
	@echo 'KPPSOLVER=___  Specifies the integrator used w/ KPP:'
	@echo '               Options: lsodes radau5 rosenbrock runge_kutta (default is lsodes)'