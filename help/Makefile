# $Id: Makefile,v 1.9 2010/03/15 19:33:18 ccarouge Exp $
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
#  23 Nov 2009 - R. Yantosca - Updated comments
#  11 Dec 2009 - R. Yantosca - Now get SHELL from Makefile_header.mk
#  21 Dec 2009 - R. Yantosca - Added info about HDF5 option
#  25 Jan 2010 - R. Yantosca - Added info about TOMAS option
#  10 Mar 2010 - C. Carouge  - Remove info about TOMAS option. Keep info
#                              about tomas target.
#EOP
#------------------------------------------------------------------------------
#BOC

# Get the Unix shell (in SHELL variable) from Makefile_header.mk
ROOTDIR = ..
include $(ROOTDIR)/Makefile_header.mk

help:
	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	@echo '%%%      GEOS-Chem Help Screen      %%%'
	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	@echo
	@echo 'Usage: make TARGET [ OPTIONAL-FLAGS ]'
	@echo ''
	@echo 'TARGET may be one of the following:'
	@echo 'all            Default target (synonym for "lib exe")'
	@echo 'lib            Builds all GEOS-Chem source code (synonym for "libkpp libutil libcore")'
	@echo 'libcore        Only builds GEOS-Chem objects & libraries in GeosCore subdir'
	@echo 'libutil        Only builds GEOS-Chem objects & libraries in GeosUtil subdir'
	@echo 'libkpp         Only builds GEOS-Chem objects & libraries in KPP subdir'
	@echo 'exe            Build GEOS-Chem executable'
	@echo 'clean          Removes *.o, *.mod files in source code subdirs only'
	@echo 'realclean      Removes all *.o, *mod, *.lib *.a, *.tex, *ps, *pdf files everywhere'
	@echo 'doc            Builds GEOS-Chem documentation (*.ps, *.pdf) with ProTeX in doc subdir'
	@echo 'docclean       Removes *.tex, *.pdf, *,ps from the doc subdir'
	@echo 'help           Displays this help screen'
	@echo ''
	@echo 'Special targets for TOMAS aerosol microphysics:'
	@echo 'tomas          Build GEOS-Chem + TOMAS code (synonym for "libtomas exetomas")'
	@echo 'libtomas       Only builds GEOS-Chem + TOMAS objects & libraries'
	@echo 'exetomas       Only builds GEOS-Chem + TOMAS executable'
	@echo 'cleantomas     Removes *.o *.mod files in GeosTomas directory only'
	@echo ''
	@echo 'Special target for mercury simulation:'
	@echo 'hg             Build GEOS-Chem + GTMM for mercury simulation'
	@echo ''
	@echo 'OPTIONAL-FLAGS may be:'
	@echo 'COMPILER=___   Options: ifort pgi sun xlf (default is ifort)'
	@echo 'HDF5=yes       Enables writing diagnostic timeseries output to HDF5 files'
	@echo 'DEBUG=yes      Compiles GEOS-Chem for use w/ a debugger'
	@echo 'BOUNDS=yes     Turns on subscript-array checking (for debugging purposes)'
	@echo 'OMP=[yes|no]   Turns OpenMP parallelization on/off (default is yes)'
	@echo 'IPO=yes        Turns on optmization options -ipo -static (ifort only)'
	@echo 'TRACEBACK=yes  Turns on -traceback option (ifort only)'
	@echo 'NONUMA=yes     Turns on -mp=nonuma option (pgi only)'
	@echo 'CHEM=___       Specifies which simulation is done. Options: standard SOA'
	@echo 'NTRAC=[43|54]  Specifies # of tracers for KPP chemical solver. Should use CHEM flag instead.'
	@echo 'KPPSOLVER=___  Specifies the integrator used w/ KPP:'
	@echo '               Options: lsodes radau5 rosenbrock runge_kutta (default is rosenbrock)'
	@echo '               (NOTE: This is set by default if you use "make tomas")'