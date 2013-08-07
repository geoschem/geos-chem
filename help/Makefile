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
#  26 Aug 2011 - R. Yantosca - Added info about APM targets
#  26 Aug 2011 - R. Yantosca - Add info about the PRECISE=no option
#  11 May 2012 - R. Yantosca - Updated to include info about new make options
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
	@echo 'lib            Builds GEOS-Chem source code'
	@echo 'libcore        Builds GEOS-Chem objs & libs only in GeosCore/'
	@echo 'libheaders     Builds GEOS-Chem objs & libs only in GeosHeaders/'
	@echo 'libiso         Builds GEOS-Chem objs & libs only in ISOROPIA/'
	@echo 'libkpp         Builds GEOS-Chem objs & libs only in KPP/'
	@echo 'libnc          Builds GEOS-Chem objs & libs only in NcdfUtil/'  
	@echo 'libutil        Builds GEOS-Chem objs & libs only in GeosUtil/'
	@echo 'ncdfcheck      Determines if the netCDF library installation works'
	@echo 'exe            Creates GEOS-Chem executable'
	@echo 'clean          Removes *.o, *.mod files in source code subdirs only'
	@echo 'realclean      Removes all *.o, *mod, *.lib *.a, *.tex, *ps, *pdf files everywhere'
	@echo 'distclean      Synonym for "make realclean"'

	@echo 'doc            Builds GEOS-Chem documentation (*.ps, *.pdf) in doc/'
	@echo 'docclean       Removes *.tex, *.pdf, *,ps from doc/'
	@echo 'help           Displays this help screen'
	@echo ''
	@echo 'Special targets for mercury simulation:'
	@echo 'allhg          Default target for Hg simulation (synonym for "libhg exehg")'
	@echo 'libhg          Builds GEOS-Chem code for Hg simulation'
	@echo 'libgtmm        Builds GEOS-Chem + GTMM code for Hg simulation'
	@echo 'exehg          Creates GEOS-Chem executable for Hg simulation'
	@echo ''
	@echo 'Special targets for TOMAS aerosol microphysics:'
	@echo 'tomas          Builds GEOS-Chem + TOMAS (synonym for "libtomas exetomas")'
	@echo 'libtomas       Builds GEOS-Chem + TOMAS objs & libs in GeosTomas/'
	@echo 'exetomas       Creates GEOS-Chem + TOMAS executable'
	@echo 'cleantomas     Removes *.o *.mod files only in GeosTomas/'
	@echo ''
	@echo 'Special targets for APM aerosol microphysics:'
	@echo 'apm            Builds GEOS-Chem + APM (synonym for libapm exeapm)'
	@echo 'libapm         Builds GEOS-Chem + APM objs & libs in GeosApm/ subdir'
	@echo 'exeapm         Creates GEOS-Chem + APM executable in GeosApm/'
	@echo 'cleanapm       Removes *.o *.mod files only in GeosApm/'
	@echo ''
	@echo 'OPTIONAL-FLAGS may be:'
	@echo 'COMPILER=___   Options: ifort pgi (default is ifort)'
	@echo 'PRECISE=no     Disable precise floating point math optimization (for speed)'
	@echo 'DEBUG=yes      Builds GEOS-Chem for a debugger (with -g -O0)'
	@echo 'BOUNDS=yes     Turns on subscript-array checking (for debug)'
	@echo 'OMP=[yes|no]   Turns OpenMP parallelization on/off (default is yes)'
	@echo 'IPO=yes        Turns on optmization options -ipo -static (ifort only)'
	@echo 'TRACEBACK=yes  Turns on -traceback option (ifort only)'
	@echo 'NONUMA=yes     Turns on -mp=nonuma option (pgi only)'
	@echo 'CHEM=___       Specifies which simulation is done. Options: standard SOA'
	@echo 'NTRAC=[43|54]  Specifies # of tracers for KPP chemical solver. Should use CHEM flag instead.'
	@echo 'KPPSOLVER=___  Specifies the integrator used w/ KPP:'
	@echo '               Options: lsodes radau5 rosenbrock runge_kutta (default is rosenbrock)'
	@echo '               (NOTE: This is set by default if you use "make tomas")'