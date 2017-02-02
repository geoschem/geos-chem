#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
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
#   make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]
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
#  08 Apr 2015 - M. Sulprizio- Updated to include make options added in v10-01
#  29 May 2015 - R. Yantosca - Updated help screen for v10-01 CHEM options
#  05 Dec 2016 - R. Yantosca - Updated for v11-01
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
	@echo 'Usage: make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]'
	@echo ''
	@echo '-jN             Compiles N files at a time (to reduce compilation time)'
	@echo ''
	@echo '--------------------------------------------------------'
	@echo 'TARGET may be one of the following:'
	@echo '--------------------------------------------------------'
	@echo 'all             Default target (synonym for "lib exe")'
	@echo 'lib             Builds GEOS-Chem source code'
	@echo 'libcore         Builds GEOS-Chem objs & libs only in GeosCore/'
	@echo 'libheaders      Builds GEOS-Chem objs & libs only in Headers/'
	@echo 'libiso          Builds GEOS-Chem objs & libs only in ISOROPIA/'
	@echo 'libkpp          Builds GEOS-Chem objs & libs only in KPP/'
	@echo 'libnc           Builds GEOS-Chem objs & libs only in NcdfUtil/'  
	@echo 'librad          Builds GEOS-Chem objs & libs only in GeosRad/'  
	@echo 'libutil         Builds GEOS-Chem objs & libs only in GeosUtil/'
	@echo 'ncdfcheck       Determines if the netCDF library installation works'
	@echo 'exe             Creates GEOS-Chem executable'
	@echo 'clean           Removes *.o, *.mod files in source code subdirs only'
	@echo 'realclean       Removes all *.o, *mod, *.lib *.a, *.tex, *ps, *pdf files everywhere'
	@echo 'distclean       Synonym for "make realclean"'
	@echo 'doc             Builds GEOS-Chem documentation (*.ps, *.pdf) in doc/'
	@echo 'docclean        Removes *.tex, *.pdf, *,ps from doc/'
	@echo 'tauclean        Removes *.pdb, *.inst, *.pp, and *.continue.* files produced by TAU'
	@echo 'help            Displays this help screen'
	@echo ''	       
	@echo 'Special targets for mercury simulation:'
	@echo '---------------------------------------'
	@echo 'allhg           Default target for Hg simulation (synonym for "libhg exehg")'
	@echo 'libhg           Builds GEOS-Chem code for Hg simulation'
	@echo 'libgtmm         Builds GEOS-Chem + GTMM code for Hg simulation'
	@echo 'exehg           Creates GEOS-Chem executable for Hg simulation'
	@echo ''	       
	@echo ''	       
	@echo '--------------------------------------------------------'
	@echo 'REQUIRED-FLAGS include:'
	@echo '--------------------------------------------------------'
	@echo 'MET=____        Specifies the met field type'
	@echo '                --> Options: merra2 merra geosfp geos5 geos4 gcap'
	@echo 'GRID=___        Specifies the horizontal grid'
	@echo '                --> Options: 4x5 2x25 05x0625 05x0666 025x03125'
	@echo ''
	@echo '--------------------------------------------------------'
	@echo 'OPTIONAL-FLAGS may be one or more of the following:'
	@echo '--------------------------------------------------------'
	@echo 'COMPILER=___    Specifies the compiler'
	@echo '                --> Options: ifort gfortran pgi (default is ifort)'
	@echo 'NEST=___        Specifies the nested-grid domain'
	@echo '                 --> Options: AS CH NA EU and CU (custom)'
	@echo 'NO_REDUCED=y    Compiles GEOS-Chem with the full vertical grid'
	@echo '                --> Default: NO_REDUCED=n (i.e use the reduced grid)'
	@echo ''	       
	@echo 'Parallelization and optimization flags:'
	@echo '---------------------------------------'
	@echo 'OMP=no          Turns off OpenMP parallelization'
	@echo '                --> Default: OMP=y'
	@echo 'NONUMA=yes      Turns on -mp=nonuma option (pgi only)'
	@echo '                --> Default: OMP=n'
	@echo 'IPO=yes         Turns on optmization options -ipo -static (ifort only)'
	@echo '                --> Default: IPO=n'
	@echo 'OPT=___         Specifies the optimization level'
	@echo '                --> Defaults OPT=-O2 (ifort and pgi)'
	@echo '                             OPT=-O3 (gfortran)'
	@echo 'M_ARCH=___      Specifies the CPU architecture (gfortran only)'
	@echo ''	       
	@echo 'Chemistry options:'
	@echo '------------------'
	@echo 'CHEM=___        Specifies which chemistry mechanism is used'
	@echo '                --> Options: Standard Tropchem UCX SOA SOA-SVPOA'
	@echo '                    Default: Standard'
	@echo 'KPPSOLVER=___   Specifies the integrator used with KPP'
	@echo '                --> Options: lsodes radau5 rosenbrock runge_kutta'
	@echo '                    Default: rosenbrock'
	@echo 'UCX=y           Turns on UCX chemistry mechanism; sets (CHEM=UCX NO_REDUCED=y)'
	@echo '                --> Default: UCX=n'
	@echo 'RRTMG=y         Turns on online radiative transfer using the RRTMG model'
	@echo '                --> Default: RRTMG=n'
	@echo ''	       
	@echo 'Debugging flags:'
	@echo '----------------'
	@echo 'BOUNDS=y        Turns on subscript-array checking'
	@echo '                --> Default: BOUNDS=n'
	@echo 'DEBUG=y         Turns on options -g -O0 for running GEOS-Chem in a debugger'
	@echo '                --> Default: DEBUG=n'
	@echo 'FPEX=y          Turns on checking for floating-point exceptions'
	@echo '                --> Default: FPEX=nn'
	@echo 'FPE=y           Synonym for FPEX=y'
	@echo '                --> Default: FPE=n'
	@echo 'TIMERS=___      Turn on GEOS-Chem timers (Use 1 for GEOS-Chem "Classic")'
	@echo '                --> Default: TIMERS=0'
	@echo 'TRACEBACK=y     Print out a list of called routines if the run dies with an error'
	@echo '                --> Default: TRACEBACK=y'
	@echo '                --> To disable: TRACEBACK=n'
	@echo 'GPROF=y         Compiles GEOS-Chem for use with the GNU profiler (aka gprof)'
	@echo '                --> Default: GPROF=n'
	@echo 'TAU_PROF=y      Compiles GEOS-Chem for use with the TAU performance profiler'
	@echo '                --> Default: TAU_PROF=n'
	@echo 'TAU_SF=___      Specifies the TAU selectfile for removing throttled files' 
	@echo ''	       
	@echo 'Aerosol microphysics flags:'
	@echo '---------------------------------------'
	@echo 'TOMAS=y         Turns on 30-bin TOMAS aerosol microphysics'
	@echo '                --> Default: TOMAS=n'
	@echo 'TOMAS12=y       Turns on 12-bin TOMAS aerosol microphysics'
	@echo '                --> Default: TOMAS12=n'
	@echo 'TOMAS15=y       Turns on 15-bin TOMAS aerosol microphysics'
	@echo '                --> Default: TOMAS15=n'
	@echo 'TOMAS40=y       Turns on 40-bin TOMAS aerosol microphysics'
	@echo '                --> Default: TOMAS40=n'
	@echo ''
	@echo 'GCHP (GEOS-Chem "High Peformance") flags:'
	@echo '---------------------------------------'
	@echo 'HPC=y           Compile GEOS-Chem for HPC environments (with ESMF & MAPL)'
	@echo '                --> Default: HPC=n'


