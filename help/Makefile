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
#  See https://github.com/geoschem/geos-chem for complete history
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
	@echo ''
	@echo 'Usage: make -jN TARGET [ OPTIONAL-FLAGS ]'
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
	@echo 'libhistory      Builds GEOS-Chem objs & libs only in History/'
	@echo 'libiso          Builds GEOS-Chem objs & libs only in ISORROPIA/'
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
#-----------------------------------------------------------------------------
# Comment out until GTMM is restored.  GTMM hasn't worked for a long time.
# (bmy, 5/16/18)
#	@echo ''
#	@echo 'Special targets for mercury simulation:'
#	@echo '---------------------------------------'
#	@echo 'allhg           Default target for Hg simulation (synonym for "libhg exehg")'
#	@echo 'libhg           Builds GEOS-Chem code for Hg simulation'
#	@echo 'libgtmm         Builds GEOS-Chem + GTMM code for Hg simulation'
#	@echo 'exehg           Creates GEOS-Chem executable for Hg simulation'
#-----------------------------------------------------------------------------
	@echo ''
	@echo ''
	@echo '--------------------------------------------------------'
	@echo 'OPTIONAL-FLAGS may be one or more of the following:'
	@echo '--------------------------------------------------------'
	@echo ''
	@echo 'Parallelization and optimization flags:'
	@echo '---------------------------------------'
	@echo 'OMP=no          Turns off OpenMP parallelization'
	@echo '                --> Default: OMP=y'
	@echo ''
	@echo 'IPO=yes         Turns on optmization options -ipo -static (ifort only)'
	@echo '                --> Default: IPO=n'
	@echo ''
	@echo 'OPT=___         Specifies the optimization level'
	@echo '                --> Defaults OPT=-O2 (ifort)'
	@echo '                             OPT=-O3 (gfortran)'
	@echo ''
	@echo 'M_ARCH=___      Specifies the CPU architecture (gfortran only)'
	@echo ''
	@echo ''
	@echo 'Chemistry options:'
	@echo '------------------'
	@echo 'CHEM=___        Specifies which chemistry mechanism is used'
	@echo '                --> Options: Standard Tropchem SOA_SVPOA'
	@echo '                --> Default: Standard'
	@echo ''
	@echo 'KPPSOLVER=___   Specifies the integrator used with KPP'
	@echo '                --> Options: lsodes radau5 rosenbrock runge_kutta'
	@echo '                --> Default: rosenbrock'
	@echo ''
	@echo 'RRTMG=y         Turns on online radiative transfer using the RRTMG model'
	@echo '                --> Default: RRTMG=n'
	@echo ''
	@echo ''
	@echo 'Diagnostics options:'
	@echo '--------------------'
	@echo 'BPCH_DIAG=n     Disable binary-punch output for diagnostics'
	@echo '                --> Default: BPCH_DIAG=y'
	@echo ''
	@echo 'BPCH_TPBC=n     Disable binary-punch I/O for nested-grid boundary'
	@echo '                condition files for TPCORE.'
	@echo '                --> Default: BPCH_TPBC=y (if NEST is specified)'
	@echo ''
	@echo 'BPCH_TIMESER=n  Disable binary-punch output for timeseries diagnostics'
	@echo '                (i.e. to maintain backwards compatibility until these'
	@echo '                diagnostics are available in netCDF format)'
	@echo '                --> Default: BPCH_TIMESER=y'
	@echo ''
	@echo 'NC_NODEFLATE=y  Disable netCDF file compression when writing'
	@echo '                diagnostic or restart files to disk.  This may'
	@echo '                be necessary for debugging, in order to ensure'
	@echo '                that file compression does not cause problems'
	@echo '                when comparing netCDF files for identicality.'
	@echo '                --> Default: NC_NODEFLATE=n'
	@echo ''
	@echo ''
	@echo 'Debugging flags:'
	@echo '----------------'
	@echo 'BOUNDS=y        Turns on subscript-array checking'
	@echo '                --> Default: BOUNDS=n'
	@echo ''
	@echo 'DEBUG=y         Turns on options -g -O0 for running GEOS-Chem in a debugger'
	@echo '                --> Default: DEBUG=n'
	@echo ''
	@echo 'FPEX=y          Turns on checking for floating-point exceptions'
	@echo '                --> Default: FPEX=nn'
	@echo ''
	@echo 'FPE=y           Synonym for FPEX=y'
	@echo '                --> Default: FPE=n'
	@echo ''
	@echo 'TIMERS=___      Turn on GEOS-Chem timers (Use 1 for GEOS-Chem "Classic")'
	@echo '                --> Default: TIMERS=0'
	@echo ''
	@echo 'TRACEBACK=y     Print out a list of called routines if the run dies with an error'
	@echo '                --> Default: TRACEBACK=y'
	@echo '                --> To disable: TRACEBACK=n'
	@echo ''
	@echo 'GPROF=y         Compiles GEOS-Chem for use with the GNU profiler (aka gprof)'
	@echo '                --> Default: GPROF=n'
	@echo ''
	@echo 'TAU_PROF=y      Compiles GEOS-Chem for use with the TAU performance profiler'
	@echo '                --> Default: TAU_PROF=n'
	@echo ''
	@echo 'TAU_SF=___      Specifies the TAU selectfile for removing throttled files'
	@echo ''
	@echo ''
	@echo 'Aerosol microphysics flags:'
	@echo '---------------------------------------'
	@echo 'TOMAS=y         Turns on 30-bin TOMAS aerosol microphysics'
	@echo '                --> Default: TOMAS=n'
	@echo ''
	@echo 'TOMAS12=y       Turns on 12-bin TOMAS aerosol microphysics'
	@echo '                --> Default: TOMAS12=n'
	@echo ''
	@echo 'TOMAS15=y       Turns on 15-bin TOMAS aerosol microphysics'
	@echo '                --> Default: TOMAS15=n'
	@echo ''
	@echo 'TOMAS40=y       Turns on 40-bin TOMAS aerosol microphysics'
	@echo '                --> Default: TOMAS40=n'
	@echo ''
	@echo ''
	@echo 'GCHP (GEOS-Chem "High Peformance") flags:'
	@echo '---------------------------------------'
	@echo 'HPC=y           Compile GEOS-Chem for HPC environments (with ESMF & MAPL)'
	@echo '                --> Default: HPC=n'


