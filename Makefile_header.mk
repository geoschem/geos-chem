#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: Makefile_header.mk
#
# !DESCRIPTION: This sub-makefile defines the variables which specify
# compilation options for the different supported compiler/platform
# combinations.  Also, the default makefile compilation rules are specified 
# here.
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#                                                                             .
# The following variables are exported to the main-level Makefile:
#                                                                             .
# Variable   Description
# --------   -----------
# CC         Contains the default C compilation commands (for PGI only)
# F90        Contains the Fortran compilation commands
# FREEFORM   Contains the command to force F90 "free format" compilation
# LD         Contains the command to link to libraries & make executable
# LINK       Contains the commands to link to GEOS-Chem built libraries
# R8         Contains the command to force REAL -> REAL*8
# SHELL      Contains the default Unix shell to use when building code
# NCL        Contains the default netCDF library link commands
#                                                                             .
# FFLAGS is a local variable that is not returned to the "outside world", 
# but is only used locally.  COMPILER, HDF5, and OMP are all input via the
# command line or via environment variables.
#                                                                             .
# NOTE: We now use SHELL :=/bin/bash as the default Unix shell.  This allows
# us to extend the Makefile ifeq statements so that we can test for more than
# one string.  The following example is used to ensure that the met field name
# selected by the user is case-insensitive:
# 
#  # %%%%% GEOS-FP %%%%%
#  REGEXP             :=(^[Gg][Ee][Oo][Ss][Ff][Pp])|(^[Gg][Ee][Oo][Ss].[Ff][Pp])
#  ifeq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
#    USER_DEFS        += -DGEOS_FP
#  endif
#                                                                             .
# The [[ ]] in bash is an evaluation.  The above ifeq statement uses regular
# expressions to test if the MET variable matches the string "GEOS" (case-
# insensitive) and either "FP" or "any character and then a FP".  This will
# return true (via the "echo true" statement) for combinations like "GEOS-FP", 
# "geosfp", "Geos-FP", "GeOs.FP", etc.  This is a robust way of evaluating
# the user's input, and will make errors less likely.
#
# !REVISION HISTORY: 
#  16 Sep 2009 - R. Yantosca - Initial version
#  22 Sep 2009 - R. Yantosca - Bug fix, added -I$(HDR) to F90 compilation lines
#  24 Sep 2009 - R. Yantosca - added NONUMA option for PGI compiler
#  07 Oct 2009 - R. Yantosca - Replaced .SUFFIXES section w/ pattern rules
#  19 Nov 2009 - R. Yantosca - Now use OMP variable to determine whether to
#                              turn on OpenMP parallelization options 
#  23 Nov 2009 - R. Yantosca - Now use -module $(MOD) instead of -I$(MOD) to 
#                              specify the directory for *.mod files on both
#                              IFORT and PGI compilers.
#  23 Nov 2009 - R. Yantosca - Now use -moddir=$(MOD) and -M$(MOD) instead of
#                              -I$(MOD) to specify the directory for *.mod 
#                              files on the SunStudio compiler.
#  23 Nov 2009 - R. Yantosca - Change DEBUG to allow for new version of 
#                              Totalview which doesn't choke when debugging
#                              parallel code (Totalview 8.6.1-1)
#  02 Dec 2009 - R. Yantosca - Added SUN32 switch for building 32-bit 
#                              executbable on the SunStudio compiler
#  11 Dec 2009 - R. Yantosca - Now define SHELL here and export to other 
#                              Makefiles, so as to have a single place where
#                              the Unix shell name is defined.
#  21 Dec 2009 - R. Yantosca - Add H5I and H5L variables to specify the
#                              HDF5 library and include paths.  Also set
#                              the default to not link to the HDF5 libraries.
#  21 Dec 2009 - R. Yantosca - Now pass LINK back to the outside world, so
#                              that the Makefile that builds the executable
#                              can reference it.
#  19 Jan 2010 - R. Yantosca - Minor fix, add -m64 if SUN32 is not defined.
#  25 Jan 2010 - R. Yantosca - Now add -DTOMAS to FFLAGS if necessary
#  28 Jan 2010 - C. Carouge  - Add -lIsoropia to LINK, for ISORROPIA II
#  16 Feb 2011 - R. Yantosca - Now add -DAPM to FFLAGS if necessary
#  25 Aug 2011 - R. Yantosca - Add "-fp-model source" to FFLAGS for IFORT 
#                              compiler.  This will prevent aggressive 
#                              optimizations from changing numerical results.
#  25 Aug 2011 - R. Yantosca - Add -CU (check for uninit'd variables) to 
#                              FFLAGS when using IFORT w/ the DEBUG option.
#  26 Aug 2011 - R. Yantosca - Allow for deactivation of the "-fp-model source"
#                              option by using the PRECISE=no env variable
#  24 Jan 2012 - R. Yantosca - If NETCDF=yes, GEOS-Chem will link and include
#                              to the netCDF dir paths that are specified
#  24 Jan 2012 - R. Yantosca - Now use := for makefile assignment statements
#  10 Feb 2012 - R. Yantosca - When compiling with NETCDF=yes or HDF5=yes,
#                              we must also add the flags -mcmodel=medium 
#                              -i-dynamic to FFLAGS in order to avoid memory 
#                              errors (for IFORT only)
#  10 Feb 2012 - R. Yantosca - Remove -CU from the DEBUG option (IFORT only)
#  19 Mar 2012 - R. Yantosca - Add optional NO_ISO switch, which will turn off
#                              the ISORROPIA ATE package for testing
#  05 Apr 2012 - R. Yantosca - Now assume netCDF is always used
#  05 Apr 2012 - R. Yantosca - Change BL_INC_NETCDF to INC_NETCDF
#  05 Apr 2012 - R. Yantosca - Change BL_INC_HDF5   to INC_HDF5
#  05 Apr 2012 - R. Yantosca - Change BL_LIB_NETCDF to LIB_NETCDF
#  05 Apr 2012 - R. Yantosca - Change BL_LIB_HDF5   to LIB_HDF5
#  30 Apr 2012 - R. Yantosca - Add NETCDF3=[yes|no] makefile option
#  30 Apr 2012 - R. Yantosca - Use separate netCDF link and include paths
#                              for netCDF3 and for netCDF4
#  30 Apr 2012 - R. Yantosca - Also add -mcmodel=medium flag for PGI compiler
#  09 May 2012 - R. Yantosca - Now try to get the proper linking sequence 
#                              for netCDF etc w/ nf-config and nc-config.
#  11 May 2012 - R. Yantosca - Now export NCL (netCDF linking sequence)
#  17 Aug 2012 - R. Yantosca - Now add RRTMG=yes option for RRTMG rad transfer
#  07 Sep 2012 - R. Yantosca - Now add OPT variable to set global opt levels
#  07 Sep 2012 - R. Yantosca - Also set TRACEBACK for PGI compiler
#  17 Apr 2013 - R. Yantosca - Add switch to set -DKPP_SOLVE_ALWAYS, which 
#                              will force KPP to get past nonconvergences
#  25 Feb 2013 - S. Farina   - Add flag for TOMAS40
#  22 Apr 2013 - R. Yantosca - TOMAS40=yes option now sets -DTOMAS -DTOMAS40
#  28 Apr 2013 - S. Farina   - Add flags for TOMAS15 and TOMAS12
#  13 Aug 2013 - R. Yantosca - Removed "define.h"; now set all GEOS-Chem
#                              user options via the Make command
#  14 Aug 2013 - R. Yantosca - Now use regular expressions to test the
#                              validity of command-line inputs
#  21 Aug 2013 - R. Yantosca - Improved error checking for command line inputs
#  26 Aug 2013 - R. Yantosca - Add -debug all as an IFORT debugging option
#  16 Sep 2013 - R. Yantosca - Now set GIGC Cpp switches first.  This allows
#                              us to skip the GRID setting if we are using
#                              EXTERNAL_GRID=yes or EXTERNAL_FORCING=yes.
#  18 Sep 2013 - M. Long     - Add edits for HPC Grid-Indpendent GEOS-Chem
#  26 Sep 2013 - R. Yantosca - MET=geosfp now sets Cpp switch w/ -DGEOS_FP
#  07 Nov 2013 - R. Yantosca - NEST=se to now sets CPP switch w/ -DNESTED_SE
#  08 Nov 2013 - R. Yantosca - Add FPEX flag to avoid conflicting with the
#                              ESMF/MAPL environment variable FPE
#  24 Feb 2014 - R. Yantosca - Add UCX=yes flag for invoking UCX strat chem
#  18 Mar 2014 - R. Yantosca - Now add TAU_PROF=y flag to invoke TAU profiler
#  19 Mar 2014 - R. Yantosca - Move library link commands after the sections
#                              that set the C-preprocessor switches
#  19 Mar 2014 - R. Yantosca - Restore GTMM compilation funcitonality
#  19 Mar 2014 - R. Yantosca - Add more visible comment section dividers
#  20 Mar 2014 - R. Yantosca - Bug fix: "+= -DDEBUG" instead of ":= -DDEBUG"
#  09 Jul 2014 - R. Yantosca - Now don't require MET or GRID if target is
#                              srcdoc, utildoc, gtmmdoc, makedoc, or hemcodoc
#  21 Jul 2014 - R. Yantosca - Update build sequence
#  03 Oct 2014 - R. Yantosca - Now turn on NO_REDUCED=y for hpc target
#  03 Oct 2014 - R. Yantosca - Now compatible with netCDF 4.1.1 or 4.2+
#  17 Oct 2014 - R. Yantosca - Don't require MET or GRID to remove ESMF etc.
#  05 Nov 2014 - R. Yantosca - Will compile w/ 8-byte precision by default
#  14 Nov 2014 - R. Yantosca - Further updates for hpc compilation
#  21 Nov 2014 - R. Yantosca - Add special compilation command for ISORROPIA
#  21 Nov 2014 - R. Yantosca - Add cosmetic changes and indentation 
#  06 Jan 2015 - R. Yantosca - Add two-way nesting options from Y. Y. Yan
#  09 Jan 2015 - M. Sulprizio- Now properly link to the RRTMG directory
#  13 Jan 2015 - R. Yantosca - Add fix for GEOS-Chem-Libraries library path
#  08 Apr 2015 - R. Yantosca - Bug fix: set RRTMG=yes if it passes the regexp
#  10 Apr 2015 - R. Yantosca - Export RRTMG_NEEDED var to be used elsewhere
#  10 Apr 2015 - R. Yantosca - Bug fix: -l rad should be -lrad in link var
#  12 May 2015 - R. Yantosca - Bug fix for PGI compiler: remove extra "-"
#                              in front of $(NC_INC_CMD) in the PGI section
#  12 May 2015 - R. Yantosca - Now use GC_BIN, GC_INCLUDE to point to the
#                              netCDF library paths and GC_F_BIN, GC_F_INCLUDE
#                              to point to netCDF-Fortran library paths.
#                              (In some cases, these are the same).
#  20 May 2015 - R. Yantosca - Test if GC_F_BIN and GC_F_INCLUDE are defined
#                              as env variables before trying to use them.
#  29 May 2015 - R. Yantosca - Now set KPP_CHEM for KPP.  We can't redefine
#                              the CHEM variable because it is an env var.
#  04 Jun 2015 - R. Yantosca - Now use RRTMG_NO_CLEAN=y or RRTMG_NOCLEAN=y to 
#                              removing RRTMG objects, modules, and libraries.
#  04 Jun 2015 - R. Yantosca - Bug fix: don't turn on UCX except for CHEM=UCX
#  15 Jun 2015 - R. Yantosca - Now define the HEMCO standalone link command
#                              separately from the GEOS-Chem link command
#  07 Jul 2015 - M. Sulprizio- Add option for CHEM=SOA_SVPOA
#  17 Jul 2015 - E. Lundgren - Remove BSTATIC option when picking pgi options 
#                              for debug run or regular run 
#  30 Jul 2015 - M. Yannetti - Added TIMERS.
#  03 Aug 2015 - M. Sulprizio- NEST=cu to now sets CPP switch w/ -DNESTED_CU for
#                              custom nested grids
#  11 Aug 2015 - R. Yantosca - Add MERRA2 as a met field option
#  24 Aug 2015 - R. Yantosca - Bug fix: Add missing | when testing USER_DEFS
#  07 Dec 2015 - R. Yantosca - Add "realclean_except_rrtmg" target that
#                              replaces the RRTMG_CLEAN variabe
#  10 Feb 2016 - E. Lundgren - Add BPCH restart file input and output switches
#  11 Feb 2016 - E. Lundgren - Change BPCH to BPCH_DIAG, NETCDF to NC_DIAG
#  12 Jul 2016 - E. Lundgren - Remove binary punch restart file option
#  19 Jul 2016 - R. Yantosca - Add more flags for enabling experimental code
#  20 Sep 2016 - M. Sulprizio- Remove NEST=se option. This grid was never fully
#                              implemented.
#  12 Dec 2016 - R. Yantosca - Allow gfortran etc. to compile with TAU_PROF=y
#  13 Dec 2016 - R. Yantosca - Add GPROF=y to compile for GNU profiler gprof
#  01 Mar 2017 - R. Yantosca - Bug fix: Make sure NO_REDUCED=no works
#  01 Mar 2017 - R. Yantosca - Set -DNC_HAS_COMPRESSION if the netCDF library
#                              can write compressed data to disk
#  07 Mar 2017 - R. Yantosca - Replace makefile variable COMPILER with
#                              COMPILER_FAMILY; also works if FC=mpif90
#  08 May 2017 - R. Yantosca - Add minor fixes to avoid Perl bareword errors
#  23 May 2017 - R. Yantosca - use -dumpversion to get the Gfortran version #
#  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
#  03 Jan 2018 - M. Sulprizio- Remove UCX flag. We now solely use Input_Opt%LUCX
#                              throughout GEOS-Chem.
#  07 Aug 2018 - R. Yantosca - For now, don't compile TOMAS/ APM when NC_DIAG=y
#  21 Aug 2018 - R. Yantosca - Simplify testing for netCDF-Fortran 
#  23 Aug 2018 - H.P. Lin    - Add NO_EXE=y to inhibit "geos" executable build
#                              and build libGeosCore.a instead for coupled
#                              models driving GEOS-Chem externally (by calling
#                              its libraries)
#  28 Aug 2018 - M. Sulprizio- Export EXE_NEEDED to be used in GeosCore/Makefile
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Set the default Unix shell and some error message variables            ###
###                                                                         ###
###############################################################################

# Set default shell to bash, for use with the Makefile conditionals
SHELL                :=/bin/bash

# Error message for bad COMPILER input
ERR_CMPLR            :="Unknown Fortran compiler!  Must be one of ifort, gfortran, pgfortran|pgi|pgf90, or mpifort|mpif90.  Check the FC environment variable in your .bashrc or .cshrc file."

# Error message for unknown compiler/OS combintation
ERR_OSCOMP           :="Makefile_header.mk not set up for this compiler/OS combination"

# Error message for bad MET input
ERR_MET              :="Select a met field: MET=geosfp, MET=merra2"

# Error message for bad GRID input
ERR_GRID             :="Select a horizontal grid: GRID=4x5. GRID=2x25, GRID=05x0625, GRID=025x03125"

# Error message for bad NEST input
ERR_NEST             :="Select a nested grid: NEST=as, NEST=ch, NEST=eu, NEST=na, NEST=cu"

# Error message for bad two-way coupled model input (yanyy,6/18/14)
ERR_COUPLECH         :="Select a coupled grid for China/SE Asia: COUPLECH=2x25ch, COUPLECH=4x5ch"
ERR_COUPLENA         :="Select a coupled grid for North America: COUPLENA=2x25na, COUPLENA=4x5na"
ERR_COUPLEEU         :="Select a coupled grid for Europe       : COUPLEEU=2x25eu, COUPLEEU=4x5eu"
ERR_COUPLE           :="Select a coupled choice: COUPLE=yes"

# Error message for bad GIGC config
ERR_GIGC             :="Unable to find the GIGC configuration file. Have you downloaded the GIGC?"

# Error message for TOMAS error message
ERR_MICPHYS           :="At present, microphysics packages (TOMAS, APM) cannot be used when NC_DIAG=y!"

###############################################################################
###                                                                         ###
###  Set C-preprocessor switches representing user options.  These are not  ###
###  specific to any compiler, but are general options for the simulation.  ###
###                                                                         ###
###  NOTE: To make the user input more robust, we use regular expression    ###
###  syntax to match characters in the various Makefile variables.  See     ###
###  this web page for more info:                                           ###
###  http://www.tldp.org/LDP/abs/html/x17046.html                           ###
###                                                                         ###
###############################################################################

#------------------------------------------------------------------------------
# Compiler settings
#------------------------------------------------------------------------------

# %%%%% OpenMP parallelization (on by default) %%%%%
ifndef OMP
  OMP                :=yes
endif

# Option to turn off OpenMP for testing
REGEXP               :=(^[Nn]|^[Oo])
ifeq ($(shell [[ "$(OMP)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DNO_OMP
endif

# %%%%% Set the HPC variable if we are building for use w/ ESMF/MPI %%%%
ifeq ($(shell [[ "$(MAKECMDGOALS)" =~ "hpc" ]] && echo true),true)
  HPC                :=yes
  export HPC
endif

# %%%%% For HPC, we disable OpenMP and turn on the full vertical grid %%%%%
REGEXP               := (^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(HPC)" =~ $(REGEXP) ]] && echo true),true)
  IS_HPC             :=1
  OMP                :=no
  NO_REDUCED         :=yes
# PRECISION          :=4
else
  IS_HPC             :=0
endif

# %%%%% Default to 8-byte precision unless specified otherwise %%%%%
ifndef PRECISION
 PRECISION           :=8
endif

# %%%%% Default to Timers disabled %%%%%
ifndef TIMERS
 TIMERS              :=0
endif

# %%%%% Turn on traceback (error stack report) by default %%%%%
ifndef TRACEBACK
 TRACEBACK           :=yes
endif

# %%%%% Set default compiler %%%%%

# %%%%% Test if mpif90/mpifort is selected (for now assume ifort) %%%%%
REGEXP               :=(^[Mm][Pp][Ii])
ifeq ($(shell [[ "$(FC)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DLINUX_IFORT
  REG_GNU            :=(^[Gg][Nn][Uu])
  REG_INTEL          :=(^[Ii][Ff][Oo][Rr][Tt])
  DISCRIM            :=$(word 1,$(shell $(FC) --version ) )
  ifeq ($(shell [[ "$(DISCRIM)" =~ $(REG_INTEL) ]] && echo true),true)
     COMPILER_FAMILY    :=Intel
     USER_DEFS          += -DLINUX_IFORT
  else ifeq ($(shell [[ "$(DISCRIM)" =~ $(REG_GNU) ]] && echo true),true)
     COMPILER_FAMILY    :=GNU
     USER_DEFS          += -DLINUX_GFORTRAN
  else
     $(error Could not determine compiler underlying mpifort/mpif90 )
  endif
endif

# %%%%% Test if Intel Fortran Compiler is selected %%%%%
REGEXP               :=(^[Ii][Ff][Oo][Rr][Tt])
ifeq ($(shell [[ "$(FC)" =~ $(REGEXP) ]] && echo true),true)

  # If we are building GCHP, then set the compile command to "mpifort",
  # which invokes the MPI magic.  Otherwise set it to $(FC). (bmy, 10/17/16)
  COMPILER_FAMILY    :=Intel
  USER_DEFS          += -DLINUX_IFORT
endif

# %%%%% Test if GNU Fortran Compiler is selected %%%%%
REGEXP               :=(^[Gg][Ff][Oo][Rr][Tt][Rr][Aa][Nn])
ifeq ($(shell [[ "$(FC)" =~ $(REGEXP) ]] && echo true),true)
  COMPILER_FAMILY    :=GNU
  USER_DEFS          += -DLINUX_GFORTRAN
endif

# %%%%% Test if PGI Fortran compiler is selected  %%%%%
REGEXP               :=(^[Pp][Gg])
ifeq ($(shell [[ "$(FC)" =~ $(REGEXP) ]] && echo true),true)
  COMPILER_FAMILY    :=PGI
  USER_DEFS          += -DLINUX_PGI
endif

# Is this GCHP?
ifeq ($(IS_HPC),1)
  COMPILE_CMD        :=mpifort
else
  COMPILE_CMD        :=$(FC)
endif

# %%%%% ERROR CHECK!  Make sure our compiler selection is valid! %%%%%
REGEXP               :=((-DLINUX_)?IFORT|PGI|GFORTRAN)
ifneq ($(shell [[ "$(USER_DEFS)" =~ $(REGEXP) ]] && echo true),true)
  $(error $(ERR_CMPLR))
endif

# Once we are sure the compiler is valid, then get the version number
COMPILER_VERSION_LONG :=$(shell $(FC) --version))
COMPILER_VERSION_LONG :=$(sort $(COMPILER_VERSION_LONG))

# For ifort, the 3rd substring of the sorted text is the version number.
# For pgfortran and gfortran, it's the 4th substring.
# NOTE: Future compiler updates may break this algorithm.
REGEXP      :=(^[Ii][Ff][Oo][Rr][Tt])
ifeq ($(shell [[ "$(FC)" =~ $(REGEXP) ]] && echo true),true)
 COMPILER_VERSION :=$(word 3, $(COMPILER_VERSION_LONG))
else
 COMPILER_VERSION :=$(word 4, $(COMPILER_VERSION_LONG))
endif

# Major version number of the compiler
# e.g. for gfortran 8.2.0, this would be "8"
COMPILER_MAJOR_VERSION   :=$(word 1,$(subst ., ,$(COMPILER_VERSION)))

#------------------------------------------------------------------------------
# Special flags for enabling experimental or development code
#------------------------------------------------------------------------------

# %%%%% DEVEL: Enable user-added experimental code %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(DEVEL)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DDEVEL
endif

# %%%%% DIAG_DEVEL: Enable experimental code specific to HEMCO %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(DIAG_DEVEL)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DDIAG_DEVEL
  BPCH_DIAG          :=no
endif

# %%%%% HCO_DEVEL: Enable experimental code specific to HEMCO %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(HCO_DEVEL)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DHCO_DEVEL
endif

# %%%%% HPC_DEVEL: Enable experimental code specific to GCHP %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(HPC_DEVEL)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DHPC_DEVEL
endif

# %%%%% Turn on tendencies computation  %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(USE_TEND)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DUSE_TEND
  BPCH_DIAG          :=no
endif

#------------------------------------------------------------------------------
# GEOS-Chem HP settings
#------------------------------------------------------------------------------

# %%%%% ESMF %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(ESMF)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DESMF_
  NO_GRID_NEEDED     :=1
endif

# %%%%% EXTERNAL_GRID %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(EXTERNAL_GRID)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DEXTERNAL_GRID
  NO_GRID_NEEDED     :=1
endif

# %%%%% EXTERNAL_FORCING %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(EXTERNAL_FORCING)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DEXTERNAL_FORCING
  NO_GRID_NEEDED     :=1
endif

#------------------------------------------------------------------------------
# Coupling GEOS-Chem to External Models settings
#------------------------------------------------------------------------------

# %%%%% NO_EXE %%%%%
# Setting NO_EXE=y will inhibit the creation of a final "geos" executable
# and create a "libGeosCore.a" in the lib folder instead.
# Used if you are linking GEOS-Chem routines to be driven by an external model.
# (hplin, 8/23/18)
EXE_NEEDED           :=1
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(NO_EXE)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DMODEL_
  EXE_NEEDED         :=0
endif

#------------------------------------------------------------------------------
# Diagnostic settings
#------------------------------------------------------------------------------

# %%%%% Use netCDF diagnostics if DEVEL=y %%%%%
ifdef DEVEL
  BPCH_DIAG          :=no
  BPCH_TBPC          :=no
endif

# %%%%% Turn on bpch code for TPCORE BC's if NEST is defined %%%%%
ifdef NEST
  BPCH_TPBC          :=yes
endif

# Turn on bpch diagnostics UNLESS specified otherwis
ifdef BPCH_DIAG
  BPCH_DIAG          :=yes
endif
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(BPCH_DIAG)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS        += -DBPCH_DIAG
endif

# %%%%% Turn netCDF diagnostics on by default to save out restart file %%%%%
REGEXP               :=(^[Nn]|^[Nn][Oo])
ifeq ($(shell [[ "$(NC_DIAG)" =~ $(REGEXP) ]] && echo true),true)

  # Set a flag to denote netCDF diagnostics are off
  IS_NC_DIAG         :=0

  # If netCDF diagnostics have not been explicitly specified, then activate
  # bpch diagnostics, bpch timeseries, AND bpch code for nested-grid BC's
  USER_DEFS          += -DBPCH_DIAG -DBPCH_TIMESER -DBPCH_TPBC

else

  # Set a flag to denote netCDF diagnostics are on
  IS_NC_DIAG         :=1

  # Turn on netCDF diagnostics if explicitly specified
  USER_DEFS          += -DNC_DIAG

  # If we are compiling GEOS-Chem "Classic", then also activate all bpch
  # timeseries diagnostics.  At this point (v11-02) there are some special
  # timeseries diagnostics that require local-time binning, which is not
  # yet available in the netCDF diagnostic output.  This will preserve
  # backwards compatibility for the time being. (bmy, 4/11/18)
  ifeq ($(IS_HPC),0)
     USER_DEFS       += -DBPCH_TIMESER
  endif

endif

#------------------------------------------------------------------------------
# KPP settings chemistry solver settings.  NOTE: We can't redefine CHEM 
# (since it is an environent variable), so define a shadow variable KPP_CHEM.
#------------------------------------------------------------------------------

# Test if the CHEM value is set
IS_CHEM_SET          :=0

# %%%%%  CHEM=Standard (aka benchmark) %%%%%
REGEXP               :=(^[Ss][Tt][Aa][Nn][Dd][Aa][Rr][Dd])
ifeq ($(shell [[ "$(CHEM)" =~ $(REGEXP) ]] && echo true),true)
  KPP_CHEM           :=Standard
  IS_CHEM_SET        :=1
  NO_REDUCED         :=yes
endif

# %%%%% Test if CHEM=SOA (same as Tropchem as of v11-02a) %%%%%
REGEXP               :=(^[Ss][Oo][Aa])
ifeq ($(shell [[ "$(CHEM)" =~ $(REGEXP) ]] && echo true),true)
  KPP_CHEM           :=Tropchem
  IS_CHEM_SET        :=1
endif

# %%%%% Test if CHEM=SOA_SVPOA %%%%%
REGEXP               :=(^[Ss][Oo][Aa]_[Ss][Vv][Pp][Oo][Aa])
ifeq ($(shell [[ "$(CHEM)" =~ $(REGEXP) ]] && echo true),true)
  KPP_CHEM           :=SOA_SVPOA
  IS_CHEM_SET        :=1
endif

# %%%%% Test if CHEM=Tropchem %%%%%
REGEXP               :=(^[Tt][Rr][Oo][Pp][Cc][Hh][Ee][Mm])
ifeq ($(shell [[ "$(CHEM)" =~ $(REGEXP) ]] && echo true),true)
  KPP_CHEM           :=Tropchem
  IS_CHEM_SET        :=1
endif

# %%%%% Test if CHEM=NOx_Ox_HC_Aer_Br (former name for Tropchem) %%%%%
REGEXP               :=(^[Nn][Oo][Xx]_[Oo][Xx]_[Hh][Cc]_[Aa][Ee][Rr]_[Bb][Rr])
ifeq ($(shell [[ "$(CHEM)" =~ $(REGEXP) ]] && echo true),true)
  KPP_CHEM           :=Tropchem
  IS_CHEM_SET        :=1
endif

# %%%%% Test if CHEM=Custom %%%%%
REGEXP               :=(^[Cc][Uu][Ss][Tt][Oo][Mm])
ifeq ($(shell [[ "$(CHEM)" =~ $(REGEXP) ]] && echo true),true)
  KPP_CHEM           :=Custom
  IS_CHEM_SET        :=1
endif

# %%%%%  Default setting %%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NOTE: For clarify in the future, the default setting should be to not set
# KPP_CHEM or IS_CHEM_SET if the CHEM compiler option is not passed. The default
# option would be reserved for specialty simulations that do not require the KPP
# code to be compiled. (mps, 4/22/16)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ifeq ($(IS_CHEM_SET),0)
  KPP_CHEM           :=Standard
  IS_CHEM_SET        :=1
endif

#------------------------------------------------------------------------------
# RRTMG radiative transfer model settings
#------------------------------------------------------------------------------

# %%%%% RRTMG %%%%%
RRTMG_NEEDED         :=0
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(RRTMG)" =~ $(REGEXP) ]] && echo true),true)
  RRTMG_NEEDED       :=1
  USER_DEFS          += -DRRTMG
endif

#------------------------------------------------------------------------------
# Met field settings
#------------------------------------------------------------------------------

# If the user has omitted MET, then throw an error UNLESS we are trying
# to compile with "clean", "distclean", "realclean", "doc", "help",
# "ncdfcheck", or "libnc".  These targets don't depend on the value of MET.
ifndef MET
  REGEXP :=($clean|^doc|^help|^libnc|^ncdfcheck|gigc_debug|the_nuclear_option|wipeout|wipeout.|baselib.)
  ifeq ($(shell [[ "$(MAKECMDGOALS)" =~ $(REGEXP) ]] && echo true),true)
    NO_MET_NEEDED    :=1
  else
    $(error $(ERR_MET))
  endif
endif

# We can skip the following checks for targets that don't require MET
ifndef NO_MET_NEEDED 

  # %%%%% MERRA-2 %%%%%
  REGEXP           :=(^[Mm][Ee][Rr][Rr][Aa]2|^[Mm][Ee][Rr][Rr][Aa].2)
  ifeq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS      += -DMERRA2
  endif

  # %%%%% GEOS-FP %%%%%
  REGEXP             :=(^[Gg][Ee][Oo][Ss][Ff][Pp])|(^[Gg][Ee][Oo][Ss].[Ff][Pp])
  ifeq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DGEOS_FP
  endif

  # %%%%% FLEXGRID %%%%%
  REGEXP             :=(^[Ff][Ll][Ee][Xx][Gg][Rr][Ii][Dd])
  ifeq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DFLEXGRID
  endif

  # %%%%% REDUCED VERTICAL GRID (default, unless specified otherwise) %%%%
  ifndef NO_REDUCED
    NO_REDUCED       :=no
  endif
  REGEXP              :=(^[Nn]|^[Nn][Oo])
  ifeq ($(shell [[ "$(NO_REDUCED)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DGRIDREDUCED
  endif

  # %%%%% ERROR CHECK!  Make sure our MET selection is valid! %%%%%
  REGEXP             :=(\-DGEOS_FP|\-DMERRA2)
  ifneq ($(shell [[ "$(USER_DEFS)" =~ $(REGEXP) ]] && echo true),true)
    $(error $(ERR_MET))
  endif

endif  # NO_MET_NEEDED

#------------------------------------------------------------------------------
# Horizontal grid settings
#------------------------------------------------------------------------------

# We can skip the following checks for targets that don't require GRID

# If the user has omitted GRID, then throw an error UNLESS we are trying
# to compile with "clean", "distclean", "realclean", "doc", "help",
# "ncdfcheck", or "libnc".  These targets don't depend on the value of GRID.
ifndef NO_GRID_NEEDED
  ifndef GRID
    REGEXP :=($clean|^doc|^help|^libnc|^ncdfcheck|gigc_debug|the_nuclear_option|wipeout|wipeout.|baselib.|^wipeout)
    ifeq ($(shell [[ $(MAKECMDGOALS) =~ $(REGEXP) ]] && echo true),true)
      NO_GRID_NEEDED :=1
    else
      $(error $(ERR_GRID))
    endif
  endif
endif

# We can skip the following checks for targets that don't require GRID
ifndef NO_GRID_NEEDED

  # %%%%% 4 x 5 %%%%%
  REGEXP             :=(^4.5|^4\.0.5\.0)
  ifeq ($(shell [[ "$(GRID)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DGRID4x5
  endif

  # %%%%% 2 x 2.5 %%%%%
  REGEXP             :=(^2.25|^2.2\.5|^2\.0.2\.5)
  ifeq ($(shell [[ "$(GRID)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DGRID2x25
  endif

  # %%%%% 0.5 x 0.625 %%%%%
  REGEXP             :=(^05.0625|^0\.5.0\.625)
  ifeq ($(shell [[ "$(GRID)" =~ $(REGEXP) ]] && echo true),true)

    # Ensure that MET=merra2
    REGEXP           :=(^[Mm][Ee][Rr][Rr][Aa]2)|(^[Mm][Ee][Rr][Rr][Aa].2)|(^[Ff][Ll][Ee][Xx][Gg][Rr][Ii][Dd])
    ifneq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
      $(error When GRID=05x0625, you can only use MET=merra2)
    endif

    # Ensure that a nested-grid option is selected
    # NOTE: For safety's sake: if a nested-grid option is selected then 
    # define the BPCH_TPBC cpp switch even if BPCH_TPBC=n was passed.
    ifndef NEST
      $(error $(ERR_NEST))
    else
      NEST_NEEDED    :=1
      USER_DEFS      += -DBPCH_TPBC -DGRID05x0625 
    endif
  endif

  # %%%%% 0.25 x 0.3125 %%%%%
  REGEXP             :=(^025.03125|^0\.25.0\.3125)
  ifeq ($(shell [[ "$(GRID)" =~ $(REGEXP) ]] && echo true),true)

    # Ensure that MET=geosfp
    REGEXP           :=(^[Gg][Ee][Oo][Ss][Ff][Pp])|(^[Gg][Ee][Oo][Ss].[Ff][Pp])|(^[Ff][Ll][Ee][Xx][Gg][Rr][Ii][Dd])
    ifneq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
      $(error When GRID=025x03125, you can only use MET=geosfp)
    endif

    # Ensure that a nested-grid option is selected
    # NOTE: For safety's sake: if a nested-grid option is selected then 
    # define the BPCH_TPBC cpp switch even if BPCH_TPBC=n was passed.
    ifndef NEST
      $(error $(ERR_NEST))
    else
      NEST_NEEDED    :=1
      USER_DEFS      += -DBPCH_TPBC -DGRID025x03125
    endif
  endif

  # %%%%% ERROR CHECK!  Make sure our GRID selection is valid! %%%%%
  REGEXP             := ((\-DGRID)?4x5|2x25|1x125|05x0666|05x0625|025x03125)
  ifneq ($(shell [[ "$(USER_DEFS)" =~ $(REGEXP) ]] && echo true),true)
    $(error $(ERR_GRID))
  endif

#------------------------------------------------------------------------------
# Nested grid settings
#------------------------------------------------------------------------------

  # %%%%% China (CH) %%%%%
  REGEXP             :=(^[Cc][Hh])
  ifeq ($(shell [[ "$(NEST)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DNESTED_CH
  endif

  # %%%%% Asia (AS) %%%%%
  REGEXP             :=(^[Aa][Ss])
  ifeq ($(shell [[ "$(NEST)" =~ $(REGEXP) ]] && echo true),true)
    # Ensure that GRID=05x0625
    REGEXP           :=(^05.0625|^0\.5.0\.625)
    ifneq ($(shell [[ "$(GRID)" =~ $(REGEXP) ]] && echo true),true)
      $(error NEST=as can only be used with GRID=05x0625)
    endif
    USER_DEFS        += -DNESTED_AS
  endif

  # %%%%% Europe (EU) %%%%%
  REGEXP             :=(^[Ee][Uu])
  ifeq ($(shell [[ "$(NEST)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DNESTED_EU
  endif

  # %%%%% North America (NA) %%%%%
  REGEXP             :=(^[Nn][Aa])
  ifeq ($(shell [[ "$(NEST)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DNESTED_NA
  endif

  # %%%%% Custom (CU) %%%%%
  REGEXP             :=(^[Cc][Uu])
  ifeq ($(shell [[ "$(NEST)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DNESTED_CU
  endif

  # %%%%% ERROR CHECK!  Make sure our NEST selection is valid! %%%%%
  ifdef NEST_NEEDED
    REGEXP           :=((\-DNESTED_)?AS|CH|EU|NA|CU)
    ifneq ($(shell [[ "$(USER_DEFS)" =~ $(REGEXP) ]] && echo true),true)
      $(error $(ERR_NEST))
    endif
  endif

endif  # NO_GRID_NEEDED

#------------------------------------------------------------------------------
# Coupled grid settings (yanyy,6/18/14)
#------------------------------------------------------------------------------

# %%%%% Couple %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(COUPLE)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DEXCHANGE
endif

# %%%%% China (CH) and 4x5 %%%%%
REGEXP               :=(^4.5[Cc][Hh]|^4\.0.5\.0[Cc][Hh])
ifeq ($(shell [[ "$(COUPLECH)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DEXCHANGE -DEXCHANGE_4x5_CH
endif

# %%%%% Europe (EU) and 4x5 %%%%%
REGEXP               :=(^4.5[Ee][Uu]|^4\.0.5\.0[Ee][Uu])
ifeq ($(shell [[ "$(COUPLEEU)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DEXCHANGE -DEXCHANGE_4x5_EU
endif

# %%%%% North America (NA) and 4x5 %%%%%
REGEXP               :=(^4.5[Nn][Aa]|^4\.0.5\.0[Nn][Aa])
ifeq ($(shell [[ "$(COUPLENA)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DEXCHANGE -DEXCHANGE_4x5_NA
endif

# %%%%% SE Asia (SE) and 4x5 %%%%%
REGEXP               :=(^4.5[Nn][Aa]|^4\.0.5\.0[Nn][Aa])
ifeq ($(shell [[ "$(COUPLESE)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DEXCHANGE -DEXCHANGE_4x5_SE
endif

# %%%%% China (CH) and 2x2.5 %%%%%
REGEXP               :=(^2.25[Cc][Hh]|^2.2\.5[Cc][Hh]|^2\.0.2\.5[Cc][Hh])
ifeq ($(shell [[ "$(COUPLECH)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DEXCHANGE -DEXCHANGE_2x25_CH
endif

# %%%%% Europe (EU) and 2x2.5 %%%%%
REGEXP               :=(^2.25[Ee][Uu]|^2.2\.5[Ee][Uu]|^2\.0.2\.5[Ee][Uu])
ifeq ($(shell [[ "$(COUPLEEU)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DEXCHANGE -DEXCHANGE_2x25_EU
endif

# %%%%% North America (NA) and 2x2.5 %%%%%
REGEXP               :=(^2.25[Nn][Aa]|^2.2\.5[Nn][Aa]|^2\.0.2\.5[Nn][Aa])
ifeq ($(shell [[ "$(COUPLENA)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DEXCHANGE -DEXCHANGE_2x25_NA
endif

# %%%%% SE Asia (SE) and 2x2.5 %%%%%
REGEXP               :=(^2.25[Nn][Aa]|^2.2\.5[Nn][Aa]|^2\.0.2\.5[Nn][Aa])
ifeq ($(shell [[ "$(COUPLESE)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DEXCHANGE -DEXCHANGE_2x25_SE
endif

# %%%%% ERROR CHECK!  Make sure our NEST selection is valid! %%%%%
ifdef COUPLE_NEEDED
  REGEXP             :=((\-DEXCHANGE_)?CH|NA|EU)
  ifneq ($(shell [[ "$(USER_DEFS)" =~ $(REGEXP) ]] && echo true),true)
    $(error $(ERR_COUPLE))
  endif
endif

#------------------------------------------------------------------------------
# Aerosol microphysics settings
# At present, TOMAS or APM cannot be compiled with NC_DIAG=y! (bmy, 8/7/18)
#------------------------------------------------------------------------------

# %%%%% TOMAS, 30 bins (default) %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(TOMAS)" =~ $(REGEXP) ]] && echo true),true)
  ifeq ($(IS_NC_DIAG),1) 
    $(error $(ERR_MICPHYS))
  else
    USER_DEFS        += -DTOMAS
  endif
endif

# %%%%% TOMAS, 40 bins %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(TOMAS40)" =~ $(REGEXP) ]] && echo true),true)
  ifeq ($(IS_NC_DIAG),1) 
    $(error $(ERR_MICPHYS))
  else
    USER_DEFS        += -DTOMAS -DTOMAS40
  endif
endif

# %%%%% TOMAS, 15 bins %%%%% 
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(TOMAS15)" =~ $(REGEXP) ]] && echo true),true)
  ifeq ($(IS_NC_DIAG),1) 
    $(error $(ERR_MICPHYS))
  else
    USER_DEFS        += -DTOMAS -DTOMAS15
  endif
endif

# %%%%% TOMAS, 12 bins %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(TOMAS12)" =~ $(REGEXP) ]] && echo true),true)
  ifeq ($(IS_NC_DIAG),1) 
    $(error $(ERR_MICPHYS))
  else
    USER_DEFS        += -DTOMAS -DTOMAS12
  endif
endif

# %%%%% APM %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(APM)" =~ $(REGEXP) ]] && echo true),true)
  ifeq ($(IS_NC_DIAG),1) 
    $(error $(ERR_MICPHYS))
  else
    USER_DEFS        += -DAPM
  endif
endif

#------------------------------------------------------------------------------
# Special chemistry settings
#------------------------------------------------------------------------------

# Activate Global Terrestrial Mercury Model (GTMM) if necessary
GTMM_NEEDED          :=0
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(GTMM_Hg)" =~ $(REGEXP) ]] && echo true),true)
  GTMM_NEEDED        :=1
  USER_DEFS          += -DGTMM_Hg
endif

# Option to turn off ISORROPIA for testing
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(NO_ISO)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DNO_ISORROPIA
endif

#------------------------------------------------------------------------------
# Performance profiling
#------------------------------------------------------------------------------

# Compile with TAU profiler (from ParaTools, Inc)
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(TAU_PROF)" =~ $(REGEXP) ]] && echo true),true)
  COMPILE_CMD        :=tau_f90.sh
endif

# Compile with GNU profiler (gprof)
IS_GPROF             :=0
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(GPROF)" =~ $(REGEXP) ]] && echo true),true)
  IS_GPROF           :=1
endif

#------------------------------------------------------------------------------
# Add test for mass conservation
#------------------------------------------------------------------------------
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(MASSCONS)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DMASSCONS
endif

###############################################################################
###                                                                         ###
###  Set linker commands for local and external libraries (incl. netCDF)    ###
###                                                                         ###
###############################################################################

# Test if we have found the nf-config file, which indicates a separate
# netCDF-Fortran build.  IS_NF_CONFIG=0 indicates that we found nf-config.
IS_NF_CONFIG         :=$(shell test -f $(GC_F_BIN)/nf-config; echo $$?)

# Test for GEOS-Chem-Libraries or onboard netCDF libraries
ifeq ($(shell [[ "$(GC_LIB)" =~ GEOS-Chem-Libraries ]] && echo true),true)

  #-----------------------------------------------------------------------
  # %%%%% We are using the GEOS-Chem-Libraries package %%%%%
  # 
  # Both netCDF-Fortran and netCDF-C library files are in the same path
  #-----------------------------------------------------------------------

  # NetCDF include command: 1 library path
  NC_INC_CMD         := -I$(GC_INCLUDE)

  # NetCDF link command: Add a workaround so that $GC_LIB will specify
  # the library path.  This should prevent any issues caused by building
  # the GEOS-Chem-Libraries in one location and moving them to another.
  NC_LINK_CMD        := $(shell $(GC_BIN)/nf-config --flibs)
  NC_LINK_CMD        += $(shell $(GC_BIN)/nc-config --libs)
  NC_LINK_CMD        := $(filter -l%,$(NC_LINK_CMD))
  NC_LINK_CMD        :=-L$(GC_LIB) $(NC_LINK_CMD)

else

  ifeq ($(IS_NF_CONFIG),0)

    #-----------------------------------------------------------------------
    # %%%%% We are using the onboard netCDF libraries %%%%%
    #
    # NetCDF-Fortran and NetCDF-C library files are in different paths,
    # which is typical of netCDF versions 4.2 and higher.
    #-----------------------------------------------------------------------

    # NetCDF include command: 2 library paths
    NC_INC_CMD       := -I$(GC_INCLUDE) -I$(GC_F_INCLUDE)

    # NetCDF link command: 2 sets of link commands
    NC_LINK_CMD      := $(shell $(GC_F_BIN)/nf-config --flibs)
    NC_LINK_CMD      += $(shell $(GC_BIN)/nc-config --libs)

  else

    #-----------------------------------------------------------------------
    # %%%%% We are using the onboard netCDF libraries %%%%%
    #
    # NetCDF-Fortran and NetCDF-C library files are in the same path,
    # which is typical netCDF versions earlier than 4.2.
    #-----------------------------------------------------------------------

    # NetCDF include command: 1 library path
    NC_INC_CMD       := -I$(GC_INCLUDE)

    # NetCDF link command: 1 set of link commands
    NC_LINK_CMD      := $(shell $(GC_BIN)/nc-config --flibs)

  endif

endif

# Save for backwards compatibility
NCL                  := $(NC_LINK_CMD)

#----------------------------
# For GEOS-Chem 
#----------------------------

# Base linker command: specify the library directory
LINK                 :=-L$(LIB)

# Append library for GTMM, if necessary
ifeq ($(GTMM_NEEDED),1)
  LINK               :=$(LINK) -lHg
endif

# Append library for RRTMG, if necessary
ifeq ($(RRTMG_NEEDED),1)
  LINK               :=$(LINK) -lrad
endif

# Create linker command to create the GEOS-Chem executable
LINK                 :=$(LINK) -lIsoropia -lHistory -lHCOI -lHCOX -lHCO 
LINK                 :=$(LINK) -lGeosUtil -lKpp -lHeaders -lNcUtils 
LINK                 :=$(LINK) $(NC_LINK_CMD)

#----------------------------
# For the HEMCO standalone
#----------------------------

# Create linker command to create the HEMCO standalone executable
LINK_HCO             :=-L$(LIB) -lHCOI -lHCOX -lHCO -lGeosUtil -lHeaders
LINK_HCO             :=$(LINK_HCO) -lNcUtils $(NC_LINK_CMD)

###############################################################################
###                                                                         ###
###  Test if the netCDF library was built with compression enabled          ###
###                                                                         ###
###  NOTE: Compressing the netCDF files will make it impossible to compare  ###
###  them for identical-ness in a unit test or diff test.  Therefore, we    ###
###  have added some extra checks to skip the compression if so desired.    ###
###                                                                         ###
###############################################################################

# Assume we will turn on netCDF compression (if present)
IS_DEFLATE           :=1

# Unless NC_NODEFLATE=y
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(NC_NODEFLATE)" =~ $(REGEXP) ]] && echo true),true)
  IS_DEFLATE         :=0
endif

# Or DEBUG=y.  This will make sure unit tests and diff tests aren't affected.
REGEXP               := (^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(DEBUG)" =~ $(REGEXP) ]] && echo true),true)
  IS_DEFLATE         :=0
endif

# Skip netCDF compression unless it's requested (or not a debug run)
ifeq ($(IS_DEFLATE),1)

  # Test if the "nf_def_var_deflate" function is defined in netcdf.inc
  # Look for netcdf.inc where the netCDF-Fortran library is located
  ifdef GC_F_INCLUDE
    GREP :=$(strip $(shell grep nf_def_var_deflate $(GC_F_INCLUDE)/netcdf.inc))
  else
    GREP :=$(strip $(shell grep nf_def_var_deflate $(GC_INCLUDE)/netcdf.inc))
  endif

  # Look for the second word of the combined search results
  WORD               :=$(word 2,"$(GREP)")

  # If it matches "nf_def_var_deflate", then define Cpp flag NC_HAS_COMPRESSION 
  ifeq ($(WORD),nf_def_var_deflate)
    USER_DEFS        += -DNC_HAS_COMPRESSION
  endif
endif

###############################################################################
###                                                                         ###
###  HPC Settings: Build & use ESMF & MAPL for Grid-Independent GEOS-Chem   ###
###                                                                         ###
###############################################################################

# If we are building w/ the HPC target, then include GIGC.mk as well
# Determine if we are building with the hpc target
ifeq ($(IS_HPC),1)
  ifneq ("$(wildcard $(CURDIR)/../GCHP/GIGC.mk)","")
    include $(CURDIR)/../GCHP/GIGC.mk
  else
  ifneq ("$(wildcard $(CURDIR)/../../GCHP/GIGC.mk)","")
    include $(CURDIR)/../../GCHP/GIGC.mk
  else
    $(error $(ERR_GIGC))
  endif
  endif
  #FFLAGS             += -double-size 32 -real-size 32 -r4
endif

###############################################################################
###                                                                         ###
###  Define settings for the GNU FORTRAN COMPILER (aka gfortran)            ###
###                                                                         ###
###############################################################################

ifeq ($(COMPILER_FAMILY),GNU) 

  # Get the GNU Fortran version
  GNU_VERSION        :=$(shell $(FC) -dumpversion)
  GNU_VERSION        :=$(subst .,,$(GNU_VERSION))
  NEWER_THAN_447     :=$(shell perl -e "print ($(GNU_VERSION) gt 447)")
  IS_GNU_8           :=$(shell perl -e "print ($(GNU_VERSION) ge 800)")

  # Base set of compiler flags
  FFLAGS             :=-cpp -w -std=legacy -fautomatic -fno-align-commons
  ifeq ($(IS_HPC),1)
    FFLAGS             += -fconvert=native
  else
    FFLAGS             += -fconvert=big-endian
  endif
  FFLAGS             += -fno-range-check

  # OPTIONAL: Add the GNU Fortran -march option, which compiles for a
  # specific computer architecture.  This may cause issues on some types
  # of CPUs (e.g. Intel), so we have left this as an optional argument.
  ifdef M_ARCH
    FFLAGS           += -march=$(M_ARCH)
  endif

  # Default optimization level for all routines (-O3)
  ifndef OPT
    # Options of interest
    #  -limf                Intel math libraries - machine must have them
    #  -O3                  Highest safe optimization level
    #  -march=native        Make the binary machine-specific. If in doubt, 
    #                        use a specific architecture, eg...
    #  -march=corei7-avx    Binary uses optimizations for 
    #                        Intel Sandy-Bridge Xeon (e.g. E5-2680)
    #  -mfpmath=sse         Use SSE extensions
    #  -funroll-loops       Enable loop unrolling
    #  -ffast-math          Enable fast math optimizations
    OPT              := -O3 -funroll-loops
    #OPT              := -O3 -march=corei7-avx -mfpmath=sse -funroll-loops
  endif

  # Pick compiler options for debug run or regular run 
  REGEXP             := (^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(DEBUG)" =~ $(REGEXP) ]] && echo true),true)
    #-fcheck=all would be more comprehensive but would force bounds checking
    FFLAGS           += -g -gdwarf-2 -gstrict-dwarf -O0
    FFLAGS           += -Wall -Wextra -Wconversion
    FFLAGS           += -Warray-temporaries -fcheck-array-temporaries
    TRACEBACK        :=yes
    USER_DEFS        += -DDEBUG
  else
    FFLAGS           += $(OPT)
  endif

  # Prevent any optimizations that would change numerical results
  #GFORTRAN_BAD#FFLAGS             += -fp-model source

  # Turn on OpenMP parallelization
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(OMP)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -fopenmp
  endif

  # Get Operating System (Linux = Linux; Darwin = MacOSX)
  ifndef UNAME
    UNAME            :=$(shell uname)
  endif

  # OSX compilation options
  ifeq ($(UNAME),Darwin)
    # This has not yet been tested
    $(error $(ERR_OSCOMP))
  #  FFLAGS           += -Wl,-stack_size,0x2cb410000  # 12 GB of stack space
  #  ifdef DEBUG
  #    FFLAGS         += -g0 -debug -save-temps -fpic -Wl,-no_pie
  #  endif
  endif

  # Add options for medium memory model.  This is to prevent G-C from 
  # running out of memory at hi-res, especially when using netCDF I/O.
  ifneq ($(UNAME),Darwin)
    #GFORTRAN_BAD#FFLAGS           += -mcmodel=medium -shared-intel
    FFLAGS           += -mcmodel=medium
  endif

  # Turn on checking for floating-point exceptions
  # These are approximately equivalent to -fpe0 -ftrapuv in IFORT
  # NOTE: GNU Fortran 4.4.7 does not allow for -finit-real-snan, so
  # we will only add this flag for versions newer than 4.4.7
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(FPE)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -ffpe-trap=invalid,zero,overflow
    ifeq ($(NEWER_THAN_447),1)
      FFLAGS           += -finit-real=snan
    endif
  endif
  ifeq ($(shell [[ "$(FPEX)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -ffpe-trap=invalid,zero,overflow
    ifeq ($(NEWER_THAN_447),1)
      FFLAGS           += -finit-real=snan
    endif
  endif

  # Add option for "array out of bounds" checking
  REGEXP             := (^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(BOUNDS)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -fbounds-check
  endif

  # Also add traceback option
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(TRACEBACK)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -fbacktrace
    ifndef DEBUG
       FFLAGS += -g
    endif
  endif

  # Compile for use with the GNU profiler (gprof), if necessary
  ifeq ($(IS_GPROF),1) 
    FFLAGS           += -pg
  endif

  # Add flexible precision declaration
  ifeq ($(PRECISION),8)
    USER_DEFS        += -DUSE_REAL8
  endif

  # Add timers declaration
  ifeq ($(TIMERS),1)
    USER_DEFS        += -DUSE_TIMERS
  endif

  # Append the user options in USER_DEFS to FFLAGS
  FFLAGS             += $(USER_DEFS)

  # Include options (i.e. for finding *.h, *.mod files)
  INCLUDE :=-J$(MOD) $(NC_INC_CMD)

  # Do not append the ESMF/MAPL/FVDYCORE includes for ISORROPIA, because it 
  # will not compile.  ISORROPIA is slated for removal shortly. (bmy, 11/21/14)
  INCLUDE_ISO        :=$(INCLUDE)

  # Append the ESMF/MAPL/FVDYCORE include commands
  ifeq ($(HPC),yes)
    INCLUDE          += $(MAPL_INC) $(ESMF_MOD) $(ESMF_INC) $(FV_INC)
  endif

  # Set the standard compiler variables
  CC                 :=
  F90                :=$(COMPILE_CMD) $(FFLAGS) $(INCLUDE)
  F90ISO             :=$(COMPILE_CMD) $(FFLAGS) $(INCLUDE_ISO)
  LD                 :=$(COMPILE_CMD) $(FFLAGS)
  FREEFORM           := -ffree-form -ffree-line-length-none
  R8                 := -fdefault-real-8 -fdefault-double-8

endif

###############################################################################
###                                                                         ###
###  Define settings for the INTEL FORTRAN COMPILER (aka ifort)             ###
###                                                                         ###
###############################################################################

ifeq ($(COMPILER_FAMILY),Intel) 

  # Base set of compiler flags
  FFLAGS             :=-cpp -w -auto -noalign -convert big_endian

  # Default optimization level for all routines (-O2)
  ifndef OPT
    OPT              := -O2
  endif

  # Pick compiler options for debug run or regular run 
  REGEXP             := (^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(DEBUG)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -g -O0 -check arg_temp_created -debug all
    TRACEBACK        :=yes
    USER_DEFS        += -DDEBUG
  else
    FFLAGS           += $(OPT) -vec-report0
  endif

  # Prevent any optimizations that would change numerical results
  FFLAGS             += -fp-model source

  # Turn on OpenMP parallelization
  # NOTE: ifort 18 and higher users -qopenmp instead of -openmp
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(OMP)" =~ $(REGEXP) ]] && echo true),true)
     REGEXP          :=^1[8-9]|^2|^3|^4|^5|^6|^7|^8|^9
     ifeq ($(shell [[ "$(COMPILER_MAJOR_VERSION)" =~ $(REGEXP) ]] && echo true),true)
      FFLAGS         += -qopenmp
    else
      FFLAGS         += -openmp
    endif
  endif

  # Get Operating System (Linux = Linux; Darwin = MacOSX)
  ifndef UNAME
    UNAME            :=$(shell uname)
  endif

  # OSX compilation options
  ifeq ($(UNAME),Darwin)
    FFLAGS           += -Wl,-stack_size,0x2cb410000  # 12 GB of stack space
    ifdef DEBUG
      FFLAGS         += -g0 -debug -save-temps -fpic -Wl,-no_pie
    endif
  endif

  # Add options for medium memory model.  This is to prevent G-C from 
  # running out of memory at hi-res, especially when using netCDF I/O.
  ifneq ($(UNAME),Darwin)
    FFLAGS           += -mcmodel=medium -shared-intel
  endif

  # Turn on checking for floating-point exceptions
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(FPE)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -fpe0 -ftrapuv
  endif
  ifeq ($(shell [[ "$(FPEX)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -fpe0 -ftrapuv
  endif

  # Add special IFORT optimization commands
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(IPO)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -ipo -static
  endif

  # Add option for "array out of bounds" checking
  REGEXP             := (^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(BOUNDS)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -check bounds
  endif

  # Also add traceback option
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(TRACEBACK)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -traceback
  endif

  # Compile for use with the GNU profiler (gprof), if necessary
  ifeq ($(IS_GPROF),1) 
    FFLAGS           += -p
  endif

  # Add flexible precision declaration
  ifeq ($(PRECISION),8)
    USER_DEFS        += -DUSE_REAL8
  endif

  # Add timers declaration
  ifeq ($(TIMERS),1)
    USER_DEFS        += -DUSE_TIMERS
  endif

  # Append the user options in USER_DEFS to FFLAGS
  FFLAGS             += $(USER_DEFS)

  # Include options (i.e. for finding *.h, *.mod files)
  INCLUDE            :=-module $(MOD) $(NC_INC_CMD)

  # Do not append the ESMF/MAPL/FVDYCORE includes for ISORROPIA, because it 
  # will not compile.  ISORROPIA is slated for removal shortly. (bmy, 11/21/14)
  INCLUDE_ISO        :=$(INCLUDE)

  # Append the ESMF/MAPL/FVDYCORE include commands
  ifeq ($(IS_HPC),1)
    INCLUDE          += $(MAPL_INC) $(ESMF_MOD) $(ESMF_INC) $(FV_INC)
  endif

  # Set the standard compiler variables
  CC                 :=
  F90                :=$(COMPILE_CMD) $(FFLAGS) $(INCLUDE)
  F90ISO             :=$(COMPILE_CMD) $(FFLAGS) $(INCLUDE_ISO)
  LD                 :=$(COMPILE_CMD) $(FFLAGS)
  FREEFORM           := -free
  #ifneq ($(shell [[ "$(HPC)" =~ $(REGEXP) ]] && echo true),true)
  #ifneq ($(HPC),yes)
    R8               := -r8
  #endif

endif

###############################################################################
###                                                                         ###
###  Define settings for the PORTLAND GROUP COMPILER (aka "pgfortran")      ###
###                                                                         ###
###############################################################################

ifeq ($(COMPILER_FAMILY),PGI) 

  # Base set of compiler flags
  FFLAGS             :=-Kieee -byteswapio -Mpreprocess -m64

  # Default optimization level for all routines (-fast)
  ifndef OPT
    OPT              :=-O2
   endif

  # Pick compiler options for debug run or regular run 
  REGEXP             := (^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(DEBUG)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -g -O0
    TRACEBACK        :=yes
    USER_DEFS        += -DDEBUG
  else
    FFLAGS           += $(OPT)
  endif

  # Add options for medium memory model.  This is to prevent G-C from 
  # running out of memory at hi-res, especially when using netCDF I/O.
  FFLAGS             += -mcmodel=medium

  # Turn on OpenMP parallelization
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(OMP)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -mp
  endif

  # Add option for suppressing PGI non-uniform memory access (numa) library 
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(NONUMA)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -mp=nonuma
  endif

  # Add option for "array out of bounds" checking
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(BOUNDS)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -Mbounds
  endif

  # Also add traceback option
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(TRACEBACK)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -traceback
  endif

  # Compile for use with the GNU profiler (gprof), if necessary
  ifeq ($(IS_GPROF),1) 
    FFLAGS           += -pg
  endif

  # Turn on checking for floating-point exceptions
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(FPE)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -Ktrap=fp
  endif

  # Switch to add detailed compiler prinotut
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(INFO)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -Minfo
  endif

  # Add flexible precision declaration
  ifeq ($(PRECISION),8)
    USER_DEFS        += -DUSE_REAL8
  endif

  # Add timers declaration
  ifeq ($(TIMERS),1)
    USER_DEFS        += -DUSE_TIMERS
  endif

  # Append the user options in USER_DEFS to FFLAGS
  FFLAGS             += $(USER_DEFS)

  # Include options (i.e. for finding *.h, *.mod files)
  INCLUDE            := -module $(MOD) $(NC_INC_CMD)

  # Do not append the ESMF/MAPL/FVDYCORE includes for ISORROPIA, because it 
  # will not compile.  ISORROPIA is slated for removal shortly. (bmy, 11/21/14)
  INCLUDE_ISO        :=$(INCLUDE)

  # Append the ESMF/MAPL/FVDYCORE include commands
  ifeq ($(IS_HPC),1)
   INCLUDE           += $(MAPL_INC) $(ESMF_MOD) $(ESMF_INC) $(FV_INC)
  endif

  # Set the standard compiler variables
  CC             :=gcc
  F90            :=$(COMPILE_CMD) $(FFLAGS) $(INCLUDE)
  F90ISO         :=$(COMPILE_CMD) $(FFLAGS) $(INCLUDE_ISO)
  LD             :=$(COMPILE_CMD) $(FFLAGS)
  FREEFORM       := -Mfree
  R8             := -Mextend -r8

endif

###############################################################################
###                                                                         ###
###  Specify pattern rules for compiliation                                 ###
###  (i.e. tell "make" how to compile files w/ different extensions)        ###
###                                                                         ###
###############################################################################

%.o : %.f
	$(F90) -c $<
%.o : %.F
	$(F90) -c $<
%.o : %.f90
	$(F90) -c $(FREEFORM) $<
%.o : %.F90
	$(F90) -c $(FREEFORM) $<
%.o : %.c
	$(CC) -c $*.c

###############################################################################
###                                                                         ###
###  Export global variables so that the main Makefile will see these       ###
###                                                                         ###
##############################6#################################################

export CC
export F90
export F90ISO
export FREEFORM
export LD
export LINK
export LINK_HCO
export R8
export SHELL
export NCL
export NC_LINK_CMD
export HPC
export PRECISION
export RRTMG_NEEDED
export RRTMG_CLEAN
export RRTMG_NO_CLEAN
export KPP_CHEM
export TIMERS
export IS_GNU_8

#EOC

###############################################################################
###                                                                         ###
###  Debug print output.  Normally you will leave the following lines       ###
###  commented out.  Uncomment these lines for testing.                     ###
###                                                                         ###
###############################################################################

#headerinfo:
#	@@echo '####### in Makefile_header.mk ########'
#	@@echo "COMPILER_FAMILY  : $(COMPILER_FAMILY)"
#	@@echo "COMPILER_VERSION : $(COMPILER_VERSION)" 
#	@@echo "COMPILE_CMD      : $(COMPILE_CMD)"
#	@@echo "DEBUG            : $(DEBUG)"
#	@@echo "BOUNDS           : $(BOUNDS)"
#	@@echo "F90              : $(F90)"
#	@@echo "CC               : $(CC)"
#	@@echo "INCLUDE          : $(INCLUDE)"
#	@@echo "LINK             : $(LINK)"
#	@@echo "USER_DEFS        : $(USER_DEFS)"
#	@@echo "IS_NC_CONFIG     : $(IS_NC_CONFIG)"
#	@@echo "NC_INC_CMD       : $(NC_INC_CMD)"
#	@@echo "NC_LINK_CMD      : $(NC_LINK_CMD)"
#	@@echo "NC_DIAG          : $(NC_DIAG)"
#	@@echo "BPCH_DIAG        : $(BPCH_DIAG)"
#	@@echo "NO_REDUCED       : $(NO_REDUCED)"

