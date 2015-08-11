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
#   # %%%%% GEOS-5 %%%%%
#   REGEXP    :=((^[Gg][Ee][Oo][Ss])?5|.5)
#   ifeq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
#   USER_DEFS += -DGEOS_5
#   endif
#                                                                             .
# The [[ ]] in bash is an evaluation.  The above ifeq statement uses regular
# expressions to test if the MET variable matches the string "GEOS" (case-
# insensitive) and either "5" or "any character and then a 5".  This will
# return true (via the "echo true" statement) for combinations like "GEOS-5", 
# "geos5", "Geos-5", "GeOs.5", etc.  This is a robust way of evaluating
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
ERR_CMPLR            :="Select a compiler: COMPILER=ifort, COMPILER=pgi"

# Error message for bad MET input
ERR_MET              :="Select a met field: MET=gcap, MET=geos4, MET=geos5, MET=merra, MET=geos-fp, MET=merra2)"

# Error message for bad GRID input
ERR_GRID             :="Select a horizontal grid: GRID=4x5. GRID=2x25, GRID=05x0666, GRID=05x0625, GRID=025x03125"

# Error message for bad NEST input
ERR_NEST             :="Select a nested grid: NEST=ch, NEST=eu, NEST=na NEST=se, NEST=cu"

# Error message for bad two-way coupled model input (yanyy,6/18/14)
ERR_COUPLECH         :="Select a coupled grid for Asia         : COUPLECH=2x25ch, COUPLECH=4x5ch"
ERR_COUPLENA         :="Select a coupled grid for North America: COUPLENA=2x25na, COUPLENA=4x5na"
ERR_COUPLEEU         :="Select a coupled grid for Europe       : COUPLEEU=2x25eu, COUPLEEU=4x5eu"
ERR_COUPLESE         :="Select a coupled grid for SE Asia      : COUPLESE=2x25se, COUPLEEU=4x5se"
ERR_COUPLE           :="Select a coupled choice: COUPLE=yes"

# Error message for bad GIGC config
ERR_GIGC             :="Unable to find the GIGC configuration file. Have you downloaded the GIGC?"

# Error message for bad GIGC config
ERR_GIGC             :="Unable to find the GIGC configuration file. Have you downloaded the GIGC?"

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

# %%%%% OpenMP parallelization default) %%%%%
ifndef OMP
  OMP                :=yes
endif

# %%%%% Set the HPC variable if we are building for use w/ ESMF/MPI %%%%
ifeq ($(shell [[ "$(MAKECMDGOALS)" =~ "hpc" ]] && echo true),true)
  HPC                :=yes
  export HPC
endif

# %%%%% For HPC, we disable OpenMP and turn on the full vertical grid %%%
ifeq ($(HPC),yes)
  OMP                :=no
  NO_REDUCED         :=yes
# PRECISION          :=4
endif

# %%%%% Default to 8-byte precision unless specified otherwise %%%%%
ifndef PRECISION
 PRECISION           :=8
endif

# %%%%% Default to Timers disabled %%%%%
ifndef TIMERS
 TIMERS              :=0
endif

# %%%%% Set default compiler %%%%%
ifndef COMPILER
  COMPILER           :=ifort
endif

# %%%%% Test if IFORT compiler is selected %%%%%
REGEXP               :=(^[Ii][Ff][Oo][Rr][Tt])
ifeq ($(shell [[ "$(COMPILER)" =~ $(REGEXP) ]] && echo true),true)
  COMPLER            :=ifort
  COMPILE_CMD        :=ifort
  USER_DEFS          += -DLINUX_IFORT
endif

# %%%%% Test if PGI compiler is selected  %%%%%
REGEXP               :=(^[Pp][Gg][Ii])
ifeq ($(shell [[ "$(COMPILER)" =~ $(REGEXP) ]] && echo true),true)
  COMPILER           :=pgi
  COMPILE_CMD        :=pgf90
  USER_DEFS          += -DLINUX_PGI
endif

# %%%%% ERROR CHECK!  Make sure our COMPILER selection is valid! %%%%%
REGEXP               :=((-DLINUX_)?IFORT|PGI)
ifneq ($(shell [[ "$(USER_DEFS)" =~ $(REGEXP) ]] && echo true),true)
  $(error $(ERR_CMPLR))
endif

#------------------------------------------------------------------------------
# GEOS-Chem HP settings
#------------------------------------------------------------------------------

# %%%%% DEVEL %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(DEVEL)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DDEVEL
endif

# %%%%% ESMF %%%%%
REGEXP               := (^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(ESMF)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DESMF_
  NO_GRID_NEEDED     :=1
endif

# %%%%% EXTERNAL_GRID %%%%%
REGEXP               := (^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(EXTERNAL_GRID)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DEXTERNAL_GRID
  NO_GRID_NEEDED     :=1
endif

# %%%%% EXTERNAL_FORCING %%%%%
REGEXP               := (^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(EXTERNAL_FORCING)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DEXTERNAL_FORCING
  NO_GRID_NEEDED     :=1
endif

# %%%%% NO_BPCH (for disabling old diagnostic arrays) %%%%%
REGEXP               := (^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(NO_BPCH)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DNO_BPCH
endif

#------------------------------------------------------------------------------
# KPP settings chemistry solver settings.  NOTE: We can't redefine CHEM 
# (since it is an environent variable), so define a shadow variable KPP_CHEM.
# %%%%% NOTE: These will probably be obsolete when FLEXCHEM is added. %%%%%
#------------------------------------------------------------------------------

# Test if the CHEM value is set
IS_CHEM_SET          :=0

# %%%%% Test if CHEM=UCX (will also turn on UCX) %%%%%
REGEXP               :=(^[Uu][Cc][Xx])
ifeq ($(shell [[ "$(CHEM)" =~ $(REGEXP) ]] && echo true),true)
  UCX                :=y
  KPP_CHEM           :=UCX
  IS_CHEM_SET        :=1
endif

# %%%%% Test if CHEM=SOA %%%%%
REGEXP               :=(^[Ss][Oo][Aa])
ifeq ($(shell [[ "$(CHEM)" =~ $(REGEXP) ]] && echo true),true)
  KPP_CHEM           :=SOA
  IS_CHEM_SET        :=1
endif

# %%%%% Test if CHEM=SOA_SVPOA %%%%%
REGEXP               :=(^[Ss][Oo][Aa]_[Ss][Vv][Pp][Oo][Aa])
ifeq ($(shell [[ "$(CHEM)" =~ $(REGEXP) ]] && echo true),true)
  KPP_CHEM           :=SOA_SVPOA
  IS_CHEM_SET        :=1
endif

# %%%%% Test if CHEM=NOx_Ox_HC_Aer_Br %%%%%
REGEXP               :=(^[Nn][Oo][Xx]_[Oo][Xx]_[Hh][Cc]_[Aa][Ee][Rr]_[Bb][Rr])
ifeq ($(shell [[ "$(CHEM)" =~ $(REGEXP) ]] && echo true),true)
  KPP_CHEM           :=NOx_Ox_HC_Aer_Br
  IS_CHEM_SET        :=1
endif

# %%%%% Test if CHEM=tropchem (synonym for NOx_Ox_HC_Aer_Br) %%%%%
REGEXP               :=(^[Tt][Rr][Oo][Pp][Cc][Hh][Ee][Mm])
ifeq ($(shell [[ "$(CHEM)" =~ $(REGEXP) ]] && echo true),true)
  KPP_CHEM           :=NOx_Ox_HC_Aer_Br
  IS_CHEM_SET        :=1
endif

# %%%%%  Default setting: CHEM=benchmark %%%%%
ifeq ($(IS_CHEM_SET),0)
  KPP_CHEM           :=benchmark
  IS_CHEM_SET        :=1
endif

#------------------------------------------------------------------------------
# UCX stratospheric-tropospheric chemistry settings
#------------------------------------------------------------------------------

# %%%%% UCX %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(UCX)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DUCX
  NO_REDUCED         :=yes
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

# %%%%% Give users the option to make realclean except for RRTMG   %%%%%
# %%%%% if they set variables  RRTMG_NOCLEAN=y or RRTMG_NO_CLEAN=y %%%%%
# %%%%% This should help reduce the amount of time to recompile.   %%%%%
RRTMG_CLEAN          :=1
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(RRTMG_NO_CLEAN)" =~ $(REGEXP) ]] && echo true),true)
  RRTMG_CLEAN        :=0
endif
ifeq ($(shell [[ "$(RRTMG_NOCLEAN)" =~ $(REGEXP) ]] && echo true),true)
  RRTMG_CLEAN        :=0
endif


#------------------------------------------------------------------------------
# Met field settings
#------------------------------------------------------------------------------

# If the user has omitted MET, then throw an error UNLESS we are trying
# to compile with "clean", "distclean", "realclean", "doc", "help",
# "ncdfcheck", or "libnc".  These targets don't depend on the value of MET.
ifndef MET
  REGEXP :=($clean|^doc|^help|^libnc|^ncdfcheck|gigc_debug|the_nuclear_option|wipeout.|baselib.)
  ifeq ($(shell [[ "$(MAKECMDGOALS)" =~ $(REGEXP) ]] && echo true),true)
    NO_MET_NEEDED    :=1
  else
    $(error $(ERR_MET))
  endif
endif

# We can skip the following checks for targets that don't require MET
ifndef NO_MET_NEEDED 

  # %%%%% GCAP %%%%%
  REGEXP             :=(^[Gg][Cc][Aa][Pp])
  ifeq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DGCAP
  endif

  # %%%%% GEOS-4 %%%%%
  REGEXP             :=((^[Gg][Ee][Oo][Ss])?4|.4)
  ifeq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DGEOS_4
  endif

  # %%%%% GEOS-5 %%%%%
  REGEXP             :=((^[Gg][Ee][Oo][Ss])?5|.5)
  ifeq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DGEOS_5
  endif

  # %%%%% MERRA or MERRA2 %%%%%
  # We have to employ a double regexp test in order to prevent
  # inadvertently setting MERRA if we want MERRA2. (bmy, 8/11/15)
  REGEXP             :=(^[Mm][Ee][Rr][Rr][Aa])
  ifeq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
    REGEXP           :=(^[Mm][Ee][Rr][Rr][Aa]2|^[Mm][Ee][Rr][Rr][Aa].2)
    ifeq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
      USER_DEFS      += -DMERRA2
    else
      USER_DEFS      += -DMERRA
    endif
  endif

  # %%%%% GEOS-FP %%%%%
  REGEXP             :=(^[Gg][Ee][Oo][Ss][Ff][Pp])|(^[Gg][Ee][Oo][Ss].[Ff][Pp])
  ifeq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DGEOS_FP
  endif

  # %%%%% REDUCED VERTICAL GRID (default, unless specified otherwise) %%%%
  ifndef NO_REDUCED
    USER_DEFS        += -DGRIDREDUCED
  else
    REGEXP           :=(^[Yy]|^[Yy][Ee][Ss])
    ifeq ($(shell [[ "$(NO_REDUCED)" =~ $(REGEXP) ]] && echo true),true)
    endif
  endif

  # %%%%% ERROR CHECK!  Make sure our MET selection is valid! %%%%%
  REGEXP             :=(\-DGCAP|\-DGEOS_4|\-DGEOS_5|\-DMERRA|\-DGEOS_FP\-DMERRA2)
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
    REGEXP :=($clean|^doc|^help|^libnc|^ncdfcheck|gigc_debug|the_nuclear_option|wipeout.|baselib.|^wipeout)
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

  # %%%%% 1 x 1.25 %%%%%
  REGEXP             :=(^1.125|^1.1\.25|^1\.0.1\.25)
  ifeq ($(shell [[ "$(GRID)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DGRID1x125
  endif

  # %%%%% 0.5 x 0.666 %%%%%
  REGEXP             :=(^05.066.|^0\.5.0\.066.)
  ifeq ($(shell [[ "$(GRID)" =~ $(REGEXP) ]] && echo true),true)

    # Ensure that MET=geos5
    REGEXP           := ((^[Gg][Ee][Oo][Ss])?5|.5)
    ifneq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
      $(error When GRID=05x0666, you can only use MET=geos5)
    endif

    # Ensure that a nested-grid option is selected
    ifndef NEST
      $(error $(ERR_NEST))
    else
      NEST_NEEDED    :=1
      USER_DEFS      += -DGRID05x0666
    endif
  endif

  # %%%%% 0.5 x 0.625 %%%%%
  REGEXP             :=(^05.0625|^0\.5.0\.625)
  ifeq ($(shell [[ "$(GRID)" =~ $(REGEXP) ]] && echo true),true)

    # Ensure that MET=geos-fp
    REGEXP           :=(^[Mm][Ee][Rr][Rr][Aa]2)|(^[Mm][Ee][Rr][Rr][Aa].2)
    ifneq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
      $(error When GRID=05x0625, you can only use MET=merra2)
    endif

    # Ensure that a nested-grid option is selected
    ifndef NEST
      $(error $(ERR_NEST))
    else
      NEST_NEEDED    :=1
      USER_DEFS      += -DGRID05x0625
    endif
  endif

  # %%%%% 0.25 x 0.3125 %%%%%
  REGEXP             :=(^025.03125|^0\.25.0\.3125)
  ifeq ($(shell [[ "$(GRID)" =~ $(REGEXP) ]] && echo true),true)

    # Ensure that MET=geos-fp
    REGEXP           :=(^[Gg][Ee][Oo][Ss][Ff][Pp])|(^[Gg][Ee][Oo][Ss].[Ff][Pp])
    ifneq ($(shell [[ "$(MET)" =~ $(REGEXP) ]] && echo true),true)
      $(error When GRID=025x03125, you can only use MET=geos-fp)
    endif

    # Ensure that a nested-grid option is selected
    ifndef NEST
      $(error $(ERR_NEST))
    else
      NEST_NEEDED    :=1
      USER_DEFS      += -DGRID025x03125
    endif
  endif

  # %%%%% ERROR CHECK!  Make sure our GRID selection is valid! %%%%%
  REGEXP             := ((\-DGRID)?4x5|2x25|1x125|05x0666|025x03125)
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

  # %%%%% SE Asia (SE) %%%%%
  REGEXP             :=(^[Ss][Ee])
  ifeq ($(shell [[ "$(NEST)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DNESTED_SE
  endif

  # %%%%% Custom (CU) %%%%%
  REGEXP             :=(^[Cc][Uu])
  ifeq ($(shell [[ "$(NEST)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DNESTED_CU
  endif

  # %%%%% ERROR CHECK!  Make sure our NEST selection is valid! %%%%%
  ifdef NEST_NEEDED
    REGEXP           :=((\-DNESTED_)?CH|NA|EU|SE|CU)
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
#------------------------------------------------------------------------------

# %%%%% TOMAS, 30 bins (default) %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(TOMAS)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DTOMAS
endif

# %%%%% TOMAS, 40 bins %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(TOMAS40)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DTOMAS -DTOMAS40
endif

# %%%%% TOMAS, 15 bins %%%%% 
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(TOMAS15)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DTOMAS -DTOMAS15
endif

# %%%%% TOMAS, 12 bins %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(TOMAS12)" =~ $(REGEXP) ]] && echo true),true)
USER_DEFS            += -DTOMAS -DTOMAS12
endif

# %%%%% APM %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(APM)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DAPM
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

# Specify year of tagged O3 prod/loss data
# NOTE: THIS IS OBSOLETE W/ HEMCO! (bmy, 11/21/14)
ifdef TAGO3YR
  USER_DEFS          += -DUSE_THIS_O3_YEAR=$(TAGO3YR)
endif

#------------------------------------------------------------------------------
# TAU Performance Profiling (only works w/ IFORT for now)
#------------------------------------------------------------------------------
ifeq ($(COMPILER),ifort)
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(TAU_PROF)" =~ $(REGEXP) ]] && echo true),true)
    COMPILE_CMD      :=tau_f90.sh
  endif
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

# netCDF Library include path.  
# Test if a separate netcdf-fortran library is specified.
ifdef GC_F_INCLUDE
  NC_INC_CMD         := -I$(GC_INCLUDE) -I$(GC_F_INCLUDE)
else
  NC_INC_CMD         := -I$(GC_INCLUDE)
endif

# Get the version number (e.g. "4130"=netCDF 4.1.3; "4200"=netCDF 4.2, etc.)
NC_VERSION           :=$(shell $(GC_BIN)/nc-config --version)
NC_VERSION           :=$(shell echo "$(NC_VERSION)" | sed 's|netCDF ||g')
NC_VERSION           :=$(shell echo "$(NC_VERSION)" | sed 's|\.||g')
NC_VERSION_LEN       :=$(shell perl -e "print length $(NC_VERSION)")
ifeq ($(NC_VERSION_LEN),3)
 NC_VERSION          :=$(NC_VERSION)0
endif
ifeq ($(NC_VERSION_LEN),2) 
 NC_VERSION          :=$(NC_VERSION)00
endif

# Test if we have at least netCDF 4.2.0.0
AT_LEAST_NC_4200     :=$(shell perl -e "print ($(NC_VERSION) ge 4200)")

ifeq ($(AT_LEAST_NC_4200),1) 

  #-------------------------------------------------------------------------
  # netCDF 4.2 and higher:
  # Use "nf-config --flibs" and "nc-config --libs"
  # Test if a separate netcdf-fortran path is specified
  #-------------------------------------------------------------------------
  ifdef GC_F_BIN 
    NC_LINK_CMD      := $(shell $(GC_F_BIN)/nf-config --flibs)
  else
    NC_LINK_CMD      := $(shell $(GC_BIN)/nf-config --flibs)
  endif
  NC_LINK_CMD        += $(shell $(GC_BIN)/nc-config --libs)

else

  #-----------------------------------------------------------------------
  # Prior to netCDF 4.2:
  # Use "nc-config --flibs" and nc-config --libs
  #-----------------------------------------------------------------------
  NC_LINK_CMD        := $(shell $(GC_BIN)/nc-config --flibs)
  NC_LINK_CMD        += $(shell $(GC_BIN)/nc-config --libs)

endif

#=============================================================================
#%%%%% FIX FOR USE WITH THE GEOS-Chem-Libraries (bmy, 1/13/15)
#%%%%% 
#%%%%% If your GEOS-Chem-Libraries netCDF/HDF5 package was built in one 
#%%%%% directory and then moved somewhere else, then nf-config and nc-config 
#%%%%% may not return the proper link directory path.  
#%%%%% 
#%%%%% To avoid this error, we shall test if the $GC_LIB environment variable 
#%%%%% contains the text "GEOS-Chem-Libraries".  (Recall that $GC_LIB is 
#%%%%% defined in either your .bashrc or .cshrc file depending on which Unix 
#%%%%% shell you use.)  If we find the text "GEOS-Chem-Libraries" in $GC_LIB, 
#%%%%% then we shall override the library path returned by nf-config and 
#%%%%% nc-config with the path specified by $GC_LIB.  This will ensure that 
#%%%%% we point to the location where the GEOS-Chem-Libraries are installed.
#%%%%%
#%%%%% NOTE: This fix should work for most users.  If it does not work, then
#%%%%% contact the GEOS-Chem Support Team (geos-chem-support@as.harvard.edu).
#%%%%%
REGEXP               :="GEOS-Chem-Libraries"
ifeq ($(shell [[ "$(GC_LIB)" =~ $(REGEXP) ]] && echo true),true)
  NC_LINK_CMD        := $(filter -l%,$(NC_LINK_CMD))
  NC_LINK_CMD        :=-L$(GC_LIB) $(NC_LINK_CMD)
endif
#=============================================================================

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
LINK                 :=$(LINK) -lIsoropia -lHCOI -lHCOX -lHCO -lGeosUtil -lKpp
LINK                 :=$(LINK) -lHeaders -lNcUtils $(NC_LINK_CMD)

#----------------------------
# For the HEMCO standalone
#----------------------------

# Create linker command to create the HEMCO standalone executable
LINK_HCO             :=-L$(LIB) -lHCOI -lHCOX -lHCO -lGeosUtil -lHeaders
LINK_HCO             :=$(LINK_HCO) -lNcUtils $(NC_LINK_CMD)

###############################################################################
###                                                                         ###
###  HPC Settings: Build & use ESMF & MAPL for Grid-Independent GEOS-Chem   ###
###                                                                         ###
###############################################################################

# If we are building w/ the HPC target, then include GIGC.mk as well
# Determine if we are building with the hpc target
ifeq ($(HPC),yes)
  ifneq ("$(wildcard $(CURDIR)/../GIGC/GIGC.mk)","")
    include $(CURDIR)/../GIGC/GIGC.mk
  else
  ifneq ("$(wildcard $(CURDIR)/../../GIGC/GIGC.mk)","")
    include $(CURDIR)/../../GIGC/GIGC.mk
  else
    $(error $(ERR_GIGC))
  endif
  endif
  #FFLAGS             += -double-size 32 -real-size 32 -r4
endif

###############################################################################
###                                                                         ###
###  IFORT compilation options.  This is the default compiler.              ###
###                                                                         ###
###############################################################################

ifeq ($(COMPILER),ifort) 

  # Default optimization level for all routines (-O2)
  ifndef OPT
    OPT              := -O2
  endif

  # Pick compiler options for debug run or regular run 
  REGEXP             := (^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(DEBUG)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           :=-cpp -w -O0 -auto -noalign -convert big_endian
    FFLAGS           += -g -check arg_temp_created -debug all
    TRACEBACK        := yes
    USER_DEFS        += -DDEBUG
  else
    FFLAGS           :=-cpp -w $(OPT) -auto -noalign -convert big_endian
    FFLAGS           += -vec-report0
  endif

  # Prevent any optimizations that would change numerical results
  FFLAGS             += -fp-model source

  # Turn on OpenMP parallelization
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(OMP)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -openmp
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

  # Loosen KPP tolerances upon non-convergence and try again
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(KPP_SOLVE_ALWAYS)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DKPP_SOLVE_ALWAYS
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
  ifeq ($(HPC),yes)
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
###  Portland Group (PGF90) compilation options                             ###
###                                                                         ###
###############################################################################

ifeq ($(COMPILER),pgi) 

  # Default optimization level for all routines (-fast)
  ifndef OPT
    OPT              :=-fast
   endif

  # Pick compiler options for debug run or regular run 
  REGEXP             := (^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(DEBUG)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           :=-byteswapio -Mpreprocess -g -O0 
    USER_DEFS        += -DDEBUG
  else
    FFLAGS           :=-byteswapio -Mpreprocess $(OPT)
  endif

  # Add options for medium memory model.  This is to prevent G-C from 
  # running out of memory at hi-res, especially when using netCDF I/O.
  FFLAGS             += -mcmodel=medium

  # Turn on OpenMP parallelization
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(OMP)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -mp -Mnosgimp -Dmultitask
  endif

  # Add option for suppressing PGI non-uniform memory access (numa) library 
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(NONUMA)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -mp=nonuma
  endif

  # Add option for "array out of bounds" checking
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(BOUNDS)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -C
  endif

  # Also add traceback option
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(TRACEBACK)" =~ $(REGEXP) ]] && echo true),true)
    FFLAGS           += -traceback
  endif

  # Loosen KPP tolerances upon non-convergence and try again
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(KPP_SOLVE_ALWAYS)" =~ $(REGEXP) ]] && echo true),true)
    USER_DEFS        += -DKPP_SOLVE_ALWAYS
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
  ifeq ($(HPC),yes)
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
###############################################################################
 
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

#EOC

###############################################################################
###                                                                         ###
###  Debug print output.  Normally you will leave the following lines       ###
###  commented out.  Uncomment these lines for testing.                     ###
###                                                                         ###
###############################################################################

headerinfo:
	@@echo '####### in Makefile_header.mk ########' 
	@@echo "COMPILER    : $(COMPILER)"
	@@echo "DEBUG       : $(DEBUG)"
	@@echo "BOUNDS      : $(BOUNDS)"
	@@echo "F90         : $(F90)"
	@@echo "CC          : $(CC)"
	@@echo "INCLUDE     : $(INCLUDE)"
	@@echo "LINK        : $(LINK)"
	@@echo "USERDEFS    : $(USER_DEFS)"
	@@echo "NC_INC_CMD  : $(NC_INC_CMD)"
	@@echo "NC_LINK_CMD : $(NC_LINK_CMD)"
