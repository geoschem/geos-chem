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
#  See https://github.com/geoschem/geos-chem for complete history
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
ERR_CMPLR            :="Unknown Fortran compiler!  Must be one of ifort, gfortran, or mpifort|mpif90.  Check the FC environment variable in your .bashrc or .cshrc file."

# Error message for unknown compiler/OS combintation
ERR_OSCOMP           :="Makefile_header.mk not set up for this compiler/OS combination"

# Error message for bad two-way coupled model input (yanyy,6/18/14)
ERR_COUPLECH         :="Select a coupled grid for China/SE Asia: COUPLECH=2x25ch, COUPLECH=4x5ch"
ERR_COUPLENA         :="Select a coupled grid for North America: COUPLENA=2x25na, COUPLENA=4x5na"
ERR_COUPLEEU         :="Select a coupled grid for Europe       : COUPLEEU=2x25eu, COUPLEEU=4x5eu"
ERR_COUPLE           :="Select a coupled choice: COUPLE=yes"

# Error message for bad GCHP config
ERR_GCHP             :="Unable to find the GCHP configuration file GIGC.mk. Make sure you have cloned the GCHP repository into the GEOS-Chem repository as subdirectory GCHP."

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
# %%%%% If not HPC, then build as GEOS-Chem Classic
REGEXP               := (^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(HPC)" =~ $(REGEXP) ]] && echo true),true)
  IS_HPC             :=1
  OMP                :=no
# PRECISION          :=4
else
  IS_HPC             :=0
  USER_DEFS          += -DMODEL_CLASSIC
endif

# %%%%% Default to 8-byte precision unless specified otherwise %%%%%
ifndef PRECISION
 PRECISION           :=8
endif

# %%%%% Turn on traceback (error stack report) by default %%%%%
ifndef TRACEBACK
 TRACEBACK           :=yes
endif

# %%%%% Set default compiler %%%%%

# %%%%% Test if mpif90/mpifort is selected (for now assume ifort) %%%%%
REGEXP               :=(^[Mm][Pp][Ii])
ifeq ($(shell [[ "$(FC)" =~ $(REGEXP) ]] && echo true),true)
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

# Is this GCHP?
ifeq ($(IS_HPC),1)
  COMPILE_CMD        :=mpifort
else
  COMPILE_CMD        :=$(FC)
endif

# %%%%% ERROR CHECK!  Make sure our compiler selection is valid! %%%%%
REGEXP               :=((-DLINUX_)?IFORT|GFORTRAN)
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

# %%%%% Turn on Luo et al (2019) wetdep scheme %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(LUO_WETDEP)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DLUO_WETDEP
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

# %%%%% MODEL_GCHP %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(MODEL_GCHP)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DMODEL_GCHP
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

# Turn OFF bpch diagnostics UNLESS specified otherwise
ifdef BPCH_DIAG
  BPCH_DIAG          :=no
endif
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(BPCH_DIAG)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DBPCH_DIAG
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
# APM radiative transfer model settings
#------------------------------------------------------------------------------

# %%%%% RRTMG %%%%%
APM_NEEDED           :=0
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(APM)" =~ $(REGEXP) ]] && echo true),true)
  APM_NEEDED         :=1
  USER_DEFS          += -DAPM
endif

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
  USER_DEFS          += -DTOMAS -DTOMAS12
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

# Append library for RRTMG, if necessary
ifeq ($(APM_NEEDED),1)
  LINK               :=$(LINK) -lApm
endif

# Append library for GCHP, if necessary
ifeq ($(IS_HPC),1)
  LINK               :=$(LINK) -lGCHPint
endif

# Create linker command to create the GEOS-Chem executable
LINK                 :=$(LINK) -lIsorropia -lObsPack -lHistory
LINK                 :=$(LINK) -lHCOI -lHCOX -lHCO
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
    $(error $(ERR_GCHP))
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

###############################################################################
###                                                                         ###
###  Export global variables so that the main Makefile will see these       ###
###                                                                         ###
##############################6#################################################

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
export APM_NEEDED
export RRTMG_NEEDED
export RRTMG_CLEAN
export RRTMG_NO_CLEAN
export KPP_CHEM
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
#	@@echo "INCLUDE          : $(INCLUDE)"
#	@@echo "LINK             : $(LINK)"
#	@@echo "USER_DEFS        : $(USER_DEFS)"
#	@@echo "IS_NC_CONFIG     : $(IS_NC_CONFIG)"
#	@@echo "NC_INC_CMD       : $(NC_INC_CMD)"
#	@@echo "NC_LINK_CMD      : $(NC_LINK_CMD)"
#	@@echo "BPCH_DIAG        : $(BPCH_DIAG)"
