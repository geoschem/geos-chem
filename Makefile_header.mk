#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
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
#   make TARGET [ OPTIONAL-FLAGS ]
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
#                                                                             .
# FFLAGS is a local variable that is not returned to the "outside world", 
# but is only used locally.  COMPILER, HDF5, and OMP are all input via the
# command line or via environment variables.
#                                                                             .
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%       NOTE: The IBM/XLF compiler has not been validated yet.         %%%
# %%%                      Beta-testers welcome!                           %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%       NOTE: GEOS-Chem has not yet been ported to GNU Fortran.        %%%
# %%%                      Beta-testers welcome!                           %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#  07 Sep 2012 - R. Yantosca - Now add OPT variable to set global opt levels
#  07 Sep 2012 - R. Yantosca - Also set TRACEBACK for PGI compiler
#  25 Feb 2013 - S. Farina   - Add flag for TOMAS40
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Default settings for Makefile options
#==============================================================================

# IFORT is default compiler
ifndef COMPILER
COMPILER  := ifort
endif

# Get Operating System (Linux = Linux; Darwin = MacOSX)
ifndef UNAME
UNAME     := $(shell uname)
endif

# OpenMP is turned on by default
ifndef OMP
OMP       := yes
endif

# HDF5 I/O is turned off by default
ifndef HDF5
HDF5      := no
endif

# Use precise FP math optimization (i.e. to avoid numerical noise)
ifndef PRECISE
PRECISE   := yes
endif 

# TOMAS runs on single processor (at least for now!)
#ifeq ($(TOMAS),yes)
#OMP       := no
#endif

#==============================================================================
# Default values for variables
#==============================================================================

# If your system uses "/bin/sh", then uncomment this line!
SHELL     := /bin/sh

# If your system uses "/bin/bash", then uncomment this line!
#SHELL     := /bin/bash

# Library include path
NCI       := -I$(GC_INCLUDE)

# Library link path: first try to get the list of proper linking flags
# for this build of netCDF with nf-config and nc-config. 
NCL       := $(shell $(GC_BIN)/nf-config --flibs)
NCL       += $(shell $(GC_BIN)/nc-config --libs)
NCL       := $(filter -l%,$(NCL))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%% NOTE TO GEOS-CHEM USERS: If you do not have netCDF-4.2 installed
#%%%% Then you can add/modify the linking sequence here.  (This sequence
#%%%% is a guess, but is probably good enough for other netCDF builds.)
ifeq ($(NCL),) 
NCL       := -lnetcdf -lhdf5_hl -lhdf5 -lz
endif
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Prepend the library directory path to the linking sequence
NCL       := -L$(GC_LIB) $(NCL)

# Command to link to the various library files (-lHeaders should be last!)
LINK      := -L$(LIB) -lKpp -lIsoropia -lGeosUtil -lHeaders
LINK      := $(LINK) -lNcUtils $(NCL)

# Commands to link to libraries, for GTMM code (-lHeaders should be last!)
LHG       := -L$(LIB) -lKpp -lIsoropia -lHg -lGeosUtil -lHeaders
LHG       := $(LINK) -lNcUtils $(NCL)

# Add the HDF5 library link commands (optional)
ifeq ($(HDF5),yes) 
LINK      := $(LINK) -L$(H5L)
LHG       := $(LINK) -L$(H5L)
endif

# For ESMF development
ifeq ($(ESMF),yes)
LINK      += -lESMF $(LIB_CHEM_BASE) $(LIB_CHEM_SHARED) $(LIB_PILGRIM)   \
                    $(LIB_MAPL_BASE) $(LIB_CFIO) $(LIB_GFIO) $(LIB_MPEU) \
                    $(LIB_ESMF) $(LIB_SDF) $(LIB_SYS) $(LIB_MPI)         \
                    $(ESMF_LDFLAGS) -lmpi_cxx -lstdc++ -limf -lrt -ldl
LHG       += -lESMF
endif

#==============================================================================
# IFORT compilation options (default)
#==============================================================================
ifeq ($(COMPILER),ifort) 

# Default optimization level for all routines (-O2)
ifndef OPT
OPT       := -O2
endif

# Turn on -traceback option by default for debugging runs
ifdef DEBUG
TRACEBACK := yes
endif

# Pick compiler options for debug run or regular run 
ifdef DEBUG
FFLAGS    := -cpp -w -O0 -auto -noalign -convert big_endian -g -DDEBUG
else
FFLAGS    := -cpp -w $(OPT) -auto -noalign -convert big_endian -vec-report0 
endif

ifdef FPE
FFLAGS    += -debug parallel -fpe3 -ftrapuv
endif

# OSX compilation options
ifeq ($(UNAME),Darwin)
FFLAGS    += -Wl,-stack_size,0x2cb410000 # Allow 12GB of stack space
ifdef DEBUG
FFLAGS    += -g0 -debug -save-temps -fpic -Wl,-no_pie
endif
endif

# Add options for medium memory model.  This is to prevent G-C from 
# running out of memory at hi-res, especially when using netCDF I/O.
ifneq ($(UNAME),Darwin)
FFLAGS    += -mcmodel=medium -i-dynamic
endif

# Prevent any optimizations that would change numerical results
# This is needed to prevent numerical noise from ISORROPIA (bmy, 8/25/11)
ifeq ($(PRECISE),yes)
FFLAGS    += -fp-model source
endif

# Turn on OpenMP parallelization
ifeq ($(OMP),yes) 
FFLAGS    += -openmp -Dmultitask
endif

# Also add TOMAS aerosol microphysics option
ifeq ($(TOMAS),yes) 
FFLAGS    += -DTOMAS
endif

ifeq ($(TOMAS40),yes) 
FFLAGS    += -DTOMAS40
endif

# Also add APM aerosol microphysics option
ifeq ($(APM),yes) 
FFLAGS    += -DAPM
endif

# Add special IFORT optimization commands
ifdef IPO
FFLAGS    += -ipo -static
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS    += -CB
endif

# Also add traceback option
ifdef TRACEBACK
FFLAGS    += -traceback
endif

# Option to turn off ISORROPIA for testing
ifdef NO_ISO
FFLAGS    += -DNO_ISORROPIA
endif

# Include options (i.e. for finding *.h, *.mod files)
INCLUDE   := -I$(HDR) -module $(MOD) $(NCI)

# Also append HDF5 include commands (optional)
ifeq ($(HDF5),yes)
INCLUDE   += -DUSE_HDF5 -I$(H5I)
endif

# Also add ESMF linking option
ifeq ($(ESMF),yes)
FFLAGS    += -DESMF_
endif

ifeq ($(ESMF_TESTBED),yes)
FFLAGS    += -DESMF_TESTBED_
INCLUDE   += -I$(HDR)
endif

#-----------------------------------------------------------------------------
# Flags for interfacing GEOS-Chem with an external GCM (mlong, bmy, 9/6/12)
#
ifeq ($(DEVEL),yes)
FFLAGS    += -DDEVEL
endif

ifeq ($(EXTERNAL_GRID),yes)
FFLAGS    += -DEXTERNAL_GRID
endif

ifeq ($(EXTERNAL_FORCING),yes)
FFLAGS    += -DEXTERNAL_FORCING
endif
#----------------------------------------------------------------------------

CC        :=
F90       := ifort $(FFLAGS) $(INCLUDE)
LD        := ifort $(FFLAGS)
FREEFORM  := -free
R8        := -r8

endif

#==============================================================================
# Portland Group (PGI) compilation options
#==============================================================================
ifeq ($(COMPILER),pgi) 

# Default optimization level for all routines (-fast)
ifndef OPT
OPT       := -fast
endif

# Pick compiler options for debug run or regular run 
ifdef DEBUG 
FFLAGS    := -byteswapio -Mpreprocess -Bstatic -g -O0 
else
FFLAGS    := -byteswapio -Mpreprocess -Bstatic $(OPT)
endif

# Add options for medium memory model.  This is to prevent G-C from 
# running out of memory at hi-res, especially when using netCDF I/O.
FFLAGS    += -mcmodel=medium

# Turn on OpenMP parallelization
ifeq ($(OMP),yes) 
FFLAGS    += -mp -Mnosgimp -Dmultitask
endif

# Add option for suppressing PGI non-uniform memory access (numa) library 
ifeq ($(NONUMA),yes) 
FFLAGS    += -mp=nonuma
endif

# Also add TOMAS aerosol microphysics option
ifeq ($(TOMAS),yes) 
FFLAGS    += -DTOMAS
endif

ifeq ($(TOMAS40),yes) 
FFLAGS    += -DTOMAS40
endif

# Also add APM aerosol microphysics option
ifeq ($(APM),yes) 
FFLAGS    += -DAPM
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS    += -C
endif

# Also add traceback option
ifdef TRACEBACK
FFLAGS    += -traceback
endif

# Option to turn off ISORROPIA for testing
ifdef NO_ISO
FFLAGS    += -DNO_ISORROPIA
endif

# Include options (i.e. for finding *.h, *.mod files)
INCLUDE   := -I$(HDR) -module $(MOD) $(NCI)

# Also append HDF5 include commands (optional)
ifeq ($(HDF5),yes)
INCLUDE   += -DUSE_HDF5 -I$(H5I)
endif

#-----------------------------------------------------------------------------
# Flags for interfacing GEOS-Chem with an external GCM (mlong, bmy, 9/6/12)
#
ifeq ($(DEVEL),yes)
FFLAGS    += -DDEVEL
endif

ifeq ($(EXTERNAL_GRID),yes)
FFLAGS    += -DEXTERNAL_GRID
endif

ifeq ($(EXTERNAL_FORCING),yes)
FFLAGS    += -DEXTERNAL_FORCING
endif
#----------------------------------------------------------------------------

CC        := gcc
F90       := pgf90 $(FFLAGS) $(INCLUDE)
LD        := pgf90 $(FFLAGS)
FREEFORM  := -Mfree
R8        := -Mextend -r8

endif

#==============================================================================
# Specify pattern rules for compiliation 
# (i.e. tell "make" how to compile different types of source code files)
#==============================================================================
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

#==============================================================================
# Export global variables so that the main Makefile will see these
#==============================================================================
export CC
export F90
export FREEFORM
export LD
export LINK
export R8
export SHELL
export NCL

#EOC
#==============================================================================
# Print variables for testing/debugging purposes (uncomment if necessary)
##=============================================================================
#headerinfo:
#	@@echo '####### in Makefile_header.mk ########' 
#	@@echo "compiler: $(COMPILER)"
#	@@echo "debug   : $(DEBUG)"
#	@@echo "bounds  : $(BOUNDS)"
#	@@echo "f90     : $(F90)"
#	@@echo "cc      : $(CC)"
#	@@echo "include : $(INCLUDE)"
#	@@echo "link    : $(LINK)"
