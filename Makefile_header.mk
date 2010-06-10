# $Id: Makefile_header.mk,v 1.13 2010/03/15 19:33:25 ccarouge Exp $
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
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Default settings for Makefile options
#==============================================================================

# IFORT is default compiler
ifndef COMPILER
COMPILER = ifort
endif

# OpenMP is turned on by default
ifndef OMP
OMP = yes
endif

# HDF5 output is turned off by defautl
ifndef HDF5
HDF5 = no
endif

# TOMAS runs on single processor (at least for now!)
ifeq ($(TOMAS),yes)
OMP = no
endif

#==============================================================================
# Default values for variables
#==============================================================================

# If your system uses "/bin/sh", then uncomment this line!
SHELL = /bin/sh

# If your system uses "/bin/bash", then uncomment this line!
#SHELL = /bin/bash

# If you have HDF5 installed on your system, then define both the include
# (H5I) and library paths (H5L) here!  Otherwise leave these blank.
H5I = /home/bmy/NASA/basedir/x86_64-unknown-linux-gnu/ifort/Linux/include/hdf5
H5L = /home/bmy/NASA/basedir/x86_64-unknown-linux-gnu/ifort/Linux/lib

# Link to library files created from code in the various subdirs
# NOTE: -lGeosUtil should always be last!
LINK  = -L$(LIB) -lKpp -lIsoropia -lGeosUtil
LHG   = -L$(LIB) -lKpp -lIsoropia -lHg -lGeosUtil

# Add the HDF5 library link commands if necessary
ifeq ($(HDF5),yes) 
LINK += -L$(H5L) -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran -lhdf5 -lsz -lz -lm
LHG  += -L$(H5L) -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran -lhdf5 -lsz -lz -lm
endif


#==============================================================================
# IFORT compilation options (default)
#==============================================================================
ifeq ($(COMPILER),ifort) 

# Turn on -traceback option by default for debugging runs
ifdef DEBUG
TRACEBACK=yes
endif

# Pick compiler options for debug run or regular run 
ifdef DEBUG
FFLAGS   = -cpp -w -O0 -auto -noalign -convert big_endian -g -p
else
FFLAGS   = -cpp -w -O2 -auto -noalign -convert big_endian -vec-report0
endif

# Turn on OpenMP parallelization
ifeq ($(OMP),yes) 
FFLAGS  += -openmp -Dmultitask
endif

# Also add TOMAS aerosol microphysics option
ifeq ($(TOMAS),yes) 
FFLAGS  += -DTOMAS
endif

# Add special IFORT optimization commands
ifdef IPO
FFLAGS  += -ipo -static
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS  += -CB
endif

# Also add traceback option
ifdef TRACEBACK
FFLAGS  += -traceback
endif

# Include options (i.e. for finding *.h, *.mod files)
INCLUDE  = -I$(HDR) -module $(MOD)

# Also append HDF5 include commands if necessary
ifeq ($(HDF5),yes)
INCLUDE += -DUSE_HDF5 -I$(H5I)
endif

CC       =
F90      = ifort $(FFLAGS) $(INCLUDE)
LD       = ifort $(FFLAGS)
FREEFORM = -free
R8       = -r8

endif

#==============================================================================
# Portland Group (PGI) compilation options
#==============================================================================
ifeq ($(COMPILER),pgi) 

# Pick compiler options for debug run or regular run 
ifdef DEBUG 
FFLAGS   = -byteswapio -Mpreprocess -Bstatic -g -O0 
else
FFLAGS   = -byteswapio -Mpreprocess -Bstatic -fast 
endif

# Turn on OpenMP parallelization
ifeq ($(OMP),yes) 
FFLAGS  += -mp -Mnosgimp -Dmultitask
endif

# Add option for suppressing PGI non-uniform memory access (numa) library 
ifeq ($(NONUMA),yes) 
FFLAGS  += -mp=nonuma
endif

# Also add TOMAS aerosol microphysics option
ifeq ($(TOMAS),yes) 
FFLAGS  += -DTOMAS
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS  += -C
endif

# Include options (i.e. for finding *.h, *.mod files)
INCLUDE  = -I$(HDR) -module $(MOD)

# Also append HDF5 include commands if necessary
ifeq ($(HDF5),yes)
INCLUDE += -DUSE_HDF5 -I$(H5I)
endif

CC       = gcc
F90      = pgf90 $(FFLAGS) $(INCLUDE)
LD       = pgf90 $(FFLAGS)
FREEFORM = -Mfree
R8       = -Mextend -r8

endif

#==============================================================================
# SunStudio compilation options
#==============================================================================
ifeq ($(COMPILER),sun) 

# Pick compiler options for debug run or regular run 
# NOTE: -native builds in proper options for whichever chipset you have!
ifdef DEBUG 
FFLAGS   = -fpp -g -O0 -stackvar -xfilebyteorder=big16:%all -native
else
FFLAGS   = -fpp -fast -stackvar -xfilebyteorder=big16:%all -native
endif

# Build Sun for 32-bit platform
ifdef SUN32
FFLAGS  += -m32
else
FFLAGS  += -m64
endif

# Turn on OpenMP parallelization
ifeq ($(OMP),yes) 
FFLAGS  += -openmp=parallel -Dmultitask
endif

# Also add TOMAS aerosol microphysics option
ifeq ($(TOMAS),yes) 
FFLAGS  += -DTOMAS
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS  += -C
endif

# Include options (i.e. for finding *.h, *.mod files)
INCLUDE  = -I$(HDR) -moddir=$(MOD) -M$(MOD)

# Also append HDF5 include commands if necessary
ifeq ($(HDF5),yes)
INCLUDE += -DUSE_HDF5 -I$(H5I)
endif

CC       =
#---------------------------------------------------------------
# If your compiler is under the name "f90", use these lines!
F90      = f90 $(FFLAGS) $(INCLUDE)
LD       = f90 $(FFLAGS)
#---------------------------------------------------------------
# If your compiler is under the name "sunf90", use these lines!
#F90      = sunf90 $(FFLAGS) $(INCLUDE)
#LD       = sunf90 $(FFLAGS)
#---------------------------------------------------------------
FREEFORM = -free
R8       = -xtypemap=real:64

endif

#==============================================================================
# IBM/XLF compilation options
# NOTE: someone who runs on IBM compiler should check this !!!
#==============================================================================
ifeq ($(COMPILER),xlf) 

# Default compilation options
FFLAGS = -bmaxdata:0x80000000 -bmaxstack:0x80000000 -qfixed -qsuffix=cpp=f -q64

# Add optimization options
FFLAGS += -O3 -qarch=auto -qtune=auto -qcache=auto -qmaxmem=-1 -qstrict 

# Turn on OpenMP parallelization
ifeq ($(OMP),yes) 
FFLAGS += -qsmp=omp:opt -WF,-Dmultitask -qthreaded
endif

# Prior to 11/19/09:
## Add more options for parallel run
#ifndef DEBUG
#FFLAGS += -qsmp=omp:opt -WF,-Dmultitask -qthreaded
#endif

# Also add TOMAS aerosol microphysics option
ifeq ($(TOMAS),yes) 
FFLAGS  += -DTOMAS
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS += -C
endif

# Include options (i.e. for finding *.h, *.mod files)
INCLUDE  = -I$(HDR) -I $(MOD)

# Also append HDF5 include commands if necessary
ifeq ($(HDF5),yes)
INCLUDE += -DUSE_HDF5 -I$(H5I)
endif

CC       =
F90      = xlf90_r $(FFLAGS) $(INCLUDE)
LD       = xlf90_r $(FFLAGS)
FREEFORM = -qrealsize=8
R8       = -r8

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

#EOC
#==============================================================================
# Print variables for testing/debugging purposes (uncomment if necessary)
#==============================================================================
#headerinfo:
#	@@echo '####### in Makefile_header.mk ########' 
#	@@echo "compiler: $(COMPILER)"
#	@@echo "debug   : $(DEBUG)"
#	@@echo "bounds  : $(BOUNDS)"
#	@@echo "f90     : $(F90)"
#	@@echo "cc      : $(CC)"
