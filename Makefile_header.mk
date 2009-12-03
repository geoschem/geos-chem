# $Id: Makefile_header.mk,v 1.9 2009/12/03 21:34:39 bmy Exp $
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
# R8         Contains the command to force REAL -> REAL*8
#                                                                             .
# FFLAGS is a local variables that is not returned to the "outside world".
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
#EOP
#------------------------------------------------------------------------------
#BOC

# Make ifort the default compiler
ifndef COMPILER
COMPILER = ifort
endif

# Turn OpenMP on by default
ifndef OMP
OMP = yes
endif

#------------------------------------------------------------------------------
# IFORT compilation options (default)
#------------------------------------------------------------------------------
ifeq ($(COMPILER),ifort) 

# Turn on -traceback option by default for debugging runs
ifdef DEBUG
TRACEBACK=yes
endif

# Pick compiler options for debug run or regular run 
ifdef DEBUG
FFLAGS = -cpp -w -O0 -auto -noalign -convert big_endian -g
else
FFLAGS = -cpp -w -O2 -auto -noalign -convert big_endian
endif

# Turn on OpenMP parallelization
ifeq ($(OMP),yes) 
FFLAGS += -openmp -Dmultitask
endif

# Add special IFORT optimization commands
ifdef IPO
FFLAGS += -ipo -static
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS += -CB
endif

# Also add traceback option
ifdef TRACEBACK
FFLAGS += -traceback
endif

CC       =
F90      = ifort $(FFLAGS) -I$(HDR) -module $(MOD)
LD       = ifort $(FFLAGS) -L$(LIB)
FREEFORM = -free
R8       = -r8

endif

#------------------------------------------------------------------------------
# Portland Group (PGI) compilation options
#------------------------------------------------------------------------------
ifeq ($(COMPILER),pgi) 

# Pick compiler options for debug run or regular run 
ifdef DEBUG 
FFLAGS = -byteswapio -Mpreprocess -Bstatic -g -O0 
else
FFLAGS = -byteswapio -Mpreprocess -Bstatic -fast 
endif

# Turn on OpenMP parallelization
ifeq ($(OMP),yes) 
FFLAGS += -mp -Mnosgimp -Dmultitask
endif

# Add option for suppressing PGI non-uniform memory access (numa) library 
ifeq ($(NONUMA),yes) 
FFLAGS += -mp=nonuma
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS += -C
endif

CC       = gcc
F90      = pgf90 $(FFLAGS) -I$(HDR) -module $(MOD)
LD       = pgf90 $(FFLAGS) -L$(LIB)
FREEFORM = -Mfree
R8       = -Mextend -r8

endif

#------------------------------------------------------------------------------
# SunStudio compilation options
#------------------------------------------------------------------------------
ifeq ($(COMPILER),sun) 

# Pick compiler options for debug run or regular run 
# NOTE: -native builds in proper options for whichever chipset you have!
ifdef DEBUG 
FFLAGS = -fpp -g -O0 -stackvar -xfilebyteorder=big16:%all -native
else
FFLAGS = -fpp -fast -stackvar -xfilebyteorder=big16:%all -native
endif

ifdef SUN32
FFLAGS += -m32
endif

# Turn on OpenMP parallelization
ifeq ($(OMP),yes) 
FFLAGS += -openmp=parallel -Dmultitask
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS += -C
endif

CC       =
#---------------------------------------------------------------
# If your compiler is under the name "f90", use these lines!
F90      = f90 $(FFLAGS) -I$(HDR) -moddir=$(MOD) -M$(MOD)
LD       = f90 $(FFLAGS) -L$(LIB)
#---------------------------------------------------------------
# If your compiler is under the name "sunf90", use these lines!
#F90      = sunf90 $(FFLAGS) -I$(HDR) -moddir=$(MOD) -M$(MOD)
#LD       = sunf90 $(FFLAGS) -L$(LIB)
#---------------------------------------------------------------
FREEFORM = -free
R8       = -xtypemap=real:64

endif

#------------------------------------------------------------------------------
# IBM/XLF compilation options
# NOTE: someone who runs on IBM compiler should check this !!!
#------------------------------------------------------------------------------
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

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS += -C
endif

CC       =
F90      = xlf90_r $(FFLAGS) -I$(HDR) -I$(MOD)
LD       = xlf90_r $(FFLAGS) -L$(LIB)
FREEFORM = -qrealsize=8
R8       = -r8

endif

#------------------------------------------------------------------------------
# Specify pattern rules for compiliation 
# (i.e. tell "make" how to compile different types of source code files)
#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
# Export global variables so that the main Makefile will see these
#------------------------------------------------------------------------------
export CC
export F90
export FREEFORM
export LD
export R8

#EOC
#------------------------------------------------------------------------------
# Print variables for testing/debugging purposes (uncomment if necessary)
#------------------------------------------------------------------------------
#headerinfo:
#	@@echo '####### in Makefile_header.mk ########' 
#	@@echo "compiler: $(COMPILER)"
#	@@echo "debug   : $(DEBUG)"
#	@@echo "bounds  : $(BOUNDS)"
#	@@echo "f90     : $(F90)"
#	@@echo "cc      : $(CC)"