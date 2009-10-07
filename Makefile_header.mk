# $Id: Makefile_header.mk,v 1.4 2009/10/07 14:49:48 bmy Exp $
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
# !REVISION HISTORY: 
#  16 Sep 2009 - R. Yantosca - Initial version
#  22 Sep 2009 - R. Yantosca - Bug fix, added -I$(HDR) to F90 compilation lines
#  24 Sep 2009 - R. Yantosca - added NONUMA option for PGI compiler
#  07 Oct 2009 - R. Yantosca - Replaced .SUFFIXES section w/ pattern rules
#EOP
#------------------------------------------------------------------------------
#BOC

# Make ifort the default compiler
ifndef COMPILER
COMPILER = ifort
endif

#------------------------------------------------------------------------------
# IFORT compilation options (default)
#------------------------------------------------------------------------------
ifeq ($(COMPILER),ifort) 

# Pick correct options for debug run or regular run 
ifdef DEBUG
FFLAGS = -cpp -w -noalign -convert big_endian -g -traceback 
else
FFLAGS = -cpp -w -O2 -auto -noalign -convert big_endian -openmp -Dmultitask
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
F90      = ifort $(FFLAGS) -I$(HDR) -I$(MOD)
LD       = ifort $(FFLAGS) -L$(LIB)
FREEFORM = -free
R8       = -r8

endif

#------------------------------------------------------------------------------
# Portland Group (PGI) compilation options
#------------------------------------------------------------------------------
ifeq ($(COMPILER),pgi) 

# Pick correct options for debug run or regular run 
ifdef DEBUG
FFLAGS = -byteswapio -Mpreprocess -fast -Bstatic
else
FFLAGS = -byteswapio -Mpreprocess -fast -mp -Mnosgimp -Dmultitask -Bstatic
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
F90      = pgf90 $(FFLAGS) -I$(HDR) -I$(MOD)
LD       = pgf90 $(FFLAGS) -L$(LIB)
FREEFORM = -Mfree
R8       = -Mextend -r8

endif

#------------------------------------------------------------------------------
# SunStudio compilation options
#------------------------------------------------------------------------------
ifeq ($(COMPILER),sun) 

# Default compilation options
# NOTE: -native builds in proper options for whichever chipset you have!
FFLAGS = -fpp -fast -stackvar -xfilebyteorder=big16:%all -native

# Additional flags for parallel run
ifndef DEBUG
FFLAGS += -openmp=parallel -Dmultitask
endif

# Add option for "array out of bounds" checking
ifdef BOUNDS
FFLAGS += -C
endif

CC       =
#---------------------------------------------------------------
# If your compiler is under the name "f90", use these lines!
F90      = f90 $(FFLAGS) -I$(HDR) -M$(MOD)
LD       = f90 $(FFLAGS) -L$(LIB)
#---------------------------------------------------------------
# If your compiler is under the name "sunf90", use these lines!
#F90      = sunf90 $(FFLAGS) -I$(HDR) -M$(MOD)
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

# Add more options for parallel run
ifndef DEBUG
FFLAGS += -qsmp=omp:opt -WF,-Dmultitask -qthreaded
endif

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
#	@echo '####### in Makefile_header.mk ########' 
#	@echo "compiler: $(COMPILER)"
#	@echo "debug   : $(DEBUG)"
#	@echo "bounds  : $(BOUNDS)"
#	@echo "f90     : $(F90)"
#	@echo "cc      : $(CC)"
