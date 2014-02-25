#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the ISOROPIA/ subdirectory)
#
# !DESCRIPTION: This makefile compiles the ISOROPIA code.  
#  Object files (*.o) are bundled into the libIsoropia.a library 
#  (located in the LIB directory).  Module files (*.mod) are copied to 
#  the MOD directory. 
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% NOTE: Normally you will not have to call this Makefile directly,     %%%
# %%% it will be called automatically from the Makefile in the directory   %%%
# %%% just above this one!                                                 %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
# ROOTDIR    Specifies the root-level directory of the GEOS-Chem code
# HDR        Specifies the directory where GEOS-Chem include files are found
# LIB        Specifies the directory where library files (*.a) are stored
# MOD        Specifies the directory where module files (*.mod) are stored
# AR         Sys var w/ name of library creator program (i.e., "ar", "ranlib")
# MAKE       Sys var w/ name of Make command (i.e, "make" or "gmake")
# R8         Specifies the command to treat "REAL" as "REAL*8"
#
# !REVISION HISTORY: 
#  21 Dec 2009 - C. Carouge  - Initial version
#  22 Aug 2011 - R. Yantosca - Add "-fp-model source" flag for IFORT compiler,
#                              which prevents random numerical noise
#  25 Aug 2011 - R. Yantosca - Remove -fp-model source flag here, as this is
#                              now added to FFLAGS in Makefile_header.mk
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
SHELL   = /bin/sh
ROOTDIR = ..
HDR     = $(ROOTDIR)/Headers
HELP    = $(ROOTDIR)/help
LIB     = $(ROOTDIR)/lib
MOD     = $(ROOTDIR)/mod

# Include header file.  This returns variables CC, F90, FREEFORM, LD, R8,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

#=============================================================================
# List of files to compile (the order is important!).  We specify these as
# a list of object files (*.o).  For each object file, the "make" utility
# will find the corresponding source code file (*.F) and compile it. 
#=============================================================================

# List of source files
SRC = $(wildcard *.F) $(wildcard *.F90)

# Replace .F and .F90 extensions with *.o
TMP = $(SRC:.F=.o)
OBJ = $(TMP:.F90=.o)

#=============================================================================
# Makefile targets: type "make help" for a complete listing!
#=============================================================================

.PHONY: clean help

lib: $(OBJ)
	$(AR) crs libIsoropia.a $(OBJ)
	mv libIsoropia.a $(LIB)

clean:
	rm -f *.o *.mod

help:
	@$(MAKE) -C $(HELP)

#=============================================================================
# Dependencies listing (grep "USE " to get the list of module references!)
#
# From this list of dependencies, the "make" utility will figure out the
# correct order of compilation (so we don't have to do that ourselves!)
#=============================================================================

isoropiaIIcode.o:  isoropiaIIcode.F  isrpia.inc
	$(F90) $(R8) -c $<
#EOC

