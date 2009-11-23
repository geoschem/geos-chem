# $Id: Makefile,v 1.2 2009/11/23 21:44:57 bmy Exp $
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the GeosUtil subdirectory)
#
# !DESCRIPTION: This makefile compiles the various GEOS-Chem utility modules,
#  which provide basic functionality for:
#
# \begin{itemize}
# \item Collapsing vertical levels in the stratosphere
# \item Date and time computations
# \item Defining data directories
# \item Defining the GEOS-Chem horizontal grid
# \item Defining the GEOS-Chem pressure coordinate grid
# \item Defining the logical units for GEOS-Chem file I/O
# \item Defining various Unix commands
# \item Platform-specific error handling
# \item Manipulating string variables
# \item Regridding data (horizontally) from fine to coarse resolution
# \end{itemize}
#
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make TARGET [ OPTIONAL-FLAGS ]
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
#
# !REVISION HISTORY: 
#  19 Nov 2009 - R. Yantosca - Initial version
#  23 Nov 2009 - R. Yantosca - Now don't copy module files; they will be
#                              automatically written to the mod directory
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
# List of files to compile.  Here the order is not important,
# as we will explicity define the dependencies listing below.
#=============================================================================

# Source file listing
SRC   = $(wildcard *.f*)

# Special files just for IFORT
ifeq ($(COMPILER),ifort) 
SRC += ifort_errmsg.f
endif

# Special files just for PGI
ifeq ($(COMPILER),pgi) 
SRC += linux_err.f
endif

# Object file listing
OBJ   = $(SRC:.f=.o)

#=============================================================================
# Makefile targets: type "make help" for a complete listing!
#=============================================================================

.PHONY: clean help

lib: $(OBJ)
	$(AR) crs libGeosUtil.a $(OBJ)
	mv libGeosUtil.a $(LIB)

clean:
	rm -f *.o *.mod

help:
	@$(MAKE) -C $(HELP)

#=============================================================================
# Dependencies listing (grep "USE " to get the list of module references!)
#
# From this list of dependencies, the "make" utility will figure out the
# correct order of compilation (so we don't have to do that ourselves).
# This also allows us to compile on multiple processors with "make -j".
#
# NOTES:
# (1) Only specify object-file dependencies that are within this directory.  
#      Object files in other directories will be referenced at link-time.
# (2) For "make -j" to work, all files in this directory must have a 
#      listed dependency.
#=============================================================================

bpch2_mod.o      : bpch2_mod.f      error_mod.o    file_mod.o   julday_mod.o
charpak_mod.o    : charpak_mod.f
directory_mod.o  : directory_mod.f
error_mod.o      : error_mod.f
file_mod.o       : file_mod.f       error_mod.o
grid_mod.o       : grid_mod.f       error_mod.o
ifort_errmsg.o   : ifort_errmsg.f
julday_mod.o     : julday_mod.f
linux_err.o      : linux_err.c
pressure_mod.o   : pressure_mod.f   error_mod.o
regrid_1x1_mod.o : regrid_1x1_mod.f charpak_mod.o  error_mod.o  grid_mod.o
time_mod.o       : time_mod.f       charpak_mod.o  error_mod.o  grid_mod.o  \
                                    julday_mod.o 
transfer_mod.o   : transfer_mod.f   error_mod.o
unix_cmds_mod.o  : unix_cmds_mod.f

#EOC
