#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
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
# \item Converting gas concentration units
# \end{itemize}
#
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
#
# !REVISION HISTORY: 
#  19 Nov 2009 - R. Yantosca - Initial version
#  23 Nov 2009 - R. Yantosca - Now don't copy module files; they will be
#                              automatically written to the mod directory
#  11 Dec 2009 - R. Yantosca - Now get SHELL from Makefile_header.mk
#  21 Dec 2009 - R. Yantosca - If HDF5=yes, then look for hdf5.mod in the
#                              HDF5 include path $(HDF5_INC).
#  01 Mar 2012 - R. Yantosca - Replace grid_mod.F with grid_mod.F90, to
#                              facilitate work on the GI model
#  03 Apr 2012 - M. Payer    - Add new module regrid_a2a_mod.F90 (M. Cooper)
#  03 Aug 2012 - R. Yantosca - Add dependency for inquireMod.F90
#  19 Mar 2014 - R. Yantosca - Add more visible comment section dividers
#  20 Jun 2014 - R. Yantosca - Removed unix_cmds_mod.F; it's now obsolete
#  20 Jun 2014 - R. Yantosca - Removed directory_mod.F; it's now obsolete
#  10 Jul 2014 - R. Yantosca - Now compile ncdf_mod.F90 in NcdfUtil/ subdir
#  21 Jul 2014 - R. Yantosca - regrid_a2a_mod.F90 no longer relies on other
#                              files in GeosUtil
#  21 Jul 2014 - R. Yantosca - Removed regrid_1x1_mod.F; it's obsolete
#  23 Jul 2014 - R. Yantosca - Removed global_grid_mod.F90; it's obsolete
#  08 Jan 2015 - E. Lundgren - Add module for unit conversion
#  03 Jun 2015 - R. Yantosca - Also remove *.mod, *.a files with "make clean"
#  07 Dec 2015 - R. Yantosca - Restore fast "clean" command; add "slowclean"
#  22 Jan 2016 - R. Yantosca - Remove linux_err.o
#  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90, etc.
#  23 Aug 2017 - R. Yantosca - Added dependencies for grid_registry_mod.F90
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Define variables
ROOT      :=..
HELP      :=$(ROOT)/help
LIB       :=$(ROOT)/lib
MOD       :=$(ROOT)/mod

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# List of source files: everything ending in .F and .F90
SOURCES   :=$(wildcard *.F) $(wildcard *.F90)

# Replace .f and .f90 extensions with *.o
TMP       :=$(SOURCES:.F=.o)
OBJECTS   :=$(TMP:.F90=.o)

# Special files just for IFORT
ifeq ($(COMPILER),ifort) 
  OBJECTS += ifort_errmsg.o
endif

# List of module files.  Convert to lowercase, then prefix directory name.
MODULES   :=$(OBJECTS:.o=.mod)
MODULES   :=$(shell echo $(MODULES) | tr A-Z a-z)
MODULES   :=$(foreach I,$(MODULES),$(MOD)/$(I))

# Library file
LIBRARY   :=libGeosUtil.a

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean help debug slowclean

lib: $(OBJECTS)
	$(AR) crs $(LIBRARY) $(OBJECTS)
	mv $(LIBRARY) $(LIB)

clean:
	@echo "===> Making clean in directory: GeosUtil <==="
	@rm -f *.o *.mod *.a *.x

slowclean:
	@echo "===> Making slowclean in directory: GeosUtil <==="
	@rm -f $(OBJECTS) $(MODULES) $(LIBRARY) $(LIB)/$(LIBRARY)

help:
	@$(MAKE) -C $(HELP)

debug:
	@echo "Targets : $(MAKECMDGOALS)"
	@echo "ROOT    : $(ROOT)"
	@echo "LIB     : $(LIB)"
	@echo "MOD     : $(MOD)"
	@echo "F90     : $(F90)"
	@echo "OBJECTS : $(OBJECTS)"
	@echo "MODULES : $(MODULES)"
	@echo "LIBRARY : $(LIBRARY)"

###############################################################################
###                                                                         ###
###  Dependencies listing                                                   ###
###  (grep "USE " to get the list of module references!)                    ###
###                                                                         ###
###  From this list of dependencies, the "make" utility will figure out     ###
###  correct order of compilation (so we don't have to do that ourselves).  ###
###  This also allows us to compile on multiple processors with "make -j".  ###
###                                                                         ###
###  NOTES:                                                                 ###
###  (1) Only specify object-file dependencies that are within this         ###
###       directory.  Object files in other directories will be referenced  ### 
###       at link-time.                                                     ###
###  (2) For "make -jN" (i.e. compile N files simultaneously), all files    ###
###       in this directory must have a listed dependency.                  ###
###                                                                         ###
###############################################################################

bpch2_mod.o         : bpch2_mod.F           \
                      error_mod.o           \
                      file_mod.o            \
                      julday_mod.o          

error_mod.o         : error_mod.F           \
                      geos_timers_mod.o

file_mod.o          : file_mod.F            \
                      error_mod.o

gc_grid_mod.o       : gc_grid_mod.F90       \
                      error_mod.o           \
                      roundoff_mod.o

grid_registry_mod.o : grid_registry_mod.F90 \
                      gc_grid_mod.o         \
                      pressure_mod.o

henry_mod.o         : henry_mod.F

ifort_errmsg.o      : ifort_errmsg.F

julday_mod.o        : julday_mod.F

pressure_mod.o      : pressure_mod.F        \
                      error_mod.o

regrid_a2a_mod.o    : regrid_a2a_mod.F90

roundoff_mod.o      : roundoff_mod.F90

time_mod.o          : time_mod.F            \
                      error_mod.o           \
                      gc_grid_mod.o         \
                      julday_mod.o 

transfer_mod.o      : transfer_mod.F        \
                      error_mod.o

unitconv_mod.o	    : unitconv_mod.F90      \
                      error_mod.o           \
                      geos_timers_mod.o

geos_timers_mod.o   : geos_timers_mod.F90

#%%% NOTE: THIS IS NOW OBSOLETE, WE SHOULD REMOVE IT %%% 
# NOTE: if HDF5=yes then we need to also look for the hdf5.mod 
# in the HDF5 include path (bmy, 12/21/09)
hdf_mod.o           : hdf_mod.F             \
                      error_mod.o           \
                      gc_grid_mod.o

ifeq ($(HDF5),yes)
	$(F90) -DUSE_HDF5 -I$(HDF_INC) -c $<
endif

#EOC
