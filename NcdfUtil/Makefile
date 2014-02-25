#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in Ncdfutil subdirectory)
#
# !DESCRIPTION: Makefile for the NcdfUtilities package.
#\\
#\\
# !REMARKS:
# To build the program, call "make" with the following syntax:
#
#   make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#
# Makefile uses the following variables:
#
# Variable   Description
# --------   -----------
# BIN        Specifies the directory where the executable will be sent
# DOC        Specifies the directory where documentation will be built
# EXE        Specifies the name of the executable file
# F90        Contains the default Fortran compilation commands
# LD         Contains the command to link to libraries & make executable
# LIB        Specifies the directory where the libNcUtils.a files will be sent
# LINK_NC    Contains the comand to link to the HDF library
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
#
# F90, LD, LINK_NC are included from "Makefile_header.mk"
#
# !REVISION HISTORY: 
#  04 Aug 2009 - R. Yantosca - Initial version
#  24 Jan 2012 - R. Yantosca - Now explicitly list all dependencies
#  25 Jan 2012 - R. Yantosca - Now use := for the assignment statement
#  05 Apr 2012 - R. Yantosca - Change BL_INC_NETCDF to INC_NETCDF
#  05 Apr 2012 - R. Yantosca - Change BL_INC_HDF5   to INC_HDF5
#  05 Apr 2012 - R. Yantosca - Change BL_LIB_NETCDF to LIB_NETCDF
#  05 Apr 2012 - R. Yantosca - Change BL_LIB_HDF5   to LIB_HDF5
#  11 May 2012 - R. Yantosca - Now use exported NCL from Makefile_header.mk
#  20 Aug 2013 - R. Yantosca - Make sure Makefile target names are consistent
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================
# Define variables
ROOTDIR := ..
BIN     := $(ROOTDIR)/bin
DOC     := $(ROOTDIR)/doc
HDR     := $(ROOTDIR)/Headers
EXE     := TestNcdfUtil.x
LIB     := $(ROOTDIR)/lib
MOD     := $(ROOTDIR)/mod

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

#==============================================================================
# Filename section
#
# NOTES:
# (1) The function $(wildcard m_*.F90) returns all files matching m_*.F90, etc.
# (2) The syntax $(SRC:.F90=.o) is a shorthand for the "patsubst" function. 
#==============================================================================

#-----------------
# Objects
#-----------------

# Files for library build (code & objects)
SRC  = $(wildcard m_*.F90)
OBJ  = $(SRC:.F90=.o)
MOD  = $(SRC:.F90=.mod)

# Files for standalone test driver (code, object, exe)
OBJx = $(EXE:.x=.o)

#-----------------
# Dependencies
#-----------------
m_do_err_out.o           : m_do_err_out.F90

m_netcdf_io_checks.o     : m_netcdf_io_checks.F90

m_netcdf_io_close.o      : m_netcdf_io_close.F90        \
                           m_do_err_out.o

m_netcdf_io_create.o     : m_netcdf_io_create.F90       \
                           m_do_err_out.o

m_netcdf_io_define.o     : m_netcdf_io_define.F90       \
                           m_do_err_out.o

m_netcdf_io_get_dimlen.o : m_netcdf_io_get_dimlen.F90   \
                           m_do_err_out.o

m_netcdf_io_handle_err.o : m_netcdf_io_handle_err.F90   \
                           m_do_err_out.o

m_netcdf_io_open.o       : m_netcdf_io_open.F90         \
                           m_do_err_out.o

m_netcdf_io_read.o       : m_netcdf_io_read.F90         \
                           m_do_err_out.o

m_netcdf_io_readattr.o   : m_netcdf_io_readattr.F90     \
                           m_do_err_out.o

m_netcdf_io_write.o      : m_netcdf_io_write.F90        \
                           m_do_err_out.o

TestNcdfUtil.o           : TestNcdfUtil.F90             \
                           m_do_err_out.o               \
                           m_netcdf_io_create.o         \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_write.o          \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_read.o           \
                           m_netcdf_io_readattr.o

#==============================================================================
# Makefile targets
#==============================================================================

.PHONY: clean help libnc ncdfcheck test

libnc: $(OBJ)
	$(AR) crs libNcUtils.a $(OBJ)
	mv libNcUtils.a $(LIB)

lib: libnc

ncdfcheck: $(OBJx)
	$(LD) $(OBJx) -L$(LIB) -lNcUtils $(NCL) -o $(EXE)
	$(EXE) 

check: ncdfcheck
test: ncdfcheck

clean:
	rm -f *.a *.mod *.o *.nc *.x 

print:
	@echo $(LINK)

#==============================================================================
# HELP Screen
#==============================================================================

help:
	@echo 'Usage: make TARGET [ OPTIONAL-FLAGS ]'
	@echo ''
	@echo 'TARGET may be one of the following:'
	@echo 'lib        Compiles code and builds the libNcUtils.a library file'
	@echo 'all        Synonym for "lib" (default target)'
	@echo 'check      Compiles and runs a test program to check the libNcUtils.a library'
	@echo 'test       Synonym for "check"'
	@echo 'clean      Removes *.o *.mod and executable files in this subdirectory'
	@echo 'help       Displays this help screen'
	@echo ''
	@echo 'OPTIONAL-FLAGS may be:'
	@echo 'COMPILER=___   Options: ifort pgi sun xlf (default is ifort)'
	@echo 'DEBUG=yes      Compiles GEOS-Chem for use w/ a debugger'
	@echo 'BOUNDS=yes     Turns on subscript-array checking (for debugging purposes)'
	@echo 'TRACEBACK=yes  Compiles with -traceback flag (ifort only)'
#
#EOC
