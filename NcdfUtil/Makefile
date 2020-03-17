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
#  See https://github.com/geoschem/geos-chem for complete history
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Directories
ROOT     :=..
BIN      :=$(ROOT)/bin
DOC      :=$(ROOT)/doc
LIB      :=$(ROOT)/lib
MOD      :=$(ROOT)/mod

# Executables
EXE      :=TestNcdfUtil.x

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# Include these files in the library libNcUtils.a
OBJlib   :=                  \
m_do_err_out.o             \
m_netcdf_io_checks.o	   \
m_netcdf_io_close.o	   \
m_netcdf_io_create.o	   \
m_netcdf_io_define.o	   \
m_netcdf_io_get_dimlen.o   \
m_netcdf_io_handle_err.o   \
m_netcdf_io_open.o	   \
m_netcdf_io_read.o	   \
m_netcdf_io_readattr.o	   \
m_netcdf_io_write.o        \
ncdf_mod.o

# But compile all of these files.  The modules charpak_mod.F and julday_mod.F
# are needed to build ncdf_mod.F90, but don't need to be added into library
# libNcUtils.a.  These modules already get packaged into libGeosUtil.a.
OBJECTS  :=charpak_mod.o julday_mod.o $(OBJlib)

# List of module files.  Convert to lowercase, then prefix directory name.
MODULES  :=$(OBJECTS:.o=.mod)
MODULES  :=$(shell echo $(MODULES) | tr A-Z a-z)
MODULES  :=$(foreach I,$(MODULES),$(MOD)/$(I))

# Compile this file to build the test driver
OBJECTSx :=$(EXE:.x=.o)

# Module file for the executable
MODULESx :=$(OBJECTSx:.o=.mod)
MODULESx :=$(shell echo $(MODULESx) | tr A-Z a-z)
MODULESx :=$(foreach I,$(MODULESx),$(MOD)/$(I))

# Library file
LIBRARY  :=libNcUtils.a

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

charpak_mod.o            : charpak_mod.F90

julday_mod.o             : julday_mod.F

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

ncdf_mod.o               : ncdf_mod.F90                 \
                           charpak_mod.F90              \
                           julday_mod.F                 \
                           m_netcdf_io_checks.o         \
                           m_netcdf_io_close.o          \
                           m_netcdf_io_create.o	        \
                           m_netcdf_io_define.o         \
                           m_netcdf_io_get_dimlen.o     \
                           m_netcdf_io_handle_err.o     \
                           m_netcdf_io_open.o           \
                           m_netcdf_io_read.o           \
                           m_netcdf_io_readattr.o       \
                           m_netcdf_io_write.o

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

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean help libnc ncdfcheck slowclean test

libnc: $(OBJECTS)
	$(AR) crs $(LIBRARY) $(OBJlib)
	mv $(LIBRARY) $(LIB)

lib: libnc

ncdfcheck: $(OBJECTSx)
	$(LD) $(OBJECTSx) -L$(LIB) -lNcUtils $(NC_LINK_CMD) -o $(EXE)
	./$(EXE)

check: ncdfcheck
test: ncdfcheck

clean:
	@echo "===> Making clean in directory: NcdfUtil <==="
	@rm -f *.o *.mod *.x *.a *.nc $(EXE)

slowclean:
	@echo "===> Making slowclean in directory: NcdfUtil <==="
	@rm -f $(OBJECTS) $(MODULES) $(LIBRARY) $(LIB)/$(LIBRARY)
	@rm -f $(OBJECTSx) $(MODULESx) $(EXE)

print:
	@echo $(LINK)

###############################################################################
###                                                                         ###
###  Help screen                                                            ###
###                                                                         ###
###############################################################################

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
	@echo 'DEBUG=yes      Compiles GEOS-Chem for use w/ a debugger'
	@echo 'BOUNDS=yes     Turns on subscript-array checking (for debugging purposes)'
	@echo 'TRACEBACK=yes  Compiles with -traceback flag (ifort only)'
#
#EOC
