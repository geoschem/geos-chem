#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the GeosCore subdirectory)
#
# !DESCRIPTION: This is the main GEOS-Chem makefile.  It compiles the 
# GEOS-Chem core source code files and bundles all of the object files (*.o) 
# into the libGeosCore.a library (located in the LIB directory).  Module files 
# (*.mod) are copied to the MOD directory.
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
# %%% NOTE: Most of the time this Makefile will be called automatically    %%%
# %%% from the router Makefile in the top-level directory.  However, if    %%%
# %%% you are in the ./GeosCore directory, then you can call this Makefile %%%
# %%% to build the GEOS-Chem source code, libraries, and executables.      %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
# ROOTDIR    Specifies the root directory for the GEOS-Chem code
# BIN        Specifies the directory where executable files are stored
# BPCH       Specifies the directory where the G-C bpch routines are stored
# DOC        Specifies the directory for generating documentation w/ ProTeX
# EXE        Specifies the name of the executable file
# HDR        Specifies the directory where include files are found
# LIB        Specifies the directory where library files (*.a) are stored
# LINK       Specifies the link commands to the GEOS-Chem library files
# KPP        Specifies the directory where th KPP solver files reside
# MOD        Specifies the directory where module files (*.mod) are stored
# NCDF       Specifies the directory where netCDF utilities are stored
# OBJ        Specifies the list of object files (*.o) to be created.
# UTIL       Specifies the directory where the G-C utility modules are found
# AR         Sys var w/ name of library creator program (i.e., "ar", "ranlib")
# MAKE       Sys var w/ name of Make command (i.e, "make" or "gmake")
# NTRAC      Cmd line argument; specifies either 43 or 54 tracer simulation
# KPPSOLVER  Cmd line argument; specifies the type of integrator to use
#                                                                             .
# NOTE: CC, F90, FREEFORM, LD, R8 are included from "Makefile_header.mk".
#                                                                             .
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% You can compile GEOS-Chem in parallel using the "make -jN" option!   %%%
# %%%                                                                      %%%
# %%% N = number of proceses that you want to run simultaneously (i.e.     %%%
# %%% (when one file is finished compiling, "make" will immediately start  %%%
# %%% on the next one).  Usually N is the # of processors on your system.  %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                             .
# GEOS-Chem routines will be compiled in the following order, by directory:
# ----------------------------------------------------------------------------
# (1) NcdfUtil/            : NetCDF I/O modules
# (2) Headers/             : Header files (i.e. CMN_SIZE_mod.F, etc.)
# (3) KPP/                 : KPP solver routines
# (4) GeosUtil/            : GEOS-Chem utility modules (i.e. pressure_mod.F)
# (5) HEMCO/Core/          : HEMCO Core modules
# (6) HEMCO/Extensions/    : HEMCO Extensions modules
# (7) HEMCO/Interfaces/    : HEMCO Interface modules
# (7) ISORROPIA/           : ISORROPIA aerosol thermodyn equilibrium module
# (8) GeosCore/            : "Core" GEOS-Chem modules
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
ROOT :=..
HELP :=$(ROOT)/help
LIB  :=$(ROOT)/lib
MOD  :=$(ROOT)/mod

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOT)/Makefile_header.mk

# List of source files: everything ending in .F or .F90
SOURCES :=$(wildcard *.F) $(wildcard *.F90)

# List of object files (replace .F and .F90 with .o)
TMP     :=$(SOURCES:.F=.o)
OBJECTS :=$(TMP:.F90=.o)

# List of module files.  Convert to lowercase, then prefix directory name.
MODULES :=$(OBJECTS:.o=.mod)
MODULES :=$(shell echo $(MODULES) | tr A-Z a-z)
MODULES :=$(foreach I,$(MODULES),$(MOD)/$(I))

# Library file
LIBRARY :=libApm.a

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean debug slowclean

all: lib

lib: $(OBJECTS) 
	$(AR) crs $(LIBRARY) $(OBJECTS)
	mv $(LIBRARY) $(LIB)

clean:
	@echo "===> Making clean in directory: APM <==="
	@rm -f *.o *.mod *.a *.x

slowclean:
	@echo "===> Making slowclean in directory: History <==="
	@rm -f $(OBJECTS) $(MODULES) $(LIBRARY) $(LIB)/$(LIBRARY)

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

apm_ATHN_mod.o              : apm_ATHN_mod.F

apm_TIMN_mod.o              : apm_TIMN_mod.F

apm_icen_mod.o              : apm_icen_mod.F90

apm_rrtmg_mods.o            : apm_rrtmg_mods.F90

apm_rrtmg_src.o             : apm_rrtmg_src.F90                              \
                              apm_rrtmg_mods.o

apm_rrtmg_sw.o              : apm_rrtmg_sw.F90                               \
                              apm_rrtmg_mods.o        apm_rrtmg_src.o

apm_mixactivate.o           : apm_mixactivate.F90

apm_nucl_mod.o              : apm_nucl_mod.F

apm_opti_mod.o              : apm_opti_mod.F

apm_radf_mod.o              : apm_radf_mod.F

apm_albd_mod.o              : apm_albd_mod.F                                 \
                              apm_init_mod.o

apm_init_mod.o              : apm_init_mod.F                                 \
                              apm_nucl_mod.o                                 \
                              apm_opti_mod.o                                 \
                              apm_ATHN_mod.o                                 \
                              apm_TIMN_mod.o

apm_grow_mod.o              : apm_grow_mod.F                                 \
                              apm_init_mod.o

apm_coag_mod.o              : apm_coag_mod.F90                               \
                              apm_init_mod.o

apm_phys_mod.o              : apm_phys_mod.F                                 \
                              apm_coag_mod.o          apm_grow_mod.o         \
                              apm_init_mod.o          apm_nucl_mod.o         \
                              apm_opti_mod.o          apm_mixactivate.o      \
                              apm_ATHN_mod.o          apm_TIMN_mod.o

module_data_mosaic_asect.o  : module_data_mosaic_asect.F90

module_data_mosaic_therm.o  : module_data_mosaic_therm.F90

module_data_mosaic_other.o  : module_data_mosaic_other.F90

module_mosaic_therm.o       : module_mosaic_therm.F90                        \
                              module_data_mosaic_asect.o                     \
                              module_data_mosaic_therm.o                     \
                              module_data_mosaic_other.o

#EOC
