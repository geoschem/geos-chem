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
# (7) ISOROPIA/            : ISORROPIA aerosol thermodyn equilibrium module
# (8) GeosCore/            : "Core" GEOS-Chem modules
#
# !REVISION HISTORY: 
#  16 Sep 2009 - R. Yantosca - Initial version
#  18 Sep 2009 - P. Le Sager - Removed -lKppInt
#  21 Sep 2009 - R. Yantosca - Now call Makefile in help directory to
#                              display the help screen options
#  19 Nov 2009 - R. Yantosca - Now compile the various GEOS-Chem utility
#                              modules in the GeosUtil subdirectory
#  19 Nov 2009 - R. Yantosca - Now compile the GEOS-Chem bpch module 
#                              separately in the GeosBpch subdirectory
#  19 Nov 2009 - R. Yantosca - Now list all object dependencies explicitly,
#                              to be able to use "make -j" (parallel make)
#  23 Nov 2009 - R. Yantosca - Remove "main.o" explicitly form the "exe"
#                              makefile target.  This will be now compiled in
#                              the proper sequence given the dependency
#                              ordering.  This allows "make -j" to work. 
#  23 Nov 2009 - R. Yantosca - Now don't copy module files; they will be
#                              automatically written to the mod directory
#  23 Nov 2009 - R. Yantosca - Removed libbpch; bundled that into libutil
#  23 Nov 2009 - R. Yantosca - Added separate target libcore.  lib is now
#                              a synonym for "libkpp libutil libcore"
#  01 Dec 2009 - R. Yantosca - Modified the "exe" target for SunStudio compiler
#                              which chokes at link time if the list of object
#                              files is not explicitly passed
#  02 Dec 2009 - R. Yantosca - Add conditional statements in dependencies list
#                              for the SunStudio compiler
#  11 Dec 2009 - R. Yantosca - Now get SHELL from Makefile_header.mk
#  21 Dec 2009 - R. Yantosca - Now get LINK from Makefile_header.mk
#  25 Jan 2010 - R. Yantosca - Whem making "realclean", also call "clean" in
#                              the $(GEOSTOM)/Makefile.  Also make sure to 
#                              remove executables in the $(BIN) directory.
#  28 Jan 2010 - C. Carouge  - Modifications for ISORROPIA II
#  08 Feb 2010 - C. Carouge  - Modifications for F. Paulot's isoprene scheme
#  10 May 2010 - R. Yantosca - Add dependency for RD_AOD.f
#  14 May 2010 - C. Carouge  - Updates for mercury simulation
#  20 Aug 2010 - R. Yantosca - Modifications for MERRA met fields
#  16 Feb 2011 - R. Yantosca - Add modifications for APM (G. Luo)
#  05 Aug 2011 - M. Long     - Now compile module files in Headers/ directory
#  04 Nov 2011 - R. Yantosca - Remove references to ESMF subdirectory
#  08 Dec 2011 - M. Payer    - Remove obsolete GEIA biogenic emissions routines
#  24 Jan 2012 - R. Yantosca - Also add libnc target to build netCDF utils
#  25 Jan 2012 - R. Yantosca - Add ncdfcheck target to check netCDF install
#  05 Apr 2012 - R. Yantosca - Now assume netCDF will always be used
#  05 Apr 2012 - R. Yantosca - Now retire rdlai.F, readlai.F
#  11 Apr 2012 - R. Yantosca - Now retire obsolete lai_mod.F
#  11 Apr 2012 - R. Yantosca - Reference modis_lai_mod.F90 in timeseries diag
#  11 Apr 2012 - R. Yantosca - Remove all references to obsolete lai_mod.F
#  12 Apr 2012 - R. Yantosca - Remove reference to findmon.F
#  19 Apr 2012 - R. Yantosca - Remove reference to rd_prof.F
#  30 Oct 2012 - R. Yantosca - Remove reference to rdsoil.F
#  27 Nov 2012 - R. Yantosca - Remove reference to getifsun.F
#  27 Mar 2013 - S.D. Eastham- Rolled photolysis functions together
#  30 May 2013 - R. Yantosca - Remove reference to GEOS-3 TPCORE routines 
#  13 Aug 2013 - M. Sulprizio- Remove reference to soaprod_mod.F (no longer
#                              needed) (SOAupdate, hotp 7/25/10)
#  20 Aug 2013 - R. Yantosca - Make sure Makefile target names are consistent
#  20 Sep 2013 - R. Yantosca - Bug fix for make realclean: Don't try to make
#                              clean unless we compile for target hpc
#  27 Jan 2014 - R. Yantosca - Add tomas_mod.o to wetscav_mod dependency list
#  03 Feb 2014 - R. Yantosca - Add linoz_mod.o to input_mod.o dependency list
#  21 Feb 2014 - M. Sulprizio- Removed reference to tropopause_mod.F. It is now
#                              obsolete because of the UCX updates.
#  18 Mar 2014 - R. Yantosca - Add tauclean target to remove TAU profile files
#  19 Mar 2014 - R. Yantosca - Now restore GTMM functionality
#  19 Mar 2014 - R. Yantosca - Add more visible comment section dividers
#  11 Apr 2014 - R. Yantosca - Add modules to input_mod.o dependencies list
#  23 Jun 2014 - R. Yantosca - Remove references to obsolete logical_mod.o
#  21 Jul 2014 - R. Yantosca - Now build HEMCO emissions modules after
#                              NcdfUtil and GeosUtil, but before GeosCore
#  22 Jul 2014 - R. Yantosca - Move hcox_paranox_mod.F90, paranox_util_mod.F90
#                              and hcox_driver_mod.F90 back to HEMCO/Extensions
#  24 Jul 2014 - R. Yantosca - Remove canopy_nox_mod.F, readfert.F, readclim.F,
#                              and readsoil.F. HEMCO makes all these obsolete.
#  15 Sep 2014 - M. Sulprizio- Remove global_oc_mod.F and global_bc_mod.F.
#                              HEMCO makes these obsolete.
#  23 Sep 2014 - M. Sulprizio- Removed global_hno3_mod.F and global_no3_mod.F.
#                              HEMCO makes these obsolete.
#  04 Jun 2015 - R. Yantosca - Add wipeout target to remove any leftover
#                              library, module, or executable files
#  03 Sep 2015 - R. Yantosca - Add gamap_mod.o dependency to input_mod.o
#  03 Dec 2015 - R. Yantosca - cleanup.o now depends on geosfp_read_mod.o
#                              and merra2_read_mod.o
#  07 Dec 2015 - R. Yantosca - Restore fast "clean" command; add "slowclean"
#  07 Dec 2015 - R. Yantosca - Add "realclean_except_rrtmg" target
#  18 May 2016 - M. Sulprizio- Remove comode_mod.F, setemis.F, and several
#                              leftover routines from SMVGEAR
#  23 Jun 2016 - R. Yantosca - Remove references to tracerid_mod
#  23 Jun 2016 - R. Yantosca - Remove references to GeosApm
#  30 Jun 2016 - M. Sulprizio- Remove readchem.F, reader.F, jsparse.F, ksparse.F
#  09 Aug 2016 - M. Sulprizio- Remove diag_pl_mod and add diag20_mod
#  13 Jul 2017 - E. Lundgren - Remove passive_species_mod.F90
#  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
#  23 Aug 2018 - H.P. Lin    - If NO_EXE is set (i.e. EXE_NEEDED=0) for coupled
#                              model driving G-C do not build "geos" but instead
#                              "libGeosCore.a"
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Directories
ROOT    :=..
BIN     :=$(ROOT)/bin
BPCH    :=$(ROOT)/GeosBpch
DOC     :=$(ROOT)/doc
GCHP    :=$(ROOT)/GCHP
GTMM    :=$(ROOT)/GTMM
HEMCO   :=$(ROOT)/HEMCO
HDR     :=$(ROOT)/Headers
HELP    :=$(ROOT)/help
HIST    :=$(ROOT)/History
ISO     :=$(ROOT)/ISOROPIA
LIB     :=$(ROOT)/lib
KPP     :=$(ROOT)/KPP
MOD     :=$(ROOT)/mod
NCDF    :=$(ROOT)/NcdfUtil
OBSP    :=$(ROOT)/ObsPack
RAD     :=$(ROOT)/GeosRad
UTIL    :=$(ROOT)/GeosUtil
DYN     :=$(GCHP)/FVdycoreCubed_GridComp

# Executables
EXE     :=geos
TOM     :=geostomas

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
LIBRARY :=libGeosCore.a

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean realclean doc docclean help distclean tauclean wipeout
.PHONY: slowclean slowrealclean

all:                                               # Build libraries
	@$(MAKE) lib                               #  and the executable
ifeq ($(UNAME),Darwin)
	ranlib -c ../lib/*
endif
ifeq ($(EXE_NEEDED),1)
	@$(MAKE) exe
endif

tomas:
	@$(MAKE) TOMAS=yes lib                     # TOMAS, 30 bins (default)
	@$(MAKE) TOMAS=yes tomexe

tomas40:
	@$(MAKE) TOMAS40=yes lib                   # TOMAS, 40 bins
	@$(MAKE) TOMAS40=yes tomexe

tomas12:
	@$(MAKE) TOMAS12=yes lib                   # TOMAS, 12 bins
	@$(MAKE) TOMAS12=yes tomexe

tomas15:
	@$(MAKE) TOMAS15=yes lib                   # TOMAS, 15 bins
	@$(MAKE) TOMAS15=yes tomexe

hpc:						   # Build ESMF & MAPL
	@$(MAKE) baselibs			   # compile for HPC environment
	@$(MAKE) lib
	@$(MAKE) libgigc
ifeq ($(UNAME),Darwin)
	ranlib -c ../lib/*
endif
	@$(MAKE) exe

lib:                                               # Build all G-C libraries
	@$(MAKE) libnc
	@$(MAKE) kppfirstpass
	@$(MAKE) libheaders         		  
	@$(MAKE) libutil 
	@$(MAKE) libkpp
	@$(MAKE) libhistory
	@$(MAKE) libobspack
	@$(MAKE) libhemco
	@$(MAKE) libiso
	@$(MAKE) librad
	@$(MAKE) libcore

libcore: $(OBJECTS)                                # Build code in GeosCore/
ifeq ($(EXE_NEEDED),0)
	$(AR) crs $(LIBRARY) $(OBJECTS)
	mv $(LIBRARY) $(LIB) # ...build a libGeosCore.a for use instead
endif


libgigc:
ifeq ($(wildcard $(FVDIR)/fvdycore.install),)
	@$(MAKE) -C $(DYN) install
	@touch $(FVDIR)/fvdycore.install
endif
	@$(MAKE) -C $(GCHP) lib

libheaders:                                        # Build code in Headers/ 
	@$(MAKE) -C $(HDR)			  

libhemco:                                          # Build code in HEMCO/*
	@$(MAKE) -C $(HEMCO) lib

libhemcosa:                                        # Build HEMCO standalone
	@$(MAKE) libnc
	@$(MAKE) kppfirstpass
	@$(MAKE) libheaders         		  
	@$(MAKE) libutil 
	@$(MAKE) libhemco

libhistory:                                        # Build code in History/
	@$(MAKE) -C $(HIST) lib

libiso:                                            # Build code in ISOROPIA/
	@$(MAKE) -C $(ISO)			  

librad:                                            # Build code in GeosRad/
ifeq ($(RRTMG_NEEDED),1)
	@$(MAKE) -C $(RAD)			  
endif

libkpp:                                            # Build code in KPP/
	@$(MAKE) -C $(KPP)			  

kppfirstpass:                                      # Compile some KPP code 1st
	@$(MAKE) -C $(KPP) firstpass

libnc:                                             # Build code in NcdfUtil/
	@$(MAKE) -C $(NCDF) libnc

libobspack:                                        # Build code in ObsPack/
	@$(MAKE) -C $(OBSP) lib

libutil:                                           # Build code in GeosUtil/
	@$(MAKE) -C $(UTIL)			  

baselibs:
	@$(MAKE) -C $(GCHP) baselibs

ncdfcheck:                                         # Check netCDF library
	@$(MAKE) libnc
	@$(MAKE) -C $(NCDF) ncdfcheck

exe:                                               # Build executable
	$(LD) $(OBJECTS) $(LINK) -o $(EXE)
	cp -f $(EXE) $(BIN)

tomexe:                                            # Build executable
	$(LD) $(OBJECTS) $(LINK) -o $(TOM)
	cp -f $(TOM) $(BIN)

clean:                                             # Remove files here
	@echo "===> Making clean in directory: GeosCore <==="
	@rm -f *.o *.mod *.a *.x $(EXE) $(TOM)

slowclean:                                         # Remove files here
	@echo "===> Making slowclean in directory: GeosCore <==="
	@rm -f $(OBJECTS) $(MODULES) $(LIBRARY) $(LIB)/$(LIBRARY)
	@rm -f $(EXE) $(BIN)/$(EXE)
	@rm -f $(TOM) $(BIN)/$(TOM)

distclean:                                         # Synonym for "realclean"
	@$(MAKE) realclean

realclean:                                         # Remove files everywhere
	@$(MAKE) clean
	@$(MAKE) -C $(GTMM)  clean
	@$(MAKE) -C $(HEMCO) clean
	@$(MAKE) -C $(HDR)   clean
	@$(MAKE) -C $(HIST)  clean
	@$(MAKE) -C $(OBSP)  clean
	@$(MAKE) -C $(ISO)   clean
	@$(MAKE) -C $(KPP)   realclean
	@$(MAKE) -C $(NCDF)  clean
	@$(MAKE) -C $(RAD)   clean
	@$(MAKE) -C $(UTIL)  clean
	@$(MAKE) docclean
ifeq ($(HPC),yes)
	@$(MAKE) -C $(GCHP)  clean
endif
	@$(MAKE) wipeout

wipeout:
	rm -f $(LIB)/*.a
	rm -f $(MOD)/*.mod
	rm -f $(BIN)/geos* $(BIN)/*.x

realclean_except_rrtmg:                            # Don't realclean RRTMG
	@$(MAKE) slowclean
	@$(MAKE) -C $(GTMM)  slowclean
	@$(MAKE) -C $(HEMCO) slowclean
	@$(MAKE) -C $(HDR)   slowclean
	@$(MAKE) -C $(HIST)  slowclean
	@$(MAKE) -C $(OBSP)  slowclean
	@$(MAKE) -C $(ISO)   slowclean
	@$(MAKE) -C $(KPP)   slowrealclean
	@$(MAKE) -C $(NCDF)  slowclean
	@$(MAKE) -C $(UTIL)  slowclean
	@$(MAKE) docclean
ifeq ($(HPC),yes)
	@$(MAKE) -C $(GCHP) clean
endif

doc:                                               # Build documentation
	@$(MAKE) -C $(DOC) doc

docclean:                                          # Remove documentation
	@$(MAKE) -C $(DOC) clean

help:                                              # Show help screen
	@$(MAKE) -C $(HELP) help

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
###  Targets for Hg simulation with Global Terrestrial Mercury Model        ###
###                                                                         ###
###############################################################################

allhg:                                             # Build Hg/GTMM executable
	@$(MAKE) GTMM_Hg=yes libhg
	@$(MAKE) GTMM_Hg=yes exehg

libhg:                                             # Compile HG/GTMM code
	@$(MAKE) GTMM_Hg=yes libgtmm           		     
	@$(MAKE) GTMM_Hg=yes lib

libgtmm:                                           # Compile GTMM model code
	@$(MAKE) -C $(GTMM) GTMM_Hg=yes

exehg:                                             # Create Hg/GTMM exec file
	$(LD) $(OBJECTS) $(LINK) -o $(EXE)
	cp -f $(EXE) $(BIN)

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

aerosol_mod.o               : aerosol_mod.F                                  \
                              dao_mod.o                                      \
                              diag_mod.o              ucx_mod.o

benchmark_mod.o             : benchmark_mod.F                              

bromocarb_mod.o             : bromocarb_mod.F         dao_mod.o              \
                              pbl_mix_mod.o

c2h6_mod.o                  : c2h6_mod.F                                     \
                              diag_mod.o              hco_interface_mod.o

carbon_mod.o                : carbon_mod.F                                   \
                              diag_mod.o              drydep_mod.o           \
                              pbl_mix_mod.o                                  \
                              tomas_mod.o             vdiff_pre_mod.o        \
                              hco_interface_mod.o     aerosol_mod.o          \
                              dao_mod.o

chemistry_mod.o             : chemistry_mod.F90                              \
                              aerosol_mod.o           isoropiaII_mod.o       \
                              c2h6_mod.o              carbon_mod.o           \
                              dust_mod.o              drydep_mod.o           \
                              global_ch4_mod.o        mercury_mod.o          \
                              pops_mod.o                                     \
                              rpmares_mod.o           RnPbBe_mod.o           \
                              seasalt_mod.o           strat_chem_mod.o       \
                              sulfate_mod.o           tagged_co_mod.o        \
                              tagged_o3_mod.o         tomas_mod.o            \
                              flexchem_mod.o          diagnostics_mod.o

cleanup.o                   : cleanup.F                                      \
                              aerosol_mod.o                                  \
                              bromocarb_mod.o         c2h6_mod.o             \
                              carbon_mod.o            co2_mod.o              \
                              dao_mod.o               depo_mercury_mod.o     \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag50_mod.o            diag51_mod.o           \
                              diag51b_mod.o           diag53_mod.o           \
                              diag_oh_mod.o           diag20_mod.o           \
                              drydep_mod.o            dust_mod.o             \
                              global_ch4_mod.o        wetscav_mod.o          \
                              isoropiaII_mod.o        land_mercury_mod.o     \
                              seasalt_mod.o           linoz_mod.o            \
                              strat_chem_mod.o        mercury_mod.o          \
                              modis_lai_mod.o         ocean_mercury_mod.o    \
                              pbl_mix_mod.o           pjc_pfix_mod.o         \
                              planeflight_mod.o       sulfate_mod.o          \
                              tagged_co_mod.o         tomas_mod.o            \
                              transport_mod.o         ucx_mod.o              \
                              uvalbedo_mod.o          tpcore_fvdas_mod.o     \
                              tpcore_window_mod.o     flexchem_mod.o         \
                              toms_mod.o              emissions_mod.o        \
                              set_global_ch4_mod.o    pops_mod.o             \
                              rrtmg_rad_transfer_mod.o

cldice_HBrHOBr_rxn.o        : cldice_HBrHOBr_rxn.F

co2_mod.o                   : co2_mod.F                                      \
                              diag04_mod.o            hco_interface_mod.o

convection_mod.o            : convection_mod.F                               \
                              diag_mod.o              depo_mercury_mod.o     \
                              wetscav_mod.o           mercury_mod.o          \
                              tendencies_mod.o        hco_interface_mod.o    \
                              diagnostics_mod.o

dao_mod.o                   : dao_mod.F

depo_mercury_mod.o          : depo_mercury_mod.F                             \
                              dao_mod.o               diag_mod.o             \
                              diag03_mod.o                     

diag03_mod.o                : diag03_mod.F

diag04_mod.o                : diag04_mod.F

diag1.o                     : diag1.F                                        \
                              hco_interface_mod.o     dao_mod.o              \
                              diag03_mod.o            diag_mod.o

diag3.o                     : diag3.F                                        \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o            diag53_mod.o           \
                              diag56_mod.o            depo_mercury_mod.o     \
                              drydep_mod.o            hco_interface_mod.o    \
                              tomas_mod.o             wetscav_mod.o

diag20_mod.o                : diag20_mod.F

diag41_mod.o                : diag41_mod.F                                   \
                              pbl_mix_mod.o

diag42_mod.o                : diag42_mod.F                                   \
                              aerosol_mod.o           carbon_mod.o            

diag48_mod.o                : diag48_mod.F                                   \
                              pbl_mix_mod.o           

diag49_mod.o                : diag49_mod.F                                   \
                              modis_lai_mod.o         pbl_mix_mod.o

diag50_mod.o                : diag50_mod.F                                   \
                              pbl_mix_mod.o           aerosol_mod.o

diag51_mod.o                : diag51_mod.F                                   \
                              pbl_mix_mod.o           

diag51b_mod.o               : diag51b_mod.F                                  \
                              pbl_mix_mod.o           

diag53_mod.o                : diag53_mod.F

diag56_mod.o                : diag56_mod.F

diag63_mod.o                : diag63_mod.F                                   \
                              diag_mod.o              pbl_mix_mod.o          \
                              hco_interface_mod.o

diag_2pm.o                  : diag_2pm.F                                     \
                              diag_mod.o

diag_mod.o                  : diag_mod.F

diag_oh_mod.o               : diag_oh_mod.F

diagnostics_mod.o           : diagnostics_mod.F90

drydep_mod.o                : drydep_mod.F                                   \
                              dao_mod.o                                      \
                              diag_mod.o              get_ndep_mod.o         \
                              pbl_mix_mod.o           tomas_mod.o            \
                              hco_interface_mod.o     diagnostics_mod.o

dust_mod.o                  : dust_mod.F                                     \
                              diag_mod.o              drydep_mod.o           \
                              hco_interface_mod.o     tomas_mod.o             

emissions_mod.o             : emissions_mod.F90                              \
                              bromocarb_mod.o         carbon_mod.o           \
                              co2_mod.o               global_ch4_mod.o       \
                              hcoi_gc_main_mod.o      mercury_mod.o          \
                              sulfate_mod.o           tomas_mod.o            \
                              ucx_mod.o               sfcvmr_mod.o           \
                              diagnostics_mod.o       pops_mod.o

exchange_mod.o              : exchange_mod.F                                 \
	                      tpcore_bc_mod.o         dao_mod.o

fast_jx_mod.o               : fast_jx_mod.F                                  \
                              dao_mod.o               toms_mod.o              

flexchem_mod.o              : flexchem_mod.F90                               \
                              aerosol_mod.o                                  \
                              diag_mod.o              diag_oh_mod.o          \
                              diag20_mod.o            dust_mod.o             \
                              fast_jx_mod.o           strat_chem_mod.o       \
                              tomas_mod.o

flexgrid_read_mod.o         : flexgrid_read_mod.F90                          \
                              dao_mod.o               diag_mod.o             \
                              get_met_mod.o

fvdas_convect_mod.o         : fvdas_convect_mod.F                            \
                              diag_mod.o              depo_mercury_mod.o

gamap_mod.o                 : gamap_mod.F                                    \
                              diag03_mod.o            diag04_mod.o           \
                              diag41_mod.o            diag42_mod.o           \
                              diag48_mod.o            diag49_mod.o           \
                              diag50_mod.o            diag51_mod.o           \
                              diag51b_mod.o           diag53_mod.o           \
                              diag56_mod.o            diag63_mod.o           \
                              drydep_mod.o            tomas_mod.o            \
                              wetscav_mod.o

gc_environment_mod.o        : gc_environment_mod.F90                         \
                              aerosol_mod.o           benchmark_mod.o        \
                              c2h6_mod.o              carbon_mod.o           \
                              co2_mod.o               diagnostics_mod.o      \
                              depo_mercury_mod.o      diag03_mod.o           \
                              diag20_mod.o            diag48_mod.o           \
                              diag49_mod.o            diag50_mod.o           \
                              diag51_mod.o            diag51b_mod.o          \
                              diag63_mod.o            diag_oh_mod.o          \
                              drydep_mod.o            dust_mod.o             \
                              gamap_mod.o             get_ndep_mod.o         \
                              global_ch4_mod.o        input_mod.o            \
                              linoz_mod.o             land_mercury_mod.o     \
                              mercury_mod.o           modis_lai_mod.o        \
                              ocean_mercury_mod.o     pops_mod.o             \
                              seasalt_mod.o           sulfate_mod.o          \
                              tagged_co_mod.o         tagged_o3_mod.o        \
                              toms_mod.o              tpcore_bc_mod.o        \
                              vdiff_pre_mod.o         wetscav_mod.o

get_met_mod.o               : get_met_mod.F90                                \
                              hco_interface_mod.o

get_ndep_mod.o              : get_ndep_mod.F                              

global_br_mod.o             : global_br_mod.F                                \
                              ocean_mercury_mod.o     hco_interface_mod.o

global_ch4_mod.o            : global_ch4_mod.F                               \
                              diag_mod.o              diag_oh_mod.o          \
                              vdiff_pre_mod.o         hco_interface_mod.o

gosat_ch4_mod.o             : gosat_ch4_mod.F

hco_interface_mod.o         : hco_interface_mod.F90

hcoi_gc_diagn_mod.o         : hcoi_gc_diagn_mod.F90                          \
                              hcoi_gc_diagn_include.H                        \
                              hco_interface_mod.o                            \
                              diag_mod.o              diag53_mod.o           \
                              diag56_mod.o 

hcoi_gc_main_mod.o          : hcoi_gc_main_mod.F90                           \
                              get_ndep_mod.o          drydep_mod.o           \
                              modis_lai_mod.o         fast_jx_mod.o          \
                              hcoi_gc_diagn_mod.o     tomas_mod.o            \
                              hco_interface_mod.o     flexgrid_read_mod.o    \
                              ocean_mercury_mod.o

initialize.o                : initialize.F                                   \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o            diag53_mod.o           \
                              diag56_mod.o

input_mod.o                 : input_mod.F                                    \
                              diag03_mod.o            diag04_mod.o           \
                              diag41_mod.o            diag42_mod.o  	     \
                              diag53_mod.o            diag56_mod.o           \
                              drydep_mod.o            planeflight_mod.o

isoropiaII_mod.o            : isoropiaII_mod.F                               \
                              hco_interface_mod.o

land_mercury_mod.o          : land_mercury_mod.F      modis_lai_mod.o        \
                              dao_mod.o               depo_mercury_mod.o     \
                              hco_interface_mod.o

linoz_mod.o                 : linoz_mod.F

main.o                      : main.F                                         \
                              benchmark_mod.o         chemistry_mod.o        \
                              carbon_mod.o            convection_mod.o       \
                              diag41_mod.o            diag42_mod.o           \
                              diag48_mod.o            diag49_mod.o           \
                              diag50_mod.o            diag51_mod.o           \
                              diag51b_mod.o           diag63_mod.o           \
                              diag_oh_mod.o           dao_mod.o              \
                              depo_mercury_mod.o      drydep_mod.o           \
                              diagnostics_mod.o       toms_mod.o             \
                              gc_environment_mod.o    global_ch4_mod.o       \
                              input_mod.o             linoz_mod.o            \
                              mapping_mod.o           toms_mod.o             \
                              mercury_mod.o           modis_lai_mod.o        \
                              ocean_mercury_mod.o     olson_landmap_mod.o    \
                              pbl_mix_mod.o           planeflight_mod.o      \
                              strat_chem_mod.o        transport_mod.o        \
                              tpcore_bc_mod.o         uvalbedo_mod.o         \
                              ucx_mod.o               vdiff_mod.o            \
                              emissions_mod.o         cleanup.o              \
                              mixing_mod.o                                   \
                              wetscav_mod.o           diag_mod.o             \
                              rrtmg_rad_transfer_mod.o                       \
                              gosat_ch4_mod.o         tccon_ch4_mod.o        \
                              aerosol_mod.o           flexgrid_read_mod.o    \
                              gc_classic_version.H

mapping_mod.o               : mapping_mod.F90                              

mercury_mod.o               : mercury_mod.F                                  \
                              dao_mod.o                                      \
                              depo_mercury_mod.o      diag03_mod.o           \
                              diag_mod.o              drydep_mod.o           \
                              global_br_mod.o         vdiff_pre_mod.o        \
                              land_mercury_mod.o      hco_interface_mod.o    \
                              ocean_mercury_mod.o     pbl_mix_mod.o

mixing_mod.o                : mixing_mod.F90                                 \
                              dao_mod.o               diagnostics_mod.o      \
                              diag_mod.o              drydep_mod.o           \
                              hcoi_gc_main_mod.o      get_ndep_mod.o         \
                              pbl_mix_mod.o           tendencies_mod.o       \
                              hco_interface_mod.o     vdiff_mod.o

modis_lai_mod.o             : modis_lai_mod.F90                              \
                              mapping_mod.o

ndxx_setup.o                : ndxx_setup.F                                   \
                              diag_mod.o                                     \
                              diag63_mod.o            diag_oh_mod.o          \
                              drydep_mod.o            planeflight_mod.o      \
                              tomas_mod.o             wetscav_mod.o

oasave.o                    : oasave.F                                       \
                              aerosol_mod.o

ocean_mercury_mod.o         : ocean_mercury_mod.F                            \
                              dao_mod.o               depo_mercury_mod.o     \
                              diag03_mod.o            toms_mod.o             \
                              hco_interface_mod.o
##############################################################################
# NOTE: For some reason gfortran 8.x.x throws an internal compiler error
# in this routine.  The error does not happen when optimization is turned
# off.  For now, lower the optimization level to get around this issue.
# The ocean mercury module might eventually be replaced later on.
#   -- Bob Yantosca, 17 Aug 2018
ifeq ($(IS_GNU_8),1)
	$(F90) -c -O1 $<
endif
##############################################################################

olson_landmap_mod.o         : olson_landmap_mod.F90                          \
                              mapping_mod.o

pbl_mix_mod.o               : pbl_mix_mod.F                                  \
                              tendencies_mod.o        diagnostics_mod.o      \
                              diag_mod.o              

pjc_pfix_window_mod.o       : pjc_pfix_window_mod.F

pjc_pfix_mod.o              : pjc_pfix_mod.F

planeflight_mod.o           : planeflight_mod.F                              \
                              diag_mod.o              ocean_mercury_mod.o     

pops_mod.o           	    : pops_mod.F    	                             \
			      diag_mod.o	      diag53_mod.o           \
                              drydep_mod.o            pbl_mix_mod.o	     \
                              vdiff_pre_mod.o         hco_interface_mod.o

RnPbBe_mod.o                : RnPbBe_mod.F                                   \
                              diag_mod.o              hco_interface_mod.o

rrtmg_rad_transfer_mod.o    : rrtmg_rad_transfer_mod.F                       \
                              set_prof_o3.o           toms_mod.o             \
                              diag_mod.o

rpmares_mod.o               : rpmares_mod.F                                  \
                              hco_interface_mod.o

seasalt_mod.o               : seasalt_mod.F           dao_mod.o              \
                              diag_mod.o              drydep_mod.o           \
                              pbl_mix_mod.o           tomas_mod.o            \
                              hco_interface_mod.o     vdiff_pre_mod.o

set_global_ch4_mod.o        : set_global_ch4_mod.F90                         \
                              pbl_mix_mod.o

set_prof_o3.o               : set_prof_o3.F                                   

sfcvmr_mod.o                : sfcvmr_mod.F90          pbl_mix_mod.o           

strat_chem_mod.o            : strat_chem_mod.F90                             \
                              linoz_mod.o             hco_interface_mod.o    \
                              dao_mod.o

sulfate_mod.o               : sulfate_mod.F                                  \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            dust_mod.o             \
                              get_ndep_mod.o                                 \
                              pbl_mix_mod.o           wetscav_mod.o          \
                              ucx_mod.o               uvalbedo_mod.o         \
                              vdiff_pre_mod.o         wetscav_mod.o          \
                              hco_interface_mod.o                            \
                              tomas_mod.o

tagged_co_mod.o             : tagged_co_mod.F                                \
                              pbl_mix_mod.o           hco_interface_mod.o    \
                              diag_mod.o

tagged_o3_mod.o             : tagged_o3_mod.F                                \
                              pbl_mix_mod.o           diag_mod.o             \
                              drydep_mod.o            hco_interface_mod.o

tccon_ch4_mod.o             : tccon_ch4_mod.F


tendencies_mod.o            : tendencies_mod.F90                             \
                              hco_interface_mod.o

toms_mod.o                  : toms_mod.F                                     \
                              hco_interface_mod.o

tpcore_bc_mod.o             : tpcore_bc_mod.F                              

tpcore_fvdas_mod.o          : tpcore_fvdas_mod.F90

tpcore_window_mod.o         : tpcore_window_mod.F90
ifeq ($(PRECISION),8)
	$(F90) -c $(FREEFORM) $(R8) $<
endif

transport_mod.o             : transport_mod.F                                \
                              diag_mod.o              dao_mod.o              \
                              pjc_pfix_mod.o          tpcore_bc_mod.o        \
                              tpcore_fvdas_mod.o      tpcore_window_mod.o    \
                              tendencies_mod.o        pjc_pfix_window_mod.o  \
                              diagnostics_mod.o

ucx_mod.o                   : ucx_mod.F                                      \
                              fast_jx_mod.o           global_ch4_mod.o       \
                              hco_interface_mod.o

uvalbedo_mod.o              : uvalbedo_mod.F90                               \
                              hco_interface_mod.o

vdiff_pre_mod.o             : vdiff_pre_mod.F90
ifeq ($(COMPILER),sun)
	$(F90) -O3 -c $<
endif

vdiff_mod.o                 : vdiff_mod.F90                                  \
                              dao_mod.o                                      \
                              depo_mercury_mod.o      diag_mod.o             \
                              drydep_mod.o            ocean_mercury_mod.o    \
                              pbl_mix_mod.o           vdiff_pre_mod.o        \
                              tendencies_mod.o        hco_interface_mod.o    \
                              global_ch4_mod.o        mercury_mod.o

wetscav_mod.o               : wetscav_mod.F                                  \
                              dao_mod.o               diag_mod.o             \
                              depo_mercury_mod.o      get_ndep_mod.o         \
                              tendencies_mod.o        hco_interface_mod.o    \
                              tomas_mod.o             diagnostics_mod.o

###############################################################################
###                                                                         ###
###  Dependencies of files specific to the TOMAS microphysics simulation    ###
###                                                                         ###
###############################################################################

tomas_mod.o                 : tomas_mod.F                                    \
                              diag_mod.o              YuIMN_Code.o           \
                              dao_mod.o                                      

tomas_tpcore_mod.o          : tomas_tpcore_mod.F90                           \
                              tomas_mod.o             


aero_drydep.o               : aero_drydep.F                                  \
                              drydep_mod.o            dust_mod.o             \
                              pbl_mix_mod.o           tomas_mod.o            \
                              diag_mod.o

YuIMN_Code.o                : YuIMN_Code.F

#EOC
