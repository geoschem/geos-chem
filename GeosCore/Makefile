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
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Define variables
ROOTDIR := ..
APM     := $(ROOTDIR)/GeosApm
BIN     := $(ROOTDIR)/bin
BPCH    := $(ROOTDIR)/GeosBpch
DOC     := $(ROOTDIR)/doc
EXE     := geos
TOM     := geostomas
HDR     := $(ROOTDIR)/Headers
HELP    := $(ROOTDIR)/help
ISO     := $(ROOTDIR)/ISOROPIA
LIB     := $(ROOTDIR)/lib
KPP     := $(ROOTDIR)/KPP
MOD     := $(ROOTDIR)/mod
NCDF    := $(ROOTDIR)/NcdfUtil
GTMM    := $(ROOTDIR)/GTMM
UTIL    := $(ROOTDIR)/GeosUtil
GIGC    := $(ROOTDIR)/GIGC

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

# List of source files
SRC := $(wildcard *.F) $(wildcard *.F90)

# Replace .f and .f90 extensions with *.o
TMP := $(SRC:.F=.o)
OBJ := $(TMP:.F90=.o)

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean realclean doc docclean help distclean tauclean

all:                                               # Build libraries
	@$(MAKE) lib                               #  and the executable
ifeq ($(UNAME),Darwin)
	ranlib -c ../lib/*
endif
	@$(MAKE) exe

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
	@$(MAKE) libheaders         		  
	@$(MAKE) libkpp             		  
	@$(MAKE) libutil 			  
	@$(MAKE) libiso 			  
	@$(MAKE) libcore			  

libcore: $(OBJ)                                    # Build code in GeosCore/

libiso:                                            # Build code in ISOROPIA/
	@$(MAKE) -C $(ISO)			  

libkpp:                                            # Build code in KPP/
	@$(MAKE) -C $(KPP)			  

libnc:                                             # Build code in NcdfUtil/
	@$(MAKE) -C $(NCDF) libnc

libgigc:
	@$(MAKE) -C $(GIGC) lib

baselibs:
	@$(MAKE) -C $(GIGC) baselibs

ncdfcheck:                                         # Check netCDF library
	@$(MAKE) libnc
	@$(MAKE) -C $(NCDF) ncdfcheck

libutil:                                           # Build code in GeosUtil/
	@$(MAKE) -C $(UTIL)			  

libheaders:                                        # Build code in Headers/ 
	@$(MAKE) -C $(HDR)			  

exe:                                               # Build executable
	$(LD) $(OBJ) $(LINK) -o $(EXE)
	cp -f $(EXE) $(BIN)

tomexe:                                               # Build executable
	$(LD) $(OBJ) $(LINK) -o $(TOM)
	cp -f $(TOM) $(BIN)

clean:                                             # Remove files here
	rm -f *.o *.mod geos geosapm

distclean:                                         # Synonym for "realclean"
	@$(MAKE) realclean

realclean:                                         # Remove files everywhere
	@$(MAKE) clean
	@$(MAKE) -C $(NCDF) clean
	@$(MAKE) -C $(APM) clean
	@$(MAKE) -C $(ISO) clean
	@$(MAKE) -C $(KPP) realclean
	@$(MAKE) -C $(GTMM) clean
	@$(MAKE) -C $(UTIL) clean
	@$(MAKE) -C $(HDR) clean
	@$(MAKE) docclean
ifeq ($(HPC),yes)
	@$(MAKE) -C $(GIGC) clean
endif
	rm -f $(LIB)/*.a
	rm -f $(MOD)/*.mod
	rm -f $(BIN)/geos*

doc:                                               # Build documentation
	@$(MAKE) -C $(DOC) doc

docclean:                                          # Remove documentation
	@$(MAKE) -C $(DOC) clean

help:                                              # Show help screen
	@$(MAKE) -C $(HELP) help

debug:
	@@echo "targets: $(MAKECMDGOALS)"

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
	$(LD) $(OBJ) $(LINK) -o $(EXE)
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

a3_read_mod.o               : a3_read_mod.F                                  \
                              diag_mod.o

a6_read_mod.o               : a6_read_mod.F                                  \
                              dao_mod.o               diag_mod.o              

acetone_mod.o               : acetone_mod.F                                  \
                              diag_mod.o              megan_mod.o 

aeic_mod.o                  : aeic_mod.F                                     \
                              chemgrid_mod.o          diag_mod.o             \
                              tracerid_mod.o

aerosol_mod.o               : aerosol_mod.F                                  \
                              chemgrid_mod.o          comode_mod.o           \
                              dao_mod.o               diag_mod.o             \
                              tracerid_mod.o          ucx_mod.o

airmas.o                    : airmas.F

anthroems.o                 : anthroems.F                                    \
                              edgar_mod.o             future_emissions_mod.o \
                              geia_mod.o              scale_anthro_mod.o     \
                              tracerid_mod.o

arctas_ship_emiss_mod.o     : arctas_ship_emiss_mod.F                        \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o

arsl1k.o                    : arsl1k.F

backsub.o                   : backsub.F

benchmark_mod.o             : benchmark_mod.F                                \
                              tracerid_mod.o

biofuel_mod.o               : biofuel_mod.F                                  \
                              bravo_mod.o             cac_anthro_mod.o       \
                              diag_mod.o              emep_mod.o             \
                              epa_nei_mod.o           future_emissions_mod.o \
                              logical_mod.o           streets_anthro_mod.o   \
                              tracerid_mod.o          tracer_mod.o

biofit.o                    : biofit.F

biomass_mod.o               : biomass_mod.F                                  \
                              diag_mod.o              gc_biomass_mod.o       \
                              gfed2_biomass_mod.o     gfed3_biomass_mod.o    \
                              logical_mod.o           tracerid_mod.o         \
                              tracer_mod.o

boxvl.o                     : boxvl.F

bravo_mod.o                 : bravo_mod.F                                    \
                              future_emissions_mod.o  scale_anthro_mod.o     \
                              tracerid_mod.o         

bromocarb_mod.o             : bromocarb_mod.F                                \
                              dao_mod.o               logical_mod.o          \
                              pbl_mix_mod.o           tracer_mod.o

c2h6_mod.o                  : c2h6_mod.F                                     \
                              biofuel_mod.o           biomass_mod.o          \
                              diag_mod.o              geia_mod.o             \
                              global_oh_mod.o         tracerid_mod.o         \
                              tracer_mod.o    

cac_anthro_mod.o            : cac_anthro_mod.F                               \
                              future_emissions_mod.o  scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o

calcrate.o                  : calcrate.F                                     \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              drydep_mod.o           \
                              pbl_mix_mod.o           planeflight_mod.o      \
                              tracerid_mod.o          cldice_HBrHOBr_rxn.o   \
                              ucx_mod.o

canopy_nox_mod.o            : canopy_nox_mod.F                               \
                              drydep_mod.o            drydep_toolbox_mod.o   \
                              modis_lai_mod.o

carbon_mod.o                : carbon_mod.F                                   \
                              biomass_mod.o           chemgrid_mod.o         \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              drydep_mod.o           \
                              future_emissions_mod.o  gfed2_biomass_mod.o    \
                              gfed3_biomass_mod.o     global_no3_mod.o       \
                              global_o3_mod.o         global_oh_mod.o        \
                              megan_mod.o             meganut_mod.o          \
                              pbl_mix_mod.o           rcp_mod.o              \
                              tomas_mod.o             tracer_mod.o           \
                              tracerid_mod.o          vdiff_pre_mod.o

ch3i_mod.o                  : ch3i_mod.F                                     \
                              biofuel_mod.o           biomass_mod.o          \
                              diag_mod.o              diag_pl_mod.o          \
                              fast_jx_mod.o           tracerid_mod.o

chemdr.o                    : chemdr.F                                       \
                              aerosol_mod.o           chemgrid_mod.o         \
                              comode_mod.o            dao_mod.o              \
                              diag_oh_mod.o           diag_pl_mod.o          \
                              dust_mod.o              future_emissions_mod.o \
                              planeflight_mod.o       restart_mod.o          \
                              tracerid_mod.o

chemgrid_mod.o              : chemgrid_mod.F                                 \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              logical_mod.o

chemistry_mod.o             : chemistry_mod.F                                \
                              acetone_mod.o           aerosol_mod.o          \
                              c2h6_mod.o              carbon_mod.o           \
                              ch3i_mod.o              comode_mod.o           \
                              drydep_mod.o            dust_mod.o             \
                              global_ch4_mod.o        h2_hd_mod.o            \
                              hcn_ch3cn_mod.o         isoropiaII_mod.o       \
                              mercury_mod.o           optdepth_mod.o         \
                              pops_mod.o              rpmares_mod.o          \
                              RnPbBe_mod.o            seasalt_mod.o          \
                              strat_chem_mod.o        sulfate_mod.o          \
                              tagged_co_mod.o         tagged_ox_mod.o        \
                              tomas_mod.o             tracerid_mod.o

cleanup.o                   : cleanup.F                                      \
                              acetone_mod.o           aerosol_mod.o          \
                              arctas_ship_emiss_mod.o aeic_mod.o             \
                              biomass_mod.o           biofuel_mod.o          \
                              bravo_mod.o             bromocarb_mod.o        \
                              c2h6_mod.o              cac_anthro_mod.o       \
                              carbon_mod.o            chemgrid_mod.o         \
                              co2_mod.o               comode_mod.o           \
                              dao_mod.o               depo_mercury_mod.o     \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag50_mod.o            diag51_mod.o           \
                              diag51b_mod.o           diag53_mod.o           \
                              diag_oh_mod.o           diag_pl_mod.o          \
                              drydep_mod.o            dust_mod.o             \
                              dust_dead_mod.o         edgar_mod.o            \
                              emep_mod.o              epa_nei_mod.o          \
                              gc_biomass_mod.o        gfed2_biomass_mod.o    \
                              gfed3_biomass_mod.o     global_bc_mod.o        \
                              global_ch4_mod.o        global_hno3_mod.o      \
                              global_no3_mod.o        global_nox_mod.o       \
                              global_o1d_mod.o        global_oc_mod.o        \
                              global_oh_mod.o         h2_hd_mod.o            \
                              hcn_ch3cn_mod.o         icoads_ship_mod.o      \
                              isoropiaII_mod.o        land_mercury_mod.o     \
                              lightning_nox_mod.o     linoz_mod.o            \
                              megan_mod.o             mercury_mod.o          \
                              modis_lai_mod.o         ocean_mercury_mod.o    \
                              pbl_mix_mod.o           pjc_pfix_mod.o         \
                              planeflight_mod.o       rcp_mod.o              \
                              retro_mod.o             seasalt_mod.o          \
                              sulfate_mod.o           strat_chem_mod.o       \
                              tagged_co_mod.o         tomas_mod.o            \
                              toms_mod.o              tracer_mod.o           \
                              transport_mod.o         ucx_mod.o              \
                              uvalbedo_mod.o          vistas_anthro_mod.o    \
                              wetscav_mod.o           tpcore_fvdas_mod.o     \
                              tpcore_geos5_window_mod.o                      \
                              tpcore_geosfp_window_mod.o                     \
                              hcoi_gc_main_mod.o

cldice_HBrHOBr_rxn.o        : cldice_HBrHOBr_rxn.F

co2_mod.o                   : co2_mod.F                                      \
                              biomass_mod.o           diag04_mod.o           \
	                      tracerid_mod.o          

comode_mod.o                : comode_mod.F

convection_mod.o            : convection_mod.F                               \
                              diag_mod.o              depo_mercury_mod.o     \
                              fvdas_convect_mod.o     wetscav_mod.o          \
                              gcap_convect_mod.o      mercury_mod.o          \
                              tracer_mod.o            tracerid_mod.o

dao_mod.o                   : dao_mod.F                                      \
                              logical_mod.o

decomp.o                    : decomp.F

depo_mercury_mod.o          : depo_mercury_mod.F                             \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracerid_mod.o         

diag03_mod.o                : diag03_mod.F                                   \
                              logical_mod.o           tracerid_mod.o

diag04_mod.o                : diag04_mod.F

diag1.o                     : diag1.F                                        \
                              chemgrid_mod.o          dao_mod.o              \
                              diag03_mod.o            diag_mod.o             \
                              tracer_mod.o            tracerid_mod.o

diag3.o                     : diag3.F                                        \
                              biofuel_mod.o           biomass_mod.o          \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o            diag53_mod.o           \
                              diag56_mod.o            diag_pl_mod.o          \
                              depo_mercury_mod.o      drydep_mod.o           \
                              logical_mod.o           tomas_mod.o            \
                              tracerid_mod.o          wetscav_mod.o

diag41_mod.o                : diag41_mod.F                                   \
                              pbl_mix_mod.o

diag42_mod.o                : diag42_mod.F                                   \
                              carbon_mod.o            logical_mod.o          \
                              tracerid_mod.o

diag48_mod.o                : diag48_mod.F                                   \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o

diag49_mod.o                : diag49_mod.F                                   \
                              modis_lai_mod.o         pbl_mix_mod.o          \
                              tracer_mod.o            tracerid_mod.o

diag50_mod.o                : diag50_mod.F                                   \
                              chemgrid_mod.o          comode_mod.o           \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracer_mod.o            tracerid_mod.o

diag51_mod.o                : diag51_mod.F                                   \
                              chemgrid_mod.o          modis_lai_mod.o        \
                              pbl_mix_mod.o           tracer_mod.o           \
                              tracerid_mod.o

diag51b_mod.o               : diag51b_mod.F                                  \
                              chemgrid_mod.o          modis_lai_mod.o        \
                              pbl_mix_mod.o           tracer_mod.o           \
                              tracerid_mod.o

diag53_mod.o                : diag53_mod.F

diag56_mod.o                : diag56_mod.F

diag63_mod.o                : diag63_mod.F                                   \
                              diag_mod.o              pbl_mix_mod.o

diag_2pm.o                  : diag_2pm.F                                     \
                              chemgrid_mod.o          diag_mod.o


diag_mod.o                  : diag_mod.F

diag_oh_mod.o               : diag_oh_mod.F                                  \
                              comode_mod.o            logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

diag_pl_mod.o               : diag_pl_mod.F                                  \
                              comode_mod.o            tracerid_mod.o         \
                              tracer_mod.o

drydep_mod.o                : drydep_mod.F                                   \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              get_ndep_mod.o         \
                              get_popsinfo_mod.o      logical_mod.o          \
                              pbl_mix_mod.o           tomas_mod.o            \
                              tracerid_mod.o          drydep_toolbox_mod.o

drydep_toolbox_mod.o        : drydep_toolbox_mod.F90

diagoh.o                    : diagoh.F                                       \
                              diag_mod.o

dust_dead_mod.o             : dust_dead_mod.F                                \
                              dao_mod.o              

dust_mod.o                  : dust_mod.F                                     \
                              comode_mod.o            diag_mod.o             \
                              drydep_mod.o            dust_dead_mod.o        \
                              tomas_mod.o             tracerid_mod.o

edgar_mod.o                 : edgar_mod.F                                    \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o

emep_mod.o                  : emep_mod.F                                     \
                              future_emissions_mod.o  scale_anthro_mod.o     \
                              tracerid_mod.o

emf_scale.o                 : emf_scale.F                                    \
                              tracerid_mod.o

emfossil.o                  : emfossil.F                                     \
                              bravo_mod.o             c2h6_mod.o             \
                              cac_anthro_mod.o        diag_mod.o             \
                              edgar_mod.o             emep_mod.o             \
                              epa_nei_mod.o           icoads_ship_mod.o      \
                              logical_mod.o           nei2005_anthro_mod.o   \
                              rcp_mod.o               retro_mod.o            \
                              streets_anthro_mod.o    tracerid_mod.o         \
                              vistas_anthro_mod.o

emissdr.o                   : emissdr.F                                      \
                              acetone_mod.o           aeic_mod.o             \
                              biofuel_mod.o           bromocarb_mod.o        \
                              canopy_nox_mod.o        diag_mod.o             \
                              emissions_mod.o         get_ndep_mod.o         \
			      lightning_nox_mod.o     megan_mod.o            \
                              meganut_mod.o           modis_lai_mod.o        \
                              rcp_mod.o               soil_nox_mod.o         \
                              tracerid_mod.o

emissions_mod.o             : emissions_mod.F                                \
                              arctas_ship_emiss_mod.o biomass_mod.o          \
                              bravo_mod.o             c2h6_mod.o             \
                              cac_anthro_mod.o        carbon_mod.o           \
                              ch3i_mod.o              co2_mod.o              \
                              dust_mod.o              edgar_mod.o            \
                              emep_mod.o              epa_nei_mod.o          \
                              global_ch4_mod.o        h2_hd_mod.o            \
                              hcn_ch3cn_mod.o         icoads_ship_mod.o      \
                              mercury_mod.o           nei2005_anthro_mod.o   \
                              paranox_mod.o           pops_mod.o             \
                              rcp_mod.o               retro_mod.o            \
                              RnPbBe_mod.o            seasalt_mod.o          \
                              streets_anthro_mod.o    sulfate_mod.o          \
                              tagged_co_mod.o         tracerid_mod.o         \
                              ucx_mod.o               vistas_anthro_mod.o

epa_nei_mod.o               : epa_nei_mod.F                                  \
                              future_emissions_mod.o  scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o

fast_jx_mod.o               : fast_jx_mod.F                                  \
                              chemgrid_mod.o          dao_mod.o              \
                              toms_mod.o

future_emissions_mod.o      : future_emissions_mod.F

fvdas_convect_mod.o         : fvdas_convect_mod.F                            \
                              diag_mod.o              depo_mercury_mod.o     \
                              logical_mod.o           tracerid_mod.o         \
                              tracer_mod.o

fyhoro.o                    : fyhoro.F

fcro2ho2.o                  : fcro2ho2.F

fyrno3.o                    : fyrno3.F

gamap_mod.o                 : gamap_mod.F                                    \
                              diag03_mod.o            diag04_mod.o           \
                              diag41_mod.o            diag42_mod.o           \
                              diag48_mod.o            diag49_mod.o           \
                              diag50_mod.o            diag51_mod.o           \
                              diag51b_mod.o           diag53_mod.o           \
                              diag56_mod.o            diag63_mod.o           \
                              diag_pl_mod.o  	      drydep_mod.o           \
                              get_popsinfo_mod.o      tomas_mod.o            \
                              tracerid_mod.o          tracer_mod.o           \
                              wetscav_mod.o

ifeq ($(COMPILER),sun)
	$(F90) -O0 -c $<
endif

gasconc.o                   : gasconc.F                                      \
                              chemgrid_mod.o          comode_mod.o           \
                              drydep_mod.o

gc_biomass_mod.o            : gc_biomass_mod.F                               \
                              future_emissions_mod.o  logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o    

gcap_convect_mod.o          : gcap_convect_mod.F                             \
                              diag_mod.o

gcap_read_mod.o             : gcap_read_mod.F                                \
                              diag_mod.o

geia_mod.o                  : geia_mod.F

geosfp_read_mod.o           : geosfp_read_mod.F90                            \
                              dao_mod.o               diag_mod.o

get_global_ch4.o            : get_global_ch4.F                               \
                              future_emissions_mod.o

get_ndep_mod.o              : get_ndep_mod.F                                 \
                              tracerid_mod.o

get_popsinfo_mod.o          : get_popsinfo_mod.F

gfed2_biomass_mod.o         : gfed2_biomass_mod.F                            \
                              future_emissions_mod.o  logical_mod.o          \
                              tracer_mod.o            tracerid_mod.o         

gfed3_biomass_mod.o         : gfed3_biomass_mod.F                            \
                              future_emissions_mod.o  logical_mod.o          \
                              tracer_mod.o            tracerid_mod.o         

gigc_environment_mod.o      : gigc_environment_mod.F90                       \
			      logical_mod.o           tracer_mod.o           \
                              vdiff_pre_mod.o

global_bc_mod.o             : global_bc_mod.F

global_br_mod.o             : global_br_mod.F                                \
                              chemgrid_mod.o          ocean_mercury_mod.o

global_ch4_mod.o            : global_ch4_mod.F                               \
                              biomass_mod.o           diag_mod.o             \
                              diag_oh_mod.o           global_oh_mod.o        \
                              tracer_mod.o            vdiff_pre_mod.o

global_hno3_mod.o           : global_hno3_mod.F                              \
                              tracer_mod.o

global_no3_mod.o            : global_no3_mod.F

global_nox_mod.o            : global_nox_mod.F

global_o1d_mod.o            : global_o1d_mod.F

global_o3_mod.o             : global_o3_mod.F

global_oc_mod.o             : global_oc_mod.F

global_oh_mod.o             : global_oh_mod.F                                \
                              dao_mod.o

h2_hd_mod.o                 : h2_hd_mod.F                                    \
                              biofuel_mod.o           biomass_mod.o          \
                              chemgrid_mod.o          diag_mod.o             \
                              drydep_mod.o            geia_mod.o             \
                              global_nox_mod.o        global_o1d_mod.o       \
                              global_oh_mod.o         megan_mod.o            \
                              meganut_mod.o           scale_anthro_mod.o     \
                              tagged_co_mod.o         tracerid_mod.o

hcn_ch3cn_mod.o             : hcn_ch3cn_mod.F                                \
                              biomass_mod.o           diag_mod.o             \
                              geia_mod.o              global_oh_mod.o        \
                              pbl_mix_mod.o           tracerid_mod.o

hco_arr_mod.o               : hco_arr_mod.F90                                \
                              hco_error_mod.o

hco_calc_mod.o              : hco_calc_mod.F90                               \
                              hco_arr_mod.o           hco_datacont_mod.o     \
                              hco_diagn_mod.o         hco_emislist_mod.o     \
                              hco_error_mod.o         hco_filedata_mod.o     \
                              hco_state_mod.o         hco_tidx_mod.o

hco_chartools_mod.o         : hco_chartools_mod.F90                          \
	                      hco_error_mod.o

hco_clock_mod.o             : hco_clock_mod.F90                              \
                              hco_error_mod.o         hco_state_mod.o

hco_config_mod.o            : hco_config_mod.F90                             \
                              hco_chartools_mod.o     hco_datacont_mod.o     \
                              hco_diagn_mod.o         hco_error_mod.o        \
                              hco_filedata_mod.o      hco_readlist_mod.o     \
                              hco_state_mod.o         hco_tidx_mod.o         \
                              hcox_extlist_mod.o

hco_datacont_mod.o          : hco_datacont_mod.F90                           \
                              hco_arr_mod.o           hco_error_mod.o        \
                              hco_filedata_mod.o      hco_state_mod.o

hco_emislist_mod.o          : hco_emislist_mod.F90                           \
                              hco_datacont_mod.o      hco_error_mod.o        \
                              hco_filedata_mod.o      hco_state_mod.o        \
                              hco_tidx_mod.o

hco_error_mod.o             : hco_error_mod.F90

hco_diagn_mod.o             : hco_diagn_mod.F90                              \
                              hco_arr_mod.o           hco_clock_mod.o        \
                              hco_error_mod.o         hco_state_mod.o        \
                              hco_unit_mod.o

hco_driver_mod.o            : hco_driver_mod.F90                             \
                              hco_calc_mod.o          hco_clock_mod.o        \
                              hco_config_mod.o        hco_diagn_mod.o        \
                              hco_emislist_mod.o      hco_readlist_mod.o     \
                              hco_datacont_mod.o      hco_error_mod.o        \
                              hco_state_mod.o         hco_tidx_mod.o

hco_filedata_mod.o          : hco_filedata_mod.F90                           \
                              hco_arr_mod.o           hco_error_mod.o        \
                              hco_tools_mod.o

hco_fluxarr_mod.o           : hco_fluxarr_mod.F90                            \
                              hco_arr_mod.o           hco_error_mod.o        \
                              hco_state_mod.o

hco_geotools_mod.o          : hco_geotools_mod.F90                           \
	                      hco_error_mod.o

hco_readlist_mod.o          : hco_readlist_mod.F90                           \
                              hco_clock_mod.o         hco_datacont_mod.o     \
                              hco_emislist_mod.o      hco_error_mod.o        \
                              hco_filedata_mod.o      hco_state_mod.o        \
                              hcoi_dataread_mod.o

hco_state_mod.o             : hco_state_mod.F90                              \
	                      hco_arr_mod.o           hco_error_mod.o

hco_tidx_mod.o              : hco_tidx_mod.F90                               \
                              hco_clock_mod.o         hco_datacont_mod.o     \
                              hco_error_mod.o         hco_filedata_mod.o     \
                              hco_state_mod.o

hco_unit_mod.o              : hco_unit_mod.F90                               \
                              hco_chartools_mod.o     hco_error_mod.o

hcoi_dataread_mod.o         : hcoi_dataread_mod.F90                          \
                              hco_datacont_mod.o      hco_error_mod.o        \
                              hco_filedata_mod.o      hco_state_mod.o        \
                              hco_tidx_mod.o          hco_unit_mod.o

hcoi_gc_main_mod.o          : hcoi_gc_main_mod.F90                           \
                              hco_arr_mod.o           hco_config_mod.o       \
                              hco_clock_mod.o         hco_chartools_mod.o    \
                              hco_driver_mod.o        hco_state_mod.o        \
                              hco_error_mod.o         hco_fluxarr_mod.o      \
                              hcoi_gc_diagn_mod.o     hco_diagn_mod.o        \
			      hcox_driver_mod.o       hcox_state_mod.o       \
                              drydep_mod.o            modis_lai_mod.o        \
                              tracerid_mod.o          bromocarb_mod.o        \
                              global_oh_mod.o         henry_mod.o

hcoi_gc_diagn_mod.o         : hcoi_gc_diagn_mod.F90                          \
                              hco_diagn_mod.o         hco_error_mod.o        \
                              hco_state_mod.o         hco_arr_mod.o

hcox_custom_mod.o           : hcox_custom_mod.F90                            \
	                      hco_diagn_mod.o         hco_fluxarr_mod.o      \
                              hco_error_mod.o         hco_geotools_mod.o     \
                              hco_state_mod.o         hcox_extlist_mod.o     \
                              hcox_state_mod.o

hcox_driver_mod.o           : hcox_driver_mod.F90                            \
                              hco_error_mod.o         hco_state_mod.o        \
                              hcox_state_mod.o        hcox_custom_mod.o      \
                              hcox_dustdead_mod.o     hcox_dustginoux_mod.o  \
                              hcox_extlist_mod.o      hcox_gfed3_mod.o       \
                              hcox_lightnox_mod.o     hcox_paranox_mod.o     \
			      hcox_seaflux_mod.o      hcox_seasalt_mod.o     \
                              hcox_soilnox_mod.o      hcox_megan_mod.o

hcox_dustdead_mod.o         : hcox_dustdead_mod.F                            \
                              hco_clock_mod.o         hco_diagn_mod.o        \
                              hco_emislist_mod.o      hco_error_mod.o        \
                              hco_fluxarr_mod.o       hco_geotools_mod.o     \
                              hcox_state_mod.o        hcox_extlist_mod.o

hcox_dustginoux_mod.o       : hcox_dustginoux_mod.F90                        \
                              hco_diagn_mod.o         hco_emislist_mod.o     \
                              hco_error_mod.o         hco_fluxarr_mod.o      \
                              hco_state_mod.o         hcox_extlist_mod.o     \
                              hcox_state_mod.o

hcox_extlist_mod.o          : hcox_extlist_mod.F90                           \
                              hco_error_mod.o         hco_state_mod.o

hcox_gfed3_mod.o            : hcox_gfed3_mod.F90                             \
                              hco_diagn_mod.o         hco_emislist_mod.o     \
                              hco_error_mod.o         hco_fluxarr_mod.o      \
                              hco_state_mod.o         hcox_extlist_mod.o     \
                              hcox_state_mod.o

hcox_lightnox_mod.o         : hcox_lightnox_mod.F90                          \
                              hco_clock_mod.o         hco_diagn_mod.o        \
                              hco_driver_mod.o        hco_error_mod.o        \
                              hco_fluxarr_mod.o       hco_geotools_mod.o     \
                              hco_state_mod.o         hcox_extlist_mod.o     \
                              hcox_state_mod.o

hcox_megan_mod.o            : hcox_megan_mod.F                               \
                              hco_clock_mod.o         hco_diagn_mod.o        \
                              hco_emislist_mod.o      hco_error_mod.o        \
                              hco_fluxarr_mod.o       hco_state_mod.o        \
                              hcox_extlist_mod.o      hcox_state_mod.o

hcox_paranox_mod.o          : hcox_paranox_mod.F90                           \
                              hco_clock_mod.o         hco_diagn_mod.o        \
                              hco_driver_mod.o        hco_error_mod.o        \
                              hco_fluxarr_mod.o       hco_state_mod.o        \
                              hcox_extlist_mod.o      hcox_state_mod.o       \
                              paranox_util_mod.o      pbl_mix_mod.o          \
                              fast_jx_mod.o

hcox_seaflux_mod.o          : hcox_seaflux_mod.F90                           \
                              hco_diagn_mod.o         hco_emislist_mod.o     \
                              hco_error_mod.o         hco_fluxarr_mod.o      \
                              hco_state_mod.o         hcox_extlist_mod.o     \
                              hcox_state_mod.o                               \
                              henry_mod.o             ocean_toolbox_mod.o

hcox_seasalt_mod.o          : hcox_seasalt_mod.F90                           \
                              hco_clock_mod.o         hco_diagn_mod.o        \
                              hco_error_mod.o         hco_fluxarr_mod.o      \
                              hco_geotools_mod.o      hco_state_mod.o        \
                              hcox_extlist_mod.o      hcox_state_mod.o

hcox_soilnox_mod.o          : hcox_soilnox_mod.F90                           \
                              hco_diagn_mod.o                                \
                              hco_emislist_mod.o      hco_error_mod.o        \
                              hco_fluxarr_mod.o       hco_state_mod.o        \
                              hcox_extlist_mod.o      hcox_state_mod.o       \
                              drydep_toolbox_mod.o

hcox_state_mod.o            : hcox_state_mod.F90                            \
                              hco_arr_mod.o           hco_error_mod.o

henry_mod.o                 : henry_mod.F                                    \
                              tracer_mod.o

i6_read_mod.o               : i6_read_mod.F                                  \
                              diag_mod.o              logical_mod.o

icoads_ship_mod.o           : icoads_ship_mod.F                              \
                              future_emissions_mod.o  scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o

initialize.o                : initialize.F                                   \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o            diag53_mod.o           \
                              diag56_mod.o            diag_pl_mod.o          \
                              logical_mod.o

input_mod.o                 : input_mod.F                                    \
                              aerosol_mod.o           benchmark_mod.o        \
                              biofuel_mod.o           biomass_mod.o          \
                              bromocarb_mod.o         carbon_mod.o           \
                              chemgrid_mod.o          depo_mercury_mod.o     \
                              diag03_mod.o            diag04_mod.o           \
                              diag41_mod.o            diag42_mod.o  	     \
                              diag48_mod.o            diag49_mod.o  	     \
                              diag50_mod.o            diag51_mod.o  	     \
                              diag51b_mod.o           diag56_mod.o  	     \
                              diag63_mod.o            diag_oh_mod.o          \
                              diag_pl_mod.o           drydep_mod.o           \
                              dust_mod.o              emissions_mod.o        \
                              future_emissions_mod.o  gamap_mod.o            \
                              get_popsinfo_mod.o      land_mercury_mod.o     \
                              linoz_mod.o             logical_mod.o          \
                              mercury_mod.o           modis_lai_mod.o        \
                              ocean_mercury_mod.o     planeflight_mod.o      \
                              pops_mod.o              rcp_mod.o              \
                              restart_mod.o           seasalt_mod.o          \
                              sulfate_mod.o           tomas_mod.o            \
                              tpcore_bc_mod.o         tracerid_mod.o         \
                              tracer_mod.o            transport_mod.o        \
                              wetscav_mod.o           tagged_co_mod.o        \
                              tagged_ox_mod.o         global_ch4_mod.o       \
                              c2h6_mod.o              hcn_ch3cn_mod.o        \
                              fvdas_convect_mod.o     henry_mod.o

isoropiaII_mod.o            : isoropiaII_mod.F                               \
                              chemgrid_mod.o          global_hno3_mod.o      \
                              tracerid_mod.o

jsparse.o                   : jsparse.F               comode_mod.o

ksparse.o                   : ksparse.F

land_mercury_mod.o          : land_mercury_mod.F                             \
                              biomass_mod.o           dao_mod.o              \
                              depo_mercury_mod.o      modis_lai_mod.o        \
                              logical_mod.o           tracerid_mod.o         

lightning_nox_mod.o         : lightning_nox_mod.F                            \
                              chemgrid_mod.o          dao_mod.o              \
                              diag_mod.o              diag56_mod.o           \
                              logical_mod.o           tracerid_mod.o

linoz_mod.o                 : linoz_mod.F                                    \
                              chemgrid_mod.o          tracerid_mod.o

logical_mod.o               : logical_mod.F

lump.o                      : lump.F                                         \
                              comode_mod.o            tracerid_mod.o

main.o                      : main.F                                         \
                              a3_read_mod.o           a6_read_mod.o          \
                              benchmark_mod.o         bromocarb_mod.o        \
                              chemgrid_mod.o          chemistry_mod.o        \
                              carbon_mod.o            convection_mod.o       \
                              comode_mod.o            diag_mod.o             \
                              diag41_mod.o            diag42_mod.o           \
                              diag48_mod.o            diag49_mod.o           \
                              diag50_mod.o            diag51_mod.o           \
                              diag51b_mod.o           diag63_mod.o           \
                              diag_oh_mod.o           dao_mod.o              \
                              depo_mercury_mod.o      drydep_mod.o           \
                              emissions_mod.o         gcap_read_mod.o        \
                              geosfp_read_mod.o       gigc_environment_mod.o \
                              global_ch4_mod.o        i6_read_mod.o          \
                              input_mod.o             lightning_nox_mod.o    \
                              linoz_mod.o             logical_mod.o          \
                              mapping_mod.o           megan_mod.o            \
                              mercury_mod.o           merra_cn_mod.o         \
                              merra_a1_mod.o          merra_a3_mod.o         \
                              merra_i6_mod.o          modis_lai_mod.o        \
                              ocean_mercury_mod.o     olson_landmap_mod.o    \
                              pbl_mix_mod.o           planeflight_mod.o      \
                              restart_mod.o           soilnox_restart_mod.o  \
                              strat_chem_mod.o        toms_mod.o             \
                              tpcore_bc_mod.o         tracer_mod.o           \
                              tracerid_mod.o          transport_mod.o        \
                              ucx_mod.o               uvalbedo_mod.o         \
                              vdiff_mod.o             wetscav_mod.o          \
                              hcoi_gc_main_mod.o

mapping_mod.o               : mapping_mod.F90                                \
                              logical_mod.o

megan_mod.o                 : megan_mod.F             a3_read_mod.o          \
                              modis_lai_mod.o         logical_mod.o          \
                              meganut_mod.o           merra_a1_mod.o         \
                              geosfp_read_mod.o

meganut_mod.o               : meganut_mod.F

mercury_mod.o               : mercury_mod.F                                  \
                              chemgrid_mod.o          dao_mod.o              \
                              depo_mercury_mod.o      diag03_mod.o           \
                              diag_mod.o              drydep_mod.o           \
                              global_br_mod.o         global_o3_mod.o        \
                              global_oh_mod.o         land_mercury_mod.o     \
                              logical_mod.o           ocean_mercury_mod.o    \
                              pbl_mix_mod.o           RnPbBe_mod.o           \
                              tracerid_mod.o          vdiff_pre_mod.o

merra_a1_mod.o              : merra_a1_mod.F                                 \
                              diag_mod.o

merra_a3_mod.o              : merra_a3_mod.F                                 \
                              dao_mod.o               diag_mod.o

merra_cn_mod.o              : merra_cn_mod.F                                 \
                              diag_mod.o

merra_i6_mod.o              : merra_i6_mod.F

modis_lai_mod.o             : modis_lai_mod.F90                              \
                              logical_mod.o           mapping_mod.o

nc_soilnox_read.o           : nc_soilnox_read.F90

ndxx_setup.o                : ndxx_setup.F                                   \
                              biofuel_mod.o           diag_mod.o             \
                              diag63_mod.o            diag_oh_mod.o          \
                              drydep_mod.o            planeflight_mod.o      \
                              tomas_mod.o             tracer_mod.o           \
                              tracerid_mod.o          wetscav_mod.o

nei2005_anthro_mod.o        : nei2005_anthro_mod.F                           \
                              future_emissions_mod.o  scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o

ocean_mercury_mod.o         : ocean_mercury_mod.F                            \
                              dao_mod.o               depo_mercury_mod.o     \
                              diag03_mod.o            logical_mod.o          \
                              tracerid_mod.o

ocean_toolbox_mod.o         : ocean_toolbox_mod.F90

ohsave.o                    : ohsave.F                                       \
                              comode_mod.o            tracerid_mod.o

olson_landmap_mod.o         : olson_landmap_mod.F90                          \
                              logical_mod.o           mapping_mod.o

optdepth_mod.o              : optdepth_mod.F          diag_mod.o  

paranox_mod.o               : paranox_mod.F

paranox_util_mod.o          : paranox_util_mod.F

partition.o                 : partition.F                                    \
                              comode_mod.o            tracerid_mod.o

pbl_mix_mod.o               : pbl_mix_mod.F                                  \
                              diag_mod.o              tracer_mod.o

pderiv.o                    : pderiv.F

physproc.o                  : physproc.F                                     \
                              comode_mod.o            chemistry_mod.o

pjc_pfix_geos5_window_mod.o : pjc_pfix_geos5_window_mod.F

pjc_pfix_geosfp_window_mod.o: pjc_pfix_geosfp_window_mod.F

pjc_pfix_mod.o              : pjc_pfix_mod.F

planeflight_mod.o           : planeflight_mod.F                              \
                              chemgrid_mod.o          comode_mod.o           \
                              diag_mod.o              ocean_mercury_mod.o    \
                              tracer_mod.o

pops_mod.o           	    : pops_mod.F    	                             \
			      diag_mod.o	      diag53_mod.o           \
                              drydep_mod.o            get_popsinfo_mod.o     \
                              global_bc_mod.o	      global_oc_mod.o        \
	        	      global_o3_mod.o         global_oh_mod.o        \
                              pbl_mix_mod.o	      tracerid_mod.o         \
                              vdiff_pre_mod.o

precipfrac.o                : precipfrac.F

rcp_mod.o                   : rcp_mod.F                                      \
                              diag_mod.o              tracerid_mod.o         \
                              tracer_mod.o

readsoil.o		    : readsoil.F

readfert.o		    : readfert.F

readclim.o		    : readclim.F

readchem.o                  : readchem.F                                     \
                              diag_pl_mod.o           drydep_mod.o

reader.o                    : reader.F

read_jv_atms_dat.o          : read_jv_atms_dat.F90

restart_mod.o               : restart_mod.F                                  \
                              comode_mod.o            logical_mod.o

retro_mod.o                 : retro_mod.F                                    \
                              future_emissions_mod.o  scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o

RnPbBe_mod.o                : RnPbBe_mod.F                                   \
                              chemgrid_mod.o          diag_mod.o             \
                              vdiff_pre_mod.o

rpmares_mod.o               : rpmares_mod.F                                  \
                              chemgrid_mod.o          global_hno3_mod.o      \
                              tracer_mod.o            tracerid_mod.o

ruralbox.o                  : ruralbox.F                                     \
                              chemgrid_mod.o          comode_mod.o  

scale_anthro_mod.o          : scale_anthro_mod.F      

seasalt_mod.o               : seasalt_mod.F                                  \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            logical_mod.o          \
                              pbl_mix_mod.o           ssa_bromine_mod.o      \
                              tomas_mod.o             tracerid_mod.o         \
                              tracer_mod.o            vdiff_pre_mod.o

setemdep.o                  : setemdep.F                                     \
                              tracerid_mod.o

setemis.o                   : setemis.F                                      \
                              biofuel_mod.o           biomass_mod.o          \
                              chemgrid_mod.o          comode_mod.o           \
                              diag_mod.o              diag63_mod.o           \
                              drydep_mod.o            edgar_mod.o            \
                              emep_mod.o              emissions_mod.o        \
                              lightning_nox_mod.o     paranox_mod.o          \
                              pbl_mix_mod.o           rcp_mod.o              \
                              tracerid_mod.o

ucx_mod.o                   : ucx_mod.F                                      \
                              chemgrid_mod.o          global_ch4_mod.o       \
                              fast_jx_mod.o           logical_mod.o          \
                              tracer_mod.o            tracerid_mod.o         \

setmodel.o                  : setmodel.F

sfcwindsqr.o                : sfcwindsqr.F

smvgear.o                   : smvgear.F                                      \
                              comode_mod.o

soil_nox_mod.o              : soil_nox_mod.F

soilnox_restart_mod.o       : soilnox_restart_mod.F

ssa_bromine_mod.o           : ssa_bromine_mod.F                              \
                              bromocarb_mod.o         chemgrid_mod.o         \
                              comode_mod.o            diag_mod.o             \
                              logical_mod.o           tracerid_mod.o

strat_chem_mod.o            : strat_chem_mod.F90                             \
                              chemgrid_mod.o          dao_mod.o              \
                              logical_mod.o           linoz_mod.o            \
                              tagged_ox_mod.o         tracer_mod.o           \
                              tracerid_mod.o

streets_anthro_mod.o        : streets_anthro_mod.F                           \
                              dao_mod.o               future_emissions_mod.o \
                              logical_mod.o           scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o            

subfun.o                    : subfun.F

sulfate_mod.o               : sulfate_mod.F                                  \
                              arctas_ship_emiss_mod.o biomass_mod.o          \
                              bravo_mod.o             cac_anthro_mod.o       \
                              chemgrid_mod.o          comode_mod.o           \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            edgar_mod.o            \
                              emep_mod.o              epa_nei_mod.o          \
                              future_emissions_mod.o  get_ndep_mod.o         \
                              gfed2_biomass_mod.o     gfed3_biomass_mod.o    \
                              global_hno3_mod.o       global_no3_mod.o       \
                              global_oh_mod.o         icoads_ship_mod.o      \
                              logical_mod.o           nei2005_anthro_mod.o   \
                              pbl_mix_mod.o           rcp_mod.o              \
                              scale_anthro_mod.o      seasalt_mod.o          \
                              streets_anthro_mod.o    tomas_mod.o            \
                              tracer_mod.o            tracerid_mod.o         \
                              ucx_mod.o               uvalbedo_mod.o         \
                              vdiff_pre_mod.o         wetscav_mod.o

sunparam.o                  : sunparam.F

tagged_co_mod.o             : tagged_co_mod.F                                \
                              biofuel_mod.o           biomass_mod.o          \
                              chemgrid_mod.o          diag_mod.o             \
                              diag_pl_mod.o           global_nox_mod.o       \
                              global_oh_mod.o         megan_mod.o            \
                              meganut_mod.o           pbl_mix_mod.o          \
                              tracer_mod.o            tracerid_mod.o

tagged_ox_mod.o             : tagged_ox_mod.F                                \
                              chemgrid_mod.o          diag_mod.o             \
                              diag_pl_mod.o           drydep_mod.o           \
                              pbl_mix_mod.o           tracerid_mod.o

tcorr.o                     : tcorr.F

toms_mod.o                  : toms_mod.F

tpcore_bc_mod.o             : tpcore_bc_mod.F                                \
                              tracer_mod.o

tpcore_geos5_window_mod.o   : tpcore_geos5_window_mod.F90
	$(F90) -c $(FREEFORM) $(R8) $<

tpcore_geosfp_window_mod.o  : tpcore_geosfp_window_mod.F90
	$(F90) -c $(FREEFORM) $(R8) $<

tracer_mod.o                : tracer_mod.F

tracerid_mod.o              : tracerid_mod.F                                 \
                              logical_mod.o           tracer_mod.o

transport_mod.o             : transport_mod.F                                \
                              diag_mod.o              pjc_pfix_mod.o         \
                              tpcore_bc_mod.o         tpcore_fvdas_mod.o     \
                              tracer_mod.o                                   \
                              tpcore_geos5_window_mod.o                      \
                              tpcore_geosfp_window_mod.o                     \
                              pjc_pfix_geos5_window_mod.o                    \
                              pjc_pfix_geosfp_window_mod.o

update.o                    : update.F

vdiff_pre_mod.o             : vdiff_pre_mod.F         tracer_mod.o
ifeq ($(COMPILER),sun)
	$(F90) -O3 -c $<
endif

vdiff_mod.o                 : vdiff_mod.F90                                  \
                              comode_mod.o            dao_mod.o              \
                              depo_mercury_mod.o      diag_mod.o             \
                              drydep_mod.o            logical_mod.o          \
                              ocean_mercury_mod.o     pbl_mix_mod.o          \
                              tracer_mod.o            tracerid_mod.o         \
                              vdiff_pre_mod.o

vistas_anthro_mod.o         : vistas_anthro_mod.F                            \
                              future_emissions_mod.o  epa_nei_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o

wetscav_mod.o               : wetscav_mod.F                                  \
                              dao_mod.o               diag_mod.o             \
                              depo_mercury_mod.o      get_ndep_mod.o         \
                              get_popsinfo_mod.o      tracerid_mod.o         \
                              tracer_mod.o            tomas_mod.o

###############################################################################
###                                                                         ###
###  Dependencies of files specific to the TOMAS microphysics simulation    ###
###                                                                         ###
###############################################################################

tomas_mod.o                 : tomas_mod.F                                    \
                              chemgrid_mod.o          dao_mod.o              \
                              diag_mod.o              diag_pl_mod.o          \
                              tracerid_mod.o          cf_nucl.o              \
                              YuIMN_Code.o            ion_nucl.o

tomas_tpcore_mod.o          : tomas_tpcore_mod.F90                           \
                              tomas_mod.o             tracerid_mod.o


aero_drydep.o               : aero_drydep.F                                  \
                              chemgrid_mod.o          diag_mod.o             \
                              drydep_mod.o            dust_mod.o             \
                              pbl_mix_mod.o           tomas_mod.o            \
                              tracerid_mod.o

ion_nucl.o                  : ion_nucl.F

YuIMN_Code.o                : YuIMN_Code.F

cf_nucl.o                   : cf_nucl.F

#EOC
