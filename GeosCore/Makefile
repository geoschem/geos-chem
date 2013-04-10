#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
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
#   make TARGET [ OPTIONAL-FLAGS ]
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
#  27 Mar 2013 - S. D. Eastham - Rolled photolysis functions together
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
ROOTDIR := ..
APM     := $(ROOTDIR)/GeosApm
BIN     := $(ROOTDIR)/bin
BPCH    := $(ROOTDIR)/GeosBpch
DOC     := $(ROOTDIR)/doc
EXE     := geos
HDR     := $(ROOTDIR)/Headers
HELP    := $(ROOTDIR)/help
ISO     := $(ROOTDIR)/ISOROPIA
LIB     := $(ROOTDIR)/lib
KPP     := $(ROOTDIR)/KPP
MOD     := $(ROOTDIR)/mod
NCDF    := $(ROOTDIR)/NcdfUtil
TOM     := $(ROOTDIR)/GeosTomas
GTMM    := $(ROOTDIR)/GTMM
UTIL    := $(ROOTDIR)/GeosUtil

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

#=============================================================================
# List of files to compile (the order is important!).  We specify these as
# a list of object files (*.o).  For each object file, the "make" utility
# will find the corresponding source code file (*.f, *.f90) and compile it. 
#=============================================================================

# List of source files
SRC := $(wildcard *.F) $(wildcard *.F90)

# Replace .f and .f90 extensions with *.o
TMP := $(SRC:.F=.o)
OBJ := $(TMP:.F90=.o)

#=============================================================================
# Makefile targets: type "make help" for a complete listing!
#=============================================================================

.PHONY: clean realclean doc docclean help distclean

all:                                               # Build libraries
	@$(MAKE) lib                               #  and the executable
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
	@$(MAKE) -C $(NCDF) lib

ncdfcheck:                                         # Check netCDF library
	@$(MAKE) libnc
	@$(MAKE) -C $(NCDF) check

libutil:                                           # Build code in GeosUtil/
	@$(MAKE) -C $(UTIL)			  

libheaders:                                        # Build code in Headers/ 
	@$(MAKE) -C $(HDR)			  

exe:                                               # Build executable
	$(LD) $(OBJ) $(LINK) -o $(EXE)
	cp -f $(EXE) $(BIN)

clean:                                             # Remove files here
	rm -f *.o *.mod geos geosapm geostomas

distclean:                                         # Synonym for "realclean"
	@$(MAKE) realclean

realclean:                                         # Remove files everywhere
	@$(MAKE) clean
	@$(MAKE) -C $(NCDF) clean
	@$(MAKE) -C $(APM) clean
	@$(MAKE) -C $(ISO) clean
	@$(MAKE) -C $(KPP) realclean
	@$(MAKE) -C $(TOM) clean
	@$(MAKE) -C $(GTMM) clean
	@$(MAKE) -C $(UTIL) clean
	@$(MAKE) -C $(HDR) clean
	@$(MAKE) docclean
	rm -f $(LIB)/*.a
	rm -f $(MOD)/*.mod
	rm -f $(BIN)/geos*

doc:                                               # Build documentation
	@$(MAKE) -C $(DOC) all

docclean:                                          # Remove documentation
	@$(MAKE) -C $(DOC) clean

help:                                              # Show help screen
	@$(MAKE) -C $(HELP)

#=============================================================================
# Targets for mercury simulation (ccc, 6/7/10)
#=============================================================================

allhg:                                             # Build Hg/GTMM executable
	@$(MAKE) libhg
	@$(MAKE) exehg

libhg:                                             # Compile HG/GTMM code
	@$(MAKE) libgtmm           		     
	@$(MAKE) lib

libgtmm:                                           # Compile GTMM model code
	@$(MAKE) -C $(GTMM)

exehg:                                             # Create Hg/GTMM exec file
	$(LD) $(OBJ) $(LHG) -o $(EXE)
	cp -f $(EXE) $(BIN)

#==============================================================================
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
#==============================================================================

a3_read_mod.o               : a3_read_mod.F           dao_mod.o              \
                              diag_mod.o              logical_mod.o

a6_read_mod.o               : a6_read_mod.F           dao_mod.o              \
                              diag_mod.o              logical_mod.o

acetone_mod.o               : acetone_mod.F           dao_mod.o              \
                              diag_mod.o              megan_mod.o 

aerosol_mod.o               : aerosol_mod.F           comode_mod.o           \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracerid_mod.o         \
                              tracer_mod.o            chemgrid_mod.o    

aeic_mod.o                  : aeic_mod.F              chemgrid_mod.o         \
                              dao_mod.o               diag_mod.o              

airmas.o                    : airmas.F

anthroems.o                 : anthroems.F                                    \
                              future_emissions_mod.o  edgar_mod.o            \
                              geia_mod.o              logical_mod.o          \
                              scale_anthro_mod.o      tracer_mod.o           \
                              tracerid_mod.o

arctas_ship_emiss_mod.o     : arctas_ship_emiss_mod.F                        \
                              logical_mod.o           scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o           

arsl1k.o                    : arsl1k.F

backsub.o                   : backsub.F

benchmark_mod.o             : benchmark_mod.F                                \
                              tracerid_mod.o          tracer_mod.o  

biofuel_mod.o               : biofuel_mod.F                                  \
                              dao_mod.o               diag_mod.o             \
                              epa_nei_mod.o           future_emissions_mod.o \
                              logical_mod.o           streets_anthro_mod.o   \
                              tracerid_mod.o          tracer_mod.o    

biofit.o                    : biofit.F

biomass_mod.o               : biomass_mod.F                                  \
                              diag_mod.o              gc_biomass_mod.o       \
                              gfed2_biomass_mod.o     gfed3_biomass_mod.o    \
                              logical_mod.o                                  \
                              tracerid_mod.o          tracer_mod.o

boxvl.o                     : boxvl.F                 dao_mod.o

bravo_mod.o                 : bravo_mod.F                                    \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         

bromocarb_mod.o             : bromocarb_mod.F                                \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracer_mod.o

c2h6_mod.o                  : c2h6_mod.F                                     \
                              biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o             \
                              geia_mod.o              global_oh_mod.o        \
                              logical_mod.o           tracerid_mod.o         \
                              tracer_mod.o    

cac_anthro_mod.o            : cac_anthro_mod.F                               \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o            

calcrate.o                  : calcrate.F                                     \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              diag63_mod.o           \
                              drydep_mod.o            emissions_mod.o        \
                              logical_mod.o           pbl_mix_mod.o          \
                              planeflight_mod.o       tracerid_mod.o         \
                              cldice_HBrHOBr_rxn.o

canopy_nox_mod.o            : canopy_nox_mod.F                               \
                              dao_mod.o               drydep_mod.o           \
			      nc_soilnox_read.o

carbon_mod.o                : carbon_mod.F                                   \
                              biomass_mod.o           comode_mod.o           \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            future_emissions_mod.o \
                              gfed2_biomass_mod.o     gfed3_biomass_mod.o    \
                              global_no3_mod.o                               \
                              global_o3_mod.o         global_oh_mod.o        \
                              logical_mod.o           megan_mod.o            \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            chemgrid_mod.o         \
                              vdiff_pre_mod.o         meganut_mod.o

ch3i_mod.o                  : ch3i_mod.F                                     \
                              biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o             \
                              diag_pl_mod.o           logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              uvalbedo_mod.o

chemdr.o                    : chemdr.F                                       \
                              aerosol_mod.o           comode_mod.o           \
                              dao_mod.o               diag_oh_mod.o          \
                              diag_pl_mod.o           future_emissions_mod.o \
                              dust_mod.o              logical_mod.o          \
                              planeflight_mod.o       restart_mod.o          \
                              tracer_mod.o            tracerid_mod.o         \
                              chemgrid_mod.o          uvalbedo_mod.o 

chemgrid_mod.o              : chemgrid_mod.F                                 \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              logical_mod.o

chemistry_mod.o             : chemistry_mod.F                                \
                              acetone_mod.o           aerosol_mod.o          \
                              c2h6_mod.o              carbon_mod.o           \
                              ch3i_mod.o              comode_mod.o           \
                              dao_mod.o               drydep_mod.o           \
                              dust_mod.o              global_ch4_mod.o       \
                              h2_hd_mod.o             hcn_ch3cn_mod.o        \
                              isoropiaII_mod.o        logical_mod.o          \
                              mercury_mod.o           optdepth_mod.o         \
                              rpmares_mod.o           RnPbBe_mod.o           \
                              seasalt_mod.o           strat_chem_mod.o       \
                              sulfate_mod.o           tagged_co_mod.o        \
                              tagged_ox_mod.o         tracerid_mod.o         \
                              tracer_mod.o            pops_mod.o

cleanup.o                   : cleanup.F                                      \
                              acetone_mod.o           aerosol_mod.o          \
                              arctas_ship_emiss_mod.o aeic_mod.o             \
                              biomass_mod.o           biofuel_mod.o          \
                              bravo_mod.o             c2h6_mod.o             \
                              cac_anthro_mod.o        carbon_mod.o           \
                              co2_mod.o               comode_mod.o           \
                              dao_mod.o               depo_mercury_mod.o     \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag50_mod.o            diag51_mod.o           \
                              diag51b_mod.o           diag_oh_mod.o          \
                              diag_pl_mod.o           drydep_mod.o           \
                              dust_mod.o              dust_dead_mod.o        \
                              edgar_mod.o             emep_mod.o             \
                              epa_nei_mod.o           gc_biomass_mod.o       \
                              isoropiaII_mod.o        gfed2_biomass_mod.o    \
                              gfed3_biomass_mod.o                            \
                              global_ch4_mod.o        global_hno3_mod.o      \
                              global_no3_mod.o        global_nox_mod.o       \
                              global_o1d_mod.o        global_oh_mod.o        \
                              gmistrat_mod.o                                 \
                              h2_hd_mod.o             hcn_ch3cn_mod.o        \
                              modis_lai_mod.o         land_mercury_mod.o     \
                              lightning_nox_mod.o     linoz_mod.o            \
                              megan_mod.o             mercury_mod.o          \
                              ocean_mercury_mod.o     pbl_mix_mod.o          \
                              pjc_pfix_mod.o          planeflight_mod.o      \
                              seasalt_mod.o           sulfate_mod.o          \
                              tagged_co_mod.o         toms_mod.o             \
                              tracer_mod.o            transport_mod.o        \
                              chemgrid_mod.o          uvalbedo_mod.o         \
                              vistas_anthro_mod.o     wetscav_mod.o          \
                              icoads_ship_mod.o       tpcore_fvdas_mod.o     \
                              tpcore_geos5_window_mod.o                      \
                              tpcore_geos57_window_mod.o                     \
                              retro_mod.o             strat_chem_mod.o

cldice_HBrHOBr_rxn.o        : cldice_HBrHOBr_rxn.F

co2_mod.o                   : co2_mod.F               biomass_mod.o          \
                              diag04_mod.o            tracer_mod.o           \
	                      tracerid_mod.o          

CO_strat_pl.o               : CO_strat_pl.F                                  \
                              dao_mod.o               tracer_mod.o           \
                              tracerid_mod.o          chemgrid_mod.o   

comode_mod.o                : comode_mod.F            tracer_mod.o

convection_mod.o            : convection_mod.F                               \
                              dao_mod.o               diag_mod.o             \
                              depo_mercury_mod.o      fvdas_convect_mod.o    \
                              gc_type_mod.o           gcap_convect_mod.o     \
                              logical_mod.o           mercury_mod.o          \
                              tracer_mod.o            tracerid_mod.o         \
                              wetscav_mod.o          

dao_mod.o                   : dao_mod.F                                      \
                              logical_mod.o           tracer_mod.o           \
			      gc_type_mod.o

decomp.o                    : decomp.F

depo_mercury_mod.o          : depo_mercury_mod.F                             \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracerid_mod.o         

diag03_mod.o                : diag03_mod.F            tracerid_mod.o

diag04_mod.o                : diag04_mod.F

diag1.o                     : diag1.F                                        \
                              dao_mod.o               diag_mod.o             \
                              tracer_mod.o            tracerid_mod.o         \
                              chemgrid_mod.o          ocean_mercury_mod.o

diag3.o                     : diag3.F                                        \
                              biofuel_mod.o           biomass_mod.o          \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o            diag56_mod.o           \
                              diag_pl_mod.o           depo_mercury_mod.o     \
                              drydep_mod.o            logical_mod.o          \
                              tracer_mod.o            tracerid_mod.o         \
                              wetscav_mod.o

diag41_mod.o                : diag41_mod.F            pbl_mix_mod.o    

diag42_mod.o                : diag42_mod.F                                   \
                              dao_mod.o               logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

diag48_mod.o                : diag48_mod.F                                   \
                              dao_mod.o               pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

diag49_mod.o                : diag49_mod.F                                   \
                              dao_mod.o               modis_lai_mod.o        \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            fast_jx_mod.o 

diag50_mod.o                : diag50_mod.F                                   \
                              comode_mod.o            dao_mod.o              \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              chemgrid_mod.o    

diag51_mod.o                : diag51_mod.F                                   \
                              dao_mod.o               modis_lai_mod.o        \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              chemgrid_mod.o    

diag51b_mod.o               : diag51b_mod.F                                  \
                              dao_mod.o               modis_lai_mod.o        \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              chemgrid_mod.o    

diag53_mod.o                : diag53_mod.F				     \
			      tracerid_mod.o               

diag56_mod.o                : diag56_mod.F

diag63_mod.o                : diag63_mod.F                                   \
                              dao_mod.o               pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o

diag_2pm.o                  : diag_2pm.F                                     \
                              diag_mod.o              chemgrid_mod.o  


diag_mod.o                  : diag_mod.F

diag_oh_mod.o               : diag_oh_mod.F                                  \
                              comode_mod.o            logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

diag_pl_mod.o               : diag_pl_mod.F                                  \
                              comode_mod.o            logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

drydep_mod.o                : drydep_mod.F            comode_mod.o          \
                              dao_mod.o               diag_mod.o            \
                              logical_mod.o           pbl_mix_mod.o         \
                              tracerid_mod.o          tracer_mod.o          \
                              meganut_mod.o           get_ndep_mod.o        \
                              get_popsinfo_mod.o

diagoh.o                    : diagoh.F                diag_mod.o

dust_dead_mod.o             : dust_dead_mod.F         dao_mod.o              

dust_mod.o                  : dust_mod.F                                     \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              drydep_mod.o           \
                              dust_dead_mod.o         logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o    

edgar_mod.o                 : edgar_mod.F                                    \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o

emep_mod.o                  : emep_mod.F                                     \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o            

emf_scale.o                 : emf_scale.F             tracerid_mod.o

emfossil.o                  : emfossil.F                                     \
                              bravo_mod.o             cac_anthro_mod.o       \
                              dao_mod.o               diag_mod.o             \
                              edgar_mod.o             emep_mod.o             \
                              epa_nei_mod.o           icoads_ship_mod.o      \
                              logical_mod.o           nei2005_anthro_mod.o   \
                              streets_anthro_mod.o    tracer_mod.o           \
                              tracerid_mod.o          vistas_anthro_mod.o    \
			      retro_mod.o             c2h6_mod.o

emissdr.o                   : emissdr.F                                      \
                              acetone_mod.o           aeic_mod.o             \
                              biofuel_mod.o           dao_mod.o              \
                              diag_mod.o              emissions_mod.o        \
			      lightning_nox_mod.o     logical_mod.o          \
                              megan_mod.o             tracer_mod.o           \
                              tracerid_mod.o          meganut_mod.o          \
                              canopy_nox_mod.o        soil_nox_mod.o         \
			      get_ndep_mod.o          bromocarb_mod.o

emissions_mod.o             : emissions_mod.F                                \
                              arctas_ship_emiss_mod.o biomass_mod.o          \
                              bravo_mod.o             c2h6_mod.o             \
                              cac_anthro_mod.o        carbon_mod.o           \
                              ch3i_mod.o              co2_mod.o              \
                              dust_mod.o              edgar_mod.o            \
                              emep_mod.o              epa_nei_mod.o          \
                              global_ch4_mod.o        h2_hd_mod.o            \
                              hcn_ch3cn_mod.o         icoads_ship_mod.o      \
                              logical_mod.o           mercury_mod.o          \
                              nei2005_anthro_mod.o    paranox_mod.o          \
                              retro_mod.o             RnPbBe_mod.o           \
                              seasalt_mod.o           streets_anthro_mod.o   \
                              sulfate_mod.o           tagged_co_mod.o        \
                              tracer_mod.o            vistas_anthro_mod.o    \
                              pops_mod.o              gmistrat_mod.o

epa_nei_mod.o               : epa_nei_mod.F                                  \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o    

fast_jx_mod.o               : fast_jx_mod.F                                  \
                              dao_mod.o               toms_mod.o

future_emissions_mod.o      : future_emissions_mod.F

fvdas_convect_mod.o         : fvdas_convect_mod.F                            \
                              dao_mod.o               diag_mod.o             \
                              depo_mercury_mod.o      logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

fyhoro.o                    : fyhoro.F

fyrno3.o                    : fyrno3.F

gamap_mod.o                 : gamap_mod.F                                    \
                              diag03_mod.o            diag04_mod.o           \
                              diag41_mod.o            diag42_mod.o           \
                              diag48_mod.o            diag49_mod.o           \
                              diag50_mod.o            diag51_mod.o           \
                              diag51b_mod.o           diag56_mod.o           \
                              diag_pl_mod.o  	      drydep_mod.o           \
                              logical_mod.o    	      tracerid_mod.o         \
                              tracer_mod.o            wetscav_mod.o          \
                              diag63_mod.o            get_popsinfo_mod.o

ifeq ($(COMPILER),sun)
	$(F90) -O0 -c $<
endif

gasconc.o                   : gasconc.F                                      \
                              comode_mod.o            dao_mod.o              \
                              drydep_mod.o            chemgrid_mod.o         \
                              logical_mod.o

gc_biomass_mod.o            : gc_biomass_mod.F                               \
                              future_emissions_mod.o  logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o    

gc_type_mod.o               : gc_type_mod.F

gcap_convect_mod.o          : gcap_convect_mod.F                             \
                              dao_mod.o               diag_mod.o

gcap_read_mod.o             : gcap_read_mod.F         dao_mod.o              \
                              diag_mod.o              logical_mod.o        

geia_mod.o                  : geia_mod.F

geos57_read_mod.o           : geos57_read_mod.F90                            \
                              dao_mod.o

get_popsinfo_mod.o          : get_popsinfo_mod.F

get_global_ch4.o            : get_global_ch4.F                               \
                              future_emissions_mod.o  logical_mod.o          

get_ndep_mod.o              : get_ndep_mod.F          tracerid_mod.o         \
			      diag_mod.o 

getifsun.o                  : getifsun.F              comode_mod.o

gfed2_biomass_mod.o         : gfed2_biomass_mod.F                            \
                              future_emissions_mod.o  logical_mod.o          \
                              tracer_mod.o            tracerid_mod.o         

gfed3_biomass_mod.o         : gfed3_biomass_mod.F                            \
                              future_emissions_mod.o  logical_mod.o          \
                              tracer_mod.o            tracerid_mod.o         

gigc_environment_mod.o      : gigc_environment_mod.F90                       \
			      gc_type_mod.o           vdiff_pre_mod.o

global_br_mod.o             : global_br_mod.F         chemgrid_mod.o         \
                              ocean_mercury_mod.o

global_ch4_mod.o            : global_ch4_mod.F                               \
                              dao_mod.o               diag_mod.o             \
                              diag_oh_mod.o           diag_pl_mod.o          \
                              global_oh_mod.o         logical_mod.o          \
                              tracer_mod.o            vdiff_pre_mod.o        \
                              biomass_mod.o  

global_hno3_mod.o           : global_hno3_mod.F                              \
                              dao_mod.o               tracer_mod.o  

global_no3_mod.o            : global_no3_mod.F

global_nox_mod.o            : global_nox_mod.F

global_o1d_mod.o            : global_o1d_mod.F

global_o3_mod.o             : global_o3_mod.F

global_oh_mod.o             : global_oh_mod.F                                \
                              dao_mod.o

global_oc_mod.o             : global_oc_mod.F

global_bc_mod.o             : global_bc_mod.F

h2_hd_mod.o                 : h2_hd_mod.F                                    \
                              biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            geia_mod.o             \
                              global_nox_mod.o        global_o1d_mod.o       \
                              global_oh_mod.o         logical_mod.o    	     \
                              scale_anthro_mod.o      tagged_co_mod.o        \
                              tracerid_mod.o          tracer_mod.o           \
                              chemgrid_mod.o          meganut_mod.o

hcn_ch3cn_mod.o             : hcn_ch3cn_mod.F                                \
                              biomass_mod.o           dao_mod.o              \
                              diag_mod.o              geia_mod.o             \
                              global_oh_mod.o         logical_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o

i6_read_mod.o               : i6_read_mod.F                                  \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           

icoads_ship_mod.o           : icoads_ship_mod.F                              \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o            

initialize.o                : initialize.F                                   \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o            diag56_mod.o           \
                              diag_pl_mod.o           logical_mod.o          \
                              diag53_mod.o

input_mod.o                 : input_mod.F                                    \
                              benchmark_mod.o         biofuel_mod.o          \
                              depo_mercury_mod.o      diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o  	      diag48_mod.o           \
                              diag49_mod.o  	      diag50_mod.o           \
			      diag51_mod.o  	      diag51b_mod.o          \
			      diag56_mod.o  	      diag_oh_mod.o          \
			      diag_pl_mod.o 	      drydep_mod.o           \
                              emissions_mod.o         future_emissions_mod.o \
                              gamap_mod.o             land_mercury_mod.o     \
                              logical_mod.o 	      mercury_mod.o          \
                              ocean_mercury_mod.o     planeflight_mod.o      \
                              restart_mod.o           tpcore_bc_mod.o        \
                              tracerid_mod.o          tracer_mod.o           \
                              transport_mod.o         wetscav_mod.o          \
                              diag63_mod.o            bromocarb_mod.o        \
                              get_popsinfo_mod.o      pops_mod.o

isoropiaII_mod.o            : isoropiaII_mod.F                               \
                              dao_mod.o               global_hno3_mod.o      \
                              logical_mod.o           tracerid_mod.o         \
                              tracer_mod.o            chemgrid_mod.o  

jsparse.o                   : jsparse.F               comode_mod.o

ksparse.o                   : ksparse.F

land_mercury_mod.o          : land_mercury_mod.F                             \
                              biomass_mod.o           dao_mod.o              \
                              depo_mercury_mod.o      modis_lai_mod.o        \
                              logical_mod.o           tracerid_mod.o         

lightning_nox_mod.o         : lightning_nox_mod.F                            \
                              dao_mod.o               diag56_mod.o           \
                              diag_mod.o              logical_mod.o        

linoz_mod.o                 : linoz_mod.F                                    \
                              dao_mod.o               tracerid_mod.o         \
                              tracer_mod.o            chemgrid_mod.o    

logical_mod.o               : logical_mod.F

lump.o                      : lump.F                                         \
                              comode_mod.o            tracerid_mod.o

main.o                      : main.F                                         \
                              a3_read_mod.o           a6_read_mod.o          \
                              benchmark_mod.o         chemistry_mod.o        \
                              convection_mod.o        comode_mod.o           \
                              diag_mod.o              diag41_mod.o           \
                              diag42_mod.o            diag48_mod.o           \
                              diag49_mod.o            diag50_mod.o           \
                              diag51_mod.o            diag51b_mod.o          \
                              diag_oh_mod.o           dao_mod.o              \
                              depo_mercury_mod.o      drydep_mod.o           \
                              emissions_mod.o         global_ch4_mod.o       \
                              gcap_read_mod.o         i6_read_mod.o          \
                              input_mod.o             modis_lai_mod.o        \
                              lightning_nox_mod.o     linoz_mod.o            \
                              logical_mod.o           megan_mod.o            \
                              pbl_mix_mod.o           ocean_mercury_mod.o    \
                              planeflight_mod.o       soaprod_mod.o          \
                              tpcore_bc_mod.o         tracer_mod.o           \
                              transport_mod.o         chemgrid_mod.o         \
                              restart_mod.o           uvalbedo_mod.o         \
                              vdiff_mod.o             wetscav_mod.o          \
                              merra_cn_mod.o          merra_a1_mod.o         \
                              merra_a3_mod.o          merra_i6_mod.o         \
                              gigc_environment_mod.o  geos57_read_mod.o      \
                              diag63_mod.o            olson_landmap_mod.o    \
                              mapping_mod.o           toms_mod.o             \
                              soilnox_restart_mod.o   gmistrat_mod.o

mapping_mod.o               : mapping_mod.F90                                \
                              logical_mod.o

megan_mod.o                 : megan_mod.F             a3_read_mod.o          \
                              modis_lai_mod.o         logical_mod.o          \
                              meganut_mod.o           merra_a1_mod.o         \
                              geos57_read_mod.o

meganut_mod.o               : meganut_mod.F           dao_mod.o     

mercury_mod.o               : mercury_mod.F                                  \
                              dao_mod.o               depo_mercury_mod.o     \
                              diag03_mod.o            diag_mod.o             \
                              drydep_mod.o            global_br_mod.o        \
                              global_o3_mod.o         global_oh_mod.o        \
                              land_mercury_mod.o                             \
                              logical_mod.o           ocean_mercury_mod.o    \
                              pbl_mix_mod.o           RnPbBe_mod.o           \
                              tracerid_mod.o          tracer_mod.o           \
                              chemgrid_mod.o          vdiff_pre_mod.o

merra_a1_mod.o              : merra_a1_mod.F                                 \
                              logical_mod.o           dao_mod.o

merra_a3_mod.o              : merra_a3_mod.F                                 \
                              logical_mod.o           dao_mod.o

merra_cn_mod.o              : merra_cn_mod.F                                 \
                              logical_mod.o           dao_mod.o

merra_i6_mod.o              : merra_i6_mod.F                                 \
                              logical_mod.o           dao_mod.o

modis_lai_mod.o             : modis_lai_mod.F90                              \
                              mapping_mod.o

nc_soilnox_read.o           : nc_soilnox_read.F90

ndxx_setup.o                : ndxx_setup.F                                   \
                              biofuel_mod.o           diag_mod.o             \
                              diag_oh_mod.o           drydep_mod.o           \
                              logical_mod.o           planeflight_mod.o      \
                              tracer_mod.o            tracerid_mod.o         \
                              wetscav_mod.o           diag63_mod.o     

nei2005_anthro_mod.o        : nei2005_anthro_mod.F                           \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o    

ocean_mercury_mod.o         : ocean_mercury_mod.F                            \
                              dao_mod.o               depo_mercury_mod.o     \
                              diag03_mod.o            logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o    

ohsave.o                    : ohsave.F                comode_mod.o           \
                              diag_mod.o              tracerid_mod.o

olson_landmap_mod.o         : olson_landmap_mod.F90                          \
                              mapping_mod.o

optdepth_mod.o              : optdepth_mod.F          diag_mod.o  

paranox_mod.o               : paranox_mod.F                                  \
                              dao_mod.o               tracer_mod.o           \
                              tracerid_mod.o

partition.o                 : partition.F                                    \
                              comode_mod.o            tracerid_mod.o

pbl_mix_mod.o               : pbl_mix_mod.F                                  \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracer_mod.o  

pderiv.o                    : pderiv.F

physproc.o                  : physproc.F                                     \
                              comode_mod.o            logical_mod.o          \
                              chemistry_mod.o

pjc_pfix_geos5_window_mod.o : pjc_pfix_geos5_window_mod.F

pjc_pfix_geos57_window_mod.o: pjc_pfix_geos57_window_mod.F

pjc_pfix_mod.o              : pjc_pfix_mod.F

planeflight_mod.o           : planeflight_mod.F                              \
                              comode_mod.o            dao_mod.o              \
                              tracer_mod.o            chemgrid_mod.o         \
                              ocean_mercury_mod.o

pops_mod.o           	    : pops_mod.F    	                             \
			      diag_mod.o	      vdiff_pre_mod.o	     \
	        	      global_oh_mod.o         pbl_mix_mod.o	     \
			      logical_mod.o	      diag53_mod.o	     \
                              global_oc_mod.o	      global_bc_mod.o	     \
			      tracer_mod.o            tracerid_mod.o         \
                              global_o3_mod.o         dao_mod.o              \
                              drydep_mod.o            get_popsinfo_mod.o

precipfrac.o                : precipfrac.F            dao_mod.o

rdsoil.o                    : rdsoil.F

readsoil.o		    : readsoil.F

readfert.o		    : readfert.F

readclim.o		    : readclim.F

readchem.o                  : readchem.F                                     \
                              diag_pl_mod.o           drydep_mod.o           \
                              logical_mod.o

reader.o                    : reader.F

read_jv_atms_dat.o          : read_jv_atms_dat.F90

restart_mod.o               : restart_mod.F                                  \
                              comode_mod.o            dao_mod.o              \
                              logical_mod.o           tracer_mod.o  

retro_mod.o                 : retro_mod.F                                    \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o

RnPbBe_mod.o                : RnPbBe_mod.F                                   \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracer_mod.o           \
                              chemgrid_mod.o    		       

rpmares_mod.o               : rpmares_mod.F                                  \
                              dao_mod.o               global_hno3_mod.o      \
                              tracerid_mod.o          tracer_mod.o           \
                              chemgrid_mod.o          

ruralbox.o                  : ruralbox.F                                     \
                              comode_mod.o            chemgrid_mod.o  

scale_anthro_mod.o          : scale_anthro_mod.F      

seasalt_mod.o               : seasalt_mod.F                                  \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            logical_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            vdiff_pre_mod.o        \
                              ssa_bromine_mod.o

setemdep.o                  : setemdep.F              drydep_mod.o           \
                              tracer_mod.o            tracerid_mod.o

setemis.o                   : setemis.F                                      \
                                                      biofuel_mod.o          \
                              biomass_mod.o           comode_mod.o           \
                              diag_mod.o              emissions_mod.o        \
                              lightning_nox_mod.o     logical_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              chemgrid_mod.o      

gmistrat_mod.o              : gmistrat_mod.F          logical_mod.o          \
                              tracer_mod.o            tracerid_mod.o         \
                              global_ch4_mod.o        chemgrid_mod.o

setmodel.o                  : setmodel.F

sfcwindsqr.o                : sfcwindsqr.F            dao_mod.o

smvgear.o                   : smvgear.F               comode_mod.o

soaprod_mod.o               : soaprod_mod.F                                  \
                              carbon_mod.o            dao_mod.o              \
                              logical_mod.o           tracer_mod.o

soil_nox_mod.o              : soil_nox_mod.F

soilnox_restart_mod.o       : soilnox_restart_mod.F   logical_mod.o

ssa_bromine_mod.o           : ssa_bromine_mod.F                              \
                              bromocarb_mod.o         comode_mod.o           \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          chemgrid_mod.o  

streets_anthro_mod.o        : streets_anthro_mod.F                           \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o            

strat_chem_mod.o            : strat_chem_mod.F90                             \
                              dao_mod.o               logical_mod.o          \
                              linoz_mod.o             tagged_ox_mod.o        \
                              tracer_mod.o            tracerid_mod.o         \
                              chemgrid_mod.o  

subfun.o                    : subfun.F

sulfate_mod.o               : sulfate_mod.F                                  \
                              arctas_ship_emiss_mod.o biomass_mod.o          \
                              bravo_mod.o             cac_anthro_mod.o       \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              drydep_mod.o           \
                              edgar_mod.o             emep_mod.o             \
                              epa_nei_mod.o           future_emissions_mod.o \
                              gfed2_biomass_mod.o     gfed3_biomass_mod.o    \
                              global_hno3_mod.o                              \
                              global_no3_mod.o        global_oh_mod.o        \
                              icoads_ship_mod.o       logical_mod.o          \
                              nei2005_anthro_mod.o    pbl_mix_mod.o          \
                              scale_anthro_mod.o      seasalt_mod.o          \
                              streets_anthro_mod.o    tracerid_mod.o         \
                              tracer_mod.o            chemgrid_mod.o         \
                              uvalbedo_mod.o          vdiff_pre_mod.o        \
                              wetscav_mod.o           get_ndep_mod.o  

sunparam.o                  : sunparam.F

tagged_co_mod.o             : tagged_co_mod.F                                \
                              biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o             \
                              diag_pl_mod.o           global_nox_mod.o       \
                              global_oh_mod.o         logical_mod.o          \
                              megan_mod.o             meganut_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            chemgrid_mod.o   

tagged_ox_mod.o             : tagged_ox_mod.F                                \
                              dao_mod.o               diag_mod.o             \
                              diag_pl_mod.o           drydep_mod.o           \
                              logical_mod.o           meganut_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            chemgrid_mod.o  

tcorr.o                     : tcorr.F

toms_mod.o                  : toms_mod.F

tpcore_bc_mod.o             : tpcore_bc_mod.F                                \
                              logical_mod.o           tracer_mod.o  

tpcore_geos5_window_mod.o   : tpcore_geos5_window_mod.F90
	$(F90) -c $(FREEFORM) $(R8) $<

tpcore_geos57_window_mod.o  : tpcore_geos57_window_mod.F90
	$(F90) -c $(FREEFORM) $(R8) $<

tpcore_mod.o                : tpcore_mod.F            dao_mod.o              \
                              diag_mod.o              global_ch4_mod.o       \
                              tracer_mod.o
	$(F90) -c $(R8) $<

tpcore_window_mod.o         : tpcore_window_mod.F     dao_mod.o              \
                              diag_mod.o              global_ch4_mod.o       \
                              tracer_mod.o
	$(F90) -c $(R8) $<

tracer_mod.o                : tracer_mod.F

tracerid_mod.o              : tracerid_mod.F                                 \
                              logical_mod.o           tracer_mod.o           \
			      gigc_environment_mod.o

transport_mod.o             : transport_mod.F                                \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           pjc_pfix_mod.o         \
                              tpcore_mod.o            tpcore_bc_mod.o        \
                              tpcore_fvdas_mod.o      tpcore_window_mod.o    \
                              tracer_mod.o                                   \
                              tpcore_geos5_window_mod.o                      \
                              tpcore_geos57_window_mod.o                     \
                              pjc_pfix_geos5_window_mod.o                    \
                              pjc_pfix_geos57_window_mod.o

tropopause.o                : tropopause.F                                   \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           chemgrid_mod.o 

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
                              logical_mod.o           scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o           

wetscav_mod.o               : wetscav_mod.F                                  \
                              dao_mod.o               diag_mod.o             \
                              depo_mercury_mod.o      logical_mod.o          \
                              mercury_mod.o           ocean_mercury_mod.o    \
                              tracerid_mod.o          tracer_mod.o           \
                              get_ndep_mod.o          get_popsinfo_mod.o

#EOC
