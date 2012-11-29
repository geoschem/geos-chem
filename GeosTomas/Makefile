#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the GeosTomas subdirectory)
#
# !DESCRIPTION: This is the main makefile for GEOS-Chem + TOMAS aerosol 
#  microphysics.  It compiles the  GEOS-Chem core source code files and 
#  into object files (*.o).  Module files (*.mod) are copied to the MOD 
#  directory.
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
#  25 Jan 2010 - R. Yantosca - Initial version for TOMAS microphysics
#  28 Jan 2010 - C. Carouge  - Modifications for ISORROPIA II
#  10 May 2010 - R. Yantosca - Add dependency for RD_AOD.f
#  26 Aug 2010 - R. Yantosca - Modifications for MERRA met fields
#  23 Sep 2010 - R. Yantosca - Removed "tropopause.f" from compile list
#  23 Sep 2010 - R. Yantosca - Updated dependencies for v9-01-01
#  06 Dec 2011 - R. Yantosca - Remove #ifdef blocks for GI model code
#  16 Feb 2012 - M. Payer    - Remove references to ESMF directory
#  05 Apr 2012 - R. Yantosca - Always assume netCDF will be used
#  05 Apr 2012 - R. Yantosca - Added c2h6_mod.o to emfossil.F dependency
#  05 Apr 2012 - R. Yantosca - Now compile the NcdfUtil routines
#  05 Apr 2012 - R. Yantosca - Remove reference to rdland.F, rdlai.F, readlai.F
#  12 Apr 2012 - R. Yantosca - Remove reference to findmon.F
#  19 Apr 2012 - R. Yantosca - Remove reference to rd_prof.F
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
ROOTDIR := ..
APM     := $(ROOTDIR)/APM
BIN     := $(ROOTDIR)/bin
BPCH    := $(ROOTDIR)/GeosBpch
CORE    := $(ROOTDIR)/GeosCore
DOC     := $(ROOTDIR)/doc
EXE     := geostomas
HDR     := $(ROOTDIR)/Headers
HELP    := $(ROOTDIR)/help
ISO     := $(ROOTDIR)/ISOROPIA
LIB     := $(ROOTDIR)/lib
KPP     := $(ROOTDIR)/KPP
MOD     := $(ROOTDIR)/mod
NCDF    := $(ROOTDIR)/NcdfUtil
UTIL    := $(ROOTDIR)/GeosUtil
GTMM    := $(ROOTDIR)/GTMM

# This directory only contains files that are different for the TOMAS
# aerosol microphysics.  For files that are the same as for the regular
# GEOS-Chem code, look in the GeosCore directory.
VPATH   := ../GeosCore

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

#=============================================================================
# List of files to compile.  We need to manually list all object files, since 
# some of these will be referenced from the ../GeosCore directory.
#=============================================================================

OBJ =                                                                        \
BLKSLV.o                 drydep_mod.o            partition.o                 \
CLDSRF.o                 dust_dead_mod.o         pbl_mix_mod.o               \
CO_strat_pl.o            dust_mod.o              pderiv.o                    \
EFOLD.o	                 edgar_mod.o             photoj.o                    \
FLINT.o                  emep_mod.o              physproc.o                  \
GAUSSP.o                 emf_scale.o             pjc_pfix_geos5_window_mod.o \
GEN.o                    emfossil.o              pjc_pfix_mod.o              \
JRATET.o                 planeflight_mod.o       geos57_read_mod.o           \
JVALUE.o                 precipfrac.o            olson_landmap_mod.o         \
LEGND0.o                 mapping_mod.o                                       \
MATIN4.o                 emissdr.o               rd_js.o                     \
MIESCT.o                 emissions_mod.o                                     \
NOABS.o                  modis_lai_mod.o                                     \
OPMIE.o                  epa_nei_mod.o                                       \
RD_TJPL.o                fast_j.o                                            \
RnPbBe_mod.o                                     SPHERE.o                    \
XSEC1D.o                 fjfunc.o                                            \
XSECO2.o                 fjx_acet_mod.o          readchem.o                  \
XSECO3.o                 future_emissions_mod.o  reader.o                    \
a3_read_mod.o            fvdas_convect_mod.o                                 \
a6_read_mod.o            fyhoro.o                restart_mod.o               \
acetone_mod.o            fyrno3.o                rpmares_mod.o               \
aerosol_mod.o            gamap_mod.o             ruralbox.o                  \
aircraft_nox_mod.o       gasconc.o               scale_anthro_mod.o          \
airmas.o                 gc_biomass_mod.o                                    \
anthroems.o              gcap_convect_mod.o      seasalt_mod.o               \
arctas_ship_emiss_mod.o  gcap_read_mod.o         set_aer.o                   \
arsl1k.o                 geia_mod.o              set_prof.o                  \
backsub.o                get_global_ch4.o        gfed3_biomass_mod.o         \
benchmark_mod.o          getifsun.o              setemdep.o                  \
biofit.o                 gfed2_biomass_mod.o     setemis.o                   \
biofuel_mod.o            global_ch4_mod.o        setmodel.o                  \
biomass_mod.o            global_hno3_mod.o       sfcwindsqr.o                \
boxvl.o                  global_no3_mod.o        smvgear.o                   \
bravo_mod.o              global_nox_mod.o        soaprod_mod.o               \
c2h6_mod.o               global_o1d_mod.o                                    \
cac_anthro_mod.o         global_o3_mod.o                                     \
calcrate.o               global_oh_mod.o                                     \
carbon_mod.o                                                                 \
ch3i_mod.o               h2_hd_mod.o                                         \
chemdr.o                 hcn_ch3cn_mod.o         streets_anthro_mod.o        \
chemistry_mod.o          i6_read_mod.o           subfun.o                    \
cleanup.o                icoads_ship_mod.o       sulfate_mod.o               \
co2_mod.o                initialize.o            sunparam.o                  \
comode_mod.o             inphot.o                tagged_co_mod.o             \
convection_mod.o         input_mod.o             tagged_ox_mod.o             \
dao_mod.o                isoropiaII_mod.o                                    \
decomp.o                 jsparse.o               toms_mod.o                  \
diag03_mod.o             jv_index.o              tpcore_bc_mod.o             \
diag04_mod.o             ksparse.o               tpcore_fvdas_mod.o          \
diag1.o                  tpcore_geos5_window_mod.o                           \
diag3.o                  lightning_nox_mod.o     tpcore_mod.o                \
diag41_mod.o             linoz_mod.o             tpcore_window_mod.o         \
diag42_mod.o             logical_mod.o           tracer_mod.o                \
diag48_mod.o             lump.o                  tracerid_mod.o              \
diag49_mod.o             main.o                  transport_mod.o             \
diag50_mod.o             megan_mod.o             retro_mod.o                 \
diag51_mod.o             meganut_mod.o           tropopause_mod.o            \
diag51b_mod.o            mercury_mod.o                                       \
diag56_mod.o             mmran_16.o              update.o                    \
diag_2pm.o               ndxx_setup.o            uvalbedo_mod.o              \
diag_mod.o               nei2005_anthro_mod.o    vdiff_mod.o                 \
diag_oh_mod.o            ocean_mercury_mod.o     vdiff_pre_mod.o             \
diag_pl_mod.o            ohsave.o                vistas_anthro_mod.o         \
diagoh.o                 optdepth_mod.o          wetscav_mod.o               \
aero_drydep.o            tomas_mod.o                                         \
tomas_tpcore_mod.o       RD_AOD.o                depo_mercury_mod.o          \
land_mercury_mod.o       global_br_mod.o         merra_a1_mod.o              \
merra_a3_mod.o           merra_cn_mod.o          merra_i6_mod.o              \
diag63_mod.o             paranox_mod.o           read_jv_atms_dat.o          \
bromocarb_mod.o          cldice_HBrHOBr_rxn.o    ssa_bromine_mod.o           \
strat_chem_mod.o         gigc_environment_mod.o                              \
canopy_nox_mod.o         get_ndep_mod.o          nc_soilnox_read.o           \
readsoil.o               readfert.o              readclim.o                  \
soil_nox_mod.o           soilnox_restart_mod.o                               \
tpcore_geos57_window_mod.o                                                   \
pjc_pfix_geos57_window_mod.o

#=============================================================================
# Makefile targets: type "make help" for a complete listing!
#=============================================================================

.PHONY: clean realclean doc docclean help

all:
	@$(MAKE) TOMAS=yes lib
	@$(MAKE) TOMAS=yes exe

lib:                                               # Build normal GEOS-Chem
	@$(MAKE) libnc				     
	@$(MAKE) libheaders         		     
	@$(MAKE) libkpp             		     
	@$(MAKE) libutil 			     
	@$(MAKE) libiso 			     
	@$(MAKE) libtomas

libtomas: $(OBJ)                                   # Build code in GeosTomas/

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

realclean:                                         # Remove files everywhere
	@$(MAKE) clean
	@$(MAKE) -C $(ISO) clean
	@$(MAKE) -C $(KPP) realclean
	@$(MAKE) -C $(UTIL) clean
	@$(MAKE) -C $(CORE) clean
	@$(MAKE) -C $(GTMM) clean
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
                              tracer_mod.o            tropopause_mod.o  

aircraft_nox_mod.o          : aircraft_nox_mod.F                             \
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

BLKSLV.o                    : BLKSLV.F

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
                              tracer_mod.o            tropopause_mod.o       \
                              vdiff_pre_mod.o         meganut_mod.o          \
                              tomas_mod.o             

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
                              tropopause_mod.o        uvalbedo_mod.o 

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
                              tracer_mod.o            tomas_mod.o

cleanup.o                   : cleanup.F                                      \
                              acetone_mod.o           aerosol_mod.o          \
                              arctas_ship_emiss_mod.o aircraft_nox_mod.o     \
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
                              h2_hd_mod.o             hcn_ch3cn_mod.o        \
                              modis_lai_mod.o         land_mercury_mod.o     \
                              lightning_nox_mod.o     linoz_mod.o            \
                              megan_mod.o             mercury_mod.o          \
                              ocean_mercury_mod.o     pbl_mix_mod.o          \
                              pjc_pfix_mod.o          planeflight_mod.o      \
                              seasalt_mod.o           sulfate_mod.o          \
                              tagged_co_mod.o         toms_mod.o             \
                              tracer_mod.o            transport_mod.o        \
                              tropopause_mod.o        uvalbedo_mod.o         \
                              vistas_anthro_mod.o     wetscav_mod.o          \
                              icoads_ship_mod.o       tpcore_fvdas_mod.o     \
                              tpcore_geos5_window_mod.o                      \
                              tpcore_geos57_window_mod.o                     \
                              retro_mod.o             strat_chem_mod.o       \
                              tomas_mod.o

cldice_HBrHOBr_rxn.o        : cldice_HBrHOBr_rxn.F

CLDSRF.o                    : CLDSRF.F

co2_mod.o                   : co2_mod.F               biomass_mod.o          \
                              diag04_mod.o            tracer_mod.o           \
	                      tracerid_mod.o          

CO_strat_pl.o               : CO_strat_pl.F                                  \
                              dao_mod.o               tracer_mod.o           \
                              tracerid_mod.o          tropopause_mod.o 

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
                              tropopause_mod.o        ocean_mercury_mod.o

diag3.o                     : diag3.F                                        \
                              biofuel_mod.o           biomass_mod.o          \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o            diag56_mod.o           \
                              diag_pl_mod.o           depo_mercury_mod.o     \
                              drydep_mod.o            logical_mod.o          \
                              tracer_mod.o            tracerid_mod.o         \
                              wetscav_mod.o           tomas_mod.o

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
                              tracer_mod.o  

diag50_mod.o                : diag50_mod.F                                   \
                              comode_mod.o            dao_mod.o              \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o  

diag51_mod.o                : diag51_mod.F                                   \
                              dao_mod.o               modis_lai_mod.o        \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o  

diag51b_mod.o               : diag51b_mod.F                                  \
                              dao_mod.o               modis_lai_mod.o        \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o  

diag56_mod.o                : diag56_mod.F

diag63_mod.o                : diag63_mod.F                                   \
                              dao_mod.o               pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o

diag_2pm.o                  : diag_2pm.F                                     \
                              diag_mod.o              tropopause_mod.o


diag_mod.o                  : diag_mod.F

diag_oh_mod.o               : diag_oh_mod.F                                  \
                              comode_mod.o            logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

diag_pl_mod.o               : diag_pl_mod.F                                  \
                              comode_mod.o            logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

drydep_mod.o                : drydep_mod.F            comode_mod.o           \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              meganut_mod.o           get_ndep_mod.o         \
                              tomas_mod.o 

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

EFOLD.o                     : EFOLD.F

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
                              acetone_mod.o           aircraft_nox_mod.o     \
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
                              tracer_mod.o            vistas_anthro_mod.o

epa_nei_mod.o               : epa_nei_mod.F                                  \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o    

fast_j.o                    : fast_j.F                                       \
                              dao_mod.o               toms_mod.o

fjfunc.o                    : fjfunc.F

fjx_acet_mod.o              : fjx_acet_mod.F

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
                              tomas_mod.o             diag63_mod.o 
ifeq ($(COMPILER),sun)
	$(F90) -O0 -c $<
endif

gasconc.o                   : gasconc.F                                      \
                              comode_mod.o            dao_mod.o              \
                              drydep_mod.o            tropopause_mod.o       \
                              logical_mod.o

GAUSSP.o                    : GAUSSP.F

gc_biomass_mod.o            : gc_biomass_mod.F                               \
                              future_emissions_mod.o  logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o    

gc_type_mod.o               : gc_type_mod.F

gcap_convect_mod.o          : gcap_convect_mod.F                             \
                              dao_mod.o               diag_mod.o

gcap_read_mod.o             : gcap_read_mod.F         dao_mod.o              \
                              diag_mod.o              logical_mod.o        

geia_mod.o                  : geia_mod.F

GEN.o                       : GEN.F

geos57_read_mod.o           : geos57_read_mod.F90                            \
                              dao_mod.o

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

global_br_mod.o             : global_br_mod.F         tropopause_mod.o 

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

h2_hd_mod.o                 : h2_hd_mod.F                                    \
                              biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            geia_mod.o             \
                              global_nox_mod.o        global_o1d_mod.o       \
                              global_oh_mod.o         logical_mod.o    	     \
                              scale_anthro_mod.o      tagged_co_mod.o        \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o        meganut_mod.o

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
                              diag_pl_mod.o           logical_mod.o 

inphot.o                    : inphot.F

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
                              tomas_mod.o             diag63_mod.o           \
                              bromocarb_mod.o

isoropiaII_mod.o            : isoropiaII_mod.F                               \
                              dao_mod.o               global_hno3_mod.o      \
                              logical_mod.o           tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o

jsparse.o                   : jsparse.F               comode_mod.o

JRATET.o                    : JRATET.F                fjx_acet_mod.o

JVALUE.o                    : JVALUE.F

jv_index.o                  : jv_index.F

ksparse.o                   : ksparse.F

land_mercury_mod.o          : land_mercury_mod.F                             \
                              biomass_mod.o           dao_mod.o              \
                              depo_mercury_mod.o      modis_lai_mod.o        \
                              logical_mod.o           tracerid_mod.o         \

LEGND0.o                    : LEGND0.F

lightning_nox_mod.o         : lightning_nox_mod.F                            \
                              dao_mod.o               diag56_mod.o           \
                              diag_mod.o              logical_mod.o        

linoz_mod.o                 : linoz_mod.F                                    \
                              dao_mod.o               tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o  

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
                              transport_mod.o         tropopause_mod.o       \
                              restart_mod.o           uvalbedo_mod.o         \
                              vdiff_mod.o             wetscav_mod.o          \
                              merra_cn_mod.o          merra_a1_mod.o         \
                              merra_a3_mod.o          merra_i6_mod.o         \
                              gigc_environment_mod.o  geos57_read_mod.o      \
                              diag63_mod.o            olson_landmap_mod.o    \
                              mapping_mod.o           toms_mod.o             \
                              soilnox_restart_mod.o 

mapping_mod.o               : mapping_mod.F90                                \
                              logical_mod.o

MATIN4.o                    : MATIN4.F

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
                              tropopause_mod.o        vdiff_pre_mod.o

merra_a1_mod.o              : merra_a1_mod.F                                 \
                              logical_mod.o           dao_mod.o

merra_a3_mod.o              : merra_a3_mod.F                                 \
                              logical_mod.o           dao_mod.o

merra_cn_mod.o              : merra_cn_mod.F                                 \
                              logical_mod.o           dao_mod.o

merra_i6_mod.o              : merra_i6_mod.F                                 \
                              logical_mod.o           dao_mod.o

MIESCT.o                    : MIESCT.F

modis_lai_mod.o             : modis_lai_mod.F90                              \
                              mapping_mod.o

nc_soilnox_read.o           : nc_soilnox_read.F90

ndxx_setup.o                : ndxx_setup.F                                   \
                              biofuel_mod.o           diag_mod.o             \
                              diag_oh_mod.o           drydep_mod.o           \
                              logical_mod.o           planeflight_mod.o      \
                              tracer_mod.o            tracerid_mod.o         \
                              wetscav_mod.o           tomas_mod.o            \
                              diag63_mod.o

nei2005_anthro_mod.o        : nei2005_anthro_mod.F                           \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o    

NOABS.o                     : NOABS.F

ocean_mercury_mod.o         : ocean_mercury_mod.F                            \
                              dao_mod.o               depo_mercury_mod.o     \
                              diag03_mod.o            logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o    

ohsave.o                    : ohsave.F                comode_mod.o           \
                              diag_mod.o              tracerid_mod.o

olson_landmap_mod.o         : olson_landmap_mod.F90                          \
                              mapping_mod.o

optdepth_mod.o              : optdepth_mod.F          diag_mod.o  

OPMIE.o                     : OPMIE.F

paranox_mod.o               : paranox_mod.F                                  \
                              dao_mod.o               tracer_mod.o           \
                              tracerid_mod.o

partition.o                 : partition.F                                    \
                              comode_mod.o            tracerid_mod.o

pbl_mix_mod.o               : pbl_mix_mod.F                                  \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracer_mod.o  

pderiv.o                    : pderiv.F

photoj.o                    : photoj.F

physproc.o                  : physproc.F                                     \
                              comode_mod.o            logical_mod.o          \
                              chemistry_mod.o

pjc_pfix_geos5_window_mod.o : pjc_pfix_geos5_window_mod.F

pjc_pfix_geos57_window_mod.o: pjc_pfix_geos57_window_mod.F

pjc_pfix_mod.o              : pjc_pfix_mod.F

planeflight_mod.o           : planeflight_mod.F                              \
                              comode_mod.o            dao_mod.o              \
                              tracer_mod.o            tropopause_mod.o       

precipfrac.o                : precipfrac.F            dao_mod.o

RD_AOD.o                    : RD_AOD.F

rd_js.o                     : rd_js.F

readsoil.o		    : readsoil.F

readfert.o		    : readfert.F

readclim.o		    : readclim.F

RD_TJPL.o                   : RD_TJPL.F

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
                              tropopause_mod.o  		       

rpmares_mod.o               : rpmares_mod.F                                  \
                              dao_mod.o               global_hno3_mod.o      \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o        

ruralbox.o                  : ruralbox.F                                     \
                              comode_mod.o            tropopause_mod.o

scale_anthro_mod.o          : scale_anthro_mod.F      

seasalt_mod.o               : seasalt_mod.F                                  \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            logical_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            vdiff_pre_mod.o        \
                              tomas_mod.o             ssa_bromine_mod.o

set_aer.o                   : set_aer.F

set_prof.o                  : set_prof.F                                     \
                              dao_mod.o               toms_mod.o

setemdep.o                  : setemdep.F              drydep_mod.o           \
                              tracer_mod.o            tracerid_mod.o

setemis.o                   : setemis.F                                      \
                              aircraft_nox_mod.o      biofuel_mod.o          \
                              biomass_mod.o           comode_mod.o           \
                              diag_mod.o              emissions_mod.o        \
                              lightning_nox_mod.o     logical_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tropopause_mod.o    

setmodel.o                  : setmodel.F

sfcwindsqr.o                : sfcwindsqr.F            dao_mod.o

smvgear.o                   : smvgear.F               comode_mod.o

soaprod_mod.o               : soaprod_mod.F                                  \
                              carbon_mod.o            dao_mod.o              \
                              logical_mod.o           tracer_mod.o

soil_nox_mod.o              : soil_nox_mod.F

soilnox_restart_mod.o       : soilnox_restart_mod.F   logical_mod.o

SPHERE.o                    : SPHERE.F

ssa_bromine_mod.o           : ssa_bromine_mod.F                              \
                              bromocarb_mod.o         comode_mod.o           \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          tropopause_mod.o

streets_anthro_mod.o        : streets_anthro_mod.F                           \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o            

strat_chem_mod.o            : strat_chem_mod.F90                             \
                              dao_mod.o               logical_mod.o          \
                              linoz_mod.o             tagged_ox_mod.o        \
                              tracer_mod.o            tracerid_mod.o         \
                              tropopause_mod.o

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
                              icoads_ship_mod.o       tomas_mod.o            \
                              logical_mod.o           nei2005_anthro_mod.o   \
                              pbl_mix_mod.o           scale_anthro_mod.o     \
                              seasalt_mod.o           streets_anthro_mod.o   \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o        uvalbedo_mod.o         \
                              vdiff_pre_mod.o         wetscav_mod.o          \
                              get_ndep_mod.o

sunparam.o                  : sunparam.F

tagged_co_mod.o             : tagged_co_mod.F                                \
                              biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o             \
                              diag_pl_mod.o           global_nox_mod.o       \
                              global_oh_mod.o         logical_mod.o          \
                              megan_mod.o             meganut_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o 

tagged_ox_mod.o             : tagged_ox_mod.F                                \
                              dao_mod.o               diag_mod.o             \
                              diag_pl_mod.o           drydep_mod.o           \
                              logical_mod.o           meganut_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o

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
                              tracer_mod.o            tomas_tpcore_mod.o
	$(F90) -c $(R8) $<

tpcore_window_mod.o         : tpcore_window_mod.F     dao_mod.o              \
                              diag_mod.o              global_ch4_mod.o       \
                              tracer_mod.o            tomas_mod.o
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
                              logical_mod.o           tropopause_mod.o 

tropopause_mod.o            : tropopause_mod.F                               \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              logical_mod.o

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
                              get_ndep_mod.o         tomas_mod.o

XSEC1D.o                    : XSEC1D.F

XSECO2.o                    : XSECO2.F

XSECO3.o                    : XSECO3.F

#------------------------------------------------------------------------------
# Dependencies of files specific to TOMAS microphysics (bmy, 1/25/10)
#------------------------------------------------------------------------------

tomas_mod.o                 : tomas_mod.F                                    \
	                      dao_mod.o               diag_mod.o             \
                              diag_pl_mod.o           logical_mod.o          \
                              tropopause_mod.o        tracer_mod.o           \
                              tracerid_mod.o

tomas_tpcore_mod.o          : tomas_tpcore_mod.F90                           \
                              tomas_mod.o             tracerid_mod.o


aero_drydep.o               : aero_drydep.F                                  \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            dust_mod.o             \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracer_mod.o            tracerid_mod.o         \
                              tropopause_mod.o        tomas_mod.o

#EOC
