#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the GeosApm/ subdirectory)
#
# !DESCRIPTION: This makefile compiles the APM microphysics core routines.
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
#
# !REVISION HISTORY: 
#  13 Jan 2011 - R. Yantosca - Initial version
#  22 Aug 2011 - R. Yantosca - Modified for v9-01-02
#  05 Aug 2011 - M. Long     - Now compile module files in Headers/ directory
#  27 Sep 2011 - M. Payer    - Modified dependencies for GFED3
#  06 Dec 2011 - R. Yantosca - Remove #ifdef blocks for GI model code
#  16 Feb 2012 - M. Payer    - Remove references to ESMF directory
#  05 Apr 2012 - R. Yantosca - Now also build code in NcdfUtil
#  05 Apr 2012 - R. Yantosca - Remove reference to rdland.F
#  05 Apr 2012 - R. Yantosca - Add c2h6_mod.o to emfossil.F dependency
#  05 Apr 2012 - R. Yantosca - Remove reference to rdlai.F
#  11 Apr 2012 - R. Yantosca - Remove reference to obsolete lai_mod.F
#  12 Apr 2012 - R. Yantosca - Remove reference to findmon.F
#  19 Apr 2012 - R. Yantosca - Remove reference to rd_prof.F
#  11 Jul 2012 - R. Yantosca - Remove references to schem.F, upbdflx_mod.F
#  30 Oct 2012 - R. Yantosca - Remove reference to rdsoil.F
#  27 Nov 2012 - R. Yantosca - Remove reference to getifsun.F
#  30 May 2013 - R. Yantosca - Remove reference to GEOS-3 TPCORE routines 
#  13 Aug 2013 - M. Sulprizio- Remove reference to soaprod_mod.F (no longer
#                              needed) (SOAupdate, hotp 7/25/10)
#  22 Aug 2013 - R. Yantosca - Bug fix: when making realclean, also do
#                              a "make clean" in the Headers directory
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
SHELL   := /bin/sh
ROOTDIR := ..
APM     := $(ROOTDIR)/APM
BIN     := $(ROOTDIR)/bin
BPCH    := $(ROOTDIR)/GeosBpch
CORE    := $(ROOTDIR)/GeosCore
DOC     := $(ROOTDIR)/doc
EXE     := geosapm
GTMM    := $(ROOTDIR)/GTMM
HDR     := $(ROOTDIR)/Headers
HELP    := $(ROOTDIR)/help
ISO     := $(ROOTDIR)/ISOROPIA
LIB     := $(ROOTDIR)/lib
KPP     := $(ROOTDIR)/KPP
MOD     := $(ROOTDIR)/mod
NCDF    := $(ROOTDIR)/NcdfUtil
UTIL    := $(ROOTDIR)/GeosUtil

# This directory only contains files that are different for the APM
# aerosol microphysics.  For files that are the same as for the regular
# GEOS-Chem code, look in the GeosCore directory.
VPATH   := ../GeosCore

# Include header file.  This returns variables CC, F90, FREEFORM, LD, R8,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

#=============================================================================
# List of files to compile (the order is important!).  We specify these as
# a list of object files (*.o).  For each object file, the "make" utility
# will find the corresponding source code file (*.f, *.f90) and compile it. 
#=============================================================================

OBJ =                                                                        \
BLKSLV.o                 CLDSRF.o                                            \
EFOLD.o                  FLINT.o                 GAUSSP.o                    \
GEN.o                    JRATET.o                JVALUE.o                    \
LEGND0.o                 MATIN4.o                MIESCT.o                    \
NOABS.o                  OPMIE.o                 RD_AOD.o                    \
RD_TJPL.o                RnPbBe_mod.o            SPHERE.o                    \
XSEC1D.o                 XSECO2.o                XSECO3.o                    \
a3_read_mod.o            a6_read_mod.o           acetone_mod.o		     \
aerosol_mod.o            aeic_mod.o              airmas.o                    \
anthroems.o              apm_coag_mod.o          apm_grow_mod.o		     \
apm_nucl_mod.o           apm_phys_mod.o          apm_driv_mod.o  	     \
apm_init_mod.o           apm_wets_mod.o		                             \
arctas_ship_emiss_mod.o  arsl1k.o                backsub.o                   \
benchmark_mod.o          biofit.o                biofuel_mod.o               \
biomass_mod.o            boxvl.o                 bravo_mod.o                 \
c2h6_mod.o               cac_anthro_mod.o        calcrate.o                  \
carbon_mod.o             ch3i_mod.o              chemdr.o	             \
chemistry_mod.o          cleanup.o               co2_mod.o                   \
comode_mod.o             convection_mod.o        dao_mod.o                   \
decomp.o                 depo_mercury_mod.o      diag03_mod.o                \
diag04_mod.o             diag1.o                 diag3.o                     \
diag41_mod.o             diag42_mod.o            diag48_mod.o                \
diag49_mod.o             diag50_mod.o            diag51_mod.o                \
diag51b_mod.o            diag56_mod.o            diag_2pm.o		     \
diag_mod.o               diag_oh_mod.o           diag_pl_mod.o		     \
diagoh.o                 drydep_mod.o            dust_dead_mod.o	     \
dust_mod.o               edgar_mod.o             emep_mod.o		     \
emf_scale.o              emfossil.o              emissdr.o		     \
emissions_mod.o          epa_nei_mod.o		 mapping_mod.o               \
fast_j.o                 fcro2ho2.o                                          \
fjfunc.o                 fjx_acet_mod.o          future_emissions_mod.o	     \
fvdas_convect_mod.o      fyhoro.o                fyrno3.o		     \
gamap_mod.o              gasconc.o               gc_biomass_mod.o	     \
gc_type_mod.o            gcap_convect_mod.o      gcap_read_mod.o	     \
geia_mod.o               get_global_ch4.o                		     \
gfed2_biomass_mod.o      gfed3_biomass_mod.o     geosfp_read_mod.o           \
global_br_mod.o          global_ch4_mod.o	 olson_landmap_mod.o         \
global_hno3_mod.o        global_no3_mod.o        global_nox_mod.o	     \
global_o1d_mod.o         global_o3_mod.o         global_oh_mod.o	     \
h2_hd_mod.o              hcn_ch3cn_mod.o	                             \
i6_read_mod.o            icoads_ship_mod.o       initialize.o		     \
inphot.o                 input_mod.o             isoropiaII_mod.o	     \
jsparse.o                jv_index.o              ksparse.o		     \
land_mercury_mod.o       lightning_nox_mod.o	                             \
linoz_mod.o              logical_mod.o           lump.o			     \
main.o                   megan_mod.o             meganut_mod.o		     \
mercury_mod.o            merra_a1_mod.o          merra_a3_mod.o		     \
merra_cn_mod.o           merra_i6_mod.o          mmran_16.o		     \
ndxx_setup.o             nei2005_anthro_mod.o    ocean_mercury_mod.o	     \
ohsave.o                 optdepth_mod.o          partition.o		     \
pbl_mix_mod.o            pderiv.o                photoj.o		     \
physproc.o               pjc_pfix_mod.o          pjc_pfix_geos5_window_mod.o \
planeflight_mod.o        precipfrac.o                   		     \
rd_js.o                  readchem.o	                                     \
reader.o                 restart_mod.o		                             \
rpmares_mod.o            ruralbox.o              scale_anthro_mod.o	     \
seasalt_mod.o            set_aer.o	                                     \
set_prof.o               setemdep.o		 modis_lai_mod.o             \
setemis.o                setmodel.o              sfcwindsqr.o		     \
smvgear.o                                                		     \
streets_anthro_mod.o     subfun.o	                         	     \
sulfate_mod.o            sunparam.o              tagged_co_mod.o	     \
tagged_ox_mod.o          toms_mod.o		                             \
tpcore_bc_mod.o          tpcore_fvdas_mod.o      tpcore_geos5_window_mod.o   \
tracer_mod.o                                                                 \
tracerid_mod.o           transport_mod.o	 tropopause_mod.o            \
update.o                 uvalbedo_mod.o                                      \
vdiff_pre_mod.o          vdiff_mod.o             vistas_anthro_mod.o         \
wetscav_mod.o            retro_mod.o             diag63_mod.o                \
paranox_mod.o            read_jv_atms_dat.o      bromocarb_mod.o             \
cldice_HBrHOBr_rxn.o     ssa_bromine_mod.o       strat_chem_mod.o            \
gigc_environment_mod.o                                                       \
canopy_nox_mod.o         get_ndep_mod.o          nc_soilnox_read.o           \
readsoil.o               readfert.o              readclim.o                  \
soil_nox_mod.o           soilnox_restart_mod.o                               \
diag53_mod.o             global_oc_mod.o         global_bc_mod.o             \
pops_mod.o               get_popsinfo_mod.o      rcp_mod.o                   \
tpcore_geosfp_window_mod.o                                                   \
pjc_pfix_geosfp_window_mod.o

#=============================================================================
# Makefile targets: type "make help" for a complete listing!
#=============================================================================

.PHONY: clean realclean doc docclean help

all:
	@$(MAKE) lib
	@$(MAKE) exe

lib:                                               # Build all G-C libraries
	@$(MAKE) libnc                             
	@$(MAKE) libheaders         		     
	@$(MAKE) libkpp             		     
	@$(MAKE) libutil 			     
	@$(MAKE) libiso 			     
	@$(MAKE) libapm

libapm: $(OBJ)                                     # Build code in GeosApm/

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
	rm -f *.o *.mod geos geosapm

realclean:                                         # Remove files everywhere
	@$(MAKE) clean
	@$(MAKE) -C $(ISO) clean
	@$(MAKE) -C $(KPP) realclean
	@$(MAKE) -C $(NCDF) clean
	@$(MAKE) -C $(UTIL) clean
	@$(MAKE) -C $(CORE) clean
	@$(MAKE) -C $(HDR)  clean
	@$(MAKE) -C $(GTMM) clean
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
# Dependencies listing (grep "USE " to get the list of module references!)
#
# From this list of dependencies, the "make" utility will figure out the
# correct order of compilation (so we don't have to do that ourselves!)
#=============================================================================

a3_read_mod.o               : a3_read_mod.F                                  \
                              dao_mod.o               diag_mod.o

a6_read_mod.o               : a6_read_mod.F                                  \
                              dao_mod.o               diag_mod.o              

acetone_mod.o               : acetone_mod.F                                  \
                              diag_mod.o              megan_mod.o 

aeic_mod.o                  : aeic_mod.F                                     \
                              diag_mod.o              tracerid_mod.o         \
                              tropopause_mod.o

aerosol_mod.o               : aerosol_mod.F                                  \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              tracerid_mod.o         \
                              tropopause_mod.o

airmas.o                    : airmas.F

anthroems.o                 : anthroems.F                                    \
                              edgar_mod.o             future_emissions_mod.o \
                              geia_mod.o              scale_anthro_mod.o     \
                              tracerid_mod.o

apm_coag_mod.o              : apm_coag_mod.F                                 \
                              apm_init_mod.o

apm_driv_mod.o              : apm_driv_mod.F                                 \
                              apm_coag_mod.o          apm_init_mod.o         \
                              apm_nucl_mod.o          apm_phys_mod.o         \
                              comode_mod.o            dao_mod.o              \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o

apm_grow_mod.o              : apm_grow_mod.F                                 \
                              apm_init_mod.o

apm_init_mod.o              : apm_init_mod.F                                 \
                              apm_nucl_mod.o

apm_nucl_mod.o              : apm_nucl_mod.F

apm_phys_mod.o              : apm_phys_mod.F                                 \
                              apm_coag_mod.o          apm_grow_mod.o         \
                              apm_init_mod.o          apm_nucl_mod.o

apm_wets_mod.o              : apm_wets_mod.F                                 \
                              apm_driv_mod.o          apm_init_mod.F         \
                              dao_mod.o               diag_mod.o             \
                              wetscav_mod.o

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

BLKSLV.o                    : BLKSLV.F

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
                              tracerid_mod.o          cldice_HBrHOBr_rxn.o

canopy_nox_mod.o            : canopy_nox_mod.F                               \
                              drydep_mod.o            modis_lai_mod.o

carbon_mod.o                : carbon_mod.F                                   \
                              biomass_mod.o           comode_mod.o           \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            future_emissions_mod.o \
                              gfed2_biomass_mod.o     gfed3_biomass_mod.o    \
                              global_no3_mod.o        global_o3_mod.o        \
                              global_oh_mod.o         megan_mod.o            \
                              meganut_mod.o           pbl_mix_mod.o          \
                              rcp_mod.o               tomas_mod.o            \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o        vdiff_pre_mod.o        \
                              apm_init_mod.F          sulfate_mod.o

ch3i_mod.o                  : ch3i_mod.F                                     \
                              biofuel_mod.o           biomass_mod.o          \
                              diag_mod.o              diag_pl_mod.o          \
                              tracerid_mod.o

chemdr.o                    : chemdr.F                                       \
                              aerosol_mod.o           comode_mod.o           \
                              dao_mod.o               diag_oh_mod.o          \
                              diag_pl_mod.o           dust_mod.o             \
                              future_emissions_mod.o  planeflight_mod.o      \
                              restart_mod.o           tracerid_mod.o         \
                              tropopause_mod.o

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
                              tomas_mod.o             tracerid_mod.o         \
                              apm_init_mod.o          apm_driv_mod.o

cleanup.o                   : cleanup.F                                      \
                              acetone_mod.o           aerosol_mod.o          \
                              arctas_ship_emiss_mod.o aeic_mod.o             \
                              biomass_mod.o           biofuel_mod.o          \
                              bravo_mod.o             bromocarb_mod.o        \
                              c2h6_mod.o              cac_anthro_mod.o       \
                              carbon_mod.o            co2_mod.o              \
                              comode_mod.o            dao_mod.o              \
                              depo_mercury_mod.o      diag_mod.o             \
                              diag03_mod.o            diag04_mod.o           \
                              diag41_mod.o            diag50_mod.o           \
                              diag51_mod.o            diag51b_mod.o          \
                              diag53_mod.o            diag_oh_mod.o          \
                              diag_pl_mod.o           drydep_mod.o           \
                              dust_mod.o              dust_dead_mod.o        \
                              edgar_mod.o             emep_mod.o             \
                              epa_nei_mod.o           gc_biomass_mod.o       \
                              gfed2_biomass_mod.o     gfed3_biomass_mod.o    \
                              global_ch4_mod.o        global_hno3_mod.o      \
                              global_no3_mod.o        global_nox_mod.o       \
                              global_o1d_mod.o        global_oh_mod.o        \
                              global_bc_mod.o         global_oc_mod.o        \
                              h2_hd_mod.o             hcn_ch3cn_mod.o        \
                              icoads_ship_mod.o       isoropiaII_mod.o       \
                              land_mercury_mod.o      lightning_nox_mod.o    \
                              linoz_mod.o             megan_mod.o            \
                              mercury_mod.o           modis_lai_mod.o        \
                              ocean_mercury_mod.o     pbl_mix_mod.o          \
                              pjc_pfix_mod.o          planeflight_mod.o      \
                              rcp_mod.o               retro_mod.o            \
                              seasalt_mod.o           sulfate_mod.o          \
                              strat_chem_mod.o        tagged_co_mod.o        \
                              tomas_mod.o             toms_mod.o             \
                              tracer_mod.o            transport_mod.o        \
                              tropopause_mod.o        uvalbedo_mod.o         \
                              vistas_anthro_mod.o     wetscav_mod.o          \
                              tpcore_fvdas_mod.o      apm_wets_mod.o         \
                              tpcore_geos5_window_mod.o                      \
                              tpcore_geosfp_window_mod.o

cldice_HBrHOBr_rxn.o        : cldice_HBrHOBr_rxn.F

CLDSRF.o                    : CLDSRF.F

co2_mod.o                   : co2_mod.F                                      \
                              biomass_mod.o           diag04_mod.o           \
	                      tracerid_mod.o          

comode_mod.o                : comode_mod.F

convection_mod.o            : convection_mod.F                               \
                              diag_mod.o              depo_mercury_mod.o     \
                              fvdas_convect_mod.o     gc_type_mod.o          \
                              gcap_convect_mod.o      mercury_mod.o          \
                              tracer_mod.o            tracerid_mod.o         \
                              wetscav_mod.o          

dao_mod.o                   : dao_mod.F                                      \
                              logical_mod.o           tracer_mod.o

decomp.o                    : decomp.F

depo_mercury_mod.o          : depo_mercury_mod.F                             \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracerid_mod.o         

diag03_mod.o                : diag03_mod.F                                   \
                              logical_mod.o           tracerid_mod.o

diag04_mod.o                : diag04_mod.F

diag1.o                     : diag1.F                                        \
                              dao_mod.o               diag03_mod.o           \
                              diag_mod.o              tracer_mod.o           \
                              tracerid_mod.o          tropopause_mod.o

diag3.o                     : diag3.F                                        \
                              biofuel_mod.o           biomass_mod.o          \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o            diag53_mod.o           \
                              diag56_mod.o            diag_pl_mod.o          \
                              depo_mercury_mod.o      drydep_mod.o           \
                              logical_mod.o           tomas_mod.o            \
                              tracerid_mod.o          wetscav_mod.o          \
                              apm_driv_mod.o

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
                              tracerid_mod.o          tracer_mod.o

diag50_mod.o                : diag50_mod.F                                   \
                              comode_mod.o            logical_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o

diag51_mod.o                : diag51_mod.F                                   \
                              modis_lai_mod.o         pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o  

diag51b_mod.o               : diag51b_mod.F                                  \
                              modis_lai_mod.o         pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o  

diag53_mod.o                : diag53_mod.F

diag56_mod.o                : diag56_mod.F

diag63_mod.o                : diag63_mod.F                                   \
                              diag_mod.o              pbl_mix_mod.o

diag_2pm.o                  : diag_2pm.F                                     \
                              diag_mod.o              tropopause_mod.o


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
                              tracerid_mod.o          apm_init_mod.o         \
                              apm_driv_mod.o

diagoh.o                    : diagoh.F                                       \
                              diag_mod.o

dust_dead_mod.o             : dust_dead_mod.F                                \
                              dao_mod.o               apm_init_mod.o

dust_mod.o                  : dust_mod.F                                     \
                              comode_mod.o            diag_mod.o             \
                              drydep_mod.o            dust_dead_mod.o        \
                              tomas_mod.o             tracerid_mod.o         \
                              apm_init_mod.o

edgar_mod.o                 : edgar_mod.F                                    \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o

EFOLD.o                     : EFOLD.F

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
                              vistas_anthro_mod.o

epa_nei_mod.o               : epa_nei_mod.F                                  \
                              future_emissions_mod.o  scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o

fast_j.o                    : fast_j.F                                       \
                              toms_mod.o

fjfunc.o                    : fjfunc.F

fjx_acet_mod.o              : fjx_acet_mod.F

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
                              comode_mod.o            drydep_mod.o           \
                              tropopause_mod.o

GAUSSP.o                    : GAUSSP.F

gc_biomass_mod.o            : gc_biomass_mod.F                               \
                              future_emissions_mod.o  logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o    

gc_type_mod.o               : gc_type_mod.F

gcap_convect_mod.o          : gcap_convect_mod.F                             \
                              diag_mod.o

gcap_read_mod.o             : gcap_read_mod.F                                \
                              diag_mod.o

geia_mod.o                  : geia_mod.F

GEN.o                       : GEN.F

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
                              ocean_mercury_mod.o     tropopause_mod.o 

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
                              diag_mod.o              drydep_mod.o           \
                              geia_mod.o              global_nox_mod.o       \
                              global_o1d_mod.o        global_oh_mod.o        \
                              megan_mod.o             meganut_mod.o          \
                              scale_anthro_mod.o      tagged_co_mod.o        \
                              tracerid_mod.o          tropopause_mod.o         

hcn_ch3cn_mod.o             : hcn_ch3cn_mod.F                                \
                              biomass_mod.o           diag_mod.o             \
                              geia_mod.o              global_oh_mod.o        \
                              pbl_mix_mod.o           tracerid_mod.o

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

inphot.o                    : inphot.F

input_mod.o                 : input_mod.F                                    \
                              aerosol_mod.o           benchmark_mod.o        \
                              biofuel_mod.o           biomass_mod.o          \
                              bromocarb_mod.o         carbon_mod.o           \
                              depo_mercury_mod.o      diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o  	      diag48_mod.o           \
                              diag49_mod.o  	      diag50_mod.o           \
                              diag51_mod.o  	      diag51b_mod.o          \
                              diag56_mod.o  	      diag63_mod.o           \
                              diag_oh_mod.o           diag_pl_mod.o          \
                              drydep_mod.o            dust_mod.o             \
                              emissions_mod.o         future_emissions_mod.o \
                              gamap_mod.o             get_popsinfo_mod.o     \
                              land_mercury_mod.o      logical_mod.o          \
                              mercury_mod.o           modis_lai_mod.o        \
                              ocean_mercury_mod.o     planeflight_mod.o      \
                              pops_mod.o              rcp_mod.o              \
                              restart_mod.o           seasalt_mod.o          \
                              sulfate_mod.o           tomas_mod.o            \
                              tpcore_bc_mod.o         tracerid_mod.o         \
                              tracer_mod.o            transport_mod.o        \
                              tropopause_mod.o        wetscav_mod.o          \
                              linoz_mod.o             apm_init_mod.o

isoropiaII_mod.o            : isoropiaII_mod.F                               \
                              global_hno3_mod.o       tracerid_mod.o         \
                              tropopause_mod.o        apm_init_mod.o

jsparse.o                   : jsparse.F               comode_mod.o

JRATET.o                    : JRATET.F                fjx_acet_mod.o

JVALUE.o                    : JVALUE.F

jv_index.o                  : jv_index.F

ksparse.o                   : ksparse.F

land_mercury_mod.o          : land_mercury_mod.F                             \
                              biomass_mod.o           dao_mod.o              \
                              depo_mercury_mod.o      modis_lai_mod.o        \
                              logical_mod.o           tracerid_mod.o         

LEGND0.o                    : LEGND0.F

lightning_nox_mod.o         : lightning_nox_mod.F                            \
                              dao_mod.o               diag56_mod.o           \
                              diag_mod.o              logical_mod.o          \
                              tracerid_mod.o          tropopause_mod.o

linoz_mod.o                 : linoz_mod.F                                    \
                              tracerid_mod.o          tropopause_mod.o  

logical_mod.o               : logical_mod.F

lump.o                      : lump.F                                         \
                              comode_mod.o            tracerid_mod.o

main.o                      : main.F                                         \
                              a3_read_mod.o           a6_read_mod.o          \
                              benchmark_mod.o         bromocarb_mod.o        \
                              chemistry_mod.o         carbon_mod.o           \
                              convection_mod.o        comode_mod.o           \
                              diag_mod.o              diag41_mod.o           \
                              diag42_mod.o            diag48_mod.o           \
                              diag49_mod.o            diag50_mod.o           \
                              diag51_mod.o            diag51b_mod.o          \
                              diag63_mod.o            diag_oh_mod.o          \
                              dao_mod.o               depo_mercury_mod.o     \
                              drydep_mod.o            emissions_mod.o        \
                              gcap_read_mod.o         geosfp_read_mod.o      \
                              gigc_environment_mod.o  global_ch4_mod.o       \
                              i6_read_mod.o           input_mod.o            \
                              lightning_nox_mod.o     linoz_mod.o            \
                              logical_mod.o           mapping_mod.o          \
                              megan_mod.o             mercury_mod.o          \
                              merra_cn_mod.o          merra_a1_mod.o         \
                              merra_a3_mod.o          merra_i6_mod.o         \
                              modis_lai_mod.o         ocean_mercury_mod.o    \
                              olson_landmap_mod.o     pbl_mix_mod.o          \
                              planeflight_mod.o       restart_mod.o          \
                              soilnox_restart_mod.o   strat_chem_mod.o       \
                              toms_mod.o              tpcore_bc_mod.o        \
                              tracer_mod.o            tracerid_mod.o         \
                              transport_mod.o         tropopause_mod.o       \
                              uvalbedo_mod.o          vdiff_mod.o            \
                              wetscav_mod.o           apm_init_mod.o         \
                              apm_driv_mod.o          apm_wets_mod.o

mapping_mod.o               : mapping_mod.F90                                \
                              logical_mod.o

MATIN4.o                    : MATIN4.F

megan_mod.o                 : megan_mod.F             a3_read_mod.o          \
                              modis_lai_mod.o         logical_mod.o          \
                              meganut_mod.o           merra_a1_mod.o         \
                              geosfp_read_mod.o

meganut_mod.o               : meganut_mod.F

mercury_mod.o               : mercury_mod.F                                  \
                              dao_mod.o               depo_mercury_mod.o     \
                              diag03_mod.o            diag_mod.o             \
                              drydep_mod.o            global_br_mod.o        \
                              global_o3_mod.o         global_oh_mod.o        \
                              land_mercury_mod.o      logical_mod.o          \
                              ocean_mercury_mod.o     pbl_mix_mod.o          \
                              RnPbBe_mod.o            tracerid_mod.o         \
                              tropopause_mod.o        vdiff_pre_mod.o

merra_a1_mod.o              : merra_a1_mod.F                                 \
                              diag_mod.o

merra_a3_mod.o              : merra_a3_mod.F                                 \
                              dao_mod.o               diag_mod.o

merra_cn_mod.o              : merra_cn_mod.F                                 \
                              diag_mod.o

merra_i6_mod.o              : merra_i6_mod.F

MIESCT.o                    : MIESCT.F

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

NOABS.o                     : NOABS.F

ocean_mercury_mod.o         : ocean_mercury_mod.F                            \
                              dao_mod.o               depo_mercury_mod.o     \
                              diag03_mod.o            logical_mod.o          \
                              tracerid_mod.o

ohsave.o                    : ohsave.F                                       \
                              comode_mod.o            tracerid_mod.o

olson_landmap_mod.o         : olson_landmap_mod.F90                          \
                              logical_mod.o           mapping_mod.o

optdepth_mod.o              : optdepth_mod.F          diag_mod.o  

OPMIE.o                     : OPMIE.F

paranox_mod.o               : paranox_mod.F

partition.o                 : partition.F                                    \
                              comode_mod.o            tracerid_mod.o

pbl_mix_mod.o               : pbl_mix_mod.F                                  \
                              diag_mod.o              tracer_mod.o

pderiv.o                    : pderiv.F

photoj.o                    : photoj.F

physproc.o                  : physproc.F                                     \
                              comode_mod.o            chemistry_mod.o

pjc_pfix_geos5_window_mod.o : pjc_pfix_geos5_window_mod.F

pjc_pfix_geosfp_window_mod.o: pjc_pfix_geosfp_window_mod.F

pjc_pfix_mod.o              : pjc_pfix_mod.F

planeflight_mod.o           : planeflight_mod.F                              \
                              comode_mod.o            diag_mod.o             \
                              ocean_mercury_mod.o     tracer_mod.o

pops_mod.o           	    : pops_mod.F    	                             \
			      diag_mod.o	      diag53_mod.o           \
                              drydep_mod.o            get_popsinfo_mod.o     \
                              global_bc_mod.o	      global_oc_mod.o        \
	        	      global_o3_mod.o         global_oh_mod.o        \
                              pbl_mix_mod.o	      tracerid_mod.o         \
                              vdiff_pre_mod.o

precipfrac.o                : precipfrac.F

RD_AOD.o                    : RD_AOD.F

rd_js.o                     : rd_js.F

rcp_mod.o                   : rcp_mod.F                                      \
                              diag_mod.o              tracerid_mod.o         \
                              tracer_mod.o

readsoil.o		    : readsoil.F

readfert.o		    : readfert.F

readclim.o		    : readclim.F

RD_TJPL.o                   : RD_TJPL.F

readchem.o                  : readchem.F                                     \
                              diag_pl_mod.o           drydep_mod.o

reader.o                    : reader.F

read_jv_atms_dat.o          : read_jv_atms_dat.F90

restart_mod.o               : restart_mod.F                                  \
                              comode_mod.o            logical_mod.o          \
                              tracer_mod.o            tracerid_mod.o         \
                              apm_init_mod.o

retro_mod.o                 : retro_mod.F                                    \
                              future_emissions_mod.o  scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o

RnPbBe_mod.o                : RnPbBe_mod.F                                   \
                              diag_mod.o              tropopause_mod.o       \
                              vdiff_pre_mod.o

rpmares_mod.o               : rpmares_mod.F                                  \
                              global_hno3_mod.o       tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o       \
                              apm_init_mod.o

ruralbox.o                  : ruralbox.F                                     \
                              comode_mod.o            tropopause_mod.o

scale_anthro_mod.o          : scale_anthro_mod.F      

seasalt_mod.o               : seasalt_mod.F                                  \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            logical_mod.o          \
                              pbl_mix_mod.o           ssa_bromine_mod.o      \
                              tomas_mod.o             tracerid_mod.o         \
                              tracer_mod.o            vdiff_pre_mod.o        \
                              apm_init_mod.o          apm_driv_mod.o

set_aer.o                   : set_aer.F

set_prof.o                  : set_prof.F

setemdep.o                  : setemdep.F                                     \
                              tracerid_mod.o

setemis.o                   : setemis.F                                      \
                              biofuel_mod.o           biomass_mod.o          \
                              comode_mod.o            diag_mod.o             \
                              diag63_mod.o            drydep_mod.o           \
                              emissions_mod.o         edgar_mod.o            \
                              emep_mod.o              lightning_nox_mod.o    \
                              paranox_mod.o           pbl_mix_mod.o          \
                              rcp_mod.o               tracerid_mod.o         \
                              tropopause_mod.o

setmodel.o                  : setmodel.F

sfcwindsqr.o                : sfcwindsqr.F

smvgear.o                   : smvgear.F                                      \
                              comode_mod.o

soil_nox_mod.o              : soil_nox_mod.F

soilnox_restart_mod.o       : soilnox_restart_mod.F

SPHERE.o                    : SPHERE.F

ssa_bromine_mod.o           : ssa_bromine_mod.F                              \
                              bromocarb_mod.o         comode_mod.o           \
                              diag_mod.o              logical_mod.o          \
                              tracerid_mod.o

strat_chem_mod.o            : strat_chem_mod.F90                             \
                              dao_mod.o               logical_mod.o          \
                              linoz_mod.o             tracer_mod.o           \
                              tracerid_mod.o          tropopause_mod.o       \
                              tagged_ox_mod.o

streets_anthro_mod.o        : streets_anthro_mod.F                           \
                              dao_mod.o               future_emissions_mod.o \
                              logical_mod.o           scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o            

subfun.o                    : subfun.F

sulfate_mod.o               : sulfate_mod.F                                  \
                              arctas_ship_emiss_mod.o biomass_mod.o          \
                              bravo_mod.o             cac_anthro_mod.o       \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              drydep_mod.o           \
                              edgar_mod.o             emep_mod.o             \
                              epa_nei_mod.o           future_emissions_mod.o \
                              get_ndep_mod.o          gfed2_biomass_mod.o    \
                              gfed3_biomass_mod.o     global_hno3_mod.o      \
                              global_no3_mod.o        global_oh_mod.o        \
                              icoads_ship_mod.o       logical_mod.o          \
                              nei2005_anthro_mod.o    pbl_mix_mod.o          \
                              rcp_mod.o               scale_anthro_mod.o     \
                              seasalt_mod.o           streets_anthro_mod.o   \
                              tomas_mod.o             tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o       \
                              uvalbedo_mod.o          vdiff_pre_mod.o        \
                              wetscav_mod.o           apm_init_mod.o         \
                              apm_driv_mod.o

sunparam.o                  : sunparam.F

tagged_co_mod.o             : tagged_co_mod.F                                \
                              biofuel_mod.o           biomass_mod.o          \
                              diag_mod.o              diag_pl_mod.o          \
                              global_nox_mod.o        global_oh_mod.o        \
                              megan_mod.o             meganut_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o 

tagged_ox_mod.o             : tagged_ox_mod.F                                \
                              diag_mod.o              diag_pl_mod.o          \
                              drydep_mod.o            pbl_mix_mod.o          \
                              tracerid_mod.o          tropopause_mod.o

tcorr.o                     : tcorr.F

toms_mod.o                  : toms_mod.F

tpcore_bc_mod.o             : tpcore_bc_mod.F                                \
                              tracer_mod.o

tpcore_geos5_window_mod.o   : tpcore_geos5_window_mod.F90
	$(F90) -c $(FREEFORM) $(R8) $<

tpcore_geosfp_window_mod.o  : tpcore_geosfp_window_mod.F90
	$(F90) -c $(FREEFORM) $(R8) $<

tracer_mod.o                : tracer_mod.F                                   \
                              apm_init_mod.o

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

tropopause_mod.o            : tropopause_mod.F                               \
                              comode_mod.o            diag_mod.o             \
                              logical_mod.o

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
                              tracer_mod.o            tomas_mod.o            \
                              apm_init_mod.o

XSEC1D.o                    : XSEC1D.F

XSECO2.o                    : XSECO2.F

XSECO3.o                    : XSECO3.F

#------------------------------------------------------------------------------
# Dependencies of files specific to TOMAS microphysics (bmy, 1/25/10)
#------------------------------------------------------------------------------

tomas_mod.o                 : tomas_mod.F                                    \
                              dao_mod.o               diag_mod.o             \
                              diag_pl_mod.o           tropopause_mod.o       \
                              tracerid_mod.o          ion_nucl.o             \
                              YuIMN_Code.o            cf_nucl.o

tomas_tpcore_mod.o          : tomas_tpcore_mod.F90                           \
                              tomas_mod.o             tracerid_mod.o


aero_drydep.o               : aero_drydep.F                                  \
                              diag_mod.o              drydep_mod.o           \
                              dust_mod.o              pbl_mix_mod.o          \
                              tomas_mod.o             tracerid_mod.o         \
                              tropopause_mod.o

ion_nucl.o                  : ion_nucl.F

YuIMN_Code.o                : YuIMN_Code.F

cf_nucl.o                   : cf_nucl.F

#EOC
