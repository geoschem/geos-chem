# $Id: Makefile,v 1.23 2010/03/15 19:33:25 ccarouge Exp $
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
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
ROOTDIR = ..
BIN     = $(ROOTDIR)/bin
BPCH    = $(ROOTDIR)/GeosBpch
DOC     = $(ROOTDIR)/doc
EXE     = geos
HDR     = $(ROOTDIR)/Headers
HELP    = $(ROOTDIR)/help
ISO     = $(ROOTDIR)/ISOROPIA
LIB     = $(ROOTDIR)/lib
KPP     = $(ROOTDIR)/KPP
MOD     = $(ROOTDIR)/mod
TOM     = $(ROOTDIR)/GeosTomas
GTMM    = $(ROOTDIR)/GTMM
UTIL    = $(ROOTDIR)/GeosUtil

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

#=============================================================================
# List of files to compile (the order is important!).  We specify these as
# a list of object files (*.o).  For each object file, the "make" utility
# will find the corresponding source code file (*.f, *.f90) and compile it. 
#=============================================================================

# List of source files
SRC = $(wildcard *.f) $(wildcard *.f90)

# Replace .f and .f90 extensions with *.o
TMP = $(SRC:.f=.o)
OBJ = $(TMP:.f90=.o)

#=============================================================================
# Makefile targets: type "make help" for a complete listing!
#=============================================================================

.PHONY: clean realclean doc docclean help

all:
	@$(MAKE) lib
	@$(MAKE) exe

lib: libkpp libutil libiso libcore

libcore: $(OBJ)

libiso:
	@$(MAKE) -C $(ISO)

libkpp: 
	@$(MAKE) -C $(KPP)

libutil:
	@$(MAKE) -C $(UTIL)

exe:
	$(LD) $(OBJ) $(LINK) -o $(EXE)
	cp -f $(EXE) $(BIN)

clean:
	rm -f *.o *.mod geos

realclean:
	@$(MAKE) clean
	@$(MAKE) -C $(ISO) clean
	@$(MAKE) -C $(KPP) realclean
	@$(MAKE) -C $(TOM) clean
	@$(MAKE) -C $(GTMM) clean
	@$(MAKE) -C $(UTIL) clean
	@$(MAKE) docclean
	rm -f $(LIB)/*.a
	rm -f $(MOD)/*.mod
	rm -f $(BIN)/geos*

doc:
	@$(MAKE) -C $(DOC) all

docclean: 
	@$(MAKE) -C $(DOC) clean

help:
	@$(MAKE) -C $(HELP)

#=============================================================================
# Targets for mercury simulation (ccc, 6/7/10)
#=============================================================================

allhg:
	@$(MAKE) libhg
	@$(MAKE) exehg

libhg: libutil libkpp libgtmm libiso libcore

libgtmm:
	@$(MAKE) -C $(GTMM)

exehg:
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

a3_read_mod.o               : a3_read_mod.f           dao_mod.o              \
                              diag_mod.o              logical_mod.o

a6_read_mod.o               : a6_read_mod.f           dao_mod.o              \
                              diag_mod.o              logical_mod.o

acetone_mod.o               : acetone_mod.f           dao_mod.o              \
                              diag_mod.o               

aerosol_mod.o               : aerosol_mod.f           comode_mod.o           \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o  

aircraft_nox_mod.o          : aircraft_nox_mod.f                             \
                              dao_mod.o               diag_mod.o              

airmas.o                    : airmas.f

anthroems.o                 : anthroems.f                                    \
                              future_emissions_mod.o  edgar_mod.o            \
                              geia_mod.o              logical_mod.o          \
                              scale_anthro_mod.o      tracer_mod.o           \
                              tracerid_mod.o

arctas_ship_emiss_mod.o     : arctas_ship_emiss_mod.f                        \
                              logical_mod.o           scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o  

arsl1k.o                    : arsl1k.f

backsub.o                   : backsub.f

benchmark_mod.o             : benchmark_mod.f                                \
                              tracerid_mod.o          tracer_mod.o  

biofuel_mod.o               : biofuel_mod.f                                  \
                              dao_mod.o               diag_mod.o             \
                              epa_nei_mod.o           future_emissions_mod.o \
                              logical_mod.o           streets_anthro_mod.o   \
                              tracerid_mod.o          tracer_mod.o    

biofit.o                    : biofit.f

biomass_mod.o               : biomass_mod.f                                  \
                              diag_mod.o              gc_biomass_mod.o       \
                              gfed2_biomass_mod.o     logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o

BLKSLV.o                    : BLKSLV.f

boxvl.o                     : boxvl.f                 dao_mod.o

bravo_mod.o                 : bravo_mod.f                                    \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o    

c2h6_mod.o                  : c2h6_mod.f                                     \
                              biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o             \
                              geia_mod.o              global_oh_mod.o        \
                              logical_mod.o           tracerid_mod.o         \
                              tracer_mod.o    

cac_anthro_mod.o            : cac_anthro_mod.f                               \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o  

calcrate.o                  : calcrate.f                                     \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              drydep_mod.o           \
                              logical_mod.o           pbl_mix_mod.o          \
                              planeflight_mod.o       tracerid_mod.o

carbon_mod.o                : carbon_mod.f                                   \
                              biomass_mod.o           comode_mod.o           \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            future_emissions_mod.o \
                              gfed2_biomass_mod.o     global_no3_mod.o       \
                              global_o3_mod.o         global_oh_mod.o        \
                              logical_mod.o           megan_mod.o            \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o       \
                              vdiff_pre_mod.o         meganut_mod.o

ch3i_mod.o                  : ch3i_mod.f                                     \
                              biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o             \
                              diag_pl_mod.o           logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              uvalbedo_mod.o

chemdr.o                    : chemdr.f                                       \
                              aerosol_mod.o           comode_mod.o           \
                              dao_mod.o               diag_oh_mod.o          \
                              diag_pl_mod.o           future_emissions_mod.o \
                              dust_mod.o              logical_mod.o          \
                              planeflight_mod.o       restart_mod.o          \
                              tracer_mod.o            tracerid_mod.o         \
                              tropopause_mod.o        uvalbedo_mod.o 

chemistry_mod.o             : chemistry_mod.f                                \
                              acetone_mod.o           aerosol_mod.o          \
                              c2h6_mod.o              carbon_mod.o           \
                              ch3i_mod.o              comode_mod.o           \
                              dao_mod.o               drydep_mod.o           \
                              dust_mod.o              global_ch4_mod.o       \
                              h2_hd_mod.o             hcn_ch3cn_mod.o        \
                              isoropiaII_mod.o        logical_mod.o          \
                              mercury_mod.o           optdepth_mod.o         \
                              rpmares_mod.o           RnPbBe_mod.o           \
                              seasalt_mod.o           sulfate_mod.o          \
                              tagged_co_mod.o         tagged_ox_mod.o        \
                              tracerid_mod.o          tracer_mod.o

cleanup.o                   : cleanup.f                                      \
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
                              global_ch4_mod.o        global_hno3_mod.o      \
                              global_no3_mod.o        global_nox_mod.o       \
                              global_o1d_mod.o        global_oh_mod.o        \
                              h2_hd_mod.o             hcn_ch3cn_mod.o        \
                              lai_mod.o               land_mercury_mod.o     \
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
                              tpcore_geos5_window_mod.o

CLDSRF.o                    : CLDSRF.f

co2_mod.o                   : co2_mod.f               biomass_mod.o          \
                              diag04_mod.o            tracer_mod.o           \
	                      tracerid_mod.o 

CO_strat_pl.o               : CO_strat_pl.f                                  \
                              dao_mod.o               tracer_mod.o           \
                              tracerid_mod.o          tropopause_mod.o 

comode_mod.o                : comode_mod.f            tracer_mod.o

convection_mod.o            : convection_mod.f                               \
                              dao_mod.o               diag_mod.o             \
                              depo_mercury_mod.o      fvdas_convect_mod.o    \
                              gcap_convect_mod.o      logical_mod.o          \
                              mercury_mod.o           tracer_mod.o           \
                              tracerid_mod.o          wetscav_mod.o

dao_mod.o                   : dao_mod.f                                      \
                              logical_mod.o           tracer_mod.o           

decomp.o                    : decomp.f

depo_mercury_mod.o          : depo_mercury_mod.f                             \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracerid_mod.o         

diag03_mod.o                : diag03_mod.f            tracerid_mod.o

diag04_mod.o                : diag04_mod.f

diag1.o                     : diag1.f                                        \
                              dao_mod.o               diag_mod.o             \
                              tracer_mod.o            tracerid_mod.o         \
                              tropopause_mod.o

diag3.o                     : diag3.f                                        \
                              biofuel_mod.o           biomass_mod.o          \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o            diag56_mod.o           \
                              diag_pl_mod.o           depo_mercury_mod.o     \
                              drydep_mod.o            logical_mod.o          \
                              tracer_mod.o            tracerid_mod.o         \
                              wetscav_mod.o

diag41_mod.o                : diag41_mod.f            pbl_mix_mod.o    

diag42_mod.o                : diag42_mod.f                                   \
                              dao_mod.o               logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

diag48_mod.o                : diag48_mod.f                                   \
                              dao_mod.o               pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

diag49_mod.o                : diag49_mod.f                                   \
                              dao_mod.o               lai_mod.o              \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o  

diag50_mod.o                : diag50_mod.f                                   \
                              comode_mod.o            dao_mod.o              \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o  

diag51_mod.o                : diag51_mod.f                                   \
                              dao_mod.o               lai_mod.o              \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o  

diag51b_mod.o               : diag51b_mod.f                                  \
                              dao_mod.o               lai_mod.o              \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o  

diag56_mod.o                : diag56_mod.f

diag_2pm.o                  : diag_2pm.f                                     \
                              diag_mod.o              tropopause_mod.o


diag_mod.o                  : diag_mod.f

diag_oh_mod.o               : diag_oh_mod.f                                  \
                              comode_mod.o            logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

diag_pl_mod.o               : diag_pl_mod.f                                  \
                              comode_mod.o            logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

drydep_mod.o                : drydep_mod.f            comode_mod.o          \
                              dao_mod.o               diag_mod.o            \
                              logical_mod.o           pbl_mix_mod.o         \
                              tracerid_mod.o          tracer_mod.o          \
                              meganut_mod.o 

diagoh.o                    : diagoh.f                diag_mod.o

dust_dead_mod.o             : dust_dead_mod.f         dao_mod.o              

dust_mod.o                  : dust_mod.f                                     \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              drydep_mod.o           \
                              dust_dead_mod.o         logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o    

edgar_mod.o                 : edgar_mod.f                                    \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o    

EFOLD.o                     : EFOLD.f

emep_mod.o                  : emep_mod.f                                     \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o  

emf_scale.o                 : emf_scale.f             tracerid_mod.o

emfossil.o                  : emfossil.f                                     \
                              bravo_mod.o             cac_anthro_mod.o       \
                              dao_mod.o               diag_mod.o             \
                              edgar_mod.o             emep_mod.o             \
                              epa_nei_mod.o           icoads_ship_mod.o      \
                              logical_mod.o           nei2005_anthro_mod.o   \
                              streets_anthro_mod.o    tracer_mod.o           \
                              tracerid_mod.o          vistas_anthro_mod.o    

emisop.o                    : emisop.f                dao_mod.o              \
                              meganut_mod.o

emisop_grass.o              : emisop_grass.f          dao_mod.o

emisop_mb.o                 : emisop_mb.f             dao_mod.o

emissdr.o                   : emissdr.f                                      \
                              acetone_mod.o           aircraft_nox_mod.o     \
                              biofuel_mod.o           dao_mod.o              \
                              diag_mod.o              emissions_mod.o        \
			      lightning_nox_mod.o     logical_mod.o          \
                              megan_mod.o             tracer_mod.o           \
                              tracerid_mod.o          meganut_mod.o

emissions_mod.o             : emissions_mod.f                                \
                              arctas_ship_emiss_mod.o biomass_mod.o          \
                              bravo_mod.o             c2h6_mod.o             \
                              cac_anthro_mod.o        carbon_mod.o           \
                              ch3i_mod.o              co2_mod.o              \
                              dust_mod.o              edgar_mod.o            \
                              emep_mod.o              epa_nei_mod.o          \
                              global_ch4_mod.o        h2_hd_mod.o            \
                              hcn_ch3cn_mod.o         icoads_ship_mod.o      \
                              logical_mod.o           mercury_mod.o          \
                              nei2005_anthro_mod.o    RnPbBe_mod.o           \
                              seasalt_mod.o           streets_anthro_mod.o   \
                              sulfate_mod.o           tagged_co_mod.o        \
                              tracer_mod.o            vistas_anthro_mod.o

emmonot.o                   : emmonot.f

epa_nei_mod.o               : epa_nei_mod.f                                  \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o    

fast_j.o                    : fast_j.f                                       \
                              dao_mod.o               toms_mod.o

fertadd.o                   : fertadd.f               logical_mod.o

findmon.o                   : findmon.f

fjfunc.o                    : fjfunc.f

fjx_acet_mod.o              : fjx_acet_mod.f

future_emissions_mod.o      : future_emissions_mod.f

fvdas_convect_mod.o         : fvdas_convect_mod.f                            \
                              dao_mod.o               diag_mod.o             \
                              depo_mercury_mod.o      logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

fyhoro.o                    : fyhoro.f

fyrno3.o                    : fyrno3.f

gamap_mod.o                 : gamap_mod.f                                    \
                              diag03_mod.o            diag04_mod.o           \
                              diag41_mod.o            diag42_mod.o           \
                              diag48_mod.o            diag49_mod.o           \
                              diag50_mod.o            diag51_mod.o           \
                              diag51b_mod.o           diag56_mod.o           \
                              diag_pl_mod.o  	      drydep_mod.o           \
                              logical_mod.o    	      tracerid_mod.o         \
                              tracer_mod.o            wetscav_mod.o  
ifeq ($(COMPILER),sun)
	$(F90) -O0 -c $<
endif

gasconc.o                   : gasconc.f                                      \
                              comode_mod.o            dao_mod.o              \
                              drydep_mod.o            tropopause_mod.o       \
                              logical_mod.o

GAUSSP.o                    : GAUSSP.f

gc_biomass_mod.o            : gc_biomass_mod.f                               \
                              future_emissions_mod.o  logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o    

gcap_convect_mod.o          : gcap_convect_mod.f                             \
                              dao_mod.o               diag_mod.o

gcap_read_mod.o             : gcap_read_mod.f         dao_mod.o              \
                              diag_mod.o              logical_mod.o        

geia_mod.o                  : geia_mod.f

GEN.o                       : GEN.f

get_global_ch4.o            : get_global_ch4.f                               \
                              future_emissions_mod.o  logical_mod.o

getifsun.o                  : getifsun.f              comode_mod.o

gfed2_biomass_mod.o         : gfed2_biomass_mod.f                            \
                              future_emissions_mod.o  logical_mod.o          \
                              tracer_mod.o            tracerid_mod.o

global_br_mod.o             : global_br_mod.f         tropopause_mod.o 

global_ch4_mod.o            : global_ch4_mod.f                               \
                              dao_mod.o               diag_mod.o             \
                              diag_oh_mod.o           diag_pl_mod.o          \
                              global_oh_mod.o         logical_mod.o          \
                              tracer_mod.o            vdiff_pre_mod.o  

global_hno3_mod.o           : global_hno3_mod.f                              \
                              dao_mod.o               tracer_mod.o  

global_no3_mod.o            : global_no3_mod.f

global_nox_mod.o            : global_nox_mod.f

global_o1d_mod.o            : global_o1d_mod.f

global_o3_mod.o             : global_o3_mod.f

global_oh_mod.o             : global_oh_mod.f

gwet_read_mod.o             : gwet_read_mod.f                                \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o 

h2_hd_mod.o                 : h2_hd_mod.f                                    \
                              biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            geia_mod.o             \
                              global_nox_mod.o        global_o1d_mod.o       \
                              global_oh_mod.o         logical_mod.o    	     \
                              scale_anthro_mod.o      tagged_co_mod.o        \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o        meganut_mod.o

hcn_ch3cn_mod.o             : hcn_ch3cn_mod.f                                \
                              biomass_mod.o           dao_mod.o              \
                              diag_mod.o              geia_mod.o             \
                              global_oh_mod.o         logical_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o

i6_read_mod.o               : i6_read_mod.f                                  \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           

icoads_ship_mod.o           : icoads_ship_mod.f                              \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o  

initialize.o                : initialize.f                                   \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o            diag56_mod.o           \
                              diag_pl_mod.o           logical_mod.o 

inphot.o                    : inphot.f

input_mod.o                 : input_mod.f                                    \
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
                              transport_mod.o         upbdflx_mod.o          \
                              wetscav_mod.o

isoropiaII_mod.o            : isoropiaII_mod.f                               \
                              dao_mod.o               global_hno3_mod.o      \
                              logical_mod.o           tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o

jsparse.o                   : jsparse.f               comode_mod.o

JRATET.o                    : JRATET.f                fjx_acet_mod.o

JVALUE.o                    : JVALUE.f

jv_index.o                  : jv_index.f

ksparse.o                   : ksparse.f

lai_mod.o                   : lai_mod.f               logical_mod.o

land_mercury_mod.o          : land_mercury_mod.f                             \
                              biomass_mod.o           dao_mod.o              \
                              depo_mercury_mod.o      lai_mod.o              \
                              logical_mod.o           tracerid_mod.o         \

LEGND0.o                    : LEGND0.f

lightning_nox_mod.o         : lightning_nox_mod.f                            \
                              dao_mod.o               diag56_mod.o           \
                              diag_mod.o              logical_mod.o        

linoz_mod.o                 : linoz_mod.f                                    \
                              dao_mod.o               tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o  

logical_mod.o               : logical_mod.f

lump.o                      : lump.f                                         \
                              comode_mod.o            tracerid_mod.o

main.o                      : main.f                                         \
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
                              gcap_read_mod.o         gwet_read_mod.o        \
                              i6_read_mod.o           input_mod.o            \
                              lai_mod.o               lightning_nox_mod.o    \
                              linoz_mod.o             logical_mod.o          \
                              megan_mod.o             pbl_mix_mod.o          \
                              ocean_mercury_mod.o     planeflight_mod.o      \
                              soaprod_mod.o           tpcore_bc_mod.o        \
                              tracer_mod.o            transport_mod.o        \
                              tropopause_mod.o        restart_mod.o          \
                              upbdflx_mod.o           uvalbedo_mod.o         \
                              vdiff_mod.o             wetscav_mod.o          \
                              xtra_read_mod.o

MATIN4.o                    : MATIN4.f

megan_mod.o                 : megan_mod.f             a3_read_mod.o          \
                              lai_mod.o               logical_mod.o          \
                              meganut_mod.o

meganut_mod.o               : meganut_mod.f           dao_mod.o     

mercury_mod.o               : mercury_mod.f                                  \
                              dao_mod.o               depo_mercury_mod.o     \
                              diag03_mod.o            diag_mod.o             \
                              drydep_mod.o            global_br_mod.o        \
                              global_o3_mod.o         global_oh_mod.o        \
                              land_mercury_mod.o      lai_mod.o              \
                              logical_mod.o           ocean_mercury_mod.o    \
                              pbl_mix_mod.o           RnPbBe_mod.o           \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o        vdiff_pre_mod.o

MIESCT.o                    : MIESCT.f

ndxx_setup.o                : ndxx_setup.f                                   \
                              biofuel_mod.o           diag_mod.o             \
                              diag_oh_mod.o           drydep_mod.o           \
                              logical_mod.o           planeflight_mod.o      \
                              tracer_mod.o            tracerid_mod.o         \
                              wetscav_mod.o     

nei2005_anthro_mod.o        : nei2005_anthro_mod.f                           \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o    

NOABS.o                     : NOABS.f

ocean_mercury_mod.o         : ocean_mercury_mod.f                            \
                              dao_mod.o               depo_mercury_mod.o     \
                              diag03_mod.o            logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o    

ohsave.o                    : ohsave.f                comode_mod.o           \
                              diag_mod.o              tracerid_mod.o

optdepth_mod.o              : optdepth_mod.f          diag_mod.o  

OPMIE.o                     : OPMIE.f

partition.o                 : partition.f                                    \
                              comode_mod.o            tracerid_mod.o

pbl_mix_mod.o               : pbl_mix_mod.f                                  \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracer_mod.o  

pderiv.o                    : pderiv.f

photoj.o                    : photoj.f

physproc.o                  : physproc.f                                     \
                              comode_mod.o            logical_mod.o          \
                              chemistry_mod.o

pjc_pfix_geos5_window_mod.o : pjc_pfix_geos5_window_mod.f

pjc_pfix_mod.o              : pjc_pfix_mod.f

planeflight_mod.o           : planeflight_mod.f                              \
                              comode_mod.o            dao_mod.o              \
                              tracer_mod.o            tropopause_mod.o  

precipfrac.o                : precipfrac.f            dao_mod.o

pulsing.o                   : pulsing.f

RD_AOD.o                    : RD_AOD.f

rd_js.o                     : rd_js.f

rd_prof.o                   : rd_prof.f

rdisopt.o                   : rdisopt.f

rdlai.o                     : rdlai.f

rdland.o                    : rdland.f                logical_mod.o

rdlight.o                   : rdlight.f

rdmonot.o                   : rdmonot.f

rdsoil.o                    : rdsoil.f

RD_TJPL.o                   : RD_TJPL.f

readchem.o                  : readchem.f                                     \
                              diag_pl_mod.o           drydep_mod.o           \
                              logical_mod.o

reader.o                    : reader.f

readlai.o                   : readlai.f               logical_mod.o

restart_mod.o               : restart_mod.f                                  \
                              comode_mod.o            dao_mod.o              \
                              logical_mod.o           tracer_mod.o  

RnPbBe_mod.o                : RnPbBe_mod.f                                   \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracer_mod.o           \
                              tropopause_mod.o  		       

rpmares_mod.o               : rpmares_mod.f                                  \
                              dao_mod.o               global_hno3_mod.o      \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o        

ruralbox.o                  : ruralbox.f                                     \
                              comode_mod.o            tropopause_mod.o

scale_anthro_mod.o          : scale_anthro_mod.f

schem.o                     : schem.f                 dao_mod.o              \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o 

seasalt_mod.o               : seasalt_mod.f                                  \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            logical_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            vdiff_pre_mod.o

set_aer.o                   : set_aer.f

set_prof.o                  : set_prof.f                                     \
                              dao_mod.o               toms_mod.o

setbase.o                   : setbase.f

setemdep.o                  : setemdep.f              drydep_mod.o           \
                              tracer_mod.o            tracerid_mod.o

setemis.o                   : setemis.f                                      \
                              aircraft_nox_mod.o      biofuel_mod.o          \
                              biomass_mod.o           comode_mod.o           \
                              diag_mod.o              emissions_mod.o        \
                              lightning_nox_mod.o     logical_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tropopause_mod.o    

setmodel.o                  : setmodel.f

sfcwindsqr.o                : sfcwindsqr.f            dao_mod.o

smvgear.o                   : smvgear.f               comode_mod.o

soaprod_mod.o               : soaprod_mod.f                                  \
                              carbon_mod.o            dao_mod.o              \
                              logical_mod.o           tracer_mod.o

soilbase.o                  : soilbase.f

soilcrf.o                   : soilcrf.f

soilnoxems.o                : soilnoxems.f                                   \
                              dao_mod.o               diag_mod.o             \
                              future_emissions_mod.o  logical_mod.o          \
                              meganut_mod.o

soiltemp.o                  : soiltemp.f

soiltype.o                  : soiltype.f

SPHERE.o                    : SPHERE.f

streets_anthro_mod.o        : streets_anthro_mod.f                           \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o  

subfun.o                    : subfun.f

sulfate_mod.o               : sulfate_mod.f                                  \
                              arctas_ship_emiss_mod.o biomass_mod.o          \
                              bravo_mod.o             cac_anthro_mod.o       \
                              comode_mod.o            dao_mod.o              \
                              diag_mod.o              drydep_mod.o           \
                              edgar_mod.o             emep_mod.o             \
                              epa_nei_mod.o           future_emissions_mod.o \
                              gfed2_biomass_mod.o     global_hno3_mod.o      \
                              global_no3_mod.o        global_oh_mod.o        \
                              icoads_ship_mod.o       logical_mod.o          \
                              nei2005_anthro_mod.o    pbl_mix_mod.o          \
                              scale_anthro_mod.o      seasalt_mod.o          \
                              streets_anthro_mod.o    tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o       \
                              uvalbedo_mod.o          vdiff_pre_mod.o        \
                              wetscav_mod.o

sunparam.o                  : sunparam.f

tagged_co_mod.o             : tagged_co_mod.f                                \
                              biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o             \
                              diag_pl_mod.o           global_nox_mod.o       \
                              global_oh_mod.o         logical_mod.o          \
                              megan_mod.o             meganut_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o 

tagged_ox_mod.o             : tagged_ox_mod.f                                \
                              dao_mod.o               diag_mod.o             \
                              diag_pl_mod.o           drydep_mod.o           \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o        meganut_mod.o

tcorr.o                     : tcorr.f

toms_mod.o                  : toms_mod.f

tpcore_bc_mod.o             : tpcore_bc_mod.f                                \
                              logical_mod.o           tracer_mod.o  

tpcore_geos5_window_mod.o   : tpcore_geos5_window_mod.f90
	$(F90) -c $(FREEFORM) $(R8) $<

tpcore_mod.o                : tpcore_mod.f            dao_mod.o              \
                              diag_mod.o              global_ch4_mod.o       \
                              tracer_mod.o
	$(F90) -c $(R8) $<

tpcore_window_mod.o         : tpcore_window_mod.f     dao_mod.o              \
                              diag_mod.o              global_ch4_mod.o       \
                              tracer_mod.o
	$(F90) -c $(R8) $<

tracer_mod.o                : tracer_mod.f

tracerid_mod.o              : tracerid_mod.f                                 \
                              logical_mod.o           tracer_mod.o

transport_mod.o             : transport_mod.f                                \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           pjc_pfix_mod.o         \
                              tpcore_mod.o            tpcore_bc_mod.o        \
                              tpcore_fvdas_mod.o      tpcore_window_mod.o    \
                              tracer_mod.o                                   \
                              tpcore_geos5_window_mod.o                      \
                              pjc_pfix_geos5_window_mod.o


tropopause.o                : tropopause.f                                   \
                              dao_mod.o               diag_mod.o             \
                              logical_mod.o           tropopause_mod.o 

tropopause_mod.o            : tropopause_mod.f                               \
                              comode_mod.o            dao_mod.o              \
                              logical_mod.o

upbdflx_mod.o               : upbdflx_mod.f                                  \
                              dao_mod.o               linoz_mod.o            \
                              logical_mod.o           tagged_ox_mod.o        \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o

update.o                    : update.f

vdiff_pre_mod.o             : vdiff_pre_mod.f         tracer_mod.o
ifeq ($(COMPILER),sun)
	$(F90) -O3 -c $<
endif

vdiff_mod.o                 : vdiff_mod.f90                                  \
                              comode_mod.o            dao_mod.o              \
                              depo_mercury_mod.o      diag_mod.o             \
                              drydep_mod.o            logical_mod.o          \
                              ocean_mercury_mod.o     pbl_mix_mod.o          \
                              tracer_mod.o            tracerid_mod.o         \
                              vdiff_pre_mod.o

vistas_anthro_mod.o         : vistas_anthro_mod.f                            \
                              future_emissions_mod.o  epa_nei_mod.o          \
                              logical_mod.o           scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o  

wetscav_mod.o               : wetscav_mod.f                                  \
                              dao_mod.o               diag_mod.o             \
                              depo_mercury_mod.o      logical_mod.o          \
                              mercury_mod.o           ocean_mercury_mod.o    \
                              tracerid_mod.o          tracer_mod.o

XSEC1D.o                    : XSEC1D.f

XSECO2.o                    : XSECO2.f

XSECO3.o                    : XSECO3.f

xtra_read_mod.o             : xtra_read_mod.f         dao_mod.o             \
                              diag_mod.o              logical_mod.o

#EOC
