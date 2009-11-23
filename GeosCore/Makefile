# $Id: Makefile,v 1.12 2009/11/23 17:15:03 bmy Exp $
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
#EOP
#------------------------------------------------------------------------------
#BOC

# Define variables
SHELL   = /bin/sh
ROOTDIR = ..
BIN     = $(ROOTDIR)/bin
BPCH    = $(ROOTDIR)/GeosBpch
DOC     = $(ROOTDIR)/doc
HDR     = $(ROOTDIR)/Headers
HELP    = $(ROOTDIR)/help
EXE     = geos
LIB     = $(ROOTDIR)/lib
KPP     = $(ROOTDIR)/KPP
MOD     = $(ROOTDIR)/mod
UTIL    = $(ROOTDIR)/GeosUtil

# Link to library files created from code in the various subdirs
# NOTE: -lGeosCore should always be first, -lGeosUtil always last!
LINK    = -lGeosCore -lGeosBpch -lKpp -lGeosUtil

# Include header file.  This returns variables CC, F90, FREEFORM, LD, R8,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

#=============================================================================
# List of files to compile (the order is important!).  We specify these as
# a list of object files (*.o).  For each object file, the "make" utility
# will find the corresponding source code file (*.f, *.f90) and compile it. 
#=============================================================================

# Source and object files
SRC   = $(wildcard *.f*)
TMP   = $(SRC:.f=.o)
OBJ   = $(TMP:.f90=.o)

#=============================================================================
# Makefile targets: type "make help" for a complete listing!
#=============================================================================

.PHONY: all clean realclean doc docclean help

all: libkpp libutil libbpch lib exe

lib: $(OBJ)
	cp *.mod $(MOD)
	$(AR) crs libGeosCore.a $(OBJ)
	mv libGeosCore.a $(LIB)

libbpch:
	@$(MAKE) -C $(BPCH)

libkpp: 
	@$(MAKE) -C $(KPP)

libutil:
	@$(MAKE) -C $(UTIL)

exe: main.o
	$(LD) main.o $(LINK) -o $(EXE)
	cp -f $(EXE) $(BIN)

clean:
	rm -f *.o *.mod geos

realclean:
	@$(MAKE) clean
	@$(MAKE) -C $(KPP) realclean
	@$(MAKE) -C $(BPCH) clean
	@$(MAKE) -C $(UTIL) clean
	@$(MAKE) docclean
	rm -f $(LIB)/*.a
	rm -f $(MOD)/*.mod

doc:
	@$(MAKE) -C $(DOC) all

docclean: 
	@$(MAKE) -C $(DOC) clean

help:
	@$(MAKE) -C $(HELP)

print:
	@echo $(OBJ)

#==============================================================================
# Dependencies listing (grep "USE " to get the list of module references!)
#
# From this list of dependencies, the "make" utility will figure out the
# correct order of compilation (so we don't have to do that ourselves).
# This also allows us to compile on multiple processors with "make -j".
#
# NOTE: Only specify object files that are within this directory.  
# Object files in other directories will be referenced at link-time.
#==============================================================================

global_ch4_mod.o            : dao_mod.o               diag_mod.o             \
                              diag_oh_mod.o           diag_pl_mod.o          \
                              logical_mod.o           tracer_mod.o           \
                              vdiff_pre_mod.o  		       

RnPbBe_mod.o                : dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracer_mod.o           \
                              tropopause_mod.o  		       

a3_read_mod.o               : dao_mod.o               diag_mod.o             \
                              logical_mod.o        		       

a6_read_mod.o               : dao_mod.o               diag_mod.o             \
                              logical_mod.o        		       

acetone_mod.o               : dao_mod.o               diag_mod.o               

aerosol_mod.o               : comode_mod.o            dao_mod.o              \
                              diag_mod.o              logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o  

aircraft_nox_mod.o          : dao_mod.o               diag_mod.o              

anthroems.o                 : future_emissions_mod.o  edgar_mod.o            \
                              geia_mod.o              logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o

arctas_ship_emiss_mod.o     : logical_mod.o           scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o  

benchmark_mod.o             : tracerid_mod.o          tracer_mod.o  

biofuel_mod.o               : dao_mod.o               diag_mod.o             \
                              epa_nei_mod.o           future_emissions_mod.o \
                              logical_mod.o           streets_anthro_mod.o   \
                              tracerid_mod.o          tracer_mod.o    

biomass_mod.o               : diag_mod.o              gc_biomass_mod.o       \
                              gfed2_biomass_mod.o     logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

boxvl.o                     : dao_mod.o

bravo_mod.o                 : future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o    

c2h6_mod.o                  : biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o             \
                              geia_mod.o              global_oh_mod.o        \
                              logical_mod.o           tracerid_mod.o         \
                              tracer_mod.o    

cac_anthro_mod.o            : future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o  

calcrate.o                  : comode_mod.o            dao_mod.o              \
                              diag_mod.o              drydep_mod.o           \
                              logical_mod.o           pbl_mix_mod.o          \
                              planeflight_mod.o       tracerid_mod.o

carbon_mod.o                : biomass_mod.o           comode_mod.o           \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            future_emissions_mod.o \
                              gfed2_biomass_mod.o     global_no3_mod.o       \
                              global_o3_mod.o         global_oh_mod.o        \
                              logical_mod.o           megan_mod.o            \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o       \
                              vdiff_pre_mod.o  

ch3i_mod.o                  : biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o             \
                              diag_pl_mod.o           logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              uvalbedo_mod.o

chemdr.o                    : aerosol_mod.o           comode_mod.o           \
                              chemistry_mod.o         dao_mod.o              \
                              diag_oh_mod.o           diag_pl_mod.o          \
                              future_emissions_mod.o  dust_mod.o             \
                              logical_mod.o           planeflight_mod.o      \
                              restart_mod.o           tracer_mod.o           \
                              tracerid_mod.o          tropopause_mod.o       \
                              uvalbedo_mod.o   

chemistry_mod.o             : acetone_mod.o           aerosol_mod.o          \
                              c2h6_mod.o              carbon_mod.o           \
                              ch3i_mod.o              comode_mod.o           \
                              dao_mod.o               drydep_mod.o           \
                              dust_mod.o              global_ch4_mod.o       \
                              h2_hd_mod.o             hcn_ch3cn_mod.o        \
	                      isoropia_mod.o                                 \
	                      logical_mod.o           mercury_mod.o          \
	                      optdepth_mod.o          rpmares_mod.o          \
                              RnPbBe_mod.o            seasalt_mod.o          \
                              sulfate_mod.o           tagged_co_mod.o        \
                              tagged_ox_mod.o         tracerid_mod.o         \
                              tracer_mod.o

cleanup.o                   : acetone_mod.o           aerosol_mod.o          \
                              arctas_ship_emiss_mod.o aircraft_nox_mod.o     \
                              biomass_mod.o           biofuel_mod.o          \
                              bravo_mod.o             c2h6_mod.o	     \
                              cac_anthro_mod.o        carbon_mod.o  	     \
                              co2_mod.o               comode_mod.o 	     \
                              dao_mod.o                	                     \
                              diag_mod.o              diag03_mod.o 	     \
                              diag04_mod.o            diag41_mod.o	     \
                              diag50_mod.o            diag51_mod.o	     \
                              diag_oh_mod.o           diag_pl_mod.o	     \
                              drydep_mod.o            dust_mod.o   	     \
                              dust_dead_mod.o         edgar_mod.o	     \
                              emep_mod.o              epa_nei_mod.o	     \
                              gc_biomass_mod.o                               \
                              gfed2_biomass_mod.o     global_ch4_mod.o       \
                              global_hno3_mod.o       global_no3_mod.o       \
                              global_nox_mod.o        global_o1d_mod.o       \
                              global_oh_mod.o         h2_hd_mod.o	     \
                              hcn_ch3cn_mod.o         lai_mod.o	             \
                              lightning_nox_mod.o     linoz_mod.o	     \
                              megan_mod.o             mercury_mod.o	     \
                              ocean_mercury_mod.o     pbl_mix_mod.o	     \
                              pjc_pfix_mod.o          planeflight_mod.o      \
                              seasalt_mod.o           sulfate_mod.o	     \
                              tagged_co_mod.o         toms_mod.o	     \
                              tracer_mod.o            transport_mod.o        \
                              tropopause_mod.o        uvalbedo_mod.o         \
                              vistas_anthro_mod.o     wetscav_mod.o          \
                              icoads_ship_mod.o       tpcore_fvdas_mod.o     \
                              tpcore_geos5_window_mod.o


co2_mod.o                   : biomass_mod.o           diag04_mod.o           \
                              tracer_mod.o    

CO_strat_pl.o               : dao_mod.o               tracer_mod.o           \
                              tracerid_mod.o          tropopause_mod.o 

comode_mod.o                : tracer_mod.o

convection_mod.o            : dao_mod.o               diag_mod.o             \
                              fvdas_convect_mod.o     gcap_convect_mod.o     \
                              logical_mod.o           ocean_mercury_mod.o    \
                              tracer_mod.o            tracerid_mod.o         \
                              wetscav_mod.o

dao_mod.o                   : logical_mod.o           tracer_mod.o           \
                              diag03_mod.o            tracerid_mod.o  

diag03_mod.o                : tracerid_mod.o

diag_2pm.o                  : tropopause_mod.o

diag1.o                     : diag_mod.o              tracer_mod.o           \
                              tracerid_mod.o          tropopause_mod.o

diag3.o                     : biofuel_mod.o           biomass_mod.o          \
                              diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o            diag56_mod.o           \
                              diag_pl_mod.o           drydep_mod.o           \
                              logical_mod.o           tracer_mod.o           \
                              tracerid_mod.o          wetscav_mod.o

diag41_mod.o                : pbl_mix_mod.o    

diag42_mod.o                : dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracerid_mod.o         \
                              tracer_mod.o  

diag48_mod.o                : dao_mod.o               pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

diag49_mod.o                : dao_mod.o               pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

diag50_mod.o                : comode_mod.o            dao_mod.o              \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o  

diag51_mod.o                : dao_mod.o               pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o  

diag_oh_mod.o               : comode_mod.o            logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o  

diag_pl_mod.o               : comode_mod.o            tracerid_mod.o         \
                              tracer_mod.o  

drydep_mod.o                : comode_mod.o            dao_mod.o              \
                              diag_mod.o              logical_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o  

dust_dead_mod.o             : dao_mod.o              

dust_mod.o                  : comode_mod.o            dao_mod.o              \
                              diag_mod.o              drydep_mod.o           \
                              dust_dead_mod.o         logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o    

edgar_mod.o                 : future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o    

emep_mod.o                  : future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o  


emf_scale.o                 : tracerid_mod.o

emfossil.o                  : bravo_mod.o             cac_anthro_mod.o       \
                              dao_mod.o               diag_mod.o             \
                              edgar_mod.o             emep_mod.o             \
                              epa_nei_mod.o           icoads_ship_mod.o      \
                              logical_mod.o           nei2005_anthro_mod.o   \
                              streets_anthro_mod.o    tracer_mod.o           \
                              tracerid_mod.o          vistas_anthro_mod.o    

emisop.o                    : dao_mod.o

emisop_grass.o              : dao_mod.o

emisop_mb.o                 : dao_mod.o

emissdr.o                   : acetone_mod.o           aircraft_nox_mod.o     \
                              biofuel_mod.o           dao_mod.o              \
                              diag_mod.o              lightning_nox_mod.o    \
                              logical_mod.o           megan_mod.o            \
                              tracer_mod.o            tracerid_mod.o  

emissions_mod.o             : arctas_ship_emiss_mod.o biomass_mod.o          \
                              bravo_mod.o             c2h6_mod.o             \
                              cac_anthro_mod.o        carbon_mod.o           \
                              ch3i_mod.o              co2_mod.o              \
                              dust_mod.o              edgar_mod.o            \
                              emep_mod.o              epa_nei_mod.o          \
                              global_ch4_mod.o        h2_hd_mod.o            \
                              hcn_ch3cn_mod.o         icoads_ship_mod.o      \
                              logical_mod.o                                  \
                              mercury_mod.o           nei2005_anthro_mod.o   \
                              RnPbBe_mod.o            seasalt_mod.o          \
                              streets_anthro_mod.o    sulfate_mod.o          \
                              tagged_co_mod.o         tracer_mod.o           \
                              vistas_anthro_mod.o  

epa_nei_mod.o               : future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o    

fast_j.o                    : dao_mod.o               toms_mod.o

fvdas_convect_mod.o         : diag_mod.o              logical_mod.o          \
                              ocean_mercury_mod.o     tracerid_mod.o         \
                              tracer_mod.o  

gamap_mod.o                 : diag03_mod.o            diag04_mod.o           \
                              diag41_mod.o            diag42_mod.o           \
                              diag48_mod.o            diag49_mod.o  	     \
                              diag50_mod.o            diag51_mod.o  	     \
                              diag56_mod.o            diag_pl_mod.o  	     \
                              drydep_mod.o            logical_mod.o    	     \
                              tracerid_mod.o          tracer_mod.o           \
                              wetscav_mod.o  

gasconc.o                   : comode_mod.o            dao_mod.o              \
                              drydep_mod.o            tropopause_mod.o       \
                              logical_mod.o

gc_biomass_mod.o            : future_emissions_mod.o  logical_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              gcap_convect_mod.o      dao_mod.o              \
                              diag_mod.o      

gcap_convect_mod.o          : dao_mod.o               diag_mod.o

gcap_read_mod.o             : dao_mod.o               diag_mod.o             \
                              logical_mod.o        

getifsun.o                  : comode_mod.o

gfed2_biomass_mod.o         : future_emissions_mod.o  logical_mod.o          \
                              tracer_mod.o  

global_ch4_mod.o            : dao_mod.o               diag_mod.o             \
                              diag_oh_mod.o           diag_pl_mod.o          \
                              logical_mod.o           tracer_mod.o           \
                              vdiff_pre_mod.o  

global_hno3_mod.o           : dao_mod.o               tracer_mod.o  

gwet_read_mod.o             : dao_mod.o               diag_mod.o             \
                              logical_mod.o           h2_hd_mod.o 

h2_hd_mod.o                 : biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o             \
                              drydep_mod.o            geia_mod.o  	     \
                              global_nox_mod.o        global_o1d_mod.o       \
                              global_oh_mod.o         logical_mod.o    	     \
                              scale_anthro_mod.o      tagged_co_mod.o        \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o  

hcn_ch3cn_mod.o             : biomass_mod.o           dao_mod.o              \
                              diag_mod.o              geia_mod.o             \
                              global_oh_mod.o         logical_mod.o          \
                              pbl_mix_mod.o      

i6_read_mod.o               : dao_mod.o               diag_mod.o             \
                              logical_mod.o           icoads_ship_mod.o      \
                              future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o  

icoads_ship_mod.o           : future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o  

initialize.o                : diag_mod.o              diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o           \
                              diag42_mod.o            diag56_mod.o           \
                              diag_pl_mod.o           logical_mod.o 

input_mod.o                 : benchmark_mod.o         biofuel_mod.o          \
                              biomass_mod.o           diag03_mod.o           \
                              diag04_mod.o            diag41_mod.o  	     \
                              diag42_mod.o            diag48_mod.o  	     \
                              diag49_mod.o            diag50_mod.o  	     \
                              diag51_mod.o            diag56_mod.o  	     \
                              diag_oh_mod.o           diag_pl_mod.o 	     \
                              drydep_mod.o            future_emissions_mod.o \
                              gamap_mod.o             logical_mod.o 	     \
                              mercury_mod.o           ocean_mercury_mod.o    \
                              planeflight_mod.o       restart_mod.o 	     \
                              tpcore_bc_mod.o         tracerid_mod.o	     \
                              tracer_mod.o            transport_mod.o	     \
                              upbdflx_mod.o           wetscav_mod.o

isoropia_mod.o              : dao_mod.o               global_hno3_mod.o      \
                              logical_mod.o           tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o	     \
                              diag04_mod.o            diag41_mod.o  	     \
                              diag42_mod.o            diag48_mod.o  	     \
                              diag49_mod.o            diag50_mod.o  	     \
                              diag51_mod.o            diag56_mod.o  	     \
                              diag_oh_mod.o           diag_pl_mod.o    	     \
                              drydep_mod.o            future_emissions_mod.o \
                              gamap_mod.o             logical_mod.o  	     \
                              mercury_mod.o           ocean_mercury_mod.o    \
                              planeflight_mod.o       restart_mod.o    	     \
                              tpcore_bc_mod.o         tracerid_mod.o  	     \
                              tracer_mod.o            transport_mod.o        \
                              upbdflx_mod.o           wetscav_mod.o 

jsparse.o                   : comode_mod.o

JRATET.o                    : fjx_acet_mod.o

lightning_nox_mod.o         : dao_mod.o               diag56_mod.o           \
                              diag_mod.o              logical_mod.o        

linoz_mod.o                 : dao_mod.o               tracerid_mod.o         \
                              tracer_mod.o            tropopause_mod.o  

lump.o                      : comode_mod.o            tracerid_mod.o

main.o                      : a3_read_mod.o           a6_read_mod.o          \
                              benchmark_mod.o         carbon_mod.o           \
                              chemistry_mod.o         convection_mod.o       \
                              comode_mod.o                                   \
                              diag_mod.o              diag41_mod.o           \
                              diag42_mod.o            diag48_mod.o           \
                              diag49_mod.o            diag50_mod.o           \
                              diag51_mod.o            diag_oh_mod.o          \
                              dao_mod.o               drydep_mod.o           \
                              emissions_mod.o         global_ch4_mod.o       \
                              gcap_read_mod.o         gwet_read_mod.o        \
                              i6_read_mod.o           input_mod.o            \
                              lai_mod.o               lightning_nox_mod.o    \
                              linoz_mod.o             logical_mod.o          \
                              megan_mod.o             pbl_mix_mod.o          \
                              ocean_mercury_mod.o     planeflight_mod.o      \
                              tracer_mod.o            transport_mod.o        \
                              tropopause_mod.o        restart_mod.o          \
                              upbdflx_mod.o           uvalbedo_mod.o         \
                              wetscav_mod.o           xtra_read_mod.o        \
                              vdiff_mod.o        

megan_mod.o                 : a3_read_mod.o           lai_mod.o              \
                              logical_mod.o      

mercury_mod.o               : dao_mod.o               diag03_mod.o           \
                              diag_mod.o              drydep_mod.o   	     \
                              global_o3_mod.o         global_oh_mod.o        \
                              logical_mod.o           ocean_mercury_mod.o    \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o    

ndxx_setup.o                : biomass_mod.o           biofuel_mod.o          \
                              diag_mod.o              diag_oh_mod.o          \
                              drydep_mod.o            logical_mod.o          \
                              planeflight_mod.o       tracer_mod.o           \
                              tracerid_mod.o          wetscav_mod.o     

nei2005_anthro_mod.o        : future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o    

ocean_mercury_mod.o         : dao_mod.o               diag03_mod.o           \
                              logical_mod.o           tracerid_mod.o         \
                              tracer_mod.o    

ohsave.o                    : comode_mod.o            diag_mod.o             \
                              tracerid_mod.o

optdepth_mod.o              : diag_mod.o  

partition.o                 : comode_mod.o            tracerid_mod.o

pbl_mix_mod.o               : dao_mod.o               diag_mod.o             \
                              logical_mod.o           tracer_mod.o  

physproc.o                  : comode_mod.o            logical_mod.o

planeflight_mod.o           : comode_mod.o            dao_mod.o              \
                              tracer_mod.o            tropopause_mod.o  

precipfrac.o                : dao_mod.o

readchem.o                  : diag_pl_mod.o           drydep_mod.o

readlai.o                   : logical_mod.o

restart_mod.o               : comode_mod.o            dao_mod.o              \
                              logical_mod.o           tracer_mod.o  

rpmares_mod.o               : dao_mod.o               global_hno3_mod.o      \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o        

ruralbox.o                  : comode_mod.o            tropopause_mod.o


schem.o                     : dao_mod.o               tracer_mod.o           \
                              tropopause_mod.o 

seasalt_mod.o               : dao_mod.o               diag_mod.o             \
                              drydep_mod.o            logical_mod.o          \
                              pbl_mix_mod.o           tracerid_mod.o         \
                              tracer_mod.o            vdiff_pre_mod.o

set_prof.o                  : dao_mod.o               toms_mod.o

setemdep.o                  : drydep_mod.o            tracer_mod.o           \
                              tracerid_mod.o

setemis.o                   : aircraft_nox_mod.o      biofuel_mod.o          \
                              biomass_mod.o           comode_mod.o           \
                              diag_mod.o              logical_mod.o          \
                              lightning_nox_mod.o     pbl_mix_mod.o          \
                              tracerid_mod.o          tropopause_mod.o    

sfcwindsqr.o                : dao_mod.o

smvgear.o                   : comode_mod.o

soilnoxems.o                : dao_mod.o               diag_mod.o             \
                              future_emissions_mod.o  logical_mod.o


streets_anthro_mod.o        : future_emissions_mod.o  logical_mod.o          \
                              scale_anthro_mod.o      tracerid_mod.o         \
                              tracer_mod.o  

sulfate_mod.o               : arctas_ship_emiss_mod.o biomass_mod.o          \
                              bravo_mod.o             cac_anthro_mod.o       \
                              comode_mod.o            dao_mod.o  	     \
                              diag_mod.o              drydep_mod.o  	     \
                              edgar_mod.o             emep_mod.o  	     \
                              epa_nei_mod.o           future_emissions_mod.o \
                              gfed2_biomass_mod.o     global_hno3_mod.o      \
                              global_no3_mod.o        global_oh_mod.o        \
                              icoads_ship_mod.o       isoropia_mod.o  	     \
                              logical_mod.o           nei2005_anthro_mod.o   \
                              pbl_mix_mod.o           scale_anthro_mod.o     \
                              seasalt_mod.o           streets_anthro_mod.o   \
                              tracerid_mod.o          tracer_mod.o    	     \
                              tropopause_mod.o        uvalbedo_mod.o  	     \
                              vdiff_pre_mod.o         wetscav_mod.o

tagged_co_mod.o             : biofuel_mod.o           biomass_mod.o          \
                              dao_mod.o               diag_mod.o  	     \
                              diag_pl_mod.o           global_nox_mod.o       \
                              global_oh_mod.o         logical_mod.o  	     \
                              megan_mod.o             pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o  

tagged_ox_mod.o             : dao_mod.o               diag_mod.o             \
                              diag_pl_mod.o           drydep_mod.o           \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o  

tpcore_bc_mod.o             : logical_mod.o           tracer_mod.o  

tracerid_mod.o              : logical_mod.o           tracer_mod.o

transport_mod.o             : dao_mod.o               diag_mod.o             \
                              logical_mod.o           pjc_pfix_mod.o         \
                              tpcore_mod.o            tpcore_bc_mod.o        \
                              tpcore_fvdas_mod.o      tpcore_window_mod.o    \
                              tracer_mod.o                                   \
                              tpcore_geos5_window_mod.o                      \
                              pjc_pfix_geos5_window_mod.o


tropopause.o                : dao_mod.o               diag_mod.o             \
                              logical_mod.o           tropopause_mod.o 

tropopause_mod.o            : comode_mod.o            dao_mod.o              \
                              logical_mod.o

upbdflx_mod.o               : dao_mod.o               linoz_mod.o            \
                              logical_mod.o           tagged_ox_mod.o        \
                              tracerid_mod.o          tracer_mod.o           \
                              tropopause_mod.o

vdiff_pre_mod.o             : tracer_mod.o

vdiff_mod.o                 : comode_mod.o            dao_mod.o              \
                              diag_mod.o              drydep_mod.o           \
                              logical_mod.o           pbl_mix_mod.o          \
                              tracer_mod.o            vdiff_pre_mod.o

vistas_anthro_mod.o         : future_emissions_mod.o  epa_nei_mod.o          \
                              logical_mod.o           scale_anthro_mod.o     \
                              tracerid_mod.o          tracer_mod.o  

wetscav_mod.o               : dao_mod.o               diag_mod.o             \
                              logical_mod.o           ocean_mercury_mod.o    \
                              tracerid_mod.o          tracer_mod.o

xltmmp.o                    : dao_mod.o

xtra_read_mod.o             : dao_mod.o               diag_mod.o             \
                              logical_mod.o

#------------------------------------------------------------------------------
# NOTE: These dependencies are special, since they need to have the $(R8) 
# option applied to force all REAL's to be REAL*8 by default.  Therefore they
# require a customized compilation rule.  For these, we also need to list
# the name of the Fortran source code file so that "make" will interpret
# them properly. (bmy, 11/20/09)
#------------------------------------------------------------------------------

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

#EOC
