#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in doc subdirectory)
#
# !DESCRIPTION: Makefile for building the documentation (in PDF and 
#  PostScript formats) for the GEOS-Chem model source code, makefiles,
#  and related scripts.
#\\
#\\
# !REMARKS:
# To build the documentation, call "make" with the following syntax:
#
#   make TARGET [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#
# You must have the LaTeX utilities (latex, dvips, dvipdf) installed
# on your system in order to build the documentation.
#
# !REVISION HISTORY: 
#  16 Sep 2009 - R. Yantosca - Initial version
#  21 Sep 2009 - R. Yantosca - Now call Makefile in help directory to
#                              display the help screen options
#  19 Nov 2009 - R. Yantosca - Now build documentation for Makefiles in the
#                              GeosUtil and GeosCore subdirectories
#  11 Dec 2009 - R. Yantosca - Now get SHELL from Makefile_header.mk
#  28 Jan 2010 - R. Yantosca - Now add the Makefiles in the GeosTomas and
#                              ISOROPIA subdirs to the documentation
#  27 Aug 2010 - R. Yantosca - Brought up to date for newly documented code
#  31 Aug 2010 - R. Yantosca - Add documentation for GTMM routines
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================

# Define variables
ROOTDIR := ..
BIN     := $(ROOTDIR)/bin
BPCH    := $(ROOTDIR)/GeosBpch
CORE    := $(ROOTDIR)/GeosCore
GTMM    := $(ROOTDIR)/GTMM
HDR     := $(ROOTDIR)/Headers
HELP    := $(ROOTDIR)/help
ISO     := $(ROOTDIR)/ISOROPIA
KPP     := $(ROOTDIR)/KPP
TOM     := $(ROOTDIR)/GeosTomas
UTIL    := $(ROOTDIR)/GeosUtil

# Get the Unix shell in SHELL from the Makefile_header.mk
include $(ROOTDIR)/Makefile_header.mk

#==============================================================================
# Makefile targets
#==============================================================================

.PHONY: all docclean help

all: srcdoc makedoc utildoc gtmmdoc

clean:
	rm -f *.tex *.ps *.pdf

help:
	@$(MAKE) -C $(HELPDIR)

#==============================================================================
# Build documentation for GEOS-Chem source code
#==============================================================================

SRC1 :=                              \
./geos_chem.intro                    \
./geos_chem.headers                  \
$(HDR)/define.h                      \
$(HDR)/CMN_SIZE                      \
$(HDR)/CMN_DIAG                      \
$(CORE)/arctas_ship_emiss_mod.f	     \
$(CORE)/bravo_mod.f		     \
$(CORE)/cac_anthro_mod.f	     \
$(CORE)/chemistry_mod.f		     \
$(CORE)/co2_mod.f		     \
$(CORE)/convection_mod.f	     \
$(CORE)/dao_mod.f		     \
$(CORE)/depo_mercury_mod.f	     \
$(CORE)/diag04_mod.f		     \
$(CORE)/diag_mod.f		     \
$(CORE)/dust_mod.f		     \
$(CORE)/emep_mod.f		     \
$(CORE)/emissions_mod.f		     \
$(CORE)/fjx_acet_mod.f		     \
$(CORE)/gamap_mod.f		     \
$(CORE)/global_o3_mod.f		     \
$(CORE)/icoads_ship_mod.f	     \
$(CORE)/input_mod.f		     \
$(CORE)/isoropiaII_mod.f	     \
$(CORE)/land_mercury_mod.f	     \
$(CORE)/linoz_mod.f		     \
$(CORE)/logical_mod.f		     \
$(CORE)/megan_mod.f	             \
$(CORE)/meganut_mod.f		     \
$(CORE)/merra_a1_mod.f		     \
$(CORE)/merra_a3_mod.f		     \
$(CORE)/merra_cn_mod.f		     \
$(CORE)/merra_i6_mod.f		     \
$(CORE)/nei2005_anthro_mod.f	     \
$(CORE)/pjc_pfix_mod.f		     \
$(CORE)/scale_anthro_mod.f	     \
$(CORE)/tpcore_fvdas_mod.f90	     \
$(CORE)/tpcore_geos5_window_mod.f90  \
$(CORE)/transport_mod.f		     \
$(CORE)/vdiff_mod.f90		     \
$(CORE)/vdiff_pre_mod.f		     \
$(CORE)/vistas_anthro_mod.f	     \
./geos_chem.subs                     \
$(CORE)/diag1.f			     \
$(CORE)/diag3.f			     \
$(CORE)/fast_j.f		     \
$(CORE)/initialize.f		     \
$(CORE)/ndxx_setup.f		     \
$(CORE)/rdlai.f			     \
$(CORE)/readlai.f                        

TEX1 := GEOS-Chem.tex
DVI1 := GEOS-Chem.dvi
PDF1 := GEOS-Chem.pdf
PS1  := GEOS-Chem.ps

srcdoc: 
	rm -f $(TEX1)
	protex -sf $(SRC1) > $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	dvipdf $(DVI1) $(PDF1)
	dvips $(DVI1) -o $(PS1)
	rm -f *.aux *.dvi *.log *.toc

#==============================================================================
# Build documentation for GEOS-Chem Makefiles
#==============================================================================

SRC2 :=                        \
./make.intro                   \
$(ROOTDIR)/Makefile            \
$(ROOTDIR)/Makefile_header.mk  \
$(UTIL)/Makefile               \
$(ISO)/Makefile                \
$(CODE)/Makefile               \
$(KPP)/Makefile                \
$(KPP)/standard/Makefile       \
$(KPP)/SOA/Makefile            \
$(TOM)/Makefile                \
$(GTMM)/Makefile               \
./Makefile                     \
$(HELP)/Makefile

TEX2 := GEOS-Chem-Makefiles.tex
DVI2 := GEOS-Chem-Makefiles.dvi
PDF2 := GEOS-Chem-Makefiles.pdf
PS2  := GEOS-Chem-Makefiles.ps

makedoc: 
	rm -f $(TEX2)
	protex -fS $(SRC2) > $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	dvipdf $(DVI2) $(PDF2)
	dvips $(DVI2) -o $(PS2)
	rm -f *.aux *.dvi *.log *.toc

#==============================================================================
# Build documentation for GEOS-Chem utility modules
#==============================================================================

SRC3 := ./util.intro $(wildcard $(UTIL)/*.f)
TEX3 := GEOS-Chem-Utility.tex
DVI3 := GEOS-Chem-Utility.dvi
PDF3 := GEOS-Chem-Utility.pdf
PS3  := GEOS-Chem-Utility.ps

utildoc: 
	rm -f $(TEX3)
	protex -sf $(SRC3) > $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	latex $(TEX3)
	dvipdf $(DVI3) $(PDF3)
	dvips $(DVI3) -o $(PS3)
	rm -f *.aux *.dvi *.log *.toc

#==============================================================================
# Build documentation for GTMM documentation utility modules
#==============================================================================

SRC4 :=                                    \
./gtmm.intro                               \
$(GTMM)/GTMM.f90                           \
$(GTMM)/CasaRegridModule.f90               \
$(GTMM)/defineArrays.f90                   \
$(GTMM)/defineConstants.f90                \
$(GTMM)/dorestart_mod.f90                  \
$(GTMM)/input_gtmm_mod.f90                 \
$(GTMM)/loadCASAinput.f90                  \
./gtmm.subs                                \
$(GTMM)/CleanupCASAarrays.f90              \
$(GTMM)/GTMM_coupled.f90                   \
$(GTMM)/HgOutForGEOS.f90                   \
$(GTMM)/assignAgeClassToRunningPool.f90    \
$(GTMM)/assignRanPoolToAgeClass.f90        \
$(GTMM)/doFPARandLAI.f90                   \
$(GTMM)/doHerbCarbon.f90                   \
$(GTMM)/doHerbCarbonHg.f90                 \
$(GTMM)/doHerbivory.f90                    \
$(GTMM)/doHgDeposition.f90                 \
$(GTMM)/doLatitude.f90                     \
$(GTMM)/doLeafRootShedding.f90             \
$(GTMM)/doMaxHg.f90                        \
$(GTMM)/doNPP.f90                          \
$(GTMM)/doOptimumTemperature.f90           \
$(GTMM)/doPET.f90                          \
$(GTMM)/doSoilMoisture.f90                 \
$(GTMM)/doTreeCarbon.f90                   \
$(GTMM)/doTreeCarbonHg.f90                 \
$(GTMM)/getAgeClassBF.f90                  \
$(GTMM)/getFireParams.f90                  \
$(GTMM)/getFuelWood.f90                    \
$(GTMM)/getSoilMoistParams.f90             \
$(GTMM)/getSoilParams.f90                  \
$(GTMM)/loadHgDeposition.f90               \
$(GTMM)/load_GC_data.f90                   \
$(GTMM)/organizeAgeClasses.f90             \
$(GTMM)/processData.f90                    \
$(GTMM)/sort_pick_veg.f90

TEX4 := GEOS-Chem-GTMM.tex
DVI4 := GEOS-Chem-GTMM.dvi
PDF4 := GEOS-Chem-GTMM.pdf
PS4  := GEOS-Chem-GTMM.ps

gtmmdoc: 
	rm -f $(TEX4)
	protex -sf $(SRC4) > $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	dvipdf $(DVI4) $(PDF4)
	dvips $(DVI4) -o $(PS4)
	rm -f *.aux *.dvi *.log *.toc

#EOC
