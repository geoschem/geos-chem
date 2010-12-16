#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: Makefile_GtmmDoc.mk (in doc subdirectory)
#
# !DESCRIPTION: This Makefile fragment contains commands to build the 
#  documentation for the Global Terrestrial Mercury Model (GTMM).  It is 
#  inlined into the Makefile (in the doc subdirectory) by an "include" 
#  command.
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
#  14 Sep 2010 - R. Yantosca - Initial version, split off from Makefile
#EOP
#------------------------------------------------------------------------------
#BOC

# List of source code files (order is important)
SRC4 :=                                    \
./intro.gtmm                               \
$(GTMM)/GTMM.f90                           \
$(GTMM)/CasaRegridModule.f90               \
$(GTMM)/defineArrays.f90                   \
$(GTMM)/defineConstants.f90                \
$(GTMM)/dorestart_mod.f90                  \
$(GTMM)/input_gtmm_mod.f90                 \
$(GTMM)/loadCASAinput.f90                  \
./subs.gtmm                                \
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


# Output file names
TEX4 := GEOS-Chem-GTMM.tex
DVI4 := GEOS-Chem-GTMM.dvi
PDF4 := GEOS-Chem-GTMM.pdf
PS4  := GEOS-Chem-GTMM.ps


# Make commands
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
