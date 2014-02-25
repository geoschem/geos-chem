#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
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
#  16 Dec 2010 - R. Yantosca - Renamed output files to "GC_Ref_Vol_4.*"
#  19 Jul 2011 - R. Yantosca - Changed *.f* to *.F* for ESMF compatibility
#  15 Jan 2014 - R. Yantosca - Now only save *.pdf output 
#EOP
#------------------------------------------------------------------------------
#BOC

# List of source code files (order is important)
SRC4 :=                                    \
./intro.gtmm                               \
$(GTMM)/GTMM.F90                           \
$(GTMM)/CasaRegridModule.F90               \
$(GTMM)/defineArrays.F90                   \
$(GTMM)/defineConstants.F90                \
$(GTMM)/dorestart_mod.F90                  \
$(GTMM)/input_gtmm_mod.F90                 \
$(GTMM)/loadCASAinput.F90                  \
./subs.gtmm                                \
$(GTMM)/CleanupCASAarrays.F90              \
$(GTMM)/GTMM_coupled.F90                   \
$(GTMM)/HgOutForGEOS.F90                   \
$(GTMM)/assignAgeClassToRunningPool.F90    \
$(GTMM)/assignRanPoolToAgeClass.F90        \
$(GTMM)/doFPARandLAI.F90                   \
$(GTMM)/doHerbCarbon.F90                   \
$(GTMM)/doHerbCarbonHg.F90                 \
$(GTMM)/doHerbivory.F90                    \
$(GTMM)/doHgDeposition.F90                 \
$(GTMM)/doLatitude.F90                     \
$(GTMM)/doLeafRootShedding.F90             \
$(GTMM)/doMaxHg.F90                        \
$(GTMM)/doNPP.F90                          \
$(GTMM)/doOptimumTemperature.F90           \
$(GTMM)/doPET.F90                          \
$(GTMM)/doSoilMoisture.F90                 \
$(GTMM)/doTreeCarbon.F90                   \
$(GTMM)/doTreeCarbonHg.F90                 \
$(GTMM)/getAgeClassBF.F90                  \
$(GTMM)/getFireParams.F90                  \
$(GTMM)/getFuelWood.F90                    \
$(GTMM)/getSoilMoistParams.F90             \
$(GTMM)/getSoilParams.F90                  \
$(GTMM)/loadHgDeposition.F90               \
$(GTMM)/load_GC_data.F90                   \
$(GTMM)/organizeAgeClasses.F90             \
$(GTMM)/processData.F90                    \
$(GTMM)/sort_pick_veg.F90


# Output file names
TEX4 := GC_Ref_Vol_4.tex
DVI4 := GC_Ref_Vol_4.dvi
PDF4 := GC_Ref_Vol_4.pdf


# Make commands
gtmmdoc: 
	rm -f $(TEX4)
	protex -sf $(SRC4) > $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	latex $(TEX4)
	dvipdf $(DVI4) $(PDF4)
	rm -f *.aux *.dvi *.log *.toc

#EOC
