#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: Makefile_SrcDoc.mk (in doc subdirectory)
#
# !DESCRIPTION: This Makefile fragment contains commands to build the
#  documentation for the GEOS-Chem Source Code.  It is inlined into
#  the Makefile (in the doc subdirectory) by an "include" command.
#\\
#\\
# !REMARKS:
# To build the documentation, call "make" with the following syntax:
#                                                                             .
#   make TARGET [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# You must have the LaTeX utilities (latex, dvips, dvipdf) installed
# on your system in order to build the documentation.
#
# !REVISION HISTORY:
#  14 Sep 2010 - R. Yantosca - Initial version, split off from Makefile
#  See https://github.com/geoschem/geos-chem for complete history
#EOP
#------------------------------------------------------------------------------
#BOC

# List of source code files (order is important)
SRC2 :=                               \
./gc_intro.P                          \
./gc_headers.P                        \
$(wildcard $(HDR)/*.F*)               \
./gc_operations.P                     \
$(CORE)/main.F                        \
$(CORE)/input_mod.F                   \
$(CORE)/gc_environment_mod.F90        \
$(CORE)/cleanup.F                     \
./gc_transport.P                      \
$(CORE)/transport_mod.F               \
$(CORE)/pjc_pfix_mod.F                \
$(CORE)/tpcore_fvdas_mod.F90          \
$(CORE)/tpcore_window_mod.F90         \
./gc_convection.P                     \
$(CORE)/convection_mod.F              \
$(CORE)/wetscav_mod.F                 \
./gc_mixing.P                         \
$(CORE)/mixing_mod.F90                \
$(CORE)/pbl_mix_mod.F                 \
$(CORE)/vdiff_mod.F90                 \
$(CORE)/vdiff_pre_mod.F90             \
./gc_drydep.P                         \
$(CORE)/drydep_mod.F                  \
$(CORE)/modis_lai_mod.F90             \
$(CORE)/olson_landmap_mod.F90         \
$(CORE)/mapping_mod.F90               \
$(CORE)/get_ndep_mod.F                \
./gc_chemistry.P                      \
$(CORE)/emissions_mod.F90             \
$(CORE)/chemistry_mod.F90             \
$(CORE)/fast_jx_mod.F                 \
$(CORE)/uvalbedo_mod.F90              \
$(CORE)/set_global_ch4_mod.F90        \
$(CORE)/set_prof_o3.F                 \
$(CORE)/toms_mod.F                    \
$(CORE)/flexchem_mod.F90              \
$(KPP)/Standard/gckpp_HetRates.F90    \
$(CORE)/ucx_mod.F                     \
$(CORE)/cldice_HBrHOBr_rxn.F          \
$(CORE)/strat_chem_mod.F90            \
$(CORE)/linoz_mod.F                   \
./gc_aerosol.P                        \
$(CORE)/aerosol_mod.F                 \
$(CORE)/carbon_mod.F                  \
$(CORE)/dust_mod.F                    \
$(CORE)/seasalt_mod.F                 \
$(CORE)/sulfate_mod.F                 \
$(CORE)/isorropiaII_mod.F             \
./gc_met.P                            \
$(CORE)/dao_mod.F                     \
$(CORE)/flexgrid_read_mod.F90         \
./gc_specialty.P                      \
$(CORE)/co2_mod.F                     \
$(CORE)/exchange_mod.F                \
$(CORE)/global_ch4_mod.F              \
$(CORE)/mercury_mod.F                 \
$(CORE)/depo_mercury_mod.F            \
$(CORE)/land_mercury_mod.F            \
$(CORE)/pops_mod.F                    \
$(CORE)/RnPbBe_mod.F                  \
$(CORE)/rrtmg_rad_transfer_mod.F      \
$(CORE)/tagged_co_mod.F               \
$(CORE)/tagged_o3_mod.F

# Output file names
TEX2 := GC_v11-02_Core_Modules.tex
DVI2 := GC_v11-02_Core_Modules.dvi
PDF2 := GC_v11-02_Core_Modules.pdf

# Make commands
srcdoc:
	rm -f $(TEX2)
	./protex -sfp $(SRC2) > $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	dvipdf $(DVI2) $(PDF2)
	rm -f *.aux *.dvi *.log *.toc

#EOC
