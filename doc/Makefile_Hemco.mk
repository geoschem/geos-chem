#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: Makefile_Hemco.mk (in doc subdirectory)
#
# !DESCRIPTION: This Makefile fragment contains commands to build the 
#  documentation for the HEMCO Source Code.  It is inlined into
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
#  08 Jul 2014 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC

# List of source code files (order is important)
SRC5 :=                               \
./intro.hemco                         \
$(CORE)/hco_arr_mod.F90               \
$(CORE)/hco_calc_mod.F90              \
$(CORE)/hco_chartools_mod.F90	      \
$(CORE)/hco_clock_mod.F90	      \
$(CORE)/hco_config_mod.F90	      \
$(CORE)/hco_datacont_mod.F90	      \
$(CORE)/hco_diagn_mod.F90	      \
$(CORE)/hco_driver_mod.F90	      \
$(CORE)/hco_emislist_mod.F90	      \
$(CORE)/hco_error_mod.F90	      \
$(CORE)/hco_filedata_mod.F90	      \
$(CORE)/hco_fluxarr_mod.F90	      \
$(CORE)/hco_geotools_mod.F90	      \
$(CORE)/hco_readlist_mod.F90	      \
$(CORE)/hco_state_mod.F90	      \
$(CORE)/hco_tidx_mod.F90	      \
$(CORE)/hco_tools_mod.F90	      \
$(CORE)/hco_unit_mod.F90	      \
$(CORE)/hcoi_dataread_mod.F90	      \
$(CORE)/hcoi_gc_diagn_mod.F90	      \
$(CORE)/hcoi_gc_main_mod.F90	      \
$(CORE)/hcox_custom_mod.F90	      \
$(CORE)/hcox_driver_mod.F90	      \
$(CORE)/hcox_dustdead_mod.F	      \
$(CORE)/hcox_dustginoux_mod.F90	      \
$(CORE)/hcox_extlist_mod.F90	      \
$(CORE)/hcox_gc_RnPbBe_mod.F90	      \
$(CORE)/hcox_gfed3_mod.F90	      \
$(CORE)/hcox_lightnox_mod.F90	      \
$(CORE)/hcox_megan_mod.F	      \
$(CORE)/hcox_paranox_mod.F90	      \
$(CORE)/hcox_seaflux_mod.F90	      \
$(CORE)/hcox_seasalt_mod.F90	      \
$(CORE)/hcox_soilnox_mod.F90          \
$(CORE)/hcox_state_mod.F90



# Output file names
TEX5 := GC_Ref_Vol_5.tex
DVI5 := GC_Ref_Vol_5.dvi
PDF5 := GC_Ref_Vol_5.pdf

# Make commands
hemcodoc: 
	rm -f $(TEX5)
	protex -sf $(SRC5) > $(TEX5)
	latex $(TEX5)
	latex $(TEX5)
	latex $(TEX5)
	dvipdf $(DVI5) $(PDF5)
	rm -f *.aux *.dvi *.log *.toc

#EOC
