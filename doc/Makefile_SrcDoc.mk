#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
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
#  14 Sep 2010 - R. Yantosca - Added optdepth_mod.f to list
#  15 Sep 2010 - R. Yantosca - Added diag_2pm, diag_56, diagoh, ohsave
#  16 Sep 2010 - R. Yantosca - Added diag_pl_mod
#  04 Nov 2010 - R. Yantosca - Added acetone_mod
#  10 Nov 2010 - R. Yantosca - Added lightning_nox_mod
#  19 Nov 2010 - R. Yantosca - Added anthroems, RnPbBe_mod, tagged_ox_mod
#  19 Nov 2010 - R. Yantosca - Added tcorr, emfossil, emf_scale
#EOP
#------------------------------------------------------------------------------
#BOC

# List of source code files (order is important)
SRC1 :=                              \
./intro.geos-chem                    \
./headers.geos-chem                  \
$(HDR)/define.h                      \
$(HDR)/CMN_SIZE                      \
$(HDR)/CMN_DIAG                      \
$(CORE)/acetone_mod.f                \
$(CORE)/arctas_ship_emiss_mod.f	     \
$(CORE)/bravo_mod.f                  \
$(CORE)/cac_anthro_mod.f             \
$(CORE)/chemistry_mod.f              \
$(CORE)/co2_mod.f                    \
$(CORE)/convection_mod.f             \
$(CORE)/dao_mod.f                    \
$(CORE)/depo_mercury_mod.f           \
$(CORE)/diag04_mod.f                 \
$(CORE)/diag56_mod.f                 \
$(CORE)/diag_pl_mod.f                \
$(CORE)/diag_oh_mod.f                \
$(CORE)/diag_mod.f                   \
$(CORE)/dust_mod.f                   \
$(CORE)/emep_mod.f                   \
$(CORE)/emissions_mod.f              \
$(CORE)/fjx_acet_mod.f               \
$(CORE)/gamap_mod.f                  \
$(CORE)/global_br_mod.f              \
$(CORE)/global_no3_mod.f             \
$(CORE)/global_nox_mod.f             \
$(CORE)/global_o1d_mod.f             \
$(CORE)/global_o3_mod.f              \
$(CORE)/global_oh_mod.f              \
$(CORE)/icoads_ship_mod.f            \
$(CORE)/input_mod.f                  \
$(CORE)/isoropiaII_mod.f             \
$(CORE)/land_mercury_mod.f           \
$(CORE)/lightning_nox_mod.f          \
$(CORE)/linoz_mod.f                  \
$(CORE)/logical_mod.f                \
$(CORE)/megan_mod.f                  \
$(CORE)/meganut_mod.f                \
$(CORE)/merra_a1_mod.f               \
$(CORE)/merra_a3_mod.f               \
$(CORE)/merra_cn_mod.f               \
$(CORE)/merra_i6_mod.f               \
$(CORE)/nei2005_anthro_mod.f         \
$(CORE)/optdepth_mod.f               \
$(CORE)/pjc_pfix_mod.f               \
$(CORE)/RnPbBe_mod.f                 \
$(CORE)/scale_anthro_mod.f           \
$(CORE)/tagged_ox_mod.f              \
$(CORE)/toms_mod.f                   \
$(CORE)/tropopause_mod.f             \
$(CORE)/tpcore_fvdas_mod.f90         \
$(CORE)/tpcore_geos5_window_mod.f90  \
$(CORE)/transport_mod.f              \
$(CORE)/vdiff_mod.f90                \
$(CORE)/vdiff_pre_mod.f              \
$(CORE)/vistas_anthro_mod.f          \
./subs.geos-chem                     \
$(CORE)/anthroems.f                  \
$(CORE)/diag1.f                      \
$(CORE)/diag3.f                      \
$(CORE)/diag_2pm.f                   \
$(CORE)/diagoh.f                     \
$(CORE)/emfossil.f                   \
$(CORE)/emf_scale.f                  \
$(CORE)/fast_j.f                     \
$(CORE)/initialize.f                 \
$(CORE)/ndxx_setup.f                 \
$(CORE)/ohsave.f                     \
$(CORE)/rdlai.f                      \
$(CORE)/readlai.f                    \
$(CORE)/ruralbox.f                   \
$(CORE)/setemis.f                    \
$(CORE)/tcorr.f            


# Output file names
TEX1 := GEOS-Chem.tex
DVI1 := GEOS-Chem.dvi
PDF1 := GEOS-Chem.pdf
PS1  := GEOS-Chem.ps


# Make commands
srcdoc: 
	rm -f $(TEX1)
	protex -sf $(SRC1) > $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	dvipdf $(DVI1) $(PDF1)
	dvips $(DVI1) -o $(PS1)
	rm -f *.aux *.dvi *.log *.toc

#EOC
