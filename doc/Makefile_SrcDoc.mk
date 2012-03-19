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
#  01 Dec 2010 - R. Yantosca - Added global_br_mod, global_no3_mod
#  01 Dec 2010 - R. Yantosca - Added global_nox_mod, global_o1d_mod
#  01 Dec 2010 - R. Yantosca - Added global_oh_mod,  toms_mod
#  02 Dec 2010 - R. Yantosca - Added upbdflx_mod, diag41_mod, diag42_mod
#  02 Dec 2010 - R. Yantosca - Added diag03_mod, diag49_mod, diag50_mod
#  02 Dec 2010 - R. Yantosca - Added diag51_mod, diag51b_mod, boxvl, rdmonot
#  02 Dec 2010 - R. Yantosca - Added rdlight, rdland, rdsoil, emmonot
#  16 Dec 2010 - R. Yantosca - Renamed output files to "GC_Ref_Vol_3.*"\
#  21 Dec 2010 - R. Yantosca - Added comode_mod
#  11 Jul 2011 - R. Yantosca - Added restart_mod
#  19 Jul 2011 - R. Yantosca - Changed *.f* to *.F* for ESMF compatibility
#  29 Jul 2011 - R. Yantosca - Added planeflight_mod
#  22 Aug 2011 - R. Yantosca - Added retro_mod
#  07 Sep 2011 - R. Yantosca - Added gfed3_biomass_mod, *jv*_mod files
#  22 Dec 2011 - M. Payer    - Added aerosol_mod, drydep_mod, seasalt_mod,
#                              and sulfate_mod
#  07 Feb 2012 - M. Payer    - Added paranox_mod, diag59_mod
#  08 Feb 2012 - R. Yantosca - Added geos57_read_mod.F90
#  19 Mar 2012 - M. Payer    - Added EF_MGN20_mod
#EOP
#------------------------------------------------------------------------------
#BOC

# List of source code files (order is important)
SRC1 :=                              \
./intro.geos-chem                    \
./headers.geos-chem                  \
$(HDR)/define.h                      \
$(HDR)/CMN_SIZE_mod.F                \
$(HDR)/CMN_DIAG_mod.F                \
$(HDR)/cmn_fj_mod.F                  \
$(HDR)/EF_MGN20_mod.F                \
$(HDR)/jv_cmn_mod.F                  \
$(HDR)/jv_mie_mod.F                  \
$(CORE)/main.F                       \
$(CORE)/acetone_mod.F                \
$(CORE)/aerosol_mod.F                \
$(CORE)/arctas_ship_emiss_mod.F	     \
$(CORE)/bravo_mod.F                  \
$(CORE)/cac_anthro_mod.F             \
$(CORE)/chemistry_mod.F              \
$(CORE)/co2_mod.F                    \
$(CORE)/comode_mod.F                 \
$(CORE)/convection_mod.F             \
$(CORE)/dao_mod.F                    \
$(CORE)/depo_mercury_mod.F           \
$(CORE)/diag03_mod.F                 \
$(CORE)/diag04_mod.F                 \
$(CORE)/diag41_mod.F                 \
$(CORE)/diag42_mod.F                 \
$(CORE)/diag49_mod.F                 \
$(CORE)/diag50_mod.F                 \
$(CORE)/diag51b_mod.F                \
$(CORE)/diag56_mod.F                 \
$(CORE)/diag59_mod.F                 \
$(CORE)/diag_pl_mod.F                \
$(CORE)/diag_oh_mod.F                \
$(CORE)/diag_mod.F                   \
$(CORE)/drydep_mod.F                 \
$(CORE)/dust_mod.F                   \
$(CORE)/emep_mod.F                   \
$(CORE)/emissions_mod.F              \
$(CORE)/fjx_acet_mod.F               \
$(CORE)/gamap_mod.F                  \
$(CORE)/geos57_read_mod.F            \
$(CORE)/gfed3_biomass_mod.F          \
$(CORE)/global_br_mod.F              \
$(CORE)/global_no3_mod.F             \
$(CORE)/global_nox_mod.F             \
$(CORE)/global_o1d_mod.F             \
$(CORE)/global_o3_mod.F              \
$(CORE)/global_oh_mod.F              \
$(CORE)/h2_hd_mod.F                  \
$(CORE)/icoads_ship_mod.F            \
$(CORE)/input_mod.F                  \
$(CORE)/isoropiaII_mod.F             \
$(CORE)/land_mercury_mod.F           \
$(CORE)/lightning_nox_mod.F          \
$(CORE)/linoz_mod.F                  \
$(CORE)/logical_mod.F                \
$(CORE)/megan_mod.F                  \
$(CORE)/meganut_mod.F                \
$(CORE)/merra_a1_mod.F               \
$(CORE)/merra_a3_mod.F               \
$(CORE)/merra_cn_mod.F               \
$(CORE)/merra_i6_mod.F               \
$(CORE)/nei2005_anthro_mod.F         \
$(CORE)/optdepth_mod.F               \
$(CORE)/paranox_mod.F                \
$(CORE)/pjc_pfix_mod.F               \
$(CORE)/planeflight_mod.F            \
$(CORE)/retro_mod.F                  \
$(CORE)/RnPbBe_mod.F                 \
$(CORE)/scale_anthro_mod.F           \
$(CORE)/seasalt_mod.F                \
$(CORE)/sulfate_mod.F                \
$(CORE)/tagged_ox_mod.F              \
$(CORE)/toms_mod.F                   \
$(CORE)/tropopause_mod.F             \
$(CORE)/tpcore_fvdas_mod.F90         \
$(CORE)/tpcore_geos5_window_mod.F90  \
$(CORE)/transport_mod.F              \
$(CORE)/upbdflx_mod.F                \
$(CORE)/vdiff_mod.F90                \
$(CORE)/vdiff_pre_mod.F              \
$(CORE)/vistas_anthro_mod.F          \
./subs.geos-chem                     \
$(CORE)/anthroems.F                  \
$(CORE)/boxvl.F                      \
$(CORE)/diag1.F                      \
$(CORE)/diag3.F                      \
$(CORE)/diag_2pm.F                   \
$(CORE)/diagoh.F                     \
$(CORE)/emfossil.F                   \
$(CORE)/emf_scale.F                  \
$(CORE)/emmonot.F                    \
$(CORE)/fast_j.F                     \
$(CORE)/findmon.F                    \
$(CORE)/initialize.F                 \
$(CORE)/ndxx_setup.F                 \
$(CORE)/ohsave.F                     \
$(CORE)/rdlai.F                      \
$(CORE)/rdland.F                     \
$(CORE)/rdsoil.F                     \
$(CORE)/rdlight.F                    \
$(CORE)/rdmonot.F                    \
$(CORE)/readlai.F                    \
$(CORE)/ruralbox.F                   \
$(CORE)/setemis.F                    \
$(CORE)/sfcwindsqr.F                 \
$(CORE)/tcorr.F


# Output file names
TEX1 := GC_Ref_Vol_3.tex
DVI1 := GC_Ref_Vol_3.dvi
PDF1 := GC_Ref_Vol_3.pdf
PS1  := GC_Ref_Vol_3.ps


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
