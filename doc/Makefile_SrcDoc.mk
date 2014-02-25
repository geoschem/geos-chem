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
#  07 Feb 2012 - M. Payer    - Added paranox_mod, diag63_mod
#  08 Feb 2012 - R. Yantosca - Added geos57_read_mod.F90
#  28 Feb 2012 - R. Yantosca - Added pbl_mix_mod
#  05 Mar 2012 - M. Payer    - Added tracer_mod
#  06 Mar 2012 - R. Yantosca - Added photoj.F and set_prof.F
#  07 Mar 2012 - M. Payer    - Added global_ch4_mod
#  19 Mar 2012 - M. Payer    - Added EF_MGN20_mod
#  22 Mar 2012 - M. Payer    - Added c2h6_mod, olson_landmap_mod
#  29 Mar 2012 - R. Yantosca - Added lai_mod
#  29 Mar 2012 - R. Yantosca - Added modis_lai_mod and mapping_mod
#  09 Apr 2012 - R. Yantosca - Added modules from Headers/ directory
#  13 Apr 2012 - R. Yantosca - Removed findmon.F, rdlai.F, lai_mod.F
#  19 Apr 2012 - R. Yantosca - Added read_jv_atms_dat.F90
#  15 May 2012 - R. Yantosca - Added tpcore_bc_mod.F
#  22 May 2012 - M. Payer    - Add bromocarb_mod.F, cldice_HBrHOBr_rxn.F,
#                              and ssa_bromine_mod.F
#  31 Jul 2012 - R. Yantosca - Added FAST-J routines etc.
#  03 Aug 2012 - R. Yantosca - Added benchmark_mod, etc.
#  06 Aug 2012 - R. Yantosca - Added gcap_read_mod, etc.
#  14 Aug 2012 - R. Yantosca - Added gc_environment_mod, etc.
#  23 Oct 2012 - R. Yantosca - Added modules in ESMF
#  23 Oct 2012 - R. Yantosca - Added tagged_co_mod
#  23 Oct 2012 - M. Payer    - Added soil NOx modules; Removed upbdflx_mod.F
#  27 Nov 2012 - M. Payer    - Added modules for POPs simulation
#  13 Dec 2012 - R. Yantosca - Added biofit, sunparam, and removed some 
#                              obsolete functions
#  22 Jul 2013 - M. Sulprizio- Added rcp_mod
#  01 Aug 2013 - M. Sulprizio- Added aeic_mod
#  20 Aug 2013 - M. Sulprizio- Added carbon_mod
#  20 Aug 2013 - R. Yantosca - Remove reference to "define.h"
#  05 Sep 2013 - M. Sulprizio- Added global_hno3_mod
#  15 Jan 2014 - R. Yantosca - Now only create *.pdf file output
#  25 Feb 2014 - M. Sulprizio- Added a3_read_mod, a6_read_mod, and i6_read_mod
#EOP
#------------------------------------------------------------------------------
#BOC

# List of source code files (order is important)
SRC1 :=                               \
./intro.geos-chem                     \
./headers.geos-chem                   \
$(HDR)/CMN_SIZE_mod.F                 \
$(HDR)/CMN_DIAG_mod.F                 \
$(HDR)/CMN_GCTM_mod.F                 \
$(HDR)/CMN_NOX_mod.F                  \
$(HDR)/CMN_O3_mod.F                   \
$(HDR)/CMN_mod.F                      \
$(HDR)/cmn_fj_mod.F                   \
$(HDR)/commsoil_mod.F                 \
$(HDR)/comode_loop_mod.F              \
$(HDR)/EF_MGN20_mod.F                 \
$(HDR)/gigc_errcode_mod.F90           \
$(HDR)/gigc_state_chm_mod.F90         \
$(HDR)/gigc_state_met_mod.F90         \
$(HDR)/jv_cmn_mod.F                   \
$(HDR)/jv_mie_mod.F                   \
$(HDR)/smv_dimension_mod.F            \
$(HDR)/smv_physconst_mod.F            \
$(CORE)/main.F                        \
$(CORE)/a3_read_mod.F                 \
$(CORE)/a6_read_mod.F                 \
$(CORE)/acetone_mod.F                 \
$(CORE)/aeic_mod.F                    \
$(CORE)/aerosol_mod.F                 \
$(CORE)/arctas_ship_emiss_mod.F	      \
$(CORE)/benchmark_mod.F               \
$(CORE)/bravo_mod.F                   \
$(CORE)/bromocarb_mod.F               \
$(CORE)/c2h6_mod.F                    \
$(CORE)/cac_anthro_mod.F              \
$(CORE)/canopy_nox_mod.F              \
$(CORE)/carbon_mod.F                  \
$(CORE)/chemistry_mod.F               \
$(CORE)/co2_mod.F                     \
$(CORE)/comode_mod.F                  \
$(CORE)/convection_mod.F              \
$(CORE)/dao_mod.F                     \
$(CORE)/depo_mercury_mod.F            \
$(CORE)/diag03_mod.F                  \
$(CORE)/diag04_mod.F                  \
$(CORE)/diag41_mod.F                  \
$(CORE)/diag42_mod.F                  \
$(CORE)/diag49_mod.F                  \
$(CORE)/diag50_mod.F                  \
$(CORE)/diag51b_mod.F                 \
$(CORE)/diag53_mod.F                  \
$(CORE)/diag56_mod.F                  \
$(CORE)/diag63_mod.F                  \
$(CORE)/diag_pl_mod.F                 \
$(CORE)/diag_oh_mod.F                 \
$(CORE)/diag_mod.F                    \
$(CORE)/drydep_mod.F                  \
$(CORE)/dust_mod.F                    \
$(CORE)/emep_mod.F                    \
$(CORE)/emissions_mod.F               \
$(CORE)/fjx_acet_mod.F                \
$(CORE)/gamap_mod.F                   \
$(CORE)/gcap_read_mod.F               \
$(CORE)/get_ndep_mod.F                \
$(CORE)/gigc_environment_mod.F90      \
$(CORE)/gc_type_mod.F                 \
$(CORE)/geosfp_read_mod.F90           \
$(CORE)/get_popsinfo_mod.F            \
$(CORE)/gfed3_biomass_mod.F           \
$(CORE)/global_bc_mod.F               \
$(CORE)/global_br_mod.F               \
$(CORE)/global_ch4_mod.F              \
$(CORE)/global_hno3_mod.F             \
$(CORE)/global_no3_mod.F              \
$(CORE)/global_nox_mod.F              \
$(CORE)/global_o1d_mod.F              \
$(CORE)/global_o3_mod.F               \
$(CORE)/global_oc_mod.F               \
$(CORE)/global_oh_mod.F               \
$(CORE)/h2_hd_mod.F                   \
$(CORE)/i6_read_mod.F                 \
$(CORE)/icoads_ship_mod.F             \
$(CORE)/input_mod.F                   \
$(CORE)/isoropiaII_mod.F              \
$(CORE)/land_mercury_mod.F            \
$(CORE)/lightning_nox_mod.F           \
$(CORE)/linoz_mod.F                   \
$(CORE)/logical_mod.F                 \
$(CORE)/mapping_mod.F90               \
$(CORE)/megan_mod.F                   \
$(CORE)/meganut_mod.F                 \
$(CORE)/merra_a1_mod.F                \
$(CORE)/merra_a3_mod.F                \
$(CORE)/merra_cn_mod.F                \
$(CORE)/merra_i6_mod.F                \
$(CORE)/modis_lai_mod.F90             \
$(CORE)/nei2005_anthro_mod.F          \
$(CORE)/olson_landmap_mod.F90         \
$(CORE)/optdepth_mod.F                \
$(CORE)/paranox_mod.F                 \
$(CORE)/pbl_mix_mod.F                 \
$(CORE)/pjc_pfix_mod.F                \
$(CORE)/planeflight_mod.F             \
$(CORE)/pops_mod.F                    \
$(CORE)/rcp_mod.F                     \
$(CORE)/retro_mod.F                   \
$(CORE)/RnPbBe_mod.F                  \
$(CORE)/scale_anthro_mod.F            \
$(CORE)/seasalt_mod.F                 \
$(CORE)/soil_nox_mod.F                \
$(CORE)/soilnox_restart_mod.F         \
$(CORE)/ssa_bromine_mod.F             \
$(CORE)/strat_chem_mod.F90            \
$(CORE)/sulfate_mod.F                 \
$(CORE)/tagged_co_mod.F               \
$(CORE)/tagged_ox_mod.F               \
$(CORE)/toms_mod.F                    \
$(CORE)/tpcore_bc_mod.F               \
$(CORE)/tracer_mod.F                  \
$(CORE)/tropopause_mod.F              \
$(CORE)/tpcore_fvdas_mod.F90          \
$(CORE)/tpcore_geos5_window_mod.F90   \
$(CORE)/transport_mod.F               \
$(CORE)/vdiff_mod.F90                 \
$(CORE)/vdiff_pre_mod.F               \
$(CORE)/vistas_anthro_mod.F           \
./subs.geos-chem                      \
$(CORE)/anthroems.F                   \
$(CORE)/biofit.F                      \
$(CORE)/boxvl.F                       \
$(CORE)/cldice_HBrHOBr_rxn.F          \
$(CORE)/diag1.F                       \
$(CORE)/diag3.F                       \
$(CORE)/diag_2pm.F                    \
$(CORE)/diagoh.F                      \
$(CORE)/emfossil.F                    \
$(CORE)/emf_scale.F                   \
$(CORE)/fast_j.F                      \
$(CORE)/gasconc.F                     \
$(CORE)/JRATET.F                      \
$(CORE)/JVALUE.F                      \
$(CORE)/jv_index.F                    \
$(CORE)/initialize.F                  \
$(CORE)/inphot.F                      \
$(CORE)/lump.F                        \
$(CORE)/ndxx_setup.F                  \
$(CORE)/ohsave.F                      \
$(CORE)/OPMIE.F                       \
$(CORE)/partition.F                   \
$(CORE)/photoj.F                      \
$(CORE)/physproc.F                    \
$(CORE)/RD_AOD.F                      \
$(CORE)/rd_js.F                       \
$(CORE)/RD_TJPL.F                     \
$(CORE)/read_jv_atms_dat.F90          \
$(CORE)/set_aer.F                     \
$(CORE)/setemdep.F                    \
$(CORE)/set_prof.F                    \
$(CORE)/SPHERE.F                      \
$(CORE)/read_jv_atms_dat.F90          \
$(CORE)/ruralbox.F                    \
$(CORE)/setemis.F                     \
$(CORE)/sfcwindsqr.F                  \
$(CORE)/sunparam.F


# Output file names
TEX1 := GC_Ref_Vol_3.tex
DVI1 := GC_Ref_Vol_3.dvi
PDF1 := GC_Ref_Vol_3.pdf


# Make commands
srcdoc: 
	rm -f $(TEX1)
	protex -sf $(SRC1) > $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	latex $(TEX1)
	dvipdf $(DVI1) $(PDF1)
	rm -f *.aux *.dvi *.log *.toc

#EOC
