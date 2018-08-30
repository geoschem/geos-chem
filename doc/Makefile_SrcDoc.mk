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
#  25 Feb 2014 - M. Sulprizio- Added chemgrid_mod, fast_jx_mod, and ucx_mod.
#                              Removed references to tropopause_mod and routines
#                              for Fast-J.
#  08 Jul 2014 - R. Yantosca - Removed obsolete routines from the list
#  06 Jan 2015 - M. Sulprizio- Remove additional routines made obsolete by HEMCO
#  07 Jan 2015 - R. Yantosca - Added exchange_mod (i.e. 2-way nesting)
#  15 Jan 2015 - M. Sulprizio- Added rrtmg_rad_transfer_mod.F and set_prof_o3.F
#  04 Mar 2015 - R. Yantosca - Add uvalbedo_mod.F
#  10 Jul 2015 - R. Yantosca - Use ./protex to avoid problems on some systems
#  10 Jul 2015 - R. Yantosca - Updated list of files as of v11-01b
#  29 Aug 2016 - M. Sulprizio- Updated list of files as of v11-01g
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
$(CORE)/restart_mod.F                 \
$(CORE)/cleanup.F                     \
./gc_transport.P                      \
$(CORE)/transport_mod.F               \
$(CORE)/pjc_pfix_mod.F                \
$(CORE)/tpcore_bc_mod.F               \
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
$(CORE)/bromocarb_mod.F               \
$(CORE)/cldice_HBrHOBr_rxn.F          \
$(CORE)/strat_chem_mod.F90            \
$(CORE)/linoz_mod.F                   \
./gc_aerosol.P                        \
$(CORE)/aerosol_mod.F                 \
$(CORE)/carbon_mod.F                  \
$(CORE)/dust_mod.F                    \
$(CORE)/seasalt_mod.F                 \
$(CORE)/sulfate_mod.F                 \
$(CORE)/isoropiaII_mod.F              \
./gc_met.P                            \
$(CORE)/dao_mod.F                     \
$(CORE)/merra2_read_mod.F90           \
$(CORE)/geosfp_read_mod.F90           \
./gc_specialty.P                      \
$(CORE)/c2h6_mod.F                    \
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
$(CORE)/tagged_o3_mod.F               \


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
