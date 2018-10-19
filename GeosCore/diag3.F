!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag3.F
!
! !DESCRIPTION: Subroutine DIAG3 prints out diagnostics to the BINARY PUNCH
!  format file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG3( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
      ! Modules from Headers directory
      USE CMN_DIAG_MOD                       ! Diagnostic switches & arrays
      USE CMN_FJX_MOD                        ! Fast-JX flux diagnostics
      USE CMN_O3_MOD                         ! FMOL
      USE CMN_SIZE_MOD                       ! Size parameters
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants,      ONLY : AVO     ! Avogadro's #
      USE PRECISION_MOD                      ! For GEOS-Chem Precision (fp)
      USE State_Met_Mod,      ONLY : MetState
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_

      ! Modules from GeosUtil directory
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE FILE_MOD                           
      USE GC_GRID_MOD,        ONLY : GET_AREA_M2
      USE GC_GRID_MOD,        ONLY : GET_YOFFSET, GET_XOFFSET
      USE TIME_MOD

      ! Modules from GeosCore directory
#if defined( BPCH_DIAG )
      USE BPCH2_MOD                          ! For binary punch I/O routines
      USE DIAG_MOD                           ! For diagnostic arrays
      USE DIAG03_MOD                         ! For Hg diagnostic
      USE DIAG04_MOD                         ! For CO2 diagnostics
      USE DIAG41_MOD                         ! For afternoon PBL diag
      USE DIAG42_MOD                         ! For SOA diag
      USE DIAG53_MOD                         ! For POPs diag
      USE DIAG56_MOD                         ! For time in tropopause diag
#endif
      USE DEPO_MERCURY_MOD                   ! For offline Hg simulation
      USE DRYDEP_MOD                         ! For dry deposition
      USE WETSCAV_MOD                        ! For wet deposition
      
#if   defined( TOMAS )
      USE TOMAS_MOD, ONLY : ICOMP, IDIAG, IBINS  !(win, 1/25/10)
#endif

      ! Interface w/ HEMCO
      USE HCO_TYPES_MOD, ONLY : DiagnCont
      USE HCO_DIAGN_MOD
      USE HCO_ERROR_MOD
      USE HCO_INTERFACE_MOD

      ! Species database
      USE Species_Mod,       ONLY : Species

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  ROUTINE diag3.F SAVES TIME-AVERAGED QUANTITIES TO THE BPCH FILE.   %%%
!  %%%  THE BPCH FILE FORMAT IS DEPRECATED.  DIAGNOSTIC OUTPUT WILL        %%%
!  %%%  EVENTUALLY BE REPLACED BY netCDF OUTPUT IN A FUTURE RELEASE.       %%%
!  %%%  FOR NOW, WE SHALL KEEP BPCH OUTPUT IN ORDER TO PRESERVE BACKWARDS  %%%
!  %%%  COMPATIBILITY WITH THE EXISTING VISUALIZATION SOFTWARE.            %%%
!  %%%                                                                     %%%
!  %%%  ALSO NOTE: MANY EMISSIONS DIAGNOSTICS ARE NOW ARCHIVED IN HEMCO,   %%%
!  %%%  AND ARE WRITTEN OUT TO BPCH FILE HERE.  THIS IS MEANT TO PRESERVE  %%%
!  %%%  THE EXISTING DIAGNOSTIC CAPABILITY FOR THE TIME BEING.             %%%
!  %%%                                                                     %%%
!  %%%     -- Bob Yantosca (14 Aug 2014)                                   %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY: 
!  (40) Bug fix: Save levels 1:LD13 for ND13 diagnostic for diagnostic
!        categories "SO2-AC-$" and "SO2-EV-$".  Now reference F90 module
!        "tracerid_mod.f".  Now reference NUMDEP from "drydep_mod.f".
!        Now save anthro, biofuel, biomass NH3 in ND13; also fixed ND13
!        tracer numbers.  For ND13, change scale factor from SCALESRCE to 1.
!        Now references "wetscav_mod.f".  Now also save true tracer numbers 
!        for ND38 and ND39 diagnostic.  Now also write out biomass SO2.
!        Now convert ND01, ND02, ND44 diagnostics for Rn/Pb/Be from kg to 
!        kg/s here. (bmy, 1/24/03)
!  (41) Now save out natural NH3 in ND13 as "NH3-NATU" (rjp, bmy, 3/23/03)
!  (42) Now replace DXYP(JREF) by routine GET_AREA_M2, GET_XOFFSET, and
!        GET_YOFFSET of "grid_mod.f".  Now references "time_mod.f".
!        DIAGb, DIAGe are now local variables.  Now remove obsolete statements
!        IF ( LBPNCH > 0 ).  Removed SCALE1, replaced with SCALEDYN. 
!        (bmy, 2/24/03)
!  (43) Added TSKIN, PARDF, PARDR, GWET to ND67 diagnostic.  For GEOS-4/fvDAS,
!        UWND, VWND, TMPU, SPHU are A-6 fields.  Adjust the ND66 scale factors 
!        accordingly.  Delete KZZ from ND66.  Updated comments. (bmy, 6/23/03)
!  (44) Bug fix: use LD68 instead of ND68 in DO-loop to avoid out-of-bounds 
!        error. (bec, bmy, 7/15/03)
!  (45) Now print out NTRACE drydep fluxes for tagged Ox.  Also tagged Ox 
!        now saves drydep in molec/cm2/s.  Now print out Kr85 prod/loss in 
!        ND03. (bmy, 8/20/03)
!  (46) Now use actual tracer number for ND37 diagnostic. (bmy, 1/21/04)
!  (47) Now loop over the actual # of soluble tracers for ND17, ND18.  
!        (bmy, 3/19/04)
!  (48) Now use the actual tracer # for ND17 and ND18 diagnostics. 
!        Rearrange ND44 code for clarity. (bmy, 3/23/04)
!  (49) Added ND06 (dust aerosol) and ND07 (carbon aerosol) diagnostics.
!        Now scale online dust optical depths by SCALECHEM in ND21 diagnostic.
!        (rjp, tdf, bmy, 4/5/04)
!  (50) Added ND08 (seasalt aerosol) diagnostic (rjp, bec, bmy, 4/20/04)
!  (51) Now save out SO2 from ships (if LSHIPSO2=T) (bec, bmy, 5/20/04)
!  (52) Added NVOC source diagnostics for ND07 (rjp, bmy, 7/13/04)
!  (53) Now reference "logical_mod.f", "tracer_mod.f", and "diag_pl_mod.f".
!        Bug fix in write to DMS_BIOG. (bmy, 7/20/04)
!  (54) Comment out ND27 for GEOS-4.  It isn't working 100% right.  If you
!        examine the flux at 200 hPa, you get the same info. (bmy, 10/15/04)
!  (55) Added biofuel SO4 to the bpch file under ND13.  Bug fix: replace ND68 
!        with LD68 in call to BPCH2 (auvray, bmy, 11/17/04)
!  (56) Now save ND03 mercury diagnostic arrays to bpch file.  Also updated
!        ND44 for tagged Hg tracers (eck, bmy, 12/14/04)
!  (57) Now print out extra ND21 diagnostics for crystalline sulfur tracers.  
!        Also now save total oceanic mass of Hg0 and Hg2.  Now call 
!        WRITE_DIAG03 from "diag03_mod.f" (bmy, 1/21/05)
!  (58) Now call WRITE_DIAG41 from "diag41_mod.f" (bmy, 2/17/05)
!  (59) Add P(SO4s) to row 8 of ND05 diagnostic.  Also remove special tracer
!        numbers for the ND67 diagnostic.  Now do not save CLDMAS for ND67
!        for GEOS-4, since GEOS-4 convection uses different met fields.
!        (bec, bmy, 5/3/05)
!  (60) Bug fix in ND68 diagnostic: use LD68 instead of ND68 in call to BPCH2.
!        Now modified for GEOS-5 and GCAP met fields.  Remove references to
!        CO-OH param simulation.  Also remove references to TRCOFFSET since
!        that is always zero now.  Now call GET_HALFPOLAR from "bpch2_mod.f" 
!        to get the HALFPOLAR value for GEOS or GCAP grids. (swu, bmy, 6/24/05)
!  (61) References ND04, WRITE_DIAG04 from "diag04_mod.f".  Also now updated
!        ND30 diagnostic for land/water/ice flags.  Also remove reference
!        to LWI array. (bmy, 8/18/05)
!  (62) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (63) Added MBO as tracer #5 in ND46 diagnostic (tmf, bmy, 10/20/05)
!  (64) Removed duplicate variable declarations.  Now remove restriction on 
!        printing out cloud mass flux in GEOS-4 for the ND66 diagnostic. 
!        (bmy, 3/14/06)
!  (65) References ND56, WRITE_DIAG56 from "diag56_mod.f" (ltm, bmy, 5/5/06)
!  (66) Now remove TRCOFFSET; it's obsolete.  References ND42, WRITE_DIAG42 
!        from "diag42_mod.f" (dkh, bmy, 5/22/06)
!  (67) Updated ND36 diagnostic for CH3I (bmy, 7/25/06)
!  (68) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (69) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (70) Now write diag 54 (time in the troposphere) if asked for (phs, 9/22/06)
!  (71) Now use new time counters for ND43 & ND45,  Also now average between
!        0 and 24 UT for ND47.  Bug fix in ND36. (phs, bmy, 3/5/07)
!  (72) Bug fix in ND65: use 3-D counter array (phs, bmy, 3/6/07)
!  (73) Bug fix in ND07: now save out IDTSOA4 tracer.  Modifications for H2/HD
!        diagnostics (ND10, ND27, ND44) (tmf, phs, bmy, 9/18/07)
!  (74) Now save out true pressure at 3-D level edges for ND31.  Change ND31
!        diagnostic category name to "PEDGE-$". Bug fix in ND28 diagnostic to 
!        allow you to print out individual biomass tracers w/o having to print 
!        all of them. (bmy, dkh, 1/24/08)
!  (75) Bug fix: Now divide ALBEDO in ND67 by SCALE_I6 for GEOS-3 met, but
!        by SCALE_A3 for all other met types (phs, bmy, 10/7/08)
!  (76) Fix ND65, ND47, and ozone case in ND45. Now only ND45 depends
!        on LD45 (phs, 11/17/08)
!  (77) Bug fix: Select the right index of AD34 to write.  Pick the right 
!         tracer field from AD22 if only a subset of tracers are requested 
!         to be printed out. (ccc, 12/15/08)
!  (78) Added ND52 for gamma(HO2) (jaegle, 02/26/09)
!  (79) Updated test on ship emissions flag for AD13 (phs, 3/3/09)     
!  (80) Add AD07_SOAGM for dicarbonyl SOA formation (tmf, 3/6/09)
!  (81) Add output in AD22 for dicarbonyl photolysis J values (tmf, 3/6/09)
!  (82) Add output in AD46 for biogenic C2H4 emissions (tmf, 3/6/09)
!  (83) Modify ND17, ND18, ND37, ND38, ND44 to output the tracers selected 
!        by the user. (ccc, 5/29/09)
!  (84) Add EFLUX output information for ND67. (lin, ccc, 5/29/09)
!  (85) Add test on ICOADS (cklee, 06/30/09)
!  (86) Add SCALE_DIAG to scale diagnostics with the number of accumulation 
!        steps. (ccc, 7/20/09)
!  (87) Add diagnostics 19, 58 and 60 for methane. (kjw, 8/18/09)
!  (88) Account for 3D AD13_NH3_an now (phs, 10/22/09)
!  (89) Added TOMAS diagnostics (win, bmy, 1/25/10)
!  (90) NBIOMAX is now in CMN_SIZE (hotp 7/31/09)
!  (91) Add SOA5 to ND07_HC, add AD57 for potential temperature. (fp, 2/3/10)
!  (92) Modify ND44 for tracers with several deposition tracers. (ccc, 2/3/10)
!  (93) Add aromatics to ND43. (dkh, 06/21/07)
!  (94) Add ND57 for potential temperature. (fp, 2/3/10)
!  (95) Re-order levels in mass fluxes diagnostics before writing them to file.
!       (ND24, 25, 26). (ccc, 3/8/10)
!  (96) Add call to update_dep for mercury simulation at the end.(ccc, 7/19/10)
!  20 Aug 2010 - R. Yantosca - Added ProTeX headers
!  20 Aug 2010 - R. Yantosca - Now pick proper scale for ND66 for MERRA
!  20 Aug 2010 - R. Yantosca - Now pick proper scale for ND67 for MERRA
!  20 Aug 2010 - R. Yantosca - Now added SCALE_A1 for hourly data
!  20 Aug 2010 - R. Yantosca - Now reference GET_A1_TIME from "time_mod.f"
!  26 May 2011 - R. Yantosca - For ND44, omit the special treatment of
!                              isoprene tracers if we are not doing fullchem
!  27 May 2011 - R. Yantosca - Now use SCALEDIAG for ND54 (time-in-trop) diag
!  08 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.x met
!  08 Feb 2012 - R. Yantosca - Restructure USE statements for clarity
!  08 Feb 2012 - R. Yantosca - Add counter for I3 (inst 3hr) met fields
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  05 Apr 2012 - R. Yantosca - Bug fix: use hourly scale for SLP in the
!                              ND67 diagnostic for GEOS-5.7.x met fields
!  14 Mar 2013 - M. Payer    - Replace NOx and Ox with NO, NO2, and O3 as part
!                              of removal of NOx-Ox partitioning
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  04 Sep 2013 - R. Yantosca - Make ND44 output consistent w/ modifications in
!                              GeosCore/gamap_mod.F.
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  03 Dec 2013 - R. Yantosca - Change unit of PBL height to meters, this used
!                              to be hPa in GEOS-1, GEOS-STRAT, GEOS-3, which
!                              are no longer supported.
!  28 Jan 2014 - R. Yantosca - Avoid array temporaries in ND60 TOMAS diagnostic
!  23 Jul 2014 - R. Yantosca - Remove reference to obsolete CMN_mod.F
!  23 Jul 2014 - R. Yantosca - Disable obsolete ND27 diagnostic
!  02 Aug 2014 - C. Keller   - Connect to HEMCO diagnostics
!  06 Aug 2014 - C. Keller   - Added wrapper subroutine for HEMCO diagnostics
!  14 Aug 2014 - R. Yantosca - Corrected units for several diagnostic outputs
!  26 Aug 2014 - M. Sulprizio- Now get ND53 POPs emissions from HEMCO
!  03 Sep 2014 - R. Yantosca - Units of AD01_Rn_SOURCE and AD01_Be7_SOURCE are
!                              now defined as kg/s in hcoi_gc_diagn_mod.F90 
!  15 Oct 2014 - C. Keller   - Updated ND37 diagnostics to write out all specs
!  10 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  24 Nov 2014 - R. Yantosca - Bug fix for ND28 BC/OC: Don't multiply by area
!  04 Dec 2014 - R. Yantosca - No longer crashes if BIOGENIC_OCPI=0 everywhere
!  12 Dec 2014 - M. Yannetti - Converted HEMCO Precision updates
!  15 Dec 2014 - M. Sulprizio- Moved radiation diagnostic from ND71 to ND72 to
!                              avoid conflicts with hourly max ppbv diagnostic.
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  23 Feb 2015 - C. Keller   - Split volcano SO2 into eruptive & degassing
!  02 Apr 2015 - M. Sulprizio- Now write out LIGHTNING_CLOUD_TOP diagnostic from
!                              HEMCO to the ctm.bpch file. Otherwise, if will
!                              save out to the HEMCO restart file.
!  22 Apr 2015 - M. Sulprizio- Remove the use of TRACER_COEFF to convert
!                              emissions of hydrocarbon species to atomsC/cm2/s.
!                              HEMCO now carries emissions for these species
!                              entirely in kgC/m2/s.
!  21 May 2015 - R. Yantosca - Bug fix: The bpch unit string for BIOMASS_NH3 
!                              should be molec/cm2/s instead of atoms C/cm2/2.
!  22 May 2015 - R. Yantosca - Remove variables made obsolete by HEMCO
!  22 May 2015 - R. Yantosca - Replace NOXEXTENT with 2 in TOMAS #if blocks
!  07 Jul 2015 - E. Lundgren - Add marine organic aerosols (B.Gantt, M.Johnson)
!  05 Jan 2016 - E. Lundgren - Use GC global param AVO instead of from HEMCO
!  21 Apr 2016 - R. Yantosca - Fixed typo that was preventing ND29 output
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  02 May 2016 - R. Yantosca - Now declare POPs and Rn tracer flags locally
!  31 May 2016 - E. Lundgren - Use species database, not XNUMOL, for molec wts
!  02 Jun 2016 - R. Yantosca - Now pass State_Met as an argument so that
!                              we can fix the units of ND62 INST-MAP 
!  16 Jun 2016 - R. Silvern  - Now define species ID flags locally witn Ind_
!  17 Jun 2016 - R. Yantosca - Renamed IDT* variables to id_* 
!  23 Jun 2016 - R. Yantosca - Bug fix: only print out CH4 diagnostics if
!                              it is an advected species
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!  18 Jul 2016 - R. Yantosca - Bug fix: look up TOMAS species ID's
!  11 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  20 Sep 2016 - R. Yantosca - Use .eqv. instead of .eq. to compare LOGICALs
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!  27 Mar 2017 - M. Sulprizio- Make anthropogenic emission diagnostics 3D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: NA, nAdvect, nDryDep, NW, NFAM
      INTEGER            :: I, IREF, J, JREF, L, M, MM, MMB, Levs
      INTEGER            :: N, NN, NMAX, NTEST, N_TOT_TRC, T
      INTEGER            :: IE, IN, IS, IW, ITEMP(3)
      INTEGER            :: NSPECOUT

      INTEGER            :: NN1    ! TOMAS tracers
      INTEGER            :: NBIN   ! TOMAS bin counter (win, 1/25/10)

      REAL(fp)           :: SCALE_TMP(IIPAR,JJPAR)
      REAL(fp)           :: SCALE_A1,   SCALE_A3  
      REAL(fp)           :: SCALED,     SCALEDYN
      REAL(fp)           :: SCALECONV,  SCALESRCE,  SCALECHEM  
      REAL(fp)           :: SCALEDIAG,  SCALE_ND66, SCALE_ND67 
      REAL(fp)           :: SCALERAD
      REAL(fp)           :: SCALEX,     SECONDS,    PMASS      
      REAL(fp)           :: PRESSX,     FDTT,       AREA_M2
      REAL(fp)           :: SCALE_I3
      REAL(f8)           :: DIAGb,      DIAGe
      
      ! For binary punch file, version 2.0
      CHARACTER (LEN=40) :: CATEGORY 
      REAL(f4)           :: ARRAY(IIPAR,JJPAR,LLPAR+1)
      REAL(f4)           :: LONRES, LATRES
      INTEGER            :: IFIRST, JFIRST, LFIRST
      INTEGER            :: HALFPOLAR
      INTEGER, PARAMETER :: CENTER180 = 1
      CHARACTER (LEN=20) :: MODELNAME 
      CHARACTER (LEN=40) :: UNIT
      CHARACTER (LEN=40) :: RESERVED = ''

#if   defined( TOMAS )
      ! For ND06 diagnostics
      CHARACTER(LEN=1)   :: ISTR1
      CHARACTER(LEN=2)   :: ISTR2

      ! For ND60 TOMAS diagnostic, avoids an array temporary (bmy, 1/28/14)
      REAL(f4)           :: ARR2D(JJPAR,LLPAR)
#endif

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER        :: Spc(:,:,:,:)
      REAL(fp), POINTER        :: AD(:,:,:)

      ! Interface w/ HEMCO diagnostics
      INTEGER                  :: FLAG
      INTEGER                  :: AFill
      CHARACTER(LEN= 63)       :: DiagnName, SrcName, SpcName, FullName
      CHARACTER(LEN=255)       :: MSG
      TYPE(DiagnCont), POINTER :: DiagnCnt
      REAL(fp)                 :: FACTOR, MW_G
      REAL(fp), PARAMETER      :: GPERKG   = 1000e+0_fp
      REAL(fp), PARAMETER      :: MWC      = 12e+0_fp ! hard-coded MW
      REAL(fp), PARAMETER      :: CM2PERM2 = 10000e+0_fp
      REAL(fp), PARAMETER      :: S_DMS    = 32e+0_fp / 62e+0_fp
      REAL(fp), PARAMETER      :: S_SO2    = 32e+0_fp / 64e+0_fp
      REAL(fp), PARAMETER      :: S_SO4    = 32e+0_fp / 96e+0_fp

      CHARACTER(LEN=255), PARAMETER :: LOC = 'DIAG3 (diag3.F)'

      ! RRTMG
      LOGICAL                  :: IS_GAS, IS_AD21, IS_DUST

      ! To point to the species database
      CHARACTER(LEN=31)        :: Name
      TYPE(Species), POINTER   :: SpcInfo

      ! Now define local tracer flags so that we can remove these from 
      ! tracerid_mod.F to facilitate FlexChem implementation (bmy, 5/2/16)
      LOGICAL, SAVE            :: FIRST = .TRUE.
      INTEGER, SAVE            :: id_Rn,       id_Pb
      INTEGER, SAVE            :: id_Be7
      INTEGER, SAVE            :: id_POPG
      INTEGER, SAVE            :: id_POPPOCPO, id_POPPOCPI
      INTEGER, SAVE            :: id_POPPBCPO, id_POPPBCPI
      INTEGER, SAVE            :: id_DST1,  id_DST2,  id_DST3,  id_DST4
      INTEGER, SAVE            :: id_DAL1,  id_DAL2,  id_DAL3,  id_DAL4
      INTEGER, SAVE            :: id_BCPI,  id_OCPI,  id_POA1,  id_MTPA
      INTEGER, SAVE            :: id_LIMO,  id_MTPO,  id_TSOA1, id_ISOA1
      INTEGER, SAVE            :: id_ASOA1, id_OPOA1, id_OPOG1
      INTEGER, SAVE            :: id_SALA,  id_SALC,  id_MOPO    
      INTEGER, SAVE            :: id_MOPI,  id_DMS,   id_SO2,   id_SO4
      INTEGER, SAVE            :: id_NH3,   id_NO,    id_CO,    id_ALK4
      INTEGER, SAVE            :: id_ACET,  id_MEK,   id_ALD2,  id_PRPE
      INTEGER, SAVE            :: id_C3H8,  id_CH2O,  id_C2H6,  id_CH4
      INTEGER, SAVE            :: id_ISOP,  id_C2H4,  id_CHBR3, id_BR2
      INTEGER, SAVE            :: id_DUST1, id_NK1,   id_SF1,   id_SS1
      INTEGER, SAVE            :: id_ECIL1, id_ECOB1, id_OCIL1, id_OCOB1
      INTEGER, SAVE            :: id_CH2BR2
      INTEGER, SAVE            :: id_PAN,   id_HNO3,  id_EOH,   id_MGLY
      INTEGER, SAVE            :: id_BENZ,  id_TOLU,  id_XYLE
      INTEGER, SAVE            :: id_NAP,   id_POG1,  id_POG2
 
!
!******************************************************************************
!  DIAG3 begins here!
!
!  Define scale factors for division.  
!  Add a small number (e.g. 1d-32) to prevent division by zero errors.
!******************************************************************************
!
      ! Assume success
      RC         = GC_SUCCESS

#if defined( BPCH_DIAG )

      ! Number of advected species
      nAdvect    = State_Chm%nAdvect

      ! Number of dry-deposited species
      nDryDep    = State_Chm%nDryDep

      ! Initialize
      SpcInfo    => NULL()
      DiagnCnt   => NULL() 

      ! Now use counter variables from "time_mod.f" (bmy, 3/27/03)
      DIAGb      = GET_DIAGb()
      DIAGe      = GET_DIAGe()
      SECONDS    = ( DIAGe - DIAGb ) * 3600e+0_fp
      SCALED     = 1e+0_fp
      SCALEDYN   = DBLE( GET_CT_DYN()  ) + 1e-32_fp
      SCALECONV  = DBLE( GET_CT_CONV() ) + 1e-32_fp
      SCALESRCE  = DBLE( GET_CT_EMIS() ) + 1e-32_fp
      SCALECHEM  = DBLE( GET_CT_CHEM() ) + 1e-32_fp
      SCALERAD   = DBLE( GET_CT_RAD()  ) + 1e-32_fp
      SCALE_A1   = DBLE( GET_CT_A1()   ) + 1e-32_fp
      SCALE_A3   = DBLE( GET_CT_A3()   ) + 1e-32_fp
      SCALE_I3   = DBLE( GET_CT_I3()   ) + 1e-32_fp
      SCALEDIAG  = DBLE( GET_CT_DIAG() ) + 1e-32_fp

!
!******************************************************************************
!  Now define local tracer flags for certain specialty simulations
!  so that we can remove them from tracerid_mod.F (bmy, 5/2/16)
!******************************************************************************
!
      IF ( FIRST ) THEN

         ! Initialize
         id_Rn        = Ind_('Rn'      )
         id_Pb        = Ind_('Pb'      )
         id_Be7       = Ind_('Be7'     )
         id_POPG      = Ind_('POPG'    )
         id_POPPOCPO  = Ind_('POPPOCPO')
         id_POPPOCPI  = Ind_('POPPOCPI')
         id_POPPBCPO  = Ind_('POPPBCPO')
         id_POPPBCPI  = Ind_('POPPBCPI')
         id_DST1      = Ind_('DST1'    )
         id_DST2      = Ind_('DST2'    )
         id_DST3      = Ind_('DST3'    )
         id_DST4      = Ind_('DST4'    )
         id_DAL1      = Ind_('DSTAL1'  )
         id_DAL2      = Ind_('DSTAL2'  )
         id_DAL3      = Ind_('DSTAL3'  )
         id_DAL4      = Ind_('DSTAL4'  )
         id_BCPI      = Ind_('BCPI'    )
         id_OCPI      = Ind_('OCPI'    )
         id_POA1      = Ind_('POA1'    )
         id_MTPA      = Ind_('MTPA'    )
         id_LIMO      = Ind_('LIMO'    )
         id_MTPO      = Ind_('MTPO'    )
         id_TSOA1     = Ind_('TSOA1'   )
         id_ISOA1     = Ind_('ISOA1'   )
         id_ASOA1     = Ind_('ASOA1'   )
         id_OPOA1     = Ind_('OPOA1'   )
         id_OPOG1     = Ind_('OPOG1'   )
         id_SALA      = Ind_('SALA'    )
         id_SALC      = Ind_('SALC'    )
         id_MOPO      = Ind_('MOPO'    ) 
         id_MOPI      = Ind_('MOPI'    )
         id_DMS       = Ind_('DMS'     )
         id_SO2       = Ind_('SO2'     )
         id_SO4       = Ind_('SO4'     )
         id_NH3       = Ind_('NH3'     )
         id_NO        = Ind_('NO'      )
         id_CO        = Ind_('CO'      )
         id_ALK4      = Ind_('ALK4'    )
         id_ACET      = Ind_('ACET'    )
         id_MEK       = Ind_('MEK'     )
         id_ALD2      = Ind_('ALD2'    )
         id_PRPE      = Ind_('PRPE'    )
         id_C3H8      = Ind_('C3H8'    )
         id_CH2O      = Ind_('CH2O'    )
         id_C2H6      = Ind_('C2H6'    )
         id_ISOP      = Ind_('ISOP'    )
         id_C2H4      = Ind_('C2H4'    )
         id_CHBR3     = Ind_('CHBR3'   )
         id_CH2BR2    = Ind_('CH2BR2'  )
         id_BR2       = Ind_('BR2'     )
         id_DUST1     = Ind_('DUST1'   )
         id_NK1       = Ind_('NK1'     )
         id_SF1       = Ind_('SF1'     )
         id_SS1       = Ind_('SS1'     )
         id_ECIL1     = Ind_('ECIL1'   )
         id_ECOB1     = Ind_('ECOB1'   )
         id_OCIL1     = Ind_('OCIL1'   )
         id_OCOB1     = Ind_('OCOB1'   )
         id_PAN       = Ind_('PAN'     )
         id_HNO3      = Ind_('HNO3'    )
         id_EOH       = Ind_('EOH'     )
         id_MGLY      = Ind_('MGLY'    )
         id_BENZ      = Ind_('BENZ'    )
         id_TOLU      = Ind_('TOLU'    )
         id_XYLE      = Ind_('XYLE'    )
         id_NAP       = Ind_('NAP'     )
         id_POG1      = Ind_('POG1'    )
         id_POG2      = Ind_('POG2'    )

         ! NOTE: CH4 can be an advected species (in CH4 or UCX-based sims),
         ! or a non-advected species (tropchem, soa, soa-svpoa).  We only
         ! want to print out diagnostics if CH4 is an advected species,
         ! so make sure to use the 'A' flag in the call to Ind_().
         ! (bmy, 6/23/16)
         id_CH4       = Ind_('CH4', 'A')

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF
!
!******************************************************************************
!  Setup for binary punch file:
!
!  IFIRST, JFIRST, LFIRST = I, J, L indices of the starting grid box 
!  LONRES                 = DISIZE, cast to REAL*4
!  LATRES                 = DJSIZE, cast to REAL*4
!******************************************************************************
!
      IFIRST = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LFIRST = 1
      LONRES = DISIZE
      LATRES = DJSIZE

      ! Get the proper model name and HALFPOLAR setting for the bpch file
      MODELNAME = GET_MODELNAME()
      HALFPOLAR = GET_HALFPOLAR()

      ! HEMCO interface: get pointer to HcoState object (of hcoi_gc_main_mod.F90)
      IF ( .NOT. ASSOCIATED(HcoState) ) THEN
         CALL ERROR_STOP ( 'HcoState not defined!', LOC ) 
      ENDIF
!
!******************************************************************************
!  ND01: Rn, Pb, Be emissions (Category: "RN--SRCE")
!
!   # : Field  : Description                    : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (1)  Rn222  : Emissions of 222Rn             : kg/s       : SCALESRCE
!  (2)  Pb210  : Emissions of 210Pb             : kg/s       : SCALECHEM
!  (3)  Be7    : Emissions of 7Be               : kg/s       : SCALESRCE
!
!  and  Rn, Pb, Be lost to radioactive decay (Category: "RN-DECAY")
!
!   # : Field  : Description                    : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (1)  Rn222  : Loss of 222Rn                  : kg/s       : SCALECHEM
!  (2)  Pb210  : Loss of 210Pb                  : kg/s       : SCALECHEM
!  (3)  Be7    : Loss of 7Be                    : kg/s       : SCALECHEM
!******************************************************************************
!
      IF ( ND01 > 0 ) THEN

         ! Diagnostic category
         CATEGORY  = 'RN--SRCE'
         UNIT      = 'kg/s'
         
         !-----------------------------------------------------------
         ! %%%%% Rn222 %%%%%
         !
         ! NOTE: This diagnostic only writes to the first level
         ! because the Rn source is only at the surface!!!
         !-----------------------------------------------------------
         IF ( id_Rn > 0 ) THEN
            DiagnName    = 'AD01_Rn_SOURCE'
            N            = id_Rn
            FACTOR       = 1.0e+0_fp
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,     N, 1, -1, .TRUE., FACTOR, RC    )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !-----------------------------------------------------------
         ! %%%%% Pb210 %%%%%
         !-----------------------------------------------------------
         IF ( id_Pb > 0 ) THEN

            ! Divide by # of chemistry timesteps
            DO L = 1, LD01
               ARRAY(:,:,L) = AD01(:,:,L) / SCALECHEM
            ENDDO

            ! Write to bpch file
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, id_Pb,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD01,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD01) )
         ENDIF

         !-----------------------------------------------------------
         ! %%%%% Be7 %%%%%
         !-----------------------------------------------------------
         IF ( id_Be7 > 0 ) THEN 
            DiagnName    = 'AD01_Be7_SOURCE'
            N            = id_Be7
            FACTOR       = 1.0e+0_fp
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT, N, 1, LD01, .TRUE., FACTOR, RC      )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF    
      ENDIF
!
!******************************************************************************
!  ND02: Rn, Pb, Be lost to radioactive decay (Category: "RN-DECAY")
!
!   # : Field  : Description                    : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (1)  Rn222  : Loss of 222Rn                  : kg/s       : SCALECHEM
!  (2)  Pb210  : Loss of 210Pb                  : kg/s       : SCALECHEM
!  (3)  Be7    : Loss of 7Be                    : kg/s       : SCALECHEM
!******************************************************************************
!
      IF ( ND02 > 0 ) THEN
         CATEGORY = 'RN-DECAY'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(2)
            N  = TINDEX(2,M)
            IF ( N > nAdvect ) CYCLE
            NN = N

            ! Divide by # of chemistry timesteps
            DO L = 1, LD02
               ARRAY(:,:,L) = AD02(:,:,L,N) / SCALECHEM
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD02,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD02) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND03: Diagnostics from Hg0/Hg2/HgP offline simulation (eck, bmy, 1/20/05)
!******************************************************************************
!
      IF ( ND03 > 0 ) THEN
         CALL WRITE_DIAG03( am_I_Root, Input_Opt, State_Chm, RC )
      ENDIF
!
!******************************************************************************
!  ND04: Diagnostics from CO2 simulation (pns, bmy, 7/26/05)
!******************************************************************************
!
      IF ( ND04 > 0 ) CALL WRITE_DIAG04
!
!******************************************************************************
!  ND05: Production/Loss for coupled fullchem/aerosol runs (NSRCX==3) or
!        offline sulfate chemistry runs (NSRCX==10).      
!
!   # : Field  : Description                     : Units        : Scale factor
!  ------------------------------------------------------------------------------
!  (1 ) SO2dms : P(SO2) from DMS + OH            : kg S         : SCALEX
!  (2 ) SO2no3 : P(SO2) from DMS + NO3           : kg S         : SCALEX
!  (3 ) SO2    : Total P(SO2)                    : kg S         : SCALEX
!  (4 ) MSAdms : P(MSA) from DMS                 : kg S         : SCALEX
!  (5 ) SO4gas : P(SO4) gas phase                : kg S         : SCALEX
!  (6 ) SO4H2O2: P(SO4) aqueous phase H2O2       : kg S         : SCALEX
!  (7 ) PSO4O3 : P(SO4) aqueous phase O3         : kg S         : SCALEX
!  (8 ) PSO4O2 : P(SO4) aqueous phase O2 (metal) : kg S         : SCALEX
!  (9 ) PSO4ss : Total P(SO4 from seasalt)       : kg S         : SCALEX
!  (10) LOH    : L(OH) by DMS                    : kg OH        : SCALEX
!  (11) LNO3   : L(NO3) by DMS                   : kg NO3       : SCALEX
!******************************************************************************
!
      IF ( ND05 > 0 ) THEN
         CATEGORY = 'PL-SUL=$'

         DO M = 1, TMAX(5)
            N = TINDEX(5,M)

            ! Tracers 9, 10 are OH, NO3
            ! and are in [kg] instead of [kg S]
            IF ( N < 10 ) THEN
               UNIT = 'kg S'
            ELSE
               UNIT = 'kg'
            ENDIF

            NN     = N
            SCALEX = 1.e+0_fp

            DO L = 1, LD05
               ARRAY(:,:,L) = AD05(:,:,L,N) / SCALEX
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD05,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD05) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND06: Dust aerosol emissions
!
!   # : Field    : Description                     : Units      : Scale factor
!  --------------------------------------------------------------------------
!  (1)  DUST     : Soil dust (4 different classes) : kg         : 1
 
!******************************************************************************
!
      IF ( ND06 > 0 .and. Input_Opt%LDUST .and. Input_Opt%LEMIS ) THEN

         ! Category & unit string
         UNIT     = 'kg'
         CATEGORY = 'DUSTSRCE'
         FACTOR   = 1.0e+0_fp

         ! Loop over # of dust bins (obtained from CMN_SIZE_mod)
         DO N = 1, NDSTBIN

#if defined( TOMAS )

            !----------------------------------------------------
            ! TOMAS simulations: many dust tracers
            !----------------------------------------------------

            ! Tracer number
            NN = id_DUST1 + ( N - 1 )

            ! Get TOMAS dust string
            IF ( N < 10 )  THEN
               WRITE( ISTR1,'(i1)' ) N
               DiagnName = 'AD06_DUST' // ISTR1
            ELSE
               WRITE( ISTR2,'(i2)' ) N
               DiagnName = 'AD06_DUST' // ISTR2
            ENDIF

#else

            !----------------------------------------------------
            ! Non-TOMAS simulations: 4 dust tracers
            !----------------------------------------------------
            SELECT CASE ( N )
               CASE ( 1 )
                  NN         = id_DST1 
                  DiagnName  = 'AD06_DST1'
               CASE ( 2 )
                  NN         = id_DST2 
                  DiagnName  = 'AD06_DST2'
               CASE ( 3 )
                  NN         = id_DST3
                  DiagnName  = 'AD06_DST3'
               CASE ( 4 )
                  NN         = id_DST4
                  DiagnName  = 'AD06_DST4'
            END SELECT
#endif

            ! Write HEMCO diagnostics to bpch file
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      NN, 1, -1, .FALSE.,  FACTOR, RC )

         ENDDO !N

         ! Include Dust alkalinity sources   tdf 6/18/2K8
         IF ( Input_Opt%LDSTUP ) THEN
            DO N = 5, NDSTBIN*2

               SELECT CASE ( N )
                  CASE ( 5 )
                     NN         = id_DAL1
                     DiagnName  = 'AD06_DSTAL1'
                  CASE ( 6 )
                     NN         = id_DAL2
                     DiagnName  = 'AD06_DSTAL2'
                  CASE ( 7 )
                     NN         = id_DAL3
                     DiagnName  = 'AD06_DSTAL3'
                  CASE ( 8 )
                     NN         = id_DAL4
                     DiagnName  = 'AD06_DSTAL4'
               END SELECT

               ! Write to HEMCO diagnostics to bpch file
               CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                         UNIT,      NN, 1, -1, .FALSE.,  FACTOR,
     &                         RC )
            ENDDO
         ENDIF

      ENDIF     
!
!******************************************************************************
!  ND07: Emissions of BC and OC aerosols
!
!   # : Field    : Description                  : Units        : Scale factor
!  --------------------------------------------------------------------------
!  (1)  Carbon   : Carbonaceous aerosols        : kg           : 1
!******************************************************************************
!
      IF ( ND07 > 0 .and. Input_Opt%LCARB ) THEN

         ! Unit
         ! updated units, OC and BC in kg C, SOA in kg (hotp 7/31/08)
         !UNIT = 'kg' 
         ! UNITS WERE CONFUSING -> everything in kgC (see latter)
         UNIT   = 'kgC' ! SOAupdate (hotp 7/29/10)
         FACTOR = SECONDS 
    
         !--------------------------------------------------------------
         ! Emission diagnostics (OCPI,OCPO,BCPI,BCPO).
         ! All carbon diagnostics are now in kgC/m2/s. Convert to kgC
         ! here. Also, HEMCO writes out diagnostics for hydrophilic and
         ! hydrophobic fractions individually. Add them together here.
         ! If the POA simulation is used, we need to replace OCPI/OCPO
         ! by POA1/POG1 because this is how those diagnostics become
         ! defined in hcoi_gc_diagn_mod.F
         !--------------------------------------------------------------

         ! Do for BC and OC
         DO J = 1, 2

            ! Select species index and name. The species name is also
            ! used to define the bpch category, so keep it OC for POA.
            IF ( J==1 ) THEN
               N       = id_BCPI
               SpcName = 'BC'
            ELSEIF ( J==2 ) THEN
               N = id_OCPI
               IF ( id_POA1 > 0 ) N = id_POA1
               SpcName = 'OC'
            ENDIF
  
            ! Do for all sources 
            DO I = 1, 3
               SELECT CASE ( I )
                  CASE ( 1 ) 
                     CATEGORY = TRIM(SpcName)//'-ANTH'
                     SRCNAME  = 'ANTHRO'
                     Levs     = LLPAR
                  CASE ( 2 ) 
                     CATEGORY = TRIM(SpcName)//'-BIOF'
                     SRCNAME  = 'BIOFUEL'
                     Levs     = 1
                  CASE ( 3 ) 
                     CATEGORY = TRIM(SpcName)//'-BIOB'
                     SRCNAME  = 'BIOMASS'
                     Levs     = 1
               END SELECT
   
               !----- (1) Hydrophilic -----
               ! Get the fullname. This is the name used in the
               ! HEMCO diagnostics.
               FullName = TRIM(SpcName)//'PI'
               IF ( N == id_POA1 ) FullName = 'POA1'

               ! Get diagnostics name
               IF ( I == 3 ) THEN
                  DiagnName='BIOMASS_'//TRIM(FullName)
               ELSE
                  DiagnName='AD07_'//TRIM(FullName)//'_'//TRIM(SRCNAME)
               ENDIF
  
               ! Get diagnostics from HEMCO
               CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                         FLAG, RC,   cName=TRIM(DiagnName), 
     &                         AutoFill=1, 
     &                         COL=HcoState%Diagn%HcoDiagnIDManual )
               IF ( RC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
                  CALL ERROR_STOP ( MSG, LOC ) 
               ENDIF
      
               ! Save into ARRAY and convert units from kg/m2/s to kg
               IF ( FLAG == HCO_SUCCESS ) THEN
                  IF ( I == 2 .or. I == 3 ) THEN
                     ARRAY(:,:,1) = DiagnCnt%Arr2D%Val(:,:) 
     &                            * FACTOR
     &                            * HcoState%Grid%AREA_M2%Val(:,:)
                  ELSE
                     ! Anthro emissions are 3D
                     DO L = 1, LLPAR
                        ARRAY(:,:,L) = DiagnCnt%Arr3D%Val(:,:,L)
     &                               * FACTOR
     &                               * HcoState%Grid%AREA_M2%Val(:,:)
                     ENDDO
                  ENDIF
               ELSE
                  MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
                  MSG = TRIM(MSG) // ' - will write zeros!'
                  CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
               ENDIF
      
               !----- (2) Add hydrophobic part -----

               ! Get the fullname. This is the name used in the
               ! HEMCO diagnostics.
               FullName = TRIM(SpcName)//'PO'
               IF ( N == id_POA1 ) FullName = 'POG1'
          
               ! Get diagnostics name
               IF ( I == 3 ) THEN
                  DiagnName='BIOMASS_'//TRIM(FullName)
               ELSE
                  DiagnName='AD07_'//TRIM(FullName)//'_'//TRIM(SRCNAME)
               ENDIF
   
               ! Get diagnostics from HEMCO
               DiagnCnt  => NULL()
               CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                         FLAG, RC,   cName=TRIM(DiagnName), 
     &                         AutoFill=1, 
     &                         COL=HcoState%Diagn%HcoDiagnIDManual )
               IF ( RC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
                  CALL ERROR_STOP ( MSG, LOC ) 
               ENDIF
      
               ! Save into ARRAY and convert units from kg/m2/s to kg
               IF ( FLAG == HCO_SUCCESS ) THEN
                  IF ( I == 2 .or. I == 3 ) THEN
                     ARRAY(:,:,1) = ARRAY(:,:,1) 
     &                            + ( DiagnCnt%Arr2D%Val(:,:)
     &                            * FACTOR 
     &                            * HcoState%Grid%AREA_M2%Val(:,:) )
                  ELSE
                     ! Anthro emissions are 3D
                     DO L = 1, LLPAR
                        ARRAY(:,:,L) = ARRAY(:,:,L) 
     &                               + ( DiagnCnt%Arr3D%Val(:,:,L)
     &                               *   FACTOR 
     &                               *   HcoState%Grid%AREA_M2%Val(:,:))
                     ENDDO 
                  ENDIF
               ELSE
                  MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
                  MSG = TRIM(MSG) // ' - will write zeros!'
                  CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
               ENDIF
      
               ! Write combined diagnostics to disk
               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, N,
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                     IIPAR,     JJPAR,     Levs,     IFIRST,     
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1:Levs) )
      
               ! Free pointer
               DiagnCnt => NULL()
            ENDDO !I
         ENDDO !J

         !--------------------------------------------------------------
         ! Other BC diagnostics
         !--------------------------------------------------------------

         !%%%%% Hydrophilic BC from Hydrophobic BC %%%%%
         CATEGORY     = 'PL-BC=$'
         N            = id_BCPI

         DO L = 1, LD07
            ARRAY(:,:,L) = AD07_BC(:,:,L) 
         ENDDO
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,
     &               IIPAR,     JJPAR,     LD07,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD07) )


         !--------------------------------------------------------------
         ! Other OC diagnostics
         !--------------------------------------------------------------

         ! These only exist with OCPI/OCPO
         IF ( id_OCPI > 0 ) THEN
            UNIT      = 'kgC'
            N         = id_OCPI

            !%%%%% Biogenic OC (kg/m2/s --> kg) %%%%%
            CATEGORY  = 'OC-BIOG'
            DiagnName = 'BIOGENIC_OCPI'
            FACTOR    = SECONDS

            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                 UNIT, N, 1, -1, .TRUE., FACTOR, RC, AreaScal=1 )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         
            !%%%%% Hydrophilic OC from hydrophobic OC %%%%%
            CATEGORY     = 'PL-OC=$'

            DO L = 1, LD07
               ARRAY(:,:,L) = AD07_OC(:,:,L) 
            ENDDO
         
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD07,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD07) )
         ENDIF

         !--------------------------------------------------------------
         ! Secondary organic aerosol (SOA) diagnostics
         !--------------------------------------------------------------
         IF ( Input_Opt%LSOA ) THEN

            ! Add units
            UNIT   = 'kg'
            FACTOR = SECONDS

            !------------------------------
            ! MTPA
            !------------------------------
            IF ( id_MTPA > 0 ) THEN
               CATEGORY  = 'OC-MTPA'
               NN        = id_MTPA

               ! (1) ALPHA-PINENE
               DiagnName = 'BIOGENIC_APIN'
               DiagnCnt  => NULL()

               ! Get diagnostics from HEMCO
               CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                         FLAG, RC,   cName=TRIM(DiagnName), 
     &                         AutoFill=1, 
     &                         COL=HcoState%Diagn%HcoDiagnIDManual )
               IF ( RC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
                  CALL ERROR_STOP ( MSG, LOC ) 
               ENDIF

               ! Save into ARRAY and convert units from kg/m2/s to kg
               IF ( FLAG == HCO_SUCCESS ) THEN
                  ARRAY(:,:,1) = DiagnCnt%Arr2D%Val(:,:)
     &                         * FACTOR * HcoState%Grid%AREA_M2%Val(:,:)
               ELSE
                  ARRAY(:,:,1) = 0.0
                  MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
                  MSG = TRIM(MSG) // ' - will write zeros!'
                  CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
               ENDIF

               ! (2) BETA-PINENE
               DiagnName = 'BIOGENIC_BPIN'
               DiagnCnt  => NULL()

               ! Get diagnostics from HEMCO
               CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                         FLAG, RC,   cName=TRIM(DiagnName), 
     &                         AutoFill=1, 
     &                         COL=HcoState%Diagn%HcoDiagnIDManual )
               IF ( RC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
                  CALL ERROR_STOP ( MSG, LOC ) 
               ENDIF

               ! Save into ARRAY and convert units from kg/m2/s to kg
               IF ( FLAG == HCO_SUCCESS ) THEN
                  ARRAY(:,:,1) = ARRAY(:,:,1)
     &                         + ( DiagnCnt%Arr2D%Val(:,:) * FACTOR
     &                         *   HcoState%Grid%AREA_M2%Val(:,:) )
               ELSE
                  MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
                  MSG = TRIM(MSG) // ' - will write zeros!'
                  CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
               ENDIF

               ! (3) Sabinene
               DiagnName = 'BIOGENIC_SABI'
               DiagnCnt  => NULL()

               ! Get diagnostics from HEMCO
               CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                         FLAG, RC,   cName=TRIM(DiagnName), 
     &                         AutoFill=1, 
     &                         COL=HcoState%Diagn%HcoDiagnIDManual )
               IF ( RC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
                  CALL ERROR_STOP ( MSG, LOC ) 
               ENDIF

               ! Save into ARRAY and convert units from kg/m2/s to kg
               IF ( FLAG == HCO_SUCCESS ) THEN
                  ARRAY(:,:,1) = ARRAY(:,:,1)
     &                         + ( DiagnCnt%Arr2D%Val(:,:) * FACTOR
     &                         *   HcoState%Grid%AREA_M2%Val(:,:) )
               ELSE
                  MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
                  MSG = TRIM(MSG) // ' - will write zeros!'
                  CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
               ENDIF

               ! (4) CARENE
               DiagnName = 'BIOGENIC_CARE'
               DiagnCnt  => NULL()

               ! Get diagnostics from HEMCO
               CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                         FLAG, RC,   cName=TRIM(DiagnName), 
     &                         AutoFill=1, 
     &                         COL=HcoState%Diagn%HcoDiagnIDManual )
               IF ( RC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
                  CALL ERROR_STOP ( MSG, LOC ) 
               ENDIF

               ! Save into ARRAY and convert units from kg/m2/s to kg
               IF ( FLAG == HCO_SUCCESS ) THEN
                  ARRAY(:,:,1) = ARRAY(:,:,1)
     &                         + ( DiagnCnt%Arr2D%Val(:,:) * FACTOR
     &                         *   HcoState%Grid%AREA_M2%Val(:,:) )
               ELSE
                  MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
                  MSG = TRIM(MSG) // ' - will write zeros!'
                  CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
               ENDIF

               ! Write combined BIOGENIC MTPA to disk
               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,
     &                     UNIT ,     DIAGb,     DIAGe,    RESERVED,   
     &                     IIPAR,     JJPAR,     1,        IFIRST,     
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1) )
            ENDIF

            !------------------------------
            ! LIMO (kg/m2/s --> kg)
            !------------------------------
            IF ( id_LIMO > 0 ) THEN
               CATEGORY  = 'OC-LIMO'
               NN        = id_LIMO
               DiagnName = 'BIOGENIC_LIMO'
               CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                         UNIT,      NN, 1, -1, .TRUE.,   FACTOR,
     &                         RC,        AreaScal=1 )
               IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
            ENDIF

            !------------------------------
            ! MTPO
            !------------------------------
            IF ( id_MTPO > 0 ) THEN
               CATEGORY  = 'OC-MTPO'
               NN        = id_MTPO

               ! (1) MYRCENE
               DiagnName = 'BIOGENIC_MYRC'
               DiagnCnt  => NULL()

               ! Get diagnostics from HEMCO
               CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                         FLAG, RC,   cName=TRIM(DiagnName), 
     &                         AutoFill=1, 
     &                         COL=HcoState%Diagn%HcoDiagnIDManual )
               IF ( RC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
                  CALL ERROR_STOP ( MSG, LOC ) 
               ENDIF

               ! Save into ARRAY and convert units from kg/m2/s to kg
               IF ( FLAG == HCO_SUCCESS ) THEN
                  ARRAY(:,:,1) = DiagnCnt%Arr2D%Val(:,:)
     &                         * FACTOR * HcoState%Grid%AREA_M2%Val(:,:)
               ELSE
                  ARRAY(:,:,1) = 0.0
                  MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
                  MSG = TRIM(MSG) // ' - will write zeros!'
                  CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
               ENDIF

               ! (2) OCIMENE
               DiagnName = 'BIOGENIC_OCIM'
               DiagnCnt  => NULL()

               ! Get diagnostics from HEMCO
               CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                         FLAG, RC,   cName=TRIM(DiagnName), 
     &                         AutoFill=1, 
     &                         COL=HcoState%Diagn%HcoDiagnIDManual )
               IF ( RC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
                  CALL ERROR_STOP ( MSG, LOC ) 
               ENDIF

               ! Save into ARRAY and convert units from kg/m2/s to kg
               IF ( FLAG == HCO_SUCCESS ) THEN
                  ARRAY(:,:,1) = ARRAY(:,:,1)
     &                         + ( DiagnCnt%Arr2D%Val(:,:) * FACTOR
     &                         *   HcoState%Grid%AREA_M2%Val(:,:) )
               ELSE
                  MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
                  MSG = TRIM(MSG) // ' - will write zeros!'
                  CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
               ENDIF

               ! (3) OTHER MONOTERPENES
               DiagnName = 'BIOGENIC_OMON'
               DiagnCnt  => NULL()

               ! Get diagnostics from HEMCO
               CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                         FLAG, RC,   cName=TRIM(DiagnName), 
     &                         AutoFill=1, 
     &                         COL=HcoState%Diagn%HcoDiagnIDManual )
               IF ( RC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
                  CALL ERROR_STOP ( MSG, LOC ) 
               ENDIF

               ! Save into ARRAY and convert units from kg/m2/s to kg
               IF ( FLAG == HCO_SUCCESS ) THEN
                  ARRAY(:,:,1) = ARRAY(:,:,1)
     &                         + ( DiagnCnt%Arr2D%Val(:,:) * FACTOR
     &                         *   HcoState%Grid%AREA_M2%Val(:,:) )
               ELSE
                  MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
                  MSG = TRIM(MSG) // ' - will write zeros!'
                  CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
               ENDIF

               ! Write combined BIOGENIC MTPO to disk
               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,
     &                     UNIT ,     DIAGb,     DIAGe,    RESERVED,   
     &                     IIPAR,     JJPAR,     1,        IFIRST,     
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1) )
            ENDIF

            !------------------------------
            ! SESQ
            !------------------------------
            IF ( id_LIMO > 0 ) THEN
               CATEGORY  = 'OC-SESQ'
               NN        = id_LIMO + 3

               ! (1) FARNESENE
               DiagnName = 'BIOGENIC_FARN'
               DiagnCnt  => NULL()

               ! Get diagnostics from HEMCO
               CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                         FLAG, RC,   cName=TRIM(DiagnName), 
     &                         AutoFill=1, 
     &                         COL=HcoState%Diagn%HcoDiagnIDManual )
               IF ( RC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
                  CALL ERROR_STOP ( MSG, LOC ) 
               ENDIF

               ! Save into ARRAY and convert units from kg/m2/s to kg
               IF ( FLAG == HCO_SUCCESS ) THEN
                  ARRAY(:,:,1) = DiagnCnt%Arr2D%Val(:,:)
     &                         * FACTOR * HcoState%Grid%AREA_M2%Val(:,:)
               ELSE
                  ARRAY(:,:,1) = 0.0
                  MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
                  MSG = TRIM(MSG) // ' - will write zeros!'
                  CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
               ENDIF

               ! (2) BETA-CARYOPHYLLENE
               DiagnName = 'BIOGENIC_BCAR'
               DiagnCnt  => NULL()

               ! Get diagnostics from HEMCO
               CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                         FLAG, RC,   cName=TRIM(DiagnName), 
     &                         AutoFill=1, 
     &                         COL=HcoState%Diagn%HcoDiagnIDManual )
               IF ( RC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
                  CALL ERROR_STOP ( MSG, LOC ) 
               ENDIF

               ! Save into ARRAY and convert units from kg/m2/s to kg
               IF ( FLAG == HCO_SUCCESS ) THEN
                  ARRAY(:,:,1) = ARRAY(:,:,1)
     &                         + ( DiagnCnt%Arr2D%Val(:,:) * FACTOR
     &                         *   HcoState%Grid%AREA_M2%Val(:,:) )
               ELSE
                  MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
                  MSG = TRIM(MSG) // ' - will write zeros!'
                  CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
               ENDIF

               ! (3) OTHER SESQUITERPENES
               DiagnName = 'BIOGENIC_OSQT'
               DiagnCnt  => NULL()

               ! Get diagnostics from HEMCO
               CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                         FLAG, RC,   cName=TRIM(DiagnName), 
     &                         AutoFill=1, 
     &                         COL=HcoState%Diagn%HcoDiagnIDManual )
               IF ( RC /= HCO_SUCCESS ) THEN
                  MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
                  CALL ERROR_STOP ( MSG, LOC ) 
               ENDIF

               ! Save into ARRAY and convert units from kg/m2/s to kg
               IF ( FLAG == HCO_SUCCESS ) THEN
                  ARRAY(:,:,1) = ARRAY(:,:,1)
     &                         + ( DiagnCnt%Arr2D%Val(:,:) * FACTOR
     &                         *   HcoState%Grid%AREA_M2%Val(:,:) ) 
               ELSE
                  MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
                  MSG = TRIM(MSG) // ' - will write zeros!'
                  CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
               ENDIF

               ! Write combined BIOGENIC SESQ to disk
               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,
     &                     UNIT ,     DIAGb,     DIAGe,    RESERVED,   
     &                     IIPAR,     JJPAR,     1,        IFIRST,     
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1) )
            ENDIF
         
            !-----------------------------------------------
            ! SOA Production from NVOC oxidation [kg]
            ! 1:ALPH+LIMO+TERP, 2:ALCO, 3:SESQ, 4:ISOP
            ! 5:AROM (Add 4 and 5 (dkh, 03/27/07)  ))
            !-----------------------------------------------
            CATEGORY = 'PL-OC=$'

            ! SOAupdate: hotp 5/24/10 (new mtp)
            DO N = 1, 6

               ! hotp 7/31/08 units
               UNIT = 'kg' 

               ! SOAupdate: semivolpoa2: units diff for POA (hotp 3/4/09)
               ! now correspond to JSV (hotp 5/24/10) (new mtp)
               IF     ( N == 1 ) THEN
                  NN = id_TSOA1
               ELSEIF ( N == 2 ) THEN
                  NN = id_ISOA1
               ELSEIF ( N == 3 ) THEN
                  NN = id_ASOA1
               ELSEIF ( N == 4 ) THEN
                  NN = id_POA1
                  UNIT = 'kgC'
               ELSEIF ( N == 5 ) THEN
                  NN = id_OPOA1
                  UNIT = 'kgC'
               ELSEIF ( N == 6 ) THEN
                  NN = id_OPOG1
                  UNIT = 'kgC'
               ENDIF

               ! Check to make sure the tracer exists (hotp 8/23/09)
               IF ( NN > 0 ) THEN
                  DO L = 1, LD07
                     ARRAY(:,:,L) = AD07_HC(:,:,L,N)
                  ENDDO

                  CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                        HALFPOLAR, CENTER180, CATEGORY, NN,
     &                        UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                        IIPAR,     JJPAR,     LD07,     IFIRST,
     &                        JFIRST,    LFIRST,    ARRAY(:,:,1:LD07) )
               ENDIF ! NN

            ENDDO

         ENDIF
      ENDIF   
!
!******************************************************************************
!  ND08: Sea salt aerosol emissions
!
!   # : Field    : Description                     : Units      : Scale factor
!  --------------------------------------------------------------------------
!  (1)  SALA     : Accumulation mode seasalt       : kg         : 1
!  (2)  SALC     : Coarse mode seasalt             : kg         : 1
!  (3)  MOPO     : Marine organic aerosol - phobo  : kg         : 1
!  (4)  MOPI     : Marine organic aerosol - philic : kg         : 1
!******************************************************************************
!
      IF ( ND08 > 0 .and. Input_Opt%LSSALT .and. Input_Opt%LEMIS ) THEN

         ! Category & unit string
         UNIT     = 'kg'
         CATEGORY = 'SALTSRCE'

         ! Loop over default seasalt tracers
         DO N = 1, 2

            ! At present we have 2 default seasalts
            IF ( N == 1 ) THEN
               NN = id_SALA
               DiagnName  = 'AD08_SALA'
            ELSEIF ( N == 2 ) THEN
               NN = id_SALC
               DiagnName  = 'AD08_SALC'
            ENDIF

            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                      UNIT, NN, 1, -1, .FALSE., 1.0e+0_fp, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         ENDDO

         ! Loop over optional marine organic aerosols
         IF ( Input_Opt%LMPOA ) THEN
            DO N = 3, 4
   
               IF ( N == 3 ) THEN
                  NN = id_MOPO
                  DiagnName  = 'AD08_MOPO'
               ELSEIF ( N == 4 ) THEN
                  NN = id_MOPI
                  DiagnName  = 'AD08_MOPI'
               ENDIF
   
               CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                         UNIT, NN, 1, -1, .FALSE., 1.0e+0_fp, RC )
               IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP (DiagnName, LOC)
   
            ENDDO
         ENDIF

      ENDIF   
!
!******************************************************************************
!  ND09: HCN/CH3CN sources/sinks (Categories: "HCN-PL-$", "HCN-SRCE")
!
!  # : Field    : Description                     : Units       : Scale factor
! ----------------------------------------------------------------------------
! (1:N) sink    : Loss of tagged tracer to OH     : kg
! (N+1) HCNbb   : HCN   from biomass burning      : molec/cm2/s : SCALESRCE
! (N+2) CH3CNbb : CH3CN from biomass burning      : molec/cm2/s : SCALESRCE
! (N+3) HCNdf   : HCN   from domestic fossil fuel : molec/cm2/s : SCALESRCE
! (N+4) CH3CNdf : CH3CN from domestic fossil fuel : molec/cm2/s : SCALESRCE
! (N+5) HCNoc   : HCN   loss to ocean uptake      : molec/cm2/s : SCALECHEM
! (N+6) CH3CNoc : CH3CN loss to ocean uptake      : molec/cm2/s : SCALECHEM
!******************************************************************************
!
      IF ( ND09 > 0 ) THEN

         ! Binary punch file
         DO M = 1, TMAX(9)
            N  = TINDEX(9,M)
            IF ( N > nAdvect+6 ) CYCLE

            ! Test tracer number
            IF ( N <= nAdvect ) THEN

               !---------------------------
               ! HCN/CH3CN sinks
               !---------------------------              
               CATEGORY  = 'HCN-PL-$'
               UNIT      = 'kg'
               NN        = N
              
               DO L = 1, LD09
                  ARRAY(:,:,L) = AD09(:,:,L,N)
               ENDDO

               ! Save to disk
               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                     IIPAR,     JJPAR,     LD09,     IFIRST,
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1:LD09) )

            ELSE 

               !---------------------------
               ! HCN/CH3CN sources
               !---------------------------
               CATEGORY     = 'HCN-SRCE'
               UNIT         = 'molec/cm2/s'
               NN           = N - nAdvect

               ! Pick proper scale
               IF ( NN <= 4 ) THEN
                  SCALEX = SCALESRCE
               ELSE
                  SCALEX = SCALECHEM
               ENDIF

               ! Scale data
               ARRAY(:,:,1) = AD09_em(:,:,NN) / SCALEX
 
               ! Write to disk
               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                     IIPAR,     JJPAR,     1,        IFIRST,
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1) )
            ENDIF
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND11: Acetone source & sink diagnostic (Category: "ACETSRCE")
!
!   # : Field  : Description                        : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (1)  ACETmo : Acetone source from MONOTERPENES   : at C/cm2/s : SCALESRCE
!  (2)  ACETmb : Acetone source from METHYL BUTENOL : at C/cm2/s : SCALESRCE 
!  (3)  ACETbg : Acetone source from DIRECT EMISSION: at C/cm2/s : SCALESRCE 
!  (4)  ACETop : Acetone source from OCEANS         : at C/cm2/s : SCALESRCE 
!  (5)  ACETol : Acetone sink   from OCEANS         : at C/cm2/s : SCALECHEM
!
!  NOTES:
!  (1 ) in HEMCO, Acetone sink from OCEAN is handled by drydep and the 
!       corresponding diagnostics. If needed, we can write a wrapper that
!       explicitly calculates the Acetone ocean sink. (ckeller, 08/04/14) 
!******************************************************************************
!
      IF ( ND11 > 0 .and. Input_Opt%LEMIS ) THEN
         CATEGORY = 'ACETSRCE'
         UNIT     = 'atoms C/cm2/s'

         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR   = GPERKG / MWC * AVO / 
     &              CM2PERM2

         ! Loop over all four emission sources
         DO N = 1, 4 
            SELECT CASE( N )
               CASE( 1 ) 
                  DiagnName = 'MEGAN_ACET_MONO'
                  AFill     = 0
               CASE( 2 )
                  DiagnName = 'MEGAN_ACET_MBO'
                  AFill     = 0
               CASE( 3 )
                  DiagnName = 'MEGAN_ACET_DIRECT'
                  AFill     = 0
               CASE( 4 )
                  DiagnName = 'AD11_OCEAN_SOURCE'
                  AFill     = 1
            END SELECT

            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                      UNIT, N, AFill, -1, .FALSE., FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         ENDDO
      ENDIF
!
!******************************************************************************
!  ND12: distribution of surface emissions in the boundry layer: [fraction]

!   # : Field   : Description                         : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1) EMDIS-BL : Fraction of BL occupied by level L  : unitless : SCALECHEM
!******************************************************************************
!
      IF ( ND12 > 0 ) THEN
         UNIT     = 'unitless'
         CATEGORY = 'EMDIS-BL'

         DO L = 1, LD12
            ARRAY(:,:,L) = AD12(:,:,L) / SCALECHEM
         ENDDO

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, 1,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LLCHEM,   IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD12) )
      ENDIF
!
!******************************************************************************
!  ND13: Sulfur emissions (for DMS/SO2/SO4/MSA/NH3/NH4/NIT chemistry)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) DMS-BIOG : Biogenic DMS emission           : kg S     : 1
!  (2 ) SO2-AC-$ : Aircraft SO2 emission           : kg S     : 1
!  (3 ) SO2-AN-$ : Anthropogenic SO2 emission      : kg S     : 1
!  (4 ) SO2-BIOB : Biomass SO2 emission            : kg S     : 1
!  (5 ) SO2-BIOF : Biofuel SO2 emission            : kg S     : 1
!  (6 ) SO2-NV-$ : Non-eruptive volcano SO2 em.    : kg S     : 1
!  (7 ) SO2-EV-$ : Eruptive volcano SO2 emissions  : kg S     : 1
!  (8 ) SO4-AN-$ : Anthropogenic SO4 emission      : kg S     : 1
!  (9 ) NH3-ANTH : Anthropogenic NH3 emission      : kg NH3   : 1
!  (10) NH3-NATU : Natural source NH3 emission     : kg NH3   : 1
!  (11) NH3-BIOB : Biomass burning NH3 emission    : kg NH3   : 1
!  (12) NH3-BIOF : Biofuel burning NH3 emission    : kg NH3   : 1
!
!  NOTES:
!  (1) ND13 diagnostics are now tracked by HEMCO. (bmy, 8/18/14)
!******************************************************************************
!
      IF ( ND13 > 0 .and. Input_Opt%LEMIS .and. 
     &   ( Input_Opt%ITS_A_FULLCHEM_SIM .or. 
     &     Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
         UNIT = 'kg S'

         !--------------------------------------------------------------
         ! DMS diagnostics
         !--------------------------------------------------------------

         !%%%%% Biogenic DMS %%%%%
         DiagnName = 'AD13_DMS_OCEAN_SOURCE'
         CATEGORY  = 'DMS-BIOG'
         FACTOR    = S_DMS
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                   UNIT, id_DMS, 1, -1, .FALSE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         !--------------------------------------------------------------
         ! SO2 diagnostics
         !--------------------------------------------------------------

         !%%%%% Aircraft SO2 %%%%%
         DiagnName = 'AD13_SO2_AIRCRAFT'
         CATEGORY = 'SO2-AC-$'
         FACTOR   = S_SO2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                   UNIT, id_SO2, 1, LLPAR, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         !%%%%% Anthro SO2 %%%%%
         DiagnName = 'AD13_SO2_ANTHROPOGENIC'
         CATEGORY  = 'SO2-AN-$'
         FACTOR    =  S_SO2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                   UNIT, id_SO2, 1, LLPAR, .FALSE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )
         
         !%%%%% Biomass SO2 (convert kg/m2/s --> kg S) %%%%%
         DiagnName = 'BIOMASS_SO2'
         CATEGORY  = 'SO2-BIOB'
         FACTOR    = S_SO2 * SECONDS
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &           UNIT, id_SO2, 1, -1, .TRUE., FACTOR, RC, AreaScal=1 )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         !%%%%% Biofuel SO2 %%%%%
         DiagnName = 'AD13_SO2_BIOFUEL'
         CATEGORY  = 'SO2-BIOF'
         FACTOR    = S_SO2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                   UNIT, id_SO2, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         !%%%%% Volcano SO2 (eruptive) %%%%%
         DiagnName = 'AD13_SO2_VOLCANO_ERUPT'
         CATEGORY  = 'SO2-EV-$'
         FACTOR    = S_SO2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                   UNIT, id_SO2, 1, LLPAR, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         !%%%%% Volcano SO2 (degassing) %%%%%
         DiagnName = 'AD13_SO2_VOLCANO_DEGAS'
         CATEGORY  = 'SO2-NV-$'
         FACTOR    = S_SO2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                   UNIT, id_SO2, 1, LLPAR, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         !%%%%% Ship SO2 %%%%%
         DiagnName = 'AD13_SO2_SHIP'
         CATEGORY  = 'SO2-SHIP'
         FACTOR    = S_SO2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                   UNIT, id_SO2, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         !--------------------------------------------------------------
         ! SO4 diagnostics
         !--------------------------------------------------------------

         !%%%%% Anthropogenic SO4 %%%%%
         DiagnName = 'AD13_SO4_ANTHROPOGENIC'
         CATEGORY  = 'SO4-AN-$'
         FACTOR    = S_SO4
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                   UNIT, id_SO4, 1, LLPAR, .FALSE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         !%%%%% Biofuel SO4 %%%%%
         DiagnName = 'AD13_SO4_BIOFUEL'
         CATEGORY  = 'SO4-BIOF'
         FACTOR    = S_SO4
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                   UNIT, id_SO4, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         !--------------------------------------------------------------
         ! NH3 diagnostics
         !--------------------------------------------------------------

         !%%%%% Anthropogenic NH3 %%%%%
         DiagnName = 'AD13_NH3_ANTHROPOGENIC'
         UNIT      = 'kg'
         CATEGORY  = 'NH3-ANTH'
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                   UNIT, id_NH3, 1, LLPAR, .FALSE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         !%%%%%  Natural source NH3 %%%%%
         DiagnName = 'AD13_NH3_NATURAL'
         CATEGORY  = 'NH3-NATU'
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                   UNIT, id_NH3, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         !%%%%% Biomass NH3 (convert kg/m2/s --> kg) %%%%%
         DiagnName = 'BIOMASS_NH3'
         CATEGORY  = 'NH3-BIOB'
         FACTOR    = SECONDS
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &          UNIT, id_NH3, 1, -1, .TRUE., FACTOR, RC, AreaScal=1 )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

         !%%%%% Biofuel NH3 %%%%%
         DiagnName = 'AD13_NH3_BIOFUEL'
         CATEGORY  = 'NH3-BIOF'
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY, 
     &                   UNIT, id_NH3, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( DiagnName, LOC )

      ENDIF    
!
!******************************************************************************
!  ND14: Upward mass flux from wet convection (NFCLDMX)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  CONVFLUP : Upward mass flux from wet conv  : kg/s     : SCALECONV
!
!  NOTES:
!  (1) Bug fix -- only write LD14 levels to the bpch file (bmy, 12/7/00)
!******************************************************************************
!
      IF ( ND14 > 0 ) THEN
         CATEGORY = 'CV-FLX-$'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(14)
            N  = TINDEX(14,M)
            IF ( N > nAdvect ) CYCLE
            NN = N
               
            ARRAY(:,:,1:LD14) = CONVFLUP(:,:,1:LD14,N) / SCALECONV

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD14,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD14) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND15: Upward mass flux from boundary layer mixing (TURBDAY)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  TURBFLUX : Upward mass flux from BL mixing : kg/s     : SCALECONV
!
!  NOTES:
!  (1) Bug fix -- only write LD15 levels to the bpch file (bmy, 12/7/00)
!******************************************************************************
!
      IF ( ND15 > 0 ) THEN
         CATEGORY = 'TURBMC-$'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(15)
            N  = TINDEX(15,M)
            IF ( N > nAdvect ) CYCLE
            NN = N

            ARRAY(:,:,1:LD15) = TURBFLUP(:,:,1:LD15,N) / SCALECONV
               
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD15,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD15) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND16: Fraction of grid box experiencing LS or convective precipitation
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WD-FLS-$ : LS precip fraction          : unitless  : CT16(:,:,:,1)
!  (2) WD-FCV-$ : Convective precip fraction  : unitless  : CT16(:,:,:,2)
!******************************************************************************
!
      IF ( ND16 > 0 ) THEN

         ! Large-scale area of precipitation
         CATEGORY = 'WD-FRC-$'
         UNIT     = 'unitless'

         ! GEOS-FP/MERRA-2 only has LS wetdep, so only need to save N=1
         ! (bmy, 11/16/17)
         DO L = 1, LD16
            SCALE_TMP(:,:) = FLOAT( CT16(:,:,L,1) ) + 1e-20_fp
            ARRAY(:,:,L)   = AD16(:,:,L,1) / SCALE_TMP(:,:)
         ENDDO

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, 1,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD16,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD16) )

      ENDIF
!
!******************************************************************************
!  ND17: Fraction of tracer lost rainout in LS and convective precip
!
!   #  Field    : Description                  : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WD-LSR-$ : Rainout fraction/LS Precip   : unitless  : CT17(:,:,:,1)
!
!  NOTES:
!  (1) Now loop over all soluble tracers (bmy, 3/19/04)
!  (2) Now use actual tracer number (bmy, 3/23/04)
!******************************************************************************
!
      IF ( ND17 > 0 ) THEN
         UNIT = 'unitless'

         ! Loop over soluble species
         DO NW = 1, State_Chm%nWetDep

            ! Get the species ID from the wetdep ID
            N = State_Chm%Map_WetDep(NW)

            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= N .AND. MM <= TMAX(17) )
               MMB = TINDEX(17,MM)
               MM  = MM + 1
            ENDDO
            
            IF ( MMB /= N ) CYCLE

            ! Large-scale rainout/washout fractions
            CATEGORY = 'WD-LSR-$'
               
            DO L = 1, LD17
               SCALE_TMP(:,:) = FLOAT( CT17(:,:,L,1) ) + 1e-20_fp
               ARRAY(:,:,L)   = AD17(:,:,L,NW,1) / SCALE_TMP(:,:) 
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD17,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD17) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND18: Fraction of tracer lost to washout in LS or convective precip
!
!   #  Field    : Description                  : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WD-LSW-$ : Washout fraction/LS precip   : unitless  : CT18(:,:,:,1)
!
!  NOTES:
!  (1) Now loop over all soluble tracers (bmy, 3/19/04)
!  (2) Now use actual tracer number (bmy, 3/23/04)
!******************************************************************************
!
      IF ( ND18 > 0 ) THEN
         UNIT = 'unitless'

         ! Loop over the soluble species
         DO NW = 1, State_Chm%nWetDep

            ! Get the species ID from the wetdep ID
            N = State_Chm%Map_WetDep(NW)

            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= N .AND. MM <= TMAX(18) )
               MMB = TINDEX(18,MM)
               MM  = MM + 1
            ENDDO
            
            IF ( MMB /= N ) CYCLE

            ! Large-scale rainout/washout fractions
            CATEGORY = 'WD-LSW-$'

            DO L = 1, LD18
               SCALE_TMP(:,:) = FLOAT( CT18(:,:,L,1) ) + 1e-20_fp
               ARRAY(:,:,L)   = AD18(:,:,L,NW,1) / SCALE_TMP(:,:) 
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD18,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD18) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND19: CH4 loss
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) CH4-LOSS : CH4 loss by reaction with OH    : kg CH4   : 1
!  (2 ) CH4-LOSS : CH4 loss in the stratosphere    : kg CH4   : 1
!  (3 ) CH4-LOSS : CH4 loss by reaction with Cl    : kg CH4   : 1
!******************************************************************************
!
      IF ( ND19 > 0 ) THEN
         UNIT     = 'kg'
         CATEGORY = 'CH4-LOSS'

         DO M = 1, TMAX(19)
            N  = TINDEX(19,M)
            IF ( N > PD19 ) CYCLE
            NN = N

            ARRAY(:,:,1:LD19) = AD19(:,:,1:LD19,N)

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD19,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD19) )

         ENDDO
      ENDIF
!
!******************************************************************************
!  ND21: Optical depth diagnostics
!
!   # : Field : Description                        : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) OPTD   Cloud Optical Depth                 : unitless : SCALECHEM
!  (2 ) CLDTOT 3-D Total Cloud fraction            : unitless : SCALECHEM
!  (3 ) CLRO   Random  Overlap Cloud Fraction      : unitless : SCALECHEM
!  (4 ) OPD    Mineral Dust Optical Depth (wvlen 1): unitless : none
!  (5 ) SD     Mineral Dust Surface Area           : cm2/cm3  : none
!  (6 ) OPSO4  Sulfate Optical Depth    (wavelen 1): unitless : SCALECHEM
!  (7 ) HGSO4  Hygroscopic growth of SO4           : unitless : SCALECHEM
!  (8 ) SSO4   Sulfate Surface Area                : cm2/cm3  : SCALECHEM
!  (9 ) OPBC   Black Carbon Optical Depth (wvlen 1): unitless : SCALECHEM
!  (10) HGBC   Hygroscopic growth of BC            : unitless : SCALECHEM
!  (11) SBC    Black Carbon Surface Area           : cm2/cm3  : SCALECHEM
!  (12) OPOC   Organic C Optical Depth  (wavelen 1): unitless : SCALECHEM
!  (13) HGOC   Hygroscopic growth of OC            : unitless : SCALECHEM
!  (14) SOC    Organic Carbon Surface Area         : cm2/cm3  : SCALECHEM
!  (15) OPSSa  Sea Salt (accum) Opt Depth (wvlen 1): unitless : SCALECHEM
!  (16) HGSSa  Hygroscopic growth of SSa           : unitless : SCALECHEM
!  (17) SSSa   Sea Salt (accum) Surface Area       : cm2/cm3  : SCALECHEM
!  (18) OPSSc  Sea Salt (coarse) Opt Depth(wvlen 1): unitless : SCALECHEM
!  (19) HGSSc  Hygroscopic growth of SSc           : unitless : SCALECHEM
!  (20) SSSc   Sea Salt (coarse) Surface Area      : cm2/cm3  : SCALECHEM  
!  (21) OPD1   Dust Optical Depth Bin 1   (wvlen 1): unitless : SCALECHEM 
!  (22) OPD2   Dust Optical Depth Bin 2   (wvlen 1): unitless : SCALECHEM 
!  (23) OPD3   Dust Optical Depth Bin 3   (wvlen 1): unitless : SCALECHEM 
!  (24) OPD4   Dust Optical Depth Bin 4   (wvlen 1): unitless : SCALECHEM 
!  (25) OPD5   Dust Optical Depth Bin 5   (wvlen 1): unitless : SCALECHEM 
!  (26) OPD6   Dust Optical Depth Bin 6   (wvlen 1): unitless : SCALECHEM 
!  (27) OPD7   Dust Optical Depth Bin 7   (wvlen 1): unitless : SCALECHEM
!  (28) OPSO4  Sulfate Optical Depth    (wavelen 2): unitless : SCALECHEM
!  (29) OPBC   Black Carbon Optical Depth (wvlen 2): unitless : SCALECHEM
!  (30) OPOC   Organic C Optical Depth  (wavelen 2): unitless : SCALECHEM
!  (31) OPSSa  Sea Salt (accum) Opt Depth (wvlen 2): unitless : SCALECHEM
!  (32) OPSSc  Sea Salt (coarse) Opt Depth(wvlen 2): unitless : SCALECHEM
!  (33) OPD1   Dust Optical Depth Bin 1   (wvlen 2): unitless : SCALECHEM 
!  (34) OPD2   Dust Optical Depth Bin 2   (wvlen 2): unitless : SCALECHEM 
!  (35) OPD3   Dust Optical Depth Bin 3   (wvlen 2): unitless : SCALECHEM 
!  (36) OPD4   Dust Optical Depth Bin 4   (wvlen 2): unitless : SCALECHEM 
!  (37) OPD5   Dust Optical Depth Bin 5   (wvlen 2): unitless : SCALECHEM 
!  (38) OPD6   Dust Optical Depth Bin 6   (wvlen 2): unitless : SCALECHEM 
!  (39) OPD7   Dust Optical Depth Bin 7   (wvlen 2): unitless : SCALECHEM
!  (40) OPSO4  Sulfate Optical Depth    (wavelen 3): unitless : SCALECHEM
!  (41) OPBC   Black Carbon Optical Depth (wvlen 3): unitless : SCALECHEM
!  (42) OPOC   Organic C Optical Depth  (wavelen 3): unitless : SCALECHEM
!  (43) OPSSa  Sea Salt (accum) Opt Depth (wvlen 3): unitless : SCALECHEM
!  (44) OPSSc  Sea Salt (coarse) Opt Depth(wvlen 3): unitless : SCALECHEM
!  (45) OPD1   Dust Optical Depth Bin 1   (wvlen 3): unitless : SCALECHEM 
!  (46) OPD2   Dust Optical Depth Bin 2   (wvlen 3): unitless : SCALECHEM 
!  (47) OPD3   Dust Optical Depth Bin 3   (wvlen 3): unitless : SCALECHEM 
!  (48) OPD4   Dust Optical Depth Bin 4   (wvlen 3): unitless : SCALECHEM 
!  (49) OPD5   Dust Optical Depth Bin 5   (wvlen 3): unitless : SCALECHEM 
!  (50) OPD6   Dust Optical Depth Bin 6   (wvlen 3): unitless : SCALECHEM 
!  (51) OPD7   Dust Optical Depth Bin 7   (wvlen 3): unitless : SCALECHEM 
!  (52) ODSLA  Strat. liquid aerosol optical depth : unitless : SCALECHEM
!  (53) SASLA  Strat. liquid aerosol surface area  : cm2/cm3  : SCALECHEM 
!  (54) NDSLA  Strat. liquid aerosol number density: num/cm3  : SCALECHEM 
!  (55) ODSPA  Strat. particulate optical depth    : unitless : SCALECHEM
!  (56) SASPA  Strat. particula. soln. surf. area  : cm2/cm3  : SCALECHEM 
!  (57) NDSPA  Strat. particula. soln. number dens,: num/cm3  : SCALECHEM 
!
!  NOTES:
!  (1 ) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2 ) Don't divide monthly mean AOD by SCALECHEM (rvm, bmy, 12/8/00)
!  (3 ) Use SCALE_A6 for GEOS-2, GEOS-3 fields, since optical depths are read
!        in from disk every 6 hours.  Use SCALECHEM for GEOS-1, GEOS-STRAT
!        fields, since optical depths are computed every chemistry timestep.
!        Use SCALEDYN for CO-OH parameterization simulation. (bmy, 4/23/01)
!  (4 ) Now GEOS-2, GEOS-3 use SCALECHEM for ND21 (bmy, 8/13/01)
!  (5 ) Updated tracers for new aerosols from Mian Chin (rvm, bmy, 3/1/02)
!  (6 ) Now scale online dust fields by SCALECHEM (bmy, 4/9/04)
!  (7 ) Also save out extra diagnostics for cryst sulfur tracers (bmy, 1/5/05)
!  (8 ) Save out extra diagnostics for dust AOD in each size bin (clh, 5/7/10)
!  (9 ) Bug fix, don't write out obsolete tracer #3 for ND21 (bmy, 12/16/11)
!  (10) Added stratospheric aerosols (SDE 04/17/13)
!******************************************************************************
!
      IF ( ND21 > 0 ) THEN
         CATEGORY = 'OD-MAP-$'

         DO M = 1, TMAX(21)
            N  = TINDEX(21,M)
            IF ( N > PD21 ) CYCLE
            NN = N 
               
            ! Select proper unit string (cf list above)
            ! Make sure OPTD and CLDF use SCALE_A3. Rest of fields are
            ! updated every chem timestep and use SCALECHEM. (mps, 3/17/16)
            SELECT CASE( N ) 
               CASE ( 1, 2 )
                  UNIT   = 'unitless'
                  SCALEX = SCALE_A3
               CASE ( 3 )
                  ! This diagnostic is for a GEOS-1/GEOS-STRAT field
                  ! so skip it! (bmy, 12/16/11)
                  CYCLE
               CASE ( 5, 8, 11, 14, 17, 20, 53, 56 )
                  UNIT   = 'cm2/cm3'
                  SCALEX = SCALECHEM
               CASE DEFAULT
                  UNIT   = 'unitless'
                  SCALEX = SCALECHEM
            END SELECT

           ! Is this something that we want to output?
           ! Only output extra wavelengths if they are requested
           IF (Input_Opt%NWVSELECT.eq.3) THEN
               IS_AD21 = .TRUE.
           ELSEIF (Input_Opt%NWVSELECT.eq.2) THEN
               ! Exclude wavelength 3
               IS_AD21 = ((N.le.(5 + 2*NRHAER + 2*(NRHAER+NDUST))).or.
     &                    (N.gt.(5 + 2*NRHAER + 3*(NRHAER+NDUST))))
           ELSEIF (Input_Opt%NWVSELECT.eq.1) THEN
               ! Exclude wavelengths 2 and 3
               IS_AD21 = ((N.le.(5 + 2*NRHAER + 1*(NRHAER+NDUST))).or.
     &                    (N.gt.(5 + 2*NRHAER + 3*(NRHAER+NDUST))))
           ENDIF
               
           ! Dust gets scaled differently
           IS_DUST = (((N.eq.4).or.(N.eq.5)).or.
     &                ((N.gt.(5 + 2*NRHAER + 1*(NRHAER))).and.
     &                 (N.lt.(5 + 2*NRHAER + 1*(NRHAER+NDUST)))).or.
     &                ((N.gt.(5 + 2*NRHAER + 2*(NRHAER))).and.
     &                 (N.le.(5 + 2*NRHAER + 2*(NRHAER+NDUST)))).or.
     &                ((N.gt.(5 + 2*NRHAER + 3*(NRHAER))).and.
     &                 (N.le.(5 + 2*NRHAER + 3*(NRHAER+NDUST)))))

           IF (IS_AD21) THEN
            IF (IS_DUST) THEN

               ! Online or offline dust fields?
               IF ( Input_Opt%LDUST ) THEN

                  ! If LDUST=T, then we are using online dust fields,
                  ! so we must scale by the chemistry timestep. (4/9/04)
                  ARRAY(:,:,1:LD21) = AD21(:,:,1:LD21,N) / SCALEX
      
               ELSE
                  
                  ! If LDUST=F, then we are using offline monthly-mean
                  ! dust fields.  These don't have to be scaled by
                  ! the chemistry timestep. (bmy, 4/9/04)
                  ARRAY(:,:,1:LD21) = AD21(:,:,1:LD21,N)

               ENDIF

            ELSE

               ! For all other types of optical depths, we need 
               ! to scale by the chemistry timestep (bmy, 4/9/04)
               ARRAY(:,:,1:LD21) = AD21(:,:,1:LD21,N) / SCALEX
               
            ENDIF
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD21,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD21) )
           ENDIF
         ENDDO    

      ENDIF
!
!******************************************************************************
!  ND22: J-value diagnostics
!
!   #  : Field : Description                   : Units : Scale factor
!  --------------------------------------------------------------------------
!  (1  ) JNO2  : NO2   J-Value                 : s-1   : SCALE_JV
!  (2  ) JHNO3 : HNO3  J-Value                 : s-1   : SCALE_JV
!  (3  ) JH2O2 : H2O2  J-Value                 : s-1   : SCALE_JV
!  (4  ) JCH2O : CH2O  J-Value                 : s-1   : SCALE_JV
!  (5  ) JO3   : O3    J-Value                 : s-1   : SCALE_JV
!  (6  ) JO3P  : O3P   J-value (Standard,UCX)  : s-1   : SCALE_JV
!        POH   : OH    from O3 (Tropchem,SOA)  : s-1   : SCALE_JV
!  (7  ) JGLYX : GLYX  J-Value                 : s-1   : SCALE
!  (8  ) JMGLY : MGLY  J-Value                 : s-1   : SCALE
!  (9  ) BrO   : BrO   J-Value                 : s-1   : SCALE_JV
!  (10 ) BHOBr : HOBr  J-Value                 : s-1   : SCALE_JV
!  (11 ) BrNO2 : BrNO2 J-Value                 : s-1   : SCALE_JV
!  (12 ) BrNO3 : BrNO3 J-Value                 : s-1   : SCALE_JV
!  (13 ) CHBr3 : CHBr3 J-Value                 : s-1   : SCALE_JV
!  (14 ) Br2   : Br2   J-Value                 : s-1   : SCALE_JV
!  (15 ) JO2   : O2    J-value                 : s-1   : SCALE_JV
!  (16 ) JN2O  : N2O   J-value                 : s-1   : SCALE_JV
!  (17 ) JNO   : NO    J-Value                 : s-1   : SCALE_JV
!  (18 ) JNO3  : NO3   J-Value                 : s-1   : SCALE_JV
!  (19 ) JCFC11: CFC11 J-value                 : s-1   : SCALE_JV
!  (20 ) JCFC12: CFC11 J-value                 : s-1   : SCALE_JV
!  (21 ) JCCl4 : CCl4  J-value                 : s-1   : SCALE_JV
!  (22 ) JCH3Cl: CH3Cl J-value                 : s-1   : SCALE_JV
!  (23 ) JACET : ACET  J-value                 : s-1   : SCALE_JV
!  (24 ) JALD2 : ALD2  J-value                 : s-1   : SCALE_JV
!  (25 ) JMVK  : MVK   J-value                 : s-1   : SCALE_JV
!  (26 ) JMACR : MACR  J-value                 : s-1   : SCALE_JV
!  (27 ) JHAC  : HAC   J-value                 : s-1   : SCALE_JV
!  (28 ) JGLYC : GLYC  J-value                 : s-1   : SCALE_JV
!  (29 ) JPIP  : PIP   J-value                 : s-1   : SCALE_JV
!  (30 ) JIPMN : IPMN  J-value                 : s-1   : SCALE_JV
!  (31 ) JETHLN: ETHLN J-value                 : s-1   : SCALE_JV
!  (32 ) JDHDN : DHDN  J-value                 : s-1   : SCALE_JV
!  (33 ) JHPALD: HPALD J-value                 : s-1   : SCALE_JV
!  (34 ) JISN1 : ISN1  J-value                 : s-1   : SCALE_JV
!  (35 ) JMONITS: MONITS J-value               : s-1   : SCALE_JV
!  (36 ) JMONITU: MONITU J-value               : s-1   : SCALE_JV
!  (37 ) JHONIT: HONIT J-value                 : s-1   : SCALE_JV
!  (71 ) JCH3I : CH3I  J-value (s^-1)          : s-1   : SCALE_JV
!  (81 ) JHCN  : HCN   J-value (s^-1)          : s-1   : SCALE_JV
!
!  NOTES:
!  (1) We must add TRCOFFSET for CH3I and HCN runs, so that GAMAP can
!       recognize those photo rates as distinct from the NO2, HNO3,
!       H2O2, CH2O, O3, and POH photo rates.
!  (2) Pick the right tracer field from AD22 if only a subset of tracers
!       are requested to be printed out. (ccc, 12/15/08)
!  (3) Add GLYX and MGLY tracers (tmf, 3/6/09)
!  (4) Replaced NOx with NO2 (mpayer, 3/14/13)
!  25 Apr 2014 - M. Sulprizio- Now include J-values for ACET
!  27 May 2014 - M. Sulprizio- Now include J-values for ALD2, MVK, MACR, HAC, 
!                              and GLYC
!******************************************************************************
!
      IF ( ND22 > 0 ) THEN
         CATEGORY  = 'JV-MAP-$'
         SCALE_TMP = FLOAT( CTJV ) + 1e-20_fp
         UNIT      = 's-1'

         DO M = 1, TMAX(22)
            N  = TINDEX(22,M)
            IF ( N > PD22 ) CYCLE
            NN = N

            DO L = 1, LD22
               ARRAY(:,:,L) = AD22(:,:,L,NN) / SCALE_TMP(:,:)
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD22,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD22) )

         ENDDO 

      ENDIF     
!
!******************************************************************************
!  ND24: Eastward mass flux from transport (TPCORE, XTP)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  MASSFLEW : Eastward mass flux - transport  : kg/s     : SCALEDYN
!
!  NOTES:
!  (1) MASSFLEW is REAL(fp)...store to ARRAY, which is REAL*4
!       before sending to BPCH or IJSCAL (bey, bmy, 4/23/99)
!  (2) Now only write LD24 levels out to the bpch file (bmy, 12/7/00)
!******************************************************************************
!
      IF ( ND24 > 0 ) THEN
         CATEGORY = 'EW-FLX-$'
         UNIT = 'kg/s'

         DO M = 1, TMAX(24)
            N  = TINDEX(24,M)
            IF ( N > nAdvect ) CYCLE
            NN = N
            
            ARRAY(:,:,1:LD24) = MASSFLEW(:,:,LLPAR:LLPAR-LD24+1:-1,N) /
     &                          SCALEDYN

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD24,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD24) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND25: Northward mass flux from transport (TPCORE, YTP)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  MASSFLNS : Northward mass flux - transport : kg/s     : SCALEDYN
!
!  NOTES:
!  (1) MASSFLNS is REAL(fp)...store to ARRAY, which is REAL*4
!       before sending to BPCH or IJSCAL (bey, bmy, 4/23/99)
!  (2) Now only write LD25 levels out to the bpch file (bmy, 12/7/00)
!******************************************************************************
!  
      IF ( ND25 > 0 ) THEN
         CATEGORY = 'NS-FLX-$'
         UNIT = 'kg/s'

         DO M = 1, TMAX(25)
            N  = TINDEX(25,M)
            IF ( N > nAdvect ) CYCLE
            NN = N

            ARRAY(:,:,1:LD25) = MASSFLNS(:,:,LLPAR:LLPAR-LD25+1:-1,N) /
     &                          SCALEDYN

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD25,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD25) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND26: Upward mass flux from transport (TPCORE, FZPPM)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  MASSFLUP : Upward mass flux - transport    : kg/s     : SCALEDYN
!
!  NOTES:
!  (1) MASSFLNS is REAL(fp)...store to ARRAY, which is REAL*4
!       before sending to BPCH or IJSCAL (bey, bmy, 4/23/99)
!  (2) Now only write LD26 levels to the bpch file (bmy, 12/7/00)
!******************************************************************************
!  
      IF ( ND26 > 0 ) THEN
         CATEGORY = 'UP-FLX-$'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(26)
            N  = TINDEX(26,M)
            IF ( N > nAdvect ) CYCLE
            NN = N
            
            ARRAY(:,:,1:LD26) = MASSFLUP(:,:,LLPAR:LLPAR-LD26+1:-1,N) /
     &                          SCALEDYN
               
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES, 
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD26,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD26) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND28: Biomass burning diagnostic 
!
!   #  : Field : Description    : Units            : Scale factor
!  --------------------------------------------------------------------------
!  (1  ) NOx   : NOx            : molec NOx /cm2/s : SCALESRCE
!  (3  ) PAN   : PAN            : molec PAN /cm2/s : SCALESRCE
!  (4  ) CO    : CO             : molec CO  /cm2/s : SCALESRCE
!  (5  ) ALK4  : >=C4 alkanes   : atoms C   /cm2/s : SCALESRCE
!  (7  ) HNO3  : HNO3           : molec HNO3/cm2/s : SCALESRCE
!  (9  ) ACET  : Acetone        : atoms C   /cm2/s : SCALESRCE
!  (10 ) MEK   : Ketones(>C3)   : atoms C   /cm2/s : SCALESRCE
!  (11 ) ALD2  : Acetaldehyde   : atoms C   /cm2/s : SCALESRCE
!  (18 ) PRPE  : Propene        : atoms C   /cm2/s : SCALESRCE
!  (19 ) C3H8  : Propane        : atoms C   /cm2/s : SCALESRCE
!  (20 ) CH2O  : Formaldehyde   : molec CH2O/cm2/s : SCALESRCE
!  (21 ) C2H6  : Ethane         : atoms C   /cm2/s : SCALESRCE
!  (26 ) SO2   : Sulfur dioxide : molec SO2 /cm2/s : SCALESRCE
!  (30 ) NH3   : Ammonia        : molec NH3 /cm2/s : SCALESRCE
!  (34 ) BC    : Black carbon   : atoms C   /cm2/s : SCALESRCE
!  (35 ) OC    : Organic carbon : atoms C   /cm2/s : SCALESRCE
!  (71 ) CH4   : Methane        : molec CH4 /cm2/s : SCALESRCE
!  (99 ) MTPA  : Lumped monoterp: atoms C   /cm2/s : SCALESRCE
!  (116) BENZ  : Benzene        : atoms C   /cm2/s : SCALESRCE
!  (117) TOLU  : Toluene        : atoms C   /cm2/s : SCALESRCE
!  (118) XYLE  : Xylene         : atoms C   /cm2/s : SCALESRCE
!  (126) EOH   : Ethanol        : atoms C   /cm2/s : SCALESRCE
!  (127) MGLY  : Methylglyoxyal : atoms C   /cm2/s : SCALESRCE
!
!  NOTES:
!  (1) Use the F90 intrinsic "ANY" function to make sure that N 
!       corresponds to actual biomass burning tracers (bmy, 4/8/99)
!  (2) ND28 now uses allocatable array AD28 instead of AIJ. (bmy, 3/16/00)
!  (3) Now write biofuel burning tracers to the punch file in the same order 
!       as they are listed in "diag.dat". (bmy, 4/17/01)
!  (4) Biomass diagnostics are now tracked by HEMCO (bmy, 8/14/14)
!  (5) NOTE: Tagged CO biomass tracers are now tracked in ND29 (bmy, 3/17/16)
!******************************************************************************
!
      IF ( ND28 > 0 ) THEN

         ! Category name for the BPCH file
         CATEGORY = 'BIOBSRCE'

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         !%%% NOTE: FACTOR is the multiplication factor that converts
         !%%% kg/m2/s (as tracked by HEMCO) to either molec/cm2/s.        
         !%%% (bmy, 8/14/14)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         !%%%%% NO %%%%%
         IF ( id_NO > 0 ) THEN
            DiagnName = 'BIOMASS_NO'
            UNIT      = 'molec/cm2/s'
            N         = id_NO
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp * CM2PERM2 )
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% PAN %%%%%
         IF ( id_PAN > 0 ) THEN
            DiagnName = 'BIOMASS_PAN'
            UNIT      = 'molec/cm2/s'
            N         = id_PAN
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp * CM2PERM2 )
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT, N, 1, -1, .TRUE., FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% CO %%%%%
         IF ( id_CO > 0 ) THEN
            DiagnName = 'BIOMASS_CO'
            UNIT      = 'molec/cm2/s'
            N         = id_CO
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% ALK4 %%%%%
         IF ( id_ALK4 > 0 ) THEN
            DiagnName = 'BIOMASS_ALK4'
            UNIT      = 'atoms C/cm2/s'
            N         = id_ALK4
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% HNO3 %%%%%
         IF ( id_HNO3 > 0 ) THEN
            DiagnName = 'BIOMASS_HNO3'
            UNIT      = 'molec/cm2/s'
            N         = id_HNO3
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT, N, 1, -1, .TRUE., FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% ACET %%%%%
         IF ( id_ACET > 0 ) THEN
            DiagnName = 'BIOMASS_ACET'
            UNIT      = 'atoms C/cm2/s'
            N         = id_ACET
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% MEK %%%%%
         IF ( id_MEK > 0 ) THEN
            DiagnName = 'BIOMASS_MEK'
            UNIT      = 'atoms C/cm2/s'
            N         = id_MEK
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% ALD2 %%%%%
         IF ( id_ALD2 > 0 ) THEN
            DiagnName = 'BIOMASS_ALD2'
            UNIT      = 'atoms C/cm2/s'
            N         = id_ALD2
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% PRPE %%%%%
         IF ( id_PRPE > 0 ) THEN
            DiagnName = 'BIOMASS_PRPE'
            UNIT      = 'atoms C/cm2/s'
            N         = id_PRPE
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% C3H8 %%%%%
         IF ( id_C3H8 > 0 ) THEN
            DiagnName = 'BIOMASS_C3H8'
            UNIT      = 'atoms C/cm2/s'
            N         = id_C3H8
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% CH2O %%%%%
         IF ( id_CH2O > 0 ) THEN
            DiagnName = 'BIOMASS_CH2O'
            UNIT      = 'molec/cm2/s'
            N         = id_CH2O
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% C2H6 %%%%%
         IF ( id_C2H6 > 0 ) THEN
            DiagnName = 'BIOMASS_C2H6'
            UNIT      = 'atoms C/cm2/s'
            N         = id_C2H6
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% SO2 %%%%%
         IF ( id_SO2 > 0 ) THEN
            DiagnName = 'BIOMASS_SO2'
            UNIT      = 'molec/cm2/s'
            N         = id_SO2
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% NH3 %%%%%
         IF ( id_NH3 > 0 ) THEN
            DiagnName = 'BIOMASS_NH3'
            UNIT      = 'molec/cm2/s'
            N         = id_NH3
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% NAP %%%%% ; added krt
         IF ( id_NAP > 0 ) THEN
            DiagnName = 'BIOMASS_NAP'
            UNIT      = 'atoms C/cm2/s'
            N         = id_NAP
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         ! Free pointer
         SpcInfo => NULL()

         !%%%%% BC %%%%%

         ! Save all as BCPI
         UNIT      = 'atoms C/cm2/s'
         N         = id_BCPI
         SpcInfo => State_Chm%SpcData(N)%Info
         FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2

         ! Get BIOMASS_BCPI 
         DiagnName = 'BIOMASS_BCPI'
         CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                   FLAG, RC,   cName=TRIM(DiagnName), 
     &                   AutoFill=1, 
     &                   COL=HcoState%Diagn%HcoDiagnIDManual )
         IF ( RC /= HCO_SUCCESS ) THEN
            MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
            CALL ERROR_STOP ( MSG, LOC ) 
         ENDIF
      
         ! Save into ARRAY. Convert to units
         IF ( FLAG == HCO_SUCCESS ) THEN
            ARRAY(:,:,1) = DiagnCnt%Arr2D%Val(:,:) 
     &                   * FACTOR 
         ELSE
            MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
            MSG = TRIM(MSG) // ' - will write zeros!'
            CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
         ENDIF
      
         ! Get BIOMASS_BCPO 
         DiagnName = 'BIOMASS_BCPO'
         DiagnCnt  => NULL()
         CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                   FLAG, RC,   cName=TRIM(DiagnName), 
     &                   AutoFill=1, 
     &                   COL=HcoState%Diagn%HcoDiagnIDManual )
         IF ( RC /= HCO_SUCCESS ) THEN
            MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
            CALL ERROR_STOP ( MSG, LOC ) 
         ENDIF
      
         ! Save into ARRAY. Convert to units
         IF ( FLAG == HCO_SUCCESS ) THEN
            ARRAY(:,:,1) = ARRAY(:,:,1) 
     &                   + ( DiagnCnt%Arr2D%Val(:,:) * FACTOR )
         ELSE
            MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
            MSG = TRIM(MSG) // ' - will write zeros!'
            CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
         ENDIF

         ! Write combined diagnostics to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT ,     DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )
      
         ! Free pointer
         DiagnCnt => NULL()
         SpcInfo => NULL()

         !%%%%% POA simulation %%%%% ; updated krt
         IF ( id_POG1 > 0 .and. id_POG2 > 0 ) THEN
            DiagnName = 'BIOMASS_POG1'
            UNIT      = 'atoms C/cm2/s'
            N         = id_POG1
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )

            DiagnName = 'BIOMASS_POG2'
            UNIT      = 'atoms C/cm2/s'
            N         = id_POG2
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )

         !%%%%% OC %%%%%
         ELSE

            ! Save all as OCPI
            UNIT      = 'atoms C/cm2/s'
            N         = id_OCPI
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
   
            ! Get BIOMASS_OCPI 
            DiagnName = 'BIOMASS_OCPI'
            CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                      FLAG, RC,   cName=TRIM(DiagnName), 
     &                      AutoFill=1, 
     &                      COL=HcoState%Diagn%HcoDiagnIDManual )
            IF ( RC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
               CALL ERROR_STOP ( MSG, LOC ) 
            ENDIF
         
            ! Save into ARRAY. Convert to units
            IF ( FLAG == HCO_SUCCESS ) THEN
               ARRAY(:,:,1) = DiagnCnt%Arr2D%Val(:,:) 
     &                      * FACTOR
            ELSE
               MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
               MSG = TRIM(MSG) // ' - will write zeros!'
               CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
            ENDIF
         
            ! Get BIOMASS_OCPO 
            DiagnName = 'BIOMASS_OCPO'
            DiagnCnt  => NULL()
            CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                         FLAG, RC, cName=TRIM(DiagnName), 
     &                         AutoFill=1, 
     &                         COL=HcoState%Diagn%HcoDiagnIDManual )
            IF ( RC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
               CALL ERROR_STOP ( MSG, LOC ) 
            ENDIF
         
            ! Save into ARRAY. Convert to units
            IF ( FLAG == HCO_SUCCESS ) THEN
               ARRAY(:,:,1) = ARRAY(:,:,1) 
     &                      + ( DiagnCnt%Arr2D%Val(:,:) * FACTOR )
            ELSE
               MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
               MSG = TRIM(MSG) // ' - will write zeros!'
               CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
            ENDIF
   
            ! Write combined diagnostics to disk
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT ,     DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         
            ! Free pointer
            DiagnCnt => NULL()
         ENDIF

         !%%%%% CH4 %%%%%
         IF ( id_CH4 > 0 ) THEN
            DiagnName = 'BIOMASS_CH4'
            UNIT      = 'molec/cm2/s'
            N         = id_CH4
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% MTPA %%%%%
         IF ( id_MTPA > 0 ) THEN
            DiagnName = 'BIOMASS_MTPA'
            UNIT      = 'atoms C/cm2/s'
            N         = id_MTPA
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT, N, 1, -1, .TRUE., FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% BENZ %%%%%
         IF ( id_BENZ > 0 ) THEN
            DiagnName = 'BIOMASS_BENZ'
            UNIT      = 'atoms C/cm2/s'
            N         = id_BENZ
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT, N, 1, -1, .TRUE., FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% TOLU %%%%%
         IF ( id_TOLU > 0 ) THEN
            DiagnName = 'BIOMASS_TOLU'
            UNIT      = 'atoms C/cm2/s'
            N         = id_TOLU
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT, N, 1, -1, .TRUE., FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% XYLE %%%%%
         IF ( id_XYLE > 0 ) THEN
            DiagnName = 'BIOMASS_XYLE'
            UNIT      = 'atoms C/cm2/s'
            N         = id_XYLE
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT, N, 1, -1, .TRUE., FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% EOH %%%%%
         IF ( id_EOH > 0 ) THEN
            DiagnName = 'BIOMASS_EOH'
            UNIT      = 'atoms C/cm2/s'
            N         = id_EOH
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT, N, 1, -1, .TRUE., FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% MGLY %%%%%
         IF ( id_MGLY > 0 ) THEN
            DiagnName = 'BIOMASS_MGLY'
            UNIT      = 'atoms C/cm2/s'
            N         = id_MGLY
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT, N, 1, -1, .TRUE., FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         ! Free pointer
         SpcInfo => NULL()

      ENDIF
!
!******************************************************************************
!  ND29: CO source diagnostics
!
!   #  Field  : Description             : Units       : Scale factor
!  --------------------------------------------------------------------------
!  (1) COanth : CO from Anthro Sources  : molec/cm2/s : SCALESRCE
!  (2) CObb   : CO from Biomass Burning : molec/cm2/s : SCALESRCE
!  (3) CObf   : CO from Biofuel Burning : molec/cm2/s : SCALESRCE
!  (4) COmeth : CO from Methanol        : molec/cm2/s : SCALESRCE
!  (5) COmono : CO from Monoterpenes    : molec/cm2/s : SCALESRCE
!
!  NOTES:
!  (1) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2) ND29 now uses allocatable array AD29 instead of AIJ. (bmy, 3/16/00)
!  (3) Added CO-sources from isoprene and monoterpenes (bnd, bmy, 1/2/01)
!  (4) Added tagged CO emissions diagnostics from HEMCO (bmy, 3/17/16)
!******************************************************************************
!
      IF ( ND29 > 0 ) THEN

         ! Diagnostic category and units
         CATEGORY = 'CO--SRCE'
         UNIT     = 'molec/cm2/s'

         SpcInfo => State_Chm%SpcData(id_CO)%Info

         ! Converts from kg/m2/s (in HEMCO) to molec/cm2/s (for bpch)
         FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2

         ! Loop over tagged tracers if necessary
         IF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN
            NN = 24
         ELSE
            NN = 6
         ENDIF

         ! Loop over all slots of ND29 
         ! 1-6  are for the fullchem simulation
         ! 7-24 are for the tagged CO simulation
         ! Also note: Tagged CO lumps biofuels with anthro emissions
         DO N = 1, NN

            ! Initalize number of levels for diagnostic
            Levs = -1

            ! Find the proper HEMCO diagnostic name 
            SELECT CASE ( N )
               CASE( 1 )
                  DiagnName = 'ANTHROPOGENIC_CO'
                  Levs      = LLPAR
               CASE( 2 )
                  DiagnName = 'BIOMASS_CO'
               CASE( 3 )
                  DiagnName = 'BIOFUEL_CO'
               CASE( 4 ) 
                  CYCLE        ! CO from methanol doesn't seem to be defined!?
               CASE ( 5 )
                  DiagnName = 'BIOGENIC_CO'
               CASE( 6 ) 
                  DiagnName = 'SHIP_CO'

               ! Tagged CO: Anthro sectors, USA
               CASE( 7 )
                  DiagnName = 'ANTHRO_BIOFUEL_TAGCO_US'
                  Levs      = LLPAR
               CASE( 8 ) 
                  DiagnName = 'AIRCRAFT_TAGCO_US'
                  Levs      = LLPAR
               CASE( 9 )
                  DiagnName = 'SHIP_TAGCO_US'

               ! Tagged CO: Anthro sectors, Europe
               CASE( 10 )
                  DiagnName = 'ANTHRO_BIOFUEL_TAGCO_EUR'
                  Levs      = LLPAR
               CASE( 11 )
                  DiagnName = 'AIRCRAFT_TAGCO_EUR'
                  Levs      = LLPAR
               CASE( 12 )
                  DiagnName = 'SHIP_TAGCO_EUR'

               ! Tagged CO: Anthro sectors, Asia
               CASE( 13 )
                  DiagnName = 'ANTHRO_BIOFUEL_TAGCO_ASIA'
                  Levs      = LLPAR
               CASE( 14 )
                  DiagnName = 'AIRCRAFT_TAGCO_ASIA'
                  Levs      = LLPAR
               CASE( 15 ) 
                  DiagnName = 'SHIP_TAGCO_ASIA'

               ! Tagged CO: Anthro sectors, rest of world
               CASE( 16 )
                  DiagnName = 'ANTHRO_BIOFUEL_TAGCO_OTHER'
                  Levs      = LLPAR
               CASE( 17 )
                  DiagnName = 'AIRCRAFT_TAGCO_OTHER'
                  Levs      = LLPAR
               CASE( 18 )
                  DiagnName = 'SHIP_TAGCO_OTHER'

               ! Tagged CO: Regional biomass tracers
               CASE( 19 )
                  DiagnName = 'BIOMASS_TAGCO_USA'
               CASE( 20 )
                  DiagnName = 'BIOMASS_TAGCO_AFRICA'
               CASE( 21 )
                  DiagnName = 'BIOMASS_TAGCO_ASIA'
               CASE( 22 )
                  DiagnName = 'BIOMASS_TAGCO_OCEANIA'
               CASE( 23 )
                  DiagnName = 'BIOMASS_TAGCO_EUROPE'
               CASE( 24 )
                  DiagnName = 'BIOMASS_TAGCO_OTHER'
                  
               ! Default
               CASE DEFAULT
                  DiagnName = 'NOTHINGTOFIND'
            END SELECT

            ! Write data from HEMCO to the bpch file
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, Levs, .TRUE., FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDDO

         ! Free pointer
         SpcInfo => NULL()
      ENDIF
!
!******************************************************************************
!  ND30: Land map diagnostic
!
!   #  Field : Description             : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) LWI   : GMAO Land-Water indices  : unitless  : SCALED 
!
!  NOTES: 
!  (1) Values are: 0=water; 1=land; 2=ice (bmy, 8/18/05)
!******************************************************************************
!
      IF ( ND30 > 0 ) THEN
         CATEGORY = 'LANDMAP'
         UNIT     = 'unitless'
            
         ARRAY(:,:,1) = AD30(:,:) / SCALEDIAG
         NN           = 1 

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, NN,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )
      ENDIF
!
!******************************************************************************
!  ND31: Surface pressure diagnostic
!
!   #  Field : Description                      : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) Pedge : Pressure at bot edge of level L  : mb        : SCALEDYN
!
!  NOTES: 
!  (1) The ASCII punch file was using SCALE2 instead of SCALE1.
!       This has now been fixed. (hyl, bmy, 12/21/99).
!  (2) Now use AD31 dynamically allocatable array (bmy, 2/17/00)
!  (3) Bug fix: write out 1 level to the bpch file (bmy, 12/7/00)
!  (4) Now remove SCALE1, replace with SCALEDYN (bmy, 2/24/03)
!  (5) Now save out true pressure at level edges.  Now   (bmy, 5/8/07)
!******************************************************************************
!   
      IF ( ND31 > 0 ) THEN
         CATEGORY          = 'PEDGE-$'
         UNIT              = 'mb'

         ARRAY(:,:,1:LD31) = AD31(:,:,1:LD31) / SCALEDIAG
         NN                = 1

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, NN,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD31,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD31) )
      ENDIF

!
!******************************************************************************
!  ND32: NOx source diagnostic
!
!  Levels        : Field                  : Units       : Scale Factor
!  -------------------------------------------------------------------------
!  1 - LLCHEM    : Aircraft NOx           : molec/cm2/s : SCALESRCE
!  1 - ?         : Anthropogenic NOx      : molec/cm2/s : SCALESRCE
!  Surface       : Biomass Burning NOx    : molec/cm2/s : SCALESRCE
!  Surface       : Biofuel Burning NOx    : molec/cm2/s : SCALESRCE
!  Surface       : Fertilizer NOx         : molec/cm2/s : SCALESRCE
!  1 - LLCONVM   : Lightning NOx          : molec/cm2/s : SCALESRCE
!  Surface       : Soil NOx               : molec/cm2/s : SCALESRCE
!  Above TP      : NOx from upper boundary: molec/cm2/s : SCALEDYN
!
!  Print out all of the types of NOx, for all levels.
!
!  NOTES:
!  (1) Only print out ND32 if for an O3 chemistry run ( NSRCX == 3 ),
!       and if NOx is a defined tracer ( id_NOX > 0 ). (bmy, 5/26/99)
!  (2) ND32 now uses allocatable arrays instead of AIJ. (bmy 3/16/00)
!  (3) Added biofuel burning to ND32 diagnostic (bmy, 9/12/00)
!  (4) Replaced NOx emissions with NO emissions (mpayer, 3/14/13)
!  (5) Aircraft NOx emissions now extend to LLPAR for AEIC inventory
!       (mpayer, 7/31/13)
!  (6) ND32 is now tracked by HEMCO (bmy, 8/15/14)
!******************************************************************************
!
      IF ( ND32  > 0                     .and. 
     &     id_NO > 0                     .and. 
     &     Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

         ! All categories of NOx are in molec/cm2/s
         UNIT   = 'molec/cm2/s'

         ! Diagnostic tracer number
         N      = id_NO
         SpcInfo => State_Chm%SpcData(N)%Info

         ! Converts from kg/m2/s (in HEMCO) to molec/cm2/s (for bpch)
         FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2

         !%%%%% Aircraft NO %%%%%
         DiagnName = 'AIRCRAFT_NO' 
         CATEGORY  = 'NO-AC-$'
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, LLPAR, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%% Anthropogenic NO %%%%%
         CATEGORY  = 'NO-AN-$'

         !----- (1) Get ANTHROPOGENIC NO from HEMCO data structure ---
         DiagnName = 'ANTHROPOGENIC_NO' 
         DiagnCnt  => NULL()

         ! Zero the diagnostic array
         ARRAY     = 0e0

         ! Get diagnostics from HEMCO
         CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                   FLAG, RC,   cName=TRIM(DiagnName), 
     &                   AutoFill=1, 
     &                   COL=HcoState%Diagn%HcoDiagnIDManual )
         IF ( RC /= HCO_SUCCESS ) THEN
            MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
            CALL ERROR_STOP ( MSG, LOC ) 
         ENDIF

         ! Save into ARRAY. Convert to units
         IF ( FLAG == HCO_SUCCESS ) THEN
            DO L = 1, LLPAR
               ARRAY(:,:,L) = DiagnCnt%Arr3D%Val(:,:,L) * FACTOR 
            ENDDO
         ELSE
            MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
            MSG = TRIM(MSG) // ' - will write zeros!'
            CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
         ENDIF

         !----- (2) Then get SHIP NO and add it to ANTHROPOGENIC NO -----
         DiagnName = 'SHIP_NO' 
         DiagnCnt  => NULL()

         ! Get diagnostics from HEMCO
         CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                   FLAG, RC,   cName=TRIM(DiagnName), 
     &                   AutoFill=1, 
     &                   COL=HcoState%Diagn%HcoDiagnIDManual )
         IF ( RC /= HCO_SUCCESS ) THEN
            MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
            CALL ERROR_STOP ( MSG, LOC ) 
         ENDIF

         ! Save into ARRAY. Convert to units
         IF ( FLAG == HCO_SUCCESS ) THEN
            ARRAY(:,:,1) = ARRAY(:,:,1) +
     &                   ( DiagnCnt%Arr2D%Val(:,:) * FACTOR )
         ELSE
            MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
            MSG = TRIM(MSG) // ' - will write zeros!'
            CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
         ENDIF

         ! Write combined ANTHROPOGENIC NO + SHIP NO to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT ,     DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LLPAR,    IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LLPAR) )

         ! Free pointer
         DiagnCnt => NULL()

         !%%%%% Biomass NO %%%%%
         DiagnName = 'BIOMASS_NO' 
         CATEGORY  = 'NO-BIOB'
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%% Biofuel NO %%%%%
         DiagnName = 'BIOFUEL_NO' 
         CATEGORY  = 'NO-BIOF'
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%% Fertilizer NO %%%%%
         DiagnName = 'FERTILIZER_NO' 
         CATEGORY  = 'NO-FERT'
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 0, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%% Lightning NO %%%%%
         DiagnName = 'LIGHTNING_NO'
         CATEGORY  = 'NO-LI-$'
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, LLPAR, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
               
         !%%%%% Soil NO %%%%%
         DiagnName = 'SOIL_NO' 
         CATEGORY  = 'NO-SOIL'
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         ! Free pointer
         SpcInfo => NULL()

      ENDIF
!
!******************************************************************************
!  ND33: Atmospheric column sum of Tracer
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) COLUMN-T : Trop. Column Sum of Tracer  : kg        : SCALEDYN
!
!  NOTES: 
!  (1) Now use dynamically allocatable array AD33 (bmy, 2/17/00)
!  (2) Rename category to COLUMN-T, since this is a column sum of tracer over
!       the entire atmosphere, not just the troposphere. (bmy, 4/3/02)
!  (3) Now replace SCALE1 with SCALEDYN (bmy, 3/27/03)
!******************************************************************************
!
      IF ( ND33 > 0 ) THEN
         CATEGORY = 'COLUMN-T'
         UNIT     = 'kg'

         DO M = 1, TMAX(33)
            N  = TINDEX(33,M)
            IF ( N > nAdvect ) CYCLE
            NN = N
            
            ARRAY(:,:,1) = AD33(:,:,N) / SCALEDIAG

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF  
!
!******************************************************************************
!  ND34: Biofuel burning diagnostic 
!
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  HEMCO LUMPS MOST BIOFUELS IN WITH THE ANTHROPOGENICS.  THUS, MANY  %%% 
!  %%%  OF THE BIOFUEL DIAGNOSTICS ARE NOW ZERO.  IT'S BETTER TO COMPARE   %%%
!  %%%  SUM OF ANTHROPOGENIC + BIOFUEL EMISSIONS. (bmy, 8/15/14)           %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   #  : Field : Description         : Units            : Scale factor
!  --------------------------------------------------------------------------
!  (1  ) NOx   : NOx                 : molec NOx /cm2/s : SCALESRCE
!  (4  ) CO    : CO                  : molec CO  /cm2/s : SCALESRCE
!  (5  ) ALK4  : Alkanes(>C4)        : atoms C   /cm2/s : SCALESRCE
!  (9  ) ACET  : Acetone             : atoms C   /cm2/s : SCALESRCE
!  (10 ) MEK   : Metyl Ethyl Ketone  : atoms C   /cm2/s : SCALESRCE
!  (11 ) ALD2  : Acetaldehyde        : atoms C   /cm2/s : SCALESRCE
!  (18 ) PRPE  : Alkenes(>=C3)       : atoms C   /cm2/s : SCALESRCE
!  (19 ) C3H8  : Propane             : atoms C   /cm2/s : SCALESRCE
!  (20 ) CH2O  : Formaldehyde        : molec CH2O/cm2/s : SCALESRCE
!  (26 ) SO2   : SO2                 : molec SO2 /cm2/s : SCALESRCE
!  (30 ) NH3   : NH3                 : molec NH3 /cm2/s : SCALESRCE
!  (116) BENZ  : Benzene             : atoms C   /cm2/s : SCALESRCE
!  (117) TOLU  : Toluene             : atoms C   /cm2/s : SCALESRCE
!  (118) XYLE  : Xylene              : atoms C   /cm2/s : SCALESRCE
!  (126) EOH   : Ethanol             : atoms C   /cm2/s : SCALESRCE
!
!  NOTES:
!  (1) Use the F90 intrinsic "ANY" function to make sure that N 
!       corresponds to actual biofuel burning tracers (bmy, 3/15/01)
!  (3) Now write biofuel burning tracers to the punch file in the same order 
!       as they are listed in "diag.dat". (bmy, 4/17/01)
!  (4) Use BFTRACE and NBFTRACE to get the right index for AD34. 
!      (ccc, 12/8/2008)
!  (5) HEMCO now tracks biofuel diagnostics.  RETRO biofuels are lumped into
!       anthropogenic emissions, so several of these entries will be zero.
!******************************************************************************
!
      IF ( ND34 > 0 ) THEN

         ! Diagnostic category for BPCH output
         CATEGORY = 'BIOFSRCE'

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         !%%% NOTE: FACTOR is the multiplication factor that converts
         !%%% kg/m2/s (as tracked by HEMCO) to either molec/cm2/s. 
         !%%% (bmy, 8/14/14)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         !%%%%% NO %%%%%
         IF ( id_NO > 0 ) THEN
            DiagnName = 'BIOFUEL_NO'
            UNIT     = 'molec/cm2/s'
            N         = id_NO
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% CO %%%%%
         IF ( id_CO > 0 ) THEN
            DiagnName = 'BIOFUEL_CO'
            UNIT      = 'molec/cm2/s'
            N         = id_CO
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% ALK4 %%%%%
         IF ( id_ALK4 > 0 ) THEN
            DiagnName = 'BIOFUEL_ALK4'
            UNIT      = 'atoms C/cm2/s'
            N         = id_ALK4
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% ACET %%%%%
         IF ( id_ACET > 0 ) THEN
            DiagnName = 'BIOFUEL_ACET'
            UNIT      = 'atoms C/cm2/s'
            N         = id_ACET
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% MEK %%%%%
         IF ( id_MEK > 0 ) THEN
            DiagnName = 'BIOFUEL_MEK'
            UNIT      = 'atoms C/cm2/s'
            N         = id_MEK
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% ALD2 %%%%%
         IF ( id_ALD2 > 0 ) THEN 
            DiagnName = 'BIOFUEL_ALD2'
            UNIT      = 'atoms C/cm2/s'
            N         = id_ALD2
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% PRPE %%%%%
         IF ( id_PRPE > 0 ) THEN
            DiagnName = 'BIOFUEL_PRPE'
            UNIT      = 'atoms C/cm2/s'
            N         = id_PRPE
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% C3H8 %%%%%
         IF ( id_C3H8 > 0 ) THEN
            DiagnName = 'BIOFUEL_C3H8'
            UNIT      = 'atoms C/cm2/s'
            N         = id_C3H8
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% CH2O %%%%%
         IF ( id_CH2O > 0 ) THEN
            DiagnName = 'BIOFUEL_CH2O'
            UNIT      = 'molec/cm2/s'
            N         = id_CH2O
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% C2H6 %%%%%
         IF ( id_C2H6 > 0 ) THEN
            DiagnName = 'BIOFUEL_C2H6'
            UNIT      = 'atoms C/cm2/s'
            N         = id_C2H6
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% SO2 %%%%%
         IF ( id_SO2 > 0 ) THEN
            DiagnName = 'BIOFUEL_SO2'
            UNIT      = 'molec/cm2/s'
            N         = id_SO2
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% NH3 %%%%%
         IF ( id_NH3 > 0 ) THEN
            DiagnName = 'BIOFUEL_NH3'
            UNIT      = 'molec/cm2/s'
            N         = id_NH3
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR,RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% BENZ %%%%%
         IF ( id_BENZ > 0 ) THEN
            DiagnName = 'BIOFUEL_BENZ'
            UNIT      = 'atoms C/cm2/s'
            N         = id_BENZ
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR,RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% TOLU %%%%%
         IF ( id_TOLU > 0 ) THEN
            DiagnName = 'BIOFUEL_TOLU'
            UNIT      = 'atoms C/cm2/s'
            N         = id_TOLU
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR,RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% XYLE %%%%%
         IF ( id_XYLE > 0 ) THEN
            DiagnName = 'BIOFUEL_XYLE'
            UNIT      = 'atoms C/cm2/s'
            N         = id_XYLE
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR,RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !%%%%% EOH %%%%%
         IF ( id_EOH > 0 ) THEN
            DiagnName = 'BIOFUEL_EOH'
            UNIT      = 'atoms C/cm2/s'
            N         = id_EOH
            SpcInfo => State_Chm%SpcData(N)%Info
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR,RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         ! Free pointer
         SpcInfo => NULL()
      ENDIF
!
!******************************************************************************
!  ND35: Tracer concentration at 500 mb 
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) 500-AVRG : Tracer at 500 mb            : v/v       : SCALEDYN
!
!  NOTES:
!  (1) Now use dynamically allocatable array AD35 (bmy, 2/17/00)
!  (2) Now replace SCALE1 with SCALEDYN (bmy, 2/24/03)
!******************************************************************************
!
      IF ( ND35 > 0 ) THEN
         CATEGORY = '500-AVRG'        
         UNIT     = ''

         DO M = 1, TMAX(35)
            N  = TINDEX(35,M)
            IF ( N > nAdvect ) CYCLE
            NN = N
               
            ARRAY(:,:,1) = AD35(:,:,N) / SCALEDIAG

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF  
!
!******************************************************************************
!  ND36: Anthropogenic source diagnostic
!
!   #    Field  : Description                     : Units         : S. Factor
!  ---------------------------------------------------------------------------
!  (1  ) NOx    : NOx                             : molec/cm2/s   : SCALE3
!  (2  ) O3     : O3                              : molec/cm2/s   : SCALE3
!  (4  ) CO     : CO                              : molec/cm2/s   : SCALE3
!  (5  ) ALK4   : Alkanes(>C4)                    : atoms C/cm2/s : SCALE3
!  (9  ) ACET   : Acetone                         : atoms C/cm2/s : SCALE3
!  (10 ) MEK    : Ketones(>C3)                    : atoms C/cm2/s : SCALE3
!  (11 ) ALD2   : Acetaldehyde                    : atoms C/cm2/s : SCALE3
!  (18 ) PRPE   : Propene                         : atoms C/cm2/s : SCALE3 
!  (19 ) C3H8   : Propane                         : atoms C/cm2/s : SCALE3 
!  (20 ) CH2O   : Formaldehyde                    : atoms C/cm2/s : SCALE3
!  (21 ) C2H6   : Ethane                          : atoms C/cm2/s : SCALE3
!  (71 ) CH4    : Methane                         : molec/cm2/s   : SCALE3
!  (116) BENZ   : Benzene                         : atoms C/cm2/s : SCALE3
!  (117) TOLU   : Toluene                         : atoms C/cm2/s : SCALE3
!  (118) XYLE   : Xylene                          : atoms C/cm2/s : SCALE3
!  (126) EOH    : Ethanol                         : atoms C/cm2/s : SCALE3
!
!  NOTES:
!  (1) ND36 is also used for CH3I emissions diagnostics when NSRCX=2.
!  (2) For an O3 run (NSRCX = 3, the "default" run) make sure that the 
!       tracer number N matches an entry in the IDEMS emission index 
!       array (bmy, 4/9/99)  
!  (3) Write the tracers out to the punch file in the same order as
!       they are listed in the IDEMS array.  Thus, we have to re-assign
!       N = IDEMS(M) after we test to make sure it is a valid tracer
!       number (bmy, 4/16/99)
!  (4) For a CH3I run, make sure that the tracer number N is not larger
!       than NTRACE (bmy, 4/9/99) 
!  (5) ND36 now uses the AD36 array instead of AIJ. (bmy, 3/16/00)
!  (6) Rewritten for clarity; also fixed for CH3I (bmy, 7/25/06)
!  (7) Bug fix: given the tracer number, now search for entry in IDEMS
!       to jive with historical baggage (bmy, 3/6/07)
!  (8) Because HEMCO now tracks ND36, print out diagnostics for every
!       emission species that is a defined tracer. (bmy, 8/29/14)
!******************************************************************************
!                     
      IF ( ND36 > 0 ) THEN

         ! Diagnostic category
         CATEGORY     = 'ANTHSRCE'

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         !%%% NOTE: FACTOR is the multiplication factor that converts
         !%%% kg/m2/s (as tracked by HEMCO) to either molec/cm2/s.        
         !%%% (bmy, 8/14/14)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! Now loop over all species and write out every defined
         ! diagnostics (ckeller, 10/15/2014) 
         DO NA = 1, nAdvect

            ! Init
            N = 0

            ! Species ID
            I = State_Chm%Map_Advect(NA)

            ! Point to the species database entry for this species name
            SpcInfo => State_Chm%SpcData(I)%Info

            ! Get variables needed for HEMCO diagnostics
            SELECT CASE( TRIM( SpcInfo%Name ) )
               CASE( 'NO'     )
                  DiagnName = 'ANTHROPOGENIC_NO'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'O3'     )
                  DiagnName = 'ANTHROPOGENIC_O3'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE ( 'CO'    ) 
                  DiagnName = 'ANTHROPOGENIC_CO'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE ( 'COus'  ) 
                  DiagnName = 'ANTHROPOGENIC_TAGCO_US'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE ( 'COeur'  ) 
                  DiagnName = 'ANTHROPOGENIC_TAGCO_EUR'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE ( 'COasia' ) 
                  DiagnName = 'ANTHROPOGENIC_TAGCO_ASIA'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE ( 'COoth'  ) 
                  DiagnName = 'ANTHROPOGENIC_TAGCO_OTHER'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'ALK4'    )
                  DiagnName = 'ANTHROPOGENIC_ALK4'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'ACET'    ) 
                  DiagnName = 'ANTHROPOGENIC_ACET'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'MEK'     )
                  DiagnName = 'ANTHROPOGENIC_MEK'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'ALD2'    )
                  DiagnName = 'ANTHROPOGENIC_ALD2'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'RCHO'    )
                  DiagnName = 'ANTHROPOGENIC_RCHO'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'MACR'    )
                  DiagnName = 'ANTHROPOGENIC_MACR'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'PRPE'    )
                  DiagnName = 'ANTHROPOGENIC_PRPE'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'C3H8'    )
                  DiagnName = 'ANTHROPOGENIC_C3H8'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'CH2O'    )
                  DiagnName = 'ANTHROPOGENIC_CH2O'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'C2H6'    )
                  DiagnName = 'ANTHROPOGENIC_C2H6'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'SO2'     )
                  DiagnName = 'ANTHROPOGENIC_SO2'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'SO4'     )
                  DiagnName = 'ANTHROPOGENIC_SO4'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'NH3'     )
                  DiagnName = 'ANTHROPOGENIC_NH3'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'BCPI'    )
                  DiagnName = 'ANTHROPOGENIC_BCPI'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'OCPI'    )
                  DiagnName = 'ANTHROPOGENIC_OCPI'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'BCPO'    )
                  DiagnName = 'ANTHROPOGENIC_BCPO'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'OCPO'    )
                  DiagnName = 'ANTHROPOGENIC_OCPO'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'NO2'     )
                  DiagnName = 'ANTHROPOGENIC_NO2'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'HNO2'    )
                  DiagnName = 'ANTHROPOGENIC_HNO2'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'CH4'     )
                  DiagnName = 'ANTHROPOGENIC_CH4'
                  UNIT      = 'molec/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'BENZ'    )
                  DiagnName = 'ANTHROPOGENIC_BENZ'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'TOLU'    )
                  DiagnName = 'ANTHROPOGENIC_TOLU'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'XYLE'    )
                  DiagnName = 'ANTHROPOGENIC_XYLE'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE( 'EOH'     )
                  DiagnName = 'ANTHROPOGENIC_EOH'
                  UNIT      = 'atoms C/cm2/s'
                  N         = SpcInfo%ModelId
                  FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) /
     &                        CM2PERM2
               CASE DEFAULT
                  ! Do nothing
            END SELECT

            ! Call only if defined!
            IF ( N > 0 ) THEN
               ! kg/m2/s to molec/cm2/s 
               CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                         UNIT, N, 1, LLPAR, .TRUE., FACTOR, RC )
               IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
            ENDIF

            ! Free pointer
            SpcInfo => NULL()

         ENDDO ! I
      ENDIF
!
!******************************************************************************
!  ND37: Fraction of tracer scavenged in convective cloud updrafts
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WETCVF-$ : Scavenging fraction         : unitless  : SCALECONV
!******************************************************************************
!
      IF ( ND37 > 0 ) THEN
         CATEGORY = 'MC-FRC-$'
         UNIT     = 'unitless'

         ! Loop over soluble tracers
         DO NW = 1, State_Chm%nWetDep

            ! Get the species ID from the wetdep ID
            N = State_Chm%Map_WetDep(NW)

            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= N .AND. MM <= TMAX(37) )
               MMB = TINDEX(37,MM)
               MM  = MM + 1
            ENDDO
            
            IF ( MMB /= N ) CYCLE

            DO L = 1, LD37
               ARRAY(:,:,L) = AD37(:,:,L,NW) / SCALECONV
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD37,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD37) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND38: Rainout loss of tracer in convective updrafts
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WETDCV-$ : Rainout loss of tracer      : kg/s      : SCALECONV
!
!  NOTES:
!  (1) Now write only LD38 levels to bpch file (bmy, 12/7/00)
!******************************************************************************
!
      IF ( ND38 > 0 ) THEN
         CATEGORY = 'WETDCV-$'
         UNIT     = 'kg/s'

         ! Loop over soluble tracers
         DO NW = 1, State_Chm%nWetDep

            ! Get the species ID from the wetdep ID
            N = State_Chm%Map_WetDep(NW)

            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= N .AND. MM <= TMAX(38) )
               MMB = TINDEX(38,MM)
               MM  = MM + 1
            ENDDO
            
            IF ( MMB /= N ) CYCLE

            ! Divide by # of convective timesteps
            DO L = 1, LD38
               ARRAY(:,:,L) = AD38(:,:,L,NW) / SCALECONV
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD38,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD38) )
         ENDDO

      ENDIF

!
!******************************************************************************
!  ND39: Rainout loss of tracer in large scale rains 
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WETDLS-$ : Large-scale loss of tracer  : kg/s      : SCALEDYN
!******************************************************************************
!
      IF ( ND39 > 0 ) THEN
         CATEGORY = 'WETDLS-$'
         UNIT     = 'kg/s'

         ! Loop over soluble species
         DO NW = 1, State_Chm%nWetDep
               
            ! Get the species ID from the wetdep ID
            N = State_Chm%Map_WetDep(NW)

            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= N .AND. MM <= TMAX(39) )
               MMB = TINDEX(39,MM)
               MM  = MM + 1
            ENDDO
               
            IF ( MMB /= N ) CYCLE

            ! Divide by # of wetdep (= dynamic) timesteps
            DO L = 1, LD39
               ARRAY(:,:,L) = AD39(:,:,L,NW) / SCALEDYN
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD39,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD39) )
         ENDDO

      ENDIF
!
!******************************************************************************
!  ND41: Afternoon boundary layer heights
!******************************************************************************
!
      IF ( ND41 > 0 ) CALL WRITE_DIAG41
!
!******************************************************************************
!  ND42: SOA concentrations [ug/m3]
!******************************************************************************
!
      IF ( ND42 > 0 ) CALL WRITE_DIAG42( Input_Opt )
!
!******************************************************************************
!  ND42: Free diagnostic as of 11/24/99
!
!  ND43: Chemical production of OH and NO
!
!   # : Field  : Description                : Units   : Scale Factor
!  ---------------------------------------------------------------------------
!  (1)  OH     : OH     Chemical Diagnostic : molec/cm3 : CTOH
!  (3)  HO2    : HO2    Chemical Diagnostic : v/v       : CTHO2
!  (4)  O1D    : O2D    Chemical Diagnostic : molec/cm3 : CTO1D
!  (5)  O3P    : O3P    Chemical Diagnostic : molec/cm3 : CTO3P
!
!
!  NOTES:
!  (1) Print output for either a NOx-Ox-HC run (NSRCX == 3), or a CO run
!       with parameterized OH (NSRCX== 5).  (bmy, 4/17/00)
!  (2) Add parentheses in IF test since .AND. has higher precedence
!       than .OR. (jsw, bmy, 12/5/00)
!  (3) Added HO2, NO2 to ND43 (rvm, bmy, 2/27/02)
!  (4) Added NO3 to ND43 (bmy, 1/16/03)
!  (5) Now uses 3D counters (phs, 1/24/07)
!  (6) Now assume that LD43 can't be higher than LD45 (phs, 1/24/07)
!  (7) Check that CTxx are not zero, instead of adding 1e-20 (phs, 11/13/07)
!  (8) Removed NO, NO2, and NO3. These are now tracers (mpayer, 3/29/13)
!  24 Feb 2014 - M. Sulprizio- Added O1D and O3P for UCX strat chem
!  03 Jan 2018 - M. Sulprizio- Replace UCX CPP switch with Input_Opt%LUCX
!******************************************************************************
!
      IF ( ND43 > 0 .and. ( Input_Opt%ITS_A_FULLCHEM_SIM .or. 
     &                      Input_Opt%ITS_A_CH4_SIM           ) ) THEN

         CATEGORY = 'CHEM-L=$' 

         ! Default units
         UNIT     = 'molec/cm3'

         DO M = 1, 2
            N  = TINDEX(43,M)

            SELECT CASE ( N )

               ! OH
               CASE ( 1 )
                  WHERE( CTOH /= 0 )
                     ARRAY(:,:,1:LD43) = AD43(:,:,1:LD43,N) /
     &                    FLOAT( CTOH )
                  ELSEWHERE
                     ARRAY(:,:,1:LD43) = 0.
                  ENDWHERE

                  UNIT = 'molec/cm3'
                  NN   = 1             ! backwards compat. w/ prior code

               ! HO2 (rvm, bmy, 2/27/02)
               CASE ( 2 )
                  WHERE( CTHO2 /= 0 )
                     ARRAY(:,:,1:LD43) = AD43(:,:,1:LD43,N) /
     &                    FLOAT( CTHO2 )
                  ELSEWHERE
                     ARRAY(:,:,1:LD43) = 0.
                  ENDWHERE

                  UNIT = 'v/v'
                  NN   = 3             ! backwards compat. w/ prior code

               CASE DEFAULT
                  ! Nothing

            END SELECT

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD43,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD43) )
         ENDDO


         IF ( Input_Opt%LUCX ) THEN

            !--------------------------------------------------------
            !           %%%%%%% UCX-based mechanisms %%%%%%%
            !--------------------------------------------------------

            DO M = 3,4
               N  = TINDEX(43,M)

               SELECT CASE ( N )

                  ! O1D
                  CASE ( 3 )
                     WHERE( CTO1D /= 0 )
                        ARRAY(:,:,1:LD43) = AD43(:,:,1:LD43,N) /
     &                       FLOAT( CTO1D )
                     ELSEWHERE
                        ARRAY(:,:,1:LD43) = 0.
                     ENDWHERE

                     UNIT = 'molec/cm3'
                     NN   = 4   ! backwards compat. w/ prior code

                  ! O3P
                  CASE ( 4 )
                     NN = 5
                     WHERE( CTO3P /= 0 )
                        ARRAY(:,:,1:LD43) = AD43(:,:,1:LD43,N) /
     &                       FLOAT( CTO3P )
                     ELSEWHERE
                        ARRAY(:,:,1:LD43) = 0.
                     ENDWHERE

                     UNIT = 'molec/cm3'
                     NN   = 5   ! backwards compat. w/ prior code

                  CASE DEFAULT
                     ! Nothing

               END SELECT

               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                     IIPAR,     JJPAR,     LD43,     IFIRST,     
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1:LD43) )
            ENDDO

         ENDIF ! LUCX

      ENDIF
!
!******************************************************************************
!  ND44: Drydep flux (molec/cm2/s) and velocity (cm/s) diagnostics
!
!   #   : Field    : Quantity           : Units               : Scale factor
!  -------------------------------------------------------------------------
!  (1 ) : DRYD-FLX : drydep fluxes      : molec/cm2/s or kg/s : SCALECHEM
!  (2 ) : DRYD-VEL : drydep velocities  : cm/s                : SCALECHEM
!
!  NOTES: 
!  (1 ) Remove diagnostics for wetdep HNO3, H2O2 from ND44.
!  (2 ) For NSRCX == 1 (Rn-Pb-Be), save the actual tracer number 
!        instead of the dry deposition index.  Add TRCOFFSET to N.
!  (3 ) For NSRCX == 6 (single tracer Ox), drydep fluxes are in kg/s.
!  (4 ) ND44 now uses allocatable array AD44 instead of AIJ. (bmy, 3/16/00)
!  (5 ) Add code from amf for multi-tracer Ox (bmy, 7/3/01)
!  (6 ) Now divide by SCALECHEM since DRYFLX is only called after the
!        chemistry routines for all relevant simulations (bmy, 1/27/03)
!  (7 ) Now print out NTRACE drydep fluxes for tagged Ox.  Also tagged Ox 
!        now saves drydep in molec/cm2/s. (bmy, 8/19/03)
!  (8 ) Rearrange ND44 code for clarity (bmy, 3/24/04)
!  (9 ) Add code for H2/HD simulation (phs, 5/8/07)
!  29 Jan 2014 - R. Yantosca - Now save DEPNAME(N) into DRYDEP_NAME, which
!                              helps avoid ambiguity for TOMAS simulations
!  12 Jun 2015 - R. Yantosca - ND44 is now always in molec/cm2/s
!  18 Jul 2016 - M. Sulprizio- Remove special handling of ISOPN and MMN. Family
!                              tracers have been eliminated.
!******************************************************************************
!
      IF ( ND44 > 0 ) THEN
         
         !==============================================================
         ! Drydep fluxes
         !==============================================================

         ! Category name
         CATEGORY = 'DRYD-FLX'

         ! # of drydep flux tracers
         IF ( Input_Opt%ITS_A_TAGO3_SIM .or. 
     &        Input_Opt%ITS_A_MERCURY_SIM ) THEN
            M = nAdvect
         ELSE
            ! Extend dry dep tracers if TOMAS aerosol is turned on
#if   defined ( TOMAS )
            IF ( id_NK1 > 0 ) THEN
               M = nDryDep + ( ( ICOMP - IDIAG )* IBINS )
            ELSE
               M = nDryDep
            ENDIF
#else
            M = nDryDep
#endif
         ENDIF

         ! Loop over drydep species
         DO N = 1, M

            IF ( Input_Opt%ITS_A_RnPbBe_SIM .or. 
     &           Input_Opt%ITS_A_H2HD_SIM ) THEN

               ! NOTE: ND44 is now archived in molec/cm2/s for all 
               ! simulations, including Rn-Pb-Be. (bmy, 6/16/15)
               ! Radon or H2/HD
               UNIT = 'molec/cm2/s'
               NN   = State_Chm%Map_DryDep(N)
 
            ELSE IF ( Input_Opt%ITS_A_TAGO3_SIM .or. 
     &                Input_Opt%ITS_A_MERCURY_SIM ) THEN

               ! Tagged O3 or Tagged Hg
               UNIT = 'molec/cm2/s'
               NN   = N
  
            ELSE 
     
               ! Other simulations
               UNIT = 'molec/cm2/s'

#if   defined( TOMAS )
               ! For extended drydep tracers, assign tracer ID of the TOMAS
               ! aerosol mass (win, 7/14/09)
               IF ( N <= nDryDep ) THEN 
                  NN  = State_Chm%Map_DryDep(N)
                  NN1 = NN
               ELSE
                  ! To calculate the id_xxx of the associated  
                  ! tracer. (ccc, 3/11/10) 
                  NN  = MOD( N - nDryDep-1, IBINS ) + id_NK1

                  ! Tracer number for bpch file
                  NN1  = ( N - nDryDep ) + ( id_NK1 + IBINS - 1 )
               ENDIF
#else
               NN   = State_Chm%Map_DryDep(N)
#endif
 
            ENDIF

            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
#if !defined( TOMAS )
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= NN .AND. MM <= TMAX(44) )
               MMB = TINDEX(44,MM)
               MM  = MM + 1
            ENDDO
#endif

            ! Save into ARRAY
            ARRAY(:,:,1) = ( AD44(:,:,N,1) / SCALECHEM )

            ! Write to file
#if   defined( TOMAS )
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN1,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
#else
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
#endif

         ENDDO

         !==============================================================
         ! Drydep velocities
         !==============================================================

         ! Category and Unit 
         CATEGORY  = 'DRYD-VEL'
         UNIT      = 'cm/s'

         ! # of drydep velocity tracers
         IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
            M = 1
         ELSE IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN  ! Ask Helen
            M = 3
         ELSE
            M = nDryDep
         ENDIF

         ! Loop over drydep tracers
         DO N = 1, M

            NN  = State_Chm%Map_DryDep(N)
            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= NN .AND. MM <= TMAX(44) )
               MMB = TINDEX(44,MM)
               MM  = MM + 1
            ENDDO

            ! Tracer number plus GAMAP offset
            ARRAY(:,:,1) = AD44(:,:,N,2) / SCALESRCE

            ! Write to file
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

         ENDDO
      ENDIF
!
!******************************************************************************
!  ND45: Tracer Mixing Ratio (v/v) for Levels L = 1, LD45
!        averaged between hours OTH_HR1 and OTH_HR2
!
!   # : Field   : Description            : Units   : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) IJ-AVG-$ : Tracer mixing ratio    : v/v     : CTOTH
!
!  NOTES:
!  (1) For NSRCX == 3 (NOx-Ox-HC run), store pure O3 with index NTRACE+1.
!  (2) Now store pure O3 as NNPAR+1 (now tracer #32). (bmy, 1/10/03)
!  (3) Now uses CTO3 instead of CTOH for pure O3 (phs, 1/24/07)
!  (4) Better handling of O3 case (phs, 11/17/08)
!  (6) Removed code for storing pure O3 as NNPAR+1 because O3 is now a tracer
!      (mpayer, 3/14/13)
!******************************************************************************
!
      IF ( ND45 > 0 ) THEN
         CATEGORY    = 'IJ-AVG-$'
         SCALE_TMP   = FLOAT( CTOTH ) + 1e-20_fp
         UNIT        = ''   

         DO M = 1, TMAX(45)
            N  = TINDEX(45,M)
            IF ( N > nAdvect ) CYCLE !(sfarina, 1/28/2013 - valid for all sims)
            NN = N

            DO L = 1, LD45
               ARRAY(:,:,L) = AD45(:,:,L,N) / SCALE_TMP(:,:)
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD45,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD45) )

         ENDDO

      ENDIF
!
!******************************************************************************
!  ND46: Biogenic source diagnostic
!
!   # : Field  : Description                : Units         : Scale Factor
!  ---------------------------------------------------------------------------
!  (1 ) ISOP   : Isoprene                   : atoms C/cm2/s : SCALESRCE
!  (2 ) ACET   : Acetone                    : atoms C/cm2/s : SCALESRCE
!  (3 ) PRPE   : Propene                    : atoms C/cm2/s : SCALESRCE
!  (4 ) MONX   : Total Monoterpenes         : atoms C/cm2/s : SCALESRCE
!  (5 ) MBOX   : Methyl Butenol             : atoms C/cm2/s : SCALESRCE
!  (6 ) C2H4   : Ethene                     : atoms C/cm2/s : SCALESRCE
!  (7 ) APIN   : Alpha pinene               : atoms C/cm2/s : SCALESRCE
!  (8 ) BPIN   : Beta pinene                : atoms C/cm2/s : SCALESRCE
!  (9 ) LIMO   : Limonene                   : atoms C/cm2/s : SCALESRCE
!  (10) SABI   : Sabinene                   : atoms C/cm2/s : SCALESRCE
!  (11) MYRC   : Mycrene                    : atoms C/cm2/s : SCALESRCE
!  (12) CARE   : 3-Carene                   : atoms C/cm2/s : SCALESRCE
!  (13) OCIM   : Ocimene                    : atoms C/cm2/s : SCALESRCE
!  (14) FAXX   : Formic acid                :   molec/cm2/s : SCALESRCE
!  (15) AAXX   : Acetic acid                :   molec/cm2/s : SCALESRCE
!  (16) ALD2   : Acetaldehyde               : atoms C/cm2/s : SCALESRCE
!  (17) OMON   : Other monoterp             : atoms C/cm2/s : SCALESRCE
!  (18) MOHX   : Methanol                   :   molec/cm2/s : SCALESRCE
!  (19) EOH    : Ethanol                    : atoms C/cm2/s : SCALESRCE
!  (20) FARN   : Farnesene                  : atoms C/cm2/s : SCALESRCE
!  (21) BCAR   : b-caryophyllene            : atoms C/cm2/s : SCALESRCE
!  (22) OSQT   : Other sesquiterpenes       : atoms C/cm2/s : SCALESRCE
!  (23) CHBr3  : Bromoform                  :       kg/m2/s : SCALESRCE
!  (24) CH2Br2 : Dibromomethane             :       kg/m2/s : SCALESRCE
!  (25) Br2    : Molec. Bromine             :       kg/m2/s : SCALESRCE
!  (26) CH4    : Methane                    :   molec/cm2/s : SCALE3
!  (27) ALD2   : ALD2 from senescing plants : atoms C/cm2/s : SCALESRCE
!  (28) EOH    : EOH  from senescing plants : atoms C/cm2/s : SCALESRCE
!  (29) ALD2   : ALD2 ocean emissions       : atoms C/cm2/s : SCALESRCE
!
!  NOTES:
!  (1) ND46 now uses allocatable array AD46 instead of AIJ (bmy, 3/16/00)
!  (2) Also write out PRPE for CO-OH run (NSRCX == 5), regardless of
!       the setting of id_PRPE.  This is to print out monterpene 
!       diagnostics. (bnd, bmy, 4/18/00)
!  (3) Added monoterpenes as tracer #4.  This requires updated versions
!       of "tracerinfo.dat" and "diaginfo.dat" for GAMAP. (bmy, 1/2/01)
!  (4) Added MBO as tracer #5. (tmf, bmy, 10/20/05)
!  (5) Added C2H4 as tracer #6. (tmf, 1/20/09)
!  13 Feb 2015 - M. Sulprizio- Add updates for MEGAN2.1 from D. Millet,
!                              based on Guenther et al. 2012
!******************************************************************************
!
      IF ( ND46 > 0 ) THEN

         ! Diagnostic category and units for BPCH output
         CATEGORY = 'BIOGSRCE'

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         !%%% NOTE: FACTOR is the multiplication factor that converts
         !%%% kg/m2/s (as tracked by HEMCO) to either molec/cm2/s.        
         !%%% (bmy, 8/14/14)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !-----------------------------------------------------------
         !%%%%% ISOPRENE %%%%%
         !-----------------------------------------------------------
         IF ( id_ISOP > 0 ) THEN
            DiagnName = 'BIOGENIC_ISOP'
            UNIT      = 'atoms C/cm2/s'
            N         = 1         
            SpcInfo => State_Chm%SpcData(id_ISOP)%Info
            ! Conversion factor from kgC/m2/s to atoms C/cm2/s
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

            ! Free pointer
            SpcInfo => NULL()
         ENDIF

         !-----------------------------------------------------------
         !%%%%% ACETONE %%%%%
         !-----------------------------------------------------------
         IF ( id_ACET > 0 ) THEN
            UNIT   = 'atoms C/cm2/s'
            N      = 2
            SpcInfo => State_Chm%SpcData(id_ACET)%Info
            ! Conversion factor from kgC/m2/s to atoms C/cm2/s
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            ! (1) BIOGENIC ACET 
            DiagnName = 'BIOGENIC_ACET'
            DiagnCnt  => NULL()

            ! Get diagnostics from HEMCO
            CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                      FLAG, RC,   cName=TRIM(DiagnName), 
     &                      AutoFill=1, 
     &                      COL=HcoState%Diagn%HcoDiagnIDManual )
            IF ( RC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
               CALL ERROR_STOP ( MSG, LOC ) 
            ENDIF

            ! Save into ARRAY. Convert to units
            IF ( FLAG == HCO_SUCCESS ) THEN
               ARRAY(:,:,1) = DiagnCnt%Arr2D%Val(:,:) * FACTOR 
            ELSE
               ARRAY(:,:,1) = 0.0
               MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
               MSG = TRIM(MSG) // ' - will write zeros!'
               CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
            ENDIF

            ! (2) Then get oceanic acetone and add to it 
            DiagnName = 'AD11_OCEAN_SOURCE'
            DiagnCnt  => NULL()

            ! Get diagnostics from HEMCO
            CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                      FLAG, RC,   cName=TRIM(DiagnName), 
     &                      AutoFill=1, 
     &                      COL=HcoState%Diagn%HcoDiagnIDManual )
            IF ( RC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
               CALL ERROR_STOP ( MSG, LOC ) 
            ENDIF

            ! Save into ARRAY. Convert to units
            IF ( FLAG == HCO_SUCCESS ) THEN
                ARRAY(:,:,1) = ARRAY(:,:,1) +
     &                      ( DiagnCnt%Arr2D%Val(:,:) * FACTOR )
            ELSE
               MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
               MSG = TRIM(MSG) // ' - will write zeros!'
               CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
            ENDIF

            ! Write combined BIOGENIC ACET and OCEAN ACET to disk
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT ,     DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

            ! Free pointer
            SpcInfo => NULL()
         ENDIF

         !-----------------------------------------------------------
         !%%%%% ALKENES %%%%%
         !-----------------------------------------------------------
         IF ( id_PRPE > 0 ) THEN
            DiagnName = 'BIOGENIC_PRPE'
            UNIT      = 'atoms C/cm2/s'
            N         = 3
            SpcInfo => State_Chm%SpcData(id_PRPE)%Info
            ! Conversion factor from kgC/m2/s to atoms C/cm2/s
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

            ! Free pointer
            SpcInfo => NULL()
         ENDIF

         !-----------------------------------------------------------
         !%%%%% TOTAL MONOTERPENES %%%%%
         !
         ! APIN + BPIN + LIMO + SABI + MYRC + CARE + OCIM + OMON
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_MONX'
         UNIT      = 'atoms C/cm2/s'
         N         = 4
         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR    = GPERKG / MWC * AVO / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% METHYL BUTENOL %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_MBOX'
         UNIT      = 'atoms C/cm2/s'
         N         = 5
         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR    = GPERKG / MWC * AVO / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% ETHENE %%%%%
         !-----------------------------------------------------------
         IF ( id_C2H4 > 0 ) THEN
            DiagnName = 'BIOGENIC_C2H4'
            UNIT      = 'atoms C/cm2/s'
            N         = 6
            SpcInfo => State_Chm%SpcData(id_C2H4)%Info
            ! Conversion factor from kgC/m2/s to atoms C/cm2/s
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

            ! Free pointer
            SpcInfo => NULL()
         ENDIF
 
         !-----------------------------------------------------------
         !%%%%% ALPHA-PINENE %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_APIN'
         UNIT      = 'atoms C/cm2/s'
         N         = 7
         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR    = GPERKG / MWC * AVO / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% BETA-PINENE %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_BPIN'
         UNIT      = 'atoms C/cm2/s'
         N         = 8
         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR    = GPERKG / MWC * AVO / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% LIMONENE %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_LIMO'
         UNIT      = 'atoms C/cm2/s'
         N         = 9
         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR    = GPERKG / MWC * AVO / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% SABINENE %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_SABI'
         UNIT      = 'atoms C/cm2/s'
         N         = 10
         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR    = GPERKG / MWC * AVO / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% MYRCENE %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_MYRC'
         UNIT      = 'atoms C/cm2/s'
         N         = 11
         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR    = GPERKG / MWC * AVO / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% 3-CARENE %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_CARE'
         UNIT      = 'atoms C/cm2/s'
         N         = 12
         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR    = GPERKG / MWC * AVO / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% OCIMENE %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_OCIM'
         UNIT      = 'atoms C/cm2/s'
         N         = 13
         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR    = GPERKG / MWC * AVO / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% FORMIC ACID %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_FAXX'
         UNIT      = 'molec/cm2/s'
         N         = 14
         ! Conversion factor from kg/m2/s to molec/cm2/s
         FACTOR    = ( AVO / 46.025e+0_fp ) / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% ACETIC ACID %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_AAXX'
         UNIT      = 'molec/cm2/s'
         N         = 15
         ! Conversion factor from kg/m2/s to molec/cm2/s
         FACTOR    = ( AVO / 60.052e+0_fp ) / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% ACETALDEHYDE %%%%%
         !-----------------------------------------------------------
         IF ( id_ALD2 > 0 ) THEN
            DiagnName = 'BIOGENIC_ALD2'
            UNIT      = 'atoms C/cm2/s'
            N         = 16
            SpcInfo => State_Chm%SpcData(id_ALD2)%Info
            ! Conversion factor from kgC/m2/s to atoms C/cm2/s
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

            ! Free pointer
            SpcInfo => NULL()
         ENDIF

         !-----------------------------------------------------------
         !%%%%% OTHER MONOTERPENES %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_OMON'
         UNIT      = 'atoms C/cm2/s'
         N         = 17
         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR    = GPERKG / MWC * AVO / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% METHANOL %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_MOHX'
         UNIT      = 'molec/cm2/s'
         N         = 18
         ! Conversion factor from kg/m2/s to molec/cm2/s
         FACTOR    =  ( AVO / 32.042e+0_fp )/ CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% ETHANOL %%%%%
         !-----------------------------------------------------------
         IF ( id_EOH > 0 ) THEN
            DiagnName = 'BIOGENIC_EOH'
            UNIT      = 'atoms C/cm2/s'
            N         = 19
            SpcInfo => State_Chm%SpcData(id_eoh)%Info
            ! Conversion factor from kgC/m2/s to atoms C/cm2/s
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

            ! Free pointer
            SpcInfo => NULL()
         ENDIF

         !-----------------------------------------------------------
         !%%%%% FARNESENE %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_FARN'
         UNIT      = 'atoms C/cm2/s'
         N         = 20
         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR    = GPERKG / MWC * AVO / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% B-CARYOPHYLLENE %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_BCAR'
         UNIT      = 'atoms C/cm2/s'
         N         = 21
         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR    = GPERKG / MWC * AVO / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% OTHER SESQUITERPENES %%%%%
         !-----------------------------------------------------------
         DiagnName = 'BIOGENIC_OSQT'
         UNIT      = 'atoms C/cm2/s'
         N         = 22
         ! Conversion factor from kgC/m2/s to atoms C/cm2/s
         FACTOR    = GPERKG / MWC * AVO / CM2PERM2
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !-----------------------------------------------------------
         !%%%%% CHBr3 %%%%%
         !-----------------------------------------------------------
         IF ( id_CHBr3 > 0 ) THEN
            DiagnName = 'BIOGENIC_CHBR3'
            UNIT      = 'kg/m2/s'
            N         = 23
            FACTOR    = 1.0e+0_fp
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF
 
         !-----------------------------------------------------------
         !%%%%% CH2Br2 %%%%%
         !-----------------------------------------------------------
         IF ( id_CH2Br2 > 0 ) THEN
            DiagnName = 'BIOGENIC_CH2BR2'
            UNIT      = 'kg/m2/s'
            N         = 24
            FACTOR    = 1.0e+0_fp
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF
 
         !-----------------------------------------------------------
         !%%%%% Br2 %%%%%
         !-----------------------------------------------------------
         IF ( id_Br2 > 0 ) THEN
            DiagnName = 'SEASALT_BR2'
            UNIT      = 'kg/m2/s'
            N         = 25
            FACTOR    = 1.0e+0_fp
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDIF

         !-----------------------------------------------------------
         !%%%%% CH4 %%%%%
         !
         ! Total biogenic CH4 = wetlands + other natural - soil absorb
         !-----------------------------------------------------------
         IF ( id_CH4 > 0 ) THEN
            UNIT      = 'molec/cm2/s'
            N         = 26
            SpcInfo => State_Chm%SpcData(id_CH4)%Info
            ! Conversion factor from kg/m2/s to molec/cm2/s
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            ! (1) Wetlands
            DiagnName = 'CH4_WETLANDS'
            DiagnCnt  => NULL()

            ! Get diagnostics from HEMCO
            CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                      FLAG, RC,   cName=TRIM(DiagnName), 
     &                      AutoFill=0, 
     &                      COL=HcoState%Diagn%HcoDiagnIDManual )
            IF ( RC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
               CALL ERROR_STOP ( MSG, LOC ) 
            ENDIF

            ! Save into ARRAY and convert units from kg/m2/s to molec/cm2/s
            IF ( FLAG == HCO_SUCCESS ) THEN
               ARRAY(:,:,1) = DiagnCnt%Arr2D%Val(:,:) * FACTOR 
            ELSE
               ARRAY(:,:,1) = 0.0
               MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
               MSG = TRIM(MSG) // ' - will write zeros!'
               CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
            ENDIF

            ! (2) Other natural
            DiagnName = 'CH4_OTHERNATUR'
            DiagnCnt  => NULL()

            ! Get diagnostics from HEMCO
            CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                      FLAG, RC,   cName=TRIM(DiagnName), 
     &                      AutoFill=0, 
     &                      COL=HcoState%Diagn%HcoDiagnIDManual )
            IF ( RC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
               CALL ERROR_STOP ( MSG, LOC ) 
            ENDIF

            ! Save into ARRAY and convert units from kg/m2/s to molec/cm2/s
            IF ( FLAG == HCO_SUCCESS ) THEN
               ARRAY(:,:,1) = ARRAY(:,:,1) +
     &                      ( DiagnCnt%Arr2D%Val(:,:) * FACTOR )
            ELSE
               MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
               MSG = TRIM(MSG) // ' - will write zeros!'
               CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
            ENDIF

            ! (3) Soil absorption (sink)
            DiagnName = 'CH4_SOILABSORB'
            DiagnCnt  => NULL()

            ! Get diagnostics from HEMCO
            CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                      FLAG, RC,   cName=TRIM(DiagnName), 
     &                      AutoFill=0, 
     &                      COL=HcoState%Diagn%HcoDiagnIDManual )
            IF ( RC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
               CALL ERROR_STOP ( MSG, LOC ) 
            ENDIF

            ! Subtract from ARRAY and convert units from kg/m2/s to molec/cm2/s
            IF ( FLAG == HCO_SUCCESS ) THEN
               ARRAY(:,:,1) = ARRAY(:,:,1) -
     &                      ( DiagnCnt%Arr2D%Val(:,:) * FACTOR )
            ELSE
               MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
               MSG = TRIM(MSG) // ' - will write zeros!'
               CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
            ENDIF

            ! Write combined biogenic CH4 emissions to disk
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT ,     DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

            ! Free pointer
            SpcInfo => NULL()
         ENDIF

         !-----------------------------------------------------------
         !%%%%% ACETALDEHYDE FROM SENESCING PLANTS %%%%%
         !-----------------------------------------------------------
         IF ( id_ALD2 > 0 ) THEN
            DiagnName = 'ALD2_SENESCING'
            UNIT      = 'atoms C/cm2/s'
            N         = 27
            SpcInfo  => State_Chm%SpcData(id_ALD2)%Info
            ! Conversion factor from kgC/m2/s to atoms C/cm2/s
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

            ! Free pointer
            SpcInfo => NULL()
         ENDIF

         !-----------------------------------------------------------
         !%%%%% ETHANOL FROM SENESCING PLANTS %%%%%
         !-----------------------------------------------------------
         IF ( id_EOH > 0 ) THEN
            DiagnName = 'EOH_SENESCING'
            UNIT      = 'atoms C/cm2/s'
            N         = 28
            SpcInfo  => State_Chm%SpcData(id_EOH)%Info
            ! Conversion factor from kgC/m2/s to atoms C/cm2/s
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

            ! Free pointer
            SpcInfo => NULL()
         ENDIF

         !-----------------------------------------------------------
         !%%%%% ACETALDEHYDE OCEAN SOURCE %%%%%
         !-----------------------------------------------------------
         IF ( id_ALD2 > 0 ) THEN
            DiagnName = 'ALD2_OCEAN_SOURCE'
            UNIT      = 'atoms C/cm2/s'
            N         = 29
            SpcInfo  => State_Chm%SpcData(id_ALD2)%Info
            ! Conversion factor from kgC/m2/s to atoms C/cm2/s
            FACTOR    = AVO / ( SpcInfo%emMW_g * 1.e-3_fp ) / CM2PERM2
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE.,    FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

            ! Free pointer
            SpcInfo => NULL()
         ENDIF

      ENDIF
!
!******************************************************************************
!  ND47: Tracer Mixing Ratio (v/v) for Levels L = 1, LD47
!        *always* averaged between 0000 and 2400 Local time.
!
!   # : Field   : Description                 : Units   : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) IJ-24H-$ : 24h avg Tracer mixing ratio : v/v     : SCALEDYN
!
!  NOTES:
!  (1) For NSRCX == 3 (NOx-Ox-HC run), store pure O3 with index NTRACE+1.
!  (2) Now store pure O3 as NNPAR+1 (now tracer #32). (bmy, 1/10/03) 
!  (3) Now replace SCALE1 with SCALEDYN
!  (4) Now averaged between 0 and 24 UT.  Replace SCALEDYN with CTOH and
!       CTO3 (phs, 1/24/07)
!  (5) Revert to SCALEDYN for all species, except O3, which uses new 
!       CTO3_24h counter (phs, 11/17/08)
!  (6) Removed code for storing pure O3 as NNPAR+1 because O3 is now a tracer
!      (mpayer, 3/14/13)
!******************************************************************************
!
      IF ( ND47 > 0 ) THEN
         CATEGORY = 'IJ-24H-$'
         UNIT     = ''   

         DO M = 1, TMAX(47)
            N  = TINDEX(47,M)
            IF ( N > nAdvect ) CYCLE
            NN = N
            
            DO L = 1, LD47
               ARRAY(:,:,L) = AD47(:,:,L,N) / SCALEDIAG
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD47,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD47) )

         ENDDO
      ENDIF

!******************************************************************************
!  ND52: Uptake coefficients (gamma)
!   # Category 
!   # : Field     : Description                      : Units     : Scale factor
!  ----------------------------------------------------------------------------
!  (1): GAMMAHO2  : Uptake coef for HO2              : unitless  : SCALECHEM
!  (2): GAMMAIEPOX: Uptake coef for IEPOX            : unitless  : SCALECHEM
!  (3): GAMMAIEPOX: Uptake coef for IMAE             : unitless  : SCALECHEM
!  (4): GAMMAIEPOX: Uptake coef for ISOPND/ISOPNB    : unitless  : SCALECHEM
!  (5): GAMMAIEPOX: Uptake coef for DHDN             : unitless  : SCALECHEM
!  (6): GAMMAIEPOX: Uptake coef for GLYX             : unitless  : SCALECHEM
!******************************************************************************
      IF ( ND52 > 0 ) THEN
        CATEGORY  = 'GAMMA-$'
        UNIT      = 'unitless'

         DO M = 1, TMAX(52)
            N = TINDEX(52,M)
            IF ( N > PD52 ) CYCLE
            NN = N

            DO L = 1, LD52
               ARRAY(:,:,L) = AD52(:,:,L,N) / SCALECHEM
            ENDDO

            ! Save to disk
            CALL BPCH2( IU_BPCH,    MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR,  CENTER180, CATEGORY, NN,
     &                  UNIT,       DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,      JJPAR,     LD52,     IFIRST,
     &                  JFIRST,     LFIRST,    ARRAY(:,:,1:LD52) )
         ENDDO

      ENDIF
!
!******************************************************************************
!  ND53: POPs emissions (eck, 10/21/10)
!
!   # : Field  : Description                    : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (1 ) POPT     : Emissions of total POPs        : kg       : 1
!  (2 ) OCPO em  : Emissions of OC-sorbed POP     : kg       : 1
!  (3 ) BCPO em  : Emissions of BC-sorbed POP     : kg       : 1
!  (4 ) G em     : Emissions of gas-phase POP     : kg       : 1
!  (5 ) kgGsoil  : Secondary emissions from soil  : kg       : 1
!  (6 ) kgGlake  : Secondary emissions from lake  : kg       : 1
!  (7 ) kgGleaf  : Secondary emissions from leaf  : kg       : 1
!  (8 ) soi2air  : Secondary positive soil flux   : ng/m2/day: 1
!  (9 ) air2soi  : Secondary negative soil flux   : ng/m2/day: 1
!  (10) lak2air  : Secondary positive lake flux   : ng/m2/day: 1
!  (11) air2lak  : Secondary negative lake flux   : ng/m2/day: 1
!  (12) lea2air  : Secondary positive leaf flux   : ng/m2/day: 1
!  (13) air2lea  : Secondary negative leaf flux   : ng/m2/day: 1
!  (14) soi/air  : Fugacity ratios soil/air       : unitless : 1
!  (15) lak/air  : Fugacity ratios lake/air       : unitless : 1
!  (16) lea/air  : Fugacity ratios leaf/air       : unitless : 1
!  (17) OCPO2G   : Gross POP OC lost to gas       : kg       : 1
!  (18) G2OCPO   : Gross POP OC formed from gas   : kg       : 1
!  (19) BCPO2G   : Gross POP BC lost to gas       : kg       : 1
!  (20) G2BCPO   : Gross POP BC formed from gas   : kg       : 1
!  (21) OHoxG    : Prod of POPG from rxn with OH  : kg       : 1
!  (22) O3oxOCPO : Prod of POPOCPO from rxn w/ O3 : kg       : 1
!  (23) O3oxOCPI : Prod of POPOCPI from rxn w/ O3 : kg       : 1
!  (24) O3oxBCPO : Prod of POPBCPO from rxn w/ O3 : kg       : 1
!  (25) O3oxBCPI : Prod of POPBCPI from rxn w/ O3 : kg       : 1
!  (26) NO3xOCPO : Prod of POPOCPO from rxn w/ NO3: kg       : 1
!  (27) NO3xOCPI : Prod of POPOCPI from rxn w/ NO3: kg       : 1
!  (28) NO3xBCPO : Prod of POPBCPO from rxn w/ NO3: kg       : 1
!  (29) NO3xBCPI : Prod of POPBCPI from rxn w/ NO3: kg       : 1
!******************************************************************************
!
      IF ( ND53 > 0 .and. Input_Opt%ITS_A_POPS_SIM ) THEN

         CATEGORY = 'PG-SRCE'

         IF ( id_POPPOCPO > 0 .and.
     &        id_POPPBCPO > 0 .and.
     &        id_POPG     > 0 ) THEN

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Total POP emissions
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            !---------------------------------------------------------
            ! (1) Get OC POP emissions from HEMCO data structure
            !---------------------------------------------------------
            DiagnName = 'GCPOPS_POPPOCPO_SOURCE'
            UNIT      = 'kg/m2/s'
            N         = 1
            DiagnCnt  => NULL()

            ! Zero the diagnostic array
            ARRAY     = 0e0

            ! Get diagnostics from HEMCO
            CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                      FLAG, RC,   cName=TRIM(DiagnName), 
     &                      AutoFill=1, 
     &                      COL=HcoState%Diagn%HcoDiagnIDManual )
            IF ( RC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
               CALL ERROR_STOP ( MSG, LOC ) 
            ENDIF

            ! Save into ARRAY
            IF ( FLAG == HCO_SUCCESS ) THEN
               ARRAY(:,:,1) = DiagnCnt%Arr2D%Val(:,:)
            ELSE
               MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
               MSG = TRIM(MSG) // ' - will write zeros!'
               CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
            ENDIF

            !---------------------------------------------------------
            ! (2) Get BC POP emissions from HEMCO data structure
            !---------------------------------------------------------
            DiagnName = 'GCPOPS_POPPBCPO_SOURCE'
            DiagnCnt  => NULL()

            ! Get diagnostics from HEMCO
            CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt, 
     &                      FLAG, RC,   cName=TRIM(DiagnName), 
     &                      AutoFill=1, 
     &                      COL=HcoState%Diagn%HcoDiagnIDManual )
            IF ( RC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
               CALL ERROR_STOP ( MSG, LOC ) 
            ENDIF

            ! Save into ARRAY
            IF ( FLAG == HCO_SUCCESS ) THEN
               ARRAY(:,:,1) = ARRAY(:,:,1) + DiagnCnt%Arr2D%Val(:,:)
            ELSE
               MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
               MSG = TRIM(MSG) // ' - will write zeros!'
               CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
            ENDIF

            !---------------------------------------------------------
            ! (3) Get gaseous POP emissions from HEMCO data structure
            !---------------------------------------------------------
            DiagnName = 'GCPOPS_POPG_SOURCE'
            DiagnCnt  => NULL()

            ! Get diagnostics from HEMCO
            CALL Diagn_Get( am_I_Root,  HcoState, .FALSE., DiagnCnt, 
     &                      FLAG, RC,   cName=TRIM(DiagnName), 
     &                      AutoFill=1, 
     &                      COL=HcoState%Diagn%HcoDiagnIDManual )
            IF ( RC /= HCO_SUCCESS ) THEN
               MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
               CALL ERROR_STOP ( MSG, LOC ) 
            ENDIF

            ! Save into ARRAY
            IF ( FLAG == HCO_SUCCESS ) THEN
               ARRAY(:,:,1) = ARRAY(:,:,1) + DiagnCnt%Arr2D%Val(:,:)
            ELSE
               MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
               MSG = TRIM(MSG) // ' - will write zeros!'
               CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
            ENDIF

            ! Write combined OC POP + BC POP + Gaseous POP to disk
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT ,     DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

            ! Free pointer
            DiagnCnt => NULL()

         ENDIF

         IF ( id_POPPOCPO > 0 ) THEN

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! OC-sorbed POP emissions
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            DiagnName = 'GCPOPS_POPPOCPO_SOURCE'
            UNIT      = 'kg/m2/s'
            N         = 2
            FACTOR    = 1.0e+0_fp
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         ENDIF

         IF ( id_POPPBCPO > 0) THEN

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! BC-sorbed POP emissions
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            DiagnName = 'GCPOPS_POPPBCPO_SOURCE'
            UNIT      = 'kg/m2/s'
            N         = 3
            FACTOR    = 1.0e+0_fp
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         ENDIF

         IF ( id_POPG > 0 ) THEN

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Gas-phase POP emissions
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            DiagnName = 'GCPOPS_POPG_SOURCE'
            N         = 4
            UNIT      = 'kg/m2/s'
            FACTOR    = 1.0e+0_fp
            CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                      UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
            IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         ENDIF

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Secondary POP emissions from soil
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         DiagnName = 'GCPOPS_POPG_SOIL'
         UNIT      = 'kg/m2/s'
         N         = 5
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Secondary POP emissions from lake
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         DiagnName = 'GCPOPS_POPG_LAKE'
         UNIT      = 'kg/m2/s'
         N         = 6
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Secondary POP emissions from leaf
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         DiagnName = 'GCPOPS_POPG_LEAF'
         UNIT      = 'kg/m2/s'
         N         = 7
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Secondary positive soil flux
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         DiagnName = 'GCPOPS_SOIL2AIR'
         UNIT      = 'ng/m2/day'
         N         = 8
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Secondary negative soil flux
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         DiagnName = 'GCPOPS_AIR2SOIL'
         UNIT      = 'ng/m2/day'
         N         = 9
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Secondary positive lake flux
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         DiagnName = 'GCPOPS_LAKE2AIR'
         UNIT      = 'ng/m2/day'
         N         = 10
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Secondary negative lake flux
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         DiagnName = 'GCPOPS_AIR2LAKE'
         UNIT      = 'ng/m2/day'
         N         = 11
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Secondary positive leaf flux
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         DiagnName = 'GCPOPS_LEAF2AIR'
         UNIT      = 'ng/m2/day'
         N         = 12
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Secondary negative leaf flux
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         DiagnName = 'GCPOPS_AIR2LEAF'
         UNIT      = 'ng/m2/day'
         N         = 13
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Fugacity ratios for soils (soil/air)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         DiagnName = 'GCPOPS_SOILAIR_FUG'
         Unit      = '1'
         N         = 14
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Fugacity ratios for lakes (water/air)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         DiagnName = 'GCPOPS_LAKEAIR_FUG'
         Unit      = '1'
         N         = 15
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Fugacity ratios for leaves (leaf sfc/air)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         DiagnName = 'GCPOPS_LEAFAIR_FUG'
         Unit      = '1'
         N         = 16
         FACTOR    = 1.0e+0_fp
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,      N, 1, -1, .TRUE., FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)

         ! Write the remaining ND53 diagnostics 
         CALL WRITE_DIAG53

      ENDIF
!
!******************************************************************************
!  ND54: Time-in-the-Troposphere diagnostic
!
!   # : Field   : Description                 : Units     : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) TIME-TPS : Time spend in troposphere   : fraction  : SCALEDYN
!******************************************************************************
!
      IF ( ND54 > 0 ) THEN
         CATEGORY = 'TIME-TPS'
         UNIT = 'unitless'

         DO L = 1, LD54
            ARRAY(:,:,L) = AD54(:,:,L) / SCALEDIAG
         ENDDO

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, 1,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD54,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD54) )
      ENDIF
!
!******************************************************************************
!  ND55: Tropopause diagnostics
!
!   # : Field   : Description                 : Units     : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) TP-LEVEL : Tropopause level            : unitless  : SCALEDYN
!  (2) TP-HGHT  : Tropopause height           : km        : SCALEDYN
!  (3) TP-PRESS : Tropopause pressure         : mb        : SCALEDYN
!******************************************************************************
!
      IF ( ND55 > 0 ) THEN
         CATEGORY = 'TR-PAUSE'

         DO M = 1, TMAX(55)
            N  = TINDEX(55,M)
            IF ( N > PD55 ) CYCLE
            NN = N

            ! Pick the appropriate unit string
            SELECT CASE ( N )
               CASE ( 1 )
                  UNIT = 'unitless'
               CASE ( 2 )
                  UNIT = 'km'
               CASE ( 3 )
                  UNIT = 'mb'
            END SELECT

            ARRAY(:,:,1) = AD55(:,:,N) / SCALEDYN

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND56: Lightning flash rate diagnostics (ltm, bmy, 5/5/06))
!******************************************************************************
!
!      IF ( ND56 > 0 ) CALL WRITE_DIAG56

       IF ( ND56 > 0 ) THEN

          CATEGORY = 'LFLASH-$'
          UNIT     = 'flashes/min/km2'

          ! Get scale factor

          ! --> If diagnostics OutOper is set to 'Mean':
          FACTOR = 1.0e+0_fp

          ! Loop over ND56 diagnostic tracers
          DO M = 1, 3

             ! Define quantities
             N = TINDEX(56,M)
            
             SELECT CASE ( N )
                CASE ( 1 )
                   DiagnName = 'LIGHTNING_TOTAL_FLASHRATE'
                CASE ( 2 )
                   DiagnName = 'LIGHTNING_INTRACLOUD_FLASHRATE'
                CASE ( 3 )
                   DiagnName = 'LIGHTNING_CLOUDGROUND_FLASHRATE'
                CASE DEFAULT
                   MSG = 'Lightning index N must not exceed 3!'
                   CALL ERROR_STOP ( MSG, LOC )
             END SELECT

             CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                       UNIT,  N,  -1, -1,   .TRUE.,    FACTOR,
     &                       RC)
             IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)
         ENDDO

         DiagnName = 'LIGHTNING_CLOUD_TOP'
         N         = 4
         CALL DIAG2BPCH( am_I_Root, HcoState, DiagnName, CATEGORY,
     &                   UNIT,  N,  -1, -1,   .TRUE.,    FACTOR, RC )
         IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( DiagnName, LOC)   

      ENDIF
!
!******************************************************************************
!  ND57: THETA
!
!  Potential temperature
!******************************************************************************
! 
! (FP 6/2009)
  
      IF ( ND57 > 0 ) THEN

         CATEGORY          = 'THETA-$'
         UNIT              = 'K'
         ARRAY(:,:,1:LD57) = AD57(:,:,1:LD57) / SCALEDIAG
         NN                = 1

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, NN,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD57,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD57) )
      ENDIF
#if   defined( TOMAS )
!
!*****************************************************************************
!  ND59: Size-resolved primary aerosol emission (win, 8/22/07)
!        Unit is amount of aerosol per time of simulation e.g. sulfate[=] kg S
!
!   # : Field   : Description                 : Units     : Scale Factor
!  ---------------------------------------------------------------------------
!  (1 )  NK-EMISS : Size-resolved number emission   : No       : 1
!  (2 )  SF-EMISS : Size-resolved sulfate emission  : kg S     : 1
!  (3 )  SS-EMISS : Size-resolved sea-salt emission : kg       : 1 
!  (4 )  ECIL-SRC : Size-resolved H-phillic EC emission   : kg : 1
!  (5 )  ECOB-SRC : Size-resolved H-phobic EC emission    : kg : 1
!  (6 )  OCIL-SRC : Size-resolved H-phillic OC emission   : kg : 1
!  (7 )  OCOB-SRC : Size-resolved H-phobic OC emission    : kg : 1
!  (8 )  DUST-SRC : Size-resolved dust emission      : kg      : 1
!
!  NOTES:
!*****************************************************************************
!
      IF ( ND59 > 0 ) THEN

         !==============================================================
         ! Size-resolved primary aerosol number emission 
         !==============================================================
         UNIT         = 'No.'
         CATEGORY     = 'NK-EMISS'
         DO NBIN =1,IBINS
            N        = id_NK1 + NBIN - 1
            DO L = 1, 2
               ARRAY(:,:,L) = AD59_NUMB(:,:,L,NBIN)
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &                  IIPAR,     JJPAR,     2,         IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     
         ENDDO

         !==============================================================
         ! Size-resolved primary aerosol sulfate mass emission 
         !==============================================================
         UNIT         = 'kg S'
         CATEGORY     = 'SF-EMISS'
         DO NBIN =1,IBINS
            N        = id_SF1 + NBIN - 1
            DO L = 1, 2
               ARRAY(:,:,L) = AD59_SULF(:,:,L,NBIN)
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &                  IIPAR,     JJPAR,     2,         IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     
         ENDDO

         !==============================================================
         ! Size-resolved primary aerosol sea-salt mass emission 
         !==============================================================
         UNIT         = 'kg'
         CATEGORY     = 'SS-EMISS'
         DO NBIN =1,IBINS
            N        = id_SS1 + NBIN - 1
            DO L = 1, 2
               ARRAY(:,:,L) = AD59_SALT(:,:,L,NBIN)
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &                  IIPAR,     JJPAR,     2,         IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     
         ENDDO

         !==============================================================
         ! Size-resolved primary Hydrophillic EC mass emission
         !==============================================================
         UNIT         = 'kg'
         CATEGORY     = 'ECIL-SRC'
         DO NBIN =1,IBINS
            N        = id_ECIL1 - 1 + NBIN
            DO L = 1, 2
               ARRAY(:,:,L) = AD59_ECIL(:,:,L,NBIN)
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &                  IIPAR,     JJPAR,     2,         IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     
         ENDDO

         !==============================================================
         ! Size-resolved primary Hydrophobic EC mass emission
         !==============================================================
         UNIT         = 'kg'
         CATEGORY     = 'ECOB-SRC'
         DO NBIN =1,IBINS
            N        = id_ECOB1 - 1 + NBIN
            DO L = 1, 2
               ARRAY(:,:,L) = AD59_ECOB(:,:,L,NBIN)
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &                  IIPAR,     JJPAR,     2,         IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     
         ENDDO

         !==============================================================
         ! Size-resolved primary Hydrophillic OC mass emission
         !==============================================================
         UNIT         = 'kg'
         CATEGORY     = 'OCIL-SRC'
         DO NBIN =1,IBINS
            N        = id_OCIL1 - 1 + NBIN
            DO L = 1, 2
               ARRAY(:,:,L) = AD59_OCIL(:,:,L,NBIN)
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &                  IIPAR,     JJPAR,     2,         IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     
         ENDDO

         !==============================================================
         ! Size-resolved primary Hydrophobic OC mass emission
         !==============================================================
         UNIT         = 'kg'
         CATEGORY     = 'OCOB-SRC'
         DO NBIN =1,IBINS
            N        = id_OCOB1 - 1 + NBIN
            DO L = 1, 2
               ARRAY(:,:,L) = AD59_OCOB(:,:,L,NBIN)
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &                  IIPAR,     JJPAR,     2,         IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     
         ENDDO

         !==============================================================
         ! Size-resolved dust mass emission
         !==============================================================
         UNIT         = 'kg'
         CATEGORY     = 'DUST-SRC'
         DO NBIN =1,IBINS
            N        = id_DUST1 - 1 + NBIN
            DO L = 1, 2
               ARRAY(:,:,L) = AD59_DUST(:,:,L,NBIN)
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &                  IIPAR,     JJPAR,     2,         IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     
         ENDDO

      ENDIF  
!
!*****************************************************************************
!  ND60: TOMAS Aerosol microphysical rate (win, 7/9/07)
!        Unit is amount of aerosol per time of simulation e.g. sulfate[=] kg S
!
!   # : Field   : Description                 : Units     : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) TMS-COND : Condensation rate           : no. or kg : NONE
!  (2) TMS-COAG : Coagulation rate            : no. or kg : NONE
!  (3) TMS-NUCL : Nucleation rate             : no. or kg : NONE
!  (4) TMS-AQOX : Aqueous oxidation rate      : no. or kg : NONE
!  (5) AERO-FIX : Accumulated error fixed     : no. or kg : NONE
!  (6) TMS-SOA  : SOA condensation rate       : no. or kg : NONE
!
!  NOTES:
!  (1 ) Add TMS-SOA diagnostic (win, 9/25/07)
!  28 Jan 2014 - R. Yantosca - Use ARR2D to avoid array temporaries
!*****************************************************************************
!
      IF ( ND60 > 0 ) THEN
         
         !==============================================================
         ! Condensation rate
         !==============================================================

         ! Category name
         CATEGORY = 'TMS-COND'
         DO M = 1,TMAX(60)
            NN = TINDEX(60,M)
            
            SCALEX = 1.e+0_fp

            IF ( ( NN .ge. id_NK1 ) .and. 
     &           ( NN .lt. id_NK1+IBINS )      ) THEN
               UNIT = 'no.'
            ELSE
               UNIT = 'kg'
            ENDIF

            DO L = 1, LD60
               ARR2D(:,L) = AD60_COND(1,:,L,M) / SCALEX
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  1,         JJPAR,     LD60,     IFIRST,
     &                  JFIRST,    LFIRST,    ARR2D(:,1:LD60)  )
         ENDDO

         !==============================================================
         ! Coagulation rate
         !==============================================================

         ! Category name
         CATEGORY = 'TMS-COAG'
         DO M = 1,TMAX(60)
            NN = TINDEX(60,M)
            
            SCALEX = 1.e+0_fp

            IF ( ( NN .ge. id_NK1 ) .and. 
     &           ( NN .lt. id_NK1+IBINS )      ) THEN
               UNIT = 'no.'
            ELSE
               UNIT = 'kg'
            ENDIF

            DO L = 1, LD60
               ARR2D(:,L) = AD60_COAG(1,:,L,M) / SCALEX
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  1,         JJPAR,     LD60,     IFIRST,
     &                  JFIRST,    LFIRST,    ARR2D(:,1:LD60)  )
         ENDDO

         !==============================================================
         ! Nucleation rate
         !==============================================================

         ! Category name
         CATEGORY = 'TMS-NUCL'
         DO M = 1,TMAX(60)
            NN = TINDEX(60,M)
            
            SCALEX = 1.e+0_fp

            IF ( ( NN .ge. id_NK1 ) .and. 
     &           ( NN .lt. id_NK1+IBINS )      ) THEN
               UNIT = 'no.'
            ELSE
               UNIT = 'kg'
            ENDIF

            DO L = 1, LD60
               ARR2D(:,L) = AD60_NUCL(1,:,L,M) / SCALEX
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  1,         JJPAR,     LD60,     IFIRST,
     &                  JFIRST,    LFIRST,    ARR2D(:,1:LD60)  )
         ENDDO

         !==============================================================
         ! Aqueous oxidation rate
         !==============================================================

         ! Category name
         CATEGORY = 'TMS-AQOX'
         DO M = 1,TMAX(60)
            NN = TINDEX(60,M)
            
            SCALEX = 1.e+0_fp

            IF ( ( NN .ge. id_NK1 ) .and. 
     &           ( NN .lt. id_NK1+IBINS )      ) THEN
               UNIT = 'no.'
            ELSE
               UNIT = 'kg'
            ENDIF

            DO L = 1, LD60
               ARR2D(:,L) = AD60_AQOX(1,:,L,M) / SCALEX
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  1,         JJPAR,     LD60,     IFIRST,
     &                  JFIRST,    LFIRST,    ARR2D(:,1:LD60)  )
         ENDDO

         !==============================================================
         ! Error fudging during aerosol microphysics
         !==============================================================

         ! Category name
         CATEGORY = 'AERO-FIX'
         DO M = 1,TMAX(60)
            NN = TINDEX(60,M)
            
            SCALEX = 1.e+0_fp

            IF ( ( NN .ge. id_NK1 ) .and. 
     &           ( NN .lt. id_NK1+IBINS )      ) THEN
               UNIT = 'no.'
            ELSE
               UNIT = 'kg'
            ENDIF

             DO L = 1, LD60
                ARR2D(:,L) = AD60_ERROR(1,:,L,M) / SCALEX
             ENDDO

             CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  1,         JJPAR,     LD60,     IFIRST,
     &                  JFIRST,    LFIRST,    ARR2D(:,1:LD60)  )
          ENDDO

         !==============================================================
         ! SOA Condensation rate  (win, 9/25/07)
         !==============================================================

         ! Category name
         CATEGORY = 'TMS-SOA'
         DO M = 1,TMAX(60)
            NN = TINDEX(60,M)
            
            SCALEX = 1.e+0_fp

            IF ( ( NN .ge. id_NK1 ) .and. 
     &           ( NN .lt. id_NK1+IBINS )      ) THEN
               UNIT = 'no.'
            ELSE
               UNIT = 'kg'
            ENDIF

            DO L = 1, LD60
               ARR2D(:,L) = AD60_SOA(1,:,L,M) / SCALEX
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  1,         JJPAR,     LD60,     IFIRST,
     &                  JFIRST,    LFIRST,    ARR2D(:,1:LD60)  )
         ENDDO

      ENDIF
!
!*****************************************************************************
!  ND61: TOMAS microphysics process rate saved in 3-D (kg/s)
!
!  Remark: for aerosol number, the unit 'kg' is used as a fake kg where
!          1 fake kg = 1 aerosol number.
!          Molecular weight of aerosol number is 1g/mol and   
!
!   # : Field   : Description                 : Units     : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) TOMAS-3D : 3-D microphysics rate       : kg/s      : SCALECHEM
!
!  NOTES: 
!  (1 ) Create ND61 to output nucleation rate in 3-D array diagnostic (win, 6/21/08)
!  (2 ) Change ND61 tracer and unit to be new list with unit of 'cm-3s-1'
!       (win, 10/6/08)
!*****************************************************************************

      IF ( ND61 > 0 ) THEN
         CATEGORY = 'TOMAS-3D'
         
         DO M = 1, TMAX(61)
            NN = TINDEX(61,M)

            SCALEX = SCALECHEM

            UNIT = 'cm-3s-1'

            DO L = 1, LD61
               ARRAY(:,:,L) = AD61(:,:,L,NN) / SCALEX
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD61,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD61)  )
         ENDDO
      ENDIF
#endif
!
!*****************************************************************************
!  ND62: I-J Instantaneous Column Maps for Tracers (molec/cm^2)  
!
!  The unit conversion is as follows:
!
!  First, convert Spc (kg/kg dry) to Spc (kg):
!
!    Spc (kg) = Spc (kg/kg dry) * AD (dry air mass in kg)
!  
!  Next, do this unit conversion:
!
!    Spc (kg) | 6.022e23 molec |   mole   | 1000 g |    1        |   m^2
!    ---------+----------------+----------+--------+-------------+----------
!             |     mole       |  MOLWT g |  kg    | AREA_M2 m^2 | 10^4 cm^2
!
!
!  which is equivalent to
!
!   ( Spc * 6.022e22 ) / ( MOLWT * AREA_M2 )
!
!*****************************************************************************
!
      IF ( ND62 > 0 ) THEN
         CATEGORY = 'INST-MAP'

         ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
         ! (mpayer, 12/6/12)
         Spc => State_Chm%Species
         AD  => State_Met%AD

         DO M = 1, TMAX(62)
            N  = TINDEX(62,M)
            IF ( N > nAdvect ) CYCLE
            NN = N

            ! Get info about this species from the database
            ! NOTE: Assumes 1:1 mapping of tracer and species indexes
            SpcInfo => State_Chm%SpcData(N)%Info

            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Grid box surface area [cm2]
               AREA_M2 = GET_AREA_M2( I, J, 1 )
               ARRAY(I,J,1) = ( SUM( Spc(I,J,:,N)   *
     &                               AD (I,J,:  ) ) * AVO    )
     &                        / ( SpcInfo%emMW_g
     &                            * AREA_M2 * 10.e+0_fp )
            ENDDO
            ENDDO

            ! Write the proper unit string
            IF ( SpcInfo%emMW_g > 12e+0_fp ) THEN
               UNIT = 'molec/cm2'
            ELSE
               UNIT = 'atoms C/cm2'
            ENDIF

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

            ! Free pointer
            SpcInfo => NULL()

         ENDDO

         ! Free pointers
         Spc => NULL()
         AD  => NULL()

      ENDIF
!
!******************************************************************************
!  ND64: Radiative flux (Category: "FJX-FLUX")
!
!   # : Field  : Description                    : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (W)  DIR-X  : Direct flux for bin X          : 1/cm2s     : 1        
!  (W)  DIF-X  : Diffuse flux for bin X         : 1/cm2s     : 1
!  (W)  NET-X  : Net flux for bin X             : 1/cm2s     : 1
!******************************************************************************
!
      IF ( ND64 > 0 ) THEN
         CATEGORY = 'FJX-FLX$'
         UNIT     = 'cm-2s-1'

         !DO M = 1, TMAX(64)
         DO MM = 1,3
            DO MMB = 1,W_
                M = ((MM-1)*W_) + MMB
                N  = TINDEX(64,M)
                NN = N
                SCALEX = SCALECHEM
   
                ! Divide by # of timesteps
                DO L = 1, LD64
                   ARRAY(:,:,L) = AD64(:,:,L,MMB,MM) / SCALEX
                ENDDO

                CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                      HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                      UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                      IIPAR,     JJPAR,     LD64,     IFIRST,     
     &                      JFIRST,    LFIRST,    ARRAY(:,:,1:LD64) )
            ENDDO
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND65: Production/Loss of specified chemical families
!
!   # : Field   : Description                 : Units       : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) PORL-L=$ : Chemical family P-L rates   : molec/cm3/s : SCALECHEM
!
!  NOTES:
!  (1 ) Make sure the units for NSRCX == 6 (single tracer O3) P-L 
!        coincide with those in "chemo3.f".  
!  (2 ) ND65 now uses allocatable array AD65 instead of AIJ. (bmy, 3/16/00)
!  (3 ) Add L(CH3I) to the ND65 diagnostic -- do not take the average 
!        but instead compute the total sum of L(CH3I) (nad, bmy, 3/20/01)
!  (4 ) Add updates for multi-tracer Ox run from amf (bmy, 7/3/01)
!  (5 ) Now account for time in troposphere for full chemistry. It is
!        assumed that LD45 >= LD65 in using CTO3 (phs, 3/6/07)
!  (6 ) Do not use CTO3 anymore, but the new CTO3_24h, which is the 3D
!        tropospheric chemistry counter (phs, 11/17/08)
!******************************************************************************
!
      IF ( ND65 > 0 ) THEN     
         CATEGORY = 'PORL-L=$'

         IF ( Input_Opt%DO_SAVE_PCO ) THEN
            NFAM = Input_Opt%NFAM + 2
         ELSE 
            NFAM = Input_Opt%NFAM
         ENDIF

         ! Loop over ND65 families
         DO N = 1, NFAM

            ! Select proper unit string
            IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
               NN          = N
               UNIT        = 'kg/s'
               
               WHERE( CTO3_24h(:,:,1:LD65) /= 0 )
                  ARRAY(:,:,1:LD65) = AD65(:,:,1:LD65,N) /
     $                                FLOAT( CTO3_24h(:,:,1:LD65) )
               ELSEWHERE
                  ARRAY(:,:,1:LD65) = 0.
               ENDWHERE

            ELSE
               NN     = N
               UNIT   = 'molec/cm3/s'

               WHERE( CTO3_24h(:,:,1:LD65) /= 0 )
                  ARRAY(:,:,1:LD65) = AD65(:,:,1:LD65,N) /
     $                                FLOAT( CTO3_24h(:,:,1:LD65) )
               ELSEWHERE
                  ARRAY(:,:,1:LD65) = 0.
               ENDWHERE
               
            ENDIF

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD65,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD65) )
         ENDDO       
      ENDIF
!
!******************************************************************************
!  ND66: GMAO 3-D fields 
!
!   # : Field  : Description                       : Units   : Scale factor
!  --------------------------------------------------------------------------
!  (1)  UWND   : GMAO Zonal Winds                  : m/s     : SCALE_A3
!  (2)  VWND   : GMAO Meridional Winds             : m/s     : SCALE_A3
!  (3)  TMPU   : GMAO Temperatures                 : K       : SCALE_I3
!  (4)  SPHU   : GMAO Specific Humidity            : g/kg    : SCALE_I3
!  (5)  CLDMAS : GMAO Cloud Mass Flux              : kg/m2/s : SCALE_A3
!  (6)  DTRAIN : GMAO Detrainment flux             : kg/m2/s : SCALE_A3
!  (7)  OMEGA  : GMAO Vertical velocity            : Pa/s    : SCALE_A3
!
!  NOTES:
!  (1) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2) Add CLDMAS to ND66 diagnostic as field #6, but with tracer index
!       #7 (for compatibility with the existing GAMAP).  (rvm, bmy, 9/8/00)
!  (3) For GEOS-4/fvDAS, UWND, VWND, TMPU, SPHU are A-6 fields.  Adjust
!       the scale factors accordingly.  Also delete KZZ. (bmy, 6/23/03)
!  (4) Modified for GEOS-5 and GCAP (bmy, 6/9/05)
!  20 Aug 2010 - R. Yantosca - Modified for MERRA
!  08 Feb 2012 - R. Yantosca - Modified for GEOS-FP
!  03 Feb 2017 - M. Sulprizio- Add OMEGA (Q. Chen)
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!******************************************************************************
!
      IF ( ND66 > 0 ) THEN
         CATEGORY = 'DAO-3D-$'

         SCALE_ND66 = SCALE_A3   ! ND66 is 3-hr time-avg data

         DO M = 1, TMAX(66)
            N  = TINDEX(66,M)
            NN = N 
            
            SELECT CASE ( N )

               ! UWND, VWND
               CASE ( 1,2 )
                  SCALEX = SCALE_ND66
                  UNIT   = 'm/s'

               ! TMPU
               CASE ( 3 )
                  SCALEX = SCALE_I3      ! T is an I3 field
                  UNIT   = 'K'

               ! SPHU
               CASE ( 4 )
                  SCALEX = SCALE_I3      ! SPHU is an I3 field
                  UNIT   = 'g/kg'

               ! CLDMAS, DTRAIN
               CASE( 5, 6 )
                  SCALEX = SCALE_ND66 
                  UNIT   = 'kg/m2/s'

               ! OMEGA
               CASE( 7 )
                  SCALEX = SCALE_ND66
                  UNIT   = 'Pa/s'

               CASE DEFAULT
                  CYCLE
            END SELECT

            ARRAY(:,:,1:LD66) = AD66(:,:,1:LD66,N) / SCALEX

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD66,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD66) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND67: GMAO surface fields
!
!   # : Field  : Description                      : Units    : Scale factor
!  -----------------------------------------------------------------------
!  (1 ) HFLUX  : GMAO Sensible Heat Flux          : W/m2     : SCALE_A1
!  (2 ) RADSWG : GMAO Insolation @ Surface        : W/m2     : SCALE_A1
!  (3 ) PREACC : GMAO Accum Precip @ Surface      : mm/day   : SCALE_A1
!  (4 ) PRECON : GMAO Conv Precip @ Surface       : mm/day   : SCALE_A1
!  (5 ) TS     : GMAO Surface Air Temperature     : K        : SCALE_A1
!  (6 ) RADSWT : GMAO Insolation @ Top of Atm     : W/m2     : SCALE_A1
!  (7 ) USTAR  : GMAO Friction Velocity           : m/s      : SCALE_A1
!  (8 ) Z0     : GMAO Roughness Height            : m        : SCALE_A1
!  (9 ) PBL    : GMAO PBL depth                   : m        : SCALE_A1
!  (10) CLDFRC : GMAO Cloud Fraction              : unitless : SCALE_A1
!  (11) U10M   : GMAO U-wind @ 10 m               : m/s      : SCALE_A1
!  (12) V10M   : GMAO V-wind @ 10 m               : m/s      : SCALE_A1
!  (13) PS-PBL : GMAO Boundary Layer Top Pressure : hPa      : SCALEDYN
!  (14) ALBD   : GMAO Surface Albedo              : unitless : SCALE_A1
!  (15) PHIS   : GMAO Geopotential Heights        : m        : SCALED 
!  (17) TROPP  : GMAO Tropopause pressure         : mb       : SCALE_A1
!  (18) SLP    : GMAO Sea Level pressure          : mb       : SCALE_A1
!  (19) TSKIN  : Ground/sea surface temp.         : hPa      : SCALE_A1
!  (20) PARDF  : Photosyn active diffuse rad.     : W/m2     : SCALE_A1
!  (21) PARDR  : Photosyn active direct  rad.     : W/m2     : SCALE_A1
!  (22) GWET   : Top soil wetness                 : unitless : SCALE_A1
!
!  NOTES:
!  (1 ) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2 ) Now use AD67 allocatable array (bmy, 2/17/00)
!  (3 ) Add TROPP as tracer #17 and SLP as tracer #18 (bmy, 10/11/00)
!  (4 ) Now replace SCALE1 with SCALEDYN (bmy, 3/27/03)
!  (5 ) Added TSKIN, PARDF, PARDR, GWET for GEOS-4 (bmy, 6/23/03)
!  (6 ) Fix SCALEX for ALBEDO: use I6 for GEOS-3 only, and A3 for other
!     models (phs, 9/3/08)
!  (7 ) add EFLUX for output. (lin, ccc 5/29/09)
!  26 Aug 2010 - R. Yantosca - Modified for MERRA
!  03 Nov 2011 - M. Payer    - Fix SCALEX for SLP: use A1 for MERRA, otherwise
!                              use I6
!  08 Feb 2012 - R. Yantosca - Modified for GEOS-5.7.x
!  03 Dec 2013 - R. Yantosca - Change unit of PBL height to meters, this used
!                              to be hPa in GEOS-1, GEOS-STRAT, GEOS-3, which
!                              are no longer supported.           
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!                              and remove CLDTOP (only available in MERRA)
!******************************************************************************
!
      IF ( ND67 > 0 ) THEN
         CATEGORY = 'DAO-FLDS'

         SCALE_ND67 = SCALE_A1    ! ND67 fields are hourly

         ! Binary punch file
         DO M = 1, TMAX(67)
            N  = TINDEX(67,M)
            NN = N 

            SELECT CASE ( N ) 
               CASE ( 1, 2, 6 )
                  SCALEX = SCALE_ND67
                  UNIT   = 'W/m2'
               CASE ( 3, 4 )
                  SCALEX = SCALE_ND67
                  UNIT   = 'mm/day'
               CASE ( 5 )
                  SCALEX = SCALE_ND67
                  UNIT   = 'K'
               CASE ( 7, 11, 12 )
                  SCALEX = SCALE_ND67
                  UNIT   = 'm/s'
               CASE ( 8 )
                  SCALEX = SCALE_ND67
                  UNIT   = 'm'
               CASE ( 9 )
                  SCALEX = SCALE_ND67
                  UNIT   = 'm'
               CASE ( 10 )
                  SCALEX = SCALE_ND67
                  UNIT   = 'unitless'
               CASE ( 13 )
                  SCALEX = SCALEDYN
                  UNIT   = 'hPa'
               CASE ( 14 ) 
                  ! Bug fix: For GEOS-3, ALBEDO is an I-6 field, but
                  ! for GEOS-4, GEOS-5, GCAP, it is an A-3 field.
                  ! (lyj, phs, bmy, 10/7/08)
                  SCALEX = SCALE_ND67
                  UNIT   = 'unitless'
               CASE ( 15 )
                  SCALEX = SCALED
                  UNIT   = 'm'
               CASE ( 16 )
                  SCALEX = SCALE_A3     ! CLDTOPS are 3-hr avg'd
                  UNIT   = 'levels'
               CASE ( 17 ) 
                  SCALEX = SCALE_ND67
                  UNIT   = 'hPa'
               CASE ( 18 ) 
                 SCALEX  = SCALE_A1     ! SLP is hourly
                  UNIT   = 'hPa'
               CASE ( 19 )
                  SCALEX = SCALE_ND67
                  UNIT   = 'K'
               CASE ( 20 )
                  SCALEX = SCALE_ND67
                  UNIT   = 'W/m2'
               CASE ( 21 ) 
                  SCALEX = SCALE_ND67
                  UNIT   = 'W/m2'
               CASE ( 22 ) 
                  SCALEX = SCALE_ND67
                  UNIT   = 'unitless'
               CASE ( 23 )
               ! add EFLUX ( Lin, 05/23/08) 
                  SCALEX = SCALE_ND67
                  UNIT   = 'W/m2'
               CASE DEFAULT
                  CYCLE
            END SELECT
                  
            ARRAY(:,:,1) = AD67(:,:,N) / SCALEX

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND68: Grid box diagnostics
!
!   # : Field   : Description                       : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1) BXHEIGHT : Grid box height                   : m        : SCALEDYN
!  (2) AD       : Air mass in grid box              : kg       : SCALEDYN
!  (3) AVGW     : Mixing ratio of water vapor       : v/v      : SCALEDYN
!  (4) AIRNUMDEN: Dry air number density            : molec/cm3: SCALEDYN
!  (5) T        : Temperature                       : K        : SCALEDYN
!  (6) PMID     : Arithmetic mean pressure          : hPa      : SCALEDYN
!  (7) PEDGE    : Level pressure (bottom edge)      : hPa      : SCALEDYN
!  (8) RH       : Relative humidity                 : %        : SCALEDYN
!
!  NOTES:
!  (1) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2) Now replaced SCALE1 with SCALEDYN (bmy, 2/24/03)
!  (3) Bug fix: replace ND68 with LD68 in call to BPCH2 (swu, bmy, 6/9/05)
!******************************************************************************
!
      IF ( ND68 > 0 ) THEN
         CATEGORY = 'BXHGHT-$'
         UNIT     = ''

         DO M = 1, TMAX(68)
            N  = TINDEX(68,M)
            IF ( N > PD68 ) CYCLE
            NN = N 

            DO L = 1, LD68
               ARRAY(:,:,L) = AD68(:,:,L,N) / SCALEDIAG
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD68,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD68) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND69: Grid Box Surface Areas
!
!   # : Field : Description                       : Units : Scale factor
!  --------------------------------------------------------------------------
!  (1) DXYP   : Surface area of grid box          : m^2   : SCALED = 1.0
!
!  NOTES:
!  (1) Only print DXYP for the first timestep, as it is an invariant field.
!  (2) We don't need to add TRCOFFSET to N.  This is not a CTM tracer.
!  (3) Now use the AD69 dynamically allocatable array (bmy, 2/17/00)
!******************************************************************************
!
      IF ( ND69 > 0 ) THEN 
         CATEGORY = 'DXYP'
         UNIT     = 'm2'

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, 1,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    AD69(:,:,1) )

         ! Set ND69 = 0 so we won't print it out again
         ND69 = 0
      ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: Bob Y. blocked out this section of code which only executes if
!%%% the GTMM Mercury Model is being used.  The bpch restart files are being 
!%%% phased out.  It is unclear anyone is still using the GTMM land Hg model 
!%%% anymore.  If so, then we can investigate a way to add the GTMM diagnostic 
!%%% quantities into the netCDF restart file. (bmy, 4/22/16)
!%%%      !==================================================================
!%%%      ! Special case for mercury simulation. We need to store AD38, AD39,
!%%%      ! AD44 to ensure that we have monthly average in GTMM restart file
!%%%      !==================================================================
!%%%      IF ( Input_Opt%LGTMM ) THEN
!%%%         N = 1
!%%%         NN = GET_WETDEP_IDWETD( N )
!%%%         DO WHILE( .NOT.(IS_Hg2( NN )) )
!%%%            
!%%%            N = N + 1
!%%%            ! Tracer number
!%%%            NN = GET_WETDEP_IDWETD( N )
!%%%            
!%%%         ENDDO
!%%%
!%%%         CALL UPDATE_DEP( N )
!%%%      ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!******************************************************************************
!  ND71: Tracer Mixing Ratio (v/v) at surface, maximum value
!
!   # : Field   : Description                 : Units   : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) IJ-MAX-$ : Peak Tracer mixing ratio    : v/v     : 1        
!
!  NOTES:
!  (1) First written (SDE 2013-11-17)
!******************************************************************************
!
      IF ( ND71 > 0 ) THEN
         CATEGORY = 'IJ-MAX-$'
         UNIT     = ''   

         DO M = 1, TMAX(71)
            N  = TINDEX(71,M)
            IF ( N > nAdvect ) CYCLE
            NN = N
            
            ARRAY(:,:,1) = AD71(:,:,N)/AD71_COUNT

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF

      ! Echo output
      WRITE( 6, '(a)' ) '     - DIAG3: Diagnostics written to bpch!'
!
!******************************************************************************
!  ND72: RRTMG radiation fields
!
!   # : Field  : Description                      : Units    : Scale
!   factor
!  -----------------------------------------------------------------------
!  (1 ) ALLTOASW  : All-sky TOA SW (Total)        : W/m2     : SCALERAD
!  (2 ) ALLSRFSW  : All-sky Surface SW (Total)    : W/m2     : SCALERAD
!  (3 ) ALLTOALW  : All-sky TOA LW (Total)        : W/m2     : SCALERAD
!  (4 ) ALLSRFLW  : All-sky Surface LW (Total)    : W/m2     : SCALERAD
!  (5 ) CLRTOASW  : Clear-sky TOA SW (Total)      : W/m2     : SCALERAD
!  (6 ) CLRSRFSW  : Clear-sky Surface SW (Total)  : W/m2     : SCALERAD
!  (7 ) CLRTOALW  : Clear-sky TOA LW (Total)      : W/m2     : SCALERAD
!  (8 ) CLRSRFLW  : Clear-sky Surface LW (Total)  : W/m2     : SCALERAD
!
!******************************************************************************
!
      !check that RRTMG is actually switched on
      IF ( Input_Opt%LRAD ) THEN
      IF ( ND72 > 0 ) THEN
         CATEGORY = 'RADMAP-$'
         write(6,*) 'Input_Opt%NWVSELECT',Input_Opt%NWVSELECT

         ! ND72 is updated every rad timestep
         SCALEX = SCALERAD

         ! Number of output species plus baseline
         IF ( Input_Opt%LUCX ) THEN
            NSPECOUT=NSPECRADMENU+1
         ELSE
            NSPECOUT=NSPECRADMENU
         ENDIF

         DO M = 1, TMAX(72)
            N  = TINDEX(72,M)
            
            IF ( N > PD72R ) CYCLE
             !only output clear-sky and all-sky if they are switched on
             !or if we're doing the optics
             !but only output 2nd and 3rd sets of optics if they're
             !requested by user in input.geos.rad
             IF (((N.LE.4).AND.(Input_Opt%LSKYRAD(2).EQV..TRUE.)).OR.
     &           ((N.GT.4).AND.(N.LE.8).AND.
     &                         (Input_Opt%LSKYRAD(1).EQV..TRUE.)).OR.
     &           ((N.GT.8).AND.(N.LE.11)).OR.                !1st set of optics
     &           ((N.GT.11).AND.(N.LE.14).AND.
     &                          (Input_Opt%NWVSELECT.GT.1)).OR.   !2nd set
     &           ((N.GT.14).AND.(N.LE.17).AND.
     &                          (Input_Opt%NWVSELECT.GT.2))) THEN !3rd set

             ! Select proper unit string (cf list above)
             IF (N.LE.8) THEN
              UNIT = 'W/m2'
             ELSE
              !AOD, SSA, ASYM
              UNIT = 'UNITLESS'
             ENDIF
               
             ! each case is for different output type, within each
             ! of those are the values for each species in the RAD input
             ! menu. Only output species that are switched on for each
           ! output type selected in the ND72 menu
           DO IS=1,NSPECOUT
            NN=(N-1)*NSPECOUT+IS

            !if output is a flux...
            IF (N.LE.8) THEN
             !always output the baseline flux
             IF (IS.EQ.1) THEN
              ARRAY(:,:,1) = AD72(:,:,NN) / SCALEX
              CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                   HALFPOLAR, CENTER180, CATEGORY, NN,
     &                   UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                   IIPAR,     JJPAR,     1,        IFIRST,
     &                   JFIRST,    LFIRST,    ARRAY(:,:,1) )
             !for other outputs check species has been selected
             ELSE IF (IS.GE.2) THEN
              IF (LSPECRADMENU(IS-1).EQ.1) THEN
               DO I=1,IIPAR
                DO J=1,JJPAR
                 !store as difference in flux from baseline
                 ARRAY(I,J,1) = (AD72(I,J,(N-1)*NSPECOUT+1) - 
     &                           AD72(I,J,NN))/SCALEX
                ENDDO
               ENDDO
             
               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                   HALFPOLAR, CENTER180, CATEGORY, NN,
     &                   UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                   IIPAR,     JJPAR,     1,        IFIRST,
     &                   JFIRST,    LFIRST,    ARRAY(:,:,1) )
              ENDIF
             ENDIF
            ! N > 8 so we are outputting optics
            ! but only for turned on species and only for aerosols
            ! i.e. skip IS=2 and IS=3
            ELSE IF (IS.GE.4) THEN
             IF (LSPECRADMENU(IS-1).EQ.1) THEN
             !output optics
              ARRAY(:,:,1) = AD72(:,:,NN) / SCALEX
              CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                   HALFPOLAR, CENTER180, CATEGORY, NN,
     &                   UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                   IIPAR,     JJPAR,     1,        IFIRST,
     &                   JFIRST,    LFIRST,    ARRAY(:,:,1) )
             ENDIF
            ENDIF !flux vs optics check
           ENDDO
          ENDIF !all-sky, clear-sky check
         ENDDO
      ENDIF
      ENDIF

!******************************************************************************
!  ND73: ISORROPIA pH and H+ (eam, 01/2015)
!   # Category 
!   # : Field     : Description                      : Units     : Scale factor
!  ----------------------------------------------------------------------------
!  (1): pH  : ISORROPIA pH                            : unitless  : SCALECHEM
!  (2): H+  : ISORROPIA H+                            : mol/L     : SCALECHEM
!  (3): H2O : ISORROPIA aerosol water                 : ug/m3 air : SCALECHEM
!
!******************************************************************************
      IF ( ND73 > 0 ) THEN
         CATEGORY  = 'ISOROP-$'
         UNIT      = 'unitless'

         DO M = 1, TMAX(73)
            N = TINDEX(73,M)
            IF ( N > PD73 ) CYCLE
            NN = N

            DO L = 1, LD73
               ARRAY(:,:,L) = AD73(:,:,L,N) / SCALECHEM
            ENDDO

            ! Get correct units:
            IF ( NN == 2 ) THEN
               UNIT = 'mol/L'
            ENDIF
            IF ( NN == 3 ) THEN
               UNIT = 'ug/m3'
            ENDIF

            ! Save to disk
            CALL BPCH2( IU_BPCH,    MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR,  CENTER180, CATEGORY, NN,
     &                  UNIT,       DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,      JJPAR,     LD73,     IFIRST,
     &                  JFIRST,     LFIRST,    ARRAY(:,:,1:LD73) )
         ENDDO
      ENDIF

!******************************************************************************

      ! Cleanup
      DiagnCnt => NULL()
      SpcInfo => NULL()

#endif

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diag2bpch
!
! !DESCRIPTION: Wrapper routine to get diagnostics from HEMCO and write 
!  them to bpch.  This will look up diagnostics 'dName' and write the 
!  corresponding diagnostics array to the bpch output file. 
!\\
!\\
!  NOTE: This is a "bridge" routine intended to provide backwards compatibility 
!  with the existing GEOS-Chem diagnostics.  Eventually we will write all 
!  GEOS-Chem diagnostics to netCDF format but we are not yet at that point.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG2BPCH( air,       HcoState, dname, 
     &                      bcat,      bUnit,    bN,
     &                      dAF,       dZ,       dOPTIONAL,  
     &                      dFACTOR,   ERR,      AreaScal  )
!
! !USES:
!
      USE HCO_STATE_MOD,  ONLY : HCO_State
!
! !INPUT PARAMETERS: 
!
      LOGICAL,          INTENT(IN)    :: air       ! Are we on the root CPU?
      TYPE(HCO_State),  POINTER       :: HcoState  ! HEMCO State object
      CHARACTER(LEN=*), INTENT(IN)    :: dname     ! Diagnostics name
      CHARACTER(LEN=*), INTENT(IN)    :: bCat      ! BPCH category
      CHARACTER(LEN=*), INTENT(IN)    :: bUnit     ! BPCH units
      INTEGER,          INTENT(IN)    :: bN        ! BPCH ID 
      INTEGER,          INTENT(IN)    :: dAF       ! AutoFill of diagnostic
      INTEGER,          INTENT(IN)    :: dZ        ! # of vertical levels
      LOGICAL,          INTENT(IN)    :: dOPTIONAL ! Optional field?
      REAL(fp),         INTENT(IN)    :: dFACTOR   ! Scale factor
      INTEGER,          OPTIONAL      :: AreaScal  ! Area scaling 
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER,          INTENT(INOUT) :: ERR       ! Return code
!
! !REMARKS:
!  (1) Data is multiplied by factor dFACTOR before writing to disk.
!  (2) dAF determines the autofill flag used for the diagnostics lookup.
!  (3) dZ are the number of vertical levels. Set to -1 for 2D arrays, 
!       otherwise a 3D array with exactly that many levels is written.
!  (4) dOptional determines whether or not this is an optional diagnostics. 
!       For optional diagnostics, zeros are written out if no corresponding 
!       HEMCO diagnostics could be found (because it is undefined or because 
!       they have never been updated). For non-optional diagnostics, the 
!       routine stops with an error if no HEMCO diagnostics is found.
!  (5) AScal can be used to multiply (1) or divide (-1) the data by the grid 
!       area.
! 
! !REVISION HISTORY: 
!  06 Aug 2014 - C. Keller   - Initial version
!  13 Aug 2014 - R. Yantosca - Added ProTeX headers; cosmetic changes
!  13 Aug 2014 - R. Yantosca - Bug fix, # of vertical levels in call to BPCH2
!                              was 1 but should have been LEVS.  Now fixed.
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      CHARACTER(LEN=255), PARAMETER :: LC = 'DIAG2BPCH (diag3.F)'
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER                       :: II, FLG, LEVS, ASCL
      CHARACTER(LEN=255)            :: MG

      ! Pointers
      TYPE(DiagnCont),    POINTER   :: DgnCont

      !=================================================================
      ! DIAG2BPCH begins here!
      !=================================================================
#if defined( BPCH_DIAG )

      ! Initialize
      ARRAY    = 0e0
      DgnCont => NULL()

      ! Set default for area scaling
      ASCL = 0
      IF ( PRESENT(AreaScal) ) ASCL = AreaScal
 
      ! Get diagnostics
      DgnCont => NULL()
      CALL Diagn_Get( air,   HcoState, .FALSE., DgnCont, 
     &                FLG,   ERR,      cName=TRIM(dname), 
     &                AutoFill=dAF, COL=HcoState%Diagn%HcoDiagnIDManual)

      ! If diagnostic not found, write error message
      IF ( ERR /= HCO_SUCCESS ) THEN
         MG = 'Cannot find diagnostics ' // TRIM(dname)
         CALL ERROR_STOP ( MG, LC ) 
      ENDIF

      ! Vertical levels
      levs = MAX(dz,1)

      ! Save into ARRAY. Apply scale factor if required
      IF ( FLG == HCO_SUCCESS ) THEN
         IF ( dz > 0 ) THEN
            ARRAY(:,:,1:levs) = 
     &           DgnCont%Arr3D%Val(:,:,1:levs) * dFACTOR
         ELSE 
            ARRAY(:,:,1) = DgnCont%Arr2D%Val(:,:) * dFACTOR
         ENDIF
      ELSE
         MG = 'No diagnostics returned: ' // TRIM(dname)
         IF ( dOptional ) THEN
            MG = TRIM(MG) // ' - will write zeros!'
            CALL HCO_WARNING ( MG, ERR, THISLOC=LC )
            ARRAY(:,:,1:levs) = 0.0
         ELSE
            CALL ERROR_STOP ( MG, LC )
         ENDIF
      ENDIF

      ! Eventually scale by area
      IF ( ASCL == 1 ) THEN
         DO II=1,levs
            ARRAY(:,:,II) = 
     &           ARRAY(:,:,II) * HcoState%Grid%AREA_M2%Val(:,:)
         ENDDO
      ELSEIF ( ASCL == -1 ) THEN
         DO II=1,levs
            ARRAY(:,:,II) = 
     &           ARRAY(:,:,II) / HcoState%Grid%AREA_M2%Val(:,:)
         ENDDO
      ENDIF
 
      ! Write to bpch file
      CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &            HALFPOLAR, CENTER180, bCat,     bN,
     &            bUnit,     DIAGb,     DIAGe,    RESERVED,   
     &            IIPAR,     JJPAR,     LEVS,     IFIRST,     
     &            JFIRST,    LFIRST,    ARRAY(:,:,1:levs) )

      ! Return w/ success
      ERR = HCO_SUCCESS

      ! Nullify pointer
      DgnCont => NULL()

#endif

      END SUBROUTINE DIAG2BPCH
!EOC
      END SUBROUTINE DIAG3    
