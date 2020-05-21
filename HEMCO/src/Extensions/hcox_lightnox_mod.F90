!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_lightnox_mod.F90
!
! !DESCRIPTION: Module HCOX\_LightNOx\_Mod contains routines to
!  compute NO lightning emissions, according to the GEOS-Chem lightning
!  algorithms.
!\\
!\\
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities. In particular, the LIS-OTD local redistribution factors are
! now read through the HEMCO framework, and the corresponding netCDF
! input file is specified in the HEMCO configuration file. The table of
! cumulative distribution functions used to vertically distribute lightning
! NOx emissions is specified in the extension switch section of the
! configuration file.
!\\
!\\
! References:
! \begin{itemize}
! \item Murray, L. T., Jacob, D. J., Logan, J. A., Hudman, R. C., and
!       Koshak, W. J.: \emph{Optimized regional and interannual variability
!       of lightning in a global chemical transport model con- strained
!       by LIS/OTD satellite data}, \underline{J. Geophys. Res.},
!       Atmospheres, 117, 2012.
! \item Ott, L. E., K. E. Pickering, G. L. Stenchikov, D. J. Allen,
!       A. J. DeCaria, B. Ridley, R.-F. Lin, S. Lang, and W.-K. Tao,
!       \emph{Production of lightning NOx and its vertical distribution
!       calculated  from three-dimensional cloud-scale chemical transport
!       model simulations}, \underline{J. Geophys. Res.}, 115, D04301, 2010.
! \end{itemize}
!
! !INTERFACE:
!
MODULE HCOX_LightNOx_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE HCOX_TOOLS_MOD
  USE HCO_State_Mod,  ONLY : HCO_State
  USE HCOX_State_MOD, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOX_LightNOX_Run
  PUBLIC  :: HCOX_LightNOX_Final
  PUBLIC  :: HCOX_LightNOX_Init
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: LIGHTNOX
  PRIVATE :: LIGHTDIST
!
! !PUBLIC DATA MEMBERS:
!
!
! !REMARKS:
!  %%% NOTE: MFLUX and PRECON methods are now deprecated (ltm, bmy, 7/9/09)
!                                                                             .
!  References:
!  ============================================================================
!  (1 ) Price & Rind (1992), JGR, vol. 97, 9919-9933.
!  (2 ) Price & Rind (1994), M. Weather Rev, vol. 122, 1930-1939.
!  (3 ) Allen & Pickering (2002), JGR, 107, D23, 4711, doi:10.1029/2002JD002066
!  (4 ) Hudman et al (2007), JGR, 112, D12S05, doi:10.1029/2006JD007912
!  (5 ) Sauvage et al, 2007, ACP,
!        http://www.atmos-chem-phys.net/7/815/2007/acp-7-815-2007.pdf
!  (6 ) Ott et al., (2010), JGR
!  (7 ) Allen et al., (2010), JGR
!  (8 ) Murray et al., (2011), in prep.
!
! !REVISION HISTORY:
!  14 Apr 2004 - L. Murray, R. Hudman - Initial version
!  (1 ) Based on "lightnox_nox_mod.f", but updated for near-land formulation
!        and for CTH, MFLUX, PRECON parameterizations (ltm, bmy, 5/10/06)
!  (2 ) Now move computation of IC/CG flash ratio out of routines FLASHES_CTH,
!        FLASHES_MFLUX, FLASHES_PRECON, and into routine GET_IC_CG_RATIO.
!        Added a fix in LIGHTDIST for pathological grid boxes.  Set E_IC_CG=1
!        according to Allen & Pickering [2002].  Rename OTDSCALE array to
!        OTD_REG_REDIST, and also add OTD_LOC_REDIST array.  Now scale
!        lightnox to 6 Tg N/yr for both 2x25 and 4x5.  Rename routine
!        GET_OTD_LIS_REDIST to GET_REGIONAL_REDIST.  Add similar routine
!        GET_LOCAL_REDIST.  Removed GET_OTD_LOCp AL_REDIST.  Bug fix: divide
!        A_M2 by 1d6 to get A_KM2. (rch, ltm, bmy, 2/22/07)
!  (3 ) Rewritten for separate treatment of LNOx emissions at tropics &
!        midlatitudes, based on Hudman et al 2007.  Removed obsolete
!        variable E_IC_CG. (rch, ltm, bmy, 3/27/07)
!  (4 ) Changes implemented in this version (ltm, bmy, 10/3/07)
!        * Revert to not classifying near-land as land
!        * Eliminate NOx emisisons per path length entirely
!        * Scale tropics to 260 mol/fl constraint from Randall Martin's
!           4.4 Tg and OTD-LIS avg ann flash rate
!        * Remove top-down scaling (remove the three functions)
!        * Allow option of mid-level scaling to match global avg ann flash
!           rate between G-C and OTD-LIS 11-year climatology (new function)
!        * Local Redist now a la Murray et al, 2007 in preparation (monthly)
!        * Replace GEMISNOX (from CMN_NOX) with module variable EMIS_LI_NOx
!  (5 ) Added MFLUX, PRECON redistribution options (ltm, bmy, 11/29/07)
!  (6 ) Updated OTD/LIS scaling for GEOS-5 to get more realistic totals
!        (ltm, bmy, 2/20/08)
!  (7 ) Now add the proper scale factors for the GEOS-5 0.5 x 0.666 grid
!        and the GEOS-3 1x1 nested N. America grid in routine
!        GET_OTD_LIS_SCALE. (yxw, dan, ltm, bmy, 11/14/08)
!  (8 ) Added quick fix for GEOS-5 reprocessed met fields (ltm, bmy, 2/18/09)
!  (9 ) Added quick fix for GEOS-5 years 2004, 2005, 2008 (ltm, bmy, 4/29/09)
!  (10) Updated OTD/LIS scaling for GEOS-5 reprocessed data (ltm, bmy, 7/10/09)
!  (11) Updated for GEOS-4 1 x 1.25 grid (lok, ltm, bmy, 1/13/10)
!  (12) Reprocessed for CLDTOPS calculation error; Updated Ott vertical
!        profiles; Removal of depreciated options, e.g., MFLUX and PRECON;
!        GEOS5 5.1.0 vs. 5.2.0 special treatment; MERRA; Other changes.
!        Please see PDF on wiki page for full description of lightnox
!        changes to v9-01-01. (ltm, 1/25/11)
!  13 Aug 2010 - R. Yantosca - Add modifications for MERRA
!  10 Nov 2010 - L. Murray   - Updated OTD/LIS local scaling for MERRA 4x5
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  02 Feb 2012 - R. Yantosca - Added modifications for GEOS-5.7.x met fields
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  22 Oct 2013 - C. Keller   - Now a HEMCO extension.
!  22 Jul 2014 - R. Yantosca - Now hardwire the Lesley Ott et al CDF's in
!                              lightning_cdf_mod.F90.  This avoids having to
!                              read an ASCII input in the ESMF environment.
!  13 Jan 2015 - L. Murray   - Add most recent lightning updates to HEMCO version
!  26 Feb 2015 - R. Yantosca - Restore reading the lightning CDF's from an
!                              ASCII file into the PROFILE array.  This helps
!                              to reduce compilation time.
!  31 Jul 2015 - C. Keller   - Added option to define scalar/gridded scale
!                              factors via HEMCO configuration file.
!  14 Oct 2016 - C. Keller   - Now use HCO_EvalFld instead of HCO_GetPtr.
!  02 Dec 2016 - M. Sulprizio- Update WEST_NS_DIV from 23d0 to 35d0 (K. Travis)
!  16 Feb 2017 - L. Murray   - Updated BETA factors for all GEOS-FP/MERRA-2
!                              products fields available by v11-01 release
!                              (through Dec. 2016), and latest version of
!                              LIS/OTD satellite climatology.
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!  17 Oct 2017 - C. Keller   - Add option to use GEOS-5 lightning flash rate
!                              (LFR). Autoselection of flash rate scale factor.
!  16 Jan 2019 - L. Murray   - Lightning flash densities and convective depths
!                              for vertical distribution are now calculated offline
!                              at the native model resolution and prescribed in
!                              HEMCO_Config.rc as the other model meteorology.
!                              Model flash climatology remains constrained to the
!                              LIS/OTD climatology, as described by Murray et al.
!                              (2012). Removed obsolete subroutines and code.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  INTEGER, PARAMETER            :: NLTYPE        = 4
  INTEGER, PARAMETER            :: NNLIGHT       = 3200
  REAL*8,  PARAMETER            :: RFLASH_MIDLAT = 3.011d26   ! 500 mol/flash applied in the N extratropics
  REAL*8,  PARAMETER            :: RFLASH_TROPIC = 1.566d26   ! 260 mol/flash applied in the tropics / S extratropics
  REAL*8,  PARAMETER            :: EAST_WEST_DIV = -30d0
  REAL*8,  PARAMETER            :: WEST_NS_DIV   =  35d0
  REAL*8,  PARAMETER            :: EAST_NS_DIV   =  35d0
!
! !PRIVATE TYPES:
!
  ! Scalars
  TYPE :: MyInst
   INTEGER                       :: Instance
   INTEGER                       :: IDTNO     ! NO tracer ID
   INTEGER                       :: ExtNr     ! HEMCO Extension ID
   LOGICAL                       :: LCNVFRC   ! Use convective fractions?
   LOGICAL                       :: LLFR      ! Use GEOS-5 flash rates

   ! Arrays
   REAL(dp), POINTER             :: PROFILE(:,:)
   REAL(hp), POINTER             :: SLBASE(:,:,:)
   REAL(sp), POINTER             :: FLASH_DENS_TOT(:,:)
   REAL(sp), POINTER             :: FLASH_DENS_IC(:,:)
   REAL(sp), POINTER             :: FLASH_DENS_CG(:,:)
   REAL(sp), POINTER             :: CONV_DEPTH(:,:)

   ! Overall scale factor to be applied to lightning NOx emissions. Must
   ! be defined in the HEMCO configuration file as extension attribute
   ! 'Scaling_NO'.
   ! SpcScalFldNme is the name of the gridded scale factor. Must be provided
   ! in the HEMCO configuration file as extension attribute 'ScaleField_NO'.
   REAL(sp), ALLOCATABLE          :: SpcScalVal(:)
   CHARACTER(LEN=61), ALLOCATABLE :: SpcScalFldNme(:)

   TYPE(MyInst), POINTER           :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to all instances
  TYPE(MyInst), POINTER            :: AllInst => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_LightNOx_Run
!
! !DESCRIPTION: Subroutine HCOX\_LIGHTNOX\_RUN is the driver routine
! to calculate lightning NOx emissions and return them to the HEMCO
! driver routine.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_LightNOx_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_FluxArr_Mod,  ONLY : HCO_EmisAdd
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER        :: ExtState   ! Module options
    TYPE(HCO_State), POINTER        :: HcoState   ! HEMCO options
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  09 Oct 1997 - R. Yantosca - Initial version
!  (1 ) Remove IOFF, JOFF from the argument list.  Also remove references
!        to header files "CMN_O3" and "comtrid.h" (bmy, 3/16/00)
!  (2 ) Now use allocatable array for ND32 diagnostic (bmy, 3/16/00)
!  (3 ) Now reference BXHEIGHT from "dao_mod.f".  Updated comments, cosmetic
!        changes.  Replace LCONVM with the parameter LLCONVM. (bmy, 9/18/02)
!  (4 ) Removed obsolete reference to "CMN".  Now bundled into
!        "lightnox_mod.f" (bmy, 4/14/04)
!  (5 ) Renamed from EMLIGHTNOX_NL to EMLIGHTNOX.  Now replace GEMISNOX
!        (from CMN_NOX) with module variable EMIS_LI_NOx. (ltm, bmy, 10/3/07)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept State_Chm
!  22 Oct 2013 - C. Keller   - Now a HEMCO extension.
!  07 Oct 2013 - C. Keller   - Now allow OTD-LIS scale factor to be set
!                              externally. Check for transition to Sep 2008.
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(MyInst), POINTER :: Inst
    INTEGER               :: Yr, Mt
    LOGICAL               :: FOUND
    CHARACTER(LEN=255)    :: MSG

    !=================================================================
    ! HCOX_LIGHTNOX_RUN begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_LightNOx_Run (hcox_lightnox_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return if extension disabled
    IF ( ExtState%LightNOx <= 0 ) THEN
       CALL HCO_LEAVE( HcoState%Config%Err,RC )
       RETURN
    ENDIF

    ! Get pointer to this instance. Varible Inst contains all module
    ! variables for the current instance. The instance number is
    ! ExtState%<yourname>.
    Inst => NULL()
    CALL InstGet ( ExtState%LightNOx, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find lightning NOx instance Nr. ', ExtState%LightNOx
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
    ENDIF

    ! Update lightnox NOx emissions (fill SLBASE)
    CALL LIGHTNOX( HcoState, ExtState, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=================================================================
    ! Pass to HEMCO State and update diagnostics
    !=================================================================
    IF ( Inst%IDTNO > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%SLBASE, Inst%IDTNO, &
                         RC, ExtNr=Inst%ExtNr)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: SLBASE', RC )
          RETURN
       ENDIF

    ENDIF

    ! Return w/ success
    Inst => NULL()
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_LightNOx_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LightNOx
!
! !DESCRIPTION: Subroutine LIGHTNOX uses Price \& Rind's formulation for
!  computing NOx emission from lightning (with various updates).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LightNOx( HcoState, ExtState, Inst, RC )
!
! !USES:
!
    USE HCO_Calc_Mod,     ONLY : HCO_EvalFld
    USE HCO_EmisList_Mod, ONLY : HCO_GetPtr
    USE HCO_GeoTools_Mod, ONLY : HCO_LANDTYPE
    USE HCO_Clock_Mod,    ONLY : HcoClock_Get
    USE HCO_Clock_Mod,    ONLY : HcoClock_First
    USE HCO_ExtList_Mod,  ONLY : GetExtOpt
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState  ! Output obj
    TYPE(Ext_State), POINTER        :: ExtState    ! Module options
    TYPE(MyInst   ), POINTER        :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  10 May 2006 - L. Murray - Initial version
!  (1 ) Now recompute the cold cloud thickness according to updated formula
!        from Lee Murray.  Rearranged argument lists to routines FLASHES_CTH,
!        FLASHES_MFLUX, FLASHES_PRECON.  Now call READ_REGIONAL_REDIST and
!        READ_LOCAL_REDIST. Updated comments accordingly.  Now apply
!        FLASH_SCALE to scale the total lightnox NOx to 6 Tg N/yr.  Now apply
!        OTD/LIS regional or local redistribution (cf. B. Sauvage) to the ND56
!        diagnostic. lightnox redistribution to the ND56 diag.  Renamed
!        REGSCALE variable to REDIST.  Bug fix: divide A_M2 by 1d6 to get
!        A_KM2. (rch, ltm, bmy, 2/14/07)
!  (2 ) Rewritten for separate treatment of LNOx emissions at tropics &
!        midlatitudes (rch, ltm, bmy, 3/27/07)
!  (3 ) Remove path-length algorithm.  Renamed from LIGHTNOX_NL to LIGHTNOX.
!        Other improvements. (ltm, bmy, 9/24/07)
!  (4 ) Remove depreciated options; Update to new Ott et al vertical profiles;
!        Reprocessed for bug in CLDTOPS calculation. See PDF on wiki for
!        full description of changes for v9-01-01. (ltm, bmy, 1/25,11)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  22 Oct 2013 - C. Keller   - Now a HEMCO extension.
!  06 Oct 2014 - C. Keller   - Now calculate pressure centers from edges.
!  16 Jan 2015 - R. Yantosca - Bug fix: TmpScale should be REAL(dp)
!  11 Mar 2015 - C. Keller   - Now determine LTOP from buoyancy for grid boxes
!                              where convection is explicitly resolved. For now,
!                              this will only work in an ESMF environment.
!  31 Jul 2015 - C. Keller   - Take into account scalar/gridded scale factors
!                              defined in HEMCO configuration file.
!  03 Mar 2016 - C. Keller   - Use buoyancy in combination with convective
!                              fraction CNV_FRC (ESMF only).
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  30 Aug 2018 - C. Keller   - Use diagnostic names for special diagnostics.
!                              This makes interface with GEOS easier.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER           :: I, J, L
    INTEGER           :: LTOP
    INTEGER           :: MONTH
    INTEGER           :: MTYPE
    REAL*8            :: A_M2
    REAL*8            :: A_KM2
    REAL*8            :: H0
    REAL*8            :: IC_CG_RATIO
    REAL*8            :: RATE
    REAL*8            :: RATE_SAVE
    REAL*8            :: TOTAL
    REAL*8            :: TOTAL_CG
    REAL*8            :: TOTAL_IC
    REAL*8            :: X
    REAL*8            :: YMID
    REAL*8            :: XMID
    REAL*8            :: VERTPROF(HcoState%NZ)
    INTEGER           :: LMAX
    INTEGER           :: LNDTYPE
    INTEGER           :: SFCTYPE
    REAL(hp)          :: TROPP
    REAL(dp)          :: TmpScale

    !=================================================================
    ! LIGHTNOX begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'LightNOx (hcox_lightnox_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Reset arrays
    Inst%SLBASE         = 0.0_hp
    Inst%FLASH_DENS_TOT = 0.0_sp
    Inst%FLASH_DENS_IC  = 0.0_sp
    Inst%FLASH_DENS_CG  = 0.0_sp
    Inst%CONV_DEPTH     = 0.0_sp

    ! LMAX: the highest L-level to look for lightning NOx (usually LLPAR-1)
    LMAX   = HcoState%NZ - 1

    ! Get current month (to be passed to LIGHTDIST)
    CALL HcoClock_Get( HcoState%Clock, cMM=MONTH, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=================================================================
    ! Compute lightning NOx emissions for each (I,J) column
    !=================================================================

!$OMP PARALLEL DO                                                     &
!$OMP DEFAULT( SHARED )                                               &
!$OMP PRIVATE( I,         J,           L,        A_M2,        A_KM2  ) &
!$OMP PRIVATE( YMID,      XMID,        LTOP,     MTYPE               ) &
!$OMP PRIVATE( LNDTYPE,   SFCTYPE,     TROPP                         ) &
!$OMP PRIVATE( RATE,      RATE_SAVE,   H0,       IC_CG_RATIO         ) &
!$OMP PRIVATE( TOTAL,     TOTAL_IC,    TOTAL_CG, VERTPROF,    X      ) &
!$OMP SCHEDULE( DYNAMIC )

    ! Loop over surface boxes
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! Grid box surface areas in [m2] and [km2]
       A_M2     = HcoState%Grid%AREA_M2%Val( I, J )
       A_KM2    = A_M2 / 1d6

       ! Grid box latitude and longitude [degrees]
       YMID     = HcoState%Grid%YMID%Val( I, J )
       XMID     = HcoState%Grid%XMID%Val( I, J )

       ! Make sure xmid is between -180 and +180
       IF ( XMID >= 180.0d0 ) XMID = XMID - 360.0d0

       ! Get surface type. Note that these types are different than
       ! the types used elsewhere: 0 = land, 1=water, 2=ice!
       LNDTYPE = HCO_LANDTYPE( ExtState%WLI%Arr%Val(I,J),  &
                               ExtState%ALBD%Arr%Val(I,J) )

       ! Adjusted SFCTYPE variable for this module:
       IF ( LNDTYPE == 2 ) THEN
          SFCTYPE = 2    ! Ice
       ELSEIF ( LNDTYPE == 1 ) THEN
          SFCTYPE = 0    ! Land
       ELSE
          SFCTYPE = 1    ! Ocean (default)
       ENDIF

       ! Tropopause pressure. Convert to Pa
       TROPP = ExtState%TROPP%Arr%Val(I,J) !* 100.0_hp

       !===========================================================
       ! Initialize
       !===========================================================
       RATE          = 0.0
       RATE_SAVE     = 0.0
       H0            = 0.0
       IC_CG_RATIO   = 1.0

       TOTAL         = 0d0
       TOTAL_IC      = 0d0
       TOTAL_CG      = 0d0

       !===========================================================
       ! (1) Get flash density [#/km2/s] from meteorology
       !===========================================================

       !-----------------------------------------------------------
       ! (1a) Prescribed in HEMCO_Config.rc
       !-----------------------------------------------------------
       RATE        = ExtState%FLASH_DENS%Arr%Val(I,J)

       !-----------------------------------------------------------
       ! (1b) From GEOS-5
       !-----------------------------------------------------------
       IF ( Inst%LLFR                                  &
            .AND. ASSOCIATED( ExtState%LFR%Arr%Val   ) &
            .AND. ASSOCIATED( ExtState%BYNCY%Arr%Val )  ) THEN
          RATE     = ExtState%LFR%Arr%Val(I,J)
       ENDIF

       ! Error check: do not continue if flash rate is zero
       IF ( RATE <= 0.0 ) CYCLE

       !===========================================================
       ! (2) Get depth of convection [m] and find associated LTOP
       !===========================================================

       !-----------------------------------------------------------
       ! (2a) Prescribed in HEMCO_Config.rc
       !-----------------------------------------------------------
       H0          = ExtState%CONV_DEPTH%Arr%Val(I,J)
       LTOP = 1
       DO L = 1, HcoState%NZ
          IF ( SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:L)) > H0 ) THEN
             LTOP = L
             EXIT
          ENDIF
       ENDDO
       ! Reset H0 to be the height of that layer in the model,
       ! to avoid negative values in the partitioning
       H0          = SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:LTOP))

       !-----------------------------------------------------------
       ! (2b) From GEOS-5
       !-----------------------------------------------------------
       IF ( Inst%LLFR                                  &
            .AND. ASSOCIATED( ExtState%LFR%Arr%Val   ) &
            .AND. ASSOCIATED( ExtState%BYNCY%Arr%Val )  ) THEN

          ! Set LTOP to top of buoyancy
          DO L = HcoState%NZ, 1, -1
             IF ( ExtState%BYNCY%Arr%Val(I,J,L) >= 0.0_sp ) THEN
                LTOP = L + 1
                EXIT
             ENDIF
          ENDDO
          !LTOP = MAX( LTOP, LMAX )
          ! H0 is the convective cloud top height [m].  This is the
          ! distance from the surface to the top edge of box (I,J,LTOP).
          H0 = SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:LTOP))
       ENDIF

       ! Save out convective cloud depth
       Inst%CONV_DEPTH(I,J) = LTOP

       !===========================================================
       ! (3) Compute ratio of CG vs total flashes
       !===========================================================

       ! Ratio of cloud-to-ground flashes to total # of flashes
       X    = 1d0 / ( 1d0 + IC_CG_RATIO )

       !-----------------------------------------------------------
       ! Store flash rates [flashes/km2/min]
       !-----------------------------------------------------------
       IF ( RATE > 0d0 ) THEN

          ! Flashes per km2 per minute
          RATE_SAVE   = RATE / 360d0

          ! Store total, IC, and CG flash rates
          Inst%FLASH_DENS_TOT(I,J) = RATE_SAVE
          Inst%FLASH_DENS_IC(I,J)  = RATE_SAVE * X
          Inst%FLASH_DENS_CG(I,J)  = H0 * 1d-3

       ENDIF

       !===========================================================
       ! (4) Compute LNOx yield for IC and CG flashes (molec/km2/s)
       !===========================================================

       ! Compute LNOx emissions for tropics or midlats
       IF ( XMID > EAST_WEST_DIV ) THEN

          !--------------------------------------------------------
          ! (4a) We are in EURASIA
          !--------------------------------------------------------
          IF ( YMID > EAST_NS_DIV ) THEN

             ! Eurasian Mid-Latitudes
             TOTAL_IC = RFLASH_MIDLAT * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_MIDLAT * RATE * X

          ELSE

             ! Eurasian Tropics
             TOTAL_IC = RFLASH_TROPIC * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_TROPIC * RATE * X

          ENDIF

       ELSE

          !--------------------------------------------------------
          ! (4b) We are in the AMERICAS
          !--------------------------------------------------------
          IF ( YMID > WEST_NS_DIV ) THEN

             ! American Mid-Latitudes
             TOTAL_IC = RFLASH_MIDLAT * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_MIDLAT * RATE * X

          ELSE

             ! American Tropics
             TOTAL_IC = RFLASH_TROPIC * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_TROPIC * RATE * X

          ENDIF
       ENDIF

       !===========================================================
       ! (5) Compute column total lightning NOx yield (kg/km2/s)
       !===========================================================

       ! Sum of IC + CG
       TOTAL = TOTAL_IC + TOTAL_CG

       ! Convert from molec km-2 s-1 to kg(NO) m-2 s-1
       TOTAL = TOTAL * ( HcoState%Spc(Inst%IDTNO)%EmMW_g / 1000.0_hp ) / &
                         HcoState%Phys%Avgdr / 1000000.0_hp

       !===========================================================
       ! (6) Distribute column LNOx vertically from surface to LTOP
       !===========================================================

       !-----------------------------------------------------------
       ! LIGHTDIST computes the lightning NOx distribution from
       ! the ground to the convective cloud top using cumulative
       ! distribution functions for ocean flashes, tropical land
       ! flashes, and non-tropical land flashes, as specified by
       ! Lesley Ott [JGR, 2010]
       !-----------------------------------------------------------

       ! If there's lightning NOx w/in the column ...
       IF ( TOTAL > 0d0 ) THEN

          ! Partition the column total NOx [kg/m2/s] from lightning
          ! into the vertical using Ott et al. [2010] PDF functions
          CALL LIGHTDIST( I, J, LTOP, H0, YMID, TOTAL, VERTPROF, &
                          ExtState, HcoState, SFCTYPE, MONTH, MTYPE, Inst )

          ! Add vertically partitioned NOx into SLBASE array
          DO L = 1, HcoState%NZ
             Inst%SLBASE(I,J,L) = VERTPROF(L)

             ! No lightning NOx emissions in the stratosphere (cdh, 4/25/2013)
             IF ( HcoState%Grid%PEDGE%Val(I,J,L) < TROPP ) THEN
                Inst%SLBASE(I,J,L) = 0.0_hp
             ENDIF

          ENDDO
       ENDIF

    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    !-----------------------------------------------------------------
    ! Eventually add scale factors
    !-----------------------------------------------------------------

    ! Eventually apply species specific scale factor
    IF ( Inst%SpcScalVal(1) /= 1.0_sp ) THEN
       Inst%SLBASE = Inst%SLBASE * Inst%SpcScalVal(1)
    ENDIF

    ! Eventually apply spatiotemporal scale factors
    CALL HCOX_SCALE( HcoState, Inst%SLBASE, TRIM(Inst%SpcScalFldNme(1)), RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE LightNOx
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LightDist
!
! !DESCRIPTION: Subroutine LightDist reads in the CDF used to partition the
!  column lightning NOx into the GEOS-Chem vertical layers.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LightDist( I, J, LTOP, H0, XLAT, TOTAL, VERTPROF, &
                        ExtState, HcoState, SFCTYPE, MONTH, MTYPE, Inst )
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN)  :: I          ! Longitude index
    INTEGER,         INTENT(IN)  :: J          ! Latitude index
    INTEGER,         INTENT(IN)  :: LTOP       ! Level of conv cloud top
    REAL*8,          INTENT(IN)  :: H0         ! Conv cloud top height [m]
    REAL*8,          INTENT(IN)  :: XLAT       ! Latitude value [degrees]
    REAL*8,          INTENT(IN)  :: TOTAL      ! Column Total # of LNOx molec
    TYPE(Ext_State), POINTER     :: ExtState   ! Module options
    TYPE(HCO_State), POINTER     :: HcoState   ! Hemco state object
    INTEGER,         INTENT(IN)  :: SFCTYPE    ! Surface type
    INTEGER,         INTENT(IN)  :: MONTH      ! Current month
    TYPE(MyInst),    POINTER     :: Inst       ! Hemco state object
!
! !OUTPUT PARAMETERS:
!
    REAL*8,          INTENT(OUT) :: VERTPROF(HcoState%NZ) ! Vertical profile
    INTEGER,         INTENT(OUT) :: MTYPE                 ! lightning type
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Pickering et al., JGR 103, 31,203 - 31,316, 1998.
!  (2 ) Ott et al., JGR, 2010
!  (3 ) Allen et al., JGR, 2010
!
! !REVISION HISTORY:
!  18 Sep 2002 - M. Evans - Initial version (based on Yuhang Wang's code)
!  (1 ) Use functions IS_LAND and IS_WATER to determine if the given grid
!        box is over land or water.  These functions work for all DAO met
!        field data sets. (bmy, 4/2/02)
!  (2 ) Renamed M2 to LTOP and THEIGHT to H0 for consistency w/ variable names
!        w/in "lightning.f".  Now read the "light_dist.dat.geos3" file for
!        GEOS-3 directly from the DATA_DIR/lightning_NOx_200203/ subdirectory.
!        Now read the "light_dist.dat" file for GEOS-1, GEOS-STRAT directly
!        from the DATA_DIR/lightning_NOx_200203/ subdirectory.  Added
!        descriptive comment header.  Now trap I/O errors across all
!        platforms with subroutine "ioerror.f".  Updated comments, cosmetic
!        changes.  Redimension FRAC(NNLIGHT) to FRAC(LLPAR). (bmy, 4/2/02)
!  (3 ) Deleted obsolete code from April 2002.  Now reference IU_FILE and
!        IOERROR from "file_mod.f".  Now use IU_FILE instead of IUNIT as the
!        file unit number. (bmy, 6/27/02)
!  (4 ) Now reference BXHEIGHT from "dao_mod.f" (bmy, 9/18/02)
!  (5 ) Bug fix: add GEOS_4 to the #if block (bmy, 3/4/04)
!  (6 ) Now bundled into "lightning_mod.f".  CDF's are now read w/in
!        routine INIT_LIGHTNOX to allow parallelization (bmy, 4/14/04)
!  (7 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (8 ) Now uses near-land formulation (ltm, bmy, 5/10/06)
!  (9 ) Added extra safety check for pathological boxes (bmy, 12/11/06)
!  (10) Remove the near-land formulation, except for PRECON (ltm, bmy, 9/24/07)
!  (11) Now use the Ott et al. [2010] profiles, and apply consistently with
!        GMI model [Allen et al., 2010] (ltm, bmy, 1/25/11).
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  15 Jun 2012 - Nielsen - INQUIRE finds free logical unit number for IU_FILE
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  22 Oct 2013 - C. Keller   - Now a HEMCO extension.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: L
    REAL*8  :: ZHEIGHT, YMID
    REAL*8  :: FRAC(HcoState%NZ)

      !=================================================================
      ! LIGHTDIST begins here!
      !=================================================================

    ! Initialize
    MTYPE    = 0
    VERTPROF = 0d0

    !%%% NOTE: Use L=1 for GRID_MOD functions.  This is OK for the
    !%%% existing GEOS-Chem with a pure cartesian grid, but may be an
    !%%% issue when interfaced with a GCM with a non-regular grid
    !%%% (bmy, 3/1/12)
    YMID     = HcoState%Grid%YMID%Val( I, J )

    !=================================================================
    ! Test whether location (I,J) is continental, marine, or snow/ice
    !
    ! Depending on the combination of land/water and latitude,
    ! assign a flag describing the type of lightning:
    !
    !   MTYPE = 1: ocean lightning
    !   MTYPE = 2: tropical continental lightning
    !   MTYPE = 3: midlatitude continental lightning
    !   MTYPE = 4: subtropical lightning
    !
    ! (ltm, bmy, 1/25/11)
    !=================================================================

    ! Assign profile kind to grid box, following Allen et al.
    ! [JGR, 2010] (ltm, 1/25,11)

    SELECT CASE ( MONTH )

       ! Southern Hemisphere Summer
       CASE ( 1,2,3,12 )

           IF ( ABS(YMID) .le. 15 ) THEN
              IF ( SFCTYPE == 0 ) THEN
                 MTYPE = 2        ! Tropical continental
              ELSE
                 MTYPE = 1        ! Tropical marine
              ENDIF
           ELSE IF ( ( YMID .gt. 15. ) .and. ( YMID .le. 30. ) ) THEN
              MTYPE = 4           ! N. Subtropics
           ELSE IF ( ( YMID .ge. -40. ) .and. ( YMID .lt. -15. ) ) THEN
              MTYPE = 4           ! S. Subtropics
           ELSE
              MTYPE = 3           ! Midlatitude
           ENDIF

        ! Equinox months
        CASE ( 4,5,10,11 )

           IF ( ABS(YMID) .le. 15 ) THEN
              IF ( SFCTYPE == 0 ) THEN
                 MTYPE = 2        ! Tropical continental
              ELSE
                 MTYPE = 1        ! Tropical marine
              ENDIF
           ELSE IF ( ABS(YMID) .le. 30 ) THEN
              MTYPE = 4           ! Subtropics
           ELSE
              MTYPE = 3           ! Midlatitude
           ENDIF

        ! Northern Hemisphere Summer
        CASE ( 6,7,8,9 )

           IF ( ABS(YMID) .le. 15 ) THEN
              IF ( SFCTYPE == 0 ) THEN
                 MTYPE = 2        ! Tropical continental
              ELSE
                 MTYPE = 1        ! Tropical marine
              ENDIF
           ELSE IF ( ( YMID .gt. 15. ) .and. ( YMID .le. 40. ) ) THEN
              MTYPE = 4           ! N. Subtropics
           ELSE IF ( ( YMID .ge. -30. ) .and. ( YMID .lt. -15. ) ) THEN
              MTYPE = 4           ! S. Subtropics
           ELSE
              MTYPE = 3           ! Midlatitude
           ENDIF

    END SELECT

    ! Extra safety check for pathological grid boxes (bmy, 11/29/06)
    IF ( MTYPE == 0 ) RETURN

    !=================================================================
    ! Use the CDF for this type of lightning to partition the total
    ! column lightning NOx into the layers
    !=================================================================
    ZHEIGHT = 0.0

    ! Compute the height [km] at the top of each vertical level.
    ! Look up the cumulative fraction of NOx for each vertical level
    DO L = 1, LTOP
       ZHEIGHT = ZHEIGHT + HcoState%Grid%BXHEIGHT_M%Val(I,J,L)
       FRAC(L) = Inst%PROFILE( NINT( ( ZHEIGHT/H0 )*3200. ), MTYPE ) *0.01
    ENDDO

    ! Convert from cumulative fraction to fraction for each level
    DO L = LTOP, 2, - 1
       FRAC(L) = FRAC(L) - FRAC(L-1)
    ENDDO

    ! Partition lightning NOx by layer into VERTPROF
    DO L = 1, LTOP
       VERTPROF(L) = ( FRAC(L) * TOTAL )
    ENDDO

  END SUBROUTINE LightDist
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_LightNOx_Init
!
! !DESCRIPTION: Subroutine HCOX\_LIGHTNOX\_INIT allocates all module arrays.
!  It also reads the lightning CDF data from disk before the first lightning
!  timestep.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_LightNOx_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_Chartools_Mod, ONLY : HCO_CharParse
    USE HCO_ExtList_Mod,   ONLY : GetExtNr
    USE HCO_ExtList_Mod,   ONLY : GetExtOpt
    USE HCO_ExtList_Mod,   ONLY : GetExtSpcVal
    USE HCO_State_Mod,     ONLY : HCO_GetHcoID
    USE HCO_State_Mod,     ONLY : HCO_GetExtHcoID
    USE HCO_ReadList_Mod,  ONLY : ReadList_Remove
    USE inquireMod,        ONLY : findfreeLUN
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState   ! Hemco options
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName    ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState     ! Module options
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)  :: RC
!
! !REVISION HISTORY:
!  14 Apr 2004 - R. Yantosca - Initial version
!  (1 ) Now reference DATA_DIR from "directory_mod.f"
!  (2 ) Now call GET_MET_FIELD_SCALE to initialize the scale factor for
!        each met field type and grid resolution (bmy, 8/25/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now get the box area at 30N for MFLUX, PRECON (lth, bmy, 5/10/06)
!  (5 ) Rename OTDSCALE to OTD_REG_REDIST.  Also add similar array
!        OTD_LOC_REDIST.  Now call GET_FLASH_SCALE_CTH, GET_FLASH_SCALE_MFLUX,
!        GET_FLASH_SCALE_PRECON depending on the type of lightning param used.
!        Updated comments.  (ltm, bmy, 1/31/07)
!  (6 ) Removed near-land stuff.  Renamed from HCOX_LightNOX_Init_NL to
!        HCOX_LightNOX_Init.  Now allocate EMIS_LI_NOx. (ltm, bmy, 10/3/07)
!  (7 ) Also update location of PDF file to lightning_NOx_200709 directory.
!        (bmy, 1/24/08)
!  (8 ) Read in new Ott profiles from lightning_NOx_201101. Remove
!        depreciated options. (ltm, bmy, 1/25/11)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Removed reference to GET_YEDGE
!  22 Oct 2013 - C. Keller   - Now a HEMCO extension.
!  26 Feb 2015 - R. Yantosca - Now re-introduce reading the CDF table from an
!                              ASCII file (reduces compilation time)
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: AS, III, IOS, JJJ, IU_FILE, nSpc
    INTEGER                        :: ExtNr
    LOGICAL                        :: FOUND, FileExists
    INTEGER, ALLOCATABLE           :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
    CHARACTER(LEN=255)             :: MSG, LOC, FILENAME, FileMsg
    TYPE(MyInst), POINTER          :: Inst

    !=======================================================================
    ! HCOX_LightNOX_Init begins here!
    !=======================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_LightNOx_Init (hcox_lightnox_mod.F90)', RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create LightNOx instance for this simulation
    Inst => NULL()
    CALL InstCreate ( ExtNr, ExtState%LightNOx, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot create LightNOx instance', RC )
       RETURN
    ENDIF

    !=======================================================================
    ! Obtain lightning CDF's from Ott et al [JGR, 2010].
    !
    ! PART 1 --- Move the file name check to the front of this routine to
    ! facilitate the GEOS-Chem dry-run and HEMCO-standalone dry-run.
    !=======================================================================

    ! Get filename from configuration file
    CALL GetExtOpt( HcoState%Config, ExtNr, 'CDF table',                     &
                    OptValChar=FILENAME, RC=RC                              )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Call HEMCO parser to replace tokens such as $ROOT, $MET, or $RES.
    ! There shouldn't be any date token in there ($YYYY, etc.), so just
    ! provide some dummy variables here
    CALL HCO_CharParse( HcoState%Config, FILENAME, -999, -1, -1, -1, -1, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------------
    ! In dry-run mode, print file path to dryrun log and exit.
    ! Otherwise, print file path to the HEMCO log file and continue.
    !-----------------------------------------------------------------------

    ! Test if the file exists
    INQUIRE( FILE=TRIM( FileName ), EXIST=FileExists )

    ! Create a display string based on whether or not the file is found
    IF ( FileExists ) THEN
       FileMsg = 'HEMCO (LIGHTNOX): Opening'
    ELSE
       FileMsg = 'HEMCO (LIGHTNOX): REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write file status to stdout and the HEMCO log
    IF ( HcoState%amIRoot ) THEN
       WRITE( 6,   300 ) TRIM( FileMsg ), TRIM( FileName )
       WRITE( MSG, 300 ) TRIM( FileMsg ), TRIM( FileName )
       CALL HCO_MSG( HcoState%Config%Err, MSG )
 300   FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulation, return to calling program.
    ! For regular simulations, throw an error if we can't find the file.
    IF ( HcoState%Options%IsDryRun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( MSG, 300 ) TRIM( FileMsg ), TRIM( FileName )
          CALL HCO_ERROR(HcoState%Config%Err, MSG, RC )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Exit if this is a GEOS-Chem or HEMCO-standalone dry-run
    !=======================================================================
    IF ( HcoState%Options%IsDryRun ) THEN
       Inst => NULL()
       CALL HCO_LEAVE( HcoState%Config%Err,RC )
       RETURN
    ENDIF

    !=======================================================================
    ! Continue for regular simulations ...
    !=======================================================================

    ! Check for usage of convective fractions. This becomes only active
    ! if both the convective fraction and the buoyancy field are available.
    CALL GetExtOpt( HcoState%Config, ExtNr, 'Use CNV_FRC', &
                     OptValBool=Inst%LCNVFRC, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) Inst%LCNVFRC = .FALSE.

    ! Check for usage of GEOS-5 lightning flash rates. If on, the GEOS-5
    ! flash rates (where available) are used instead of the computed flash
    ! rates. This is off by default.
    CALL GetExtOpt( HcoState%Config, ExtNr, 'GEOS-5 flash rates', &
                     OptValBool=Inst%LLFR, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) Inst%LLFR = .FALSE.

    ! Get species ID
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( nSpc /= 1 ) THEN
       MSG = 'Lightning NOx module must have exactly one species!'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF
    Inst%IDTNO = HcoIDs(1)

    ! Get species scale factor
    CALL GetExtSpcVal( HcoState%Config, ExtNr, nSpc, &
                       SpcNames, 'Scaling', 1.0_sp, Inst%SpcScalVal, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL GetExtSpcVal( HcoState%Config, ExtNr, nSpc, &
                       SpcNames, 'ScaleField', HCOX_NOSCALE, Inst%SpcScalFldNme, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Echo info about this extension
    IF ( HcoState%amIRoot ) THEN
       MSG = 'Use lightning NOx emissions (extension module)'
       CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
       WRITE(MSG,*) ' - Use species ', TRIM(SpcNames(1)), '->', Inst%IDTNO
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - Use GEOS-5 flash rates: ', Inst%LLFR
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - Use scalar scale factor: ', Inst%SpcScalVal(1)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - Use gridded scale field: ', TRIM(Inst%SpcScalFldNme(1))
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    !=======================================================================
    ! Allocate arrays
    !=======================================================================

    ALLOCATE( Inst%PROFILE( NNLIGHT, NLTYPE ), STAT=AS )
    IF( AS /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'PROFILE', RC )
       RETURN
    ENDIF
    Inst%PROFILE = 0.0_hp

    ALLOCATE( Inst%SLBASE(HcoState%NX,HcoState%NY,HcoState%NZ), STAT=AS )
    IF( AS /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'SLBASE', RC )
       RETURN
    ENDIF
    Inst%SLBASE = 0.0_hp

    ALLOCATE ( Inst%FLASH_DENS_TOT( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLASH_DENS_TOT', RC )
       RETURN
    ENDIF
    Inst%FLASH_DENS_TOT = 0.0_sp

    ALLOCATE ( Inst%FLASH_DENS_IC( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLASH_DENS_IC', RC )
       RETURN
    ENDIF
    Inst%FLASH_DENS_IC = 0.0_sp

    ALLOCATE ( Inst%FLASH_DENS_CG( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLASH_DENS_CG', RC )
       RETURN
    ENDIF
    Inst%FLASH_DENS_CG = 0.0_sp

    ALLOCATE ( Inst%CONV_DEPTH( HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'CONV_DEPTH', RC )
       RETURN
    ENDIF
    Inst%CONV_DEPTH = 0.0_sp

    !=======================================================================
    ! Obtain lightning CDF's from Ott et al [JGR, 2010].
    !
    ! PART 2 --- Read the data!
    !=======================================================================

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open file containing lightning NOx PDF data
    OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       MSG = 'IOERROR: LightDist: 1'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! Read 12 header lines
    DO III = 1, 12
       READ( IU_FILE, '(a)', IOSTAT=IOS )
       IF ( IOS /= 0 ) THEN
          MSG = 'IOERROR: LightDist: 2'
          CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
          RETURN
       ENDIF
    ENDDO

    ! Read NNLIGHT types of lightning profiles
    DO III = 1, NNLIGHT
       READ( IU_FILE,*,IOSTAT=IOS) (Inst%PROFILE(III,JJJ),JJJ=1,NLTYPE)
       IF ( IOS /= 0 ) THEN
          MSG = 'IOERROR: LightDist: 3'
          CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
          RETURN
       ENDIF
    ENDDO

    ! Close file
    CLOSE( IU_FILE )

    !=======================================================================
    ! Create diagnostics for lightning flash rates and convective cloud height
    !=======================================================================
    CALL Diagn_Create( HcoState  = HcoState,                        &
                       cName     = 'LightningFlashRate_Total' ,     &
                       ExtNr     = ExtNr,                           &
                       Cat       = -1,                              &
                       Hier      = -1,                              &
                       HcoID     = -1,                              &
                       SpaceDim  = 2,                               &
                       OutUnit   = 'flashes/min/km2',               &
                       AutoFill  = 0,                               &
                       Trgt2D    = Inst%FLASH_DENS_TOT,             &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,                        &
                       cName     = 'LightningFlashRate_IntraCloud', &
                       ExtNr     = ExtNr,                           &
                       Cat       = -1,                              &
                       Hier      = -1,                              &
                       HcoID     = -1,                              &
                       SpaceDim  = 2,                               &
                       OutUnit   = 'flashes/min/km2',               &
                       AutoFill  = 0,                               &
                       Trgt2D    = Inst%FLASH_DENS_IC,              &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,                         &
                       cName     = 'LightningFlashRate_CloudGround', &
                       ExtNr     = ExtNr,                            &
                       Cat       = -1,                               &
                       Hier      = -1,                               &
                       HcoID     = -1,                               &
                       SpaceDim  = 2,                                &
                       OutUnit   = 'flashes/min/km2',                &
                       AutoFill  = 0,                                &
                       Trgt2D    = Inst%FLASH_DENS_CG,               &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,                         &
                       cName     = 'ConvectiveCloudTopHeight',       &
                       ExtNr     = ExtNr,                            &
                       Cat       = -1,                               &
                       Hier      = -1,                               &
                       HcoID     = -1,                               &
                       SpaceDim  = 2,                                &
                       OutUnit   = '1',                              &
                       AutoFill  = 0,                                &
                       Trgt2D    = Inst%CONV_DEPTH,                  &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=======================================================================
    ! Activate met fields required by this module
    !=======================================================================
    ExtState%TK%DoUse         = .TRUE.
    ExtState%TROPP%DoUse      = .TRUE.
    ExtState%CNV_MFC%DoUse    = .TRUE.
    ExtState%CNV_FRC%DoUse    = .TRUE.
    ExtState%ALBD%DoUse       = .TRUE.
    ExtState%WLI%DoUse        = .TRUE.
    ExtState%LFR%DoUse        = .TRUE.
    ExtState%FLASH_DENS%DoUse = .TRUE.
    ExtState%CONV_DEPTH%DoUse = .TRUE.

    ! Only activate BYNCY and LFR if they are needed
    IF ( Inst%LCNVFRC .OR. Inst%LLFR ) ExtState%BYNCY%DoUse = .TRUE.
    IF ( Inst%LLFR ) ExtState%LFR%DoUse = .TRUE.

    ! Cleanup
    Inst => NULL()

    ! Leave w/ success
    IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_LightNOx_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hcox_lightnox_final
!
! !DESCRIPTION: Subroutine HCOX\_LIGHTNOX\_FINAL deallocates all module
!  arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_LightNOx_Final( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !REVISION HISTORY:
!  14 Apr 2004 - R. Yantosca - Initial version
!  (1 ) Now deallocates OTDSCALE (ltm, bmy, 5/10/06)
!  (2 ) Rename OTDSCALE to OTD_REG_REDIST.  Now deallocate OTD_LOC_REDIST.
!        (bmy, 1/31/07)
!  (3 ) Renamed from HCOX_LightNOX_Final_NL to HCOX_LightNOX_Final.
!        Now deallocate EMIS_LI_NOx. (ltm, bmy, 10/3/07)
!  (4 ) Remove depreciated options. (ltm, bmy, 1/25/11)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  22 Oct 2013 - C. Keller   - Now a HEMCO extension.
!  22 Jul 2014 - R. Yantosca - PROFILE is now set in lightning_cdf_mod.F90
!  26 Feb 2015 - R. Yantosca - Now re-introduce PROFILE, as we read the CDF
!                              table from an ASCII file (reduces compile time)
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! Cleanup module arrays
    !=================================================================
    CALL InstRemove ( ExtState%LightNOx )

  END SUBROUTINE HCOX_LightNOx_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a pointer to the desired instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstGet ( Instance, Inst, RC, PrevInst )
!
! !INPUT PARAMETERS:
!
    INTEGER                             :: Instance
    TYPE(MyInst),     POINTER           :: Inst
    INTEGER                             :: RC
    TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst),     POINTER    :: PrvInst

    !=================================================================
    ! InstGet begins here!
    !=================================================================

    ! Get instance. Also archive previous instance.
    PrvInst => NULL()
    Inst    => AllInst
    DO WHILE ( ASSOCIATED(Inst) )
       IF ( Inst%Instance == Instance ) EXIT
       PrvInst => Inst
       Inst    => Inst%NextInst
    END DO
    IF ( .NOT. ASSOCIATED( Inst ) ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Pass output arguments
    IF ( PRESENT(PrevInst) ) PrevInst => PrvInst

    ! Cleanup & Return
    PrvInst => NULL()
    RC = HCO_SUCCESS

  END SUBROUTINE InstGet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate
!
! !DESCRIPTION: Subroutine InstCreate adds a new instance to the list of
!  instances, assigns a unique instance number to this new instance, and
!  archives this instance number to output argument Instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstCreate ( ExtNr, Instance, Inst, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,       INTENT(IN)       :: ExtNr
!
! !OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(  OUT)    :: Instance
    TYPE(MyInst),  POINTER          :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(INOUT)    :: RC
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst), POINTER          :: TmpInst
    INTEGER                        :: nnInst

    !=================================================================
    ! InstCreate begins here!
    !=================================================================

    ! ----------------------------------------------------------------
    ! Generic instance initialization
    ! ----------------------------------------------------------------

    ! Initialize
    Inst => NULL()

    ! Get number of already existing instances
    TmpInst => AllInst
    nnInst = 0
    DO WHILE ( ASSOCIATED(TmpInst) )
       nnInst  =  nnInst + 1
       TmpInst => TmpInst%NextInst
    END DO

    ! Create new instance
    ALLOCATE(Inst)
    Inst%Instance = nnInst + 1
    Inst%ExtNr    = ExtNr

    ! Attach to instance list
    Inst%NextInst => AllInst
    AllInst       => Inst

    ! Update output instance
    Instance = Inst%Instance

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstRemove
!
! !DESCRIPTION: Subroutine InstRemove removes an instance from the list of
! instances.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstRemove ( Instance )
!
! !INPUT PARAMETERS:
!
    INTEGER                         :: Instance
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                     :: RC
    TYPE(MyInst), POINTER       :: PrevInst => NULL()
    TYPE(MyInst), POINTER       :: Inst     => NULL()

    !=================================================================
    ! InstRemove begins here!
    !=================================================================

    ! Get instance. Also archive previous instance.
    PrevInst => NULL()
    Inst     => NULL()
    CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

    ! Instance-specific deallocation
    IF ( ASSOCIATED(Inst) ) THEN
       ! Free pointer
       IF ( ASSOCIATED( Inst%PROFILE       ) ) DEALLOCATE ( Inst%PROFILE       )
       IF ( ASSOCIATED( Inst%SLBASE        ) ) DEALLOCATE ( Inst%SLBASE        )
       IF ( ASSOCIATED( Inst%FLASH_DENS_TOT) ) DEALLOCATE ( Inst%FLASH_DENS_TOT)
       IF ( ASSOCIATED( Inst%FLASH_DENS_IC ) ) DEALLOCATE ( Inst%FLASH_DENS_IC )
       IF ( ASSOCIATED( Inst%FLASH_DENS_CG ) ) DEALLOCATE ( Inst%FLASH_DENS_CG )
       IF ( ASSOCIATED( Inst%CONV_DEPTH    ) ) DEALLOCATE ( Inst%CONV_DEPTH    )
       IF ( ALLOCATED ( Inst%SpcScalVal    ) ) DEALLOCATE ( Inst%SpcScalVal    )
       IF ( ALLOCATED ( Inst%SpcScalFldNme ) ) DEALLOCATE ( Inst%SpcScalFldNme )

       ! Pop off instance from list
       IF ( ASSOCIATED(PrevInst) ) THEN
          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL()
    ENDIF

   END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_LightNOx_Mod
