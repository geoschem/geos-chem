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
  PRIVATE :: FLASHES_CTH              
  PRIVATE :: GET_IC_CG_RATIO          
  PRIVATE :: GET_OTD_LIS_SCALE        
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  INTEGER, PARAMETER            :: NLTYPE        = 4
  INTEGER, PARAMETER            :: NNLIGHT       = 3200
  REAL*8,  PARAMETER            :: RFLASH_MIDLAT = 3.011d26   ! 500 mol/flash
  REAL*8,  PARAMETER            :: RFLASH_TROPIC = 1.566d26   ! 260 mol/flash
  REAL*8,  PARAMETER            :: EAST_WEST_DIV = -30d0
  REAL*8,  PARAMETER            :: WEST_NS_DIV   =  35d0
  REAL*8,  PARAMETER            :: EAST_NS_DIV   =  35d0
  REAL*8,  PARAMETER            :: T_NEG_BOT     = 273.0d0    !   0 C 
  REAL*8,  PARAMETER            :: T_NEG_CTR     = 258.0d0    ! -15 C
  REAL*8,  PARAMETER            :: T_NEG_TOP     = 233.0d0    ! -40 C
                               
  ! testing only               
  integer, parameter            :: ix = -1 !30 !19 
  integer, parameter            :: iy = -1 !6  !33 
  integer, parameter            :: iz = -1 !9  !9
!
! !PRIVATE TYPES:
!
  ! Scalars
  TYPE :: MyInst
   INTEGER                       :: Instance
   INTEGER                       :: IDTNO     ! NO tracer ID
   INTEGER                       :: ExtNr     ! HEMCO Extension ID
   LOGICAL                       :: DoDiagn 
!   REAL*8                        :: AREA_30N
   REAL*8                        :: OTD_LIS_SCALE
   LOGICAL                       :: OTD_LIS_PRESC ! Is OTD_LIS_SCALE prescribed?
   LOGICAL                       :: LOTDLOC       ! Use OTD-LIS dist factors?
   LOGICAL                       :: LCNVFRC       ! Use convective fractions? 
   LOGICAL                       :: LLFR          ! Use GEOS-5 flash rates 

   ! Arrays
   REAL(dp), POINTER             :: PROFILE(:,:)
   REAL(hp), POINTER             :: SLBASE(:,:,:)

   ! OTD scale factors read through configuration file
   REAL(hp), POINTER :: OTDLIS(:,:) => NULL()

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
! to calculate lightnox NOx emissions and return them to the HEMCO
! driver routine.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_LightNOx_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_FluxArr_Mod,  ONLY : HCO_EmisAdd 
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root
    TYPE(Ext_State), POINTER        :: ExtState    ! Module options
    TYPE(HCO_State), POINTER        :: HcoState   ! Hemco options 
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
    CALL LIGHTNOX ( am_I_Root, HcoState, ExtState, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=================================================================
    ! Pass to HEMCO State and update diagnostics 
    !=================================================================
    IF ( Inst%IDTNO > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( am_I_Root, HcoState, Inst%SLBASE, Inst%IDTNO, & 
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
!  computing NOx emission from lightnox (with various updates).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LightNOx( am_I_Root, HcoState, ExtState, Inst, RC )
!
! !USES:
!
    USE HCO_Calc_Mod,     ONLY : HCO_EvalFld
    USE HCO_EmisList_Mod, ONLY : HCO_GetPtr      
    USE HCO_GeoTools_Mod, ONLY : HCO_LANDTYPE
    USE HCO_Clock_Mod,    ONLY : HcoClock_Get
    USE HCO_Clock_Mod,    ONLY : HcoClock_First
    USE HCO_ExtList_Mod,  ONLY : GetExtOpt
    USE HCO_Types_Mod,    ONLY : DiagnCont
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root
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
    INTEGER           :: I,         J,           L,        LCHARGE
    INTEGER           :: LMAX,      LTOP,        LBOTTOM,  L_MFLUX
    INTEGER           :: cMt,       MTYPE,       LTOP1,    LTOP2 
    INTEGER           :: LTOP_LFR
    REAL*8            :: A_KM2,     A_M2,        CC,       DLNP     
    REAL*8            :: DZ,        FLASHRATE,   H0,       HBOTTOM
    REAL*8            :: HCHARGE,   IC_CG_RATIO, MFLUX,    P1
    REAL*8            :: P2,        P3,          RAIN,     RATE
    REAL*8            :: RATE_SAVE, REDIST,      T1,       T2
    REAL*8            :: TOTAL,     TOTAL_CG,    TOTAL_IC, X       
    REAL*8            :: YMID,      Z_IC,        Z_CG,     ZUP
    REAL*8            :: XMID
    REAL*8            :: VERTPROF(HcoState%NZ)
    REAL*8            :: RATE_LFR, IC_CG_RATIO_LFR, H0_LFR
    INTEGER           :: LNDTYPE, SFCTYPE
    INTEGER           :: DiagnID
    REAL(hp), TARGET  :: DIAGN(HcoState%NX,HcoState%NY,3)
    REAL(hp), POINTER :: Arr2D(:,:)
    TYPE(DiagnCont), POINTER :: TmpCnt
    REAL(hp)          :: TROPP
    REAL(dp)          :: TmpScale

    ! Cloud top height
    REAL(hp), TARGET  :: TOPDIAGN(HcoState%NX,HcoState%NY)

    !=================================================================
    ! LIGHTNOX begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'LightNOx (hcox_lightnox_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Init
    Arr2D  => NULL() 
    TmpCnt => NULL()

    ! ----------------------------------------------------------------
    ! First call routines
    ! ----------------------------------------------------------------
    IF ( HcoClock_First( HcoState%Clock, .TRUE. ) ) THEN

       ! ckeller, 8/30/18: Always write out diagnostics
       Inst%DoDiagn = .TRUE.

       ! Get scale factor. 
       ! - Try to read from configuration file first.
       CALL GetExtOpt( HcoState%Config, Inst%ExtNr, 'OTD-LIS scaling', &
                        OptValDp = TmpScale, FOUND=Inst%OTD_LIS_PRESC, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( Inst%OTD_LIS_PRESC ) THEN
          Inst%OTD_LIS_SCALE = TmpScale
       ! - Get according to compiler switches otherwise
       ELSE
          CALL GET_OTD_LIS_SCALE( am_I_Root, HcoState, Inst%OTD_LIS_SCALE, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

    ENDIF

    ! Eventually get OTD-LIS local redistribution factors from HEMCO.
    IF ( Inst%LOTDLOC ) THEN
       CALL HCO_EvalFld( am_I_Root, HcoState, 'LIGHTNOX_OTDLIS', Inst%OTDLIS, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Reset arrays 
    Inst%SLBASE = 0.0_hp
    IF (Inst%DoDiagn) THEN
       DIAGN    = 0.0_hp
       TOPDIAGN = 0.0_hp
    ENDIF

    ! LMAX: the highest L-level to look for lightnox (usually LLPAR-1)
    LMAX   = HcoState%NZ - 1

    ! Get current month (to be passed to LIGHTDIST)
    CALL HcoClock_Get( am_I_Root, HcoState%Clock, cMM=cMt, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=================================================================
    ! Compute lightnox emissions for each (I,J) column
    !=================================================================

!$OMP PARALLEL DO                                                   &
!$OMP DEFAULT( SHARED )                                             &
!$OMP PRIVATE( I,           J,        L,        A_M2,   A_KM2     ) &
!$OMP PRIVATE( YMID,        XMID,     LCHARGE,  P1,     P2        ) &
!$OMP PRIVATE( T1,          T2,       DLNP,     DZ,     P3        ) &
!$OMP PRIVATE( ZUP,         HCHARGE,  LTOP,     H0,     Z_CG      ) &
!$OMP PRIVATE( Z_IC,        LBOTTOM,  HBOTTOM,  CC,     FLASHRATE ) &
!$OMP PRIVATE( IC_CG_RATIO, L_MFLUX,  MFLUX,    RAIN,   RATE      ) &
!$OMP PRIVATE( X,           TOTAL_IC, TOTAL_CG, TOTAL,  REDIST    ) &
!$OMP PRIVATE( RATE_SAVE,   VERTPROF, SFCTYPE,  LNDTYPE, TROPP    ) &
!$OMP PRIVATE( MTYPE,       LTOP1,    LTOP2                       ) &
!$OMP PRIVATE( RATE_LFR,    H0_LFR,   IC_CG_RATIO_LFR, LTOP_LFR   ) &
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

       ! Initialize
       LBOTTOM       = 0 
       LCHARGE       = 0
       CC            = 0d0
       HCHARGE       = 0d0
       HBOTTOM       = 0d0
       TOTAL         = 0d0
       TOTAL_IC      = 0d0
       TOTAL_CG      = 0d0
       
       ! Get factors for OTD-LIS local redistribution or none.
       ! This constrains the seasonality and spatial distribution
       ! of the parameterized lightnox to match the HRMC v2.2
       ! product from LIS/OTD, while still allowing the model to
       ! place lightnox locally within deep convective events.
       ! (ltm, bmy, 1/31/07)
       IF ( Inst%LOTDLOC ) THEN
          REDIST = Inst%OTDLIS(I,J)
       ELSE
          REDIST = 1.0d0
       ENDIF

       !===========================================================
       ! Use GEOS-5 LFR if option is selected
       ! LFR is in flashes s-1 km-2
       ! Use arbitrary value of 1.0 for IC/CG ratio. This value has
       ! no impact on the total lightning emissions. 
       !===========================================================
       IF ( Inst%LLFR                                  &
            .AND. ASSOCIATED( ExtState%LFR%Arr%Val   ) &
            .AND. ASSOCIATED( ExtState%BYNCY%Arr%Val )  ) THEN
          RATE_LFR        = ExtState%LFR%Arr%Val(I,J) * A_KM2 * 60.0d0 * 360.0d0
          IC_CG_RATIO_LFR = 1.0

          ! Apply scaling factor to make sure annual average flash rate 
          ! equals that of the climatology. (ltm, 09/24/07)
          RATE_LFR = RATE_LFR * Inst%OTD_LIS_SCALE

          ! Set LTOP to top of buoyancy
          DO L = HcoState%NZ, 1, -1
             IF ( ExtState%BYNCY%Arr%Val(I,J,L) >= 0.0_sp ) THEN
                LTOP_LFR = L + 1
                EXIT
             ENDIF
          ENDDO
          !LTOP = MAX( LTOP, LMAX )

          ! H0 is the convective cloud top height [m].  This is the
          ! distance from the surface to the top edge of box (I,J,LTOP).
          H0_LFR = SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:LTOP_LFR))

       ELSE
          IC_CG_RATIO_LFR = -999.0
          RATE_LFR        = -999.0
          H0_LFR          = -999.0
          LTOP_LFR        = -999.0

       ENDIF

       ! Init
       RATE        = -999.0
       H0          = -999.0
       IC_CG_RATIO = -999.0

       !===========================================================
       ! (1) FIND NEGATIVE CHARGE LAYER
       !
       ! LCHARGE is the L-value where the negative charge layer is
       ! found.  According to Williams (1985), the negative charge
       ! layer occurs where T is between 0 C and -40 C.  The 
       ! original model code set this at -10 C, but according to 
       ! Houze (1993), a good proxy for the negative charge layer
       ! maximum density is at -15 C.
       !
       ! Also of interest for later, will be the bottom of the
       ! negative charge layer (i.e., temp = 0 C) in calculating 
       ! the cold cloud depth.
       !
       ! If LCHARGE=1, then it is too cold to have water droplets
       ! in the column, so there will be no lightnox events,
       ! and we go to the next (I,J) box.
       !
       ! (ltm, bmy, 5/10/06, 12/11/06)
       !===========================================================

       ! Find negative charge layer
       DO L = 1, LMAX
          IF ( ExtState%TK%Arr%Val(I,J,L) <= T_NEG_CTR ) THEN
             LCHARGE = L
             EXIT
          ENDIF
       ENDDO

       ! Error check LCHARGE
       IF ( LCHARGE >= LMAX ) LCHARGE = LMAX
       !IF ( LCHARGE <= 1    ) CYCLE
       IF ( LCHARGE > 1 ) THEN

       !-----------------------------------------------------------
       ! (1a) Define more quantities
       !-----------------------------------------------------------
           
       ! Pressure [Pa] at the centers of grid
       ! boxes (I,J,LCHARGE-1) and (I,J,LCHARGE)
       ! Now calculate from grid edges (ckeller, 10/06/2014)
       P1 = ( HcoState%Grid%PEDGE%Val(I,J,LCHARGE-1) &
            + HcoState%Grid%PEDGE%Val(I,J,LCHARGE  ) ) / 2.0_hp
       P2 = ( HcoState%Grid%PEDGE%Val(I,J,LCHARGE  ) &
            + HcoState%Grid%PEDGE%Val(I,J,LCHARGE+1) ) / 2.0_hp

       ! Temperatures [K] at the centers of grid
       ! boxes (I,J,LCHARGE-1) and (I,J,LCHARGE)
       T1   = ExtState%TK%Arr%Val(I,J,LCHARGE-1)
       T2   = ExtState%TK%Arr%Val(I,J,LCHARGE  )
 
       ! DZ is the height [m] from the center of box (I,J,LCHARGE-1)
       ! to the negative charge layer.  It may be found in either
       ! the (LCHARGE)th sigma layer or the (LCHARGE-1)th layer.
       ! We use the hypsometric eqn to find the distance between
       ! the center of (LCHARGE)th and (LCHARGE-1)th boxes, then
       ! assume a linear temp distribution to scale between the two.
       DLNP = LOG( P1 / P2 ) / ( T1 - T2 ) * ( T1 - T_NEG_CTR )
       DZ   = HcoState%Phys%Rdg0 * ( ( T1 + T2 ) / 2d0 ) * DLNP

       ! Pressure [Pa] at the bottom edge of box (I,J,LCHARGE),
       ! or, equivalently, the top edge of box (I,J,LCHARGE-1).
       P3   = HcoState%Grid%PEDGE%Val( I, J, LCHARGE )

       ! Height [m] from the center of grid box (I,J,LCHARGE-1) 
       ! to the top edge of grid box (I,J,LCHARGE-1)
       ZUP  = HcoState%Phys%Rdg0 * T1 * LOG( P1 / P3 )

       !-----------------------------------------------------------
       ! (1b) HCHARGE is the height of the negative charge layer 
       ! above the bottom edge of box (I,J,LCHARGE).  
       ! 
       ! If DZ < ZUP, then DZ is in grid box (I,J,LCHARGE-1);
       ! therefore subtract 1 from LCHARGE and compute HCHARGE 
       ! accordingly.  
       !
       ! In this case, please note that BXHEIGHT(I,J,LCHARGE)-ZUP 
       ! is the distance from the bottom edge of the grid box to 
       ! the center of the newly defined (LCHARGE)th layer.
       !-----------------------------------------------------------
       IF ( DZ >= ZUP ) THEN
          HCHARGE = DZ - ZUP
       ELSE
          LCHARGE = LCHARGE - 1
          HCHARGE = (HcoState%Grid%BXHEIGHT_M%Val(I,J,LCHARGE)-ZUP) + DZ
       ENDIF
 
       !===========================================================
       ! (2) COMPUTE CONVECTIVE CLOUD TOP HEIGHT
       !
       ! LTOP is the L-layer where the convective cloud top is 
       ! found.  The cloud top is located at the highest sigma 
       ! level for which the cloud mass flux is nonzero.  Since 
       ! GMAO cloud mass flux is defined at the top of each sigma 
       ! level, the convective cloud top is located at the top 
       ! edge of layer LTOP.
       !
       ! For lightnox to exist, the cloud must straddle the 
       ! negative charge layer (in other words, at the very 
       ! minimum, the cloud bottom must occur in the LCHARGEth 
       ! layer).  If LTOP < LCHARGE go to the next (I,J) location.
       ! 
       ! Additionally, because the negative charge layer extends 
       ! from 0 C to around -40 C (Williams 1985), any cloud type 
       ! heights that are not colder than -40 C will be considered 
       ! unable to create the necessary dipole.  Therefore, if 
       ! T(I,J,LTOP) >= -40 C, go to the next (I,J) location. 
       !
       ! To be easily translatable to an ESMF environment, we now 
       ! use the convective cloud mass flux to determine LTOP.
       ! Use the same definition as used in GEOS-Chem.
       !
       ! (ltm, bmy, 5/10/06, 12/11/06)
       !
       ! GEOS-FP turns off convection in grid boxes where vertical
       ! transport is explicitly resolved. The convective mass flux
       ! (which is computed within the convection code) is then zero
       ! in these grid boxes, even though convection did occur at 
       ! these places. 
       ! This may become increasingly relevant as GEOS-FP operates 
       ! at even higher resolutions. 
       ! If available, also determine cloud top height from 
       ! buoyancy and the convective fraction. Define it as the
       ! highest level with non-negative buoyancy and for columns 
       ! with non-zero convective fraction (ckeller, 3/04/16).
       !===========================================================

       ! 'Traditional definition of cloud top level
       LTOP1 = 1
       DO L = HcoState%NZ, 1, -1
          IF ( ExtState%CNV_MFC%Arr%Val(I,J,L) > 0.0_hp ) THEN
             LTOP1 = L + 1
             EXIT
          ENDIF
       ENDDO 

       ! To determine cloud top height from buoyancy for all grid
       ! boxes with non-zero convective fraction (define cloud top
       ! as top level with positive buoyancy).
       LTOP2 = 0
       IF ( Inst%LCNVFRC                         .AND. &
            ASSOCIATED(ExtState%BYNCY%Arr%Val  ) .AND. &
            ASSOCIATED(ExtState%CNV_FRC%Arr%Val)        ) THEN
          IF ( ExtState%CNV_FRC%Arr%Val(I,J) > 0.0_sp ) THEN
             DO L = HcoState%NZ, 1, -1
                IF ( ExtState%BYNCY%Arr%Val(I,J,L) >= 0.0_sp ) THEN 
                   LTOP2= L + 1
                   EXIT
                ENDIF
             ENDDO
          ENDIF 
       ENDIF 

       ! Take whichever value is higher
       LTOP = MAX(LTOP1,LTOP2)

       !----------------------------------------------------------------
       ! Error checks for LTOP 
       !----------------------------------------------------------------

       ! Error check LTOP
       !IF ( LTOP == 0 ) CYCLE
       IF ( LTOP > 0 ) THEN

       ! Error check LTOP as described above
       IF ( LTOP        >  LMAX      ) LTOP = LMAX
       !IF ( LTOP        <  LCHARGE   ) CYCLE
       IF ( LTOP        >= LCHARGE   ) THEN

       ! Diagnose used LTOP
       IF ( Inst%DoDiagn ) THEN
          TOPDIAGN(I,J) = LTOP
       ENDIF

       ! H0 is the convective cloud top height [m].  This is the
       ! distance from the surface to the top edge of box (I,J,LTOP).
       H0   = SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:LTOP))

       ! Z_CG is the cloud-ground path (ground --> HCHARGE) [m]
       Z_CG = SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:LCHARGE-1)) + HCHARGE

       ! Z_IC is the intra-cloud path (HCHARGE --> cloud top) [m]
       Z_IC = SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,LCHARGE:LTOP)) - HCHARGE

       !===========================================================
       ! (3) COMPUTE COLD CLOUD THICKNESS
       ! 
       ! Find the cold cloud thickness (CC) -- the distance from 
       ! where the temperature is 0 C up to the top of the cloud.  
       ! This is necessary for calculating the f_CG/f_IC ratio as 
       ! per Price and Rind 1993.  
       !
       ! This is a clone of the method above to find height to 
       ! HCHARGE, and we can recycle many of the same variables 
       ! that aren't used again.
       ! 
       ! Grid box (I,J,LBOTTOM) is the model layer where the 
       ! temperature of the cloud is 0C.
       !
       ! NOTE: If no temperature in the column is above 0 C, it 
       ! moves on to the next (I,J) box as before with the -15 C.
       !
       ! (ltm, bmy, 5/10/06, 12/11/06)
       !===========================================================

       ! Find the level where T = 0 C
       DO L = 1, LMAX
          IF ( ExtState%TK%Arr%Val(I,J,L) <= T_NEG_BOT ) THEN
             LBOTTOM = L
             EXIT
          ENDIF
       ENDDO

       ! Error check LBOTTOM as described above
       IF ( LBOTTOM >= LMAX ) LBOTTOM = LMAX
       !IF ( LBOTTOM <= 1    ) CYCLE 
       IF ( LBOTTOM > 1    ) THEN 

       !-----------------------------------------------------------
       ! (3a) Define more quantities
       !-----------------------------------------------------------

       ! Pressure [Pa] at the centers of grid
       ! boxes (I,J,LBOTTOM-1) and (I,J,LBOTTOM)
       ! Now calculate from grid edges (ckeller, 10/06/2014)
       P1 = ( HcoState%Grid%PEDGE%Val(I,J,LBOTTOM-1) &
            + HcoState%Grid%PEDGE%Val(I,J,LBOTTOM  ) ) / 2.0_hp
       P2 = ( HcoState%Grid%PEDGE%Val(I,J,LBOTTOM  ) &
            + HcoState%Grid%PEDGE%Val(I,J,LBOTTOM+1) ) / 2.0_hp

       ! Temperature [K] at the centers of grid
       ! boxes (I,J,LBOTTOM-1) and (I,J,LBOTTOM)
       T1   = ExtState%TK%Arr%Val(I,J,LBOTTOM-1)
       T2   = ExtState%TK%Arr%Val(I,J,LBOTTOM  )
       
       ! DZ is the height [m] from the center of box (I,J,LCHARGE-1)
       ! to the negative charge layer.  It may be found in either
       ! the (LCHARGE)th sigma layer or the (LCHARGE-1)th layer.
       ! We use the hypsometric eqn to find the distance between
       ! the center of (LCHARGE)th and (LCHARGE-1)th boxes, then
       ! assume a linear temp distribution to scale between the two.
       DLNP = LOG( P1 / P2 ) / ( T1 - T2 ) * ( T1 - T_NEG_BOT )
       DZ   = HcoState%Phys%Rdg0 * ( ( T1 + T2 ) / 2d0 ) * DLNP

       ! Pressure [Pa] at the bottom edge of box (I,J,LBOTTOM),
       ! or, equivalently, the top edge of box (I,J,BOTTOM-1).
       P3   = HcoState%Grid%PEDGE%Val( I, J, LBOTTOM )

       ! Height [m] from the center of grid box (I,J,LBOTTOM-1) 
       ! to the top edge of grid box (I,J,LBOTTOM-1)
       ZUP  = HcoState%Phys%Rdg0 * T1 * LOG( P1 / P3 )

       !-----------------------------------------------------------
       ! (3b) HBOTTOM is the height of the 0 C layer above the 
       ! bottom edge of box (I,J,LBOTTOM).  
       ! 
       ! If DZ < ZUP, then DZ is in grid box (I,J,LBOTTOM-1);
       ! therefore subtract 1 from LBOTTOM and compute HBOTTOM
       ! accordingly.
       !
       ! In this case, please note that BXHEIGHT(I,J,LBOTTOM)-ZUP 
       ! is the distance from the bottom edge of the grid box to 
       ! the center of the newly defined (LBOTTOM)th layer.
       !-----------------------------------------------------------
       IF ( DZ >= ZUP ) THEN
          HBOTTOM = DZ - ZUP
       ELSE
          LBOTTOM = LBOTTOM - 1
          HBOTTOM = (HcoState%Grid%BXHEIGHT_M%Val(I,J,LBOTTOM) - ZUP) + DZ
       ENDIF
  
       ! Cold cloud thickness is difference of cloud top 
       ! height (H0) and the height to the bottom.
       CC = H0 - SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:LBOTTOM-1) ) - &
                     HBOTTOM 

       !===========================================================
       ! (4) COMPUTE IC/CG FLASH_RATIO FROM COLD-CLOUD DEPTH
       !
       ! This is necessary as an input for the MFLUX and PRECON
       ! parameterizations, as well as for determining the fraction 
       ! of LNOX generated by either type of flash, and will
       ! eventually be used for separate vertical distributions
       ! when they become available.  (ltm, bmy, 12/11/06)
       !===========================================================

       ! Get Inter-Cloud/Cloud-Ground flash ratio [unitless]
       IC_CG_RATIO = GET_IC_CG_RATIO( CC )

       !===========================================================
       ! (5) COMPUTE LIGHTNOX FLASH RATES
       !
       ! Now that we have computed the the ratio of intra-cloud
       ! flashes to cloud-ground flashes, compute the lightnox
       ! flash rate via one of these parameterizations:
       !
       ! (a) Cloud top height (CTH)
       ! (b) Mass flux (MFLUX)
       ! (c) Convective Precpitation (PRECON)
       ! 
       ! (ltm, bmy, 5/10/06, 12/11/06)
       !===========================================================

       !--------------------------------------------------------
       ! (5a) CLOUD TOP HEIGHT PARAMETERIZATION (all met fields)
       !
       ! Based on Price & Rind [1992,1993,1994].
       !--------------------------------------------------------

       ! Get lightnox flash rate per minute and IC/CG ratio
       CALL FLASHES_CTH( I, J, H0, FLASHRATE, SFCTYPE ) 

       !===========================================================
       ! (6) COMPUTE TOTAL LNOx AND PARTITION INTO VERTICAL LAYERS
       ! 
       ! (6a) We convert FLASHRATE (computed above) to units of
       ! [flashes/6h] and store in the RATE variable.
       !
       ! We then multiply RATE by a scale factor based on 
       ! OTD/LIS observations.  This is necessary in order to make 
       ! sure that the lightnox flashes happen in the correct 
       ! locations as diagnosed by OTD/LIS satellite observations.  
       ! There are two redistribution options:
       !
       !   (1) Apply regional scale factors based on OTD/LIS
       !        observations (method of L. Jourdain et al)
       !
       !   (2) Apply box-by-box scale scale factors based on
       !        OTD/LIS observations (method of B. Sauvage)
       !
       ! NOTE: As of 3/27/07, only method (1) is implemented.
       ! 
       ! (6b) We then compute X, which is the ratio
       !   [cloud-ground flashes / total flashes].
       !
       ! The amount of lightnox released will depend whether we
       ! are in the tropics or in mid-latitudes.
       !
       !
       ! (6c) LIGHTNOX NOx EMISSIONS IN THE TROPICS:
       ! ----------------------------------------------------------
       ! N. American / S. American     tropics: lat <= 23 N 
       ! African / Oceanian / Eurasian tropics: lat <= 35 N
       !
       ! The lightnox NOx released in the inter-cloud (IC) and 
       ! cloud-ground (CG) paths are given by:
       !
       !   TOTAL_IC = RFLASH_TROPIC * RATE * (1-X) * Z_IC
       !   TOTAL_CG = RFLASH_TROPIC * RATE * (  X) * Z_CG
       !
       ! where:
       !   RFLASH_TROPIC = # of NOx molecules released per flash 
       !                    per meter (same as in previous code)
       !   RATE          = lightnox flashes / 6h computed above
       !   Z_IC          = IC pathway in meters (from the negative    
       !                    cloud layer to the cloud top)
       !   Z_CG          = CG pathway in meters (from the negative
       !                    cloud layer to the ground surface)
       !   
       ! We also apply a top-down final global scaling factor, 
       ! calculated by previously bringing total global LNOx to 
       ! 6 Tg N/yr  (2x2.5: 0.3683, 4x5: 0.8996).  In 2004, the 
       ! tropics-only contribution to LNOx was 4.5379 Tg N.
       !
       ! 
       ! (6d) LIGHTING NOx EMISSIONS AT MIDLATITUDES:
       ! ----------------------------------------------------------
       ! N. American midlatitudes : lat > 23N
       ! Eurasian    midlatitudes : lat > 35N
       !
       ! The lightnox NOx released at midlatitudes is independent
       ! of path length.  Thus:
       !
       !   TOTAL_IC = RFLASH_MIDLAT * RATE * (1-X) * MID_LAT_SCALE
       !   TOTAL_CG = RFLASH_MIDLAT * RATE *    X  * MID_LAT_SCALE
       !
       ! where 
       !   RFLASH_MIDLAT = # of NOx molecules released per flash 
       !                    per meter (based on 500 mol/flash)
       !   RATE          = lightnox flashes / 6h computed above
       !   Z_IC          = IC pathway in meters (from the negative    
       !                    cloud layer to the cloud top)
       !   Z_CG          = CG pathway in meters (from the negative
       !                    cloud layer to the ground surface)
       !
       ! We now emit at the Northern Mid-latitudes using an RFLASH
       ! value of 500 mol/flash.  This is independent of path 
       ! length.  
       !
       ! NOTE: The OTD-LIS local redistribution method was expanded
       ! upon from Sauvage et al, 2007, ACP.
       ! http://www.atmos-chem-phys.net/7/815/2007/acp-7-815-2007.pdf
       !
       ! (6e) The total lightnox is the sum of IC+CG paths:
       !   TOTAL = TOTAL_IC + TOTAL_CG
       !
       ! (6g) We then partition the NOx into each of the vertical 
       ! grid boxes within the column with a Ken Pickering PDF
       ! (see comments below).
       !
       ! (ltm, rch, bmy, 5/10/06, 3/27/07)
       !===========================================================

       !-----------------------------------------------------------
       ! (6a) Compute flash rate and apply OTD/LIS redistribution
       !-----------------------------------------------------------

       ! Convert [flashes/min] to [flashes/6h]
       RATE     = FLASHRATE * 360.0d0

       ! Apply regional or local OTD-LIS redistribution so that the 
       ! flashes occur in the right place. 
       RATE = RATE * REDIST

       ! Apply scaling factor to make sure annual average flash rate 
       ! equals that of the climatology. (ltm, 09/24/07)
       RATE = RATE * Inst%OTD_LIS_SCALE

       END IF ! LBOTTOM >  1
       END IF ! LTOP    >= LCHARGE 
       END IF ! LTOP    >  0 
       END IF ! LCHARGE >  1 

       ! Do not allow flash density to become unrealistically high
       ! Globally limit the flash rate to its highest observed value
       ! of 4.2e-3 flashes / km2 / s from the ENTLN global product
       ! (ltm, mps, 8/9/18)
       IF ( ( RATE / 21600. / A_KM2 ) > 0.004177159 ) THEN
          RATE = 0.004177159 * 21600. * A_KM2
       END IF
       
       !-----------------------------------------------------------
       ! Eventually overwrite with values determined from imported
       ! flash rate 
       !-----------------------------------------------------------
       IF ( RATE_LFR > RATE ) THEN
          RATE        = RATE_LFR
          H0          = H0_LFR
          LTOP        = LTOP_LFR
          IC_CG_RATIO = IC_CG_RATIO_LFR
       END IF

       ! Error check: continue if flash rate is zero
       IF ( RATE <= 0.0 ) CYCLE

       !-----------------------------------------------------------
       ! (6b) Compute cloud-ground/total flash ratio
       !-----------------------------------------------------------

       ! Ratio of cloud-to-ground flashes to total # of flashes
       X    = 1d0 / ( 1d0 + IC_CG_RATIO )

       ! Compute LNOx emissions for tropics or midlats 
       IF ( XMID > EAST_WEST_DIV ) THEN 
          
          !--------------------------------------------------------
          ! (6c,6d) We are in EURASIA
          !--------------------------------------------------------
          IF ( YMID > EAST_NS_DIV ) THEN 


             ! 6d: Eurasian Mid-Latitudes
             TOTAL_IC = RFLASH_MIDLAT * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_MIDLAT * RATE * X
          ELSE                           

             ! 6c: Eurasian Tropics
             TOTAL_IC = RFLASH_TROPIC * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_TROPIC * RATE * X

          ENDIF

       ELSE   
          
          !--------------------------------------------------------
          ! (6c,6d) We are in the AMERICAS
          !--------------------------------------------------------
          IF ( YMID > WEST_NS_DIV ) THEN 

             ! 6d: American Mid-Latitudes
             TOTAL_IC = RFLASH_MIDLAT * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_MIDLAT * RATE * X

          ELSE                           

             ! 6c: American Tropics
             TOTAL_IC = RFLASH_TROPIC * RATE * ( 1d0 - X )
             TOTAL_CG = RFLASH_TROPIC * RATE * X

          ENDIF
       ENDIF

       !-----------------------------------------------------------
       ! (6e) Compute total lightnox
       !-----------------------------------------------------------

       ! Sum of IC + CG [molec/6h]
       TOTAL = TOTAL_IC + TOTAL_CG

       !-----------------------------------------------------------
       ! Store flash rates [flashes/min/km2] for diagnostics
       !-----------------------------------------------------------
       IF ( Inst%DoDiagn .AND. RATE > 0d0 ) THEN

          ! LightNOX flashes per minute per km2
          RATE_SAVE   = RATE / A_KM2 / 360d0

          ! Store total, IC, and CG flash rates in AD56
          DIAGN(I,J,1) = RATE_SAVE
          DIAGN(I,J,3) = RATE_SAVE * X 
          DIAGN(I,J,2) = H0 * 1d-3

       ENDIF

       !-----------------------------------------------------------
       ! (6g) LIGHTDIST computes the lightnox distribution from 
       ! the ground to the convective cloud top using cumulative
       ! distribution functions for ocean flashes, tropical land
       ! flashes, and non-tropical land flashes, as specified by
       ! Lesley Ott [JGR, 2010]
       !-----------------------------------------------------------

       ! If there's lightnox w/in the column ...
       IF ( TOTAL > 0d0 ) THEN

          ! Partition the column total NOx [molec/6h] from lightnox 
          ! into the vertical using Pickering PDF functions
          CALL LIGHTDIST( I, J, LTOP, H0, YMID, TOTAL, VERTPROF, &
                          ExtState, HcoState, SFCTYPE, cMt, MTYPE, Inst )

          ! Add vertically partitioned NOx into SLBASE array
          DO L = 1, HcoState%NZ
             Inst%SLBASE(I,J,L) = Inst%SLBASE(I,J,L) + VERTPROF(L) 

             ! No lightnox emissions in the stratosphere (cdh, 4/25/2013)
             IF ( HcoState%Grid%PEDGE%Val(I,J,L) < TROPP ) THEN 
                Inst%SLBASE(I,J,L) = 0.0_hp

             ELSE
                ! Convert to kg/m2/s
                ! SLBASE(I,J,L) has units [molec NOx/6h/box], convert units:
                ! [molec/6h/box] * [6h/21600s] * [area/AREA_M2 m2] /
                ! [MW/(g/mol)] / [Avgrd/(molec/mol)] * [1kg/1000g] = [kg/m2/s]
                Inst%SLBASE(I,J,L) = Inst%SLBASE(I,J,L)                       &
                              / (21600.d0*HcoState%Grid%AREA_M2%Val(I,J))     &
                              * HcoState%Spc(Inst%IDTNO)%EmMW_g               &
                              / HcoState%Phys%Avgdr / 1000.0d0
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
    CALL HCOX_SCALE ( am_I_Root, HcoState, Inst%SLBASE, &
                      TRIM(Inst%SpcScalFldNme(1)), RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! Eventually add diagnostics 
    !-----------------------------------------------------------------

    ! Eventually add individual diagnostics. These go by names!
    IF ( Inst%DoDiagn ) THEN
       ! Total flash rates
       Arr2D   => DIAGN(:,:,1)
        CALL Diagn_Update( am_I_Root, HcoState,               &
                          cName='LIGHTNING_TOTAL_FLASHRATE',  &
                          ExtNr=Inst%ExtNr, Array2D=Arr2D, RC=RC         ) 

       IF ( RC /= HCO_SUCCESS ) RETURN 
       Arr2D => NULL() 
   
       ! Intracloud flash rates 
       Arr2D     => DIAGN(:,:,2)
       CALL Diagn_Update( am_I_Root, HcoState,                    &
                          cName='LIGHTNING_INTRACLOUD_FLASHRATE', &
                          ExtNr=Inst%ExtNr, Array2D=Arr2D, RC=RC         ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
       Arr2D => NULL() 
   
       ! Cloud to ground flash rates
       Arr2D     => DIAGN(:,:,3)
       CALL Diagn_Update( am_I_Root, HcoState,                     &
                          cName='LIGHTNING_CLOUDGROUND_FLASHRATE', &
                          ExtNr=Inst%ExtNr, Array2D=Arr2D, RC=RC         ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
       Arr2D => NULL() 

       ! Cloud top height
       Arr2D     => TOPDIAGN(:,:)
       CALL Diagn_Update( am_I_Root,   HcoState,       &
                          cName='LIGHTNING_CLOUD_TOP', &
                          ExtNr=Inst%ExtNr, Array2D=Arr2D, RC=RC         ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
       Arr2D => NULL() 

    ENDIF

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
!  column lightnox NOx into the GEOS-Chem vertical layers. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LightDist( I, J, LTOP, H0, XLAT, TOTAL, VERTPROF, &
                        ExtState, HcoState, SFCTYPE, cMt, MTYPE, Inst )
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
    INTEGER,         INTENT(IN)  :: cMt        ! Current month 
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
!        w/in "lightnox.f".  Now read the "light_dist.dat.geos3" file for 
!        GEOS-3 directly from the DATA_DIR/lightnox_NOx_200203/ subdirectory.
!        Now read the "light_dist.dat" file for GEOS-1, GEOS-STRAT directly 
!        from the DATA_DIR/lightnox_NOx_200203/ subdirectory.  Added 
!        descriptive comment header.  Now trap I/O errors across all 
!        platforms with subroutine "ioerror.f".  Updated comments, cosmetic 
!        changes.  Redimension FRAC(NNLIGHT) to FRAC(LLPAR). (bmy, 4/2/02)
!  (3 ) Deleted obsolete code from April 2002.  Now reference IU_FILE and
!        IOERROR from "file_mod.f".  Now use IU_FILE instead of IUNIT as the
!        file unit number. (bmy, 6/27/02)
!  (4 ) Now reference BXHEIGHT from "dao_mod.f" (bmy, 9/18/02)
!  (5 ) Bug fix: add GEOS_4 to the #if block (bmy, 3/4/04)
!  (6 ) Now bundled into "lightnox_mod.f".  CDF's are now read w/in
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
    ! assign a flag describing the type of lightnox:
    !
    !   MTYPE = 1: ocean lightnox
    !   MTYPE = 2: tropical continental lightnox
    !   MTYPE = 3: midlatitude continental lightnox 
    !   MTYPE = 4: subtropical lightnox
    !             
    ! (ltm, bmy, 1/25/11)
    !=================================================================

    ! Assign profile kind to grid box, following Allen et al. 
    ! [JGR, 2010] (ltm, 1/25,11)
!    MONTH = cMt

    SELECT CASE (cMt)

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
    ! Use the CDF for this type of lightnox to partition the total
    ! column lightnox into the layers
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
      
    ! Partition lightnox NOx by layer into VERTPROF
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
! !IROUTINE: Flashes_CTH
!
! !DESCRIPTION: Subroutine Flashes\_CTH determines the rate of lightnox 
!  flashes per minute based on the height of convective cloud tops, and the 
!  intra-cloud to cloud-ground strike ratio.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Flashes_CTH( I, J, HEIGHT, FLASHRATE, SFCTYPE ) 
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN)  :: I           ! Longitude index
    INTEGER, INTENT(IN)  :: J           ! Latitude index
    REAL*8,  INTENT(IN)  :: HEIGHT      ! Height of conv cloud top [m]
    INTEGER, INTENT(IN)  :: SFCTYPE     ! Surface type 
!
! !OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(OUT) :: FLASHRATE   ! LightNOX flash rate [flashes/min]
!
! !REVISION HISTORY: 
!  10 May 2006 - L. Murray - Initial version
!  (1  ) Subroutine renamed from FLASHES (ltm, bmy, 5/10/06)
!  (2  ) Remove CCTHICK, IC_CG_RATIO as arguments.  Remove computation of
!         IC_CG_RATIO and move that to GET_IC_CG_RATIO. (ltm, bmy, 12/11/06)
!  (3  ) Remove the near-land formulation (i.e. use function IS_LAND 
!         instead of IS_NEAR).(ltm, bmy, 9/24/07)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  22 Oct 2013 - C. Keller   - Now a HEMCO extension.
!EOP
!------------------------------------------------------------------------------
!BOC
    !================================================================
    ! FLASHES_CTH begins here!
    !
    ! COMPUTE LIGHTNOX FLASH RATE / MINUTE
    !
    ! Price & Rind (1992) give the following parameterizations for
    ! lightnox flash rates as a function of convective cloud top
    ! height [km]:
    !
    !    FLAND  = 3.44e-5 * ( CLDTOP HEIGHT [km] ^ 4.9  )
    !    FOCEAN = 6.4e-4  * ( CLDTOP HEIGHT [km] ^ 1.73 )
    !
    ! LightNOX will therefore occur much more often on land.  It 
    ! goes as approx. the 5th power of height, as opposed to approx. 
    ! the 2nd power of height over oceans.
    !
    ! We suppress lightnox where the surface is mostly ice.  
    !
    ! (ltm, bmy, 5/10/06, 12/11/06)
    !================================================================

    ! Test for land type
    IF ( SFCTYPE == 0 ) THEN

       ! Flashes/min over land boxes
       FLASHRATE   = 3.44d-5 * ( ( HEIGHT * 1d-3 )**4.9d0  )

    ELSE IF ( SFCTYPE == 1 ) THEN

       ! Flahes/min over water
       FLASHRATE   = 6.4d-4  * ( ( HEIGHT * 1d-3 )**1.73d0 )

    ELSE IF ( SFCTYPE == 2 ) THEN

       ! Suppress lightnox over snow/ice
       FLASHRATE   = 0d0

    ENDIF

  END SUBROUTINE Flashes_CTH
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_IC_CG_Ratio
!
! !DESCRIPTION: Function Get\_IC\_CG\_Ratio calculates the Intra-Cloud (IC) 
!  and Cloud-to-Ground (CG) lightnox flash ratio based on the method of 
!  Price and Rind 1993, which is calculated from the cold-cloud depth 
!  (CCTHICK).
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_IC_CG_ratio( CCTHICK ) RESULT( IC_CG_RATIO )
!
! !INPUT PARAMETERS: 
!
    REAL*8,  INTENT(IN) :: CCTHICK       ! Cold cloud thickness [m]
!
! !RETURN VALUE:
!
    REAL*8              :: IC_CG_RATIO   ! Intra-cloud/cloud-ground ratio
!
! !REVISION HISTORY: 
!  11 Dec 2006 - R. Yantosca - Initial version
!  (1 ) Split off from FLASHES_CTH, FLASHES_MFLUX, FLASHES_PRECON into this
!        separate function (ltm, bmy, 12/11/06)
!  (2 ) Bug fix for XLF compiler (morin, bmy, 7/8/09)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  22 Oct 2013 - C. Keller   - Now a HEMCO extension.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8 :: CC, F_CG

    !=================================================================
    ! GET_IC_CG_RATIO begins here!
    !
    ! COMPUTE INTRA-CLOUD / CLOUD-GROUND FLASH RATIO 
    !
    ! Price & Rind (1993) compute the ratio of Cloud-Ground 
    ! to Total Flashes by the parameterization:
    !
    ! For 5.5 < dz < 14:
    !
    !     f_CG = 1 / (( A*dz^4 + B*dz^3 + C*dz^2 + D*dz + E ) + 1 )
    !
    ! For dz > 14:
    !
    !     f_CG = 0.02
    !                                        
    ! Where:
    !
    ! (1) dz is the depth [km] of the cloud above the freezing 
    !     level.  The cold-cloud thickness (dz) is the depth of 
    !     the layer between the cloud top and the center of the 
    !     highest layer for which the temperature exceeds 273 K. 
    !     The cold-cloud thickness is set to 5.5 km at grid points 
    !     where it is less than 5.5 km.
    !
    ! (2) The polynomial coefficients are:
    !        A=0.021,  B=-0.648,  C=7.493,  D=-36.54,  E=63.09
    !
    ! 
    ! Note: f_IC = 1 - f_CG
    ! 
    ! And hence,
    !
    !     IC_CG_RATIO = ( 1 - f_CG ) / f_CG
    !
    !
    ! IC_CG_RATIO is passed back to routine the LIGHTNOX_NL, where
    ! it is passed to FLASHES_MFLUX and FLASHES_PRECON.  In these
    ! routines, the fraction of total lightnox flashes that are 
    ! cloud-ground (CG) flashes is computed by:
    ! 
    !     F_CG        = 1d0 / ( 1d0 + IC_CG_RATIO )
    !
    ! and the fraction of the total lightnox flashes that are
    ! intra-cloud (IC) flashes is computed by:
    !
    !     F_IC        = 1d0 - 1d0 / ( 1d0 + IC_CG_RATIO )
    !=====================================================================

    ! Convert cold cloud thickness from [m] to [km] (min value: 5.5 km)
    CC = MAX( CCTHICK * 1d-3, 5.5d0 )

    ! Compute cloud-ground flash ratio as described above
    IF ( CC > 14d0 ) THEN

       ! Constant value above 14 km
       F_CG = 0.02d0

    ELSE

       ! First create the polynomial expression
       F_CG = 63.09d0 + CC * ( -36.54d0  + &
                        CC * (   7.493d0 + &
                        CC * (  -0.648d0 + &
                        CC * (   0.021d0 ) ) ) )

       ! Then put it in the denominator
       F_CG = 1d0 / ( F_CG + 1d0 )
                  
    ENDIF

    ! Intra-Cloud / Cloud-Ground flash ratio
    IC_CG_RATIO = ( 1d0 - F_CG ) / F_CG

  END FUNCTION Get_IC_CG_Ratio
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_OTD_LIS_Scale
!
! !DESCRIPTION: Function GET\_OTD\_LIS\_SCALE returns a met-field dependent 
!  scale factor which is to be applied to the lightnox flash rate to bring 
!  the annual average flash rate to match that of the OTD-LIS climatology 
!  ( ~ 45.9 flashes/sec ). Computed by running the model over the 11-year 
!  OTD-LIS campaign window and comparing the average flash rates, or as 
!  many years as are available.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_OTD_LIS_Scale( am_I_Root, HcoState, BETA, RC ) 
!
! !USES:
!
    USE HCO_Clock_Mod, ONLY : HcoClock_Get
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN   )   :: am_I_Root  ! Root CPU?
    TYPE(HCO_State), POINTER :: HcoState   ! HEMCO state obj
    REAL*8,  INTENT(  OUT)   :: BETA       ! Scale factor
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT)   :: RC         ! Suc
!
! !REVISION HISTORY: 
!  24 Sep 2007 - L. Murray - Initial version
!  (1 ) Added MFLUX, PRECON scaling for GEOS-4.  Also write messages for met
!        field types/grids where scaling is not defined. (ltm, bmy, 11/29/07)
!  (2 ) Now use different divisor for local redist (ltm, bmy, 2/20/08)
!  (3 ) Now compute the proper scale factor for GEOS-5 0.5 x 0.666 grids
!        and the GEOS-3 1x1 nested NA grid (yxw, dan, ltm, bmy, 11/14/08)
!  (4 ) Added "quick fix" for reprocessed GEOS-5 met fields to be used when 
!        the IN_CLOUD_OD switch is turned on. (ltm, bmy, 2/18/09)
!  (5 ) Added "quick fix" for 2004, 2005, 2008 OTD/LIS (ltm, bmy, 4/29/09)
!  (6 ) Updated scale factors for GEOS-5 based on 4+ years of data.  Remove
!        temporary fixes. (bmy, 7/10/09)
!  (7 ) Modification for GEOS-4 1 x 1.25 grid (lok, ltm, bmy, 1/13/10)
!  (8 ) Reprocessed for error in CLDTOPS field; Updated for GEOS
!        5.1.0 vs. 5.2.0; MERRA added; (ltm, bmy, 1/25/11)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  02 Feb 2012 - R. Yantosca - Compute BETA for MERRA 2 x 2.5
!  02 Feb 2012 - R. Yantosca - Compute BETA for GEOS-5.7.x
!  22 Oct 2013 - C. Keller   - Now a HEMCO extension.
!  04 Nov 2014 - Y. X. Wang  - Define BETA, ANN_AVG_FLASHRATE for the
!                              GEOS-FP 025x03125 NESTED_CH grid
!  14 Jan 2015 - L. Murray   - Updated GEOS-FP files through Oct 2014
!  01 Apr 2015 - R. Yantosca - Cosmetic changes
!  01 Apr 2015 - R. Yantosca - Bug fix: GRID025x0325 should be GRID025x03125
!  01 Mar 2016 - L. Murray   - Add preliminary values for MERRA-2 4x5, NA, CH
!  19 Jul 2016 - L. Murray   - Add preliminary values for MERRA-2 2x2.5
!  24 Sep 2017 - L. Murray   - Removed legacy resolutions. Updated LIS/OTD
!                              HRMC climatology. Final global MERRA-2 values.
!                              Updated GEOS-FP and regional MERRA-2 values.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! lun of error log
    INTEGER :: cYr, cMt
    REAL*8  :: DY

    !=================================================================
    ! Define the average annual flash rate (flashes per second), as
    ! calculated from the LIS/OTD High Resolution Monthly Climatology
    ! (LISOTD_HRMC_V2.3.2015.hdf)
    ! 
    ! doi: 10.5067/LIS/LIS-OTD/DATA303
    ! 
    ! The climatology contains data from May 1995 through Dec 2014.
    ! Slight difference in global mean total when averaging over different
    ! GEOS-Chem horizontal resolutions.
    !=================================================================
#if   defined( GRID2x25      ) 
    REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 46.019893d0

#elif defined( GRID4x5       )
    REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 46.019893d0

#elif defined( GRID025x03125 ) && defined( NESTED_CH )
    REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 4.8334473d0    

#elif defined( GRID025x03125 ) && defined( NESTED_EU )
    REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 1.3384335d0

#elif defined( GRID025x03125 ) && defined( NESTED_NA )
    REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 6.4666451d0

#elif defined( GRID05x0625   ) && defined( NESTED_AS )
    REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 9.2583674d0

#elif defined( GRID05x0625   ) && defined( NESTED_EU )
    REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 1.7925457d0

#elif defined( GRID05x0625   ) && defined( NESTED_NA )
    REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 6.7011175d0
    
#endif

    !=================================================================
    ! GET_OTD_LIS_SCALE begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'Get_OTD_LIS_Scale (hcox_lightnox_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Extract current year and month
    CALL HcoClock_Get( am_I_Root, HcoState%Clock, cYYYY=cYr, cMM=cMt, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! The lightnox flash rate equations are sensitive to model resolution
    ! and convection scheme used in the data assimilation.
    ! We know from the LIS/OTD satellite products that the global annual
    ! average flash rate is 46 fl s-1. We determine a single scaling 
    ! factor, beta, to be applied uniformly 
    !
    ! beta =  ( Annual Average Flash Rate Observed ) / 
    !           ( Annual Average Flash Rate Unconstr Parameterization )
    !
    ! This is equivalent to modifying the first coefficient of the 
    ! Price and Rind [1992] formulation to get the right magnitude 
    ! for a given model framework.
    ! 
    ! Beta corresponds to beta in Murray et al. [2011]
    !
    ! (ltm, 1/25/11)

    ! Initialize
    BETA = 1d0

#if   defined( GEOS_FP ) && defined( GRID4x5 ) 

    !---------------------------------------
    ! GEOS-FP: 4 x 5 global simulation
    !---------------------------------------
    
    ! Constrained with simulated "climatology" for
    ! April 2012 - Jul 2017. Will need to be updated as more
    ! met fields become available (ltm, 2017-09-24).
    IF ( ( cYr .eq. 2017 .and. cMt .le. 7 ) .or. cYr .le. 2016 ) THEN
       BETA = ANN_AVG_FLASHRATE / 85.362449d0
    ENDIF

#elif defined( GEOS_FP ) && defined( GRID2x25 )

    !---------------------------------------
    ! GEOS-FP: 2 x 2.5 global simulation
    !---------------------------------------
    
    ! Constrained with simulated "climatology" for
    ! April 2012 - Jul 2017. Will need to be updated as more
    ! met fields become available (ltm, 2017-09-24).    
    IF ( ( cYr .eq. 2017 .and. cMt .le. 7 ) .or. cYr .le. 2016 ) THEN
       BETA = ANN_AVG_FLASHRATE / 269.13945d0
    ENDIF
    
#elif defined( GEOS_FP ) && defined( GRID025x03125 ) && defined( NESTED_CH )

    !---------------------------------------
    ! GEOS-FP: Nested China simulation
    !---------------------------------------

    ! Constrained with simulated "climatology" for
    ! April 2012 - Jul 2017. Will need to be updated as more
    ! met fields become available (ltm, 2017-09-28).
    IF ( ( cYr .eq. 2017 .and. cMt .le. 7 ) .or. cYr .le. 2016 ) THEN
       BETA = ANN_AVG_FLASHRATE / 1030.6438d0
    ENDIF

#elif defined( GEOS_FP ) && defined( GRID025x03125 ) && defined( NESTED_EU )

    !---------------------------------------
    ! GEOS-FP: Nested Europe simulation
    !---------------------------------------

    ! Constrained with simulated "climatology" for
    ! April 2012 - Jul 2017. Will need to be updated as more
    ! met fields become available (ltm, 2017-09-28).
    IF ( ( cYr .eq. 2017 .and. cMt .le. 7 ) .or. cYr .le. 2016 ) THEN
       BETA = ANN_AVG_FLASHRATE / 90.403359d0
    ENDIF

#elif defined( GEOS_FP ) && defined( GRID025x03125 ) && defined( NESTED_NA )

    !---------------------------------------
    ! GEOS-FP: Nested North America simulation
    !---------------------------------------

    ! Constrained with simulated "climatology" for
    ! April 2012 - Jan 2017. Will need to be updated as more
    ! met fields become available (ltm, 2017-09-28).
    IF ( ( cYr .eq. 2017 .and. cMt .le. 7 ) .or. cYr .le. 2016 ) THEN    
       BETA = ANN_AVG_FLASHRATE / 770.29078d0
    ENDIF

#elif defined( MERRA2 ) && defined( GRID4x5 )

    !---------------------------------------
    ! MERRA2: 4 x 5 global simulation
    !---------------------------------------

    ! Constrained with simulated "climatology" for
    ! full LIS/OTD observational period (May 1995-Dec 2014). 
    ! Does not need to be updated (ltm, 2017-09-24).
    BETA = ANN_AVG_FLASHRATE / 102.38173d0

#elif defined( MERRA2 ) && defined( GRID2x25 )

    !---------------------------------------
    ! MERRA2: 2 x 2.5 global simulation
    !---------------------------------------

    ! Constrained with simulated "climatology" for
    ! full LIS/OTD observational period (May 1995-Dec 2014). 
    ! Does not need to be updated (ltm, 2017-09-24).
    BETA = ANN_AVG_FLASHRATE / 322.83040d0

#elif defined( MERRA2 ) && defined( GRID05x0625  ) && defined( NESTED_AS )
    
    !---------------------------------------
    ! MERRA-2: Nested Asia simulation
    !---------------------------------------

    ! Constrained with simulated "climatology" for
    ! full LIS/OTD observational period (May 1995-Dec 2014). 
    ! Does not need to be updated (ltm, 2017-09-28).
    BETA = ANN_AVG_FLASHRATE / 1177.4952d0
    
#elif defined( MERRA2 ) && defined( GRID05x0625  ) && defined( NESTED_EU )

    !------------------------------------------
    ! MERRA-2: Nested Europe simulation
    !------------------------------------------

    ! Constrained with simulated "climatology" for
    ! full LIS/OTD observational period (May 1995-Dec 2014). 
    ! Does not need to be updated (ltm, 2017-09-23).
    BETA = ANN_AVG_FLASHRATE / 47.187111d0
    
#elif defined( MERRA2 ) && defined( GRID05x0625  ) && defined( NESTED_NA )

    !------------------------------------------
    ! MERRA-2: Nested North America simulation
    !------------------------------------------

    ! Constrained with simulated "climatology" for
    ! full LIS/OTD observational period (May 1995-Dec 2014). 
    ! Does not need to be updated (ltm, 2017-09-28).
    BETA = ANN_AVG_FLASHRATE / 305.06467d0
    
#endif

    ! If not defined yet, try to estimate it from grid spacing
    ! ckeller, 6/6/2017
    IF ( BETA == 1d0 ) THEN
       ! average latitude spacing
       DY = ABS(MAXVAL(HcoState%Grid%YMID%Val) - MINVAL(HcoState%Grid%YMID%Val)) / ( HcoState%NY - 1 )

       ! 0.125 degrees / C720
       IF ( DY < 0.175d0 ) THEN
          BETA = 1.4152d-3
       ! 0.25 degrees / C360
       ELSEIF ( DY < 0.375d0 ) THEN
          BETA = 7.024d-3
       ! 0.5 degrees / C180
       ELSEIF ( DY < 0.75d0 ) THEN
          BETA = 1.527d-2
       ! 1.0 degrees / C90
       ELSEIF ( DY < 1.5d0 ) THEN
          BETA = 0.10d0
       ! 2.0 degrees / C48
       ELSEIF ( DY < 3.0d0 ) THEN
          BETA = 0.355d0
       ENDIF
    ENDIF

    IF ( BETA == 1d0 ) THEN

       WRITE( *,* ) 'WARNING:'       
       WRITE( *,* ) ''
       WRITE( *,* ) 'Your model configuration has not had lightning NOx'       
       WRITE( *,* ) 'emissions processed, or you are outside the period'
       WRITE( *,* ) 'for which lightning has been constrained to observations.'
       WRITE( *,* ) 'As a precaution, the model has gracefully stopped.'
       WRITE( *,* ) ''
       WRITE( *,* ) 'Please contact Lee Murray (lee.murray@rochester.edu),'
       WRITE( *,* ) 'who can help you pepare the necessary modifications'
       WRITE( *,* ) 'and files. You may also find what you are looking for at'
       WRITE( *,* ) 'http://ees.rochester.edu/atmos/data.html'
       WRITE( *,* ) ''
       WRITE( *,* ) 'You may remove this error trap at your peril by either'
       WRITE( *,* ) 'commenting out the call to HCO_ERROR in'
       WRITE( *,* ) 'HEMCO/Extensions/hcox_lightnox_mod.F90, or by manually'
       WRITE( *,* ) 'setting BETA in the HEMCO configuration file to a'
       WRITE( *,* ) 'value other than 1.0 in the Lightning NOx settings:'
       WRITE( *,* )
       WRITE( *,* ) '# ExtNr ExtName            on/off Species'
       WRITE( *,* ) '103     LightNOx         : on     NO'
       WRITE( *,* ) '    --> OTD-LIS scaling  :        0.55'
       WRITE( *,* ) ''
       WRITE( *,* ) 'For a description of BETA, see Murray et al. [2012].'
       WRITE( *,* ) 'It is the scaling factor that scales the incorrect'
       WRITE( *,* ) 'total mean flash rate in the model to the correct'
       WRITE( *,* ) 'climatological global mean of 46 fl/s.'
       WRITE( *,* ) ''
       WRITE( *,* ) 'However, if you disable the error trap, be aware that'
       WRITE( *,* ) 'the magnitude and distribution of lightning NOx may'
       WRITE( *,* ) 'be wildly unrealistic. This is especially possible for the'
       WRITE( *,* ) 'regional simulations, especially GEOS-FP.'
       WRITE( *,* ) 'You should always check to make sure that you have'
       WRITE( *,* ) '~6 Tg N a-1 globally, or typical regional values.'
       
       CALL HCO_ERROR( HcoState%Config%Err, 'No beta - see information in standard output', RC )
       RETURN
 
    ENDIF

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE Get_OTD_LIS_Scale
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_LightNOx_Init
!
! !DESCRIPTION: Subroutine HCOX\_LIGHTNOX\_INIT allocates all module arrays.  
!  It also reads the lightnox CDF data from disk before the first lightnox 
!  timestep. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_LightNOx_Init( am_I_Root, HcoState, ExtName, ExtState, RC ) 
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
    LOGICAL,          INTENT(IN   )  :: am_I_Root
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
!        GET_FLASH_SCALE_PRECON depending on the type of lightnox param used.
!        Updated comments.  (ltm, bmy, 1/31/07)
!  (6 ) Removed near-land stuff.  Renamed from HCOX_LightNOX_Init_NL to
!        HCOX_LightNOX_Init.  Now allocate EMIS_LI_NOx. (ltm, bmy, 10/3/07)
!  (7 ) Also update location of PDF file to lightnox_NOx_200709 directory. 
!        (bmy, 1/24/08)
!  (8 ) Read in new Ott profiles from lightnox_NOx_201101. Remove
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
    LOGICAL                        :: FOUND
    INTEGER, ALLOCATABLE           :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
    CHARACTER(LEN=255)             :: MSG, LOC, FILENAME
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

    ! Create AeroCom instance for this simulation
    Inst => NULL()
    CALL InstCreate ( ExtNr, ExtState%LightNOx, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot create AeroCom instance', RC )
       RETURN
    ENDIF

    ! Read settings specified in configuration file
    ! Note: the specified strings have to match those in 
    !       the config. file!
    CALL GetExtOpt( HcoState%Config, ExtNr, 'OTD-LIS factors', &
                     OptValBool=Inst%LOTDLOC, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If OTD-LIS factor are not being used, make sure that the corresponding
    ! gridded data will be ignored (e.g. not read) by HEMCO.
    IF ( .NOT. Inst%LOTDLOC ) THEN
       CALL ReadList_Remove ( am_I_Root, HcoState, 'LIGHTNOX_OTDLIS', RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Note: the OTD-LIS scale factor will be determined during run time
    ! as it requires the current time information.

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
    IF ( am_I_Root ) THEN
       MSG = 'Use lightning NOx emissions (extension module)'
       CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
       WRITE(MSG,*) ' - Use species ', TRIM(SpcNames(1)), '->', Inst%IDTNO 
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - Use OTD-LIS factors from file? ', Inst%LOTDLOC 
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - Use GEOS-5 flash rates: ', Inst%LLFR 
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - Use scalar scale factor: ', Inst%SpcScalVal(1)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) ' - Use gridded scale field: ', TRIM(Inst%SpcScalFldNme(1))
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    !-----------------
    ! Allocate arrays
    !-----------------

    ! Allocate PROFILE (holds the CDF table)
    ALLOCATE( Inst%PROFILE( NNLIGHT, NLTYPE ), STAT=AS )
    IF( AS /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'PROFILE', RC )
       RETURN
    ENDIF
    Inst%PROFILE = 0d0

    ! Allocate SLBASE (holds NO emissins from lightning)
    ALLOCATE( Inst%SLBASE(HcoState%NX,HcoState%NY,HcoState%NZ), STAT=AS )
    IF( AS /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'SLBASE', RC )
       RETURN
    ENDIF
    Inst%SLBASE = 0d0

    ! Allocate SLBASE (holds NO emissins from lightning)
    IF ( Inst%LOTDLOC ) THEN 
       ALLOCATE( Inst%OTDLIS(HcoState%NX,HcoState%NY), STAT=AS )
       IF( AS /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, 'OTDLIS', RC )
          RETURN
       ENDIF
       Inst%OTDLIS = 0d0
    ENDIF

    !=======================================================================
    ! Obtain lightning CDF's from Ott et al [JGR, 2010]. (ltm, 1/25/11)
    !=======================================================================

    ! Get filename from configuration file
    CALL GetExtOpt( HcoState%Config, ExtNr, 'CDF table', &
                     OptValChar=FILENAME, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Call HEMCO parser to replace tokens such as $ROOT, $MET, or $RES.
    ! There shouldn't be any date token in there ($YYYY, etc.), so just
    ! provide some dummy variables here
    CALL HCO_CharParse( HcoState%Config, FILENAME, -999, -1, -1, -1, -1, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Echo info
    IF ( am_I_Root ) THEN
       WRITE( MSG, 100 ) TRIM( FILENAME )
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF
100 FORMAT( '     - INIT_LIGHTNOX: Reading ', a )

    ! Find a free file LUN
    IU_FILE = findFreeLUN()
      
    ! Open file containing lightnox PDF data
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
         
    ! Read NNLIGHT types of lightnox profiles
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
    ! Further HEMCO setup
    !=======================================================================

    ! Activate met. fields required by this module
    ExtState%TK%DoUse      = .TRUE.
    ExtState%TROPP%DoUse   = .TRUE.
    ExtState%CNV_MFC%DoUse = .TRUE.
    ExtState%CNV_FRC%DoUse = .TRUE.
    ExtState%ALBD%DoUse    = .TRUE.
    ExtState%WLI%DoUse     = .TRUE.
    ExtState%LFR%DoUse     = .TRUE.

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

    ! ----------------------------------------------------------------
    ! Type specific initialization statements follow below
    ! ----------------------------------------------------------------
    Inst%DoDiagn = .FALSE.

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
       IF ( ASSOCIATED( Inst%OTDLIS        ) ) DEALLOCATE ( Inst%OTDLIS        )

       IF ( ASSOCIATED( Inst%PROFILE       ) ) DEALLOCATE ( Inst%PROFILE       )
       IF ( ASSOCIATED( Inst%SLBASE        ) ) DEALLOCATE ( Inst%SLBASE        )
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
