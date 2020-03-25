!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_soilnox_mod.F90
!
! !DESCRIPTION: Module HCOX\_SoilNOx\_Mod contains routines to compute soil
!  NOx emissions.  We follow the implementation in GEOS-Chem by Hudman
!  et al 2012.
!\\
!\\
! !INTERFACE:
!
MODULE HCOX_SoilNOx_Mod
!
! !USES:
!
  USE HCO_ERROR_Mod
  USE HCO_CHARTOOLS_MOD
  USE HCO_DIAGN_Mod
  USE HCOX_TOOLS_MOD
  USE HCOX_State_Mod,     ONLY : Ext_State
  USE HCO_STATE_Mod,      ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_SoilNOx_Run
  PUBLIC :: HCOX_SoilNOx_Init
  PUBLIC :: HCOX_SoilNOx_Final
!
! !REMARKS:
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities.
!                                                                             .
!  Original codes from:
!    HARVARD ATMOSPHERIC CHEMISTRY MODELING GROUP
!    MODULE FOR SOIL NOX EMISSIONS
!    by Yuhang Wang, Gerry Gardner, and Prof. Daniel Jacob
!  Updated model code:
!    by  Rynda Hudman, Neil Moore, Randall Martin, and Bram Maasakkers
!                                                                             .
!  The soil NOx code has been updated from the original implementation
!  of Yienger & Levy [1995] from  Wang et al., [1998] as summarized below.
!                                                                             .
!  Old:
!  ENOx   = f( T, biome, w/d)  x Pulse(precip) x canopy uptake + FERT
!                                                                             .
!  New:
!  ENOx   = f( T, biome, WFPS, Fert)  x Pulse(dryspell) x canopy uptake
!                                                                             .
!  1 - Update moisture treatment: soil moisture as a continuous variable
!  using WFPS rather than discrete wet/dry states and purely exponential
!  T impact (impact = -1. Tg N/yr)
!                                                                             .
!  2 - Update to Fertilizer:  new fertilizer maps including chemical and
!  manure fertilizer from Potter et al., [2010] distributed using MODIS EVI
!  seasonality, online-N deposition as a fertilizer source, and N-fertilizer
!  source subject to T, WFPS, and pulsing like other N (impact = +1.3 Tg N/yr)
!                                                                             .
!  3- Update Pulsing Scheme: Yan et al., [2005] (shorter, stronger pulses)
!  (impact = +1. Tg N/yr). Also added restart file containing dry spell
!  information to properly account for dry spell length in continuing runs.
!                                                                             .
!  References:
!  ============================================================================
!  (1 ) Wang, Y., D.J. Jacob, and J.A. Logan,  Global simulation of
!        tropospheric O3-NOx-hydrocarbon chemistry, 1. Model formulation,
!        J. Geophys. Res., 103/D9, 10, 713-10,726, 1998.
!  (2 ) Yienger, J.J, and H. Levy, Empirical model of global soil-biogenic
!        NOx emissions, J. Geophys. Res., 100, D6, 11,447-11464, June 20, 1995.
!  (3 ) Yan, X., T. Ohara, and H. Akimoto, Statistical modeling of global
!        soil NOx emissions, Global Biogeochem. Cycles, 19, GB3019,
!        doi:10.1029/2004GB002276, 2005.
!  (4 ) Potter, P., Ramankutty, N., Bennett, E., and Donner, S.:
!        Characterizing the Spatial Patterns of Global Fertilizer Application
!        and Manure Production, Earth Interactions, 14, 1-22,
!        10.1175/2009EI288.1, 2010.
!  (5 ) Moore, N.E., Improving global bottom-up biogenic soil NOx inventories,
!        Master's Thesis, Dalhousie University, 2007.
!  (6 ) Hudman, R.C., N.E. Moore, A.K. Mebust, R.V. Martin, A.R. Russell,
!        L.C. Valin, and R.C Cohen, Steps toward a mechanistic model of global
!        soil nitric oxide emissions: implementation and space
!        based-constraints, Atmos. Chem. Phys., 12, 7779-7795,
!        doi:10.5194/acp-12-7779-2012, 2012.
!
! !REVISION HISTORY:
!
!  17 Aug 2009 - R. Yantosca     - Columnized and cleaned up
!  17 Aug 2009 - R. Yantosca     - Added ProTeX headers
!  31 Jan 2011 - R. Hudman       - Added new code12259.perceus-ucb0
!  31 Jan 2011 - R. Hudman       - Updated headers
!  29 Aug 2012 - J.D. Maasakkers - Implemented Jacob and Bakwin CRF
!  29 Aug 2012 - J.D. Maasakkers - Adapted code to work with new (online
!                                  regridded) landfraction, climate and
!                                  fertilizer data
!  29 Aug 2012 - J.D. Maasakkers - Removed all unused Wang et al. code
!                                  (comments)
!  04 Nov 2013 - C. Keller       - Moved all soil NOx routines into one
!                                  module. Now a HEMCO extension.
!  28 Jul 2014 - C. Keller       - Now allow DRYCOEFF to be read through
!                                  configuration file (as setting)
!  11 Dec 2014 - M. Yannetti     - Changed REAL*8 to REAL(hp)
!  14 Oct 2016 - C. Keller       - Now use HCO_EvalFld instead of HCO_GetPtr.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
  ! Derived type to hold MODIS land types
  TYPE MODL
     REAL(hp), POINTER             :: VAL              (:,:)
  ENDTYPE MODL

  TYPE :: MyInst
     INTEGER                        :: Instance
     INTEGER                        :: ExtNr          ! Extension number
     INTEGER                        :: IDTNO          ! NO tracer ID
     LOGICAL                        :: LFERTILIZERNOX ! Use fertilizer NOx?
     REAL(hp)                       :: FERT_SCALE     ! fertilizer scale factor

     ! Dry period length (from restart)
     REAL(sp), ALLOCATABLE          :: DRYPERIOD    (:,:  )

     ! Pulse factors (from restart)
     REAL(sp), ALLOCATABLE          :: PFACTOR      (:,:  )
     REAL(sp), ALLOCATABLE          :: GWET_PREV    (:,:  )

     ! Deposition reservoir (from restart)
     REAL(sp), ALLOCATABLE          :: DEP_RESERVOIR(:,:  )

     ! NOx in the canopy
     REAL(hp), ALLOCATABLE          :: CANOPYNOX        (:,:,:)

     ! MODIS landtype
     TYPE(MODL), POINTER            :: LANDTYPE         (:    ) => NULL()

     ! Soil fertilizer (kg/m3)
     REAL(hp), POINTER              :: SOILFERT         (:,:  ) => NULL()

     ! Fraction of arid and non-arid land
     REAL(hp), POINTER              :: CLIMARID         (:,:  ) => NULL()
     REAL(hp), POINTER              :: CLIMNARID        (:,:  ) => NULL()

     ! DRYCOEFF (if read from settings in configuration file)
     REAL(hp), POINTER              :: DRYCOEFF(:)

     ! Overall scale factor to be applied to total soil NOx emissions. Must
     ! be defined in the HEMCO configuration file as extension attribute
     ! 'Scaling_NO'
     REAL(sp),          ALLOCATABLE :: SpcScalVal(:)
     CHARACTER(LEN=61), ALLOCATABLE :: SpcScalFldNme(:)

     TYPE(MyInst), POINTER          :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to all instances
  TYPE(MyInst), POINTER             :: AllInst => NULL()
!
! !DEFINED PARAMETERS:
!
  ! # of MODIS/Koppen biome types
  INTEGER, PARAMETER               :: NBIOM = 24

  ! Max. # of allowed drycoeff vars
  INTEGER, PARAMETER               :: MaxDryCoeff = 50

  ! Canopy wind extinction coefficients
  ! (cf. Yienger & Levy [1995], Sec 5), now a function of the
  ! MODIS/KOPPEN biometype (J.D. Maasakkers)
  REAL(hp),  PARAMETER, PRIVATE :: SOILEXC(NBIOM) =                 (/ &
        0.10, 0.50, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 1.00,    &
        1.00, 1.00, 1.00, 2.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00,    &
        4.00, 2.00, 0.10, 2.00                                       /)

  ! Steinkamp and Lawrence, 2011 A values, wet biome coefficients
  ! for each of the 24 soil biomes [ng N/m2/s].
  REAL(hp),  PARAMETER, PRIVATE  :: A_BIOME(NBIOM) =                (/ &
        0.00, 0.00, 0.00, 0.00, 0.00, 0.06, 0.09, 0.09, 0.01, 0.84,    &
        0.84, 0.24, 0.42, 0.62, 0.03, 0.36, 0.36, 0.35, 1.66, 0.08,    &
        0.44, 0.57, 0.57, 0.57                                       /)

  ! "A" coefficients for converting surface temp to soil temp
  ! for each of the 24 soil biomes
  REAL(hp),  PARAMETER, PRIVATE :: SOILTA(NBIOM)  =                 (/ &
        0.00, 0.92, 0.00, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66, 0.66,    &
        0.66, 0.66, 0.66, 0.66, 0.84, 0.84, 0.84, 0.84, 0.84, 0.84,    &
        0.84, 1.03, 1.03, 1.03                                       /)

  ! "B" coefficients for converting surface temp to soil temp
  ! for each of the 24 soil biomes
  REAL(hp),  PARAMETER, PRIVATE :: SOILTB(NBIOM)  =                 (/ &
        0.00, 4.40, 0.00, 8.80, 8.80, 8.80, 8.80, 8.80, 8.80, 8.80,    &
        8.80, 8.80, 8.80, 8.80, 3.60, 3.60, 3.60, 3.60, 3.60, 3.60,    &
        3.60, 2.90, 2.90, 2.90                                       /)

  ! MODIS/Koppen resistance values
  INTEGER, PARAMETER, PRIVATE :: SNIMODIS(NBIOM) =                  (/ &
           1,    2,    3,    4,    5,    6,    7,    8,    9,   10,    &
          11,   12,   13,   14,   15,   16,   17,   18,   19,   20,    &
          21,   22,   23,   24                                       /)

  INTEGER, PARAMETER, PRIVATE :: SNIRI(NBIOM) =                     (/ &
        9999,  200, 9999, 9999, 9999, 9999,  200,  200,  200,  200,    &
         200,  200,  200,  200,  200,  200,  200,  400,  400,  200,    &
         200,  200, 9999,  200                                       /)

  INTEGER, PARAMETER, PRIVATE :: SNIRLU(NBIOM) =                    (/ &
        9999, 9000, 9999, 9999, 9999, 9999, 9000, 9000, 9000, 9000,    &
        9000, 9000, 9000, 9000, 9000, 1000, 9000, 9000, 9000, 9000,    &
        1000, 9000, 9999, 9000                                       /)

  INTEGER, PARAMETER, PRIVATE :: SNIRAC(NBIOM) =                    (/ &
           0,  300,    0,    0,    0,    0,  100,  100,  100,  100,    &
         100,  100,  100,  100, 2000, 2000, 2000, 2000, 2000, 2000,    &
        2000,  200,  100,  200                                       /)

  INTEGER, PARAMETER, PRIVATE :: SNIRGSS(NBIOM) =                   (/ &
           0,    0,  100, 1000,  100, 1000,  350,  350,  350,  350,    &
         350,  350,  350,  350,  500,  200,  500,  500,  500,  500,    &
         200,  150,  400,  150                                       /)

  INTEGER, PARAMETER, PRIVATE :: SNIRGSO(NBIOM) =                   (/ &
        2000, 1000, 3500,  400, 3500,  400,  200,  200,  200,  200,    &
         200,  200,  200,  200,  200,  200,  200,  200,  200,  200,    &
         200,  150,  300,  150                                       /)

  INTEGER, PARAMETER, PRIVATE :: SNIRCLS(NBIOM) =                   (/ &
       9999, 2500, 9999, 9999, 9999, 9999, 2000, 2000, 2000, 2000,     &
       2000, 2000, 2000, 2000, 2000, 9999, 2000, 2000, 2000, 2000,     &
       9999, 2000, 9999, 2000                                        /)

  INTEGER, PARAMETER, PRIVATE :: SNIRCLO(NBIOM) =                   (/ &
       9999, 1000, 1000, 9999, 1000, 9999, 1000, 1000, 1000, 1000,     &
       1000, 1000, 1000, 1000, 1000, 9999, 1000, 1000, 1000, 1000,     &
       9999, 1000, 9999, 1000                                        /)

  INTEGER, PARAMETER, PRIVATE :: SNIVSMAX(NBIOM) =                  (/ &
         10,  100,  100,   10,  100,   10,  100,  100,  100,  100,     &
        100,  100,  100,  100,  100,  100,  100,  100,  100,  100,     &
        100,  100,  100,  100                                        /)

!  ! Conversion factor from kg NO to ng N
!  REAL(hp),  PARAMETER, PRIVATE :: kgNO_to_ngN = 4.666d11 !(14/30 * 1e12)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_SoilNOx_Run
!
! !DESCRIPTION: Subroutine HcoX\_SoilNox\_Run is the driver routine to
! calculate ship NOx emissions for the current time step. Emissions in
! [kg/m2/s] are added to the emissions array of the passed
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_SoilNOx_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_Types_Mod,      ONLY : DiagnCont
    USE HCO_CLOCK_MOD,      ONLY : HcoClock_First
    USE HCO_CLOCK_MOD,      ONLY : HcoClock_Rewind
    USE HCO_FLuxArr_Mod,    ONLY : HCO_EmisAdd
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE HCO_Calc_Mod,       ONLY : HCO_EvalFld
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt
    USE HCO_ExtList_Mod,    ONLY : HCO_GetOpt
    USE HCO_Restart_Mod,    ONLY : HCO_RestartGet
    USE HCO_Restart_Mod,    ONLY : HCO_RestartWrite
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER        :: ExtState    ! Module options

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState   ! Output obj
    INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  05 Nov 2013 - C. Keller - Initial Version
!  08 May 2015 - C. Keller - Now read/write restart variables from here to
!                            accomodate replay runs in GEOS-5.
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  29 Mar 2017 - M. Sulprizio- Read DEP_RESERVOIR_DEFAULT field from file for
!                              use when when DEP_RESERVOIR is not found in the
!                              HEMCO restart file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                  :: I, J, N
    REAL(hp), TARGET         :: FLUX_2D(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET         :: DIAG   (HcoState%NX,HcoState%NY)
    REAL(hp), TARGET         :: Tmp2D  (HcoState%NX,HcoState%NY)
    REAL(sp)                 :: Def2D  (HcoState%NX,HcoState%NY)
    REAL(hp)                 :: FERTDIAG, DEP_FERT, SOILFRT
    REAL*4                   :: TSEMIS
    REAL(hp)                 :: UNITCONV, IJFLUX
    REAL(dp), ALLOCATABLE    :: VecDp(:)
    LOGICAL                  :: FIRST
    LOGICAL                  :: aIR, FOUND
    CHARACTER(LEN= 31)       :: DiagnName
    CHARACTER(LEN=255)       :: MSG, DMY
    TYPE(MyInst),    POINTER :: Inst

    ! For manual diagnostics
    LOGICAL, SAVE            :: DoDiagn = .FALSE.
    TYPE(DiagnCont), POINTER :: TmpCnt

    !=================================================================
    ! HCOX_SoilNOx_RUN begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_SoilNox_Run (hcox_soilnox_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return if extension disabled
    IF ( ExtState%SoilNOx < 0 ) RETURN

    ! Nullify
    Inst   => NULL()
    TmpCnt => NULL()

    ! Get Instance
    CALL InstGet ( ExtState%SoilNox, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find soil NOx instance Nr. ', ExtState%SoilNOx
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
    ENDIF

    ! Conversion factor from ng N to kg NO
    UNITCONV = 1.0e-12_hp / 14.0e+0_hp * HcoState%Spc(Inst%IDTNO)%EmMW_g

    !-----------------------------------------------------------------
    ! On first call, set pointers to all arrays needed by SoilNOx
    !-----------------------------------------------------------------
    FIRST = HcoClock_First ( HcoState%Clock, .TRUE. )

    !IF ( FIRST ) THEN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK1',  Inst%LANDTYPE(1)%VAL,  RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK2',  Inst%LANDTYPE(2)%VAL,  RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK3',  Inst%LANDTYPE(3)%VAL,  RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK4',  Inst%LANDTYPE(4)%VAL,  RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK5',  Inst%LANDTYPE(5)%VAL,  RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK6',  Inst%LANDTYPE(6)%VAL,  RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK7',  Inst%LANDTYPE(7)%VAL,  RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK8',  Inst%LANDTYPE(8)%VAL,  RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK9',  Inst%LANDTYPE(9)%VAL,  RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK10', Inst%LANDTYPE(10)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK11', Inst%LANDTYPE(11)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK12', Inst%LANDTYPE(12)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK13', Inst%LANDTYPE(13)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK14', Inst%LANDTYPE(14)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK15', Inst%LANDTYPE(15)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK16', Inst%LANDTYPE(16)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK17', Inst%LANDTYPE(17)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK18', Inst%LANDTYPE(18)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK19', Inst%LANDTYPE(19)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK20', Inst%LANDTYPE(20)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK21', Inst%LANDTYPE(21)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK22', Inst%LANDTYPE(22)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK23', Inst%LANDTYPE(23)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_LANDK24', Inst%LANDTYPE(24)%VAL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_FERT',    Inst%SOILFERT,         RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_ARID',    Inst%CLIMARID,         RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL HCO_EvalFld( HcoState, 'SOILNOX_NONARID', Inst%CLIMNARID,        RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

    IF ( FIRST ) THEN
       ! Check if ExtState variables DRYCOEFF is defined. Otherwise, try to
       ! read it from settings.
       IF ( .NOT. ASSOCIATED(ExtState%DRYCOEFF) ) THEN
          CALL GetExtOpt( HcoState%Config, Inst%ExtNr, 'DRYCOEFF', &
                           OptValChar=DMY, FOUND=FOUND, RC=RC )
          IF ( .NOT. FOUND ) THEN
             CALL HCO_ERROR( HcoState%Config%Err, 'DRYCOEFF not defined', RC )
             RETURN
          ENDIF
          ALLOCATE(VecDp(MaxDryCoeff))
          CALL HCO_CharSplit( DMY, HCO_GetOpt(HcoState%Config%ExtList,'Separator'), &
                              HCO_GetOpt(HcoState%Config%ExtList,'Wildcard'), VecDp, N, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          ALLOCATE(Inst%DRYCOEFF(N))
          Inst%DRYCOEFF(1:N) = VecDp(1:N)
          ExtState%DRYCOEFF => Inst%DRYCOEFF
          DEALLOCATE(VecDp)
       ENDIF

       ! Check if we need to write manual fertilizer NO diagnostics
       DiagnName = 'EmisNO_Fert'
       CALL DiagnCont_Find ( HcoState%Diagn, -1, -1, -1, -1, -1, &
                             DiagnName, 0, DoDiagn, TmpCnt )
       TmpCnt => NULL()
    ENDIF

    !---------------------------------------------------------------
    ! Fill restart variables. Restart variables can be specified
    ! in the HEMCO configuration file. In an ESMF environment, they
    ! can also be defined as internal state fields. The internal
    ! state fields take precedence over fields read through the
    ! HEMCO interface.
    ! The restart variables must be read on the first time step or
    ! after the rewinding the clock.
    !---------------------------------------------------------------

    IF ( FIRST .OR. HcoClock_Rewind( HcoState%Clock, .TRUE. ) ) THEN

       ! DEP_RESERVOIR. Read in kg NO/m3
       CALL HCO_EvalFld( HcoState, 'DEP_RESERVOIR_DEFAULT', &
                         Tmp2D, RC, FOUND=FOUND )
       IF ( FOUND ) THEN
          Def2D = Tmp2D
       ELSE
          Def2D = 1.0e-4_sp
       ENDIF
       CALL HCO_RestartGet( HcoState, 'DEP_RESERVOIR', &
                            Inst%DEP_RESERVOIR, RC, Def2D=Def2D )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! GWET_PREV [unitless]
       CALL HCO_RestartGet( HcoState, 'GWET_PREV', &
                            Inst%GWET_PREV, RC, FILLED=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( .NOT. FOUND ) THEN
          Inst%GWET_PREV = 0.0_sp
          IF ( HcoState%amIRoot ) THEN
             MSG = 'Cannot find GWET_PREV restart variable - initialized to 0.0!'
             CALL HCO_WARNING(HcoState%Config%Err,MSG,RC)
          ENDIF
       ENDIF

       ! PFACTOR [unitless]
       CALL HCO_RestartGet( HcoState, 'PFACTOR', &
                            Inst%PFACTOR, RC, FILLED=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( .NOT. FOUND ) THEN
          Inst%PFACTOR = 1.0_sp
          IF ( HcoState%amIRoot ) THEN
             MSG = 'Cannot find PFACTOR restart variable - initialized to 1.0!'
             CALL HCO_WARNING(HcoState%Config%Err,MSG,RC)
          ENDIF
       ENDIF

       ! DRYPERIOD [unitless]
       CALL HCO_RestartGet( HcoState, 'DRYPERIOD', &
                            Inst%DRYPERIOD, RC, FILLED=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( .NOT. FOUND ) THEN
          Inst%DRYPERIOD = 0.0_sp
          IF ( HcoState%amIRoot ) THEN
             MSG = 'Cannot find DRYPERIOD restart variable - initialized to 0.0!'
             CALL HCO_WARNING(HcoState%Config%Err,MSG,RC)
          ENDIF
       ENDIF

    ENDIF

    !---------------------------------------------------------------
    ! Now need to call GET_CANOPY_NOX to break ugly dependency between
    ! drydep and soil NOx emissions. (bmy, 6/22/09)
    ! Now a function of the new MODIS/Koppen biome map (J.D. Maasakkers)
    CALL GET_CANOPY_NOX( HcoState, ExtState, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Init
    TSEMIS  = HcoState%TS_EMIS
    FLUX_2D = 0e+0_hp
    IF(DoDiagn) DIAG = 0e+0_hp

    ! Loop over each land grid-box, removed loop over landpoints
!$OMP PARALLEL DO                                               &
!$OMP DEFAULT( SHARED )                                         &
!$OMP PRIVATE( I, J, DEP_FERT, SOILFRT, FERTDIAG, IJFLUX )
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! Do not calculate soil NOx emissions if there is no soil in
       ! the gridbox
       IF ( Inst%LANDTYPE(1)%VAL(I,J) == 1.0_hp ) CYCLE

       ! Get Deposited Fertilizer DEP_FERT [kg NO/m2]
       CALL GET_DEP_N( I, J, ExtState, HcoState, Inst, DEP_FERT )

       ! Get N fertilizer reservoir associated with chemical and
       ! manure fertilizer [kg NO/m2]
       IF ( Inst%LFERTILIZERNOX ) THEN
          SOILFRT =  Inst%SOILFERT( I, J )
       ELSE
          SOILFRT = 0e+0_hp
       ENDIF

       ! Put in constraint if dry period gt 1 yr, keep at 1yr to
       ! avoid unrealistic pulse
       IF ( Inst%DRYPERIOD(I,J) > 8760e+0_sp ) &
          Inst%DRYPERIOD(I,J) = 8760e+0_sp

       ! Return NO emissions from soils [kg NO/m2/s]
       CALL SOIL_NOX_EMISSION( ExtState,  Inst,                &
                               TSEMIS, I, J,                   &
                               SOILFRT,                        &
                               Inst%GWET_PREV(I,J), Inst%DRYPERIOD(I,J), &
                               Inst%PFACTOR(I,J),   IJFLUX,         &
                               DEP_FERT,       FERTDIAG,       &
                               UNITCONV,                       &
                               Inst%CANOPYNOX(I,J,:)                 )

       ! Write out
       FLUX_2D(I,J) = MAX(IJFLUX,0.0_hp)
       IF (DoDiagn) DIAG(I,J) = FERTDIAG

    ENDDO !J
    ENDDO !I
!$OMP END PARALLEL DO

    !-----------------------------------------------------------------
    ! EVENTUALLY ADD SCALE FACTORS
    !-----------------------------------------------------------------

    ! Eventually apply species specific scale factor
    IF ( Inst%SpcScalVal(1) /= 1.0_sp ) THEN
       FLUX_2D = FLUX_2D * Inst%SpcScalVal(1)
    ENDIF

    ! Eventually apply spatiotemporal scale factors
    CALL HCOX_SCALE ( HcoState, FLUX_2D, TRIM(Inst%SpcScalFldNme(1)), RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS
    !-----------------------------------------------------------------

    ! Add flux to emission array
    CALL HCO_EmisAdd( HcoState, FLUX_2D, Inst%IDTNO, &
                      RC,       ExtNr=Inst%ExtNr )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error', RC )
       RETURN
    ENDIF

    ! 'EmisNO_Fert' is the fertilizer NO emissions.
    ! This is a manual diagnostic created in GeosCore/hcoi_gc_diagn_mod.F90.
    ! If an empty pointer (i.e. not associated) is passed to Diagn_Update,
    ! diagnostics are treated as zeros!
    IF ( DoDiagn ) THEN
       DiagnName = 'EmisNO_Fert'
       CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                          cName=TRIM(DiagnName), Array2D=DIAG, RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! Eventually copy internal values to ESMF internal state object
    ! ----------------------------------------------------------------

    ! DEP_RESERVOIR [kg/m3]
    CALL HCO_RestartWrite( HcoState, 'DEP_RESERVOIR', Inst%DEP_RESERVOIR, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! GWET_PREV [unitless]
    CALL HCO_RestartWrite( HcoState, 'GWET_PREV', Inst%GWET_PREV, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! PFACTOR [unitless]
    CALL HCO_RestartWrite( HcoState, 'PFACTOR', Inst%PFACTOR, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! DRYPERIOD [unitless]
    CALL HCO_RestartWrite( HcoState, 'DRYPERIOD', Inst%DRYPERIOD, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave w/ success
    Inst => NULL()
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_SoilNox_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_SoilNOx_Init
!
! !DESCRIPTION: Subroutine HcoX\_SoilNox\_Init initializes the HEMCO
! SOILNOX extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_SoilNOx_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod,        ONLY : GetExtNr, GetExtOpt
    USE HCO_ExtList_Mod,        ONLY : GetExtSpcVal
    USE HCO_STATE_MOD,          ONLY : HCO_GetExtHcoID
    USE HCO_Restart_Mod,        ONLY : HCO_RestartDefine
!
! !ARGUMENTS:
!
    TYPE(HCO_State),  POINTER        :: HcoState   ! Output obj
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName    ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState   ! Module options
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  05 Nov 2013 - C. Keller   - Initial Version
!  12 May 2015 - R. Yantosca - Cosmetic changes
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: ExtNr
    CHARACTER(LEN=255)             :: MSG, LOC
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
    INTEGER, ALLOCATABLE           :: HcoIDs(:)
    INTEGER                        :: nSpc, I, J, II, AS
    TYPE(MyInst), POINTER          :: Inst

    !=================================================================
    ! HCOX_SoilNOx_INIT begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_SoilNOx_Init (hcox_soilnox_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create instance
    Inst => NULL()
    CALL InstCreate ( ExtNr, ExtState%SoilNox, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot create soil NOx instance', RC )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------------
    ! Get species IDs and settings
    ! ----------------------------------------------------------------------

    ! Read settings specified in configuration file
    ! Note: the specified strings have to match those in
    !       the config. file!
    CALL GetExtOpt( HcoState%Config, ExtNr, 'Use fertilizer NOx', &
                     OptValBool=Inst%LFERTILIZERNOX, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get global scale factor
    Inst%FERT_SCALE = HCOX_SoilNOx_GetFertScale()

    ! Get HEMCO species IDs
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( nSpc /= 1 ) THEN
       MSG = 'Module soil NOx accepts only one species!'
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

    ! Verbose mode
    IF ( HcoState%amIRoot ) THEN
       MSG = 'Use soil NOx emissions (extension module)'
       CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )

       WRITE(MSG,*) '   - NOx species            : ', TRIM(SpcNames(1)), Inst%IDTNO
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) '   - NOx scale factor       : ', Inst%SpcScalVal(1)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) '   - NOx scale field        : ', TRIM(Inst%SpcScalFldNme(1))
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) '   - Use fertilizer NOx     : ', Inst%LFERTILIZERNOX
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) '   - Fertilizer scale factor: ', Inst%FERT_SCALE
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP2='-')
    ENDIF

    ! ----------------------------------------------------------------------
    ! Set module variables
    ! ----------------------------------------------------------------------

    ! horizontal dimensions
    I = HcoState%NX
    J = HcoState%NY

    ALLOCATE( Inst%DRYPERIOD( I, J ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'DRYPERIOD', RC )
       RETURN
    ENDIF

    ALLOCATE( Inst%PFACTOR( I, J ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'PFACTOR', RC )
       RETURN
    ENDIF

    ALLOCATE( Inst%GWET_PREV( I, J ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'GWET_PREV', RC )
       RETURN
    ENDIF

    ALLOCATE( Inst%DEP_RESERVOIR( I, J ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'DEP_RESERVOIR', RC )
       RETURN
    ENDIF

    ALLOCATE( Inst%CANOPYNOX( I, J, NBIOM ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'CANOPYNOX', RC )
       RETURN
    ENDIF

    ! Reserve 24 pointers for land fractions for each Koppen category
    ALLOCATE ( Inst%LANDTYPE(NBIOM), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'LANDTYPE', RC )
       RETURN
    ENDIF
    DO II = 1,NBIOM
       ALLOCATE( Inst%LANDTYPE(II)%VAL( I, J ), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'LANDTYPE array', RC )
          RETURN
       ENDIF
       Inst%LANDTYPE(II)%Val = 0.0_hp
    ENDDO

    ALLOCATE ( Inst%SOILFERT  ( I, J ), &
               Inst%CLIMARID  ( I, J ), &
               Inst%CLIMNARID ( I, J ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'SOILFERT', RC )
       RETURN
    ENDIF
    Inst%SOILFERT  = 0.0_hp
    Inst%CLIMARID  = 0.0_hp
    Inst%CLIMNARID = 0.0_hp

    ! Zero arrays
    Inst%DRYPERIOD     = 0.0_sp
    Inst%PFACTOR       = 0.0_sp
    Inst%GWET_PREV     = 0.0_sp
    Inst%DEP_RESERVOIR = 0.0_sp
    Inst%CANOPYNOX     = 0e+0_hp

    ! Initialize pointers
    !Inst%CLIMARID  => NULL()
    !Inst%CLIMNARID => NULL()
    !Inst%SOILFERT  => NULL()

    ! ----------------------------------------------------------------------
    ! Set diagnostics
    ! ----------------------------------------------------------------------
    CALL HCO_RestartDefine( HcoState, 'PFACTOR', &
                            Inst%PFACTOR, '1',  RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_RestartDefine( HcoState, 'DRYPERIOD', &
                            Inst%DRYPERIOD, '1',  RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_RestartDefine( HcoState, 'GWET_PREV', &
                            Inst%GWET_PREV, '1',  RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_RestartDefine( HcoState, '   DEP_RESERVOIR', &
                            Inst%DEP_RESERVOIR, 'kg/m3', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ----------------------------------------------------------------------
    ! Set HEMCO extensions variables
    ! ----------------------------------------------------------------------

    ! Activate required met fields
    ExtState%T2M%DoUse       = .TRUE.
    ExtState%GWETTOP%DoUse   = .TRUE.
    ExtState%SUNCOS%DoUse    = .TRUE.
    ExtState%U10M%DoUse      = .TRUE.
    ExtState%V10M%DoUse      = .TRUE.
    ExtState%LAI%DoUse       = .TRUE.
    ExtState%ALBD%DoUse      = .TRUE.
    ExtState%RADSWG%DoUse    = .TRUE.
    ExtState%CLDFRC%DoUse    = .TRUE.

    ! Activate required deposition parameter
    ExtState%DRY_TOTN%DoUse = .TRUE.
    ExtState%WET_TOTN%DoUse = .TRUE.

    ! Leave w/ success
    Inst => NULL()
    IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_SoilNOx_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_SoilNOx_Final
!
! !DESCRIPTION: Subroutine HcoX\_SoilNOx\_Final finalizes the HEMCO
! SOILNOX extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_SoilNOx_Final( HcoState, ExtState, RC )
!
! !USES
!
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState      ! HEMCO State obj
    TYPE(Ext_State),  POINTER       :: ExtState      ! Extension state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  05 Nov 2013 - C. Keller - Initial Version
!
! !NOTES:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    !=================================================================
    ! HCOX_SoilNOx_FINAL begins here!
    !=================================================================
    CALL InstRemove ( ExtState )

  END SUBROUTINE HCOX_SoilNox_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_SoilNOx_GetFertScale
!
! !DESCRIPTION: Function HCOX\_SoilNOx\_GETFERTSCALE returns the scale factor
! applied to fertilizer NOx emissions.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCOX_SoilNOx_GetFertScale() RESULT ( FERT_SCALE )
!
! !ARGUMENTS:
!
    REAL(hp) :: FERT_SCALE
!
! !REMARKS:
!
! !REVISION HISTORY:
!  11 Dec 2013 - C. Keller   - Initial version
!  12 May 2015 - R. Yantosca - Bug fix: PGI expects routine name to end w/ ()
!!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    ! Scale factor so that fertilizer emission = 1.8 Tg N/yr
    ! (Stehfest and Bouwman, 2006)
    ! before canopy reduction
    FERT_SCALE = 0.0068
    ! Value calculated by running the 2x2.5 model
    ! For now, use this value for all resolutions since regular soil NOx
    ! emissions change with resolution as well (J.D. Maasakkers)

  END FUNCTION HCOX_SoilNOx_GetFertScale
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Soil_NOx_Emission
!
! !DESCRIPTION: Subroutine Soil\_NOx\_Emission computes the emission of soil and
!  fertilizer NOx for the GEOS-Chem model.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Soil_NOx_Emission( ExtState,  Inst, &
                                TS_EMIS,   I, J, &
                                SOILFRT,   &
                                GWET_PREV, DRYPERIOD, &
                                PFACTOR,   SOILNOx,   &
                                DEPN,      FERTDIAG,  &
                                UNITCONV,  R_CANOPY )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER :: ExtState      ! Module options
    TYPE(MyInst),    POINTER :: Inst          ! Instance object
    REAL*4,   INTENT(IN)  :: TS_EMIS          ! Emission timestep [s]
    INTEGER,  INTENT(IN)  :: I                ! grid box lon index
    INTEGER,  INTENT(IN)  :: J                ! grid box lat index
    REAL(hp), INTENT(IN)  :: DEPN             ! Dry Dep Fert term [kg/m2]
    REAL(hp), INTENT(IN)  :: SOILFRT          ! Fertilizer emissions [kg/m2]
    REAL(hp), INTENT(IN)  :: UNITCONV         ! ng N to kg NO

    !Input parameters for the canopy reduction factor
    REAL(hp), INTENT(IN)  :: R_CANOPY(:)      ! Resist. of canopy to NOx [1/s]
!
! !OUTPUT PARAMETERS:
!
    REAL(hp), INTENT(OUT) :: SOILNOx          ! Soil NOx emissions [kg/m2/s]
    REAL(sp), INTENT(OUT) :: GWET_PREV        ! Soil Moisture Prev timestep
    REAL(sp), INTENT(OUT) :: DRYPERIOD        ! Dry period length in hours
    REAL(sp), INTENT(OUT) :: PFACTOR          ! Pulsing Factor
    REAL(hp), INTENT(OUT) :: FERTDIAG         ! Fert emissions [kg/m2/s]
!
! !REMARKS:
!  R_CANOPY is computed in routine GET_CANOPY_NOX of "canopy_nox_mod.f".
!  This was originally in the GEOS-Chem dry deposition code, but was split
!  off in order to avoid an ugly code dependency between the dry deposition
!  and soil NOx codes.
!
!  As of v9-02, this module uses the MODIS/Koppen biome types instead
!  of the Olson land type / biome type, making it different from the original
!  dry deposition code (J.D. Maasakkers)
!
! !REVISION HISTORY:
!  17 Aug 2009 - R. Yantosca - Columnized and cleaned up
!  17 Aug 2009 - R. Yantosca - Added ProTeX headers
!  31 Jan 2011 - R. Hudman   - New Model added
!  23 Oct 2012 - M. Payer    - Now reference Headers/gigc_errcode_mod.F90
!  12 May 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: K
    REAL(hp)  :: BASE_TERM, CRF_TERM,  PULSE
    REAL(hp)  :: TC,        TEMP_TERM, WINDSQR
    REAL(hp)  :: WET_TERM,  A_FERT,    A_BIOM
    REAL(hp)  :: LAI,       SUNCOS,    GWET
    REAL(hp)  :: ARID,      NARID

    !=================================================================
    ! Initialize
    !=================================================================

    ! Initialize
    SOILNOX        = 0e+0_hp
    FERTDIAG       = 0e+0_hp

    ! Surface temperature [C]
    TC             = ExtState%T2M%Arr%Val(I,J) - 273.15e+0_hp

    ! Surface wind speed, squared
    WINDSQR        = ExtState%U10M%Arr%Val(I,J)**2 + &
                     ExtState%V10M%Arr%Val(I,J)**2

    ! Leaf area index
    LAI = ExtState%LAI%Arr%Val(I,J)

    ! Cosine of Solar Zenit Angle
    SUNCOS = ExtState%SUNCOS%Arr%Val(I,J)

    ! Top soil wetness [unitless]
    GWET = ExtState%GWETTOP%Arr%Val(I,J)

    !=================================================================
    ! Compute soil NOx emissions
    !=================================================================

    ! Cumulative multiplication factor (over baseline emissions)
    ! that accounts for soil pulsing
    PULSE = PULSING( GWET, TS_EMIS, GWET_PREV, PFACTOR, DRYPERIOD )

    ! ------Loop Over MODIS/Koppen  Landtypes
    DO K = 1, 24

       ! Temperature-dependent term of soil NOx emissions [unitless]
       ! Use GWET instead of climo wet/dry
       TEMP_TERM = SOILTEMP( K , TC, GWET)

       ! Soil moisture scaling of soil NOx emissions
       ARID     = Inst%CLIMARID(I,J)
       NARID    = Inst%CLIMNARID(I,J)
       WET_TERM = SOILWET( GWET , ARID, NARID )

       ! Fertilizer emission [kg/m2/s]
       A_FERT = FERTADD( SOILFRT , DEPN)

       ! Scale fertilizer emissions as specified
       ! (scale needed to force fert emiss of 1.8 Tg N/yr w/o canopy uptake)
       A_FERT = A_FERT * Inst%FERT_SCALE

       ! Canopy reduction factor
       CRF_TERM  = SOILCRF( K, LAI, R_CANOPY(K), WINDSQR, SUNCOS )

       ! Base emission. ng N/m2/s --> kg NO/m2/s
       A_BIOM = A_BIOME(K) * UNITCONV

       ! SOILNOX includes fertilizer
       SOILNOX   = (SOILNOX                            &
                 + ( A_BIOM + A_FERT )                 &
                 * ( TEMP_TERM * WET_TERM * PULSE )    &
                 * Inst%LANDTYPE(K)%VAL(I,J)           &
                 * ( 1.e+0_hp - CRF_TERM  )                 )

       ! FERTDIAG, only used for the fertilizer diagnostic
       FERTDIAG  = (FERTDIAG                           &
                 + ( A_FERT )                          &
                 * ( TEMP_TERM * WET_TERM * PULSE )    &
                 * Inst%LANDTYPE(K)%VAL(I,J)           &
                 * ( 1.e+0_hp - CRF_TERM  )                 )

    ENDDO

  END SUBROUTINE Soil_NOx_Emission
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Canopy_NOx
!
! !DESCRIPTION: Subroutine Get\_Canopy\_NOx computes the bulk surface
!  resistance of the canopy to NOx.  This computation was originally done
!  within legacy routine DEPVEL (in "drydep\_mod.f").  Moving this computation
!  to Get\_Canopy\_NOx now allows for a totally clean separation between
!  dry deposition routines and emissions routines in GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Canopy_NOx( HcoState, ExtState, Inst, RC )
!
! !USES:
!
    USE Drydep_Toolbox_Mod, ONLY : BIOFIT
!
! !ARGUMENTS:
!
    TYPE(HCO_State), POINTER        :: HcoState
    TYPE(Ext_State), POINTER        :: ExtState
    TYPE(MyInst),    POINTER        :: Inst
    INTEGER,         INTENT(INOUT)  :: RC
!
! !REMARKS:
!  For backwards compatibility, the bulk surface resistance is stored
!  in common block array CANOPYNOX in "commsoil.h".  Leave it like this
!  for the time being...we'll clean it up when we fix all of the soil
!  NOx routines.
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca     - Split off from "drydep_mod.f"
!  14 Jun 2012 - J.D. Maasakkers - Rewritten as a function of the
!                                     MODIS/Koppen biometype
!  09 Nov 2012 - M. Payer        - Replaced all met field arrays with State_Met
!                                   derived type object
!  13 Dec 2012 - R. Yantosca     - Removed ref to obsolete CMN_DEP_mod.F
!  28 Jul 2014 - C. Keller       - Added error trap for DRYCOEFF
!  11 Dec 2014 - M. Yannetti     - Added BIO_RESULT
!  12 May 2015 - R. Yantosca     - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Molecular weight of water [kg]
    REAL(hp), PARAMETER :: XMWH2O = 18e-3_hp

    ! Surface pressure??? [Pa]
    REAL(hp), PARAMETER :: PRESS  = 1.5e+5_hp
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: I,     J,     K,      KK
    INTEGER             :: DCSZ
    REAL(hp)            :: F0,    HSTAR, XMW
    REAL(hp)            :: DTMP1, DTMP2, DTMP3,  DTMP4, GFACT, GFACI
    REAL(hp)            :: RT,    RAD0,  RIX,    RIXX,  RDC,   RLUXX
    REAL(hp)            :: RGSX,  RCLX,  TEMPK,  TEMPC
    REAL(hp)            :: LAI,   SUNCOS, CLDFRC
    REAL(hp)            :: BIO_RESULT

    ! Arrays
    REAL(hp)            :: RI  (NBIOM)
    REAL(hp)            :: RLU (NBIOM)
    REAL(hp)            :: RAC (NBIOM)
    REAL(hp)            :: RGSS(NBIOM)
    REAL(hp)            :: RGSO(NBIOM)
    REAL(hp)            :: RCLS(NBIOM)
    REAL(hp)            :: RCLO(NBIOM)

    !=================================================================
    ! GET_CANOPY_NOX begins here!
    !=================================================================

    ! Set physical parameters
    HSTAR = 0.01e+0_hp              ! Henry's law constant
    F0    = 0.1e+0_hp               ! Reactivity factor for biological oxidation
    XMW   = 46e-3_hp               ! Molecular wt of NO2 (kg)

    ! Get size of DRYCOEFF (will be passed to routine BIOFIT)
    DCSZ = SIZE( ExtState%DRYCOEFF )

    ! Loop over surface boxes
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! Surface temperature [K] and [C]
       TEMPK = ExtState%T2M%Arr%Val(I,J)
       TEMPC = ExtState%T2M%Arr%Val(I,J) - 273.15e+0_hp

       ! Compute bulk surface resistance for gases.
       !
       !  Adjust external surface resistances for temperature;
       !  from Wesely [1989], expression given in text on p. 1296.
       RT = 1000.0e+0_hp * EXP( -TEMPC - 4.0e+0_hp )

       !--------------------------------------------------------------
       ! Get surface resistances - loop over biome types K
       !
       ! The land types within each grid square are defined using the
       ! Olson land-type database.  Each of the Olson land types is
       ! assigned a corresponding "deposition land type" with
       ! characteristic values of surface resistance components.
       ! There are 74 Olson land-types but only 11 deposition
       ! land-types (i.e., many of the Olson land types share the
       ! same deposition characteristics).  Surface resistance
       ! components for the "deposition land types" are from Wesely
       ! [1989] except for tropical forests [Jacob and Wofsy, 1990]
       ! and for tundra [Jacob et al., 1992].  All surface resistance
       ! components are normalized to a leaf area index of unity.
       !--------------------------------------------------------------
       !Loop over all biometypes
       DO K = 1, NBIOM

          ! Skip if not present
          IF ( Inst%LANDTYPE(K)%VAL(I,J) == 0.0_hp ) CYCLE

          ! Set second loop variable to K to allow snow/ice correction
          KK = K

          ! If the surface is snow or ice, then set K=3
          IF ( ExtState%ALBD%Arr%Val(I,J) > 0.4 ) KK = 3

          ! USE new MODIS/KOPPEN Biometypes to read data

          ! Read the internal resistance RI (minimum stomatal resistance
          ! for water vapor, per unit area of leaf) from the IRI array;
          ! a '9999' value means no deposition to stomata so we impose a
          ! very large value for RI.
          RI(K) = DBLE( SNIRI(KK) )
          IF ( RI(K) >= 9999.e+0_hp ) RI(K)= 1.e+12_hp

          ! Cuticular resistances IRLU read in from 'drydep.table'
          ! are per unit area of leaf; divide them by the leaf area index
          ! to get a cuticular resistance for the bulk canopy.  If IRLU is
          !'9999' it means there are no cuticular surfaces on which to
          ! deposit so we impose a very large value for RLU.
          IF ( SNIRLU(KK) >= 9999 .OR. &
               ExtState%LAI%Arr%Val(I,J) <= 0e+0_hp ) THEN
             RLU(K)  = 1.e+6_hp
          ELSE
             RLU(K)= DBLE(SNIRLU(KK)) / ExtState%LAI%Arr%Val(I,J) + RT
          ENDIF

          ! The following are the remaining resistances for the Wesely
          ! resistance-in-series model for a surface canopy
          ! (see Atmos. Environ. paper, Fig.1).
          RAC(K)  = MAX( DBLE( SNIRAC(KK)  ),      1e+0_hp )
          RGSS(K) = MAX( DBLE( SNIRGSS(KK) ) + RT, 1e+0_hp )
          RGSO(K) = MAX( DBLE( SNIRGSO(KK) ) + RT, 1e+0_hp )
          RCLS(K) =      DBLE( SNIRCLS(KK) ) + RT
          RCLO(K) =      DBLE( SNIRCLO(KK) ) + RT

          IF (  RAC(K) >= 9999.e+0_hp ) RAC(K)  = 1e+12_hp
          IF ( RGSS(K) >= 9999.e+0_hp ) RGSS(K) = 1e+12_hp
          IF ( RGSO(K) >= 9999.e+0_hp ) RGSO(K) = 1e+12_hp
          IF ( RCLS(K) >= 9999.e+0_hp ) RCLS(K) = 1e+12_hp
          IF ( RCLO(K) >= 9999.e+0_hp ) RCLO(K) = 1e+12_hp

          !-------------------------------------------------------------
          ! Adjust stomatal resistances for insolation and temperature:
          !
          ! Temperature adjustment is from Wesely [1989], equation (3).
          !
          ! Light adjustment by the function BIOFIT is described by Wang
          ! [1996].  It combines:
          !
          ! - Local dependence of stomal resistance on the intensity I
          !   of light impinging the leaf; this is expressed as a
          !   multiplicative factor I/(I+b) to the stomatal resistance
          !   where b = 50 W m-2
          !   (equation (7) of Baldocchi et al. [1987])
          ! - Radiative transfer of direct and diffuse radiation in the
          !   canopy using equations (12)-(16) from Guenther et al.
          !   [1995]
          ! - Separate accounting of sunlit and shaded leaves using
          !   equation (12) of Guenther et al. [1995]
          ! - Partitioning of the radiation at the top of the canopy
          !   into direct and diffuse components using a
          !   parameterization to results from an atmospheric radiative
          !   transfer model [Wang, 1996]
          !
          ! The dependent variables of the function BIOFIT are the leaf
          ! area index (XYLAI), the cosine of zenith angle (SUNCOS) and
          ! the fractional cloud cover (CFRAC).  The factor GFACI
          ! integrates the light dependence over the canopy depth; so
          ! be scaled by LAI to yield a bulk canopy value because that's
          ! already done in the GFACI formulation.
          !-------------------------------------------------------------

          ! Radiation @ sfc [W/m2]
          RAD0 = ExtState%RADSWG%Arr%Val(I,J)

          ! Internal resistance
          RIX  = RI(K)

          ! Skip the following block if the resistance RIX is high
          IF ( RIX < 9999e+0_hp ) THEN
             GFACT = 100.0e+0_hp

             IF ( TEMPC > 0.e+0_hp .AND. TEMPC < 40.e+0_hp) THEN
                GFACT = 400.e+0_hp / TEMPC / ( 40.0e+0_hp - TEMPC )
             ENDIF

             GFACI = 100.e+0_hp

             IF ( RAD0 > 0e+0_hp .AND. ExtState%LAI%Arr%Val(I,J) > 0e+0_hp ) THEN

                LAI    = ExtState%LAI%Arr%Val(I,J)
                SUNCOS = ExtState%SUNCOS%Arr%Val(I,J)
                CLDFRC = ExtState%CLDFRC%Arr%Val(I,J)

                BIO_RESULT = BIOFIT( ExtState%DRYCOEFF, LAI, &
                              SUNCOS, CLDFRC, SIZE(ExtState%DRYCOEFF) )

                GFACI= 1e+0_hp / BIO_RESULT
             ENDIF

             RIX = RIX * GFACT * GFACI
          ENDIF

          ! Compute aerodynamic resistance to lower elements in lower
          ! part of the canopy or structure, assuming level terrain -
          ! equation (5) of Wesely [1989].
          RDC = 100.e+0_hp*(1.0e+0_hp+1000.0e+0_hp/(RAD0 + 10.e+0_hp))

          ! Loop over species; species-dependent corrections to resistances
          ! are from equations (6)-(9) of Wesely [1989].
          !
          ! NOTE: here we only consider NO2 (bmy, 6/22/09)
          RIXX   = RIX * DIFFG( TEMPK, PRESS, XMWH2O ) / &
                         DIFFG( TEMPK, PRESS, XMW    )   &
                 + 1.e+0_hp / ( HSTAR/3000.e+0_hp + 100.e+0_hp*F0  )

          RLUXX  = 1.e+12_hp

          IF ( RLU(K) < 9999.e+0_hp ) THEN
             RLUXX = RLU(K) / ( HSTAR / 1.0e+05_hp + F0 )
          ENDIF

          ! To prevent virtually zero resistance to species with huge HSTAR,
          ! such as HNO3, a minimum value of RLUXX needs to be set.
          ! The rationality of the existence of such a minimum is
          ! demonstrated by the observed relationship between Vd(NOy-NOx)
          ! and Ustar in Munger et al.[1996]; Vd(HNO3) never exceeds 2 cm/s
          ! in observations. The corresponding minimum resistance is 50 s/m.
          ! was introduced by J.Y. Liang on 7/9/95.
          RGSX = 1e+0_hp / ( HSTAR/1e+5_hp/RGSS(K) + F0/RGSO(K) )
          RCLX = 1e+0_hp / ( HSTAR/1e+5_hp/RCLS(K) + F0/RCLO(K) )

          ! Get the bulk surface resistance of the canopy
          ! from the network of resistances in parallel and in series
          ! (Fig. 1 of Wesely [1989])
          DTMP1 = 1.e+0_hp / RIXX
          DTMP2 = 1.e+0_hp / RLUXX
          DTMP3 = 1.e+0_hp / ( RAC(K) + RGSX )
          DTMP4 = 1.e+0_hp / ( RDC      + RCLX )

          ! Save the within canopy depvel of NOx, used in calculating
          ! the canopy reduction factor for soil emissions [1/s]
          Inst%CANOPYNOX(I,J,K) = DTMP1 + DTMP2 + DTMP3 + DTMP4

       ENDDO !K
    ENDDO !I
    ENDDO !J

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Get_Canopy_NOx
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DiffG
!
! !DESCRIPTION: Function DiffG calculates the molecular diffusivity [m2/s] in
!  air for a gas X of molecular weight XM [kg] at temperature TK [K] and
!  pressure PRESS [Pa].
!\\
!\\
! !INTERFACE:
!
  FUNCTION DiffG( TK, PRESS, XM ) RESULT( DIFF_G )
!
! !INPUT PARAMETERS:
!
    REAL(hp), INTENT(IN) :: TK      ! Temperature [K]
    REAL(hp), INTENT(IN) :: PRESS   ! Pressure [hPa]
    REAL(hp), INTENT(IN) :: XM      ! Molecular weight of gas [kg]
!
! !RETURN VALUE:
!
    REAL(hp)             :: DIFF_G  ! Molecular diffusivity [m2/s]
!
! !REMARKS:
!  We specify the molecular weight of air (XMAIR) and the hard-sphere molecular
!  radii of air (RADAIR) and of the diffusing gas (RADX).  The molecular
!  radius of air is given in a Table on p. 479 of Levine [1988].  The Table
!  also gives radii for some other molecules.  Rather than requesting the user
!  to supply a molecular radius we specify here a generic value of 2.E-10 m for
!  all molecules, which is good enough in terms of calculating the diffusivity
!  as long as molecule is not too big.
!
! !REVISION HISTORY:
!     22 Jun 2009 - R. Yantosca - Copied from "drydep_mod.f"
!     07 Jan 2016 - E. Lundgren - Update Avogadro's # to NIST 2014 value
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(hp), PARAMETER  :: XMAIR  = 28.8e-3_hp
    REAL(hp), PARAMETER  :: RADAIR = 1.2e-10_hp
    REAL(hp), PARAMETER  :: PI     = 3.14159265358979323e+0_hp
    REAL(hp), PARAMETER  :: RADX   = 1.5e-10_hp
    REAL(hp), PARAMETER  :: RGAS   = 8.32e+0_hp
    REAL(hp), PARAMETER  :: AVOGAD = 6.022140857e+23_hp
!
! !LOCAL VARIABLES:
!
    REAL(hp)             :: AIRDEN, Z, DIAM, FRPATH, SPEED

    !=================================================================
    ! DIFFG begins here!
    !=================================================================

    ! Air density
    AIRDEN = ( PRESS * AVOGAD ) / ( RGAS * TK )

    ! DIAM is the collision diameter for gas X with air.
    DIAM   = RADX + RADAIR

    ! Calculate the mean free path for gas X in air:
    ! eq. 8.5 of Seinfeld [1986];
    Z      = XM  / XMAIR
    FRPATH = 1e+0_hp /( PI * SQRT( 1e+0_hp + Z ) * AIRDEN*( DIAM**2 ) )

    ! Calculate average speed of gas X; eq. 15.47 of Levine [1988]
    SPEED  = SQRT( 8e+0_hp * RGAS * TK / ( PI * XM ) )

    ! Calculate diffusion coefficient of gas X in air;
    ! eq. 8.9 of Seinfeld [1986]
    DIFF_G = ( 3e+0_hp * PI / 32e+0_hp ) * ( 1e+0_hp + Z ) * FRPATH * SPEED

  END FUNCTION DiffG
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Dep_N
!
! !DESCRIPTION: Subroutine GET\_DEP\_N sums dry and wet deposition since prev.
! timestep and calculates contribution to fertilizer N source. Output is in
! kg NO/m2.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Dep_N( I, J, ExtState, HcoState, Inst, DEP_FERT )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)         :: I
    INTEGER,  INTENT(IN)         :: J
    TYPE(Ext_State), POINTER     :: ExtState
    TYPE(HCO_State), POINTER     :: HcoState
    TYPE(MyInst),    POINTER     :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Dep emitted as Fert [kgNO/m2]
    REAL(hp) ,  INTENT(INOUT) :: DEP_FERT
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(hp),  PARAMETER :: TAU_MONTHS   = 6. ! Decay rate of dep. N [months]
    REAL(hp),  PARAMETER :: SECPERDAY    = 86400.e+0_hp
    REAL(hp),  PARAMETER :: DAYSPERMONTH = 30.
!
! !LOCAL VARIABLES:
!
    REAL(hp)             :: DRYN  ! Dry dep. N since prev timestep
                                  ! Units ng N/m2/s
    REAL(hp)             :: WETN  ! Wet dep. N since prev timestep
    REAL(hp)             :: DEPN  ! dep. N since prev timestep

    REAL(hp)             :: C1
    REAL(hp)             :: C2
    REAL(hp)             :: TAU_SEC
    REAL(hp)             :: TS_SEC

    !Total all N species & convert molec/cm2/s --> kg NO/m2/s
    DRYN = SOURCE_DRYN( I, J, ExtState, HcoState, Inst )

    !Total all N species & convert kg/s --> kg NO/m2/s
    WETN = SOURCE_WETN( I, J, ExtState, HcoState )

    ! Sum wet and dry deposition [kg NO/m2/s]
    DEPN = DRYN + WETN

    ! Emission Timestep in seconds
    TS_SEC = HcoState%TS_EMIS

    !Do mass balance (see Intro to Atm Chem Chap. 3)
    !m(t) = m(0) * exp(-t/tau) + Source * tau * (1 - exp(-t/tau))

    !convert months -->  seconds (assume 30 days months)
    TAU_SEC = TAU_MONTHS * DAYSPERMONTH * SECPERDAY

    C1 = EXP( - TS_SEC / TAU_SEC)
    C2 = 1.e+0_hp - C1

    ! kg NO/m2
    ! NOTE: DEP_RESERVOIR is stored in kg NO/m3, but we just assume
    ! that this is kg NO/m2.
    Inst%DEP_RESERVOIR(I,J) = ( Inst%DEP_RESERVOIR (I,J) * C1 ) &
                            + DEPN * TAU_SEC * C2

    ! 40% runoff.
    DEP_FERT = Inst%DEP_RESERVOIR(I,J) * 0.6e+0_hp

  END SUBROUTINE Get_Dep_N
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Source_DryN
!
! !DESCRIPTION: Subroutine SOURCE\_DRYN gets dry deposited Nitrogen since
!               last emission time step, converts to kg NO/m2/s.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Source_Dryn( I, J, ExtState, HcoState, Inst ) RESULT( DRYN )
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN) :: I
    INTEGER,         INTENT(IN) :: J
    TYPE(Ext_State), POINTER    :: ExtState   ! Module options
    TYPE(HCO_State), POINTER    :: HcoState   ! Output obj
    TYPE(MyInst),    POINTER    :: Inst
!
! !RETURN VALUE:
!
    REAL(hp)                      :: DRYN       ! Dry dep. N since prev timestep
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp),  PARAMETER   :: CM2_PER_M2  = 1.e+4_hp
    REAL(hp)               :: NTS

    ! Divide through by number of chemistry timesteps
    ! because DRY_TOTN is summed over chemistry timesteps
    ! need to get average

    !Molecules/cm2/s --> kg NO/m2/s
    NTS  = HcoState%TS_EMIS / HcoState%TS_CHEM
    DRYN = ExtState%DRY_TOTN%Arr%Val(I,J) * CM2_PER_M2 / NTS / &
           HcoState%Phys%Avgdr * HcoState%Spc(Inst%IDTNO)%EmMW_g / 1000.0e+0_hp

  END FUNCTION Source_DryN
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Source_WetN
!
! !DESCRIPTION: Subroutine Source\_WetN gets wet deposited Nitrogen since
!  last emission time step, converts to kg NO/m2/s.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Source_WetN( I, J, ExtState, HcoState ) RESULT(WETN )
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN) :: I
    INTEGER,         INTENT(IN) :: J
    TYPE(Ext_State), POINTER    :: ExtState   ! Module options
    TYPE(HCO_State), POINTER    :: HcoState   ! Output obj
!
! !RETURN VALUE:
!
    REAL(hp)                      :: WETN       ! Dry dep. N since prev timestep
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp)               :: NTS,    AREA_M2

    ! Divide through by number of transport timesteps
    ! because WET_TOTN is summed over transport timesteps
    ! need to get average

    NTS     = HcoState%TS_EMIS / HcoState%TS_DYN
    AREA_M2 = HcoState%Grid%AREA_M2%Val(I,J)

    ! Total N wet dep
    WETN = ExtState%WET_TOTN%Arr%Val(I,J) / AREA_M2 / NTS

  END FUNCTION Source_WetN
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SoilTemp
!
! !DESCRIPTION: Function SoilTemp computes the temperature-dependent term
!  of the soil NOx emissions in ng N/m2/s and converts to molec/cm2/s
!\\
!\\
! !INTERFACE:
!
  FUNCTION SoilTemp( NN, TC, GWET ) RESULT( SOIL_TEMP )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: NN           ! Soil biome type
    REAL(hp), INTENT(IN) :: TC           ! Surface air temperature [C]
    REAL(hp), INTENT(IN) :: GWET         ! Top soil moisture
!
! !RETURN VALUE:
!
    REAL(hp)             :: SOIL_TEMP    ! Temperature-dependent term of
                                         ! soil NOx emissions [unitless]
!
! !REMARKS:
!    Based on Ormeci et al., [1999] and Otter et al., [1999]
!    there exists and entirely exponential relationship between
!    temperature and soil NOx emissions at constant soil moisture
!    Therefore we use the following relationship based
!    on Yienger and Levy et al., [1995] for temperatures 0-30C:
!
!
!         f(T) =  exp( 0.103+/-0.04 * T )
!           in ng N/m2/s
!
!
!     where T is the temperature in degrees Celsius....Below
!     0 C, we assume emissions are zero because they are insignificant
!     for the purposes of this global source. ...
!
!  References:
!  ============================================================================
!  (1 ) Ormeci, B., S. L. Sanin, and J. J. Pierce, Laboratory study of
!        NO flux from agricultural soil: Effects of soil moisture, pH,
!        and temperature, J. Geophys. Res., 104 ,16211629, 1999.
!  (2 ) Otter, L. B., W. X. Yang, M. C. Scholes, and F. X. Meixner,
!        Nitric oxide emissions from a southern African savanna, J.
!        Geophys. Res., 105 , 20,69720,706, 1999.
!  (3 ) Yienger, J.J, and H. Levy, Empirical model of global soil-biogenic
!        NOx emissions, J. Geophys. Res., 100, D6, 11,447-11464, June 20, 1995.
!
! !REVISION HISTORY:
!  17 Aug 2009 - R. Yantosca - Initial Version
!  17 Aug 2009 - R. Yantosca - Added ProTeX headers
!  31 Jan 2011 - R. Hudman   - Added new soil T dependance
!  31 Jan 2011 - R. Hudman   - Updated headers
!  12 May 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    REAL(hp)  :: TMMP

    !==============================================================
    ! 1) Convert from Surface Temp  --> Soil Temp
    !==============================================================

    ! Save surface air temp in shadow variable TMMP
    TMMP   = TC

    ! DRY
    IF ( GWET < 0.3e+0_hp ) THEN

       ! Convert surface air temperature to model temperature
       ! by adding 5 degrees C to model temperature
       TMMP = TMMP + 5e+0_hp

    ! WET
    ELSE

       TMMP = SOILTA(NN) * TMMP + SOILTB(NN)

    ENDIF

    !==============================================================
    ! 2) Compute Temperature Dependence
    !==============================================================

    ! Compute the soil temperature dependence term according
    ! to equations 9b, 9a of Yienger & Levy [1995].
    ! We now assume that soil response is exponential 0-30C
    ! based on latest observations, caps at 30C

    IF ( TMMP <= 0e+0_hp ) THEN

       ! No soil emissions if temp below freezing
       SOIL_TEMP = 0e+0_hp

    ELSE

       ! Caps temperature response at 30C
       IF ( TMMP >= 30.e+0_hp ) TMMP = 30.e+0_hp

       SOIL_TEMP =  EXP( 0.103 * TMMP )

    ENDIF

  END FUNCTION SoilTemp
!EOC
!----------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SoilWet
!
! !DESCRIPTION: Function SoilWet returns the soil moisture scaling
!  of soil NOx emissions (values from 0-1).
!\\
!\\
! !INTERFACE:
!
  FUNCTION SoilWet( GWET, ARID, NONARID ) RESULT( WETSCALE )
!
! !INPUT PARAMETERS:
!
    ! Top soil wetness [unitless]
    REAL(hp), INTENT(IN) :: GWET

    ! Fraction of arid & non-arid soil in the gridbox
    REAL(hp), INTENT(IN) :: ARID
    REAL(hp), INTENT(IN) :: NONARID
!
! !RETURN_VALUE:
!
    ! A scaling term between 0-1 based on soil moisture
    REAL(hp)             :: WETSCALE
!
! !REMARKS:
!  Soil moisture and temperature and now decoupled, the temperature
!  term is scaled with a value from 0-1 based on water filled pore space
!  WFPS in top-soil.
!
!  From N.E. Moore thesis:
!  The response of SNOx is not monotonic to WFPS. SNOx are low for the
!  extreme values of WFPS (0 and 1). For low values, emissions are
!  substrate-limited. For high values, emissions are trapped and cannot
!  diffuse to the surface [Yan et al., 2005]. SNOx dependence on soil
!  moisture is best described as a Poisson function [Parsons et al., 1996;
!  Otter et al., 1999; Pierce and Aneja, 2000; Kirkman et al., 2001;
!  van Dijk and Meixner, 2001; van Dijk et al., 2002]:
!
!     scaling = a*x*exp(-b*x^2)
!
!  where the values of a and b are chosen such that the maximum value
!  (unity) occurs for WFPS=0.3, which laboratory and field measurements have
!  found to be the optimal value for emissions in most soils. The typical
!  range of values are 0.2 (arid) up to 0.45 (floodplain)
!  [Yang and Meixner, 1997; Ormeci et al., 1999].
!
!  Rice paddies no longer have to be scaled as in the Yienger & Levy model.
!
!  References:
!  ============================================================================
!  (1 ) Galbally, I. E., and R. Roy, Loss of fixed nitrogen from soils
!        by nitric oxide exhalation, Nature, 275 , 734735, 1978.
!  (2 ) Kirkman, G. A., W. X. Yang, and F. X. Meixner, Biogenic nitric
!        oxide emissions upscaling: An approach for Zimbabwe, Global
!        Biogeochemical Cycles, 15 ,1005 1020, 2001.
!  (3 ) Ormeci, B., S. L. Sanin, and J. J. Pierce, Laboratory study of NO
!        flux from agricultural soil: Effects of soil moisture, pH, and
!        temperature, J. Geophys. Res., 104 , 16211629, 1999.
!  (4 ) Otter, L. B., W. X. Yang, M. C. Scholes, and F. X. Meixner,
!        Nitric oxide emissions from a southern African savanna, J.
!        Geophys. Res., 105 , 20,69720,706, 1999.
!  (5 ) Parsons, D. A., M. C. Scholes, R. J. Scholes, and J. S. Levine,
!        Biogenic NO emissions from savanna soils as a function of fire
!        regime, soil type, soil nitrogen, and water status, J. Geophys.
!        Res., 101 , 23,68323,688, 1996.
!  (6 ) Pierce, J. J., and V. P. Aneja, Nitric oxide emissions from
!        engineered soil systems, Journal of Environmental Engineering,
!        pp. 225232, 2000.
!  (7 ) van Dijk, S. M., and J. H. Duyzer, Nitric oxide emissions from
!        forest soils, J. Geophys. Res., 104 , 15,95515,961, 1999.
!  (8 ) van Dijk, S. M., and F. X. Meixner, Production and consumption of
!        NO in forest and pasture soils from the Amazon basin, Water, Air,
!        and Soil Pollution: Focus 1 , pp. 119130, 2001.
!  (9 ) Yang, W. X., and F. X. Meixner, Gaseous Nitrogen Emissions from
!        Grasslands, CAB Int., Wallingford, UK, 1997, 67-71.
!
! !REVISION HISTORY:
!  17 Aug 2009 - R. Yantosca - Columnized and cleaned up
!  17 Aug 2009 - R. Yantosca - Added ProTeX headers
!  31 Jan 2011 - R. Hudman   - Rewrote scaling scheme
!  31 Jan 2011 - R.Hudman    - Updated ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !Scale by soil moisture
    IF ( ARID .GE. NONARID .AND. ARID .NE. 0) THEN
       !Arid, max Poisson = 0.2
       WETSCALE = 8.24 * GWET * EXP(-12.5*GWET*GWET)
    ELSE
       !Non-arid, max Poisson = 0.3
       WETSCALE = 5.5 * GWET * EXP( -5.55 * GWET * GWET)
    ENDIF

  END FUNCTION SoilWet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SoilCrf
!
! !DESCRIPTION: Computes the canopy reduction factor for the soil NOx
!  emissions according to Jacob \% Bakwin [1991] (and as used in Wang
!  et al [1998]).
!\\
!\\
! !INTERFACE:
!
  FUNCTION SoilCrf( K, LAI, CPYNOX, WINDSQR, SUNCOS ) RESULT( SOIL_CRF )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: K          ! Soil biome type
    REAL(hp),  INTENT(IN) :: LAI        ! Leaf area index [cm2/cm2]
    REAL(hp),  INTENT(IN) :: CPYNOX     ! Bulk sfc resistance to NOx [1/s]
    REAL(hp),  INTENT(IN) :: WINDSQR    ! Square of sfc wind speed [m2/s2]
    REAL(hp),  INTENT(IN) :: SUNCOS     ! Cosine of solar zenith angle
!
! !RETURN_VALUE:
!
    REAL(hp)              :: SOIL_CRF   ! Canopy reduction factor (see below)
!
! !REMARKS:
!  Also note, CANOPYNOX (the bulk surface resistance to NOx) is computed
!  in routine GET_CANOPY_NOx (in "canopy_nox_mod.f") and is passed here
!  as an argument.
!
! !REVISION HISTORY:
!  17 Aug 2009 - R. Yantosca - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Ventilation velocity for NOx, day & night values [m/s]
    REAL(hp),  PARAMETER :: VFDAY   = 1.0e-2_hp
    REAL(hp),  PARAMETER :: VFNIGHT = 0.2e-2_hp
!
! !LOCAL VARIABLES:
!
    REAL(hp) :: VFNEW

    ! Pick proper ventilation velocity for day or night
    IF ( SUNCOS > 0e+0_hp ) THEN
       VFNEW = VFDAY
    ELSE
       VFNEW = VFNIGHT
    ENDIF

    ! If the leaf area index and the bulk surface resistance
    ! of the canopy to NOx deposition are both nonzero ...
    IF ( LAI > 0e+0_hp .and. CPYNOX > 0e+0_hp ) THEN

       ! Adjust the ventilation velocity.
       ! NOTE: SOILEXC(21) is the canopy wind extinction
       ! coefficient for the tropical rainforest biome.
       VFNEW    = (VFNEW * SQRT( WINDSQR/9e+0_hp * 7e+0_hp/LAI     ) * &
                  ( SOILEXC(21)  / SOILEXC(K) ))

       ! Soil canopy reduction factor
       SOIL_CRF = CPYNOX / ( CPYNOX + VFNEW )

    ELSE

       ! Otherwise set the soil canopy reduction factor to zero
       SOIL_CRF = 0e+0_hp

    ENDIF

  END FUNCTION SoilCrf
!EOC
!-----------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FertAdd
!
! !DESCRIPTION: Function FertAdd computes fertilizer emissions
!\\
!\\
! !INTERFACE:
!
  FUNCTION FertAdd( SOILFRT, DEPN ) RESULT( FERT_ADD )
!
! !INPUT PARAMETERS:
!
    REAL(hp), INTENT(IN) :: DEPN       ! N emissions from deposition
    REAL(hp), INTENT(IN) :: SOILFRT    ! N emissions from fertilizers
                                       !  read in from disk and passed
                                       !  here as an argument [ng N/m2/s]
!
! !RETURN_VALUE:
!
    REAL(hp)             :: FERT_ADD   ! Total Fert emissions
!
! !REMARKS:
!  We use a new spatially explicit data set of chemical and manure fert
!  (native resolution 0.5\B0x0.5\B0) from Potter et al., [2010]
!  distributed using MODIS EVI seasonality as described in
!  N.E. Moore thesis, and Hudman et al., in prep.
!
!  In previous model, fertilizer emissions were emitted instantaneously as
!  2.5% of applied fertilizer, independent of soil moisture/soil temperature,
!  so that they were constant over the growing season.
!
!  Similar to the YL  parameterization, we now treat fertilizer emissions
!  as part of the Aw. If we treat the wet biome coefficient as a measure of
!  available N multiplied by a mean emission rate, we can treat fertilizer
!  N in the same manner.
!
!  AW = SOILAW(BinewsoilAWS_08112011_emissonlyome) + N available in soil
!       x mean emission rate
!
!  Instead of choosing an emission rate for each box equivalent to 2.5%
!  of applied N yearly as done in the YL scheme, we chose the mean emission
!  rate so that the total global above canopy SNOx due to fertilizer matches
!  observed estimates of fertilizer emissions of 1.8 Tg N yr-1 from Stehfest
!  and Bouman [2006].  This treatment allows for interannual and daily
!  variability in the strength  of response to temperature and precipitation.
!  Note: this scaling must be set for each resolution.
!
!  References:
!  ============================================================================
!  (1 ) Potter, P., Ramankutty, N., Bennett, E.,  and Donner, S.:
!        Characterizing the Spatial Patterns of Global Fertilizer
!        Application and Manure Production, Earth Interactions,
!        in press, 2010.
!  (2 ) Stehfest, E. and L. Bouwman, N2O and NO emission from
!        agricultural fields and soils under natural vegetation:
!        summarizing available measurement data and modeling
!        of global annual emissions, Nutrient Cycling in Agroecosystems
!        (2006), 74:207-228 DOI 10.1007/s10705-006-9000-7.
!
! !REVISION HISTORY:
!  17 Aug 2009 - R. Yantosca - Columnized and cleaned up
!  17 Aug 2009 - R. Yantosca - Added ProTeX headers
!  31 Jan 2011 - R. Hudman   - Rewrote pulsing scheme
!  31 Jan 2011 - R. Hudman   - Updated ProTex headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Seconds per year
    REAL(hp), PARAMETER      :: SEC_PER_YEAR = 3.1536e+7_hp


    ! Soil fert and dep [ kg/m2 ], a measure of N avail. in soil
    FERT_ADD = ( SOILFRT + DEPN ) / SEC_PER_YEAR

    ! Convert [ng N/m2] --> [kg/m2/s]
    ! (scale needed to force fert emiss of 1.8 Tg N/yr w/o canopy uptake)
!    FERT_ADD = FERT_ADD * FERT_SCALE

  END FUNCTION FERTADD
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Pulsing
!
! !DESCRIPTION: Function Pulsing calculates the increase (or "pulse") of
!  soil NOx emission that happens after preciptiation falls on dry soil.
!                                                                             .
!  According to  Yan et al., [2005] , this pulsing process is thought to
!  be due to a release of inorganic nitrogen trapped on top of the dry soil
!  and a subsequent reactivation of water-stressed bacteria, which then
!  metabolize the excess nitrogen. This can happen in seasonally dry
!  grasslands and savannahs or over freshly fertilized fields.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Pulsing( GWET,      TS_EMIS,             &
                    GWET_PREV, PFACTOR, DRYPERIOD ) &
                    RESULT( THE_PULSING )
!
! !INPUT PARAMETERS:
!
    REAL(hp), INTENT(IN)    :: GWET        ! Soil Moisture
    REAL*4,   INTENT(IN)    :: TS_EMIS     ! Emissions timestep [s]

! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp), INTENT(INOUT) :: GWET_PREV   ! Soil Moisture Prev timestep
    REAL(sp), INTENT(INOUT) :: PFACTOR     ! Pulsing Factor
    REAL(sp), INTENT(INOUT) :: DRYPERIOD   ! Dry period length in hours
!
! !RETURN VALUE:
!
    REAL(hp)                :: THE_PULSING ! Factor to multiply baseline
                                           ! emissions by to account for
                                           ! soil pulsing of all types
!
! !REMARKS:
!  Soil NOx emissions consist of baseline emissions plus discrete "pulsing"
!  episodes.  We follow thw Yan et al., [2005] algorithm, where the pulse
!  (relative to the flux prewetting) is determined by the antecedent dry
!  period, with a simple logarithmic relationship,
!
!     PFACTOR = 13.01 ln ( DRYPERIOD ) -  53.6
!
!  where PFACTOR is the magnitude of peak flux relative to prewetting flux,
!  and DRYPERIOD  is the length of the antecedent dry period in hours.
!
!  The pulse decays with
!
!     PFACTOR = PFACTOR * EXP( -0.068e+0_hp * DTSRCE )
!
!  References:
!  ============================================================================
!  (1 ) Yan, X., T. Ohara, and H. Akimoto (2005), Statistical modeling of
!        global soil NOx emissions, Global Biogeochem. Cycles, 19, GB3019,
!        doi:10.1029/2004GB002276.Section 2.3.3
!
! !REVISION HISTORY:
!  17 Aug 2009 - R. Yantosca - Columnized and cleaned up
!  17 Aug 2009 - R. Yantosca - Added ProTeX headers
!  31 Jan 2011 - R. Hudman   - Rewrote pulsing scheme
!  31 Jan 2011 - R. Hudman   - Updated ProTex header
!  29 May 2013 - R. Yantosca - Bug fix: prevent log(0) from happening
!  21 Oct 2014 - C. Keller   - Limit PFACTOR to 1.
!  12 May 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp)  :: DTSRCE, GDIFF

    !=================================================================
    ! PULSING begins here!
    !=================================================================

    ! Emission timestep [s --> hours]
    DTSRCE = TS_EMIS / 3600e+0_hp

    ! If soil moisture less than 0.3 and no pulse is taking place
    IF ( GWET < 0.3e+0_hp .and. PFACTOR == 1.e+0_hp) THEN

       ! Get change in soil moisture since previous timestep
       GDIFF = ( GWET - GWET_PREV )

       ! If change in soil moisture is > 0.01 (rains)
       IF ( GDIFF > 0.01 ) THEN

          ! Initialize new pulse factor (dry period hours)
          IF ( DRYPERIOD > 0e+0_hp ) THEN
             PFACTOR = 13.01e+0_hp * LOG( DRYPERIOD ) - 53.6e+0_hp
          ELSE
             PFACTOR = -53.6e+0_hp
          ENDIF

          ! If dry period < ~3 days then no pulse
          IF ( PFACTOR < 1.0 ) PFACTOR = 1.0

          ! Reinitialize dry period
          DRYPERIOD = 0e+0_hp

       ! If no rain (i.e.,  change in soil moisture is < 0.01)
       ELSE

          ! Add one timestep to dry period
          DRYPERIOD = DRYPERIOD + DTSRCE

       ENDIF

    ! If box is already pulsing , then decay pulse one timestep
    ELSEIF ( PFACTOR /= 1.e+0_hp) THEN

       ! Decay pulse
       PFACTOR   = PFACTOR * EXP( -0.068e+0_hp * DTSRCE )

       ! Update dry period
       IF ( GWET < 0.3e+0_hp ) DRYPERIOD = DRYPERIOD + DTSRCE

       ! If end of pulse
       IF ( PFACTOR < 1.e+0_hp ) PFACTOR = 1.e+0_hp

    ENDIF

    ! Update soil moisture holder for previous timestep
    GWET_PREV = GWET

    ! Return the pulsing factor
    THE_PULSING = PFACTOR

  END FUNCTION Pulsing
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

    ! Make sure pointers are not dangling
    Inst%DRYCOEFF  => NULL()
    Inst%CLIMARID  => NULL()
    Inst%CLIMNARID => NULL()
    Inst%SOILFERT  => NULL()
    Inst%LANDTYPE  => NULL()

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
  SUBROUTINE InstRemove ( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState      ! Extension state
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                     :: RC
    INTEGER                     :: I
    TYPE(MyInst), POINTER       :: PrevInst
    TYPE(MyInst), POINTER       :: Inst

    !=================================================================
    ! InstRemove begins here!
    !=================================================================

    ! Get instance. Also archive previous instance.
    PrevInst => NULL()
    Inst     => NULL()
    CALL InstGet ( ExtState%SoilNOx, Inst, RC, PrevInst=PrevInst )

    ! Instance-specific deallocation
    IF ( ASSOCIATED(Inst) ) THEN

       ! Deallocate arrays
       IF ( ALLOCATED  ( Inst%DRYPERIOD     ) ) DEALLOCATE ( Inst%DRYPERIOD     )
       IF ( ALLOCATED  ( Inst%PFACTOR       ) ) DEALLOCATE ( Inst%PFACTOR       )
       IF ( ALLOCATED  ( Inst%GWET_PREV     ) ) DEALLOCATE ( Inst%GWET_PREV     )
       IF ( ALLOCATED  ( Inst%CANOPYNOX     ) ) DEALLOCATE ( Inst%CANOPYNOX     )
       IF ( ALLOCATED  ( Inst%DEP_RESERVOIR ) ) DEALLOCATE ( Inst%DEP_RESERVOIR )
       IF ( ALLOCATED  ( Inst%SpcScalVal    ) ) DEALLOCATE ( Inst%SpcScalVal    )
       IF ( ALLOCATED  ( Inst%SpcScalFldNme ) ) DEALLOCATE ( Inst%SpcScalFldNme )

       ! Deallocate LANDTYPE vector
       IF ( ASSOCIATED(Inst%LANDTYPE) ) THEN
          DO I = 1,NBIOM
             IF ( ASSOCIATED(Inst%LANDTYPE(I)%VAL) ) &
                DEALLOCATE(Inst%LANDTYPE(I)%Val)
          ENDDO
          DEALLOCATE ( Inst%LANDTYPE )
       ENDIF

       ! Eventually deallocate DRYCOEFF. Make sure ExtState DRYCOEFF pointer is
       ! not dangling!
       IF ( ASSOCIATED ( Inst%DRYCOEFF ) ) THEN
          DEALLOCATE ( Inst%DRYCOEFF )
          ExtState%DRYCOEFF => NULL()
       ENDIF

       ! Free pointers
       IF ( ASSOCIATED( Inst%CLIMARID  ) ) DEALLOCATE ( Inst%CLIMARID  )
       IF ( ASSOCIATED( Inst%CLIMNARID ) ) DEALLOCATE ( Inst%CLIMNARID )
       IF ( ASSOCIATED( Inst%SOILFERT  ) ) DEALLOCATE ( Inst%SOILFERT  )

       ! ----------------------------------------------------------------
       ! Pop off instance from list
       ! ----------------------------------------------------------------
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
END MODULE HCOX_SoilNOx_Mod
!EOM
