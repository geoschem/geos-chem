!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoi_gc_diagn_mod.F90
!
! !DESCRIPTION: Module HCOI\_GC\_Diagn\_Mod.F90 is the GEOS-Chem interface 
! module for the HEMCO diagnostics. For every GEOS-Chem emissions diagnostics,
! a corresponding HEMCO diagnostics is created. The HEMCO diagnostics become
! (automatically) filled and updated when calling HEMCO. They are passed
! back to GEOS-Chem when writing the diagnostics (e.g. in diag3.F).
!\\
!\\
! Notes:
! \begin{itemize}
! \item The category specific diagnostics (anthropogenic, aircraft, etc.)
!  explicitly assume certain category numbers in the HEMCO configuration 
!  file (e.g. Cat=1 for anthropogenic, Cat=20 for aircraft, etc.).
!  Diagnostics will not represent what they should if these category numbers
!  get changed!
! \item Most biofuel emissions are included in the anthropogenic inventories
!  and hence not distinguishable from those. For most compounds, no biofuel 
!  diagnostics are written.
! \item In HEMCO, ocean sinks are treated as drydep and the calculated 
!  deposition velocities are passed to drydep\_mod.F. Hence, no Acetone or ALD2
!  ocean sink is calculated by HEMCO and the DMS diagnostics only includes
!  the ocean flux (this is NOT the net flux!!). 
!  If needed, we can build a simple wrapper in hcoi\_gc\_main\_mod.F90 that
!  explicitly calculates oceanic fluxes.
! \end{itemize}
!
! !INTERFACE:
!
MODULE HCOI_GC_Diagn_Mod
!
! !USES:
!
  ! GEOS-Chem diagnostic switches and arrays
  USE CMN_SIZE_Mod
#if defined( BPCH_DIAG )
  USE CMN_DIAG_Mod
  USE DIAG_Mod
  USE DIAG53_Mod
  USE DIAG56_Mod
#endif
  USE HCO_Diagn_Mod
  USE HCO_Error_Mod
  USE HCO_Interface_Mod

  IMPLICIT NONE
  PRIVATE

  ! Get parameters that define the different categories
#include "hcoi_gc_diagn_include.H"
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOI_GC_Diagn_Init
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !REMARKS:
!  This is currently a "bridge" module to provide backwards compatibility
!  with existing GEOS-Chem diagnostics.  We will eventually write all
!  diagnostics to netCDF format, but we are not quite there yet.
!
! !REVISION HISTORY:
!  04 May 2014 - C. Keller   - Initial version. 
!  11 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  11 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  28 Jul 2014 - C. Keller   - Split off from hcoio_diagn_mod.F90 and moved
!                              from HEMCO/Interface to GeosCore.
!  20 Aug 2014 - R. Yantosca - Add wrapper function GetHemcoId to simplify
!                              the process of getting the HEMCO species ID
!  20 Aug 2014 - R. Yantosca - Split code into several routines, for clarity
!  26 Aug 2014 - M. Sulprizio- Add modifications for POPs emissions diagnostics
!  23 Sep 2014 - C. Keller   - Added Hg diagnostics
!  11 Nov 2014 - C. Keller   - Added call to ESMF diagnostics.
!  22 Apr 2015 - M. Sulprizio- Now save out hydrocarbons in units kgC/m2/s
!  27 Feb 2016 - C. Keller   - Update to HEMCO v2.0
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_GC_Diagn_Init
!
! !DESCRIPTION: Subroutine HCOI\_GC\_Diagn\_Init initializes the HEMCO 
! diagnostics in GEOS-Chem. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_GC_Diagn_Init( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt
    USE HCO_State_Mod,      ONLY : HCO_GetHcoID
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  The category numbers must correspond to those in the HEMCO_Config.rc file.
!  We will have to come up with a better way of making sure that these
!  are consistent in the future.

!  CO emissions (ND29) 
!  --> Anthropogenic, biogenic, biomass and biofuel emissions are 
!      all covered in the respective sections. 
!  --> CO produced from methanol doesn't seem to be written anymore?!
!      Not filled for now.
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller   - Initial version 
!  11 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  11 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  13 Aug 2014 - R. Yantosca - Cosmetic changes
!  20 Aug 2014 - R. Yantosca - Now call wrapper function GetHemcoId to get 
!                              the HEMCO ID # for each species name
!  29 Aug 2014 - R. Yantosca - Now exit if any of the subroutines come
!                              back with RC = HCO_FAIL
!  22 Apr 2015 - M. Sulprizio- Now save out hydrocarbons in units kgC/m2/s
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: YesOrNo
    INTEGER            :: I, J,  HcoID, N,    AS
    INTEGER            :: ExtNr, Cat, Hier
    CHARACTER(LEN=31)  :: SpcName, DiagnName, Unit
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
 
    !=======================================================================
    ! HCOI_GC_DIAGN_INIT begins here!
    !=======================================================================

    ! Initialize
    RC      = HCO_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
       ' -> at HCOI_GC_Diagn_Init (in module GeosCore/hcoi_gc_diagn_mod.F90)'

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%  NOTE: Emissions for CH4 specialty simulations are passed to  %%%%
    !%%%%  global_ch4_mod.F90 via HEMCO diagnostics, and not directly   %%%%
    !%%%%  from the HEMCO state pointer.  Therefore, we need to make    %%%%
    !%%%%  sure that routine DIAGN_CH4 is outside the BPCH_DIAG #if     %%%%
    !%%%%  block.  -- Bob Yantosca (25 Jan 2018)                        %%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL Diagn_CH4( am_I_Root, Input_Opt, HcoState, ExtState, RC )

    ! Trap potential errors
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Diagn_CH4"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%  NOTE: Emissions for Hg specialty simulations are passed to   %%%%
    !%%%%  mercury_mod.F90 via HEMCO diagnostics, and not directly      %%%%
    !%%%%  from the HEMCO state pointer.  Therefore, we need to make    %%%%
    !%%%%  sure that routine DIAGN_Hg is outside the BPCH_DIAG #if      %%%%
    !%%%%  block.  -- Bob Yantosca (25 Jan 2018)                        %%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL Diagn_Hg( am_I_Root, Input_Opt, HcoState, ExtState, RC )

    ! Trap potential errors
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Diagn_Hg"'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%  NOTE: Some diagnostics for the POPs simulation do not have   %%%%
    !%%%%  any species associated with them, and thus need to be        %%%%
    !%%%%  declared as manual diagnostics.  We need to move the call    %%%%
    !%%%%  to Diagn_POPs outside of the #ifdef BPCH_DIAG block.         %%%%
    !%%%%    -- Bob Yantosca (09 Oct 2018)                              %%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL Diagn_POPs( am_I_Root, Input_Opt, HcoState, ExtState, RC )

    ! Trap potential errors
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Diagn_POPs"'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

#if defined( BPCH_DIAG )

    !=======================================================================
    ! Define manual diagnostics
    !
    !      CATEGORY                    : GEOS-CHEM DIAGNOSTICS
    ! (1 ) Rn-Pb-Be emissions          : ND01
    ! (2 ) Dust emissions              : ND06
    ! (3 ) Carbon aerosols             : ND07
    ! (4 ) Sea salt emissions          : ND08
    ! (5 ) Acetone ocean source        : ND11
    ! (6 ) Sulfur emisisons            : ND13
    ! (6 ) Biomass emissions           : ND07, ND13, ND28, ND29, ND32
    ! (7 ) NO emissions                : ND07, ND13, ND28, ND29, ND32
    ! (8 ) Biofuel emissions           : ND29, ND32, ND34
    ! (9 ) Anthropogenic emissions     : ND29, ND32, ND34
    ! (10) Biogenic emissions          : ND46
    ! (11) Lightning flash diagnostics : ND56
    ! (12) PARANOX diagnostics         : ND63
    ! (13) POPs emissions              : ND53
    !=======================================================================
    CALL Diagn_Radon   ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Dust    ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Carbon  ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_SeaSalt ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_AcetSrc ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Sulfur  ( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Biomass ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_NOsrc   ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Biofuel ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Anthro  ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Biogenic( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_LFlash  ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_ParaNOx ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

#if defined( TOMAS )
    CALL Diagn_TOMAS   ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
#endif

    !=======================================================================
    ! Define automatic diagnostics (AutoFill)
    !=======================================================================

    ! This is for testing only. Only activate if needed.
    IF ( .FALSE. ) THEN    ! Deactivated

!       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN 

          !-------------------------------------
          ! Hourly emissions
          !-------------------------------------

          ! Do for all emission species
          DO I = 1,43

             ! Get species name
             SELECT CASE ( I )
                CASE ( 1 )
                   SpcName = 'NO'
                   Unit    = 'kg/m2/s'
                CASE ( 2 )
                   SpcName = 'CO'
                   Unit    = 'kg/m2/s'
                CASE ( 3 )
                   SpcName = 'NH3'
                   Unit    = 'kg/m2/s'
                CASE ( 4 )
                   SpcName = 'SO2'
                   Unit    = 'kg/m2/s'
                CASE ( 5 )
                   SpcName = 'SO4'
                   Unit    = 'kg/m2/s'
                CASE ( 6 )
                   SpcName = 'GLYX'
                   Unit    = 'kg/m2/s'
                CASE ( 7 )
                   SpcName = 'MGLY'
                   Unit    = 'kg/m2/s'
                CASE ( 8 )
                   SpcName = 'C2H6'
                   Unit    = 'kgC/m2/s'
                CASE ( 9 )
                   SpcName = 'ALK4'
                   Unit    = 'kgC/m2/s'
                CASE ( 10 )
                   SpcName = 'ACET'
                   Unit    = 'kgC/m2/s'
                CASE ( 11 )
                   SpcName = 'MEK'
                   Unit    = 'kgC/m2/s'
                CASE ( 12 )
                   SpcName = 'ALD2'
                   Unit    = 'kgC/m2/s'
                CASE ( 13 )
                   SpcName = 'PRPE'
                   Unit    = 'kgC/m2/s'
                CASE ( 14 )
                   SpcName = 'C3H8'
                   Unit    = 'kgC/m2/s'
7                CASE ( 15 )
                   SpcName = 'CH2O'
                   Unit    = 'kg/m2/s'
                CASE ( 16 )
                   SpcName = 'BENZ'
                   Unit    = 'kgC/m2/s'
                CASE ( 17 )
                   SpcName = 'TOLU'
                   Unit    = 'kgC/m2/s'
                CASE ( 18 )
                   SpcName = 'XYLE'
                   Unit    = 'kgC/m2/s'
                CASE ( 19 )
                   SpcName = 'C2H4'
                   Unit    = 'kgC/m2/s'
                CASE ( 20 )
                   SpcName = 'C2H2'
                   Unit    = 'kgC/m2/s'
                CASE ( 21 )
                   SpcName = 'CHBr3'
                   Unit    = 'kg/m2/s'
                CASE ( 22 )
                   SpcName = 'CH2Br2'
                   Unit    = 'kg/m2/s'
                CASE ( 23 )
                   SpcName = 'BCPI'
                   Unit    = 'kg/m2/s'
                CASE ( 24 )
                   SpcName = 'BCPO'
                   Unit    = 'kg/m2/s'
                CASE ( 25 )
                   SpcName = 'OCPI'
                   Unit    = 'kg/m2/s'
                CASE ( 26 )
                   SpcName = 'OCPO'
                   Unit    = 'kg/m2/s'
                CASE ( 27 )
                   SpcName = 'RCHO'
                   Unit    = 'kg/m2/s'
                CASE ( 28 )
                   SpcName = 'MACR'
                   Unit    = 'kg/m2/s'
                CASE ( 29 )
                   SpcName = 'DMS'
                   Unit    = 'kg/m2/s'
                CASE ( 30 )
                   SpcName = 'DST1'
                   Unit    = 'kg/m2/s'
                CASE ( 31 )
                   SpcName = 'DST2'
                   Unit    = 'kg/m2/s'
                CASE ( 32 )
                   SpcName = 'DST3'
                   Unit    = 'kg/m2/s'
                CASE ( 33 )
                   SpcName = 'DST4'
                   Unit    = 'kg/m2/s'
                CASE ( 34 )
                   SpcName = 'SALA'
                   Unit    = 'kg/m2/s'
                CASE ( 35 )
                   SpcName = 'SALC'
                   Unit    = 'kg/m2/s'
                CASE ( 36 )
                   SpcName = 'Br2'
                   Unit    = 'kg/m2/s'
                CASE ( 37 )
                   SpcName = 'ISOP'
                   Unit    = 'kgC/m2/s'
                CASE ( 38 )
                   SpcName = 'Rn'
                   Unit    = 'kg/m2/s'
                CASE ( 39 )
                   SpcName = 'O3'
                   Unit    = 'kg/m2/s'
                CASE ( 40 )
                   SpcName = 'PAN'
                   Unit    = 'kg/m2/s'
                CASE ( 41 )
                   SpcName = 'HNO3'
                   Unit    = 'kg/m2/s'
                CASE ( 42 )
                   SpcName = 'MTPA'
                   Unit    = 'kgC/m2/s'
                CASE ( 43 )
                   SpcName = 'EOH'
                   Unit    = 'kgC/m2/s'
                CASE DEFAULT
                   SpcName = 'DUMMY'
             END SELECT

             HcoID = HCO_GetHcoID( TRIM(SpcName), HcoState )
             IF ( HcoID > 0 ) THEN
                CALL Diagn_Create ( am_I_Root,                          &
                                    cName     = 'EMIS_'//TRIM(SpcName), &
                                    HcoState  = HcoState,               &
                                    ExtNr     = -1,                     &
                                    Cat       = -1,                     &
                                    Hier      = -1,                     &
                                    HcoID     = HcoID,                  &
                                    SpaceDim  = 2,                      &
                                    LevIDx    = -1,                     &
                                    OutUnit   = TRIM(Unit),             &
                                    AutoFill  = 1,                      &
                                    COL       = HcoState%Diagn%HcoDiagnIDDefault,      &
                                    RC        = RC                       ) 
                IF ( RC /= HCO_SUCCESS ) RETURN
             ENDIF

             !-------------------------------------
             ! Emissions per category (NO only)
             !-------------------------------------
             IF ( TRIM(SpcName) == 'NO' .and. HcoID > 0 ) THEN

                ! There are 3 different categories
                DO J = 1, 7
                   SELECT CASE ( J )
                      CASE ( 1 )
                         DiagnName = 'EMIS_NO_ANTHRO'
                         ExtNr     = 0
                         Cat       = 1
                      CASE ( 2 )
                         DiagnName = 'EMIS_NO_AVIATION'
                         ExtNr     = 0
                         Cat       = 20
                      CASE ( 3 )
                         DiagnName = 'EMIS_NO_PARANOX'
                         ExtNr     = 102
                         Cat       = -1
                      CASE ( 4 )
                         DiagnName = 'EMIS_NO_LIGHTNING'
                         ExtNr     = 103
                         Cat       = -1
                      CASE ( 5 )
                         DiagnName = 'EMIS_NO_SOIL'
                         ExtNr     = 104
                         Cat       = -1
                      CASE ( 6 )
                         DiagnName = 'EMIS_NO_BIOMASS'
                         ExtNr     = 111
                         Cat       = -1
                      CASE ( 7 )
                         DiagnName = 'EMIS_NO_BIOFUEL'
                         ExtNr     = 0
                         Cat       = 2
                      CASE DEFAULT
                         DiagnName = 'EMIS_NO_DUMMY'
                         ExtNr     = 999
                         Cat       = 999
                   END SELECT

                   CALL Diagn_Create ( am_I_Root,                     &
                                       cName     = DiagnName,         &
                                       HcoState  = HcoState,          &
                                       ExtNr     = ExtNr,             &
                                       Cat       = Cat,               &
                                       Hier      = -1,                &
                                       HcoID     = HcoID,             &
                                       SpaceDim  = 2,                 &
                                       LevIDx    = -1,                &
                                       OutUnit   = 'kg/m2/s',         & 
                                       AutoFill  = 1,                 &
                                       COL       = HcoState%Diagn%HcoDiagnIDDefault, &
                                       RC        = RC ) 
                   IF ( RC /= HCO_SUCCESS ) RETURN

                ENDDO ! J
             ENDIF    ! NO
          ENDDO       ! I
!       ENDIF          ! fullchem
    ENDIF             ! testing toggle

    ! Leave w/ success
    RC = HCO_SUCCESS 
#endif

  END SUBROUTINE HCOI_GC_Diagn_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_Radon
!
! !DESCRIPTION: Subroutine Diagn\_Radon initializes diagnostics for the
!  Rn-Pb-Be simulation (ND01).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Radon( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Split off code from HCOI_GC_Diagn_Init into smaller routines in order to
!  make the code more manageable.
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use Rn-Pb-Be
!  03 Sep 2014 - R. Yantosca - Don't define diagnostic container for Pb
!  03 Sep 2014 - R. Yantosca - Change units from kg/m2/s to kg; also in diag3.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, N
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_RADON (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! Define ND01 diagnostics (Rn-Pb-Be emissions)
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if the Rn-Pb-Be simulation is not selected
    IF ( .not. Input_Opt%ITS_A_RnPbBe_SIM ) RETURN

    ! Define diagnostics
    IF ( ( ExtState%GC_RnPbBe > 0 ) .and. ( ND01 > 0 ) ) THEN

       ! HEMCO extension # for Rn-Pb-Be
       ExtNr = GetExtNr( HcoState%Config%ExtList, 'GC_Rn-Pb-Be' )
       IF ( ExtNr <= 0 ) THEN
          CALL HCO_Error ( 'Cannot find Rn-Pb-Be extension', RC, THISLOC=LOC )
          RETURN      
       ENDIF

       !-------------------------------------------
       ! %%%%% Rn222 %%%%%
       !-------------------------------------------
 
       ! HEMCO species ID
       HcoID = GetHemcoId( 'Rn', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'AD01_Rn_SOURCE'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState = HcoState,           &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = -1,                &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/s',            &
                          AutoFill  = 1,                 &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       !-------------------------------------------
       ! %%%%% Be7 %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'Be7', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'AD01_Be7_SOURCE'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = -1,                &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/s',            &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF
#endif

  END SUBROUTINE Diagn_Radon
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_Dust
!
! !DESCRIPTION: Subroutine Diagn\_Dust initializes diagnostics for the
!  mineral dust aerosols (ND06).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Dust( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE CMN_SIZE_Mod,       ONLY : NDSTBIN
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Split off code from HCOI_GC_Diagn_Init into smaller routines in order to
!  make the code more manageable.
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use dust
!  30 Sep 2014 - R. Yantosca - Update for TOMAS dust species
!  25 Oct 2016 - R. Yantosca - Make sure to cast INTEGER to LOGICAL values
!                              before comparing them in an IF statement
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: Is_DustDead
    LOGICAL            :: Is_DustGinoux
    INTEGER            :: ExtNr, HcoID, I, N
    CHARACTER(LEN=1)   :: ISTR1
    CHARACTER(LEN=2)   :: ISTR2
    CHARACTER(LEN=15)  :: SpcName
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_DUST (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_DUST begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if we are doing a specialty simulation w/o dust
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )   .and. &
         ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
       RETURN
    ENDIF

    ! Now use local LOGICAL variables to save ExtState%DustDead and
    ! ExtState%DustGinoux.  This will make sure these variables are 
    ! cast to LOGICAL, so that we can compare them in the same IF
    ! statement.  Otherwise GNU Fortran will choke. (bmy, 10/25/16)
    Is_DustDead   = ( ExtState%DustDead   )
    Is_DustGinoux = ( ExtState%DustGinoux )

    ! Define diagnostics if dust is used
    IF ( ( Is_DustDead .OR. Is_DustGinoux )  .AND. &
         ( ND06 > 0                       ) ) THEN

       ! Get Ext. Nr of used extension
       ExtNr = GetExtNr( HcoState%Config%ExtList, 'DustDead' )
       IF ( ExtNr <= 0 ) ExtNr = GetExtNr( HcoState%Config%ExtList, 'DustGinoux' )
       IF ( ExtNr <= 0 ) THEN
          CALL HCO_Error( 'Cannot find dust extension', RC, THISLOC=LOC )
          RETURN      
       ENDIF

       ! Do for each dust bin
       DO I = 1, NDSTBIN

#if defined( TOMAS )

          ! Get species name (i.e. DUST1 .. DUST40) for TOMAS simulatiosn
          IF ( I < 10 )  THEN
             WRITE( ISTR1,'(i1)' ) I
             SpcName = 'DUST'   // ISTR1
          ELSE
             WRITE( ISTR2,'(i2)' ) I
             SpcName = 'DUST'   // ISTR2
          ENDIF
#else

          ! Get species name (i.e. DST1 .. DST4) for non TOMAS simualtions
          WRITE( ISTR1,'(i1)' ) I
          SpcName   = 'DST'   // ISTR1

#endif

          DiagnName = 'AD06_' // TRIM( SpcName )         

          ! HEMCO species ID 
          HcoID = GetHemcoId( TRIM( SpcName ), HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create diagnostic container
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg',              &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 

       ENDDO 

       ! Add diagnostics for dust alkalinity
       IF ( Input_Opt%LDSTUP ) THEN

          ! Get Ext. Nr of used extension
          ExtNr = GetExtNr( HcoState%Config%ExtList, 'DustAlk' )
          IF ( ExtNr <= 0 ) THEN
             CALL HCO_Error( 'Cannot find dust alk extension', RC, &
                              THISLOC=LOC )
             RETURN      
          ENDIF

          ! Do for each dust bin
          DO I = 1, NDSTBIN

             ! Get species name (i.e. DSTAL1 .. DSTAL4)
             WRITE( ISTR1,'(i1)' ) I
             SpcName   = 'DSTAL' // ISTR1
             DiagnName = 'AD06_' // TRIM( SpcName )

             ! HEMCO species ID 
             HcoID = GetHemcoId( TRIM( SpcName ), HcoState, LOC, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN

             ! Create diagnostic container
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = -1,                &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg',              &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDDO
       ENDIF

    ENDIF
#endif

  END SUBROUTINE Diagn_Dust
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_Carbon
!
! !DESCRIPTION: Subroutine Diagn\_Carbon initializes diagnostics for the
!  carbon aerosols (ND07).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Carbon( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      only : Ind_
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Biomass diagnostics (ND28) are defined in routine Diagn_Biomass.
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use carbon
!  16 Jun 2016 - C. Miller   - Now define species ID's with Ind_ function 
!  27 Mar 2017 - M. Sulprizio- Make anthropogenic emissions diagnostics 3D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, Cat, N, I, J
    CHARACTER(LEN=31)  :: SpcName, SrcName, DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_CARBON (hcoi_gc_diagn_mod.F90)'
    INTEGER            :: id_POA1, id_POG1
    INTEGER            :: SpaceDim

    !=======================================================================
    ! DIAGN_CARBON begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if we are doing a specialty simulation w/o carbon aerosols
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )   .and. &
         ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
       RETURN
    ENDIF
    
    ! Define advected species ID's
    id_POA1 = Ind_('POA1','A')
    id_POG1 = Ind_('POG1','A')

    ! Define diagnostics
    IF ( ND07 > 0 .AND. Input_Opt%LCARB ) THEN

       ! Do for all species
       DO I = 1, 4
          
          ! Get species name
          SELECT CASE ( I ) 
             CASE ( 1 ) 
                SpcName = 'BCPI'
             CASE ( 2 ) 
                SpcName = 'BCPO'
             CASE ( 3 )
                SpcName = 'OCPI'
                IF ( id_POA1 > 0 ) SpcName = 'POA1'
             CASE ( 4 ) 
                SpcName = 'OCPO'
                IF ( id_POG1 > 0 ) SpcName = 'POG1'
          END SELECT

          ! HEMCO species ID
          HcoID = GetHemcoId( SpcName, HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Do for all sources
          DO J = 1, 2

             SELECT CASE ( J )
                CASE ( 1 )
                   SrcName = 'ANTHRO'
                   Cat     = CATEGORY_ANTHRO
                   SpaceDim= 3
                CASE ( 2 )
                   SrcName = 'BIOFUEL'
                   Cat     = CATEGORY_BIOFUEL 
                   SpaceDim= 2
             END SELECT

             !-------------------------------------------
             ! %%%%% DEFINE DIAGNOSTICS %%%%%% 
             ! -------------------------------------------

             ! Set DiagnName
             DiagnName = "AD07_"//TRIM(SpcName)//"_"//TRIM(SrcName)

             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = 0,                 &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = SpaceDim,          &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 

          ENDDO
       ENDDO
    ENDIF 
#endif

  END SUBROUTINE Diagn_Carbon
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_SeaSalt
!
! !DESCRIPTION: Subroutine Diagn\_SeaSalt initializes diagnostics for the
!  sea salt aerosols (ND08).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_SeaSalt( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Split off code from HCOI_GC_Diagn_Init into smaller routines in order to
!  make the code more manageable.
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use sea salt
!  09 Jul 2015 - E. Lundgren - Add marine organic aerosols
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, ExtNrSS, ExtNrMPOA
    INTEGER            :: HcoID, I, N, NSALT
    CHARACTER(LEN=15)  :: SpcName
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_SEASALT (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_SEASALT begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if we are doing a specialty simulation w/o sea salt
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )   .and. &
         ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
       RETURN
    ENDIF

    ! Define diagnostics
    IF ( ND08 > 0 .AND. Input_Opt%LSSALT .AND. ( ExtState%SeaSalt > 0 ) ) THEN

       ! Get HEMCO extension # for SeaSalt
       ExtNrSS = GetExtNr( HcoState%Config%ExtList, 'SeaSalt' )
       IF ( ExtNrSS <= 0 ) THEN
          CALL HCO_Error( 'Cannot find extension SeaSalt', RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Get HEMCO extension # for marine organic aerosols and
       ! set number of seasalt tracers
       IF ( Input_Opt%LMPOA ) THEN
          ExtNrMPOA = GetExtNr( HcoState%Config%ExtList, 'MarinePOA' )
          IF ( ExtNrMPOA <= 0 ) THEN
             CALL HCO_Error( 'Cannot find extension MarinePOA', RC,  &
                             THISLOC=LOC )
             RETURN
          ENDIF
          NSALT = 4
       ELSE
          NSALT = 2
       ENDIF          
       
       ! I=1 is SALA, I=2 is SALC, I=3 is MOPO, I=4 is MOPI
       DO I = 1,NSALT

          ! Pick the proper species, diagnostic name, and extension #
          SELECT CASE( I )
             CASE( 1 )
                ExtNr     = ExtNrSS
                SpcName   = 'SALA'
                DiagnName = 'AD08_SALA'
             CASE( 2 )
                ExtNr     = ExtNrSS
                SpcName   = 'SALC'
                DiagnName = 'AD08_SALC'
             CASE( 3 )
                ExtNr     = ExtNrMPOA
                SpcName   = 'MOPO'
                DiagnName = 'AD08_MOPO'
             CASE( 4 )
                ExtNr     = ExtNrMPOA
                SpcName   = 'MOPI'
                DiagnName = 'AD08_MOPI'
          END SELECT

          ! HEMCO species ID 
          HcoID = GetHemcoId( TRIM( SpcName ), HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create the diagnostic
          CALL Diagn_Create ( am_I_Root,                     & 
                              HcoState  = HcoState,          &
                              cName     = TRIM( DiagnName ), &
                              ExtNr     = ExtNr,             &
                              Cat       = -1,                &
                              Hier      = -1,                &
                              HcoID     = HcoID,             &
                              SpaceDim  = 2,                 &
                              LevIDx    = -1,                &
                              OutUnit   = 'kg',              &
                              COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                              AutoFill  = 1,                 &
                              RC        = RC                  ) 

          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDDO 
    ENDIF
#endif

  END SUBROUTINE Diagn_SeaSalt
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_AcetSrc
!
! !DESCRIPTION: Subroutine Diagn\_SeaSalt initializes diagnostics for the
!  acetone ocean sources (ND11).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_AcetSrc( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  ACETONE (AD11)
!  --> 3 manually defined diagnostics in MEGAN (defined in ND46)
!  --> 1 automatically filled diagnostics in Seaflux
!  --> Ocean sink is passed to drydep and not explicitly written out!!
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use acetone
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, N
    CHARACTER(LEN=15)  :: SpcName
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_ACETSRC (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_ACETSRC begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if we are doing a specialty simulation w/o acetone
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN 
       RETURN
    ENDIF

    ! Define diagnostics
    IF ( ND11 > 0 .OR. ND46 > 0 ) THEN 

       ! Get extension # for SeaFlux 
       ExtNr = GetExtNr( HcoState%Config%ExtList, 'SeaFlux' )

       IF ( ExtNr <= 0 ) THEN

          ! SeaFlux not used, so print warning
          CALL HCO_Warning( 'Cannot find extension SeaFlux', RC, THISLOC=LOC )

       ELSE

          ! HEMCO species ID
          HcoID = GetHemcoId( 'ACET', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create diagnostics. Set AutoFill to on for this diagnostics.
          DiagnName = 'AD11_OCEAN_SOURCE'
          CALL Diagn_Create ( am_I_Root,                     & 
                              HcoState  = HcoState,          &
                              cName     = TRIM( DiagnName ), &
                              ExtNr     = ExtNr,             &
                              Cat       = -1,                &
                              Hier      = -1,                &
                              HcoID     = HcoID,             &
                              SpaceDim  = 2,                 &
                              LevIDx    = -1,                &
                              OutUnit   = 'kgC/m2/s',         &
                              COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                              AutoFill  = 1,                 &
                              RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF
    ENDIF
#endif

  END SUBROUTINE Diagn_AcetSrc
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_Sulfur
!
! !DESCRIPTION: Subroutine Diagn\_Sulfur initializes diagnostics for the
!  sulfur sources (ND13).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Sulfur( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Sulfur emissions (ND13) 
!  --> For DMS, only positive flux is diagnosed
!  --> Don't diagnose biofuel as most inventory include it w/ anthro
!  --> Volcano emissions are lumped (eruptive + noneruptive)
!  ==> BIOMASS diagnostics (ND28) are defined in routine DIAGN_BIOMASS
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use sulfur
!  23 Feb 2015 - C. Keller   - Split volcano into eruptive and degassing.
!  27 Mar 2017 - M. Sulprizio- Make anthropogenic emissions diagnostics 3D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, N
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_SULFUR (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_SULFUR begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if we are doing a specialty simulation w/o carbon aerosols
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )   .and. &
         ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
       RETURN
    ENDIF
    
    ! Define diagnostics
    IF ( ND13 > 0 ) THEN

       !-------------------------------------------
       ! %%%%% DMS %%%%%
       !
       ! As we did for acetone, only keep track 
       ! of flux from ocean.  Deposition from 
       ! atmosphere is handled by drydep.
       !-------------------------------------------

       ! HEMCO extension # for SeaFlux
       ExtNr = GetExtNr( HcoState%Config%ExtList, 'SeaFlux' )

       IF ( ExtNr <= 0 ) THEN

          ! SeaFlux not found, print warning
          MSG = 'Cannot find SeaFlux extension - ' // &
                 'no DMS diagnostics will be written!'
          CALL HCO_Warning( MSG, RC, THISLOC=LOC )

       ELSE

          ! HEMCO species ID
          HcoID = GetHemcoId( 'DMS', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create diagnostic container
          DiagnName = 'AD13_DMS_OCEAN_SOURCE'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg',              &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% SO2 %%%%%
       !-------------------------------------------
       ExtNr = 0

       ! HEMCO species ID
       HcoID = GetHemcoId( 'SO2', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! ... from aircrafts ...
       DiagnName = 'AD13_SO2_AIRCRAFT'
       CALL Diagn_Create( am_I_Root,                     &   
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_AIRCRAFT, &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! ... anthropogenic ...
       DiagnName = 'AD13_SO2_ANTHROPOGENIC'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! ... biofuel ...
       DiagnName = 'AD13_SO2_BIOFUEL'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_BIOFUEL,  &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! ... from volcanoes (eruptive) ...
       DiagnName = 'AD13_SO2_VOLCANO_ERUPT'
       CALL Diagn_Create( am_I_Root,                           & 
                          HcoState  = HcoState,                &
                          cName     = TRIM( DiagnName ),       &
                          ExtNr     = ExtNr,                   &
                          Cat       = CATEGORY_VOLCANO_ERUPT,  &
                          Hier      = -1,                      &
                          HcoID     = HcoID,                   &
                          SpaceDim  = 3,                       &
                          LevIDx    = -1,                      &
                          OutUnit   = 'kg',                    &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,        &
                          AutoFill  = 1,                       &
                          RC        = RC                        )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! ... from volcanoes (non-eruptive / degassing) ...
       DiagnName = 'AD13_SO2_VOLCANO_DEGAS'
       CALL Diagn_Create( am_I_Root,                           & 
                          HcoState  = HcoState,                &
                          cName     = TRIM( DiagnName ),       &
                          ExtNr     = ExtNr,                   &
                          Cat       = CATEGORY_VOLCANO_DEGAS,  &
                          Hier      = -1,                      &
                          HcoID     = HcoID,                   &
                          SpaceDim  = 3,                       &
                          LevIDx    = -1,                      &
                          OutUnit   = 'kg',                    &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,        &
                          AutoFill  = 1,                       &
                          RC        = RC                        )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! ... from ships ...
       DiagnName = 'AD13_SO2_SHIP'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_SHIP,     &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       !-------------------------------------------
       ! %%%%% NH3 %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'NH3', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! ... anthropogenic ...
       ExtNr     = 0
       DiagnName = 'AD13_NH3_ANTHROPOGENIC'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! ... biofuel ...
       ExtNr     = 0
       DiagnName = 'AD13_NH3_BIOFUEL'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_BIOFUEL,  &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! ... natural
       ExtNr     = 0
       DiagnName = 'AD13_NH3_NATURAL'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_NATURAL,  &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       !-------------------------------------------
       ! %%%%% SO4 %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'SO4', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! ... anthropogenic ...
       ExtNr     = 0
       DiagnName = 'AD13_SO4_ANTHROPOGENIC' 
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! ... biofuel ...
       ExtNr     = 0
       DiagnName = 'AD13_SO4_BIOFUEL' 
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_BIOFUEL,  &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF
#endif

  END SUBROUTINE Diagn_Sulfur
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_Biomass
!
! !DESCRIPTION: Subroutine Diagn\_Biomass initializes diagnostics for the
!  biomass emissions species (ND13, ND28, ND29, ND32).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Biomass( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCO_State_Mod,      ONLY : HCO_GetHcoID
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : Ind_
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Biomass burning emissions (ND28, ND07, ND13, ND29, ND32) 
!  ==> write one single biomass burning diagnostics per species.
!  ==> Biomass buring comes from GFED or FINN inventory. If none
!      of those inventories is used, it is assumed that biomass
!      burning has category 3.
!       --> NOTE: 3 is the same category as natural source NH3!
!  ==> Diagnostics are returned in kg/m2/s.
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use biomass
!  17 Jun 2016 - R. Yantosca - Now define species ID's with Ind_ function
!  10 Mar 2017 - M. Sulprizio- Change SpaceDim from 2 to 3 for emitting 35% of
!                              biomass burning source into the free troposphere
!  24 Apr 2017 - M. Sulprizio- Change SpaceDim back to 2 for now

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: Cat, ExtNr, HcoID, N, N_CO
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_BIOMASS (hcoi_gc_diagn_mod.F90)'

    ! CO tracer names
    INTEGER, PARAMETER :: N_BIOM_CO             = 7
    CHARACTER(LEN=7)   :: CO_Tracers(N_BIOM_CO) =          &
         (/ 'CO     ', 'CObbam ', 'CObbaf ', 'CObbas ' ,   &
            'CObboc ', 'CObbeu ', 'CObboth'              /) 

    !=======================================================================
    ! DIAGN_BIOMASS begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if we are doing a specialty simulation w/o biomass
    IF ( Input_Opt%ITS_A_POPS_SIM   ) RETURN
    IF ( Input_Opt%ITS_A_RnPbBe_SIM ) RETURN
    IF ( Input_Opt%ITS_A_TAGO3_SIM  ) RETURN

    ! First test if GFED is used.  If not, then test if FINN is used.
    ! If not, then use extension # 0 and the default biomass category.
    Cat   = -1
    ExtNr = GetExtNr( HcoState%Config%ExtList, 'GFED' )
    IF ( ExtNr <= 0 ) ExtNr = GetExtNr( HcoState%Config%ExtList, 'FINN' )
    IF ( ExtNr <= 0 ) THEN
       ExtNr = 0
       Cat   = CATEGORY_BIOMASS
    ENDIF
    
    ! ND28 only, for VOC species
    IF ( ND28 > 0 ) THEN

       !-------------------------------------------
       ! %%%%% Biomass ALK4 %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'ALK4', HcoState )
       IF ( HcoID > 0 ) THEN  

          ! Create diagnostic container
          DiagnName = 'BIOMASS_ALK4'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kgC/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass ACET %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'ACET', HcoState )
       IF ( HcoID > 0 ) THEN  

          ! Create diagnostic container
          DiagnName = 'BIOMASS_ACET'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kgC/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass MEK %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'MEK', HcoState )
       IF ( HcoID > 0 ) THEN  

          ! Create diagnostic container
          DiagnName = 'BIOMASS_MEK'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kgC/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass ALD2 %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'ALD2', HcoState )
       IF ( HcoID > 0 ) THEN  

          ! Create diagnostic container
          DiagnName = 'BIOMASS_ALD2'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kgC/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass PRPE %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'PRPE', HcoState )
       IF ( HcoID > 0 ) THEN  

          ! Create diagnostic container
          DiagnName = 'BIOMASS_PRPE'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kgC/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass C3H8 %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'C3H8', HcoState )
       IF ( HcoID > 0 ) THEN  

          ! Create diagnostic container
          DiagnName = 'BIOMASS_C3H8'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kgC/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass CH2O %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'CH2O', HcoState )
       IF ( HcoID > 0 ) THEN  

          ! Create diagnostic container
          DiagnName = 'BIOMASS_CH2O'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass C2H6 %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'C2H6', HcoState )
       IF ( HcoID > 0 ) THEN  
          ! Create diagnostic container
          DiagnName = 'BIOMASS_C2H6'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kgC/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
  
       !-------------------------------------------
       ! %%%%% Biomass CH4 %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'CH4', HcoState )
       IF ( HcoID > 0 ) THEN  
          ! Create diagnostic container
          DiagnName = 'BIOMASS_CH4'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kgC/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass NH3 %%%%%
       !-------------------------------------------
       HcoID = HCO_GetHcoID( 'NH3', HcoState )
       IF ( HcoID > 0 ) THEN  
          DiagnName = 'BIOMASS_NH3'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass PAN %%%%%
       !-------------------------------------------
       HcoID = HCO_GetHcoID( 'PAN', HcoState )
       IF ( HcoID > 0 ) THEN  
          DiagnName = 'BIOMASS_PAN'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass HNO3 %%%%%
       !-------------------------------------------
       HcoID = HCO_GetHcoID( 'HNO3', HcoState )
       IF ( HcoID > 0 ) THEN  
          DiagnName = 'BIOMASS_HNO3'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass MTPA %%%%%
       !-------------------------------------------
       HcoID = HCO_GetHcoID( 'MTPA', HcoState )
       IF ( HcoID > 0 ) THEN  
          DiagnName = 'BIOMASS_MTPA'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kgC/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass BENZ %%%%%
       !-------------------------------------------
       HcoID = HCO_GetHcoID( 'BENZ', HcoState )
       IF ( HcoID > 0 ) THEN  
          DiagnName = 'BIOMASS_BENZ'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kgC/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass TOLU %%%%%
       !-------------------------------------------
       HcoID = HCO_GetHcoID( 'TOLU', HcoState )
       IF ( HcoID > 0 ) THEN  
          DiagnName = 'BIOMASS_TOLU'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kgC/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass XYLE %%%%%
       !-------------------------------------------
       HcoID = HCO_GetHcoID( 'XYLE', HcoState )
       IF ( HcoID > 0 ) THEN  
          DiagnName = 'BIOMASS_XYLE'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kgC/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass EOH %%%%%
       !-------------------------------------------
       HcoID = HCO_GetHcoID( 'EOH', HcoState )
       IF ( HcoID > 0 ) THEN  
          DiagnName = 'BIOMASS_EOH'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kgC/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass MGLY %%%%%
       !-------------------------------------------
       HcoID = HCO_GetHcoID( 'MGLY', HcoState )
       IF ( HcoID > 0 ) THEN  
          DiagnName = 'BIOMASS_MGLY'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       !-------------------------------------------
       ! %%%%% Biomass NAP %%%%%
       !-------------------------------------------
       HcoID = HCO_GetHcoID( 'NAP', HcoState )
       IF ( HcoID > 0 ) THEN  
          DiagnName = 'BIOMASS_NAP'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

    ENDIF ! ND28 only

    
    !-------------------------------------------
    ! %%%%% Biomass SO2 %%%%%
    !-------------------------------------------
    IF ( ND28 > 0 .OR. ND13 > 0 ) THEN

       HcoID = HCO_GetHcoID( 'SO2', HcoState )
       IF ( HcoID > 0 ) THEN
          ! Create diagnostic container
          DiagnName = 'BIOMASS_SO2'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
    ENDIF

    ! BC and OC
    IF ( ND28 > 0 .OR. ( ND07 > 0 .AND. Input_Opt%LCARB ) ) THEN

       ! BC and OC are only defined for either
       ! full-chemistry or offline aerosol simulations
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or.             &
            Input_Opt%ITS_AN_AEROSOL_SIM ) THEN 

          !-------------------------------------------
          ! %%%%% Biomass POG1 and POG2 %%%%% krt, 8/24/17
          !-------------------------------------------
          HcoID = HCO_GetHcoID( 'POG1', HcoState )
          IF ( HcoID > 0 ) THEN  
             DiagnName = 'BIOMASS_POG1'
             CALL Diagn_Create( am_I_Root,                     & 
                  HcoState  = HcoState,          &
                  cName     = TRIM( DiagnName ), &
                  ExtNr     = ExtNr,             &
                  Cat       = Cat,               &
                  Hier      = -1,                &
                  HcoID     = HcoID,             &
                  SpaceDim  = 2,                 &
                  LevIDx    = -1,                &
                  OutUnit   = 'kg/m2/s',         &
                  COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                  AutoFill  = 1,                 &
                  RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          HcoID = HCO_GetHcoID( 'POG2', HcoState )
          IF ( HcoID > 0 ) THEN  
             DiagnName = 'BIOMASS_POG2'
             CALL Diagn_Create( am_I_Root,                     & 
                  HcoState  = HcoState,          &
                  cName     = TRIM( DiagnName ), &
                  ExtNr     = ExtNr,             &
                  Cat       = Cat,               &
                  Hier      = -1,                &
                  HcoID     = HcoID,             &
                  SpaceDim  = 2,                 &
                  LevIDx    = -1,                &
                  OutUnit   = 'kg/m2/s',         &
                  COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                  AutoFill  = 1,                 &
                  RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          !----------------------------------------
          ! %%%%% Biomass OC %%%%%
          !----------------------------------------
          HcoID = GetHemcoId( 'OCPI', HcoState, LOC, RC )
          IF ( HcoID > 0 ) THEN
             DiagnName = 'BIOMASS_OCPI'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          HcoID = GetHemcoId( 'OCPO', HcoState, LOC, RC )
          IF ( HcoID > 0 ) THEN
             DiagnName = 'BIOMASS_OCPO'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          !----------------------------------------
          ! %%%%% Biomass BC %%%%%
          !----------------------------------------
          
          ! HEMCO species ID
          HcoID = GetHemcoId( 'BCPI', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          
          ! Create diagnostic container
          DiagnName = 'BIOMASS_BCPI'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! HEMCO species ID
          HcoID = GetHemcoId( 'BCPO', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          
          ! Create diagnostic container
          DiagnName = 'BIOMASS_BCPO'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
    ENDIF

    !----------------------------------------------
    ! %%%%% Biomass CO %%%%%
    !----------------------------------------------
    IF ( ND28 > 0 .OR. ND29 > 0 ) THEN

       ! CO is only defined for full chemistry and tagged CO simulations
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or. &
            Input_Opt%ITS_A_TAGCO_SIM    ) THEN

          ! Loop over tagged CO tracers if necessary
          IF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN
             N_CO = N_BIOM_CO
          ELSE
             N_CO = 1
          ENDIF

          ! Loop over all CO tracers
          DO N = 1, N_CO

             ! Pick the various category names
             SELECT CASE( TRIM( CO_Tracers(N) ) )
                CASE( 'CO'      )
                   HcoId     =  GetHemcoId( 'CO', HcoState, LOC, RC )
                   DiagnName = 'BIOMASS_CO'
                CASE( 'CObbam'  ) 
                   HcoId     =  GetHemcoId( 'CObbam', HcoState, LOC, RC )
                   DiagnName = 'BIOMASS_TAGCO_USA'
                CASE( 'CObbaf'  ) 
                   HcoId     =  GetHemcoId( 'CObbaf', HcoState, LOC, RC )
                   DiagnName = 'BIOMASS_TAGCO_AFRICA'
                CASE( 'CObbas'  ) 
                   HcoId     =  GetHemcoId( 'CObbas', HcoState, LOC, RC )
                   DiagnName = 'BIOMASS_TAGCO_ASIA'
                CASE( 'CObboc'  ) 
                   HcoId     =  GetHemcoId( 'CObboc', HcoState, LOC, RC )
                   DiagnName = 'BIOMASS_TAGCO_OCEANIA'
                CASE( 'CObbeu'  ) 
                   HcoId     =  GetHemcoId( 'CObbeu', HcoState, LOC, RC )
                   DiagnName = 'BIOMASS_TAGCO_EUROPE'
                CASE( 'CObboth' ) 
                   HcoId     =  GetHemcoId( 'CObboth', HcoState, LOC, RC )
                   DiagnName = 'BIOMASS_TAGCO_OTHER'
                CASE DEFAULT
                   HcoId     = -1
                   DiagnName = ''
             END SELECT
             
             ! Define the diagnostic catetory if the HEMCO id is found
             IF ( HcoId > 0 ) THEN
                CALL Diagn_Create( am_I_Root,                       & 
                                   HcoState  = HcoState,            &
                                   cName     = TRIM( DiagnName ),   &
                                   ExtNr     = ExtNr,               &
                                   Cat       = Cat,                 &
                                   Hier      = -1,                  &
                                   HcoID     = HcoID,               &
                                   SpaceDim  = 2,                   &
                                   LevIDx    = -1,                  &
                                   OutUnit   = 'kg/m2/s',           &
                                   COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                   AutoFill  = 1,                   &
                                   RC        = RC                  ) 
             ENDIF
          ENDDO
       ENDIF
    ENDIF

    !----------------------------------------------
    ! %%%%% Biomass NO %%%%%
    !----------------------------------------------
    IF ( ND28 > 0 .OR. ND32 > 0 ) THEN

       ! NO is only defined for the full chemistry simulation
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

          ! HEMCO species ID
          HcoID = GetHemcoId( 'NO', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create diagnostic container
          DiagnName = 'BIOMASS_NO'
          CALL Diagn_Create( am_I_Root,                   & 
                             HcoState  = HcoState,        &
                             cName     = TRIM(DiagnName), &
                             ExtNr     = ExtNr,           &
                             Cat       = Cat,             &
                             Hier      = -1,              &
                             HcoID     = HcoID,           &
                             SpaceDim  = 2,               &
                             LevIDx    = -1,              &
                             OutUnit   = 'kg/m2/s',       &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,&
                             AutoFill  = 1,               &
                             RC        = RC                ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
    ENDIF

! Now moved to Diagn_Hg
!    !----------------------------------------------
!    ! %%%%% Biomass Hg0 %%%%%
!    !----------------------------------------------
!
!    ! Do only if Hg0 is defined ... 
!    HcoID = HCO_GetHcoID( 'Hg0', HcoState )
!    IF ( HcoID > 0 ) THEN
!
!       ! Create diagnostic container
!       DiagnName = 'BIOMASS_HG0'
!       CALL Diagn_Create( am_I_Root,                   & 
!                          HcoState  = HcoState,        &
!                          cName     = TRIM(DiagnName), &
!                          ExtNr     = ExtNr,           &
!                          Cat       = Cat,             &
!                          Hier      = -1,              &
!                          HcoID     = HcoID,           &
!                          SpaceDim  = 2,               &
!                          LevIDx    = -1,              &
!                          OutUnit   = 'kg/m2/s',       &
!                          COL       = HcoState%Diagn%HcoDiagnIDManual,&
!                          AutoFill  = 1,               &
!                          RC        = RC                ) 
!       IF ( RC /= HCO_SUCCESS ) RETURN
!    ENDIF
#endif

  END SUBROUTINE Diagn_Biomass
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_NOsrc
!
! !DESCRIPTION: Subroutine Diagn\_NOsrc initializes diagnostics for the
!  NO emissions (ND32)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_NOsrc( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  NO emissions (ND32) 
!  --> Anthropogenic, biogenic, biomass and biofuel emissions are 
!      all covered in the respective sections. 
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use NO
!  22 Jan 2015 - M. Yannetti - Corrected typo in LOC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: YesOrNo
    INTEGER            :: Cat, ExtNr, HcoID, N
    CHARACTER(LEN=1)   :: ISTR
    CHARACTER(LEN=15)  :: SpcName
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_NOSRC (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_NOSRC begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if we are doing a specialty simulation w/o NO
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

    ! Only if diagnostics are defined ...
    IF ( ND32 > 0 ) THEN

       ! Extension number
       ExtNr = 0
   
       ! HEMCO species ID
       HcoID = GetHemcoId( 'NO', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
   
       !----------------------------------------------
       ! %%%%% Aircraft NO %%%%%
       !----------------------------------------------
       DiagnName = 'AIRCRAFT_NO'
       CALL Diagn_Create( am_I_Root,                      & 
                          HcoState  = HcoState,           &
                          cName     = TRIM( DiagnName ),  &
                          ExtNr     = ExtNr,              &
                          Cat       = CATEGORY_AIRCRAFT,  &
                          Hier      = -1,                 &
                          HcoID     = HcoID,              &
                          SpaceDim  = 3,                  &
                          LevIDx    = -1,                 &
                          OutUnit   = 'kg/m2/s',          &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,   &
                          AutoFill  = 1,                  &
                          RC        = RC                   )
       IF ( RC /= HCO_SUCCESS ) RETURN 
   
   
       !----------------------------------------------
       ! %%%%% Ship NO %%%%%
       !
       ! ==> Only define if ParaNOx is not used. 
       !     SHIP_NO from ParaNOx is defined in ND63.
       !----------------------------------------------
       Cat   = -1
       ExtNr = GetExtNr( HcoState%Config%ExtList, 'ParaNOx' )
   
       IF ( ExtNr <= 0 ) THEN
          ExtNr     = 0
          Cat       = CATEGORY_SHIP
          DiagnName = 'SHIP_NO'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF
   
       !----------------------------------------------
       ! %%%%% Lightning NO %%%%%
       !
       ! ==> Only define if LightNox is turned on
       !----------------------------------------------
       Cat   = -1
       ExtNr = GetExtNr( HcoState%Config%ExtList, 'LightNOx')
       IF ( ExtNr > 0 ) THEN
          DiagnName = 'LIGHTNING_NO'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 3,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF
         
      
       !----------------------------------------------
       ! %%%%% Soil and Fertilizer NO %%%%%
       !
       ! ==> Only define if SoilNox is turned on
       !----------------------------------------------
       Cat   = -1
       ExtNr = GetExtNr( HcoState%Config%ExtList, 'SoilNOx')
       IF ( ExtNr > 0 ) THEN
   
          ! %%%%%% Soil NO %%%%%%
          DiagnName = 'SOIL_NO'
          CALL Diagn_Create ( am_I_Root,                     & 
                              HcoState  = HcoState,          &
                              cName     = TRIM( DiagnName ), &
                              ExtNr     = ExtNr,             &
                              Cat       = Cat,               &
                              Hier      = -1,                &
                              HcoID     = HcoID,             &
                              SpaceDim  = 2,                 &
                              LevIDx    = -1,                &
                              OutUnit   = 'kg/m2/s',         &
                              COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                              AutoFill  = 1,                 &
                              RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN 
   
          ! %%%%%% Fertilizer NO %%%%%%
          CALL GetExtOpt( HcoState%Config, ExtNr, &
             'Use fertilizer NOx', OptValBool=YesOrNo, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
   
          IF ( YesOrNo .eqv. .FALSE. ) THEN
             MSG = 'Fertilizer NOx disabled - diagnostics will be zero!'
             CALL HCO_Warning( MSG, RC, THISLOC=LOC )
          ENDIF
   
          DiagnName = 'FERTILIZER_NO'
          CALL Diagn_Create ( am_I_Root,                     & 
                              HcoState  = HcoState,          &
                              cName     = TRIM( DiagnName ), &
                              ExtNr     = ExtNr,             &
                              Cat       = Cat,               &
                              Hier      = -1,                &
                              HcoID     = HcoID,             &
                              SpaceDim  = 2,                 &
                              LevIDx    = -1,                &
                              OutUnit   = 'kg/m2/s',         &
                              COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                              AutoFill  = 0,                 &
                              RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF
  
    ENDIF !ND32
#endif

  END SUBROUTINE Diagn_NOsrc
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_Biofuel
!
! !DESCRIPTION: Subroutine Diagn\_Biofuel initializes diagnostics for the
!  biofuel sources (ND34).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Biofuel( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCO_State_Mod,      ONLY : HCO_GetHcoId
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Biofuel emissions (ND34, ND29, ND32)
!  ==> write one single biofuel emissions diagnostics per species.
!  ==> most inventories include biofuel emissions in the anthrop.
!      sector. For explicit biofuel emissions, assume they are
!      assigned category 3 in the HEMCO configuration file.
!  ==> Diagnostics are returned in kg/m2/s.
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use biofuels
!  22 Apr 2015 - M. Sulprizio- Now save out hydrocarbons in units kgC/m2/s
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, N, I
    CHARACTER(LEN=31)  :: DiagnName, SpcName, Unit
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_BIOFUEL (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_BIOFUEL begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if we are doing a specialty simulation w/o biofuels
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) RETURN
    IF ( Input_Opt%ITS_A_POPS_SIM    ) RETURN
    IF ( Input_Opt%ITS_A_TAGO3_SIM   ) RETURN
    IF ( Input_Opt%ITS_A_RnPbBe_SIM  ) RETURN

    ! Extension number
    ExtNr = 0
 
    ! ND34 only
    IF ( ND34 > 0 ) THEN

       ! Loop over speices
       DO I = 1, 14

          ! Select species
          SELECT CASE ( I ) 
             CASE ( 1 )
                SpcName = 'SO2'
                Unit    = 'kg/m2/s'
             CASE ( 2 )
                SpcName = 'NH3'
                Unit    = 'kg/m2/s'
             CASE ( 3 )
                SpcName = 'ALK4'
                Unit    = 'kgC/m2/s'
             CASE ( 4 )
                SpcName = 'ALD2'
                Unit    = 'kgC/m2/s'
             CASE ( 5 )
                SpcName = 'ACET'
                Unit    = 'kgC/m2/s'
             CASE ( 6 )
                SpcName = 'MEK'
                Unit    = 'kgC/m2/s'
             CASE ( 7 )
                SpcName = 'PRPE'
                Unit    = 'kgC/m2/s'
             CASE ( 8 )
                SpcName = 'C2H6'
                Unit    = 'kgC/m2/s'
             CASE ( 9 )
                SpcName = 'C3H8'
                Unit    = 'kgC/m2/s'
             CASE ( 10)
                SpcName = 'CH2O'
                Unit    = 'kg/m2/s'
             CASE ( 11)
                SpcName = 'BENZ'
                Unit    = 'kgC/m2/s'
             CASE ( 12)
                SpcName = 'TOLU'
                Unit    = 'kgC/m2/s'
             CASE ( 13)
                SpcName = 'XYLE'
                Unit    = 'kgC/m2/s'
             CASE ( 14)
                SpcName = 'EOH'
                Unit    = 'kgC/m2/s'
             CASE DEFAULT
                SpcName = 'DUMMY'
          END SELECT

          ! Check for species ID
          HcoID = HCO_GetHcoId( TRIM(SpcName), HcoState )
          IF ( HcoID > 0 ) THEN 
             ! Create diagnostic container
             DiagnName = 'BIOFUEL_'//TRIM(SpcName)
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = CATEGORY_BIOFUEL,  &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = TRIM(Unit),        &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF
       ENDDO !I
    ENDIF !ND34   

    !----------------------------------------------
    ! %%%%% Biofuel NO %%%%%
    !----------------------------------------------
    IF ( ND34 > 0 .OR. ND32 > 0 ) THEN

       ! NO is only defined for the fullchem simulation
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

         ! HEMCO species ID
         HcoID = GetHemcoId( 'NO', HcoState, LOC, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
  
         ! Create diagnostic container
         DiagnName = 'BIOFUEL_NO'
         CALL Diagn_Create( am_I_Root,                     & 
                            HcoState  = HcoState,          &
                            cName     = TRIM( DiagnName ), &
                            ExtNr     = ExtNr,             &
                            Cat       = CATEGORY_BIOFUEL,  &
                            Hier      = -1,                &
                            HcoID     = HcoID,             &
                            SpaceDim  = 2,                 &
                            LevIDx    = -1,                &
                            OutUnit   = 'kg/m2/s',         &
                            COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                            AutoFill  = 1,                 &
                            RC        = RC                  ) 
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF
   ENDIF
      
   !-------------------------------------------
   ! %%%%% Biofuel CO %%%%%
   !-------------------------------------------
   IF ( ND34 > 0 .OR. ND29 > 0 ) THEN

      ! CO is only defined for fullchem and tagged CO simulations
      IF ( Input_Opt%ITS_A_FULLCHEM_SIM  .or.              &
           Input_Opt%ITS_A_TAGCO_SIM    ) THEN

         ! HEMCO species ID
         HcoID = GetHemcoId( 'CO', HcoState, LOC, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
   
         ! Create diagnostic container
         DiagnName = 'BIOFUEL_CO'
         CALL Diagn_Create( am_I_Root,                     & 
                            HcoState  = HcoState,          &
                            cName     = TRIM( DiagnName ), &
                            ExtNr     = ExtNr,             &
                            Cat       = CATEGORY_BIOFUEL,  &
                            Hier      = -1,                &
                            HcoID     = HcoID,             &
                            SpaceDim  = 2,                 &
                            LevIDx    = -1,                &
                            OutUnit   = 'kg/m2/s',         &
                            COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                            AutoFill  = 1,                 &
                            RC        = RC                  ) 
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF
   ENDIF
#endif

  END SUBROUTINE Diagn_Biofuel
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_Anthro
!
! !DESCRIPTION: Subroutine Diagn\_Anthro initializes diagnostics for the
!  anthropogenic emissions (ND07, ND28, ND29, ND32, ND36)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Anthro( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      only : Ind_
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Anthropogenic emissions (ND36, ND32, ND29) 
!  ==> write one single anthropogenic emissions diagnostics per species.
!  ==> it is assumed that anthropogenic emissions are given category 1
!      in the HEMCO configuration file
!  ==> Diagnostics are returned in kg/m2/s.
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use anthro
!  28 Aug 2014 - R. Yantosca - Add IF statements to prevent defining diags
!                              for species that aren't present in a given sim
!  22 Apr 2015 - M. Sulprizio- Now save out hydrocarbons in units kgC/m2/s
!  27 Mar 2017 - M. Sulprizio- Make anthropogenic emissions diagnostics 3D;
!                              Add anthropogenic emissions diagnostics for
!                              remaining NEI2011 species
!EOC
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, I, N, N_CO
    INTEGER            :: id_OCPI, id_OCPO
    CHARACTER(LEN=15)  :: SpcName
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=31)  :: DiagnName_AN
    CHARACTER(LEN=31)  :: DiagnName_AC
    CHARACTER(LEN=31)  :: DiagnName_BF
    CHARACTER(LEN=31)  :: DiagnName_SH
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_ANTHRO (hcoi_gc_diagn_mod.F90)'

    ! CO tracer names
    INTEGER, PARAMETER :: N_ANTH_CO             = 5
    CHARACTER(LEN=6)   :: CO_Tracers(N_ANTH_CO) =               &
         (/ 'CO    ', 'COus  ', 'COeur ', 'COasia' , 'COoth ' /) 

    !=======================================================================
    ! DIAGN_Anthro begins here!
    !=======================================================================

    ! Extension number
    ExtNr = 0

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if we are doing a specialty simulation w/o the anthro species below
    IF ( Input_Opt%ITS_A_C2H6_SIM    ) RETURN
    IF ( Input_Opt%ITS_A_CH4_SIM     ) RETURN
    IF ( Input_Opt%ITS_A_CO2_SIM     ) RETURN
    IF ( Input_Opt%ITS_A_HCN_SIM     ) RETURN
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) RETURN
    IF ( Input_Opt%ITS_A_POPS_SIM    ) RETURN
    IF ( Input_Opt%ITS_A_RnPbBe_SIM  ) RETURN
    IF ( Input_Opt%ITS_A_TAGO3_SIM   ) RETURN

    ! Define advected species ID's
    id_OCPI = Ind_('OCPI','A')
    id_OCPO = Ind_('OCPO','A')

    ! ND36 only: VOC's are only defined for fullchem (not tagged CO)
    IF ( ND36 > 0 .and. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       !----------------------------------------
       ! %%%%% Anthropogenic O3 %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'O3', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_O3'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic ALK4 %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'ALK4', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_ALK4'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kgC/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic ACET %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'ACET', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_ACET'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kgC/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic MEK %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'MEK', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_MEK'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kgC/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic ALD2 %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'ALD2', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_ALD2'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kgC/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic RCHO %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'RCHO', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_RCHO'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic MACR %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'MACR', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_MACR'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic PRPE %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'PRPE', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_PRPE'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kgC/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic C3H8 %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'C3H8', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_C3H8'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kgC/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic CH2O %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'CH2O', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_CH2O'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic C2H6 %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'C2H6', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_C2H6'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kgC/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic SO2 %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'SO2', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_SO2'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic SO4 %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'SO4', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_SO4'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic NH3 %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'NH3', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_NH3'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic BCPI %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'BCPI', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_BCPI'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic OCPI %%%%%
       !----------------------------------------

!sfarinaTMP
       ! Only create diagnostic if OCPI is a defined species
       IF ( id_OCPI > 0 ) THEN

          ! HEMCO species ID
          HcoID = GetHemcoId( 'OCPI', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create diagnostic container
          DiagnName = 'ANTHROPOGENIC_OCPI'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = CATEGORY_ANTHRO,   &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 3,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN

       ENDIF

       !----------------------------------------
       ! %%%%% Anthropogenic BCPO %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'BCPO', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_BCPO'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic OCPO %%%%%
       !----------------------------------------

       ! Only create diagnostic if OCPO is a defined species
       IF ( id_OCPO > 0 ) THEN

          ! HEMCO species ID
          HcoID = GetHemcoId( 'OCPO', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create diagnostic container
          DiagnName = 'ANTHROPOGENIC_OCPO'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = CATEGORY_ANTHRO,   &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 3,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',        &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN

       ENDIF

       !----------------------------------------
       ! %%%%% Anthropogenic NO2 %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'NO2', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_NO2'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic HNO2 %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'HNO2', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_HNO2'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN



       !----------------------------------------
       ! %%%%% Anthropogenic BENZ %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'BENZ', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_BENZ'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kgC/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic TOLU %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'TOLU', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_TOLU'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kgC/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic XYLE %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'XYLE', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_XYLE'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kgC/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       !----------------------------------------
       ! %%%%% Anthropogenic EOH %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'EOH', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_EOH'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kgC/m2/s',        &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

    ENDIF

    !-------------------------------------------
    ! %%%%% Anthropogenic NO %%%%%
    !-------------------------------------------
    IF ( ND36 > 0 .OR. ND32 > 0 ) THEN

       ! NO is only defined for the full-chemistry simulation
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

          ! HEMCO species ID
          HcoID = GetHemcoId( 'NO', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create diagnostic container
          DiagnName = 'ANTHROPOGENIC_NO'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = CATEGORY_ANTHRO,   &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 3,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
    ENDIF

    !-------------------------------------------
    ! %%%%% Anthropogenic CO sectors %%%%%
    !-------------------------------------------
    IF ( ND36 > 0 .OR. ND29 > 0 ) THEN

       ! CO is only defined for the full-chemistry and tagged CO simulations
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or.   &
            Input_Opt%ITS_A_TAGCO_SIM         ) THEN

          ! Loop over tagged CO tracers if necessary
          IF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN
             N_CO = N_ANTH_CO
          ELSE
             N_CO = 1
          ENDIF

          ! Loop over all CO tracers
          DO N = 1, N_CO

             ! Pick the various category names
             SELECT CASE( TRIM( CO_Tracers(N) ) )
                CASE( 'CO'   )
                   HcoId        =  GetHemcoId( 'CO', HcoState, LOC, RC )
                   DiagnName_AN = 'ANTHROPOGENIC_CO'
                   DiagnName_AC = 'AIRCRAFT_CO'
                   DiagnName_BF = 'BIOFUEL_CO'    
                   DiagnName_SH = 'SHIP_CO'
                CASE( 'COus' )                 
                   HcoId        =  GetHemcoId( 'COus', HcoState, LOC, RC )
                   DiagnName_AN = 'ANTHRO_BIOFUEL_TAGCO_US'
                   DiagnName_AC = 'AIRCRAFT_TAGCO_US'
                   DiagnName_BF = ''
                   DiagnName_SH = 'SHIP_TAGCO_US'
                CASE( 'COeur'  )                 
                   HcoId        =  GetHemcoId( 'COeur', HcoState, LOC, RC )
                   DiagnName_AN = 'ANTHRO_BIOFUEL_TAGCO_EUR'
                   DiagnName_AC = 'AIRCRAFT_TAGCO_EUR'
                   DiagnName_BF = ''
                   DiagnName_SH = 'SHIP_TAGCO_EUR'
                CASE( 'COasia' )
                   HcoId        =  GetHemcoId( 'COasia', HcoState, LOC, RC )
                   DiagnName_AN = 'ANTHRO_BIOFUEL_TAGCO_ASIA'
                   DiagnName_AC = 'AIRCRAFT_TAGCO_ASIA'
                   DiagnName_BF = ''
                   DiagnName_SH = 'SHIP_TAGCO_ASIA'
                CASE( 'COoth'  )
                   HcoId        =  GetHemcoId( 'COoth', HcoState, LOC, RC )
                   DiagnName_AN = 'ANTHRO_BIOFUEL_TAGCO_OTHER'
                   DiagnName_AC = 'AIRCRAFT_TAGCO_OTHER'
                   DiagnName_BF = ''
                   DiagnName_SH = 'SHIP_TAGCO_OTHER'
                CASE DEFAULT
                   HcoId        = -1
                   DiagnName_AN = ''
                   DiagnName_AC = ''
                   DiagnName_BF = ''
                   DiagnName_SH = ''
             END SELECT
             
             ! If a valid tracer
             IF ( HcoId > 0 ) THEN
      
                ! Anthropogenic
                CALL Diagn_Create( am_I_Root,                             & 
                                   HcoState  = HcoState,                  &
                                   cName     = TRIM( DiagnName_AN ),      &
                                   ExtNr     = ExtNr,                     &
                                   Cat       = CATEGORY_ANTHRO,           &
                                   Hier      = -1,                        &
                                   HcoID     = HcoID,                     &
                                   SpaceDim  = 3,                         &
                                   LevIDx    = -1,                        &
                                   OutUnit   = 'kg/m2/s',                 &
                                   COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                   AutoFill  = 1,                         &
                                   RC        = RC                        ) 

                ! Aircraft
                CALL Diagn_Create( am_I_Root,                             & 
                                   HcoState  = HcoState,                  &
                                   cName     = TRIM( DiagnName_AC ),      &
                                   ExtNr     = ExtNr,                     &
                                   Cat       = CATEGORY_AIRCRAFT,         &
                                   Hier      = -1,                        &
                                   HcoID     = HcoID,                     &
                                   SpaceDim  = 3,                         &
                                   LevIDx    = -1,                        &
                                   OutUnit   = 'kg/m2/s',                 &
                                   COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                   AutoFill  = 1,                         &
                                   RC        = RC                        ) 

                ! Biofuel
                ! (NOTE: For tagged CO, biofuel is lumped in w/ anthro)
                IF ( LEN_TRIM( DiagnName_BF ) > 0 ) THEN 
                   CALL Diagn_Create( am_I_Root,                          & 
                                      HcoState  = HcoState,               &
                                      cName     = TRIM( DiagnName_BF ),   &
                                      ExtNr     = ExtNr,                  &
                                      Cat       = CATEGORY_BIOFUEL,       &
                                      Hier      = -1,                     &
                                      HcoID     = HcoID,                  &
                                      SpaceDim  = 2,                      &
                                      LevIDx    = -1,                     &
                                      OutUnit   = 'kg/m2/s',              &
                                      COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                      AutoFill  = 1,                      &
                                      RC        = RC                     ) 
                ENDIF

                ! Ship
                CALL Diagn_Create( am_I_Root,                             & 
                                   HcoState  = HcoState,                  &
                                   cName     = TRIM( DiagnName_SH ),      &
                                   ExtNr     = ExtNr,                     &
                                   Cat       = CATEGORY_SHIP,             &
                                   Hier      = -1,                        &
                                   HcoID     = HcoID,                     &
                                   SpaceDim  = 2,                         &
                                   LevIDx    = -1,                        &
                                   OutUnit   = 'kg/m2/s',                 &
                                   COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                   AutoFill  = 1,                         &
                                   RC        = RC                        ) 

             ENDIF
          ENDDO
       ENDIF
    ENDIF
#endif

  END SUBROUTINE Diagn_Anthro
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_Biogenic
!
! !DESCRIPTION: Subroutine Diagn\_Biogenic initializes diagnostics for the
!  mineral dust aerosols (ND06).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Biogenic( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCO_State_Mod,      ONLY : HCO_GetHcoId
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Biogenic emissions (ND46, ND29, ND11) 
!  ==> Biogenic emissions are taken from MEGAN inventory
!  ==> write one single biogenic emissions diagnostics per species.
!  ==> Diagnostics are returned in kg/m2/s.
!  ==> Oceanic acetone is defined in ND11.
!  ==> For now, only GC species totals are diagnosed. To add 
!      individual MEGAN species (Sabinene, Limonene, ...), those
!      have to be explicitly added to hcox_megan_mod (as for acetone) 
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use MEGAN
!  18 Feb 2015 - M. Sulprizio- Add manual diagnostics for individual MEGAN
!                              species (MBOX, APIN, BPIN, etc.)
!  10 Mar 2015 - R. Yantosca - Remove double-definition of BIOGENIC_LIMO
!  30 Mar 2015 - R. Yantosca - Bug fix: Now test if Br2 is a HEMCO species
!  22 Apr 2015 - M. Sulprizio- Now save out hydrocarbons in units kgC/m2/s
!  02 Jun 2016 - R. Yantosca - Bug fix: only save seasalt Br2 diagnostics
!                              for full-chemistry or aerosol-only simulations
!  01 Mar 2017 - M. Sulprizio- Add ALD2 senescing, EOH senescing, and ALD2 from
!                              ocean source
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: YesOrNo
    INTEGER            :: Cat, ExtNr, HcoID, I, N
    CHARACTER(LEN=1)   :: ISTR
    CHARACTER(LEN=15)  :: SpcName
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_BIOGENIC (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_BIOGENIC begins here!
    !=======================================================================
    
    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if we are doing a specialty simulation w/o biofuels
    IF ( Input_Opt%ITS_A_HCN_SIM     ) RETURN
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) RETURN
    IF ( Input_Opt%ITS_A_POPS_SIM    ) RETURN
    IF ( Input_Opt%ITS_A_RnPbBe_SIM  ) RETURN
    IF ( Input_Opt%ITS_A_TAGO3_SIM   ) RETURN

    ! Extension and category #'s for MEGAN
    ExtNr = GetExtNr( HcoState%Config%ExtList, 'MEGAN')
    Cat   = -1

    ! Make sure MEGAN is on if ND46 is used
    IF ( ExtNr <= 0 .AND. ND46 > 0 ) THEN
       MSG = 'MEGAN is not enabled - cannot write biogenic diagnostics!'
       CALL HCO_Error ( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
 
    ! Only if MEGAN is on ... 
    IF ( ExtNr > 0 ) THEN


#if ! defined ( TOMAS )
       IF ( ND07 > 0 ) THEN
#endif
         !-------------------------------------------------------------------
         ! %%%%% diag for direct emission of SOAS OC for TOMAS          %%%%%
         ! %%%%% this is optional for non-tomas simulations             %%%%%
         ! %%%%% this is not optional for tomas simulations             %%%%%
         !-------------------------------------------------------------------
          ! HEMCO species ID
          HcoID = HCO_GetHcoID( 'SOAS', HcoState )
      
          ! Create diagnostic container
          IF ( HcoID > 0 ) THEN
             DiagnName = 'BIOGENIC_SOAS'
             CALL Diagn_Create( am_I_Root,                     &  
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &  
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 
          ENDIF

#if ! defined ( TOMAS )
       ENDIF
#endif

       ! ND46 only
       IF ( ND46 > 0 ) THEN

          !----------------------------------------
          ! %%%%% Biogenic ISOP %%%%%
          !----------------------------------------

          ! HEMCO species ID
          HcoID = HCO_GetHcoID( 'ISOP', HcoState )
          
          ! Create diagnostic container
          IF ( HcoID > 0 ) THEN
             DiagnName = 'BIOGENIC_ISOP'
             CALL Diagn_Create ( am_I_Root,                     & 
                                 HcoState  = HcoState,          &
                                 cName     = TRIM( DiagnName ), &
                                 ExtNr     = ExtNr,             &
                                 Cat       = Cat,               &
                                 Hier      = -1,                &
                                 HcoID     = HcoID,             &
                                 SpaceDim  = 2,                 &
                                 LevIDx    = -1,                &
                                 OutUnit   = 'kgC/m2/s',        &
                                 COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                 AutoFill  = 1,                 &
                                 RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          !----------------------------------------
          ! %%%%% Biogenic ALD2 %%%%%
          !----------------------------------------

          ! HEMCO species ID
          HcoID = HCO_GetHcoID( 'ALD2', HcoState )

          ! Create diagnostic container (if ALD2 is defined)
          IF ( HcoID > 0 ) THEN
             DiagnName = 'BIOGENIC_ALD2'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kgC/m2/s',        &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          !----------------------------------------
          ! %%%%% Biogenic EOH %%%%%
          !----------------------------------------

          ! HEMCO species ID
          HcoID = HCO_GetHcoID( 'EOH', HcoState )

          ! Create diagnostic container (if EOH is defined)
          IF ( HcoID > 0 ) THEN
             DiagnName = 'BIOGENIC_EOH'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kgC/m2/s',        &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          !----------------------------------------
          ! %%%%% Biogenic PRPE %%%%%
          !----------------------------------------

          ! HEMCO species ID
          HcoID = HCO_GetHcoID( 'PRPE', HcoState )

          ! Create diagnostic container
          IF ( HcoID > 0 ) THEN
             DiagnName = 'BIOGENIC_PRPE'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kgC/m2/s',        &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          !----------------------------------------
          ! %%%%% Biogenic C2H4 %%%%%
          !----------------------------------------

          ! HEMCO species ID
          HcoID = HCO_GetHcoID( 'C2H4', HcoState )

          ! Create diagnostic container (if C2H4 is defined)
          IF ( HcoID > 0 ) THEN
             DiagnName = 'BIOGENIC_C2H4'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kgC/m2/s',        &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          !----------------------------------------
          ! %%%%% Biogenic CHBr3 %%%%%
          !----------------------------------------

          ! HEMCO species ID
          HcoID = HCO_GetHcoID( 'CHBr3', HcoState )

          ! Create diagnostic container
          ! NOTE: CHBr3 and CH2Br2 are emitted through 
          ! HEMCO core, i.e. extension number is 0!
          IF ( HcoID > 0 ) THEN
             DiagnName = 'BIOGENIC_CHBR3'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = 0,                 &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          !----------------------------------------
          ! %%%%% Biogenic CH2Br2 %%%%%
          !----------------------------------------

          ! HEMCO species ID
          HcoID = HCO_GetHcoId( 'CH2Br2', HcoState )

          ! Create diagnostic container
          ! NOTE: CHBr3 and CH2Br2 are emitted through 
          ! HEMCO core, i.e. extension number is 0!
          IF ( HcoID > 0 ) THEN
             DiagnName = 'BIOGENIC_CH2BR2'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = 0,                 &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%% These diagnostics are explicitly filled in hcox_megan_mod.F 
          !%%% The diagnostics  name defined below must match the names 
          !%%% used in the MEGAN extension!
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          ! There are 3 manual diagnostics in MEGAN
          DO I = 1,3
   
             ! Define diagnostics names. These names have to match the
             ! names called in hcox_megan_mod.F90.
             IF ( I == 1 ) THEN
                DiagnName = 'BIOGENIC_FAXX'
             ELSEIF ( I == 2 ) THEN
                DiagnName = 'BIOGENIC_AAXX'
             ELSEIF ( I == 3 ) THEN
                DiagnName = 'BIOGENIC_MOHX'
             ENDIF
      
             ! Create diagnostics. Don't use AutoFill here since the 
             ! diagnostics update calls are explicitly called in 
             ! hcox_megan_mod.F90.
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = -1,                &
                                Hier      = -1,                &
                                HcoID     = -1,                &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                OutOper   = 'Mean',            &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 
          ENDDO

       ENDIF

       !----------------------------------------
       ! %%%%% Biogenic ACET %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'ACET', HcoState )

       IF ( HcoID > 0 ) THEN
          !%%% For ND46 or ND11 %%%
          IF ( ND46 > 0 .OR. ND11 > 0 ) THEN
   
             ! Create diagnostic container
             DiagnName = 'BIOGENIC_ACET'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kgC/m2/s',        &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF
   
          !%%% For ND11 only %%%
          !%%% These diagnostics are explicitly filled in hcox_megan_mod. 
          !%%% The diagnostics  name defined below must match the names 
          !%%% used in the MEGAN extension!
          IF ( ND11 > 0 ) THEN
   
             ! There are three manual acetone diagnostics in MEGAN
             DO I = 1,3
   
                ! Define diagnostics names. These names have to match the
                ! names called in hcox_megan_mod.F90.
                IF ( I == 1 ) THEN
                   DiagnName = 'MEGAN_ACET_MONO'
                ELSEIF ( I == 2 ) THEN
                   DiagnName = 'MEGAN_ACET_MBO'
                ELSEIF ( I == 3 ) THEN
                   DiagnName = 'MEGAN_ACET_DIRECT'
                ENDIF
      
                ! Create diagnostics. Don't use AutoFill here since the 
                ! diagnostics update calls are explicitly called in 
                ! hcox_megan_mod.F90.
                CALL Diagn_Create( am_I_Root,                     & 
                                   HcoState  = HcoState,          &
                                   cName     = TRIM( DiagnName ), &
                                   ExtNr     = ExtNr,             &
                                   Cat       = -1,                &
                                   Hier      = -1,                &
                                   HcoID     = HcoID,             &
                                   SpaceDim  = 2,                 &
                                   LevIDx    = -1,                &
                                   OutUnit   = 'kgC/m2/s',        &
                                   COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                   AutoFill  = 0,                 &
                                   RC        = RC                  ) 
                IF ( RC /= HCO_SUCCESS ) RETURN 
             ENDDO
          ENDIF
       ENDIF ! ACET
    ENDIF !MEGAN

    !=======================================================================
    ! These diagnostics use the MEGAN Monoterpenes extension
    !=======================================================================


    ! Extension # of MEGAN monoterpenes
    ExtNr = GetExtNr( HcoState%Config%ExtList, 'MEGAN_Mono')
    IF ( ExtNr > 0 ) THEN

       !%%% For ND46 diagnostic %%%
       IF ( ND46 > 0 ) THEN

          !----------------------------------------
          ! %%%%% Biogenic monoterpene species %%%%%
          !----------------------------------------

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%% These diagnostics are explicitly filled in hcox_megan_mod.F 
          !%%% The diagnostics  name defined below must match the names 
          !%%% used in the MEGAN extension!
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          ! There are 9 manual monoterpene diagnostics in MEGAN
          DO I = 1,13
   
             ! Define diagnostics names. These names have to match the
             ! names called in hcox_megan_mod.F90.
             IF ( I == 1 ) THEN
                DiagnName = 'BIOGENIC_MBOX'
             ELSEIF ( I ==  2 ) THEN
                DiagnName = 'BIOGENIC_APIN'
             ELSEIF ( I ==  3 ) THEN
                DiagnName = 'BIOGENIC_BPIN'
             ELSEIF ( I ==  4 ) THEN
                DiagnName = 'BIOGENIC_LIMO'
             ELSEIF ( I ==  5 ) THEN
                DiagnName = 'BIOGENIC_SABI'
             ELSEIF ( I ==  6 ) THEN
                DiagnName = 'BIOGENIC_MYRC'
             ELSEIF ( I ==  7 ) THEN
                DiagnName = 'BIOGENIC_CARE'
             ELSEIF ( I ==  8 ) THEN
                DiagnName = 'BIOGENIC_OCIM'
             ELSEIF ( I ==  9 ) THEN
                DiagnName = 'BIOGENIC_OMON'
             ! Now include diagnostic for total monoterpenes (MONX)
             ! even if tracer is not defined (mps, 2/23/15)
             ELSEIF ( I == 10 ) THEN
                DiagnName = 'BIOGENIC_MONX'
             ! Make diagnostic containers for SOA species to avoid errors
             ! in diag3.F. These diagnostics will be zero if MEGAN_SOA is
             ! turned off. (mps, 2/23/15)
             ELSEIF ( I == 11 ) THEN
                DiagnName = 'BIOGENIC_FARN'
             ELSEIF ( I == 12 ) THEN
                DiagnName = 'BIOGENIC_BCAR'
             ELSEIF ( I == 13 ) THEN
                DiagnName = 'BIOGENIC_OSQT'
             ENDIF

             ! Create diagnostics. Don't use AutoFill here since the 
             ! diagnostics update calls are explicitly called in 
             ! hcox_megan_mod.F90.
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = -1,                &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                OutOper   = 'Mean',            &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 
          ENDDO

       ENDIF

       !%%% For ND07 diagnostic %%%
       IF ( ND07 > 0 ) THEN

          !----------------------------------------
          ! %%%%% Biogenic OC %%%%%
          !----------------------------------------

          ! HEMCO species ID
          HcoID = HCO_GetHcoID( 'OCPI', HcoState )

          ! Create diagnostic container
          IF ( HcoID > 0 ) THEN
             DiagnName = 'BIOGENIC_OCPI'
             CALL Diagn_Create( am_I_Root,                     &  
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &  
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF
       ENDIF

       !%%% For ND29 diag %%%
       IF ( ND29 > 0 ) THEN

          !----------------------------------------
          ! %%%%% Biogenic CO %%%%%
          !----------------------------------------

          ! HEMCO species ID
          HcoID = GetHemcoId( 'CO', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create diagnostic container
          DiagnName = 'BIOGENIC_CO'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF
    ENDIF ! Megan mono

    !=======================================================================
    ! These diagnostics use the SeaSalt extension
    !=======================================================================
    IF ( ( ND46 > 0                     )   .and. &
         ( Input_Opt%ITS_A_FULLCHEM_SIM ) ) THEN

       ! Extension # of SeaSalt
       ExtNr = GetExtNr( HcoState%Config%ExtList, 'SeaSalt')
       IF ( ExtNr <= 0 ) THEN
          CALL HCO_Error ( 'SeaSalt extension not enabled', RC, THISLOC=LOC )
          RETURN      
       ENDIF

       ! Find out if SeaSalt Br2 is enabled
       CALL GetExtOpt ( HcoState%Config, ExtNr, &
                       'Emit Br2', OptValBool=YesOrNo, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Only save out SeaSalt Br2 diagnostic if the Br2 option is enabled
       IF ( YesOrNo ) THEN 

          !----------------------------------------
          ! %%%%% Biogenic Br2 %%%%%
          !----------------------------------------

          ! HEMCO species ID
          HcoID = GetHemcoId( 'Br2', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create diagnostic container
          IF ( HcoID > 0 ) THEN
             CALL Diagn_Create( am_I_Root,                    & 
                                HcoState  = HcoState,         &
                                cName     = 'SEASALT_BR2',    &
                                ExtNr     = ExtNr,            &
                                Cat       = -1,               &
                                Hier      = -1,               &
                                HcoID     = HcoID,            &
                                SpaceDim  = 2,                &
                                LevIDx    = -1,               &
                                OutUnit   = 'kg/m2/s',        &
                                COL       = HcoState%Diagn%HcoDiagnIDManual, &
                                AutoFill  = 1,                &
                                RC        = RC                 ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 
          ENDIF
       ENDIF
    ENDIF

    !=======================================================================
    ! These diagnostics are from the base emissions
    !=======================================================================
    IF ( ( ND46 > 0                     )   .and. &
         ( Input_Opt%ITS_A_FULLCHEM_SIM ) ) THEN

       !----------------------------------------
       ! %%%%% ALD2 from senescing plants %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'ALD2', HcoState )

       ! Create diagnostic container
       IF ( HcoID > 0 ) THEN
          ExtNr     = 0
          DiagnName = 'ALD2_SENESCING'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = CATEGORY_BIOGENIC, &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',   &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF

       !----------------------------------------
       ! %%%%% EOH from senescing plants %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'EOH', HcoState )

       ! Create diagnostic container
       IF ( HcoID > 0 ) THEN
          ExtNr     = 0
          DiagnName = 'EOH_SENESCING'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = CATEGORY_BIOGENIC, &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',   &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF

    ENDIF

    !=======================================================================
    ! These diagnostics use the SeaFlux extension
    !
    ! As we did for acetone, only keep track of flux from ocean. Deposition
    ! from the atmosphere is handled by drydep.
    !=======================================================================
    IF ( ( ND46 > 0                     )   .and. &
         ( Input_Opt%ITS_A_FULLCHEM_SIM     .or.  &
           Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN

       !----------------------------------------
       ! %%%%% ALD2 from ocean source %%%%%
       !----------------------------------------
       ExtNr = GetExtNr( HcoState%Config%ExtList, 'SeaFlux')
       IF ( ExtNr <= 0 ) THEN
          CALL HCO_Error ( 'SeaFlux extension not enabled', RC, THISLOC=LOC )
          RETURN      
       ENDIF

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'ALD2', HcoState )

       ! Create diagnostic container
       IF ( HcoID > 0 ) THEN
          DiagnName = 'ALD2_OCEAN_SOURCE'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',   &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF

    ENDIF
#endif

  END SUBROUTINE Diagn_Biogenic
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_LFlash
!
! !DESCRIPTION: Subroutine Diagn\_LFlash initializes diagnostics for the
!  lightning flash diagnostics (ND56).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_LFlash( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Lightning flashes (ND56)
!  ==> Manually set in hcox_lightning_mod.F90
!  ==> The original code diagnoses the 'cumulative mean', i.e. it
!      sums all quantities and simply divides by the number of 
!      met. time steps when writing the diagnostics. Here, we only 
!      use the average since the last writeout. To simulate the
!      behavior of the original ND56 diagnostics, set 'OutOper' to
!      'Cumsum' (cumulative sum) and manually divide by the # of
!      timesteps when writing data to the output file!
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use lightning=
!  19 Feb 2015 - C. Keller   - Added entry for MTYPE
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, I, N, COL
    CHARACTER(LEN=1)   :: ISTR
    CHARACTER(LEN=15)  :: SpcName
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_LFLASH (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_LFLASH begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if we are doing a specialty simulation w/o lightning
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

    ! Define diagnostics
    IF ( Input_Opt%ND56 > 0 ) THEN

       ! Extension number
       ExtNr = GetExtNr( HcoState%Config%ExtList, 'LightNOx')

       ! Exit if LightNOx is not turned on - the lightning NOx extension
       ! must be enabled in the HEMCO configuration file.
       IF ( ExtNr <= 0 ) THEN
          MSG = 'Lightning NOx is not enabled - cannot write diagnostics ND56!'
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! HEMCO species ID
       HcoID = GetHemcoId( 'NO', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Define collection: in development mode or if netCDF is enabled,
       ! add it to the default HEMCO collection. Otherwise, add it to the
       ! manual collection and the diagnostics will be written to the
       ! bpch file in diag3.F.
#if defined( NC_DIAG )
       COL = HcoState%Diagn%HcoDiagnIDDefault
#else
       COL = HcoState%Diagn%HcoDiagnIDManual
#endif

       ! Loop over lighthing flash quantities
       DO I = 1, 3

          ! Pick the proper diagnostic name
          SELECT CASE( I )
             CASE( 1 )
                DiagnName = 'LIGHTNING_TOTAL_FLASHRATE'
             CASE( 2 )
                DiagnName = 'LIGHTNING_INTRACLOUD_FLASHRATE'
             CASE( 3 )
                DiagnName = 'LIGHTNING_CLOUDGROUND_FLASHRATE'
          END SELECT

          ! Define diagnostics ID
          N = 56000 + I

          ! Create diagnostic container
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             cID       = N,                 &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'flashes/min/km2', &
                             OutOper   = 'Mean',            &
                             COL       = COL,               &
                             AutoFill  = 0,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDDO
 
       ! ---------------------------------------------------------- 
       ! Diagnostics for convective cloud top height.
       ! ---------------------------------------------------------- 

       ! Define diagnostics name and ID
       DiagnName = 'LIGHTNING_CLOUD_TOP'
       N         = 56004 

       ! Create diagnostic container
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          cID       = N,                 & 
                          ExtNr     = ExtNr,             &
                          Cat       = -1,                &
                          Hier      = -1,                &
                          HcoID     = -1,                &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = '1',               &
                          OutOper   = 'Mean',            &
                          COL       = COL,               &
                          AutoFill  = 0,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

    ENDIF ! ND56 
#endif

  END SUBROUTINE Diagn_LFlash
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_ParaNOx
!
! !DESCRIPTION: Subroutine Diagn\_ParaNox initializes diagnostics for the
!  ParaNOx ship plume emissions  (ND32, ND63).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_ParaNOx( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  SHIP DIAGNOSTICS (ND63, ND32)
!   ==> Manually set in hcox_paranox_mod.F90
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use PARANOX
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, N
    CHARACTER(LEN=15)  :: SpcName
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_PARANOX (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_PARANOX begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if we are doing a specialty simulation w/o lightning
    IF ( .NOT. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

    ! Extension number
    ExtNr = GetExtNr( HcoState%Config%ExtList, 'ParaNOx')

    ! Exit if PARANOX extension was turned off
    IF ( ExtNr <= 0 .AND. Input_Opt%DO_ND63 ) THEN
       MSG = 'ParaNOx is not enabled - cannot write diagnostics ND63!'
       CALL HCO_Error( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    
    IF ( ExtNr > 0 ) THEN

       ! HEMCO species ID
       HcoID = GetHemcoId( 'NO', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       !-------------------------------------------
       ! %%%%% Ship NO %%%%%
       !
       ! These are the final NO emissions
       !-------------------------------------------       
       IF ( Input_Opt%DO_ND63 .OR. ND32 > 0 ) THEN
          
          ! Create diagnostic container
          DiagnName = 'SHIP_NO'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
 
       !-------------------------------------------
       ! %%%%% PARANOX specific diagnostics %%%%%
       !-------------------------------------------
       IF ( Input_Opt%DO_ND63 ) THEN
 
          ! These are the ship NO emissions before PARANOX chemistry
          DiagnName = 'PARANOX_TOTAL_SHIPNOX'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 0,                 &
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN
  
          DiagnName = 'PARANOX_NOXFRAC_REMAINING'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = '1',               &
                             OutOper   = 'Mean',            &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 0,                 &
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN
  
          DiagnName = 'PARANOX_OPE'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = '1',               &
                             OutOper   = 'Mean',            &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 0,                 &
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! This is the O3 production through PARANOX 
          HcoID = GetHemcoId( 'O3', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( HcoID > 0 ) THEN
             DiagnName = 'PARANOX_O3_PRODUCTION'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = -1,                &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 0,                 &
                                RC        = RC                  )
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF 
       ENDIF
    ENDIF
#endif

  END SUBROUTINE Diagn_ParaNOx
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_POPs
!
! !DESCRIPTION: Subroutine Diagn\_POPs initializes several HEMCO manual 
!  diagnostics for the POPs simulation.  These diagnostics are updated
!  in the HEMCO extensions module hcox\_gc\_POPs\_mod.F90.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_POPs( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Split off code from HCOI_GC_Diagn_Init into smaller routines in order to
!  make the code more manageable.
!
! !REVISION HISTORY: 
!  26 Aug 2014 - M. Sulprizio- Initial version
!  15 Oct 2018 - R. Yantosca - Changed "AD53_*" to "GCPOPS_*"
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, N, I
    CHARACTER(LEN=31)  :: DiagnName, OutOper, OutUnit
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_POPs (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! Define ND53 diagnostics (POPs emissions)
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Exit if the POPs simulation is not selected
    IF ( .not. Input_Opt%ITS_A_POPS_SIM ) RETURN

#if defined( NC_DIAG ) 
    ! For the HISTORY netCDF diagnostics, we want to get the instantaneous
    ! values archived by HEMCO and then let HISTORY do the averaging.
    OutOper = 'Instantaneous'
#endif

#if defined( BPCH_DIAG )
    ! Exit if ND53 diagnostics aren't turned on
    IF ( ND53 <= 0 ) RETURN

    ! For the bpch diagnostics, change units to kg/s to help in validating
    ! the netCDF diagnostics.  But select "Mean" diagnostics since these
    ! are only ouptut at the end.
    OutOper = 'Mean'
#endif

    ! Define diagnostics
    IF ( ExtState%GC_POPs > 0 ) THEN

       ! HEMCO extension # for POPs
       ExtNr = GetExtNr( HcoState%Config%ExtList, 'GC_POPs' )
       IF ( ExtNr <= 0 ) THEN
          MSG = 'Cannot find the POPs extension for HEMCO!'
          CALL HCO_Error( MSG, RC, THISLOC=LOC )
          RETURN      
       ENDIF

       !-------------------------------------------
       ! %%%%% Gas-phase POP emissions %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'POPG', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'GCPOPS_POPG_SOURCE'
       CALL Diagn_Create( am_I_Root,                                         & 
                          HcoState  = HcoState,                              &
                          cName     = TRIM( DiagnName ),                     &
                          ExtNr     = ExtNr,                                 &
                          Cat       = -1,                                    &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &
                          OutOper   = OutOper,                               &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                    ) 

       ! Trap potential errors
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Could not define POPs diagnostic: '// TRIM( DiagnName )
          CALL HCO_Error( MSG, RC, THISLOC=LOC )
          RETURN 
       ENDIF

       !-------------------------------------------
       ! %%%%% OC-phase POP emissions%%%%%
       !-------------------------------------------
 
       ! HEMCO species ID
       HcoID = GetHemcoId( 'POPPOCPO', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'GCPOPS_POPPOCPO_SOURCE'
       CALL Diagn_Create( am_I_Root,                                         & 
                          HcoState  = HcoState,                              &
                          cName     = TRIM( DiagnName ),                     &
                          ExtNr     = ExtNr,                                 &
                          Cat       = -1,                                    &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &  
                          OutOper   = OutOper,                               &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                    ) 

       ! Trap potential errors
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Could not define POPs diagnostic: '// TRIM( DiagnName )
          CALL HCO_Error( MSG, RC, THISLOC=LOC )
          RETURN 
       ENDIF

       !-------------------------------------------
       ! %%%%% BC-phase POP emissions %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'POPPBCPO', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'GCPOPS_POPPBCPO_SOURCE'
       CALL Diagn_Create( am_I_Root,                                         &
                          HcoState  = HcoState,                              &
                          cName     = TRIM(DiagnName),                       &
                          ExtNr     = ExtNr,                                 &
                          Cat       = -1,                                    &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &
                          OutOper   = OutOper,                               &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                    ) 

       ! Trap potential errors
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Could not define POPs diagnostic '// TRIM( DiagnName )
          CALL HCO_Error( MSG, RC, THISLOC=LOC )
          RETURN 
       ENDIF

       !-------------------------------------------
       ! %%%%% Manual diagnostics %%%%%
       !-------------------------------------------

       DO I = 1,12

          ! Define diagnostic names. These have to match the names
          ! in module HEMCO/Extensions/hcox_gc_POPs_mod.F90.
          IF ( I == 1 ) THEN
             DiagnName = 'GCPOPS_POPG_SOIL'
             OutUnit   = 'kg/m2/s'
          ELSEIF ( I == 2  ) THEN
             DiagnName = 'GCPOPS_POPG_LAKE'
             OutUnit   = 'kg/m2/s'
          ELSEIF ( I == 3  ) THEN
             DiagnName = 'GCPOPS_POPG_LEAF'
             OutUnit   = 'kg/m2/s'
          ELSEIF ( I == 4  ) THEN
             DiagnName = 'GCPOPS_SOIL2AIR'
             OutUnit   = 'ng/m2/day'
          ELSEIF ( I == 5  ) THEN
             DiagnName = 'GCPOPS_AIR2SOIL'
             OutUnit   = 'ng/m2/day'
          ELSEIF ( I == 6  ) THEN
             DiagnName = 'GCPOPS_LAKE2AIR'
             OutUnit   = 'ng/m2/day'
          ELSEIF ( I == 7  ) THEN
             DiagnName = 'GCPOPS_AIR2LAKE'
             OutUnit   = 'ng/m2/day'
          ELSEIF ( I == 8  ) THEN
             DiagnName = 'GCPOPS_LEAF2AIR'
             OutUnit   = 'ng/m2/day'
          ELSEIF ( I == 9  ) THEN
             DiagnName = 'GCPOPS_AIR2LEAF'
             OutUnit   = 'ng/m2/day'
          ELSEIF ( I == 10 ) THEN
             DiagnName = 'GCPOPS_SOILAIR_FUG'
             OutUnit   = '1'
          ELSEIF ( I == 11 ) THEN
             DiagnName = 'GCPOPS_LAKEAIR_FUG'
             OutUnit   = '1'
          ELSEIF ( I == 12 ) THEN
             DiagnName = 'GCPOPS_LEAFAIR_FUG'
             OutUnit   = '1'
          ENDIF

          ! Create manual diagnostics
          CALL Diagn_Create( am_I_Root,                                      & 
                             HcoState  = HcoState,                           &
                             cName     = TRIM( DiagnName ),                  &
                             ExtNr     = ExtNr,                              &
                             Cat       = -1,                                 &
                             Hier      = -1,                                 &
                             HcoID     = -1,                                 &
                             SpaceDim  = 2,                                  &
                             LevIDx    = -1,                                 &
                             OutUnit   = OutUnit,                            &
                             OutOper   = OutOper,                            &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,    &
                             AutoFill  = 1,                                  &
                             RC        = RC                                 ) 

          ! Trap potential errors
          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Could not define POPs diagnostic: '// TRIM( DiagnName )
             CALL HCO_Error( Msg, RC, THISLOC=LOC )
             RETURN 
          ENDIF
       ENDDO

    ENDIF

  END SUBROUTINE Diagn_POPs
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_CH4
!
! !DESCRIPTION: Subroutine Diagn\_CH4 initializes diagnostics for the
!  CH4 simulation (ND58).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_CH4( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCO_State_Mod,      ONLY : HCO_GetHcoID
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : Ind_
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Split off code from HCOI_GC_Diagn_Init into smaller routines in order to
!  make the code more manageable.
!
!  Biomass diagnostics are defined in routine Diagn_Biomass.
!\\
!\\
!  CH4 diagnostics need to be defined even if ND58 is turned off because
!  the diagnostics are also being used to write CH4 emissions from the 
!  individual sources (gas, coal, etc.) into STT (in global\_ch4\_mod.F).
!  The categories defined here need to match the ones specified in the 
!  HEMCO configuration file.
!
! !REVISION HISTORY: 
!  13 Sep 2014 - C. Keller   - Initial version
!  16 Jun 2016 - C. Miller   - Now define species ID's with Ind_ funciton
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, id_CH4, Cat, HcoID, N
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_CH4 (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! Define ND58 diagnostics (CH4 emissions)
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Exit if the CH4 simulation is not selected
    IF ( .NOT. ( Input_Opt%ITS_A_CH4_SIM .OR. Ind_('CH4','A') > 0 ) ) RETURN

    ! Get default HEMCO species ID for CH4 
    id_CH4 = HCO_GetHcoID( 'CH4', HcoState )

    ! Extension number is zero (HEMCO core) until defined otherwise
    ExtNr = 0

    !-----------------------------------------------------------------
    ! %%%%% CH4 from oil (Category 1 or species CH4_OIL)  %%%%%
    !-----------------------------------------------------------------

    ! Check if tagged CH4 simulation
    ! Otherwise, use CH4 category 1 emissions
    Cat   = 1
    HcoID = HCO_GetHcoID( 'CH4_OIL', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_OIL'
       CALL Diagn_Create( am_I_Root,                                         & 
                          HcoState  = HcoState,                              &
                          cName     = TRIM( DiagnName ),                     &
                          ExtNr     = ExtNr,                                 &
                          Cat       = Cat,                                   &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                    ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------
    ! %%%%% CH4 from natural gas (Category 2 or species CH4_GAS)  %%%%%
    !-----------------------------------------------------------------

    ! Check if tagged CH4 simulation
    ! Otherwise, use CH4 category 1 emissions
    Cat   = 2
    HcoID = HCO_GetHcoID( 'CH4_GAS', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_GAS'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !----------------------------------------------------------
    ! %%%%% CH4 from coal (Category 3 or species CH4_COL)  %%%%%
    !----------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 3
    HcoID = HCO_GetHcoID( 'CH4_COL', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_COAL'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !---------------------------------------------------------------
    ! %%%%% CH4 from livestock (Category 4 or species CH4_LIV)  %%%%%
    !---------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 4
    HcoID = HCO_GetHcoID( 'CH4_LIV', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_LIVESTOCK'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF   
 
    !---------------------------------------------------------------
    ! %%%%% CH4 from landfills (Category 5 or species CH4_LDF)  %%%%%
    !---------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 5
    HcoID = HCO_GetHcoID( 'CH4_LDF', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_LANDFILLS'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF  

    !---------------------------------------------------------------
    ! %%%%% CH4 from wastewater (Category 6 or species CH4_WST)  %%%%%
    !---------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 6
    HcoID = HCO_GetHcoID( 'CH4_WST', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_WASTEWATER'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF   
 
    !---------------------------------------------------------------
    ! %%%%% CH4 from rice (Category 7 or species CH4_RIC)  %%%%%
    !---------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 7
    HcoID = HCO_GetHcoID( 'CH4_RIC', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_RICE'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from other anth. sources (Category 8 or species CH4_OTA)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 8
    HcoID = HCO_GetHcoID( 'CH4_OTA', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_ANTHROTHER'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from biomass burning (Category 9 or species CH4_BBN)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 9
    HcoID = HCO_GetHcoID( 'CH4_BBN', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_BIOMASS'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from wetlands (Category 10 or species CH4_WTL)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 10
    HcoID = HCO_GetHcoID( 'CH4_WTL', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_WETLAND'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from seeps (Category 11 or species CH4_SEE)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 11
    HcoID = HCO_GetHcoID( 'CH4_SEE', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_SEEPS'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from lakes (Category 12 or species CH4_LAK)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 12
    HcoID = HCO_GetHcoID( 'CH4_LAK', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_LAKES'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from termites (Category 13 or species CH4_TER)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 13
    HcoID = HCO_GetHcoID( 'CH4_TER', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_TERMITES'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from soil absorption (Category 14 or species CH4_SAB)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 14
    HcoID = HCO_GetHcoID( 'CH4_SAB', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_SOILABSORB'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

  END SUBROUTINE Diagn_CH4
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_Hg
!
! !DESCRIPTION: Subroutine Diagn\_Hg initializes diagnostics for the
!  Hg specialty simulation. For now, this is just a placeholder. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Hg( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCO_State_Mod,      ONLY : HCO_GetHcoID
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  23 Sep 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, Cat, N
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_Hg (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! Define diagnostics (POPs emissions)
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Exit if the mercury simulation is not selected
    IF ( .not. Input_Opt%ITS_A_MERCURY_SIM ) RETURN

    ! HEMCO species ID
    HcoID = GetHemcoId( 'Hg0', HcoState, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Extension number is always zero (HEMCO core)
    ExtNr = 0

    !-------------------------------------------
    ! %%%%% ARTISANAL HG (Hg0) %%%%%
    !-------------------------------------------
 
    ! Create diagnostic container
    DiagnName = 'HG0_ARTISANAL'
    Cat       = 8
    CALL Diagn_Create( am_I_Root,                                            & 
                       HcoState  = HcoState,                                 &
                       cName     = TRIM(DiagnName),                          &
                       ExtNr     = ExtNr,                                    &
                       Cat       = Cat,                                      &
                       Hier      = -1,                                       &
                       HcoID     = HcoID,                                    &
                       SpaceDim  = 2,                                        &
                       LevIDx    = -1,                                       &
                       OutUnit   = 'kg/m2/s',                                &
                       COL       = HcoState%Diagn%HcoDiagnIDManual,          &
                       AutoFill  = 1,                                        &
                       RC        = RC                                       ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-------------------------------------------
    ! %%%%% NATURAL HG (Hg0) %%%%%
    !-------------------------------------------
 
    ! Create diagnostic container
    DiagnName = 'HG0_NATURAL'
    Cat       = CATEGORY_NATURAL
    CALL Diagn_Create( am_I_Root,                                            & 
                       HcoState  = HcoState,                                 &
                       cName     = TRIM(DiagnName),                          &
                       ExtNr     = ExtNr,                                    &
                       Cat       = Cat,                                      &
                       Hier      = -1,                                       &
                       HcoID     = HcoID,                                    &
                       SpaceDim  = 2,                                        &
                       LevIDx    = -1,                                       &
                       OutUnit   = 'kg/m2/s',                                &
                       COL       = HcoState%Diagn%HcoDiagnIDManual,          &
                       AutoFill  = 1,                                        &
                       RC        = RC                                       )
    IF ( RC /= HCO_SUCCESS ) RETURN


    !---------------------------------------------
    ! %%%%% ANTHROPOGENIC HG (Hg0, Hg2, HgP) %%%%%
    !---------------------------------------------
 
    Cat = CATEGORY_ANTHRO

    ! Create diagnostic container
    DiagnName = 'HG0_ANTHRO'
    CALL Diagn_Create( am_I_Root,                                            & 
                       HcoState  = HcoState,                                 &
                       cName     = TRIM(DiagnName),                          &
                       ExtNr     = ExtNr,                                    &
                       Cat       = Cat,                                      &
                       Hier      = -1,                                       &
                       HcoID     = HcoID,                                    &
                       SpaceDim  = 2,                                        &
                       LevIDx    = -1,                                       &
                       OutUnit   = 'kg/m2/s',                                &
                       COL       = HcoState%Diagn%HcoDiagnIDManual,          &
                       AutoFill  = 1,                                        &
                       RC        = RC                                       ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Hg2
    HcoID = GetHemcoId( 'Hg2', HcoState, LOC, RC, ERR=.FALSE. )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create diagnostic container
    IF ( HcoID > 0 ) THEN
       DiagnName = 'HG2_ANTHRO'
       CALL Diagn_Create( am_I_Root,                                         & 
                          HcoState  = HcoState,                              &
                          cName     = TRIM(DiagnName),                       &
                          ExtNr     = ExtNr,                                 &
                          Cat       = Cat,                                   &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                    )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! HgP
    HcoID = GetHemcoId( 'HgP', HcoState, LOC, RC, ERR=.FALSE. )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create diagnostic container
    IF ( HcoID > 0 ) THEN
       DiagnName = 'HGP_ANTHRO'
       CALL Diagn_Create( am_I_Root,                                         & 
                          HcoState  = HcoState,                              &
                          cName     = TRIM(DiagnName),                       &
                          ExtNr     = ExtNr,                                 &
                          Cat       = Cat,                                   &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                     ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    !-------------------------------------------
    ! %%%%% BIOMASS BURNING HG %%%%%
    ! ==> defined in Diagn_Biomass
    !-------------------------------------------

    ! Do only if Hg0 is defined ... 
    HcoID = HCO_GetHcoID( 'Hg0', HcoState )
    IF ( HcoID > 0 ) THEN

        Cat   = -1
        ExtNr = GetExtNr( HcoState%Config%ExtList, 'GFED' )
        IF ( ExtNr <= 0 ) ExtNr = GetExtNr( HcoState%Config%ExtList, 'FINN' )
        IF ( ExtNr <= 0 ) THEN
           ExtNr = 0
           Cat   = CATEGORY_BIOMASS
        ENDIF


       ! Create diagnostic container
       DiagnName = 'BIOMASS_HG0'
       CALL Diagn_Create( am_I_Root,                                         & 
                          HcoState  = HcoState,                              &
                          cName     = TRIM(DiagnName),                       &
                          ExtNr     = ExtNr,                                 &
                          Cat       = Cat,                                   &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                    ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Diagn_Hg
!EOC
#if defined( TOMAS )
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_TOMAS
!
! !DESCRIPTION: This creates diagnostics for bulk emissions that will be called
! to scale into TOMAS bins. May not even be necessary. (JKodros 6/2/15)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_TOMAS( am_I_Root, Input_Opt, HcoState, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCO_State_Mod,      ONLY : HCO_GetHcoID
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success

!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, Cat, HcoID, N
    INTEGER            :: id_BCPI, id_BCPO, id_OCPI, id_OCPO
    INTEGER            :: IDSO4
    INTEGER            :: IDCO
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_TOMAS (hcoi_gc_diagn_mod.F90)'
    !=======================================================================
    ! Define ND?? diagnostics (TOMAS-related emissions)
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if the CH4 simulation is not selected
    !IF ( .NOT. ( Input_Opt%ITS_A_CH4_SIM .OR. id_CH4 > 0 ) ) RETURN
    ! SOME SORT OF IF DEFINED TOMAS HERE

    ! Get default HEMCO species ID for BC/OC
    id_BCPI = HCO_GetHcoID( 'BCPI', HcoState )
    id_BCPO = HCO_GetHcoID( 'BCPO', HcoState )
    id_OCPI = HCO_GetHcoID( 'OCPI', HcoState )
    id_OCPO = HCO_GetHcoID( 'OCPO', HcoState )

    ! Extension number is zero (HEMCO core) until defined otherwise
    ExtNr = 0
    !-----------------------------------------------------------------
    ! %%%%% BCPI from anthro (Category 1 or species BCPI_ANTH)  %%%%%
    !-----------------------------------------------------------------

       Cat = CATEGORY_ANTHRO
       HcoID = id_BCPI
       ! Create diagnostic container
       DiagnName = 'BCPI_ANTH'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% BCPO from anthro (Category 1 or species BCPO_ANTH)  %%%%%
    !-----------------------------------------------------------------

       Cat = CATEGORY_ANTHRO
       HcoID = id_BCPO
       ! Create diagnostic container
       DiagnName = 'BCPO_ANTH'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% OCPI from anthro (Category 1 or species OCPI_ANTH)  %%%%%
    !-----------------------------------------------------------------

       Cat = CATEGORY_ANTHRO
       HcoID = id_OCPI
       ! Create diagnostic container
       DiagnName = 'OCPI_ANTH'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% OCPO from anthro (Category 1 or species OCPO_ANTH)  %%%%%
    !-----------------------------------------------------------------

       Cat = CATEGORY_ANTHRO
       HcoID = id_OCPO
       ! Create diagnostic container
       DiagnName = 'OCPO_ANTH'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% BCPI from BIOFUEL (Category 1 or species BCPI_BF)  %%%%%
    !-----------------------------------------------------------------

       Cat = CATEGORY_BIOFUEL
       HcoID = id_BCPI
       ! Create diagnostic container
       DiagnName = 'BCPI_BF'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% BCPO from BIOFUEL (Category 2 or species BCPO_BF)  %%%%%
    !-----------------------------------------------------------------

       Cat = CATEGORY_BIOFUEL
       HcoID = id_BCPO
       ! Create diagnostic container
       DiagnName = 'BCPO_BF'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% OCPI from BIOFUEL (Category 1 or species OCPI_BF)  %%%%%
    !-----------------------------------------------------------------

       Cat = CATEGORY_BIOFUEL
       HcoID = id_OCPI
       ! Create diagnostic container
       DiagnName = 'OCPI_BF'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% OCPO from BIOFUEL (Category 2 or species OCPO_BF)  %%%%%
    !-----------------------------------------------------------------

       Cat = CATEGORY_BIOFUEL
       HcoID = id_OCPO
       ! Create diagnostic container
       DiagnName = 'OCPO_BF'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    ! ------------ NOW DEAL WITH BIOMASS BURNING --------------------
    ! First test if GFED is used.  If not, then test if FINN is used.
    ! If not, then use extension # 0 and the default biomass category.
      Cat   = -1
      ExtNr = GetExtNr( HcoState%Config%ExtList, 'GFED' )
      IF ( ExtNr <= 0 ) ExtNr = GetExtNr( HcoState%Config%ExtList, 'FINN' )
      IF ( ExtNr <= 0 ) THEN
         ExtNr = 0
         Cat   = CATEGORY_BIOMASS
      ENDIF

    !-----------------------------------------------------------------
    ! %%%%% BPCI from BIOB (Category ? or species BCPI_bb)  %%%%%
    !-----------------------------------------------------------------

       HcoID = id_BCPI
       ! Create diagnostic container
       DiagnName = 'BCPI_BB'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% BPCO from BIOB (Category ? or species BCPO_bb)  %%%%%
    !-----------------------------------------------------------------

       HcoID = id_BCPO
       ! Create diagnostic container
       DiagnName = 'BCPO_BB'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% OCPI from BIOB (Category ? or species OCPI_bb)  %%%%%
    !-----------------------------------------------------------------

       HcoID = id_OCPI
       ! Create diagnostic container
       DiagnName = 'OCPI_BB'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% OCPO from BIOB (Category ? or species OCPI_bb)  %%%%%
    !-----------------------------------------------------------------


       HcoID = id_OCPO
       ! Create diagnostic container
       DiagnName = 'OCPO_BB'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% SO4 from ANTRHO (Category ? or species SO4_ANTH)  %%%%%
    !-----------------------------------------------------------------
    ! ID for SO4
    IDSO4 = HCO_GetHcoID( 'SO4', HcoState )

    ! Extension number is zero (HEMCO core) until defined otherwise
    ExtNr = 0

    Cat = CATEGORY_ANTHRO
    HcoID = IDSO4
    ! Create diagnostic container
    DiagnName = 'SO4_ANTH'
    CALL Diagn_Create( am_I_Root,                     &
                       HcoState  = HcoState,          &
                       cName     = TRIM( DiagnName ), &
                       ExtNr     = ExtNr,             &
                       Cat       = Cat,               &
                       Hier      = -1,                &
                       HcoID     = HcoID,             &
                       SpaceDim  = 3,                 &
                       LevIDx    = -1,                &
                       OutUnit   = 'kg/m2/s',         &
                       COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                       AutoFill  = 1,                 &
                       RC        = RC                  )
    IF ( RC /= HCO_SUCCESS ) RETURN

    Cat = CATEGORY_BIOFUEL
    HcoID = IDSO4
    ! Create diagnostic container
    DiagnName = 'SO4_BIOF'
    CALL Diagn_Create( am_I_Root,                     &
                       HcoState  = HcoState,          &
                       cName     = TRIM( DiagnName ), &
                       ExtNr     = ExtNr,             &
                       Cat       = Cat,               &
                       Hier      = -1,                &
                       HcoID     = HcoID,             &
                       SpaceDim  = 2,                 &
                       LevIDx    = -1,                &
                       OutUnit   = 'kg/m2/s',         &
                       COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                       AutoFill  = 1,                 &
                       RC        = RC                  )
    IF ( RC /= HCO_SUCCESS ) RETURN


    !-----------------------------------------------------------------
    ! %%%%% CO from ANTRHO (Category ? or species CO_ANTH)  %%%%%
    !-----------------------------------------------------------------
    ! ID for CO
    IDCO = HCO_GetHcoID( 'CO', HcoState )

    ! Extension number is zero (HEMCO core) until defined otherwise
    ExtNr = 0

    Cat = CATEGORY_ANTHRO
    HcoID = IDCO
    ! Create diagnostic container
    DiagnName = 'CO_ANTH'
    CALL Diagn_Create( am_I_Root,                     &
                       HcoState  = HcoState,          &
                       cName     = TRIM( DiagnName ), &
                       ExtNr     = ExtNr,             &
                       Cat       = Cat,               &
                       Hier      = -1,                &
                       HcoID     = HcoID,             &
                       SpaceDim  = 3,                 &
                       LevIDx    = -1,                &
                       OutUnit   = 'kg/m2/s',         &
                       COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                       AutoFill  = 1,                 &
                       RC        = RC                  )
    IF ( RC /= HCO_SUCCESS ) RETURN

#endif

  END SUBROUTINE Diagn_TOMAS
!EOC
#endif
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetHemcoId
!
! !DESCRIPTION: Function GetHemcoId returns the HEMCO species ID number 
!  corresponding to a given HEMCO species name.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GetHemcoId( HcoName, HcoState, Loc, RC, ERR ) RESULT( HcoID )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
    USE HCO_State_Mod, ONLY : HCO_GetHcoID
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*),  INTENT(IN)  :: HcoName    ! HEMCO species name
    TYPE(HCO_State),   POINTER     :: HcoState   ! HEMCO State object
    CHARACTER(LEN=*),  INTENT(IN)  :: Loc        ! Calling routine
    LOGICAL, OPTIONAL, INTENT(IN)  :: Err        ! Return error if not found
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT) :: RC         ! Success or failure?
!
! !RETURN VALUE:
!
    INTEGER                        :: HcoID      ! HEMCO species ID #
!
! !REMARKS:
!  This is a wrapper function to simplify the code above.   Calls to 
!  HCO_GetHcoId and HCO_Error are made from here. 
! 
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: MSG
    LOGICAL            :: ERROR = .TRUE.

    !=======================================================================
    ! GetHemcoId begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Prompt error?
    IF ( PRESENT(ERR) ) ERROR = ERR

    ! Get the HEMCO species ID from the name
    HcoID = HCO_GetHcoID( HcoName, HcoState )

    ! Exit with error if the species is not valid
    ! (HCO_Error will set RC = HCO_FAIL)
    IF ( HcoID <= 0 .AND. ERROR ) THEN
       MSG = 'This is not a HEMCO species: ' // HcoName
       CALL HCO_Error( MSG, RC, THISLOC=Loc )
    ENDIF

  END FUNCTION GetHemcoId
!EOC
END MODULE HCOI_GC_Diagn_Mod
