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
!  deposition velocities are passed to drydep\_mod.F. Hence, no Acetone
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
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
 
  ! GEOS-Chem diagnostic switches and arrays
  USE CMN_SIZE_Mod
  USE CMN_DIAG_Mod
  USE DIAG_Mod
  USE DIAG53_Mod
  USE DIAG56_Mod

  IMPLICIT NONE
  PRIVATE

  ! Get parameters that define the different categories
#include "hcoi_gc_diagn_include.H"

  ! Define default output frequency. In an ESMF mode, this must be 'Manual'.
#if defined(ESMF_)
  CHARACTER(LEN=16), PARAMETER :: Default_WriteFreq = "Manual"
#else
  CHARACTER(LEN=16), PARAMETER :: Default_WriteFreq = "Daily"
#endif
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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_GetHcoID
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: YesOrNo
    INTEGER            :: I, J,  HcoID, N,    AS
    INTEGER            :: ExtNr, Cat, Hier
    CHARACTER(LEN=15)  :: SpcName
    CHARACTER(LEN=31)  :: DiagnName, WriteFreq
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'HCOI_GC_DIAGN_INIT (hcoi_gc_diagn_mod.F90)'
 
    !=======================================================================
    ! HCOI_GC_DIAGN_INIT begins here!
    !=======================================================================

    ! Assume success
    RC  = HCO_SUCCESS

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

    CALL Diagn_NOSrc   ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
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

    CALL Diagn_POPs    ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_CH4     ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Hg      ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=======================================================================
    ! Define automatic diagnostics (AutoFill)
    !=======================================================================

    ! This is for testing only. Only activate if needed.
   !IF ( .TRUE.  ) THEN     ! Activated
    IF ( .FALSE. ) THEN     ! Deactivated

       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN 

          !-------------------------------------
          ! Hourly emissions
          !-------------------------------------

          ! Do for all emission species
          DO I = 1,37

             ! Get species name
             SELECT CASE ( I )
                CASE ( 1 )
                   SpcName = 'NO'
                CASE ( 2 )
                   SpcName = 'CO'
                CASE ( 3 )
                   SpcName = 'NH3'
                CASE ( 4 )
                   SpcName = 'SO2'
                CASE ( 5 )
                   SpcName = 'SO4'
                CASE ( 6 )
                   SpcName = 'GLYX'
                CASE ( 7 )
                   SpcName = 'MGLY'
                CASE ( 8 )
                   SpcName = 'C2H6'
                CASE ( 9 )
                   SpcName = 'ALK4'
                CASE ( 10 )
                   SpcName = 'ACET'
                CASE ( 11 )
                   SpcName = 'MEK'
                CASE ( 12 )
                   SpcName = 'ALD2'
                CASE ( 13 )
                   SpcName = 'PRPE'
                CASE ( 14 )
                   SpcName = 'C3H8'
                CASE ( 15 )
                   SpcName = 'CH2O'
                CASE ( 16 )
                   SpcName = 'BENZ'
                CASE ( 17 )
                   SpcName = 'TOLU'
                CASE ( 18 )
                   SpcName = 'XYLE'
                CASE ( 19 )
                   SpcName = 'C2H4'
                CASE ( 20 )
                   SpcName = 'C2H2'
                CASE ( 21 )
                   SpcName = 'CHBr3'
                CASE ( 22 )
                   SpcName = 'CH2Br2'
                CASE ( 23 )
                   SpcName = 'BCPI'
                CASE ( 24 )
                   SpcName = 'BCPO'
                CASE ( 25 )
                   SpcName = 'OCPI'
                CASE ( 26 )
                   SpcName = 'OCPO'
                CASE ( 27 )
                   SpcName = 'RCHO'
                CASE ( 28 )
                   SpcName = 'MACR'
                CASE ( 29 )
                   SpcName = 'DMS'
                CASE ( 30 )
                   SpcName = 'DST1'
                CASE ( 31 )
                   SpcName = 'DST2'
                CASE ( 32 )
                   SpcName = 'DST3'
                CASE ( 33 )
                   SpcName = 'DST4'
                CASE ( 34 )
                   SpcName = 'SALA'
                CASE ( 35 )
                   SpcName = 'SALC'
                CASE ( 36 )
                   SpcName = 'Br2'
                CASE ( 37 )
                   SpcName = 'ISOP'
                CASE DEFAULT
                   SpcName = 'DUMMY'
             END SELECT

             HcoID = HCO_GetHcoID( TRIM(SpcName), HcoState )
             IF ( HcoID > 0 ) THEN
                CALL Diagn_Create ( am_I_Root,                          &
                                    HcoState,                           &
                                    cName     = 'EMIS_'//TRIM(SpcName), &
                                    ExtNr     = -1,                     &
                                    Cat       = -1,                     &
                                    Hier      = -1,                     &
                                    HcoID     = HcoID,                  &
                                    SpaceDim  = 2,                      &
                                    LevIDx    = -1,                     &
                                    OutUnit   = 'kg/m2/s',              &
                                    WriteFreq = Default_WriteFreq,      &
                                    AutoFill  = 1,                      &
                                    cID       = N,                      & 
                                 RC        = RC ) 
                IF ( RC /= HCO_SUCCESS ) RETURN
             ENDIF

             !-------------------------------------
             ! Emissions per category (NO only)
             !-------------------------------------
             IF ( TRIM(SpcName) == 'NO' .and. HcoID > 0 ) THEN

                ! There are 3 different categories
                DO J = 1, 6
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
                      CASE DEFAULT
                         DiagnName = 'EMIS_NO_DUMMY'
                         ExtNr     = 999
                         Cat       = 999
                   END SELECT

                   CALL Diagn_Create ( am_I_Root, &
                                       HcoState,  &
                                       cName    = DiagnName, &
                                       ExtNr    = ExtNr, &
                                       Cat      = Cat, &
                                       Hier     = -1, &
                                       HcoID    = HcoID, &
                                       SpaceDim = 2, &
                                       LevIDx   = -1, &
                                       OutUnit  = 'kg/m2/s', &
                                       WriteFreq = Default_WriteFreq, &
                                       AutoFill  = 1,                 &
                                       cID       = N,                 & 
                                       RC        = RC ) 
                   IF ( RC /= HCO_SUCCESS ) RETURN

                ENDDO ! J
             ENDIF    ! NO
          ENDDO       ! I
       ENDIF          ! fullchem
    ENDIF             ! testing toggle

    ! Leave w/ success
    RC = HCO_SUCCESS 

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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
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

    ! Exit if the Rn-Pb-Be simulation is not selected
    IF ( .not. Input_Opt%ITS_A_RnPbBe_SIM ) RETURN

    ! Define diagnostics
    IF ( ExtState%GC_RnPbBe .and. ( ND01 > 0 ) ) THEN

       ! HEMCO extension # for Rn-Pb-Be
       ExtNr = GetExtNr( 'GC_Rn-Pb-Be' )
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
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = -1,                &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/s',            &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
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
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = -1,                &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/s',            &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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

    ! Exit if we are doing a specialty simulation w/o dust
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )   .and. &
         ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
       RETURN
    ENDIF

    ! Define diagnostics if dust is used
    IF ( ( ExtState%DustDead .OR. ExtState%DustGinoux )   .AND. &
         ( ND06 > 0                                   ) ) THEN

       ! Get Ext. Nr of used extension
       ExtNr = GetExtNr( 'DustDead' )
       IF ( ExtNr <= 0 ) ExtNr = GetExtNr( 'DustGinoux' )
       IF ( ExtNr <= 0 ) THEN
          CALL HCO_Error( 'Cannot find dust extension', RC, THISLOC=LOC )
          RETURN      
       ENDIF

       ! Do for each dust bin
       DO I = 1, Input_Opt%N_DUST_BINS

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

          ! Get species name (i.e. DST1 .. DST4 for non TOMAS simualtions
          WRITE( ISTR1,'(i1)' ) I
          SpcName   = 'DST'   // ISTR1

#endif

          DiagnName = 'AD06_' // TRIM( SpcName )         

          ! HEMCO species ID 
          HcoID = GetHemcoId( TRIM( SpcName ), HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create diagnostic container
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg',              &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDDO 
    ENDIF

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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE TRACERID_MOD,       ONLY : IDTPOA1, IDTPOG1
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
!  20 Jan 2015 - M. Yannetti - Added LEMIS Catch Exit
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

    !=======================================================================
    ! DIAGN_CARBON begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Exit if we are doing a specialty simulation w/o carbon aerosols
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )   .and. &
         ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
       RETURN
    ENDIF

    ! Exit if emissions are turned off
    IF ( .not. Input_Opt%LEMIS ) RETURN
   
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
                IF ( IDTPOA1 > 0 ) SpcName = 'POA1'
             CASE ( 4 ) 
                SpcName = 'OCPO'
                IF ( IDTPOG1 > 0 ) SpcName = 'POG1'
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
                CASE ( 2 )
                   SrcName = 'BIOFUEL'
                   Cat     = CATEGORY_BIOFUEL 
             END SELECT

             !-------------------------------------------
             ! %%%%% DEFINE DIAGNOSTICS %%%%%% 
             ! -------------------------------------------

             ! Set DiagnName
             DiagnName = "AD07_"//TRIM(SpcName)//"_"//TRIM(SrcName)

             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState,                      &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = 0,                 &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                WriteFreq = 'Manual',          &
                                AutoFill  = 1,                 &
                                cID       = N,                 & 
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 

          ENDDO
       ENDDO
    ENDIF 

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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, I, N
    CHARACTER(LEN=15)  :: SpcName
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_SEASALT (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_SEASALT begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Exit if we are doing a specialty simulation w/o sea salt
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )   .and. &
         ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
       RETURN
    ENDIF

    ! Define diagnostics
    IF ( ND08 > 0 .AND. Input_Opt%LSSALT .AND. ExtState%SeaSalt ) THEN

       ! Get HEMCO extension # for SeaSalt
       ExtNr = GetExtNr( 'SeaSalt' )
       IF ( ExtNr <= 0 ) THEN
          CALL HCO_Error( 'Cannot find extension SeaSalt', RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! I=1 is SALA, I=2 is SALC
       DO I = 1, 2

          ! Pick the proper species & diagnostic name
          SELECT CASE( I )
             CASE( 1 )
                SpcName   = 'SALA'
                DiagnName = 'AD08_SALA'
             CASE( 2 )
                SpcName   = 'SALC'
                DiagnName = 'AD08_SALC'
          END SELECT

          ! HEMCO species ID 
          HcoID = GetHemcoId( TRIM( SpcName ), HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create the diagnostic
          CALL Diagn_Create ( am_I_Root,                     & 
                              HcoState,                      &
                              cName     = TRIM( DiagnName ), &
                              ExtNr     = ExtNr,             &
                              Cat       = -1,                &
                              Hier      = -1,                &
                              HcoID     = HcoID,             &
                              SpaceDim  = 2,                 &
                              LevIDx    = -1,                &
                              OutUnit   = 'kg',              &
                              WriteFreq = 'Manual',          &
                              AutoFill  = 1,                 &
                              cID       = N,                 & 
                              RC        = RC                  ) 

          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDDO 
    ENDIF

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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
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
!  21 Jan 2015 - M. Yannetti - Added catch for emissions being turned off
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

    ! Exit if emissions are turned off
    IF ( .not. Input_Opt%LEMIS ) RETURN


    ! Exit if we are doing a specialty simulation w/o acetone
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN 
       RETURN
    ENDIF

    ! Define diagnostics
    IF ( ND11 > 0 .OR. ND46 > 0 ) THEN 

       ! Get extension # for SeaFlux 
       ExtNr = GetExtNr( 'SeaFlux' )

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
                              HcoState,                      &
                              cName     = TRIM( DiagnName ), &
                              ExtNr     = ExtNr,             &
                              Cat       = -1,                &
                              Hier      = -1,                &
                              HcoID     = HcoID,             &
                              SpaceDim  = 2,                 &
                              LevIDx    = -1,                &
                              OutUnit   = 'kg/m2/s',         &
                              WriteFreq = 'Manual',          &
                              AutoFill  = 1,                 &
                              cID       = N,                 & 
                              RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF
    ENDIF


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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
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
!  21 Jan 2015 - M. Yannetti - Added catch for emissions being turned off
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

    ! Exit if emissions are turned off
    IF ( .not. Input_Opt%LEMIS ) RETURN


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
       ExtNr = GetExtNr( 'SeaFlux' )

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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg',              &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
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
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_AIRCRAFT, &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! ... anthropogenic ...
       DiagnName = 'AD13_SO2_ANTHROPOGENIC'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! ... biofuel ...
       DiagnName = 'AD13_SO2_BIOFUEL'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_BIOFUEL,  &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! ... from volcanoes ...
       DiagnName = 'AD13_SO2_VOLCANO'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_VOLCANO,  &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! ... from ships ...
       DiagnName = 'AD13_SO2_SHIP'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_SHIP,     &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
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
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! ... biofuel ...
       ExtNr     = 0
       DiagnName = 'AD13_NH3_BIOFUEL'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_BIOFUEL,  &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! ... natural
       ExtNr     = 0
       DiagnName = 'AD13_NH3_NATURAL'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_NATURAL,  &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
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
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! ... biofuel ...
       ExtNr     = 0
       DiagnName = 'AD13_SO4_BIOFUEL' 
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_BIOFUEL,  &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCO_State_Mod,      ONLY : HCO_GetHcoID
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE TRACERID_MOD,       ONLY : IDTPOA1
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
!  21 Jan 2015 - M. Yannetti - Added catch for emissions being turned off
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: Cat, ExtNr, HcoID, N
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_BIOMASS (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_BIOMASS begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Exit if emissions are turned off
    IF ( .not. Input_Opt%LEMIS ) RETURN

    ! Exit if we are doing a specialty simulation w/o biomass
    IF ( Input_Opt%ITS_A_POPS_SIM   ) RETURN
    IF ( Input_Opt%ITS_A_RnPbBe_SIM ) RETURN
    IF ( Input_Opt%ITS_A_TAGOX_SIM  ) RETURN

    ! First test if GFED3 is used.  If not, then test if FINN is used.
    ! If not, then use extension # 0 and the default biomass category.
    Cat   = -1
    ExtNr = GetExtNr( 'GFED3' )
    IF ( ExtNr <= 0 ) ExtNr = GetExtNr( 'FINN' )
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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
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

          !----------------------------------------
          ! %%%%% FOR POA SIMULATION %%%%%
          !----------------------------------------
          IF ( IDTPOA1 > 0 ) THEN
             ! HEMCO species ID
             HcoID = GetHemcoId( 'POA1', HcoState, LOC, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
             
             ! Create diagnostic container
             DiagnName = 'BIOMASS_POA1'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState,                      &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                WriteFreq = 'Manual',          &
                                AutoFill  = 1,                 &
                                cID       = N,                 & 
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN

          !----------------------------------------
          ! %%%%% Biomass OC %%%%%
          !----------------------------------------
          ELSE
          
             ! HEMCO species ID
             HcoID = GetHemcoId( 'OCPI', HcoState, LOC, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
             
             ! Create diagnostic container
             DiagnName = 'BIOMASS_OCPI'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState,                      &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                WriteFreq = 'Manual',          &
                                AutoFill  = 1,                 &
                                cID       = N,                 & 
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
             
             ! HEMCO species ID
             HcoID = GetHemcoId( 'OCPO', HcoState, LOC, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
             
             ! Create diagnostic container
             DiagnName = 'BIOMASS_OCPO'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState,                      &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                WriteFreq = 'Manual',          &
                                AutoFill  = 1,                 &
                                cID       = N,                 & 
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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! HEMCO species ID
          HcoID = GetHemcoId( 'BCPO', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          
          ! Create diagnostic container
          DiagnName = 'BIOMASS_BCPO'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
    ENDIF

    !----------------------------------------------
    ! %%%%% Biomass CO %%%%%
    !----------------------------------------------
    IF ( ND28 > 0 .OR. ND29 > 0 ) THEN

       ! CO is only defined for full chemistry and tagged CO simulations
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or.             &
            Input_Opt%ITS_A_TAGCO_SIM    ) THEN

          ! HEMCO species ID
          HcoID = GetHemcoId( 'CO', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create diagnostic container
          DiagnName = 'BIOMASS_CO'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
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
                             HcoState,                    &
                             cName     = TRIM(DiagnName), &
                             ExtNr     = ExtNr,           &
                             Cat       = Cat,             &
                             Hier      = -1,              &
                             HcoID     = HcoID,           &
                             SpaceDim  = 2,               &
                             LevIDx    = -1,              &
                             OutUnit   = 'kg/m2/s',       &
                             WriteFreq = 'Manual',        &
                             AutoFill  = 1,               &
                             cID       = N,               & 
                             RC        = RC                ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
    ENDIF

    !----------------------------------------------
    ! %%%%% Biomass Hg0 %%%%%
    !----------------------------------------------

    ! Do only if Hg0 is defined ... 
    HcoID = HCO_GetHcoID( 'Hg0', HcoState )
    IF ( HcoID > 0 ) THEN

       ! Create diagnostic container
       DiagnName = 'BIOMASS_HG0'
       CALL Diagn_Create( am_I_Root,                   & 
                          HcoState,                    &
                          cName     = TRIM(DiagnName), &
                          ExtNr     = ExtNr,           &
                          Cat       = Cat,             &
                          Hier      = -1,              &
                          HcoID     = HcoID,           &
                          SpaceDim  = 2,               &
                          LevIDx    = -1,              &
                          OutUnit   = 'kg/m2/s',       &
                          WriteFreq = 'Manual',        &
                          AutoFill  = 1,               &
                          cID       = N,               & 
                          RC        = RC                ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt
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
!  22 Jan 2015 - M. Yannetti - Corrected typo in LOC, exit for emissions off
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

    ! Exit if emissions are turned off
    IF ( .not. Input_Opt%LEMIS ) RETURN


    ! Exit if we are doing a specialty simulation w/o NO
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

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
                       HcoState,                       &
                       cName     = TRIM( DiagnName ),  &
                       ExtNr     = ExtNr,              &
                       Cat       = CATEGORY_AIRCRAFT,  &
                       Hier      = -1,                 &
                       HcoID     = HcoID,              &
                       SpaceDim  = 3,                  &
                       LevIDx    = -1,                 &
                       OutUnit   = 'kg/m2/s',          &
                       WriteFreq = 'Manual',           &
                       AutoFill  = 1,                  &
                       cID       = N,                  & 
                       RC        = RC                   )
    IF ( RC /= HCO_SUCCESS ) RETURN 


    !----------------------------------------------
    ! %%%%% Ship NO %%%%%
    !
    ! ==> Only define if ParaNOx is not used. 
    !     SHIP_NO from ParaNOx is defined in ND63.
    !----------------------------------------------
    Cat   = -1
    ExtNr = GetExtNr( 'ParaNOx' )

    IF ( ExtNr <= 0 ) THEN
       ExtNr     = 0
       Cat       = CATEGORY_SHIP
       DiagnName = 'SHIP_NO'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !----------------------------------------------
    ! %%%%% Lightning NO %%%%%
    !
    ! ==> Only define if LightNox is turned on
    !----------------------------------------------
    Cat   = -1
    ExtNr = GetExtNr('LightNOx')
    IF ( ExtNr > 0 ) THEN
       DiagnName = 'LIGHTNING_NO'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF
      

    !----------------------------------------------
    ! %%%%% Soil and Fertilizer NO %%%%%
    !
    ! ==> Only define if SoilNox is turned on
    !----------------------------------------------
    Cat   = -1
    ExtNr = GetExtNr('SoilNOx')
    IF ( ExtNr > 0 ) THEN

       ! %%%%%% Soil NO %%%%%%
       DiagnName = 'SOIL_NO'
       CALL Diagn_Create ( am_I_Root,                     & 
                           HcoState,                      &
                           cName     = TRIM( DiagnName ), &
                           ExtNr     = ExtNr,             &
                           Cat       = Cat,               &
                           Hier      = -1,                &
                           HcoID     = HcoID,             &
                           SpaceDim  = 2,                 &
                           LevIDx    = -1,                &
                           OutUnit   = 'kg/m2/s',         &
                           WriteFreq = 'Manual',          &
                           AutoFill  = 1,                 &
                           cID       = N,                 & 
                           RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! %%%%%% Fertilizer NO %%%%%%
       CALL GetExtOpt( ExtNr, 'Use fertilizer', OptValBool=YesOrNo, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       IF ( YesOrNo == .FALSE. ) THEN
          MSG = 'Fertilizer NOx disabled - diagnostics will be zero!'
          CALL HCO_Warning( MSG, RC, THISLOC=LOC )
       ENDIF

       DiagnName = 'FERTILIZER_NO'
       CALL Diagn_Create ( am_I_Root,                     & 
                           HcoState,                      &
                           cName     = TRIM( DiagnName ), &
                           ExtNr     = ExtNr,             &
                           Cat       = Cat,               &
                           Hier      = -1,                &
                           HcoID     = HcoID,             &
                           SpaceDim  = 2,                 &
                           LevIDx    = -1,                &
                           OutUnit   = 'kg/m2/s',         &
                           WriteFreq = 'Manual',          &
                           AutoFill  = 0,                 &
                           cID       = N,                 & 
                           RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_GetHcoId
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
!  22 Jan 2015 - M. Yannetti - Added LEMIS Catch Exit
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, N, I
    CHARACTER(LEN=31)  :: DiagnName, SpcName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_BIOFUEL (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_BIOFUEL begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Exit if emissions are turned off
    IF ( .not. Input_Opt%LEMIS ) RETURN


    ! Exit if we are doing a specialty simulation w/o biofuels
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) RETURN
    IF ( Input_Opt%ITS_A_POPS_SIM    ) RETURN
    IF ( Input_Opt%ITS_A_TAGOX_SIM   ) RETURN
    IF ( Input_Opt%ITS_A_RnPbBe_SIM  ) RETURN

    ! Extension number
    ExtNr = 0
 
    ! ND34 only
    IF ( ND34 > 0 ) THEN

       ! Loop over speices
       DO I = 1, 10

          ! Select species
          SELECT CASE ( I ) 
             CASE ( 1 )
                SpcName = 'SO2'
             CASE ( 2 )
                SpcName = 'NH3'
             CASE ( 3 )
                SpcName = 'ALK4'
             CASE ( 4 )
                SpcName = 'ALD2'
             CASE ( 5 )
                SpcName = 'ACET'
             CASE ( 6 )
                SpcName = 'MEK'
             CASE ( 7 )
                SpcName = 'PRPE'
             CASE ( 8 )
                SpcName = 'C2H6'
             CASE ( 9 )
                SpcName = 'C3H8'
             CASE ( 10)
                SpcName = 'CH2O'

             CASE DEFAULT
                SpcName = 'DUMMY'
          END SELECT

          ! Check for species ID
          HcoID = HCO_GetHcoId( TRIM(SpcName), HcoState )
          IF ( HcoID > 0 ) THEN 
             ! Create diagnostic container
             DiagnName = 'BIOFUEL_'//TRIM(SpcName)
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState,                      &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = CATEGORY_BIOFUEL,  &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                WriteFreq = 'Manual',          &
                                AutoFill  = 1,                 &
                                cID       = N,                 & 
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
                            HcoState,                      &
                            cName     = TRIM( DiagnName ), &
                            ExtNr     = ExtNr,             &
                            Cat       = CATEGORY_BIOFUEL,  &
                            Hier      = -1,                &
                            HcoID     = HcoID,             &
                            SpaceDim  = 2,                 &
                            LevIDx    = -1,                &
                            OutUnit   = 'kg/m2/s',         &
                            WriteFreq = 'Manual',          &
                            AutoFill  = 1,                 &
                            cID       = N,                 & 
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
                            HcoState,                      &
                            cName     = TRIM( DiagnName ), &
                            ExtNr     = ExtNr,             &
                            Cat       = CATEGORY_BIOFUEL,  &
                            Hier      = -1,                &
                            HcoID     = HcoID,             &
                            SpaceDim  = 2,                 &
                            LevIDx    = -1,                &
                            OutUnit   = 'kg/m2/s',         &
                            WriteFreq = 'Manual',          &
                            AutoFill  = 1,                 &
                            cID       = N,                 & 
                            RC        = RC                  ) 
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF
   ENDIF

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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
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
!  22 Jan 2015 - M. Yannetti - Added LEMIS Catch Exit
!EOC
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, I, N
    CHARACTER(LEN=15)  :: SpcName
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_ANTHRO (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_Anthro begins here!
    !=======================================================================

    ! Extension number
    ExtNr = 0

    ! Assume success
    RC = HCO_SUCCESS

    ! Exit if emissions are turned off
    IF ( .not. Input_Opt%LEMIS ) RETURN

    ! Exit if we are doing a specialty simulation w/o the anthro species below
    IF ( Input_Opt%ITS_A_C2H6_SIM    ) RETURN
    IF ( Input_Opt%ITS_A_CH4_SIM     ) RETURN
    IF ( Input_Opt%ITS_A_CO2_SIM     ) RETURN
    IF ( Input_Opt%ITS_A_HCN_SIM     ) RETURN
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) RETURN
    IF ( Input_Opt%ITS_A_POPS_SIM    ) RETURN
    IF ( Input_Opt%ITS_A_RnPbBe_SIM  ) RETURN
    IF ( Input_Opt%ITS_A_TAGOX_SIM   ) RETURN

    ! ND36 only: VOC's are only defined for fullchem (not tagged CO)
    IF ( ND36 > 0 .and. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       !----------------------------------------
       ! %%%%% Anthropogenic ALK4 %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'ALK4', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'ANTHROPOGENIC_ALK4'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
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
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
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
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
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
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
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
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
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
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
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
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
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
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = CATEGORY_ANTHRO,   &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = CATEGORY_ANTHRO,   &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
    ENDIF

    !-------------------------------------------
    ! %%%%% Anthropogenic CO %%%%%
    !-------------------------------------------
    IF ( ND36 > 0 .OR. ND29 > 0 ) THEN

       ! CO is only defined for the full-chemistry simulation
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM    .or.            &
            Input_Opt%ITS_A_TAGCO_SIM     ) THEN

          ! HEMCO species ID
          HcoID = GetHemcoId( 'CO', HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          
          ! Create diagnostic container
          DiagnName = 'ANTHROPOGENIC_CO'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = CATEGORY_ANTHRO,   &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
    ENDIF

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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCO_State_Mod,      ONLY : HCO_GetHcoId
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt
    USE TRACERID_MOD,       ONLY : IDTPOA1
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
!  22 Jan 2015 - M. Yannetti - Added LEMIS Catch Exit
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

    ! Exit if emissions are turned off
    IF ( .not. Input_Opt%LEMIS ) RETURN

    ! Exit if we are doing a specialty simulation w/o biofuels
    IF ( Input_Opt%ITS_A_HCN_SIM     ) RETURN
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) RETURN
    IF ( Input_Opt%ITS_A_POPS_SIM    ) RETURN
    IF ( Input_Opt%ITS_A_RnPbBe_SIM  ) RETURN
    IF ( Input_Opt%ITS_A_TAGOX_SIM   ) RETURN

    ! Extension and category #'s for MEGAN
    ExtNr = GetExtNr('MEGAN')
    Cat   = -1

    ! Make sure MEGAN is on if ND46 is used
    IF ( ExtNr <= 0 .AND. ND46 > 0 ) THEN
       MSG = 'MEGAN is not enabled - cannot write biogenic diagnostics!'
       CALL HCO_Error ( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
 
    ! Only if MEGAN is on ... 
    IF ( ExtNr > 0 ) THEN

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
                                 HcoState,                      &
                                 cName     = TRIM( DiagnName ), &
                                 ExtNr     = ExtNr,             &
                                 Cat       = Cat,               &
                                 Hier      = -1,                &
                                 HcoID     = HcoID,             &
                                 SpaceDim  = 2,                 &
                                 LevIDx    = -1,                &
                                 OutUnit   = 'kg/m2/s',         &
                                 WriteFreq = 'Manual',          &
                                 AutoFill  = 1,                 &
                                 cID       = N,                 & 
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
                                HcoState,                      &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                WriteFreq = 'Manual',          &
                                AutoFill  = 1,                 &
                                cID       = N,                 & 
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
                                HcoState,                      &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                WriteFreq = 'Manual',          &
                                AutoFill  = 1,                 &
                                cID       = N,                 & 
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
                                HcoState,                      &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = 0,                 &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                WriteFreq = 'Manual',          &
                                AutoFill  = 1,                 &
                                cID       = N,                 & 
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
                                HcoState,                      &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = 0,                 &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                WriteFreq = 'Manual',          &
                                AutoFill  = 1,                 &
                                cID       = N,                 & 
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF
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
                                HcoState,                      &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                WriteFreq = 'Manual',          &
                                AutoFill  = 1,                 &
                                cID       = N,                 & 
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF
   
          !%%% For ND11 only %%%
          !%%% These diagnostics are explicitly filled in hcox_megan_mod. 
          !%%% The diagnostics  name defined below must match the names 
          !%%% used in the MEGAN extension!
          IF ( ND11 > 0 ) THEN
   
             ! There are three manual diagnostics in MEGAN
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
                                   HcoState,                      &
                                   cName     = TRIM( DiagnName ), &
                                   ExtNr     = ExtNr,             &
                                   Cat       = -1,                &
                                   Hier      = -1,                &
                                   HcoID     = HcoID,             &
                                   SpaceDim  = 2,                 &
                                   LevIDx    = -1,                &
                                   OutUnit   = 'kg/m2/s',         &
                                   WriteFreq = 'Manual',          &
                                   AutoFill  = 0,                 &
                                   cID       = N,                 & 
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
    ExtNr = GetExtNr('MEGAN_Mono')
    IF ( ExtNr > 0 ) THEN

       !%%% For ND46 diagnostic %%%
       IF ( ND46 > 0 ) THEN

          !----------------------------------------
          ! %%%%% Biogenic MONX %%%%%
          !----------------------------------------

          ! HEMCO species ID
          HcoID = HCO_GetHcoID( 'MONX', HcoState )

          ! Create diagnostic container (if MONX is defined)
          IF ( HcoID > 0 ) THEN
             DiagnName = 'BIOGENIC_MONX'
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState,                      &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                WriteFreq = 'Manual',          &
                                AutoFill  = 1,                 &
                                cID       = N,                 & 
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN 
          ENDIF
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
                                HcoState,                      &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &  
                                Cat       = Cat,               &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg/m2/s',         &
                                WriteFreq = 'Manual',          &
                                AutoFill  = 1,                 &
                                cID       = N,                 & 
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
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF
    ENDIF ! Megan mono

    !=======================================================================
    ! These diagnostics use the MEGAN SOA extension
    !=======================================================================
    IF ( Input_Opt%LSOA .AND. ( ND46 > 0 .OR. ND07 > 0 ) ) THEN

       ! Extension # of MEGAN SOA
       ExtNr = GetExtNr('MEGAN_SOA')
       IF ( ExtNr < 0 ) THEN
          MSG = 'MEGAN SOA emissions are not turned on!'
          CALL HCO_Error( MSG, RC, THISLOC=LOC )
          RETURN      
       ENDIF
         
       !----------------------------------------
       ! %%%%% Biogenic MTPA %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'MTPA', HcoState )

       IF ( HcoID > 0 ) THEN
          ! Create diagnostic container
          CALL Diagn_Create( am_I_Root,                    & 
                             HcoState,                     &
                             cName     = 'BIOGENIC_MTPA',  &
                             ExtNr     = ExtNr,            &
                             Cat       = -1,               &
                             Hier      = -1,               &
                             HcoID     = HcoID,            &
                             SpaceDim  = 2,                &
                             LevIDx    = -1,               &
                             OutUnit   = 'kg/m2/s',        &
                             WriteFreq = 'Manual',         &
                             AutoFill  = 1,                &
                             cID       = N,                & 
                             RC        = RC                 ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF

       !----------------------------------------
       ! %%%%% Biogenic MTPO %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoId( 'MTPO', HcoState )

       IF ( HcoID > 0 ) THEN
          ! Create diagnostic container
          CALL Diagn_Create( am_I_Root,                    & 
                             HcoState,                     &
                             cName     = 'BIOGENIC_MTPO',  &
                             ExtNr     = ExtNr,            &
                             Cat       = -1,               &
                             Hier      = -1,               &
                             HcoID     = HcoID,            &
                             SpaceDim  = 2,                &
                             LevIDx    = -1,               &
                             OutUnit   = 'kg/m2/s',        &
                             WriteFreq = 'Manual',         &
                             AutoFill  = 1,                &
                             cID       = N,                & 
                             RC        = RC                 ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF
       !----------------------------------------
       ! %%%%% Biogenic LIMO %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoId( 'LIMO', HcoState )

       IF ( HcoID > 0 ) THEN
          ! Create diagnostic container
          CALL Diagn_Create( am_I_Root,                    & 
                             HcoState,                     &
                             cName     = 'BIOGENIC_LIMO',  &
                             ExtNr     = ExtNr,            &
                             Cat       = -1,               &
                             Hier      = -1,               &
                             HcoID     = HcoID,            &
                             SpaceDim  = 2,                &
                             LevIDx    = -1,               &
                             OutUnit   = 'kg/m2/s',        &
                             WriteFreq = 'Manual',         &
                             AutoFill  = 1,                &
                             cID       = N,                & 
                             RC        = RC                 ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF

       !----------------------------------------
       ! %%%%% Biogenic C2H6 %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( 'SESQ', HcoState ) 

       IF ( HcoID > 0 ) THEN
          ! Create diagnostic container
          CALL Diagn_Create( am_I_Root,                    & 
                             HcoState,                     &
                             cName     = 'BIOGENIC_SESQ',  &
                             ExtNr     = ExtNr,            &
                             Cat       = -1,               &
                             Hier      = -1,               &
                             HcoID     = HcoID,            &
                             SpaceDim  = 2,                &
                             LevIDx    = -1,               &
                             OutUnit   = 'kg/m2/s',        &
                             WriteFreq = 'Manual',         &
                             AutoFill  = 1,                &
                             cID       = N,                & 
                             RC        = RC                 ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF
    ENDIF ! SOA simulation

    !=======================================================================
    ! These diagnostics use the SeaSalt extension
    !=======================================================================
    IF ( ND46 > 0 ) THEN

       ! Extension # of SeaSalt
       ExtNr = GetExtNr('SeaSalt')
       IF ( ExtNr <= 0 ) THEN
          CALL HCO_Error ( 'SeaSalt extension not enabled', RC, THISLOC=LOC )
          RETURN      
       ENDIF

       ! Find out if SeaSalt Br2 is enabled
       CALL GetExtOpt ( ExtNr, 'Emit Br2', OptValBool=YesOrNo, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( YesOrNo == .FALSE. ) THEN
          CALL HCO_Error ( 'SeaSalt Br2 not enabled', RC, THISLOC=LOC )
          RETURN      
       ENDIF

       !----------------------------------------
       ! %%%%% Biogenic Br2 %%%%%
       !----------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'Br2', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       CALL Diagn_Create( am_I_Root,                    & 
                          HcoState,                     &
                          cName     = 'SEASALT_BR2',    &
                          ExtNr     = ExtNr,            &
                          Cat       = -1,               &
                          Hier      = -1,               &
                          HcoID     = HcoID,            &
                          SpaceDim  = 2,                &
                          LevIDx    = -1,               &
                          OutUnit   = 'kg/m2/s',        &
                          WriteFreq = 'Manual',         &
                          AutoFill  = 1,                &
                          cID       = N,                & 
                          RC        = RC                 ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
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
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use lightning
!  22 Jan 2015 - M. Yannetti - Added LEMIS Catch Exit
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, I, N
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

    ! Exit if emissions are turned off
    IF ( .not. Input_Opt%LEMIS ) RETURN

    ! Exit if we are doing a specialty simulation w/o lightning
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

    ! Define diagnostics
    IF ( ND56 > 0 ) THEN

       ! Extension number
       ExtNr = GetExtNr('LightNOx')

       ! Exit if LightNOx was not turned on
       IF ( ExtNr <= 0 ) THEN
          MSG = 'Lightning NOx is not enabled - cannot write diagnostics ND56!'
          CALL HCO_Error( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! HEMCO species ID
       HcoID = GetHemcoId( 'NO', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

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

          ! Create diagnostic container
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'flashes/min/km2', &
                             OutOper   = 'Mean',            &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 0,                 &
                             cID       = N,                 & 
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDDO
    ENDIF 

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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
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

    ! Exit if we are doing a specialty simulation w/o lightning
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

    ! Extension number
    ExtNr = GetExtNr('ParaNOx')

    ! Exit if PARANOX extension was turned off
    IF ( ExtNr <= 0 .AND. ND63 > 0 ) THEN
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
       IF ( ND63 > 0 .OR. ND32 > 0 ) THEN
          
          ! Create diagnostic container
          DiagnName = 'SHIP_NO'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 &
                             cID       = N,                 & 
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
 
       !-------------------------------------------
       ! %%%%% PARANOX specific diagnostics %%%%%
       !-------------------------------------------
       IF ( ND63 > 0 ) THEN
 
          ! These are the ship NO emissions before PARANOX chemistry
          DiagnName = 'PARANOX_TOTAL_SHIPNOX'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 0,                 &
                             cID       = N,                 & 
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN
  
          ! These is the O3 production/loss through PARANOX 
          DiagnName = 'PARANOX_O3_PRODUCTION'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 0,                 &
                             cID       = N,                 & 
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN
 
          DiagnName = 'PARANOX_NOXFRAC_REMAINING'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'unitless',        &
                             OutOper   = 'Mean',            &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 0,                 &
                             cID       = N,                 & 
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN
  
          DiagnName = 'PARANOX_OPE'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'unitless',        &
                             OutOper   = 'Mean',            &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 0,                 &
                             cID       = N,                 & 
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN
 
       ENDIF
    ENDIF

  END SUBROUTINE Diagn_ParaNOx
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_POPs
!
! !DESCRIPTION: Subroutine Diagn\_POPs initializes diagnostics for the
!  POPs simulation (ND53).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_POPs( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, N
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_POPs (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! Define ND53 diagnostics (POPs emissions)
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Exit if the POPs simulation is not selected
    IF ( .not. Input_Opt%ITS_A_POPS_SIM ) RETURN

    ! Define diagnostics
    IF ( ExtState%GC_POPs .and. ( ND53 > 0 ) ) THEN

       ! HEMCO extension # for POPs
       ExtNr = GetExtNr( 'GC_POPs' )
       IF ( ExtNr <= 0 ) THEN
          CALL HCO_Error ( 'Cannot find POPs extension', RC, THISLOC=LOC )
          RETURN      
       ENDIF

!------------------------------------------------------------------------------
! Prior to 8/27/14:
! Comment out for now -- Need to figure out how to track total POPs emissions
! in HEMCO (mps/8/27/14)
!       !-------------------------------------------
!       ! %%%%% Total POP %%%%%
!       !-------------------------------------------
! 
!       ! HEMCO species ID
!       HcoID = GetHemcoId( 'POPG', HcoState, LOC, RC )
!       IF ( RC /= HCO_SUCCESS ) RETURN
!
!       ! Create diagnostic container
!       DiagnName = 'AD01_POPT_SOURCE'
!       CALL Diagn_Create( am_I_Root,                     & 
!                          HcoState,                      &
!                          cName     = TRIM( DiagnName ), &
!                          ExtNr     = ExtNr,             &
!                          Cat       = -1,                &
!                          Hier      = -1,                &
!                          HcoID     = HcoID,             &
!                          SpaceDim  = 2,                 &
!                          LevIDx    = -1,                &
!                          OutUnit   = 'kg',              &
!                          WriteFreq = 'Manual',          &
!                          AutoFill  = 1,                 &
!                          cID       = N,                 & 
!                          RC        = RC                  ) 
!       IF ( RC /= HCO_SUCCESS ) RETURN 
!------------------------------------------------------------------------------

       !-------------------------------------------
       ! %%%%% OC-phase POP %%%%%
       !-------------------------------------------
 
       ! HEMCO species ID
       HcoID = GetHemcoId( 'POPPOC', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'AD01_POPPOC_SOURCE'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = -1,                &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       !-------------------------------------------
       ! %%%%% BC-phase POP %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'POPPBC', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'AD01_POPPBC_SOURCE'
       CALL Diagn_Create( am_I_Root,                   & 
                          HcoState,                    &
                          cName     = TRIM(DiagnName), &
                          ExtNr     = ExtNr,           &
                          Cat       = -1,              &
                          Hier      = -1,              &
                          HcoID     = HcoID,           &
                          SpaceDim  = 2,               &
                          LevIDx    = -1,              &
                          OutUnit   = 'kg',            &
                          WriteFreq = 'Manual',        &
                          AutoFill  = 1,               &
                          cID       = N,               & 
                          RC        = RC                ) 
       IF ( RC /= HCO_SUCCESS ) RETURN


       !-------------------------------------------
       ! %%%%% Gas-phase POP %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'POPG', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'AD01_POPG_SOURCE'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = -1,                &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg',              &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCO_State_Mod,      ONLY : HCO_GetHcoID
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE TRACERID_MOD,       ONLY : IDTCH4
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, IDCH4, Cat, HcoID, N
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_CH4 (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! Define ND58 diagnostics (CH4 emissions)
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Exit if the CH4 simulation is not selected
    IF ( .NOT. ( Input_Opt%ITS_A_CH4_SIM .OR. IDTCH4 > 0 ) ) RETURN

    ! Get default HEMCO species ID for CH4 
    IDCH4 = HCO_GetHcoID( 'CH4', HcoState )

    ! Extension number is zero (HEMCO core) until defined otherwise
    ExtNr = 0

    !-----------------------------------------------------------------
    ! %%%%% CH4 from gas and oil (Category 1 or species CH4_ga)  %%%%%
    !-----------------------------------------------------------------

    ! Check if there is a specific HEMCO species defined for this 
    ! category, in which case we use the total of this species.
    ! Otherwise, use CH4 category 1 emissions.
    HcoID = HCO_GetHcoID( 'CH4_ga', HcoState )
    IF ( HcoID > 0 ) THEN
       Cat   = -1
    ELSE
       HcoID = IDCH4
       Cat   = 1
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_GAS'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !----------------------------------------------------------
    ! %%%%% CH4 from coal (Category 2 or species CH4_co)  %%%%%
    !----------------------------------------------------------

    HcoID = HCO_GetHcoID( 'CH4_co', HcoState )
    IF ( HcoID > 0 ) THEN
       Cat   = -1
    ELSE
       HcoID = IDCH4
       Cat   = 2
    ENDIF
    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_COAL'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !---------------------------------------------------------------
    ! %%%%% CH4 from livestock (Category 3 or species CH4_ef)  %%%%%
    !---------------------------------------------------------------

    HcoID = HCO_GetHcoID( 'CH4_ef', HcoState )
    IF ( HcoID > 0 ) THEN
       Cat   = -1
    ELSE
       HcoID = IDCH4
       Cat   = 3
    ENDIF
    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_LIVESTOCK'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF   
 
    !---------------------------------------------------------------
    ! %%%%% CH4 from waste (Category 4 or species CH4_wa)  %%%%%
    !---------------------------------------------------------------

    HcoID = HCO_GetHcoID( 'CH4_wa', HcoState )
    IF ( HcoID > 0 ) THEN
       Cat   = -1
    ELSE
       HcoID = IDCH4
       Cat   = 4
    ENDIF
    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_WASTE'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF   
 
    !---------------------------------------------------------------
    ! %%%%% CH4 from biofuel (Category 5 or species CH4_bf)  %%%%%
    !---------------------------------------------------------------

    HcoID = HCO_GetHcoID( 'CH4_bf', HcoState )
    IF ( HcoID > 0 ) THEN
       Cat   = -1
    ELSE
       HcoID = IDCH4
       Cat   = 5
    ENDIF
    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_BIOFUEL'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from other anth. sources (Category 6 or species CH4_oa)  %%%%%
    !-------------------------------------------------------------------------

    HcoID = HCO_GetHcoID( 'CH4_oa', HcoState )
    IF ( HcoID > 0 ) THEN
       Cat   = -1
    ELSE
       HcoID = IDCH4
       Cat   = 6
    ENDIF
    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_ANTHROTHER'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from soil absorption (Category 7 or species CH4_sa)  %%%%%
    !-------------------------------------------------------------------------

    HcoID = HCO_GetHcoID( 'CH4_sa', HcoState )
    IF ( HcoID > 0 ) THEN
       Cat   = -1
    ELSE
       HcoID = IDCH4
       Cat   = 7
    ENDIF
    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_SOILABSORB'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !---------------------------------------------------------------------------
    ! %%%%% CH4 from other natural sources (Category 8 or species CH4_on)  %%%%%
    !---------------------------------------------------------------------------

    HcoID = HCO_GetHcoID( 'CH4_on', HcoState )
    IF ( HcoID > 0 ) THEN
       Cat   = -1
    ELSE
       HcoID = IDCH4
       Cat   = 8
    ENDIF
    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_OTHERNATUR'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 1,                 &
                          cID       = N,                 & 
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !--------------------------------------------------------------------------
    ! %%%%% CH4 from biomass burning (automatically filled in extension)  %%%%%
    !--------------------------------------------------------------------------

    ! HEMCO extension # for wetland ch4 
    ExtNr = GetExtNr( 'GFED3' )
    IF ( ExtNr <= 0 ) ExtNr = GetExtNr( 'FINN' )
    IF ( ExtNr <= 0 ) THEN
       CALL HCO_Warning ( 'Biomass burning emissions not turned on!!', RC, THISLOC=LOC )
    ENDIF
    IF ( ExtNr > 0 ) THEN
       IF ( IDCH4 < 0 ) THEN
          HcoID = HCO_GetHcoID( 'CH4_tot', HcoState )
       ELSE
          HcoID = IDCH4
       ENDIF
       IF ( HcoID > 0 ) THEN 
   
          ! Create diagnostic container
          DiagnName = 'CH4_BIOMASS'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 1,                 & 
                             cID       = N,                 & 
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF     
    ENDIF     
   
    !----------------------------------------------------------------------
    ! %%%%% CH4 from rice (manual diagnostics in wetlands extension)  %%%%%
    !----------------------------------------------------------------------

    ! HEMCO extension # for wetland ch4 
    ExtNr = GetExtNr( 'CH4_WETLANDS' )
    IF ( ExtNr <= 0 ) THEN
       CALL HCO_Warning ( 'Wetland emissions not turned on!!', RC, THISLOC=LOC )
    ENDIF
    IF ( ExtNr > 0 ) THEN
       IF ( IDCH4 < 0 ) THEN
          HcoID = HCO_GetHcoID( 'CH4_tot', HcoState )
       ELSE
          HcoID = IDCH4
       ENDIF
       IF ( HcoID > 0 ) THEN 
   
          ! Create diagnostic container
          DiagnName = 'CH4_RICE'
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState,                      &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             WriteFreq = 'Manual',          &
                             AutoFill  = 0,                 &  ! Manually filled !!
                             cID       = N,                 & 
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 
       ENDIF     
    ENDIF     
   
    !--------------------------------------------------------------------------
    ! %%%%% CH4 from wetlands (manual diagnostics in wetlands extension)  %%%%%
    !--------------------------------------------------------------------------

    IF ( ExtNr > 0 .AND. HcoID > 0 ) THEN
       DiagnName = 'CH4_WETLAND'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState,                      &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = -1,                &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          WriteFreq = 'Manual',          &
                          AutoFill  = 0,                 &  ! Manually filled !!
                          cID       = N,                 & 
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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
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
    CALL Diagn_Create( am_I_Root,                   & 
                       HcoState,                    &
                       cName     = TRIM(DiagnName), &
                       ExtNr     = ExtNr,           &
                       Cat       = Cat,             &
                       Hier      = -1,              &
                       HcoID     = HcoID,           &
                       SpaceDim  = 2,               &
                       LevIDx    = -1,              &
                       OutUnit   = 'kg/m2/s',       &
                       WriteFreq = 'Manual',        &
                       AutoFill  = 1,               &
                       cID       = N,               & 
                       RC        = RC                ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-------------------------------------------
    ! %%%%% NATURAL HG (Hg0) %%%%%
    !-------------------------------------------
 
    ! Create diagnostic container
    DiagnName = 'HG0_NATURAL'
    Cat       = CATEGORY_NATURAL
    CALL Diagn_Create( am_I_Root,                   & 
                       HcoState,                    &
                       cName     = TRIM(DiagnName), &
                       ExtNr     = ExtNr,           &
                       Cat       = Cat,             &
                       Hier      = -1,              &
                       HcoID     = HcoID,           &
                       SpaceDim  = 2,               &
                       LevIDx    = -1,              &
                       OutUnit   = 'kg/m2/s',       &
                       WriteFreq = 'Manual',        &
                       AutoFill  = 1,               &
                       cID       = N,               & 
                       RC        = RC                ) 
    IF ( RC /= HCO_SUCCESS ) RETURN


    !---------------------------------------------
    ! %%%%% ANTHROPOGENIC HG (Hg0, Hg2, HgP) %%%%%
    !---------------------------------------------
 
    Cat = CATEGORY_ANTHRO

    ! Create diagnostic container
    DiagnName = 'HG0_ANTHRO'
    CALL Diagn_Create( am_I_Root,                   & 
                       HcoState,                    &
                       cName     = TRIM(DiagnName), &
                       ExtNr     = ExtNr,           &
                       Cat       = Cat,             &
                       Hier      = -1,              &
                       HcoID     = HcoID,           &
                       SpaceDim  = 2,               &
                       LevIDx    = -1,              &
                       OutUnit   = 'kg/m2/s',       &
                       WriteFreq = 'Manual',        &
                       AutoFill  = 1,               &
                       cID       = N,               & 
                       RC        = RC                ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Hg2
    HcoID = GetHemcoId( 'Hg2', HcoState, LOC, RC, ERR=.FALSE. )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create diagnostic container
    IF ( HcoID > 0 ) THEN
       DiagnName = 'HG2_ANTHRO'
       CALL Diagn_Create( am_I_Root,                   & 
                          HcoState,                    &
                          cName     = TRIM(DiagnName), &
                          ExtNr     = ExtNr,           &
                          Cat       = Cat,             &
                          Hier      = -1,              &
                          HcoID     = HcoID,           &
                          SpaceDim  = 2,               &
                          LevIDx    = -1,              &
                          OutUnit   = 'kg/m2/s',       &
                          WriteFreq = 'Manual',        &
                          AutoFill  = 1,               &
                          cID       = N,               & 
                          RC        = RC                ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! HgP
    HcoID = GetHemcoId( 'HgP', HcoState, LOC, RC, ERR=.FALSE. )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create diagnostic container
    IF ( HcoID > 0 ) THEN
       DiagnName = 'HGP_ANTHRO'
       CALL Diagn_Create( am_I_Root,                   & 
                          HcoState,                    &
                          cName     = TRIM(DiagnName), &
                          ExtNr     = ExtNr,           &
                          Cat       = Cat,             &
                          Hier      = -1,              &
                          HcoID     = HcoID,           &
                          SpaceDim  = 2,               &
                          LevIDx    = -1,              &
                          OutUnit   = 'kg/m2/s',       &
                          WriteFreq = 'Manual',        &
                          AutoFill  = 1,               &
                          cID       = N,               & 
                          RC        = RC                ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    !-------------------------------------------
    ! %%%%% BIOMASS BURNING HG %%%%%
    ! ==> defined in Diagn_Biomass
    !-------------------------------------------

  END SUBROUTINE Diagn_Hg
!EOC
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
    USE HCO_State_Mod, ONLY : HCO_GetHcoID
    USE HCO_State_Mod, ONLY : HCO_State
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
