!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_gfed3_mod.F90
!
! !DESCRIPTION: Module HCOX\_GFED3\_MOD contains routines to calculate
! GFED3 biomass burning emissions in HEMCO. 
!
! !NOTES:
!
! !REFERENCES:
!
! !INTERFACE: 
!
MODULE HCOX_GFED3_MOD
!
! !USES:
! 
  USE HCO_ERROR_MOD
  USE HCO_DIAGN_MOD
  USE HCO_STATE_MOD,  ONLY : HCO_State
  USE HCOX_State_MOD, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_GFED3_Init
  PUBLIC :: HCOX_GFED3_Run
  PUBLIC :: HCOX_GFED3_Final
!
! !REMARKS:
!  Monthly emissions of DM are read from disk,
!  multiplied by daily and 3hourly fractions (if necessary), and then
!  multiplied by the appropriate emission factors to produce biomass
!  burning emissions.
!                                                                             .
!  References:
!  ============================================================================
!  (1 ) Original GFED3 database from Guido van der Werf 
!        http://www.falw.vu/~gwerf/GFED/GFED3/emissions/
!  (2 ) Giglio, L., Randerson, J. T., van der Werf, G. R., Kasibhatla, P. S.,
!       Collatz, G. J., Morton, D. C., and DeFries, R. S.: Assessing
!       variability and long-term trends in burned area by merging multiple 
!       satellite fire products, Biogeosciences, 7, 1171-1186, 
!       doi:10.5194/bg-7-1171-2010, 2010.
!  (3 ) van der Werf, G. R., Randerson, J. T., Giglio, L., Collatz, G. J.,
!       Mu, M., Kasibhatla, P. S., Morton, D. C., DeFries, R. S., Jin, Y., 
!       and van Leeuwen, T. T.: Global fire emissions and the contribution of 
!       deforestation, savanna, forest, agricultural, and peat fires 
!       (1997â~@~S2009), Atmos. Chem. Phys., 10, 11707-11735, 
!       doi:10.5194/acp-10-11707-2010, 2010.
!
! !REVISION HISTORY: 
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  07 Sep 2011 - R. Yantosca   - Added ProTeX headers 
!  14 Feb 2012 - M. Payer      - Add modifications for CH4 (K. Wecht)
!  01 Mar 2012 - R. Yantosca   - Now reference new grid_mod.F90
!  06 Mar 2012 - P. Kasibhatla - Final version
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  15 Dec 2013 - C. Keller   - Now a HEMCO extension. Emissions in kg/m2/s,
!                              emission factors in kg/kgDM.
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  08 Aug 2014 - R. Yantosca - Now avoid ASCII file reads for ESMF
!EOP
!------------------------------------------------------------------------------
!
! !DEFINED PARAMETERS:
!
  !=================================================================
  ! MODULE PARAMETERS
  !
  ! N_EMFAC : Number of emission factors per species
  ! N_SPEC  : Max. number of species
  !=================================================================
  INTEGER,           PARAMETER :: N_EMFAC = 6
  INTEGER,           PARAMETER :: N_SPEC  = 25
!
! !PRIVATE TYPES:
!
  !=================================================================
  ! HEMCO VARIABLES 
  !
  ! ExtNr   : Extension number 
  ! DoDay   : TRUE if dialy scale factors are used 
  ! Do3Hr   : TRUE if 3-hourly scale factors are used 
  !=================================================================
  INTEGER                       :: ExtNr
  LOGICAL                       :: DoDay
  LOGICAL                       :: Do3Hr

  !=================================================================
  ! SPECIES VARIABLES 
  !
  ! nSpc     : Number of GFED3 species (specified in config. file)
  ! SpcNames : Names of all used GFED3 species
  ! HcoIDs   : HEMCO species IDs of all used GFED3 species 
  ! gfedIDs  : Index of used GFED3 species in scale factor table 
  !=================================================================
  INTEGER                        :: nSpc
  CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
  INTEGER,           ALLOCATABLE :: HcoIDs(:)
  INTEGER,           ALLOCATABLE :: GfedIDs(:)

  !=================================================================
  ! SCALE FACTORS 
  !
  ! GFED3_EMFAC: emission scale factors for each species and 
  !              emission factor type. The filename of the emissions
  !              emissions factor table is specified in the HEMCO
  !              configuration file. All scale factors in kg/kgDM.
  ! COScale    : CO scale factor to account for porduction from 
  !              VOCs. Read from HEMCO configuration file.
  !=================================================================
  REAL(hp),          ALLOCATABLE :: GFED3_EMFAC(:,:)
  REAL(sp)                       :: COScale

  !=================================================================
  ! DATA ARRAY POINTERS 
  !
  ! These are the pointers to the 6 input data specified in the 
  ! the configuration file
  !=================================================================
  REAL(hp), POINTER   :: GFED3_WDL(:,:) => NULL()
  REAL(hp), POINTER   :: GFED3_SAV(:,:) => NULL()
  REAL(hp), POINTER   :: GFED3_PET(:,:) => NULL()
  REAL(hp), POINTER   :: GFED3_FOR(:,:) => NULL()
  REAL(hp), POINTER   :: GFED3_AGW(:,:) => NULL()
  REAL(hp), POINTER   :: GFED3_DEF(:,:) => NULL()
  REAL(hp), POINTER   :: HUMTROP  (:,:) => NULL()
  REAL(hp), POINTER   :: DAYSCAL  (:,:) => NULL()
  REAL(hp), POINTER   :: HRSCAL   (:,:) => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_GFED3_Run 
!
! !DESCRIPTION: Subroutine HcoX\_GFED3\_Run is the driver run routine to 
! calculate seasalt emissions in HEMCO.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_GFED3_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_EmisList_Mod, ONLY : EmisList_GetDataArr
    USE HCO_FluxArr_MOD,  ONLY : HCO_EmisAdd
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState   ! Output obj
    TYPE(Ext_State), POINTER        :: ExtState  ! Module options  
    INTEGER,         INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  15 Dec 2013 - C. Keller   - Now a HEMCO extension 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL             :: IsCO 
    LOGICAL, SAVE       :: FIRST = .TRUE.
    INTEGER             :: N, M
    REAL(hp), POINTER   :: Arr2D (:,:) => NULL()
    REAL(hp), POINTER   :: TMPPTR(:,:) => NULL()

    REAL(hp), TARGET    :: SpcArr(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET    :: TypArr(HcoState%NX,HcoState%NY)
   
    !=================================================================
    ! HCOX_GFED3_Run begins here!
    !=================================================================

    ! Return if extension disabled 
    IF ( .NOT. ExtState%GFED3 ) RETURN

    ! Enter 
    CALL HCO_ENTER ( 'HCOX_GFED3_Run (hcox_gfed3_mod.F90)', RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! Get pointers to data arrays 
    !-----------------------------------------------------------------
    IF ( FIRST ) THEN
       CALL EmisList_GetDataArr ( am_I_Root, 'GFED3_WDL', GFED3_WDL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL EmisList_GetDataArr ( am_I_Root, 'GFED3_SAV', GFED3_SAV, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL EmisList_GetDataArr ( am_I_Root, 'GFED3_PET', GFED3_PET, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL EmisList_GetDataArr ( am_I_Root, 'GFED3_FOR', GFED3_FOR, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL EmisList_GetDataArr ( am_I_Root, 'GFED3_AGW', GFED3_AGW, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL EmisList_GetDataArr ( am_I_Root, 'GFED3_DEF', GFED3_DEF, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       CALL EmisList_GetDataArr ( am_I_Root, 'HUMTROP', HUMTROP, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Also point to scale factors if needed
       IF ( DoDay ) THEN
          CALL EmisList_GetDataArr ( am_I_Root, 'GFED3_FRAC_DAY', &
                                     DAYSCAL,   RC               )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
       IF ( Do3Hr ) THEN
          CALL EmisList_GetDataArr ( am_I_Root, 'GFED3_FRAC_3HOUR', &
                                     HRSCAL,    RC                 )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       FIRST = .FALSE.
    ENDIF

    !-----------------------------------------------------------------
    ! Calculate emissions for defined species
    !-----------------------------------------------------------------

    DO N = 1, nSpc

       ! Continue if species not defined
       IF ( HcoIDs(N)  < 0 ) CYCLE
       IF ( GfedIDs(N) < 0 ) CYCLE

       ! SpcArr are the total biomass burning emissions for this
       ! species. TypArr are the emissions from a given source type. 
       SpcArr = 0.0_hp

       ! Is this CO?
       IF ( TRIM(SpcNames(N)) == "CO" ) THEN
          IsCO = .TRUE.
       ELSE
          IsCO = .FALSE.
       ENDIF

       ! Calculate emissions for all source types
       DO M = 1, N_EMFAC
         
          ! Point to the emission factor array for each source type
          SELECT CASE ( M ) 
             CASE( 1 )
                TMPPTR => GFED3_AGW
             CASE( 2 )
                TMPPTR => GFED3_DEF
             CASE( 3 )
                TMPPTR => GFED3_FOR
             CASE( 4 )
                TMPPTR => GFED3_PET
             CASE( 5 )
                TMPPTR => GFED3_SAV
             CASE( 6 )
                TMPPTR => GFED3_WDL
             CASE DEFAULT
                CALL HCO_ERROR ( 'Undefined emission factor', RC )
                RETURN
          END SELECT

          ! Calculate emissions for this type. The emission factors 
          ! per type are in kgDM/m2/s, and the GFED3_EMFAC scale factors
          ! are in kg/kgDM (or kgC/kgDM for VOCs). This gives us TypArr
          ! in kg/m2/s.
          TypArr = TmpPtr * GFED3_EMFAC(GfedIDs(N),M)

          ! Use woodland emission factors for 'deforestation' outside
          ! humid tropical forest
          IF ( M == 2 ) THEN
             WHERE ( HUMTROP == 0.0_hp ) 
                TypArr = TmpPtr * GFED3_EMFAC(GfedIDs(N),6)
             END WHERE
          ENDIF

          ! Eventually add daily / 3-hourly scale factors. These scale
          ! factors are unitless.
          IF ( DoDay ) TypArr = TypArr * DAYSCAL
          IF ( Do3Hr ) TypArr = TypArr * HRSCAL

          ! For CO, multiply with CO scale factor
          IF ( IsCo ) TypArr = TypArr * COScale        

          ! Add to output array
          SpcArr = SpcArr + TypArr

          ! Nullify pointer
          TmpPtr  => NULL()

       ENDDO !M

       ! Add flux to HEMCO emission array
       CALL HCO_EmisAdd( HcoState, SpcArr, HcoIDs(N), RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Eventually update diagnostics
       IF ( Diagn_AutoFillLevelDefined(2) ) THEN
          Arr2D => SpcArr
          CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr, &
                             Cat=-1, Hier=-1, HcoID=HcoIDs(N), &
                             AutoFill=1, Array2D=Arr2D, RC=RC   )
          IF ( RC /= HCO_SUCCESS ) RETURN
          Arr2D => NULL()
       ENDIF
    ENDDO !N

    ! Nullify pointers for safety's sake
    TmpPtr  => NULL()
    Arr2D   => NULL()

    ! Leave w/ success
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE HCOX_GFED3_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_GFED3_Init
!
! !DESCRIPTION: Subroutine HcoX\_GFED3\_Init initializes all
!  extension variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_GFED3_Init ( am_I_Root, HcoState, ExtName, &
                               ExtState,  RC                  ) 
!
! !USES:
!
    USE HCO_STATE_MOD,          ONLY : HCO_GetHcoID
    USE HCO_STATE_MOD,          ONLY : HCO_GetExtHcoID
    USE HCO_ExtList_Mod,        ONLY : GetExtNr, GetExtOpt
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root   ! root CPU?
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName     ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState    ! Options object
!                                                   
! !INPUT/OUTPUT PARAMETERS:                         
!                                                   
    TYPE(HCO_State),  POINTER        :: HcoState    ! HEMCO state object 
    INTEGER,          INTENT(INOUT)  :: RC          ! Return status
!
! !REVISION HISTORY:
!  15 Dec 2013 - C. Keller   - Initial version
!  08 Aug 2014 - R. Yantosca - Now include hcox_gfed3_include.H, which defines
!                              GFED3_SPEC_NAME and GFED3_EMFAC arrays
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: MSG, ScalFile
    CHARACTER(LEN=255) :: GFED3_SPEC_NAME(N_SPEC)
    INTEGER            :: tmpNr, AS, IU_FILE, IOS
    INTEGER            :: N, M, NDUM
    CHARACTER(LEN=31)  :: tmpName
    LOGICAL            :: Matched

    !=================================================================
    ! HCOX_GFED3_Init begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter 
    CALL HCO_ENTER ( 'HCOX_GFED3_Init (hcox_gfed3_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ---------------------------------------------------------------------- 
    ! Get settings 
    ! ---------------------------------------------------------------------- 
 
    ! Get CO scale factor 
    CALL GetExtOpt ( ExtNr, 'CO scale factor', &
                     OptValSp=COScale, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Use daily scale factors?
    tmpName = TRIM(ExtName) // "_monthly"
    tmpNr   = GetExtNr( TRIM(tmpName) )
    IF ( tmpNr > 0 ) THEN
       DoDay = .TRUE.
    ELSE
       DoDay = .FALSE.
    ENDIF

    ! Use 3-hourly scale factors?
    tmpName = TRIM(ExtName) // "_3hourly"
    tmpNr   = GetExtNr( TRIM(tmpName) )
    IF ( tmpNr > 0 ) THEN
       Do3Hr = .TRUE.
    ELSE
       Do3Hr = .FALSE.
    ENDIF

    !----------------------------------------------------------------------- 
    ! Initialize GFED3 scale factors
    !----------------------------------------------------------------------- 

    ! Allocate scale factors table
    ALLOCATE ( GFED3_EMFAC ( N_SPEC, N_EMFAC ), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( 'Cannot allocate GFED3_EMFAC', RC )
       RETURN
    ENDIF
    GFED3_EMFAC = 0.0_hp

    ! Now get definitions for GFED3_EMFAC and GFED3_SPEC_NAME from an include 
    ! file.  This avoids ASCII file reads in the ESMF environment.  To update
    ! the emission factors, one just needs to modify the include file.
    ! This can be done with the script HEMCO/Extensions/Preprocess/gfed3.pl,
    ! (bmy, 8/14/14)
#include "hcox_gfed3_include.H"

    !----------------------------------------------------------------------- 
    ! Match specified species with GFED3 species
    ! The species to be used are specified in the HEMCO configuration file.
    ! Match these species with the ones found in the scale factors table.
    !----------------------------------------------------------------------- 

    ! Prompt to log file
    MSG = 'Use GFED3 extension'
    CALL HCO_MSG( MSG, SEP1='-' )
    WRITE(MSG,*) '   - Use daily scale factors : ', DoDay 
    CALL HCO_MSG( MSG )
    WRITE(MSG,*) '   - Use hourly scale factors: ', Do3Hr
    CALL HCO_MSG( MSG )
    WRITE(MSG,*) '   - CO scale factor         : ', COScale
    CALL HCO_MSG( MSG )

    ! Get HEMCO species IDs of all species specified in configuration file
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( nSpc == 0 ) THEN
       MSG = 'No GFED3 species specified'
       CALL HCO_ERROR ( MSG, RC ) 
       RETURN
    ENDIF

    ! GFEDIDS are the matching indeces of the HEMCO species in GFED3_EMFAC.
    ALLOCATE ( GfedIDs(nSpc), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( 'Cannot allocate GfedIDs', RC )
       RETURN
    ENDIF
    GfedIDs = -1

    ! Find matching GFED3 index for each specified species
    DO N = 1, nSpc
       IF ( HcoIDs(N) < 0 ) CYCLE
       Matched = .FALSE.
       DO M = 1, N_SPEC 
          IF ( TRIM(SpcNames(N)) == TRIM(GFED3_SPEC_NAME(M)) ) THEN
             GfedIDs(N) = M
             Matched    = .TRUE.

             MSG = '   - Use GFED3 species ' // TRIM(SpcNames(N))
             CALL HCO_MSG( MSG )
             EXIT ! go to next species
          ENDIF
       ENDDO
       IF ( .NOT. Matched ) THEN
          MSG = 'Species '// TRIM(SpcNames(N)) //' not found in GFED3'
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF
    ENDDO

    ! Enable module
    ExtState%GFED3 = .TRUE.

    ! Return w/ success
    CALL HCO_LEAVE ( RC ) 
 
  END SUBROUTINE HCOX_GFED3_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_GFED3_Final 
!
! !DESCRIPTION: Subroutine HcoX\_GFED3\_Final deallocates 
!  all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_GFED3_Final
!
! !REVISION HISTORY:
!  15 Dec 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! HCOX_GFED3_Final begins here!
    !=================================================================

    ! Free pointers
    GFED3_WDL => NULL()
    GFED3_SAV => NULL()
    GFED3_PET => NULL()
    GFED3_FOR => NULL()
    GFED3_AGW => NULL()
    GFED3_DEF => NULL()
    HUMTROP   => NULL()
    DAYSCAL   => NULL()
    HRSCAL    => NULL()

    ! Cleanup module arrays
    IF ( ALLOCATED( GFED3_EMFAC ) ) DEALLOCATE( GFED3_EMFAC )
    IF ( ALLOCATED( GfedIDs     ) ) DEALLOCATE( GfedIds     )
    IF ( ALLOCATED( HcoIDs      ) ) DEALLOCATE( HcoIDs      )
    IF ( ALLOCATED( SpcNames    ) ) DEALLOCATE( SpcNames    )

  END SUBROUTINE HCOX_GFED3_Final
!EOC
END MODULE HCOX_GFED3_MOD
