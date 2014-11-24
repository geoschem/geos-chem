!BOM
#if defined(ESMF_) 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoi_esmf_mod.F90 
!
! !DESCRIPTION: Module HCOI\_ESMF\_MOD.F90 contains routines to
! interface HEMCO and ESMF. This code is run "outside" of HEMCO and always
! w/in the ESMF interface, so we can use the ESMF/MAPL specific error syntax.
! \\
! !INTERFACE:
!
      MODULE HCOI_ESMF_MOD
!
! !USES:
!
      USE ESMF
      USE MAPL_Mod

      IMPLICIT NONE

#     include "MAPL_Generic.h"
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: HCOI_ESMF_SetServices
      PUBLIC :: HCOI_ESMF_INIT
      PUBLIC :: HCOI_ESMF_RUN
      PUBLIC :: HCOI_ESMF_FINAL
!
! !PRIVATE MODULE FUNCTIONS:
!
      PRIVATE :: Hemco2FluxOut
!
! !PRIVATE MODULE VARIABLES:
!
      REAL*8, POINTER                 :: XMD   (:,:,:) => NULL()
      REAL*8, POINTER                 :: XDG   (:,:,:) => NULL()
      REAL*8, POINTER                 :: YMD   (:,:,:) => NULL()
      REAL*8, POINTER                 :: YDG   (:,:,:) => NULL()
      REAL*8, POINTER                 :: YSN   (:,:,:) => NULL()
      REAL*8, POINTER                 :: AM2   (:,:,:) => NULL()
      REAL*8, POINTER                 :: BXHGHT(:,:,:) => NULL()
      LOGICAL                         :: OnSimGrid
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller: Initial version 
!
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCOI_ESMF_SetServices
!
! !DESCRIPTION: Subroutine HCOI\_ESMF\_SetServices registers all required HEMCO 
! data so that it can be imported through the ESMF import state. 
! This routine is called at the beginning of a simulation - even ahead of 
! the initialization routines. Since this routine is called from outside of 
! the HEMCO environment, use the MAPL specific error codes! 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOI_ESMF_SetServices( am_I_Root, GC, ConfigFile, RC ) 
!
! !USES:
!
      USE HCO_ERROR_MOD
      USE HCO_CONFIG_MOD,  ONLY : Read2Buffer
      USE HCO_CONFIG_MOD,  ONLY : HcoCfgLin, GetNextLine
!
! !ARGUMENTS:
!
      LOGICAL,             INTENT(IN   )   :: am_I_Root
      TYPE(ESMF_GridComp), INTENT(INOUT)   :: GC
      CHARACTER(LEN=*),    INTENT(IN   )   :: ConfigFile
      INTEGER,             INTENT(  OUT)   :: RC
!
! !REVISION HISTORY:
!  29 Aug 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                     :: FLAG
      TYPE(HcoCfgLin), POINTER    :: ThisLine => NULL()

      ! ================================================================
      ! HCOI_ESMF_SetServices begins here
      ! ================================================================

      ! For MAPL/ESMF error handling (defined Iam and STATUS)
      __Iam__('HCOI_ESMF_SetServices (hcoi_esmf_mod.F90)') 

      ! ---------------------------------------------------------------------
      ! Read file into buffer
      ! ---------------------------------------------------------------------

      CALL Read2Buffer( am_I_Root, TRIM(ConfigFile), STATUS )
      IF ( STATUS /= HCO_SUCCESS ) THEN
         ASSERT_(.FALSE.)
      ENDIF

      ! Loop over all lines and set services according to inputi file content
      CALL GetNextLine ( ThisLine, FLAG )
      DO WHILE ( FLAG == HCO_SUCCESS ) 

         ! Skip data that is not read from external file but direct from
         ! the HEMCO input file
         IF ( TRIM(ThisLine%srcFile) == '-' ) THEN 
            CALL GetNextLine ( ThisLine, FLAG )
            CYCLE
         ENDIF

         ! 2D data
         IF ( TRIM(ThisLine%srcDim) == 'xy' ) THEN

            CALL MAPL_AddImportSpec(GC,            &
               SHORT_NAME = TRIM(ThisLine%cName),   &
               LONG_NAME  = TRIM(ThisLine%cName),   &
               UNITS      = TRIM(ThisLine%srcUnit), &
               DIMS       = MAPL_DimsHorzOnly,      &
               VLOCATION  = MAPL_VLocationNone,     &
               RC         = STATUS                   )
            VERIFY_(STATUS)

         ! 3D data: Assume central location in vertical dimension!
         ELSEIF ( TRIM(ThisLine%srcDim) == 'xyz' ) THEN
 
            CALL MAPL_AddImportSpec(GC,            &
               SHORT_NAME = TRIM(ThisLine%cName),   &
               LONG_NAME  = TRIM(ThisLine%cName),   &
               UNITS      = TRIM(ThisLine%srcUnit), &
               DIMS       = MAPL_DimsHorzVert,      &
               VLOCATION  = MAPL_VLocationCenter,   &
               RC         = STATUS                   )
            VERIFY_(STATUS)

         ! Return w/ error if not 2D or 3D data 
         ELSE
            ASSERT_(.FALSE.) 
         ENDIF

         ! Advance to next line
         CALL GetNextLine ( ThisLine, FLAG ) 

      ENDDO

      ! Free pointer
      ThisLine => NULL()

      ! Return success
      RETURN_(ESMF_SUCCESS)

      END SUBROUTINE HCOI_ESMF_SetServices
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_ESMF_INIT
!
! !DESCRIPTION: Subroutine HCOI\_ESMF\_INIT initializes the HEMCO derived
! types and arrays within the ESMF environment. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOI_ESMF_INIT( am_I_Root, Input_Opt,    &
                                 State_Met, State_Chm,    & 
                                 IMPORT,    ConfigFile, RC ) 
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod,    ONLY : OptInput
      USE GIGC_State_Met_Mod,    ONLY : MetState
      USE GIGC_State_Chm_Mod,    ONLY : ChmState
      USE GRID_MOD,              ONLY : XMID,  YMID
      USE GRID_MOD,              ONLY : XEDGE, YEDGE, YSIN
      USE GRID_MOD,              ONLY : AREA_M2
      USE TRACERID_MOD,          ONLY : IDEMIS, CTRMB
      USE CMN_SIZE_MOD,          ONLY : IIPAR, JJPAR, LLPAR
      USE TIME_MOD,              ONLY : GET_TS_EMIS, GET_TS_DYN

      ! HEMCO core routines
      USE HCO_ERROR_MOD
      USE HCO_MAIN_MOD,        ONLY : HCO_INIT
!
! !INPUT/OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
      TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input opts
      TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state 
      TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chem state 
      TYPE(ESMF_State), INTENT(IN   )  :: IMPORT     ! ESMF Import obj.
      CHARACTER(LEN=*), INTENT(IN   )  :: ConfigFile ! HEMCO config file
      INTEGER,          INTENT(  OUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8, ALLOCATABLE             :: SpecMW(:)
      REAL*8, ALLOCATABLE             :: EmSpecMW(:)
      REAL*8, ALLOCATABLE             :: MOLEC_RATIO(:)
      REAL*8                          :: TS_EMIS, TS_DYN
      INTEGER, ALLOCATABLE            :: SpecIDs(:)
      CHARACTER(LEN=255), ALLOCATABLE :: SpecNames(:)
      INTEGER                         :: I, N, tmpID
      INTEGER                         :: nLon, nLat, nLev
      INTEGER                         :: JSP, JNP 
      INTEGER                         :: HMRC

      !=================================================================
      ! HCOI_ESMF_INIT begins here!
      !=================================================================

      ! Assume success until otherwise
      __Iam__('HCOI_ESMF_INIT (HCOI_ESMF_MOD.F90)')

      !-----------------------------------------------------------------
      ! Molecules of tracer per molecules of species
      ! For now, pass default value of 1 and calculate the correct value 
      ! from CTRMB of tracerid_mod.F when used (hemco_dataread_mod).
      !-----------------------------------------------------------------

      ! Number of defined species
      N = Input_Opt%N_TRACERS

      ! Allocate temporary arrays
      ALLOCATE ( Molec_Ratio(N))
      ALLOCATE ( SpecMW(N))
      ALLOCATE ( EmSpecMW(N))
      ALLOCATE ( SpecIDs(N))
      ALLOCATE ( SpecNames(N))

      ! Set species IDs for more efficient identification 
      SpecIDs(1:N) = State_Chm%Trac_ID(1:N)

      ! Set species names
      !SpecNames(1:N) = Input_Opt%TRACER_FULLNAME(1:N)
      SpecNames(1:N) = State_Chm%Trac_Name(1:N)

      ! testing only
      if ( am_I_Root ) then
         write(6,*) 'Init HEMCO tracers:'
         write(6,*) 'N tracers: ', N
         write(6,*) 'IDs: ', SpecIDs
         write(6,*) 'SpecNames: ', SpecNames
      endif

!***********************************************************************
      ! Molecular weights of species
      ! For now, set to tracer molecular weight!! This is only correct
      ! for species where the emitted species is equal to the GC
      ! species, but not for VOCs! 
      ! As long as we use this (wrong) MW, the emission fields MUST be
      ! in units of mass carbon and NOT mass species, otherwise these
      ! emissions won't be properly calculated!!
      SpecMW(1:N) = Input_Opt%Tracer_MW_G(1:N)
!***********************************************************************

      ! Molecular weights of emitted species
      EmSpecMW(1:N) = Input_Opt%Tracer_MW_G(1:N)

      ! Emitted molecules per molecule of species. 
      DO I = 1, Input_Opt%N_TRACERS
         tmpID    = Input_Opt%ID_EMITTED(I) 
         IF ( tmpID <= 0 ) THEN
            Molec_Ratio(I) = 1.0d0
         ELSE 
            Molec_Ratio(I) = Input_Opt%TRACER_COEFF(I,tmpID)
         ENDIF
      ENDDO

      !-----------------------------------------------------------------
      ! Set grid settings
      !-----------------------------------------------------------------

      ! Always use simulation grid in ESMF! 

         ! Set grid dimensions
         nLon = IIPAR
         nLat = JJPAR
         nLev = LLPAR

         ! Point to all grid variables
         XMD    => XMID
         YMD    => YMID
         XDG    => XEDGE
         YDG    => YEDGE
         YSN    => YSIN
         AM2    => AREA_M2
         BXHGHT => State_Met%BXHEIGHT

         ! Set flag
         OnSimGrid = .TRUE.

      ! Get emission timestep in seconds
      TS_EMIS = GET_TS_EMIS() * 60d0
      TS_DYN  = GET_TS_DYN() * 60d0

      ! Set return code flag to HEMCO success. This value should be
      ! preserved throughout all HEMCO calls, otherwise an error
      ! will be returned!
      HMRC = HCO_SUCCESS

      !-----------------------------------------------------------------
      ! Initialize HEMCO
      !-----------------------------------------------------------------
      CALL HCO_INIT (   am_I_Root  = am_I_Root,                 &
                        nSpecies   = N,                         &
                        SpecIDs    = SpecIDs,                   &
                        SpecNames  = SpecNames,                 &
                        SpecMW     = SpecMW,                    &
                        EmSpecMW   = EmSpecMW,                  &
                        MolecRatio = Molec_Ratio,               &
                        nLon       = nLon,                      & 
                        nLat       = nLat,                      &
                        nLev       = nLev,                      &
                        XMID       = XMD,                       &
                        YMID       = YMD,                       &
                        XEDGE      = XDG,                       &
                        YEDGE      = YDG,                       &
                        YSIN       = YSN,                       &
                        AREA_M2    = AM2,                       &
                        BXHEIGHT_M = BXHGHT,                    & 
                        OnSimGrid  = OnSimGrid,                 &
                        ConfigFile = ConfigFile,                &
                        IMPORT     = IMPORT,                    &
                        TS_EMIS    = TS_EMIS,                   &
                        TS_DYN     = TS_DYN,                    &
                        RC         = HMRC                     )
      IF ( HMRC /= HCO_SUCCESS ) THEN
         ASSERT_(.FALSE.)
      ENDIF

      !-----------------------------------------------------------------
      ! Leave
      !-----------------------------------------------------------------
      IF ( ALLOCATED ( Molec_Ratio ) ) DEALLOCATE ( Molec_Ratio )
      IF ( ALLOCATED ( EmSpecMW    ) ) DEALLOCATE ( EmSpecMW    )
      IF ( ALLOCATED ( SpecMW      ) ) DEALLOCATE ( SpecMW      )
      IF ( ALLOCATED ( SpecIDs     ) ) DEALLOCATE ( SpecIDs     )
      IF ( ALLOCATED ( SpecNames   ) ) DEALLOCATE ( SpecNames   )

      IF ( OnSimGrid ) THEN
         XMD    => NULL()
         YMD    => NULL()
         XDG    => NULL()
         YDG    => NULL()
         YSN    => NULL()
         AM2    => NULL()
         BXHGHT => NULL()
      ENDIF

      ! Return w/ (ESMF) success
      RETURN_(ESMF_SUCCESS)

      END SUBROUTINE HCOI_ESMF_INIT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_ESMF_RUN
!
! !DESCRIPTION: Subroutine HCOI\_ESMF\_RUN runs HEMCO within ESMF. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOI_ESMF_RUN( am_I_Root, FluxOut, RC ) 
!
! !USES:
!
      USE HCO_ERROR_MOD
      USE HCO_MAIN_MOD,        ONLY : HCO_RUN, HCO_GetState
      USE HCO_TYPE_MOD,        ONLY : HCO_State, ResetArrays 

      USE TIME_MOD,              ONLY : GET_YEAR, GET_MONTH,  GET_DAY
      USE TIME_MOD,              ONLY : GET_HOUR, GET_MINUTE, GET_SECOND
      USE TIME_MOD,              ONLY : GET_DAY_OF_YEAR, GET_DAY_OF_WEEK
      USE CMN_O3_MOD,            ONLY : FSCALYR
!
! !INPUT/OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root        ! root CPU?
      REAL*8,           POINTER        :: FluxOut(:,:,:,:) ! output array
      INTEGER,          INTENT(  OUT)  :: RC               ! Failure or success
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE(HCO_State), POINTER  :: HcoState => NULL()
      INTEGER                   :: HMRC

      !=================================================================
      ! HCOI_ESMF_RUN begins here!
      !=================================================================

      ! Enter
      __Iam__('HCOI_ESMF_RUN (HCOI_ESMF_MOD.F90)')

      ! Get HEMCO state object (from hemco_main)
      CALL HCO_GetState ( HcoState )

      ! ================================================================
      ! Reset all emission and deposition values
      ! ================================================================
      CALL ResetArrays ( HcoState, HMRC )
      IF ( HMRC /= HCO_SUCCESS ) THEN
         ASSERT_(.FALSE.)
      ENDIF

      !=================================================================
      ! Set Hemco options and define all arrays needed by the core 
      ! module and extensions 
      !=================================================================

      ! Set range of tracers and emissions categories. 
      ! Set extension number to 0 (indicating that the core module 
      ! is executed).
      HcoState%TrcMin = 1 
      HcoState%TrcMax = HcoState%nSpecies 
      HcoState%CatMin = 1 
      HcoState%CatMax = -1 
      HcoState%ExtNr  = 0 

      ! Current simulation time 
      IF ( FSCALYR < 0 ) THEN
         HcoState%sYear = GET_YEAR()
      ELSE
         HcoState%sYear = FSCALYR
      ENDIF
      HcoState%sMonth     = GET_MONTH()
      HcoState%sDay       = GET_DAY()
      HcoState%sHour      = GET_HOUR()
      HcoState%sMin       = GET_MINUTE()
      HcoState%sSec       = GET_SECOND()
      HcoState%sDayOfYear = GET_DAY_OF_YEAR()
      HcoState%sWeekDAy   = GET_DAY_OF_WEEK() !0=Sun,...,6=Sat

      ! Use temporary array?
      HcoState%FillTemp3D = .FALSE. 

      ! Set return code flag to HEMCO success. This value should be
      ! preserved throughout all HEMCO calls, otherwise an error
      ! will be returned!
      HMRC = HCO_SUCCESS

      ! ================================================================
      ! Run HEMCO core module 
      ! ================================================================
      CALL HCO_RUN ( am_I_Root, HcoState, HMRC )
      IF ( HMRC /= HCO_SUCCESS ) THEN
         ASSERT_(.FALSE.)
      ENDIF

      ! ================================================================
      ! Run HEMCO extensions
      ! ================================================================

      ! ... no extensions at the moment ...

      ! ================================================================
      ! Translate emissions array from HEMCO state onto output array
      ! Values are retained in [kg/m2/s] 
      ! ================================================================

      ! Pass emissions from HEMCO state object to output array. Reset
      ! output emissions first!
      FluxOut = 0d0
      CALL Hemco2FluxOut ( HcoState, FluxOut, HMRC )
      IF ( HMRC /= HCO_SUCCESS ) THEN
         ASSERT_(.FALSE.)
      ENDIF

      ! Free pointers
      HcoState => NULL()

      ! Return w/ ESMF success
      RETURN_(ESMF_SUCCESS)

      END SUBROUTINE HCOI_ESMF_RUN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Hemco2FluxOut
!
! !DESCRIPTION: Function Hemco2FluxOut fills the passed output array 
! with the corresponding emission values of the HEMCO state object. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Hemco2FluxOut ( HcoState, FluxOut, HMRC )
!
! !USES:
!
      USE HCO_ERROR_MOD 
      USE HCO_TYPE_MOD,       ONLY : HCO_State
!
! !ARGUMENTS:
!
      TYPE(HCO_State), POINTER        :: HcoState         ! Hemco options
      REAL*8,          POINTER        :: FluxOut(:,:,:,:) ! Output array 
      INTEGER,         INTENT(INOUT)  :: HMRC          ! Failure?
!
! !REVISION HISTORY:
!  01 May 2012 - C. Keller - Initial Version
!  20 Aug 2013 - C. Keller - Now pass from HEMCO to chemistry state
!  12 Sep 2013 - C. Keller - Now pass FluxOut instead of chemistry state
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
      INTEGER                  :: N, nTrc, thisTrc 
      CHARACTER(LEN=255)       :: MSG, LOC

      !=================================================================
      ! Hemco2FluxOut begins here
      !=================================================================

      ! For (HEMCO) error handling
      LOC = 'Hemco2FluxOut (HCOI_ESMF_MOD.F90)'            
  
      ! Get number of specified tracers 
      nTrc = SIZE( FluxOut, 4 ) 

      ! Loop over all specified tracers 
      DO N = 1, nTrc

         ! if emissions are defined, copy array from Emsr3D into
         ! output array  
        IF ( ASSOCIATED(HcoState%Emsr3D(N)%Arr3D) ) THEN
            IF ( HcoState%OnSimGrid ) THEN
               FluxOut(:,:,:,N) = HcoState%Emsr3D(N)%Arr3D(:,:,:)
            ELSE
               MSG = 'Emission and simulation grid do not match!'
               CALL HCO_WARNING ( MSG, LOC, HMRC )
               RETURN
            ENDIF
         ENDIF
      ENDDO !N 

      ! Return w/ success
      HMRC = HCO_SUCCESS

      END SUBROUTINE Hemco2FluxOut 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_ESMF_FINAL
!
! !DESCRIPTION: Subroutine HCOI\_ESMF\_FINAL cleans up HEMCO. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOI_ESMF_FINAL
!
! !USES:
!
      USE HCO_MAIN_MOD,        ONLY : HCO_FINAL
!
! !INPUT/OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! HCOI_ESMF_FINAL begins here!
      !=================================================================

      ! Cleanup HEMCO
      CALL HCO_FINAL

      END SUBROUTINE HCOI_ESMF_FINAL
!EOC
      END MODULE HCOI_ESMF_MOD
#endif
!EOM
