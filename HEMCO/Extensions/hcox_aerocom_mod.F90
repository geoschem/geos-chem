!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_aerocom_mod.F90
!
! !DESCRIPTION: Module HCOX\_AeroCom\_Mod.F90 is a HEMCO extension to use
! AeroCom volcano emissions from ascii tables. This module reads the daily
! AeroCom tables and emits the emissions according to the information in 
! this file.
!\\
!\\
! Each AeroCom table is expected to list the volcano location, sulfur 
! emissions (in kg S/s), and the volcano elevation as well as the 
! volcano plume column height. These entries need be separated by space
! characters. For example:
!###  LAT (-90,90), LON (-180,180), SULFUR [kg S/s], ELEVATION [m], CLOUD_COLUMN_HEIGHT [m]
!### If elevation=cloud_column_height, emit in layer of elevation
!### else, emit in top 1/3 of cloud_column_height
!volcano::
!50.170 6.850 3.587963e-03 600. 600. 
!::
!\\
!\\
! The sulfur read from table is emitted as the species defined in the
! AeroCom settings section. More than one species can be provided. Mass
! sulfur is automatically converted to mass of emitted species (using the
! emitted molecular weight and molecular ratio of the corresponding HEMCO
! species). Additional scale factors can be defined in the settings section
! by using the (optional) setting 'Scaling_<SpecName>'.
! For example, to emit SO2 and BrO from volcanoes, with an additional scale
! factor of 1e-4 kg BrO / kgS for BrO, use the following setting:
!115     AeroCom_Volcano   : on    SO2/BrO
!    --> Scaling_BrO       :       1.0e-4
!    --> AeroCom_Table     :       /path/to/Aerocom.so2_volcanic.$YYYY$MM$DD.rc
!\\
!\\
! !REMARKS:
! This extension was primarily added for usage within GEOS-5. 
! When using this extension, you should turn off any other volcano emission
! inventories!
!\\
!\\
! !INTERFACE:
!
MODULE HCOX_AeroCom_Mod 
!
! !USES:
!
  USE HCO_Error_MOD
  USE HCO_Diagn_MOD
  USE HCOX_State_MOD, ONLY : Ext_State
  USE HCO_State_MOD,  ONLY : HCO_State 

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_AeroCom_Run
  PUBLIC :: HCOX_AeroCom_Init
  PUBLIC :: HCOX_AeroCom_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: ReadVolcTable
  PRIVATE :: EmitVolc
!
! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
  INTEGER                         :: ExtNr     = -1  ! Extension number
  INTEGER                         :: CatErupt  = -1  ! Category of eruptive emissions
  INTEGER                         :: CatDegas  = -1  ! Category of degassing emissions
  INTEGER                         :: nSpc      =  0  ! # of species
  INTEGER                         :: nVolc     =  0  ! # of volcanoes in buffer 
  INTEGER,  ALLOCATABLE           :: SpcIDs(:)       ! HEMCO species IDs
  REAL(sp), ALLOCATABLE           :: SpcScl(:)       ! Species scale factors
  REAL(sp), ALLOCATABLE           :: VolcSlf(:)      ! Sulface emissions [kg S/s]
  REAL(sp), ALLOCATABLE           :: VolcElv(:)      ! Elevation [m]
  REAL(sp), ALLOCATABLE           :: VolcCld(:)      ! Cloud column height [m]
  INTEGER,  ALLOCATABLE           :: VolcIdx(:)      ! Lon grid index 
  INTEGER,  ALLOCATABLE           :: VolcJdx(:)      ! Lat grid index 
  CHARACTER(LEN=255)              :: FileName        ! Volcano file name

  ! AeroCom data is in kgS. Will be converted to kg emitted species.
  ! MW_S is the molecular weight of sulfur 
  REAL(hp), PARAMETER             :: MW_S = 32.0_hp
CONTAINS
!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_AeroCom_Run
!
! !DESCRIPTION: Subroutine HCOX\_AeroCom\_Run is the driver routine 
! for the customizable HEMCO extension. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_AeroCom_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_FluxArr_Mod,  ONLY : HCO_EmisAdd
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root   ! Are we on the root CPU?
    TYPE(Ext_State), POINTER       :: ExtState    ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState    ! Hemco state 
    INTEGER,         INTENT(INOUT) :: RC          ! Success or failure
!
! !REMARKS:
!  
!
! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: N
    REAL(sp)            :: SO2degas(HcoState%NX,HcoState%NY,HcoState%NZ) ! degassing
    REAL(sp)            :: SO2erupt(HcoState%NX,HcoState%NY,HcoState%NZ) ! eruptive
    REAL(sp)            :: iFlx    (HcoState%NX,HcoState%NY,HcoState%NZ)
    LOGICAL             :: ERR
    CHARACTER(LEN=255)  :: MSG

    !=================================================================
    ! HCOX_AEROCOM_RUN begins here!
    !=================================================================

    ! Sanity check: return if extension not turned on
    IF ( .NOT. ExtState%AeroCom ) RETURN

    ! Enter
    CALL HCO_ENTER ( 'HCOX_AeroCom_Run (hcox_aerocom_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Read/update the volcano data (will be done only if this is a new
    ! day) 
    CALL ReadVolcTable ( am_I_Root, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN 

    ! Emit volcanos into SO2degas and SO2erupt
    CALL EmitVolc ( am_I_Root, HcoState, ExtState, SO2degas, SO2erupt, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Add eruptive and degassing emissions to emission arrays & diagnostics 
    DO N = 1, nSpc
       iFlx = SO2degas * SpcScl(N) 
       CALL HCO_EmisAdd( am_I_Root, HcoState,    iFlx, SpcIDs(N), &
                         RC,        ExtNr=ExtNr, Cat=CatDegas )
       IF ( RC /= HCO_SUCCESS ) RETURN

       iFlx = SO2erupt * SpcScl(N) 
       CALL HCO_EmisAdd( am_I_Root, HcoState,    iFlx, SpcIDs(N), &
                         RC,        ExtNr=ExtNr, Cat=CatErupt )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDDO !N

    ! Return w/ success
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE HCOX_AeroCom_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_AeroCom_Init
!
! !DESCRIPTION: Subroutine HCOX\_AeroCom\_Init initializes the HEMCO
! CUSTOM extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_AeroCom_Init( am_I_Root, HcoState, ExtName, &
                                ExtState,  RC                  )
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt
    USE HCO_STATE_MOD,      ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root
    CHARACTER(LEN=*), INTENT(IN   ) :: ExtName    ! Extension name
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options      
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState   ! Hemco state 
    INTEGER,          INTENT(INOUT) :: RC 

! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(sp)                       :: ValSp
    INTEGER                        :: N, Dum
    LOGICAL                        :: FOUND
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
    CHARACTER(LEN=255)             :: MSG 

    !=================================================================
    ! HCOX_AEROCOM_INIT begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER ( 'HCOX_AeroCom_Init (hcox_aerocom_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get species IDs. 
    CALL HCO_GetExtHcoID( HcoState, ExtNr, SpcIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! There must be at least one species
    IF ( nSpc == 0 ) THEN
       CALL HCO_ERROR ( 'No AeroCom species specified', RC )
       RETURN
    ENDIF

    ! Determine scale factor to be applied to each species. This is 1.00 
    ! by default, but can be set in the HEMCO configuration file via setting
    ! Scaling_<SpcName>. 
    ALLOCATE ( SpcScl(nSpc) ) 
    DO N = 1, nSpc
       CALL GetExtOpt( ExtNr, 'Scaling_'//TRIM(SpcNames(N)), OptValSp=ValSp, &
                       FOUND=FOUND, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( FOUND ) THEN
          SpcScl(N) = ValSp
       ELSE
          SpcScl(N) = 1.0_sp
       ENDIF

       ! Add conversion factor from kg S to kg of emitted species
       SpcScl(N) = SpcScl(N) * HcoState%Spc(SpcIDs(N))%EmMW_g &
                 * HcoState%Spc(SpcIDs(N))%MolecRatio / MW_S

    ENDDO

    ! Get location of volcano table. This must be provided.
    CALL GetExtOpt( ExtNr, 'AeroCom_Table', OptValChar=FileName, &
                    FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS .OR. .NOT. FOUND ) THEN
       MSG = 'Cannot properly read AeroCom table file name. Please provide ' // &
             'the AeroCom table as a setting to the AeroCom extension. The ' // &
             'name of this setting must be `AeroCom_Table`.'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF

    ! See if eruptive and degassing hierarchies are given
    CatErupt = 51
    CatDegas = 52
    CALL GetExtOpt( ExtNr, 'Cat_Degassing', OptValInt=Dum, &
                    FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( FOUND ) CatDegas = Dum
    CALL GetExtOpt( ExtNr, 'Cat_Eruptive', OptValInt=Dum, &
                    FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( FOUND ) CatErupt = Dum

    ! Verbose mode
    IF ( am_I_Root ) THEN
       MSG = 'Use emissions extension `AeroCom_Volcano`:'
       CALL HCO_MSG( MSG )

       MSG = ' - use the following species (Name, HcoID, Scaling relative to kgS):'
       CALL HCO_MSG(MSG)
       DO N = 1, nSpc
          WRITE(MSG,*) TRIM(SpcNames(N)), ', ', SpcIDs(N), ', ', SpcScl(N)
          CALL HCO_MSG(MSG)
       ENDDO
       WRITE(MSG,*) ' - Emit eruptive emissions as category ', CatErupt
       CALL HCO_MSG( MSG )
       WRITE(MSG,*) ' - Emit degassing emissions as category ', CatDegas
       CALL HCO_MSG( MSG )
    ENDIF

    ! Activate this extension
    ExtState%AeroCom = .TRUE.

    ! Cleanup 
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)

    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE HCOX_AeroCom_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_AeroCom_Final
!
! !DESCRIPTION: Subroutine HCOX\_AeroCom\_Final finalizes the HEMCO
!  AeroCom extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_AeroCom_Final
!
! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! HCOX_AEROCOM_FINAL begins here!
    !=================================================================
    IF ( ALLOCATED(VolcSlf  ) ) DEALLOCATE ( VolcSlf   )
    IF ( ALLOCATED(VolcElv  ) ) DEALLOCATE ( VolcElv   )
    IF ( ALLOCATED(VolcCld  ) ) DEALLOCATE ( VolcCld   )
    IF ( ALLOCATED(VolcIdx  ) ) DEALLOCATE ( VolcIdx   )
    IF ( ALLOCATED(VolcJdx  ) ) DEALLOCATE ( VolcJdx   )
    IF ( ALLOCATED(SpcIDs   ) ) DEALLOCATE ( SpcIDs    )
    IF ( ALLOCATED(SpcScl   ) ) DEALLOCATE ( SpcScl    )

  END SUBROUTINE HCOX_AeroCom_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadVolcTable 
!
! !DESCRIPTION: Subroutine ReadVolcTable reads the AeroCom volcano table of the
! current day. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadVolcTable ( am_I_Root, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_CharTools_Mod
    USE inquireMod,         ONLY : findfreeLun
    USE HCO_CLOCK_MOD,      ONLY : HcoClock_NewDay
    USE HCO_CLOCK_MOD,      ONLY : HcoClock_Get
    USE HCO_GeoTools_MOD,   ONLY : HCO_GetHorzIJIndex 
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options      
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState   ! Hemco state 
    INTEGER,          INTENT(INOUT) :: RC 

! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: YYYY, MM, DD
    INTEGER               :: N, LUN, IOS, AS
    INTEGER               :: nCol
    REAL(sp)              :: Dum(10)
    REAL(hp), ALLOCATABLE :: VolcLon(:)      ! Volcano longitude [deg E]
    REAL(hp), ALLOCATABLE :: VolcLat(:)      ! Volcano latitude  [deg N]
    LOGICAL               :: FileExist, EOF
    CHARACTER(LEN=255)    :: ThisFile, ThisLine
    CHARACTER(LEN=255)    :: MSG 
    CHARACTER(LEN=255)    :: LOC = 'ReadVolcTable (hcox_aerocom_mod.F90)' 

    !=================================================================
    ! ReadVolcTable begins here!
    !=================================================================

    ! Do only if it's a new day...
    IF ( HcoClock_NewDay( EmisTime=.TRUE. ) ) THEN

       ! Get current year, month, day
       CALL HcoClock_Get ( cYYYY=YYYY, cMM=MM, cDD=DD, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Get file name
       ThisFile = FileName
       CALL HCO_CharParse( ThisFile, YYYY, MM, DD, 0, 0, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
 
       ! Verbose
       IF ( am_I_Root ) THEN
          MSG = 'AeroCom: reading ' // TRIM(ThisFile)
          CALL HCO_MSG(MSG)
       ENDIF
 
       ! Check if file exists
       INQUIRE ( FILE=TRIM(ThisFile), EXIST=FileExist )
       IF ( .NOT. FileExist ) THEN
          MSG = 'Cannot find ' // TRIM(ThisFile)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Open file
       LUN = findFreeLun()
       OPEN ( LUN, FILE=TRIM(ThisFile), STATUS='OLD', IOSTAT=IOS )
       IF ( IOS /= 0 ) THEN
          MSG = 'Error reading ' // TRIM(ThisFile)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Get number of volcano records 
       nVolc = 0
       DO 
          CALL GetNextLine( am_I_Root, LUN, ThisLine, EOF, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( EOF ) EXIT
          
          ! Skip any entries that contain '::'
          IF ( INDEX( TRIM(ThisLine), '::') > 0 ) CYCLE

          ! If we make it to here, this is a valid entry
          nVolc = nVolc + 1
       ENDDO 

       ! Verbose
       IF ( HCO_IsVerb(2) ) THEN
          WRITE(MSG,*) 'Number of volcanoes: ', nVolc 
          CALL HCO_MSG(MSG)
       ENDIF
 
       ! Allocate arrays
       IF ( nVolc > 0 ) THEN
          ! Eventually deallocate previously allocated data
          IF ( ALLOCATED(VolcSlf) ) DEALLOCATE(VolcSlf)
          IF ( ALLOCATED(VolcElv) ) DEALLOCATE(VolcElv)
          IF ( ALLOCATED(VolcCld) ) DEALLOCATE(VolcCld)
          IF ( ALLOCATED(VolcIdx) ) DEALLOCATE(VolcIdx)
          IF ( ALLOCATED(VolcJdx) ) DEALLOCATE(VolcJdx)

          ALLOCATE(VolcLon(nVolc), VolcLat(nVolc), VolcSlf(nVolc), &
                   VolcElv(nVolc), VolcCld(nVolc), VolcIdx(nVolc), &
                   VolcJdx(nVolc), STAT=AS )
          IF ( AS /= 0 ) THEN
             CALL HCO_ERROR ( 'Volc allocation error', RC, THISLOC=LOC )
             RETURN
          ENDIF
          VolcLon = 0.0_hp
          VolcLat = 0.0_hp
          VolcSlf = 0.0_sp
          VolcElv = 0.0_sp
          VolcCld = 0.0_sp

       ELSE
          WRITE(MSG,*) 'No volcano data found for year/mm/dd: ', YYYY, MM, DD
          CALL HCO_WARNING(MSG,RC,WARNLEV=1,THISLOC=LOC)
       ENDIF
    
       ! Now read records
       IF ( nVolc > 0 ) THEN
          REWIND( LUN )
  
          N = 0 
          DO 
             CALL GetNextLine( am_I_Root, LUN, ThisLine, EOF, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
             IF ( EOF ) EXIT
          
             ! Skip any entries that contain '::'
             IF ( INDEX( TRIM(ThisLine), '::') > 0 ) CYCLE

             ! Write this data into the following vector element
             N = N + 1
             IF ( N > nVolc ) THEN
                WRITE(MSG,*) 'N exceeds nVolc: ', N, nVolc, &
                             ' - This error occurred when reading ', &
                             TRIM(ThisFile), '. This line: ', TRIM(ThisLine)
                CALL HCO_ERROR ( MSG, RC, THISLOC = LOC )
                RETURN
             ENDIF   
     
             CALL HCO_CharSplit( TRIM(ThisLine), ' ', HCO_WCD(), Dum, nCol, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN

             ! Expect 5 values
             IF ( nCol /= 5 ) THEN
                WRITE(MSG,*) 'Cannot parse line ', TRIM(ThisLine), &
                             'Expected five entries, separated by space ', &
                             'character, instead found ', nCol
                CALL HCO_ERROR ( MSG, RC, THISLOC = LOC )
                RETURN
             ENDIF

             ! Now pass to vectors 
             VolcLat(N) = Dum(1)
             VolcLon(N) = Dum(2) 
             VolcSlf(N) = Dum(3)
             VolcElv(N) = Dum(4)
             VolcCld(N) = Dum(5)
          ENDDO

          ! At this point, we should have read exactly nVolc entries!
          IF ( N /= nVolc ) THEN
             WRITE(MSG,*) 'N /= nVolc: ', N, nVolc, &
                          ' - This error occurred when reading ', TRIM(ThisFile)
             CALL HCO_ERROR ( MSG, RC, THISLOC = LOC )
             RETURN
          ENDIF 

       ENDIF

       ! All done
       CLOSE ( LUN )

       ! Get grid box indeces for each location
       IF ( nVolc > 0 ) THEN
          CALL HCO_GetHorzIJIndex( am_I_Root, HcoState, nVolc, &
                                   VolcLon,   VolcLat,  VolcIdx, VolcJdx, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

    ENDIF ! new day

    ! Cleanup
    IF ( ALLOCATED(VolcLon) ) DEALLOCATE(VolcLon)
    IF ( ALLOCATED(VolcLat) ) DEALLOCATE(VolcLat)

    ! Return w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE ReadVolcTable 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EmitVolc 
!
! !DESCRIPTION: Subroutine EmitVolc reads the AeroCom volcano table of the
! current day. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EmitVolc ( am_I_Root, HcoState, ExtState, SO2d, SO2e, RC ) 
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options      
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState   ! Hemco state 
    INTEGER,          INTENT(INOUT) :: RC 
!
! !OUTPUT PARAMETERS:
!
    REAL(sp),         INTENT(  OUT) :: SO2e(HcoState%NX,HcoState%NY,HcoState%NZ)
    REAL(sp),         INTENT(  OUT) :: SO2d(HcoState%NX,HcoState%NY,HcoState%NZ)
!
! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N
    LOGICAL            :: Erupt
    REAL(sp)           :: nSO2, zTop, zBot, PlumeHgt
    REAL(sp)           :: z1,   z2
    REAL(sp)           :: tmp1, tmp2, Frac
    REAL(sp)           :: totE, totD, volcE, volcD
    CHARACTER(LEN=255) :: MSG 
    CHARACTER(LEN=255) :: LOC = 'EmitVolc (hcox_aerocom_mod.F90)'

    !=================================================================
    ! EmitVolc begins here!
    !=================================================================

    ! Init
    SO2e = 0.0_sp
    SO2d = 0.0_sp
    totE = 0.0_sp
    totD = 0.0_sp

    ! Make sure all required grid quantities are defined
    IF ( .NOT. ASSOCIATED(HcoState%Grid%AREA_M2%Val) ) THEN
       CALL HCO_ERROR ( 'Grid box areas not defined', RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( .NOT. ASSOCIATED(HcoState%Grid%ZSFC%Val) ) THEN
       CALL HCO_ERROR ( 'surface heights not defined', RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( .NOT. ASSOCIATED(HcoState%Grid%BXHEIGHT_M%Val) ) THEN
       CALL HCO_ERROR ( 'Grid box heights not defined', RC, THISLOC=LOC )
       RETURN
    ENDIF
 
    ! Do for every volcano
    IF ( nVolc > 0 ) THEN
       DO N = 1, nVolc

          ! Grid box index for this volcano
          I = VolcIdx(N)
          J = VolcJdx(N)

          ! Skip if outside of domain
          IF( I < 1 .OR. J < 1 ) CYCLE

          ! total emissions of this volcano
          volcE = 0.0_sp
          volcD = 0.0_sp

          z1 = HcoState%Grid%ZSFC%Val(I,J)

          ! Get total emitted kgS/m2/s. Data in table is in kgS/s.
          nSo2 = VolcSlf(N) / HcoState%Grid%AREA_M2%Val(I,J)

          ! Elevation of volcano base and volcano cloud top height [m]
          ! Make sure that the bottom / top are at least at surface level
          zBot = MAX(VolcElv(N),z1)
          zTop = MAX(VolcCld(N),z1)

          ! If volcano is eruptive, zBot /= zTop. In this case, evenly 
          ! distribute emissions in top 1/3 of the plume
          IF ( zBot /= zTop ) THEN
             zBot  = zTop - ( ( zTop - zBot ) / 3.0_sp )
             Erupt = .TRUE.
          ELSE
             Erupt = .FALSE.
          ENDIF       

          ! Volcano plume height
          PlumeHgt = zTop - zBot

          ! Distribute emissions into emission arrays. The volcano plume
          ! ranges from zBot to zTop. 
          DO L = 1, HcoState%NZ

             ! Get top height of this box
             z2 = z1 + HcoState%Grid%BXHEIGHT_M%Val(I,J,L)

             ! Skip if the plume bottom is above this grid box top
             IF ( zBot >= z2 ) THEN
                z1 = z2
                CYCLE
             ENDIF

             ! If the plume top is below this grid box bottom, we can exit
             ! since there will be no more emissions to distribute.
             IF ( zTop < z1 ) EXIT

             ! If we make it to here, the volcano plume is at least partly
             ! within this level. Determine the fraction of the plume that
             ! is within heights z1 to z2. 
            
             ! Get the bottom and top height of the plume within this layer. 
             tmp1 = MAX(z1,zBot)  ! this layer's plume bottom
             tmp2 = MIN(z2,zTop)  ! this layer's plume top

             ! Special case that zTop is heigher than the highest level: make
             ! sure that all emissions are going to be used.
             IF ( ( L == HcoState%NZ ) .AND. ( zTop > z2 ) ) THEN
                tmp2 = zTop
             ENDIF

             ! Fraction of total plume that is within this layer
             IF ( PlumeHgt == 0.0_sp ) THEN
                Frac = 1.0_sp
             ELSE
                Frac = (tmp2-tmp1) / PlumeHgt 
             ENDIF              
 
             ! Distribute emissions
             IF ( Erupt ) THEN
                SO2e(I,J,L) = SO2e(I,J,L) + ( Frac * nSo2 ) 
                volcE       = volcE &
                            + ( Frac * nSo2 * HcoState%Grid%AREA_M2%Val(I,J) )
             ELSE
                SO2d(I,J,L) = SO2d(I,J,L) + ( Frac * nSo2 ) 
                volcD       = volcD &
                            + ( Frac * nSo2 * HcoState%Grid%AREA_M2%Val(I,J) )
             ENDIF 

             ! The top height is the new bottom
             z1 = z2
          ENDDO

          ! testing 
          !IF ( HCO_IsVerb(3) ) THEN
          !   WRITE(MSG,*) 'Total eruptive  emissions of volcano ', N, ' [kgS/s]: ', volcE
          !   CALL HCO_MSG(MSG)
          !   WRITE(MSG,*) 'Total degassing emissions of volcano ', N, ' [kgS/s]: ', volcD
          !   CALL HCO_MSG(MSG)
          !ENDIF

          ! total 
          totE = totE + volcE
          totD = totD + volcD
 
       ENDDO 
    ENDIF

    ! verbose
    IF ( HCO_IsVerb(3) ) THEN
       WRITE(MSG,*) 'Total eruptive  emissions [kgS/s]: ', totE 
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) 'Total degassing emissions [kgS/s]: ', totD
       CALL HCO_MSG(MSG)
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE EmitVolc 
!EOC
END MODULE HCOX_AeroCom_Mod
