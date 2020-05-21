!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_calc_mod.F90
!
! !DESCRIPTION: Module HCO\_Calc\_Mod contains routines to calculate
! HEMCO core emissions based on the content of the HEMCO EmisList
! object. All emissions are in [kg/m2/s].
!\\
!\\
! Emissions for the current datetime are calculated by multiplying base
! emissions fields with the associated scale factors. Different
! inventories are merged/overlayed based upon the category and hierarchy
! attributes assigned to the individual base fields. Within the same
! category, fields of higher hierarchy overwrite lower-hierarchy fields.
! Emissions of different categories are always added.
!\\
!\\
! The assembled emission array is written into the corresponding emission
! rates array of the HEMCO state object: HcoState%Spc(HcoID)%Emis, where
! HcoID denotes the corresponding species ID. Emis covers dimension lon,
! lat, lev on the HEMCO grid, i.e. unlike the emission arrays in EmisList
! that only cover the levels defined in the source files, Emis extends
! over all vertical model levels.
!\\
!\\
! Negative emissions are not supported and are ignored. Surface
! deposition velocities are stored in HcoState%Spc(HcoID)%Depv and can
! be added therein.
!\\
!\\
! In addition to emissions and surface deposition rates, HEMCO also
! supports concentrations (kg/m3). Data is automatically written into
! the concentration array HcoState%Spc(HcoID)%Conc if the source data
! is marked as being concentration data (i.e. if Dta%IsConc is .TRUE.).
! The latter is automatically determined by HEMCO based upon the data
! units.
!\\
!\\
! All emission calculation settings are passed through the HcoState
! options object (HcoState%Options). These include:
!
! \begin{itemize}
!  \item ExtNr: extension number to be considered.
!  \item SpcMin: lower species ID (HEMCO ID) to be considered.
!  \item SpcMax: upper species ID (HEMCO ID) to be considered. If set
!        to -1, all species above or equal to SpcMin are considered.
!  \item CatMin: lower emission category to be considered.
!  \item CatMax: upper emission category to be considered. If set to
!        -1, all categories above or equal to CatMin are considered.
!  \item FillBuffer: if set to TRUE, the emissions will be written into
!        buffer array HcoState%Buffer3D instead of HcoState%Spc(ID)%Emis.
!        If this option is enabled, only one species can be calculated at
!        a time (by setting SpcMin/SpcMax, CatMin/CatMax and/or ExtNr
!        accordingly). This option is useful for extensions, e.g. if
!        additional scalings need to be done on some emission fields
!        assembled by HEMCO (e.g. PARANOX extension).
! \end{itemize}
!
! !INTERFACE:
!
MODULE HCO_Calc_Mod
!
! !USES:
!
  USE HCO_Diagn_Mod
  USE HCO_Error_Mod
  USE HCO_Types_Mod
  USE HCO_DataCont_Mod, ONLY : Pnt2DataCont

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCO_CalcEmis
  PUBLIC  :: HCO_CheckDepv
  PUBLIC  :: HCO_EvalFld
  PUBLIC  :: HCO_MaskFld
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: GET_CURRENT_EMISSIONS
  PRIVATE :: GetMaskVal
  PRIVATE :: GetDilFact
  PRIVATE :: GetVertIndx
  PRIVATE :: GetIdx
!
! !PARAMETER
!
  ! Mask threshold. All mask values below this value will be evaluated
  ! as zero (= outside of mask), and all values including and above this
  ! value as inside the mask. This is only of relevance if the MaskFractions
  ! option is false. If MaskFractions is true, the fractional mask values are
  ! considered, e.g. a grid box can contribute 40% to a mask region, etc.
  ! The MaskFractions toggle can be set in the settings section of the HEMCO
  ! configuration file (Use mask fractions: true/false). It defaults to false.
  REAL(sp), PARAMETER  :: MASK_THRESHOLD = 0.5_sp
!
! ============================================================================
!
! !REVISION HISTORY:
!  25 Aug 2012 - C. Keller   - Initial version.
!  06 Jun 2014 - R. Yantosca - Add cosmetic changes in ProTeX headers
!  08 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  29 Dec 2014 - C. Keller   - Added MASK_THRESHOLD parameter. Added option to
!                              apply scale factors only over masked area.
!  08 Apr 2015 - C. Keller   - Added option for fractional masks.
!  11 May 2015 - C. Keller   - Added HCO_EvalFld interface.
!EOP
!------------------------------------------------------------------------------
!BOC
  INTERFACE HCO_EvalFld
     MODULE PROCEDURE HCO_EvalFld_2D
     MODULE PROCEDURE HCO_EvalFld_3D
  END INTERFACE HCO_EvalFld

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_CalcEmis
!
! !DESCRIPTION: Subroutine HCO\_CalcEmis calculates the 3D emission
! fields at current datetime for the specified species, categories, and
! extension number.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_CalcEmis( HcoState, UseConc, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
    USE HCO_ARR_MOD,      ONLY : HCO_ArrAssert
    USE HCO_DATACONT_MOD, ONLY : ListCont_NextCont
    USE HCO_FILEDATA_MOD, ONLY : FileData_ArrIsDefined
    USE HCO_Scale_Mod,    ONLY : HCO_ScaleArr
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: UseConc    ! Use concentration fields?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState   ! HEMCO state object
    INTEGER,         INTENT(INOUT)  :: RC         ! Return code
!
! !REVISION HISTORY:
!  25 Aug 2012 - C. Keller   - Initial Version
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX header
!  03 Aug 2014 - C. Keller   - Bug fix for adding data to diagnostics. Now
!                              explicitly check for new species OR category.
!  21 Aug 2014 - C. Keller   - Added concentration.
!  14 Apr 2016 - C. Keller   - Bug fix: avoid double-counting if multiple
!                              regional inventories have the same hierarchy.
!  19 Sep 2016 - R. Yantosca - Use .neqv. for LOGICAL comparisons
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Working pointers: list and data container
    TYPE(ListCont), POINTER :: Lct
    TYPE(DataCont), POINTER :: Dct

    ! Temporary emission arrays
    REAL(hp), POINTER       :: OutArr(:,:,:) => NULL()
    REAL(hp), TARGET        :: SpcFlx( HcoState%NX, &
                                       HcoState%NY, &
                                       HcoState%NZ   )
    REAL(hp), TARGET        :: CatFlx( HcoState%NX, &
                                       HcoState%NY, &
                                       HcoState%NZ   )
    REAL(hp), TARGET        :: TmpFlx( HcoState%NX, &
                                       HcoState%NY, &
                                       HcoState%NZ   )
    REAL(hp)                :: Mask  ( HcoState%NX, &
                                       HcoState%NY, &
                                       HcoState%NZ   )
    REAL(hp)                :: HirFlx( HcoState%NX, &
                                       HcoState%NY, &
                                       HcoState%NZ   )
    REAL(hp)                :: HirMsk( HcoState%NX, &
                                       HcoState%NY, &
                                       HcoState%NZ   )

    ! Integers
    INTEGER             :: ThisSpc, PrevSpc ! current and previous species ID
    INTEGER             :: ThisCat, PrevCat ! current and previous category
    INTEGER             :: ThisHir, PrevHir ! current and previous hierarchy
    INTEGER             :: SpcMin,  SpcMax  ! range of species to be considered
    INTEGER             :: CatMin,  CatMax  ! range of categories to be considered
    INTEGER             :: ExtNr            ! Extension Nr to be used
    INTEGER             :: nI, nJ, nL
    INTEGER             :: nnSpec, FLAG

    LOGICAL             :: Found, DoDiagn, EOL, UpdateCat

    ! For error handling & verbose mode
    CHARACTER(LEN=255)  :: MSG

    ! testing / debugging
    integer :: ix,iy

    !=================================================================
    ! HCO_CalcEmis begins here!
    !=================================================================

    ! testing only
    ix = 30
    iy = 34

    ! Initialize
    Lct => NULL()
    Dct => NULL()

    ! Enter routine
    CALL HCO_ENTER (HcoState%Config%Err,'HCO_CalcEmis (HCO_CALC_MOD.F90)', RC )
    IF(RC /= HCO_SUCCESS) RETURN

    !-----------------------------------------------------------------
    ! Initialize variables
    !-----------------------------------------------------------------

    ! Initialize
    SpcFlx(:,:,:)    = 0.0_hp
    CatFlx(:,:,:)    = 0.0_hp
    HirFlx(:,:,:)    = 0.0_hp
    HirMsk(:,:,:)    = 0.0_hp
    PrevSpc          = -1
    PrevHir          = -1
    PrevCat          = -1
    nnSpec           = 0

    ! Pass emission grid dimensions
    nI = HcoState%NX
    nJ = HcoState%NY
    nL = HcoState%NZ

    ! Pass calculation options
    SpcMin  = HcoState%Options%SpcMin        !Lower species ID
    SpcMax  = HcoState%Options%SpcMax        !Upper species ID
    CatMin  = HcoState%Options%CatMin        !Lower emission category
    CatMax  = HcoState%Options%CatMax        !Upper emission category
    ExtNr   = HcoState%Options%ExtNr         !Extension number
    DoDiagn = HcoState%Options%AutoFillDiagn !Write AutoFill diagnostics?

    ! Verbose mode
    IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
       WRITE (MSG, *) 'Run HEMCO calculation w/ following options:'
       CALL HCO_MSG ( HcoState%Config%Err, MSG )
       WRITE (MSG, "(A20,I5)")    'Extension number:', ExtNr
       CALL HCO_MSG ( HcoState%Config%Err, MSG )
       WRITE (MSG, "(A20,I5,I5)") 'Tracer range:', SpcMin, SpcMax
       CALL HCO_MSG ( HcoState%Config%Err, MSG )
       WRITE (MSG, "(A20,I5,I5)") 'Category range:', CatMin, CatMax
       CALL HCO_MSG ( HcoState%Config%Err, MSG )
       WRITE (MSG, *) 'Auto diagnostics: ', DoDiagn
       CALL HCO_MSG ( HcoState%Config%Err, MSG )
    ENDIF

    !=================================================================
    ! Walk through all containers of EmisList and determine the
    ! emissions for all containers that qualify for calculation.
    ! The containers in EmisList are sorted by species, category and
    ! hierarchy. This enables a straightforward, piece-by-piece
    ! assembly of the final emission array (start with lowest
    ! hierarchy emissions, then overwrite piece-by-piece with higher
    ! hierarchy values).
    !=================================================================

    ! Point to the head of the emissions linked list
    EOL = .FALSE. ! End of list
    Lct => NULL()
    CALL ListCont_NextCont ( HcoState%EmisList, Lct, FLAG )

    ! Do until end of EmisList (==> loop over all emission containers)
    DO
       ! Have we reached the end of the list?
       IF ( FLAG /= HCO_SUCCESS ) THEN
          EOL = .TRUE.
       ELSE
          EOL = .FALSE.
       ENDIF

       ! ------------------------------------------------------------
       ! Select container and update all working variables & arrays.
       ! ------------------------------------------------------------
       IF ( .NOT. EOL ) THEN

          ! Dct is the current data container
          Dct => Lct%Dct

          ! Check if this is a base field
          IF ( Dct%DctType /= HCO_DCTTYPE_BASE ) THEN
             CALL ListCont_NextCont ( HcoState%EmisList, Lct, FLAG )
             CYCLE
          ENDIF

          ! Sanity check: Make sure this container holds data.
          ! 'Empty' containers are possible if the simulation time
          ! is outside of the specified data time range and time
          ! slice cycling is deactivated (CycleFlag > 1).
          IF( .NOT. FileData_ArrIsDefined(Lct%Dct%Dta) ) THEN
             CALL ListCont_NextCont ( HcoState%EmisList, Lct, FLAG )
             CYCLE
          ENDIF

          ! Check if this is the specified extension number
          IF ( Dct%ExtNr /= ExtNr ) THEN
             CALL ListCont_NextCont ( HcoState%EmisList, Lct, FLAG )
             CYCLE
          ENDIF

          ! Advance to next container if the species ID is outside
          ! the specified species range (SpcMin - SpcMax). Consider
          ! all species above SpcMin if SpcMax is negative!
          IF( (  Dct%HcoID < SpcMin                     ) .OR. &
              ( (Dct%HcoID > SpcMax) .AND. (SpcMax > 0) ) ) THEN
             CALL ListCont_NextCont ( HcoState%EmisList, Lct, FLAG )
             CYCLE
          ENDIF

          ! Advance to next emission field if the emission category of
          ! the current container is outside of the specified species
          ! range (CatMin - CatMax). Consider all categories above CatMin
          ! if CatMax is negative!
          IF( (  Dct%Cat < CatMin                     ) .OR. &
              ( (Dct%Cat > CatMax) .AND. (CatMax > 0) ) ) THEN
             CALL ListCont_NextCont ( HcoState%EmisList, Lct, FLAG )
             CYCLE
          ENDIF

          ! Check if this container holds data in the desired unit format,
          ! i.e. concentration data if UseConc is enabled, emission data
          ! otherwise.
          IF ( UseConc .NEQV. Dct%Dta%IsConc ) THEN
             CALL ListCont_NextCont ( HcoState%EmisList, Lct, FLAG )
             CYCLE
          ENDIF

          ! Update working variables
          ThisSpc = Dct%HcoID
          ThisCat = Dct%Cat
          ThisHir = Dct%Hier

       ! If end of list, use dummy values for ThisSpc, ThisCat and ThisHir
       ! to make sure that emissions are added to HEMCO in the section
       ! below!
       ELSE
          ThisSpc = -1
          ThisCat = -1
          ThisHir = -1
       ENDIF

       !--------------------------------------------------------------------
       ! Before computing emissions of current data container make sure that
       ! emissions of previous container are properly archived.
       !--------------------------------------------------------------------

       ! Add emissions on hierarchy level to the category flux array. Do
       ! this only if this is a new species, a new category or a new
       ! hierarchy level.
       ! Note: no need to add to diagnostics because hierarchy level
       ! diagnostics are filled right after computing the emissions of
       ! a given data container (towards the end of the DO loop).
       IF ( (ThisHir /= PrevHir) .OR. &
            (ThisSpc /= PrevSpc) .OR. &
            (ThisCat /= PrevCat)        ) THEN

          ! Add hierarchy level emissions to category array over the
          ! covered regions.
          CatFlx = ( (1.0_hp - HirMsk) * CatFlx ) + HirFlx

          ! Reset
          HirFlx = 0.0_hp
          HirMsk = 0.0_hp
       ENDIF

       !--------------------------------------------------------------------
       ! If this is a new species or category, pass the previously collected
       ! emissions to the species array. Update diagnostics at category level.
       ! Skip this step for first species, i.e. if PrevSpc is still -1.
       !--------------------------------------------------------------------
       UpdateCat = .FALSE.
       IF ( ThisCat /= PrevCat ) UpdateCat = .TRUE.
       IF ( ThisSpc /= PrevSpc ) UpdateCat = .TRUE.
       IF ( PrevCat <= 0 .OR. PrevSpc <= 0 ) UpdateCat = .FALSE.
       IF ( UpdateCat ) THEN

          ! CatFlx holds the emissions for this category. Pass this to
          ! the species array SpcFlx.
          SpcFlx(:,:,:) = SpcFlx(:,:,:) + CatFlx(:,:,:)

          ! verbose
          IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
             WRITE(MSG,*) 'Added category emissions to species array: '
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             WRITE(MSG,*) 'Species       : ', PrevSpc
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             WRITE(MSG,*) 'Category      : ', PrevCat
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             WRITE(MSG,*) 'Cat. emissions: ', SUM(CatFlx)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             WRITE(MSG,*) 'Spc. emissions: ', SUM(SpcFlx)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF

          ! Add category emissions to diagnostics at category level
          ! (only if defined in the diagnostics list).
          IF ( Diagn_AutoFillLevelDefined(HcoState%Diagn,3) .AND. DoDiagn ) THEN
             CALL Diagn_Update( HcoState,    ExtNr=ExtNr,   &
                                Cat=PrevCat, Hier=-1,  HcoID=PrevSpc, &
                                AutoFill=1,  Array3D=CatFlx, COL=-1, RC=RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          ! Reset CatFlx array and the previously used hierarchy
          ! ==> Emission hierarchies are only important within the
          ! same category, hence always start over at lowest hierarchy
          ! when entering a new category.
          CatFlx(:,:,:)  = 0.0_hp
          PrevHir        = -1
       ENDIF

       !--------------------------------------------------------------------
       ! If this is a new species, pass previously calculated emissions
       ! to the final emissions array in HcoState.
       ! Update diagnostics at extension number level.
       ! Don't do before first emission calculation, i.e. if PrevSpc
       ! is still the initialized value of -1!
       !--------------------------------------------------------------------
       IF ( ThisSpc /= PrevSpc .AND. PrevSpc > 0 ) THEN

          ! Add to OutArr
          OutArr(:,:,:) = OutArr(:,:,:) + SpcFlx(:,:,:)

          ! testing only
          IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
             WRITE(MSG,*) 'Added total emissions to output array: '
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             WRITE(MSG,*) 'Species: ', PrevSpc
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             WRITE(MSG,*) 'SpcFlx : ', SUM(SpcFlx)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             WRITE(MSG,*) 'OutArr : ', SUM(OutArr)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF

          ! Add to diagnostics at extension number level.
          ! The same diagnostics may be updated multiple times during
          ! the same time step, continuously adding emissions to it.
          IF ( Diagn_AutoFillLevelDefined(HcoState%Diagn,2) .AND. DoDiagn ) THEN
             CALL Diagn_Update( HcoState,  ExtNr=ExtNr,  &
                                Cat=-1,    Hier=-1,  HcoID=PrevSpc, &
                                AutoFill=1,Array3D=SpcFlx, COL=-1, RC=RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          ! Reset arrays and previous hierarchy.
          SpcFlx(:,:,:)  =  0.0_hp
          PrevCat        =  -1
          PrevHir        =  -1
          OutArr         => NULL()
       ENDIF

       !--------------------------------------------------------------------
       ! Exit DO loop here if end of list
       !--------------------------------------------------------------------
       IF ( EOL ) EXIT

       !--------------------------------------------------------------------
       ! Update/archive information on species level if needed
       !--------------------------------------------------------------------
       IF ( ThisSpc /= PrevSpc .AND. ThisSpc > 0 ) THEN

          ! Update number of species for which emissions have been
          ! calculated.
          nnSpec = nnSpec + 1

          ! To write emissions into temporary array, make OutArr point
          ! to the buffer array HcoState%Buffer3D.
          IF ( HcoState%Options%FillBuffer ) THEN

             ! Cannot use temporary array for more than one species!
             IF ( nnSpec > 1 ) THEN
                MSG = 'Cannot fill buffer for more than one species!'
                CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
                RETURN
             ENDIF

             ! Point to array and check allocation status as well as
             ! array size.
             OutArr => HcoState%Buffer3D%Val
             IF ( .NOT. ASSOCIATED( OutArr ) ) THEN
                MSG = 'Buffer array is not associated'
                CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
                RETURN
             ENDIF
             IF ( (SIZE(OutArr,1) /= nI) .OR. &
                  (SIZE(OutArr,2) /= nJ) .OR. &
                  (SIZE(OutArr,3) /= nL)       ) THEN
                MSG = 'Buffer array has wrong dimension!'
                CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
                RETURN
             ENDIF

          ! To write emissions directly into HcoState, make OutArr
          ! point to current species' array in HcoState. Use emission
          ! array for emissions, and concentration array for concentrations.
          ELSE

             ! For concentrations:
             IF ( UseConc ) THEN
                CALL HCO_ArrAssert( HcoState%Spc(ThisSpc)%Conc, &
                                    nI, nJ, nL, RC             )
                IF ( RC /= HCO_SUCCESS ) RETURN
                OutArr => HcoState%Spc(ThisSpc)%Conc%Val

             ! For emissions:
             ELSE
                CALL HCO_ArrAssert( HcoState%Spc(ThisSpc)%Emis, &
                                    nI, nJ, nL, RC             )
                IF ( RC /= HCO_SUCCESS ) RETURN
                OutArr => HcoState%Spc(ThisSpc)%Emis%Val
             ENDIF

          ENDIF

          ! verbose mode
          IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
             WRITE(MSG,*) 'Calculating emissions for species ', &
                           TRIM(HcoState%Spc(ThisSpc)%SpcName)
             CALL HCO_MSG( HcoState%Config%Err, MSG, SEP1='-', SEP2='-' )
          ENDIF
       ENDIF

       !--------------------------------------------------------------------
       ! Get current emissions and write into TmpFlx array. The array Mask
       ! denotes all valid grid boxes for this inventory.
       !--------------------------------------------------------------------
       TmpFlx(:,:,:) = 0.0_hp
       CALL GET_CURRENT_EMISSIONS( HcoState, Dct, nI, nJ, nL, TmpFlx, Mask, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Eventually add universal scale factor
       CALL HCO_ScaleArr( HcoState, ThisSpc, TmpFlx, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Check for negative values according to the corresponding setting
       ! in the configuration file: 2 means allow negative values, 1 means
       ! set to zero and prompt a warning, else return with error.
       IF ( HcoState%Options%NegFlag /= 2 ) THEN

          IF ( ANY(TmpFlx < 0.0_hp) ) THEN

             ! Set to zero and prompt warning
             IF ( HcoState%Options%NegFlag == 1 ) THEN
                WHERE ( TmpFlx < 0.0_hp ) TmpFlx = 0.0_hp
                MSG = 'Negative emissions set to zero: '// TRIM(Dct%cName)
                CALL HCO_WARNING( HcoState%Config%Err, MSG, RC )

             ! Return with error
             ELSE
                MSG = 'Negative emissions in: '// TRIM(Dct%cName) // '. ' // &
                'To allow negatives, edit settings in the configuration file.'
                CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
                RETURN
             ENDIF
          ENDIF
       ENDIF

       ! ------------------------------------------------------------
       ! Collect all emissions of the same category (and species) on
       ! the hierarchy level into array HirFlx. HirMsk contains the
       ! combined covered region. That is, if there are two regional
       ! inventories with the same hierarchy HirMsk will cover both
       ! of these regions.
       ! The specified field hierarchies determine whether the
       ! temporary emissions are added (if hierarchy is the same
       ! as the previously used hierarchy), or if they overwrite the
       ! previous values in HirFlx (if hierarchy is higher than the
       ! previous hierarchy).
       ! ------------------------------------------------------------

       ! Add emissions to the hierarchy array HirFlx if this hierarchy
       ! is the same as previous hierarchy
       IF ( ThisHir == PrevHir ) THEN
          HirFlx = HirFlx + TmpFlx
          HirMsk = HirMsk + Mask

          ! Make sure mask values do not exceed 1.0
          WHERE(HirMsk > 1.0 ) HirMsk = 1.0

       ! If hierarchy is larger than those of the previously used
       ! fields, overwrite HirFlx with new values.
       ELSE

          HirFlx = TmpFlx
          HirMsk = Mask

       ENDIF

       ! Update diagnostics at hierarchy level. Make sure that only
       ! positive values are used.
       ! The same diagnostics may be updated multiple times
       ! during the same time step, continuously adding
       ! emissions to it.
       ! Now remove PosOnly flag. TmpFlx is initialized to zero, so it's
       ! ok to keep negative values (ckeller, 7/12/15).
       IF ( Diagn_AutoFillLevelDefined(HcoState%Diagn,4) .AND. DoDiagn ) THEN
          CALL Diagn_Update( HcoState,       ExtNr=ExtNr,   &
                             Cat=ThisCat,Hier=ThisHir,   HcoID=ThisSpc, &
                             !AutoFill=1, Array3D=TmpFlx, PosOnly=.TRUE.,&
                             AutoFill=1, Array3D=TmpFlx, &
                             COL=-1, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! Update previously used species, category and hierarchy
       PrevSpc = ThisSpc
       PrevCat = ThisCat
       PrevHir = ThisHir

       ! Advance to next emission container
       CALL ListCont_NextCont( HcoState%EmisList, Lct, FLAG )

    ENDDO ! Loop over EmisList

    ! Make sure internal pointers are nullified
    Lct    => NULL()
    Dct    => NULL()
    OutArr => NULL()

    ! verbose
    IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
       WRITE (MSG, *) 'HEMCO emissions successfully calculated!'
       CALL HCO_MSG ( HcoState%Config%Err, MSG )
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE HCO_CalcEmis
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_CheckDepv
!
! !DESCRIPTION: Subroutine HCO\_CheckDepv is a simple routine to check the
! dry deposition frequency value. This is to avoid unrealistically high
! deposition frequencies that may occur if grid box concentrations are very
! low. The deposition frequency is limited to a value that will make sure
! that the drydep exponent ( exp( -depfreq * dt ) ) is still small enough to
! remove all species mass. The maximum limit of depfreq * dt can be defined
! as a HEMCO option (MaxDepExp). Its default value is 20.0.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_CheckDepv( HcoState, Depv, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState   ! HEMCO state object
    REAL(hp),        INTENT(INOUT)  :: Depv       ! Deposition velocity
    INTEGER,         INTENT(INOUT)  :: RC         ! Return code
!
! !REVISION HISTORY:
!  11 Mar 2015 - C. Keller   - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp)  :: ExpVal

    !=================================================================
    ! HCO_CheckDepv begins here!
    !=================================================================

    ExpVal = Depv * HcoState%TS_EMIS
    IF ( ExpVal > HcoState%Options%MaxDepExp ) THEN
       Depv = HcoState%Options%MaxDepExp / HcoState%TS_EMIS
    ENDIF

  END SUBROUTINE HCO_CheckDepv
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Current_Emissions
!
! !DESCRIPTION: Subroutine Get\_Current\_Emissions calculates the current
!  emissions for the specified emission container.
!  This subroutine is only called by HCO\_CalcEmis and for base emission
!  containers, i.e. containers of type 1.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Current_Emissions( HcoState,   BaseDct,   &
                                    nI, nJ, nL, OUTARR_3D, MASK, RC, UseLL )
!
! !USES:
!
    USE HCO_State_Mod,    ONLY : HCO_State
    USE HCO_tIdx_MOD,     ONLY : tIDx_GetIndx
    USE HCO_FileData_Mod, ONLY : FileData_ArrIsDefined
!
! !INPUT PARAMETERS:
!
    INTEGER,           INTENT(IN)  :: nI                  ! # of lons
    INTEGER,           INTENT(IN)  :: nJ                  ! # of lats
    INTEGER,           INTENT(IN)  :: nL                  ! # of levs
!
! !INPUT/OUTPUT PARAMETERS:
!

    TYPE(HCO_State), POINTER       :: HcoState            ! HEMCO state object
    TYPE(DataCont),  POINTER       :: BaseDct             ! base emission
                                                          !  container
    REAL(hp),        INTENT(INOUT) :: OUTARR_3D(nI,nJ,nL) ! output array
    REAL(hp),        INTENT(INOUT) :: MASK     (nI,nJ,nL) ! mask array
    INTEGER,         INTENT(INOUT) :: RC
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(  OUT), OPTIONAL :: UseLL
!
! !REMARKS:
!  This routine uses multiple loops over all grid boxes (base emissions
!  and scale factors use separate loops). In an OMP environment, this approach
!  seems to be faster than using only one single loop (but repeated calls to
!  point to containers, etc.). The alternative approach is used in routine
!  Get\_Current\_Emissions\_B at the end of this module and may be employed
!  on request.
!
! !REVISION HISTORY:
!  25 Aug 2012 - C. Keller   - Initial Version
!  09 Nov 2012 - C. Keller   - MASK update. Masks are now treated
!                              separately so that multiple masks can be
!                              added.
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  07 Sep 2014 - C. Keller   - Mask update. Now set mask to zero as soon as
!                              on of the applied masks is zero.
!  03 Dec 2014 - C. Keller   - Now calculate time slice index on-the-fly.
!  29 Dec 2014 - C. Keller   - Added scale factor masks.
!  02 Mar 2015 - C. Keller   - Now check for missing values. Missing values are
!                              excluded from emission calculation.
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  11 May 2017 - C. Keller   - Added universal scaling
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DataCont), POINTER :: ScalDct
    TYPE(DataCont), POINTER :: MaskDct
    TYPE(DataCont), POINTER :: LevDct1
    TYPE(DataCont), POINTER :: LevDct2

    ! Scalars
    REAL(sp)                :: TMPVAL, MaskScale
    REAL(hp)                :: DilFact
    REAL(hp)                :: ScalFact
    INTEGER                 :: tIDx, IDX
    INTEGER                 :: I, J, L, N
    INTEGER                 :: LowLL, UppLL, ScalLL, TmpLL
    INTEGER                 :: ERROR
    INTEGER                 :: TotLL, nnLL
    CHARACTER(LEN=255)      :: MSG, LOC
    LOGICAL                 :: NegScalExist
    LOGICAL                 :: MaskFractions

    ! testing only
    INTEGER                 :: IX, IY

    !=================================================================
    ! GET_CURRENT_EMISSIONS begins here
    !=================================================================

    ! Initialize
    ScalDct => NULL()
    MaskDct => NULL()

    ! Enter
    CALL HCO_ENTER(HcoState%Config%Err,'GET_CURRENT_EMISSIONS', RC )
    IF(RC /= HCO_SUCCESS) RETURN

    ! testing only:
    IX = 25 !-1
    IY = 25 !-1

    ! Check if container contains data
    IF ( .NOT. FileData_ArrIsDefined(BaseDct%Dta) ) THEN
       MSG = 'Array not defined: ' // TRIM(BaseDct%cName)
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF

    ! Initialize mask. By default, assume that we use all grid boxes.
    MASK(:,:,:)  = 1.0_hp
    MaskFractions = HcoState%Options%MaskFractions

    ! Verbose
    IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
       WRITE(MSG,*) 'Evaluate field ', TRIM(BaseDct%cName)
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1=' ')
    ENDIF

    ! ----------------------------------------------------------------
    ! Set base emissions
    ! ----------------------------------------------------------------

    ! Initialize ERROR. Will be set to 1 if error occurs below
    ERROR = 0

    ! Initialize variables to compute average vertical level index
    totLL = 0
    nnLL  = 0

    ! Check for level index containers
    IF ( BaseDct%levScalID1 > 0 ) THEN
       CALL Pnt2DataCont( HcoState, BaseDct%levScalID1, LevDct1, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ELSE
       LevDct1 => NULL()
    ENDIF
    IF ( BaseDct%levScalID2 > 0 ) THEN
       CALL Pnt2DataCont( HcoState, BaseDct%levScalID2, LevDct2, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ELSE
       LevDct2 => NULL()
    ENDIF

    ! Loop over all latitudes and longitudes
!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT( SHARED )                                                &
!$OMP PRIVATE( I, J, L, tIdx, TMPVAL, DilFact, LowLL, UppLL          ) &
!$OMP SCHEDULE( DYNAMIC )
    DO J = 1, nJ
    DO I = 1, nI

       ! Get current time index for this container and at this location
       tIDx = tIDx_GetIndx( HcoState, BaseDct%Dta, I, J )
       IF ( tIDx < 1 ) THEN
          WRITE(MSG,*) 'Cannot get time slice index at location ',I,J,&
                       ': ', TRIM(BaseDct%cName), tIDx
          ERROR = 1
          EXIT
       ENDIF

       ! Get lower and upper vertical index
       CALL GetVertIndx ( HcoState, BaseDct, LevDct1, LevDct2, &
                          I, J, LowLL, UppLL, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          WRITE(MSG,*) 'Error getting vertical index at location ',I,J,&
                       ': ', TRIM(BaseDct%cName)
          ERROR = 1 ! Will cause error
          EXIT
       ENDIF

       ! average upper level
       totLL = totLL + UppLL
       nnLL  = nnLL + 1

       ! Loop over all levels
       DO L = LowLL, UppLL

          ! Get base value. Use uniform value if scalar field.
          IF ( BaseDct%Dta%SpaceDim == 1 ) THEN
             TMPVAL = BaseDct%Dta%V2(tIDx)%Val(1,1)
          ELSEIF ( BaseDct%Dta%SpaceDim == 2 ) THEN
             TMPVAL = BaseDct%Dta%V2(tIDx)%Val(I,J)
          ELSE
             TMPVAL = BaseDct%Dta%V3(tIDx)%Val(I,J,L)
          ENDIF

          ! If it's a missing value, mask box as unused and set value to
          ! zero
#if defined( ESMF_ )
          ! SDE 2017-01-07: Temporary kludge. MAPL ExtData sets missing
          ! data to 1e15, but HEMCO uses a different value!
          IF ( ( TMPVAL == HCO_MISSVAL ) .or. ( TMPVAL > 1.0e+14 ) ) THEN
#else
          IF ( TMPVAL == HCO_MISSVAL ) THEN
#endif
             MASK(I,J,:)      = 0.0_hp
             OUTARR_3D(I,J,L) = 0.0_hp

          ! Pass base value to output array
          ELSE

             ! Get dilution factor. Never dilute 3D emissions.
             IF ( BaseDct%Dta%SpaceDim == 3 ) THEN
                DilFact = 1.0

             ! 2D dilution factor
             ELSE
                CALL GetDilFact ( HcoState,    BaseDct%Dta%EmisL1, &
                                  BaseDct%Dta%EmisL1Unit, BaseDct%Dta%EmisL2,  &
                                  BaseDct%Dta%EmisL2Unit, I, J, L, LowLL,  &
                                  UppLL, DilFact, RC )
                IF ( RC /= HCO_SUCCESS ) THEN
                   WRITE(MSG,*) 'Error getting dilution factor at ',I,J,&
                                ': ', TRIM(BaseDct%cName)
                   ERROR = 1
                   EXIT
                ENDIF
             ENDIF

             ! Scale base emission by dilution factor
             OUTARR_3D(I,J,L) = DilFact * TMPVAL
          ENDIF
       ENDDO !L

    ENDDO !I
    ENDDO !J
!$OMP END PARALLEL DO

    ! Check for error
    IF ( ERROR == 1 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! Apply scale factors
    ! The container IDs of all scale factors associated with this base
    ! container are stored in vector Scal_cID.
    ! ----------------------------------------------------------------

    ! Loop over scale factors
    IF ( BaseDct%nScalID > 0 ) THEN

    DO N = 1, BaseDct%nScalID

       ! Get the scale factor container ID for the current slot
       IDX = BaseDct%Scal_cID(N)

       ! Point to data container with the given container ID
       CALL Pnt2DataCont( HcoState, IDX, ScalDct, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Sanity check: scale field cannot be a base field
       IF ( (ScalDct%DctType == HCO_DCTTYPE_BASE) ) THEN
          MSG = 'Wrong scale field type: ' // TRIM(ScalDct%cName)
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
          RETURN
       ENDIF

       ! Skip this scale factor if no data defined. This is possible
       ! if scale factors are only defined for a given time range and
       ! the simulation datetime is outside of this range.
       IF ( .NOT. FileData_ArrIsDefined(ScalDct%Dta) ) THEN
          IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
             MSG = 'Skip scale factor '//TRIM(ScalDct%cName)// &
                   ' because it is not defined for this datetime.'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
          CYCLE
       ENDIF

       ! Verbose mode
       IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
          MSG = 'Applying scale factor ' // TRIM(ScalDct%cName)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       ! Get vertical extension of this scale factor array.
       IF( (ScalDct%Dta%SpaceDim<=2) ) THEN
          ScalLL = 1
       ELSE
          ScalLL = SIZE(ScalDct%Dta%V3(1)%Val,3)
       ENDIF

       ! Check if there is a mask field associated with this scale
       ! factor. In this case, get a pointer to the corresponding
       ! mask field and evaluate scale factors only inside the mask
       ! region.
       IF ( ASSOCIATED(ScalDct%Scal_cID) ) THEN
          CALL Pnt2DataCont( HcoState, ScalDct%Scal_cID(1), MaskDct, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Must be mask field
          IF ( MaskDct%DctType /= HCO_DCTTYPE_MASK ) THEN
             MSG = 'Invalid mask for scale factor: '//TRIM(ScalDct%cName)
             MSG = TRIM(MSG) // '; mask: '//TRIM(MaskDct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF
       ENDIF

       ! Reinitialize error flag. Will be set to 1 or 2 if error occurs,
       ! and to -1 if negative scale factor is ignored.
       ERROR = 0

       ! Loop over all latitudes and longitudes
!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT( SHARED )                                                &
!$OMP PRIVATE( I, J, tIdx, TMPVAL, L, LowLL, UppLL, tmpLL, MaskScale ) &
!$OMP SCHEDULE( DYNAMIC )
       DO J = 1, nJ
       DO I = 1, nI

          ! ------------------------------------------------------------
          ! If there is a mask associated with this scale factors, check
          ! if this grid box is within or outside of the mask region.
          ! Values that partially fall into the mask region are either
          ! treated as binary (100% inside or outside), or partially
          ! (using the real grid area fractions), depending on the
          ! HEMCO options.
          ! ------------------------------------------------------------

          ! Default mask scaling is 1.0 (no mask applied)
          MaskScale = 1.0_sp

          ! If there is a mask applied to this scale factor ...
          IF ( ASSOCIATED(MaskDct) ) THEN
             CALL GetMaskVal ( MaskDct, I, J, MaskScale, MaskFractions, RC )
             IF ( RC /= HCO_SUCCESS ) THEN
                ERROR = 4
                EXIT
             ENDIF
          ENDIF

          ! We can skip this grid box if mask is completely zero
          IF ( MaskScale <= 0.0_sp ) CYCLE

          ! Get current time index for this container and at this location
          tIDx = tIDx_GetIndx( HcoState, ScalDct%Dta, I, J )
          IF ( tIDx < 1 ) THEN
             WRITE(*,*) 'Cannot get time slice index at location ',I,J,&
                          ': ', TRIM(ScalDct%cName), tIDx
             ERROR = 3
             EXIT
          ENDIF

          ! Check if this is a mask. If so, add mask values to the MASK
          ! array. For now, we assume masks to be binary, i.e. 0 or 1.
          ! We may want to change that in future to also support values
          ! in between. This is especially important when regridding
          ! high resolution masks onto coarser grids!
          ! ------------------------------------------------------------
          IF ( ScalDct%DctType == HCO_DCTTYPE_MASK ) THEN

             ! Get mask value
             CALL GetMaskVal ( ScalDct, I, J, TMPVAL, MaskFractions, RC )
             IF ( RC /= HCO_SUCCESS ) THEN
                ERROR = 4
                EXIT
             ENDIF

             ! Pass to output mask
             MASK(I,J,:) = MASK(I,J,:) * TMPVAL

             ! testing only
             IF ( HCO_IsVerb(HcoState%Config%Err,2) .AND. I==1 .AND. J==1 ) THEN
                write(MSG,*) 'Mask field ', TRIM(ScalDct%cName),   &
                     ' found and added to temporary mask.'
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF

             ! Advance to next grid box
             CYCLE
          ENDIF! DctType=MASK

          ! ------------------------------------------------------------
          ! For non-mask fields, apply scale factors to all levels
          ! of the base field individually. If the scale factor
          ! field has more than one vertical level, use the
          ! vertical level closest to the corresponding vertical
          ! level of the base emission field
          ! ------------------------------------------------------------

          ! Get lower and upper vertical index
          CALL GetVertIndx ( HcoState, BaseDct, &
                             LevDct1, LevDct2, I, J, LowLL, UppLL, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ERROR = 1 ! Will cause error
             EXIT
          ENDIF

          ! Loop over all vertical levels of the base field
          DO L = LowLL,UppLL
             ! If the vertical level exceeds the number of available
             ! scale factor levels, use the highest available level.
             IF ( L > ScalLL ) THEN
                TmpLL = ScalLL
             ! Otherwise use the same vertical level index.
             ELSE
                TmpLL = L
             ENDIF

             ! Get scale factor for this grid box. Use same uniform
             ! value if it's a scalar field
             IF ( ScalDct%Dta%SpaceDim == 1 ) THEN
                TMPVAL = ScalDct%Dta%V2(tidx)%Val(1,1)
             ELSEIF ( ScalDct%Dta%SpaceDim == 2 ) THEN
                TMPVAL = ScalDct%Dta%V2(tidx)%Val(I,J)
             ELSE
                TMPVAL = ScalDct%Dta%V3(tidx)%Val(I,J,TmpLL)
             ENDIF

             ! Set missing value to one
             IF ( TMPVAL == HCO_MISSVAL ) TMPVAL = 1.0_sp

             ! Eventually apply mask scaling
             IF ( MaskScale /= 1.0_sp ) THEN
                TMPVAL = TMPVAL * MaskScale
             ENDIF

             ! For negative scale factor, proceed according to the
             ! negative value setting specified in the configuration
             ! file (NegFlag = 2: use this value):
             IF ( TMPVAL < 0.0_sp .AND. HcoState%Options%NegFlag /= 2 ) THEN

                ! NegFlag = 1: ignore and show warning
                IF ( HcoState%Options%NegFlag == 1 ) THEN
                   ERROR = -1 ! Will prompt warning
                   CYCLE

                ! Return w/ error otherwise
                ELSE
                   WRITE(*,*) 'Negative scale factor at ',I,J,TmpLL,tidx,&
                              ': ', TRIM(ScalDct%cName), TMPVAL
                   ERROR = 1 ! Will cause error
                   EXIT
                ENDIF
             ENDIF

             ! -------------------------------------------------------
             ! Apply scale factor in accordance to field operator
             ! -------------------------------------------------------

             ! Oper 1: multiply
             IF ( ScalDct%Oper == 1 ) THEN
                OUTARR_3D(I,J,L) = OUTARR_3D(I,J,L) * TMPVAL

             ! Oper -1: divide
             ELSEIF ( ScalDct%Oper == -1 ) THEN
                ! Ignore zeros to avoid NaN
                IF ( TMPVAL /= 0.0_sp ) THEN
                   OUTARR_3D(I,J,L) = OUTARR_3D(I,J,L) / TMPVAL
                ENDIF

             ! Oper 2: square
             ELSEIF ( ScalDct%Oper == 2 ) THEN
                OUTARR_3D(I,J,L) = OUTARR_3D(I,J,L) * TMPVAL * TMPVAL

             ! Return w/ error otherwise (Oper 3 is only allowed for masks!)
             ELSE
                WRITE(*,*) 'Illegal operator for ', TRIM(ScalDct%cName), ScalDct%Oper
                ERROR = 2 ! Will cause error
                EXIT
             ENDIF

          ENDDO !LL

          ! Verbose mode
          if ( HCO_IsVerb(HcoState%Config%Err,3) .and. i == ix .and. j == iy ) then
             write(MSG,*) 'Scale field ', TRIM(ScalDct%cName)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             write(MSG,*) 'Time slice: ', tIdx
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             write(MSG,*) 'IX, IY: ', IX, IY
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             write(MSG,*) 'Scale factor (IX,IY,L1): ', TMPVAL
             CALL HCO_MSG(HcoState%Config%Err,MSG)
             write(MSG,*) 'Mathematical operation : ', ScalDct%Oper
             CALL HCO_MSG(HcoState%Config%Err,MSG)
!             write(lun,*) 'Updt (IX,IY,L1): ', OUTARR_3D(IX,IY,1)
          endif

       ENDDO !I
       ENDDO !J
!$OMP END PARALLEL DO

       ! error check
       IF ( ERROR > 0 ) THEN
          IF ( ERROR == 1 ) THEN
             MSG = 'Negative scale factor found (aborted): ' // TRIM(ScalDct%cName)
          ELSEIF ( ERROR == 2 ) THEN
             MSG = 'Illegal mathematical operator for scale factor: ' // TRIM(ScalDct%cName)
          ELSEIF ( ERROR == 3 ) THEN
             MSG = 'Encountered negative time index for scale factor: ' // TRIM(ScalDct%cName)
          ELSEIF ( ERROR == 3 ) THEN
             MSG = 'Mask error in ' // TRIM(ScalDct%cName)
          ELSE
             MSG = 'Error when applying scale factor: ' // TRIM(ScalDct%cName)
          ENDIF
          ScalDct => NULL()
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
          RETURN
       ENDIF

       ! eventually prompt warning for negative values
       IF ( ERROR == -1 ) THEN
          MSG = 'Negative scale factor found (ignored): ' // TRIM(ScalDct%cName)
          CALL HCO_WARNING( HcoState%Config%Err, MSG, RC )
       ENDIF

       ! Free pointer
       MaskDct => NULL()

    ENDDO ! N
    ENDIF ! N > 0

    ! Update optional variables
    IF ( PRESENT(UseLL) ) THEN
       UseLL = 1
       IF ( nnLL > 0 ) UseLL = NINT(REAL(TotLL,4)/REAL(nnLL,4))
    ENDIF

    ! Weight output emissions by mask
    OUTARR_3D = OUTARR_3D * MASK

    ! Cleanup and leave w/ success
    ScalDct => NULL()
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE Get_Current_Emissions
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Current_Emissions_b (NOT USED!!)
!
! !DESCRIPTION: Subroutine Get\_Current\_Emissions\_B calculates the current
!  emissions for the specified emission field and passes the result to
!  OUTARR\_3D.
!\\
!\\
!  This subroutine is only called by HCO\_CalcEmis and for fields with a valid
!  species ID, i.e. for base emission fields.
!
! !!! WARNING: this routine is not actively developed any more and may lag
! !!! behind Get\_Current\_Emissions
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Current_Emissions_B( HcoState, BaseDct, &
                                      nI, nJ, nL, OUTARR_3D, MASK, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
    USE HCO_TIDX_MOD,     ONLY : tIDx_GetIndx
    USE HCO_FILEDATA_MOD, ONLY : FileData_ArrIsDefined
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN)    :: nI                  ! # of lons
    INTEGER,         INTENT(IN)    :: nJ                  ! # of lats
    INTEGER,         INTENT(IN)    :: nL                  ! # of levs
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState            ! HEMCO state object
    TYPE(DataCont),  POINTER       :: BaseDct             ! base emission
                                                          !  container
    REAL(hp),        INTENT(INOUT) :: OUTARR_3D(nI,nJ,nL) ! output array
    REAL(hp),        INTENT(INOUT) :: MASK     (nI,nJ,nL) ! mask array
    INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  25 Aug 2012 - C. Keller   - Initial Version
!  09 Nov 2012 - C. Keller   - MASK update. Masks are now treated
!                              separately so that multiple masks can be
!                              added.
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX header
!  07 Sep 2014 - C. Keller   - Mask update. Now set mask to zero as soon as
!                              on of the applied masks is zero.
!  02 Mar 2015 - C. Keller   - Now check for missing values
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DataCont), POINTER :: ScalDct
    TYPE(DataCont), POINTER :: MaskDct
    REAL(sp)                :: TMPVAL, MaskScale
    INTEGER                 :: tIdx, IDX
    INTEGER                 :: I, J, L, N
    INTEGER                 :: LowLL, UppLL, ScalLL, TmpLL
    INTEGER                 :: IJFILLED
    INTEGER                 :: ERROR
    CHARACTER(LEN=255)      :: MSG, LOC
    LOGICAL                 :: MaskFractions

    ! testing only
    INTEGER                 :: IX, IY
    LOGICAL                 :: verb

    !=================================================================
    ! GET_CURRENT_EMISSIONS_B begins here
    !=================================================================

    ! Initialize
    ScalDct => NULL()
    MaskDct => NULL()

    ! Enter
    CALL HCO_ENTER(HcoState%Config%Err,'GET_CURRENT_EMISSIONS_B', RC )
    IF(RC /= HCO_SUCCESS) RETURN

    ! testing only
    verb = HCO_IsVerb(HcoState%Config%Err,1)
    IX = 60 !40 !19 43 61
    IY = 32 !36 !33 26 37

    ! Check if field data is defined
    IF ( .NOT. FileData_ArrIsDefined(BaseDct%Dta) ) THEN
       MSG = 'Array not defined: ' // TRIM(BaseDct%cName)
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF

    ! Testing only:
    IF ( verb ) THEN
       write(MSG,*) '--> GET EMISSIONS FOR ', TRIM(BaseDct%cName)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    ! Initialize mask values
    MASK(:,:,:)  = 1.0_hp
    MaskFractions = HcoState%Options%MaskFractions

    ! Initialize ERROR. Will be set to 1 if error occurs below
    ERROR = 0

    ! Loop over all grid boxes
!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT( SHARED )                                                &
!$OMP PRIVATE( I, J, LowLL, UppLL, tIdx, IJFILLED, L                 ) &
!$OMP PRIVATE( TMPVAL, N, IDX, ScalDct, ScalLL, tmpLL, MaskScale     ) &
!$OMP SCHEDULE( DYNAMIC )
    DO J = 1, nJ
    DO I = 1, nI

       ! -------------------------------------------------------------
       ! Set base emissions
       ! -------------------------------------------------------------

       ! Get vertical extension of base emission array.
       ! Unlike the output array OUTARR_3D, the data containers do not
       ! necessarily extent over the entire troposphere but only cover
       ! the effectively filled vertical levels. For most inventories,
       ! this is only the first model level.
       IF ( BaseDct%Dta%SpaceDim==3 ) THEN
          LowLL = 1
          UppLL = SIZE(BaseDct%Dta%V3(1)%Val,3)
       ELSE
          !LowLL = BaseDct%Dta%Lev2D
          !UppLL = BaseDct%Dta%Lev2D
          LowLL = 1
          UppLL = 1
       ENDIF

       ! Precalculate timeslice index. The data containers can
       ! carry 2D/3D arrays for multiple time steps (i.e. for
       ! every hour of the day), stored in a vector.
       ! tIdxVec contains the vector index to be used at the current
       ! datetime. This parameter may vary with longitude due to time
       ! zone shifts!
       tIDx = tIDx_GetIndx( HcoState, BaseDct%Dta, I, J )
       IF ( tIDx < 0 ) THEN
          write(MSG,*) 'Cannot get time slice index at location ',I,J,&
                       ': ', TRIM(BaseDct%cName)
          ERROR = 3
          EXIT
       ENDIF

       ! # of levels w/ defined emissions
       IJFILLED = 0

       ! Loop over all levels
       DO L = LowLL, UppLL

          ! Get base value. Use uniform value if scalar field.
          IF ( BaseDct%Dta%SpaceDim == 1 ) THEN
             TMPVAL = BaseDct%Dta%V2(tIDx)%Val(1,1)
          ELSEIF ( BaseDct%Dta%SpaceDim == 2 ) THEN
             TMPVAL = BaseDct%Dta%V2(tIDx)%Val(I,J)
          ELSE
             TMPVAL = BaseDct%Dta%V3(tIDx)%Val(I,J,L)
          ENDIF

          ! Check for missing value
          IF ( TMPVAL == HCO_MISSVAL ) THEN
             OUTARR_3D(I,J,L) = 0.0_hp
             MASK(I,J,:)      = 0.0_hp

          ! Pass base value to output array
          ELSE
             OUTARR_3D(I,J,L) = TMPVAL
          ENDIF

          ! Update IJFILLED
          IJFILLED = IJFILLED + 1

       ENDDO !L

       ! -------------------------------------------------------------
       ! Apply scale factors
       ! The container IDs of all scale factors associated with this base
       ! container are stored in vector Scal_cID.
       ! -------------------------------------------------------------

       ! Loop over maximum number of scale factors
       IF ( BaseDct%nScalID > 0 ) THEN
       DO N = 1, BaseDct%nScalID

          ! Get the scale factor container ID for the current slot
          IDX = BaseDct%Scal_cID(N)

          ! Point to emission container with the given container ID
          CALL Pnt2DataCont( HcoState, IDX, ScalDct, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ERROR = 4
             EXIT
          ENDIF

          ! Scale field cannot be a base field
          IF ( (ScalDct%DctType == HCO_DCTTYPE_BASE) ) THEN
             ERROR = 4
             EXIT
          ENDIF

          ! Skip this scale factor if no data defined. This is possible
          ! if scale factors are only defined for a given time range and
          ! the simulation datetime is outside of this range.
          IF ( .NOT. FileData_ArrIsDefined(ScalDct%Dta) ) THEN
             MSG = 'Array not defined: ' // TRIM(ScalDct%cName)
             CALL HCO_WARNING( HcoState%Config%Err, MSG, RC )
             CYCLE
          ENDIF

          ! Check if there is a mask field associated with this scale
          ! factor. In this case, get a pointer to the corresponding
          ! mask field and evaluate scale factors only inside the mask
          ! region.
          IF ( ASSOCIATED(ScalDct%Scal_cID) ) THEN
             CALL Pnt2DataCont( HcoState, ScalDct%Scal_cID(1), MaskDct, RC )
             IF ( RC /= HCO_SUCCESS ) THEN
                ERROR = 5
                EXIT
             ENDIF

             ! Must be mask field
             IF ( MaskDct%DctType /= HCO_DCTTYPE_MASK ) THEN
                MSG = 'Invalid mask for scale factor: '//TRIM(ScalDct%cName)
                MSG = TRIM(MSG) // '; mask: '//TRIM(MaskDct%cName)
                CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
                ERROR = 5
                EXIT
             ENDIF

             ! Get mask value
             CALL GetMaskVal( ScalDct, I, J, TMPVAL, MaskFractions, RC )
             IF ( RC /= HCO_SUCCESS ) THEN
                ERROR = 6
                EXIT
             ENDIF

          ENDIF

          ! Get vertical extension of this scale factor array.
          IF( (ScalDct%Dta%SpaceDim<=2) ) THEN
             ScalLL = 1
          ELSE
             ScalLL = SIZE(ScalDct%Dta%V3(1)%Val,3)
          ENDIF

          ! Get current time index
          tIDx = tIDx_GetIndx( HcoState, ScalDct%Dta, I, J )
          IF ( tIDx < 0 ) THEN
             write(MSG,*) 'Cannot get time slice index at location ',I,J,&
                          ': ', TRIM(ScalDct%cName)
             ERROR = 3
             EXIT
          ENDIF

          ! ------------------------------------------------------------
          ! Check if this is a mask. If so, add mask values to the MASK
          ! array. For now, we assume masks to be binary, i.e. 0 or 1.
          ! We may want to change that in future to also support values
          ! in between. This is especially important when regridding
          ! high resolution masks onto coarser grids!
          ! ------------------------------------------------------------
          IF ( ScalDct%DctType == HCO_DCTTYPE_MASK ) THEN

             ! Get mask value
             CALL GetMaskVal( ScalDct, I, J, TMPVAL, MaskFractions, RC )
             IF ( RC /= HCO_SUCCESS ) THEN
                ERROR = 6
                EXIT
             ENDIF

             ! Pass to mask
             MASK(I,J,:) = MASK(I,J,:) * TMPVAL

             ! testing only
             if ( verb .and. i == ix .and. j == iy ) then
                write(*,*) 'Mask field ', TRIM(ScalDct%cName),   &
                     ' found and added to temporary mask.'
             ENDIF

             ! Advance to next scale factor
             CYCLE
          ENDIF! DctType=MASK

          ! ------------------------------------------------------------
          ! For non-mask fields, apply scale factors to all levels
          ! of the base field individually. If the scale factor
          ! field has more than one vertical level, use the
          ! vertical level closest to the corresponding vertical
          ! level in the base emission field
          ! ------------------------------------------------------------

          ! Loop over all vertical levels of the base field
          DO L = LowLL,UppLL
             ! If the vertical level exceeds the number of available
             ! scale factor levels, use the highest available level.
             IF ( L > ScalLL ) THEN
                TmpLL = ScalLL
             ! Otherwise use the same vertical level index.
             ELSE
                TmpLL = L
             ENDIF

             ! Get scale factor for this grid box. Use same uniform
             ! value if it's a scalar field
             IF ( ScalDct%Dta%SpaceDim == 1 ) THEN
                TMPVAL = ScalDct%Dta%V2(tidx)%Val(1,1)
             ELSEIF ( ScalDct%Dta%SpaceDim == 2 ) THEN
                TMPVAL = ScalDct%Dta%V2(tidx)%Val(I,J)
             ELSE
                TMPVAL = ScalDct%Dta%V3(tidx)%Val(I,J,TmpLL)
             ENDIF

             ! Check for missing value
             IF ( TMPVAL == HCO_MISSVAL ) TMPVAL = 1.0_sp

             ! For negative scale factor, proceed according to the
             ! negative value setting specified in the configuration
             ! file (NegFlag = 2: use this value):
             IF ( TMPVAL < 0.0_sp .AND. HcoState%Options%NegFlag /= 2 ) THEN

                ! NegFlag = 1: ignore and show warning
                IF ( HcoState%Options%NegFlag == 1 ) THEN
                   ERROR = -1 ! Will prompt warning
                   CYCLE

                ! Return w/ error otherwise
                ELSE
                   WRITE(*,*) 'Negative scale factor at ',I,J,TmpLL,tidx,&
                              ': ', TRIM(ScalDct%cName), TMPVAL
                   ERROR = 1 ! Will cause error
                   EXIT
                ENDIF
             ENDIF

             ! -------------------------------------------------------
             ! Apply scale factor according to field operator
             ! -------------------------------------------------------

             ! Oper 1: multiply
             IF ( ScalDct%Oper == 1 ) THEN
                OUTARR_3D(I,J,L) = OUTARR_3D(I,J,L) * TMPVAL

             ! Oper -1: divide
             ELSEIF ( ScalDct%Oper == -1 ) THEN
                ! Ignore zeros to avoid NaN
                IF ( TMPVAL /= 0.0_sp ) THEN
                   OUTARR_3D(I,J,L) = OUTARR_3D(I,J,L) / TMPVAL
                ENDIF

             ! Oper 2: square
             ELSEIF ( ScalDct%Oper == 2 ) THEN
                OUTARR_3D(I,J,L) = OUTARR_3D(I,J,L) * TMPVAL * TMPVAL

             ! Return w/ error otherwise (Oper 3 only allowed for masks!)
             ELSE
                MSG = 'Illegal data operator: ' // TRIM(ScalDct%cName)
                CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
                ERROR = 2
                EXIT
             ENDIF
          ENDDO !LL
       ENDDO ! N
       ENDIF ! N > 0

       ! ----------------------------
       ! Masks
       ! ----------------------------

       ! Apply the mask. Make sure that emissions become negative
       ! outside the mask region. This is to make sure that these
       ! grid boxes will be ignored when calculating the final
       ! emissions.
       WHERE ( MASK(I,J,:) == 0 )
          OUTARR_3D(I,J,:) = 0.0_hp
       ENDWHERE

    ENDDO !I
    ENDDO !J
!$OMP END PARALLEL DO

    ! Error check
    IF ( ERROR > 0 ) THEN
       IF ( ERROR == 1 ) THEN
          MSG = 'Negative scale factor found (aborted): ' // TRIM(ScalDct%cName)
       ELSEIF ( ERROR == 2 ) THEN
          MSG = 'Illegal mathematical operator for scale factor: ' // TRIM(ScalDct%cName)
       ELSEIF ( ERROR == 3 ) THEN
          MSG = 'Encountered negative time index for scale factor: ' // TRIM(ScalDct%cName)
       ELSE
          MSG = 'Error when applying scale factor: ' // TRIM(ScalDct%cName)
       ENDIF
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
       ScalDct => NULL()
       RETURN
    ENDIF

    ! eventually prompt warning for negative values
    IF ( ERROR == -1 ) THEN
       MSG = 'Negative scale factor found (ignored): ' // TRIM(ScalDct%cName)
       CALL HCO_WARNING( HcoState%Config%Err, MSG, RC )
    ENDIF

    ! Leave
    ScalDct => NULL()
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE Get_Current_Emissions_B
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EvalFld_3D
!
! !DESCRIPTION: Subroutine HCO\_EvalFld\_3D returns the 3D data field belonging
!  to the emissions list data container with field name 'cName'. The returned
!  data field is the completely evaluated field, e.g. the base field multiplied
!  by all scale factors and with all masking being applied (as specified in the
!  HEMCO configuration file). This distinguished this routine from HCO\_GetPtr
!  in hco\_emislist\_mod.F90, which returns a reference to the unevaluated data
!  field.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EvalFld_3D( HcoState, cName, Arr3D, RC, FOUND )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
    USE HCO_DATACONT_MOD, ONLY : ListCont_Find
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   )  :: cName
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState     ! HEMCO state object
    REAL(hp),         INTENT(INOUT)  :: Arr3D(:,:,:) ! 3D array
    INTEGER,          INTENT(INOUT)  :: RC           ! Return code
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(  OUT), OPTIONAL  :: FOUND
!
! !REVISION HISTORY:
!  11 May 2015 - C. Keller   - Initial Version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                 :: FND
    INTEGER                 :: AS, nI, nJ, nL, FLAG

    ! Arrays
    REAL(hp), ALLOCATABLE   :: Mask(:,:,:)

    ! Working pointers: list and data container
    TYPE(ListCont), POINTER :: Lct

    ! For error handling & verbose mode
    CHARACTER(LEN=255)  :: MSG
    CHARACTER(LEN=255)  :: LOC = "HCO_EvalFld_3d (HCO_calc_mod.F90)"

    !=================================================================
    ! HCO_EvalFld_3D begins here!
    !=================================================================

    ! Init
    RC    = HCO_SUCCESS
    Lct   => NULL()
    IF ( PRESENT(FOUND) ) FOUND = .FALSE.

    ! Search for base container
    CALL ListCont_Find ( HcoState%EmisList, TRIM(cName), FND, Lct )
    IF ( PRESENT(FOUND) ) FOUND = FND

    ! If not found, return here
    IF ( .NOT. FND ) THEN
       IF ( PRESENT(FOUND) ) THEN
          RETURN
       ELSE
          MSG = 'Cannot find in EmisList: ' // TRIM(cName)
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDIF

    ! Init
    Arr3D = 0.0_hp

    ! Define output dimensions
    nI = SIZE(Arr3D,1)
    nJ = SIZE(Arr3D,2)
    nL = SIZE(Arr3D,3)

    ! Sanity check: horizontal grid dimensions are expected to be on HEMCO grid
    IF ( nI /= HcoState%NX .OR. nJ /= HcoState%nY ) THEN
       WRITE(MSG,*) "Horizontal dimension error: ", TRIM(cName), nI, nJ
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Make sure mask array is defined
    ALLOCATE(MASK(nI,nJ,nL),STAT=AS)
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate MASK', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Calculate emissions for base container
    CALL GET_CURRENT_EMISSIONS( HcoState, Lct%Dct, &
                                nI, nJ, nL, Arr3D,    Mask, RC  )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! All done
    IF (ALLOCATED(MASK) ) DEALLOCATE(MASK)
    Lct => NULL()

  END SUBROUTINE HCO_EvalFld_3D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EvalFld_2D
!
! !DESCRIPTION: Subroutine HCO\_EvalFld\_2D returns the 2D data field belonging
!  to the emissions list data container with field name 'cName'. The returned
!  data field is the completely evaluated field, e.g. the base field multiplied
!  by all scale factors and with all masking being applied (as specified in the
!  HEMCO configuration file). This distinguished this routine from HCO\_GetPtr
!  in hco\_emislist\_mod.F90, which returns a reference to the unevaluated data
!  field.
!\\
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EvalFld_2D( HcoState, cName, Arr2D, RC, FOUND )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
    USE HCO_DATACONT_MOD, ONLY : ListCont_Find
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   )  :: cName
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState     ! HEMCO state object
    REAL(hp),         INTENT(INOUT)  :: Arr2D(:,:)   ! 2D array
    INTEGER,          INTENT(INOUT)  :: RC           ! Return code
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(  OUT), OPTIONAL  :: FOUND
!
! !REVISION HISTORY:
!  11 May 2015 - C. Keller   - Initial Version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  01 Nov 2016 - C. Keller   - Added error trap for UseLL
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                 :: FND
    INTEGER                 :: AS, nI, nJ, nL, UseLL, FLAG

    ! Arrays
    REAL(hp), ALLOCATABLE   :: Mask (:,:,:)
    REAL(hp), ALLOCATABLE   :: Arr3D(:,:,:)

    ! Working pointers: list and data container
    TYPE(ListCont), POINTER :: Lct

    ! For error handling & verbose mode
    CHARACTER(LEN=255)  :: MSG
    CHARACTER(LEN=255)  :: LOC = "HCO_EvalFld_2d (HCO_calc_mod.F90)"

    !=================================================================
    ! HCO_EvalFld_2D begins here!
    !=================================================================

    ! Init
    RC    = HCO_SUCCESS
    Lct   => NULL()
    IF ( PRESENT(FOUND) ) FOUND = .FALSE.

    ! Search for base container
    CALL ListCont_Find ( HcoState%EmisList, TRIM(cName), FND, Lct )
    IF ( PRESENT(FOUND) ) FOUND = FND

    ! If not found, return here
    IF ( .NOT. FND ) THEN
       IF ( PRESENT(FOUND) ) THEN
          RETURN
       ELSE
          MSG = 'Cannot find in EmisList: ' // TRIM(cName)
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDIF

    ! Init Arr2D
    Arr2D = 0.0_hp

    ! Define output dimensions
    nI = SIZE(Arr2D,1)
    nJ = SIZE(Arr2D,2)
    nL = 1

    ! Sanity check: horizontal grid dimensions are expected to be on HEMCO grid
    IF ( nI /= HcoState%NX .OR. nJ /= HcoState%nY ) THEN
       WRITE(MSG,*) "Horizontal dimension error: ", TRIM(cName), nI, nJ
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Make sure mask array is defined
    ALLOCATE(MASK(nI,nJ,nL),Arr3D(nI,nJ,nL),STAT=AS)
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate MASK', RC, THISLOC=LOC )
       RETURN
    ENDIF
    Arr3D = 0.0_hp
    Mask  = 0.0_hp

    ! Calculate emissions for base container
    CALL GET_CURRENT_EMISSIONS( HcoState, Lct%Dct, &
                                nI, nJ, nL, Arr3D, Mask, RC, UseLL=UseLL )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Place 3D array into 2D array. UseLL returns the vertical level into which
    ! emissions have been added within GET_CURRENT_EMISSIONS. This should be
    ! level 1 for most cases but it can be another level if specified so. Return
    ! a warning if level is not 1 (ckeller, 11/1/16).
    UseLL = MIN( MAX(useLL,1), SIZE(Arr3D,3) )
    IF ( UseLL /= 1 ) THEN
       WRITE(MSG,*) "2D data was emitted above surface - this information might be lost: " , TRIM(cName), UseLL
       CALL HCO_WARNING( HcoState%Config%Err, MSG, RC, THISLOC=LOC, WARNLEV=2 )
    ENDIF

    ! Pass 3D data to 2D array
    Arr2D(:,:) = Arr3D(:,:,UseLL)

    ! All done
    IF (ALLOCATED(MASK ) ) DEALLOCATE(MASK )
    IF (ALLOCATED(Arr3D) ) DEALLOCATE(Arr3D)
    Lct => NULL()

  END SUBROUTINE HCO_EvalFld_2D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetMaskVal
!
! !DESCRIPTION: Subroutine GetMaskVal is a helper routine to get the mask
!  value at a given location.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetMaskVal ( Dct, I, J, MaskVal, Fractions, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN   ) :: I                   ! # of lons
    INTEGER,         INTENT(IN   ) :: J                   ! # of lats
    LOGICAL,         INTENT(IN   ) :: Fractions           ! Use fractions?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DataCont),  POINTER       :: Dct                 ! Mask container
    REAL(sp),        INTENT(INOUT) :: MaskVal
    INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  09 Apr 2015 - C. Keller   - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    !=================================================================
    ! GetMaskVal begins here
    !=================================================================

    ! Mask value over this grid box
    MaskVal = Dct%Dta%V2(1)%Val(I,J)

    ! Negative mask values are treated as zero (exclude).
    IF ( (MaskVal <= 0.0_sp) .OR. (MaskVal == HCO_MISSVAL) ) THEN
       MaskVal = 0.0_sp
    ELSEIF ( MaskVal > 1.0_sp ) THEN
       MaskVal = 1.0_sp
    ENDIF

    ! For operator set to 3, mirror value
    ! MaskVal=1 becomes 0 and MaskVal=0/missing becomes 1
    IF ( Dct%Oper == 3 ) THEN
       IF ( (MaskVal == 0.0_sp) .OR. (MaskVal == HCO_MISSVAL) ) THEN
          MaskVal = 1.0_sp
       ELSEIF ( MaskVal == 1.0_sp ) THEN
          MaskVal = 1.0_sp - MaskVal
       ENDIF
    ENDIF

    ! Treat as binary?
    IF ( .NOT. Fractions ) THEN
       IF ( MaskVal < MASK_THRESHOLD ) THEN
          MaskVal = 0.0_sp
       ELSE
          MaskVal = 1.0_sp
       ENDIF
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE GetMaskVal
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_MaskFld
!
! !DESCRIPTION: Subroutine HCO\_MaskFld is a helper routine to get the mask
! field with the given name. The returned mask field is fully evaluated,
! e.g. the data operation flag associated with this mask field is already
! taken into account. For instance, if the data operator of a mask field is
! set to 3, the returned array contains already the mirrored mask values.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_MaskFld ( HcoState, MaskName, Mask, RC, FOUND )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
    USE HCO_DATACONT_MOD, ONLY : ListCont_Find
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_STATE), POINTER                 :: HcoState
    CHARACTER(LEN=*),INTENT(IN   )           :: MaskName
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp),        INTENT(INOUT)           :: Mask(:,:)
    INTEGER,         INTENT(INOUT)           :: RC
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,         INTENT(  OUT), OPTIONAL :: FOUND
!
! !REVISION HISTORY:
!  11 Jun 2015 - C. Keller   - Initial Version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                 :: I, J, FLAG

    LOGICAL                 :: FND, ERR
    LOGICAL                 :: Fractions

    TYPE(ListCont), POINTER :: MaskLct

    CHARACTER(LEN=255)      :: MSG
    CHARACTER(LEN=255)      :: LOC = 'HCO_MaskFld (hco_calc_mod.F90)'

    !=================================================================
    ! HCO_MaskFld begins here
    !=================================================================

    ! Nullify
    MaskLct  => NULL()

    ! Init: default is mask value of 1
    MASK = 1.0_sp
    ERR  = .FALSE.
    FND  = .FALSE.

    ! Search for mask field within EmisList
    CALL ListCont_Find ( HcoState%EmisList, TRIM(MaskName), FND, MaskLct )

    IF ( .NOT. FND .AND. .NOT. PRESENT(FOUND) ) THEN
       MSG = 'Cannot find mask field ' // TRIM(MaskName)
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='!')
       MSG = 'Make sure this field is listed in the mask section '  // &
           'of the HEMCO configuration file. You may also need to ' // &
           'set the optional attribute `ReadAlways` to `yes`, e.g.'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       MSG = '5000 TESTMASK     -140/10/-40/90 - - - xy 1 1 -140/10/-40/90 yes'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       CALL HCO_ERROR ( HcoState%Config%Err, &
                        'Error reading mask '//TRIM(MaskName), RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( PRESENT(FOUND) ) FOUND = FND

    ! Do only if found
    IF ( FND ) THEN

       ! Use mask fractions?
       Fractions = HcoState%Options%MaskFractions

       ! Make sure mask array has correct dimensions
       IF ( SIZE(MASK,1) /= HcoState%NX .OR. SIZE(MASK,2) /= HcoState%NY ) THEN
          WRITE(MSG,*) 'Input mask array has wrong dimensions. Must be ', &
             HcoState%NX, HcoState%NY, ' but found ', SIZE(MASK,1), SIZE(MASK,2)
          CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Do for every grid box
!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT( SHARED )                                                &
!$OMP PRIVATE( I, J                                                  ) &
!$OMP SCHEDULE( DYNAMIC )
       DO J = 1, HcoState%NY
       DO I = 1, HcoState%NX
          CALL GetMaskVal( MaskLct%Dct, I, J, Mask(I,J), Fractions, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ERR = .TRUE.
             EXIT
          ENDIF
       ENDDO
       ENDDO
!$OMP END PARALLEL DO

       ! Error check
       IF ( ERR ) THEN
          MSG = 'Error in GetMaskVal'
          CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

    ENDIF

    ! Free pointer
    MaskLct  => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_MaskFld
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetVertIndx
!
! !DESCRIPTION: Subroutine GetVertIndx is a helper routine to get the vertical
!  index range of the given data field.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetVertIndx ( HcoState, Dct, &
                           LevDct1, LevDct2, I, J, LowLL, UppLL, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER          :: HcoState    ! HEMCO state object
    TYPE(DataCont),  POINTER          :: LevDct1     ! Level index 1 container
    TYPE(DataCont),  POINTER          :: LevDct2     ! Level index 2 container
    INTEGER,  INTENT(IN   )           :: I           ! lon index
    INTEGER,  INTENT(IN   )           :: J           ! lat index
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DataCont),  POINTER          :: Dct         ! Mask container
    INTEGER,  INTENT(INOUT)           :: LowLL       ! lower level index
    INTEGER,  INTENT(INOUT)           :: UppLL       ! upper level index
    INTEGER,  INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  06 May 2016 - C. Keller   - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: EmisLUnit
    REAL(hp) :: EmisL

    !=================================================================
    ! GetVertIndx begins here
    !=================================================================

    ! Get vertical extension of base emission array.
    ! Unlike the output array OUTARR_3D, the data containers do not
    ! necessarily extent over the entire troposphere but only cover
    ! the effectively filled vertical levels. For most inventories,
    ! this is only the first model level.
    IF ( Dct%Dta%SpaceDim==3 ) THEN
       LowLL = 1
       UppLL = SIZE(Dct%Dta%V3(1)%Val,3)

    ! For 2D field, check if it shall be spread out over multiple
    ! levels. Possible to go from PBL to max. specified level.
    ELSE
       ! Lower level
       ! --> Check if scale factor is used to determine lower and/or
       !     upper level
       IF ( ASSOCIATED(LevDct1) ) THEN
          EmisL = GetEmisL( HcoState, LevDct1, I, J )
          IF ( EmisL < 0.0_hp ) THEN
             RC = HCO_FAIL
             RETURN
          ENDIF
          EmisLUnit = GetEmisLUnit( HcoState, LevDct1 )
          IF ( EmisLUnit < 0 ) THEN
             RC = HCO_FAIL
             RETURN
          ENDIF
       ELSE
          EmisL     = Dct%Dta%EmisL1
          EmisLUnit = Dct%Dta%EmisL1Unit
       ENDIF
       CALL GetIdx( HcoState, I, J, EmisL, EmisLUnit, LowLL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Upper level
       IF ( ASSOCIATED(LevDct2) ) THEN
          EmisL = GetEmisL( HcoState, LevDct2, I, J )
          IF ( EmisL < 0.0_hp ) THEN
             RC = HCO_FAIL
             RETURN
          ENDIF
          EmisLUnit = GetEmisLUnit( HcoState, LevDct2 )
          IF ( EmisLUnit < 0 ) THEN
             RC = HCO_FAIL
             RETURN
          ENDIF
       ELSE
          EmisL     = Dct%Dta%EmisL2
          EmisLUnit = Dct%Dta%EmisL2Unit
       ENDIF
       CALL GetIdx( HcoState, I, J, EmisL, EmisLUnit, UppLL, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Upper level must not be lower than lower level
       UppLL = MAX(LowLL, UppLL)
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE GetVertIndx
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: GetEmisL
!
! !DESCRIPTION: Returns the emission level read from a scale factor.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GetEmisL ( HcoState, LevDct, I, J ) RESULT ( EmisL )
!
! !USES:
!
    USE HCO_TYPES_MOD
    USE HCO_STATE_MOD,    ONLY : HCO_State
    USE HCO_tIdx_MOD,     ONLY : tIDx_GetIndx
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState       ! HEMCO state object
    TYPE(DataCont),  POINTER        :: LevDct         ! Level index 1 container
    INTEGER,         INTENT(IN   )  :: I, J           ! horizontal index
!
! !RETURN VALUE:
!
    REAL(hp)                        :: EmisL
!
! !REVISION HISTORY:
!  26 Jan 2018 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: levtidx

    !=================================================================
    ! GetEmisL begins here
    !=================================================================
    levtidx = tIDx_GetIndx( HcoState, LevDct%Dta, I, J )
    IF ( levtidx <= 0 ) THEN
       WRITE(*,*)' Cannot get time slice for field '//&
       TRIM(LevDct%cName)//': GetEmisL (hco_calc_mod.F90)'
       EmisL = -1.0
       RETURN
    ENDIF

    IF ( LevDct%Dta%SpaceDim == 1 ) THEN
       EmisL = LevDct%Dta%V2(levtidx)%Val(1,1)
    ELSEIF ( LevDct%Dta%SpaceDim == 2 ) THEN
       EmisL = LevDct%Dta%V2(levtidx)%Val(I,J)
    ELSEIF ( LevDct%Dta%SpaceDim == 3 ) THEN
       EmisL = LevDct%Dta%V3(levtidx)%Val(I,J,1)
    ENDIF

    IF ( EmisL == HCO_MISSVAL ) EmisL = 0.0_hp

END FUNCTION GetEmisL
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: GetEmisLUnit
!
! !DESCRIPTION: Returns the emission level unit read from a scale factor.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GetEmisLUnit ( HcoState, LevDct ) RESULT( EmisLUnit )
!
! !USES:
!
    USE HCO_TYPES_MOD
    USE HCO_STATE_MOD,    ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState       ! HEMCO state object
    TYPE(DataCont),  POINTER        :: LevDct         ! Level index 1 container
!
! !RETURN VALUE:
!
    INTEGER                         :: EmisLUnit
!
! !REVISION HISTORY:
!  26 Jan 2018 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !=================================================================
    ! GetEmisLUnit begins here
    !=================================================================

    ! For now, only meters are supported
    EmisLUnit = HCO_EMISL_M

    ! Dummy check that units on field are actually in meters
    IF ( TRIM(LevDct%Dta%OrigUnit) /= 'm' .AND. &
         TRIM(LevDct%Dta%OrigUnit) /= '1'        ) THEN
       WRITE(*,*) TRIM(LevDct%cName)// &
       ' must have units of `m`, instead found '//&
       TRIM(LevDct%Dta%OrigUnit)//': GetEmisLUnit (hco_calc_mod.F90)'
       EmisLUnit = -1
    ENDIF

END FUNCTION GetEmisLUnit
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetIdx
!
! !DESCRIPTION: Subroutine GetIdx is a helper routine to return the vertical
!  level index for a given altitude. The altitude can be provided in level
!  coordinates, in units of meters or as the 'PBL mixing height'.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetIdx( HcoState, I, J, alt, altu, lidx, RC )
!
! !USES:
!
    USE HCO_TYPES_MOD
    USE HCO_STATE_MOD,   ONLY : HCO_STATE
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState       ! HEMCO state object
    INTEGER,         INTENT(IN   )  :: I, J           ! horizontal index
    INTEGER,         INTENT(IN   )  :: altu           ! altitude unit
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(  OUT)  :: lidx           ! level index
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(hp),        INTENT(INOUT)  :: alt            ! altitude
    INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  09 May 2016 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                 :: L
    REAL(hp)                :: altb, altt
    CHARACTER(LEN=255)      :: MSG
    CHARACTER(LEN=255)      :: LOC = 'GetIdx (hco_geotools_mod.F90)'

    !=================================================================
    ! HCO_GetVertIndx begins here
    !=================================================================

    ! Init
    RC = HCO_SUCCESS

    ! Simple case: data is already on level unit
    IF ( altu == HCO_EMISL_LEV ) THEN
       lidx = INT(alt)

    ELSEIF ( altu == HCO_EMISL_M .OR. altu == HCO_EMISL_PBL ) THEN

       ! Eventually get altitude from PBL height
       IF ( altu == HCO_EMISL_PBL ) THEN
          IF ( .NOT. ASSOCIATED(HcoState%Grid%PBLHEIGHT%Val) ) THEN
             MSG = 'PBL field missing in HEMCO state'
             CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
          alt = HcoState%Grid%PBLHEIGHT%Val(I,J)
       ENDIF

       ! Special case of negative height
       IF ( alt <= 0.0_hp ) THEN
          lidx = 1
          RETURN
       ENDIF

       ! Error check
       IF ( .NOT. ASSOCIATED(HcoState%Grid%BXHEIGHT_M%Val) ) THEN
          MSG = 'Boxheight missing in HEMCO state'
          CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Loop over data until we are within desired level
       altt = 0.0_hp
       altb = 0.0_hp
       lidx = -1
       DO L = 1, HcoState%NZ
          altt = altb + HcoState%Grid%BXHEIGHT_M%Val(I,J,L)
          IF ( alt >= altb .AND. alt < altt ) THEN
             lidx = L
             RETURN
          ENDIF
          altb = altt
       ENDDO

       ! If altitude is above maximum level
       IF ( lidx == -1 .AND. alt >= altt ) THEN
          lidx = HcoState%NZ
          WRITE(MSG,*)  'Level is above max. grid box level - use top level ', alt
          CALL HCO_WARNING ( HcoState%Config%Err, MSG, RC, THISLOC=LOC, WARNLEV=2 )
          RETURN
       ENDIF

    ELSE
       MSG = 'Illegal altitude unit'
       CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE GetIdx
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetDilFact
!
! !DESCRIPTION: Subroutine GetDilFact returns the vertical dilution factor,
! that is the factor that is to be applied to distribute emissions into
! multiple vertical levels. If grid box height information are available,
! these are used to compute the distribution factor. Otherwise, equal weight
! is given to all vertical levels.
!\\
!\\
! !TODO: Dilution factors are currently only weighted by grid box heights
! (if these information are available) but any pressure information are
! ignored.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetDilFact ( HcoState, EmisL1, EmisL1Unit, &
                          EmisL2, EmisL2Unit, I, J, L, LowLL, UppLL, &
                          DilFact, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER          :: HcoState    ! HEMCO state object
    INTEGER,  INTENT(IN   )           :: I           ! lon index
    INTEGER,  INTENT(IN   )           :: J           ! lat index
    INTEGER,  INTENT(IN   )           :: L           ! lev index
    INTEGER,  INTENT(IN   )           :: LowLL       ! lower level index
    INTEGER,  INTENT(IN   )           :: UppLL       ! upper level index
!
! !OUTPUT PARAMETERS:
!
    REAL(hp), INTENT(  OUT)           :: DilFact     ! Dilution factor
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(hp), INTENT(INOUT)           :: EmisL1
    INTEGER,  INTENT(INOUT)           :: EmisL1Unit
    REAL(hp), INTENT(INOUT)           :: EmisL2
    INTEGER,  INTENT(INOUT)           :: EmisL2Unit
    INTEGER,  INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  06 May 2016 - C. Keller   - Initial Version
!  16 Jul 2018 - C. Keller   - Bug fix: get PBL height in m, properly compute
!                              fractions in lowest / highest level.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER    :: L1
    REAL(hp)   :: h1, h2, lh, dh

    !=================================================================
    ! GetDilFact begins here
    !=================================================================

    ! Init
    DilFact = 1.0_hp
    RC = HCO_SUCCESS

    ! Nothing to do if it's only one level
    IF ( LowLL == UppLL ) RETURN

    ! Compute 'accurate' dilution factor if boxheights are available
    IF ( HcoState%Options%VertWeight .AND. &
         ASSOCIATED( HcoState%Grid%BXHEIGHT_M%Val ) ) THEN

       ! Height of grid box of interest (in m)
       dh = HcoState%Grid%BXHEIGHT_M%Val(I,J,L)

       ! Get height of bottom level LowLL (in m)
       IF ( EmisL1Unit == HCO_EMISL_M ) THEN
          h1 = EmisL1
       ELSEIF ( EmisL1Unit == HCO_EMISL_PBL ) THEN
          h1 = HcoState%Grid%PBLHEIGHT%Val(I,J)
       ELSE
          IF ( LowLL > 1 ) THEN
             h1 = SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:(LowLL-1)))
          ELSE
             h1 = 0.0_hp
          ENDIF
       ENDIF

       ! Get height of top level UppLL (in m)
       IF ( EmisL2Unit == HCO_EMISL_M ) THEN
          h2 = EmisL2
       ELSEIF ( EmisL2Unit == HCO_EMISL_PBL ) THEN
          h2 = HcoState%Grid%PBLHEIGHT%Val(I,J)
       ELSE
          h2 = SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:UppLL))
       ENDIF

       ! Adjust dh if we are in lowest level
       IF ( L == LowLL ) THEN
          dh = SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:LowLL)) - h1
       ENDIF

       ! Adjust dh if we are in top level
       IF ( L == UppLL ) THEN
          dh = h2 - SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:(UppLL-1)))
       ENDIF

       ! compute dilution factor: the new flux should emit the same mass per
       ! volume, i.e. flux_total/column_total = flux_level/column_level
       ! --> flux_level = fluxtotal * column_level / column_total.
       IF ( h2 > h1 ) THEN
          DilFact = dh / ( h2 - h1 )
       ELSE
          !write(*,*) 'Warning: GetDilFact h2 not greater than h1!!! ',LowLL,UppLL,L,EmisL1,EmisL2,EmisL1Unit,EmisL2Unit
          DilFact = 1.0_hp
       ENDIF

    ! Approxiate dilution factor otherwise
    ELSE

       DilFact = 1.0_hp / REAL(UppLL-LowLL+1,hp)
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE GetDilFact
!EOC
END MODULE HCO_Calc_Mod
