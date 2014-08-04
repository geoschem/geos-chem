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
  USE HCO_DataCont_Mod, ONLY : DataCont, ListCont, SclMax
  USE HCO_DataCont_Mod, ONLY : Pnt2DataCont

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCO_CalcEmis
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: GET_CURRENT_EMISSIONS
!
! ============================================================================
!
! !REVISION HISTORY:
!  25 Aug 2012 - C. Keller   - Initial version.
!  06 Jun 2014 - R. Yantosca - Add cosmetic changes in ProTeX headers
!  08 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
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
! !IROUTINE: HCO_CalcEmis
!
! !DESCRIPTION: Subroutine HCO\_CalcEmis calculates the 3D emission 
! fields at current datetime for the specified species, categories, and 
! extension number. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_CalcEmis( am_I_Root, HcoState, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
    USE HCO_ARR_MOD,      ONLY : HCO_ArrAssert
    USE HCO_EMISLIST_MOD, ONLY : EmisList_NextCont
    USE HCO_FILEDATA_MOD, ONLY : FileData_ArrIsDefined
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root  ! Root CPU?
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Working pointers: list and data container 
    TYPE(ListCont), POINTER :: Lct => NULL()
    TYPE(DataCont), POINTER :: Dct => NULL()

    ! Temporary emission arrays
    REAL(hp), POINTER       :: Diag3D(:,:,:) => NULL()
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

    ! Integers
    INTEGER             :: ThisSpc, PrevSpc ! current and previous species ID
    INTEGER             :: ThisCat, PrevCat ! current and previous category 
    INTEGER             :: ThisHir, PrevHir ! current and previous hierarchy 
    INTEGER             :: SpcMin,  SpcMax  ! range of species to be considered 
    INTEGER             :: CatMin,  CatMax  ! range of categories to be considered 
    INTEGER             :: ExtNr            ! Extension Nr to be used 
    INTEGER             :: nI, nJ, nL 
    INTEGER             :: nnSpec, FLAG

    LOGICAL             :: DoDiagn

    ! For error handling & verbose mode
    LOGICAL             :: verb
    CHARACTER(LEN=255)  :: MSG

    ! testing / debugging
    integer :: ix,iy

    ! testing only
    integer :: modid

    !=================================================================
    ! HCO_CalcEmis begins here!
    !=================================================================

    ! testing only
    ix = 30 
    iy = 34

    ! Enter routine 
    CALL HCO_ENTER ('HCO_CalcEmis (HCO_CALC_MOD.F90)', RC )
    IF(RC /= HCO_SUCCESS) RETURN

    ! verb mode? 
    verb = HCO_VERBOSE_CHECK() .AND. am_I_Root

    !-----------------------------------------------------------------
    ! Initialize variables 
    !-----------------------------------------------------------------

    ! Initialize
    SpcFlx(:,:,:)    = 0.0_hp
    CatFlx(:,:,:)    = 0.0_hp
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
    IF ( verb ) THEN
       WRITE (MSG, *) 'Run HEMCO calculation w/ following options:'
       CALL HCO_MSG ( MSG )
       WRITE (MSG, "(A20,I5)")    'Extension number:', ExtNr 
       CALL HCO_MSG ( MSG )
       WRITE (MSG, "(A20,I5,I5)") 'Tracer range:', SpcMin, SpcMax
       CALL HCO_MSG ( MSG )
       WRITE (MSG, "(A20,I5,I5)") 'Category range:', CatMin, CatMax
       CALL HCO_MSG ( MSG )
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
    Lct => NULL()
    CALL EmisList_NextCont ( Lct, FLAG ) 

    ! Do until end of EmisList (==> loop over all emission containers) 
    DO WHILE ( FLAG == HCO_SUCCESS )

       ! ------------------------------------------------------------
       ! Select container and update all working variables & arrays.
       ! ------------------------------------------------------------

       ! Dct is the current data container 
       Dct => Lct%Dct

       ! Check if this is a base field (type = 1).
       IF ( Dct%DctType /= 1 ) THEN
          CALL EmisList_NextCont ( Lct, FLAG )
          CYCLE
       ENDIF

       ! Sanity check: Make sure this container holds data.
       ! 'Empty' containers are possible if the simulation time
       ! is outside of the specified data time range and time
       ! slice cycling is deactivated (CycleFlag > 1). 
       IF( .NOT. FileData_ArrIsDefined(Lct%Dct%Dta) ) THEN
          CALL EmisList_NextCont ( Lct, FLAG )
          CYCLE
       ENDIF

       ! Check if this is the specified extension number
       IF ( Dct%ExtNr /= ExtNr ) THEN 
          CALL EmisList_NextCont ( Lct, FLAG )
          CYCLE
       ENDIF

       ! Advance to next container if the species ID is outside 
       ! the specified species range (SpcMin - SpcMax). Consider 
       ! all species above SpcMin if SpcMax is negative!
       IF( (  Dct%HcoID < SpcMin                     ) .OR. &
           ( (Dct%HcoID > SpcMax) .AND. (SpcMax > 0) ) ) THEN
          CALL EmisList_NextCont ( Lct, FLAG )
          CYCLE
       ENDIF

       ! Advance to next emission field if the emission category of 
       ! the current container is outside of the specified species 
       ! range (CatMin - CatMax). Consider all categories above CatMin
       ! if CatMax is negative!
       IF( (  Dct%Cat < CatMin                     ) .OR. &
           ( (Dct%Cat > CatMax) .AND. (CatMax > 0) ) ) THEN
          CALL EmisList_NextCont ( Lct, FLAG )
          CYCLE
       ENDIF

       ! Update working variables
       ThisSpc = Dct%HcoID
       ThisCat = Dct%Cat
       ThisHir = Dct%Hier

       !--------------------------------------------------------------------
       ! If this is a new species or category, pass the previously collected 
       ! emissions to the species array. Update diagnostics at category level.
       ! Skip this step for first species, i.e. if PrevSpc is still -1. 
       !--------------------------------------------------------------------
       IF ( (PrevSpc>0) .AND. (ThisCat/=PrevCat .OR. ThisSpc/=PrevSpc) ) THEN

          ! CatFlx holds the emissions for this category. Pass this to 
          ! the species array SpcFlx.
          SpcFlx(:,:,:) = SpcFlx(:,:,:) + CatFlx(:,:,:)

          ! Add category emissions to diagnostics at category level
          ! (only if defined in the diagnostics list).
          IF ( Diagn_AutoFillLevelDefined(3) .AND. DoDiagn ) THEN 
             Diag3D => CatFlx
             CALL Diagn_Update( am_I_Root,   HcoState, ExtNr=ExtNr,   &
                                Cat=PrevCat, Hier=-1,  HcoID=PrevSpc, &
                                AutoFill=1,  Array3D=Diag3D, RC = RC ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
             Diag3D => NULL() 
          ENDIF

          ! Reset CatFlx array and the previously used hierarchy 
          ! ==> Emission hierarchies are only important within the 
          ! same category, hence always start over at lowest hierarchy
          ! when entering a new category.
          CatFlx(:,:,:)  = 0.0_df
          PrevHir        = -1
          ThisCat        = Dct%Cat

       ENDIF

       !--------------------------------------------------------------------
       ! If this is a new species, pass previously calculated emissions
       ! to the final emissions array in HcoState. 
       ! Update diagnostics at extension number level. 
       !--------------------------------------------------------------------
       IF ( ThisSpc /= PrevSpc ) THEN

          ! Don't do before first emission calculation, i.e. if PrevSpc 
          ! is still the initialized value of -1!
          IF ( PrevSpc > 0 ) THEN

             ! Add to OutArr
             OutArr(:,:,:) = OutArr(:,:,:) + SpcFlx(:,:,:)

             ! Add to diagnostics at extension number level. 
             ! The same diagnostics may be updated multiple times during 
             ! the same time step, continuously adding emissions to it.
             IF ( Diagn_AutoFillLevelDefined(2) .AND. DoDiagn ) THEN 
                Diag3D => SpcFlx
                CALL Diagn_Update(am_I_Root, HcoState, ExtNr=ExtNr,   &
                                  Cat=-1,    Hier=-1,  HcoID=PrevSpc, &
                                  AutoFill=1,Array3D=Diag3D, RC = RC ) 
                IF ( RC /= HCO_SUCCESS ) RETURN
                Diag3D => NULL()
             ENDIF
 
             ! Reset arrays and previous hierarchy. 
             SpcFlx(:,:,:)  =  0.0_df
             PrevCat        =  -1
             OutArr         => NULL()
          ENDIF

          ! Update number of species for which emissions have been
          ! calculated. 
          nnSpec = nnSpec + 1

          ! To write emissions into temporary array, make OutArr point
          ! to the buffer array HcoState%Buffer3D. 
          IF ( HcoState%Options%FillBuffer ) THEN

             ! Cannot use temporary array for more than one species!
             IF ( nnSpec > 1 ) THEN
                MSG = 'Cannot fill buffer for more than one species!'
                CALL HCO_ERROR( MSG, RC ) 
                RETURN
             ENDIF

             ! Point to array and check allocation status as well as 
             ! array size.
             OutArr => HcoState%Buffer3D%Val
             IF ( .NOT. ASSOCIATED( OutArr ) ) THEN
                MSG = 'Buffer array is not associated'
                CALL HCO_ERROR( MSG, RC )
                RETURN
             ENDIF
             IF ( (SIZE(OutArr,1) /= nI) .OR. &
                  (SIZE(OutArr,2) /= nJ) .OR. &
                  (SIZE(OutArr,3) /= nL)       ) THEN
                MSG = 'Buffer array has wrong dimension!'
                CALL HCO_ERROR( MSG, RC )
                RETURN
             ENDIF

          ! To write emissions directly into HcoState%Emsr3D, make OutArr
          ! point to current species' emission array in HcoState.
          ELSE

             ! Check allocation status of emission array in HcoState and
             ! allocate if necessary.
             CALL HCO_ArrAssert( HcoState%Spc(ThisSpc)%Emis, &
                                 nI, nJ, nL, RC             )
             IF ( RC /= HCO_SUCCESS ) RETURN
             OutArr => HcoState%Spc(ThisSpc)%Emis%Val
          ENDIF

          ! verbose mode
          IF ( verb ) THEN
             write(MSG,*) 'Now calculating emissions for species ', &
                           TRIM(HcoState%Spc(ThisSpc)%SpcName)
             CALL HCO_MSG( MSG, SEP1='-', SEP2='-' )
          ENDIF
 
       ENDIF

       !--------------------------------------------------------------------
       ! Define TmpFlx. This is the array containing the emissions 
       ! for the current container. Set initial values to -999 and 
       ! not to 0 in order to be able to distinguish between untouched 
       ! grid boxes (e.g. for regional emissions or invetories that
       ! don't extent through the entire troposphere) and boxes with 
       ! defined but zero emissions.
       !--------------------------------------------------------------------
       TmpFlx(:,:,:) = -999.0_df
       CALL GET_CURRENT_EMISSIONS( am_I_Root, HcoState, & 
                                   Dct,    nI, nJ, nL, TmpFlx, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! ------------------------------------------------------------
       ! Collect all emissions of the same category (and species) in 
       ! array CatFlx.
       ! The specified field hierarchies determine whether the
       ! temporary emissions are added to CatFlx (if hierarchy is
       ! the same as the previously used hierarchy), or if they 
       ! overwrite the previous values in CatFlx (if hierarchy is 
       ! higher than the previous hierarchy).
       ! ------------------------------------------------------------

       ! Add emissions to the category array CatFlx if this hierarchy
       ! is the same as previous hierarchy
       IF ( ThisHir == PrevHir ) THEN

          ! Ignore negative emissions!
          WHERE ( TmpFlx >= 0.0_df )
             CatFlx = CatFlx + TmpFlx
          END WHERE

!          ! testing only
!          IF ( verb ) THEN
!             write(lun,*) 'Field ', TRIM(Dct%cName),              &
!                        ' added to emissions (tracer ', ThisSpc,     &
!                        '; Category = ', ThisCat, ')' 
!          ENDIF
 
       ! If hierarchy is larger than those of the previously used
       ! fields, overwrite CatFlx w/ new values. 
       ELSEIF ( ThisHir > PrevHir ) THEN
        
          ! Ignore negative emissions!
          WHERE ( TmpFlx >= 0.0_df )
             CatFlx = TmpFlx
          END WHERE

!          ! testing only
!          IF ( verb ) THEN
!             write(lun,*) 'Field ', TRIM(Dct%cName),              &
!                        ' replaced old emissions (tracer ', ThisSpc, &
!                        '; Category = ', ThisCat, ')' 
!          ENDIF

       ELSE
          MSG = 'Hierarchy error in calc_emis: ' // TRIM(Dct%cName)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Update diagnostics at hierarchy level. Make sure that only 
       ! positive values are used.
       ! The same diagnostics may be updated multiple times 
       ! during the same time step, continuously adding
       ! emissions to it. 
       IF ( Diagn_AutoFillLevelDefined(4) .AND. DoDiagn ) THEN 
          Diag3D => TmpFlx
          CALL Diagn_Update( am_I_Root,  HcoState,       ExtNr=ExtNr,   &
                             Cat=ThisCat,Hier=ThisHir,   HcoID=ThisSpc, &
                             AutoFill=1, Array3D=Diag3D, PosOnly=.TRUE.,&
                             RC=RC ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
          Diag3D => NULL()
       ENDIF

       ! Update previously used species, category and hierarchy
       PrevSpc = ThisSpc
       PrevCat = ThisCat
       PrevHir = ThisHir
 
       ! Advance to next emission container
       CALL EmisList_NextCont( Lct, FLAG )

    ENDDO ! Loop over EmisList

    !=======================================================================
    ! Also pass emissions of the last category to output array 
    !=======================================================================
    IF ( nnSpec > 0 ) THEN
       SpcFlx(:,:,:) = SpcFlx(:,:,:) + CatFlx(:,:,:)
       OutArr(:,:,:) = OutArr(:,:,:) + SpcFlx(:,:,:)

       ! Diagnostics at category level
       IF ( Diagn_AutoFillLevelDefined(3) .AND. DoDiagn ) THEN 
          Diag3D => CatFlx
          CALL Diagn_Update( am_I_Root,   HcoState, ExtNr=ExtNr,   &
                             Cat=PrevCat, Hier=-1,  HcoID=PrevSpc, &
                             AutoFill=1,  Array3D=Diag3D, RC = RC ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
          Diag3D => NULL() 
       ENDIF

       ! Diagnostics at extension number level
       IF ( Diagn_AutoFillLevelDefined(2) .AND. DoDiagn ) THEN 
          Diag3D => SpcFlx
          CALL Diagn_Update( am_I_Root,  HcoState,       ExtNr=ExtNr,   &
                             Cat=-1,     Hier=-1,        HcoID=PrevSpc, &
                             AutoFill=1, Array3D=Diag3D, RC = RC         )
          IF ( RC /= HCO_SUCCESS ) RETURN
          Diag3D => NULL()
       ENDIF
    ENDIF ! nnSpec > 0

    ! Make sure internal pointers are nullified 
    Lct    => NULL()
    Dct    => NULL()
    OutArr => NULL()

    ! verbose
    IF ( verb ) THEN
       WRITE (MSG, *) 'HEMCO emissions successfully calculated!'
       CALL HCO_MSG ( MSG )
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE HCO_CalcEmis
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
  SUBROUTINE Get_Current_Emissions( am_I_Root, HcoState, BaseDct, &
                                    nI, nJ, nL, OUTARR_3D, RC )
!
! !USES:
!
    USE HCO_State_Mod,    ONLY : HCO_State
    USE HCO_tIdx_MOD,     ONLY : tIDx_GetIndxVec
    USE HCO_FileData_Mod, ONLY : FileData_ArrIsDefined
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN )   :: am_I_Root           ! Root CPU?
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
    INTEGER,         INTENT(INOUT) :: RC
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
!  25 Aug 2012 - C. Keller   -  Initial Version
!  09 Nov 2012 - C. Keller   -  MASK update. Masks are now treated
!                               separately so that multiple masks can be 
!                               added.
!  06 Jun 2014 - R. Yantosca -  Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DataCont), POINTER :: ScalDct => NULL()

    ! Scalars
    REAL(hp)                :: MASK(nI,nJ,1)
    REAL(hp)                :: TMPVAL
    INTEGER                 :: tIdxVec(nI), tIDx
    INTEGER                 :: IDX
    INTEGER                 :: I, J, L, N
    INTEGER                 :: BaseLL, ScalLL, TmpLL
    INTEGER                 :: IJFILLED
    LOGICAL                 :: DO_MASK, ERR
    CHARACTER(LEN=255)      :: MSG, LOC
 
    ! testing only
    INTEGER                 :: IX, IY
    LOGICAL                 :: verb

    !=================================================================
    ! GET_CURRENT_EMISSIONS begins here
    !=================================================================

    ! Enter
    CALL HCO_ENTER('GET_CURRENT_EMISSIONS', RC )
    IF(RC /= HCO_SUCCESS) RETURN
    ERR = .FALSE.

    ! Verbose mode 
    verb = HCO_VERBOSE_CHECK() .AND. am_I_Root 

    ! testing only:
    IX = 25 !-1 
    IY = 25 !-1 

    ! Check if container contains data
    IF ( .NOT. FileData_ArrIsDefined(BaseDct%Dta) ) THEN
       MSG = 'Array not defined: ' // TRIM(BaseDct%cName)
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF

    ! Initialize mask
    MASK(:,:,:) = 0.0_hp
    DO_MASK     = .FALSE.

    ! Verbose 
    IF ( verb ) THEN
       write(MSG,*) '--> GET EMISSIONS FOR ', TRIM(BaseDct%cName)
       CALL HCO_MSG(MSG)
    ENDIF

    ! ----------------------------------------------------------------
    ! Set base emissions
    ! ----------------------------------------------------------------

    ! Get vertical extension of base emission array.
    ! Unlike the output array OUTARR_3D, the data containers do not
    ! necessarily extent over the entire troposphere but only cover
    ! the effectively filled vertical levels. For most inventories, 
    ! this is only the first model level.
    IF ( BaseDct%Dta%SpaceDim==3 ) THEN 
       BaseLL = SIZE(BaseDct%Dta%V3(1)%Val,3) 
    ELSE
       BaseLL = 1
    ENDIF

    ! Precalculate timeslice index. The data containers can 
    ! carry 2D/3D arrays for multiple time steps (i.e. for 
    ! every hour of the day), stored in a vector.
    ! tIdxVec contains the vector index to be used at the current
    ! datetime. This parameter may vary with longitude due to time
    ! zone shifts! 
    tIdxVec = tIDx_GetIndxVec( BaseDct%Dta, nI ) 

    ! Loop over all latitudes and longitudes
!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT( SHARED )                                                &
!$OMP PRIVATE( I, J, L, tIdx, IJFILLED, TMPVAL                       ) & 
!$OMP SCHEDULE( DYNAMIC )
    DO J = 1, nJ
    DO I = 1, nI
 
       ! # of levels w/ defined emissions 
       IJFILLED = 0

       ! Time slice index for this lon
       tIdx = tIdxVec(I)

       ! Loop over all levels
       DO L = 1, BaseLL

          ! Get base value. Use uniform value if scalar field.
          IF ( BaseDct%Dta%SpaceDim == 1 ) THEN
             TMPVAL = BaseDct%Dta%V2(tIDx)%Val(1,1)
          ELSEIF ( BaseDct%Dta%SpaceDim == 2 ) THEN
             TMPVAL = BaseDct%Dta%V2(tIDx)%Val(I,J)
          ELSE
             TMPVAL = BaseDct%Dta%V3(tIDx)%Val(I,J,L)
          ENDIF
          
          ! Advance to next grid box if base value is negative, 
          ! indicating that this inventory is not defined over 
          ! this grid box.
          IF ( TMPVAL < 0.0_hp ) CYCLE 

          ! Pass base value to output array
          OUTARR_3D(I,J,L) = TMPVAL

          ! Update IJFILLED
          IJFILLED = IJFILLED + 1

       ENDDO !L

       ! If emissions are defined for at least one level, make 
       ! sure that emissions in all other levels are set to zero. 
       ! This is to make sure that higher hierarchy emissions entirely
       ! overwrite lower hierarchy emissions (emissions are only over-
       ! written where updated emissions are zero or higher).
       IF ( IJFILLED > 0 ) THEN
          WHERE ( OUTARR_3D(I,J,:) < 0.0_hp ) 
             OUTARR_3D(I,J,:) = 0.0_hp
          ENDWHERE
       ENDIF

    ENDDO !I
    ENDDO !J
!$OMP END PARALLEL DO

    ! ----------------------------------------------------------------
    ! Apply scale factors
    ! The container IDs of all scale factors associated with this base 
    ! container are stored in vector Scal_cID.
    ! ----------------------------------------------------------------

    ! Loop over maximum number of scale factors 
    DO N = 1, SclMax

       ! Get the scale factor container ID for the current slot
       IDX = BaseDct%Scal_cID(N)

       ! Leave if container ID is not defined. The container IDs
       ! were filled from left, i.e. beginning with the first element
       ! of Scal_cID. 
       IF ( IDX <= 0 ) THEN
          EXIT
       ENDIF

       ! Point to data container with the given container ID
       CALL Pnt2DataCont( IDX, ScalDct, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Sanity check: scale field cannot be a base field 
       IF ( (ScalDct%DctType == 1) ) THEN
          MSG = 'Wrong scale field type: ' // TRIM(ScalDct%cName)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

       ! Skip this scale factor if no data defined. This is possible
       ! if scale factors are only defined for a given time range and
       ! the simulation datetime is outside of this range.
       IF ( .NOT. FileData_ArrIsDefined(ScalDct%Dta) ) THEN
          MSG = 'Scale factor not defined: ' // TRIM(ScalDct%cName)
          CALL HCO_WARNING( MSG, RC )
          CYCLE
       ENDIF

       ! Get vertical extension of this scale factor array.
       IF( (ScalDct%Dta%SpaceDim<=2) ) THEN
          ScalLL = 1
       ELSE
          ScalLL = SIZE(ScalDct%Dta%V3(1)%Val,3)
       ENDIF

       ! Get vector of time slice indeces
       tIDxVec = tIDx_GetIndxVec( ScalDct%Dta, nI ) 

       ! Loop over all latitudes and longitudes
!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT( SHARED )                                                &
!$OMP PRIVATE( I, J, tIdx, TMPVAL, L, tmpLL                          ) & 
!$OMP SCHEDULE( DYNAMIC )
       DO J = 1, nJ
       DO I = 1, nI

          ! Get current time index
          tIdx = tIdxVec(I)
            
          ! ------------------------------------------------------------ 
          ! Check if this is a mask. If so, add mask values to the MASK
          ! array. For now, we assume masks to be binary, i.e. 0 or 1.
          ! We may want to change that in future to also support values
          ! in between. This is especially important when regridding 
          ! high resolution masks onto coarser grids! 
          ! ------------------------------------------------------------ 
          IF ( ScalDct%DctType == 3 ) THEN  

             ! Mask value over this grid box
             TMPVAL = ScalDct%Dta%V2(1)%Val(I,J)
 
             ! Mask values must not be negative 
             IF ( TMPVAL <= 0.0_hp ) THEN
                TMPVAL = 0.0_hp
             ELSE
                TMPVAL = 1.0_hp
             ENDIF

             ! For operator set to 3, mirror value
             IF ( ScalDct%Oper == 3 ) THEN
                TMPVAL = 1.0_hp - TMPVAL 
             ENDIF

             ! Add to mask and set mask flag to TRUE
             MASK(I,J,1) = MASK(I,J,1) + TMPVAL
             DO_MASK     = .TRUE.

             ! testing only
             IF ( verb .AND. I==1 .AND. J==1 ) THEN
                write(MSG,*) 'Mask field ', TRIM(ScalDct%cName),   &
                     ' found and added to temporary mask.'
                CALL HCO_MSG( MSG ) 
             ENDIF

             ! Advance to next grid box 
             CYCLE 
          ENDIF! DctType=3 

          ! ------------------------------------------------------------ 
          ! For non-mask fields, apply scale factors to all levels
          ! of the base field individually. If the scale factor
          ! field has more than one vertical level, use the
          ! vertical level closest to the corresponding vertical
          ! level of the base emission field
          ! ------------------------------------------------------------ 

          ! Loop over all vertical levels of the base field
          DO L = 1,BaseLL
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

             ! Advance to next grid box if scale factor is negative
             IF ( TMPVAL < 0.0_hp ) CYCLE

             ! -------------------------------------------------------
             ! Apply scale factor in accordance to field operator
             ! -------------------------------------------------------

             ! Oper 1: multiply
             IF ( ScalDct%Oper == 1 ) THEN
                OUTARR_3D(I,J,L) = OUTARR_3D(I,J,L) * TMPVAL

             ! Oper -1: divide 
             ELSEIF ( ScalDct%Oper == -1 ) THEN
                ! Ignore zeros to avoid NaN
                IF ( TMPVAL /= 0.0_hp ) THEN
                   OUTARR_3D(I,J,L) = OUTARR_3D(I,J,L) / TMPVAL
                ENDIF

             ! Oper 2: square
             ELSEIF ( ScalDct%Oper == 2 ) THEN
                OUTARR_3D(I,J,L) = OUTARR_3D(I,J,L) * TMPVAL * TMPVAL

             ! Return w/ error otherwise (Oper 3 is only allowed for masks!)
             ELSE
                MSG = 'Illegal data operator: ' // TRIM(ScalDct%cName)
                CALL HCO_ERROR( MSG, RC )
                ERR = .TRUE.
                EXIT
             ENDIF
          ENDDO !LL

          ! testing only
          if ( verb .and. i == ix .and. j == iy ) then
             write(MSG,*) 'Scale field ', TRIM(ScalDct%cName)
             CALL HCO_MSG( MSG )
             write(MSG,*) 'Time slice: ', tIdx
             CALL HCO_MSG( MSG )
             write(MSG,*) 'IX, IY: ', IX, IY
             CALL HCO_MSG( MSG )
             write(MSG,*) 'Scale factor (IX,IY,L1): ', TMPVAL
             CALL HCO_MSG( MSG )
             write(MSG,*) 'Mathematical operation : ', ScalDct%Oper 
             CALL HCO_MSG( MSG )
!             write(lun,*) 'Updt (IX,IY,L1): ', OUTARR_3D(IX,IY,1)
          endif

       ENDDO !I
       ENDDO !J
!$OMP END PARALLEL DO

       ! error check
       IF ( ERR ) THEN
          ScalDct => NULL()
          RC = HCO_FAIL
          RETURN
       ENDIF

    ENDDO ! N

    ! ----------------------------
    ! Masks 
    ! ----------------------------
    IF ( DO_MASK ) THEN

       ! Restrict mask values to a maximum of 1. 
       WHERE ( MASK > 1.0_hp ) MASK = 1.0_hp

       ! Apply mask. Make sure that emissions become negative
       ! outside the mask region. This is to make sure that these 
       ! grid boxes will be ignored when calculating the final  
       ! emissions. 
       DO L = 1, BaseLL
          WHERE ( MASK(:,:,1) <= 0.0_hp )
             OUTARR_3D(:,:,L) = -999.0_hp
          ELSEWHERE
             OUTARR_3D(:,:,L) = OUTARR_3D(:,:,L) * MASK(:,:,1) 
          ENDWHERE
       ENDDO

    ENDIF

    ! Cleanup and leave w/ success
    ScalDct => NULL()
    CALL HCO_LEAVE ( RC )
 
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
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Current_Emissions_B( am_I_Root, HcoState, BaseDct, &
                                      nI, nJ, nL, OUTARR_3D, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
    USE HCO_TIDX_MOD,     ONLY : tIDx_GetIndx
    USE HCO_FILEDATA_MOD, ONLY : FileData_ArrIsDefined
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)    :: am_I_Root           ! Root CPU?
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
    INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  25 Aug 2012 - C. Keller   - Initial Version
!  09 Nov 2012 - C. Keller   - MASK update. Masks are now treated
!                              separately so that multiple masks can be 
!                              added.
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DataCont), POINTER :: ScalDct => NULL()
    REAL(hp)                :: MASK
    REAL(hp)                :: TMPVAL
    INTEGER                 :: tIdx, IDX
    INTEGER                 :: I, J, L, N
    INTEGER                 :: BaseLL, ScalLL, TmpLL
    INTEGER                 :: IJFILLED
    LOGICAL                 :: DO_MASK, ERR
    CHARACTER(LEN=255)      :: MSG, LOC
 
    ! testing only
    INTEGER                 :: IX, IY
    LOGICAL                 :: verb

    !=================================================================
    ! GET_CURRENT_EMISSIONS_B begins here
    !=================================================================

    ! Enter
    CALL HCO_ENTER('GET_CURRENT_EMISSIONS_B', RC )
    IF(RC /= HCO_SUCCESS) RETURN

    ERR = .FALSE.

    ! testing only
    verb = HCO_VERBOSE_CHECK() .AND. am_I_Root
    IX = 60 !40 !19 43 61
    IY = 32 !36 !33 26 37

    ! Check if field data is defined
    IF ( .NOT. FileData_ArrIsDefined(BaseDct%Dta) ) THEN
       MSG = 'Array not defined: ' // TRIM(BaseDct%cName)
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF

    ! Testing only:
    IF ( verb ) THEN
       write(MSG,*) '--> GET EMISSIONS FOR ', TRIM(BaseDct%cName)
       CALL HCO_MSG(MSG)
    ENDIF

    ! Loop over all grid boxes
!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT( SHARED )                                                &
!$OMP PRIVATE( I, J, MASK, DO_MASK, BaseLL, tIdx, IJFILLED, L        ) & 
!$OMP PRIVATE( TMPVAL, N, IDX, ScalDct, ScalLL, tmpLL                ) & 
!$OMP SCHEDULE( DYNAMIC )
    DO J = 1, nJ
    DO I = 1, nI

       ! Initialize mask
       MASK    = 0.0_hp
       DO_MASK = .FALSE.

       ! -------------------------------------------------------------
       ! Set base emissions
       ! -------------------------------------------------------------

       ! Get vertical extension of base emission array.
       ! Unlike the output array OUTARR_3D, the data containers do not
       ! necessarily extent over the entire troposphere but only cover
       ! the effectively filled vertical levels. For most inventories, 
       ! this is only the first model level.
       IF ( BaseDct%Dta%SpaceDim==3 ) THEN 
          BaseLL = SIZE(BaseDct%Dta%V3(1)%Val,3) 
       ELSE
          BaseLL = 1
       ENDIF

       ! Precalculate timeslice index. The data containers can 
       ! carry 2D/3D arrays for multiple time steps (i.e. for 
       ! every hour of the day), stored in a vector.
       ! tIdxVec contains the vector index to be used at the current
       ! datetime. This parameter may vary with longitude due to time
       ! zone shifts! 
       tIdx = tIDx_GetIndx ( BaseDct%Dta, I )

       ! # of levels w/ defined emissions 
       IJFILLED = 0

       ! Loop over all levels
       DO L = 1, BaseLL

          ! Get base value. Use uniform value if scalar field.
          IF ( BaseDct%Dta%SpaceDim == 1 ) THEN
             TMPVAL = BaseDct%Dta%V2(tIDx)%Val(1,1)
          ELSEIF ( BaseDct%Dta%SpaceDim == 2 ) THEN
             TMPVAL = BaseDct%Dta%V2(tIDx)%Val(I,J)
          ELSE
             TMPVAL = BaseDct%Dta%V3(tIDx)%Val(I,J,L)
          ENDIF

          ! Advance to next grid box if base value is negative, 
          ! indicating that this inventory is not defined over 
          ! this grid box.
          IF ( TMPVAL < 0.0_hp ) CYCLE 

          ! Pass base value to output array
          OUTARR_3D(I,J,L) = TMPVAL

          ! Update IJFILLED
          IJFILLED = IJFILLED + 1

       ENDDO !L

       ! If emissions are defined for at least one level, make 
       ! sure that emissions in all other levels are set to zero.
       ! This is to make sure that higher hierarchy emissions entirely
       ! overwrite lower hierarchy emissions (emissions are only over-
       ! written where updated emissions are zero or higher).
       IF ( IJFILLED > 0 ) THEN
          WHERE ( OUTARR_3D(I,J,:) < 0.0_hp ) 
             OUTARR_3D(I,J,:) = 0.0_hp
          ENDWHERE
       ENDIF

       ! -------------------------------------------------------------
       ! Apply scale factors
       ! The container IDs of all scale factors associated with this base 
       ! container are stored in vector Scal_cID.
       ! -------------------------------------------------------------
   
       ! Loop over maximum number of scale factors 
       DO N = 1, SclMax
   
          ! Get the scale factor container ID for the current slot
          IDX = BaseDct%Scal_cID(N)
   
          ! Leave if field ID is not defined
          IF ( IDX <= 0 ) THEN
             EXIT
          ENDIF
   
          ! Point to emission container with the given container ID
          CALL Pnt2DataCont( IDX, ScalDct, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ERR = .TRUE.; EXIT
          ENDIF
   
          ! Scale field cannot be a base field 
          IF ( (ScalDct%DctType == 1) ) THEN
             MSG = 'Wrong scale field type: ' // TRIM(ScalDct%cName)
             CALL HCO_ERROR( MSG, RC )
             ERR = .TRUE.
             EXIT 
          ENDIF
 
          ! Skip this scale factor if no data defined. This is possible
          ! if scale factors are only defined for a given time range and
          ! the simulation datetime is outside of this range.
          IF ( .NOT. FileData_ArrIsDefined(ScalDct%Dta) ) THEN
             MSG = 'Array not defined: ' // TRIM(ScalDct%cName)
             CALL HCO_WARNING( MSG, RC )
             CYCLE
          ENDIF

          ! Get vertical extension of this scale factor array.
          IF( (ScalDct%Dta%SpaceDim<=2) ) THEN
             ScalLL = 1
          ELSE
             ScalLL = SIZE(ScalDct%Dta%V3(1)%Val,3)
          ENDIF
   
          ! Get current time index
          tIdx = tIDx_GetIndx ( ScalDct%Dta, I )
            
          ! ------------------------------------------------------------ 
          ! Check if this is a mask. If so, add mask values to the MASK
          ! array. For now, we assume masks to be binary, i.e. 0 or 1.
          ! We may want to change that in future to also support values
          ! in between. This is especially important when regridding
          ! high resolution masks onto coarser grids! 
          ! ------------------------------------------------------------ 
          IF ( ScalDct%DctType == 3 ) THEN  

             TMPVAL = ScalDct%Dta%V2(1)%Val(I,J)

             ! Mask values are restricted to 0 or 1
             IF ( TMPVAL <= 0.0_hp ) THEN
                TMPVAL = 0.0_hp
             ELSE
                TMPVAL = 1.0_hp
             ENDIF

             ! For operator set to 3, mirror value
             IF ( ScalDct%Oper == 3 ) THEN
                TMPVAL = 1.0_hp - TMPVAL 
             ENDIF

             ! Add to mask and set mask flag to TRUE
             MASK        = MASK + TMPVAL
             DO_MASK     = .TRUE.

             ! testing only
             if ( verb .and. i == ix .and. j == iy ) then
                write(*,*) 'Mask field ', TRIM(ScalDct%cName),   &
                     ' found and added to temporary mask.' 
             ENDIF

             ! Advance to next scale factor 
             CYCLE 
          ENDIF! DctType=3 

          ! ------------------------------------------------------------ 
          ! For non-mask fields, apply scale factors to all levels
          ! of the base field individually. If the scale factor
          ! field has more than one vertical level, use the
          ! vertical level closest to the corresponding vertical
          ! level in the base emission field
          ! ------------------------------------------------------------ 

          ! Loop over all vertical levels of the base field
          DO L = 1,BaseLL
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

             ! Advance to next grid box if scale factor is negative
             IF ( TMPVAL < 0.0_hp ) CYCLE

             ! -------------------------------------------------------
             ! Apply scale factor according to field operator
             ! -------------------------------------------------------

             ! Oper 1: multiply
             IF ( ScalDct%Oper == 1 ) THEN
                OUTARR_3D(I,J,L) = OUTARR_3D(I,J,L) * TMPVAL

             ! Oper -1: divide 
             ELSEIF ( ScalDct%Oper == -1 ) THEN
                ! Ignore zeros to avoid NaN
                IF ( TMPVAL /= 0.0_hp ) THEN
                   OUTARR_3D(I,J,L) = OUTARR_3D(I,J,L) / TMPVAL
                ENDIF

             ! Oper 2: square
             ELSEIF ( ScalDct%Oper == 2 ) THEN
                OUTARR_3D(I,J,L) = OUTARR_3D(I,J,L) * TMPVAL * TMPVAL

             ! Return w/ error otherwise (Oper 3 only allowed for masks!)
             ELSE
                MSG = 'Illegal data operator: ' // TRIM(ScalDct%cName)
                CALL HCO_ERROR( MSG, RC )
                ERR = .TRUE.
                EXIT 
             ENDIF
          ENDDO !LL
       ENDDO ! N

       ! ----------------------------
       ! Masks 
       ! ----------------------------
       IF ( DO_MASK ) THEN
   
          ! Restrict mask values to a maximum of 1. 
          IF ( MASK > 1.0_hp ) MASK = 1.0_hp
   
          ! Apply the mask. Make sure that emissions become negative
          ! outside the mask region. This is to make sure that these 
          ! grid boxes will be ignored when calculating the final  
          ! emissions. 
          IF ( MASK <= 0.0_hp ) THEN
             OUTARR_3D(I,J,:) = -999.0_hp
          ELSE
             OUTARR_3D(I,J,:) = OUTARR_3D(I,J,:) * MASK
          ENDIF
       ENDIF
       
    ENDDO !I
    ENDDO !J
!$OMP END PARALLEL DO

    ! Cleanup 
    ScalDct => NULL()

    ! Check exit status
    IF ( ERR ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Leave
    CALL HCO_LEAVE( RC )

  END SUBROUTINE Get_Current_Emissions_B
!EOC
END MODULE HCO_Calc_Mod
