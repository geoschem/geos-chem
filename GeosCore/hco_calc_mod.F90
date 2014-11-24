!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: HCO_CALC_MOD.F90
!
! !DESCRIPTION: Module HCO\_CALC\_MOD contains the HEMCO calculation
!  routines.\\
!  The emissions are calculated using the EMISSIONS linked list (see
!  NG\_EMIS\_MOD.F90). In detail, the current emissions are obtained by 
!  multiplying the base emissions with the associated scaling factors.
!  Both the base emissions and the scaling factors can vary with time
!  (depending on their 4th dimension and the definition of the temporal 
!  resolution, see also NG\_EMIS\_TIME\_MOD.F90).\\
!  Once the current emissions are calculated, they are stored in the
!  emission array State_Chm%Ems_Nomix, which contains the emissions for 
!  each tracer separately.\\
!  The field parameters CATEGORY and HIERARCHY defines the order in
!  which emissions are handled. Within the same category, higher hierarchy 
!  emissions overwrite lower hierarchy emission fields. Fields with 
!  the same hierarchy are added up, as are fields with hierarchy 0.
!  Emissions are fully calculated for a given category and only then 
!  passed to the emission array.
!  Note that negative emissions are not supported and are ignored.
!  See NG\_EMIS\_MOD.F90 for more details on the emissions structure and
!  calculations.
!
!  All emissions are in kg/m2/s. 
! ============================================================================
! \\
! !INTERFACE: 
!
      MODULE HCO_CALC_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD
      USE HCO_EMISLL_MOD,      ONLY : EmisCont
      USE HCO_EMISLL_MOD,      ONLY : SclMax
      USE HCO_EMISLL_MOD,      ONLY : Pnt2EmisCont

      IMPLICIT NONE

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CALC_EMIS
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: GET_CURRENT_EMISSIONS
!
! !REVISION HISTORY:
!  25 Aug 2012 - C. Keller - Initialization (stripped from ng_emis_mod.F90)
!
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CALC_EMIS
!
! !DESCRIPTION: Subroutine CALC_EMIS calculates the 3D emission fields 
!  for all tracers and for the current (emission) time step. Vertical
!  mixing into the boundary layer is not yet accounted for, and the
!  emission fields are saved into HcoState%Emsr3D(TrcID)%Arr3D.
!  All emissions in [kg/m2/s]. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_EMIS ( am_I_Root, HcoState, EmisLL, Clock, RC )
!
! !USES:
!
      USE HCO_TYPE_MOD,     ONLY : HCO_State, AllocCheck
      USE HCO_EMISLL_MOD,   ONLY : Get_nnEmisCont
      USE HCO_EMISLL_MOD,   ONLY : Get_nnIDList
      USE HCO_EMISLL_MOD,   ONLY : Create_IDList
      USE HCO_EMISLL_MOD,   ONLY : Show_EmisLL
      USE HCO_TIME_MOD,     ONLY : Update_tSlc
      USE HCO_TIME_MOD,     ONLY : HcoClock
!
! !ARGUMENTS:
!
      LOGICAL,         INTENT(IN   )  :: am_I_Root
      TYPE(HCO_State), POINTER        :: HcoState   ! HEMCO state
      TYPE(EmisCont),  POINTER        :: EmisLL     ! Emissions linked list
      TYPE(HcoClock),  POINTER        :: Clock
      INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!     13 March 2012 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Current field to work on
      TYPE(EmisCont), POINTER :: ThisCnt => NULL()

      ! Accumulated emissions per category (on emission grid)
      REAL*8, TARGET      :: CatFlx (HcoState%ISIZE, &
                                     HcoState%JSIZE, &
                                     HcoState%LSIZE )
      REAL*8, TARGET      :: TmpFlx (HcoState%ISIZE, &
                                     HcoState%JSIZE, &
                                     HcoState%LSIZE )

      ! Integers
      INTEGER             :: ThisTrc, LastTrc ! this and last used tracer ID
      INTEGER             :: ThisCat, LastCat ! this and last used category 
      INTEGER             :: ThisHir, LastHir ! this and last used hierarchy 
      INTEGER             :: TrcMin,  TrcMax  ! tracers to be considered 
      INTEGER             :: CatMin,  CatMax  ! categories to be considered 
      INTEGER             :: ExtNr 
      INTEGER             :: nI, nJ, nL 
      INTEGER             :: nnSpec

      ! Working array
      REAL*8, POINTER     :: TrcArr(:,:,:) => NULL()

      ! For error handling
      CHARACTER(LEN=255)  :: LOC, MSG

      ! testing / debugging
      logical :: debug

      !=================================================================
      ! CALC_EMIS begins here!
      !=================================================================

      ! Enter routine 
      LOC = 'CALC_EMIS (HCO_CALC_MOD.F90)'

      ! testing
      debug = HcoState%verbose .AND. am_I_Root

      ! Initialize
      CatFlx(:,:,:)    = 0d0
      LastTrc          = -1
      LastHir          = -1
      LastCat          = -1
      nnSpec           = 0

      ! Pass emission grid dimensions
      nI = HcoState%ISIZE
      nJ = HcoState%JSIZE
      nL = HcoState%LSIZE

      ! Pass calculation options
      TrcMin = HcoState%TrcMin
      TrcMax = HcoState%TrcMax
      CatMin = HcoState%CatMin
      CatMax = HcoState%CatMax
      ExtNr  = HcoState%ExtNr

      ! for debugging
      IF ( debug ) THEN
         WRITE ( 6, * ) 'Run HEMCO calculation w/ following options:'
         WRITE ( 6, * ) 'Extension number: ', ExtNr 
         WRITE ( 6, * ) 'Tracer range    : ', TrcMin, TrcMax
         WRITE ( 6, * ) 'Category range  : ', CatMin, CatMax
      ENDIF

      ! Return if no emission fields are defined
      IF ( Get_nnEmisCont() == 0 ) THEN
         IF ( debug ) THEN
            WRITE ( 6, * ) 'No emission fields defined - Leave!'
         ENDIF
         RC = HCO_SUCCESS
         RETURN
      ENDIF

      ! Create the IDLIST which allows a quick access to all containers 
      ! of EmisLL based upon their field IDs. Do this whenever the total
      ! number of containers has changed. Typically, this is only the
      ! case after the first emission time step.
      IF ( Get_nnEmisCont() /= Get_nnIDList() ) THEN
         CALL Create_IDList ( am_I_Root, EmisLL, HcoState, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! Print emissions linked list
         IF ( debug ) CALL Show_EmisLL( am_I_Root, EmisLL )
      ENDIF

      ! Update pointers to time slices 
      CALL Update_tSlc ( Clock, HcoState, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Point to the head of the emissions linked list
      ThisCnt => EmisLL

      ! Loop over all emission containers 
      DO WHILE ( ASSOCIATED ( ThisCnt ) )

         ! Check if this is a base field (type = 1).
         IF ( ThisCnt%DataType /= 1 ) THEN
            ThisCnt => ThisCnt%NEXTCONT
            CYCLE
         ENDIF

         ! Check if this is the correct extension number
         IF ( ThisCnt%ExtNr /= ExtNr ) THEN 
            ThisCnt => ThisCnt%NEXTCONT
            CYCLE
         ENDIF

         ! Advance to next emission field if no (valid) tracer ID is defined
         ! for the current emission field (which is the case for all
         ! scale factors), or if the tracer ID is outside of the desired
         ! range. Ignore if maximum tracer ID is negative!
         IF ( TrcMax > 0 ) THEN
            IF ( ThisCnt%TrcID < TrcMin .OR. ThisCnt%TrcID > TrcMax ) THEN
               ThisCnt => ThisCnt%NEXTCONT
               CYCLE
            ENDIF
         ENDIF

         ! Skip if category of this container is outside of the
         ! desired range. Ignore if maximum category is negative!
         IF ( CatMax > 0 ) THEN
            IF ( ThisCnt%Cat < CatMin .OR. &
                 ThisCnt%Cat > CatMax       ) THEN
               ThisCnt => ThisCnt%NEXTCONT
               CYCLE
            ENDIF
         ENDIF

         ! If a valid tracer ID is found, update current tracer index,
         ! category and hierarchy. 
         ThisTrc = ThisCnt%TrcID
         ThisCat = ThisCnt%Cat
         ThisHir = ThisCnt%Hier

         ! If this is a new tracer or new category for the same tracer,
         ! pass the collected emissions of the former tracer/category 
         ! to the tracer array 
         IF ( (ThisTrc /= LastTrc .AND. LastTrc > 0      ) .OR.        &
              (ThisTrc == LastTrc .AND. ThisCat /= LastCat) ) THEN

            ! Pass CatFlx to the array.
            TrcArr(:,:,:) = TrcArr(:,:,:) + CatFlx(:,:,:)

            ! Reset CatFlux array and LastHir.
            CatFlx(:,:,:)  = 0d0
            LastHir        = -1
         ENDIF

         ! If we deal with a new tracer, then make sure that the working
         ! array TrcArr points to the corresponding array in HcoState.
         IF ( ThisTrc /= LastTrc ) THEN

            ! Update number of species for which emissions are
            ! calculated during this call
            nnSpec = nnSpec + 1

            ! To write emissions into temporary array:
            IF ( HcoState%FillTemp3D ) THEN

               ! Cannot use temporary array for more than one species!
               IF ( nnSpec > 1 ) THEN
                  MSG = 'Cannot fill Temp3D for more than one species!'
                  CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
               ENDIF

               ! Point to array and check allocation status as well as
               ! array size 
               TrcArr => HcoState%Temp3D
               IF ( .NOT. ASSOCIATED( TrcArr ) ) THEN
                  MSG = 'Temp3D array not associated'
                  CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
               ENDIF
               IF ( (SIZE(TrcArr,1) /= nI) .OR. &
                    (SIZE(TrcArr,2) /= nJ) .OR. &
                    (SIZE(TrcArr,3) /= nL)       ) THEN
                  MSG = 'TmpArr has wrong size!'
                  CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
               ENDIF

            ! To write emissions directly into HcoState%Emsr3D:
            ELSE

               ! Check allocation status of emission array in HcoState and
               ! allocate if necessary.
               CALL AllocCheck( HcoState, ThisTrc, 1, RC )
               IF ( RC /= HCO_SUCCESS ) RETURN

               ! Make working pointer TrcArr point to corresponding flux
               ! array in HEMCO state object. 
               TrcArr => HcoState%Emsr3d(ThisTrc)%Arr3d
            ENDIF

            ! Reset flux before filling 
            TrcArr(:,:,:) = 0d0
         ENDIF

         ! Define TmpFlx. This is the array containing the current
         ! emissions (lon,lat,lev). Set initial values to -999 and not 
         ! to 0 in order to be able to distinguish between untouched 
         ! grid boxes and boxes with defined but zero emissions.
         TmpFlx(:,:,:) = -999d0

         ! Multiply base emissions with the associated emission factors
         ! Grid boxes with no defined emissions remain -999. 
         CALL GET_CURRENT_EMISSIONS( am_I_Root, HcoState, & 
                                     ThisCnt, nI, nJ, nL, TmpFlx, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! ------------------------------------------------------------
         ! Collect all emissions of the same category and tracer in 
         ! the CatFlx array.
         ! The specified field hierarchies determine whether the
         ! temporary emissions are added to CatFlx (if hierarchy is
         ! the same as the last used hierarchy), or if they overwrite
         ! the previous values in CatFlx (if hierarchy is higher than 
         ! the previous hierarchy).
         ! ------------------------------------------------------------

         ! Add emissions to the category array CatFlx if this hierarchy
         ! is the same as last hierarchy
         IF ( ThisHir == LastHir ) THEN

            ! Ignore negative emissions!
            WHERE ( TmpFlx >= 0d0 )
               CatFlx = CatFlx + TmpFlx
            END WHERE
        
            ! testing only
            IF ( debug ) THEN
               write(*,*) 'Field ', TRIM(ThisCnt%cName),               &
                          ' added to emissions (tracer ', ThisTrc,     &
                          '; Category = ', ThisCat, ')' 
            ENDIF
 
         ! If hierarchy is larger than those of the previously used
         ! fields, overwrite CatFlx w/ new values. 
         ELSEIF ( ThisHir > LastHir ) THEN
        
            ! Ignore negative emissions!
            WHERE ( TmpFlx >= 0d0 )
               CatFlx = TmpFlx
            END WHERE

            ! testing only
            IF ( debug ) THEN
               write(*,*) 'Field ', TRIM(ThisCnt%cName),               &
                          ' replaced old emissions (tracer ', ThisTrc, &
                          '; Category = ', ThisCat, ')' 
            ENDIF

         ELSE
            MSG = 'Hierarchy error in calc_emis: ' // TRIM(ThisCnt%cName)
            CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
         ENDIF

         ! Update last used tracer, category and hierarchy
         LastTrc = ThisTrc
         LastCat = ThisCat
         LastHir = ThisHir

         ! Advance to next emission container
         ThisCnt => ThisCnt%NEXTCONT

      ENDDO

      ! Also pass the emissions of the last category to the final 
      ! emission array
      TrcArr(:,:,:) = TrcArr(:,:,:) + CatFlx(:,:,:)

      ! Free pointer
      ThisCnt => NULL()
      TrcArr  => NULL()

      ! debugging
      IF ( debug ) THEN
         write(*,*) 'leaving CALC_EMIS now!'
      ENDIF

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE CALC_EMIS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_CURRENT_EMISSIONS
!
! !DESCRIPTION: Subroutine GET_CURRENT_EMISSIONS calculates the current 
!  emissions for the specified emission field and passes the result to 
!  OUTARR_3D.\\
!  This subroutine is only called (by calc_emis) for fields with a valid
!  tracer ID, i.e. for base emission fields. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_CURRENT_EMISSIONS( am_I_Root, HcoState, BaseCnt, &
                                        nI, nJ, nL, OUTARR_3D, RC )
!
! !USES:
!
      USE HCO_TYPE_MOD,    ONLY : HCO_State
!
! !ARGUMENTS:
!
      LOGICAL,         INTENT(IN )   :: am_I_Root ! Root CPU?
      TYPE(HCO_State), POINTER       :: HcoState ! HEMCO state
      TYPE(EmisCont),  POINTER       :: BaseCnt   ! base emission field
      INTEGER,         INTENT(IN)    :: nI
      INTEGER,         INTENT(IN)    :: nJ
      INTEGER,         INTENT(IN)    :: nL
      REAL*8,          INTENT(INOUT) :: OUTARR_3D(:,:,:)
      INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!     25 May 2012 - C. Keller - Initial Version
!     09 Nov 2012 - C. Keller - MASK update. Masks are now treated
!                               separately so that multiple masks can be 
!                               added.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE(EMISCONT), POINTER     :: ScalCnt     => NULL()
      REAL*8                      :: MASK(nI,nJ,1)
      REAL*8                      :: TMPVAL
      INTEGER                     :: tSlc, IDX
      INTEGER                     :: I, J, L, N
      INTEGER                     :: BaseLL, ScalLL, TmpLL
      INTEGER                     :: IJFILLED
      LOGICAL                     :: DO_MASK
      CHARACTER(LEN=255)          :: MSG, LOC
 
      ! testing only
      INTEGER                     :: IX, IY
      LOGICAL                     :: debug

      !=================================================================
      ! GET_CURRENT_EMISSIONS begins here
      !=================================================================

      ! Enter
      LOC = 'GET_CURRENT_EMISSIONS'

      ! testing only
      debug = HcoState%Verbose

      ! testing only:
      IX = 60 !40 !19 43 61
      IY = 32 !36 !33 26 37

      ! Check if emission field is indeed defined
      IF ( .NOT. ASSOCIATED ( BaseCnt ) ) THEN
         CALL HCO_ERROR ( 'Base field not defined', LOC, RC ); RETURN
      ENDIF
      
      ! Check if field data is defined
      IF ( .NOT. ASSOCIATED ( BaseCnt%Array ) ) THEN
         CALL HCO_ERROR ( 'Base field data not defined', LOC, RC )
         RETURN
      ENDIF

      ! Initialize mask
      MASK(:,:,:) = 0d0
      DO_MASK     = .FALSE.

      ! Testing only:
      IF ( debug .and. am_I_Root ) THEN
         write(*,*) '##################################################'
         write(*,*) ' GET EMISSIONS FOR ', TRIM(BaseCnt%cName)
         write(*,*) '##################################################'
      ENDIF

      ! ----------------------------
      ! Base emissions
      ! ----------------------------

      ! Base field cannot be a mask field!
      IF ( BaseCnt%DataType /= 1 ) THEN
         MSG = 'Wrong base field type: ' // TRIM(BaseCnt%cName)
         CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
      ENDIF

      ! Get vertical extension of base emission array
      BaseLL = SIZE(BaseCnt%Array,3) 

      ! Loop over all latitudes and longitudes
      DO J = 1, nJ
      DO I = 1, nI
 
         ! Get current index of the 4th dimension. This parameter
         ! may vary with longitude due to the time shifts. 
         tSlc = BaseCnt%TimeSlc(I)

         ! Logical if any i,j box is filled
         IJFILLED = 0

         DO L = 1, BaseLL

            ! Get base value
            IF ( BaseCnt%IsScalar ) THEN
               TMPVAL = BaseCnt%Array(1,1,1,tSlc)
            ELSE
               TMPVAL = BaseCnt%Array(I,J,L,tSlc)
            ENDIF

            ! Advance to next grid box if base value is negative
            IF ( TMPVAL < 0d0 ) CYCLE 

            ! Set basevalue in output array
            OUTARR_3D(I,J,L) = TMPVAL

            ! Update IJFILLED
            IJFILLED = IJFILLED + 1

         ENDDO !L

         ! If any grid box (I,J,*) is filled, make sure that all
         ! emissions in this column are at least zero!
         IF ( IJFILLED > 0 ) THEN
            WHERE ( OUTARR_3D(I,J,:) < 0d0 ) 
               OUTARR_3D(I,J,:) = 0d0
            ENDWHERE
         ENDIF

         ! Testing only:
         if ( debug .and. i == ix .and. j == iy .and. am_I_Root ) THEN
         write(*,*) 'Base field ', TRIM(BaseCnt%cName)
         write(*,*) 'Time index: ', tSlc
         write(*,*) 'IX, IY: ', IX, IY
         write(*,*) 'Value (IX,IY,L1): ', BaseCnt%Array(IX,IY,1,tSlc)
         write(*,*) 'Updated emissions (IX,IY): ', OUTARR_3D(IX,IY,1)
         if ( BaseLL > 1 ) then
         write(*,*) 'Value (IX,IY,L2): ', BaseCnt%Array(IX,IY,2,tSlc)
         write(*,*) 'Updated emissions (IX,IY): ', OUTARR_3D(IX,IY,2)
         endif
         endif

      ENDDO !I
      ENDDO !J

      ! Testing only:
      IF ( debug .and. am_I_Root ) THEN
         write(*,*) 'Total level 1: ', SUM(OUTARR_3D(:,:,1))
         if ( BaseLL > 1 ) then
            write(*,*) 'Total level 2: ', SUM(OUTARR_3D(:,:,2))
         endif
      ENDIF

      ! ----------------------------
      ! Scale factors
      ! ----------------------------

      ! Loop over all scale factor slots
      DO N = 1, SclMax

         ! Get the scale factor container ID for the current slot
         IDX = BaseCnt%Scal_cID(N)

         ! Leave if field ID is not defined
         IF ( IDX <= 0 ) THEN
            EXIT
         ENDIF

         ! Point to emission container with the given field ID
         CALL Pnt2EmisCont ( IDX, ScalCnt, RC )

         ! Base field cannot be a mask field!
         IF ( (ScalCnt%DataType == 1) .OR. (ScalCnt%DataType == 4) ) THEN
            MSG = 'Wrong scale field type: ' // TRIM(ScalCnt%cName)
            CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
         ENDIF

         ! Check if fetched data array is valid
         IF ( .NOT. ASSOCIATED ( ScalCnt%Array ) ) THEN
            MSG = 'Field data must be defined!'
            CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
         ENDIF

         ! Get vertical extension of this scale factor array.  
         ScalLL = SIZE(ScalCnt%Array,3)

         ! Loop over all latitudes and longitudes
         DO J = 1, nJ
         DO I = 1, nI

            ! Get current time index
            tSlc = ScalCnt%TimeSlc(I) 
            
            ! ------------------------------------------------------------ 
            ! Check if this is a mask. If so, add mask values to the MASK
            ! array. Only consider valid (positive) mask values.
            ! Don't just set the value to 1 but rather add the values
            ! together. This is useful if masks are not binary but
            ! contain values between 0 and 1 as well. Mask values 
            ! will be restricted to a total of 1 lateron.
            ! ------------------------------------------------------------ 
            IF ( ScalCnt%DataType == 3 ) THEN  

               TMPVAL = ScalCnt%Array(I,J,1,1)

               ! Add to mask and set mask flag to TRUE
               MASK(I,J,1) = MASK(I,J,1) + TMPVAL
               DO_MASK     = .TRUE.

               ! testing only
               IF ( debug .and. am_I_Root .AND. I==1 .AND. J==1 ) THEN
                  write(*,*) 'Mask field ', TRIM(ScalCnt%cName),   &
                          ' found and added to temporary mask.' 
               ENDIF

               ! Advance to next grid box 
               CYCLE 
            ENDIF! DataType=3 

            ! ------------------------------------------------------------ 
            ! For no-mask fields, apply scale factors to all levels
            ! of the base field individually. If the scale factor
            ! field has more than one vertical level, use the
            ! vertical level closest to the corresponding vertical
            ! level in the base emission field
            ! ------------------------------------------------------------ 

            ! Loop over all vertical levels of the base field
            DO L = 1,BaseLL
               IF ( L > ScalLL ) THEN 
                  ! If the vertical level exceeds the number of available 
                  ! scale factor levels, use the highest available level.
                  TmpLL = ScalLL
            ELSE 
                  ! Otherwise use the same vertical level index.
                  TmpLL = L
               ENDIF

               ! Get scale factor for this grid box
               IF ( ScalCnt%IsScalar ) THEN
                  TMPVAL = ScalCnt%Array(1,1,1,tSlc)
               ELSE
                  TMPVAL = ScalCnt%Array(I,J,TmpLL,tSlc)
               ENDIF

               ! Advance to next grid box if scale factor is negative
               IF ( TMPVAL < 0d0 ) CYCLE

               ! Apply scale factor
               OUTARR_3D(I,J,L) = OUTARR_3D(I,J,L) * TMPVAL
            ENDDO !LL

            ! testing only
            if ( debug .and. i == ix .and. j == iy ) then
               write(*,*) 'Scale field ', TRIM(ScalCnt%cName)
               write(*,*) 'Time slice: ', tSlc
               write(*,*) 'IX, IY: ', IX, IY
               write(*,*) 'Scale factor (IX,IY,L1): ', TMPVAL
               write(*,*) 'Updt (IX,IY,L1): ', OUTARR_3D(IX,IY,1)
            endif

         ENDDO !I
         ENDDO !J

         ! testing only
         IF ( debug .and. am_I_Root ) THEN
            write(*,*) 'Total level 1: ', SUM(OUTARR_3D(:,:,1))
            if ( BaseLL > 1 ) then
            write(*,*) 'Updated emiss (IX,IY,L2): ', OUTARR_3D(IX,IY,2)
            write(*,*) 'Total level 2: ', SUM(OUTARR_3D(:,:,2))
            endif
         ENDIF

      ENDDO ! N

      ! Apply mask
      IF ( DO_MASK ) THEN

         ! testing only:
         IF ( debug .and. am_I_Root ) THEN
            write(*,*) 'Apply masks...'
         ENDIF

         ! Restrict mask values to a maximum of 1. Higher values are
         ! possible if multiple masks are used which overlap, so make
         ! sure that emissions from these grid boxes are not
         ! artificially increased.
         WHERE ( MASK > 1d0 ) MASK = 1d0

         ! Apply the mask. Make sure that emissions become negative
         ! for all non-mask areas. This is required to make sure that
         ! these grid boxes will be ignored when calculating the
         ! final emissions. 
         DO L = 1, BaseLL
            WHERE ( MASK(:,:,1) <= 0d0 )
               OUTARR_3D(:,:,L) = -999d0
            ELSEWHERE
               OUTARR_3D(:,:,L) = OUTARR_3D(:,:,L) * MASK(:,:,1) 
            ENDWHERE
         ENDDO

         ! testing only:
         IF ( debug .and. am_I_Root ) THEN
            write(*,*) 'Scale factor (IX,IY,L1): ', MASK(IX,IY,1)
            write(*,*) 'Updated emiss (IX,IY,L1): ', OUTARR_3D(IX,IY,1)
         ENDIF

      ENDIF

      ! Cleanup and leave w/ success
      ScalCnt  => NULL()
      RC = HCO_SUCCESS
 
      END SUBROUTINE GET_CURRENT_EMISSIONS
! EOC
      END MODULE HCO_CALC_MOD
!EOF
