!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: unitconv_mod
!
! !DESCRIPTION: Module UNITCONV\_MOD contains routines to convert data
! to GEOS-Chem units.
!
! \\
! !INTERFACE:
!
      MODULE UNITCONV_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CHANGE_UNITS
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: IsInWord
!
! !PRIVATE MODULE VARIABLES:
!
      REAL*8,  PARAMETER   :: N_0             = 6.022d+23
      REAL*8,  PARAMETER   :: SEC_IN_DAY      = 86400d0
      REAL*8,  PARAMETER   :: SEC_IN_LEAPYEAR = SEC_IN_DAY * 366d0
      REAL*8,  PARAMETER   :: SEC_IN_REGYEAR  = SEC_IN_DAY * 365d0
!
! !REVISION HISTORY:
!  15 May 2012 - C. Keller: Initialization
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
! !IROUTINE: change_units
!
! !DESCRIPTION: Subroutine CHANGE\_UNITS converts the values of the 
! passed array to GEOS-Chem tracer units. This is molec/cm2/s for fluxes and
! molec/cm3 for concentrations. Hydrocarbons are treated as molec carbon
! (molec C).\\ 
! The mass and area/volume conversion is always performed, but the time 
! conversion is only done if a valid time string is provided.
! For example, if the input unit is kg/m3 it will be converted to molec/cm3, 
! and ug/m2/year will be converted to molec/cm2/s. 
! If no (valid) area/volume is given in the unit string, the data is 
! normalized by the grid box areas calculated from the specified ncfile,
! i.e. an input unit of kg will be converted to molec/cm2!!
! The following masses are supported: ng, ug, mg, g, kg, atom(s),
! molec(ules), nmol(s), umol(s), mmol(s), mol(s),  
! The following times are supported: /s, s-1, /d, d-1, /day, /mt, /month, 
! /y, y-1, /yr, /year.
! The following volumes/areas are supported: /cm2, /m2, /cm3, /m3, /L.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHANGE_UNITS ( ARRAY,  UNITS,  XNUMOL,    & 
                                Tr_Sp,  NCFILE, YYYY,  MM )
!
! !USES:
!
      USE ERROR_MOD,                ONLY : ERROR_STOP
      USE CHARPAK_MOD,              ONLY : TRANLC
      USE CMN_GCTM_MOD,             ONLY : PI, Re
!<<MSL>>      USE NCDF_MOD,                 ONLY : NC_READ_GRID
      USE TIME_MOD,                 ONLY : ITS_A_LEAPYEAR
!
! !ARGUMENTS:
!
      ! Emission field 
      REAL*8,           POINTER         :: ARRAY(:,:,:,:)

      ! Current units
      CHARACTER(LEN=*), INTENT(IN )     :: UNITS 

      ! Molecules of tracer per kg tracer
      REAL*8,           INTENT(IN )     :: XNUMOL

      ! Molecules of tracer per molecules of species
      ! (e.g. 3 for VOC with 3 carbon atoms if tracer is carbon) 
      REAL*8,           INTENT(IN )     :: Tr_Sp

      ! Data year and month
      INTEGER,          INTENT(IN )     :: YYYY
      INTEGER,          INTENT(IN )     :: MM

      ! ncdf source file (for grid box areas)
      CHARACTER(LEN=*), INTENT(IN )     :: NCFILE 
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller - Initial Version
!  23 May 2013 - C. Keller - Now use additive method
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
      REAL*8                :: Tr_Kg, Coef1 
      REAL*8                :: Fact
      REAL*8                :: RLAT, AREA
      CHARACTER(LEN=255)    :: MSG, LOC
      CHARACTER(LEN=31 )    :: unt
      INTEGER               :: Year, Month
      INTEGER               :: J, tmpID
      INTEGER               :: MONTHDAYS(12) = (/ 31, 28, 31, 30, &
                                                  31, 30, 31, 31, & 
                                                  30, 31, 30, 31 /)

      INTEGER               :: NLON,  NLAT
      REAL*8, ALLOCATABLE   :: XEDGE(:), YSIN(:)

      !=================================================================
      ! CHANGE_UNIT begins here
      !=================================================================

      ! Init
      LOC    = 'CHANGE_UNIT (unitconv_mod.F90)'
      Fact   = 1.0d0
      Coef1  = 1.0d0

      ! Get field unit and translate to lower case letters
      unt = TRIM(UNITS)
      CALL TRANLC( unt )        
      
      ! Tracer specific conversion factors:
      Month = MAX(MM,1)
      Year  = MAX(YYYY,1951)

      ! Molecules of tracer per kg
      Tr_Kg = XNUMOL

      !=================================================================
      ! Convert mass to molec tracer
      !=================================================================

      ! molecules atoms of carbon 
      IF ( IsInWord(unt,'molecc'  ) .OR. IsInWord(unt,'molec c') .OR. &
           IsInWord(unt,'atomc'   ) .OR. IsInWord(unt,'atom c' ) .OR. &
           IsInWord(unt,'atomsc'  ) .OR. IsInWord(unt,'atoms c') .OR. &
           IsInWord(unt,'molec(c)') .OR. IsInWord(unt,'molec tracer') &
         ) THEN
         Fact = Fact

      ! molecules / atoms of species
      ELSEIF ( IsInWord(unt,'molec') .OR. IsInWord(unt,'atom') ) THEN
         Fact = Fact / Tr_Sp 

      ! mols of species
      ELSEIF ( IsInWord(unt,'nmol') ) THEN
         Fact = Fact * 1.0d-9 * N_0 * Tr_Sp 
      ELSEIF ( IsInWord(unt,'umol') ) THEN
         Fact = Fact * 1.0d-6 * N_0 * Tr_Sp 
      ELSEIF ( IsInWord(unt,'mmol') ) THEN
         Fact = Fact * 1.0d-3 * N_0 * Tr_Sp 
      ELSEIF ( IsInWord(unt,'mol') ) THEN
         Fact = Fact * N_0 * Tr_Sp 

      ! mass Carbon of species
      ELSEIF ( IsInWord(unt,'kgc') .OR. IsInWord(unt,'kg c') .OR. &
               IsInWord(unt,'kg(c)') ) THEN
         Fact = Fact * Tr_Kg
      ELSEIF ( IsInWord(unt,'gc') .OR. IsInWord(unt,'g c') .OR. &
               IsInWord(unt,'g(c)') ) THEN
         Fact = Fact * 1.0d-3 * Tr_Kg

      ! mass of species
      ELSEIF ( IsInWord(unt,'ng') ) THEN
         Fact = Fact * 1.0d-12 * Tr_Kg * Tr_Sp 
      ELSEIF ( IsInWord(unt,'ug') ) THEN
         Fact = Fact * 1.0d-9 * Tr_Kg * Tr_Sp 
      ELSEIF ( IsInWord(unt,'mg') ) THEN
         Fact = Fact * 1.0d-6 * Tr_Kg * Tr_Sp 
      ELSEIF ( IsInWord(unt,'kg') ) THEN
         Fact = Fact * Tr_Kg * Tr_Sp 
      ELSEIF ( IsInWord(unt,'g') ) THEN
         Fact = Fact * 1.0d-3 * Tr_Kg * Tr_Sp 

      ! Stop run with error if none of the above specified mass 
      ! units is provided
      ELSE
         MSG = 'unrecognized mass unit: ' // TRIM(unt)
         CALL ERROR_STOP ( MSG, LOC )
      ENDIF

      !=================================================================
      ! Convert time to sec 
      !=================================================================

      ! second
      IF ( IsInWord(unt,'/s') .OR. IsInWord(unt,'s-1') .OR. &
           IsInWord(unt,'sec') ) THEN
         Fact = Fact

      ! hour
      ELSEIF ( IsInWord(unt,'hr') .OR. IsInWord(unt,'hour') ) THEN
         Fact = Fact / 3600.0d0

      ! day
      ELSEIF ( IsInWord(unt,'/d') .OR. IsInWord(unt,'-1d') .OR. &
               IsInWord(unt,'day') ) THEN
         Fact = Fact / SEC_IN_DAY

      ! month
      ELSEIF ( IsInWord(unt,'mt') .OR. IsInWord(unt,'month') ) THEN

         IF ( ITS_A_LEAPYEAR ( YEAR_IN=Year, FORCE=.TRUE. ) ) THEN
            MONTHDAYS(2) = 29
         ENDIF
         Coef1 = MONTHDAYS(Month) * SEC_IN_DAY 
         Fact  = Fact / Coef1
      
      ! year
      ELSEIF ( IsInWord(unt,'/y') .OR. IsInWord(unt,'y-1') .OR. &
               IsInWord(unt,'yr') .OR. IsInWord(unt,'year') ) THEN

         IF ( ITS_A_LEAPYEAR ( YEAR_IN=Year, FORCE=.TRUE. ) ) THEN
            Coef1 = SEC_IN_LEAPYEAR
         ELSE
            Coef1 = SEC_IN_REGYEAR
         ENDIF
         Fact = Fact / Coef1
      ENDIF

      !=================================================================
      ! Convert area to cm2 
      !=================================================================

      ! cm2
      IF ( IsInWord(unt,'/cm2' ) .OR. IsInWord(unt,'cm-2') .OR. &
           IsInWord(unt,'/cm^2') ) THEN 
         Fact = Fact

      ! m2
      ELSEIF ( IsInWord(unt,'/m2' ) .OR. IsInWord(unt,'m-2') .OR. & 
               IsInWord(unt,'/m^2') ) THEN 
         Fact = Fact * 1.0d-4

      !=================================================================
      ! Convert volume to cm3 
      !=================================================================

      ! cm3
      ELSEIF ( IsInWord(unt,'/cm3') .OR. IsInWord(unt,'cm-3') ) THEN 
         Fact = Fact

      ! dm3
      ELSEIF ( IsInWord(unt,'/dm3') .OR. IsInWord(unt,'dm-3') ) THEN 
         Fact = Fact * 1.0d-3

      ! m3
      ELSEIF ( IsInWord(unt,'/m3') .OR. IsInWord(unt,'m-3') ) THEN 
         Fact = Fact * 1.0d-6

      ! L
      ELSEIF ( IsInWord(unt,'/l') .OR. IsInWord(unt,'l-1') ) THEN 
         Fact = Fact * 1.0d-3

      ! If no area/volume conversion done up to here, normalize data by
      ! each grid box unit!!
      ELSE
         ! Get grid dimensions
         NLON = SIZE(ARRAY,1)
         NLAT = SIZE(ARRAY,2)

         ! Allocate grid edges
         ALLOCATE ( XEDGE(NLON+1) )
         ALLOCATE ( YSIN (NLAT+1) )
         XEDGE = 0d0
         YSIN  = 0d0

         ! Read input grid specifications
!<<MSL>>         CALL NC_READ_GRID( NLON, NLAT, TRIM(NCFILE), XEDGE, YSIN )
 
         ! Loop over all latitudes
         DO J = 1, NLAT

            ! get grid box area in cm2
            RLAT = ABS( YSIN(J+1) - YSIN(J) )
            AREA = ( 2d0 * PI * Re * RLAT * 1d4 * Re ) / DBLE( NLON )

            ! convert to cm-2
            ARRAY(:,J,:,:) = ARRAY(:,J,:,:) / AREA 
         ENDDO

         ! Prompt a warning
         WRITE(6, '(a  )' ) REPEAT( '#', 79 )
         WRITE(6, 100     ) TRIM( unt ), TRIM( NCFILE )
100      FORMAT ( 'WARNING: No area unit found - convert unit ', a, &
                  ' of file ', a, ' to /cm2!!' )

         ! Deallocate arrays
         DEALLOCATE ( XEDGE, YSIN )

      ENDIF

      ! Now apply correction factor
      ARRAY(:,:,:,:) = ARRAY(:,:,:,:) * Fact

      !=================================================================

      END SUBROUTINE CHANGE_UNITS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IsInWord
!
! !DESCRIPTION: Function IsInWord checks if the word InString
! contains the sequence of SearchString. 
!\\
! !INTERFACE:
!
      FUNCTION IsInWord ( InString, SearchString ) RESULT ( Cnt ) 
!
! !USES:
!
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*), INTENT(IN   )    :: InString
      CHARACTER(LEN=*), INTENT(IN   )    :: SearchString
      LOGICAL                            :: Cnt
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
      INTEGER  :: L1, L2, I

      !=================================================================
      ! IsInWord begins here
      !=================================================================

      ! Init
      L1 = MAX(1,LEN_TRIM(InString    ))
      L2 = MAX(1,LEN_TRIM(SearchString))
      Cnt = .FALSE.

      ! Return here if search string is larger than input word
      IF ( L2 > L1 ) RETURN

      ! Loop over all characters of input word and check if search string 
      ! matches the input word from the current position onwards.
      DO I = 1, L1 

         ! Check remaining word length
         IF ( (I+L2-1) > L1 ) EXIT 

         ! Compare strings
         IF ( InString(I:(I+L2-1)) == SearchString ) THEN
            Cnt = .TRUE.
            EXIT
         ENDIF
      ENDDO !I

      END FUNCTION IsInWord
!EOC
      END MODULE UNITCONV_MOD
!EOM
