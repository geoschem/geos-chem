!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: merra_a1_mod
!
! !DESCRIPTION: Module MERRA\_A1\_MOD contains subroutines for reading the 
!  1-hour time averaged (aka "A1") fields from the MERRA data archive.
!\\
!\\
! !INTERFACE: 
!
      MODULE MERRA_A1_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC  :: GET_MERRA_A1_FIELDS
      PUBLIC  :: OPEN_MERRA_A1_FIELDS
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: A1_CHECK
      PRIVATE :: DO_OPEN_A1
      PRIVATE :: READ_A1
!
! !REMARKS:
!  Don't bother with the file unzipping anymore.
!
! !REVISION HISTORY:
!  19 Aug 2010 - R. Yantosca - Initial version, based on a3_read_mod.f
!  03 Aug 2012 - R. Yantosca - Now make IU_A1 a private module variable
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      INTEGER :: N_A1_FIELDS    ! # of fields in the file
      INTEGER :: IU_A1          ! File LUN for A-1 fields

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_open_a1
!
! !DESCRIPTION: Function DO\_OPEN\_A1 returns TRUE if is time to open the A1 
!  met field file or FALSE otherwise.  This prevents us from opening a file 
!  which has already been opened. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION DO_OPEN_A1( NYMD, NHMS, RESET ) RESULT( DO_OPEN )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)           :: NYMD     ! YYYYMMDD and hhmmss to test
      INTEGER, INTENT(IN)           :: NHMS     !  if it's time to open file
      LOGICAL, INTENT(IN), OPTIONAL :: RESET    ! Reset the 
!
! !RETURN VALUE:
!
      LOGICAL                       :: DO_OPEN  ! =T if it's time to open file
!
! !REVISION HISTORY: 
!  19 Aug 2010 - R. Yantosca - Initial version, based on a3_read_mod.f
!  21 Sep 2010 - R. Yantosca - Add RESET via the argument list to reset
!                              the FIRST flag if so desired.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST    = .TRUE.
      INTEGER, SAVE :: LASTNYMD = -1
      INTEGER, SAVE :: LASTNHMS = -1
      
      !=================================================================
      ! DO_OPEN_A1 begins here!
      !=================================================================

      ! Reset the FIRST flag if necessary (i.e. if we have been
      ! reading A1 fields for MEGAN, then FIRST=.FALSE.).  This will
      ! allow us to start a simulation at hours other than 0 GMT
      ! (bmy, 9/21/10)
      IF ( PRESENT( RESET ) ) THEN
         IF ( RESET ) FIRST = .TRUE.
      ENDIF

      ! Initialize
      DO_OPEN = .FALSE.

      ! Return if we have already opened the file
      IF ( NYMD == LASTNYMD .and. NHMS == LASTNHMS ) THEN
         DO_OPEN = .FALSE. 
         GOTO 999
      ENDIF

      ! Open A1 file if it's 00:30 GMT,  or on the first call
      IF ( NHMS == 003000 .or. FIRST ) THEN
         DO_OPEN = .TRUE. 
         GOTO 999
      ENDIF

      !=================================================================
      ! Reset quantities for next call
      !=================================================================
 999  CONTINUE
      LASTNYMD = NYMD
      LASTNHMS = NHMS
      FIRST    = .FALSE.

      END FUNCTION DO_OPEN_A1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: open_merra_a1_fields
!
! !DESCRIPTION: Subroutine OPEN\_MERRA\_A1\_FIELDS opens the A1 met fields 
!  file for date NYMD and time NHMS. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OPEN_MERRA_A1_FIELDS( NYMD,  NHMS, Input_Opt, 
     &                                 RESET, IUNIT            )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_RES_EXT
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE FILE_MOD,           ONLY : FILE_EXISTS
      USE FILE_MOD,           ONLY : IOERROR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,           ONLY : EXPAND_DATE
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)            :: NYMD      ! YYYYMMDD date
      INTEGER,        INTENT(IN)            :: NHMS      ! hhmmss time
      TYPE(OptInput), INTENT(IN)            :: Input_Opt ! Input Options
      LOGICAL,        INTENT(IN),  OPTIONAL :: RESET     ! Reset first A1 flag?
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT), OPTIONAL :: IUNIT     ! Returns IU_A1
! 
! !REVISION HISTORY:
!  19 Aug 2010 - R. Yantosca - Initial version, based on a3_read_mod.f
!  03 Aug 2012 - R. Yantosca - Now use findFreeLUN to define IU_A1 locally
!  06 Aug 2012 - R. Yantosca - Add optional IUNIT to pass LUN to calling routine
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  11 Apr 2013 - R. Yantosca - Now pass directory fields with Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: DO_RESET
      LOGICAL            :: IT_EXISTS, IS_OPEN
      INTEGER            :: IOS
      CHARACTER(LEN=2)   :: DUM
      CHARACTER(LEN=8)   :: IDENT
      CHARACTER(LEN=255) :: A1_FILE
      CHARACTER(LEN=255) :: GEOS_DIR
      CHARACTER(LEN=255) :: PATH

      !=================================================================
      ! OPEN_MERRA_A1_FIELDS begins here!
      !=================================================================
      
      ! Define shadow variable for optional RESET switch
      IF ( PRESENT( RESET ) ) THEN
         DO_RESET = RESET
      ELSE
         DO_RESET = .FALSE.
      ENDIF
            
      ! Check if it's time to open file
      IF ( DO_OPEN_A1( NYMD, NHMS, DO_RESET ) ) THEN

         !---------------------------
         ! Initialization
         !---------------------------

         ! Strings for directory & filename
         GEOS_DIR = TRIM( Input_Opt%MERRA_DIR )
         A1_FILE  = 'YYYYMMDD.a1.' // GET_RES_EXT()

         ! Replace date tokens
         CALL EXPAND_DATE( A1_FILE,  NYMD, NHMS )
         CALL EXPAND_DATE( GEOS_DIR, NYMD, NHMS )

         ! Full file path
         PATH = TRIM( Input_Opt%DATA_DIR ) // 
     &          TRIM( GEOS_DIR           ) // 
     &          TRIM( A1_FILE            )

         ! Close previously opened A-1 file
	 INQUIRE( IU_A1, OPENED=IS_OPEN )
         IF ( IS_OPEN ) CLOSE( IU_A1 )

         ! Make sure the file unit is valid before we open the file
         IF ( .not. FILE_EXISTS( IU_A1 ) ) THEN
            CALL ERROR_STOP( 'Could not find file!', 
     &                       'OPEN_MERRA_A1_FIELDS (merra_a1_mod.f)' )
         ENDIF

         !---------------------------
         ! Open the A1 file
         !---------------------------

         ! Find a free file LUN
         IU_A1 = findFreeLUN()

         ! If IUNIT is specified, then pass the value of
         ! IU_A1 back to the calling program
         IF ( PRESENT( IUNIT ) ) THEN
            IUNIT = IU_A1
         ENDIF

         ! Open the file
         OPEN( UNIT   = IU_A1,         FILE   = TRIM( PATH ),
     &         STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &         FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_A1, 'open_merra_a1_fields:1' )
         ENDIF

         ! Echo info
         WRITE( 6, 100 ) TRIM( PATH ), IU_A1
 100     FORMAT( '     - Opening: ', a, ' on unit ', i4 )
         
         !---------------------------
         ! Get # of fields in file
         !---------------------------
         
         ! Read the IDENT string
         READ( IU_A1, IOSTAT=IOS ) IDENT

         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_A1, 'open_merra_a1_fields:2' )
         ENDIF

         ! The last 2 digits of the ident string
         ! is the # of fields contained in the file
         READ( IDENT(7:8), '(i2.2)' ) N_A1_FIELDS        
        
      ENDIF

      END SUBROUTINE OPEN_MERRA_A1_FIELDS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_merra_a1_fields 
!
! !DESCRIPTION: Subroutine GET\_MERRA\_A1\_FIELDS is a wrapper for routine
!  READ\_A1.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_MERRA_A1_FIELDS( NYMD, NHMS, State_Met )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: NYMD        ! YYYYMMDD and 
      INTEGER,        INTENT(IN)    :: NHMS        !  hhmmss of data to read
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
! 
! !REVISION HISTORY: 
!  19 Aug 2010 - R. Yantosca - Initial version, based on a3_read_mod.f
!  25 Aug 2010 - R. Yantosca - Now pass LWI down to READ_A1
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Remove reference to dao_mod.F module arrays
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, SAVE :: LASTNYMD = -1
      INTEGER, SAVE :: LASTNHMS = -1

      !=================================================================
      ! Initialization
      !=================================================================

      ! Skip over previously-read A-3 fields
      IF ( NYMD == LASTNYMD .and. NHMS == LASTNHMS ) THEN
         WRITE( 6, 100 ) NYMD, NHMS
 100     FORMAT( '     - A-3 met fields for NYMD, NHMS = ', 
     &           i8.8, 1x, i6.6, ' have been read already' ) 
         RETURN
      ENDIF

      !=================================================================      
      ! Read data from disk
      !=================================================================
      CALL READ_A1( NYMD     = NYMD, 
     &              NHMS     = NHMS,
     &              ALBEDO   = State_Met%ALBD,
     &              CLDTOT   = State_Met%CLDFRC,   
     &              EFLUX    = State_Met%EFLUX,    
     &              EVAP     = State_Met%EVAP,
     &              FRSEAICE = State_Met%FRSEAICE,
     &              FRSNO    = State_Met%FRSNO, 
     &              GRN      = State_Met%GRN,      
     &              GWETROOT = State_Met%GWETROOT,
     &              GWETTOP  = State_Met%GWETTOP, 
     &              HFLUX    = State_Met%HFLUX,  
     &              LAI      = State_Met%LAI,
     &              LWI      = State_Met%LWI,
     &              LWGNT    = State_Met%RADLWG,
     &              PARDF    = State_Met%PARDF,
     &              PARDR    = State_Met%PARDR,
     &              PBLH     = State_Met%PBLH,
     &              PRECANV  = State_Met%PRECANV,   
     &              PRECTOT  = State_Met%PRECTOT,  
     &              PRECCON  = State_Met%PRECCON, 
     &              PRECLSC  = State_Met%PRECLSC,   
     &              PRECSNO  = State_Met%PRECSNO,
     &              SEAICE00 = State_Met%SEAICE00,
     &              SEAICE10 = State_Met%SEAICE10,
     &              SEAICE20 = State_Met%SEAICE20,
     &              SEAICE30 = State_Met%SEAICE30,
     &              SEAICE40 = State_Met%SEAICE40,
     &              SEAICE50 = State_Met%SEAICE50,
     &              SEAICE60 = State_Met%SEAICE60,
     &              SEAICE70 = State_Met%SEAICE70,
     &              SEAICE80 = State_Met%SEAICE80,
     &              SEAICE90 = State_Met%SEAICE90,
     &              SLP      = State_Met%SLP,
     &              SNODP    = State_Met%SNODP,
     &              SNOMAS   = State_Met%SNOMAS,
     &              SWGNT    = State_Met%RADSWG,
     &              TROPPT   = State_Met%TROPP,
     &              T2M      = State_Met%TS,
     &              TS       = State_Met%TSKIN, 
     &              U10M     = State_Met%U10M,
     &              USTAR    = State_Met%USTAR,
     &              V10M     = State_Met%V10M,
     &              Z0M      = State_Met%Z0        )

      ! Save NYMD, NHMS for next call
      LASTNYMD = NYMD
      LASTNHMS = NHMS

      END SUBROUTINE GET_MERRA_A1_FIELDS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_a1
!
! !DESCRIPTION: Subroutine READ\_A1 reads MERRA 1-hour time averaged ("A1") 
!  met fields from disk.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_A1( NYMD,     NHMS, 
     &                    ALBEDO,   CLDTOT,   EFLUX,    EVAP,    
     &                    FRSEAICE, FRSNO,    GRN,      GWETROOT, 
     &                    GWETTOP,  HFLUX,    LAI,      LWGNT,   
     &                    LWI,      PARDF,    PARDR,    PBLH,     
     &                    PRECANV,  PRECTOT,  PRECCON,  PRECLSC,  
     &                    PRECSNO,  SEAICE00, SEAICE10, SEAICE20, 
     &                    SEAICE30, SEAICE40, SEAICE50, SEAICE60, 
     &                    SEAICE70, SEAICE80, SEAICE90, SLP,      
     &                    SNODP,    SNOMAS,   SWGNT,    TROPPT,   
     &                    T2M,      TS,       U10M,     USTAR,    
     &                    V10M,     Z0M                           )
!
! !USES:
!
      USE DIAG_MOD,     ONLY : AD67
      USE FILE_MOD,     ONLY : IOERROR
      USE TIME_MOD,     ONLY : SET_CT_A1
      USE TIME_MOD,     ONLY : TIMESTAMP_STRING
      USE TRANSFER_MOD, ONLY : TRANSFER_2D
      USE TRANSFER_MOD, ONLY : TRANSFER_TO_1D

      USE CMN_SIZE_MOD                             ! Size parameters
      USE CMN_DIAG_MOD                             ! ND67 flag
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: NYMD                   ! YYYYMMDD and hhmmss     
      INTEGER, INTENT(IN)  :: NHMS                   !  of data to read
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: ALBEDO  (IIPAR,JJPAR)  ! Sfc albedo [unitless]
      REAL*8,  INTENT(OUT) :: CLDTOT  (IIPAR,JJPAR)  ! Column cld fraction
      REAL*8,  INTENT(OUT) :: EFLUX   (IIPAR,JJPAR)  ! Latent heat flux [W/m2]
      REAL*8,  INTENT(OUT) :: EVAP    (IIPAR,JJPAR)  ! Surface evap [kg/m2/s]
      REAL*8,  INTENT(OUT) :: FRSEAICE(IIPAR,JJPAR)  ! Sfc sea ice fraction
      REAL*8,  INTENT(OUT) :: FRSNO   (IIPAR,JJPAR)  ! Sfc snow fraction
      REAL*8,  INTENT(OUT) :: GRN     (IIPAR,JJPAR)  ! Greenness fraction
      REAL*8,  INTENT(OUT) :: GWETROOT(IIPAR,JJPAR)  ! Root soil wetness [frac]
      REAL*8,  INTENT(OUT) :: GWETTOP (IIPAR,JJPAR)  ! Topsoil wetness [frac]
      REAL*8,  INTENT(OUT) :: HFLUX   (IIPAR,JJPAR)  ! Sensible H-flux [W/m2]
      REAL*8,  INTENT(OUT) :: LAI     (IIPAR,JJPAR)  ! Leaf area index [m2/m2]
      REAL*8,  INTENT(OUT) :: LWI     (IIPAR,JJPAR)  ! Leaf area index [m2/m2]
      REAL*8,  INTENT(OUT) :: LWGNT   (IIPAR,JJPAR)  ! Net LW rad @ sfc [W/m2]
      REAL*8,  INTENT(OUT) :: PARDF   (IIPAR,JJPAR)  ! Diffuse PAR [W/m2]
      REAL*8,  INTENT(OUT) :: PARDR   (IIPAR,JJPAR)  ! Direct PAR [W/m2]
      REAL*8,  INTENT(OUT) :: PBLH    (IIPAR,JJPAR)  ! PBL height [m]
      REAL*8,  INTENT(OUT) :: PRECANV (IIPAR,JJPAR)  ! Anv prec @ sfc [kg/m2/s]
      REAL*8,  INTENT(OUT) :: PRECTOT (IIPAR,JJPAR)  ! Tot prec @ sfc [kg/m2/s]
      REAL*8,  INTENT(OUT) :: PRECCON (IIPAR,JJPAR)  ! CV prec @ sfc [kg/m2/s]
      REAL*8,  INTENT(OUT) :: PRECLSC (IIPAR,JJPAR)  ! LS prec @ sfc [kg/m2/s]
      REAL*8,  INTENT(OUT) :: PRECSNO (IIPAR,JJPAR)  ! Snow precip [kg/m2/s]
      REAL*8,  INTENT(OUT) :: SEAICE00(IIPAR,JJPAR)  ! Sea ice coverage 00-10%
      REAL*8,  INTENT(OUT) :: SEAICE10(IIPAR,JJPAR)  ! Sea ice coverage 10-20%
      REAL*8,  INTENT(OUT) :: SEAICE20(IIPAR,JJPAR)  ! Sea ice coverage 20-30%
      REAL*8,  INTENT(OUT) :: SEAICE30(IIPAR,JJPAR)  ! Sea ice coverage 30-40%
      REAL*8,  INTENT(OUT) :: SEAICE40(IIPAR,JJPAR)  ! Sea ice coverage 40-50%
      REAL*8,  INTENT(OUT) :: SEAICE50(IIPAR,JJPAR)  ! Sea ice coverage 50-60%
      REAL*8,  INTENT(OUT) :: SEAICE60(IIPAR,JJPAR)  ! Sea ice coverage 60-70%
      REAL*8,  INTENT(OUT) :: SEAICE70(IIPAR,JJPAR)  ! Sea ice coverage 70-80%
      REAL*8,  INTENT(OUT) :: SEAICE80(IIPAR,JJPAR)  ! Sea ice coverage 80-90%
      REAL*8,  INTENT(OUT) :: SEAICE90(IIPAR,JJPAR)  ! Sea ice coverage 90-100%
      REAL*8,  INTENT(OUT) :: SLP     (IIPAR,JJPAR)  ! Sea level pressure [hPa]
      REAL*8,  INTENT(OUT) :: SNODP   (IIPAR,JJPAR)  ! Snow depth [m]
      REAL*8,  INTENT(OUT) :: SNOMAS  (IIPAR,JJPAR)  ! Snow mass [kg/m2]
      REAL*8,  INTENT(OUT) :: SWGNT   (IIPAR,JJPAR)  ! SW rad @ sfc [W/m2] 
      REAL*8,  INTENT(OUT) :: TROPPT  (IIPAR,JJPAR)  ! T'pause pressure [hPa]
      REAL*8,  INTENT(OUT) :: T2M     (IIPAR,JJPAR)  ! T @ 2m height [K]
      REAL*8,  INTENT(OUT) :: TS      (IIPAR,JJPAR)  ! Sfc skin T [K]
      REAL*8,  INTENT(OUT) :: U10M    (IIPAR,JJPAR)  ! U-wind @ 10m [m/s]
      REAL*8,  INTENT(OUT) :: USTAR   (IIPAR,JJPAR)  ! Friction velocity [m/s]
      REAL*8,  INTENT(OUT) :: V10M    (IIPAR,JJPAR)  ! V-wind @ 10m [m/s]
      REAL*8,  INTENT(OUT) :: Z0M     (IIPAR,JJPAR)  ! Roughness height [m]
! 
! !REVISION HISTORY: 
!  19 Aug 2010 - R. Yantosca - Initial version, based on a3_read_mod.f
!  25 Aug 2010 - R. Yantosca - Now read LWI (land/water/ice) from disk
!  15 Aug 2011 - R. Yantosca - Now save SWGDN in 2nd slot of ND67 diagnostic
!  25 Mar 2011 - R. Yantosca - Bug fix: make local SWGDN array for ND67 diag
!  08 Jun 2012 - S. Philip   - Correction for MERRA boundary layer height
!  03 Aug 2012 - R. Yantosca - Now use locally-defined IU_A1 file LUN
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER           :: I, IOS, J, NFOUND, XYMD, XHMS
      CHARACTER(LEN=8)  :: NAME
      CHARACTER(LEN=16) :: STAMP

      ! Arrays
      REAL*4            :: Q2   (IIPAR,JJPAR)
      REAL*4            :: SWGDN(IIPAR,JJPAR)

      !=================================================================
      ! READ_A1 begins here!      
      !=================================================================

      ! Zero the number of A1 fields that we have found
      NFOUND = 0

      !=================================================================
      ! Read the A1 fields from disk
      !=================================================================
      DO

         ! Read the A1 field name
         READ( IU_A1, IOSTAT=IOS ) NAME

         ! End of file test -- make sure we have found all fields
         IF ( IOS < 0 ) THEN
            CALL A1_CHECK( NFOUND, N_A1_FIELDS )
            EXIT
         ENDIF

         ! IOS > 0: True I/O error; stop w/ err msg
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:1' )

         ! CASE statement for A1 fields
         SELECT CASE ( TRIM( NAME ) )

            !-------------------------------------
            ! ALBEDO: surface albedo
            !-------------------------------------
            CASE ( 'ALBEDO' ) 
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:2' )

               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, ALBEDO )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! CLDTOT: column cloud fraction
            !--------------------------------------
            CASE ( 'CLDTOT' )
               READ ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:3' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, CLDTOT )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! EFLUX: latent heat flux
            !--------------------------------------
            CASE ( 'EFLUX' )
               READ ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:4' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, EFLUX )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! EVAP: surface evaporation
            !--------------------------------------
            CASE ( 'EVAP' )
               READ ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:5' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, EVAP )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! FRSEAICE: sea ice fraction
            !--------------------------------------
            CASE ( 'FRSEAICE' )
               READ ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:6' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, FRSEAICE )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! FRSNO: snow fraction
            !--------------------------------------
            CASE ( 'FRSNO' )
               READ ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:7' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, FRSNO )
                  NFOUND = NFOUND + 1
               ENDIF

            !--------------------------------------
            ! GRN: evapotranspiration flux
            !--------------------------------------
            CASE ( 'GRN' )
               READ ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:8' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, GRN )
                  NFOUND = NFOUND + 1
               ENDIF

            !-------------------------------------
            ! GWETROOT: root soil wetness
            !-------------------------------------
            CASE ( 'GWETROOT' ) 
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:9' )

               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, GWETROOT )
                  NFOUND = NFOUND + 1
               ENDIF

            !-------------------------------------
            ! GWETTOP: topsoil wetness 
            !-------------------------------------
            CASE ( 'GWETTOP' ) 
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:10' )

               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, GWETTOP )
                  NFOUND = NFOUND + 1
               ENDIF

            !-------------------------------------
            ! HFLUX: sensible heat flux 
            !-------------------------------------
            CASE ( 'HFLUX' ) 
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:11' )

               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, HFLUX )
                  NFOUND = NFOUND + 1
               ENDIF

            !-------------------------------------
            ! LAI: GMAO leaf area index
            !-------------------------------------
            CASE ( 'LAI' ) 
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:12' )

               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, LAI )
                  NFOUND = NFOUND + 1
               ENDIF

            !-------------------------------------
            ! LWI: land/water/ice flags
            ! (used for backwards compatibility)
            !-------------------------------------
            CASE ( 'LWI' ) 
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:13' )

               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, LWI )
                  NFOUND = NFOUND + 1
               ENDIF

            !-------------------------------------
            ! LWGNT: net LW radiation @ ground
            !-------------------------------------
            CASE ( 'LWGNT' ) 
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:14' )

               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, LWGNT )
                  NFOUND = NFOUND + 1
               ENDIF

            !-------------------------------------
            ! PARDF: photosyn active diff rad
            !-------------------------------------
            CASE ( 'PARDF' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:15' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, PARDF )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! PARDR: photosyn active direct rad
            !-------------------------------------
            CASE ( 'PARDR' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:16' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, PARDR )
                  NFOUND = NFOUND + 1
               ENDIF 

            !-------------------------------------
            ! PBLH: boundary layer height
            !-------------------------------------
            CASE ( 'PBLH' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:17' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, PBLH )
                  NFOUND = NFOUND + 1
               ENDIF       

            !-------------------------------------
            ! PRECANV: anvil precip @ ground
            !-------------------------------------
            CASE ( 'PRECANV' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:18' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, PRECANV )
                  NFOUND = NFOUND + 1
               ENDIF       

            !-------------------------------------
            ! PRECCON: convective precip @ ground
            !-------------------------------------
            CASE ( 'PRECCON' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:19' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, PRECCON )
                  NFOUND = NFOUND + 1
               ENDIF   

            !-------------------------------------
            ! PRECLSC: LS precip @ ground
            !-------------------------------------
            CASE ( 'PRECLSC' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:20' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, PRECLSC )
                  NFOUND = NFOUND + 1
               ENDIF   

            !-------------------------------------
            ! PRECTOT: total precip @ ground
            !-------------------------------------
            CASE ( 'PRECTOT' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:21' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, PRECTOT )
                  NFOUND = NFOUND + 1
               ENDIF   

            !-------------------------------------
            ! PRECSNO: snow precip @ ground
            !-------------------------------------
            CASE ( 'PRECSNO' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:22' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, PRECSNO )
                  NFOUND = NFOUND + 1
               ENDIF   

            !-------------------------------------
            ! SEAICE00: Sea ice bin 0-10%
            !-------------------------------------
            CASE ( 'SEAICE00' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:23' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SEAICE00 )
                  NFOUND = NFOUND + 1
               ENDIF   

            !-------------------------------------
            ! SEAICE10: Sea ice bin 10-20%
            !-------------------------------------
            CASE ( 'SEAICE10' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:24' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SEAICE10 )
                  NFOUND = NFOUND + 1
               ENDIF   

            !-------------------------------------
            ! SEAICE20: Sea ice bin 20-30%
            !-------------------------------------
            CASE ( 'SEAICE20' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:25' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SEAICE20 )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! SEAICE30: Sea ice bin 30-40%
            !-------------------------------------
            CASE ( 'SEAICE30' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:26' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SEAICE30 )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! SEAICE40: Sea ice bin 40-50%
            !-------------------------------------
            CASE ( 'SEAICE40' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:27' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SEAICE40 )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! SEAICE50: Sea ice bin 50-60%
            !-------------------------------------
            CASE ( 'SEAICE50' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:28' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SEAICE50 )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! SEAICE60: Sea ice bin 60-70%
            !-------------------------------------
            CASE ( 'SEAICE60' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:29' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SEAICE60 )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! SEAICE70: Sea ice bin 70-80%
            !-------------------------------------
            CASE ( 'SEAICE70' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:30' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SEAICE70 )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! SEAICE80: Sea ice bin 80-90%
            !-------------------------------------
            CASE ( 'SEAICE80' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:31' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SEAICE80 )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! SEAICE90: Sea ice bin 90-100%
            !-------------------------------------
            CASE ( 'SEAICE90' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:32' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SEAICE90 )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! SLP: sea level pressure
            !-------------------------------------
            CASE ( 'SLP' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:33' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SLP )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! SNODP: snow depth
            !-------------------------------------
            CASE ( 'SNODP' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:34' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SNODP )
                  NFOUND = NFOUND + 1
               ENDIF 
 
            !-------------------------------------
            ! SNOMAS: snow mass
            !-------------------------------------
            CASE ( 'SNOMAS' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:35' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SNOMAS )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! SWGNT: Net SW radiation @ ground
            !-------------------------------------
            CASE ( 'SWGNT' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:36' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, SWGNT )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! TROPPT: Troopause pressure, T-based
            !-------------------------------------
            CASE ( 'TROPPT' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:37' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, TROPPT )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! TS: Surface skin temperature
            !-------------------------------------
            CASE ( 'TS' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:38' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, TS )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! T2M: Temp @ 2m altitude
            !-------------------------------------
            CASE ( 'T2M' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:39' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, T2M )
                  NFOUND = NFOUND + 1
               ENDIF  

            !-------------------------------------
            ! U10M: U-wind @ 10m altitude
            !-------------------------------------
            CASE ( 'U10M' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:40' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, U10M )
                  NFOUND = NFOUND + 1
               ENDIF 

            !-------------------------------------
            ! USTAR: Friction velocity
            !-------------------------------------
            CASE ( 'USTAR' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:41' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, USTAR )
                  NFOUND = NFOUND + 1
               ENDIF 

            !-------------------------------------
            ! V10M: V-wind @ 10m altitude
            !-------------------------------------
            CASE ( 'V10M' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:42' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, V10M )
                  NFOUND = NFOUND + 1
               ENDIF 

            !-------------------------------------
            ! Z0M: roughness height
            !-------------------------------------
            CASE ( 'Z0M' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:43' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, Z0M )
                  NFOUND = NFOUND + 1
               ENDIF 

            !-------------------------------------
            ! Save SWGDN locally for ND67 diag
            !-------------------------------------
            CASE ( 'SWGDN' )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:44' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  SWGDN  = Q2
                  NFOUND = NFOUND + 1
               ENDIF

            !-------------------------------------
            ! Skip over these fields for now:
            ! LWTUP, QV2M
            !-------------------------------------
            CASE ( 'LWTUP',  'QV2M'  )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:45' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  NFOUND = NFOUND + 1
               ENDIF

            !------------------------------------------------
            ! Field not found -- skip over
            !------------------------------------------------
            CASE DEFAULT
               WRITE ( 6, 200 )
               READ( IU_A1, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_A1, 'read_a1:46' )

         END SELECT
            
         !==============================================================
         ! If we have found all the fields for this time, then exit 
         ! the loop.  Otherwise, go on to the next iteration.
         !==============================================================
         IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
            IF ( NFOUND == N_A1_FIELDS ) THEN 
               STAMP = TIMESTAMP_STRING( NYMD, NHMS )
               WRITE( 6, 210 ) NFOUND, STAMP, IU_A1
               EXIT
            ENDIF
         ENDIF
      ENDDO

      ! FORMATs
 200  FORMAT( 'Searching for next MERRA A1 field!'                    )
 210  FORMAT( '     - Found all ', i3, ' MERRA A1 met fields for ', 
     &             a, ' on unit ', i4                                 )

      !=================================================================
      !        %%%%% SPECIAL HANDLING FOR CERTAIN FIELDS %%%%% 
      !
      ! In MERRA, the PRECTOT etc. surface precipitation met fields
      ! fields have units of [kg/m2/s].  In all other GEOS versions, 
      ! PREACC and PRECON have units of [mm/day].  
      !
      ! Therefore, for backwards compatibility with existing code, 
      ! apply the following unit conversion to the GEOS-5 PRECTOT and
      ! PRECCON fields:
      !
      !
      !     kg  |    m3    | 86400 s | 1000 mm
      !   ------+----------+---------+--------- = 86400 
      !    m2 s |  1000 kg |  day    |   m
      !              ^
      !              |
      !       1 / density of water 
      !=================================================================
      
      ! Convert from [kg/m2/s] --> [mm/day]
      PRECANV = PRECANV * 86400d0
      PRECCON = PRECCON * 86400d0
      PRECLSC = PRECLSC * 86400d0
      PRECTOT = PRECTOT * 86400d0

      ! Correction for MERRA PBLH:
      ! PBLH must be greater than a minimum mechanical mixing depth,
      ! defined as 700*friction velocity
      ! (Koracin and Berkowicz, 1988; Lin and McElroy, 2010)
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         PBLH(I,J) = max( PBLH(I,J), 700d0*USTAR(I,J) )
      ENDDO
      ENDDO

      !=================================================================
      ! ND67 diagnostic: A1 surface fields
      !=================================================================
      IF ( ND67 > 0 ) THEN
         AD67(:,:,1 ) = AD67(:,:,1 ) + HFLUX    ! Sensible heat flux [W/m2]
         AD67(:,:,2 ) = AD67(:,:,2 ) + SWGDN    ! Incident SW rad @ sfc [W/m2]
         AD67(:,:,3 ) = AD67(:,:,3 ) + PRECTOT  ! Tot prec @ sfc [kg/m2/s]
         AD67(:,:,4 ) = AD67(:,:,4 ) + PRECCON  ! CV prec @ sfc [kg/m2/s]
         AD67(:,:,5 ) = AD67(:,:,5 ) + T2M      ! T @ 2m height [K]
         AD67(:,:,6 ) = AD67(:,:,6 ) + 0e0      !
         AD67(:,:,7 ) = AD67(:,:,7 ) + USTAR    ! Friction velocity [m/s]
         AD67(:,:,8 ) = AD67(:,:,8 ) + Z0M      ! Roughness height [m]
         AD67(:,:,9 ) = AD67(:,:,9 ) + PBLH     ! PBL height [m]
         AD67(:,:,10) = AD67(:,:,10) + CLDTOT   ! Column cld fraction
         AD67(:,:,11) = AD67(:,:,11) + U10M     ! U-wind @ 10m [m/s]
         AD67(:,:,12) = AD67(:,:,12) + V10M     ! V-wind @ 10m [m/s]
         AD67(:,:,14) = AD67(:,:,14) + ALBEDO   ! Sfc albedo [unitless]
         AD67(:,:,17) = AD67(:,:,17) + TROPPT   ! T'pause pressure [hPa]
         AD67(:,:,18) = AD67(:,:,18) + SLP      ! Sea level pressure [hPa]
         AD67(:,:,19) = AD67(:,:,19) + TS       ! Sfc skin temperature [K]
         AD67(:,:,20) = AD67(:,:,20) + PARDF    ! Diffuse PAR [W/m2]
         AD67(:,:,21) = AD67(:,:,21) + PARDR    ! Direct PAR [W/m2]
         AD67(:,:,22) = AD67(:,:,22) + GWETTOP  ! Topsoil wetness [frac]
         AD67(:,:,23) = AD67(:,:,23) + EFLUX    ! Latent heat flux [W/m2]
      ENDIF

      !=================================================================
      ! Cleanup and quit
      !=================================================================

      ! Increment the # of times READ_A1 is called
      CALL SET_CT_A1( INCREMENT=.TRUE. )

      END SUBROUTINE READ_A1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: a1_check
!
! !DESCRIPTION: Subroutine A1\_CHECK prints an error message if not all of 
!  the A-3 met fields are found.  The run is also terminated. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE A1_CHECK( NFOUND, N_A1 )
!
! !USES:
!
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: NFOUND   ! Number of met fields read in from disk
      INTEGER, INTENT(IN) :: N_A1     ! Number of expected met fields
! 
! !REVISION HISTORY: 
!  19 Aug 2010 - R. Yantosca - Initial version, based on a3_read_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Test if N_FOUND == N_A1
      IF ( NFOUND /= N_A1 ) THEN
         
         ! Write error msg
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100   ) 
         WRITE( 6, 110   ) N_A1, NFOUND
         WRITE( 6, 120   )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! FORMATs
 100     FORMAT( 'ERROR -- not enough MERRA A1 fields found!' )
 110     FORMAT( 'There are ', i2, ' fields but only ', i2 ,
     &           ' were found!'                               )
 120     FORMAT( '### STOP in A1_CHECK (merra_a1_mod.f)'      )

         ! Deallocate arrays and stop
         CALL GEOS_CHEM_STOP

      ENDIF

      END SUBROUTINE A1_CHECK
!EOC
      END MODULE MERRA_A1_MOD
