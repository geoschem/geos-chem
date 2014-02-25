!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: merra_i6_mod
!
! !DESCRIPTION: Module MERRA\_I6\_MOD contains subroutines for reading the 
!  6-hour instantaneous (aka "I6") fields from the MERRA data archive.
!\\
!\\
! !INTERFACE: 
!
      MODULE MERRA_I6_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC  :: GET_MERRA_I6_FIELDS_1
      PUBLIC  :: GET_MERRA_I6_FIELDS_2
      PUBLIC  :: OPEN_MERRA_I6_FIELDS
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: I6_CHECK
      PRIVATE :: READ_I6
!
! !REMARKS:
!  Don't bother with the file unzipping anymore.
!
! !REVISION HISTORY:
!  19 Aug 2010 - R. Yantosca - Initial version, based on i6_read_mod.f
!  03 Aug 2012 - R. Yantosca - Now make IU_I6 a private module variable
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      INTEGER :: N_I6_FIELDS    ! # of fields in the file
      INTEGER :: IU_I6          ! File LUN for I-6 fields

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: open_merra_i6_fields
!
! !DESCRIPTION: Subroutine OPEN\_MERRA\_I6\_FIELDS opens the MERRA "I6" 
!  met fields file for date NYMD and time NHMS.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OPEN_MERRA_I6_FIELDS( NYMD, NHMS, Input_Opt )
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
      INTEGER,        INTENT(IN) :: NYMD        ! YYYYMMDD date
      INTEGER,        INTENT(IN) :: NHMS        ! hhmmss time
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
!
! !REVISION HISTORY: 
!  19 Aug 2010 - R. Yantosca - Initial version, based on i6_read_mod.f
!  03 Aug 2012 - R. Yantosca - Now use findFreeLUN to define IU_I6 locally
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  11 Apr 2013 - R. Yantosca - Now pass directory fields with Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST = .TRUE.
      LOGICAL            :: IT_EXISTS
      LOGICAL            :: IS_OPEN
      INTEGER            :: IOS, IUNIT
      CHARACTER(LEN=8)   :: IDENT
      CHARACTER(LEN=255) :: GEOS_DIR
      CHARACTER(LEN=255) :: I6_FILE, TP_FILE
      CHARACTER(LEN=255) :: PATH

      !=================================================================
      ! OPEN_MERRA_I6_FIELDS begins here!
      !=================================================================

      ! Check if it's time to open file
      IF ( NHMS == 000000 .or. FIRST ) THEN

         !---------------------------
         ! Initialization
         !---------------------------

         ! Strings for directory & filename
         GEOS_DIR = TRIM( Input_Opt%MERRA_DIR )
         I6_FILE  = 'YYYYMMDD.i6.' // GET_RES_EXT()

         ! Replace date tokens
         CALL EXPAND_DATE( GEOS_DIR, NYMD, NHMS )
         CALL EXPAND_DATE( I6_FILE,  NYMD, NHMS )

         ! Full file path
         PATH = TRIM( Input_Opt%DATA_DIR ) // 
     &          TRIM( GEOS_DIR           ) // 
     &          TRIM( I6_FILE            )

         ! Close previously opened A-3 file
	 INQUIRE( IU_I6, OPENED=IS_OPEN )
         IF ( IS_OPEN ) CLOSE( IU_I6 )

         ! Make sure the file unit is valid before we open it 
         IF ( .not. FILE_EXISTS( IU_I6 ) ) THEN 
            CALL ERROR_STOP( 'Could not find file!', 
     &                       'OPEN_MERRA_I6_FIELDS (merra_i6_mod.f)' )
         ENDIF

         !---------------------------
         ! Open the I6 file
         !---------------------------

         ! Find a free file LUN
         IU_I6 = findFreeLUN()

         ! Open the file
         OPEN( UNIT   = IU_I6,         FILE   = TRIM( PATH ),
     &         STATUS = 'OLD',         ACCESS = 'SEQUENTIAL',  
     &         FORM   = 'UNFORMATTED', IOSTAT = IOS )
               
         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_I6, 'open_merra_i6_fields:1' )
         ENDIF

         ! Echo info
         WRITE( 6, 100 ) TRIM( PATH ), IU_I6
 100     FORMAT( '     - Opening: ', a, ' on unit ', i4 )
         
         ! Set the proper first-time-flag false
         FIRST = .FALSE.

         !---------------------------
         ! Get # of fields in file
         !---------------------------

         ! Read the IDENT string
         READ( IU_I6, IOSTAT=IOS ) IDENT

         IF ( IOS /= 0 ) THEN
            CALL IOERROR( IOS, IU_I6, 'open_merra_i6_fields:2' )
         ENDIF
         
         ! The last 2 digits of the ident string
         ! is the # of fields contained in the file
         READ( IDENT(7:8), '(i2.2)' ) N_I6_FIELDS    

      ENDIF

      END SUBROUTINE OPEN_MERRA_I6_FIELDS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_merra_i6_fields_1
!
! !DESCRIPTION: Subroutine GET\_MERRA\_I6\_FIELDS\_1 is a wrapper for routine 
!  READ\_I6.  It reads the initial data at the start of a GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_MERRA_I6_FIELDS_1( NYMD, NHMS, State_Met )
!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: NYMD   ! YYYYMMDD date 
      INTEGER,        INTENT(IN)    :: NHMS   ! hhmmss time of desired data
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
! 
! !REVISION HISTORY: 
!  19 Aug 2010 - R. Yantosca - Initial version, based on i6_read_mod.f
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================      
      ! Read data from disk
      !=================================================================
      CALL READ_I6( NYMD = NYMD, 
     &              NHMS = NHMS, 
     &              PS   = State_Met%PS1, 
     &              RH   = State_Met%RH1   )

      END SUBROUTINE GET_MERRA_I6_FIELDS_1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_merra_i6_fields_2
!
! !DESCRIPTION: Subroutine GET\_MERRA\_I6\_FIELDS\_2 is a wrapper for routine 
!  READ\_I6.  It reads the data every 6 hours during a GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_MERRA_I6_FIELDS_2( NYMD, NHMS, State_Met )
!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: NYMD   ! YYYYMMDD date 
      INTEGER,        INTENT(IN)    :: NHMS   ! hhmmss time of desired data
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
! 
! !REVISION HISTORY: 
!  19 Aug 2010 - R. Yantosca - Initial version, based on i6_read_mod.f
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================      
      ! Read data from disk
      !=================================================================
      CALL READ_I6( NYMD = NYMD, 
     &              NHMS = NHMS,
     &              PS   = State_Met%PS2, 
     &              RH   = State_Met%RH2  )

      END SUBROUTINE GET_MERRA_I6_FIELDS_2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_i6
!
! !DESCRIPTION: Subroutine READ\_I6 reads GEOS-Chem I-6 (instantaneous 6-hour)
!  met fields from disk.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_I6( NYMD, NHMS, PS, RH )
!
! !USES:
!
      USE FILE_MOD,     ONLY : IOERROR
      USE TIME_MOD,     ONLY : SET_CT_I6
      USE TIME_MOD,     ONLY : TIMESTAMP_STRING
      USE TRANSFER_MOD, ONLY : TRANSFER_2D
      USE TRANSFER_MOD, ONLY : TRANSFER_3D

      USE CMN_SIZE_MOD                               ! Size parameters
      USE CMN_DIAG_MOD                               ! NDxx flags
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: NYMD                   ! YYYYMMDD and hhmmss
      INTEGER, INTENT(IN)  :: NHMS                   !  time of desired data
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: PS(IIPAR,JJPAR      )  ! Surface pressure [hPa]
      REAL*8,  INTENT(OUT) :: RH(IIPAR,JJPAR,LLPAR)  ! Rel. humidity    [1]
!
! !REVISION HISTORY: 
!  19 Aug 2010 - R. Yantosca - Initial version, based on i6_read_mod.f
!  03 Aug 2012 - R. Yantosca - Now use locally-defined IU_I6 file LUN
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  15 Nov 2013 - R. Yantosca - Now convert RH from [1] to [%], in order
!                              to be consistent with GEOS-Chem convention
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: IOS,  NFOUND, N_I6      
      INTEGER            :: XYMD, XHMS
      CHARACTER(LEN=8)   :: NAME
      CHARACTER(LEN=16)  :: STAMP

      ! Arrays
      REAL*4             :: Q2(IIPAR,JJPAR)
      REAL*4             :: Q3(IIPAR,JJPAR,LGLOB)

      !=================================================================
      ! READ_I6 begins here!
      !=================================================================

      ! Zero the number of I-6 fields we have already found
      NFOUND = 0

      !=================================================================
      ! Read the I-6 fields from disk
      !=================================================================
      DO 

         ! I-6 field name
         READ( IU_I6, IOSTAT=IOS ) NAME

         ! IOS < 0: End-of-file, but make sure we have 
         ! found all I-6 fields before exiting loop!
         IF ( IOS < 0 ) THEN
            CALL I6_CHECK( NFOUND, N_I6_FIELDS )
            EXIT
         ENDIF

         ! IOS > 0: True I/O error, stop w/ error msg
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:1' )

         ! CASE statement for met fields
         SELECT CASE ( TRIM( NAME ) )

            !------------------------------------
            ! PS: Surface pressure [hPa]
            !------------------------------------
            CASE ( 'PS' )
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, Q2
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:2' )

               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_2D( Q2, PS )
                  NFOUND = NFOUND + 1
               ENDIF

            !------------------------------------
            ! RH: Relative humidity [1]
            !     and converted to  [%]
            !------------------------------------
            CASE ( 'RH' )
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:3' )

               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  CALL TRANSFER_3D( Q3, RH )
                  NFOUND = NFOUND + 1

                  ! Convert RH from [1] to [%] 
                  RH = RH * 100d0
               ENDIF

            !------------------------------------
            ! PV, OMEGA: just skip over these
            !------------------------------------
            CASE ( 'PV', 'OMEGA' ) 
               READ( IU_I6, IOSTAT=IOS ) XYMD, XHMS, Q3
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_I6, 'read_i6:14' )
             
               IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
                  NFOUND = NFOUND + 1
               ENDIF

            !------------------------------------
            ! Field not found
            !------------------------------------
            CASE DEFAULT
               WRITE ( 6, 200 ) 
               
         END SELECT

         !==============================================================
         ! If we have found all the fields for this time, then exit 
         ! the loop and return to the calling program.  Otherwise, 
         ! go to the next iteration.
         !==============================================================
         IF ( XYMD == NYMD .and. XHMS == NHMS ) THEN
            IF ( NFOUND == N_I6_FIELDS ) THEN
               STAMP = TIMESTAMP_STRING( NYMD, NHMS )
               WRITE( 6, 210 ) NFOUND, STAMP, IU_I6
               EXIT
            ENDIF
         ENDIF
      ENDDO

      ! FORMATs
 200  FORMAT( 'Searching for next MERRA I6 field!'                    )
 210  FORMAT( '     - Found all ', i3, ' MERRA I6 met fields for ', 
     &             a, ' on unit ', i4                                 )

      !=================================================================
      ! Cleanup and quit
      !=================================================================

      ! Increment the # of times READ_I6 is called. 
      CALL SET_CT_I6( INCREMENT=.TRUE. )

      END SUBROUTINE READ_I6
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: i6_check
!
! !DESCRIPTION: Subroutine I6\_CHECK prints an error message if not all of 
!  the I6 met fields are found.  The run is also terminated.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE I6_CHECK( NFOUND, N_I6 )
!
! !USES:
!
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: NFOUND   ! Number of met fields read in from disk
      INTEGER, INTENT(IN) :: N_I6     ! Number of expected met fields
!
! !REVISION HISTORY: 
!  19 Aug 2010 - R. Yantosca - Initial version, based on i6_read_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Test if NFOUND == N_I6
      IF ( NFOUND /= N_I6 ) THEN

         ! Write error msg
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100   ) 
         WRITE( 6, 110   ) N_I6, NFOUND
         WRITE( 6, 120   )
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! FORMATs
 100     FORMAT( 'ERROR -- not enough MERRA I6 fields found!' )
 110     FORMAT( 'There are ', i2, ' fields but only ', i2 ,
     &           ' were found!'                               )
 120     FORMAT( '### STOP in I6_CHECK (merra_i6_mod.f)'      )

         ! Deallocate arrays and stop
         CALL GEOS_CHEM_STOP
      ENDIF

      END SUBROUTINE I6_CHECK
!EOC
      END MODULE MERRA_I6_MOD
