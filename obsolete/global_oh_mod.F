!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_oh_mod
!
! !DESCRIPTION: Module GLOBAL\_OH\_MOD contains variables and routines for 
!  reading the global monthly mean OH concentration from disk. 
!\\
!\\
! !INTERFACE: 
!
      MODULE GLOBAL_OH_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      ! Array to store global monthly mean OH field [molec/cm3]
      REAL(fp), PUBLIC, ALLOCATABLE :: OH(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_GLOBAL_OH
      PUBLIC :: GET_GLOBAL_OH
      PUBLIC :: GET_DIURNAL_OH
      PUBLIC :: GET_SZAFACT
      PUBLIC :: INIT_GLOBAL_OH
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: CALC_SUMCOSZA
!
! !REVISION HISTORY:
!  28 Jul 2000 - R. Yantosca - Initial version
!  (1 ) Updated comments (bmy, 9/4/01)
!  (2 ) Now use routines from "transfer_mod.f" to regrid OH to 30 levels
!        for reduced GEOS-3 grid.  Also size OH array properly. (bmy, 1/14/02)
!  (3 ) Eliminate obsolete code from 11/01 (bmy, 2/27/02)
!  (4 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (5 ) Now use updated OH fields (bmy, 10/2/02)
!  (6 ) Now references "error_mod.f" (bmy, 10/15/02)
!  (7 ) Minor bug fixes in FORMAT statements (bmy, 3/23/03)
!  (8 ) Cosmetic changes to simplify output (bmy, 3/27/03)
!  (9 ) Bug fix: OH should be (IIPAR,JJPAR,LLPAR) (bmy, 5/4/04)
!  (10) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  14 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Sum of cosine of the solar zenith angle. Used to impose a
      ! diurnal variability on OH concentrations
      REAL(fp), ALLOCATABLE         :: SUMCOSZA(:,:)

      ! Allocation status
      LOGICAL                     :: Arrays_Allocated = .FALSE.

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_global_oh
!
! !DESCRIPTION: Subroutine GET\_GLOBAL\_OH reads global OH from binary 
!  punch files stored in the /data/ctm/GEOS\_MEAN directory.  This OH data 
!  is needed as oxidant for various offline chemistry mechanisms.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GLOBAL_OH( am_I_Root, Input_Opt, 
     &                          State_Met, THISMONTH, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT
      USE BPCH2_MOD,            ONLY : GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0
      USE BPCH2_MOD,            ONLY : READ_BPCH2
      USE DAO_MOD,              ONLY : AIRDEN_FULLGRID
      USE DAO_MOD,              ONLY : AIRQNT_FULLGRID
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod,   ONLY : OptInput
      USE GIGC_State_Met_Mod,   ONLY : MetState
      USE TRANSFER_MOD,         ONLY : TRANSFER_3D
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Obj for input options
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Obj for met fields
      INTEGER,        INTENT(IN)    :: THISMONTH   ! Current month
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  28 Jul 2000 - R. Yantosca - Initial version
!  (1 ) GET_GLOBAL_OH assumes that we are reading global OH data that occupies
!        all CTM levels.  Contact Bob Yantosca (bmy@io.harvard.edu) for IDL
!        regridding code which will produce the appropriate OH files.
!  (2 ) Now use version of GET_TAU0 with 3 arguments.  Now call READ_BPCH2
!        with IIPAR,JJPAR,LGLOB.  Call TRANSFER_3D to cast from REAL*4 to
!        REAL(fp) and to regrid to 30 levels for GEOS-3 (if necessary).
!        ARRAY should now be of size (IIPAR,JJPAR,LGLOB). (bmy, 1/11/02)
!  (3 ) Now point to new OH files in the v4-26 subdirectory.  Also eliminated
!        obsolete code from 11/01. (bmy, 2/27/02)
!  (4 ) Now point to OH files in the v4-33 subdirectory. (bmy, 10/2/02)
!  (5 ) Replace missing commas in the FORMAT statement (bmy, 3/23/03)
!  (6 ) Cosmetic changes to simplify output (bmy, 3/27/03)
!  (7 ) Add Mat's OH as an option.  Also read bpch file quietly (bmy, 5/4/04)
!  (8 ) Now use OH_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!  23 Oct 2013 - R. Yantosca - Now accept Input_Opt argument
!  29 Oct 2013 - R. Yantosca - Remove TRANSFER_3D_NOLUMP routine, we can just
!                              instead do a direct cast assignment
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, L
      REAL(f4)           :: ARRAY(IIPAR,JJPAR,LGLOB)
      REAL(f8)           :: XTAU
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! GET_GLOBAL_OH begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Allocate OH array, if this is the first call
      IF ( .NOT. Arrays_Allocated ) THEN
         CALL INIT_GLOBAL_OH
      ENDIF

#if defined( NESTED_NA )
      ! Filename
      FILENAME = TRIM( Input_Opt%OH_DIR ) // 'OH_3Dglobal.' // 
     &           GET_NAME_EXT()           // '.'            //
     &           GET_RES_EXT()  // '_NA'
#else
      ! Filename
      FILENAME = TRIM( Input_Opt%OH_DIR ) // 'OH_3Dglobal.' // 
     &           GET_NAME_EXT()           // '.'            // 
     &           GET_RES_EXT()
#endif

      ! Echo some information to the standard output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_OH: Reading OH from: ', a )

      ! Get the TAU0 value for the start of the given month
      ! Assume "generic" year 1985 (TAU0 = [0, 744, ... 8016])
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Read OH data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 1,     
     &                 XTAU,      IIPAR,     JJPAR,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! Do not use routines AIRQNT_FULLGRID or TRANSFER_3D_NO_LUMP when
      ! we are connecting to the GEOS-5 GCM via ESMF. (mpayer, 3/14/13)
      !-----------------------------------------------------------------

      ! Assign data from ARRAY to the module variable OH
      CALL TRANSFER_3D( ARRAY, OH )

#else
      !-----------------------------------------------------------------
      !                %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Include fix for vertical regridding of OH for offline simulations
      ! provided by Chris Holmes (mpayer, 3/14/13)
      !-----------------------------------------------------------------
#if defined( GRIDREDUCED )

      ! Calculate the number density of air on the full vertical grid
      CALL AIRQNT_FULLGRID( am_I_Root, Input_Opt, RC )

      ! Convert molec/cm3 -> mixing ratio (molec m3 kg-1 cm-3)
      ! because regridding requires mixing ratio quantity
      ! Note that peculiar units don't matter because we immediately
      ! convert back after regridding
      DO L=1, LGLOB
         ARRAY(:,:,L) = ARRAY(:,:,L) / AIRDEN_FULLGRID(:,:,L)
      ENDDO
         
      ! Assign data from ARRAY to the module variable OH
      CALL TRANSFER_3D( ARRAY, OH )
      
      ! convert mixing ratio (molec m3 kg-1 cm-3) -> molec/cm3
      DO L=1, LLPAR
         OH(:,:,L) = OH(:,:,L) * State_Met%AIRDEN(:,:,L)
      ENDDO

#else

      ! Assign data from ARRAY to the module variable OH
      OH = ARRAY

#endif
#endif

      END SUBROUTINE GET_GLOBAL_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_szafact
!
! !DESCRIPTION:
!  Subroutine GET\_SZAFACT returns diurnal scale factors from dividing
!  the sza by the sum of the total sza per day. These factors are mainly
!  imposed to the monthly OH climatology. 
!  However, the same scale factors are dimensionless and can hence be 
!  applied to other compounds too (e.g. O3).
!\\
! !INTERFACE:
!
      FUNCTION GET_SZAFACT ( I, J, State_Met ) RESULT( FACT )
!
! !USES:
!
      ! References to F90 modules
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !ARGUMENTS:
!
      INTEGER,        INTENT(IN) :: I, J
      TYPE(MetState), INTENT(IN) :: State_Met

      ! Local variables
      REAL(fp)              :: FACT
!
! !REVISION HISTORY: 
! 01 Mar 2013 - C. Keller - Imported from carbon_mod.F, where these
! calculations are done w/in GET_OH
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! GET_SZAFACT begins here!
      !=================================================================

      ! Check for allocation status of arrays 
      IF ( .NOT. Arrays_Allocated ) CALL INIT_GLOBAL_OH

      ! Update SUMCOSZA
      CALL CALC_SUMCOSZA

      ! Test for sunlight...
      IF ( State_Met%SUNCOSmid(I,J) > 0e+0_fp  .AND. 
     &     SUMCOSZA(I,J)            > 0e+0_fp ) THEN

         ! Impose a diurnal variation on OH during the day
         FACT = ( State_Met%SUNCOSmid(I,J) / SUMCOSZA(I,J) ) *
     &          ( 1440e+0_fp                   / GET_TS_CHEM() )

      ELSE

         ! At night, OH goes to zero
         FACT = 0e+0_fp

      ENDIF

      END FUNCTION GET_SZAFACT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_diurnal_oh
!
! !DESCRIPTION:
!  Subroutine GET\_DIURNAL\_OH returns the OH concentration of a given
!  field with the diurnal scaling factor applied to it. 
!\\
! !INTERFACE:
!
      FUNCTION GET_DIURNAL_OH ( I, J, L, State_Met ) 
     &  RESULT( OH_MOLEC_CM3 )
!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !ARGUMENTS:
!
      INTEGER,        INTENT(IN) :: I, J, L
      TYPE(MetState), INTENT(IN) :: State_Met
!
! !REVISION HISTORY: 
! 01 Mar 2013 - C. Keller - Imported from carbon_mod.F, where this
! calculations are done w/in GET_OH
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local variables
      REAL(fp)              :: FACT
      REAL(fp)              :: OH_MOLEC_CM3

      !=================================================================
      ! GET_DIURNAL_OH begins here!
      !=================================================================

      ! Calculate fraction to be applied
      FACT = GET_SZAFACT( I, J, State_Met )

      ! Multiply monthly average OH field with diurnal scale factor
      OH_MOLEC_CM3 = OH(I,J,L) * FACT 

      ! Make sure OH is not negative
      OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0e+0_fp )

      END FUNCTION GET_DIURNAL_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_sumcosza
!
! !DESCRIPTION:
!  Subroutine CALC\_SUMCOSZA computes the sum of cosine of the solar zenith
!  angle over a 24 hour day, as well as the total length of daylight. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_SUMCOSZA
!
! !USES:
!
      USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC
      USE TIME_MOD, ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_GCTM_MOD
!
! !REVISION HISTORY: 
! 01 Mar 2013 - C. Keller - Imported from carbon_mod.F, where it's
! called OHNO3TIME
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, SAVE       :: SAVEDOY = -1
      INTEGER             :: I, IJLOOP, J, L, N, NT, NDYSTEP
      REAL(fp)              :: A0, A1, A2, A3, B1, B2, B3
      REAL(fp)              :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
      REAL(fp)              :: SUNTMP(MAXIJ)

      !=================================================================
      ! CALC_SUMCOSZA begins here!
      !=================================================================

      !  Solar declination angle (low precision formula, good enough for us):
      A0 = 0.006918
      A1 = 0.399912
      A2 = 0.006758
      A3 = 0.002697
      B1 = 0.070257
      B2 = 0.000907
      B3 = 0.000148
      R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

      DEC = A0 - A1*cos(  R) + B1*sin(  R)
     &         - A2*cos(2*R) + B2*sin(2*R)
     &         - A3*cos(3*R) + B3*sin(3*R)

      LHR0 = int(float( GET_NHMSb() )/10000.)

      ! Only do the following at the start of a new day
      IF ( SAVEDOY /= GET_DAY_OF_YEAR() ) THEN 

         ! Zero arrays
         SUMCOSZA(:,:) = 0e+0_fp

         ! NDYSTEP is # of chemistry time steps in this day
         NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 60 / GET_TS_CHEM()      

         ! NT is the elapsed time [s] since the beginning of the run
         NT = GET_ELAPSED_SEC()

         ! Loop forward through NDYSTEP "fake" timesteps for this day 
         DO N = 1, NDYSTEP
            
            ! Zero SUNTMP array
            SUNTMP(:) = 0e+0_fp

            ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YMID_R, IJLOOP, TIMLOC, AHR )
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Grid box latitude center [radians]
               YMID_R = GET_YMID_R( I, J, 1 )

               ! Increment IJLOOP
               IJLOOP = ( (J-1) * IIPAR ) + I
               TIMLOC = real(LHR0) + real(NT)/3600.0 + 
     &                  GET_XMID( I, J, 1 ) / 15.0
         
               DO WHILE (TIMLOC .lt. 0)
                  TIMLOC = TIMLOC + 24.0
               ENDDO

               DO WHILE (TIMLOC .gt. 24.0)
                  TIMLOC = TIMLOC - 24.0
               ENDDO

               AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

            !===========================================================
            ! The cosine of the solar zenith angle (SZA) is given by:
            !     
            !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR) 
            !                   
            ! where LAT = the latitude angle, 
            !       DEC = the solar declination angle,  
            !       AHR = the hour angle, all in radians. 
            !
            ! If SUNCOS < 0, then the sun is below the horizon, and 
            ! therefore does not contribute to any solar heating.  
            !===========================================================

               ! Compute Cos(SZA)
               SUNTMP(IJLOOP) = sin(YMID_R) * sin(DEC) +
     &                          cos(YMID_R) * cos(DEC) * cos(AHR)

               ! SUMCOSZA is the sum of SUNTMP at location (I,J)
               ! Do not include negative values of SUNTMP
               SUMCOSZA(I,J) = SUMCOSZA(I,J) + 
     &                         MAX(SUNTMP(IJLOOP),0e+0_fp)

            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Increment elapsed time [sec]
            NT = NT + ( GET_TS_CHEM() * 60 )             
         ENDDO

         ! Set saved day of year to current day of year 
         SAVEDOY = GET_DAY_OF_YEAR()

      ENDIF

      ! Return to calling program
      END SUBROUTINE CALC_SUMCOSZA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_oh
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_OH allocates and zeroes
!  all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GLOBAL_OH
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD 
! 
! !REVISION HISTORY: 
!  28 Jul 2000 - R. Yantosca - Initial version
!  (1 ) OH array now needs to be sized (IIPAR,JJPAR,LGLOB) (bmy, 1/14/02)
!  (2 ) Also eliminated obsolete code from 11/01 (bmy, 2/27/02)
!  (3 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (4 ) OH should be (IIPAR,JJPAR,LLPAR): avoid subscript errors (bmy, 5/4/04)
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_GLOBAL_OH begins here!
      !=================================================================

      ! Allocate OH array
      ALLOCATE( OH( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OH' )

      ! Zero OH array
      OH = 0e+0_fp

      ! Allocate SUMCOSZA array
      ALLOCATE( SUMCOSZA( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUMCOSZA' )
      SUMCOSZA = 0e+0_fp

      ! Arrays are now allocated
      Arrays_Allocated = .TRUE.

      END SUBROUTINE INIT_GLOBAL_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_global_oh
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_OH deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GLOBAL_OH
! 
! !REVISION HISTORY: 
!  28 Jul 2000 - R. Yantosca - Initial version
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_GLOBAL_OH begins here!
      !=================================================================
      IF ( ALLOCATED( OH       ) ) DEALLOCATE( OH       ) 
      IF ( ALLOCATED( SUMCOSZA ) ) DEALLOCATE( SUMCOSZA ) 
     
      END SUBROUTINE CLEANUP_GLOBAL_OH
!EOC
      END MODULE GLOBAL_OH_MOD
