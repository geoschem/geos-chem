!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tropopause_mod
!
! !DESCRIPTION: Module TROPOPAUSE\_MOD contains routines and variables for 
!  reading and returning the value of the annual mean tropopause.
!\\
!\\
! !INTERFACE: 
!
      MODULE TROPOPAUSE_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_TROPOPAUSE
      PUBLIC :: CHECK_VAR_TROP
      PUBLIC :: COPY_FULL_TROP
      PUBLIC :: DIAG_TROPOPAUSE
      PUBLIC :: GET_MIN_TPAUSE_LEVEL
      PUBLIC :: GET_MAX_TPAUSE_LEVEL
      PUBLIC :: GET_TPAUSE_LEVEL
      PUBLIC :: INIT_TROPOPAUSE
      PUBLIC :: ITS_IN_THE_TROP
      PUBLIC :: ITS_IN_THE_STRAT
      PUBLIC :: READ_TROPOPAUSE
      PUBLIC :: SAVE_FULL_TROP
!
! !PUBLIC DATA MEMBERS:
!
      PUBLIC :: LMIN
      PUBLIC :: LMAX
!
! !REVISION HISTORY:
!  22 Aug 2005 - R. Yantosca - Initial version
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Simplify counting of tropospheric boxes (bmy, 11/1/05)
!  (3 ) Added case of variable tropopause.
!        The definition of the tropopause boxes is different in the two cases.
!        They are part of the troposphere in the case of a variable 
!        troposphere. LMAX, LMIN are the min and max extent of the troposphere
!        in that case.  (bdf, phs, 1/19/07)
!  (4 ) Bug fix: set NCS=NCSURBAN for safety's sake (bmy, 4/25/07)
!  (5 ) Updated comments (bmy, 9/18/07)
!  (6 ) Bug fix: make ITS_IN_THE_STRAT more robust. (phs, 11/14/08)
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!  29 Mar 2013 - R. Yantosca - Now make INIT_TROPOPAUSE a public function
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!      
      ! Scalars
      INTEGER              :: LMIN, LMAX

      ! Arrays
      INTEGER, ALLOCATABLE :: TROPOPAUSE(:,:)
      INTEGER, ALLOCATABLE :: IFLX(:,:)

      ! Shadow variables for values in Input_Opt
      LOGICAL              :: LVARTROP

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: copy_full_trop
!
! !DESCRIPTION: Subroutine COPY\_FULL\_TROP takes the saved full troposphere 
!  and copies chemical species into the current troposphere that will be used 
!  in SMVGEAR for this timestep.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE COPY_FULL_TROP( State_Chm )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_MOD,         ONLY : CSPEC
      USE COMODE_MOD,         ONLY : IXSAVE, IYSAVE, IZSAVE
      USE COMODE_LOOP_MOD
      USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
!
! !REMARKS:
!  ROUTINE NEEDED BECAUSE WITH VARIABLE TROPOPAUSE 
!  JLOOP WILL NOT ALWAYS REFER TO THE SAME (I,J,L) BOX
! 
! !REVISION HISTORY: 
!  14 Sep 2006 - P. Le Sager - Initial version
!  (1 ) Very similar to a get_properties of an object. Should probably
!        be in COMODE_MOD.F, and called GET_SPECIES_CONCENTRATION (phs)
!  (2 ) Bug fix: set NCS=NCSURBAN for safety's sake (bmy, 4/25/07)
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: JGAS, JLOOP, IX, IY, IZ
      INTEGER :: LOCATION(4)

      !=================================================================
      ! COPY_FULL_TROP begins here!
      !=================================================================

      ! Reset NCS to NCSURBAN for safety's sake (bmy, 4/25/07)
      NCS = NCSURBAN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( JGAS, JLOOP, IX, IY, IZ )

      ! Loop over species (active + inactive)
      DO JGAS  = 1, NTSPEC(NCS)

         ! Loop over 1-D grid boxes
         DO JLOOP = 1, NTLOOP

            ! 3-D array indices
            IX = IXSAVE(JLOOP)
            IY = IYSAVE(JLOOP)
            IZ = IZSAVE(JLOOP)

            ! Copy from 3-D array
            CSPEC(JLOOP,JGAS) = State_Chm%Species(IX,IY,IZ,JGAS)

         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE COPY_FULL_TROP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: save_full_trop
!
! !DESCRIPTION: Subroutine SAVE\_FULL\_TROP takes the current troposphere and 
!  copies chemical species into the full troposphere that will be used in 
!  SMVGEAR for this timestep.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SAVE_FULL_TROP( State_Chm )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_MOD,         ONLY : CSPEC
      USE COMODE_MOD,         ONLY : IXSAVE, IYSAVE, IZSAVE
      USE COMODE_LOOP_MOD
      USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REMARKS:
!  ROUTINE NEEDED BECAUSE WITH VARIABLE TROPOPAUSE 
!  JLOOP WILL NOT ALWAYS REFER TO THE SAME (I,J,L) BOX
! 
! !REVISION HISTORY: 
!  14 Sep 2006 - P. Le Sager - Initial version
!  (1 ) Very similar to a set_properties of an object. Should probably
!        be in COMODE_MOD.F, and called SAVE_SPECIES_CONCENTRATION (phs)
!  (2 ) Bug fix: set NCS=NCSURBAN for safety's sake! (bmy, 4/25/07)
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: JGAS, JLOOP, IX, IY, IZ

      !=================================================================
      ! SAVE_FULL_TROP begins here!
      !=================================================================

      ! Reset NCS to NCSURBAN for safety's sake (bmy, 4/25/07)
      NCS = NCSURBAN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( JGAS, JLOOP, IX, IY, IZ )

      ! Loop over species (active + inactive)
      DO JGAS = 1, NTSPEC(NCS)

         ! Loop over 1-D grid boxes
         DO JLOOP = 1, NTLOOP

            ! 3-D array indices
            IX = IXSAVE(JLOOP)
            IY = IYSAVE(JLOOP)
            IZ = IZSAVE(JLOOP)

            ! Save in 3-D array
            State_Chm%Species(IX,IY,IZ,JGAS) = CSPEC(JLOOP,JGAS)

         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE SAVE_FULL_TROP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_var_trop
!
! !DESCRIPTION: Subroutine CHECK\_VAR\_TROP checks that the entire variable 
!  troposphere is included in the 1..LLTROP range, and set the LMIN and LMAX
!  to current min and max tropopause. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_VAR_TROP( State_Met )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE GIGC_State_Met_Mod, ONLY : MetState

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_MOD                  ! LPAUSE, for backwards compatibility
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
! 
! !REVISION HISTORY: 
!  24 Aug 2006 - P. Le Sager - Initial version
!  (1 ) LLTROP is set at the first level entirely above 20 km (phs, 9/29/06)
!  (2 ) Fix LPAUSE for CH4 chemistry (phs, 1/19/07)
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J
      REAL*8  :: TPAUSE_LEV(IIPAR,JJPAR)

      !=================================================================
      ! CHECK_VAR_TROP begins here!
      !=================================================================

      ! set LMIN and LMAX to current min and max tropopause
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         TPAUSE_LEV(I,J) = GET_TPAUSE_LEVEL( I, J, State_Met )
      ENDDO
      ENDDO

      LMIN = MINVAL( TPAUSE_LEV )
      LMAX = MAXVAL( TPAUSE_LEV )

      !### For backwards compatibility during transition (still needed??)
      !### LPAUSE is still used by CH4 chemistry and ND27 (phs, 1/19/07)
      LPAUSE = TPAUSE_LEV - 1

      ! check to be sure LLTROP is large enough.
      IF ( LLTROP < LMAX ) THEN
         WRITE( 6, '(a)' ) 'CHECK_VAR_TROP: LLTROP is set too low!' 
         WRITE( 6, 10   ) LMAX, LLTROP
 10      FORMAT( 'MAX TROPOSPHERE LEVEL = ', i3, ' and LLTROP = ', i3 )
         WRITE( 6, '(a)' ) 'STOP in TROPOPAUSE_MOD.F!!!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      END SUBROUTINE CHECK_VAR_TROP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_tropopause
!
! !DESCRIPTION: Subroutine READ\_TROPOPAUSE reads in the annual mean 
!  tropopause. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_TROPOPAUSE
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD      ! Size parameters
      USE CMN_MOD           ! LPAUSE, for backwards compatibility
! 
! !REVISION HISTORY: 
!  13 Dec 1999 - Q. Li, R. Yantosca - Initial version
!  (1 ) Call READ_BPCH2 to read in the annual mean tropopause data
!        which is stored in binary punch file format. (bmy, 12/13/99)
!  (2 ) Now also read integer flags for ND27 diagnostic -- these determine
!        how to sum fluxes from boxes adjacent to the annual mean tropoause.
!        (qli, bmy, 1/7/00)
!  (3 ) Cosmetic changes (bmy, 3/17/00)
!  (4 ) Reference F90 module "bpch2_mod" which contains routine "read_bpch2"
!        for reading data from binary punch files (bmy, 6/28/00)
!  (5 ) Call TRANSFER_2D from "transfer_mod.f" to cast data from REAL*4 to
!        INTEGER and also to resize to (IIPAR,JJPAR).  ARRAY needs to be of 
!        size (IIPAR,JJPAR).  Also updated comments and made cosmetic changes. 
!        Removed obsolete variables.(bmy, 9/26/01)
!  (6 ) Removed obsolete code from 9/01 (bmy, 10/26/01)
!  (7 ) Now read annual mean tropopause files from the ann_mean_trop_200202/
!        subdirectory of DATA_DIR (bmy, 1/24/02)
!  (8 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (9 ) Now write file name to stdout (bmy, 4/3/02)
!  (10) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  (11) Now call READ_BPCH2 with QUIET=.TRUE. to suppress printing of extra
!        info to stdout.  Also updated FORMAT strings. (bmy, 3/14/03)
!  (12) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (13) Now bundled into "tropopause_mod.f' (bmy, 2/10/05)
!  (14) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (15) Simplify counting of # of tropospheric boxes (bmy, 11/1/05)
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!  29 Mar 2013 - R. Yantosca - Remove call to INIT_TROPOPAUSE, this is now
!                              done in the init stage
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, COUNT
      REAL*4             :: ARRAY(IIPAR,JJPAR,1)
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! READ_TROPOPAUSE begins here!
      !
      ! Read the annual mean tropopause from disk (binary punch file 
      ! format).  Transfer data into an array of size (IIPAR,JJPAR).
      !=================================================================
      
      ! Create filename
      FILENAME = TRIM( DATA_DIR )                      // 
     &           'ann_mean_trop_200202/ann_mean_trop.' //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT()

      ! Write file name to stdout
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - READ_TROPOPAUSE: Reading ', a )

      ! Annual mean tropopause is tracer #1  
      CALL READ_BPCH2( FILENAME, 'TR-PAUSE', 1, 
     &                 0d0,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Copy from REAL*4 to INTEGER and resize to (IIPAR,JJPAR)
      CALL TRANSFER_2D( ARRAY(:,:,1), TROPOPAUSE )

      !### For backwards compatibility during transition
      LPAUSE = TROPOPAUSE

      !=================================================================
      ! L <  TROPOPAUSE(I,J) are tropospheric boxes  
      ! L >= TROPOPAUSE(I,J) are stratospheric boxes
      !
      ! LMIN   = level where minimum extent of the TROPOPAUSE occurs
      ! LMAX   = level where maximum extent of the TROPOPAUSE occurs
      !
      ! LMIN-1 = level where minimum extent of the TROPOSPHERE occurs
      ! LMAX-1 = level where maximum extent of the TROPOSPHERE occurs
      !
      ! Write LMAX-1 and LMIN-1 to the standard output.
      !
      ! Also make sure that LMAX-1 does not exceed LLTROP, since LLTROP 
      ! is used to dimension the chemistry arrays in "comode.h". 
      !=================================================================
      LMIN = MINVAL( TROPOPAUSE )
      LMAX = MAXVAL( TROPOPAUSE )
      
      WRITE( 6, 120 ) LMIN-1
 120  FORMAT( '     - READ_TROPOPAUSE: Minimum tropospheric extent,',
     &        ' L=1 to L=', i3 )

      WRITE( 6, 130 ) LMAX-1
 130  FORMAT( '     - READ_TROPOPAUSE: Maximum tropospheric extent,',
     &        ' L=1 to L=', i3 )
    
      IF ( LMAX-1 > LLTROP ) THEN
         WRITE( 6, '(a)' ) 'READ_TROPOPAUSE: LLTROP is set too low!' 
         WRITE( 6, 131   ) LMAX-1, LLTROP
 131     FORMAT( 'LMAX = ', i3, '  LLTROP = ', i3 )
         WRITE( 6, '(a)' ) 'STOP in READ_TROPOPAUSE.F!!!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! Write the number of tropopsheric and stratospheric boxes.
      ! Recall that tropospheric boxes extend up to TROPOPAUSE - 1.
      !=================================================================
      COUNT = SUM( TROPOPAUSE - 1 )

      WRITE( 6, 140 ) COUNT
 140  FORMAT( '     - READ_TROPOPAUSE: # of tropopsheric boxes:  ', i8 )
      
      WRITE( 6, 150 ) ( IIPAR * JJPAR * LLPAR ) - COUNT
 150  FORMAT( '     - READ_TROPOPAUSE: # of stratospheric boxes: ', i8 )

      END SUBROUTINE READ_TROPOPAUSE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_max_tpause_level
!
! !DESCRIPTION: Function GET\_MAX\_TPAUSE\_LEVEL returns GEOS-Chem level at 
!  the highest extent of the annual mean tropopause.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_MAX_TPAUSE_LEVEL() RESULT( L_MAX )
!
! !RETURN VALUE:
!
      INTEGER :: L_MAX    ! Maximum tropopause level
!
! !REVISION HISTORY: 
!  10 Feb 2005 - R. Yantosca - Initial version
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_MAX_TPAUSE_LEVEL begins here!
      !=================================================================
      L_MAX = LMAX

      END FUNCTION GET_MAX_TPAUSE_LEVEL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_min_tpause_level
!
! !DESCRIPTION: Function GET\_MIN\_TPAUSE\_LEVEL returns GEOS-Chem level 
!  at the lowest extent of the annual mean tropopause.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_MIN_TPAUSE_LEVEL() RESULT( L_MIN )
!
! !RETURN VALUE:
!
      INTEGER :: L_MIN   ! Minimum tropopause level
! 
! !REVISION HISTORY: 
!  10 Feb 2005 - R. Yantosca - Initial version
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_MIN_TPAUSE_LEVEL begins here!
      !=================================================================
      L_MIN = LMIN

      END FUNCTION GET_MIN_TPAUSE_LEVEL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_tpause_level
!
! !DESCRIPTION: Function GET\_TPAUSE\_LEVEL returns the tropopause level L\_TP
!  at surface location (I,J).  Therefore, grid box (I,J,L\_TP) is partially
!  in the troposphere and partially in the stratosphere.  The grid box below
!  this, (I,J,L\_TP-1), is the last totally tropospheric box in the column.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_TPAUSE_LEVEL( I, J, State_Met ) RESULT( L_TP )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE LOGICAL_MOD,        ONLY : LVARTROP
      USE PRESSURE_MOD,       ONLY : GET_PEDGE

      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: I           ! Longitude index
      INTEGER,        INTENT(IN) :: J           ! Latitude index
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
      INTEGER                    :: L_TP   ! Tropopause level at (I,J)
!
! !REVISION HISTORY: 
!  22 Aug 2005 - R. Yantosca - Initial version
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!  10 Sep 2010 - R. Yantosca - Update comments, remove obsolete documentation
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: L
      REAL*8  :: PRESS_BEDGE

      !=================================================================
      ! GET_TPAUSE_LEVEL begins here!
      !=================================================================
      IF ( LVARTROP ) THEN

         !--------------------------
         ! Dynamic tropopause
         !--------------------------

         ! Start at the surface level
         L = 1

         ! Loop over vertical levels in the (I,J) column
         DO

            ! Pressure [hPa] at the bottom edge of grid box (I,J,L) 
            PRESS_BEDGE = GET_PEDGE( I, J, L )

            ! Break out of this loop if we encounter the box (I,J,L_TP)
            ! where the tropopause occurs.  This box is partially in the
            ! trop and partially in the strat.
            IF ( State_Met%TROPP(I,J) >= PRESS_BEDGE ) THEN
               L_TP = L - 1       
               EXIT
            ENDIF

            ! Increment L for next iteration
            L = L + 1

            ! Stop w/ error if tropopause not found
            ! (i.e. in case TROPP value is bad)
            IF ( L .GT. LLPAR ) THEN
               WRITE( 6, '(a)' ) 'GET_TPAUSE_LEVEL: CANNOT ' //
     &              'FIND T-PAUSE !'
               WRITE( 6, 160   ) L
 160           FORMAT( 'L reaches ', i3 )
               WRITE( 6, '(a)' ) 'STOP in GET_TPAUSE_LEVEL'
               WRITE( 6, '(a)' ) REPEAT( '=', 79 )
               CALL GEOS_CHEM_STOP
            ENDIF

         ENDDO

      ELSE

         !--------------------------
         ! Annual mean tropopause 
         !--------------------------

         ! Otherwise, if we are using the annual mean tropopause,
         ! set L_TP to the value read in from disk.
         L_TP = TROPOPAUSE(I,J)

      ENDIF

      END FUNCTION GET_TPAUSE_LEVEL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_in_the_trop
!
! !DESCRIPTION: Function ITS\_IN\_THE\_TROP returns TRUE if grid box (I,J,L) 
!  lies within the troposphere, or FALSE otherwise. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_IN_THE_TROP( I, J, L, State_Met ) RESULT ( IS_TROP )
!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE LOGICAL_MOD,  ONLY : LVARTROP
      USE PRESSURE_MOD, ONLY : GET_PEDGE
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: I           ! Longitude index
      INTEGER,        INTENT(IN) :: J           ! Latitude index
      INTEGER,        INTENT(IN) :: L           ! Level index
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
      LOGICAL                    :: IS_TROP   ! =T if we are in the troposphere 
!
! !REMARKS:
! 
! 
! !REVISION HISTORY: 
!  10 Feb 2005 - P. Le Sager - Initial version
!  (1 ) Modified for variable tropopause (phs, 9/14/06)
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8 :: PRESS_BEDGE

      !=================================================================
      ! ITS_IN_THE_TROP begins here
      !=================================================================
      IF ( LVARTROP ) THEN

         ! Get bottom pressure edge
         PRESS_BEDGE = GET_PEDGE(I,J,L)

         ! Check against actual tropopause pressure
         IS_TROP     = ( PRESS_BEDGE > State_Met%TROPP(I,J) )

      ELSE
         
         ! Check against annual mean tropopause
         IS_TROP     = ( L < TROPOPAUSE(I,J) ) 

      ENDIF

      END FUNCTION ITS_IN_THE_TROP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_in_the_strat
!
! !DESCRIPTION: Function ITS\_IN\_THE\_STRAT returns TRUE if grid box (I,J,L) 
!  lies within the stratosphere, or FALSE otherwise. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_IN_THE_STRAT( I, J, L, State_Met ) RESULT( IS_STRAT )
!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: I           ! Longitude index
      INTEGER,        INTENT(IN) :: J           ! Latitude index
      INTEGER,        INTENT(IN) :: L           ! Level index
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
      LOGICAL                    :: IS_STRAT  ! =T if we are in the stratosphere
!
! !REVISION HISTORY: 
!  10 Feb 2005 - P. Le Sager - Initial version
!  (1 ) Modified for variable tropopause (phs, 9/14/06)
!  (2 ) Now return the opposite value of ITS_IN_THE_TROP.  This should help
!        to avoid numerical issues. (phs, 11/14/08)
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! ITS_IN_THE_STRAT begins here
      !=================================================================
      ! Make the algorithm more robust by making ITS_IN_THE_STRAT be the 
      ! exact opposite of function ITS_IN_THE_TROP.  This should avoid
      ! numerical issues. (phs, 11/14/08)
      IS_STRAT = ( .not. ITS_IN_THE_TROP( I, J, L, State_Met ) )

      END FUNCTION ITS_IN_THE_STRAT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diag_tropopause
!
! !DESCRIPTION: Subroutine TROPOPAUSE archives the ND55 tropopause diagnostic.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG_TROPOPAUSE( State_Met )
!
! !USES:
!
      USE DIAG_MOD,           ONLY : AD55
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE LOGICAL_MOD,        ONLY : LVARTROP
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE PRESSURE_MOD,       ONLY : GET_PEDGE

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! Diagnostic switches
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REMARKS:
!  For GEOS-4, GEOS-5, 'MERRA', we use the tropopause pressure from the met 
!  field archive to determine if we are in the tropopause or not.  Therefore, 
!  the 3rd slot of AD55 should be archived with the tropopause pressure from 
!  the met fields.
!                                                                             .
!  For other met fields, we have to estimate the tropopause pressure from the
!  tropopause level.  Archive the pressure at the midpoint of the level in 
!  which the tropopause occurs.  NOTE: this may result in lower minimum 
!  tropopause pressure than reality. 
!
! !REVISION HISTORY:
!  30 Nov 1999 - H. Liu, R. Yantosca - Initial version
!  (1 ) Make sure the DO-loops go in the order L-J-I, wherever possible.
!  (2 ) Now archive ND55 diagnostic here rather than in DIAG1.F.  Also,
!        use an allocatable array (AD55) to archive tropopause heights.
!  (3 ) HTPAUSE is now a local variable, since it is only used here.
!  (4 ) Make LTPAUSE a local variable, since LPAUSE is used to store
!        the annual mean tropopause. (bmy, 4/17/00)
!  (5 ) Replace PW(I,J) with P(I,J).  Also updated comments. (bmy, 10/3/01)
!  (6 ) Removed obsolete code from 9/01 and 10/01 (bmy, 10/24/01)
!  (7 ) Added polar tropopause for GEOS-3 in #if defined( GEOS_3 ) block 
!        (bmy, 5/20/02) 
!  (8 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (9 ) Now use GET_PCENTER from "pressure_mod.f" to compute the pressure
!        at the midpoint of box (I,J,L).  Also deleted obsolete, commented-out
!        code. (dsa, bdf, bmy, 8/21/02)
!  (10) Now reference BXHEIGHT and T from "dao_mod.f".  Also reference routine
!        ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!  (11) Now uses routine GET_YMID from "grid_mod.f" to compute grid box 
!        latitude. (bmy, 2/3/03)
!  (12) Add proper polar tropopause level for GEOS-4 (bmy, 6/18/03)
!  (13) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (14) Get tropopause level from TROPOPAUSE_MOD.F routines (phs, 10/17/06)
!  10 Sep 2010 - R. Yantosca - Added ProTeX headers
!  10 Sep 2010 - R. Yantosca - For GEOS-4, GEOS-5, MERRA met fields, take the
!                              the tropopause pressure directly from the
!                              met fields rather than computing it here.
!  10 Sep 2010 - R. Yantosca - Remove reference to LPAUSE, it's obsolete
!  10 Sep 2010 - R. Yantosca - Reorganize #if blocks for clarity
!  10 Sep 2010 - R. Yantosca - Renamed to DIAG_TROPOPAUSE and bundled into
!                              tropopause_mod.f
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( GEOS_4 ) || defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
!
! !LOCAL VARIABLES:
! 
      ! Scalars
      INTEGER :: I, J,    L,  L_TP
      REAL*8  :: H, FRAC, Pb, Pt

      !=================================================================
      ! %%%%% GEOS-4, GEOS-5, MERRA met fields %%%%%
      !
      ! We get tropopause pressure directly from the met field archive
      ! Compute tropopause height to be consistent w/ the pressure
      !=================================================================
      IF ( ND55 > 0 ) THEN

         ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L_TP, H, Pb, Pt, FRAC )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            !---------------------------
            ! Compute quantities
            !---------------------------
    
            ! For this (I,J) column, get the level where the t'pause occurs
            L_TP = GET_TPAUSE_LEVEL( I, J, State_Met )

            ! Get height (from surface to top edge) of all boxes that lie
            ! totally w/in the troposphere.  NOTE: Grid box (I,J,L_TP-1)
            ! is the highest purely tropospheric grid box in the column.
            H    = SUM( State_Met%BXHEIGHT( I, J, 1:L_TP-1 ) )

            ! Get the pressures [hPa] at the bottom and top edges
            ! of the grid box in which the tropopause occurs
            Pb   = GET_PEDGE( I, J, L_TP   )  
            Pt   = GET_PEDGE( I, J, L_TP+1 )

            ! FRAC is the fraction of the grid box (I,J,L_TP) 
            ! that lies totally within the troposphere
            FRAC = ( Pb - State_Met%TROPP(I,J) ) / ( Pb - Pt ) 

            ! Add to H the height [m] of the purely tropospheric 
            ! fraction of grid box (I,J,L_TP)
            H    = H + ( FRAC * State_Met%BXHEIGHT(I,J,L_TP) )

            !---------------------------
            ! Archive into ND55 array
            !---------------------------
            AD55(I,J,1) = AD55(I,J,1) + L_TP        ! T'pause level
            AD55(I,J,2) = AD55(I,J,2) + H/1.0d3     ! T'pause height [km]
            AD55(I,J,3) = AD55(I,J,3) + 
     &                    State_Met%TROPP(I,J)      ! T'pause pressure [hPa]

         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

#else

!
! !LOCAL VARIABLES:
! 
      ! Scalars
      INTEGER :: I, J, L

      ! Arrays
      REAL*8  :: H(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! %%%%% ALL OTHER MET FIELDS %%%%%
      !
      ! We compute tropopause pressure from the tropopause level (which 
      ! is taken from the thermally-derived annual mean tropopause data 
      ! read from disk).
      !
      ! NOTE: Keep the existing algorithm for backwards compatibility.
      !=================================================================

      ! Find height of the midpoint of the first level
      ! H (in m) is the height of the midpoint of layer L (hyl, 03/28/99)
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         H(I,J,1) = State_Met%BXHEIGHT(I,J,1) / 2.d0
      ENDDO
      ENDDO

      ! Add to H 1/2 of the sum of the two adjacent boxheights
      DO L = 1, LLPAR-1
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         H(I,J,L+1) = H(I,J,L) + 
     &               ( State_Met%BXHEIGHT(I,J,L) +
     &                 State_Met%BXHEIGHT(I,J,L+1) ) / 2.d0
      ENDDO
      ENDDO
      ENDDO

      !=================================================================
      ! ND55: Tropopause level, height [ km ], and pressure [ mb ]
      !=================================================================
      IF ( ND55 > 0 ) THEN
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Get the tropopause level
            L           = GET_TPAUSE_LEVEL( I, J, State_Met )

            ! If we are using the variable tropopause, then (I,J,L) is the
            ! highest purely tropospheric grid box.  The grid box in which
            ! the tropopause actually occurs is then (I,J,L+1).
            IF ( LVARTROP ) L = L + 1

            ! Archive level at which tropopause occurs
            AD55(I,J,1) = AD55(I,J,1) + L

            ! Archive tropopause height [km]
            AD55(I,J,2) = AD55(I,J,2) + H(I,J,L) / 1.0d3 ! m --> km

            ! We have to estimate the tropopause pressure from the 
            ! tropopause level.  Archive the pressure at the midpoint
            ! of the level in which the tropopause occurs.  NOTE: this may
            ! result in lower minimum tropopause pressure than reality.
            AD55(I,J,3) = AD55(I,J,3) + GET_PCENTER(I,J,L)

         ENDDO
         ENDDO
      ENDIF

#endif

      END SUBROUTINE DIAG_TROPOPAUSE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tropopause
!
! !DESCRIPTION: Subroutine INIT\_TROPOPAUSE allocates and zeroes module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TROPOPAUSE( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE ERROR_MOD,          ONLY : ALLOC_ERR
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REMARKS:
!  Now 
!
! !REVISION HISTORY: 
!  10 Feb 2005 - R. Yantosca - Initial version
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!  29 Mar 2013 - R. Yantosca - Now made public so we can shadow LVARTROP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS 

      !=================================================================
      ! INIT_TROPOPAUSE
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Shadow Input_Opt%LVARTROP in a module variable
      LVARTROP = Input_Opt%LVARTROP

      ! Only allocate arrays if we are not using the variable tropopause
      IF ( .not. LVARTROP ) THEN 

         ! Array for the ann mean tropopause
         ALLOCATE( TROPOPAUSE( IIPAR, JJPAR ), STAT=AS ) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TROPOPAUSE' )
         TROPOPAUSE = 0

         ! For now don't allocate IFLX
         !ALLOCATE( IFLX( IIPAR, JJPAR ), STAT=AS ) 
         !IF ( AS /= 0 ) CALL ALLOC_ERR( 'IFLX' )
         !IFLX = 0
         
      ENDIF

      END SUBROUTINE INIT_TROPOPAUSE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_tropopause
!
! !DESCRIPTION: Subroutine CLEANUP\_TROPOPAUSE deallocates module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_TROPOPAUSE
! 
! !REVISION HISTORY: 
!  10 Feb 2005 - R. Yantosca - Initial version
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_TROPOPAUSE begins here!
      !=================================================================
      IF ( ALLOCATED( TROPOPAUSE ) ) DEALLOCATE( TROPOPAUSE )
      IF ( ALLOCATED( IFLX       ) ) DEALLOCATE( IFLX       ) 

      END SUBROUTINE CLEANUP_TROPOPAUSE 
!EOC
      END MODULE TROPOPAUSE_MOD
