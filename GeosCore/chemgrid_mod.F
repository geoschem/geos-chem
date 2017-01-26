!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: chemgrid_mod.F
!
! !DESCRIPTION: Module CHEMGRID\_MOD contains routines and variables for 
!  reading and returning vertical layer limits.
!\\
!\\
! !INTERFACE: 
!
      MODULE CHEMGRID_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : LLSTRAT

      USE PRECISION_MOD   ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_CHEMGRID
      PUBLIC  :: CHECK_VAR_TROP
      PUBLIC  :: DIAG_TROPOPAUSE
      PUBLIC  :: GET_MIN_CHEMGRID_LEVEL
      PUBLIC  :: GET_MAX_CHEMGRID_LEVEL
      PUBLIC  :: GET_MIN_TPAUSE_LEVEL
      PUBLIC  :: GET_MAX_TPAUSE_LEVEL
      PUBLIC  :: GET_TPAUSE_LEVEL
      PUBLIC  :: GET_CHEMGRID_LEVEL
      PUBLIC  :: INIT_CHEMGRID
      PUBLIC  :: ITS_IN_THE_TROP
      PUBLIC  :: ITS_IN_THE_STRAT
      PUBLIC  :: ITS_IN_THE_MESO
      PUBLIC  :: ITS_IN_THE_STRATMESO
      PUBLIC  :: ITS_IN_THE_CHEMGRID
      PUBLIC  :: ITS_IN_THE_NOCHEMGRID
      PUBLIC  :: READ_TROPOPAUSE
!
! !PUBLIC DATA MEMBERS:
!
      ! Scalars
      INTEGER, PUBLIC :: LMIN   ! Minimum level where chemistry happens
      INTEGER, PUBLIC :: LMAX   ! Maximum levle where chemistry happens
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
!  09 Sep 2010 - R. Yantosca   - Added ProTeX headers
!  24 Mar 2013 - S. D. Eastham - Adapted from tropopause_mod.F. Can now
!                                identify atmospheric layers, chemistry grids
!                                and aerosol grids independently of one another
!  20 Feb 2014 - M. Sulprizio  - Removed "define.h", this is now obsolete
!  21 Feb 2014 - M. Sulprizio  - Now make INIT_CHEMGRID a public function
!  06 Nov 2014 - M. Yannetti   - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!      
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
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE State_Met_Mod,      ONLY : MetState
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
!  21 Feb 2014 - M. Sulprizio- Now pass State_Met object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J
      REAL(fp)  :: TPAUSE_LEV(IIPAR,JJPAR)

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

      ! check to be sure LLTROP is large enough.
      IF ( LLTROP < LMAX ) THEN
         WRITE( 6, '(a)' ) 'CHECK_VAR_TROP: LLTROP is set too low!' 
         WRITE( 6, 10   ) LMAX, LLTROP
 10      FORMAT( 'MAX TROPOSPHERE LEVEL = ', i3, ' and LLTROP = ', i3 )
         WRITE( 6, '(a)' ) 'STOP in CHECK_VAR_TROP.F!!!'
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
      SUBROUTINE READ_TROPOPAUSE( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_NAME_EXT, GET_RES_EXT
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE Input_Opt_Mod,      ONLY : OptInput
      USE m_netcdf_io_open
      USE m_netcdf_io_read
      USE m_netcdf_io_readattr
      USE m_netcdf_io_close
!
! !INPUT PARAMETERS:
!     
      LOGICAL,        INTENT(IN)  :: am_I_Root  ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt  ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC         ! Success or failure?
! 
! !REMARKS:
!  ##########################################################################
!  #####    NOTE: BINARY PUNCH INPUT IS BEING PHASED OUT.  THIS DATA    #####
!  #####    WILL EVENTUALLY BE READ IN FROM netCDF FILES VIA HEMCO!     #####
!  #####       -- Bob Yantosca (05 Mar 2015)                            #####
!  ##########################################################################
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
!  21 Feb 2014 - M. Sulprizio- Remove call to INIT_TROPOPAUSE, this is now
!                              done in the init stage
!  20 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  06 Nov 2014 - R. Yantosca - Replace TRANSFER_2D with direct casts
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  05 Mar 2015 - R. Yantosca - Add Input_Opt%RES_DIR to data path
!  06 Mar 2015 - R. Yantosca - Now read ann mean trop from netCDF 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, COUNT
      REAL(f4)           :: ARRAY(IIPAR,JJPAR,1)
      CHARACTER(LEN=255) :: FILENAME

      ! Variables for netCDF
      INTEGER            :: fId                 ! netCDF file ID
      INTEGER            :: st3d(3), ct3d(3)    ! For 3D arrays 
      CHARACTER(LEN=255) :: v_name              ! netCDF variable name 
      CHARACTER(LEN=255) :: a_name              ! netCDF attribute name
      CHARACTER(LEN=255) :: a_val               ! netCDF attribute value

      !=================================================================
      ! READ_TROPOPAUSE begins here!
      !
      ! Read the annual mean tropopause from disk (binary punch file 
      ! format).  Transfer data into an array of size (IIPAR,JJPAR).
      !=================================================================
      
      ! Assume success
      RC = GC_SUCCESS

      !=================================================================
      ! Open and read data from the netCDF file
      !=================================================================

      ! Create filename
      FILENAME = TRIM( Input_Opt%CHEM_INPUTS_DIR )      // 
     &           'ann_mean_trop_200202/ann_mean_trop.'  //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT() // '.nc'

      ! Open netCDF file
      CALL Ncop_Rd( fId, TRIM( FILENAME ) )
      
      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 200 ) REPEAT( '%', 79 )
         WRITE( 6, 210 ) TRIM( FILENAME )
      ENDIF
    
      !----------------------------------------
      ! VARIABLE: TP_LEVEL
      !----------------------------------------

      ! Variable and attribute names
      v_name = "TP_LEVEL"
      a_name = "units"

      ! Read TP_LEVEL from file
      st3d   = (/ 1,     1,     1 /)
      ct3d   = (/ IIPAR, JJPAR, 1 /)
      CALL NcRd( ARRAY, fId, TRIM( v_name ), st3d, ct3d )
      
      ! Read the TP_LEVEL:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes(fId, TRIM(v_name), TRIM(a_name), a_val)

      ! Echo info to stdout
      IF ( am_I_Root ) THEN 
         WRITE( 6, 220 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      ! Close netCDF file
      CALL NcCl( fId )
    
      ! Echo info to stdout
      IF ( am_I_Root ) THEN 
         WRITE( 6, 230 )
         WRITE( 6, 200 ) REPEAT( '%', 79 )
      ENDIF

      ! Copy from REAL*4 to INTEGER and resize to (IIPAR,JJPAR)
      TROPOPAUSE = ARRAY(:,:,1)

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
      
      IF ( am_I_Root ) THEN
         WRITE( 6, 120 ) LMIN-1
 120     FORMAT( '     - READ_TROPOPAUSE: Minimum tropospheric extent,',
     &           ' L=1 to L=', i3 )

         WRITE( 6, 130 ) LMAX-1
 130     FORMAT( '     - READ_TROPOPAUSE: Maximum tropospheric extent,',
     &            ' L=1 to L=', i3 )
      ENDIF

      IF ( LMAX-1 > LLTROP ) THEN
         IF ( am_I_Root ) THEN
            WRITE( 6, '(a)' ) 'READ_TROPOPAUSE: LLTROP is set too low!' 
            WRITE( 6, 131   ) LMAX-1, LLTROP
 131        FORMAT( 'LMAX = ', i3, '  LLTROP = ', i3 )
            WRITE( 6, '(a)' ) 'STOP in READ_TROPOPAUSE.F!!!'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         ENDIF
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! Write the number of tropopsheric and stratospheric boxes.
      ! Recall that tropospheric boxes extend up to TROPOPAUSE - 1.
      !=================================================================
      COUNT = SUM( TROPOPAUSE - 1 )

      IF ( am_I_Root ) THEN
       WRITE( 6, 140 ) COUNT
 140   FORMAT( '     - READ_TROPOPAUSE: # of tropopsheric boxes:  ', i8)
      
       WRITE( 6, 150 ) ( IIPAR * JJPAR * LLPAR ) - COUNT
 150   FORMAT( '     - READ_TROPOPAUSE: # of stratospheric boxes: ', i8)
      ENDIF

 200  FORMAT( a                                              )
 210  FORMAT( '%% Opening file  : ',         a               )
 220  FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
 230  FORMAT( '%% Successfully closed file!'                 )

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
! !IROUTINE: get_min_chemgrid_level
!
! !DESCRIPTION: Function GET\_MIN\_CHEMGRID\_LEVEL returns GEOS-Chem level 
!  at the lowest extent of the chemistry grid
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_MIN_CHEMGRID_LEVEL() RESULT( L_MIN )
!
! !USES:
!
!
! !RETURN VALUE:
!
      INTEGER :: L_MIN   ! Minimum chemistry grid level
! 
! !REVISION HISTORY: 
!  10 Feb 2005 - R. Yantosca   - Initial version
!  09 Sep 2010 - R. Yantosca   - Added ProTeX headers
!  25 Mar 2013 - S. D. Eastham - Adapted from GET\_MIN\_TPAUSE\_LEVEL
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! GET_MIN_CHEMGRID_LEVEL begins here!
      !=================================================================

#if defined( UCX )
      ! Minimum level is fixed at stratopause
      L_MIN = LLSTRAT
#else
      L_MIN = LMIN
#endif

      END FUNCTION GET_MIN_CHEMGRID_LEVEL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_max_chemgrid_level
!
! !DESCRIPTION: Function GET\_MAX\_CHEMGRID\_LEVEL returns GEOS-Chem level at 
!  the highest extent of the chemistry grid.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_MAX_CHEMGRID_LEVEL() RESULT( L_MAX )
!
! !USES:
!
!
! !RETURN VALUE:
!
      INTEGER :: L_MAX    ! Maximum chemistry grid level
!
! !REVISION HISTORY: 
!  10 Feb 2005 - R. Yantosca   - Initial version
!  09 Sep 2010 - R. Yantosca   - Added ProTeX headers
!  28 Aug 2013 - S. D. Eastham - Adapted from GET\_MAX\_TPAUSE\_LEVEL
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! GET_MAX_CHEMGRID_LEVEL begins here!
      !=================================================================

#if defined( UCX )
      ! Maximum level is fixed at stratopause
      L_MAX = LLSTRAT
#else
      L_MAX = LMAX
#endif

      END FUNCTION GET_MAX_CHEMGRID_LEVEL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_chemgrid_level
!
! !DESCRIPTION: Function GET\_CHEMGRID\_LEVEL returns the chemistry grid limit
!  L\_CHEM at surface location (I,J). Therefore, grid box (I,J,L\_CHEM) is
!  partially within the chemistry grid and partially outside of it. The grid
!  box below this, (I,J,L\_CHEM-1), is the last fully legitimate chemistry
!  grid box in the column.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CHEMGRID_LEVEL( I, J, State_Met ) RESULT( L_CHEM )
!
! !USES:
!
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN) :: I           ! Longitude index
      INTEGER,        INTENT(IN) :: J           ! Latitude index
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
      INTEGER             :: L_CHEM ! Chemgrid limit at (I,J)
!
! !REVISION HISTORY:
!  25 Mar 2013 - S. D. Eastham - Initial version
!  21 Feb 2014 - M. Sulprizio  - Now pass State_Met object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! GET_CHEMGRID_LEVEL begins here!
      !=================================================================

#if defined( UCX )
      ! Chemistry grid level is fixed at stratopause
      ! None of the vertical grids currently employed have terrain
      ! following coordinates at this altitude, so the level is fixed
      L_CHEM = LLSTRAT
#else
      ! Use GET_TPAUSE_LEVEL otherwise - no point in duplicating code
      L_CHEM = GET_TPAUSE_LEVEL( I, J, State_Met )
#endif
      END FUNCTION GET_CHEMGRID_LEVEL
!
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
      USE State_Met_Mod,      ONLY : MetState

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
!  21 Feb 2014 - M. Sulprizio- Now pass State_Met object via the arg list
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!                              Remove dependency on pressure_mod.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: L
      REAL(fp)  :: PRESS_BEDGE

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
            PRESS_BEDGE = State_Met%PEDGE( I, J, L )

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
               WRITE( 6, '(a)' ) 'CGM: GET_TPAUSE_LEVEL: CANNOT ' //
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
      USE State_Met_Mod,      ONLY : MetState
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
!  21 Feb 2014 - M. Sulprizio- Now pass State_Met object via the arg list
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!                              Remove dependency on pressure_mod.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: PRESS_BEDGE

      !=================================================================
      ! ITS_IN_THE_TROP begins here
      !=================================================================
      IF ( LVARTROP ) THEN

         ! Get bottom pressure edge
         PRESS_BEDGE = State_Met%PEDGE(I,J,L)

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
! !IROUTINE: its_in_the_stratmeso
!
! !DESCRIPTION: Function ITS\_IN\_THE\_STRATMESO returns TRUE if grid box 
!  (I,J,L) lies without the troposphere, or FALSE otherwise. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_IN_THE_STRATMESO( I, J, L, State_Met )
     &                             RESULT( IS_SMESO )
!
! !USES:
!
      USE State_Met_Mod,      ONLY : MetState
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
      LOGICAL                    :: IS_SMESO    ! =F if we are in the trop
!
! !REVISION HISTORY: 
!  10 Feb 2005 - P. Le Sager - Initial version
!  (1 ) Modified for variable tropopause (phs, 9/14/06)
!  (2 ) Now return the opposite value of ITS_IN_THE_TROP.  This should help
!        to avoid numerical issues. (phs, 11/14/08)
!  09 Sep 2010 - R. Yantosca   - Added ProTeX headers
!  25 Mar 2013 - S. D. Eastham - Adapted from ITS_IN_THE_STRAT
!  21 Feb 2014 - M. Sulprizio  - Now pass State_Met object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! ITS_IN_THE_STRATMESO begins here
      !=================================================================
      ! Make the algorithm more robust by making ITS_IN_THE_STRATMESO be 
      ! exact opposite of function ITS_IN_THE_TROP.  This should avoid
      ! numerical issues. (phs, 11/14/08)
      IS_SMESO = ( .not. ITS_IN_THE_TROP( I, J, L, State_Met ) )

      END FUNCTION ITS_IN_THE_STRATMESO
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
      USE State_Met_Mod,      ONLY : MetState
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
!  24 Mar 2013 - S. D. Eastham - Now distinguish between strat and mesosphere
!  21 Feb 2014 - M. Sulprizio  - Now pass State_Met object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! ITS_IN_THE_STRAT begins here
      !=================================================================

      IS_STRAT = ( ( L .le. LLSTRAT ) .and.
     &             ( ITS_IN_THE_STRATMESO( I, J, L, State_Met ) ) )

      END FUNCTION ITS_IN_THE_STRAT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_in_the_meso
!
! !DESCRIPTION: Function ITS\_IN\_THE\_MESO returns TRUE if grid box (I,J,L) 
!  lies within the mesosphere, or FALSE otherwise. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_IN_THE_MESO( I, J, L ) RESULT( IS_MESO )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I          ! Longitude index
      INTEGER, INTENT(IN) :: J          ! Latitude index
      INTEGER, INTENT(IN) :: L          ! Level index
!
! !RETURN VALUE:
!
      LOGICAL             :: IS_MESO    ! =T if we are in the mesosphere
!
! !REVISION HISTORY: 
!  10 Feb 2005 - P. Le Sager - Initial version
!  (1 ) Modified for variable tropopause (phs, 9/14/06)
!  (2 ) Now return the opposite value of ITS_IN_THE_TROP.  This should help
!        to avoid numerical issues. (phs, 11/14/08)
!  09 Sep 2010 - R. Yantosca   - Added ProTeX headers
!  24 Mar 2013 - S. D. Eastham - Now distinguish between strat and mesosphere
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! ITS_IN_THE_MESO begins here
      !=================================================================

      IS_MESO = ( L .gt. LLSTRAT )

      END FUNCTION ITS_IN_THE_MESO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_in_the_chemgrid
!
! !DESCRIPTION: Function ITS\_IN\_THE\_CHEMGRID returns TRUE if grid box
!  (I,J,L) lies within the chemistry grid, or FALSE otherwise. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_IN_THE_CHEMGRID( I, J, L, State_Met ) 
     &                              RESULT( IS_CHEM )
!
! !USES:
!
      USE State_Met_Mod,      ONLY : MetState
      USE CMN_SIZE_MOD,       Only : LLCHEM
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
      LOGICAL                    :: IS_CHEM     ! =T if we are in the chem grid
!
! !REVISION HISTORY: 
!  24 Mar 2013 - S. D. Eastham - Initial version
!  21 Feb 2014 - M. Sulprizio  - Now pass State_Met object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! ITS_IN_THE_CHEMGRID begins here
      !=================================================================

      If (L.gt.LLCHEM) Then
         IS_CHEM = .False.
      Else
#if defined( UCX )
         IS_CHEM = ( .not. ITS_IN_THE_MESO( I, J, L ) )
#else
         IS_CHEM = ITS_IN_THE_TROP( I, J, L, State_Met )
#endif
      End If

      END FUNCTION ITS_IN_THE_CHEMGRID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_in_the_nochemgrid
!
! !DESCRIPTION: Function ITS\_IN\_THE\_NOCHEMGRID returns TRUE if grid box
!  (I,J,L) lies outside the chemistry grid, or FALSE otherwise. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_IN_THE_NOCHEMGRID( I, J, L, State_Met )
     &                              RESULT( IS_NOCHEM )
!
! !USES:
!
      USE State_Met_Mod,      ONLY : MetState
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
      LOGICAL                    :: IS_NOCHEM    ! =F if we are in the chem grid
!
! !REVISION HISTORY: 
!  24 Mar 2013 - S. D. Eastham - Initial version
!  21 Feb 2014 - M. Sulprizio  - Now pass State_Met object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! ITS_IN_THE_NOCHEMGRID begins here
      !=================================================================

      ! Function simply complement of ITS_IN_THE_CHEMGRID
      ! Provided for consistency with old approach
      IS_NOCHEM = ( .not. ITS_IN_THE_CHEMGRID( I, J, L, State_Met ) )

      END FUNCTION ITS_IN_THE_NOCHEMGRID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diag_tropopause
!
! !DESCRIPTION: Subroutine DIAG\_TROPOPAUSE archives the ND55 tropopause
!  diagnostic.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG_TROPOPAUSE( am_I_Root, Input_Opt, State_Met, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! Diagnostic switches
#if defined( BPCH_DIAG )
      USE DIAG_MOD,           ONLY : AD55
#endif
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Met_Mod,      ONLY : MetState    
#if defined( NC_DIAG )
      USE HCO_INTERFACE_MOD,  ONLY : HcoState
      USE Error_Mod,          ONLY : Error_Stop
      USE HCO_Diagn_Mod,      ONLY : Diagn_Update
      USE HCO_Error_Mod
#endif
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root    ! Is this the root CPU?
      TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
      TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(INOUT) :: RC           ! Success or failure
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
!  21 Feb 2014 - M. Sulprizio- Now pass State_Met object via the arg list
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE and GET_PCENTER with
!                              State_Met%PEDGE and State_Met%PMID.
!                              Remove dependency on pressure_mod.
!  21 Jan 2016 - E. Lundgren - Update diagnostic for netcdf output using HEMCO

!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC

!
! !LOCAL VARIABLES:
!
      !Scalars
      INTEGER :: I, J, L

      ! Met-specific variables
#if defined( GCAP )
      REAL(fp)  :: H(IIPAR,JJPAR,LLPAR)
#else
      INTEGER   :: L_TP
      REAL(fp)  :: H, FRAC, Pb, Pt
#endif

#if defined( NC_DIAG )
      ! NetCDF diagnostic variables
      INTEGER           :: M, HCRC
      CHARACTER(LEN=60) :: LOC = 'Diag_Tropopause (chemgrid_mod.F)' 
      CHARACTER(LEN=60) :: DiagnName, MSG
      REAL(fp), TARGET  :: TRPAUSE_L(IIPAR,JJPAR)
      REAL(fp), TARGET  :: TRPAUSE_P(IIPAR,JJPAR)
      REAL(fp), TARGET  :: TRPAUSE_H(IIPAR,JJPAR)
      REAL(fp), POINTER :: Ptr2D(:,:)
#endif

      !=======================================================================
      ! DIAG_TROPOPAUSE begins here!
      !=======================================================================

      ! Assume successful return
      RC = GC_SUCCESS

      ! Exit if diagnostic is not turned on
#if defined( BPCH_DIAG )
      IF ( Input_Opt%ND55 <= 0 ) RETURN
#endif
#if defined( NC_DIAG )
      ! If outputting NetCDF, temporarily using ND68 as indicator since
      ! tropopause will be output as part of Met diagnostic group
      IF ( Input_Opt%ND68 <= 0 ) RETURN
      
      ! Nullify pointers
      Ptr2d => NULL()
#endif
      
#if defined( GCAP)
      !=================================================================
      ! %%%%% GCAP MET FIELDS ONLY %%%%%
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
         H(I,J,1) = State_Met%BXHEIGHT(I,J,1) / 2.e+0_fp
      ENDDO
      ENDDO

      ! Add to H 1/2 of the sum of the two adjacent boxheights
      DO L = 1, LLPAR-1
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         H(I,J,L+1) = H(I,J,L) + 
     &               ( State_Met%BXHEIGHT(I,J,L) +
     &                 State_Met%BXHEIGHT(I,J,L+1) ) / 2.e+0_fp
      ENDDO
      ENDDO
      ENDDO

      !=================================================================
      ! ND55: Tropopause level, height [ km ], and pressure [ mb ]
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Get the tropopause level
         L           = GET_TPAUSE_LEVEL( I, J, State_Met )

         ! If we are using the variable tropopause, then (I,J,L) is the
         ! highest purely tropospheric grid box.  The grid box in which
         ! the tropopause actually occurs is then (I,J,L+1).
         IF ( LVARTROP ) L = L + 1

#if defined( BPCH_DIAG )
         ! Archive level at which tropopause occurs
         AD55(I,J,1) = AD55(I,J,1) + L

         ! Archive tropopause height [km]
         AD55(I,J,2) = AD55(I,J,2) + H(I,J,L) / 1.0e+3_fp ! m --> km

         ! We have to estimate the tropopause pressure from the 
         ! tropopause level.  Archive the pressure at the midpoint
         ! of the level in which the tropopause occurs.  NOTE: this may
         ! result in lower minimum tropopause pressure than reality.
         AD55(I,J,3) = AD55(I,J,3) + State_Met%PMID(I,J,L)
#endif
#if defined( NC_DIAG )
         TRPAUSE_L(I,J) = L                       ! tropopause level
         TRPAUSE_H(I,J) = H(I,J,L) / 1.0e+3_fp    ! tropopause height [km]
         TRPAUSE_P(I,J) = State_Met%PMID(I,J,L)   ! tropopause press [hPa]
#endif
      ENDDO
      ENDDO
#else
      !=================================================================
      ! %%%%% GEOS-4, GEOS-5, GEOS-FP, MERRA, MERRA2 met fields %%%%%
      !
      ! We get tropopause pressure directly from the met field archive
      ! Compute tropopause height to be consistent w/ the pressure
      !=================================================================

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
         Pb   = State_Met%PEDGE( I, J, L_TP   )  
         Pt   = State_Met%PEDGE( I, J, L_TP+1 )

         ! FRAC is the fraction of the grid box (I,J,L_TP) 
         ! that lies totally within the troposphere
         FRAC = ( Pb - State_Met%TROPP(I,J) ) / ( Pb - Pt ) 

         ! Add to H the height [m] of the purely tropospheric 
         ! fraction of grid box (I,J,L_TP)
         H    = H + ( FRAC * State_Met%BXHEIGHT(I,J,L_TP) )

#if defined( BPCH_DIAG )
         !---------------------------
         ! Archive into ND55 array
         !---------------------------
         AD55(I,J,1) = AD55(I,J,1) + L_TP        ! T'pause level
         AD55(I,J,2) = AD55(I,J,2) + H/1.0e+3_fp     ! T'pause height [km]
         AD55(I,J,3) = AD55(I,J,3) +             ! T'pause pressure [hPa]
     &                    State_Met%TROPP(I,J)
#endif
#if defined( NC_DIAG )
         TRPAUSE_L(I,J) = FLOAT( L_TP )          ! tropopause level
         TRPAUSE_H(I,J) = H / 1.0e+3_fp          ! tropopause height [km]
         TRPAUSE_P(I,J) = State_Met%TROPP(I,J)   ! tropopause press [hPa]
#endif
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#if defined( NC_DIAG )
      ! Assume success
      HCRC = HCO_SUCCESS

      ! Loop over tropopause diagnostics
      DO M = 1, 3
   
         ! Define diagnostic name and update diagnostic container
         SELECT CASE ( M )
            CASE ( 1 )
               DiagnName = 'MET_2D_TRPAUSE_LVL'
               Ptr2D => TRPAUSE_L
            CASE ( 2 )
               DiagnName = 'MET_2D_TRPAUSE_HGHT'
               Ptr2D => TRPAUSE_H
            CASE ( 3 )
               DiagnName = 'MET_2D_TRPAUSE_PRESS'
               Ptr2D => TRPAUSE_P
         END SELECT

         CALL Diagn_Update( am_I_Root, HcoState, 
     &                      cName     = TRIM( DiagnName ),            
     &                      Array2D   = Ptr2D,                        
     &                      RC        = HCRC )
   
         ! Free the pointers
         Ptr2D => NULL()
   
         ! Stop with error if the diagnostics update was unsuccessful
         IF ( HCRC /= HCO_SUCCESS ) THEN
            MSG = 'Cannot update diagnostic: ' // TRIM( DiagnName )
            CALL ERROR_STOP( MSG, LOC ) 
         ENDIF 
 
      ENDDO
#endif

#endif

      END SUBROUTINE DIAG_TROPOPAUSE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_chemgrid
!
! !DESCRIPTION: Subroutine INIT\_CHEMGRID allocates and zeroes module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_CHEMGRID( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE Input_Opt_Mod,      ONLY : OptInput
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
      ! INIT_CHEMGRID begins here!
      !=================================================================

      ! Assume success
      RC        =  GC_SUCCESS

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

      END SUBROUTINE INIT_CHEMGRID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_chemgrid
!
! !DESCRIPTION: Subroutine CLEANUP\_CHEMGRID deallocates module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_CHEMGRID
! 
! !REVISION HISTORY: 
!  10 Feb 2005 - R. Yantosca - Initial version
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_CHEMGRID begins here!
      !=================================================================
      IF ( ALLOCATED( TROPOPAUSE ) ) DEALLOCATE( TROPOPAUSE )
      IF ( ALLOCATED( IFLX       ) ) DEALLOCATE( IFLX       ) 

      END SUBROUTINE CLEANUP_CHEMGRID 
!EOC
      END MODULE CHEMGRID_MOD
