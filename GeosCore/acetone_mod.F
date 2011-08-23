!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: acetone_mod
!
! !DESCRIPTION: Module ACETONE\_MOD contains subroutines to emit the biogenic 
!  flux of acetone into the full chemistry simulation.
!\\
!\\
! !INTERFACE:
!
      MODULE ACETONE_MOD
!
! !USES:
!
      IMPLICIT NONE
#     include "define.h"
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_ACETONE   
      PUBLIC :: EMISS_BIOACET
      PUBLIC :: OCEAN_SOURCE_ACET
      PUBLIC :: OCEAN_SINK_ACET
      PUBLIC :: READ_JO1D
      PUBLIC :: READ_RESP
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Jacob, D.J., B.D. Field, E. Jin, I. Bey, Q. Li, J.A. Logan, and 
!        R.M. Yantosca, "Atmospheric budget of acetone", Geophys. Res. Lett., 
!        107(D11), 4100, 2002. 
!  (2 ) Nightingale et al [2000a], J. Geophys. Res, 14, 373-387
!  (3 ) Nightingale et al [2000b], Geophys. Res. Lett, 27, 2117-2120
!
! !REVISION HISTORY:
!  NOTES:
!  18 Sep 2001 - B. Field, R. Yantosca - Initial version  
!  (1 ) Added changes from bdf and updated comments (bmy, 9/5/01)
!  (2 ) Updated comments (bmy, 9/12/01)
!  (3 ) Removed VERBOSE flag and all "print-to-log-file" diagnostics.  The
!        ND11 diagnostic produces the same totals. (bdf, bmy, 9/18/01)
!  (4 ) Now cal GET_TAU0 w/ 3 arguments instead of 2.  Also minor bug
!        fix in READ_RESP (bmy, 11/15/01)
!  (5 ) Implemented fix for ocean source/sink from Mat Evans.  Also deleted 
!        obsolete code from 11/01. (bmy, 11/26/01)
!  (6 ) Eliminated more obsolete code from 11/01 (bmy, 2/27/02)
!  (7 ) Removed duplicate variable definitions (bmy, 3/20/02)
!  (8 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (9 ) Bug fix: Now apply true exponential loss in OCEAN_SINK_ACET, instead
!        of just the 1st order approximation. (bdf, bmy, 7/11/02)
!  (10) Scale the ocean source of acetone for GEOS-3 meteorology in order to
!        match the total listed in Jacob et al 2002. (bdf, bmy, 9/16/02)
!  (11) Now references "error_mod.f" (bmy, 10/15/02)
!  (12) Minor modifications to READ_JO1D, READ_RESP (bmy, 3/14/03)
!  (13) Add surface area scale factor for ocean source for 1x1 nested
!        grids.  (yxw, bmy, 5/16/03)
!  (14) Scale ACET ocean source to Jacob et al 2002 for GEOS-4, and now
!        account for surface area ratio for all GEOS grids. (bmy, 3/15/04)
!  (15) Now references "directory_mod.f" (bmy, 7/19/04)
!  (16) Now can read data from GEOS and GCAP grids.  Also now use Nightingale
!        et al 2000b formulation for piston velocity KL. (swu, bmy, 8/16/05)
!  (17) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (18) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (19) Updates for nested EU and NA grids (amv, bmy, 12/18/09)
!  (20) Updates for GEOS-4 1 x 1.25 grid (lok, bmy, 1/13/10)
!  13 Aug 2010 - R. Yantosca - Add modifications for MERRA (treat like GEOS-5)
!  04 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
! 
      ! Array for J-O1D values
      REAL*8,  ALLOCATABLE :: JO1D(:,:)

      ! Array for heterotrophic respiration values
      REAL*8,  ALLOCATABLE :: XRESP(:,:) 
!
! !DEFINED PARAMETERS:
!
      ! Avogadro's number
      REAL*8,  PARAMETER   :: AVO      = 6.022d23 

      ! Molecules C / kg C 
      REAL*8,  PARAMETER   :: XNUMOL_C = AVO / 12d-3 

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_jo1d
!
! !DESCRIPTION: Subroutine READ\_JO1D reads in the J-Values for O1D from disk 
!  that are needed for the biogenic acetone fluxes,
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_JO1D( THISMONTH )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD                 ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: THISMONTH   ! Current month
!
! !REMARKS:
!  J-values for O1D are are stored in the JO1D module array in [s^-1].
! 
! !REVISION HISTORY: 
!  14 Sep 2001 - B. Field    - Initial version
!  (1 ) Now use TRANSFER_2D from "transfer_mod" to cast from REAL*4 to REAL*8
!        and to resize data block to (IIPAR,JJPAR).  Also use 3-argument
!        form of GET_TAU0 (bmy, 11/15/01)
!  (2 ) Removed obsolete code from 11/01 (bmy, 11/26/01, bmy, 2/27/02)
!  (3 ) Now reference routines ALLOC_ERR and ERROR_STOP from "error_mod.f" 
!        (bmy, 10/15/02)
!  (4 ) Now call READ_BPCH2 with QUIET=.TRUE. to suppress printing of extra
!        info to stdout.  Also made cosmetic changes. (bmy, 3/14/03)
!  (5 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/19/04)
!  (6 ) Now can read data from GEOS and GCAP grids (bmy, 8/16/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  04 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST = .TRUE.
      INTEGER            :: AS
      REAL*4             :: ARRAY(IIPAR,JJPAR,1)
      REAL*8             :: TAU0
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! READ_JO1D begins here!
      !=================================================================

      ! Allocate JO1D only on the first call
      IF ( FIRST ) THEN
         ALLOCATE( JO1D( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'JO1D' )
         JO1D  = 0d0
         FIRST = .FALSE.
      ENDIF

      ! Construct filename
      FILENAME = TRIM( DATA_DIR )       // 
     &           'acetone_200108/JO1D.' // GET_NAME_EXT_2D() //
     &           '.'                    // GET_RES_EXT()

      ! Echo filename 
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_JO1D: Reading ', a )  

      ! Get TAU0 value -- use "generic" year 1985
      TAU0 = GET_TAU0( THISMONTH, 1, 1985 )

      ! Read JO1D data 
      CALL READ_BPCH2( FILENAME, 'JV-MAP-$',    51, 
     &                 TAU0,      IIPAR,        JJPAR,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )
   
      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), JO1D )

      ! Make sure JO1D is not zero everywhere 
      IF ( .not. ( SUM( JO1D ) > 0d0 ) ) THEN
         CALL ERROR_STOP( 'JO1D is zero everywhere!', 'READ_JO1D' )
      ENDIF
         
      ! Echo information
      WRITE( 6, 110 ) MAXVAL( JO1D )
 110  FORMAT( '     - READ_JO1D: Max value of JO1D: ', es13.6 )

      END SUBROUTINE READ_JO1D
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_resp
!
! !DESCRIPTION: Subroutine READ\_RESP reads in the monthly heterotrophic 
!  respiration measured in g of plant matter/m\^2 flowing out of the biosphere.
!  This is needed for the biogenic acetone fluxes.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_RESP( THISMONTH )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD                  ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: THISMONTH   ! Current month
!
! !REMARKS:
!  Respiration values are stored in the RESP module array in [g C/m2/s].
! 
! !REVISION HISTORY: 
!  14 Sep 2001 - B. Field    - Initial version
!  (1 ) Now use TRANSFER_2D from "transfer_mod" to cast from REAL*4 to REAL*8
!        and to resize data block to (IIPAR,JJPAR).  Bug fix: THISMONTH > 12
!        is never valid.  Also use 3-argument form of GET_TAU0 (bmy, 11/15/01)
!  (2 ) Removed obsolete code from 11/01 (bmy, 11/26/01, bmy, 2/27/02)
!  (3 ) Now reference ALLOC_ERR from "error_mod.f".  Also use version of
!        GET_TAU0 w/ 3 arguments. (bmy, 10/15/02)
!  (4 ) Now call READ_BPCH2 with QUIET=.TRUE. to suppress printing of extra
!        info to stdout.  Also made cosmetic changes. (bmy, 3/14/03)
!  (5 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (6 ) Now can read files for both GEOS and GCAP grids (bmy, 8/16/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  04 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST = .TRUE.
      INTEGER            :: AS
      REAL*4             :: ARRAY(IIPAR,JJPAR,1)
      REAL*8             :: TAU0, TAU1
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! READ_RESP begins here!
      !=================================================================
      
      ! Allocate XRESP only on the first call
      IF ( FIRST ) THEN
         ALLOCATE( XRESP( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'XRESP' )
         XRESP = 0d0
         FIRST = .FALSE.
      ENDIF

      ! Construct filename
      FILENAME = TRIM( DATA_DIR )       // 
     &           'acetone_200108/resp.' // GET_NAME_EXT_2D() //
     &           '.'                    // GET_RES_EXT()
      
      ! Echo filename 
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_RESP: Reading ', a )  

      ! Get TAU0 value -- use "generic" year 1985
      TAU0 = GET_TAU0( THISMONTH, 1, 1985 )

      ! Get TAU1 value -- use generic year 1985 
      ! Be careful -- if Dec then nextmonth is Jan 1986
      IF ( THISMONTH == 12 ) THEN
         TAU1 = GET_TAU0( 1, 1, 1986 )
      ELSE
         TAU1 = GET_TAU0( THISMONTH+1, 1, 1985 )
      ENDIF

      ! Read heterotrophic respiration data in [g C/m2/month]
      CALL READ_BPCH2( FILENAME, 'HET-RESP',    1, 
     &                 TAU0,      IIPAR,        JJPAR,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )
   
      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), XRESP )

      ! Convert RESP field to [g C/m2/s]
      XRESP = XRESP / ( ( TAU1 - TAU0 ) * 3600d0 )

      END SUBROUTINE READ_RESP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ocean_source_acet
!
! !DESCRIPTION: Subroutine OCEAN\_SOURCE\_ACET specifies the ocean source of 
!  acetone.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OCEAN_SOURCE_ACET( I, J, ACETONE )
!
! !USES:
!
      USE ERROR_MOD, ONLY : CHECK_VALUE
      USE DAO_MOD,   ONLY : ALBD, TS
      USE DIAG_MOD,  ONLY : AD11
      USE GRID_MOD,  ONLY : GET_AREA_CM2
      USE TIME_MOD,  ONLY : GET_TS_EMIS

      USE CMN_SIZE_MOD                  ! Size parameters
      USE CMN_DIAG_MOD                  ! ND11
      USE CMN_DEP_MOD                   ! FRCLND
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: I         ! GEOS-Chem longitude index
      INTEGER, INTENT(IN)    :: J         ! GEOS-Chem latitude index
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(INOUT) :: ACETONE   ! Acetone emission [atoms C/s]
!
! !REMARKS:

! 
! !REVISION HISTORY: 
!  14 Sep 2001 - B. Field    - Initial version
!  (1 ) Now compute u = SQRT( U10M^2 + V10M^2 ) as SQRT( SFCWINDSQR(I,J) ).
!        This is necessary since U10M and V10M are missing for 1996, and
!        need to be computed from UWND and VWND.  (bmy, 9/5/01)
!  (2 ) Bug fixes: multiply kg by 360000 and use exponent to the -0.5 power
!        in the expression for Kl.  Also update value of the OCEAN_SCALE
!        factor to 3.63e11.  Also updated comments. (bdf, bmy, 9/5/01)
!  (3 ) Bug fix: ACETONE has units of [atoms C/box/s], to match those of
!        EMISRR.  This involves an extra division by DTSRCE. (bmy, 9/14/01)
!  (4 ) Removed diagnostic variable OCEAN_SOURCE (bmy, 9/18/01)
!  (5 ) JO1D(IREF,JREF) is now JO1D(I,J).  Bug fix: Zero the ocean source
!        of acetone in grid boxes that are covered by less than 50% ocean.  
!        Bug fix: make sure -5 <= TC <= 30, in order to prevent the power
!        series for Schmidt # from going negative.  Also eliminate IREF,
!        JREF, we don't need them anymore. (mje, rvm, bmy, 11/26/01)
!  (6 ) Eliminated obsolete code from 11/01 (bmy, 2/27/02)
!  (7 ) Scale the ocean source of acetone for GEOS-3 meteorology in order to
!        match the total listed in Jacob et al 2002. (bdf, bmy, 9/16/02)
!  (8 ) Now use function GET_AREA_CM2 of "grid_mod.f" to return the
!        grid box area in cm2.  Use function GET_TS_EMIS from "time_mod.f".
!        Remove reference to CMN header file. (bmy, 2/11/03)
!  (9 ) Apply surface area scale factor for 1x1 nested grids, in order to
!        make the total ocean source the same as for 4x5. (yxw, bmy, 5/16/03)
!  (10) Scale the ocean source to Jacob et al 2002 for GEOS-4.  Also account
!        for surface area ratio for all GEOS grids. (bmy, 3/15/04)
!  (11) Added space in #ifdef block for GEOS-4 x 1x125 grid (bmy, 12/1/04)
!  (12) Now use Nightingale et al 2000b formulation for piston velocity KL.
!        (swu, bmy, 8/16/05)
!  (13) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (14) Adjust SCALE_FACTOR for 0.5 x 0.667 grid (dan, bmy, 11/6/08)
!  (15) Additional scale factors for NESTED_NA and NESTED_EU calculated and 
!        included (amv, bmy, 12/18/09)
!  (16) Added scale factor for GEOS-4 1 x 1.25 grid (lok, bmy, 1/13/10)
!  13 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  04 Nov 2010 - R. Yantosca - Added ProTeX headers
!  04 Nov 2010 - R. Yantosca - Cleaned up #if statements for clarity
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8                 :: KG, U, TC, SC, kl, KKL, HSTAR, KL600
      REAL*8                 :: DTSRCE, OCEAN_ACET, AREA_CM2, FOCEAN
      REAL*8,  EXTERNAL      :: SFCWINDSQR
!
! !DEFINED PARAMETERS:
!
      ! ALPHA scale factor in kg*s^2/cm
      REAL*8,  PARAMETER     :: OCEAN_SCALE = 3.63d11     

      ! Coefficients for fitting the Schmidt number for acetone
      REAL*8,  PARAMETER     :: A0 =  3287.687d0
      REAL*8,  PARAMETER     :: A1 = -136.2176d0
      REAL*8,  PARAMETER     :: A2 =  2.20642d0
      REAL*8,  PARAMETER     :: A3 = -0.01410642d0

      !=================================================================
      ! Since the algorithm below was developed for the 4x5 GEOS-STRAT
      ! model, we need to scale the emissions to the A Posterioris in 
      ! Jacob et al 2002 = 16.74 Tg C/yr.
      !
      ! This scaling is necessary in order (1) to account for variations
      ! in surface temperature and wind speed between the different
      ! GEOS met field versions and (2) to account for the different
      ! surface areas between the 1x1, 2x2.5 and 4x5 grid boxes.
      !=================================================================
#if   defined( MERRA ) && defined( GRID4x5 ) 

      !-----------------------------------------------------------------
      ! %%% MERRA, 4 x 5, global %%%
      !
      ! Acetone produced : 17.0919 Tg C
      ! Scale factor     : 16.74 / 17.0919 = 0.979411
      !-----------------------------------------------------------------
      REAL*8, PARAMETER :: SCALE_FACTOR = 0.979411d0

#elif defined( MERRA ) && defined( GRID2x25 )

      !-----------------------------------------------------------------
      ! %%% MERRA, 2 x 2.5, global %%%
      !
      ! %%% SCALING NEEDS TO BE DONE WHEN MET IS FINALLY AVAILABLE %%%
      !-----------------------------------------------------------------
      REAL*8, PARAMETER :: SCALE_FACTOR = 1d0

#elif defined( GEOS_5 ) && defined( GRID4x5 )
 
      !-----------------------------------------------------------------
      ! %%% GEOS-5, 4 x 5, global %%%
      !
      ! Acetone produced : 17.5272 Tg C
      ! Scale factor     : 16.74 / 17.5272 = 0.9551
      !-----------------------------------------------------------------
      REAL*8, PARAMETER :: SCALE_FACTOR =   0.9551d0

#elif defined( GEOS_5 ) && defined( GRID2x25 ) 

      !-----------------------------------------------------------------
      ! %%% GEOS-5, 2 x 2.5, global %%%
      !
      ! Acetone produced : 75.6781 Tg C
      ! Scale factor     : 16.74 / 75.6781 = 0.2212
      !-----------------------------------------------------------------
      REAL*8, PARAMETER :: SCALE_FACTOR = 0.2212d0

#elif defined( GEOS_5 ) && defined( GRID05x0666 ) && defined ( NESTED_CH )

      !-----------------------------------------------------------------
      ! %%% GEOS-5, 0.5 x 0.666, NESTED CHINA %%%     (May Fu, Dan Chen)
      !
      ! ACET; 4x5 BC region      : 2.16838 Tg C
      ! ACET; 0.5 x 0.666 region : 134.757 Tg C
      ! Scale factor             : 0.9551*(2.16838/134.757) = 0.01537
      !-----------------------------------------------------------------
      REAL*8, PARAMETER :: SCALE_FACTOR = 0.015369d0      

#elif defined( GEOS_5 ) && defined( GRID05x0666 ) && defined ( NESTED_NA )

      !-----------------------------------------------------------------
      ! %%% GEOS-5, 0.5 x 0.666, NESTED N. AMERICA %%%   (van Donkelaar)
      !
      ! ACET; 2x25 BC region     : 1.6794 Tg C
      ! ACET; 0.5 x 0.666 region : 1.724  Tg C
      ! Scale factor             : 0.015369*(1.6794/1.724) = 0.01497
      !-----------------------------------------------------------------
      REAL*8, PARAMETER :: SCALE_FACTOR = 0.01497d0
 
#elif defined( GEOS_5 ) && defined( GRID05x0666 ) && defined ( NESTED_EU )

      !-----------------------------------------------------------------
      ! %%% GEOS-5, 0.5 x 0.666, NESTED EUROPE %%%       (van Donkelaar)
      !
      ! ACET; 2x25 BC region     : 0.2868 Tg C
      ! ACET; 0.5 x 0.666 region : 0.2689 Tg C
      ! Scale factor             : 0.015369*(0.2868/0.2689) = 0.01639
      !-----------------------------------------------------------------
      REAL*8, PARAMETER :: SCALE_FACTOR = 0.01639d0

#elif defined( GEOS_4 ) && defined( GRID4x5 )
 
      !-----------------------------------------------------------------
      ! %%% GEOS-4, 4 x 5, global %%%
      !
      ! Acetone produced : 19.0973 Tg C
      ! Scale factor     : 16.74 / 19.0973 = 0.8765
      !-----------------------------------------------------------------
      REAL*8, PARAMETER :: SCALE_FACTOR = 0.8765d0
 
#elif defined( GEOS_4 ) && defined( GRID2x25 )
 
      !-----------------------------------------------------------------
      ! %%% GEOS-4, 2 x 2.5, global %%%
      !
      ! Acetone produced : 80.2888 Tg C
      ! Scale factor     : 16.74 / 80.2888 = 0.2085
      !-----------------------------------------------------------------
      REAL*8, PARAMETER :: SCALE_FACTOR = 0.2085d0
 
#elif defined( GEOS_4 ) && defined( GRID1x125 )
 
      !-----------------------------------------------------------------
      ! %%% GEOS-4, 1 x 1.25, global %%%                    (Lok Lamsal)
      !
      ! Acetone produced : 290.4233 Tg C
      ! Scale factor     : 16.74 / 290.4233 = 0.05764
      !-----------------------------------------------------------------
      REAL*8, PARAMETER :: SCALE_FACTOR = 0.05764d0
 
#elif  defined( GEOS_3 ) && defined( GRID4x5 )
 
      !-----------------------------------------------------------------
      ! %%% GEOS-3, 4 x 5, global %%%
      !
      ! Acetone produced : 20.16 Tg C
      ! Scale factor     : 16.74 / 20.16 =  0.83
      !-----------------------------------------------------------------
      REAL*8, PARAMETER :: SCALE_FACTOR = 0.83d0
 
#elif defined( GEOS_3 ) && defined( GRID2x25 )
 
      !-----------------------------------------------------------------
      ! %%% GEOS-3, 2 x 2.5, global %%%
      !
      ! Acetone produced : 80.76 Tg C
      ! Scale factor     : 16.74 / 80.76  = 0.2075
      !-----------------------------------------------------------------
      REAL*8, PARAMETER :: SCALE_FACTOR = 0.2075d0
 
#elif defined( GEOS_3 ) && defined( GRID1x1 )
 
      !-----------------------------------------------------------------
      ! %%% GEOS-3, 1 x 1, NESTED ?? %%%
      ! Acetone produced : 334.8 Tg C
      ! Scale factor     : 16.74 / 334.8  = 0.05
      !-----------------------------------------------------------------
      REAL*8, PARAMETER :: SCALE_FACTOR = 0.05d0
 
#else
      
      ! Otherwise set to 1
      REAL*8, PARAMETER :: SCALE_FACTOR = 1d0
 
#endif

      !=================================================================
      ! OCEAN_SOURCE_ACET begins here!
      !=================================================================

      ! Emission timestep in seconds
      DTSRCE   = GET_TS_EMIS() * 60d0 

      ! Fraction of (I,J) that is ocean
      FOCEAN   = 1d0 - FRCLND(I,J)    

      ! Area of grid box (I,J) in cm^2
      AREA_CM2 = GET_AREA_CM2( J )    

      !=================================================================
      ! Compute ocean source by Henry's law
      ! Now make sure only to compute the ocean source if more than
      ! 50% of the box is covered by water (mje, bdf, bmy, 11/26/01)
      !=================================================================
      IF ( FOCEAN > 0.5d0 ) THEN 

         ! Henry's law [unitless] using 32 M/atm Henry's law constant
         HSTAR    = ( 1d0 / 730d0 ) *
     &        EXP( -5500d0 * ( 298d0 - TS(I,J) ) / ( TS(I,J) * 298d0 ) )

         ! Magnitude of resultant wind [m/s]
         ! SFCWINDSQR(I,J) is needed since this will compute the square
         ! of the surface wind correctly for all GEOS models (bmy, 9/5/01)
         U        = SQRT( SFCWINDSQR( I, J ) )

         ! TC is temperature in Celsius
         ! Also make sure -5 <= TC <= 30 (mje, rvm, bmy, 11/26/01)
         TC       = TS(I,J) - 273.15d0                           
         TC       = MIN( MAX( TC, -5d0 ), 30d0 )

         ! SC is Schmidt # for acetone [unitless]
         SC       = A0 + A1*TC + A2*TC**2 + A3*TC**3 

         ! KL is conductance for mass transfer in liquid phase 
         ! (Nightingale et al 2000b), which has units of [cm/hr]
         KL       = ( 0.24d0*U*U + 0.061d0*U ) * SQRT( 600d0/Sc )  

         ! KG is conductance for mass transfer in gas phase (Asher 1997)
         ! Multiply KG by 360000 to convert from [m/s] to [cm/hr]
         KG       = SQRT( 18d0 / 58d0 ) * ( 5.2d-5 + 3.2d-3 * U )   
         KG       = KG * 360000d0

         ! KKL is the air-to-sea transfer velocity (Liss and Slater 1974)
         ! Multiply KKL by 3600 to convert from [cm/hr] to [cm/s]
         KKL      = 1d0 / ( 1d0/KL + 1d0/( HSTAR * KG ) ) 
         KKL      = KKL / 3600d0

         ! Turn off ocean uptake over snow/ice (ALBEDO > 0.4)
         IF ( ALBD(I,J) > 0.4d0 ) KKL = 0d0
      
         !==============================================================
         ! The acetone ocean source [kg C] includes the following terms:
         !
         !   (1) OCEAN_SCALE = 3.63e11 [kg C s^2/cm].  
         !
         !   (2) J(O1D), which is a proxy for the UV-B flux [s^-1].  
         !
         !       The product of J(O1D) and OCEAN_SCALE reflect the 
         !       photoproduction of acetone by microbes in the ocean
         !       microlayer.
         !
         !   (3) KKL, the air-to-sea transfer velocity [cm/s], which
         !       represents uptake of acetone by microbes in the ocean
         !       microlayer.  This term is a net loss.
         !
         !   (4) FOCEAN, the fraction of the grid box that is ocean.
         !
         ! The units of the ocean source (stored in OCEAN_ACET) are:
         !
         !      1     |   kg C   | s^2 | 1 | cm       kg C   | 
         !  ----------+----------+-----+---+---- = ----------+----------
         !  emission  | grid box | cm  | s | s      grid box | emission 
         !  time step                                          time step
         !==============================================================
         OCEAN_ACET = OCEAN_SCALE * JO1D(I,J) * KKL * FOCEAN 
         
         ! Apply further scale factor to account for variations in surface 
         ! temperature wind speed between GEOS met fields -- and also 
         ! surface area between 1x1, 2x2.5, and 4x5 grids. (bmy, 3/15/04)
         OCEAN_ACET = OCEAN_ACET * SCALE_FACTOR

      ELSE

         ! If there is less than 50% water in the grid box, zero 
         ! the ocean source from acetone (mje, rvm, bmy, 11/26/01)
         OCEAN_ACET = 0d0
         
      ENDIF

      ! Add ocean source to total biogenic source in [atoms C/box/s]
      ACETONE = ACETONE + ( OCEAN_ACET * XNUMOL_C / DTSRCE )
            
      !=================================================================
      ! ND11 diag -- save ACETONE from the ocean in [atoms C/cm2/s]
      !=================================================================
      IF ( ND11 > 0 ) THEN
         AD11(I,J,6) = AD11(I,J,6) + ( OCEAN_ACET * XNUMOL_C ) / 
     &                               ( AREA_CM2   * DTSRCE   )
      ENDIF

      END SUBROUTINE OCEAN_SOURCE_ACET
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ocean_sink_acet
!
! !DESCRIPTION: Subroutine OCEAN\_SINK\_ACET applies the ocean sink to global
!  acetone concentrations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OCEAN_SINK_ACET( ACETONE )
!
! !USES:
!
      USE DAO_MOD,  ONLY : ALBD
      USE DAO_MOD,  ONLY : TS
      USE DIAG_MOD, ONLY : AD11
      USE GRID_MOD, ONLY : GET_AREA_CM2
      USE TIME_MOD, ONLY : GET_TS_CHEM

      USE CMN_SIZE_MOD                              ! Size parameters 
      USE CMN_DIAG_MOD                              ! ND11
      USE CMN_DEP_MOD                               ! FRCLND 
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8, INTENT(INOUT) :: ACETONE(IIPAR,JJPAR)   ! Acetone mass [kg C]
! 
! !REVISION HISTORY: 
!  14 Sep 2001 - B. Field    - Initial version
!  (1 ) Remove references to F77_CMN_UV10M and F77_CMN_LWI -- these are now
!        obsolete in GEOS-CHEM versions 4.18 and higher (bmy, 9/5/01)
!  (2 ) Now compute u = SQRT( U10M^2 + V10M^2 ) as SQRT( SFCWINDSQR(I,J) ).
!        This is necessary since U10M and V10M are missing for 1996, and
!        need to be computed from UWND and VWND.  (bmy, 8/2/01)
!  (3 ) Now declare OCEANSINK_SCALE = 0.15 as a parameter.  This is the
!        optimized value of BETA from Emily Jin's analysis.  Also updated
!        comments. (bdf, bmy, 9/5/01)
!  (4 ) Updated comments.  Also parallellized DO loops. (bmy, 9/14/01)
!  (5 ) Removed diagnostic variable OCEAN_LOSS (bmy, 9/18/01)
!  (6 ) Bug fix: Zero the ocean sink of acetone in grid boxes where there
!        is less than 50% of ocean, and where there is ice on the surface.
!        Bug fix: Make sure -5 <= TC <= 30, in order to prevent the power
!        series for Schmidt # from going negative.  Also eliminate IREF,
!        JREF, we don't need them. (mje, rvm, bmy, 11/26/01)
!  (7 ) Eliminated obsolete code from 11/01 (bmy, 2/27/02)
!  (8 ) Bug fix: now use true exponential for loss instead of just 1st
!        order term.  Also added PRE_ACET variable to save previous acetone
!        mass for diagnostic, before applying loss.  (bdf, bmy, 7/11/02)
!  (9 ) Now use function GET_AREA_CM2 of "grid_mod.f" to return the
!        grid box area in cm2.  Now use function GET_TS_CHEM from
!        "time_mod.f".  Remove reference to CMN header file. (bmy, 2/11/03)
!  (12) Now use Nightingale et al 2000b formulation for piston velocity KL.
!        (swu, bmy, 8/16/05)
!  04 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I, IREF, J, JREF
      REAL*8            :: KH298, DHR, KH, U, TC, SC, KL, KG 
      REAL*8            :: KKL, CG, F, T1L, H, KL600, FLUX, HSTAR
      REAL*8            :: AREA_CM2, DTCHEM, FOCEAN, OCEAN_ACET
      REAL*8            :: PRE_ACET

      ! External functions
      REAL*8, EXTERNAL  :: BOXVL, SFCWINDSQR
!
! !DEFINED PARAMETERS:
!
      ! Optimized value of BETA for ocean sink found by
      ! Emily Jin from inverse modeling analysis
      REAL*8, PARAMETER :: OCEANSINK_SCALE = 0.15d0

      ! Coefficients for fitting the Schmidt number for acetone
      REAL*8, PARAMETER :: A0 =  3287.687d0
      REAL*8, PARAMETER :: A1 = -136.2176d0
      REAL*8, PARAMETER :: A2 =  2.20642d0
      REAL*8, PARAMETER :: A3 = -0.01410642d0

      !=================================================================
      ! OCEAN_SINK_ACET begins here! 
      !
      ! Compute acetone lost to ocean sink and subtract from STT
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,  J, AREA_CM2, FOCEAN, HSTAR, U, TC, SC )   
!$OMP+PRIVATE( KL600, KL, KG, KKL, CG, FLUX, OCEAN_ACET, PRE_ACET )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR

         ! Grid box area in cm2
         AREA_CM2 = GET_AREA_CM2( J )

         DO I = 1, IIPAR

            ! Fraction of grid box that is ocean
            FOCEAN = 1d0 - FRCLND(I,J)

            !===========================================================
            ! Only compute ocean sink if there is more than 50% ocean
            ! in the grid box, and if it is not ice (albedo > 0.4)
            ! (mje, rvm, bmy, 11/26/01)
            !===========================================================
            IF ( FOCEAN > 0.5d0 .and. ALBD(I,J) <= 0.4d0 ) THEN

               ! Henry's law [unitless] using 32 M/atm Henry's law constant
               HSTAR = ( 1d0 / 792d0 ) * 
     &        EXP( -5500d0 * ( 298d0 - TS(I,J) ) / ( TS(I,J) * 298d0 ) )

               ! Magnitude of surface wind [m/s]
               ! SFCWINDSQR(I,J) is needed since this will compute the 
               ! square of the surface wind correctly for all GEOS models 
               U     = SQRT( SFCWINDSQR( I, J ) )
            
               ! TC is temperature in Celsius
               ! Bug fix: make sure -5 <= TC <= 30 (mje, rvm, bmy, 11/26/01) 
               TC    = TS(I,J) - 273.15d0  
               TC    = MIN( MAX( TC, -5d0 ), 30d0 )

               ! SC is Schmidt # for acetone [unitless]
               SC    = A0 + A1*TC + A2*TC**2 + A3*TC**3 

               ! KL is conductance for mass transfer in liquid phase 
               ! (Nightingale et al 2000b), which has units of [cm/hr]
               KL    = ( 0.24d0*U*U + 0.061d0*U ) * SQRT( 600d0/Sc )  

               ! KG is conductance for mass transfer in gas phase (Asher 1997)
               ! Multiply KG by 360000 to convert from [m/s] to [cm/hr]
               KG    = SQRT( 18d0 / 58d0 ) * ( 5.2d-5 + 3.2d-3 * U ) 
               KG    = KG * ( 360000d0 )     

               ! KKL is the air-to-sea transfer velocity (Liss and Slater 1974)
               ! Multiply KKL by 3600 to convert from [cm/hr] to [cm/s]
               KKL   = 1d0 / ( 1d0/KL + 1d0/( HSTAR * KG ) ) 
               KKL   = KKL / 3600d0       

               ! CG is the gas concentration of acetone [kg C/cm3]
               CG    = ACETONE(I,J) / BOXVL(I,J,1)

               ! FLUX is the air-to-sea flux of acetone in [kg C/cm2/s].
               FLUX  = KKL * CG / HSTAR    

               ! Multiply FLUX by OCEANSINK_SCALE, which is the optimized 
               ! BETA value (= 0.15) found from Emily Jin's analysis.
               FLUX  = FLUX * OCEANSINK_SCALE

               !========================================================
               ! Ocean loss of acetone consists of the following terms:
               !
               ! (1) FLUX, the air-to-sea flux of acetone in [kg C/cm2/s]
               !
               ! (2) AREA_CM2, the grid box surface area
               !
               ! (3) DTCHEM, the number of seconds per chemistry timestep
               !
               ! (4) FOCEAN, the fraction of the grid box that is ocean.
               !
               ! The units of the resultant ocean loss (in OCEAN_ACET) are:
               !
               !     kg C  | AREA_CM2 cm2 |  DTCHEM s           kg C
               !  ---------+--------------+--------------- = ------------
               !   cm2 * s |   grid box   | chem timestep     box * step
               !========================================================
               OCEAN_ACET = ( FLUX * AREA_CM2 * DTCHEM * FOCEAN )
               
            ELSE

               ! If there is less than 50% water in the grid box, or  
               ! if there is ice on the ocean, then zero the ocean sink
               ! for acetone (mje, rvm, bmy, 11/26/01)
               OCEAN_ACET = 0d0

            ENDIF

            ! Save mass of acetone in tmp variable for diagnostic
            PRE_ACET = ACETONE(I,J) 

            ! Apply exponential loss to acetone mass
            ACETONE(I,J) = ACETONE(I,J) * EXP(-OCEAN_ACET/ACETONE(I,J))

            !===========================================================
            ! Diagnostics: save ACETONE lost to ocean in [atoms C/cm2/s]
            !===========================================================
            IF ( ND11 > 0 ) THEN
               AD11(I,J,7) = AD11(I,J,7) + 
     &              ( ( PRE_ACET - ACETONE(I,J) ) * XNUMOL_C ) / 
     &              ( AREA_CM2 * DTCHEM   )
            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE OCEAN_SINK_ACET
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_bioacet
!
! !DESCRIPTION: Subroutine EMISS\_BIOACET computes the biogenic emissions of 
!  ACETONE from monoterpenes, isoprene, methyl butenol, dry leaf matter, and 
!  grasslands.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_BIOACET( I,    J,    TMMP,  EMMO, 
     &                          EMIS, EMMB, GRASS, ACETONE )
!
! !USES:
!
      USE DIAG_MOD, ONLY : AD11
      USE GRID_MOD, ONLY : GET_AREA_M2
      USE GRID_MOD, ONLY : GET_XMID
      USE GRID_MOD, ONLY : GET_YMID
      USE TIME_MOD, ONLY : GET_TS_EMIS

      USE CMN_SIZE_MOD                  ! Size parameters
      USE CMN_DIAG_MOD                  ! ND11
      USE CMN_MONOT_MOD                 ! BASEMONOT
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: I         ! Grid box longitude index
      INTEGER, INTENT(IN)    :: J         ! Grid box latitude index
      REAL*8,  INTENT(IN)    :: TMMP      ! Surface temperature [K]
      REAL*8,  INTENT(IN)    :: EMMO      ! Monoterpene emission [atoms C]
      REAL*8,  INTENT(IN)    :: EMIS      ! Isoprene emission [atoms C]
      REAL*8,  INTENT(IN)    :: EMMB      ! Methylbutenol emission  [atoms C]
      REAL*8,  INTENT(IN)    :: GRASS     ! Isoprene from grasslands [atoms C]
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(INOUT) :: ACETONE   ! Biogenic acetone [atoms C/s]
! 
! !REVISION HISTORY: 
!  18 Sep 2001 - B. Field    - Initial version
!  (1 ) Now pass acetone array (e.g. from STT) thru the argument list, since
!        this avoids dependence on IDTACET in this program (bmy, 8/1/01)
!  (2 ) Updated scale factors (bdf, bmy, 9/5/01)
!  (3 ) Updated comments (bmy, 9/14/01)
!  (4 ) Removed diagnostic variables: MONOTERPENES, ISOPRENE, ISOP_TOTAL,
!        MONO_TOTAL, NA_TOT, RESP_TOT, GRASS_TOT.  These have now been 
!        supplanted by the ND11 acetone source diagnostic. (bdf, bmy, 9/18/01)
!  (5 ) XRESP(I+I0,J+J0) is now XRESP(I,J) (bmy, 11/26/01)
!  (6 ) Eliminated obsolete code from 11/01 (bmy, 2/27/02)
!  (7 ) Removed duplicate definitions of EMMB and GRASS (bmy, 3/20/02)
!  (8 ) Now use functions from "grid_mod.f" to get surface area, lon, and
!        lat of grid box (I,J).   Use function GET_AREA_M2 to get the grid
!        box surface area in m2, then convert to cm2.   Now use function
!        GET_TS_EMIS from "time_mod.f".  Remove reference to CMN header
!        file. (bmy, 2/11/03)
!  04 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8                 :: EMMO_MOL,    YIELD_MO   
      REAL*8                 :: YIELD_RESP,  ACET_MOL,    ACET_MB    
      REAL*8                 :: ACET_MO,     DTSRCE,      AREA_CM2 
      REAL*8                 :: YIELD_ISOP,  MB_MOL,      ACET_RESP 
      REAL*8                 :: ACET_C,      ACET_ISOP,   YIELD_GRASS 
      REAL*8                 :: ACET_GRASS,  ACETSCAL,    AREA_M2
      REAL*8                 :: X,           Y

      ! Scale factors for a posteriori
      REAL*8                 :: MONO_SCALE, DIRECT_SCALE, MB_SCALE 
      REAL*8                 :: DP_SCALE,   GRASS_SCALE

      !=================================================================
      ! EMISS_BIOACET begins here!
      !
      ! The yield for acetone from biogenic sources comes from 
      ! experimental yields from monoterpenes, (Reissell et. al. 1999), 
      ! methyl butenol (Alvarado et. al. 1999) and emissions for 
      ! monoterpenes, methyl butenol, and acetone from Guenther et. al. 
      ! 1999.  Guenther's emissions are for North America, and have 
      ! been scaled to the entire globe
      !=================================================================

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Grid box areas in [m2] and [cm2]
      AREA_M2  = GET_AREA_M2( J )
      AREA_CM2 = AREA_M2 * 1d4

      !=================================================================
      ! (1) BIOGENIC EMISSIONS OF ACETONE FROM MONOTERPENES
      !
      ! Monoterpenes has same # molecules/kg of carbon as isoprene
      ! The yield for monoterpenes is .12 mol/mol from Reisell et.al. 
      ! 1999 (this does not includes direct acetone emissions)
      !=================================================================

      ! Convert [atoms C/box/step] to [molec MONOTERPENE/box/step]
      ! There are 10 C atoms per molecule of MONOTERPENE
      EMMO_MOL   = EMMO / 10d0      

      ! Apply yield from monoterpenes to get [molec ACET/box/step]
      YIELD_MO   = 0.116d0
      ACET_MOL   = EMMO_MOL * YIELD_MO

      ! Convert acetoneemissions back into [atoms C/box/step] 
      ACET_MO    = ACET_MOL * 3.d0

      ! Scale to a posteriori source from Jacob et al 2001 (bdf, 9/5/01)
      MONO_SCALE = 0.89d0
      ACET_MO    = ACET_MO * MONO_SCALE

      ! Convert monoterpene yield to [atoms C/box/s] and 
      ! add to the total biogenic acetone emissions
      ACETONE    = ACETONE + ( ACET_MO / DTSRCE )

      ! Diagnostics -- save ACETONE from MONOTERPENES in [atoms C/cm2/s]
      IF ( ND11 > 0 ) THEN
         AD11(I,J,1) = AD11(I,J,1) + 
     &                 ( ACET_MO / ( AREA_CM2 * DTSRCE ) )
      ENDIF
 
      !=================================================================
      ! (2) BIOGENIC ACETONE FROM METHYL BUTENOL -- NORTH AMERICA
      !
      ! Methyl Butenol (a.k.a. MBO) produces acetone with a molar yield 
      ! of 0.6 [Alvarado (1999)].  The biogenic source of MBO is thought 
      ! to be restricted to North America.  According to Guenther (1999) 
      ! North america emits 3.2Tg-C of MBO, producing 1.15 Tg-C of 
      ! Acetone in North America.
      !=================================================================
      ACET_MB = 0D0

      ! Lon and lat of grid box (I,J) in degrees
      X = GET_XMID( I )
      Y = GET_YMID( J )

      ! Methyl butenol is emitted only in North America, where
      ! ( -167.5 <= lon <= -52.5 ) and ( 16.0 <= lat <= 72.0 ) 
      IF ( ( X >= -167.5d0 .and. X <= -52.5d0 ) .AND. 
     &     ( Y >=   16.0d0 .and. Y <=  72.0d0 ) ) THEN

         ! Convert from [atoms C/box/step] to [molec MBO/box/step] 
         ! There are 5 C atoms per molecule MBO
         MB_MOL   = EMMB / 5.d0      

         ! Apply yield from MBO to get [molec ACET/box/step]
         MB_SCALE = 0.6d0                
         ACET_MOL = MB_MOL * MB_SCALE

         ! Convert from [molec ACET/box/step] to [atoms C/box/step]
         ! There are 3 C atoms per acetone molecule
         ACET_MB  = ACET_MOL * 3.d0      

         ! Scale to a posteriori source from Jacob et al 2001 (bdf, 9/5/01)
         MB_SCALE = 0.76d0 
         ACET_MB  = ACET_MB * MB_SCALE

         ! Convert MBO yield to [atoms C/box/s] and add 
         ! to the total biogenic acetone emissions
         ACETONE  = ACETONE + ( ACET_MB / DTSRCE )
            
         ! Diagnostics -- save ACETONE from MBO in [atoms C/cm2/s]
         IF ( ND11 > 0 ) THEN
            AD11(I,J,2) = AD11(I,J,2) + 
     &                   ( ACET_MB / ( AREA_CM2 * DTSRCE ) )
         ENDIF
      ENDIF

      !=================================================================
      ! (3) BIOGENIC ACETONE -- DIRECT EMISSION 
      ! 
      ! With communication with Singh we have a direct acetone emission 
      ! source of 18 Tg acet/yr that scales to the isoprene emissions.
      !=================================================================

      ! Compute [atoms C/box/step] for ACET, using yield from ISOP
      YIELD_ISOP   = 0.0282d0
      ACET_ISOP    = EMIS * YIELD_ISOP

      ! Scale to a posteriori source from Jacob et al 2001 (bdf, 9/5/01)
      DIRECT_SCALE = 1.06d0
      ACET_ISOP    = ACET_ISOP * DIRECT_SCALE

      ! Convert isoprene yield to [atoms C/box/s] and 
      ! add to the total biogenic acetone emissions
      ACETONE      = ACETONE + ( ACET_ISOP / DTSRCE )

      ! Diagnostics -- save ACETONE from DIRECT EMISSION [atoms C/cm2/s]
      IF ( ND11 > 0 ) THEN
         AD11(I,J,3) = AD11(I,J,3) + 
     &                 ( ACET_ISOP / ( AREA_CM2 * DTSRCE ) )
      ENDIF

      !=================================================================
      ! (4) BIOGENIC ACETONE FROM DRY LEAF MATTER / DEAD PLANTS
      !
      ! According to Warneke et al. 1999, 1 g C of dry leaf matter 
      ! produces at least 10^-4 g C in acetone.  We use this lower limit 
      ! as our g C -> g C yield.  The monthly values of dry leaf matter 
      ! comes from estimates of resp, heterotrophic respiration, from 
      ! Parvada Suntharalingham.  XRESP is in units of g C/m2/s.
      !=================================================================

      ! Convert from [g C in dry leaves/m2/s] to [g C in ACETONE/m2/s]
      YIELD_RESP = 1d-4
      ACET_RESP  = XRESP(I,J) * YIELD_RESP

      ! Scale to a posteriori source from Jacob et al 2001 (bdf, 9/5/01)
      DP_SCALE   = 0.23d0
      ACET_RESP  = ACET_RESP * DP_SCALE

      ! Convert [g C in acetone/m2/s] to [atoms C/box/s] 
      ! and add to the total biogenic ACETONE emissions
      ! The 1000 is for changing [g C] to [kg C]
      ACETONE    = ACETONE + 
     &             ( ACET_RESP * AREA_M2 * XNUMOL_C / 1000d0 )

      ! Diagnostics -- save ACETONE from DRY LEAF MATTER in [atoms C/cm2/s]
      ! the 1000 is for [g C] to [kg C], the 1d4 is for [m-2] to [cm-2]
      IF ( ND11 > 0 ) THEN
         AD11(I,J,4) = AD11(I,J,4) + 
     &                 ( ACET_RESP * XNUMOL_C / ( 1000 * 1D4 ) )
      ENDIF

      !=================================================================
      ! (5) BIOGENIC ACETONE FROM GRASSLANDS 
      !
      ! Direct grass emissions should be about 5 TgC/yr from 
      ! Kirstine et al 1998
      !=================================================================

      ! Compute from [atoms C/box/step] of acetone from grass
      ! for all VOC emissions acetone is about 15%
      YIELD_GRASS = 0.15d0  
      ACET_GRASS  = GRASS * YIELD_GRASS

      ! Scale to a posteriori source from Jacob et al 2001 (bdf, 9/5/01)
      GRASS_SCALE = 1.61d0 
      ACET_GRASS  = ACET_GRASS * GRASS_SCALE

      ! Convert grassland acetone yield to [atoms C/box/s] 
      ! and add to the total biogenic acetone emissions
      ACETONE     = ACETONE + ( ACET_GRASS / DTSRCE )

      ! Diagnostics -- save ACETONE from GRASSLANDS in [atoms C/cm2/s]
      IF ( ND11 > 0 ) THEN
         AD11(I,J,5) = AD11(I,J,5) + 
     &                 ( ACET_GRASS / ( AREA_CM2 * DTSRCE ) )
      ENDIF

      END SUBROUTINE EMISS_BIOACET
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_acetone
!
! !DESCRIPTION: Subroutine CLEANUP\_ACETONE deallocates module arrays
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_ACETONE
! 
! !REVISION HISTORY: 
!  14 Sep 2001 - R. Yantosca - Initial version
!  04 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( ALLOCATED( JO1D  ) ) DEALLOCATE( JO1D  )
      IF ( ALLOCATED( XRESP ) ) DEALLOCATE( XRESP )

      END SUBROUTINE CLEANUP_ACETONE   
!EOC
      END MODULE ACETONE_MOD
