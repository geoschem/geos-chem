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
!  24 Aug 2017 - M. Sulprizio  - Remove support for GCAP, GEOS-4, GEOS-5 and
!                                MERRA
!  06 Nov 2017 - R. Yantosca   - Remove the fixed-tropopause arrays and related
!                                routines INIT_CHEMGRID, CLEANUP_CHEMGRID
!  06 Nov 2017 - R. Yantosca   - Remove LVARTROP, we will always use the 
!                                dynamic tropopause for GEOS-FP, MERRA-2
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!      
      ! Diagnostic flags
      LOGICAL :: Archive_TropopauseLevel
      LOGICAL :: Archive_TropopauseHeight

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
      SUBROUTINE CHECK_VAR_TROP( am_I_Root, State_Met, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
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
      INTEGER  :: I, J
      REAL(fp) :: TPAUSE_LEV(IIPAR,JJPAR)

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
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         
         ! Exit with error!
         RC = GC_FAILURE
         RETURN
      ENDIF

      END SUBROUTINE CHECK_VAR_TROP
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
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE State_Met_Mod, ONLY : MetState

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
!  06 Nov 2017 - R. Yantosca - Remove fixed-tropopause IF block
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
     &           'FIND T-PAUSE !'
            WRITE( 6, 160   ) L
 160        FORMAT( 'L reaches ', i3 )
            WRITE( 6, '(a)' ) 'STOP in GET_TPAUSE_LEVEL'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            CALL GEOS_CHEM_STOP
         ENDIF

      ENDDO

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
!  06 Nov 2017 - R. Yantosca - Remove fixed-tropopause IF block
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

      ! Get bottom pressure edge
      PRESS_BEDGE = State_Met%PEDGE(I,J,L)

      ! Check against actual tropopause pressure
      IS_TROP     = ( PRESS_BEDGE > State_Met%TROPP(I,J) )

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
      SUBROUTINE DIAG_TROPOPAUSE( am_I_Root, Input_Opt, 
     &                            State_Met, State_Diag, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
#if defined( BPCH_DIAG )
      USE DIAG_MOD,       ONLY : AD55
#endif
      USE ErrCode_Mod
      USE Input_Opt_Mod,  ONLY : OptInput
      USE State_Met_Mod,  ONLY : MetState    
      USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
      INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
!
! !REMARKS:
!  For GEOS-FP and MERRA-2, we use the tropopause pressure from the met 
!  field archive to determine if we are in the tropopause or not.
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
!  06 Nov 2017 - R. Yantosca - Now attach netCDF diagnostics (level, height)
!EOP
!------------------------------------------------------------------------------
!BOC

!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL   :: Do_Diag
      INTEGER   :: I, J, L
      INTEGER   :: L_TP
      REAL(fp)  :: H, FRAC, Pb, Pt

      !=======================================================================
      ! DIAG_TROPOPAUSE begins here!
      !=======================================================================

      ! Initialize
      RC      = GC_SUCCESS
      Do_Diag = ( Input_Opt%ND55 > 0      .or.
     &            Archive_TropopauseLevel .or. 
     &            Archive_TropopauseHeight )

      ! Exit if diagnostic is not turned on
      IF ( .not. Do_Diag ) RETURN

      !=================================================================
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
         !-------------------------------------------------------------
         ! ND55 (bpch) DIAGNOSTIC
         !
         ! Tropopause level, height [km], and pressure [hPa]
         !-------------------------------------------------------------
         AD55(I,J,1) = AD55(I,J,1) + L_TP        ! T'pause level
         AD55(I,J,2) = AD55(I,J,2) + H/1.0e+3_fp ! T'pause height [km]
         AD55(I,J,3) = AD55(I,J,3) +             ! T'pause pressure [hPa]
     &                    State_Met%TROPP(I,J)
#endif


#if defined( NC_DIAG )
         !-------------------------------------------------------------
         ! HISTORY (aka netCDF diagnostics)
         !
         ! Tropopause level and height [km]
         ! Tropopuase pressure is obtained with MET_TROPP
         !-------------------------------------------------------------
         IF ( Archive_TropopauseLevel ) THEN
            State_Diag%TropopauseLevel(I,J) = L_TP
         ENDIF

         IF ( Archive_TropopauseHeight ) THEN
            State_Diag%TropopauseHeight(I,J) = H / 1.0e+3_fp
         ENDIF
#endif

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

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
      SUBROUTINE INIT_CHEMGRID( am_I_Root, Input_Opt, State_Diag, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,  ONLY : OptInput
      USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  10 Feb 2005 - R. Yantosca - Initial version
!  09 Sep 2010 - R. Yantosca - Added ProTeX headers
!  29 Mar 2013 - R. Yantosca - Now made public so we can shadow LVARTROP
!  06 Nov 2017 - R. Yantosca - Remove allocations, initialize diagnostic flags
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

      ! Test if the tropopause level diagnostic is activated
      Archive_TropopauseLevel  = 
     &     ASSOCIATED( State_Diag%TropopauseLevel )

      ! Test if the tropopause height diagnostic is activated
      Archive_TropopauseHeight = 
     &     ASSOCIATED( State_Diag%TropopauseHeight )

      END SUBROUTINE INIT_CHEMGRID
!EOC
      END MODULE CHEMGRID_MOD
