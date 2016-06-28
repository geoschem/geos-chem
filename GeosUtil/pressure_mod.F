!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: pressure_mod
!
! !DESCRIPTION: Module PRESSURE\_MOD contains variables and routines which 
!  specify the grid box pressures for both hybrid or pure-sigma models.  
!  This is necessary for running GEOS-Chem with the GEOS-4 or GEOS-5 hybrid
!  grids.
!\\
!\\
! !INTERFACE: 
!
      MODULE PRESSURE_MOD
! 
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: INIT_PRESSURE
      PUBLIC  :: GET_AP
      PUBLIC  :: GET_BP
      PUBLIC  :: GET_PCENTER              ! wet air P at grid center
      PUBLIC  :: GET_PEDGE                ! wet air P at lower grid edge
      PUBLIC  :: GET_PEDGE_FULLGRID       ! wet air P at lower grid edge
      PUBLIC  :: CLEANUP_PRESSURE
      PUBLIC  :: GET_PEDGE_DRY                ! MOISTURE FIX (ewl, 6/16/16)
      PUBLIC  :: GET_PEDGE_TEST               ! MOISTURE FIX (ewl, 6/16/16)
      PUBLIC  :: SET_FLOATING_PRESSURES       ! MOISTURE FIX (ewl, 5/4/16)
!      PUBLIC  :: GET_PFLT_DRY                 ! MOISTURE FIX (ewl, 5/4/16)

#if defined( ESMF_ )
      PUBLIC  :: Accept_External_Pedge
#endif
!
! !REMARKS:
!
!  Hybrid Grid Coordinate Definition: (dsa, bmy, 8/27/02, 2/2/12)
!  ============================================================================
!                                                                             .
!  GEOS-4, GEOS-5, GEOS-5.7, and MERRA (hybrid grids):
!  ----------------------------------------------------------------------------
!  For GEOS-4/GEOS-5/MERRA met data products, the pressure at the bottom edge 
!  of grid box (I,J,L) is defined as follows:
!                                                                             .
!     Pedge(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
!                                                                             .
!  where
!                                                                             .
!     Psurface(I,J) is  the "true" surface pressure at lon,lat (I,J)
!     Ap(L)         has the same units as surface pressure [hPa]
!     Bp(L)         is  a unitless constant given at level edges
!                                                                             .
!  Ap(L) and Bp(L) are given to us by GMAO.
!                                                                             .
!                                                                             .
!  GEOS-3 (pure-sigma) and GCAP (hybrid grid):
!  ----------------------------------------------------------------------------
!  GEOS-3 is a pure-sigma grid.  GCAP is a hybrid grid, but its grid is
!  defined as if it were a pure sigma grid (i.e. PTOP=150 hPa, and negative
!  sigma edges at higher levels).  For these grids, can stil use the same
!  formula as for GEOS-4/GEOS-5/MERRA, with one modification:
!                                                                             .
!     Pedge(I,J,L) = Ap(L) + [ Bp(L) * ( Psurface(I,J) - PTOP ) ]
!                                                                             .
!  where
!                                                                             .
!     Psurface(I,J) = the "true" surface pressure at lon,lat (I,J)
!     Ap(L)         = PTOP    = model top pressure
!     Bp(L)         = SIGE(L) = bottom sigma edge of level L
!                                                                             .
!                                                                             .
!  The following are true for GCAP, GEOS-3, GEOS-4, GEOS-5, MERRA:
!  ----------------------------------------------------------------------------
!  (1) Bp(LLPAR+1) = 0.0          (L=LLPAR+1 is the atmosphere top)
!  (2) Bp(1)       = 1.0          (L=1       is the surface       )
!  (3) PTOP        = Ap(LLPAR+1)  (L=LLPAR+1 is the atmosphere top) 
!
! !REVISION HISTORY:
!  27 Aug 2002 - D. Abbot & R. Yantosca - Initial version 
!  (1 ) Be sure to check PFLT for NaN or Infinities (bmy, 8/27/02)
!  (2 ) Updated comments (bmy, 5/8/03)
!  (3 ) Updated format string for fvDAS (bmy, 6/19/03)
!  (4 ) Bug fix: use PFLT instead of PFLT-PTOP for GEOS-4 (bmy, 10/24/03)
!  (5 ) Modifications for 30L and 55L GEOS-4 grids (bmy, 11/3/03)
!  (6 ) Added parallel DO-loop in SET_FLOATING_PRESSURE (bmy, 4/14/04)
!  (7 ) Modified for GCAP and GEOS-5 grids (swu, bmy, 5/24/05)
!  (8 ) Removed obsolete reference to "CMN" (bmy, 4/25/06)
!  (9 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (10) Added Ap and Bp for GEOS-5 met fields (bmy, 10/30/07)
!  20 Nov 2009 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Added modifications for MERRA met fields
!  30 Aug 2010 - R. Yantosca - Updated comments
!  02 Feb 2012 - R. Yantosca - Added modifications for GEOS-5.7.x met fields
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  31 Jul 2012 - R. Yantosca - Modifications for grid-independence
!  10 Aug 2012 - R. Yantosca - Remove DEVEL from #ifdef for EXTERNAL_PEDGE
!  11 Dec 2012 - R. Yantosca - Now make EXTERNAL_PEDGE private
!  11 Dec 2012 - R. Yantosca - Add new routine ACCEPT_PEDGE_FROM_ESMF to set
!                              EXTERNAL_PEDGE from the ESMF environment
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  18 Sep 2013 - M. Long     - Now use #if defined( ESMF_ ) for HPC code
!  02 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  11 Aug 2015 - R. Yantosca - Added support for MERRA2 data
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Module variables
      REAL(fp), ALLOCATABLE :: AP(:)                  ! "A" term for hybrid grid
      REAL(fp), ALLOCATABLE :: BP(:)                  ! "B" term for hybrid grid
! Prior to 04 May 2016: old variable PFLT includes moisture
!      REAL(fp), ALLOCATABLE :: PFLT(:,:)              ! "Floating" sfc pressure
      REAL(fp), ALLOCATABLE :: PFLT_DRY(:,:)          ! "Floating" dry sfc pres
      REAL(fp), ALLOCATABLE :: PFLT_WET(:,:)          ! "Floating" wet sfc pres
      REAL(fp), ALLOCATABLE :: AP_FULLGRID(:)         ! "A" term for full grid
      REAL(fp), ALLOCATABLE :: BP_FULLGRID(:)         ! "B" term for full grid
#if defined( ESMF_ )
      REAL(fp), ALLOCATABLE :: EXTERNAL_PEDGE(:,:,:)  ! Pressure edges from 
                                                      !  external grid
#endif
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ap
!
! !DESCRIPTION: Function GET\_AP returns the "A" term [hPa] for the 
!  hybrid ETA coordinate.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_AP( L ) RESULT( AP_TEMP )
!
! !USES:
!
      USE CMN_SIZE_MOD                ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: L        ! GEOS-Chem level index
!
! !RETURN VALUE: 
!
      REAL(fp)              :: AP_TEMP  ! Corresponding "A" value [hPa]
                                      !  at bottom edge of level L
!
! !REVISION HISTORY:
!  20 Aug 2002 - D. Abbot & R. Yantosca - Initial version  
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      AP_TEMP = AP(L)

      END FUNCTION GET_AP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_bp
!
! !DESCRIPTION: Function GET\_BP returns the "B" term [unitless] for the 
!  hybrid ETA coordinate.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_BP( L ) RESULT( BP_TEMP )
!
! !USES:
!
      USE CMN_SIZE_MOD                ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: L        ! GEOS-Chem level index
!
! !RETURN VALUE: 
!
      REAL(fp)              :: BP_TEMP  ! Corresponding "B" value [unitless]
                                      !  at bottom edge of level L
!
! !REVISION HISTORY:
!  20 Aug 2002 - D. Abbot & R. Yantosca - Initial version  
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      BP_TEMP = BP(L)

      END FUNCTION GET_BP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_floating_pressures
!
! !DESCRIPTION: Subroutine SET\_FLOATING\_PRESSURE initializes the 
!  dry floating pressure field PFLT_DRY with the State_Met field
!  PSC2_DRY (interpolated dry surface pressure).  It then calculates 
!  the wet floating pressure from the dry floating pressure by 
!  incorporating the current state specific humidity.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_FLOATING_PRESSURES( am_I_Root, State_Met, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD  ! Size parameters
      USE ERROR_MOD, ONLY : CHECK_VALUE
      USE GIGC_ErrCode_Mod
      USE GIGC_State_Met_Mod, ONLY : MetState
!     
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on root CPU?
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology state object
!
! !OUTPUT ARGUMENTS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Jun 2016 - E. Lundgren- Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER            :: I, J, L
      INTEGER            :: ERR_LOC(4)
      CHARACTER(LEN=255) :: ERR_VAR
      CHARACTER(LEN=255) :: ERR_MSG
      REAL(fp)           :: PEDGE1, PEDGE2, SPHU_KGKG

      !=================================================================
      ! SET_FLOATING_PRESSURES begins here!
      !=================================================================

      ! DEBUGGING (ewl)
      PRINT *, " "
      PRINT *, "In SET_FLOATING_PRESSURES"
      PRINT *, "   Old PFLT_DRY(56,20): ", PFLT_DRY(56,20)
      PRINT *, "   Old PFLT_WET(56,20): ", PFLT_WET(56,20)
      ! END DEBUGGING

      ! Set PFLT_DRY equal to input value PS
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, ERR_LOC, ERR_VAR, ERR_MSG )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! DEBUGGING (ewl)
         ! NOTE: If using gcap, need to add PTOP - come back to this later
         ! END DEBUGGING

         ! Set the floating dry pressure, PFLT_DRY, to the 
         ! interpolate dry pressure, PSC2_DRY
         PFLT_DRY(I,J) = State_Met%PSC2_DRY(I,J)

         ! Check for NaN or Infinities
         ERR_LOC = (/ I, J, 0, 0 /)
         ERR_VAR = 'PFLT_DRY'
         ERR_MSG = 'set_floating_pressures:1'
         CALL CHECK_VALUE( PFLT_DRY(I,J), ERR_LOC, ERR_VAR, ERR_MSG )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! for the wet pressure, use the old way and set PFLT_WET to PSC2 
      ! (the interpolated GMAO pressure). The new way is under 
      ! construction

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, ERR_LOC, ERR_VAR, ERR_MSG )
      DO J = 1, JJPAR
      DO I = 1, IIPAR    
         PFLT_WET(I,J) = State_Met%PSC2(I,J)

         ! Check for NaN or Infinities
         ERR_LOC = (/ I, J, 0, 0 /)
         ERR_VAR = 'PFLT_WET'
         ERR_MSG = 'set_floating_pressures:1'
         CALL CHECK_VALUE( PFLT_WET(I,J), ERR_LOC, ERR_VAR, ERR_MSG )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

!      ! New way: calculate PFLT_WET using PFLT_DRY and specific humidity
!      PFLT_WET = 0.e+0_fp
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L, PEDGE1, PEDGE2, SPHU_KgKg )
!!$OMP+PRIVATE( ERR_LOC, ERR_VAR, ERR_MSG )
!      DO L = 1, LLPAR
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Convert specific humidity from g/kg to kg/kg locally
!         SPHU_KgKg = State_Met%SPHU(I,J,L) * 1.0e-3_fp
!
!         ! Calculate top and bottom dry edge pressures for this 
!         ! level by applying AP and BP params to dry surface pressure
!         PEDGE1    = GET_AP(L)   + GET_BP(L)   * PFLT_DRY(I,J)
!         PEDGE2    = GET_AP(L+1) + GET_BP(L+1) * PFLT_DRY(I,J)
!
!         ! Add the wet delta pressure for this level to the
!         ! cumulative wet pressure in this column
!         PFLT_WET(I,J) = PFLT_WET(I,J) + ( PEDGE2 - PEDGE1 ) / 
!     &                   ( 1.e+0_fp - SPHU_KgKg )
!      ENDDO
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO

      ! NOTE: To accurately get the wet pressure back from the dry
      !       pressure using this method, must also add in a PTOP
      !       bias factor...this is under construction

      ! DEBUGGING (ewl)
      PRINT *, "   New PFLT_DRY(56,20): ", PFLT_DRY(56,20)
      PRINT *, "   New PFLT_WET(56,20): ", PFLT_WET(56,20)
      ! END DEBUGGING

      ! Return to calling program
      END SUBROUTINE SET_FLOATING_PRESSURES
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_pedge
!
! !DESCRIPTION: Function GET\_PEDGE returns the pressure at the bottom edge 
!  of level L.  L=1 is the surface, L=LLPAR+1 is the atm top.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_PEDGE( I, J, L ) RESULT( PEDGE )
!
! !USES:
!
      USE CMN_SIZE_MOD   ! PTOP
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I        ! GEOS-Chem lon   index
      INTEGER, INTENT(IN)   :: J        ! GEOS-Chem lat   index
      INTEGER, INTENT(IN)   :: L        ! GEOS-Chem level index
!
! !RETURN VALUE:
!
      REAL(f8)              :: PEDGE  ! Pressure @ bottom edge of (I,J,L) [hPa]
!
! !REVISION HISTORY:
!  20 Aug 2002 - D. Abbot & R. Yantosca - Initial version  
!  (1 ) Bug fix: use PFLT instead of PFLT-PTOP for GEOS-4 (bmy, 10/24/03)
!  (2 ) Now treat GEOS-5 the same way as GEOS-4 (bmy, 10/30/07)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  13 Aug 2010 - R. Yantosca - Compute PEDGE for MERRA the same as for GEOS-5
!  02 Feb 2012 - R. Yantosca - Compute PEDGE for GEOS-5.7.2 the same as MERRA
!  10 Aug 2012 - R. Yantosca - Need to put #ifdef for EXTERNAL_PEDGE in the
!                              section for GEOS-4, GEOS-5, MERRA, GEOS-5.7.x
!  10 Aug 2012 - R. Yantosca - Now only use Cpp switches EXTERNAL_GRID or 
!                              EXTERNAL_FORCING to use the GCM pressures.
!                              This prevents problems when compiling G-C with
!                              the DEVEL tag when using traditional main.F.
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  23 Dec 2014 - M. Yannetti - Changed output to REAL(f8)
!  11 Aug 2015 - R. Yantosca - Compute PEDGE for MERRA2 the same as for GEOS-FP
!  04 May 2016 - E. Lundgren - Replace PFLT with new variable name PFLT_WET
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( GEOS_4 ) || defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )

      !-------------------------------------------------------
      ! GEOS-4, GEOS-5, GEOS-FP, MERRA, and MERRA-2 met fields
      !-------------------------------------------------------

#if defined( ESMF_ )
      ! Pressure [hPa] at bottom edge of level L (see documentation header)
      ! Taken from the GCM fields
      PEDGE = EXTERNAL_PEDGE(I,J,L)
#else
      ! MOISTURE FIX (ewl, 5/4/16) - Use new variable name PFLT_WET
      PEDGE = AP(L) + ( BP(L) * PFLT_WET(I,J) )

      ! Prior to 04 May 2016: Use old variable name PFLT
      !! Pressure [hPa] at bottom edge of level L (see documentation header)
      !! Computed for use w/in GEOS-Chem
      !PEDGE = AP(L) + ( BP(L) * PFLT(I,J) )
#endif

#else 

      !-----------------------------
      ! GEOS-3 & GCAP met fields
      !-----------------------------

#if defined( ESMF_ )
      ! Pressure [hPa] at bottom edge of level L (see documentation header)
      ! Taken from the GCM fields
      PEDGE = EXTERNAL_PEDGE(I,J,L)
#else
      ! MOISTURE FIX (ewl, 5/4/16) - Use new variable name PFLT_WET
      PEDGE = AP(L) + ( BP(L) * ( PFLT_WET(I,J) - PTOP ) ) 

      ! Prior to 04 May 2016: Use old variable name PFLT
      !! Pressure [hPa] at bottom edge of level L (see documentation header)
      !! Computed for use within GEOS-Chem
      !PEDGE = AP(L) + ( BP(L) * ( PFLT(I,J) - PTOP ) )
#endif

#endif     

      END FUNCTION GET_PEDGE 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_pedge_test
!
! !DESCRIPTION: Function GET\_PEDGE\_TEST returns the pressure at the 
!  bottom edge of level L by calculating it from the dry surface pressure
!  and state humidity. This function is a test to allow passing State_Met
!  to pedge in dao_mod only as part of the moisture fix debugging.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_PEDGE_TEST( II, JJ, LL, State_Met ) RESULT( PEDGE )
!
! !USES:
!
      USE CMN_SIZE_MOD   ! PTOP
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)        :: II        ! GEOS-Chem lon   index
      INTEGER, INTENT(IN)        :: JJ        ! GEOS-Chem lat   index
      INTEGER, INTENT(IN)        :: LL        ! GEOS-Chem level index
      TYPE(MetState), INTENT(IN) :: State_Met ! Meteorology state object
!
! !RETURN VALUE:
!
      INTEGER               :: L, P1, P2
      REAL(f8)              :: PEDGE  ! Pressure @ bottom edge of (I,J,L) [hPa]
!
! !REVISION HISTORY:
!  21 Jun 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      ! NOTE: This routine is for testing only!!!

      ! Calculate PEDGE from PFLT_DRY and SPHU
      PEDGE = 0.e+0_fp
      DO L = LL, LLPAR
        P1 = GET_AP(L) + GET_BP(L) * PFLT_DRY(II,JJ)
        P2 = GET_AP(L+1) + GET_BP(L+1) * PFLT_DRY(II,JJ)
        PEDGE = PEDGE + ( P1 - P2 ) / 
     &                ( 1.e+0_fp - State_Met%SPHU(II,JJ,L) * 1.0e-3_fp )
      ENDDO   

      END FUNCTION GET_PEDGE_TEST 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_pedge_dry
!
! !DESCRIPTION: Function GET\_PEDGE\_DRY returns the pressure at the 
!  bottom edge of level L, reconstructed using the dry surface pressure. 
!  L=1 is the surface, L=LLPAR+1 is the atm top.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_PEDGE_DRY( I, J, L ) RESULT( PEDGE_DRY )
!
! !USES:
!
      USE CMN_SIZE_MOD   ! PTOP
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I      ! GEOS-Chem lon   index
      INTEGER, INTENT(IN) :: J      ! GEOS-Chem lat   index
      INTEGER, INTENT(IN) :: L      ! GEOS-Chem level index
!
! !RETURN VALUE:
!
      REAL(f8)              :: PEDGE_DRY  ! Prssr @ bottom edge of (I,J,L) [hPa]
!
! !REVISION HISTORY:
!  16 Jun 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

! MOISTURE FIX NOTE: Will this break GCHP? Took out external grid...
! Should be okay as long as this is used only in transport.
! But... it is in AIRQNT and interp, etc... Might need to bracket
! those with if transport is on. 

#if   defined( GEOS_4 ) || defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )

      !-------------------------------------------------------
      ! GEOS-4, GEOS-5, GEOS-FP, MERRA, and MERRA-2 met fields
      !-------------------------------------------------------
      PEDGE_DRY = AP(L) + ( BP(L) * PFLT_DRY(I,J) )

#else 

      !-----------------------------
      ! GEOS-3 & GCAP met fields
      !-----------------------------
      PEDGE_DRY = AP(L) + ( BP(L) * ( PFLT_DRY(I,J) - PTOP ) ) 
#endif    

      END FUNCTION GET_PEDGE_DRY 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_pedge_fullgrid
!
! !DESCRIPTION: Function GET\_PEDGE\_FULLGRID returns the pressure at the
!  bottom edge of level L of the unreduced vertical grid.  L=1 is the surface, 
!  L=LLLPAR+1 is the atm top.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_PEDGE_FULLGRID( I, J, L ) RESULT( PEDGE )
!
! !USES:
!
      USE CMN_SIZE_MOD   ! PTOP
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I      ! GEOS-Chem lon   index
      INTEGER, INTENT(IN) :: J      ! GEOS-Chem lat   index
      INTEGER, INTENT(IN) :: L      ! GEOS-Chem level index
!
! !RETURN VALUE:
!
      REAL(fp)              :: PEDGE  ! Pressure @ bottom edge of (I,J,L) [hPa]
!
! !REVISION HISTORY:
!  (1 ) Modified from GET_PEDGE (cdh, 1/22/09)
!  02 Feb 2012 - R. Yantosca - Compute PEDGE for GEOS-5.7.2 the same as MERRA
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  11 Aug 2015 - R. Yantosca - Compute PEDGE for MERRA2 the same as for GEOS-FP
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_PEDGE_FULLGRID begins here!
      !=================================================================

#if   defined( GEOS_4 ) || defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )

      !-----------------------------
      ! GEOS-4 & GEOS-5 met fields
      !-----------------------------

      ! MOISTURE FIX (ewl, 5/4/16) - Use new variable name PFLT_WET
      PEDGE = AP_FULLGRID(L) + ( BP_FULLGRID(L) * PFLT_WET(I,J) )

      ! Prior to 04 May 2016: Use old variable name PFLT
      !! Pressure [hPa] at bottom edge of level L (see documentation header)
      !PEDGE = AP_FULLGRID(L) + ( BP_FULLGRID(L) * PFLT(I,J) )

#else 

      !-----------------------------
      ! GEOS-3 & GCAP met fields
      !-----------------------------

      ! MOISTURE FIX (ewl, 5/4/16) - Use new variable name PFLT_WET
      PEDGE = AP_FULLGRID(L) + ( BP_FULLGRID(L) * 
     &        ( PFLT_WET(I,J) - PTOP ) ) 

      ! Prior to 04 May 2016: Use old variable name PFLT
      !! Pressure [hPa] at bottom edge of level L (see documentation header)
      !PEDGE = AP_FULLGRID(L) + ( BP_FULLGRID(L) * ( PFLT(I,J) - PTOP ) )

#endif     
 
      ! Return to calling program
      END FUNCTION GET_PEDGE_FULLGRID
!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: get_pflt_dry
!!
!! !DESCRIPTION: Function GET\_PFLT\_DRY returns the dry surface pressure at 
!!   (lat,lon) grid point (I,J).
!!\\
!!\\
!! !INTERFACE:
!!
!      FUNCTION GET_PFLT_DRY( I, J ) RESULT( PS )
!!
!! !INPUT PARAMETERS: 
!!
!      INTEGER, INTENT(IN) :: I        ! GEOS-Chem lon   index
!      INTEGER, INTENT(IN) :: J        ! GEOS-Chem lat   index
!!
!! !RETURN VALUE:
!!
!      REAL(fp)            :: PS  ! Floating dry pressure @ (I,J,1) [hPa]
!!
!! !REVISION HISTORY:
!!  04 May 2016 - E. Lundgren - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!      PS = PFLT_DRY(I,J)
!
!      END FUNCTION GET_PFLT_DRY
!!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_pcenter
!
! !DESCRIPTION: Function GET\_PCENTER returns the pressure at the vertical
!  midpoint of level L.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_PCENTER( I, J, L ) RESULT( PCENTER )
!
! !USES:
!
      USE CMN_SIZE_MOD   ! PTOP
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I        ! GEOS-Chem lon   index
      INTEGER, INTENT(IN) :: J        ! GEOS-Chem lat   index
      INTEGER, INTENT(IN) :: L        ! GEOS-Chem level index
!
! !RETURN VALUE:
!
      REAL(fp)              :: PCENTER  ! Pressure @ center of (I,J,L) [hPa]
!
! !REVISION HISTORY:
!  20 Aug 2002 - D. Abbot & R. Yantosca - Initial version  
!  (1 ) Updated format string for fvDAS (bmy, 6/19/03)
!  (2 ) Removed reference to "CMN", it's obsolete (bmy, 4/25/06)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC

      ! The pressure at the center of a grid-box is found
      ! by averaging the pressures at the box's two edges
      PCENTER = 0.5e+0_fp * ( GET_PEDGE(I,J,L) + GET_PEDGE(I,J,L+1) )

      END FUNCTION GET_PCENTER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pressure
!
! !DESCRIPTION: Subroutine INIT\_PRESSURE allocates and initializes the AP 
!  and BP arrays.  It must be called in "main.f", after SIGE is defined.  
!  GEOS-4 and GEOS-5 requires the hybrid pressure system specified by 
!  the listed values of AP and BP, while earlier versions of GEOS use a pure 
!  sigma pressure system.  GCAP met fields (based on GISS) also use a hybrid 
!  system. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_PRESSURE( am_I_Root )
!
! !USES:
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD  ! LLPAR, PTOP
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !REVISION HISTORY:
!  27 Aug 2002 - D. Abbot, S. Wu, & R. Yantosca - Initial version 
!  (1 ) Now reference ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Now echo Ap, Bp to std output (bmy, 3/14/03)
!  (3 ) Now print LLPAR+1 levels for Ap, Bp.  Remove reference to SIGE, it's
!        obsolete.  Also now use C-preprocessor switch GRID30LEV instead of
!        IF statements to define vertical coordinates. (bmy, 11/3/03)
!  (4 ) Now modified for both GCAP & GEOS-5 vertical grids (swu, bmy, 5/24/05)
!  (5 ) Renamed GRID30LEV to GRIDREDUCED (bmy, 10/30/07)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  13 Aug 2010 - R. Yantosca - Compute Ap and Bp for MERRA the same way as for
!                              GEOS-5.  The vertical grids are identical.
!  30 Aug 2010 - R. Yantosca - Updated comments
!  30 Nov 2010 - R. Yantosca - Further improved comments about how GEOS-4 and
!                              GEOS-5 vertical levels are lumped together.\
!  02 Feb 2012 - R. Yantosca - Compute Ap and Bp for GEOS-5.7.x in the same way
!                              as for GEOS-5 and MERRA (grids are identical)
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  11 Aug 2015 - R. Yantosca - Init MERRA2 Ap & Bp the same way as for GEOS-FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS
      INTEGER :: L

      ! Temporary array for GCAP met fields (swu, bmy, 5/24/05)
      REAL(fp)  :: SIGE_GCAP(LLPAR+1)  

      !=================================================================
      ! INIT_PRESSURE begins here!
      !=================================================================
! Prior to 04 May 2016:
!      ALLOCATE( PFLT( IIPAR, JJPAR ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PFLT' )
!      PFLT = 0e+0_fp

      ALLOCATE( PFLT_DRY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PFLT_DRY' )
      PFLT_DRY = 0e+0_fp

      ALLOCATE( PFLT_WET( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PFLT_WET' )
      PFLT_WET = 0e+0_fp

      ALLOCATE( AP( LLPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AP' )
      AP = 1e+0_fp

      ALLOCATE( BP( LLPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BP' )
      BP = 0e+0_fp

      ALLOCATE( AP_FULLGRID( LGLOB+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AP_FULLGRID' )
      AP = 1e+0_fp

      ALLOCATE( BP_FULLGRID( LGLOB+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BP_FULLGRID' )
      BP = 0e+0_fp

#if defined( ESMF_ )
      ALLOCATE( EXTERNAL_PEDGE( IIPAR, JJPAR, LLPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EXTERNAL_PEDGE' )
      EXTERNAL_PEDGE = 0e+0_fp
#endif

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )

      !=================================================================
      ! GEOS-5/MERRA vertical coordinates (47 or 72 levels)
      !=================================================================

#if   defined( GRIDREDUCED )

      !-----------------------------------------------------------------
      ! GEOS-5/MERRA/GEOS-FP/MERRA2 47-level reduced vertical grid
      !  
      !  Bottom   Bottom    # levels
      !  edge of  edge prs  lumped 
      !  level    (hPa)     together
      !
      !   PTOP       0.010   
      !    47        0.066     4
      !    46        0.211     4
      !    45        0.617     4
      !    44        1.651     4
      !    43        4.077     4
      !    42        9.293     4
      !    41       19.792     4
      !    40       28.368     2
      !    39       40.175     2
      !    38       56.388     2
      !    37       78.512     2
      ! %%%% START LUMPING LEVELS ABOVE HERE %%%%%
      !    36       92.366       
      !    35      108.663
      !    34      127.837
      !    33      150.393
      !    32      176.930
      ! %%%% FIXED-PRESSURE LEVELS BEGIN HERE %%%%
      !-----------------------------------------------------------------

      ! Ap [hPa] for 47 levels (48 edges)
      AP = (/ 0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01,
     &        1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01,
     &        4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01,
     &        7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01,
     &        1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02,
     &        1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02,
     &        2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02,
     &        2.243630d+02, 2.168650d+02, 2.011920d+02, 1.769300d+02,
     &        1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01,
     &        7.851231d+01, 5.638791d+01, 4.017541d+01, 2.836781d+01, 
     &        1.979160d+01, 9.292942d+00, 4.076571d+00, 1.650790d+00, 
     &        6.167791d-01, 2.113490d-01, 6.600001d-02, 1.000000d-02 /)

      ! Bp [unitless] for 47 levels (48 edges)
      BP = (/ 1.000000d+00, 9.849520d-01, 9.634060d-01, 9.418650d-01,
     &        9.203870d-01, 8.989080d-01, 8.774290d-01, 8.560180d-01,
     &        8.346609d-01, 8.133039d-01, 7.919469d-01, 7.706375d-01,
     &        7.493782d-01, 7.211660d-01, 6.858999d-01, 6.506349d-01,
     &        6.158184d-01, 5.810415d-01, 5.463042d-01, 4.945902d-01,
     &        4.437402d-01, 3.928911d-01, 3.433811d-01, 2.944031d-01,
     &        2.467411d-01, 2.003501d-01, 1.562241d-01, 1.136021d-01,
     &        6.372006d-02, 2.801004d-02, 6.960025d-03, 8.175413d-09,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00 /)

      !--------------------------------
      ! GEOS-5 72 level grid
      !--------------------------------

      ! Ap [hPa] for 72 levels (73 edges)
      AP_FULLGRID =
     &     (/ 0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01,
     &        1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01,
     &        4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01,
     &        7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01,
     &        1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02,
     &        1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02,
     &        2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02,
     &        2.243630d+02, 2.168650d+02, 2.011920d+02, 1.769300d+02,
     &        1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01,
     &        7.851231d+01, 6.660341d+01, 5.638791d+01, 4.764391d+01,
     &        4.017541d+01, 3.381001d+01, 2.836781d+01, 2.373041d+01,
     &        1.979160d+01, 1.645710d+01, 1.364340d+01, 1.127690d+01,
     &        9.292942d+00, 7.619842d+00, 6.216801d+00, 5.046801d+00,
     &        4.076571d+00, 3.276431d+00, 2.620211d+00, 2.084970d+00,
     &        1.650790d+00, 1.300510d+00, 1.019440d+00, 7.951341d-01,
     &        6.167791d-01, 4.758061d-01, 3.650411d-01, 2.785261d-01,
     &        2.113490d-01, 1.594950d-01, 1.197030d-01, 8.934502d-02,
     &        6.600001d-02, 4.758501d-02, 3.270000d-02, 2.000000d-02,
     &        1.000000d-02 /)

      ! Bp [unitless] for 72 levels (73 edges)
      BP_FULLGRID =
     &     (/ 1.000000d+00, 9.849520d-01, 9.634060d-01, 9.418650d-01,
     &        9.203870d-01, 8.989080d-01, 8.774290d-01, 8.560180d-01,
     &        8.346609d-01, 8.133039d-01, 7.919469d-01, 7.706375d-01,
     &        7.493782d-01, 7.211660d-01, 6.858999d-01, 6.506349d-01,
     &        6.158184d-01, 5.810415d-01, 5.463042d-01, 4.945902d-01,
     &        4.437402d-01, 3.928911d-01, 3.433811d-01, 2.944031d-01,
     &        2.467411d-01, 2.003501d-01, 1.562241d-01, 1.136021d-01,
     &        6.372006d-02, 2.801004d-02, 6.960025d-03, 8.175413d-09,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00 /)

#else

      !-----------------------------------------------------------------
      ! GEOS-5 72 level grid
      !-----------------------------------------------------------------

      ! Ap [hPa] for 72 levels (73 edges)
      AP = (/ 0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01,
     &        1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01,
     &        4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01,
     &        7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01,
     &        1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02,
     &        1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02,
     &        2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02,
     &        2.243630d+02, 2.168650d+02, 2.011920d+02, 1.769300d+02,
     &        1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01,
     &        7.851231d+01, 6.660341d+01, 5.638791d+01, 4.764391d+01,
     &        4.017541d+01, 3.381001d+01, 2.836781d+01, 2.373041d+01,
     &        1.979160d+01, 1.645710d+01, 1.364340d+01, 1.127690d+01,
     &        9.292942d+00, 7.619842d+00, 6.216801d+00, 5.046801d+00,
     &        4.076571d+00, 3.276431d+00, 2.620211d+00, 2.084970d+00,
     &        1.650790d+00, 1.300510d+00, 1.019440d+00, 7.951341d-01,
     &        6.167791d-01, 4.758061d-01, 3.650411d-01, 2.785261d-01,
     &        2.113490d-01, 1.594950d-01, 1.197030d-01, 8.934502d-02,
     &        6.600001d-02, 4.758501d-02, 3.270000d-02, 2.000000d-02,
     &        1.000000d-02 /)

      ! Bp [unitless] for 72 levels (73 edges)
      BP = (/ 1.000000d+00, 9.849520d-01, 9.634060d-01, 9.418650d-01,
     &        9.203870d-01, 8.989080d-01, 8.774290d-01, 8.560180d-01,
     &        8.346609d-01, 8.133039d-01, 7.919469d-01, 7.706375d-01,
     &        7.493782d-01, 7.211660d-01, 6.858999d-01, 6.506349d-01,
     &        6.158184d-01, 5.810415d-01, 5.463042d-01, 4.945902d-01,
     &        4.437402d-01, 3.928911d-01, 3.433811d-01, 2.944031d-01,
     &        2.467411d-01, 2.003501d-01, 1.562241d-01, 1.136021d-01,
     &        6.372006d-02, 2.801004d-02, 6.960025d-03, 8.175413d-09,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00,
     &        0.000000d+00 /)
      
#endif

#elif defined( GEOS_4 )
      
      !=================================================================
      ! GEOS-4 vertical coordinates (30 or 55 levels)
      !=================================================================

#if   defined( GRIDREDUCED )

      !-----------------------------------------------------------------
      ! GEOS-4 30-level reduced vertical grid
      !  
      !  Bottom    Bottom    # levels
      !  edge of   edge prs  lumped 
      !  level     (hPa)     together
      !  
      !   PTOP       0.010   
      !    30        0.066      4
      !    29        0.211      4
      !    28        0.617      4
      !    27        1.651      4
      !    26        4.077      4
      !    25        9.293      4
      !    24       19.792      4
      !    23       28.368      2
      !    22       40.175      2
      !    21       56.388      2
      !    20       78.512      2
      ! %%%% START LUMPING LEVELS ABOVE HERE %%%%%
      !    19       92.366       
      !    18      108.663
      !    17      127.837
      !    16      150.393
      !    15      176.930
      ! %%%% FIXED-PRESSURE LEVELS BEGIN HERE %%%%
      !-----------------------------------------------------------------

      ! Ap [hPa] for 30 levels (31 edges)
      AP = (/  0.000000e+0_fp,   0.000000e+0_fp,
     &        12.704939e+0_fp,  35.465965e+0_fp, 
     &        66.098427e+0_fp, 101.671654e+0_fp,
     &       138.744400e+0_fp, 173.403183e+0_fp, 
     &       198.737839e+0_fp, 215.417526e+0_fp,
     &       223.884689e+0_fp, 224.362869e+0_fp, 
     &       216.864929e+0_fp, 201.192093e+0_fp,
     &       176.929993e+0_fp, 150.393005e+0_fp, 
     &       127.837006e+0_fp, 108.663429e+0_fp,
     &        92.365662e+0_fp,  78.512299e+0_fp, 
     &        56.387939e+0_fp,  40.175419e+0_fp,
     &        28.367815e+0_fp,  19.791553e+0_fp, 
     &         9.292943e+0_fp,   4.076567e+0_fp,
     &         1.650792e+0_fp,   0.616779e+0_fp, 
     &         0.211349e+0_fp,   0.066000e+0_fp,
     &         0.010000e+0_fp /)

      ! Bp [unitless] for 30 levels (31 edges)
      BP = (/  1.000000e+0_fp,   0.985110e+0_fp,
     &         0.943290e+0_fp,   0.867830e+0_fp, 
     &         0.764920e+0_fp,   0.642710e+0_fp,
     &         0.510460e+0_fp,   0.378440e+0_fp, 
     &         0.270330e+0_fp,   0.183300e+0_fp,
     &         0.115030e+0_fp,   0.063720e+0_fp, 
     &         0.028010e+0_fp,   0.006960e+0_fp,
     &         0.000000e+0_fp,   0.000000e+0_fp, 
     &         0.000000e+0_fp,   0.000000e+0_fp,
     &         0.000000e+0_fp,   0.000000e+0_fp, 
     &         0.000000e+0_fp,   0.000000e+0_fp,
     &         0.000000e+0_fp,   0.000000e+0_fp, 
     &         0.000000e+0_fp,   0.000000e+0_fp,
     &         0.000000e+0_fp,   0.000000e+0_fp, 
     &         0.000000e+0_fp,   0.000000e+0_fp,
     &         0.000000e+0_fp /)

      !-----------------------------------
      ! GEOS-4 55 level grid
      !-----------------------------------

      ! AP [hPa] for 55 levels (56 edges)
      AP_FULLGRID =
     &     (/ 0.000000e+0_fp,   0.000000e+0_fp,
     &       12.704939e+0_fp,  35.465965e+0_fp, 
     &       66.098427e+0_fp, 101.671654e+0_fp,
     &      138.744400e+0_fp, 173.403183e+0_fp,
     &      198.737839e+0_fp, 215.417526e+0_fp,
     &      223.884689e+0_fp, 224.362869e+0_fp,
     &      216.864929e+0_fp, 201.192093e+0_fp,
     &      176.929993e+0_fp, 150.393005e+0_fp,
     &      127.837006e+0_fp, 108.663429e+0_fp,
     &       92.365662e+0_fp,  78.512299e+0_fp, 
     &       66.603378e+0_fp,  56.387939e+0_fp,
     &       47.643932e+0_fp,  40.175419e+0_fp, 
     &       33.809956e+0_fp,  28.367815e+0_fp,
     &       23.730362e+0_fp,  19.791553e+0_fp, 
     &       16.457071e+0_fp,  13.643393e+0_fp,
     &       11.276889e+0_fp,   9.292943e+0_fp,
     &        7.619839e+0_fp,   6.216800e+0_fp,
     &        5.046805e+0_fp,   4.076567e+0_fp, 
     &        3.276433e+0_fp,   2.620212e+0_fp,
     &        2.084972e+0_fp,   1.650792e+0_fp,
     &        1.300508e+0_fp,   1.019442e+0_fp,
     &        0.795134e+0_fp,   0.616779e+0_fp, 
     &        0.475806e+0_fp,   0.365041e+0_fp,
     &        0.278526e+0_fp,   0.211349e+0_fp, 
     &        0.159495e+0_fp,   0.119703e+0_fp,
     &        0.089345e+0_fp,   0.066000e+0_fp, 
     &        0.047585e+0_fp,   0.032700e+0_fp,
     &        0.020000e+0_fp,   0.010000e+0_fp /)

      ! BP [unitless] for 55 levels (56 edges)
      BP_FULLGRID =
     &     (/  1.000000e+0_fp,  0.985110e+0_fp,
     &         0.943290e+0_fp,  0.867830e+0_fp,
     &         0.764920e+0_fp,  0.642710e+0_fp,
     &         0.510460e+0_fp,  0.378440e+0_fp,
     &         0.270330e+0_fp,  0.183300e+0_fp,
     &         0.115030e+0_fp,  0.063720e+0_fp,
     &         0.028010e+0_fp,  0.006960e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp /)

#else

      !-----------------------------------------------------------------
      ! GEOS-4 55 level grid
      !-----------------------------------------------------------------

      ! AP [hPa] for 55 levels (56 edges)
      AP = (/ 0.000000e+0_fp,   0.000000e+0_fp,
     &       12.704939e+0_fp,  35.465965e+0_fp, 
     &       66.098427e+0_fp, 101.671654e+0_fp,
     &      138.744400e+0_fp, 173.403183e+0_fp,
     &      198.737839e+0_fp, 215.417526e+0_fp,
     &      223.884689e+0_fp, 224.362869e+0_fp,
     &      216.864929e+0_fp, 201.192093e+0_fp,
     &      176.929993e+0_fp, 150.393005e+0_fp,
     &      127.837006e+0_fp, 108.663429e+0_fp,
     &       92.365662e+0_fp,  78.512299e+0_fp, 
     &       66.603378e+0_fp,  56.387939e+0_fp,
     &       47.643932e+0_fp,  40.175419e+0_fp, 
     &       33.809956e+0_fp,  28.367815e+0_fp,
     &       23.730362e+0_fp,  19.791553e+0_fp, 
     &       16.457071e+0_fp,  13.643393e+0_fp,
     &       11.276889e+0_fp,   9.292943e+0_fp,
     &        7.619839e+0_fp,   6.216800e+0_fp,
     &        5.046805e+0_fp,   4.076567e+0_fp, 
     &        3.276433e+0_fp,   2.620212e+0_fp,
     &        2.084972e+0_fp,   1.650792e+0_fp,
     &        1.300508e+0_fp,   1.019442e+0_fp,
     &        0.795134e+0_fp,   0.616779e+0_fp, 
     &        0.475806e+0_fp,   0.365041e+0_fp,
     &        0.278526e+0_fp,   0.211349e+0_fp, 
     &        0.159495e+0_fp,   0.119703e+0_fp,
     &        0.089345e+0_fp,   0.066000e+0_fp, 
     &        0.047585e+0_fp,   0.032700e+0_fp,
     &        0.020000e+0_fp,   0.010000e+0_fp /)

      ! BP [unitless] for 55 levels (56 edges)
      BP = (/  1.000000e+0_fp,  0.985110e+0_fp,
     &         0.943290e+0_fp,  0.867830e+0_fp,
     &         0.764920e+0_fp,  0.642710e+0_fp,
     &         0.510460e+0_fp,  0.378440e+0_fp,
     &         0.270330e+0_fp,  0.183300e+0_fp,
     &         0.115030e+0_fp,  0.063720e+0_fp,
     &         0.028010e+0_fp,  0.006960e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp,
     &         0.000000e+0_fp,  0.000000e+0_fp /)


#endif

#elif defined( GCAP )

      !=================================================================
      ! GCAP vertical coordinates (23 levels)
      !=================================================================

      ! Define SIGMA edges from GISS 23L model
      SIGE_GCAP = (/ 1.0e+0_fp,           0.9712230e+0_fp,   
     &               0.9340528e+0_fp,     
     &               0.8800959e+0_fp,     0.8021583e+0_fp,   
     &               0.6714628e+0_fp, 
     &               0.5035971403e+0_fp,  0.3297362030e+0_fp,
     &               0.1966426820e+0_fp,  
     &               0.1139088720e+0_fp,  0.0503597111e+0_fp,
     &               0.0000000000e+0_fp, 
     &              -0.0395683460e+0_fp, -0.0764988065e+0_fp,
     &              -0.1124700233e+0_fp, 
     &              -0.1419664323e+0_fp, -0.1585131884e+0_fp,
     &              -0.1678657085e+0_fp, 
     &              -0.1743045598e+0_fp, -0.1781055182e+0_fp,
     &              -0.1793033630e+0_fp, 
     &              -0.1796822548e+0_fp, -0.1798187047e+0_fp,
     &              -0.1798536479e+0_fp /)

      ! Convert SIGMA to AP and BP coordinates
      DO  L = 1, LLCHEM
         AP(L) = 150e+0_fp + ( SIGE_GCAP(L) * ( PTOP - 150e+0_fp ) )
         BP(L) = SIGE_GCAP(L)
      ENDDO

      DO L = LLCHEM+1, LLPAR+1  
         AP(L) = 150e+0_fp + ( SIGE_GCAP(L) * (984e+0_fp - 150e+0_fp) )
         BP(L) = 0e+0_fp
      ENDDO

      AP_FULLGRID = AP
      BP_FULLGRID = BP

#endif
      
#if !defined( ESMF_ ) 
      ! Echo info to std output (skip if using ESMF interface to GEOS-5 GCM)
      IF ( am_I_Root ) THEN
         WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
         WRITE( 6, '(a,/)' ) 'V E R T I C A L   G R I D   S E T U P'
         WRITE( 6, '(a,/)' ) 'INIT_PRESSURE: Vertical coordinates!'
         WRITE( 6, '( ''Ap '', /, 6(f11.6,1x) )' ) AP(1:LLPAR+1)
         WRITE( 6, '(a)'   )
         WRITE( 6, '( ''Bp '', /, 6(f11.6,1x) )' ) BP(1:LLPAR+1)
         WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      ENDIF
#endif

      END SUBROUTINE INIT_PRESSURE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_pressure
!
! !DESCRIPTION: Subroutine CLEANUP\_PRESSURE deallocates all allocated arrays 
!  at the end of a GEOS-Chem model run.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_PRESSURE
!
! !REVISION HISTORY:
!  20 Aug 2002 - D. Abbot & R. Yantosca - Initial version  
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( ALLOCATED( AP          ) ) DEALLOCATE( AP          )
      IF ( ALLOCATED( BP          ) ) DEALLOCATE( BP          )
      IF ( ALLOCATED( AP_FULLGRID ) ) DEALLOCATE( AP_FULLGRID )
      IF ( ALLOCATED( BP_FULLGRID ) ) DEALLOCATE( BP_FULLGRID )
! Prior to 04 May 2016:
!      IF ( ALLOCATED( PFLT        ) ) DEALLOCATE( PFLT        )
      IF ( ALLOCATED( PFLT_DRY    ) ) DEALLOCATE( PFLT_DRY    )
      IF ( ALLOCATED( PFLT_WET    ) ) DEALLOCATE( PFLT_WET    )
#if defined( ESMF_ )
      IF ( ALLOCATED( EXTERNAL_PEDGE ) ) DEALLOCATE( EXTERNAL_PEDGE )
#endif

      END SUBROUTINE CLEANUP_PRESSURE
!EOC
!#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
#if defined( ESMF_ )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: accept_external_pedge
!
! !DESCRIPTION: Subroutine ACCEPT\_EXTERNAL\_PEDGE sets the GEOS-Chem
!  pressure edge variable with the values obtained from an external GCM
!  (such as the NASA GEOS-5 GCM).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Accept_External_Pedge( am_I_Root, State_Met, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE GIGC_State_Met_Mod, ONLY : MetState
!     
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on root CPU?
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology state object
!
! !OUTPUT ARGUMENTS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine is a setter for EXTERNAL_PEDGE.  It allows us to keep the 
!  EXTERNAL_PEDGE array PRIVATE to this module, which is good programming
!  practice.
! 
! !REVISION HISTORY: 
!  06 Dec 2012 - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Set EXTERNAL_PEDGE to the pressure edges [hPa] carried in the
      ! State_Met object, which were obtained from the external GCM
      EXTERNAL_PEDGE = State_Met%PEDGE

      ! Return successfully
      RC             = GIGC_SUCCESS
      
      END SUBROUTINE Accept_External_Pedge
!EOC
#endif
      END MODULE PRESSURE_MOD

