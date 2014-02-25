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
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_PRESSURE
      PUBLIC  :: GET_AP
      PUBLIC  :: GET_BP
      PUBLIC  :: GET_PCENTER
      PUBLIC  :: GET_PEDGE
      PUBLIC  :: GET_PEDGE_FULLGRID
      PUBLIC  :: INIT_PRESSURE
      PUBLIC  :: SET_FLOATING_PRESSURE
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Module variables
      REAL*8, ALLOCATABLE :: AP(:)                  ! "A" term for hybrid grid
      REAL*8, ALLOCATABLE :: BP(:)                  ! "B" term for hybrid grid
      REAL*8, ALLOCATABLE :: PFLT(:,:)              ! "Floating" sfc pressure
      REAL*8, ALLOCATABLE :: AP_FULLGRID(:)         ! "A" term for full grid
      REAL*8, ALLOCATABLE :: BP_FULLGRID(:)         ! "B" term for full grid
#if defined( ESMF_ )
      REAL*8, ALLOCATABLE :: EXTERNAL_PEDGE(:,:,:)  ! Pressure edges from 
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
      REAL*8              :: AP_TEMP  ! Corresponding "A" value [hPa]
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
      REAL*8              :: BP_TEMP  ! Corresponding "B" value [unitless]
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
! !IROUTINE: set_floating_pressure
!
! !DESCRIPTION: Subroutine SET\_FLOATING\_PRESSURE initializes the floating 
!  pressure field PFLT with a pressure from the main program.  This is needed
!  to initialize and reset PFLT after transport.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_FLOATING_PRESSURE( PS )
!
! !USES:
!
      USE ERROR_MOD, ONLY : CHECK_VALUE

      USE CMN_SIZE_MOD  ! Size parameters
!
! !INPUT PARAMETERS: 
!
!------------------------------------------------------------------------------
! Prior to 11/4/13:
! Eliminate array temporary by accepting assumed-shape argument into this
! routine. (mpayer, 11/4/13)
!      ! Array containing pressure with which to initialize PFLT [hPa]
!      REAL*8, INTENT(IN) :: PS(IIPAR,JJPAR)  
!------------------------------------------------------------------------------
      ! Array containing pressure with which to initialize PFLT [hPa]
      REAL*8, INTENT(IN) :: PS(:,:)  
!
! !REVISION HISTORY:
!  27 Aug 2002 - D. Abbot, B. Field,  R. Yantosca - Initial version  
!  (1 ) Now check PFLT for NaN or Infinities (bmy, 8/27/02)
!  (2 ) Added parallel DO-loop (bmy, 4/14/04)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  10 Jun 2013 - R. Yantosca - Avoid array temporaries in CHECK_VALUE
!  04 Nov 2013 - M. Sulprizio- Eliminate array temporary by accepting assumed-
!                              shape dummy argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER            :: I, J
      INTEGER            :: ERR_LOC(4)
      CHARACTER(LEN=255) :: ERR_VAR
      CHARACTER(LEN=255) :: ERR_MSG

      !=================================================================
      ! SET_FLOATING_PRESSURE begins here!
      !=================================================================

      ! Define quantities for CHECK_VALUE
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, ERR_LOC, ERR_VAR, ERR_MSG )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Assign into PFLT array
         PFLT(I,J) = PS(I,J) !  12 Jun 2013 - R. Yantosca -  


         ! Check for NaN or Infinities
         ERR_LOC = (/ I, J, 0, 0 /)
         ERR_VAR = 'PFLT'
         ERR_MSG = 'set_floating_pressure:1'
         CALL CHECK_VALUE( PFLT(I,J), ERR_LOC, ERR_VAR, ERR_MSG )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SET_FLOATING_PRESSURE
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
      INTEGER, INTENT(IN) :: I      ! GEOS-Chem lon   index
      INTEGER, INTENT(IN) :: J      ! GEOS-Chem lat   index
      INTEGER, INTENT(IN) :: L      ! GEOS-Chem level index
!
! !RETURN VALUE:
!
      REAL*8              :: PEDGE  ! Pressure @ bottom edge of (I,J,L) [hPa]
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
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( GEOS_4 ) || defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

      !-----------------------------
      ! GEOS-4 & GEOS-5 met fields
      !-----------------------------

#if defined( ESMF_ )
      ! Pressure [hPa] at bottom edge of level L (see documentation header)
      ! Taken from the GCM fields
      PEDGE = EXTERNAL_PEDGE(I,J,L)
#else
      ! Pressure [hPa] at bottom edge of level L (see documentation header)
      ! Computed for use w/in GEOS-Chem
      PEDGE = AP(L) + ( BP(L) * PFLT(I,J) )
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
      ! Pressure [hPa] at bottom edge of level L (see documentation header)
      ! Computed for use within GEOS-Chem
      PEDGE = AP(L) + ( BP(L) * ( PFLT(I,J) - PTOP ) )
#endif

#endif     
 
      END FUNCTION GET_PEDGE 
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
      REAL*8              :: PEDGE  ! Pressure @ bottom edge of (I,J,L) [hPa]
!
! !REVISION HISTORY:
!  (1 ) Modified from GET_PEDGE (cdh, 1/22/09)
!  02 Feb 2012 - R. Yantosca - Compute PEDGE for GEOS-5.7.2 the same as MERRA
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_PEDGE_FULLGRID begins here!
      !=================================================================

#if   defined( GEOS_4 ) || defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

      !-----------------------------
      ! GEOS-4 & GEOS-5 met fields
      !-----------------------------

      ! Pressure [hPa] at bottom edge of level L (see documentation header)
      PEDGE = AP_FULLGRID(L) + ( BP_FULLGRID(L) * PFLT(I,J) )

#else 

      !-----------------------------
      ! GEOS-3 & GCAP met fields
      !-----------------------------

      ! Pressure [hPa] at bottom edge of level L (see documentation header)
      PEDGE = AP_FULLGRID(L) + ( BP_FULLGRID(L) * ( PFLT(I,J) - PTOP ) )

#endif     
 
      ! Return to calling program
      END FUNCTION GET_PEDGE_FULLGRID
!EOC
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
      REAL*8              :: PCENTER  ! Pressure @ center of (I,J,L) [hPa]
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
      PCENTER = 0.5d0 * ( GET_PEDGE(I,J,L) + GET_PEDGE(I,J,L+1) )

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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS
      INTEGER :: L

      ! Temporary array for GCAP met fields (swu, bmy, 5/24/05)
      REAL*8  :: SIGE_GCAP(LLPAR+1)  

      !=================================================================
      ! INIT_PRESSURE begins here!
      !=================================================================
      ALLOCATE( PFLT( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PFLT' )
      PFLT = 0d0

      ALLOCATE( AP( LLPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AP' )
      AP = 1d0

      ALLOCATE( BP( LLPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BP' )
      BP = 0d0

      ALLOCATE( AP_FULLGRID( LGLOB+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AP_FULLGRID' )
      AP = 1d0

      ALLOCATE( BP_FULLGRID( LGLOB+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BP_FULLGRID' )
      BP = 0d0

#if defined( ESMF_ )
      ALLOCATE( EXTERNAL_PEDGE( IIPAR, JJPAR, LLPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EXTERNAL_PEDGE' )
      EXTERNAL_PEDGE = 0d0
#endif

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

      !=================================================================
      ! GEOS-5/MERRA vertical coordinates (47 or 72 levels)
      !=================================================================

#if   defined( GRIDREDUCED )

      !-----------------------------------------------------------------
      ! GEOS-5/MERRA 47-level reduced vertical grid
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
      AP = (/  0.000000d0,   0.000000d0,  12.704939d0,  35.465965d0, 
     &        66.098427d0, 101.671654d0, 138.744400d0, 173.403183d0, 
     &       198.737839d0, 215.417526d0, 223.884689d0, 224.362869d0, 
     &       216.864929d0, 201.192093d0, 176.929993d0, 150.393005d0, 
     &       127.837006d0, 108.663429d0,  92.365662d0,  78.512299d0, 
     &        56.387939d0,  40.175419d0,  28.367815d0,  19.791553d0, 
     &         9.292943d0,   4.076567d0,   1.650792d0,   0.616779d0, 
     &         0.211349d0,   0.066000d0,   0.010000d0 /)

      ! Bp [unitless] for 30 levels (31 edges)
      BP = (/  1.000000d0,   0.985110d0,   0.943290d0,   0.867830d0, 
     &         0.764920d0,   0.642710d0,   0.510460d0,   0.378440d0, 
     &         0.270330d0,   0.183300d0,   0.115030d0,   0.063720d0, 
     &         0.028010d0,   0.006960d0,   0.000000d0,   0.000000d0, 
     &         0.000000d0,   0.000000d0,   0.000000d0,   0.000000d0, 
     &         0.000000d0,   0.000000d0,   0.000000d0,   0.000000d0, 
     &         0.000000d0,   0.000000d0,   0.000000d0,   0.000000d0, 
     &         0.000000d0,   0.000000d0,   0.000000d0 /)

      !-----------------------------------
      ! GEOS-4 55 level grid
      !-----------------------------------

      ! AP [hPa] for 55 levels (56 edges)
      AP_FULLGRID =
     &     (/ 0.000000d0,   0.000000d0,  12.704939d0,  35.465965d0, 
     &       66.098427d0, 101.671654d0, 138.744400d0, 173.403183d0,
     &      198.737839d0, 215.417526d0, 223.884689d0, 224.362869d0,
     &      216.864929d0, 201.192093d0, 176.929993d0, 150.393005d0,
     &      127.837006d0, 108.663429d0,  92.365662d0,  78.512299d0, 
     &       66.603378d0,  56.387939d0,  47.643932d0,  40.175419d0, 
     &       33.809956d0,  28.367815d0,  23.730362d0,  19.791553d0, 
     &       16.457071d0,  13.643393d0,  11.276889d0,   9.292943d0,
     &        7.619839d0,   6.216800d0,   5.046805d0,   4.076567d0, 
     &        3.276433d0,   2.620212d0,   2.084972d0,   1.650792d0,
     &        1.300508d0,   1.019442d0,   0.795134d0,   0.616779d0, 
     &        0.475806d0,   0.365041d0,   0.278526d0,   0.211349d0, 
     &        0.159495d0,   0.119703d0,   0.089345d0,   0.066000d0, 
     &        0.047585d0,   0.032700d0,   0.020000d0,   0.010000d0 /)

      ! BP [unitless] for 55 levels (56 edges)
      BP_FULLGRID =
     &     (/  1.000000d0,  0.985110d0,   0.943290d0,   0.867830d0,
     &         0.764920d0,  0.642710d0,   0.510460d0,   0.378440d0,
     &         0.270330d0,  0.183300d0,   0.115030d0,   0.063720d0,
     &         0.028010d0,  0.006960d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0 /)

#else

      !-----------------------------------------------------------------
      ! GEOS-4 55 level grid
      !-----------------------------------------------------------------

      ! AP [hPa] for 55 levels (56 edges)
      AP = (/ 0.000000d0,   0.000000d0,  12.704939d0,  35.465965d0, 
     &       66.098427d0, 101.671654d0, 138.744400d0, 173.403183d0,
     &      198.737839d0, 215.417526d0, 223.884689d0, 224.362869d0,
     &      216.864929d0, 201.192093d0, 176.929993d0, 150.393005d0,
     &      127.837006d0, 108.663429d0,  92.365662d0,  78.512299d0, 
     &       66.603378d0,  56.387939d0,  47.643932d0,  40.175419d0, 
     &       33.809956d0,  28.367815d0,  23.730362d0,  19.791553d0, 
     &       16.457071d0,  13.643393d0,  11.276889d0,   9.292943d0,
     &        7.619839d0,   6.216800d0,   5.046805d0,   4.076567d0, 
     &        3.276433d0,   2.620212d0,   2.084972d0,   1.650792d0,
     &        1.300508d0,   1.019442d0,   0.795134d0,   0.616779d0, 
     &        0.475806d0,   0.365041d0,   0.278526d0,   0.211349d0, 
     &        0.159495d0,   0.119703d0,   0.089345d0,   0.066000d0, 
     &        0.047585d0,   0.032700d0,   0.020000d0,   0.010000d0 /)

      ! BP [unitless] for 55 levels (56 edges)
      BP = (/  1.000000d0,  0.985110d0,   0.943290d0,   0.867830d0,
     &         0.764920d0,  0.642710d0,   0.510460d0,   0.378440d0,
     &         0.270330d0,  0.183300d0,   0.115030d0,   0.063720d0,
     &         0.028010d0,  0.006960d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0,
     &         0.000000d0,  0.000000d0,   0.000000d0,   0.000000d0 /)


#endif

#elif defined( GCAP )

      !=================================================================
      ! GCAP vertical coordinates (23 levels)
      !=================================================================

      ! Define SIGMA edges from GISS 23L model
      SIGE_GCAP = (/ 1.0d0,           0.9712230d0,     0.9340528d0,     
     &               0.8800959d0,     0.8021583d0,     0.6714628d0, 
     &               0.5035971403d0,  0.3297362030d0,  0.1966426820d0,  
     &               0.1139088720d0,  0.0503597111d0,  0.0000000000d0, 
     &              -0.0395683460d0, -0.0764988065d0, -0.1124700233d0, 
     &              -0.1419664323d0, -0.1585131884d0, -0.1678657085d0, 
     &              -0.1743045598d0, -0.1781055182d0, -0.1793033630d0, 
     &              -0.1796822548d0, -0.1798187047d0, -0.1798536479d0 /)

      ! Convert SIGMA to AP and BP coordinates
      DO  L = 1, LLTROP
         AP(L) = 150d0 + ( SIGE_GCAP(L) * ( PTOP - 150d0 ) )
         BP(L) = SIGE_GCAP(L)
      ENDDO

      DO L = LLTROP+1, LLPAR+1  
         AP(L) = 150d0 + ( SIGE_GCAP(L) * ( 984d0 - 150d0 ) )
         BP(L) = 0d0
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
      IF ( ALLOCATED( PFLT        ) ) DEALLOCATE( PFLT        )
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

