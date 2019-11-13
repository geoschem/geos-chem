!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: transfer_mod
!
! !DESCRIPTION: Module TRANSFER\_MOD contains routines used to copy data
!  from REAL*4 to REAL(fp) arrays after being read from disk.  Also, vertical
!  levels will be collapsed in the stratosphere if necessary.  This will help
!  us to gain computational advantage.
!\\
!\\
! !INTERFACE:
!
      MODULE TRANSFER_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,      ONLY : ALLOC_ERR
      USE ERROR_MOD,      ONLY : GEOS_CHEM_STOP
      USE PRECISION_MOD
      USE State_Grid_Mod, ONLY : GrdState

      IMPLICIT NONE

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: TRANSFER_3D_yan
      PUBLIC  :: INIT_TRANSFER
      PUBLIC  :: CLEANUP_TRANSFER
!
! !REMARKS:
!
!  Hybrid Grid Coordinate Definition: (dsa, bmy, 8/27/02, 8/11/15)
!  ============================================================================
!                                                                             .
!  GEOS-4, GEOS-5, GEOS-FP, MERRA, and MERRA-2 (hybrid grids):
!  ----------------------------------------------------------------------------
!  For GEOS-4 and GEOS-5, the pressure at the bottom edge of grid box (I,J,L)
!  is defined as follows:
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
!  formula as for GEOS-4, with one modification:
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
!  The following are true for GCAP, GEOS-3, GEOS-4:
!  ----------------------------------------------------------------------------
!  (1) Bp(NZ+1) = 0.0       (L=NZ+1 is the atmosphere top)
!  (2) Bp(1)    = 1.0       (L=1    is the surface       )
!  (3) PTOP     = Ap(NZ+1)  (L=NZ+1 is the atmosphere top)
!
! !REVISION HISTORY:
!  21 Sep 2010 - M. Evans    - Initial version
!  (1 ) GEOS-3 Output levels were determined by Mat Evans.  Groups of 2 levels
!        and groups of 4 levels on the original grid are merged together into
!        thick levels for the output grid. (mje, bmy, 9/26/01)
!  (2 ) Assumes that LLPAR == LGLOB for GEOS-1, GEOS-STRAT (bmy, 9/26/01)
!  (3 ) EDGE_IN needs to be provided for each model type, within an #ifdef
!        block, in order to ensure compilation.  However, EDGE_IN is currently
!        only used for regridding GEOS-3 data (and probably also GEOS-4 when
!        that becomes available). (bmy, 9/26/01)
!  (4 ) Add interfaces TRANSFER_2D and TRANSFER_ZONAL (bmy, 9/27/01)
!  (5 ) Added routine TRANSFER_2D_R4.  Added TRANSFER_2D_R4 to the generic
!        TRANSFER_2D interface. (bmy, 1/25/02)
!  (6 ) Updated comments, cosmetic changes (bmy, 2/28/02)
!  (7 ) Bug fix: remove extraneous "," in GEOS-1 definition of EDGE_IN array.
!        (bmy, 3/25/02)
!  (8 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Also add MODULE INTERFACES section,
!        since we have an interface here. (bmy, 5/28/02)
!  (9 ) Now references "pressure_mod.f" (dsa, bdf, bmy, 8/22/02)
!  (10) Bug fix in "init_transfer", declare variable L.  Also reference
!        GEOS_CHEM_STOP from "error_mod.f" for safe stop (bmy, 10/15/02)
!  (11) Added routine TRANSFER_3D_TROP.  Also updated comments. (bmy, 10/31/02)
!  (12) Now uses functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        (bmy, 3/11/03)
!  (13) Added code to regrid GEOS-4 from 55 --> 30 levels.  Renamed module
!        variable SIGE_IN to EDGE_IN. (mje, bmy, 10/31/03)
!  (14) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/24/05)
!  (15) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (16) Modified for GEOS-5.  Rewritten for clarity. (bmy, 10/30/07)
!  13 Aug 2010 - R. Yantosca - Added modifications for MERRA met fields
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  02 Feb 2012 - R. Yantosca - Added modifications for GEOS-5.7.x met fields
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Updated to use grid_mod.F90 for the GI model
!  20 Jul 2012 - R. Yantosca - Add routine TRANSFER_3D_Bry, which takes
!                              data sized (144,91,:) as inputs & outputs
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  29 Oct 2013 - R. Yantosca - Remove TRANSFER_3D_NOLUMP routine, we can just
!                              instead do a direct cast assignment
!  03 Apr 2014 - R. Yantosca - Add TRANSFER_3D_R4 and TRANSFER_3D_R8 routines
!                              so that they can be overloaded w/ an interface
!  06 Nov 2014 - R. Yantosca - Remove obsolete TRANSFER_A6 function
!  06 Nov 2014 - R. Yantosca - Remove obsolete TRANSFER_ZONAL* functions
!  06 Nov 2014 - R. Yantosca - Remove obsolete TRANSFER_TO_1D function
!  06 Nov 2014 - R. Yantosca - Remove obsolete TRANSFER_2D* functions
!  06 Nov 2014 - R. Yantosca - Remove obsolete TRANSFER_3D_TROP function
!  04 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  11 Aug 2015 - R. Yantosca - Add support for MERRA2 data
!  14 Aug 2018 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and
!                              MERRA;
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Scalars
      INTEGER             :: I0
      INTEGER             :: J0
      INTEGER             :: L_COPY

      ! Arrays
      REAL(fp), ALLOCATABLE :: EDGE_IN(:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Transfer_3d_yan
!
! !DESCRIPTION: Subroutine TRANSFER\_3D\_YAN
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_3D_yan( NI, NJ, NK, IN, OUT )
!
! !INPUT PARAMETERS:
!
      INTEGER,  INTENT(IN)  :: NI, NJ, NK
      REAL*4,   INTENT(IN)  :: IN(NI,NJ,NK)    ! Input data
!
! !OUTPUT PARAMETERS:
!
      REAL(fp), INTENT(OUT) :: OUT(NI,NJ,NK)   ! Output data
!
! !REVISION HISTORY:
!  08 Feb 2007 - R. Yantosca - Initial version
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA the same way as GEOS-5, because
!                              the vertical grids are identical
!  02 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x the same way as MERRA
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  11 Aug 2015 - R. Yantosca - MERRA2 behaves as GEOS-5, MERRA, GEOS-FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER   :: I,J,K

      !=================================================================
      ! TRANSFER_3D_Lp1 begins here!
      !=================================================================

      ! Copy the first L_COPY+1 levels
      OUT(:,:,:) = IN(:,:,:)

      END SUBROUTINE TRANSFER_3D_yan
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Transfer
!
! !DESCRIPTION: Subroutine INIT\_TRANSFER initializes and zeroes
!  all module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TRANSFER( State_Grid, THIS_I0, THIS_J0 )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      TYPE(GrdState), INTENT(IN) :: State_Grid ! Grid State object
      INTEGER,        INTENT(IN) :: THIS_I0    ! Global X (longitude) offset
      INTEGER,        INTENT(IN) :: THIS_J0    ! Global Y (latitude)  offset
!
! !REVISION HISTORY:
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Removed additional "," for GEOS-1 definition of EDGE_IN (bmy, 3/25/02)
!  (2 ) Now use GET_BP from "pressure_mod.f" to get sigma edges for all
!        grids except GEOS-3 (dsa, bdf, bmy, 8/22/02)
!  (3 ) Declare L as a local variable.  Also reference ALLOC_ERR from module
!        "error_mod.f" (bmy, 10/15/02)
!  (4 ) Renamed SIGE_IN to EDGE_IN to denote that it is not always a sigma
!        coordinate (as for GEOS-4).  Now assign original Ap coordinates from
!        the GEOS-4 grid to the EDGE_IN array (bmy, 10/31/03)
!  (5 ) Now modified for GEOS-5 met fields (bmy, 5/24/05)
!  (6 ) Rewritten for clarity.  Remove references to "grid_mod.f" and
!        "pressure_mod.f".  Now pass I0, J0 from "grid_mod.f" via the arg list.
!         (bmy, 2/8/07)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA the same way as GEOS-5, because
!                              the vertical grids are identical
!  02 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x the same way as MERRA
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  12 Aug 2015 - R. Yantosca - Treat MERRA2 in the same way as GEOS-FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS, L

      !=================================================================
      ! INIT_TRANSFER begins here!
      !=================================================================

      ! Return if we have already initialized
      IF ( IS_INIT ) RETURN

      !-----------------------------------------------------------------
      ! Get global X and Y offsets (usually =0, even for nested grid)
      !-----------------------------------------------------------------
      I0 = THIS_I0
      J0 = THIS_J0

      !-----------------------------------------------------------------
      ! Get the # of levels to copy in the vertical
      !-----------------------------------------------------------------
      IF ( State_Grid%NZ == State_Grid%NativeNZ ) THEN

         ! Full vertical resolution; copy all levels!
         L_COPY = State_Grid%NativeNZ

      ELSE

         ! Copy up to L=36 (GEOS-FP, MERRA-2)
         L_COPY = 36

      ENDIF

      !=================================================================
      ! Define vertical edges for collapsing stratospheric levels
      !=================================================================

      ! Allocate the EDGE_IN array
      ALLOCATE( EDGE_IN( State_Grid%NativeNZ + 1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EDGE_IN' )
      EDGE_IN = 0e+0_fp

      !-----------------------------------------------------------------
      ! Levels 1-31 are "terrain-following" coordinates
      ! (i.e. vary with location), and levels 32-72 are
      ! fixed pressure levels.  The transition pressure is 176.93 hPa,
      ! which is the edge between L=31 and L=32.
      !
      ! Initialize EDGE_IN with the original 73 Ap values for GEOS-5.
      !-----------------------------------------------------------------
      EDGE_IN = (/
     &  0.000000e+00_fp, 4.804826e-02_fp,
     &  6.593752e+00_fp, 1.313480e+01_fp,
     &  1.961311e+01_fp, 2.609201e+01_fp,
     &  3.257081e+01_fp, 3.898201e+01_fp,
     &  4.533901e+01_fp, 5.169611e+01_fp,
     &  5.805321e+01_fp, 6.436264e+01_fp,
     &  7.062198e+01_fp, 7.883422e+01_fp,
     &  8.909992e+01_fp, 9.936521e+01_fp,
     &  1.091817e+02_fp, 1.189586e+02_fp,
     &  1.286959e+02_fp, 1.429100e+02_fp,
     &  1.562600e+02_fp, 1.696090e+02_fp,
     &  1.816190e+02_fp, 1.930970e+02_fp,
     &  2.032590e+02_fp, 2.121500e+02_fp,
     &  2.187760e+02_fp, 2.238980e+02_fp,
     &  2.243630e+02_fp, 2.168650e+02_fp,
     &  2.011920e+02_fp,
!------- EDGES OF GEOS-5 FIXED PRESSURE LEVELS OCCUR BELOW THIS LINE ------
     &  1.769300e+02_fp,
     &  1.503930e+02_fp, 1.278370e+02_fp,
     &  1.086630e+02_fp, 9.236572e+01_fp,
     &  7.851231e+01_fp, 6.660341e+01_fp,
     &  5.638791e+01_fp, 4.764391e+01_fp,
     &  4.017541e+01_fp, 3.381001e+01_fp,
     &  2.836781e+01_fp, 2.373041e+01_fp,
     &  1.979160e+01_fp, 1.645710e+01_fp,
     &  1.364340e+01_fp, 1.127690e+01_fp,
     &  9.292942e+00_fp, 7.619842e+00_fp,
     &  6.216801e+00_fp, 5.046801e+00_fp,
     &  4.076571e+00_fp, 3.276431e+00_fp,
     &  2.620211e+00_fp, 2.084970e+00_fp,
     &  1.650790e+00_fp, 1.300510e+00_fp,
     &  1.019440e+00_fp, 7.951341e-01_fp,
     &  6.167791e-01_fp, 4.758061e-01_fp,
     &  3.650411e-01_fp, 2.785261e-01_fp,
     &  2.113490e-01_fp, 1.594950e-01_fp,
     &  1.197030e-01_fp, 8.934502e-02_fp,
     &  6.600001e-02_fp, 4.758501e-02_fp,
     &  3.270000e-02_fp, 2.000000e-02_fp,
     &  1.000000e-02_fp /)

      ! We have now initialized everything
      IS_INIT = .TRUE.

      END SUBROUTINE INIT_TRANSFER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Transfer
!
! !DESCRIPTION: Subroutine CLEANUP\_TRANSFER deallocates all module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_TRANSFER
!
! !REVISION HISTORY:
!  19 Sep 2001 - R. Yantosca - Initial version
!  31 Oct 2003 - R. Yantosca - Renamed SIGE_IN to EDGE_IN to denote that it
!                              is not always a sigma coordinate (as for GEOS-4)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_TRANSFER begins here!
      !=================================================================
      IF ( ALLOCATED( EDGE_IN ) ) DEALLOCATE( EDGE_IN )

      END SUBROUTINE CLEANUP_TRANSFER
!EOC
      END MODULE TRANSFER_MOD
