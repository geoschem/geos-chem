!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: transfer_mod.F90
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
#ifdef EXCHANGE
!
! !USES:
!
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
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE WILL BE A STUB UNLESS GEOS-Chem IS COMPILED    %%%
!  %%%  WITH THE EXCHANGE=y OPTION. (bmy, 10/4/19)                       %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
!  See https://github.com/geoschem/geos-chem for complete history
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
!  See https://github.com/geoschem/geos-chem for complete history
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    EDGE_IN = (/ 0.000000e+00_fp, 4.804826e-02_fp, &
                 6.593752e+00_fp, 1.313480e+01_fp, &
                 1.961311e+01_fp, 2.609201e+01_fp, &
                 3.257081e+01_fp, 3.898201e+01_fp, &
                 4.533901e+01_fp, 5.169611e+01_fp, &
                 5.805321e+01_fp, 6.436264e+01_fp, &
                 7.062198e+01_fp, 7.883422e+01_fp, &
                 8.909992e+01_fp, 9.936521e+01_fp, &
                 1.091817e+02_fp, 1.189586e+02_fp, &
                 1.286959e+02_fp, 1.429100e+02_fp, &
                 1.562600e+02_fp, 1.696090e+02_fp, &
                 1.816190e+02_fp, 1.930970e+02_fp, &
                 2.032590e+02_fp, 2.121500e+02_fp, &
                 2.187760e+02_fp, 2.238980e+02_fp, &
                 2.243630e+02_fp, 2.168650e+02_fp, &
                 2.011920e+02_fp,
    !---- EDGES OF GEOS-5 FIXED PRESSURE LEVELS OCCUR BELOW THIS LINE ------
                 1.769300e+02_fp, &
                 1.503930e+02_fp, 1.278370e+02_fp, &
                 1.086630e+02_fp, 9.236572e+01_fp, &
                 7.851231e+01_fp, 6.660341e+01_fp, &
                 5.638791e+01_fp, 4.764391e+01_fp, &
                 4.017541e+01_fp, 3.381001e+01_fp, &
                 2.836781e+01_fp, 2.373041e+01_fp, &
                 1.979160e+01_fp, 1.645710e+01_fp, &
                 1.364340e+01_fp, 1.127690e+01_fp, &
                 9.292942e+00_fp, 7.619842e+00_fp, &
                 6.216801e+00_fp, 5.046801e+00_fp, &
                 4.076571e+00_fp, 3.276431e+00_fp, &
                 2.620211e+00_fp, 2.084970e+00_fp, &
                 1.650790e+00_fp, 1.300510e+00_fp, &
                 1.019440e+00_fp, 7.951341e-01_fp, &
                 6.167791e-01_fp, 4.758061e-01_fp, &
                 3.650411e-01_fp, 2.785261e-01_fp, &
                 2.113490e-01_fp, 1.594950e-01_fp, &
                 1.197030e-01_fp, 8.934502e-02_fp, &
                 6.600001e-02_fp, 4.758501e-02_fp, &
                 3.270000e-02_fp, 2.000000e-02_fp, &
                 1.000000e-02_fp /)

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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( ALLOCATED( EDGE_IN ) ) DEALLOCATE( EDGE_IN )

  END SUBROUTINE CLEANUP_TRANSFER
#endif
!EOC
END MODULE TRANSFER_MOD
