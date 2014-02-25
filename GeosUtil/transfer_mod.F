!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: transfer_mod
!
! !DESCRIPTION: Module TRANSFER\_MOD contains routines used to copy data 
!  from REAL*4 to REAL*8 arrays after being read from disk.  Also, vertical 
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
      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      IMPLICIT NONE

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: TRANSFER_A6
      PUBLIC  :: TRANSFER_2D
      PUBLIC  :: TRANSFER_3D
      PUBLIC  :: TRANSFER_3D_Bry
      PUBLIC  :: TRANSFER_3D_Lp1
      PUBLIC  :: TRANSFER_3D_TROP
      PUBLIC  :: TRANSFER_G5_PLE
      PUBLIC  :: TRANSFER_ZONAL
      PUBLIC  :: TRANSFER_TO_1D
      PUBLIC  :: INIT_TRANSFER
      PUBLIC  :: CLEANUP_TRANSFER

      INTERFACE TRANSFER_2D
         MODULE PROCEDURE TRANSFER_2D_INT
         MODULE PROCEDURE TRANSFER_2D_R4
         MODULE PROCEDURE TRANSFER_2D_R8
      END INTERFACE

      INTERFACE TRANSFER_ZONAL
         MODULE PROCEDURE TRANSFER_ZONAL_R4
         MODULE PROCEDURE TRANSFER_ZONAL_R8
      END INTERFACE
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: LUMP_2
      PRIVATE :: LUMP_2_R4
      PRIVATE :: LUMP_2_R8
      PRIVATE :: LUMP_4
      PRIVATE :: LUMP_4_R4
      PRIVATE :: LUMP_4_R8
      PRIVATE :: TRANSFER_2D_INT
      PRIVATE :: TRANSFER_2D_R4
      PRIVATE :: TRANSFER_2D_R8
      PRIVATE :: TRANSFER_ZONAL_R4
      PRIVATE :: TRANSFER_ZONAL_R8

      INTERFACE LUMP_2
         MODULE PROCEDURE LUMP_2_R4
         MODULE PROCEDURE LUMP_2_R8
      END INTERFACE

      INTERFACE LUMP_4
         MODULE PROCEDURE LUMP_4_R4
         MODULE PROCEDURE LUMP_4_R8
      END INTERFACE
!
! !REMARKS:
!
!  Hybrid Grid Coordinate Definition: (dsa, bmy, 8/27/02, 2/2/12)
!  ============================================================================
!                                                                             .
!  GEOS-4, GEOS-5, GEOS-5.7, and MERRA (hybrid grids):
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
!  (1) Bp(LLPAR+1) = 0.0          (L=LLPAR+1 is the atmosphere top)
!  (2) Bp(1)       = 1.0          (L=1       is the surface       )
!  (3) PTOP        = Ap(LLPAR+1)  (L=LLPAR+1 is the atmosphere top) 
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
      REAL*8, ALLOCATABLE :: EDGE_IN(:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: transfer_a6
!
! !DESCRIPTION: Subroutine TRANSFER\_A6 transfers A-6 data from a REAL*4 
!  array to a REAL*8 array.  Vertical layers are collapsed (from LGLOB to 
!  LLPAR) if necessary.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_A6( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR,LGLOB)     ! Input data
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: OUT(LLPAR,IIPAR,JJPAR)    ! Output data
!
! !REMARKS:
! 
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) A-6 fields are dimensioned (LLPAR,IIPAR,JJPAR) since for Fortran
!        efficiency, since the code loops over vertical layers L in a column
!        located above a certain surface location (I,J). (bmy, 9/21/01)
!  (2 ) Assumes that LLPAR == LGLOB for GEOS-1, GEOS-STRAT (bmy, 9/21/01)
!  (3 ) Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        Now I0, J0 are local variables. (bmy, 3/11/03)
!  (4 ) Added code to regrid GEOS-4 from 55 --> 30 levels (mje, bmy, 10/31/03)
!  (5 ) Now modified for GEOS-5 met fields (bmy, 5/24/05)
!  (6 ) Rewritten for clarity (bmy, 2/8/07)
!  (7 ) Now get nested-grid offsets (dan, bmy, 11/6/08)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA the same way as GEOS-5, because
!                              the vertical grids are identical
!  02 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x the same way as MERRA
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L
      REAL*4  :: INCOL(LGLOB)
      
      !================================================================
      ! TRANSFER_A6 begins here!
      !================================================================

      ! Copy the first L_COPY levels
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      DO L = 1, L_COPY
         OUT(L,I,J) = IN(I+I0,J+J0,L)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Exit if we are running at full vertical resolution
      IF ( LLPAR == LGLOB ) RETURN

      !=================================================================
      ! Collapse levels in the stratosphere 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, INCOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Store vertical column at (IREF,JREF) in a 1-D vector
         INCOL = IN( I+I0, J+J0, 1:LGLOB )
            
#if   defined( GEOS_4 )

         !--------------------------------------------------------------
         ! GEOS-4: Lump 55 levels into 30 levels, starting above L=20
         ! Lump levels in groups of 2, then 4. (cf. Mat Evans)
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time
         OUT(20,I,J) = LUMP_2( INCOL, LGLOB, 20 )
         OUT(21,I,J) = LUMP_2( INCOL, LGLOB, 22 )
         OUT(22,I,J) = LUMP_2( INCOL, LGLOB, 24 )
         OUT(23,I,J) = LUMP_2( INCOL, LGLOB, 26 )

         ! Lump 4 levels together at a time
         OUT(24,I,J) = LUMP_4( INCOL, LGLOB, 28 )
         OUT(25,I,J) = LUMP_4( INCOL, LGLOB, 32 )
         OUT(26,I,J) = LUMP_4( INCOL, LGLOB, 36 ) 
         OUT(27,I,J) = LUMP_4( INCOL, LGLOB, 40 ) 
         OUT(28,I,J) = LUMP_4( INCOL, LGLOB, 44 ) 
         OUT(29,I,J) = LUMP_4( INCOL, LGLOB, 48 ) 
         OUT(30,I,J) = LUMP_4( INCOL, LGLOB, 52 ) 

#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

         !--------------------------------------------------------------
         ! GEOS-5/MERRA: Lump 72 levels into 47 levels, starting above 
         ! L=36.  Lump levels in groups of 2, then 4. (cf. Bob Yantosca)
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time
         OUT(37,I,J) = LUMP_2( INCOL, LGLOB, 37 )
         OUT(38,I,J) = LUMP_2( INCOL, LGLOB, 39 )
         OUT(39,I,J) = LUMP_2( INCOL, LGLOB, 41 )
         OUT(40,I,J) = LUMP_2( INCOL, LGLOB, 43 )

         ! Lump 4 levels together at a time
         OUT(41,I,J) = LUMP_4( INCOL, LGLOB, 45 )
         OUT(42,I,J) = LUMP_4( INCOL, LGLOB, 49 )
         OUT(43,I,J) = LUMP_4( INCOL, LGLOB, 53 ) 
         OUT(44,I,J) = LUMP_4( INCOL, LGLOB, 57 ) 
         OUT(45,I,J) = LUMP_4( INCOL, LGLOB, 61 ) 
         OUT(46,I,J) = LUMP_4( INCOL, LGLOB, 65 ) 
         OUT(47,I,J) = LUMP_4( INCOL, LGLOB, 69 ) 

#endif

         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE TRANSFER_A6
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: transfer_3d
!
! !DESCRIPTION: Subroutine TRANSFER\_3D transfers 3-dimensional data from a 
!  REAL*4  array to a REAL*8 array.  Vertical layers are collapsed (from LGLOB
!  to LLPAR) if necessary.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_3D( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR,LGLOB)    ! Input data
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLPAR)   ! Output data
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Lump levels together in groups of 2 or 4, as dictated by Mat Evans.
!        (bmy, 9/21/01)
!  (2 ) Assumes that LLPAR == LGLOB for GEOS-1, GEOS-STRAT (bmy, 9/21/01)
!  (3 ) Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        Now I0, J0 are local variables. (bmy, 3/11/03)
!  (4 ) Added code to regrid GEOS-4 from 55 --> 30 levels (mje, bmy, 10/31/03)
!  (5 ) Now modified for GEOS-5 met fields (bmy, 5/24/05)
!  (6 ) Rewritten for clarity (bmy, 2/8/07)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA the same way as GEOS-5, because
!                              the vertical grids are identical
!  02 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x the same way as MERRA
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J
      REAL*4  :: INCOL(LGLOB)
     
      !================================================================
      ! TRANSFER_3D begins here!
      !================================================================

      ! Copy the first L_COPY levels
      OUT(:,:,1:L_COPY) = IN( 1+I0:IIPAR+I0, 1+J0:JJPAR+J0, 1:L_COPY )

      ! Exit if we are running at full vertical resolution
      IF ( LLPAR == LGLOB ) RETURN

      !================================================================
      ! Collapse levels in the stratosphere
      !================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, INCOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Copy a vertical column into INCOL
         INCOL = IN( I+I0, J+J0, 1:LGLOB )

#if   defined( GEOS_4 )

         !--------------------------------------------------------------
         ! GEOS-4: Lump 55 levels into 30 levels, starting above L=20
         ! Lump levels in groups of 2, then 4. (cf. Mat Evans)
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time, starting at L=20
         OUT(I,J,20) = LUMP_2( INCOL, LGLOB, 20 )
         OUT(I,J,21) = LUMP_2( INCOL, LGLOB, 22 )
         OUT(I,J,22) = LUMP_2( INCOL, LGLOB, 24 )
         OUT(I,J,23) = LUMP_2( INCOL, LGLOB, 26 )

         ! Lump 4 levels together at a time, starting at L=28
         OUT(I,J,24) = LUMP_4( INCOL, LGLOB, 28 )
         OUT(I,J,25) = LUMP_4( INCOL, LGLOB, 32 )
         OUT(I,J,26) = LUMP_4( INCOL, LGLOB, 36 ) 
         OUT(I,J,27) = LUMP_4( INCOL, LGLOB, 40 ) 
         OUT(I,J,28) = LUMP_4( INCOL, LGLOB, 44 ) 
         OUT(I,J,29) = LUMP_4( INCOL, LGLOB, 48 ) 
         OUT(I,J,30) = LUMP_4( INCOL, LGLOB, 52 ) 

#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

         !--------------------------------------------------------------
         ! GEOS-5/MERRA Lump 72 levels into 47 levels, starting above 
         ! L=36.  Lump levels in groups of 2, then 4. (cf. Bob Yantosca)
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time
         OUT(I,J,37) = LUMP_2( INCOL, LGLOB, 37 )
         OUT(I,J,38) = LUMP_2( INCOL, LGLOB, 39 )
         OUT(I,J,39) = LUMP_2( INCOL, LGLOB, 41 )
         OUT(I,J,40) = LUMP_2( INCOL, LGLOB, 43 )

         ! Lump 4 levels together at a time
         OUT(I,J,41) = LUMP_4( INCOL, LGLOB, 45 )
         OUT(I,J,42) = LUMP_4( INCOL, LGLOB, 49 )
         OUT(I,J,43) = LUMP_4( INCOL, LGLOB, 53 ) 
         OUT(I,J,44) = LUMP_4( INCOL, LGLOB, 57 ) 
         OUT(I,J,45) = LUMP_4( INCOL, LGLOB, 61 ) 
         OUT(I,J,46) = LUMP_4( INCOL, LGLOB, 65 ) 
         OUT(I,J,47) = LUMP_4( INCOL, LGLOB, 69 ) 

#endif

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE TRANSFER_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: transfer_3d_bry
!
! !DESCRIPTION: Subroutine TRANSFER\_3D\_Bry transfers 3-dimensional data 
!  from a REAL*4 array to a REAL*8 array.  Vertical layers are collapsed 
!  (from LGLOB to LLPAR) if necessary.  This routine was specially created 
!  for 2 x 2.5 Bry data in strat\_chem\_mod.F90.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_3D_Bry( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(144,91,LGLOB)    ! Input data
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: OUT(144,91,LLPAR)   ! Output data
! 
! !REVISION HISTORY: 
!  20 Jul 2012 - R. Yantosca - Initial version, based on TRANSFER_3D
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J
      REAL*4  :: INCOL(LGLOB)
     
      !================================================================
      ! TRANSFER_3D begins here!
      !================================================================

      ! Copy the first L_COPY levels
      OUT(:,:,1:L_COPY) = IN(:,:,1:L_COPY)

      ! Exit if we are running at full vertical resolution
      IF ( LLPAR == LGLOB ) RETURN

      !================================================================
      ! Collapse levels in the stratosphere
      !================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, INCOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, 91
      DO I = 1, 144

         ! Copy a vertical column into INCOL
         INCOL = IN( I, J, 1:LGLOB )

#if   defined( GEOS_4 )

         !--------------------------------------------------------------
         ! GEOS-4: Lump 55 levels into 30 levels, starting above L=20
         ! Lump levels in groups of 2, then 4. (cf. Mat Evans)
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time, starting at L=20
         OUT(I,J,20) = LUMP_2( INCOL, LGLOB, 20 )
         OUT(I,J,21) = LUMP_2( INCOL, LGLOB, 22 )
         OUT(I,J,22) = LUMP_2( INCOL, LGLOB, 24 )
         OUT(I,J,23) = LUMP_2( INCOL, LGLOB, 26 )

         ! Lump 4 levels together at a time, starting at L=28
         OUT(I,J,24) = LUMP_4( INCOL, LGLOB, 28 )
         OUT(I,J,25) = LUMP_4( INCOL, LGLOB, 32 )
         OUT(I,J,26) = LUMP_4( INCOL, LGLOB, 36 ) 
         OUT(I,J,27) = LUMP_4( INCOL, LGLOB, 40 ) 
         OUT(I,J,28) = LUMP_4( INCOL, LGLOB, 44 ) 
         OUT(I,J,29) = LUMP_4( INCOL, LGLOB, 48 ) 
         OUT(I,J,30) = LUMP_4( INCOL, LGLOB, 52 ) 

#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

         !--------------------------------------------------------------
         ! GEOS-5/MERRA Lump 72 levels into 47 levels, starting above 
         ! L=36.  Lump levels in groups of 2, then 4. (cf. Bob Yantosca)
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time
         OUT(I,J,37) = LUMP_2( INCOL, LGLOB, 37 )
         OUT(I,J,38) = LUMP_2( INCOL, LGLOB, 39 )
         OUT(I,J,39) = LUMP_2( INCOL, LGLOB, 41 )
         OUT(I,J,40) = LUMP_2( INCOL, LGLOB, 43 )

         ! Lump 4 levels together at a time
         OUT(I,J,41) = LUMP_4( INCOL, LGLOB, 45 )
         OUT(I,J,42) = LUMP_4( INCOL, LGLOB, 49 )
         OUT(I,J,43) = LUMP_4( INCOL, LGLOB, 53 ) 
         OUT(I,J,44) = LUMP_4( INCOL, LGLOB, 57 ) 
         OUT(I,J,45) = LUMP_4( INCOL, LGLOB, 61 ) 
         OUT(I,J,46) = LUMP_4( INCOL, LGLOB, 65 ) 
         OUT(I,J,47) = LUMP_4( INCOL, LGLOB, 69 ) 

#endif

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE TRANSFER_3D_Bry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: transfer_g5_ple
!
! !DESCRIPTION: Subroutine TRANSFER\_G5\_PLE transfers GEOS-5/MERRA pressure 
!  edge data from the native 72-level grid to the reduced 47-level grid.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_G5_PLE( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR,LGLOB+1)    ! Input data
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLPAR+1)   ! Output data
!
! !REVISION HISTORY: 
!  08 Feb 2007 - R. Yantosca - Initial version
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA the same way as GEOS-5, because
!                              the vertical grids are identical
!  02 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x the same way as MERRA
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J
     
      !================================================================
      ! TRANSFER_PLE begins here!
      !================================================================

      ! Copy the first L_COPY+1 edges (which define L_COPY levels)
      OUT(:,:,1:L_COPY+1) = IN( 1+I0:IIPAR+I0,1+J0:JJPAR+J0,1:L_COPY+1 )

      ! Exit if we are running at full vertical resolution
      IF ( LLPAR == LGLOB ) RETURN

      !================================================================
      ! Return GEOS-5 pressure edges for reduced grid
      !================================================================

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J)
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Top edges of levels lumped by 2's
         OUT(I,J,38) = IN(I+I0,J+J0,39)
         OUT(I,J,39) = IN(I+I0,J+J0,41)
         OUT(I,J,40) = IN(I+I0,J+J0,43)
         OUT(I,J,41) = IN(I+I0,J+J0,45)

         ! Top edges of levels lumped by 4's
         OUT(I,J,42) = IN(I+I0,J+J0,49)
         OUT(I,J,43) = IN(I+I0,J+J0,53)
         OUT(I,J,44) = IN(I+I0,J+J0,57)
         OUT(I,J,45) = IN(I+I0,J+J0,61)
         OUT(I,J,46) = IN(I+I0,J+J0,65)
         OUT(I,J,47) = IN(I+I0,J+J0,69)
         OUT(I,J,48) = IN(I+I0,J+J0,73)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#endif

      END SUBROUTINE TRANSFER_G5_PLE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: transfer_3d_lp1
!
! !DESCRIPTION: Subroutine TRANSFER\_3D\_Lp1 transfers 3-D data from a REAL*4 
!  array of dimension (IIPAR,JJPAR,LGLOB+1) to a REAL*8 array of dimension 
!  (IIPAR,JJPAR,LLPAR+1).  Regrid in the vertical if needed.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_3D_Lp1( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR,LGLOB+1)    ! Input data
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLPAR+1)   ! Output data
! 
! !REVISION HISTORY: 
!  08 Feb 2007 - R. Yantosca - Initial version
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA the same way as GEOS-5, because
!                              the vertical grids are identical
!  02 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x the same way as MERRA
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.  
      INTEGER       :: I, J
      REAL*4        :: INCOL(LGLOB)
     
      !=================================================================
      ! TRANSFER_3D_Lp1 begins here!
      !=================================================================

      ! Copy the first L_COPY+1 levels
      OUT(:,:,1:L_COPY+1) = IN( 1+I0:IIPAR+I0, 1+J0:JJPAR+J0,1:L_COPY+1)

      ! Exit if we are running full vertical resolution
      IF ( LLPAR == LGLOB ) RETURN

      !=================================================================
      ! Collapse levels in the stratosphere
      !
      ! %%% TEMPORARY KLUDGE!!!!
      ! %%% NOTE: For now do the same thing as in TRANSFER_G5_PLE, i.e.
      ! %%% return the values at the edges.  The only other field than 
      ! %%% PLE defined on the edges is CMFMC and that is always zero 
      ! %%% above about 120 hPa. (bmy, 2/8/07)
      !=================================================================
      
#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Top edges of levels lumped by 2's
         OUT(I,J,38) = IN(I+I0,J+J0,39)
         OUT(I,J,39) = IN(I+I0,J+J0,41)
         OUT(I,J,40) = IN(I+I0,J+J0,43)
         OUT(I,J,41) = IN(I+I0,J+J0,45)

         ! Top edges of levels lumped by 4's
         OUT(I,J,42) = IN(I+I0,J+J0,49)
         OUT(I,J,43) = IN(I+I0,J+J0,53)
         OUT(I,J,44) = IN(I+I0,J+J0,57)
         OUT(I,J,45) = IN(I+I0,J+J0,61)
         OUT(I,J,46) = IN(I+I0,J+J0,65)
         OUT(I,J,47) = IN(I+I0,J+J0,69)
         OUT(I,J,48) = IN(I+I0,J+J0,73)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#endif

      END SUBROUTINE TRANSFER_3D_Lp1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: transfer_3d_trop
!
! !DESCRIPTION: Subroutine TRANSFER\_3D\_TROP transfers tropospheric 3-D 
!  data from a REAL*4  array to a REAL*8 array. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_3D_TROP( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR,LLTROP_FIX)     ! Input data
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLTROP_FIX)    ! Output data
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - M. Evans    - Initial version
!  08 Feb 2007 - R. Yantosca - Now use LLTROP_FIX instead of LLTROP, since 
!                              most of the offline simulations use the annual 
!                              mean tropopause
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: L
     
      !=================================================================
      ! TRANSFER_3D_TROP begins here!
      !=================================================================

      ! Cast to REAL*8 and resize up to LLTROP
      DO L = 1, LLTROP_FIX
         CALL TRANSFER_2D( IN(:,:,L), OUT(:,:,L) )
      ENDDO
      
      END SUBROUTINE TRANSFER_3D_TROP 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: transfer_zonal_r4
!
! !DESCRIPTION: Subroutine TRANSFER\_ZONAL\_R4 transfers zonal-mean data from 
!  a REAL*4 array to a REAL*8 array.  Vertical levels are collapsed (from 
!  LGLOB to LLPAR) if necessary. (mje, bmy, 9/21/01, 2/8/07)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_ZONAL_R4( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(JJPAR,LGLOB)     ! Input data
!
! !OUTPUT PARAMETERS:
!
      REAL*4,  INTENT(OUT) :: OUT(JJPAR,LLPAR)    ! Output data
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - M. Evans    - Initial version
!  (1 ) Lump levels together in groups of 2 or 4, as dictated by Mat Evans.
!        (bmy, 9/21/01)
!  (2 ) Assumes that LLPAR == LGLOB for GEOS-1, GEOS-STRAT (bmy, 9/21/01)
!  (3 ) Now use function GET_YOFFSET from "grid_mod.f".  Now I0 and J0 are 
!        local variables (bmy, 3/11/03)
!  (4 ) Added code to regrid GEOS-4 from 55 --> 30 levels (mje, bmy, 10/31/03)
!  (5 ) Rewritten for clarity (bmy, 2/8/07)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA the same way as GEOS-5, because
!                              the vertical grids are identical
!  02 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x the same way as MERRA
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: J
      REAL*4  :: INCOL(LGLOB)
     
      !================================================================
      ! TRANSFER_ZONAL_R4 begins here!
      !================================================================

      ! Copy the first L_COPY levels
      OUT(:,1:L_COPY) = IN( 1+J0:JJPAR+J0, 1:L_COPY )

      ! Exit if we are running at full vertical resolution
      IF ( LLPAR == LGLOB ) RETURN

      !================================================================
      ! Collapse levels in the stratosphere
      !================================================================    

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J, INCOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR

         ! Store vertical column at (I,J) in a 1-D vector
         INCOL = IN( J+J0, 1:LGLOB )

#if   defined( GEOS_4 )
      
         !--------------------------------------------------------------
         ! GEOS-4: Lump 55 levels into 30 levels, starting above L=20
         ! Lump levels in groups of 2, then 4. (cf. Mat Evans)
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time
         OUT(J,20) = LUMP_2( INCOL, LGLOB, 20 )
         OUT(J,21) = LUMP_2( INCOL, LGLOB, 22 )
         OUT(J,22) = LUMP_2( INCOL, LGLOB, 24 )
         OUT(J,23) = LUMP_2( INCOL, LGLOB, 26 )

         ! Lump 4 levels together at a time
         OUT(J,24) = LUMP_4( INCOL, LGLOB, 28 )
         OUT(J,25) = LUMP_4( INCOL, LGLOB, 32 )
         OUT(J,26) = LUMP_4( INCOL, LGLOB, 36 ) 
         OUT(J,27) = LUMP_4( INCOL, LGLOB, 40 ) 
         OUT(J,28) = LUMP_4( INCOL, LGLOB, 44 ) 
         OUT(J,29) = LUMP_4( INCOL, LGLOB, 48 ) 
         OUT(J,30) = LUMP_4( INCOL, LGLOB, 52 ) 

#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

         !--------------------------------------------------------------
         ! GEOS-5/MERRA: Lump 72 levels into 47 levels, starting above 
         ! L=36.  Lump levels in groups of 2, then 4.  
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time
         OUT(J,37) = LUMP_2( INCOL, LGLOB, 37 )
         OUT(J,38) = LUMP_2( INCOL, LGLOB, 39 )
         OUT(J,39) = LUMP_2( INCOL, LGLOB, 41 )
         OUT(J,40) = LUMP_2( INCOL, LGLOB, 43 )

         ! Lump 4 levels together at a time
         OUT(J,41) = LUMP_4( INCOL, LGLOB, 45 )
         OUT(J,42) = LUMP_4( INCOL, LGLOB, 49 )
         OUT(J,43) = LUMP_4( INCOL, LGLOB, 53 ) 
         OUT(J,44) = LUMP_4( INCOL, LGLOB, 57 ) 
         OUT(J,45) = LUMP_4( INCOL, LGLOB, 61 ) 
         OUT(J,46) = LUMP_4( INCOL, LGLOB, 65 ) 
         OUT(J,47) = LUMP_4( INCOL, LGLOB, 69 ) 

#endif

      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE TRANSFER_ZONAL_R4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: transfer_zonal_r8
!
! !DESCRIPTION: Subroutine TRANSFER\_ZONAL\_R8 transfers zonal mean or lat-alt 
!  data from a REAL*4 array of dimension (JJPAR,LGLOB) to a REAL*8 array of 
!  dimension (JJPAR,LLPAR).   Regrid data in the vertical if necessary by
!  lumping levels.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_ZONAL_R8( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(JJPAR,LGLOB)    ! Input data
!
! !OUTPUT PARAMETERS: 
! 
      REAL*8,  INTENT(OUT) :: OUT(JJPAR,LLPAR)   ! Output data
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Lump levels together in groups of 2 or 4, as dictated by Mat Evans.
!        (bmy, 9/21/01)
!  (2 ) Assumes that LLPAR == LGLOB for GEOS-1, GEOS-STRAT (bmy, 9/21/01)
!  (3 ) Now use function GET_YOFFSET from "grid_mod.f".  Now I0 and J0 are 
!        local variables (bmy, 3/11/03)
!  (4 ) Added code to regrid GEOS-4 from 55 --> 30 levels (mje, bmy, 10/31/03)
!  (5 ) Now modified for GEOS-5 met fields (bmy, 5/24/05)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA the same way as GEOS-5, because
!                              the vertical grids are identical
!  02 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x the same way as MERRA
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: J
      REAL*4  :: INCOL(LGLOB)
     
      !================================================================
      ! TRANSFER_ZONAL_R8 begins here!
      !================================================================

      ! Copy the first L_COPY levels
      OUT(:,1:L_COPY) = IN( 1+J0:JJPAR+J0, 1:L_COPY )

      ! Exit if we are running at full vertical resolution
      IF ( LLPAR == LGLOB ) RETURN

      !================================================================
      ! Collapse levels in the stratosphere
      !================================================================    

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J, INCOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR

         ! Store vertical column at (I,J) in a 1-D vector
         INCOL = IN( J+J0, 1:LGLOB )

#if   defined( GEOS_4 )

         !--------------------------------------------------------------
         ! GEOS-4: Lump 55 levels into 30 levels, starting above L=20
         ! Lump levels in groups of 2, then 4. (cf. Mat Evans)
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time
         OUT(J,20) = LUMP_2( INCOL, LGLOB, 20 )
         OUT(J,21) = LUMP_2( INCOL, LGLOB, 22 )
         OUT(J,22) = LUMP_2( INCOL, LGLOB, 24 )
         OUT(J,23) = LUMP_2( INCOL, LGLOB, 26 )

         ! Lump 4 levels together at a time
         OUT(J,24) = LUMP_4( INCOL, LGLOB, 28 )
         OUT(J,25) = LUMP_4( INCOL, LGLOB, 32 )
         OUT(J,26) = LUMP_4( INCOL, LGLOB, 36 ) 
         OUT(J,27) = LUMP_4( INCOL, LGLOB, 40 ) 
         OUT(J,28) = LUMP_4( INCOL, LGLOB, 44 ) 
         OUT(J,29) = LUMP_4( INCOL, LGLOB, 48 ) 
         OUT(J,30) = LUMP_4( INCOL, LGLOB, 52 ) 

#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

         !--------------------------------------------------------------
         ! GEOS-5/MERRA: Lump 72 levels into 47 levels, starting above 
         ! L=36.  Lump levels in groups of 2, then 4.  
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time
         OUT(J,37) = LUMP_2( INCOL, LGLOB, 37 )
         OUT(J,38) = LUMP_2( INCOL, LGLOB, 39 )
         OUT(J,39) = LUMP_2( INCOL, LGLOB, 41 )
         OUT(J,40) = LUMP_2( INCOL, LGLOB, 43 )

         ! Lump 4 levels together at a time
         OUT(J,41) = LUMP_4( INCOL, LGLOB, 45 )
         OUT(J,42) = LUMP_4( INCOL, LGLOB, 49 )
         OUT(J,43) = LUMP_4( INCOL, LGLOB, 53 ) 
         OUT(J,44) = LUMP_4( INCOL, LGLOB, 57 ) 
         OUT(J,45) = LUMP_4( INCOL, LGLOB, 61 ) 
         OUT(J,46) = LUMP_4( INCOL, LGLOB, 65 ) 
         OUT(J,47) = LUMP_4( INCOL, LGLOB, 69 ) 
#endif

      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE TRANSFER_ZONAL_R8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: transfer_2d_int
!
! !DESCRIPTION: Subroutine TRANSFER\_2D\_INT transfers 2-D data from a REAL*4 
!  array of dimension (IIPAR,JJPAR) to an INTEGER array of dimension
!  (IIPAR,JJPAR). 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_2D_INT( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR)    ! Input data
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: OUT(IIPAR,JJPAR)   ! Output data
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Use parallel DO loops to speed things up (bmy, 9/21/01)!
!  (2 ) Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        Now I0 and J0 are local variables. (bmy, 3/11/03)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! TRANSFER_2D_INT begins here!
      !=================================================================

      ! Copy and cast array
      OUT = IN( 1+I0:IIPAR+I0, 1+J0:JJPAR+J0 )

      END SUBROUTINE TRANSFER_2D_INT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: transfer_2d_r4
!
! !DESCRIPTION: Subroutine TRANSFER\_2D\_R4 transfers 2-D data from a REAL*4 
!  array of dimension (IIPAR,JJPAR) to a REAL*4 array of dimension 
!  (IIPAR,JJPAR). 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_2D_R4( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR)
!
! !OUTPUT PARAMETERS:
!
      REAL*4,  INTENT(OUT) :: OUT(IIPAR,JJPAR)
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Use parallel DO loops to speed things up (bmy, 9/21/01)
!  (2 ) Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f"
!        Now I0 and J0 are local variables (bmy, 3/11/03)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! TRANSFER_2D_R4 begins here!
      !=================================================================
      
      ! Copy and cast array
      OUT = IN( 1+I0:IIPAR+I0, 1+J0:JJPAR+J0 )

      END SUBROUTINE TRANSFER_2D_R4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: transfer_2d_r8
!
! !DESCRIPTION: Subroutine TRANSFER\_2D\_R8 transfers 2-D data from a REAL*4 
!  array of dimension (IIPAR,JJPAR) to a REAL*8 array of dimension
!  (IIPAR,JJPAR). 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_2D_R8( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR)    ! Input data
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR)   ! Output data
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Use parallel DO loops to speed things up (bmy, 9/21/01)
!  (2 ) Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f"
!        Now I0 and J0 are local variables. (bmy, 3/11/03)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! TRANSFER_2D_R8 begins here!
      !=================================================================

      ! Copy and cast array
      OUT = IN( 1+I0:IIPAR+I0, 1+J0:JJPAR+J0 )
      
      END SUBROUTINE TRANSFER_2D_R8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: transfer_to_1d
!
! !DESCRIPTION: Subroutine TRANSFER\_TO\_1D transfers 2-D data from a REAL*4 
!  array of dimension (IIPAR,JJPAR) to 1-D a REAL*8 array of dimension (MAXIJ),
!  where MAXIJ = IIPAR * JJPAR.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_TO_1D( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR)    ! Input data
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: OUT(MAXIJ)         ! Output data
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Use single-processor DO-loops for now (bmy, 9/21/01)
!  (2 ) Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        Now I0 and J0 are local variables. (bmy, 3/11/03)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, IREF, J, JREF, IJLOOP

      !=================================================================
      ! TRANSFER_TO_1D begins here!
      !=================================================================

      ! 1-D counter
      IJLOOP = 0

      ! IJLOOP = ( (J-1) * IIPAR ) + I

      DO J = 1, JJPAR
         JREF = J + J0
         DO I = 1, IIPAR
            IREF        = I + I0
            IJLOOP      = IJLOOP + 1
            OUT(IJLOOP) = IN(IREF,JREF)
         ENDDO
      ENDDO
      
      END SUBROUTINE TRANSFER_TO_1D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lump_2_r4
!
! !DESCRIPTION: Function LUMP\_2\_R4 lumps 2 sigma levels into one thick 
!  level.  Input arguments must be REAL*4. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION LUMP_2_R4( IN, L_IN, L ) RESULT( OUT )
!
! !USES:
!
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN) :: IN(L_IN)   ! Column of data on input grid
      INTEGER, INTENT(IN) :: L_IN       ! Vertical dimension of the IN array
      INTEGER, INTENT(IN) :: L          ! Level on input grid from which 
                                        !  to start regridding
!
! !RETURN VALUE:
!
      REAL*4              :: OUT        ! Data on output grid: 4 lumped levels
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Renamed SIGE_IN to EDGE_IN to denote that it is not always a sigma
!        coordinate (as for GEOS-4).  Also updated comments (bmy, 10/31/03)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! LUMP_2_R4 begins here!
      !=================================================================

      ! Error check: prevent array out of bounds error
      IF ( L < 1 .or. L > L_IN .or. L+2 > L_IN+1 ) THEN 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'Error: L < 1 or L > L_IN or L+2 > L_IN+1!'
         WRITE( 6, '(a)' ) 'STOP in LUMP_2 ("regrid_mod.f")'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF
      
      !=================================================================
      ! When lumping the levels together, we need to perform a weighted
      ! average by air mass.  The air mass in a grid box is given by:
      !
      !  Air Mass = Delta-P [hPa] * 100/g * Grid box sfc area [cm2]
      !    
      ! Where Delta-P is the difference in pressure between the bottom 
      ! and top edges of the grid box (Delta-P is positive), 100/g is a 
      ! constant, and the grid box surface area is also constant w/in 
      ! the same vertical column.  Therefore, for a vertical column,
      ! the air mass in a grid box really only depends on Delta-P.
      !
      ! Because of this, we may compute the quantity Q(L') on the new 
      ! merged sigma according to the weighted average (EQUATION 1):
      !
      !             [ ( Q(L  ) * ( PEDGE(L  ) - PEDGE(L+1) ) ) + 
      !               ( Q(L+1) * ( PEDGE(L+1) - PEDGE(L+2) ) ) ]
      !  Q(L') = -------------------------------------------------
      !                       PEDGE(L) - PEDGE(L+2)
      !
      ! where PEDGE(L) is the pressure at the bottom edge of layer L.
      !
      ! GEOS-4/GEOS-5/MERRA are a hybrid sigma-pressure grid, with 
      ! all of the levels above level a certain level being pure 
      ! pressure levels.  Therefore, for these grids, we may just 
      ! use EQUATION 1 exactly as written above.
      !
      ! However, GEOS-3 is a pure sigma grid.  The pressure at the 
      ! edge of a grid box is given by (EQUATION 2):
      ! 
      !  PEDGE(I,J,L) = PTOP + ( SIG_EDGE(L) * ( Psurf(I,J) - PTOP) )
      !
      ! In a vertical column, then ( Psurf(I,J) - PTOP ) will be the 
      ! same for all vertical levels, and will divide out of the
      ! equation.  Also the PTOP's will cancel each other out.  Thus
      ! for GEOS-3, the above equation reduces to (EQUATION 3):
      !
      !           [ ( Q(L  ) * ( SIG_EDGE(L  ) - SIG_EDGE(L+1) ) ) + 
      !             ( Q(L+1) * ( SIG_EDGE(L+1) - SIG_EDGE(L+2) ) ) ]
      !  Q(L') = ----------------------------------------------------
      !                     SIG_EDGE(L) - SIG_EDGE(L+2)
      !=================================================================     

      ! For GEOS-3, EDGE_IN are the sigma values at grid box edges
      ! Otherwise,  EDGE_IN are the pressures at grid box edges
      OUT   = ( IN(L  ) * ( EDGE_IN(L  ) - EDGE_IN(L+1) ) ) +
     &        ( IN(L+1) * ( EDGE_IN(L+1) - EDGE_IN(L+2) ) )

      ! Divde by sigma thickness of new thick level
      OUT   = OUT / ( EDGE_IN(L) - EDGE_IN(L+2) )
       
      END FUNCTION LUMP_2_R4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lump_2_r8
!
! !DESCRIPTION: Function LUMP\_2\_R8 lumps 2 sigma levels into one thick 
!  level.  Input arguments must be REAL*8. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION LUMP_2_R8( IN, L_IN, L ) RESULT( OUT )
!
! !USES:
!
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: IN(L_IN)   ! Column of data on input grid
      INTEGER, INTENT(IN) :: L_IN       ! Vertical dimension of the IN array
      INTEGER, INTENT(IN) :: L          ! Level on input grid from which 
                                        !  to start regridding
!
! !RETURN VALUE:
!
      REAL*8              :: OUT        ! Data on output grid: 2 lumped levels
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Renamed SIGE_IN to EDGE_IN to denote that it is not always a sigma
!        coordinate (as for GEOS-4).  Also updated comments (bmy, 10/31/03)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! LUMP_2_R8 begins here!
      !=================================================================      

      ! Error check: prevent array out of bounds error
      IF ( L < 1 .or. L > L_IN .or. L+2 > L_IN+1 ) THEN 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'ERROR: L < 1 or L > L_IN or L+2 > L_IN+1!'
         WRITE( 6, '(a)' ) 'STOP in LUMP_2 ("regrid_mod.f")'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! When lumping the levels together, we need to perform a weighted
      ! average by air mass.  The air mass in a grid box is given by:
      !
      !  Air Mass = Delta-P [hPa] * 100/g * Grid box sfc area [cm2]
      !    
      ! Where Delta-P is the difference in pressure between the bottom 
      ! and top edges of the grid box (Delta-P is positive), 100/g is a 
      ! constant, and the grid box surface area is also constant w/in 
      ! the same vertical column.  Therefore, for a vertical column,
      ! the air mass in a grid box really only depends on Delta-P.
      !
      ! Because of this, we may compute the quantity Q(L') on the new 
      ! merged sigma according to the weighted average (EQUATION 1):
      !
      !             [ ( Q(L  ) * ( PEDGE(L  ) - PEDGE(L+1) ) ) + 
      !               ( Q(L+1) * ( PEDGE(L+1) - PEDGE(L+2) ) ) +
      !               ( Q(L+2) * ( PEDGE(L+2) - PEDGE(L+3) ) ) +
      !               ( Q(L+3) * ( PEDGE(L+3) - PEDGE(L+4) ) ) ]
      !  Q(L') = ------------------------------------------------
      !                       PEDGE(L) - PEDGE(L+4)
      !
      ! where PEDGE(L) is the pressure at the bottom edge of layer L.
      !
      ! GEOS-4/GEOS-5/MERRA are a hybrid sigma-pressure grid, with 
      ! all of the levels above level a certain level being pure 
      ! pressure levels.  Therefore, for these grids, we may just 
      ! use EQUATION 1 exactly as written above.
      !
      ! However, GEOS-3 is a pure sigma grid.  The pressure at the 
      ! edge of a grid box is given by EQUATION 2:
      ! 
      !  PEDGE(I,J,L) = PTOP + ( SIG_EDGE(L) * ( Psurf(I,J) - PTOP) )
      !
      ! In a vertical column, then ( Psurf(I,J) - PTOP ) will be the 
      ! same for all vertical levels, and will divide out of the
      ! equation.  Also the PTOP's will cancel each other out.  Thus
      ! for GEOS-3, the above equation reduces to (EQUATION 3):
      !
      !           [ ( Q(L  ) * ( SIG_EDGE(L  ) - SIG_EDGE(L+1) ) ) + 
      !             ( Q(L+1) * ( SIG_EDGE(L+1) - SIG_EDGE(L+2) ) ) ]
      !  Q(L') = ----------------------------------------------------
      !                     SIG_EDGE(L) - SIG_EDGE(L+2)
      !=================================================================     

      ! For GEOS-3, EDGE_IN are the sigma values at grid box edges
      ! Othewise,   EDGE_IN are the pressures at grid box edges
      OUT   = ( IN(L  ) * ( EDGE_IN(L  ) - EDGE_IN(L+1) ) ) +
     &        ( IN(L+1) * ( EDGE_IN(L+1) - EDGE_IN(L+2) ) )

      ! Divde by thickness of new lumped level
      OUT   = OUT / ( EDGE_IN(L) - EDGE_IN(L+2) )
       
      END FUNCTION LUMP_2_R8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lump_4_r4
!
! !DESCRIPTION: Function LUMP\_4\_R4 lumps 4 sigma levels into one thick 
!  level.  Input arguments must be REAL*4.
!\\
!\\
! !INTERFACE:
!
      FUNCTION LUMP_4_R4( IN, L_IN, L ) RESULT( OUT )
!
! !USES:
!
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN) :: IN(L_IN)   ! Column of data on input grid
      INTEGER, INTENT(IN) :: L_IN       ! Vertical dimension of the IN array
      INTEGER, INTENT(IN) :: L          ! Level on input grid from which 
                                        !  to start regridding
!
! !RETURN VALUE:
!
      REAL*4              :: OUT        ! Data on output grid: 4 lumped levels
!
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Renamed SIGE_IN to EDGE_IN to denote that it is not always a sigma
!        coordinate (as for GEOS-4).  Also updated comments (bmy, 10/31/03)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! LUMP_4_R4 begins here!
      !=================================================================      

      ! Error check: prevent array out of bounds error
      IF ( L < 1 .or. L > L_IN .or. L+4 > L_IN+1 ) THEN 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 ) 
         WRITE( 6, '(a)' ) 'ERROR: L < 1 or L > L_IN or L+4 > L_IN+1!'
         WRITE( 6, '(a)' ) 'STOP in LUMP_4 ("regrid_mod.f")'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 ) 
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! When lumping the levels together, we need to perform a weighted
      ! average by air mass.  The air mass in a grid box is given by:
      !
      !  Air Mass = Delta-P [hPa] * 100/g * Grid box sfc area [cm2]
      !    
      ! Where Delta-P is the difference in pressure between the bottom 
      ! and top edges of the grid box (Delta-P is positive), 100/g is a 
      ! constant, and the grid box surface area is also constant w/in 
      ! the same vertical column.  Therefore, for a vertical column,
      ! the air mass in a grid box really only depends on Delta-P.
      !
      ! Because of this, we may compute the quantity Q(L') on the new 
      ! merged sigma according to the weighted average (EQUATION 1):
      !
      !             [ ( Q(L  ) * ( PEDGE(L  ) - PEDGE(L+1) ) ) + 
      !               ( Q(L+1) * ( PEDGE(L+1) - PEDGE(L+2) ) ) +
      !               ( Q(L+2) * ( PEDGE(L+2) - PEDGE(L+3) ) ) +
      !               ( Q(L+3) * ( PEDGE(L+3) - PEDGE(L+4) ) ) ]
      !  Q(L') = --------------------------------------------------
      !                       PEDGE(L) - PEDGE(L+4)
      !
      ! where PEDGE(L) is the pressure at the bottom edge of layer L.
      !
      ! GEOS-4/GEOS-5/MERRA are a hybrid sigma-pressure grid, with 
      ! all of the levels above level a certain level being pure 
      ! pressure levels.  Therefore, for these grids, we may just 
      ! use EQUATION 1 exactly as written above.
      !
      ! However, GEOS-3 is a pure sigma grid.  The pressure at the 
      ! edge of a grid box is given by EQUATION 2:
      ! 
      !  PEDGE(I,J,L) = PTOP + ( SIG_EDGE(L) * ( Psurf(I,J) - PTOP) )
      !
      ! In a vertical column, then ( Psurf(I,J) - PTOP ) will be the 
      ! same for all vertical levels, and will divide out of the
      ! equation.  Also the PTOP's will cancel each other out.  Thus
      ! for GEOS-3, the above equation reduces to (EQUATION 3):
      !
      !           [ ( Q(L  ) * ( SIG_EDGE(L  ) - SIG_EDGE(L+1) ) ) + 
      !             ( Q(L+1) * ( SIG_EDGE(L+1) - SIG_EDGE(L+2) ) ) +
      !             ( Q(L+2) * ( SIG_EDGE(L+2) - SIG_EDGE(L+3) ) ) +
      !             ( Q(L+3) * ( SIG_EDGE(L+3) - SIG_EDGE(L+4) ) ) ]
      !  Q(L') = ----------------------------------------------------
      !                     SIG_EDGE(L) - SIG_EDGE(L+4)
      !=================================================================     

      ! For GEOS-3, EDGE_IN are the sigma values at grid box edges
      ! Otherwise,  EDGE_IN are the pressures at grid box edges
      OUT   = ( IN(L  ) * ( EDGE_IN(L  ) - EDGE_IN(L+1) ) ) +
     &        ( IN(L+1) * ( EDGE_IN(L+1) - EDGE_IN(L+2) ) ) +
     &        ( IN(L+2) * ( EDGE_IN(L+2) - EDGE_IN(L+3) ) ) +
     &        ( IN(L+3) * ( EDGE_IN(L+3) - EDGE_IN(L+4) ) ) 

      ! Divde by thickness of new lumped level
      OUT   = OUT / ( EDGE_IN(L) - EDGE_IN(L+4) )
       
      END FUNCTION LUMP_4_R4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lump_4_r8
!
! !DESCRIPTION: Function LUMP\_4\_R8 lumps 4 sigma levels into one thick 
!  level.  Input arguments must be REAL*8.
!\\
!\\
! !INTERFACE:
!
      FUNCTION LUMP_4_R8( IN, L_IN, L ) RESULT( OUT )
!
! !USES:
!
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: IN(L_IN)   ! Column of data on input grid
      INTEGER, INTENT(IN) :: L_IN       ! Vertical dimension of the IN array
      INTEGER, INTENT(IN) :: L          ! Level on input grid from which 
                                        !  to start regridding
!
! !RETURN VALUE:
!
      REAL*8              :: OUT        ! Data on output grid: 4 lumped levels
! 
! !REVISION HISTORY: 
!  19 Sep 2001 - R. Yantosca - Initial version
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Renamed SIGE_IN to EDGE_IN to denote that it is not always a sigma
!        coordinate (as for GEOS-4).  Also updated comments (bmy, 10/31/03)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! LUMP_4_R8 begins here!
      !=================================================================      

      ! Error check: prevent array out of bounds error
      IF ( L < 1 .or. L > L_IN .or. L+4 > L_IN+1 ) THEN 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 ) 
         WRITE( 6, '(a)' ) 'ERROR: L < 1 or L > L_IN or L+4 > L_IN+1!'
         WRITE( 6, '(a)' ) 'STOP in LUMP_4 ("regrid_mod.f")'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 ) 
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! When lumping the levels together, we need to perform a weighted
      ! average by air mass.  The air mass in a grid box is given by:
      !
      !  Air Mass = Delta-P [hPa] * 100/g * Grid box sfc area [cm2]
      !    
      ! Where Delta-P is the difference in pressure between the bottom 
      ! and top edges of the grid box (Delta-P is positive), 100/g is a 
      ! constant, and the grid box surface area is also constant w/in 
      ! the same vertical column.  Therefore, for a vertical column,
      ! the air mass in a grid box really only depends on Delta-P.
      !
      ! Because of this, we may compute the quantity Q(L') on the new 
      ! merged sigma according to the weighted average (EQUATION 1):
      !
      !             [ ( Q(L  ) * ( PEDGE(L  ) - PEDGE(L+1) ) ) + 
      !               ( Q(L+1) * ( PEDGE(L+1) - PEDGE(L+2) ) ) +
      !               ( Q(L+2) * ( PEDGE(L+2) - PEDGE(L+3) ) ) +
      !               ( Q(L+3) * ( PEDGE(L+3) - PEDGE(L+4) ) ) ]
      !  Q(L') = ------------------------------------------------
      !                       PEDGE(L) - PEDGE(L+4)
      !
      ! where PEDGE(L) is the pressure at the bottom edge of layer L.
      !
      ! GEOS-4/GEOS-5/MERRA are a hybrid sigma-pressure grid, with 
      ! all of the levels above level a certain level being pure 
      ! pressure levels.  Therefore, for these grids, we may just 
      ! use EQUATION 1 exactly as written above.
      !
      ! However, GEOS-3 is a pure sigma grid.  The pressure at the 
      ! edge of a grid box is given by EQUATION 2:
      ! 
      !  PEDGE(I,J,L) = PTOP + ( SIG_EDGE(L) * ( Psurf(I,J) - PTOP) )
      !
      ! In a vertical column, then ( Psurf(I,J) - PTOP ) will be the 
      ! same for all vertical levels, and will divide out of the
      ! equation.  Also the PTOP's will cancel each other out.  Thus
      ! for GEOS-3, the above equation reduces to (EQUATION 3):
      !
      !           [ ( Q(L  ) * ( SIG_EDGE(L  ) - SIG_EDGE(L+1) ) ) + 
      !             ( Q(L+1) * ( SIG_EDGE(L+1) - SIG_EDGE(L+2) ) ) +
      !             ( Q(L+2) * ( SIG_EDGE(L+2) - SIG_EDGE(L+3) ) ) +
      !             ( Q(L+3) * ( SIG_EDGE(L+3) - SIG_EDGE(L+4) ) ) ]
      !  Q(L') = ------------------------------------------------
      !                     SIG_EDGE(L) - SIG_EDGE(L+4)
      !================================================================= 

      ! For GEOS-3, EDGE_IN are the sigma values at grid box edges
      ! Otherwise,  EDGE_IN are the pressures at grid box edges
      OUT   = ( IN(L  ) * ( EDGE_IN(L  ) - EDGE_IN(L+1) ) ) +
     &        ( IN(L+1) * ( EDGE_IN(L+1) - EDGE_IN(L+2) ) ) +
     &        ( IN(L+2) * ( EDGE_IN(L+2) - EDGE_IN(L+3) ) ) +
     &        ( IN(L+3) * ( EDGE_IN(L+3) - EDGE_IN(L+4) ) ) 

      ! Divde by thickness of new lumped level
      OUT   = OUT / ( EDGE_IN(L) - EDGE_IN(L+4) )
       
      END FUNCTION LUMP_4_R8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_transfer
!
! !DESCRIPTION: Subroutine INIT\_TRANSFER initializes and zeroes 
!  all module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TRANSFER( THIS_I0, THIS_J0 )
!
! !USES:
!
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: THIS_I0    ! Global X (longitude) offset
      INTEGER, INTENT(IN) :: THIS_J0    ! Global Y (latitude)  offset
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
      IF ( LLPAR == LGLOB ) THEN

         ! Full vertical resolution; copy all levels!
         L_COPY = LGLOB 

      ELSE

#if   defined( GEOS_4  )
         L_COPY = 19       ! GEOS-4: Copy up to L=19
#elif defined( GEOS_5  )
         L_COPY = 36       ! GEOS-5: Copy up to L=36
#elif defined( GEOS_FP )
         L_COPY = 36       ! GEOS-FP: Copy up to L=36
#elif defined( MERRA   )
         L_COPY = 36       ! MERRA: Copy up to L=36
#elif defined( GCAP    )
         L_COPY = LGLOB    ! GCAP: Copy all levels
#endif

      ENDIF

      !=================================================================      
      ! Define vertical edges for collapsing stratospheric levels
      !=================================================================

      ! Allocate the EDGE_IN array
      ALLOCATE( EDGE_IN( LGLOB + 1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EDGE_IN' )
      EDGE_IN = 0d0

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

      !-----------------------------------------------------------------
      ! For GEOS-5/MERRA, levels 1-31 are "terrain-following" 
      ! coordinates (i.e. vary with location), and levels 32-72 are 
      ! fixed pressure levels.  The transition pressure is 176.93 hPa, 
      ! which is the edge between L=31 and L=32.
      !
      ! Initialize EDGE_IN with the original 73 Ap values for GEOS-5.
      !-----------------------------------------------------------------
      EDGE_IN = (/ 
     &        0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01,
     &        1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01,
     &        4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01,
     &        7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01,
     &        1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02,
     &        1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02,
     &        2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02,
     &        2.243630d+02, 2.168650d+02, 2.011920d+02, 
!------- EDGES OF GEOS-5 FIXED PRESSURE LEVELS OCCUR BELOW THIS LINE ------
     &                                                  1.769300d+02,
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

#elif defined( GEOS_4 )

      !-----------------------------------------------------------------
      ! For GEOS-4, levels 1-14 are "terrain-following" coordinates
      ! (i.e. vary with location), and levels 15-55 are fixed pressure 
      ! levels.  The transition pressure is 176.93 hPa, which is the
      ! edge between L=14 and L=15.
      !
      ! Initialize EDGE_IN with the original 56 Ap values for GEOS-4.
      !-----------------------------------------------------------------
      EDGE_IN = (/  0.000000d0,   0.000000d0,  12.704939d0,  
     &             35.465965d0,  66.098427d0, 101.671654d0, 
     &            138.744400d0, 173.403183d0, 198.737839d0, 
     &            215.417526d0, 223.884689d0, 224.362869d0,
     &            216.864929d0, 201.192093d0, 176.929993d0, 
     &            150.393005d0, 127.837006d0, 108.663429d0,  
     &             92.365662d0,  78.512299d0,  66.603378d0,  
     &             56.387939d0,  47.643932d0,  40.175419d0, 
     &             33.809956d0,  28.367815d0,  23.730362d0,  
     &             19.791553d0,  16.457071d0,  13.643393d0,  
     &             11.276889d0,   9.292943d0,   7.619839d0,   
     &              6.216800d0,   5.046805d0,   4.076567d0, 
     &              3.276433d0,   2.620212d0,   2.084972d0,   
     &              1.650792d0,   1.300508d0,   1.019442d0,   
     &              0.795134d0,   0.616779d0,   0.475806d0,   
     &              0.365041d0,   0.278526d0,   0.211349d0, 
     &              0.159495d0,   0.119703d0,   0.089345d0,   
     &              0.066000d0,   0.047585d0,   0.032700d0,   
     &              0.020000d0,   0.010000d0 /)

#endif

      ! We have now initialized everything
      IS_INIT = .TRUE.

      END SUBROUTINE INIT_TRANSFER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_transfer
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
