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
      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
      USE PRECISION_MOD


      IMPLICIT NONE

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: TRANSFER_3D
      PUBLIC  :: TRANSFER_3D_Lp1
      PUBLIC  :: TRANSFER_3D_yan
      PUBLIC  :: TRANSFER_G5_PLE
      PUBLIC  :: INIT_TRANSFER
      PUBLIC  :: CLEANUP_TRANSFER

      INTERFACE TRANSFER_3D
         MODULE PROCEDURE TRANSFER_3D_R4
         MODULE PROCEDURE TRANSFER_3D_R8
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
      PRIVATE :: TRANSFER_3D_R4
      PRIVATE :: TRANSFER_3D_R8

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
! !IROUTINE: Transfer_3d_r4
!
! !DESCRIPTION: Subroutine TRANSFER\_3D\_R8 transfers 3-dimensional data from a 
!  REAL*4  array to a REAL*4 array.  Vertical layers are collapsed (from LGLOB
!  to LLPAR) if necessary.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_3D_R4( IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR,LGLOB)    ! Input data
!
! !OUTPUT PARAMETERS:
!
      REAL*4,  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLPAR)   ! Output data
! 
! !REVISION HISTORY: 
!  03 Apr 2014 - R. Yantosca - Initial version, based on TRANSFER_3D_R8
!  11 Aug 2015 - R. Yantosca - MERRA2 behaves as GEOS-5, MERRA, GEOS-FP
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

         !--------------------------------------------------------------
         ! Lump 72 levels into 47 levels, starting above 
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

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE TRANSFER_3D_R4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Transfer_3d_r8
!
! !DESCRIPTION: Subroutine TRANSFER\_3D\_R8 transfers 3-dimensional data from a 
!  REAL*4  array to a REAL(fp) array.  Vertical layers are collapsed (from LGLOB
!  to LLPAR) if necessary.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TRANSFER_3D_R8( IN, OUT )
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
!  03 Apr 2014 - R. Yantosca - Renamed to TRANSFER_3D_R8 so that it can
!                              be overloaded with an interface
!  11 Aug 2015 - R. Yantosca - MERRA2 behaves as GEOS-5, MERRA, GEOS-FP
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

         !--------------------------------------------------------------
         ! Lump 72 levels into 47 levels, starting above 
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

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE TRANSFER_3D_R8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Transfer_G5_Ple
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
      REAL(fp),  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLPAR+1)   ! Output data
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

      END SUBROUTINE TRANSFER_G5_PLE
      
      SUBROUTINE TRANSFER_3D_yan( NI,NJ,NK, IN, OUT )
!
! !INPUT PARAMETERS: 
!
      REAL*4,  INTENT(IN)  :: IN(NI,NJ,NK)    ! Input data
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
      INTEGER,          INTENT(IN)     :: NI, NJ, NK
      INTEGER   :: I,J,K
     
      !=================================================================
      ! TRANSFER_3D_Lp1 begins here!
      !=================================================================

      ! Copy the first L_COPY+1 levels

      OUT(:,:,:) = IN(:,:,:)
      


      END SUBROUTINE TRANSFER_3D_yan
!EOC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Transfer_3d_Lp1
!
! !DESCRIPTION: Subroutine TRANSFER\_3D\_Lp1 transfers 3-D data from a REAL*4 
!  array of dimension (IIPAR,JJPAR,LGLOB+1) to a REAL(fp) array of dimension 
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
      REAL(fp),  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLPAR+1)   ! Output data
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

      END SUBROUTINE TRANSFER_3D_Lp1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Lump_2_r4
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
! !IROUTINE: Lump_2_r8
!
! !DESCRIPTION: Function LUMP\_2\_R8 lumps 2 sigma levels into one thick 
!  level.  Input arguments must be REAL(fp). 
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
! !IROUTINE: Lump_4_r4
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
!  level.  Input arguments must be REAL(fp).
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
! !IROUTINE: Init_Transfer
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
      IF ( LLPAR == LGLOB ) THEN

         ! Full vertical resolution; copy all levels!
         L_COPY = LGLOB 

      ELSE

         ! Copy up to L=36 (GEOS-FP, MERRA-2)
         L_COPY = 36

      ENDIF

      !=================================================================      
      ! Define vertical edges for collapsing stratospheric levels
      !=================================================================

      ! Allocate the EDGE_IN array
      ALLOCATE( EDGE_IN( LGLOB + 1 ), STAT=AS )
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
