! $Id: transfer_mod.f,v 1.1 2009/11/20 21:43:02 bmy Exp $
      MODULE TRANSFER_MOD
!
!******************************************************************************
!  Module TRANSFER_MOD contains routines used to copy data from REAL*4 to
!  REAL*8 arrays after being read from disk.  Also, vertical levels will be
!  collapsed in the stratosphere if necessary.  This will help us to gain 
!  computational advantage. (mje, bmy, 9/27/01, 10/3/07)
!
!  NOTE: The level above which we start collapsing layers is ~78 hPa.  
!
!  Module Variables:
!  ============================================================================
!  (1 ) EDGE_IN           : Input sigma edges (for pure sigma models)
!  (2 ) I0                : Global longitude offset (%%% NOTE: usually=0 %%%)
!  (3 ) J0                : Global latitude  offset (%%% NOTE: usually=0 %%%)
!  (4 ) L_COPY            : # of levels to copy (before stratosphere lumping)
!
!  Module Routines:
!  ============================================================================
!  (1 ) TRANSFER_A6       : Transfers GEOS A-6   fields, regrids if necessary
!  (2 ) TRANSFER_3D       : Transfers GEOS 3-D   fields, regrids if necessary
!  (3 ) TRANSFER_3D_TROP  : Transfers GEOS 3-D   fields up to tropopause level
!  (4 ) TRANSFER_G5_PLE   : Transfers GEOS-5 3-D pressure edges, regrids
!  (5 ) TRANSFER_3D_Lp1   : Transfers GEOS-5 3-D fields defined on level edges
!  (6 ) TRANSFER_ZONAL_R4 : Transfers GEOS zonal fields, regrids (REAL*4)
!  (7 ) TRANSFER_ZONAL_R8 : Transfers GEOS zonal fields, regrids (REAL*8) 
!  (8 ) TRANSFER_ZONAL    : Transfers GEOS zonal fields, regrids if necessary 
!  (9 ) TRANSFER_2D_INT   : Transfers GEOS 2-D   fields (INTEGER argument)
!  (10) TRANSFER_2D_R4    : Transfers GEOS 2-D   fields (REAL*4 argument)
!  (11) TRANSFER_2D_R8    : Transfers GEOS 2-D   fields (REAL*8 argument)
!  (12) TRANSFER_TO_1D    : Transfers GEOS 2-D   fields to a 1-D array
!  (13) LUMP_2_R4         : Combines 2 levels into 1 thick level (REAL*4) 
!  (14) LUMP_2_R8         : Combines 2 levels into 1 thick level (REAL*8) 
!  (15) LUMP_4_R4         : Combines 4 levels into 1 thick level (REAL*4)
!  (16) LUMP_4_R8         : Combines 4 levels into 1 thick level (REAL*8)
!  (17) INIT_TRANSFER     : Allocates and initializes the EDGE_IN array
!  (18) CLEANUP_TRANSFER  : Deallocates the EDGE_IN array
!
!  Module Interfaces:
!  ============================================================================
!  (1 ) LUMP_2            : Overloads LUMP_2_*         module routines
!  (2 ) LUMP_4            : Overloads LUMP_4_*         module routines
!  (3 ) TRANSFER_2D       : Overloads TRANSFER_2D_*    module routines
!  (4 ) TRANSFER_ZONAL    : Overloads TRANSFER_ZONAL_* module routines
!
!  GEOS-Chem modules referenced by "transfer_mod.f"
!  ============================================================================
!  (1 ) error_mod.f      : Module w/ NaN and other error check routines
!
!  NOTES:
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
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "transfer_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: TRANSFER_A6
      PUBLIC :: TRANSFER_2D
      PUBLIC :: TRANSFER_3D
      PUBLIC :: TRANSFER_3D_Lp1
      PUBLIC :: TRANSFER_3D_TROP
      PUBLIC :: TRANSFER_G5_PLE
      PUBLIC :: TRANSFER_ZONAL
      PUBLIC :: TRANSFER_TO_1D
      PUBLIC :: INIT_TRANSFER
      PUBLIC :: CLEANUP_TRANSFER

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      INTEGER             :: I0
      INTEGER             :: J0
      INTEGER             :: L_COPY

      ! Arrays
      REAL*8, ALLOCATABLE :: EDGE_IN(:)

      !=================================================================
      ! MODULE INTERFACES -- "bind" two or more routines with different
      ! argument types or # of arguments under one unique name
      !================================================================= 

      ! Interface for routines to lump 2 levels together
      INTERFACE LUMP_2
         MODULE PROCEDURE LUMP_2_R4
         MODULE PROCEDURE LUMP_2_R8
      END INTERFACE

      ! Interface for routines to lump 4 levels together
      INTERFACE LUMP_4
         MODULE PROCEDURE LUMP_4_R4
         MODULE PROCEDURE LUMP_4_R8
      END INTERFACE

      ! Interface for routines which copy 2-D data 
      INTERFACE TRANSFER_2D
         MODULE PROCEDURE TRANSFER_2D_INT
         MODULE PROCEDURE TRANSFER_2D_R4
         MODULE PROCEDURE TRANSFER_2D_R8
      END INTERFACE

      ! Interface for routines which copy zonal data 
      INTERFACE TRANSFER_ZONAL
         MODULE PROCEDURE TRANSFER_ZONAL_R4
         MODULE PROCEDURE TRANSFER_ZONAL_R8
      END INTERFACE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_A6( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_A6 transfers A-6 data from a REAL*4 array to a REAL*8
!  array.  Vertical layers are collapsed (from LGLOB to LLPAR) if necessary.
!  (mje, bmy, 9/21/01, 11/6/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN  (REAL*4) : Input field,  dimensioned (IGLOB,JGLOB,LGLOB)
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) OUT (REAL*4) : Output field, dimensioned (LLPAR,IIPAR,JJPAR)
!
!  NOTES:
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
!******************************************************************************
!
#     include "CMN_SIZE"    ! Size parameters

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IGLOB,JGLOB,LGLOB)
      REAL*8,  INTENT(OUT) :: OUT(LLPAR,IIPAR,JJPAR)

      ! Local variables
      INTEGER              :: I, J, L
      REAL*4               :: INCOL(LGLOB)
      
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
            
#if   defined( GEOS_3 ) 

         !--------------------------------------------------------------
         ! GEOS-3: Lump 48 levels into 30 levels, starting above L=22.  
         ! Lump levels in groups of 2, then 4. (cf. Mat Evans)
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time, starting at L=23
         OUT(23,I,J) = LUMP_2( INCOL, LGLOB, 23 )
         OUT(24,I,J) = LUMP_2( INCOL, LGLOB, 25 )
         OUT(25,I,J) = LUMP_2( INCOL, LGLOB, 27 )

         ! Lump 4 levels together at a time, starting at L=29
         OUT(26,I,J) = LUMP_4( INCOL, LGLOB, 29 )
         OUT(27,I,J) = LUMP_4( INCOL, LGLOB, 33 )
         OUT(28,I,J) = LUMP_4( INCOL, LGLOB, 37 ) 
         OUT(29,I,J) = LUMP_4( INCOL, LGLOB, 41 ) 
         OUT(30,I,J) = LUMP_4( INCOL, LGLOB, 45 ) 

#elif defined( GEOS_4 )

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

#elif defined( GEOS_5 )

         !--------------------------------------------------------------
         ! GEOS-5: Lump 72 levels into 47 levels, starting above L=36
         ! Lump levels in groups of 2, then 4. (cf. Bob Yantosca)
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

      ! Return to calling program
      END SUBROUTINE TRANSFER_A6

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_3D( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_3D transfers A-6 data from a REAL*4 array to a REAL*8
!  array.  Vertical layers are collapsed (from LGLOB to LLPAR) if necessary.
!  (mje, bmy, 9/21/01, 2/8/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN  (REAL*4) : Input field,  of dimension (IGLOB,JGLOB,LGLOB)
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) OUT (REAL*8) : Output field, of dimension (IIPAR,JJPAR,LLPAR)
!
!  NOTES:
!  (1 ) Lump levels together in groups of 2 or 4, as dictated by Mat Evans.
!        (bmy, 9/21/01)
!  (2 ) Assumes that LLPAR == LGLOB for GEOS-1, GEOS-STRAT (bmy, 9/21/01)
!  (3 ) Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        Now I0, J0 are local variables. (bmy, 3/11/03)
!  (4 ) Added code to regrid GEOS-4 from 55 --> 30 levels (mje, bmy, 10/31/03)
!  (5 ) Now modified for GEOS-5 met fields (bmy, 5/24/05)
!  (6 ) Rewritten for clarity (bmy, 2/8/07)
!******************************************************************************
!
#     include "CMN_SIZE"    ! Size parameters

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR,LGLOB)
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER              :: I, J
      REAL*4               :: INCOL(LGLOB)
     
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

#if   defined( GEOS_3 ) 

         !--------------------------------------------------------------
         ! GEOS-3: Lump 48 levels into 30 levels, starting above L=22.  
         ! Lump levels in groups of 2, then 4. (cf. Mat Evans)
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time, starting at L=23
         OUT(I,J,23) = LUMP_2( INCOL, LGLOB, 23 )
         OUT(I,J,24) = LUMP_2( INCOL, LGLOB, 25 )
         OUT(I,J,25) = LUMP_2( INCOL, LGLOB, 27 )

         ! Lump 4 levels together at a time, starting at L=29
         OUT(I,J,26) = LUMP_4( INCOL, LGLOB, 29 )
         OUT(I,J,27) = LUMP_4( INCOL, LGLOB, 33 )
         OUT(I,J,28) = LUMP_4( INCOL, LGLOB, 37 ) 
         OUT(I,J,29) = LUMP_4( INCOL, LGLOB, 41 ) 
         OUT(I,J,30) = LUMP_4( INCOL, LGLOB, 45 ) 

#elif defined( GEOS_4 )

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

#elif defined( GEOS_5 )

         !--------------------------------------------------------------
         ! GEOS-5: Lump 72 levels into 47 levels, starting above L=36
         ! Lump levels in groups of 2, then 4. (cf. Bob Yantosca)
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

      ! Return to calling program
      END SUBROUTINE TRANSFER_3D

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_G5_PLE( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_G5_PLE transfers GEOS-5 pressure edge data from the
!  native 72-level grid to the reduced 47-level grid.  (bmy, 2/8/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN  (REAL*4) : Input field,  of dimension (IGLOB,JGLOB,LGLOB)
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) OUT (REAL*8) : Output field, of dimension (IIPAR,JJPAR,LLPAR)
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE"    ! Size parameters

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR,LGLOB+1)
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLPAR+1)

      ! Local variables
      INTEGER              :: I, J
     
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

#if   defined( GEOS_5 )

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

      ! Return to calling program
      END SUBROUTINE TRANSFER_G5_PLE

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_3D_Lp1( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_3D_Lp1 transfers 3-D data from a REAL*4 array of 
!  dimension (IGLOB,JGLOB,LGLOB+1) to a REAL*8 array of dimension 
!  (IIPAR,JJPAR,LLPAR+1).  Regrid in the vertical if needed.
!  (bmy, 9/21/01, 2/8/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN  (REAL*4) : Input field,  of dimension (IGLOB,JGLOB,LGLOB+1)
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) OUT (REAL*8) : Output field, of dimension (IIPAR,JJPAR,LLPAR+1)
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE"

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IGLOB,JGLOB,LGLOB+1)
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLPAR+1)

      ! Local variables
      LOGICAL, SAVE        :: FIRST = .TRUE.  
      INTEGER              :: I, J
      REAL*4               :: INCOL(LGLOB)
     
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
      
#if   defined( GEOS_5 )

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

      ! Return to calling program
      END SUBROUTINE TRANSFER_3D_Lp1

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_3D_TROP( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_3D_TROP transfers tropospheric 3-D data from a REAL*4 
!  array to a REAL*8 array. (mje, bmy, 9/21/01, 2/8/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN  (REAL*4) : Input field,  of dimension (IGLOB,JGLOB,LLTROP)
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) OUT (REAL*8) : Output field, of dimension (IIPAR,LLPAR,LLTROP)
!
!  NOTES:
!  (1 ) Now use LLTROP_FIX instead of LLTROP, since most of the offline
!        simulations use the annual mean tropopause (bmy, 2/8/07)
!******************************************************************************
!
#     include "CMN_SIZE"    ! Size parameters

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IGLOB,JGLOB,LLTROP_FIX)
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLTROP_FIX)

      ! Local variables
      INTEGER              :: L
     
      !=================================================================
      ! TRANSFER_3D_TROP
      !=================================================================

      ! Cast to REAL*8 abd resize up to LLTROP
      DO L = 1, LLTROP_FIX
         CALL TRANSFER_2D( IN(:,:,L), OUT(:,:,L) )
      ENDDO
      
      ! Return to calling program
      END SUBROUTINE TRANSFER_3D_TROP 

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_ZONAL_R4( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_ZOJAL_R4 transfers zonal-mean data from a REAL*4 array 
!  to a REAL*8 array.  Vertical levels are collapsed (from LGLOB to LLPAR) if 
!  necessary. (mje, bmy, 9/21/01, 2/8/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN  (REAL*4) : Input  field, of dimension (JGLOB,LGLOB)
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) OUT (REAL*8) : Output field, of dimension (JJPAR,LLPAR)
!
!  NOTES:
!  (1 ) Lump levels together in groups of 2 or 4, as dictated by Mat Evans.
!        (bmy, 9/21/01)
!  (2 ) Assumes that LLPAR == LGLOB for GEOS-1, GEOS-STRAT (bmy, 9/21/01)
!  (3 ) Now use function GET_YOFFSET from "grid_mod.f".  Now I0 and J0 are 
!        local variables (bmy, 3/11/03)
!  (4 ) Added code to regrid GEOS-4 from 55 --> 30 levels (mje, bmy, 10/31/03)
!  (5 ) Rewritten for clarity (bmy, 2/8/07)
!******************************************************************************
!
#     include "CMN_SIZE"    ! Size parameters

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(JGLOB,LGLOB)
      REAL*4,  INTENT(OUT) :: OUT(JJPAR,LLPAR)

      ! Local variables
      INTEGER              :: J
      REAL*4               :: INCOL(LGLOB)
     
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

#if   defined( GEOS_3 )

         !--------------------------------------------------------------
         ! GEOS-3: Lump 48 levels into 30 levels, starting above L=22.  
         ! Lump levels in groups of 2, then 4. (cf. Mat Evans)
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time
         OUT(J,23) = LUMP_2( INCOL, LGLOB, 23 )
         OUT(J,24) = LUMP_2( INCOL, LGLOB, 25 )
         OUT(J,25) = LUMP_2( INCOL, LGLOB, 27 )

         ! Lump 4 levels together at a time
         OUT(J,26) = LUMP_4( INCOL, LGLOB, 29 )
         OUT(J,27) = LUMP_4( INCOL, LGLOB, 33 )
         OUT(J,28) = LUMP_4( INCOL, LGLOB, 37 ) 
         OUT(J,29) = LUMP_4( INCOL, LGLOB, 41 ) 
         OUT(J,30) = LUMP_4( INCOL, LGLOB, 45 ) 

#elif defined( GEOS_4 )
      
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

#elif defined( GEOS_5 )

         !--------------------------------------------------------------
         ! GEOS-5: Lump 72 levels into 47 levels, starting above L=36
         ! Lump levels in groups of 2, then 4.  
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

      ! Return to calling program
      END SUBROUTINE TRANSFER_ZONAL_R4

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_ZONAL_R8( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_ZONAL_R4 transfers zonal mean or lat-alt data from a 
!  REAL*4 array of dimension (JGLOB,LGLOB) to a REAL*8 array of dimension 
!  (JJPAR,LLPAR).   Regrid GEOS-3 data from 48 --> 30 levels or GEOS-4 data
!  from 55 --> 30 levels if necessary. (bmy, 9/21/01, 5/24/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN  (REAL*4) : Input  field, of dimension (JGLOB,LGLOB)
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) OUT (REAL*8) : Output field, of dimension (JJPAR,LLPAR)
!
!  NOTES:
!  (1 ) Lump levels together in groups of 2 or 4, as dictated by Mat Evans.
!        (bmy, 9/21/01)
!  (2 ) Assumes that LLPAR == LGLOB for GEOS-1, GEOS-STRAT (bmy, 9/21/01)
!  (3 ) Now use function GET_YOFFSET from "grid_mod.f".  Now I0 and J0 are 
!        local variables (bmy, 3/11/03)
!  (4 ) Added code to regrid GEOS-4 from 55 --> 30 levels (mje, bmy, 10/31/03)
!  (5 ) Now modified for GEOS-5 met fields (bmy, 5/24/05)
!******************************************************************************
!
#     include "CMN_SIZE"    ! Size parameters

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(JGLOB,LGLOB)
      REAL*8,  INTENT(OUT) :: OUT(JJPAR,LLPAR)

      ! Local variables
      INTEGER              :: J
      REAL*4               :: INCOL(LGLOB)
     
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

#if   defined( GEOS_3 )

         !--------------------------------------------------------------
         ! GEOS-3: Lump 48 levels into 30 levels, starting above L=22.  
         ! Lump levels in groups of 2, then 4. (cf. Mat Evans)
         !--------------------------------------------------------------

         ! Lump 2 levels together at a time
         OUT(J,23) = LUMP_2( INCOL, LGLOB, 23 )
         OUT(J,24) = LUMP_2( INCOL, LGLOB, 25 )
         OUT(J,25) = LUMP_2( INCOL, LGLOB, 27 )

         ! Lump 4 levels together at a time
         OUT(J,26) = LUMP_4( INCOL, LGLOB, 29 )
         OUT(J,27) = LUMP_4( INCOL, LGLOB, 33 )
         OUT(J,28) = LUMP_4( INCOL, LGLOB, 37 ) 
         OUT(J,29) = LUMP_4( INCOL, LGLOB, 41 ) 
         OUT(J,30) = LUMP_4( INCOL, LGLOB, 45 ) 

#elif defined( GEOS_4 )

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

#elif defined( GEOS_5 )

         !--------------------------------------------------------------
         ! GEOS-5: Lump 72 levels into 47 levels, starting above L=36
         ! Lump levels in groups of 2, then 4.  
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

      ! Return to calling program
      END SUBROUTINE TRANSFER_ZONAL_R8

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_2D_INT( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_2D_INT transfers 2-D data from a REAL*4 array of 
!  dimension (IGLOB,JGLOB) to an INTEGER array of dimension (IIPAR,JJPAR). 
!  (bmy, 9/21/01, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN  (REAL*4 ) : Input field,  of dimension (IGLOB,JGLOB)
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) OUT (INTEGER) : Output field, of dimension (IIPAR,JJPAR)
!
!  NOTES:
!  (1 ) Use parallel DO loops to speed things up (bmy, 9/21/01)!
!  (2 ) Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        Now I0 and J0 are local variables. (bmy, 3/11/03)
!******************************************************************************
!
#     include "CMN_SIZE"

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IGLOB,JGLOB)
      INTEGER, INTENT(OUT) :: OUT(IIPAR,JJPAR)

      !=================================================================
      ! TRANSFER_2D_INT begins here!
      !=================================================================

      ! Copy and cast array
      OUT = IN( 1+I0:IIPAR+I0, 1+J0:JJPAR+J0 )

      ! Return to calling program
      END SUBROUTINE TRANSFER_2D_INT

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_2D_R4( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_2D_R4 transfers 2-D data from a REAL*4 array of 
!  dimension (IGLOB,JGLOB) to a REAL*4 array of dimension (IIPAR,JJPAR). 
!  (bmy, 1/25/02, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN  (REAL*4) : Input field,  of dimension (IGLOB,JGLOB)
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) OUT (REAL*4) : Output field, of dimension (IIPAR,JJPAR)
!
!  NOTES:
!  (1 ) Use parallel DO loops to speed things up (bmy, 9/21/01)
!  (2 ) Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f"
!        Now I0 and J0 are local variables (bmy, 3/11/03)
!******************************************************************************
!
#     include "CMN_SIZE"

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IGLOB,JGLOB)
      REAL*4,  INTENT(OUT) :: OUT(IIPAR,JJPAR)

      !=================================================================
      ! TRANSFER_2D_R4 begins here!
      !=================================================================
      
      ! Copy and cast array
      OUT = IN( 1+I0:IIPAR+I0, 1+J0:JJPAR+J0 )

      ! Return to calling program
      END SUBROUTINE TRANSFER_2D_R4

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_2D_R8( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_2D_R8 transfers 2-D data from a REAL*4 array of 
!  dimension (IGLOB,JGLOB) to a REAL*8 array of dimension (IIPAR,JJPAR). 
!  (bmy, 9/21/01, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN  (REAL*4) : Input field,  of dimension (IGLOB,JGLOB)
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) OUT (REAL*8) : Output field, of dimension (IIPAR,JJPAR)
!
!  NOTES:
!  (1 ) Use parallel DO loops to speed things up (bmy, 9/21/01)
!  (2 ) Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f"
!        Now I0 and J0 are local variables. (bmy, 3/11/03)
!******************************************************************************
!
#     include "CMN_SIZE"

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IGLOB,JGLOB)
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR)

      !=================================================================
      ! TRANSFER_2D_R8 begins here!
      !=================================================================

      ! Copy and cast array
      OUT = IN( 1+I0:IIPAR+I0, 1+J0:JJPAR+J0 )
      
      ! Return to calling program
      END SUBROUTINE TRANSFER_2D_R8

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_TO_1D( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_TO_1D transfers 2-D data from a REAL*4 array of 
!  dimension (IGLOB,JGLOB) to 1-D a REAL*8 array of dimension (MAXIJ),
!  where MAXIJ = IIPAR * JJPAR. (bmy, 9/21/01, 3/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN  (REAL*4) : Input field,  of dimension (IGLOB,JGLOB)
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) OUT (REAL*8) : Output field, of dimension (MAXIJ)
!
!  NOTES:
!  (1 ) Use single-processor DO-loops for now (bmy, 9/21/01)
!  (2 ) Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        Now I0 and J0 are local variables. (bmy, 3/11/03)
!******************************************************************************
!
#     include "CMN_SIZE"

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IGLOB,JGLOB)
      REAL*8,  INTENT(OUT) :: OUT(MAXIJ)

      ! Local variables
      INTEGER              :: I, IREF, J, JREF, IJLOOP

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
      
      ! Return to calling program
      END SUBROUTINE TRANSFER_TO_1D

!------------------------------------------------------------------------------

      FUNCTION LUMP_2_R4( IN, L_IN, L ) RESULT( OUT )
!
!******************************************************************************
!  Function LUMP_2_R4 lumps 2 sigma levels into one thick level. 
!  Input arguments must be REAL*4. (bmy, 9/18/01, 10/31/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN   (REAL*4 ) : Column of data on input vertical grid
!  (2 ) L_IN (INTEGER) : Vertical dimension of the IN array
!  (3 ) L    (INTEGER) : Level on input grid from which to start regridding
!
!  Function Value:
!  ============================================================================
!  (4 ) OUT  (REAL*4 ) : Data on output grid -- 2 levels merged together
! 
!  NOTES:
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Renamed SIGE_IN to EDGE_IN to denote that it is not always a sigma
!        coordinate (as for GEOS-4).  Also updated comments (bmy, 10/31/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: L_IN, L
      REAL*4,  INTENT(IN) :: IN(L_IN)

      ! Function value
      REAL*4              :: OUT

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
      ! GEOS-4 is a hybrid sigma-pressure grid, with all of the levels
      ! above level 14 being pure pressure levels.  Therefore, for 
      ! GEOS-4, we may just use EQUATION 1 exactly as written above.
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
      ! For GEOS-4, EDGE_IN are the pressures at grid box edges
      OUT   = ( IN(L  ) * ( EDGE_IN(L  ) - EDGE_IN(L+1) ) ) +
     &        ( IN(L+1) * ( EDGE_IN(L+1) - EDGE_IN(L+2) ) )

      ! Divde by sigma thickness of new thick level
      OUT   = OUT / ( EDGE_IN(L) - EDGE_IN(L+2) )
       
      ! Return to calling routine
      END FUNCTION LUMP_2_R4

!------------------------------------------------------------------------------

      FUNCTION LUMP_2_R8( IN, L_IN, L ) RESULT( OUT )
!
!******************************************************************************
!  Function LUMP_2_R8 lumps 2 sigma levels into one thick level. 
!  Input arguments must be REAL*8. (bmy, 9/18/01, 10/31/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN   (REAL*8 ) : Column of data on input vertical grid
!  (2 ) L_IN (INTEGER) : Vertical dimension of the IN array
!  (3 ) L    (INTEGER) : Level on input grid from which to start regridding
!
!  Function Value:
!  ============================================================================
!  (4 ) OUT  (REAL*8 ) : Data on output grid -- 2 levels merged together
!
!  NOTES:
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Renamed SIGE_IN to EDGE_IN to denote that it is not always a sigma
!        coordinate (as for GEOS-4).  Also updated comments (bmy, 10/31/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: L_IN, L
      REAL*8,  INTENT(IN) :: IN(L_IN)

      ! Function value
      REAL*8              :: OUT

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
      ! GEOS-4 is a hybrid sigma-pressure grid, with all of the levels
      ! above level 14 being pure pressure levels.  Therefore, for 
      ! GEOS-4, we may just use EQUATION 1 exactly as written above.
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
      ! For GEOS-4, EDGE_IN are the pressures at grid box edges
      OUT   = ( IN(L  ) * ( EDGE_IN(L  ) - EDGE_IN(L+1) ) ) +
     &        ( IN(L+1) * ( EDGE_IN(L+1) - EDGE_IN(L+2) ) )

      ! Divde by thickness of new lumped level
      OUT   = OUT / ( EDGE_IN(L) - EDGE_IN(L+2) )
       
      ! Return to calling routine
      END FUNCTION LUMP_2_R8

!------------------------------------------------------------------------------

      FUNCTION LUMP_4_R4( IN, L_IN, L ) RESULT( OUT )
!
!******************************************************************************
!  Function LUMP_4_R4 lumps 4 sigma levels into one thick level. 
!  Input arguments must be REAL*4. (bmy, 9/18/01, 10/31/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN   (REAL*4 ) : Column of data on input vertical grid
!  (2 ) L_IN (INTEGER) : Vertical dimension of the IN array
!  (3 ) L    (INTEGER) : Level on input grid from which to start regridding
!
!  Function Value:
!  ============================================================================
!  (4 ) OUT  (REAL*4 ) : Data on output grid -- 4 levels merged together
!
!  NOTES:
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Renamed SIGE_IN to EDGE_IN to denote that it is not always a sigma
!        coordinate (as for GEOS-4).  Also updated comments (bmy, 10/31/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: L_IN, L
      REAL*4,  INTENT(IN) :: IN(L_IN)

      ! Function value
      REAL*4              :: OUT

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
      ! GEOS-4 is a hybrid sigma-pressure grid, with all of the levels
      ! above level 14 being pure pressure levels.  Therefore, for 
      ! GEOS-4, we may just use EQUATION 1 exactly as written above.
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
      ! For GEOS-4, EDGE_IN are the pressures at grid box edges
      OUT   = ( IN(L  ) * ( EDGE_IN(L  ) - EDGE_IN(L+1) ) ) +
     &        ( IN(L+1) * ( EDGE_IN(L+1) - EDGE_IN(L+2) ) ) +
     &        ( IN(L+2) * ( EDGE_IN(L+2) - EDGE_IN(L+3) ) ) +
     &        ( IN(L+3) * ( EDGE_IN(L+3) - EDGE_IN(L+4) ) ) 

      ! Divde by thickness of new lumped level
      OUT   = OUT / ( EDGE_IN(L) - EDGE_IN(L+4) )
       
      ! Return to calling routine
      END FUNCTION LUMP_4_R4

!------------------------------------------------------------------------------

      FUNCTION LUMP_4_R8( IN, L_IN, L ) RESULT( OUT )
!
!******************************************************************************
!  Function LUMP_4_R8 lumps 4 sigma levels into one thick level. 
!  Input arguments must be REAL*8. (bmy, 9/18/01, 10/31/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN   (REAL*8 ) : Column of data on input vertical grid
!  (2 ) L_IN (INTEGER) : Vertical dimension of the IN array
!  (3 ) L    (INTEGER) : Level on input grid from which to start regridding
!
!  Function Value:
!  ============================================================================
!  (4 ) OUT  (REAL*8 ) : Data on output grid -- 4 levels merged together
!
!  NOTES:
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Renamed SIGE_IN to EDGE_IN to denote that it is not always a sigma
!        coordinate (as for GEOS-4).  Also updated comments (bmy, 10/31/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

#     include "CMN_SIZE"

      ! Arguments
      INTEGER, INTENT(IN) :: L_IN, L
      REAL*8,  INTENT(IN) :: IN(L_IN)

      ! Function value
      REAL*8              :: OUT

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
      ! GEOS-4 is a hybrid sigma-pressure grid, with all of the levels
      ! above level 14 being pure pressure levels.  Therefore, for 
      ! GEOS-4, we may just use EQUATION 1 exactly as written above.
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
      ! For GEOS-4, EDGE_IN are the pressures at grid box edges
      OUT   = ( IN(L  ) * ( EDGE_IN(L  ) - EDGE_IN(L+1) ) ) +
     &        ( IN(L+1) * ( EDGE_IN(L+1) - EDGE_IN(L+2) ) ) +
     &        ( IN(L+2) * ( EDGE_IN(L+2) - EDGE_IN(L+3) ) ) +
     &        ( IN(L+3) * ( EDGE_IN(L+3) - EDGE_IN(L+4) ) ) 

      ! Divde by thickness of new lumped level
      OUT   = OUT / ( EDGE_IN(L) - EDGE_IN(L+4) )
       
      ! Return to calling routine
      END FUNCTION LUMP_4_R8

!------------------------------------------------------------------------------

      SUBROUTINE INIT_TRANSFER( THIS_I0, THIS_J0 )
!
!******************************************************************************
!  Subroutine INIT_TRANSFER initializes and zeroes the EDGE_IN array.
!  (bmy, 9/19/01, 2/8/07)
!
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,  ONLY : ALLOC_ERR

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: THIS_I0, THIS_J0

      ! Local variables
      LOGICAL, SAVE       :: IS_INIT = .FALSE.
      INTEGER             :: AS, L
       
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

#if   defined( GEOS_3 )
         L_COPY = 22       ! GEOS-3: Copy up to L=22
#elif defined( GEOS_4 )
         L_COPY = 19       ! GEOS-4: Copy up to L=19
#elif defined( GEOS_5 )
         L_COPY = 36       ! GEOS-5: Copy up to L=36
#elif defined( GCAP   )
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

#if   defined( GEOS_5 )

      !-----------------------------------------------------------------
      ! For GEOS-5, levels 1-31 are "terrain-following" coordinates
      ! (i.e. vary with location), and levels 32-72 are fixed pressure 
      ! levels.  The transition pressure is 176.93 hPa, which is the
      ! edge between L=31 and L=32.
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

#elif defined( GEOS_3 ) 

      !-----------------------------------------------------------------
      ! For GEOS-3, this is a pure-sigma grid.  
      ! Initialize EDGE_IN with the original 49 sigma edges.
      !-----------------------------------------------------------------
      EDGE_IN = (/ 1.000000d0, 0.997095d0, 0.991200d0, 0.981500d0,    
     &             0.967100d0, 0.946800d0, 0.919500d0, 0.884000d0,    
     &             0.839000d0, 0.783000d0, 0.718200d0, 0.647600d0,    
     &             0.574100d0, 0.500000d0, 0.427800d0, 0.359500d0,    
     &             0.297050d0, 0.241950d0, 0.194640d0, 0.155000d0,    
     &             0.122680d0, 0.096900d0, 0.076480d0, 0.060350d0,   
     &             0.047610d0, 0.037540d0, 0.029600d0, 0.023330d0,   
     &             0.018380d0, 0.014480d0, 0.011405d0, 0.008975d0,  
     &             0.007040d0, 0.005500d0, 0.004280d0, 0.003300d0,  
     &             0.002530d0, 0.001900d0, 0.001440d0, 0.001060d0,  
     &             0.000765d0, 0.000540d0, 0.000370d0, 0.000245d0, 
     &             0.000155d0, 9.20000d-5, 4.75000d-5, 1.76800d-5, 
     &             0.000000d0 /)

#endif

      ! We have now initialized everything
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_TRANSFER
      
!------------------------------------------------------------------------------
      
      SUBROUTINE CLEANUP_TRANSFER 
!
!******************************************************************************
!  Subroutine CLEANUP_TRANSFER deallocates the EDGE_IN array.
!  (bmy, 9/19/01, 10/31/03)
!
!  NOTES:
!  (1 ) Renamed SIGE_IN to EDGE_IN to denote that it is not always a sigma
!        coordinate (as for GEOS-4). (bmy, 10/31/03)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_TRANSFER begins here!
      !=================================================================
      IF ( ALLOCATED( EDGE_IN ) ) DEALLOCATE( EDGE_IN )

      ! Return to calling program
      END SUBROUTINE CLEANUP_TRANSFER
      
!------------------------------------------------------------------------------

      ! End of module
      END MODULE TRANSFER_MOD
