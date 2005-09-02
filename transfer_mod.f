! $Id: transfer_mod.f,v 1.4 2005/09/02 15:17:28 bmy Exp $
      MODULE TRANSFER_MOD
!
!******************************************************************************
!  Module TRANSFER_MOD contains routines used to copy data from REAL*4 to
!  REAL*8 arrays after being read from disk.  Also, regridding of GEOS-3
!  vertical levels from 48 levels to 30 levels will be done if necessary.
!  (mje, bmy, 9/27/01, 6/7/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) EDGE_IN          : Input sigma edges (for pure sigma models)
!
!  Module Routines:
!  ============================================================================
!  (1 ) GET_L_COPY       : Returns # of vertical levels to be copied
!  (2 ) TRANSFER_A6      : Transfers GEOS-3 A-6   fields, regrids if necessary
!  (3 ) TRANSFER_3D      : Transfers GEOS-3 3-D   fields, regrids if necessary
!  (4 ) TRANSFER_3D_TROP : Transfers GEOS-3 3-D   fields up to tropopause level
!  (5 ) TRANSFER_ZONAL_R4: Transfers GEOS-3 zonal fields, regrids (REAL*4)
!  (6 ) TRANSFER_ZONAL_R8: Transfers GEOS-3 zonal fields, regrids (REAL*8) 
!  (7 ) TRANSFER_ZONAL   : Transfers GEOS-3 zonal fields, regrids if necessary 
!  (8 ) TRANSFER_2D_INT  : Transfers GEOS   2-D   fields (INTEGER argument)
!  (9 ) TRANSFER_2D_R4   : Transfers GEOS   2-D   fields (REAL*4 argument)
!  (10) TRANSFER_2D_R8   : Transfers GEOS   2-D   fields (REAL*8 argument)
!  (11) TRANSFER_TO_1D   : Transfers GEOS   2-D   fields to a 1-D array
!  (12) LUMP_2_R4        : Combines 2 sigma levels into 1 thick level (REAL*4) 
!  (13) LUMP_2_R8        : Combines 2 sigma levels into 1 thick level (REAL*8) 
!  (14) LUMP_4_R4        : Combines 4 sigma levels into 1 thick level (REAL*4)
!  (15) LUMP_4_R8        : Combines 4 sigma levels into 1 thick level (REAL*8)
!  (16) INIT_TRANSFER    : Allocates and initializes the EDGE_IN array
!  (17) CLEANUP_TRANSFER : Deallocates the EDGE_IN array
!
!  Module Interfaces:
!  ============================================================================
!  (1 ) LUMP_2           : Overloads routines LUMP_2_R4 and LUMP_2_R8
!  (2 ) LUMP_4           : Overloads routines LUMP_4_R4 and LUMP_4_R8
!  (3 ) TRANSFER_2D      : Overloads TRANSFER_2D_INT, TRANSFER_2D_R4
!                           and TRANSFER_2D_R8
!  (4 ) TRANSFER_ZONAL   : Overloads TRANSFER_ZONAL_R4 and TRANSFER_ZONAL_R8
!
!  GEOS-CHEM modules referenced by regrid_mod.f
!  ============================================================================
!  (1 ) error_mod.f      : Module containing NaN and other error check routines
!  (2 ) pressure_mod.f   : Module containing routines to compute P(I,J,L)
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
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "transfer_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE EDGE_IN

      ! PRIVATE module routines
      PRIVATE LUMP_2_R4,         LUMP_2_R8
      PRIVATE LUMP_4_R4,         LUMP_4_R8
      PRIVATE TRANSFER_2D_INT,   TRANSFER_2D_R4,    TRANSFER_2D_R8
      PRIVATE TRANSFER_ZONAL_R4, TRANSFER_ZONAL_R8, GET_L_COPY

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      REAL*8, ALLOCATABLE :: EDGE_IN(:)

      !=================================================================
      ! MODULE INTERFACES -- "bind" two or more routines with different
      ! argument types or # of arguments under one unique name
      !================================================================= 

      ! Interface for routines to lump 2 levels together (REAL*4 and REAL*8)
      INTERFACE LUMP_2
         MODULE PROCEDURE LUMP_2_R4
         MODULE PROCEDURE LUMP_2_R8
      END INTERFACE

      ! Interface for routines to lump 2 levels together (REAL*4 and REAL*8)
      INTERFACE LUMP_4
         MODULE PROCEDURE LUMP_4_R4
         MODULE PROCEDURE LUMP_4_R8
      END INTERFACE

      ! Interface for routines which copy 2-D data 
      ! (INTEGER, REAL*4, and REAL*8)
      INTERFACE TRANSFER_2D
         MODULE PROCEDURE TRANSFER_2D_INT
         MODULE PROCEDURE TRANSFER_2D_R4
         MODULE PROCEDURE TRANSFER_2D_R8
      END INTERFACE

      ! Interface for routines which copy lat-alt data (REAL*4 and REAL*8)
      INTERFACE TRANSFER_ZONAL
         MODULE PROCEDURE TRANSFER_ZONAL_R4
         MODULE PROCEDURE TRANSFER_ZONAL_R8
      END INTERFACE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION GET_L_COPY() RESULT( L_COPY )
!
!******************************************************************************
!  Function GET_L_COPY returns the value L_COPY, which is the number
!  of vertical levels to copy from a REAL*4 array to a REAL*8 array.  
!  Levels above L_COPY will be vertically regridded, if necessary. 
!  (bmy, 6/18/03, 5/24/05)
!
!               { LGLOB, for GCAP
!               { LGLOB, for GEOS-1
!               { LGLOB, for GEOS-STRAT
!      L_INT =  { LGLOB, for GEOS-3            (if LLPAR == LGLOB)
!               { 22   , for GEOS-3            (if LLPAR /= LGLOB)
!               { LGLOB, for GEOS-4 or GEOS-5  (if LLPAR == LGLOB)
!               { 19   , for GEOS-4 or GEOS-5  (if LLPAR /= LGLOB)
!
!  NOTES:
!  (1 ) Add value of L_COPY for GEOS-4/fvDAS (6/18/03)
!  (2 ) Add LCOPY = 19 for GEOS-4 if regridding (bmy, 10/31/03)
!  (3 ) Now modified for GEOS-5 and GCAP met fields (bmy, 5/24/05)
!******************************************************************************
!
#     include "CMN_SIZE"

      ! Function return variable
      INTEGER :: L_COPY

#if   defined( GEOS_1 ) 

      ! Copy all vertical levels for GEOS-1
      L_COPY = LGLOB

#elif defined( GEOS_STRAT )

      ! Copy all vertical levels for GEOS-STRAT
      L_COPY = LGLOB

#elif defined( GEOS_3 )

      ! For GEOS-3, only regrid if LLPAR does not equal LGLOB 
      IF ( LLPAR == LGLOB ) THEN
         L_COPY = LGLOB
      ELSE
         L_COPY = 22
      ENDIF
      
#elif defined( GEOS_4 ) || defined( GEOS_5 )

      ! For GEOS-4 & GEOS-5, only regrid if LLPAR does not equal LGLOB 
      IF ( LLPAR == LGLOB ) THEN
         L_COPY = LGLOB
      ELSE
         L_COPY = 19
      ENDIF

#elif defined( GCAP )

      ! Copy all vertical levels for GCAP met fields (swu, bmy, 5/24/05)
      L_COPY = LGLOB

#endif

      ! Return to calling program
      END FUNCTION GET_L_COPY

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_A6( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_A6 transfers A-6 data from a REAL*4 array of dimension
!  (IGLOB,JGLOB,LGLOB) to a REAL*8 array of dimension (LLPAR,IIPAR,JJPAR).
!  Regrid GEOS-3 data from 48 --> 30 levels or GEOS-4 data from 55 --> 30 
!  levels if necessary. (bmy, 9/21/01, 5/24/05)
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
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XOFFSET, GET_YOFFSET

#     include "CMN_SIZE"

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IGLOB,JGLOB,LGLOB)
      REAL*8,  INTENT(OUT) :: OUT(LLPAR,IIPAR,JJPAR)

      ! Local variables
      LOGICAL, SAVE        :: FIRST = .TRUE.  
      INTEGER              :: I, I0, IREF, J, J0, JREF, L, L_COPY
      REAL*4               :: INCOL(LGLOB)
      
      !================================================================
      ! TRANSFER_A6 begins here!
      !
      ! Also convert from global horizontal indices (IREF,JREF) to
      ! window indices (I,J), for consistency with existing code
      !================================================================

      ! Allocate EDGE_IN on the first call, if necessary
      IF ( FIRST ) THEN
         CALL INIT_TRANSFER
         FIRST = .FALSE.
      ENDIF

      ! Get nested-grid offsets 
      I0     = GET_XOFFSET()
      J0     = GET_YOFFSET()

      ! L_COPY is the # of levels to copy from REAL*4 to REAL*8
      L_COPY = GET_L_COPY()

      ! Copy levels up to L_COPY, for all model types
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, IREF, J, JREF, L )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
         JREF = J + J0
      DO I = 1, IIPAR
         IREF = I + I0
      DO L = 1, L_COPY
         OUT(L,I,J) = IN(IREF,JREF,L)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#if   defined( GEOS_3 )

      !================================================================
      ! For GEOS-3, when LLPAR /= LGLOB, then lump 48 levels into 
      ! 30 levels, according to the formulation proposed by Mat Evans.  
      ! Lump levels above L_COPY = 22 in groups of 2 or groups of 4.
      !================================================================ 

      ! Return if LLPAR = LGLOB
      IF ( LLPAR == LGLOB ) RETURN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, IREF, J, JREF, L, INCOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
         JREF = J + J0

         DO I = 1, IIPAR
            IREF = I + I0

            ! Store vertical column at (IREF,JREF) in a 1-D vector
            DO L = 1, LGLOB
               INCOL(L) = IN(IREF,JREF,L)
            ENDDO

            ! Lump 2 levels together at a time, starting at the given level
            OUT(23,I,J) = LUMP_2( INCOL, LGLOB, 23 )
            OUT(24,I,J) = LUMP_2( INCOL, LGLOB, 25 )
            OUT(25,I,J) = LUMP_2( INCOL, LGLOB, 27 )

            ! Lump 4 levels together at a time, starting at the given level
            OUT(26,I,J) = LUMP_4( INCOL, LGLOB, 29 )
            OUT(27,I,J) = LUMP_4( INCOL, LGLOB, 33 )
            OUT(28,I,J) = LUMP_4( INCOL, LGLOB, 37 ) 
            OUT(29,I,J) = LUMP_4( INCOL, LGLOB, 41 ) 
            OUT(30,I,J) = LUMP_4( INCOL, LGLOB, 45 ) 
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

#elif defined( GEOS_4 ) || defined( GEOS_5 )

      !================================================================
      ! GEOS-4/GEOS-5: When LLPAR /= LGLOB, then lump 55 levels into 
      ! 30 levels, according to the formulation proposed by Mat Evans.  
      ! Lump levels above L_COPY = 19 in groups of 2 or groups of 4.
      !================================================================ 

      ! Return if LLPAR = LGLOB
      IF ( LLPAR == LGLOB ) RETURN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, IREF, J, JREF, L, INCOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
         JREF = J + J0

         DO I = 1, IIPAR
            IREF = I + I0

            ! Store vertical column at (IREF,JREF) in a 1-D vector
            DO L = 1, LGLOB
               INCOL(L) = IN(IREF,JREF,L)
            ENDDO

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
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

#endif

      ! Return to calling program
      END SUBROUTINE TRANSFER_A6

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_3D( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_3D transfers 3-D data from a REAL*4 array of dimension
!  (IGLOB,JGLOB,LGLOB) to a REAL*8 array of dimension (IIPAR,JJPAR,LLPAR).
!  Regrid GEOS-3 data from 48 --> 30 levels or GEOS-4 data from 55 --> 30 
!  levels if necessary. (bmy, 9/21/01, 5/24/05)
!
!  NOTE: TRANSFER_3D can be used w/ I-6 met fields, since these are of
!  dimensions (IIPAR,JJPAR,LLPAR).
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
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XOFFSET, GET_YOFFSET

#     include "CMN_SIZE"

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IIPAR,JJPAR,LGLOB)
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLPAR)

      ! Local variables
      LOGICAL, SAVE        :: FIRST = .TRUE.  
      INTEGER              :: I, I0, IREF, J, J0, JREF, L, L_COPY
      REAL*4               :: INCOL(LGLOB)
     
      !================================================================
      ! TRANSFER_3D begins here!
      !================================================================

      ! Allocate EDGE_IN on the first call, if necessary
      IF ( FIRST ) THEN
         CALL INIT_TRANSFER
         FIRST = .FALSE.
      ENDIF

      ! Get nested-grid offsets
      I0     = GET_XOFFSET()
      J0     = GET_YOFFSET()

      ! L_COPY is the # of levels to copy from REAL*4 to REAL*8
      L_COPY = GET_L_COPY()

      ! Copy levels up to L_COPY, for all model types
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, IREF, J, JREF, L )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, L_COPY
      DO J = 1, JJPAR
         JREF = J + J0
      DO I = 1, IIPAR
         IREF = I + I0
         OUT(I,J,L) = IN(IREF,JREF,L)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#if   defined( GEOS_3 )

      !================================================================
      ! For GEOS-3, when LLPAR /= LGLOB, then lump 48 levels into 
      ! 30 levels, according to the formulation proposed by Mat Evans.  
      ! Lump levels above L_INT = 22 in groups of 2 or groups of 4.
      !================================================================    

      ! Return if LLPAR = LGLOB
      IF ( LLPAR == LGLOB ) RETURN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, IREF, J, JREF, L, INCOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
         JREF = J + J0

         DO I = 1, IIPAR
            IREF = I + I0

            ! Store vertical column at (I,J) in a 1-D vector
            DO L = 1, LGLOB
               INCOL(L) = IN(IREF,JREF,L)
            ENDDO

            ! Lump 2 levels together at a time, starting at the given level
            OUT(I,J,23) = LUMP_2( INCOL, LGLOB, 23 )
            OUT(I,J,24) = LUMP_2( INCOL, LGLOB, 25 )
            OUT(I,J,25) = LUMP_2( INCOL, LGLOB, 27 )

            ! Lump 4 levels together at a time, starting at the given level
            OUT(I,J,26) = LUMP_4( INCOL, LGLOB, 29 )
            OUT(I,J,27) = LUMP_4( INCOL, LGLOB, 33 )
            OUT(I,J,28) = LUMP_4( INCOL, LGLOB, 37 ) 
            OUT(I,J,29) = LUMP_4( INCOL, LGLOB, 41 ) 
            OUT(I,J,30) = LUMP_4( INCOL, LGLOB, 45 ) 
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

#elif defined( GEOS_4 ) || defined( GEOS_5 )

      !================================================================
      ! GEOS-4/GEOS-5: When LLPAR /= LGLOB, then lump 55 levels into 
      ! 30 levels, according to the formulation proposed by Mat Evans.  
      ! Lump levels above L_COPY = 19 in groups of 2 or groups of 4.
      !================================================================ 

      ! Return if LLPAR = LGLOB
      IF ( LLPAR == LGLOB ) RETURN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, IREF, J, JREF, L, INCOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
         JREF = J + J0

         DO I = 1, IIPAR
            IREF = I + I0

            ! Store vertical column at (I,J) in a 1-D vector
            DO L = 1, LGLOB
               INCOL(L) = IN(IREF,JREF,L)
            ENDDO

            ! Lump 2 levels together at a time
            OUT(I,J,20) = LUMP_2( INCOL, LGLOB, 20 )
            OUT(I,J,21) = LUMP_2( INCOL, LGLOB, 22 )
            OUT(I,J,22) = LUMP_2( INCOL, LGLOB, 24 )
            OUT(I,J,23) = LUMP_2( INCOL, LGLOB, 26 )

            ! Lump 4 levels together at a time
            OUT(I,J,24) = LUMP_4( INCOL, LGLOB, 28 )
            OUT(I,J,25) = LUMP_4( INCOL, LGLOB, 32 )
            OUT(I,J,26) = LUMP_4( INCOL, LGLOB, 36 ) 
            OUT(I,J,27) = LUMP_4( INCOL, LGLOB, 40 ) 
            OUT(I,J,28) = LUMP_4( INCOL, LGLOB, 44 ) 
            OUT(I,J,29) = LUMP_4( INCOL, LGLOB, 48 ) 
            OUT(I,J,30) = LUMP_4( INCOL, LGLOB, 52 ) 
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

#endif

      ! Return to calling program
      END SUBROUTINE TRANSFER_3D

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_3D_TROP( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_3D_TROP transfers 3-D data from a REAL*4 array of 
!  dimension (IGLOB,JGLOB,LLTROP) to a REAL*8 array of dimension 
!  (IIPAR,JJPAR,LLTROP).  Use this routine to cast & resize data that are
!  only defined up to the max tropopause level LLTROP. (bmy, 10/31/02, 3/11/03)
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
!******************************************************************************
!
#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IGLOB,JGLOB,LLTROP)
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR,LLTROP)

      ! Local variables
      INTEGER              :: L
      
      !=================================================================
      ! TRANSFER_3D_TROP
      !=================================================================

      ! Cast to REAL*8 abd resize up to LLTROP
      DO L = 1, LLTROP
         CALL TRANSFER_2D( IN(:,:,L), OUT(:,:,L) )
      ENDDO
      
      ! Return to calling program
      END SUBROUTINE TRANSFER_3D_TROP 

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_ZONAL_R4( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_ZONAL_R4 transfers zonal mean or lat-alt data from a 
!  REAL*4 array of dimension (JGLOB,LGLOB) to a REAL*4 array of dimension 
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
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_YOFFSET

#     include "CMN_SIZE"

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(JGLOB,LGLOB)
      REAL*4,  INTENT(OUT) :: OUT(JJPAR,LLPAR)

      ! Local variables
      LOGICAL, SAVE        :: FIRST = .TRUE.  
      INTEGER              :: J, J0, JREF, L, L_COPY
      REAL*4               :: INCOL(LGLOB)
     
      !================================================================
      ! TRANSFER_ZONAL_R4 begins here!
      !================================================================

      ! Allocate EDGE_IN on the first call, if necessary
      IF ( FIRST ) THEN
         CALL INIT_TRANSFER
         FIRST = .FALSE.
      ENDIF

      ! Get nested grid offsets
      J0     = GET_YOFFSET()

      ! L_COPY is the # of levels to copy from REAL*4 to REAL*8
      L_COPY = GET_L_COPY()

      ! Copy levels up to L_COPY, for all model types
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J, JREF, L )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, L_COPY
      DO J = 1, JJPAR
         JREF = J + J0
         OUT(J,L) = IN(JREF,L)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#if   defined( GEOS_3 )

      !================================================================
      ! For GEOS-3, when LLPAR /= LGLOB, then lump 48 levels into 
      ! 30 levels, according to the formulation proposed by Mat Evans.  
      ! Lump levels above L_INT = 22 in groups of 2 or groups of 4.
      !================================================================    

      ! Return if LLPAR = LGLOB
      IF ( LLPAR == LGLOB ) RETURN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J, JREF, L, INCOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
         JREF = J + J0

         ! Store vertical column at (I,J) in a 1-D vector
         DO L = 1, LGLOB
            INCOL(L) = IN(JREF,L)
         ENDDO

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
      ENDDO
!$OMP END PARALLEL DO

#elif defined( GEOS_4 ) || defined( GEOS_5 )
      
      !================================================================
      ! GEOS-4/GEOS-5: When LLPAR /= LGLOB, then lump 55 levels into 
      ! 30 levels, according to the formulation proposed by Mat Evans.  
      ! Lump levels above L_COPY = 19 in groups of 2 or groups of 4.
      !================================================================ 

      ! Return if LLPAR = LGLOB
      IF ( LLPAR == LGLOB ) RETURN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J, JREF, L, INCOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
         JREF = J + J0

         ! Store vertical column at (I,J) in a 1-D vector
         DO L = 1, LGLOB
            INCOL(L) = IN(JREF,L)
         ENDDO

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
      ENDDO
!$OMP END PARALLEL DO

#endif

      ! Return to calling program
      END SUBROUTINE TRANSFER_ZONAL_R4

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_ZONAL_R8( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_ZONAL_R8 transfers zonal mean or lat-alt data from a 
!  REAL*4 array of dimension (JGLOB,LGLOB) to a REAL*8 array of dimension 
!  (JJPAR,LLPAR).  Regrid GEOS-3 data from 48 --> 30 levels or GEOS-4 data
!  from 55 --> 30 levels if necessary. (bmy, 9/21/01, 5/24/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IN  (REAL*4) : Input  field, of dimension (JGLOB,LGLOB)
!  (2 ) OUT (REAL*8) : Output field, of dimension (JJPAR,LLPAR)
!
!  NOTES:
!  (1 ) Lump levels together in groups of 2 or 4, as dictated by Mat Evans.
!        (bmy, 9/21/01)
!  (2 ) Assumes that LLPAR == LGLOB for GEOS-1, GEOS-STRAT (bmy, 9/21/01)
!  (3 ) Now use functions GET_YOFFSET from "grid_mod.f".  Now J0 is a local 
!        variable. (bmy, 3/11/03)
!  (4 ) Added code to regrid GEOS-4 from 55 --> 30 levels (mje, bmy, 10/31/03)
!  (5 ) Now modified for GEOS-5 met fields (bmy, 5/24/05)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_YOFFSET

#     include "CMN_SIZE"

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(JGLOB,LGLOB)
      REAL*8,  INTENT(OUT) :: OUT(JJPAR,LLPAR)

      ! Local variables
      LOGICAL, SAVE        :: FIRST = .TRUE.  
      INTEGER              :: J, J0, JREF, L, L_COPY
      REAL*4               :: INCOL(LGLOB)
     
      !================================================================
      ! TRANSFER_ZONAL begins here!
      !================================================================

      ! Allocate EDGE_IN on the first call, if necessary
      IF ( FIRST ) THEN
         CALL INIT_TRANSFER
         FIRST = .FALSE.
      ENDIF

      ! Get nested-grid offset
      J0     = GET_YOFFSET()

      ! L_COPY is the # of levels to copy from REAL*4 to REAL*8
      L_COPY = GET_L_COPY()

      ! Copy levels up to L_COPY, for all model types
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J, JREF, L )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, L_COPY
      DO J = 1, JJPAR
         JREF = J + J0
         OUT(J,L) = IN(JREF,L)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#if   defined( GEOS_3 )

      !================================================================
      ! For GEOS-3, when LLPAR /= LGLOB, then lump 48 levels into 
      ! 30 levels, according to the formulation proposed by Mat Evans.  
      ! Lump levels above L_COPY = 22 in groups of 2 or groups of 4.
      !================================================================    

      ! Return if LLPAR = LGLOB
      IF ( LLPAR == LGLOB ) RETURN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J, JREF, L, INCOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
         JREF = J + J0

         ! Store vertical column at (I,J) in a 1-D vector
         DO L = 1, LGLOB
            INCOL(L) = IN(JREF,L)
         ENDDO

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
      ENDDO
!$OMP END PARALLEL DO

#elif defined( GEOS_4 ) || defined( GEOS_5 )
      
      !================================================================
      ! GEOS-4/GEOS-5: When LLPAR /= LGLOB, then lump 55 levels into 
      ! 30 levels, according to the formulation proposed by Mat Evans.  
      ! Lump levels above L_COPY = 19 in groups of 2 or groups of 4.
      !================================================================  

      ! Return if LLPAR = LGLOB
      IF ( LLPAR == LGLOB ) RETURN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( J, JREF, L, INCOL )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
         JREF = J + J0

         ! Store vertical column at (I,J) in a 1-D vector
         DO L = 1, LGLOB
            INCOL(L) = IN(JREF,L)
         ENDDO

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
      ENDDO
!$OMP END PARALLEL DO

#endif

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
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XOFFSET, GET_YOFFSET

#     include "CMN_SIZE"

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IGLOB,JGLOB)
      INTEGER, INTENT(OUT) :: OUT(IIPAR,JJPAR)

      ! Local variables
      INTEGER              :: I, I0, IREF, J, J0, JREF

      !=================================================================
      ! TRANSFER_2D_INT begins here!
      !=================================================================

      ! Get nested-grid offsets
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, IREF, J, JREF )
      DO J = 1, JJPAR
         JREF = J + J0
      DO I = 1, IIPAR
         IREF = I + I0
         OUT(I,J) = IN(IREF,JREF)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
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
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XOFFSET, GET_YOFFSET

#     include "CMN_SIZE"

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IGLOB,JGLOB)
      REAL*4,  INTENT(OUT) :: OUT(IIPAR,JJPAR)

      ! Local variables
      INTEGER              :: I, I0, IREF, J, J0, JREF

      !=================================================================
      ! TRANSFER_2D_R4 begins here!
      !=================================================================
      
      ! Get nested-grid offsets
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, IREF, J, JREF )
      DO J = 1, JJPAR
         JREF = J + J0
      DO I = 1, IIPAR
         IREF = I + I0
         OUT(I,J) = IN(IREF,JREF)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
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
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XOFFSET, GET_YOFFSET

#     include "CMN_SIZE"

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IGLOB,JGLOB)
      REAL*8,  INTENT(OUT) :: OUT(IIPAR,JJPAR)

      ! Local variables
      INTEGER              :: I, I0, IREF, J, J0, JREF

      !=================================================================
      ! TRANSFER_2D_R8 begins here!
      !=================================================================

      ! Get nested-grid offsets
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, IREF, J, JREF )
      DO J = 1, JJPAR
         JREF = J + J0
      DO I = 1, IIPAR
         IREF = I + I0
         OUT(I,J) = IN(IREF,JREF)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
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
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XOFFSET, GET_YOFFSET

#     include "CMN_SIZE"

      ! Arguments
      REAL*4,  INTENT(IN)  :: IN(IGLOB,JGLOB)
      REAL*8,  INTENT(OUT) :: OUT(MAXIJ)

      ! Local variables
      INTEGER              :: I, I0, IREF, J, J0, JREF, IJLOOP

      !=================================================================
      ! TRANSFER_TO_1D begins here!
      !=================================================================

      ! Get nested-grid offsets
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

      ! 1-D counter
      IJLOOP = 0

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

      SUBROUTINE INIT_TRANSFER
!
!******************************************************************************
!  Subroutine INIT_TRANSFER initializes and zeroes the EDGE_IN array.
!  (bmy, 9/19/01, 5/24/05)
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
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE PRESSURE_MOD, ONLY : GET_BP

#     include "CMN_SIZE" 

      ! Local variables
      INTEGER :: AS, L
      
      !=================================================================
      ! INIT_TRANSFER begins here!
      !=================================================================

      ! Test if EDGE_IN has already been allocated
      IF ( .not. ALLOCATED( EDGE_IN ) ) THEN

         ! Allocate the EDGE_IN array
         ALLOCATE( EDGE_IN( LGLOB + 1 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EDGE_IN' )

#if   defined( GEOS_4 ) || defined( GEOS_5 )

         ! For GEOS-4, this is a hybrid grid.  Assign the original 56 
         ! Ap values (which for stratospheric levels are pressures
         ! at grid box edges in [hPa]) to the EDGE_IN array.
         EDGE_IN = (/  0.000000d0,   0.000000d0,  12.704939d0,  
     &                35.465965d0,  66.098427d0, 101.671654d0, 
     &               138.744400d0, 173.403183d0, 198.737839d0, 
     &               215.417526d0, 223.884689d0, 224.362869d0,
     &               216.864929d0, 201.192093d0, 176.929993d0, 
     &               150.393005d0, 127.837006d0, 108.663429d0,  
     &                92.365662d0,  78.512299d0,  66.603378d0,  
     &                56.387939d0,  47.643932d0,  40.175419d0, 
     &                33.809956d0,  28.367815d0,  23.730362d0,  
     &                19.791553d0,  16.457071d0,  13.643393d0,  
     &                11.276889d0,   9.292943d0,   7.619839d0,   
     &                 6.216800d0,   5.046805d0,   4.076567d0, 
     &                 3.276433d0,   2.620212d0,   2.084972d0,   
     &                 1.650792d0,   1.300508d0,   1.019442d0,   
     &                 0.795134d0,   0.616779d0,   0.475806d0,   
     &                 0.365041d0,   0.278526d0,   0.211349d0, 
     &                 0.159495d0,   0.119703d0,   0.089345d0,   
     &                 0.066000d0,   0.047585d0,   0.032700d0,   
     &                 0.020000d0,   0.010000d0 /)

#elif defined( GEOS_3 ) 

         ! For GEOS-3, this is a pure-sigma grid.  Assign 
         ! the original 49 sigma edges to the EDGE_IN array.
         EDGE_IN = (/ 1.000000d0, 0.997095d0, 0.991200d0, 0.981500d0,    
     &                0.967100d0, 0.946800d0, 0.919500d0, 0.884000d0,    
     &                0.839000d0, 0.783000d0, 0.718200d0, 0.647600d0,    
     &                0.574100d0, 0.500000d0, 0.427800d0, 0.359500d0,    
     &                0.297050d0, 0.241950d0, 0.194640d0, 0.155000d0,    
     &                0.122680d0, 0.096900d0, 0.076480d0, 0.060350d0,   
     &                0.047610d0, 0.037540d0, 0.029600d0, 0.023330d0,   
     &                0.018380d0, 0.014480d0, 0.011405d0, 0.008975d0,  
     &                0.007040d0, 0.005500d0, 0.004280d0, 0.003300d0,  
     &                0.002530d0, 0.001900d0, 0.001440d0, 0.001060d0,  
     &                0.000765d0, 0.000540d0, 0.000370d0, 0.000245d0, 
     &                0.000155d0, 9.20000d-5, 4.75000d-5, 1.76800d-5, 
     &                0.000000d0 /)

#else

         ! We are not reducing the layers for GEOS-1 and GEOS-STRAT,
         ! thus LGLOB = LLPAR.  Use GET_BP to initialize EDGE_IN,
         ! since for a pure sigma grid, BP are just the sigma edges.
         DO L = 1, LGLOB+1
            EDGE_IN(L) = GET_BP(L)
         ENDDO

#endif

      ENDIF

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
      IF ( ALLOCATED( EDGE_IN ) ) DEALLOCATE( EDGE_IN )

      ! Return to calling program
      END SUBROUTINE CLEANUP_TRANSFER
      
!------------------------------------------------------------------------------

      END MODULE TRANSFER_MOD
