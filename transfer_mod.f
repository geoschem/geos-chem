! $Id: transfer_mod.f,v 1.1 2003/06/30 20:26:05 bmy Exp $
      MODULE TRANSFER_MOD
!
!******************************************************************************
!  Module TRANSFER_MOD contains routines used to copy data from REAL*4 to
!  REAL*8 arrays after being read from disk.  Also, regridding of GEOS-3
!  vertical levels from 48 levels to 30 levels will be done if necessary.
!  (mje, bmy, 9/27/01, 3/11/03)
!
!  Module Variables:
!  ============================================================================
!  (1 ) SIGE_IN          : Input sigma edges for the current model grid
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
!  (16) INIT_TRANSFER    : Allocates and initializes the SIGE_IN array
!  (17) CLEANUP_TRANSFER : Deallocates the SIGE_IN array
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
!  (3 ) SIGE_IN needs to be provided for each model type, within an #ifdef
!        block, in order to ensure compilation.  However, SIGE_IN is currently
!        only used for regridding GEOS-3 data (and probably also GEOS-4 when 
!        that becomes available). (bmy, 9/26/01)
!  (4 ) Add interfaces TRANSFER_2D and TRANSFER_ZONAL (bmy, 9/27/01)
!  (5 ) Added routine TRANSFER_2D_R4.  Added TRANSFER_2D_R4 to the generic
!        TRANSFER_2D interface. (bmy, 1/25/02)
!  (6 ) Updated comments, cosmetic changes (bmy, 2/28/02) 
!  (7 ) Bug fix: remove extraneous "," in GEOS-1 definition of SIGE_IN array.
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
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "transfer_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE SIGE_IN

      ! PRIVATE module routines
      PRIVATE LUMP_2_R4,         LUMP_2_R8
      PRIVATE LUMP_4_R4,         LUMP_4_R8
      PRIVATE TRANSFER_2D_INT,   TRANSFER_2D_R4,    TRANSFER_2D_R8
      PRIVATE TRANSFER_ZONAL_R4, TRANSFER_ZONAL_R8, GET_L_COPY

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      REAL*8, ALLOCATABLE :: SIGE_IN(:)

      !=================================================================
      ! MODULE INTERFACES -- "bind" two or more routines with different
      ! argument types or # of arguments under one unique name
      !================================================================= 

      ! Interface for routines to lump 2 levels together (REAL*4 and REAL*8)
      INTERFACE LUMP_2
         MODULE PROCEDURE LUMP_2_R4, LUMP_2_R8
      END INTERFACE

      ! Interface for routines to lump 2 levels together (REAL*4 and REAL*8)
      INTERFACE LUMP_4
         MODULE PROCEDURE LUMP_4_R4, LUMP_4_R8
      END INTERFACE

      ! Interface for routines which copy 2-D data 
      ! (INTEGER, REAL*4, and REAL*8)
      INTERFACE TRANSFER_2D
         MODULE PROCEDURE TRANSFER_2D_INT,
     &                    TRANSFER_2D_R4, 
     &                    TRANSFER_2D_R8
      END INTERFACE

      ! Interface for routines which copy lat-alt data (REAL*4 and REAL*8)
      INTERFACE TRANSFER_ZONAL
         MODULE PROCEDURE TRANSFER_ZONAL_R4, TRANSFER_ZONAL_R8
      END INTERFACE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION GET_L_COPY() RESULT( L_COPY )
!
!******************************************************************************
!  Function GET_L_COPY returns the value L_COPY, which is the number of 
!  vertical levels to copy from a REAL*4 array to a REAL*8 array.  Levels
!  above L_COPY will be vertically regridded, if necessary. (bmy, 6/18/03)
!
!               { LGLOB,            for GEOS-1
!               { LGLOB,            for GEOS-STRAT
!      L_INT =  { LGLOB,            for GEOS-3 if LLPAR == LGLOB 
!               { 22   ,            for GEOS-3 if LLPAR /= LGLOB
!               { to be deterimned, for GEOS-4/fvDAS
!
!  NOTES:
!  (1 ) Add value of L_COPY for GEOS-4/fvDAS (bmy, 6/18/03)
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
      
#elif defined( GEOS_4 )

      ! For GEOS-4, only regrid if LLPAR does not equal LGLOB 
      IF ( LLPAR == LGLOB ) THEN
         L_COPY = LGLOB
      ELSE
         L_COPY = 22  ! change later (bmy, 6/18/03)
      ENDIF
#endif

      ! Return to calling program
      END FUNCTION GET_L_COPY

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_A6( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_A6 transfers A-6 data from a REAL*4 array of dimension
!  (IGLOB,JGLOB,LGLOB) to a REAL*8 array of dimension (LLPAR,IIPAR,JJPAR).
!  Also, GEOS-3 data will be regridded from 48 --> 30 levels if necessary.
!  (bmy, 9/21/01, 3/11/03)
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

      ! Allocate SIGE_IN on the first call, if necessary
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

#elif defined( GEOS_4 )

      !================================================================
      ! Leave room for future expansion to GEOS-4 
      !================================================================    

#endif

      ! Return to calling program
      END SUBROUTINE TRANSFER_A6

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_3D( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_3D transfers 3-D data from a REAL*4 array of dimension
!  (IGLOB,JGLOB,LGLOB) to a REAL*8 array of dimension (IIPAR,JJPAR,LLPAR).
!  Also, GEOS-3 data will be regridded from 48 --> 30 levels if necessary.
!  (bmy, 9/21/01, 3/11/03)
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
      ! TRANSFER_I6 begins here!
      !
      ! Also convert from global horizontal indices (IREF,JREF) to
      ! window indices (I,J), for consistency with existing code
      !================================================================

      ! Allocate SIGE_IN on the first call, if necessary
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

#elif defined( GEOS_4 )

      !================================================================
      ! Leave room for future expansion to GEOS-4 
      !================================================================    

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
!  (JJPAR,LLPAR).  Also, GEOS-3 data will be regridded from  48 --> 30 levels 
!  if necessary. (bmy, 9/21/01, 3/11/03)
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
      !
      ! Also convert from global horizontal indices (IREF,JREF) to
      ! window indices (I,J), for consistency with existing code
      !================================================================

      ! Allocate SIGE_IN on the first call, if necessary
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

         ! Lump 2 levels together at a time, starting at the given level
         OUT(J,23) = LUMP_2( INCOL, LGLOB, 23 )
         OUT(J,24) = LUMP_2( INCOL, LGLOB, 25 )
         OUT(J,25) = LUMP_2( INCOL, LGLOB, 27 )

         ! Lump 4 levels together at a time, starting at the given level
         OUT(J,26) = LUMP_4( INCOL, LGLOB, 29 )
         OUT(J,27) = LUMP_4( INCOL, LGLOB, 33 )
         OUT(J,28) = LUMP_4( INCOL, LGLOB, 37 ) 
         OUT(J,29) = LUMP_4( INCOL, LGLOB, 41 ) 
         OUT(J,30) = LUMP_4( INCOL, LGLOB, 45 ) 
      ENDDO
!$OMP END PARALLEL DO

#elif defined( GEOS_4 )
      
      !================================================================
      ! Leave room for future expansion to GEOS-4 
      !================================================================  

#endif

      ! Return to calling program
      END SUBROUTINE TRANSFER_ZONAL_R4

!------------------------------------------------------------------------------

      SUBROUTINE TRANSFER_ZONAL_R8( IN, OUT )
!
!******************************************************************************
!  Subroutine TRANSFER_ZONAL_R8 transfers zonal mean or lat-alt data from a 
!  REAL*4 array of dimension (JGLOB,LGLOB) to a REAL*8 array of dimension 
!  (JJPAR,LLPAR).  Also, GEOS-3 data will be regridded from  48 --> 30 levels 
!  if necessary. (bmy, 9/21/01)
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
      !
      ! Also convert from global horizontal indices (IREF,JREF) to
      ! window indices (I,J), for consistency with existing code
      !================================================================

      ! Allocate SIGE_IN on the first call, if necessary
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

         ! Lump 2 levels together at a time, starting at the given level
         OUT(J,23) = LUMP_2( INCOL, LGLOB, 23 )
         OUT(J,24) = LUMP_2( INCOL, LGLOB, 25 )
         OUT(J,25) = LUMP_2( INCOL, LGLOB, 27 )

         ! Lump 4 levels together at a time, starting at the given level
         OUT(J,26) = LUMP_4( INCOL, LGLOB, 29 )
         OUT(J,27) = LUMP_4( INCOL, LGLOB, 33 )
         OUT(J,28) = LUMP_4( INCOL, LGLOB, 37 ) 
         OUT(J,29) = LUMP_4( INCOL, LGLOB, 41 ) 
         OUT(J,30) = LUMP_4( INCOL, LGLOB, 45 ) 
      ENDDO
!$OMP END PARALLEL DO

#elif defined( GEOS_4 )
      
      !================================================================
      ! Leave room for future expansion to GEOS-4 
      !================================================================  

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
!  Input arguments must be REAL*4. (bmy, 9/18/01)
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
      ! The quantity Q(L') on the new merged sigma level is given 
      ! by the weighted average:
      !
      !             [ ( Q(L  ) * ( SIGE(L  ) - SIGE(L+1) ) ) + 
      !               ( Q(L+1) * ( SIGE(L+1) - SIGE(L+2) ) ) ]
      !  Q(L') = ------------------------------------------------
      !                       SIGE(L) - SIGE(L+2)
      !=================================================================     

      ! Compute weighted sum over 2 sigma levels
      OUT   = ( IN(L  ) * ( SIGE_IN(L  ) - SIGE_IN(L+1) ) ) +
     &        ( IN(L+1) * ( SIGE_IN(L+1) - SIGE_IN(L+2) ) )

      ! Divde by sigma thickness of new thick level
      OUT   = OUT / ( SIGE_IN(L) - SIGE_IN(L+2) )
       
      ! Return to calling routine
      END FUNCTION LUMP_2_R4

!------------------------------------------------------------------------------

      FUNCTION LUMP_2_R8( IN, L_IN, L ) RESULT( OUT )
!
!******************************************************************************
!  Function LUMP_2_R8 lumps 2 sigma levels into one thick level. 
!  Input arguments must be REAL*8. (bmy, 9/18/01, 10/15/02)
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
      ! The quantity Q(L') on the new merged sigma level is given by
      ! the weighted average:
      !
      !             [ ( Q(L  ) * ( SIGE(L  ) - SIGE(L+1) ) ) + 
      !               ( Q(L+1) * ( SIGE(L+1) - SIGE(L+2) ) ) ]
      !  Q(L') = ------------------------------------------------
      !                       SIGE(L) - SIGE(L+2)
      !=================================================================     

      ! Compute weighted sum over 2 sigma levels
      OUT   = ( IN(L  ) * ( SIGE_IN(L  ) - SIGE_IN(L+1) ) ) +
     &        ( IN(L+1) * ( SIGE_IN(L+1) - SIGE_IN(L+2) ) )

      ! Divde by sigma thickness of new thick level
      OUT   = OUT / ( SIGE_IN(L) - SIGE_IN(L+2) )
       
      ! Return to calling routine
      END FUNCTION LUMP_2_R8

!------------------------------------------------------------------------------

      FUNCTION LUMP_4_R4( IN, L_IN, L ) RESULT( OUT )
!
!******************************************************************************
!  Function LUMP_4_R4 lumps 4 sigma levels into one thick level. 
!  Input arguments must be REAL*4. (bmy, 9/18/01, 10/15/02)
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
      ! The quantity Q(L') on the new merged sigma level is given by
      ! the weighted average:
      !
      !             [ ( Q(L  ) * ( SIGE(L  ) - SIGE(L+1) ) ) + 
      !               ( Q(L+1) * ( SIGE(L+1) - SIGE(L+2) ) ) +
      !               ( Q(L+2) * ( SIGE(L+2) - SIGE(L+3) ) ) +
      !               ( Q(L+3) * ( SIGE(L+3) - SIGE(L+4) ) ) ]
      !  Q(L') = ------------------------------------------------
      !                       SIGE(L) - SIGE(L+4)
      !=================================================================     

      ! Compute weighted sum over 4 sigma levels
      OUT   = ( IN(L  ) * ( SIGE_IN(L  ) - SIGE_IN(L+1) ) ) +
     &        ( IN(L+1) * ( SIGE_IN(L+1) - SIGE_IN(L+2) ) ) +
     &        ( IN(L+2) * ( SIGE_IN(L+2) - SIGE_IN(L+3) ) ) +
     &        ( IN(L+3) * ( SIGE_IN(L+3) - SIGE_IN(L+4) ) ) 

      ! Divde by sigma thickness of new thick level
      OUT   = OUT / ( SIGE_IN(L) - SIGE_IN(L+4) )
       
      ! Return to calling routine
      END FUNCTION LUMP_4_R4

!------------------------------------------------------------------------------

      FUNCTION LUMP_4_R8( IN, L_IN, L ) RESULT( OUT )
!
!******************************************************************************
!  Function LUMP_4_R8 lumps 4 sigma levels into one thick level. 
!  Input arguments must be REAL*8. (bmy, 9/18/01, 10/15/02)
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
      ! The quantity Q(L') on the new merged sigma level is given by
      ! the weighted average:
      !
      !             [ ( Q(L  ) * ( SIGE(L  ) - SIGE(L+1) ) ) + 
      !               ( Q(L+1) * ( SIGE(L+1) - SIGE(L+2) ) ) +
      !               ( Q(L+2) * ( SIGE(L+2) - SIGE(L+3) ) ) +
      !               ( Q(L+3) * ( SIGE(L+3) - SIGE(L+4) ) ) ]
      !  Q(L') = ------------------------------------------------
      !                       SIGE(L) - SIGE(L+4)
      !=================================================================     

      ! Compute weighted sum over 4 sigma levels
      OUT   = ( IN(L  ) * ( SIGE_IN(L  ) - SIGE_IN(L+1) ) ) +
     &        ( IN(L+1) * ( SIGE_IN(L+1) - SIGE_IN(L+2) ) ) +
     &        ( IN(L+2) * ( SIGE_IN(L+2) - SIGE_IN(L+3) ) ) +
     &        ( IN(L+3) * ( SIGE_IN(L+3) - SIGE_IN(L+4) ) ) 

      ! Divde by sigma thickness of new thick level
      OUT   = OUT / ( SIGE_IN(L) - SIGE_IN(L+4) )
       
      ! Return to calling routine
      END FUNCTION LUMP_4_R8

!------------------------------------------------------------------------------

      SUBROUTINE INIT_TRANSFER
!
!******************************************************************************
!  Subroutine INIT_TRANSFER initializes and zeroes the SIGE_IN array.
!  (bmy, 9/19/01, 10/15/02)
!
!  NOTES:
!  (1 ) Removed additional "," for GEOS-1 definition of SIGE_IN (bmy, 3/25/02)
!  (2 ) Now use GET_BP from "pressure_mod.f" to get sigma edges for all
!        grids except GEOS-3 (dsa, bdf, bmy, 8/22/02)
!  (3 ) Declare L as a local variable.  Also reference ALLOC_ERR from module
!        "error_mod.f" (bmy, 10/15/02)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE PRESSURE_MOD, ONLY : GET_BP

#     include "CMN_SIZE" 

      ! Local variables
      INTEGER       :: AS, L
      
      !=================================================================
      ! INIT_TRANSFER begins here!
      !
      ! Allocate, initialize, and assign values to SIGE_IN
      !=================================================================
      IF ( .not. ALLOCATED( SIGE_IN ) ) THEN
         ALLOCATE( SIGE_IN( LGLOB + 1 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SIGE_IN' )

#if   defined( GEOS_3 ) 

         ! SIGE_IN = 49 sigma edges 
         SIGE_IN = (/ 1.000000d0, 0.997095d0, 0.991200d0, 0.981500d0,    
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
         ! thus LGLOB = LLPAR.  Use GET_BP to initialize SIGE_IN,
         ! since for a pure sigma grid, BP are just the sigma edges.
         DO L = 1, LGLOB+1
            SIGE_IN(L) = GET_BP(L)
         ENDDO

#endif

      ENDIF

      ! Return to calling program
      END SUBROUTINE INIT_TRANSFER
      
!------------------------------------------------------------------------------
      
      SUBROUTINE CLEANUP_TRANSFER 
!
!******************************************************************************
!  Subroutine CLEANUP_TRANSFER deallocates the SIGE_IN array (bmy, 9/19/01)
!******************************************************************************
!
      IF ( ALLOCATED( SIGE_IN ) ) DEALLOCATE( SIGE_IN )

      END SUBROUTINE CLEANUP_TRANSFER
      
!------------------------------------------------------------------------------

      END MODULE TRANSFER_MOD
