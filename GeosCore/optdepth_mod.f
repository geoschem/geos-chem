! $Id: optdepth_mod.f,v 1.1 2009/09/16 14:06:17 bmy Exp $
      MODULE OPTDEPTH_MOD
!
!******************************************************************************
!  Module OPTDEPTH_MOD contains routines to compute optical depths for GEOS-3
!  GEOS-4, and GCAP met data sets. (bmy, 8/15/01, 8/4/06)
! 
!  Module Routines:
!  ============================================================================
!  (1 ) OD_GEOS3_GEOS4 : Computes optical depths for GEOS-2 or GEOS-3 
!
!  Module Interfaces:
!  ============================================================================
!  (1 ) OPTDEPTH       : Connects routines OD_GEOS1_GEOSS, OD_GEOS2_GEOS3
!
!  GEOS-CHEM modules referenced by optdepth_mod.f
!  ============================================================================
!  (1 ) diag_mod.f     : Module containing GEOS-CHEM diagnostic arrays
!
!  NOTES: 
!  (1 ) Now add parallel DO-loops (bmy, 8/15/01)
!  (2 ) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (3 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Also add MODULE INTERFACES section,
!        since we have an interface here. (bmy, 5/28/02)
!  (4 ) Renamed OD_GEOS2_GEOS_3 to OD_GEOS3_GEOS4.  (bmy, 4/20/05)
!  (5 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "optdepth_mod.f"
      !=================================================================

      ! PRIVATE module routines
      PRIVATE OD_GEOS3_GEOS4

      !=================================================================
      ! MODULE INTERFACES -- "bind" two or more routines with different
      ! argument types or # of arguments under one unique name
      !================================================================= 
      INTERFACE OPTDEPTH
         MODULE PROCEDURE OD_GEOS3_GEOS4
      END INTERFACE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE OD_GEOS3_GEOS4( NVERT, CLDF, OPTDEP, OPTD )
!
!******************************************************************************
!  Subroutine OD_GEOS3_GEOS4 copies the DAO grid box optical depth from
!  the OPTDEP met field array into the OPTD array.  Diagnostics are also
!  archived. (bmy, 8/15/01, 4/20/05)
!   
!  Arguments as input: 
!  ===========================================================================
!  (1 ) NVERT  (INTEGER) : Number of levels to compute Optical Depth fo
!  (2 ) CLDF   (REAL*8 ) : GEOS-3/GEOS-4 3/D cloud fraction [unitless]
!  (3 ) OPTDEP (REAL*8 ) : GEOS-3/GEOS-4 grid box optical depths [unitless]
!
!  Arguments as output:
!  ===========================================================================
!  (4 ) OPTD   (REAL*8 ) : DAO optical depth at grid box (I,J,L) [unitless]
!
!  NOTES:
!  (1 ) Now parallelize I-J DO loops (bmy, 8/15/01)
!  (2 ) Renamed to OD_GEOS3_GEOS4.  Also now saves CLDF in AD21(I,J,L,2)
!        for the ND21 diagnostic (bmy, 4/20/05)
!******************************************************************************
! 
      ! References to F90 modules
      USE DIAG_MOD, ONLY: AD21

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_DIAG"   ! ND21

      ! Arguments
      INTEGER, INTENT(IN)  :: NVERT
      REAL*8,  INTENT(IN)  :: CLDF  (LLPAR,IIPAR,JJPAR)
      REAL*8,  INTENT(IN)  :: OPTDEP(LLPAR,IIPAR,JJPAR)  
      REAL*8,  INTENT(OUT) :: OPTD  (LLPAR,IIPAR,JJPAR)  
      
      ! Local Variables
      INTEGER              :: I, J, L
      
      !=================================================================
      ! OD_GEOS3_GEOS4 begins here!
      !
      ! GEOS-3/GEOS-4 optical depth is stored in the OPTDEP array,
      ! which is read in routine "read_a6" of "dao_read_mod.f".
      !
      ! OPTDEP is archived every 6 hours, nevertheless, each chemistry
      ! timestep we copy this into the OPTD array and archive for the
      ! ND21 diagnostic.  This way the ND21 diagnostic is consistent
      ! with GEOS-1/GEOS-STRAT.
      !
      ! OPTDEP and OPTD are dimensioned (LLPAR,IIPAR,JJPAR) to maximize
      ! loop efficiency for processing an (I,J) column layer by layer.
      !
      ! Now also save CLDTOT to the ND21 diagnostic (bmy, 4/20/05)
      !================================================================= 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      DO L = 1, NVERT   
         
         ! Copy optical depth over from OPTDEP array
         OPTD(L,I,J) = OPTDEP(L,I,J) 
         
         ! Save to AD21 array only if ND21 is turned on
         IF ( ND21 > 0 .and. L <= LD21 ) THEN
            AD21(I,J,L,1) = AD21(I,J,L,1) + OPTD(L,I,J) 
            AD21(I,J,L,2) = AD21(I,J,L,2) + CLDF(L,I,J)
         ENDIF 
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE OD_GEOS3_GEOS4

!------------------------------------------------------------------------------

      ! End of module
      END MODULE OPTDEPTH_MOD





