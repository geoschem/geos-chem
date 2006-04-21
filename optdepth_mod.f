! $Id: optdepth_mod.f,v 1.5 2006/04/21 15:40:04 bmy Exp $
      MODULE OPTDEPTH_MOD
!
!******************************************************************************
!  Module OPTDEPTH_MOD contains routines to compute optical depths for GEOS-1
!  GEOS-STRAT, GEOS-2, and GEOS-3 met data sets. (bmy, 8/15/01, 5/28/02)
! 
!  Module Routines:
!  ============================================================================
!  (1 ) OD_GEOS1_GEOSS : Computes optical depths for GEOS-1 or GEOS-STRAT
!  (2 ) OD_GEOS3_GEOS4 : Computes optical depths for GEOS-2 or GEOS-3 
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
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "optdepth_mod.f"
      !=================================================================

      ! PRIVATE module routines
      PRIVATE OD_GEOS1_GEOSS
      PRIVATE OD_GEOS3_GEOS4

      !=================================================================
      ! MODULE INTERFACES -- "bind" two or more routines with different
      ! argument types or # of arguments under one unique name
      !================================================================= 
      INTERFACE OPTDEPTH
         MODULE PROCEDURE OD_GEOS1_GEOSS
         MODULE PROCEDURE OD_GEOS3_GEOS4
      END INTERFACE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE OD_GEOS1_GEOSS( NVERT, CLMO, CLRO, DELP, T, OPTD ) 
!
!******************************************************************************
!  Subroutine OD_GEOS1_GEOSS computes the DAO cloud optical depth based 
!  on the temperature and cloud fractions of grid box (I,J,L).  
!  (bmy, 6/5/98, 10/24/01)
!    
!  Arguments as input: 
!  ===========================================================================
!  (1 ) NVERT (INTEGER) : number of sigma levels to use      [unitless]
!  (2 ) CLMO  (REAL*8 ) : DAO maximum overlap cloud fraction [unitless]
!  (3 ) CLRO  (REAL*8 ) : DAO random  overlap cloud fraction [unitless]
!  (4 ) DELP  (REAL*8 ) : Delta-P thickness of grid box      [mb      ]
!  (5 ) T     (REAL*8 ) : DAO temperature of grid box        [K       ] 
!
!  Arguments as output:
!  ===========================================================================
!  (6 ) OPTD  (REAL*8 ) : DAO optical depth at (I,J,L)       [unitless]
!
!  Reference:
!  =======================================================================
!  (1 ) L. Takacs, A. Molod, T. Wang, "Technical Report Series on Global
!        Modeling and Data Assimilation": NASA Technical Memorandum
!        104606, Vol 1., Sept. 1994, pp. 14-15.
!
!  NOTES:
!  (1 ) OPTDEPTH is written in Fixed-Form Fortran 90.  Also use F90 syntax
!       for declarations (bmy, 3/29/99)
!  (2 ) MAIN now passes the Harvard CTM variable for temperature
!        T(IGCMPAR,JGCMPAR,LLPAR) to OPTDEPTH.  Use window offsets I+I0,
!        J+J0 when accessing T.  This eliminates the need of having two
!        separate arrays to contain temperature values.
!  (3 ) CLMO, CLRO, DELP, TAUCLD, CLDTOT are dimensioned (LLPAR,IIPAR,JJPAR) 
!        to maximize loop efficiency for processing an (I,J) column layer by 
!        layer.
!  (4 ) Eliminated variables LN21 and X -- these are obsolete (bmy, 8/9/01)
!  (5 ) Now parallelize I-J do loops (bmy, 8/15/01)
!  (6 ) T(IREF,JREF,L) is now T(I,J,L).  Also, T is now dimensioned 
!        (IIPAR,JJPAR,LLPAR).  Remove IREF,JREF -- these are obsolete.
!        (bmy, 9/27/01)
!  (7 ) Removed obsolete code from 9/01 (bmy, 10/24/01)
!******************************************************************************
! 
      ! References to F90 modules
      USE DIAG_MOD, ONLY: AD21

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_DIAG"   ! ND21

      ! Arguments
      INTEGER, INTENT(IN)  :: NVERT
      REAL*8,  INTENT(IN)  :: CLMO(LLPAR,IIPAR,JJPAR)  
      REAL*8,  INTENT(IN)  :: CLRO(LLPAR,IIPAR,JJPAR) 
      REAL*8,  INTENT(IN)  :: DELP(LLPAR,IIPAR,JJPAR)  
      REAL*8,  INTENT(IN)  :: T   (IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: OPTD(LLPAR,IIPAR,JJPAR)  
      
      ! Local Variables
      INTEGER              :: L, I, J
      REAL*8               :: T1, TAURO, TAUMO
      
      !=================================================================
      ! OD_GEOS1_GEOSS begins here!
      !
      ! Optical depth has to be computed as a function of 
      ! temperature and cloud fractions.  See the reference above.
      !================================================================= 
      
      ! TAUMO = the optical depth for maximum overlap cloudiness CLMO (mb^-1) 
      ! Note: TAUMO is the same for all grid boxes (I,J,L) = 0.16
      TAUMO = 0.16d0

      ! Loop over grid boxes in the following order: J - I - L.  This 
      ! allows us to process an (I,J) column layer by layer.  
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, T1, TAURO )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
         DO I = 1, IIPAR

            !=========================================================== 
            ! TAURO = optical depth for random overlap cloudiness 
            ! CLRO (mb^-1), as given by the following formula:
            !
            ! TAURO = { 0.0,                                T <  190.66
            !         { (T - 190.66)^2 * (2.0e-6), 190.66 < T <= 263.16
            !         { (6.95e-3 * T) - 1.82,      263.16 < T <= 273.38 
            !         { 0.08,                      273.38 < T
            !
            ! Where T is the temperature of grid box (I,J,L).
            !
            ! Overall optical depth TAU (km^-1) is given by
            !
            !  TAU(L,I,J) = ( TAUMO        * CLMO(L,I,J) ) + 
            !               ( TAURO(L,I,J) * CLRO(L,I,J) )
            !
            ! Multiply TAU(L,I,J) by DELP(L,I,J) to obtain unitless 
            ! optical depths.
            !
            ! Attach ND21 diagnostics -- map of optical depths and 
            ! cloud fractions (bmy, 9/11/98)
            !
            ! NVERT should be the highest level where we are doing 
            ! chemistry.
            !=========================================================== 
            DO L = 1, NVERT
               T1 = T(I,J,L)

               IF ( T1 <= 190.66 ) THEN
                  TAURO = 0.0d0
                  
               ELSE IF ( T1 >  190.66  .and. 
     &                   T1 <= 263.16 ) THEN
                  TAURO = 2.0d-6 * ( T1 - 190.66d0 )**2
            
               ELSE IF ( T1 >  263.16  .and. 
     &                   T1 <= 273.38 ) THEN
                  TAURO = ( 6.95d-3 * T1 ) - 1.82d0
                  
               ELSE
                  TAURO = 0.08d0
                  
               ENDIF             

               ! Compute unitless optical depths
               OPTD(L,I,J) = DELP(L,I,J) * 
     &            ( ( TAUMO * CLMO(L,I,J) ) + ( TAURO * CLRO(L,I,J) ) )

               ! Save to AD21 array only if ND21 is turned on
               IF ( ND21 > 0 .and. L <= LD21 ) THEN
                  AD21(I,J,L,1) = AD21(I,J,L,1) + OPTD(L,I,J)
                  AD21(I,J,L,2) = AD21(I,J,L,2) + CLMO(L,I,J)
                  AD21(I,J,L,3) = AD21(I,J,L,3) + CLRO(L,I,J)
               ENDIF 
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE OD_GEOS1_GEOSS

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

      END MODULE OPTDEPTH_MOD





