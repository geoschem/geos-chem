! $Id: ohsave.f,v 1.4 2004/12/02 21:48:38 bmy Exp $
      SUBROUTINE OHSAVE( N_TRACERS, XNUMOL,  STT,    FRACO3,  
     &                   FRACNO,    FRACNO2, SAVEOH, SAVEHO2, 
     &                   SAVENO,    SAVENO2, SAVENO3 )
!
!******************************************************************************
!  Subroutine OHSAVE stores the concentrations of OH, HO2, NO, NO2, and NO3
!  for the ND43 diagnostic.  Also the O3/Ox, NO/NOx and NO2/NOx fractions
!  are computed and returned to the calling program. (bmy, 2/27/02, 7/20/04) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N_TRACERS (INTEGER) : Number of tracers in XNUMOL and STT
!  (2 ) XNUMOL    (REAL*8 ) : Array of molec/kg for each tracer
!  (3 ) STT       (REAL*8 ) : Array containing CTM tracers
!
!  Arguments as Output:
!  ============================================================================
!  (4 ) FRACO3    (REAL*8 ) : Array of O3/Ox   fractions 
!  (5 ) FRACNO    (REAL*8 ) : Array of NO/NOx  fractions 
!  (6 ) FRACNO2   (REAL*8 ) : Array of NO2/NOx fractions
!  (7 ) SAVEOH    (REAL*8 ) : Array of OH  concentrations [molec/cm3]
!  (8 ) SAVEHO2   (REAL*8 ) : Array of HO2 concentrations [v/v]
!  (9 ) SAVENO    (REAL*8 ) : Array of NO  concentrations [v/v]
!  (10) SAVENO2   (REAL*8 ) : Array of NO2 concentrations [v/v]
!  (11) SAVENO3   (REAL*8 ) : Array of NO3 concentrations [v/v]
!
!  NOTES:
!  (1 ) Original code from lwh, gmg, djj, jyl, etc, 1990's.  Modified for
!        GEOS-CHEM by Bob Yantosca et al.
!  (2 ) Added comment header and F90 declaration syntax.  Also now specify
!        the units of each variable for clarity. 
!  (3 ) Deleted NTRACER, it is not used.  Also added FRACNO2 and SAVEHO2
!        variables.  Updated comments, cosmetic changes (rvm, bmy, 2/27/02)
!  (4 ) Bug fix: swap the order of the lines where TMPNOX is computed.
!        Also deleted obsolete code from 2/02. (bmy, 7/31/02)
!  (5 ) Now reference IDTOX, IDTNOX, etc from "tracerid_mod.f". (1/13/03)
!  (6 ) Added OpenMP parallelization commands (bmy, 8/1/03)
!  (7 ) Now compute quantities for mean OH in "diag_oh_mod.f".  Now also
!        references STT from "tracer_mod.f".  Added N_TRACERS to the arg list.
!        Now dimension args XNUMOL, STT w/ N_TRACERS and not NNPAR. 
!        (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD, ONLY : AIRDENS, CSPEC, JLOP, T3, VOLUME
      USE DIAG_MOD,   ONLY : DIAGCHLORO
      USE TRACERID_MOD

      IMPLICIT NONE

#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"   ! VOLUME, CSPEC, NPVERT, NLAT, NLONG

      ! Arguments
      INTEGER, INTENT(IN) :: N_TRACERS
      REAL*8, INTENT(IN)  :: XNUMOL(N_TRACERS)
      REAL*8, INTENT(IN)  :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS) 
      REAL*8, INTENT(OUT) :: FRACO3(IIPAR,JJPAR,LLPAR)
      REAL*8, INTENT(OUT) :: FRACNO(IIPAR,JJPAR,LLPAR)
      REAL*8, INTENT(OUT) :: FRACNO2(IIPAR,JJPAR,LLPAR)
      REAL*8, INTENT(OUT) :: SAVEOH(IIPAR,JJPAR,LLPAR)
      REAL*8, INTENT(OUT) :: SAVEHO2(IIPAR,JJPAR,LLPAR)
      REAL*8, INTENT(OUT) :: SAVENO(IIPAR,JJPAR,LLPAR)     
      REAL*8, INTENT(OUT) :: SAVENO2(IIPAR,JJPAR,LLPAR)
      REAL*8, INTENT(OUT) :: SAVENO3(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER             :: I, J, L, JLOOP   ! (bmy, 7/20/04)
      REAL*8              :: TEMPOX, TEMPNOX  !, KCLO, XLOSS, XOHMASS

      !=================================================================
      ! DIAGOH begins here!
      !
      ! Save info on ozone, OH, and NO concentrations
      ! for consistency with the old method of doing O3, we'll archive
      ! the fraction O3/Ox, and the fraction NO/NOx
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, JLOOP, TEMPOX, TEMPNOX )
!$OMP+SCHEDULE( DYNAMIC )
      DO 370 L = 1, NPVERT
      DO 360 J = 1, NLAT
      DO 350 I = 1, NLONG

         ! 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Cycle if this isn't a valid SMVGEAR gridbox
         IF ( JLOOP == 0 ) GOTO 350

         ! Total Ox concentration, convert from [kg] to [molec/cm3]
         TEMPOX         = STT(I,J,L,IDTOX)
         TEMPOX         = TEMPOX  * XNUMOL(IDTOX) /VOLUME(JLOOP)

         ! Total NOx concentration, convert from [kg] to [molec/cm3]
         TEMPNOX        = STT(I,J,L,IDTNOX)
         TEMPNOX        = TEMPNOX * XNUMOL(IDTNOX)/VOLUME(JLOOP)

         ! Ox/O3 fraction [unitless]
         FRACO3(I,J,L)  = CSPEC(JLOOP,IDO3)  / TEMPOX

         ! NO/NOx fraction [unitless]
         FRACNO(I,J,L)  = CSPEC(JLOOP,IDNO)  / TEMPNOX

         ! NO2/NOx fraction [unitless]
         FRACNO2(I,J,L) = CSPEC(JLOOP,IDNO2) / TEMPNOX    

         ! OH concentration [molec/cm3]
         SAVEOH(I,J,L)  = CSPEC(JLOOP,IDOH)

         ! HO2 concentration [v/v] 
         SAVEHO2(I,J,L) = CSPEC(JLOOP,IDHO2) / AIRDENS(JLOOP)  

         ! NO concentration [v/v]
         SAVENO(I,J,L)  = CSPEC(JLOOP,IDNO)  / AIRDENS(JLOOP)

         ! NO2 concentration [v/v]
         SAVENO2(I,J,L) = CSPEC(JLOOP,IDNO2) / AIRDENS(JLOOP)  

         ! NO3 concentration [v/v]
         SAVENO3(I,J,L) = CSPEC(JLOOP,IDNO3) / AIRDENS(JLOOP)

 350  CONTINUE
 360  CONTINUE
 370  CONTINUE
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE OHSAVE
