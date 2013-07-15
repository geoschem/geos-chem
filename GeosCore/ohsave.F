!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ohsave
!
! !DESCRIPTION: Subroutine OHSAVE stores the concentrations of OH, HO2, 
!  NO, NO2, and NO3 for the ND43 diagnostic.  Also the O3/Ox, NO/NOx and 
!  NO2/NOx fractions are computed and returned to the calling program. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OHSAVE( N_TRACERS, XNUMOL,  STT,    FRACO3,  
     &                   FRACNO,    FRACNO2, SAVEOH, SAVEHO2, 
     &                   SAVENO,    SAVENO2, SAVENO3,SAVEO1D )
!
! !USES:
!
      USE COMODE_MOD,   ONLY : AIRDENS, CSPEC, JLOP, T3, VOLUME
      USE DIAG_MOD,     ONLY : DIAGCHLORO
      USE TRACERID_MOD, ONLY : IDTOX, IDTNOX, IDO3,  IDNO
      USE TRACERID_MOD, ONLY : IDNO2, IDOH,   IDHO2, IDNO3  
      USE TRACERID_MOD, ONLY : IDO1D

      USE CMN_SIZE_MOD     ! Size parameters
      USE COMODE_LOOP_MOD     ! VOLUME, CSPEC, NPVERT, NLAT, NLONG

      IMPLICIT NONE
#     include "define.h"

!
! !INPUT PARAMETERS: 
!
      ! Number of tracers in XNUMOL and STT
      INTEGER, INTENT(IN) :: N_TRACERS

      ! Array of molec/kg for each tracer
      REAL*8, INTENT(IN)  :: XNUMOL(N_TRACERS)

      ! Array containing CTM tracers
      REAL*8, INTENT(IN)  :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS) 
!
! !OUTPUT PARAMETERS:
!
      ! Array of O3/Ox fractions 
      REAL*8, INTENT(OUT) :: FRACO3(IIPAR,JJPAR,LLPAR)

      ! Array of NO/NOx fractions 
      REAL*8, INTENT(OUT) :: FRACNO(IIPAR,JJPAR,LLPAR)

      ! Array of NO2/NOx fractions
      REAL*8, INTENT(OUT) :: FRACNO2(IIPAR,JJPAR,LLPAR)

      ! Array of OH  concentrations [molec/cm3]
      REAL*8, INTENT(OUT) :: SAVEOH(IIPAR,JJPAR,LLPAR)

      ! Array of HO2 concentrations [v/v]
      REAL*8, INTENT(OUT) :: SAVEHO2(IIPAR,JJPAR,LLPAR)

      ! Array of NO  concentrations [v/v]
      REAL*8, INTENT(OUT) :: SAVENO(IIPAR,JJPAR,LLPAR)     

      ! Array of NO2 concentrations [v/v]
      REAL*8, INTENT(OUT) :: SAVENO2(IIPAR,JJPAR,LLPAR)

      ! Array of NO3 concentrations [v/v]
      REAL*8, INTENT(OUT) :: SAVENO3(IIPAR,JJPAR,LLPAR)

      ! Array of O1D concentrations (molec/cm3)
      REAL*8, INTENT(OUT) :: SAVEO1D(IIPAR,JJPAR,LLPAR)
!
! !REVISION HISTORY: 
!  27 Feb 2002 - R. Yantosca - Initial version
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
!  (8 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (9 ) Reset FRAC* and SAVE* arrays, so that we don't carry dubious data
!        over from boxes that used to be in the tropopause but aren't anymore.
!        (phs, 1/19/07)
!  15 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L, JLOOP   ! (bmy, 7/20/04)
      REAL*8  :: TEMPOX, TEMPNOX  !, KCLO, XLOSS, XOHMASS

      !=================================================================
      ! OHSAVE begins here!
      !
      ! Save info on ozone, OH, and NO concentrations
      ! for consistency with the old method of doing O3, we'll archive
      ! the fraction O3/Ox, and the fraction NO/NOx
      !=================================================================

      ! Reset because of variable tropopause.  Ensure that data for boxes 
      ! that once were in the troposphere, and are not anymore, do not
      ! carry dubious data over. (phs, 1/19/07)
      FRACO3  = 0d0
      FRACNO  = 0d0
      FRACNO2 = 0d0
      SAVEOH  = 0d0
      SAVEHO2 = 0d0
      SAVENO  = 0d0
      SAVENO2 = 0d0
      SAVENO3 = 0d0
      SAVEO1D = 0d0

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

#if defined( GMISTRAT )
         ! O1D concentration [molec/cm3]
         SAVEO1D(I,J,L) = CSPEC(JLOOP,IDO1D)
#endif

 350  CONTINUE
 360  CONTINUE
 370  CONTINUE
!$OMP END PARALLEL DO

      END SUBROUTINE OHSAVE
!EOC
