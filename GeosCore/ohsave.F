!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ohsave
!
! !DESCRIPTION: Subroutine OHSAVE stores the concentrations of OH and HO2 
!  for the ND43 diagnostic.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OHSAVE( SAVEOH, SAVEHO2 )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_MOD,      ONLY : AIRDENS, CSPEC, JLOP
      USE COMODE_LOOP_MOD
      USE TRACERID_MOD,    ONLY : IDOH, IDHO2

      IMPLICIT NONE
!
! !OUTPUT PARAMETERS:
!
      ! Array of OH  concentrations [molec/cm3]
      REAL*8, INTENT(OUT) :: SAVEOH(IIPAR,JJPAR,LLPAR)

      ! Array of HO2 concentrations [v/v]
      REAL*8, INTENT(OUT) :: SAVEHO2(IIPAR,JJPAR,LLPAR)

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
!  14 Mar 2013 - M. Payer    - FRACO3, FRACNO, and FRACNO2 are no longer needed
!                              because O3, NO, and NO2 are now tracers.
!  29 Mar 2013 - M. Payer    - Removed SAVENO, SAVENO2, SAVENO3.
!  31 May 2013 - R. Yantosca - Extra cleanup, remove N_TRACERS, XNUMOL, STT
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L, JLOOP   ! (bmy, 7/20/04)

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
      SAVEOH  = 0d0
      SAVEHO2 = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, JLOOP )
!$OMP+SCHEDULE( DYNAMIC )
      DO 370 L = 1, NPVERT
      DO 360 J = 1, NLAT
      DO 350 I = 1, NLONG

         ! 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Cycle if this isn't a valid SMVGEAR gridbox
         IF ( JLOOP == 0 ) GOTO 350

         ! OH concentration [molec/cm3]
         SAVEOH(I,J,L)  = CSPEC(JLOOP,IDOH)

         ! HO2 concentration [v/v] 
         SAVEHO2(I,J,L) = CSPEC(JLOOP,IDHO2) / AIRDENS(JLOOP)  

 350  CONTINUE
 360  CONTINUE
 370  CONTINUE
!$OMP END PARALLEL DO

      END SUBROUTINE OHSAVE
!EOC
