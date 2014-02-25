!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: diagoh
!
! !DESCRIPTION: Subroutine DIAGOH saves chemical diagnostic quantities for 
!  the ND43 chemical diagnostics.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAGOH
!
! !USES:
!
      USE DIAG_MOD, ONLY: AD43, LTOH, LTHO2

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_O3_MOD    ! SAVEOH
      USE CMN_DIAG_MOD  ! Diagnostic switches & arrays

      IMPLICIT NONE
! 
! !REVISION HISTORY: 
!  01 May 1998 - R. Yantosca - Initial version
!  (1 ) Now use F90 syntax for declarations (bmy, 3/29/99)
!  (2 ) Cosmetic changes (bmy, 3/29/99)
!  (3 ) AD43 and DIAGCHLORO are now declared allocatable in "diag_mod.f". 
!        Also eliminate obsolete code. (bmy, 11/29/99)
!  (4 ) LTNO, LTOH are now allocatable arrays in "diag_mod.f" (bmy, 3/17/00)
!  (5 ) Don't save OH into STT(:,:,:NTRACER+2) anymore.  The SAVEOH 
!        array is now used to save OH concentrations for diagnostics.  
!        Also revised out-of-date comments. (bmy, 4/24/00)
!  (6 ) Also save out NO2 and HO2 for use w/ the ND43 diagnostic.  
!        Now also reference LTNO2, LTHO2 arrays from "diag_mod.f".
!        Updated comments, cosmetic changes. (rvm, bmy, 2/27/02)
!  (7 ) Removed obsolete reference to DIAGCHLORO (bmy, 8/2/02)
!  (8 ) Now save NO3 [molec/cm3] as AD43(:,:,:,5) (bmy, 1/13/03)
!  (9 ) Corrected typo in comments (bmy, 8/10/09)
!  15 Sep 2010 - R. Yantosca - Added ProTeX headers
!  29 Mar 2013 - M. Payer    - Removed NO, NO2, and NO3 from ND43. These are
!                              now tracers.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: I,  J,  L 
      REAL*8   :: OH, HO2

      !=================================================================
      ! DIAGOH begins here!
      !
      ! ND43 diagnostic: Save OH, HO2 between HR1_OH and HR2_OH
      !
      ! Store the following chemical diagnostics into the AD43 array:
      !    AD43(:,:,:,1) = OH    [molec/cm3]
      !    AD43(:,:,:,3) = HO2   [v/v]
      !=================================================================
      IF ( ND43 > 0 ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, OH, HO2 )  
         DO L = 1, LD43
         DO J = 1, JJPAR
         DO I = 1, IIPAR
  
            ! Save OH as AD43(:,:,:,1)
            OH            = SAVEOH(I,J,L)  * LTOH(I,J)
            AD43(I,J,L,1) = AD43(I,J,L,1)  + OH

            ! Save HO2 as AD43(:,:,:,3)
            HO2           = SAVEHO2(I,J,L) * LTHO2(I,J)
            AD43(I,J,L,3) = AD43(I,J,L,3)  + HO2

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF  

      END SUBROUTINE DIAGOH
!EOC
