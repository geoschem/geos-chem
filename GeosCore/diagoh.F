!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      USE DIAG_MOD, ONLY: AD43, LTNO, LTOH, LTNO2, LTHO2, LTNO3
      USE DIAG_MOD, ONLY: LTO1D, LTO3P

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_O3_MOD    ! SAVEOH, SAVENO
      USE CMN_DIAG_MOD  ! Diagnostic switches & arrays

      IMPLICIT NONE
#     include "define.h"

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
!  15 Aug 2013 - S. D. Eastham - Added O1D (for strat-chem)
!  29 Aug 2013 - S. D. Eastham - Added O3P (for strat-chem)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: I,  J,  L 
      REAL*8   :: OH, NO, HO2, NO2, NO3
      REAL*8   :: O1D,O3P

      !=================================================================
      ! DIAGOH begins here!
      !
      ! ND43 diagnostic: Save OH, HO2, NO3  between HR1_OH and HR2_OH
      !                  Save NO, NO2 between times HR1_NO and HR2_NO
      !
      ! Store the following chemical diagnostics into the AD43 array:
      !    AD43(:,:,:,1) = OH    [molec/cm3]
      !    AD43(:,:,:,2) = NO    [v/v]
      !    AD43(:,:,:,3) = HO2   [v/v]
      !    AD43(:,:,:,4) = NO2   [v/v]
      !    AD43(:,:,:,5) = NO3   [v/v]
      !    AD43(:,:,:,6) = O1D   [molec/cm3]
      !    AD43(:,:,:,7) = O3P   [molec/cm3]
      !=================================================================
      IF ( ND43 > 0 ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, OH, NO, HO2, NO2, NO3, O1D, O3P )  
         DO L = 1, LD43
         DO J = 1, JJPAR
         DO I = 1, IIPAR
  
            ! Save OH as AD43(:,:,:,1)
            OH            = SAVEOH(I,J,L)  * LTOH(I,J)
            AD43(I,J,L,1) = AD43(I,J,L,1)  + OH

            ! Save NO as AD43(:,:,:,2)
            NO            = SAVENO(I,J,L)  * LTNO(I,J)
            AD43(I,J,L,2) = AD43(I,J,L,2)  + NO

            ! Save HO2 as AD43(:,:,:,3)
            HO2           = SAVEHO2(I,J,L) * LTHO2(I,J)
            AD43(I,J,L,3) = AD43(I,J,L,3)  + HO2

            ! Save NO2 as AD43(:,:,:,4)
            NO2           = SAVENO2(I,J,L) * LTNO2(I,J)
            AD43(I,J,L,4) = AD43(I,J,L,4)  + NO2

            ! Save NO3 as AD43(:,:,:,5)
            NO3           = SAVENO3(I,J,L) * LTNO3(I,J)
            AD43(I,J,L,5) = AD43(I,J,L,5)  + NO3

#if defined( UCX )
            ! Save O1D as AD43(:,:,:,6)
            O1D           = SAVEO1D(I,J,L) * LTO1D(I,J)
            AD43(I,J,L,6) = AD43(I,J,L,6)  + O1D

            ! Save O3P as AD43(:,:,:,7)
            O3P           = SAVEO3P(I,J,L) * LTO3P(I,J)
            AD43(I,J,L,7) = AD43(I,J,L,7)  + O3P
#endif

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF  

      END SUBROUTINE DIAGOH
!EOC
