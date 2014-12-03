!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: henry_coeffs.F
!
! !DESCRIPTION: Module HENRY\_COEFFS contains Henry coefficients for a
! suite of species. 
!
! !INTERFACE: 
!
      MODULE HENRY_COEFFS
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)
 
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: GET_HENRY_CONSTANT
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller   - Initial version
!  02 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!

      ! Constants to be used for each tracer
      INTEGER,           PARAMETER :: N = 4

      CHARACTER(LEN=31), PARAMETER :: NAME(N) = (/ 'H2O2',   ! Jacob et al. 2000 
     &                                             'CH3I',   ! Moore et al. 1995
     &                                             'DMS',    ! De Bruyn et al. 1995
     &                                             'ACET' /) ! Benkelberg et al. 1995

      REAL(fp),          PARAMETER :: K0(N)   = (/ 8.3e+4_fp,  
     &                                             1.4e-1_fp,
     &                                             4.8e-1_fp, 
     &                                             27e+0_fp   /)

      REAL(fp),          PARAMETER :: CR(N)   = (/ 7.4e+3_fp,  
     &                                             4.3e+3_fp,
     &                                             3.1e+3_fp,  
     &                                             5300e+0_fp /) 

      REAL(fp),          PARAMETER :: PKa(N)  = (/ -999e+0_fp, 
     &                                             -999e+0_fp,
     &                                             -999e+0_fp, 
     &                                             -999e+0_fp /) 

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_HENRY_CONSTANT
!
! !DESCRIPTION: Subroutine GET\_HENRY\_CONSTANT returns the Henry
! constants for the given species by name.
! K0 is the Henry constant @ 298K in [M/atm, liquid over gas], and CR is 
! the temperature dependency -d ln(kH) / d(1/T), following Sander (1999):
! http://www.mpch-mainz.mpg.de/~sander/res/henry.html
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_HENRY_CONSTANT ( NM, thisK0, thisCR, thispKa, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
      CHARACTER(LEN=*), INTENT(IN   )  :: NM
      REAL(fp),           INTENT(  OUT)  :: thisK0      
      REAL(fp),           INTENT(  OUT)  :: thisCR     
      REAL(fp),           INTENT(  OUT)  :: thispKa    
      INTEGER,          INTENT(INOUT)  :: RC    
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      INTEGER :: I

      !=================================================================
      ! GET_HENRY_CONSTANT begins here!
      !=================================================================

      ! Init
      RC      = 0
      thisK0  = -999.0e+0_fp
      thisCR  = -999.0e+0_fp
      thisPKa = -999.0e+0_fp

      ! Loop over all species
      DO I = 1,N
         IF ( TRIM(NAME(I)) == TRIM(NM) ) THEN
            thisK0  = K0(I)
            thisCR  = CR(I)
            thisPKa = pKA(I)
            RETURN
         ENDIF
      ENDDO
 
      END SUBROUTINE GET_HENRY_CONSTANT 
!EOC
      END MODULE HENRY_COEFFS
!EOM
