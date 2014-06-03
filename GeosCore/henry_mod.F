!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: HENRY_MOD
!
! !DESCRIPTION: Module HENRY\_MOD contains routines to calculate
! The dimensionless liquid-over-gas Henry constant KH as well as the
! effective Henry constant HEFF, which accounts for hydrolysis.
!
! KH   = K0 * exp ( CR * (1/T - 1/Tref)  ) * R * T / 101.325 
! HEFF = KH * ( 1 + 10^(pH-pKa) )
!
! where K0 is the value of KH at standard conditions [M/atm], CR is the
! temperature dependency of KH [K], T is the temperature in Kelvin, Tref 
! is the reference temperature (298.15 K), and R is the universal gas
! constant R = 8.314 JK-1mol-1.
!
! Note that K0, CR and pKa are defined below.
!
! !REFERENCES:
! Sander, R: Compilation of Henry's law constant for inorganic and
! organic species of potential importance in environmental chemistry,
! 1999. http://www.mpch-mainz.mpg.de/~sander/res/henry.html.
!\\
!\\
! !INTERFACE: 
!
      MODULE HENRY_MOD
!
! !USES:
! 
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CALC_KH
      PUBLIC :: CALC_HEFF
      PUBLIC :: GET_HENRY_CONSTANTS
      PUBLIC :: CLEANUP_HENRY 
      PUBLIC :: INIT_HENRY
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE VARIABLES:
!
      ! Constants to be used for each tracer
      REAL*8, ALLOCATABLE, TARGET  :: K0(:) 
      REAL*8, ALLOCATABLE, TARGET  :: CR(:) 
      REAL*8, ALLOCATABLE, TARGET  :: PKa(:) 

      ! Parameter
      REAL*8,  PARAMETER    :: TREF = 298.15d0    ![K]
      REAL*8,  PARAMETER    :: R    = 8.3144621d0 ![JK-1mol-1]
      REAL*8,  PARAMETER    :: ATM  = 101.325     ![mPa (!)]

      CONTAINS
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CALC_KH
!
! !DESCRIPTION: Function CALC\_KH calculates the liquid over gas
! dimensionless Henry constant for the given tracer and temperature.
!
! Reference: http://www.mpch-mainz.mpg.de/~sander/res/henry.html
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_KH ( thisK0, thisCR, TK, KH, RC ) 
!
! !USES:
! 
!
! !INPUT ARGUMENTS:
!
      REAL*8,  INTENT(IN)         :: thisK0     ! [M/atm] 
      REAL*8,  INTENT(IN)         :: thisCR     ! [-d ln kH / d(1/T) ]
      REAL*8,  INTENT(IN)         :: TK         ! Temperature [K]
!
! !OUTPUT ARGUMENTS:
!
      REAL*8,  INTENT(OUT)        :: KH      ! Henry liquid over gas constant [-]
      INTEGER, INTENT(INOUT)      :: RC      ! Error handling 
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!

      !=================================================================
      ! CALC_HENRY begins here!
      !=================================================================

      ! Assume success
      RC = 0

      ! Error if not defined
      IF ( thisK0 == 0d0 ) THEN
         RC = -999
         KH = -999
         RETURN
      ENDIF 

      ! Calculate Henry coefficient for given temperature
      KH = thisK0 * exp( thisCR * (1/TK - 1/TREF) ) * R * TK / ATM

      END SUBROUTINE CALC_KH
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CALC_HEFF
!
! !DESCRIPTION: Function CALC_HEFF calculates the effective Henry
! constant taking into account hydrolysis effects. For instance, aqueous 
! HBr will dissociate to Br- (and H+), and often we are interested in the 
! total Br concentration in water, i.e. Brx = HBr + Br-, with the equilibrium
! concentration of Brx and gaseous HBr being the effective Henry
! constant. This function provides the correction factor to calculate
! the effective Henry constant from the 'regular' Henry constant
! (above):
! Heff = KH * CORR
! 
! The correction term is derived as following:
!
! The regular Henry constant is:   H    = HA(liq) / HA(g)
! The effective Henry constant is: Heff = ( HA(liq) + A(liq)) / HA(g)
! Equilibrium between HA and A is: pH = pK + log ( A/HA )
! A(liq) hence becomes: HA(liq) * 10**(pH-pK)
!
! ==> Heff = ( HA(liq) * ( 1 + 10**(pH-pK) ) / HA(g)
!          = HA(liq) / HA(g) * ( 1 + 10**(pH-pK) )
!          = H * ( 1 + 10**(pH-pK) )
!
! ==> CORR = 1 + 10**(pH-pK)
!
! !NOTES:
! It should be noted that the correction term calculated here is from a
! 'acid perspective', i.e. for compounds with the acid being in the
! gaseous phase. The correction term reads 1 + 10**(-pH+pK) for
! compounds with the base in the gas phase (e.g. ammonia).
!
! The correction term becomes more complicated for compounds with more
! than two equilibrium compounds that are relevant under the current
! conditions (e.g. CO2).
!
! We ignore any temperature dependencies of pKa for now.
!
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_HEFF ( thispKa, PH, KH, HEFF, RC ) 
!
! !USES:
! 
!
! !INPUT ARGUMENTS:
!
      REAL*8,  INTENT(IN)         :: thispKa    ! pKa value [-] 
      REAL*8,  INTENT(IN)         :: PH         ! PH value [-] 
      REAL*8,  INTENT(IN)         :: KH         ! gas/aq Henry constant [-] 
!
! !OUTPUT ARGUMENTS:
!
      REAL*8,  INTENT(OUT)        :: HEFF  ! effective gas/aq constant [-] 
      INTEGER, INTENT(INOUT)      :: RC    ! for error handling 
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!

      !=================================================================
      ! PH_CORR begins here!
      !=================================================================

      ! Assume success
      RC = 0

      ! Calculate correction term.
      IF ( thispKa > 0d0 ) THEN
         HEFF = KH * ( 1d0 + 10d0**( pH - thispKa ) )
      ELSE
         HEFF = KH
      ENDIF

      END SUBROUTINE CALC_HEFF
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_HENRY_CONSTANTS
!
! !DESCRIPTION: Subroutine GET\_HENRY\_CONSTANTS returns the Henry
! constants for the given GC tracer. 
! K0 is the Henry constant @ 298K in [M/atm, liquid over gas], and CR is 
! the temperature dependency -d ln(kH) / d(1/T), following Sander (1999):
! http://www.mpch-mainz.mpg.de/~sander/res/henry.html
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_HENRY_CONSTANTS ( ID, thisK0, thisCR, thispKa, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(IN   )  :: ID       
      REAL*8,           INTENT(  OUT)  :: thisK0      
      REAL*8,           INTENT(  OUT)  :: thisCR     
      REAL*8,           INTENT(  OUT)  :: thispKa    
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

      !=================================================================
      ! GET_HENRY_CONSTANTS begins here!
      !=================================================================

      ! Init
      thisK0  = 0.0d0
      thisCR  = 0.0d0
      thisPKa = 0.0d0

      ! Tracer ID not defined:
      IF ( ID < 0 ) THEN
         RC = 0
 
      ! Tracer ID defined:
      ELSEIF ( ID <= SIZE(K0) ) THEN
         thisK0  = K0(ID)
         thisCR  = CR(ID)
         thisPKa = pKA(ID)
         RC = 0
      ELSE
         RC = -999
      ENDIF

      END SUBROUTINE GET_HENRY_CONSTANTS 
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_HENRY
!
! !DESCRIPTION: Subroutine INIT\_HENRY initializes all module variables.
! K0 is the Henry constant @ 298K in [M/atm, liquid over gas], and CR is 
! the temperature dependency -d ln(kH) / d(1/T), following Sander (1999):
! http://www.mpch-mainz.mpg.de/~sander/res/henry.html
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_HENRY ( SpecNames, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
      CHARACTER(LEN=*), INTENT(IN   )  :: SpecNames(:)
      INTEGER,          INTENT(INOUT)  :: RC         ! Return status
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      INTEGER  :: I, N, AS

      !=================================================================
      ! INIT_HENRY begins here!
      !=================================================================

      ! Set return status to wrong before having done all allocations
      RC = -999

      ! Get number of tracers
      N = SIZE( SpecNames ) 

      ! Allocate vectors and initialize values
      ALLOCATE( K0(N), STAT=AS )
      IF ( AS /= 0 ) RETURN 
      K0(:) = 0d0

      ALLOCATE( CR(N), STAT=AS )
      IF ( AS /= 0 ) RETURN 
      CR(:) = 0d0

      ALLOCATE( PKa(N), STAT=AS )
      IF ( AS /= 0 ) RETURN 
      PKa(:) = 0d0

      ! Loop over all tracers and assign constants as provided below
      DO I = 1,N 

         ! -------------------------------------------------
         ! H2O2
         ! -------------------------------------------------
         IF (TRIM(SpecNames(I)) == 'H2O2' ) THEN
            K0(I) = 8.3d4    ! Jacob et al., 2000
            CR(I) = 7.4d3    ! Jacob et al., 2000
            CYCLE
         ENDIF      

         ! -------------------------------------------------
         ! CH3I
         ! -------------------------------------------------
         IF (TRIM(SpecNames(I)) == 'CH3I' ) THEN
            K0(I) = 1.4d-1  ! Moore et al., 1995
            CR(I) = 4.3d3   ! Moore et al., 1995
            CYCLE
         ENDIF

         ! -------------------------------------------------
         ! DMS
         ! -------------------------------------------------
         IF (TRIM(SpecNames(I)) == 'DMS' ) THEN
            K0(I) = 4.8d-1  ! De Bruyn et al. [1995] 
            CR(I) = 3.1d3   ! De Bruyn et al. [1995]
            CYCLE
         ENDIF

         ! -------------------------------------------------
         ! Acetone
         ! -------------------------------------------------
         IF (TRIM(SpecNames(I)) == 'ACET' ) THEN
            K0(I) = 27d0     ! Benkelberg et al. [1995] 
            CR(I) = 5300d0   ! Benkelberg et al. [1995]
            CYCLE
         ENDIF

      ENDDO !I

      ! Return w/ success
      RC = 0

      END SUBROUTINE INIT_HENRY
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP_HENRY
!
! !DESCRIPTION: Subroutine CLEANUP\_HENRY deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_HENRY
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! CLEANUP_HENRY begins here!
      !=================================================================

      ! Nullify all pointer
      IF ( ALLOCATED ( K0 ) ) DEALLOCATE( K0 )
      IF ( ALLOCATED ( CR ) ) DEALLOCATE( CR )

      END SUBROUTINE CLEANUP_HENRY
!EOC
      END MODULE HENRY_MOD
!EOM
