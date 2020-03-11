!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: henry_mod.F90
!
! !DESCRIPTION: Module HENRY\_MOD contains routines to calculate
! The dimensionless liquid-over-gas Henry constant KH as well as the
! effective Henry constant HEFF, which accounts for hydrolysis.
!\\
!\\
! KH   = K0 * exp ( CR * (1/T - 1/Tref)  ) * R * T / 101.325
!\\
!\\
! HEFF = KH * ( 1 + 10\^(pH-pKa) )
!\\
!\\
! where K0 is the value of KH at standard conditions [M/atm], CR is the
! temperature dependency of KH [K], T is the temperature in Kelvin, Tref
! is the reference temperature (298.15 K), and R is the universal gas
! constant R = 8.314 JK-1mol-1.
!\\
!\\
! References:
! \begin{itemize}
! \item Sander, R: Compilation of Henry's law constant for inorganic and
! organic species of potential importance in environmental chemistry, 1999.
! \item  http://www.mpch-mainz.mpg.de/~sander/res/henry.html.
! \end{itemize}
!
! !INTERFACE:
!
MODULE HENRY_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: CALC_KH
  PUBLIC :: CALC_HEFF
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  REAL*8, PARAMETER :: TREF = 298.15d0      ! [K          ]
  REAL*8, PARAMETER :: R    = 8.3144598d0   ! [J K-1 mol-1]
  REAL*8, PARAMETER :: ATM  = 101.325d0     ! [mPa (!)    ]

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_kh
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
! !INPUT PARAMETERS:
!
    REAL*8,  INTENT(IN)    :: thisK0  ! [M/atm]
    REAL*8,  INTENT(IN)    :: thisCR  ! [-d ln kH / d(1/T) ]
    REAL*8,  INTENT(IN)    :: TK      ! Temperature [K]
!
! !OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(OUT)   :: KH      ! Henry liquid over gas constant [-]
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT) :: RC      ! Error handling
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! CALC_KH begins here!
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
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_heff
!
! !DESCRIPTION: Function CALC\_HEFF calculates the effective Henry
! constant taking into account hydrolysis effects. For instance, aqueous
! HBr will dissociate to Br- (and H+), and often we are interested in the
! total Br concentration in water, i.e. Brx = HBr + Br-, with the equilibrium
! concentration of Brx and gaseous HBr being the effective Henry
! constant. This function provides the correction factor to calculate
! the effective Henry constant from the 'regular' Henry constant
! (above):
!\\
!\\
! Heff = KH * CORR
!\\
!\\
! The correction term is derived as following:
!\\
!\\
! The regular Henry constant is:   H    = HA(liq) / HA(g)
! The effective Henry constant is: Heff = ( HA(liq) + A(liq)) / HA(g)
! Equilibrium between HA and A is: pH = pK + log ( A/HA )
! A(liq) hence becomes: HA(liq) * 10**(pH-pK)
!\\
!\\
! ==> Heff = ( HA(liq) * ( 1 + 10**(pH-pK) ) / HA(g)
!          = HA(liq) / HA(g) * ( 1 + 10**(pH-pK) )
!          = H * ( 1 + 10**(pH-pK) )
!\\
!\\
! ==> CORR = 1 + 10**(pH-pK)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_HEFF ( thispKa, PH, KH, HEFF, RC )
!
! !INPUT PARAMETERS:
!
    REAL*8,  INTENT(IN)    :: thispKa  ! pKa value [-]
    REAL*8,  INTENT(IN)    :: PH       ! PH value [-]
    REAL*8,  INTENT(IN)    :: KH       ! gas/aq Henry constant [-]
!
! !OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(OUT)   :: HEFF     ! effective gas/aq constant [-]
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT) :: RC       ! for error handling
!
! !REMARKS:
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
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! CALC_HEFF begins here!
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
END MODULE HENRY_MOD
