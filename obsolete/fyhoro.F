!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: FYHORO
!
! !DESCRIPTION: \subsection*{Overview}
!  Function FYHORO returns returns the branching ratio between 
!  HOC2H4O oxidation and dissociation:
!  (1) HOC2H4 --O2--> HO2 + GLYC
!  (2) HOC2H4 ------> HO2 + 2CH2O

!\subsection*{References}
!  \begin{enumerate}
!  \item Orlando et al., 1998: \emph{Laboratory and theoretical study of the 
!         oxyradicals in the OH- and Cl-initiated oxidation of ethene}, 
!        \underline{J. Phys. Chem. A}, \textbf{102}, 8116-8123.
!  \item Orlando et al., 2003: \emph{The atmospheric chemistry of alkoxy 
!         radicals}, \underline{Chem. Rev.}, \textbf{103}, 4657-4689.
!  \end{enumerate}
!
!\\
!\\
! !INTERFACE: 
!
      REAL*8 FUNCTION FYHORO( ZDNUM, TT )
! 
! !USES:
!
      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      ! Air density   [molec/cm3 ]
      REAL*8, INTENT(IN) :: ZDNUM

      ! Temperature   [K         ]
      REAL*8, INTENT(IN) :: TT

!
! !REVISION HISTORY:
!  (1 ) Branching ratio calculation (tmf, 2/6/05).
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8             :: K1, K2, O2DNUM

      !=================================================================
      ! FYHORO begins here!
      !=================================================================
      O2DNUM = ZDNUM * 0.21D0
      K1     = 6.0D-14 * EXP(-550.D0/TT) * O2DNUM
      K2     = 9.5D+13 * EXP(-5988.D0/TT) 

      FYHORO = K1 / (K1 + K2)
     
      ! Return to calling program
      END FUNCTION FYHORO
!EOC 
