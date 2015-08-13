!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: get_popsinfo_mod
!
! !DESCRIPTION: Module GET\_POPSINFO\_MOD contains variables and routines for
!  the GEOS-Chem peristent organic pollutants (POPs) simulation. 
!\\
!\\
! !INTERFACE: 
!
      MODULE GET_POPSINFO_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: GET_POP_TYPE
      PUBLIC :: GET_EMISSFILE
      PUBLIC :: GET_POP_XMW
      PUBLIC :: GET_POP_HSTAR
      PUBLIC :: GET_POP_DEL_Hw
      PUBLIC :: GET_POP_DEL_H
      PUBLIC :: GET_POP_KBC
      PUBLIC :: GET_POP_K_POPP_O3A
      PUBLIC :: GET_POP_K_POPP_O3B
      PUBLIC :: GET_POP_K_POPG_OH
      PUBLIC :: GET_POP_KOA
      PUBLIC :: INIT_POP_PARAMS
!
! !REVISION HISTORY:
!  30 Sep 2012 - C. Pike Thackray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_pop_type
!
! !DESCRIPTION: Function GET\_POP\_TYPE is used to retrieve type of POP.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_POP_TYPE( IN_TYPE )
!
! !INPUT PARAMETERS:
!
      CHARACTER(LEN=3), INTENT(IN) :: IN_TYPE
!
! !REVISION HISTORY:
!  30 Sep 2012 - C. Pike Thackray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=3)             :: GET_POP_TYPE
      CHARACTER(LEN=3), SAVE       :: POP_TYPE
      LOGICAL, SAVE                :: IS_SET

      !=================================================================
      ! GET_POP_TYPE begins here
      !=================================================================

      IF ( IS_SET ) THEN
         GET_POP_TYPE = POP_TYPE
         RETURN
      ENDIF

      POP_TYPE = IN_TYPE
      IS_SET = .TRUE.
      GET_POP_TYPE = POP_TYPE
      RETURN

      END FUNCTION GET_POP_TYPE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_pop_xmw
!
! !DESCRIPTION: Function GET\_POP\_XMW returns POP molecular weight in kg/mol.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_POP_XMW( IN_XMW )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: IN_XMW
!
! !REVISION HISTORY:
!  30 Sep 2012 - C. Pike Thackray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8             :: GET_POP_XMW
      REAL*8, SAVE       :: POP_XMW
      LOGICAL, SAVE      :: IS_SET

      !=================================================================
      ! GET_POP_XMW begins here
      !=================================================================

      IF ( IS_SET ) THEN
         GET_POP_XMW = POP_XMW
         RETURN
      ENDIF

      POP_XMW = IN_XMW
      IS_SET = .TRUE.
      GET_POP_XMW = POP_XMW
      RETURN

      END FUNCTION GET_POP_XMW
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_pop_koa
!
! !DESCRIPTION: Function GET\_POP\_KOA returns the POP octanol-water
!  partition coefficient [unitless].
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_POP_KOA( IN_KOA )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: IN_KOA
!
! !REVISION HISTORY:
!  30 Sep 2012 - C. Pike Thackray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8             :: GET_POP_KOA
      REAL*8, SAVE       :: POP_KOA
      LOGICAL, SAVE      :: IS_SET

      !=================================================================
      ! GET_POP_KOA begins here
      !=================================================================

      IF ( IS_SET ) THEN
         GET_POP_KOA = POP_KOA
         RETURN
      ENDIF

      POP_KOA = IN_KOA
      IS_SET = .TRUE.
      GET_POP_KOA = POP_KOA
      RETURN

      END FUNCTION GET_POP_KOA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_pop_kbc
!
! !DESCRIPTION: Function GET\_POP\_KBC returns the POP black carbon-air
!  partition coefficient [unitless].
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_POP_KBC( IN_KBC )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: IN_KBC
!
! !REVISION HISTORY:
!  30 Sep 2012 - C. Pike Thackray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8             :: GET_POP_KBC
      REAL*8, SAVE       :: POP_KBC
      LOGICAL, SAVE      :: IS_SET

      !=================================================================
      ! GET_POP_KBC begins here
      !=================================================================

      IF ( IS_SET ) THEN
         GET_POP_KBC = POP_KBC
         RETURN
      ENDIF

      POP_KBC = IN_KBC
      IS_SET = .TRUE.
      GET_POP_KBC = POP_KBC
      RETURN

      END FUNCTION GET_POP_KBC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_pop_k_popg_oh
!
! !DESCRIPTION: Function GET\_POP\_K\_POPG\_OH returns the POP reaction rate
!  constant for reaction of gas phase POP with hydroxyl radical [cm3/molecule/s]
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_POP_K_POPG_OH( IN_K_POPG_OH )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: IN_K_POPG_OH
!
! !REVISION HISTORY:
!  30 Sep 2012 - C. Pike Thackray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8             :: GET_POP_K_POPG_OH
      REAL*8,  SAVE      :: POP_K_POPG_OH
      LOGICAL, SAVE      :: IS_SET

      !=================================================================
      ! GET_POP_K_POPG_OH begins here
      !=================================================================

      IF ( IS_SET ) THEN
         GET_POP_K_POPG_OH = POP_K_POPG_OH
         RETURN
      ENDIF

      POP_K_POPG_OH = IN_K_POPG_OH
      IS_SET = .TRUE.
      GET_POP_K_POPG_OH = POP_K_POPG_OH
      RETURN

      END FUNCTION GET_POP_K_POPG_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_pop_k_popp_o3a
!
! !DESCRIPTION: Function GET\_POP\_K\_POPP\_O3A returns the POP reaction rate
!  constant for reaction of particle phase POP with ozone [s\^-1].
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_POP_K_POPP_O3A( IN_K_POPP_O3A )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: IN_K_POPP_O3A
!
! !REVISION HISTORY:
!  30 Sep 2012 - C. Pike Thackray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8             :: GET_POP_K_POPP_O3A
      REAL*8, SAVE       :: POP_K_POPP_O3A
      LOGICAL, SAVE      :: IS_SET

      !=================================================================
      ! GET_POP_K_POPP_O3A begins here
      !=================================================================

      IF ( IS_SET ) THEN
         GET_POP_K_POPP_O3A = POP_K_POPP_O3A
         RETURN
      ENDIF

      POP_K_POPP_O3A = IN_K_POPP_O3A
      IS_SET = .TRUE.
      GET_POP_K_POPP_O3A = POP_K_POPP_O3A
      RETURN

      END FUNCTION GET_POP_K_POPP_O3A
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_pop_k_popp_o3b
!
! !DESCRIPTION: Function GET\_POP\_K\_POPP\_O3B returns the POP reaction rate
!  constant for reaction of particle phase POP with ozone [molec/cm3].
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_POP_K_POPP_O3B( IN_K_POPP_O3B )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: IN_K_POPP_O3B
!
! !REVISION HISTORY:
!  30 Sep 2012 - C. Pike Thackray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8             :: GET_POP_K_POPP_O3B
      REAL*8, SAVE       :: POP_K_POPP_O3B
      LOGICAL, SAVE      :: IS_SET

      !=================================================================
      ! GET_POP_K_POPP_O3B begins here
      !=================================================================

      IF ( IS_SET ) THEN
         GET_POP_K_POPP_O3B = POP_K_POPP_O3B
         RETURN
      ENDIF

      POP_K_POPP_O3B = IN_K_POPP_O3B
      IS_SET = .TRUE.
      GET_POP_K_POPP_O3B = POP_K_POPP_O3B
      RETURN

      END FUNCTION GET_POP_K_POPP_O3B
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_pop_hstar
!
! !DESCRIPTION: Function GET\_POP\_HSTAR returns the POP Henry's Law constant in
!  atm/M/K.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_POP_HSTAR( IN_HSTAR )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: IN_HSTAR
!
! !REVISION HISTORY:
!  30 Sep 2012 - C. Pike Thackray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8             :: GET_POP_HSTAR
      REAL*8, SAVE       :: POP_HSTAR
      LOGICAL, SAVE      :: IS_SET

      !=================================================================
      ! GET_POP_HSTAR begins here
      !=================================================================

      IF ( IS_SET ) THEN
         GET_POP_HSTAR = POP_HSTAR
         RETURN
      ENDIF

      POP_HSTAR = IN_HSTAR
      IS_SET = .TRUE.
      GET_POP_HSTAR = POP_HSTAR
      WRITE ( 6, '( a)') 'THIS is the first HSTAR'
      RETURN

      END FUNCTION GET_POP_HSTAR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_pop_del_H
!
! !DESCRIPTION: Function GET\_POP\_DEL\_H returns the enthalpy of air-water
!  exchange (K).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_POP_DEL_H( IN_DEL_H )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: IN_DEL_H
!
! !REVISION HISTORY:
!  30 Sep 2012 - C. Pike Thackray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8             :: GET_POP_DEL_H
      REAL*8, SAVE       :: POP_DEL_H
      LOGICAL, SAVE      :: IS_SET

      !=================================================================
      ! GET_POP_DEL_H begins here
      !=================================================================

      IF ( IS_SET ) THEN
         GET_POP_DEL_H = POP_DEL_H
         RETURN
      ENDIF

      POP_DEL_H = IN_DEL_H
      IS_SET = .TRUE.
      GET_POP_DEL_H = POP_DEL_H
      RETURN

      END FUNCTION GET_POP_DEL_H
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_pop_del_Hw
!
! !DESCRIPTION: Function GET\_POP\_DEL\_Hw returns the enthalpy of phase
!  transfer from gas phaseto particle phase.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_POP_DEL_Hw( IN_DEL_Hw )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: IN_DEL_Hw
!
! !REVISION HISTORY:
!  30 Sep 2012 - C. Pike Thackray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8             :: GET_POP_DEL_Hw
      REAL*8, SAVE       :: POP_DEL_Hw
      LOGICAL, SAVE      :: IS_SET

      !=================================================================
      ! GET_POP_DEL_Hw begins here
      !=================================================================

      IF ( IS_SET ) THEN
         GET_POP_DEL_Hw = POP_DEL_Hw
         RETURN
      ENDIF

      POP_DEL_Hw = IN_DEL_Hw
      IS_SET = .TRUE.
      GET_POP_DEL_Hw = POP_DEL_Hw
      RETURN

      END FUNCTION GET_POP_DEL_Hw
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  init_pop_params
!
! !DESCRIPTION: Subroutine INIT\_POP\_PARAMS initalizes POP parameters.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_POP_PARAMS( POP_XMW,        POP_KOA, 
     &                            POP_KBC,        POP_K_POPG_OH,
     &                            POP_K_POPP_O3A, POP_K_POPP_O3B, 
     &                            POP_HSTAR,      POP_DEL_H,
     &                            POP_DEL_Hw      )
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8 :: POP_XMW,        POP_KOA
      REAL*8 :: POP_KBC,        POP_K_POPG_OH
      REAL*8 :: POP_K_POPP_O3A, POP_K_POPP_O3B
      REAL*8 :: POP_HSTAR,      POP_DEL_H
      REAL*8 :: POP_DEL_Hw
!
! !REVISION HISTORY:
!  30 Sep 2012 - C. Pike Thackray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! INIT_POP_PARAMS begins here
      !=================================================================

      POP_XMW        = GET_POP_XMW( POP_XMW )
      POP_KOA        = GET_POP_KOA( POP_KOA )
      POP_KBC        = GET_POP_KBC( POP_KBC )
      POP_K_POPG_OH  = GET_POP_K_POPG_OH( POP_K_POPG_OH )
      POP_K_POPP_O3A = GET_POP_K_POPP_O3A( POP_K_POPP_O3A )
      POP_K_POPP_O3B = GET_POP_K_POPP_O3B( POP_K_POPP_O3B )
      POP_HSTAR      = GET_POP_HSTAR( POP_HSTAR )
      POP_DEL_H      = GET_POP_DEL_H( POP_DEL_H )
      POP_DEL_Hw     = GET_POP_DEL_Hw( POP_DEL_Hw )

      END SUBROUTINE INIT_POP_PARAMS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_emissfile
!
! !DESCRIPTION: Function GET\_EMISSFILE returns the emissions file for
!  particular POP.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_EMISSFILE( IN_FILE )
!
! !INPUT PARAMETERS: 
!
      CHARACTER, INTENT(IN)    :: IN_FILE
!
! !REVISION HISTORY:
!  30 Sep 2012 - C. Pike Thackray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=225)       :: GET_EMISSFILE
      CHARACTER(LEN=225), SAVE :: EMISSFILE
      LOGICAL, SAVE            :: IS_SET

      !=================================================================
      ! GET_EMISSFILE begins here
      !=================================================================

      IF ( IS_SET ) THEN
         GET_EMISSFILE = EMISSFILE
         RETURN
      ENDIF

      EMISSFILE = IN_FILE
      IS_SET = .TRUE.
      GET_EMISSFILE = EMISSFILE
      RETURN

      END FUNCTION GET_EMISSFILE
!EOC
      END MODULE GET_POPSINFO_MOD

