!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: get_ndep_mod
!     
! !DESCRIPTION: Module GET\_NDEP\_MOD contains routines for computing the
!  nitrogen dry and wet deposition. This quantity is needed by 
!  GEOS-Chem soil emissions "FERT\_AW"
!\\
!\\ 
! !INTERFACE: 
!     
      MODULE GET_NDEP_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
! 
! !PUBLIC MEMBER FUNCTIONS:
!      
      PUBLIC  :: SOIL_DRYDEP
      PUBLIC  :: SOIL_WETDEP
      PUBLIC  :: GET_DEP_N
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers    
!EOP
!------------------------------------------------------------------------------
!BOC
!
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: soil_drydep
!
! !DESCRIPTION: Subroutine SOIL\_DRY\_DEP holds dry deposited species
!               [molec/cm2/s]. This is called from dry\_dep\_mod.F.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SOIL_DRYDEP( I, J, L, NN, TDRYFX )
!
! !USES:
!
      USE TRACERID_MOD, ONLY : IDTNH4,    IDTNIT, IDTNH4aq 
      USE TRACERID_MOD, ONLY : IDTHNO3,   IDTNH3, IDTNITs
      USE TRACERID_MOD, ONLY : IDTNO2,     IDTPAN
      USE COMMSOIL_MOD
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I          ! I
      INTEGER, INTENT(IN)   :: J          ! J
      INTEGER, INTENT(IN)   :: L          ! Level
      INTEGER, INTENT(IN)   :: NN         ! Dry Dep Tracer # 
      REAL*8 , INTENT(IN)   :: TDRYFX     ! Dry dep flux [molec/cm2/s] 
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!  14 Mar 2013 - M. Payer    - Replace NOx with NO2 as part of removal of
!                              NOx-Ox partitioning
!EOP
!------------------------------------------------------------------------------
!BOC
!
      IF ( NN == IDTNO2   ) DRY_NO2(I,J)  = DRY_NO2(I,J)  + TDRYFX 
      IF ( NN == IDTPAN   ) DRY_PAN(I,J)  = DRY_PAN(I,J)  + TDRYFX 
      IF ( NN == IDTHNO3  ) DRY_HNO3(I,J) = DRY_HNO3(I,J) + TDRYFX 
      IF ( NN == IDTNH3   ) DRY_NH3(I,J)  = DRY_NH3(I,J)  + TDRYFX 
      IF ( NN == IDTNH4   ) DRY_NH4(I,J)  = DRY_NH4(I,J)  + TDRYFX 
      IF ( NN == IDTNH4aq ) DRY_NH4(I,J)  = DRY_NH4(I,J)  + TDRYFX
      IF ( NN == IDTNIT   ) DRY_NIT(I,J)  = DRY_NIT(I,J)  + TDRYFX
      IF ( NN == IDTNITs  ) DRY_NIT(I,J)  = DRY_NIT(I,J)  + TDRYFX

      END SUBROUTINE   SOIL_DRYDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: soil_wetdep
!
! !DESCRIPTION: Subroutine SOIL\_WETDEP holds wet deposited species
!               [molec/cm2/s]. This is called from wetscav\_mod.F.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SOIL_WETDEP( I, J, L, NN, TWETFX )
!
! !USES:
!
      USE TRACERID_MOD, ONLY : IDTNH4,    IDTNIT, IDTNH4aq 
      USE TRACERID_MOD, ONLY : IDTHNO3,   IDTNH3, IDTNITs
      USE COMMSOIL_MOD
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I          ! I
      INTEGER, INTENT(IN)   :: J          ! J
      INTEGER, INTENT(IN)   :: L          ! Level
      INTEGER, INTENT(IN)   :: NN         ! Wet Dep Tracer # 
      REAL*8 , INTENT(IN)   :: TWETFX     ! Wet dep flux [kg/s] 
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
      IF ( NN == IDTHNO3  ) WET_HNO3(I,J) = WET_HNO3(I,J) + TWETFX 
      IF ( NN == IDTNH3   ) WET_NH3(I,J)  = WET_NH3(I,J)  + TWETFX 
      IF ( NN == IDTNH4   ) WET_NH4(I,J)  = WET_NH4(I,J)  + TWETFX 
      IF ( NN == IDTNH4aq ) WET_NH4(I,J)  = WET_NH4(I,J)  + TWETFX
      IF ( NN == IDTNIT   ) WET_NIT(I,J)  = WET_NIT(I,J)  + TWETFX
      IF ( NN == IDTNITs  ) WET_NIT(I,J)  = WET_NIT(I,J)  + TWETFX

      END SUBROUTINE   SOIL_WETDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_dep_N
!
! !DESCRIPTION: Subroutine GET\_DEP\_N sums dry and wet deposition since prev.
!               timestep and calculates contribution to fertilizer N source.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_DEP_N ( I, J, TS_EMIS, DEP_FERT )
!
! !USES:
!
      USE COMMSOIL_MOD
!
! !INPUT PARAMETERS: 
!
      INTEGER,  INTENT(IN)  :: I
      INTEGER,  INTENT(IN)  :: J
      REAL*8 ,  INTENT(IN)  :: TS_EMIS     ! Emission timestep [min]
!
! !INPUT/OUTPUT PARAMETERS: 
!
      ! Dep emitted as Fert [ng N/m2/s]
      REAL*8 ,  INTENT(INOUT) :: DEP_FERT  
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER :: TAU_MONTHS   = 6. ! Decay rate of dep. N [months]
      REAL*8,  PARAMETER :: SECPERDAY    = 86400.d0
      REAL*8,  PARAMETER :: DAYSPERMONTH = 30.
!
! !LOCAL VARIABLES:
!
      REAL*8             :: DRYN  ! Dry dep. N since prev timestep
                                  ! Units ng N/m2/s     
      REAL*8             :: WETN  ! Wet dep. N since prev timestep 
      REAL*8             :: DEPN  ! dep. N since prev timestep 

      REAL*8             :: C1
      REAL*8             :: C2 
      REAL*8             :: TAU_SEC
      REAL*8             :: TS_SEC  

      !Total all N species & convert molec/cm2/s --> ng N/m2/s
      DRYN = SOURCE_DRYN( I,  J )

      !Total all N species & convert kg/s --> ng N/m2/s
      WETN = SOURCE_WETN( I,  J )

      ! Sum wet and dry deposition
      DEPN = DRYN + WETN 

      !Emission Timestep in seconds
      TS_SEC = TS_EMIS * 60.d0

      !Do mass balance (see Intro to Atm Chem Chap. 3)
      !m(t) = m(0) * exp(-t/tau) + Source * tau * (1 - exp(-t/tau))

      !convert months -->  seconds (assume 30 days months)
      TAU_SEC = TAU_MONTHS * DAYSPERMONTH * SECPERDAY

      C1 = EXP( - TS_SEC / TAU_SEC)
      C2 = 1.d0 - C1        

      ! ng N/m2
      DEP_RESERVOIR(I,J) = DEP_RESERVOIR (I,J) * C1
     &                        + DEPN * TAU_SEC * C2

      ! 40% runoff 
      DEP_FERT = DEP_RESERVOIR(I,J) * 0.6d0   

      END SUBROUTINE  GET_DEP_N
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: source_dryN
!
! !DESCRIPTION: Subroutine SOURCE\_DRYN gets dry deposited Nitrogen since
!               last emission time step, converts to ng N/m2/s.
!\\
!\\
! !INTERFACE:
!
      FUNCTION SOURCE_DRYN( I, J  ) RESULT( DRYN )
!
! !USES:
!
      USE COMMSOIL_MOD   
      USE CMN_DIAG_MOD        
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: I           
      INTEGER, INTENT(IN)  :: J           
!
! !RETURN VALUE:
!
      REAL*8               :: DRYN         !Dry dep. N since prev timestep
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8,  PARAMETER   :: AVO         = 6.022d23 ! Avogadro's # (molec/mol)
      REAL*8,  PARAMETER   :: MW          = 14.d0    ! Molec weight N (14gN/mol)
      REAL*8,  PARAMETER   :: NG_PER_GRAM = 1.d9
      REAL*8,  PARAMETER   :: CM2_PER_M2  = 1.d4
      
      !Total
      DRYN = DRY_HNO3(I,J) + DRY_NH3(I,J) + DRY_NH4(I,J) + DRY_NIT(I,J)
     &     + DRY_NO2(I,J)  + DRY_PAN(I,J)  
 
      !Molecules/cm2/s --> ngN/m2/s 
      DRYN = DRYN * CM2_PER_M2 / AVO * MW * NG_PER_GRAM 

      END FUNCTION SOURCE_DRYN

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: source_wetN
!
! !DESCRIPTION: Subroutine SOURCE\_WETN gets wet deposited Nitrogen since
!               last emission time step, converts to ng N/m2/s.
!\\
!\\
! !INTERFACE:
!
      FUNCTION SOURCE_WETN( I, J  ) RESULT(WETN )
!
! !USES:
!
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TIME_MOD,     ONLY : GET_TS_DYN
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE COMMSOIL_MOD
      USE CMN_DIAG_MOD
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: I           
      INTEGER, INTENT(IN)  :: J           
!
! !RETURN VALUE:
!
      REAL*8               :: WETN         !Dry dep. N since prev timestep
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!  
      REAL*8,  PARAMETER   :: CONVHNO3 = 0.222 ! MWN/MWHNO3
      REAL*8,  PARAMETER   :: CONVNH4  = 0.777 ! MWN/MWNH4
      REAL*8,  PARAMETER   :: CONVNH3  = 0.823 ! MWN/MWNH3
      REAL*8,  PARAMETER   :: CONVNIT  = 0.226 ! MWN/MWNIT
!
! !LOCAL VARIABLES:
!
      REAL*8               :: HNO3NG, NITNG,   NH3NG,  NH4NG
      REAL*8               :: NTS,    TS_EMIS, TS_DYN, AREA_M2     
 
      ! Divide through by number of transport timesteps
      ! because WET_xxxx is summed over transport timesteps
      ! need to get average

      TS_EMIS = GET_TS_EMIS()
      TS_DYN  = GET_TS_DYN() 
      NTS     = TS_EMIS/TS_DYN 
      AREA_M2 = GET_AREA_M2( 1, J, 1 )

      ! In my code and I think as it should be in everyones 
      ! EMIS TS = 2 x DYN TS 
      ! Comment out for test simulations where transport is not done
      ! (skim, 1/14/13)
      ! IF (NTS /= 2.) WRITE(6,*) 'ERROR! CHECK THIS'

      ! Convert kg species/s * xx dyn t.s --> ng N/m2/s
      HNO3NG  = WET_HNO3(I,J) * CONVHNO3 * 1.d12 / 
     &          AREA_M2 / NTS
      NH4NG   = WET_NH4(I,J)  * CONVNH4  * 1.d12 / 
     &          AREA_M2 / NTS  
      NH3NG   = WET_NH3(I,J)  * CONVNH3  * 1.d12 / 
     &          AREA_M2 / NTS 
      NITNG   = WET_NIT(I,J)  * CONVNIT  * 1.d12 / 
     &          AREA_M2 / NTS  

      ! Total N wet dep
      WETN = HNO3NG + NH3NG  + NH4NG + NITNG  

      END FUNCTION SOURCE_WETN
!EOC
      END MODULE GET_NDEP_MOD
