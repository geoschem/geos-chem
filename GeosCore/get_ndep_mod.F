!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group         !
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
      PUBLIC  :: RESET_DEP_N
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

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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

      ! ckeller (14/04/02)
      IF ( NN == IDTNO2 .OR. NN == IDTPAN .OR. NN == IDTHNO3  .OR.
     &     NN == IDTNH3 .OR. NN == IDTNH4 .OR. NN == IDTNH4aq .OR. 
     &     NN == IDTNIT .OR. NN == IDTNITs                        ) THEN
         DRY_TOTN(I,J) = DRY_TOTN(I,J) + TDRYFX
      ENDIF

      END SUBROUTINE   SOIL_DRYDEP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      ! ckeller (14/04/02)
      REAL*8        :: SCAL

      ! ckeller (14/04/02)
      SCAL = 0d0
      IF ( NN == IDTHNO3  ) SCAL = CONVHNO3 
      IF ( NN == IDTNH3   ) SCAL = CONVNH3
      IF ( NN == IDTNH4   ) SCAL = CONVNH4
      IF ( NN == IDTNH4aq ) SCAL = CONVNH4
      IF ( NN == IDTNIT   ) SCAL = CONVNIT
      IF ( NN == IDTNITs  ) SCAL = CONVNIT

      IF ( SCAL > 0d0 ) THEN
         WET_TOTN(I,J) = WET_TOTN(I,J) + ( TWETFX * SCAL )
      ENDIF

      END SUBROUTINE SOIL_WETDEP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TIME_MOD,     ONLY : GET_TS_CHEM
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
      REAL*8,  PARAMETER   :: NG_P_GRAM = 1.d9
      REAL*8,  PARAMETER   :: CM2_P_M2  = 1.d4
      REAL*8               :: NTS, TS_EMIS, TS_CHEM

      ! ckeller (14/04/02)
      TS_EMIS = GET_TS_EMIS()
      TS_CHEM = GET_TS_CHEM() 
      NTS     = TS_EMIS/TS_CHEM
      DRYN    = DRY_TOTN(I,J) / NTS * CM2_P_M2 / AVO * MW * NG_P_GRAM

      END FUNCTION SOURCE_DRYN

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      IF (NTS /= 2.) WRITE(6,*) 'ERROR! CHECK THIS'

      ! ckeller (14/04/02)
      WETN = WET_TOTN(I,J) / NTS * 1.d12 / AREA_M2

      END FUNCTION SOURCE_WETN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: reset_dep_N
!
! !DESCRIPTION: Subroutine RESET\_DEP\_N resets the dry and wet deposition 
!               arrays and variables so that they can be refilled. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RESET_DEP_N 
!
! !USES:
!
      USE COMMSOIL_MOD
!
! !REVISION HISTORY:
!  03 Apr 2014 - C. Keller      - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Reset all variables
      DRY_TOTN = 0d0
      WET_TOTN = 0d0

      END SUBROUTINE  RESET_DEP_N
!EOC
      END MODULE GET_NDEP_MOD
