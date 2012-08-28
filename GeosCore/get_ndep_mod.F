! $Id:get_ndep_mod.F,v 1.1.1.1 2009/10/24 20:34:56 rch Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group         !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: get_ndep_mod
!     
! !DESCRIPTION: Module GET\_NDEP\_MOD contains routines for computing the
!  nitrogen dry and wet deposition. This quantity is needed by 
!  GEOS-Chem soil emissions "FERT_AW"
!\\
!\\ 
! !INTERFACE: 
!     
      MODULE GET_NDEP_MOD

      IMPLICIT NONE
      PRIVATE
  
! !PUBLIC MEMBER FUNCTIONS:
!      PUBLIC  :: GET_NDEP
      PUBLIC  :: SOIL_DRYDEP
      PUBLIC  :: SOIL_WETDEP
      PUBLIC  :: GET_DEP_N


!
! !USES
!
!      USE TRACERID_MOD, ONLY : IDTNH3, IDTNH4, IDTHNO3, IDTNIT
!      USE COMMSOIL_MOD
!      USE CMN_SIZE_MOD

!
!  !REMARKS:
!  References:
!  ===========================================================================
! !REVISION HISTORY:
!     
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: soil_drydep
!
! !DESCRIPTION: Subroutine SOIL_DRY_DEP holds dry deposited 
!               species [molec/cm2/s]. This is
!               called from dry_dep_mod.f
!\\

      !###Rynda
      SUBROUTINE SOIL_DRYDEP( I, J, L, NN, TDRYFX )

      ! Reference to F90 modules
      USE TRACERID_MOD,      ONLY : IDTNH4,    IDTNIT, IDTNH4aq 
      USE TRACERID_MOD,      ONLY : IDTHNO3,   IDTNH3, IDTNITs
      USE TRACERID_MOD,      ONLY : IDTNOx,    IDTPAN
      USE COMMSOIL_MOD

      INTEGER, INTENT(IN)  :: I           ! I
      INTEGER, INTENT(IN)  :: J           ! J
      INTEGER, INTENT(IN)  :: L           ! Level
      INTEGER, INTENT(IN)  :: NN          ! Dry Dep Tracer # 
      REAL*8 , INTENT(IN)  :: TDRYFX      ! Dry dep flux [molec/cm2/s] 


      IF ( NN == IDTNOx   ) DRY_NO2(I,J)  = DRY_NO2(I,J)  + TDRYFX 
      IF ( NN == IDTPAN   ) DRY_PAN(I,J)  = DRY_PAN(I,J)  + TDRYFX 
      IF ( NN == IDTHNO3  ) DRY_HNO3(I,J) = DRY_HNO3(I,J) + TDRYFX 
      IF ( NN == IDTNH3   ) DRY_NH3(I,J)  = DRY_NH3(I,J)  + TDRYFX 
      IF ( NN == IDTNH4   ) DRY_NH4(I,J)  = DRY_NH4(I,J)  + TDRYFX 
      IF ( NN == IDTNH4aq ) DRY_NH4(I,J)  = DRY_NH4(I,J)  + TDRYFX
      IF ( NN == IDTNIT   ) DRY_NIT(I,J)  = DRY_NIT(I,J)  + TDRYFX
      IF ( NN == IDTNITs  ) DRY_NIT(I,J)  = DRY_NIT(I,J)  + TDRYFX

      END SUBROUTINE   SOIL_DRYDEP



!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: soil_drydep
!
! !DESCRIPTION: Subroutine SOIL_WETDEP holds wet deposited 
!               species [molec/cm2/s]. This is
!               called from dry_dep_mod.f
!\\

      SUBROUTINE SOIL_WETDEP( I, J, L, NN, TWETFX )

      ! Reference to F90 modules
      USE TRACERID_MOD,      ONLY : IDTNH4,    IDTNIT, IDTNH4aq 
      USE TRACERID_MOD,      ONLY : IDTHNO3,   IDTNH3, IDTNITs
      USE COMMSOIL_MOD

      INTEGER, INTENT(IN)  :: I           ! I
      INTEGER, INTENT(IN)  :: J           ! J
      INTEGER, INTENT(IN)  :: L           ! Level
      INTEGER, INTENT(IN)  :: NN          ! Wet Dep Tracer # 
      REAL*8 , INTENT(IN)  :: TWETFX      ! Wet dep flux [kg/s] 


      IF ( NN == IDTHNO3  ) WET_HNO3(I,J) = WET_HNO3(I,J) + TWETFX 
      IF ( NN == IDTNH3   ) WET_NH3(I,J)  = WET_NH3(I,J)  + TWETFX 
      IF ( NN == IDTNH4   ) WET_NH4(I,J)  = WET_NH4(I,J)  + TWETFX 
      IF ( NN == IDTNH4aq ) WET_NH4(I,J)  = WET_NH4(I,J)  + TWETFX
      IF ( NN == IDTNIT   ) WET_NIT(I,J)  = WET_NIT(I,J)  + TWETFX
      IF ( NN == IDTNITs  ) WET_NIT(I,J)  = WET_NIT(I,J)  + TWETFX

      END SUBROUTINE   SOIL_WETDEP



!------------------------------------------------------------------------------

!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GET_DEP_N
!
! !DESCRIPTION: Subroutine sums dry and wet deposition since prev.
!               timestep and calculates contribution to fertilizer
!               N source.       
!\\

      SUBROUTINE GET_DEP_N ( I, J, TS_EMIS, DEP_FERT )

      !Just for testing
	!Rynda Testing !maasatemp
      USE COMMSOIL_MOD

      INTEGER,  INTENT(IN)  :: I
      INTEGER,  INTENT(IN)  :: J
      REAL*8 ,  INTENT(IN)  :: TS_EMIS     ! Emission timestep [min]

      ! Return value to emissdr.f Dep emitted as Fert [ng N/m2/s]
      REAL*8 ,  INTENT(INOUT) :: DEP_FERT  


      REAL*8,  PARAMETER :: TAU_MONTHS   = 6. !Decay rate of dep. N [months]
      REAL*8,  PARAMETER :: SECPERDAY    = 86400.d0
      REAL*8,  PARAMETER :: DAYSPERMONTH = 30.

      !Testing
      REAL*8  ::  LAT, LON

      REAL*8  :: DRYN  !Dry dep. N since prev timestep
                                    ! Units ng N/m2/s     
      REAL*8  :: WETN  !Wet dep. N since prev timestep 
      REAL*8  :: DEPN  !dep. N since prev timestep 

      REAL*8  :: C1
      REAL*8  :: C2 
      REAL*8  :: TAU_SEC
      REAL*8  :: TS_SEC  

      !Total all N species & convert molec/cm2/s --> ng N/m2/s
      DRYN = SOURCE_DRYN( I,  J )

      !Total all N species & convert kg/s --> ng N/m2/s
      WETN = SOURCE_WETN( I,  J )

      ! Sum wet and dry deposition
      DEPN = DRYN + WETN 


      !Emission Timestep in seconds
      TS_SEC = TS_EMIS * 60.d0

      !Do mass balance (see Intro to Atm Chem Chap. 3)
      ! m(t) = m(0) * exp(-t/tau) + Source * tau * (1 - exp(-t/tau))

      !convert months -->  seconds (assume 30 days months)
      TAU_SEC = TAU_MONTHS * DAYSPERMONTH * SECPERDAY

      C1 = EXP( - TS_SEC / TAU_SEC)
      C2 = 1.d0 - C1        

      ! ng N/m2
      DEP_RESERVOIR(I,J) = DEP_RESERVOIR (I,J) * C1
     &                        + DEPN * TAU_SEC * C2

      ! 40% runoff 
      DEP_FERT = DEP_RESERVOIR(I,J) * 0.6d0

      IF ( I == 41 .and. J == 65) THEN
         WRITE(6,*) 'DEPFERT',  DEP_FERT
         WRITE(6,*) 'RES' ,     DEP_RESERVOIR(I,J)
	 WRITE(6,*) 'DRYN (ng N /m2/s)', DRYN
         WRITE(6,*) 'WETN (ng N /m2/s)', WETN
         WRITE(6,*) 'DRYN + WETN (ng N /m2/s)', DRYN + WETN
         CALL FLUSH(6)

      ENDIF
    

      END SUBROUTINE  GET_DEP_N
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_Ndry
!
! !DESCRIPTION: This subroutine gets dry deposited Nitrogen since
!               last emission time step, converts to ng N/m2/s
!\\
!\\
! !INTERFACE:
!
      FUNCTION SOURCE_DRYN( I, J  ) RESULT( DRYN )

      USE COMMSOIL_MOD   
      USE CMN_DIAG_MOD    
      !Just for testing
      USE GRID_MOD,          ONLY : GET_YMID !Maasatemp
      USE GRID_MOD,          ONLY : GET_XMID !Maasatemp
      USE DIAG_MOD,          ONLY : AD44 !TESTING
       
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: I           
      INTEGER, INTENT(IN)  :: J           

! !RETURN VALUE:
!
      REAL*8               :: DRYN         !Dry dep. N since prev timestep
      REAL*8               :: LON, LAT , CONV !Testing

! !LOCAL VARIABLES:

      REAL*8,  PARAMETER   :: AVO = 6.022d23 !Avogadro's # (molec/mol)
      REAL*8,  PARAMETER   :: MW  = 14. ! molec. weight nitrogen (14gN/mol)
      REAL*8,  PARAMETER   :: NG_PER_GRAM  = 1.d9
      REAL*8,  PARAMETER   :: CM2_PER_M2  = 1.d4

      
      !Total
      DRYN = DRY_HNO3(I,J) + DRY_NH3(I,J) + DRY_NH4(I,J) + DRY_NIT(I,J)
     &       + DRY_NO2(I,J) + DRY_PAN(I,J)  
 
      !Molecules/cm2/s --> ngN/m2/s 
      DRYN = DRYN * CM2_PER_M2 / AVO * MW * NG_PER_GRAM 

      IF ( I == 41 .and. J == 65) THEN
      CONV = CM2_PER_M2 / AVO * MW * NG_PER_GRAM 

      !Rynda Testing !maasatemp
      LAT     = GET_YMID( I,J,1 ) !maasatemp
      LON     = GET_XMID( I,J,1 ) !maasatemp


      WRITE(6,*)   LON, LAT     
      WRITE(6,*)  'DRY NO2_SOURCE(ng)', DRY_NO2(I,J) * CONV 
      WRITE(6,*)  'DRY PAN_SOURCE(ng)', DRY_PAN(I,J) * CONV 
      WRITE(6,*)  'DRY HNO3_SOURCE(ng)', DRY_HNO3(I,J) * CONV 
      WRITE(6,*)  'DRY HNO3_SOURCE(molec/cm2/s)', DRY_HNO3(I,J) 
      !Testing
      	IF ( ND44 > 0 ) THEN
      	WRITE(6,*)  'AD44 HNO3', AD44(I,J,4,1)  - TOTAL_DRYOLD
      	ENDIF
      WRITE(6,*)  'DRY NIT_SOURCE(ng)',  DRY_NIT(I,J)  * CONV
      WRITE(6,*)  'DRY NH4_SOURCE(ng)',  DRY_NH4(I,J)  * CONV
      WRITE(6,*)  'DRY NH3_SOURCE(ng)',  DRY_NH3(I,J)  * CONV    
      CALL FLUSH(6)
      
	IF ( ND44 > 0 ) THEN
      	TOTAL_DRYOLD =  AD44(I,J,4,1)
	ENDIF 

      ENDIF


      END FUNCTION SOURCE_DRYN

 !------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_Ndry
!
! !DESCRIPTION: This subroutine gets wet deposited Nitrogen since
!               last emission time step, converts to ng N/m2/s
!\\
!\\
! !INTERFACE:
!
      FUNCTION SOURCE_WETN( I, J  ) RESULT(WETN )


      !Rynda Testing
      USE GRID_MOD,          ONLY : GET_YMID !Maasatemp
      USE GRID_MOD,          ONLY : GET_XMID !Maasatemp
      USE TIME_MOD,          ONLY : GET_TS_EMIS
      USE TIME_MOD,          ONLY : GET_TS_DYN
      USE DIAG_MOD,          ONLY : AD38, AD39 !TESTING
      USE COMMSOIL_MOD
      USE GRID_MOD,          ONLY : GET_AREA_M2
      USE CMN_DIAG_MOD

! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: I           
      INTEGER, INTENT(IN)  :: J           

! !RETURN VALUE:
!
      REAL*8               :: WETN         !Dry dep. N since prev timestep

! !LOCAL VARIABLES:
  
      REAL*8,  PARAMETER   :: CONVHNO3 = 0.222 !MWN/MWHNO3
      REAL*8,  PARAMETER   :: CONVNH4  = 0.777 !MWN/MWNH4
      REAL*8,  PARAMETER   :: CONVNH3  = 0.823 !MWN/MWNH3
      REAL*8,  PARAMETER   :: CONVNIT  = 0.226 !MWN/MWNIT
      REAL*8               :: LON, LAT !Testing
      REAL*8               :: HNO3NG, NITNG, NH3NG, NH4NG
      REAL*8               :: NTS, TS_EMIS, TS_DYN, AREA_M2     
      REAL*8               :: TOTAL_HNO3 !Testing
      INTEGER              :: L !Testing


 
      !divide through by number of transport timesteps
      !because WET_xxxx is summed over transport timesteps
      ! need to get average

      TS_EMIS = GET_TS_EMIS()
      TS_DYN  = GET_TS_DYN() 
      NTS     = TS_EMIS/TS_DYN 
      AREA_M2 = GET_AREA_M2( 1, J, 1 )

      !WRITE(6,*) 'NTS', NTS
      !In my code and I think as it should be in everyones 
      ! EMIS TS = 2 x DYN TS 
      IF (NTS /= 2.) WRITE(6,*) 'ERROR! CHECK THIS'

      !convert kg species/s * xx dyn t.s --> ng N/m2/s
      HNO3NG  = WET_HNO3(I,J) * CONVHNO3 * 1.d12 / 
     &          AREA_M2 / NTS
      NH4NG   = WET_NH4(I,J)  * CONVNH4  * 1.d12 / 
     &          AREA_M2 / NTS  
      NH3NG   = WET_NH3(I,J)  * CONVNH3  * 1.d12 / 
     &          AREA_M2 / NTS 
      NITNG   = WET_NIT(I,J)  * CONVNIT  * 1.d12 / 
     &          AREA_M2 / NTS  

      !Total N wet dep
      WETN = HNO3NG + NH3NG  + NH4NG + NITNG  

	IF ( ND38 > 0 .and. ND39 > 0 ) THEN
      		IF ( I == 41 .and. J == 65) THEN
	!Rynda Testing Maasa
      LAT     = GET_YMID( I,J,1 ) !maasatemp
      LON     = GET_XMID( I,J,1 ) !maasatemp

      WRITE(6,*)   LON, LAT
      WRITE(6,*)  '# DYN TS in getdep', NTS
      WRITE(6,*)  'WET HNO3_SOURCE(ng)', HNO3NG
      WRITE(6,*)  'WET HNO3_SOURCE(kg/s)', WET_HNO3(I,J)  
      WRITE(6,*)  'WET NIT_SOURCE(ng)',  NITNG
      WRITE(6,*)  'WET NH4_SOURCE(ng)',  NH4NG 
      WRITE(6,*)  'WET NH3_SOURCE(ng)',  NH3NG
      CALL FLUSH(6)

	
      TOTAL_HNO3 = 0.
      DO L = 1,30          
         TOTAL_HNO3 = TOTAL_HNO3 + AD39(I,J,L,1) + AD38(I,J,L,1)
      ENDDO

      WRITE(6,*)  'WET DIAG HNO3(kg/s)', TOTAL_HNO3 - TOTAL_WETOLD

      TOTAL_WETOLD = TOTAL_HNO3

		ENDIF
      	ENDIF


      END FUNCTION SOURCE_WETN

 

      END MODULE GET_NDEP_MOD


