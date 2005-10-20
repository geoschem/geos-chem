! $Id: epa_nei_mod.f,v 1.4 2005/10/20 14:03:25 bmy Exp $
      MODULE EPA_NEI_MOD
!
!******************************************************************************
!  Module EPA_NEI_MOD contains variables and routines to read the
!  weekday/weekend emissions from the EPA/NEI emissions inventory.
!  (rch, bmy, 11/10/04, 10/3/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) USA_MASK        : Array used to mask out the continental USA
!  (2 ) EPA_WD_AN_NOX   : Avg weekday NOX  anthro  emissions [  molec/cm2/s]
!  (3 ) EPA_WD_AN_CO    : Avg weekday CO   anthro  emissions [  molec/cm2/s]
!  (4 ) EPA_WD_AN_ALK4  : Avg weekday ALK4 anthro  emissions [atoms C/cm2/s]
!  (5 ) EPA_WD_AN_ACET  : Avg weekday ACET anthro  emissions [atoms C/cm2/s]
!  (6 ) EPA_WD_AN_MEK   : Avg weekday MEK  anthro  emissions [atoms C/cm2/s]
!  (7 ) EPA_WD_AN_PRPE  : Avg weekday PRPE anthro  emissions [atoms C/cm2/s]
!  (8 ) EPA_WD_AN_C2H6  : Avg weekday C2H6 anthro  emissions [atoms C/cm2/s]
!  (9 ) EPA_WD_AN_C3H8  : Avg weekday C3H8 anthro  emissions [atoms C/cm2/s] 
!  (10) EPA_WD_AN_CH2O  : Avg weekday CH2O anthro  emissions [  molec/cm2/s]
!  (11) EPA_WD_AN_NH3   : Avg weekday NH3  anthro  emissions [  molec/cm2/s]
!  (12) EPA_WD_AN_SO2   : Avg weekday SO2  anthro  emissions [  molec/cm2/s]   
!  (13) EPA_WD_AN_SO4   : Avg weekday SO4  anthro  emissions [  molec/cm2/s]   
!  (14) EPA_WE_AN_NOX   : Avg weekend NOX  anthro  emissions [  molec/cm2/s]
!  (15) EPA_WE_AN_CO    : Avg weekend CO   anthro  emissions [atoms C/cm2/s]
!  (16) EPA_WE_AN_ALK4  : Avg weekend ALK4 anthro  emissions [atoms C/cm2/s]
!  (17) EPA_WE_AN_ACET  : Avg weekend ACET anthro  emissions [atoms C/cm2/s]
!  (18) EPA_WE_AN_MEK   : Avg weekend MEK  anthro  emissions [atoms C/cm2/s] 
!  (19) EPA_WE_AN_PRPE  : Avg weekend PRPE anthro  emissions [atoms C/cm2/s]
!  (20) EPA_WE_AN_C2H6  : Avg weekend C2H6 anthro  emissions [atoms C/cm2/s]
!  (21) EPA_WE_AN_C3H8  : Avg weekend C3H8 anthro  emissions [atoms C/cm2/s]
!  (22) EPA_WE_AN_CH2O  : Avg weekend CH2O anthro  emissions [  molec/cm2/s]
!  (23) EPA_WE_AN_NH3   : Avg weekend NH3  anthro  emissions [  molec/cm2/s]
!  (24) EPA_WE_AN_SO2   : Avg weekend SO2  anthro  emissions [  molec/cm2/s]
!  (25) EPA_WE_AN_SO4   : Avg weekend SO4  anthro  emissions [  molec/cm2/s]
!  (26) EPA_WD_BF_NOX   : Avg weekday NOX  anthro  emissions [  molec/cm2/s]
!  (27) EPA_WD_BF_CO    : Avg weekday CO   anthro  emissions [  molec/cm2/s]
!  (28) EPA_WD_BF_ALK4  : Avg weekday ALK4 anthro  emissions [atoms C/cm2/s]
!  (29) EPA_WD_BF_ACET  : Avg weekday ACET anthro  emissions [atoms C/cm2/s]
!  (30) EPA_WD_BF_MEK   : Avg weekday MEK  anthro  emissions [atoms C/cm2/s]
!  (31) EPA_WD_BF_PRPE  : Avg weekday PRPE anthro  emissions [atoms C/cm2/s]
!  (32) EPA_WD_BF_C2H6  : Avg weekday C2H6 anthro  emissions [atoms C/cm2/s]
!  (33) EPA_WD_BF_C3H8  : Avg weekday C3H8 anthro  emissions [atoms C/cm2/s] 
!  (34) EPA_WD_BF_CH2O  : Avg weekday CH2O anthro  emissions [  molec/cm2/s]
!  (35) EPA_WD_BF_NH3   : Avg weekday NH3  anthro  emissions [  molec/cm2/s]
!  (36) EPA_WD_BF_SO2   : Avg weekday SO2  anthro  emissions [  molec/cm2/s]   
!  (37) EPA_WD_BF_SO4   : Avg weekday SO4  anthro  emissions [  molec/cm2/s]   
!  (38) EPA_WE_BF_NOX   : Avg weekend NOX  anthro  emissions [  molec/cm2/s]
!  (39) EPA_WE_BF_CO    : Avg weekend CO   anthro  emissions [atoms C/cm2/s]
!  (40) EPA_WE_BF_ALK4  : Avg weekend ALK4 anthro  emissions [atoms C/cm2/s]
!  (41) EPA_WE_BF_ACET  : Avg weekend ACET anthro  emissions [atoms C/cm2/s]
!  (42) EPA_WE_BF_MEK   : Avg weekend MEK  anthro  emissions [atoms C/cm2/s] 
!  (43) EPA_WE_BF_PRPE  : Avg weekend PRPE anthro  emissions [atoms C/cm2/s]
!  (44) EPA_WE_BF_C2H6  : Avg weekend C2H6 anthro  emissions [atoms C/cm2/s]
!  (45) EPA_WE_BF_C3H8  : Avg weekend C3H8 anthro  emissions [atoms C/cm2/s]
!  (46) EPA_WE_BF_CH2O  : Avg weekend CH2O anthro  emissions [  molec/cm2/s]
!  (47) EPA_WE_BF_NH3   : Avg weekend NH3  anthro  emissions [  molec/cm2/s]
!  (48) EPA_WE_BF_SO2   : Avg weekend SO2  anthro  emissions [  molec/cm2/s]
!  (49) EPA_WE_BF_SO4   : Avg weekend SO4  anthro  emissions [  molec/cm2/s]
!
!  Module Routines:
!  ============================================================================
!  (1 ) EMISS_EPA       : Driver routine for EPA emissions
!  (2 ) READ_EPA        : Reads EPA emissions from disk
!  (3 ) READ_USA_MASK   : Reads USA Mask from disk
!  (4 ) TOTAL_ANTHRO_TG : Prints monthly anthro  emission sums in Tg or Tg C
!  (5 ) TOTAL_BIOFUEL_TG: Prints monthly biofuel emission sums in Tg or Tg C
!  (6 ) GET_USA_MASK    : Returns value of USA mask  at a (I,J) location
!  (7 ) GET_EPA_ANTHRO  : Gets EPA anthro  emissions at a (I,J) location
!  (8 ) GET_EPA_BIOFUEL : Gets EPA biofuel emissions at a (I,J) location 
!  (9 ) INIT_EPA_NEI    : Allocates and zeroes module arrays 
!  (10) CLEANUP_EPA_NEI : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by epa_nei_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) directory_mod.f : Module containing GEOS-CHEM met field and data dirs
!  (3 ) error_mod.f     : Module containing I/O error and NaN check routines
!  (3 ) file_mod.f      : Module containing file unit numbers and error checks
!  (4 ) logical_mod.f   : 
!  (5 ) time_mod.f      : Module containing routines for computing time & date
!  (6 ) transfer_mod.f  : Module containing routines to cast & resize arrays
!
!  NOTES:
!  (1 ) Prevent out of bounds errors in routines TOTAL_ANTHRO_TG and 
!        TOTAL_BIOFUEL_TG (bmy, 1/26/05)
!  (2 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "epa_nei_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CLEANUP_EPA_NEI
      PUBLIC :: EMISS_EPA_NEI
      PUBLIC :: GET_USA_MASK
      PUBLIC :: GET_EPA_ANTHRO
      PUBLIC :: GET_EPA_BIOFUEL 

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! USA Mask
      REAL*8,  ALLOCATABLE :: USA_MASK(:,:)

      ! Fossil fuel arrays -- weekday
      REAL*4,  ALLOCATABLE :: EPA_WD_AN_NOX(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_AN_CO(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_AN_ALK4(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_AN_ACET(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_AN_MEK(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_AN_PRPE(:,:)      
      REAL*4,  ALLOCATABLE :: EPA_WD_AN_C2H6(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_AN_C3H8(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_AN_CH2O(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_AN_NH3(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_AN_SO2(:,:)      
      REAL*4,  ALLOCATABLE :: EPA_WD_AN_SO4(:,:)      
      
      ! Fossil fuel arrays -- weekend
      REAL*4,  ALLOCATABLE :: EPA_WE_AN_NOX(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_AN_CO(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_AN_ALK4(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_AN_ACET(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_AN_MEK(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_AN_PRPE(:,:)      
      REAL*4,  ALLOCATABLE :: EPA_WE_AN_C2H6(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_AN_C3H8(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_AN_CH2O(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_AN_NH3(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_AN_SO2(:,:) 
      REAL*4,  ALLOCATABLE :: EPA_WE_AN_SO4(:,:) 

      ! Biofuel arrays -- weekday
      REAL*4,  ALLOCATABLE :: EPA_WD_BF_NOX(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_BF_CO(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_BF_ALK4(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_BF_ACET(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_BF_MEK(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_BF_PRPE(:,:)      
      REAL*4,  ALLOCATABLE :: EPA_WD_BF_C2H6(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_BF_C3H8(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_BF_CH2O(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_BF_NH3(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WD_BF_SO2(:,:)  
      REAL*4,  ALLOCATABLE :: EPA_WD_BF_SO4(:,:)  

      ! Biofuel arrays -- weekend
      REAL*4,  ALLOCATABLE :: EPA_WE_BF_NOX(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_BF_CO(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_BF_ALK4(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_BF_ACET(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_BF_MEK(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_BF_PRPE(:,:)      
      REAL*4,  ALLOCATABLE :: EPA_WE_BF_C2H6(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_BF_C3H8(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_BF_CH2O(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_BF_NH3(:,:)
      REAL*4,  ALLOCATABLE :: EPA_WE_BF_SO2(:,:)  
      REAL*4,  ALLOCATABLE :: EPA_WE_BF_SO4(:,:)  

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
      
      SUBROUTINE EMISS_EPA_NEI
!
!******************************************************************************
!  Subroutine EMISS_EPA_NEI reads all EPA emissions from disk at the start
!  of a new month. (rch, bmy, 11/10/04, 8/16/05)
!
!  NOTES:
!  (1 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!******************************************************************************
!                               
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TIME_MOD,      ONLY : EXPAND_DATE,     GET_MONTH

#     include "CMN_SIZE"   ! Size parameters

      ! Local variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I, J, THISMONTH, YYYYMMDD
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! EMISS_EPA_NEI begins here!
      !=================================================================
      
      ! First-time initialization
      IF ( FIRST ) THEN

         ! Allocate arrays
         CALL INIT_EPA_NEI

         ! Read mask over the USA
         CALL READ_USA_MASK

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Get the current month
      THISMONTH = GET_MONTH()

      ! Get date for 1999 emissions
      YYYYMMDD = 19990000 + ( THISMONTH * 100 ) + 01
      
      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100   ) 
 100  FORMAT( 'E P A  /  N E I   E M I S S I O N S',
     &       '  -- Baseline Year: 1999', / )
    
      !=================================================================
      ! Read EPA weekday average anthropogenic emissions
      !=================================================================

      ! Weekday anthro file name
      FILENAME = TRIM( DATA_DIR )                         // 
     &           'EPA_NEI_200411/wkday_avg_an.YYYYMM.'    //
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()     

      ! Replace date in filename
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

      ! Read weekday data 
      CALL READ_EPA( FILENAME,       
     &               EPA_WD_AN_NOX,  EPA_WD_AN_CO,   EPA_WD_AN_ALK4, 
     &               EPA_WD_AN_ACET, EPA_WD_AN_MEK,  EPA_WD_AN_PRPE, 
     &               EPA_WD_AN_C3H8, EPA_WD_AN_CH2O, EPA_WD_AN_C2H6, 
     &               EPA_WD_AN_SO2,  EPA_WD_AN_SO4,  EPA_WD_AN_NH3 )

      !=================================================================
      ! Read EPA weekend average anthropogenic emissions
      !=================================================================

      ! Weekend anthro file name
      FILENAME = TRIM( DATA_DIR )                         // 
     &           'EPA_NEI_200411/wkend_avg_an.YYYYMM.'    //
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

      ! Replace date in filename
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

      ! Read weekday data 
      CALL READ_EPA( FILENAME,       
     &               EPA_WE_AN_NOX,  EPA_WE_AN_CO,   EPA_WE_AN_ALK4, 
     &               EPA_WE_AN_ACET, EPA_WE_AN_MEK,  EPA_WE_AN_PRPE, 
     &               EPA_WE_AN_C3H8, EPA_WE_AN_CH2O, EPA_WE_AN_C2H6, 
     &               EPA_WE_AN_SO2,  EPA_WE_AN_SO4,  EPA_WE_AN_NH3 )

      !=================================================================
      ! Read EPA weekday average biofuel emissions
      !=================================================================

      ! Weekday biofuel file name
      FILENAME = TRIM( DATA_DIR )                         // 
     &           'EPA_NEI_200411/wkday_avg_bf.YYYYMM.'    //
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

      ! Replace date in filename
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

      ! Read weekday data 
      CALL READ_EPA( FILENAME,       
     &               EPA_WD_BF_NOX,  EPA_WD_BF_CO,   EPA_WD_BF_ALK4, 
     &               EPA_WD_BF_ACET, EPA_WD_BF_MEK,  EPA_WD_BF_PRPE, 
     &               EPA_WD_BF_C3H8, EPA_WD_BF_CH2O, EPA_WD_BF_C2H6, 
     &               EPA_WD_BF_SO2,  EPA_WD_BF_SO4,  EPA_WD_BF_NH3 )

      !=================================================================
      ! Read EPA weekend average biofuel emissions
      !=================================================================

      ! Weekend biofuel file name
      FILENAME = TRIM( DATA_DIR )                         // 
     &           'EPA_NEI_200411/wkend_avg_bf.YYYYMM.'    //
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT()  

      ! Replace date in filename
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

      ! Read weekend data 
      CALL READ_EPA( FILENAME,       
     &               EPA_WE_BF_NOX,  EPA_WE_BF_CO,   EPA_WE_BF_ALK4, 
     &               EPA_WE_BF_ACET, EPA_WE_BF_MEK,  EPA_WE_BF_PRPE, 
     &               EPA_WE_BF_C3H8, EPA_WE_BF_CH2O, EPA_WE_BF_C2H6, 
     &               EPA_WE_BF_SO2,  EPA_WE_BF_SO4,  EPA_WE_BF_NH3 )

      !=================================================================
      ! Apply USA Mask (keep emissions over US, zero elsewhere)
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Weekday avg anthro
         EPA_WD_AN_NOX (I,J) = EPA_WD_AN_NOX (I,J) * USA_MASK(I,J)
         EPA_WD_AN_CO  (I,J) = EPA_WD_AN_CO  (I,J) * USA_MASK(I,J) 
         EPA_WD_AN_ALK4(I,J) = EPA_WD_AN_ALK4(I,J) * USA_MASK(I,J) 
         EPA_WD_AN_ACET(I,J) = EPA_WD_AN_ACET(I,J) * USA_MASK(I,J) 
         EPA_WD_AN_MEK (I,J) = EPA_WD_AN_MEK (I,J) * USA_MASK(I,J) 
         EPA_WD_AN_PRPE(I,J) = EPA_WD_AN_PRPE(I,J) * USA_MASK(I,J)      
         EPA_WD_AN_C3H8(I,J) = EPA_WD_AN_C3H8(I,J) * USA_MASK(I,J) 
         EPA_WD_AN_CH2O(I,J) = EPA_WD_AN_CH2O(I,J) * USA_MASK(I,J) 
         EPA_WD_AN_C2H6(I,J) = EPA_WD_AN_C2H6(I,J) * USA_MASK(I,J)
         EPA_WD_AN_SO2 (I,J) = EPA_WD_AN_SO2 (I,J) * USA_MASK(I,J)   
         EPA_WD_AN_SO4 (I,J) = EPA_WD_AN_SO4 (I,J) * USA_MASK(I,J)  
         EPA_WD_AN_NH3 (I,J) = EPA_WD_AN_NH3 (I,J) * USA_MASK(I,J) 
                             
         ! Weekend avg anthro     
         EPA_WE_AN_NOX (I,J) = EPA_WE_AN_NOX (I,J) * USA_MASK(I,J)
         EPA_WE_AN_CO  (I,J) = EPA_WE_AN_CO  (I,J) * USA_MASK(I,J)
         EPA_WE_AN_ALK4(I,J) = EPA_WE_AN_ALK4(I,J) * USA_MASK(I,J)
         EPA_WE_AN_ACET(I,J) = EPA_WE_AN_ACET(I,J) * USA_MASK(I,J)
         EPA_WE_AN_MEK (I,J) = EPA_WE_AN_MEK (I,J) * USA_MASK(I,J)
         EPA_WE_AN_PRPE(I,J) = EPA_WE_AN_PRPE(I,J) * USA_MASK(I,J)   
         EPA_WE_AN_C3H8(I,J) = EPA_WE_AN_C3H8(I,J) * USA_MASK(I,J)
         EPA_WE_AN_CH2O(I,J) = EPA_WE_AN_CH2O(I,J) * USA_MASK(I,J)
         EPA_WE_AN_C2H6(I,J) = EPA_WE_AN_C2H6(I,J) * USA_MASK(I,J)
         EPA_WE_AN_SO2 (I,J) = EPA_WE_AN_SO2 (I,J) * USA_MASK(I,J)
         EPA_WE_AN_SO4 (I,J) = EPA_WE_AN_SO4 (I,J) * USA_MASK(I,J)
         EPA_WE_AN_NH3 (I,J) = EPA_WE_AN_NH3 (I,J) * USA_MASK(I,J)
                              
         ! Weekday avg biofuel    
         EPA_WD_BF_NOX (I,J) = EPA_WD_BF_NOX (I,J) * USA_MASK(I,J)
         EPA_WD_BF_CO  (I,J) = EPA_WD_BF_CO  (I,J) * USA_MASK(I,J)
         EPA_WD_BF_ALK4(I,J) = EPA_WD_BF_ALK4(I,J) * USA_MASK(I,J)
         EPA_WD_BF_ACET(I,J) = EPA_WD_BF_ACET(I,J) * USA_MASK(I,J)
         EPA_WD_BF_MEK (I,J) = EPA_WD_BF_MEK (I,J) * USA_MASK(I,J)
         EPA_WD_BF_PRPE(I,J) = EPA_WD_BF_PRPE(I,J) * USA_MASK(I,J)
         EPA_WD_BF_C3H8(I,J) = EPA_WD_BF_C3H8(I,J) * USA_MASK(I,J)
         EPA_WD_BF_CH2O(I,J) = EPA_WD_BF_CH2O(I,J) * USA_MASK(I,J)    
         EPA_WD_BF_C2H6(I,J) = EPA_WD_BF_C2H6(I,J) * USA_MASK(I,J)
         EPA_WD_BF_SO2 (I,J) = EPA_WD_BF_SO2 (I,J) * USA_MASK(I,J)
         EPA_WD_BF_SO4 (I,J) = EPA_WD_BF_SO4 (I,J) * USA_MASK(I,J)
         EPA_WD_BF_NH3 (I,J) = EPA_WD_BF_NH3 (I,J) * USA_MASK(I,J)

         ! Weekend avg biofuel    
         EPA_WE_BF_NOX (I,J) = EPA_WE_BF_NOX (I,J) * USA_MASK(I,J)
         EPA_WE_BF_CO  (I,J) = EPA_WE_BF_CO  (I,J) * USA_MASK(I,J)
         EPA_WE_BF_ALK4(I,J) = EPA_WE_BF_ALK4(I,J) * USA_MASK(I,J)
         EPA_WE_BF_ACET(I,J) = EPA_WE_BF_ACET(I,J) * USA_MASK(I,J)
         EPA_WE_BF_MEK (I,J) = EPA_WE_BF_MEK (I,J) * USA_MASK(I,J)
         EPA_WE_BF_PRPE(I,J) = EPA_WE_BF_PRPE(I,J) * USA_MASK(I,J)   
         EPA_WE_BF_C3H8(I,J) = EPA_WE_BF_C3H8(I,J) * USA_MASK(I,J)
         EPA_WE_BF_CH2O(I,J) = EPA_WE_BF_CH2O(I,J) * USA_MASK(I,J)
         EPA_WE_BF_C2H6(I,J) = EPA_WE_BF_C2H6(I,J) * USA_MASK(I,J)
         EPA_WE_BF_SO2 (I,J) = EPA_WE_BF_SO2 (I,J) * USA_MASK(I,J)
         EPA_WE_BF_SO4 (I,J) = EPA_WE_BF_SO4 (I,J) * USA_MASK(I,J)
         EPA_WE_BF_NH3 (I,J) = EPA_WE_BF_NH3 (I,J) * USA_MASK(I,J)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Print totals to the log file
      CALL TOTAL_ANTHRO_TG( THISMONTH )
      CALL TOTAL_BIOFUEL_TG( THISMONTH )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE EMISS_EPA_NEI

!------------------------------------------------------------------------------

      SUBROUTINE READ_EPA( FILENAME, NOX,  CO,   ALK4, ACET, MEK,       
     &                     PRPE,     C3H8, CH2O, C2H6, SO2,  SO4,  NH3 )
!
!******************************************************************************
!  Subroutine READ_EPA reads an EPA data file (biomass or anthro) from disk.
!  The entire file is read through on one pass for better I/O optimization.
!  (rch, bmy, 7/1/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Name of anthro or biomass file to read
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) NOx      (REAL*4   ) : Array for NOx  anthro or biomass data
!  (3 ) CO       (REAL*4   ) : Array for CO   anthro or biomass data
!  (4 ) ALK4     (REAL*4   ) : Array for ALK4 anthro or biomass data
!  (5 ) ACET     (REAL*4   ) : Array for ACET anthro or biomass data
!  (6 ) MEK      (REAL*4   ) : Array for MEK  anthro or biomass data
!  (7 ) PRPE     (REAL*4   ) : Array for PRPE anthro or biomass data 
!  (8 ) C3H8     (REAL*4   ) : Array for C3H8 anthro or biomass data 
!  (9 ) CH2O     (REAL*4   ) : Array for CH2O anthro or biomass data 
!  (10) C2H6     (REAL*4   ) : Array for C2H6 anthro or biomass data
!  (11) NH3      (REAL*4   ) : Array for NH3  anthro or biomass data
!  (12) SO2      (REAL*4   ) : Array for SO4  anthro or biomass data 
!  (13) SO4      (REAL*4   ) : Array for SO4  anthro or biomass data 
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,    ONLY : OPEN_BPCH2_FOR_READ
      USE FILE_MOD,     ONLY : IU_FILE, IOERROR
      USE TRANSFER_MOD, ONLY : TRANSFER_2D
   
#     include "CMN_SIZE"             ! Size parameters
  
      ! Arguments
      CHARACTER(LEN=*), INTENT(IN)    :: FILENAME
      REAL*4,           INTENT(INOUT) :: NOX(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: CO(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: ALK4(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: ACET(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: MEK(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: PRPE(IIPAR,JJPAR)      
      REAL*4,           INTENT(INOUT) :: C3H8(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: CH2O(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: C2H6(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: SO2(IIPAR,JJPAR)  
      REAL*4,           INTENT(INOUT) :: SO4(IIPAR,JJPAR)  
      REAL*4,           INTENT(INOUT) :: NH3(IIPAR,JJPAR)

      ! Local variables
      INTEGER                         :: I,  J,  L,  N,  IOS
      INTEGER                         :: NTRACER,   NSKIP
      INTEGER                         :: HALFPOLAR, CENTER180
      INTEGER                         :: NI,        NJ,        NL
      INTEGER                         :: IFIRST,    JFIRST,    LFIRST
      REAL*4                          :: ARRAY(IGLOB,JGLOB,1)
      REAL*4                          :: LONRES,    LATRES
      REAL*8                          :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)               :: MODELNAME
      CHARACTER(LEN=40)               :: CATEGORY
      CHARACTER(LEN=40)               :: UNIT     
      CHARACTER(LEN=40)               :: RESERVED

      !=================================================================
      ! READ_EPA begins here!
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_EPA: Reading ', a )

      ! Open file
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )

      ! Read the entire file in one pass (for I/O optimization)
      DO 

         ! Read 1st data block header line
         READ( IU_FILE, IOSTAT=IOS ) 
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

         ! Check for EOF or errors
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_data:2' )

         ! Read 2nd data block header line
         READ( IU_FILE, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP

         ! Error check
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_data:3' )

         ! Read data
         READ( IU_FILE, IOSTAT=IOS ) 
     &        ( ( ( ARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         ! Error check
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_data:4' )

         !==============================================================
         ! Save into tracer arrays
         !==============================================================
         SELECT CASE ( NTRACER )
            CASE( 1  )
               CALL TRANSFER_2D( ARRAY(:,:,1), NOx  )
            CASE( 4  )
               CALL TRANSFER_2D( ARRAY(:,:,1), CO   )
            CASE( 5  )
               CALL TRANSFER_2D( ARRAY(:,:,1), ALK4 )
            CASE( 9  )
               CALL TRANSFER_2D( ARRAY(:,:,1), ACET )
            CASE( 10 )
               CALL TRANSFER_2D( ARRAY(:,:,1), MEK  )
            CASE( 18 )
               CALL TRANSFER_2D( ARRAY(:,:,1), PRPE )
            CASE( 19 )
               CALL TRANSFER_2D( ARRAY(:,:,1), C3H8 )
            CASE( 20 )
               CALL TRANSFER_2D( ARRAY(:,:,1), CH2O )
            CASE( 21 )
               CALL TRANSFER_2D( ARRAY(:,:,1), C2H6 )
            CASE( 26 )
               CALL TRANSFER_2D( ARRAY(:,:,1), SO2  )
            CASE( 27 )
               CALL TRANSFER_2D( ARRAY(:,:,1), SO4  )
            CASE( 29 )
               CALL TRANSFER_2D( ARRAY(:,:,1), NH3  )
            CASE DEFAULT
               ! Nothing
         END SELECT

      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE READ_EPA

!------------------------------------------------------------------------------

      SUBROUTINE READ_USA_MASK
!
!******************************************************************************
!  Subroutine READ_USA_MASK reads the USA mask from disk.   The USA mask is
!  the fraction of the grid box (I,J) which lies w/in the continental USA.
!  (rch, bmy, 11/10/04, 10/3/05)
!
!  NOTES:
!  (1 ) Now can read data for GEOS and GCAP grids (bmy, 8/16/05)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! Reference to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      REAL*4             :: ARRAY(IGLOB,JGLOB,1)
      REAL*8             :: XTAU
      CHARACTER(LEN=255) :: FILENAME 

      !=================================================================
      ! READ_USA_MASK begins here!
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR )           //
     &           'EPA_NEI_200411/usa_mask.' // GET_NAME_EXT_2D() //
     &           '.'                        // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_USA_MASK: Reading ', a )

      ! Get TAU0 for Jan 1985
      XTAU  = GET_TAU0( 1, 1, 1985 )

      ! USA mask is stored in the bpch file as #2
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2, 
     &                 XTAU,      IGLOB,    JGLOB,     
     &                 1,         ARRAY,    QUIET=.TRUE. ) 

      ! Cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), USA_MASK )

      ! Return to calling program
      END SUBROUTINE READ_USA_MASK

!------------------------------------------------------------------------------

      SUBROUTINE TOTAL_ANTHRO_TG( THISMONTH )
!
!******************************************************************************
!  Subroutine TOTAL_ANTHRO_TG prints the amount of EPA/NEI anthropogenic
!  emissions that are emitted each month in Tg or Tg C. 
!  (rch, bmy, 11/10/04, 1/25/05)
!  
!  Arguments as Input:
!  ============================================================================
!  (1  ) FFARRAY  (REAL*8 ) : Fossil Fuel CO emissions [molec (C)/cm2/month]
!  (2-4) IX,JX,LX (INTEGER) : Dimensions of FFARRAY 
!  (5  ) MOLWT    (REAL*8 ) : Molecular wt [kg/mole] for the given tracer
!  (6  ) NAME     (REAL*8 ) : Tracer name
!  (7  ) NSEASON  (INTEGER) : Number of the season, for seasonal NOx/SOX
!
!  NOTES:
!  (1) Scale factors were determined by Jennifer Logan (jal@io.harvard.edu),
!      Bryan Duncan (bnd@io.harvard.edu), and Daniel Jacob (djj@io.harvard.edu)
!  (2) Now replace DXYP(J)*1d4 with routine GET_AREA_CM2 from "grid_mod.f".
!       (bmy, 2/4/03)
!  (3) Prevent out of bounds error when tracers are undefined (bmy, 1/25/05)
!  (4) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACERID_MOD, ONLY : IDTACET, IDTALK4, IDTC2H6, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTCH2O, IDTCO,   IDTMEK,  IDTNOX
      USE TRACERID_MOD, ONLY : IDTNH3,  IDTPRPE, IDTSO2,  IDTSO4  

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_O3"    ! FMOL

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, J
      REAL*8              :: WD_NOX,  WD_CO,   WD_ALK4, WD_ACET
      REAL*8              :: WD_MEK,  WD_PRPE, WD_C2H6, WD_C3H8
      REAL*8              :: WD_CH2O, WD_NH3,  WD_SO2,  WD_SO4,  A
      REAL*8              :: WE_NOX,  WE_CO,   WE_ALK4, WE_ACET
      REAL*8              :: WE_MEK,  WE_PRPE, WE_C2H6, WE_C3H8
      REAL*8              :: WE_CH2O, WE_NH3,  WE_SO2,  WE_SO4
      REAL*8              :: F_NOX,   F_CO,    F_ALK4,  F_ACET
      REAL*8              :: F_MEK,   F_PRPE,  F_C2H6,  F_C3H8
      REAL*8              :: F_CH2O,  F_SO2,   F_SO4,   F_NH3
      CHARACTER(LEN=6)    :: UNIT

      ! Days per month
      INTEGER             :: D(12) = (/ 31, 28, 31, 30, 31, 30,
     &                                  31, 31, 30, 31, 30, 31 /)

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Summing variables for weekday avg anthro
      WD_NOX  = 0d0 
      WD_CO   = 0d0 
      WD_ALK4 = 0d0 
      WD_ACET = 0d0 
      WD_MEK  = 0d0 
      WD_PRPE = 0d0 
      WD_C2H6 = 0d0 
      WD_C3H8 = 0d0 
      WD_CH2O = 0d0 
      WD_NH3  = 0d0 
      WD_SO2  = 0d0 
      WD_SO4  = 0d0 

      ! Summing variables for weekend avg anthro
      WE_NOX  = 0d0 
      WE_CO   = 0d0 
      WE_ALK4 = 0d0 
      WE_ACET = 0d0 
      WE_MEK  = 0d0 
      WE_PRPE = 0d0 
      WE_C2H6 = 0d0 
      WE_C3H8 = 0d0 
      WE_CH2O = 0d0 
      WE_NH3  = 0d0 
      WE_SO2  = 0d0 
      WE_SO4  = 0d0 

      ! Molecular weights
      F_NOX   = 0d0   
      F_CO    = 0d0 
      F_ALK4  = 0d0 
      F_ACET  = 0d0 
      F_MEK   = 0d0 
      F_PRPE  = 0d0 
      F_C2H6  = 0d0 
      F_C3H8  = 0d0 
      F_CH2O  = 0d0 
      F_SO2   = 0d0 
      F_SO4   = 0d0 
      F_NH3   = 0d0 

      ! Prevent array out of bounds error for undefined tracers
      IF ( IDTNOX  > 0 ) F_NOX  = FMOL(IDTNOX )
      IF ( IDTCO   > 0 ) F_CO   = FMOL(IDTCO  )
      IF ( IDTALK4 > 0 ) F_ALK4 = FMOL(IDTALK4)
      IF ( IDTACET > 0 ) F_ACET = FMOL(IDTACET)
      IF ( IDTMEK  > 0 ) F_MEK  = FMOL(IDTMEK )
      IF ( IDTPRPE > 0 ) F_PRPE = FMOL(IDTPRPE)
      IF ( IDTC2H6 > 0 ) F_C2H6 = FMOL(IDTC2H6)
      IF ( IDTC3H8 > 0 ) F_C3H8 = FMOL(IDTC3H8)
      IF ( IDTCH2O > 0 ) F_CH2O = FMOL(IDTCH2O)
      IF ( IDTSO2  > 0 ) F_SO2  = FMOL(IDTSO2 )
      IF ( IDTSO4  > 0 ) F_SO4  = FMOL(IDTSO4 )
      IF ( IDTNH3  > 0 ) F_NH3  = FMOL(IDTNH3 )

      !=================================================================
      ! Sum anthropogenic emissions
      !=================================================================

      ! Loop over latitudes
      DO J = 1, JJPAR
            
         ! Surface area [cm2] * seconds in this month / AVOGADRO's number
         ! Also multiply by the factor 1d-9 to convert kg to Tg
         A = GET_AREA_CM2( J ) * ( D(THISMONTH) * 86400d-9 ) / 6.0225d23
         
         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Weekday avg emissions
            WD_NOX  = WD_NOX  + EPA_WD_AN_NOX (I,J) * A * F_NOX 
            WD_CO   = WD_CO   + EPA_WD_AN_CO  (I,J) * A * F_CO  
            WD_ALK4 = WD_ALK4 + EPA_WD_AN_ALK4(I,J) * A * F_ALK4
            WD_ACET = WD_ACET + EPA_WD_AN_ACET(I,J) * A * F_ACET
            WD_MEK  = WD_MEK  + EPA_WD_AN_MEK (I,J) * A * F_MEK 
            WD_PRPE = WD_PRPE + EPA_WD_AN_PRPE(I,J) * A * F_PRPE
            WD_C2H6 = WD_C2H6 + EPA_WD_AN_C2H6(I,J) * A * F_C2H6
            WD_C3H8 = WD_C3H8 + EPA_WD_AN_C3H8(I,J) * A * F_C3H8
            WD_CH2O = WD_CH2O + EPA_WD_AN_CH2O(I,J) * A * F_CH2O
            WD_SO2  = WD_SO2  + EPA_WD_AN_SO2 (I,J) * A * F_SO2 
            WD_SO4  = WD_SO4  + EPA_WD_AN_SO4 (I,J) * A * F_SO4 
            WD_NH3  = WD_NH3  + EPA_WD_AN_NH3 (I,J) * A * F_NH3 

            ! Weekend avg emissions
            WE_NOX  = WE_NOX  + EPA_WE_AN_NOX (I,J) * A * F_NOX 
            WE_CO   = WE_CO   + EPA_WE_AN_CO  (I,J) * A * F_CO  
            WE_ALK4 = WE_ALK4 + EPA_WE_AN_ALK4(I,J) * A * F_ALK4
            WE_ACET = WE_ACET + EPA_WE_AN_ACET(I,J) * A * F_ACET
            WE_MEK  = WE_MEK  + EPA_WE_AN_MEK (I,J) * A * F_MEK 
            WE_PRPE = WE_PRPE + EPA_WE_AN_PRPE(I,J) * A * F_PRPE
            WE_C2H6 = WE_C2H6 + EPA_WE_AN_C2H6(I,J) * A * F_C2H6
            WE_C3H8 = WE_C3H8 + EPA_WE_AN_C3H8(I,J) * A * F_C3H8
            WE_CH2O = WE_CH2O + EPA_WE_AN_CH2O(I,J) * A * F_CH2O
            WE_SO2  = WE_SO2  + EPA_WE_AN_SO2 (I,J) * A * F_SO2 
            WE_SO4  = WE_SO4  + EPA_WE_AN_SO4 (I,J) * A * F_SO4 
            WE_NH3  = WE_NH3  + EPA_WE_AN_NH3 (I,J) * A * F_NH3 

         ENDDO
      ENDDO
 
      !=================================================================
      ! Print info
      !=================================================================
      
      ! Weekday avg anthro
      WRITE( 6, '(a)' )
      WRITE( 6, 100   ) 'NOx ', THISMONTH, WD_NOX,  '  '
      WRITE( 6, 100   ) 'CO  ', THISMONTH, WD_CO,   '  '
      WRITE( 6, 100   ) 'ALK4', THISMONTH, WD_ALK4, ' C'
      WRITE( 6, 100   ) 'ACET', THISMONTH, WD_ACET, ' C'
      WRITE( 6, 100   ) 'MEK ', THISMONTH, WD_MEK,  ' C'
      WRITE( 6, 100   ) 'PRPE', THISMONTH, WD_PRPE, ' C'
      WRITE( 6, 100   ) 'C3H8', THISMONTH, WD_C3H8, ' C'
      WRITE( 6, 100   ) 'CH2O', THISMONTH, WD_CH2O, '  '
      WRITE( 6, 100   ) 'C2H6', THISMONTH, WD_C2H6, ' C'
      WRITE( 6, 100   ) 'SO2 ', THISMONTH, WD_SO2,  '  '
      WRITE( 6, 100   ) 'SO4 ', THISMONTH, WD_SO4,  '  '
      WRITE( 6, 100   ) 'NH3 ', THISMONTH, WD_NH3,  '  '
 100  FORMAT( 'Total weekday avg anthro ', a4, ' for 1999/', 
     &         i2.2, ': ', f13.6, ' Tg', a2 )

      ! Weekend avg anthro
      WRITE( 6, '(a)' )
      WRITE( 6, 110   ) 'NOx ', THISMONTH, WE_NOX,  '  '
      WRITE( 6, 110   ) 'CO  ', THISMONTH, WE_CO,   '  '
      WRITE( 6, 110   ) 'ALK4', THISMONTH, WE_ALK4, ' C'
      WRITE( 6, 110   ) 'ACET', THISMONTH, WE_ACET, ' C'
      WRITE( 6, 110   ) 'MEK ', THISMONTH, WE_MEK,  ' C'
      WRITE( 6, 110   ) 'PRPE', THISMONTH, WE_PRPE, ' C'
      WRITE( 6, 110   ) 'C3H8', THISMONTH, WE_C3H8, ' C'
      WRITE( 6, 110   ) 'CH2O', THISMONTH, WE_CH2O, '  '
      WRITE( 6, 110   ) 'C2H6', THISMONTH, WE_C2H6, ' C'
      WRITE( 6, 110   ) 'SO2 ', THISMONTH, WE_SO2,  '  '
      WRITE( 6, 110   ) 'SO4 ', THISMONTH, WE_SO4,  '  '
      WRITE( 6, 110   ) 'NH3 ', THISMONTH, WE_NH3,  '  '
 110  FORMAT( 'Total weekend avg anthro ', a4, ' for 1999/', 
     &         i2.2, ': ', f13.6, ' Tg', a2 )

      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_TG

!------------------------------------------------------------------------------

      SUBROUTINE TOTAL_BIOFUEL_TG( THISMONTH )
!
!******************************************************************************
!  Subroutine TOTAL_BIOFUEL_TG prints the amount of EPA/NEI biofuel emissions
!  that are emitted each month in Tg or Tg C. (rch, bmy, 11/10/04, 1/26/05)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number
!
!  NOTES:
!  (1 ) Prevent out of bounds error when tracers are undefined (bmy, 1/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACERID_MOD, ONLY : IDTACET, IDTALK4, IDTC2H6, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTCH2O, IDTCO,   IDTMEK,  IDTNOX
      USE TRACERID_MOD, ONLY : IDTNH3,  IDTPRPE, IDTSO2,  IDTSO4  

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_O3"    ! FMOL

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, J
      REAL*8              :: WD_NOX,  WD_CO,   WD_ALK4, WD_ACET
      REAL*8              :: WD_MEK,  WD_PRPE, WD_C2H6, WD_C3H8
      REAL*8              :: WD_CH2O, WD_NH3,  WD_SO2,  WD_SO4,  A
      REAL*8              :: WE_NOX,  WE_CO,   WE_ALK4, WE_ACET
      REAL*8              :: WE_MEK,  WE_PRPE, WE_C2H6, WE_C3H8
      REAL*8              :: WE_CH2O, WE_NH3,  WE_SO2,  WE_SO4
      REAL*8              :: F_NOX,   F_CO,    F_ALK4,  F_ACET
      REAL*8              :: F_MEK,   F_PRPE,  F_C2H6,  F_C3H8
      REAL*8              :: F_CH2O,  F_SO2,   F_SO4,   F_NH3
      CHARACTER(LEN=6)    :: UNIT

      ! Days per month
      INTEGER             :: D(12) = (/ 31, 28, 31, 30, 31, 30,
     &                                  31, 31, 30, 31, 30, 31 /)

      !=================================================================
      ! TOTAL_BIOFUEL_TG begins here!
      !=================================================================

      ! Summing variables for weekday avg anthro
      WD_NOX  = 0d0 
      WD_CO   = 0d0 
      WD_ALK4 = 0d0 
      WD_ACET = 0d0 
      WD_MEK  = 0d0 
      WD_PRPE = 0d0 
      WD_C2H6 = 0d0 
      WD_C3H8 = 0d0 
      WD_CH2O = 0d0 
      WD_NH3  = 0d0 
      WD_SO2  = 0d0 
      WD_SO4  = 0d0 

      ! Summing variables for weekend avg anthro
      WE_NOX  = 0d0 
      WE_CO   = 0d0 
      WE_ALK4 = 0d0 
      WE_ACET = 0d0 
      WE_MEK  = 0d0 
      WE_PRPE = 0d0 
      WE_C2H6 = 0d0 
      WE_C3H8 = 0d0 
      WE_CH2O = 0d0 
      WE_NH3  = 0d0 
      WE_SO2  = 0d0 
      WE_SO4  = 0d0 

      ! Molecular weights
      F_NOX   = 0d0   
      F_CO    = 0d0 
      F_ALK4  = 0d0 
      F_ACET  = 0d0 
      F_MEK   = 0d0 
      F_PRPE  = 0d0 
      F_C2H6  = 0d0 
      F_C3H8  = 0d0 
      F_CH2O  = 0d0 
      F_SO2   = 0d0 
      F_SO4   = 0d0 
      F_NH3   = 0d0 

      ! Prevent array out of bounds error for undefined tracers
      IF ( IDTNOX  > 0 ) F_NOX  = FMOL(IDTNOX )
      IF ( IDTCO   > 0 ) F_CO   = FMOL(IDTCO  )
      IF ( IDTALK4 > 0 ) F_ALK4 = FMOL(IDTALK4)
      IF ( IDTACET > 0 ) F_ACET = FMOL(IDTACET)
      IF ( IDTMEK  > 0 ) F_MEK  = FMOL(IDTMEK )
      IF ( IDTPRPE > 0 ) F_PRPE = FMOL(IDTPRPE)
      IF ( IDTC2H6 > 0 ) F_C2H6 = FMOL(IDTC2H6)
      IF ( IDTC3H8 > 0 ) F_C3H8 = FMOL(IDTC3H8)
      IF ( IDTCH2O > 0 ) F_CH2O = FMOL(IDTCH2O)
      IF ( IDTSO2  > 0 ) F_SO2  = FMOL(IDTSO2 )
      IF ( IDTSO4  > 0 ) F_SO4  = FMOL(IDTSO4 )
      IF ( IDTNH3  > 0 ) F_NH3  = FMOL(IDTNH3 )

      !=================================================================
      ! Sum biofuel emissions
      !=================================================================

      ! Loop over surface boxes
      DO J = 1, JJPAR

         ! Surface area [cm2] * seconds in this month / AVOGADRO's number
         ! Also multiply by the factor 1d-9 to convert kg to Tg
         A = GET_AREA_CM2( J ) * ( D(THISMONTH) * 86400d-9 ) / 6.0225d23

         DO I = 1, IIPAR

            ! Weekday avg emissions
            WD_NOX  = WD_NOX  + EPA_WD_BF_NOX (I,J) * A * F_NOX
            WD_CO   = WD_CO   + EPA_WD_BF_CO  (I,J) * A * F_CO  
            WD_ALK4 = WD_ALK4 + EPA_WD_BF_ALK4(I,J) * A * F_ALK4
            WD_ACET = WD_ACET + EPA_WD_BF_ACET(I,J) * A * F_ACET
            WD_MEK  = WD_MEK  + EPA_WD_BF_MEK (I,J) * A * F_MEK 
            WD_PRPE = WD_PRPE + EPA_WD_BF_PRPE(I,J) * A * F_PRPE
            WD_C2H6 = WD_C2H6 + EPA_WD_BF_C2H6(I,J) * A * F_C2H6
            WD_C3H8 = WD_C3H8 + EPA_WD_BF_C3H8(I,J) * A * F_C3H8
            WD_CH2O = WD_CH2O + EPA_WD_BF_CH2O(I,J) * A * F_CH2O
            WD_SO2  = WD_SO2  + EPA_WD_BF_SO2 (I,J) * A * F_SO2 
            WD_SO4  = WD_SO4  + EPA_WD_BF_SO4 (I,J) * A * F_SO4 
            WD_NH3  = WD_NH3  + EPA_WD_BF_NH3 (I,J) * A * F_NH3 
            
            ! Weekend avg emissions
            WE_NOX  = WE_NOX  + EPA_WE_BF_NOX (I,J) * A * F_NOX 
            WE_CO   = WE_CO   + EPA_WE_BF_CO  (I,J) * A * F_CO  
            WE_ALK4 = WE_ALK4 + EPA_WE_BF_ALK4(I,J) * A * F_ALK4
            WE_ACET = WE_ACET + EPA_WE_BF_ACET(I,J) * A * F_ACET
            WE_MEK  = WE_MEK  + EPA_WE_BF_MEK (I,J) * A * F_MEK 
            WE_PRPE = WE_PRPE + EPA_WE_BF_PRPE(I,J) * A * F_PRPE
            WE_C2H6 = WE_C2H6 + EPA_WE_BF_C2H6(I,J) * A * F_C2H6
            WE_C3H8 = WE_C3H8 + EPA_WE_BF_C3H8(I,J) * A * F_C3H8
            WE_CH2O = WE_CH2O + EPA_WE_BF_CH2O(I,J) * A * F_CH2O
            WE_SO2  = WE_SO2  + EPA_WE_BF_SO2 (I,J) * A * F_SO2 
            WE_SO4  = WE_SO4  + EPA_WE_BF_SO4 (I,J) * A * F_SO4 
            WE_NH3  = WE_NH3  + EPA_WE_BF_NH3 (I,J) * A * F_NH3 
            
         ENDDO
      ENDDO
 
      !=================================================================
      ! Print info
      !=================================================================
      
      ! Weekday avg biofuel
      WRITE( 6, '(a)' )
      WRITE( 6, 100   ) 'NOx ', THISMONTH, WD_NOX,  '  '
      WRITE( 6, 100   ) 'CO  ', THISMONTH, WD_CO,   '  '
      WRITE( 6, 100   ) 'ALK4', THISMONTH, WD_ALK4, ' C'
      WRITE( 6, 100   ) 'ACET', THISMONTH, WD_ACET, ' C'
      WRITE( 6, 100   ) 'MEK ', THISMONTH, WD_MEK,  ' C'
      WRITE( 6, 100   ) 'PRPE', THISMONTH, WD_PRPE, ' C'
      WRITE( 6, 100   ) 'C3H8', THISMONTH, WD_C3H8, ' C'
      WRITE( 6, 100   ) 'CH2O', THISMONTH, WD_CH2O, '  '
      WRITE( 6, 100   ) 'C2H6', THISMONTH, WD_C2H6, ' C'
      WRITE( 6, 100   ) 'SO2 ', THISMONTH, WD_SO2,  '  '
      WRITE( 6, 100   ) 'SO4 ', THISMONTH, WD_SO4,  '  '
      WRITE( 6, 100   ) 'NH3 ', THISMONTH, WD_NH3,  '  '
 100  FORMAT( 'Total weekday avg biofuel ', a4, ' for 1999/', 
     &         i2.2, ': ', f13.6, ' Tg', a2 )

      ! Weekend avg biofuel
      WRITE( 6, '(a)' )
      WRITE( 6, 110   ) 'NOx ', THISMONTH, WE_NOX,  '  '
      WRITE( 6, 110   ) 'CO  ', THISMONTH, WE_CO,   '  '
      WRITE( 6, 110   ) 'ALK4', THISMONTH, WE_ALK4, ' C'
      WRITE( 6, 110   ) 'ACET', THISMONTH, WE_ACET, ' C'
      WRITE( 6, 110   ) 'MEK ', THISMONTH, WE_MEK,  ' C'
      WRITE( 6, 110   ) 'PRPE', THISMONTH, WE_PRPE, ' C'
      WRITE( 6, 110   ) 'C3H8', THISMONTH, WE_C3H8, ' C'
      WRITE( 6, 110   ) 'CH2O', THISMONTH, WE_CH2O, '  '
      WRITE( 6, 110   ) 'C2H6', THISMONTH, WE_C2H6, ' C'
      WRITE( 6, 110   ) 'SO2 ', THISMONTH, WE_SO2,  '  '
      WRITE( 6, 110   ) 'SO4 ', THISMONTH, WE_SO4,  '  '
      WRITE( 6, 110   ) 'NH3 ', THISMONTH, WE_NH3,  '  '
 110  FORMAT( 'Total weekend avg biofuel ', a4, ' for 1999/', 
     &         i2.2, ': ', f13.6, ' Tg', a2 )

      ! Return to calling program
      END SUBROUTINE TOTAL_BIOFUEL_TG

!------------------------------------------------------------------------------

      FUNCTION GET_USA_MASK( I, J ) RESULT( USA )
!
!******************************************************************************
!  Function GET_USA_MASK returns the value of the USA mask (i.e. the fraction
!  of a grid box which lies w/in the continental USA) at a given (I,J)
!  location. (rch, bmy, 11/10/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) : Longitude index
!  (2 ) J (INTEGER) : Latitude index
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Local variables
      REAL*8              :: USA 

      !=================================================================
      ! GET_USA_MASK begins here!
      !=================================================================
      USA = USA_MASK(I,J)
      
      ! Return to calling program
      END FUNCTION GET_USA_MASK

!------------------------------------------------------------------------------

      FUNCTION GET_EPA_ANTHRO( I, J, N, WEEKDAY ) RESULT( EPA_NEI )
!
!******************************************************************************
!  Function GET_EPA_ANTHRO returns the EPA weekday avg or weekend avg 
!  anthropogenic emissions at a (I,J) location. (rch, bmy, 11/10/04, 10/3/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I       (INTEGER) : Longitude index
!  (2 ) J       (INTEGER) : Latitude index
!  (3 ) N       (INTEGER) : Tracer index (i.e. as listed in "inptr.ctm")
!  (4 ) WEEKDAY (LOGICAL) : Fla for weekday (=T) or weekend (=F) emissions
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY.  Also remove reference
!        to BPCH2_MOD and TRACERID_MOD, they're not needed.  (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD, ONLY : IDTACET, IDTALK4, IDTC2H6, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTCH2O, IDTCO,   IDTMEK,  IDTNOX
      USE TRACERID_MOD, ONLY : IDTNH3,  IDTPRPE, IDTSO2,  IDTSO4  

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      LOGICAL, INTENT(IN) :: WEEKDAY
      INTEGER, INTENT(IN) :: I, J, N
      
      ! Local variables
      REAL*8              :: EPA_NEI

      !=================================================================
      ! GET_EPA_ANTHRO begins here!
      !=================================================================

      ! Return either weekday or weekend avg emissions
      IF ( WEEKDAY ) THEN 

         !--------------------
         ! Weekday avg anthro
         !--------------------
         IF ( N == IDTNOX ) THEN
            EPA_NEI = EPA_WD_AN_NOX(I,J)

         ELSE IF ( N == IDTCO ) THEN 
            EPA_NEI = EPA_WD_AN_CO(I,J)
            
         ELSE IF ( N == IDTALK4 ) THEN
            EPA_NEI = EPA_WD_AN_ALK4(I,J)

         ELSE IF ( N == IDTACET ) THEN
            EPA_NEI = EPA_WD_AN_ACET(I,J)

         ELSE IF ( N == IDTMEK ) THEN
            EPA_NEI = EPA_WD_AN_MEK(I,J)

         ELSE IF ( N == IDTPRPE ) THEN
            EPA_NEI = EPA_WD_AN_PRPE(I,J)

         ELSE IF ( N == IDTC3H8 ) THEN
            EPA_NEI = EPA_WD_AN_C3H8(I,J)

         ELSE IF ( N == IDTCH2O ) THEN
            EPA_NEI = EPA_WD_AN_CH2O(I,J)

         ELSE IF ( N == IDTC2H6 ) THEN
            EPA_NEI = EPA_WD_AN_C2H6(I,J)

         ELSE IF ( N == IDTSO2 ) THEN
            EPA_NEI = EPA_WD_AN_SO2(I,J)

         ELSE IF ( N == IDTSO4 ) THEN
            EPA_NEI = EPA_WD_AN_SO4(I,J)

         ELSE IF ( N == IDTNH3 ) THEN
            EPA_NEI = EPA_WD_AN_NH3(I,J)

         ELSE
            EPA_NEI = 0d0
            
         ENDIF

      ELSE

         !--------------------
         ! Weekend avg anthro
         !--------------------
         IF (  N == IDTNOX  ) THEN
            EPA_NEI = EPA_WE_AN_NOX(I,J)

         ELSE IF ( N == IDTCO ) THEN 
            EPA_NEI = EPA_WE_AN_CO(I,J)
            
         ELSE IF ( N == IDTALK4 ) THEN
            EPA_NEI = EPA_WE_AN_ALK4(I,J)

         ELSE IF ( N == IDTACET ) THEN
            EPA_NEI = EPA_WE_AN_ACET(I,J)

         ELSE IF ( N == IDTMEK ) THEN
            EPA_NEI = EPA_WE_AN_MEK(I,J)

         ELSE IF ( N == IDTPRPE ) THEN
            EPA_NEI = EPA_WE_AN_PRPE(I,J)

         ELSE IF ( N == IDTC3H8 ) THEN
            EPA_NEI = EPA_WE_AN_C3H8(I,J)

         ELSE IF ( N == IDTCH2O ) THEN
            EPA_NEI = EPA_WE_AN_CH2O(I,J)

         ELSE IF ( N == IDTC2H6 ) THEN
            EPA_NEI = EPA_WE_AN_C2H6(I,J)

         ELSE IF ( N == IDTSO2 ) THEN
            EPA_NEI = EPA_WE_AN_SO2(I,J)

         ELSE IF ( N == IDTSO4 ) THEN
            EPA_NEI = EPA_WE_AN_SO4(I,J)

         ELSE IF ( N == IDTNH3 ) THEN
            EPA_NEI = EPA_WE_AN_NH3(I,J)

         ELSE
            EPA_NEI = 0d0
            
         ENDIF

      ENDIF

      ! Return to calling program
      END FUNCTION GET_EPA_ANTHRO

!------------------------------------------------------------------------------

      FUNCTION GET_EPA_BIOFUEL( I, J, N, WEEKDAY ) RESULT( EPA_NEI )
!
!******************************************************************************
!  Function GET_EPA_BIOFUEL returns the EPA weekday avg or weekend avg 
!  biofuel emissions at a (I,J) location. (rch, bmy, 11/10/04, 10/3/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I       (INTEGER) : Longitude index
!  (2 ) J       (INTEGER) : Latitude index
!  (3 ) N       (INTEGER) : Tracer index (i.e. as listed in "inptr.ctm")
!  (4 ) WEEKDAY (LOGICAL) : Fla for weekday (=T) or weekend (=F) emissions
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY.  Also remove reference
!        to BPCH2_MOD and TRACERID_MOD, they're not needed.  (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD, ONLY : IDTACET, IDTALK4, IDTC2H6, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTCH2O, IDTCO,   IDTMEK,  IDTNOX
      USE TRACERID_MOD, ONLY : IDTNH3,  IDTPRPE, IDTSO2,  IDTSO4  

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      LOGICAL, INTENT(IN) :: WEEKDAY
      INTEGER, INTENT(IN) :: I, J, N
      
      ! Local variables
      REAL*8              :: EPA_NEI

      !=================================================================
      ! GET_EPA_BIOFUEL begins here!
      !=================================================================

      ! Return either weekday or weekend avg emissions
      IF ( WEEKDAY ) THEN 

         !---------------------
         ! Weekday avg biofuel
         !---------------------
         IF ( N == IDTNOX ) THEN
            EPA_NEI = EPA_WD_BF_NOX(I,J)

         ELSE IF ( N == IDTCO ) THEN 
            EPA_NEI = EPA_WD_BF_CO(I,J)
            
         ELSE IF ( N == IDTALK4 ) THEN
            EPA_NEI = EPA_WD_BF_ALK4(I,J)

         ELSE IF ( N == IDTACET ) THEN
            EPA_NEI = EPA_WD_BF_ACET(I,J)

         ELSE IF ( N == IDTMEK ) THEN
            EPA_NEI = EPA_WD_BF_MEK(I,J)

         ELSE IF ( N == IDTPRPE ) THEN
            EPA_NEI = EPA_WD_BF_PRPE(I,J)

         ELSE IF ( N == IDTC3H8 ) THEN
            EPA_NEI = EPA_WD_BF_C3H8(I,J)

         ELSE IF ( N == IDTCH2O ) THEN
            EPA_NEI = EPA_WD_BF_CH2O(I,J)

         ELSE IF ( N == IDTC2H6 ) THEN
            EPA_NEI = EPA_WD_BF_C2H6(I,J)

         ELSE IF ( N == IDTSO2 ) THEN
            EPA_NEI = EPA_WD_BF_SO2(I,J)

         ELSE IF ( N == IDTSO4 ) THEN
            EPA_NEI = EPA_WD_BF_SO4(I,J)

         ELSE IF ( N == IDTNH3 ) THEN
            EPA_NEI = EPA_WD_BF_NH3(I,J)

         ELSE
            EPA_NEI = 0d0
            
         ENDIF

      ELSE

         !---------------------
         ! Weekend avg biofuel
         !---------------------
         IF (  N == IDTNOX  ) THEN
            EPA_NEI = EPA_WE_BF_NOX(I,J)

         ELSE IF ( N == IDTCO ) THEN 
            EPA_NEI = EPA_WE_BF_CO(I,J)
            
         ELSE IF ( N == IDTALK4 ) THEN
            EPA_NEI = EPA_WE_BF_ALK4(I,J)

         ELSE IF ( N == IDTACET ) THEN
            EPA_NEI = EPA_WE_BF_ACET(I,J)

         ELSE IF ( N == IDTMEK ) THEN
            EPA_NEI = EPA_WE_BF_MEK(I,J)

         ELSE IF ( N == IDTPRPE ) THEN
            EPA_NEI = EPA_WE_BF_PRPE(I,J)

         ELSE IF ( N == IDTC3H8 ) THEN
            EPA_NEI = EPA_WE_BF_C3H8(I,J)

         ELSE IF ( N == IDTCH2O ) THEN
            EPA_NEI = EPA_WE_BF_CH2O(I,J)

         ELSE IF ( N == IDTC2H6 ) THEN
            EPA_NEI = EPA_WE_BF_C2H6(I,J)

         ELSE IF ( N == IDTSO2 ) THEN
            EPA_NEI = EPA_WE_BF_SO2(I,J)

         ELSE IF ( N == IDTSO4 ) THEN
            EPA_NEI = EPA_WE_BF_SO4(I,J)

         ELSE IF ( N == IDTNH3 ) THEN
            EPA_NEI = EPA_WE_BF_NH3(I,J)

         ELSE
            EPA_NEI = 0d0
            
         ENDIF

      ENDIF

      ! Return to calling program
      END FUNCTION GET_EPA_BIOFUEL

!------------------------------------------------------------------------------

      SUBROUTINE INIT_EPA_NEI
!
!******************************************************************************
!  Subroutine INIT_EPA_NEI allocates and zeroes all module arrays.
!  (rch, bmy, 11/10/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LNEI99

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER :: AS
      
      !=================================================================
      ! INIT_EPA_NEI begins here!
      !=================================================================

      ! Return if we LNEI99 = .FALSE.
      IF ( .not. LNEI99 ) RETURN

      ! USA Mask
      ALLOCATE( USA_MASK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'USA_MASK' )
      USA_MASK = 0d0

      !-----------------------
      ! Anthro - weekday avg
      !-----------------------
      ALLOCATE( EPA_WD_AN_NOX( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_NOX' )
      EPA_WD_AN_NOX = 0e0

      ALLOCATE( EPA_WD_AN_CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_CO' )
      EPA_WD_AN_CO = 0e0

      ALLOCATE( EPA_WD_AN_ALK4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_ALK4' )
      EPA_WD_AN_ALK4 = 0e0

      ALLOCATE( EPA_WD_AN_ACET( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_ACET' )
      EPA_WD_AN_ACET = 0e0

      ALLOCATE( EPA_WD_AN_MEK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_MEK' )
      EPA_WD_AN_MEK = 0e0

      ALLOCATE( EPA_WD_AN_PRPE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_PRPE' )
      EPA_WD_AN_PRPE = 0e0

      ALLOCATE( EPA_WD_AN_C2H6( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_C2H6' )
      EPA_WD_AN_C2H6 = 0e0

      ALLOCATE( EPA_WD_AN_C3H8( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_C3H8' )
      EPA_WD_AN_C3H8 = 0e0

      ALLOCATE( EPA_WD_AN_CH2O( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_CH2O' )
      EPA_WD_AN_CH2O = 0e0

      ALLOCATE( EPA_WD_AN_NH3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_NH3' )
      EPA_WD_AN_NH3 = 0e0

      ALLOCATE( EPA_WD_AN_SO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_SO2' )
      EPA_WD_AN_SO2 = 0e0

      ALLOCATE( EPA_WD_AN_SO4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_AN_SO4' )
      EPA_WD_AN_SO4 = 0e0

      !-----------------------
      ! Anthro - weekend avg
      !-----------------------
      ALLOCATE( EPA_WE_AN_NOX( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_NOX' )
      EPA_WE_AN_NOX = 0e0

      ALLOCATE( EPA_WE_AN_CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_CO' )
      EPA_WE_AN_CO = 0e0

      ALLOCATE( EPA_WE_AN_ALK4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_ALK4' )
      EPA_WE_AN_ALK4 = 0e0

      ALLOCATE( EPA_WE_AN_ACET( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_ACET' )
      EPA_WE_AN_ACET = 0e0

      ALLOCATE( EPA_WE_AN_MEK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_MEK' )
      EPA_WE_AN_MEK = 0e0

      ALLOCATE( EPA_WE_AN_PRPE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_PRPE' )
      EPA_WE_AN_PRPE = 0e0

      ALLOCATE( EPA_WE_AN_C2H6( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_C2H6' )
      EPA_WE_AN_C2H6 = 0e0

      ALLOCATE( EPA_WE_AN_C3H8( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_C3H8' )
      EPA_WE_AN_C3H8 = 0e0

      ALLOCATE( EPA_WE_AN_CH2O( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_CH2O' )
      EPA_WE_AN_CH2O = 0e0

      ALLOCATE( EPA_WE_AN_NH3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_NH3' )
      EPA_WE_AN_NH3 = 0e0

      ALLOCATE( EPA_WE_AN_SO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_SO2' )
      EPA_WE_AN_SO2 = 0e0

      ALLOCATE( EPA_WE_AN_SO4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_AN_SO4' )
      EPA_WE_AN_SO4 = 0e0

      !-----------------------
      ! Biofuel - weekday avg
      !-----------------------
      ALLOCATE( EPA_WD_BF_NOX(  IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_NOX' )
      EPA_WD_BF_NOX = 0e0

      ALLOCATE( EPA_WD_BF_CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_CO' )
      EPA_WD_BF_CO = 0e0

      ALLOCATE( EPA_WD_BF_ALK4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_ALK4' )
      EPA_WD_BF_ALK4 = 0e0

      ALLOCATE( EPA_WD_BF_ACET( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_ACET' )
      EPA_WD_BF_ACET = 0e0

      ALLOCATE( EPA_WD_BF_MEK(  IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_MEK' )
      EPA_WD_BF_MEK = 0e0

      ALLOCATE( EPA_WD_BF_PRPE( IIPAR, JJPAR ), STAT=AS )     
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_PRPE' )
      EPA_WD_BF_PRPE = 0e0

      ALLOCATE( EPA_WD_BF_C2H6( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_C2H6' )
      EPA_WD_BF_C2H6 = 0e0

      ALLOCATE( EPA_WD_BF_C3H8( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_C3H8' )
      EPA_WD_BF_C3H8 = 0e0

      ALLOCATE( EPA_WD_BF_CH2O( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_CH2O' )
      EPA_WD_BF_CH2O = 0e0

      ALLOCATE( EPA_WD_BF_NH3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_NH3' )
      EPA_WD_BF_NH3 = 0e0

      ALLOCATE( EPA_WD_BF_SO2( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_SO2' )
      EPA_WD_BF_SO2 = 0e0

      ALLOCATE( EPA_WD_BF_SO4( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WD_BF_SO4' )
      EPA_WD_BF_SO4 = 0e0

      !-----------------------
      ! Biofuel - weekend avg
      !-----------------------
      ALLOCATE( EPA_WE_BF_NOX(  IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_NOX' )
      EPA_WE_BF_NOX = 0e0

      ALLOCATE( EPA_WE_BF_CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_CO' )
      EPA_WE_BF_CO = 0e0

      ALLOCATE( EPA_WE_BF_ALK4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_ALK4' )
      EPA_WE_BF_ALK4 = 0e0

      ALLOCATE( EPA_WE_BF_ACET( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_ACET' )
      EPA_WE_BF_ACET = 0e0

      ALLOCATE( EPA_WE_BF_MEK(  IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_MEK' )
      EPA_WE_BF_MEK = 0e0

      ALLOCATE( EPA_WE_BF_PRPE( IIPAR, JJPAR ), STAT=AS )     
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_PRPE' )
      EPA_WE_BF_PRPE = 0e0

      ALLOCATE( EPA_WE_BF_C2H6( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_C2H6' )
      EPA_WE_BF_C2H6 = 0e0

      ALLOCATE( EPA_WE_BF_C3H8( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_C3H8' )
      EPA_WE_BF_C3H8 = 0e0

      ALLOCATE( EPA_WE_BF_CH2O( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_CH2O' )
      EPA_WE_BF_CH2O = 0e0

      ALLOCATE( EPA_WE_BF_NH3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_NH3' )
      EPA_WE_BF_NH3 = 0e0

      ALLOCATE( EPA_WE_BF_SO2( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_SO2' )
      EPA_WE_BF_SO2 = 0e0

      ALLOCATE( EPA_WE_BF_SO4( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPA_WE_BF_SO4' )
      EPA_WE_BF_SO4 = 0e0

      ! Return to calling program
      END SUBROUTINE INIT_EPA_NEI

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_EPA_NEI
!
!******************************************************************************
!  Subroutine CLEANUP_EPA_NEI deallocates all module arrays 
!  (rch, bmy, 11/10/04)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_EPA_NEI begins here!
      !=================================================================

      ! USA mask
      IF ( ALLOCATED( USA_MASK       ) ) DEALLOCATE( USA_MASK       )

      ! Fossil fuel -- weekday
      IF ( ALLOCATED( EPA_WD_AN_NOX  ) ) DEALLOCATE( EPA_WD_AN_NOX  )
      IF ( ALLOCATED( EPA_WD_AN_CO   ) ) DEALLOCATE( EPA_WD_AN_CO   )
      IF ( ALLOCATED( EPA_WD_AN_ALK4 ) ) DEALLOCATE( EPA_WD_AN_ALK4 )
      IF ( ALLOCATED( EPA_WD_AN_ACET ) ) DEALLOCATE( EPA_WD_AN_ACET )
      IF ( ALLOCATED( EPA_WD_AN_MEK  ) ) DEALLOCATE( EPA_WD_AN_MEK  )
      IF ( ALLOCATED( EPA_WD_AN_PRPE ) ) DEALLOCATE( EPA_WD_AN_PRPE )
      IF ( ALLOCATED( EPA_WD_AN_C2H6 ) ) DEALLOCATE( EPA_WD_AN_C2H6 )
      IF ( ALLOCATED( EPA_WD_AN_C3H8 ) ) DEALLOCATE( EPA_WD_AN_C3H8 )
      IF ( ALLOCATED( EPA_WD_AN_C2H6 ) ) DEALLOCATE( EPA_WD_AN_C2H6 )
      IF ( ALLOCATED( EPA_WD_AN_NH3  ) ) DEALLOCATE( EPA_WD_AN_NH3  )
      IF ( ALLOCATED( EPA_WD_AN_SO2  ) ) DEALLOCATE( EPA_WD_AN_SO2  )
      IF ( ALLOCATED( EPA_WD_AN_SO4  ) ) DEALLOCATE( EPA_WD_AN_SO4  )

      ! Fossil fuel -- weekend
      IF ( ALLOCATED( EPA_WE_AN_NOX  ) ) DEALLOCATE( EPA_WE_AN_NOX  )
      IF ( ALLOCATED( EPA_WE_AN_CO   ) ) DEALLOCATE( EPA_WE_AN_CO   )
      IF ( ALLOCATED( EPA_WE_AN_ALK4 ) ) DEALLOCATE( EPA_WE_AN_ALK4 )
      IF ( ALLOCATED( EPA_WE_AN_ACET ) ) DEALLOCATE( EPA_WE_AN_ACET )
      IF ( ALLOCATED( EPA_WE_AN_MEK  ) ) DEALLOCATE( EPA_WE_AN_MEK  )
      IF ( ALLOCATED( EPA_WE_AN_PRPE ) ) DEALLOCATE( EPA_WE_AN_PRPE )
      IF ( ALLOCATED( EPA_WE_AN_C2H6 ) ) DEALLOCATE( EPA_WE_AN_C2H6 )
      IF ( ALLOCATED( EPA_WE_AN_C3H8 ) ) DEALLOCATE( EPA_WE_AN_C3H8 )
      IF ( ALLOCATED( EPA_WE_AN_C2H6 ) ) DEALLOCATE( EPA_WE_AN_C2H6 )
      IF ( ALLOCATED( EPA_WE_AN_NH3  ) ) DEALLOCATE( EPA_WE_AN_NH3  )
      IF ( ALLOCATED( EPA_WE_AN_SO2  ) ) DEALLOCATE( EPA_WE_AN_SO2  )
      IF ( ALLOCATED( EPA_WE_AN_SO4  ) ) DEALLOCATE( EPA_WE_AN_SO4  )

      ! Biofuel -- weekday
      IF ( ALLOCATED( EPA_WD_BF_NOX  ) ) DEALLOCATE( EPA_WD_BF_NOX  )
      IF ( ALLOCATED( EPA_WD_BF_CO   ) ) DEALLOCATE( EPA_WD_BF_CO   )
      IF ( ALLOCATED( EPA_WD_BF_ALK4 ) ) DEALLOCATE( EPA_WD_BF_ALK4 )
      IF ( ALLOCATED( EPA_WD_BF_ACET ) ) DEALLOCATE( EPA_WD_BF_ACET )
      IF ( ALLOCATED( EPA_WD_BF_MEK  ) ) DEALLOCATE( EPA_WD_BF_MEK  )
      IF ( ALLOCATED( EPA_WD_BF_PRPE ) ) DEALLOCATE( EPA_WD_BF_PRPE )
      IF ( ALLOCATED( EPA_WD_BF_C2H6 ) ) DEALLOCATE( EPA_WD_BF_C2H6 )
      IF ( ALLOCATED( EPA_WD_BF_C3H8 ) ) DEALLOCATE( EPA_WD_BF_C3H8 )
      IF ( ALLOCATED( EPA_WD_BF_C2H6 ) ) DEALLOCATE( EPA_WD_BF_C2H6 )
      IF ( ALLOCATED( EPA_WD_BF_NH3  ) ) DEALLOCATE( EPA_WD_BF_NH3  )
      IF ( ALLOCATED( EPA_WD_BF_SO2  ) ) DEALLOCATE( EPA_WD_BF_SO2  )
      IF ( ALLOCATED( EPA_WD_BF_SO4  ) ) DEALLOCATE( EPA_WD_BF_SO4  )

      ! Biofuel -- weekend
      IF ( ALLOCATED( EPA_WE_BF_NOX  ) ) DEALLOCATE( EPA_WE_BF_NOX  )
      IF ( ALLOCATED( EPA_WE_BF_CO   ) ) DEALLOCATE( EPA_WE_BF_CO   )
      IF ( ALLOCATED( EPA_WE_BF_ALK4 ) ) DEALLOCATE( EPA_WE_BF_ALK4 )
      IF ( ALLOCATED( EPA_WE_BF_ACET ) ) DEALLOCATE( EPA_WE_BF_ACET )
      IF ( ALLOCATED( EPA_WE_BF_MEK  ) ) DEALLOCATE( EPA_WE_BF_MEK  )
      IF ( ALLOCATED( EPA_WE_BF_PRPE ) ) DEALLOCATE( EPA_WE_BF_PRPE )
      IF ( ALLOCATED( EPA_WE_BF_C2H6 ) ) DEALLOCATE( EPA_WE_BF_C2H6 )
      IF ( ALLOCATED( EPA_WE_BF_C3H8 ) ) DEALLOCATE( EPA_WE_BF_C3H8 )
      IF ( ALLOCATED( EPA_WE_BF_C2H6 ) ) DEALLOCATE( EPA_WE_BF_C2H6 )
      IF ( ALLOCATED( EPA_WE_BF_NH3  ) ) DEALLOCATE( EPA_WE_BF_NH3  )
      IF ( ALLOCATED( EPA_WE_BF_SO2  ) ) DEALLOCATE( EPA_WE_BF_SO2  )
      IF ( ALLOCATED( EPA_WE_BF_SO4  ) ) DEALLOCATE( EPA_WE_BF_SO4  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_EPA_NEI

!------------------------------------------------------------------------------

      END MODULE EPA_NEI_MOD
