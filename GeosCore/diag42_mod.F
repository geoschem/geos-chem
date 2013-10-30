!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag42_mod
!
! !DESCRIPTION: Module DIAG42\_MOD contains arrays and routines for archiving 
!  the ND42 diagnostic -- secondary organic aerosols [ug/m3]. 
!\\
!\\
! !INTERFACE: 
!
      MODULE DIAG42_MOD
!
! !USES:
!
      IMPLICIT NONE
#     include "define.h"
      PRIVATE
!
! !DEFINED PARAMETERS:
!
      ! Maximum number of output:
      ! SOA1, SOA2, SOA3, SOA4, SOA5, SUM(SOA1-3), SUM(SOA1-4), SUM(SOA1-5),
      ! SUM(SOA1-5+OC), SUM(SOA1-5+OC), SUM(SOA1-5+OC), OC, BC, SOA4, NH4, NIT,
      ! SSALT, SUM(aerosols), SOAG, SOAM, SUM(SOA1-5+SOAG+SOAM),
      ! SUM(SOA1-5+SOAG+SOAM+OC), SUM(SOA1-5+SOAG+SOAM), 
      ! SUM(SOA1-5+SOAG+SOAM+OC)
      INTEGER, PUBLIC, PARAMETER   :: PD42 = 24
!
! !PUBLIC DATA MEMBERS:
!
      INTEGER, PUBLIC              :: ND42            ! ND42 on/off flag
      INTEGER, PUBLIC              :: LD42            ! # of levels for ND42

      ! Arrays
      REAL*4,  PUBLIC, ALLOCATABLE :: AD42(:,:,:,:)   ! Array for SOA [ug/m3]
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC :: DIAG42
      PUBLIC :: ZERO_DIAG42
      PUBLIC :: WRITE_DIAG42
      PUBLIC :: INIT_DIAG42
      PUBLIC :: CLEANUP_DIAG42
!
! !REVISION HISTORY:
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  (1 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (2 ) Now use ratio of 2.1 instead of 1.4 for SOA4 (dkh, bmy, 3/29/07)
!  (3 ) Add diagnostics for SOAG and SOAM (tmf, 1/7/09)
!  (4 ) Increase PD42 to 24. (fp, hotp, 2/3/10)

!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diag42
!
! !DESCRIPTION: Subroutine DIAG42 archives SOA concentrations [ug/m3] 
!  for the ND42 diagnostic.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG42( State_Met, State_Chm )
!
! !USES:
!
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TRACERID_MOD,       ONLY : IDTSOA1, IDTSOA2, IDTSOA3, IDTSOA4
      USE TRACERID_MOD,       ONLY : IDTSOA5
      USE TRACERID_MOD,       ONLY : IDTOCPI, IDTOCPO
      USE TRACERID_MOD,       ONLY : IDTSOAG, IDTSOAM
      USE TRACERID_MOD,       ONLY : IDTSO4,  IDTNIT, IDTNH4
      USE TRACERID_MOD,       ONLY : IDTSALA, IDTSALC
      USE TRACERID_MOD,       ONLY : IDTBCPI, IDTBCPO

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! NDxx flags
!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
! 
! !REVISION HISTORY: 
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  (1 ) Now use ratio of 2.1 instead of 1.4 for SOA4 (dkh, bmy, 3/29/07)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I,      J,    L
      REAL*8            :: FACTOR, PRES
!
! !DEFINED PARAMETERS:
!
      ! Factor for computing standard volume
      REAL*8, PARAMETER :: STD_VOL_FAC = 1013.25d0 / 273.15d0
     
      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !================================================================= 
      ! DIAG42 begins here! 
      !================================================================= 

      ! Error check
      IF ( IDTSOA1 == 0 ) RETURN
      IF ( IDTSOA2 == 0 ) RETURN
      IF ( IDTSOA3 == 0 ) RETURN
      IF ( IDTSOA4 == 0 ) RETURN
      IF ( IDTSOA5 == 0 ) RETURN
      IF ( IDTOCPO == 0 ) RETURN
      IF ( IDTOCPI == 0 ) RETURN

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Loop over grid boxes     
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, FACTOR, PRES )  
      DO L = 1, LD42  
      DO J = 1, JJPAR 
      DO I = 1, IIPAR

         ! Conversion factor from [kg] --> [ug/m3]
         ! (LTOTH=1 if between OTH_HR1 and OTH_HR2, LTOTH=0 otherwise)
         !FACTOR        = 1d9 / State_Met%AIRVOL(I,J,L) * LTOTH(I,J) 

         ! Conversion factor from [kg] --> [ug/m3]
         FACTOR        = 1d9 / State_Met%AIRVOL(I,J,L)

         ! SOA1 [ug/m3]
         AD42(I,J,L,1) = AD42(I,J,L,1)        + 
     &                   ( STT(I,J,L,IDTSOA1) * FACTOR )
 
         ! SOA2 [ug/m3]
         AD42(I,J,L,2) = AD42(I,J,L,2)        + 
     &                   ( STT(I,J,L,IDTSOA2) * FACTOR )

         ! SOA3 [ug/m3]
         AD42(I,J,L,3) = AD42(I,J,L,3)        + 
     &                   ( STT(I,J,L,IDTSOA3) * FACTOR )

         ! SOA4 [ug/m3]
         AD42(I,J,L,4) = AD42(I,J,L,4)        + 
     &                   ( STT(I,J,L,IDTSOA4) * FACTOR )

         ! SOA5 [ug/m3]
         AD42(I,J,L,5) = AD42(I,J,L,5)        + 
     &                   ( STT(I,J,L,IDTSOA5) * FACTOR )

         ! Sum of original 3 SOA types [ug/m3]
         AD42(I,J,L,6) = AD42(I,J,L,6)        + 
     &                   ( STT(I,J,L,IDTSOA1) + 
     &                     STT(I,J,L,IDTSOA2) +  
     &                     STT(I,J,L,IDTSOA3) ) * FACTOR

         ! Sum of all biogenic SOA [ug/m3] 
         AD42(I,J,L,7) = AD42(I,J,L,7)        + 
     &                   ( STT(I,J,L,IDTSOA1) + 
     &                     STT(I,J,L,IDTSOA2) + 
     &                     STT(I,J,L,IDTSOA3) + 
     &                     STT(I,J,L,IDTSOA4) ) * FACTOR

         ! Sum of all SOA [ug/m3] 
         AD42(I,J,L,8) = AD42(I,J,L,8)        +
     &                   ( STT(I,J,L,IDTSOA1) + 
     &                     STT(I,J,L,IDTSOA2) + 
     &                     STT(I,J,L,IDTSOA3) + 
     &                     STT(I,J,L,IDTSOA4) + 
     &                     STT(I,J,L,IDTSOA5) ) * FACTOR


         ! Sum of primary OC + SOA1 to SOA4 [ug C/m3] 
         ! Use higher ratio (2.1) of molecular weight of
         ! organic mass per carbon mass accounting for non-carbon
         ! components attached to OC [Turpin and Lim, 2001] 
         AD42(I,J,L,9) = AD42(I,J,L,9)          +
     &                   ( ( STT(I,J,L,IDTSOA1) + 
     &                       STT(I,J,L,IDTSOA2) + 
     &                       STT(I,J,L,IDTSOA3) + 
     &                       STT(I,J,L,IDTSOA4) +
     &                       STT(I,J,L,IDTSOA5))   / 2.1d0
     &                   + ( STT(I,J,L,IDTOCPO) + 
     &                       STT(I,J,L,IDTOCPI) ) ) * FACTOR

         ! Sum of PRIMARY OC + SOA1 to SOA4 [ug C/sm3] at STP
         PRES          = GET_PCENTER( I, J, L )
         AD42(I,J,L,10) = AD42(I,J,L,9)      * STD_VOL_FAC * 
     &                    State_Met%T(I,J,L) / PRES

         ! Sum of all OA in ug/m3
         AD42(I,J,L,11) = AD42(I,J,L,11)        +
     &                    ( STT(I,J,L,IDTSOA1) + 
     &                      STT(I,J,L,IDTSOA2) + 
     &                      STT(I,J,L,IDTSOA3) + 
     &                      STT(I,J,L,IDTSOA4) + 
     &                      STT(I,J,L,IDTSOA5) +
     &                    ( STT(I,J,L,IDTOCPO) +
     &                      STT(I,J,L,IDTOCPI) ) * 2.1d0 )
     &                    * FACTOR

!-------------------------------------------
! additional aerosol tracers (hotp 10/26/07)
!-------------------------------------------
         ! OC [ugC/m3]
         AD42(I,J,L,12) = AD42(I,J,L,12)       + 
     &                    ( STT(I,J,L,IDTOCPI) + 
     &                      STT(I,J,L,IDTOCPO) ) * FACTOR

         ! BC [ugC/m3]
         AD42(I,J,L,13) = AD42(I,J,L,13)        + 
     &                   ( STT(I,J,L,IDTBCPI)   +
     &                     STT(I,J,L,IDTBCPO)   ) * FACTOR 

         ! SO4 [ug/m3]
         AD42(I,J,L,14) = AD42(I,J,L,14)        + 
     &                   ( STT(I,J,L,IDTSO4)   * FACTOR )

         ! NH4 [ug/m3]
         AD42(I,J,L,15) = AD42(I,J,L,15)        + 
     &                   ( STT(I,J,L,IDTNH4)   * FACTOR )

         ! NIT [ug/m3]
         AD42(I,J,L,16) = AD42(I,J,L,16)        + 
     &                   ( STT(I,J,L,IDTNIT)   * FACTOR )

         ! SAL [ug/m3]
         AD42(I,J,L,17) = AD42(I,J,L,17)        + 
     &                   ( STT(I,J,L,IDTSALA)   +
     &                     STT(I,J,L,IDTSALC)   ) * FACTOR 

         ! total aerosol [ug/m3]
         AD42(I,J,L,18) = AD42(I,J,L,18)        +
     &                    ( STT(I,J,L,IDTSOA1)  +
     &                     STT(I,J,L,IDTSOA2)   +
     &                     STT(I,J,L,IDTSOA3)   +
     &                     STT(I,J,L,IDTSOA4)   +
     &                     STT(I,J,L,IDTSOA5)   +
     &                     STT(I,J,L,IDTSO4)    +
     &                     STT(I,J,L,IDTNH4)    +
     &                     STT(I,J,L,IDTNIT)    +
     &                     STT(I,J,L,IDTBCPI)   +
     &                     STT(I,J,L,IDTBCPO)   +
     &                    ( STT(I,J,L,IDTOCPO)  +
     &                      STT(I,J,L,IDTOCPI) ) * 2.1 )
     &                    * FACTOR

!--------------------------------------------------------
! Additional diagnostics for SOAG, SOAM (tmf, 12/8/07) 
! Assume SOAG mass = GLYX mass, SOAM mass = MGLY mass
! Test if SOAG and SOAM are simulated (ccc, 12/18/08)
!--------------------------------------------------------
         IF ( IDTSOAG /= 0 .AND. IDTSOAM /=0 ) THEN
            ! SOAG [ug total mass /m3]
            AD42(I,J,L,19) = AD42(I,J,L,19)        + 
     &                      ( STT(I,J,L,IDTSOAG) * 1.d0 * FACTOR )

            ! SOAM [ug total mass /m3]
            AD42(I,J,L,20) = AD42(I,J,L,20)        + 
     &                      ( STT(I,J,L,IDTSOAM) * 1.d0 * FACTOR )


            ! Sum of SOA1 to SOA4, SOAG, SOAM (tmf, 1/31/07)
            AD42(I,J,L,21) = AD42(I,J,L,21)        + 
     &                      ( STT(I,J,L,IDTSOA1) + 
     &                        STT(I,J,L,IDTSOA2) + 
     &                        STT(I,J,L,IDTSOA3) + 
     &                        STT(I,J,L,IDTSOA4)  +
     &                        STT(I,J,L,IDTSOA5)  +
     &                      ( STT(I,J,L,IDTSOAG) * 1.d0 ) +
     &                      ( STT(I,J,L,IDTSOAM) * 1.d0 )) * FACTOR 

            ! Sum of SOA1 to SOA4, SOAG, SOAM in carbon (tmf, 1/31/07) 
            ! Except SOAG is 0.41 carbon, SOAM is 0.5 carbon
            AD42(I,J,L,22) = AD42(I,J,L,22)          +
     &                      ( ( STT(I,J,L,IDTSOA1) + 
     &                          STT(I,J,L,IDTSOA2) + 
     &                          STT(I,J,L,IDTSOA3) + 
     &                          STT(I,J,L,IDTSOA4) + 
     &                          STT(I,J,L,IDTSOA5) )   / 2.1d0 +
     &                         ( STT(I,J,L,IDTSOAG) * 0.41D0 ) +
     &                         ( STT(I,J,L,IDTSOAM) * 0.50D0 ) +
     &                        ( STT(I,J,L,IDTOCPO) + 
     &                          STT(I,J,L,IDTOCPI) ) ) * FACTOR

            ! Sum of SOA1 to SOA4, SOAG, SOAM at STP [ug/sm3 STP] (tmf, 1/31/07)  
            PRES          = GET_PCENTER( I, J, L )
            AD42(I,J,L,23) = AD42(I,J,L,21)     * STD_VOL_FAC *
     &                       State_Met%T(I,J,L) / PRES

            ! Sum of all OC [ug C/sm3] at STP (including SOAG, SOAM)
            AD42(I,J,L,24) = AD42(I,J,L,22)     * STD_VOL_FAC *
     &                       State_Met%T(I,J,L) / PRES
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO 

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE DIAG42
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: zero_diag42
!
! !DESCRIPTION: Subroutine ZERO\_DIAG42 zeroes all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ZERO_DIAG42
! 
! !REVISION HISTORY: 
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! ZERO_DIAG42 begins here!
      !=================================================================

      ! Exit if ND42 is turned off
      IF ( ND42 == 0 ) RETURN

      ! Zero arrays
      AD42(:,:,:,:) = 0e0

      END SUBROUTINE ZERO_DIAG42
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diag42
!
! !DESCRIPTION: Subroutine WRITE\_DIAG42 writes the ND42 diagnostic arrays 
!  to the binary punch file at the proper time.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_DIAG42
!
! !USES:
!
      USE BPCH2_MOD,    ONLY : BPCH2
      USE BPCH2_MOD,    ONLY : GET_MODELNAME
      USE BPCH2_MOD,    ONLY : GET_HALFPOLAR
      !USE DIAG_MOD,     ONLY : CTOTH
      USE FILE_MOD,     ONLY : IU_BPCH
      USE GRID_MOD,     ONLY : GET_XOFFSET
      USE GRID_MOD,     ONLY : GET_YOFFSET
      USE TIME_MOD,     ONLY : GET_CT_DIAG
      USE TIME_MOD,     ONLY : GET_DIAGb
      USE TIME_MOD,     ONLY : GET_DIAGe

      USE CMN_SIZE_MOD     ! Size parameters
      USE CMN_DIAG_MOD     ! TINDEX
!
! !REMARKS:
!   # : Field    : Description                 : Units    : Scale factor
!  -----------------------------------------------------------------------
!  (1 ) IJ-SOA-$ : SOA1                        : ug/m3    : SCALE_OTH
!  (2 ) IJ-SOA-$ : SOA2                        : ug/m3    : SCALE_OTH
!  (3 ) IJ-SOA-$ : SOA3                        : ug/m3    : SCALE_OTH
!  (4 ) IJ-SOA-$ : SOA4                        : ug/m3    : SCALE_OTH
!  (5 ) IJ-SOA-$ : SOA1 + SOA2 + SOA3          : ug/m3    : SCALE_OTH
!  (6 ) IJ-SOA-$ : SOA1 + SOA2 + SOA3 + SOA4   : ug/m3    : SCALE_OTH
!  (7 ) IJ-SOA-$ : Sum of all Org Carbon       : ug C/m3  : SCALE_OTH
!  (8 ) IJ-SOA-$ : Sum of all Org Carbon @ STP : ug C/sm3 : SCALE_OTH
! 
! !REVISION HISTORY: 
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  (1 ) Replace TINY(1d0) with 1d-32 to avoid problems  on SUN 4100 platform
!        (bmy, 9/5/06)
!  (2 ) Use TS_DIAG for scaling instead of TS_DYN. (ccc, 8/18/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: CENTER180, HALFPOLAR
      INTEGER           :: L,         M,         N
      INTEGER           :: IFIRST,    JFIRST,    LFIRST        
      REAL*4            :: LONRES,    LATRES
      REAL*4            :: ARRAY(IIPAR,JJPAR,LLPAR)
      !REAL*8            :: SCALE(IIPAR,JJPAR)
      REAL*8            :: SCALE
      REAL*8            :: DIAGb,     DIAGe
      CHARACTER(LEN=20) :: MODELNAME 
      CHARACTER(LEN=40) :: CATEGORY
      CHARACTER(LEN=40) :: RESERVED
      CHARACTER(LEN=40) :: UNIT

      !=================================================================
      ! WRITE_DIAG42 begins here!
      !=================================================================

      ! Exit if ND03 is turned off
      IF ( ND42 == 0 ) RETURN

      ! Initialize
      CENTER180 = 1
      DIAGb     = GET_DIAGb()
      DIAGe     = GET_DIAGe()
      HALFPOLAR = GET_HALFPOLAR()
      IFIRST    = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST    = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LATRES    = DJSIZE
      LFIRST    = 1
      LONRES    = DISIZE
      MODELNAME = GET_MODELNAME()
      RESERVED  = ''
      SCALE     = DBLE( GET_CT_DIAG() ) + TINY( 1e0 )

      !=================================================================
      ! Write data to the bpch file
      !=================================================================

      ! Loop over ND03 diagnostic tracers
      DO M = 1, TMAX(42)

         ! Define quantities
         N        = TINDEX(42,M)
         CATEGORY = 'IJ-SOA-$'

         ! Pick proper unit
         SELECT CASE ( N )
            CASE( 10, 24 )
               UNIT = 'ug C/sm3'
            CASE( 9, 12, 13, 22 )
               UNIT = 'ug C/m3'
            CASE( 23 )
               UNIT = 'ug/sm3'
            CASE DEFAULT
               UNIT = 'ug/m3'
         END SELECT

         ! Apply scale factor
         DO L = 1, LD42
            !ARRAY(:,:,L) = AD42(:,:,L,N) / SCALE(:,:)
            ARRAY(:,:,L) = AD42(:,:,L,N) / SCALE
         ENDDO

         ! Write data to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD42,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD42) )
      ENDDO

      END SUBROUTINE WRITE_DIAG42
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diag42
!
! !DESCRIPTION: Subroutine INIT\_DIAG42 allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAG42
!
! !USES:
!
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LSOA

      USE CMN_SIZE_MOD    ! Size parameters 
! 
! !REVISION HISTORY: 
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS
      
      !=================================================================
      ! INIT_DIAG42 begins here!
      !=================================================================

      ! Turn off ND42 if SOA tracers are not used
      IF ( .not. LSOA ) THEN
         ND42 = 0
         RETURN
      ENDIF

      ! Exit if ND42 is turned off
      IF ( ND42 == 0 ) RETURN

      ! Number of levels to save for this diagnostic
      LD42 = MIN( ND42, LLPAR )

      ! 2-D array ("LFLASH-$")
      ALLOCATE( AD42( IIPAR, JJPAR, LD42, PD42 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD42' )

      ! Zero arrays
      CALL ZERO_DIAG42

      END SUBROUTINE INIT_DIAG42
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_diag42
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG42 deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_DIAG42
! 
! !REVISION HISTORY: 
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_DIAG42 begins here!
      !=================================================================
      IF ( ALLOCATED( AD42 ) ) DEALLOCATE( AD42 ) 

      END SUBROUTINE CLEANUP_DIAG42
!EOC
      END MODULE DIAG42_MOD
