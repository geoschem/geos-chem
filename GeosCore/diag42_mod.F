!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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

      ! SOAupdate: ratio of OM/OC (hotp 6/12/09)
      REAL*8, PARAMETER    :: OCFPOA  = 1.4d0
      REAL*8, PARAMETER    :: OCFOPOA = 1.4d0*1.5d0  ! 2.1

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
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
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
! !IROUTINE: diag42
!
! !DESCRIPTION: Subroutine DIAG42 archives SOA concentrations [ug/m3] 
!  for the ND42 diagnostic.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG42( Input_Opt, State_Met, State_Chm )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TRACERID_MOD,       ONLY : IDTOCPI, IDTOCPO
      USE TRACERID_MOD,       ONLY : IDTSOAG, IDTSOAM
      USE TRACERID_MOD,       ONLY : IDTSO4,  IDTNIT, IDTNH4
      USE TRACERID_MOD,       ONLY : IDTSALA, IDTSALC
      USE TRACERID_MOD,       ONLY : IDTBCPI, IDTBCPO
      ! SOAupdate: consider additional species (hotp 10/26/07) MERGE1
      ! semivolpoa: replace OCPO and OCPI with POA (hotp 2/17/09)
      USE TRACERID_MOD,       ONLY : IDTPOA1,  IDTPOA2
      ! semivolpoa4: add OPOA (hotp 3/27/09)
      USE TRACERID_MOD,       ONLY : IDTOPOA1, IDTOPOA2
      ! add lumped arom/ivoc (hotp 5/15/10)
      USE TRACERID_MOD,       ONLY : IDTASOAN, IDTASOA1
      USE TRACERID_MOD,       ONLY : IDTASOA2, IDTASOA3
      ! monoterpene + sesquiterpene SOA (hotp 5/24/10)
      USE TRACERID_MOD,       ONLY : IDTTSOA1, IDTTSOA2
      USE TRACERID_MOD,       ONLY : IDTTSOA3, IDTTSOA0
      ! isoprene SOA (hotp 5/24/10)
      USE TRACERID_MOD,       ONLY : IDTISOA1, IDTISOA2, IDTISOA3
      ! NO branching ratio diagnostic (hotp 5/24/10)
      USE CARBON_MOD,         ONLY : BETANOSAVE

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! NDxx flags
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I,      J,    L
      REAL*8            :: FACTOR, PRES
      ! SOAupdate: Temporary variables for diagnostics (hotp 8/24/09)
      REAL*8            :: TEMP6, TEMP7

      ! SOAupdate: Logical SOA tracer flags (mpayer, 7/12/11)
      LOGICAL           :: IS_SOA1,  IS_SOA2,  IS_SOA3,  IS_SOA4
      LOGICAL           :: IS_SOA5,  IS_SOA1to5
      LOGICAL           :: IS_OC,    IS_BC,    IS_SO4,   IS_NH4
      LOGICAL           :: IS_NIT,   IS_SAL,   IS_SOAG,  IS_SOAM
      LOGICAL           :: IS_TSOA,  IS_ISOA,  IS_ASOA
      LOGICAL           :: IS_POA,   IS_OPOA 
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

      ! SOAupdate: Define logical flags to decide whether or not to archive
      ! into AD42 array.  This will prevent out-of-bounds errors. 
      ! (mpayer, 7/12/11)
      IS_OC      = ( IDTOCPI  > 0 .AND. IDTOCPO  > 0 )
      IS_BC      = ( IDTBCPI  > 0 .AND. IDTBCPO  > 0 )
      IS_SO4     = ( IDTSO4   > 0 )
      IS_NH4     = ( IDTNH4   > 0 )
      IS_NIT     = ( IDTNIT   > 0 )
      IS_SAL     = ( IDTSALA  > 0 .AND. IDTSALC  > 0 )
      IS_SOAG    = ( IDTSOAG  > 0 )
      IS_SOAM    = ( IDTSOAM  > 0 )
      IS_TSOA    = ( IDTTSOA1 > 0 .AND. IDTTSOA2 > 0 .AND. IDTTSOA3 > 0
     &         .AND. IDTTSOA0 > 0 )
      IS_ISOA    = ( IDTISOA1 > 0 .AND. IDTISOA2 > 0 .AND. IDTISOA3 > 0)
      IS_ASOA    = ( IDTASOAN > 0 .AND. IDTASOA1 > 0 .AND. IDTASOA2 > 0 
     &         .AND. IDTASOA3 > 0 )
      IS_POA     = ( IDTPOA1  > 0 .AND. IDTPOA2  > 0 )
      IS_OPOA    = ( IDTOPOA1 > 0 .AND. IDTOPOA2 > 0 )

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! SOAupdate: TEMPS need to be private!!!(hotp 10/14/09)
      ! Loop over grid boxes     
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, FACTOR, PRES )  
!$OMP+PRIVATE( TEMP6, TEMP7 )
      DO L = 1, LD42  
      DO J = 1, JJPAR 
      DO I = 1, IIPAR

         ! Conversion factor from [kg] --> [ug/m3]
         ! (LTOTH=1 if between OTH_HR1 and OTH_HR2, LTOTH=0 otherwise)
         !FACTOR        = 1d9 / State_Met%AIRVOL(I,J,L) * LTOTH(I,J) 

         ! Conversion factor from [kg] --> [ug/m3]
         FACTOR        = 1d9 / State_Met%AIRVOL(I,J,L)

         ! SOAupdate: AD42(1-11) updated (hotp 7/29/10)
         ! TSOA (terpene SOA) (hotp 5/24/10)
         IF ( IS_TSOA ) THEN
            AD42(I,J,L,1) = AD42(I,J,L,1)       + 
     &                  ( ( STT(I,J,L,IDTTSOA1) +
     &                      STT(I,J,L,IDTTSOA2) +
     &                      STT(I,J,L,IDTTSOA3) +
     &                      STT(I,J,L,IDTTSOA0) ) * FACTOR )
         ENDIF

         ! ISOA (isoprene SOA) (hotp 5/24/10)
         IF ( IS_ISOA ) THEN
            AD42(I,J,L,2) = AD42(I,J,L,2)       + 
     &                  ( ( STT(I,J,L,IDTISOA1) +
     &                      STT(I,J,L,IDTISOA2) +
     &                      STT(I,J,L,IDTISOA3) ) * FACTOR )
         ENDIF

         ! ASOA (benz, tolu, xyle, + NAP/IVOC SOA) (hotp 5/24/10)
         IF ( IS_ASOA ) THEN
            AD42(I,J,L,3) = AD42(I,J,L,3)       + 
     &                  ( ( STT(I,J,L,IDTASOAN) +
     &                      STT(I,J,L,IDTASOA1) +
     &                      STT(I,J,L,IDTASOA2) +
     &                      STT(I,J,L,IDTASOA3) ) * FACTOR )
         ENDIF

         ! POA hotp
         IF ( IS_POA ) THEN
            AD42(I,J,L,4) = AD42(I,J,L,4)      + 
     &                  ( ( STT(I,J,L,IDTPOA1) +
     &                      STT(I,J,L,IDTPOA2) ) * OCFPOA  * FACTOR )
         ELSEIF ( IS_OC ) THEN
            AD42(I,J,L,4) = AD42(I,J,L,4)      + 
     &                  ( ( STT(I,J,L,IDTOCPI) +
     &                      STT(I,J,L,IDTOCPO) ) * OCFOPOA * FACTOR )
         ENDIF

         ! OPOA hotp
         IF ( IS_OPOA ) THEN
            AD42(I,J,L,5) = AD42(I,J,L,5)       + 
     &                  ( ( STT(I,J,L,IDTOPOA1) +
     &                      STT(I,J,L,IDTOPOA2) ) * OCFOPOA 
     &                    * FACTOR )
         ENDIF

         ! sum of all OA (hotp 5/24/10)
         IF ( IS_TSOA .AND. IS_ISOA .AND. IS_ASOA ) THEN
            TEMP6 = STT(I,J,L,IDTTSOA1) +
     &              STT(I,J,L,IDTTSOA2) + 
     &              STT(I,J,L,IDTTSOA3) + 
     &              STT(I,J,L,IDTTSOA0) + 
     &              STT(I,J,L,IDTISOA1) + 
     &              STT(I,J,L,IDTISOA2) + 
     &              STT(I,J,L,IDTISOA3) + 
     &              STT(I,J,L,IDTASOAN) + 
     &              STT(I,J,L,IDTASOA1) + 
     &              STT(I,J,L,IDTASOA2) + 
     &              STT(I,J,L,IDTASOA3) 
 
            IF ( IS_POA ) THEN
               TEMP6 = TEMP6 + STT(I,J,L,IDTPOA1) * OCFPOA +
     &                         STT(I,J,L,IDTPOA2) * OCFPOA
            ENDIF
          
            IF ( IS_OPOA ) THEN
               TEMP6 = TEMP6 + STT(I,J,L,IDTOPOA1) * OCFOPOA +
     &                         STT(I,J,L,IDTOPOA2) * OCFOPOA
            ENDIF

            IF ( IS_OC ) THEN
               TEMP6 = TEMP6 + STT(I,J,L,IDTOCPI) * OCFOPOA +
     &                         STT(I,J,L,IDTOCPO) * OCFOPOA
            ENDIF

            AD42(I,J,L,6) = AD42(I,J,L,6)  + ( TEMP6 * FACTOR )

            ! sum of all OC (hotp 5/24/10)
            TEMP7 = (  STT(I,J,L,IDTTSOA1) +
     &                 STT(I,J,L,IDTTSOA2) + 
     &                 STT(I,J,L,IDTTSOA3) + 
     &                 STT(I,J,L,IDTTSOA0) + 
     &                 STT(I,J,L,IDTISOA1) + 
     &                 STT(I,J,L,IDTISOA2) + 
     &                 STT(I,J,L,IDTISOA3) + 
     &                 STT(I,J,L,IDTASOAN) + 
     &                 STT(I,J,L,IDTASOA1) + 
     &                 STT(I,J,L,IDTASOA2) + 
     &                 STT(I,J,L,IDTASOA3)  ) / 2.1d0

            IF ( IS_POA ) THEN
               TEMP7 = TEMP7 + STT(I,J,L,IDTPOA1) +
     &                         STT(I,J,L,IDTPOA2) 
            ENDIF
          
            IF ( IS_OPOA ) THEN
               TEMP7 = TEMP7 + STT(I,J,L,IDTOPOA1) +
     &                         STT(I,J,L,IDTOPOA2)
            ENDIF

            IF ( IS_OC ) THEN
               TEMP7 = TEMP7 + STT(I,J,L,IDTOCPI) +
     &                         STT(I,J,L,IDTOCPO)
            ENDIF

            AD42(I,J,L,7) = AD42(I,J,L,7)  +
     &                    ( TEMP7 * FACTOR )

         ENDIF

         ! sum of biogenic OC (hotp 5/24/10)
         IF ( IS_TSOA .AND. IS_ISOA ) THEN
            AD42(I,J,L,8) = AD42(I,J,L,8)       +
     &                  ( ( STT(I,J,L,IDTTSOA1) +
     &                      STT(I,J,L,IDTTSOA2) +
     &                      STT(I,J,L,IDTTSOA3) +
     &                      STT(I,J,L,IDTTSOA0) +
     &                      STT(I,J,L,IDTISOA1) +
     &                      STT(I,J,L,IDTISOA2) +
     &                      STT(I,J,L,IDTISOA3) ) * FACTOR )
         ENDIF
          
         ! NO branching ratio (hotp 5/24/10)
         ! will have zero or junk values if not in troposphere
         AD42(I,J,L,9) = AD42(I,J,L,9) + BETANOSAVE(I,J,L)

         ! POA (ugC/m3) hotp 7/28/10
         IF ( IS_POA ) THEN
            AD42(I,J,L,10) = AD42(I,J,L,10)     + 
     &                   ( ( STT(I,J,L,IDTPOA1) +
     &                       STT(I,J,L,IDTPOA2) )  * FACTOR )
         ELSEIF ( IS_OC ) THEN
            AD42(I,J,L,10) = AD42(I,J,L,10)     + 
     &                   ( ( STT(I,J,L,IDTOCPI) +
     &                       STT(I,J,L,IDTOCPO) )  * FACTOR )
         ENDIF

         ! OPOA (ugC/m3) hotp 7/28/10
         IF ( IS_OPOA ) THEN
            AD42(I,J,L,11) = AD42(I,J,L,11)      + 
     &                   ( ( STT(I,J,L,IDTOPOA1) +
     &                       STT(I,J,L,IDTOPOA2) ) * FACTOR )
         ENDIF

         !--------------------------------------------------------------
         ! Additional aerosol tracers (hotp 10/26/07)
         !--------------------------------------------------------------

         ! OC [ugC/m3]
         IF ( IS_OC ) THEN
            AD42(I,J,L,12) = AD42(I,J,L,12)       + 
     &                       ( STT(I,J,L,IDTOCPI) + 
     &                         STT(I,J,L,IDTOCPO) ) * FACTOR
         ENDIF

         ! BC [ugC/m3]
         IF ( IS_BC ) THEN
            AD42(I,J,L,13) = AD42(I,J,L,13)       + 
     &                       ( STT(I,J,L,IDTBCPI) +
     &                         STT(I,J,L,IDTBCPO) ) * FACTOR 
         ENDIF

         ! SO4 [ug/m3]
         IF ( IS_SO4 ) THEN
            AD42(I,J,L,14) = AD42(I,J,L,14)      + 
     &                       ( STT(I,J,L,IDTSO4) * FACTOR )
         ENDIF

         ! NH4 [ug/m3]
         IF ( IS_NH4 ) THEN
            AD42(I,J,L,15) = AD42(I,J,L,15)      + 
     &                       ( STT(I,J,L,IDTNH4) * FACTOR )
         ENDIF

         ! NIT [ug/m3]
         IF ( IS_NIT ) THEN
            AD42(I,J,L,16) = AD42(I,J,L,16)      + 
     &                       ( STT(I,J,L,IDTNIT) * FACTOR )
         ENDIF

         ! SAL [ug/m3]
         IF ( IS_SAL ) THEN
            AD42(I,J,L,17) = AD42(I,J,L,17)       + 
     &                       ( STT(I,J,L,IDTSALA) +
     &                        STT(I,J,L,IDTSALC) ) * FACTOR 
         ENDIF

         !--------------------------------------------------------------
         ! Additional diagnostics for SOAG, SOAM (tmf, 12/8/07) 
         ! Assume SOAG mass = GLYX mass, SOAM mass = MGLY mass
         ! Test if SOAG and SOAM are simulated (ccc, 12/18/08)
         !--------------------------------------------------------------
         IF ( IS_SOAG .AND. IS_SOAM ) THEN

            ! SOAG [ug total mass /m3]
            AD42(I,J,L,19) = AD42(I,J,L,19)        + 
     &                      ( STT(I,J,L,IDTSOAG) * 1.d0 * FACTOR )


            ! SOAM [ug total mass /m3]
            AD42(I,J,L,20) = AD42(I,J,L,20)        + 
     &                      ( STT(I,J,L,IDTSOAM) * 1.d0 * FACTOR )

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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      SUBROUTINE WRITE_DIAG42( Input_Opt )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : BPCH2
      USE BPCH2_MOD,          ONLY : GET_MODELNAME
      USE BPCH2_MOD,          ONLY : GET_HALFPOLAR
      !USE DIAG_MOD,          ONLY : CTOTH
      USE FILE_MOD,           ONLY : IU_BPCH
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GRID_MOD,           ONLY : GET_XOFFSET
      USE GRID_MOD,           ONLY : GET_YOFFSET
      USE TIME_MOD,           ONLY : GET_CT_DIAG
      USE TIME_MOD,           ONLY : GET_DIAGb
      USE TIME_MOD,           ONLY : GET_DIAGe

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! TINDEX
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
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
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
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
            ! SOAupdate: update units (hotp 5/24/10)
            CASE( 7, 10, 11, 12, 13, 22 )
               UNIT = 'ug C/m3'
            CASE( 9 )
               UNIT = 'dimless'
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
