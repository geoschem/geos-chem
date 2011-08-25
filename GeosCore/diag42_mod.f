! $Id: diag42_mod.f,v 1.3 2010/03/15 19:33:24 ccarouge Exp $
      MODULE DIAG42_MOD
!
!******************************************************************************
!  Module DIAG42_MOD contains arrays and routines for archiving the ND42
!  diagnostic -- secondary organic aerosols [ug/m3]. (dkh,bmy,5/22/06,3/29/07)
!
!  Module Variables:
!  ============================================================================
!  (1 ) AD42 (REAL*4)  : Array for SOA concentrations [ug/m3]
!
!  Module Routines:
!  ============================================================================
!  (1 ) DIAG42         : Archives quantities for diagnostic
!  (2 ) ZERO_DIAG42    : Sets all module arrays to zero
!  (3 ) WRITE_DIAG42   : Writes data in module arrays to bpch file
!  (4 ) INIT_DIAG42    : Allocates all module arrays
!  (5 ) CLEANUP_DIAG42 : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by diag03_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f    : Module w/ routines for binary pch file I/O
!  (2 ) error_mod.f    : Module w/ NaN and other error check routines
!  (3 ) file_mod.f     : Module w/ file unit numbers and error checks
!  (4 ) grid_mod.f     : Module w/ horizontal grid information
!  (5 ) pressure_mod.f : Module w/ routines to compute P(I,J,L)
!  (6 ) time_mod.f     : Module w/ routines to compute date & time
!
!  References:
!  ============================================================================
!  (1 ) Pye, H.O.T., and J.H. Seinfeld, "A global perspective on aerosol from
!        low-volatility organic compounds", Atmos. Chem. & Phys., Vol 10, pp
!        4377-4401, 2010.
!
!  NOTES:
!  (1 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (2 ) Now use ratio of 2.1 instead of 1.4 for SOA4 (dkh, bmy, 3/29/07)
!  (3 ) Add diagnostics for SOAG and SOAM (tmf, 1/7/09)
!  (4 ) Increase PD42 to 24. (fp, hotp, 2/3/10)
!  08 Jul 2011 - M. Payer    - Add modifications for SOA + semivol POA (H. Pye)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag42_mod.f"
      !=================================================================

      ! Make everything PUBLIC
      PUBLIC 

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Scalars
      INTEGER              :: ND42, LD42

      ! Parameters
      ! Maximum number of output:
      ! SOA1, SOA2, SOA3, SOA4, SOA5, SUM(SOA1-3), SUM(SOA1-4), SUM(SOA1-5),
      ! SUM(SOA1-5+OC), SUM(SOA1-5+OC), SUM(SOA1-5+OC), OC, BC, SOA4, NH4, NIT,
      ! SSALT, SUM(aerosols), SOAG, SOAM, SUM(SOA1-5+SOAG+SOAM),
      ! SUM(SOA1-5+SOAG+SOAM+OC), SUM(SOA1-5+SOAG+SOAM), 
      ! SUM(SOA1-5+SOAG+SOAM+OC)
      !INTEGER, PARAMETER   :: PD42 = 14
      INTEGER, PARAMETER   :: PD42 = 24
      ! OM/OC Ratios for POA & OPOA (hotp, mpayer, 7/13/11)
      ! (see pp 4382 & 4386, Pye & Seinfeld, 2010)
      REAL*8, PARAMETER    :: OCFPOA  = 1.4d0
      REAL*8, PARAMETER    :: OCFOPOA = 1.4d0*1.5d0  ! 2.1

      ! Arrays
      REAL*4,  ALLOCATABLE :: AD42(:,:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
     
!------------------------------------------------------------------------------

      SUBROUTINE DIAG42
!
!******************************************************************************
!  Subroutine DIAG42 archives SOA concentrations [ug/m3] for the ND42
!  diagnostic. (dkh, bmy, 5/22/06, 3/29/07)
!
!  NOTES:
!  (1 ) Now use ratio of 2.1 instead of 1.4 for SOA4 (dkh, bmy, 3/29/07)
!  12 Jul 2011 - M. Payer    - Add modifications for SOA + semivol POA (H. Pye)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AIRVOL, T
      !USE DIAG_MOD,     ONLY : LTOTH
      USE PRESSURE_MOD, ONLY : GET_PCENTER
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTSOA1, IDTSOA2, IDTSOA3, IDTSOA4
      USE TRACERID_MOD, ONLY : IDTSOA5
      USE TRACERID_MOD, ONLY : IDTOCPI, IDTOCPO
      USE TRACERID_MOD, ONLY : IDTSOAG, IDTSOAM
      USE TRACERID_MOD, ONLY : IDTSO4,  IDTNIT, IDTNH4, IDTSALA, IDTSALC
      USE TRACERID_MOD, ONLY : IDTBCPI, IDTBCPO
      ! For SOA + semivolatile POA (hotp, mpayer, 7/12/11)
      USE TRACERID_MOD, ONLY : IDTPOA1,  IDTPOA2
      USE TRACERID_MOD, ONLY : IDTOPOA1, IDTOPOA2
      USE TRACERID_MOD, ONLY : IDTASOAN, IDTASOA1, IDTASOA2, IDTASOA3
      USE TRACERID_MOD, ONLY : IDTTSOA1, IDTTSOA2, IDTTSOA3, IDTTSOA0
      USE TRACERID_MOD, ONLY : IDTISOA1, IDTISOA2, IDTISOA3
      USE CARBON_MOD,   ONLY : BETANOSAVE
      USE LOGICAL_MOD,  ONLY : LSVPOA

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! NDxx flags

      ! Local variables
      INTEGER               :: I,      J,    L
      REAL*8                :: FACTOR, PRES
      REAL*8                :: TEMP6, TEMP7 ! (hotp, mpayer, 7/12/11)

      ! Factor for computing standard volume
      REAL*8, PARAMETER     :: STD_VOL_FAC = 1013.25d0 / 273.15d0

      ! Logical SOA tracer flags (mpayer, 7/12/11)
      LOGICAL               :: IS_SOA1,  IS_SOA2,  IS_SOA3,  IS_SOA4
      LOGICAL               :: IS_SOA5,  IS_SOA1to5
      LOGICAL               :: IS_OC,    IS_BC,    IS_SO4,   IS_NH4
      LOGICAL               :: IS_NIT,   IS_SAL,   IS_SOAG,  IS_SOAM
      LOGICAL               :: IS_TSOA,  IS_ISOA,  IS_ASOA
      LOGICAL               :: IS_POA,   IS_OPOA 
     
      !================================================================= 
      ! DIAG42 begins here! 
      !================================================================= 

!-----------------------------------------------------------------------
! Prior to 7/12/11:
! Removed for SOA + semivol POA, now use logical flags (mpayer, 7/12/11)
!      ! Error check
!!      IF ( IDTSOA1 == 0 ) RETURN
!!      IF ( IDTSOA2 == 0 ) RETURN
!!      IF ( IDTSOA3 == 0 ) RETURN
!!      IF ( IDTSOA4 == 0 ) RETURN
!!      IF ( IDTSOA5 == 0 ) RETURN
!!      IF ( IDTOCPO == 0 ) RETURN
!!      IF ( IDTOCPI == 0 ) RETURN
!-----------------------------------------------------------------------
        
      ! Define logical flags to decide whether or not to archive
      ! into AD42 array.  This will prevent out-of-bounds errors. 
      ! (mpayer, 7/12/11)
      IS_SOA1    = ( IDTSOA1  > 0 )
      IS_SOA2    = ( IDTSOA2  > 0 )
      IS_SOA3    = ( IDTSOA3  > 0 )
      IS_SOA4    = ( IDTSOA4  > 0 )
      IS_SOA5    = ( IDTSOA5  > 0 )
      IS_SOA1to5 = ( IDTSOA1  > 0 .AND. IDTSOA2  > 0 .AND. IDTSOA3  > 0
     &         .AND. IDTSOA4  > 0 .AND. IDTSOA5  > 0 )
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
      IS_ISOA   = ( IDTISOA1 > 0 .AND. IDTISOA2 > 0 .AND. IDTISOA3 > 0 )
      IS_ASOA    = ( IDTASOAN > 0 .AND. IDTASOA1 > 0 .AND. IDTASOA2 > 0 
     &         .AND. IDTASOA3 > 0 )
      IS_POA     = ( IDTPOA1  > 0 .AND. IDTPOA2  > 0 )
      IS_OPOA    = ( IDTOPOA1 > 0 .AND. IDTOPOA2 > 0 )

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
         !FACTOR        = 1d9 / AIRVOL(I,J,L) * LTOTH(I,J) 

         ! Conversion factor from [kg] --> [ug/m3]
         FACTOR        = 1d9 / AIRVOL(I,J,L)

         !--------------------------------------------------------------
         ! Traditional SOA tracers
         !--------------------------------------------------------------
         ! Add logical flags to all tracers (mpayer, 7/12/11)

         ! SOA1 [ug/m3]
         IF ( IS_SOA1 ) THEN 
            AD42(I,J,L,1) = AD42(I,J,L,1)        + 
     &                      ( STT(I,J,L,IDTSOA1) * FACTOR )
         ENDIF
 
         ! SOA2 [ug/m3]
         IF ( IS_SOA2 ) THEN
            AD42(I,J,L,2) = AD42(I,J,L,2)        + 
     &                      ( STT(I,J,L,IDTSOA2) * FACTOR )
         ENDIF

         ! SOA3 [ug/m3]
         IF ( IS_SOA3 ) THEN
            AD42(I,J,L,3) = AD42(I,J,L,3)        + 
     &                      ( STT(I,J,L,IDTSOA3) * FACTOR )
         ENDIF

         ! SOA4 [ug/m3]
         IF ( IS_SOA4 ) THEN
            AD42(I,J,L,4) = AD42(I,J,L,4)        + 
     &                      ( STT(I,J,L,IDTSOA4) * FACTOR )
         ENDIF

         ! SOA5 [ug/m3]
         IF ( IS_SOA5 ) THEN
            AD42(I,J,L,5) = AD42(I,J,L,5)        + 
     &                      ( STT(I,J,L,IDTSOA5) * FACTOR )
         ENDIF

         ! Sum of original 3 SOA types [ug/m3]
         IF ( IS_SOA1 .AND. IS_SOA2 .AND. IS_SOA3 ) THEN
            AD42(I,J,L,6) = AD42(I,J,L,6)        + 
     &                      ( STT(I,J,L,IDTSOA1) + 
     &                        STT(I,J,L,IDTSOA2) +  
     &                        STT(I,J,L,IDTSOA3) ) * FACTOR
         ENDIF

         ! Sum of all biogenic SOA [ug/m3] 
         IF ( IS_SOA1 .AND. IS_SOA2 .AND. IS_SOA3 .AND. IS_SOA4 ) THEN
            AD42(I,J,L,7) = AD42(I,J,L,7)        + 
     &                      ( STT(I,J,L,IDTSOA1) + 
     &                        STT(I,J,L,IDTSOA2) + 
     &                        STT(I,J,L,IDTSOA3) + 
     &                        STT(I,J,L,IDTSOA4) ) * FACTOR
         ENDIF

         ! Sum of all SOA [ug/m3]
         IF ( IS_SOA1to5 ) THEN
            AD42(I,J,L,8) = AD42(I,J,L,8)        +
     &                      ( STT(I,J,L,IDTSOA1) + 
     &                        STT(I,J,L,IDTSOA2) + 
     &                        STT(I,J,L,IDTSOA3) + 
     &                        STT(I,J,L,IDTSOA4) + 
     &                        STT(I,J,L,IDTSOA5) ) * FACTOR
         ENDIF


         ! Sum of primary OC + SOA1 to SOA4 [ug C/m3] 
         ! Use higher ratio (2.1) of molecular weight of
         ! organic mass per carbon mass accounting for non-carbon
         ! components attached to OC [Turpin and Lim, 2001] 
         IF ( IS_SOA1to5 .AND. IS_OC ) THEN
            AD42(I,J,L,9) = AD42(I,J,L,9)          +
     &                      ( ( STT(I,J,L,IDTSOA1) + 
     &                          STT(I,J,L,IDTSOA2) + 
     &                          STT(I,J,L,IDTSOA3) + 
     &                          STT(I,J,L,IDTSOA4) +
     &                          STT(I,J,L,IDTSOA5))   / 2.1d0
     &                      + ( STT(I,J,L,IDTOCPO) + 
     &                          STT(I,J,L,IDTOCPI) ) ) * FACTOR

            ! Sum of PRIMARY OC + SOA1 to SOA4 [ug C/sm3] at STP
            PRES           = GET_PCENTER( I, J, L )
            AD42(I,J,L,10) = AD42(I,J,L,9) * STD_VOL_FAC * T(I,J,L) 
     &                       / PRES

            ! Sum of all OA in ug/m3
            AD42(I,J,L,11) = AD42(I,J,L,11)       +
     &                       ( STT(I,J,L,IDTSOA1) + 
     &                         STT(I,J,L,IDTSOA2) + 
     &                         STT(I,J,L,IDTSOA3) + 
     &                         STT(I,J,L,IDTSOA4) + 
     &                         STT(I,J,L,IDTSOA5) +
     &                       ( STT(I,J,L,IDTOCPO) +
     &                         STT(I,J,L,IDTOCPI) ) * 2.1d0 )
     &                       * FACTOR
         ENDIF

         !--------------------------------------------------------------
         ! additional aerosol tracers (hotp 10/26/07)
         !--------------------------------------------------------------
         ! OC [ugC/m3]
         IF ( IS_OC ) THEN
            AD42(I,J,L,12) = AD42(I,J,L,12)       + 
     &                       ( STT(I,J,L,IDTOCPI) + 
     &                         STT(I,J,L,IDTOCPO) ) * FACTOR
         ENDIF

         ! Add option for POA (hotp, mpayer, 7/12/11)
         IF ( IS_POA ) THEN
            AD42(I,J,L,12) = AD42(I,J,L,12)       + 
     &                     ( ( STT(I,J,L,IDTPOA1) +
     &                         STT(I,J,L,IDTPOA2) ) * FACTOR )
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

         ! total aerosol [ug/m3]
         IF ( IS_SOA1to5 .AND. IS_SO4 .AND. IS_NH4 .AND. IS_NIT .AND. 
     &        IS_BC      .AND. IS_OC ) THEN
            AD42(I,J,L,18) = AD42(I,J,L,18)        +
     &                       ( STT(I,J,L,IDTSOA1)  +
     &                        STT(I,J,L,IDTSOA2)   +
     &                        STT(I,J,L,IDTSOA3)   +
     &                        STT(I,J,L,IDTSOA4)   +
     &                        STT(I,J,L,IDTSOA5)   +
     &                        STT(I,J,L,IDTSO4)    +
     &                        STT(I,J,L,IDTNH4)    +
     &                        STT(I,J,L,IDTNIT)    +
     &                        STT(I,J,L,IDTBCPI)   +
     &                        STT(I,J,L,IDTBCPO)   +
     &                       ( STT(I,J,L,IDTOCPO)  +
     &                         STT(I,J,L,IDTOCPI) ) * 2.1 )
     &                       * FACTOR
         ENDIF

         !--------------------------------------------------------------
         ! Additional diagnostics for SOAG, SOAM (tmf, 12/8/07) 
         ! Assume SOAG mass = GLYX mass, SOAM mass = MGLY mass
         ! Test if SOAG and SOAM are simulated (ccc, 12/18/08)
         !--------------------------------------------------------------
!-----------------------------------------------------------------------
! Prior to 7/12/11:
! Now use logical flags (mpayer, 7/12/11)
!         IF ( IDTSOAG /= 0 .AND. IDTSOAM /=0 ) THEN
!-----------------------------------------------------------------------
         IF ( IS_SOAG .AND. IS_SOAM ) THEN

            ! SOAG [ug total mass /m3]
            AD42(I,J,L,19) = AD42(I,J,L,19)        + 
     &                      ( STT(I,J,L,IDTSOAG) * 1.d0 * FACTOR )


            ! SOAM [ug total mass /m3]
            AD42(I,J,L,20) = AD42(I,J,L,20)        + 
     &                      ( STT(I,J,L,IDTSOAM) * 1.d0 * FACTOR )


            ! Sum of SOA1 to SOA4, SOAG, SOAM (tmf, 1/31/07)
            IF ( IS_SOA1to5 ) THEN
               AD42(I,J,L,21) = AD42(I,J,L,21)       + 
     &                          ( STT(I,J,L,IDTSOA1) + 
     &                            STT(I,J,L,IDTSOA2) + 
     &                            STT(I,J,L,IDTSOA3) + 
     &                            STT(I,J,L,IDTSOA4) +
     &                            STT(I,J,L,IDTSOA5) +
     &                          ( STT(I,J,L,IDTSOAG) * 1.d0 ) +
     &                          ( STT(I,J,L,IDTSOAM) * 1.d0 )) * FACTOR
            ENDIF


            ! Sum of SOA1 to SOA4, SOAG, SOAM in carbon (tmf, 1/31/07) 
            ! Except SOAG is 0.41 carbon, SOAM is 0.5 carbon
            IF ( IS_SOA1to5 .AND. IS_OC ) THEN
               AD42(I,J,L,22) = AD42(I,J,L,22)         +
     &                          ( ( STT(I,J,L,IDTSOA1) + 
     &                              STT(I,J,L,IDTSOA2) + 
     &                              STT(I,J,L,IDTSOA3) + 
     &                              STT(I,J,L,IDTSOA4) + 
     &                              STT(I,J,L,IDTSOA5) )   / 2.1d0 +
     &                            ( STT(I,J,L,IDTSOAG) * 0.41D0 ) +
     &                            ( STT(I,J,L,IDTSOAM) * 0.50D0 ) +
     &                            ( STT(I,J,L,IDTOCPO) + 
     &                              STT(I,J,L,IDTOCPI) ) ) * FACTOR 
            ENDIF
          

            ! Sum of SOA1 to SOA4, SOAG, SOAM at STP [ug/sm3 STP] (tmf, 1/31/07)  
            IF ( IS_SOA1to5 ) THEN            
               PRES           = GET_PCENTER( I, J, L )
               AD42(I,J,L,23) = AD42(I,J,L,21) * STD_VOL_FAC * T(I,J,L) 
     &                          / PRES
            ENDIF

            ! Sum of all OC [ug C/sm3] at STP (including SOAG, SOAM)
            IF ( IS_SOA1to5 .AND. IS_OC ) THEN
            AD42(I,J,L,24) = AD42(I,J,L,22) * STD_VOL_FAC * T(I,J,L) 
     &                       / PRES
            ENDIF

         ENDIF

         !--------------------------------------------------------------
         ! SOA + semivolatile POA tracers (hotp, mpayer, 7/12/11)
         !--------------------------------------------------------------
         IF ( LSVPOA ) THEN

            ! TSOA (terpene SOA) [ug/m3]
            IF ( IS_TSOA ) THEN
               AD42(I,J,L,25) = AD42(I,J,L,25) + 
     &                    ( ( STT(I,J,L,IDTTSOA1) +
     &                        STT(I,J,L,IDTTSOA2) +
     &                        STT(I,J,L,IDTTSOA3) +
     &                        STT(I,J,L,IDTTSOA0)  ) * FACTOR )
            ENDIF

            ! ISOA (isoprene SOA) [ug/m3]
            IF ( IS_ISOA ) THEN
               AD42(I,J,L,26) = AD42(I,J,L,26) + 
     &                    ( ( STT(I,J,L,IDTISOA1) +
     &                        STT(I,J,L,IDTISOA2) +
     &                        STT(I,J,L,IDTISOA3)  ) * FACTOR )
            ENDIF

            ! ASOA (aromatic SOA: BENZ, TOLU, XYLE, NAP/IVOC) [ug/m3]
            IF ( IS_ASOA ) THEN
               AD42(I,J,L,27) = AD42(I,J,L,27) + 
     &                    ( ( STT(I,J,L,IDTASOAN) +
     &                        STT(I,J,L,IDTASOA1) +
     &                        STT(I,J,L,IDTASOA2) +
     &                        STT(I,J,L,IDTASOA3)  ) * FACTOR )
            ENDIF

            ! POA or OC [ug/m3]
            IF ( IS_POA ) THEN
               AD42(I,J,L,28) = AD42(I,J,L,28) + 
     &                    ( ( STT(I,J,L,IDTPOA1) +
     &                        STT(I,J,L,IDTPOA2)  ) * OCFPOA  * FACTOR )
            ENDIF

            IF ( IS_OC ) THEN
               AD42(I,J,L,28) = AD42(I,J,L,28) + 
     &                    ( ( STT(I,J,L,IDTOCPI) +
     &                        STT(I,J,L,IDTOCPO)  ) * OCFOPOA * FACTOR )
            ENDIF

            ! OPOA [ug/m3]
            IF ( IS_OPOA ) THEN
               AD42(I,J,L,29) = AD42(I,J,L,29) + 
     &                   ( ( STT(I,J,L,IDTOPOA1) +
     &                       STT(I,J,L,IDTOPOA2)  ) * OCFOPOA * FACTOR )
            ENDIF

            ! sum of all OA [ug/m3]
            IF ( IS_TSOA .AND. IS_ISOA .AND. IS_ASOA ) THEN
               TEMP6 = STT(I,J,L,IDTTSOA1) +
     &                 STT(I,J,L,IDTTSOA2) + 
     &                 STT(I,J,L,IDTTSOA3) + 
     &                 STT(I,J,L,IDTTSOA0) + 
     &                 STT(I,J,L,IDTISOA1) + 
     &                 STT(I,J,L,IDTISOA2) + 
     &                 STT(I,J,L,IDTISOA3) + 
     &                 STT(I,J,L,IDTASOAN) + 
     &                 STT(I,J,L,IDTASOA1) + 
     &                 STT(I,J,L,IDTASOA2) + 
     &                 STT(I,J,L,IDTASOA3) 
 
               IF ( IS_POA ) THEN
                  TEMP6 = TEMP6 + STT(I,J,L,IDTPOA1) * OCFPOA +
     &                            STT(I,J,L,IDTPOA2) * OCFPOA
               ENDIF
          
               IF ( IS_OPOA ) THEN
                  TEMP6 = TEMP6 + STT(I,J,L,IDTOPOA1) * OCFOPOA +
     &                            STT(I,J,L,IDTOPOA2) * OCFOPOA
               ENDIF

               IF ( IS_OC ) THEN
                  TEMP6 = TEMP6 + STT(I,J,L,IDTOCPI) * OCFOPOA +
     &                            STT(I,J,L,IDTOCPO) * OCFOPOA
               ENDIF

               AD42(I,J,L,30) = AD42(I,J,L,30)  + ( TEMP6 * FACTOR )

               ! sum of all OC [ugC/m3]
               TEMP7 = (  STT(I,J,L,IDTTSOA1) +
     &                    STT(I,J,L,IDTTSOA2) + 
     &                    STT(I,J,L,IDTTSOA3) + 
     &                    STT(I,J,L,IDTTSOA0) + 
     &                    STT(I,J,L,IDTISOA1) + 
     &                    STT(I,J,L,IDTISOA2) + 
     &                    STT(I,J,L,IDTISOA3) + 
     &                    STT(I,J,L,IDTASOAN) + 
     &                    STT(I,J,L,IDTASOA1) + 
     &                    STT(I,J,L,IDTASOA2) + 
     &                    STT(I,J,L,IDTASOA3)  ) / 2.1d0
               IF ( IS_POA ) THEN
                  TEMP7 = TEMP7 + STT(I,J,L,IDTPOA1) +
     &                            STT(I,J,L,IDTPOA2) 
               ENDIF
          
               IF ( IS_OPOA ) THEN
                  TEMP7 = TEMP7 + STT(I,J,L,IDTOPOA1) +
     &                            STT(I,J,L,IDTOPOA2)
               ENDIF

               IF ( IS_OC ) THEN
                  TEMP7 = TEMP7 + STT(I,J,L,IDTOCPI) +
     &                            STT(I,J,L,IDTOCPO)
               ENDIF

               AD42(I,J,L,31) = AD42(I,J,L,31)  +
     &                    ( TEMP7 * FACTOR )

            ENDIF

            ! sum of biogenic OC [ug/m3]
            IF ( IS_TSOA .AND. IS_ISOA ) THEN
               AD42(I,J,L,32) = AD42(I,J,L,32)  +
     &                    ( ( STT(I,J,L,IDTTSOA1) +
     &                        STT(I,J,L,IDTTSOA2) +
     &                        STT(I,J,L,IDTTSOA3) +
     &                        STT(I,J,L,IDTTSOA0) +
     &                        STT(I,J,L,IDTISOA1) +
     &                        STT(I,J,L,IDTISOA2) +
     &                        STT(I,J,L,IDTISOA3)  ) * FACTOR )
            ENDIF
          
            ! NO branching ratio
            ! will have zero or junk values if not in troposphere
            AD42(I,J,L,33) = AD42(I,J,L,33) + BETANOSAVE(I,J,L)


            ! OPOA [ug C/m3]
            IF ( IS_OPOA ) THEN
               AD42(I,J,L,34) = AD42(I,J,L,34) + 
     &                    ( ( STT(I,J,L,IDTOPOA1) +
     &                        STT(I,J,L,IDTOPOA2)  ) * FACTOR )
            ENDIF

         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


      ! Return to calling program
      END SUBROUTINE DIAG42

!------------------------------------------------------------------------------

      SUBROUTINE ZERO_DIAG42
!
!******************************************************************************
!  Subroutine ZERO_DIAG42 zeroes the ND03 diagnostic arrays. 
!  (dkh, bmy, 5/22/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! ZERO_DIAG42 begins here!
      !=================================================================

      ! Exit if ND42 is turned off
      IF ( ND42 == 0 ) RETURN

      ! Zero arrays
      AD42(:,:,:,:) = 0e0

      ! Return to calling program
      END SUBROUTINE ZERO_DIAG42

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_DIAG42
!
!******************************************************************************
!  Subroutine WRITE_DIAG03 writes the ND03 diagnostic arrays to the binary
!  punch file at the proper time. (bmy, 5/22/06, 9/5/06)
!
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
!  NOTES:
!  (1 ) Replace TINY(1d0) with 1d-32 to avoid problems  on SUN 4100 platform
!        (bmy, 9/5/06)
!  (2 ) Use TS_DIAG for scaling instead of TS_DYN. (ccc, 8/18/09)
!  12 Jul 2011 - M. Payer    - Add modifications for SOA + semivol POA (H. Pye)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,    ONLY : BPCH2, GET_MODELNAME, GET_HALFPOLAR
      !USE DIAG_MOD,     ONLY : CTOTH
      USE FILE_MOD,     ONLY : IU_BPCH
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,     ONLY : GET_CT_DIAG,  GET_DIAGb,  GET_DIAGe

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! TINDEX

      ! Local variables
      INTEGER               :: CENTER180, HALFPOLAR
      INTEGER               :: L,         M,         N
      INTEGER               :: IFIRST,    JFIRST,    LFIRST        
      REAL*4                :: LONRES,    LATRES
      REAL*4                :: ARRAY(IIPAR,JJPAR,LLPAR)
      !REAL*8                :: SCALE(IIPAR,JJPAR)
      REAL*8                :: SCALE
      REAL*8                :: DIAGb,     DIAGe
      CHARACTER(LEN=20)     :: MODELNAME 
      CHARACTER(LEN=40)     :: CATEGORY
      CHARACTER(LEN=40)     :: RESERVED
      CHARACTER(LEN=40)     :: UNIT

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
!-----------------------------------------------------------------------
! Prior to 7/12/11:
! Update for SOA + semivolatile POA (hotp, mpayer, 7/12/11)
!         SELECT CASE ( N )
!            CASE( 10, 24 )
!               UNIT = 'ug C/sm3'
!            CASE( 9, 12, 13, 22 )
!               UNIT = 'ug C/m3'
!            CASE( 23 )
!               UNIT = 'ug/sm3'
!            CASE DEFAULT
!               UNIT = 'ug/m3'
!         END SELECT
!-----------------------------------------------------------------------
         SELECT CASE ( N )
            CASE( 10, 24 )
               UNIT = 'ug C/sm3'
            CASE( 9, 12, 13, 22, 31, 34 )
               UNIT = 'ug C/m3'
            CASE( 23 )
               UNIT = 'ug/sm3'
            CASE( 33 )
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

      ! Return to calling program
      END SUBROUTINE WRITE_DIAG42

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG42
!
!******************************************************************************
!  Subroutine INIT_DIAG42 allocates all module arrays (bmy, 5/22/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LSOA

#     include "CMN_SIZE"    ! Size parameters 

      ! Local variables
      INTEGER              :: AS
      
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

      ! Return to calling program
      END SUBROUTINE INIT_DIAG42

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG42
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG42 deallocates all module arrays (bmy, 5/22/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG42 begins here!
      !=================================================================
      IF ( ALLOCATED( AD42 ) ) DEALLOCATE( AD42 ) 

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG42

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG42_MOD
