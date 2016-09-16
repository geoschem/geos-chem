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
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !DEFINED PARAMETERS:
!
      ! Maximum number of output:
      INTEGER, PUBLIC, PARAMETER   :: PD42 = 20
!
! !PUBLIC DATA MEMBERS:
!
      INTEGER, PUBLIC              :: ND42            ! ND42 on/off flag
      INTEGER, PUBLIC              :: LD42            ! # of levels for ND42

      ! SOAupdate: ratio of OM/OC (hotp 6/12/09)
      REAL(fp), PARAMETER    :: OCFPOA  = 1.4e+0_fp
      REAL(fp), PARAMETER    :: OCFOPOA = 1.4e+0_fp*1.5e+0_fp  ! 2.1

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
!  10 Nov 2014 - M. Yannetti - Added PRECISION_MOD
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
      USE CARBON_MOD,         ONLY : BETANOSAVE ! NO branching ratio diagnostic
      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! NDxx flags
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants,      ONLY : ATM
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Met_Mod,      ONLY : MetState
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
!  26 Feb 2015 - E. Lundgren - Remove dependency on pressure_mod (not used)
!  25 Mar 2015 - E. Lundgren - Change tracer units from kg to kg/kg
!  06 Jan 2016 - E. Lundgren - Use global physical parameter ATM
!  16 Jun 2016 - K. Yu       - Now define species ID's with the Ind_ function
!  17 Jun 2016 - R. Yantosca - Now only define species ID's on the first call
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  11 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  01 Sep 2016 - M. Sulprizio- Add PM2.5 concentrations to index 18
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I,       J,   L
      ! SOAupdate: Temporary variables for diagnostics (hotp 8/24/09)
      REAL(fp)            :: FACTOR, PRES
      REAL(fp)            :: TEMP6, TEMP7

      ! Species ID flags
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER, SAVE       :: id_OCPI,  id_OCPO
      INTEGER, SAVE       :: id_SOAG,  id_SOAM
      INTEGER, SAVE       :: id_SO4,   id_NIT,   id_NH4
      INTEGER, SAVE       :: id_SALA,  id_SALC
      INTEGER, SAVE       :: id_BCPI,  id_BCPO
      INTEGER, SAVE       :: id_POA1,  id_POA2
      INTEGER, SAVE       :: id_OPOA1, id_OPOA2
      INTEGER, SAVE       :: id_ASOAN, id_ASOA1
      INTEGER, SAVE       :: id_ASOA2, id_ASOA3
      INTEGER, SAVE       :: id_TSOA1, id_TSOA2
      INTEGER, SAVE       :: id_TSOA3, id_TSOA0
      INTEGER, SAVE       :: id_ISOA1, id_ISOA2, id_ISOA3
      INTEGER, SAVE       :: id_DST1,  id_DST2

      ! SOAupdate: Logical SOA tracer flags (mpayer, 7/12/11)
      LOGICAL             :: IS_SOA1,  IS_SOA2,  IS_SOA3,  IS_SOA4
      LOGICAL             :: IS_SOA5,  IS_SOA1to5
      LOGICAL             :: IS_OC,    IS_BC,    IS_SO4,   IS_NH4
      LOGICAL             :: IS_NIT,   IS_SAL,   IS_SOAG,  IS_SOAM
      LOGICAL             :: IS_TSOA,  IS_ISOA,  IS_ASOA
      LOGICAL             :: IS_POA,   IS_OPOA,  IS_DST
!
! !DEFINED PARAMETERS:
!
      ! Factor for computing standard volume
      REAL(fp), PARAMETER :: STD_VOL_FAC = ATM * 1e-2_fp / 273.15e+0_fp
     
      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)

      !================================================================= 
      ! DIAG42 begins here! 
      !================================================================= 

      ! First-time setup
      IF ( FIRST ) THEN

         ! Define species ID flags
         id_OCPI  = Ind_('OCPI' )
         id_OCPO  = Ind_('OCPO' )
         id_SOAG  = Ind_('SOAG' )
         id_SOAM  = Ind_('SOAM' )
         id_SO4   = Ind_('SO4'  )
         id_NIT   = Ind_('NIT'  )
         id_NH4   = Ind_('NH4'  )
         id_DST1  = Ind_('DST1' )
         id_DST2  = Ind_('DST2' )
         id_SALA  = Ind_('SALA' )
         id_SALC  = Ind_('SALC' )
         id_BCPI  = Ind_('BCPI' )
         id_BCPO  = Ind_('BCPO' )
         id_POA1  = Ind_('POA1' )
         id_POA2  = Ind_('POA2' )
         id_OPOA1 = Ind_('OPOA1')
         id_OPOA2 = Ind_('OPOA2')
         id_ASOAN = Ind_('ASOAN')
         id_ASOA1 = Ind_('ASOA1')
         id_ASOA2 = Ind_('ASOA2')
         id_ASOA3 = Ind_('ASOA3')
         id_TSOA1 = Ind_('TSOA1')
         id_TSOA2 = Ind_('TSOA2')
         id_TSOA3 = Ind_('TSOA3')
         id_TSOA0 = Ind_('TSOA0')
         id_ISOA1 = Ind_('ISOA1')
         id_ISOA2 = Ind_('ISOA2') 
         id_ISOA3 = Ind_('ISOA3')

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! SOAupdate: Define logical flags to decide whether or not to archive
      ! into AD42 array.  This will prevent out-of-bounds errors. 
      ! (mpayer, 7/12/11)
      IS_OC      = ( id_OCPI  > 0 .AND. id_OCPO  > 0 )
      IS_BC      = ( id_BCPI  > 0 .AND. id_BCPO  > 0 )
      IS_SO4     = ( id_SO4   > 0 )
      IS_NH4     = ( id_NH4   > 0 )
      IS_NIT     = ( id_NIT   > 0 )
      IS_DST     = ( id_DST1  > 0 .AND. id_DST2  > 0 )
      IS_SAL     = ( id_SALA  > 0 .AND. id_SALC  > 0 )
      IS_SOAG    = ( id_SOAG  > 0 )
      IS_SOAM    = ( id_SOAM  > 0 )
      IS_TSOA    = ( id_TSOA1 > 0 .AND. id_TSOA2 > 0 
     &         .AND. id_TSOA3 > 0 .AND. id_TSOA0 > 0 )
      IS_ISOA    = ( id_ISOA1 > 0 .AND. id_ISOA2 > 0 
     &         .AND. id_ISOA3 > 0)
      IS_ASOA    = ( id_ASOAN > 0 .AND. id_ASOA1 > 0 
     &         .AND. id_ASOA2 > 0 .AND. id_ASOA3 > 0 )
      IS_POA     = ( id_POA1  > 0 .AND. id_POA2  > 0 )
      IS_OPOA    = ( id_OPOA1 > 0 .AND. id_OPOA2 > 0 )

      ! Point to the 
      Spc => State_Chm%Species

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
         ! (removed by ewl, 3/25/15)
         ! FACTOR        = 1e+9_fp / State_Met%AIRVOL(I,J,L)

         ! Conversion factor from [kg/kg] --> [ug/m3] 
         ! (ewl, 3/25/15)
         FACTOR        = 1e+9_fp * State_Met%AIRDEN(I,J,L)

         !--------------------------------------------------------------
         ! TSOA [ug/m3]
         ! (terpene SOA)
         !--------------------------------------------------------------
         IF ( IS_TSOA ) THEN
            AD42(I,J,L,1) = AD42(I,J,L,1)       + 
     &                  ( ( Spc(I,J,L,id_TSOA1) +
     &                      Spc(I,J,L,id_TSOA2) +
     &                      Spc(I,J,L,id_TSOA3) +
     &                      Spc(I,J,L,id_TSOA0) ) * FACTOR )
         ENDIF

         !--------------------------------------------------------------
         ! ISOA [ug/m3]
         ! (isoprene SOA)
         !--------------------------------------------------------------
         IF ( IS_ISOA ) THEN
            AD42(I,J,L,2) = AD42(I,J,L,2)       + 
     &                  ( ( Spc(I,J,L,id_ISOA1) +
     &                      Spc(I,J,L,id_ISOA2) +
     &                      Spc(I,J,L,id_ISOA3) ) * FACTOR )
         ENDIF

         !--------------------------------------------------------------
         ! ASOA [ug/m3]
         ! (benz, tolu, xyle, + NAP/IVOC SOA)
         !--------------------------------------------------------------
         IF ( IS_ASOA ) THEN
            AD42(I,J,L,3) = AD42(I,J,L,3)       + 
     &                  ( ( Spc(I,J,L,id_ASOAN) +
     &                      Spc(I,J,L,id_ASOA1) +
     &                      Spc(I,J,L,id_ASOA2) +
     &                      Spc(I,J,L,id_ASOA3) ) * FACTOR )
         ENDIF

         !--------------------------------------------------------------
         ! POA [ug/m3]
         !--------------------------------------------------------------
         IF ( IS_POA ) THEN
            AD42(I,J,L,4) = AD42(I,J,L,4)      + 
     &                  ( ( Spc(I,J,L,id_POA1) +
     &                      Spc(I,J,L,id_POA2) ) * OCFPOA  * FACTOR )
         ELSEIF ( IS_OC ) THEN
            AD42(I,J,L,4) = AD42(I,J,L,4)      + 
     &                  ( ( Spc(I,J,L,id_OCPI) +
     &                      Spc(I,J,L,id_OCPO) ) * OCFOPOA * FACTOR )
         ENDIF

         !--------------------------------------------------------------
         ! OPOA [ug/m3]
         !--------------------------------------------------------------
         IF ( IS_OPOA ) THEN
            AD42(I,J,L,5) = AD42(I,J,L,5)       + 
     &                  ( ( Spc(I,J,L,id_OPOA1) +
     &                      Spc(I,J,L,id_OPOA2) ) * OCFOPOA 
     &                    * FACTOR )
         ENDIF

         IF ( IS_TSOA .AND. IS_ISOA .AND. IS_ASOA ) THEN
            !-----------------------------------------------------------
            ! Sum of all organic aerosol [ug/m3]
            !-----------------------------------------------------------
            TEMP6 = Spc(I,J,L,id_TSOA1) +
     &              Spc(I,J,L,id_TSOA2) + 
     &              Spc(I,J,L,id_TSOA3) + 
     &              Spc(I,J,L,id_TSOA0) + 
     &              Spc(I,J,L,id_ISOA1) + 
     &              Spc(I,J,L,id_ISOA2) + 
     &              Spc(I,J,L,id_ISOA3) + 
     &              Spc(I,J,L,id_ASOAN) + 
     &              Spc(I,J,L,id_ASOA1) + 
     &              Spc(I,J,L,id_ASOA2) + 
     &              Spc(I,J,L,id_ASOA3) 
 
            IF ( IS_POA ) THEN
               TEMP6 = TEMP6 + Spc(I,J,L,id_POA1) * OCFPOA +
     &                         Spc(I,J,L,id_POA2) * OCFPOA
            ENDIF
          
            IF ( IS_OPOA ) THEN
               TEMP6 = TEMP6 + Spc(I,J,L,id_OPOA1) * OCFOPOA +
     &                         Spc(I,J,L,id_OPOA2) * OCFOPOA
            ENDIF

            IF ( IS_OC ) THEN
               TEMP6 = TEMP6 + Spc(I,J,L,id_OCPI) * OCFOPOA +
     &                         Spc(I,J,L,id_OCPO) * OCFOPOA
            ENDIF

            AD42(I,J,L,6) = AD42(I,J,L,6)  + ( TEMP6 * FACTOR )

            !-----------------------------------------------------------
            ! Sum of all organic carbon [ugC/m3]
            !-----------------------------------------------------------
            TEMP7 = (  Spc(I,J,L,id_TSOA1) +
     &                 Spc(I,J,L,id_TSOA2) + 
     &                 Spc(I,J,L,id_TSOA3) + 
     &                 Spc(I,J,L,id_TSOA0) + 
     &                 Spc(I,J,L,id_ISOA1) + 
     &                 Spc(I,J,L,id_ISOA2) + 
     &                 Spc(I,J,L,id_ISOA3) + 
     &                 Spc(I,J,L,id_ASOAN) + 
     &                 Spc(I,J,L,id_ASOA1) + 
     &                 Spc(I,J,L,id_ASOA2) + 
     &                 Spc(I,J,L,id_ASOA3)  ) / 2.1e+0_fp

            IF ( IS_POA ) THEN
               TEMP7 = TEMP7 + Spc(I,J,L,id_POA1) +
     &                         Spc(I,J,L,id_POA2) 
            ENDIF
          
            IF ( IS_OPOA ) THEN
               TEMP7 = TEMP7 + Spc(I,J,L,id_OPOA1) +
     &                         Spc(I,J,L,id_OPOA2)
            ENDIF

            IF ( IS_OC ) THEN
               TEMP7 = TEMP7 + Spc(I,J,L,id_OCPI) +
     &                         Spc(I,J,L,id_OCPO)
            ENDIF

            AD42(I,J,L,7) = AD42(I,J,L,7)  +
     &                    ( TEMP7 * FACTOR )

         ENDIF

         !--------------------------------------------------------------
         ! Sum of biogenic aerosol [ug/m3]
         !--------------------------------------------------------------
         IF ( IS_TSOA .AND. IS_ISOA ) THEN
            AD42(I,J,L,8) = AD42(I,J,L,8)       +
     &                  ( ( Spc(I,J,L,id_TSOA1) +
     &                      Spc(I,J,L,id_TSOA2) +
     &                      Spc(I,J,L,id_TSOA3) +
     &                      Spc(I,J,L,id_TSOA0) +
     &                      Spc(I,J,L,id_ISOA1) +
     &                      Spc(I,J,L,id_ISOA2) +
     &                      Spc(I,J,L,id_ISOA3) ) * FACTOR )
         ENDIF
          
         !--------------------------------------------------------------
         ! NO branching ratio [dimless]
         !--------------------------------------------------------------
         ! will have zero or junk values if not in troposphere
         AD42(I,J,L,9) = AD42(I,J,L,9) + BETANOSAVE(I,J,L)

         !--------------------------------------------------------------
         ! POA [ugC/m3]
         !--------------------------------------------------------------
         IF ( IS_POA ) THEN
            AD42(I,J,L,10) = AD42(I,J,L,10)     + 
     &                   ( ( Spc(I,J,L,id_POA1) +
     &                       Spc(I,J,L,id_POA2) )  * FACTOR )
         ELSEIF ( IS_OC ) THEN
            AD42(I,J,L,10) = AD42(I,J,L,10)     + 
     &                   ( ( Spc(I,J,L,id_OCPI) +
     &                       Spc(I,J,L,id_OCPO) )  * FACTOR )
         ENDIF

         !--------------------------------------------------------------
         ! OPOA [ugC/m3]
         !--------------------------------------------------------------
         IF ( IS_OPOA ) THEN
            AD42(I,J,L,11) = AD42(I,J,L,11)      + 
     &                   ( ( Spc(I,J,L,id_OPOA1) +
     &                       Spc(I,J,L,id_OPOA2) ) * FACTOR )
         ENDIF

         !--------------------------------------------------------------
         ! Additional aerosol tracers
         !--------------------------------------------------------------

         ! OC [ugC/m3]
         IF ( IS_OC ) THEN
            AD42(I,J,L,12) = AD42(I,J,L,12)       + 
     &                       ( Spc(I,J,L,id_OCPI) + 
     &                         Spc(I,J,L,id_OCPO) ) * FACTOR
         ENDIF

         ! BC [ugC/m3]
         IF ( IS_BC ) THEN
            AD42(I,J,L,13) = AD42(I,J,L,13)       + 
     &                       ( Spc(I,J,L,id_BCPI) +
     &                         Spc(I,J,L,id_BCPO) ) * FACTOR 
         ENDIF

         ! SO4 [ug/m3]
         IF ( IS_SO4 ) THEN
            AD42(I,J,L,14) = AD42(I,J,L,14)      + 
     &                       ( Spc(I,J,L,id_SO4) * FACTOR )
         ENDIF

         ! NH4 [ug/m3]
         IF ( IS_NH4 ) THEN
            AD42(I,J,L,15) = AD42(I,J,L,15)      + 
     &                       ( Spc(I,J,L,id_NH4) * FACTOR )
         ENDIF

         ! NIT [ug/m3]
         IF ( IS_NIT ) THEN
            AD42(I,J,L,16) = AD42(I,J,L,16)      + 
     &                       ( Spc(I,J,L,id_NIT) * FACTOR )
         ENDIF

         ! SAL [ug/m3]
         IF ( IS_SAL ) THEN
            AD42(I,J,L,17) = AD42(I,J,L,17)       + 
     &                       ( Spc(I,J,L,id_SALA) +
     &                        Spc(I,J,L,id_SALC) ) * FACTOR 
         ENDIF

         !--------------------------------------------------------------
         ! PM2.5 [ug/m3]
         !
         ! Randall Martin wrote on 6/29/16:
         !
         ! The PM2.5 calculation should account for aerosol water of
         ! PM2.5 in a way that is consistent with the PM2.5 measurements.
         ! The water content of PM2.5 is operationally defined. In the
         ! US an RH of 35% is often used and would be appropriate for
         ! the benchmark. In Europe an RH of 50% is often used.
         !
         ! Below are the factors that account for aerosol water in a
         ! consistent manner as used in jv_spec.dat, as calculated by
         ! Aaron van Donkelaar and Sajeev Philip.
         !   At 35% RH:
         !    1.33 for SO4, NIT, and NH4
         !    1.12 for OCPI and SOA
         !    1.86 for SALA
         !
         ! The DST2 bin includes aerosols with diameter both smaller
         ! and larger than 2.5 um. [The Aerosols WG decided] that 38% of
         ! the DST2 bin should be included in the calculation of PM2.5.
         !
         ! In summary, the PM2.5 benchmark should be
         !
         !   PM25 = 1.33 (NH4 + NIT  + SO4) + BCPI + BCPO +
         !          2.1  (OCPO + 1.12 OCPI) + 1.12 SOA +
         !          DST1 + 0.38 DST2 + 1.86 SALA
         !
         ! where the tracer units are in ug/m3 at STP, and the value of
         ! 2.1 is the global mean OM/OC as recommended by the Aerosol WG. 
         !--------------------------------------------------------------
         IF ( IS_OC   .AND. IS_BC   .AND. IS_SO4  .AND. IS_NH4  .AND.
     &        IS_NIT  .AND. IS_DST  .AND. IS_SAL  .AND. IS_TSOA .AND.
     &        IS_ISOA .AND. IS_ASOA ) THEN
            AD42(I,J,L,18) = AD42(I,J,L,18)                           + 
     &                     ( 1.33e+0_fp*( Spc(I,J,L,id_NH4  )         + 
     &                                    Spc(I,J,L,id_NIT  )         +
     &                                    Spc(I,J,L,id_SO4  ) )       +
     &                                    Spc(I,J,L,id_BCPI )         + 
     &                                    Spc(I,J,L,id_BCPO )         +
     &                       2.10e+0_fp*( Spc(I,J,L,id_OCPO )         +
     &                       1.12e+0_fp * Spc(I,J,L,id_OCPI ) )       +
     &                       1.12e+0_fp*( Spc(I,J,L,id_TSOA0)         +
     &                                    Spc(I,J,L,id_TSOA1)         +
     &                                    Spc(I,J,L,id_TSOA2)         +
     &                                    Spc(I,J,L,id_TSOA3)         +
     &                                    Spc(I,J,L,id_ISOA1)         + 
     &                                    Spc(I,J,L,id_ISOA2)         + 
     &                                    Spc(I,J,L,id_ISOA3)         + 
     &                                    Spc(I,J,L,id_ASOAN)         + 
     &                                    Spc(I,J,L,id_ASOA1)         + 
     &                                    Spc(I,J,L,id_ASOA2)         + 
     &                                    Spc(I,J,L,id_ASOA3) )       +
     &                                    Spc(I,J,L,id_DST1)          +
     &                       0.38e+0_fp * Spc(I,J,L,id_DST2)          +
     &                       1.86e+0_fp * Spc(I,J,L,id_SALA) ) * FACTOR 
         ENDIF

         !--------------------------------------------------------------
         ! Additional diagnostics for SOAG, SOAM (tmf, 12/8/07) 
         !
         ! Assume SOAG mass = GLYX mass, SOAM mass = MGLY mass
         !--------------------------------------------------------------
         IF ( IS_SOAG .AND. IS_SOAM ) THEN

            ! SOAG [ug total mass /m3]
            AD42(I,J,L,19) = AD42(I,J,L,19)        + 
     &                      ( Spc(I,J,L,id_SOAG) * 1.e+0_fp * FACTOR )


            ! SOAM [ug total mass /m3]
            AD42(I,J,L,20) = AD42(I,J,L,20)        + 
     &                      ( Spc(I,J,L,id_SOAM) * 1.e+0_fp * FACTOR )

         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


      ! Free pointer
      Spc => NULL()

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
      USE CMN_DIAG_MOD             ! TINDEX
      USE CMN_SIZE_MOD             ! Size parameters
      USE FILE_MOD,           ONLY : IU_BPCH
      USE GRID_MOD,           ONLY : GET_XOFFSET
      USE GRID_MOD,           ONLY : GET_YOFFSET
      USE Input_Opt_Mod,      ONLY : OptInput
      USE TIME_MOD,           ONLY : GET_CT_DIAG
      USE TIME_MOD,           ONLY : GET_DIAGb
      USE TIME_MOD,           ONLY : GET_DIAGe
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
! 
! !REVISION HISTORY: 
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  (1 ) Replace TINY(1d0) with 1d-32 to avoid problems  on SUN 4100 platform
!        (bmy, 9/5/06)
!  (2 ) Use TS_DIAG for scaling instead of TS_DYN. (ccc, 8/18/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
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
      !REAL(fp)         :: SCALE(IIPAR,JJPAR)
      REAL(fp)          :: SCALE
      REAL(f8)          :: DIAGb,     DIAGe
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
            CASE( 7, 10, 11, 12, 13 )
               UNIT = 'ug C/m3'
            CASE( 9 )
               UNIT = 'dimless'
            CASE DEFAULT
               UNIT = 'ug/m3'
         END SELECT

         ! Apply scale factor
         DO L = 1, LD42
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
      SUBROUTINE INIT_DIAG42( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD   
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE Input_Opt_Mod,      ONLY : OptInput

!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  22 May 2006 - D. Henze, R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  23 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! INIT_DIAG42 begins here!
      !=================================================================

      ! Turn off ND42 if SOA tracers are not used
      IF ( .not. Input_Opt%LSOA ) THEN
         ND42 = 0
         RETURN
      ENDIF

      ! Exit if ND42 is turned off
      IF ( ND42 == 0 ) RETURN

      ! Number of levels to save for this diagnostic
      LD42 = MIN( ND42, LLPAR )

      ! 2-D array ("LFLASH-$")
      ALLOCATE( AD42( IIPAR, JJPAR, LD42, PD42 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'AD42' )

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
