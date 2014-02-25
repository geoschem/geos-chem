!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: initialize
!
! !DESCRIPTION: Subroutine INITIALIZE does the following:
!
! \begin{enumerate}
! \item Zeroes globally defined GEOS-CHEM variables.
! \item Zeroes accumulating diagnostic arrays.
! \item Resets certain year/month/day and counter variables used 
!       in GEOS-Chem diagnostic subroutines. 
! \end{enumerate}
!
! !INTERFACE:
!
      SUBROUTINE INITIALIZE( IFLAG, am_I_Root )
!
! !USES:
!
      ! Modules from Headers subdirectory
      USE CMN_SIZE_MOD          
      USE CMN_DIAG_MOD
      USE ERROR_MOD
      USE TIME_MOD

      ! Modules from GeosCore subdirectory
      USE DIAG_MOD
      USE DIAG03_MOD
      USE DIAG04_MOD
      USE DIAG41_MOD
      USE DIAG42_MOD
      USE DIAG53_MOD
      USE DIAG56_MOD
      USE DIAG_PL_MOD
      USE LOGICAL_MOD

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      ! If IFLAG=1, zero global CTM arrays 
      ! If IFLAG=2, zero accumulating diagnostic arrays
      ! If IFLAG=3, zero accumulating diagnostic counters
      INTEGER, INTENT(IN)  :: IFLAG

      ! Is this the root CPU?
      LOGICAL, INTENT(IN) :: am_I_Root   
!
! !REMARKS:
!  Eventually we will fold this into "diag_mod.f" in a cleaner,
!   more consistent fashion.  Think about this later (bmy, 11/14/02)
! 
! !REVISION HISTORY: 
!  15 Jun 1998 - M. Prather  - Initial version
!  (1 ) INITIALIZE is written in Fixed-Form Fortran 90.
!  (2 ) To ensure double precision accuracy, use 0d0 instead of 0.0.
!  (3 ) Also zero the mass flux arrays from TPCORE (bmy, 4/26/99)
!  (4 ) Only zero allocatable arrays that are turned on. (bmy, 11/29/99)
!  (5 ) Added arrays for ND13 diagnostic -- sulfur emissions.
!        Also updated comments (bmy, 6/21/00)
!  (6 ) Remove SAVEJ and SAVEL -- we don't call DIAG0 anymore (bmy, 9/8/00)
!  (7 ) Add array AD32_bf for ND32 NOx biofuel diagnostic (bmy, 9/12/00)
!  (8 ) Also zero the FAMPL array for ND65 (bmy, 12/5/00)
!  (9 ) Now initialize AD34 array for biofuel emissions (bmy, 3/15/01)
!  (10) Now initialize AD12 array for boundary layer emissions in "setemis.f".
!        Also made cosmetic changes & updated comments. (bdf, bmy, 6/15/01)
!  (11) Now initialize AD11 array for acetone diagnostic (bmy, 8/1/01)
!  (12) Remove reference to AVGF -- it is obsolete.  Also, AVGW is now
!        included in "dao_mod.f", and is initialized there. (bmy, 9/25/01)
!  (13) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (14) Make sure FAMPL is allocated before we reference it (bmy, 1/15/02)
!  (15) Eliminated obsolete code from 1/02.  Now also zero CTNO2, CTHO2
!        counter arrays. (bmy, 2/27/02)
!  (16) Bug fix: CTHO2 and CTNO2 should be zeroed if ND43 > 0, not if
!        ND45 > 0.  Fix this typo. (bmy, 4/19/02) 
!  (17) Now also zero AD01, AD02 arrays (bmy, 8/7/02)
!  (18) Remove reference to arrays P, SIG, SIGE from "CMN", since we now
!        use floating pressure + the hybrid grid. (dsa, bdf, bmy, 8/21/02)
!  (19) Now zero the AD05 array for sulfate P-L (rjp, bdf, bmy, 9/20/02)
!  (20) Now we no longer have to zero the T array.  Also reference ERROR_STOP
!        from "error_mod.f". Now also initialize AD13_NH3_an, AD13_NH3_bb,
!        AD13_NH3_bf. (bmy, 12/13/02)
!  (21) Now also zero AD13_NH3_na array for ND13 (rjp, bmy, 3/23/03)
!  (22) Now references "time_mod.f" (bmy, 3/27/03)
!  (23) Now zeroes AD03 array for Kr85 prod/loss diag. (jsw, bmy, 8/20/03)
!  (24) Now also zeroes AD06 and AD07* arrays (rjp, tdf, bmy, 4/5/04)
!  (25) Now also zeroes AD08 array (rjp, bec, bmy, 4/20/04)
!  (26) Now also initialize AD13_SO2_sh array (bec, bmy, 5/20/04)
!  (27) Now also initialize AD07_HC array (rjp, bmy, 7/13/04)
!  (28) Now references AD65 & FAM_PL from "diag_pl_mod.f".  Now remove
!        reference to DIAGCHLORO, it's obsolete. (bmy, 7/20/04)
!  (29) Now initialize extra arrays for ND03 mercury diag.  Also remove
!        reference to obsolete TOFDY0 variable. (eck, bmy, 12/7/04)
!  (30) Now initialize AD21_cr array for ND21 diag.  Also references 
!        LCRYST from "logical_mod.f"  Now call ZERO_DIAG03 from "diag03_mod.f"
!        to zero ND03 arrays (bmy, 1/21/05) 
!  (31) Now call ZERO_DIAG41 from "diag41_mod.f".  Also removed references
!        to AD41 and AFTTOT. (bmy, 2/17/05)
!  (32) Now zero AD09 and AD09_em for HCN simulation (xyp, bmy, 6/27/05)
!  (33) Now references ND04, ZERO_DIAG04 from "diag04_mod.f".  Also remove
!        reference to "CMN" and XTRA2.  Now zeroes AD30 array (bmy, 8/18/05)
!  (34) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (35) Now resets SET_CT_XTRA at the beginning of the run. (tmf, 10/20/05)
!  (36) Now references ND56, ZERO_DIAG56 from "diag56_mod.f" (ltm, bmy, 5/5/06)
!  (37) Now references ND42, ZERO_DIAG42 from "diag42_mod.f" (dkh, bmy,5/22/06)
!  (38) take care of AD54 (time in the troposphere diagnostic) (phs, 10/17/06)
!  (39) Now also zero CTO3 array.  Bug fix: ZERO_DIAG42 is now called when
!        ND42 is turned on. (phs, bmy, 1/30/07)
!  (40) Now zero AD10 and AD10em for H2HD simulation (phs, 9/18/07)
!  (41) Now zero CTO3_24h (phs, 11/17/08)
!  (42) Now zero AD52 for Gamma HO2 diag. (ccc, jaegle, 2/26/09)
!  (43) Updated to diagnose GLYX production of SOAG in ND07. (tmf, 1/7/09)
!  (44) Add initialization of counter for diag time steps. (ccc, 7/20/09)
!  (45) Define new diagnostics, ND19, ND58, ND60 for methane 
!       (kjw, 8/18/09)
!  (46) Add ND59 and ND60 for initialization (win, 7/28/09)
!  (47) Add potential temperature diagnostic. (fp, 06/09)
!  (48) Add TOMAS diags using ifdefs. (sfarina, 01/13)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  25 Aug 2010 - R. Yantosca - Now also reset the counter for A1 timesteps
!  08 Feb 2012 - R. Yantosca - Rewrote USE statements, for clarity
!  08 Feb 2012 - R. Yantosca - Now also reset the counter for I3 timesteps
!  15 Oct 2012 - R. Yantosca - Bug fix, make sure Counter arrays CTLBRO2H etc. 
!                              are allocated before we use them
!  02 Apr 2013 - M. Payer    - Remove code for CTNO, CTNO2, and CTNO3. These
!                              are no onger needed because NO, NO2, and NO3
!                              are now tracers.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  08 Nov 2013 - M. Sulprizio- Removed CTO3. It is no longer used because O3
!                              is now a tracer.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! INITIALIZE begins here!
      !=================================================================

      ! Error condition if IFLAG does not equal 2, or 3!
      IF ( IFLAG < 2 .or. IFLAG > 3 ) THEN
         CALL ERROR_STOP( 'Invalid IFLAG!', 'initialize.f' )
      ENDIF  

      !=================================================================
      ! If IFLAG=2 then zero the accumulating arrays
      !=================================================================
      IF ( IFLAG == 2 ) THEN

         ! Allocatable arrays are zeroed only if their
         ! respective diagnostics are turned on (bmy, 2/17/00)
         IF ( ND01 > 0 ) AD01     = 0e0
         IF ( ND02 > 0 ) AD02     = 0e0
         IF ( ND05 > 0 ) AD05     = 0e0
         IF ( ND06 > 0 ) AD06     = 0e0
         IF ( ND08 > 0 ) AD08     = 0e0
         IF ( ND11 > 0 ) AD11     = 0e0
         IF ( ND12 > 0 ) AD12     = 0e0
         IF ( ND14 > 0 ) CONVFLUP = 0d0
         IF ( ND15 > 0 ) TURBFLUP = 0d0
         IF ( ND16 > 0 ) AD16     = 0e0
         IF ( ND17 > 0 ) AD17     = 0e0
         IF ( ND18 > 0 ) AD18     = 0e0         
         IF ( ND19 > 0 ) AD19     = 0e0
         IF ( ND22 > 0 ) AD22     = 0e0
         IF ( ND24 > 0 ) MASSFLEW = 0d0
         IF ( ND25 > 0 ) MASSFLNS = 0d0
         IF ( ND26 > 0 ) MASSFLUP = 0d0
         IF ( ND28 > 0 ) AD28     = 0e0
         IF ( ND29 > 0 ) AD29     = 0e0
         IF ( ND30 > 0 ) AD30     = 0e0
         IF ( ND31 > 0 ) AD31     = 0e0
         IF ( ND33 > 0 ) AD33     = 0e0
         IF ( ND34 > 0 ) AD34     = 0e0
         IF ( ND35 > 0 ) AD35     = 0e0
         IF ( ND37 > 0 ) AD37     = 0e0
         IF ( ND38 > 0 ) AD38     = 0e0
         IF ( ND39 > 0 ) AD39     = 0e0
         IF ( ND43 > 0 ) AD43     = 0e0
         IF ( ND44 > 0 ) AD44     = 0e0
         IF ( ND45 > 0 ) AD45     = 0e0
         IF ( ND46 > 0 ) AD46     = 0e0
         IF ( ND47 > 0 ) AD47     = 0e0
         IF ( ND52 > 0 ) AD52     = 0e0
         IF ( ND54 > 0 ) AD54     = 0e0
         IF ( ND55 > 0 ) AD55     = 0e0
         !potential temp (fp, 6/2009)
         IF ( ND57 > 0 ) AD57     = 0e0
         IF ( ND58 > 0 ) AD58     = 0e0
         IF ( ND60 > 0 ) AD60     = 0e0
         IF ( ND66 > 0 ) AD66     = 0e0
         IF ( ND67 > 0 ) AD67     = 0e0
         IF ( ND68 > 0 ) AD68     = 0e0
         IF ( ND69 > 0 ) AD69     = 0e0

         ! For ND03 - mercury simulations (eck, sas, bmy, 1/20/05)
         IF ( ND03 > 0 ) THEN
            CALL ZERO_DIAG03
         ENDIF

         ! For ND53: POPS simulations (eck, 9/20/10)
         IF ( ND53 > 0 ) THEN
            CALL ZERO_DIAG53
         ENDIF

         ! For ND04 - CO2 simulation (pns, bmy, 7/26/05)
         IF ( ND04 > 0 ) THEN
            CALL ZERO_DIAG04
         ENDIF

         ! ND07 -- carbon aerosol emissions (rjp, tdf, bmy, 4/5/04)
         IF ( ND07 > 0 ) THEN
            AD07       = 0e0
            AD07_BC    = 0e0
            AD07_OC    = 0e0
            AD07_HC    = 0e0
            AD07_SOAGM = 0e0
#if   defined( APM )
            AD07_OM    = 0e0
#endif
         ENDIF

         ! ND09 -- HCN & CH3CN simulation (xyp, bmy, 6/27/05)
         IF ( ND09 > 0 ) THEN
            AD09     = 0e0
            AD09_em  = 0e0
         ENDIF
         
         ! For ND10 -- H2 & HD simulation (phs, 9/18/07)
         IF ( ND10 > 0 ) THEN 
            AD10    = 0e0
            AD10em  = 0e0
         ENDIF

         ! For ND13 - sulfur emissions (bmy, 6/6/00, 5/20/04)
         IF ( ND13 > 0 ) THEN
            AD13_DMS    = 0e0  
            AD13_SO2_ac = 0e0
            AD13_SO2_an = 0e0
            AD13_SO2_bb = 0e0
            AD13_SO2_bf = 0e0
            AD13_SO2_nv = 0e0
            AD13_SO2_ev = 0e0
            AD13_SO2_sh = 0e0
            AD13_SO4_an = 0e0
            AD13_SO4_bf = 0e0
            AD13_NH3_an = 0e0
            AD13_NH3_na = 0e0
            AD13_NH3_bb = 0e0
            AD13_NH3_bf = 0e0
         ENDIF

         ! ND21 -- optical depths
         IF ( ND21 > 0 ) THEN
            AD21 = 0e0
            IF ( LCRYST ) AD21_cr = 0e0
         ENDIF

         ! For ND32 -- NOx source diagnostics (bmy, 3/28/00)
         IF ( ND32 > 0 ) THEN 
            AD32_ac = 0e0
            AD32_an = 0e0
            AD32_bb = 0e0
            AD32_bf = 0e0
            AD32_fe = 0e0
            AD32_li = 0e0
            AD32_so = 0e0
            AD32_ub = 0e0
            AD32_ship       = 0e0
            AD32_ship_count = 0
         ENDIF

         IF ( ND36 > 0 ) THEN
            AD36            = 0e0
            AD36_SHIP       = 0e0
            AD36_SHIP_COUNT = 0
         ENDIF

         ! For ND41 - afternoon PBL heights (bmy, 2/17/05)
         IF ( ND41 > 0 ) THEN
            CALL ZERO_DIAG41
         ENDIF

         ! For ND42 - SOA concentrations (dkh, bmy, 5/22/06)
         IF ( ND42 > 0 ) THEN
            CALL ZERO_DIAG42
         ENDIF

         ! For ND56 - lightning flash rates (bmy, 5/5/06)
         IF ( ND56 > 0 ) THEN
            CALL ZERO_DIAG56
         ENDIF

#if   defined( TOMAS )
         ! For ND59 - size-resolved primary aerosol emissions (win, 7/28/09)
         IF ( ND59 > 0 ) THEN 
            AD59_NUMB = 0e0
            AD59_SULF = 0e0
            AD59_SALT = 0e0
            AD59_ECIL = 0e0
            AD59_ECOB = 0e0
            AD59_OCIL = 0e0
            AD59_OCOB = 0e0
            AD59_DUST = 0e0
         ENDIF

         ! For ND60 - aerosol microphysics (win, 7/28/09)
         IF ( ND60 > 0 ) THEN
            AD60_COND = 0e0
            AD60_COAG = 0e0
            AD60_NUCL = 0e0
            AD60_AQOX = 0e0
            AD60_ERROR = 0e0
            AD60_SOA   = 0e0
         ENDIF
#endif

         ! For ND65 -- Chemical production & loss (bmy, 12/5/00)
         IF ( ND65 > 0 ) THEN
            AD65  = 0d0
            IF ( ALLOCATED( FAM_PL ) ) FAM_PL = 0d0      
         ENDIF

         ! Echo output
         IF ( am_I_Root ) THEN
            WRITE( 6, '(a)' ) '     - INITIALIZE: Diag arrays zeroed!'
         ENDIF
      ENDIF  

      !================================================================= 
      ! If IFLAG=3 then zero the counter variables & arrays
      !=================================================================
      IF ( IFLAG == 3 ) THEN

         ! Now reset timesteps here for now 
         CALL SET_CT_A1(   RESET=.TRUE. )
         CALL SET_CT_A3(   RESET=.TRUE. )
         CALL SET_CT_A6(   RESET=.TRUE. )
         CALL SET_CT_CHEM( RESET=.TRUE. )
         CALL SET_CT_CONV( RESET=.TRUE. )
         CALL SET_CT_DYN(  RESET=.TRUE. )
         CALL SET_CT_EMIS( RESET=.TRUE. )
         CALL SET_CT_I3(   RESET=.TRUE. )
         CALL SET_CT_I6(   RESET=.TRUE. )
         CALL SET_CT_XTRA( RESET=.TRUE. )
         CALL SET_CT_DIAG( RESET=.TRUE. )

         ! Leave the ND48 counter for now
         KDA48     = 0

         ! Allocatable counter arrays
         IF ( ND16 >               0 ) CT16       = 0
         IF ( ND17 >               0 ) CT17       = 0
         IF ( ND18 >               0 ) CT18       = 0
         IF ( ND22 >               0 ) CTJV       = 0
         IF ( ND43 >               0 ) CTOH       = 0
         IF ( ND45 >               0 ) CTOTH      = 0
         IF ( ND47 > 0 .OR. ND65 > 0 ) CTO3_24h   = 0
         IF ( ND43 >               0 ) CTHO2      = 0
         ! update for arom (dkh, 06/21/07)  
         IF ( ND43 > 0 ) THEN
            IF ( ALLOCATED( CTLBRO2H ) ) CTLBRO2H   = 0
            IF ( ALLOCATED( CTLBRO2N ) ) CTLBRO2N   = 0
            IF ( ALLOCATED( CTLTRO2H ) ) CTLTRO2H   = 0
            IF ( ALLOCATED( CTLTRO2N ) ) CTLTRO2N   = 0
            IF ( ALLOCATED( CTLXRO2H ) ) CTLXRO2H   = 0
            IF ( ALLOCATED( CTLXRO2N ) ) CTLXRO2N   = 0
         ENDIF

         ! Echo output
         IF ( am_I_Root ) THEN
            WRITE( 6, '(a)' ) '     - INITIALIZE: Diag counters zeroed!'
         ENDIF
      ENDIF

      END SUBROUTINE INITIALIZE
!EOC
