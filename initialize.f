! $Id: initialize.f,v 1.6 2004/07/15 18:17:46 bmy Exp $
      SUBROUTINE INITIALIZE( IFLAG )
!
!******************************************************************************
!  Subroutine INITIALIZE (bmy, 6/15/98, 7/13/04) does the following:
!     (1) Zeroes globally defined GEOS-CHEM variables.
!     (2) Zeroes accumulating diagnostic arrays.
!     (3) Resets certain year/month/day and counter variables used 
!         in GEOS-CHEM diagnostic subroutines.
!
!  NOTE: Eventually we will fold this into "diag_mod.f" in a cleaner,
!        more consistent fashion.  Think about this later (bmy, 11/14/02)
!
!  Arguments as Input/Output
!  ============================================================================
!  (1 ) IFLAG : IFLAG=1, zero global CTM arrays 
!             : IFLAG=2, zero accumulating diagnostic arrays
!             : IFLAG=3, zero accumulating diagnostic counters
!
!  CTM arrays passed via COMMON blocks:
!  ============================================================================
!  (2 ) XTRA2       : Contains global boundary layer height in # of layers
!
!  Allocatable arrays passed via F90 module "diag_mod.f"
!  ============================================================================
!  (-1) AD11        : ND11 array -- acetone source diagnostic
!  (0 ) AD12        : ND12 array -- boundary layer emissions in "setemis.f"
!  (1 ) AD13_DMS    : ND13 array -- DMS emissions
!  (2 ) AD13_SO2_ac : ND13 array -- SO2 aircraft emissions
!  (3 ) AD13_SO2_an : ND13 array -- SO2 anthro emissions
!  (4 ) AD13_SO2_bb : ND13 array -- SO2 biomass emissions
!  (4a) AD13_SO2_bf : ND13 array -- SO2 biofuel emissions
!  (5 ) AD13_SO2_nv : ND13 array -- SO2 non-eruptive volcano emissions
!  (6 ) AD13_SO2_ev : ND13 array -- SO2 eruptive volcano emissions
!  (6a) AD13_SO2_sh : ND13 array -- SO2 ship emissions
!  (7 ) AD13_SO4_an : ND13 array -- SO4 anthro emissions
!  (8 ) AD13_NH3_an : ND13 array -- NH3 anthro emissions
!  (8a) AD13_NH3_na : ND13 array -- NH3 natural source emissions
!  (9 ) AD13_NH3_bb : ND13 array -- NH3 biomass emissions
!  (10) AD13_NH3_bf : ND13 array -- NH3 biofuel emissions
!  (11) CONVFLUP    : ND14 array -- cloud convection fluxes
!  (12) TURBFLUP    : ND15 array -- mass change in BL mixing
!  (13) AD16        : ND16 array -- precip fractions for wetdep
!  (14) AD17        : ND17 array -- rainout fractions
!  (15) AD18        : ND18 array -- washout fractions
!  (16) AD21        : ND21 array -- optical depths, cloud fractions
!  (17) AD22        : ND22 array -- J-values
!  (18) DIAGCHLORO  : ND23 array -- CH3CCl3 lifetime 
!  (19) MASSFLEW    : ND24 array -- E-W transport fluxes
!  (20) MASSFLNS    : ND25 array -- N-S transport fluxes
!  (21) MASSFLUP    : ND26 array -- vertical transport fluxes
!  (22) AD31        : ND31 array -- Psurface - PTOP
!  (23) AD33        : ND33 array -- tropopsheric sum of tracer
!  (24) AD32_ac     : ND32 array -- NOx source from aircraft 
!  (25) AD32_an     : ND32 array -- NOx source from anthro emissions 
!  (26) AD32_bb     : ND32 array -- NOx source from biomass burning
!  (27) AD32_bf     : ND32 array -- NOx source from biofuel burning
!  (28) AD32_fe     : ND32 array -- NOx source from fertilizers
!  (29) AD32_li     : ND32 array -- NOx source from lightning
!  (30) AD32_so     : ND32 array -- NOx source from soils
!  (31) AD32_ub     : ND32 array -- NOx source from upper boundary
!  (32) AD34        : ND34 array -- biofuel burning emissions
!  (33) AD35        : ND35 array -- tracer at 500 mb
!  (34) AD37        : ND37 array -- wet scavenging fraction
!  (35) AD38        : ND38 array -- rainout in wet conv
!  (36) AD39        : ND39 array -- washout in aerosol deposition
!  (37) AD41        : ND41 array -- afternoon PBL depths
!  (38) AD43        : ND43 array -- OH, NO concentrations
!  (39) AD45        : ND45 array -- tracer concentrations
!  (40) AD47        : ND47 array -- 24-h avg'd tracer conc.
!  (41) TCOBOX      : ND48 array -- station time series
!  (42) AD55        : ND55 array -- tropopause quantities
!  (43) AD65        : ND65 array -- chemical prod & loss
!  (44) FAMPL       : ND65 array -- accumulator for chemical prod & loss
!  (45) AD66        : ND66 array -- DAO 3-D fields
!  (46) AD67        : ND67 array -- DAO surface fields
!  (47) AD68        : ND68 array -- boxheights, air mass, water vapor, 
!                                   Air number density 
!  (48) AD69        : ND69 array -- surface areas
!
!  Scalars & Counter variables passed via COMMON blocks
!  ============================================================================
!  (1 ) TAU0        : beginning of diagnostic interval
!  (2 ) NTAU0       : integer representation of TAU0
!  (3 ) IDAY0       : day at beginning of diagnostic interval
!  (4 ) TOFDY0      : GMT at beginning of diagnostic interval
!  (5 ) JDATE0      : day of month  at beginning of diagnostic interval
!  (6 ) JMNTH0      : month of year at beginning of diagnostic interval
!  (7 ) JYEAR0      : year          at beginning of diagnostic interval
!  (8 ) KDA48       : Counter for timeseries accumulation (ND48 diagnostic)
!  (9 ) KDACC       : Counter for DIAG1
!  (10) KDADYN      : Counter of dynamic timesteps
!  (11) KDACONV     : Counter of convective timesteps
!  (12) KDASRCE     : Counter of emission timesteps
!  (13) KDACHEM     : Counter of chemistry timesteps
!  (14) KDA3FLDS    : Counter for # of times A-3 fields are read
!  (15) KDA6FLDS    : Counter for # of times A-6 fields are read
!  (16) KDI6FLDS    : Counter for # of times I-6 fields are read
!  (17) KDKZZFLDS   : Counter for # of times KZZ fields are read
!
!  Dynamically allocatable counter variables passed via F90 Modules
!  ============================================================================
!  (1 ) CT16        : ND16 counter array 
!  (2 ) CT17        : ND17 counter array
!  (3 ) CT18        : ND18 counter array
!  (4 ) CTJV        : ND22 counter array
!  (5 ) AFTTOT      : ND41 counter array
!  (6 ) CTNO        : ND43 counter array -- NO
!  (7 ) CTOH        : ND43 counter array -- OH
!  (8 ) CTOTH       : ND45 counter array 
!
!  NOTES:
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
!******************************************************************************
! 
      ! References to F90 modules
      USE DIAG_MOD
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE TIME_MOD

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! Many arrays...
#     include "CMN_DIAG"  ! NDxx flags

      ! Arguments 
      INTEGER, INTENT(IN) :: IFLAG
 
      !=================================================================
      ! INITIALIZE begins here!
      !
      ! Error condition if IFLAG does not equal 1, 2, or 3!
      !=================================================================
      IF ( IFLAG < 1 .or. IFLAG > 3 ) THEN
         CALL ERROR_STOP( 'Invalid IFLAG!', 'initialize.f' )
      ENDIF  

      !=================================================================
      ! If IFLAG=1 then zero the following CTM variables:
      !    T, XTRA2, SIG, SIGE, DSIG, DIAGCHLORO 
      !=================================================================
      IF ( IFLAG == 1 ) THEN
         XTRA2 = 0d0

         ! Only zero DIAGCHLORO if ND23 is turned on.  DIAGCHLORO only 
         ! needs to be zeroed at the start of the run. (bmy, 11/30/99)
         IF ( ND23 > 0 ) DIAGCHLORO = 0d0
      ENDIF  

      !=================================================================
      ! If IFLAG=2 then zero the accumulating arrays
      !=================================================================
      IF ( IFLAG == 2 ) THEN
         TOFDY0 = 0e0

         ! Allocatable arrays are zeroed only if their
         ! respective diagnostics are turned on (bmy, 2/17/00)
         IF ( ND01 > 0 ) AD01     = 0e0
         IF ( ND02 > 0 ) AD02     = 0e0
         IF ( ND03 > 0 ) AD03     = 0e0
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
         IF ( ND21 > 0 ) AD21     = 0e0
         IF ( ND22 > 0 ) AD22     = 0e0
         IF ( ND24 > 0 ) MASSFLEW = 0d0
         IF ( ND25 > 0 ) MASSFLNS = 0d0
         IF ( ND26 > 0 ) MASSFLUP = 0d0
         IF ( ND28 > 0 ) AD28     = 0e0
         IF ( ND29 > 0 ) AD29     = 0e0
         IF ( ND31 > 0 ) AD31     = 0e0
         IF ( ND33 > 0 ) AD33     = 0e0
         IF ( ND34 > 0 ) AD34     = 0e0
         IF ( ND35 > 0 ) AD35     = 0e0
         IF ( ND36 > 0 ) AD36     = 0e0
         IF ( ND37 > 0 ) AD37     = 0e0
         IF ( ND38 > 0 ) AD38     = 0e0
         IF ( ND39 > 0 ) AD39     = 0e0
         IF ( ND41 > 0 ) AD41     = 0e0
         IF ( ND43 > 0 ) AD43     = 0e0
         IF ( ND44 > 0 ) AD44     = 0e0
         IF ( ND45 > 0 ) AD45     = 0e0
         IF ( ND46 > 0 ) AD46     = 0e0
         IF ( ND47 > 0 ) AD47     = 0e0
         IF ( ND48 > 0 ) TCOBOX   = 0d0
         IF ( ND55 > 0 ) AD55     = 0e0
         IF ( ND66 > 0 ) AD66     = 0e0
         IF ( ND67 > 0 ) AD67     = 0e0
         IF ( ND68 > 0 ) AD68     = 0e0
         IF ( ND69 > 0 ) AD69     = 0e0

         ! ND07 -- carbon aerosol emissions (rjp, tdf, bmy, 4/5/04)
         IF ( ND07 > 0 ) THEN
            AD07     = 0e0
            AD07_BC  = 0e0
            AD07_OC  = 0e0
            AD07_HC  = 0e0
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
            AD13_NH3_an = 0e0
            AD13_NH3_na = 0e0
            AD13_NH3_bb = 0e0
            AD13_NH3_bf = 0e0
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
         ENDIF

         ! For ND65 -- Chemical production & loss (bmy, 12/5/00)
         IF ( ND65 > 0 ) THEN
            AD65  = 0e0
            IF ( ALLOCATED( FAMPL ) ) FAMPL = 0d0      
         ENDIF

         ! Echo output
         WRITE( 6, '(a)' ) '     - INITIALIZE: Diag arrays zeroed!'
      ENDIF  

      !================================================================= 
      ! If IFLAG=3 then zero the counter variables & arrays
      !=================================================================
      IF ( IFLAG == 3 ) THEN

         ! Now reset timesteps here for now 
         CALL SET_CT_A3(   RESET=.TRUE. )
         CALL SET_CT_A6(   RESET=.TRUE. )
         CALL SET_CT_CHEM( RESET=.TRUE. )
         CALL SET_CT_CONV( RESET=.TRUE. )
         CALL SET_CT_DYN(  RESET=.TRUE. )
         CALL SET_CT_EMIS( RESET=.TRUE. )
         CALL SET_CT_I6(   RESET=.TRUE. )

         ! Leave the ND48 counter for now
         KDA48     = 0

         ! Allocatable counter arrays
         IF ( ND16 > 0 ) CT16   = 0
         IF ( ND17 > 0 ) CT17   = 0
         IF ( ND18 > 0 ) CT18   = 0
         IF ( ND22 > 0 ) CTJV   = 0
         IF ( ND41 > 0 ) AFTTOT = 0
         IF ( ND43 > 0 ) CTNO   = 0
         IF ( ND43 > 0 ) CTOH   = 0
         IF ( ND45 > 0 ) CTOTH  = 0
         IF ( ND43 > 0 ) CTNO2  = 0
         IF ( ND43 > 0 ) CTHO2  = 0
         IF ( ND43 > 0 ) CTNO3  = 0

         ! Echo output
         WRITE( 6, '(a)' ) '     - INITIALIZE: Diag counters zeroed!'
      ENDIF
 
      ! Return to calling program
      END SUBROUTINE INITIALIZE
