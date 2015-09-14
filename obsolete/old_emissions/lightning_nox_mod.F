!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: lightning_nox_mod
!
! !DESCRIPTION: Module LIGHTNING\_NOx\_MOD contains variables and routines for 
!  emitting NOx from lightning into the atmosphere.  Original code comes from 
!  the old GISS-II CTM's of Yuhang Wang, Gerry Gardner, \& Larry Horowitz.  
!\\
!\\
! !INTERFACE:
!
      MODULE LIGHTNING_NOx_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: LIGHTNING
      PUBLIC  :: EMLIGHTNING
      PUBLIC  :: CLEANUP_LIGHTNING_NOX
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: LIGHTDIST                
      PRIVATE :: FLASHES_CTH              
      PRIVATE :: GET_IC_CG_RATIO          
      PRIVATE :: READ_LOCAL_REDIST        
      PRIVATE :: GET_OTD_LIS_SCALE        
      PRIVATE :: INIT_LIGHTNING_NOX       
!
! !PUBLIC DATA MEMBERS:
!
      ! Lightning NOx emissions [molec/cm3/s]
      REAL*8, ALLOCATABLE, PUBLIC :: EMIS_LI_NOx(:,:,:)
!
! !REMARKS:
!  (1) MFLUX and PRECON methods are now deprecated (ltm, bmy, 7/9/09)
!  (2) Starting w/ GEOS-Chem v9-02, we read OTD-LIS local redistribution
!       data files contained in subdirectory lightning_NOx_201311/.
!                                                                             .
!  References:
!  ============================================================================
!  (1 ) Price & Rind (1992), JGR, vol. 97, 9919-9933.
!  (2 ) Price & Rind (1994), M. Weather Rev, vol. 122, 1930-1939.
!  (3 ) Allen & Pickering (2002), JGR, 107, D23, 4711, doi:10.1029/2002JD002066
!  (4 ) Hudman et al (2007), JGR, 112, D12S05, doi:10.1029/2006JD007912
!  (5 ) Sauvage et al, 2007, ACP, 
!        http://www.atmos-chem-phys.net/7/815/2007/acp-7-815-2007.pdf
!  (6 ) Ott et al., (2010), JGR
!  (7 ) Allen et al., (2010), JGR
!  (8 ) Murray et al., (2011), in prep.
!
! !REVISION HISTORY:
!  14 Apr 2004 - L. Murray, R. Hudman - Initial version
!  (1 ) Based on "lightning_nox_mod.f", but updated for near-land formulation
!        and for CTH, MFLUX, PRECON parameterizations (ltm, bmy, 5/10/06)
!  (2 ) Now move computation of IC/CG flash ratio out of routines FLASHES_CTH, 
!        FLASHES_MFLUX, FLASHES_PRECON, and into routine GET_IC_CG_RATIO.
!        Added a fix in LIGHTDIST for pathological grid boxes.  Set E_IC_CG=1 
!        according to Allen & Pickering [2002].  Rename OTDSCALE array to
!        OTD_REG_REDIST, and also add OTD_LOC_REDIST array.  Now scale 
!        lightning to 6 Tg N/yr for both 2x25 and 4x5.  Rename routine
!        GET_OTD_LIS_REDIST to GET_REGIONAL_REDIST.  Add similar routine
!        GET_LOCAL_REDIST.  Removed GET_OTD_LOCp AL_REDIST.  Bug fix: divide 
!        A_M2 by 1d6 to get A_KM2. (rch, ltm, bmy, 2/22/07)
!  (3 ) Rewritten for separate treatment of LNOx emissions at tropics & 
!        midlatitudes, based on Hudman et al 2007.  Removed obsolete
!        variable E_IC_CG. (rch, ltm, bmy, 3/27/07)
!  (4 ) Changes implemented in this version (ltm, bmy, 10/3/07)
!        * Revert to not classifying near-land as land
!        * Eliminate NOx emisisons per path length entirely
!        * Scale tropics to 260 mol/fl constraint from Randall Martin's 
!           4.4 Tg and OTD-LIS avg ann flash rate
!        * Remove top-down scaling (remove the three functions)
!        * Allow option of mid-level scaling to match global avg ann flash
!           rate between G-C and OTD-LIS 11-year climatology (new function)
!        * Local Redist now a la Murray et al, 2007 in preparation (monthly)
!        * Replace GEMISNOX (from CMN_NOX) with module variable EMIS_LI_NOx
!  (5 ) Added MFLUX, PRECON redistribution options (ltm, bmy, 11/29/07)
!  (6 ) Updated OTD/LIS scaling for GEOS-5 to get more realistic totals
!        (ltm, bmy, 2/20/08)
!  (7 ) Now add the proper scale factors for the GEOS-5 0.5 x 0.666 grid
!        and the GEOS-3 1x1 nested N. America grid in routine 
!        GET_OTD_LIS_SCALE. (yxw, dan, ltm, bmy, 11/14/08)
!  (8 ) Added quick fix for GEOS-5 reprocessed met fields (ltm, bmy, 2/18/09)
!  (9 ) Added quick fix for GEOS-5 years 2004, 2005, 2008 (ltm, bmy, 4/29/09)
!  (10) Updated OTD/LIS scaling for GEOS-5 reprocessed data (ltm, bmy, 7/10/09)
!  (11) Updated for GEOS-4 1 x 1.25 grid (lok, ltm, bmy, 1/13/10)
!  (12) Reprocessed for CLDTOPS calculation error; Updated Ott vertical
!        profiles; Removal of depreciated options, e.g., MFLUX and PRECON;
!        GEOS5 5.1.0 vs. 5.2.0 special treatment; MERRA; Other changes.
!        Please see PDF on wiki page for full description of lightning
!        changes to v9-01-01. (ltm, 1/25/11)
!  13 Aug 2010 - R. Yantosca - Add modifications for MERRA
!  10 Nov 2010 - L. Murray   - Updated OTD/LIS local scaling for MERRA 4x5
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  02 Feb 2012 - R. Yantosca - Added modifications for GEOS-5.7.x met fields
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Scalars
      INTEGER              :: NNLIGHT
      REAL*8               :: AREA_30N
      REAL*8               :: OTD_LIS_SCALE

      ! Parameters
      INTEGER, PARAMETER   :: NLTYPE        = 4
      REAL*8,  PARAMETER   :: RFLASH_MIDLAT = 3.011d26   ! 500 mol/flash
      REAL*8,  PARAMETER   :: RFLASH_TROPIC = 1.566d26   ! 260 mol/flash
      REAL*8,  PARAMETER   :: EAST_WEST_DIV = -30d0
      REAL*8,  PARAMETER   :: WEST_NS_DIV   =  23d0
      REAL*8,  PARAMETER   :: EAST_NS_DIV   =  35d0
      REAL*8,  PARAMETER   :: T_NEG_BOT     = 273.0d0    !   0 C 
      REAL*8,  PARAMETER   :: T_NEG_CTR     = 258.0d0    ! -15 C
      REAL*8,  PARAMETER   :: T_NEG_TOP     = 233.0d0    ! -40 C

      ! Arrays
      REAL*8,  ALLOCATABLE :: PROFILE(:,:)
      REAL*8,  ALLOCATABLE :: SLBASE(:,:,:)
      REAL*8,  ALLOCATABLE :: OTD_REG_REDIST(:,:)
      REAL*8,  ALLOCATABLE :: OTD_LOC_REDIST(:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lightning
!
! !DESCRIPTION: Subroutine LIGHTNING uses Price \& Rind's formulation for 
!  computing NOx emission from lightning (with various updates).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LIGHTNING( State_Met )
!
! !USES:
!
      USE DIAG56_MOD,         ONLY : AD56,      ND56
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_YMID,  GET_XMID,   GET_AREA_M2
      USE LOGICAL_MOD,        ONLY : LOTDLOC
      USE PRESSURE_MOD,       ONLY : GET_PEDGE, GET_PCENTER
      USE TIME_MOD,           ONLY : GET_MONTH, GET_YEAR

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_GCTM_MOD             ! Physical constants
!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REMARKS:
!  Output Lightning NOX [molec/cm3/s] is stored in the EMIS_NOX_LI array.
! 
! !REVISION HISTORY: 
!  10 May 2006 - L. Murray - Initial version  
!  (1 ) Now recompute the cold cloud thickness according to updated formula 
!        from Lee Murray.  Rearranged argument lists to routines FLASHES_CTH, 
!        FLASHES_MFLUX, FLASHES_PRECON.  Now call READ_REGIONAL_REDIST and
!        READ_LOCAL_REDIST. Updated comments accordingly.  Now apply 
!        FLASH_SCALE to scale the total lightning NOx to 6 Tg N/yr.  Now apply
!        OTD/LIS regional or local redistribution (cf. B. Sauvage) to the ND56 
!        diagnostic. lightning redistribution to the ND56 diag.  Renamed
!        REGSCALE variable to REDIST.  Bug fix: divide A_M2 by 1d6 to get
!        A_KM2. (rch, ltm, bmy, 2/14/07)
!  (2 ) Rewritten for separate treatment of LNOx emissions at tropics & 
!        midlatitudes (rch, ltm, bmy, 3/27/07)
!  (3 ) Remove path-length algorithm.  Renamed from LIGHTNING_NL to LIGHTNING.
!        Other improvements. (ltm, bmy, 9/24/07)
!  (4 ) Remove depreciated options; Update to new Ott et al vertical profiles;
!        Reprocessed for bug in CLDTOPS calculation. See PDF on wiki for
!        full description of changes for v9-01-01. (ltm, bmy, 1/25,11)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST      = .TRUE.
      INTEGER, SAVE :: LASTMONTH  = -1
      INTEGER, SAVE :: LASTSEASON = -1
      INTEGER       :: I,         J,           L,        LCHARGE
      INTEGER       :: LMAX,      LTOP,        LBOTTOM,  L_MFLUX
      INTEGER       :: MONTH,     YEAR
      REAL*8        :: A_KM2,     A_M2,        CC,       DLNP     
      REAL*8        :: DZ,        FLASHRATE,   H0,       HBOTTOM
      REAL*8        :: HCHARGE,   IC_CG_RATIO, MFLUX,    P1
      REAL*8        :: P2,        P3,          RAIN,     RATE
      REAL*8        :: RATE_SAVE, REDIST,      T1,       T2
      REAL*8        :: TOTAL,     TOTAL_CG,    TOTAL_IC, X       
      REAL*8        :: YMID,      Z_IC,        Z_CG,     ZUP
      REAL*8        :: XMID
      REAL*8        :: VERTPROF(LLPAR)

      ! Pointers
      INTEGER, POINTER :: CLDTOPS(:,:)
      REAL*8,  POINTER :: BXHEIGHT(:,:,:)
      REAL*8,  POINTER :: T(:,:,:)

      !=================================================================
      ! LIGHTNING begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_LIGHTNING_NOX
         FIRST = .FALSE.
      ENDIF

      ! Initialize pointers
      BXHEIGHT => State_Met%BXHEIGHT
      CLDTOPS  => State_Met%CLDTOPS
      T        => State_Met%T

      ! LMAX: the highest L-level to look for lightning (usually LLPAR-1)
      LMAX   = LLPAR - 1

      ! Get current month
      MONTH  = GET_MONTH()

      ! Check if it's a new month
      IF ( MONTH /= LASTMONTH ) THEN
         
         ! OTD-LIS local redistribution: read monthly redist factors
         IF ( LOTDLOC ) THEN
            CALL READ_LOCAL_REDIST( MONTH )
         ENDIF

#if defined( GEOS_5 ) 
      ! Because of different convection in GEOS 5.1.0 and GEOS 5.2.0,
      ! this value is different before and after Sept 1, 2008. 
      ! So reset value at start of each month, just in case it's
      ! a 2008 simulation. (ltm,1/26/11)
      OTD_LIS_SCALE = GET_OTD_LIS_SCALE()
#endif

         ! Reset for next month
         LASTMONTH = MONTH
      ENDIF

      ! Array containing molecules NOx / grid box / 6h. 
      SLBASE = 0d0

      !=================================================================
      ! Compute lightning emissions for each (I,J) column
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,           J,        L,        A_M2,   A_KM2     )
!$OMP+PRIVATE( YMID,        XMID,     LCHARGE,  P1,     P2        )
!$OMP+PRIVATE( T1,          T2,       DLNP,     DZ,     P3        )
!$OMP+PRIVATE( ZUP,         HCHARGE,  LTOP,     H0,     Z_CG      )
!$OMP+PRIVATE( Z_IC,        LBOTTOM,  HBOTTOM,  CC,     FLASHRATE )
!$OMP+PRIVATE( IC_CG_RATIO, L_MFLUX,  MFLUX,    RAIN,   RATE      )
!$OMP+PRIVATE( X,           TOTAL_IC, TOTAL_CG, TOTAL,  REDIST    )
!$OMP+PRIVATE( RATE_SAVE,   VERTPROF                              )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over surface boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

            !%%% NOTE: Use L=1 for GRID_MOD functions.  This is OK for the
            !%%% existing GEOS-Chem with a pure cartesian grid, but may be an
            !%%% issue when interfaced with a GCM with a non-regular grid
            !%%% (bmy, 3/1/12)

            ! Grid box surface areas in [m2] and [km2]
            A_M2     = GET_AREA_M2( I, J, 1 )
            A_KM2    = A_M2 / 1d6

            ! Grid box latitude [degrees]
            YMID     = GET_YMID( I, J, 1 )

            ! Grid box longitude [degrees]
            XMID     = GET_XMID( I, J, 1 )

            ! Initialize
            LBOTTOM  = 0
            LCHARGE  = 0
            CC       = 0d0
            HCHARGE  = 0d0
            HBOTTOM  = 0d0
            REDIST   = 0d0
            TOTAL    = 0d0
            TOTAL_IC = 0d0
            TOTAL_CG = 0d0

            !===========================================================
            ! (1) FIND NEGATIVE CHARGE LAYER
            !
            ! LCHARGE is the L-value where the negative charge layer is
            ! found.  According to Williams (1985), the negative charge
            ! layer occurs where T is between 0 C and -40 C.  The 
            ! original model code set this at -10 C, but according to 
            ! Houze (1993), a good proxy for the negative charge layer
            ! maximum density is at -15 C.
            !
            ! Also of interest for later, will be the bottom of the
            ! negative charge layer (i.e., temp = 0 C) in calculating 
            ! the cold cloud depth.
            !
            ! If LCHARGE=1, then it is too cold to have water droplets
            ! in the column, so there will be no lightning events,
            ! and we go to the next (I,J) box.
            !
            ! (ltm, bmy, 5/10/06, 12/11/06)
            !===========================================================

            ! Find negative charge layer
            DO L = 1, LMAX
               IF ( T(I,J,L) <= T_NEG_CTR ) THEN
                  LCHARGE = L
                  EXIT
               ENDIF
            ENDDO

            ! Error check LCHARGE
            IF ( LCHARGE >= LMAX ) LCHARGE = LMAX
            IF ( LCHARGE <= 1    ) CYCLE

            !-----------------------------------------------------------
            ! (1a) Define more quantities
            !-----------------------------------------------------------
           
            ! Pressure [hPa] at the centers of grid
            ! boxes (I,J,LCHARGE-1) and (I,J,LCHARGE)
            P1   = GET_PCENTER( I, J, LCHARGE-1 )
            P2   = GET_PCENTER( I, J, LCHARGE   )

            ! Temperatures [K] at the centers of grid
            ! boxes (I,J,LCHARGE-1) and (I,J,LCHARGE)
            T1   = T(I,J,LCHARGE-1)
            T2   = T(I,J,LCHARGE  )
 
            ! DZ is the height [m] from the center of box (I,J,LCHARGE-1)
            ! to the negative charge layer.  It may be found in either
            ! the (LCHARGE)th sigma layer or the (LCHARGE-1)th layer.
            ! We use the hypsometric eqn to find the distance between
            ! the center of (LCHARGE)th and (LCHARGE-1)th boxes, then
            ! assume a linear temp distribution to scale between the two.
            DLNP = LOG( P1 / P2 ) / ( T1 - T2 ) * ( T1 - T_NEG_CTR )
            DZ   = Rdg0 * ( ( T1 + T2 ) / 2d0 ) * DLNP

            ! Pressure [hPa] at the bottom edge of box (I,J,LCHARGE),
            ! or, equivalently, the top edge of box (I,J,LCHARGE-1).
            P3   = GET_PEDGE( I, J, LCHARGE )

            ! Height [m] from the center of grid box (I,J,LCHARGE-1) 
            ! to the top edge of grid box (I,J,LCHARGE-1)
            ZUP  = Rdg0 * T1 * LOG( P1 / P3 )

            !-----------------------------------------------------------
            ! (1b) HCHARGE is the height of the negative charge layer 
            ! above the bottom edge of box (I,J,LCHARGE).  
            ! 
            ! If DZ < ZUP, then DZ is in grid box (I,J,LCHARGE-1);
            ! therefore subtract 1 from LCHARGE and compute HCHARGE 
            ! accordingly.  
            !
            ! In this case, please note that BXHEIGHT(I,J,LCHARGE)-ZUP 
            ! is the distance from the bottom edge of the grid box to 
            ! the center of the newly defined (LCHARGE)th layer.
            !-----------------------------------------------------------
            IF ( DZ >= ZUP ) THEN
               HCHARGE = DZ - ZUP
            ELSE
               LCHARGE = LCHARGE - 1
               HCHARGE = ( BXHEIGHT(I,J,LCHARGE) - ZUP ) + DZ
            ENDIF

            !===========================================================
            ! (2) COMPUTE CONVECTIVE CLOUD TOP HEIGHT
            !
            ! LTOP is the L-layer where the convective cloud top is 
            ! found.  The cloud top is located at the highest sigma 
            ! level for which the cloud mass flux is nonzero.  Since 
            ! GMAO cloud mass flux is defined at the top of each sigma 
            ! level, the convective cloud top is located at the top 
            ! edge of layer LTOP.
            !
            ! For lightning to exist, the cloud must straddle the 
            ! negative charge layer (in other words, at the very 
            ! minimum, the cloud bottom must occur in the LCHARGEth 
            ! layer).  If LTOP < LCHARGE go to the next (I,J) location.
            ! 
            ! Additionally, because the negative charge layer extends 
            ! from 0 C to around -40 C (Williams 1985), any cloud type 
            ! heights that are not colder than -40 C will be considered 
            ! unable to create the necessary dipole.  Therefore, if 
            ! T(I,J,LTOP) >= -40 C, go to the next (I,J) location. 
            !
            ! (ltm, bmy, 5/10/06, 12/11/06)
            !===========================================================

            ! Cloud top level
            LTOP = CLDTOPS(I,J)
            
            ! Error check LTOP
            IF ( LTOP == 0 ) CYCLE

            ! Error check LTOP as described above
            IF ( LTOP        >  LMAX      ) LTOP = LMAX
            IF ( LTOP        <  LCHARGE   ) CYCLE

#if    defined( GEOS_4 )

            !--------------------------
            ! GEOS-4 only
            !--------------------------
            ! Shallow-cloud inhibition trap (see Murray et al. [2011])
            IF ( T(I,J,LTOP) >= T_NEG_TOP ) CYCLE

#endif

            ! H0 is the convective cloud top height [m].  This is the
            ! distance from the surface to the top edge of box (I,J,LTOP).
            H0   = SUM( BXHEIGHT(I,J,1:LTOP) )

            ! Z_CG is the cloud-ground path (ground --> HCHARGE) [m]
            Z_CG = SUM( BXHEIGHT(I,J,1:LCHARGE-1)  ) + HCHARGE

            ! Z_IC is the intra-cloud path (HCHARGE --> cloud top) [m]
            Z_IC = SUM( BXHEIGHT(I,J,LCHARGE:LTOP) ) - HCHARGE

            !===========================================================
            ! (3) COMPUTE COLD CLOUD THICKNESS
            ! 
            ! Find the cold cloud thickness (CC) -- the distance from 
            ! where the temperature is 0 C up to the top of the cloud.  
            ! This is necessary for calculating the f_CG/f_IC ratio as 
            ! per Price and Rind 1993.  
            !
            ! This is a clone of the method above to find height to 
            ! HCHARGE, and we can recycle many of the same variables 
            ! that aren't used again.
            ! 
            ! Grid box (I,J,LBOTTOM) is the model layer where the 
            ! temperature of the cloud is 0C.
            !
            ! NOTE: If no temperature in the column is above 0 C, it 
            ! moves on to the next (I,J) box as before with the -15 C.
            !
            ! (ltm, bmy, 5/10/06, 12/11/06)
            !===========================================================

            ! Find the level where T = 0 C
            DO L = 1, LMAX
               IF ( T(I,J,L) <= T_NEG_BOT ) THEN
                  LBOTTOM = L
                  EXIT
               ENDIF
            ENDDO

            ! Error check LBOTTOM as described above
            IF ( LBOTTOM >= LMAX ) LBOTTOM = LMAX
            IF ( LBOTTOM <= 1    ) CYCLE 

            !-----------------------------------------------------------
            ! (3a) Define more quantities
            !-----------------------------------------------------------

            ! Pressure [hPa] at the centers of grid
            ! boxes (I,J,LBOTTOM-1) and (I,J,LBOTTOM)
            P1   = GET_PCENTER( I, J, LBOTTOM-1 )
            P2   = GET_PCENTER( I, J, LBOTTOM   )

            ! Temperature [K] at the centers of grid
            ! boxes (I,J,LBOTTOM-1) and (I,J,LBOTTOM)
            T1   = T(I,J,LBOTTOM-1)
            T2   = T(I,J,LBOTTOM  )
       
            ! DZ is the height [m] from the center of box (I,J,LCHARGE-1)
            ! to the negative charge layer.  It may be found in either
            ! the (LCHARGE)th sigma layer or the (LCHARGE-1)th layer.
            ! We use the hypsometric eqn to find the distance between
            ! the center of (LCHARGE)th and (LCHARGE-1)th boxes, then
            ! assume a linear temp distribution to scale between the two.
            DLNP = LOG( P1 / P2 ) / ( T1 - T2 ) * ( T1 - T_NEG_BOT )
            DZ   = Rdg0 * ( ( T1 + T2 ) / 2d0 ) * DLNP

            ! Pressure [hPa] at the bottom edge of box (I,J,LBOTTOM),
            ! or, equivalently, the top edge of box (I,J,BOTTOM-1).
            P3   = GET_PEDGE( I, J, LBOTTOM )

            ! Height [m] from the center of grid box (I,J,LBOTTOM-1) 
            ! to the top edge of grid box (I,J,LBOTTOM-1)
            ZUP  = Rdg0 * T1 * LOG( P1 / P3 )

            !-----------------------------------------------------------
            ! (3b) HBOTTOM is the height of the 0 C layer above the 
            ! bottom edge of box (I,J,LBOTTOM).  
            ! 
            ! If DZ < ZUP, then DZ is in grid box (I,J,LBOTTOM-1);
            ! therefore subtract 1 from LBOTTOM and compute HBOTTOM
            ! accordingly.
            !
            ! In this case, please note that BXHEIGHT(I,J,LBOTTOM)-ZUP 
            ! is the distance from the bottom edge of the grid box to 
            ! the center of the newly defined (LBOTTOM)th layer.
            !-----------------------------------------------------------
            IF ( DZ >= ZUP ) THEN
               HBOTTOM = DZ - ZUP
            ELSE
               LBOTTOM = LBOTTOM - 1
               HBOTTOM = ( BXHEIGHT(I,J,LBOTTOM) - ZUP ) + DZ
            ENDIF
  
            ! Cold cloud thickness is difference of cloud top 
            ! height (H0) and the height to the bottom.
            CC = H0 - SUM( BXHEIGHT(I,J,1:LBOTTOM-1) ) -
     &           HBOTTOM 

            !===========================================================
            ! (4) COMPUTE IC/CG FLASH_RATIO FROM COLD-CLOUD DEPTH
            !
            ! This is necessary as an input for the MFLUX and PRECON
            ! parameterizations, as well as for determining the fraction 
            ! of LNOX generated by either type of flash, and will
            ! eventually be used for separate vertical distributions
            ! when they become available.  (ltm, bmy, 12/11/06)
            !===========================================================
			
            ! Get Inter-Cloud/Cloud-Ground flash ratio [unitless]
            IC_CG_RATIO = GET_IC_CG_RATIO( CC )

            !===========================================================
            ! (5) COMPUTE LIGHTNING FLASH RATES
            !
            ! Now that we have computed the the ratio of intra-cloud
            ! flashes to cloud-ground flashes, compute the lightning
            ! flash rate via one of these parameterizations:
            !
            ! (a) Cloud top height (CTH)
            ! (b) Mass flux (MFLUX)
            ! (c) Convective Precpitation (PRECON)
            ! 
            ! (ltm, bmy, 5/10/06, 12/11/06)
            !===========================================================

            !--------------------------------------------------------
            ! (5a) CLOUD TOP HEIGHT PARAMETERIZATION (all met fields)
            !
            ! Based on Price & Rind [1992,1993,1994].
            !--------------------------------------------------------

            ! Get lightning flash rate per minute and IC/CG ratio
            CALL FLASHES_CTH( I, J, H0, FLASHRATE, State_Met )

            !===========================================================
            ! (6) COMPUTE TOTAL LNOx AND PARTITION INTO VERTICAL LAYERS
            ! 
            ! (6a) We convert FLASHRATE (computed above) to units of
            ! [flashes/6h] and store in the RATE variable.
            !
            ! We then multiply RATE by a scale factor based on 
            ! OTD/LIS observations.  This is necessary in order to make 
            ! sure that the lightning flashes happen in the correct 
            ! locations as diagnosed by OTD/LIS satellite observations.  
            ! There are two redistribution options:
            !
            !   (1) Apply regional scale factors based on OTD/LIS
            !        observations (method of L. Jourdain et al)
            !
            !   (2) Apply box-by-box scale scale factors based on
            !        OTD/LIS observations (method of B. Sauvage)
            !
            ! NOTE: As of 3/27/07, only method (1) is implemented.
            ! 
            ! (6b) We then compute X, which is the ratio
            !   [cloud-ground flashes / total flashes].
            !
            ! The amount of lightning released will depend whether we
            ! are in the tropics or in mid-latitudes.
            !
            !
            ! (6c) LIGHTNING NOx EMISSIONS IN THE TROPICS:
            ! ----------------------------------------------------------
            ! N. American / S. American     tropics: lat <= 23 N 
            ! African / Oceanian / Eurasian tropics: lat <= 35 N
            !
            ! The lightning NOx released in the inter-cloud (IC) and 
            ! cloud-ground (CG) paths are given by:
            !
            !   TOTAL_IC = RFLASH_TROPIC * RATE * (1-X) * Z_IC
            !   TOTAL_CG = RFLASH_TROPIC * RATE * (  X) * Z_CG
            !
            ! where:
            !   RFLASH_TROPIC = # of NOx molecules released per flash 
            !                    per meter (same as in previous code)
            !   RATE          = lightning flashes / 6h computed above
            !   Z_IC          = IC pathway in meters (from the negative    
            !                    cloud layer to the cloud top)
            !   Z_CG          = CG pathway in meters (from the negative
            !                    cloud layer to the ground surface)
            !   
            ! We also apply a top-down final global scaling factor, 
            ! calculated by previously bringing total global LNOx to 
            ! 6 Tg N/yr  (2x2.5: 0.3683, 4x5: 0.8996).  In 2004, the 
            ! tropics-only contribution to LNOx was 4.5379 Tg N.
            !
            ! 
            ! (6d) LIGHTING NOx EMISSIONS AT MIDLATITUDES:
            ! ----------------------------------------------------------
            ! N. American midlatitudes : lat > 23N
            ! Eurasian    midlatitudes : lat > 35N
            !
            ! The lightning NOx released at midlatitudes is independent
            ! of path length.  Thus:
            !
            !   TOTAL_IC = RFLASH_MIDLAT * RATE * (1-X) * MID_LAT_SCALE
            !   TOTAL_CG = RFLASH_MIDLAT * RATE *    X  * MID_LAT_SCALE
            !
            ! where 
            !   RFLASH_MIDLAT = # of NOx molecules released per flash 
            !                    per meter (based on 500 mol/flash)
            !   RATE          = lightning flashes / 6h computed above
            !   Z_IC          = IC pathway in meters (from the negative    
            !                    cloud layer to the cloud top)
            !   Z_CG          = CG pathway in meters (from the negative
            !                    cloud layer to the ground surface)
            !
            ! We now emit at the Northern Mid-latitudes using an RFLASH
            ! value of 500 mol/flash.  This is independent of path 
            ! length.  
            !
            ! NOTE: The OTD-LIS local redistribution method was expanded
            ! upon from Sauvage et al, 2007, ACP.
            ! http://www.atmos-chem-phys.net/7/815/2007/acp-7-815-2007.pdf
            !
            ! (6e) The total lightning is the sum of IC+CG paths:
            !   TOTAL = TOTAL_IC + TOTAL_CG
            !
            ! (6g) We then partition the NOx into each of the vertical 
            ! grid boxes within the column with a Ken Pickering PDF
            ! (see comments below).
            !
            ! (ltm, rch, bmy, 5/10/06, 3/27/07)
            !===========================================================

            !-----------------------------------------------------------
            ! (6a) Compute flash rate and apply OTD/LIS redistribution
            !-----------------------------------------------------------

            ! Convert [flashes/min] to [flashes/6h]
            RATE     = FLASHRATE * 360.0d0

            ! Get factors for OTD-LIS local redistribution or none.
            ! This constrains the seasonality and spatial distribution
            ! of the parameterized lightning to match the HRMC v2.2
            ! product from LIS/OTD, while still allowing the model to
            ! place lightning locally within deep convective events.
            ! (ltm, bmy, 1/31/07)
            IF ( LOTDLOC ) THEN
               REDIST = OTD_LOC_REDIST(I,J)
            ELSE
               REDIST = 1d0
            ENDIF

            ! Apply regional or local OTD-LIS redistribution so that the 
            ! flashes occur in the right place. 
            RATE = RATE * REDIST

            ! Apply scaling factor to make sure annual average flash rate 
            ! equals that of the climatology. (ltm, 09/24/07)
            RATE = RATE * OTD_LIS_SCALE

            !-----------------------------------------------------------
            ! (6b) Compute cloud-ground/total flash ratio
            !-----------------------------------------------------------

            ! Ratio of cloud-to-ground flashes to total # of flashes
            X    = 1d0 / ( 1d0 + IC_CG_RATIO )

            ! Compute LNOx emissions for tropics or midlats 

            IF ( XMID > EAST_WEST_DIV ) THEN 
               
               !--------------------------------------------------------
               ! (6c,6d) We are in EURASIA
               !--------------------------------------------------------
               IF ( YMID > EAST_NS_DIV ) THEN 


                  ! 6d: Eurasian Mid-Latitudes
                  TOTAL_IC = RFLASH_MIDLAT * RATE * ( 1d0 - X )
                  TOTAL_CG = RFLASH_MIDLAT * RATE * X
               ELSE                           

                  ! 6c: Eurasian Tropics
                  TOTAL_IC = RFLASH_TROPIC * RATE * ( 1d0 - X )
                  TOTAL_CG = RFLASH_TROPIC * RATE * X

               ENDIF

            ELSE   
               
               !--------------------------------------------------------
               ! (6c,6d) We are in the AMERICAS
               !--------------------------------------------------------
               IF ( YMID > WEST_NS_DIV ) THEN 

                  ! 6d: American Mid-Latitudes
                  TOTAL_IC = RFLASH_MIDLAT * RATE * ( 1d0 - X )
                  TOTAL_CG = RFLASH_MIDLAT * RATE * X

               ELSE                           

                  ! 6c: American Tropics
                  TOTAL_IC = RFLASH_TROPIC * RATE * ( 1d0 - X )
                  TOTAL_CG = RFLASH_TROPIC * RATE * X

               ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! (6e) Compute total lightning
            !-----------------------------------------------------------

            ! Sum of IC + CG [molec/6h]
            TOTAL = TOTAL_IC + TOTAL_CG

            !-----------------------------------------------------------
            ! (6f) ND56 diagnostic: store flash rates [flashes/min/km2]
            !-----------------------------------------------------------
            IF ( ND56 > 0 .and. RATE > 0d0 ) THEN

               ! Lightning flashes per minute per km2
               RATE_SAVE   = RATE / A_KM2 / 360d0

               ! Store total, IC, and CG flash rates in AD56
               AD56(I,J,1) = AD56(I,J,1) +   RATE_SAVE
               !AD56(I,J,2) = AD56(I,J,2) + ( RATE_SAVE * ( 1d0 - X ) )
               AD56(I,J,3) = AD56(I,J,3) + ( RATE_SAVE *         X   )

               AD56(I,J,2) = AD56(I,J,2) + H0 * 1d-3

            ENDIF

            !-----------------------------------------------------------
            ! (6g) LIGHTDIST computes the lightning distribution from 
            ! the ground to the convective cloud top using cumulative
            ! distribution functions for ocean flashes, tropical land
            ! flashes, and non-tropical land flashes, as specified by
            ! Lesley Ott [JGR, 2010]
            !-----------------------------------------------------------

            ! If there's lightning w/in the column ...
            IF ( TOTAL > 0d0 ) THEN

               ! Partition the column total NOx [molec/6h] from lightning 
               ! into the vertical using Pickering PDF functions
               CALL LIGHTDIST( I, J, LTOP, H0, YMID, TOTAL, VERTPROF,
     &                         State_Met )

               ! Add vertically partitioned NOx into SLBASE array
               DO L = 1, LLPAR
                  SLBASE(I,J,L) = SLBASE(I,J,L) + VERTPROF(L) 
               ENDDO
            ENDIF

         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointers
      NULLIFY( BXHEIGHT )
      NULLIFY( CLDTOPS  )
      NULLIFY( T        )

      END SUBROUTINE LIGHTNING
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lightdist
!
! !DESCRIPTION: Subroutine LIGHTDIST reads in the CDF used to partition the 
!  column lightning NOx into the GEOS-Chem vertical layers. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LIGHTDIST( I, J, LTOP, H0, XLAT, TOTAL, VERTPROF,
     &                      State_Met )
!
! !USES:
!
      USE DAO_MOD,            ONLY : IS_LAND,  IS_WATER
      USE DIRECTORY_MOD,      ONLY : DATA_DIR
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE FILE_MOD,           ONLY : IOERROR
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_YMID
      USE TIME_MOD,           ONLY : GET_MONTH

      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I          ! Longitude index
      INTEGER,        INTENT(IN)  :: J          ! Latitude index 
      INTEGER,        INTENT(IN)  :: LTOP       ! Level of conv cloud top
      REAL*8,         INTENT(IN)  :: H0         ! Conv cloud top height [m]
      REAL*8,         INTENT(IN)  :: XLAT       ! Latitude value [degrees]
      REAL*8,         INTENT(IN)  :: TOTAL      ! Column Total # of LNOx molec 
      TYPE(MetState), INTENT(IN)  :: State_Met  ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(OUT) :: VERTPROF(LLPAR) ! Vertical profile of LNOx
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Pickering et al., JGR 103, 31,203 - 31,316, 1998.
!  (2 ) Ott et al., JGR, 2010
!  (3 ) Allen et al., JGR, 2010
! 
! !REVISION HISTORY: 
!  18 Sep 2002 - M. Evans - Initial version (based on Yuhang Wang's code)
!  (1 ) Use functions IS_LAND and IS_WATER to determine if the given grid
!        box is over land or water.  These functions work for all DAO met
!        field data sets. (bmy, 4/2/02)
!  (2 ) Renamed M2 to LTOP and THEIGHT to H0 for consistency w/ variable names
!        w/in "lightning.f".  Now read the "light_dist.dat.geos3" file for 
!        GEOS-3 directly from the DATA_DIR/lightning_NOx_200203/ subdirectory.
!        Now read the "light_dist.dat" file for GEOS-1, GEOS-STRAT directly 
!        from the DATA_DIR/lightning_NOx_200203/ subdirectory.  Added 
!        descriptive comment header.  Now trap I/O errors across all 
!        platforms with subroutine "ioerror.f".  Updated comments, cosmetic 
!        changes.  Redimension FRAC(NNLIGHT) to FRAC(LLPAR). (bmy, 4/2/02)
!  (3 ) Deleted obsolete code from April 2002.  Now reference IU_FILE and
!        IOERROR from "file_mod.f".  Now use IU_FILE instead of IUNIT as the
!        file unit number. (bmy, 6/27/02)
!  (4 ) Now reference BXHEIGHT from "dao_mod.f" (bmy, 9/18/02)
!  (5 ) Bug fix: add GEOS_4 to the #if block (bmy, 3/4/04)
!  (6 ) Now bundled into "lightning_mod.f".  CDF's are now read w/in
!        routine INIT_LIGHTNING to allow parallelization (bmy, 4/14/04)
!  (7 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (8 ) Now uses near-land formulation (ltm, bmy, 5/10/06)
!  (9 ) Added extra safety check for pathological boxes (bmy, 12/11/06)
!  (10) Remove the near-land formulation, except for PRECON (ltm, bmy, 9/24/07)
!  (11) Now use the Ott et al. [2010] profiles, and apply consistently with
!        GMI model [Allen et al., 2010] (ltm, bmy, 1/25/11).
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  15 Jun 2012 - Nielsen - INQUIRE finds free logical unit number for IU_FILE
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: M, MTYPE, L, III, IOS, IUNIT, JJJ, MONTH
      REAL*8             :: ZHEIGHT, YMID
      REAL*8             :: FRAC(LLPAR)
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! LIGHTDIST begins here!
      !=================================================================

      ! Initialize 
      MTYPE    = 0
      VERTPROF = 0d0
      MONTH    = GET_MONTH()

      !%%% NOTE: Use L=1 for GRID_MOD functions.  This is OK for the
      !%%% existing GEOS-Chem with a pure cartesian grid, but may be an
      !%%% issue when interfaced with a GCM with a non-regular grid
      !%%% (bmy, 3/1/12)
      YMID     = GET_YMID( I, J, 1 )

      !=================================================================
      ! Test whether location (I,J) is continental, marine, or snow/ice
      !
      ! Depending on the combination of land/water and latitude, 
      ! assign a flag describing the type of lightning:
      !
      !   MTYPE = 1: ocean lightning
      !   MTYPE = 2: tropical continental lightning
      !   MTYPE = 3: midlatitude continental lightning 
      !   MTYPE = 4: subtropical lightning
      !             
      ! (ltm, bmy, 1/25/11)
      !=================================================================

      ! Assign profile kind to grid box, following Allen et al. 
      ! [JGR, 2010] (ltm, 1/25,11)
      SELECT CASE (MONTH)

      ! Southern Hemisphere Summer
      CASE ( 1,2,3,12 )

         IF ( ABS(YMID) .le. 15 ) THEN
            IF ( IS_LAND( I, J, State_Met ) ) THEN
               MTYPE = 2        ! Tropical continental
            ELSE
               MTYPE = 1        ! Tropical marine
            ENDIF
         ELSE IF ( ( YMID .gt. 15. ) .and. ( YMID .le. 30. ) ) THEN
            MTYPE = 4           ! N. Subtropics
         ELSE IF ( ( YMID .ge. -40. ) .and. ( YMID .lt. -15. ) ) THEN
            MTYPE = 4           ! S. Subtropics
         ELSE
            MTYPE = 3           ! Midlatitude
         ENDIF

      ! Equinox months
      CASE ( 4,5,10,11 )

         IF ( ABS(YMID) .le. 15 ) THEN
            IF ( IS_LAND( I, J, State_Met ) ) THEN
               MTYPE = 2        ! Tropical continental
            ELSE
               MTYPE = 1        ! Tropical marine
            ENDIF
         ELSE IF ( ABS(YMID) .le. 30 ) THEN
            MTYPE = 4           ! Subtropics
         ELSE
            MTYPE = 3           ! Midlatitude
         ENDIF

      ! Northern Hemisphere Summer
      CASE ( 6,7,8,9 )

         IF ( ABS(YMID) .le. 15 ) THEN
            IF ( IS_LAND( I, J, State_Met ) ) THEN
               MTYPE = 2        ! Tropical continental
            ELSE
               MTYPE = 1        ! Tropical marine
            ENDIF
         ELSE IF ( ( YMID .gt. 15. ) .and. ( YMID .le. 40. ) ) THEN
            MTYPE = 4           ! N. Subtropics
         ELSE IF ( ( YMID .ge. -30. ) .and. ( YMID .lt. -15. ) ) THEN
            MTYPE = 4           ! S. Subtropics
         ELSE
            MTYPE = 3           ! Midlatitude
         ENDIF
         
      END SELECT

      ! Extra safety check for pathological grid boxes (bmy, 11/29/06)
      IF ( MTYPE == 0 ) RETURN

      !=================================================================
      ! Use the CDF for this type of lightning to partition the total
      ! column lightning into the GEOS-3, GEOS-4, or GEOS-5 layers
      !=================================================================
      ZHEIGHT = 0.0

      ! Compute the height [km] at the top of each vertical level.
      ! Look up the cumulative fraction of NOx for each vertical level
      DO L = 1, LTOP
         ZHEIGHT = ZHEIGHT + State_Met%BXHEIGHT(I,J,L)
         FRAC(L) = PROFILE( NINT( ( ZHEIGHT/H0 )*3200. ), MTYPE ) *0.01
      ENDDO

      ! Convert from cumulative fraction to fraction for each level
      DO L = LTOP, 2, - 1
         FRAC(L) = FRAC(L) - FRAC(L-1)
      ENDDO 
      
      ! Partition lightning NOx by layer into VERTPROF
      DO L = 1, LTOP
         VERTPROF(L) = ( FRAC(L) * TOTAL )
      ENDDO

      END SUBROUTINE LIGHTDIST
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: flashes_cth
!
! !DESCRIPTION: Subroutine FLASHES\_CTH determines the rate of lightning 
!  flashes per minute based on the height of convective cloud tops, and the 
!  intra-cloud to cloud-ground strike ratio.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE FLASHES_CTH( I, J, HEIGHT, FLASHRATE, State_Met )
!
! !USES:
!

      USE DAO_MOD,            ONLY : IS_ICE
      USE DAO_MOD,            ONLY : IS_LAND
      USE DAO_MOD,            ONLY : IS_WATER
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I           ! Longitude index
      INTEGER,        INTENT(IN)  :: J           ! Latitude index
      REAL*8,         INTENT(IN)  :: HEIGHT      ! Height of conv cloud top [m]
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(OUT) :: FLASHRATE   ! Lightning flast rate
                                                 !  [flashes/min]
!
! !REVISION HISTORY: 
!  10 May 2006 - L. Murray - Initial version
!  (1  ) Subroutine renamed from FLASHES (ltm, bmy, 5/10/06)
!  (2  ) Remove CCTHICK, IC_CG_RATIO as arguments.  Remove computation of
!         IC_CG_RATIO and move that to GET_IC_CG_RATIO. (ltm, bmy, 12/11/06)
!  (3  ) Remove the near-land formulation (i.e. use function IS_LAND 
!         instead of IS_NEAR).(ltm, bmy, 9/24/07)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      !================================================================
      ! FLASHES_CTH begins here!
      !
      ! COMPUTE LIGHTNING FLASH RATE / MINUTE
      !
      ! Price & Rind (1992) give the following parameterizations for
      ! lightning flash rates as a function of convective cloud top
      ! height [km]:
      !
      !    FLAND  = 3.44e-5 * ( CLDTOP HEIGHT [km] ^ 4.9  )
      !    FOCEAN = 6.4e-4  * ( CLDTOP HEIGHT [km] ^ 1.73 )
      !
      ! Lightning will therefore occur much more often on land.  It 
      ! goes as approx. the 5th power of height, as opposed to approx. 
      ! the 2nd power of height over oceans.
      !
      ! We suppress lightning where the surface is mostly ice.  
      !
      ! (ltm, bmy, 5/10/06, 12/11/06)
      !================================================================

      ! Test for land type
      IF ( IS_LAND( I, J, State_Met ) ) THEN

         ! Flashes/min over land boxes
         FLASHRATE   = 3.44d-5 * ( ( HEIGHT * 1d-3 )**4.9d0  )

      ELSE IF ( IS_WATER( I, J, State_Met ) ) THEN

         ! Flahes/min over water
         FLASHRATE   = 6.4d-4  * ( ( HEIGHT * 1d-3 )**1.73d0 )

      ELSE IF ( IS_ICE( I, J, State_Met ) ) THEN

         ! Suppress lightning over snow/ice
         FLASHRATE   = 0d0

      ENDIF

      END SUBROUTINE FLASHES_CTH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ic_cg_ratio
!
! !DESCRIPTION: Function GET\_IC\_CG\_RATIO calculates the Intra-Cloud (IC) 
!  and Cloud-to-Ground (CG) lightning flash ratio based on the method of 
!  Price and Rind 1993, which is calculated from the cold-cloud depth 
!  (CCTHICK).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_IC_CG_RATIO( CCTHICK ) RESULT( IC_CG_RATIO )
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN)  :: CCTHICK       ! Cold cloud thickness [m]
!
! !RETURN VALUE:
!
      REAL*8               :: IC_CG_RATIO   ! Intra-cloud/cloud-ground ratio
!
! !REVISION HISTORY: 
!  11 Dec 2006 - R. Yantosca - Initial version
!  (1 ) Split off from FLASHES_CTH, FLASHES_MFLUX, FLASHES_PRECON into this
!        separate function (ltm, bmy, 12/11/06)
!  (2 ) Bug fix for XLF compiler (morin, bmy, 7/8/09)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8 :: CC, F_CG

      !=================================================================
      ! GET_IC_CG_RATIO begins here!
      !
      ! COMPUTE INTRA-CLOUD / CLOUD-GROUND FLASH RATIO 
      !
      ! Price & Rind (1993) compute the ratio of Cloud-Ground 
      ! to Total Flashes by the parameterization:
      !
      ! For 5.5 < dz < 14:
      !
      !     f_CG = 1 / (( A*dz^4 + B*dz^3 + C*dz^2 + D*dz + E ) + 1 )
      !
      ! For dz > 14:
      !
      !     f_CG = 0.02
      !                                        
      ! Where:
      !
      ! (1) dz is the depth [km] of the cloud above the freezing 
      !     level.  The cold-cloud thickness (dz) is the depth of 
      !     the layer between the cloud top and the center of the 
      !     highest layer for which the temperature exceeds 273 K. 
      !     The cold-cloud thickness is set to 5.5 km at grid points 
      !     where it is less than 5.5 km.
      !
      ! (2) The polynomial coefficients are:
      !        A=0.021,  B=-0.648,  C=7.493,  D=-36.54,  E=63.09
      !
      ! 
      ! Note: f_IC = 1 - f_CG
      ! 
      ! And hence,
      !
      !     IC_CG_RATIO = ( 1 - f_CG ) / f_CG
      !
      !
      ! IC_CG_RATIO is passed back to routine the LIGHTNING_NL, where
      ! it is passed to FLASHES_MFLUX and FLASHES_PRECON.  In these
      ! routines, the fraction of total lightning flashes that are 
      ! cloud-ground (CG) flashes is computed by:
      ! 
      !     F_CG        = 1d0 / ( 1d0 + IC_CG_RATIO )
      !
      ! and the fraction of the total lightning flashes that are
      ! intra-cloud (IC) flashes is computed by:
      !
      !     F_IC        = 1d0 - 1d0 / ( 1d0 + IC_CG_RATIO )
      !=====================================================================

      ! Convert cold cloud thickness from [m] to [km] (min value: 5.5 km)
      CC = MAX( CCTHICK * 1d-3, 5.5d0 )

      ! Compute cloud-ground flash ratio as described above
      IF ( CC > 14d0 ) THEN

         ! Constant value above 14 km
         F_CG = 0.02d0

      ELSE

         ! First create the polynomial expression
         F_CG = 63.09d0 + CC * ( -36.54d0  +
     &                    CC * (   7.493d0 + 
     &                    CC * (  -0.648d0 +
     &                    CC * (   0.021d0 ) ) ) )

         ! Then put it in the denominator
         F_CG = 1d0 / ( F_CG + 1d0 )
                  
      ENDIF

      ! Intra-Cloud / Cloud-Ground flash ratio
      IC_CG_RATIO = ( 1d0 - F_CG ) / F_CG

      END FUNCTION GET_IC_CG_RATIO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_local_redist
!
! !DESCRIPTION: Subroutine READ\_LOCAL\_REDIST reads in seasonal factors 
!  in order to redistribute GEOS-Chem flash rates according the "local 
!  redistribution" method of Bastien Sauvage.  This helps to make sure 
!  that the lightning flashes occur according to the distribution of 
!  observed convection.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_LOCAL_REDIST( MONTH )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE TIME_MOD,      ONLY : GET_TAU
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD                ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: MONTH   ! Current month
! 
! !REVISION HISTORY: 
!  26 Jan 2007 - B. Sauvage - Initial version
!  (1 ) Change from seasonal to monthly.  Rename all filenames from "v2"
!        to "v3". (ltm, bmy, 9/24/07)
!  (2 ) Change all filenames from "v2" to "v3".  Also now read from the 
!        directory lightning_NOx_200709. (ltm, bmy, 9/24/07)
!  (3 ) Added "quick fix" for reprocessed GEOS-5 met fields to be used when
!        the IN_CLOUD_OD switch is turned on. (ltm, bmy, 2/18/09)
!  (4 ) Now read from lightning_NOx_200907 directory for GEOS-4 and
!        GEOS-5 CTH parameterizations.  Updated OTD/LIS for GEOS-5 based on
!        4+ years of data; removed temporary fixes. (ltm, bmy, 7/10/09)
!  (5 ) Remove depreciated options and update to v5 of redist files in
!       new data directory. Special handling for GEOS5.1.0 and 5.2.0 added. 
!       (ltm, bmy, 1/25/11)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  02 Feb 2012 - R. Yantosca - Added modifications for GEOS-5.7.x met
!  18 Apr 2013 - R. Yantosca - Bug fix, prefix DATA_DIR to GEOS-5.7.x file
!  26 Sep 2013 - R. Yantosca - Remove SEAC4RS C-preprocessor switch
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  07 Nov 2013 - R. Yantosca - Now read files from lightning_NOx_201311 dir
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4             :: ARRAY(IIPAR,JJPAR,1)
      REAL*8             :: TAU0
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=9)   :: MODELNAME

      !=================================================================
      ! READ_LOCAL_REDIST begins here!
      !=================================================================

      ! Get model name
#if   defined( MERRA )
      MODELNAME = 'merra'
#elif defined( GEOS_FP )
      MODELNAME = 'geosfp'
#elif defined( GEOS_5 )
      ! GEOS-5 convection scheme changed in GCM post 9/1/2008
      IF ( GET_TAU() .ge. GET_TAU0( 9, 1, 2008 ) ) THEN
         MODELNAME = 'geos5.2.0'
      ELSE
         MODELNAME = 'geos5.1.0'
      ENDIF
#else
      MODELNAME = GET_NAME_EXT()
#endif

      ! OTD-LIS local redist filename for CTH parameterization
      ! Now read from lightning_NOx_201301 for GEOS-Chem v9-02+ (bmy, 11/7/13)
#if defined( GEOS_5 ) && defined( NESTED_EU ) 
      FILENAME =
     &        'lightning_NOx_201101/OTD-LIS-Local-Redist.CTH.v5.' //
     &        TRIM( MODELNAME ) // '.' // GET_RES_EXT()
#else
      FILENAME = 
     &        'lightning_NOx_201311/OTD-LIS-Local-Redist.CTH.v5.' // 
     &        TRIM( MODELNAME ) // '.' // GET_RES_EXT()
#endif

      ! Prefix directory to file name
      FILENAME = TRIM( DATA_DIR ) // TRIM( FILENAME ) 

      ! Append suffix for GEOS-FP nested grids
#if   defined( GEOS_FP ) && defined( NESTED_CH )
      FILENAME = TRIM( FILENAME ) // '.CH'
#elif defined( GEOS_FP ) && defined( NESTED_NA )
      FILENAME = TRIM( FILENAME ) // '.NA'
#elif defined( GEOS_FP ) && defined( NESTED_EU ) 
      FILENAME = TRIM( FILENAME ) // '.EU'
#elif defined( GEOS_FP ) && defined( NESTED_SE ) 
      FILENAME = TRIM( FILENAME ) // '.SE'
#endif

      ! Append suffix for GEOS-5 nested grids
#if   defined( GEOS_5 ) && defined( NESTED_CH ) 
      FILENAME = TRIM( FILENAME ) // '.CH'
#elif defined( GEOS_5 ) && defined( NESTED_NA ) 
      FILENAME = TRIM( FILENAME ) // '.NA'
#elif defined( GEOS_5 ) && defined( NESTED_EU ) 
      FILENAME = TRIM( FILENAME ) // '.EU'
#endif

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_LOCAL_REDIST: Reading ', a )

      ! Use "generic" year 1985 for time indexing
      TAU0 = GET_TAU0( MONTH, 1, 1985 )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'OTD-LOC',  1, 
     &                 TAU0,      IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )  
   
      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), OTD_LOC_REDIST )

      END SUBROUTINE READ_LOCAL_REDIST
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emlightning
!
! !DESCRIPTION: Subroutine EMLIGHTNING converts lightning emissions to 
!  [molec/cm3/s] and stores them in the GEMISNOX array, which gets passed 
!  to SMVGEAR.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMLIGHTNING( State_Met, State_Chm )

!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_STRAT
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD32_li
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TRACERID_MOD,       ONLY : IDTNO
!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
! 
! !REVISION HISTORY: 
!  09 Oct 1997 - R. Yantosca - Initial version
!  (1 ) Remove IOFF, JOFF from the argument list.  Also remove references
!        to header files "CMN_O3" and "comtrid.h" (bmy, 3/16/00)
!  (2 ) Now use allocatable array for ND32 diagnostic (bmy, 3/16/00)  
!  (3 ) Now reference BXHEIGHT from "dao_mod.f".  Updated comments, cosmetic
!        changes.  Replace LCONVM with the parameter LLCONVM. (bmy, 9/18/02)
!  (4 ) Removed obsolete reference to "CMN".  Now bundled into 
!        "lightning_mod.f" (bmy, 4/14/04)
!  (5 ) Renamed from EMLIGHTNING_NL to EMLIGHTNING.  Now replace GEMISNOX
!        (from CMN_NOX) with module variable EMIS_LI_NOx. (ltm, bmy, 10/3/07)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept State_Chm
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
      INTEGER             :: I,J,L
      REAL*8              :: TMP

      ! External functions
      REAL*8, EXTERNAL    :: BOXVL

      !=================================================================
      ! EMLIGHTNING begins here!
      !=================================================================
      DO L = 1, LLCONVM 
      DO J = 1, JJPAR
      DO I = 1, IIPAR
                
         ! No lightning emissions in the stratosphere (cdh, 4/25/2013)
         IF ( ITS_IN_THE_STRAT(I, J, L, State_Met) ) EXIT

         ! SLBASE(I,J,L) has units [molec NOx/6h/box], convert units:
         ! [molec/6h/box] * [6h/21600s] * [box/BOXVL cm3] = [molec/cm3/s]
         TMP = SLBASE(I,J,L) / ( 21600.d0 * BOXVL(I,J,L,State_Met) )

         EMIS_LI_NOx(I,J,L)                = TMP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% KLUDGE FOR SEAC4RS (bmy, 9/25/13)
!%%% For SEAC4RS, knock the emissions down by a factor of 5. Keep code here.
!%%%         ! ( skim, 8/15/13)
!%%%         EMIS_LI_NOx(I,J,L)                = TMP / 5.0
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! ND32 Diagnostic: Lightning NOx [molec NOx/cm2/s]
         IF ( ND32 > 0 ) THEN
            AD32_li(I,J,L) = AD32_li(I,J,L) + 
     &                     ( TMP * State_Met%BXHEIGHT(I,J,L) * 1d2 )
         ENDIF
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE EMLIGHTNING
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_otd_lis_scale
!
! !DESCRIPTION: Function GET\_OTD\_LIS\_SCALE returns a met-field dependent 
!  scale factor which is to be applied to the lightning flash rate to bring 
!  the annual average flash rate to match that of the OTD-LIS climatology 
!  ( ~ 45.9 flashes/sec ). Computed by running the model over the 11-year 
!  OTD-LIS campaign window and comparing the average flash rates, or as 
!  many years as are available.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_OTD_LIS_SCALE() RESULT( BETA )
!
! !USES:
!
      USE BPCH2_MOD,   ONLY : GET_TAU0
      USE ERROR_MOD,   ONLY : GEOS_CHEM_STOP
      USE TIME_MOD,    ONLY : GET_TAU, GET_MONTH, GET_YEAR
!
! !RETURN VALUE:
!
      REAL*8 :: BETA    ! Scale factor
!
! !REMARKS:
!  (1) Starting in G-C v9-02, we now read data from lightning_NOx_201311.
! 
! !REVISION HISTORY: 
!  24 Sep 2007 - L. Murray - Initial version
!  (1 ) Added MFLUX, PRECON scaling for GEOS-4.  Also write messages for met
!        field types/grids where scaling is not defined. (ltm, bmy, 11/29/07)
!  (2 ) Now use different divisor for local redist (ltm, bmy, 2/20/08)
!  (3 ) Now compute the proper scale factor for GEOS-5 0.5 x 0.666 grids
!        and the GEOS-3 1x1 nested NA grid (yxw, dan, ltm, bmy, 11/14/08)
!  (4 ) Added "quick fix" for reprocessed GEOS-5 met fields to be used when 
!        the IN_CLOUD_OD switch is turned on. (ltm, bmy, 2/18/09)
!  (5 ) Added "quick fix" for 2004, 2005, 2008 OTD/LIS (ltm, bmy, 4/29/09)
!  (6 ) Updated scale factors for GEOS-5 based on 4+ years of data.  Remove
!        temporary fixes. (bmy, 7/10/09)
!  (7 ) Modification for GEOS-4 1 x 1.25 grid (lok, ltm, bmy, 1/13/10)
!  (8 ) Reprocessed for error in CLDTOPS field; Updated for GEOS
!        5.1.0 vs. 5.2.0; MERRA added; (ltm, bmy, 1/25/11)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  02 Feb 2012 - R. Yantosca - Compute BETA for MERRA 2 x 2.5
!  02 Feb 2012 - R. Yantosca - Compute BETA for GEOS-5.7.x
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! Define the average annual flash rate (flashes per second), as
      ! calculated from the OTD-LIS HR Monthly Climatology observations
      ! from May 1995 through Dec 2005.  Slight difference when
      ! averaging over different resolutions. (ltm, 09/24/07, 11/14/08)
      !=================================================================
#if   defined( GRID2x25 ) 
      REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 45.8650d0
#elif defined( GRID4x5  )
      REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 45.8658d0
#elif defined( GRID1x125 )
      REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 45.8655d0
#elif defined( GRID05x0666 ) && defined( NESTED_CH )
      REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 8.7549280d0
#elif defined( GRID05x0666 ) && defined( NESTED_NA )
      REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 6.9685368d0
#elif defined( GRID025x03125 ) && defined( NESTED_NA )
      REAL*8, PARAMETER     :: ANN_AVG_FLASHRATE = 6.7167603d0
#endif

#if   defined( GEOS_5 )

      ! Are we using GEOS 5.2.0 or GEOS 5.1.0?
      LOGICAL               :: GEOS_520

      ! Lightning is sensitive to which convection scheme
      ! is used in the GCM used for the data assimilation.
      ! GEOS-5 changed its scheme in met fields following 9/1/2008,
      ! and requires special treatment. (ltm, 1/25/11)
      if ( GET_TAU() .GE. GET_TAU0( 9, 1, 2008 ) ) then
         GEOS_520 = .TRUE.      ! Using GEOS 5.2.0
      else
         GEOS_520 = .FALSE.     ! Using GEOS 5.1.0
      endif

#endif

      !=================================================================
      ! GET_OTD_LIS_SCALE begins here!
      !=================================================================

      ! The lightning flash rate equations are sensitive to model resolution
      ! and convection scheme used in the data assimilation.
      ! We know from the LIS/OTD satellite products that the global annual
      ! average flash rate is 46 fl s-1. We determine a single scaling 
      ! factor, beta, to be applied uniformly 
      !
      ! beta =  ( Annual Average Flash Rate Observed ) / 
      !           ( Annual Average Flash Rate Unconstr Parameterization )
      !
      ! This is equivalent to modifying the first coefficient of the 
      ! Price and Rind [1992] formulation to get the right magnitude 
      ! for a given model framework.
      ! 
      ! Beta corresponds to beta in Murray et al. [2011]
      !
      ! Note: GEOS-5 requires separate factors for GEOS 5.2.0 and 5.1.0.
      ! (ltm, 1/25/11)

      ! Initialize
      BETA = 1d0

#if   defined( GEOS_FP ) && defined( GRID4x5 ) 

      !---------------------------------------
      ! GEOS-FP: 4 x 5 global simulation
      !---------------------------------------

      ! Constrained with simulated "climatology" for
      ! April 2012 - Sept 2013. Will need to be updated as more
      ! met fields become available (ltm, 11/07/13).
      IF ( ( GET_YEAR() .eq. 2012 .and. GET_MONTH() .ge. 4 ) .or.
     &     ( GET_YEAR() .eq. 2013 .and. GET_MONTH() .le. 9 ) ) THEN
         BETA = ANN_AVG_FLASHRATE / 82.003230d0
      ENDIF

#elif defined( GEOS_FP ) && defined( GRID2x25 )

      !---------------------------------------
      ! GEOS-FP: 2 x 2.5 global simulation
      !---------------------------------------

      ! Constrained with simulated "climatology" for
      ! April 2012 - Sept 2013. Will need to be updated as more
      ! met fields become available (ltm, 01/15/14).
      IF ( ( GET_YEAR() .eq. 2012 .and. GET_MONTH() .ge. 4 ) .or.
     &     ( GET_YEAR() .eq. 2013 .and. GET_MONTH() .le. 9 ) ) THEN
         BETA = ANN_AVG_FLASHRATE / 257.93269d0
      ENDIF

#elif defined( GEOS_FP ) && defined( GRID025x0325 ) && defined( NESTED_CH )

      !---------------------------------------
      ! GEOS-FP: Nested China simulation
      !---------------------------------------
      
      ! ltm: Will need to be determined when met fields become available.

#elif defined( GEOS_FP ) && defined( GRID025x03125 ) && defined( NESTED_NA )

      !---------------------------------------
      ! GEOS-FP: Nested SEAC4RS simulation
      !---------------------------------------

      ! Constrained with simulated "climatology" for
      ! April 2012 - Sept 2013. Will need to be updated as more
      ! met fields become available (ltm, 11/14/13).
      IF ( ( GET_YEAR() .eq. 2012 .and. GET_MONTH() .ge. 4 ) .or. 
     &     ( GET_YEAR() .eq. 2013 .and. GET_MONTH() .le. 9 ) ) THEN
         BETA = ANN_AVG_FLASHRATE / 652.44105d0
      ENDIF

#elif defined( MERRA ) && defined( GRID2x25 )

      !---------------------------------------
      ! MERRA: 2 x 2.5 global simulation
      !---------------------------------------
      BETA = ANN_AVG_FLASHRATE / 253.55888d0

#elif defined( MERRA ) && defined( GRID4x5 )

      !---------------------------------------
      ! MERRA: 4 x 5 global simulation
      !---------------------------------------
      BETA = ANN_AVG_FLASHRATE / 76.019042d0

#elif defined( GEOS_5 ) && defined( GRID05x0666 ) && defined( NESTED_NA)

      !---------------------------------------
      ! GEOS 5: 0.5 x 0.666
      ! Nested grid simulation: North America
      !---------------------------------------
      if ( GEOS_520 ) then
         ! Constrained with simulated climatology for
         ! Sept 2009 - May 2013 (ltm, 11/07/13)
         BETA = ANN_AVG_FLASHRATE / 170.05559d0
      else
         BETA = ANN_AVG_FLASHRATE / 160.51908d0
      endif

      ! Discourage users from using lightning outside the constraint period.
      ! You may comment out these lines, but should verify that lightning
      ! doesn't become unreasonably high anywere in the domain. (ltm, 11/07/13)
      IF (   GET_YEAR() .ge. 2014 .or.
     &    ( GET_YEAR() .eq. 2013 .and. GET_MONTH() .gt. 5 ) ) BETA = 1d0

#elif defined( GEOS_5 ) && defined( GRID05x0666 ) && defined( NESTED_CH)

      !---------------------------------------
      ! GEOS 5: 0.5 x 0.666
      ! Nested grid simulation: China
      !---------------------------------------
      if ( GEOS_520 ) then
         BETA = ANN_AVG_FLASHRATE / 573.24835d0
      else
         BETA = ANN_AVG_FLASHRATE / 546.56367d0
      endif  

#elif defined( GEOS_5 ) && defined( GRID2x25 )

      !---------------------------------------
      ! GEOS 5: 2 x 2.5 global simulation
      !---------------------------------------
      if ( GEOS_520 ) then
         BETA = ANN_AVG_FLASHRATE / 221.72962d0
      else
         BETA = ANN_AVG_FLASHRATE / 199.54964d0
      endif

#elif defined( GEOS_5 ) && defined( GRID4x5 )
      
      !---------------------------------------
      ! GEOS 5: 4 x 5 global simulation
      !---------------------------------------
      if ( GEOS_520 ) then
         BETA = ANN_AVG_FLASHRATE / 70.236997d0
      else
         BETA = ANN_AVG_FLASHRATE / 64.167893d0
      endif

#elif defined( GEOS_4 ) && defined( GRID2x25 )
      
      !---------------------------------------
      ! GEOS 4: 2 x 2.5 global simulation
      !---------------------------------------
      BETA = ANN_AVG_FLASHRATE / 83.522403d0

#elif defined( GEOS_4 ) && defined( GRID4x5 )

      !---------------------------------------
      ! GEOS 4: 4 x 5 global simulation
      !---------------------------------------
      BETA = ANN_AVG_FLASHRATE / 29.359449d0

#elif   defined( GCAP )
      
      !---------------------------------------
      ! GCAP: 4 x 5 global simulation
      !---------------------------------------
      BETA = ANN_AVG_FLASHRATE / 48.681763d0

#endif

      IF ( BETA .eq. 1d0 ) THEN

         WRITE( 6,* ) 'Your model framework has not had its'
         WRITE( 6,* ) 'lightning code reprocessed for the correction'
         WRITE( 6,* ) 'to how CLDTOPS are calculated, probably due to'
         WRITE( 6,* ) 'the lack of your met fields at Harvard.'
         WRITE( 6,* ) ''
         WRITE( 6,* ) 'Please contact Lee Murray'
         WRITE( 6,* ) '(ltmurray@post.harvard.edu), who can help you'
         WRITE( 6,* ) 'prepare the necessary modifications and files'
         WRITE( 6,* ) 'to get lightning working for you.'
         WRITE( 6,* ) ''
         WRITE( 6,* ) 'You may remove this trap in lightning_nox_mod.f'
         WRITE( 6,* ) 'at your own peril, but be aware that the'
         WRITE( 6,* ) 'magnitude and distribution of lightning may be'
         WRITE( 6,* ) 'unrealistic.'
         
        ! CALL GEOS_CHEM_STOP
         
      ENDIF

      END FUNCTION GET_OTD_LIS_SCALE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lightning_NOx
!
! !DESCRIPTION: Subroutine INIT\_LIGHTNING\_NOx allocates all module arrays.  
!  It also reads the lightning CDF data from disk before the first lightning 
!  timestep. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_LIGHTNING_NOx
!
! !USES:
!
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE FILE_MOD,      ONLY : IOERROR
      USE GRID_MOD,      ONLY : GET_AREA_M2
      USE LOGICAL_MOD,   ONLY : LOTDLOC

      USE CMN_SIZE_MOD      ! Size parameters
! 
! !REVISION HISTORY:
!  14 Apr 2004 - R. Yantosca - Initial version
!  (1 ) Now reference DATA_DIR from "directory_mod.f"
!  (2 ) Now call GET_MET_FIELD_SCALE to initialize the scale factor for
!        each met field type and grid resolution (bmy, 8/25/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now get the box area at 30N for MFLUX, PRECON (lth, bmy, 5/10/06)
!  (5 ) Rename OTDSCALE to OTD_REG_REDIST.  Also add similar array 
!        OTD_LOC_REDIST.  Now call GET_FLASH_SCALE_CTH, GET_FLASH_SCALE_MFLUX,
!        GET_FLASH_SCALE_PRECON depending on the type of lightning param used.
!        Updated comments.  (ltm, bmy, 1/31/07)
!  (6 ) Removed near-land stuff.  Renamed from INIT_LIGHTNING_NOX_NL to
!        INIT_LIGHTNING_NOX.  Now allocate EMIS_LI_NOx. (ltm, bmy, 10/3/07)
!  (7 ) Also update location of PDF file to lightning_NOx_200709 directory. 
!        (bmy, 1/24/08)
!  (8 ) Read in new Ott profiles from lightning_NOx_201101. Remove
!        depreciated options. (ltm, bmy, 1/25/11)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Removed reference to GET_YEDGE
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: AS, III, IOS, JJJ, IU_FILE
      REAL*8              :: Y0, Y1
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! INIT_LIGHTNING_NOX begins here!
      !=================================================================

      !------------------
      ! Define variables
      !------------------

      ! Get scaling factor to match annual average global flash rate
      ! (ltm, 09/24/07)
      OTD_LIS_SCALE = GET_OTD_LIS_SCALE()

      ! NNLIGHT is the number of points for the lightning CDF's
      NNLIGHT = 3200

      !-----------------
      ! Allocate arrays
      !-----------------

      ! Allocate EMIS_LI_NOx
      ALLOCATE( EMIS_LI_NOx( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMIS_LI_NOx' )
      EMIS_LI_NOx = 0d0

      ! Allocate PROFILE
      ALLOCATE( PROFILE( NNLIGHT, NLTYPE ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROFILE' )
      PROFILE = 0d0

      ! Allocate SLBASE
      ALLOCATE( SLBASE( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SLBASE' )
      SLBASE = 0d0

      IF ( LOTDLOC ) THEN

         ! Array for OTD-LIS local redistribution factors
         ALLOCATE( OTD_LOC_REDIST( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'OTD_LOC_REDIST' )
         OTD_LOC_REDIST = 0d0

      ENDIF

      !=================================================================
      ! Read lightning CDF from Ott et al [JGR, 2010]. (ltm, 1/25/11)
      !=================================================================

      ! Define filename
      FILENAME = 'lightning_NOx_201101/light_dist.ott2010.dat'
      FILENAME = TRIM( DATA_DIR ) // TRIM( FILENAME ) 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - INIT_LIGHTNING: Reading ', a )

      ! Find a free file LUN
      IU_FILE = findFreeLUN()
      
      ! Open file containing lightning PDF data
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'lightdist:1' )
         
      ! Read 12 header lines
      DO III = 1, 12
         READ( IU_FILE, '(a)', IOSTAT=IOS ) 
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'lightdist:2' )
      ENDDO
         
      ! Read NNLIGHT types of lightning profiles
      DO III = 1, NNLIGHT
         READ( IU_FILE,*,IOSTAT=IOS) (PROFILE(III,JJJ),JJJ=1,NLTYPE)
      ENDDO
         
      ! Close file
      CLOSE( IU_FILE )

      END SUBROUTINE INIT_LIGHTNING_NOx
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_lightning_NOx
!
! !DESCRIPTION: Subroutine CLEANUP\_LIGHTNING\_NOx deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_LIGHTNING_NOx
! 
! !REVISION HISTORY: 
!  14 Apr 2004 - R. Yantosca - Initial version
!  (1 ) Now deallocates OTDSCALE (ltm, bmy, 5/10/06)
!  (2 ) Rename OTDSCALE to OTD_REG_REDIST.  Now deallocate OTD_LOC_REDIST.
!        (bmy, 1/31/07)
!  (3 ) Renamed from CLEANUP_LIGHTNING_NOX_NL to CLEANUP_LIGHTNING_NOX.
!        Now deallocate EMIS_LI_NOx. (ltm, bmy, 10/3/07)
!  (4 ) Remove depreciated options. (ltm, bmy, 1/25/11)
!  10 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_LIGHTNING_NOX begins here!
      !=================================================================
      IF ( ALLOCATED( EMIS_LI_NOx    ) ) DEALLOCATE( EMIS_LI_NOx    )
      IF ( ALLOCATED( PROFILE        ) ) DEALLOCATE( PROFILE        )
      IF ( ALLOCATED( SLBASE         ) ) DEALLOCATE( SLBASE         )
      IF ( ALLOCATED( OTD_LOC_REDIST ) ) DEALLOCATE( OTD_LOC_REDIST )

      END SUBROUTINE CLEANUP_LIGHTNING_NOx
!EOC
      END MODULE LIGHTNING_NOx_MOD
