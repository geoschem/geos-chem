!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: canopy_nox_mod
!     
! !DESCRIPTION: Module CANOPY\_NOX\_MOD contains routines for computing the
!  bulk surface resistance of the canopy to NOx.  This quantity is needed by 
!  GEOS-Chem soil emissions routine "soilcrf.f".
!\\   
!\\   
! !INTERFACE: 
!     
      MODULE CANOPY_NOX_MOD
!
! !USES
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: GET_CANOPY_NOX
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: DIFFG
!
!  !REMARKS:
!  This computation was originally done in legacy routine DEPVEL (located in 
!  "drydep\_mod.f"), and the bulk surface resistance was stored in common block 
!  array CANOPYNOX (in "commsoil.h").  However, the legacy code caused an ugly 
!  dependency.  Drydep routine DEPVEL had to be called BEFORE the soil NOx 
!  emissions routines.  Routine GET\_CANOPY\_NOX in this module performs the 
!  same computation of NOx from tree canopies independent of "drydep\_mod.f", 
!  thus allowing for a totally clean separation between dry deposition routines
!  and emisisons routines.
! 
!  References:
!  ============================================================================
!  (1 ) Baldocchi, D.D., B.B. Hicks, and P. Camara, "A canopy stomatal
!        resistance model for gaseous deposition to vegetated surfaces",
!        Atmos. Environ. 21, 91-101, 1987.
!  (2 ) Brutsaert, W., "Evaporation into the Atmosphere", Reidel, 1982.
!  (3 ) Businger, J.A., et al., "Flux-profile relationships in the atmospheric 
!        surface layer", J. Atmos. Sci., 28, 181-189, 1971.
!  (4 ) Dwight, H.B., "Tables of integrals and other mathematical data",
!        MacMillan, 1957.
!  (5 ) Guenther, A., and 15 others, A global model of natural volatile
!         organic compound emissions, J. Geophys. Res., 100, 8873-8892, 1995.
!  (6 ) Hicks, B.B., and P.S. Liss, "Transfer of SO2 and other reactive
!        gases across the air-sea interface", Tellus, 28, 348-354, 1976.
!  (7 ) Jacob, D.J., and S.C. Wofsy, "Budgets of reactive nitrogen,
!        hydrocarbons, and ozone over the Amazon forest during the wet season",
!        J.  Geophys. Res., 95, 16737-16754, 1990.
!  (8 ) Jacob, D.J., et al, "Deposition of ozone to tundra", J. Geophys. Res., 
!        97, 16473-16479, 1992.
!  (9 ) Levine, I.N., "Physical Chemistry, 3rd ed.", McGraw-Hill, 
!        New York, 1988.
!  (10) Munger, J.W., et al, "Atmospheric deposition of reactive nitrogen 
!        oxides and ozone in a temperate deciduous forest and a sub-arctic 
!        woodland", J. Geophys. Res., in press, 1996.
!  (11) Walcek, C.J., R.A. Brost, J.S. Chang, and M.L. Wesely, "SO2, sulfate, 
!        and HNO3 deposition velocities computed using regional landuse and
!        meteorological data", Atmos. Environ., 20, 949-964, 1986.
!  (12) Wang, Y.H., paper in preparation, 1996.
!  (13) Wesely, M.L, "Improved parameterizations for surface resistance to
!        gaseous dry deposition in regional-scale numerical models", 
!        Environmental Protection Agency Report EPA/600/3-88/025,
!        Research Triangle Park (NC), 1988.
!  (14) Wesely, M. L., Parameterization of surface resistance to gaseous dry 
!        deposition in regional-scale numerical models.  Atmos. Environ., 23
!        1293-1304, 1989. 
!  (15) Price, H., L. Jaeglé, A. Rice, P. Quay, P.C. Novelli, R. Gammon, 
!        Global Budget of Molecular Hydrogen and its Deuterium Content: 
!        Constraints from Ground Station, Cruise, and Aircraft Observations,
!        submitted to J. Geophys. Res., 2007.
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca - Split off from "drydep_mod.f"
!EOP
!------------------------------------------------------------------------------
!BOC
!
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_canopy_nox
!
! !DESCRIPTION: Subroutine GET\_CANOPY\_NOX computes the bulk surface 
!  resistance of the canopy to NOx.  This computation was originally done 
!  within legacy routine DEPVEL (in "drydep\_mod.f").  Moving this computation 
!  to GET\_CANOPY\_NOX now allows for a totally clean separation between 
!  dry deposition routines and emissions routines in GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_CANOPY_NOX( State_Met )
!
! !USES:
!
      USE DRYDEP_MOD,         ONLY : DRYCOEFF	!todo
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE MODIS_LAI_MOD

      USE COMMSOIL_MOD
      USE CMN_SIZE_MOD
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !REMARKS:
!     For backwards compatibility, the bulk surface resistance is stored
!     in common block array CANOPYNOX in "commsoil.h".  Leave it like this
!     for the time being...we'll clean it up when we fix all of the soil
!     NOx routines.
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca     - Split off from "drydep_mod.f"
!  14 Jun 2012 - J.D. Maasakkers - Rewritten as a function of the
!                                     MODIS/Koppen biometype
!  09 Nov 2012 - M. Payer        - Replaced all met field arrays with State_Met
!                                   derived type object
!  13 Dec 2012 - R. Yantosca     - Removed ref to obsolete CMN_DEP_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !EXTERNAL FUNCTIONS:
!
      REAL*8,  EXTERNAL:: BIOFIT
!
! !DEFINED PARAMETERS:
!
      ! Molecular weight of water [kg]
      REAL*8, PARAMETER :: XMWH2O = 18d-3

      ! Surface pressure??? [Pa]
      REAL*8, PARAMETER :: PRESS  = 1.5d5
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER          :: I,      IJLOOP, J, K, KK
      REAL*8           :: F0,     HSTAR, XMW              
      REAL*8           :: DTMP1,  DTMP2, DTMP3,  DTMP4, GFACT, GFACI
      REAL*8           :: RT,     RAD0,  RIX,    RIXX,  RDC,   RLUXX
      REAL*8           :: RGSX,   RCLX,  TEMPK,  TEMPC

      LOGICAL, SAVE    :: FIRSTCANOPY = .TRUE. 

      ! Arrays
      REAL*8           :: RI  (NSOILB)       
      REAL*8           :: RLU (NSOILB)      
      REAL*8           :: RAC (NSOILB)      
      REAL*8           :: RGSS(NSOILB)     
      REAL*8           :: RGSO(NSOILB)     
      REAL*8           :: RCLS(NSOILB)     
      REAL*8           :: RCLO(NSOILB)  

      !=================================================================
      ! GET_CANOPY_NOX begins here!
      !=================================================================

      ! Set physical parameters
      HSTAR = 0.01d0              ! Henry's law constant
      F0    = 0.1d0               ! Reactivity factor for biological oxidation 
      XMW   = 46d-3               ! Molecular wt of NO2 (kg)

      IF ( FIRSTCANOPY ) THEN
         CALL READ_RESISTANCES
         FIRSTCANOPY = .FALSE.    ! Now only read resistances once
      ENDIF

      ! Loop over surface boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! 1-D grid box index
         IJLOOP = ( (J-1) * IIPAR ) + I

         ! Surface temperature [K] and [C]
         TEMPK = State_Met%TS(I,J)
         TEMPC = State_Met%TS(I,J) - 273.15d0

         ! Compute bulk surface resistance for gases.    
         !                                  
         !  Adjust external surface resistances for temperature; 
         !  from Wesely [1989], expression given in text on p. 1296.        
         RT = 1000.0D0 * EXP( -TEMPC - 4.0d0 )

         !--------------------------------------------------------------
         ! Get surface resistances - loop over biome types K
         !
         ! The land types within each grid square are defined using the 
         ! Olson land-type database.  Each of the Olson land types is 
         ! assigned a corresponding "deposition land type" with 
         ! characteristic values of surface resistance components.  
         ! There are 74 Olson land-types but only 11 deposition 
         ! land-types (i.e., many of the Olson land types share the 
         ! same deposition characteristics).  Surface resistance 
         ! components for the "deposition land types" are from Wesely 
         ! [1989] except for tropical forests [Jacob and Wofsy, 1990] 
         ! and for tundra [Jacob et al., 1992].  All surface resistance 
         ! components are normalized to a leaf area index of unity.
         !--------------------------------------------------------------
 	 !Loop over all biometypes
         DO K = 1, 24

            ! Skip if not present
            IF ( LAND2(I,J,K) == 0 ) CYCLE

            ! Set second loop variable to K to allow snow/ice correction
	    KK = K

            ! If the surface is snow or ice, then set K=3
            IF ( State_Met%ALBD(I,J) > 0.4 ) KK = 3

            ! USE new MODIS/KOPPEN Biometypes to read data

            ! Read the internal resistance RI (minimum stomatal resistance 
            ! for water vapor, per unit area of leaf) from the IRI array; 
            ! a '9999' value means no deposition to stomata so we impose a 
            ! very large value for RI.
            RI(K) = DBLE( SNIRI(KK) )
            IF ( RI(K) >= 9999.D0 ) RI(K)= 1.D12

            ! Cuticular resistances IRLU read in from 'drydep.table'
            ! are per unit area of leaf; divide them by the leaf area index 
            ! to get a cuticular resistance for the bulk canopy.  If IRLU is 
            !'9999' it means there are no cuticular surfaces on which to 
            ! deposit so we impose a very large value for RLU.
            IF ( SNIRLU(KK) >= 9999 .OR. GC_LAI(I,J) <= 0d0 ) THEN
               RLU(K)  = 1.D6
            ELSE
               RLU(K)= DBLE( SNIRLU(KK) ) / GC_LAI(I,J) + RT
            ENDIF

            ! The following are the remaining resistances for the Wesely
            ! resistance-in-series model for a surface canopy
            ! (see Atmos. Environ. paper, Fig.1).  
            RAC(K)  = MAX( DBLE( SNIRAC(KK)  ),      1d0 )
            RGSS(K) = MAX( DBLE( SNIRGSS(KK) ) + RT, 1d0 )
            RGSO(K) = MAX( DBLE( SNIRGSO(KK) ) + RT, 1d0 ) 
            RCLS(K) =      DBLE( SNIRCLS(KK) ) + RT           
            RCLO(K) =      DBLE( SNIRCLO(KK) ) + RT 

            IF (  RAC(K) >= 9999.D0 ) RAC(K)  = 1d12
            IF ( RGSS(K) >= 9999.D0 ) RGSS(K) = 1d12
            IF ( RGSO(K) >= 9999.D0 ) RGSO(K) = 1d12
            IF ( RCLS(K) >= 9999.D0 ) RCLS(K) = 1d12         
            IF ( RCLO(K) >= 9999.D0 ) RCLO(K) = 1d12

            !-------------------------------------------------------------
            ! Adjust stomatal resistances for insolation and temperature:  
            ! 
            ! Temperature adjustment is from Wesely [1989], equation (3).
            ! 
            ! Light adjustment by the function BIOFIT is described by Wang 
            ! [1996].  It combines:
            !
            ! - Local dependence of stomal resistance on the intensity I 
            !   of light impinging the leaf; this is expressed as a 
            !   multiplicative factor I/(I+b) to the stomatal resistance 
            !   where b = 50 W m-2
            !   (equation (7) of Baldocchi et al. [1987])
            ! - Radiative transfer of direct and diffuse radiation in the 
            !   canopy using equations (12)-(16) from Guenther et al. 
            !   [1995]
            ! - Separate accounting of sunlit and shaded leaves using
            !   equation (12) of Guenther et al. [1995]
            ! - Partitioning of the radiation at the top of the canopy 
            !   into direct and diffuse components using a 
            !   parameterization to results from an atmospheric radiative 
            !   transfer model [Wang, 1996]
            !
            ! The dependent variables of the function BIOFIT are the leaf 
            ! area index (XYLAI), the cosine of zenith angle (SUNCOS) and 
            ! the fractional cloud cover (CFRAC).  The factor GFACI 
            ! integrates the light dependence over the canopy depth; so
            ! be scaled by LAI to yield a bulk canopy value because that's 
            ! already done in the GFACI formulation.
            !-------------------------------------------------------------

            ! Radiation @ sfc [W/m2]
            RAD0 = State_Met%RADSWG(I,J)
            
            ! Internal resistance
            RIX  = RI(K)

            ! Skip the following block if the resistance RIX is high
            IF ( RIX < 9999d0 ) THEN
               GFACT = 100.0D0

               IF ( TEMPC > 0.D0 .AND. TEMPC < 40.D0) THEN
                  GFACT = 400.D0 / TEMPC / ( 40.0D0 - TEMPC )
               ENDIF

               GFACI = 100.D0

               IF ( RAD0 > 0d0 .AND. GC_LAI(I,J) > 0d0 ) THEN
                  GFACI= 1d0 / 
     &                   BIOFIT( DRYCOEFF,       
     &                           GC_LAI(I,J),
     &                           State_Met%SUNCOS(I,J), 
     &                           State_Met%CLDFRC(I,J)  )
               ENDIF
            
               RIX = RIX * GFACT * GFACI
            ENDIF

            ! Compute aerodynamic resistance to lower elements in lower 
            ! part of the canopy or structure, assuming level terrain - 
            ! equation (5) of Wesely [1989].                     
            RDC = 100.D0*(1.0D0+1000.0D0/(RAD0 + 10.D0))

            ! Loop over species; species-dependent corrections to resistances
            ! are from equations (6)-(9) of Wesely [1989].
            !
            ! NOTE: here we only consider NO2 (bmy, 6/22/09)
            RIXX   = RIX * DIFFG( TEMPK, PRESS, XMWH2O ) /
     &                     DIFFG( TEMPK, PRESS, XMW    )
     &             + 1.D0 / ( HSTAR/3000.D0 + 100.D0*F0  )

            RLUXX  = 1.D12

            IF ( RLU(K) < 9999.D0 ) THEN
               RLUXX = RLU(K) / ( HSTAR / 1.0D+05 + F0 )
            ENDIF

            ! To prevent virtually zero resistance to species with huge HSTAR, 
            ! such as HNO3, a minimum value of RLUXX needs to be set. 
            ! The rationality of the existence of such a minimum is 
            ! demonstrated by the observed relationship between Vd(NOy-NOx) 
            ! and Ustar in Munger et al.[1996]; Vd(HNO3) never exceeds 2 cm/s 
            ! in observations. The corresponding minimum resistance is 50 s/m.
            ! was introduced by J.Y. Liang on 7/9/95.
            RGSX = 1d0 / ( HSTAR/1d5/RGSS(K) + F0/RGSO(K) )
            RCLX = 1d0 / ( HSTAR/1d5/RCLS(K) + F0/RCLO(K) )

            ! Get the bulk surface resistance of the canopy
            ! from the network of resistances in parallel and in series 
            ! (Fig. 1 of Wesely [1989])
            DTMP1 = 1.D0 / RIXX
            DTMP2 = 1.D0 / RLUXX
            DTMP3 = 1.D0 / ( RAC(K) + RGSX )
            DTMP4 = 1.D0 / ( RDC      + RCLX )

            ! Save the within canopy depvel of NOx, used in calculating 
            ! the canopy reduction factor for soil emissions [1/s]
            CANOPYNOX(IJLOOP,K) = DTMP1 + DTMP2 + DTMP3 + DTMP4

         ENDDO
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE GET_CANOPY_NOx
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diffg
!
! !DESCRIPTION: Function DIFFG calculates the molecular diffusivity [m2/s] in 
!  air for a gas X of molecular weight XM [kg] at temperature TK [K] and 
!  pressure PRESS [Pa].
!\\
!\\
! !INTERFACE:
!
      FUNCTION DIFFG( TK, PRESS, XM ) RESULT( DIFF_G )
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN) :: TK      ! Temperature [K]
      REAL*8, INTENT(IN) :: PRESS   ! Pressure [hPa]
      REAL*8, INTENT(IN) :: XM      ! Molecular weight of gas [kg]
!
! !RETURN VALUE:
!
      REAL*8             :: DIFF_G  ! Molecular diffusivity [m2/s]
!
! !REMARKS:
!  We specify the molecular weight of air (XMAIR) and the hard-sphere molecular
!  radii of air (RADAIR) and of the diffusing gas (RADX).  The molecular
!  radius of air is given in a Table on p. 479 of Levine [1988].  The Table
!  also gives radii for some other molecules.  Rather than requesting the user
!  to supply a molecular radius we specify here a generic value of 2.E-10 m for
!  all molecules, which is good enough in terms of calculating the diffusivity
!  as long as molecule is not too big.
! 
! !REVISION HISTORY:
!     22 Jun 2009 - R. Yantosca - Copied from "drydep_mod.f"
!EOP
!------------------------------------------------------------------------------
!BOC      
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER  :: XMAIR  = 28.8d-3 
      REAL*8, PARAMETER  :: RADAIR = 1.2d-10
      REAL*8, PARAMETER  :: PI     = 3.1415926535897932d0
      REAL*8, PARAMETER  :: RADX   = 1.5d-10
      REAL*8, PARAMETER  :: RGAS   = 8.32d0
      REAL*8, PARAMETER  :: AVOGAD = 6.023d23
!
! !LOCAL VARIABLES:
!
      REAL*8             :: AIRDEN, Z, DIAM, FRPATH, SPEED      

      !=================================================================
      ! DIFFG begins here!
      !=================================================================

      ! Air density
      AIRDEN = ( PRESS * AVOGAD ) / ( RGAS * TK )

      ! DIAM is the collision diameter for gas X with air.
      DIAM   = RADX + RADAIR

      ! Calculate the mean free path for gas X in air: 
      ! eq. 8.5 of Seinfeld [1986];
      Z      = XM  / XMAIR
      FRPATH = 1d0 /( PI * SQRT( 1d0 + Z ) * AIRDEN*( DIAM**2 ) )

      ! Calculate average speed of gas X; eq. 15.47 of Levine [1988]
      SPEED  = SQRT( 8d0 * RGAS * TK / ( PI * XM ) )

      ! Calculate diffusion coefficient of gas X in air; 
      ! eq. 8.9 of Seinfeld [1986]
      DIFF_G = ( 3d0 * PI / 32d0 ) * ( 1d0 + Z ) * FRPATH * SPEED

      ! Return to calling program
      END FUNCTION DIFFG
!EOC
      END MODULE CANOPY_NOX_MOD
