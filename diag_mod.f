! $Id: diag_mod.f,v 1.4 2004/03/24 20:52:29 bmy Exp $
      MODULE DIAG_MOD 
!
!******************************************************************************
!  Module DIAG_MOD contains declarations for allocatable arrays for use with 
!  GEOS-CHEM diagnostics. (amf, bdf, bmy, 11/30/99, 8/20/03)
!
!  Module Routines:
!  ============================================================================
!  (1 ) CLEANUP_DIAG : Deallocates all module arrays  
!
!  GEOS-CHEM modules referenced by diag_mod.f
!  ============================================================================
!  none
!
!  NOTES:
!  (1 ) DIAG_MOD is written in Fixed-Format F90.
!  (2 ) Call subroutine CLEANUP at the end of the MAIN program to deallocate
!        the memory before the run stops.  It is always good style to free
!        any memory we have dynamically allocated when we don't need it
!        anymore
!  (3 ) Added ND13 arrays for sulfur emissions (bmy, 6/6/00)
!  (4 ) Moved ND51 arrays to "diag51_mod.f" (bmy, 11/29/00)
!  (5 ) Added AD34 array for biofuel burning emissions (bmy, 3/15/01)
!  (6 ) Eliminated old commented-out code (bmy, 4/20/01)
!  (7 ) Added AD12 array for boundary layer emissions in routine "setemis.f".
!        (bdf, bmy, 6/15/01)
!  (8 ) Added CHEML24, DRYDL24, CTCHDD for archiving daily mean chemical
!        and drydep loss in chemo3 and chemo3.f (amf, bmy, 7/2/01)
!  (9 ) Add ND43 arrays LTNO2, CTNO2, LTHO2, CTHO2 (rvm, bmy, 2/27/02)
!  (10) Add AD01, AD02 arrays for Rn-Pb-Be simulation (hyl, bmy, 8/7/02)
!  (11) Add AD05 array for sulfate P-L diagnostic (rjp, bdf, bmy, 9/20/02)
!  (12) Added subroutine CLEANUP_DIAG...moved code here from "cleanup.f", 
!        so that it is internal to "diag_mod.f".  Added arrays AD13_NH3_bb,
!        AD13_NH3_bf, AD13_NH3_an for NH3 emissons in ND13.  Deleted obsolete
!        allocatable arrays CHEML24, DRYDL24, CTCHDD.  Now also added LTNO3
!        and CTNO3 arrays for ND43 diagnostic.  Added AD13_SO2_bf array for
!        SO2 biofuel. (bmy, 1/16/03)
!  (13) Added array AD13_NH3_na for ND13 diagnostic (rjp, bmy, 3/23/03)
!  (14) Removed P24H and L24H -- these are now defined w/in "tagged_ox_mod.f"
!        Also added AD03 array for Kr85 prod/loss diag. (jsw, bmy, 8/20/03)
!******************************************************************************
!     
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! For ND01 -- Rn, Pb, Be emissions
      REAL*4,  ALLOCATABLE :: AD01(:,:,:,:)

      ! For ND02 -- Rn, Pb, Be decay
      REAL*4,  ALLOCATABLE :: AD02(:,:,:,:)

      ! For ND03 -- Kr85 prod/loss
      REAL*4,  ALLOCATABLE :: AD03(:,:,:,:)

      ! For ND05 -- Sulfate prod/loss diagnostics
      REAL*4,  ALLOCATABLE :: AD05(:,:,:,:)

      ! For ND12 -- boundary layer multiplication factor
      REAL*4,  ALLOCATABLE :: AD11(:,:,:)

      ! For ND12 -- boundary layer multiplication factor
      REAL*4,  ALLOCATABLE :: AD12(:,:,:)

      ! For ND13 -- Sulfur emissions
      REAL*4,  ALLOCATABLE :: AD13_DMS(:,:)
      REAL*4,  ALLOCATABLE :: AD13_SO2_ac(:,:,:)
      REAL*4,  ALLOCATABLE :: AD13_SO2_an(:,:,:)
      REAL*4,  ALLOCATABLE :: AD13_SO2_bb(:,:)
      REAL*4,  ALLOCATABLE :: AD13_SO2_bf(:,:)
      REAL*4,  ALLOCATABLE :: AD13_SO2_nv(:,:,:)
      REAL*4,  ALLOCATABLE :: AD13_SO2_ev(:,:,:)
      REAL*4,  ALLOCATABLE :: AD13_SO4_an(:,:,:)
      REAL*4,  ALLOCATABLE :: AD13_NH3_an(:,:)
      REAL*4,  ALLOCATABLE :: AD13_NH3_na(:,:)
      REAL*4,  ALLOCATABLE :: AD13_NH3_bb(:,:)
      REAL*4,  ALLOCATABLE :: AD13_NH3_bf(:,:)

      ! For ND14 -- wet convection mass flux diagnostic
      REAL*8,  ALLOCATABLE :: CONVFLUP(:,:,:,:)

      ! For ND15 -- BL mixing mass flux diagnostic
      REAL*8,  ALLOCATABLE :: TURBFLUP(:,:,:,:)

      ! For ND16 -- Fraction of grid box that is precipitating
      REAL*4,  ALLOCATABLE :: AD16(:,:,:,:)  
      INTEGER, ALLOCATABLE :: CT16(:,:,:,:)
      
      ! For ND17 -- Fraction of tracer lost to rainout 
      REAL*4,  ALLOCATABLE :: AD17(:,:,:,:,:)   
      INTEGER, ALLOCATABLE :: CT17(:,:,:,:)

      ! For ND18 -- Fraction of tracer lost to washout
      REAL*4,  ALLOCATABLE :: AD18(:,:,:,:,:)   
      INTEGER, ALLOCATABLE :: CT18(:,:,:,:)

      ! For ND20 -- saves P,L rates for Ox
      REAL*8,  ALLOCATABLE :: PL24H(:,:,:,:)

      ! For ND21 -- Optical Depth diagnostic
      REAL*4,  ALLOCATABLE :: AD21(:,:,:,:)

      ! For ND22 -- J-value diagnostic
      REAL*4,  ALLOCATABLE :: AD22(:,:,:,:)      
      INTEGER, ALLOCATABLE :: LTJV(:,:)
      INTEGER, ALLOCATABLE :: CTJV(:,:) 

      ! For ND23 -- CH3CCl3 lifetime diagnostic
      REAL*8,  ALLOCATABLE :: DIAGCHLORO(:,:,:,:)

      ! For ND24 -- E/W transport mass flux diagnostic
      REAL*8,  ALLOCATABLE :: MASSFLEW(:,:,:,:)

      ! For ND25 -- N/S transport mass flux diagnostic
      REAL*8,  ALLOCATABLE :: MASSFLNS(:,:,:,:)

      ! For ND26 -- UP/DOWN transport mass flux diagnostic
      REAL*8,  ALLOCATABLE :: MASSFLUP(:,:,:,:)

      ! For ND28 -- Biomass burning diagnostic
      REAL*4,  ALLOCATABLE :: AD28(:,:,:)

      ! For ND29 -- CO source diagnostic
      REAL*4,  ALLOCATABLE :: AD29(:,:,:)

      ! For ND31 -- surface pressures
      REAL*4,  ALLOCATABLE :: AD31(:,:,:)

      ! For ND32 -- NOx sources 
      REAL*4,  ALLOCATABLE :: AD32_ac(:,:,:)
      REAL*4,  ALLOCATABLE :: AD32_an(:,:,:)
      REAL*4,  ALLOCATABLE :: AD32_bb(:,:)
      REAL*4,  ALLOCATABLE :: AD32_bf(:,:)
      REAL*4,  ALLOCATABLE :: AD32_fe(:,:)
      REAL*4,  ALLOCATABLE :: AD32_li(:,:,:)
      REAL*4,  ALLOCATABLE :: AD32_so(:,:)
      REAL*4,  ALLOCATABLE :: AD32_ub(:,:)

      ! For ND33 -- tropopsheric sum of tracer
      REAL*4,  ALLOCATABLE :: AD33(:,:,:)

      ! For ND34 -- biofuel emissions
      REAL*4,  ALLOCATABLE :: AD34(:,:,:)

      ! For ND35 -- 500 mb tracer
      REAL*4,  ALLOCATABLE :: AD35(:,:,:)

      ! For ND36 -- Anthropogenic source diagnostic
      REAL*4,  ALLOCATABLE :: AD36(:,:,:)

      ! For ND37 -- Fraction of tracer scavenged in cloud updrafts
      REAL*4,  ALLOCATABLE :: AD37(:,:,:,:)      

      ! For ND38 -- Rainout in moist convection diagnostic
      REAL*4,  ALLOCATABLE :: AD38(:,:,:,:)      

      ! For ND39 -- Washout in aerosol wet deposition diagnostic
      REAL*4,  ALLOCATABLE :: AD39(:,:,:,:)      

      ! For ND41 -- afternoon PBL depths
      REAL*4,  ALLOCATABLE :: AD41(:,:,:)
      INTEGER, ALLOCATABLE :: AFTTOT(:,:)

      ! For ND43 -- OH, NO, NO2, HO2 chemical diagnostics
      REAL*4,  ALLOCATABLE :: AD43(:,:,:,:)      
      INTEGER, ALLOCATABLE :: LTNO(:,:)
      INTEGER, ALLOCATABLE :: CTNO(:,:)
      INTEGER, ALLOCATABLE :: LTOH(:,:)
      INTEGER, ALLOCATABLE :: CTOH(:,:)
      INTEGER, ALLOCATABLE :: LTNO2(:,:)
      INTEGER, ALLOCATABLE :: CTNO2(:,:)
      INTEGER, ALLOCATABLE :: LTHO2(:,:)
      INTEGER, ALLOCATABLE :: CTHO2(:,:)
      INTEGER, ALLOCATABLE :: LTNO3(:,:)
      INTEGER, ALLOCATABLE :: CTNO3(:,:)

      ! For ND44 -- Dry deposition fluxes & velocities
      REAL*4,  ALLOCATABLE :: AD44(:,:,:,:)

      ! For ND45 -- Tracer concentration diagnostic
      REAL*4,  ALLOCATABLE :: AD45(:,:,:,:)      
      INTEGER, ALLOCATABLE :: LTOTH(:,:)
      INTEGER, ALLOCATABLE :: CTOTH(:,:)

      ! For ND45 -- Tracer concentration diagnostic
      REAL*4,  ALLOCATABLE :: AD46(:,:,:)      

      ! For ND47 -- 24-h tracer concentration diagnostic
      REAL*4,  ALLOCATABLE :: AD47(:,:,:,:)      

      ! For ND48 -- station timeseries diagnostic
      REAL*8,  ALLOCATABLE :: TCOBOX(:,:)

      ! Dynamically allocatable array -- local only to DIAG50.F
      REAL*8,  ALLOCATABLE :: STT_TEMPO2(:,:,:,:)

      ! For ND55 -- tropopause diagnostics
      REAL*4,  ALLOCATABLE :: AD55(:,:,:)
      
      ! For ND65 -- Chemical family P-L diagnostic
      REAL*4,  ALLOCATABLE :: AD65(:,:,:,:)
      REAL*8,  ALLOCATABLE :: FAMPL(:,:,:,:)

      ! For ND66 -- I-6 fields diagnostic
      REAL*4,  ALLOCATABLE :: AD66(:,:,:,:)      

      ! For ND67 -- DAO surface fields diagnostic
      REAL*4,  ALLOCATABLE :: AD67(:,:,:)      

      ! For ND68 -- BXHEIGHT, AD, AVGW diagnostic
      REAL*4,  ALLOCATABLE :: AD68(:,:,:,:)      

      ! For ND69 -- DXYP diagnostic
      REAL*4,  ALLOCATABLE :: AD69(:,:,:)      

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG deallocates all module arrays.
!  (bmy, 12/13/02, 8/20/03)
!
!  NOTES:
!  (1 ) Now also deallocate AD13_NH3_an, AD13_NH3_bb, AD13_NH3_bf arrays
!        for the ND13 diagnostic.  (bmy, 12/13/02)
!  (2 ) Now also deallocate AD13_NH3_na array for ND13 (rjp, bmy, 3/23/03)
!  (3 ) Removed P24H and L24H, these are now defined within "tagged_ox_mod.f".
!       Now also deallocate AD03 array for Kr85 prod/loss (jsw, bmy, 8/20/03)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG begins here!
      !=================================================================
      IF ( ALLOCATED( AD01        ) ) DEALLOCATE( AD01        )
      IF ( ALLOCATED( AD02        ) ) DEALLOCATE( AD02        )
      IF ( ALLOCATED( AD03        ) ) DEALLOCATE( AD03        )
      IF ( ALLOCATED( AD11        ) ) DEALLOCATE( AD11        )
      IF ( ALLOCATED( AD12        ) ) DEALLOCATE( AD12        )
      IF ( ALLOCATED( AD13_DMS    ) ) DEALLOCATE( AD13_DMS    )
      IF ( ALLOCATED( AD13_SO2_ac ) ) DEALLOCATE( AD13_SO2_ac )
      IF ( ALLOCATED( AD13_SO2_an ) ) DEALLOCATE( AD13_SO2_an )
      IF ( ALLOCATED( AD13_SO2_bb ) ) DEALLOCATE( AD13_SO2_bb )
      IF ( ALLOCATED( AD13_SO2_bf ) ) DEALLOCATE( AD13_SO2_bf )
      IF ( ALLOCATED( AD13_SO2_nv ) ) DEALLOCATE( AD13_SO2_nv )
      IF ( ALLOCATED( AD13_SO2_ev ) ) DEALLOCATE( AD13_SO2_ev )
      IF ( ALLOCATED( AD13_SO4_an ) ) DEALLOCATE( AD13_SO4_an )
      IF ( ALLOCATED( AD13_NH3_an ) ) DEALLOCATE( AD13_NH3_an )
      IF ( ALLOCATED( AD13_NH3_na ) ) DEALLOCATE( AD13_NH3_na )
      IF ( ALLOCATED( AD13_NH3_bb ) ) DEALLOCATE( AD13_NH3_bb )
      IF ( ALLOCATED( AD13_NH3_bf ) ) DEALLOCATE( AD13_NH3_bf )
      IF ( ALLOCATED( AD16        ) ) DEALLOCATE( AD16        )
      IF ( ALLOCATED( AD17        ) ) DEALLOCATE( AD17        )
      IF ( ALLOCATED( AD18        ) ) DEALLOCATE( AD18        )
      IF ( ALLOCATED( AD21        ) ) DEALLOCATE( AD21        )
      IF ( ALLOCATED( AD22        ) ) DEALLOCATE( AD22        ) 
      IF ( ALLOCATED( AD28        ) ) DEALLOCATE( AD28        ) 
      IF ( ALLOCATED( AD29        ) ) DEALLOCATE( AD29        ) 
      IF ( ALLOCATED( AD31        ) ) DEALLOCATE( AD31        ) 
      IF ( ALLOCATED( AD32_ac     ) ) DEALLOCATE( AD32_ac     ) 
      IF ( ALLOCATED( AD32_an     ) ) DEALLOCATE( AD32_an     )
      IF ( ALLOCATED( AD32_bb     ) ) DEALLOCATE( AD32_bb     )
      IF ( ALLOCATED( AD32_bf     ) ) DEALLOCATE( AD32_bf     )
      IF ( ALLOCATED( AD32_fe     ) ) DEALLOCATE( AD32_fe     )
      IF ( ALLOCATED( AD32_li     ) ) DEALLOCATE( AD32_li     )
      IF ( ALLOCATED( AD32_so     ) ) DEALLOCATE( AD32_so     )
      IF ( ALLOCATED( AD32_ub     ) ) DEALLOCATE( AD32_ub     )
      IF ( ALLOCATED( AD33        ) ) DEALLOCATE( AD33        )
      IF ( ALLOCATED( AD34        ) ) DEALLOCATE( AD34        )
      IF ( ALLOCATED( AD35        ) ) DEALLOCATE( AD35        )
      IF ( ALLOCATED( AD36        ) ) DEALLOCATE( AD36        )
      IF ( ALLOCATED( AD37        ) ) DEALLOCATE( AD37        )
      IF ( ALLOCATED( AD38        ) ) DEALLOCATE( AD38        )  
      IF ( ALLOCATED( AD39        ) ) DEALLOCATE( AD39        )
      IF ( ALLOCATED( AD41        ) ) DEALLOCATE( AD41        )
      IF ( ALLOCATED( AD43        ) ) DEALLOCATE( AD43        )
      IF ( ALLOCATED( AD44        ) ) DEALLOCATE( AD44        )
      IF ( ALLOCATED( AD45        ) ) DEALLOCATE( AD45        )
      IF ( ALLOCATED( AD46        ) ) DEALLOCATE( AD46        )
      IF ( ALLOCATED( AD47        ) ) DEALLOCATE( AD47        )
      IF ( ALLOCATED( AD55        ) ) DEALLOCATE( AD55        )
      IF ( ALLOCATED( AD65        ) ) DEALLOCATE( AD65        )
      IF ( ALLOCATED( AD66        ) ) DEALLOCATE( AD66        )
      IF ( ALLOCATED( AD68        ) ) DEALLOCATE( AD68        )
      IF ( ALLOCATED( AD69        ) ) DEALLOCATE( AD69        )
      IF ( ALLOCATED( AFTTOT      ) ) DEALLOCATE( AFTTOT      )
      IF ( ALLOCATED( CONVFLUP    ) ) DEALLOCATE( CONVFLUP    )
      IF ( ALLOCATED( CT16        ) ) DEALLOCATE( CT16        )
      IF ( ALLOCATED( CT17        ) ) DEALLOCATE( CT17        )
      IF ( ALLOCATED( CT18        ) ) DEALLOCATE( CT18        )
      IF ( ALLOCATED( CTJV        ) ) DEALLOCATE( CTJV        )
      IF ( ALLOCATED( CTNO        ) ) DEALLOCATE( CTNO        )
      IF ( ALLOCATED( CTOH        ) ) DEALLOCATE( CTOH        )
      IF ( ALLOCATED( CTNO2       ) ) DEALLOCATE( CTNO2       )
      IF ( ALLOCATED( CTNO3       ) ) DEALLOCATE( CTNO3       )
      IF ( ALLOCATED( CTHO2       ) ) DEALLOCATE( CTHO2       )
      IF ( ALLOCATED( CTOTH       ) ) DEALLOCATE( CTOTH       )
      IF ( ALLOCATED( DIAGCHLORO  ) ) DEALLOCATE( DIAGCHLORO  )
      IF ( ALLOCATED( FAMPL       ) ) DEALLOCATE( FAMPL       )
      IF ( ALLOCATED( LTJV        ) ) DEALLOCATE( LTJV        )
      IF ( ALLOCATED( LTNO        ) ) DEALLOCATE( LTNO        )
      IF ( ALLOCATED( LTOH        ) ) DEALLOCATE( LTOH        )
      IF ( ALLOCATED( LTNO2       ) ) DEALLOCATE( LTNO2       )
      IF ( ALLOCATED( LTNO3       ) ) DEALLOCATE( LTNO3       )
      IF ( ALLOCATED( LTHO2       ) ) DEALLOCATE( LTHO2       )
      IF ( ALLOCATED( LTOTH       ) ) DEALLOCATE( LTOTH       )
      IF ( ALLOCATED( MASSFLEW    ) ) DEALLOCATE( MASSFLEW    )
      IF ( ALLOCATED( MASSFLNS    ) ) DEALLOCATE( MASSFLNS    )
      IF ( ALLOCATED( MASSFLUP    ) ) DEALLOCATE( MASSFLUP    )
      IF ( ALLOCATED( PL24H       ) ) DEALLOCATE( PL24H       )
      IF ( ALLOCATED( TCOBOX      ) ) DEALLOCATE( TCOBOX      )
      IF ( ALLOCATED( TURBFLUP    ) ) DEALLOCATE( TURBFLUP    )
      IF ( ALLOCATED( STT_TEMPO2  ) ) DEALLOCATE( STT_TEMPO2  ) 

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG
      
!------------------------------------------------------------------------------

      END MODULE DIAG_MOD 

