!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag_mod.F
!
! !DESCRIPTION: Module DIAG\_MOD contains declarations for allocatable arrays 
!  for use with GEOS-CHEM diagnostics. 
!\\
!\\
! !INTERFACE: 
!
      MODULE DIAG_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_DIAG
!
! !PUBLIC DATA MEMBERS:
!
#if defined( BPCH_DIAG ) || defined( MODEL_GEOS )

      ! For ND01 -- Pb emissions (Rn, Be now tracked thru HEMCO)
      REAL*4,  ALLOCATABLE :: AD01(:,:,:)

      ! For ND02 -- Rn, Pb, Be decay
      REAL*4,  ALLOCATABLE :: AD02(:,:,:,:)

      ! For ND05 -- Sulfate prod/loss diagnostics
      REAL*4,  ALLOCATABLE :: AD05(:,:,:,:)

      ! For ND07 -- Carbon aerosol emission
      REAL*4,  ALLOCATABLE :: AD07_BC(:,:,:)
      REAL*4,  ALLOCATABLE :: AD07_OC(:,:,:)
      REAL*4,  ALLOCATABLE :: AD07_HC(:,:,:,:)

      ! For ND09 -- HCN / CH3CN simulation
      REAL*4,  ALLOCATABLE :: AD09(:,:,:,:)
      REAL*4,  ALLOCATABLE :: AD09_em(:,:,:)

      ! For ND12 -- boundary layer multiplication factor
      REAL*4,  ALLOCATABLE :: AD12(:,:,:)

      ! For ND14 -- wet convection mass flux diagnostic
      REAL(fp),  ALLOCATABLE :: CONVFLUP(:,:,:,:)

      ! For ND15 -- BL mixing mass flux diagnostic
      REAL(fp),  ALLOCATABLE :: TURBFLUP(:,:,:,:)

      ! For ND16 -- Fraction of grid box that is precipitating
      REAL*4,  ALLOCATABLE :: AD16(:,:,:,:)  
      INTEGER, ALLOCATABLE :: CT16(:,:,:,:)
      
      ! For ND17 -- Fraction of tracer lost to rainout 
      REAL*4,  ALLOCATABLE :: AD17(:,:,:,:,:)   
      INTEGER, ALLOCATABLE :: CT17(:,:,:,:)

      ! For ND18 -- Fraction of tracer lost to washout
      REAL*4,  ALLOCATABLE :: AD18(:,:,:,:,:)   
      INTEGER, ALLOCATABLE :: CT18(:,:,:,:)

      ! For ND21 -- Optical Depth diagnostic
      REAL*4,  ALLOCATABLE :: AD21(:,:,:,:)

      ! For ND22 -- J-value diagnostic
      REAL*4,  ALLOCATABLE :: AD22(:,:,:,:)      
      INTEGER, ALLOCATABLE :: LTJV(:,:)
      INTEGER, ALLOCATABLE :: CTJV(:,:) 

      ! For ND24 -- E/W transport mass flux diagnostic
      REAL(fp),  ALLOCATABLE, TARGET :: MASSFLEW(:,:,:,:)

      ! For ND25 -- N/S transport mass flux diagnostic
      REAL(fp),  ALLOCATABLE, TARGET :: MASSFLNS(:,:,:,:)

      ! For ND26 -- UP/DOWN transport mass flux diagnostic
      REAL(fp),  ALLOCATABLE, TARGET :: MASSFLUP(:,:,:,:)

      ! For ND30 -- land / water / ice flags
      REAL*4,  ALLOCATABLE :: AD30(:,:)

      ! For ND31 -- surface pressures
      REAL*4,  ALLOCATABLE :: AD31(:,:,:)

      ! For ND33 -- tropopsheric sum of tracer
      REAL*4,  ALLOCATABLE :: AD33(:,:,:)

      ! For ND35 -- 500 mb tracer
      REAL*4,  ALLOCATABLE :: AD35(:,:,:)

      ! For ND37 -- Fraction of tracer scavenged in cloud updrafts
      REAL*4,  ALLOCATABLE :: AD37(:,:,:,:)      

      ! For ND38 -- Rainout in moist convection diagnostic
      REAL*4,  ALLOCATABLE :: AD38(:,:,:,:)      

      ! For ND39 -- Washout in aerosol wet deposition diagnostic
      REAL*4,  ALLOCATABLE :: AD39(:,:,:,:)      

      ! For ND43 -- OH and HO2 chemical diagnostics
      REAL*4,  ALLOCATABLE :: AD43(:,:,:,:)      
      INTEGER, ALLOCATABLE :: LTOH(:,:)
      INTEGER, ALLOCATABLE :: CTOH(:,:,:)
      INTEGER, ALLOCATABLE :: LTHO2(:,:)
      INTEGER, ALLOCATABLE :: CTHO2(:,:,:)
      ! Added for stratospheric chem. (SDE 2013-07-15)
      INTEGER, ALLOCATABLE :: LTO1D(:,:)
      INTEGER, ALLOCATABLE :: CTO1D(:,:,:)
      INTEGER, ALLOCATABLE :: LTO3P(:,:)
      INTEGER, ALLOCATABLE :: CTO3P(:,:,:)

      ! For ND44 -- Dry deposition fluxes & velocities
      REAL*4,  ALLOCATABLE :: AD44(:,:,:,:)
!------------------------------------------------------------------------------
! This diagnostic requires some fixes. Comment out for now (mps,3/2/18)
!      ! Instantaneous dry dep velocity (krt, 11/2017)
!      REAL*4,  ALLOCATABLE :: AD44b(:,:,:)
!      ! For aerodynamic resistance (krt, 11/2017)
!      REAL*4,  ALLOCATABLE :: AD_RA(:,:)
!      ! For test of M/O theory applicability (krt, 02/2018)
!      REAL*4,  ALLOCATABLE :: AD_ZL(:,:)
!------------------------------------------------------------------------------
      ! For ND45 -- Tracer concentration diagnostic
      REAL*4,  ALLOCATABLE :: AD45(:,:,:,:)
      INTEGER, ALLOCATABLE :: LTOTH(:,:)
      INTEGER, ALLOCATABLE :: CTOTH(:,:)

      ! For ND47 -- 24-h tracer concentration diagnostic
      REAL*4,  ALLOCATABLE :: AD47(:,:,:,:)      

      ! For ND47(O3) / ND65 -- 24-h tracer diagnostic
      INTEGER, ALLOCATABLE :: CTO3_24h(:,:,:)

      ! Dynamically allocatable array -- local only to DIAG50.F
      REAL(fp),  ALLOCATABLE :: STT_TEMPO2(:,:,:,:)

      ! For ND52 -- gamma HO2 diagnostic
      REAL*4,  ALLOCATABLE :: AD52(:,:,:,:)

      ! For ND54 -- tropopause diagnostics
      REAL*4,  ALLOCATABLE :: AD54(:,:,:)

      ! For ND55 -- tropopause diagnostics
      REAL*4,  ALLOCATABLE :: AD55(:,:,:)

      ! For ND57 -- theta, potential temp (FP 6/2009)
      REAL*4,  ALLOCATABLE :: AD57(:,:,:)

      ! -- for methane simulation diagnostics
      REAL*4,  ALLOCATABLE :: AD19(:,:,:,:)
      REAL*4,  ALLOCATABLE :: AD58(:,:,:)
      REAL*4,  ALLOCATABLE :: AD60(:,:)

#if   defined( TOMAS )
      ! For ND59 -- Size-resolved primary aerosol emissions
      REAL*4,  ALLOCATABLE :: AD59_NUMB(:,:,:,:)
      REAL*4,  ALLOCATABLE :: AD59_SULF(:,:,:,:)
      REAL*4,  ALLOCATABLE :: AD59_SALT(:,:,:,:)
      REAL*4,  ALLOCATABLE :: AD59_ECIL(:,:,:,:)
      REAL*4,  ALLOCATABLE :: AD59_ECOB(:,:,:,:)
      REAL*4,  ALLOCATABLE :: AD59_OCIL(:,:,:,:)
      REAL*4,  ALLOCATABLE :: AD59_OCOB(:,:,:,:)
      REAL*4,  ALLOCATABLE :: AD59_DUST(:,:,:,:)

      ! For ND60 -- TOMAS condensation rate diagnostic
      REAL*4,  ALLOCATABLE :: AD60_COND(:,:,:,:)      
      REAL*4,  ALLOCATABLE :: AD60_COAG(:,:,:,:)      
      REAL*4,  ALLOCATABLE :: AD60_NUCL(:,:,:,:)      
      REAL*4,  ALLOCATABLE :: AD60_AQOX(:,:,:,:)      
      REAL*4,  ALLOCATABLE :: AD60_ERROR(:,:,:,:)      
      REAL*4,  ALLOCATABLE :: AD60_SOA(:,:,:,:)

      ! For ND61 -- 3D TOMAS rate diagnostic
      REAL*4,  ALLOCATABLE :: AD61(:,:,:,:)
      REAL*4,  ALLOCATABLE :: AD61_inst(:,:,:,:)
#endif

      ! For ND64 -- radiative flux from FAST-JX
      REAL*4,  ALLOCATABLE :: AD64(:,:,:,:,:)

      ! For ND65 -- production/loss rates
      REAL*4,  ALLOCATABLE :: AD65(:,:,:,:)

      ! For ND66 -- I-6 fields diagnostic
      REAL*4,  ALLOCATABLE :: AD66(:,:,:,:)      

      ! For ND67 -- DAO surface fields diagnostic
      REAL*4,  ALLOCATABLE :: AD67(:,:,:)      

      ! For ND68 -- BXHEIGHT, AD, AVGW diagnostic
      REAL*4,  ALLOCATABLE :: AD68(:,:,:,:)      

      ! For ND69 -- DXYP diagnostic
      REAL*4,  ALLOCATABLE :: AD69(:,:,:)      

      ! For ND71 -- Surface peak tracer conc. diagnostic
      REAL*4,  ALLOCATABLE :: AD71(:,:,:)      
      REAL*4,  ALLOCATABLE :: AD71_HR(:,:,:)
      REAL*4,  ALLOCATABLE :: AD71_DAY(:,:,:)
      INTEGER              :: AD71_COUNT
      INTEGER              :: AD71_HRCT
      INTEGER              :: AD71_LHR
      INTEGER              :: AD71_LDAY

      ! For ND72 -- Radiation output diagnostic
      REAL*4,  ALLOCATABLE :: AD72(:,:,:)

      ! For ND64 -- ISORROPIA pH, H+, and aerosol H2O diagnostic:
      REAL*4,  ALLOCATABLE :: AD73(:,:,:,:)
#endif

      ! For ND23 -- CH3CCl3 lifetime diagnostic
      REAL(fp),  ALLOCATABLE :: DIAGCHLORO(:,:,:,:)
!
! !REVISION HISTORY:
!  30 Nov 1999 - A. Fiore - Initial version
!  (1 ) DIAG_MOD is written in Fixed-Format F90.
!  (2 ) Call subroutine CLEANUP at the end of the MAIN program to deallocate
!        the memory before the run stops.  It is always good style to free
!        any memory we have dynamically allocated when we don't need it
!        anymoren
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
!  (15) Added ND06 (dust emission) and ND07 (carbon aerosol emission) 
!        diagnostic arrays (rjp, tdf, bmy, 4/5/04)
!  (16) Added AD13_SO2_sh diagnostic array for ND13 (bec, bmy, 5/20/04)
!  (17) Added AD07_HC diagnostic array for ND07 (rjp, bmy, 7/13/04)
!  (18) Moved AD65 & FAMPL to "diag65_mod.f" (bmy, 7/20/04)
!  (19) Added array AD13_SO4_bf (bmy, 11/17/04)!
!  (20) Added extra arrays for ND03 mercury diagnostics (eck, bmy, 12/7/04)
!  (21) Added extra ND21 array for crystalline sulfur tracers.  Also remove
!        ND03 and ND48 arrays; they are obsolete (bmy, 1/21/05)
!  (22) Removed AD41 and AFTTOT arrays; they're obsolete (bmy, 2/17/05)
!  (23) Added AD09, AD09_em arrays for HCN/CH3CN simulation (xyp, bmy, 6/27/05)
!  (24) Added AD30 array for land/water/ice output (bmy, 8/18/05)
!  (25) Added AD54 array for time spend in the troposphere (phs, 9/22/06)
!  (26) Added CTO3 counter. Convert ND43 counter arrays from 2D to 3D, for
!        the variable tropopause. (phs, 1/19/07)
!  (27) Added AD10 and AD10em arrays for ND10 H2-HD-sim diag (phs, 9/18/07)
!  (28) Added CTO3_24h to account for time in the troposphere for O3 in
!        ND47 (phs, 11/17/08)
!  (29) Added AD52 for Gamma HO2 diagnostic. (jaegle, ccc, 2/26/09)
!  (30) Updated to save out GLYX production of SOAG in ND07.
!       (tmf, 3/6/09)
!  (31) Add LTO3 for ND45 diag. (ccc, 7/20/09)
!  (32) Add AD19, AD58, AD60 for CH4 (kjw, 8/18/09)
!  (33) AD13_NH3_an is 3D now (phs, 10/22/09)  
!  (34) Add AD59_NUMB, AD59_SULF, AD59_SALT, AD59_ECOB, AD59_ECIL, AD59_OCOB, 
!        AD59_OCIL, and AD59_DUST for size-resolved emission  (win, 1/25/10)
!  (35) Add AD60_COND, AD60_COAG, AD60_NUCL, AD60_AQOX, AD60_SOA, and 
!        AD60_ERROR for TOMAS process rate diagnostics (win, 1/25/10)
!  (36) Add AD61 and AD61_INST for saving 3-D TOMAS rate (win, 1/25/10)   
!  (37) Add counter for aromatics SOA and add AD57 diagnostic for potential
!       temperature. (fp, 2/3/10)
!  26 Aug 2010 - R. Yantosca - Added ProTeX headers
!  02 Apr 2013 - M. Payer    - Removed *NO, *NO2, and *NO3 arrays for ND43
!                              diagnostic. These are no longer needed because
!                              NO, NO2, and NO3 are now tracers.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  08 Nov 2013 - M. Sulprizio- Removed CTO3 and LTO3. They are no longer used
!                              because O3 is now a tracer.
!  15 Aug 2014 - R. Yantosca - Remove diagnostics now tracked by HEMCO
!  21 Aug 2014 - R. Yantosca - AD01 is now only a 3-D array (for Pb emissions)
!  12 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  20 Jan 2015 - R. Yantosca - Added DIAGINIT_TRACER_CONC routine
!  15 Apr 2015 - R. Yantosca - Add TARGET spec to MASSFLEW, MASSFLNS, MASSFLUP
!  20 Jan 2016 - M. Sulprizio- Remove AD06, AD07, and AD08 arrays since these
!                              diagnostics are tracked by HEMCO
!  23 Jun 2016 - R. Yantosca - Remove references to APM code; it is no longer
!                              compatible with the FlexChem implementation
!  02 Aug 2016 - M. Sulprizio- Move AD65 array here from diag_pl_mod.F
!  03 Oct 2016 - R. Yantosca - Bug fix: scalars cannot be ALLOCATABLE
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
! !IROUTINE: cleanup_diag
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_DIAG
!
! !REVISION HISTORY:
!  13 Dec 2002 - R. Yantosca - Initial version
!  (1 ) Now also deallocate AD13_NH3_an, AD13_NH3_bb, AD13_NH3_bf arrays
!        for the ND13 diagnostic.  (bmy, 12/13/02)
!  (2 ) Now also deallocate AD13_NH3_na array for ND13 (rjp, bmy, 3/23/03)
!  (3 ) Removed P24H and L24H, these are now defined within "tagged_ox_mod.f".
!       Now also deallocate AD03 array for Kr85 prod/loss (jsw, bmy, 8/20/03)
!  (4 ) Now also deallocate AD06 and AD07* arrays (rjp, bdf, bmy, 4/5/04)
!  (5 ) Now also deallocate AD08 array (rjp, bec, bmy, 4/20/04)
!  (6 ) Now also deallocaes AD13_SO2_sh array (bec, bmy, 5/20/04)
!  (7 ) Now also deallocates AD07_HC array (rjp, bmy, 7/13/04)
!  (8 ) Now also deallocate AD13_SO4_bf array (bmy, 11/17/04)
!  (9 ) Now deallocate extra arrays for ND03 diagnostics (eck, bmy, 12/7/04)
!  (10) Now deallocates AD21_cr array.  Remove reference to arrays for ND03
!        and ND48 diagnostics, they're obsolete. (cas, sas, bmy, 1/21/05)
!  (11) Removed AD41 and AFTTOT arrays; they're obsolete (bmy, 2/17/05)
!  (12) Now also deallocate AD09 and AD09_em (bmy, 6/27/05)
!  (13) Now deallocate AD30 (bmy, 8/18/05)
!  (14) Now deallocate CTO3, AD10, AD10em arrays (phs, 9/18/07)
!  (15) Now deallocate TOMAS related arrays (win, bmy, 1/25/10)
!  15 Feb 2011 - R. Yantosca - Add modifications for APM microphysics
!  15 Aug 2014 - R. Yantosca - Remove diagnostics now tracked by HEMCO
!   3 Oct 2016 - R. Yantosca - Bug fix: Don't deallocate ND71 scalar variables
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_DIAG begins here!
      !=================================================================
#if defined( BPCH_DIAG ) || defined( MODEL_GEOS )
      IF ( ALLOCATED( AD01        ) ) DEALLOCATE( AD01        )
      IF ( ALLOCATED( AD02        ) ) DEALLOCATE( AD02        )
      IF ( ALLOCATED( AD07_BC     ) ) DEALLOCATE( AD07_BC     )
      IF ( ALLOCATED( AD07_OC     ) ) DEALLOCATE( AD07_OC     )
      IF ( ALLOCATED( AD07_HC     ) ) DEALLOCATE( AD07_HC     )
      IF ( ALLOCATED( AD09        ) ) DEALLOCATE( AD09        )
      IF ( ALLOCATED( AD09_em     ) ) DEALLOCATE( AD09_em     )
      IF ( ALLOCATED( AD12        ) ) DEALLOCATE( AD12        )
      IF ( ALLOCATED( AD16        ) ) DEALLOCATE( AD16        )
      IF ( ALLOCATED( AD17        ) ) DEALLOCATE( AD17        )
      IF ( ALLOCATED( AD18        ) ) DEALLOCATE( AD18        )
      IF ( ALLOCATED( AD21        ) ) DEALLOCATE( AD21        )
      IF ( ALLOCATED( AD22        ) ) DEALLOCATE( AD22        ) 
      IF ( ALLOCATED( AD30        ) ) DEALLOCATE( AD30        ) 
      IF ( ALLOCATED( AD31        ) ) DEALLOCATE( AD31        ) 
      IF ( ALLOCATED( AD33        ) ) DEALLOCATE( AD33        )
      IF ( ALLOCATED( AD35        ) ) DEALLOCATE( AD35        )
      IF ( ALLOCATED( AD37        ) ) DEALLOCATE( AD37        )
      IF ( ALLOCATED( AD38        ) ) DEALLOCATE( AD38        )  
      IF ( ALLOCATED( AD39        ) ) DEALLOCATE( AD39        )
      IF ( ALLOCATED( AD43        ) ) DEALLOCATE( AD43        )
      IF ( ALLOCATED( AD44        ) ) DEALLOCATE( AD44        )
!------------------------------------------------------------------------------
! This diagnostic requires some fixes. Comment out for now (mps,3/2/18)
!      IF ( ALLOCATED( AD44b       ) ) DEALLOCATE( AD44b       )
!      IF ( ALLOCATED( AD_RA       ) ) DEALLOCATE( AD_RA       )
!      IF ( ALLOCATED( AD_ZL       ) ) DEALLOCATE( AD_ZL       )
!------------------------------------------------------------------------------
      IF ( ALLOCATED( AD45        ) ) DEALLOCATE( AD45        )
      IF ( ALLOCATED( AD47        ) ) DEALLOCATE( AD47        )
      IF ( ALLOCATED( AD52        ) ) DEALLOCATE( AD52        )
      IF ( ALLOCATED( AD54        ) ) DEALLOCATE( AD54        )
      IF ( ALLOCATED( AD55        ) ) DEALLOCATE( AD55        )
      IF ( ALLOCATED( AD57        ) ) DEALLOCATE( AD57        )
      IF ( ALLOCATED( AD19        ) ) DEALLOCATE( AD19        )
      IF ( ALLOCATED( AD58        ) ) DEALLOCATE( AD58        )
#if   defined( TOMAS )
      IF ( ALLOCATED( AD59_NUMB   ) ) DEALLOCATE( AD59_NUMB   ) 
      IF ( ALLOCATED( AD59_SULF   ) ) DEALLOCATE( AD59_SULF   )
      IF ( ALLOCATED( AD59_SALT   ) ) DEALLOCATE( AD59_SALT   )
      IF ( ALLOCATED( AD59_ECIL   ) ) DEALLOCATE( AD59_ECIL   )
      IF ( ALLOCATED( AD59_ECOB   ) ) DEALLOCATE( AD59_ECOB   )
      IF ( ALLOCATED( AD59_OCIL   ) ) DEALLOCATE( AD59_OCIL   )
      IF ( ALLOCATED( AD59_OCOB   ) ) DEALLOCATE( AD59_OCOB   )
      IF ( ALLOCATED( AD59_DUST   ) ) DEALLOCATE( AD59_DUST   )
      IF ( ALLOCATED( AD60_COND   ) ) DEALLOCATE( AD60_COND   )
      IF ( ALLOCATED( AD60_COAG   ) ) DEALLOCATE( AD60_COAG   )
      IF ( ALLOCATED( AD60_NUCL   ) ) DEALLOCATE( AD60_NUCL   )
      IF ( ALLOCATED( AD60_AQOX   ) ) DEALLOCATE( AD60_AQOX   )
      IF ( ALLOCATED( AD60_ERROR  ) ) DEALLOCATE( AD60_ERROR  )
      IF ( ALLOCATED( AD60_SOA    ) ) DEALLOCATE( AD60_SOA    )
      IF ( ALLOCATED( AD61        ) ) DEALLOCATE( AD61        )
      IF ( ALLOCATED( AD61_inst   ) ) DEALLOCATE( AD61_inst   ) 
#else
      IF ( ALLOCATED( AD60        ) ) DEALLOCATE( AD60        )
#endif
      IF ( ALLOCATED( AD64        ) ) DEALLOCATE( AD64        )
      IF ( ALLOCATED( AD65        ) ) DEALLOCATE( AD65        )
      IF ( ALLOCATED( AD66        ) ) DEALLOCATE( AD66        )
      IF ( ALLOCATED( AD68        ) ) DEALLOCATE( AD68        )
      IF ( ALLOCATED( AD69        ) ) DEALLOCATE( AD69        )
      IF ( ALLOCATED( AD71        ) ) DEALLOCATE( AD71        )
      IF ( ALLOCATED( AD71_HR     ) ) DEALLOCATE( AD71_HR     )
      IF ( ALLOCATED( AD71_DAY    ) ) DEALLOCATE( AD71_DAY    )
      IF ( ALLOCATED( AD72        ) ) DEALLOCATE( AD72        )
      IF ( ALLOCATED( AD73        ) ) DEALLOCATE( AD73        )
      IF ( ALLOCATED( CONVFLUP    ) ) DEALLOCATE( CONVFLUP    )
      IF ( ALLOCATED( CT16        ) ) DEALLOCATE( CT16        )
      IF ( ALLOCATED( CT17        ) ) DEALLOCATE( CT17        )
      IF ( ALLOCATED( CT18        ) ) DEALLOCATE( CT18        )
      IF ( ALLOCATED( CTJV        ) ) DEALLOCATE( CTJV        )
      IF ( ALLOCATED( CTO3_24h    ) ) DEALLOCATE( CTO3_24h    )
      IF ( ALLOCATED( CTOH        ) ) DEALLOCATE( CTOH        )
      IF ( ALLOCATED( CTO1D       ) ) DEALLOCATE( CTO1D       )
      IF ( ALLOCATED( LTO1D       ) ) DEALLOCATE( LTO1D       )
      IF ( ALLOCATED( CTO3P       ) ) DEALLOCATE( CTO3P       )
      IF ( ALLOCATED( LTO3P       ) ) DEALLOCATE( LTO3P       )
      IF ( ALLOCATED( CTHO2       ) ) DEALLOCATE( CTHO2       )
      IF ( ALLOCATED( CTOTH       ) ) DEALLOCATE( CTOTH       )
      IF ( ALLOCATED( DIAGCHLORO  ) ) DEALLOCATE( DIAGCHLORO  )
      IF ( ALLOCATED( LTJV        ) ) DEALLOCATE( LTJV        )
      IF ( ALLOCATED( LTOH        ) ) DEALLOCATE( LTOH        )
      IF ( ALLOCATED( LTHO2       ) ) DEALLOCATE( LTHO2       )
      IF ( ALLOCATED( LTOTH       ) ) DEALLOCATE( LTOTH       )
      IF ( ALLOCATED( MASSFLEW    ) ) DEALLOCATE( MASSFLEW    )
      IF ( ALLOCATED( MASSFLNS    ) ) DEALLOCATE( MASSFLNS    )
      IF ( ALLOCATED( MASSFLUP    ) ) DEALLOCATE( MASSFLUP    )
      IF ( ALLOCATED( TURBFLUP    ) ) DEALLOCATE( TURBFLUP    )
      IF ( ALLOCATED( STT_TEMPO2  ) ) DEALLOCATE( STT_TEMPO2  )
#endif

      END SUBROUTINE CLEANUP_DIAG
!EOC
      END MODULE DIAG_MOD 


