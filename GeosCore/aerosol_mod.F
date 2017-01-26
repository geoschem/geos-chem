!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: aerosol_mod.F
!
! !DESCRIPTION: Module AEROSOL\_MOD contains variables and routines for
!  computing optical properties for aerosols which are needed for both the
!  FAST-J photolysis and ND21 optical depth diagnostics. (bmy, 7/20/04,
!  2/10/09)
!\\
!\\
! !INTERFACE: 
!
      MODULE AEROSOL_MOD
! 
! !USES:
!
      USE PRECISION_MOD

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: AEROSOL_CONC
      PUBLIC :: CLEANUP_AEROSOL
      PUBLIC :: INIT_AEROSOL
      PUBLIC :: RDAER
!
! !PUBLIC DATA MEMBERS:
!
      !========================================================================
      ! BCPI        : Hydrophilic black carbon aerosol   [kg/m3]
      ! BCPO        : Hydrophobic black carbon aerosol   [kg/m3]
      ! OCPI        : Hydrophilic organic carbon aerosol [kg/m3]
      ! OCPO        : Hydrophilic organic carbon aerosol [kg/m3]
      ! OCPISOA     : Hydrophilic OC + SOA aerosol       [kg/m3]
      ! SALA        : Accumulation mode seasalt aerosol  [kg/m3] 
      ! SALC        : Coarse mode seasalt aerosol        [kg/m3]
      ! SO4_NH4_NIT : Lumped SO4-NH4-NIT aerosol         [kg/m3]
      ! SO4         : Sulfate aerosol                    [kg/m3]
      ! NH4         : Ammonium aerosol                   [kg/m3]
      ! NIT         : Inorganic nitrate aerosol          [kg/m3]
      ! SOILDUST    : Mineral dust aerosol from soils    [kg/m3]
      ! SLA         : Stratospheric liquid aerosol       [kg/m3]
      ! SPA         : Stratospheric particulate aerosol  [kg/m3]
      ! TSOA        : Terpene SOA                        [kg/m3]
      ! ISOA        : Isoprene SOA                       [kg/m3]
      ! ASOA        : Aromatic + IVOC SOA                [kg/m3]
      ! OPOA        : Aerosol product of SVOC oxidation  [kg/m3]
      ! SOAG        : SOA product of GLYX                [kg/m3]
      ! SOAM        : SOA product of MYLY                [kg/m3]
      ! PM25        : Particulate matter < 2.5 um        [kg/m3]
      !========================================================================
      REAL(fp), ALLOCATABLE, PUBLIC :: BCPI(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: BCPO(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: OCPI(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: OCPO(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: OCPISOA(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: SALA(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: SALC(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: SO4_NH4_NIT(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: SO4(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: NH4(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: NIT(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: FRAC_SNA(:,:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: SOILDUST(:,:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: SLA(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: SPA(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: TSOA(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: ISOA(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: ASOA(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: OPOA(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: SOAG(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: SOAM(:,:,:)
      REAL(fp), ALLOCATABLE, PUBLIC :: PM25(:,:,:)
!
! !DEFINED PARAMETERS:
!
      REAL(fp), PARAMETER, PUBLIC :: OCFPOA  = 1.4e+0_fp ! OM/OC for POA
      REAL(fp), PARAMETER, PUBLIC :: OCFOPOA = 2.1e+0_fp ! OM/OC for OPOA, OCPI,
                                                         !  and OCPO

      ! For SOAG, assume the total aerosol mass/glyoxal mass = 1.d0 
      ! for now (tmf, 1/7/09)
      REAL(fp), PARAMETER, PUBLIC :: OCFG = 1.e+0_fp

      ! For SOAM, assume the total aerosol mass/methylglyoxal mass = 1.d0 
      ! for now (tmf, 1/7/09)
      REAL(fp), PARAMETER, PUBLIC :: OCFM = 1.e+0_fp

!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Pye, H.O.T., and J.H. Seinfeld, "A global perspective on aerosol from
!        low-volatility organic compounds", Atmos. Chem. & Phys., Vol 10, pp
!        4377-4401, 2010.
!
! !REVISION HISTORY:
!  (1 ) Added AEROSOL_RURALBOX routine (bmy, 9/28/04)
!  (2 ) Now convert ABSHUM from absolute humidity to relative humidity in
!         AEROSOL_RURALBOX, using the same algorithm as in "gasconc.f".
!         (bmy, 1/27/05)
!  (3 ) Now references "tropopause_mod.f" (bmy, 8/22/05)
!  (4 ) Now add contribution of SOA4 into Hydrophilic OC (dkh, bmy, 5/18/06)
!  (5 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (6 ) Add support for variable tropopause (bdf, phs, 9/14/06)
!  (7 ) Now set OCF=2.1 in AEROSOL_CONC for consistency w/ carbon_mod.f
!       (tmf, 2/10/09)
!  (8 ) Add WTAREA and WERADIUS for dicarbonyl SOA production.  
!       WTAREA is the same as TAREA, but excludes dry dust, BCPO and OCPO; 
!       use same units as TAREA.
!       WERADIUS is same as ERADIUS, but excludes dry dust, BCPO and OCPO;
!       use same units as ERADIUS. (tmf, 3/2/09)
!  (9 ) Add SOAG and SOAM species. (tmf, ccc, 3/2/09)
!  (10) Modify AOD output to wavelength specified in jv_spec_aod.dat 
!       (clh, 05/07/10)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  05 Mar 2013 - R. Yantosca - Now make INIT_AEROSOL a public routine
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!  05 Nov 2014 - M. Yannetti - PRECISION_MOD Changed REAL*8 to REAL(fp)
!  05 Mar 2015 - R. Yantosca - Remove DO_READ_DATA option.  We can no longer
!                              have binary punch input in HPC environments.
!  16 May 2016 - M. Sulprizio- Remove AEROSOL_RURALBOX. The implementation of
!                              FlexChem has rendered the routine obsolete.
!  17 Jun 2016 - R. Yantosca - Now use locally-defined species ID flags
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Mass of hydrophobic aerosol from Mian Chin
      REAL(fp), ALLOCATABLE, SAVE   :: DAERSL(:,:,:,:)

      ! Mass of hydrophilic aerosol from Mian Chin
      REAL(fp), ALLOCATABLE, SAVE   :: WAERSL(:,:,:,:)

      ! Add tracer ID flags as module variables (bmy, 6/16/16)
      INTEGER :: id_BCPI,  id_BCPO,  id_DST1,  id_DST2
      INTEGER :: id_DST3,  id_DST4,  id_NH4,   id_NIT  
      INTEGER :: id_OCPO,  id_OCPI,  id_SALA,  id_SALC 
      INTEGER :: id_SO4,   id_SOAG,  id_SOAM,  id_SO4s
      INTEGER :: id_NITs,  id_POA1,  id_POA2,  id_OPOA1
      INTEGER :: id_OPOA2, id_TSOA1, id_TSOA2, id_TSOA3 
      INTEGER :: id_TSOA0, id_ISOA1, id_ISOA2, id_ISOA3
      INTEGER :: id_ASOAN, id_ASOA1, id_ASOA2, id_ASOA3
      INTEGER :: id_DUST1

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
! !IROUTINE: aerosol_conc
!
! !DESCRIPTION: Subroutine AEROSOL\_CONC computes aerosol concentrations in
!  kg/m3 from the tracer mass in kg in the Species array.  These are needed to
!  compute optical properties for photolysis, for the optical depth diagnostics,
!  and for the SOA concentration diagnostics. (bmy, 7/20/04, 2/10/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AEROSOL_CONC( am_I_Root, Input_Opt, 
     &                         State_Met, State_Chm, RC )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_TROP
      USE CMN_FJX_MOD,        ONLY : REAA
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Met_Mod,      ONLY : MetState
      USE State_Chm_Mod,      ONLY : ChmState
      USE UCX_MOD,            ONLY : KG_STRAT_AER
      USE UnitConv_Mod
#if   defined( TOMAS ) 
      USE TOMAS_MOD,          ONLY : IBINS
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  This code was originally included in "chemdr.f", but the same computation 
!  also needs to be done for offline aerosol simulations.  Therefore, we have 
!  split this code off into a separate subroutine which can be called by both 
!  fullchem and offline aerosol simulations.
!
! !REVISION HISTORY: 
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Now add contribution from SOA4 into Hydrophilic OC (dkh, bmy, 5/18/06)
!  (3 ) Now set OCF=2.1 to be consistent w/ "carbon_mod.f" (tmf, 2/10/09)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  13 Nov 2012 - R. Yantosca - Now pass am_I_Root, Input_Opt, RC as arguments
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  05 Mar 2013 - R. Yantosca - Remove call to INIT_AEROSOL, this is now done
!                              in the initialization stage
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  17 Jun 2014 - J. Pierce   - Now allow bulk dust in TOMAS simulations
!  07 Apr 2015 - R. Yantosca - Add check to prevent div-by-zero errors
!  25 Jun 2015 - E. Lundgren - Use L. Zhang dust size distribution params for 
!                              mineral dust aerosols (non-TOMAS)
!  13 Aug 2015 - E. Lundgren - Convert tracer units to kg if input units are
!                              mixing ratio
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  10 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  18 Nov 2016 - M. Sulprizio- Move calculation of SOA and PM2.5 concentrations
!                              here from diag42_mod.F
!  22 Nov 2016 - M. Sulprizio- Add online calculation of aerosol growth factors
!                              at 35% RH
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd variables
      LOGICAL,  SAVE      :: FIRST = .TRUE.
      REAL(fp), SAVE      :: SIA_GROWTH
      REAL(fp), SAVE      :: ORG_GROWTH
      REAL(fp), SAVE      :: SSA_GROWTH

      ! Non-SAVEd variables
      INTEGER             :: I, J, L, N, K
      INTEGER             :: k_SO4
      INTEGER             :: k_ORG
      INTEGER             :: k_SSA
      REAL(fp)            :: Rad_wet, Rad_dry
      REAL(fp)            :: Rho_wet, Rho_dry

      ! Logical flags
      LOGICAL             :: prtDebug
      LOGICAL             :: UNITCONVERSION
      LOGICAL             :: LCARB
      LOGICAL             :: LDUST
      LOGICAL             :: LSOA
      LOGICAL             :: LSSALT
      LOGICAL             :: LSULF
      LOGICAL             :: LPRT
      LOGICAL             :: IS_OCPO,  IS_OCPI,  IS_BC
      LOGICAL             :: IS_SO4,   IS_NH4,   IS_NIT
      LOGICAL             :: IS_SAL,   IS_DST
      LOGICAL             :: IS_TSOA,  IS_ISOA,  IS_ASOA
      LOGICAL             :: IS_POA,   IS_OPOA
      LOGICAL             :: IS_SOAG,  IS_SOAM

      ! Pointers
      REAL(fp), POINTER   :: Spc(:,:,:,:)
      REAL(fp), POINTER   :: AIRVOL(:,:,:)

      !=================================================================
      ! AEROSOL_CONC begins here!
      !=================================================================

      ! Assume success
      RC     = GC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      LCARB  = Input_Opt%LCARB
      LDUST  = Input_Opt%LDUST
      LSOA   = Input_Opt%LSOA
      LSSALT = Input_Opt%LSSALT
      LSULF  = Input_Opt%LSULF
      LPRT   = Input_Opt%LPRT

      ! Do we have to print debug output?
      prtDebug   = ( LPRT .and. am_I_Root )

      ! Define logical flags
      IS_OCPI    = ( id_OCPI  > 0 )
      IS_OCPO    = ( id_OCPO  > 0 )
      IS_BC      = ( id_BCPI  > 0 .AND. id_BCPO  > 0 )
      IS_SO4     = ( id_SO4   > 0 )
      IS_NH4     = ( id_NH4   > 0 )
      IS_NIT     = ( id_NIT   > 0 )
      IS_DST     = ( id_DST1  > 0 .AND. id_DST2  > 0 )
      IS_SAL     = ( id_SALA  > 0 .AND. id_SALC  > 0 )
      IS_TSOA    = ( id_TSOA1 > 0 .AND. id_TSOA2 > 0 
     &         .AND. id_TSOA3 > 0 .AND. id_TSOA0 > 0 )
      IS_ISOA    = ( id_ISOA1 > 0 .AND. id_ISOA2 > 0 
     &         .AND. id_ISOA3 > 0)
      IS_ASOA    = ( id_ASOAN > 0 .AND. id_ASOA1 > 0 
     &         .AND. id_ASOA2 > 0 .AND. id_ASOA3 > 0 )
      IS_POA     = ( id_POA1  > 0 .AND. id_POA2  > 0 )
      IS_OPOA    = ( id_OPOA1 > 0 .AND. id_OPOA2 > 0 )
      IS_SOAG    = ( id_SOAG  > 0 )
      IS_SOAM    = ( id_SOAM  > 0 )

      ! If species units are 'kg/kg dry', convert to [kg] for this routine.
      ! Otherwise, if tracer units are not 'kg', throw an error to the 
      ! calling routine (ewl, 8/12/15)
      UNITCONVERSION = .FALSE.
      IF ( TRIM( State_Chm%Spc_Units ) .eq. 'kg/kg dry' ) THEN
         CALL ConvertSpc_KgKgDry_to_Kg( am_I_Root, State_Met, 
     &                                  State_Chm, RC         )  
         UNITCONVERSION = .TRUE.
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GC_Error('Unit conversion error', RC, 
     &                    'AEROSOL_CONC in aerosol_mod.F')
            RETURN
         ENDIF
      ELSE IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
         CALL GC_Error('Unexpected species units' // 
     &                 State_Chm%Spc_Units, RC, 
     &                 'AEROSOL_CONC in aerosol_mod')
         RETURN
      ENDIF

      ! Initialize pointers
      Spc    => State_Chm%Species
      AIRVOL => State_Met%AIRVOL

      !=================================================================
      ! Compute growth factors at 35% RH
      !
      ! GF = 1 + [ ( r_wet / r_dry )^3 -1 ] * [ rho_wet - rho_dry ]
      !
      ! and use rho_wet = 1000 kg/m3 
      !=================================================================
      IF ( FIRST ) THEN

         ! Species index of REAA from RD_AOD (in fast_jx_mod.F)
         k_SO4      = 1
         k_ORG      = 3
         k_SSA      = 4

         ! Density of H2O [kg/m3]
         Rho_wet    = 1000e+0_fp

         ! Growth factor for SO4 + NIT + NH4
         Rad_dry    = REAA(1,k_SO4)
         Rad_wet    = REAA(1,k_SO4) + 35e+0_fp *
     &              ( REAA(2,k_SO4) - REAA(1,k_SO4) ) / 50e+0_fp
         Rho_dry    = State_Chm%SpcData(id_SO4)%Info%Density
         SIA_GROWTH = 1 + ( ( ( Rad_wet / Rad_dry ) ** 3 - 1 ) *
     &                        ( Rho_wet / Rho_dry ) )

         ! Growth factor for OCPI + SOA
         Rad_dry    = REAA(1,k_ORG)
         Rad_wet    = REAA(1,k_ORG) + 35e+0_fp *
     &              ( REAA(2,k_ORG) - REAA(1,k_ORG) ) / 50e+0_fp
         IF ( IS_POA ) THEN
            Rho_dry    = State_Chm%SpcData(id_POA1)%Info%Density
         ELSE IF ( IS_OCPI ) THEN
            Rho_dry    = State_Chm%SpcData(id_OCPI)%Info%Density
         ENDIF
         ORG_GROWTH = 1 + ( ( ( Rad_wet / Rad_dry ) ** 3 - 1 ) *
     &                        ( Rho_wet / Rho_dry ) )

         ! Growth factor for SALA
         Rad_dry    = REAA(1,k_SSA)
         Rad_wet    = REAA(1,k_SSA) + 35e+0_fp *
     &              ( REAA(2,k_SSA) - REAA(1,k_SSA) ) / 50e+0_fp
         Rho_dry    = State_Chm%SpcData(id_SALA)%Info%Density
         SSA_GROWTH = 1 + ( ( ( Rad_wet / Rad_dry ) ** 3 - 1 ) *
     &                        ( Rho_wet / Rho_dry ) )

         ! Print debug info
         IF ( prtDebug ) THEN
            WRITE( 6,'(a)') 'Growth factors at 35% RH:'
            WRITE( 6, 100 ) SIA_GROWTH, ' for SO4, NIT, and NH4'
            WRITE( 6, 100 ) ORG_GROWTH, ' for OCPI and SOA'
            WRITE( 6, 100 ) SSA_GROWTH, ' for SALA'
 100        FORMAT(F5.2,A)
         ENDIF

      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, K )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      
         !==============================================================
         ! S U L F A T E   A E R O S O L S
         !
         ! Dump hydrophilic aerosols into one array that will be passed 
         ! to RDAER and then used for heterogeneous chemistry as well 
         ! as photolysis rate calculations interatively. 
         !
         ! For the full-chemistry run, If LSULF=F, then we read these 
         ! aerosol data from Mian's simulation.  If LSULF=T then we use 
         ! the online tracers.
         !
         ! Now assume that all sulfate, ammonium, and nitrate are 
         ! hydrophilic but sooner or later we can pass only hydrophilic 
         ! aerosols from the thermodynamic calculations for this 
         ! purpose.  This dumping should be done before calling INITGAS, 
         ! which converts the unit of Spc from kg/box to molec/cm3.
         !
         ! Units of SO4_NH4_NIT are [kg/m3].  (rjp, bmy, 3/23/03)
         !==============================================================
         IF ( LSULF ) THEN
#if defined( UCX )
            ! If we are using the full stratospheric chemistry mechanism,
            ! stratospheric NH4 is ignored, stratospheric NIT is taken
            ! as available for NAT formation and stratospheric SO4 is
            ! taken as sulfuric acid
            IF ( ITS_IN_THE_TROP( I, J, L, State_Met ) ) THEN

               ! Tropospheric - keep as normal
               ! now including sulfate and nitrate associated with sea-salt
               ! NOTE: these should be treated as having a sea-salt size
               ! distribution but are currently simply treated in the same
               ! way (size and optics) as all other sulfate aerosol (DAR
               ! 2013)
               SO4_NH4_NIT(I,J,L) = ( Spc(I,J,L,id_SO4)    + 
     &                                Spc(I,J,L,id_NH4)    +
     &                                Spc(I,J,L,id_NIT)    +
     &                                Spc(I,J,L,id_SO4s)   +
     &                                Spc(I,J,L,id_NITs) ) /
     &                                AIRVOL(I,J,L)
               SO4(I,J,L) = Spc(I,J,L,id_SO4) / AIRVOL(I,J,L)
               NH4(I,J,L) = Spc(I,J,L,id_NH4) / AIRVOL(I,J,L)
               NIT(I,J,L) = Spc(I,J,L,id_NIT) / AIRVOL(I,J,L)
               SLA(I,J,L) = 0e+0_fp
               SPA(I,J,L) = 0e+0_fp

            ELSE

               ! Tropospheric sulfate is zero in stratosphere
               SO4_NH4_NIT(I,J,L) = 0e+0_fp
               SO4(I,J,L) = 0e+0_fp
               NH4(I,J,L) = 0e+0_fp
               NIT(I,J,L) = 0e+0_fp
               SLA(I,J,L) = KG_STRAT_AER(I,J,L,1) / AIRVOL(I,J,L)
               SPA(I,J,L) = KG_STRAT_AER(I,J,L,2) / AIRVOL(I,J,L)

            ENDIF
#else
            ! Compute SO4 aerosol concentration [kg/m3]
            ! now including sulfate and nitrate associated with sea-salt
            ! NOTE: these should be treated as having a sea-salt size
            ! distribution but are currently simply treated in the same
            ! way (size and optics) as all other sulfate aerosol (DAR
            ! 2013)
            SO4_NH4_NIT(I,J,L) = ( Spc(I,J,L,id_SO4)    + 
     &                             Spc(I,J,L,id_NH4)    +
     &                             Spc(I,J,L,id_NIT)    +
     &                             Spc(I,J,L,id_SO4s)   +
     &                             Spc(I,J,L,id_NITs) ) /
     &                             AIRVOL(I,J,L)

            SO4(I,J,L) = Spc(I,J,L,id_SO4) / AIRVOL(I,J,L)
            NH4(I,J,L) = Spc(I,J,L,id_NH4) / AIRVOL(I,J,L)
            NIT(I,J,L) = Spc(I,J,L,id_NIT) / AIRVOL(I,J,L)

#endif
            
            ! Add error check for safe division (bmy, 4/7/15)
            IF ( SO4_NH4_NIT(I,J,L) > 0e+0_fp ) THEN

               ! Save these fractions for partitioning of optics
               ! until later when these may be treated independently
               FRAC_SNA(I,J,L,1) = ( ( Spc(I,J,L,id_SO4 ) 
     &                           +     Spc(I,J,L,id_SO4s) )
     &                           /   AIRVOL(I,J,L)          )
     &                           / SO4_NH4_NIT(I,J,L)

               FRAC_SNA(I,J,L,2) = ( ( Spc(I,J,L,id_NIT )
     &                           +     Spc(I,J,L,id_NITs) )
     &                           /   AIRVOL(I,J,L)          )
     &                           / SO4_NH4_NIT(I,J,L)

               FRAC_SNA(I,J,L,3) = ( Spc(I,J,L,id_NH4)
     &                           /   AIRVOL(I,J,L)        ) 
     &                           / SO4_NH4_NIT(I,J,L)
               
            ELSE
               
               ! If SO4_NH4_NIT(I,J,L) is zero, then avoid a div-by-zero
               ! error.  Set all of these to zero because the division
               ! cannot be done.
               FRAC_SNA(I,J,L,1) = 0e+0_fp
               FRAC_SNA(I,J,L,2) = 0e+0_fp
               FRAC_SNA(I,J,L,3) = 0e+0_fp

            ENDIF

         ENDIF

         !==============================================================
         ! C A R B O N  &  2 n d A R Y   O R G A N I C   A E R O S O L S
         !
         ! Compute hydrophilic and hydrophobic BC and OC in [kg/m3]
         ! Also add online 2ndary organics if necessary
         !==============================================================
         IF ( LCARB ) THEN

            ! Hydrophilic BC [kg/m3]
            BCPI(I,J,L) = Spc(I,J,L,id_BCPI) / AIRVOL(I,J,L)

            ! Hydrophobic BC [kg/m3]
            BCPO(I,J,L) = Spc(I,J,L,id_BCPO) / AIRVOL(I,J,L)

            ! Hydrophobic OC [kg/m3]
            ! SOAupdate: Treat either OCPO (x2.1) or POA (x1.4) 
            IF ( IS_POA ) THEN
               OCPO(I,J,L) = ( Spc(I,J,L,id_POA1) + Spc(I,J,L,id_POA2) )
     &                       * OCFPOA / AIRVOL(I,J,L)
            ELSE IF ( IS_OCPO ) THEN
               OCPO(I,J,L) = Spc(I,J,L,id_OCPO)
     &                       * OCFOPOA / AIRVOL(I,J,L)
            ENDIF

            ! Hydrophilic OC [kg/m3]
            IF ( IS_OCPI ) THEN
               OCPI(I,J,L) = Spc(I,J,L,id_OCPI)
     &                       * OCFOPOA / AIRVOL(I,J,L)
            ENDIF

            IF ( LSOA ) THEN

               ! Hydrophilic primary OC plus SOA [kg/m3].
               !
               ! We need to multiply by OCF to account for the mass of
               ! other components which are attached to the OC aerosol.
               ! (rjp, bmy, 7/15/04)
               !
               ! SOAupdate: Update traditional SOA (hotp 7/21/10)
               ! for new mtp + isop + lumparomivoc (hotp 5/20/10)
               IF ( IS_TSOA .and. IS_ISOA .and. IS_ASOA ) THEN
                  OCPISOA(I,J,L) = ( Spc(I,J,L,id_TSOA1)
     &                             + Spc(I,J,L,id_TSOA2)
     &                             + Spc(I,J,L,id_TSOA3)
     &                             + Spc(I,J,L,id_TSOA0)
     &                             + Spc(I,J,L,id_ISOA1)
     &                             + Spc(I,J,L,id_ISOA2)
     &                             + Spc(I,J,L,id_ISOA3) 
     &                             + Spc(I,J,L,id_ASOAN) 
     &                             + Spc(I,J,L,id_ASOA1)
     &                             + Spc(I,J,L,id_ASOA2)
     &                             + Spc(I,J,L,id_ASOA3) )
     &                             / AIRVOL(I,J,L) 
               ENDIF

               IF ( IS_OPOA ) THEN  ! hotp 7/28/10
                  OCPISOA(I,J,L) = OCPISOA(I,J,L) +
     &                             ( Spc(I,J,L,id_OPOA1)
     &                             + Spc(I,J,L,id_OPOA2) )
     &                             * OCFOPOA / AIRVOL(I,J,L)
               ENDIF

               IF ( IS_OCPI ) THEN  ! hotp 7/28/10
                  OCPISOA(I,J,L) = OCPISOA(I,J,L)
     &                             + Spc(I,J,L,id_OCPI)
     &                             * OCFOPOA / AIRVOL(I,J,L) 
               ENDIF

               ! Check to see if we are simulating SOAG and SOAM (tmf, 1/7/09)
               IF ( IS_SOAG ) THEN
                 OCPISOA(I,J,L) = OCPISOA(I,J,L)
     &                            + Spc(I,J,L,id_SOAG)
     &                            * OCFG / AIRVOL(I,J,L)
               ENDIF

               IF ( IS_SOAM ) THEN
                 OCPISOA(I,J,L) = OCPISOA(I,J,L)
     &                            + Spc(I,J,L,id_SOAM)
     &                            * OCFM / AIRVOL(I,J,L)
               ENDIF

            ELSE

               ! Hydrophilic primary and SOA OC [kg/m3].   We need
               ! to multiply by OCF to account for the mass of other 
               ! components which are attached to the OC aerosol.
               ! (rjp, bmy, 7/15/04)
               ! SOAupdate: use 2.1 (OCFOPOA) (hotp 7/21/10)
               OCPISOA(I,J,L) = Spc(I,J,L,id_OCPI)
     &                          * OCFOPOA / AIRVOL(I,J,L)
                  
            ENDIF ! LSOA

            ! Now avoid division by zero (bmy, 4/20/04)
            BCPI(I,J,L)    = MAX( BCPI(I,J,L), 1e-35_fp )
            OCPI(I,J,L)    = MAX( OCPI(I,J,L), 1e-35_fp )
            BCPO(I,J,L)    = MAX( BCPO(I,J,L), 1e-35_fp )
            OCPO(I,J,L)    = MAX( OCPO(I,J,L), 1e-35_fp )
            OCPISOA(I,J,L) = MAX( OCPISOA(I,J,L), 1e-35_fp )
         
         ENDIF ! LCARB
            
         !===========================================================
         ! M I N E R A L   D U S T   A E R O S O L S
         !
         ! NOTE: We can do better than this! Currently we carry 4 
         ! dust tracers...but het. chem and fast-j use 7 dust bins 
         ! hardwired from Ginoux.
         !
         ! Now, I apportion the first dust tracer into four smallest 
         ! dust bins equally in mass for het. chem and fast-j. 
         ! 
         ! Maybe we need to think about chaning our fast-j and het. 
         ! chem to use just four dust bins or more flexible 
         ! calculations depending on the number of dust bins. 
         ! (rjp, 03/27/04)
         !
         ! Now splitting mass into bins in fractions derived from
         ! Highwood et al. (2003).  Data is from log-normal fit of
         ! PCASP measurements of Saharan dust (Solid line in Fig.4b)
         ! (dar, 04/25/10) [see Ridley et al., 2012, JGR]
         !
         ! Updated for TOMAS (Jeffrey Pierce, 6/17/14)
         !===========================================================
#if defined( TOMAS )

         !-----------------------------------------------------------
         ! TOMAS simulations only 
         !-----------------------------------------------------------
         IF ( LDUST ) THEN

            ! Zero SOILDUST
            SOILDUST(I,J,L,:) = 0.e0_fp

            ! Loop over the # of TOMAS dust bins
            DO K = 1, IBINS

               ! Bin #1
               IF ( Input_Opt%DUSTREFF(K) < 0.2e-6_fp ) THEN
                  SOILDUST(I,J,L,1) = SOILDUST(I,J,L,1)
     &                              + Spc(I,J,L,id_DUST1+K-1)  
     &                              / AIRVOL(I,J,L)


               ! Bin #2
               ELSE IF ( Input_Opt%DUSTREFF(K) < 0.325e-6_fp ) THEN
                  SOILDUST(I,J,L,2) = SOILDUST(I,J,L,2)
     &                              + Spc(I,J,L,id_DUST1+K-1) 
     &                              / AIRVOL(I,J,L)

               ! Bin #3
               ELSE IF ( Input_Opt%DUSTREFF(K) < 0.6e-6_fp ) THEN
                  SOILDUST(I,J,L,3) = SOILDUST(I,J,L,3)
     &                              + Spc(I,J,L,id_DUST1+K-1) 
     &                              / AIRVOL(I,J,L)

               ! Bin #4
               ELSE IF ( Input_Opt%DUSTREFF(K) < 1.15e-6_fp ) THEN
                  SOILDUST(I,J,L,4) = SOILDUST(I,J,L,4)
     &                              + Spc(I,J,L,id_DUST1+K-1) 
     &                              / AIRVOL(I,J,L)

               ! Bin #5
               ELSE IF ( Input_Opt%DUSTREFF(K) < 2.0e-6_fp ) THEN
                  SOILDUST(I,J,L,5) = SOILDUST(I,J,L,5)
     &                              + Spc(I,J,L,id_DUST1+K-1) 
     &                              / AIRVOL(I,J,L)

               ! Bin #6
               ELSE IF ( Input_Opt%DUSTREFF(K) < 3.25e-6_fp ) THEN
                  SOILDUST(I,J,L,6) = SOILDUST(I,J,L,6)
     &                              + Spc(I,J,L,id_DUST1+K-1) 
     &                              / AIRVOL(I,J,L)

               ! Bin #7
               ELSE 
                  SOILDUST(I,J,L,7) = SOILDUST(I,J,L,7)
     &                              + Spc(I,J,L,id_DUST1+K-1) 
     &                              / AIRVOL(I,J,L)

               ENDIF
            ENDDO

         ENDIF
#else

         !-----------------------------------------------------------
         ! Preserve original code for non-TOMAS simulations
         !-----------------------------------------------------------
         IF ( LDUST ) THEN

            ! Lump 1st dust tracer for het chem
            ! Now use dust size distribution scheme to improve PM2.5
            ! surface dust conc over western U.S. (L. Zhang, 6/25/15)
            SOILDUST(I,J,L,1) = 
     &            0.007e+0_fp * Spc(I,J,L,id_DST1) / AIRVOL(I,J,L)
            SOILDUST(I,J,L,2) =
     &            0.0332e+0_fp * Spc(I,J,L,id_DST1) / AIRVOL(I,J,L)
            SOILDUST(I,J,L,3) =
     &            0.2487e+0_fp * Spc(I,J,L,id_DST1) / AIRVOL(I,J,L)
            SOILDUST(I,J,L,4) =
     &            0.7111e+0_fp * Spc(I,J,L,id_DST1) / AIRVOL(I,J,L)

            ! Other hetchem bins
            SOILDUST(I,J,L,5) = Spc(I,J,L,id_DST2) / AIRVOL(I,J,L)
            SOILDUST(I,J,L,6) = Spc(I,J,L,id_DST3) / AIRVOL(I,J,L)
            SOILDUST(I,J,L,7) = Spc(I,J,L,id_DST4) / AIRVOL(I,J,L)
            
         ENDIF

#endif

         !===========================================================
         ! S E A S A L T   A E R O S O L S
         !
         ! Compute accumulation & coarse mode concentration [kg/m3]
         !===========================================================
         IF ( LSSALT ) THEN

            ! Accumulation mode seasalt aerosol [kg/m3]
            SALA(I,J,L) = Spc(I,J,L,id_SALA) / AIRVOL(I,J,L)
            
            ! Coarse mode seasalt aerosol [kg/m3]
            SALC(I,J,L) = Spc(I,J,L,id_SALC) / AIRVOL(I,J,L)

            ! Avoid division by zero
            SALA(I,J,L) = MAX( SALA(I,J,L), 1e-35_fp )
            SALC(I,J,L) = MAX( SALC(I,J,L), 1e-35_fp )

         ENDIF

         !===========================================================
         ! S E C O N D A R Y   O R G A N I C   A E R O S O L S
         !
         ! Compute SOA concentration [kg/m3]
         !===========================================================
         IF ( LSOA ) THEN

            ! TSOA (terpene SOA) [kg/m3]
            IF ( IS_TSOA ) THEN
               TSOA(I,J,L) = ( Spc(I,J,L,id_TSOA1)
     &                       + Spc(I,J,L,id_TSOA2)
     &                       + Spc(I,J,L,id_TSOA3)
     &                       + Spc(I,J,L,id_TSOA0) ) / AIRVOL(I,J,L) 
            ENDIF

            ! ISOA (isoprene SOA) [kg/m3]
            IF ( IS_ISOA ) THEN
               ISOA(I,J,L) = ( Spc(I,J,L,id_ISOA1)
     &                       + Spc(I,J,L,id_ISOA2)
     &                       + Spc(I,J,L,id_ISOA3) ) / AIRVOL(I,J,L) 
            ENDIF

            ! ASOA (benz, tolu, xyle, + NAP/IVOC SOA) [kg/m3]
            IF ( IS_ASOA ) THEN
               ASOA(I,J,L) = ( Spc(I,J,L,id_ASOAN)
     &                       + Spc(I,J,L,id_ASOA1)
     &                       + Spc(I,J,L,id_ASOA2)
     &                       + Spc(I,J,L,id_ASOA3) ) / AIRVOL(I,J,L) 
            ENDIF

            ! OPOA [kg/m3]
            IF ( IS_OPOA ) THEN
               OPOA(I,J,L) = ( Spc(I,J,L,id_OPOA1)
     &                       + Spc(I,J,L,id_OPOA2) )
     &                       * OCFOPOA / AIRVOL(I,J,L) 
            ENDIF

         ENDIF

         !===========================================================
         ! SOAG and SOAM [kg/m3]
         !===========================================================
         IF ( IS_SOAG ) THEN
            SOAG(I,J,L) = Spc(I,J,L,id_SOAG) * OCFG / AIRVOL(I,J,L)
         ENDIF

         IF ( IS_SOAM ) THEN
            SOAM(I,J,L) = Spc(I,J,L,id_SOAM) * OCFM / AIRVOL(I,J,L)
         ENDIF

         !===========================================================
         ! P A R T I C U L A T E   M A T T E R
         !
         ! Compute PM2.5 concentration [kg/m3]
         !
         ! PM25 = 1.33 (NH4 + NIT  + SO4) + BCPI + BCPO +
         !        2.10 (OCPO + 1.16 OCPI) + 1.16 SOA    +
         !        DST1 + 0.38 DST2 + 1.86 SALA + SOAG + SOAM
         !
         ! NOTES:
         ! - We apply growth factors at 35% RH (computed above): 
         !    1.33 for SO4, NIT, and NH4
         !    1.16 for OCPI and SOA
         !    1.86 for SALA
         ! - Ratio of OM/OC = 2.1 is applied to OCPI and OCPO above
         ! - Aerosol WG recommends including 38% of DST2 in PM2.5
         !===========================================================

         PM25(I,J,L) = NH4(I,J,L)  * SIA_GROWTH
     &               + NIT(I,J,L)  * SIA_GROWTH
     &               + SO4(I,J,L)  * SIA_GROWTH
     &               + BCPI(I,J,L)
     &               + BCPO(I,J,L)
     &               + OCPO(I,J,L)
     &               + OCPI(I,J,L) * ORG_GROWTH
     &               + TSOA(I,J,L) * ORG_GROWTH
     &               + ISOA(I,J,L) * ORG_GROWTH
     &               + ASOA(I,J,L) * ORG_GROWTH
     &               + SALA(I,J,L) * SSA_GROWTH
     &               + SOILDUST(I,J,L,1) + SOILDUST(I,J,L,2) ! DST1
     &               + SOILDUST(I,J,L,3) + SOILDUST(I,J,L,4) ! DST1
     &               + SOILDUST(I,J,L,5) * 0.38              ! 38% of DST2
     &               + SOAG(I,J,L)
     &               + SOAM(I,J,L)

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! If units were originally 'kg/kg dry', convert [kg] back to mixing 
      ! ratio (ewl, 8/12/15)
      IF ( UNITCONVERSION ) THEN
         CALL ConvertSpc_Kg_to_KgKgDry( am_I_Root, State_Met, 
     &                                  State_Chm, RC         )  
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GC_Error('Unit conversion error', RC, 
     &                    'AEROSOL_CONC in aerosol_mod.F')
            RETURN
         ENDIF
      ENDIF

      ! Free pointers
      Spc    => NULL()
      AIRVOL => NULL()

      END SUBROUTINE AEROSOL_CONC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rdaer
!
! !DESCRIPTION: Subroutine RDAER reads global aerosol concentrations as
!  determined by Mian Chin.  Calculates optical depth at each level for
!  "set\_prof.f". Also calculates surface area for heterogeneous chemistry. It
!  uses aerosol parameters in FAST-J input file "jv\_spec.dat" for these
!  calculations. (rvm, rjp, tdf, bmy, 11/04/01, 7/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RDAER( am_I_Root, Input_Opt, State_Met,  State_Chm,
     &                  RC,        MONTH,     YEAR,       ODSWITCH   )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_CHEMGRID
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_NOCHEMGRID
      USE CMN_DIAG_MOD
      USE CMN_FJX_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD21
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants,      ONLY : CONSVAP
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      USE TIME_MOD,           ONLY : SYSTEM_TIMESTAMP
      USE TRANSFER_MOD,       ONLY : TRANSFER_3D
      USE UCX_MOD,            ONLY : GET_STRAT_OPT
      USE UCX_MOD,            ONLY : NDENS_AER

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root     ! Are we on root CPU?
      TYPE(MetState), INTENT(IN)  :: State_Met     ! Meteorology State
      INTEGER,        OPTIONAL    :: MONTH         ! # of current month
      INTEGER,        OPTIONAL    :: YEAR          ! 4-digit year
      INTEGER,        OPTIONAL    :: ODSWITCH      ! Logical indicator
                                                   !  = 0: AOD computed 
                                                   !       at 999 nm
                                                   !  = 1: AOD computed 
                                                   !       at wavelength set 
                                                   !       in Radiation Menu
      TYPE(OptInput), INTENT(IN)  :: Input_Opt     ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:                             
!
      INTEGER,        INTENT(OUT) :: RC            ! Success or failure?
!
! !REVISION HISTORY: 
!  (1 ) At the point in which "rdaer.f" is called, ABSHUM is actually
!        absolute humidity and not relative humidity (rvm, bmy, 2/28/02)
!  (2 ) Now force double-precision arithmetic by using the "D" exponent.
!        (bmy, 2/28/02)
!  (3 ) At present aerosol growth is capped at 90% RH.  The data
!        in jv_spec.dat could be used to allow a particle to grow to
!        99% RH if desired. (rvm, 3/15/02)
!  (4 ) Bug fix: TEMP2 needs to be sized (IIPAR,JJPAR,LLPAR) (bmy, 5/30/02)
!  (5 ) Now reference BXHEIGHT from "dao_mod.f".  Also references ERROR_STOP
!        from "error_mod.f".  Delete local declaration of TIME, since that
!        is also declared w/in comode.h -- this causes compile-time errors
!        on the ALPHA platform. (gcc, bmy, 11/6/02)
!  (6 ) Now use the online SO4, NH4, NIT aerosol, taken from the STT array, 
!        and passed via SO4_NH4_NIT argument if sulfate chemistry is turned on.
!        Otherwise, read monthly mean sulfate from disk.  (rjp, bmy, 3/23/03)
!  (7 ) Now call READ_BPCH2 with QUIET=.TRUE., which prevents info from being
!        printed to stdout.  Also made cosmetic changes. (bmy, 3/27/03)
!  (8 ) Add BCPI, BCPO, OCPI, OCPO to the arg list.  Bug fix: for online
!        sulfate & carbon aerosol tracers, now make sure these get updated
!        every timestep.  Now references "time_mod.f".  Now echo info about
!        which online/offline aerosols we are using.  Updated comments.
!        (bmy, 4/9/04)
!  (9 ) Add SALA, SALC to the arg list (rjp, bec, bmy, 4/20/04)
!  (10) Now references DATA_DIR from "directory_mod.f".  Now references LSULF,
!        LCARB, LSSALT from "logical_mod.f".  Added minor bug fix for 
!        conducting the appropriate scaling for optical depth for ND21
!        diagnostic.  Now make MONTH and YEAR optional arguments.  Now bundled
!        into "aerosol_mod.f".  (rvm, aad, clh, bmy, 7/20/04)
!  (11) Now remove FWET from extinction efficiency computation (avd, 8/3/10)
!  (12) Include third input argument to determine the wavelength at which
!        the AOD should be computed. This will set the optical properties
!        that are used for the calculation of the AOD. The ND21 diagnostic
!        should only be updated when WAVELENGTH = 1. (skim, 02/03/11)
!  09 Mar 2011 - R. Yantosca - Set MSDENS(2) = 1800 for APM (G. Luo)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  13 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC arguments for GIGC
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  27 Mar 2013 - S.D. Eastham- Upgraded from Fast-J to Fast-JX
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  05 Mar 2015 - R. Yantosca - Remove DO_READ_DATA option.  We can no longer
!                              have binary punch input in HPC environments.
!  31 Mar 2015 - R. Yantosca - Bug fix: Declare LINTERP as !$OMP PRIVATE
!  07 Apr 2015 - R. Yantosca - Bug fix: Add error check to JLOOP do-loop
!                              so that we cycle if JLOOP <= 0
!  24 Jun 2015 - M. Sulprizio- Update organic aerosol density MSDENS(3) from
!                              1800 to 1300 kg/m3 (E. Marais, M. Hammer)
!  15 Oct 2015 - C. Keller   - Empty ODAER before refilling.
!  12 May 2016 - M. Sulprizio- Remove 1D arrays that depend on JLOOP. ABSHUM
!                              is now a 3D field in State_Met. ERADIUS, TAREA
!                              WERADIUS, WTAREA are now pointers that point to
!                              3D fields in State_Chm.
!  16 May 2016 - M. Sulprizio- Remove JLOOP entirely and loop over LLPAR, JJPAR,
!                              IIPAR instead.
!  28 Jun 2016 - M. Sulprizio- Make IOUT a local variable to remove dependence
!                              on comode_loop_mod.F
!  22 Nov 2016 - M. Sulprizio- Now define aerosol densities in species database
!                              and use SpcData%Info%Density to populate MSDENS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL             :: FIRST = .TRUE.
      LOGICAL             :: LINTERP
      CHARACTER(LEN=16)   :: STAMP
      INTEGER             :: I, J, L, N, R, IRH, W, IRHN
      INTEGER             :: AA, IWV, IIWV, NWVS, IR, NRT
      INTEGER             :: IOUT
      REAL*4              :: TEMP(IIPAR,JJPAR,LGLOB)
      REAL(fp)            :: TEMP2(IIPAR,JJPAR,LLPAR)
      REAL(fp)            :: MSDENS(NAER), DRYAREA
      REAL(f8)            :: XTAU

      ! Variables for speed diagnostics
      INTEGER             :: ITIMEVALS(8)
      REAL*8              :: OLDSECS, NEWSECS
      REAL*8              :: OLDSECST, NEWSECST

      ! Effective radius at RH bins read in from "FJX_spec.dat"
      REAL(fp)		  :: RW(NRH)	

      ! Effective radius at RH after interpolation
      REAL(fp)		  :: REFF       

      ! Q at different RH bins read in from "FJX_spec.dat"
      REAL(fp)            :: AW(NRH)
      REAL(fp)            :: QW(NRH)	
      REAL(fp)            :: SSW(NRH)
      REAL(fp)            :: ASYW(NRH)	
      REAL(fp)            :: AW0
      REAL(fp)            :: QW0
      REAL(fp)            :: SSW0
      REAL(fp)            :: ASYW0
      REAL(fp)            :: DENOM,NUMER

      ! Used to interpolate between sizes
      REAL(fp)		  :: FRAC       
 
      ! Change in Q (extinction efficiency)
      REAL(fp)		  :: SCALEQ     

      ! Change in Radius with RH
      REAL(fp)		  :: SCALER     

      ! Change in SSA with RH
      REAL(fp)            :: SCALESSA

      ! Change in asym parameter with RH
      REAL(fp)            :: SCALEASY

      ! Chnge in Optical properties vs RH
      REAL(fp)            :: SCALEA
      REAL(fp)            :: SCALEOD 

      ! Change in Vol vs RH 
      REAL(fp)            :: SCALEVOL 

      ! Convert AbsHum to RelHum
      REAL(fp)  :: TK,CONSEXP,VPRESH2O,RELHUM

      ! Relative Humidities
      REAL(fp),  SAVE     :: RH(NRH)   = (/0e+0_fp,0.5e+0_fp, 
     &   0.7e+0_fp,0.8e+0_fp,0.9e+0_fp/)

      ! Temporary varaibles
      REAL(fp)            :: RAER, SADSTRAT, RHOSTRAT, XSASTRAT
      INTEGER             :: ISTRAT

      ! Index to aerosol types in FJX_spec.dat
      ! The following are ordered according to the mass densities below
#if defined( UCX )
      INTEGER, SAVE       :: IND(NAER) = (/22, 29, 36, 43, 50, 4, 14/)
#else
      INTEGER, SAVE	  :: IND(NAER) = (/22, 29, 36, 43, 50/)
#endif

      ! Local variables for quantities from Input_Opt
      LOGICAL             :: LCARB
      LOGICAL             :: LSSALT
      LOGICAL             :: LSULF
      LOGICAL             :: LSTRATOD
      LOGICAL             :: LRAD
      LOGICAL             :: IS_POA, IS_OCPI
      REAL(fp)            :: GF_RH

      ! Pointers
      REAL(fp), POINTER   :: BXHEIGHT(:,:,:)
      REAL(fp), POINTER   :: ERADIUS(:,:,:,:)
      REAL(fp), POINTER   :: TAREA(:,:,:,:)
      REAL(fp), POINTER   :: WERADIUS(:,:,:,:)
      REAL(fp), POINTER   :: WTAREA(:,:,:,:)

      !=================================================================
      ! RDAER begins here!
      !=================================================================
!     speed diagnostic
!            CALL DATE_AND_TIME( VALUES=ITIMEVALS )
!            OLDSECS=real(ITIMEVALS(5))*3600.0+real(ITIMEVALS(6))*60.0+
!     &              real(ITIMEVALS(7))+real(ITIMEVALS(8))/1000.0
!

      ! Assume success
      RC                   = GC_SUCCESS

      ! Copy fields from INPUT_OPT to local variables for use below
      LCARB                = Input_Opt%LCARB
      LSSALT               = Input_Opt%LSSALT
      LSULF                = Input_Opt%LSULF
      LSTRATOD             = Input_Opt%LSTRATOD
      LRAD                 = Input_Opt%LRAD

      ! Define logical flags
      IS_OCPI              = ( id_OCPI > 0 )
      IS_POA               = ( id_POA1 > 0 .AND. id_POA2 > 0 )

      ! Initialize pointers
      BXHEIGHT            => State_Met%BXHEIGHT    ! Grid box height [m]
      ERADIUS             => State_Chm%AeroRadi    ! Aerosol Radius [cm]
      TAREA               => State_Chm%AeroArea    ! Aerosol Area [cm2/cm3]
      WERADIUS            => State_Chm%WetAeroRadi ! Wet Aerosol Radius [cm]
      WTAREA              => State_Chm%WetAeroArea ! Wet Aerosol Area [cm2/cm3]

      !=================================================================
      ! S U L F A T E   A E R O S O L S
      !
      ! If LSULF = TRUE, then take the lumped SO4, NH4, NIT 
      ! concentrations [kg/m3] computed by AEROSOL_CONC, and save 
      ! into WAERSL(:,:,:,1) for use w/ FAST-J and hetchem.  This is 
      ! updated every timestep.  (For fullchem and offline runs)
      !
      ! If LSULF = FALSE, then read monthly mean offline sulfate aerosol   
      ! concentrations [kg/m3] from disk at the start of each month.
      ! (For fullchem simulations only)
      !=================================================================
      IF ( LSULF ) THEN 

         !-----------------------------------
         ! Use online aerosol concentrations
         !-----------------------------------
         IF ( FIRST ) THEN
            IF ( am_I_Root ) WRITE( 6, 100 )
 100        FORMAT( '     - RDAER: Using online SO4 NH4 NIT!' )
         ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            WAERSL(I,J,L,1) = SO4_NH4_NIT(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !=================================================================
      ! C A R B O N  &  2 n d A R Y   O R G A N I C   A E R O S O L S
      !
      ! If LCARB = TRUE, then take Hydrophilic OC, Hydrophobic OC,
      ! Hydropilic BC, and Hydrophobic BC, and 2ndary organic aerosol
      ! concentrations [kg/m3] that have been computed by AEROSOL_CONC.   
      ! Save these into DAERSL and WAERSL for use w/ FAST-J and hetchem.  
      ! These fields are updated every chemistry timestep.
      ! (For both fullchem and offline simulations)
      !
      ! If LCARB = FALSE, then read monthly mean carbon aerosol
      ! concentrations [kg/m3] from disk at the start of each month.
      ! (For full chemistry simulations only)
      !=================================================================
      IF ( LCARB ) THEN

         !-----------------------------------
         ! Use online aerosol concentrations
         !-----------------------------------
         IF ( FIRST ) THEN
            IF ( am_I_Root ) WRITE( 6, 110 )
 110        FORMAT( '     - RDAER: Using online BCPI OCPI BCPO OCPO!' )
         ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Hydrophilic BC (a.k.a EC) [kg/m3]
            WAERSL(I,J,L,2) = BCPI(I,J,L)

            ! Hydrophilic OC [kg/m3]
            WAERSL(I,J,L,3) = OCPISOA(I,J,L)

            ! Hydrophobic BC (a.k.a EC) [kg/m3]
            DAERSL(I,J,L,1) = BCPO(I,J,L)

            ! Hydrophobic OC [kg/m3]
            DAERSL(I,J,L,2) = OCPO(I,J,L)

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !=================================================================
      ! S E A S A L T   A E R O S O L S
      !
      ! If LSSALT = TRUE, then take accumulation and coarse mode
      ! seasalt aerosol concentrations [kg/m3] that are passed from
      ! "chemdr.f".  Save these into WAERSL for use w/ FAST-J and
      ! hetchem.  These fields are updated every chemistry timestep.
      ! (For both fullchem and offline simulations)
      !
      ! If LSSALT = FALSE, then read monthly-mean coarse sea-salt 
      ! aerosol concentrations [kg/m3] from the binary punch file.  
      ! Also merge the coarse sea salt aerosols into a combined bin 
      ! rather than carrying them separately.
      ! (For fullchem simulations only)
      !=================================================================
      IF ( LSSALT ) THEN

         !-----------------------------------
         ! Use online aerosol concentrations
         !-----------------------------------
         IF ( FIRST ) THEN
            IF ( am_I_Root ) WRITE( 6, 120 )
 120        FORMAT( '     - RDAER: Using online SALA SALC' )
         ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Accumulation mode seasalt aerosol [kg/m3]
            WAERSL(I,J,L,4) = SALA(I,J,L)

            ! Coarse mode seasalt aerosol [kg/m3]
            WAERSL(I,J,L,5) = SALC(I,J,L)

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF 

#if defined( UCX )
      ! SDE 04/17/13 - Transfer stratospheric aerosol data
      WAERSL(:,:,:,NRHAER+1) = SLA
      WAERSL(:,:,:,NRHAER+2) = SPA
#endif

      !=================================================================
      ! Calculate optical depth and surface area at each timestep
      ! to account for the change in relative humidity
      !
      ! For the optical depth calculation, this involves carrying the 
      ! optical depth at each RH as separate aerosols since OPMIE.f 
      ! treats the phase functions and single scattering albedos 
      ! separately. (An alternative would be to rewrite OPMIE.f)
      !
      ! Scaling is sufficient for the surface area calculation
      !=================================================================
      MSDENS(1) = State_Chm%SpcData(id_SO4)%Info%Density
      MSDENS(2) = State_Chm%SpcData(id_BCPI)%Info%Density
      IF ( IS_POA ) THEN
         MSDENS(3) = State_Chm%SpcData(id_POA1)%Info%Density
      ELSE IF ( IS_OCPI ) THEN
         MSDENS(3) = State_Chm%SpcData(id_OCPI)%Info%Density
      ENDIF
      MSDENS(4) = State_Chm%SpcData(id_SALA)%Info%Density
      MSDENS(5) = State_Chm%SpcData(id_SALC)%Info%Density
#if defined( UCX )
      ! These default values unused (actively retrieved from ucx_mod)
      MSDENS(NRHAER+1) = 1700.0d0! SSA/STS
      MSDENS(NRHAER+2) = 1000.0d0! NAT/ice PSC
#endif

      !set default values for RH index array
      IRHARR(:,:,:)=1

      ! Empty ODAER before refilling. This is required to make sure that
      ! all gridboxes outside the chemistry grid are zero and do not carry
      ! over values from previous time steps. (ckeller, 10/15/15)
      ODAER(:,:,:,:,:) = 0.0d0

      ! DAR 09/2013
      ! There are two ways RDAER can be called:
      ! (1) When Fast-J requires aerosol optics at 1000nm (ODSWITCH=0)
      ! (2) Before diags are accumulated, from RECOMPUTE_AOD (ODSWITCH=1) 
      ! for (1) we just need the optics stored at a single wavelength
      !     not for all user specified wavelengths, hence NWVS=1, IWV=IWV1000
      ! for (2) we need to determine if RRTMG is switched on
      !     if LRAD=true, calculation is over total minus standard wavelengths 
      !     in optics dat files (NWVAA-NWVAA0)
      !     if LRAD=false, calculation is for the wavelengths required for
      !     user-requested wavelength output. These are determined in CALC_AOD
      !     within RD_AOD.F and stored in IWVREQUIRED. The coefficients to 
      !     interpolate from the LUT wavelengths to the user-requested
      !     waveelenths (in CALC_AOD) are used here.

      ! Select number of wavelengths required to loop over
      IF (ODSWITCH .EQ. 0) THEN !this is the call for Fast_JX at 1000nm
         NWVS   = 1
      ELSE
         IF ( LRAD ) THEN       !Loop over all RT wavelengths (30)
            ! plus any required for calculating the AOD
            NWVS = NWVAA-NWVAA0+NWVREQUIRED
         ELSE                   !Loop over wavelengths needed for 
                                !interpolation to those requested in input.geos 
                                !(determined in RD_AOD)
            NWVS = NWVREQUIRED
         ENDIF                 
      ENDIF

      DO IIWV = 1, NWVS
         !now select the correct LUT wavelength
         IF (ODSWITCH .EQ. 0) THEN
            ! only doing for 1000nm (IWV1000 is set in RD_AOD)
            ! N.B. NWVS is fixed to 1 above - only one wavelength
            IWV=IWV1000
         ELSE
            IF ( LRAD ) THEN
               ! RRTMG wavelengths begin after NWVAA0 standard wavelengths
               ! but add on any others required 
               IF (IIWV.LE.30) THEN
                  !index of RRTMG wavelengths starts after the standard NWVAA0
                  !(currently NWVAA0=11, set in CMN_FJX_mod based on the .dat LUT)
                  IWV = IIWV+NWVAA0
               ELSE
                  !now we calculate at wvs for the requested AOD
                  IWV = IWVREQUIRED(IIWV-30)
               ENDIF
            ELSE
               ! IWVREQUIRED lists the index of requires standard wavelengths
               IWV = IWVREQUIRED(IIWV)
            ENDIF               
         ENDIF

      ! Loop over types of aerosol with hygroscopic growth
      ! (this will include strat aerosol if UCX=yes)
      DO N = 1, NAER 
         !index for strat aerosol (only >0 for strat aero)
         ISTRAT=N-NRHAER

         !NRT is subscript for RT arrays that contain SNA separately
         !so their optics can be treated separately in future
         IF (N.GT.1) THEN
            NRT=N+2
         ELSE
            NRT=N
         ENDIF

         !==============================================================
         ! Determine aerosol growth rates from the relative 
         ! humidity in each box
         !
         ! The optical depth scales with the radius and Q alone
         ! since SCALEDENS cancels as follows
         ! 
         !    SCALER 	= RW / RDRY
         !    SCALEDENS = DENSWET / DENSDRY
         !    SCALEM 	= SCALEDENS * SCALER**3
         !    SCALEOD 	= (SCALEQ * SCALEM) / (SCALEDENS * SCALER)
         !          	= SCALEQ * SCALER**2
         !
         ! Cap aerosol values at 90% relative humidity since
         ! aerosol growth at that point becomes highly nonlinear and 
         ! relative humidities above this value essentially mean
         ! there is a cloud in that grid box
         !
         ! Q is the extinction efficiency
         !
         ! Each grid box (I,J,L) will fall into one of the RH bins, 
         ! since each grid box will have a different RH value.  So,
         ! for SCALEOD(I,J,L,:), only one of the IRH bins will contain
         ! nonzero data, while the other IRH bins will all be zero.
         !==============================================================
         ! We loop over all selected wavelengths, referenced by IWV
         ! if FAST_J is calling then IWV will be for 1000nm only
         ! if RRTMG is on then IWV will be 30 wavelengths + AOD wavs,
         ! otherwise IWV will be at user input specified wavelengths

         ! Loop over relative humidity bins
         DO R = 1, NRH

            ! Wet radius in aerosol LUT files
            RW(R) = REAA(R,N)

            ! Extinction efficiency for Q for each RH bin
            QW(R)   = QQAA(IWV,R,N)
            AW(R)   = ALPHAA(IWV,R,N)
            SSW(R)  = SSAA(IWV,R,N)
            ASYW(R) = ASYMAA(IWV,R,N)

         ENDDO

         ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+PRIVATE( I,        J,       L,        IRH                 )
!$OMP+PRIVATE( AW0,      QW0,     SSW0,     ASYW0,    REFF      )
!$OMP+PRIVATE( SCALEA,   SCALEQ,  SCALESSA, SCALEASY, FRAC      )
!$OMP+PRIVATE( SCALER,   SCALEOD, SCALEVOL, DRYAREA             )
!$OMP+PRIVATE( TK,       CONSEXP, VPRESH2O, RELHUM              )
#if defined( RRTMG )
!$OMP+PRIVATE( IR                                               )
#endif
#if defined( UCX )
!$OMP+PRIVATE( RHOSTRAT, RAER,    SADSTRAT, XSASTRAT            )
#endif
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Skip non-chemistry boxes
            IF ( ITS_IN_THE_NOCHEMGRID( I, J, L, State_Met ) ) CYCLE

            ! Calculate RH. Not clear why the result of this calc is 
            ! slightly different than State_Met%RH
            RELHUM   = State_Met%AVGW(I,J,L) *
     &                 State_Met%AIRNUMDEN(I,J,L)
            TK       = State_Met%T(I,J,L)
            CONSEXP  = 17.2693882e+0_fp * (TK - 273.16e+0_fp) /
     &                 (TK - 35.86e+0_fp)
            VPRESH2O = CONSVAP * EXP(CONSEXP) / TK 
            RELHUM   = RELHUM / VPRESH2O 

            ! Sort into relative humidity bins
            ! Currently uses the 0, 50, 70, 80, 90% RH bins
            ! (95 and 99% also stored in data files but not used)
            IF (      RELHUM <= RH(2) ) THEN   
               IRH = 1
            ELSE IF ( RELHUM <= RH(3) ) THEN
               IRH = 2
            ELSE IF ( RELHUM <= RH(4) ) THEN
               IRH = 3
            ELSE IF ( RELHUM <= RH(5) ) THEN
               IRH = 4
            ELSE 
               IRH = 5
            ENDIF

            ! save the index of the relative humidity into array that is
            ! used by the photolysis to pull the right optics
            ! information from the FJX_spec.dat
            ! Previously the ODAER was organized by aerosol species and
            ! each RH bin, but this leaves a lot of the array redundant
            ! because only one RH bin is used at once
            ! This becomes a waste of memory with multiple wavelengths
            IRHARR(I,J,L) = IRH

            ! For the NRHth bin, we don't have to interpolate
            ! For the other bins, we have to interpolate 
            ! Sometimes values can be zero, don't want to divide by 0...
            IF (AW(1).EQ.0.0) THEN
               AW0=1.0
            ELSE
               AW0=AW(1)
            ENDIF
            IF (QW(1).EQ.0.0) THEN
               QW0=1.0
            ELSE
               QW0=QW(1)
            ENDIF
            IF (SSW(1).EQ.0.0) THEN
               SSW0=1.0
            ELSE
               SSW0=SSW(1)
            ENDIF
            IF (ASYW(1).EQ.0.0) THEN
               ASYW0=1.0
            ELSE
               ASYW0=ASYW(1)
            ENDIF

            IF ( IRH == NRH ) THEN
               REFF     = RW(NRH)
               SCALEA   = AW(NRH)   / AW0  !QW(1) is dry extinction eff.
               SCALEQ   = QW(NRH)   / QW0 
               SCALESSA = SSW(NRH)  / SSW0
               SCALEASY = ASYW(NRH) / ASYW0
            ELSE                
               ! Interpolate between different RH
               FRAC = (RELHUM-RH(IRH)) / (RH(IRH+1)-RH(IRH))
               IF ( FRAC > 1.0d0 ) FRAC = 1.0d0
               REFF    = FRAC*RW(IRH+1)  + (1.d0-FRAC)*RW(IRH)
               SCALEA  =(FRAC*AW(IRH+1)  + (1.d0-FRAC)*AW(IRH))/AW0
               SCALEQ  =(FRAC*QW(IRH+1)  + (1.d0-FRAC)*QW(IRH))/QW0
               SCALESSA=(FRAC*SSW(IRH+1) + (1.d0-FRAC)*SSW(IRH))/SSW0
               SCALEASY=(FRAC*ASYW(IRH+1)+ (1.d0-FRAC)*ASYW(IRH))/ASYW0
            ENDIF

            SCALER  = REFF / RW(1)
            SCALEOD = SCALEQ * SCALER * SCALER

#if defined( UCX )
            IF (N.LE.NRHAER) THEN
#endif
            !calculate optics for hyrdophillic aerosol here
            !However MDENS in LUT was in g/cm3 not kg/m3 so x1e3
            ODAER(I,J,L,IWV,N) = SCALEOD*BXHEIGHT(I,J,L)*0.75d0 
     &                         * WAERSL(I,J,L,N) * QQAA(IWV,1,N) / 
     &                         ( MSDENS(N) * REAA(1,N) * 1.0D-6 )

            IF ((N.eq.2).or.(N.eq.3)) THEN
               !now combine with hydrophilic OD as before
               ODAER(I,J,L,IWV,N)= ODAER(I,J,L,IWV,N) +
     &                             0.75d0 * BXHEIGHT(I,J,L) *
     &                             DAERSL(I,J,L,N-1) * QQAA(IWV,1,N)  /
     &                             ( MSDENS(N) * REAA(1,N) * 1.0D-6 )
            ENDIF

#if defined( UCX )
            ELSE           
            ! Get aerosol effective radius
            CALL GET_STRAT_OPT(I,J,L,ISTRAT,RAER,REFF,SADSTRAT,XSASTRAT)

            ! SDE 2014-02-04
            ! The calculation used for the aerosols above
            ! is essentially a roundabout way of deriving
            ! the cross-sectional area. For log-normally
            ! distributed aerosols, this is much easier,
            ! and a direct query prevents the possibility
            ! of dividing a small mass by a small calculated
            ! radius and blowing up

            ! Aerosol optical depth
            ODAER(I,J,L,IWV,N) = BXHEIGHT(I,J,L) * XSASTRAT
     &                           * QQAA(IWV,1,N)
            ENDIF

#endif

#if defined( RRTMG )            
            !SNA currently treated as one with optics but considered
            !separately for RT, so we split them by mass here
            IF (N.EQ.1) THEN
             DO IR=1,3
              RTODAER(I,J,L,IWV,N+IR-1)= ODAER(I,J,L,IWV,N)*
     &                                      FRAC_SNA(I,J,L,IR)
              RTSSAER(I,J,L,IWV,N+IR-1)   = SCALESSA*SSAA(IWV,1,N)
              RTASYMAER(I,J,L,IWV,N+IR-1) = SCALEASY*ASYMAA(IWV,1,N)
             ENDDO
            ELSE
             !RT arrays now offset from NAER by 2 (NRT=N+2 for N>1) 
             !If strat aerosol switched on (UCX) then this will
             !automatically be added after the standard aerosol
             !(NRHAER+1,2) but before dust
             RTODAER(I,J,L,IWV,NRT)     = ODAER(I,J,L,IWV,N)
             RTSSAER(I,J,L,IWV,NRT)     = SCALESSA*SSAA(IWV,1,N)
             RTASYMAER(I,J,L,IWV,NRT)   = SCALEASY*ASYMAA(IWV,1,N)
            ENDIF
#endif

            !==============================================================
            ! ND21 Diagnostic: 
            !
            ! Computed here:
            ! --------------
            ! #7  Hygroscopic growth of SO4                [unitless]
            ! #10 Hygroscopic growth of Black Carbon       [unitless]
            ! #13 Hygroscopic growth of Organic Carbon     [unitless]
            ! #16 Hygroscopic growth of Sea Salt (accum)   [unitless]
            ! #19 Hygroscopic growth of Sea Salt (coarse)  [unitless]
            !==============================================================
            !only need to do hyg once, not for each wavelength             
            IF ((IIWV.EQ.1).and.(ISTRAT.le.0)) THEN

               !save increase in AOD resulting from water uptake
               IF( ND21 > 0 .AND. L <= LD21 .AND.ODSWITCH.EQ.1) THEN
                  AD21(I,J,L,7+3*(N-1)) = AD21(I,J,L,7+3*(N-1)) +
     &                                    SCALEOD
               ENDIF
            
               !now calulate the surface areas
               !==============================================================
               !  Calculate Aerosol Surface Area
               !
               !  Units ==> AERSL    [ kg aerosol m^-3 air ]
               !            MSDENS   [ kg aerosol m^-3 aerosol ]
               !            ERADIUS  [ cm      ]
               !            TAREA    [ cm^2 dry aerosol/cm^3 air ]
               !
               !  Note: first find volume of aerosol (cm^3 arsl/cm^3 air), then
               !        multiply by 3/radius to convert to surface area in cm^2
               !
               !  Wet Volume = AERSL * SCALER**3 / MSDENS
               !  Wet Surface Area = 3 * (Wet Volume) / ERADIUS 
               !
               !  Effective radius for surface area and optical depths 
               !  are identical.
               !==============================================================
               !========================================================
               ! NOTES:          
               !    WAERSL   [ kg dry mass of wet aerosol m^-3 air ]
               !    ERADIUS  [ cm wet aerosol radius ]
               !    MSDENS   [ kg dry mass of aerosol m^-3 dry volume of aerosol ]
               !    TAREA    [ cm^2 wet sfc area of aerosol cm^-3 air ]   
               !    WTAREA   : same as TAREA, but excludes dry dust, BCPO, OCPO 
               !               use same units as TAREA    (tmf, 4/18/07) 
               !    WERADIUS : same as ERADIUS, but excludes dry dust, BCPO,OCPO
               !               use same units as ERADIUS  (tmf, 4/18/07)
               ! Wet dust WTAREA and WERADIUS are archived in dust_mod.f.
               !========================================================
            
               !get scaling for R and VOL
               SCALER                 = REFF / RW(1)
               SCALEVOL               = SCALER**3
               ERADIUS(I,J,L,N+NDUST) = 1.0D-4 * REFF

               ! Store aerosol surface areas in TAREA, and be sure
               ! to list them following the dust surface areas
               TAREA(I,J,L,N+NDUST)   = 3.D0                     * 
     &                                  WAERSL(I,J,L,N)          *  
     &                                  SCALEVOL                 / 
     &                                  ( ERADIUS(I,J,L,N+NDUST) * 
     &                                    MSDENS(N) )  

               WTAREA(I,J,L,N+NDUST)   = TAREA(I,J,L,N+NDUST)
               WERADIUS(I,J,L,N+NDUST) = ERADIUS(I,J,L,N+NDUST)     

               !include hydrophobic BC and OC
               !stored separate to hydrophillic in RT variables
               IF ((N.eq.2).or.(N.eq.3)) THEN
             
                  ! Dry surface area
                  ! SDE 2015-10-27: RW is in um, but everything
                  ! else is in terms of cm. Correct with 10^-4 factor
                  DRYAREA = 3.D0 * DAERSL(I,J,L,N-1) / ( RW(1) *
     &                      1.0D-4 * MSDENS(N) )  

                  ! Add surface area to TAREA array
                  TAREA(I,J,L,N+NDUST) = WTAREA(I,J,L,N+NDUST) + DRYAREA

                  ! Define a new effective radius that accounts 
                  ! for the hydrophobic aerosol 
                  ERADIUS(I,J,L,NDUST+N) = ( WERADIUS(I,J,L,NDUST+N) *
     &                                       WTAREA(I,J,L,N+NDUST)   +
     &                                       RW(1) * 1.0D-4 * DRYAREA )/
     &                                     ( WTAREA(I,J,L,N+NDUST)   +
     &                                       DRYAREA )
               ENDIF !Hydrophobic aerosol surface area
            ENDIF    !Surface area calcs
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDDO !Loop over NAER
      ENDDO !Loop over NWVSELECT

#if defined( UCX )
      !==============================================================
      ! Account for stratospheric aerosols (SDE 04/17/13)
      !==============================================================
      DO ISTRAT = 1,NSTRATAER

         ! Index for combination of aerosol type and RH
         N = NRHAER + ISTRAT
         ! For ND21
         IOUT = 5 + 2*NRHAER + (3*(NDUST+NRHAER)) + (ISTRAT-1)*3

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, RAER, REFF )
!$OMP+PRIVATE( SADSTRAT, XSASTRAT )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Get aerosol effective radius
            CALL GET_STRAT_OPT(I,J,L,ISTRAT,RAER,REFF,SADSTRAT,XSASTRAT)

            ! Moved this from a separate loop for clarity
            IF ( ITS_IN_THE_CHEMGRID( I, J, L, State_Met ) ) THEN

               ! Add surface area to TAREA array
               TAREA(I,J,L,NDUST+N)  = SADSTRAT
               WTAREA(I,J,L,NDUST+N) = SADSTRAT

               ! Store radius
               ERADIUS(I,J,L,NDUST+N)  = RAER
               WERADIUS(I,J,L,NDUST+N) = RAER 

            ENDIF

            ! Store SAD
            IF ((ND21.gt.0).and.(L.le.LD21)) THEN
                AD21(I,J,L,IOUT+2) = AD21(I,J,L,IOUT+2) + SADSTRAT
            ENDIF

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDDO
#endif
       
      !==============================================================
      ! ND21 Diagnostic: Aerosol OD's, Growth Rates, Surface Areas
      !
      ! Computed in other routines:
      ! ---------------------------------
      ! #1: Cloud optical depths (1000 nm) --> from "optdepth_mod.f"       
      ! #2: Max Overlap Cld Frac           --> from "optdepth_mod.f" 
      ! #3: Random Overlap Cld Frac        --> from "optdepth_mod.f" 
      ! #4: Dust optical depths            --> from "rdust.f"
      ! #5: Dust surface areas             --> from "rdust.f"
      !
      ! Computed previously in "rdaer.f":
      ! ---------------------------------
      ! #7  Hygroscopic growth of SO4                [unitless]
      ! #10 Hygroscopic growth of Black Carbon       [unitless]
      ! #13 Hygroscopic growth of Organic Carbon     [unitless]
      ! #16 Hygroscopic growth of Sea Salt (accum)   [unitless]
      ! #19 Hygroscopic growth of Sea Salt (coarse)  [unitless]
      !
      ! Computed here:
      ! ---------------------------------
      ! #6  Sulfate Optical Depth (400 nm)           [unitless]
      ! #8  Sulfate Surface Area                     [cm2/cm3 ]
      ! #9  Black Carbon Optical Depth (400 nm)      [unitless]
      ! #11 Black Carbon Surface Area                [cm2/cm3 ]
      ! #12 Organic Carbon Optical Depth (400 nm)    [unitless]
      ! #14 Organic Carbon Surface Area              [cm2/cm3 ]
      ! #15 Sea Salt (accum) Opt Depth (400 nm)      [unitless]
      ! #17 Sea Salt (accum) Surface Area            [cm2/cm3 ]
      ! #18 Sea Salt (coarse) Opt Depth(400 nm)      [unitless]
      ! #20 Sea Salt (coarse) Surface Area           [cm2/cm3 ]
      !
      ! #21: Dust optical depths (0.15 um)     --> from "rdust.f"
      ! #22: Dust optical depths (0.25 um)     --> from "rdust.f"
      ! #23: Dust optical depths (0.4 um)     --> from "rdust.f"
      ! #24: Dust optical depths (0.8 um)     --> from "rdust.f"
      ! #25: Dust optical depths (1.5 um)     --> from "rdust.f"
      ! #26: Dust optical depths (2.5 um)     --> from "rdust.f"
      ! #27: Dust optical depths (4.0 um)     --> from "rdust.f"
      !
      ! #28: Strat. liquid aerosol optical depth     [unitless]
      ! #29: Strat. liquid aerosol surface area      [cm2/cm3 ]
      ! #30: Strat. liquid aerosol number density    [#/cm3   ]
      ! #31: Strat. particulate aerosol opt. depth   [unitless]
      ! #32: Strat. particulate aerosol surface area [cm2/cm3 ]
      ! #33: Strat. particulate aerosol num. density [#/cm3   ]
      !
      ! NOTE: The cloud optical depths are actually recorded at
      !       1000 nm, but vary little with wavelength.
      !==============================================================

      IF ( ND21 > 0 .and. ODSWITCH .EQ. 1 ) THEN

#if defined( UCX )
         ! Get stratospheric aerosol properties from 1st WL (for now)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, W, IOUT, ISTRAT )
!$OMP+SCHEDULE( DYNAMIC )
         !---------------------------------------------------
         ! Strat Aerosol OD [-] and ND [#/cm3]
         ! SAD was stored in previous loop
         !---------------------------------------------------
         DO ISTRAT=1,NSTRATAER
            W = 1
            N = NRHAER + ISTRAT
            IOUT = 5 + 5*NRHAER + 3*NDUST + (ISTRAT-1)*3
            DO L=1,LD21
            DO J=1,JJPAR
            DO I=1,IIPAR
               AD21(I,J,L,IOUT+1) = AD21(I,J,L,IOUT+1) +
     &                              ODAER(I,J,L,IWVSELECT(1,W),N)
               AD21(I,J,L,IOUT+3) = AD21(I,J,L,IOUT+3) +
     &                              NDENS_AER(I,J,L,ISTRAT)*1.d-6
            ENDDO
            ENDDO
            ENDDO
         ENDDO
!$OMP END PARALLEL DO
#endif
         ! Loop over aerosol types (dust handled in dust_mod.f)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, W, LINTERP, IOUT )
!$OMP+SCHEDULE( DYNAMIC )
         DO N = 1, NRHAER 
            !------------------------------------
            ! Aerosol Optical Depths [unitless]
            ! Scale of optical depths w/ RH
            !------------------------------------
            DO W = 1, NWVSELECT 
               !no interpolation required if these are the same
               IF (IWVSELECT(1,W).EQ.IWVSELECT(2,W)) THEN
                  LINTERP=.FALSE.
               ELSE
                  LINTERP=.TRUE.
               ENDIF

               ! Determine subscript for ND21
               ! Store first wavelength in usual position and any
               ! subsequent wavelengths at end
               ! There are always 5 entries first, but separation then
               ! depends on W.
               SELECT CASE(W)
                  CASE(1)
                     IOUT = 5 + (3*(N-1)) + 1
                  CASE(2)
                     IOUT = 5 + (3*NRHAER) + (  NDUST) + N
                  CASE(3)
                     IOUT = 5 + (4*NRHAER) + (2*NDUST) + N
                  CASE DEFAULT
                     IOUT = 5 + (3*(N-1)) + 1
               END SELECT

               DO L = 1, LD21
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  ! Optical Depths
                  IF ( .not. LINTERP ) THEN
                     AD21(I,J,L,IOUT) = AD21(I,J,L,IOUT) +
     &                                  ODAER(I,J,L,IWVSELECT(1,W),N)
                  ELSE
                     ! Interpolated using angstrom exponent between
                     ! Closest available wavelengths
                     ! (coefs pre-calculated in CALC_AOD (RD_AOD.F)
                     !catch any zero values before interpolation
                     IF ((ODAER(I,J,L,IWVSELECT(2,W),N).GT.0).AND.
     &                   (ODAER(I,J,L,IWVSELECT(1,W),N).GT.0)) THEN
                        AD21(I,J,L,IOUT) = AD21(I,J,L,IOUT) +
     &                   ODAER(I,J,L,IWVSELECT(2,W),N)*ACOEF_WV(W)**
     &                   (BCOEF_WV(W)*LOG(ODAER(I,J,L,IWVSELECT(1,W),N)/
     &                   ODAER(I,J,L,IWVSELECT(2,W),N)))
                     ENDIF
                  ENDIF

               ENDDO
               ENDDO
               ENDDO

            ENDDO
         
            !------------------------------------
            ! Aerosol Surface Areas [cm2/cm3]
            !------------------------------------
            DO L = 1, LD21
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Get 3-D grid box indices
               IOUT = 5 + (3*(N-1)) + 3

               ! Add aerosol surface areas
               IF ( L <= LD21 ) THEN
                  AD21(I,J,L,IOUT) = AD21(I,J,L,IOUT) +
     &                               TAREA(I,J,L,N+NDUST)
               ENDIF

            ENDDO
            ENDDO
            ENDDO

         ENDDO
!$OMP END PARALLEL DO
      
      ENDIF 
#if defined( UCX )
      ! Turn off radiative effects of stratospheric aerosols?
      IF (.not.(LSTRATOD)) THEN
         ODAER(:,:,:,:,NRH+1) = 0.d0
         ODAER(:,:,:,:,NRH+2) = 0.d0
      ENDIF
#endif
      !=================================================================
      ! To turn off the radiative effects of different aerososl
      ! uncomment the following lines
      !=================================================================
      !DO R = 1,NRH
      !  ODAER(:,:,:,R)       = 0.d0  !sulfate
      !  ODAER(:,:,:,R+NRH)   = 0.d0  !BC
      !  ODAER(:,:,:,R+2*NRH) = 0.d0  !OC
      !  ODAER(:,:,:,R+3*NRH) = 0.d0  !SS(accum)
      !  ODAER(:,:,:,R+4*NRH) = 0.d0  !SS(coarse)
      !ENDDO
      !ODAER(:,:,:,NRHAER*NRH+1) = 0.d0   !SLA
      !ODAER(:,:,:,NRHAER*NRH+2) = 0.d0   !SPA

      !=================================================================
      ! To turn off heterogeneous chemistry on different aerosols
      ! uncomment the following lines
      !=================================================================
      !TAREA(:,NDUST+1) = 0.d0	!Sulfate
      !TAREA(:,NDUST+2) = 0.d0	!BC 
      !TAREA(:,NDUST+3) = 0.d0	!OC 
      !TAREA(:,NDUST+4) = 0.d0	!SS (accum)
      !TAREA(:,NDUST+5) = 0.d0	!SS (coarse)
      !TAREA(:,NDUST+NRHAER+1) = 0.d0 !SLA
      !TAREA(:,NDUST+NRHAER+2) = 0.d0 !SPA

      ! Free pointers
      NULLIFY( BXHEIGHT, ERADIUS, TAREA, WERADIUS, WTAREA )

      ! Reset first-time flag
      FIRST = .FALSE.

      END SUBROUTINE RDAER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_aerosol
!
! !DESCRIPTION: Subroutine INIT\_AEROSOL allocates and zeroes module arrays
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_AEROSOL( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : Ind_
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  05 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC arguments
!  17 Jun 2016 - R. Yantosca - Now use locally-defined species ID flags,
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_AEROSOL begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

      ! Add tracer ID flags as module variables (bmy, 6/16/16)
      id_BCPI  = Ind_('BCPI' )
      id_BCPO  = Ind_('BCPO' )
      id_DST1  = Ind_('DST1' )
      id_DST2  = Ind_('DST2' )
      id_DST3  = Ind_('DST3' ) 
      id_DST4  = Ind_('DST4' )
      id_DUST1 = Ind_('DUST1')
      id_NH4   = Ind_('NH4'  )
      id_NIT   = Ind_('NIT'  ) 
      id_OCPO  = Ind_('OCPO' ) 
      id_OCPI  = Ind_('OCPI' )  
      id_SALA  = Ind_('SALA' )  
      id_SALC  = Ind_('SALC' )   
      id_SO4   = Ind_('SO4'  )         
      id_SOAG  = Ind_('SOAG' ) 
      id_SOAM  = Ind_('SOAM' ) 
      id_SO4s  = Ind_('SO4s' )  
      id_NITs  = Ind_('NITs' ) 
      id_POA1  = Ind_('POA1' )   
      id_POA2  = Ind_('POA2' )    
      id_OPOA1 = Ind_('OPOA1')   
      id_OPOA2 = Ind_('OPOA2')    
      id_TSOA1 = Ind_('TSOA1') 
      id_TSOA2 = Ind_('TSOA2') 
      id_TSOA3 = Ind_('TSOA3')
      id_TSOA0 = Ind_('TSOA0') 
      id_ISOA1 = Ind_('ISOA1')  
      id_ISOA2 = Ind_('ISOA2')  
      id_ISOA3 = Ind_('ISOA3')
      id_ASOAN = Ind_('ASOAN') 
      id_ASOA1 = Ind_('ASOA1')  
      id_ASOA2 = Ind_('ASOA2') 
      id_ASOA3 = Ind_('ASOA3')

      ! Allocate arrays
      ALLOCATE( BCPI( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCPI' )
      BCPI = 0d0

      ALLOCATE( BCPO( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCPO' )
      BCPO = 0d0

      ALLOCATE( OCPI( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCPI' )
      OCPI = 0d0

      ALLOCATE( OCPO( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCPO' )
      OCPO = 0d0

      ALLOCATE( OCPISOA( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCPISOA' )
      OCPISOA = 0d0

      ALLOCATE( SALA( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SALA' )
      SALA = 0d0

      ALLOCATE( SALC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SALC' )
      SALC = 0d0

      ALLOCATE( SO4_NH4_NIT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4_NH4_NIT' )
      SO4_NH4_NIT = 0d0

      ALLOCATE( SO4( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4' )
      SO4 = 0d0

      ALLOCATE( NH4( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NH4' )
      NH4 = 0d0

      ALLOCATE( NIT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NIT' )
      NIT = 0d0

      ALLOCATE( FRAC_SNA( IIPAR, JJPAR, LLPAR, 3 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRAC_SNA' )
      FRAC_SNA = 0d0

      ALLOCATE( SOILDUST( IIPAR, JJPAR, LLPAR, NDUST ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SOILDUST' )
      SOILDUST = 0d0

#if defined( UCX )
      ALLOCATE( SLA( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SLA' )
      SLA = 0d0

      ALLOCATE( SPA( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SPA' )
      SPA   = 0d0
#endif

      ALLOCATE( TSOA( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TSOA' )
      TSOA = 0d0

      ALLOCATE( ISOA( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ISOA' )
      ISOA = 0d0

      ALLOCATE( ASOA( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ASOA' )
      ASOA = 0d0

      ALLOCATE( OPOA( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OPOA' )
      OPOA = 0d0

      ALLOCATE( SOAG( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SOAG' )
      SOAG = 0d0

      ALLOCATE( SOAM( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SOAM' )
      SOAM = 0d0

      ALLOCATE( PM25( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PM25' )
      PM25 = 0d0

      ! Mass of hydrophobic aerosol from Mian Chin
      ALLOCATE( DAERSL(IIPAR,JJPAR,LLPAR,2), STAT=AS)       
      IF ( AS /= 0) CALL ALLOC_ERR( 'DAERSL' )
      DAERSL = 0d0

      ! Mass of hydrophilic aerosol from Mian Chin
      ALLOCATE( WAERSL(IIPAR,JJPAR,LLPAR,NAER), STAT=AS)       
      IF ( AS /= 0) CALL ALLOC_ERR( 'WAERSL' )
      DAERSL = 0d0

      END SUBROUTINE INIT_AEROSOL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_aerosol
!
! !DESCRIPTION: Subroutine CLEANUP\_AEROSOL deallocates all module arrays
!  (bmy, 7/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_AEROSOL
! 
! !REVISION HISTORY:
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!

      !=================================================================
      ! CLEANUP_AEROSOL begins here!
      !=================================================================
      IF ( ALLOCATED( BCPI        ) ) DEALLOCATE( BCPI        )
      IF ( ALLOCATED( BCPO        ) ) DEALLOCATE( BCPO        )
      IF ( ALLOCATED( OCPI        ) ) DEALLOCATE( OCPI        )
      IF ( ALLOCATED( OCPO        ) ) DEALLOCATE( OCPO        )
      IF ( ALLOCATED( OCPISOA     ) ) DEALLOCATE( OCPISOA     )
      IF ( ALLOCATED( SALA        ) ) DEALLOCATE( SALA        )
      IF ( ALLOCATED( SALC        ) ) DEALLOCATE( SALC        )
      IF ( ALLOCATED( SO4_NH4_NIT ) ) DEALLOCATE( SO4_NH4_NIT )
      IF ( ALLOCATED( SO4         ) ) DEALLOCATE( SO4         )
      IF ( ALLOCATED( NH4         ) ) DEALLOCATE( NH4         )
      IF ( ALLOCATED( NIT         ) ) DEALLOCATE( NIT         )
      IF ( ALLOCATED( FRAC_SNA    ) ) DEALLOCATE( FRAC_SNA    )
      IF ( ALLOCATED( SOILDUST    ) ) DEALLOCATE( SOILDUST    )
      IF ( ALLOCATED( SLA         ) ) DEALLOCATE( SLA         )
      IF ( ALLOCATED( SPA         ) ) DEALLOCATE( SPA         )
      IF ( ALLOCATED( TSOA        ) ) DEALLOCATE( TSOA        )
      IF ( ALLOCATED( ISOA        ) ) DEALLOCATE( ISOA        )
      IF ( ALLOCATED( ASOA        ) ) DEALLOCATE( ASOA        )
      IF ( ALLOCATED( OPOA        ) ) DEALLOCATE( OPOA        )
      IF ( ALLOCATED( SOAG        ) ) DEALLOCATE( SOAG        )
      IF ( ALLOCATED( SOAM        ) ) DEALLOCATE( SOAM        )
      IF ( ALLOCATED( PM25        ) ) DEALLOCATE( PM25        )

      END SUBROUTINE CLEANUP_AEROSOL
!EOC
      END MODULE AEROSOL_MOD
