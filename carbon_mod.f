! $Id: carbon_mod.f,v 1.1 2004/04/13 14:52:57 bmy Exp $
      MODULE CARBON_MOD
!
!******************************************************************************
!  Module CARBON_MOD contains arrays and routines for performing an offline 
!  carbonaceous aerosol simulation.  Original code taken from Mian Chin's 
!  GOCART model and modified accordingly. (rjp, bmy, 4/2/04)
!
!  4 Aerosol species : Organic and Black carbon 
!                    : hydrophilic (soluble) and hydrophobic of each
!
!  Module Variables:
!  ============================================================================
!  (1 ) ANTH_BLKC        (REAL*8 ) : BC anthropogenic emissions       [kg C ]
!  (2 ) ANTH_ORGC        (REAL*8 ) : OC anthropogenic emissions       [kg C ]
!  (3 ) BIOB_BLKC        (REAL*8 ) : BC biomass emissions             [kg C ]
!  (4 ) BIOB_ORGC        (REAL*8 ) : OC biomass emissions             [kg C ]
!  (5 ) BIOF_BLKC        (REAL*8 ) : BC biofuel emissions             [kg C ]
!  (6 ) BIOF_ORGC        (REAL*8 ) : OC biofuel emissions             [kg C ]
!  (7 ) EF_BLKC          (REAL*8 ) : Emission factors for BC          [kg/kg]
!  (8 ) EF_ORGC          (REAL*8 ) : Emission factors for OC          [kg/kg]
!  (9 ) TERP_ORGC        (REAL*8 ) : Secodary organic emissions       [kg C ]
!  (10) BCCONV           (REAL*8 ) : Hydrophilic BC from Hydrophobic  [kg C ]
!  (11) OCCONV           (REAL*8 ) : Hydrophilic OC from Hydrophobic  [kg C ]
!  (13) USE_MONTHLY_ANTH (LOGICAL) : Toggles monthly or annual anthro emissions
!  (14) USE_MONTHLY_BIOB (LOGICAL) : Toggles monthly or annual biomass emiss.
!
!  Module Routines:
!  ============================================================================
!  (1 ) CHEMCARBON         : Driver program for carbon aerosol chemistry
!  (2 ) CHEM_BCPO          : Chemistry routine for hydrophobic BC (aka EC)
!  (3 ) CHEM_BCPI          : Chemistry routine for hydrophilic BC (aka EC)
!  (4 ) CHEM_OCPO          : Chemistry routine for hydrophobic OC
!  (5 ) CHEM_OCPI          : Chemistry routine for hydrophilic OC
!  (6 ) EMISSCARBON        : Driver routine for carbon aerosol emissions
!  (7 ) BIOGENIC_OC        : Computes biogenic OC [each time step]
!  (8 ) ANTHRO_CARB_TBOND  : Computes anthropogenic OC/EC [annual data]
!  (9 ) ANTHRO_CARB_COOKE  : Computes anthropogenic OC/EC [monthly data]
!  (10) BIOMASS_CARB_TBOND : Computes biomass burning OC/EC [annual data]
!  (11) BIOMASS_CARB_GEOS  : Computes biomass burning OC/EC [monthly data]
!  (12) EMITHIGH           : Computes complete mixing of emission within PBL
!
!  NOTE: Choose either (8 ) or (9 ) for ANTHROPOGENIC emission
!        Choose either (11) or (12) for BIOMASS BURNING emission.
!
!  GEOS-CHEM modules referenced by carbon_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f    : Module containing routines for binary punch file I/O
!  (2 ) dao_mod.f      : Module containing arrays for DAO met fields
!  (3 ) diag_mod.f     : Module containing GEOS-CHEM diagnostic arrays
!  (4 ) drydep_mod.f   : 
!  (4 ) error_mod.f    : Module containing I/O error and NaN check routines
!  (5 ) grid_mod.f     : Module containing horizontal grid information
!  (  ) pressure_mod.f : Module containing routines to compute P(I,J,L)
!  (6 ) time_mod.f     : Module containing routines for computing time & date
!  (7 ) tracerid_mod.f : Module containing pointers to tracers & emissions
!  (8 ) transfer_mod.f : Module containing routines to cast & resize arrays
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "biomass_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE             :: NDAYS,     ANTH_BLKC, ANTH_ORGC
      PRIVATE             :: BIOB_BLKC, BIOB_ORGC, BIOF_BLKC 
      PRIVATE             :: BIOF_ORGC, EF_BLKC,   EF_ORGC   
      PRIVATE             :: TERP_ORGC, BCCONV,    OCCONV 
      PRIVATE             :: SMALLNUM,  DRYBCPI,   DRYOCPI
      PRIVATE             :: DRYBCPO,   DRYOCPO
      PRIVATE             :: USE_MONTHLY_ANTH
      PRIVATE             :: USE_MONTHLY_BIOB

      ! PRIVATE module routines
      PRIVATE             :: CHEM_BCPO          
      PRIVATE             :: CHEM_BCPI           
      PRIVATE             :: CHEM_OCPO           
      PRIVATE             :: CHEM_OCPI           
      PRIVATE             :: BIOGENIC_OC         
      PRIVATE             :: ANTHRO_CARB_TBOND   
      PRIVATE             :: ANTHRO_CARB_COOKE   
      PRIVATE             :: BIOMASS_CARB_TBOND  
      PRIVATE             :: BIOMASS_CARB_GEOS   
      PRIVATE             :: EMITHIGH            

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      LOGICAL             :: USE_MONTHLY_ANTH = .FALSE.
      LOGICAL             :: USE_MONTHLY_BIOB = .TRUE.
      INTEGER             :: DRYBCPI, DRYOCPI, DRYBCPO, DRYOCPO
      REAL*8, PARAMETER   :: SMALLNUM = 1d-20
      REAL*8, ALLOCATABLE :: ANTH_BLKC(:,:,:)
      REAL*8, ALLOCATABLE :: ANTH_ORGC(:,:,:)
      REAL*8, ALLOCATABLE :: BIOB_BLKC(:,:,:)
      REAL*8, ALLOCATABLE :: BIOB_ORGC(:,:,:)
      REAL*8, ALLOCATABLE :: BIOF_BLKC(:,:,:)
      REAL*8, ALLOCATABLE :: BIOF_ORGC(:,:,:)
      REAL*8, ALLOCATABLE :: EF_BLKC(:,:)
      REAL*8, ALLOCATABLE :: EF_ORGC(:,:)
      REAL*8, ALLOCATABLE :: TERP_ORGC(:,:)
      REAL*8, ALLOCATABLE :: BCCONV(:,:,:)
      REAL*8, ALLOCATABLE :: OCCONV(:,:,:)

      ! Days per month (based on 1998)
      INTEGER             :: NDAYS(12) = (/ 31, 28, 31, 30, 31, 30, 
     &                                      31, 31, 30, 31, 30, 31 /)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE CHEMCARBON
!
!******************************************************************************
!  Subroutine CHEMCARBON is the interface between the GEOS-CHEM main 
!  program and the carbon aerosol chemistry routines that calculates
!  dry deposition and chemical conversion between hydrophilic and 
!  hydrophobic.
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DRYDEP_MOD,   ONLY : DEPNAME, NUMDEP
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE TRACERID_MOD

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! STT, etc

      ! Local variables
      LOGICAL, SAVE :: FIRSTCHEM = .TRUE.
      INTEGER       :: N

      !=================================================================
      ! CHEMCARBON begins here!
      !=================================================================

      ! Establish indices w/in DEPSAV array
      IF ( FIRSTCHEM ) THEN

         ! Initialize arrays (if not already done before)
         CALL INIT_CARBON

         ! Find drydep species in DEPSAV
         DO N = 1, NUMDEP
            SELECT CASE ( TRIM( DEPNAME(N) ) )
               CASE ( 'BCPI' )
                  DRYBCPI = N
               CASE ( 'OCPI' )
                  DRYOCPI = N
               CASE ( 'BCPO' )
                  DRYBCPO = N
               CASE ( 'OCPO' )
                  DRYOCPO = N
               CASE DEFAULT
                  ! Nothing
            END SELECT        
         ENDDO

         ! Reset first-time flag
         FIRSTCHEM = .FALSE.
      ENDIF

      ! Chemistry for hydrophobic BC
      IF ( IDTBCPO > 0 ) THEN
         CALL CHEM_BCPO( STT(:,:,:,IDTBCPO) )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_BCPO' )
      ENDIF

      ! Chemistry for hydrophilic BC
      IF ( IDTBCPI > 0 ) THEN
         CALL CHEM_BCPI( STT(:,:,:,IDTBCPI) )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_BCPI' )
      ENDIF

      ! Chemistry for hydrophobic OC
      IF ( IDTOCPO > 0 ) THEN
         CALL CHEM_OCPO( STT(:,:,:,IDTOCPO) )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_OCPO' )
      ENDIF

      ! Chemistry for hydrophilic OC
      IF ( IDTOCPI > 0 ) THEN 
         CALL CHEM_OCPI( STT(:,:,:,IDTOCPI) )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMCARBON: a CHEM_OCPI' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEMCARBON

!-----------------------------------------------------------------------------

      SUBROUTINE CHEM_BCPO( TC )
!
!******************************************************************************
!  Subroutine CHEM_BCPO converts hydrophobic BC to hydrophilic BC and
!  calculates the dry deposition of hydrophobic BC. (rjp, bmy, 4/1/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8) : Array of hydrophobic BC tracer 
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44, AD07_BC 
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACERID_MOD, ONLY : IDTBCPO
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! NCHEM, DXYP
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44, ND07, LD07

      ! Arguments
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I,       J,   L
      REAL*8                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)     
      REAL*8                 :: DTCHEM, FLUX, KBC, FREQ
      REAL*8                 :: TC0,    CNEW, RKT, AREA_CM2, BL_FRAC
      REAL*8,  PARAMETER     :: BC_LIFE = 1.15D0


      !=================================================================
      ! CHEM_BCPO begins here!
      !=================================================================

      ! Return if BCPO isn't defined
      IF ( IDTBCPO == 0 .or. DRYBCPO == 0 ) RETURN

      ! Initialize
      KBC    = 1.D0 / ( 86400d0 * BC_LIFE )
      DTCHEM = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         BCCONV(I,J,L) = 0d0

         ! Initialize for drydep diagnostic
         IF ( ND44 > 0 ) ND44_TMP(I,J,L) = 0d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! For tracers with dry deposition, the loss rate of dry dep is 
      ! combined in chem loss term.
      !
      ! Conversion from hydrophobic to hydrophilic:  
      ! e-folding time 1.15 days 
      ! ----------------------------------------
      ! Use an e-folding time of 1.15 days or a convertion rate 
      ! of 1.0e-5 /sec. 
      !
      ! Hydrophobic(2) --> Hydrophilic(1) ,  k  = 1.0e-5          
      ! Both aerosols are dry-deposited,     kd = Dvel/DELZ (sec-1)      
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, FREQ, BL_FRAC, RKT, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial BC mass [kg]
         TC0  = TC(I,J,L)

         ! Zero drydep freq
         FREQ = 0d0

         ! PBLFRAC is only defined up to the tropopause, but we need to 
         ! do the conversion from H-philic to H-phobic at all levels
         IF ( L > LLTROP ) THEN
            BL_FRAC = 0d0
         ELSE
            BL_FRAC = PBLFRAC(I,J,L)
         ENDIF

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! BC drydep frequency [1/s] -- PBLFRAC accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYBCPO) * BL_FRAC

         ENDIF

         ! Amount of BCPO left after chemistry and drydep [kg]
         RKT  = ( KBC + FREQ ) * DTCHEM
         CNEW = TC0 * EXP( -RKT )

         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Amount of BCPO converted to BCPI [kg/timestep]
         BCCONV(I,J,L) = ( TC0 - CNEW ) * KBC / ( KBC + FREQ )

         !==============================================================
         ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
         !==============================================================
         IF ( ND44 > 0 .AND. FREQ > 0d0 ) THEN

             ! Surface area [cm2]
             AREA_CM2 = GET_AREA_CM2( J )

             ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]  
             ! XNUMOL is the ratio [molec tracer/kg tracer]   
             FLUX     = TC0 - CNEW - BCCONV(I,J,L) 
             FLUX     = FLUX * XNUMOL(IDTBCPO) / ( DTCHEM * AREA_CM2 )

             ! Store in ND44_TMP as a placeholder
             ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
         ENDIF

         !==============================================================
         ! ND07 diagnostic: H-philic BC from H_phobic BC [kg/timestep]
         !==============================================================
         IF ( ND07 > 0 .and. L <= LD07 ) THEN
             AD07_BC(I,J,L) = AD07_BC(I,J,L) + BCCONV(I,J,L)
         ENDIF

         ! Store new concentration back into tracer array
         TC(I,J,L) = CNEW
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !===============================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !===============================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYBCPO,1) = AD44(I,J,DRYBCPO,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF      

      ! Return to calling program
      END SUBROUTINE CHEM_BCPO

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_BCPI( TC )
!
!******************************************************************************
!  Subroutine CHEM_BCPI calculates dry deposition of hydrophilic BC.
!  (rjp, bmy, 4/1/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8) : Array of hydrophilic BC tracer 
! 
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44 
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACERID_MOD, ONLY : IDTBCPI
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! NCHEM, DXYP
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I, J, L
      REAL*8                 :: DTCHEM, FLUX, BL_FRAC, AREA_CM2
      REAL*8                 :: TC0,    CNEW, CCV,     FREQ
      REAL*8                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! CHEM_BCPI begins here!
      !=================================================================

      ! Return if BCPI isn't defined
      IF ( IDTBCPI == 0 .or. DRYBCPI == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize for ND44 diagnostic
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ND44_TMP(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, CCV, FREQ, BL_FRAC, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial H-philic BC [kg]
         TC0 = TC(I,J,L)

         ! H-philic BC that used to be H-phobic BC [kg]
         CCV = BCCONV(I,J,L)
         
         ! PBLFRAC is only defined up to the tropopause, but we need to 
         ! do the conversion from H-philic to H-phobic everywhere
         IF ( L > LLTROP ) THEN
            BL_FRAC = 0d0
         ELSE
            BL_FRAC = PBLFRAC(I,J,L)
         ENDIF

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! Drydep frequency
            FREQ = DEPSAV(I,J,DRYBCPI) * BL_FRAC
            
            !===========================================================
            ! Note, This is an analytical solution of first order 
            ! partial differential equations (w/ 2 solutions):
            !
            ! #1) CNEW = Cphi * exp(-RKT) + Cconv/RKT * (1.-exp(-RKT)) 
            ! #2) CNEW = ( Cphi + Cconv ) * exp(-RKT)
            !===========================================================

            ! Comment out for now
            !CNEW = TC0 * EXP( -FREQ * DTCHEM ) 
            !     + CCV / FREQ * ( 1.D0 - EXP( -FREQ * DTCHEM ) )

            ! Amount of BCPI left after drydep [kg]
            CNEW = ( TC0 + CCV ) * EXP( -FREQ * DTCHEM )

            !===========================================================
            ! ND44 diagnostic: drydep flux [atoms C/cm2/s]
            !===========================================================
            IF ( ND44 > 0 .and. FREQ > 0d0 ) THEN
  
               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( J )

               ! Convert drydep loss from [kg/timestep] to [molec/cm2/s]
               FLUX = ( TC0 + CCV - CNEW ) 
               FLUX = FLUX * XNUMOL(IDTBCPI) / ( AREA_CM2 * DTCHEM )
             
               ! Store in ND44_TMP as a placeholder
               ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
            ENDIF

         ELSE

            ! Otherwise, omit the exponential to save on clock cycles
            CNEW = TC0 + CCV

         ENDIF
      
         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Save new concentration of H-philic IC in tracer array
         TC(I,J,L) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !=================================================================
      ! Zero out the BCCONV array for the next iteration
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         BCCONV(I,J,L) = 0.d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !=================================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYBCPI,1) = AD44(I,J,DRYBCPI,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF      

      ! Return to calling program
      END SUBROUTINE CHEM_BCPI

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_OCPO( TC )
!
!******************************************************************************
!  Subroutine CHEM_OCPO converts hydrophobic OC to hydrophilic OC and
!  calculates the dry deposition of hydrophobic OC. (rjp, bmy, 4/1/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8) : Array of hydrophobic OC tracer [kg]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44, AD07_OC 
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACERID_MOD, ONLY : IDTOCPO
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! NCHEM, DXYP
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44, ND07, LD07

      ! Arguments
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variable
      INTEGER                :: I, J, L
      REAL*8                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)
      REAL*8                 :: DTCHEM, FLUX, KOC,  BL_FRAC
      REAL*8                 :: TC0,    FREQ, CNEW, RKT, AREA_CM2
      REAL*8,  PARAMETER     :: OC_LIFE = 1.15D0

      !=================================================================
      ! CHEM_OCPO begins here!
      !=================================================================

      ! Return if OCPO isn't defined
      IF ( IDTOCPO == 0 .or. DRYOCPO == 0 ) RETURN

      ! Initialize
      KOC    = 1.D0 / ( 86400d0 * OC_LIFE )
      DTCHEM = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         OCCONV(I,J,L) = 0d0

         ! Initialize for drydep diagnostic
         IF ( ND44 > 0 ) ND44_TMP(I,J,L) = 0d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! For tracers with dry deposition, the loss rate of dry dep is 
      ! combined in chem loss term.
      !
      ! Conversion from hydrophobic to hydrophilic:  
      ! e-folding time 1.15 days 
      ! ----------------------------------------
      ! Use an e-folding time of 1.15 days or a convertion rate 
      ! of 1.0e-5 /sec. 
      !    Hydrophobic --> Hydrophilic,  k  = 1.0e-5          
      !    Aerosols are dry-deposited,   kd = DEPSAV (sec-1)      
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, FREQ, BL_FRAC, RKT, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial OC [kg]
         TC0  = TC(I,J,L)

         ! Zero drydep freq 
         FREQ = 0d0

         ! PBLFRAC is only defined up to the tropopause, but we need to 
         ! do the conversion from H-philic to H-phobic everywhere
         IF ( L > LLTROP ) THEN
            BL_FRAC = 0d0
         ELSE
            BL_FRAC = PBLFRAC(I,J,L)
         ENDIF

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! OC drydep frequency [1/s] -- PBLFRAC accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYOCPO) * BL_FRAC

         ENDIF

         ! Amount of OCPO left after chemistry and drydep [kg]
         RKT  = ( KOC + FREQ ) * DTCHEM
         CNEW = TC0 * EXP( -RKT )

         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Amount of OCPO converted to OCPI [kg/timestep]
         OCCONV(I,J,L) = ( TC0 - CNEW ) * KOC / ( KOC + FREQ )

         !==============================================================
         ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
         !==============================================================
         IF ( ND44 > 0 .AND. FREQ > 0d0 ) THEN

             ! Surface area [cm2]
             AREA_CM2 = GET_AREA_CM2( J )

             ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]
             ! XNUMOL is the ratio [molec tracer/kg tracer]     
             FLUX     = TC0 - CNEW - OCCONV(I,J,L)
             FLUX     = FLUX * XNUMOL(IDTOCPO) / ( DTCHEM * AREA_CM2 )

             ! Store in ND44_TMP as a placeholder
             ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
         ENDIF

         !==============================================================
         ! ND07 diagnostic: H-Philic OC from H-phobic [kg/timestep]
         !==============================================================
         IF ( ND07 > 0 .and. L <= LD07 ) THEN
            AD07_OC(I,J,L) = AD07_OC(I,J,L) + OCCONV(I,J,L) 
         ENDIF

         ! Store modified OC concentration back in tracer array
         TC(I,J,L) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !=================================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !=================================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYOCPO,1) = AD44(I,J,DRYOCPO,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF   

      ! Return to calling program
      END SUBROUTINE CHEM_OCPO

!-----------------------------------------------------------------------

      SUBROUTINE CHEM_OCPI( TC )
!
!******************************************************************************
!  Subroutine CHEM_BCPI calculates dry deposition of hydrophilic OC.
!  (rjp, bmy, 4/1/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC (REAL*8) : Array of hydrophilic BC tracer 
! 
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44 
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACERID_MOD, ONLY : IDTOCPI
      USE TIME_MOD,     ONLY : GET_TS_CHEM

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! NCHEM, DXYP
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variable
      INTEGER                :: I, J, L
      REAL*8                 :: DTCHEM, FLUX, BL_FRAC
      REAL*8                 :: TC0, CNEW, CCV, FREQ, AREA_CM2
      REAL*8                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! CHEM_OCPI begins here!
      !=================================================================
      IF ( IDTOCPI == 0 .or. DRYOCPI == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize for drydep diagnostic
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ND44_TMP(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, CCV, FREQ, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initial H-philic OC [kg]
         TC0 = TC(I,J,L)

         ! H-philic OC that used to be H-phobic OC [kg]
         CCV = OCCONV(I,J,L)

         ! PBLFRAC is only defined up to the tropopause, but we need to 
         ! do the conversion from H-philic to H-phobic everywhere
         IF ( L > LLTROP ) THEN
            BL_FRAC = 0d0
         ELSE
            BL_FRAC = PBLFRAC(I,J,L)
         ENDIF

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! Drydep frequency [1/s]
            FREQ = DEPSAV(I,J,DRYOCPI) * BL_FRAC

            !===========================================================
            ! Note, This is an analytical solution of first order 
            ! partial differential equations (w/ 2 solutions):
            !
            ! #1) CNEW = Cphi * exp(-RKT) + Cconv/RKT * (1.-exp(-RKT))
            ! #2) CNEW = ( Cphi + Cconv ) * exp(-RKT)
            !===========================================================

            ! CNEW = TC0 * EXP( -FREQ * DTCHEM ) 
            !       + CCV / FREQ * ( 1.D0 - EXP( -FREQ * DTCHEM ) )

            ! Amount of BCPI left after drydep [kg]
            CNEW = ( TC0 + CCV ) * EXP( -FREQ * DTCHEM )

            !===========================================================
            ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
            !===========================================================
            IF ( ND44 > 0 ) THEN

               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( J )

               ! Convert drydep loss from [kg/timestep] to [atoms C/cm2/s]
               FLUX = ( TC0 + CCV - CNEW ) 
               FLUX = FLUX * XNUMOL(IDTOCPI) / ( AREA_CM2 * DTCHEM )
             
               ! Store in ND44_TMP as a placeholder
               ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
            ENDIF

         ELSE

            ! Otherwise, avoid doing the exponential
            ! to preserve precision and clock cycles
            CNEW = TC0 + CCV

         ENDIF
      
         ! Prevent underflow condition
         IF ( CNEW < SMALLNUM ) CNEW = 0d0

         ! Store modified concentration back in tracer array [kg]
         TC(I,J,L) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  

      !=================================================================
      ! Zero OCCONV array for next timestep
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         OCCONV(I,J,L) = 0d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


      !=================================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !=================================================================
      IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DRYOCPI,1) = AD44(I,J,DRYOCPI,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF    

      ! Return to calling program
      END SUBROUTINE CHEM_OCPI

!-----------------------------------------------------------------------

      SUBROUTINE EMISSCARBON
!
!=========================================================================
!  Subroutine EMISSCARBON is the interface between the GEOS-CHEM model
!  and the CARBONACEOUS AEROSOL emissions (rjp, 1/24/02)
!
!==========================================================================
!
      ! References to F90 modules
      USE DIAG_MOD,    ONLY : AD07
      USE DAO_MOD,     ONLY : PBL
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE TIME_MOD,    ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TRACERID_MOD

#     include "CMN_SIZE"    ! Size paramters
#     include "CMN"         ! STT
#     include "CMN_DIAG"    ! ND07

      ! Local variables
      LOGICAL, SAVE        :: FIRST = .TRUE.
      INTEGER              :: I, J, MONTH, N
      REAL*8               :: BCSRC(IIPAR,JJPAR,2)
      REAL*8               :: OCSRC(IIPAR,JJPAR,2)

      !=================================================================
      ! EMISSCARBON begins here!
      !
      ! Read carbonaceous aerosols from disk and compute hydrophilic 
      ! and hydrophobic fractions. NOTE, CARBON AEROSOLS HAVE TO BE 
      ! ORDERED AS Hydrophilic(BC[1], OC[2]) Hydrophobic(BC[3], OC[4]).
      !=================================================================      

      !--------------------------
      ! Read time-invariant data
      !--------------------------
      IF ( FIRST ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 100 )

         ! Monthly or annual ANTHRO emissions?
         IF ( USE_MONTHLY_ANTH ) THEN
            WRITE( 6, 110 )
            WRITE( 6, 111 )
         ELSE
            WRITE( 6, 120 )
         ENDIF

         ! Monthly or annual BIOMASS emissions?
         IF ( USE_MONTHLY_BIOB ) THEN
            WRITE( 6, 130 )
         ELSE
            WRITE( 6, 140 )
         ENDIF
         
         ! Write spacer
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! FORMAT strings
 100     FORMAT( 'C A R B O N   A E R O S O L   E M I S S I O N S'    )
 110     FORMAT( 'w/ ANTHROPOGENIC emissions from Cooke et al [1999]' )
 111     FORMAT( ' and seasonality following Park et al [2003]',      )
 120     FORMAT( 'w/ ANTHROPOGENIC emissions from Bond et al [2004]'  )
 130     FORMAT( 'w/ BIOMASS emissions from GEOS-CHEM inventory'      )
 140     FORMAT( 'w/ BIOMASS emissions from Bond et al [2004]'        )

         ! Initialize arrays
         CALL INIT_CARBON       

         ! Read annual mean anthro emissions
         IF ( .not. USE_MONTHLY_ANTH ) CALL ANTHRO_CARB_TBOND

         ! Read annual mean biomass emissions
         IF ( .not. USE_MONTHLY_BIOB ) CALL BIOMASS_CARB_TBOND

         ! Reset flag
         FIRST = .FALSE.
      ENDIF
      
      !--------------------------
      ! Read monthly-mean data
      !--------------------------
      IF ( ITS_A_NEW_MONTH() ) THEN
      
         ! Current month
         MONTH = GET_MONTH()

         ! Read monthly mean anthro emissions
         IF ( USE_MONTHLY_ANTH ) CALL ANTHRO_CARB_COOKE( MONTH )

         ! Read monthly mean biomass emissions
         IF ( USE_MONTHLY_BIOB ) CALL BIOMASS_CARB_GEOS( MONTH )
      ENDIF
      
      !--------------------------
      ! Compute biogenic OC
      !--------------------------
      CALL BIOGENIC_OC
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISCARBON: after data' )

      !=================================================================
      ! Sum up BC and OC sources. 
      ! N=1 is HYDROPHILIC; N=2 is HYDROPHOBIC.
      !
      ! COMMENT: Maybe someday we'll want to play with the different 
      ! emission height for different source type.  For example the
      ! carbon from biomass burning could be emitted to the higher 
      ! altitude due to the thermal bouyancy and shallow convection.
      ! The current setting to use EMITHIGH seems rather inefficient 
      ! but robust for sensitivity studies for emission height 
      ! variation on carbon concentrations, so please keep using the 
      ! current setup until we decide otherwise. (rjp, 4/2/02)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Total HYDROPHILIC BC source [kg]
         BCSRC(I,J,1) = ANTH_BLKC(I,J,1) + 
     &                  BIOF_BLKC(I,J,1) + 
     &                  BIOB_BLKC(I,J,1)   

         ! Total HYDROPHOBIC BC source [kg]
         BCSRC(I,J,2) = ANTH_BLKC(I,J,2) +
     &                  BIOF_BLKC(I,J,2) +
     &                  BIOB_BLKC(I,J,2)  
 
         ! Total HYDROPHILIC OC source [kg]
         OCSRC(I,J,1) = ANTH_ORGC(I,J,1) + 
     &                  BIOF_ORGC(I,J,1) + 
     &                  BIOB_ORGC(I,J,1) + 
     &                  TERP_ORGC(I,J)

         ! Total HYDROPHOBIC OC source [kg]
         OCSRC(I,J,2) = ANTH_ORGC(I,J,2) + 
     &                  BIOF_ORGC(I,J,2) + 
     &                  BIOB_ORGC(I,J,2) 
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Sum up all carbon tracers throughout the boundary layer
      CALL EMITHIGH( BCSRC, OCSRC )
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISCARBON: after EMITHIGH' )

      !=================================================================
      ! ND07 diagnostic: Carbon aerosol emissions [kg/timestep]
      !=================================================================
      IF ( ND07 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Anthropogenic BC source
            AD07(I,J,1) = AD07(I,J,1)        +
     &                    ( ANTH_BLKC(I,J,1) + ANTH_BLKC(I,J,2) )
            
            ! Biogenic BC source
            AD07(I,J,2) = AD07(I,J,2)        +
     &                    ( BIOB_BLKC(I,J,1) + BIOB_BLKC(I,J,2) )

            ! Biofuel BC source
            AD07(I,J,3) = AD07(I,J,3)        +
     &                    ( BIOF_BLKC(I,J,1) + BIOF_BLKC(I,J,2) )

            ! Anthropogenic OC source
            AD07(I,J,4) = AD07(I,J,4)        +
     &                    ( ANTH_ORGC(I,J,1) + ANTH_ORGC(I,J,2) )

            ! Biomass OC source
            AD07(I,J,5) = AD07(I,J,5)        +
     &                    ( BIOB_ORGC(I,J,1) + BIOB_ORGC(I,J,2) )

            ! Biofuel OC source
            AD07(I,J,6) = AD07(I,J,6)        + 
     &                    ( BIOF_ORGC(I,J,1) + BIOF_ORGC(I,J,2) )

            ! Terpene source
            AD07(I,J,7) = AD07(I,J,7)        + TERP_ORGC(I,J)
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISCARBON: after ND07' )


      ! Return to calling program
      END SUBROUTINE EMISSCARBON

!------------------------------------------------------------------------------

      SUBROUTINE BIOGENIC_OC
!
!******************************************************************************
!  Subroutine BIOGENIC_OC emits secondary organic carbon aerosols.
!  (rjp, bmy, 4/1/04)
!
!  Terpene emissions as a source of OC:  TERP.GEIA90.a1.2x2.5.*
!  Assuming 10% yield of OC(hydrophilic) from terpene emission.
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_SETUP"    ! DATA_DIR

      ! Local variables
      LOGICAL, SAVE         :: FIRSTEMIS = .TRUE.
      INTEGER               :: I, J, IJLOOP
      REAL*8                :: CONVERT(NVEGTYPE)
      REAL*8                :: GMONOT(NVEGTYPE)
      REAL*8                :: TMMP, EMMO

      ! Fraction of yield of OC (hydrophilic) from terpene emission
      REAL*8, PARAMETER     :: FBIOG = 1.0d-1

      ! External functions
      REAL*8,  EXTERNAL     :: XLTMMP
      REAL*8,  EXTERNAL     :: EMMONOT

      !=================================================================
      ! BIOGENIC_OC begins here!
      !=================================================================

      ! Get ISOPRENE baseline emissions (first-time only)
      IF ( FIRSTEMIS ) THEN
         CALL RDISOPT ( CONVERT )
         CALL RDMONOT ( GMONOT  )
         CALL SETBASE ( CONVERT, GMONOT )

         ! Reset first-time flag
         FIRSTEMIS = .FALSE.
      ENDIF

      ! 1-D loop index
      IJLOOP = 0

      ! Loop over surface gridboxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! 1-D loop index
         IJLOOP         = IJLOOP + 1

         ! Surface temperature [K]
         TMMP           = XLTMMP(I,J,IJLOOP)

         ! EMMO = [kg C/box/time-step] from monoterpenes
         EMMO           = EMMONOT( IJLOOP, TMMP, 1.d0 )

         ! Fraction of EMMO that converts into OC [kg/box/timestep]
         TERP_ORGC(I,J) = EMMO * FBIOG
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE BIOGENIC_OC

!------------------------------------------------------------------------------

      SUBROUTINE ANTHRO_CARB_TBOND
!
!******************************************************************************
!  Subroutine ANTHRO_CARB_TBOND computes annual mean anthropogenic and 
!  biofuel emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC 
!  CARBON.  It also separates these into HYDROPHILIC and HYDROPHOBIC 
!  fractions. (rjp, bmy, 4/2/04)
!
!  Emissions data comes from the Bond et al [2004] inventory and has units
!  of [kg C/yr].  This will be converted to [kg C/timestep] below.
!
!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_SETUP"  ! DATA_DIR

      ! Local variables
      INTEGER             :: I, J
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU, STEPS_PER_YR
      REAL*8              :: FD2D(IIPAR,JJPAR)
      CHARACTER(LEN=255)  :: FILENAME

      ! Hydrophilic fraction of BLACK CARBON (aka ELEMENTAL CARBON)
      REAL*8, PARAMETER   :: FHB = 0.2d0

      ! Hydrophilic fraction of ORGANIC CARBON 
      REAL*8, PARAMETER   :: FHO = 0.5d0 

      !=================================================================
      ! ANTHRO_CARB_TBOND begins here!
      !=================================================================

      ! Number of emission timesteps per year
      STEPS_PER_YR = ( ( 1440 * 365 ) / GET_TS_EMIS() )

      ! Get TAU0 value to index the punch file
      XTAU         = GET_TAU0( 1, 1, 2001 )

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from 
      ! anthropogenic sources as tracer #34 in [kg C/year].  
      ! Then convert to [kg C/timestep] and store in ANTH_BLKC.
      !=================================================================

      ! Filename for carbon aerosol from fossil fuel use
      FILENAME = TRIM( DATA_DIR )                        // 
     &           'carbon_200404/BCOC_TBond_fossil.geos.' // 
     &           GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - ANTHRO_CARB_TBOND: Reading ', a )

      ! Read BLCK emission
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Hydrophilic BLACK CARBON from anthropogenics [kg C/timestep]
         ANTH_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_YR
         
         ! Hydrophobic BLACK CARBON from anthropogenics [kg C/timestep]
         ANTH_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_YR
        
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read ORGANIC CARBON from anthropogenic sources as tracer #35
      ! in [kg C/year].  Then Convert to [kg C/timestep] and store in 
      ! ANTH_ORGC.
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Hydrophilic ORGANIC CARBON from anthropogenics [kg C/timestep]
         ANTH_ORGC(I,J,1) =          FHO *   FD2D(I,J) / STEPS_PER_YR

         ! Hydrophobic ORGANIC CARBON from anthropogenics [kgC/timestep]
         ANTH_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_YR

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from biofuel 
      ! combustion as tracer #34 in [kg C/year].  Then convert to 
      ! [kg C/timestep] and store in BIOF_BLKC.
      !=================================================================

      ! Filename
      FILENAME = TRIM( DATA_DIR )                         // 
     &           'carbon_200404/BCOC_TBond_biofuel.geos.' // 
     &           GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Hydrophilic BLACK CARBON from biofuels [kg C /timestep]
         BIOF_BLKC(I,J,1) =          FHB *   FD2D(I,J) / STEPS_PER_YR
         
         ! Hydrophobic BLACK CARBON from biofuels [kg C/timestep]
         BIOF_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_YR

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read ORGANIC CARBON from biofuel combustion as tracer #35 in
      ! [kg C/year].  Convert to [kg C/timestep] and store in BIOF_BLKC.
      !=================================================================
      
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Hydrophilic ORGANIC CARBON from biofuels [kg C/timestep]
         BIOF_ORGC(I,J,1) =          FHO   * FD2D(I,J) / STEPS_PER_YR

         ! Hydrophobic ORGANIC CARBON from biofuels [kg C/timestep]
         BIOF_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_YR

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE ANTHRO_CARB_TBOND

!------------------------------------------------------------------------------

      SUBROUTINE ANTHRO_CARB_COOKE( THISMONTH )
!
!******************************************************************************
!  Subroutine ANTHRO_CARB_COOKE computes monthly mean anthropogenic and 
!  biofuel emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC 
!  CARBON.  It also separates these into HYDROPHILIC and HYDROPHOBIC 
!  fractions. (rjp, bmy, 4/2/04)
!
!  Emissions data comes from the Cooke et al. [1999] inventory and 
!  seasonality imposed by Park et al. [2003].  The data has units of 
!  [kg C/month].  This will be converted to [kg C/timestep] below.
!
!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

#     include "CMN_SIZE"
#     include "CMN_SETUP"  ! DATA_DIR

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, J
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU, STEPS_PER_MON
      REAL*8              :: FD2D(IIPAR,JJPAR)
      CHARACTER(LEN=255)  :: FILENAME

      ! Hydrophilic fraction of BLACK CARBON aerosol
      REAL*8, PARAMETER   :: FHB = 0.2d0

      ! Hydrophilic fraction of ORGANIC CARBON aerosol
      REAL*8, PARAMETER   :: FHO = 0.5d0

      !=================================================================
      ! ANTHRO_CARB_COOKE begins here!
      !=================================================================

      ! Number of emission timesteps per month
      STEPS_PER_MON = ( ( 1440 * NDAYS( THISMONTH ) ) / GET_TS_EMIS() )
      
      ! Get TAU0 value to index the punch file
      XTAU = GET_TAU0( THISMONTH, 1, 1998 )

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from 
      ! anthropogenic sources as tracer #34 in [kg C/month].  
      ! Then convert to [kg C/timestep] and store in ANTH_BLKC.
      !=================================================================

      ! Filename
      FILENAME = TRIM( DATA_DIR )                    //
     &           'carbon_200404/BCOC_anthsrce.geos.' // 
     &            GET_RES_EXT()
       
      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - ANTHRO_CARB_COOKE: Reading ', a )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Hydrophilic BLACK CARBON from anthropogenics [kg C/timestep]
         ANTH_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic BLACK CARBON from anthropogenics [kg C/timestep]
         ANTH_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_MON
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read ORGANIC CARBON from anthropogenic sources as tracer #35
      ! in [kg C/month].  Then Convert to [kg C/timestep] and store in 
      ! ANTH_ORGC.
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Hydrophilic ORGANIC CARBON from anthropogenics [kg C/timestep]
         ANTH_ORGC(I,J,1) = FHO * FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic ORGANIC CARBON from anthropogenics [kg C/timestep]
         ANTH_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_MON
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from biofuel 
      ! combustion over Canada and the US as tracer #34 in [kg C/year].  
      ! Then convert to [kg C/timestep] and store in BIOF_BLKC.
      !
      ! Seasonality has been imposed using the heating degree approach 
      ! for year 1998 [Park et al., 2003].
      !=================================================================

      ! Filename
      FILENAME = TRIM( DATA_DIR )                   //
     &           'carbon_200404/BCOC_biofuel.geos.' // 
     &           GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Hydrophilic BLACK CARBON from biofuels [kg C/timestep]
         BIOF_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic BLACK CARBON from biofuels [kg C/timestep]
         BIOF_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_MON

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read ORGANIC CARBON emission from biofuel combustion over 
      ! Canada and the US as tracer #35 in [kg C/year].  Then convert 
      ! to [kg C/timestep] and store in BIOF_ORGC.
      !
      ! Seasonality has been imposed using the heating degree approach 
      ! for year 1998 [Park et al., 2003].
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Hydrophilic ORGANIC CARBON from biofuels [kg C/timestep]
         BIOF_ORGC(I,J,1) =          FHO   * FD2D(I,J) / STEPS_PER_MON

         ! Hydrophobic ORGANIC CARBON from biofuels [kg C/timestep]
         BIOF_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_MON

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE ANTHRO_CARB_COOKE

!------------------------------------------------------------------------------

      SUBROUTINE BIOMASS_CARB_TBOND
!
!******************************************************************************
!  Subroutine BIOMASS_CARB_TBOND computes annual mean biomass burning 
!  emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC CARBON.  
!  It also separates these into HYDROPHILIC and HYDROPHOBIC fractions. 
!  (rjp, bmy, 4/2/04)
!
!  Emissions data comes from the Bond et al [2004] inventory and has units
!  of [kg C/yr].  This will be converted to [kg C/timestep] below.
!
!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

#     include "CMN_SIZE"   ! Size parameters 
#     include "CMN_SETUP"  ! DATA_DIR

      ! Local variables
      INTEGER             :: I, J
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU, STEPS_PER_YR     
      REAL*8              :: FD2D(IIPAR,JJPAR)
      CHARACTER(LEN=255)  :: FILENAME

      ! Hydrophilic fraction of carbonaceous aerosols
      REAL*8, PARAMETER   :: FHB = 0.2d0
      REAL*8, PARAMETER   :: FHO = 0.5d0

      !=================================================================
      ! BIOMASS_CARB_TBOND begins here!
      !=================================================================

      ! Number of emission timesteps per year
      STEPS_PER_YR = ( ( 1440 * 365 ) / GET_TS_EMIS() )

      ! Filename containing biomass emissions
      FILENAME = TRIM( DATA_DIR )                         //
     &           'carbon_200404/BCOC_TBond_biomass.geos.' // 
     &            GET_RES_EXT()

      ! Get TAU0 value to index the punch file
      XTAU = GET_TAU0( 1, 1, 2001 )

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - BIOMASS_CARB_TBOND: Reading ', a )

      !=================================================================
      ! Read BLACK CARBON (aka ELEMENTAL CARBON) emission from  
      ! biomass burning as tracer #34 in [kg C/year].  Then 
      ! convert to [kg C/timestep] and store in BIOB_BLKC.
      !=================================================================  
      CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 34, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Hydrophilic BLACK CARBON from biomass [kg C/timestep]
         BIOB_BLKC(I,J,1) =          FHB   * FD2D(I,J) / STEPS_PER_YR

         ! Hydrophobic BLACK CARBON from biomass [kg C/timestep]
         BIOB_BLKC(I,J,2) = ( 1.d0 - FHB ) * FD2D(I,J) / STEPS_PER_YR

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Read ORGANIC CARBON from biomass burning as tracer #35 in 
      ! [kg C/year].  Then convert to [kg C/timestep] and store in 
      ! BIOF_BLKC.
      !=================================================================  
      CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 35, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Hydrophilic ORGANIC CARBON from biomass [kg C/timestep]
         BIOB_ORGC(I,J,1) =          FHO   * FD2D(I,J) / STEPS_PER_YR

         ! Hydrophobic ORGANIC CARBON from biomass [kg C/timestep]
         BIOB_ORGC(I,J,2) = ( 1.d0 - FHO ) * FD2D(I,J) / STEPS_PER_YR

      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      ! Return to calling program
      END SUBROUTINE BIOMASS_CARB_TBOND

!------------------------------------------------------------------------------

      SUBROUTINE BIOMASS_CARB_GEOS( THISMONTH )
!
!******************************************************************************
!  Subroutine BIOMASS_CARB_TBOND computes annual mean biomass burning 
!  emissions of BLACK CARBON (aka ELEMENTAL CARBON) and ORGANIC CARBON.  
!  It also separates these into HYDROPHILIC and HYDROPHOBIC fractions. 
!  (rjp, bmy, 4/2/04)
!
!  Emissions data comes from the Bond et al [2004] inventory and has units
!  of [kg C/yr].  This will be converted to [kg C/timestep] below.
!
!  We also assume that 20% of BC and 50% of OC from anthropogenic 
!  emissions are hydrophilic (soluble) and the rest are hydrophobic.
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! NSRCE, DXYP
#     include "CMN_SETUP"    ! DATA_DIR

      !-------------------
      ! Arguments
      !-------------------
      INTEGER, INTENT(IN) :: THISMONTH

      !-------------------
      ! Local variables
      !-------------------

      ! Hydrophilic fraction of BLACK CARBON
      REAL*8, PARAMETER   :: FHB = 0.2d0

      ! Hydrophilic fraction of ORGANIC CARBON
      REAL*8, PARAMETER   :: FHO = 0.5d0 

      ! Black Carbon aerosols emission factor 
      ! from biomass burning (0.002 kgC/kg)
      REAL*8, PARAMETER   :: FEC = 2.d-3

      ! Organic Carbon aerosols emission factor 
      ! from biomass burning (0.014 kgC/kg)
      REAL*8, PARAMETER   :: FOC = 1.4d-2

      INTEGER             :: I, J
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: FD2D(IIPAR,JJPAR)
      REAL*8              :: XTAU, BIOCARB
      REAL*8              :: STEPS_PER_MON, AREA_M2
      CHARACTER(LEN=255)  :: FILENAME
      LOGICAL, SAVE       :: FIRSTEMIS = .TRUE.
      LOGICAL             :: EMFAC_VEG = .TRUE.

      !=================================================================
      ! BIOMASS_CARB_GEOS begins here!
      !=================================================================

      ! Number of emission timesteps per month
      STEPS_PER_MON = ( 1440 * NDAYS( THISMONTH ) ) / GET_TS_EMIS()

      ! Only do the following on the first timestep
      IF ( FIRSTEMIS ) THEN

         ! Read vegetation factors?
         IF ( EMFAC_VEG ) THEN

            !===========================================================
            ! If EMFAC_VEG=T, then read carbon aerosol emission factor 
            ! [kg/kg] from biomass burning sources depedning on surface 
            ! vegetation type compiled by [rjp, 2003].  Emission factors 
            ! are from Andreae and Merlet [2001].
            !===========================================================
            FILENAME = TRIM( DATA_DIR )               //
     &                'carbon_200404/emis_fac.EC-OC.' // GET_RES_EXT()

            ! TAU value for reading from the bpch files
            XTAU = GET_TAU0( 1, 1, 1985 )

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )
 100        FORMAT( '     - BIOMASS_CARB_GEOS: Reading ', a )

            !------------------
            ! BLACK CARBON
            !------------------
            CALL READ_BPCH2( FILENAME, 'EMISFAC', 81,  
     &                       XTAU,      IGLOB,    JGLOB,
     &                       1,         ARRAY,    QUIET=.TRUE. ) 

            ! Cast to REAL*8 and resize
            CALL TRANSFER_2D( ARRAY(:,:,1), EF_BLKC )

            !------------------
            ! ORGANIC CARBON
            !------------------
            CALL READ_BPCH2( FILENAME, 'EMISFAC', 82,  
     &                       XTAU,      IGLOB,    JGLOB,     
     &                       1,         ARRAY,    QUIET=.TRUE. ) 
            
            ! Cast to REAL*8 and resize
            CALL TRANSFER_2D( ARRAY(:,:,1), EF_ORGC )

         ELSE

            !===========================================================
            ! If EMFAC_VEG=F, then use these emission factors:
            !   BLACK CARBON   : 0.002 [kg C/kg]
            !   ORGANIC CARBON : 0.014 [kg C/kg] 
            ! Thus OC/BC = 7. [Chin et al., 2000]
            !============================================================
            EF_BLKC(:,:) = FEC
            EF_ORGC(:,:) = FOC
           
         ENDIF
        
         ! Reset first-time flag
         FIRSTEMIS = .FALSE.          
      ENDIF

      !=================================================================
      ! Read TOTAL biomass burning [g/cm2/month] as tracer #33
      ! Convert to [kg C/box/timestep] and store in BIOCARB.  
      ! 
      ! Then compute HYDROPHILIC and HYDROPHOBIC fractions of
      ! BLACK CARBON and ORGANIC CARBON.
      !=================================================================

      ! Filename
      FILENAME = TRIM( DATA_DIR )                        //
     &           'biomass_200110/bioburn.seasonal.geos.' //
     &           GET_RES_EXT()

      ! Get TAU value for reading the punch file
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 33, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D ( ARRAY(:,:,1), FD2D )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, AREA_M2, BIOCARB )
        DO J = 1, JJPAR

           ! Surface area [m2]
           AREA_M2 = GET_AREA_M2( J )

           DO I = 1, IIPAR

              ! Convert total biomass burned from [g/cm2/month] to
              ! [kg C/box/timestep].  The factor 10 comes from 
              ! 1e4 cm2/m2 divided by 1000 g/kg.  
              BIOCARB = FD2D(I,J) * AREA_M2 * 10.d0 / STEPS_PER_MON

              ! Hydrophilic BLACK CARBON from biomass [kg C/timestep]
              BIOB_BLKC(I,J,1) =          FHB   * EF_BLKC(I,J) * BIOCARB

              ! Hydrophobic BLACK CARBON from biomass [kg C/timestep]
              BIOB_BLKC(I,J,2) = ( 1.D0 - FHB ) * EF_BLKC(I,J) * BIOCARB

              ! Hydrophilic ORGANIC CARBON from biomass [kg C/timestep]
              BIOB_ORGC(I,J,1) =          FHO   * EF_ORGC(I,J) * BIOCARB

              ! Hydrophobic ORGANIC CARBON from biomass [kg C/timestep]
              BIOB_ORGC(I,J,2) = ( 1.D0 - FHO ) * EF_ORGC(I,J) * BIOCARB

           ENDDO
        ENDDO
!$OMP END PARALLEL DO  

        ! Return to calling program
        END SUBROUTINE BIOMASS_CARB_GEOS

!------------------------------------------------------------------------------

      SUBROUTINE EMITHIGH( BCSRC, OCSRC )
!
!******************************************************************************
!  Subroutine EMITHIGH mixes tracer completely from the surface to the PBL
!  top. (rjp, bmy, 4/2/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) BCSRC (REAL*8) : Array which holds Total BC (H-phobic & H-philic)
!  (2 ) OCSRC (REAL*8) : Array which holds Total OC (H-phobic & H-philic)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : PBL
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TRACERID_MOD, ONLY : IDTBCPI, IDTBCPO, IDTOCPI, IDTOCPO
      USE PRESSURE_MOD, ONLY : GET_PEDGE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! STT
#     include "CMN_GCTM"  ! SCALE_HEIGHT

      ! Arguments
      REAL*8, INTENT(IN) :: BCSRC(IIPAR,JJPAR,2)
      REAL*8, INTENT(IN) :: OCSRC(IIPAR,JJPAR,2)

      ! Local variables
      LOGICAL            :: IS_BCPO, IS_OCPO, IS_BCPI, IS_OCPI
      INTEGER            :: I,       J,       L
      REAL*8             :: BLTOP,   BLTHIK,  FTOT,    PB
      REAL*8             :: PT,      DELP,    FEMIS

      !=================================================================
      ! EMITHIGH begins here!
      !=================================================================

      ! Define logical flags for expediency
      IS_BCPI = ( IDTBCPI > 0 )
      IS_OCPI = ( IDTOCPI > 0 ) 
      IS_BCPO = ( IDTBCPO > 0 )
      IS_OCPO = ( IDTOCPO > 0 )

      !=================================================================
      ! Compute FEMIS -- fraction of box (I,J,L) w/in the PBL
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, FTOT, BLTOP, BLTHIK, L, PB, PT, DELP, FEMIS )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize
         FTOT  = 0d0
         FEMIS = 0d0

#if   defined( GEOS_4 )

         ! BLTOP = pressure at PBL top [hPa]
         ! Use barometric law since PBL is in [m]
         BLTOP  = GET_PEDGE(I,J,1) * EXP( -PBL(I,J) / SCALE_HEIGHT )

         ! BLTHIK is PBL thickness [hPa]
         BLTHIK = GET_PEDGE(I,J,1) - BLTOP

#else

         ! BLTOP = pressure of PBL top [hPa]
         BLTOP  = GET_PEDGE(I,J,1) - MAX( PBL(I,J), 1.D0 )

         ! BLTHIK is PBL thickness [hPa]
         BLTHIK = PBL(I,J)

#endif

         !==============================================================
         ! Loop thru tropospheric levels
         !==============================================================
         DO L = 1, LLTROP

            ! Pressure at edges of grid box(I,J,L) [hPa]
            PB   = GET_PEDGE(I,J,L)
            PT   = GET_PEDGE(I,J,L+1)

            ! Thickness of grid box (I,J,L) [hPa]
            DELP = PB - PT

            ! FEMIS is the fraction of the PBL 
            ! which is occupied by this level L
            IF ( BLTOP <= PT )  THEN
               FEMIS = DELP / BLTHIK

            ELSEIF ( BLTOP > PT .AND. BLTOP <= PB ) THEN
               FEMIS = ( PB - BLTOP ) / BLTHIK

            ELSEIF ( BLTOP > PB ) THEN
               CYCLE      
 
            ENDIF
            
            ! Fraction of data partitioned into 
            ! each level should sum to 1.0
            FTOT = FTOT + FEMIS

            !===========================================================
            ! Partition organic tracers equally throughout the 
            ! boundary layer -- store back into STT array
            !===========================================================

            ! Hydrophilic BLACK CARBON
            IF ( IS_BCPI ) THEN
               STT(I,J,L,IDTBCPI) = STT(I,J,L,IDTBCPI) + 
     &                              ( FEMIS * BCSRC(I,J,1) )
            ENDIF

            ! Hydrophilic ORGANIC CARBON
            IF ( IS_OCPI ) THEN
               STT(I,J,L,IDTOCPI) = STT(I,J,L,IDTOCPI) + 
     &                              ( FEMIS * OCSRC(I,J,1) )
            ENDIF
            
            ! Hydrophobic BLACK CARBON
            IF ( IS_BCPO ) THEN
               STT(I,J,L,IDTBCPO) = STT(I,J,L,IDTBCPO) + 
     &                              ( FEMIS * BCSRC(I,J,2) )
            ENDIF

            ! Hydrophobic ORGANIC CARBON
            IF ( IS_OCPO ) THEN
               STT(I,J,L,IDTOCPO) = STT(I,J,L,IDTOCPO) + 
     &                              ( FEMIS * OCSRC(I,J,2) )
            ENDIF

         ENDDO

         ! Error check
         IF ( ABS( FTOT - 1.d0 ) > 1.d-3 ) THEN
            CALL ERROR_STOP( 'Check vertical. distribution!',
     &                       'EMITHIGH ("carbon_mod.f")' )
         ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE EMITHIGH

!------------------------------------------------------------------------------

      SUBROUTINE INIT_CARBON
!
!******************************************************************************
!  Subroutine INIT_CARBON initializes all module arrays (rjp, bmy, 4/1/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE" ! Size parameters

      ! Local variables
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS

      !=================================================================
      ! INIT_CARBON begins here!
      !=================================================================
      
      ! Return if we already allocated arrays
      IF ( IS_INIT ) RETURN

      ALLOCATE( ANTH_BLKC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ANTH_BLKC' )
      ANTH_BLKC = 0d0

      ALLOCATE( ANTH_ORGC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ANTH_ORGC' )
      ANTH_ORGC = 0d0

      ALLOCATE( BIOB_BLKC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOB_BLKC' )
      BIOB_BLKC = 0d0

      ALLOCATE( BIOB_ORGC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOB_ORGC' )
      BIOB_ORGC = 0d0

      ALLOCATE( BIOF_BLKC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOF_BLKC' )
      BIOF_BLKC = 0d0

      ALLOCATE( BIOF_ORGC( IIPAR, JJPAR, 2 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOF_ORGC' )
      BIOF_ORGC = 0d0

      ALLOCATE( TERP_ORGC( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TERP_ORGC' )
      TERP_ORGC = 0d0

      ALLOCATE( BCCONV( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCCONV' )
      BCCONV = 0d0

      ALLOCATE( OCCONV( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCCONV' )
      OCCONV = 0d0

      ! These only have to be allocated if we are
      ! reading in monthly mean biomass burning
      IF ( USE_MONTHLY_BIOB ) THEN

         ALLOCATE( EF_BLKC( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EF_BLKC' )
         EF_BLKC = 0d0

         ALLOCATE( EF_ORGC( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EF_ORGC' )
         EF_ORGC = 0d0

      ENDIF

      ! Reset IS_INIT
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_CARBON

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_CARBON
!
!******************************************************************************
!  Subroutine CLEANUP_CARBON deallocates all module arrays (rjp, bmy, 4/1/04)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_CARBON begins here!
      !=================================================================
      IF ( ALLOCATED( ANTH_BLKC ) ) DEALLOCATE( ANTH_BLKC )
      IF ( ALLOCATED( ANTH_ORGC ) ) DEALLOCATE( ANTH_ORGC )
      IF ( ALLOCATED( BIOB_BLKC ) ) DEALLOCATE( BIOB_BLKC )
      IF ( ALLOCATED( BIOB_ORGC ) ) DEALLOCATE( BIOB_ORGC )
      IF ( ALLOCATED( BIOF_BLKC ) ) DEALLOCATE( BIOF_BLKC )
      IF ( ALLOCATED( BIOF_ORGC ) ) DEALLOCATE( BIOF_ORGC )
      IF ( ALLOCATED( TERP_ORGC ) ) DEALLOCATE( TERP_ORGC )
      IF ( ALLOCATED( BCCONV    ) ) DEALLOCATE( BCCONV    )
      IF ( ALLOCATED( OCCONV    ) ) DEALLOCATE( OCCONV    )
      IF ( ALLOCATED( EF_BLKC   ) ) DEALLOCATE( EF_BLKC   )
      IF ( ALLOCATED( EF_ORGC   ) ) DEALLOCATE( EF_ORGC   )

      ! Return to calling program
      END SUBROUTINE CLEANUP_CARBON

!------------------------------------------------------------------------------

      END MODULE CARBON_MOD
