
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: pops_mod.f
!
! !DESCRIPTION: This module contains variables and routines for the 
!  GEOS-Chem peristent organic pollutants (POPs) simulation. 
!\\
!\\
! !INTERFACE: 
!
      MODULE POPS_MOD
! 
! !USES:
!
      IMPLICIT NONE
! Make everything Private ...
      PRIVATE
!
! !PUBLIC TYPES:
!
      PUBLIC :: EMISSPOPS
      PUBLIC :: POPS_READ_FOC

! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CHEMPOPS
      PUBLIC :: INIT_POPS
!      PUBLIC :: GET_POP_TYPE
!      PUBLIC :: GET_EMISSFILE
!      PUBLIC :: GET_POP_XMW
!      PUBLIC :: GET_POP_HSTAR
!      PUBLIC :: GET_POP_DEL_Hw
!      PUBLIC :: GET_POP_KOA
!
! !PUBLIC DATA MEMBERS:
!

! !REVISION HISTORY:
!  20 September 2010 N.E. Selin - Initial Version
!
! !REMARKS:
! Under construction
!
!EOP

!******************************************************************************
!Comment header:
!
!  Module Variables:
!  ===========================================================================
!  (1 ) TCOSZ    (REAL*8) : Sum of COS(Solar Zenith Angle ) [unitless]
!  (2 ) TTDAY    (REAL*8) : Total daylight time at location (I,J) [minutes]
!  (3 ) ZERO_DVEL(REAL*8) : Array with zero dry deposition velocity [cm/s]
!  (4 ) COSZM    (REAL*8) : Max daily value of COS(S.Z. angle) [unitless]  
!
!  Module Routines:
!  ===========================================================================
!  (1 ) CHEMPOPS
!  (2 ) INIT_POPS
!  (3 ) CHEM_POPGP
!  (4 ) EMISSPOPS
!  (5 ) EMITPOP
!  (6 ) OHNO3TIME
!  (7 ) CLEANUP_POPS
!
!  Module Functions:
!  ===========================================================================
!
!
!  GEOS-CHEM modules referenced by pops_mod.f
!  ===========================================================================
!
!
!  Nomenclature: 
!  ============================================================================
!
!
!  POPs Tracers
!  ============================================================================
!  (1 ) POPG               : Gaseous POP - total tracer  
!  (2 ) POPPOCPO           : Hydrophobic OC-sorbed POP  - total tracer
!  (3 ) POPPBCPO           : Hydrophobic BC-sorbed POP  - total tracer
!  (4 ) POPPOCPI           : Hydrophilic OC-sorbed POP  - total tracer
!  (5 ) POPPBCPI           : Hydrophilic BC-sorbed POP  - total tracer
!
!
!  References:
!  ============================================================================
!
!
!  Notes:
!  ============================================================================
!  (1) 20 September 2010 N.E. Selin - Initial version
!  (2) 4 January 2011 C.L. Friedman - Expansion on initial version
!
!
!******************************************************************************
!
      ! References to F90 modules

      

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
      REAL*8,  PARAMETER   :: SMALLNUM = 1D-20
      ! Arrays
      REAL*8,  ALLOCATABLE :: TCOSZ(:,:)
      REAL*8,  ALLOCATABLE :: TTDAY(:,:)
      REAL*8,  ALLOCATABLE :: ZERO_DVEL(:,:)
      REAL*8,  ALLOCATABLE :: COSZM(:,:)
      REAL*8,  ALLOCATABLE :: EPOP_G(:,:,:)
      REAL*8,  ALLOCATABLE :: EPOP_OC(:,:,:)
      REAL*8,  ALLOCATABLE :: EPOP_BC(:,:,:)
      REAL*8,  ALLOCATABLE :: EPOP_P_TOT(:,:,:)
      REAL*8,  ALLOCATABLE :: POP_TOT_EM(:,:)
      REAL*8,  ALLOCATABLE :: POP_SURF(:,:)
      REAL*8,  ALLOCATABLE :: EPOP_VEG(:,:)
      REAL*8,  ALLOCATABLE :: EPOP_LAKE(:,:)
      REAL*8,  ALLOCATABLE :: EPOP_SOIL(:,:)
      REAL*8,  ALLOCATABLE :: EPOP_OCEAN(:,:)
      REAL*8,  ALLOCATABLE :: EPOP_SNOW(:,:)
      REAL*8,  ALLOCATABLE :: F_OC_SOIL(:,:)
      REAL*8,  ALLOCATABLE :: C_OC(:,:,:),        C_BC(:,:,:)
      REAL*8,  ALLOCATABLE :: SUM_OC_EM(:,:), SUM_BC_EM(:,:)
      REAL*8,  ALLOCATABLE :: SUM_G_EM(:,:), SUM_OF_ALL(:,:)


      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CHEMPOPS
!
! !DESCRIPTION: This routine is the driver routine for POPs chemistry 
!  (eck, 9/20/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEMPOPS
!
! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
! 
      ! References to F90 modules
      USE DRYDEP_MOD,    ONLY : DEPSAV
      USE ERROR_MOD,     ONLY : DEBUG_MSG
      USE GLOBAL_OH_MOD, ONLY : GET_GLOBAL_OH
      USE GLOBAL_O3_MOD, ONLY : GET_GLOBAL_O3
      USE GLOBAL_NO3_MOD, ONLY: GET_GLOBAL_NO3
      USE GLOBAL_OC_MOD, ONLY : GET_GLOBAL_OC !clf, 1/20/2011
      USE GLOBAL_BC_MOD, ONLY : GET_GLOBAL_BC !clf, 1/20/2011
      USE PBL_MIX_MOD,   ONLY : GET_PBL_MAX_L
      USE LOGICAL_MOD,   ONLY : LPRT, LGTMM, LNLPBL !CDH added LNLPBL
      USE TIME_MOD,      ONLY : GET_MONTH, ITS_A_NEW_MONTH, GET_YEAR
      USE TRACER_MOD,    ONLY : N_TRACERS
      USE DRYDEP_MOD,    ONLY : DRYPOPG, DRYPOPP_OCPO, DRYPOPP_BCPO
      USE DRYDEP_MOD,    ONLY : DRYPOPP_OCPI, DRYPOPP_BCPI  !clf, 2/13/2012

      USE CMN_SIZE_MOD             ! Size parameters ! hma, 07Mar2013

!#     include "CMN_SIZE"      ! Size parameters

!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on CHEMMERCURY from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine CHEMPOPS is the driver routine for POPs chemistry
!  in the GEOS-CHEM module. (eck, clf, 1/4/2011)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) 
!
!
!  Local variables:
!  ============================================================================
!  (1 )


!  NOTES:
!  (1 )
!******************************************************************************

!BOC

      ! Local variables
      ! LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I, J, L, MONTH, YEAR, N, PBL_MAX

      !=================================================================
      ! CHEMPOPS begins here!
      !
      ! Read monthly mean OH fields for oxidation and monthly OC and BC
      ! fields for gas-particle partioning
      !=================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN 

         ! Get the current month
         MONTH = GET_MONTH()

         ! Get the current year
         YEAR = GET_YEAR()

         ! Read monthly mean oxidants from disk
         
         CALL GET_GLOBAL_OH( MONTH)
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMPOPS: a GET_GLOBAL_OH' )

         ! Read monthly mean O3 from disk
         CALL GET_GLOBAL_O3( MONTH )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMPOPS: a GET_GLOBAL_O3' )

         ! Read monthly mean NO3 from disk 
         CALL GET_GLOBAL_NO3( MONTH )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMPOPS: a GET_GLOBAL_NO3' )

         ! Read monthly OC from disk
         ! Do this in EMISSPOPS now
         !CALL GET_GLOBAL_OC( MONTH )
         !IF ( LPRT ) CALL DEBUG_MSG( '### CHEMPOPS: a GET_GLOBAL_OC' )

         ! Read monthly BC from disk
         ! Do this in EMISSPOPS now
         !CALL GET_GLOBAL_BC( MONTH )
         !IF ( LPRT ) CALL DEBUG_MSG( '### CHEMPOPS: a GET_GLOBAL_BC' ) 

      ENDIF
     
      ! If it's a new 6-hr mean, then get the current average 3-D temperature    

      !=================================================================
      ! Perform chemistry on POPs tracers
      !=================================================================
      
      ! Compute diurnal scaling for OH and NO3
      CALL OHNO3TIME
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMPOPS: a OHNO3TIME' )

      !-------------------------
      ! GAS AND PARTICLE PHASE chemistry
      !-------------------------
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMPOPS: b CHEM_GASPART' )
      
      ! Add option for non-local PBL (cdh, 08/27/09)
      IF ( LNLPBL ) THEN

         ! Dry deposition occurs with PBL mixing,
         ! pass zero deposition frequency
         CALL CHEM_POPGP( ZERO_DVEL, ZERO_DVEL, ZERO_DVEL, ZERO_DVEL,
     &      ZERO_DVEL)
         
      ELSE

         ! For addition of hydrophilic OC/BC POP tracers: for now, assume that 
         ! if dry deposition of hydrophoblic tracers is active, then dry deposition
         ! of hydrophilic tracers is also active (ie., they're coupled). clf, 2/12/2012
         IF ( DRYPOPG > 0 .and. DRYPOPP_OCPO > 0 .and.
     &      DRYPOPP_BCPO > 0) THEN
         
            ! Dry deposition active for both POP-Gas and POP-Particle; 
            ! pass drydep frequency to CHEM_POPGP (NOTE: DEPSAV has units 1/s)
            CALL CHEM_POPGP(DEPSAV(:,:,DRYPOPG), 
     &          DEPSAV(:,:,DRYPOPP_OCPO), DEPSAV(:,:,DRYPOPP_OCPI), 
     &          DEPSAV(:,:,DRYPOPP_BCPO), DEPSAV(:,:,DRYPOPP_BCPI)) 

           ELSEIF (DRYPOPG > 0 .and. DRYPOPP_OCPO > 0 .and. 
     &          DRYPOPP_BCPO .le. 0 ) THEN

            ! Only POPG and POPP_OC dry deposition are active
            CALL CHEM_POPGP(DEPSAV(:,:,DRYPOPG), 
     &           DEPSAV(:,:,DRYPOPP_OCPO), DEPSAV(:,:,DRYPOPP_OCPI),
     &           ZERO_DVEL, ZERO_DVEL) 

           ELSEIF (DRYPOPG > 0 .and. DRYPOPP_OCPO .le. 0 .and. 
     &          DRYPOPP_BCPO > 0 ) THEN

            ! Only POPG and POPP_BC dry deposition are active
            CALL CHEM_POPGP(DEPSAV(:,:,DRYPOPG), ZERO_DVEL, ZERO_DVEL, 
     &          DEPSAV(:,:,DRYPOPP_BCPO), DEPSAV(:,:,DRYPOPP_BCPI)) 
         
           ELSEIF (DRYPOPG > 0 .and. DRYPOPP_OCPO .le. 0 .and. 
     &          DRYPOPP_BCPO .le. 0 ) THEN

            ! Only POPG dry deposition is active
            CALL CHEM_POPGP( DEPSAV(:,:,DRYPOPG), ZERO_DVEL, ZERO_DVEL,
     &          ZERO_DVEL, ZERO_DVEL) 
            
           ELSEIF (DRYPOPG <= 0 .and. DRYPOPP_OCPO > 0 .and. 
     &          DRYPOPP_BCPO > 0) THEN

            ! Only POPP dry deposition is active
            CALL CHEM_POPGP( ZERO_DVEL , DEPSAV(:,:,DRYPOPP_OCPO), 
     &           DEPSAV(:,:,DRYPOPP_OCPI), DEPSAV(:,:,DRYPOPP_BCPO),
     &           DEPSAV(:,:,DRYPOPP_BCPI))

           ELSEIF (DRYPOPG <= 0 .and. DRYPOPP_OCPO > 0 .and. 
     &          DRYPOPP_BCPO <= 0) THEN

            ! Only POPP_OC dry deposition is active
            CALL CHEM_POPGP( ZERO_DVEL , DEPSAV(:,:,DRYPOPP_OCPO), 
     &           DEPSAV(:,:,DRYPOPP_OCPI), ZERO_DVEL, ZERO_DVEL)

           ELSEIF (DRYPOPG <= 0 .and. DRYPOPP_OCPO <= 0 .and. 
     &          DRYPOPP_BCPO > 0) THEN

            ! Only POPP_OC dry deposition is active
            CALL CHEM_POPGP( ZERO_DVEL , ZERO_DVEL, ZERO_DVEL,
     &           DEPSAV(:,:,DRYPOPP_BCPO), DEPSAV(:,:,DRYPOPP_BCPI))            
         ELSE

            ! No dry deposition, pass zero deposition frequency
            CALL CHEM_POPGP( ZERO_DVEL, ZERO_DVEL, ZERO_DVEL,
     &            ZERO_DVEL, ZERO_DVEL)

         ENDIF

      ENDIF      

      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMPOPS: a CHEM_GASPART' )
   
    
      ! Return to calling program
      END SUBROUTINE CHEMPOPS

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CHEM_POPGP
!
! !DESCRIPTION: This routine does chemistry for POPs gas and particles
!  (eck, 9/20/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_POPGP (V_DEP_G, V_DEP_P_OCPO, V_DEP_P_OCPI,
     &                       V_DEP_P_BCPO, V_DEP_P_BCPI)

!     References to F90 modules
      USE TRACER_MOD,   ONLY : STT,        XNUMOL
      USE TRACERID_MOD, ONLY : IDTPOPG,    IDTPOPPOCPO,  IDTPOPPBCPO
      USE TRACERID_MOD, ONLY : IDTPOPPOCPI,  IDTPOPPBCPI
      USE DIAG53_MOD,   ONLY : AD53_PG_OC_NEG, AD53_PG_BC_NEG
      USE DIAG53_MOD,   ONLY : AD53_PG_OC_POS, AD53_PG_BC_POS
      USE DIAG53_MOD,   ONLY : AD53_POPP_OCPO_O3, AD53_POPP_OCPI_O3
      USE DIAG53_MOD,   ONLY : AD53_POPP_BCPO_O3, AD53_POPP_BCPI_O3
      USE DIAG53_MOD,   ONLY : AD53_POPP_OCPO_NO3, AD53_POPP_OCPI_NO3
      USE DIAG53_MOD,   ONLY : AD53_POPP_BCPO_NO3, AD53_POPP_BCPI_NO3
      !USE DIAG53_MOD,   ONLY : ND53,       LD53,       AD53_POPG_OH
      USE DIAG53_MOD,   ONLY :  AD53_POPG_OH !hma, 8mar2013
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE DIAG_MOD,     ONLY : AD44
      USE LOGICAL_MOD,  ONLY : LNLPBL,     LGTMM
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE DAO_MOD,      ONLY : T,          AIRVOL
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE GET_POPSINFO_MOD, ONLY: GET_POP_DEL_H, GET_POP_KOA
      USE GET_POPSINFO_MOD, ONLY: GET_POP_KBC, GET_POP_K_POPG_OH
      USE GET_POPSINFO_MOD, ONLY: GET_POP_K_POPP_O3A
      USE GET_POPSINFO_MOD, ONLY: GET_POP_K_POPP_O3B

      ! H. Amos, 7Mar2013
      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND44
!#     include "CMN_SIZE" ! Size parameters
!#     include "CMN_DIAG" ! ND44

!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN)    :: V_DEP_G(IIPAR,JJPAR)
      REAL*8, INTENT(IN)    :: V_DEP_P_OCPO(IIPAR,JJPAR)
      REAL*8, INTENT(IN)    :: V_DEP_P_BCPO(IIPAR,JJPAR)
      REAL*8, INTENT(IN)    :: V_DEP_P_OCPI(IIPAR,JJPAR)
      REAL*8, INTENT(IN)    :: V_DEP_P_BCPI(IIPAR,JJPAR)

!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
    
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on CHEM_HG0_HG2 from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine CHEM_POPGP is the chemistry subroutine for the oxidation,
!  gas-particle partitioning, and deposition of POPs.
!  (eck, clf, 1/4/2011)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) V_DEP_G (REAL*8)      : Dry deposition frequency for gaseous POP [/s]
!  (2 ) V_DEP_P_OCPO (REAL*8) : Dry deposition frequency for OCPO-POP [/s]
!  (3 ) V_DEP_P_BCPO (REAL*8) : Dry deposition frequency for BCPO-POP [/s]
!  (4 ) V_DEP_P_OCPI (REAL*8) : Dry deposition frequency for OCPI-POP [/s]
!  (5 ) V_DEP_P_BCPI (REAL*8) : Dry deposition frequency for BC-PIPOP [/s]
!
!  Local variables:
!  ============================================================================
!  (1 )
!
!
!     
!  NOTES:
!  (1 ) 
!  
!  REFS:
!  (1 ) For OH rate constant: Brubaker & Hites. 1998. OH reaction kinetics of
!  PAHs and PCDD/Fs. J. Phys. Chem. A. 102:915-921. 
!
!******************************************************************************
!BOC
      ! Local variables
      INTEGER               :: I, J, L
      REAL*8                :: DTCHEM,       SUM_F
      REAL*8                :: KOA_T,        KBC_T
      REAL*8                :: KOC_BC_T,     KBC_OC_T
      REAL*8                :: TK
      REAL*8                :: AREA_CM2
      REAL*8                :: F_PBL
      REAL*8                :: C_OH,         C_OC_CHEM,   C_BC_CHEM
      REAL*8                :: C_OC_CHEM1,   C_BC_CHEM1
      REAL*8                :: K_OH,         AIR_VOL
      REAL*8                :: K_OX,         C_O3
      REAL*8                :: E_KOX_T
      REAL*8                :: K_DEPG,       K_DEPP_OCPO,   K_DEPP_BCPO
      REAL*8                :: K_DEPP_OCPI,  K_DEPP_BCPI
      REAL*8                :: OLD_POPG,    OLD_POPP_OCPO, OLD_POPP_BCPO
      REAL*8                :: OLD_POPP_OCPI, OLD_POPP_BCPI
      REAL*8                :: OLD_POPP_OC,   OLD_POPP_BC
      REAL*8                :: NEW_POPG,    NEW_POPP_OCPO, NEW_POPP_BCPO
      REAL*8                :: NEW_POPP_OCPI, NEW_POPP_BCPI
      REAL*8                :: POPG_BL,     POPP_OCPO_BL,  POPP_BCPO_BL
      REAL*8                :: POPP_OCPI_BL, POPP_BCPI_BL
      REAL*8                :: POPG_FT,     POPP_OCPO_FT,  POPP_BCPO_FT
      REAL*8                :: POPP_OCPI_FT, POPP_BCPI_FT
      REAL*8                :: TMP_POPG,      TMP_OX
      REAL*8                :: GROSS_OX,      GROSS_OX_OH, NET_OX
      REAL*8                :: DEP_POPG,    DEP_POPP_OCPO, DEP_POPP_BCPO
      REAL*8                :: DEP_POPP_OCPI, DEP_POPP_BCPI
      REAL*8                :: DEP_POPG_DRY,  DEP_POPP_OCPO_DRY
      REAL*8                :: DEP_POPP_BCPO_DRY, DEP_POPP_OCPI_DRY
      REAL*8                :: DEP_POPP_BCPI_DRY
      REAL*8                :: DEP_DRY_FLXG,  DEP_DRY_FLXP_OCPO
      REAL*8                :: DEP_DRY_FLXP_BCPO, DEP_DRY_FLXP_OCPI
      REAL*8                :: DEP_DRY_FLXP_BCPI
      REAL*8                :: OLD_POP_T
      REAL*8                :: VR_OC_AIR,     VR_BC_AIR
      REAL*8                :: VR_OC_BC,      VR_BC_OC
      REAL*8                :: F_POP_OC,      F_POP_BC
      REAL*8                :: F_POP_G
      REAL*8                :: MPOP_OCPO,     MPOP_BCPO,     MPOP_G
      REAL*8                :: MPOP_OCPI,     MPOP_BCPI
      REAL*8                :: MPOP_OC, MPOP_BC
      REAL*8                :: DIFF_G,        DIFF_OC,     DIFF_BC
      REAL*8                :: OC_AIR_RATIO,  OC_BC_RATIO, BC_AIR_RATIO
      REAL*8                :: BC_OC_RATIO,   SUM_DIFF
      REAL*8                :: FOLD, KOCBC, NEW_OCPI, NEW_BCPI
      REAL*8                :: GROSS_OX_OCPO, GROSS_OX_OCPI
      REAL*8                :: GROSS_OX_BCPO, GROSS_OX_BCPI
      REAL*8                :: GROSS_OX_O3_OCPO, GROSS_OX_O3_OCPI
      REAL*8                :: GROSS_OX_O3_BCPO, GROSS_OX_O3_BCPI
      REAL*8                :: GROSS_OX_NO3_OCPO, GROSS_OX_NO3_OCPI
      REAL*8                :: GROSS_OX_NO3_BCPO, GROSS_OX_NO3_BCPI
      REAL*8                :: TMP_POPP_OCPO, TMP_POPP_OCPI
      REAL*8                :: E_KOX_T_BC
      REAL*8                :: TMP_OX_P_OCPO, TMP_OX_P_OCPI
      REAL*8                :: TMP_OX_P_BCPO, TMP_OX_P_BCPI
      REAL*8                :: TMP_POPP_BCPO, TMP_POPP_BCPI
      REAL*8                :: NET_OX_OCPO, NET_OX_OCPI
      REAL*8                :: NET_OX_BCPO, NET_OX_BCPI
      REAL*8                :: K_O3_BC, C_NO3, K_OX_P, K_NO3_BC

      ! Delta H for POP [kJ/mol]. Delta H is enthalpy of phase transfer
      ! from gas phase to OC. For now we use Delta H for phase transfer 
      ! from the gas phase to the pure liquid state. 
      ! For PHENANTHRENE: 
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -74000 [J/mol].
      ! For PYRENE:
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -87000 [J/mol].    
      ! For BENZO[a]PYRENE:
      ! this is also taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 452, Prob 11.1), or -110,000 [J/mol]
      REAL*8     :: DEL_H

      ! R = universal gas constant for adjusting KOA for temp: 8.3145 [J/mol/K]
      REAL*8, PARAMETER     :: R          = 8.31d0  

      ! KOA_298 for partitioning of gas phase POP to atmospheric OC
      ! KOA_298 = Cpop in octanol/Cpop in atmosphere at 298 K 
      ! For PHENANTHRENE:
      ! log KOA_298 = 7.64, or 4.37*10^7 [unitless]
      ! For PYRENE:
      ! log KOA_298 = 8.86, or 7.24*10^8 [unitless]
      ! For BENZO[a]PYRENE:
      ! log KOA_298 = 11.48, or 3.02*10^11 [unitless]
      ! (Ma et al., J. Chem. Eng. Data, 2010, 55:819-825).
      REAL*8     :: KOA_298

      ! KBC_298 for partitioning of gas phase POP to atmospheric BC
      ! KBC_298 = Cpop in black carbon/Cpop in atmosphere at 298 K
      ! For PHENANTHRENE:
      ! log KBC_298 = 10.0, or 1.0*10^10 [unitless]
      ! For PYRENE:
      ! log KBC_298 = 11.0, or 1.0*10^11 [unitless]
      ! For BENZO[a]PYRENE:
      ! log KBC_298 = 13.9, or 7.94*10^13 [unitless]
      ! (Lohmann and Lammel, EST, 2004, 38:3793-3802)
      REAL*8     :: KBC_298

      ! DENS_OCT = density of octanol, needed for partitioning into OC
      ! 820 [kg/m^3]
      REAL*8, PARAMETER     :: DENS_OCT   = 82d1

      ! DENS_BC = density of BC, needed for partitioning onto BC
      ! 1 [kg/L] or 1000 [kg/m^3] 
      ! From Lohmann and Lammel, Environ. Sci. Technol., 2004, 38:3793-3803.
      REAL*8, PARAMETER     :: DENS_BC    = 1d3

      ! K for reaction POPG + OH  [cm3 /molecule /s]
      ! For PHENANTHRENE: 2.70d-11
      ! (Source: Brubaker & Hites, J. Phys Chem A 1998)
      ! For PYRENE: 5.00d-11
      ! Calculated by finding the ratio between kOH of phenanthrene and kOH of pyrene
      ! using structure-activity relationships (Schwarzenback, Gschwend, Imboden, 
      ! pg 680) and scaling the experimental kOH for phenanthrene from Brubaker and Hites
      ! For BENZO[a]PYRENE: 5.68d-11
      ! Calculated by finding the ratio between kOH of phenanthrene and kOH of pyrene
      ! using structure-activity relationships (Schwarzenback, Gschwend, Imboden, 
      ! pg 680) and scaling the experimental kOH for phenanthrene from Brubaker and Hites
      ! Could potentially set this to change with temmperature 

      ! Using EPA AOPWIN values:
      ! For PHENANTHRENE: 13d-12
      ! For PYRENE: 50d-12
      ! For BaP: 50d-12
      REAL*8     :: K_POPG_OH !(Gas phase)

      ! k for reaction POPP + O3 [/s] depends on fitting parameters A and B. 
      ! A represents the maximum number of surface sites available to O3, and B 
      ! represents the ratio of desorption/adsorption rate coefficients for both bulk
      ! phases (Ref: Kahan et al Atm Env 2006, 40:3448)
      ! k(obs) = A x [O3(g)] / (B + [O3(g)])
      ! For PHENANTHRENE: A = 0.5 x 10^-3 s^-1, B = 2.15 x 10^15 molec/cm3
      ! For PYRENE: A = 0.7 x 10^-3 s^-1, B = 3 x 10^15 molec/cm3
      ! for BaP: A = 5.5 x 10^-3 s^-1, B = 2.8 x 10^15 molec/cm3
      REAL*8     :: AK  ! s^-1
      REAL*8     :: BK  ! molec/cm3

      ! k for reaction POPP + NO3 taken from Liu et al. EST 2012 "Kinetic studies of heterogeneous
      ! reactions of PAH aerosols with NO3 radicals", for now. CLF, 1/24/2012
      ! For PYR, 6.4 x 10^-12 [cm3 / molec / s]
      REAL*8, PARAMETER     :: K_POPP_NO3 = 0d0 ! 6.4d-12

      ! OC/BC hydrophobic POP lifetime before folding to hydrophilic
      REAL*8, PARAMETER     :: OCBC_LIFE = 1.15D0

      REAL*8                :: DUM

      !=================================================================
      ! CHEM_POPGP begins here!
      !=================================================================
      DUM = 1.0
      DEL_H = GET_POP_DEL_H(DUM)
      KOA_298 = GET_POP_KOA(DUM)
      KBC_298 = GET_POP_KBC(DUM)
      K_POPG_OH = GET_POP_K_POPG_OH(DUM)
      AK = GET_POP_K_POPP_O3A(DUM)
      BK = GET_POP_K_POPP_O3B(DUM)

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Zero concentrations in loop
         MPOP_G = 0d0
         MPOP_OCPO = 0d0
         MPOP_BCPO = 0d0
         MPOP_OCPI = 0d0
         MPOP_BCPI = 0d0
         OLD_POPG = 0d0
         OLD_POPP_OCPO = 0d0
         OLD_POPP_BCPO = 0d0
         OLD_POPP_OCPI = 0d0
         OLD_POPP_BCPI = 0d0
         OLD_POP_T = 0d0
         NEW_POPG = 0d0
         NEW_POPP_OCPO = 0d0
         NEW_POPP_BCPO = 0d0
         NEW_POPP_OCPI = 0d0
         NEW_POPP_BCPI = 0d0
         POPG_BL = 0d0
         POPP_OCPO_BL = 0d0
         POPP_BCPO_BL = 0d0
         POPP_OCPI_BL = 0d0
         POPP_BCPI_BL = 0d0
         POPG_FT = 0d0
         POPP_OCPO_FT = 0d0
         POPP_BCPO_FT = 0d0
         POPP_OCPI_FT = 0d0
         POPP_BCPI_FT = 0d0
         DIFF_G = 0d0
         DIFF_OC = 0d0
         DIFF_BC = 0d0
         NET_OX = 0d0 
         TMP_POPG = 0d0
         TMP_OX = 0d0      
         GROSS_OX = 0d0
         GROSS_OX_OH = 0d0 
         DEP_POPG = 0d0
         DEP_POPP_OCPO = 0d0
         DEP_POPP_BCPO = 0d0
         DEP_POPP_OCPI = 0d0
         DEP_POPP_BCPI = 0d0
         DEP_POPG_DRY = 0d0
         DEP_POPP_OCPO_DRY = 0d0
         DEP_POPP_BCPO_DRY = 0d0
         DEP_POPP_OCPI_DRY = 0d0
         DEP_POPP_BCPI_DRY = 0d0
         DEP_DRY_FLXG = 0d0
         DEP_DRY_FLXP_OCPO = 0d0
         DEP_DRY_FLXP_BCPO = 0d0
         DEP_DRY_FLXP_OCPI = 0d0
         DEP_DRY_FLXP_BCPI = 0d0
         E_KOX_T = 0d0
         K_OX = 0d0
!         K_O3_OC = 0d0
         K_O3_BC = 0d0
         GROSS_OX_OCPO = 0d0
         GROSS_OX_BCPO = 0d0
         GROSS_OX_OCPI = 0d0
         GROSS_OX_BCPI = 0d0
         GROSS_OX_O3_OCPO = 0d0
         GROSS_OX_O3_BCPO = 0d0
         GROSS_OX_O3_OCPI = 0d0
         GROSS_OX_O3_BCPI = 0d0
         GROSS_OX_NO3_OCPO = 0d0
         GROSS_OX_NO3_BCPO = 0d0
         GROSS_OX_NO3_OCPI = 0d0
         GROSS_OX_NO3_BCPI = 0d0
         TMP_POPP_OCPO = 0d0 
         TMP_POPP_OCPI = 0d0
         TMP_POPP_BCPO = 0d0 
         TMP_POPP_BCPI = 0d0
         E_KOX_T_BC    = 0d0
         TMP_OX_P_OCPO = 0d0
         TMP_OX_P_OCPI = 0d0
         TMP_OX_P_BCPO = 0d0
         TMP_OX_P_BCPI = 0d0
         NET_OX_OCPO   = 0d0
         NET_OX_OCPI   = 0d0
         NET_OX_BCPO   = 0d0
         NET_OX_BCPI   = 0d0

         ! Save local temperature in TK for convenience [K]
         TK = T(I,J,L)

         ! Get monthly mean OH concentrations (molec/cm3)
         C_OH        = GET_OH( I, J, L )

         ! Get monthly mean O3 concentrations (molec/cm3)
         C_O3        = GET_O3( I, J, L )

         ! Get monthly mean NO3 concentration (molec/cm3)
! ----> this is causing a segfault because it's looking for 
!       L = 39 and NO3 is (72, 46, 38). (H. Amos, 10April2013)
         C_NO3       = GET_NO3( I, J, L )

         ! Fraction of box (I,J,L) underneath the PBL top [dimensionless]
         F_PBL = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Define K for the oxidation reaction with POPG [/s]
         K_OH        = K_POPG_OH * C_OH

         ! Define K for the oxidation reaction with POPPOC and POPPBC [/s]
         ! Kahan:
         K_O3_BC        = (AK * C_O3) / (BK + C_O3)
!         WRITE (6,*) 'K_O3_BC=', K_O3_BC
         
         ! Define K for oxidation of POPP by NO3  [/s]
         K_NO3_BC    = K_POPP_NO3 * C_NO3 
!         WRITE (6,*) 'K_NO3_BC=', K_NO3_BC

         ! Total K for gas phase oxidation [/s]
         K_OX        = K_OH !+ ...

         ! Total K for particle phase oxidation [/s]
         K_OX_P      = K_O3_BC + K_NO3_BC

         ! Define Ks for dry deposition of gas phase POP [/s]
         K_DEPG = V_DEP_G(I,J)

         ! Define Ks for dry deposition of hydrophoblic particle phase POP [/s]
         K_DEPP_OCPO = V_DEP_P_OCPO(I,J)
         K_DEPP_BCPO = V_DEP_P_BCPO(I,J)

         ! Define Ks for dry deposition of hydrophilic particle phase POP [/s]
         K_DEPP_OCPI = V_DEP_P_OCPI(I,J)
         K_DEPP_BCPI = V_DEP_P_BCPI(I,J)

         ! Precompute exponential factors [dimensionless]
         ! For gas phase
         E_KOX_T  = EXP( -K_OX  * DTCHEM )
         ! For particle phase
         E_KOX_T_BC = EXP(-K_OX_P  * DTCHEM)


         !==============================================================
         ! GAS-PARTICLE PARTITIONING
         !==============================================================

         OLD_POPG = MAX( STT(I,J,L,IDTPOPG), SMALLNUM )  ![kg]
         OLD_POPP_OCPO = MAX( STT(I,J,L,IDTPOPPOCPO), SMALLNUM )  ![kg]
         OLD_POPP_BCPO = MAX( STT(I,J,L,IDTPOPPBCPO), SMALLNUM )  ![kg]
         OLD_POPP_OCPI = MAX( STT(I,J,L,IDTPOPPOCPI), SMALLNUM ) ![kg]
         OLD_POPP_BCPI = MAX( STT(I,J,L,IDTPOPPBCPI), SMALLNUM ) ![kg]

         ! Total POPs in box I,J,L 
         OLD_POP_T = OLD_POPG + OLD_POPP_OCPO + OLD_POPP_BCPO +
     &                  OLD_POPP_OCPI + OLD_POPP_BCPI

         ! Define temperature-dependent partition coefficients:
         KOA_T = KOA_298 * EXP((-DEL_H/R) * ((1d0/TK) - 
     &                  (1d0/298d0))) 

         ! Define KBC_T, the BC-air partition coeff at temp T [unitless]
         ! TURN OFF TEMPERATURE DEPENDENCY FOR SENSITIVITY ANALYSIS
         KBC_T = KBC_298 * EXP((-DEL_H/R) * ((1d0/TK) - 
     &                  (1d0/298d0)))

         ! Define KOC_BC_T, the theoretical OC-BC part coeff at temp T [unitless]
         KOC_BC_T = KOA_T / KBC_T

         ! Define KBC_OC_T, the theoretical BC_OC part coeff at temp T [unitless]
         KBC_OC_T = 1d0 / KOC_BC_T

         ! Get monthly mean OC and BC concentrations [kg/box] 
         C_OC_CHEM = GET_OC(I,J,L)
         C_BC_CHEM = GET_BC(I,J,L)

         ! Convert to units of volume per box [m^3 OC or BC/box]
         C_OC_CHEM1        = C_OC_CHEM / DENS_OCT
         C_BC_CHEM1        = C_BC_CHEM / DENS_BC

         ! Get AIRVOL
         AIR_VOL = AIRVOL(I,J,L)

         ! Define volume ratios:
         ! VR_OC_AIR = volume ratio of OC to air [unitless]     
         VR_OC_AIR   = C_OC_CHEM1 / AIR_VOL ! could be zero

         ! VR_OC_BC  = volume ratio of OC to BC [unitless]
         VR_OC_BC    = C_OC_CHEM1 / C_BC_CHEM1 ! could be zero or undefined

         ! VR_BC_AIR = volume ratio of BC to air [unitless]
         VR_BC_AIR   = VR_OC_AIR / VR_OC_BC ! could be zero or undefined

         ! VR_BC_OC  = volume ratio of BC to OC [unitless]
         VR_BC_OC    = 1d0 / VR_OC_BC ! could be zero or undefined

         ! Redefine fractions of total POPs in box (I,J,L) that are OC-phase, 
         ! BC-phase, and gas phase with new time step (should only change if 
         ! temp changes or OC/BC concentrations change) 
         OC_AIR_RATIO = 1d0 / (KOA_T * VR_OC_AIR) 
         OC_BC_RATIO = 1d0 / (KOC_BC_T * VR_OC_BC) 

         BC_AIR_RATIO = 1d0 / (KBC_T * VR_BC_AIR) 
         BC_OC_RATIO = 1d0 / (KBC_OC_T * VR_BC_OC)

         ! If there are zeros in OC or BC concentrations, make sure they
         ! don't cause problems with phase fractions
         IF ( C_OC_CHEM > SMALLNUM .and. C_BC_CHEM > SMALLNUM ) THEN
            F_POP_OC  = 1d0 / (1d0 + OC_AIR_RATIO + OC_BC_RATIO) 
            F_POP_BC  = 1d0 / (1d0 + BC_AIR_RATIO + BC_OC_RATIO)
         
           ELSE IF (C_OC_CHEM > SMALLNUM .and.
     &             C_BC_CHEM .le. SMALLNUM ) THEN
           F_POP_OC  = 1d0 / (1d0 + OC_AIR_RATIO)
           F_POP_BC  = SMALLNUM           

           ELSE IF ( C_OC_CHEM .le. SMALLNUM .and.
     &              C_BC_CHEM > SMALLNUM ) THEN
           F_POP_OC  = SMALLNUM
           F_POP_BC  = 1d0 / (1d0 + BC_AIR_RATIO)

           ELSE IF ( C_OC_CHEM .le. SMALLNUM .and. 
     &              C_BC_CHEM .le. SMALLNUM) THEN
           F_POP_OC = SMALLNUM
           F_POP_BC = SMALLNUM
        ENDIF

         ! Gas-phase:
         F_POP_G   = 1d0 - F_POP_OC - F_POP_BC

         ! Check that sum equals 1
         SUM_F = F_POP_OC + F_POP_BC + F_POP_G
!         WRITE (6,*) 'SUM_F =', SUM_F
         
         ! Calculate new masses of POP in each phase [kg]
         ! OC-phase:
         MPOP_OC    = F_POP_OC * OLD_POP_T

         ! BC-phase
         MPOP_BC     = F_POP_BC * OLD_POP_T

         ! Gas-phase
         MPOP_G     = F_POP_G  * OLD_POP_T

         ! Ensure new masses of POP in each phase are positive
         MPOP_OC = MAX(MPOP_OC, SMALLNUM)
         MPOP_BC = MAX(MPOP_BC, SMALLNUM)
         MPOP_G  = MAX(MPOP_G,  SMALLNUM)     

         ! Calculate differences in masses in each phase from previous time
         ! step for storage in ND53 diagnostic

            DIFF_G = MPOP_G - OLD_POPG
            DIFF_OC = MPOP_OC - OLD_POPP_OCPO - OLD_POPP_OCPI
            DIFF_BC = MPOP_BC - OLD_POPP_BCPO - OLD_POPP_BCPI

          ! Sum of differences should equal zero
            SUM_DIFF = MAX(DIFF_G + DIFF_OC + DIFF_BC, SMALLNUM)
!            WRITE (6,*) 'SUM_DIFF =', SUM_DIFF

            !==============================================================
            ! ND53 diagnostic: Differences in distribution of gas and
            ! particle phases between time steps [kg]
            !==============================================================

            IF ( ND53 > 0 .AND. L <= LD53 ) THEN ! LD53 is max level
               
               IF (DIFF_OC .lt. 0) THEN
               AD53_PG_OC_NEG(I,J,L) = AD53_PG_OC_NEG(I,J,L)  + 
     &              DIFF_OC

               ELSE IF (DIFF_OC .eq. 0 .or. DIFF_OC .gt. 0) THEN
               AD53_PG_OC_POS(I,J,L) = AD53_PG_OC_POS(I,J,L)  + 
     &              DIFF_OC

               ENDIF

               IF (DIFF_BC .lt. 0) THEN
               AD53_PG_BC_NEG(I,J,L) = AD53_PG_BC_NEG(I,J,L)  + 
     &              DIFF_BC
               
               ELSE IF (DIFF_BC .eq. 0 .or. DIFF_BC .gt. 0) THEN
               AD53_PG_BC_POS(I,J,L) = AD53_PG_BC_POS(I,J,L)  + 
     &              DIFF_BC

               ENDIF
            ENDIF

         !==============================================================
         ! HYDROPHOBIC PARTICULATE POPS DECAY TO HYDROPHILIC (clf, 2/12/2012)
         !==============================================================

         ! Define the lifetime and e-folding time to be the same as hydrophobic to 
         ! hydrophilic aerosols
         KOCBC = 1.D0 / (86400d0 * OCBC_LIFE)

         ! Send hydrophobic to hydrophilic
         FOLD = KOCBC * DTCHEM

         ! Amount of hydrophobic particulate POP left after folding to hydrophilic
         MPOP_OCPO = MPOP_OC * EXP( -FOLD)
         MPOP_BCPO = MPOP_BC * EXP( -FOLD)

         ! Hydrophilic particulate POP already existing
         MPOP_OCPI = MAX( STT(I,J,L,IDTPOPPOCPI), SMALLNUM )  ![kg]
         MPOP_BCPI = MAX( STT(I,J,L,IDTPOPPBCPI), SMALLNUM )  ![kg]

         ! Hydrophilic POP that used to be hydrophobic
         NEW_OCPI = MPOP_OC - MPOP_OCPO
         NEW_BCPI = MPOP_BC - MPOP_BCPO

         ! Add new hydrophilic to old hydrophilic
         ! Don't do this - already added into total for redistribution (clf 3/27/2012)
         ! MPOP_OCPI = MPOP_OCPI + NEW_OCPI
         ! MPOP_BCPI = MPOP_BCPI + NEW_BCPI
         MPOP_OCPI = NEW_OCPI
         MPOP_BCPI = NEW_BCPI


         !==============================================================
         ! CHEMISTRY AND DEPOSITION REACTIONS
         !==============================================================

         IF ( F_PBL < 0.05D0 .OR. 
     &           K_DEPG < SMALLNUM ) THEN

               !==============================================================
               ! Entire box is in the free troposphere
               ! or deposition is turned off, so use RXN without deposition
               ! for gas phase POPs
               ! For particle POPs, no rxn and no deposition
               !==============================================================

            ! OH:
               CALL RXN_OX_NODEP( MPOP_G, K_OX,
     &              E_KOX_T, NEW_POPG, GROSS_OX )

            ! O3 and NO3:
               CALL RXN_OX_NODEP( MPOP_OCPO, K_OX_P, 
     &              E_KOX_T_BC, NEW_POPP_OCPO, GROSS_OX_OCPO)

               CALL RXN_OX_NODEP( MPOP_OCPI, K_OX_P, 
     &              E_KOX_T_BC, NEW_POPP_OCPI, GROSS_OX_OCPI)

               CALL RXN_OX_NODEP( MPOP_BCPO, K_OX_P, 
     &              E_KOX_T_BC, NEW_POPP_BCPO, GROSS_OX_BCPO)

               CALL RXN_OX_NODEP( MPOP_BCPI, K_OX_P, 
     &              E_KOX_T_BC, NEW_POPP_BCPI, GROSS_OX_BCPI)

               ! No deposition occurs [kg]
               DEP_POPG = 0D0
               DEP_POPP_OCPO = 0D0
               DEP_POPP_BCPO = 0D0
               DEP_POPP_OCPI = 0D0
               DEP_POPP_BCPI = 0D0
               

            ELSE IF ( F_PBL > 0.95D0 ) THEN 

               !==============================================================
               ! Entire box is in the boundary layer
               ! so use RXN with deposition for gas phase POPs
               ! Deposition only (no rxn) for particle phase POPs
               !==============================================================

               CALL RXN_OX_WITHDEP( MPOP_G,   K_OX,
     &              K_DEPG,   DTCHEM,  E_KOX_T, NEW_POPG,
     &              GROSS_OX,  DEP_POPG )

!               CALL NO_RXN_WITHDEP( MPOP_OCPO, K_DEPP_OCPO, DTCHEM,
!     &              NEW_POPP_OCPO, DEP_POPP_OCPO )

!               CALL NO_RXN_WITHDEP( MPOP_BCPO, K_DEPP_BCPO, DTCHEM,
!     &              NEW_POPP_BCPO, DEP_POPP_BCPO )

!               CALL NO_RXN_WITHDEP( MPOP_OCPI, K_DEPP_OCPI, DTCHEM,
!     &              NEW_POPP_OCPI, DEP_POPP_OCPI )

!               CALL NO_RXN_WITHDEP( MPOP_BCPI, K_DEPP_BCPI, DTCHEM,
!     &              NEW_POPP_BCPI, DEP_POPP_BCPI )

               CALL RXN_OX_WITHDEP( MPOP_OCPO,   K_OX_P,
     &              K_DEPP_OCPO,   DTCHEM,  E_KOX_T_BC, NEW_POPP_OCPO,
     &              GROSS_OX_OCPO,  DEP_POPP_OCPO )

               CALL RXN_OX_WITHDEP( MPOP_OCPI,   K_OX_P,
     &              K_DEPP_OCPI,   DTCHEM,  E_KOX_T_BC, NEW_POPP_OCPI,
     &              GROSS_OX_OCPI,  DEP_POPP_OCPI )

               CALL RXN_OX_WITHDEP( MPOP_BCPO,   K_OX_P,
     &              K_DEPP_BCPO,   DTCHEM,  E_KOX_T_BC, NEW_POPP_BCPO,
     &              GROSS_OX_BCPO,  DEP_POPP_BCPO )

               CALL RXN_OX_WITHDEP( MPOP_BCPI,   K_OX_P,
     &              K_DEPP_BCPI,   DTCHEM,  E_KOX_T_BC, NEW_POPP_BCPI,
     &              GROSS_OX_BCPI,  DEP_POPP_BCPI )

            ELSE

               !==============================================================
               ! Box spans the top of the boundary layer
               ! Part of the mass is in the boundary layer and subject to 
               ! deposition while part is in the free troposphere and
               ! experiences no deposition.
               !
               ! We apportion the mass between the BL and FT according to the
               ! volume fraction of the box in the boundary layer.
               ! Arguably we should assume uniform mixing ratio, instead of
               ! uniform density but if the boxes are short, the air density
               ! doesn't change much.
               ! But assuming uniform mixing ratio across the inversion layer
               ! is a poor assumption anyway, so we are just using the
               ! simplest approach.
               !==============================================================

               ! Boundary layer portion of POPG [kg]
               POPG_BL = MPOP_G * F_PBL 

               ! Boundary layer portion of POPP_OCPO [kg]
               POPP_OCPO_BL = MPOP_OCPO * F_PBL

               ! Boundary layer portion of POPP_BCPO [kg]
               POPP_BCPO_BL = MPOP_BCPO * F_PBL

               ! Boundary layer portion of POPP_OCPI [kg]
               POPP_OCPI_BL = MPOP_OCPI * F_PBL

               ! Boundary layer portion of POPP_BCPI [kg]
               POPP_BCPI_BL = MPOP_BCPI * F_PBL

               ! Free troposphere portion of POPG [kg]
               POPG_FT = MPOP_G - POPG_BL

               ! Free troposphere portion of POPP_OCPO [kg]
               POPP_OCPO_FT = MPOP_OCPO - POPP_OCPO_BL

               ! Free troposphere portion of POPP_BCPO [kg]
               POPP_BCPO_FT = MPOP_BCPO - POPP_BCPO_BL

               ! Free troposphere portion of POPP_OCPI [kg]
               POPP_OCPI_FT = MPOP_OCPI - POPP_OCPI_BL

               ! Free troposphere portion of POPP_BCPI [kg]
               POPP_BCPI_FT = MPOP_BCPI - POPP_BCPI_BL
               
               ! Do chemistry with deposition on BL fraction for gas phase
               CALL RXN_OX_WITHDEP( POPG_BL,  K_OX,
     &              K_DEPG,   DTCHEM, E_KOX_T,
     &              NEW_POPG, GROSS_OX,  DEP_POPG )           

               ! Do chemistry without deposition on the FT fraction for gas phase
               CALL RXN_OX_NODEP( POPG_FT, K_OX,
     &              E_KOX_T, TMP_POPG, TMP_OX )            

               ! Now do the same with the OC and BC phase:

               ! Do chemistry with deposition on BL fraction for OCPO phase
               CALL RXN_OX_WITHDEP( POPP_OCPO_BL,  K_OX_P,
     &              K_DEPP_OCPO,   DTCHEM, E_KOX_T_BC,
     &              NEW_POPP_OCPO, GROSS_OX_OCPO,  DEP_POPP_OCPO )           

               ! Do chemistry without deposition on the FT fraction for OCPO phase
               CALL RXN_OX_NODEP( POPP_OCPO_FT, K_OX_P,
     &              E_KOX_T_BC, TMP_POPP_OCPO, TMP_OX_P_OCPO )

               CALL RXN_OX_WITHDEP( POPP_OCPI_BL,  K_OX_P,
     &              K_DEPP_OCPI,   DTCHEM, E_KOX_T_BC,
     &              NEW_POPP_OCPI, GROSS_OX_OCPI,  DEP_POPP_OCPI )           

               ! Do chemistry without deposition on the FT fraction for OCPI phase
               CALL RXN_OX_NODEP( POPP_OCPI_FT, K_OX_P,
     &              E_KOX_T_BC, TMP_POPP_OCPI, TMP_OX_P_OCPI )

               ! Do chemistry with deposition on BL fraction for BCPO phase
               CALL RXN_OX_WITHDEP( POPP_BCPO_BL,  K_OX_P,
     &              K_DEPP_BCPO,   DTCHEM, E_KOX_T_BC,
     &              NEW_POPP_BCPO, GROSS_OX_BCPO,  DEP_POPP_BCPO )           

               ! Do chemistry without deposition on the FT fraction for BCPO phase
               CALL RXN_OX_NODEP( POPP_BCPO_FT, K_OX_P,
     &              E_KOX_T_BC, TMP_POPP_BCPO, TMP_OX_P_BCPO ) 

               ! Do chemistry with deposition on BL fraction for BCPI phase
               CALL RXN_OX_WITHDEP( POPP_BCPI_BL,  K_OX_P,
     &              K_DEPP_BCPI,   DTCHEM, E_KOX_T_BC,
     &              NEW_POPP_BCPI, GROSS_OX_BCPI,  DEP_POPP_BCPI )           

               ! Do chemistry without deposition on the FT fraction for BCPI phase
               CALL RXN_OX_NODEP( POPP_BCPI_FT, K_OX_P,
     &              E_KOX_T_BC, TMP_POPP_BCPI, TMP_OX_P_BCPI )    

               ! Do deposition (no chemistry) on BL fraction for particulate phase
               ! No deposition (and no chem) on the FT fraction
               ! for the particulate phase
!               CALL NO_RXN_WITHDEP(POPP_OCPO_BL, K_DEPP_OCPO, DTCHEM,  
!     &              NEW_POPP_OCPO, DEP_POPP_OCPO)

!               CALL NO_RXN_WITHDEP(POPP_BCPO_BL, K_DEPP_BCPO, DTCHEM,  
!     &              NEW_POPP_BCPO, DEP_POPP_BCPO)

!               CALL NO_RXN_WITHDEP(POPP_OCPI_BL, K_DEPP_OCPI, DTCHEM,  
!     &              NEW_POPP_OCPI, DEP_POPP_OCPI)

!               CALL NO_RXN_WITHDEP(POPP_BCPI_BL, K_DEPP_BCPI, DTCHEM,  
!     &              NEW_POPP_BCPI, DEP_POPP_BCPI)
               
               ! Recombine the boundary layer and free troposphere parts [kg]
               NEW_POPG    = NEW_POPG + TMP_POPG
               NEW_POPP_OCPO = NEW_POPP_OCPO + TMP_POPP_OCPO
               NEW_POPP_BCPO = NEW_POPP_BCPO + TMP_POPP_BCPO          
               NEW_POPP_OCPI = NEW_POPP_OCPI + TMP_POPP_OCPI
               NEW_POPP_BCPI = NEW_POPP_BCPI + TMP_POPP_BCPI
               
               ! Total gross oxidation of gas phase in the BL and FT [kg]
               GROSS_OX = GROSS_OX + TMP_OX

               ! Total gross oxidation of particulate phase in the BL and FT [kg]
               GROSS_OX_OCPO = GROSS_OX_OCPO + TMP_OX_P_OCPO
               GROSS_OX_OCPI = GROSS_OX_OCPI + TMP_OX_P_OCPI
               GROSS_OX_BCPO = GROSS_OX_BCPO + TMP_OX_P_BCPO
               GROSS_OX_BCPI = GROSS_OX_BCPI + TMP_OX_P_BCPI

            ENDIF

            ! Ensure positive concentration [kg]
            NEW_POPG    = MAX( NEW_POPG, SMALLNUM )
            NEW_POPP_OCPO = MAX( NEW_POPP_OCPO, SMALLNUM )
            NEW_POPP_BCPO = MAX( NEW_POPP_BCPO, SMALLNUM )
            NEW_POPP_OCPI = MAX( NEW_POPP_OCPI, SMALLNUM )
            NEW_POPP_BCPI = MAX( NEW_POPP_BCPI, SMALLNUM )

            ! Archive new POPG and POPP values [kg]
            STT(I,J,L,IDTPOPG)   = NEW_POPG
            STT(I,J,L,IDTPOPPOCPO) = NEW_POPP_OCPO
            STT(I,J,L,IDTPOPPBCPO) = NEW_POPP_BCPO
            STT(I,J,L,IDTPOPPOCPI) = NEW_POPP_OCPI
            STT(I,J,L,IDTPOPPBCPI) = NEW_POPP_BCPI

            ! Net oxidation [kg] (equal to gross ox for now)
            NET_OX = MPOP_G - NEW_POPG - DEP_POPG
            NET_OX_OCPO = MPOP_OCPO - NEW_POPP_OCPO - DEP_POPP_OCPO
            NET_OX_OCPI = MPOP_OCPI - NEW_POPP_OCPI - DEP_POPP_OCPI
            NET_OX_BCPO = MPOP_BCPO - NEW_POPP_BCPO - DEP_POPP_BCPO 
            NET_OX_BCPI = MPOP_BCPI - NEW_POPP_BCPI - DEP_POPP_BCPI                    

            ! Error check on gross oxidation [kg]
            IF ( GROSS_OX < 0D0 ) 
     &          CALL DEBUG_MSG('CHEM_POPGP: negative gross gas oxid')

            IF ( GROSS_OX_OCPO < 0D0 ) 
     &          CALL DEBUG_MSG('CHEM_POPGP: negative gross OCPO oxid')

            IF ( GROSS_OX_OCPI < 0D0 ) 
     &          CALL DEBUG_MSG('CHEM_POPGP: negative gross OCPI oxid')

            IF ( GROSS_OX_BCPO < 0D0 ) 
     &          CALL DEBUG_MSG('CHEM_POPGP: negative gross BCPO oxid')

            IF ( GROSS_OX_BCPI < 0D0 ) 
     &          CALL DEBUG_MSG('CHEM_POPGP: negative gross BCPI oxid')

            ! Apportion gross oxidation between OH (and no other gas-phase oxidants considered now) [kg]
            IF ( (K_OX     < SMALLNUM) .OR. 
     &           (GROSS_OX < SMALLNUM) ) THEN
               GROSS_OX_OH = 0D0
            ELSE
               GROSS_OX_OH = GROSS_OX * K_OH / K_OX
            ENDIF

            ! Small number check for particulate O3 oxidation
            ! Now apportion total particulate oxidation between O3 and NO3
            IF ( (K_OX_P < SMALLNUM) .OR.
     &           (GROSS_OX_OCPO < SMALLNUM) ) THEN
               GROSS_OX_OCPO = 0D0
            ELSE
               GROSS_OX_O3_OCPO  = GROSS_OX_OCPO * K_O3_BC / K_OX_P
               GROSS_OX_NO3_OCPO = GROSS_OX_OCPO * K_NO3_BC / K_OX_P
            ENDIF

            IF ( (K_OX_P < SMALLNUM) .OR.
     &           (GROSS_OX_OCPI < SMALLNUM) ) THEN
               GROSS_OX_OCPI = 0D0
            ELSE
               GROSS_OX_O3_OCPI  = GROSS_OX_OCPI * K_O3_BC / K_OX_P
               GROSS_OX_NO3_OCPI = GROSS_OX_OCPI * K_NO3_BC / K_OX_P
            ENDIF

            IF ( (K_OX_P < SMALLNUM) .OR.
     &           (GROSS_OX_BCPO < SMALLNUM) ) THEN
               GROSS_OX_BCPO = 0D0
            ELSE
               GROSS_OX_O3_BCPO  = GROSS_OX_BCPO * K_O3_BC / K_OX_P
               GROSS_OX_NO3_BCPO = GROSS_OX_BCPO * K_NO3_BC / K_OX_P
            ENDIF

            IF ( (K_OX_P < SMALLNUM) .OR.
     &           (GROSS_OX_BCPI < SMALLNUM) ) THEN
               GROSS_OX_BCPI = 0D0
            ELSE
               GROSS_OX_O3_BCPI  = GROSS_OX_BCPI * K_O3_BC / K_OX_P
               GROSS_OX_NO3_BCPI = GROSS_OX_BCPI * K_NO3_BC / K_OX_P
            ENDIF

            ! Apportion deposition [kg]
            ! Right now only using dry deposition (no sea salt) (clf, 1/27/11)
            ! If ever use dep with sea salt aerosols,
            ! will need to multiply DEP_POPG by the ratio 
            ! of K_DRYG (rate of dry dep) to K_DEPG (total dep rate).
            IF ( (K_DEPG  < SMALLNUM) .OR. 
     &           (DEP_POPG < SMALLNUM) ) THEN
               DEP_POPG_DRY  = 0D0
            ELSE 
               DEP_POPG_DRY  = DEP_POPG   
            ENDIF

            IF ( (K_DEPP_OCPO  < SMALLNUM) .OR. 
     &           (DEP_POPP_OCPO < SMALLNUM) ) THEN
               DEP_POPP_OCPO_DRY  = 0D0
            ELSE
               DEP_POPP_OCPO_DRY  = DEP_POPP_OCPO
            ENDIF

            IF ( (K_DEPP_BCPO  < SMALLNUM) .OR. 
     &           (DEP_POPP_BCPO < SMALLNUM) ) THEN
               DEP_POPP_BCPO_DRY  = 0D0
            ELSE
               DEP_POPP_BCPO_DRY  = DEP_POPP_BCPO 
            ENDIF

            IF ( (K_DEPP_OCPI  < SMALLNUM) .OR. 
     &           (DEP_POPP_OCPI < SMALLNUM) ) THEN
               DEP_POPP_OCPI_DRY  = 0D0
            ELSE
               DEP_POPP_OCPI_DRY  = DEP_POPP_OCPI
            ENDIF

            IF ( (K_DEPP_BCPI  < SMALLNUM) .OR. 
     &           (DEP_POPP_BCPI < SMALLNUM) ) THEN
               DEP_POPP_BCPI_DRY  = 0D0
            ELSE
               DEP_POPP_BCPI_DRY  = DEP_POPP_BCPI 
            ENDIF

            !=================================================================
            ! ND44 diagnostic: drydep flux of POPG and POPP [molec/cm2/s]
            !=================================================================

            IF ( ( ND44 > 0 .OR. LGTMM ) .AND. (.NOT. LNLPBL) ) THEN
            ! Not using LGTMM right now (logical switch for using GTMM soil model)
            ! Also not using non-local PBL mode yet (clf, 1/27/2011)

               ! Grid box surface area [cm2]
               !AREA_CM2 = GET_AREA_CM2( J )      ! bug fix
               AREA_CM2 = GET_AREA_CM2( I, J, 1 ) ! hma, 8mar2013

               ! Amt of POPG lost to drydep [molec/cm2/s]
               DEP_DRY_FLXG  = DEP_POPG_DRY * XNUMOL(IDTPOPG) / 
     &              ( AREA_CM2 * DTCHEM )

               ! Archive POPG drydep flux in AD44 array [molec/cm2/s]
               AD44(I,J,IDTPOPG,1) = AD44(I,J,IDTPOPG,1) +
     &              DEP_DRY_FLXG

               ! Amt of POPPOC lost to drydep [molec/cm2/s]
               DEP_DRY_FLXP_OCPO = DEP_POPP_OCPO_DRY * 
     &                 XNUMOL(IDTPOPPOCPO)/( AREA_CM2 * DTCHEM )        

               ! Archive POPPOC drydep flux in AD44 array [molec/cm2/s]
               AD44(I,J,IDTPOPPOCPO,1) = 
     &              AD44(I,J,IDTPOPPOCPO,1) + DEP_DRY_FLXP_OCPO

               ! Amt of POPPBC lost to drydep [molec/cm2/s] 
               DEP_DRY_FLXP_BCPO = DEP_POPP_BCPO_DRY * 
     &                 XNUMOL(IDTPOPPBCPO)/( AREA_CM2 * DTCHEM )        

               ! Archive POPPBC drydep flux in AD44 array [molec/cm2/s]
               AD44(I,J,IDTPOPPBCPO,1) = 
     &              AD44(I,J,IDTPOPPBCPO,1) + DEP_DRY_FLXP_BCPO

              ! Amt of POPPOC lost to drydep [molec/cm2/s]
               DEP_DRY_FLXP_OCPI = DEP_POPP_OCPI_DRY * 
     &                 XNUMOL(IDTPOPPOCPI)/( AREA_CM2 * DTCHEM )        

               ! Archive POPPOC drydep flux in AD44 array [molec/cm2/s]
               AD44(I,J,IDTPOPPOCPI,1) = 
     &              AD44(I,J,IDTPOPPOCPI,1) + DEP_DRY_FLXP_OCPI

               ! Amt of POPPBC lost to drydep [molec/cm2/s] 
               DEP_DRY_FLXP_BCPI = DEP_POPP_BCPI_DRY * 
     &                 XNUMOL(IDTPOPPBCPI)/( AREA_CM2 * DTCHEM )        

               ! Archive POPPBC drydep flux in AD44 array [molec/cm2/s]
               AD44(I,J,IDTPOPPBCPI,1) = 
     &              AD44(I,J,IDTPOPPBCPI,1) + DEP_DRY_FLXP_BCPI


            ENDIF
           

            !==============================================================
            ! ND53 diagnostic: Oxidized POP production [kg]
            !==============================================================

            IF ( ND53 > 0 .AND. L <= LD53 ) THEN ! LD53 is max level

               ! OH:
               AD53_POPG_OH(I,J,L)= AD53_POPG_OH(I,J,L) + GROSS_OX_OH
               
               ! O3:
               AD53_POPP_OCPO_O3(I,J,L)=AD53_POPP_OCPO_O3(I,J,L)
     &             + GROSS_OX_O3_OCPO
               AD53_POPP_OCPI_O3(I,J,L)=AD53_POPP_OCPI_O3(I,J,L) 
     &             + GROSS_OX_O3_OCPI
               AD53_POPP_BCPO_O3(I,J,L)=AD53_POPP_BCPO_O3(I,J,L) 
     &             + GROSS_OX_O3_BCPO
               AD53_POPP_BCPI_O3(I,J,L)=AD53_POPP_BCPI_O3(I,J,L) 
     &             + GROSS_OX_O3_BCPI

               ! NO3:
               AD53_POPP_OCPO_NO3(I,J,L)=AD53_POPP_OCPO_NO3(I,J,L)
     &             + GROSS_OX_NO3_OCPO
               AD53_POPP_OCPI_NO3(I,J,L)=AD53_POPP_OCPI_NO3(I,J,L) 
     &             + GROSS_OX_NO3_OCPI
               AD53_POPP_BCPO_NO3(I,J,L)=AD53_POPP_BCPO_NO3(I,J,L) 
     &             + GROSS_OX_NO3_BCPO
               AD53_POPP_BCPI_NO3(I,J,L)=AD53_POPP_BCPI_NO3(I,J,L) 
     &             + GROSS_OX_NO3_BCPI

            ENDIF

      ENDDO
      ENDDO
      ENDDO



! END OMP stuff here if added

      END SUBROUTINE CHEM_POPGP 

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  RXN_OX_NODEP
!
! !DESCRIPTION: Subroutine RXN_OX_NODEP calculates new mass of POPG for given
! oxidation rates, without any deposition. This is for the free troposphere, or
! simulations with deposition turned off. (clf, 1/27/11, based on RXN_REDOX_NODEP
! in mercury_mod.f).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RXN_OX_NODEP( OLD_POPG, K_OX, E_KOX_T,
     &     NEW_POPG, GROSS_OX )
!

! !INPUT PARAMETERS: 
      REAL*8,  INTENT(IN)  :: OLD_POPG
      REAL*8,  INTENT(IN)  :: K_OX
      REAL*8,  INTENT(IN)  :: E_KOX_T
      
!
! !INPUT/OUTPUT PARAMETERS:   
!
!
! !OUTPUT PARAMETERS:
      REAL*8,  INTENT(OUT) :: NEW_POPG,  GROSS_OX
!
! !REVISION HISTORY: 
!  27 January 2011 - CL Friedman - Initial Version
!
! !REMARKS:
! (1) Based on RXN_REDOX_NODEP in mercury_mod.f
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine RXN_OX_NODEP calculates new mass of POPG for given
! oxidation rates, without any deposition. This is for the free troposphere, or
! simulations with deposition turned off. (clf, 1/27/11, based on RXN_REDOX_NODEP
! in mercury_mod.f).
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) OLD_POPG (REAL*8) : 
!  (2 ) DT       (REAL*8) : 
!  (3 ) K_OX     (REAL*8) :
!  (4 ) E_KOX_T  (REAL*8) :
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) NEW_POPG (REAL*8) :
!  (2 ) GROSS_OX (REAL*8) :
!
!  Local variables:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  (1 ) 
!  
!  REFS:
!  (1 )  
!
!******************************************************************************
!BOC
      
      ! Local variables
      ! None

      !=================================================================
      ! RXN_OX_NODEP begins here!
      !=================================================================

         !=================================================================
         ! Oxidation
         !=================================================================

         IF (K_OX < SMALLNUM ) THEN

            GROSS_OX = 0d0
            NEW_POPG = OLD_POPG

         ELSE 

         ! New concentration of POPG
         NEW_POPG = OLD_POPG * E_KOX_T

         ! Gross oxidation 
         GROSS_OX = OLD_POPG - NEW_POPG
         GROSS_OX = MAX( GROSS_OX, 0D0 )

         ENDIF

      END SUBROUTINE RXN_OX_NODEP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  RXN_OX_WITHDEP
!
! !DESCRIPTION: Subroutine RXN_OX_WITHDEP calculates new mass of POPG for given
! rates of oxidation and deposition. This is for the boundary layer.
! (clf, 1/27/11, based on RXN_REDOX_NODEP in mercury_mod.f).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RXN_OX_WITHDEP( OLD_POPG, K_OX, K_DEPG, DT, E_KOX_T,
     &     NEW_POPG, GROSS_OX, DEP_POPG )
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP


! !INPUT PARAMETERS: 
      REAL*8,  INTENT(IN)  :: OLD_POPG,  DT
      REAL*8,  INTENT(IN)  :: K_OX, K_DEPG
      REAL*8,  INTENT(IN)  :: E_KOX_T
      
!
! !INPUT/OUTPUT PARAMETERS:   
!
!
! !OUTPUT PARAMETERS:
      REAL*8,  INTENT(OUT) :: NEW_POPG,  GROSS_OX
      REAL*8,  INTENT(OUT) :: DEP_POPG
!
! !REVISION HISTORY: 
!  27 January 2011 - CL Friedman - Initial Version
!
! !REMARKS:
! (1) Based on RXN_REDOX_WITHDEP in mercury_mod.f
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine RXN_OX_WITHDEP calculates new mass of POPG for given
! rates of oxidation and deposition. This is for the boundary layer.
! (clf, 1/27/11, based on RXN_REDOX_WITHDEP in mercury_mod.f).
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) OLD_POPG (REAL*8) : 
!  (2 ) DT       (REAL*8) : 
!  (3 ) K_OX     (REAL*8) :
!  (4 ) K_DEPG   (REAL*8) :
!  (5 ) E_KOX_T  (REAL*8) :
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) NEW_POPG (REAL*8) :
!  (2 ) GROSS_OX (REAL*8) :
!  (3 ) DEP_POPG (REAL*8) :
!
!  Local variables:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  (1 ) 
!  
!  REFS:
!  (1 )  
!
!******************************************************************************
!BOC
      
      ! Local variables
      REAL*8               :: E_KDEPG_T
      REAL*8               :: NEWPOPG_OX
      REAL*8               :: NEWPOPG_DEP

      !=================================================================
      ! RXN_OX_WITHDEP begins here!
      !=================================================================

      ! Precompute exponential factor for deposition [dimensionless]
      E_KDEPG_T = EXP( -K_DEPG * DT )

      IF (K_OX < SMALLNUM) THEN     
         
         !=================================================================
         ! No Chemistry, Deposition only
         !=================================================================

         ! New mass of POPG [kg]
         NEW_POPG = OLD_POPG * E_KDEPG_T
         
         ! Oxidation of POPG [kg]
         GROSS_OX = 0D0

         ! Deposited POPG [kg]
         DEP_POPG = OLD_POPG - NEW_POPG

      ELSE

         !=================================================================
         ! Oxidation and Deposition 
         !=================================================================

         ![POPG](t) = [POPG](0) exp( -(kOx + kDPOPG) t)
         !Ox(t)     = ( [POPG](0) - [POPG](t) ) * kOx / ( kOx + kDPOPG )
         !Dep_POPG(t)   = ( [POPG](0) - [POPG](t) - Ox(t) ) 

         ! New concentration of POPG [kg]
         NEW_POPG = OLD_POPG * E_KOX_T * E_KDEPG_T

         ! Gross oxidized gas phase mass [kg]
         GROSS_OX = ( OLD_POPG - NEW_POPG ) * K_OX / ( K_OX + K_DEPG )
         GROSS_OX = MAX( GROSS_OX, 0D0 )

         ! POPG deposition [kg]
         DEP_POPG = ( OLD_POPG - NEW_POPG - GROSS_OX )       
         DEP_POPG = MAX( DEP_POPG, 0D0 )

      ENDIF

      END SUBROUTINE RXN_OX_WITHDEP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  NO_RXN_WITHDEP
!
! !DESCRIPTION: Subroutine NO_RXN_WITHDEP calculates new mass of POPP for given
! rate of deposition. No oxidation of POPP. This is for the boundary layer.
! (clf, 2/9/11)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NO_RXN_WITHDEP( OLD_POPP, K_DEPP, DT,
     &     NEW_POPP, DEP_POPP )
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP


! !INPUT PARAMETERS: 
      REAL*8,  INTENT(IN)  :: OLD_POPP
      REAL*8,  INTENT(IN)  :: K_DEPP
      REAL*8,  INTENT(IN)  :: DT
      
!
! !INPUT/OUTPUT PARAMETERS:   
!
!
! !OUTPUT PARAMETERS:
      REAL*8,  INTENT(OUT) :: NEW_POPP
      REAL*8,  INTENT(OUT) :: DEP_POPP
!
! !REVISION HISTORY: 
!  9 February 2011 - CL Friedman - Initial Version
!
! !REMARKS:
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine NO_RXN_WITHDEP calculates new mass of POPP for given
! rate of deposition. This is for the boundary layer.
! (clf, 1/27/11, based on RXN_REDOX_NODEP in mercury_mod.f).
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) OLD_POPP (REAL*8) : 
!  (2 ) K_DEPP   (REAL*8) :
!  (3 ) DT       (REAL*8) :
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) NEW_POPP (REAL*8) :
!  (2 ) DEP_POPP (REAL*8) :
!
!  Local variables:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  (1 ) 
!  
!  REFS:
!  (1 )  
!
!******************************************************************************
!BOC
      
      ! Local variables
      REAL*8               :: E_KDEPP_T

      !=================================================================
      ! NO_RXN_WITHDEP begins here!
      !=================================================================

      ! Precompute exponential factors [dimensionless]
      E_KDEPP_T = EXP( -K_DEPP * DT )     

      !=================================================================
      ! No Chemistry, Deposition only
      !=================================================================

      ! New mass of POPP [kg]
      NEW_POPP = OLD_POPP * E_KDEPP_T

      ! POPP deposition [kg]
      DEP_POPP = OLD_POPP - NEW_POPP
      DEP_POPP = MAX( DEP_POPP, 0D0 )


      END SUBROUTINE NO_RXN_WITHDEP
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  EMISSPOPS
!
! !DESCRIPTION: This routine is the driver routine for POPs emissions
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISSPOPS
!
! !INPUT PARAMETERS: 
!

!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on EMISSMERCURY from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine EMISSPOPS is the driver subroutine for POPs emissions.
!  !
!  Arguments as Input:
!  ============================================================================
!  (1 )
!
!  Local variables:
!  ============================================================================
!  (1 ) I, J, L   (INTEGER) : Long, lat, level
!  (2 ) N         (INTEGER) : Tracer ID
!  (3 ) PBL_MAX   (INTEGER) : Maximum extent of boundary layer [level]
!  (4 ) DTSRCE    (REAL*8)  : Emissions time step  [s]
!  (5 ) T_POP     (REAL*8)  : POP emission rate [kg/s]
!  (6 ) E_POP     (REAL*8)  : POPs emitted into box [kg]
!  (7 ) F_OF_PBL  (REAL*8)  : Fraction of box within boundary layer [unitless]
!     
!  NOTES:
!  (1 ) 
!  
!  REFS:
!  (1 )
!******************************************************************************
!BOC
      ! References to F90 modules
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP
      USE LOGICAL_MOD,       ONLY : LPRT, LNLPBL !CDH added LNLPBL
      USE TIME_MOD,          ONLY : GET_MONTH, ITS_A_NEW_MONTH, GET_YEAR
      USE TRACER_MOD,        ONLY : STT
      USE VDIFF_PRE_MOD,     ONLY : EMIS_SAVE !cdh for LNLPBL
      USE GRID_MOD,          ONLY : GET_XMID, GET_YMID
      USE DAO_MOD,           ONLY : T, AIRVOL, SNOMAS
      USE DAO_MOD,           ONLY : IS_ICE, IS_LAND, IS_WATER
      USE DAO_MOD,           ONLY : FRLAKE
!LAI calls obsolete?
! Update to new LAI naming convention. (H.Amos, 7Mar2013)
!      USE MODIS_LAI_MOD,     ONLY : GC_LAI, GC_LAI_PM, GC_LAI_CM !hma, 7Mar2013
!      USE LAI_MOD,           ONLY : INIT_LAI, DAYS_BTW_MON !hma, 7Mar2013
!      USE LAI_MOD,           ONLY : ISOLAI, MISOLAI, PMISOLAI
!      USE LAI_MOD,           ONLY : INIT_LAI, DAYS_BTW_M
      USE GLOBAL_OC_MOD,     ONLY : GET_GLOBAL_OC
      USE GLOBAL_BC_MOD,     ONLY : GET_GLOBAL_BC
      ! Reference to diagnostic arrays
      !USE DIAG53_MOD,   ONLY : AD53, ND53
      USE DIAG53_MOD,   ONLY : AD53 !hma, 8mar2013
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRACERID_MOD, ONLY : IDTPOPG, IDTPOPPOCPO,  IDTPOPPBCPO
      USE LAND_POPS_MOD, ONLY: SOILEMISPOP, LAKEEMISPOP, VEGEMISPOP
      USE GET_POPSINFO_MOD, ONLY: GET_POP_DEL_H, GET_POP_KOA
      USE GET_POPSINFO_MOD, ONLY: GET_POP_KBC, GET_POP_K_POPG_OH
      USE GET_POPSINFO_MOD, ONLY: GET_POP_K_POPP_O3A
      USE GET_POPSINFO_MOD, ONLY: GET_POP_K_POPP_O3B

      ! Replace include statements. I also think CMN_DEP is
      ! obsolete (H. Amos, 7 Mar 2013)
      USE CMN_SIZE_MOD        ! Size parameters  
      USE CMN_DIAG_MOD        ! Diagnostic flags    
!#     include "CMN_SIZE"     ! Size parameters
!#     include "CMN_DEP"      ! FRCLND

      ! Local variables
      INTEGER               :: I,   J,    L,    N,    PBL_MAX
      INTEGER               :: MONTH, YEAR
      REAL*8                :: DTSRCE, F_OF_PBL, TK
      REAL*8                :: E_POP, T_POP
      REAL*8                :: C_OC1,   C_BC1,  AIR_VOL
      REAL*8                :: C_OC2,   C_BC2
      REAL*8                :: F_POP_OC, F_POP_BC 
      REAL*8                :: F_POP_G
      REAL*8                :: KOA_T, KBC_T, KOC_BC_T, KBC_OC_T
      REAL*8                :: VR_OC_AIR, VR_OC_BC
      REAL*8                :: VR_BC_AIR, VR_BC_OC, SUM_F
      REAL*8                :: OC_AIR_RATIO, OC_BC_RATIO
      REAL*8                :: BC_AIR_RATIO, BC_OC_RATIO
      REAL*8                :: MAXVAL_EMISSPOPS
      REAL*8                :: MINVAL_EMISSPOPS
      REAL*8                :: FRAC_SNOW_OR_ICE, FRAC_SNOWFREE_LAND
      REAL*8                :: FRAC_LEAF, FRAC_LAKE, FRAC_SOIL
      LOGICAL, SAVE         :: FIRST = .TRUE.
      LOGICAL               :: IS_SNOW_OR_ICE, IS_LAND_OR_ICE

      REAL*8,    SAVE       :: POP_XMW, POP_KOA, POP_KBC, POP_K_POPG_OH
      REAL*8,    SAVE       :: POP_K_POPP_O3A, POP_K_POPP_O3B
      REAL*8,    SAVE       :: POP_HSTAR, POP_DEL_H, POP_DEL_Hw

      ! Delta H for POP [kJ/mol]. Delta H is enthalpy of phase transfer
      ! from gas phase to OC. For now we use Delta H for phase transfer 
      ! from the gas phase to the pure liquid state. 
      ! For PHENANTHRENE: 
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -74000 [J/mol].
      ! For PYRENE:
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -87000 [J/mol].    
      ! For BENZO[a]PYRENE:
      ! this is also taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 452, Prob 11.1), or -110,000 [J/mol]
      REAL*8     :: DEL_H

      ! R = universal gas constant for adjusting KOA for temp: 8.3145 [J/mol/K]
      REAL*8, PARAMETER     :: R          = 8.31d0  

      ! KOA_298 for partitioning of gas phase POP to atmospheric OC
      ! KOA_298 = Cpop in octanol/Cpop in atmosphere at 298 K 
      ! For PHENANTHRENE:
      ! log KOA_298 = 7.64, or 4.37*10^7 [unitless]
      ! For PYRENE:
      ! log KOA_298 = 8.86, or 7.24*10^8 [unitless]
      ! For BENZO[a]PYRENE:
      ! log KOA_298 = 11.48, or 3.02*10^11 [unitless]
      ! (Ma et al., J. Chem. Eng. Data, 2010, 55:819-825).
      REAL*8     :: KOA_298

      ! KBC_298 for partitioning of gas phase POP to atmospheric BC
      ! KBC_298 = Cpop in black carbon/Cpop in atmosphere at 298 K
      ! For PHENANTHRENE:
      ! log KBC_298 = 10.0, or 1.0*10^10 [unitless]
      ! For PYRENE:
      ! log KBC_298 = 11.0, or 1.0*10^11 [unitless]
      ! For BENZO[a]PYRENE:
      ! log KBC_298 = 13.9, or 7.94*10^13 [unitless]
      ! (Lohmann and Lammel, EST, 2004, 38:3793-3802)
      REAL*8     :: KBC_298

      ! DENS_OCT = density of octanol, needed for partitioning into OC
      ! 820 [kg/m^3]
      REAL*8, PARAMETER     :: DENS_OCT   = 82d1

      ! DENS_BC = density of BC, needed for partitioning onto BC
      ! 1 [kg/L] or 1000 [kg/m^3]
      ! From Lohmann and Lammel, Environ. Sci. Technol., 2004, 38:3793-3803.
      REAL*8, PARAMETER     :: DENS_BC    = 1d3

      REAL*8                :: DUM

      !=================================================================
      ! EMISSPOPS begins here!
      !=================================================================

      DUM = 1.0
      DEL_H = GET_POP_DEL_H(DUM)
      KOA_298 = GET_POP_KOA(DUM)
      KBC_298 = GET_POP_KBC(DUM)


      CALL INIT_POPS(DUM,DUM,DUM,DUM,DUM,DUM,DUM,DUM,DUM)

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Read primary emissions from disk
         !CALL INIT_POPS 
         CALL POPS_READYR

         ! Read land concentrations from disk
         CALL POPS_READSURFCONC

         ! Read soil organic carbon fractions from disk
         CALL POPS_READ_FOC

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Read monthly OC and BC fields for gas-particle partitioning
      !=================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN 

         ! Get the current month
         MONTH = GET_MONTH()

         ! Get the current month
         YEAR = GET_YEAR()

         ! Read monthly OC and BC from disk
         CALL GET_GLOBAL_OC( MONTH, YEAR )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMPOPS: a GET_GLOBAL_OC' )

         CALL GET_GLOBAL_BC( MONTH, YEAR )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMPOPS: b GET_GLOBAL_BC' )

      ENDIF

      ! If we are using the non-local PBL mixing,
      ! we need to initialize the EMIS_SAVE array (cdh, 08/27/09)
      IF (LNLPBL) EMIS_SAVE = 0d0

      ! Emission timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0

      ! Maximum extent of the PBL [model level]
      PBL_MAX = GET_PBL_MAX_L() 
 
      !=================================================================
      ! Call emissions routines for revolatilization fluxes from surfaces
      ! Assume all re-emited POPs are in the gas phase until partitioning
      ! to ambient OC and BC in the boundary layer.
      ! Re-emission flux/mass depends on type of POP
      ! First draft, CLF, 28 Aug 2012
      !=================================================================

       CALL SOILEMISPOP( POP_SURF, F_OC_SOIL, EPOP_SOIL)
       CALL LAKEEMISPOP( POP_SURF, EPOP_LAKE )
       CALL VEGEMISPOP(  POP_SURF, EPOP_VEG )

       EPOP_SNOW = 0d0
       EPOP_OCEAN = 0d0


      ! Loop over grid boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR          

!            CALL FLUSH (6)  

         F_OF_PBL = 0d0 
         T_POP = 0d0       

         !Here, save the total from the emissions array
         !into the T_POP variable [kg/s]
         T_POP = POP_TOT_EM(I,J)

         ! Now add revolatilization (secondary) emissions to primary [kg/s]
         T_POP = T_POP + EPOP_VEG(I,J) + EPOP_LAKE(I,J) + EPOP_SOIL(I,J) 
!     &         +  EPOP_SNOW(I,J) + EPOP_OCEAN(I,J)

         !==============================================================
         ! Apportion total POPs emitted to gas phase, OC-bound, and BC-bound
         ! emissions (clf, 2/1/2011)         
         ! Then partition POP throughout PBL; store into STT [kg]
         ! Now make sure STT does not underflow (cdh, bmy, 4/6/06; eck 9/20/10)
         !==============================================================

         ! Loop up to max PBL level
         DO L = 1, PBL_MAX

            !Get temp [K]
            TK = T(I,J,L)

            ! Define temperature-dependent partition coefficients:
            ! KOA_T, the octanol-air partition coeff at temp T [unitless]
            KOA_T = KOA_298 * EXP((-DEL_H/R) * ((1d0/TK) - 
     &              (1d0/298d0)))

            ! Define KBC_T, the BC-air partition coeff at temp T [unitless]
            ! TURN OFF TEMPERATURE DEPENDENCY FOR SENSITIVITY ANALYSIS
            KBC_T = KBC_298 * EXP((-DEL_H/R) * ((1d0/TK) - 
     &              (1d0/298d0)))

            ! Define KOC_BC_T, the theoretical OC-BC part coeff at temp T [unitless]
            KOC_BC_T = KOA_T / KBC_T

            ! Define KBC_OC_T, the theoretical BC_OC part coeff at temp T [unitless]
            KBC_OC_T = 1d0 / KOC_BC_T

           ! Get monthly mean OC and BC concentrations [kg/box]
            C_OC1        = GET_OC( I, J, L )
            C_BC1        = GET_BC( I, J, L )
           
            ! Convert C_OC and C_BC units to volume per box 
            ! [m^3 OC or BC/box]
            C_OC2        = C_OC1 / DENS_OCT
            C_BC2        = C_BC1 / DENS_BC

            ! Get air volume (m^3)
            AIR_VOL     = AIRVOL(I,J,L) 

            ! Define volume ratios:
            ! VR_OC_AIR = volume ratio of OC to air [unitless]    
            VR_OC_AIR = C_OC2 / AIR_VOL

            ! VR_OC_BC  = volume ratio of OC to BC [unitless]
            VR_OC_BC    = C_OC2 / C_BC2

            ! VR_BC_AIR = volume ratio of BC to air [unitless]
            VR_BC_AIR   = VR_OC_AIR / VR_OC_BC

            ! VR_BC_OC  = volume ratio of BC to OC [unitless]
            !VR_BC_OC(I,J,L)    = 1d0 / VR_OC_BC(I,J,L)
            VR_BC_OC    = 1d0 / VR_OC_BC 

            ! Redefine fractions of total POPs in box (I,J,L) that are OC-phase, 
            ! BC-phase, and gas phase with new time step (should only change if 
            ! temp changes or OC/BC concentrations change) 
            OC_AIR_RATIO = 1d0 / (KOA_T * VR_OC_AIR) 
            OC_BC_RATIO = 1d0 / (KOC_BC_T * VR_OC_BC) 
  
            BC_AIR_RATIO = 1d0 / (KBC_T * VR_BC_AIR) 
            BC_OC_RATIO = 1d0 / (KBC_OC_T * VR_BC_OC)

            ! If there are zeros in OC or BC concentrations, make sure they
            ! don't cause problems with phase fractions
            IF ( C_OC1 > SMALLNUM .and. C_BC1 > SMALLNUM ) THEN
               F_POP_OC  = 1d0 / (1d0 + OC_AIR_RATIO + OC_BC_RATIO) 
               F_POP_BC  = 1d0 / (1d0 + BC_AIR_RATIO + BC_OC_RATIO)
         
             ELSE IF (C_OC1 > SMALLNUM .and.
     &             C_BC1 .le. SMALLNUM ) THEN
             F_POP_OC  = 1d0 / (1d0 + OC_AIR_RATIO)
             F_POP_BC  = SMALLNUM           

             ELSE IF ( C_OC1 .le. SMALLNUM .and.
     &             C_BC1 > SMALLNUM ) THEN
             F_POP_OC  = SMALLNUM
             F_POP_BC  = 1d0 / (1d0 + BC_AIR_RATIO)

             ELSE IF ( C_OC1 .le. SMALLNUM .and. 
     &             C_BC1 .le. SMALLNUM) THEN
             F_POP_OC = SMALLNUM
             F_POP_BC = SMALLNUM
            ENDIF

            ! Gas-phase:
            F_POP_G   = 1d0 - F_POP_OC - F_POP_BC

            ! Check that sum of fractions equals 1
            SUM_F = F_POP_OC + F_POP_BC + F_POP_G                
            
            ! Fraction of PBL that box (I,J,L) makes up [unitless]
            F_OF_PBL    = GET_FRAC_OF_PBL(I,J,L)

            ! Calculate rates of POP emissions in each phase [kg/s]
            ! OC-phase:
            EPOP_OC(I,J,L) = F_POP_OC * F_OF_PBL * T_POP                        

            ! BC-phase
            EPOP_BC(I,J,L) = F_POP_BC * F_OF_PBL * T_POP                         

            ! Gas-phase
            EPOP_G(I,J,L)  = F_POP_G * F_OF_PBL * T_POP

            ! Toggle for turning off emissions in specific regions

            ! USA: 
!            IF ( GET_XMID(I) > -125 .AND. GET_XMID(I) < -65 .AND.
!     &           GET_YMID(J) >  25  .AND. GET_YMID(J) < 50 ) THEN
!
!               EPOP_OC(I,J,L) = 0d0
!               EPOP_BC(I,J,L) = 0d0
!               EPOP_G(I,J,L)  = 0d0
!
!            ENDIF


            ! Canada (includes Alaska): 
!            IF ( GET_XMID(I) > -170 .AND. GET_XMID(I) < -57 .AND.
!     &           GET_YMID(J) >=  50  .AND. GET_YMID(J) < 75 ) THEN
!
!               EPOP_OC(I,J,L) = 0d0
!               EPOP_BC(I,J,L) = 0d0
!               EPOP_G(I,J,L)  = 0d0
!
!            ENDIF


            ! Europe: 
!            IF ( GET_XMID(I) > -12.5 .AND. GET_XMID(I) < 40 .AND.
!     &           GET_YMID(J) >  36  .AND. GET_YMID(J) < 78 ) THEN
!
!               EPOP_OC(I,J,L) = 0d0
!               EPOP_BC(I,J,L) = 0d0
!               EPOP_G(I,J,L)  = 0d0
!
!            ENDIF

            ! Asia:
!            IF ( GET_XMID(I) >= 70.0 .and. GET_XMID(I) < 152.5 .and.
!     &           GET_YMID(J) >=  8.0 .and. GET_YMID(J) < 51.0 ) THEN
!
!               EPOP_OC(I,J,L) = 0d0
!               EPOP_BC(I,J,L) = 0d0
!               EPOP_G(I,J,L)  = 0d0
!
!            ENDIF         
        

            ! Africa:
!            IF ( GET_XMID(I) >= -18 .and. GET_XMID(I) < 49 .and.
!     &           GET_YMID(J) >= -41 .and. GET_YMID(J) < 33 ) THEN
!
!               EPOP_OC(I,J,L) = 0d0
!               EPOP_BC(I,J,L) = 0d0
!               EPOP_G(I,J,L)  = 0d0
!
!            ENDIF 

            ! Russia:
!            IF ( GET_XMID(I) >= 40 .and. GET_XMID(I) < 180 .and.
!     &           GET_YMID(J) >= 51 .and. GET_YMID(J) < 75 ) THEN
!
!              EPOP_OC(I,J,L) = 0d0
!               EPOP_BC(I,J,L) = 0d0
!               EPOP_G(I,J,L)  = 0d0
!
!            ENDIF       

     
!     TURNING OFF EMISSIONS IN LRTAP SOURCE REGIONS:
            ! Europe including North Africa:
!            IF ( GET_XMID(I) >= -10 .and. GET_XMID(I) <= 50 .and.
!     &           GET_YMID(J) >= 25 .and. GET_YMID(J) <= 65 ) THEN
!
!               EPOP_OC(I,J,L) = 0d0
!               EPOP_BC(I,J,L) = 0d0
!               EPOP_G(I,J,L)  = 0d0
!
!            ENDIF 

            ! North America:
!            IF ( GET_XMID(I) >= -125 .and. GET_XMID(I) <= -60 .and.
!     &           GET_YMID(J) >= 15 .and. GET_YMID(J) <= 55 ) THEN

!               EPOP_OC(I,J,L) = 0d0
!               EPOP_BC(I,J,L) = 0d0
!               EPOP_G(I,J,L)  = 0d0

!            ENDIF 

            ! East Asia:
!            IF ( GET_XMID(I) >= 95 .and. GET_XMID(I) <= 160 .and.
!     &           GET_YMID(J) >= 15 .and. GET_YMID(J) <= 50 ) THEN
!
!               EPOP_OC(I,J,L) = 0d0
!               EPOP_BC(I,J,L) = 0d0
!               EPOP_G(I,J,L)  = 0d0
!
!            ENDIF  

            ! South Asia:
!            IF ( GET_XMID(I) >= 50 .and. GET_XMID(I) <= 95 .and.
!     &           GET_YMID(J) >= 5 .and. GET_YMID(J) <= 35 ) THEN
!
!               EPOP_OC(I,J,L) = 0d0
!               EPOP_BC(I,J,L) = 0d0
!               EPOP_G(I,J,L)  = 0d0
!
!            ENDIF 

            ! Russia:
!            IF ( GET_XMID(I) >= 50 .and. GET_XMID(I) <= 180 .and.
!     &           GET_YMID(J) >= 51 .and. GET_YMID(J) <= 75 ) THEN
!
!              EPOP_OC(I,J,L) = 0d0
!              EPOP_BC(I,J,L) = 0d0
!               EPOP_G(I,J,L)  = 0d0
!
!            ENDIF  

            !-----------------
            ! OC-PHASE EMISSIONS
            !-----------------
            N           = IDTPOPPOCPO
            E_POP       = EPOP_OC(I,J,L) * DTSRCE
            CALL EMITPOP( I, J, L, N, E_POP )

            !-----------------
            ! BC-PHASE EMISSIONS
            !-----------------
            N           = IDTPOPPBCPO
            E_POP       = EPOP_BC(I,J,L) * DTSRCE
            CALL EMITPOP( I, J, L, N, E_POP )

            !-----------------
            ! GASEOUS EMISSIONS
            !-----------------
            N           = IDTPOPG
            E_POP       = EPOP_G(I,J,L) * DTSRCE
            CALL EMITPOP( I, J, L, N, E_POP )
             
            ENDDO


         !==============================================================
         ! Sum different POPs emissions phases (OC, BC, and gas phase)
         ! through bottom layer to top of PBL for storage in ND53 diagnostic
         !==============================================================

           SUM_OC_EM(I,J) =  SUM(EPOP_OC(I,J,1:PBL_MAX))  
           SUM_BC_EM(I,J) =  SUM(EPOP_BC(I,J,1:PBL_MAX))
           SUM_G_EM(I,J)  =  SUM(EPOP_G(I,J,1:PBL_MAX))           
       
         SUM_OF_ALL(I,J) = SUM_OC_EM(I,J) + SUM_BC_EM(I,J) + 
     &                      SUM_G_EM(I,J)

         ! Check that sum thru PBL is equal to original emissions array
         SUM_OF_ALL(I,J) = POP_TOT_EM(I,J) / SUM_OF_ALL(I,J)
        

         !==============================================================
         ! ND53 diagnostic: POP emissions [kg]
         ! 1 = total;  2 = OC;  3 = BC;  4 = gas phase
         !==============================================================

         IF ( ND53 > 0 ) THEN
            AD53(I,J,1) = AD53(I,J,1) + (T_POP * DTSRCE)
            AD53(I,J,2) = AD53(I,J,2) + (SUM_OC_EM(I,J) * DTSRCE)
            AD53(I,J,3) = AD53(I,J,3) + (SUM_BC_EM(I,J) * DTSRCE)
            AD53(I,J,4) = AD53(I,J,4) + (SUM_G_EM(I,J) * DTSRCE)
         ENDIF
         
      ENDDO
      ENDDO


      ! Return to calling program
      END SUBROUTINE EMISSPOPS

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  EMITPOP
!
! !DESCRIPTION: This routine directs emission either to STT directly or to EMIS_SAVE
!  for use by the non-local PBL mixing. This is a programming convenience.
!  (cdh, 08/27/09, modified for pops by eck, 9/20/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMITPOP( I, J, L, ID, E_POP )
!! 
! !USES:
      ! Reference to diagnostic arrays
      USE TRACER_MOD,   ONLY : STT
      USE LOGICAL_MOD,  ONLY : LNLPBL
      USE VDIFF_PRE_MOD,ONLY : EMIS_SAVE
! !INPUT PARAMETERS: 
      INTEGER, INTENT(IN)   :: I, J, L, ID
      REAL*8,  INTENT(IN)   :: E_POP
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on EMITHG from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!******************************************************************************
!Comment header
!  Subroutine EMITPOP directs emission either to STT directly or to EMIS_SAVE
!  for use by the non-local PBL mixing. This is a programming convenience.
!  (cdh, 08/27/09, modified for pops by eck, 9/20/10)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) I, J, L            INTEGERS  Grid box dimensions
!  (2 ) ID                 INTEGER   Tracer ID
!  (3 ) E_POP              REAL*8    POP emissions [kg/s]
!
!  Local variables:
!  ============================================================================
!  (1 ) 
!     
!  NOTES:
!  (1 ) Based on EMITHG in mercury_mod.f
!  
!  REFS:
!  (1 )

!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! EMITPOP begins here!
      !=================================================================

      ! Save emissions [kg/s] for non-local PBL mixing or emit directly.
      ! Make sure that emitted mass is non-negative
      ! This is here only for consistency with old code which warned of
      ! underflow error (cdh, 08/27/09, modified for POPs 9/20/10)
      IF (LNLPBL) THEN
         EMIS_SAVE(I,J,ID) = EMIS_SAVE(I,J,ID) + MAX( E_POP, 0D0 )
      ELSE
          STT(I,J,L,ID) = STT(I,J,L,ID) + E_POP
      ENDIF

      END SUBROUTINE EMITPOP
!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  POPS_READYR
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE POPS_READYR
!! 
! !USES:
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : READ_BPCH2, GET_TAU0
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1, POP_EMISDIR
      USE REGRID_1x1_MOD,    ONLY : DO_REGRID_1x1
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE GET_POPSINFO_MOD,  ONLY : GET_EMISSFILE
    
! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  3 February 2011 - CL Friedman - Initial Version
!
! !REMARKS:
! (1) Based initially on MERCURY_READYR from MERCURY_MOD (clf, 2/3/2011)
!
!EOP
!******************************************************************************
!Comment header
!  Subroutine POP_READYR read the year-invariant emissions for POPs (PAHs) from 
!  all sources combined
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) 
!
!  Local variables:
!  ============================================================================
!  (1 ) 
!     
!  NOTES:
!  (1 ) Based on MERCURY_READYR in mercury_mod.f
!  
!  REFS:
!  (1 ) Zhang, Y. and Tao, S. 2009. Global atmospheric emission inventory
!  of polycyclic aromatic hydrocarbons (PAHs) for 2004. Atm Env. 43:812-819.

!------------------------------------------------------------------------------
!BOC
      USE CMN_SIZE_MOD             ! Size parameters ! hma 7Mar2013
!#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      REAL*4               :: ARRAY(I1x1,J1x1,1) 
      REAL*8               :: XTAU, MAX_A, MIN_A
      REAL*8               :: MAX_B, MIN_B
      REAL*8, PARAMETER    :: SEC_PER_YR = 365.35d0 * 86400d0  
      CHARACTER(LEN=225)   :: FILENAME 
      INTEGER              :: NYMD

      !=================================================================
      ! POP_READYR begins here!
      !=================================================================

      !IF (LPAH) THEN
      
      ! POLYCYCLIC AROMATIC HYDROCARBONS (PAHS):
      ! PAH emissions are for the year 2004
      ! Each PAH congener is emitted individually and contained in separate
      ! files
 
      ! Filename for congener you wish to model:
      !FILENAME = TRIM( DATA_DIR_1x1 )       // 
!     &           'PAHs_2004/PHE_EM_4x5.bpch' 
!      FILENAME = '/net/fs03/d0/geosdata/data/GEOS_4x5/PAHs_2004/' //
!     &           '/1x1/updated060911/PYR_EM_1x1.bpch'
      FILENAME = POP_EMISDIR
      
      ! Timestamp for emissions
      ! All PAH emissions are for the year 2004
      XTAU = GET_TAU0( 1, 1, 2004) 

      ! Echo info
      WRITE( 6, 100 )
100        FORMAT( '     - POPS_READYR: Reading ', a )

       ! Read data in [Mg/yr]
       CALL READ_BPCH2( FILENAME, 'PG-SRCE', 1, 
     &           XTAU,      I1x1,     J1x1,    
     &           1,         ARRAY,   QUIET=.FALSE. )

       ! Cast to REAL*8 and resize       
       CALL DO_REGRID_1x1( 'kg', ARRAY, POP_TOT_EM )  

       ! Convert from [Mg/yr] to [kg/s]
       POP_TOT_EM = POP_TOT_EM * 1000d0 / SEC_PER_YR
      
      ! Return to calling program
      END SUBROUTINE POPS_READYR

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  POPS_READSURFCONC
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE POPS_READSURFCONC
!! 
! !USES:
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : READ_BPCH2, GET_TAU0
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1
      USE TRANSFER_MOD,      ONLY : TRANSFER_2D
      USE TIME_MOD,          ONLY : EXPAND_DATE
    
! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  3 September 2012 - CL Friedman - Initial Version
!
! !REMARKS:
! (1) 
!
!EOP
!******************************************************************************
!Comment header
!  Subroutine POPS_READSURFCONC reads archived surface concentrations
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) 
!
!  Local variables:
!  ============================================================================
!  (1 ) 
!     
!  NOTES:
!  (1 ) Based on MERCURY_READYR in mercury_mod.f
!  
!  REFS:
!  (1 ) 

!------------------------------------------------------------------------------
!BOC
      USE CMN_SIZE_MOD             ! Size parameters ! hma 7Mar2013
!#     include "CMN_SIZE"       ! Size parameters
#     include "define.h"

      ! Local variables
      REAL*4               :: ARRAY(IGLOB,JGLOB,1) 
      REAL*8               :: XTAU, MAX_A, MIN_A
      REAL*8               :: MAX_B, MIN_B
      REAL*8, PARAMETER    :: SEC_PER_YR = 365.35d0 * 86400d0  
      CHARACTER(LEN=225)   :: FILENAME 
      INTEGER              :: NYMD

      !=================================================================
      ! POP_READSURFCONC begins here!
      !=================================================================
 
      !IF (LPAH) THEN

      ! For PAHs:
      ! Read in file which contains 1-yr worth of deposition based on
      ! annual average deposition rates [kg]
      ! Filename for PAHs:
!Default: GEOS 5 @ 4x5
!      FILENAME = '/net/fs03/d0/geosdata/data/GEOS_4x5/PAHs_2004/' //
!     &           '4x5/soil_conc/pyr_soil_conc_ann.bpch'

#if defined ( GEOS_5 ) && defined ( GRID4x5 )
      FILENAME = '/net/fs03/d0/geosdata/data/GEOS_4x5/PAHs_2004' //
     &           '/4x5/soil_conc/pyr_soil_conc_ann.bpch'

#elif defined ( GCAP ) && defined ( GRID4X5 )
      FILENAME = '/net/fs03/d0/geosdata/data/GCAP_4x5/PAHs_2004' //
     &           '/soil_conc/pyr_soil_conc_ann_GCAP.bpch'

#endif
     
         ! Timestamp for emissions
         ! All PAH soil conc are for the year 2009
         XTAU = GET_TAU0( 1, 1, 2009) 

         ! Echo info
         WRITE( 6, 100 )
100        FORMAT( '     - POPS_READSURFCONC: Reading ', a )

          PRINT*,'FILENAME = ', FILENAME
          ! Read data in [kg in surface layer]
           CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1, 
     &           XTAU,      IGLOB,     JGLOB,    
     &           1,         ARRAY,   QUIET=.FALSE. )
          
           ! Cast to REAL*8 and resize       
           CALL TRANSFER_2D( ARRAY(:,:,1), POP_SURF )

      
      ! Return to calling program
      END SUBROUTINE POPS_READSURFCONC

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  POPS_READ_FOC
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE POPS_READ_FOC
!! 
! !USES:
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : READ_BPCH2, GET_TAU0
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD,    ONLY : DO_REGRID_1x1
      USE TIME_MOD,          ONLY : EXPAND_DATE
    
! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  3 September 2012 - CL Friedman - Initial Version
!
! !REMARKS:
! (1) 
!
!EOP
!******************************************************************************
!Comment header
!  Subroutine POP_READ_FOC reads the mass of carbon in soils calculated with CASA
!  via the GTMM (files provided by Bess Sturges Corbitt)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) 
!
!  Local variables:
!  ============================================================================
!  (1 ) 
!     
!  NOTES:
!  (1 ) Based on MERCURY_READYR in mercury_mod.f
!  
!  REFS:
!  (1 ) 

!------------------------------------------------------------------------------
!BOC
      USE CMN_SIZE_MOD             ! Size parameters ! hma 7Mar2013
!#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      REAL*4               :: ARRAY(I1x1,J1x1,1) 
      REAL*8               :: XTAU, MAX_A, MIN_A
      REAL*8               :: MAX_B, MIN_B!, F_OC_SOIL
      REAL*8, PARAMETER    :: SEC_PER_YR = 365.35d0 * 86400d0  
      CHARACTER(LEN=255)   :: FILENAME 
      INTEGER              :: NYMD, I, J

      !=================================================================
      ! POP_READ_FOC begins here!
      !=================================================================

      ! Read in 1x1 file which contains sum of all soil carbon pools calculated
      ! by CASA/GTTM in g/m2
      ! Filename for PAHs
         FILENAME = '/net/fs03/d0/geosdata/data/GEOS_4x5/PAHs_2004' //
     &              '/1x1/soil_carbon/from_GTMM_09102012.bpch'
!      FILENAME = '/home/hamos/GC/input_files/POPs/PAHs_2004_GEOS' //
!     &           '/1x1/soil_carbon/from_GTMM_09102012.bpch'

         ! Timestamp for emissions
         ! All PAH soil conc are for the generic year 1985
         XTAU = GET_TAU0( 1, 1, 1985) 

         ! Echo info
         WRITE( 6, 100 )
100        FORMAT( '     - POPS_READ_FOC: Reading ', a )

         ! Read data in [g/m2 in top 30 cm of surface soil, according to Potter et al. 1993]
         ! NaNs and zeros have been replaced with SMALLNUM already 
         CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1, 
     &           XTAU,      I1x1,     J1x1,    
     &           1,         ARRAY,   QUIET=.FALSE. )

          ! Cast to REAL*8 and resize (does DO_REGRID_1x1 regrid 2-D arrays? -Yes!)      
          CALL DO_REGRID_1x1( 'g soil C/m2', ARRAY, F_OC_SOIL ) 

      DO J=1, JJPAR
      DO I=1, IIPAR

         ! Assume most of carbon mass extends to 5 cm and calculate concentration in g/g
         ! For now, assume a mean soil bulk density of 1300 kg/m3 similar to McLachlan 2002
         ! to calculate a dry weight fraction
         F_OC_SOIL(I,J) = F_OC_SOIL(I,J) / 30d-2 / 13d2 / 1d3


      ENDDO
      ENDDO

      
      ! Return to calling program
      END SUBROUTINE POPS_READ_FOC

!EOC
!------------------------------------------------------------------------------

      FUNCTION GET_NO3( I, J, L ) RESULT( NO3_MOLEC_CM3 ) 
!
!******************************************************************************
!  Function GET_NO3 returns monthly mean NO3 from archives. For offline runs, the 
!  conc of NO3 is set to zero during the day. (rjp, bmy, 12/16/02, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) Now references ERROR_STOP from "error_mod.f".  We also assume that
!        SETTRACE has been called to define IDNO3.  Now also set NO3 to 
!        zero during the day. (rjp, bmy, 12/16/02)
!  (2 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,     ONLY : CSPEC, JLOP
      USE DAO_MOD,        ONLY : AD,    SUNCOS
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE GLOBAL_NO3_MOD, ONLY : NO3
      USE TRACERID_MOD,   ONLY : IDNO3

      ! Replace include statements (H. Amos, 7 Mar 2013)
      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_MOD        ! NSRCX 
!#     include "CMN_SIZE"  ! Size parameters
!#     include "CMN"       ! NSRCX

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP
      REAL*8              :: NO3_MOLEC_CM3
      REAL*8,  PARAMETER  :: XNUMOL_NO3 = 6.022d23 / 62d-3
 
      ! External functions
      REAL*8,  EXTERNAL   :: BOXVL

      !=================================================================
      ! GET_NO3 begins here!
      !=================================================================
c$$$      IF ( ITS_A_FULLCHEM_SIM() ) THEN
c$$$
c$$$         !----------------------
c$$$         ! Fullchem simulation
c$$$         !----------------------
c$$$
c$$$         ! 1-D SMVGEAR grid box index
c$$$         JLOOP = JLOP(I,J,L)
c$$$
c$$$         ! Take NO3 from the SMVGEAR array CSPEC
c$$$         ! NO3 is defined only in the troposphere
c$$$         IF ( JLOOP > 0 ) THEN
c$$$            NO3_MOLEC_CM3 = CSPEC(JLOOP,IDNO3)
c$$$         ELSE
c$$$            NO3_MOLEC_CM3 = 0d0
c$$$         ENDIF
c$$$
c$$$      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

         !==============================================================  
         ! Offline simulation: Read monthly mean GEOS-CHEM NO3 fields
         ! in [v/v].  Convert these to [molec/cm3] as follows:
         !
         !  vol NO3   moles NO3    kg air     kg NO3/mole NO3
         !  ------- = --------- * -------- * ---------------- =  kg NO3 
         !  vol air   moles air      1        kg air/mole air
         !
         ! And then we convert [kg NO3] to [molec NO3/cm3] by:
         !  
         !  kg NO3   molec NO3   mole NO3     1     molec NO3
         !  ------ * --------- * -------- * ----- = --------- 
         !     1     mole NO3     kg NO3     cm3       cm3
         !          ^                    ^
         !          |____________________|  
         !            this is XNUMOL_NO3
         !
         ! If at nighttime, use the monthly mean NO3 concentration from
         ! the NO3 array of "global_no3_mod.f".  If during the daytime,
         ! set the NO3 concentration to zero.  We don't have to relax to 
         ! the monthly mean concentration every 3 hours (as for HNO3) 
         ! since NO3 has a very short lifetime. (rjp, bmy, 12/16/02) 
         !==============================================================

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! Test if daylight
         IF ( SUNCOS(JLOOP) > 0d0 ) THEN

            ! NO3 goes to zero during the day
            NO3_MOLEC_CM3 = 0d0
              
         ELSE

            ! At night: Get NO3 [v/v] and convert it to [kg]
            NO3_MOLEC_CM3 = NO3(I,J,L) * AD(I,J,L) * ( 62d0/28.97d0 ) 
!            WRITE (6,*) 'NO3 v/v =', NO3(I,J,L)
!            WRITE (6,*) 'NO3 kg =', NO3_MOLEC_CM3

            ! Convert NO3 from [kg] to [molec/cm3]
            NO3_MOLEC_CM3 = NO3_MOLEC_CM3 * XNUMOL_NO3 / BOXVL(I,J,L)
!            WRITE (6,*) 'XNUMOL_NO3 =', XNUMOL_NO3
!            WRITE (6,*) 'BOXVL =', BOXVL(I,J,L)
!            WRITE (6,*) 'NO3 molec/cm3 =', NO3_MOLEC_CM3
                  
         ENDIF
            
         ! Make sure NO3 is not negative
         NO3_MOLEC_CM3  = MAX( NO3_MOLEC_CM3, 0d0 )

c$$$      ELSE
c$$$
c$$$         !----------------------
c$$$         ! Invalid sim type!
c$$$         !----------------------       
c$$$         CALL ERROR_STOP( 'Invalid Simulation Type!',
c$$$     &                    'GET_NO3 ("carbon_mod.f")' )
c$$$
c$$$      ENDIF

      ! Return to calling program
      END FUNCTION GET_NO3

!------------------------------------------------------------------------------

      FUNCTION GET_O3( I, J, L ) RESULT( O3_MOLEC_CM3 )
!
!******************************************************************************
!  Function GET_O3 returns monthly mean O3 for offline sulfate aerosol
!  simulations. (bmy, 12/16/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L   (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDO3. (bmy, 12/16/02)
!  (2 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : AD
      USE GLOBAL_O3_MOD, ONLY : O3

      ! Replace include statement (H. Amos, 7Mar2013)
      USE CMN_SIZE_MOD ! Size parameters
!#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      REAL*8              :: O3_MOLEC_CM3

      ! External functions
      REAL*8, EXTERNAL    :: BOXVL
      
      !=================================================================
      ! GET_O3 begins here!
      !=================================================================

      ! Get ozone [v/v] for this gridbox & month
      ! and convert to [molec/cm3] (eck, 12/2/04)
      O3_MOLEC_CM3 = O3(I,J,L) * ( 6.022d23 / 28.97d-3 ) * 
     &               AD(I,J,L)  /  BOXVL(I,J,L)

      ! Return to calling program
      END FUNCTION GET_O3

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GET_OH
!
! !DESCRIPTION: Function GET_OH returns monthly mean OH and imposes a diurnal
! variation. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_OH( I, J, L ) RESULT( OH_MOLEC_CM3 )

      ! References to F90 modules
      USE DAO_MOD,       ONLY : SUNCOS 
      USE GLOBAL_OH_MOD, ONLY : OH
      USE TIME_MOD,      ONLY : GET_TS_CHEM

      ! Replace include statement (H. Amos, 7Mar2013)
      USE CMN_SIZE_MOD ! Size parameters
!#     include "CMN_SIZE"  ! Size parameters

! !INPUT PARAMETERS: 

      INTEGER, INTENT(IN) :: I, J, L
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  03 February 2011 - CL Friedman - Initial Version
!
! !REMARKS:
! Copied GET_OH function from mercury_mod.f - CLF
!
!EOP
!------------------------------------------------------------------------------
!BOC

       ! Local variables
       INTEGER       :: JLOOP
       REAL*8        :: OH_MOLEC_CM3
    
       !=================================================================
       ! GET_OH begins here!
       !=================================================================

       ! 1-D grid box index for SUNCOS
       JLOOP = ( (J-1) * IIPAR ) + I

       ! Test for sunlight...
       IF ( SUNCOS(JLOOP) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN

         ! Impose a diurnal variation on OH during the day
         OH_MOLEC_CM3 = OH(I,J,L)                      *           
     &                  ( SUNCOS(JLOOP) / TCOSZ(I,J) ) *
     &                  ( 1440d0        / GET_TS_CHEM() )

         ! Make sure OH is not negative
         OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0d0 )
               
       ELSE

         ! At night, OH goes to zero
         OH_MOLEC_CM3 = 0d0

       ENDIF

       ! Return to calling program
       END FUNCTION GET_OH
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GET_OC
!
! !DESCRIPTION: Function GET_OC returns monthly mean organic carbon 
! concentrations [kg/box]
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_OC( I, J, L) RESULT( C_OC )
!
! !INPUT PARAMETERS:

!     References to F90 modules
      USE GLOBAL_OC_MOD, ONLY : OC 

      INTEGER, INTENT(IN) :: I, J, L 
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  03 February 2011 - CL Friedman - Initial Version
!
! !REMARKS:
! Test
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Local variables
      REAL*8            :: C_OC

      !=================================================================
      ! GET_OC begins here!
      !=================================================================

      ! Get organic carbon concentration [kg/box] for this gridbox and month
      C_OC = OC(I,J,L)

      ! Make sure OC is not negative
      C_OC = MAX( C_OC, 0d0 )

      ! Return to calling program
      END FUNCTION GET_OC
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GET_BC
!
! !DESCRIPTION: Function GET_BC returns monthly mean black carbon concentrations
! [kg/box]
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_BC( I, J, L) RESULT( C_BC )

!
! !INPUT PARAMETERS: 
!
      ! References to F90 modules
      USE GLOBAL_BC_MOD, ONLY : BC
   
      INTEGER, INTENT(IN) :: I, J, L   
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  03 February 2011 - CL Friedman - Initial Version
!
! !REMARKS:
! Test
!
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Local variables
      REAL*8      :: C_BC    
    
      !=================================================================
      ! GET_BC begins here!
      !=================================================================

      ! Get black carbon concentration [kg/box] for this gridbox and month
      C_BC = BC(I,J,L)

      ! Return to calling program

      END FUNCTION GET_BC
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  OHNO3TIME
!
! !DESCRIPTION: Subroutine OHNO3TIME computes the sum of cosine of the solar zenith
!  angle over a 24 hour day, as well as the total length of daylight. 
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 12/8/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OHNO3TIME
!! 
! !USES:
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC
      USE TIME_MOD, ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT

      ! Replace include statements (H. Amos, 7Mar2013)
      USE CMN_SIZE_MOD ! Size parameters
      USE CMN_GCTM_MOD ! Physical constants
!#     include "CMN_SIZE"  ! Size parameters
!#     include "CMN_GCTM"  ! Physical constants

! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version for POPS_MOD
!
! !REMARKS:
!  (1 ) Copy code from COSSZA directly for now, so that we don't get NaN
!        values.  Figure this out later (rjp, bmy, 1/10/03)
!  (2 ) Now replace XMID(I) with routine GET_XMID from "grid_mod.f".  
!        Now replace RLAT(J) with routine GET_YMID_R from "grid_mod.f". 
!        Removed NTIME, NHMSb from the arg list.  Now use GET_NHMSb,
!        GET_ELAPSED_SEC, GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT from 
!        "time_mod.f". (bmy, 3/27/03)
!  (3 ) Now store the peak SUNCOS value for each surface grid box (I,J) in 
!        the COSZM array. (rjp, bmy, 3/30/04)
!  (4 ) Also added parallel loop over grid boxes (eck, bmy, 12/8/04)
!  (5 ) copied from mercury_mod by eck (9/20/10)
!******************************************************************************
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Local variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I, IJLOOP, J, L, N, NT, NDYSTEP
      REAL*8              :: A0, A1, A2, A3, B1, B2, B3
      REAL*8              :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
      REAL*8              :: SUNTMP(MAXIJ)
      
      !=================================================================
      ! OHNO3TIME begins here!
      !=================================================================

      !  Solar declination angle (low precision formula, good enough for us):
      A0 = 0.006918
      A1 = 0.399912
      A2 = 0.006758
      A3 = 0.002697
      B1 = 0.070257
      B2 = 0.000907
      B3 = 0.000148
      R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

      DEC = A0 - A1*cos(  R) + B1*sin(  R)
     &         - A2*cos(2*R) + B2*sin(2*R)
     &         - A3*cos(3*R) + B3*sin(3*R)

      LHR0 = int(float( GET_NHMSb() )/10000.)

      ! Only do the following at the start of a new day
      IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN 
      
         ! Zero arrays
         TTDAY(:,:) = 0d0
         TCOSZ(:,:) = 0d0
         COSZM(:,:) = 0d0

         ! NDYSTEP is # of chemistry time steps in this day
         NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 60 / GET_TS_CHEM()         

         ! NT is the elapsed time [s] since the beginning of the run
         NT = GET_ELAPSED_SEC()

         ! Loop forward through NDYSTEP "fake" timesteps for this day 
         DO N = 1, NDYSTEP
            
            ! Zero SUNTMP array
            SUNTMP(:) = 0d0

            ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YMID_R, IJLOOP, TIMLOC, AHR )
            DO J = 1, JJPAR

c$$$               ! Grid box latitude center [radians]
c$$$               !YMID_R = GET_YMID_R( J )
c$$$               YMID_R = GET_YMID_R( I, J, 1 ) ! hma, 8mar2013

            DO I = 1, IIPAR

               ! Grid box latitude center [radians]
               ! (moved inside I,J loops - H. Amos, 10Apr2013)
               YMID_R = GET_YMID_R( I, J, 1)


               ! Increment IJLOOP
               IJLOOP = ( (J-1) * IIPAR ) + I
               !TIMLOC = real(LHR0) + real(NT)/3600.0 + GET_XMID(I)/15.0
               TIMLOC = real(LHR0) + real(NT)/3600.0 +
     &                  GET_XMID( I, J, 1 )/15.0   ! bug fix (hma, 8mar2013)
         
               DO WHILE (TIMLOC .lt. 0)
                  TIMLOC = TIMLOC + 24.0
               ENDDO

               DO WHILE (TIMLOC .gt. 24.0)
                  TIMLOC = TIMLOC - 24.0
               ENDDO

               AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

            !===========================================================
            ! The cosine of the solar zenith angle (SZA) is given by:
            !     
            !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR) 
            !                   
            ! where LAT = the latitude angle, 
            !       DEC = the solar declination angle,  
            !       AHR = the hour angle, all in radians. 
            !
            ! If SUNCOS < 0, then the sun is below the horizon, and 
            ! therefore does not contribute to any solar heating.  
            !===========================================================

               ! Compute Cos(SZA)
               SUNTMP(IJLOOP) = sin(YMID_R) * sin(DEC) +
     &                          cos(YMID_R) * cos(DEC) * cos(AHR)

               ! TCOSZ is the sum of SUNTMP at location (I,J)
               ! Do not include negative values of SUNTMP
               TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(IJLOOP), 0d0 )

               ! COSZM is the peak value of SUMTMP during a day at (I,J)
               ! (rjp, bmy, 3/30/04)
               COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(IJLOOP) )

               ! TTDAY is the total daylight time at location (I,J)
               IF ( SUNTMP(IJLOOP) > 0d0 ) THEN
                  TTDAY(I,J) = TTDAY(I,J) + DBLE( GET_TS_CHEM() )
               ENDIF
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Increment elapsed time [sec]
            NT = NT + ( GET_TS_CHEM() * 60 )             
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Return to calling program
      END SUBROUTINE OHNO3TIME

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  INIT_POPS
!
! !DESCRIPTION: Subroutine INIT_POPS allocates and zeroes all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_POPS(POP_XMW, POP_KOA, POP_KBC, POP_K_POPG_OH,
     &                         POP_K_POPP_O3A, POP_K_POPP_O3B, 
     &                         POP_HSTAR, POP_DEL_H, POP_DEL_Hw )
!
      ! References to F90 modules
      USE DRYDEP_MOD,   ONLY : DEPNAME,   NUMDEP
      USE ERROR_MOD,    ONLY : ALLOC_ERR, ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT,    LDRYD,     LNLPBL
      USE TRACER_MOD,   ONLY : N_TRACERS
      USE PBL_MIX_MOD,  ONLY : GET_PBL_MAX_L
      ! Call to INIT_LAI obsolete (H. Amos, 7Mar2013)
!      USE LAI_MOD,      ONLY : INIT_LAI
      USE GET_POPSINFO_MOD, ONLY : INIT_POP_PARAMS


      ! Replace include statements (H. Amos, 7Mar2013)
      USE CMN_SIZE_MOD ! Size parameters
      USE CMN_DIAG_MOD ! ND44
!#     include "CMN_SIZE"     ! Size parameters
!#     include "CMN_DIAG"     ! ND44

! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on INIT_MERCURY from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Local variables
      LOGICAL, SAVE         :: IS_INIT = .FALSE. 
      INTEGER               :: AS, N!, PBL_MAX
      REAL*8                :: MAX_A, MIN_A


      REAL*8                :: POP_XMW, POP_KOA, POP_KBC, POP_K_POPG_OH
      REAL*8                :: POP_K_POPP_O3A, POP_K_POPP_O3B
      REAL*8                :: POP_HSTAR, POP_DEL_H, POP_DEL_Hw
      CHARACTER             :: POP_TYPE


      !=================================================================
      ! INIT_POPS begins here!
      !=================================================================

      ! Maximum extent of the PBL
      !PBL_MAX = GET_PBL_MAX_L()

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

      CALL INIT_POP_PARAMS( POP_XMW, POP_KOA, POP_KBC, POP_K_POPG_OH,
     &                         POP_K_POPP_O3A, POP_K_POPP_O3B, 
     &                         POP_HSTAR, POP_DEL_H, POP_DEL_Hw )

      !=================================================================
      ! Allocate and initialize arrays
      !=================================================================
      ALLOCATE( COSZM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSZM' )
      COSZM = 0d0

      ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
      TCOSZ = 0d0

      ALLOCATE( TTDAY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TTDAY' )
      TTDAY = 0d0

      ALLOCATE( C_OC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'C_OC' )
      C_OC = 0d0

      ALLOCATE( C_BC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'C_BC' )
      C_BC = 0d0

      ALLOCATE( SUM_OC_EM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUM_OC_EM' )
      SUM_OC_EM = 0d0

      ALLOCATE( SUM_BC_EM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUM_BC_EM' )
      SUM_BC_EM = 0d0

      ALLOCATE( SUM_G_EM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUM_G_EM' )
      SUM_G_EM = 0d0

      ALLOCATE( SUM_OF_ALL( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUM_OF_ALL' )
      SUM_OF_ALL = 0d0

      ALLOCATE( EPOP_G( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_G' )
      EPOP_G = 0d0

      ALLOCATE( EPOP_OC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_OC' )
      EPOP_OC = 0d0

      ALLOCATE( EPOP_BC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_BC' )
      EPOP_BC = 0d0

      ALLOCATE( EPOP_P_TOT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_P_TOT' )
      EPOP_P_TOT = 0d0

      ALLOCATE( POP_TOT_EM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'POP_TOT_EM' )
      POP_TOT_EM = 0d0

      ALLOCATE( POP_SURF( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'POP_SURF' )
      POP_SURF = 0d0

      ALLOCATE( EPOP_VEG( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_VEG' )
      EPOP_VEG = 0d0

      ALLOCATE( EPOP_LAKE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_LAKE' )
      EPOP_LAKE = 0d0

      ALLOCATE( EPOP_SOIL( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_SOIL' )
      EPOP_SOIL = 0d0

      ALLOCATE( EPOP_OCEAN( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_OCEAN' )
      EPOP_OCEAN = 0d0

      ALLOCATE( EPOP_SNOW( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_SNOW' )
      EPOP_SNOW = 0d0

      ALLOCATE( F_OC_SOIL( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'F_OC_SOIL' )
      F_OC_SOIL = 0d0

      ! Allocate ZERO_DVEL if we use non-local PBL mixing or
      ! if drydep is turned off 
      IF ( LNLPBL .OR. (.not. LDRYD) ) THEN
         ALLOCATE( ZERO_DVEL( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ZERO_DVEL' )
         ZERO_DVEL = 0d0
      ENDIF

      ! Initialze LAI arrays
!      CALL INIT_LAI

      !=================================================================
      ! Done
      !=================================================================

      ! Reset IS_INIT, since we have already allocated arrays
      IS_INIT = .TRUE.
      
      ! Return to calling program
      END SUBROUTINE INIT_POPS

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CLEANUP_POPS
!
! !DESCRIPTION: Subroutine CLEANUP_POPS deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_POPS
!

! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on INIT_MERCURY from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!BOC

      IF ( ALLOCATED( COSZM     ) ) DEALLOCATE( COSZM   )     
      IF ( ALLOCATED( TCOSZ     ) ) DEALLOCATE( TCOSZ   )
      IF ( ALLOCATED( TTDAY     ) ) DEALLOCATE( TTDAY   )
      IF ( ALLOCATED( ZERO_DVEL ) ) DEALLOCATE( ZERO_DVEL )
      IF ( ALLOCATED( EPOP_G    ) ) DEALLOCATE( EPOP_G   )
      IF ( ALLOCATED( EPOP_OC   ) ) DEALLOCATE( EPOP_OC  )
      IF ( ALLOCATED( EPOP_BC   ) ) DEALLOCATE( EPOP_BC )
      IF ( ALLOCATED( EPOP_P_TOT) ) DEALLOCATE( EPOP_P_TOT )
      IF ( ALLOCATED( POP_TOT_EM) ) DEALLOCATE( POP_TOT_EM )
      IF ( ALLOCATED( POP_SURF  ) ) DEALLOCATE( POP_SURF )
      IF ( ALLOCATED( EPOP_VEG  ) ) DEALLOCATE( EPOP_VEG )
      IF ( ALLOCATED( EPOP_LAKE ) ) DEALLOCATE( EPOP_LAKE )
      IF ( ALLOCATED( EPOP_SOIL ) ) DEALLOCATE( EPOP_SOIL )
      IF ( ALLOCATED( EPOP_SNOW ) ) DEALLOCATE( EPOP_SNOW )
      IF ( ALLOCATED( EPOP_OCEAN) ) DEALLOCATE( EPOP_OCEAN )
      IF ( ALLOCATED( F_OC_SOIL ) ) DEALLOCATE( F_OC_SOIL )
      IF ( ALLOCATED( C_OC      ) ) DEALLOCATE( C_OC )
      IF ( ALLOCATED( C_BC      ) ) DEALLOCATE( C_BC )
      IF ( ALLOCATED( SUM_OC_EM  ) ) DEALLOCATE( SUM_OC_EM )
      IF ( ALLOCATED( SUM_BC_EM  ) ) DEALLOCATE( SUM_BC_EM )
      IF ( ALLOCATED( SUM_G_EM   ) ) DEALLOCATE( SUM_G_EM ) 
      IF ( ALLOCATED( SUM_OF_ALL ) ) DEALLOCATE( SUM_OF_ALL ) 

      END SUBROUTINE CLEANUP_POPS
!EOC
!------------------------------------------------------------------------------
      END MODULE POPS_MOD

