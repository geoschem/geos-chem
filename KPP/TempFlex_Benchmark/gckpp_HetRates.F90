!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gckpp_HetRates
!
! !DESCRIPTION: FlexChem module for heterogeneous chemistry, via KPP.
!\\
!\\
! !INTERFACE:
!
MODULE GCKPP_HETRATES
!
! !USES:
!
  USE CMN_FJX_MOD,        ONLY : NDUST
  USE CMN_FJX_MOD,        ONLY : NAER
  USE CMN_SIZE_MOD,       ONLY : LLSTRAT
  USE COMODE_LOOP_MOD,    ONLY : CONSVAP
  USE ERROR_MOD,          ONLY : ERROR_STOP
  USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
  USE ERROR_MOD,          ONLY : IS_SAFE_DIV
  USE gckpp_Precision
  USE gckpp_Parameters
  USE gckpp_Global,       ONLY : HET
  USE GIGC_State_Chm_Mod, ONLY : ChmState
  USE GIGC_State_Chm_Mod, ONLY : Get_Indx
  USE GIGC_State_Met_Mod, ONLY : MetState
  USE GIGC_Input_Opt_Mod, ONLY : OptInput
  USE Precision_Mod,      ONLY : fp

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: SET_HET
!
! !PRIVATE MEMBER FUNCTIONS:
!
  ! These functions are used for all mechanisms
  PRIVATE :: HETNO3
  PRIVATE :: HETNO2
  PRIVATE :: HETHO2
  PRIVATE :: HETHBr
  PRIVATE :: HETN2O5
  PRIVATE :: HETBrNO3
  PRIVATE :: HETHOBr
  PRIVATE :: HETHOBr_ice
  PRIVATE :: HETHBr_ice
  PRIVATE :: N2O5
  PRIVATE :: HO2
  PRIVATE :: CLD1K_BrNO3
  PRIVATE :: FCRO2HO2
  PRIVATE :: FYHORO
  PRIVATE :: FYRNO3
  PRIVATE :: ARSL1K

#if defined( UCX )
  ! These functions are only used for UCX-based mechanisms
  PRIVATE :: HETClNO3_PSC1
  PRIVATE :: HETClNO3_PSC2
  PRIVATE :: HETClNO3_PSC3
  PRIVATE :: HETBrNO3_PSC
  PRIVATE :: HETHOCl_PSC1
  PRIVATE :: HETHOCl_PSC2
  PRIVATE :: HETHOBr_PSC
  PRIVATE :: HETN2O5_PSC
#endif
!
! !PRIVATE DATA MEMBERS:
!
  !%%% NOTE: SOME OF THESE VARIABLE DEFINITIONS MAY BE BETTER IMPLEMENTED %%%
  !%%% AS LOCAL VARIABLES WITHIN EACH OF THE FUNCTIONS.  DO THIS LATER.   %%%

  ! Scalars
  INTEGER  :: II,JJ,LL
  INTEGER  :: N
  INTEGER  :: NAERO
  LOGICAL  :: NATSURFACE,   SAFEDIV
  LOGICAL  :: KII_KI,       PSCBOX,    STRATBOX
  REAL(fp) :: TEMPK,        RELHUM,    XSTKCF
  REAL(fp) :: VPRESH2O,     CONSEXP
  REAL(fp) :: TRC_NIT,      TRC_SO4,   XNM_SO4,   XNM_NIT
  REAL(fp) :: GAMMA_HO2,    XTEMP,     XDENA,     ADJUSTEDRATE
  REAL(fp) :: CLD_BRNO3_RC, KI_HBR,    KI_HOBr,   QLIQ,        QICE
  REAL(fp) :: DUMMY
  REAL(fp) :: hobr_rtemp,   hbr_rtemp
  REAL(fp) :: SPC_HBr,      SPC_HOBr
#if defined( UCX )
  INTEGER  :: PSCIDX
  REAL(fp) :: SPC_N2O5,     SPC_H2O,   SPC_HCl
  REAL(fp) :: SPC_HOCl,     SPC_ClNO3, SPC_BrNO3
  REAL(fp) :: EDUCTCONC,    LIMITCONC
#endif

  ! Arrays
  REAL(fp) :: SCF2(3)
  REAL(fp) :: XAREA(25)
  REAL(fp) :: XRADI(25)
  REAL(fp) :: KHETI_SLA(11)
#if defined( UCX )
  REAL(fp) :: PSCEDUCTCONC(11,2)
#endif

!$OMP THREADPRIVATE( KII_KI,    NAERO,        N                        )
!$OMP THREADPRIVATE( RELHUM,    CONSEXP,      VPRESH2O                 )
!$OMP THREADPRIVATE( PSCBOX,    STRATBOX,     NATSURFACE               )
!$OMP THREADPRIVATE( TRC_SO4,   TRC_NIT,      XNM_SO4,    XNM_NIT      )
!$OMP THREADPRIVATE( XAREA,     XRADI,        TEMPK,      XTEMP        )
!$OMP THREADPRIVATE( XDENA,     GAMMA_HO2,    QICE,       QLIQ         )
!$OMP THREADPRIVATE( DUMMY,     ki_hbr,       ki_hobr,    cld_brno3_rc )
!$OMP THREADPRIVATE( SAFEDIV,   hbr_rtemp,    hobr_rtemp               )
!$OMP THREADPRIVATE( XSTKCF,    ADJUSTEDRATE, SPC_HBr,    SPC_HOBr     )
#if defined( UCX )
!$OMP THREADPRIVATE( SPC_N2O5,  SPC_H2O,      SPC_HCl,    SPC_HOCl     )
!$OMP THREADPRIVATE( SPC_ClNO3, SPC_BrNO3,    KHETI_SLA,  PSCEDUCTCONC )
!$OMP THREADPRIVATE( PSCIDX,    EDUCTCONC,    LIMITCONC                ) 
#endif
!
! !DEFINED PARAMETERS:
!
  REAL(fp), PARAMETER :: PSCMINLIFE = 1.e-3_fp
!
! !REMARKS:
!  Need 
!  - TOTAREA (previously used for archiving N2O5 hydrolysis in the planeflight
!             diagnostic only)
!  - Air NUM. DENSITY
!  - TEMPERATURE
!  - Aerosol Surface Area
!  - Aerosol Type
!  - Gamma (XSTKCF; sticking factor)
!  - ARR
!  - Species num density (mcl cm-3)
!  - Continental PBL or no?
!  - In stratosphere or no?
!  - Reaction index (e.g. NK1HBr, NK2HBr)
!
!  According to S. Eastham, we should also include
!  cloud and ice area explicitly, in addition to
!  aerosol area
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - NOTE: SPC_HBR and SPC_HOBR are defined
!                              for trop-only mechanisms
!  29 Mar 2016 - R. Yantosca - Added ProTeX headers
!  29 Mar 2016 - R. Yantosca - Moved all the UCX-based functions to the
!                              end of the module, for clarity
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
! !IROUTINE: set_het
!
! !DESCRIPTION: Main heterogenous chemistry driver routine.  Sets up the
!  vector of heterogeneous chemistry rates for the KPP chemistry solver.
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE SET_HET( I, J, L, SC, SM, IO, SCF )
!
! !INPUT PARAMETERS: 
!
      INTEGER        :: I, J, L   ! Lon, lat, level indices
      TYPE(MetState) :: SM        ! Meteorology State object
      TYPE(OptInput) :: IO        ! Input Options object 

!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState) :: SC        ! Chemistry Sate object
      REAL(fp)       :: SCF(3)    ! Coefficients (Need help documenting this)
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: IND

      !====================================================================
      ! SET_HET begins here!
      !====================================================================

      ! For Debugging
      II = I
      JJ = J
      LL = L

      ! Divide by educt concentration
      KII_KI     = .false.

      ! Initialize logicals for UCX
      PSCBOX     = .false.
      STRATBOX   = .false.
      NATSURFACE = .false.

      NAERO = SC%nAero

#if defined( UCX )
      ! Copy sticking coefficients for PSC reactions on SLA
      KHETI_SLA(:) = SC%KHETI_SLA(I,J,L,:)
#endif

      ! Calculate RH. Not clear why the result of this calc is 
      ! slightly different than SM%RH
      RELHUM        = SM%AVGW(I,J,L) * SM%AIRNUMDEN(I,J,L)
      CONSEXP       = 17.2693882e+0_fp * (SM%T(I,J,L) - 273.16e+0_fp) / &
                      (SM%T(I,J,L) - 35.86e+0_fp)
      VPRESH2O      = CONSVAP * EXP(CONSEXP) / SM%T(I,J,L) 
      RELHUM = RELHUM / VPRESH2O 

      ! Get tracer concentrations [kg]
      IND = get_indx('SO4',IO%ID_TRACER,IO%TRACER_NAME)
      IF (IND .le. 0) THEN
         TRC_SO4    = 0.0e+0_fp
         XNM_SO4    = 0.0e+0_fp
      ELSE
         TRC_SO4    = SC%Tracers(I,J,L,IND)
         XNM_SO4    = IO%XNUMOL(IND)
      ENDIF

      IND = get_indx('NIT',IO%ID_TRACER,IO%TRACER_NAME)
      IF (IND .le. 0) THEN
         TRC_NIT    = 0.0e+0_fp
         XNM_NIT    = 0.0e+0_fp
      ELSE
         TRC_NIT    = SC%Tracers(I,J,L,IND)
         XNM_NIT    = IO%XNUMOL(IND)
      ENDIF

      IND = get_indx('HBr',SC%Spec_ID,SC%Spec_Name)
      IF (IND .le. 0) THEN
         SPC_HBr    = 0.0e+0_fp
      ELSE
         SPC_HBr    = SC%Species(I,J,L,IND)
      ENDIF

      IND = get_indx('HOBr',SC%Spec_ID,SC%Spec_Name)
      IF (IND .le. 0) THEN
         SPC_HOBr   = 0.0e+0_fp
      ELSE
         SPC_HOBr   = SC%Species(I,J,L,IND)
      ENDIF

#if defined( UCX )
      ! Get species concentrations [molec/cm3]
      IND = get_indx('N2O5',SC%Spec_ID,SC%Spec_Name)
      IF (IND .le. 0) THEN
         SPC_N2O5   = 0.0e+0_fp
      ELSE
         SPC_N2O5   = SC%Species(I,J,L,IND)
      ENDIF

      IND = get_indx('H2O',SC%Spec_ID,SC%Spec_Name)
      IF (IND .le. 0) THEN
         SPC_H2O    = 0.0e+0_fp
      ELSE
         SPC_H2O    = SC%Species(I,J,L,IND)
      ENDIF

      IND = get_indx('HCl',SC%Spec_ID,SC%Spec_Name)
      IF (IND .le. 0) THEN
         SPC_HCl    = 0.0e+0_fp
      ELSE
         SPC_HCl    = SC%Species(I,J,L,IND)
      ENDIF

      IND = get_indx('ClNO3',SC%Spec_ID,SC%Spec_Name)
      IF (IND .le. 0) THEN
         SPC_ClNO3  = 0.0e+0_fp
      ELSE
         SPC_ClNO3  = SC%Species(I,J,L,IND)
      ENDIF

      IND = get_indx('BrNO3',SC%Spec_ID,SC%Spec_Name)
      IF (IND .le. 0) THEN
         SPC_BrNO3  = 0.0e+0_fp
      ELSE
         SPC_BrNO3  = SC%Species(I,J,L,IND)
      ENDIF

      IND = get_indx('HOCl',SC%Spec_ID,SC%Spec_Name)
      IF (IND .le. 0) THEN
         SPC_HOCl   = 0.0e+0_fp
      ELSE
         SPC_HOCl   = SC%Species(I,J,L,IND)
      ENDIF

      ! Set PSC educt concentrations (SDE 04/24/13)
      PSCEDUCTCONC( 1,1) = SPC_N2O5
      PSCEDUCTCONC( 1,2) = SPC_H2O

      PSCEDUCTCONC( 2,1) = SPC_N2O5
      PSCEDUCTCONC( 2,2) = SPC_HCl

      PSCEDUCTCONC( 3,1) = SPC_ClNO3
      PSCEDUCTCONC( 3,2) = SPC_H2O

      PSCEDUCTCONC( 4,1) = SPC_ClNO3
      PSCEDUCTCONC( 4,2) = SPC_HCl

      PSCEDUCTCONC( 5,1) = SPC_ClNO3
      PSCEDUCTCONC( 5,2) = SPC_HBr

      PSCEDUCTCONC( 6,1) = SPC_BrNO3
      PSCEDUCTCONC( 6,2) = SPC_H2O

      PSCEDUCTCONC( 7,1) = SPC_BrNO3
      PSCEDUCTCONC( 7,2) = SPC_HCl

      PSCEDUCTCONC( 8,1) = SPC_HOCl
      PSCEDUCTCONC( 8,2) = SPC_HCl

      PSCEDUCTCONC( 9,1) = SPC_HOCl
      PSCEDUCTCONC( 9,2) = SPC_HBr

      PSCEDUCTCONC(10,1) = SPC_HOBr
      PSCEDUCTCONC(10,2) = SPC_HCl

      ! This is still pseudo-first-order - ignore
      PSCEDUCTCONC(11,1) = SPC_HOBr
      PSCEDUCTCONC(11,2) = SPC_HBr
#endif

      XAREA(1:SC%nAero) = SC%AeroArea(I,J,L,:)
      XRADI(1:SC%nAero) = SC%AeroRadi(I,J,L,:)

      TEMPK = SM%T(I,J,L)
      XTEMP = sqrt(SM%T(I,J,L))
      XDENA = SM%AIRNUMDEN(I,J,L)

      GAMMA_HO2 = IO%GAMMA_HO2

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )
            
      ! GEOS-5 / MERRA / GEOS-FP / MERRA-2 have QI and QL defined as 
      ! met fields so use these to define the QICE, QLIQ arrays. 
      QICE       = SM%QI(I,J,L)
      QLIQ       = SM%QL(I,J,L)
      
#else
      
      ! Otherwise, compute QLIQ as a function of temperature ...
      IF ( SM%T(I,J,L) .LE. 248e+0_fp ) THEN
         QLIQ  = 0e+0_fp
      ELSE IF ( SM%T(I,J,L) .GE. 268e+0_fp ) THEN
         QLIQ  = 1e-6_fp
      ELSE
         QLIQ  = 1e-6_fp * ( ( SM%T(I,J,L) - 248e+0_fp ) / 20e+0_fp)
      ENDIF
      
      ! ... and compute QICE from QLIQ (bmy, 9/24/12)
      QICE     = 1e-6_fp - QLIQ
      
#endif

#if defined( UCX )
      ! Check surface type of PSCs (SDE 04/17/13)
      CALL CHECK_NAT( I,  J,  L, NATSURFACE, PSCBOX, STRATBOX, &
                      IO, SM, SC )
#endif

      ! ----------------------------------------------
      !  Calculate rate for cloud heterogeneous
      !  chemistry (jpp, 2/28/2011)
      ! ----------------------------------------------
      IF (.not.PSCBOX) THEN
         cld_brno3_rc = CLD1K_BrNO3(I,J,L,XDENA,QLIQ, SM )
      END IF

      ! ----------------------------------------------
      !  Calculate rates for HOBr + HBr + ice --> Br2
      !  for cold and mixed clouds. (jpp, 6/16/2011)
      ! ----------------------------------------------
      IF (.not.PSCBOX) THEN
         DUMMY = 0.0e+0_fp
         CALL cldice_hbrhobr_rxn( I,J,L,XDENA,QICE,SPC_HBr,SPC_HOBr, &
              ki_hbr, ki_hobr, DUMMY, SM )
      ELSE
         ! For PSCs, het chem already accounted for in
         ! aerosol code <-- IS THIS STILL TRUE? (MSL)
         ki_hbr = 0e+0_fp
         ki_hobr = 0e+0_fp
      ENDIF

      HET = 0.

      !----------------------------------------------------------------
      ! Calculate and pass het rates to the KPP rate array
      !----------------------------------------------------------------
      HET(ind_HO2,1)   = HETHO2(        3.30E1_fp, 2E-1_fp)
      HET(ind_NO2,1)   = HETNO2(        4.60E1_fp, 1E-4_fp)
      HET(ind_NO3,1)   = HETNO3(        6.20E1_fp, 1E-1_fp)
      HET(ind_N2O5,1)  = HETN2O5(       1.08E2_fp, 1E-1_fp)
      HET(ind_BrNO3,1) = HETBrNO3(      1.42E2_fp, 3E-1_fp)
      HET(ind_HOBr,1)  = HETHOBr(       0.97E2_fp, 2E-1_fp)
      HET(ind_HBr,1)   = HETHBr(        0.81E2_fp, 2E-1_fp)
      HET(ind_HOBr,2)  = HETHOBr_ice(   0.97E2_fp, 1E-1_fp)
      HET(ind_HBr,2)   = HETHBr_ice(    0.81E2_fp, 1E-1_fp)
#if defined( UCX )
      HET(ind_N2O5,2)  = HETN2O5_PSC(   1.08E2_fp, 0E+0_fp)
      HET(ind_ClNO3,1) = HETClNO3_PSC1( 0.97E2_fp, 0E+0_fp)
      HET(ind_ClNO3,2) = HETClNO3_PSC2( 0.97E2_fp, 0E+0_fp)
      HET(ind_ClNO3,3) = HETClNO3_PSC3( 0.97E2_fp, 0E+0_fp)
      HET(ind_BrNO3,2) = HETBrNO3_PSC(  1.42E2_fp, 0E+0_fp)
      HET(ind_HOCl,1)  = HETHOCl_PSC1(  0.52E2_fp, 0E+0_fp)
      HET(ind_HOCl,2)  = HETHOCl_PSC2(  0.52E2_fp, 0E+0_fp)
      HET(ind_HOBr,3)  = HETHOBr_PSC(   0.97E2_fp, 0E+0_fp)
#endif

      !----------------------------------------------------------------
      ! Kludging the rates to be equal to one another to avoid having
      ! to keep setting equality in solver. (jpp, 5/10/2011)
      !----------------------------------------------------------------
      IF ( ( HET(ind_HBr,1) > 0 ) .and. ( HET(ind_HOBr,1) > 0 ) ) THEN

         ! select the min of the two rates
         hbr_rtemp  = HET(ind_HBr,1)  * SPC_HBr
         hobr_rtemp = HET(ind_HOBr,1) * SPC_HOBr

         ! if HBr rate is larger than HOBr rate
         IF ( hbr_rtemp > hobr_rtemp ) THEN

            SAFEDIV = IS_SAFE_DIV( HET(ind_HOBr,1) * SPC_HOBr, SPC_HBr )

            IF (SAFEDIV) THEN
               ! 2. if it is safe, then go ahead
               HET(ind_HBr,1) = HET(ind_HOBr,1) * SPC_HOBr / SPC_HBr
            ELSE
               ! if not, then set rates really small...
               ! b/c the largest contributor is very small.
               HET(ind_HBr,1)  = TINY(1.e+0_fp)
               HET(ind_HOBr,1) = TINY(1.e+0_fp)
            ENDIF

            ! if HOBr rate is larger than HBr rate
         ELSE

            ! 1. is it safe to divide?
            SAFEDIV = IS_SAFE_DIV( HET(ind_HBr,1) * SPC_HBr, SPC_HOBr )

            IF (SAFEDIV) THEN
               ! 2. if it is safe, then go ahead
               HET(ind_HOBr,1) = HET(ind_HBr,1) * SPC_HBr / SPC_HOBr
            ELSE
               ! if not, then set rates really small...
               ! b/c the largest contributor is very small.
               HET(ind_HBr,1)  = TINY(1.e+0_fp)
               HET(ind_HOBr,1) = TINY(1.e+0_fp)
            ENDIF
         ENDIF
      ENDIF

      !----------------------------------------------------------------
      ! SDE 05/30/13: Limit rates to prevent solver failure for PSC
      ! het. chem.
      !----------------------------------------------------------------
#if defined( UCX )
      DO PSCIDX=1,10

         ! Pseudo-first-order reactions - divide by number-conc
         ! of aerosol-phase educt to yield 2nd-order constant
         EDUCTCONC = PSCEDUCTCONC(PSCIDX,2)
         LIMITCONC = PSCEDUCTCONC(PSCIDX,1)

         ! Initialize adjusted rates
         IF     ( PSCIDX .eq. 1 ) THEN
            ! N2O5 + H2O
            ADJUSTEDRATE = HET(ind_N2O5,1)
         ELSEIF ( PSCIDX .eq. 2 ) THEN
            ! N2O5 + HCl
            ADJUSTEDRATE = HET(ind_N2O5,2)
         ELSEIF ( PSCIDX .eq. 3 ) THEN
            ! ClNO3 + H2O
            ADJUSTEDRATE = HET(ind_ClNO3,1)
         ELSEIF ( PSCIDX .eq. 4 ) THEN
            ! ClNO3 + HCl
            ADJUSTEDRATE = HET(ind_ClNO3,2)
         ELSEIF ( PSCIDX .eq. 5 ) THEN
            ! ClNO3 + HBr
            ADJUSTEDRATE = HET(ind_ClNO3,3)
         ELSEIF ( PSCIDX .eq. 6 ) THEN
            ! BrNO3 + H2O
            ADJUSTEDRATE = HET(ind_BrNO3,1)
         ELSEIF ( PSCIDX .eq. 7 ) THEN
            ! BrNO3 + HCl
            ADJUSTEDRATE = HET(ind_BrNO3,2)
         ELSEIF ( PSCIDX .eq. 8 ) THEN
            ! HOCl + HCl
            ADJUSTEDRATE = HET(ind_HOCl,1)
         ELSEIF ( PSCIDX .eq. 9 ) THEN
            ! HOCl + HBr
            ADJUSTEDRATE = HET(ind_HOCl,2)
         ELSEIF ( PSCIDX .eq. 10) THEN
            ! HOBr + HCl
            ADJUSTEDRATE = HET(ind_HOBr,3)
         ENDIF

         ! ---SAFETY-CHECK REACTION---
         ! Definition of 2nd order reaction rate:
         ! k[A][B] = -d[A]/dt = -d[B]/dt
         !
         ! However, here we are using a pseudo-first order
         ! reaction rate, ki, and assuming that [B] is
         ! abundant. To get k, we will therefore perform:
         ! k = ki/[B]
         !
         ! This will yield the following when solved:
         ! -d[A]/dt = ki[A] = -d[B]/dt
         !
         ! This has some problems, especially for small [B]!
         ! To get around this, we run the following checks:
         !
         ! 1. The lifetime of [A] is 1/ki. If this is below
         !    PSCMINLIFE, limit reaction rate to yield the
         !    specified lifetimedepletion (ki = 1/60)
         ! 2. The depletion time of [B] is [B]/(ki[A]). If
         !    this is below PSCMINLIFE, limit reaction rate
         !    (ki = [B]/(T*[A])
         ! 3. If [B] is < 100 molec/cm3, or ki/[B] yields
         !    a Nan, force k = 0.
         !
         ! If all these checks are passed, we set k = ki/[B].
         ! Rxn 11 is first-order - ignore
         IF (PSCIDX.eq.1) THEN
            ! Convert from 1st-order to 2nd-order
            SAFEDIV = IS_SAFE_DIV(EDUCTCONC,LIMITCONC)
            IF (SAFEDIV) THEN
               ! Temporarily store [B]/(T*[A])
               LIMITCONC = EDUCTCONC/(PSCMINLIFE*LIMITCONC)
               IF (ADJUSTEDRATE.gt.LIMITCONC) THEN
                  ADJUSTEDRATE = LIMITCONC
               ENDIF
            ELSE
               ADJUSTEDRATE = 0e+0_fp
            ENDIF
            SAFEDIV = IS_SAFE_DIV(ADJUSTEDRATE,EDUCTCONC)
            IF ((EDUCTCONC.gt.1.e+2_fp).and. (SAFEDIV)) THEN
               ADJUSTEDRATE = ADJUSTEDRATE/EDUCTCONC
            ELSE
               ADJUSTEDRATE = 0e+0_fp
            ENDIF
         ELSEIF (PSCIDX.ne.11) THEN
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
               ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
            ENDIF
            ! Convert from 1st-order to 2nd-order
            SAFEDIV = IS_SAFE_DIV(EDUCTCONC,LIMITCONC)
            IF (SAFEDIV) THEN
               ! Temporarily store [B]/(T*[A])
               LIMITCONC = EDUCTCONC/(PSCMINLIFE*LIMITCONC)
               IF (ADJUSTEDRATE.gt.LIMITCONC) THEN
                  ADJUSTEDRATE = LIMITCONC
               ENDIF
            ELSE
               ADJUSTEDRATE = 0e+0_fp
            ENDIF
            SAFEDIV = IS_SAFE_DIV(ADJUSTEDRATE,EDUCTCONC)
            IF ((EDUCTCONC.gt.1.e+2_fp).and. (SAFEDIV)) THEN
               ADJUSTEDRATE = ADJUSTEDRATE/EDUCTCONC
            ELSE
               ADJUSTEDRATE = 0e+0_fp
            ENDIF
         ENDIF

         ! Copy adjusted rates to HET
         IF     ( PSCIDX .eq. 1 ) THEN
            ! N2O5 + H2O
            HET(ind_N2O5,1) = ADJUSTEDRATE
         ELSEIF ( PSCIDX .eq. 2 ) THEN
            ! N2O5 + HCl
            HET(ind_N2O5,2) = ADJUSTEDRATE
         ELSEIF ( PSCIDX .eq. 3 ) THEN
            ! ClNO3 + H2O
            HET(ind_ClNO3,1) = ADJUSTEDRATE
         ELSEIF ( PSCIDX .eq. 4 ) THEN
            ! ClNO3 + HCl
            HET(ind_ClNO3,2) = ADJUSTEDRATE
         ELSEIF ( PSCIDX .eq. 5 ) THEN
            ! ClNO3 + HBr
            HET(ind_ClNO3,3) = ADJUSTEDRATE
         ELSEIF ( PSCIDX .eq. 6 ) THEN
            ! BrNO3 + H2O
            HET(ind_BrNO3,1) = ADJUSTEDRATE
         ELSEIF ( PSCIDX .eq. 7 ) THEN
            ! BrNO3 + HCl
            HET(ind_BrNO3,2) = ADJUSTEDRATE
         ELSEIF ( PSCIDX .eq. 8 ) THEN
            ! HOCl + HCl
            HET(ind_HOCl,1) = ADJUSTEDRATE
         ELSEIF ( PSCIDX .eq. 9 ) THEN
            ! HOCl + HBr
            HET(ind_HOCl,2) = ADJUSTEDRATE
         ELSEIF ( PSCIDX .eq. 10) THEN
            ! HOBr + HCl
            HET(ind_HOBr,3) = ADJUSTEDRATE
         ENDIF

      ENDDO
#endif

      SCF = SCF2

      RETURN

    END SUBROUTINE SET_HET
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hetno3
!
! !DESCRIPTION: Set the heterogenous chemistry rate for NO3.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETNO3( A, B ) RESULT( HET_NO3 )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_NO3
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Initialize
      HET_NO3      = 0.0_fp
      ADJUSTEDRATE = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAERO

         XSTKCF = B
         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF
         
         IF (KII_KI .and. N.gt.12) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
               ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
            ENDIF
         ENDIF
         
         ! Add to overall reaction rate
         HET_NO3 = HET_NO3 + ADJUSTEDRATE

      ENDDO
      
    END FUNCTION HETNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hetno2
!
! !DESCRIPTION: Set the heterogenous chemistry rate for NO2.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETNO2( A, B ) RESULT( HET_NO2 )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
! 
      REAL(fp)             :: HET_NO2
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Initialize
      HET_NO2      = 0.0_fp
      ADJUSTEDRATE = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAERO

         XSTKCF = B
         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF
         
         IF (KII_KI .and. N.gt.12) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
               ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
            ENDIF
         ENDIF
         
         ! Add to overall reaction rate
         HET_NO2 = HET_NO2 + ADJUSTEDRATE

      END DO
      
    END FUNCTION HETNO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hetho2
!
! !DESCRIPTION: Set the heterogenous chemistry rate for HO2.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHO2( A, B ) RESULT( HET_HO2 )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_HO2
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Initialize
      HET_HO2      = 0.0_fp
      ADJUSTEDRATE = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAERO

         IF (N.gt.12) THEN
            XSTKCF = TINY(1e+0_fp)
         ELSE
            XSTKCF = GAMMA_HO2
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF
         
         IF (KII_KI .and. N.gt.12) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
               ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
            ENDIF
         ENDIF
         
         ! Add to overall reaction rate
         HET_HO2 = HET_HO2 + ADJUSTEDRATE

      ENDDO

    END FUNCTION HETHO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hethbr
!
! !DESCRIPTION: Set the heterogeneous rate for HBr.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHBr( A, B ) RESULT( HET_HBr )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_HBr
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Initialize
      HET_HBr      = 0.0_fp
      ADJUSTEDRATE = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAERO

         ! jpp, 3/22/11: set the sticking coefficient to 
         !  ~0 for aerosol types we don't want reactions on
         !  for the HBr and HOBr surface reaction
         
         ! Only apply adjustment if at high altitude
         KII_KI = STRATBOX
         
         ! Select proper aerosol type
         IF ( (N == 8) .OR. (N == 11) .OR. (N == 12)) THEN
            ! sulfate, 2 modes of sea-salt
            XSTKCF = B
         ELSEIF ( N == 13 ) THEN
            XSTKCF = KHETI_SLA(11)
         ELSEIF ( N == 14 ) THEN
            XSTKCF = 0.1e+0_fp
         ELSE
            XSTKCF = 0e+0_fp
         ENDIF
         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         scf2(2) = xstkcf
         
         IF (KII_KI .and. N.gt.12) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
               ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_HBr = HET_HBr + ADJUSTEDRATE

      END DO

    END FUNCTION HETHBr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hetn2o5
!
! !DESCRIPTION: Set heterogenous chemistry rate for N2O5.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETN2O5( A, B ) RESULT( HET_N2O5 )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_N2O5
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: TMP1, TMP2

      ! Initialize
      HET_N2O5     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp

      ! Always apply adjustment
      KII_KI = .TRUE.

      ! Loop over aerosol types
      DO N = 1, NAERO

         ! Get GAMMA for N2O5 hydrolysis, which is
         ! a function of aerosol type, temp, and RH
         IF (N.eq.14) THEN
            IF (NATSURFACE) THEN
               XSTKCF = 4.0e-4_fp ! NAT
            ELSE
               XSTKCF = 0.02e+0_fp ! Ice
            ENDIF
         ELSEIF (N.eq.13) THEN
            ! Stratospheric aerosol
            XSTKCF = KHETI_SLA(1)
         ELSE
            ! In UCX, ABSHUMK will have been set by
            ! STT(I,J,L,IDTH2O)
            XSTKCF = N2O5( N, TEMPK, RELHUM )
         ENDIF
         ! Nitrate effect; reduce the gamma on nitrate by a
         ! factor of 10 (lzh, 10/25/2011)
         IF ( N == 8 ) THEN
            ! WARNING! It appears that these should be in units of
            !          mcl/cc. This is discerned from output in
            !          the old calcrate.F routine which gets the
            !          tracer concentrations from the state_chm object.
            !          When the values are dumped in calcrate.F, they were
            !          in mcl/cc. Still, comments later in calcrate.F indicate
            !          that Tracers should be in units of kg/box.
            ! -- As a fix, here, we simply impose the equivalent of a kg to
            !    mcl/cc conversion using the SO4 and NIT molecular weights.
            !    This should be investigated and the proper units applied.
            !    It will and does have a large impct on heterogenous
            !    N chemistry and on NOx in remote regions.
            ! -- In any case, the result from below yields the same
            !    ratio of TMP2/TMP1 as calcrate does with its
            !    current settings.
            !    MSL - Feb. 16, 2016
            TMP1 = (TRC_NIT*XNM_NIT)+(TRC_SO4*XNM_SO4)
            TMP2 = TRC_NIT*XNM_NIT
            IF ( TMP1 .GT. 0.0 ) THEN
               XSTKCF = XSTKCF * ( 1.0e+0_fp - 0.9e+0_fp &
                                   *TMP2/TMP1 )
            ENDIF
         ENDIF
         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF
         
         IF (KII_KI .and. N.gt.12) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
               ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
            ENDIF
         ENDIF
         
         ! Add to overall reaction rate
         HET_N2O5 = HET_N2O5 + ADJUSTEDRATE

      END DO

    END FUNCTIOn HETN2O5
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hetbrno3
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for BrNO3.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETBrNO3( A, B ) RESULT( HET_BrNO3 )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_BrNO3
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Initialize
      HET_BrNO3    = 0.0_fp
      ADJUSTEDRATE = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAERO

         ! Only apply adjustment if at high altitude
         KII_KI = STRATBOX
         
         ! Get the aerosol type
         ! If it's sulfate then use 0.8 for alpha, following
         !  JPL 2006 kinetics evaluation... holds for many
         !  temperatures and percent weights of sulfate.
         ! If not, then use the IUPAC recommendation of
         !  0.3, which is an input in globchem.dat
         ! (jpp, 5/4/10)
         IF ( N == 8 ) THEN
            ! sulfate aerosol
            XSTKCF = 0.8e+0_fp
         ELSE IF ( (N == 11) .OR. ( N == 12) ) THEN
            ! 2 modes of sea-salt
            XSTKCF = B
         ELSE IF ( N == 13 ) THEN
            ! SSA/STS
            XSTKCF = KHETI_SLA(6)
         ELSE IF ( N == 14 ) THEN 
            ! Ice/NAT PSC
            IF (NATSURFACE) THEN 
               XSTKCF = 0.001e+0_fp
            ELSE
               XSTKCF = 0.3e+0_fp
            ENDIF
         ELSE
            XSTKCF = 0e+0_fp
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF
         
         IF (KII_KI .and. N.gt.12) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
               ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
            ENDIF
         ENDIF
         
         ! Add to overall reaction rate
         HET_BrNO3 = HET_BrNO3 + ADJUSTEDRATE
      END DO
      IF (.not.PSCBOX) THEN
         HET_BrNO3 = HET_BrNO3 + cld_brno3_rc
      ENDIF

    END FUNCTIOn HETBrNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hethobr
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for HOBr.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOBr( A, B ) RESULT( HET_HOBr )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_HOBr
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Initialize
      HET_HOBr     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAERO

         ! jpp, 3/22/11: set the sticking coefficient to 
         !  ~0 for aerosol types we don't want reactions on
         !  for the HBr and HOBr surface reaction

         ! Only apply adjustment if at high altitude
         KII_KI = STRATBOX
         
         ! Select proper aerosol type
         IF ( (N == 8) .OR. (N == 11) .OR. (N == 12)) THEN
            ! sulfate, 2 modes of sea-salt
            XSTKCF = B
         ELSEIF ( N == 13 ) THEN
            XSTKCF = KHETI_SLA(11)
         ELSEIF ( N == 14 ) THEN
            XSTKCF = 0.1e+0_fp
         ELSE
            XSTKCF = 0e+0_fp
         ENDIF
         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF
         
         IF (KII_KI .and. N.gt.12) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
               ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
            ENDIF
         ENDIF
         
         ! Add to overall reaction rate
         HET_HOBr = HET_HOBr + ADJUSTEDRATE

      END DO

    END FUNCTIOn HETHOBr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hethobr_ice
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for HOBr (on ice).
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOBr_ice( A, B ) RESULT( HET_HObr_ice )

!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
! 
      REAL(fp)             :: HET_HObr_ice
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      HET_HOBr_ice = KI_HOBr

    END FUNCTIOn HETHOBr_ice
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hetbr_ice
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for HBr (on ice).
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHBr_ice( A, B ) RESULT( HET_HBr_ice )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE
!
      REAL(fp)             :: HET_HBr_ice
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      HET_HBr_ice = KI_HBr
      scf2(3)     = KI_HBr

    END FUNCTIOn HETHBr_ice
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: n2o5
!
! !DESCRIPTION: Internal function N2O5 computes the GAMMA sticking factor
!  for N2O5 hydrolysis. (mje, bmy, 8/7/03)
!\\
!\\
! !INTERFACE:
!
    FUNCTION N2O5( AEROTYPE, TEMP, RH ) RESULT( GAMMA )
!
! !INPUT PARAMETERS: 
!
      INTEGER,   INTENT(IN) :: AEROTYPE  ! Denoting aerosol type (cf FAST_JX)
      REAL(fp),  INTENT(IN) :: TEMP      ! Temperature [K]
      REAL(fp),  INTENT(IN) :: RH        ! Relative humidity [1]
!
! !RETURN VALUE:
!
      REAL(fp)              :: GAMMA

!
! !REMARKS:
!  Taken from the old SMVGEAR function calcrate.F.
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local variables
      REAL(fp) :: RH_P, FACT, TTEMP

      !=================================================================
      ! N2O5 begins here!
      !=================================================================

      ! Convert RH to % (max = 100%)
      RH_P  = MIN( RH * 100e+0_fp, 100e0_fp )

      ! Default value
      GAMMA = 0.01e+0_fp

      ! Special handling for various aerosols
      SELECT CASE ( AEROTYPE )

         !----------------
         ! Dust 
         !----------------
         CASE ( 1, 2, 3, 4, 5, 6, 7 )      
                                
            ! Based on unpublished Crowley work
            GAMMA = 0.01e+0_fp

         !----------------
         ! Sulfate
         !----------------
         CASE ( 8 )            
    
            !===========================================================
            ! RH dependence from Kane et al., Heterogenous uptake of 
            ! gaseous N2O5 by (NH4)2SO4, NH4HSO4 and H2SO4 aerosols
            ! J. Phys. Chem. A , 2001, 105, 6465-6470 
            !===========================================================
            ! No RH dependence above 50.0% (lzh, 10/26/2011)
            ! According to Bertram and Thornton, ACP, 9, 8351-8363, 2009
            RH_P  = MIN( RH_P, 50e+0_fp )

            GAMMA = 2.79e-4_fp + RH_P*(  1.30e-4_fp +    &
                              RH_P*( -3.43e-6_fp +       &
                              RH_P*(  7.52e-8_fp ) ) )

            !===========================================================
            ! Temperature dependence factor (Cox et al, Cambridge UK) 
            ! is of the form:
            !
            !          10^( LOG10( G294 ) - 0.04 * ( TTEMP - 294 ) )
            ! FACT = -------------------------------------------------
            !                     10^( LOG10( G294 ) )
            !
            ! Where G294 = 1e-2 and TTEMP is MAX( TEMP, 282 ).
            ! 
            ! For computational speed, replace LOG10( 1e-2 ) with -2
            ! and replace 10^( LOG10( G294 ) ) with G294 
            !===========================================================
            TTEMP = MAX( TEMP, 282e0_fp )
            FACT  = 10.e0_fp**( -2e+0_fp - 4e-2_fp       &
                  *( TTEMP - 294.e+0_fp ) ) / 1e-2_fp

            ! Apply temperature dependence
            GAMMA = GAMMA * FACT

         !----------------
         ! Black Carbon
         !----------------
         CASE ( 9 )  

             ! From IUPAC
             GAMMA = 0.005e+0_fp

         !----------------
         ! Organic Carbon
         !----------------           
         CASE ( 10 )          

            !===========================================================
            ! Based on Thornton, Braban and Abbatt, 2003
            ! N2O5 hydrolysis on sub-micron organic aerosol: the effect
            ! of relative humidity, particle phase and particle size
            !===========================================================
            IF ( RH_P >= 57e+0_fp ) THEN
               GAMMA = 0.03e+0_fp
            ELSE
               GAMMA = RH_P * 5.2e-4_fp
            ENDIF

         !----------------
         ! Sea salt
         ! accum & coarse
         !----------------
         CASE ( 11, 12 )        
          
            ! Based on IUPAC recomendation
            IF ( RH_P >= 62 ) THEN 
               GAMMA = 0.03e+0_fp
            ELSE
               GAMMA = 0.005e+0_fp
            ENDIF

         !----------------
         ! Strat. aerosols
         !----------------
         CASE ( 13, 14 )
       
            ! These are handled outside this routine - something
            ! is wrong if AEROTYPE=13 or 14 reaches this point
            WRITE (6,*) 'Stratospheric aerosols should not '
            WRITE (6,*) 'be passed to general N2O5 het. '
            WRITE (6,*) 'chem. subroutine'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            CALL GEOS_CHEM_STOP

         !----------------         
         ! Default
         !----------------
         CASE DEFAULT
            WRITE (6,*) 'Not a suitable aerosol surface '
            WRITE (6,*) 'for N2O5 hydrolysis'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            CALL GEOS_CHEM_STOP

      END SELECT   
         
    END FUNCTION N2O5
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ho2
!
! !DESCRIPTION: Function HO2 computes the GAMMA reaction probability
!  for HO2 loss in aerosols based on the recommendation of 
!  Thornton, Jaegle, and McNeill, "Assessing Known Pathways For HO2 Loss in 
!  Aqueous Atmospheric Aerosols: Regional and Global Impacts on Tropospheric 
!  Oxidants" J. Geophys. Res.,  doi:10.1029/2007JD009236, 2008  
!\\
!\\
! !INTERFACE:
!
      FUNCTION HO2( RADIUS,          TEMP,     DENAIR,       &
                    SQM,             HO2DENS,  AEROTYPE,     &
                    CONTINENTAL_PBL, Input_Opt            )  &
                    RESULT( GAMMA )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS: 
!
      ! Arguments
      REAL(fp),       INTENT(IN) :: RADIUS          ! Aerosol radius [cm]
      REAL(fp),       INTENT(IN) :: TEMP            ! Temperature [K]
      REAL(fp),       INTENT(IN) :: DENAIR          ! Air density [molec/cm3]
      REAL(fp),       INTENT(IN) :: HO2DENS         ! HO2 density [molec/cm3]
      REAL(fp),       INTENT(IN) :: SQM             ! Square root of MW [g/mol]
      INTEGER,        INTENT(IN) :: AEROTYPE        ! Aerosol type (cf FAST-JX)
      INTEGER,        INTENT(IN) :: CONTINENTAL_PBL ! Flag set to 1 if the box
                                                    !  box is located in the 
                                                    !  continenal boundary 
                                                    !  layer, otherwise 0.
                                                    !  Also check for ICE/SNOW
                                                    !  (to disable this at 
                                                    !  high latitudes).
      TYPE(OptInput), INTENT(IN) :: Input_Opt       ! Input Options object
!
! !RETURN VALUE:
!
      REAL(fp)                   :: GAMMA           ! Reaction probability
      
! !REMARKS:
!  Taken from the old SMVGEAR routine calcrate.F.
!  Gamma(HO2) is a function of aerosol type, radius, temperature.

!  eferences:
!  ---------------------------------------------------------------
!  (1) Jacob, D.J., Heterogeneous chemistry and tropospheric ozone,
!       Atmos. Environ., 34, 2131-2159, 2000. [full text (pdf)]
!  (2) J. Mao, Fan, S., Jacob, D. J., and Travis, K. R.: Radical
!       loss in the atmosphere from Cu-Fe redox coupling in aerosols,
!       Atmos. Chem. Phys., 13, 509-519, doi:10.5194/acp-13-509-2013,
!       2013.
!
! !REVISION HISTORY:
!  17 May 2013 - M. Payer    - Add improved HO2 uptake (J. Mao)
!  22 May 2013 - M. Payer    - Added option to read GAMMA_HO2 from
!                              input.geos. Recommended value is 0.2
!                              based on Jacob et al (2000) and Mao
!                              et al. (2013).
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  29 Mar 2016 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)             :: ALPHA
      REAL(fp)             :: delG, Keq, w, H_eff
      REAL(fp)             :: A1, B1, k1, k2, A, B, C
      REAL(fp)             :: kaq, kmt, o2_ss, fluxrxn, DFKG
      REAL(fp)             :: TEST
!
! !DEFINED PARAMETERS:
!
      !%%% NOTE: WE SHOULD EVENTUALLY USE THE VALUES FROM physconsts.F %%%

      ! Avogadro's number
      REAL(fp),  PARAMETER :: Na = 6.022e+23_fp

      ! Ideal gas constant [atm cm3/mol/K], Raq
      REAL(fp),  PARAMETER :: Raq=82.e+0_fp

      !=================================================================
      ! HO2 begins here!
      !=================================================================

      ! Default value
      GAMMA = 0.0e+0_fp

      ! Error check
      IF (RADIUS.le.1e-30_fp) THEN
         RETURN
      ENDIF

      ! Special handling for various aerosols
      SELECT CASE ( AEROTYPE )

         !----------------
         ! Dust 
         !----------------
         CASE ( 1, 2, 3, 4, 5, 6, 7 )      
                                
            ! Assume default gamma=0.1 on dust aerosols
            ! This is tentative as no lab measurements presently exist
            ! for gamma(HO2) on dust aerosols. We assume the rate to
            ! be fast on dust aerosols as transition metal ion induced
            ! chemistry is likely to occur in a thin aqueous surface layer.
            GAMMA = 0.1e+0_fp

         !----------------
         ! For Sulfate(8), Black Carbon (9), Organic Carbon (10),
         ! Sea-salt accum & coarse (11,12) calculate the 
         ! reaction probability due to self reaction 
         ! by using the algebraic expression in Thornton et al.  (2008)
         ! (equation 7) which is a function of temperature, aerosol radius,
         ! air density and HO2 concentration. 
         !
         ! Transition metal ions (such as copper and iron) in sea-salt and 
         ! carbonaceous aerosols are complexed to ligands and/or exist at 
         ! a concentration too low to catalyze HO2 loss efficiently, so we 
         ! apply the HO2 self reaction expression directly for these aerosols.
         ! 
         ! In the case of sulfate aerosol, the aerosols likely
         ! contain copper in the continental boundary layer and
         ! HO2 uptake proceeds rapidly. To account for the metal catalyzed
         ! uptake, we assume gamma(HO2)=0.07 (in the mid-range of the recommended
         ! 0.04-0.1 by Thornton et al, based on observed copper concentrations
         ! in the US boundary layer). Outside the continental boundary layer, we
         ! use the HO2-only algebraic expression.
         !
         ! SDE 04/18/13: Added stratospheric sulfur aerosols
         !
         !----------------
         CASE ( 8, 9, 10, 11, 12, 13 )  

            ! Mean molecular speed [cm/s]
            w = 14550.5e+0_fp * sqrt(TEMP/(SQM*SQM))

            ! DFKG = Gas phase diffusion coeff [cm2/s]
            DFKG  = 9.45E+17_fp/DENAIR * SQRT(TEMP) *      &
                    SQRT(3.472E-2_fp + 1.E+0_fp/(SQM*SQM))

            !calculate T-dependent solubility and aq. reaction rate constants
            ! hydronium ion concentration
            ! A1 = 1.+(Keq/hplus) 
            ! with Keq = 2.1d-5 [M], Equilibrium constant for 
            ! HO2aq = H+ + O2- (Jacob, 2000)
            !      hplus=10.e+0_fp^(-pH), with pH = 5
            ! B1 = Req * TEMP
            ! with Req = 1.987d-3 [kcal/K/mol], Ideal gas constant
            ! Note that we assume a constant pH of 5.
            A1 = 1.+ (2.1e-5_fp / (10.e+0_fp**(-5) ) )
            B1 = 1.987e-3_fp * TEMP

            ! Free energy change for solvation of HO2 (Hanson 1992, Golden 1991)
            ! in [kcal/mol]:
            ! delG = -4.9-(TEMP-298e+0_fp)*delS
            ! with delS=-0.023  [kcal/mol/K],  Entropy change for solvation of HO2
            delG  = -4.9e+0_fp - (TEMP-298.e+0_fp) * (-0.023)
            H_eff = exp( -delG / B1 ) * A1

            ! Estimated temp dependent value for HO2 + O2- (k1) and 
            ! HO2+HO2 (see Jacob 1989)
            k1  =   1.58e+10_fp * exp( -3. / B1 )
            k2  =   2.4e+9_fp   * exp( -4.7 / B1 )
            kaq = ( k1 * (A1 - 1.e+0_fp) + k2) / (A1**2)

            ! Calculate the mass transfer rate constant and s.s. conc. of 
            ! total HO2 in the aqueous phase:
            ! kmt = (RADIUS/DFKG + 4e+0_fp/w/alpha)^(-1)
            ! with alpha = mass accomodation coefficient, assumed 
            ! to be 1 (Thornton et al.)
            kmt = 1.e+0_fp/( RADIUS/DFKG + 4e+0_fp/w/1. )

            !use quadratic formula to obtain [O2-] in particle of radius RADIUS
            A = -2e+0_fp * kaq
            B = -3e+0_fp * kmt / RADIUS / (H_eff * 0.082 * TEMP)
            C =  3e+0_fp * kmt * HO2DENS * 1000e+0_fp / RADIUS / Na

            ! Error check that B^2-(4e+0_fp*A*C) is not negative
            TEST= B**2-(4e+0_fp*A*C)
            IF ( TEST < 0e+0_fp ) THEN
                GAMMA = 0e+0_fp
            ELSE
                ! Calculate the concentration of O2- in the aerosol
                o2_ss= ( -B  -sqrt(B**2-(4e+0_fp*A*C)) )/(2e+0_fp*A)

                ! Calculate the reactive flux
                fluxrxn = kmt*HO2DENS - o2_ss*Na*kmt/H_eff/Raq/TEMP

                IF ( fluxrxn <= 0e0_fp ) THEN
                   GAMMA = 0e+0_fp
                ELSE
                   ! Gamma for HO2 at TEMP, ho2, and RADIUS given
                   GAMMA = 1./( ( ( HO2DENS/fluxrxn ) -              &
                                  ( RADIUS/DFKG ) ) * w / 4.e+0_fp )
                ENDIF
            ENDIF
            ! For sulfate aerosols, check whether we are in
            ! the continental boundary layer, in which case
            ! copper catalyzed HO2 uptake likely dominates and
            ! speeds up the reaction: we assume gamma=0.07,
            ! which is in the middle of the 0.04-0.1 range recommended
            ! by Thornton et al. (2008)
            !
            IF ( AEROTYPE == 8 .and. CONTINENTAL_PBL == 1) THEN
                GAMMA = 0.07
            ENDIF 

         !----------------
         ! NAT/ice (SDE 04/18/13)
         !----------------
         CASE ( 14 )
       
            GAMMA = 0.e+0_fp

         !----------------
         ! Default
         !----------------
         CASE DEFAULT
            WRITE (6,*) 'Not a suitable aerosol surface '
            WRITE (6,*) 'for HO2 uptake'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            CALL GEOS_CHEM_STOP

      END SELECT
     
      ! If negative value is calculated, set it to zero
      IF ( GAMMA  <= 0e+0_fp ) GAMMA = 0e+0_fp

      ! This is for the improved HO2 uptake (J. Mao)
      GAMMA = Input_Opt%GAMMA_HO2

    END FUNCTION HO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cld1k_brno3
!
! !DESCRIPTION: Function CLD1K_BrNO3 calculates the rate constant for
!  heterogeneous cycling of BrNO3 off of cloud particles.
!\\
!\\
! !INTERFACE:
!
    FUNCTION CLD1K_BrNO3( I,      J,        L,   DENAIR,         &
                          QL,     State_Met )    RESULT( cld1k )

!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN) :: I         ! Longitude index
      INTEGER,        INTENT(IN) :: J         ! Latitude  index
      INTEGER,        INTENT(IN) :: L         ! Altitude  index
      REAL(fp),       INTENT(IN) :: DENAIR    ! Density of air [#/cm3]
      REAL(fp),       INTENT(IN) :: QL        ! Cloud water mixing ratio [kg/kg]
      TYPE(MetState), INTENT(IN) :: State_Met ! Meteorology State object
!
! !RETURN VALUE:
!
      REAL(fp)              :: cld1k          ! Rate constant for 
                                              ! heterogeneous cycling
                                              ! of BrNO3 off of cloud 
!                                             ! particles
!
! !REMARKS:
!  The rate constant for heterogeneous cycling of BrNO3 off of cloud particles
!  is calculated assuming:
!                                                                             .
!    1. A sticking coefficient of 0.3 [Yang et al. 2005]
!    2. uniform cloud droplet size for 2 types of clouds
!       - continental warm clouds: r =  6d-4 [cm]
!       - marine warm clouds:      r = 10d-4 [cm]
!       * no distributions are assumed
!
!  ** Calculation of a 1st order rate constent barrowed from the
!     subroutine arsl1k.f. Below are comments from that code:
!                                                                             .
!       The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
!       is computed as:
!                                                                             .
!         ARSL1K [1/s] = area / [ radius/dfkg + 4./(stkcf * nu) ]        
!                                                                             .
!       where nu   = Mean molecular speed [cm/s] = sqrt(8R*TK/pi/M) for Maxwell
!             DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
!
! !REVISION HISTORY:
!  27 Feb 2011 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  06 Nov 2014 - R. Yantosca - Now use State_Met%CLDF(I,J,L)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Cloud droplet radius in continental warm clouds [cm]
      REAL(fp), PARAMETER :: XCLDR_CONT =  6.e-4_fp

      ! Cloud droplet radius in marine warm clouds [cm]
      REAL(fp), PARAMETER :: XCLDR_MARI = 10.e-4_fp

      !%%% NOTE: WE SHOULD EVENTUALLY USE THE VALUES FROM physconsts.F %%%

      REAL(fp), PARAMETER :: R = 8.314472                  ! [J/mol/K]
      REAL(fp), PARAMETER :: mw_brno3 = 0.142              ! [kg/mol]
      rEAL(fp), PARAMETER :: pi = 3.14159265358979323846e+0_fp ! [unitless]
      REAL(fp), PARAMETER :: alpha = 0.3                   ! sticking coefficient
      REAL(fp), PARAMETER :: dens_h2o = 0.001e+0_fp            ! [kg/cm3]
!
! !LOCAL VARIABLES:
!
      REAL(fp)            :: nu         ! Mean molecular speed
      REAL(fp)            :: RADIUS     ! Radius of cloud droplet      [cm]
      REAL(fp)            :: SQM        ! Square root of molec. weight [g/mol]
      REAL(fp)            :: STK        ! Square root of temperature   [K]
      REAL(fp)            :: AREA       ! Surface area                 [m2]
      REAL(fp)            :: DFKG       ! Gas diffusion coefficient    [cm2/s]
      REAL(fp)            :: Vc         ! Volume of the cloud          [cm3]
      REAL(fp)            :: XAIRM3     ! Volume of air                [m3]
      LOGICAL             :: yn_continue, IS_LAND, IS_ICE
   
      ! Pointers
      REAL(fp), POINTER   :: AD(:,:,:)
      REAL(fp), POINTER   :: AIRVOL(:,:,:)
      REAL(fp), POINTER   :: CLDF(:,:,:)
      REAL(fp), POINTER   :: FRLAND(:,:)
      REAL(fp), POINTER   :: FROCEAN(:,:)
      REAL(fp), POINTER   :: T(:,:,:)

      !=================================================================
      ! CLD1K_BrNO3 begins here!
      !=================================================================

      ! Initialize pointers
      AD      => State_Met%AD
      AIRVOL  => State_Met%AIRVOL
      CLDF    => State_Met%CLDF
      FRLAND  => State_Met%FRLAND
      FROCEAN => State_Met%FROCEAN
      T       => State_Met%T

      ! -- IS THIS LAND? -- (Adapted from DAO_MOD function)
#if   defined( GCAP )

      !--------------------------
      ! GCAP
      !--------------------------

      ! It's a land box if 50% or more of the box is covered by 
      ! land and less than 50% of the box is covered by ice
      IS_LAND = ( State_Met%LWI_GISS(I,J) >= 0.5e+0_fp .and. &
                  State_Met%SNICE(I,J)    <  0.5e+0_fp )

#else

      !--------------------------
      ! GEOS-4 / GEOS-5 / MERRA
      !--------------------------

      ! LWI=1 and ALBEDO less than 69.5% is a LAND box 
      IS_LAND = ( NINT( State_Met%LWI(I,J) ) == 1       .and. &
                     State_Met%ALBD(I,J)  <  0.695e+0_fp )

#endif
      ! Done with 'Is this land' ---------------------
      ! -- IS THIS ICE? -- (Adapted from DAO_MOD function)
#if   defined( GCAP )

      !--------------------------
      ! GCAP
      !--------------------------

      ! It's an ice box if 50% or more of the box is covered by ice
      IS_ICE = ( State_Met%SNICE(I,J) >= 0.5e+0_fp )

#else

      !--------------------------
      ! GEOS-4 / GEOS-5 / MERRA
      !--------------------------

      ! LWI=2 or ALBEDO > 69.5% is ice
      IS_ICE = ( NINT( State_Met%LWI(I,J) ) == 2       .or. &
                    State_Met%ALBD(I,J)  >= 0.695e+0_fp )

#endif
      ! Done with 'Is this ice' ---------------------

      ! ----------------------------------------------
      ! 1.
      !   calculate the mean molecular speed of the
      !   molecules given the temperature.
      ! ----------------------------------------------
      nu   = sqrt( 8.e+0_fp * R * T(I,J,L) / (mw_brno3 * pi) )

      ! ----------------------------------------------
      ! Test conditions to see if we want to continue
      ! or set the cloud rate equal to zero.
      ! ----------------------------------------------

      ! continental or marine clouds only...
#if defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
      IF ( (FRLAND (I,J) > 0) .or. (FROCEAN(I,J) > 0) ) THEN
#else
      ! Above line is to skip over land ice (Greenland and Antartica). This
      ! should do the same (and also work for GEOS-5, but leave above for now).
      IF ( IS_LAND .and. .not. IS_ICE  ) THEN
#endif
         ! do we have clouds? and do we have warm temperatures?
         IF ( ( CLDF(I,J,L) > 0    )   .and.           &
              ( T(I,J,L)    > 258.0) ) THEN
            yn_continue = .TRUE.
         ELSE
            yn_continue = .FALSE.
         ENDIF
      ELSE
         yn_continue = .FALSE.
      ENDIF

      ! test
      IF ( .not. yn_continue ) THEN
         ! nothing to calculate...
         cld1k = 0.e+0_fp
         RETURN
      ENDIF


      ! ----------------------------------------------
      ! 2.
      !   calculate the surface area of cloud droplets
      !   in the given grid box, assuming 1 of 2
      !   conditions:
      !     a. marine warm cloud
      !       or
      !     b. continental warm cloud
      !
      !
      !   * Calculation for area is derived follows,
      !     assuming that RADIUS is constant:
      !
      !                         4/3 (pi) (RADIUS)**3
      !  1) FC = Vc / Vb = N  -------------------------
      !                                  Vb
      !
      !
      !       where N      = number of cloud droplets
      !             RADIUS = radius of cloud droplet
      !             Vc     = volumn of the cloud
      !             Vb     = volumn of the box = AIRVOL (in GEOS-Chem)
      !
      !
      !                     Vb
      !  2) N = FC --------------------
      !            4/3 (pi) (RADIUS)**3
      !
      !
      !  So the surface area [m2] is calculated as
      !
      !  3) total surface A = N * 4 * (pi) * (RADIUS)**2
      !
      !                  3*Vb
      !          = FC ----------
      !                 RADIUS
      !
      !  4) for this routine though we want
      !     AREA in [cm2/cm3], surface area to volume air:
      !
      !                   3
      !     AREA = FC ---------
      !                RADIUS (in cm)
      !
      !
      !    or    
      !                   3 x Vc
      !     AREA =  -----------------
      !              AIRVOL x RADIUS      (in cm)
      ! ----------------------------------------------
#if defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
      IF ( FRLAND(I,J) > FROCEAN(I,J) ) THEN
#else
      IF ( IS_LAND ) THEN
#endif
         ! Continental cloud droplet radius [cm]
         RADIUS = XCLDR_CONT
      ELSE
         ! Marine cloud droplet radius [cm]
         RADIUS = XCLDR_MARI
      ENDIF

      ! store the volume of air [m3]
      XAIRM3 = AIRVOL(I,J,L)
      ! convert to [cm3]
      XAIRM3 = XAIRM3 * (100.e+0_fp)**3

      ! get the volume of cloud [cm3]
#if defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
      ! QL is [g/g]
      Vc = CLDF(I,J,L) * QL * AD(I,J,L) / dens_h2o
#else
      ! QL is [cm3/cm3]
      Vc = CLDF(I,J,L) * QL * XAIRM3
#endif

      ! now calculate the cloud droplet surface area
      AREA    = 3.e+0_fp * (Vc/XAIRM3) / (RADIUS) ! keep Radius in [cm]

      ! ----------------------------------------------------
      ! 3.
      !   Now finish calculating the 1st order rate
      !   constant for BrNO3 hydrolysis.
      !
      !   (a) calculate the gas phase diffusion coefficient;
      !
      !   (b) calculate the hydrolysis rxn rate.
      ! ----------------------------------------------------
      SQM = sqrt(mw_brno3 * 1.e+3_fp)    ! square root of molar mass [g/mole]
      STK = sqrt(T(I,J,L)) ! square root of temperature [K]

      ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
      DFKG  = 9.45E+17_fp/DENAIR * STK * SQRT(3.472E-2_fp     &
           + 1.E+0_fp/(SQM*SQM))

      ! Compute ARSL1K according to the formula listed above
      cld1k = AREA / ( RADIUS/DFKG + 2.749064E-4              &
           * SQM/(alpha*STK) )

      ! Free Pointers
      NULLIFY( AD      )
      NULLIFY( AIRVOL  )
      NULLIFY( CLDF    )
      NULLIFY( FRLAND  )
      NULLIFY( FROCEAN )
      NULLIFY( T       )

    END FUNCTION CLD1K_BrNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fcro2ho2
!
! !DESCRIPTION: !fgap, based on saunder 2003 k14
!\\
!\\
! !INTERFACE:
!
    FUNCTION FCRO2HO2( XCARBN ) RESULT( FC_RO2HO2 )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: XCARBN
!
! !RETURN VALUE:
!
      REAL(fp)             :: FC_RO2HO2
!
! !REVISION HISTORY:
!  24 Jul 2014 - R. Yantosca - Now inlined to calcrate.F
!EOP
!------------------------------------------------------------------------------
!BOC

      FC_RO2HO2 = 1E+0_fp - EXP( -0.245E+0_fp * XCARBN )
     
    END FUNCTION FCRO2HO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: FYHORO
!
! !DESCRIPTION: \subsection*{Overview}
!  Function FYHORO returns returns the branching ratio between 
!  HOC2H4O oxidation and dissociation:
!  (1) HOC2H4 --O2--> HO2 + GLYC
!  (2) HOC2H4 ------> HO2 + 2CH2O
!
!\subsection*{References}
!  \begin{enumerate}
!  \item Orlando et al., 1998: \emph{Laboratory and theoretical study of the 
!         oxyradicals in the OH- and Cl-initiated oxidation of ethene}, 
!        \underline{J. Phys. Chem. A}, \textbf{102}, 8116-8123.
!  \item Orlando et al., 2003: \emph{The atmospheric chemistry of alkoxy 
!         radicals}, \underline{Chem. Rev.}, \textbf{103}, 4657-4689.
!  \end{enumerate}
!
!\\
!\\
! !INTERFACE: 
!
    FUNCTION FYHORO( ZDNUM, TT ) RESULT( FY_HORO )
!
! !INPUT PARAMETERS: 
!
      REAL(fp), INTENT(IN) :: ZDNUM   ! Air density   [molec/cm3 ]
      REAL(fp), INTENT(IN) :: TT      ! Temperature   [K         ]
!
! !RETURN VALUE:
!
      REAL(fp)             :: FY_HORO
!
! !REVISION HISTORY:
!  (1 ) Branching ratio calculation (tmf, 2/6/05).
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  25 Jul 2014 - R. Yantosca - Now inlined into calcrate.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: K1, K2, O2DNUM

      !=================================================================
      ! FYHORO begins here!
      !=================================================================
      O2DNUM  = ZDNUM * 0.21E+0_fp
      K1      = 6.0E-14_fp * EXP(-550.E+0_fp/TT) * O2DNUM
      K2      = 9.5E+13_fp * EXP(-5988.E+0_fp/TT) 

      FY_HORO = K1 / (K1 + K2)
        
    END FUNCTION FYHORO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fyrno3
!
! !DESCRIPTION: Function FYRNO3 returns organic nitrate yields 
!  YN = RKA/(RKA+RKB) from RO2+NO reactions as a function of the number 
!  N of carbon atoms.
!\\
!\\
! !INTERFACE:
!
    FUNCTION FYRNO3( XCARBN, ZDNUM, TT ) RESULT( FYR_NO3 )
!
! !INPUT PARAMETERS: 
!
      REAL(fp), INTENT(IN) :: XCARBN   ! Number of C atoms in RO2
      REAL(fp), INTENT(IN) :: ZDNUM    ! Air density   [molec/cm3 ]
      REAL(fp), INTENT(IN) :: TT       ! Temperature   [K         ]
!
! !RETURN VALUE:
!
      REAL(fp)             :: FYR_NO3 
! 
! !REVISION HISTORY: 
!  (1 ) Original code from Larry Horowitz, Jinyou Liang, Gerry Gardner,
!        and Daniel Jacob circa 1989/1990.
!  (2 ) Updated following Atkinson 1990.
!  (3 ) Change yield from Isoprene Nitrate (ISN2) from 0.44% to 12%,
!        according to Sprengnether et al., 2002. (amf, bmy, 1/7/02)
!  (4 ) Eliminate obsolete code from 1/02 (bmy, 2/27/02)
!  (5 ) Updated comment description of XCARBN (bmy, 6/26/03)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  25 Jul 2014 - R. Yantosca - Now inlined into calcrate.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: YYYN, XXYN,  AAA,  RARB, ZZYN
      REAL(fp) :: XF,   ALPHA, Y300, BETA, XMINF, XM0

      ! Initialize static variables
      DATA   Y300,ALPHA,BETA,XM0,XMINF,XF/.826,1.94E-22,.97,0.,8.1,.411/

      !=================================================================
      ! FYRNO3 begins here!
      !=================================================================
      XXYN    = ALPHA*EXP(BETA*XCARBN)*ZDNUM*((300./TT)**XM0)
      YYYN    = Y300*((300./TT)**XMINF)
      AAA     = LOG10(XXYN/YYYN)
      ZZYN    = 1./(1.+ AAA*AAA )
      RARB    = (XXYN/(1.+ (XXYN/YYYN)))*(XF**ZZYN)
      FYR_NO3 = RARB/(1. + RARB)
     
    END FUNCTION FYRNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: arsl1k
!
! !DESCRIPTION: Function ARSL1K calculates the 1st-order loss rate of species 
!  on wet aerosol surface.
!\\
!\\
! !INTERFACE:
!
    FUNCTION ARSL1K( AREA, RADIUS, DENAIR, STKCF, STK, SQM ) &
         RESULT( ARS_L1K )
!
! !INPUT PARAMETERS: 
!
      ! Surface  area of wet aerosols/volume of air [cm2/cm3]
      REAL(fp), INTENT(IN) :: AREA     

      ! Radius of wet aerosol [cm], order of 0.01-10 um;
      ! Note that radius here is Rd, not Ro
      REAL(fp), INTENT(IN) :: RADIUS 
  
      ! Density of air [#/cm3]
      REAL(fp), INTENT(IN) :: DENAIR  
 
      ! Sticking coefficient [unitless], order of 0.1
      REAL(fp), INTENT(IN) :: STKCF  
  
      ! Square root of temperature [K]
      REAL(fp), INTENT(IN) :: STK  
    
      ! Square root of molecular weight [g/mole]
      REAL(fp), INTENT(IN) :: SQM      
!
! !RETURN VALUE:
!
      REAL(fp)             :: ARS_L1K
!
! !REMARKS:
!  The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
!  is computed as:
!                                                                             .
!      ARSL1K [1/s] = area / [ radius/dfkg + 4./(stkcf * xmms) ]        
!                                                                             .
!  where XMMS = Mean molecular speed [cm/s] = sqrt(8R*TK/pi/M) for Maxwell 
!        DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)

! !REVISION HISTORY:
!  01 Jul 1994 - lwh, jyl, gmg, djj - Initial version 
!  04 Apr 2003 - R. Yantosca - Updated comments, cosmetic changes
!  07 Apr 2004 - R. Yantosca - Now return w/ default value if RADIUS is zero 
!                              (i.e. is smaller than a very small number)
!  03 Dec 2009 - R. Yantosca - Prevent div-by-zero errors by returning the
!                              default value if any of the args are zero 
!  03 Dec 2009 - R. Yantosca - Added ProTeX Header
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: DFKG

      !=================================================================
      ! ARSL1K begins here!
      !=================================================================
      IF ( AREA < 0e0_fp      .or. DENAIR < 1e-30_fp    .or.  & 
           RADIUS < 1e-30_fp  .or. SQM  < 1e-30_fp      .or.  &
           STK    < 1e-30_fp  .or.  STKCF  < 1e-30_fp ) THEN

         ! Use default value if any of the above values are zero
         ! This will prevent div-by-zero errors in the eqns below
         ! Value changed from 1d-3 to 1d-30 (bhh, jmao, eam, 7/18/2011)
         ARS_L1K = 1.E-30_fp

      ELSE

         ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
         DFKG  = 9.45E+17_fp/DENAIR * STK * SQRT(3.472E-2_fp + 1.E0_fp/ &
          (SQM*SQM))

         ! Compute ARSL1K according to the formula listed above
         ARS_L1K = AREA / ( RADIUS/DFKG + 2.749064E-4*SQM/(STKCF*STK) )

      ENDIF

    END FUNCTION ARSL1K
!EOC
!##############################################################################
!###                                                                        ###
!###   THE FOLLOWING FUNCTIONS ARE ONLY DEFINED FOR UCX-BASED MECHANISMS    ### 
!###                                                                        ###
!##############################################################################
#if defined( UCX )
!BOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hetn2o5_psc
!
! !DESCRIPTION: Set heterogenous chemistry rate for N2O5(g) + HCl(l,s)
!  in polar stratospheric clouds.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETN2O5_PSC( A, B ) RESULT( HET_N2O5_PSC )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_N2O5_PSC
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms.
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Initialize
      HET_N2O5_PSC = 0.0_fp
      ADJUSTEDRATE = 0.0_fp

      ! Always apply adjustment
      KII_KI = .TRUE.

      ! Only consider PSC reactions in strat
      IF (STRATBOX) THEN

         ! Loop over aerosol types
         DO N = 1, NAERO

            IF (N.eq.8) THEN
               XSTKCF = 0.1e-4_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(2)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.003e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.03e+0_fp ! Ice
               ENDIF
            ELSE
               XSTKCF = 0e+0_fp
            ENDIF

            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            IF (KII_KI .and. N.gt.12) THEN
               ! PSC reaction - prevent excessive reaction rate
               IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
                  ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
               ENDIF
            ENDIF

            ! Add to overall reaction rate
            HET_N2O5_PSC = HET_N2O5_PSC + ADJUSTEDRATE

         END DO

      ENDIF

    END FUNCTION HETN2O5_PSC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hetclno3_psc1
!
! !DESCRIPTION: Set heterogenous chemistry rate for ClNO3(g) + H2O(l,s)
!  in polar stratopsheric clouds.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETClNO3_PSC1( A, B ) RESULT( HET_ClNO3_PSC1 )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)            :: HET_ClNO3_PSC1
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms.
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Initialize
      HET_ClNO3_PSC1 = 0.0_fp
      ADJUSTEDRATE   = 0.0_fp

      ! Always apply adjustment
      KII_KI = .TRUE.

      ! Only consider PSC reactions in strat
      IF (STRATBOX) THEN

         ! Loop over aerosol types
         DO N = 1, NAERO

            IF (N.eq.8) THEN
               XSTKCF = 0.1e-3_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(3)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.004e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
            ELSE
               XSTKCF = 0e+0_fp
            ENDIF

            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            IF (KII_KI .and. N.gt.12) THEN
               ! PSC reaction - prevent excessive reaction rate
               IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
                  ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
               ENDIF
            ENDIF

            ! Add to overall reaction rate
            HET_ClNO3_PSC1 = HET_ClNO3_PSC1 + ADJUSTEDRATE

         END DO

      ENDIF

    END FUNCTION HETClNO3_PSC1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hetclno3_psc2
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for ClNO3(g) + HCl(l,s)
! in polar stratospheric clouds.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETClNO3_PSC2( A, B ) RESULT( HET_ClNO3_PSC2 )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_ClNO3_PSC2
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms.
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Initialize
      HET_ClNO3_PSC2 = 0.0_fp
      ADJUSTEDRATE   = 0.0_fp

      ! Always apply adjustment
      KII_KI = .TRUE.

      ! Only consider PSC reactions in strat
      IF (STRATBOX) THEN

         ! Loop over aerosol types
         DO N = 1, NAERO

            IF (N.eq.8) THEN
               XSTKCF = 0.1e-4_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(4)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.2e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
            ELSE
               XSTKCF = 0e+0_fp
            ENDIF

            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            IF (KII_KI .and. N.gt.12) THEN
               ! PSC reaction - prevent excessive reaction rate
               IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
                  ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
               ENDIF
            ENDIF

            ! Add to overall reaction rate
            HET_ClNO3_PSC2 = HET_ClNO3_PSC2 + ADJUSTEDRATE

         END DO

      ENDIF

    END FUNCTION HETClNO3_PSC2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hetclno3_psc3
!
! !DESCRIPTION: Set heterogenous chemistry rate for ClNO3(g) + HBr(l,s)
!  in polar stratospheric clouds.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETClNO3_PSC3( A, B ) RESULT( HET_ClNO3_PSC3 )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_ClNO3_PSC3
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms.
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Initialize
      HET_ClNO3_PSC3 = 0.0_fp
      ADJUSTEDRATE   = 0.0_fp

      ! Always apply adjustment
      KII_KI = .TRUE.

      ! Only consider PSC reactions in strat
      IF (STRATBOX) THEN

         ! Loop over aerosol types
         DO N = 1, NAERO

            IF (N.eq.8) THEN
               XSTKCF = 0.e+0_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(5)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.3e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
            ELSE
               XSTKCF = 0e+0_fp
            ENDIF

            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            IF (KII_KI .and. N.gt.12) THEN
               ! PSC reaction - prevent excessive reaction rate
               IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
                  ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
               ENDIF
            ENDIF

            ! Add to overall reaction rate
            HET_ClNO3_PSC3 = HET_ClNO3_PSC3 + ADJUSTEDRATE

         END DO

      ENDIF

    END FUNCTION HETClNO3_PSC3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hetbrno3_psc
!
! !DESCRIPTION: Set heterogenous chemistry rate for BrNO3(g) + HCl(l,s)
!  in polar stratospheric clouds.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETBrNO3_PSC( A, B ) RESULT( HET_BrNO3_PSC )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE: 
!
      REAL(fp)             :: HET_BrNO3_PSC
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms.
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Initialize
      HET_BrNO3_PSC = 0.0_fp
      ADJUSTEDRATE  = 0.0_fp

      ! Always apply adjustment
      KII_KI = .TRUE.

      ! Only consider PSC reactions in strat
      IF (STRATBOX) THEN

         ! Loop over aerosol types
         DO N = 1, NAERO

            IF (N.eq.8) THEN
               XSTKCF = 0.9e+0_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(7)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.3e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
            ELSE
               XSTKCF = 0e+0_fp
            ENDIF

            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            IF (KII_KI .and. N.gt.12) THEN
               ! PSC reaction - prevent excessive reaction rate
               IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
                  ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
               ENDIF
            ENDIF

            ! Add to overall reaction rate
            HET_BrNO3_PSC = HET_BrNO3_PSC + ADJUSTEDRATE

         END DO

      ENDIF

    END FUNCTION HETBrNO3_PSC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hethocl_psc1
!
! !DESCRIPTION: Set heterogenous chemistry rate for HOCl(g) + HCl(l,s)
!  in polar stratopsheric clouds.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOCl_PSC1( A, B ) RESULT( HET_HOCl_PSC1 )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)            :: HET_HOCl_PSC1
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms.
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Initialize
      HET_HOCl_PSC1 = 0.0_fp
      ADJUSTEDRATE  = 0.0_fp

      ! Always apply adjustment
      KII_KI = .TRUE.

      ! Only consider PSC reactions in strat
      IF (STRATBOX) THEN

         ! Loop over aerosol types
         DO N = 1, NAERO

            IF (N.eq.8) THEN
               XSTKCF = 0.8e+0_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(8)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.1e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.2e+0_fp ! Ice
               ENDIF
            ELSE
               XSTKCF = 0e+0_fp
            ENDIF

            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            IF (KII_KI .and. N.gt.12) THEN
               ! PSC reaction - prevent excessive reaction rate
               IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
                  ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
               ENDIF
            ENDIF

            ! Add to overall reaction rate
            HET_HOCl_PSC1 = HET_HOCl_PSC1 + ADJUSTEDRATE

         END DO

      ENDIF

    END FUNCTION HETHOCl_PSC1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hethocl_psc2
!
! !DESCRIPTION: Set heterogenous chemistry rate for HOCl(g) + HBr(l,s)
!  in polar stratospheric clouds.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOCl_PSC2( A, B ) RESULT( HET_HOCl_PSC2 )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_HOCl_PSC2
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms.
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC


      ! Initialize
      HET_HOCl_PSC2 = 0.0_fp
      ADJUSTEDRATE  = 0.0_fp

      ! Always apply adjustment
      KII_KI = .TRUE.

      ! Only consider PSC reactions in strat
      IF (STRATBOX) THEN

         ! Loop over aerosol types
         DO N = 1, NAERO

            IF (N.eq.8) THEN
               XSTKCF = 0.8e+0_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(9)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.3e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
            ELSE
               XSTKCF = 0e+0_fp
            ENDIF

            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            IF (KII_KI .and. N.gt.12) THEN
               ! PSC reaction - prevent excessive reaction rate
               IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
                  ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
               ENDIF
            ENDIF

            ! Add to overall reaction rate
            HET_HOCl_PSC2 = HET_HOCl_PSC2 + ADJUSTEDRATE

         END DO

      ENDIF

    END FUNCTION HETHOCl_PSC2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hethobr_psc
!
! !DESCRIPTION: Set heterogenous chemistry rate for HOBr(g) + HCl(l,s)
!  in polar stratospheric clouds.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOBr_PSC( A, B ) RESULT( HET_HOBr_PSC )
!
! !INPUT PARAMETERS: 
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_HOBr_PSC
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms.
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Initialize
      HET_HOBr_PSC = 0.0_fp
      ADJUSTEDRATE = 0.0_fp

      ! Always apply adjustment
      KII_KI = .TRUE.

      ! Only consider PSC reactions in strat
      IF (STRATBOX) THEN

         ! Loop over aerosol types
         DO N = 1, NAERO

            IF (N.eq.8) THEN
               XSTKCF = 0.8e+0_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(10)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.1e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
            ELSE
               XSTKCF = 0e+0_fp
            ENDIF

            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            IF (KII_KI .and. N.gt.12) THEN
               ! PSC reaction - prevent excessive reaction rate
               IF (ADJUSTEDRATE.gt.(1.e+0_fp/PSCMINLIFE)) THEN
                  ADJUSTEDRATE = 1.e+0_fp/PSCMINLIFE
               ENDIF
            ENDIF

            ! Add to overall reaction rate
            HET_HOBr_PSC = HET_HOBr_PSC + ADJUSTEDRATE

         ENDDO

      ENDIF

    END FUNCTION HETHOBr_PSC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_nat
!
! !DESCRIPTION: Subroutine CHECK\_NAT determines whether the solid PSC is 
!  composed of ice or NAT (nitric acid trihydrate) (needed for heterogeneous
!  chemistry), or indeed if there is any direct PSC calculation at all. This
!  is important for determining whether to use the JPP or Kirner scheme for
!  ice cloud radii.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_NAT( I, J, L, IS_NAT, IS_PSC, IS_STRAT, &
                            Input_Opt, State_Met, State_Chm )
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)  :: I,J,L      ! Grid indices
      TYPE(OptInput), INTENT(IN)  :: Input_Opt  ! Input options
      TYPE(MetState), INTENT(IN)  :: State_Met  ! Meteorology State object
      TYPE(ChmState), INTENT(IN)  :: State_Chm  ! Chemistry State object
!
! !OUTPUT VARIABLES:
!
      LOGICAL,        INTENT(OUT) :: IS_NAT     ! Is surface NAT?
      LOGICAL,        INTENT(OUT) :: IS_PSC     ! Are there solid PSCs?
      LOGICAL,        INTENT(OUT) :: IS_STRAT   ! Are we in the strat?
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms
!
! !REVISION HISTORY:
!  17 Apr 2013 - S. D. Eastham - Initial version
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm
!                                objects via the arg list
!  08 Apr 2015 - R. Yantosca   - Remove call to READ_PSC_FILE, this is
!                                now done from DO_CHEMISTRY (chemistry_mod.F)
!  28 Jan 2016 - M. Sulprizio  - Moved this routine from ucx_mod.F to
!                                gckpp_HetRates.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: IS_TROP

      !=================================================================
      ! CHECK_NAT begins here!
      !=================================================================

      ! Check if box is in the troposphere
      IS_TROP  = ( State_Met%PEDGE(I,J,L) > State_Met%TROPP(I,J) )

      ! Check if box is in the stratosphere
      IS_STRAT = ( ( L .le. LLSTRAT ) .and. ( .not. IS_TROP) )

      ! Check if there are solid PSCs
      IS_PSC   = ( ( Input_Opt%LPSCCHEM ) .and. &
                 ( State_Chm%STATE_PSC(I,J,L) >= 2.0 ) .and. ( IS_STRAT ) )

      ! Check if there is surface NAT
      IS_NAT   = ( ( IS_PSC ) .and. ( TRC_NIT .gt. TINY(1e+0_fp) ) )

    END SUBROUTINE CHECK_NAT
!EOC
#endif
  END MODULE GCKPP_HETRATES
