!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: aciduptake_DustChemFuncs
!
! !DESCRIPTION: Module containing rate-law functions for the dust acid
!  uptake species, in the aciduptake mechanism.
!\\
!\\
! !INTERFACE:

MODULE aciduptake_DustChemFuncs
!
! !USES:
!
  USE GcKpp_Precision

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: aciduptake_DustChem
  PUBLIC :: aciduptake_InitDustChem
!
! !PRIVATE DATA MEMBERS:
!
  INTEGER, PRIVATE :: id_DSTAL1, id_DSTAL2, id_DSTAL3, id_DSTAL4
  INTEGER, PRIVATE :: id_DST1,   id_DST2,   id_DST3,   id_DST4

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aciduptake_DustChem
!
! !DESCRIPTION: Computes the reaction rates [1/s] for acid uptake on dust
!  species. 
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE aciduptake_DustChem( I,         J,         L,                   &
                                  Input_Opt, State_Chm, State_Met, RC       )
!
! !USES:
!
    USE CMN_SIZE_Mod,     ONLY : NDSTBIN
    USE ErrCode_Mod
    USE GcKpp_Global,     ONLY : C, K_DST
    USE GcKpp_Parameters
    USE GcKpp_Precision
    USE Input_Opt_Mod,    ONLY : OptInput
    USE PhysConstants,    ONLY : AIRMW, AVO
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Met_Mod,    ONLY : MetState
    USE rateLawUtilFuncs, ONLY : KIIR1Ltd
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L     ! X, Y, Z grid indices
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: prtDebug
    INTEGER            :: N
    REAL(dp)           :: K0, Ki, F, M, KK, F1

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    REAL(dp)           :: ALK_d(NDSTBIN)
    REAL(dp)           :: KTS(NDSTBIN)
    REAL(dp)           :: KTN(NDSTBIN)
    REAL(dp)           :: KTH(NDSTBIN)

    !========================================================================
    ! aciduptake_DustChem begins here!
    !========================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at SET_CLD_S (in module GeosCore/sulfate_mod.F90)'

    ! Updated to match JPL 2006 + full chem (jaf, 10/14/09)

    !========================================================================
    ! Get dust alkalinity ALK_d (NDSTBIN) [v/v], Uptake rates for
    ! sulfate, KTS(NDSTBIN), and nitrate, KTN(NDSTBIN) on dust [s-1]
    !========================================================================
    CALL Get_Dust_Alk(  I         = I,                                       &
                        J         = J,                                       &
                        L         = L,                                       &
                        ALK_d     = ALK_d,                                   &
                        KTS       = KTS,                                     &
                        KTN       = KTN,                                     &
                        KTH       = KTH,                                     &
                        Input_Opt = Input_Opt,                               &
                        State_Met = State_Met,                               &
                        State_Chm = State_Chm                               )

    !========================================================================
    ! Get dust alkalinity ALK_d (NDSTBIN) [v/v], Uptake rates for
    ! sulfate = KTS(NDSTBIN) on dust [s-1]
    ! nitrate = KTN(NDSTBIN) on dust [s-1]
    ! H2SO4   = KTH(NDSTBIN) on dust [s-1]
    !========================================================================

    K_DST(1)  = KIIR1Ltd( C(ind_SO2),   2.0_dp*C(ind_DSTAL1), KTS(1) )
    K_DST(2)  = KIIR1Ltd( C(ind_SO2),   2.0_dp*C(ind_DSTAL2), KTS(2) )
    K_DST(3)  = KIIR1Ltd( C(ind_SO2),   2.0_dp*C(ind_DSTAL3), KTS(3) )
    K_DST(4)  = KIIR1Ltd( C(ind_SO2),   2.0_dp*C(ind_DSTAL4), KTS(4) ) 

    K_DST(5)  = KIIR1Ltd( C(ind_HNO3),  C(ind_DSTAL1),        KTN(1) )
    K_DST(6)  = KIIR1Ltd( C(ind_HNO3),  C(ind_DSTAL2),        KTN(2) )
    K_DST(7)  = KIIR1Ltd( C(ind_HNO3),  C(ind_DSTAL3),        KTN(3) )
    K_DST(8)  = KIIR1Ltd( C(ind_HNO3),  C(ind_DSTAL4),        KTN(4) )

    K_DST(9)  = KIIR1Ltd( C(ind_H2SO4), 2.0_dp*C(ind_DSTAL1), KTH(1) )
    K_DST(10) = KIIR1Ltd( C(ind_H2SO4), 2.0_dp*C(ind_DSTAL2), KTH(2) )
    K_DST(11) = KIIR1Ltd( C(ind_H2SO4), 2.0_dp*C(ind_DSTAL3), KTH(3) )
    K_DST(12) = KIIR1Ltd( C(ind_H2SO4), 2.0_dp*C(ind_DSTAL4), KTH(4) )

    !========================================================================
    ! Gas phase SO4 production is done here in offline run only
    !========================================================================
    K0        = 3.3e-31_dp * ( 300.0_dp / State_Met%T(I,J,L) )**4.3_dp
    Ki        = 1.6e-12_dp
    F         = 1000.0_dp / AIRMW * AVO * 1.e-6_dp
    M         = State_Met%AIRDEN(I,J,L) * F
    KK        = K0 * M / Ki
    F1        = 1.0_dp / ( 1.0_dp + ( LOG10( KK ) )**2 )
    K_DST(13) = ( K0 * M / ( 1.0_dp + KK ) ) * 0.6_dp**F1

  END SUBROUTINE aciduptake_DustChem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_dust_alk
!
! !DESCRIPTION: Subroutine GET\_DUST\_ALK returns: (1) dust alkalinity,
!  ALK\_d(NDSTBIN) [v/v], (2) rate coefficients, KTS(NDSTBIN), KTN(NDSTBIN),
!  for uptake of SO2 and HNO3 on dust for use in sulfate\_mod.f for chemistry
!  on dust aerosols, (3) fraction, KTH(NDSTBIN), of the size-weighted total
!  area of aerosols in the grid box. GET\_DUST\_ALK is analogous to GET\_ALK
!  for seasalt (bec, 12/7/04; tdf 04/08/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_DUST_ALK( I,   J,   L,         ALK_d,     KTS,              &
                           KTN, KTH, Input_Opt, State_Met, State_Chm        )
!
! !USES:
!
    USE CMN_SIZE_MOD,    ONLY : NDUST, NDSTBIN
    USE ERROR_MOD,       ONLY : IT_IS_NAN
    USE GcKpp_Parameters
    USE Input_Opt_Mod,   ONLY : OptInput
    USE PhysConstants,   ONLY : PI, AIRMW
    USE Species_Mod,     ONLY : SpcConc
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Met_Mod,   ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L        ! Grid box indices
    TYPE(OptInput), INTENT(IN)    :: Input_Opt      ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met      ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm      ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    REAL(dp),       INTENT(OUT)   :: ALK_d(NDSTBIN) ! Dust alkalinity [v/v]
    REAL(dp),       INTENT(OUT)   :: KTS  (NDSTBIN) ! Rate coef for uptake of
                                                    ! SO2 on dust [s-1]
    REAL(dp),       INTENT(OUT)   :: KTN  (NDSTBIN) ! Rate coef for uptake of
                                                    ! HNO3 on dust [s-1]
    REAL(dp),       INTENT(OUT)   :: KTH  (NDSTBIN) ! Fraction of the size-
                                                    ! weighted total area
                                                    ! of aerosols in grid box
!
! !REVISION HISTORY:
!  08 Apr 2008 - T.D. Fairlie- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(dp), PARAMETER :: MINDAT = 1.d-20

    ! Need this for dust
    !REAL(dp), PARAMETER :: GAMMA_SO2 = 0.05d0  !(Song & Carmichael, 2001)
    ! 200 times smaller 8/28/2K9
    REAL(dp), PARAMETER :: GAMMA_SO2 = 2.5d-4

    !tdf V9 4/1/2K9 Applying Song et al.(2007) reduced value
    REAL(dp), PARAMETER :: GAMMA_H2SO4 = 1.d0

    !  Need this for dust
    !REAL(dp), PARAMETER :: GAMMA_HNO3 = 0.1d0 ! (Song & Carmichael, 2001)
    ! 200 times smaller 8/28/2K9
    REAL(dp), PARAMETER :: GAMMA_HNO3 = 5.0d-4

    REAL(dp), PARAMETER :: DG = 0.2d0 ! gas phase diffusion coeff. [cm2/s]
    REAL(dp), PARAMETER :: v = 3.0d4  ! cm/s
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: IRH
    INTEGER             :: IBIN, ISBIN
    REAL(dp)            :: N1, KT1, KT1N
    REAL(dp)            :: AREA, HGF
    REAL(dp)            :: CONST1, CONST2, CONST
    REAL(dp)            :: AIR_DENS
    REAL(dp)            :: A1, B1, A1N, B1N, T1, R1
    REAL(dp)            :: DD, RD, DM
    REAL(dp)            :: TOTAL_AREA, DF_TOTAL_AREA
    REAL(dp)            :: DST_d (NDSTBIN), ALK, GAMS, GAMN
    REAL(dp)            :: SULF_AREA, BC_AREA, OC_AREA
    REAL(dp)            :: SULF_RAD, BC_RAD, OC_RAD
    REAL(dp)            :: SULF_FAC, BC_FAC, OC_FAC
    REAL(dp)            :: SSA_AREA, SSC_AREA
    REAL(dp)            :: SSA_RAD, SSC_RAD
    REAL(dp)            :: SSA_FAC, SSC_FAC
    REAL(dp)            :: term
    LOGICAL, SAVE       :: FIRST = .TRUE.

    ! Arrays
    ! Dust Surface Areas                                ! tdf 08/20/09
    REAL(dp)            :: AREA_d(NDSTBIN)              ! [cm^2/cm^3]

    ! Dust Surface Areas within sub-bins 1-4 of BIN 1   ! tdf 08/20/09
    REAL(dp)            :: AREA_sd1(4)                  ! [cm^2/cm^3]

    ! Dust Effective Radius                             ! tdf 08/20/09
    REAL(dp)            :: RD_d(NDSTBIN)                ! [cm]

    ! Dust Effective Radii for sub-bins 1-4 of BIN 1    ! tdf 08/20/09
    REAL(dp)            :: RD_sd1(4)                    ! [cm]

    ! Dust size-weighted Surface Areas                  ! tdf 08/20/09
    REAL(dp)            :: DF_AREA_d(NDSTBIN)           ! [1/s]

    ! Dust size-weighted Surface Areas for sub-bins 1-4 ! tdf 08/20/09
    REAL(dp)            :: DF_AREA_sd1(4)               ! [1/s]

    ! Molecular weights
    REAL(dp)            :: MW_DST1, MW_DST2, MW_DST3, MW_DST4

    ! Pointers
    TYPE(SpcConc), POINTER   :: Spc(:)
    REAL(dp),      POINTER   :: ERADIUS(:,:,:,:)
    REAL(dp),      POINTER   :: TAREA(:,:,:,:)

    !=================================================================
    ! GET_DUST_ALK begins here!
    !=================================================================

    ! Initialize pointers
    Spc       => State_Chm%Species  ! GEOS-Chem species array [v/v dry]
    ERADIUS   => State_Chm%AeroRadi    ! Aerosol Radius [cm]
    TAREA     => State_Chm%AeroArea    ! Aerosol Area [cm2/cm3]

    ! Get MWs from species database
    MW_DST1   = State_Chm%SpcData(id_DST1)%Info%MW_g
    MW_DST2   = State_Chm%SpcData(id_DST2)%Info%MW_g
    MW_DST3   = State_Chm%SpcData(id_DST3)%Info%MW_g
    MW_DST4   = State_Chm%SpcData(id_DST4)%Info%MW_g

    ! Zero variables
    ALK_d     = 0.0_dp
    KTS       = 0.0_dp
    KTN       = 0.0_dp
    KTH       = 0.0_dp
    AREA_d    = 0.0_dp
    RD_d      = 0.0_dp

    ! Air density [kg/m3]
    AIR_DENS  = State_Met%AD(I,J,L) / State_Met%AIRVOL(I,J,L)

    ! Retrieve Dust Alkalinity [v/v dry from Spc array
    ALK_d(1)  = Spc(id_DSTAL1)%Conc(I,J,L)
    ALK_d(2)  = Spc(id_DSTAL2)%Conc(I,J,L)
    ALK_d(3)  = Spc(id_DSTAL3)%Conc(I,J,L)
    ALK_d(4)  = Spc(id_DSTAL4)%Conc(I,J,L)

    ! Dust [kg/m3] from Spc, used to compute dust surface area
    ! Units: (moles/mole).(kg(air)/m3).(kg(dust)/mole)/(kg(air)/mole)
    DST_d(1)  = Spc(id_DST1)%Conc(I,J,L) * AIR_DENS / ( AIRMW / MW_DST1 )
    DST_d(2)  = Spc(id_DST2)%Conc(I,J,L) * AIR_DENS / ( AIRMW / MW_DST2 )
    DST_d(3)  = Spc(id_DST3)%Conc(I,J,L) * AIR_DENS / ( AIRMW / MW_DST3 )
    DST_d(4)  = Spc(id_DST4)%Conc(I,J,L) * AIR_DENS / ( AIRMW / MW_DST4 )

    ! tdf Now get aerosol surface area from TAREA (cm2/cm3)
    SULF_AREA = TAREA(I,J,L,NDUST+1)
    BC_AREA   = TAREA(I,J,L,NDUST+2)
    OC_AREA   = TAREA(I,J,L,NDUST+3)
    SSA_AREA  = TAREA(I,J,L,NDUST+4)
    SSC_AREA  = TAREA(I,J,L,NDUST+5)

    ! tdf Now get aerosol effective radius from ERADIUS (cm)
    SULF_RAD  = ERADIUS(I,J,L,NDUST+1)
    BC_RAD    = ERADIUS(I,J,L,NDUST+2)
    OC_RAD    = ERADIUS(I,J,L,NDUST+3)
    SSA_RAD   = ERADIUS(I,J,L,NDUST+4)
    SSC_RAD   = ERADIUS(I,J,L,NDUST+5)

    ! tdf Quotients [s/cm] used to weight surface area for H2SO4 uptake
    term      = 4.0_dp / ( V * GAMMA_H2SO4 ) 
    SULF_FAC  = ( SULF_RAD / DG + term )
    BC_FAC    = (   BC_RAD / DG + term )
    OC_FAC    = (   OC_RAD / DG + term )
    SSA_FAC   = (  SSA_RAD / DG + term )
    SSC_FAC   = (  SSC_RAD / DG + term )

    !tdf Surface areas and effective radii for sub-bins 1-4 of dust bin 1
    DO ISBIN = 1, 4
       T1 = TAREA  (I,J,L,ISBIN)
       R1 = ERADIUS(I,J,L,ISBIN)
       AREA_sd1    (ISBIN) = T1
       RD_sd1      (ISBIN) = R1
       !tdf surface area for sub bins 1-4 in bin 1, weighted by gas-phase
       !tdf diffusion and collision limitations
       !tdf used to compute proportionate uptake of H2SO4 only  [1/s]
       DF_AREA_sd1 (ISBIN) = T1 / (R1/DG + 4.0e+0_dp/(V*GAMMA_H2SO4))
    END DO

    !-----------------------------------------------------------------------
    ! Very Simple Formulation: For each size bin (i)   ! tdf 8/20/09
    ! Dust Area density = 3 * Dust Mass density  / (REFF(i) * DUSTDEN)
    ! TAREA computed   in RDUST_ONLINE - Units: cm^2(dust) / cm^3(air)
    ! ERADIUS computed in RDUST_ONLINE - Units: cm
    ! NB: I am now subdividing the submicron dust size bin
    !     using TAREA (I,J,L,1->4), and ERADIUS (I,J,L,1->4).
    !-----------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !  Find Dust surface area density in grid-box, AREA_d [cm^2/cm^3].
    !  Also find the size-weighted surface area density, DF_AREA_d [1/s].
    !  The latter represents the gas-phase diffusion and surface
    !  limited weighting and is used to determine the fraction of H2SO4
    !  taken up on dust, versus taken up on other aerosols.
    !                                                 tdf 08/21/09
    !-------------------------------------------------------------------------

    ! tdf Loop over size bins  (NDSTBIN = 4)
    DO IBIN = 1, NDSTBIN

       ! Dust Area density in grid box,      AREA_d [cm^2/cm^3]    tdf 8/21/09
       ! Dust weighted surface area density, DF_AREA_d [1/s]       tdf 8/21/09

       IF (IBIN .EQ. 1) THEN
          ! For Dust size bin 1, sum over the 4 size sub bins  tdf 8/21/09
          AREA_d   (IBIN) = AREA_sd1(1) + AREA_sd1(2) &       ![cm^2/cm^3]
                          + AREA_sd1(3) + AREA_sd1(4)
          DF_AREA_d(IBIN) = DF_AREA_sd1(1) + DF_AREA_sd1(3) &  ! [1/s]
                          + DF_AREA_sd1(2) + DF_AREA_sd1(4)
       ELSE
          T1 = TAREA(I,J,L,3+IBIN)      ! [cm^2/cm^3]
          R1 = ERADIUS(I,J,L,3+IBIN)    ! [cm]
          RD_d     (IBIN) = R1
          AREA_d   (IBIN) = T1          ! [cm^2/cm^3]
          DF_AREA_d(IBIN) = T1 / (R1/DG + 4.0D0/(V*GAMMA_H2SO4)) ! [1/s]
       ENDIF

    END DO

    ! tdf total aerosol surface area  [cm^2/cm^3]
    TOTAL_AREA = SULF_AREA + BC_AREA + OC_AREA + SSA_AREA  + SSC_AREA + &
                 AREA_d(1) + AREA_d(2) + AREA_d(3) + AREA_d(4)

    ! tdf total surface area weighted by gas-phase diffusion limitation [1/s]
    DF_TOTAL_AREA = SULF_AREA / SULF_FAC + &
                    BC_AREA   / BC_FAC   + &
                    OC_AREA   / OC_FAC   + &
                    SSA_AREA  / SSA_FAC  + &
                    SSC_AREA  / SSC_FAC  + &
                    DF_AREA_d(1)         + &
                    DF_AREA_d(2)         + &
                    DF_AREA_d(3)         + &
                    DF_AREA_d(4)

    ! tdf Total Dust Alkalinity
    ALK = ALK_d(1) + ALK_d(2) + ALK_d(3) + ALK_d(4)  ! [v/v]

    ! set humidity index IRH as a percent
    IRH = State_Met%RH(I,J,L)
    IRH = MAX(  1, IRH )
    IRH = MIN( 99, IRH )

    ! hygroscopic growth factor for dust: Set to NO GROWTH for now
    IF ( IRH < 100 ) HGF = 1.0e+0_dp

    ! tdf Loop over size bins (NDSTBIN = 4)
    DO IBIN = 1, NDSTBIN

       !----------------------------------
       ! SO2 uptake onto particles
       !----------------------------------

       !tdf 2/11/2K9
       !tdf Following relative uptake rates of Preszler-Prince et al.(2007)
       IF ( IRH >= 90.0_dp ) THEN
          GAMS = GAMMA_SO2 * 2.0_dp
       ELSE IF ( IRH >= 84.0_dp ) THEN
          GAMS = GAMMA_SO2                                                   &
               * ( 0.5_dp  + 1.5_dp*(IRH - 84.0_dp)  / (90.0_dp - 84.0_dp)  )
       ELSE IF ( IRH >= 76.0_dp ) THEN
          GAMS = GAMMA_SO2                                                   &
               * ( 0.16_dp + 0.34_dp*(IRH - 76.0_dp) / (84.0_dp - 76.0_dp)  )
       ELSE IF ( IRH >= 33.0_dp ) THEN
          GAMS = GAMMA_SO2                                                   &
               * ( 0.03_dp + 0.13_dp*(IRH-33.e+0_dp) / (76.0_dp - 33.0_dp)  )
       ELSE IF ( IRH >= 20.0_dp ) THEN
          GAMS = GAMMA_SO2 * 0.03_dp
       ELSE                                       ! 0.0 below 20%
          GAMS = 0.0_dp
       ENDIF

       ! Check for sufficient alkalinity          tdf 3/28/2K8
       IF ( ALK > MINDAT ) THEN

          ! calculate gas-to-particle rate constant for uptake of
          ! SO2 onto dust aerosols [Jacob, 2000] analytical solution
          ! Corrected based on discussions with Becky     tdf 07/14/08
          KT1    = 0.0_dp

          IF (IBIN .EQ. 1) THEN

             ! tdf Sum over the 1-4 sub-bins for bin 1      ! tdf 08/21/2K9
             DO ISBIN = 1, 4
                RD     = RD_sd1 (ISBIN)        ! effective radius [cm]
                AREA   = AREA_sd1 (ISBIN)      ! Dust Surface Area [cm^2/cm^3]

                ! Prevent divide by zero if GAMS = 0 (tdf, mps, 11/14/13)
                IF ( GAMS > 0.0_dp ) THEN
                   CONST1 = 4.0_dp/(V*GAMS)    ! Collision [s/cm]
                   CONST2 = RD/DG              ! Diffusion [s/cm]
                   CONST  = CONST1 + CONST2
                   KT1    = KT1 + AREA / CONST ! [cm^2/cm^3] * [cm/s] = [1/s]
                ELSE
                   KT1    = KT1                ! [cm^2/cm^3] * [cm/s] = [1/s]
                ENDIF
             END DO

          ELSE

             RD     = RD_d (IBIN)              ! effective radius [cm]
             AREA   = AREA_d (IBIN)            ! Dust Surface Area [cm^2/cm^3]

             ! Prevent divide by zero if GAMS = 0 (tdf, mps, 11/14/13)
             IF ( GAMS > 0.0_dp ) THEN
                CONST1 = 4.0_dp/(V*GAMS)       ! Collision [s/cm]
                CONST2 = RD/DG                 ! Diffusion [s/cm]
                CONST  = CONST1 + CONST2
                KT1    = AREA / CONST          ! [cm^2/cm^3] * [cm/s] = [1/s]
             ELSE
                KT1    = 0.0_dp                ! [cm^2/cm^3] * [cm/s] = [1/s]
             ENDIF

          ENDIF

          KTS(IBIN) = KT1

       ELSE

          ! If no alkalinity, set rate coefficients to zero
          !tdf
          KTS(IBIN)  = 0.0_dp

       ENDIF

       !----------------------------------
       ! HNO3 uptake onto particles
       !----------------------------------

       !tdf 2/11/2K9
       !tdf Following uptake coefficients of Liu et al.(2007)
       IF (IRH >= 80.0_dp ) THEN
          GAMN = GAMMA_HNO3 * 2.1_dp
       ELSE IF (IRH >= 70.0_dp ) THEN
          GAMN = GAMMA_HNO3                                                  &
               * ( 1.3_dp  + 0.7_dp   * (IRH - 70.0_dp) / 10.0_dp )
       ELSE IF (IRH >= 60.0_dp ) THEN
          GAMN = GAMMA_HNO3                                                  &
               * ( 1.0_dp  + 0.3_dp   * (IRH - 60.0_dp) / 10.0_dp )
       ELSE IF ( IRH >= 50.0_dp ) THEN
          GAMN = GAMMA_HNO3                                                  &
               * ( 0.7_dp  + 0.3_dp   * (IRH - 50.0_dp) / 10.0_dp )
       ELSE IF ( IRH >= 30.0_dp ) THEN
          GAMN = GAMMA_HNO3                                                  &
               * ( 0.19_dp + 0.255_dp * (IRH - 30.0_dp) / 10.0_dp )
       ELSE IF ( IRH >= 10.0_dp ) THEN
          GAMN = GAMMA_HNO3                                                  &
               * ( 0.03_dp + 0.08_dp  * (IRH - 10.0_dp) / 10.0_dp )
       ELSE
          ! 0.0 below 10%
          GAMN = 0.0_dp
       ENDIF

       ! Check for sufficient alkalinity      tdf 3/28/2K8
       IF ( ALK > MINDAT ) THEN

          ! calculate gas-to-particle rate constant for uptake of
          ! HNO3 onto dust aerosols [Jacob, 2000] analytical solution
          ! Corrected based on discussions with Becky     tdf 07/14/08
          KT1    = 0.0e+0_dp

          IF (IBIN .EQ. 1) THEN

             ! tdf Sum over the 1-4 sub-bins for bin 1      ! tdf 08/21/2K9
             DO ISBIN = 1, 4
                RD = RD_sd1 (ISBIN)            ! effective radius [cm]
                AREA = AREA_sd1 (ISBIN)        ! Dust Surface Area [cm^2/cm^3]

                ! Prevent divide by zero if GAMN = 0 (tdf, mps, 11/14/13)
                IF ( GAMN > 0.0_dp ) THEN
                   CONST1 = 4.0_dp/(V*GAMN)    ! Collision [s/cm]
                   CONST2 = RD/DG              ! Diffusion [s/cm]
                   CONST  = CONST1 + CONST2
                   KT1    = KT1 + AREA / CONST ! [cm^2/cm^3] * [cm/s] = [1/s]
                ELSE
                   KT1    = KT1                ! [cm^2/cm^3] * [cm/s] = [1/s]
                ENDIF
             END DO

          ELSE

             RD     = RD_d (IBIN)              ! effective radius [cm]
             AREA   = AREA_d (IBIN)            ! Dust Surface Area [cm^2/cm^3]

             ! Prevent divide by zero if GAMN = 0 (tdf, mps, 11/14/13)
             IF ( GAMN > 0.0_dp ) THEN
                CONST1 = 4.0_dp/(V*GAMN)       ! Collision [s/cm]
                CONST2 = RD/DG                 ! Diffusion [s/cm]
                CONST  = CONST1 + CONST2
                KT1    = AREA / CONST          ! [cm^2/cm^3] * [cm/s] = [1/s]
             ELSE
                KT1    = 0.0_dp               ! [cm^2/cm^3] * [cm/s] = [1/s]
             ENDIF

          ENDIF

          KTN(IBIN) = KT1

       ELSE

          ! If no alkalinity, set rate coefficients to zero
          !tdf
          KTN(IBIN)  = 0.0_dp

       ENDIF

       !----------------------------------
       ! H2SO4 uptake onto particles
       !----------------------------------

       ! Uptake not limited by dust alkalinity      tdf 3/02/2K9

       !tdf As of 08/20/09, we use AREA and size weighted uptake
       !tdf where now KTH is a fractional uptake for each size bin
       !tdf with respect to total aerosol surface area.

       KT1    = DF_AREA_d(IBIN) / DF_TOTAL_AREA   ! Fraction

       KTH(IBIN) = KT1

    END DO ! tdf End Loop over size bins

    ! Free pointers
    Spc     => NULL()
    ERADIUS => NULL()
    TAREA   => NULL()

  END SUBROUTINE Get_Dust_Alk
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aciduptake_InitDustChem
!
! !DESCRIPTION: Defines species indices used for the dust acid uptake 
!  rate-law functions.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE aciduptake_InitDustChem( RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE State_Chm_Mod, ONLY : Ind_
!
! !INPUT/OUTPUT PARAMETERS: 
!
    INTEGER, INTENT(INOUT) :: RC   ! Success or failure
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Initialize
    RC        = GC_SUCCESS

    ! Dust alkalinity species
    id_DSTAL1 = Ind_( 'DSTAL1' )
    id_DSTAL2 = Ind_( 'DSTAL2' )
    id_DSTAL3 = Ind_( 'DSTAL3' )
    id_DSTAL4 = Ind_( 'DSTAL4' )

    ! Dust species
    id_DST1   = Ind_( 'DST1'   )
    id_DST2   = Ind_( 'DST2'   )    
    id_DST3   = Ind_( 'DST3'   ) 
    id_DST4   = Ind_( 'DST4'   ) 

  END SUBROUTINE aciduptake_InitDustChem
!EOC
END MODULE aciduptake_DustChemFuncs
