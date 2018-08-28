!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ras_mod.F90
!
! !DESCRIPTION: Computes convective mass flux and detrainment for
! offline GEOS-Chem simulations at the resolution of the simulation
! (rather than the resolution of the met). 
!\\
!\\
! !INTERFACE:
!
MODULE RAS_Mod
!
! !USES:
!
  USE Errcode_Mod
  USE PhysConstants
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Do_RAS
!
! !PRIVATE MEMBER FUNCTIONS:
! 
  PRIVATE :: RasE
  PRIVATE :: Qsat
  PRIVATE :: Rearrange_Grid
  PRIVATE :: Rearrange_Grid_Int
  PRIVATE :: Rearrange_2d
! 
! !DEFINED PARAMETERS:
!
  ! Parameters for use below (private to this module)
  REAL(fp), PARAMETER :: MAPL_GRAV    = 9.80665_fp                 !m^2/s
  REAL(fp), PARAMETER :: MAPL_RUNIV   = 8314.47_fp                 ! J/Kmol K
  REAL(fp), PARAMETER :: MAPL_AIRMW   = 28.965_fp                  ! kg/Kmol
  REAL(fp), PARAMETER :: MAPL_H2OMW   = 18.015_fp
  REAL(fp), PARAMETER :: MAPL_RDRY    = MAPL_RUNIV / MAPL_AIRMW
  REAL(fp), PARAMETER :: MAPL_CPDRY   = 3.5_fp     * MAPL_RDRY
  REAL(fp), PARAMETER :: MAPL_KAPPA   = MAPL_RDRY  / MAPL_CPDRY
  REAL(fp), PARAMETER :: MAPL_RGAS    = MAPL_RDRY
  REAL(fp), PARAMETER :: MAPL_CP      = MAPL_RGAS  / MAPL_KAPPA
  REAL(fp), PARAMETER :: MAPL_ALHL    = 2.4665e6_fp                ! J/kg @15C
  REAL(fp), PARAMETER :: MAPL_EPSILON = MAPL_H2OMW / MAPL_AIRMW
  REAL(fp), PARAMETER :: MAPL_VIREPS  = 1.0_fp / MAPL_EPSILON - 1.0_fp
!
! !REMARKS:
! 
! !REVISION HISTORY:
!  23 Feb 2017 - K. Yu       - Initial version taken from strippedras.F90
!                              and qsat.F90 from Andrea Molod. 
!  27 Aug 2018 - R. Yantosca - Add updates for GEOS-Chem 12.  Move the
!                              RASPARAMS type into State_Met.
!EOP
!-----------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
! 
! !IROUTINE: Do_RAS
! 
! !DESCRTION: Wrapper for RASE (the Relaxed Arakawa-Schubert convection code).
!  Computes the convection fields at the resolution of the simulation and
!  then stores those in State\_Met.
!\\
!\\
!
  SUBROUTINE Do_RAS( am_I_Root, Input_Opt,  State_Met,                       &
                     State_Chm, State_Diag, RC                              )
!
! !USES:
! 
    USE CMN_SIZE_Mod
    USE Physconsts
    USE ErrCode_Mod
    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Met_Mod,   ONLY : MetState
    USE State_Diag_Mod,  ONLY : DgnState
    USE Time_Mod,        ONLY : GET_TS_DYN
    USE Pressure_Mod,    ONLY : GET_AP,   GET_BP
    USE Grid_Mod,        ONLY : GET_XMID, GET_YMID
#if defined( BPCH_DIAG )
    USE CMN_DIAG_MOD
    USE DIAG_MOD,        ONLY : AD66
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
!
! !INPUT/OUTPUT PARAMETERS:
! 
    
    TYPE(MetState), INTENT(INOUT) :: State_Met
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
    TYPE(DgnState), INTENT(INOUT) :: State_Diag
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  23 Feb 2017 - K. Yu       - Initial version
!  27 Aug 2018 - R. Yantosca - RasParams is now folded into the State_Met obj
!EOP
!-----------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL         :: prtDebug
    INTEGER         :: IDIM, IRUN
    INTEGER         :: CBL_METHOD
    INTEGER         :: I, J, L
    INTEGER         :: ICMIN
    REAL(fp)        :: PMIN_DET, PMIN_CBL
    REAL(fp)        :: AUTOC_CN_OCN, AUTOC_CN_LAND, AUTOC_CN_ZDEP
    INTEGER         :: KCBLMIN
    REAL(fp)        :: CBL_QPERT, CBL_TPERT, CBL_TPERT_MXOCN, CBL_TPERT_MXLND
    REAL(fp)        :: CNV_Q600_MIN, CNV_Q600_MAX
    REAL(fp)        :: DT, LON, part_tmp
    REAL(fp)        :: COLUMN_RAS, COLUMN_FP, REPARTITION
    REAL(fp)        :: REPARTITION_DQR, REPARTITION_PFI
   
    ! Arrays
    REAL(fp)        :: PREF(LLPAR+1)
    INTEGER         :: TMP_INT(IIPAR, JJPAR)
    REAL(fp)        :: TMP_REAL(IIPAR, JJPAR)
    INTEGER         :: KPBL(IIPAR, JJPAR)
    REAL(fp)        :: EDGE_PRES(IIPAR, JJPAR)
    REAL(fp)        :: QV600(IIPAR, JJPAR)
    REAL(fp)        :: TH(IIPAR,JJPAR,LLPAR)
    REAL(fp)        :: CMFMC(IIPAR, JJPAR, LLPAR+1)
    REAL(fp)        :: DTRAIN(IIPAR, JJPAR, LLPAR)
    REAL(fp)        :: PART(IIPAR,JJPAR)
    REAL(fp)        :: CONV_FRAC(IIPAR,JJPAR)

    ! reshaped arrays
    REAL(fp)        :: TEMPERATURE(IIPAR*JJPAR, LLPAR)
    REAL(fp)        :: SURF_TEMP(IIPAR*JJPAR)
    REAL(fp)        :: FRLAND(IIPAR*JJPAR)
    INTEGER         :: SEEDRAS(IIPAR*JJPAR, 2)
    INTEGER         :: IRAS(IIPAR*JJPAR)
    INTEGER         :: JRAS(IIPAR*JJPAR)
    REAL(fp)        :: SIGE(LLPAR+1)
    REAL(fp)        :: PK(IIPAR*JJPAR,LLPAR)
    REAL(fp)        :: PLO(IIPAR*JJPAR, LLPAR)
    REAL(fp)        :: ZLO(IIPAR*JJPAR, LLPAR)
    REAL(fp)        :: ZLE(IIPAR*JJPAR, LLPAR+1)  
    REAL(fp)        :: ZCBLx(IIPAR*JJPAR)
    REAL(fp)        :: WGT0(IIPAR*JJPAR, LLPAR)
    REAL(fp)        :: WGT1(IIPAR*JJPAR, LLPAR)
    INTEGER         :: KCBL(IIPAR*JJPAR)
    REAL(fp)        :: CNV_PLE(IIPAR*JJPAR,LLPAR+1)
    REAL(fp)        :: TPERT(IIPAR*JJPAR)
    REAL(fp)        :: QPERT(IIPAR*JJPAR)
    REAL(fp)        :: QSSFC(IIPAR*JJPAR)
    REAL(fp)        :: PKE(IIPAR*JJPAR, LLPAR+1) !exner fcn at edge
    REAL(fp)        :: QSS(IIPAR*JJPAR, LLPAR)  !saturation spec humid
    REAL(fp)        :: DQS(IIPAR*JJPAR, LLPAR) !d(qss)/dt
    REAL(fp)        :: CNV_FRAC(IIPAR*JJPAR)
    REAL(fp)        :: THO(IIPAR*JJPAR, LLPAR)  ! potentail temp
    REAL(fp)        :: QHO(IIPAR*JJPAR, LLPAR) ! specific humidity
    REAL(fp)        :: CNV_FLX(IIPAR*JJPAR, LLPAR)
    REAL(fp)        :: CNV_FLXD(IIPAR*JJPAR, LLPAR)
    REAL(fp)        :: CNV_FLXC(IIPAR*JJPAR, LLPAR+1)
    REAL(fp)        :: RASAL2_2d(IIPAR*JJPAR)
    REAL(fp)        :: MXDIAMx(IIPAR*JJPAR)
    INTEGER         :: IRC(IIPAR*JJPAR, LLPAR)
    
    ! Pointers 
    REAL(fp),       POINTER :: TEMP(:,:,:) 
    REAL(fp),       POINTER :: PMID_DRY(:,:,:)
    REAL(fp),       POINTER :: PLE(:,:,:)
    REAL(fp),       POINTER :: Q(:,:,:)
    REAL(fp),       POINTER :: DELP(:,:,:)
    TYPE(ParamRas), POINTER :: RasParams

    ! Strings

    !=======================================================================
    ! Do_RAS begins here!
    !=======================================================================

    ! Initialize
    RC        =  GIGC_SUCCESS
    prtDebug  =  ( am_I_Root .and. Input_Opt%LPRT )
    ErrMsg    = ''
    ThisLoc   = ' -> at Do_Ras (in module GeosCore/ras_mod.F90'
    IDIM      =  IIPAR * JJPAR
    IRUN      =  IIPAR * JJPAR
    DT        =  GET_TS_DYN() !* 60.0_fp  !### Timesteps are now in seconds!

    ! Point to fields of State_Met (flip in vertical where necessary)
    RasParams => State_Met%RasParams
    TEMP      => State_Met%T       (:,:,LLPAR:1:-1  )
    PMID_DRY  => State_Met%PMID_DRY(:,:,LLPAR:1:-1  )
    PLE       => State_Met%PEDGE   (:,:,LLPAR+1:1:-1)
    Q         => State_Met%SPHU    (:,:,LLPAR:1:-1  )
    DELP      => State_Met%DELP    (:,:,LLPAR:1:-1  )

    !-----------------------------------------------------------------------
    ! Compute ICMIN and associated values
    !-----------------------------------------------------------------------
    
    PMIN_DET      = 3000.0_fp
    PMIN_CBL      = 50000.0_fp
    AUTOC_CN_OCN  = 2.5e-3_fp
    AUTOC_CN_LAND = AUTOC_CN_OCN

    DO L = 1, LLPAR+1
       PREF(L) = (GET_AP(LLPAR-L+2) + (GET_BP(LLPAR-L+2) * 1000.0)) * 100.0
    ENDDO

    ! USE PBL Height as subcloud layer for RAS
    ICMIN = MAX( 1, COUNT( PREF < PMIN_DET ) )

    !-----------------------------------------------------------------------
    ! SEEDRAS: set random see for RAS
    !-----------------------------------------------------------------------
    TMP_INT = 0
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       TMP_INT(I,J) = 1000000 * ( 100*TEMP(I,J,LLPAR)   -                   &
                             INT( 100*TEMP(I,J,LLPAR) ) ) 
    ENDDO
    ENDDO

    ! Rearrange to 1-d
    CALL REARRANGE_GRID_INT(TMP_INT, SEEDRAS(:,1))

    TMP_INT = 0
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       TMP_INT(I,J) = 1000000 * ( 100*TEMP(I,J,LLPAR-1) -                   &
                             INT( 100*TEMP(I,J,LLPAR-1) ) ) 
    ENDDO
    ENDDO

    ! Rearrange to 1-d
    CALL REARRANGE_GRID_INT(TMP_INT, SEEDRAS(:,2))
      
    !-----------------------------------------------------------------------
    ! Compute IRAS, JRAS, and SIGE
    !----------------------------------------------------------------------
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       LON = GET_XMID(I,J,1)
       IF ( LON < 0.0_fp ) LON = 360.0_fp + LON
       IRAS((J-1)*IIPAR + I) = NINT( LON / 360.0_fp * 2.0_fp * PI * 100.0_fp)
       JRAS((J-1)*IIPAR + I) = NINT( GET_YMID(I,J,1) / 360.0_fp * 2.0_fp * PI * 100.0_fp)
    ENDDO
    ENDDO

    DO L = 1, LLPAR+1
       SIGE(L) = PREF(L) / PREF(LLPAR+1)
    ENDDO

    !-----------------------------------------------------------------------
    ! Compute KCBL
    !-----------------------------------------------------------------------
    CBL_METHOD = 6
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       IF(State_Met%PBL_TOP_L(I,J).ne.0) THEN
          KPBL(I,J) = MIN(LLPAR-State_Met%PBL_TOP_L(i,j)+1,LLPAR-1)
       ELSE
          KPBL(I,J) = LLPAR-1
       ENDIF
    ENDDO
    ENDDO

    KCBLMIN  = COUNT( PREF < PMIN_CBL ) 

    DO J=1,JJPAR
    DO I=1,IIPAR
       TMP_INT(I,J) = KPBL(I,J)
       TMP_INT(I,J) = MAX( TMP_INT (I,J), KCBLMIN )
    ENDDO
    ENDDO

    CALL REARRANGE_GRID_INT(TMP_INT, KCBL)

    DO I = 1, IDIM
       WGT0(I,:)              = 0.0_fp
       WGT0(I, KCBL(I):LLPAR) = 1.0_fp
       WGT1(I,:)              = 0.0_fp
       WGT1(I, KCBL(I):LLPAR) = 1.0_fp
    ENDDO

    !-----------------------------------------------------------------------
    ! Compute ZLE and ZLO
    !-----------------------------------------------------------------------

    ! Compute potential temperature
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
      TH(I,J,L) = TEMP(I,J,L) *                                             &
                  (State_Met%PSC2(I,J) / PMID_DRY(I,J,L)) ** MAPL_KAPPA
    ENDDO
    ENDDO
    ENDDO

    ! Rearrange into 2d grid
    DO L = 1, LLPAR
       CALL REARRANGE_GRID( TH(:,:,L), THO(:,L) )
    ENDDO

    DO L = 1, LLPAR+1
       CALL REARRANGE_GRID( PLE(:,:,L), CNV_PLE(:,L) )
    ENDDO

    PLO = 0.5_fp* ( CNV_PLE(:,1:LLPAR) + CNV_PLE(:,2:LLPAR+1) )
    PK  = ( PLO     / 1000.0_fp ) ** ( MAPL_RGAS / MAPL_CP )
    PKE = ( CNV_PLE / 1000.0_fp ) ** ( MAPL_RGAS / MAPL_CP )

    !### Debug
    IF ( prtDebug ) THEN
       PRINT*, 'PLO min max mean',                                          &
                minval(PLO), maxval(PLO),sum(PLO)/(IDIM*LLPAR)
       PRINT*, 'PKE min max mean',                                          &
                minval(PKE), maxval(PKE),sum(PKE)/(IDIM*LLPAR+IDIM)
    ENDIF

    DO L = 1, LLPAR
      CALL REARRANGE_GRID( Q(:,:,L), QHO(:,L) )
    ENDDO
    QHO = QHO * 0.001_fp    ! convert from g/kg to kg/kg

    ZLE(:,LLPAR+1) = 0.0_fp

    DO L = LLPAR,1,-1
    DO I = 1, IDIM
       ZLE(I,L) = THO(I,L)   * ( 1.0_fp + MAPL_VIREPS * QHO(I,L) )

       ZLO(I,L) = ZLE(I,L+1) + ( MAPL_CP    / MAPL_GRAV )                    &
                             * ( PKE(I,L+1) - PK(I,L)   )                    &
                             * ZLE(I,L)

       ZLE(I,L) = ZLO(I,L)   + ( MAPL_CP  / MAPL_GRAV )                      &
                             * ( PK (I,L) - PKE(I,L)  )                      &
                             * ZLE(I,L)
    ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    ! Compute ZCBL
    !-----------------------------------------------------------------------
    DO I = 1, IDIM
       ZCBLx(I) = ZLE(I, KCBL(I))
    ENDDO

    !### Debug
    IF ( prtDebug ) THEN
       print*, 'ZLE 1 min max mean', minval(ZLE(:,LLPAR)),                   &
                                     maxval(ZLE(:,LLPAR)),                   &
                                     sum(ZLE(:,LLPAR))/IDIM

       print*, 'ZLE 50min max mean', minval(ZLE(:,50)),                      &
                                     maxval(ZLE(:,50)),                      &
                                     sum(ZLE(:,50))/IDIM

       print*, 'KCBL  min max mean', minval(KCBL ),                          &
                                     maxval(KCBL ),                          &
                                     sum(KCBL)/IDIM

       print*, 'ZCBLx min max mean', minval(ZCBLx),                          &
                                     maxval(ZCBLx),                          &
                                     sum(ZCBLx)/IDIM
    ENDIF

    !------------------------------------------------------------------------
    ! Compute CNV_FRACTION
    ! NOTE: This block seems to be omitted, can probably comment out.
    !------------------------------------------------------------------------
    IF ( .FALSE. ) THEN 
       CNV_Q600_MIN = 0.00250_fp
       CNV_Q600_MAX = 0.00600_fp

       TMP_REAL = 0.0_fp
       
       ! Look for 600 mb grid cell
       DO J = 1, JJPAR
       DO I = 1, IIPAR
       DO L = 1, LLPAR
          QV600(I,J) = Q(I,J,L)
          EDGE_PRES(I,J) = PLE(I,J,L)
          IF ( EDGE_PRES(I,J) > 600.0_fp ) EXIT
       ENDDO
       ENDDO
       ENDDO
       QV600 = QV600 / 1000.0_fp

       IF ( CNV_Q600_MAX > CNV_Q600_MIN ) THEN
          DO J = 1, JJPAR
          DO I = 1, IIPAR
             TMP_REAL(I,J) = MAX(0.0_fp, MIN(1.0_fp, (QV600(I,J) - &
                             CNV_Q600_MIN) / (CNV_Q600_MAX - CNV_Q600_MIN)))
          ENDDO
          ENDDO
       ENDIF

       CONV_FRAC = TMP_REAL
   
!%%% BMY NOTE: Attach netCDF diagnostics here
!   AD66(:,:,1,9) = AD66(:,:,1,9) + TMP_REAL

       ! Reshape to 1-d
       CALL REARRANGE_GRID(TMP_REAL, CNV_FRAC)
    ENDIF

    CNV_FRAC = 1.0_fp

    !---------------------------------------------------------------------
    ! Compute TPERT and QPERT
    !---------------------------------------------------------------------
    CBL_QPERT        =  0.0_fp
    CBL_TPERT        = -2.0_fp
    CBL_TPERT_MXOCN  =  3.0_fp
    CBL_TPERT_MXLND  =  3.0_fp
   !CBL_TPERT_MXLND  =  0.0_fp
    QSSFC            =  0.0_fp ! for now
       
    ! Compute TPERT
    DO L = 1, LLPAR
       CALL REARRANGE_GRID(TEMP(:,:,L), TEMPERATURE(:,L))
    ENDDO
    CALL REARRANGE_GRID( State_Met%TSKIN, SURF_TEMP )

    TPERT  = ABS(CBL_TPERT) * ( SURF_TEMP - ( TEMPERATURE(:,LLPAR) + MAPL_GRAV*ZLO(:,LLPAR)/MAPL_CP )  ) 
    CALL REARRANGE_2D(TPERT,TMP_REAL)

!%%% BMY NOTE: Attach netCDF diagnostics here
!       AD66(:,:,6,9) = AD66(:,:,6,9) + TMP_REAL

    !if (CBL_TPERT < 0) then
    ! Make TPERT 0 in areas of deep convection
    !TPERT = TPERT*(1.0-CNV_FRAC)
    !endif
    QPERT  = CBL_QPERT * ( QSSFC - QHO(:,LLPAR) )          
    TPERT  = MAX( TPERT, 0.0_fp )
    QPERT  = MAX( QPERT, 0.0_fp )

    CALL REARRANGE_GRID(State_Met%FRLAND, FRLAND)

    WHERE ( FRLAND < 0.1_fp ) 
       TPERT = MIN( TPERT , CBL_TPERT_MXOCN ) ! ocean
    ELSEWHERE
       TPERT = MIN( TPERT , CBL_TPERT_MXLND ) ! land
    END WHERE

    IF ( prtDebug ) THEN
       print*, 'max tpert', maxval(TPERT)
       print*, 'min tpert', minval(TPERT)
    ENDIF

    !-----------------------------------------------------------------------
    ! Compute RASAL_2d
    !-----------------------------------------------------------------------
    IF ( RASPARAMS%RASAL2 > 0.0_fp ) THEN
       RASAL2_2d(:) = RASPARAMS%RASAL2
    ELSE
       ! include CNV dependence
       DO I=1, IDIM
          RASAL2_2d(I) = CNV_FRAC(I)                                         &
               * ABS( RASPARAMS%RASAL2 )                                     &
               + ( 1.0_fp - CNV_FRAC(I) )                                    &
               * RASPARAMS%RASAL1
       ENDDO
    ENDIF

    !------------------------------------------------------------------------
    ! Compute QSS And DQS
    !------------------------------------------------------------------------
    DO L = 1, LLPAR
    DO I = 1, IDIM
       CALL QSAT(TEMPERATURE(I,L), PLO(I,L), QSS(I,L), DQS(I,L), .TRUE.)
    ENDDO
    ENDDO

    IF ( prtDebug ) THEN
       print*, 'QSS', minval(QSS), maxval(QSS), sum(QSS)/(IDIM*LLPAR)
    ENDIF

    !------------------------------------------------------------------------
    ! Call RAS 
    !------------------------------------------------------------------------
    CALL RASE( IDIM          = IDIM,                                         &
               IRUN          = IRUN,                                         &
               K0            = LLPAR,                                        &
               ICMIN         = ICMIN,                                        &
               DT            = DT,                                           &
               CP0           = MAPL_CP,                                      &
               ALHLO         = MAPL_ALHL,                                    &
               GRAVO         = MAPL_GRAV,                                    &
               SEEDRAS       = SEEDRAS,                                      & 
               IRAS          = IRAS,                                         &
               JRAS          = JRAS,                                         &
               SIGE          = SIGE,                                         &
               !----------------------
               ! Inputs for CBL
               !----------------------
               KCBL          = KCBL,                                         &
               WGT0          = WGT0,                                         &
               WGT1          = WGT1,                                         &
               ZCBL          = ZCBLx,                                        &
               MXDIAM        = MXDIAMx,                                      &
               TPERT         = TPERT,                                        &
               QPERT         = QPERT                                         &
               !----------------------
               ! Inputs
               !----------------------
               TH0          = THO,                                           &
               QH0          = QHO,                                           & 
               QSS          = QSS,                                           & 
               DQS          = DQS,                                           &
               CNV_FRACTION = CNV_FRAC,                                      &
               RASAL2_2d    = RASAL2_2d,                                     &
               PLE          = CNV_PLE,                                       &
               PKE          = PKE,                                           &
               !----------------------
               ! Outputs
               !----------------------
               FLX          = CNV_FLX,                                       &
               FLXD         = CNV_FLXD,                                      &
               FLXC         = CNV_FLXC,                                      &
               RASPARAMS    = RASPARAMS,                                     &
               IRC          = IRC                                           )

    !------------------------------------------------------------------------
    ! Remap outputs to 2D grid
    !------------------------------------------------------------------------
    DO L = 1, LLPAR
       CALL REARRANGE_2D( CNV_FLXC(:,L), CMFMC(:,:,L)  ) 
       CALL REARRANGE_2D( CNV_FLXD(:,L), DTRAIN(:,:,L) )
    ENDDO
      
    CALL REARRANGE_2D( CNV_FLXC(:,LLPAR+1), CMFMC(:,:,LLPAR+1) )

    !------------------------------------------------------------------
    ! Save to State_Met
    !------------------------------------------------------------------
    State_Met%RAS_CMFMC  = CMFMC (:,:,LLPAR+1:1:-1)
    State_Met%RAS_DTRAIN = DTRAIN(:,:,LLPAR  :1:-1)

    !------------------------------------------------------------------
    ! Compute repartitioning
    !------------------------------------------------------------------
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       COLUMN_RAS = SUM( State_Met%RAS_CMFMC(I,J,:) )
       COLUMN_FP  = SUM( State_Met%CMFMC    (I,J,:) )
       !IF (CONV_FRAC(I,J) < 0.5) THEN
       IF ( COLUMN_RAS < 1e-2_fp                                  .or.        &
            ABS( COLUMN_RAS - COLUMN_FP) / COLUMN_RAS < 0.25_fp ) THEN
          PART(I,J) = 0.0_fp
       ELSE
          PART(I,J) = 1.0_fp - ( COLUMN_FP / COLUMN_RAS )
       ENDIF
    ENDDO
    ENDDO

    IF ( prtDebug ) THEN
       print*, 'max part', maxval(part)
       print*, 'min part', minval(part)
    ENDIF

    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       part_tmp = PART(I,J)

       IF ( part_tmp >= 0.0_fp ) THEN
          repartition_dqr = MIN( part_tmp * State_Met%DQRLSAN(I,J,L),        &
                                            State_Met%DQRLSAN(I,J,L)        )

          repartition_pfi = MIN( part_tmp * State_Met%PFILSAN(I,J,L),        &
                                            State_Met%PFILSAN(I,J,L)        )
       ELSE
          repartition_dqr = MAX( part_tmp * State_Met%DQRLSAN(I,J,L),        &
                                  -1.0_fp * State_Met%DQRCU(I,J,L)          )

          repartition_pfi = MAX( part_tmp * State_Met%PFILSAN(I,J,L),        &
                                 0.0 )!-1.0*State_Met%PFICU(I,J,L))
       ENDIF

       ! Repartition DQR
       State_Met%RAS_DQRCU(I,J,L)   = State_Met%DQRCU(I,J,L)   !+ repart'd_dqr 
       State_Met%RAS_DQRLSAN(I,J,L) = State_Met%DQRLSAN(I,J,L) !- repart'd_dqr

       ! Repartition PFI
       State_Met%RAS_PFICU(I,J,L)   = State_Met%PFICU(I,J,L)   !+ repart'd_pfi
       State_Met%RAS_PFILSAN(I,J,L) = State_Met%PFILSAN(I,J,L) !- repart'd_pfi
       
       IF ( part_tmp >= 0.0_fp ) THEN
          repartition = MIN( part_tmp * State_Met%PFLLSAN(I,J,L),            &
                                       State_Met%PFLLSAN(I,J,L)             )
       ELSE
          repartition = MAX( part_tmp * State_Met%PFLLSAN(I,J,L),            &
                             0.0 ) !-1.0*State_Met%PFLCU(I,J,L))
       ENDIF

       ! Repartition PFL
       State_Met%RAS_PFLCU(I,J,L)   = State_Met%PFLCU(I,J,L)   !+ repartition
       State_Met%RAS_PFLLSAN(I,J,L) = State_Met%PFLLSAN(I,J,L) !- repartition

       IF ( part_tmp >= 0.0_fp ) THEN
          repartition = MIN( part_tmp * State_Met%REEVAPLS(I,J,L),           &
                             0.5_fp   * State_Met%REEVAPLS(I,J,L)           )
       ELSE
          repartition = MAX( part_tmp * 0.5_fp*State_Met%REEVAPLS(I,J,L),    &
                             0.0) ! -1.0*State_Met%REEVAPCN(I,J,L))
       ENDIF

       ! Repartition REEVAP
       !IF (State_Met%REEVAPCN(I,J,L) <0 .and. repartition <0) THEN
       !    repartition = 0.0
       !ENDIF
       State_Met%RAS_REEVAPCN(I,J,L) = State_Met%REEVAPCN(I,J,L) !+ repartition
       State_Met%RAS_REEVAPLS(I,J,L) = State_Met%REEVAPLS(I,J,L) !- repartition

    ENDDO
    ENDDO
    ENDDO

    !### Debug
    IF ( prtDebug ) THEN
       PRINT*, 'max DQR',    MAXVAL( State_Met%PFICU        ),               &
                             MAXVAL( State_Met%RAS_PFICU    )
       PRINT*, 'min DQR',    MINVAL( State_Met%PFICU        ),               &
                             MINVAL( State_Met%RAS_PFICU    )
       print*, 'max REEVAP', MAXVAL( State_Met%REEVAPCN     ),               &
                             MAXVAL( State_Met%RAS_REEVAPCN )
       print*, 'min REEVAP', MINVAL( State_Met%REEVAPCN     ),               &
                             MINVAL( State_Met%RAS_REEVAPCN )
    ENDIF

    !------------------------------------------------------------------
    !%%% BMY NOTE: Update netCDF diagnostics here
    ! Write to diagnostic
    !------------------------------------------------------------------
    !AD66(:,:,:,7) = AD66(:,:,:,7) + State_Met%RAS_CMFMC
    !AD66(:,:,:,8) = AD66(:,:,:,8) + State_Met%RAS_DTRAIN
    !AD66(:,:,2,9) = AD66(:,:,2,9) + part
    !AD66(:,:,3,9) = AD66(:,:,3,9) + State_Met%PBLH_MAX
    
    ! Free pointers
    RasParams => NULL()
    TEMP      => NULL() 
    PMID_DRY  => NULL()
    PLE       => NULL()
    Q         => NULL() 
    DELP      => NULL() 

  END SUBROUTINE DO_RAS
!EOC
!------------------------------------------------------------------------
!               GEOS-Chem Global Chemical Transport Model
!-------------------------------------------------------------------------
!BOP
!
      SUBROUTINE RASE(IDIM, IRUN, K0, ICMIN, DT ,            &
         CPO,ALHLO,GRAVO,                                 &
         SEEDRAS,IRAS,JRAS,SIGE,                          &
         KCBL,WGT0,WGT1,ZCBL,MXDIAM,TPERT,QPERT,          &
         THO, QHO,                                        & 
         QSS, DQS, CNV_FRACTION, RASAL2_2d,               &
         PLE, PKE, FLX, FLXD, FLXC,                       &
         RASPARAMS,                                       &
         IRC                                              )


      !*********************************************************************
      !*********************************************************************
      !******************** Relaxed Arakawa-Schubert ***********************
      !************************ Parameterization ***************************
      !********************** SCALAR RAS-1 VERSION  ************************
      !************************* 31 DECEMBER 1999 **************************
      !*********************************************************************
      !************************** Developed By *****************************
      !*********************************************************************
      !************************ Shrinivas Moorthi **************************
      !******************************* and *********************************
      !************************** Max J. Suarez ****************************
      !*********************************************************************
      !******************** Laboratory for Atmospheres *********************
      !****************** NASA/GSFC, Greenbelt, MD 20771 *******************
      !*********************************************************************
      !*********************************************************************

      !  Input:
      !  ------
      ! 
      !     K0      : Number of vertical levels (increasing downwards)
      !
      !     DT      : Time step in seconds
      !
      !     RASAL   : Array of dimension K-1 containing relaxation parameters
      !               for cloud-types detraining at those levels
      !
      !     CPO     : Specific heat at constant pressure (J/kg/K)
      !
      !     ALHLO   : Latent Heat of condensation (J/kg)
      !
      !     GRAVO   : Acceleration due to gravity (m/s^2)
      !
      !     PLE     : 2D array of dimension (IDIM,K0+1) containing pressure
      !               in hPa at the interfaces of K-layers from top of the 
      !               atmosphere to the bottom  (mb)
      !
      !     PKE     : 2D array of dimension (IDIM,K0+1) containing (PRS/P00) **
      !               RKAP.  i.e. Exner function at layer edges.
      !
      !     PKL     : 2D array of dimension (IDIM,K0) ) containing the
      !               Exner function at the layers.
      !
      !     QSS     : 2D array of dimension (IDIM,K0  ) containing the
      !               saturation specific humidity at the layers. (kg/kg)
      !
      !     DQS     : 2D array of dimension (IDIM,K0  ) containing
      !               d(qss)/dt at the layers.  (1/K)
      !   
      !     CNV_FRACTION    : 1D array of dimension (IDIM) containing
      !               fraction of grid cell considered to be convective
      !   
      !  Update:
      !  -------
      !
      !     THO     : 2D array of dimension (IDIM,K0) containing potential
      !               temperature (K)
      !
      !     QHO     : 2D array of dimension (IDIM,K0) containing specific
      !               humidity (kg/kg)
      !
      !  Output:
      !  -------
      !!
      !     FLX     : 2D array of dimension (IDIM,K0) containing the
      !               cloud-base mass flux for each cloud type ordered by
      !               detrainment level.   (kg/m^2/s) 
      !
      !     FLXD    : 2D array of dimension (IDIM,K0) containing the
      !               detrained  mass flux for each cloud type ordered by
      !               detrainment level.   (kg/m^2/s) 
      !
      !     FLXC    : 2D array of dimension (IDIM,K0+1) containing the
      !               total cloud mass flux for all cloud types through
      !               the top of each level. (e.g., FLXC(K)=SUM(FLX(ICMIN:K))
      !               and  FLXD(L) = FLXC(L+1)-FLSD(L) )
      !                (kg/m^2/s) 
      !
      !
      !************************************************************************

    USE PRECISION_MOD

      !  ARGUMENTS

      INTEGER,                     INTENT(IN   ) ::  IDIM, IRUN, K0, ICMIN
      REAL(fp), DIMENSION (IDIM,K0  ), INTENT(INOUT) ::  THO, QHO
      REAL(fp), DIMENSION (IDIM,K0+1), INTENT(IN   ) ::  PLE, PKE
      REAL(fp), DIMENSION (IDIM,K0  ), INTENT(IN   ) ::  QSS, DQS
      REAL(fp), DIMENSION (IDIM     ), INTENT(IN   ) ::  CNV_FRACTION, RASAL2_2d
      REAL(fp), DIMENSION (     K0+1), INTENT(IN   ) ::  SIGE
      REAL(fp), DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  FLX, FLXD
      REAL(fp), DIMENSION (IDIM,K0+1), INTENT(  OUT) ::  FLXC
      REAL(fP),                        INTENT(IN   ) ::  DT,  CPO, ALHLO, GRAVO
      INTEGER, DIMENSION (IDIM,2), INTENT(IN   ) ::  SEEDRAS
      INTEGER, DIMENSION (IDIM  ), INTENT(IN   ) ::  IRAS,JRAS,KCBL
      REAL(fp), DIMENSION (IDIM     ), INTENT(IN   ) ::  ZCBL,TPERT,QPERT
      REAL(fp), DIMENSION (IDIM     ), INTENT(  OUT) ::  MXDIAM
      REAL(fp), DIMENSION (IDIM,K0  ), INTENT(IN   ) ::  WGT0,WGT1

      !     REAL, DIMENSION(:),          INTENT(IN   ) ::  RASPARAMS
      type (RASPARAM_TYPE),        INTENT(IN   ) ::  RASPARAMS

      INTEGER, DIMENSION(IDIM,K0), INTENT(  OUT) ::  IRC

      !  LOCALS


      REAL(fp),  DIMENSION(K0) :: POI_SV, QOI_SV
      REAL(fp),  DIMENSION(K0) :: POI, QOI, DQQ, BET, GAM
      REAL(fp),  DIMENSION(K0) :: POI_c, QOI_c
      REAL(fp),  DIMENSION(K0) :: PRH,  PRI,  GHT, DPT, DPB, PKI
      REAL(fp),  DIMENSION(K0) :: TCU, QCU, UCU, VCU,  CLN, RNS, POL,DM
      REAL(fp),  DIMENSION(K0) :: QST, SSL,  RMF, RNN, RN1, RMFC, RMFP
      REAL(fp),  DIMENSION(K0) :: GMS, ETA, GMH, EHT,  HCC, RMFD
      REAL(fp),  DIMENSION(K0) :: HOL, HST, QOL, ZOL, HCLD
      REAL(fp),  DIMENSION(K0) :: LAMBDSV, BKE , CVW, UPDFRC
      REAL(fp),  DIMENSION(K0) :: RASAL, MTKWI, UPDFRP,BK2,BK3



      REAL(fp),  DIMENSION(K0+1) :: PRJ, PRS, QHT, SHT ,ZET, XYD, XYD0

      INTEGER,  DIMENSION(K0-1) :: RC

      INTEGER :: K,MY_PE

      REAL(fp), DIMENSION(IDIM,K0) :: LAMBDSV2


      REAL(fp) TX2, TX3, UHT, VHT, AKM, ACR, ALM, TTH, QQH, SHTRG, DQX
      REAL(fp) WFN, TEM, TRG, TRGEXP, EVP, WLQ, QCC,MTKW_MAX !, BKE
      REAL(fp) SHTRG_FAC, SIGE_MINHOL, WFNOG

      INTEGER I, IC, L, ICL , ITR , ICL_C, N_DTL
      INTEGER NDTLEXPON

      INTEGER,  ALLOCATABLE, DIMENSION(:) :: ICL_V

      !  RASE GLOBAL CONSTANTS


      REAL(fp) GRAV, CP, ALHL, CPBG, ALHI, CPI, GRAVI, DDT, LBCP, OBG, AFC 

!!!!!!!!!
      REAL(fp) PBLFRAC,AUTORAMPB,CO_ZDEP
      REAL(fp) RASAL1, RASAL2, RASAL2i, CO_T, RASNCL,FRICLAMBDA,SDQVT1,SDQV2 
      REAL(fp) LAMBDA_FAC,STRAPPING,ACRITFAC,HMINTRIGGER,LLDISAGGXP
      REAL(fp) LAMBMX_FAC, DIAMMN_MIN,RDTLEXPON, CLI_CRIT,SDQV3, MAXDALLOWED_D, MAXDALLOWED_S
      REAL(fp) RHMN, RHMX
      INTEGER RASAL_EXP

      real(fp) cld_radius, areal_frac, spect_mflx, cvw_cbase
!!!!!!!!!

      REAL(fp), PARAMETER :: ONEPKAP = 1.+ 2./7., DAYLEN = 86400.0
      !      REAL, PARAMETER :: PBLFRAC = 0.5
      REAL(fp), PARAMETER :: RHMAX   = 0.9999

      !  LAMBDA LIMITS
      REAL(fp)            :: LAMBDA_MIN
      REAL(fp)            :: LAMBDA_MAX

      !  TRIGGER PARAMETERS
      real(fp), parameter :: RHO_W  =  1.0e3  ! Density of liquid water in kg/m^3


      character(len=255)          :: CBL_STYLE

      real*8, dimension(k0) :: tcu8, qcu8, pcu, flx8
      real*8    :: cup  !, dpd, tla
      logical :: revap, wrkfun, calkpb, crtfun, lprnt, dndrft

      real*8, dimension(k0) :: toi8, qoi8, prsm8, phil8, qli8, qii8, trcfac
      real*8, dimension(k0) :: ALFIND,ALFINT,ALFINQ,RHC_LS
      real*8, dimension(k0+1) :: prs8, phih8
      real*8 :: FRACBL, dt8, rasalf
      integer :: KPBL



      real*8 :: MAX_NEG_BOUY = 1.0 ! no inhibition for =1.0
      !!real*8 :: ALFINT = 0.5
      !!real*8 :: ALFINQ = 0.5
      real*8 :: RHFACL = 0.0 ! not used
      real*8 :: RHFACS = 0.0 ! no inhibition
      !!real*8 :: ALFIND = 1.0
      !!real*8 :: RHC_LS = 0.80
      real*8 :: dsfc   = 0.001
      real*8 :: cd     = 1.e-3
      real*8 :: wfnc   = 0.0
      real*8 :: tla    = -1.0
      real*8 :: dpd    = 300.


      ! *********************************************************************


      IF(IRUN <= 0) RETURN

      IRC = -2

      RASAL1        = RASPARAMS%RASAL1
      RASAL2        = RASPARAMS%RASAL2
      RASNCL        = RASPARAMS%RASNCL
      LAMBDA_FAC    = RASPARAMS%LAMBDA_FAC
      LAMBMX_FAC    = RASPARAMS%LAMBMX_FAC
      DIAMMN_MIN    = RASPARAMS%MIN_DIAMETER
      RDTLEXPON     = RASPARAMS%RDTLEXPON
      ACRITFAC      = RASPARAMS%ACRITFAC
      PBLFRAC       = RASPARAMS%PBLFRAC
      AUTORAMPB     = RASPARAMS%RASAUTORAMPB
      CO_ZDEP       = RASPARAMS%AUTOC_CN_ZDEP
      MAXDALLOWED_S = RASPARAMS%MAXDALLOWED_S
      MAXDALLOWED_D = RASPARAMS%MAXDALLOWED_D
      RASAL_EXP     = RASPARAMS%RASAL_EXP
      RHMN          = RASPARAMS%RAS_RHMIN
      RHMX          = RASPARAMS%RAS_RHFULL

      GRAV  = GRAVO
      ALHL  = ALHLO
      CP    = CPO
      CPI   = 1.0/CP      
      ALHI  = 1.0/ALHL
      GRAVI = 1.0/GRAV
      CPBG  = CP*GRAVI
      DDT   = DAYLEN/DT
      AFC   = -1.04E-4*SQRT(DT*113.84)
      LBCP  = ALHL*CPI
      OBG   = 100.*GRAVI

      DO I=1,IRUN

         K = KCBL(I)

         rc(icmin) = 0

         CALL FINDDTLS

         IF ( K > 0 ) THEN 
            CALL STRAP( FINAL=0 )
            call HTEST

            DO ICL_C = 1,N_DTL
               ICL   = ICL_V( ICL_C )

               IF( ICL > ICMIN ) THEN
                  CALL CLOUDE(ICL)
               ENDIF

            ENDDO

            IF ( SUM( RMF(ICMIN:K) ) > 0.0 ) THEN

               CALL STRAP( FINAL=1 )

            ELSE

               CALL STRAP( FINAL=2 )

            ENDIF

         ELSE 

            CALL STRAP( FINAL=2 )

         ENDIF

      ENDDO

      IF(ALLOCATED(ICL_V)) DEALLOCATE(ICL_V)

      RETURN

   CONTAINS

      !*********************************************************************

      SUBROUTINE CLOUDE(IC)

         INTEGER, INTENT(IN ) :: IC
         REAL :: DEEP_FACT,WSCALE

         REAL :: CLI , TE_A, C00_X, CLI_CRIT_X, PETE, TOKI, GMHx, HSTx  !, dQx
         REAL :: DT_LYR, RATE, CLOSS, F2, F3, F4, F5
         INTEGER :: K700

         !=============================AER_CLOUD local variables ====================
         REAL :: WBASE, NDROP, NICE, FP_D, FF_A, FP_I, FICE, &
               NDROP_AMB, NSOOT_AMB, NSOOT, NIN, INSOOT, dCVW2, QICE, &
               dQICE, dQIG, FPICE, dNICE, dNDROP, DSOOT_AMB, DSOOT, QLIQ, dQLIQ, FPRECIP, AUX, QT, &
               MAXNICE, MAXNDROP, MINNICE, MINNDROP, NDROP_ACT, RIMM, FNDRIM, TminusTa, Tparcel , &
	      alph_e, beta_e, RH_AMB, ECRIT

         INTEGER :: INX, naux

         ALM   = 0.
         TRG   = AMIN1(1.,(QOI(K)/QST(K)-RHMN)/(RHMX-RHMN))

         F4  = MIN(   1.0,  MAX( 0.0 , (AUTORAMPB-SIGE(IC))/0.2 )  )  ! F4 should ramp from 0 at SIG=AUTORAMPB
         ! to 1 at SIG=AUTORAMPB-0.2

         if ( SIGE(IC) >= 0.5 ) then
            F5 = 1.0
         else
            F5 = 1.0 - 2.*CO_ZDEP *( 0.5 - SIGE(IC) )
            F5 = MAX( F5 , 0.0 )
         endif


         IF(TRG <= 1.0E-5) THEN    ! TRIGGER  =========>>
            RC(IC) = 7
            RETURN
         ENDIF

         !  RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL

         POI_c = POI
         QOI_c = QOI
         POI_c(K) =  POI_c(K) + TPERT(I)
         QOI_c(K) =  QOI_c(K) + QPERT(I)

         ZET(K+1) = 0.
         SHT(K+1) = CP*POI_c(K)*PRJ(K+1)
         DO L=K,IC,-1
            QOL(L)  = AMIN1(QST(L)*RHMAX,QOI_c(L))
            QOL(L)  = AMAX1( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
            SSL(L)  = CP*PRJ(L+1)*POI_c(L) + GRAV*ZET(L+1)
            HOL(L)  = SSL(L) + QOL(L)*ALHL
            HST(L)  = SSL(L) + QST(L)*ALHL
            TEM     = POI_c(L)*(PRJ(L+1)-PRJ(L))*CPBG
            ZET(L)  = ZET(L+1) + TEM
            ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI_c(L)*CPBG
         ENDDO

         DO L=IC+1,K
            TEM  = (PRJ(L)-PRH(L-1))/(PRH(L)-PRH(L-1))
            SHT(L)  = SSL(L-1) + TEM*(SSL(L)-SSL(L-1)) 
            QHT(L)  = .5*(QOL(L)+QOL(L-1))
         ENDDO

         !  CALCULATE LAMBDA, ETA, AND WORKFUNCTION

         LAMBDA_MIN = .2/MXDIAM(I)
         LAMBDA_MAX = .2/  200. 


         IF (HOL(K) <= HST(IC)) THEN   ! CANNOT REACH IC LEVEL  ======>>
            RC(IC) = 1
            RETURN
         ENDIF

         !  LAMBDA CALCULATION: MS-A18

         TEM  =       (HST(IC)-HOL(IC))*(ZOL(IC)-ZET(IC+1)) 
         DO L=IC+1,K-1
            TEM = TEM + (HST(IC)-HOL(L ))*(ZET(L )-ZET(L +1))
         ENDDO

         IF(TEM <= 0.0) THEN         ! NO VALID LAMBDA  ============>>
            RC(IC) = 2
            RETURN
         ENDIF

         ALM     = (HOL(K)-HST(IC)) / TEM

         IF(ALM > LAMBDA_MAX) THEN
            RC(IC) = 3
            RETURN
         ENDIF

         TOKI=1.0

         IF(ALM < LAMBDA_MIN) THEN
       TOKI = ( ALM/LAMBDA_MIN )**2  !we can probably replace this by a actual distribution based on grid cell size

            !RC(IC) = 6
            !RETURN
         ENDIF

         !  ETA CALCULATION: MS-A2

         DO L=IC+1,K
            ETA(L) = 1.0 + ALM * (ZET(L )-ZET(K))
         ENDDO
         ETA(IC) = 1.0 + ALM * (ZOL(IC)-ZET(K))

         !  WORKFUNCTION CALCULATION:  MS-A22

         WFN     = 0.0
         HCC(K)  = HOL(K)
         DO L=K-1,IC+1,-1
            HCC(L) = HCC(L+1) + (ETA(L) - ETA(L+1))*HOL(L)
            TEM    = HCC(L+1)*DPB(L) + HCC(L)*DPT(L)
            EHT(L) = ETA(L+1)*DPB(L) + ETA(L)*DPT(L)
            WFN    = WFN + (TEM - EHT(L)*HST(L))*GAM(L)
         ENDDO
         HCC(IC) = HST(IC)*ETA(IC)
         WFN     = WFN + (HCC(IC+1)-HST(IC)*ETA(IC+1))*GAM(IC)*DPB(IC)

         !   ALPHA CALCULATION 
         RASAL2i = RASAL2_2d(I)
         IF ( ZET(IC) <  2000. ) RASAL(IC) = RASAL1
         IF ( ZET(IC)  >= 2000. ) THEN 
   !WMP     RASAL(IC) = RASAL1 + (RASAL2i-RASAL1)*(ZET(IC) - 2000.)/8000.
            RASAL(IC) = RASAL1 + (RASAL2i-RASAL1)*MIN(1.0,(ZET(IC) - 2000.)/8000.)**RASAL_EXP
         ENDIF
   !WMP  RASAL(IC) = MIN( RASAL(IC) , 1.0e5 )

         RASAL(IC) = DT / RASAL(IC)

         !  TEST FOR CRITICAL WORK FUNCTION

         CALL ACRITN(POL(IC), PRS(K), ACR)

         IF(WFN <= ACR) THEN   ! SUB-CRITICAL WORK FUNCTION ======>>
            RC(IC) = 4
            RETURN
         ENDIF


         !print *, '========================================='
         !DO L=K-1,IC,-1

         !     CALCULATE GAMMAS AND KERNEL

         GMS(K) =          (SHT(K)-SSL(K))*PRI(K)          ! MS-A30 (W/O GRAV)
         GMH(K) = GMS(K) + (QHT(K)-QOL(K))*PRI(K)*ALHL     ! MS-A31 (W/O GRAV)
         AKM    = GMH(K)*GAM(K-1)*DPB(K-1)                 ! MS-A37 (W/O GRAV)

         TX2     = GMH(K)
         DO L=K-1,IC+1,-1
            GMS(L) = ( ETA(L  )*(SHT(L)-SSL(L  ))   & 
                  + ETA(L+1)*(SSL(L)-SHT(L+1)) )     *PRI(L)
            GMH(L) = GMS(L)                         &
                  + ( ETA(L  )*(QHT(L)-QOL(L  ))   &
                  + ETA(L+1)*(QOL(L)-QHT(L+1)) )*ALHL*PRI(L)
            TX2 = TX2 + (ETA(L) - ETA(L+1)) * GMH(L)
            AKM = AKM - GMS(L)*EHT(L)*PKI(L) + TX2*GHT(L)

         ENDDO

         GMS(IC) = ETA(IC+1)*(SSL(IC)-SHT(IC+1))*PRI(IC)
         AKM     = AKM - GMS(IC)*ETA(IC+1)*DPB(IC)*PKI(IC)

         GMH(IC) =   GMS(IC) + ( ETA(IC+1)*(QOL(IC)-QHT(IC+1))*ALHL &
               + ETA(IC  )*(HST(IC)-HOL(IC  ))     )*PRI(IC)

         !    CLOUD BASE MASS FLUX

         IF (AKM >= 0.0 .OR. WLQ < 0.0)  THEN  !  =========>
            RC(IC) = 5
            RETURN
         ENDIF

         WFN = - (WFN-ACR)/AKM ! MS-A39 MASS-FLUX IN Pa/step
         WFN = MIN( ( RASAL(IC)*TRG*TOKI )*WFN  ,   (PRS(K+1)-PRS(K) )*(100.*PBLFRAC))



         !    CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT

         WFNOG    = WFN*GRAVI
         TEM      = WFN*GRAVI
         RMF (IC) = RMF (IC) +     TEM           ! (kg/m^2/step)
         RMFD(IC) = RMFD(IC) +     TEM * ETA(IC) ! (kg/m^2/step)

         DO L=IC+1,K
            RMFP(L) = TEM * ETA(L)                 ! (kg/m^2/step)
            RMFC(L) = RMFC(L)  +  RMFP(L)          ! (kg/m^2/step)

         ENDDO



         !    THETA AND Q CHANGE DUE TO CLOUD TYPE IC

         DO L=IC,K
            GMH(L) = GMH(L) * WFN
            GMS(L) = GMS(L) * WFN
            QOI(L) = QOI(L) + (GMH(L) - GMS(L)) * ALHI
            POI(L) = POI(L) + GMS(L)*PKI(L)*CPI
            QST(L) = QST(L) + GMS(L)*BET(L)*CPI
         ENDDO

         RC(IC) = 0

         RETURN

      END SUBROUTINE CLOUDE

      SUBROUTINE ACRITN( PL, PLB, ACR)

         REAL(fp), INTENT(IN ) :: PL, PLB
         REAL(fp), INTENT(OUT) :: ACR

         INTEGER IWK

         !!REAL, PARAMETER :: FACM=0.5

         REAL(fp), PARAMETER :: &
               PH(15)=(/150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, &
               550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0/)

         REAL(fp), PARAMETER :: & 
               A(15)=(/ 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677, &
               0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664, &
               0.0553, 0.0445, 0.0633     /)   !!*FACM

         IWK = PL * 0.02 - 0.999999999

         IF (IWK .GT. 1 .AND. IWK .LE. 15) THEN
            ACR = A(IWK-1) + (PL-PH(IWK-1))*.02*(A(IWK)-A(IWK-1))
         ELSEIF(IWK > 15) THEN
            ACR = A(15)
         ELSE
            ACR = A(1)
         ENDIF

         ACR = ACRITFAC  * ACR * (PLB - PL)

         RETURN

      END SUBROUTINE ACRITN



      SUBROUTINE HTEST
         REAL,  DIMENSION(K0) :: HOL1
         integer  :: lminhol
         real     :: minhol


         hol=0.            ! HOL initialized here in order not to confuse Valgrind debugger
         lminhol  = K+1
         MINHOL   = -999999.
         ZET(K+1) = 0
         SHT(K+1) = CP*POI(K)*PRJ(K+1)
         DO L=K,ICMIN,-1
            QOL(L)  = AMIN1(QST(L)*RHMAX,QOI(L))
            QOL(L)  = AMAX1( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
            SSL(L)  = CP*PRJ(L+1)*POI(L) + GRAV*ZET(L+1)
            HOL(L)  = SSL(L) + QOL(L)*ALHL
            HST(L)  = SSL(L) + QST(L)*ALHL
            TEM     = POI(L)*(PRJ(L+1)-PRJ(L))*CPBG
            ZET(L)  = ZET(L+1) + TEM
            ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI(L)*CPBG
         ENDDO

         HOL1=HOL
         DO L=K-1,ICMIN+1,-1
            HOL1(L)  = (0.25*HOL(L+1)+0.50*HOL(L)+0.25*HOL(L-1))
            if ( ( MINHOL>=HOL1(L) ) .OR. (MINHOL<0.) ) THEN
               MINHOL   =  HOL1(L)
               LMINHOL  =  L
            end if
         ENDDO

         SIGE_MINHOL = SIGE(LMINHOL)

      end subroutine HTEST


      subroutine FINDDTLS
         real :: SIGDT0,sigmax,sigmin
         integer :: LL
#ifndef __GFORTRAN__
         integer :: THE_SEED(2)
#else
         integer :: THE_SEED(12)
#endif

         THE_SEED(1)=SEEDRAS(I,1)*IRAS(I) + SEEDRAS(I,2)*JRAS(I)
         THE_SEED(2)=SEEDRAS(I,1)*JRAS(I) + SEEDRAS(I,2)*IRAS(I)
         THE_SEED(1)=THE_SEED(1)*SEEDRAS(I,1)/( SEEDRAS(I,2) + 10)
         THE_SEED(2)=THE_SEED(2)*SEEDRAS(I,1)/( SEEDRAS(I,2) + 10)
         if(THE_SEED(1) == 0) THE_SEED(1) =  5
         if(THE_SEED(2) == 0) THE_SEED(2) = -5

#ifdef __GFORTRAN__
         THE_SEED(3:12) = 0
#endif

         call random_seed(PUT=THE_SEED)

         SIGMAX=SIGE(K)
         SIGMIN=SIGE(ICMIN)

         if ( RASNCL < 0.0 ) then 
            !! NO SHALLOW CONV   N_DTL = 56 - ICMIN 
            N_DTL = K - ICMIN 
         else 
            N_DTL = min( int( RASNCL ) , K-ICMIN )
         endif

         if(allocated(ICL_V)) deallocate(ICL_V)
         allocate(ICL_V(N_DTL))

         if ( ( RASNCL < 0.0 ) .and. ( RASNCL >=-100.) ) then 
            do L=1,N_DTL
               ICL_V(L) = ICMIN + L - 1
            enddo
         else if ( RASNCL < -100.0 ) then 
            do L=1,N_DTL
               ICL_V(L) = K - L
               !! NO SHALLOW CONV           ICL_V(L) = 56 - L
            enddo
         else
            do L=1,N_DTL
               call random_number ( SIGDT0 )
               SIGDT0 = 1.00 - ( SIGDT0**RDTLEXPON )
               SIGDT0 = SIGMIN+SIGDT0*(SIGMAX-SIGMIN)

               do LL=ICMIN,K
                  if ( (SIGE(LL+1)>=SIGDT0) .and. (SIGE(LL)<SIGDT0 ) ) ICL_V(L)=LL
               enddo
            end do
         endif

      end subroutine FINDDTLS

      SUBROUTINE STRAP(FINAL)

         INTEGER :: FINAL
         REAL , DIMENSION(K0)  :: WGHT, MASSF

         REAL :: WGHT0, PRCBL

         integer, parameter :: nrands=1
         real ::    rndu(nrands)
         integer :: seedcbl(nrands)

         ! !DESCRIPTION: 
         !   {\tt STRAP} is called: FINAL=0, to compute cloud base layer CBL properties
         !   given a value K for the index of the upper {\em EDGE} of the CBL; FINAL=1
         !   to redistribute convective tendencies within CBL

         integer :: KK

         !  LOCAL VARIABLES FOR USE IN CLOUDE

         !!IF (.NOT. PRESENT(FINAL)) THEN
         IF (FINAL==0) THEN

            do kk=icmin,k+1
               PRJ(kk) = PKE(I,kk)
            enddo

            poi=0.        ! These initialized here in order not to confuse Valgrind debugger
            qoi=0.        ! Do not believe it actually makes any difference.

            PRS(ICMIN:K0+1) = PLE(I,ICMIN:K0+1)
            POI(ICMIN:K)   = THO(I,ICMIN:K)
            QOI(ICMIN:K)   = QHO(I,ICMIN:K)

            QST(ICMIN:K) = QSS(I,ICMIN:K)
            DQQ(ICMIN:K) = DQS(I,ICMIN:K)

!!! Mass fraction of each layer below cloud base
!!! contributed to aggregate cloudbase layer (CBL) 
            MASSF(:) = WGT0(I,:)

!!! RESET PRESSURE at bottom edge of CBL 
            PRCBL = PRS(K)
            do l= K,K0
               PRCBL = PRCBL + MASSF(l)*( PRS(l+1)-PRS(l) )
            end do
            PRS(K+1) = PRCBL
            PRJ(K+1) = (PRS(K+1)/1000.)**(MAPL_RGAS/CP)

            DO L=K,ICMIN,-1
               POL(L)  = 0.5*(PRS(L)+PRS(L+1))
               PRH(L)  = (PRS(L+1)*PRJ(L+1)-PRS(L)*PRJ(L)) &
                     / (ONEPKAP*(PRS(L+1)-PRS(L)))
               PKI(L)  = 1.0 / PRH(L)
               DPT(L)  = PRH(L  ) - PRJ(L)
               DPB(L)  = PRJ(L+1) - PRH(L)
               PRI(L)  = .01 / (PRS(L+1)-PRS(L))
            ENDDO

!!!!! RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
            if( K<=K0) then
               POI(K) = 0.
               QOI(K) = 0.

               !! SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
               WGHT = 0.
               DO L=K,K0
                  WGHT(L)   = MASSF(L) *                &
                        ( PLE(I,L+1) - PLE(I,L) ) &
                        /( PRS(K+1)   - PRS(K)  )
               END DO

               DO L=K,K0
                  POI(K) = POI(K) + WGHT(L)*THO(I,L)
                  QOI(K) = QOI(K) + WGHT(L)*QHO(I,L)
               ENDDO

               CALL QSAT( POI(K)*PRH(K) , POL(K), QST(K), DQQ(K), .FALSE. )

            endif

            rndu(:) = max( seedras(I,1)/1000000., 1e-6 )

            MXDIAM(I) = CNV_FRACTION(I)*ABS(MAXDALLOWED_D) + (1-CNV_FRACTION(I))*ABS(MAXDALLOWED_S)
            if (MAXDALLOWED_D > 0) then
               ! Make MXDIAM stochastic
               MXDIAM(I) = MXDIAM(I)*( rndu(1)**(-1./2.) )
            endif

            DO L=K,ICMIN,-1
               BET(L)  = DQQ(L)*PKI(L)  !*
               GAM(L)  = PKI(L)/(1.0+LBCP*DQQ(L)) !*
               IF(L<K) THEN
                  GHT(L+1) = GAM(L)*DPB(L) + GAM(L+1)*DPT(L+1)
               ENDIF
            ENDDO

            TCU(ICMIN:K) = -POI(ICMIN:K)*PRH(ICMIN:K)
            QCU(ICMIN:K) = -QOI(ICMIN:K)

            RMF  = 0.
            RMFD = 0.
            RMFC = 0.
            RMFP = 0.

            POI_SV = POI
            QOI_SV = QOI

         END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (FINAL==1) THEN 

            THO(I,ICMIN:K-1) = POI(ICMIN:K-1)
            QHO(I,ICMIN:K-1) = QOI(ICMIN:K-1)

            !! De-strap tendencies from RAS
            !! specify weighting "SHAPE"
            WGHT   = WGT1(I,:)

            !! Scale properly by layer masses
            wght0 = 0.
            DO L=K,K0 
               wght0 = wght0 + WGHT(L)* ( PLE(I,L+1) - PLE(I,L) )
            END DO

            wght0 = ( PRS(K+1)   - PRS(K)  )/wght0

            WGHT  = wght0 * WGHT


            DO L=K,K0 
               THO(I,L) =  THO(I,L) + WGHT(L)*(POI(K) - POI_SV(K))
               QHO(I,L) =  QHO(I,L) + WGHT(L)*(QOI(K) - QOI_SV(K))
            END DO

            FLX (I,ICMIN:K) = RMF (ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD BASE)
            FLXD(I,ICMIN:K) = RMFD(ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD TOP)
            FLXC(I,ICMIN:K) = RMFC(ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD TOP)

            FLX (I,1:ICMIN-1) = 0.
            FLXD(I,1:ICMIN-1) = 0.
            FLXC(I,1:ICMIN-1) = 0.


            IF ( K < K0 ) THEN 
               FLX (I,K:K0) = 0.
               FLXD(I,K:K0) = 0.
               FLXC(I,K:K0) = 0.
            END IF

            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)

            IF(ALLOCATED(ICL_V)) DEALLOCATE( ICL_V )

         ENDIF

         IF (FINAL==2) THEN 


            FLX (I,:) = 0.
            FLXD(I,:) = 0.
            FLXC(I,:) = 0.

            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)

         ENDIF

         RETURN

      END SUBROUTINE STRAP


      END SUBROUTINE RASE
!EOC
!--------------------------------------------------------------------------
!               GEOS-Chem Global Chemical Transport Model
!--------------------------------------------------------------------------
      SUBROUTINE qsat (TT,P,Q,DQDT,LDQDT)
!***********************************************************************
!
!  PURPOSE:
!  ========
!    Compute Saturation Specific Humidity
!
!  INPUT:
!  ======
!    TT ......... Temperature (Kelvin)
!    P .......... Pressure (mb)
!    LDQDT ...... Logical Flag to compute QSAT Derivative
!
!  OUTPUT:
!  =======
!    Q .......... Saturation Specific Humidity
!    DQDT ....... Saturation Specific Humidity Derivative wrt Temperature
!
!
!***********************************************************************
!*                  GODDARD LABORATORY FOR ATMOSPHERES                 *
!***********************************************************************

      USE PRECISION_MOD

      IMPLICIT NONE
      REAL(fp)   TT, P, Q, DQDT
      LOGICAL LDQDT
      REAL(fp)   AIRMW, H2OMW
      
      PARAMETER ( AIRMW  = 28.97      )                                         
      PARAMETER ( H2OMW  = 18.01      )                                         

      REAL ESFAC, ERFAC
      PARAMETER ( ESFAC = H2OMW/AIRMW       )
      PARAMETER ( ERFAC = (1.0-ESFAC)/ESFAC )

      real aw0, aw1, aw2, aw3, aw4, aw5, aw6
      real bw0, bw1, bw2, bw3, bw4, bw5, bw6
      real ai0, ai1, ai2, ai3, ai4, ai5, ai6
      real bi0, bi1, bi2, bi3, bi4, bi5, bi6

      real d0, d1, d2, d3, d4, d5, d6
      real e0, e1, e2, e3, e4, e5, e6
      real f0, f1, f2, f3, f4, f5, f6
      real g0, g1, g2, g3, g4, g5, g6

! ********************************************************
! ***  Polynomial Coefficients WRT Water (Lowe, 1977) ****
! ***              (Valid +50 C to -50 C)             ****
! ********************************************************

      parameter ( aw0 =  6.107799961e+00 * esfac )
      parameter ( aw1 =  4.436518521e-01 * esfac )
      parameter ( aw2 =  1.428945805e-02 * esfac )
      parameter ( aw3 =  2.650648471e-04 * esfac )
      parameter ( aw4 =  3.031240396e-06 * esfac )
      parameter ( aw5 =  2.034080948e-08 * esfac )
      parameter ( aw6 =  6.136820929e-11 * esfac )

      parameter ( bw0 = +4.438099984e-01 * esfac )
      parameter ( bw1 = +2.857002636e-02 * esfac )
      parameter ( bw2 = +7.938054040e-04 * esfac )
      parameter ( bw3 = +1.215215065e-05 * esfac )
      parameter ( bw4 = +1.036561403e-07 * esfac )
      parameter ( bw5 = +3.532421810e-10 * esfac )
      parameter ( bw6 = -7.090244804e-13 * esfac )


! ********************************************************
! ***   Polynomial Coefficients WRT Ice  (Lowe, 1977) ****
! ***              (Valid  +0 C to -50 C)             ****
! ********************************************************

      parameter ( ai0 = +6.109177956e+00 * esfac )
      parameter ( ai1 = +5.034698970e-01 * esfac )
      parameter ( ai2 = +1.886013408e-02 * esfac )
      parameter ( ai3 = +4.176223716e-04 * esfac )
      parameter ( ai4 = +5.824720280e-06 * esfac )
      parameter ( ai5 = +4.838803174e-08 * esfac )
      parameter ( ai6 = +1.838826904e-10 * esfac )

      parameter ( bi0 = +5.030305237e-01 * esfac )
      parameter ( bi1 = +3.773255020e-02 * esfac )
      parameter ( bi2 = +1.267995369e-03 * esfac )
      parameter ( bi3 = +2.477563108e-05 * esfac )
      parameter ( bi4 = +3.005693132e-07 * esfac )
      parameter ( bi5 = +2.158542548e-09 * esfac )
      parameter ( bi6 = +7.131097725e-12 * esfac )


! ********************************************************
! ***         Polynomial Coefficients WRT Ice         ****
! ***   Starr and Cox (1985) (Valid -40 C to -70 C)   ****
! ********************************************************


      parameter ( d0 = 0.535098336e+01 * esfac )
      parameter ( d1 = 0.401390832e+00 * esfac )
      parameter ( d2 = 0.129690326e-01 * esfac )
      parameter ( d3 = 0.230325039e-03 * esfac )
      parameter ( d4 = 0.236279781e-05 * esfac )
      parameter ( d5 = 0.132243858e-07 * esfac )
      parameter ( d6 = 0.314296723e-10 * esfac )

      parameter ( e0 = 0.469290530e+00 * esfac )
      parameter ( e1 = 0.333092511e-01 * esfac )
      parameter ( e2 = 0.102164528e-02 * esfac )
      parameter ( e3 = 0.172979242e-04 * esfac )
      parameter ( e4 = 0.170017544e-06 * esfac )
      parameter ( e5 = 0.916466531e-09 * esfac )
      parameter ( e6 = 0.210844486e-11 * esfac )


! ********************************************************
! ***         Polynomial Coefficients WRT Ice         ****
! ***   Starr and Cox (1985) (Valid -65 C to -95 C)   ****
! ********************************************************

      parameter ( f0 = 0.298152339e+01 * esfac )
      parameter ( f1 = 0.191372282e+00 * esfac )
      parameter ( f2 = 0.517609116e-02 * esfac )
      parameter ( f3 = 0.754129933e-04 * esfac )
      parameter ( f4 = 0.623439266e-06 * esfac )
      parameter ( f5 = 0.276961083e-08 * esfac )
      parameter ( f6 = 0.516000335e-11 * esfac )

      parameter ( g0 = 0.312654072e+00 * esfac )
      parameter ( g1 = 0.195789002e-01 * esfac )
      parameter ( g2 = 0.517837908e-03 * esfac )
      parameter ( g3 = 0.739410547e-05 * esfac )
      parameter ( g4 = 0.600331350e-07 * esfac )
      parameter ( g5 = 0.262430726e-09 * esfac )
      parameter ( g6 = 0.481960676e-12 * esfac )

      REAL        TMAX, TICE
      PARAMETER ( TMAX=323.15, TICE=273.16)
      
      REAL T, D, W, QX, DQX
      T = MIN(TT,TMAX) - TICE
      DQX = 0.
      QX  = 0.

! Fitting for temperatures above 0 degrees centigrade
! ---------------------------------------------------
      if(t.gt.0.) then
       qx = aw0+T*(aw1+T*(aw2+T*(aw3+T*(aw4+T*(aw5+T*aw6)))))
      if (ldqdt)  then
      dqx = bw0+T*(bw1+T*(bw2+T*(bw3+T*(bw4+T*(bw5+T*bw6)))))
      endif
      endif

! Fitting for temperatures between 0 and -40
! ------------------------------------------
      if( t.le.0. .and. t.gt.-40.0 ) then
        w = (40.0 + t)/40.0
       qx =     w *(aw0+T*(aw1+T*(aw2+T*(aw3+T*(aw4+T*(aw5+T*aw6))))))  &
          + (1.-w)*(ai0+T*(ai1+T*(ai2+T*(ai3+T*(ai4+T*(ai5+T*ai6))))))
      if (ldqdt)  then
      dqx =     w *(bw0+T*(bw1+T*(bw2+T*(bw3+T*(bw4+T*(bw5+T*bw6))))))  &
          + (1.-w)*(bi0+T*(bi1+T*(bi2+T*(bi3+T*(bi4+T*(bi5+T*bi6))))))
      endif
      endif

! Fitting for temperatures between -40 and -70
! --------------------------------------------
      if( t.le.-40.0 .and. t.ge.-70.0 ) then
       qx = d0+T*(d1+T*(d2+T*(d3+T*(d4+T*(d5+T*d6)))))
      if (ldqdt) then
      dqx = e0+T*(e1+T*(e2+T*(e3+T*(e4+T*(e5+T*e6)))))
      endif
      endif

! Fitting for temperatures less than -70
! --------------------------------------
      if(t.lt.-70.0) then
       qx = f0+t*(f1+t*(f2+t*(f3+t*(f4+t*(f5+t*f6)))))
      if (ldqdt) then
      dqx = g0+t*(g1+t*(g2+t*(g3+t*(g4+t*(g5+t*g6)))))
      endif
      endif

! Compute Saturation Specific Humidity
! ------------------------------------
      D = (P-ERFAC*QX)
      IF(D.LT.0.) THEN
       Q = 1.0
       IF (LDQDT)  DQDT = 0.
      ELSE
       D = 1.0 / D
       Q = MIN(QX * D,1.0)
       IF (LDQDT)  DQDT = (1.0 + ERFAC*Q) * D * DQX
      ENDIF
      RETURN

      END SUBROUTINE qsat
!EOC
!------------------------------------------------------------------------------
!               GEOS-Chem Global Chemical Transport Model
!------------------------------------------------------------------------------
!BOP
!
! !DESCRIPTION: Rearaanges a 2-d grid into 1-d for RAS
!\\
!\\
!
  SUBROUTINE REARRANGE_GRID(INGRID, OUTGRID)
!
! !USES:
!
    USE PRECISION_MOD
    USE CMN_SIZE_MOD
      ! 
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: INGRID(:,:)
! 
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: OUTGRID(:)
!
! !REVISION HISTORY:
! 28 Feb 2017 - K. Yu - Initial Version
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER :: I, J

    DO J = 1, JJPAR
    DO I = 1, IIPAR
       OUTGRID(((J-1)*IIPAR) + I) = INGRID(I,J)
    ENDDO
    ENDDO

  END SUBROUTINE REARRANGE_GRID
!EOC
!------------------------------------------------------------------------------
!               GEOS-Chem Global Chemical Transport Model
!------------------------------------------------------------------------------
!BOP
!
! !DESCRIPTION: Rearaanges a 2-d grid into 1-d for RAS
!\\
!\\
!
  SUBROUTINE REARRANGE_GRID_INT(INGRID, OUTGRID)
!
! !USES:
!
    USE PRECISION_MOD
    USE CMN_SIZE_MOD
! 
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)   :: INGRID(:,:)
! 
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT)  :: OUTGRID(:)
!
! !REVISION HISTORY:
! 28 Feb 2017 - K. Yu - Initial Version
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER               :: I, J

    DO J = 1, JJPAR
    DO I = 1, IIPAR
       OUTGRID(((J-1)*IIPAR) + I) = INGRID(I,J)
    ENDDO
    ENDDO

  END SUBROUTINE REARRANGE_GRID_INT
!EOC
!------------------------------------------------------------------------------
!               GEOS-Chem Global Chemical Transport Model
!------------------------------------------------------------------------------
!BOP
!
! !DESCRIPTION: Rearaanges a 1-d grid into 2-d for RAS
!\\
!\\
!
  SUBROUTINE REARRANGE_2D(INGRID, OUTGRID)
!
! !USES:
!
    USE PRECISION_MOD
    USE CMN_SIZE_MOD
! 
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: INGRID(:)
! 
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: OUTGRID(:,:)
!
! !REVISION HISTORY:
! 28 Feb 2017 - K. Yu - Initial Version
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER :: I, J

    DO J = 1, JJPAR
    DO I = 1, IIPAR
       OUTGRID(I,J) = INGRID(((J-1)*IIPAR) + I)
    ENDDO
    ENDDO

  END SUBROUTINE REARRANGE_2D
!EOC
END MODULE RAS_MOD 
