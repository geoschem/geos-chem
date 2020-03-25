#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag3.F90
!
! !DESCRIPTION: Subroutine DIAG3 prints out diagnostics to the BINARY PUNCH
!  format file.
!\\
!\\
! !INTERFACE:
!
SUBROUTINE DIAG3( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
  USE BPCH2_MOD                          ! For binary punch I/O routines
  USE CMN_DIAG_MOD                       ! Diagnostic switches & arrays
  USE CMN_FJX_MOD                        ! Fast-JX flux diagnostics
  USE CMN_O3_MOD                         ! FMOL
  USE CMN_SIZE_MOD,   ONLY : NDSTBIN
  USE DEPO_MERCURY_MOD                   ! For offline Hg simulation
  USE DIAG_MOD                           ! For diagnostic arrays
  USE DIAG03_MOD                         ! For Hg diagnostic
  USE DIAG53_MOD                         ! For POPs diag
  USE DRYDEP_MOD                         ! For dry deposition
  USE ErrCode_Mod
  USE ERROR_MOD,      ONLY : ERROR_STOP
  USE FILE_MOD
  USE HCO_TYPES_MOD, ONLY : DiagnCont
  USE HCO_DIAGN_MOD
  USE HCO_ERROR_MOD
  USE HCO_INTERFACE_MOD
  USE Input_Opt_Mod,  ONLY : OptInput
  USE TIME_MOD
  USE PhysConstants,  ONLY : AVO         ! Avogadro's #
  USE Precision_Mod                      ! For GEOS-Chem Precision (fp)
  USE Species_Mod,    ONLY : Species
  USE State_Chm_Mod,  ONLY : ChmState
  USE State_Chm_Mod,  ONLY : Ind_
  USE State_Grid_Mod, ONLY : GrdState
  USE State_Met_Mod,  ONLY : MetState
  USE WETSCAV_MOD                        ! For wet deposition
#ifdef APM
  ! Modules from GeosApm directory
  USE APM_DRIV_MOD, ONLY : IFTEMPOUT
  USE APM_DRIV_MOD, ONLY : TEMPOUT
  USE APM_DRIV_MOD, ONLY : NTEMPOUT
  USE APM_DRIV_MOD, ONLY : NPOUTSTEPS
  USE APM_DRIV_MOD, ONLY : APM_RADFOUT
  USE APM_INIT_MOD, ONLY : IFRADF
#endif
#ifdef TOMAS
  USE TOMAS_MOD, ONLY : ICOMP, IDIAG, IBINS  !(win, 1/25/10)
#endif

  IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
  TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
  TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
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
  INTEGER            :: NA, nAdvect, nDryDep, NW, NFAM
  INTEGER            :: I, IREF, J, JREF, L, M, MM, MMB, Levs
  INTEGER            :: N, NN, NMAX, NTEST, N_TOT_TRC, T
  INTEGER            :: IE, IN, IS, IW, ITEMP(3)
  INTEGER            :: NSPECOUT

  INTEGER            :: NN1    ! TOMAS tracers
  INTEGER            :: NBIN   ! TOMAS bin counter (win, 1/25/10)

  REAL(fp)           :: SCALE_TMP(State_Grid%NX,State_Grid%NY)
  REAL(fp)           :: SCALE_A3
  REAL(fp)           :: SCALED,     SCALEDYN
  REAL(fp)           :: SCALECONV,  SCALESRCE,  SCALECHEM
  REAL(fp)           :: SCALEDIAG,  SCALE_ND66, SCALE_ND67
  REAL(fp)           :: SCALERAD
  REAL(fp)           :: SCALEX,     SECONDS,    PMASS
  REAL(fp)           :: PRESSX,     FDTT,       AREA_M2
  REAL(fp)           :: SCALE_I3
  REAL(f8)           :: DIAGb,      DIAGe

  ! For binary punch file, version 2.0
  CHARACTER (LEN=40) :: CATEGORY
  REAL(f4)           :: ARRAY(State_Grid%NX,State_Grid%NY,State_Grid%NX+1)
  REAL(f4)           :: LONRES, LATRES
  INTEGER            :: IFIRST, JFIRST, LFIRST
  INTEGER            :: HALFPOLAR
  INTEGER, PARAMETER :: CENTER180 = 1
  CHARACTER (LEN=20) :: MODELNAME
  CHARACTER (LEN=40) :: UNIT
  CHARACTER (LEN=40) :: RESERVED = ''

#ifdef TOMAS
  ! For ND06 diagnostics
  CHARACTER(LEN=1)   :: ISTR1
  CHARACTER(LEN=2)   :: ISTR2

  ! For ND60 TOMAS diagnostic, avoids an array temporary (bmy, 1/28/14)
  REAL(f4)           :: ARR2D(State_Grid%NY,State_Grid%NZ)
#endif

  ! Pointers
  ! We need to define local arrays to hold corresponding values
  ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
  REAL(fp), POINTER        :: Spc(:,:,:,:)
  REAL(fp), POINTER        :: AD(:,:,:)

  ! Interface w/ HEMCO diagnostics
  INTEGER                  :: FLAG
  INTEGER                  :: AFill
  CHARACTER(LEN= 63)       :: DiagnName, SrcName, SpcName, FullName
  CHARACTER(LEN=255)       :: MSG
  TYPE(DiagnCont), POINTER :: DiagnCnt
  REAL(fp)                 :: FACTOR, MW_G
  REAL(fp), PARAMETER      :: GPERKG   = 1000e+0_fp
  REAL(fp), PARAMETER      :: MWC      = 12e+0_fp ! hard-coded MW
  REAL(fp), PARAMETER      :: CM2PERM2 = 10000e+0_fp
  REAL(fp), PARAMETER      :: S_DMS    = 32e+0_fp / 62e+0_fp
  REAL(fp), PARAMETER      :: S_SO2    = 32e+0_fp / 64e+0_fp
  REAL(fp), PARAMETER      :: S_SO4    = 32e+0_fp / 96e+0_fp

  CHARACTER(LEN=255), PARAMETER :: LOC = 'DIAG3 (diag3.F)'

  ! To point to the species database
  CHARACTER(LEN=31)        :: Name
  TYPE(Species), POINTER   :: SpcInfo

  ! Now define local tracer flags so that we can remove these from
  ! tracerid_mod.F to facilitate FlexChem implementation (bmy, 5/2/16)
  LOGICAL, SAVE            :: FIRST = .TRUE.
  INTEGER, SAVE            :: id_Rn222
  INTEGER, SAVE            :: id_Pb210,    id_Pb210Strat
  INTEGER, SAVE            :: id_Be7,      id_Be7Strat
  INTEGER, SAVE            :: id_Be10,     id_Be10Strat
  INTEGER, SAVE            :: id_POPG
  INTEGER, SAVE            :: id_POPPOCPO, id_POPPOCPI
  INTEGER, SAVE            :: id_POPPBCPO, id_POPPBCPI
  INTEGER, SAVE            :: id_DST1,  id_DST2,  id_DST3,  id_DST4
  INTEGER, SAVE            :: id_DAL1,  id_DAL2,  id_DAL3,  id_DAL4
  INTEGER, SAVE            :: id_BCPI,  id_OCPI,  id_POA1,  id_MTPA
  INTEGER, SAVE            :: id_LIMO,  id_MTPO,  id_TSOA1
  INTEGER, SAVE            :: id_ASOA1, id_OPOA1, id_OPOG1
  INTEGER, SAVE            :: id_SALA,  id_SALC,  id_MOPO
  INTEGER, SAVE            :: id_MOPI,  id_DMS,   id_SO2,   id_SO4
  INTEGER, SAVE            :: id_NH3,   id_NO,    id_CO,    id_ALK4
  INTEGER, SAVE            :: id_ACET,  id_MEK,   id_ALD2,  id_PRPE
  INTEGER, SAVE            :: id_C3H8,  id_CH2O,  id_C2H6,  id_CH4
  INTEGER, SAVE            :: id_ISOP,  id_C2H4,  id_CHBR3, id_BR2
  INTEGER, SAVE            :: id_DUST1, id_NK1,   id_SF1,   id_SS1
  INTEGER, SAVE            :: id_ECIL1, id_ECOB1, id_OCIL1, id_OCOB1
  INTEGER, SAVE            :: id_CH2BR2
  INTEGER, SAVE            :: id_PAN,   id_HNO3,  id_EOH,   id_MGLY
  INTEGER, SAVE            :: id_BENZ,  id_TOLU,  id_XYLE,  id_MOH
  INTEGER, SAVE            :: id_NAP,   id_POG1,  id_POG2
!
!******************************************************************************
!  DIAG3 begins here!
!
!  Define scale factors for division.
!  Add a small number (e.g. 1d-32) to prevent division by zero errors.
!******************************************************************************
!
  ! Assume success
  RC         = GC_SUCCESS

  ! Number of advected species
  nAdvect    = State_Chm%nAdvect

  ! Number of dry-deposited species
  nDryDep    = State_Chm%nDryDep

  ! Initialize
  SpcInfo    => NULL()
  DiagnCnt   => NULL()

  ! Now use counter variables from "time_mod.f" (bmy, 3/27/03)
  DIAGb      = GET_DIAGb()
  DIAGe      = GET_DIAGe()
  SECONDS    = ( DIAGe - DIAGb ) * 3600e+0_fp
  SCALED     = 1e+0_fp
  SCALEDYN   = DBLE( GET_CT_DYN()  ) + 1e-32_fp
  SCALECONV  = DBLE( GET_CT_CONV() ) + 1e-32_fp
  SCALESRCE  = DBLE( GET_CT_EMIS() ) + 1e-32_fp
  SCALECHEM  = DBLE( GET_CT_CHEM() ) + 1e-32_fp
  SCALERAD   = DBLE( GET_CT_RAD()  ) + 1e-32_fp
  SCALE_A3   = DBLE( GET_CT_A3()   ) + 1e-32_fp
  SCALE_I3   = DBLE( GET_CT_I3()   ) + 1e-32_fp
  SCALEDIAG  = DBLE( GET_CT_DIAG() ) + 1e-32_fp
!
!******************************************************************************
!  Now define local tracer flags for certain specialty simulations
!  so that we can remove them from tracerid_mod.F (bmy, 5/2/16)
!******************************************************************************
!
  IF ( FIRST ) THEN

     ! Initialize
     id_POPG      = Ind_('POPG'    )
     id_POPPOCPO  = Ind_('POPPOCPO')
     id_POPPOCPI  = Ind_('POPPOCPI')
     id_POPPBCPO  = Ind_('POPPBCPO')
     id_POPPBCPI  = Ind_('POPPBCPI')
     id_DST1      = Ind_('DST1'    )
     id_DST2      = Ind_('DST2'    )
     id_DST3      = Ind_('DST3'    )
     id_DST4      = Ind_('DST4'    )
     id_DAL1      = Ind_('DSTAL1'  )
     id_DAL2      = Ind_('DSTAL2'  )
     id_DAL3      = Ind_('DSTAL3'  )
     id_DAL4      = Ind_('DSTAL4'  )
     id_BCPI      = Ind_('BCPI'    )
     id_OCPI      = Ind_('OCPI'    )
     id_POA1      = Ind_('POA1'    )
     id_MTPA      = Ind_('MTPA'    )
     id_LIMO      = Ind_('LIMO'    )
     id_MTPO      = Ind_('MTPO'    )
     id_TSOA1     = Ind_('TSOA1'   )
     id_ASOA1     = Ind_('ASOA1'   )
     id_OPOA1     = Ind_('OPOA1'   )
     id_OPOG1     = Ind_('OPOG1'   )
     id_SALA      = Ind_('SALA'    )
     id_SALC      = Ind_('SALC'    )
     id_MOPO      = Ind_('MOPO'    )
     id_MOPI      = Ind_('MOPI'    )
     id_DMS       = Ind_('DMS'     )
     id_SO2       = Ind_('SO2'     )
     id_SO4       = Ind_('SO4'     )
     id_NH3       = Ind_('NH3'     )
     id_NO        = Ind_('NO'      )
     id_CO        = Ind_('CO'      )
     id_ALK4      = Ind_('ALK4'    )
     id_ACET      = Ind_('ACET'    )
     id_MEK       = Ind_('MEK'     )
     id_ALD2      = Ind_('ALD2'    )
     id_PRPE      = Ind_('PRPE'    )
     id_C3H8      = Ind_('C3H8'    )
     id_CH2O      = Ind_('CH2O'    )
     id_C2H6      = Ind_('C2H6'    )
     id_ISOP      = Ind_('ISOP'    )
     id_C2H4      = Ind_('C2H4'    )
     id_CHBR3     = Ind_('CHBR3'   )
     id_CH2BR2    = Ind_('CH2BR2'  )
     id_BR2       = Ind_('BR2'     )
     id_DUST1     = Ind_('DUST1'   )
     id_NK1       = Ind_('NK1'     )
     id_SF1       = Ind_('SF1'     )
     id_SS1       = Ind_('SS1'     )
     id_ECIL1     = Ind_('ECIL1'   )
     id_ECOB1     = Ind_('ECOB1'   )
     id_OCIL1     = Ind_('OCIL1'   )
     id_OCOB1     = Ind_('OCOB1'   )
     id_PAN       = Ind_('PAN'     )
     id_HNO3      = Ind_('HNO3'    )
     id_MOH       = Ind_('MOH'     )
     id_EOH       = Ind_('EOH'     )
     id_MGLY      = Ind_('MGLY'    )
     id_BENZ      = Ind_('BENZ'    )
     id_TOLU      = Ind_('TOLU'    )
     id_XYLE      = Ind_('XYLE'    )
     id_NAP       = Ind_('NAP'     )
     id_POG1      = Ind_('POG1'    )
     id_POG2      = Ind_('POG2'    )

     ! NOTE: CH4 can be an advected species (in CH4 or UCX-based sims),
     ! or a non-advected species (tropchem, soa, soa-svpoa).  We only
     ! want to print out diagnostics if CH4 is an advected species,
     ! so make sure to use the 'A' flag in the call to Ind_().
     ! (bmy, 6/23/16)
     id_CH4       = Ind_('CH4', 'A')

     ! Reset first-time flag
     FIRST = .FALSE.
  ENDIF
!
!******************************************************************************
!  Setup for binary punch file:
!
!  IFIRST, JFIRST, LFIRST = I, J, L indices of the starting grid box
!  LONRES                 = State_Grid%DX, cast to REAL*4
!  LATRES                 = State_Grid%DY, cast to REAL*4
!******************************************************************************
!
  IFIRST = State_Grid%XMinOffset + 1
  JFIRST = State_Grid%YMinOffset + 1
  LFIRST = 1
  LONRES = State_Grid%DX
  LATRES = State_Grid%DY

  ! Get the proper model name and HALFPOLAR setting for the bpch file
  MODELNAME = GET_MODELNAME( Input_Opt, State_Grid )
  HALFPOLAR = GET_HALFPOLAR()

  ! HEMCO interface: get pointer to HcoState object (of hcoi_gc_main_mod.F90)
  IF ( .NOT. ASSOCIATED(HcoState) ) THEN
     CALL ERROR_STOP ( 'HcoState not defined!', LOC )
  ENDIF

  !****************************************************************************
  !  ND03: Diagnostics from Hg0/Hg2/HgP offline simulation (eck, bmy, 1/20/05)
  !****************************************************************************
  IF ( ND03 > 0 ) THEN
     CALL WRITE_DIAG03( Input_Opt, State_Chm, State_Grid, RC )
  ENDIF

#ifdef TOMAS

  !****************************************************************************
  !  ND44: Drydep flux (molec/cm2/s) and velocity (cm/s) diagnostics
  !
  !   #   : Field    : Quantity           : Units               : Scale factor
  !  -------------------------------------------------------------------------
  !  (1 ) : DRYD-FLX : drydep fluxes      : molec/cm2/s or kg/s : SCALECHEM
  !  (2 ) : DRYD-VEL : drydep velocities  : cm/s                : SCALECHEM
  !****************************************************************************
  IF ( ND44 > 0 ) THEN

     !==============================================================
     ! Drydep fluxes
     !==============================================================

     ! Category name
     CATEGORY = 'DRYD-FLX'

     ! # of drydep flux tracers
     IF ( Input_Opt%ITS_A_TAGO3_SIM .or. &
          Input_Opt%ITS_A_MERCURY_SIM ) THEN
        M = nAdvect
     ELSE
        ! Extend dry dep tracers if TOMAS aerosol is turned on
        IF ( id_NK1 > 0 ) THEN
           M = nDryDep + ( ( ICOMP - IDIAG )* IBINS )
        ELSE
           M = nDryDep
        ENDIF
     ENDIF

     ! Loop over drydep species
     DO N = 1, M

        IF ( Input_Opt%ITS_A_RnPbBe_SIM ) THEN

           ! NOTE: ND44 is now archived in molec/cm2/s for all
           ! simulations, including Rn-Pb-Be. (bmy, 6/16/15)
           UNIT = 'molec/cm2/s'
           NN   = State_Chm%Map_DryDep(N)

        ELSE IF ( Input_Opt%ITS_A_TAGO3_SIM .or. &
                  Input_Opt%ITS_A_MERCURY_SIM ) THEN

           ! Tagged O3 or Tagged Hg
           UNIT = 'molec/cm2/s'
           NN   = N

        ELSE

           ! Other simulations
           UNIT = 'molec/cm2/s'

           ! For extended drydep tracers, assign tracer ID of the TOMAS
           ! aerosol mass (win, 7/14/09)
           IF ( N <= nDryDep ) THEN
              NN  = State_Chm%Map_DryDep(N)
              NN1 = NN
           ELSE
              ! To calculate the id_xxx of the associated
              ! tracer. (ccc, 3/11/10)
              NN  = MOD( N - nDryDep-1, IBINS ) + id_NK1

              ! Tracer number for bpch file
              NN1  = ( N - nDryDep ) + ( id_NK1 + IBINS - 1 )
           ENDIF

        ENDIF

        ! To output only the species asked in input.geos
        ! (ccc, 5/15/09)
        MM  = 1
        MMB = 0
        DO WHILE ( MMB /= NN .AND. MM <= TMAX(44) )
           MMB = TINDEX(44,MM)
           MM  = MM + 1
        ENDDO

        ! Save into ARRAY
        ARRAY(:,:,1) = ( AD44(:,:,N,1) / SCALECHEM )

        ! Write to file
        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY, NN1,      &
                    UNIT,      DIAGb,     DIAGe,    RESERVED, &
                    State_Grid%NX, State_Grid%NY, 1, IFIRST, &
                    JFIRST,    LFIRST,    ARRAY(:,:,1) )

     ENDDO

     !==============================================================
     ! Drydep velocities
     !==============================================================

     ! Category and Unit
     CATEGORY  = 'DRYD-VEL'
     UNIT      = 'cm/s'

     ! # of drydep velocity tracers
     IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
        M = 1
     ELSE IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN  ! Ask Helen
        M = 3
     ELSE
        M = nDryDep
     ENDIF

     ! Loop over drydep tracers
     DO N = 1, M

        NN  = State_Chm%Map_DryDep(N)
        ! To output only the species asked in input.geos
        ! (ccc, 5/15/09)
        MM  = 1
        MMB = 0
        DO WHILE ( MMB /= NN .AND. MM <= TMAX(44) )
           MMB = TINDEX(44,MM)
           MM  = MM + 1
        ENDDO

        ! Tracer number plus GAMAP offset
        ARRAY(:,:,1) = AD44(:,:,N,2) / SCALESRCE

        ! Write to file
        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY, NN,       &
                    UNIT,      DIAGb,     DIAGe,    RESERVED, &
                    State_Grid%NX, State_Grid%NY, 1, IFIRST,  &
                    JFIRST,    LFIRST,    ARRAY(:,:,1) )

     ENDDO
  ENDIF
  
  !****************************************************************************
  !  ND06: Dust aerosol emissions
  !
  !   # : Field    : Description                     : Units      : Scale factor
  !  --------------------------------------------------------------------------
  !  (1)  DUST     : Soil dust (4 different classes) : kg         : 1
  !****************************************************************************
  IF ( ND06 > 0 .and. Input_Opt%LDUST .and. Input_Opt%LEMIS ) THEN

     ! Category & unit string
     UNIT     = 'kg'
     CATEGORY = 'DUSTSRCE'
     FACTOR   = 1.0e+0_fp

     ! Loop over # of dust bins (obtained from CMN_SIZE_mod)
     DO N = 1, NDSTBIN

        !----------------------------------------------------
        ! TOMAS simulations: many dust tracers
        !----------------------------------------------------

        ! Tracer number
        NN = id_DUST1 + ( N - 1 )

        ! Get TOMAS dust string
        IF ( N < 10 )  THEN
           WRITE( ISTR1,'(i1)' ) N
           DiagnName = 'AD06_DUST' // ISTR1
        ELSE
           WRITE( ISTR2,'(i2)' ) N
           DiagnName = 'AD06_DUST' // ISTR2
        ENDIF

        ! Write HEMCO diagnostics to bpch file
        CALL DIAG2BPCH( HcoState, DiagnName, CATEGORY, &
                        UNIT,     NN, 1, -1, .FALSE.,  FACTOR, RC )

     ENDDO !N

     ! Include Dust alkalinity sources   tdf 6/18/2K8
     IF ( Input_Opt%LDSTUP ) THEN
        DO N = 5, NDSTBIN*2

           SELECT CASE ( N )
           CASE ( 5 )
              NN         = id_DAL1
              DiagnName  = 'AD06_DSTAL1'
           CASE ( 6 )
              NN         = id_DAL2
              DiagnName  = 'AD06_DSTAL2'
           CASE ( 7 )
              NN         = id_DAL3
              DiagnName  = 'AD06_DSTAL3'
           CASE ( 8 )
              NN         = id_DAL4
              DiagnName  = 'AD06_DSTAL4'
           END SELECT

           ! Write to HEMCO diagnostics to bpch file
           CALL DIAG2BPCH( HcoState, DiagnName, CATEGORY, &
                           UNIT,     NN, 1, -1, .FALSE.,  FACTOR, RC )
        ENDDO
     ENDIF

  ENDIF

  !****************************************************************************
  !  ND59: Size-resolved primary aerosol emission (win, 8/22/07)
  !        Unit is amount of aerosol per time of simulation e.g. sulfate[=] kg S
  !
  !   # : Field   : Description                 : Units     : Scale Factor
  !  ---------------------------------------------------------------------------
  !  (1 )  NK-EMISS : Size-resolved number emission   : No       : 1
  !  (2 )  SF-EMISS : Size-resolved sulfate emission  : kg S     : 1
  !  (3 )  SS-EMISS : Size-resolved sea-salt emission : kg       : 1
  !  (4 )  ECIL-SRC : Size-resolved H-phillic EC emission   : kg : 1
  !  (5 )  ECOB-SRC : Size-resolved H-phobic EC emission    : kg : 1
  !  (6 )  OCIL-SRC : Size-resolved H-phillic OC emission   : kg : 1
  !  (7 )  OCOB-SRC : Size-resolved H-phobic OC emission    : kg : 1
  !  (8 )  DUST-SRC : Size-resolved dust emission      : kg      : 1
  !****************************************************************************
  IF ( ND59 > 0 ) THEN

     !==============================================================
     ! Size-resolved primary aerosol number emission
     !==============================================================
     UNIT     = 'No.'
     CATEGORY = 'NK-EMISS'
     DO NBIN = 1,IBINS
        N = id_NK1 + NBIN - 1
        DO L = 1, 2
           ARRAY(:,:,L) = AD59_NUMB(:,:,L,NBIN)
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY,  N,        &
                    UNIT,      DIAGb,     DIAGe,     RESERVED, &
                    State_Grid%NX, State_Grid%NY, 2, IFIRST,   &
                    JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     ENDDO

     !==============================================================
     ! Size-resolved primary aerosol sulfate mass emission
     !==============================================================
     UNIT     = 'kg S'
     CATEGORY = 'SF-EMISS'
     DO NBIN = 1,IBINS
        N = id_SF1 + NBIN - 1
        DO L = 1, 2
           ARRAY(:,:,L) = AD59_SULF(:,:,L,NBIN)
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY,  N,        &
                    UNIT,      DIAGb,     DIAGe,     RESERVED, &
                    State_Grid%NX, State_Grid%NY, 2, IFIRST,   &
                    JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     ENDDO

     !==============================================================
     ! Size-resolved primary aerosol sea-salt mass emission
     !==============================================================
     UNIT     = 'kg'
     CATEGORY = 'SS-EMISS'
     DO NBIN = 1,IBINS
        N = id_SS1 + NBIN - 1
        DO L = 1, 2
           ARRAY(:,:,L) = AD59_SALT(:,:,L,NBIN)
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,    &
                     HALFPOLAR, CENTER180, CATEGORY,  N,        &
                     UNIT,      DIAGb,     DIAGe,     RESERVED, &
                     State_Grid%NX, State_Grid%NY, 2, IFIRST,   &
                     JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     ENDDO

     !==============================================================
     ! Size-resolved primary Hydrophillic EC mass emission
     !==============================================================
     UNIT     = 'kg'
     CATEGORY = 'ECIL-SRC'
     DO NBIN = 1,IBINS
        N = id_ECIL1 - 1 + NBIN
        DO L = 1, 2
           ARRAY(:,:,L) = AD59_ECIL(:,:,L,NBIN)
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY,  N,        &
                    UNIT,      DIAGb,     DIAGe,     RESERVED, &
                    State_Grid%NX, State_Grid%NY, 2, IFIRST,   &
                    JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     ENDDO

     !==============================================================
     ! Size-resolved primary Hydrophobic EC mass emission
     !==============================================================
     UNIT     = 'kg'
     CATEGORY = 'ECOB-SRC'
     DO NBIN = 1,IBINS
        N = id_ECOB1 - 1 + NBIN
        DO L = 1, 2
           ARRAY(:,:,L) = AD59_ECOB(:,:,L,NBIN)
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY,  N,        &
                    UNIT,      DIAGb,     DIAGe,     RESERVED, &
                    State_Grid%NX, State_Grid%NY, 2, IFIRST,   &
                    JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     ENDDO

     !==============================================================
     ! Size-resolved primary Hydrophillic OC mass emission
     !==============================================================
     UNIT     = 'kg'
     CATEGORY = 'OCIL-SRC'
     DO NBIN = 1,IBINS
        N = id_OCIL1 - 1 + NBIN
        DO L = 1, 2
           ARRAY(:,:,L) = AD59_OCIL(:,:,L,NBIN)
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY,  N,        &
                    UNIT,      DIAGb,     DIAGe,     RESERVED, &
                    State_Grid%NX, State_Grid%NY, 2, IFIRST,   &
                    JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     ENDDO

     !==============================================================
     ! Size-resolved primary Hydrophobic OC mass emission
     !==============================================================
     UNIT     = 'kg'
     CATEGORY = 'OCOB-SRC'
     DO NBIN = 1,IBINS
        N = id_OCOB1 - 1 + NBIN
        DO L = 1, 2
           ARRAY(:,:,L) = AD59_OCOB(:,:,L,NBIN)
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY,  N,        &
                    UNIT,      DIAGb,     DIAGe,     RESERVED, &
                    State_Grid%NX, State_Grid%NY, 2, IFIRST,   &
                    JFIRST,    LFIRST,    ARRAY(:,:,1:2) )
     ENDDO

     !==============================================================
     ! Size-resolved dust mass emission
     !==============================================================
     UNIT     = 'kg'
     CATEGORY = 'DUST-SRC'
     DO NBIN = 1,IBINS
        N = id_DUST1 - 1 + NBIN
        DO L = 1, 2
           ARRAY(:,:,L) = AD59_DUST(:,:,L,NBIN)
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY,  N,        &
                    UNIT,      DIAGb,     DIAGe,     RESERVED, &
                    State_Grid%NX, State_Grid%NY, 2, IFIRST,   &
                    JFIRST,    LFIRST,    ARRAY(:,:,1:2) )

     ENDDO

  ENDIF

  !****************************************************************************
  !  ND60: TOMAS Aerosol microphysical rate (win, 7/9/07)
  !        Unit is amount of aerosol per time of simulation e.g. sulfate[=] kg S
  !
  !   # : Field   : Description                 : Units     : Scale Factor
  !  ---------------------------------------------------------------------------
  !  (1) TMS-COND : Condensation rate           : no. or kg : NONE
  !  (2) TMS-COAG : Coagulation rate            : no. or kg : NONE
  !  (3) TMS-NUCL : Nucleation rate             : no. or kg : NONE
  !  (4) TMS-AQOX : Aqueous oxidation rate      : no. or kg : NONE
  !  (5) AERO-FIX : Accumulated error fixed     : no. or kg : NONE
  !  (6) TMS-SOA  : SOA condensation rate       : no. or kg : NONE
  !****************************************************************************
  IF ( ND60 > 0 ) THEN

     !==============================================================
     ! Condensation rate
     !==============================================================

     ! Category name
     CATEGORY = 'TMS-COND'
     DO M = 1,TMAX(60)
        NN = TINDEX(60,M)

        SCALEX = 1.e+0_fp

        IF ( ( NN .ge. id_NK1 ) .and. &
             ( NN .lt. id_NK1+IBINS )      ) THEN
           UNIT = 'no.'
        ELSE
           UNIT = 'kg'
        ENDIF

        DO L = 1, LD60
           ARR2D(:,L) = AD60_COND(1,:,L,M) / SCALEX
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY, NN,       &
                    UNIT,      DIAGb,     DIAGe,    RESERVED, &
                    1,         State_Grid%NY,     LD60,     IFIRST, &
                    JFIRST,    LFIRST,    ARR2D(:,1:LD60)  )
     ENDDO

     !==============================================================
     ! Coagulation rate
     !==============================================================

     ! Category name
     CATEGORY = 'TMS-COAG'
     DO M = 1,TMAX(60)
        NN = TINDEX(60,M)

        SCALEX = 1.e+0_fp

        IF ( ( NN .ge. id_NK1 ) .and. &
             ( NN .lt. id_NK1+IBINS )      ) THEN
           UNIT = 'no.'
        ELSE
           UNIT = 'kg'
        ENDIF

        DO L = 1, LD60
           ARR2D(:,L) = AD60_COAG(1,:,L,M) / SCALEX
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY, NN,       &
                    UNIT,      DIAGb,     DIAGe,    RESERVED, &
                    1,         State_Grid%NY,     LD60,     IFIRST, &
                    JFIRST,    LFIRST,    ARR2D(:,1:LD60)  )
     ENDDO

     !==============================================================
     ! Nucleation rate
     !==============================================================

     ! Category name
     CATEGORY = 'TMS-NUCL'
     DO M = 1,TMAX(60)
        NN = TINDEX(60,M)

        SCALEX = 1.e+0_fp

        IF ( ( NN .ge. id_NK1 ) .and. &
             ( NN .lt. id_NK1+IBINS )      ) THEN
           UNIT = 'no.'
        ELSE
           UNIT = 'kg'
        ENDIF

        DO L = 1, LD60
           ARR2D(:,L) = AD60_NUCL(1,:,L,M) / SCALEX
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY, NN,       &
                    UNIT,      DIAGb,     DIAGe,    RESERVED, &
                    1,         State_Grid%NY,     LD60,     IFIRST, &
                    JFIRST,    LFIRST,    ARR2D(:,1:LD60)  )
     ENDDO

     !==============================================================
     ! Aqueous oxidation rate
     !==============================================================

     ! Category name
     CATEGORY = 'TMS-AQOX'
     DO M = 1,TMAX(60)
        NN = TINDEX(60,M)

        SCALEX = 1.e+0_fp

        IF ( ( NN .ge. id_NK1 ) .and. &
             ( NN .lt. id_NK1+IBINS )      ) THEN
           UNIT = 'no.'
        ELSE
           UNIT = 'kg'
        ENDIF

        DO L = 1, LD60
           ARR2D(:,L) = AD60_AQOX(1,:,L,M) / SCALEX
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY, NN,       &
                    UNIT,      DIAGb,     DIAGe,    RESERVED, &
                    1,         State_Grid%NY,     LD60,     IFIRST, &
                    JFIRST,    LFIRST,    ARR2D(:,1:LD60)  )
     ENDDO

     !==============================================================
     ! Error fudging during aerosol microphysics
     !==============================================================

     ! Category name
     CATEGORY = 'AERO-FIX'
     DO M = 1,TMAX(60)
        NN = TINDEX(60,M)

        SCALEX = 1.e+0_fp

        IF ( ( NN .ge. id_NK1 ) .and. &
             ( NN .lt. id_NK1+IBINS )      ) THEN
           UNIT = 'no.'
        ELSE
           UNIT = 'kg'
        ENDIF

        DO L = 1, LD60
           ARR2D(:,L) = AD60_ERROR(1,:,L,M) / SCALEX
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY, NN,       &
                    UNIT,      DIAGb,     DIAGe,    RESERVED, &
                    1,         State_Grid%NY,     LD60,     IFIRST, &
                    JFIRST,    LFIRST,    ARR2D(:,1:LD60)  )
     ENDDO

     !==============================================================
     ! SOA Condensation rate  (win, 9/25/07)
     !==============================================================

     ! Category name
     CATEGORY = 'TMS-SOA'
     DO M = 1,TMAX(60)
        NN = TINDEX(60,M)

        SCALEX = 1.e+0_fp

        IF ( ( NN .ge. id_NK1 ) .and. &
             ( NN .lt. id_NK1+IBINS )      ) THEN
           UNIT = 'no.'
        ELSE
           UNIT = 'kg'
        ENDIF

        DO L = 1, LD60
           ARR2D(:,L) = AD60_SOA(1,:,L,M) / SCALEX
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,   &
                    HALFPOLAR, CENTER180, CATEGORY, NN,       &
                    UNIT,      DIAGb,     DIAGe,    RESERVED, &
                    1,         State_Grid%NY,     LD60,     IFIRST, &
                    JFIRST,    LFIRST,    ARR2D(:,1:LD60)  )
     ENDDO

  ENDIF

  !****************************************************************************
  !  ND61: TOMAS microphysics process rate saved in 3-D (kg/s)
  !
  !  Remark: for aerosol number, the unit 'kg' is used as a fake kg where
  !          1 fake kg = 1 aerosol number.
  !          Molecular weight of aerosol number is 1g/mol and
  !
  !   # : Field   : Description                 : Units     : Scale Factor
  !  --------------------------------------------------------------------------
  !  (1) TOMAS-3D : 3-D microphysics rate       : kg/s      : SCALECHEM
  !****************************************************************************
  IF ( ND61 > 0 ) THEN
     CATEGORY = 'TOMAS-3D'

     DO M = 1, TMAX(61)
        NN = TINDEX(61,M)

        SCALEX = SCALECHEM

        UNIT = 'cm-3s-1'

        DO L = 1, LD61
           ARRAY(:,:,L) = AD61(:,:,L,NN) / SCALEX
        ENDDO

        CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     &
                    HALFPOLAR, CENTER180, CATEGORY, NN,         &
                    UNIT,      DIAGb,     DIAGe,    RESERVED,   &
                    State_Grid%NX, State_Grid%NY, LD61, IFIRST, &
                    JFIRST,    LFIRST,    ARRAY(:,:,1:LD61)  )
     ENDDO
  ENDIF
#endif

#ifdef RRTMG

  !****************************************************************************
  !  ND72: RRTMG radiation fields
  !
  !   # : Field  : Description                      : Units    : Scale
  !   factor
  !  -----------------------------------------------------------------------
  !  (1 ) ALLTOASW  : All-sky TOA SW (Total)        : W/m2     : SCALERAD
  !  (2 ) ALLSRFSW  : All-sky Surface SW (Total)    : W/m2     : SCALERAD
  !  (3 ) ALLTOALW  : All-sky TOA LW (Total)        : W/m2     : SCALERAD
  !  (4 ) ALLSRFLW  : All-sky Surface LW (Total)    : W/m2     : SCALERAD
  !  (5 ) CLRTOASW  : Clear-sky TOA SW (Total)      : W/m2     : SCALERAD
  !  (6 ) CLRSRFSW  : Clear-sky Surface SW (Total)  : W/m2     : SCALERAD
  !  (7 ) CLRTOALW  : Clear-sky TOA LW (Total)      : W/m2     : SCALERAD
  !  (8 ) CLRSRFLW  : Clear-sky Surface LW (Total)  : W/m2     : SCALERAD
  !****************************************************************************
  IF ( Input_Opt%LRAD .and. ND72 > 0 ) THEN
     CATEGORY = 'RADMAP-$'
     IF ( Input_Opt%amIRoot ) &
          WRITE(6,*) 'Input_Opt%NWVSELECT',Input_Opt%NWVSELECT

     ! ND72 is updated every rad timestep
     SCALEX = SCALERAD

     ! Number of output species plus baseline
     IF ( Input_Opt%LUCX ) THEN
        NSPECOUT=Input_Opt%NSPECRADMENU+1
     ELSE
        NSPECOUT=Input_Opt%NSPECRADMENU
     ENDIF

     DO M = 1, TMAX(72)
        N  = TINDEX(72,M)

        IF ( N > PD72R ) CYCLE
        !only output clear-sky and all-sky if they are switched on
        !or if we're doing the optics
        !but only output 2nd and 3rd sets of optics if they're
        !requested by user in input.geos.rad
        IF (((N.LE.4 ).AND.(Input_Opt%LSKYRAD(2).EQV..TRUE.)).OR. &
            ((N.GT.4 ).AND.(N.LE.8).AND. &
                          (Input_Opt%LSKYRAD(1).EQV..TRUE.)).OR. &
            ((N.GT.8 ).AND.(N.LE.11)).OR. &                  !1st set of optics
            ((N.GT.11).AND.(N.LE.14).AND. &
                           (Input_Opt%NWVSELECT.GT.1)).OR. &   !2nd set
            ((N.GT.14).AND.(N.LE.17).AND. &
                           (Input_Opt%NWVSELECT.GT.2))) THEN !3rd set

           ! Select proper unit string (cf list above)
           IF (N.LE.8) THEN
              UNIT = 'W/m2'
           ELSE
              !AOD, SSA, ASYM
              UNIT = 'UNITLESS'
           ENDIF

           ! each case is for different output type, within each
           ! of those are the values for each species in the RAD input
           ! menu. Only output species that are switched on for each
           ! output type selected in the ND72 menu
           DO IS=1,NSPECOUT
              NN=(N-1)*NSPECOUT+IS

              !if output is a flux...
              IF (N.LE.8) THEN
                 !always output the baseline flux
                 IF (IS.EQ.1) THEN
                    ARRAY(:,:,1) = AD72(:,:,NN) / SCALEX
                    CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,   &
                                HALFPOLAR, CENTER180, CATEGORY, NN,       &
                                UNIT,      DIAGb,     DIAGe,    RESERVED, &
                                State_Grid%NX, State_Grid%NY, 1, IFIRST,  &
                                JFIRST,    LFIRST,    ARRAY(:,:,1) )
                 !for other outputs check species has been selected
                 ELSE IF (IS.GE.2) THEN
                    IF (Input_Opt%LSPECRADMENU(IS-1).EQ.1) THEN
                       DO I=1,State_Grid%NX
                       DO J=1,State_Grid%NY
                          !store as difference in flux from baseline
                          ARRAY(I,J,1) = (AD72(I,J,(N-1)*NSPECOUT+1) - &
                                          AD72(I,J,NN))/SCALEX
                       ENDDO
                       ENDDO

                       CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,   &
                                   HALFPOLAR, CENTER180, CATEGORY, NN,       &
                                   UNIT,      DIAGb,     DIAGe,    RESERVED, &
                                   State_Grid%NX, State_Grid%NY, 1, IFIRST,  &
                                   JFIRST,    LFIRST,    ARRAY(:,:,1) )
                    ENDIF
                 ENDIF
                 ! N > 8 so we are outputting optics
                 ! but only for turned on species and only for aerosols
                 ! i.e. skip IS=2 and IS=3
              ELSE IF (IS.GE.4) THEN
                 IF (Input_Opt%LSPECRADMENU(IS-1).EQ.1) THEN
                    !output optics
                    ARRAY(:,:,1) = AD72(:,:,NN) / SCALEX
                    CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,   &
                                HALFPOLAR, CENTER180, CATEGORY, NN,       &
                                UNIT,      DIAGb,     DIAGe,    RESERVED, &
                                State_Grid%NX, State_Grid%NY, 1, IFIRST,  &
                                JFIRST,    LFIRST,    ARRAY(:,:,1) )
                 ENDIF
              ENDIF !flux vs optics check
           ENDDO
        ENDIF !all-sky, clear-sky check
     ENDDO
  ENDIF
#endif

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diag2bpch
!
! !DESCRIPTION: Wrapper routine to get diagnostics from HEMCO and write
!  them to bpch.  This will look up diagnostics 'dName' and write the
!  corresponding diagnostics array to the bpch output file.
!\\
!\\
!  NOTE: This is a "bridge" routine intended to provide backwards compatibility
!  with the existing GEOS-Chem diagnostics.  Eventually we will write all
!  GEOS-Chem diagnostics to netCDF format but we are not yet at that point.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DIAG2BPCH( HcoState, dname, bcat, bUnit, bN, &
                        dAF, dZ, dOPTIONAL, dFACTOR, ERR, AreaScal )
!
! !USES:
!
    USE HCO_STATE_MOD,  ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState  ! HEMCO State object
    CHARACTER(LEN=*), INTENT(IN)    :: dname     ! Diagnostics name
    CHARACTER(LEN=*), INTENT(IN)    :: bCat      ! BPCH category
    CHARACTER(LEN=*), INTENT(IN)    :: bUnit     ! BPCH units
    INTEGER,          INTENT(IN)    :: bN        ! BPCH ID
    INTEGER,          INTENT(IN)    :: dAF       ! AutoFill of diagnostic
    INTEGER,          INTENT(IN)    :: dZ        ! # of vertical levels
    LOGICAL,          INTENT(IN)    :: dOPTIONAL ! Optional field?
    REAL(fp),         INTENT(IN)    :: dFACTOR   ! Scale factor
    INTEGER,          OPTIONAL      :: AreaScal  ! Area scaling
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: ERR       ! Return code
!
! !REMARKS:
!  (1) Data is multiplied by factor dFACTOR before writing to disk.
!  (2) dAF determines the autofill flag used for the diagnostics lookup.
!  (3) dZ are the number of vertical levels. Set to -1 for 2D arrays,
!       otherwise a 3D array with exactly that many levels is written.
!  (4) dOptional determines whether or not this is an optional diagnostics.
!       For optional diagnostics, zeros are written out if no corresponding
!       HEMCO diagnostics could be found (because it is undefined or because
!       they have never been updated). For non-optional diagnostics, the
!       routine stops with an error if no HEMCO diagnostics is found.
!  (5) AScal can be used to multiply (1) or divide (-1) the data by the grid
!       area.
!
! !REVISION HISTORY:
!  06 Aug 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    CHARACTER(LEN=255), PARAMETER :: LC = 'DIAG2BPCH (diag3.F)'
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                       :: II, FLG, LEVS, ASCL
    CHARACTER(LEN=255)            :: MG

    ! Pointers
    TYPE(DiagnCont),    POINTER   :: DgnCont

    !=================================================================
    ! DIAG2BPCH begins here!
    !=================================================================

    ! Initialize
    ARRAY    = 0e0
    DgnCont => NULL()

    ! Set default for area scaling
    ASCL = 0
    IF ( PRESENT(AreaScal) ) ASCL = AreaScal

    ! Get diagnostics
    DgnCont => NULL()
    CALL Diagn_Get( HcoState, .FALSE., DgnCont, &
                    FLG,   ERR,      cName=TRIM(dname), &
                    AutoFill=dAF, COL=HcoState%Diagn%HcoDiagnIDManual)

    ! If diagnostic not found, write error message
    IF ( ERR /= HCO_SUCCESS ) THEN
       MG = 'Cannot find diagnostics ' // TRIM(dname)
       CALL ERROR_STOP ( MG, LC )
    ENDIF

    ! Vertical levels
    levs = MAX(dz,1)

    ! Save into ARRAY. Apply scale factor if required
    IF ( FLG == HCO_SUCCESS ) THEN
       IF ( dz > 0 ) THEN
          ARRAY(:,:,1:levs) = DgnCont%Arr3D%Val(:,:,1:levs) * dFACTOR
       ELSE
          ARRAY(:,:,1) = DgnCont%Arr2D%Val(:,:) * dFACTOR
       ENDIF
    ELSE
       MG = 'No diagnostics returned: ' // TRIM(dname)
       IF ( dOptional ) THEN
          MG = TRIM(MG) // ' - will write zeros!'
          CALL HCO_WARNING ( MG, ERR, THISLOC=LC )
          ARRAY(:,:,1:levs) = 0.0
       ELSE
          CALL ERROR_STOP ( MG, LC )
       ENDIF
    ENDIF

    ! Eventually scale by area
    IF ( ASCL == 1 ) THEN
       DO II=1,levs
          ARRAY(:,:,II) = ARRAY(:,:,II) * HcoState%Grid%AREA_M2%Val(:,:)
       ENDDO
    ELSEIF ( ASCL == -1 ) THEN
       DO II=1,levs
          ARRAY(:,:,II) = ARRAY(:,:,II) / HcoState%Grid%AREA_M2%Val(:,:)
       ENDDO
    ENDIF

    ! Write to bpch file
    CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,         &
                HALFPOLAR, CENTER180, bCat,     bN,             &
                bUnit,     DIAGb,     DIAGe,    RESERVED,       &
                State_Grid%NX, State_Grid%NY, LEVS,     IFIRST, &
                JFIRST,    LFIRST,    ARRAY(:,:,1:levs) )

    ! Return w/ success
    ERR = HCO_SUCCESS

    ! Nullify pointer
    DgnCont => NULL()

  END SUBROUTINE DIAG2BPCH
!EOC
END SUBROUTINE DIAG3
#endif
