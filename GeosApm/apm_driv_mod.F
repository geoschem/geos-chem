!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: apm_driv_mod
!
! !DESCRIPTION: Module APM\_DRIV\_MOD contains variables and routines to drive 
!  the Advanced Particle Microphysics (APM) model.  It serves as the 
!  interface between APM module and the 3D model.
!\\
!\\
! !INTERFACE:
!
      MODULE APM_DRIV_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !DEFINED PARAMETERS:
!
      INTEGER, PARAMETER,   PUBLIC :: NTEMPOUT1 = 1
!
! !PUBLIC DATA MEMBERS:
!
      INTEGER,              PUBLIC :: NPOUTSTEPS
      INTEGER,              PUBLIC :: IFTEMPOUT
      INTEGER,              PUBLIC :: NTEMPOUT
      REAL*8,  ALLOCATABLE, PUBLIC :: FCLOUD(:,:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: GFTOT3D(:,:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: DENWET3D(:,:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PLVSOG(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: TEMPOUT(:,:,:,:)

      ! Bin index for cloud act diameters taking into account uptake of
      ! NIT, NH4, SOA
      INTEGER, ALLOCATABLE, PUBLIC :: IACT1(:,:,:)
      INTEGER, ALLOCATABLE, PUBLIC :: IACT2(:,:,:)
      INTEGER, ALLOCATABLE, PUBLIC :: IACT3(:,:,:) 

      ! H2SO4 gas production rate
      REAL*8,  ALLOCATABLE, PUBLIC :: PSO4GAS(:,:,:) 
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC                       :: APM_DRIV
      PUBLIC                       :: AERONUM
      PUBLIC                       :: INIT_APM3D
      PUBLIC                       :: CLEANUP_APM3D
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE                      :: APM3DQGEOS
      PRIVATE                      :: GET_OH
!
! !REMARKS:
!  The APM model was designed and developed for implementation into GEOS-Chem 
!  by Fangqun Yu and Gan Luo at State University of New York (SUNY) at Albany.
!  (Emails: yfq@asrc.cestm.albany.edu; ganluo@asrc.cestm.albany.edu) 
!                                                                             .
!  Fangqun focused on overall strategy, aerosol structure (bins, compositions, 
!  types, mixing states) design, computationally efficient schemes, particle 
!  microphysics (nucleation, growth, coagulation), model evaluation and 
!  improvement, and model application.
!  Gan focused on overall strategy, integration of APM with GEOS-Chem, model 
!  input, emission, transport, size-resolved dry deposition and wet scavenging,
!  output visualization, and model application.
!                                                                             .
!  The major reference of the APM model implemented in GEOS-Chem is:
!                                                                             .
!  1. Yu, F., and G. Luo, Simulation of particle size distribution with a 
!        global aerosol model: Contribution of nucleation to aerosol and CCN 
!        number concentrations, Atmos. Chem. Phys., 9, 7691-7710, 2009. 
!                                                                             .
!  The APM model is an advanced multi-type, multi-component, size-resolved 
!  microphysics model developed for a wide range of applications. The current 
!  APM model is the result of the past developments and validation efforts 
!  aimed at explaining atmospheric particle observations, reported in a 
!  number of previous publications including:
!                                                                             .
!  2. Turco, R. P., Hamill, P., Toon, O. B., Whitten, R. C., and Kiang, C. S.: 
!        A one-dimensional model describing aerosol formation and evolution in 
!        the stratosphere, Part I, Physical processes and mathematical analogs,
!        J. Atmos. Sci., 36, 699-717, 1979.
!  3. Toon, O. B., Turco, R. P., Westphal, D., Malone, R., and Liu, M. S.: A 
!        multidimensional model for aerosols: Description of computational 
!        analogs, J. Atmos. Sci., 45, 2123-2143, 1988.
!  4. Jacobson, M., Turco, R., Jensen, E. and Toon O.: Modeling coagulation 
!        among particles of different composition and size, Atmos. Environ., 
!        28, 1327-1338, 1994.
!  5. Jacobson, M. Z., and Turco, R. P.: Simulating condensational growth, 
!        evaporation and coagulation of aerosols using a combined moving and 
!        stationary size grid, Aerosol Sci. Tech., 22, 73-92, 1995.
!  6. Yu, F., and R. P. Turco: The role of ions in the formation and evolution 
!        of particles in aircraft plumes, Geophys. Res. Lett., 24, 1927-1930, 
!        1997.
!  7. Yu, F.: A Study of the Formation and Evolution of Aerosols and Contrails 
!        in Aircraft Wakes: Development, Validation and Application of an 
!        Advanced Particle Microphysics (APM) Model, Doctoral Dissertation, 
!        UCLA, 1998.
!  8. Yu, F., From molecular clusters to nanoparticles: Second-generation 
!        ion-mediated nucleation model, Atmos. Chem. Phys.,6, 5193-5211, 2006.
!  9. Yu, F., and R. P. Turco, Case studies of particle formation events 
!        observed in boreal forests: Implications for nucleation mechanisms, 
!        Atmos. Chem. Phys., 8, 6085-6102, 2008.
!                                                                             .
!  While the core components of the present APM model implemented in GEOS-Chem 
!  was written from scratch using fortran 95, some algorithms and ideas were 
!  inherited from the above mentioned references. 
!                                                                             .
!  For more information and updates, please check the APM aerosol microphysics 
!  wiki page 
!  http://wiki.seas.harvard.edu/geos-chem/index.php/APM_aerosol_microphysics.
!
! !REVISION HISTORY: 
!  8/2008 - 10/2010 - F. Yu & G. Luo  - Initial versions
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      REAL*8,        ALLOCATABLE :: XN4D(:,:,:,:)
      REAL*8,  SAVE, ALLOCATABLE :: XQ3D(:,:,:)
      INTEGER,       ALLOCATABLE :: NCOAG3D(:,:,:,:)
      INTEGER,       ALLOCATABLE :: IFOUTIJ(:,:)
      INTEGER,       ALLOCATABLE :: SITEID(:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm_driv
!
! !DESCRIPTION: Subroutine APM\_DRIV is the interface between APM and
!  the GEOS-Chem model.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE APM_DRIV( State_Met, State_Chm )
!
! !USES:
!
      USE APM_INIT_MOD,       ONLY : NGCOND,   NSO4,       NSEA
      USE APM_INIT_MOD,       ONLY : NDSTB,    NTYP
      USE APM_INIT_MOD,       ONLY : NCTSO4,   NCTBCOC,    NCTDST
      USE APM_INIT_MOD,       ONLY : NCTSEA,   NBCOCT
      USE APM_INIT_MOD,       ONLY : NBCPIF,   NBCPIO
      USE APM_INIT_MOD,       ONLY : NOCPIF,   NOCPIO
      USE APM_INIT_MOD,       ONLY : NBCPOF,   NBCPOO
      USE APM_INIT_MOD,       ONLY : NOCPOF,   NOCPOO !luogan
      USE APM_INIT_MOD,       ONLY : IDTSO4G,  IDTSO4BIN1
      USE APM_INIT_MOD,       ONLY : IDTCTSO4, IDTCTBCOC,  IDTCTDST
      USE APM_INIT_MOD,       ONLY : IDTCTSEA, IDTSEABIN1, IDTDSTBIN1
      USE APM_INIT_MOD,       ONLY : IDTBCPIFF
      USE APM_INIT_MOD,       ONLY : IFNUCL,   IFAG
      USE APM_INIT_MOD,       ONLY : XMACID,   XMLVSOG
      USE APM_INIT_MOD,       ONLY : M1ACID,   M1LVSOG
      USE APM_INIT_MOD,       ONLY : MAXSITE,  MSITE
      USE APM_INIT_MOD,       ONLY : ISITES,   JSITES,     LOUT
      USE APM_INIT_MOD,       ONLY : IFSITE,   IFSITEOUT,  IFQANN
      USE APM_COAG_MOD,       ONLY : READCK6DTABLE
      USE APM_PHYS_MOD,       ONLY : APM_PHYS
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GRID_MOD,           ONLY : GET_XMID,  GET_YMID
      USE PBL_MIX_MOD,        ONLY : GET_PBL_TOP_L,GET_PBL_TOP_hPa
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TIME_MOD,           ONLY : GET_TS_CHEM, ITS_A_NEW_MONTH
      USE TIME_MOD,           ONLY : GET_YEAR, GET_MONTH
      USE TIME_MOD,           ONLY : GET_DAY,GET_HOUR,GET_MINUTE
      USE TRACER_MOD,         ONLY : TRACER_MW_KG
      USE TRACERID_MOD,       ONLY : IDTNH4,  IDTNIT,  IDTSO4, IDTNH3
      USE TRACERID_MOD,       ONLY : IDTBCPO, IDTBCPI, IDTOCPO, IDTOCPI
      USE TRACERID_MOD,       ONLY : IDTTSOA1, IDTTSOA2
      USE TRACERID_MOD,       ONLY : IDTTSOA3, IDTTSOA0
      USE TRACERID_MOD,       ONLY : IDTISOA1, IDTISOA2, IDTISOA3
      USE TRACERID_MOD,       ONLY : IDTASOA1,  IDTASOA2, IDTASOA3
      USE TRACERID_MOD,       ONLY : IDTMSA
      USE TRACERID_MOD,       ONLY : IDTDST1,IDTDST2,IDTDST3,IDTDST4
      USE TRACERID_MOD,       ONLY : IDTSO2

      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY: 
!  8/2008 - 10/2010 - F. Yu & G. Luo  - Initial versions
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  17 Oct 2013 - M. Sulprizio- Update SOAT so that it is consistent with tracers
!                              in Havala Pye's updated SOA simulation (G. Luo)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I,J,L,N,SIZENUM,IY,MDAY

      REAL*8  :: PRESS,TK,RHIN,XQ,CACID,PACID,DTAPM,DT
      REAL*8  :: MSO4,MNIT,MNH4,MMSA,SOAT
      REAL*8  :: MBCS, MOCS   ! mass of sulfate attached to BC, OC
      REAL*8  :: MDSTS, MSALTS   ! mass of sulfate attached to dust,sea salt
      REAL*8  :: MSULFT   ! total sulfate
      REAL*8  :: MBCOC8(8)
      REAL*8  :: XM1D(NSO4+NSEA), XN1D(NSO4),TEMPOUT1(NTEMPOUT1)
      REAL*8  :: XMDST(NDSTB)

      REAL*8  :: MASS1, MASS2

      REAL*8  :: VOL, CLVSOG, PLVSOG1
      REAL*8  :: MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV

      REAL*8  :: XOH, XU, XV
      INTEGER :: KYEAR,KMON,KDAY,KHOUR,KMIN,ISITE,JSITE,NSITE
      REAL*8  :: TOP, TOPP
      INTEGER :: KKOUT
      REAL*8  :: XLAT, XLON

      REAL*8  :: CSO2,CNH3,XN0
      REAL*8  :: CSOG1,CSOG2,CSOG3,CSOG4,CSOA1,CSOA2,CSOA3,CSOA4
      REAL*8  :: CSOG5,CSOA5
      REAL*8  :: GFTOT1,GFTOT2,DENWET1,DENWET2
      INTEGER :: IACT10, IACT20, IACT30   ! bin index for cloud act 
                                          ! diameters corresponding to RDRY
      REAL*8  :: FCLOUD1(NSO4+4)
      INTEGER :: NCOAG1,NCOAG2
      REAL*8  :: SITEOUT(MAXSITE,LOUT,35),YTOP(MAXSITE,2)
      REAL*8  :: SITEOUT1(MAXSITE,LOUT,NSO4+NSEA+NDSTB+8)
      REAL*8  :: SITEOUT2(MAXSITE,LOUT,NTEMPOUT1)
         
      LOGICAL, SAVE    :: FIRST = .TRUE.
      LOGICAL, SAVE    :: FIRST1 = .TRUE.

      ! Make a pointer to the tracer array
      REAL*8, POINTER :: STT(:,:,:,:)

      WRITE(6,*)'    - APM calculation '

      ! Point to STT
      STT => State_Chm%TRACERS

      ! Chemistry timestep [s]
      DT = GET_TS_CHEM() * 60d0
      DTAPM = DT  !s

! Calculate 3-D ionization rate
      IF(IFNUCL.EQ.1) THEN
        IF(IFQANN.EQ.1) THEN
         IF(FIRST) THEN   ! only calculate Q once 
          WRITE(6,*)"CALCULATE 3-D ANNUAL MEAN IONIZATION RATES"
          CALL APM3DQGEOS( State_Met )
          FIRST = .FALSE.
         ENDIF
        ENDIF
      ENDIF

! TEMPOUT
      IF(FIRST1)THEN
        CALL  READCK6DTABLE
! Intitialize for first step useage
        GFTOT3D = 1.0
        DENWET3D = 2.0
! init FCLOUD
        FCLOUD=0.D0
        DO N=26,NSO4
        DO L=1,LLPAR
        DO J=1,JJPAR
        DO I=1,IIPAR
          FCLOUD(I,J,L,N)=1/15.D0
        ENDDO
        ENDDO
        ENDDO
        ENDDO

!define the variables you want to output
        IFTEMPOUT=.TRUE.
        IF(IFTEMPOUT)THEN
         NTEMPOUT=NTEMPOUT1 !number of species you want to output
         ALLOCATE( TEMPOUT( IIPAR, JJPAR, LLPAR, NTEMPOUT ) )
         TEMPOUT = 0d0
         NPOUTSTEPS = 0
         FIRST1=.FALSE.
        ENDIF
      ENDIF

      SITEOUT =0.
      SITEOUT1 =0.
      SITEOUT2 =0.

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
!$OMP+PRIVATE( SIZENUM, PRESS, TK, RHIN )
!$OMP+PRIVATE( CACID,PACID )
!$OMP+PRIVATE( MSO4,MNIT,MNH4,SOAT)
!$OMP+PRIVATE( MBCS, MOCS, MSULFT,MDSTS,MSALTS)
!$OMP+PRIVATE( MBCOC8, MMSA)
!$OMP+PRIVATE( XMDST)
!$OMP+PRIVATE( MASS1, MASS2)
!$OMP+PRIVATE( CSO2,CNH3,XN0)
!$OMP+PRIVATE( CSOG1,CSOG2,CSOG3,CSOG4,CSOG5)
!$OMP+PRIVATE( CSOA1,CSOA2,CSOA3,CSOA4,CSOA5)
!$OMP+PRIVATE( VOL)
!$OMP+PRIVATE( CLVSOG,MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV)
!$OMP+PRIVATE( XM1D,XN1D,TEMPOUT1)
!$OMP+PRIVATE( XQ,PLVSOG1,GFTOT1,GFTOT2,DENWET1,DENWET2)
!$OMP+PRIVATE( IACT10,IACT20,IACT30,FCLOUD1)
!$OMP+PRIVATE( NCOAG1,NCOAG2)
!$OMP+PRIVATE( XOH, XU, XV,XLAT, XLON)
!$OMP+PRIVATE( KYEAR,KMON,KDAY,KHOUR,KMIN,ISITE,JSITE,NSITE)
!$OMP+PRIVATE( TOP, TOPP)
!$OMP+PRIVATE( KKOUT)
!$OMP+SCHEDULE( DYNAMIC )

      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        PRESS = GET_PCENTER(I,J,L) * 1.d2         ! P at level center [Pa]
        TK    = State_Met%T(I,J,L)                ! Temperature [K]
        VOL   = State_Met%AIRVOL(I,J,L)           ! Volume of grid box [m3]
        RHIN  = MIN( 99.d0, State_Met%RH(I,J,L) ) ! Cap RH at 99%
        RHIN  = MAX( 0.d0, RHIN )                 ! Safety check

        IF(IFNUCL.EQ.1) THEN   ! IMN
         XQ = XQ3D(I,J,L)   !ion-pairs/cm3s
        ELSE
         XQ = 1.d-20   !ion-pairs/cm3s
        ENDIF

        !kg/(box*timestep) to #/(cm3*s)
        PACID = XMACID/96.* PSO4GAS(I,J,L)/(VOL*M1ACID*DT)*1.d-6

        !kg/box to #/cm3
        CACID = STT(I,J,L,IDTSO4G)/(VOL*M1ACID)*1.d-6

        IF(IFAG.EQ.1) THEN
          PLVSOG1 = PLVSOG(I,J,L)/(M1LVSOG*1.d+15)  !#/cm3s, PLVSOG in ug/m3s, M1LVSOG in kg
          CLVSOG  = STT(I,J,L,IDTSO4G+1)/(VOL*M1LVSOG)*1.d-6 
! coatd LV-SOA
          MSULFLV = STT(I,J,L,IDTCTSO4)/VOL   !kg/m3
          MBCLV   = STT(I,J,L,IDTCTBCOC+2)/VOL !kg/m3
          MOCLV   = STT(I,J,L,(IDTCTBCOC+3))/VOL !kg/m3
          MDSTLV  = STT(I,J,L,IDTCTDST+1)/VOL !kg/m3
          MSALTLV = STT(I,J,L,IDTCTSEA+1)/VOL !kg/m3
        ELSE
          PLVSOG1 = 0.
          CLVSOG  = 0.
          MSULFLV = 0.
          MBCLV   = 0.
          MOCLV   = 0.
          MDSTLV  = 0.
          MSALTLV = 0.
        ENDIF

!------------------------------------------------------------------------------
! Prior to 10/17/13:
!        SOAT= (STT(I,J,L,IDTSOA1)+STT(I,J,L,IDTSOA2)+STT(I,J,L,IDTSOA3)
!     &         + STT(I,J,L,IDTSOA4) + STT(I,J,L,IDTSOA5))/VOL   ! Total SV-&MV-SOA
!------------------------------------------------------------------------------
        SOAT= (STT(I,J,L,IDTTSOA1)+STT(I,J,L,IDTTSOA2)+
     &         STT(I,J,L,IDTTSOA3)+STT(I,J,L,IDTTSOA0)+
     &         STT(I,J,L,IDTISOA1)+STT(I,J,L,IDTISOA2)+
     &         STT(I,J,L,IDTISOA3)+STT(I,J,L,IDTASOA1)+
     &         STT(I,J,L,IDTASOA2)+STT(I,J,L,IDTASOA3) )/VOL ! Total SV-&MV-SOA



        MSO4 = 0.d0
        DO N=1,NSO4
          SIZENUM=IDTSO4BIN1+N-1
          XM1D(N)=STT(I,J,L,SIZENUM)/VOL  !kg/m3
          MSO4 = MSO4 + XM1D(N)   ! total bin sulfate mass
          XN1D(N)=XN4D(I,J,L,N)
        ENDDO

        DO N=1,NSEA        ! sea salt
          SIZENUM=IDTSEABIN1+N-1
          XM1D(NSO4+N)=STT(I,J,L,SIZENUM)/VOL  !kg/m3
        ENDDO

        DO N=1,NDSTB      ! dust
          XMDST(N) = STT(I,J,L,(IDTDSTBIN1+N-1))/VOL   !kg/m3
          IF(XMDST(N).LT.1.d-20) THEN   !debug force
             XMDST(N) = 1.d-20
             STT(I,J,L,(IDTDSTBIN1+N-1)) = XMDST(N) * VOL
          ENDIF
        ENDDO

! BCOC associated with fossil fuel and biomass/biofuel seperated
        DO N= 1, NBCOCT
           SIZENUM=IDTBCPIFF+N-1
           MBCOC8(N)=STT(I,J,L,SIZENUM)/VOL  !kg/m3
        ENDDO

        MNIT = STT(I,J,L,IDTNIT)/VOL   !kg/m3
        MNH4 = STT(I,J,L,IDTNH4)/VOL   !kg/m3
        MMSA = STT(I,J,L,IDTMSA)/VOL   !kg/m3

! coated sulfate
        MBCS = STT(I,J,L,(IDTCTBCOC-1+1))/VOL ! SULF on BC kg/m3
        MOCS = STT(I,J,L,(IDTCTBCOC-1+2))/VOL ! SULF on OC kg/m3
        MDSTS = STT(I,J,L,IDTCTDST)/VOL ! SULF on DUST kg/m3
        MSALTS = STT(I,J,L,IDTCTSEA)/VOL ! SULF on SEA-SALT kg/m3

        NCOAG1 = NCOAG3D(I,J,L,1)
        NCOAG2 = NCOAG3D(I,J,L,2)

        CALL APM_PHYS(I,J,L,
     &       NCOAG1,NCOAG2,IACT10,IACT20,IACT30,NTEMPOUT1,
     &       PRESS,TK,RHIN,XQ,PLVSOG1,CACID,PACID,
     &       DTAPM,MMSA,MNIT,MNH4,MBCS, MOCS,
     &       MDSTS, MSALTS,MBCOC8, SOAT,
     &       CLVSOG,MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV,
     &       GFTOT1,GFTOT2,DENWET1,DENWET2,
     &       XM1D,XN1D,TEMPOUT1,XMDST,FCLOUD1)


        GFTOT3D(I,J,L,1) = GFTOT1
        GFTOT3D(I,J,L,2) = GFTOT2
        DENWET3D(I,J,L,1) = DENWET1
        DENWET3D(I,J,L,2) = DENWET2
        IACT1(I,J,L) = IACT10
        IACT2(I,J,L) = IACT20
        IACT3(I,J,L) = IACT30
        DO N=1,NSO4+4
           FCLOUD(I,J,L,N) = FCLOUD1(N)
        ENDDO
        NCOAG3D(I,J,L,1)=NCOAG1   ! count step that coag is not call in the grid
        NCOAG3D(I,J,L,2)=NCOAG2   ! count step that coag is not call in the grid

        STT(I,J,L,IDTSO4G)=  CACID*VOL*M1ACID*1.d6

        DO N=1,NSO4
          SIZENUM=IDTSO4BIN1+N-1
          STT(I,J,L,SIZENUM)=XM1D(N)*VOL
        ENDDO
        STT(I,J,L,IDTCTBCOC) = MBCS * VOL 
        STT(I,J,L,(IDTCTBCOC+1)) = MOCS * VOL
        STT(I,J,L,IDTCTDST) = MDSTS * VOL
        STT(I,J,L,IDTCTSEA) = MSALTS * VOL

        IF(IFAG.EQ.1) THEN
          STT(I,J,L,IDTSO4G+1)=  CLVSOG*VOL*M1LVSOG*1.d6
          STT(I,J,L,IDTCTSO4)= MSULFLV * VOL   !kg
          STT(I,J,L,IDTCTBCOC+2) = MBCLV * VOL   !kg
          STT(I,J,L,IDTCTBCOC+3) = MOCLV * VOL   !kg
          STT(I,J,L,IDTCTDST+1) = MDSTLV * VOL   !kg
          STT(I,J,L,IDTCTSEA+1) = MSALTLV * VOL   !kg
        ENDIF

        DO N=1,NSEA        ! sea salt
          SIZENUM=IDTSEABIN1+N-1
          STT(I,J,L,SIZENUM)=XM1D(NSO4+N)*VOL  !kg
        ENDDO

! no update for dust needed for now because of no coag yet (need update if coag)

      !luogan temp output
       IF(IFTEMPOUT)THEN
        DO N=1,NTEMPOUT1
         TEMPOUT(I,J,L,N)=TEMPOUT(I,J,L,N)+ TEMPOUT1(N)
        ENDDO
       ENDIF

      ENDDO
      ENDDO
      ENDDO

!$OMP END PARALLEL DO

      NPOUTSTEPS = NPOUTSTEPS + 1

      ! Clear the pointer
      NULLIFY( STT )

      END SUBROUTINE APM_DRIV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aeronum
!
! !DESCRIPTION: Subroutine AERONUM calculates aerosol number concentration
!  based on the mass.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AERONUM( State_Met, State_Chm )
!
! !USES:
!
      USE APM_INIT_MOD,       ONLY : NGCOND
      USE APM_INIT_MOD,       ONLY : NSO4
      USE APM_INIT_MOD,       ONLY : VDRY
      USE APM_INIT_MOD,       ONLY : DENSULF
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TRACER_MOD,         ONLY : N_TRACERS
                        
      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
! 
! !REVISION HISTORY: 
!  28 Aug 2008 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I,J,L,N,SIZENUM
      REAL*8  :: XMATEMP

      ! Make a pointer to the tracer array
      REAL*8, POINTER :: STT(:,:,:,:)
!
!******************************************************************************
      ! Point to STT
      STT => State_Chm%TRACERS

      XN4D = 0.D0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
!$OMP+PRIVATE( SIZENUM )
!$OMP+PRIVATE( XMATEMP )
!$OMP+SCHEDULE( DYNAMIC )

      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Should we put the N loop on the outside?
         DO N=1,NSO4
            SIZENUM =N_TRACERS+NGCOND+N
            XMATEMP=STT(I,J,L,SIZENUM) / State_Met%AIRVOL(I,J,L)
            XMATEMP= MAX( 1.d-40, XMATEMP)
            XN4D(I,J,L,N)=XMATEMP/(DENSULF*VDRY(N))*1.E-9 !XN4D in #/cm3, VDRY in m3
         ENDDO

      ENDDO
      ENDDO
      ENDDO

!$OMP END PARALLEL DO

      ! Clear the pointer
      NULLIFY( STT )

      END SUBROUTINE AERONUM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm3dqgeos
!
! !DESCRIPTION: Subroutine APM3DQGEOS ...
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE APM3DQGEOS( State_Met )
!
! !USES:
!
      USE DAO_MOD,            ONLY : IS_WATER
      USE DAO_MOD,            ONLY : IS_ICE
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_XMID
      USE GRID_MOD,           ONLY : GET_YMID
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE PRESSURE_MOD,       ONLY : GET_PEDGE
      USE APM_NUCL_MOD,       ONLY : IONRATE0

      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY: 
!  17 Mar 2010 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: ISURF, I, J, L
      REAL*8  :: XLON,XLAT,YPR,YQ,YPSURF

      DO J = 1, JJPAR
        DO I = 1, IIPAR


          YPSURF = GET_PEDGE(I,J,1)  !surface P in mb
          IF( IS_WATER( I, J, State_Met ) .or.
     &        IS_ICE  ( I, J, State_Met )) THEN
             ISURF = 0
          ELSE
             ISURF = 1
          ENDIF
          DO L = 1, LLPAR
             XLON = GET_XMID( I, J, L  )  ! Grid box longitude [degrees]
             XLAT = GET_YMID( I, J, L  )  ! Grid box latitude center [degree]
             YPR  = GET_PCENTER(I,J,L)    ! in mb
             CALL IONRATE0(ISURF, YPSURF, XLON, XLAT, YPR, YQ)
             XQ3D(I,J,L)= YQ
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE APM3DQGEOS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: Function GET\_OH returns OH from SMVGEAR's CSPEC array 
!  (for coupled runs) or monthly mean OH (for offline runs).  Imposes a 
!  diurnal variation on OH for offline simulations.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_OH( I, J, L, State_Met ) RESULT( OH_MOLEC_CM3 )
!
! !USES:
!
      USE COMODE_MOD,         ONLY : CSPEC, JLOP
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACER_MOD,         ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,         ONLY : ITS_AN_AEROSOL_SIM
      USE TRACERID_MOD,       ONLY : IDOH

      USE CMN_SIZE_MOD                          ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I           ! Longitude index
      INTEGER,        INTENT(IN)  :: J           ! Latitude index
      INTEGER,        INTENT(IN)  :: L           ! Level index
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
      REAL*8                      :: OH_MOLEC_CM3 ! OH [molec/cm3]
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version  
!  (1 ) We assume SETTRACE has been called to define IDOH (bmy, 11/1/02)
!  (2 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (3 ) Now reference ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM from
!        "tracer_mod.f".  Also replace IJSURF w/ an analytic function.
!        (bmy, 7/20/04)
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!  14 Dec 2012 - M. Payer    - Replace SUNCOS with State_Met%SUNCOS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: JLOOP

      !=================================================================
      ! GET_OH begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !---------------------
         ! Coupled simulation
         !---------------------

         ! JLOOP = SMVGEAR 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Take OH from the SMVGEAR array CSPEC
         ! OH is defined only in the troposphere
         IF ( JLOOP > 0 ) THEN
            OH_MOLEC_CM3 = CSPEC(JLOOP,IDOH)
         ELSE
            OH_MOLEC_CM3 = 0d0
         ENDIF

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

         !---------------------
         ! Offline simulation
         !---------------------

!         ! Test for sunlight...
!         IF ( State_Met%SUNCOS(I,J) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN
!
!            ! Impose a diurnal variation on OH during the day
!            OH_MOLEC_CM3 = OH(I,J,L)                      *
!     &                     ( State_Met%SUNCOS(I,J) / TCOSZ(I,J) ) *
!     &                     ( 1440d0        / GET_TS_CHEM() )
!
!            ! Make sure OH is not negative
!            OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0d0 )
!
!         ELSE

            ! At night, OH goes to zero
            OH_MOLEC_CM3 = 0d0

!         ENDIF

      ELSE

         !---------------------
         ! Invalid simulation
         !---------------------
         CALL ERROR_STOP( 'Invalid NSRCX!', 'GET_OH (sulfate_mod.f)')

      ENDIF

      END FUNCTION GET_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_apm3d
!
! !DESCRIPTION: Subroutine INIT_APM3D allocates and zeroes module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_APM3D
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE APM_INIT_MOD,ONLY : NSO4 , NSEA

      USE CMN_SIZE_MOD      ! Size parameters
! 
! !REVISION HISTORY: 
!  28 Aug 2008 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
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

      ALLOCATE( FCLOUD(IIPAR,JJPAR,LLPAR,NSO4+4), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FCLOUD' )
      FCLOUD = 0d0

      ALLOCATE( XN4D(IIPAR,JJPAR,LLPAR,NSO4), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'XN4D' )
      XN4D = 0d0

      ALLOCATE( XQ3D(IIPAR,JJPAR,LLPAR), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'XQ3D' )
      XQ3D = 0d0

      ALLOCATE( GFTOT3D(IIPAR,JJPAR,LLPAR,2), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFTOT3D' )
      GFTOT3D = 1d0

      ALLOCATE( DENWET3D(IIPAR,JJPAR,LLPAR,2), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DENWET3D' )
      DENWET3D = 2d0

      ALLOCATE( NCOAG3D(IIPAR,JJPAR,LLPAR,2), STAT=AS )
      NCOAG3D = 0d0

      ALLOCATE( IACT1(IIPAR,JJPAR,LLPAR), STAT=AS )
      IACT1 = 15

      ALLOCATE( IACT2(IIPAR,JJPAR,LLPAR), STAT=AS )
      IACT2 = 15

      ALLOCATE( IACT3(IIPAR,JJPAR,LLPAR), STAT=AS )
      IACT3 = 15

      ALLOCATE( PSO4GAS( IIPAR, JJPAR, LLPAR ), STAT=AS )
      PSO4GAS = 0.d0

      ALLOCATE( PLVSOG( IIPAR, JJPAR, LLPAR ), STAT=AS )
      PLVSOG = 0d0

      ALLOCATE( IFOUTIJ( IIPAR, JJPAR ) )
      IFOUTIJ = 0

      ALLOCATE( SITEID( IIPAR, JJPAR) )
      SITEID = 0

      END SUBROUTINE INIT_APM3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_apm3d
!
! !DESCRIPTION: Subroutine CLEANUP\_APM3D deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_APM3D
! 
! !REVISION HISTORY: 
!  28 Aug 2008 - F. Yu       - Initial version  
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_APM3D begins here!
      !=================================================================
      IF ( ALLOCATED( FCLOUD      ) ) DEALLOCATE( FCLOUD      )
      IF ( ALLOCATED( XN4D        ) ) DEALLOCATE( XN4D        )
      IF ( ALLOCATED( XQ3D        ) ) DEALLOCATE( XQ3D        )
      IF ( ALLOCATED( GFTOT3D     ) ) DEALLOCATE( GFTOT3D     )
      IF ( ALLOCATED( DENWET3D    ) ) DEALLOCATE( DENWET3D    )
      IF ( ALLOCATED( NCOAG3D     ) ) DEALLOCATE( NCOAG3D     )
      IF ( ALLOCATED( IACT1       ) ) DEALLOCATE( IACT1       )
      IF ( ALLOCATED( IACT2       ) ) DEALLOCATE( IACT2       )
      IF ( ALLOCATED( IACT3       ) ) DEALLOCATE( IACT3       )
      IF ( ALLOCATED( PSO4GAS     ) ) DEALLOCATE( PSO4GAS     )  
      IF ( ALLOCATED( PLVSOG      ) ) DEALLOCATE( PLVSOG      ) 
      IF ( ALLOCATED( TEMPOUT     ) ) DEALLOCATE( TEMPOUT     )
      IF ( ALLOCATED( IFOUTIJ     ) ) DEALLOCATE( IFOUTIJ     )
      IF ( ALLOCATED( SITEID      ) ) DEALLOCATE( SITEID      )

      END SUBROUTINE CLEANUP_APM3D
!EOC
      END MODULE APM_DRIV_MOD


