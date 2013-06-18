#if   defined( RRTMG )
        MODULE RRTMG_RAD_TRANSFER_MOD

        use parrrtm, only : nbndlw
        use parrrsw, only : nbndsw

        USE CMN_SIZE_MOD, ONLY : IIPAR,JJPAR,LLPAR, NDUST, NAER
        USE JV_CMN_MOD, ONLY: RTODAER, RTSSAER, RTASYMAER, 
     &                        WVAA, SPECMASK, LSPECRADMENU,
     &                        NSPECRADMENU
        USE DIAG_MOD,      ONLY : AD71 !rad output diagnostic array
!        IMPLICIT NONE
!        PRIVATE


        !PUBLIC MEMBER FUNCTIONS
        PUBLIC  :: CLEANUP_SURFACE_RAD
        PUBLIC  :: INIT_SURFACE_RAD
        PUBLIC  :: READ_SURFACE_RAD

!
! !PUBLIC DATA MEMBERS:
!

      REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: ALBDIR(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: ALBDIF(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: EMISS (:,:,:)

      REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CH4CLIM(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: N2OCLIM(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CFC11CLIM(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CFC12CLIM(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CCL4CLIM(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CFC22CLIM(:,:,:)

        REAL*8  :: LW_UFLUX(IIPAR*JJPAR,LLPAR+1)
        REAL*8  :: LW_DFLUX(IIPAR*JJPAR,LLPAR+1)
        REAL*8  :: SW_UFLUX(IIPAR*JJPAR,LLPAR+1)
        REAL*8  :: SW_DFLUX(IIPAR*JJPAR,LLPAR+1)
        REAL*8  :: LW_UFLUXC(IIPAR*JJPAR,LLPAR+1)
        REAL*8  :: LW_DFLUXC(IIPAR*JJPAR,LLPAR+1)
        REAL*8  :: SW_UFLUXC(IIPAR*JJPAR,LLPAR+1)
        REAL*8  :: SW_DFLUXC(IIPAR*JJPAR,LLPAR+1)

        real*8   :: rrtmg_lmb(nbndlw+nbndsw)

        INTEGER  :: id_aer_lmb0 (nbndlw+nbndsw)
        INTEGER  :: id_aer_lmb1 (nbndlw+nbndsw)
        CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_rrtmg_rad_transfer
!
! !DESCRIPTION: Internal subroutine 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_RRTMG_RAD_TRANSFER( THISDAY, THISMONTH, ICLD, 
     &                        ISPECMENU, ISEED)
!
! !USES:
!
      USE parkind,        ONLY : im=>kind_im, rb=>kind_rb
      USE rrtmg_lw_init,  ONLY : rrtmg_lw_ini
      USE rrtmg_sw_init,  ONLY : rrtmg_sw_ini

      
      !USE STRAT_CHEM_MOD, ONLY : GET_SPECIES
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP
      USE PRESSURE_MOD,   ONLY : GET_PCENTER,      GET_PEDGE
      USE COMODE_MOD,   ONLY : AIRDENS,        CSPEC,       JLOP
      USE COMODE_LOOP_MOD,   ONLY : ICH4, AVG, WTAIR
      USE DAO_MOD,      ONLY : BXHEIGHT, QI, QL, TAUCLI, TAUCLW
      USE DAO_MOD,      ONLY : CLDF, SUNCOS, AVGW, TSKIN, T, ALBD
      USE DAO_MOD,      ONLY : AIRVOL, AIRDEN
      USE TIME_MOD,     ONLY : GET_DAY_OF_YEAR, GET_HOUR
      USE TRACERID_MOD, ONLY : IDO3
      USE DIAG_MOD,     ONLY : AD71
      USE JV_CMN_MOD,   ONLY : NSPECRAD  ! number of species for RT
      USE JV_CMN_MOD,   ONLY : NASPECRAD ! number of aerosol species
      USE JV_CMN_MOD,   ONLY : SPECMASK, IWVSELECT, NWVSELECT
      USE JV_CMN_MOD,   ONLY : ACOEF_WV, BCOEF_WV, CCOEF_WV
      USE LOGICAL_MOD,  ONLY : LSUMRAD,LLWRAD,LSWRAD,LSKYRAD

      USE JV_CMN_MOD,   ONLY : DO3
      USE TOMS_MOD,     ONLY: GET_OVERHEAD_O3
      USE GRID_MOD,     ONLY : GET_YMID

 
!
      parameter (nwv=37)
      INTEGER, INTENT(IN) :: THISDAY      ! Current day
      INTEGER, INTENT(IN) :: THISMONTH    ! Current month
      INTEGER, INTENT(IN) :: ICLD         ! Cloud flag for RRTMG: 0-nocloud, 1-grey cloud
      INTEGER, INTENT(IN) :: ISPECMENU    ! The species being included
                                          ! Needed for output purposes
      INTEGER, INTENT(IN) :: ISEED  
      !LOGICAL, INTENT(IN) :: SPECMASK(NSPECRAD) ! Mask for number of species considered 
      LOGICAL             :: DOAERAD      ! Include aerosol determined from SPECMASK
      LOGICAL             :: LOUTPUTAERO  ! Output aerosol diagnostics?
      INTEGER             :: NAD71        ! Number of outputs per field
      INTEGER             :: ITIMEVALS(8)
      INTEGER             :: IDIAGOUT     ! Index of the species optics for output
      REAL*8              :: OLDSECS, NEWSECS

!      REAL*8,  INTENT(IN) :: ALBEDO(IIPAR,JJPAR)
! 

! !REMARKS:
!  This routine passes inputs to the RRTMG driver routine "RAD_DRIVER"
!  via the argument list.  This prevents circular references.
!
! !REVISION HISTORY: 
!  17 Aug 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

! Temporary aerosol variables
      real*8    :: lambda(nwv)


      ! Scalars
      INTEGER :: I, J ,L, LL, N, NCOL, JLOOP, IJLOOP
      INTEGER :: OUTIDX,IOUTWV
      integer :: ib,ib_sw,is,nbnds,nspec
      INTEGER :: IS_ON,NASPECRAD_ON
      INTEGER :: IASPECRAD_ON(NASPECRAD)
      INTEGER :: IRTWVSELECT(2) !required SW wavelengths for AOD output
      REAL*8  :: rhoice=0.9167, rholiq=1.    ! g/cm3
      ! rel and rei from personal communication from Lazaros Oreopoulos
      ! (GSFC) 12/12/12
      REAL*8  :: rel_def = 14.2, rei_def=24.8    ! microns
      REAL*8  :: DOY

      integer           :: ihr
      character(len=2)  :: chr 

      ! Arrays from GC
      REAL*8  :: PCENTER(IIPAR,JJPAR,LLPAR)
      REAL*8  :: PEDGE  (IIPAR,JJPAR,LLPAR)
      LOGICAL :: INTROP (IIPAR,JJPAR,LLPAR)

      !Arrays for RRTMG
      REAL(kind=rb) :: O3VMR (IIPAR*JJPAR,LLPAR)
      REAL(kind=rb) :: CH4VMR(IIPAR*JJPAR,LLPAR)
      REAL(kind=rb) :: N2OVMR(IIPAR*JJPAR,LLPAR)
      real(kind=rb) :: CFC11VMR(IIPAR*JJPAR,LLPAR)
      real(kind=rb) :: CFC12VMR(IIPAR*JJPAR,LLPAR)
      real(kind=rb) :: CFC22VMR(IIPAR*JJPAR,LLPAR)
      real(kind=rb) :: CCL4VMR(IIPAR*JJPAR,LLPAR)

      REAL(kind=rb) :: taucld(IIPAR*JJPAR,LLPAR)
      REAL(kind=rb) :: cldfr (IIPAR*JJPAR,LLPAR)
      REAL(kind=rb) :: reliq (IIPAR*JJPAR,LLPAR)
      REAL(kind=rb) :: reice (IIPAR*JJPAR,LLPAR)
      REAL(kind=rb) :: cliqwp (IIPAR*JJPAR,LLPAR)
      REAL(kind=rb) :: cicewp (IIPAR*JJPAR,LLPAR)
      REAL(kind=rb) :: albvis (IIPAR*JJPAR)

      REAL(kind=rb) :: tauaer_lw (IIPAR*JJPAR,LLPAR,nbndlw)
      REAL(kind=rb) :: tauaer_sw (IIPAR*JJPAR,LLPAR,nbndsw)
      REAL(kind=rb) :: ssaaer (IIPAR*JJPAR,LLPAR,nbndsw)
      REAL(kind=rb) :: asmaer (IIPAR*JJPAR,LLPAR,nbndsw)
      !TO STORE THE OPTICS FOR THE AEROSOLS WE ARE INTERESTED IN
      REAL*8        :: tauaerDIAG (IIPAR*JJPAR,LLPAR,nbndsw)
      REAL*8        :: ssaaerDIAG (IIPAR*JJPAR,LLPAR,nbndsw)
      REAL*8        :: asmaerDIAG (IIPAR*JJPAR,LLPAR,nbndsw)

      !Variables for obtaining stratospheric variables
      real*8       :: o3col, ylat, air_tmp

      !surface
      REAL(kind=rb) :: albdirvis (IIPAR*JJPAR)
      REAL(kind=rb) :: albdifvis (IIPAR*JJPAR)
      REAL(kind=rb) :: albdirnir (IIPAR*JJPAR)
      REAL(kind=rb) :: albdifnir (IIPAR*JJPAR)
      REAL(kind=rb) :: rtemiss (IIPAR*JJPAR,nbndlw)

      REAL*8        :: NUMER,DENOM
      REAL*4        :: AODTMP,AODOUT,SSATMP,SSAOUT
      REAL*4        :: ASYMTMP,ASYMOUT
!      INTEGER       :: flg_first_strat, flg_calc_strat_scl
      INTEGER       :: flg_first_strat(iipar*jjpar)
      REAL*4        :: ch4scl(iipar*jjpar)

! Local aerosol arrays and variables for interpolating to RRTMG bands ! and combining species
!      real*8    :: od_aer(ndust+naer+2)
!      real*8    :: ssa_aer(ndust+naer+2)
!      real*8    :: asym_aer(ndust+naer+2)
!      real*8    :: angst


      !=================================================================
      ! Prepare inputs for RAD_DRIVER
      !=================================================================

      !CALL GET_SPECIES( .TRUE., THISMONTH, 'O3',  GMI_O3  )
      !CALL GET_SPECIES( .TRUE., THISMONTH, 'CH4', GMI_CH4 )

      !=================================================================
      ! Prepare inputs for RRTMG
      !=================================================================

      ! Initialize
      !print *,'in rrtmg_rad'
      nspec = naer+ndust+4
      flg_first_strat(:) = 0 !flag to determine if the first stratospheric level has been reached
 
      !determine if we are running with aerosol
      !create index for aerosols required
      LOUTPUTAERO=.TRUE. !set aerosol diag output to true initially
      DOAERAD = .FALSE.
      NASPECRAD_ON = 0
      IASPECRAD_ON(:) = 0
      DO N=1,NASPECRAD
       IF (SPECMASK(N).GT.0) THEN
        DOAERAD = .TRUE.
        NASPECRAD_ON = NASPECRAD_ON +1
        IASPECRAD_ON(NASPECRAD_ON) = N
        IDIAGOUT = MAX(IDIAGOUT,SPECMASK(N))
        
       ENDIF
      ENDDO
!      write(6,*) 'CURRENT SPEC:',NASPECRAD_ON
!      write(6,*) IASPECRAD_ON
      !set up the SW wavelength required for AOD output
      !convert from LUT index to IB_SW
      DO N=1,NWVSELECT
       IRTWVSELECT(N)=IWVSELECT(N)-23
      ENDDO 
      taucld(:,:) = 0.0
      cldfr(:,:) = 0.0
      reliq(:,:) = 0.0
      reice(:,:) = 0.0
      cliqwp(:,:) = 0.0
      cicewp(:,:) = 0.0
      tauaer_lw(:,:,:) = 0.0
      tauaer_sw(:,:,:) = 0.0
      ssaaer(:,:,:) = 0.0
      asmaer(:,:,:) = 0.0
      TAUAERDIAG(:,:,:) = 0.0d0
      SSAAERDIAG(:,:,:) = 0.0d0
      ASMAERDIAG(:,:,:) = 0.0d0
      nbnds = nbndlw+nbndsw

!!X$OMP PARALLEL DO
!!X$OMP+DEFAULT( SHARED )
!!X$OMP+PRIVATE( I, J, IJLOOP, IB, IB_SW)
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ! 1-D index
         ijloop         = (j-1)*iipar+i
       !broadband albedo
       albdirvis(ijloop) = albdir(i,j,1)
       albdifvis(ijloop) = albdif(i,j,1)
       albdirnir(ijloop) = albdir(i,j,2)
       albdifnir(ijloop) = albdif(i,j,2)
       if ((albdirvis(ijloop).gt.0.999).or.
     &  (albdirvis(ijloop).lt.0.001)) then
        write(6,*) 'ALBEDO DRVIS out of range',i,j,
     &              albdirvis(ijloop)
       endif
       if ((albdifvis(ijloop).gt.0.999).or.
     &  (albdifvis(ijloop).lt.0.001)) then
        write(6,*) 'ALBEDO DFVIS out of range',i,j,
     &    albdifvis(ijloop)
       endif
       if ((albdirnir(ijloop).gt.0.999).or.
     &  (albdirnir(ijloop).lt.0.001)) then
        write(6,*) 'ALBEDO DRNIR out of range',i,j,
     &    albdirnir(ijloop)
       endif
       if ((albdifnir(ijloop).gt.0.999).or.
     &  (albdifnir(ijloop).lt.0.001)) then
        write(6,*) 'ALBEDO DFNIR out of range',i,j,
     &   albdifnir(ijloop)
       endif

      DO IB= 1, nbnds
         ib_sw = ib-nbndlw
         IF ( ib .le. 16) THEN
          rtemiss(ijloop,ib)     = emiss(i,j,ib)
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!!X$OMP END PARALLEL DO

      !I DONT THINK FUNCTION CALLS CAN BE PARALLELIZED
      !get PCENTER, PEDGE and determine if in TROP before parallel loop
      DO I = 1, IIPAR
      DO J = 1, JJPAR
      DO L = 1, LLPAR
         PCENTER(I,J,L) = GET_PCENTER( I, J, L )
         PEDGE  (I,J,L) = GET_PEDGE  ( I, J, L )
         INTROP (I,J,L) = ITS_IN_THE_TROP( I, J, L )
      ENDDO
!         o3col(I,J) = get_overhead_o3(I,J)
!         call set_prof_o3 (i,j,ylat,thismonth,thisday,t(i,j,:),o3col)
      ENDDO
      ENDDO
!!X$OMP PARALLEL DO
!!X$OMP+DEFAULT( SHARED )
!!X$OMP+PRIVATE( I, J, L, IJLOOP, JLOOP, AIR_TMP )
!!X$OMP+PRIVATE( YLAT, O3COL)
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
       
       ! 1-D index
       ijloop         = (j-1)*iipar+i
       !PCENTER(I,J,L) = GET_PCENTER( I, J, L )
       !PEDGE  (I,J,L) = GET_PEDGE  ( I, J, L )
       !INTROP (I,J,L) = ITS_IN_THE_TROP(I,J,L)      
       !if (.false.) then
       IF (INTROP(I,J,L).EQ..TRUE.) THEN 
            !-----------------------------
            ! We are in the troposphere
            !-----------------------------

            ! SET O3, CH4, N2O AND CFC PROFILES
            ! G-C chemistry is only done in the trop
            ! Therefore CSPEC will only be defined in the trop
            ! 1-D index for CSPEC array
            JLOOP = JLOP(I,J,L)
            !if ozone not selected then SPECMASK will be zero
            o3vmr(ijloop,l)  = (cspec(jloop,ido3)/airdens(jloop))*
     &                         SPECMASK(NASPECRAD+1) !gases come after aerosol
            ch4vmr(ijloop,l) = cspec(jloop,ich4)/airdens(jloop)*
     &                         SPECMASK(NASPECRAD+2)
            n2ovmr(ijloop,l) = n2oclim(i,j,l)/1e9

            !CFC climatology from UARS and MIPAS
            cfc11vmr(ijloop,l) = cfc11clim(i,j,l)/1e9
            cfc12vmr(ijloop,l) = cfc12clim(i,j,l)/1e9
            ccl4vmr(ijloop,l)  = ccl4clim(i,j,l)/1e9
            cfc22vmr(ijloop,l) = cfc22clim(i,j,l)/1e9


            taucld(ijloop,l) = taucli(i,j,l)+tauclw(i,j,l)
            cldfr(ijloop,l) =  cldf(i,j,l)
            albvis(ijloop)   = albd(i,j)
                         

            !-----------------------------
            ! Get cloud properties by setting reasonable values for Rel
            ! and Rei in microns and calculating LWP and IWP from
            ! visible optical depth  (in g/m2)
            ! This is a diferent approach than that used up through
            ! October, where we tried to used the cloud liquid and ice
            ! mixing ratios to get LWP and IWP. Note that rho in g/cm3
            ! and re in microns make the cxxxwp variables come out in
            ! g/m2.
            !-----------------------------

            !-----------------------------
            ! Set Re values to avoid problems in cloud code

            IF (ICLD.NE.0) THEN
             ! Liquid
             cliqwp(ijloop,l) = 0.667*tauclw(i,j,l)*rholiq*rel_def
             reliq(ijloop,l) = rel_def
             ! Ice
             cicewp(ijloop,l) = 0.667*taucli(i,j,l)*rhoice*rei_def
             reice(ijloop,l) = rei_def
            ENDIF !clouds
       ELSE
           !-----------------------------
           ! We are in the stratosphere
           !-----------------------------

           !! density of air in g/cm2
           air_tmp = airden(l,i,j)*bxheight(i,j,l)*1.0e-1
           !! density of air in molec/cm2
           air_tmp = avg*air_tmp/wtair
           o3col = get_overhead_o3(I,J)
           call set_prof_o3 (i,j,ylat,thismonth,thisday,t(i,j,:),o3col)
           o3vmr(ijloop,l) = do3(l)/air_tmp

            !get scalings if this is the first level in the strat
            if (flg_first_strat(ijloop).eq.0) then
                flg_first_strat(ijloop) = 1
                ch4scl(ijloop) = (cspec(jloop,ich4)/airdens(jloop))/
     &             (ch4clim(i,j,l)/1e9)
            endif
            
            !TES profiles interpolated to GC grid when saved
            !so we just need to scale to current conc at top of trop
              
            ch4vmr(ijloop,l) = ch4scl(ijloop)*ch4clim(i,j,l)/1e9
            n2ovmr(ijloop,l) = n2oclim(i,j,l)/1e9

            !CFC climatology from UARS and MIPAS
            cfc11vmr(ijloop,l) = cfc11clim(i,j,l)/1e9
            cfc12vmr(ijloop,l) = cfc12clim(i,j,l)/1e9
            ccl4vmr(ijloop,l)  = ccl4clim(i,j,l)/1e9
            cfc22vmr(ijloop,l) = cfc22clim(i,j,l)/1e9
       ENDIF
       !endif
      ENDDO
      ENDDO
      ENDDO
!!X$OMP END PARALLEL DO
      IF (DOAERAD) THEN
       DO IB = 1,NBNDS
        !rrtmg wavebands start at i=8 in GC arrays
        !in same order, so just offset
        IBX=IB+7
        IB_SW = IB-NBNDLW
        DO IS = 1,NASPECRAD
         !THE AEROSOL SPECIES WE ARE CURRENTLY CALCULATING FOR WILL BE
         !SET TO THE LSPECRADMENU VALUE FOR THAT SPECIES.
         !THIS MEANS THAT RRTMG REQUIRES *ALL OTHER* SPECIES SO THAT THE
         !FLUX IN ABSENCE OF THE SPECIES CAN BE CALCULATED
         !THEREFORE WE COMPILE TWO SETS OF AEROSOL PROPERTIES:
         !(1) ALL BUT THE CURRENT SPECIES TO SEND TO RRTMG
         !(2) THE CURRENT SPECIES FOR OUTPUT TO THE RT DIAGNOSTICS
         ! We must merge aerosol properties for the species to be output
         ! (i.e. combine hydrophilic/phobic and multiple sizes for the same) 
!!X$OMP PARALLEL DO
!!X$OMP+DEFAULT( SHARED )
!!X$OMP+PRIVATE( I, J, L, IJLOOP)
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
          ! 1-D index
          IJLOOP = (J-1)*IIPAR+I       
          IF (INTROP(I,J,L).EQ..TRUE.) THEN 
           IF (IB.LE.16) THEN !LW
            IF (SPECMASK(IS).EQ.1) THEN
                  tauaer_lw(ijloop,l,ib) = tauaer_lw(ijloop,l,ib) +
     &                           rtodaer(i,j,l,ibx,is)
            ENDIF
           ELSE !SW
            !IF SPECMASK(IS)=1 THEN WE AGGREGATE THAT SPECIES FOR RRTMG
            !IF SPECMASK(IS)>1 THEN WE SAVE THAT SPECIES FOR DIAG OUTPUT
            IF (SPECMASK(IS).EQ.1) THEN
             tauaer_sw(ijloop,l,ib_sw)=tauaer_sw(ijloop,l,ib_sw)+
     &                       rtodaer(i,j,l,ibx,is)
             ssaaer(ijloop,l,ib_sw) =  ssaaer(ijloop,l,ib_sw) + 
     &                       rtssaer(i,j,l,ibx,is)*rtodaer(i,j,l,ibx,is)
             asmaer(ijloop,l,ib_sw) = asmaer(ijloop,l,ib_sw) +
     &                     rtasymaer(i,j,l,ibx,is) *
     &                     rtodaer(i,j,l,ibx,is)*rtssaer(i,j,l,ibx,is)
            ENDIF 
            IF (SPECMASK(IS).GT.1) THEN
             TAUAERDIAG(ijloop,l,ib_sw)=TAUAERDIAG(ijloop,l,ib_sw)+
     &                       rtodaer(i,j,l,ibx,is)
             SSAAERDIAG(ijloop,l,ib_sw) = SSAAERDIAG(ijloop,l,ib_sw) +
     &                       rtssaer(i,j,l,ibx,is)*rtodaer(i,j,l,ibx,is)
             ASMAERDIAG(ijloop,l,ib_sw) = ASMAERDIAG(ijloop,l,ib_sw) +
     &                     rtasymaer(i,j,l,ibx,is) *
     &                     rtodaer(i,j,l,ibx,is)*rtssaer(i,j,l,ibx,is)
              !TAUAERDIAG(ijloop,l,ib_sw) = 1.0d0
              !ASMAERDIAG(IJLOOP,L,IB_SW) = 0.02d0
              !SSAAERDIAG(IJLOOP,L,IB_SW) = 0.99d0*0.02d0
            ENDIF
           ENDIF 
          ENDIF
         ENDDO
         ENDDO
         ENDDO
!!X$OMP END PARALLEL DO 
        ENDDO !SPECIES
        !NOW AEROSOL HAVE BEEN SUMMED AND WEIGHTED BY AOD AND SSA
        !DIVIDE THROUGH BY TOTAL AOD (FOR SSA) AND AOD*SSA (FOR ASYM)
        IF (IB.GT.16) THEN !SW

!switching this OMP on halves the optics values for some unknown reason.
!The benefit of having it switched on is ~1% reduction in processing
!time so currently doing without parallelizing this. DAR 02/2013
!!!!!$OMP PARALLEL DO
!!!!!$OMP+DEFAULT( SHARED )
!!!!!$OMP+PRIVATE( I, J, L, IJLOOP )

         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         IJLOOP = (J-1)*IIPAR+I       
          IF (INTROP(I,J,L).EQ..TRUE.) THEN 
           IF ((tauaer_sw(ijloop,l,ib_sw).GT.0).AND.
     $        (    ssaaer(ijloop,l,ib_sw).GT.0)) THEN
              !divide sum(asym*ssa*od) by sum(ssa*od) to get 
              !od*ssa weighted asym
              asmaer(ijloop,l,ib_sw) = asmaer(ijloop,l,ib_sw) /
     &                                 ssaaer(ijloop,l,ib_sw)
              !divide sum(ssa*od) by sum(od) to get od weighted ssa
              SSAAER(ijloop,l,ib_sw) = ssaaer(ijloop,l,ib_sw) / 
     &                                 tauaer_sw(ijloop,l,ib_sw)
           ENDIF
                !AND DO THE SAME FOR THE SPECIES WE'RE INTERESTED IN
           IF ((tauaerDIAG(ijloop,l,ib_sw).GT.0).AND.
     $        ( ssaaerDIAG(ijloop,l,ib_sw).GT.0)) THEN
              !divide sum(asym*ssa*od) by sum(ssa*od) to get
              !od*ssa weighted asym
              asmaerDIAG(ijloop,l,ib_sw) = asmaerDIAG(ijloop,l,ib_sw) /
     &                                     ssaaerDIAG(ijloop,l,ib_sw)
              !divide sum(ssa*od) by sum(od) to get od weighted ssa
              ssaaerDIAG(ijloop,l,ib_sw) = ssaaerDIAG(ijloop,l,ib_sw) /
     &                                     tauaerDIAG(ijloop,l,ib_sw)
           ENDIF
          ENDIF
         ENDDO
         ENDDO
         ENDDO
!!!!!$OMP END PARALLEL DO
        ENDIF
       ENDDO !BAND
      ELSE
       !no aerosol, set all to safe values
!!X$OMP PARALLEL DO
!!X$OMP+DEFAULT( SHARED )
!!X$OMP+PRIVATE( I, J, L, IB, IJLOOP, JLOOP, IB_SW)
       DO IB= 1, NBNDS
        IB_SW = IB-NBNDLW
        DO L = 1, LLPAR
        DO J = 1, JJPAR
        DO I = 1, IIPAR
         ! 1-D index
         IJLOOP = (J-1)*IIPAR+I       
         IF (INTROP(I,J,L).EQ..TRUE.) THEN 
          IF (IB.LE.16) THEN
           TAUAER_LW(IJLOOP,L,IB)    = 0.0
          ELSE
           TAUAER_SW(IJLOOP,L,IB_SW) = 0.0d0
           SSAAER(IJLOOP,L,IB_SW)    = 0.99d0
           ASMAER(IJLOOP,L,IB_SW)    = 0.02d0
           TAUAERDIAG(IJLOOP,L,IB_SW) = 0.0d0
           SSAAERDIAG(IJLOOP,L,IB_SW) = 0.99d0
           ASMAERDIAG(IJLOOP,L,IB_SW) = 0.02d0
          ENDIF
         ENDIF
        ENDDO
        ENDDO
        ENDDO
       ENDDO
!!X$OMP END PARALLEL DO 
      ENDIF
      ! Initialization
      !IF (LLWRAD) THEN
       CALL RRTMG_LW_INI
      !ENDIF
      IF (LSWRAD) THEN
       CALL RRTMG_SW_INI
      ENDIF

      ! RAD_DRIVER is the interface between G-C and RRTMG
      doy = get_day_of_year()
      ncol = iipar*jjpar
      !print *,'will call rad_driver'
      CALL RAD_DRIVER( IIPAR, JJPAR,   LLPAR, NCOL, ICLD,ISEED,T,TSKIN,
     &                 AVGW,  SUNCOS,  DOY,   PEDGE, PCENTER,
     &                 O3VMR, CH4VMR, N2OVMR,
     &                 cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr,
     &                 TAUCLD,CLDFR,CLIQWP,CICEWP,
     &                 RELIQ,REICE,
     &                 ALBDIRVIS,ALBDIFVIS,ALBDIRNIR,ALBDIFNIR,
     &                 RTEMISS,LLWRAD,LSWRAD,
     &                 TAUAER_LW,TAUAER_SW,SSAAER,ASMAER,
     &                 LW_UFLUX,LW_DFLUX,SW_UFLUX,SW_DFLUX,  
     &                 LW_UFLUXC,LW_DFLUXC,SW_UFLUXC,SW_DFLUXC)

      ! Output radiation variables to diagnostic
      ! if calc with aerosols and gases completed
      ! USE ISPECMENU to determine which flux has been calculated
      ! if =0 then we just did the baseline flux
      ! else if LSUMRAD is switched on then this is the combined
      IF (ISPECMENU.EQ.0) THEN
       OUTIDX = 1
      ELSE IF ( LSUMRAD ) THEN
        ! store in combined (2nd set)
        OUTIDX = 2
      ELSE
        ! DIAGNOSTIC INDEX IS ISPECMENU+2 (BASELINE=1 SUM=2)
        OUTIDX = ISPECMENU + 2
        !OUTIDX = IDIAGOUT + 2
      ENDIF

      !the number of ND71 outputs per field
      NAD71 = NSPECRADMENU + 2

      !first check if we have already output aerosol diagnostics
      !(i.e. if both ALL-SKY and CLEAR-SKY are switched on)
      IF ((LSKYRAD(1)).AND.(LSKYRAD(2))) THEN
      !we only need to output during one of these, so dont when ICLD=0
       IF (ICLD.EQ.0) THEN
        LOUTPUTAERO=.FALSE.
       ENDIF
      ENDIF


!!X$OMP PARALLEL DO
!!X$OMP+DEFAULT( SHARED )
!!X$OMP+PRIVATE( I, J, LL, IJLOOP )
!!X$OMP+PRIVATE( AODTMP, SSATMP, ASYMTMP)
!!X$OMP+PRIVATE( AODOUT, SSAOUT, ASYMOUT)
      DO J=1,JJPAR
      DO I=1,IIPAR
       IJLOOP   = (J-1)*IIPAR+I
       IF (ICLD.GT.0) THEN
        !all-sky (we get clear-sky with this too)
        !N.B. Upwelling should be negative as down is +ve
        AD71(I,J,OUTIDX) = AD71(I,J,OUTIDX) -  
     &                     SNGL(SW_UFLUX(IJLOOP,LLPAR+1))
        AD71(I,J,OUTIDX+NAD71) = AD71(I,J,OUTIDX+NAD71) +
     &                     SNGL(SW_DFLUX(IJLOOP,1))
        AD71(I,J,OUTIDX+2*NAD71) = AD71(I,J,OUTIDX+2*NAD71) -
     &                     SNGL(LW_UFLUX(IJLOOP,LLPAR+1))
        AD71(I,J,OUTIDX+3*NAD71) = AD71(I,J,OUTIDX+3*NAD71) +
     &                     SNGL(LW_DFLUX(IJLOOP,1))
        AD71(I,J,OUTIDX+4*NAD71) = AD71(I,J,OUTIDX+4*NAD71) -
     &                     SNGL(SW_UFLUXC(IJLOOP,LLPAR+1))
        AD71(I,J,OUTIDX+5*NAD71) = AD71(I,J,OUTIDX+5*NAD71) +
     &                     SNGL(SW_DFLUXC(IJLOOP,1))
        AD71(I,J,OUTIDX+6*NAD71) = AD71(I,J,OUTIDX+6*NAD71) -
     &                     SNGL(LW_UFLUXC(IJLOOP,LLPAR+1))
        AD71(I,J,OUTIDX+7*NAD71) = AD71(I,J,OUTIDX+7*NAD71) +
     &                     SNGL(LW_DFLUXC(IJLOOP,1))
       ELSE
        !clear-sky (running with clouds off)
        AD71(I,J,OUTIDX+4*NAD71) = AD71(I,J,OUTIDX+4*NAD71) -
     &                     SNGL(SW_UFLUX(IJLOOP,LLPAR+1))
        AD71(I,J,OUTIDX+5*NAD71) = AD71(I,J,OUTIDX+5*NAD71) +
     &                     SNGL(SW_DFLUX(IJLOOP,1))
        AD71(I,J,OUTIDX+6*NAD71) = AD71(I,J,OUTIDX+6*NAD71) -
     &                     SNGL(LW_UFLUX(IJLOOP,LLPAR+1))
        AD71(I,J,OUTIDX+7*NAD71) = AD71(I,J,OUTIDX+7*NAD71) +
     &                     SNGL(LW_DFLUX(IJLOOP,1))
       ENDIF
       !output optics for each
       !check that we have some aerosol to output 
       !IF OUTIDX=1 THEN IT IS THE BASELINE WHICH NOW CONTAINS ALL
       !AEROSOL - IDEALLY THIS NEEDS DUMPING INTO THE PM TRACER
       !SKIP OUTIDX=3,4 (OZONE AND CH4)
       IF (((OUTIDX.LE.2).OR.(OUTIDX.GE.5)).AND.(LOUTPUTAERO)) THEN
       !interpolate to the requested wavelength       
      ! IF (NWVSELECT.GT.1) THEN
        AODTMP  = 0.0d0
        SSATMP  = 0.0d0
        ASYMTMP = 0.0d0
        AODOUT  = 0.0d0
        SSAOUT  = 0.0d0
        ASYMOUT = 0.0d0
        !! Test impact on RTAOD of only going to LEV=30
        DO LL=1,30 !!LLPAR
         !check AOD is non-zero before LOG...
         IF((TAUAERDIAG(IJLOOP,LL,IRTWVSELECT(2)).GT.0).AND.
     &      (TAUAERDIAG(IJLOOP,LL,IRTWVSELECT(1)).GT.0)) THEN
          AODTMP=SNGL(TAUAERDIAG(IJLOOP,LL,IRTWVSELECT(2))*ACOEF_WV**
     &           (BCOEF_WV*LOG(TAUAERDIAG(IJLOOP,LL,IRTWVSELECT(1))/
     &            TAUAERDIAG(IJLOOP,LL,IRTWVSELECT(2)))))
          SSATMP=SNGL( CCOEF_WV*SSAAERDIAG(IJLOOP,LL,IRTWVSELECT(2))+
     &     (1.0d0-CCOEF_WV)*SSAAERDIAG(IJLOOP,LL,IRTWVSELECT(1)))*AODTMP
          
          ASYMTMP=SNGL( CCOEF_WV*ASMAERDIAG(IJLOOP,LL,IRTWVSELECT(2))+
     &     (1.0d0-CCOEF_WV)*ASMAERDIAG(IJLOOP,LL,IRTWVSELECT(1)))*SSATMP
          AODOUT=AODOUT+AODTMP
          SSAOUT=SSAOUT+SSATMP
          ASYMOUT=ASYMOUT+ASYMTMP
         ENDIF
        ENDDO
        !WE ARE SAVING COLUMN AVERAGED VALUES FOR EACH SPECIES
        !divide through by number of gridboxes and aod*ssa
        ASYMOUT=ASYMOUT/SSAOUT
        !divide through by number of gridboxes and aod
        SSAOUT=SSAOUT/AODOUT
!     $              OUTIDX,AODOUT,SSAOUT,ASYMOUT
        AD71(I,J,OUTIDX+8*NAD71) = AD71(I,J,OUTIDX+8*NAD71) +
     &                             AODOUT
        AD71(I,J,OUTIDX+9*NAD71) = AD71(I,J,OUTIDX+9*NAD71) +
     &                             SSAOUT
        AD71(I,J,OUTIDX+10*NAD71)= AD71(I,J,OUTIDX+10*NAD71) +
     &                             ASYMOUT
        
       ENDIF
      ENDDO 
      ENDDO
!!X$OMP END PARALLEL DO

!EOC
      END SUBROUTINE DO_RRTMG_RAD_TRANSFER

!******************************************************************************

      SUBROUTINE SET_SPECMASK(ISPECRADMENU)
!     This routine converts the species switches in the input.mod
!     radiation section into the list of species that should be passed
!     through to RRTMG. This must be done in a subtractive way, e.g. if
!     we require the DRE of sulfate then the baseline will contain all
!     species and the sulfate run will contain everything but sulfate,
!     this way the contribution of sulfate can be inferred. Therefore,
!     all species are initially set to 1 and their inclusion results in
!     SPECMASK for the particular species being set to zero.

      USE JV_CMN_MOD,   ONLY :  SPECMASK,NASPECRAD,
     &                          LSPECRADMENU,NSPECRADMENU
      USE LOGICAL_MOD, ONLY  :  LSUMRAD
      INTEGER, INTENT(IN)    :: ISPECRADMENU
      INTEGER                :: N0,N

      !convert the current species selection from the input menu into
      !the required species to be included in the RRTMG calculation
      SPECMASK(:)=1
      !If ISPECRADMENU is zero, we just wanted SPECMASK zeroing, so skip
      IF (ISPECRADMENU.GT.0) THEN
      !need to loop over all if sum requested
      IF (LSUMRAD) THEN
       N0 = 1
       N  = NSPECRADMENU
      ELSE
       !otherwise just set single species
       N0 = ISPECRADMENU
       N  = ISPECRADMENU
      ENDIF
      DO I = N0, N
       ! for the benefit of LSUMRAD=T, just combine switched on species
       IF(LSPECRADMENU(I).EQ.1) THEN 
       SELECT CASE( I )
       !ozone
       CASE( 1 ) 
        SPECMASK(17)=0
        !write(6,*) 'RT:OZONE'
       !methane
       CASE( 2 )
        SPECMASK(18)=0
        !write(6,*) 'RT:CH4'
       !SU
       CASE( 3 )
        SPECMASK(1)=3
        !write(6,*) 'RT:SU'
       !NI
       CASE( 4 )
        SPECMASK(2)=4
        !write(6,*) 'RT:NI'
       !AM
       CASE( 5 )
        SPECMASK(3)=5
        !write(6,*) 'RT:AM'
       !BC
       CASE( 6 )
        SPECMASK(4)=6
        SPECMASK(8)=6
        !write(6,*) 'RT:BC'
       !OA
       CASE( 7 )
        SPECMASK(5)=7
        SPECMASK(9)=7
        !write(6,*) 'RT:OA'
       !SS
       CASE( 8 )
        SPECMASK(6)=8
        SPECMASK(7)=8
        !write(6,*) 'RT:SS'
       !Dust
       CASE( 9 )
        SPECMASK(10)=9
        SPECMASK(11)=9
        SPECMASK(12)=9
        SPECMASK(13)=9
        SPECMASK(14)=9
        SPECMASK(15)=9
        SPECMASK(16)=9
        !write(6,*) 'RT:DU'
       !PM - all aerosol
       CASE( 10 )
        !set all but gases to 0
        DO II = 1, NASPECRAD
         SPECMASK(II)=10
        ENDDO
       ! write(6,*) 'RT:PM'
       END SELECT
       ENDIF
      ENDDO      
      ENDIF

      END SUBROUTINE SET_SPECMASK

!******************************************************************************

      SUBROUTINE READ_SURFACE_RAD(FORCEREAD)

!
!******************************************************************************
!  Subroutine READ_SURFACE_RAD gets the surface albedo and emissivity
!  from data files processed from MODIS MCD43C3.5 and MOD11C2
!  Albedo is direct (black sky) and diffuse (white sky) and interpolated
!  to the first 21 wavelengths of GADS, used in the RRTMG code
!  emissivity has been interpolated to last 40 wavelengths of GADS
!
!  This routine is called from main.F when ITS_TIME_FOR_SURFACE_RAD() is
!  true (every 8 days)
!  DAR (10/2012)
!******************************************************************************
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TIME_MOD,      ONLY : GET_YEAR, GET_MONTH 
      USE TIME_MOD,      ONLY : GET_DAY,  GET_DAY_OF_YEAR  
      USE TIME_MOD,      ONLY : YMD_EXTRACT
      USE JULDAY_MOD,    ONLY : JULDAY, CALDATE

      USE CMN_SIZE_MOD        ! Size parameters
      USE JV_CMN_MOD

      LOGICAL, INTENT(IN), OPTIONAL :: FORCEREAD   ! Reset first-time flag?

      ! Local variables
      LOGICAL               :: DOREAD
      INTEGER               :: THISDAY, THISMONTH, THISYEAR
      INTEGER               :: THISDOFY, I, J
      INTEGER, SAVE         :: LASTDAY = -1
      REAL*4                :: ARRAY1(IIPAR,JJPAR)
      REAL*4                :: ARRAY2(IIPAR,JJPAR,NEMISS)
      REAL*8                :: XTAU
      CHARACTER(LEN=255)    :: FILENAME
      CHARACTER(LEN=4)      :: STHISYEAR
      CHARACTER(LEN=3)      :: STHISDOFY
      INTEGER               :: NYMD, NHMS
      REAL*8                :: JDNEW


      ! READ_SURFACE_RAD begins here!
      !=================================================================
      !Flag to force read from closest day
      DOREAD=.FALSE.
      ! Reset the FIRST flag if
      IF ( PRESENT( FORCEREAD ) ) THEN
         IF ( FORCEREAD ) DOREAD = .TRUE.
      ENDIF


      ! TAU for reading the bpch files
      ! Fix year as 2002 for climatology
      THISYEAR  = 2002 !GET_YEAR()
      THISMONTH = GET_MONTH()
      THISDOFY  = GET_DAY_OF_YEAR()
      THISDAY   = GET_DAY() 
      
      !Force read of closest file
      IF ( DOREAD ) THEN
       !Number of days until next file
       DAYSAWAY = 9 - MOD(THISDOFY,8)
       !Read closest 8-day file
       IF ( DAYSAWAY .GE. 4) THEN
        THISDOFY = THISDOFY - (8 - DAYSAWAY)
       ELSE
        THISDOFY = THISDOFY + DAYSAWAY
       ENDIF
      ENDIF

       ! Calculate new month and day from Julian Day (clh)
      JDNEW = THISDOFY + JULDAY(THISYEAR, 1, 0d0)
      CALL CALDATE(JDNEW, NYMD, NHMS)
      CALL YMD_EXTRACT( NYMD, THISYEAR, THISMONTH, THISDAY )

      IF ( THISDOFY == LASTDAY ) THEN
       WRITE(6,*) 'SURFACE RAD ALREADY READ', 
     &            THISDAY,THISMONTH,THISYEAR
       RETURN
      ENDIF
      WRITE(STHISDOFY,101) THISDOFY
101   FORMAT( I3.3 )
      WRITE(STHISYEAR,102) THISYEAR
102   FORMAT( I4.4 )
      XTAU      = GET_TAU0( THISMONTH, THISDAY, THISYEAR )


      ! Filename and time
      FILENAME  = TRIM( DATA_DIR )         //
     &            'modis_surf_201210/' //
     &            'modis_surf_02-07_' //
     &            STHISDOFY//'.'// GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - GET_SURFACE_RAD_DATA: Reading ', a )
      !WRITE( 6,*) 'XTAU: ',XTAU
      !-----------------------
      ! Direct surface albedo (VIS)
      !-----------------------
      CALL READ_BPCH2( FILENAME, 'MODAVG-$', 71,
     &                 XTAU,      IIPAR,    JJPAR,
     &                 1,     ARRAY1,    QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      DO I=1,IIPAR
      DO J=1,JJPAR
      ALBDIR(I,J,1) = ARRAY1(I,J)
      ENDDO
      ENDDO

      !write(6,*) 'ALBDIR1 read'
      !-----------------------
      ! Diffuse surface albedo (VIS)
      !-----------------------
      CALL READ_BPCH2( FILENAME, 'MODAVG-$', 73,
     &                 XTAU,      IIPAR,    JJPAR,
     &                 1,     ARRAY1,    QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      DO I=1,IIPAR
      DO J=1,JJPAR
      ALBDIF(I,J,1) = ARRAY1(I,J)
      ENDDO
      ENDDO

      !write(6,*) 'ALBDIF1 read'
      !-----------------------
      ! Direct surface albedo (NIR)
      !-----------------------
      CALL READ_BPCH2( FILENAME, 'MODAVG-$', 72,
     &                 XTAU,      IIPAR,    JJPAR,
     &                 1,     ARRAY1,    QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      DO I=1,IIPAR
      DO J=1,JJPAR
      ALBDIR(I,J,2) = ARRAY1(I,J)
      ENDDO
      ENDDO

      !write(6,*) 'ALBDIR2 read'
      !-----------------------
      ! Diffuse surface albedo (NIR)
      !-----------------------
      CALL READ_BPCH2( FILENAME, 'MODAVG-$', 74,
     &                 XTAU,      IIPAR,    JJPAR,
     &                 1,     ARRAY1,    QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      DO I=1,IIPAR
      DO J=1,JJPAR
      ALBDIF(I,J,2) = ARRAY1(I,J)
      ENDDO
      ENDDO

      !write(6,*) 'ALBDIF2 read'
      !-----------------------
      ! Surface emissivity
      !-----------------------
      CALL READ_BPCH2( FILENAME, 'MODAVG-$', 75,
     &                 XTAU,      IIPAR,    JJPAR,
     &                 NEMISS,     ARRAY2,    QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      EMISS(:,:,:) = ARRAY2(:,:,:)

      !save this day to prevent re-read of same file
      LASTDAY = THISDOFY
      WRITE(6,*) 'SURFACE REFLECTIVITY READ SUCCESSFULLY'
      ! Return to calling program
      END SUBROUTINE

!******************************************************************************

      SUBROUTINE READ_STRAT_CLIM()

!
!******************************************************************************
!  Subroutine READ_STRAT_CLIM gets the stored N2O and CH4 profiles
!  created based on TES climatology.
!  These are global July mean profiles from TES that are 
!  interpolated to the 3D GEOS-chem grid as a simple solution to
!  account for stratospheric contribution of these species to radiative
!  balance.
!
!  This routine is called from main.F one at the start of the model run
!  DAR (12/2012)
!******************************************************************************
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TIME_MOD,      ONLY : GET_YEAR, GET_MONTH
      USE TIME_MOD,      ONLY : GET_DAY,  GET_DAY_OF_YEAR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D

      USE CMN_SIZE_MOD        ! Size parameters
      USE JV_CMN_MOD


      ! Local variables
      LOGICAL               :: DOREAD
      REAL*4                :: ARRAY2(IIPAR,JJPAR,LGLOB)
      REAL*8                :: XTAU
      CHARACTER(LEN=255)    :: FILENAME


      ! READ_STRAT_CHEM begins here!
      !=================================================================

      ! Filename and time
      FILENAME  = TRIM( DATA_DIR )         //
!     &            'strat_chem_201206/' //
     &            'species_clim_profiles' //
     &            '.'// GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - GET_TES_STRAT_CLIM: Reading ', a )

      XTAU      = GET_TAU0( 1, 1, 2000 )
      WRITE(6,*) XTAU,LLPAR
      !-----------------------
      ! N2O
      !-----------------------
      CALL READ_BPCH2( FILENAME, 'O3CLIM-$', 71,
     &                 XTAU,      IIPAR,    JJPAR,
     &                 LGLOB,     ARRAY2,    QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_3D( ARRAY2, N2OCLIM )

      !-----------------------
      ! CH4
      !-----------------------
      CALL READ_BPCH2( FILENAME, 'O3CLIM-$', 72,
     &                 XTAU,      IIPAR,    JJPAR,
     &                 LGLOB,     ARRAY2,    QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_3D( ARRAY2, CH4CLIM)

      WRITE(6,*) 'CFC MIPAS & UARS CLIMATOLOGY READ SUCCESSFULLY'

      !-----------------------
      ! CFC11
      !-----------------------
      CALL READ_BPCH2( FILENAME, 'O3CLIM-$', 73,
     &                 XTAU,      IIPAR,    JJPAR,
     &                 LGLOB,     ARRAY2,    QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_3D( ARRAY2, CFC11CLIM )

      !-----------------------
      ! CFC12
      !-----------------------
      CALL READ_BPCH2( FILENAME, 'O3CLIM-$', 74,
     &                 XTAU,      IIPAR,    JJPAR,
     &                 LGLOB,     ARRAY2,    QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_3D( ARRAY2, CFC12CLIM)

      !-----------------------
      ! CCL4
      !-----------------------
      CALL READ_BPCH2( FILENAME, 'O3CLIM-$', 75,
     &                 XTAU,      IIPAR,    JJPAR,
     &                 LGLOB,     ARRAY2,    QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_3D( ARRAY2, CCL4CLIM )

      !-----------------------
      ! CFC22
      !-----------------------
      CALL READ_BPCH2( FILENAME, 'O3CLIM-$', 76,
     &                 XTAU,      IIPAR,    JJPAR,
     &                 LGLOB,     ARRAY2,    QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_3D( ARRAY2, CFC22CLIM)

     

      ! Return to calling program
      END SUBROUTINE

!******************************************************************************


!------------------------------------------------------------------------------

      SUBROUTINE INIT_SURFACE_RAD
!
!******************************************************************************
!  Subroutine INIT_SURFACE_RAD initializes all allocatable module arrays.
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE JV_CMN_MOD
      USE CMN_SIZE_MOD

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ALLOCATE( ALBDIR( IIPAR, JJPAR, NALBD ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALBDIR' )
      ALBDIR = 0d0

      ALLOCATE( ALBDIF( IIPAR, JJPAR, NALBD ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALBDIF' )
      ALBDIF = 0d0

      ALLOCATE( EMISS( IIPAR, JJPAR, NEMISS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMISS' )
      EMISS = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_SURFACE_RAD

!------------------------------------------------------------------------------

      SUBROUTINE INIT_STRAT_CLIM
!
!******************************************************************************
!  Subroutine INIT_STRAT_CLIM initializes all allocatable module arrays.
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE JV_CMN_MOD
      USE CMN_SIZE_MOD

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ALLOCATE( N2OCLIM( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'N2OCLIM' )
      N2OCLIM = 0d0

      ALLOCATE( CH4CLIM( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4CLIM' )
      CH4CLIM = 0d0

      ALLOCATE( CFC11CLIM( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CFC11CLIM' )
      CFC11CLIM = 0d0

      ALLOCATE( CFC12CLIM( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CFC12CLIM' )
      CFC12CLIM = 0d0

      ALLOCATE( CCL4CLIM( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CCL4CLIM' )
      CCL4CLIM = 0d0

      ALLOCATE( CFC22CLIM( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CFC22CLIM' )
      CFC22CLIM = 0d0


      ! Return to calling program
      END SUBROUTINE INIT_STRAT_CLIM

!------------------------------------------------------------------------------


!
      SUBROUTINE CLEANUP_SURFACE_RAD
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! CLEANUP_DAO begins here!
      !=================================================================
      IF ( ALLOCATED( ALBDIR          ) ) DEALLOCATE( ALBDIR          )
      IF ( ALLOCATED( ALBDIF          ) ) DEALLOCATE( ALBDIF          )
      IF ( ALLOCATED( EMISS           ) ) DEALLOCATE( EMISS           )

      END SUBROUTINE CLEANUP_SURFACE_RAD

!------------------------------------------------------------------------------


!
      SUBROUTINE CLEANUP_STRAT_CLIM
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! CLEANUP_DAO begins here!
      !=================================================================
      IF ( ALLOCATED( CH4CLIM         ) ) DEALLOCATE( CH4CLIM         )
      IF ( ALLOCATED( N2OCLIM         ) ) DEALLOCATE( N2OCLIM         )
      IF ( ALLOCATED( CFC11CLIM       ) ) DEALLOCATE( CFC11CLIM       )
      IF ( ALLOCATED( CFC12CLIM       ) ) DEALLOCATE( CFC12CLIM       )
      IF ( ALLOCATED( CCL4CLIM        ) ) DEALLOCATE( CCL4CLIM        )
      IF ( ALLOCATED( CFC22CLIM       ) ) DEALLOCATE( CFC22CLIM       )


      END SUBROUTINE CLEANUP_STRAT_CLIM


      END MODULE RRTMG_RAD_TRANSFER_MOD
#endif

