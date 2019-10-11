#ifdef APM
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: apm_coag_mod
!
! !DESCRIPTION: Module APM\_COAG\_MOD contains variables and routines for
!  coagulation calculation.
!\\
!\\
! !INTERFACE:
!
      MODULE APM_COAG_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: APM_COAG
      PUBLIC :: APM_COAGSCAV
      PUBLIC :: READCK6DTABLE
      PUBLIC :: YCOAGKERN_TABLE
!
! !PRIVATE MEMBER FUNCTIONS:
!

! !REVISION HISTORY:
!  28 Nov 2008 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! MR1  : Number of points in 1st particle radius dimension
      ! MD1  : Number of points in 1st particle density dimension
      ! MR2  : Number of points in 2nd particle radius dimension
      ! MD2  : Number of points in 2nd particle density dimension
      ! MT   : Number of points in temperature dimension
      ! MP   : Number of points in pressure dimension
      INTEGER, PRIVATE, PARAMETER :: MR1=128
!      INTEGER, PRIVATE, PARAMETER :: MD1=5
      INTEGER, PRIVATE, PARAMETER :: MD1=3
      INTEGER, PRIVATE, PARAMETER :: MR2=128
!      INTEGER, PRIVATE, PARAMETER :: MD2=5
      INTEGER, PRIVATE, PARAMETER :: MD2=3
      INTEGER, PRIVATE, PARAMETER :: MT =8
      INTEGER, PRIVATE, PARAMETER :: MP =21
!
! !LOCAL VARIABLES:
!
      ! R1   : Values at points in 1st particle radius dimension
      ! D1   : Values at points in 1st particle density dimension
      ! R2   : Values at points in 2nd particle radius dimension
      ! D2   : Values at points in 2nd particle density dimension
      ! T    : Values at points in temperature dimension
      ! P    : Values at points in pressure dimension
      ! XCK  : Coagulation kernel (cm3/s) AT ALL POINTS IN 6-D SPACE
      REAL*8,  PRIVATE            :: R1(MR1)
      REAL*8,  PRIVATE            :: DEN1(MD1)
      REAL*8,  PRIVATE            :: R2(MR2)
      REAL*8,  PRIVATE            :: DEN2(MD2)
      REAL*8,  PRIVATE            :: T(MT)
      REAL*8,  PRIVATE            :: P(MP)
      REAL*8,  PRIVATE            :: XCK8(MP)
      REAL,    PRIVATE            :: XCK(MR1,MD1,MR2,MD2,MT,MP)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm_coag
!
! !DESCRIPTION: Coagulation solver (by Fangqun Yu, UAlbany, 2006, updated 2008)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE APM_COAG( ITYPE, NMAX, YT, PMB, DT, R, DEN, XN, XVA )
!
! !USES:
!
      USE APM_INIT_MOD, ONLY : VDRY
      USE APM_INIT_MOD, ONLY : VSALT
      USE APM_INIT_MOD, ONLY : VBCOC
      USE APM_INIT_MOD, ONLY : COAGPAR
      USE APM_INIT_MOD, ONLY : COAGPARSS
      USE APM_INIT_MOD, ONLY : COAGPBCOC
!
! !INPUT PARAMETERS:
!
      ! ITYPE = 1: Coagulation among sulfate particles
      !       = 2: Coagulation among sea salt particles
      INTEGER :: ITYPE

      ! Number of bins resolved for the specified type of particles
      INTEGER :: NMAX

      ! Temperature (K)
      REAL*8  :: YT

      ! Pressure (mb)
      REAL*8  :: PMB

      ! Time step (s)
      REAL*8  :: DT
!
! !INPUT/OUTPUT PARAMETERS:
!
      ! Particle wet/total radius (cm)
      REAL*8  :: R(NMAX)

      ! Density of particles (g/cm3)
      REAL*8  :: DEN

      ! Number size distribution (#/cm3)
      REAL*8  :: XN(NMAX)

      ! Volume distribution of all components (cm3/cm3)
      REAL*8  :: XVA(NMAX)
!
! !REMARKS:
!  COAGPAR = (Partition Fractions)
!  BEITA   = Coagulation Kernels
!
! !REVISION HISTORY:
!  28 Nov 2008 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      INTEGER, PARAMETER :: NC = 1
!
! !LOCAL VARIABLES:
!
      INTEGER :: N, IC, I, J, K
      REAL*8  :: XV(NMAX,NC)
      REAL*8  :: XN0(NMAX),XV0(NMAX,NC)
      REAL*8  :: BEITA(NMAX,NMAX),YCOAGPAR(NMAX,NMAX,NMAX),VDRY1(NMAX)
      REAL*8  :: YLOSS, YPRODVA(NC),TOTC
      REAL*8  :: R1,R2,YCK1
!
      IF(ITYPE.EQ.1) THEN                  ! sulfate
         YCOAGPAR(:,:,:) = COAGPAR(:,:,:)
         VDRY1(:) = VDRY
      ELSEIF(ITYPE.EQ.2) THEN   ! seasalt
         YCOAGPAR(:,:,:) = COAGPARSS(:,:,:)
         VDRY1(:) = VSALT
      ELSEIF(ITYPE.EQ.4) THEN   ! BC
         YCOAGPAR(:,:,:) = COAGPBCOC(:,:,:)
         VDRY1(:) = VBCOC
      ELSEIF(ITYPE.EQ.5) THEN   ! OC
         YCOAGPAR(:,:,:) = COAGPBCOC(:,:,:)
         VDRY1(:) = VBCOC
      ELSE
         WRITE(6,*)"STOP AT COAG, check ITYPE"
         STOP
      ENDIF

! Find coagulation Kernels
      DO I = 1, NMAX
         R1 = R(I)
         DO J = I, NMAX
            R2 = R(J)
            CALL YCOAGKERN_TABLE(R1,DEN,R2,DEN,YT,PMB,YCK1) !Use lookup table to find YCK

            BEITA(I,J) = YCK1
            BEITA(J,I) = BEITA(I,J)
         ENDDO
      ENDDO

      DO N=1,NMAX
         XV(N,1) = XVA(N)       ! to be modified if NC>1
      ENDDO
!
      DO N = 1, NMAX
         XN0(N) = XN(N)
         DO IC=1,NC
            XV0(N,IC) = XV(N,IC)
         ENDDO
      ENDDO
!
      DO i = 1, NMAX
         YLOSS = 0.
         DO j = 1, NMAX
            YLOSS=YLOSS + XN0(j)*BEITA(j,i)*(1.-YCOAGPAR(i,j,i))
         ENDDO

         DO IC=1,NC
            YPRODVA(IC) = 0.
            DO j = 1, i
               DO k = 1, i-1
                  IF(YCOAGPAR(j,k,i).GT.0.) THEN
                     YPRODVA(IC)=YPRODVA(IC) + BEITA(j,k)*XN0(j) &
                                             *XV(k,IC)*YCOAGPAR(j,k,i)
                  ENDIF
               ENDDO
            ENDDO

            XV(i,IC)=(XV0(i,IC)+YPRODVA(IC)*DT)/(1.+YLOSS*DT)
         ENDDO
      ENDDO

!
      DO N=1,NMAX
!         XVA(N) = XV(N,1) ! to be modified if NC>1
         XVA(N) = MAX(XV(N,1),1.d-40) ! to be modified if NC>1
      ENDDO
!
! Recalculate XN from XV
!
      DO I=1,NMAX
         TOTC = 0.
         DO IC = 1, NC
            TOTC = TOTC + XV(I,IC)
         ENDDO
         XN(I) = TOTC/(1.E6*VDRY1(I)) ! TOTC in cm3/cm3, VDRY in m3,XN in #/cm3
      ENDDO
!
!Mass cons. check
!
!      DO IC=1,NC
!         ZPRODVA(IC)=0.
!         ZTOTVA(IC)=0.
!      ENDDO
!      DO i = 1, NMAX
!         DO IC=1,NC
!            ZTOTVA(IC)=ZTOTVA(IC)+XV(i,IC)
!            ZPRODVA(IC)=ZPRODVA(IC)+(XV(i,IC)-XV0(i,IC))
!         ENDDO
!      ENDDO
!
!      DO IC=1,NC
!         write(95+IC,100)IC,ZTOTVA(IC),ZPRODVA(IC)
!      ENDDO
! 100  FORMAT(I3,5(1PE12.4))
 110  FORMAT(100(1E9.2))
!
      END SUBROUTINE APM_COAG
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apm_coagscav
!
! !DESCRIPTION: Subroutine to calculate the scavenging of secondary
!  particles by primary aerosols
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE APM_COAGSCAV( TK,         PMB,    DT,       DENWET, &
                               RWETCM,     YCS,    RLOSULF,  DENSALTWET, &
                               RSALTWETCM, XNSALT, & 
                               DENBCWET,RBCWET,XNBC, &
                               DENOCWET,ROCWET,XNOC, &
                               DENAER,     XNDST,  XVA,      MBCS, &
                               MOCS,       MDUSTS, MSALTS,   MBCLV, &
                               MOCLV,      MDSTLV, MSALTLV ) !inout
!
! !USES:
!
      USE APM_INIT_MOD, ONLY : DENSULF
      USE APM_INIT_MOD, ONLY : RDST
      USE APM_INIT_MOD, ONLY : NSO4
      USE APM_INIT_MOD, ONLY : NSEA
      USE APM_INIT_MOD, ONLY : NDSTB
      USE APM_INIT_MOD, ONLY : NTYP
      USE APM_INIT_MOD, ONLY : NBCOC
!
! !INPUT PARAMETERS:
!
      REAL*8  :: TK
      REAL*8  :: PMB
      REAL*8  :: DT
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8  :: DENWET
      REAL*8  :: RWETCM(NSO4)
      REAL*8  :: YCS(NTYP)
      REAL*8  :: RLOSULF
      REAL*8  :: DENSALTWET
      REAL*8  :: RSALTWETCM(NSEA)
      REAL*8  :: XNSALT(NSEA)

      REAL*8  :: DENBCWET
      REAL*8  :: RBCWET(NBCOC)
      REAL*8  :: XNBC(NBCOC)
      REAL*8  :: DENOCWET
      REAL*8  :: ROCWET(NBCOC)
      REAL*8  :: XNOC(NBCOC)

      REAL*8  :: DENAER(NTYP)
      REAL*8  :: XNDST(NDSTB)
      REAL*8  :: XVA(NSO4)
      REAL*8  :: MBCS         ! mass of sulfate coated on BC
      REAL*8  :: MOCS         ! mass of sulfate coated on OC
      REAL*8  :: MDUSTS       ! mass of sulfate coated on dust
      REAL*8  :: MSALTS       ! mass of sulfate coated on sea salt
      REAL*8  :: MSULFLV      ! mass of LV-SOA coated on primary particles
      REAL*8  :: MBCLV        ! mass of LV-SOA coated on primary particles
      REAL*8  :: MOCLV        ! mass of LV-SOA coated on primary particles
      REAL*8  :: MDSTLV       ! mass of LV-SOA coated on primary particles
      REAL*8  :: MSALTLV      ! mass of LV-SOA coated on primary particles
!
! !REVISION HISTORY:
!  28 Nov 2008 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: NI,NA,NJ
      REAL*8  :: YR1,YR2,YDEN2,YCK1
      REAL*8  :: AERLOSS(NSO4)
      REAL*8  :: AERLOSSBC(NSO4), AERLOSSOC(NSO4)
      REAL*8  :: AERLOSSDST(NSO4), AERLOSSSALT(NSO4)
      REAL*8  :: TOTV,TOTVLOSS,XDV,LOSSRATIO
      REAL*8  :: TOTVLOSSBC, TOTVLOSSOC, TOTVLOSSDST, TOTVLOSSSALT

      DO NI=1,NSO4
         AERLOSS(NI) =1.E-20
         AERLOSSBC(NI) =1.E-30
         AERLOSSOC(NI) =1.E-30
         AERLOSSDST(NI) =1.E-30
         AERLOSSSALT(NI) =1.E-30
      ENDDO
      DO NA = 2,NTYP            ! first type is secondary particles
       IF(YCS(NA).GT.1.E-5) THEN ! only when primary particles is substantial
        IF(NA.EQ.2) THEN    ! sea salt
         YDEN2 = DENSALTWET
          DO NJ=1,NSEA
            YR2 = RSALTWETCM(NJ)
             DO NI=1,NSO4
               YR1 = RWETCM(NI)
                IF(YR1.LT.YR2) THEN ! only allow scavenging of sulfate particles by bigger particles
                 CALL YCOAGKERN_TABLE(YR1,DENWET,YR2,YDEN2,TK,PMB,YCK1) !Use lookup table to find YCK
                 AERLOSSSALT(NI) = AERLOSSSALT(NI) + YCK1*XNSALT(NJ) !s-1
               ENDIF
             ENDDO
            ENDDO
            DO NI=1,NSO4
               AERLOSS(NI) = AERLOSS(NI) + AERLOSSSALT(NI)
            ENDDO
           ELSEIF(NA.EQ.3) THEN !dust
            YDEN2 = DENAER(NA)  !for now, didn't consider effect of coating (typically small)
            DO NJ=1,NDSTB
             YR2 = RDST(NJ)*100.  ! RDST in m, YR2 in cm
             DO NI=1,NSO4
               YR1 = RWETCM(NI)
               IF(YR1.LT.YR2) THEN ! only allow scavenging of sulfate particles by bigger particles
                CALL YCOAGKERN_TABLE(YR1,DENWET,YR2,YDEN2,TK,PMB,YCK1) !Use lookup table to find YCK
                AERLOSSDST(NI) = AERLOSSDST(NI) + YCK1*XNDST(NJ)    !s-1
               ENDIF
             ENDDO
            ENDDO
            DO NI=1,NSO4
               AERLOSS(NI) = AERLOSS(NI) + AERLOSSDST(NI)
            ENDDO
          ELSEIF(NA.EQ.4) THEN !BC
            YDEN2 = DENBCWET
            DO NJ=1,NBCOC
             YR2 = RBCWET(NJ)
             DO NI=1,NSO4
              YR1 = RWETCM(NI)
              IF(YR1.LT.YR2) THEN ! only allow scavenging of sulfate particles by bigger BC particles
               CALL YCOAGKERN_TABLE(YR1,DENWET,YR2,YDEN2,TK,PMB,YCK1)
               AERLOSSBC(NI) = AERLOSSBC(NI) + YCK1*XNBC(NJ) !s-1
              ENDIF
             ENDDO
            ENDDO
            DO NI=1,NSO4
               AERLOSS(NI) = AERLOSS(NI) + AERLOSSBC(NI)
            ENDDO
          ELSEIF(NA.EQ.5) THEN !OC
            YDEN2 = DENOCWET
            DO NJ=1,NBCOC
             YR2 = ROCWET(NJ)
             DO NI=1,NSO4
              YR1 = RWETCM(NI)
              IF(YR1.LT.YR2) THEN ! only allow scavenging of sulfate particles by bigger OC particles
               CALL YCOAGKERN_TABLE(YR1,DENWET,YR2,YDEN2,TK,PMB,YCK1)
               AERLOSSOC(NI) = AERLOSSOC(NI) + YCK1*XNOC(NJ) !s-1
              ENDIF
             ENDDO
            ENDDO
            DO NI=1,NSO4
               AERLOSS(NI) = AERLOSS(NI) + AERLOSSOC(NI)
            ENDDO
          ENDIF
         ENDIF
      ENDDO
      TOTV = 1.E-50
      TOTVLOSS = 1.d-50
      TOTVLOSSBC = 1.d-50
      TOTVLOSSOC = 1.d-50
      TOTVLOSSDST = 1.d-50
      TOTVLOSSSALT =  1.d-50

      DO NI=1,NSO4
         TOTV = TOTV + XVA(NI)
         XDV = XVA(NI)*(1.-exp(-AERLOSS(NI)*DT))
         XVA(NI) = XVA(NI) - XDV
         IF(XVA(NI).LT.1.d-40)XVA(NI)=1.d-40
         TOTVLOSS = TOTVLOSS + XDV

         TOTVLOSSBC = TOTVLOSSBC + XDV*AERLOSSBC(NI)/AERLOSS(NI)
         TOTVLOSSOC = TOTVLOSSOC + XDV*AERLOSSOC(NI)/AERLOSS(NI)
         TOTVLOSSDST = TOTVLOSSDST + XDV*AERLOSSDST(NI)/AERLOSS(NI)
         TOTVLOSSSALT = TOTVLOSSSALT + XDV*AERLOSSSALT(NI)/AERLOSS(NI)

      ENDDO

      LOSSRATIO = (TOTVLOSSBC + TOTVLOSSOC)/TOTVLOSS

      ! SULFATE and LVSOA MASS scavenged by BC and OC,
      ! RLOSULF is ratio of LVSOA to (LVSOA+SULF)

      ! TOTVLOSSBC in cm3/cm3  MBCS in kg/m3
      MBCS = MBCS + (1.-RLOSULF)*TOTVLOSSBC * DENSULF*1.d3

      ! MBCLV in kg/m3
      MBCLV = MBCLV + RLOSULF*TOTVLOSSBC * DENSULF*1.d3

      ! TOTVLOSSOC in cm3/cm3  MOCS in kg/m3
      MOCS = MOCS + (1.-RLOSULF)*TOTVLOSSOC * DENSULF*1.d3

      ! TOTVLOSSOC in cm3/cm3  MOCLV in kg/m3
      MOCLV = MOCLV + RLOSULF*TOTVLOSSOC * DENSULF*1.d3

      ! TOTVLOSSDUST in cm3/cm3  MDUSTS in kg/m3
      MDUSTS = MDUSTS + (1.-RLOSULF)*TOTVLOSSDST * DENSULF*1.d3

      ! TOTVLOSSDUST in cm3/cm3  MDSTLV in kg/m3
      MDSTLV = MDSTLV + RLOSULF*TOTVLOSSDST * DENSULF*1.d3

      ! TOTVLOSSSALT in cm3/cm3  MSALTS in kg/m3
      MSALTS = MSALTS + (1.-RLOSULF)*TOTVLOSSSALT * DENSULF*1.d3

      ! TOTVLOSSSALT in cm3/cm3  MSALTLV in kg/m3
      MSALTLV = MSALTLV + RLOSULF*TOTVLOSSSALT * DENSULF*1.d3

!      IF(MDSTLV.GT.1.d-8.or.MSALTLV.GT.1.d-8) THEN
      IF(MDSTLV.GT.1.d-7.or.MSALTLV.GT.1.d-7) THEN
         WRITE(6,*)"11101",RLOSULF,TOTVLOSSDST, &
              TOTVLOSSSALT,MDSTLV,MSALTLV
      ENDIF

      END SUBROUTINE APM_COAGSCAV
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ycoagkern_table
!
! !DESCRIPTION: This subroutine is to calculate coagulation kernels of two
!  particles from lookup tables (No interpolation is needed because the table
!  has high enough resolution in each dimension)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE YCOAGKERN_TABLE( X, Y, Z, U, V, W, YCK )
!
! !INPUT PARAMETERS:
!
      REAL*8  :: X      ! Radius of first particle in cm  (6.0E-8 - 1.0E-3)
      REAL*8  :: Y      ! Density of first particle in g/cm3 (1.0-2.8)
      REAL*8  :: Z      ! Radius of second particle in cm  (6.0E-8 - 1.0E-3)
      REAL*8  :: U      ! Density of second particle in g/cm3 (1.0-2.8)
      REAL*8  :: V      ! T (in K) (180-320)
      REAL*8  :: W      ! P (in mb) (10-1020)
!
! !OUTPUT PARAMETERS:
!
      REAL*8  :: YCK    ! Coagulation kernel (cm3/s)
!
! !REMARKS:
!  WRITTEN by Fangqun Yu, SUNY-Albany, 2008
!
! !REVISION HISTORY:
!  28 Aug 2008 - F. Yu       - Initial version
!  08 Nov 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER :: DRN = 30.   ! bin per decade in R space
!
! !LOCAL VARIABLES:
!
      INTEGER           :: IX1, IY1, IZ1, IU1, IV1, IW1

!
! The lookup tables should cover most of possible values in the atmosphere.
!
      ! DRN points per decade of R space
      IX1 =MAX0(INT(1.5+DRN*DLOG10(X/R1(1))),1)
      IX1 = MIN0(IX1,MR1)
!     WRITE(6,100)IX1,X,R1(IX1-1),R1(IX1),R1(IX1+1)

!      IY1 = MAX0(INT(0.5+(Y-DEN1(1))/0.4)+1,1)
      IY1 = MAX0(INT(0.5+(Y-DEN1(1))/0.8)+1,1)
      IY1 = MIN0(IY1,MD1)
!     WRITE(6,100)IY1,Y,DEN1(IY1-1),DEN1(IY1),DEN1(IY1+1)

      ! DRN points per decade of R space
      IZ1 =MAX0(INT(1.5+DRN*DLOG10(Z/R2(1))),1)
      IZ1 = MIN0(IZ1,MR2)
!     WRITE(6,100)IZ1,Z,R2(IZ1-1),R2(IZ1),R2(IZ1+1)

!      IU1 = MAX0(INT(0.5+(U-DEN2(1))/0.4)+1,1)
      IU1 = MAX0(INT(0.5+(U-DEN2(1))/0.8)+1,1)
      IU1 = MIN0(IU1,MD2)
!     WRITE(6,100)IU1,U,DEN2(IU1-1),DEN2(IU1),DEN2(IU1+1)

      IV1 = MAX0(INT(0.5+(V-T(1))/20.0)+1,1)
      IV1 = MIN0(IV1,MT)
!     WRITE(6,100)IV1,V,T(IV1-1),T(IV1),T(IV1+1)

      ! 10 per decade of P space
      IW1 = MAX0(INT(1.5+10.*DLOG10(W/P(1))),1)
      IW1 = MIN0(IW1,MP)
!     WRITE(6,100)IW1,W,P(IW1-1),P(IW1),P(IW1+1)

      YCK = DBLE(XCK(IX1,IY1,IZ1,IU1,IV1,IW1))

 100  FORMAT(I4,10(1PE10.3))

      END SUBROUTINE YCOAGKERN_TABLE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readck6dtable
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READCK6DTABLE
!
! !USES:
!
      USE APM_INIT_MOD, ONLY : DATA_DIR_1x1
!
! !REMARKS:
!  WRITTEN by Fangqun Yu, SUNY-Albany, 2008
!
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
      INTEGER :: IR1,ID1,IR2,ID2,IT,IP
      REAL*8  :: DRN, YY

      CHARACTER*999 YPATH
      YPATH = TRIM(DATA_DIR_1x1)//'/APM_data_201906/YCK140609/'
      WRITE(6,*)"Read coagulation kernel look-up tables"

      open(41,file=TRIM(YPATH)//'Yu_CK_1R1.txt',form='formatted')
      open(42,file=TRIM(YPATH)//'Yu_CK_2D1.txt',form='formatted')
      open(43,file=TRIM(YPATH)//'Yu_CK_3R2.txt',form='formatted')
      open(44,file=TRIM(YPATH)//'Yu_CK_4D2.txt',form='formatted')
      open(45,file=TRIM(YPATH)//'Yu_CK_5T.txt',form='formatted')
      open(46,file=TRIM(YPATH)//'Yu_CK_6P.txt',form='formatted')
      open(47,file=TRIM(YPATH)//'Yu_CK_7CK.txt',form='formatted')

!
      READ(41,100)(R1(IR1),IR1=1,MR1)
      WRITE(6,*)"R1(I), I=1, ", MR1, ":"
      WRITE(6,100)(R1(IR1),IR1=1,MR1)
!
      READ(42,100)(DEN1(ID1),ID1=1,MD1)
      WRITE(6,*)"DEN1(I), I=1, ", MD1, ":"
      WRITE(6,100)(DEN1(ID1),ID1=1,MD1)
!
      READ(43,100)(R2(IR2),IR2=1,MR2)
      WRITE(6,*)"R2(I), I=1, ", MR2, ":"
      WRITE(6,100)(R2(IR2),IR2=1,MR2)
!
      READ(44,100)(DEN2(ID2),ID2=1,MD2)
      WRITE(6,*)"DEN2(I), I=1, ", MD2, ":"
      WRITE(6,100)(DEN2(ID2),ID2=1,MD2)
!
      READ(45,100)(T(IT),IT=1,MT)
      WRITE(6,*)"T(I), I=1, ", MT, ":"
      WRITE(6,100)(T(IT),IT=1,MT)
!
      READ(46,100)(P(IP),IP=1,MP)
      WRITE(6,*)"P(I), I=1, ", MP, ":"
      WRITE(6,100)(P(IP),IP=1,MP)
!
! Use the formula to calculate parameters to get values with more digits,
! otherwise may cause problem when input values are very clsoe to the
! parameter values as formula is used to determine the location of inputted
! values. Also serve as a double check to make sure the consisistency in
! parameter spaces.
!
      DRN = 30.
      R1(1) = 6.0E-8            ! cm
      DO IR1 = 2, MR1
         YY = R1(IR1)
         R1(IR1)=R1(1)*10.**(float(IR1-1)/DRN)
         IF(abs(1.-YY/R1(IR1)).GT.0.02) THEN
            write(6,*)"need check CK look-up table R1 inputs"
            stop
         ENDIF
      ENDDO

      DEN1(1) = 1.0             !g/cm3
      DO ID1 = 2, MD1
         YY = DEN1(ID1)
         DEN1(ID1) = DEN1(1)+0.8*float(ID1-1)
         IF(abs(1.-YY/DEN1(ID1)).GT.0.02) THEN
            write(6,*)"need check CK look-up table DEN1 inputs"
            stop
         ENDIF
      ENDDO

      R2(1) = 6.0E-8            ! cm
      DO IR2 = 2, MR2
         YY = R2(IR2)
         R2(IR2)=R2(1)*10.**(float(IR2-1)/DRN)
         IF(abs(1.-YY/R2(IR2)).GT.0.02) THEN
            write(6,*)"need check CK look-up table R2 inputs"
            stop
         ENDIF
      ENDDO

      DEN2(1) = 1.0             !g/cm3
      DO ID2 = 2, MD2
         YY = DEN2(ID2)
         DEN2(ID2) = DEN2(1)+0.8*float(ID2-1)
         IF(abs(1.-YY/DEN2(ID2)).GT.0.02) THEN
            write(6,*)"need check CK look-up table DEN2 inputs"
            stop
         ENDIF
      ENDDO

      T(1) = 180.               !k
      DO IT = 2, MT
         YY = T(IT)
         T(IT) = T(1)+20.*float(IT-1)
         IF(abs(1.-YY/T(IT)).GT.0.02) THEN
            write(6,*)"need check CK look-up table T inputs"
            stop
         ENDIF
      ENDDO

      P(1) = 10.                ! mb
      DO IP = 2, MP
         YY = P(IP)
         P(IP) = P(1) * 10.**(0.1*float(IP-1))
         IF(abs(1.-YY/P(IP)).GT.0.02) THEN
            write(6,*)"need check CK look-up table P inputs"
            stop
         ENDIF
      ENDDO

!
! Read in the 6-D Table
!
      DO IR1 = 1, MR1
      DO ID1 = 1, MD1
      DO IR2 = 1, MR2
      DO ID2 = 1, MD2
      DO IT  = 1, MT
         READ(47,100)(XCK8(IP),IP=1,MP)
         DO IP=1,MP
            XCK(IR1,ID1,IR2,ID2,IT,IP)=REAL(XCK8(IP))
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!
 100  FORMAT(1PE8.2,200(1PE9.2))

      CLOSE(41)
      CLOSE(42)
      CLOSE(43)
      CLOSE(44)
      CLOSE(45)
      CLOSE(46)
      CLOSE(47)
!
      END SUBROUTINE READCK6DTABLE
!EOC
      END MODULE APM_COAG_MOD
#endif
