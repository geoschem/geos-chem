C $Id: diffuse.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
       SUBROUTINE DIFFUSE(ISTRIP,NPCS,NDTURB,
     &                    ITRTRB,PZ,TMP,SPHU,
     &                    IM,JM,LM,NM,SIGE,SIG,DSIG,PTOP,
     &                    KDIFF,TC,XTRA2)

c-----------------------------------------------------------------------
c      In this subroutine the geometrical height
c      used in the diffusion equation is derived from the
c       sigma-coordinate system.
c      In the first part the transformation is
c      performed, then the diffusion eq is solved is called in TRBDIF1.
c
c      istrip  - number of horizontal points to be handled at a time on
c      npcs    - number of regions in which the domain is diveded
c      ndturb  - turbulence time step integer in HHMMSS format
c      pz      - surface pressure minus ptop in mb - real[lon,lat]
c      tz      - temperature [deg K]  - level 1 at top
c      qz      - specific humidity in g/kg - level 1 at top
c      im      - number of points in the longitude direction
c      jm      - number of points in the latitude direction
c      lm      - number of vertical levels
c      nm      - number of species
c      sig     - model mid-sigma values - level =1 at top
c      sige    - model upper edge-sigma values - level =1 at top
c      dsig    - model edge-to-edge delta sigma values - level =1 at top
c      ptop    - model top pressure - rigid lid assumed - real
c
c      kdiff     - diffusion coefficients - level =1 at top
c updated:
c      tc     - quantity  to be diffused - level =1 at top
c
c
C  Here is the mapping of input arguments from main.f (bdf, bmy, 2/24/99)
C     ISTRIP  <-- JJPAR or JM
C     NPCS    <-- IIPAR or IM
C     NTDT    <-- NTDT
C     ITRTRB  <-- 1
C     PSC     <-- PSC (or PW) = surface pressure - PTOP
C     TMP     <-- T
C     SPHU    <-- SPHU
C     IM      <-- IIPAR or IM
C     JM      <-- JJPAR or JM
C     LM      <-- LLPAR or LM
C     NM      <-- NTRACE
C     SIGEf   <-- G_SIGEf
C     SIGf    <-- G_SIGf
C     DSIG    <-- G_DSIGEf
C     PTOP    <-- PTOP
C     KDIFF   <-- KZZ
C     TC      <-- STT
C     XTRA2   <-- XTRA2
c----------------------------------------------------------------------------
C  CHANGES IN DIFFUSE.F
C
C 1)  XTRA2 gives the height of the boundry layer.  The kzz values are not
C     trusted above the boundry layer, as the generation of the values comes
C     from a boundry layer turbulence code.  Therefore we zero all kzz values
C     above the boundry layer.  (bdf 3/3/99)
C 
C 2)  change 1-d loops over 2-d arrays to 2-d loops.  This vectorization 
C     works in cray fortran, but not with our compiler.
C
C 3)  GETCON is no longer called...define parameters directly (bmy, 3/23/99)
C
C 4)  Now use double precision exponents (e.g. 1D0, to ensure that we 
C     don't lose precision to roundoff errors (bmy, 3/23/99)
C
C 5)  Make the following changes (bmy, 4/6/99):
C     (1) TH(J,L) is replaced by TZ(NN,J,L)
C     (2) SH(J,L) is replaced by QZ(NN,J,L)
C     (3) KH(J,L) is replaced by KDIFF2(NN,J,L)
C     (4) Arrays PK, PKZ, PLZ, TEMP are eliminated
C     (5) Call to PKAPPA.F is eliminated -- PKZ isn't used anywhere!!!
C     (6) Define new parameters to avoid repeated mult/div in do loops
C 
C 6) Fixed bug in do loop for FLXFAC (bdf,  4/8/99)
C   
C 7) XTRA2(IREF,JREF,5) is now XTRA2(I,J) (bmy, 9/25/01)
C
C 8) Removed obsolete code from 9/01 (bmy, 10/23/01)
C
C----------------------------------------------------------------------------
C  Subroutines needed for use with DIFFUSE:
C  ========================================
C  (1) diffuse.f  (2) TRBDIF.f  (3) VTRI0.f
C----------------------------------------------------------------------------
      IMPLICIT NONE
      
c Input Variables
c ---------------
      INTEGER, INTENT(IN)    :: NDTURB,ISTRIP,NPCS

      INTEGER, INTENT(IN)    :: IM,JM,LM,NM       ! Physics Grid
      REAL*8,  INTENT(IN)    :: PTOP              ! Physics Grid
      REAL*8,  INTENT(IN)    :: SIGE(LM+1)        ! Physics Grid
      REAL*8,  INTENT(IN)    :: SIG(LM)           ! Physics Grid
      REAL*8,  INTENT(IN)    :: DSIG(LM)          ! Physics Grid

      REAL*8,  INTENT(IN)    :: PZ(IM,JM)         ! Dynamics State
      REAL*8,  INTENT(IN)    :: TMP(IM,JM,LM)
      REAL*8,  INTENT(IN)    :: SPHU(IM,JM,LM)
      REAL*8,  INTENT(IN)    :: KDIFF(IM,JM,LM)   ! Diffusion coefficients
      REAL*8,  INTENT(IN)    :: XTRA2(IM,JM)   
      REAL*8,  INTENT(INOUT) :: TC(IM,JM,LM,NM)   ! quantity to be diffused

c Local Variables
c ---------------
!**** Work arrays
      REAL*8                 :: TZ(IM,JM,LM)      ! Dynamics State
      REAL*8                 :: QZ(IM,JM,LM)      ! Dynamics State
      REAL*8                 :: PKHT(IM,JM,LM)    ! Dynamics State
      REAL*8                 :: KDIFF2(IM,JM,LM)  ! zeroed kzz values
      REAL*8                 :: XX(JM,LM+1)       ! Temp storage for TC

!**** The following arrays are not necessary (bmy, 4/5/99)
!****      REAL*8                 :: TH(ISTRIP,lm)
!****      REAL*8                 :: SH(ISTRIP,lm) 
!****      REAL*8                 :: PKZ(IM,JM,lm)
!****      REAL*8                 :: KH(ISTRIP,lm)  
!****      REAL*8                 :: TEMP(ISTRIP,LM)
!****      REAL*8                 :: P(ISTRIP,LM)
!****      REAL*8                 :: PK(ISTRIP,LM) 
!****      REAL*8                 :: PLZ(IM,JM,LM)
      REAL*8                 :: THV(ISTRIP,LM) 
      REAL*8                 :: PE(ISTRIP,LM+1)
      REAL*8                 :: PKE(ISTRIP,LM+1)

c-----Dimensions for the variables used to convert the grid
      REAL*8                 :: ADZ1(ISTRIP,LM )
      REAL*8                 :: ADZ2(ISTRIP,LM-1)
      REAL*8                 :: RHODZ2(ISTRIP,LM-1)
      REAL*8                 :: RHOKDZ(ISTRIP,LM)
      REAL*8                 :: FLXFAC(ISTRIP,LM+1)

!**** Scalar variables
      LOGICAL                :: LTYPE
      
      REAL*8                 :: DTTRB,TIMSTP,ATIMSTP
      REAL*8                 :: FAC1,FAC2

!**** bmy added GD1 to save on repeated operations (bmy, 4/6/99)
      REAL*8                 :: GD1

      INTEGER                :: IT,ITRC
      INTEGER                :: LP1,ITER
      INTEGER                :: I,J,L,NN,ITRTRB
      INTEGER                :: NLEVM1,NLEVP1
!**** These variables are loop bounds for Cray 1-D loops (bmy, 4/6/99)
!****      INTEGER                :: IMJM,ISTNP1,ISTNLAY,ISTNM1

!**** Parameters -- GETCON is obsolete (bmy, 4/6/99)
      REAL*8, PARAMETER      :: CP      = 1004.16d0
      REAL*8, PARAMETER      :: GRAV    = 9.81d0
      REAL*8, PARAMETER      :: VIRTCON = 0.609d0
      REAL*8, PARAMETER      :: RGAS    = 8314.3d0 / 28.97d0
      REAL*8, PARAMETER      :: AKAP    = RGAS     / CP
      REAL*8, PARAMETER      :: EPS     = 1d-30

!**** Add extra parameters to eliminate repeated operations (bmy, 4/5/99)
      REAL*8, PARAMETER      :: CP_GRAV = CP / GRAV
      REAL*8, PARAMETER      :: RGAS01  = RGAS * 0.01d0
C
C   SET VARIABLES THAT DO NOT CHANGE
C
      TIMSTP = NDTURB

!**** These variables are loop bounds for Cray 1-D loops (bmy, 4/6/99)
!****      ISTNP1   =   ISTRIP * (lm+1)
!****      ISTNLAY  =   ISTRIP * lm
!****      ISTNM1=istrip*NLEVM1

      LTYPE  = .TRUE.
      NLEVM1 = LM-1
      NLEVP1 = LM+1
c     
c   Convert and flip up-side-down to have edge=lm+1 at earth surface
c

!**** zero all kzz values above the boundry layer. 
!**** (kzz values have been flipped)  (bdf, 4/5/99)
      DO J=1,JM
         DO I=1,IM
            L = CEILING( XTRA2( I, J ) )
            KDIFF2( I, J, (LM-L) : LM  ) = KDIFF( I, J, (LM-L) : LM )
            KDIFF2( I, J, 1 : (LM-L-1) ) = 0d0
         ENDDO
      ENDDO

      DO L=1,LM
         DO J=1,JM
            DO I=1,IM
c
c Specific humidity, convert g/kg -> kg/kg
c
               QZ(I,J,L) = SPHU(I,J,L)*0.001d0
c
c Potential temperature T*(p0/p)**kappa
c with  p0=1 mb and kappa=Rd/cp
c --> theta=T/p**kappa (about 40 K at the earth surface)
c
               TZ(I,J,L) = TMP(I,J,L)/(SIG(L)*PZ(I,J) + PTOP )**AKAP
c
c Model pressure at upper model edges
c
               PKHT(I,J,L) = (SIGE(L+1)*PZ(I,J) + PTOP )**AKAP
            ENDDO
         ENDDO
      ENDDO
c
      DTTRB   = TIMSTP / FLOAT(ITRTRB)
      ATIMSTP = 1.d0 / TIMSTP

!**** Save GRAV * DTTRB * 0.1d0 in a variable (bmy, 4/5/99)
      GD1     = GRAV * DTTRB * 0.01d0

!**** We don't need to call PKAPPA anymore (bmy, 4/8/99)
!****C **********************************************************************
!****C   COMPUTE P ** KAPPA AT LAYER CENTERS
!****C **********************************************************************
!****
!****      call pkappa ( pz,pkht,pkz,ptop,sige,dsig,im,jm,lm )

c********************************************************
c     LOOP OVER THE REGIONS
!**** Recall that NN is the longitude index (bmy, 4/5/99)
c*********************************************************
      DO 2000 NN = 1, NPCS
!**** TH(J,L) is not needed...replace by TZ(NN,J,L) below AND
!**** SH(J,L) is not needed...replace by QZ(NN,J,L) below AND
!**** PK is not used anywhere!!! AND
!**** KH(J,L) is replaced by KDIFF2(NN,J,L) below (bmy, 4/6/99)
!****         TH (:,:)      = TZ    (NN,:,:)
!****         SH (:,:)      = QZ    (NN,:,:)
!****         PK (:,:)      = PKZ   (NN,:,:)
!****         KH (:,:)      = KDIFF2(NN,:,:)
!**** Surface pressure PZ(I,J) is stored in PE(J,LM+1) (bmy, 4/5/99)
         PE (:,LM+1)   = PZ    (NN,:  )
!**** PKHT(I,J,L) is stored in PKE(J,L,2:LM+1) (bmy, 4/5/99)
         PKE(:,2:LM+1) = PKHT  (NN,:,:)  

c*****

! PE  is pressure at sigma edges
! PKE is P**Kappa at sigma edges
         DO I =1,ISTRIP
            PE(I,1)     = ( SIGE(1)*PE(I,LM+1)+PTOP )
            PKE(I,1)    = ( SIGE(1)*PE(I,LM+1)+PTOP )**AKAP
         ENDDO   

         DO L = 1, LM
         DO I = 1, ISTRIP
!**** P is not used anywhere! (bmy, 4/6/99)
!****               P(I,L)    = SIG (L)   * PE(I,LM+1) + PTOP
            PE(I,L+1) = SIGE(L+1) * PE(I,LM+1) + PTOP
         ENDDO   
         ENDDO   


C COMPUTE VIRTUAL POTENTIAL TEMPERATURES
C --------------------------------------
         DO L = 1, LM
         DO J = 1, JM
!**** Replace SH(J,L) with QZ(NN,J,L) below and also
!**** replace TH(J,L) by TZ(NN,J,L) below (bmy, 4/5/99)
!****               THV(J,L) = 1d0 + VIRTCON * SH(J,L)
!****               THV(J,L) = TH(J,L) * THV(J,L)
            THV(J,L) = 1d0 + VIRTCON * QZ(NN,J,L) 
            THV(J,L) = TZ(NN,J,L) * THV(J,L)
         ENDDO
         ENDDO

C ---------
c*******************************************************************
c    Convert the vertical coordinate from a sigma system to 
c    cartesian one

C COMPUTE VERTICAL GRID
C ---------------------
         DO L = 1, LM
         DO J = 1, JM
!**** Replace (CP/GRAV) by CP_GRAV, to avoid repeated divisions (bmy, 4/5/99)
!****               ADZ1(J,L) = (CP/GRAV)*(PKE(J,L+1)-PKE(J,L))
            ADZ1(J,L) = CP_GRAV  * (PKE(J,L+1)-PKE(J,L))
            ADZ1(J,L) = THV(J,L) * ADZ1(J,L)
         ENDDO
         ENDDO

         DO L = 1, NLEVM1
         DO J = 1, JM
            ADZ2(J,L) = 0.5d0 * (ADZ1(J,L)+ADZ1(J,L+1))
         ENDDO
         ENDDO

C COMPUTE RHO BY DZ AT MID AND EDGE LEVELS
c (but only for the levels 1,lm-1
C ----------------------------------------
         DO 200 L = 1,NLEVM1
            LP1 = L + 1
            FAC1 = DSIG(L) / ( DSIG(L) + DSIG(LP1) )
            FAC2 = 1.d0 - FAC1
            DO 9058 I =1,ISTRIP
               RHODZ2(I,L) = FAC1 * THV(I,LP1)
               RHODZ2(I,L) = RHODZ2(I,L) + FAC2 * THV(I,L)
 9058       CONTINUE
 200     CONTINUE
c***
         DO L = 1, NLEVM1
         DO J = 1, JM
!**** Replace RGAS*0.01d0 by RGAS01 to eliminate repetitions (bmy, 4/5/99)
!****               RHODZ2(J,L) = (RGAS*0.01d0) * RHODZ2(J,L)
            RHODZ2(J,L) = RGAS01 * RHODZ2(J,L)
!**** We don't need TEMP(J,L) in this expression (bmy, 4/6/99)
!****               TEMP(J,L) = PKE(J,L+1) * ADZ2(J,L)
!****               RHODZ2(J,L) = TEMP(J,L) * RHODZ2(J,L)
            RHODZ2(J,L) = ( PKE(J,L+1) * ADZ2(J,L) ) * RHODZ2(J,L)
            RHODZ2(J,L) = PE(J,L+1)  / RHODZ2(J,L)
         ENDDO
         ENDDO

C COMPUTE FLXFAC FOR LAYERS AND EDGES
C -----------------------------------------------------------------
         DO L = 1, LM
         DO J = 1, JM
!**** Replace GRAV * DTTRB * 0.01d0 with parameter GD1 (bmy, 4/5/99)
!****               FLXFAC(J,L) = (GRAV*DTTRB*0.01d0) / (PE(J,L+1)-PE(J,L))
            FLXFAC(J,L) = GD1 / (PE(J,L+1)-PE(J,L))
         ENDDO
         ENDDO

         DO I =1,ISTRIP
            FLXFAC(I,NLEVP1) = 0.d0
         ENDDO

c-----COMPUTE RHOKDZ
c*********************************************************************
c     In this version RHODZ2 for levels 1,NLEV-1 is calculated
c     above, while RHODZ2 for level NLEV is not needed
c     as rhokdz=rhodz2*kh is set = 0 for level NLEV
c     from the file rhodz2.dat
c*******************************************************************
         DO L = 1, NLEVM1
         DO J = 1, JM
!**** Replace KH(J,L) by KDIFF2(NN,J,L) (bmy, 4/5/99)
!****               RHOKDZ(J,L) = RHODZ2(J,L) * KH(J,L)
            RHOKDZ(J,L) = RHODZ2(J,L) * KDIFF2(NN,J,L)
         ENDDO
         ENDDO

C
         DO J = 1, JM
            RHOKDZ(J,LM) = 0.0d0
         ENDDO

         DO 3010 ITRC = 1, NM

!**** Store a Lat-Alt slice of TC for this longitude and tracer (bmy, 4/6/99)
            XX(1:JM, 1:LM) = TC(NN, 1:JM, 1:LM, ITRC)
            XX(1:JM, LM+1) = 0d0

C LOOP ON TIME
C
            DO 3000 ITER = 1, ITRTRB
C
c-----Solves diffusion equation using a backward implicit scheme
C
               CALL TRBDIF(XX,RHOKDZ,FLXFAC,LM,LTYPE,EPS,ISTRIP)

 3000       CONTINUE

!**** Store the new concentrations back into TC (bmy, 4/6/99)
            TC(NN, 1:JM, 1:LM, ITRC) = XX(1:JM, 1:LM)

 3010    CONTINUE

 2000 CONTINUE

C**********************************************************************
C   END REGIONS LOOP and exit this subroutine with a new
c   value of XX
c**********************************************************************
      RETURN
      END SUBROUTINE DIFFUSE 
