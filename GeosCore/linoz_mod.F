!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: linoz_mod
!
! !DESCRIPTION: Module LINOZ\_MOD contains routines to perform the Linoz 
!  stratospheric ozone chemistry.
!\\
!\\
! !INTERFACE:
!
      MODULE LINOZ_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !DEFINED PARAMETERS:
!
      INTEGER, PARAMETER  :: NFIELDS_LINOZ = 7   ! # of Linoz fields
      INTEGER, PARAMETER  :: NLEVELS_LINOZ = 25  ! # of levels in Linoz fields
      INTEGER, PARAMETER  :: NLAT_LINOZ    = 18  ! # latitudes in Linoz fields
      INTEGER, PARAMETER  :: NMONTHS_LINOZ = 12  ! # of months in Linoz fields
!
! !PRIVATE DATA MEMBERS:
!
      REAL*8, ALLOCATABLE :: TPARM(:,:,:,:)
      REAL*8, ALLOCATABLE :: TLSTT(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_LINOZ
      PUBLIC  :: DO_LINOZ    
      PUBLIC  :: LINOZ_READ
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: INIT_LINOZ
      PRIVATE :: LINOZ_CHEM3
      PRIVATE :: LINOZ_STRATL
      PRIVATE :: LINOZ_STRT2M
      PRIVATE :: LINOZ_SOMLFQ
      PRIVATE :: LINOZ_INTPL
      PRIVATE :: STRAT_INIT
!
! !REMARKS:
!  Dylan Jones (dbj@atmosp.physics.utoronto.ca) wrote:
!                                                                             .
!    Testing this code [in v8-02-04] was more difficult that I thought. 
!    I began by trying to compare the output of v8-02-04 with our previous 
!    runs with v8-02-01.  I accounted for the changes in the transport_mod.f 
!    and I tried to undo the changes in when the diagnostics are archived in 
!    v8-02-04, but I was still getting large differences between v8-02-04 
!    and v8-02-01. I finally gave up on this since I may have made a mistake
!    in reverting to the old way of doing the diagnostics in v8-02-04.  In 
!    the end I took the new linoz code from v8-02-04 and used it in v8-02-01. 
!    I ran two GEOS-5 full chemistry simulations for 2007 and the output 
!    were consistent over the full year.  
!                                                                             .
!    I think that it is safe to release [Linoz in v8-02-04].  However, we 
!    should acknowledge that it was [only] tested in v8-02-01, since I was 
!    not able to assess the quality of the output in v8-02-04. 
!
! REVISION HISTORY:
!  23 Mar 2000 - P. Cameron-Smith    - Initial version adapted heavily 
!                                      from McLinden's original file.
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!  28 May 2009 - D. Jones            - Further modifications
!  18 Nov 2009 - D. Jones            - Further modifications
!EOP
!------------------------------------------------------------------------------
!BOC

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_linoz
!
! !DESCRIPTION: Subroutine DO\_LINOZ is the main driver for the Linoz
!  stratospheric Ozone chemistry package.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_LINOZ   
!
! !USES:
!
      USE TIME_MOD

#     include "F77_CMN_SIZE"
!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST     = .TRUE.
      INTEGER, SAVE :: LASTMONTH = -99
      REAL*8        :: NSCHEM

      ! if new month, get new parameters?
      IF ( GET_MONTH() /= LASTMONTH ) THEN
         CALL LINOZ_STRATL
         LASTMONTH =  GET_MONTH()
      ENDIF

      ! Linoz needs time step in seconds
      NSCHEM = GET_TS_CHEM() * 60D0 

      ! Call the Linoz chemistry
      CALL LINOZ_CHEM3( NSCHEM )

      END SUBROUTINE DO_LINOZ    
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: linoz_chem3
!
! !DESCRIPTION: Subroutine LINOZ\_CHEM3 applies linearized chemistry based on 
!  tables from PRATMO model using climatological T, O3, time of year
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LINOZ_CHEM3( DTCHEM )
!
! !USES:
!
      USE DAO_MOD                               
      USE GRID_MOD,       ONLY : GET_AREA_CM2   
      USE TRACER_MOD
      USE TRACERID_MOD                          
      USE TROPOPAUSE_MOD, ONLY : GET_TPAUSE_LEVEL 
      USE TROPOPAUSE_MOD, ONLY : GET_MAX_TPAUSE_LEVEL
      USE PRESSURE_MOD,   ONLY : GET_PEDGE
      USE PRESSURE_MOD,   ONLY : GET_PCENTER
  
#     include "F77_CMN_SIZE"
#     include "CMN"
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN) :: DTCHEM    ! Time step [seconds]
!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!  18 Nov 2009 - D. Jones            - For now, set tagged stratospheric 
!                                      tracer to total O3 in the overworld 
!                                      to avoid issues with spin ups
!  08 Feb 2010 - R. Yantosca         - Deleted obsolete local variables
!  22 Oct 2010 - R. Yantosca         - Added OMP parallel loop
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER :: IM,      JM,         LM
      INTEGER :: I,       J,          L,    N
      INTEGER :: K,       M,          LBOT, L_OVERWRLD
      INTEGER :: NTRACER, NUM_TRACER, LPOS, ITRC
      REAL*8  :: CLIMO3,  CLIMPML,    DCO3, DERO3, DERTMP
      REAL*8  :: DERCO3,  DMASS,      DTMP, SSO3

      ! Arrays
      REAL*8  :: DCOLO3(IIPAR,JJPAR,LLPAR)
      REAL*8  :: COLO3(IIPAR,JJPAR,LLPAR)
      REAL*8  :: OUT_DATA(IIPAR,JJPAR,LLPAR)
      
      ! Assign values for local IM and JM   (dbj 6/24/03) 
      IM     = IIPAR
      JM     = JJPAR
      LM     = LLPAR   
      L_OVERWRLD = GET_MAX_TPAUSE_LEVEL()

      ! Stratospheric Chemistry Tables for O3:
      ! ======================================
      !   7 tables, each a function of month (12), latitude 
      !   (18, -85 to 85 in 10 deg. increments) and altitude 
      !   (25, z*=10-58 km in 2 km increments).
      !  1- ozone (Logan climatology), v/v
      !  2- Temperature climatology, K
      !  3- Column ozone climatology, Logan ozone integrated above box, DU
      !  4- ozone (P-L) for climatological ozone, v/v/s
      !  5- d(P-L) / dO3, 1/s
      !  6- d(P-L) / dT, v/v/s/K      
      !  7- d(P-L) / d(column O3), v/v/s/DU
      !
      ! zero storage arrays
      !      do n=1,ntrace
      !         sttold(n)=0.d0
      !      enddo

      !=================================================================
      ! Select the proper tracer number to store O3 into, depending on
      ! whether this is a full chemistry run or a tagged Ox run.
      ! If tagged Ox, tracer 2 should be the stratospheric tracer.  (dbj)
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN
          NUM_TRACER = 1
      ELSE
         IF ( ITS_A_TAGOX_SIM() ) THEN
           IF (N_TRACERS > 1) THEN
             NUM_TRACER = 2
           ELSE
             NUM_TRACER = 1
           ENDIF
         ELSE
            ! All other simulations don't use O3...print error message
            WRITE( 6, '(a)' ) 'This simulation does not use O3!!'
            WRITE( 6, '(a)' ) 'STOP in linoz_chem3.f!'
            STOP
         ENDIF
      ENDIF

      ! Echo info
      WRITE( 6, 100 ) 
 100  FORMAT('     - LINOZ_CHEM3: Doing LINOZ stratospheric chemistry')

      ! **** note dbj: check STT(I,J,20:LLPAR,NTRACER) = with trop level
      ! ****         : check DMASS

      DO ITRC=1,NUM_TRACER     ! dbj loop for tagged

         IF ( ITS_A_FULLCHEM_SIM() ) THEN
            NTRACER = IDTOX
         ELSE
            NTRACER = ITRC
         ENDIF

         ! Start at top layer and continue to lowest layer for strat. chem
         OUT_DATA = 0d0

!------------------------------------------------------------------------------
! Prior to 2/7/11:
! Un-parallelize until further notice (ccarouge, bmy, 2/7/11)
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I,       J,     LBOT,   LPOS,   L    )
!!$OMP+PRIVATE( CLIMPML, DERO3, CLIMO3, DERCO3, DCO3 )
!!$OMP+PRIVATE( DERTMP,  DTMP,  SSO3,   DMASS        )
!------------------------------------------------------------------------------
         DO J = 1, JM
         DO I = 1, IM
            LBOT = GET_TPAUSE_LEVEL(I,J)+1
            LPOS = 1
            DO WHILE (GET_PEDGE(I,J,LPOS+1) .GE. 0.3D0)
              LPOS = LPOS +1
            ENDDO
            LPOS = LPOS-1

            !---------------------------------------------------------
            ! dbj: for now, set tagged stratospheric tracer to total 
            ! O3 in the overworld to avoid issues with spin ups
            !---------------------------------------------------------
            IF ( ITS_A_TAGOX_SIM() ) THEN
              STT(I,J,(L_OVERWRLD+1):LLPAR,NTRACER) =
     &                STT(I,J,(L_OVERWRLD+1):LLPAR,1)
            ENDIF


            DO L = LM,LBOT,-1        

               !IF (L .LT. LPAUSE(I,J)) GOTO 31

               IF (STT(I,J,L,NTRACER) .LE. 0.D0) CYCLE     

               ! calculate ozone column above box (and save)
               ! dcolo3 = ozone column (in DU) in given layer
               ! colo3 =  ozone column above layer + half of
               ! column in layer

               ! bdf stt is in v/v, make conversion to DU
               if (l.eq.lm) then !top model layer
                  dcolo3(i,j,l) = (stt(i,j,l,NTRACER)*AD(I,J,L)/
     &                 TCVV(NTRACER))/ GET_AREA_CM2(J) *
     &                 6.022d23/(28.97/TCVV(NTRACER)*1d-3)/ 2.687d16
                  colo3(i,j,l) = dcolo3(i,j,l)*0.5
               else
                  dcolo3(i,j,l) = (stt(i,j,l,NTRACER)*AD(I,J,L)/
     &                 TCVV(NTRACER))/ GET_AREA_CM2(J) *
     &                 6.022d23/(28.97/TCVV(NTRACER)*1d-3)/ 2.687d16
                  colo3(i,j,l) = colo3(i,j,l+1) +
     &                 (dcolo3(i,j,l)+dcolo3(i,j,l+1))*0.5
               endif
               out_data(i,j,l) = colo3(i,j,l)

               ! ++++++ climatological P-L:   ++++++          
               climpml=tlstt(j,l,4)      ! Climatological P-L = (P-L)^o

               ! ++++++ local ozone feedback: ++++++ 
               dero3=tlstt(j,l,5)               ! Derivative w.r.t. O3. 
                                                !  dero3=-1/(time constant)
               IF (dero3.EQ.0) CYCLE            ! Skip Linoz if lifetime 
                                                !  is infinite.
               climo3=tlstt(j,l,1)              ! Climatological O3 = f^o
               derco3=tlstt(j,l,7)              ! Derivative w.r.t. Column O3
               dco3=(colo3(i,j,l)-tlstt(j,l,3)) ! deviation from o3 climatology.
               ! ++++++ temperature feedback: ++++++ 
               dertmp=tlstt(j,l,6)              ! Derivative w.r.t. Temperature
               dtmp=(T(I,J,L)-tlstt(j,l,2))     ! Deviation in Temperature 
                                                !  from climatology.

               ! ++++++ calculate steady-state ozone: ++++++ 
               sso3=climo3 - (climpml+dtmp*dertmp+dco3*derco3)/dero3 

               ! ++++++ change in ozone mass due to chemistry: ++++++
               !ssO3 = f^*
               dmass=(sso3-STT(I,J,L,NTRACER))*(1.0-exp(dero3*dtchem))

               ! ++++++ update ozone mass ++++++ 
               ! LINOX valid only up to 58 km, so do not use above 0.3 hPa
               ! dbj: impose exponential fall off of mixing ratio 
               ! between 0.3 and 0.01 hPa (with fall off of a scale height)
               IF (GET_PEDGE(I,J,L) .LE. 0.3D0) THEN
                  STT(I,J,L,NTRACER) = 
     &                (GET_PCENTER(I,J,L)/GET_PCENTER(I,J,LPOS-1)) * 
     &                STT(I,J,LPOS-1,NTRACER)
               ELSE
                  STT(I,J,L,NTRACER) = STT(I,J,L,NTRACER)+DMASS 
               ENDIF

            ENDDO       ! loop over L

         ENDDO          ! loop over I

      ENDDO             ! loop pver J
!------------------------------------------------------------------------------
! Prior to 2/7/11:
! Un-parallelize until further notice (ccarouge, bmy, 2/7/11)
!!$OMP END PARALLEL DO
!------------------------------------------------------------------------------

      !write our calculated column o3 maximum
      !write(6,*) 'max of columns= ',maxval(out_data)

      ENDDO   ! loop over ntracers

      END SUBROUTINE LINOZ_CHEM3
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: linoz_stratl
!
! !DESCRIPTION: Subroutine LINOZ\_STRATL performs a monthly fixup of chemistry 
!  parameters for the Linoz stratospheric ozone chemistry.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LINOZ_STRATL
!
! !USES:
!
      USE GRID_MOD,       ONLY : GET_YMID 
      USE TIME_MOD,       ONLY : GET_MONTH
      USE PRESSURE_MOD                   

#     include "F77_CMN_SIZE"
#     include "CMN"
!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer J,K,L,N,JLATMD(jjpar),JXXX,LR,JJ,i,im1,im2 !,je
!         integer jdofm(nmonths_linoz+1),jdmc(nmonths_linoz)
!      parameter (je=18) !number of latitudes in look-up table

      ! Now declare IM, JM as local variables
      ! since we have removed them from the common block (dbj 6/24/03)
      INTEGER IM, JM, MONTH

      real*8 STRTX(nlevels_linoz),YSTRT(nlat_linoz)
      real*8 P0L(llpar+1)
      real*8 STRT0L(llpar+1),STRT1L(llpar+1),STRT2L(llpar+1)
      real*8, PARAMETER :: PSF=1010D0

      !Define Month names locally      (dbj 6/25/03)
      CHARACTER(LEN=3)    :: CMONTH(12) = (/'jan', 'feb', 'mar', 'apr',
     &                                      'may', 'jun', 'jul', 'aug',
     &                                      'sep', 'oct', 'nov', 'dec'/)

!      data JDOFM/0,31,59,90,120,151,181,212,243,273,304,334,365/
c-----------------------------------------------------------------------

      ! Assign values for local IM and JM   (dbj 6/24/03)
      IM     = IIPAR
      JM     = JJPAR

      ! added call to GET_MONTH (dbj 6/25/03) 
      WRITE(6,*)'#####################################################'
      WRITE(6,*)'#    Interpolating Linoz fields for ',
     &                CMONTH( GET_MONTH() ),
     &                                                 '             #'
      WRITE(6,*)'#####################################################'


! ***** Linear interpolation between months is not currently used {PJC} *****
!c get weights for month interpolation
!      do i=1,nmonths_linoz
!        jdmc(i) = jdofm(i+1) - (jdofm(i+1)-jdofm(i))/2
!      enddo
!
!      im1=0
!      do i=1,nmonths_linoz
!        if (jdmc(i).lt.jday) then
!          im1=i
!        endif
!      enddo
!      if (im1.eq.0) then
!        im1=nmonths_linoz
!        im2=1
!        wm1=(jdmc(im2)-jday)*1.0/(jdmc(im2)-(jdmc(im1)-365.0))
!      elseif (im1.eq.nmonths_linoz) then
!        im2=1
!        wm1=(jdmc(im2)+365.0-jday)/(jdmc(im2)+365.0-jdmc(im1))
!      else
!        im2=im1+1
!        wm1=(jdmc(im2)-jday)*1.0/(jdmc(im2)-jdmc(im1))
!      endif
!      wm2=1.0-wm1
!
!c      write(6,*)iday,jday,' weights: ',wm1,wm2
!c      write(6,*)'months: ',im1,im2,month
!c      write(6,*)'between: ',jdmc(im1),jdmc(im2)
! ***************************************************************************

c latitude interpolation   

      YSTRT(1) = -85.d0   !Latitudes = -85, -75, -65, .... +75, +85.
      do J = 2,NLAT_LINOZ
        YSTRT(J) = YSTRT(J-1) + 10.d0
      enddo


      DO J = 1,JJPAR
        JXXX     = int(0.1d0 * GET_YMID(J) +10.d0)   ! (dbj 6/25/03)
        JLATMD(J) = MIN(18,MAX(1,JXXX)) !index of nearest Linoz data column
      ENDDO

      DO L = 1,LLPAR+1
         P0L(L) = GET_AP(LLPAR+2-L) + (GET_BP(LLPAR+2-L)*PSF)  ! dbj
      ENDDO

c-------- TPARM(25,18,12,N) defined for --------------------------------
c            25 layers from 58 km to 10 km by 2 km intervals
c            18 LATS (85S, 75S, ...85N) 
c            12 months
c            N tables = NTBLS
c-------- skip interpolating, pick nearest latitude --------------------

      DO N = 1,nfields_linoz

! ***** Interpolation between latitudes is not currently used {PJC} *****
!c----- interpolating along latitude, from TPAR2 to STRTXY
!        do K = 1,nlevels_linoz
!          do J = 1,nlat_linoz
!c            TPAR2(K,J) = TPARM(K,J,MONTH,N)
!            TPAR2(K,J) = TPARM(K,J,im1,N)
!          enddo
!        enddo
!        call LINOZ_INTPL(nlevels_linoz,NLAT_LINOZ,JPAR,JM,YSTRT,YDGRD,
!     &       TPAR2,STRTXY1)
!        do K = 1,nlevels_linoz
!          do J = 1,nlat_linoz
!            TPAR2(K,J) = TPARM(K,J,im2,N)
!          enddo
!        enddo
!        call LINOZ_INTPL(nlevels_linoz,NLAT_LINOZ,JPAR,JM,YSTRT,YDGRD,
!     &       TPAR2,STRTXY2)
! ***********************************************************************

         DO J = 1,JM
            JJ         = JLATMD(J)
            DO K = 1,nlevels_linoz
! linearly interpolate in latitude and month
!            STRTX(K) = STRTXY1(K,J)*wm1 + STRTXY2(K,J)*wm2
! linearly interpolate in latitude, single month
!            STRTX(K) = STRTXY2(K,J)
! nearest latitude, linearly interpolate in month
!            STRTX(K) = TPARM(K,JJ,im1,N)*wm1 + TPARM(K,JJ,im2,N)*wm2
! nearest latitude, single month
               STRTX(K) = TPARM(K,JJ,GET_MONTH(),N) ! (dbj 6/25/03)
            ENDDO   ! loop over K


            ! *PJC* Interpolate and calculate moments of column distribution 
            CALL LINOZ_STRT2M(STRTX,nlevels_linoz,STRT0L,STRT1L,STRT2L,
     &           P0L,LLPAR)

            ! Store loss freq/yields & moments in TLSTT/SWT/SWW 
            ! for exact CTM layers LM down
            ! Order reversed from C.McLinden version {PJC}
            DO LR = 1,LLPAR  
               TLSTT(J,LR,N) = STRT0L(LLPAR+1-LR)
            ENDDO
         ENDDO     ! loop over J
      ENDDO      ! loop over N

      END SUBROUTINE LINOZ_STRATL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: linoz_strt2m
!
! !DESCRIPTION: Subroutine LINOZ\_STRT2M sets up a std z* atmosphere: 
!  p = 1000 * 10**(-z*/16 km).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LINOZ_STRT2M(STRTX,NX,STRT0L,STRT1L,STRT2L,P0L,NSTRT)
!
! !USES:
!
#     include "F77_CMN_SIZE"
!
! !DEFINED PARAMETERS:
!
      ! Parameter (ncbox=25)  
      ! Now use nlevels_linoz for all routines. {PJC}
      INTEGER, PARAMETER   :: NL = NLEVELS_LINOZ+5
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)  :: NX
      INTEGER, INTENT(IN)  :: NSTRT
      REAL*8,  INTENT(IN)  :: STRTX(NLEVELS_LINOZ)
      REAL*8,  INTENT(IN)  :: P0L(LLPAR+1)
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: STRT0L(LLPAR+1)
      REAL*8,  INTENT(OUT) :: STRT1L(LLPAR+1)
      REAL*8,  INTENT(OUT) :: STRT2L(LLPAR+1)
!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer ncbox,l,k

      real*8  P1,P2,F0,F1,F2,PS(NL+1),F(NL)
      real*8  XPSD,XPSLM1,XPSL
c-----------------------------------------------------------------------
c         set up std z* atmosphere: p = 1000 * 10**(-z*/16 km)
c         assume that stratospheric chemical parameters always start at
cc         52 km (N=27) scan downward from 52 km to 14 km (NX=20) by 2 km
c         58 km (N=30) scan downward from 58 km to 10 km (NX=25) by 2 km      
c         intervals, constant >58km
c-------- N.B. F(@30km) assumed to be constant from 29-31 km (by mass) 
c
C======== Comments from Chris McLinden by Email ={PJC}==================
C     CALL SOMLFQ(P1,P2,F0,F1,F2,PS,F,NL)
C   - P1,P2 are the pressure EDGES for the CTM layer onto which the
C     coefficients will be mapped. [P1>P2 I believe {PJC}]
C   - F0,F1,F2 are the CTM layer vertical moments determined in SOMLFQ
C   - PS are the pressure layer edges of the original [ie Linox] grid
C   - F is the column of coefficients (on the original grid); note
C     F is flipped relative to STRTX and since the coefficients begin
C     at z*=10, F(1)=F(2)=...=F(5)=0
C   - NL is 30; size of F()
C      
C     The box model calculations were performed at z*=10km, 12km, ... and
C     so these would represent the centres with the corresponding edges at
C     9,11km ; 11,13km; ...
C     PS() represents the edges (although PS(1) is set to 1000mb).
C     The first few values are:
C     PS(1)=1000
C     PS(2)=874.947105    (note PS(2) is not quite 1000 exp(-1/16) as the
C     PS(3)=656.117767     the average pressure is used - not the pressure
C     PS(4)=492.018914     at the average z*)
C     PS(5)=368.96213
C     PS(6)=276.68257
C     PS(7)=207.48266
C     ...
C     PS(30)=0.276682568
C     PS(31)=0.0
C     
C     F(1) spans PS(1)-PS(2)
C     F(2) spans PS(2)-PS(3)
C     ...
C     F(30) spans PS(30)-PS(31)
C=======================================================================


      XPSD       = 10.D0 **(-0.125D0)
      XPSLM1     = 1000.D0
      PS(1)      = 1000.D0
      DO L = 2,NL
        XPSL     = XPSLM1 *XPSD
        PS(L)    = 0.5D0 *(XPSLM1 +XPSL)
        XPSLM1   = XPSL
      ENDDO
      PS(NL+1)   = 0.D0
      DO L = 1,NL-NX
        F(L)     = 0.D0
      ENDDO
c-------- K=1 is at the top of atmosphere ------------------------------
      DO K = 1,NX
        F(NL+1-K)= STRTX(K) !STRTX has increasing preasure. {PJC}
      ENDDO
      DO K = 1,NSTRT
        P1       = P0L(K+1)
        P2       = P0L(K)
        CALL LINOZ_SOMLFQ(P1,P2,F0,F1,F2,PS,F,NL)
        STRT0L(K)= F0
        STRT1L(K)= F1
        STRT2L(K)= F2
      ENDDO

      END SUBROUTINE LINOZ_STRT2M
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: linoz_somlfq
!
! !DESCRIPTION: subroutine LINOZ\_SOMLFQ calculates loss freq moments from a 
!  set of loss frequencies at std z*, given a CTM model interval pressure
!  range: P1 > P2 (decreasing up)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LINOZ_SOMLFQ(P1,P2,F0,F1,F2,PS,F,NL)
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)  :: NL
      REAL*8,  INTENT(IN)  :: F(NL)
      REAL*8,  INTENT(IN)  :: PS(NL+1)
      REAL*8,  INTENT(OUT) :: P1
      REAL*8,  INTENT(OUT) :: P2
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: F0
      REAL*8,  INTENT(OUT) :: F1
      REAL*8,  INTENT(OUT) :: F2
!
! REMARKS:
! The pressure levels BETWEEN z* values are:
!      PS(i) > PS(i+1) bounds z*(i)
!                                                                             .
! NL:  z* levels, ==> PS(NL+1) = 0  (extrapolate chemical loss to top)
!      Z1 = 16.D0*LOG10(1000.D0/P1)
!      Z2 = 16.D0*LOG10(1000.D0/P2)
!                                                                             .
! The MOMENTS for a square-wave or 'bar': F(x)=f0  b<=x<=c, =0.0 else
!      S0 =   f0 (x)                      [from x=b to x=c]
!      S1 = 3 f0 (x^2 - x)                [from x=b to x=c]
!      S2 = 5 f0 (2x^3 - 3x^2 + x)        [from x=b to x=c]
!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer  I
      real*8   XB,XC,PC,PB,THIRD,sgnf0

      F0     = 0.D0
      F1     = 0.D0
      F2     = 0.D0
      DO I = 1,NL
         PC   = MIN(P1,PS(I))
         PB   = MAX(P2,PS(I+1))
         IF (PC .GT. PB)  THEN

            ! have condition:  P1>=PC > PB>=P2, 0<=XB < XC<=1
            XC = (PC-P2)/(P1-P2)
            XB = (PB-P2)/(P1-P2)

            ! assume that the loss freq, F, is constant over interval [XB,XC],
            ! F0: (c-b),   
            ! F1: 6((c2-c)-(b2-b)),  
            ! F2: 5((2c3-3c2+c)-(2b3-3b2+b))
            ! calculate its contribution to the moments in the interval [0,1]
            F0 = F0 +F(I) *(XC -XB)
            F1 = F1 +F(I) *3.D0 *((XC *XC -XC) - (XB *XB -XB))
            F2 = F2 +F(I) *5.D0 *
     &           ((XC+XC-1.D0)*(XC*XC -XC) - (XB+XB-1.D0)*(XB*XB -XB))
         ENDIF
      ENDDO

      ! RESTRAIN moments: force monotonicity & positive at min end pt

      ! cam: tables can be + or -
      if (f0.ne.0.0) then
         sgnf0=f0 / abs(f0)
      else
         sgnf0=1.0
      endif
      f0=abs(f0)
      
      !F0 = MAX(F0, 0.D0)
      THIRD = 1.D0/3.D0
      IF (F2 .GT. 0.D0)  THEN


         ! do not allow reversal of curvature: F2 > 0
         F2   = MIN(F2, ABS(F1)*THIRD, 5.D-1*F0)
         IF (F1 .LT .0.D0)  THEN
            F1 = MAX(-(F0+F2), F1)
         ELSE
            F1 = MIN(+(F0+F2), F1)
         ENDIF
      ELSE

         ! F2 < 0 = curved down at ends, allow if F1 < F0
         F1  = MIN(F0,MAX(-F0,F1))
         F2  = MAX(F2,(ABS(F1)-F0),(-ABS(F1)*THIRD))
      ENDIF

      ! cam: apply sign
      f0=sgnf0 * f0
      f1=sgnf0 * f1
      f2=sgnf0 * f2

      END SUBROUTINE LINOZ_SOMLFQ
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: linoz_read
!
! !DESCRIPTION: Subroutine LINOZ\_READ reads the input data file for the
!  Linoz stratospheric ozone chemistry.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LINOZ_READ
!
! !USES:
!
      USE FILE_MOD,      ONLY : IU_FILE        ! Logical unit #
      USE FILE_MOD,      ONLY : IOERROR        ! I/O error subroutine
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1   ! Data directory path

#     include "F77_CMN_SIZE"
!
! !REMARKS:
!  LINOZ_READ is called from "main.f" at the start of the simulation.
!  LINOZ_READ will also call INIT_LINOZ to initialize the arrays.
!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!  16 Oct 2009 - R. Yantosca         - Now use IU_FILE instead of IU_LINOZ
!  16 Oct 2009 - R. Yantosca         - Read file from DATA_DIR_1x1
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      LOGICAL, SAVE      :: FIRST = .TRUE.
      INTEGER            :: K,       J,   M, NTBLS, IOS
      REAL*8             :: TMAX,    TMIN
      CHARACTER(LEN=80)  :: HEADING, TITL1
      CHARACTER(LEN=255) :: FILENAME

      ! Only initialize arrays on first timestep
      IF ( FIRST ) THEN 
         CALL INIT_LINOZ
         FIRST = .FALSE.
      ENDIF

      ! Filename
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'Linoz_200910/Linoz_March2007.dat'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - LINOZ_READ: Reading ', a )

      ! new std z*=2km levels from model:  z*=10,12,...(25*2)+8 km
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', 
     &               FORM='FORMATTED',      IOSTAT=IOS )

      ! 
      IF ( IOS /= 0 ) THEN
         CALL IOERROR( IOS, IU_FILE, 'read_linoz_coeff_file' )
      ENDIF

      ! Reade header
      READ ( IU_FILE, 901 ) HEADING
      WRITE(6,*)            TRIM( HEADING )

      ! Loop thru file
      DO NTBLS = 1,nfields_linoz
         TMIN = +1.d30
         TMAX = -1.d30
         READ (IU_FILE,901) TITL1
         do M = 1,nmonths_linoz !Month
            do J = 1,nlat_linoz !Latitudes
               READ (IU_FILE,907)  
     &               (TPARM(K,J,M,NTBLS),K=nlevels_linoz,1,-1)
               do K=1,nlevels_linoz
                  TMAX = max (TMAX, TPARM(K,J,M,ntbls))
                  TMIN = min (TMIN, TPARM(K,J,M,ntbls))
               enddo
            enddo
         enddo
         write (6,912) TITL1,TMIN,TMAX
      enddo

      WRITE(6,*)'$$ Finished Reading Linoz Data $$'
      WRITE(6,*)

      GOTO 999

      ! If error has occurred
 101  CONTINUE
      WRITE(6,*)'**** STOP: Error reading Linoz Coefficients {PJC} ****'
      write(6,*)'NTBLS =',ntbls,',  M =',m,',  J =',j,',  K =',k
      write(6,*)'TPARM(K,J,M,NTBLS) =',TPARM(K,J,M,NTBLS)
      STOP

      ! Format strings
 901  FORMAT(A)
 907  FORMAT(20X,6E11.4/(8E11.4))
c907  FORMAT(20X,6E10.3/(8E10.3))
 912  FORMAT('  Linoz Data:  ',a80,1p,2e10.3)
      
 999  CONTINUE

      ! Close the files
      CLOSE( IU_FILE ) 

      END SUBROUTINE LINOZ_READ
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: linoz_intpl
!
! !DESCRIPTION: Subroutine LINOZ\_INTPL does some kind of interpolation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LINOZ_INTPL(KE,IE,ND,NE,XI,XN,YI,YN)
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)  :: KE
      INTEGER, INTENT(IN)  :: IE
      INTEGER, INTENT(IN)  :: ND
      INTEGER, INTENT(IN)  :: NE
      REAL*8,  INTENT(IN)  :: XI(IE) 
      REAL*8,  INTENT(IN)  :: XN(ND)
      REAL*8,  INTENT(IN)  :: YI(KE,IE)
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: YN(KE,ND)

!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer  I,II,J,K
      real*8   CNST1,CNST2

      ! k=height; i=lat
      J       = 2
      do I = 1,NE
        if (XN(I) .gt. XI(1        ))  then
          if (XN(I) .lt. XI(IE))  then
            CNST1     = (XI(J) - XN(I)) / (XI(J) - XI(J-1))
            CNST2     = (XN(I) - XI(J-1)) / (XI(J) - XI(J-1))
            do K = 1,KE
              YN(K,I) = CNST1 * YI(K,J-1) + CNST2 * YI(K,J)
            enddo
            II    = min(I+1,NE)
            if (XN(II) .gt. XI(J))  J = min(IE,J+1)
          else
            do K = 1 ,KE
              YN(K,I) = YI(K,IE)
            enddo
          endif
        else
          do K = 1,KE
            YN(K,I)   = YI(K,1)
          enddo
        endif
        !write(6,*)i,(yn(k,i),k=1,ke)
      enddo

      END SUBROUTINE LINOZ_INTPL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: strat_init
!
! !DESCRIPTION: Subroutine STRAT\_INIT copies the ozone computed by the
!  Linoz stratospheric chemistry algorithm back into the GEOS-Chem
!  tracer array.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE STRAT_INIT
!
! !USES:
!
      USE TRACERID_MOD 
      USE TRACER_MOD   

#     include "F77_CMN_SIZE"
#     include "CMN"
!
! !REVISION HISTORY: 
!  24 Jun 2003 - B. Field & D. Jones - Further updates for GEOS-Chem
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER I, J, L

      CALL LINOZ_STRATL

      DO J = 1,JJPAR
      DO I = 1,IIPAR
      DO L = MINVAL(LPAUSE),LLPAR
         IF (L .LT. LPAUSE(I,J)) CYCLE
         STT(I,J,L,IDTOX) = TLSTT(J,L,1) / TCVV(IDTOX)
      ENDDO
      ENDDO
      ENDDO

!      call write_fields2(7)
!      call flush(12)

      END SUBROUTINE STRAT_INIT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_linoz
!
! !DESCRIPTION: Subroutine INIT\_LINOZ allocates and zeroes the module arrays 
!  used in the Linoz stratospheric ozone algorithm.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_LINOZ
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "F77_CMN_SIZE" 
!
! !REVISION HISTORY: 
!  16 Oct 2009 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: AS

      ! For safety's sake, only allocate arrays on first call
      IF ( FIRST ) THEN 

         ! Allocate TPARM array
         ALLOCATE( TPARM( NLEVELS_LINOZ, NLAT_LINOZ,
     &                    NMONTHS_LINOZ, NFIELDS_LINOZ ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TPARM' )
         TPARM = 0d0

         ! Allocate TLSTT array
         ALLOCATE( TLSTT( JJPAR, LLPAR, NFIELDS_LINOZ ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TPARM' )
         TLSTT = 0d0

         ! Reset FIRST
         FIRST = .FALSE.
      ENDIF

      END SUBROUTINE INIT_LINOZ
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_linoz
!
! !DESCRIPTION: Subroutine CLEANUP\_LINOZ deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_LINOZ
!
! !REVISION HISTORY: 
!  16 Oct 2009 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Deallocate arrays
      IF ( ALLOCATED( TPARM ) ) DEALLOCATE( TPARM )
      IF ( ALLOCATED( TLSTT ) ) DEALLOCATE( TLSTT )
      
      END SUBROUTINE CLEANUP_LINOZ
!EOC

      ! End of module
      END MODULE LINOZ_MOD
