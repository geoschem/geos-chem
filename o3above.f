C $Id: o3above.f,v 1.1 2003/06/30 20:26:05 bmy Exp $
C*****************************************************************************
      SUBROUTINE o3above(NSEASON,NCLIMATOLOGY,NCLIMATOLOGY2)
C*****************************************************************************
      
      ! References to F90 modules
      USE DAO_MOD, ONLY : BXHEIGHT

      IMPLICIT NONE

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! LPAUSE
#     include "CMN_SETUP"  ! REMOVE_CMD, SPACE, ZIP_SUFFIX, etc.
#     include "CMN_CO"     ! CO arrays
!----------------------------------------------------------------------------
! Prior to 9/18/02:
!#     include "CMN_NOX"    ! BXHEIGHT
!----------------------------------------------------------------------------

      INTEGER NSEASON,NSEAS,NLAYO32,MVERT,ML,NLATS,JJ
      PARAMETER (NSEAS=4,NLATS=8,NLAYO32=9,NPROF=40)
      INTEGER MLAT,MLON,M,K,L,J,I,MM,LATDEG,NPROF,NDO1
      REAL*8 AALT,AA,HENHT,AB,AC,Pjal,Pdao,RAT,P100,AP
      REAL*8 COL,R,DENSO3(IIPAR,JJPAR,LLPAR)
      REAL*8 O3TROP(IIPAR,JJPAR,NLAYO32),Pclim(NLAYO32)
      REAL*8 O3TROPC(IIPAR,JJPAR,LLPAR),
     *     O3TROPD(IIPAR,JJPAR,LLPAR),sum,
     *     SO3(NPROF),SALT(NPROF),sumht,SPRESS(NPROF),
     *     delZ(IIPAR,JJPAR,LLPAR),totZ(IIPAR,JJPAR,LLPAR)
      REAL*8 BAND1(NLATS),BAND2(NLATS)
      CHARACTER ( LEN=6 )   :: MERGE
      CHARACTER ( LEN=18 )  :: MERGE5
      CHARACTER ( LEN=255 ) :: O3COL_DIR,CHARO3,CHARRM,CHART
      CHARACTER ( LEN=5 )   :: TEMPO
      CHARACTER ( LEN=5 )   :: TEMPO2
      CHARACTER ( LEN=3 )   :: XSEASONS(NSEAS)
      CHARACTER ( LEN=3 )   :: XSEASONS2(NSEAS)
      CHARACTER ( LEN=6 )   :: XLATS(NLATS)
      DATA XSEASONS  /'win','spr','sum','aut'/
      DATA XSEASONS2 /'WIN','SPR','SUM','AUT'/
      DATA Pclim /1000.,900.,800.,700.,500.,300.,200.,150.,100./
      DATA XLATS /'75S65S','55S45S','00035S','25S05S',
     &     '05N25N','00035N','45N55N','65N75N'/
      DATA BAND1/-90.,-60.,-40.,-30.,0.,30.,40.,60./
      DATA BAND2/-60.,-40.,-30.,0.,30.,40.,60.,90./
      
C*****************************************************************************
C*****************************************************************************
C This SR calculates the O3 COLumn (DU) above a certain
C     grid box.
C*****************************************************************************
C*****************************************************************************
C The annual mean tropopause is stored in the LPAUSE array 
C (from header file "CMN").  LPAUSE is defined such that: 
C
C Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
C         LPAUSE(I,J) <= L <= LLPAR           are stratospheric
C
C We now use LPAUSE instead of NSKIPL to denote the strat/trop boundary. 
C (bmy, 4/18/00)  
C*****************************************************************************
C
      print*,'Calculate the overhead O3 column in SR o3above.'
C
C o3up = O3 column (DU) above
C
      o3up(:,:,:)=0.
C
C*****************************************************************************
C Loop over season.
C
      DO MM=NSEASON,NSEASON
C
C*****************************************************************************
C First, calculate the tropospheric portion (i.e., > 100 mb)
C   of the column using jal's tropospheric O3 profiles.
C*****************************************************************************
C
C A) Read in O3 (ppbv) from the O3/Temperature climatologies.
C    Store these values in O3TROP(:,:,:).
C    The pressure levels correpsonding to O3TROP are
C    stored in Pclim().
C
C      
cbnd         O3COL_DIR = '/users/ctm/bnd/o3clim/'
C
cbnd         TEMPO='tempo'
C
C Write appropriate profile name into MERGE.
C
         WRITE(MERGE,2)XSEASONS(MM)
 2       FORMAT('o3.',A3)
C
C unzip and copy file to temporary directory.
C
cbnd         CHARO3 =
cbnd     *        TRIM(UNZIP_CMD)//SPACE(1:1)//
cbnd     *        TRIM(O3COL_DIR)//TRIM(MERGE)// !TRIM(ZIP_SUFFIX)//
cbnd     *        TRIM(REDIRECT)//SPACE(1:1)//
cbnd     *        TRIM(TEMP_DIR)//TEMPO
C
cbnd         CALL SYSTEM(TRIM(CHARO3))
C
         CHART = TRIM(DATA_DIR)//MERGE
C
         CLOSE(NO3READ)
         OPEN(NO3READ,FILE=CHART,STATUS='OLD')
         REWIND(NO3READ)
C     
         DO I=1,IIPAR
         DO J=1,JJPAR
         DO L=1,NLAYO32
            READ(NO3READ,*)O3TROP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
C
         CLOSE(NO3READ)
C     
C B) Interpolate using ln(pressure) to put O3 values
C    on DAO model sigma layers.  Store the interpolated
C    climatological O3 values in O3TROPC(:,:,:).
C
         P100=DLOG(Pclim(NLAYO32))
         O3TROPC(:,:,:)=0.

         DO I=1,IIPAR
         DO J=1,JJPAR
         ! M now goes up to the annual mean tropopause (bmy, 4/18/00)
         !DO M=1,NSKIPL-1
         DO M = 1, LPAUSE(I,J) - 1

            Pdao=0.
            Pdao=DLOG(Pavg(I,J,M))
         
            DO L=1,NLAYO32

               Pjal=0.
               Pjal=DLOG(Pclim(L))
C
C  If > 100. then fill in with stratospheric info.
C
               IF(Pdao.LT.P100) THEN
                  O3TROPC(I,J,M)=0.
                  GOTO 299
               ENDIF

               IF(L.EQ.1.AND.Pjal.LT.Pdao) THEN
                  O3TROPC(I,J,M)=O3TROP(I,J,L)
                  GOTO 299
               ENDIF

               IF(Pjal.EQ.Pdao) THEN
                  O3TROPC(I,J,M)=O3TROP(I,J,L)
                  GOTO 299
               ENDIF

               IF(Pjal.GT.Pdao) GOTO 298
               
               IF(Pjal.LT.Pdao) THEN
                  RAT=(DLOG(Pclim(L-1))-Pdao)/(DLOG(Pclim(L-1))-
     *                 Pjal)
                  O3TROPC(I,J,M)=O3TROP(I,J,L-1)-RAT*
     *                 (O3TROP(I,J,L-1)-O3TROP(I,J,L))
                  GOTO 299
               ENDIF
               
               PRINT*,'Problem in SR O3above!'
               PRINT*,'Grid interpolation failed!'
               STOP
           
 298           CONTINUE
            
            ENDDO
        
 299        CONTINUE
        
         ENDDO
         ENDDO
         ENDDO
C
C C) Convert O3(ppbv) to O3(molec/cc).
C
         R=0.08204*1000.*1000.
         DENSO3(:,:,:)=0.
C
         ! Make sure DO loops are ordered L-J-I (bmy, 4/16/00)
         DO L=1,LLPAR
         DO J=1,JJPAR
         DO I=1,IIPAR
            DENSO3(I,J,L)=Pavg(I,J,L)*AVGNO/(R*Tavg(I,J,L))
            O3TROPD(I,J,L)=O3TROPC(I,J,L)*1.E-9
            O3TROPC(I,J,L)=O3TROPC(I,J,L)*1.E-9*DENSO3(I,J,L)
         ENDDO
         ENDDO
         ENDDO
C
C If NCLIMATOLOGY = 1 then fill O3 fields with 3-D climatology.
C
         IF(NCLIMATOLOGY.EQ.1) THEN

            PRINT*,'In SR O3ABOVE.f'
            PRINT*,'O3TROPD only goes to NSKIPL-1 so may not cover'
            PRINT*,'entire height of troposphere!'
            PRINT*,'Need to fix this!'
            STOP

            ! Make sure DO-loops are ordered L-J-I (bmy, 4/16/00) 
            ! L goes from 1 up to the annual mean tropopause (bmy, 4/18/00)
            !DO L=1,NSKIPL-1
            DO L=1,MAXVAL( LPAUSE )
            DO J=1,JJPAR
            DO I=1,IIPAR

               ! Only process tropospheric boxes (bmy, 4/18/00)
               IF ( L < LPAUSE(I,J) ) THEN
                  BBIJ(I,J,L,8)=O3TROPD(I,J,L)
               ENDIF
            ENDDO
            ENDDO
            ENDDO

         ENDIF
C
C D) Integrate "tropospheric part" of column above.
C    
         delZ(:,:,:)=0.
         totZ(:,:,:)=0.
C
C ****************************************************************************
         DO MLON=1,IIPAR
C ****************************************************************************
C
            DO MLAT=1,JJPAR
C
               sumht=0.
C
               ! ML goes up to the annual mean tropopause (bmy, 4/18/00)
               !DO ML=1,NSKIPL-1
               DO ML = 1, LPAUSE(MLON,MLAT) - 1
C
                  delZ(MLON,MLAT,ML)=BXHEIGHT(MLON,MLAT,ML)/1000.
                  sumht=sumht+delZ(MLON,MLAT,ML)
                  totZ(MLON,MLAT,ML)=sumht
C     
               ENDDO
C
C ****************************************************************************
               ! MVERT goes up to the annual mean tropopause (bmy, 4/18/00)
               !DO MVERT=1,NSKIPL-1
               DO MVERT = 1, LPAUSE(MLON,MLAT) - 1
C ****************************************************************************
C
                  AALT=totZ(MLON,MLAT,MVERT)
C
                  COL=0.
                  NDO1=0
C
                  ! J goes up to the annual mean tropopause (bmy, 4/18/00)
                  !DO 87 J=2,NSKIPL-1
                  DO 87 J = 2, LPAUSE(MLON,MLAT) - 1
                     IF(AALT.LT.totZ(MLON,MLAT,J)) GOTO 80
 87               CONTINUE
C
                  COL=0.
                  NDO1=1

 80               L=J-1
C
                  IF(NDO1.EQ.0) THEN
                     COL=(totZ(MLON,MLAT,L)-AALT)*
     &                    ((O3TROPC(MLON,MLAT,L+1)+
     *                    O3TROPC(MLON,MLAT,L))/2.)*1.E05
C     
                     ! I goes up to the annual mean tropopause (bmy, 4/18/00)
                     !DO 81 I=L,NSKIPL-1-1
                     DO 81 I = 1, LPAUSE(MLON,MLAT) - 1
                        HENHT=(totZ(MLON,MLAT,I+1)-
     &                       totZ(MLON,MLAT,I))*1.E05
                        COL=COL+HENHT*((O3TROPC(MLON,MLAT,I)+
     *                       O3TROPC(MLON,MLAT,I+1))/2.)
 81                  CONTINUE
C
                  ENDIF
C
                  COL=COL/2.69E16
C     
                  o3up(MLON,MLAT,MVERT)=COL
C
C ****************************************************************************
               ENDDO
C ****************************************************************************
            ENDDO
C ****************************************************************************
         ENDDO
C ****************************************************************************
C
C remove file from temporary directory.
C
         CHARRM =
     *        TRIM(REMOVE_CMD)//SPACE(1:1)//
     *        TRIM(TEMP_DIR)//TEMPO
C     
         CALL SYSTEM(TRIM(CHARRM))
C
C*****************************************************************************
C Second, fill the boxes with new.* O3/T profiles. Only
C     the stratospheric portion (i.e., < P100) of these 
C     files will be used here. The files are 2-d (e.g.
C     30-40N zonally averaged latitude band).
C*****************************************************************************
C
cbnd         O3COL_DIR='/users/ctm/bnd/chem1d/input/innaatms/trop/'
         TEMPO2='tempo2'
C
         DO K=1,NLATS
C
C Write appropriate profile name into MERGE.
C
            WRITE(MERGE5,22)XSEASONS2(MM),XLATS(K)
 22         FORMAT('new.atm.',A3,'.',A6)
C
C unzip and copy file to temporary directory.
C
            CHARO3 =
     *           TRIM(UNZIP_CMD)//SPACE(1:1)//
cbnd     *           TRIM(O3COL_DIR)//TRIM(MERGE5)//TRIM(ZIP_SUFFIX)//
     *           TRIM(MERGE5)//TRIM(ZIP_SUFFIX)//
     *           TRIM(REDIRECT)//SPACE(1:1)//
     *           TRIM(TEMP_DIR)//TEMPO2
C
            CALL SYSTEM(TRIM(CHARO3))
C     
            CHART = TRIM(TEMP_DIR)//TEMPO2
C     
            CLOSE(NO3READ)
            OPEN(NO3READ,FILE=CHART,STATUS='UNKNOWN')
            REWIND(NO3READ)
C
            READ(NO3READ,903)
 903        FORMAT(////)
C
            DO 89 J=1,NPROF
               READ(NO3READ,*)AA,SALT(J),AB,SO3(J),SPRESS(J)
 89         CONTINUE
C
            CLOSE(NO3READ)
C
            LATDEG=-90.
C
C ****************************************************************************
            DO MLAT=1,JJPAR
C ****************************************************************************
C
               IF(MLAT.EQ.1.OR.MLAT.EQ.46) THEN
                  LATDEG=LATDEG+2.
               ELSE
                  LATDEG=LATDEG+4.
               ENDIF
C
C 
C ****************************************************************************
               IF(LATDEG.GT.BAND1(K).AND.LATDEG.LE.BAND2(K)) THEN
C ****************************************************************************
C
C ****************************************************************************
                  DO MLON=1,IIPAR
C ****************************************************************************
C
C Add stratospheric portion of O3 column to tropopspheric portion.
C
                     ! Now use LPAUSE instead of NSKIPL (bmy, 4/18/00)
                     !AALT=totZ(MLON,MLAT,NSKIPL-1)
                     AALT=totZ(MLON,MLAT,LPAUSE(MLON,MLAT)-1)
C
                     COL=0.
C
                     DO 870 J=2,NPROF
                        IF(AALT.LT.SALT(J)) GOTO 800
 870                 CONTINUE

 800                 L=J-1

                     COL=(SALT(L)-AALT)*((SO3(L+1)+SO3(L))/2.)*1.E05
                     DO 810 I=L,NPROF-1
                        HENHT=(SALT(I+1)-SALT(I))*1.E05
                        COL=COL+HENHT*((SO3(I)+SO3(I+1))/2.)
 810                 CONTINUE
C
                     COL=COL/2.69E16
C
C ****************************************************************************
                     ! MVERT goes up to the annual mean tropopause 
                     ! (bmy, 4/18/00)
                     !DO MVERT=1,NSKIPL-1
                     DO MVERT = 1, LPAUSE(MLON,MLAT) - 1
C ****************************************************************************
C
                        o3up(MLON,MLAT,MVERT)=o3up(MLON,MLAT,MVERT)+COL
C
C ****************************************************************************
                     ENDDO
C ****************************************************************************
C
C ****************************************************************************
                     ! MVERT goes from the annual mean tropopause up to
                     ! level LLPAR-1 (bmy, 4/18/00)
                     !DO MVERT=NSKIPL,LLPAR-1
                     DO MVERT = LPAUSE(MLON,MLAT), LLPAR-1
C ****************************************************************************
C 
                        Pdao=0.
                        Pdao=Pavg(MLON,MLAT,MVERT)
C     
C ****************************************************************************
C  Above NSKIPL-1 fill in with stratospheric info
C ****************************************************************************
                        DO JJ = 1,NPROF
C ****************************************************************************
                           IF(SPRESS(JJ).LT.Pdao) THEN
C     
                              AALT=SALT(JJ)
C
                              COL=0.
C
                              DO 8700 J=JJ-2,NPROF
                                 IF(AALT.LT.SALT(J)) GOTO 8000
 8700                         CONTINUE

 8000                         L=J-1
C
                              COL=(SALT(L)-AALT)*((SO3(L+1)+SO3(L))/2.)*
     &                             1.E05
                              DO 8100 I=L,NPROF-1
                                 HENHT=(SALT(I+1)-SALT(I))*1.E05
                                 COL=COL+HENHT*((SO3(I)+SO3(I+1))/2.)
 8100                         CONTINUE
C     
                              COL=COL/2.69E16
C     
                              GOTO 388
C
                           ENDIF ! IF(SPRESS(JJ).LT.Pdao) THEN

C ****************************************************************************
                        ENDDO   ! DO JJ = 1,NPROF 
C ****************************************************************************
C
 388                    CONTINUE
C
                        o3up(MLON,MLAT,MVERT)=COL
C     
C ****************************************************************************
                     ENDDO ! DO MVERT=NSKIPL,LLPAR-1
C ****************************************************************************
                  ENDDO    ! DO MLON=1,IIPAR
C ****************************************************************************
               ENDIF       ! IF(LATDEG.GT.BAND1(K).AND.LATDEG.LE.BAND2(K)) THEN
C ****************************************************************************
            ENDDO          ! DO MLAT=1,JJPAR
C ****************************************************************************
         ENDDO             ! DO K=1,NLATS
C ****************************************************************************

         CHARRM =
     *        TRIM(REMOVE_CMD)//SPACE(1:1)//
     *        TRIM(TEMP_DIR)//TEMPO2
C     
         CALL SYSTEM(TRIM(CHARRM))
C
C*****************************************************************************
C End loop over season.
C
      ENDDO
C*****************************************************************************

      RETURN
      END
