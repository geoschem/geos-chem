! $Id: diag48.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      SUBROUTINE DIAG48 
!
!******************************************************************************
!  Subroutine DIAG48 does time series diagnostics.
!  (bmy, bey, amf, 6/1/99, 2/11/03)
!------------------------------------------------------------------------------
!  **** ND48 must be 1 in input.ctm 
!  Information about stations are given in inptr.ctm
!  Example: 
!      I    J    L  MS   NS
!     53   27    5   1   10     
!      
!  I, J, L, MS, NS are read as '5I5' in the subroutine inptr.f
!
!  if MS   = 0 -> Store a time series of the tracer NS ( in [v/v] )
!                 for the box (I,J) from the level 1 to L  
!
!  if MS   = 1 -> Store a time series of the quantity defined by NS
!                 for the box (I,J) from the level 1 to L  
!
!                 NS =  1 : O3   concentration    [v/v]
!                 NS =  2 : OH   concentration    [molec/cm3]
!                 NS =  3 : NOY  concentration    [v/v]
!                 NS =  4 : NO2  drydep velocity  [cm/s]
!                 NS =  5 : O3   drydep velocity  [cm/s]
!                 NS =  6 : PAN  drydep velocity  [cm/s]    
!                 NS =  7 : HNO3 drydep velocity  [cm/s]   
!                 NS =  8 : H2O2 drydep velocity  [cm/s]   
!                 NS =  9 : NO   concentration    [v/v
!                 NS = 10 : Boundary layer height [mb]
!                 NS = 11 : Local Time at station [hours]
!                 NS = 12 : ISOP emission flux    [molec/cm2/s]
!                 NS = 13 : Monin-Obhukov Length  [m]
!                 NS = 14 : NO3  concentration    [v/v]
!                 NS = 15 : Surface pressure      [mb]
!                 NS = 16 : Water vapor           [v/v]
!                 NS = 17 : Relative Humidity     [%]
!                 NS = 18 : HCN Column density    [molec/cm2]
!                 NS = 19 : Temperature           [K]
!
!------------------------------------------------------------------------------
!  Time series info is stored in TCOBOX(MAXACC,MAXSTA), which is now
!  an allocatable array in "diag_mod.f". (bmy, 4/24/0)
! 
!  The maximum number accumulating timesteps is set in CMN_DIAG
!  and is currently set to 1200. (bmy, 4/24/00)
!  .
!------------------------------------------------------------------------------
!  NOTES:
!  (1 ) KDA48 is now incremented in "main.f" (bmy, 11/5/99)
!  (2 ) TCOBOX is now declared allocatable in "diag_mod.f". (bmy, 11/29/99)
!  (3 ) It is advised to store local time in the range 0-24 hours to 
!        avoid precision problems when reading the "ctm.ts" file in IDL. 
!        (bmy, 2/4/00)
!  (4 ) Now get NO3 concentrations from array SAVENO3.  Also updated comments.
!        (bmy, 4/20/00)
!  (5 ) Now reference AVGW, OBK, PBL, PS from "dao_mod.f".  Removed PBL, PS,
!        OBK as arguments.  Make sure AVGW and OBK are allocated before
!        referencing them (bmy, 9/25/01)
!  (6 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (7 ) Now use P(I,J)+PTOP instead of PS(I,J), since that is the way to
!        make sure the P is consistent with the airmass AD. (bmy, 4/11/02)
!  (8 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE.  Also removed
!        obsolete, commented-out code. (bmy, 6/25/02)
!  (9 ) Replaced references to P(I,J) with call to GET_PEDGE(I,J,1) 
!        from "pressure_mod.f".  Also eliminated obsolete, commented-out
!        code. (dsa, bdf, bmy, 8/20/02)
!  (10) Now reference AD, T from "dao_mod.f".  Now call GEOS_CHEM_STOP to
!        free memory before stopping the run.  Also references IDTNOX, IDTOX,
!        etc. from "tracerid_mod.f".  Remove reference to OBK, since that is
!        now a local variable in "drydep_mod.f". (bmy, 11/20/02)
!  (11) Now replace DXYP(JREF) with routine GET_AREA_M2 of "grid_mod.f"
!        Now use function GET_LOCALTIME from "time_mod.f".  Now removed
!        reference to obsolete STH and RH2 arrays. (bmy, 2/11/03)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD, AVGW, PBL, T
      USE DIAG_MOD,     ONLY : TCOBOX
      USE ERROR_MOD,    ONLY : GEOS_CHEM_STOP
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : GET_LOCALTIME
      USE TRACERID_MOD 
      
      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "CMN"
#     include "CMN_DIAG"
#     include "CMN_O3"
#     include "CMN_DEP"

      ! Local variables
      INTEGER            :: I, J, L, N, K, MS, IPL

      REAL*8             :: FDTT, TCBXXX
      REAL*8             :: XLOCTM, TOPBL
      REAL*8,  SAVE      :: XXLOCTM(NNSTA) 

      LOGICAL, SAVE      :: FIRST(NNSTA)
      LOGICAL, SAVE      :: FFIRST = .TRUE.
      LOGICAL            :: AVGW_ALLOCATED  !, OBK_ALLOCATED (bmy, 11/20/02)
!
!*****************************************************************************
!  Initialize the FIRST array.  This is used for computing the local time
!  at each station.
!*****************************************************************************
!
      IF ( FFIRST ) THEN 
         DO K = 1, MAXSTA 
            FIRST(K) = .TRUE.
         ENDDO
         FFIRST = .FALSE. 
      ENDIF

      ! Also test if AVGW and OBK are allocated (bmy, 9/25/01)
      AVGW_ALLOCATED = ALLOCATED( AVGW )
!
!*****************************************************************************
!  STATISTICS FOR TIME SERIES: 
!
!  Get I, J, L, N, MS from arrays read in from the 'inptr.ctm' file 
!  for each station.  Exit loop and return I, J, L are less than 1.
!*****************************************************************************
!
      IF ( ND48 > 0 .and. KDA48 <= MAXACC ) THEN

         DO K = 1, MAXSTA 
            I = ISAVTC(K)             
            J = JSAVTC(K)  
            L = LSAVTC(K) 
            N = NSAVTC(K)
            MS= MSAVTC(K)

            IF ( I <= 0 .or. J <= 0 .or. L <= 0 ) EXIT

            ! Save air mass in FDTT
            FDTT = AD(I,J,L)
!
!*****************************************************************************
!  If MS = 0 then store STT time series for tracer (V/V) at each station
!*****************************************************************************
!
            IF ( MS == 0 ) THEN
               TCOBOX(KDA48,K) = TCVV(N) * STT(I,J,L,N) / FDTT
!
!*****************************************************************************
!  If MS = 1, then store other time series quantities (as listed above)
!*****************************************************************************!
!
            ELSE IF ( MS == 1 ) THEN 

               !==============================================================
               ! MS = 1, Tracer 1:  Store O3 time series [v/v] 
               !
               ! Multiply Ox (use IDTOX) by FO3(5) to get O3 concentration.
               ! Only do this for an NOx-Ox-HC chemistry run (NSRCX == 3)
               !==============================================================
               IF ( N == 1 ) THEN
                  IF ( NSRCX /= 3 ) CYCLE

                  IF ( IDTOX > 0 ) THEN 
                     TCBXXX = TCVV(IDTOX) * STT(I,J,L,IDTOX) / FDTT

                     CALL O3COMP(I,J,L)
                     
                     TCOBOX(KDA48,K) = TCBXXX * FO3(5)
                  ENDIF

               !==============================================================
               ! MS == 1, Tracer == 2: Store OH time series [molec/cm3/s] 
               !==============================================================
               ELSE IF ( N == 2 ) THEN
                  IF ( NSRCX /= 3 ) CYCLE

                  TCOBOX(KDA48,K) = SAVEOH(I,J,L)

               !==============================================================
               ! MS == 1, Tracer == 3:  Store NOy time series [v/v] 
               !
               ! Add up all of the NOy tracers, if they are switched on
               ! NOy = { NOX, PAN, HNO3, PMN, PPN, ISN2, R4N2, N2O5, HNO4 }
               !
               ! Only do this for an NOx-Ox-HC chemistry run (NSRCX == 3) 
               !==============================================================
               ELSE IF ( N == 3 ) THEN
                  IF ( NSRCX /= 3 ) CYCLE

                  TCBXXX = 0d0
                 
                  IF ( IDTNOX > 0 ) THEN
                     TCBXXX = TCBXXX + 
     &                    TCVV(IDTNOX) * STT(I,J,L,IDTNOX) / FDTT 
                  ENDIF

                  IF ( IDTPAN > 0 ) THEN
                     TCBXXX = TCBXXX + 
     &                    TCVV(IDTPAN) * STT(I,J,L,IDTPAN) / FDTT 
                  ENDIF
                  
                  IF ( IDTHNO3 > 0 ) THEN
                     TCBXXX = TCBXXX + 
     &                    TCVV(IDTHNO3) * STT(I,J,L,IDTHNO3) / FDTT 
                  ENDIF   
   
                  IF ( IDTPMN > 0 ) THEN
                     TCBXXX = TCBXXX + 
     &                    TCVV(IDTPMN) * STT(I,J,L,IDTPMN) / FDTT 
                  ENDIF

                  IF ( IDTPPN > 0 ) THEN
                     TCBXXX = TCBXXX + 
     &                    TCVV(IDTPPN) * STT(I,J,L,IDTPPN) / FDTT 
                  ENDIF

                  IF ( IDTISN2 > 0 ) THEN
                     TCBXXX = TCBXXX + 
     &                    TCVV(IDTISN2) * STT(I,J,L,IDTISN2) / FDTT 
                  ENDIF

                  IF ( IDTR4N2 > 0 ) THEN
                     TCBXXX = TCBXXX + 
     &                    TCVV(IDTR4N2) * STT(I,J,L,IDTR4N2) / FDTT 
                  ENDIF

                  IF ( IDTN2O5 > 0 ) THEN
                     TCBXXX = TCBXXX + 
     &                    2 * TCVV(IDTN2O5) * STT(I,J,L,IDTN2O5) / FDTT 
                  ENDIF
                     
                  IF ( IDTHNO4 > 0 ) THEN
                     TCBXXX = TCBXXX + 
     &                    TCVV(IDTHNO4) * STT(I,J,L,IDTHNO4) / FDTT 
                  ENDIF

                  TCOBOX(KDA48,K) = TCBXXX

               !==============================================================
               ! MS == 1, Tracer == 9: Store NO concentration [v/v]
               !
               ! Multiply NOx concentration by FNO to get NO concentration.
               ! Only do this for an O3 chemistry run (NSRCX == 3) 
               !==============================================================
               ELSE IF ( N == 9 ) THEN
                  IF ( NSRCX /= 3 ) CYCLE

                  IF ( IDTNOX > 0 ) THEN
                     TCBXXX = TCVV(IDTNOX) * STT(I,J,L,IDTNOX) / FDTT
                     
                     CALL O3COMP(I,J,L)
                     
                     TCOBOX(KDA48,K) = TCBXXX * FNO(5)
                  ENDIF

               !==============================================================
               ! MS = 1, Tracer = 10: Store PBL top (PS - PBL) [mb] 
               !==============================================================
               ELSE IF ( N == 10 ) THEN
                  TOPBL = GET_PEDGE(I,J,1) - PBL(I,J)
                  TCOBOX(KDA48,K) = TOPBL 

               !==============================================================
               ! MS == 1, Tracer == 11: Store Local Time [hours]
               ! 
               ! NOTE: It is recommended to save local time in the range 
               ! 0-24 hours, in order to avoid problems indexing local time
               ! in the IDL programs that read and plot timeseries data.
               !==============================================================
               ELSE IF ( N == 11 ) THEN
                  
                  ! Get local time
                  TCOBOX(KDA48,K) = GET_LOCALTIME( I )
                  
               !==============================================================
               ! MS == 1, Tracer == 12: Store ISOP Emissions [molec/cm2/s]
               !   
               ! Only do this for an NOx-Ox-HC chemistry run (NSRCX == 3)
               !==============================================================
               ELSE IF ( N == 12 ) THEN
                  IF ( NSRCX /= 3 ) CYCLE

                  TCOBOX(KDA48,K) = EMISRR48(I,J)

               !==============================================================
               ! MS == 1, Tracer == 14: Store NO3 mixing ratios [v/v] (bey)
               !   
               ! Only do this for a NOx-Ox-HC chemistry run (NSRCX == 3)
               ! Now get NO3 from the SAVENO3 array (amf, bmy, 4/25/00)
               !==============================================================
               ELSE IF ( N == 14 ) THEN
                  IF ( NSRCX /= 3 ) CYCLE
                  
                  TCOBOX(KDA48,K) = SAVENO3(I,J,L)

               !==============================================================
               ! MS == 1, Tracer == 15: Store Surface Pressure [mb]
               !==============================================================
               ELSE IF ( N == 15 ) THEN
                  TCOBOX(KDA48,K) = GET_PEDGE(I,J,1)

               !==============================================================
               ! MS == 1, Tracer == 16: Store Water Vapor [v/v]
               !==============================================================
               ELSE IF ( N == 16 ) THEN
                  IF ( AVGW_ALLOCATED ) THEN
                     TCOBOX(KDA48,K) = AVGW(I,J,L)
                  ELSE
                     TCOBOX(KDA48,K) = 0e0
                  ENDIF

               !==============================================================
               ! MS == 1, Tracer == 18: Store HCN column density [molec/cm2]
               !   
               ! NOTE: Assume HCN is tracer 1 (qli, bmy, 6/1/99)
               !==============================================================
               ELSE IF ( N == 18 .and. NSRCX == 4 ) THEN
                  TCBXXX =  ( SUM( STT(I,J,:,1) ) * 6.022d22         ) / 
     &                      ( TCMASS(1)           * GET_AREA_M2( J ) )

                  TCOBOX(KDA48,K) = TCBXXX

               !==============================================================
               ! MS == 1; Tracer == 19: Store temperature [K]                  
               !==============================================================
               ELSE IF ( N == 19 ) THEN
                  TCOBOX(KDA48,K) = T(I,J,L)
!
!*****************************************************************************
!  Any other tracer is an invalid diagnostic!!
!*****************************************************************************
!
               ELSE
                  PRINT*, 'Invalid TS station (I,J,L,MS,N) : ', 
     &                 I, J, L, MS, N
                  PRINT*, 'STOP in diag48.f'
                  CALL GEOS_CHEM_STOP
               ENDIF  

            ENDIF               ! MS = 0/1
         ENDDO
      ENDIF
!
!*****************************************************************************
!  Return to calling program
!*****************************************************************************
!
      END SUBROUTINE DIAG48
