! $Id: diag5.f,v 1.1 2003/06/30 20:26:01 bmy Exp $
      SUBROUTINE DIAG5                                                         
!
!*****************************************************************************
!  Subroutine DIAG5 (bmy, mgs, 5/8/98, 2/11/03) prints 
!  Time Series info to the "ctm.ts" file.
!
!  See comments to DIAG48 for an explanation of which quantities
!  are saved as time series.
!
!  NOTES:
!  (1 ) Use F90 syntax (bmy, 3/22/99).
!  (2 ) ### mgs, 24 Nov 1998: also uses TRCOFFSET (60, 70) to indicate 
!       Radon or CH3I runs.  TRCOFFSET is defined in NDXX_SETUP.
!  (3 ) Eliminate the ASCII map file ( LPRNT > 0 ), as well as
!       statistics output. (bmy, 3/22/99)
!  (4 ) TCOBOX is now declared allocatable in "diag_mod.f". (bmy, 11/29/99)
!  (5 ) Now reference file unit IU_TS from "file_mod.f" (bmy, 6/27/02)
!  (6 ) Now references GET_DIAGb, GET_DIAGe from "time_mod.f".  Now make 
!        NTAU0, NTAU local variables. (bmy, 2/11/03)
!*****************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD, ONLY : TCOBOX
      USE FILE_MOD, ONLY : IU_TS
      USE TIME_MOD, ONLY : GET_DIAGb, GET_DIAGe

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "CMN"
#     include "CMN_DIAG"
#     include "CMN_O3"

      ! Local variables
      REAL*8               :: ASCALB, TCBMAX, TCTMBX(MAXACC)

      INTEGER              :: I, ISCALB, J, L, N, MS, K, KD, KDMAX, NL
      INTEGER              :: NTAU0, NTAU
      CHARACTER ( LEN=16 ) :: TITLESTR
!
!*****************************************************************************
!  Print out time series
!*****************************************************************************
!
      IF ( ND48 > 0 .and. KDA48 > 0 ) THEN

         ! NTAU0, NTAU are just integers of DIAGb, DIAGe
         NTAU0 = GET_DIAGb()
         NTAU  = GET_DIAGe()
!
!*****************************************************************************
!  KDMAX is the number of accumulating time steps ( max = MAXACC )
!
!  Loop over each station K.  Convert I and J to global reference points.
!*****************************************************************************
!
         KDMAX = MIN( MAXACC, KDA48 )

         DO K = 1, MAXSTA 

            ! Get station coordinates...no longer use modulus for IMX
            I  = ISAVTC(K)
            J  = JSAVTC(K)
            L  = LSAVTC(K)
            MS = MSAVTC(K)
            N  = NSAVTC(K) 
            
            ! Ending line is denoted by I < 0 or J < 0
            IF ( I <= 0 .or. J <= 0 ) EXIT
            IF ( N <= 0 ) CYCLE
!
!*****************************************************************************
!  Loop over all the stations and find the maximum value of 
!  the RURAL BOX time series = TCBMAX.  
!
!  If the TCBMAX = 0 then all values are zero, so go to the next 
!  time series station.
!
!  Compute the RURAL BOX scale factor based on the value of TCBMAX
!
!  TCTMBX = TIME SERIES scaled data values 
!*****************************************************************************
!
            TCBMAX            = MAXVAL( TCOBOX( 1:KDMAX, K )  ) + 1d-30
            ISCALB            = INT( LOG10( TCBMAX ) + 100.0 ) - 100 
            ASCALB            = 10.0**( ISCALB - 1 )  

            TCTMBX( 1:KDMAX ) = TCOBOX( 1:KDMAX, K ) / ASCALB
!
!*****************************************************************************
!  Write Diagnostic output of full time series on unit IUNIT (ctm.ts)
!  Do NOT write out time series statistics anymore
!
!  NL is the number of lines for each data block, where each line 
!  contains 12 data values.
!*****************************************************************************
!
            NL = ( KDA48 + 11 ) / 12                                       

            WRITE ( IU_TS, '(2X, ''TB'', 6I5, 2I10, 1x, 1E10.3)' )
     &           I,J,L,MS,N,NL,NTAU0,NTAU,ASCALB
            
            WRITE ( IU_TS, '(12F6.2)' ) ( TCTMBX(KD), KD=1,KDMAX )  
         ENDDO

         ! Echo output
         WRITE( 6, '(a)' ) '     - DIAG48: Wrote output to ctm.ts!'
      ENDIF       

      ! Return to calling program
      END SUBROUTINE DIAG5
