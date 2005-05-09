! $Id: emisshcn.f,v 1.5 2005/05/09 14:33:58 bmy Exp $
      SUBROUTINE EMISSHCN
!
!******************************************************************************
!  Subroutine EMISSHCN specifies hydrogen cyanide (HCN) emissions 
!  (qli, bmy, 1/10/99, 7/20/04)
!
!  NOTES:
!  (1 ) The following sources of HCN are now used:
!        (a) HCN from biomass burning (scaled from CO values)
!  (2 ) Use F90 syntax (bmy, 3/24/99)
!  (3 ) Now use double-precision exponents, e.g. 1d0 (bmy, 3/18/99)
!  (4 ) Now pass SUNCOS and FIRSTEMISS as arguments (qli, 5/28/99)
!  (5 ) Now use allocatable arrays for ND28, ND46 diagnostics.
!        Eliminate reference to "CMN_O3" since we don't need that anymore.
!        Added a few minor cosmetic changes (bmy, 3/16/00)
!  (6 ) Now reference TS from "dao_mod.f" instead of from common block
!        header file "CMN_TS". (bmy, 6/23/00)
!  (7 ) Added reference to F90 modules "biomass_mod.f". Also, BURNEMIS is 
!        now referenced with global offsets IREF = I + I0 and JREF = J + J0. 
!        (bmy, 9/12/00)
!  (8 ) Removed obsolete code from 9/12/00 (bmy, 12/21/00)
!  (9 ) Now make sure IDBCO is not zero -- to avoid subscript range
!        errors when indexing BURNEMIS. (bmy, 3/20/01)
!  (10) Now prompt user to check IDBCO in "tracer.dat" if this switch
!        is turned off (bmy, 6/19/01)
!  (11) BURNEMIS(IDBCO,IREF,JREF) is now BURNEMIS(IDBCO,I,J).  Also DXYP
!        is dimensioned JGLOB, so use J+J0 to reference it. (bmy, 9/28/01)
!  (12) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (13) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (14) Removed obsolete code from 6/02.  Also commented out the sensitivity
!        cases from Larry Horowitz. (bmy, 8/26/02)
!  (15) Now reference BXHEIGHT from "dao_mod.f".  Now references ERROR_STOP 
!        from "error_mod.f".  Updated comments.  Now references IDBCO from
!        F90 module "tracerid_mod.f".  Now references SUNCOS from "dao_mod.f".
!        Now make FIRSTEMISS a local SAVEd variable. (bmy, 11/15/02)
!  (16) Now replace DXYP(J+J0)*1d4 with routine GET_AREA_CM2 from "grid_mod.f"
!        Now use function GET_TS_EMIS from "time_mod.f".  Removed MONTH from
!        call to BIOBURN. (bmy, 2/11/03)
!  (17) Now pass I, J to EMISOP (bmy, 12/9/03)
!  (18) Now references STT from "tracer_mod.f".  Now references LEMIS from
!        "logical_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BIOMASS_MOD,  ONLY : BURNEMIS, BIOBURN
      USE DAO_MOD,      ONLY : BXHEIGHT, TS,  SUNCOS
      USE DIAG_MOD,     ONLY : AD28,     AD46
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,  ONLY : LEMIS
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDBCO

      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! DIAGNOSTIC
#     include "CMN_MONOT"    ! Monoterpine variables
#     include "CMN_HCN"      ! HCN variables
#     include "hcn.h"        ! HCN header file w/ switches

      ! Local variables
      LOGICAL, SAVE          :: FIRSTEMISS = .TRUE.
      REAL*8                 :: DTSRCE, FLUX, T1L, AREA_CM2
      INTEGER                :: I, IREF, J, JREF, L, N

!-----------------------------------------------------------------------------
!  Define the following quantities for Biogenic HCN source 
!
#if   defined( BIOHCN )
      REAL*8                 :: TMMP, EMXX, EMX, EMXKG
      REAL*8                 :: TLAI,EMBIO,CLIGHT,XLTMMP,BIOFIT,TCORR
      REAL*8                 :: CONVERT(NVEGTYPE),GMONOT(NVEGTYPE)
      INTEGER                :: IJLOOP,INVEG

      ! HCNSCAL: Scale factor mol HCN / mol ISOP
      REAL*8, PARAMETER      :: HCNSCAL = 0.1d-2

      ! XNUMOL_ISOP: molecules ISOP / kg ISOP
      REAL*8, PARAMETER      :: XNUMOL_ISOP = 6.022d23 / 0.068d0

      ! External functions for Biogenic HCN
      REAL*8, EXTERNAL       ::  XLTMMP, EMISOP
#endif
!-----------------------------------------------------------------------------

      ! EHCN: Scale factor mol HCN / mol CO
      REAL*8, PARAMETER      :: EHCN = 1.1d-2

      ! XNUMOL_HCN: molecules HCN / kg HCN
      REAL*8, PARAMETER      :: XNUMOL_HCN = 6.022d23 / 0.027d0

      ! External functions
      REAL*8, EXTERNAL       :: BOXVL
!
!*****************************************************************************
!  EMISSHCN begins here!
!*****************************************************************************
!

      ! Return if we are not doing emissions ( if LSRCE = .FALSE. )
      IF ( .not. LEMIS ) RETURN

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0
      
!
!*****************************************************************************
!  HCN from Biomass burning
!
!  Biomass burning CO is stored in BURNEMIS(IDBCO,:,:) in molec/cm^3/s
!  Convert to kg HCN using the scale factors below.
!  DTSRCE   = # of seconds in this time interval
!
!  BURNEMIS is now referenced with global offset IREF = I + I0 and 
!  JREF = J + J0 (unlike previously).  (bmy, 9/12/00)
!*****************************************************************************
!
      ! Need to make sure CO is a turned-on biomass tracer
      IF ( IDBCO == 0 ) THEN
         CALL ERROR_STOP( 'IDBCO = 0 in "tracer.dat"!', 'emisshcn.f' )
      ENDIF

      ! Get biomass burning emissions 
      CALL BIOBURN
      
      N = 1
      L = 1
      DO J = 1, JJPAR

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )

      DO I = 1, IIPAR

         ! Emission flux in molec/cm^2/time step 
         FLUX = EHCN * BURNEMIS(IDBCO,I,J) * 
     &                 ( BXHEIGHT(I,J,L) * 100d0 )         

         ! Save to ND28 diagnostic as [molec/cm^2/time step]
         IF ( ND28 > 0 ) THEN
            AD28(I,J,N) = AD28(I,J,N) + FLUX
         ENDIF
         
         ! Convert to kg/grid box  
         T1L = FLUX * DTSRCE * AREA_CM2 / XNUMOL_HCN

         ! Add to tracer mass and global emission sum
         STT(I,J,L,N) = STT(I,J,L,N) + T1L 
         HCN_EMS = HCN_EMS + T1L  
         
      ENDDO
      ENDDO


!-----------------------------------------------------------------------------
!  A simple scale of biogenic HCN emission to that of Isoprene (qli,5/1/99)
!
#if   defined( BIOHCN )

      IF( FIRSTEMISS ) THEN
          CALL RDLIGHT
          CALL RDISOPT ( CONVERT )
          CALL RDMONOT ( GMONOT  )
          CALL SETBASE ( CONVERT, GMONOT )
      ENDIF

      FIRSTEMISS = .FALSE.

      N = 1
      IJLOOP = 0

      DO J = 1, JJPAR

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )

      DO I = 1, IIPAR
         IJLOOP = IJLOOP + 1
         TMMP = XLTMMP(I,J, IJLOOP) 
         EMXX = EMISOP(I,J,IJLOOP,SUNCOS,TMMP,XNUMOL_ISOP) 
  
         EMX          = EMXX
         EMXKG        = EMXX / XNUMOL_HCN * HCNSCAL
         STT(I,J,L,N) = STT(I,J,L,N) + EMXKG 

         ! ND46 -- archive as [atoms C/cm2/s] here (bmy, 9/13/01)
         IF ( ND46 > 0 ) THEN
            AD46(I,J,N) = AD46(I,J,N) + 
     &           ( EMXX * HCNSCAL ) / AREA_CM2 / DTSRCE
         ENDIF

      ENDDO
      ENDDO    

#endif
!-----------------------------------------------------------------------------

      ! Set FIRSTEMISS to .FALSE. since we have gone through one step
      FIRSTEMISS = .FALSE.

      ! Return to calling program
      END SUBROUTINE EMISSHCN

