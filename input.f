! $Id: input.f,v 1.2 2003/07/08 15:32:39 bmy Exp $
      SUBROUTINE INPUT
!
! NOTE: INPUT is mostly historical baggage, but it works for now.
!       We will eventually get rid of it (bmy, 11/6/02)
!
      ! References to F90 modules
      USE BIOMASS_MOD,  ONLY : SET_BIOTRCE 
      USE BIOFUEL_MOD,  ONLY : SET_BFTRACE
      USE DRYDEP_MOD,   ONLY : INIT_DRYDEP
      USE GRID_MOD,     ONLY : COMPUTE_GRID, SET_XOFFSET, SET_YOFFSET
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TIME_MOD,     ONLY : SET_TIMESTEPS
      USE TRACERID_MOD, ONLY : NEMANTHRO, NEMBIOG, TRACERID

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "CMN"
#     include "CMN_DIAG"
#     include "CMN_O3"
#     include "CMN_DEP"
#     include "comode.h"
#     include "CMN_GCTM"
#     include "CMN_SETUP" ! GEOS-CHEM flags

      ! Local variables
      INTEGER I, J, L, N, K, DUM1
      REAL*8  FJMID, FJBND, CPI180
      
      ! Now declare IM, JM, IMX, JMX as local variables
      ! since we have removed them from the common block
      INTEGER IM, JM, IMX, JMX

      !### Kludge: for now make NREAD, NWRITE, NDYN, NCONV, NSRCE,
      !### NCHEM, NDIAG local variables.  Call internal routine
      !### CHECK_TIMESTEPS to set the timesteps in time_mod.f
      !### 
      !### NOTE: We will replace "input.f" very soon, this is just
      !### a quick fix in order for us to proceed (bmy, 3/11/03)
      INTEGER NREAD, NWRITE, NDYN, NCONV, NSRCE, NCHEM, NDIAG, NUNIT
 
      !### Kludge: need to declare these as local variables temporarily
      INTEGER I0, J0, NTAU0, JMSIZE, IYEAR, I00, J00

      !=================================================================
      ! INPUT begins here!
      !=================================================================

      ! Initialize variables (bmy, 8/27/02)
      STT    = 0d0
      KDA48  = 0
      IYEAR  = 0 ! (bmy, 4/29/03)

      ! Now open unit #5 w/in "input.f" instead of w/in "open_files.f" 
      ! (bmy, 6/27/02)
      OPEN( 5, FILE='input.ctm', STATUS='OLD' )

      ! Now read XLABEL in as a CHAR*80 instead of 20 CHAR*4's (bmy, 5/27/99)
      READ (5,'(a)') XLABEL
      WRITE(6,960) XLABEL
      READ (5,902) LCONT,LSAVE,LINIT,LWINDO,LYRRD
      READ (5,902) LSTRAT,LTROP,LSRCE,LZONE,LPUFF
      READ (5,902) LSOM,LPFILT,LCONV4,LPLU,L2PM
      READ (5,902) LBIONOX,LAIRNOX,LLIGHTNOX,LSOILNOX,LFFNOX
      READ (5,902) LFOSSIL,LWOODCO,LUPBD,LBLMX,LEMBED
      READ (5,904) IEBD1,IEBD2,JEBD1,JEBD2

      ! For GEOS-CHEM, TAUI, TAUE are computed in MAIN.F
      ! Read a blank line from unit 5 (bmy, 10/9/97)
      READ(5,*)

C-------INPUT BASIC GRID DATA + SPECIFIC CTM DATA
      ! Use placeholders for NINST, NINIT, ICASE, LTM (bmy, 8/26/02)
      READ(5,904) NREAD,NWRITE,NDYN,NCONV,NSRCE,NCHEM,NDIAG,DUM1
      READ(5,904) JMSIZE,IMX,JMX,LM,DUM1,LCONVM,DUM1,DUM1

      ! IMRSLV and LMSTRT are GISS-CTM specific!
      ! Read 2 blank lines from unit 5 (bmy, 1/30/98)
      READ(5,*)
      READ(5,*)
      READ(5,*)

      ! For GEOS-CTM, SIGE is now read in SETUP.F
      ! Read 2 blank lines from unit 5 (bmy, 1/30/98)
      READ(5,*)
      READ(5,*)
      READ(5,904) NTRACE,I0,IM,J0,JM
      
      ! BMY added ISCALYR here for scaling of fossil fuels from 1985 data!
      READ(5,904) DUM1, NSRCX, DUM1, DUM1, FSCALYR

      ! If FSCALYR = 0 then turn off emissions by setting to 1985!
      IF ( FSCALYR .eq. 0 ) FSCALYR = 1985

      ! We now have scalefoss* data up to 1998 (bmy, 1/23/02)
      IF ( FSCALYR > 1998 ) THEN
         CALL ERROR_STOP( 'FSCALYR > 1998!', 'input.f' )
      ENDIF
      
      ! Read in NSKIPL -- use placeholders (bmy, 8/26/02)
      READ(5,904) DUM1, DUM1, DUM1, DUM1, NSKIPL

      ! Set LSOM, LPFILT, LCONV4, LPLU, LBLMX, LSTRAT = FALSE 
      ! Set LDIFF1, LDIFFM, KDIFFU, NSTRTC = 0
      ! This prevents errors for GEOS-CTM (bmy, 1/30/98)
      LSTRAT = .FALSE.
      LSOM   = .FALSE.
      LPLU   = .FALSE.
      LPFILT = .FALSE.
      LCONV4 = .FALSE.
      LBLMX  = .FALSE.

      ! Hardwire IM, JM, LM, IMX, JMX, LCONVM for safety's sake (bmy, 10/10/01)
      IM     = IIPAR
      JM     = JJPAR
      LM     = LLPAR
      IMX    = IGLOB
      JMX    = JGLOB
      LCONVM = LLCONVM

      ! Also make sure that NTRACE does not exceed NNPAR (bmy, 4/25/00)
      IF ( NTRACE > NNPAR ) THEN 
         CALL ERROR_STOP( 'NTRACE > NNPAR!', 'input.f' )
      ENDIF
      
      ! Read diagnostic dates andd NDxx diagnostic flags
      READ (5,901)
      READ (5,906) NJDAY
      READ (5,901)
      READ (5,905) ND01,ND02,ND03,ND04,ND05,ND06,ND07,ND08,ND09,ND10
      READ (5,905) ND11,ND12,ND13,ND14,ND15,ND16,ND17,ND18,ND19,ND20
      READ (5,905) ND21,ND22,ND23,ND24,ND25,ND26,ND27,ND28,ND29,ND30
      READ (5,905) ND31,ND32,ND33,ND34,ND35,ND36,ND37,ND38,ND39,ND40
      READ (5,905) ND41,ND42,ND43,ND44,ND45,ND46,ND47,ND48,ND49,ND50
      READ (5,905) ND51,ND52,ND53,ND54,ND55,ND56,ND57,ND58,ND59,ND60
      READ (5,905) ND61,ND62,ND63,ND64,ND65,ND66,ND67,ND68,ND69,ND70

      ! NDYN  now holds the value of NTDT (dyn step) as assigned in MAIN.F
      ! NREAD now holds the value of NDT  (6h  step) as assigned in MAIN.F
      ! NCHEM,  NSRCE, NCONV  must be multiples of NDYN
      ! NWRITE, NDIAG,        must be multiples of NREAD 
      ! Also replace label 26 with a new CONTINUE statement (bmy, 1/30/98)
 26   CONTINUE
      IF (mod(NCHEM,NDYN) .ne. 0) then
         CALL ERROR_STOP( 'NCHEM not a multiple of NDYN', 'input.f' )
      ENDIF
      IF (mod(NSRCE,NDYN) .ne. 0) THEN
         CALL ERROR_STOP( 'NSRCE not a multiple of NDYN', 'input.f' )
      ENDIF
      IF (mod(NCONV,NDYN) .ne. 0) THEN
         CALL ERROR_STOP( 'NCONV not a multiple of NDYN', 'input.f' )
      ENDIF
      IF (mod(NDIAG,NDYN) .ne. 0) THEN
         CALL ERROR_STOP( 'NDIAG not a multiple of NDYN', 'input.f' )
      ENDIF 
      IF (mod(NWRITE,NREAD) .ne. 0) THEN
         CALL ERROR_STOP( 'NWRITE not a multiple of NDYN', 'input.f' )
      ENDIF

 50   CONTINUE

      ! Set the following GCM specific variables = 0 (bmy, 1/30/98)
      LHEADR = .FALSE.
      I00    = 0
      J00    = 0

C---------------------------LIST INPUT----------------------------------
      WRITE(6,964) LCONT,LSAVE,LINIT,LWINDO,LYRRD,
     1             LSTRAT,LTROP,LSRCE,LZONE,LPUFF,
     2             LSOM,LPFILT,LCONV4,LPLU,L2PM,
     3             LBIONOX,LAIRNOX,LLIGHTNOX,LSOILNOX,LFFNOX,
     4             LFOSSIL,LWOODCO,LUPBD,LBLMX,LEMBED,
     5             IEBD1,IEBD2,JEBD1,JEBD2
      WRITE(6,965) NREAD,NWRITE,NDYN,NCONV,NSRCE,NCHEM,NDIAG,0
      WRITE(6,966) JMSIZE,IMX,JMX,LM,0,LCONVM,0,0
      WRITE(6,969) NTRACE,I0,IM,J0,JM,IYEAR,NSRCX,0, 0
      WRITE(6,970) 0, 0, 0, 0, NSKIPL
      WRITE(6,971)
      WRITE(6,906) NJDAY
C------------------CALCULATE SPHERICAL GEOMETRY---FIXED TO GLOBE--------
C--------JMSIZE=46 = OLD 7.83X10,  JMSIZE=45 = 8X10,  JMSIZE=90 = 4X5 GRID
 60   CONTINUE
      WRITE(6,980) DJSIZE,DISIZE

      !### Temporary: we will eventually replace "input.f" with 
      !### a new input file reader (bmy, 2/3/03)
      CALL SET_XOFFSET( I0 )
      CALL SET_YOFFSET( J0 )
      CALL COMPUTE_GRID
      
#if   defined( LGEOSCO ) 
      !  Redefine NSRCX to 5 if we are doing a CO run. (bmy, 4/16/00)
      NSRCX = 5
#endif

! Make sure photolysis is defined for certain simulations (bmy, 9/12/00)
#if   !defined( LFASTJ ) && !defined( LSLOWJ )
      SELECT CASE ( NSRCX )

         ! CH3I chemistry
         CASE ( 2, 3, 4 )
            CALL ERROR_STOP( 'LFASTJ and LSLOWJ are both undefined!',
     &                       'input.f' )
      
         CASE DEFAULT
            ! Nothing

      END SELECT
#endif

      !=================================================================
      ! Initialize timesteps in "time_mod.f" for now...we will 
      ! write new code for reading in the input files (bmy, 3/11/03)
      !=================================================================
      CALL CHECK_TIME_STEPS

      !=================================================================
      ! Call INPTR to read zonal mean concentrations from "inptr.ctm"
      ! This is mostly baggage, we'll replace this soon (bmy, 11/6/02)
      !=================================================================
      CALL INPTR

      !=================================================================
      ! Call TRACERID to set up ids for tracers, emissions, and biomass
      ! burning species.  Moved some code outside of TRACERID so as to
      ! avoid making circular references w/ F90 modules (bmy, 11/6/02)
      !
      ! NOTE: Now need to call INPTR before calling TRACERID, since
      !       it will read the TCNAME array (bmy, 11/12/02)
      !=================================================================
      CALL TRACERID

      ! Echo anthro & biogenic emitted tracers
      WRITE( 6, 110 ) IDEMS ( 1:NEMANTHRO+NEMBIOG )
 110  FORMAT( 'TRACERID: Emitted tracers (anthro & bio) :', 20i3 )

      ! Set NBIOTRCE in "biomass_mod.f"
      CALL SET_BIOTRCE

      ! Set NBFTRACE in "biofuel_mod.f"
      CALL SET_BFTRACE

      !=================================================================
      ! We have to call INIT_DRYDEP after TRACERID so that we can
      ! initialize the order of the drydep species.  Also, we have to
      ! do this before calling NDXX_SETUP, since NDXX_SETUP needs
      ! NUMDEP to be defined so that it can allocate the AD44 array.
      !=================================================================
      CALL INIT_DRYDEP

      !=================================================================
      ! Call NDXX_SETUP here...since NEMANTHRO, NEMBIOG, NBIOTRCE, and
      ! NBFTRACE will be defined.  At some point we should fold this 
      ! into an expanded "diag_mod.f". (bmy, 11/6/02) 
      !=================================================================
      CALL NDXX_SETUP

      !=================================================================
      ! Call READ_TINDEX to read the "diag.dat" file (bmy, 3/30/99)
      !=================================================================
      CALL READ_TINDEX

      ! Close the INPUT.CTM file (bmy, 3/29/99)
      CLOSE(5)
      RETURN
C-----------------------------------------------------------------------
C-----TERMINATE BECAUSE OF IMPROPER STARTUP CONDITIONS
 800  WRITE (6,990)
      CALL ERROR_STOP( 'STOP 800', 'input.f' )

C-----TERMINATE BECAUSE OF INCONSISTENCY WITH THE GCM WIND TAPE
 810  CONTINUE
      CALL ERROR_STOP( 'STOP 810', 'input.f' )
C-----------------------------------------------------------------------
  901 FORMAT (20A4)
  902 FORMAT (5L2)
  903 FORMAT (5E10.3)
  904 FORMAT (16I5)
  905 FORMAT (5X,25I3)
  906 FORMAT(5X,31I1/5X,29I1/5X,31I1/5X,30I1/5X,31I1/5X,30I1/5X,31I1/
     5 5X,31I1/5X,30I1/5X,31I1/5X,30I1/5X,31I1)
  907 FORMAT(8X,2L1,10I5)
  908 FORMAT(5X,7E10.3)
  909 FORMAT(5X,15F5.2)
C-----------------------------------------------------------------------
  960 FORMAT(' INPUT: ',20A4)
  961 FORMAT('  FULL TRACER CONTINUATION FROM TAU=',F10.2,I5/1X,20A4)
  962 FORMAT('  USE OLD TRACER TO INITIALIZE  TAU=',F10.2/1X,20A4)
  963 FORMAT(' BEGIN/END (HRS)=',2F12.2)
  964 FORMAT(' LCONT/LSAVE/LINIT/LWINDO/LYRRD:',5L2/
     F ' LSTRAT/LTROP/LSRCE/LZONE/LPUFF:',5L2/
     F ' LSOM/LPFILT/LCONV4/LPLU/L2PM:',5L2/
     F ' LBIONOX/LAIRNOX/LLIGHTNOX/LSOILNOX/LFFNOX:',5L2/
     F ' LFOSSIL/LWOODCO/LUPBD/LBLMX/LEMBED:',5L2/
     F ' IEBD1-2,JEBD1-2:',4I4/
     F ' LSRC1-2-3-4:',4L2)
  965 FORMAT(' NREAD/NWRITE/NDYN/NCONV/NSRCE/NCHEM/NDIAG/NINST',16I5)
  966 FORMAT(' JMSIZE/IMX/JMX/LM/LTM/LCONVM/NINIT/ICASE/:',17I5)
  967 FORMAT(' IMRSLV:',24I3)
  968 FORMAT(' PTOP/PSF=',2F10.2,'    SIGE='/(12F9.6))
  969 FORMAT(' NTRACE/I0-IM/J0-JM:',5I5/' IYEAR/NSRCX/II1/II2:',4I5)
  970 FORMAT(' LDIFF1/LDIFFM:',2I3,'  KDIFFU:',I5,'  NSTRTC:',I3,
     x       ' NSKIPL:',I3)
  971 FORMAT(' MON/',3('123456789.'),'1X')
  972 FORMAT(' ND..  1  2  3  4  5  6  7  8  9  .')
  974 FORMAT(' DIAG PRESSURE LEVELS:'/(10F10.4))
  975 FORMAT(' ND',I1,'-',10I3)
  980 FORMAT('   ----GRID SIZE =',F8.4,' DEG LAT BY',F5.1,' DEG LONG')
  990 FORMAT (' *****ERROR READING RESTART FILE ON UNIT=1')
  991 FORMAT (' *****ERROR ILLOGICAL WIND TAPE/TRACER RUN',/,
     1 '   NTAU/NTBEG/NTEND',10I10)

      CONTAINS

!-----------------------------------------------------------------------------

      SUBROUTINE CHECK_TIME_STEPS

      !=================================================================
      ! Internal subroutine CHECK_TIME_STEPS computes the smallest 
      ! dynamic time step for the model, based on which operations 
      ! are turned on.  
      !
      ! NOTE: moved this here from "main.f" as a stopgap measure.  We
      ! will totally rewrite the input file reader soon. (bmy, 3/11/03)
      !=================================================================

      ! Local variables
      INTEGER              :: I, J, K, L, NSMALLEST

      ! NUNIT is time step in minutes for unit conversion
      NUNIT = -1

      IF ( LTRAN .or. LCONV .or. LTURB ) THEN
         NUNIT = MAX( NDYN, NCONV )
      ENDIF

      ! Compute NSMALLEST as the minimum of NDYN, NCONV, NSRCE, NCHEM
      I = NDYN
      J = NCONV
      K = NSRCE
      L = NCHEM

      IF ( .not. LTRAN                  ) I = 999999 
      IF ( .not. LCONV .and..not. LTURB ) J = 999999
      IF ( .not. LDRYD .and..not. LEMIS ) K = 999999
      IF ( .not. LCHEM                  ) L = 999999

      NSMALLEST = MIN( I, J, K, L )


      ! If all of the operators above are turned off, 
      ! then set NSMALLEST to NDYN.
      IF ( NSMALLEST == 999999 ) THEN 
         NSMALLEST = NDYN
      ENDIF
      
      
      ! If NDYN is smaller than NSMALLEST, reset NSMALLEST
      ! to NDYN.  Also reset NTDT and NSTP accordingly.
      ! This is useful for runs where transport is turned off,
      ! but where chemistry is turned on. 
      IF ( NDYN < NSMALLEST ) THEN
         NDYN = NSMALLEST
      ENDIF

      ! Initialize timesteps in "time_mod.f"
      CALL SET_TIMESTEPS( CHEMISTRY=NCHEM, EMISSION=NSRCE, 
     &                    DYNAMICS=NDYN,   UNIT_CONV=NUNIT,
     &                    CONVECTION=NCONV )

      ! Return to MAIN program
      END SUBROUTINE CHECK_TIME_STEPS

!------------------------------------------------------------------------------

      ! Return to calling program
      END SUBROUTINE INPUT




