! $Id: chemo3.f,v 1.1 2003/06/30 20:26:08 bmy Exp $
      SUBROUTINE CHEMO3
!
!******************************************************************************
!  Subroutine CHEMO3 performs ozone chemistry with specified production
!  and loss rates (amf, bey, bmy, 6/9/99, 2/11/03)
! 
!  Arguments as input:
!  ==========================================================================
!  (1 ) NYMD (INTEGER) : Current YYYYMMDD value
!  (2 ) NHMS (INTEGER) : Current HHMMSS value
!
!  NOTES:
!  (1 ) This subroutine should be called once every NCHEM minutes.
!  (2 ) Now use IOERROR to trap I/O errors (bmy, 6/9/99)
!  (3 ) Now remove dry deposition losses for NO2 and O3 from Ox.
!        Also DO-loops and IF statements for more efficient execution.
!        (bmy, 11/23/99) 
!  (4 ) LLTROP is now the highest tropospheric level underneath the
!        annual mean tropopause.  
!  (5 ) Prior to 12/13/99, P(O3) and L(O3) rates were only saved up to 
!        level 14.  With the annual mean tropopause code now in existence,
!        the troposphere can now extend as high as level 16.  For now, just
!        copy the 14th level of P24H and L24H into levels 15 and 16.  This
!        is necessary in order to read existing P(O3) and L(O3) rates.
!        Replace this when another full-year run is done.  (bmy, 12/13/99)
!  (6 ) Now declare P24H and L24H as allocatable in "DIAG_MOD.F".  Also
!        now only read the P(O3) and L(O3) rates from the binary file
!        once per day.  This will save time. (bmy, 12/13/99)
!  (7 ) Add changes from amf for multi-tracer Ox -- now call 
!        subroutine CHEMO3_SPLIT (amf, bmy, 7/5/01)
!  (8 ) Removed hardwired stuff and daily loss files -- amf says that these 
!        shouldn't be needed for other users. (amf, bmy, 8/28/01)
!  (9 ) Bug fix: now use true exponential loss for drydep instead of just 
!        a first-order approximation. (bdf, bmy, 7/11/02)
!  (10) Bug fix: now use true exponential loss for drydep instead of just 
!        a first-order approximation. (bdf, bmy, 7/11/02)
!  (11) Now read rate files in binary punch file format.  Also deleted 
!        obsolete code. (bmy, 7/31/02)
!  (12) Now references ALLOC_ERR from "error_mod.f".  Also now make FIRSTCHEM 
!        a local SAVEd variable.  Now references DEPSAV from "drydep_mod.f".
!        Eliminate extraneous dimension in DEPSAV.  Now eliminate reference
!        to CMN_SAV header file. (bmy, 11/19/02)
!  (13) Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2 of "grid_mod.f"
!        Now removed NYMD, NHMS from the arg list.  Now references functions
!        DATE_STRING, GET_NHMS, GET_MONTH, GET_DAY, GET_YEAR, GET_TS_CHEM,
!        GET_TAU from the new "time_mod.f"(bmy, 2/11/03)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIAG_MOD,   ONLY : AD44, AD65, P24H, L24H
      USE DRYDEP_MOD, ONLY : DEPSAV
      USE ERROR_MOD,  ONLY : ALLOC_ERR
      USE GRID_MOD,   ONLY : GET_AREA_M2
      USE TIME_MOD,   ONLY : DATE_STRING, GET_NHMS, GET_TS_CHEM,
     &                       GET_MONTH,   GET_DAY,  GET_YEAR,   GET_TAU

      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters 
#     include "CMN"          ! STT, NSRCX, MONTH, JDATE, YEAR, etc
#     include "CMN_DIAG"     ! Diagnostic switches
#     include "CMN_O3"       ! 
#     include "CMN_SETUP"    ! DATA_DIR
#     include "CMN_TIMES"    ! for dynamic allocation for US array.

      ! Local variables
      LOGICAL, SAVE          :: FIRSTCHEM = .TRUE.
      INTEGER                :: I, J, L, N, NN,  IOS, GMTRC, GMNL
      REAL*4                 :: ARRAY(IGLOB,JGLOB,LLTROP) 
      REAL*8                 :: DTCHEM, XTAU, PP, LL, PL, DTC
      CHARACTER(LEN=80)      :: TOPTITLE
      CHARACTER(LEN=255)     :: FILENAME

      ! External routines
      REAL*8, EXTERNAL       :: BOXVL

      !=================================================================
      ! CHEMO3 begins here!
      !
      ! Convert NCHEM from mn to sec, store in DTCHEM
      !=================================================================
      WRITE( 6, '( '' --- CHEMO3 - at : '', f10.2 )' ) GET_TAU()

      DTCHEM   = GET_TS_CHEM() * 60d0
      TOPTITLE = 'GEOS-CHEM time series for geographical domain'

      !=================================================================
      ! If this is the first time through the subroutine, allocate the
      ! P24H and L24H arrays.  This will allow us not to declare these 
      ! large arrays if we are not doing a single tracer O3 run.    
      ! 
      ! Allocate CHEML24, DRYDL24, CTCHDD only for area specified in 
      ! timeseries.dat -- need to include CMN_TIMES for this!
      ! (amf, 6/28/01)
      !=================================================================
      IF ( FIRSTCHEM ) THEN
         ALLOCATE( P24H( IIPAR, JJPAR, LLTROP ), STAT=IOS )
         IF ( IOS > 0 ) CALL ALLOC_ERR( 'P24H' )

         ALLOCATE( L24H( IIPAR, JJPAR, LLTROP ), STAT=IOS )
         IF ( IOS > 0 ) CALL ALLOC_ERR( 'L24H' )         
      ENDIF

      !=================================================================
      ! If it's a new day, read P(O3) and L(O3) from disk
      !=================================================================
      IF ( FIRSTCHEM .or. GET_NHMS() == 000000 ) THEN

         ! Define the data directory (uncomment the one you want)
         FILENAME = '/data/ctm/GEOS_MEAN/O3_PROD_LOSS/'
         !FILENAME = TRIM( FILENAME ) // '1991v4.6/'
         !FILENAME = TRIM( FILENAME ) // '1994/'
         !FILENAME = TRIM( FILENAME ) // '1994v4.6/'
         !FILENAME = TRIM( FILENAME ) // '1995/'
         FILENAME = TRIM( FILENAME ) // '1996-97r2x25v4.26.aer/'
         !FILENAME = TRIM( FILENAME ) // '1996v4.6/'
         !FILENAME = TRIM( FILENAME ) // '1997v4.11/'
         !FILENAME = TRIM( FILENAME ) // '1998/'
         !FILENAME = TRIM( FILENAME ) // '2000v4.26/'
         !FILENAME = TRIM( FILENAME ) // '2001v4.26/'

         ! Pick the filename
         FILENAME  = TRIM( FILENAME ) // 'rate.' // DATE_STRING()

         ! Echo information
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( ' --- CHEMO3 - reading ', a )

         ! Get the TAU0 value for today
         XTAU = GET_TAU0( GET_MONTH(), GET_DAY(), GET_YEAR() )

         !==============================================================
         ! Read P(O3)
         !==============================================================
         CALL READ_BPCH2( FILENAME, 'PORL-L=$', 1,      XTAU,  
     &                    IGLOB,    JGLOB,      LLTROP, ARRAY )

         ! Cast from REAL*4 to REAL*8 
         P24H(:,:,1:LLTROP) = ARRAY(:,:,1:LLTROP)

         !==============================================================
         ! Read L(O3)
         !==============================================================
         CALL READ_BPCH2( FILENAME, 'PORL-L=$', 2,      XTAU,  
     &                    IGLOB,    JGLOB,      LLTROP, ARRAY )

         ! Cast from REAL*4 to REAL*8 
         L24H(:,:,1:LLTROP) = ARRAY(:,:,1:LLTROP)

         ! We have gone through one iteration
         FIRSTCHEM = .FALSE.

      ENDIF

      !=================================================================
      ! Change tracer array: New STT = Old STT + ( Production - Loss )
      !
      ! BOXVL  = Volume of grid box            [ cm3             ]
      ! DTCHEM = Length of chemistry timestep  [ sec             ]
      ! PP     = Production of O3              [ kg/box/timestep ]
      ! LL     = Loss       of O3              [ kg/box/timestep ]
      !
      ! amf - for 12 tracer geographically split O3 run, call 
      !       subroutine CHEMO3_SPLIT, modeled on code below, but kept 
      !       separately to permit single tracer O3 run to be used. 
      !=================================================================
      IF ( LSPLIT ) THEN

         ! For geographcially tagged Ox
         CALL CHEMO3_SPLIT
      ELSE

         ! For single tracer Ox
         DO L = 1, LLTROP
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Skip over stratospheric boxes
            IF ( L >= LPAUSE(I,J) ) CYCLE

            PP = P24H(I,J,L) * BOXVL(I,J,L) * DTCHEM

            !-----------------------------------------------------------
            ! Debug...change PP for testing for different levels
            !PP = 0D0                  ! strato only
            !IF (L .LE. 11) PP = 0d0   ! upper trop
            !IF (L .LE. 7 ) PP = 0D0   ! mid   trop
            !IF (L .GE. 12) PP = 0D0   ! mid   trop
            !IF (L .GE. 8)  PP = 0D0   ! lower trop
            !-----------------------------------------------------------

            LL = STT(I,J,L,1) * L24H(I,J,L) * BOXVL(I,J,L) * DTCHEM 

            !===========================================================
            ! Archive ND65 diagnostic 
            ! (You have a choice of units, see below)
            !
            ! Save O3 Production in Ad65(I,J,L,1)
            ! Save O3 Loss       in AD65(I,J,L,2)
            !
            ! NOTES: 
            ! (1) Make sure that you fix the unit string for the binary
            !      punch file in "diag3.f" so that it coincides with 
            !      the unit you have selected.            
            !===========================================================
            IF ( ND65 > 0 ) THEN

               !--------------------------------------------------------
               ! To save O3 production in ppb/s, uncomment this line:
               !PL = P24H(I,J,L) * BOXVL(I,J,L) * 
               !     TCVV(1)     / AD(I,J,L)
               !--------------------------------------------------------
               ! To save O3 production in kg/s, uncomment this line:
               !PL = P24H(I,J,L) * BOXVL(I,J,L) 
               !--------------------------------------------------------
               ! To save O3 production in molec/s, uncomment this line:
               !PL = P24H(I,J,L) * BOXVL(I,J,L) * XNUMOL(1)
               !--------------------------------------------------------
               ! To save O3 production in molec/cm3/s, uncomment this line:
               PL = P24H(I,J,L) * XNUMOL(1)
               !---------------------------------------------------------
               !amf check canadian box
               !if (i.eq. 20 .and. j.eq.36 .and.l.eq.1) then
               !   write(444,*) 'read = ', P24H(I,J,L)
               !   write(444,*) 'XNUMOL=', XNUMOL(1)
               !   write(444,*) 'PROD = ', PL
               !endif 
               !
               !if (i.eq. 20 .and. j.eq.33 .and.l.eq.1) then
               !   write(445,*) 'read = ', P24H(I,J,L)
               !   write(445,*) 'XNUMOL=', XNUMOL(1)
               !   write(445,*) 'PROD = ', PL
               !endif 
               !---------------------------------------------------------

               AD65(I,J,L,1) = AD65(I,J,L,1) + PL

!-----------------------------------------------------------------------------
! To save O3 loss in ppb/s, uncomment this line:
!               PL = STT(I,J,L,1) * L24H(I,J,L) * BOXVL(I,J,L) *
!     &              TCVV(1)      / AD(I,J,L)
!-----------------------------------------------------------------------------
! To save O3 loss in kg/s, uncomment this line
!               PL = STT(I,J,L,1) * L24H(I,J,L) * BOXVL(I,J,L) 
!-----------------------------------------------------------------------------
! To save O3 loss in molec/s, uncomment this line:
!               PL = STT(I,J,L,1) * L24H(I,J,L) * 
!     &              BOXVL(I,J,L) * XNUMOL(1)
!-----------------------------------------------------------------------------
! To save O3 loss in molec/cm3/s, uncomment this line:
               PL = STT(I,J,L,1) * L24H(I,J,L) * XNUMOL(1)
!-----------------------------------------------------------------------------

               !### Debug
               !if (i.eq.20.and.j.eq.33.and.l.eq.1) then
               !  write(*,*) 'PL = ', PL, ND65
               !endif

               AD65(I,J,L,2) = AD65(I,J,L,2) + PL

            ENDIF               !nd65

            !===========================================================
            ! Perform dry deposition losses and also archive them to 
            ! the ND44 diagnostic (if necessary).  Drydep losses are 
            ! in kg/box/sec. 
            !
            ! N      = index for DEPSAV array (N=1 is O3, N=2 is NO2)
            ! DEPSAV = Drydep loss frequency of O3 [s^-1]
            ! DTCHEM = Chemistry time interval     [s  ]   
            !===========================================================
            IF ( L == 1 ) THEN 
               IF ( LDRYD ) THEN

                  ! Loop over drydep species (for now, just do O3)
                  DO N = 1, 1

                     ! Save fluxes to ND44
                     IF ( ND44 > 0 ) THEN

!-----------------------------------------------------------------------
! TO SAVE In kg/sec.  Don't forget to change units in diag3.f!!
!                        DTC = DEPSAV(I,J,N) * STT(I,J,1,N)
!------------------------------------------------------------------------
! TO SAVE IN molec/cm^2/s for comparison with full chemistry run, 
! uncomment this line: (convert dxyp to cm^2)
                        DTC = DEPSAV(I,J,N) * STT(I,J,1,N) * 
     &                        XNUMOL(1)     * 1e-4 / GET_AREA_M2(J)
!------------------------------------------------------------------------

                        AD44(I,J,N,1) = AD44(I,J,N,1) + DTC
                     ENDIF
                  ENDDO         

                  ! Apply drydep loss to STT tracer array --
                  ! Only if drydep is turned on and for 1st layer only
                  STT(I,J,L,1) = STT(I,J,L,1) * 
     &                           EXP( -DEPSAV(I,J,N) * DTCHEM )

               ENDIF
            ENDIF               

            ! Apply production & loss to STT -- all layers
            STT(I,J,L,1) = STT(I,J,L,1) + PP - LL

            ! Warn if STT is negative
            IF ( STT(I,J,L,1) < 0d0 ) THEN
               WRITE( 6, 110 ) I, J, L, GET_TAU()
               !STOP
            ENDIF
         ENDDO
         ENDDO
         ENDDO

      ENDIF !endif else of LSPLIT -- single tracer O3 run amf. 

      RETURN

      ! FORMAT statements
 110  FORMAT( '### CHEMO3: STT < 0 AT (I,J,L) = ', 3i5, 
     &        ' at TAU = ', f10.2, ' (before drydep)' )
 120  FORMAT( '### CHEMO3: STT < 0 AT (I,J,L) = ', 3i5,
     &        ' at TAU = ', f10.2, ' (after drydep)' )

      ! Return to calling program
      END SUBROUTINE CHEMO3
