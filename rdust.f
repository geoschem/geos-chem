! $Id: rdust.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      SUBROUTINE RDUST( THISMONTH, THISYEAR )
!
!******************************************************************************
!  Subroutine RDUST reads global mineral dust concentrations as determined 
!  by P. Ginoux.  Calculates dust optical depth at each level for the
!  FAST-J routine "set_prof.f". (rvm, bmy, 9/30/00, 6/30/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Number of the current month (1-12)
!  (2 ) THISYEAR  (INTEGER) : 4-digit year number (e.g. 1996, 2001)
!
!  NOTES:
!  (1 ) RDUST was patterned after rdaerosol.f (rvm, 9/30/00)
!  (2 ) Don't worry about rewinding the binary file...reading from
!        binary files is pretty fast.  And it's only done once a month.
!  (3 ) Now references punch file utility routines from F90 module
!        "bpch2_mod.f".  Also reference variable DATA_DIR from the
!         header file "CMN_SETUP". (bmy, 9/30/00) 
!  (4 ) Now selects proper GEOS-STRAT dust field for 1996 or 1997.
!        Also need to pass THISYEAR thru the arg list. (rvm, bmy, 11/21/00)
!  (5 ) CONC is now declared as REAL*8 (rvm, bmy, 12/15/00)
!  (6 ) Removed obsolete code from 12/15/00 (bmy, 12/21/00)
!  (7 ) CONC(IGLOB,JGLOB,LGLOB,NDUST) is now CONC(IIPAR,JJPAR,LLPAR,NDUST).
!        Now use routine TRANSFER_3D from "transfer_mod.f" to cast from REAL*4
!        to REAL*8 and also to convert from {IJL}GLOB to IIPAR,JJPAR,LLPAR 
!        space.  Use 3 arguments in call to GET_TAU0.  Updated comments.
!        (bmy, 9/26/01)
!  (8 ) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (9 ) Now reference ERADIUS, IXSAVE, IYSAVE, IZSAVE, TAREA from 
!        "comode_mod.f".  Compute ERADIUS and TAREA for the NDUST dust
!        size bins from FAST-J.  Renamed CONC to DUST to avoid conflicts.
!        Also reference NTTLOOP from "comode.h".  Also added parallel
!        DO-loops.  Also renamed MONTH and YEAR to THISMONTH and THISYEAR
!        to avoid conflicts w/ other variables. (bmy, 11/15/01)
!  (10) Bug fix: Make sure to use 1996 dust data for Dec 1995 for the
!        GEOS-STRAT met field dataset.  Set off CASE statement with an
!        #if defined( GEOS_STRAT ) block. (rvm, bmy, 1/2/02)
!  (11) Eliminate obsolete code from 1/02 (bmy, 2/27/02)
!  (12) Now report dust optical depths in ND21 diagnostic at 400 nm.  Now
!       report dust optical depths as one combined diagnostic field instead 
!        of 7 separate fields.  Now reference JLOP from "comode_mod.f".  
!        Now save aerosol surface areas as tracer #5 of the ND21 diagnostic.  
!        (rvm, bmy, 2/28/02)
!  (13) Remove declaration for TIME, since that is also defined in the
!        header file "comode.h" (bmy, 3/20/02)
!  (14) Now read mineral dust files directly from the DATA_DIR/dust_200203/
!        subdirectory (bmy, 4/2/02)
!  (15) Now reference BXHEIGHT from "dao_mod.f".  Also reference ERROR_STOP
!        from "error_mod.f". (bmy, 10/15/02)
!  (16) Now call READ_BPCH2 with QUIET=TRUE to suppress extra informational
!        output from being printed.  Added cosmetic changes. (bmy, 3/14/03)
!  (17) Since December 1997 dust data does not exist, use November 1997 dust
!        data as a proxy. (bnd, bmy, 6/30/03)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE COMODE_MOD,   ONLY : ERADIUS, IXSAVE, IYSAVE, 
     &                         IZSAVE,  JLOP,   TAREA
      USE DAO_MOD,      ONLY : BXHEIGHT
      USE DIAG_MOD,     ONLY : AD21
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TRANSFER_MOD, ONLY : TRANSFER_3D

      IMPLICIT NONE

#     include "cmn_fj.h"   ! LPAR, CMN_SIZE
#     include "jv_cmn.h"   ! ODMDUST, QAA, RAA
#     include "CMN_DIAG"   ! ND21, LD21
#     include "CMN_SETUP"  ! DATA_DIR
#     include "comode.h"   ! NTTLOOP

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH, THISYEAR

      ! Local variables
      INTEGER             :: I, J, JLOOP, L, N
      INTEGER, SAVE       :: MONTH_LAST = -999
      REAL*4              :: TEMP(IGLOB,JGLOB,LGLOB)
      REAL*8              :: DUST(IIPAR,JJPAR,LLPAR,NDUST)
      REAL*8              :: MSDENS(NDUST), XTAU
      CHARACTER (LEN=255) :: FILENAME

      !=================================================================
      ! RDUST begins here!
      !
      ! Read aerosol data from the binary punch file during the first 
      ! chemistry timestep and, after that, at the start of each month.
      !=================================================================
      IF ( THISMONTH /= MONTH_LAST ) THEN   
         
         ! Save the current month
         MONTH_LAST = THISMONTH

         ! Get TAU0 value used to index the punch file
         ! Use the "generic" year 1985
         XTAU = GET_TAU0( THISMONTH, 1, 1985 )
         
#if   defined( GEOS_STRAT )

         ! Select proper dust file name for GEOS-STRAT (1996 or 1997 data)
         SELECT CASE ( THISYEAR )

            ! GEOS-STRAT -- 1996 dust fields from P. Ginoux
            ! Since GEOS-STRAT covers Dec 1995, use the 1996 file as a
            ! proxy.  1995 dust data doesn't exist (rvm, bmy, 1/2/02)
            CASE ( 1995, 1996 )
               FILENAME = TRIM( DATA_DIR ) // 'dust_200203/dust.' //
     &                    GET_NAME_EXT()   // '.'                 // 
     &                    GET_RES_EXT()    // '.1996'

            ! GEOS-STRAT -- 1997 dust fields from P. Ginoux
            CASE ( 1997 )
               FILENAME = TRIM( DATA_DIR ) // 'dust_200203/dust.' //
     &                    GET_NAME_EXT()   // '.'                 // 
     &                    GET_RES_EXT()    // '.1997'

               ! KLUDGE -- there isn't dust data for December 1997, so 
               ! just use November's data for December (bnd, bmy, 6/30/03)
               IF ( THISMONTH == 12 ) THEN
                  XTAU = GET_TAU0( 11, 1, 1985 )
               ENDIF

            ! Error: THISYEAR is outside valid range for GEOS-STRAT
            CASE DEFAULT
               CALL ERROR_STOP( 'Invalid GEOS-STRAT year!', 'rdust.f' )

         END SELECT

#else
         ! Select proper dust file name for GEOS-1, GEOS-3, or GEOS-4
         FILENAME = TRIM( DATA_DIR ) // 'dust_200203/dust.' //
     &              GET_NAME_EXT()   // '.'                 // 
     &              GET_RES_EXT()

#endif

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - RDUST: Reading ', a )

         ! Read aerosol concentrations [kg/m3] for each 
         ! dust type from the binary punch file
         DO N = 1, NDUST 
            CALL READ_BPCH2( FILENAME, 'MDUST-$', N,     XTAU,
     &                       IGLOB,     JGLOB,    LGLOB, TEMP, 
     &                       QUIET=.TRUE. )

            CALL TRANSFER_3D( TEMP, DUST(:,:,:,N) )
         ENDDO

         !==============================================================
         ! Convert concentration [kg/m3] to optical depth [unitless].
         !
         ! ODMDUST = ( 0.75 * BXHEIGHT * CONC * QAA ) / 
         !           ( MSDENS * RAA * 1e-6 )
         ! (see Tegen and Lacis, JGR, 1996, 19237-19244, eq. 1)
         !
         !  Units ==> DUST     [ kg/m3    ]
         !            MSDENS   [ kg/m3    ]
         !            RAA      [ um       ]
         !            BXHEIGHT [ m        ]
         !            QAA      [ unitless ]
         !            ODMDUST  [ unitless ]
         !
         ! NOTES: 
         ! (1) Do the calculation at QAA(4,:) (i.e. 999 nm).          
         !==============================================================
         MSDENS(1) = 2500.0
         MSDENS(2) = 2500.0
         MSDENS(3) = 2500.0
         MSDENS(4) = 2500.0
         MSDENS(5) = 2650.0
         MSDENS(6) = 2650.0
         MSDENS(7) = 2650.0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, NDUST
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ODMDUST(I,J,L,N) = 0.75d0        * BXHEIGHT(I,J,L) * 
     &                         DUST(I,J,L,N) * QAA(4,14+N)     / 
     &                        ( MSDENS(N) * RAA(4,14+N) * 1.0D-6 )
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Echo information
         WRITE( 6, 110 ) 
 110     FORMAT( '     - RDUST: Finished computing optical depths' )

         !==============================================================
         ! Calculate Dust Surface Area
         !
         ! Units ==> DUST     [ kg dust/m^3 air    ]
         !           MSDENS   [ kg dust/m^3 dust   ]
         !           RAA      [ um                 ]
         !           TAREA    [ cm^2 dust/cm^3 air ]
         !           ERADIUS  [ cm                 ]
         !
         ! NOTE: first find volume of dust (cm3 dust/cm3 air), then 
         !       multiply by 3/radius to convert to surface area in cm2
         !  
         ! TAREA(:,1:NDUST) and ERADIUS(:,1:NDUST) are for 
         ! the NDUST FAST-J dust wavelength bins (read into DUST)
         !==============================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L, N )
         DO N     = 1, NDUST
         DO JLOOP = 1, NTTLOOP

            ! Compute 3-D grid box indices
            I = IXSAVE(JLOOP)
            J = IYSAVE(JLOOP)
            L = IZSAVE(JLOOP)

            ERADIUS(JLOOP,N) = RAA(4,14+N) * 1.0D-4

            TAREA(JLOOP,N)   = 3.D0 / ERADIUS(JLOOP,N) *
     &                         DUST(I,J,L,N) / MSDENS(N)  
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         !==============================================================
         ! ND21 Diagnostic: 
         !
         ! Tracer #1: Cloud optical depths    (from "optdepth_mod.f")
         ! Tracer #2: Max Overlap Cld Frac    (from "optdepth_mod.f")
         ! Tracer #3: Random Overlap Cld Frac (from "optdepth_mod.f")
         ! Tracer #4: Dust optical depths at 400 nm (from all size bins)
         ! Tracer #5: Dust surface areas (from all size bins)
         !==============================================================
         IF ( ND21 > 0 ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JLOOP, L, N ) 
            DO N = 1, NDUST
            DO L = 1, LD21
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               !--------------------------------------
               ! ND21 tracer #4: Dust optical depths
               !--------------------------------------
               AD21(I,J,L,4) = AD21(I,J,L,4) + 
     &            ( ODMDUST(I,J,L,N) * QAA(2,14+N) / QAA(4,14+N) )

               !--------------------------------------
               ! ND21 tracer #5: Dust surface areas
               !--------------------------------------
               IF ( L <= LLTROP ) THEN

                  ! Convert 3-D indices to 1-D index
                  ! JLOP is only defined in the tropopause
                  JLOOP = JLOP(I,J,L)
             
                  ! Add to AD21
                  IF ( JLOOP > 0 ) THEN
                     AD21(I,J,L,5) = AD21(I,J,L,5) + TAREA(JLOOP,N)
                  ENDIF
               ENDIF
            ENDDO
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ENDIF 
      ENDIF

      ! Return to calling program
      END SUBROUTINE RDUST
