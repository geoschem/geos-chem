! $Id: inptr.f,v 1.1 2003/06/30 20:26:01 bmy Exp $
      SUBROUTINE INPTR                                                       
      
      !=================================================================
      ! NOTE: The original INPTR was historical baggage from the old
      ! GISS-II 9-layer model.  Most of the features in the "inptr.ctm"
      ! file are now obsolete for GEOS-CHEM.  As a stopgap measure we
      ! have rewritten the "inptr.ctm" and changed it so that tracers
      ! are listed vertically instead of horizontally.  At some point
      ! we will totally redesign the input file reading subroutines
      ! so that INPTR will no longer be necessary. (bmy, 11/20/02)
      !=================================================================

      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IOERROR

      IMPLICIT NONE
        
#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! STT
#     include "CMN_DIAG"   ! ISAVTC, JSAVTC, LSAVTC, MSAVTC, NSAVTC
#     include "CMN_O3"     ! XNUMOL

      ! Local variables
      INTEGER :: I, J, L, N, K, IS, JS, LS, MS, NS, IOS

      !=================================================================     
      ! INPTR begins here!
      !=================================================================
      OPEN( 9, FILE='inptr.ctm', STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 9, 'inptr:1' )
      
      ! Read header line
      READ( 9, '(a)', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 9, 'inptr:2' )       

      ! Loop over tracers
      DO N = 1, NTRACE

         ! Read tracer name, and molecular weight [g/mole]
         READ( 9, '(6x,a4,1x,f7.1)', IOSTAT=IOS ) TCNAME(N), TCMASS(N) 
         IF ( IOS /= 0 ) CALL IOERROR( IOS, 9, 'inptr:3' )

         ! FMOL are the molecular weights [kg/mole]
         FMOL(N)   = TCMASS(N) * 1d-3

         ! XNUMOL is the ratio [molec tracer/kg tracer]
         XNUMOL(N) = 6.022d+23 / FMOL(N)  
 
         ! TCVV is the ratio MW air / MW tracer
         TCVV(N)   = 28.97d0   / TCMASS(N)
      ENDDO

      ! Also compute the ratio [molec air/kg air]
      XNUMOLAIR = 6.022d+23 * 1.0d+03 / 28.9644d0

      ! Read an extra line from "inptr.ctm"
      READ( 9, '(a)', IOSTAT=IOS ) 
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 9, 'inptr:4' )

      !=================================================================
      ! Print out tracer information in a better way (bmy, 11/20/02)
      !=================================================================
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'INPTR: Tracer Quantities'
      WRITE( 6, '(a)'   ) '  #   Name   Molec.    TCVV     XNUMOL'
      WRITE( 6, '(a)'   ) '             Weight             molec/kg'
      WRITE( 6, '(a)'   ) REPEAT( '-', 42 )

      DO N = 1, NTRACE
         WRITE( 6, '(i3,3x,a4,2(3x,f6.3),3x,es10.3)' ) 
     &        N, TRIM( TCNAME(N) ), FMOL(N), TCVV(N), XNUMOL(N)
      ENDDO
    
      !=================================================================
      ! Read information about stations for time series diagnostic (48). 
      ! The diagnostic 48 has to be turn on in input.ctm.
      ! See also diag48.f for description of the different quantities 
      ! available
      !
      ! ### mgs, 24 Nov 1998: changed to automatically assign stations
      ! ### from level 1 to level indicated in inptr.cmt
      ! 
      ! ### Replaced obsolete Prather code.  Also halt if the 
      ! ### station coordinates exceed array boundaries. (bmy, 1/5/00)
      !=================================================================

      ISAVTC(:) = 0
      JSAVTC(:) = 0
      LSAVTC(:) = 0
      MSAVTC(:) = 0
      NSAVTC(:) = 0

      WRITE( 6, '(/,a)' ) 'ND48 Timeseries Stations:'

      ! Loop over each station
      K = 1
      DO WHILE ( K <= NNSTA ) 

         ! Read station coordinates from disk
         READ( 9, '(5i5)', IOSTAT=IOS ) IS, JS, LS, MS, NS
         IF ( IOS /= 0 ) CALL IOERROR( IOS, 9, 'inptr:5' )

         WRITE( 6, '(3i5,''-'',i2,2I5)' ) IS, JS, 1, LS, MS, NS

         ! Exit the loop if IS and JS are less than 1
         IF ( IS < 1 .and. JS < 1 ) EXIT

         ! Make sure IS is not larger than IIPAR
         IF ( IS > IIPAR ) THEN 
            CALL ERROR_STOP( 'IS > IIPAR!', 'inptr.f' )
         ENDIF

         ! Make sure JS is not larger than JJPAR
         IF ( JS > JJPAR ) THEN
            CALL ERROR_STOP( 'JS > JJPAR!', 'inptr.f' )
         ENDIF

         ! Make sure LS is not larger than LLPAR
         IF ( LS > LLPAR ) THEN
            CALL ERROR_STOP( 'LS > LLPAR!', 'inptr.f' )
         ENDIF

         ! Make sure MS is either 0 or 1
         IF ( MS < 0 .or. MS > 1 ) THEN 
            CALL ERROR_STOP( 'MS must be 0 or 1!', 'inptr.f' )
         ENDIF

         ! If MS == 0, then we are saving tracer concentrations,
         ! so make sure that NS is not larger than NTRACE!!
         IF ( MS == 0 .and. NS > NTRACE ) THEN
            CALL ERROR_STOP( 'NS > NTRACE!', 'inptr.f' )
         ENDIF

         ! Loop over each level and save station coordinates
         DO L = 1, LS
            ISAVTC(K) = IS
            JSAVTC(K) = JS
            LSAVTC(K) = L
            MSAVTC(K) = MS
            NSAVTC(K) = NS
            
            ! Treat each level as a new station
            ! Halt if we exceed the max # of stations, NNSTA
            K = K + 1
            IF ( K > NNSTA ) THEN
               CALL ERROR_STOP( 'Too many ND48 stations!', 'inptr.f' )
            ENDIF
         ENDDO
      ENDDO

      ! Added the CLOSE command w/in "inptr.f" (bmy, 6/27/02)
      CLOSE( 9 )
      
      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE INPTR       
                                                              
