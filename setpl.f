! $Id: setpl.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
      SUBROUTINE SETPL
!
!******************************************************************************
!  Subroutine SETPL flags the reactions and species which contribute to 
!  production or loss for a given ND65 prodloss diagnostic family.  
!  (ljm, bey, 1999; bmy, 5/1/03)
!
!  NOTES:
!  (1 ) Now references "file_mod.f" and "error_mod.f".  Also now use IOERROR 
!        to trap I/O errors, and ERROR_STOP to stop the run and deallocate
!        all module arrays.  NAMESPEC is now NAMEGAS for SMVGEAR II. Now 
!        uses F90 declaration syntax.  Set NCS = NCSURBAN for now, since we 
!        have defined our GEOS-CHEM mechanism in the urban slot of SMVGEAR II
!        Updated comments.  (bmy, 5/1/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP, GEOS_CHEM_STOP
      USE FILE_MOD,  ONLY : IU_FILE,    IOERROR

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "comode.h"
      
      ! Parameters
      INTEGER, PARAMETER :: MAXPL=100, MAXMEM=10

      ! Local variables
      INTEGER            :: ICOUNT, I, J, INDEX, IOS
      INTEGER            :: K, M, N, NK, NREAC, NPROD, NPOS
      INTEGER            :: IREAC1, IREAC2, IREAC3
      INTEGER            :: IPROD1, IPROD2, IPROD3
      INTEGER            :: NFAMMEM(MAXFAM)
      INTEGER            :: IFAMMEM(MAXMEM,MAXFAM)
      INTEGER            :: ITEMPREAC(NMRPROD) 
      INTEGER            :: NNPL(MAXFAM)
      INTEGER            :: NKPL(MAXPL,MAXFAM)
      INTEGER            :: IPLREAC(NMRPROD,MAXPL,MAXFAM)
      REAL*8             :: PL
      REAL*8             :: COEFMEM(MAXMEM,MAXFAM)
      REAL*8             :: COEFPL(MAXPL,MAXFAM)
      CHARACTER(LEN=5)   :: MOLNAME, FAMNAME, EXTRACHAR
      CHARACTER(LEN=7)   :: JUNK

      !=================================================================
      ! SETPL begins here!
      !=================================================================

      ! Set NCS = NCSURBAN for now, since we have defined our GEOS-CHEM
      ! mechanism in the urban slot of SMVGEAR II. (bmy, 4/21/03)
      NCS = NCSURBAN
      
      ! Initialize
      ICOUNT = 0

      ! Open the "prodloss.dat" file for input
      OPEN( IU_FILE, FILE='prodloss.dat', STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'setpl:1' )

      !=================================================================
      ! Read info from "prodloss.dat"
      !=================================================================
      DO
         
         ! Read line
         READ( IU_FILE, '(a7)', IOSTAT=IOS ) JUNK
         
         ! IOS < 0 is EOF; otherwise it's an I/O error
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'setpl:2' )

         ! A new family is denoted by "*family" in the file
         IF ( JUNK == '*family' ) THEN

            !-------------------------------------------
            ! Increment # of families
            !-------------------------------------------
            ICOUNT = ICOUNT + 1

            ! Make sure we don't exceed MAXFAM families
            IF ( ICOUNT > MAXFAM ) THEN
               CALL ERROR_STOP( 'Too many ND65 families!', 'setpl.f' )
            ENDIF

            !-------------------------------------------
            ! Read family name
            !-------------------------------------------
            READ( IU_FILE, '(22x,a5)', IOSTAT=IOS ) FAMNAME
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'setpl:3' )

            ! Store family index in IFAM
            DO J = 1, NSPEC(NCS)
               IF ( NAMEGAS(J) == FAMNAME ) IFAM(ICOUNT) = J
            ENDDO

            !-------------------------------------------
            ! Read family type: production or loss
            !-------------------------------------------
            READ( IU_FILE, '(22x,a4)', IOSTAT=IOS ) PORL(ICOUNT)
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'setpl:4' )

            ! Convert PORL to lower case if necessary
            IF ( PORL(ICOUNT) == 'PROD' ) PORL(ICOUNT) = 'prod'
            IF ( PORL(ICOUNT) == 'LOSS' ) PORL(ICOUNT) = 'loss'

            !-------------------------------------------
            ! Read number of members for this family
            !-------------------------------------------
            READ( IU_FILE, '(22x,i2)', IOSTAT=IOS ) NFAMMEM(ICOUNT)
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'setpl:5' )

            ! Make sure we don't exceed MAXMEM
            IF ( NFAMMEM(ICOUNT) > MAXMEM ) THEN
               CALL ERROR_STOP( 'Too many family members!', 'setpl.f' )
            ENDIF 

            ! Write to "smv2.log"
            WRITE( IO93, 104 ) 
     &           ICOUNT, FAMNAME, PORL(ICOUNT), NFAMMEM(ICOUNT) 
 104        FORMAT(/, 'Family ', i2, ' is ' ,a5, ' ', a4,
     &                ' with ',  i2, ' members' )

            WRITE( IO93, 105 ) 
 105        FORMAT( 'ind', 2x, 'species', 1x, 'jnum', 2x, 'coef' )

            !-------------------------------------------
            ! Read info for each family member
            !-------------------------------------------
            DO I = 1, NFAMMEM(ICOUNT) 
               READ( IU_FILE, '(22X,A5,f11.0)', IOSTAT=IOS ) 
     &              MOLNAME, COEFMEM(I,ICOUNT)
               IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'setpl:6' )

               ! Store each family  member in IFAMMEM
               DO J = 1, NSPEC(NCS)
                  IF ( NAMEGAS(J) == MOLNAME ) IFAMMEM(I,ICOUNT) = J
               ENDDO

               ! Write to "smv2.log"
               WRITE( IO93, '(i2,3x,a5,2x,i3,2x,f5.1 )') 
     &              ICOUNT,            MOLNAME, 
     &              IFAMMEM(I,ICOUNT), COEFMEM(I,ICOUNT)
            ENDDO

         ENDIF

      ENDDO

      ! Close "prodloss.dat"
      CLOSE( IU_FILE )

      ! Make sure ICOUNT is the same as NFAMILIES
      IF ( ICOUNT /= NFAMILIES ) THEN
         CALL ERROR_STOP( '# of families not as expected', 'setpl.f' )
      ENDIF

      !=================================================================
      ! Now determine which reactions are sources or sinks of the
      ! specified families. Amend the IRM array accordingly.
      !=================================================================
      DO N = 1, NFAMILIES
         NNPL(N) = 0
      ENDDO

      ! Loop over all rxns (NTRATES = # of kinetic + photo rxns)
      DO NK = 1, NTRATES(NCS)

         ! If this rxn hasn't been turned off...
         IF ( LSKIP(NK,NCS) == 0 ) THEN

            ! Index of first reactant
            IREAC1 = IRM(1,NK,NCS)

            ! Index of first product
            IPROD1 = IRM(NPRODLO,NK,NCS)

            ! Skip emission rxns
            IF ( NAMEGAS(IREAC1) == 'EMISSION' ) GOTO 150  

            ! Skip drydep rxns
            DO N = 1, NDRYDEP(NCS)
               IF ( NK ==  NKDRY(N,NCS) ) GOTO 150
            ENDDO

            !===========================================================
            ! For this rxn, loop over all prod/loss diagnostic families
            !===========================================================
            DO N = 1, NFAMILIES

               ! Initialize for each family
               PL        = 0
               NPROD     = 0
               ICOUNT    = 0
               ITEMPREAC = 0

               !========================================================
               ! For each rxn, loop over reactants and products
               ! and compute how many moles are gained and lost 
               !========================================================
               DO I = 1, NPRODHI

                  ! Increment product count (1st 4 slots are reactants)
                  IF ( I > 4 ) ICOUNT = ICOUNT + 1

                  ! Skip blank entries
                  IF ( IRM(I,NK,NCS) /= 0 ) THEN

                     ! Store reactant index for later use
                     ITEMPREAC(I) = IRM(I,NK,NCS)

                     ! Ensure NPROD skips over the reactant slots of IRM
                     IF ( I     > 4      ) NPROD = NPROD + 1
                     IF ( NPROD < ICOUNT ) NPROD = ICOUNT

                     ! Loop over all family members
                     DO J = 1, NFAMMEM(N)

                        ! Test for product or reactant
                        IF ( IRM(I,NK,NCS) == IFAMMEM(J,N) ) THEN

                           !============================================
                           ! PRODUCT: The # of moles that prodloss 
                           ! family N gains is the # of moles that 
                           ! species M contributes to family N (i.e. 
                           ! COEFMEM(J,N) ) times the # of moles of 
                           ! species M gained in the reaction (i.e. 
                           ! FKOEF(I,NK,NCS) ).
                           !============================================
                           IF ( I >= NPRODLO ) THEN
                              PL = PL + COEFMEM(J,N) * FKOEF(I,NK,NCS)
                           ENDIF

                           !============================================
                           ! REACTANT: The # of moles that prodloss 
                           ! family N loses is the # of moles that 
                           ! species M contributes to family N (i.e. 
                           ! COEFMEM(J,N) ).  Here FKOEF is almost 
                           ! always 1 for reactants.
                           !============================================
                           IF ( I < NPRODLO ) THEN 
                              PL = PL - COEFMEM(J,N)
                           ENDIF
                        ENDIF   
                     ENDDO      
                  ENDIF         
               ENDDO

               !========================================================
               ! If there is a production or loss for prodloss family 
               ! N, then update IRM and the other arrays
               !========================================================
               IF ( ( PL > 0 .AND. PORL(N) == 'prod' )  .OR.
     &              ( PL < 0 .AND. PORL(N) == 'loss' ) ) THEN

                  ! # of prod or loss rxns for family N
                  NNPL(N)            = NNPL(N) + 1

                  ! Error check
                  IF ( NNPL(N) .GT. MAXPL ) THEN
                     CALL ERROR_STOP( 'Number of rxns exceeds MAXPL!', 
     &                                'setpl.f' )
                  ENDIF   

                  ! Index of IRM for one beyond the next product
                  NPOS               = NPRODLO + NPROD 

                  ! Store # of each rxn in NKPL for output below
                  NKPL(NNPL(N),N)    = NK

                  ! Store P/L coeff for each rxn in COEFPL for output below
                  COEFPL(NNPL(N),N)  = PL

                  ! Store the family name as the "last" product of the
                  ! of the rxn -- in the (NPRODLO+NPROD)th slot of IRM
                  IRM(NPOS,NK,NCS)   = IFAM(N)

                  ! Also store the total prod/loss of family N 
                  ! in the (NPRODLO+NPROD)th of the FKOEF array
                  FKOEF(NPOS,NK,NCS) = ABS( PL )

                  ! Loop over all reactants and products
                  DO I = 1, NMRPROD

                     ! Zero any negative reactant/product indices
                     IF ( ITEMPREAC(I) < 0 ) ITEMPREAC(I) = 0

                     ! 3-body rxn???
                     IF ( ITEMPREAC(3) > 0 ) THEN
                        WRITE( 6, 1190 ) NK
 1190                   FORMAT( 'SETPL: Problem with rxn # ',i4 )
                        CALL GEOS_CHEM_STOP
                     ENDIF

                     ! Save reactants and products for this
                     ! reaction in IPLREAC for output below
                     IPLREAC(I,NNPL(N),N) = ITEMPREAC(I)
                  ENDDO
               ENDIF    
            ENDDO
         ENDIF

         !-------------------------------
         ! Skip emission & drydep rxns 
         !-------------------------------
 150     CONTINUE 
      ENDDO

      !=================================================================
      ! Write out prod or loss reactions to "smv2.log"
      !=================================================================
      WRITE( IO93, '(/,a)' ) REPEAT( '=', 79 )
      WRITE( IO93, '(a)'   ) 'Here are the prod and loss reactions'
      WRITE( IO93, '(a)'   ) REPEAT( '=', 79 )

      ! Loop over P/L diagnostic families
      DO N = 1, NFAMILIES

         ! Write family header
         WRITE( IO93, 587 ) NAMEGAS(IFAM(N)), PORL(N), NNPL(N)
 587     FORMAT( /, 'Family ',a5,' ',a4,' -- no of rxns is  ',i3, 5x,
     &           'coefficient')

         ! Loop over prod/loss reactions
         DO I = 1, NNPL(N)

            ! Rxn number
            NK        = NKPL(I,N)

            ! Reactant indices
            IREAC1    = IPLREAC(1,I,N)
            IREAC2    = IPLREAC(2,I,N)

            ! Product indices
            IPROD1    = IPLREAC( NPRODLO,   I,N)
            IPROD2    = IPLREAC((NPRODLO+1),I,N)
            IPROD3    = IPLREAC((NPRODLO+2),I,N)

            ! Character to denote 3 or more products
            EXTRACHAR = '     '
            IF ( IPROD3 .GT. 0 ) EXTRACHAR = '+ ...'

            ! Test for kinetic or photo rxns
            IF ( NK .LE. NRATES(NCS) ) THEN

               !----------------------
               ! Write kinetic rxns
               !----------------------
               WRITE(IO93,588) I, NK, NAMEGAS(IREAC1),
     &              NAMEGAS(IREAC2),  NAMEGAS(IPROD1),
     &              NAMEGAS(IPROD2),  EXTRACHAR, COEFPL(I,N)

 588           FORMAT(I3,1X,I3,1X,A5,' + ',A5,' = ',A5,' + ',A5,
     &                A5,1X,ES13.6)

            ELSE

               !----------------------
               ! Write photo rxns
               !----------------------
               WRITE(IO93,589) I, NK, NAMEGAS(IREAC1),
     &              NAMEGAS(IPROD1),  NAMEGAS(IPROD2), 
     &              EXTRACHAR,        COEFPL(I,N)

 589           FORMAT(I3,1X,I3,1X,A5,' +  hv   = ',A5,' + ',A5,
     &                A5,1X,1P1E13.6)

            ENDIF        
            
!###            !### Debug
!###            WRITE( 6, '(i4,1x,16(a,'':'')))' ) 
!###     &           NK, ( TRIM(NAMEGAS(IRM(J,NK,NCS))), J=1,16 )
!###            WRITE( 6, '(i4,1x,4f4.1,''/'',12f4.1)' ) 
!###     &           NK, ( FKOEF(J,NK,NCS), J=1,16 )
         ENDDO             
      ENDDO

      ! Return to calling program
      END SUBROUTINE SETPL
