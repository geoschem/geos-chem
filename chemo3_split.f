! $Id: chemo3_split.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      SUBROUTINE CHEMO3_SPLIT
!
!******************************************************************************
!  Subroutine CHEMO3_split performs ozone chemistry with specified production
!  and loss rates for 13 geographical regions. (amf, bey, bmy, 6/9/99, 2/11/03)
! 
!  NOTES:
!  (1 ) This subroutine should be called once every NCHEM minutes.
!  (2 ) Now use IOERROR to trap I/O errors (bmy, 6/9/99)
!  (3 ) Now remove dry deposition losses for NO2 and O3 from Ox.
!        Also DO-loops and IF statements for more efficient execution.
!        (bmy, 11/23/99) 
!  (4 ) LLTROP is now the highest tropospheric level underneath the
!       annual mean tropopause.
!  (5 ) Need to turn on ND65 and ND44 to use 24hr avg chem & drydep 
!        fluxes (amf, 8/9/00)
!  (6 ) Updated comments, cosmetic changes (bmy, 7/2/01)
!  (7 ) Removed hardwired stuff and daily loss files -- amf says that
!        we shouldn't need these anymore. (amf, bmy, 8/28/01)
!  (8 ) Bug fix: now use true exponential loss for drydep, instead of 
!        first-order approximation. (bdf, bmy, 7/11/02)
!  (9 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 7/30/02)
!  (10) Now reference DEPSAV from "drydep_mod.f".  Also remove extraneous
!        dimension from DEPSAV.  Remove reference to "CMN_SAV".(bmy, 11/19/02)
!  (11) Now replace DXYP(J) with routine GET_AREA_M2 from "grid_mod.f"
!        Now use functions GET_TS_CHEM and GET_TAU from "time_mod.f",  Now
!        define parameters for 1x1 grid for quick fix.  (bmy, 2/11/03)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,   ONLY : AD44, AD65, P24H, L24H
      USE DRYDEP_MOD, ONLY : DEPSAV
      USE GRID_MOD,   ONLY : GET_AREA_M2
      USE TIME_MOD,   ONLY : GET_TS_CHEM, GET_TAU

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "CMN"
#     include "CMN_DIAG"
#     include "CMN_O3"
#     include "CMN_SETUP"
#     include "CMN_TIMES" !amf - for saving dryd & chem loss fluxes

!------------------------------------------------------------------------------
!  based on emisi_co.f: Define boundaries for 
!  geographic CO regions for the appropriate grid
!
      INTEGER                :: LDIV(2) = (/ 6, 11 /)

#if   defined( GRID4x5 )
      INTEGER                :: IDIV(5) = (/ 11, 24, 34, 48, 66 /)
      INTEGER                :: SHDIV   = 27
      INTEGER                :: JDIV    = 33
      INTEGER                :: USDIV(2)= (/ 29, 35 /)

#elif defined( GRID2x25 )
      INTEGER                :: IDIV(5) = (/ 22, 47, 67, 95, 131 /)
      INTEGER                :: SHDIV   = 53
      INTEGER                :: JDIV    = 64
      INTEGER                :: USDIV(2)= (/ 58, 70 /)

#elif defined( GRID1x1 )

      ! NOTE: These are not correct, but they will allow the code
      !       to compile for now.  fix these later. (bmy, 3/11/03)
      INTEGER                :: IDIV(5) = (/ 22, 47, 67, 95, 131 /)
      INTEGER                :: SHDIV   = 53
      INTEGER                :: JDIV    = 64
      INTEGER                :: USDIV(2)= (/ 58, 70 /)

#endif
!------------------------------------------------------------------------------


      ! Local variables
      INTEGER                :: I, J, L, N, LN65
      INTEGER                :: LLEV, M

      ! amf: indices for correctly accessing STT array
      INTEGER                :: ISET, JSET, LSET, II, JJ, LL2, NN

      REAL*8                 :: DTCHEM
      REAL*8                 :: PP, PPROD, LL, PL, DTC

      ! External routines
      REAL*8, EXTERNAL       :: BOXVL

      !=================================================================
      ! CHEMO3_SPLIT begins here!
      !=================================================================
      WRITE( 6, '( '' --- CHEMO3_SPLIT - at : '', f10.2 )' ) GET_TAU()

      !### Debug
      !write(222,*) P24H(24:26,30:38,14:16)
      !write(322,*) L24H(20:25,30:38,14:16)
      !write(*,*) 'LLTROP = ', LLTROP
      !CALL FLUSH(6)
      !CALL FLUSH(222)
      !CALL FLUSH(322)

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! calculate offset to correctly access where STT_TEMP51
      ! is in the global STT array:
      ISET = IMIN_AREA - 1
      JSET = JMIN_AREA - 1
      LSET = LMIN_AREA - 1

      !### Debug
      !print*, 'IMIN... ', IMIN_AREA, JMIN_AREA, LMIN_AREA
      !print*, 'ISET, JSET, LSET : ', ISET, JSET, LSET

      !=================================================================
      ! Change tracer array: New STT = Old STT + ( Production - Loss )
      ! 
      ! BOXVL  = Volume of grid box            [ cm3             ]
      ! DTCHEM = Length of chemistry timestep  [ sec             ]
      ! PP     = Production of O3              [ kg/box/timestep ]
      ! LL     = Loss       of O3              [ kg/box/timestep ]
      !
      ! amf -- now add loop over 13 tracers & partition accordingly. 
      !
      ! TRACER 1  = TOTAL O3
      ! TRACER 2  = Upper Trop (above layer 11)
      ! TRACER 3  = Middle Trop (layers 7-11)
      ! TRACER 4  = Southern hemisphere BL (levs 1-6)
      ! TRACER 5  = Pacific BL (1-6)
      ! TRACER 6  = N America BL (1-6)
      ! TRACER 7  = Atlantic BL (1-6)
      ! TRACER 8  = Europe BL (1-6)
      ! TRACER 9  = N Africa BL (1-6)
      ! TRACER 10 = Asia BL (1-6)
      ! TRACER 11 = flux in from stratosphere 
      !              (specified in upbdflx_O3.f)
      ! TRACER 12 = initial O3 read in from restart file 
      !              in read_restart_file.f
      ! TRACER 13 = USA production
      !=================================================================

      DO L = 1, LLTROP
         LL2 = L - LSET
         
         IF ( L > LDIV(2) ) THEN
            LLEV = 1
         ELSE IF (L > LDIV(1) .and. L <= LDIV(2) ) THEN 
            LLEV = 2
         ELSE 
            LLEV = 3
         ENDIF
         
         DO J = 1, JJPAR
            JJ = J - JSET
            DO I = 1, IIPAR
               II = I - ISET

               ! Skip over stratospheric boxes
               IF ( L >= LPAUSE(I,J) ) CYCLE

               PPROD = P24H(I,J,L) * BOXVL(I,J,L) * DTCHEM

               !### Dbug
               !write(555,*) 'BOXVL, DTCHEM', BOXVL(I,J,L), DTCHEM
               !write(555,*)  'I,J,L,P24,L24', I,J,L, P24H(I,J,L), L24H(I,J,L)


               !========================================================
               ! Loop to fill appropriate regions with production rates.
               ! there is never any production in 11 (strat) and 12 
               ! (init cond), so keep PP = 0 for these, and loop only 
               ! through tracers 1-10 and 13. 
               !========================================================
               DO N = 1, NTRACE 
                  PP = 0d0
                  LL = STT(I,J,L,N) * L24H(I,J,L) * 
     &                 BOXVL(I,J,L) * DTCHEM

                  !### Debug
                  !write(*,*) 'I,J,L,N,LLEV', I,J,L,N,LLEV
                  !call flush(6)
                  !write(555,*) 'PP, PROD ', PP, PPROD
                  !write(555,*) 'LL ', LOSS

                  SELECT CASE (N) 
                                
                     ! Total O3, always include 
                     CASE ( 1 )
                        PP = PPROD
                        !write(555,*) 'added to case 1'

                     ! This is UT -  LLEV = 1
                     CASE ( 2 )
                        IF ( LLEV == 1 ) THEN
                           PP = PPROD
                           !write(555,*) 'added to case 2'
                        ENDIF

                     !THis is MT - llev = 2
                     CASE ( 3 ) 
                        IF ( LLEV == 2 )  THEN
                           PP = PPROD
                           !write(555,*) 'added to case 3'
                        ENDIF

                     ! This is all boundary layer so if 
                     ! LLEV /= 3, leave PP at 0d0
                     CASE ( 4:10, 13 ) 
                        IF ( LLEV == 3 ) THEN
                           !write(555,*) 'case(4:10)'

                           !TEST FOR SH -- if so, zero out SH tracer. 
                           IF ( J < SHDIV ) THEN
                              !write(555,*) 'southern hemisphere'

                              ! add to southern hemispheric tracer 
                              IF ( N == 4 ) THEN
                                 PP = PPROD 
                                 !write(555,*) 'added to SH'
                              ENDIF 

                           ! If not southern hemisphere, than in NH 
                           ! so start testing for I values in NH to 
                           ! separate remaining regions
                           ! PACIFIC - zero all other tracers except 
                           ELSE IF ( I >  IDIV(5)  .or. 
     &                               I <= IDIV(1) ) THEN
                              !write(555,*) 'pacific'

                              IF (N .eq. 5) THEN
                                 PP = PPROD 
                                 !write(555,*) 'added to pacific'
                              ENDIF

                           ! North America
                           ELSE IF ( I > IDIV(1)   .and. 
     &                               I <= IDIV(2) ) THEN
                              !write(555,*) 'north america'
                              
                              IF ( N == 6 ) then
                                 PP = PPROD
                                 !write(555,*) 'added to n amer'
                              ENDIF

                              ! Test for USA
                              IF ( J >  USDIV(1)  .and. 
     &                             J <= USDIV(2) ) THEN
                                 IF (N .eq. 13) then
                                    PP = PPROD
                                 ENDIF
                              ENDIF

                           ! Atlantic Ocean
                           ELSE IF ( I >  IDIV(2)  .and. 
     &                               I <= IDIV(3) ) THEN
                              !write(555,*) 'atlantic'

                              IF ( N == 7 ) then
                                 PP = PPROD
                                 !write(555,*) 'added to atlantic'
                              ENDIF
                  
                           ! Europe & N Africa
                           ELSE IF ( I >  IDIV(3)  .and. 
     &                               I <= IDIV(4) ) THEN

                              ! Europe
                              IF ( J >= JDIV ) THEN
                                 !write(555,*) 'EUROPE'

                                 IF ( N == 8 ) then
                                    PP = PPROD 
                                    !write(555,*) 'added to europe'
                                 ENDIF

                              ELSE

                                 ! N Africa
                                 IF (N .eq. 9)  then
                                    PP = PPROD 
                                    !write(555,*) 'added to N africa'
                                 ENDIF
                              ENDIF  !J .lt. JDIV
                
                           !ASIA : I .ge. 49 .and. I .le. 66        
                           ELSE  
                              IF ( N == 10) then
                                 PP = PPROD
                                 !write(555,*) 'added to ASIA'
                              ENDIF
                              
                           ENDIF !end if J NH   
                        ENDIF   ! end if boundary layer
                  END SELECT

                  !=====================================================
                  ! archive prod/loss in ND65
                  !=====================================================
                  IF ( ND65 > 0) THEN

                     ! only record production if production 
                     ! is appropriate for this tracer !
                     IF ( PP > 0d0 ) THEN
!-----------------------------------------------------------------------------
! To save O3 production in ppb/s, uncomment this line:
!                        PL = P24H(I,J,L) * BOXVL(I,J,L) *
!                        TCVV(1)     / AD(I,J,L)
!-----------------------------------------------------------------------------
! To save O3 production in kg/s, uncomment this line:
                        PL = P24H(I,J,L) * BOXVL(I,J,L)
!-----------------------------------------------------------------------------
! To save O3 production in molec/s, uncomment this line:
!                        PL = P24H(I,J,L) * BOXVL(I,J,L) * XNUMOL(1)
!-----------------------------------------------------------------------------
! To save O3 production in molec/cm3/s, uncomment this line:
!                        PL = P24H(I,J,L) * XNUMOL(1)
!-----------------------------------------------------------------------------

                        AD65(I,J,L,N) = AD65(I,J,L,N) + PL

                     ENDIF      !end if production in this box
                              
                     !=====================================================
                     ! Archive 24h avg loss in CHEMPL24 
                     !=====================================================
!-----------------------------------------------------------------------------
! To save O3 loss in ppb/s, uncomment this line:
!                     PL = STT(I,J,L,N) * L24H(I,J,L) * BOXVL(I,J,L) *
!     &                    TCVV(1)      / AD(I,J,L)
!-----------------------------------------------------------------------------
! To save O3 loss in kg/s, uncomment this line
                     PL = STT(I,J,L,N) * L24H(I,J,L) * BOXVL(I,J,L)
!-----------------------------------------------------------------------------
! To save O3 loss in molec/s, uncomment this line:
!                     PL = STT(I,J,L,N) * L24H(I,J,L) *
!     &                    BOXVL(I,J,L) * XNUMOL(1)
!-----------------------------------------------------------------------------
! To save O3 loss in molec/cm3/s, uncomment this line:
!                     PL = STT(I,J,L,N) * L24H(I,J,L) * XNUMOL(1)
!-----------------------------------------------------------------------------

                     AD65(I,J,L,NTRACE+N) = AD65(I,J,L,NTRACE+N) + PL
                  ENDIF

                  !=====================================================
                  ! Perform dry deposition losses and also archive them 
                  ! to the ND44 diagnostic (if necessary).  Drydep 
                  ! losses are in kg/box/sec.
                  !
                  ! M        = index for the DEPSAV array
                  !             amf: 1 is actually Ox if only 1 
                  !             tracer is declared.
                  ! DEPSAV   = Drydep loss frequency of O3 [s^-1] 
                  !             calc in drydep.f
                  ! DTCHEM   = Chemistry time interval [s]
                  !
                  ! Do DRY DEP for level 1 from all tracers: - 
                  ! all tracers deposit !
                  !=====================================================
                  IF ( L == 1 ) THEN
                     IF ( LDRYD ) THEN 

                        ! Save fluxes to ND44 -- 
                        ! NO2 first 12, then O3 (NTRACE should be 12)
                        IF ( ND44 > 0 ) THEN

!----------------------------------------------------------------------------
!TO SAVE In kg/sec.  Don't forget to change units in diag3.f!!
!                           DTC = DEPSAV(I,J,1) * STT(I,J,1,N)
!----------------------------------------------------------------------------
!TO SAVE IN molec/cm^2/s for comparison with full chemistry run, uncomment this
!line:
                           DTC = DEPSAV(I,J,1) * STT(I,J,1,N) *
     &                           XNUMOL(1)     / GET_AREA_M2(J)
!-----------------------------------------------------------------------------
                           AD44(I,J,N,1) = AD44(I,J,N,1) + DTC
                        ENDIF

                        ! Apply drydep loss to STT tracer array --
                        ! Only if drydep is turned on and for 1st layer only
                        STT(I,J,L,1) = STT(I,J,L,1) * 
     &                                 EXP( -DEPSAV(I,J,N) * DTCHEM )

                     ENDIF       
                  ENDIF

                  ! Add appropriate prod/loss for tracer N
                  STT(I,J,L,N) = STT(I,J,L,N) + PP - LL

                  ! Test for negative STTs: - reset to zero.  
                  ! monitor for strat/init cond decay.
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     STT(I,J,L,N) = 0d0
                     WRITE( 110, 110 ) I, J, L, N, GET_TAU()
                     !STOP
                  ENDIF

               ENDDO            !end N tracer loop 
            ENDDO               !end I 
         ENDDO                  !end J loop
      ENDDO                     !end L loop
          
      !### Debug amf 6/9
      !write (444,*) 'SUMS of different tracers'
      !write (444,*) '1: ', sum(STT(:,:,:,1))
      !write (444,*) '2: ', sum(STT(:,:,:,2))
      !write (444,*) '3: ', sum(STT(:,:,:,3))
      !write (444,*) '4: ', sum(STT(:,:,:,4))
      !write (444,*) '5: ', sum(STT(:,:,:,5))
      !write (444,*) '6: ', sum(STT(:,:,:,6))
      !write (444,*) '7: ', sum(STT(:,:,:,7))
      !write (444,*) '8: ', sum(STT(:,:,:,8))
      !write (444,*) '9: ', sum(STT(:,:,:,9))
      !write (444,*) '10: ', sum(STT(:,:,:,10))
      !write (444,*) '11: ', sum(STT(:,:,:,11))
      !write (444,*) '12: ', sum(STT(:,:,:,12))
      !write (444,*) '2-12: ', sum(STT(:,:,:,2:12))
      !call flush(444)
      !STOP 

      !  FORMAT statements
 110  FORMAT( '### CHEMO3: STT < 0 AT (I,J,L,N) = ', 4i5, 
     &        ' at TAU = ', f10.2, ' (before drydep)' )
 120  FORMAT( '### CHEMO3: STT < 0 AT (I,J,L,M) = ', 5i5,
     &        ' at TAU = ', f10.2, ' (after drydep)' )

      ! Return to calling program
      END SUBROUTINE CHEMO3_SPLIT
