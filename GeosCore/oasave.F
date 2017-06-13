!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: oasave
!
! !DESCRIPTION: Subroutine OASAVE stores the concentrations of organic aerosols
!  for the ND42 diagnostic and for the timeseries and satellite diagnostics
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OASAVE( SAVE_TMP6, am_I_Root, 
     &                   Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE CMN_SIZE_MOD
      USE COMODE_MOD,      ONLY : AIRDENS, CSPEC, JLOP
      USE COMODE_LOOP_MOD 
      USE TRACERID_MOD,       ONLY : IDTOCPI, IDTOCPO
      ! monoterpene + sesquiterpene SOA (hotp 5/24/10)
      USE TRACERID_MOD,       ONLY : IDTTSOA1, IDTTSOA2
      USE TRACERID_MOD,       ONLY : IDTTSOA3, IDTTSOA0
      ! isoprene SOA (hotp 5/24/10)
      USE TRACERID_MOD,       ONLY : IDTISOA1, IDTISOA2, IDTISOA3
      ! add lumped arom/ivoc (hotp 5/15/10)
      USE TRACERID_MOD,       ONLY : IDTASOAN, IDTASOA1
      USE TRACERID_MOD,       ONLY : IDTASOA2, IDTASOA3
      ! Irreversible uptake of isoprene SOA precursors (eam, 2014):
      USE TRACERID_MOD,       ONLY : IDTSOAIE, IDTSOAME, IDTINDIOL      
      ! semivolpoa: replace OCPO and OCPI with POA (hotp 2/17/09)
      USE TRACERID_MOD,       ONLY : IDTPOA1,  IDTPOA2
      ! semivolpoa4: add OPOA (hotp 3/27/09)
      USE TRACERID_MOD,       ONLY : IDTOPOA1, IDTOPOA2
      ! Glyoxal and methylglyoxal (eam, 2014):
      USE TRACERID_MOD,       ONLY : IDTSOAGX, IDTSOAMG, IDTSOAGC
      ! Low volatility SOA precursor from ISOPOOH (eam, 02/2015):
      USE TRACERID_MOD,       ONLY : IDTLVOCOA
      ! SOA from ISOP+NO3 (eam, 02/2015):
      USE TRACERID_MOD,       ONLY : IDTISN1OA
      !USE CMN_O3_MOD            ! SAVE_TMP6  (total organic aerosol)

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      ! Array of total aerosol concentration [ug/m3]
      REAL*8,   INTENT(OUT) :: SAVE_TMP6(IIPAR,JJPAR,LLPAR)
      INTEGER,  INTENT(OUT) :: RC          ! Success or failure?

!
! !REVISION HISTORY: 
!  10 Jul 2014 - E. A. Marais - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L, JLOOP   
      REAL*8  :: TEMP6
      REAL*8  :: FACTOR

      ! SOAupdate: ratio of OM/OC (hotp 6/12/09)
      REAL*8, PARAMETER    :: OCFPOA  = 1.4d0
      REAL*8, PARAMETER    :: OCFOPOA = 1.4d0*1.5d0  ! 2.1

      ! SOAupdate: Logical SOA tracer flags (mpayer, 7/12/11)
      LOGICAL           :: IS_OC
      LOGICAL           :: IS_SOAGX,  IS_SOAMG,  IS_SOAGC
      LOGICAL           :: IS_TSOA,   IS_ISVOC,  IS_ASOA
      LOGICAL           :: IS_POA,    IS_OPOA     
      LOGICAL           :: IS_SOAIE,  IS_SOAME,  IS_INDIOL
      LOGICAL           :: IS_LVOCOA, IS_ISN1OA
      LOGICAL           :: LISOPOA

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! OASAVE begins here!
      !
      ! Save organic aerosol concentrations (TEMP6 is moved from 
      ! DIAG42_MOD (eam, 2014)
      !=================================================================

      ! SOAupdate: Define logical flags to decide whether or not to add to
      ! array.  This will prevent out-of-bounds errors. 
      ! (mpayer, 7/12/11)
      IS_OC      = ( IDTOCPI  > 0 .AND. IDTOCPO  > 0 )
      IS_SOAGX   = ( IDTSOAGX > 0 )
      IS_SOAMG   = ( IDTSOAMG > 0 )
      IS_SOAGC   = ( IDTSOAGC > 0 )
      IS_TSOA    = ( IDTTSOA1 > 0 .AND. IDTTSOA2 > 0 .AND. IDTTSOA3 > 0
     &         .AND. IDTTSOA0 > 0 )
      IS_ISVOC   = ( IDTISOA1 > 0 .AND. IDTISOA2 > 0 
     &         .AND. IDTISOA3 > 0 )
      IS_ASOA    = ( IDTASOAN > 0 .AND. IDTASOA1 > 0 .AND. IDTASOA2 > 0 
     &         .AND. IDTASOA3 > 0 )
      IS_POA     = ( IDTPOA1  > 0 .AND. IDTPOA2  > 0 )
      IS_OPOA    = ( IDTOPOA1 > 0 .AND. IDTOPOA2 > 0 )
      IS_SOAIE   = ( IDTSOAIE > 0 )
      IS_SOAME   = ( IDTSOAME > 0 )
      IS_INDIOL  = ( IDTINDIOL> 0 )
      IS_LVOCOA  = ( IDTLVOCOA> 0 )
      IS_ISN1OA  = ( IDTISN1OA> 0 )

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Copy fields from INPUT_OPT to local variables for use below
      LISOPOA = Input_Opt%LISOPOA

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, FACTOR )
!$OMP+PRIVATE( TEMP6 )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Conversion factor from [kg] --> [ug/m3]
         FACTOR        = 1d9 / State_Met%AIRVOL(I,J,L)

         ! sum of all OA (hotp 5/24/10)
         ! Non-isoprene semivolatile partitioning terpenes:
         IF ( IS_TSOA .AND. IS_ASOA ) THEN
            TEMP6 = STT(I,J,L,IDTTSOA1 ) +
     &              STT(I,J,L,IDTTSOA2 ) + 
     &              STT(I,J,L,IDTTSOA3 ) + 
     &              STT(I,J,L,IDTTSOA0 ) + 
     &              STT(I,J,L,IDTASOAN ) + 
     &              STT(I,J,L,IDTASOA1 ) + 
     &              STT(I,J,L,IDTASOA2 ) + 
     &              STT(I,J,L,IDTASOA3 )

         ENDIF

         If ( LISOPOA ) THEN

            ! Reactive uptake to form isoprene SOA:
            IF ( IS_SOAIE  ) TEMP6 = TEMP6 + STT(I,J,L,IDTSOAIE  )
            IF ( IS_SOAME  ) TEMP6 = TEMP6 + STT(I,J,L,IDTSOAME  )
            IF ( IS_INDIOL ) TEMP6 = TEMP6 + STT(I,J,L,IDTINDIOL )
            IF ( IS_LVOCOA ) TEMP6 = TEMP6 + STT(I,J,L,IDTLVOCOA )
            IF ( IS_ISN1OA ) TEMP6 = TEMP6 + STT(I,J,L,IDTISN1OA )
            IF ( IS_SOAGX  ) TEMP6 = TEMP6 + STT(I,J,L,IDTSOAGX  )
            IF ( IS_SOAMG  ) TEMP6 = TEMP6 + STT(I,J,L,IDTSOAMG  )
            IF ( IS_SOAGC  ) TEMP6 = TEMP6 + STT(I,J,L,IDTSOAGC  )

         ELSE 

            ! Semivolatile partitioning to form isoprene SOA:
            IF ( IS_ISVOC ) THEN
               TEMP6 = TEMP6 +
     &              STT(I,J,L,IDTISOA1 ) + 
     &              STT(I,J,L,IDTISOA2 ) + 
     &              STT(I,J,L,IDTISOA3 )
            ENDIF

         ENDIF

         IF ( IS_POA ) THEN
            TEMP6 = TEMP6 + STT(I,J,L,IDTPOA1) * OCFPOA +
     &                      STT(I,J,L,IDTPOA2) * OCFPOA
         ENDIF
          
         IF ( IS_OPOA ) THEN
            TEMP6 = TEMP6 + STT(I,J,L,IDTOPOA1) * OCFOPOA +
     &                      STT(I,J,L,IDTOPOA2) * OCFOPOA
         ENDIF

         IF ( IS_OC ) THEN
            TEMP6 = TEMP6 + STT(I,J,L,IDTOCPI) * OCFOPOA +
     &                      STT(I,J,L,IDTOCPO) * OCFOPOA
         ENDIF

         ! Convert from [kg] --> [ug/m3]
         SAVE_TMP6(I,J,L) = TEMP6 * FACTOR

      ENDDO
      ENDDO
      ENDDO

!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE OASAVE
!EOC
