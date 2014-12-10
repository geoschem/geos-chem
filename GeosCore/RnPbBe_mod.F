!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: RnPbBe_mod
!
! !DESCRIPTION: Module RnPbBe\_MOD contains variables and routines used 
!  for the 222Rn-210Pb-7Be simulation. (hyl, swu, bmy, 6/14/01, 8/4/06)
!\\
!\\
! !INTERFACE: 
!
      MODULE RnPbBe_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC  :: CHEMRnPbBe

      ! NOTE: WE NEED TO KEEP THIS HERE FOR NOW SINCE IT IS USED
      ! BY CODE IN mercury_mod.F.  REMOVE IT LATER. (bmy, 7/7/14)
      PUBLIC  :: SLQ
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Liu,H., D.Jacob, I.Bey, and R.M.Yantosca, Constraints from 210Pb 
!        and 7Be on wet deposition and transport in a global three-dimensional
!        chemical tracer model driven by assimilated meteorological fields, 
!        JGR, 106, D11, 12,109-12,128, 2001.
!  (2 ) Jacob et al.,Evaluation and intercomparison of global atmospheric 
!        transport models using Rn-222 and other short-lived tracers, 
!        JGR, 1997 (102):5953-5970
!  (3 ) Dorothy Koch, JGR 101, D13, 18651, 1996.
!  (4 ) Lal, D., and B. Peters, Cosmic ray produced radioactivity on the 
!        Earth. Handbuch der Physik, 46/2, 551-612, edited by K. Sitte, 
!        Springer-Verlag, New York, 1967. 
!
! !REVISION HISTORY:
!  14 Jun 2001 - H. Liu      - Initial version  
!  (1 ) Added existing routines to this module (bmy, 6/14/01)
!  (2 ) Updated comments (bmy, 9/4/01)
!  (3 ) Eliminate AVGF; redimensioned XTRA2 (bmy, 9/25/01)
!  (4 ) Replace references to PW(I,J) with P(I,J) (bmy, 10/3/01)
!  (5 ) Remove obsolete code from 9/01 and 10/01 (bmy, 10/23/01)
!  (6 ) Removed duplicate variable declarations (bmy, 11/15/01)
!  (7 ) Now read files from DATA_DIR/RnPbBe_200203/ directory.  
!        Also updated comments. (bmy, 3/29/02)
!  (8 ) Incorporated latest changes from Hongyu Liu.  Also split off the
!        code to read in the 7Be emissions into a separate routine. 
!        Add parallel DO-loops in several places.  Cleaned up DRYFLXRnPbBe,
!        and now make sure ND44 accurately represents the drydep fluxes
!        of 210Pb and 7Be. (hyl, bmy, 8/7/02)
!  (9 ) Now reference AD from "dao_mod.f".  Now references "error_mod.f".
!        Moved routine DRYFLXRnPbBe into "drydep_mod.f".  (bmy, 1/27/03)
!  (10) Now references the new "time_mod.f" (bmy, 2/11/03)
!  (11) Bug fix in EMISSRnPbBe -- take abs( lat) for 7Be emiss. (bmy, 6/10/03)
!  (12) Bug fix in EMISSRnPbBe -- shut off 222Rn emissions in polar regions
!        (swu, bmy, 10/28/03)
!  (13) Now references "directory_mod.f", "logical_mod.f", and "tracer_mod.f"
!        (bmy, 7/20/04)
!  (14) Now modified for GCAP and GEOS-5 met fields (swu, bmy, 5/24/05)
!  (15) Now references "tropopause_mod.f"
!  (16) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  19 Nov 2010 - R. Yantosca - Added ProTeX headers
!  08 Nov 2011 - R. Yantosca - Prevent out-of-bounds errors in diagnostics
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now use routines from the new grid_mod.F90
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  07 Jul 2014 - R. Yantosca - Removed routines now orphaned by HEMCO
!  22 Aug 2014 - R. Yantosca - Removed LATSOU, PRESOU, BESOU arrays, these
!                              are now defined in the HEMCO code.
!  22 Aug 2014 - R. Yantosca - Remove XNUMOL_* parameters; these are obsolete
!  04 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC   
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemRnPbBe
!
! !DESCRIPTION: Subroutine CHEMRnPbBe performs loss chemistry on 222Rn, 
!  210Pb, and 7Be.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEMRnPbBe( am_I_Root, Input_Opt, 
     &                       State_Met, State_Chm, RC )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_STRATMESO
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD01
      USE DIAG_MOD,           ONLY : AD02
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACERID_MOD,       ONLY : IDTRn
      USE TRACERID_MOD,       ONLY : IDTPb
      USE TRACERID_MOD,       ONLY : IDTBe7    
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoVal
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
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY:
!  31 Oct 1999 - H. Liu - Initial version
!  (1 ) Now use F90 syntax (bmy, hyl, 3/22/99)
!  (2 ) Add FIRSTCHEM as an argument.  Only compute the exponential terms
!        when FIRSTCHEM = .TRUE., and save the values for later use
!        (bmy, 3/24/99)
!  (3 ) Cosmetic changes (bmy, 10/13/99)
!  (4 ) Eliminate obsolete code and ND63 diagnostic (bmy, 4/12/00)
!  (5 ) Cosmetic changes (bmy, 7/12/00)
!  (6 ) Added to module "RnPbBe_mod.f".  Also updated comments 
!        and made cosmetic changes. (bmy, 6/14/01)
!  (7 ) Add diagnostics for Rn/Be emissions.  Also cleaned up some old code
!        and added parallel DO-loops.  Updated comments. (hyl, 8/6/02)
!  (8 ) Now make FIRSTCHEM a local SAVEd variable.  (bmy, 1/27/03)
!  (9 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 2/11/03)
!  (10) Now references STT and N_TRACERS from "tracer_mod.f" (bmy, 7/20/04)
!  (11) Remove reference to CMN; it's obsolete.  Now use inquiry functions 
!        from "tropopause_mod.f" to diagnose strat boxes. (bmy, 8/15/05)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  08 Nov 2011 - R. Yantosca - Prevent out-of-bounds errors in diagnostics
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  22 Aug 2014 - R. Yantosca - Copy emissions above the PBL to tracer array
!  22 Aug 2014 - R. Yantosca - Cosmetic changes, for clarity
!  04 Sep 2014 - R. Yantosca - Add minor changes for efficiency
!  04 Nov 2014 - M. Yannetti - Changed REAL*8 to REAL(fp) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER             :: I,        J,         L
      INTEGER             :: LSTART,   N,         N_TRACERS            
      REAL(fp)            :: ADD_Pb,   Be_LOST,   DTCHEM   
      REAL(fp)            :: Pb_LOST,  XNUMOL_Rn, XNUMOL_Pb
      REAL(fp)            :: XNUMOL_Be

      ! Static scalars
      LOGICAL,   SAVE     :: FIRSTCHEM = .TRUE.
      REAL(fp),  SAVE     :: EXP_Rn,   EXP_Pb,    EXP_Be

      ! Arrays
      REAL(fp)            :: Rn_LOST(IIPAR,JJPAR,LLPAR)

      ! To calculate decay rate
      REAL(fp)            :: Decay

      ! To get HEMCO values
      REAL*8              :: Emis
      LOGICAL             :: FOUND
!
! !DEFINED PARAMETERS
!
      ! Ratio of molecular weights of 210Pb/222Rn
      REAL(fp), PARAMETER :: Pb_Rn_RATIO = 210e+0_fp / 222e+0_fp

      ! Ln 2
      REAL(fp), PARAMETER :: ln2   = 0.693147181E+00_fp

      ! Lifetime of Radon-222 in days
      REAL(fp), PARAMETER :: RnTau = 3.83E+00_fp

      !=================================================================
      ! CHEMRnPbBe begins here!
      !=================================================================

      ! Assume success
      RC           =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      N_TRACERS    = Input_Opt%N_TRACERS

      ! Chemistry timestep [s]
      DTCHEM       =  GET_TS_CHEM() * 60e+0_fp

      ! Pre-compute exponential terms only on first timestep
      IF ( FIRSTCHEM ) THEN 

         ! Fraction of (222Rn, 210Pb, 7Be) left after radioactive decay
!-----------------------------------------------------------------------------
! Prior to 11/12/14:
! New algorithm from Christoph Keller (bmy, 11/12/14)
!         EXP_Rn    = EXP( -DTCHEM * 2.097e-6_fp  )
!         EXP_Pb    = EXP( -DTCHEM * 9.725e-10_fp ) 
!         EXP_Be    = EXP( -DTCHEM * 1.506e-7_fp  )
!-----------------------------------------------------------------------------
         Decay     = ln2 / ( RnTau * 24.E+00_fp * 3600.E+00_fp ) 
         EXP_Rn    = EXP( -DTCHEM * Decay        )
         EXP_Pb    = EXP( -DTCHEM * 9.725E-10_fp ) 
         EXP_Be    = EXP( -DTCHEM * 1.506E-7_fp  )

         ! Reset FIRSTCHEM flag
         FIRSTCHEM = .FALSE.

         ! testing only 
         IF(am_I_Root) THEN
            write(*,*) ''
            write(*,*) '### GEOS-Chem Radon simulation ###'
            write(*,*) '    Timestep (secs)   : ', DTCHEM
            write(*,*) '    Rn lifetime (days): ', RnTau
            write(*,*) '    Rn decadence      : ', EXP_Rn
            write(*,*) ''
         ENDIF

      ENDIF

      !=================================================================
      ! Add emissions above the PBL to the tracer array. 
      ! For more information about emissions and PBL mixing, see:
      ! http://wiki.geos-chem.org/Distributing_emissions_in_the_PBL
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, Emis, FOUND )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, N_TRACERS

         ! Skip Pb210, it's source is from Rn decay
         IF ( N == IDTPb ) CYCLE

         ! Loop over all grid boxes
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! If we are using the non-local PBL mixing scheme, then
            ! only add emissions above the PBL top.  All wmissions within
            ! the PBL are handled in routine VDIFFDR in vdiff_mod.F90.
            ! Skip if we have PBL mixing completely turned off, or if
            ! we are using the TURBDAY full-mixing scheme.
            IF ( Input_Opt%LNLPBL ) THEN
               IF ( L < State_Met%PBL_TOP_L(I,J) + 1 ) CYCLE
            ENDIF

            ! Add emissions into the tracer array
            ! Convert from [kg/m2/s] to [kg]
            ! Now use GetHcoVal (ckeller, 10/21/2014)
            CALL GetHcoVal( N, I, J, L, FOUND, Emis=Emis )
            IF ( .NOT. FOUND ) CYCLE
            State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)
     &                                 + ( Emis
     &                                 *   State_Met%AREA_M2  (I,J,L  )
     &                                 *   DTCHEM                      )

         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Radioactive decay of Rn222
      !=================================================================
      
      ! Tracer index
      N = IDTRn

      ! Make sure Rn222 is a defined tracer
      IF ( N > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Rn_LOST = amount of 222Rn lost to decay [kg]
            Rn_LOST(I,J,L)             = State_Chm%TRACERS(I,J,L,N) 
     &                                 * ( 1e+0_fp - EXP_Rn )

            ! ND02 diag: 222Rn lost to decay [kg/s]
            IF ( ND02 > 0 .and. L <= LD02 ) THEN
               AD02(I,J,L,1)           = AD02(I,J,L,1) 
     &                                 + ( Rn_LOST(I,J,L) / DTCHEM )
            ENDIF

            ! Subtract Rn_LOST from STT [kg]
            State_Chm%TRACERS(I,J,L,N) = 
     &      State_Chm%TRACERS(I,J,L,N) - Rn_LOST(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !=================================================================
      ! Radioactive decay of Pb210
      !=================================================================
      
      ! Tracer index
      N = IDTPb

      ! Make sure Pb210 is a defined tracer
      IF ( N > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, ADD_Pb, Pb_LOST )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
           
            ! ADD_Pb = Amount of 210Pb gained by decay from 222Rn [kg]
            ADD_Pb                     = Rn_LOST(I,J,L) * Pb_Rn_RATIO 

#if defined( GCAP )
            !%%%%% FOR GCAP SIMULATION: Divide the 210Pb source by 3.5
            !%%%%% to correct for the strat-trop exchange.  A similar
            !%%%%% line of code is done in the HEMCO extension where
            !%%%%% the Be7 source is computed.  This replicates the
            !%%%%% pre-HEMCO code. (bmy, 9/4/14)
            IF ( ITS_IN_THE_STRATMESO( I, J, L, State_Met ) ) THEN
               ADD_Pb                  = ADD_Pb / 3.5e+0_fp
            ENDIF
#endif

            ! ND01 diag: 210Pb emission from 222Rn decay [kg/s]
            IF ( ND01 > 0 .and. L <= LD01 ) THEN
               AD01(I,J,L)             = AD01(I,J,L) + 
     &                                   ( ADD_Pb / DTCHEM )
            ENDIF

            ! Add 210Pb gained by decay from 222Rn into STT [kg]
            State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) 
     &                                 + ADD_Pb          

            ! Amount of 210Pb lost to radioactive decay [kg]
            ! NOTE: we've already added in the 210Pb gained from 222Rn
            Pb_LOST                    = State_Chm%TRACERS(I,J,L,N) 
     &                                 * ( 1e+0_fp - EXP_Pb )

            ! ND02 diag: 210Pb lost to decay [kg/s]
            IF ( ND02 > 0 .and. L <= LD02 ) THEN
               AD02(I,J,L,2)           = AD02(I,J,L,2) 
     &                                 + ( Pb_LOST / DTCHEM )
            ENDIF

            ! Subtract 210Pb lost to decay from STT [kg]
            State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) 
     &                                 - Pb_LOST
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !=================================================================
      ! Radioactive decay of Be7
      !=================================================================

      ! Tracer index
      N = IDTBe7

      ! Make sure Be7 is a defined tracer
      IF ( N > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, Be_LOST )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Amount of 7Be lost to decay [kg]
            Be_LOST                    = State_Chm%TRACERS(I,J,L,N) 
     &                                 * ( 1d0 - EXP_Be )

            ! ND02 diag: 7Be lost to decay [kg/s]
            IF ( ND02 > 0 .and. L <= LD02 ) THEN
               AD02(I,J,L,3)           = AD02(I,J,L,3) 
     &                                 + ( Be_LOST / DTCHEM )
            ENDIF

            ! Subtract amount of 7Be lost to decay from STT [kg]
            State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) 
     &                                 - Be_LOST
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      END SUBROUTINE CHEMRnPbBe
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: slq
!
! !DESCRIPTION: Subroutine SLQ is an interpolation subroutine from a 
!  Chinese reference book (says Hongyu Liu).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SLQ( X, Y, Z, N, M, U, V, W )
!
! !INPUT PARAMETERS: 
!
      INTEGER :: N        ! First dimension of Z
      INTEGER :: M        ! Second dimension of Z
      REAL(fp)  :: X(N)     ! X-axis coordinate on original grid
      REAL(fp)  :: Y(M)     ! Y-axis coordinate on original grid
      REAL(fp)  :: Z(N,M)   ! Array of data on original grid
      REAL(fp)  :: U        ! X-axis coordinate for desired interpolated value
      REAL(fp)  :: V        ! Y-axis coordinate for desired interpolated value
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)  :: W        ! Interpolated value of Z array, at coords (U,V) 
!
! !REMARKS:
! 
! 
! !REVISION HISTORY: 
!  17 Mar 1998 - H. Liu      - Initial version
!  (1 ) Added to "RnPbBe_mod.f" (bmy, 7/16/01)
!  (2 ) Removed duplicate definition of IQ.  Added comments. (bmy, 11/15/01)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  04 Nov 2014 - M. Yannetti - Changed REAL*8 to REAL(fp) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)  :: B(3), HH
      INTEGER :: NN,   IP, I, J, L, IQ, K, MM

      !=================================================================
      ! SLQ begins here!
      !=================================================================
      NN=3
      IF(N.LE.3) THEN
         IP=1
         NN=N
      ELSE IF (U.LE.X(2)) THEN
         IP=1
      ELSE IF (U.GE.X(N-1)) THEN
         IP=N-2
      ELSE
         I=1
         J=N
 10      IF (IABS(I-J).NE.1) THEN
            L=(I+J)/2
            IF (U.LT.X(L)) THEN
               J=L
            ELSE
               I=L
            END IF
            GOTO 10
         END IF
         IF (ABS(U-X(I)).LT.ABS(U-X(J))) THEN
            IP=I-1
         ELSE
            IP=I
         END IF
      END IF
      MM=3
      IF (M.LE.3) THEN
         IQ=1
         MM=N
      ELSE IF (V.LE.Y(2)) THEN
         IQ=1
      ELSE IF (V.GE.Y(M-1)) THEN
         IQ=M-2
      ELSE
         I=1
         J=M
 20      IF (IABS(J-I).NE.1) THEN
            L=(I+J)/2
            IF (V.LT.Y(L)) THEN
               J=L
            ELSE
               I=L
            END IF
            GOTO 20
         END IF
         IF (ABS(V-Y(I)).LT.ABS(V-Y(J))) THEN
            IQ=I-1
         ELSE
            IQ=I
         END IF
      END IF
      DO 50 I=1,NN
         B(I)=0.0
         DO 40 J=1,MM
            HH=Z(IP+I-1,IQ+J-1)
            DO 30 K=1,MM
               IF (K.NE.J) THEN
                  HH=HH*(V-Y(IQ+K-1))/(Y(IQ+J-1)-Y(IQ+K-1))
               END IF
 30         CONTINUE
            B(I)=B(I)+HH
 40      CONTINUE
 50   CONTINUE
      W=0.0
      DO 70 I=1,NN
         HH=B(I)
         DO 60 J=1,NN
            IF (J.NE.I) THEN
               HH=HH*(U-X(IP+J-1))/(X(IP+I-1)-X(IP+J-1))
            END IF
 60      CONTINUE
         W=W+HH
 70   CONTINUE

      END SUBROUTINE SLQ
!EOC
      END MODULE RnPbBe_MOD


