#if   defined ( TOMAS )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: aero_drydep
!
! !DESCRIPTION: Subroutine AERO\_DRYDEP removes size-resolved aerosol number 
!  and mass by dry deposition.  The deposition velocities are calcualted from 
!  drydep_mod.f and only aerosol number NK1-NK30 are really treated as dry 
!  depositing tracers while each of the mass species are depositing accordingly
!  with number.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AERO_DRYDEP( am_I_Root, Input_Opt, 
     &                        State_Met, State_Chm, RC )
!
! !USES:
!
      USE DIAG_MOD,           ONLY : AD44
      USE DRYDEP_MOD,         ONLY : DEPSAV, NUMDEP, NTRAIND 
      USE DUST_MOD,           ONLY : SETTLEDUST   
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTNK1, IDTH2SO4
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_STRAT
      USE TOMAS_MOD

      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD

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
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  22 Jul 2007 - Win T.      - Initial version
!  19 Mar 2008 - Win T.      - Calculate removal by gravitational settling for
!                              the internally-mixed aerosol population based
!                              on the same manner implemented by Rokjin Park
!  25 Jan 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  13 Dec 2012 - M. Payer    - Add am_I_Root, Input_Opt, RC as arguments
!  31 May 2013 - R. Yantosca - Now pass Input_Opt, State_Chm to routines
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, L, N, JC, BIN, ID, IJLOOP
      REAL*8                 :: DTCHEM,  AREA_CM2, FLUX
      REAL*8                 :: X, X0(IBINS, ICOMP-IDIAG+1)
      REAL*8                 :: Y, Y0, RKT
      REAL*8                 :: NU0(IIPAR,JJPAR,LLPAR,IBINS)
      REAL*8                 :: DU0(IIPAR,JJPAR,LLPAR,IBINS)
      REAL*8                 :: SIZ_DIA(IIPAR*JJPAR,IBINS)
      REAL*8                 :: SIZ_DEN(IIPAR*JJPAR,IBINS)
      REAL*8                 :: DEN,   DP,    PDP,  TEMP,  P
      REAL*8                 :: CONST, SLIP,  VISC, DELZ,  DELZ1
      REAL*8                 :: TOT1,  TOT2
      REAL*8,  PARAMETER     :: Go = 9.8d0  !Gravitational accel [m/s2]
      REAL*8                 :: TC(LLPAR),  TC0(LLPAR)
      INTEGER,SAVE           :: DRYD(IBINS)
      LOGICAL,SAVE           :: FIRST = .TRUE.
      INTEGER,SAVE           :: H2SO4ID
      LOGICAL                :: DOSETTLING = .TRUE.

      ! Settling velocity of particle (m/s)
      REAL*8                :: VTS(LLPAR)  
      
      !debug
      integer   :: ii, jj , ix, jx, bb, ll
      data ii,jj, ix, jx, bb, ll /55, 29, 55, 29, 1, 1 /

      ! For fields from Input_Opt
      LOGICAL         :: LPRT
      INTEGER         :: N_TRACERS
      REAL*8          :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      REAL*8, POINTER :: STT(:,:,:,:)
      REAL*8, POINTER :: BXHEIGHT(:,:,:)
      REAL*8, POINTER :: T(:,:,:)

      !=================================================================
      ! AERO_DRYDEP begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM    = GET_TS_CHEM() * 60d0

      ! Get logical values from Input_Opt
      LPRT      = ( Input_Opt%LPRT .and. am_I_Root )
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Initialize pointers
      STT      => State_Chm%TRACERS
      BXHEIGHT => State_Met%BXHEIGHT
      T        => State_Met%T

      ! First identify if the size-resolved aerosol tracers have their
      ! deposition velocity calculated.
      ! dryd is an array that keeps the drydep tracer ID.  So if the 
      ! aerosol component has dryd = 0, that means it was not included 
      ! as a dry depositting tracer.
      IF(FIRST) THEN
      DRYD = 0
      DO BIN = 1, IBINS
         DO N = 1, NUMDEP
            !just want to match only once (win, 5/24/06)
            IF ( BIN == 1 .and. NTRAIND(N)==IDTH2SO4) H2SO4ID = N
cdbg            print *, 'DRYDEP Species:',N
            IF ( NTRAIND(N) == ( IDTNK1-1+BIN ) )THEN
cdbg               print *,'Match tracer:',IDTNK1-1+bin,'Bin',bin
               DRYD( BIN ) = N
               GOTO 100
            ENDIF            
         ENDDO
 100     CONTINUE
      ENDDO
      FIRST = .FALSE.
      ENDIF !First


cdbg      print *,'dryd(30)'
cdbg      print *, dryd(:)

      !---------- GRAVITATIONAL SETTLING -------------
      !
      ! First calculate vertical movement and removal by 
      ! gravitational settling
      !
      ! Clarify units:
      !
      !      v_settling = rho   * Dp**2  *  g    *  C 
      !                  -----------------------------
      !                   18    *  visc
      ! [units]
      !         m/s    = kg/m^3 *  m^2   * m/s^2  * -
      !                  -----------------------------
      !                    -    * kg/m/s
      !
      ! NOTES: 
      ! (1 ) Pa s = kg/m/s
      ! (2 ) Slip correction factor is unitless, however, the
      !      equation from Hinds' Aerosol Technology that is 
      !      a function of P and Dp needs the correct units
      !      P [=] kPa and Dp [=] um



      IF ( DOSETTLING ) THEN

      ! SIZ_DIA [=] m  and SIZ_DEN [=] kg/m3
      CALL AERO_DIADEN( SIZ_DIA(:,:), SIZ_DEN(:,:), 1, 
     &                  Input_Opt,    State_Met,    State_Chm )

!$OMP PARALLEL DO
!$OMP+PRIVATE( BIN, I, J, IJLOOP, DP, DEN, CONST, L, P, TEMP )
!$OMP+PRIVATE( PDP, SLIP, VISC, VTS, JC, ID, TC0, TC )
!$OMP+PRIVATE( DELZ, DELZ1, AREA_CM2, TOT1, TOT2, FLUX )
!$OMP+DEFAULT( SHARED )
!$OMP+SCHEDULE( DYNAMIC )
      DO I = 1, IIPAR
      DO J = 1, JJPAR
      DO BIN = 1, IBINS

         ! 1-D grid box index
         IJLOOP      = ( (J-1) * IIPAR ) + I

         DP  = SIZ_DIA(IJLOOP, BIN) * 1.d6 ![=] um
         DEN = SIZ_DEN(IJLOOP, BIN) ![=] kg/m3
         CONST = DEN *  (DP*1.d-6)**2.d0 * Go / 18.d0  

!debug            if(i==ix .and. j==jx .and. bin==bb) then
!               print *, 'L, P, Dp, DEN, SLIP, VISC, VTS(L)'
!            endif
!debug------------

         DO L = 1, LLPAR
            
            ! Get P [kPa], T [K], and P*DP
            P    = GET_PCENTER(I,J,L) * 0.1d0  ![=] kPa
            TEMP = T(I,J,L)          ![=] K
            PDP  = P * DP

            !=====================================================
            ! # air molecule number density
            ! num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
            !
            ! # gas mean free path
            ! lamda = 1.d6 / 
            !     &   ( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
            !
            ! # Slip correction
            ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * 
            !      &  exp( -1.1 * Dp / (2. * lamda))) / Dp
            !=====================================================
            ! NOTE, Slip correction factor calculations following 
            !       Seinfeld, pp464 which is thought to be more 
            !       accurate but more computation required.
            !=====================================================

            ! Slip correction factor as function of (P*dp)
            SLIP = 1d0 + 
     &             ( 15.60d0 + 7.0d0 * EXP(-0.059d0*PDP) ) / PDP
            
            !=====================================================
            ! NOTE, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
            ! which produce slip correction factor with small 
            ! error compared to the above with less computation.
            !=====================================================

            ! Viscosity [Pa s] of air as a function of temp (K)
            ! Sutherland eqn. (ref. pp 25 in Hinds (Aerosol Technology)
            VISC = 1.458d-6 * (TEMP)**(1.5d0) / ( TEMP + 110.4d0 )

            ! Settling velocity [m/s]
            VTS(L) = CONST * SLIP / VISC
!debug            if(i==ix .and. j==jx .and. bin==bb ) then
!               print *, L,P, Dp, DEN, SLIP, VISC, VTS(L)
!            endif
!debug---------

         ENDDO  ! L-loop


         DO JC = 1, ICOMP-IDIAG+1
            ID = IDTNK1 - 1 + BIN + ( IBINS * (JC-1) )

!debug            if(i==ix .and. j==jx .and. l==ll) then
!               write(200,*)'BIN , TC0,  TC, , VTS(L), JC=',JC
!            endif
!debug------------

            ! Method is to solve bidiagonal matrix 
            ! which is implicit and first order accurate in Z
            DO L = 1, LLPAR
               TC0(L) = STT(I,J,L,ID)
               TC(L)  = TC0(L)
            ENDDO

            ! We know the boundary condition at the model top
            L     = LLTROP
            DELZ  = BXHEIGHT(I,J,L)           ![=] meter
            TC(L) = TC(L) / 
     &                   ( 1.d0 + DTCHEM * VTS(L) / DELZ )

            DO L = LLTROP-1, 1, -1
               DELZ  = BXHEIGHT(I,J,L)
               DELZ1 = BXHEIGHT(I,J,L+1)
               TC(L) = 1.d0 / 
     &                       ( 1.d0 + DTCHEM * VTS(L) / DELZ )
     &                     * ( TC(L)  + DTCHEM * VTS(L+1) /
     &                         DELZ1 *  TC(L+1)              )
            ENDDO

            DO L = 1, LLPAR
               STT(I,J,L,ID) = TC(L)

!debug            if(i==ix .and. j==jx .and. l==ll ) 
!     &              print *, BIN, TC0(L), TC(L), VTS(L)              
!debug---------
            ENDDO


            !========================================================      
            ! ND44: Dry deposition diagnostic [#/cm2/s]
            !========================================================
            IF ( ND44 > 0 ) THEN 

               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( I, J, 1 )

               ! Initialize
               TOT1 = 0d0
               TOT2 = 0d0
            
               ! Compute column totals of TCO(:) and TC(:)
               DO L = 1, LLPAR
                  TOT1 = TOT1 + TC0(L)
                  TOT2 = TOT2 + TC(L)
               ENDDO

               ! Convert dust flux from [kg/s] to [#/cm2/s]
               FLUX = ( TOT1 - TOT2 ) / DTCHEM
               FLUX = FLUX * XNUMOL(ID) / AREA_CM2 
           
               ! Save in AD44
               IF( JC == 1 ) THEN 
                  AD44(I,J,DRYD(BIN),1) = AD44(I,J,DRYD(BIN),1) + FLUX
               ELSE
                  AD44(I,J,NUMDEP+BIN+(JC-2)*IBINS,1) =  
     &                 AD44(I,J,NUMDEP+BIN+(JC-2)*IBINS,1) + FLUX
               ENDIF
                 
            ENDIF

         ENDDO ! JC-loop
 


      ENDDO  ! I-loop
      ENDDO  ! J-loop
      ENDDO  ! Bin-loop
!$OMP END PARALLEL DO

      ENDIF  ! DOSETTLING 

      ! Dust gravitational settling 
      IF ( .not. DOSETTLING ) THEN
         CALL SETTLEDUST( am_I_Root, Input_Opt,
     &                    State_Met, State_Chm, RC )
      ENDIF


      !---------- DRY DEPOSITION ----------
      !


      ! Initialize array
      X = 0d0
      X0(:,:) = 0d0

!$OMP PARALLEL DO
!$OMP+PRIVATE( L, J, I, AREA_CM2, RKT, flux, JC, BIN )
!$OMP+PRIVATE( ID, X0, X, Y0, Y )
!$OMP+DEFAULT( SHARED )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLTROP ! Loop over tropospheric grid boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize for safety's sake
         AREA_CM2 = 0d0
         RKT    = 0d0  
         flux   = 0d0

         ! Skip stratospheric boxes
         IF (ITS_IN_THE_STRAT( I, J, L, State_Met ) ) CYCLE

         ! Save the initial 30-bin number and icomp-1 mass component
         DO JC = 1, ICOMP-IDIAG+1
            DO BIN = 1, IBINS
               ID = IDTNK1 - 1 + BIN + ( IBINS * (JC-1) )
               X0(BIN,JC) = STT(I,J,L,ID)
            ENDDO
         ENDDO
!debug            if(i==ii .and. j==jj .and. L==1) 
c     &          print *,'L    STT(',I,J,'L',bb,')   DIF    ',
c     &           'FLUX  AD44' 
c            if(i==ix .and. j==jx .and. L==1) 
c     &          print *,'L    STT(',I,J,'L',bb,')   DIF    ',
c     &           'FLUX  AD44' 
!debug-----
 
         ! Dry deposit 1 aerosol component at a time, start looping from
         ! number and then the icomp-1 mass.
         DO JC = 1, ICOMP-IDIAG+1
         DO BIN = 1, IBINS
            X = 0d0
            ID = IDTNK1 - 1 + BIN + (( JC-1 )* IBINS)

            ! RKT is drydep frequency [1/s] -- PBLFRAC accounts for the 
            ! fraction of each vertical level that is located below the PBL top
            RKT = DEPSAV(I,J,DRYD(BIN))* GET_FRAC_UNDER_PBLTOP( I, J, L)
!debug            if(i==ii .and. j==jj .and. L==1) 
c     &           print *,'JC=',JC,'BIN=',BIN,'ID=',ID,'RKT',RKT
            IF (RKT > 0d0) THEN 
               RKT = RKT * DTCHEM
               ! Remaining amount after drydep
               X  = X0(BIN,JC)* EXP(-RKT)
            ELSE
               X = X0(BIN,JC)
            ENDIF
               
            !==============================================================
            ! ND44 Diagnostic: Drydep flux of bin1,..,bin30 [molec/cm2/s]
            !==============================================================
            IF ( ND44 > 0 .AND. RKT > 0d0 ) THEN

               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( I, J, 1 )
               
               ! Convert from [kg/timestep] to [molec/cm2/s] 
               ! Store in AD44
               FLUX = X0(BIN,JC) - X
               FLUX = FLUX * XNUMOL(ID)/ AREA_CM2/ DTCHEM
               IF ( JC == 1 ) THEN 
                  AD44(I,J,DRYD(BIN),1) = AD44(I,J,DRYD(BIN),1)+ FLUX   
               ELSE
                  AD44(I,J,NUMDEP+BIN+(JC-2)*IBINS,1) =  
     &                 AD44(I,J,NUMDEP+BIN+(JC-2)*IBINS,1) + FLUX
               ENDIF
               
            ENDIF
!debug
c            if(i==ii .and. j==jj .and. 
c     &           bin==bb .and. JC==1) 
c     &           print *,'>',L, STT(I,J,L,ID), X0(BIN,JC) - X, FLUX, 
c     &           AD44(I,J,DRYD(BIN),1)
c            if(i==ii .and. j==jj .and. 
c     &           bin==bb .and. JC==2) 
c     &           print *,'>',L, STT(I,J,L,ID), X0(BIN,JC) - X, FLUX, 
c     &           AD44(I,J,NUMDEP+BIN+(JC-2)*IBINS,1)
!debug            if(i==ix .and. j==jx .and. 
c     &           bin==bb .and. JC==ICOMP) 
c     &          print *, L, STT(I,J,L,ID), X0(BIN,JC) - X, FLUX, 
c     &           AD44(I,J,NUMDEP+BIN+(JC-2)*IBINS,1)
!debug------


            ! Swap X back into STT array
            STT(I,J,L,ID) = X


         ENDDO
         ENDDO
         
         !Dry deposit H2SO4 gas (win, 5/24/06)
         Y0 = STT(I,J,L,IDTH2SO4)
         RKT = DEPSAV(I,J,H2SO4ID)* GET_FRAC_UNDER_PBLTOP( I, J, L)
         Y = Y0 * EXP(-RKT)
         !==============================================================
         ! ND44 Diagnostic: Drydep flux of H2SO4 [molec/cm2/s]
         !==============================================================
         IF ( ND44 > 0 .AND. RKT > 0d0 ) THEN

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )

            ! Convert from [kg/timestep] to [molec/cm2/s] 
            ! Store in AD44
            FLUX = Y0 - Y
            FLUX = FLUX * XNUMOL(IDTH2SO4)/ AREA_CM2/ DTCHEM
            AD44(I,J,H2SO4ID,1) = AD44(I,J,H2SO4ID,1) + FLUX   

         ENDIF

         !Swap final H2SO4 back into STT array
         STT(I,J,L,IDTH2SO4) = Y
         

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


      IF ( LPRT ) PRINT *,'### Finish AERO_DRYDEP'

      ! Free pointers
      NULLIFY( STT      )
      NULLIFY( BXHEIGHT )
      NULLIFY( T        )

      ! Return to calling program
      END SUBROUTINE AERO_DRYDEP
!EOC
#endif
