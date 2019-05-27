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
!  depositing species while each of the mass species are depositing accordingly
!  with number.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AERO_DRYDEP( am_I_Root, Input_Opt,  State_Met,
     &                        State_Chm, State_Diag, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DUST_MOD,           ONLY : SETTLEDUST   
      USE ErrCode_Mod
      USE ERROR_MOD
      USE GC_GRID_MOD,        ONLY : GET_AREA_CM2
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants,      ONLY : g0
      USE PhysConstants,      ONLY : AVO
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE PRECISION_MOD
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TOMAS_MOD
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
      USE DIAG_MOD,           ONLY : AD44
#endif

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
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
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
!  05 Nov 2014 - M. Yannetti - PRECISION_MOD convert REAL*8 to REAL(fp)
!  26 Feb 2015 - E. Lundgren - Replace GET_PCENTER with State_Met%PMID.
!                              Remove dependency on pressure_mod.
!  12 Aug 2015 - E. Lundgren - Now pass am_I_Root and RC to AERO_DIADEN 
!                              to enable unit conversion if needed
!  12 May 2016 - M. Sulprizio- Remove IJLOOP and changed SIZ_DIA and SIZ_DEN
!                              array dimensions from (IIPAR*JJPAR,IBINS) to
!                              (IIPAR,JJPAR,IBINS)
!  07 Jun 2016 - E. Lundgren - Remove XNUMOL; use species database MW instead
!  16 Jun 2016 - R. Silvern  - Now use Ind_ function to define species IDs
!  17 Jun 2016 - R. Yantosca - Now define species ID's on the first call
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!  04 Jan 2018 - M. Sulprizio- Now accept State_Diag as an argument
!  21 Sep 2018 - P.F. Liu    - Fix the bug in unit conversion
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      ! SAVEd scalars 
      LOGICAL,  SAVE     :: DOSETTLING = .TRUE.
      LOGICAL,  SAVE     :: FIRST      = .TRUE.
      INTEGER,  SAVE     :: H2SO4ID
      INTEGER,  SAVE     :: id_H2SO4
      INTEGER,  SAVE     :: id_NK1

      ! Scalars
      LOGICAL            :: LPRT
      INTEGER            :: nDryDep
      INTEGER            :: I,      J,        L
      INTEGER            :: N,      JC,       BIN,   ID
      REAL(fp)           :: DTCHEM, AREA_CM2, FLUX,  X,    Y
      REAL(fp)           :: Y0,     RKT,      DEN,   DP,   PDP
      REAL(fp)           :: TEMP,   P,        CONST, SLIP, VISC
      REAL(fp)           :: DELZ,   DELZ1,    TOT1,  TOT2

      ! Strings
      CHARACTER(LEN=255) :: LOC, MSG

      ! Arrays
      REAL(fp)           :: TC     (            LLPAR      )
      REAL(fp)           :: TC0    (            LLPAR      )
      REAL(fp)           :: VTS    (            LLPAR      ) ! Settling V [m/s]
      REAL(fp)           :: NU0    (IIPAR,JJPAR,LLPAR,IBINS)
      REAL(fp)           :: DU0    (IIPAR,JJPAR,LLPAR,IBINS)
      REAL(fp)           :: SIZ_DIA(IIPAR,JJPAR,      IBINS)
      REAL(fp)           :: SIZ_DEN(IIPAR,JJPAR,      IBINS)
      REAL(fp)           :: X0     (IBINS,ICOMP-IDIAG+1    )

      ! Pointers
      REAL(fp), POINTER  :: Spc     (:,:,:,:)
      REAL(fp), POINTER  :: BXHEIGHT(:,:,:  )
      REAL(fp), POINTER  :: T       (:,:,:  )
      REAL(fp), POINTER  :: DEPSAV  (:,:,:  )                ! IM, JM, nDryDep

      ! SAVEd arrays
      INTEGER,  SAVE     :: DRYD(IBINS)

!      !debug
!     integer   :: ii, jj , ix, jx, bb, ll
!     data ii,jj, ix, jx, bb, ll /55, 29, 55, 29, 1, 1 /

      !=================================================================
      ! AERO_DRYDEP begins here!
      !=================================================================

      ! Assume success
      RC        = GC_SUCCESS

      ! Number of dry-deposited species
      nDryDep   = State_Chm%nDryDep

      ! Check that species units are in [kg] (ewl, 8/13/15)
      IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
         MSG = 'Incorrect species units: ' // TRIM(State_Chm%Spc_Units)
         LOC = 'Routine AERO_DRYDEP in aero_drydep.F'
         CALL GC_Error( MSG, RC, LOC )
      ENDIF

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM    = GET_TS_CHEM()

      ! Get logical values from Input_Opt
      LPRT      = ( Input_Opt%LPRT .and. am_I_Root )

      ! Initialize pointers
      Spc      => State_Chm%Species
      BXHEIGHT => State_Met%BXHEIGHT
      T        => State_Met%T
      DEPSAV   => State_Chm%DryDepSav

      ! First-time setup
      IF ( FIRST ) THEN

         ! Define species ID flags
         id_H2SO4 = Ind_('H2SO4')
         id_NK1   = Ind_('NK1'  )

         ! Make sure species are defined
         IF ( id_H2SO4 < 0 ) THEN
            MSG = 'H2SO4 is not a defined species!'
            LOC = 'Routine AERO_DRYDEP in aero_drydep.F'
            CALL ERROR_STOP( MSG, LOC )
         ENDIF
         IF ( id_NK1 < 0 ) THEN
            MSG = 'NK1 is not a defined species!'
            LOC = 'Routine AERO_DRYDEP in aero_drydep.F'
            CALL ERROR_STOP( MSG, LOC )
         ENDIF

         ! First identify if the size-resolved aerosol species have their
         ! deposition velocity calculated.
         ! dryd is an array that keeps the drydep species ID.  So if the 
         ! aerosol component has dryd = 0, that means it was not included 
         ! as a dry depositting species.
         DRYD = 0
         DO BIN = 1, IBINS
            DO N   = 1, nDryDep
               !just want to match only once (win, 5/24/06)
               IF ( BIN == 1 .and.
     &              State_Chm%Map_DryDep(N)==id_H2SO4 ) THEN
                  H2SO4ID = N
cdbg              print *, 'DRYDEP Species:',N
               ENDIF
               IF ( State_Chm%Map_DryDep(N) == ( id_NK1-1+BIN ) )THEN
cdbg              print *,'Match species:',IDTNK1-1+bin,'Bin',bin
                  DRYD( BIN ) = N
                  GOTO 100
               ENDIF
            ENDDO
 100        CONTINUE
         ENDDO

         ! Reset first-time falg
         FIRST = .FALSE.

      ENDIF

cdbg  print *,'dryd(30)'
cdbg  print *, dryd(:)

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
      CALL AERO_DIADEN( am_I_Root, 1,       Input_Opt, State_Met,
     &                  State_Chm, SIZ_DIA, SIZ_DEN,   RC         )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL ERROR_STOP('CALL AERO_DIADEN',
     &                   'AERO_DRYDEP in aero_drydep.F')
      ENDIF

!$OMP PARALLEL DO
!$OMP+PRIVATE( BIN, I, J, DP, DEN, CONST, L, P, TEMP )
!$OMP+PRIVATE( PDP, SLIP, VISC, VTS, JC, ID, TC0, TC )
!$OMP+PRIVATE( DELZ, DELZ1, AREA_CM2, TOT1, TOT2, FLUX )
!$OMP+DEFAULT( SHARED )
!$OMP+SCHEDULE( DYNAMIC )
      DO I = 1, IIPAR
      DO J = 1, JJPAR
      DO BIN = 1, IBINS

         DP    = SIZ_DIA(I,J,BIN) * 1.d6 ![=] um
         DEN   = SIZ_DEN(I,J,BIN) ![=] kg/m3
         CONST = DEN *  (DP*1.d-6)**2.d0 * g0 / 18.d0  

!debug            if(i==ix .and. j==jx .and. bin==bb) then
!               print *, 'L, P, Dp, DEN, SLIP, VISC, VTS(L)'
!            endif
!debug------------

         DO L = 1, LLPAR
            
            ! Get P [kPa], T [K], and P*DP
            ! Use moist pressure for mean free path (ewl, 3/2/2015)
            P    = State_Met%PMID(I,J,L) * 0.1d0  ![=] kPa
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
            ID = id_NK1 - 1 + BIN + ( IBINS * (JC-1) )

!debug            if(i==ix .and. j==jx .and. l==ll) then
!               write(200,*)'BIN , TC0,  TC, , VTS(L), JC=',JC
!            endif
!debug------------

            ! Method is to solve bidiagonal matrix 
            ! which is implicit and first order accurate in Z
            DO L = 1, LLPAR
               TC0(L) = Spc(I,J,L,ID)
               TC(L)  = TC0(L)
            ENDDO

            ! We know the boundary condition at the model top
            L     = LLCHEM
            DELZ  = BXHEIGHT(I,J,L)           ![=] meter
            TC(L) = TC(L) / ( 1.d0 + DTCHEM * VTS(L) / DELZ )

            DO L = LLCHEM-1, 1, -1
               DELZ  = BXHEIGHT(I,J,L)
               DELZ1 = BXHEIGHT(I,J,L+1)
               TC(L) = 1.d0 / 
     &               ( 1.d0   + DTCHEM * VTS(L)   / DELZ ) *
     &               ( TC(L)  + DTCHEM * VTS(L+1) / DELZ1  *  TC(L+1) )
            ENDDO

            DO L = 1, LLPAR
               Spc(I,J,L,ID) = TC(L)

!debug            if(i==ix .and. j==jx .and. l==ll ) 
!     &              print *, BIN, TC0(L), TC(L), VTS(L)              
!debug---------
            ENDDO


#if defined( BPCH_DIAG ) 
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
               FLUX = FLUX * State_Chm%SpcData(ID)%Info%emMW_g * 
     &                1.e-3_fp / AREA_CM2 
           
               ! Save in AD44
               IF( JC == 1 ) THEN 
                  AD44(I,J,DRYD(BIN),1) = AD44(I,J,DRYD(BIN),1) + FLUX
               ELSE
                  AD44(I,J,nDryDep+BIN+(JC-2)*IBINS,1) =  
     &                 AD44(I,J,nDryDep+BIN+(JC-2)*IBINS,1) + FLUX
               ENDIF
                 
            ENDIF
#endif

         ENDDO ! JC-loop
 


      ENDDO  ! I-loop
      ENDDO  ! J-loop
      ENDDO  ! Bin-loop
!$OMP END PARALLEL DO

      ENDIF  ! DOSETTLING 

      ! Dust gravitational settling 
      IF ( .not. DOSETTLING ) THEN
         CALL SETTLEDUST( am_I_Root, Input_Opt,  State_Met,
     &                    State_Chm, State_Diag, RC )
      ENDIF


      !---------- DRY DEPOSITION ----------
      !


      ! Initialize array
      X = 0d0
      X0(:,:) = 0d0

      ! Loop over chemically-active grid boxes
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%% TEMPORARY FIX: REVERSE ORDER OF LOOPS IN ORDER TO PASS UNIT TESTS    %%%
! %%%                                                                      %%%
! %%% Sal Farina wrote: Change the loop order from LJI to IJL or JIL.      %%%
! %%% This will make the MP code add up the diagnostic in the same order   %%%
! %%% as SP mode (only the outermost loop gets parallelized). Yes looping  %%%
! %%% over LJI should be faster than IJL, but (and correct me if i'm       %%%
! %%% wrong) if we are looping over species inside that loop, all benefits %%%
! %%% are totally lost anyway. The way Spc / STT is defined, tracerid      %%%
! %%% should always be the outermost loop...                               %%%
! %$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!$OMP PARALLEL DO
!$OMP+PRIVATE( L, J, I, AREA_CM2, RKT, flux, JC, BIN )
!$OMP+PRIVATE( ID, X0, X, Y0, Y )
!$OMP+DEFAULT( SHARED )
!$OMP+SCHEDULE( DYNAMIC )
      DO I = 1, IIPAR
      DO J = 1, JJPAR
      DO L = 1, LLCHEM 

         ! Initialize for safety's sake
         AREA_CM2 = 0d0
         RKT    = 0d0  
         flux   = 0d0

         ! Skip non-chemsitry boxes
         IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

         ! Save the initial 30-bin number and icomp-1 mass component
         DO JC = 1, ICOMP-IDIAG+1
            DO BIN = 1, IBINS
               ID = id_NK1 - 1 + BIN + ( IBINS * (JC-1) )
               X0(BIN,JC) = Spc(I,J,L,ID)
            ENDDO
         ENDDO
!debug            if(i==ii .and. j==jj .and. L==1) 
c     &          print *,'L    Spc(',I,J,'L',bb,')   DIF    ',
c     &           'FLUX  AD44' 
c            if(i==ix .and. j==jx .and. L==1) 
c     &          print *,'L    Spc(',I,J,'L',bb,')   DIF    ',
c     &           'FLUX  AD44' 
!debug-----
 
         ! Dry deposit 1 aerosol component at a time, start looping from
         ! number and then the icomp-1 mass.
         DO JC = 1, ICOMP-IDIAG+1
         DO BIN = 1, IBINS
            X = 0d0
            ID = id_NK1 - 1 + BIN + (( JC-1 )* IBINS)

! ******************************************************************************
! NOTE: I'm not sure if this is now covered by dry-deposition in mixing_mod.F90
! (ckeller, 3/5/15)
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
! ******************************************************************************
               
#if defined( BPCH_DIAG )
            !==============================================================
            ! ND44 Diagnostic: Drydep flux of bin1,..,bin30 [molec/cm2/s]
            !==============================================================
            IF ( ND44 > 0 .AND. RKT > 0d0 ) THEN

               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( I, J, 1 )
               
               ! Convert from [kg/timestep] to [molec/cm2/s] 
               ! Store in AD44
               FLUX = X0(BIN,JC) - X
               FLUX = FLUX / (State_Chm%SpcData(ID)%Info%emMW_g * 
     &                1.e-3_fp) / AREA_CM2 / DTCHEM * AVO

               IF ( JC == 1 ) THEN 
                  AD44(I,J,DRYD(BIN),1) = AD44(I,J,DRYD(BIN),1)+ FLUX   
               ELSE
                  AD44(I,J,nDryDep+BIN+(JC-2)*IBINS,1) =  
     &                 AD44(I,J,nDryDep+BIN+(JC-2)*IBINS,1) + FLUX
               ENDIF
               
            ENDIF
!debug
c            if(i==ii .and. j==jj .and. 
c     &           bin==bb .and. JC==1) 
c     &           print *,'>',L, Spc(I,J,L,ID), X0(BIN,JC) - X, FLUX, 
c     &           AD44(I,J,DRYD(BIN),1)
c            if(i==ii .and. j==jj .and. 
c     &           bin==bb .and. JC==2) 
c     &           print *,'>',L, Spc(I,J,L,ID), X0(BIN,JC) - X, FLUX, 
c     &           AD44(I,J,nDryDep+BIN+(JC-2)*IBINS,1)
!debug            if(i==ix .and. j==jx .and. 
c     &           bin==bb .and. JC==ICOMP) 
c     &          print *, L, Spc(I,J,L,ID), X0(BIN,JC) - X, FLUX, 
c     &           AD44(I,J,nDryDep+BIN+(JC-2)*IBINS,1)
!debug------
#endif


            ! Swap X back into Spc array
            Spc(I,J,L,ID) = X


         ENDDO
         ENDDO
         
! ******************************************************************************
! NOTE: I'm not sure if this is now covered by dry-deposition in mixing_mod.F90
! (ckeller, 3/5/15)
         !Dry deposit H2SO4 gas (win, 5/24/06)
         Y0 = Spc(I,J,L,id_H2SO4)
         RKT = DEPSAV(I,J,H2SO4ID)* GET_FRAC_UNDER_PBLTOP( I, J, L)
         Y = Y0 * EXP(-RKT)
! ******************************************************************************
#if defined( BPCH_DIAG )
         !==============================================================
         ! ND44 Diagnostic: Drydep flux of H2SO4 [molec/cm2/s]
         !==============================================================
         IF ( ND44 > 0 .AND. RKT > 0d0 ) THEN

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )

            ! Convert from [kg/timestep] to [molec/cm2/s] 
            ! Store in AD44
            FLUX = Y0 - Y
            FLUX = FLUX / (State_Chm%SpcData(id_H2SO4)%Info%emMW_g *
     &             1.e-3_fp) / AREA_CM2 / DTCHEM * AVO

            AD44(I,J,H2SO4ID,1) = AD44(I,J,H2SO4ID,1) + FLUX   

         ENDIF
#endif

         !Swap final H2SO4 back into Spc array
         Spc(I,J,L,id_H2SO4) = Y

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF ( LPRT ) PRINT *,'### Finish AERO_DRYDEP'

      ! Free pointers
      Spc      => NULL()
      BXHEIGHT => NULL()
      T        => NULL()
      DEPSAV   => NULL()

      ! Check that species units are still in [kg] (ewl, 8/13/15)
      IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
         MSG = 'Incorrect species units at end of routine: ' 
     &            // TRIM(State_Chm%Spc_Units)
         LOC = 'Routine AERO_DRYDEP in aero_drydep.F'
         CALL GC_Error( MSG, RC, LOC )
      ENDIF

      ! Return to calling program
      END SUBROUTINE AERO_DRYDEP
!EOC
#endif
