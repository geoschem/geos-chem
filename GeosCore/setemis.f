! $Id: setemis.f,v 1.3 2009/11/05 15:35:30 phs Exp $
      SUBROUTINE SETEMIS( EMISRR, EMISRRN )
!
!******************************************************************************
!  Subroutine SETEMIS places emissions computed from GEOS-Chem
!  subroutines into arrays for SMVGEAR II chemistry. 
!  (lwh, jyl, gmg, djj, bdf, bmy, 6/8/98, 6/11/08)
!
!  SETEMIS converts from units of [molec tracer/box/s] to units of
!  [molec chemical species/cm3/s], and stores in the REMIS array.  For
!  hydrocarbons that are carried through the GEOS-CHEM model as [molec C], 
!  these are converted back to [molec hydrocarbon], and then stored in REMIS.  
!
!  REMIS(JLOOP,N) = emis. rate of species corr. to tracer N in box JLOOP
!                   (reaction number NTEMIS(N))
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) EMISRR   (REAL*8 ) : CO, hydrocarbon emission   [molec tracer/box/s ]
!  (2 ) EMISRRN  (REAL*8 ) : Multi-level NOx emissions  [molec NOx/box/s    ]
!
!  Variables taken from F90 Modules:
!  ============================================================================
!  (1 ) BIOFUEL  (REAL*8 ) : Biofuel burning emissions  [molec (C)/cm3/s    ]
!  (2 ) BFTRACE  (INTEGER) : Index array for biofuels   [CTM tracer #       ]
!  (3 ) NBFTRACE (INTEGER) : Number of biofuel species  [unitless           ]
!  (4 ) BURNEMIS (REAL*8 ) : Biomass burning emissions  [molec (C)/cm3/s    ] 
!  (5 ) BIOTRCE  (INTEGER) : Index array for bioburn    [CTM tracer #       ] 
!  (6 ) NBIOTRCE (INTEGER) : Number of bioburn species  [unitless           ]
!  (7 ) JLOP     (INTEGER) : SMVGEAR grid box index     [unitless           ]
!  (8 ) REMIS    (REAL*8 ) : SMVGEAR emissions array    [molec species/cm3/s]
!  (9 ) VOLUME   (REAL*8 ) : SMVGEAR volume array       [cm3                ]
!
!  NOTES: 
!  (1 ) Original code from Harvard Tropospheric Chemistry Module for 3-D 
!        applications by Larry Horowitz, Jinyou Liang, Gerry Gardner, 
!        Prof. Daniel Jacob of Harvard University (Release V2.0)  
!  (2 ) New version 3.0 by Bob Yantosca to place NOx emissions into boxes  
!        above the surface. (bmy, 6/8/98)     
!  (3 ) Also now do chemistry up to the location of the annual mean         
!         tropopause (bmy, 12/9/99)                                         
!  (4 ) BURNEMIS is now dynamically allocatable and is contained in F90       
!        module "biomass_mod.f".  BIOTRCE and NBIOTRCE are also contained
!        in "biomass_mod.f".  (bmy, 9/12/00)                     
!  (5 ) BIOFUEL is now dynamically allocatable and is contained in F90
!        module "biofuel_mod.f".  BFTRACE and NBFTRACE are also contained
!        in "biofuel_mod.f" (bmy, 9/12/00, 4/17/01)
!  (6 ) BURNEMIS and BIOFUEL are now treated as true global arrays,  
!        and need to be referenced by the global offset variables          
!        IREF = I + I0 and JREF = J + J0 (bmy, 9/12/00)                    
!  (7 ) Now reference JLOP, REMIS, VOLUME from F90 module "comode_mod.f",  
!        in order to save memory (bmy, 10/19/00)                           
!  (8 ) Now add in up to NBFTRACE biofuel species (bmy, 4/17/01) 
!  (9 ) Add new subroutine header, updated comments, cosmetic changes.
!        (bmy, 4/17/01)  
!  (10) Updated comments -- GEMISNOX is [molec/cm3/s]. (bdf, bmy, 6/7/01)     
!  (11) For GEOS-3, we now distributes surface emissions throughout the 
!        boundary layer.  This is necessary since the first couple of GEOS-3 
!        surface layers are very thin.  Piling up of emissions into a small 
!        layer will cause SMVGEAR to choke.  (bdf, bmy, 6/15/01)
!  (12) Also now reference BFTRACE and NBFTRACE from "biofuel_mod.f", 
!        and reference AD12 from "diag_mod.f". (bdf, bmy, 6/15/01)
!  (13) For GEOS-1, GEOS-STRAT, emit into the surface layer, as we did
!        in prior versions. (bmy, 6/26/01)
!  (14) Bug fix: corrected a typo for the biofuel emissions (bmy, 7/10/01)
!  (15) Bug fix: make sure BIOMASS and BIOFUEL, and SOIL NOx emissions have 
!        units of [molec/box/s] before distributing thru the boundary layer.  
!        This involves multiplication by VOLUME(JLOOP1) and division by
!        VOLUME(JLOOP). (bmy, 7/16/01)
!  (16) XTRA2(IREF,JREF,5) is now XTRA2(I,J).  BIOFUEL(:,IREF,JREF) is now
!        BIOFUEL(:,I,J). BURNEMIS(:,IREF,JREF) is now BURNEMIS(:,I,J).
!        Replace PW(I,J) with P(I,J). (bmy, 9/28/01)
!  (17) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (18) Now references GET_PEDGE from "pressure_mod.f", to compute P at 
!        the bottom edge of grid box (I,J,L).  (dsa, bdf, bmy, 8/21/02)
!  (19) Now reference IDTNOX, IDENOX, etc from "tracerid_mod.f" (bmy, 11/6/02)
!  (20) Remove references to IREF, JREF (bmy, 2/11/03)
!  (21) NEMIS is now NEMIS(NCS) for SMVGEAR II (gcc, bdf, bmy, 4/1/03)
!  (22) Added parallel loop over N.  Also directly substituted JLOP(I,J,1) 
!        for all instances of JLOOP1.  Updated comments. (hamid, bmy, 3/19/04)
!  (23) Bug fix for COMPAQ compiler...do not use EXIT from w/in parallel loop.
!        (auvray, bmy, 11/29/04)
!  (24) Now replace XTRA2 with GET_PBL_TOP_L in "pbl_mix_mod.f".  Now remove
!        reference to CMN, it's obsolete.  Now references GET_TPAUSE_LEVEL
!        from "tropopause_mod.f" (bmy, 8/22/05)
!  (25) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (26) Now updated for new "biomass_mod.f" (bmy, 4/5/06)
!  (27) Now account for the different definition of tropopause in case 
!        of variable tropopause.   The BIOMASS array from "biomass_mod.f" is 
!        now in units of [molec CO/cm2/s].  Adjust unit conversion accordingly.
!        Also replace NBIOMAX with NBIOMAX_GAS, since aerosol biomass is
!        handled elsewhere.  (bdf, phs, bmy, 9/27/06)
!  (28) Now replace GEMISNOX array (from CMN_NOX) with module arrays
!        EMIS_LI_NOx and EMIS_AC_NOx (ltm, bmy, 10/3/07)
!  (29) Bug fix: resize EMISRR to be consistent w/ CMN_O3 (bmy, jaf, 6/11/08)
!  (30) Limit emissions into the surface level only (lin, 5/29/09)
!  (31) Bug fix: cycle if IDEMIS(NN) <= 0 to avoid array-out-of-bounds
!        errors (bmy, 8/6/09)
!  (32) Check for emissions above PBL -anthro NOx only for now- (phs, 10/27/09)
!******************************************************************************
!
      ! References to F90 modules 
      USE AIRCRAFT_NOX_MOD,  ONLY : EMIS_AC_NOx
      USE BIOFUEL_MOD,       ONLY : BIOFUEL,   BFTRACE, NBFTRACE
      USE BIOMASS_MOD,       ONLY : BIOMASS,   BIOTRCE, NBIOMAX_GAS
      USE COMODE_MOD,        ONLY : JLOP,      REMIS,   VOLUME
      USE COMODE_MOD,        ONLY : IYSAVE
      USE DIAG_MOD,          ONLY : AD12
      USE GRID_MOD,          ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,       ONLY : LVARTROP
      USE LIGHTNING_NOX_MOD, ONLY : EMIS_LI_NOx
      USE PBL_MIX_MOD,       ONLY : GET_PBL_TOP_L
      USE PRESSURE_MOD,      ONLY : GET_PEDGE
      USE TRACERID_MOD,      ONLY : CTRMB,     IDEMIS,  IDENOX
      USE TROPOPAUSE_MOD,    ONLY : GET_TPAUSE_LEVEL
      USE LOGICAL_MOD,       ONLY : LNLPBL ! (Lin, 03/31/09)

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_NOX"   ! GEMISNOX2
#     include "CMN_DIAG"  ! Diagnostic flags
#     include "comode.h"  ! IDEMS, NEMIS

      ! Arguments
      REAL*8,  INTENT(IN) :: EMISRR(IIPAR,JJPAR,NEMPARA+NEMPARB)
      REAL*8,  INTENT(IN) :: EMISRRN(IIPAR,JJPAR,NOXEXTENT)  

      ! Local variables
      LOGICAL             :: IS_LI_NOx, IS_AC_NOx
      INTEGER             :: I, J,  JLOOP, JLOOP1, LTROP
      INTEGER             :: L, LL, N, NN,  NBB, NBF, TOP, TOPMIX
      REAL*8              :: COEF1,   TOTPRES, DELTPRES
      REAL*8              :: EMIS_BL, NOXTOT,  TOTAL, A_CM2

      !=================================================================
      ! SETEMIS begins here!
      !=================================================================

      ! some ajdustments for non-local PBL (Lin, 03/31/09)
      call flush(6)
      IF (NCS == 0) THEN
        REMIS(:,:)=0.
        RETURN
      ENDIF

      ! Test if the EMIS_LI_NOx and EMIS_AC_NOx arrays are allocated
      IS_LI_NOx = ALLOCATED( EMIS_LI_NOx )
      IS_AC_NOX = ALLOCATED( EMIS_AC_NOX )

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( N,     NN,  NBB,     NBF,    I,        J,     L, JLOOP )
!$OMP+PRIVATE( COEF1, TOP, TOTPRES, NOXTOT, DELTPRES, EMIS_BL,  A_CM2 )
!$OMP+PRIVATE( TOPMIX )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over emission species
      DO N = 1, NEMIS(NCS)

         ! Get CTM tracer number NN corresponding to emission species N
         NN = IDEMS(N)
         IF ( NN == 0 ) CYCLE

         ! We have to search for the biomass burning species in 
         ! BIOTRCE with the same CTM tracer number NN as in IDEMS
         NBB = 0
         IF ( ALLOCATED( BIOMASS ) ) THEN 
            DO I = 1, NBIOMAX_GAS
               IF ( BIOTRCE(I) == NN ) THEN 
                 NBB = I
                  EXIT
               ENDIF
            ENDDO
         ENDIF

         ! We have to search for the biofuel burning species in 
         ! BFTRACE with the same CTM tracer number NN as in IDEMS
         NBF = 0
         IF ( ALLOCATED( BIOFUEL ) ) THEN
            DO I = 1, NBFTRACE
               IF ( BFTRACE(I) == NN ) THEN
                  NBF = I
                  EXIT
              ENDIF
            ENDDO
         ENDIF            

         ! Initialize REMIS(*,N) -- the emission rate array
         DO JLOOP = 1, NTTLOOP
            REMIS(JLOOP,N) = 0d0
         ENDDO       

         ! Skip to next tracer if IDEMIS(NN) is not defined in 
         ! order to avoid array-out-of-bounds errors (bmy, 8/6/04)
         IF ( IDEMIS(NN) <= 0 ) THEN
            PRINT*, 'Tracer ', NN, ' is not an emitted species!'
            CYCLE
         ENDIF

         ! COEF1 = molecules of emission species / molecules of tracer
         COEF1 = 1.0 + CTRMB(NN, IDEMIS(NN))         

         ! Loop over Lat and Long boxes
         DO J = 1, NLAT
         DO I = 1, NLONG

            !===========================================================
            ! For GEOS-3: distribute surface emissions throughout the
            ! entire boundary layer.  Define some variables here.
            ! (bdf, 6/15/01)
            !===========================================================

            ! Top level of the boundary layer
            ! guard for b.l. being in first level.
            TOP = FLOOR( GET_PBL_TOP_L( I, J ) )
            IF ( TOP == 0 ) TOP = 1

            ! Check for possible emissions above PBL (phs, 27/10/09)
            TOPMIX = MIN( TOP, NOXEXTENT )
            
            ! Add option for non-local PBL (Lin, 03/31/09)
            IF (LNLPBL) TOP = 1

            ! Pressure thickness of entire boundary layer [hPa]
            TOTPRES = GET_PEDGE(I,J,1) - GET_PEDGE(I,J,TOP+1)

            ! For NOx only....
            IF ( N == IDENOX ) THEN

               !========================================================
               ! Anthropogenic NOx emissions [molec/box/s]
               ! Distribute emissions thru the entire boundary layer
               !========================================================

               ! Sum anthro NOx emissions over all levels [molec NOx/box/s]
               NOXTOT = 0d0
!--prior to 27/10/09 (phs)
!               DO L = 1, NOXEXTENT
               DO L = 1, TOPMIX
                  NOXTOT = NOXTOT + EMISRRN(I,J,L)
               ENDDO

               ! Loop over the boundary layer
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0

                  IF ( JLOOP /= 0 ) THEN

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

                     ! Add option for non-local PBL (Lin, 03/31/09)
                     IF (LNLPBL) DELTPRES = TOTPRES

                     ! Of the total anthro NOx, the fraction DELTPRES/TOTPRES
                     ! goes into level L, since that is the fraction of the
                     ! boundary layer occupied by level L.  Also divide NOx 
                     ! by COEF1 to convert from [molec NOx/box/s] to 
                     ! [molec NO/box/s], which is really what gets emitted.
                     EMIS_BL        = ( NOXTOT   / COEF1   ) *
     &                                ( DELTPRES / TOTPRES )

                     ! Convert anthro NOx emissions from [molec NO/box/s]
                     ! to [molec NO/cm3/s] and store in the REMIS array
                     REMIS(JLOOP,N) = EMIS_BL / VOLUME(JLOOP)
                  ENDIF
               ENDDO

               ! put any emissions above PBL into corresponding
               ! REMIS (phs, 27/10/09)
               IF ( NOXEXTENT > TOPMIX ) THEN
                  DO L = TOPMIX+1, NOXEXTENT
                     JLOOP          = JLOP(I,J,L)
                     REMIS(JLOOP,N) = EMISRRN(I,J,L) /
     $                    ( VOLUME(JLOOP) * COEF1 )
                  ENDDO
               ENDIF   
                     
               
               !========================================================
               ! Soil Nox emissions [molec/cm3/s] 
               ! Distribute emissions thru the entire boundary layer
               !========================================================
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0

                  IF ( JLOOP /= 0 ) THEN

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)
                     
                     ! Add option for non-local PBL (Lin, 03/31/09)
                     IF (LNLPBL) DELTPRES = TOTPRES

                     ! Soil NOx is in [molec/cm3/s], so we need to multiply
                     ! by VOLUME(JLOP(I,J,1)) to convert it to [molec/box/s],
                     ! VOLUME(JLOP(I,J,1)) is the volume in cm3 of the surface
                     ! grid box (I,J,1).  Then we need to divide that by 
                     ! COEF1 to convert from [molec NOx/box/s] to 
                     ! [molec NO/box/s], which is really what gets emitted.  
                     ! Of the total soil NOx, the fraction DELTPRES/TOTPRES 
                     ! goes into level L, since that is the fraction of the 
                     ! boundary layer occupied by level L.  Store in EMIS_BL.
                     EMIS_BL    = ( GEMISNOX2(I,J)                    *
     &                              VOLUME( JLOP(I,J,1) ) / COEF1   ) * 
     &                              ( DELTPRES            / TOTPRES )

                     ! Since EMIS_BL is in [molec NO/box/s], we have to
                     ! divide by VOLUME(JLOOP), which is the volume of the
                     ! grid box (I,J,L) to convert back to [molec/cm3/s].
                     ! Store in the REMIS array for SMVGEAR.
                     REMIS(JLOOP,N) = REMIS(JLOOP,N) + 
     &                                ( EMIS_BL / VOLUME(JLOOP) ) 
                  ENDIF
               ENDDO

               !========================================================
               ! Aircraft and Lightning NOx [molec/cm3/s]
               ! Distribute emissions in the troposphere
               !========================================================

               ! bdf - variable tropopause is a tropospheric box
               IF ( LVARTROP ) THEN 
                  LTROP = GET_TPAUSE_LEVEL( I, J ) 
               ELSE
                  LTROP = GET_TPAUSE_LEVEL( I, J ) - 1
               ENDIF


               DO L = 1, LTROP 
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0

                  IF ( JLOOP /= 0 ) THEN

                     !-----------------
                     ! Lightning NOx
                     !-----------------
                     IF ( IS_LI_NOx ) THEN
 
                        ! Divide lightning NOx by COEF1 to convert
                        ! from [molec NOx/cm3/s] to [molec NO/cm3/s], since
                        ! NO is the actual emission species for NOx.
                        EMIS_BL        = EMIS_LI_NOx(I,J,L) / COEF1 

                        ! Save lightning NOx [molec NO/cm3/s] in REMIS
                        REMIS(JLOOP,N) = REMIS(JLOOP,N) + EMIS_BL
                     ENDIF

                     !-----------------
                     ! Aircraft NOx
                     !-----------------
                     IF ( IS_AC_NOx ) THEN

                        ! Divide aircraft NOx by COEF1 to convert
                        ! from [molec NOx/cm3/s] to [molec NO/cm3/s], since
                        ! NO is the actual emission species for NOx.
                        EMIS_BL        = EMIS_AC_NOx(I,J,L) / COEF1 

                        ! Save aircraft NOx [molec NO/cm3/s] in REMIS
                        REMIS(JLOOP,N) = REMIS(JLOOP,N) + EMIS_BL
                     ENDIF

                  ENDIF
               ENDDO

            ELSE

               !========================================================
               ! Anthropogenic tracers other than NOx [molec/box/s]
               ! Distribute emissions thru the entire boundary layer
               !========================================================
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0

                  IF ( JLOOP /= 0 ) THEN 

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

                     ! Add option for non-local PBL (Lin, 03/31/09)
                     IF (LNLPBL) DELTPRES = TOTPRES

                     ! Of the total tracer, the fraction DELTPRES/TOTPRES
                     ! goes into level L, since that is the fraction of the
                     ! boundary layer occupied by level L.  Also divide the
                     ! tracer by COEF1 to convert from [molec tracer/box/s] 
                     ! to [molec species/box/s].  For example, ISOPRENE is
                     ! carried by GEOS-CHEM as 5 carbons, so you would divide
                     ! by 5 to get "effective molecules" of ISOPRENE.
                     EMIS_BL        = ( EMISRR(I,J,N) / COEF1   ) *
     &                                ( DELTPRES      / TOTPRES )

                     ! Convert emissions from [molec species/box/s] to 
                     ! [molec species/cm3/s] and store in the REMIS array
                     REMIS(JLOOP,N) = EMIS_BL / VOLUME(JLOOP)
                  ENDIF
               ENDDO
            ENDIF

            !===========================================================
            ! Add biomass burning source [molec/cm3/s]
            ! Distribute emissions thru the entire boundary layer
            !===========================================================
            IF ( NBB /= 0 ) THEN
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0

                  IF ( JLOOP /= 0 ) THEN

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

                     ! Add option for non-local PBL (Lin, 03/31/09)
                     IF (LNLPBL) DELTPRES = TOTPRES
 
                     ! Grid box area [cm2]
                     A_CM2    = GET_AREA_CM2( IYSAVE(JLOOP) )
               
                     ! Biomass burning is in [molec/cm2/s], so we need to 
                     ! multiply by A_CM2 to convert it to [molec/box/s].
                     ! Then we need to divide that by COEF1 to convert from 
                     ! [molec tracer/box/s] to [molec species/box/s].
                     ! Of the total biomass burning emissions, the fraction 
                     ! DELTPRES/TOTPRES goes into level L, since that is the 
                     ! fraction of the boundary layer occupied by level L.  
                     ! Store in EMIS_BL.
                     EMIS_BL  = ( BIOMASS(I,J,NBB) * A_CM2 / COEF1   ) *
     &                          ( DELTPRES                 / TOTPRES )

                     ! Since EMIS_BL is in [molec species/box/s], we 
                     ! have to divide by VOLUME(JLOOP), which is the 
                     ! volume of the grid box (I,J,L) to convert back to 
                     ! [molec species/cm3/s].  Store in the REMIS array.
                     REMIS(JLOOP,N) = REMIS(JLOOP,N) + 
     &                                ( EMIS_BL / VOLUME(JLOOP) )
                  ENDIF
               ENDDO
            ENDIF

            !===========================================================
            ! Add biofuel burning source [molec/cm3/s]
            ! Distribute emissions thru the entire boundary layer 
            !===========================================================
            IF ( NBF /= 0 ) THEN
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0

                  IF ( JLOOP /= 0 ) THEN

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

                     ! Add option for non-local PBL (Lin, 03/31/09)
                     IF (LNLPBL) DELTPRES = TOTPRES

                     ! Biofuel burning is in [molec/cm3/s], so we need to 
                     ! multiply by VOLUME(JLOP(I,J,1)) to convert it to 
                     ! [molec/box/s], VOLUME(JLOP(I,J,1)) is the volume in cm3 
                     ! of the surface grid box (I,J,1).  Then we need to 
                     ! divide that by COEF1 to convert from 
                     ! [molec tracer/box/s] to [molec species/box/s].
                     ! Of the total biofuel burning emissions, the fraction 
                     ! DELTPRES/TOTPRES goes into level L, since that is the 
                     ! fraction of the boundary layer occupied by level L.  
                     ! Store in EMIS_BL.
                     EMIS_BL  = ( BIOFUEL(NBF,I,J) *
     &                            VOLUME( JLOP(I,J,1) ) / COEF1 ) *
     &                            ( DELTPRES / TOTPRES )

                     ! Since EMIS_BL is in [molec species/box/s], we 
                     ! have to divide by VOLUME(JLOOP), which is the 
                     ! volume of the grid box (I,J,L) to convert back to 
                     ! [molec species/cm3/s].  Store in the REMIS array.
                     REMIS(JLOOP,N) = REMIS(JLOOP,N) + 
     &                                ( EMIS_BL / VOLUME(JLOOP) )
                  ENDIF
               ENDDO
            ENDIF

            !===========================================================
            ! ND12 Diagnostic: Save the fraction of the boundary layer
            ! occupied by level L into the AD12 diagnostic array.
            !===========================================================
            IF ( N == 1 .and. ND12 > 0 ) THEN
               DO L = 1, MIN( FLOOR( GET_PBL_TOP_L( I, J ) ), LD12 )

                  ! Thickness of layer L [mb]
                  DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)
 
                 ! Save boundary layer fraction into AD12
                  AD12(I,J,L) = AD12(I,J,L) + ( DELTPRES / TOTPRES )
               ENDDO
            ENDIF

         ENDDO  ! I
         ENDDO  ! J
      ENDDO     ! N
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SETEMIS
