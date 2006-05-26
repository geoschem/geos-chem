! $Id: setemis.f,v 1.9 2006/05/26 17:45:26 bmy Exp $
      SUBROUTINE SETEMIS( EMISRR, EMISRRN )
!
!******************************************************************************
!  Subroutine SETEMIS places emissions computed from GEOS-CHEM 
!  subroutines into arrays for SMVGEAR II chemistry. 
!  (lwh, jyl, gmg, djj, bdf, bmy, 6/8/98, 4/5/06)
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
!******************************************************************************
!
      ! References to F90 modules 
      USE BIOFUEL_MOD,    ONLY : BIOFUEL,     BFTRACE, NBFTRACE
      USE BIOMASS_MOD,    ONLY : BIOMASS,     BIOTRCE, NBIOMAX
      USE COMODE_MOD,     ONLY : JLOP,        REMIS,   VOLUME
      USE DIAG_MOD,       ONLY : AD12
      USE PBL_MIX_MOD,    ONLY : GET_PBL_TOP_L
      USE PRESSURE_MOD,   ONLY : GET_PEDGE
      USE TRACERID_MOD,   ONLY : CTRMB,       IDEMIS,  IDENOX
      USE TROPOPAUSE_MOD, ONLY : GET_TPAUSE_LEVEL

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_NOX"   ! GEMISNOX, GEMISNOX2
#     include "CMN_DIAG"  ! Diagnostic flags
#     include "comode.h"  ! IDEMS, NEMIS

      ! Arguments
      REAL*8,  INTENT(IN) :: EMISRR(IIPAR,JJPAR,2:NEMPARA+NEMPARB)
      REAL*8,  INTENT(IN) :: EMISRRN(IIPAR,JJPAR,NOXEXTENT)  

      ! Local variables
      INTEGER             :: I, J,  JLOOP, JLOOP1, LTROP
      INTEGER             :: L, LL, N, NN,  NBB, NBF, TOP
      REAL*8              :: COEF1,   TOTPRES, DELTPRES
      REAL*8              :: EMIS_BL, NOXTOT,  TOTAL

      !=================================================================
      ! SETEMIS begins here!
      !
      ! Loop over emission species
      !=================================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( N,     NN,  NBB,     NBF,    I,        J,     L, JLOOP )
!$OMP+PRIVATE( COEF1, TOP, TOTPRES, NOXTOT, DELTPRES, EMIS_BL         )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, NEMIS(NCS)

         ! Get CTM tracer number NN corresponding to emission species N
         NN = IDEMS(N)
         IF ( NN == 0 ) CYCLE

         ! We have to search for the biomass burning species in 
         ! BIOTRCE with the same CTM tracer number NN as in IDEMS
         NBB = 0
         IF ( ALLOCATED( BIOMASS ) ) THEN 
            DO I = 1, NBIOMAX
               IF ( BIOTRCE(I) == NN ) THEN 
                  NBB = I
#if   defined( COMPAQ )
                  ! COMPAQ has an issue with EXIT from w/in parallel loop
                  ! (auvray, bmy, 11/29/04)
#else
                  EXIT
#endif
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
#if   defined( COMPAQ )
                  ! COMPAQ has an issue with EXIT from w/in parallel loop
                  ! (auvray, bmy, 11/29/04)
#else
                  EXIT
#endif 
              ENDIF
            ENDDO
         ENDIF            

         ! Initialize REMIS(*,N) -- the emission rate array
         DO JLOOP = 1, NTTLOOP
            REMIS(JLOOP,N) = 0d0
         ENDDO       

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
               DO L = 1, NOXEXTENT
                  NOXTOT = NOXTOT + EMISRRN(I,J,L)
               ENDDO

               ! Loop over the boundary layer
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0

                  IF ( JLOOP /= 0 ) THEN

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

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
               LTROP = GET_TPAUSE_LEVEL( I, J ) - 1
               DO L = 1, LTROP 
                  JLOOP   = JLOP(I,J,L)
                  EMIS_BL = 0d0

                  IF ( JLOOP /= 0 ) THEN

                     ! Divide aircraft & lightning NOx by COEF1 to convert
                     ! from [molec NOx/cm3/s] to [molec NO/cm3/s], since
                     ! NO is the actual emission species for NOx.
                     EMIS_BL        = GEMISNOX(I,J,L) / COEF1

                     ! Save aircraft/lightning NOx [molec NO/cm3/s] in REMIS
                     REMIS(JLOOP,N) = REMIS(JLOOP,N) + EMIS_BL
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
                                          
                     ! Biomass burning is in [molec/cm3/s], so we need to 
                     ! multiply by VOLUME(JLOP(I,J,1)) to convert it to 
                     ! [molec/box/s], VOLUME(JLOP(I,J,1)) is the volume in cm3 
                     ! of the surface grid box (I,J,1).  Then we need to 
                     ! divide that by COEF1 to convert from 
                     ! [molec tracer/box/s] to [molec species/cm3/s].
                     ! Of the total biomass burning emissions, the fraction 
                     ! DELTPRES/TOTPRES goes into level L, since that is the 
                     ! fraction of the boundary layer occupied by level L.  
                     ! Store in EMIS_BL.
                     EMIS_BL  = ( BIOMASS(I,J,NBB)  *
     &                            VOLUME( JLOP(I,J,1) ) / COEF1  ) *
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

                     ! Biofuel burning is in [molec/cm3/s], so we need to 
                     ! multiply by VOLUME(JLOP(I,J,1)) to convert it to 
                     ! [molec/box/s], VOLUME(JLOP(I,J,1)) is the volume in cm3 
                     ! of the surface grid box (I,J,1).  Then we need to 
                     ! divide that by COEF1 to convert from 
                     ! [molec tracer/box/s] to [molec species/cm3/s].
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
               DO L = 1, MIN( TOP, LD12 )

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
