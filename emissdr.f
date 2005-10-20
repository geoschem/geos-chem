! $Id: emissdr.f,v 1.9 2005/10/20 14:03:24 bmy Exp $
      SUBROUTINE EMISSDR
!
!******************************************************************************
!  Subroutine EMISSDR computes emissions for the full chemistry simulation
!  (NSRCX == 3).  Emissions are stored in GEMISNOX and EMISRR arrays, 
!  which are then passed to the SMVGEAR subroutines (bmy, 10/8/98, 10/3/05)
!
!  NOTES:
!  (1 ) Now accounts for seasonal NOx emissions, and multi-level NOx 
!        emissions. (bmy, 10/8/98)
!  (2 ) Surface NOx and 100m NOx are now placed into the correct sigma level.
!        (bmy, 10/8/98)
!  (3 ) Eliminate GISS-Specific code (bmy, 3/15/99)
!  (4 ) Now includes monoterpenes for ACETONE emissions (bdf, 4/8/99)
!  (5 ) Now uses allocatable arrays for ND29, ND36, and ND46 diagnostics.
!        Also made some cosmetic changes, and updated comments (bmy, 3/16/00)
!  (6 ) Eliminate obsolete code and ND63 diagnostic (bmy, 4/12/00)
!  (7 ) Add reference to BIOBURN in "biomass_mod.f" and to BIOFUEL_BURN
!        in "biofuel_mod.f".  Also remove references to BURNEMIS and TWOODIJ.
!        (bmy, 9/12/00)
!  (8 ) Remove reference to "biomass.h" -- that is replaced by F90 module
!        "biomass_mod.f". (bmy, 9/25/00)
!  (9 ) Remove obsolete code from 9/12/00 and 9/25/00 (bmy, 12/21/00)
!  (10) Add CO source from monoterpenes (bnd, bmy, 12/21/00)
!  (11) Add monoterpene source to the ND46 diagnostic.  Renamed EMXX to
!        EMIS (for Isoprene), and commented out Larry Horowitz's "special
!        cases".  Also made some cosmetic changes. (bmy, 1/2/01)
!  (12) Added CO source from CH3OH oxidation (bmy, 1/3/01)
!  (13) Removed obsolete code from 1/2/01 (bmy, 3/15/01)
!  (14) Now initialize GEMISNOX2.  Also updated comments. (bdf, bmy, 6/15/01)
!  (15) Now references routines from "acetone_mod.f" for the biogenic
!        emission of acetone into the SMVGEAR arrays.  Now use 
!        EMISRR(I,J,IDEACET) to archive ND46, since the biogenic
!        acetone emissions are now computed in this array.  Also define
!        XNUMOL_C so as not to rely on IDTISOP being defined.  Also add
!        LASTMONTH variable to flag when we change month. (bmy, 9/4/01)
!  (16) Now reference AIREMISS from "aircraft_nox_mod.f" (bmy, 2/14/02)
!  (17) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE.  Also removed
!        obsolete, commented-out code. (bmy, 6/25/02)
!  (18) Now references IDTNOX, etc. from "tracerid_mod.f".  Now references
!        SUNCOS from "dao_mod.f".  Now make FIRSTEMISS a local SAVEd variable
!        instead of an argument. (bmy, 11/15/02)
!  (19) Now replaced DXYP(JREF)*1d4 with GET_AREA_CM2 from "grid_mod.f".
!        Now remove MONTH from call to BIOBURN.  Now use functions GET_MONTH,
!        GET_LOCALTIME, GET_ELAPSED_MIN, GET_TS_EMIS, GET_LOCALTIME from 
!        "time_mod.f".  Now use functions GET_XOFFSET and GET_YOFFSET from 
!        "grid_mod.f". (bmy, 2/11/03)
!  (20) Now pass I, J to EMISOP, EMISOP_GRASS, EMISOP_MB (bmy, 12/9/03)
!  (21) Now references EMLIGHTNING from "lightning_nox_mod.f" (bmy, 4/14/04)
!  (22) Now references "logical_mod.f".  Now replaced LFOSSIL with LANTHRO.
!        (bmy, 7/20/04)
!  (23) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ACETONE_MOD,       ONLY : EMISS_BIOACET, OCEAN_SOURCE_ACET
      USE ACETONE_MOD,       ONLY : READ_JO1D,     READ_RESP
      USE AIRCRAFT_NOX_MOD,  ONLY : AIREMISS
      USE BIOFUEL_MOD,       ONLY : BIOFUEL_BURN
      USE BIOMASS_MOD,       ONLY : BIOBURN
      USE DAO_MOD,           ONLY : SUNCOS
      USE DIAG_MOD,          ONLY : AD29,          AD46
      USE GRID_MOD,          ONLY : GET_AREA_CM2
      USE GRID_MOD,          ONLY : GET_XOFFSET,   GET_YOFFSET
      USE LIGHTNING_NOX_MOD, ONLY : EMLIGHTNING
      USE LOGICAL_MOD,       ONLY : LANTHRO,       LLIGHTNOX, LSOILNOX  
      USE LOGICAL_MOD,       ONLY : LAIRNOX,       LBIONOX,   LWOODCO   
      USE TIME_MOD,          ONLY : GET_MONTH,     GET_TAU
      USE TIME_MOD,          ONLY : GET_TS_EMIS,   GET_LOCALTIME
      USE TRACERID_MOD,      ONLY : IDEACET,       IDTISOP,   IDEISOP   
      USE TRACERID_MOD,      ONLY : IDECO,         IDEPRPE,   NEMANTHRO 

      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! STT, many other variables
#     include "CMN_DIAG"     ! Diagnostic arrays and switches
#     include "CMN_O3"       ! Emissions arrays, XNUMOL
#     include "CMN_NOX"      ! GEMISNOX, GEMISNOX2
#     include "CMN_MONOT"    ! Monoterpenes
#     include "comode.h"     ! IVERT?

      ! Local variables
      LOGICAL, SAVE          :: FIRSTEMISS = .TRUE. 
      INTEGER                :: I, J, L, N, IJLOOP 
      INTEGER                :: I0, J0, IOFF, JOFF, IREF, JREF
      INTEGER                :: NDAY,JSCEN,NN
      INTEGER, SAVE          :: LASTMONTH = -99
      REAL*8                 :: DTSRCE,  XLOCTM,   EMIS,  TMMP
      REAL*8                 :: BIOSCAL, AREA_CM2, EMMO,  ACETSCAL
      REAL*8                 :: TMPVAL,  EMMB,     GRASS, BIO_ACET
      REAL*8                 :: CONVERT(NVEGTYPE), GMONOT(NVEGTYPE)
      
      ! Molecules C / kg C
      REAL*8,  PARAMETER     :: XNUMOL_C = 6.022d+23 / 12d-3 

      ! External functions
      REAL*8, EXTERNAL       :: BOXVL,   XLTMMP,    EMISOP 
      REAL*8, EXTERNAL       :: EMMONOT, EMISOP_MB, EMISOP_GRASS
!
!******************************************************************************
!  EMISSDR begins here!
!
!  DTSRCE = emission timestep in seconds
!
!  Call subroutines to set up ISOP and monoterpene emission (first time only!)
!******************************************************************************
!
      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Get nested-grid offsets
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

      IF ( FIRSTEMISS ) THEN
          CALL RDLIGHT
          CALL RDISOPT( CONVERT )
          CALL RDMONOT( GMONOT  )
          CALL SETBASE( CONVERT, GMONOT )
          FIRSTEMISS = .FALSE.
      ENDIF
!
!******************************************************************************
!  BE CAREFUL WITH USING A WINDOW RELATIVE TO THE EMISSIONS WINDOW
!  NEED SPECIFY THE OFFSET OF THE SUB-WINDOW
!
!  first zero the arrays in which emissions will be stored
!
!  EMISRRN(I,J,L)   = Emission rate of NOx (N = IDENOX ) into box (I,J,L).
!                     Units are [molec NOx/box/s].
!
!  EMISRR(I,J,N)    = Emission rate of tracer N into surface box (I,J,1). 
!                     Units are [molec tracer/box/s].
!
!  GEMISNOX(I,J,L)  = Array which stores NOx emissions from aircraft,
!                     and lightning.  Units are [molec NOx/cm3/s].
!
!  GEMISNOX2(I,J)   = Array which stores NOx emissions from soils.
!                     Units are [molec NOx/cm3/s].
!
!  NOTE: Now use F90 array initialization syntax (bmy, 3/15/99)
!******************************************************************************
!
      ! These need to be initialized on every call
      EMISRRN   = 0d0
      EMISRR    = 0d0
      GEMISNOX  = 0d0
      GEMISNOX2 = 0d0
      
      ! Loop over latitudes
      IJLOOP = 0
      DO J = 1, JJPAR
         JREF = J + J0

         ! Compute surface area of grid boxes in cm^2 
         AREA_CM2 = GET_AREA_CM2( J )

         ! Loop over longitues
         DO I = 1, IIPAR
            IREF   = I + I0
            IJLOOP = IJLOOP + 1
         
            ! Zero biogenic acetone (bmy, 9/14/01)
            BIO_ACET  = 0d0

            ! Use function GET_LOCALTIME to get the local time at lon I
            ! Middle of time step is between 10pm-2am when IHOUR = 1
            IHOUR = INT( ( GET_LOCALTIME( I ) ) / 4 ) + 1
            IF ( IHOUR .EQ. 7 ) IHOUR = 1

            !=================================================================
            ! attenuate emissions on the weekend ---
            ! scale factors for Saturday/Sunday/Weekday must average out to 1!
            !    JSCEN = 1 Saturday 
            !    JSCEN = 2 Sunday
            !    JSCEN = 3 Weekday 
            !
            ! 1 Jan 1980 and 1 Jan 1985 were both Tuesdays, so NDAY mod 7 = 4 
            ! is a Saturday and NDAY mod 7 = 5 is a Sunday (bmy, 3/23/98)
            !=================================================================
            NDAY = ( GET_TAU() / 24d0 ) 
            IF ( MOD( NDAY, 7 ) .eq. 4 ) THEN
               JSCEN = 1
            ELSE IF ( MOD( NDAY, 7 ) .eq. 5 ) THEN
               JSCEN = 2
            ELSE
               JSCEN = 3
            ENDIF

            ! Fossil Fuel emissions (kg / Grid-Box / Time-Step)
            ! NN = tracer number corresponding to emission species N
            IF ( LANTHRO ) THEN 
               DO N = 1, NEMANTHRO
                  NN = IDEMS(N)
                  IF ( NN /= 0 ) THEN
                     CALL EMFOSSIL( I, J, N, NN, IREF, JREF, JSCEN )
                  ENDIF
               ENDDO
            ENDIF

!-----------------------------------------------------------------------------
! LIGHTNING EMISSIONS NOX [molecules/cm3/s]
!
            IF ( LLIGHTNOX ) CALL EMLIGHTNING( I, J )
!-----------------------------------------------------------------------------
! SOIL EMISSIONS NOX [molecules/cm3/s]
! Now have to pass SUNCOS to SOILNOXEMS and SOILCRF (bmy, 10/20/99)
!
            IF ( LSOILNOX .AND. I == 1 .AND. J == 1 ) 
     &         CALL SOILNOXEMS( SUNCOS )
!-----------------------------------------------------------------------------
! AIRCRAFT emissions NOx [molecules/cm3/s]
!
            IF ( LAIRNOX .AND. I == 1 .AND. J == 1 ) CALL AIREMISS
!-----------------------------------------------------------------------------
! BIOMASS BURNING emissions NOx [molecules/cm3/s]
!
            IF ( LBIONOX .AND. I == 1 .AND. J == 1 ) CALL BIOBURN
!-----------------------------------------------------------------------------
! NOx AND CO from biofuel combustion [kg/box]  
!
            IF ( LWOODCO .AND. I == 1 .AND. J == 1 ) CALL BIOFUEL_BURN
!----------------------------------------------------------------------------
! BIOGENIC EMISSIONS OF VARIOUS QUANTITIES [Atoms C/box/time step]
!
            ! Temperature
            TMMP  = XLTMMP(I,J,IJLOOP)

            ! Monoterpenes 
            EMMO  = EMMONOT( IJLOOP, TMMP, XNUMOL_C )

            ! Isoprene
            EMIS  = EMISOP( I, J, IJLOOP, SUNCOS, TMMP, XNUMOL_C )
 
            ! Methyl Butenol (MBO)
            EMMB  = EMISOP_MB( I, J, IJLOOP, SUNCOS, TMMP, XNUMOL_C )

            ! Isoprene emissions from grasslands
            GRASS = EMISOP_GRASS( I, J, IJLOOP, SUNCOS, TMMP, XNUMOL_C) 
!-----------------------------------------------------------------------------
! BIOGENIC ACETONE EMISSIONS
!
            IF ( IDEACET /= 0 ) THEN

               ! Read monthly mean JO1D and leaf respiration values
               ! These will be stored internally in "acetone_mod.f"
               IF ( I==1 .and. J==1 ) THEN
                  IF ( GET_MONTH() /= LASTMONTH ) THEN
                     CALL READ_JO1D( GET_MONTH() )
                     CALL READ_RESP( GET_MONTH() )
                     LASTMONTH = GET_MONTH()
                  ENDIF
               ENDIF

               ! Compute biogenic acetone emissions [atoms C/box/s]
               CALL EMISS_BIOACET( I,    J,    TMMP,  EMMO, 
     &                             EMIS, EMMB, GRASS, BIO_ACET )
               
               ! Also add ocean source of acetone [atoms C/box/s]
               CALL OCEAN_SOURCE_ACET( I, J, BIO_ACET )

               ! Add biogenic acetone to anthro source [atoms C/box/s]
               EMISRR(I,J,IDEACET) = EMISRR(I,J,IDEACET) + BIO_ACET
            ENDIF
!-----------------------------------------------------------------------------

            !=================================================================
            ! save biogenic isoprene emission for later use
            ! EMISRR has units [atoms C/box/s] 
            !=================================================================
            IF ( IDTISOP /= 0 ) THEN
               EMISRR(I,J,IDEISOP) = EMISRR(I,J,IDEISOP) + 
     &                               ( EMIS / DTSRCE )

               !- bey - for ND48 time serie archive 
               ! EMISRR48 has units [atoms C/cm2/s]
               EMISRR48(I,J) = EMIS / ( DTSRCE * AREA_CM2 )
            ENDIF
!------------------------------------------------------------------------------
!
!******************************************************************************
!  Biogenic source of CO -- from oxidation of METHANOL and MONOTERPENES
!
!  CO from METHANOL oxidation -- scaled from ISOPRENE (bnd, 1/2/01)
!
!    We need to scale the Isoprene flux to get the CH3OH (methanol) flux.
!    Currently, the annual isoprene flux in GEOS-CHEM is ~ 397 Tg C.
!
!    Daniel Jacob recommends a flux of 100 Tg/yr CO from CH3OH oxidation  
!    based on Singh et al. 2000 [JGR 105, 3795-3805] who estimate a global 
!    methanol source of 122 Tg yr-1, of which most (75 Tg yr-1) is 
!    "primary biogenic".  He also recommends for now that the CO flux 
!    from CH3OH oxidation be scaled to monthly mean isoprene flux.
!
!    To get CO from METHANOL oxidation, we must therefore multiply
!    the ISOPRENE flux by the following scale factor:
!      ( 100 Tg CO from CH3OH Oxidation  / 397 Tg C from Isoprene Flux ) *
!      (  12 g C/mole                    / 28 g CO/mole                )
!
!  CO from MONOTERPENE oxidation (bnd, bmy, 1/2/01)
!
!    Assume the production of CO from monoterpenes is instantaneous even 
!    though the lifetime of intermediate species may be on the order of hours 
!    or days.  This assumption will likely cause CO from monoterpene oxidation
!    to be too high in the box in which the monoterpene is emitted.
!
!    The CO yield here is taken from:
!
!    Hatakeyama et al. JGR, Vol. 96, p. 947-958 (1991)
!      "The ultimate yield of CO from the tropospheric oxidation of terpenes 
!       (including both O3 and OH reactions) was estimated to be 20% on the 
!       carbon number basis."  They studied ALPHA- & BETA-pinene.
!
!    Vinckier et al. Fresenius Env. Bull., Vol. 7, p.361-368 (1998)
!      "R(CO)=1.8+/-0.3" : 1.8/10 is about 20%.
!******************************************************************************
!
            !========================================================
            ! CO from MONOTERPENE oxidation [molec CO/box/s] 
            !========================================================
            TMPVAL            = ( EMMO / DTSRCE ) * 0.2d0
            EMISRR(I,J,IDECO) = EMISRR(I,J,IDECO) + TMPVAL
         
            ! ND29: CO-source from monoterpenes [molec/cm2/s]
            IF ( ND29 > 0 ) THEN
               AD29(I,J,5) = AD29(I,J,5) + ( TMPVAL / AREA_CM2 )
            ENDIF
!
!******************************************************************************
!  Biogenic source of PRPE -- scaled to ISOPRENE
!
!  Also, add biogenic emissions of alkenes. We do this by scaling to
!  isoprene emissions (probably OK for summertime conditions). The
!  scaling factor is based on work by Allen Goldstein. His values
!  indicate emission ratios of ethene:propene:butene=4:2:1 (on a
!  per molecule basis), with total emissions approx. equal to
!  10% of isoprene emissions (again, on molecule basis).
!  BIOSCAL is in units of atoms C (alkenes) / atoms C (isoprene)
!******************************************************************************
! Change this factor to exclude ethene (bey, ljm)
!    (10 molec alkenes / 100 molec isop) * (1 molec isop / 5 atoms C isop)
!    *(3 molec butene + propene / 7 molec total alkenes)
!    *(3.3333 atoms C but+prop mix/ 1 molec but+prop mix)
!    = 0.0286 atoms C butene+propene / atom C isop
! Note that 3.3333 atoms C/molecule is the weighted average for this mix.
!******************************************************************************
!
            BIOSCAL = 0.0286d0  ! new factor, (ljm, bey, 9/28/98)

            IF ( IDEPRPE /= 0 ) THEN
               EMISRR(I,J,IDEPRPE) = EMISRR(I,J,IDEPRPE) +
     &                               ( EMIS / DTSRCE ) * BIOSCAL
            ENDIF
!
!******************************************************************************
!  ND46 diagnostic: Biogenic emissions 
!
!     AD46(:,:,1) = Total biogenic ISOP  emissions [atoms C/cm2/s]
!     AD46(:,:,2) = Total biogenic ACET  emissions [atoms C/cm2/s]
!     AD46(:,:,3) = Total biogenic PRPE  emissions [atoms C/cm2/s]
!     AD46(:,:,4) = Total biogenic MONOT emissions [atoms C/cm2/s]
!
!  NOTES: 
!  (1 ) Now make ACET tracer #2 and PRPE tracer #3 (bmy, 9/13/01)
!  (2 ) Now archive ND46 as [atoms C/cm2/s] here (bmy, 9/13/01)
!******************************************************************************
!
            IF ( ND46 > 0 ) THEN

               ! ISOP emissions [atoms C/cm2/s] -- tracer #1
               AD46(I,J,1) = AD46(I,J,1) + ( EMIS / AREA_CM2 / DTSRCE )

               ! ACET emissions [atoms C/cm2/s] -- tracer #2
               AD46(I,J,2) = AD46(I,J,2) + ( BIO_ACET / AREA_CM2 )

               ! PRPE emissions [atoms C/cm2/s] -- tracer #3
               AD46(I,J,3) = AD46(I,J,3) + 
     &                      ( EMIS * BIOSCAL / AREA_CM2 / DTSRCE )

               ! Monoterpene emissions [atoms C/cm2/s] -- tracer #4
               AD46(I,J,4) = AD46(I,J,4) + ( EMMO / AREA_CM2 / DTSRCE ) 
            ENDIF           
         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE EMISSDR                                                

