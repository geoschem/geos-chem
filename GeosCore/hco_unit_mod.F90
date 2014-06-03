!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_unit_mod
!
! !DESCRIPTION: Module HCO_UNIT\_MOD contains routines to check/convert
! units. 
!
! \\
! !INTERFACE:
!
      MODULE HCO_UNIT_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD    

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: HCO_UNIT_CHANGE
      PUBLIC :: HCO_UNIT_GetMassScal
      PUBLIC :: HCO_UNIT_GetAreaScal
      PUBLIC :: HCO_UNIT_GetTimeScal
      PUBLIC :: HCO_UNIT_SCALCHECK 
!
! !PRIVATE MODULE VARIABLES:
!
      REAL(hp),  PARAMETER   :: N_0             = 6.022e+23_hp
      REAL(hp),  PARAMETER   :: SEC_IN_DAY      = 86400_hp
      REAL(hp),  PARAMETER   :: SEC_IN_LEAPYEAR = SEC_IN_DAY * 366_hp 
      REAL(hp),  PARAMETER   :: SEC_IN_REGYEAR  = SEC_IN_DAY * 365_hp
!
! !REVISION HISTORY:
!  15 May 2012 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_UNIT_CHANGE 
!
! !DESCRIPTION: Subroutine HCO_UNIT_CHANGE converts the values of the 
! passed array to units of (emitted) kg/m2/s. 
! The mass is in units of kg and refers to mass of emitted species. For
! most compounds, this corresponds to the molecular mass of the species,
! but e.g. for VOCs this can be mass of carbon instead.
! The mass and area/volume conversion is always performed, but the time 
! conversion is only done if a valid time string is provided. For example, 
! if the input unit is kg/cm3 it will be converted to kg/m3, while 
! ug/m2/year is converted to kg/m2/s. If no (valid) area/volume is given in 
! the unit string, the return flag PerArea is set to False (True 
! otherwise).\\ 
! The input argument UNITS refers to the unit of the input data.
! Argument MW\_IN denotes the molecular weight of the input unit
! (g/mol), while MW\_OUT is the molecular weight of the output unit.
! They can differ e.g. for VOCs whose output units are in mass carbon
! instead of mass species. The argument MOLEC\_RATIO is the coefficient 
! used for the conversion of molecules of species to molecules of carbon. 
! For example, MOLEC\_RATIO should be set to 2 for C2H6, which will 
! properly convert kg species (or molec species) to kg C. If the input 
! unit is already on a per carbon basis (e.g. kgC, or molecC), no species
! coefficients will be applied!
!\\
! Supported unit values:
!\\
! MASSES:
! - molec (includes molecC; molec(C); molec tracer; molecN, molec(N))
! - atom or atoms (incl. atomC, etc.)
! - kg, g, mg, ug, ng (incl. kgC, etc.)
!\\
! TIMES:
! s, sec, hr, hour, d, day, mt, month, y, year.
! Valid formats: /s, s-1, s^-1. 
!\\
! VOLUMES/AREAS: 
! cm2, m2, km2, cm3, dm3, m3, l.
! Valid formats: /cm3, cm-3, /cm^3, cm^-3.
!\\
! The following units will be ignored (no unit conversion is applied):
! - unitless, fraction, factor, hours, degC, 1
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_UNIT_CHANGE( ARRAY,       UNITS, MW_IN, MW_OUT, &
                                  MOLEC_RATIO, YYYY,  MM,    IsPerArea, RC )
!
! !USES:
!
      USE CHARPAK_MOD,              ONLY : CSTRIP
!
! !ARGUMENTS:
!
      REAL(sp),         POINTER         :: ARRAY(:,:,:,:) ! Data
      CHARACTER(LEN=*), INTENT(IN )     :: UNITS          ! Data unit
      REAL(hp),         INTENT(IN )     :: MW_IN          ! MW g/mol 
      REAL(hp),         INTENT(IN )     :: MW_OUT         ! MW g/mol
      REAL(hp),         INTENT(IN )     :: MOLEC_RATIO    ! molec. ratio
      INTEGER,          INTENT(IN )     :: YYYY           ! Data year 
      INTEGER,          INTENT(IN )     :: MM             ! Data month
      LOGICAL,          INTENT(OUT)     :: IsPerArea      ! Is per area? 
      INTEGER,          INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller - Initial Version
!  23 May 2013 - C. Keller - Now use additive method
!  01 Oct 2013 - C. Keller - Now convert to kg/m2/s instead of
!                            molec/cm2/s
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
      REAL(hp)              :: Fact, Coef1
      CHARACTER(LEN=31 )    :: unt

      !=================================================================
      ! HCO_UNIT_CHANGE begins here
      !=================================================================

      ! Init
      RC        = 0  
      Fact      = 1.0_hp
      Coef1     = 1.0_hp
      IsPerArea = .TRUE.

      ! Get input data unit and strip all blanks. 
      unt = TRIM(UNITS)
      CALL CSTRIP( unt )

      !=================================================================
      ! For special case that data is unitless, a fraction or any other
      ! quantity that shall not be converted - or if it's already in 
      ! units of kg/m2/s.
      !=================================================================
      IF ( HCO_UNIT_SCALCHECK(unt) < 2 ) RETURN

      !=================================================================
      ! Get scale factor for mass. Force to be a valid factor.
      !=================================================================
      Coef1 = HCO_UNIT_GetMassScal ( unt, MW_IN, MW_OUT, MOLEC_RATIO )
      IF ( Coef1 < 0.0_hp ) THEN
         WRITE(6,*) 'unrecognized mass unit: ', TRIM(unt)
         RC = -999
         RETURN 
      ENDIF
      Fact = Fact * Coef1

      !=================================================================
      ! Get scale factor for time. Skip if invalid factor. This makes
      ! sure that concentrations (e.g. kg/m3) are supported! 
      !=================================================================
      Coef1 = HCO_UNIT_GetTimeScal ( unt, MM, YYYY ) 
      IF ( Coef1 > 0.0_hp ) Fact = Fact * Coef1

      !=================================================================
      ! Get scale factor for area/volume. If no area conversion
      ! factor can be determined, set PerArea flag to False. 
      !=================================================================
      Coef1 = HCO_UNIT_GetAreaScal ( unt ) 
      IF ( Coef1 < 0.0_hp ) THEN
         IsPerArea = .FALSE.
      ELSE
         Fact = Fact * Coef1
      ENDIF

      ! Apply correction factor
      ARRAY(:,:,:,:) = ARRAY(:,:,:,:) * Fact

      ! Leave
      RC = 0 

      END SUBROUTINE HCO_UNIT_CHANGE
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_UNIT_GetMassScal
!
! !DESCRIPTION: Returns the mass scale factors for the given unit.
! This is the scale factor required to convert from unit 'Unit' to
! HEMCO units (i.e. kg).
!\\
!\\
! !INTERFACE:
!
      FUNCTION HCO_UNIT_GetMassScal ( unt, MW_IN, MW_OUT, MOLEC_RATIO ) &
         RESULT ( Scal ) 
!
! !USES:
!
      USE HCO_CHARTOOLS_MOD
!
! !INPUT PARAMETERS:
!
      CHARACTER(LEN=*), INTENT(IN   )   :: unt 
      REAL(hp),         INTENT(IN   )   :: MW_IN          ! MW g/mol 
      REAL(hp),         INTENT(IN   )   :: MW_OUT         ! MW g/mol
      REAL(hp),         INTENT(IN   )   :: MOLEC_RATIO    ! molec. ratio
!
! !OUTPUT PARAMETERS:
!
      REAL(hp)                         :: Scal
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!

      !=================================================================
      ! HCO_UNIT_GetMassScal begins here 
      !=================================================================

      ! Init
      Scal = -999.0_hp

      ! Molecules / atoms of carbon: convert to kg carbon
      ! Molecules / atoms of nitrogen: convert to kg output tracer
      IF ( IsInWord(unt,'molecC'  ) .OR. IsInWord(unt,'molecN')   .OR. &
           IsInWord(unt,'atomC'   ) .OR. IsInWord(unt,'atomsC' )  .OR. &
           IsInWord(unt,'molec(C)') .OR. IsInWord(unt,'atoms(C)') .OR. &
           IsInWord(unt,'atomN'   ) .OR. IsInWord(unt,'atomsN' )  .OR. &
           IsInWord(unt,'molec(N)') .OR. IsInWord(unt,'atoms(N)') .OR. &
           IsInWord(unt,'molectracer') ) THEN
         Scal = MW_OUT * 1e-3_hp / N_0

      ! Molecules / atoms of species: convert to kg output species.
      ELSEIF ( IsInWord(unt,'molec') .OR. IsInWord(unt,'atom') ) THEN
         Scal = MOLEC_RATIO * MW_OUT * 1e-3_hp / N_0

      ! Mols carbon / nitrogen of species
      ELSEIF ( IsInWord(unt,'nmolC') .OR. IsInWord(unt,'nmol(C)') .OR. &
               IsInWord(unt,'nmolN') .OR. IsInWord(unt,'nmol(N)') ) THEN
         Scal = 1e-9_hp * MW_OUT * 1e-3_hp
      ELSEIF ( IsInWord(unt,'umolC') .OR. IsInWord(unt,'umol(C)') .OR. &
               IsInWord(unt,'umolN') .OR. IsInWord(unt,'umol(N)') ) THEN
         Scal = 1e-6_hp * MW_OUT * 1e-3_hp 
      ELSEIF ( IsInWord(unt,'mmolC') .OR. IsInWord(unt,'mmol(C)') .OR. &
               IsInWord(unt,'mmolN') .OR. IsInWord(unt,'mmol(N)') ) THEN
         Scal = 1e-3_hp * MW_OUT * 1e-3_hp 
      ELSEIF ( IsInWord(unt,'molC') .OR. IsInWord(unt,'mol(C)')   .OR. &
               IsInWord(unt,'molN') .OR. IsInWord(unt,'mol(N)') ) THEN
         Scal = MW_OUT * 1e-3_hp 

      ! Mols of species
      ELSEIF ( IsInWord(unt,'nmol') ) THEN
         Scal = 1e-9_hp * MOLEC_RATIO * MW_OUT * 1e-3_hp
      ELSEIF ( IsInWord(unt,'umol') ) THEN
         Scal = 1e-6_hp * MOLEC_RATIO * MW_OUT * 1e-3_hp
      ELSEIF ( IsInWord(unt,'mmol') ) THEN
         Scal = 1e-3_hp * MOLEC_RATIO * MW_OUT * 1e-3_hp
      ELSEIF ( IsInWord(unt,'mol') ) THEN
         Scal = MOLEC_RATIO * MW_OUT * 1e-3_hp

      ! Mass Carbon of species
      ELSEIF ( IsInWord(unt,'ngC') .OR. IsInWord(unt,'ng(C)') ) THEN 
         Scal = 1e-12_hp / 12_hp * MW_OUT 
      ELSEIF ( IsInWord(unt,'ugC') .OR. IsInWord(unt,'ug(C)') ) THEN 
         Scal = 1e-9_hp / 12_hp * MW_OUT 
      ELSEIF ( IsInWord(unt,'mgC') .OR. IsInWord(unt,'mg(C)') ) THEN 
         Scal = 1e-6_hp / 12_hp * MW_OUT 
      ELSEIF ( IsInWord(unt,'kgC') .OR. IsInWord(unt,'kg(C)') ) THEN 
         Scal = 1.0_hp / 12_hp * MW_OUT 
      ELSEIF ( IsInWord(unt,'gC') .OR. IsInWord(unt,'g(C)') ) THEN 
         Scal = 1e-3_hp / 12_hp * MW_OUT

      ! Mass Nitrogen of species
      ELSEIF ( IsInWord(unt,'ngN') .OR. IsInWord(unt,'ng(N)') ) THEN 
         Scal = 1e-12_hp / 14_hp * MW_OUT 
      ELSEIF ( IsInWord(unt,'ugN') .OR. IsInWord(unt,'ug(N)') ) THEN 
         Scal = 1e-9_hp / 14_hp * MW_OUT 
      ELSEIF ( IsInWord(unt,'mgN') .OR. IsInWord(unt,'mg(N)') ) THEN 
         Scal = 1e-6_hp / 14_hp * MW_OUT 
      ELSEIF ( IsInWord(unt,'kgN') .OR. IsInWord(unt,'kg(N)') ) THEN 
         Scal = 1.0_hp / 14_hp * MW_OUT 
      ELSEIF ( IsInWord(unt,'gN') .OR. IsInWord(unt,'g(N)') ) THEN 
         Scal = 1e-3_hp / 14_hp * MW_OUT

      ! Mass of species
      ELSEIF ( IsInWord(unt,'ng') ) THEN
         Scal = 1e-12_hp * MOLEC_RATIO * MW_OUT / MW_IN 
      ELSEIF ( IsInWord(unt,'ug') ) THEN
         Scal = 1e-9_hp * MOLEC_RATIO * MW_OUT / MW_IN 
      ELSEIF ( IsInWord(unt,'mg') ) THEN
         Scal = 1e-6_hp * MOLEC_RATIO * MW_OUT / MW_IN 
      ELSEIF ( IsInWord(unt,'kg') ) THEN
         Scal = MOLEC_RATIO * MW_OUT / MW_IN 
      ELSEIF ( IsInWord(unt,'g') ) THEN
         Scal = 1e-3_hp * MOLEC_RATIO * MW_OUT / MW_IN 
      ENDIF

      END FUNCTION HCO_UNIT_GetMassScal 
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_UNIT_GetTimeScal
!
! !DESCRIPTION: Returns the time scale factors for the given unit.
! This is the scale factor required to convert from unit 'Unit' to
! HEMCO units (i.e. per second).
!\\
!\\
! !INTERFACE:
!
      FUNCTION HCO_UNIT_GetTimeScal ( unt, MM, YYYY ) &
         RESULT ( Scal ) 
!
! !USES:
!
      USE HCO_CHARTOOLS_MOD
!
! !INPUT PARAMETERS:
!
      CHARACTER(LEN=*), INTENT(IN   )   :: unt   ! This unit
      INTEGER,          INTENT(IN   )   :: MM    ! Current month
      INTEGER,          INTENT(IN   )   :: YYYY  ! Current year 
!
! !OUTPUT PARAMETERS:
!
      REAL(hp)                         :: Scal
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!
      INTEGER  :: Month, Year
      LOGICAL  :: IS_LEAPYEAR
      INTEGER  :: MONTHDAYS(12) = (/ 31, 28, 31, 30, 31, 30, &
                                     31, 31, 30, 31, 30, 31   /)

      !=================================================================
      ! HCO_UNIT_GetTimeScal begins here 
      !=================================================================

      ! Init
      Scal = -999.0_hp

      ! Is this a leap year?
      Year  = MAX(YYYY,1)
      IS_LEAPYEAR = ( (MOD(Year,4) == 0) .AND. (MOD(Year,400) /= 0) )
      IF ( IS_LEAPYEAR ) MONTHDAYS(2) = 29

      ! second
      IF ( IsInWord(unt,'/s')   .OR. IsInWord(unt,'s-1') .OR. &
           IsInWord(unt,'s^-1') .OR. IsInWord(unt,'sec') ) THEN
         Scal = 1.0_hp

      ! hour
      ELSEIF ( IsInWord(unt,'hr') .OR. IsInWord(unt,'hour') ) THEN
         Scal = 1.0_hp / 3600_hp

      ! day
      ELSEIF ( IsInWord(unt,'/d')   .OR. IsInWord(unt,'d-1') .OR. &
               IsInWord(unt,'d^-1') .OR. IsInWord(unt,'day') ) THEN
         Scal = 1.0_hp / SEC_IN_DAY

      ! month
      ELSEIF ( IsInWord(unt,'mt') .OR. IsInWord(unt,'month') ) THEN
         Month = MAX(MM,1)
         Scal  = 1.0_hp / MONTHDAYS(Month) / SEC_IN_DAY 
      
      ! year
      ELSEIF ( IsInWord(unt,'/y')   .OR. IsInWord(unt,'y-1') .OR. &
               IsInWord(unt,'y^-1') .OR. IsInWord(unt,'yr')  .OR. &
               IsInWord(unt,'year') ) THEN

         IF ( IS_LEAPYEAR ) THEN 
            Scal = 1.0_hp / SEC_IN_LEAPYEAR 
         ELSE
            Scal = 1.0_hp / SEC_IN_REGYEAR 
         ENDIF
      ENDIF

      END FUNCTION HCO_UNIT_GetTimeScal 
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_UNIT_GetAreaScal
!
! !DESCRIPTION: Returns the area/volume scale factors for the given unit.
! This is the scale factor required to convert from unit 'Unit' to
! HEMCO units (i.e. per m2 or per m3).
!\\
!\\
! !INTERFACE:
!
      FUNCTION HCO_UNIT_GetAreaScal ( unt ) & 
         RESULT ( Scal ) 
!
! !USES:
!
      USE HCO_CHARTOOLS_MOD
!
! !INPUT PARAMETERS:
!
      CHARACTER(LEN=*), INTENT(IN   )   :: unt   ! This unit
!
! !OUTPUT PARAMETERS:
!
      REAL(hp)                         :: Scal
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!

      !=================================================================
      ! HCO_UNIT_GetAreaScal begins here 
      !=================================================================

      ! Init
      Scal = -999.0_hp

      ! cm2
      IF ( IsInWord(unt,'/cm2' ) .OR. IsInWord(unt,'cm-2' ) .OR. &
           IsInWord(unt,'/cm^2') .OR. IsInWord(unt,'cm^-2')       ) THEN
         Scal = 1.0_hp / 1e-4_hp

      ! km2
      ELSEIF ( IsInWord(unt,'/km2' ) .OR. IsInWord(unt,'km-2' ) .OR. & 
               IsInWord(unt,'/km^2') .OR. IsInWord(unt,'km^-2')       ) THEN 
         Scal = 1e6_hp

      ! m2
      ELSEIF ( IsInWord(unt,'/m2' ) .OR. IsInWord(unt,'m-2' ) .OR. & 
               IsInWord(unt,'/m^2') .OR. IsInWord(unt,'m^-2')       ) THEN 
         Scal = 1.0_hp

      !=================================================================
      ! Convert volume to m3 
      !=================================================================

      ! cm3
      ELSEIF ( IsInWord(unt,'/cm3')  .OR. IsInWord(unt,'cm-3' ) .OR. &
               IsInWord(unt,'/cm^3') .OR. IsInWord(unt,'cm^-3')       ) THEN 
         Scal = 1e-6_hp

      ! dm3
      ELSEIF ( IsInWord(unt,'/dm3')  .OR. IsInWord(unt,'dm-3' ) .OR. & 
               IsInWord(unt,'/dm^3') .OR. IsInWord(unt,'dm^-3')       ) THEN 
         Scal = 1e-3_hp

      ! m3
      ELSEIF ( IsInWord(unt,'/m3')  .OR. IsInWord(unt,'m-3' ) .OR. & 
               IsInWord(unt,'/m^3') .OR. IsInWord(unt,'m^-3')       ) THEN 
         Scal = 1.0_hp

      ! L
      ELSEIF ( IsInWord(unt,'/l') .OR. IsInWord(unt,'l-1') .OR. & 
               IsInWord(unt,'l^-1')                              ) THEN 
         Scal = 1e-3_hp
      ENDIF

      END FUNCTION HCO_UNIT_GetAreaScal 
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_UNIT_SCALCHECK 
!
! !DESCRIPTION: Check if the provided unit is unitless. Returns
! 0 if Unit is unitless, 1 if it's not unitless but in correct
! HEMCO units (i.e. kg/m2/s), 2 otherwise. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION HCO_UNIT_SCALCHECK ( Unit ) Result ( Flag ) 
!
! !USES:
!
      USE CHARPAK_MOD, ONLY : TRANLC 
!
! !INPUT PARAMETERS:
!
      CHARACTER(LEN=*), INTENT(IN   )  :: Unit 
!
! !OUTPUT PARAMETERS:
!
      INTEGER                          :: Flag  ! 0=ok, 1=warning, 2=error
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!
      CHARACTER(LEN=31) :: tmpU

      !=================================================================
      ! HCO_UNIT_SCALCHECK begins here 
      !=================================================================

      ! Mirror
      tmpU = Unit

      ! lower case
      CALL TRANLC( tmpU )

      ! Error (default):
      Flag = 2

      ! Ok:
      IF ( TRIM(tmpU) == 'unitless' .OR. &
           TRIM(tmpU) == 'fraction' .OR. &
           TRIM(tmpU) == 'factor'   .OR. &
           TRIM(tmpU) == 'scale'    .OR. &
           TRIM(tmpU) == 'hours'    .OR. &
           TRIM(tmpU) == 'm2/m2'    .OR. &
           TRIM(tmpU) == '1'            ) THEN
         Flag = 0

      ! Warning:
      ELSEIF ( TRIM(tmpU) == 'kg/m2/s'     .OR. &
               TRIM(tmpU) == 'kgc/m2/s'    .OR. &
               TRIM(tmpU) == 'kg(c)/m2/s'        ) THEN 
         Flag = 1
      ENDIF

      END FUNCTION HCO_UNIT_SCALCHECK 
!EOC
      END MODULE HCO_UNIT_MOD
!EOM
