!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !!
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_unit_mod.F90
!
! !DESCRIPTION: Module HCO\_Unit\_Mod contains routines to check/convert
! units.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_Unit_Mod
!
! !USES:
!
  USE HCO_Error_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCO_Unit_Change
  PUBLIC :: HCO_Unit_GetMassScal
  PUBLIC :: HCO_Unit_GetAreaScal
  PUBLIC :: HCO_Unit_GetTimeScal
  PUBLIC :: HCO_Unit_ScalCheck
  PUBLIC :: HCO_IsUnitLess
  PUBLIC :: HCO_IsIndexData
  PUBLIC :: HCO_UnitTolerance
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: HCO_Unit_Change_SP
  PRIVATE :: HCO_Unit_Change_DP
!
! !REVISION HISTORY:
!  15 May 2012 - C. Keller   - Initialization
!  08 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  08 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  24 Jul 2014 - C. Keller   - Now define unitless & 'standard' units as
!                              parameter
!  13 Aug 2014 - C. Keller   - Interface for sp & dp arrays
!  13 Mar 2015 - R. Yantosca - Add m and m2 to the "unitless" list
!  16 Mar 2015 - R. Yantosca - Also allow "kg m-2 s-1" and similar units
!  16 Mar 2015 - R. Yantosca - Add dobsons and dobsons/day units
!  16 Jun 2015 - R. Yantosca - Add % and percent to the unitless list
!  07 Jan 2016 - E. Lundgren - Update Avogadro's # to NIST 2014 value
!  19 Sep 2016 - R. Yantosca - Make sure all strings are the same length in
!                              the array constructor or Gfortran will choke
!  10 Jun 2018 - C. Keller   - Add mol/mol to unitless quantities
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  REAL(dp),  PARAMETER :: N_0             = 6.022140857e+23_dp
  REAL(hp),  PARAMETER :: SEC_IN_DAY      = 86400_hp
  REAL(hp),  PARAMETER :: SEC_IN_LEAPYEAR = SEC_IN_DAY * 366_hp
  REAL(hp),  PARAMETER :: SEC_IN_REGYEAR  = SEC_IN_DAY * 365_hp

  ! Accepted units for unitless data. No unit conversion is applied to
  ! data with any of these units. The first entry represents the
  ! character that denotes unitless data in the HEMCO configuration
  ! file. The second entry represents the character that denotes
  ! index data. Different regridding algorithms are applied to
  ! index data compared to unitless data.
  ! All units listed below will not be converted by HEMCO. You can
  ! add more units if you don't want HEMCO to attempt to convert data
  ! in these units.
  ! All characters in this list should be lower case!
  INTEGER,           PARAMETER :: NUL = 38
  CHARACTER(LEN=15), PARAMETER :: UL(NUL) = (/ '1          ',   &
                                               'count      ',   &
                                               'unitless   ',   &
                                               'fraction   ',   &
                                               'factor     ',   &
                                               'scale      ',   &
                                               'hours      ',   &
                                               'mol/mol    ',   &
                                               'v/v        ',   &
                                               'v/v/s      ',   &
                                               's-1        ',   &
                                               's^-1       ',   &
                                               'm2/m2      ',   &
                                               'm2m-2      ',   &
                                               'kg/kg      ',   &
                                               'kgkg-1     ',   &
                                               'mg/m3      ',   &
                                               'mg/m2/d    ',   &
                                               'k          ',   &
                                               'w/m2       ',   &
                                               'wm-2       ',   &
                                               'pptv       ',   &
                                               'ppt        ',   &
                                               'ppbv       ',   &
                                               'ppb        ',   &
                                               'ppmv       ',   &
                                               'ppm        ',   &
                                               'm/s        ',   &
                                               'ms-1       ',   &
                                               'm          ',   &
                                               'cm2cm-2    ',   &
                                               'dobsons    ',   &
                                               'dobsons/day',   &
                                               'DU         ',   &
                                               'pa         ',   &
                                               'hpa        ',   &
                                               '%          ',   &
                                               'percent    '    /)

  ! Accepted units for data on HEMCO standard units. No unit conversion
  ! is applied to data with any of these units.
  ! All characters in this list should be lower case!

  ! Emission units
  INTEGER,           PARAMETER :: NHE = 6
  CHARACTER(LEN=15), PARAMETER :: HE(NHE) = (/ 'kg/m2/s    ',   &
                                               'kgc/m2/s   ',   &
                                               'kg(c)/m2/s ',   &
                                               'kgm-2s-1   ',   &
                                               'kgcm-2s-1  ',   &
                                               'kg(c)m-2s-1'  /)
  ! Concentration units
  INTEGER,           PARAMETER :: NHC = 3
  CHARACTER(LEN=15), PARAMETER :: HC(NHC) = (/ 'kg/m3 ',        &
                                               'kgm-3 ',        &
                                               'kgm^-3'       /)

  ! Interfaces:
  INTERFACE HCO_UNIT_CHANGE
     MODULE PROCEDURE HCO_UNIT_CHANGE_SP
     MODULE PROCEDURE HCO_UNIT_CHANGE_DP
  END INTERFACE HCO_UNIT_CHANGE

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Unit_Change_sp
!
! !DESCRIPTION: Subroutine HCO\_UNIT\_CHANGE\_SP is a wrapper routine to
! convert the values of the passed single precision array to units of
! (emitted) kg/m2/s.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Unit_Change_SP( HcoConfig,   ARRAY, UNITS, MW_IN, MW_OUT,   &
                                 MOLEC_RATIO, YYYY,  MM,    AreaFlag, &
                                 TimeFlag,    FACT,  RC,    KeepSpec   )
!
! !USES:
!
    USE HCO_TYPES_MOD,    ONLY : ConfigObj
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER                 :: HcoConfig
    CHARACTER(LEN=*), INTENT(IN )             :: UNITS          ! Data unit
    REAL(hp),         INTENT(IN )             :: MW_IN          ! MW g/mol
    REAL(hp),         INTENT(IN )             :: MW_OUT         ! MW g/mol
    REAL(hp),         INTENT(IN )             :: MOLEC_RATIO    ! molec. ratio
    INTEGER,          INTENT(IN )             :: YYYY           ! Data year
    INTEGER,          INTENT(IN )             :: MM             ! Data month
    LOGICAL,          INTENT(IN ), OPTIONAL   :: KeepSpec       ! Keep input species?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: AreaFlag       ! 2 if per area, 3 if per volume, 0 otherwise
    INTEGER,          INTENT(INOUT)           :: TimeFlag       ! 1 if per time, 0 otherwise
    REAL(sp),         POINTER                 :: ARRAY(:,:,:,:) ! Data
    INTEGER,          INTENT(INOUT)           :: RC
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),         INTENT(  OUT), OPTIONAL :: FACT           ! Applied factor
!
! !REVISION HISTORY:
!  13 Aug 2014 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    REAL(hp)     :: Factor

    !=================================================================
    ! HCO_UNIT_CHANGE_SP begins here
    !=================================================================

    CALL HCO_Unit_Factor( HcoConfig,UNITS, MW_IN, MW_OUT, MOLEC_RATIO, YYYY, &
                          MM, AreaFlag, TimeFlag, Factor, RC, KeepSpec=KeepSpec   )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Apply correction factor
    IF ( Factor /= 1.0_hp ) THEN
       ARRAY(:,:,:,:) = ARRAY(:,:,:,:) * Factor
    ENDIF

    ! Eventually return factor
    IF ( PRESENT(FACT) ) FACT = Factor

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_Unit_Change_SP
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Unit_Change_dp
!
! !DESCRIPTION: Subroutine HCO\_UNIT\_CHANGE\_DP is a wrapper routine to
! convert the values of the passed double precision array to units of
! (emitted) kg/m2/s.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Unit_Change_DP( HcoConfig,   ARRAY, UNITS, MW_IN, MW_OUT,   &
                                 MOLEC_RATIO, YYYY,  MM,    AreaFlag, &
                                 TimeFlag,    FACT,  RC,    KeepSpec   )
!
! !USES:
!
    USE HCO_TYPES_MOD,    ONLY : ConfigObj
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER                 :: HcoConfig
    CHARACTER(LEN=*), INTENT(IN )             :: UNITS          ! Data unit
    REAL(hp),         INTENT(IN )             :: MW_IN          ! MW g/mol
    REAL(hp),         INTENT(IN )             :: MW_OUT         ! MW g/mol
    REAL(hp),         INTENT(IN )             :: MOLEC_RATIO    ! molec. ratio
    INTEGER,          INTENT(IN )             :: YYYY           ! Data year
    INTEGER,          INTENT(IN )             :: MM             ! Data month
    LOGICAL,          INTENT(IN ), OPTIONAL   :: KeepSpec       ! Keep input species?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: AreaFlag       ! 2 if per area, 3 if per volume, 0 otherwise
    INTEGER,          INTENT(INOUT)           :: TimeFlag       ! 1 if per time, 0 otherwise
    REAL(dp),         POINTER                 :: ARRAY(:,:,:,:) ! Data
    INTEGER,          INTENT(INOUT)           :: RC
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),         INTENT(  OUT), OPTIONAL :: FACT           ! Applied factor
!
! !REVISION HISTORY:
!  13 Aug 2014 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    REAL(hp)     :: Factor

    !=================================================================
    ! HCO_UNIT_CHANGE_DP begins here
    !=================================================================

    CALL HCO_Unit_Factor( HcoConfig, UNITS, MW_IN, MW_OUT, MOLEC_RATIO, YYYY, &
                          MM, AreaFlag, TimeFlag, Factor, RC, KeepSpec=KeepSpec )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Apply correction factor
    IF ( Factor /= 1.0_hp ) THEN
       ARRAY(:,:,:,:) = ARRAY(:,:,:,:) * Factor
    ENDIF

    ! Eventually return factor
    IF ( PRESENT(FACT) ) FACT = Factor

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_Unit_Change_DP
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Unit_Factor
!
! !DESCRIPTION: Subroutine HCO\_UNIT\_Factor calculates the conversion
! factor needed to conver unit UNIT to HEMCO units.
!\\
!\\
! The mass is in units of kg and refers to mass of emitted species. For
! most compounds, this corresponds to the molecular mass of the species,
! but e.g. for VOCs this can be mass of carbon instead.
! The mass and area/volume conversion is always performed, but the time
! conversion is only done if a valid time string is provided. For example,
! if the input unit is kg/cm3 it will be converted to kg/m3, while
! ug/m2/year is converted to kg/m2/s. If no (valid) area/volume is given in
! the unit string, the return flag PerArea is set to False (True
! otherwise).
!\\
!\\
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
!\\
! Supported unit values:
!\\
!\\
! MASSES:
! \begin{itemize}
! \item molec (includes molecC; molec(C); molec tracer; molecN, molec(N))
! \item atom or atoms (incl. atomC, etc.)
! \item kg, g, mg, ug, ng (incl. kgC, etc.)
! \end{itemize}
!
! TIMES:
! \begin{itemize}
! \item s, sec, hr, hour, d, day, mt, month, y, year.
! \item Valid formats: \/s, s-1, s\^-1.
! \end{itemize}
!
! VOLUMES/AREAS:
! \begin{itemize}
! \item cm2, m2, km2, cm3, dm3, m3, l.
! \item Valid formats: \/cm3, cm-3, \/cm\^3, cm\^-3.
! \end{itemize}
!
! The following units will be ignored (no unit conversion is applied):
! \begin{itemize}
! \item unitless
! \item fraction
! \item factor
! \item hours
! \item degC
! \item 1
! \end{itemize}
!
! !INTERFACE:
!
  SUBROUTINE HCO_Unit_Factor( HcoConfig, UNITS, MW_IN,    MW_OUT,   MOLEC_RATIO, YYYY, &
                              MM,       AreaFlag, TimeFlag, Factor,      RC, KeepSpec )
!
! !USES:
!
    USE HCO_TYPES_MOD,    ONLY : ConfigObj
    USE CharPak_Mod,      ONLY : CStrip
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER               :: HcoConfig
    CHARACTER(LEN=*), INTENT(IN )           :: UNITS       ! Data unit
    REAL(hp),         INTENT(IN )           :: MW_IN       ! MW g/mol
    REAL(hp),         INTENT(IN )           :: MW_OUT      ! MW g/mol
    REAL(hp),         INTENT(IN )           :: MOLEC_RATIO ! molec. ratio
    INTEGER,          INTENT(IN )           :: YYYY        ! Data year
    INTEGER,          INTENT(IN )           :: MM          ! Data month
    LOGICAL,          INTENT(IN ), OPTIONAL :: KeepSpec    ! Keep input species?
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)         :: AreaFlag
    INTEGER,          INTENT(  OUT)         :: TimeFlag
    REAL(hp),         INTENT(  OUT)         :: Factor         ! Conversion factor
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)         :: RC
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller - Initial Version
!  23 May 2013 - C. Keller - Now use additive method
!  01 Oct 2013 - C. Keller - Now convert to kg/m2/s instead of
!                            molec/cm2/s
!  13 Aug 2014 - C. Keller - Split off from subroutine HCO\_UNIT\_CHANGE
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER                       :: FLAG, CHECK
    REAL(hp)                      :: Coef1
    CHARACTER(LEN=31 )            :: unt
    CHARACTER(LEN=255)            :: MSG
    CHARACTER(LEN=255), PARAMETER :: LOC = 'HCO_UNIT_FACTOR (HCO_UNIT_MOD.F90)'

    !=================================================================
    ! HCO_UNIT_FACTOR begins here
    !=================================================================

    ! Init
    RC        = 0
    Factor    = 1.0_hp
    Coef1     = 1.0_hp

    ! Get input data unit and strip all blanks.
    unt = TRIM(UNITS)
    CALL CSTRIP( unt )

    !=================================================================
    ! For special case that data is unitless, a fraction or any other
    ! quantity that shall not be converted - or if it's already in
    ! units of kg/m2/s.
    !=================================================================
    CHECK = HCO_UNIT_SCALCHECK(unt)

    ! unitless data
    IF ( CHECK == 0 ) THEN

       AreaFlag = -1
       TimeFlag = -1
       RETURN

!    ! emissions
!    ELSEIF ( CHECK == 1 ) THEN
!       AreaFlag = 2
!       TimeFlag = 1
!       RETURN
!
!    ! concentrations
!    ELSEIF ( CHECK == 2 ) THEN
!       AreaFlag = 3
!       TimeFlag = 0
!       RETURN

    ENDIF

    !=================================================================
    ! Get scale factor for mass. Force to be a valid factor.
    !=================================================================
    CALL HCO_UNIT_GetMassScal ( HcoConfig, unt, MW_IN, MW_OUT, &
                                MOLEC_RATIO, Coef1, KeepSpec=KeepSpec )

    IF ( Coef1 < 0.0_hp ) THEN
       MSG = 'cannot do unit conversion. Mass unit: ' // TRIM(unt)
       CALL HCO_ERROR ( HcoConfig%Err, MSG, RC, ThisLoc = LOC )
       RETURN
    ENDIF
    Factor = Factor * Coef1

    !=================================================================
    ! Get scale factor for time. Skip if invalid factor. This makes
    ! sure that concentrations (e.g. kg/m3) are supported!
    !=================================================================
    CALL HCO_UNIT_GetTimeScal ( unt, MM, YYYY, Coef1, Flag )
    IF ( Flag > 0 ) Factor   = Factor * Coef1
    TimeFlag = Flag

    !=================================================================
    ! Get scale factor for area/volume. If no area conversion
    ! factor can be determined, set PerArea flag to False.
    !=================================================================
    CALL HCO_UNIT_GetAreaScal ( unt, Coef1, Flag )
    IF ( Flag > 0 ) Factor = Factor * Coef1
    AreaFlag = Flag

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_Unit_Factor
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Unit_GetMassScal
!
! !DESCRIPTION: Returns the mass scale factors for the given unit.
! This is the scale factor required to convert from input units to
! HEMCO units (i.e. kg). If KeepSpec is set to true, the molecular
! weight of the input data is preserved, e.g. data in kgC is kept in
! kgC, data in molecC is converted to kgC, etc.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_UNIT_GetMassScal( HcoConfig, unt, MW_IN, MW_OUT, &
                                   MOLEC_RATIO, Scal, KeepSpec )
!
! !USES:
!
    USE HCO_CharTools_Mod
    USE HCO_TYPES_MOD,    ONLY : ConfigObj
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER               :: HcoConfig
    CHARACTER(LEN=*), INTENT(IN)            :: unt         ! Input units
    REAL(hp),         INTENT(IN)            :: MW_IN       ! MW g/mol
    REAL(hp),         INTENT(IN)            :: MW_OUT      ! MW g/mol
    REAL(hp),         INTENT(IN)            :: MOLEC_RATIO ! molec. ratio
    LOGICAL,          INTENT(IN ), OPTIONAL :: KeepSpec    ! Keep input species?
!
! !OUTPUT PARAMETER:
!
    REAL(hp),         INTENT(OUT)           :: Scal        ! Scale factor
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!  17 Sep 2014 - C. Keller - Now accept '1' as mass unit (e.g. 1/m3/s)
!EOP
!------------------------------------------------------------------------------
!BOC
    LOGICAL             :: KEEP
    REAL(hp)            :: MW_N, MW_C
    REAL(hp)            :: MR
    CHARACTER(LEN=255)  :: MSG

    !=================================================================
    ! HCO_UNIT_GetMassScal begins here
    !=================================================================

    ! Init
    Scal = -999.0_hp

    ! Keep same species
    IF ( PRESENT(KeepSpec) ) THEN
       KEEP = KeepSpec
    ELSE
       KEEP = .FALSE.
    ENDIF

    ! Set molecular weights of N and C
    MW_N = MW_OUT
    MW_C = MW_OUT

    ! Mass unit of '1': keep as is. Note: this has to be the first
    ! entry of the unit string (e.g. 1/m2/s or 1 cm^-3). Don't
    ! accept any other ones (e.g. kg m^-2 s^-1)!
    IF ( unt(1:1) == '1' ) THEN
       Scal = 1.0_hp
       RETURN
    ENDIF

    ! Error checks: all passed variables must be defined! If this is not
    ! the case, only accept kg as valid unit!
    IF ( MW_IN       <= 0.0_hp .OR. &
         MW_OUT      <= 0.0_hp .OR. &
         MOLEC_RATIO <= 0.0_hp       ) THEN
       IF ( IsInWord(unt,'kg') .OR. KEEP ) THEN
          Scal = 1.0_hp
          RETURN
       ELSE
          MSG = 'Cannot determine unit conversion factor for mass - ' // &
                'not all species parameter are defined!'
          CALL HCO_MSG(HcoConfig%Err,MSG)
          RETURN
       ENDIF
    ENDIF

    ! Molecules / atoms of carbon: convert to kg carbon
    ! Molecules / atoms of nitrogen: convert to kg output tracer
    IF ( IsInWord(unt,'molecC'  ) .OR. IsInWord(unt,'atomC' )   .OR. &
         IsInWord(unt,'atomsC'  ) .OR. IsInWord(unt,'molec(C)') .OR. &
         IsInWord(unt,'atoms(C)') ) THEN
       IF ( KEEP ) THEN
          Scal = MW_C * 1e-3_hp / N_0
       ELSE
          Scal = MW_OUT * 1e-3_hp / N_0
       ENDIF

    ELSEIF ( IsInWord(unt,'molecN') .OR. IsInWord(unt,'atomN') .OR. &
             IsInWord(unt,'atomsN') .OR. IsInWord(unt,'molec(N)') .OR. &
             IsInWord(unt,'atoms(N)') ) THEN
       IF ( KEEP ) THEN
          Scal = MW_N * 1e-3_hp / N_0
       ELSE
          Scal = MW_OUT * 1e-3_hp / N_0
       ENDIF

    ! Molecules / atoms of species: convert to kg output species.
    ELSEIF ( IsInWord(unt,'molec') .OR. IsInWord(unt,'atom') ) THEN
       Scal = MOLEC_RATIO * MW_OUT * 1e-3_hp / N_0

    ! Mols carbon / nitrogen of species
    ELSEIF ( IsInWord(unt,'nmolC') .OR. IsInWord(unt,'nmol(C)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-9_hp * MW_C * 1e-3_hp
       ELSE
          Scal = 1e-9_hp * MW_OUT * 1e-3_hp
       ENDIF
    ELSEIF ( IsInWord(unt,'nmolN') .OR. IsInWord(unt,'nmol(N)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-9_hp * MW_N * 1e-3_hp
       ELSE
          Scal = 1e-9_hp * MW_OUT * 1e-3_hp
       ENDIF
    ELSEIF ( IsInWord(unt,'umolC') .OR. IsInWord(unt,'umol(C)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-6_hp * MW_C * 1e-3_hp
       ELSE
          Scal = 1e-6_hp * MW_OUT * 1e-3_hp
       ENDIF
    ELSEIF ( IsInWord(unt,'umolN') .OR. IsInWord(unt,'umol(N)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-6_hp * MW_N * 1e-3_hp
       ELSE
          Scal = 1e-6_hp * MW_OUT * 1e-3_hp
       ENDIF
    ELSEIF ( IsInWord(unt,'mmolC') .OR. IsInWord(unt,'mmol(C)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-3_hp * MW_C * 1e-3_hp
       ELSE
          Scal = 1e-3_hp * MW_OUT * 1e-3_hp
       ENDIF
    ELSEIF ( IsInWord(unt,'mmolN') .OR. IsInWord(unt,'mmol(N)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-3_hp * MW_N * 1e-3_hp
       ELSE
          Scal = 1e-3_hp * MW_OUT * 1e-3_hp
       ENDIF
    ELSEIF ( IsInWord(unt,'molC') .OR. IsInWord(unt,'mol(C)') ) THEN
       IF ( KEEP ) THEN
          Scal = MW_C * 1e-3_hp
       ELSE
          Scal = MW_OUT * 1e-3_hp
       ENDIF
    ELSEIF ( IsInWord(unt,'molN') .OR. IsInWord(unt,'mol(N)') ) THEN
       IF ( KEEP ) THEN
          Scal = MW_N * 1e-3_hp
       ELSE
          Scal = MW_OUT * 1e-3_hp
       ENDIF

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
       IF ( KEEP ) THEN
          Scal = 1e-12_hp
       ELSE
          Scal = 1e-12_hp / 12_hp * MW_OUT
       ENDIF
    ELSEIF ( IsInWord(unt,'ugC') .OR. IsInWord(unt,'ug(C)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-9_hp
       ELSE
          Scal = 1e-9_hp / 12_hp * MW_OUT
       ENDIF
    ELSEIF ( IsInWord(unt,'mgC') .OR. IsInWord(unt,'mg(C)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-6_hp
       ELSE
          Scal = 1e-6_hp / 12_hp * MW_OUT
       ENDIF
    ELSEIF ( IsInWord(unt,'kgC') .OR. IsInWord(unt,'kg(C)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1.0_hp
       ELSE
          Scal = 1.0_hp / 12_hp * MW_OUT
       ENDIF
    ELSEIF ( IsInWord(unt,'gC') .OR. IsInWord(unt,'g(C)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-3_hp
       ELSE
          Scal = 1e-3_hp / 12_hp * MW_OUT
       ENDIF

    ! Mass Nitrogen of species
    ELSEIF ( IsInWord(unt,'ngN') .OR. IsInWord(unt,'ng(N)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-12_hp
       ELSE
          Scal = 1e-12_hp / 14_hp * MW_OUT
       ENDIF
    ELSEIF ( IsInWord(unt,'ugN') .OR. IsInWord(unt,'ug(N)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-9_hp
       ELSE
          Scal = 1e-9_hp / 14_hp * MW_OUT
       ENDIF
    ELSEIF ( IsInWord(unt,'mgN') .OR. IsInWord(unt,'mg(N)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-6_hp
       ELSE
          Scal = 1e-6_hp / 14_hp * MW_OUT
       ENDIF
    ELSEIF ( IsInWord(unt,'kgN') .OR. IsInWord(unt,'kg(N)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1.0_hp
       ELSE
          Scal = 1.0_hp / 14_hp * MW_OUT
       ENDIF
    ELSEIF ( IsInWord(unt,'gN') .OR. IsInWord(unt,'g(N)') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-3_hp
       ELSE
          Scal = 1e-3_hp / 14_hp * MW_OUT
       ENDIF

    ! Mass of species
    ELSEIF ( IsInWord(unt,'ng') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-12_hp
       ELSE
          Scal = 1e-12_hp * MOLEC_RATIO * MW_OUT / MW_IN
       ENDIF
    ELSEIF ( IsInWord(unt,'ug') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-9_hp
       ELSE
          Scal = 1e-9_hp * MOLEC_RATIO * MW_OUT / MW_IN
       ENDIF
    ELSEIF ( IsInWord(unt,'mg') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-6_hp
       ELSE
          Scal = 1e-6_hp * MOLEC_RATIO * MW_OUT / MW_IN
       ENDIF
    ELSEIF ( IsInWord(unt,'kg') ) THEN
       IF ( KEEP ) THEN
          Scal = 1.0_hp
       ELSE
          Scal = MOLEC_RATIO * MW_OUT / MW_IN
       ENDIF
    ELSEIF ( IsInWord(unt,'g') ) THEN
       IF ( KEEP ) THEN
          Scal = 1e-3_hp
       ELSE
          Scal = 1e-3_hp * MOLEC_RATIO * MW_OUT / MW_IN
       ENDIF
    ENDIF

  END SUBROUTINE HCO_Unit_GetMassScal
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Unit_GetTimeScal
!
! !DESCRIPTION: Returns the time scale factors for the given unit.
! This is the scale factor required to convert from unit 'Unit' to
! HEMCO units (i.e. per second).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Unit_GetTimeScal( unt, MM, YYYY, Scal, Flag )
!
! !USES:
!
    USE HCO_CharTools_Mod
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: unt   ! This unit
    INTEGER,          INTENT(IN)  :: MM    ! Current month
    INTEGER,          INTENT(IN)  :: YYYY  ! Current year
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),         INTENT(OUT) :: Scal  ! Scale factor
    INTEGER,          INTENT(OUT) :: Flag  ! 1=per time, 0 otherwise
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
    Flag = 0

    ! Is this a leap year?
    Year  = MAX(YYYY,1)
    IS_LEAPYEAR = ( (MOD(Year,4) == 0) .AND. (MOD(Year,400) /= 0) )
    IF ( IS_LEAPYEAR ) MONTHDAYS(2) = 29

    ! second
    IF ( IsInWord(unt,'/s')   .OR. IsInWord(unt,'s-1') .OR. &
         IsInWord(unt,'s^-1') .OR. IsInWord(unt,'sec') ) THEN
       Scal = 1.0_hp
       Flag = 1

    ! hour
    ELSEIF ( IsInWord(unt,'hr') .OR. IsInWord(unt,'hour') ) THEN
       Scal = 1.0_hp / 3600_hp
       Flag = 1

    ! day
    ELSEIF ( IsInWord(unt,'/d')   .OR. IsInWord(unt,'d-1') .OR. &
             IsInWord(unt,'d^-1') .OR. IsInWord(unt,'day') ) THEN
       Scal = 1.0_hp / SEC_IN_DAY
       Flag = 1

      ! month
    ELSEIF ( IsInWord(unt,'mt') .OR. IsInWord(unt,'month') ) THEN
       Month = MAX(MM,1)
       Scal  = 1.0_hp / MONTHDAYS(Month) / SEC_IN_DAY
       Flag  = 1

    ! year
    ELSEIF ( IsInWord(unt,'/y')   .OR. IsInWord(unt,'y-1') .OR. &
             IsInWord(unt,'y^-1') .OR. IsInWord(unt,'yr')  .OR. &
             IsInWord(unt,'year') ) THEN

       IF ( IS_LEAPYEAR ) THEN
          Scal = 1.0_hp / SEC_IN_LEAPYEAR
       ELSE
          Scal = 1.0_hp / SEC_IN_REGYEAR
       ENDIF
       Flag = 1
    ENDIF

  END SUBROUTINE HCO_UNIT_GetTimeScal
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Unit_GetAreaScal
!
! !DESCRIPTION: Returns the area/volume scale factors for the given unit.
! This is the scale factor required to convert from unit 'Unit' to
! HEMCO units (i.e. per m2 or per m3).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Unit_GetAreaScal( unt, Scal, Flag )
!
! !USES:
!
    USE HCO_CharTools_Mod
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: unt   ! This unit
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),         INTENT(OUT) :: Scal  ! scale factor
    INTEGER,          INTENT(OUT) :: Flag  ! 2=per area, 3= per volume, 0 otherwise
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! HCO_UNIT_GetAreaScal begins here
    !=================================================================

    ! Init
    Scal = -999.0_hp
    Flag = 0

    ! cm2
    IF ( IsInWord(unt,'/cm2' ) .OR. IsInWord(unt,'cm-2' ) .OR. &
         IsInWord(unt,'/cm^2') .OR. IsInWord(unt,'cm^-2')       ) THEN
       Scal = 1.0_hp / 1e-4_hp
       Flag = 2

    ! km2
    ELSEIF ( IsInWord(unt,'/km2' ) .OR. IsInWord(unt,'km-2' ) .OR. &
             IsInWord(unt,'/km^2') .OR. IsInWord(unt,'km^-2')       ) THEN
       Scal = 1.0_hp / 1e6_hp
       Flag = 2

    ! m2
    ELSEIF ( IsInWord(unt,'/m2' ) .OR. IsInWord(unt,'m-2' ) .OR. &
             IsInWord(unt,'/m^2') .OR. IsInWord(unt,'m^-2')       ) THEN
       Scal = 1.0_hp
       Flag = 2

    !=================================================================
    ! Convert volume to m3
    !=================================================================

    ! cm3
    ELSEIF ( IsInWord(unt,'/cm3')  .OR. IsInWord(unt,'cm-3' ) .OR. &
             IsInWord(unt,'/cm^3') .OR. IsInWord(unt,'cm^-3')       ) THEN
       Scal = 1.0_hp / 1e-6_hp
       Flag = 3

    ! dm3
    ELSEIF ( IsInWord(unt,'/dm3')  .OR. IsInWord(unt,'dm-3' ) .OR. &
             IsInWord(unt,'/dm^3') .OR. IsInWord(unt,'dm^-3')       ) THEN
       Scal = 1.0_hp / 1e-3_hp
       Flag = 3

    ! m3
    ELSEIF ( IsInWord(unt,'/m3')  .OR. IsInWord(unt,'m-3' ) .OR. &
             IsInWord(unt,'/m^3') .OR. IsInWord(unt,'m^-3')       ) THEN
       Scal = 1.0_hp
       Flag = 3

    ! L
    ELSEIF ( IsInWord(unt,'/l') .OR. IsInWord(unt,'l-1') .OR. &
             IsInWord(unt,'l^-1')                              ) THEN
       Scal = 1.0_hp / 1e-3_hp
       Flag = 3
    ENDIF

  END SUBROUTINE HCO_Unit_GetAreaScal
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Unit_ScalCheck
!
! !DESCRIPTION: Check if the provided unit is unitless. Returns
! 0 if Unit is unitless, 1 if it's not unitless but in correct
! HEMCO emission units (i.e. kg/m2/s), 2 if it's in HEMCO concentration
! units (kg/m3), -1 otherwise.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_Unit_ScalCheck( Unit ) Result ( Flag )
!
! !USES:
!
    USE CharPak_Mod, ONLY : TRANLC
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: Unit
!
! !OUTPUT PARAMETERS:
!
    INTEGER                       :: Flag  ! 0=ok, 1=warning, 2=error
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=31) :: tmpU
    INTEGER           :: I

    !=================================================================
    ! HCO_UNIT_SCALCHECK begins here
    !=================================================================

    ! Mirror
    tmpU = Unit

    ! lower case
    CALL TRANLC( tmpU )

    ! No recognized unit (default):
    Flag = -1

    ! Check for unitless factors
    DO I = 1, NUL
       IF ( TRIM(tmpU) == TRIM(UL(I)) ) THEN
          Flag = 0
          RETURN
       ENDIF
    ENDDO

    ! Check for HEMCO emissions units
    DO I = 1, NHE
       IF ( TRIM(tmpU) == TRIM(HE(I)) ) THEN
          Flag = 1
          RETURN
       ENDIF
    ENDDO

    ! Check for HEMCO concentration units
    DO I = 1, NHC
       IF ( TRIM(tmpU) == TRIM(HC(I)) ) THEN
          Flag = 2
          RETURN
       ENDIF
    ENDDO

  END FUNCTION HCO_Unit_ScalCheck
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_IsUnitless
!
! !DESCRIPTION: Returns TRUE if the passed string corresponds to the HEMCO
! unitless data unit string.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_IsUnitless( Unt ) Result ( Bool )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)   :: Unt
!
! !OUTPUT PARAMETERS:
!
    LOGICAL                        :: Bool
!
! !REVISION HISTORY:
!  24 Jul 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    Bool = ( TRIM(Unt) == TRIM(UL(1)) )

  END FUNCTION HCO_IsUnitless
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_IsIndexData
!
! !DESCRIPTION: Returns TRUE if the passed string corresponds to the HEMCO
! index data unit string.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_IsIndexData( Unt ) Result ( Bool )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)   :: Unt
!
! !OUTPUT PARAMETERS:
!
    LOGICAL                        :: Bool
!
! !REVISION HISTORY:
!  24 Jul 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    Bool = ( TRIM(Unt) == TRIM(UL(2)) )

  END FUNCTION HCO_IsIndexData
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_UnitTolerance
!
! !DESCRIPTION: Returns the HEMCO unit tolerance as defined in the HEMCO
! configuration file (under settings). Returns a default value of 0 (no
! tolerance) if no value is set in the configuration file.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_UnitTolerance( HcoConfig ) Result ( UnitTolerance )
!
! !USES:
!
    USE HCO_TYPES_MOD,    ONLY : ConfigObj
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt, CoreNr
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj), POINTER       :: HcoConfig
!
! !OUTPUT PARAMETERS:
!
    INTEGER                        :: UnitTolerance
!
! !REVISION HISTORY:
!  24 Jul 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, SAVE       :: Tolerance = -1
    LOGICAL             :: FOUND
    INTEGER             :: RC
    CHARACTER(LEN=255)  :: MSG

    !=================================================================
    ! HCO_UnitTolerance begins here
    !=================================================================

    ! On first call, try to get unit tolerance value from core settings.
    ! Use value of zero if not specified in the configuration file.
    IF ( Tolerance < 0 ) THEN
       CALL GetExtOpt ( HcoConfig, CoreNr, 'Unit tolerance', &
                        OptValInt=Tolerance, FOUND=FOUND, RC=RC )
       IF ( .NOT. FOUND ) Tolerance = 0
    ENDIF

    ! Return
    UnitTolerance = Tolerance

  END FUNCTION HCO_UnitTolerance
!EOC
END MODULE HCO_Unit_Mod
