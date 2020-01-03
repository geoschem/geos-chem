#if defined( TOMAS )
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_tomas_jeagle_mod.F90
!
! !DESCRIPTION: Module HCOX\_TOMAS\_JEAGLE\_Mod contains routines to
!  calculate sea salt aerosol emissions for the TOMAS aerosol microphysics
!  package. JKODROS - This is an update of hcox\_tomas\_seasalt\_mod.F90 to
!  use Jeagle emissions. Should bring TOMAS emissions in line with bulk sea
!  salt.
!\\
!\\
!  This is a HEMCO extension module that uses many of the HEMCO core
!  utilities.
!\\
!\\
!  References:
!  \begin{itemize}
!  \item Clarke, A.D., Owens, S., Zhou, J. \emph{An ultrafine sea-salt flux
!        from breaking waves: Implications for CCN in the remote marine
!        atmosphere}, \underline{J. Geophys. Res.}, 2006.
!  \end{itemize}
!
! !INTERFACE:
!
MODULE HCOX_TOMAS_Jeagle_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE HCO_State_Mod,  ONLY : HCO_State
  USE HCOX_State_Mod, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_TOMAS_Jeagle_Init
  PUBLIC :: HCOX_TOMAS_Jeagle_Run
  PUBLIC :: HCOX_TOMAS_Jeagle_Final
!
! !REVISION HISTORY:
!  01 Oct 2014 - R. Yantosca - Initial version, based on TOMAS code
!  20 May 2015 - J. Kodros   - Added fixes to integrate TOMAS with HEMCO
!  02 JUL 2015 - J. Kodros   - Updating to use scale factors from Jeagle
!                              et al. (2011)
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE TYPES:
!
  TYPE :: MyInst

   ! Scalars
   INTEGER               :: Instance
   INTEGER               :: ExtNr                  ! HEMCO extension #
   REAL(dp)              :: TOMAS_COEF             ! Seasalt emiss coeff.

   ! Arrays
   INTEGER,  ALLOCATABLE :: HcoIDs    (:      )    ! HEMCO species ID's
   REAL(dp), POINTER     :: TOMAS_DBIN(:      )    ! TOMAS bin width
   REAL(dp), POINTER     :: DRFAC     (:      )    ! TOMAS area?
   REAL(dp), POINTER     :: TC1       (:,:,:,:)    ! Aerosol mass
   REAL(dp), POINTER     :: TC2       (:,:,:,:)    ! Aerosol number

   TYPE(MyInst), POINTER :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to instances
  TYPE(MyInst), POINTER  :: AllInst => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_TOMAS_Jeagle_Run
!
! !DESCRIPTION: Subroutine HCOX\_TOMAS\_Jeagle\_Run emits sea-salt into the
!  TOMAS sectional sea-salt mass and aerosol number arrays.  Sea-salt emission
!  parameterization of Jeagle et al. (2011).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_TOMAS_Jeagle_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_GeoTools_Mod, ONLY : HCO_LandType
    USE HCO_FluxArr_mod,  ONLY : HCO_EmisAdd
    USE HCO_State_Mod,    ONLY : HCO_GetHcoID
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState    ! Extension Options object
    TYPE(HCO_State),  POINTER       :: HcoState    ! HEMCO state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC          ! Success or failure?
!
! !REMARKS:
!
!
! !REVISION HISTORY:
!  01 Oct 2014 - R. Yantosca - Initial version, based on TOMAS SRCSALT30 code
!  20 May 2015 - J. Kodros   - Add seasalt number & mass to HEMCO state
!  20 May 2015 - R. Yantosca - Pass am_I_Root to HCO_EMISADD routine
!  22 May 2015 - R. Yantosca - Extend up to 40 size bins
!  10 Jul 2015 - R. Yantosca - Fixed minor issues in the ProTeX headers
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER           :: I,      J,    L,      K, HcoID
    REAL(sp)          :: FOCEAN, W10M, DTEMIS
    REAL(dp)          :: F100,   W, A_M2, FEMIS, NUMBER, MASS, NUMBER_TOT
    REAL(dp)          :: rwet, dfo, B, A, SST, SCALE
    CHARACTER(LEN=255):: MSG

    REAL*8, PARAMETER :: BETHA=2.22d0   !wet diameter (80% Rel Hum) to dry diam

    ! Strings
    CHARACTER(LEN=31) :: SpcName

    ! Pointers
    TYPE(MyInst), POINTER :: Inst
    REAL(dp), POINTER :: ptr3D(:,:,:)

    ! For debugging
    !INTEGER            :: ii=50, jj=10

    ! Error handling
    LOGICAL                :: ERR

    !=================================================================
    ! SRCSALT30 begins here!
    !=================================================================

    ! Return if extension disabled
    IF ( ExtState%TOMAS_Jeagle <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER ( HcoState%Config%Err, 'HCOX_TOMAS_Jeagle_Run (hcox_TOMAS_Jeagle_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get instance
    Inst   => NULL()
    CALL InstGet ( ExtState%TOMAS_Jeagle, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find TOMAS_Jeagle instance Nr. ', ExtState%TOMAS_Jeagle
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
    ENDIF

    !INIT VALUES
    Inst%TC1 = 0.0_hp
    Inst%TC2 = 0.0_hp

    ! Depending on the grid resolution. 4x5 (default) doesn't need
    ! adjusting coeff

    !### Debug
    !print*, 'JACK IN HCOX TOMAS Jeagle'

    ! Init
    ptr3D => NULL()

    ! Emission timestep [s]
    DTEMIS = HcoState%TS_EMIS

    ! Loop over grid cells
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! Grid box surface area [m2]
       A_M2  = HcoState%Grid%AREA_M2%Val(I,J)

       ! Get the fraction of the box that is over water
       IF ( HCO_LandType( ExtState%WLI%Arr%Val(I,J),              &
                          ExtState%ALBD%Arr%Val(I,J) ) == 0 ) THEN
          FOCEAN = 1d0 - ExtState%FRCLND%Arr%Val(I,J)
       ELSE
          FOCEAN = 0.d0
       ENDIF

       ! Skip boxes that are not at least 50% water
       IF ( FOCEAN > 0.5d0 ) THEN

          ! Wind speed at 10 m altitude [m/s]
          W10M = SQRT( ExtState%U10M%Arr%Val(I,J)**2  &
               +       ExtState%V10M%Arr%Val(I,J)**2 )

      ! Sea surface temperature in Celcius
      SST = ExtState%TSKIN%Arr%Val(I,J) - 273.15d0

      ! Limit SST to 0-30C range
          SST = MAX( SST , 0d0 )  ! limit to  0C
          SST = MIN( SST , 30d0 ) ! limit to 30C

          ! Empirical SST scaling factor (jaegle 5/11/11)
          SCALE = 0.329d0 + 0.0904d0*SST -  &
                  0.00717d0*SST**2d0 + 0.000207d0*SST**3d0

          !---------------------------------------------------------------
          ! Partition TOMAS_Jeagle emissions w/in the boundary layer
          !---------------------------------------------------------------
          DO K = 1, HcoState%MicroPhys%nBins
             rwet=Inst%TOMAS_DBIN(k)*1.0E6*BETHA/2. ! convert from dry diameter [m] to wet (80% RH) radius [um]
         ! jkodros - testing out BETHA 7/29/15
             if (rwet > 0.d0) then
                  A=4.7*(1.+30.*rwet)**(-0.017*rwet**(-1.44))
                  B=(0.433-log10(rwet))/0.433

                  dfo=1.373*W10M**3.41*rwet**(-1.*A)  & !m-2 um-1 s-1
                    *(1.+0.057*rwet**3.45)*10.**(1.607*exp(-1.*B**2))

             else

        dfo=0.d0
         endif

             dfo=dfo*Inst%DRFAC(k)*BETHA  !hemco units???? jkodros
         dfo=dfo*focean*SCALE

             ! Loop thru the boundary layer
             DO L = 1, HcoState%Nz

                ! Fraction of the PBL spanned by box (I,J,L) [unitless]
                FEMIS = ExtState%FRAC_OF_PBL%Arr%Val(I,J,L)

                ! Only handle grid boxes w/in the boundary layer
                IF ( FEMIS > 0d0 ) THEN

                   ! Number
                   NUMBER = dfo * FEMIS

                   ! Mass
                   MASS   = NUMBER                                      &
                          * SQRT( HcoState%MicroPhys%BinBound(K  ) *    &
                                  HcoState%MicroPhys%BinBound(K+1)   )

                   ! Store number & mass
                   Inst%TC1(I,J,L,K) = NUMBER
                   Inst%TC2(I,J,L,K) = MASS

                ENDIF
             ENDDO
          ENDDO
       ELSE
          Inst%TC1(I,J,:,:) = 0d0
          Inst%TC2(I,J,:,:) = 0d0
       ENDIF
    ENDDO
    ENDDO

    !### Debug
    !print*, 'JACK SEASALT EMISSIONS AT 50, 10,7: ', TC2(ii,jj,1,7)
    !print*, 'BINS: ', HcoState%MicroPhys%nBins

    ! Loop over # of microphysics bins
    DO K = 1, HcoState%MicroPhys%nBins

       ! Add mass to the HEMCO data structure (jkodros)
       CALL HCO_EmisAdd( HcoState, Inst%TC2(:,:,:,K), Inst%HcoIDs(K), RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: FLUXSALT', RC )
          RETURN
       ENDIF

       ! Get the proper species name
       IF ( K<10  ) THEN
         WRITE(SpcName,'(A2,I1)') 'NK', K
       ELSE
         WRITE(SpcName,'(A2,I2)') 'NK', K
       ENDIF

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( TRIM(SpcName), HcoState )

       !### Debug
       !print*, 'JACK SEASALT EMISSIONS AT 50, 10,: ', TC1(ii,jj,1,k)
       !print*, 'JACK HCO ID: ', HcoID

       ! Add number to the HEMCO data structure
       CALL HCO_EmisAdd( HcoState, Inst%TC1(:,:,:,K), HcoID, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: FLUXSALT', RC )
          RETURN
       ENDIF
    ENDDO

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Nullify pointers
    Inst    => NULL()

    ! Leave w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE HCOX_TOMAS_Jeagle_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_TOMAS_Jeagle_Init
!
! !DESCRIPTION: Subroutine HcoX\_TOMAS\_Jeagle\_Init initializes all
!  extension variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_TOMAS_Jeagle_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_State_Mod,   ONLY : HCO_GetHcoID
    USE HCO_STATE_MOD,   ONLY : HCO_GetExtHcoID
    USE HCO_ExtList_Mod, ONLY : GetExtNr
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState    ! HEMCO state object
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName     ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState    ! Extension options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  15 Dec 2013 - C. Keller   - Initial version
!  10 Jul 2015 - R. Yantosca - Fixed minor issues in ProTeX header
!  24 Aug 2017 - M. Sulprizio- Remove support for GRID1x1
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                        :: N, R, AS, ExtNr
    REAL*8                         :: A, B, R0, R1
    REAL*8                         :: CONST_N
    INTEGER                        :: nSpc, minLen
    CHARACTER(LEN=255)             :: MSG

    ! Arrays
!    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)

    ! Pointers
    TYPE(MyInst), POINTER          :: Inst

    !=================================================================
    ! HCOX_TOMAS_Jeagle_Init begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_TOMAS_Jeagle_Init (hcox_tomas_jeagle_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create Instance
    Inst => NULL()
    CALL InstCreate ( ExtNr, ExtState%TOMAS_Jeagle, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot create TOMAS_Jeagle instance', RC )
       RETURN
    ENDIF
    ! Also fill Inst%ExtNr
    Inst%ExtNr = ExtNr

    ! ----------------------------------------------------------------------
    ! Get species IDs and settings
    ! ----------------------------------------------------------------------

    ! Get HEMCO species IDs
    CALL HCO_GetExtHcoID( HcoState, Inst%ExtNr, Inst%HcoIDs, SpcNames, &
                          nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( nSpc < HcoState%MicroPhys%nBins ) THEN
       MSG = 'Not enough sea salt emission species set'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! Allocate TOMAS_DBIN
    ALLOCATE ( Inst%TOMAS_DBIN( HcoState%MicroPhys%nBins ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot allocate TOMAS_DBIN array (hcox_tomas_jeagle_mod.F90)'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! Allocate TOMAS_A
    ALLOCATE ( Inst%DRFAC( HcoState%MicroPhys%nBins ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot allocate DRFAC array (hcox_tomas_jeagle_mod.F90)'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! JKODROS - ALLOCATE TC1 and TC2
    ALLOCATE ( Inst%TC1( HcoState%NX, HcoState%NY,&
               HcoState%NZ, HcoState%MicroPhys%nBins ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot allocate TC1 array (hcox_tomas_jeagle_mod.F90)'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ELSE
    Inst%TC1 = 0d0
    ENDIF

    ! JKODROS - ALLOCATE TC1 and TC2
    ALLOCATE ( Inst%TC2( HcoState%NX, HcoState%NY,&
               HcoState%NZ, HcoState%MicroPhys%nBins ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot allocate TC2 array (hcox_tomas_jeagle_mod.F90)'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ELSE
    Inst%TC2 = 0d0
    ENDIF

! ----- IMPORTANT BINS ONLY CORRECTLY SET UP FOR TOMAS 15 PLEASE ADJUST OTHERS -jkodros (7/21/15)
! ----  6/24/16 - JKodros - I have updated the DRFAC. They should (hopefully) be
! ----  correct now. DRFAC is the bin width (radius not diameter) for DRY SS
#if defined( TOMAS12 )

    !-----------------------------------------------------------------------
    ! TOMAS simulation w/ 12 size-resolved bins
    !-----------------------------------------------------------------------

    Inst%TOMAS_DBIN = (/                                                &
          9.68859d-09,   1.53797d-08,    2.44137d-08,    3.87544d-08,   &
          6.15187d-08,   9.76549d-08,    1.55017d-07,    2.46075d-07,   &
          3.90620d-07,   6.20070d-07,    9.84300d-07,    3.12500d-06  /)

    Inst%DRFAC     = (/                                                 &
          2.84132d-03,   4.51031d-03,    7.15968d-03,    1.13653d-02,   &
          1.80413d-02,   2.86387d-02,    4.54612d-02,    7.21651d-02,   &
          1.14555d-01,   1.81845d-01,    1.06874d+00,    3.39304d+00 /)
#elif defined( TOMAS15 )

    !-----------------------------------------------------------------------
    ! TOMAS simulation w/ 15 size-resolved bins
    !-----------------------------------------------------------------------

    Inst%TOMAS_DBIN = (/         0d0,            0d0,             0d0,   &
          1.22069d-08,   1.93772d-08,    3.07594d-08,     4.88274d-08,   &
          7.75087d-08,   1.23037d-07,    1.95310d-07,     3.10035d-07,   &
          4.92150d-07,   7.81239d-07,    1.74054d-06,     5.52588d-06 /)

    Inst%DRFAC      = (/         0d0,            0d0,             0d0,   &
          2.84132d-03,   4.51031d-03,    7.15968d-03,     1.13653d-02,   &
          1.80413d-02,   2.86387d-02,    4.54612d-02,     7.21651d-02,   &
          1.14555d-01,   1.81845d-01,    1.06874d+00,     3.39304d+00 /)

#elif defined( TOMAS40 )

    !-----------------------------------------------------------------------
    ! TOMAS simulation w/ 40 size-resolved bins
    !-----------------------------------------------------------------------

    Inst%TOMAS_DBIN = (/                                                 &
       0.0d0      , 0.0d0      , 0.0d0     ,  0.0d0      , 0.0d0      ,  &
       0.0d0      , 0.0d0      , 0.0d0     ,  0.0d0      , 0.0d0      ,  &
       9.68859d-09, 1.22069d-08, 1.53797d-08, 1.93772d-08, 2.44137d-08,  &
       3.07594d-08, 3.87544d-08, 4.88274d-08, 6.15187d-08, 7.75087d-08,  &
       9.76549d-08, 1.23037d-07, 1.55017d-07, 1.95310d-07, 2.46075d-07,  &
       3.10035d-07, 3.90620d-07, 4.92150d-07, 6.20070d-07, 7.81239d-07,  &
       9.84300d-07, 1.24014d-06, 1.56248d-06, 1.96860d-06, 2.48028d-06,  &
       3.12496d-06, 3.93720d-06, 4.96056d-06, 6.24991d-06, 7.87440d-06 /)

    Inst%DRFAC     = (/                                                  &
       0.0d0      , 0.0d0      , 0.0d0      , 0.0d0      , 0.0d0      ,  &
       0.0d0      , 0.0d0      , 0.0d0      , 0.0d0      , 0.0d0      ,  &
       1.24737d-03, 1.57158d-03, 1.98007d-03, 2.49473d-03, 3.14317d-03,  &
       3.96014d-03, 4.98947d-03, 6.28633d-03, 7.92028d-03, 9.97893d-03,  &
       1.25727d-02, 1.58406d-02, 1.99579d-02, 2.51453d-02, 3.16811d-02,  &
       3.99157d-02, 5.02906d-02, 6.33623d-02, 7.98314d-02, 1.00581d-01,  &
       1.26725d-01, 1.59663d-01, 2.01163d-01, 2.53449d-01, 3.19326d-01,  &
       4.02325d-01, 5.06898d-01, 6.38652d-01, 8.04651d-01, 1.01380d+00 /)
#else

    !-----------------------------------------------------------------------
    ! TOMAS simulation w/ 30 size-resolved bins (default)
    !-----------------------------------------------------------------------

    Inst%TOMAS_DBIN = (/                                                 &
       9.68859d-09, 1.22069d-08, 1.53797d-08, 1.93772d-08, 2.44137d-08,  &
       3.07594d-08, 3.87544d-08, 4.88274d-08, 6.15187d-08, 7.75087d-08,  &
       9.76549d-08, 1.23037d-07, 1.55017d-07, 1.95310d-07, 2.46075d-07,  &
       3.10035d-07, 3.90620d-07, 4.92150d-07, 6.20070d-07, 7.81239d-07,  &
       9.84300d-07, 1.24014d-06, 1.56248d-06, 1.96860d-06, 2.48028d-06,  &
       3.12496d-06, 3.93720d-06, 4.96056d-06, 6.24991d-06, 7.87440d-06 /)

    Inst%DRFAC     = (/                                                  &
       1.24737d-03, 1.57158d-03, 1.98007d-03, 2.49473d-03, 3.14317d-03,  &
       3.96014d-03, 4.98947d-03, 6.28633d-03, 7.92028d-03, 9.97893d-03,  &
       1.25727d-02, 1.58406d-02, 1.99579d-02, 2.51453d-02, 3.16811d-02,  &
       3.99157d-02, 5.02906d-02, 6.33623d-02, 7.98314d-02, 1.00581d-01,  &
       1.26725d-01, 1.59663d-01, 2.01163d-01, 2.53449d-01, 3.19326d-01,  &
       4.02325d-01, 5.06898d-01, 6.38652d-01, 8.04651d-01, 1.01380d+00 /)
#endif

    !=======================================================================
    ! Allocate quantities depending on horizontal resolution
    !=======================================================================
    IF ( TRIM( HcoState%Config%GridRes) == '4.0x5.0' ) THEN

       !-----------------------------------------------------------------------
       ! TOMAS simulations at 4 x 5 global resolution
       !-----------------------------------------------------------------------
       Inst%TOMAS_COEF = 1.d0

    ELSE IF ( TRIM( HcoState%Config%GridRes) == '2.0x2.5' ) THEN

       !-----------------------------------------------------------------------
       ! TOMAS simulations at 2 x 2.5 global resolution
       !-----------------------------------------------------------------------
       Inst%TOMAS_COEF = 1.d0

    ELSE

       MSG = 'Adjust TOMAS_Jeagle emiss coeff (TOMAS_COEF) for your model res: SRCSALT30: hcox_TOMAS_jeagle_mod.F90'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )

    ENDIF

    !=======================================================================
    ! Activate this module and the fields of ExtState that it uses
    !=======================================================================

    ! Activate met fields
    ExtState%WLI%DoUse         = .TRUE.
    ExtState%ALBD%DoUse        = .TRUE.
    ExtState%TSKIN%DoUse       = .TRUE.
    ExtState%U10M%DoUse        = .TRUE.
    ExtState%V10M%DoUse        = .TRUE.
    ExtState%FRAC_OF_PBL%DoUse = .TRUE.
    ExtState%FRCLND%DoUse      = .TRUE.

    !=======================================================================
    ! Leave w/ success
    !=======================================================================
    IF ( ALLOCATED( SpcNames ) ) DEALLOCATE( SpcNames )

    ! Nullify pointers
    Inst    => NULL()

    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE HCOX_TOMAS_Jeagle_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_TOMAS_Jeagle_Final
!
! !DESCRIPTION: Subroutine HcoX\_TOMAS\_Jeagle\_Final deallocates
!  all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_TOMAS_Jeagle_Final( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !REVISION HISTORY:
!  15 Dec 2013 - C. Keller   - Initial version
!  20 May 2015 - J. Kodros   - Deallocate HcoIDs, TC1, TC2 arrays
!EOP
!------------------------------------------------------------------------------
!BOC
!
    CALL InstRemove ( ExtState%TOMAS_Jeagle )

  END SUBROUTINE HCOX_TOMAS_Jeagle_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a poiner to the desired instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstGet ( Instance, Inst, RC, PrevInst )
!
! !INPUT PARAMETERS:
!
    INTEGER                             :: Instance
    TYPE(MyInst),     POINTER           :: Inst
    INTEGER                             :: RC
    TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst),     POINTER    :: PrvInst

    !=================================================================
    ! InstGet begins here!
    !=================================================================

    ! Get instance. Also archive previous instance.
    PrvInst => NULL()
    Inst    => AllInst
    DO WHILE ( ASSOCIATED(Inst) )
       IF ( Inst%Instance == Instance ) EXIT
       PrvInst => Inst
       Inst    => Inst%NextInst
    END DO
    IF ( .NOT. ASSOCIATED( Inst ) ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Pass output arguments
    IF ( PRESENT(PrevInst) ) PrevInst => PrvInst

    ! Cleanup & Return
    PrvInst => NULL()
    RC = HCO_SUCCESS

  END SUBROUTINE InstGet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate
!
! !DESCRIPTION: Subroutine InstCreate creates a new instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstCreate ( ExtNr, Instance, Inst, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,       INTENT(IN)       :: ExtNr
!
! !OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(  OUT)    :: Instance
    TYPE(MyInst),  POINTER          :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(INOUT)    :: RC
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst), POINTER          :: TmpInst
    INTEGER                        :: nnInst

    !=================================================================
    ! InstCreate begins here!
    !=================================================================

    ! ----------------------------------------------------------------
    ! Generic instance initialization
    ! ----------------------------------------------------------------

    ! Initialize
    Inst => NULL()

    ! Get number of already existing instances
    TmpInst => AllInst
    nnInst = 0
    DO WHILE ( ASSOCIATED(TmpInst) )
       nnInst  =  nnInst + 1
       TmpInst => TmpInst%NextInst
    END DO

    ! Create new instance
    ALLOCATE(Inst)
    Inst%Instance = nnInst + 1
    Inst%ExtNr    = ExtNr

    ! Attach to instance list
    Inst%NextInst => AllInst
    AllInst       => Inst

    ! Update output instance
    Instance = Inst%Instance

    ! ----------------------------------------------------------------
    ! Type specific initialization statements follow below
    ! ----------------------------------------------------------------

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!BOP
!
! !IROUTINE: InstRemove
!
! !DESCRIPTION: Subroutine InstRemove creates a new instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstRemove ( Instance )
!
! !INPUT PARAMETERS:
!
    INTEGER                         :: Instance
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                     :: RC
    TYPE(MyInst), POINTER       :: PrevInst
    TYPE(MyInst), POINTER       :: Inst

    !=================================================================
    ! InstRemove begins here!
    !=================================================================

    ! Init
    PrevInst => NULL()
    Inst     => NULL()

    ! Get instance. Also archive previous instance.
    CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

    ! Instance-specific deallocation
    IF ( ASSOCIATED(Inst) ) THEN

       ! Pop off instance from list
       IF ( ASSOCIATED(PrevInst) ) THEN

          IF ( ASSOCIATED( Inst%TOMAS_DBIN ) ) DEALLOCATE( Inst%TOMAS_DBIN )
          IF ( ASSOCIATED( Inst%DRFAC      ) ) DEALLOCATE( Inst%DRFAC      )
          IF ( ASSOCIATED( Inst%TC1        ) ) DEALLOCATE( Inst%TC1        )
          IF ( ASSOCIATED( Inst%TC2        ) ) DEALLOCATE( Inst%TC2        )
          IF ( ALLOCATED ( Inst%HcoIDs     ) ) DEALLOCATE( Inst%HcoIDs     )

          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL()
    ENDIF

   END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_TOMAS_Jeagle_Mod
#endif
