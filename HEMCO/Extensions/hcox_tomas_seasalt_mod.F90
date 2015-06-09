#if defined( TOMAS )
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_tomas_seasalt_mod.F90
!
! !DESCRIPTION: Module HCOX\_TOMAS\_SeaSalt\_Mod contains routines to 
!  calculate sea salt aerosol emissions for the TOMAS aerosol microphysics
!  package. 
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
!
! !INTERFACE: 
!
MODULE HCOX_TOMAS_SeaSalt_Mod
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
  PUBLIC :: HCOX_TOMAS_SeaSalt_Init
  PUBLIC :: HCOX_TOMAS_SeaSalt_Run
  PUBLIC :: HCOX_TOMAS_SeaSalt_Final
!
! !REVISION HISTORY:
!  01 Oct 2014 - R. Yantosca - Initial version, based on TOMAS code
!  20 May 2015 - J. Kodros   - Added fixes to integrate TOMAS with HEMCO
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE TYPES:
!
  ! Scalars
  INTEGER               :: ExtNr                  ! HEMCO extension #   
  REAL(dp)              :: TOMAS_COEF             ! Seasalt emiss coeff.

  ! Arrays
  INTEGER,  ALLOCATABLE :: HcoIDs    (:      )    ! HEMCO species ID's
  REAL(dp), ALLOCATABLE :: TOMAS_DBIN(:      )    ! TOMAS bin width
  REAL(dp), ALLOCATABLE :: TOMAS_A   (:      )    ! TOMAS area? 
  REAL(dp), ALLOCATABLE :: TC1       (:,:,:,:)    ! Aerosol mass
  REAL(dp), ALLOCATABLE :: TC2       (:,:,:,:)    ! Aerosol number

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_TOMAS_SeaSalt_Run
!
! !DESCRIPTION: Subroutine SRCSALT_TOMAS emits sea-salt into the TOMAS
!  sectional sea-salt mass and aerosol number arrays.  Sea-salt emission 
!  parameterization of Clarke et al. [2006].  Formerly named SRCSALT30.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_TOMAS_SeaSalt_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_GeoTools_Mod, ONLY : HCO_LandType
    USE HCO_FluxArr_mod,  ONLY : HCO_EmisAdd
    USE HCO_State_Mod,    ONLY : HCO_GetHcoID
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root   ! root CPU?
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

    ! Strings
    CHARACTER(LEN=31) :: SpcName
       
    ! Pointers
    REAL(dp), POINTER :: ptr3D(:,:,:) => NULL()

    ! For debugging
    !INTEGER            :: ii=50, jj=10

    !=================================================================
    ! SRCSALT30 begins here!
    !=================================================================

    ! Depending on the grid resolution. 4x5 (default) doesn't need 
    ! adjusting coeff

    !### Debug
    !print*, 'JACK IN HCOX TOMAS SEASALT'

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
                
          ! in ocean area - calc wind speed/eqm conc
          ! calculate the fraction of whitecap coverage
          W = 3.84E-6 * W10M ** (3.41)

          !---------------------------------------------------------------
          ! Partition TOMAS_SeaSalt emissions w/in the boundary layer
          !---------------------------------------------------------------
          DO K = 1, HcoState%MicroPhys%nBins
               
             ! Sea salt number
             F100   = TOMAS_A(K)

             ! JKODROS get number as flux
             NUMBER_TOT = F100 * W * FOCEAN * TOMAS_COEF

             ! Loop thru the boundary layer
             DO L = 1, HcoState%Nz

                ! Fraction of the PBL spanned by box (I,J,L) [unitless]
                FEMIS = ExtState%FRAC_OF_PBL%Arr%Val(I,J,L)

                ! Only handle grid boxes w/in the boundary layer
                IF ( FEMIS > 0d0 ) THEN

                   ! Number
                   NUMBER = NUMBER_TOT * FEMIS

                   ! Mass
                   MASS   = NUMBER                                      &
                          * SQRT( HcoState%MicroPhys%BinBound(K  ) *    &
                                  HcoState%MicroPhys%BinBound(K+1)   )
                   
                   ! Store number & mass
                   TC1(I,J,L,K) = NUMBER
                   TC2(I,J,L,K) = MASS

                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    ENDDO

    !### Debug
    !print*, 'JACK SEASALT EMISSIONS AT 50, 10,7: ', TC2(ii,jj,1,7)
    !print*, 'BINS: ', HcoState%MicroPhys%nBins

    ! Loop over # of microphysics bins
    DO K = 1, HcoState%MicroPhys%nBins

       ! Add mass to the HEMCO data structure (jkodros)
       CALL HCO_EmisAdd( am_I_Root, HcoState, TC2(:,:,:,K), HcoIDs(K), RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXSALT', RC )
          RETURN
       ENDIF

       ! Get the proper species name
       !%%% This is a kludge, do better later %%%
       IF ( K==1  ) SpcName = 'NK1'
       IF ( K==2  ) SpcName = 'NK2'
       IF ( K==3  ) SpcName = 'NK3'
       IF ( K==4  ) SpcName = 'NK4'
       IF ( K==5  ) SpcName = 'NK5'
       IF ( K==6  ) SpcName = 'NK6'
       IF ( K==7  ) SpcName = 'NK7'
       IF ( K==8  ) SpcName = 'NK8'
       IF ( K==9  ) SpcName = 'NK9'
       IF ( K==10 ) SpcName = 'NK10'
       IF ( K==11 ) SpcName = 'NK11'
       IF ( K==12 ) SpcName = 'NK12'
       IF ( K==13 ) SpcName = 'NK13'
       IF ( K==14 ) SpcName = 'NK14'
       IF ( K==15 ) SpcName = 'NK15'
       IF ( K==16 ) SpcName = 'NK16'
       IF ( K==17 ) SpcName = 'NK17'
       IF ( K==18 ) SpcName = 'NK18'
       IF ( K==19 ) SpcName = 'NK19'
       IF ( K==20 ) SpcName = 'NK20'
       IF ( K==21 ) SpcName = 'NK21'
       IF ( K==22 ) SpcName = 'NK22'
       IF ( K==23 ) SpcName = 'NK23'
       IF ( K==24 ) SpcName = 'NK24'
       IF ( K==25 ) SpcName = 'NK25'
       IF ( K==26 ) SpcName = 'NK26'
       IF ( K==27 ) SpcName = 'NK27'
       IF ( K==28 ) SpcName = 'NK28'
       IF ( K==29 ) SpcName = 'NK29'
       IF ( K==30 ) SpcName = 'NK30'
       IF ( K==31 ) SpcName = 'NK31'
       IF ( K==32 ) SpcName = 'NK32'
       IF ( K==33 ) SpcName = 'NK33'
       IF ( K==34 ) SpcName = 'NK34'
       IF ( K==35 ) SpcName = 'NK35'
       IF ( K==36 ) SpcName = 'NK36'
       IF ( K==37 ) SpcName = 'NK37'
       IF ( K==38 ) SpcName = 'NK38'
       IF ( K==39 ) SpcName = 'NK39'
       IF ( K==40 ) SpcName = 'NK40'

       ! HEMCO species ID
       HcoID = HCO_GetHcoID( TRIM(SpcName), HcoState ) 
        
       !### Debug
       !print*, 'JACK SEASALT EMISSIONS AT 50, 10,: ', TC1(ii,jj,1,k)
       !print*, 'JACK HCO ID: ', HcoID
    
       ! Add number to the HEMCO data structure
       CALL HCO_EmisAdd( am_I_Root, HcoState, TC1(:,:,:,K), HcoID, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXSALT', RC )
          RETURN
       ENDIF
    ENDDO
                
  END SUBROUTINE HCOX_TOMAS_SeaSalt_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_TOMAS_SeaSalt_Init
!
! !DESCRIPTION: Subroutine HcoX\_TOMAS_SeaSalt\_Init initializes all
!  extension variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_TOMAS_SeaSalt_Init( am_I_Root, HcoState,     &
                                      ExtName,   ExtState, RC ) 
!
! !USES:
!
    USE HCO_State_Mod,   ONLY : HCO_GetHcoID
    USE HCO_STATE_MOD,   ONLY : HCO_GetExtHcoID
    USE HCO_ExtList_Mod, ONLY : GetExtNr
    USE HCO_ExtList_Mod, ONLY : GetExtOpt
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root   ! root CPU?
    TYPE(HCO_State),  POINTER        :: HcoState    ! HEMCO state object 
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName     ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState    ! Extension options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  15 Dec 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                        :: N, R, AS
    REAL*8                         :: A, B, R0, R1
    REAL*8                         :: CONST_N
    INTEGER                        :: nSpc, minLen
    CHARACTER(LEN=255)             :: MSG

    ! Arrays
!    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)

    !=================================================================
    ! HCOX_TOMAS_SeaSalt_Init begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN
 
    ! Enter 
    CALL HCO_ENTER( 'HCOX_TOMAS_SeaSalt_Init (hcox_tomas_seasalt_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ---------------------------------------------------------------------- 
    ! Get species IDs and settings 
    ! ---------------------------------------------------------------------- 
  
    ! Get HEMCO species IDs
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( nSpc < HcoState%MicroPhys%nBins ) THEN
       MSG = 'Not enough sea salt emission species set' 
       CALL HCO_ERROR ( MSG, RC ) 
       RETURN
    ENDIF

    ! Allocate TOMAS_DBIN
    ALLOCATE ( TOMAS_DBIN( HcoState%MicroPhys%nBins ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot allocate TOMAS_DBIN array (hcox_tomas_seasalt_mod.F90)'
       CALL HCO_ERROR( MSG, RC )      
       RETURN
    ENDIF

    ! Allocate TOMAS_A
    ALLOCATE ( TOMAS_A( HcoState%MicroPhys%nBins ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot allocate TOMAS_A array (hcox_tomas_seasalt_mod.F90)'
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF

   ! JKODROS - ALLOCATE TC1 and TC2
    ALLOCATE ( TC1( HcoState%NX, HcoState%NY,& 
           HcoState%NZ, HcoState%MicroPhys%nBins ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot allocate TC1 array (hcox_tomas_seasalt_mod.F90)'
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ELSE
 TC1 = 0d0
    ENDIF

   ! JKODROS - ALLOCATE TC1 and TC2
    ALLOCATE ( TC2( HcoState%NX, HcoState%NY,& 
           HcoState%NZ, HcoState%MicroPhys%nBins ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot allocate TC2 array (hcox_tomas_seasalt_mod.F90)'
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ELSE
 TC2 = 0d0
    ENDIF

#if defined( TOMAS12 ) 

    !-----------------------------------------------------------------------
    ! TOMAS simulation w/ 12 size-resolved bins
    !-----------------------------------------------------------------------

    TOMAS_DBIN = (/                                                     &
          9.68859d-09,   1.53797d-08,    2.44137d-08,    3.87544d-08,   &
          6.15187d-08,   9.76549d-08,    1.55017d-07,    2.46075d-07,   &
          3.90620d-07,   6.20070d-07,    9.84300d-07,    3.12500d-06  /)

    TOMAS_A   = (/                                                      &
        4607513.229d0, 9309031.200d0, 12961629.010d0, 13602132.943d0,   & 
       11441451.509d0, 9387934.311d0,  8559624.313d0,  7165322.549d0,   &
        4648135.263d0, 2447035.933d0,  3885009.997d0,  1006980.679d0  /)

#elif defined( TOMAS15 )

    !-----------------------------------------------------------------------
    ! TOMAS simulation w/ 15 size-resolved bins
    !-----------------------------------------------------------------------

    TOMAS_DBIN = (/              0d0,            0d0,            0d0,  &
          9.68859d-09,   1.53797d-08,    2.44137d-08,    3.87544d-08,  &
          6.15187d-08,   9.76549d-08,    1.55017d-07,    2.46075d-07,  &
          3.90620d-07,   6.20070d-07,    9.84300d-07,    3.12500d-06 /)

    TOMAS_A    = (/              0d0,            0d0,            0d0,  & 
        4607513.229d0, 9309031.200d0, 12961629.010d0, 13602132.943d0,  & 
       11441451.509d0, 9387934.311d0,  8559624.313d0,  7165322.549d0,  &
        4648135.263d0, 2447035.933d0,  3885009.997d0,  1006980.679d0 /)

#elif defined( TOMAS40 )

    !-----------------------------------------------------------------------
    ! TOMAS simulation w/ 40 size-resolved bins
    !-----------------------------------------------------------------------

    TOMAS_DBIN = (/                                                      &
       0.0d0      , 0.0d0      , 0.0d0     ,  0.0d0      , 0.0d0      ,  &
       0.0d0      , 0.0d0      , 0.0d0     ,  0.0d0      , 0.0d0      ,  &
       9.68859d-09, 1.22069d-08, 1.53797d-08, 1.93772d-08, 2.44137d-08,  &
       3.07594d-08, 3.87544d-08, 4.88274d-08, 6.15187d-08, 7.75087d-08,  &
       9.76549d-08, 1.23037d-07, 1.55017d-07, 1.95310d-07, 2.46075d-07,  &
       3.10035d-07, 3.90620d-07, 4.92150d-07, 6.20070d-07, 7.81239d-07,  &
       9.84300d-07, 1.24014d-06, 1.56248d-06, 1.96860d-06, 2.48028d-06,  &
       3.12496d-06, 3.93720d-06, 4.96056d-06, 6.24991d-06, 7.87440d-06 /)

    TOMAS_A    = (/                                                          &
  0.0d0        , 0.0d0        , 0.0d0        , 0.0d0       ,  0.0d0   ,      &
  0.0d0        , 0.0d0        , 0.0d0        , 0.0d0       ,  0.0d0   ,      &
  1719793.975d0, 2887719.254d0, 4086059.079d0, 5222972.121d0, 6172287.155d0, &
  6789341.855d0, 6954290.435d0, 6647842.508d0, 6030292.470d0, 5411159.039d0, &
  4920485.633d0, 4467448.678d0, 4379031.834d0, 4180592.479d0, 3836983.331d0, &
  3328339.218d0, 2675909.440d0, 1972225.823d0, 1384692.112d0, 1062343.821d0, &
  913194.1118d0, 859176.8257d0, 812688.4300d0, 719215.3301d0, 580735.2991d0, &
  418247.5535d0, 273217.6572d0, 183340.5653d0, 132174.9032d0,      0.0000d0/) 


#else

    !-----------------------------------------------------------------------
    ! TOMAS simulation w/ 30 size-resolved bins (default)
    !-----------------------------------------------------------------------

    TOMAS_DBIN = (/                                                      &
       9.68859d-09, 1.22069d-08, 1.53797d-08, 1.93772d-08, 2.44137d-08,  &
       3.07594d-08, 3.87544d-08, 4.88274d-08, 6.15187d-08, 7.75087d-08,  &
       9.76549d-08, 1.23037d-07, 1.55017d-07, 1.95310d-07, 2.46075d-07,  &
       3.10035d-07, 3.90620d-07, 4.92150d-07, 6.20070d-07, 7.81239d-07,  &
       9.84300d-07, 1.24014d-06, 1.56248d-06, 1.96860d-06, 2.48028d-06,  &
       3.12496d-06, 3.93720d-06, 4.96056d-06, 6.24991d-06, 7.87440d-06 /)

    TOMAS_A    = (/                                                          &
  1719793.975d0, 2887719.254d0, 4086059.079d0, 5222972.121d0, 6172287.155d0, &
  6789341.855d0, 6954290.435d0, 6647842.508d0, 6030292.470d0, 5411159.039d0, &
  4920485.633d0, 4467448.678d0, 4379031.834d0, 4180592.479d0, 3836983.331d0, &
  3328339.218d0, 2675909.440d0, 1972225.823d0, 1384692.112d0, 1062343.821d0, &
  913194.1118d0, 859176.8257d0, 812688.4300d0, 719215.3301d0, 580735.2991d0, &
  418247.5535d0, 273217.6572d0, 183340.5653d0, 132174.9032d0,      0.0000d0/) 

#endif

    !=======================================================================
    ! Allocate quantities depending on horizontal resolution
    !=======================================================================
#if defined( GRID4x5  )

    !-----------------------------------------------------------------------
    ! TOMAS simulations at 4 x 5 global resolution
    !-----------------------------------------------------------------------
    TOMAS_COEF = 1.d0

#elif defined( GRID2x25 )

    !-----------------------------------------------------------------------
    ! TOMAS simulations at 2 x 2.5 global resolution
    !-----------------------------------------------------------------------
    TOMAS_COEF = 1.d0

#elif defined( GRID1x1  ) && defined( NESTED_CH )

    !-----------------------------------------------------------------------
    ! TOMAS simulations at 1 x 1 nested-grid resolution
    !
    ! NOTE: Monthly emission over the China nested-grid domain is about
    !       1.3 times higher compared to the same domain in 4x5 resolution
    !       Thus applying 1/1.30 = 0.77 factor to correct the emission.
    !
    ! NOTE: May need to update this for other nested-grids
    !-----------------------------------------------------------------------
    TOMAS_COEF = 0.77d0

#else
    MSG = 'Adjust TOMAS_SeaSalt emiss coeff (TOMAS_COEF) for your model res: SRCSALT30: hcox_TOMAS_SeaSalt_mod.F90'
      call HCO_ERROR( MSG, RC )
#endif

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

!    print*, 'ENABLE MODULE'
    ! Enable module
    ExtState%TOMAS_SeaSalt     = .TRUE.
    

    ! Return w/ success
!    IF ( ALLOCATED( HcoIDs   ) ) DEALLOCATE( HcoIDs   )
    IF ( ALLOCATED( SpcNames ) ) DEALLOCATE( SpcNames )
    CALL HCO_LEAVE ( RC ) 
 
  END SUBROUTINE HCOX_TOMAS_SeaSalt_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_TOMAS_SeaSalt_Final 
!
! !DESCRIPTION: Subroutine HcoX\_TOMAS_SeaSalt\_Final deallocates 
!  all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_TOMAS_SeaSalt_Final
!
! !REVISION HISTORY:
!  15 Dec 2013 - C. Keller   - Initial version
!  20 May 2015 - J. Kodros   - Deallocate HcoIDs, TC1, TC2 arrays
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! HCOX_TOMAS_SeaSalt_Final begins here!
    !=================================================================
    IF ( ALLOCATED( TOMAS_DBIN ) ) DEALLOCATE( TOMAS_DBIN )
    IF ( ALLOCATED( TOMAS_A    ) ) DEALLOCATE( TOMAS_A    )
    IF ( ALLOCATED( HcoIDs     ) ) DEALLOCATE( HcoIDs     )
    IF ( ALLOCATED( TC1        ) ) DEALLOCATE( TC1        )
    IF ( ALLOCATED( TC2        ) ) DEALLOCATE( TC2        )

  END SUBROUTINE HCOX_TOMAS_SeaSalt_Final
!EOC
END MODULE HCOX_TOMAS_SeaSalt_Mod
#endif
