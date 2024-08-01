#include "MAPL_Generic.h"
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: geos_TaggedSpecies
!
! !DESCRIPTION: Module with routines and variables to handle tagged species. 
! For now, this module only covers tagged NO(x). It can handle up to 100 
! individual tagged NO species. The number of tagged species is specified
! in the input file 'GEOSChem_TaggedNOx.rc'. If partner_tag is enabled, a
! tagged NO2 species is also used and tagging will be done as NOx (NO+NO2).

! The tagged species names are hardcoded to 'NOTAG<XX>' and NO2TAG<XX>',
! respectively. To run with the tagged species, they need to be added to 
! the list of advected species in geoschem_config.yaml. Also, for each 
! tagged species, there needs to be corresponding entry in the species 
! database file (species_database.yml), e.g.:
!
!NOTAG1:
!  Background_VV: 1.0e-30
!  Formula: 'NO'
!  FullName: Nitrogen oxide tag 1
!  Is_Advected: true
!  Is_Gas: true
!  Is_Photolysis: false
!  MW_g: 30.01
!
!NO2TAG1:
!  Background_VV: 1.0e-30
!  DD_F0: 0.1
!  DD_Hstar: 1.0e-2
!  Formula: NO2
!  FullName: Nitrogen dioxide tag 1
!  Is_Advected: true
!  Is_DryDep: true
!  Is_Gas: true
!  Is_Photolysis: false
!  MW_g: 46.01
!
! Chemistry prod/loss rates are inherited from the parent species (NO/NO2),
! all other processes are performed at the tagged species level. Emissions 
! for each tagged species need to be assigned in the HEMCO configuration 
! file. For example, to assign NO biomass burning emissions to tagged species
! #1 (NOTAG1): 
!0 QFED_NO_TF_SFC   $ROOT/QFED/v2014-09/$YYYY/$MM/qfed2.emis_no.006.$YYYY$MM$DD.nc4 biomass_tf 2000-2018/1-12/1-31/0 C xyL=1:PBL kg/m2/s NO     75/311/545/592 8 1
!0 QFED_NOT1_TF_SFC -                                                               -          -                     - -         -       NOTAG1 75/311/545/592 8 1 
! 
!\\
!\\
! !INTERFACE:
!
MODULE GEOS_TaggedSpecies
!
! !USES:
!
  ! MAPL/ESMF
  USE ESMF     
  USE MAPL_Mod 
  ! GEOS-Chem stuff
  USE Precision_Mod
  USE ErrCode_Mod                                    ! Error numbers
  USE PHYSCONSTANTS
  USE Input_Opt_Mod,         ONLY : OptInput
  USE State_Chm_Mod,         ONLY : ChmState, Ind_   ! Chemistry State obj
  USE State_Diag_Mod,        ONLY : DgnState         ! Diagnostics State obj

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC   :: Init_TaggedSpecies
  PUBLIC   :: Run_TaggedSpecies
  PUBLIC   :: Finalize_TaggedSpecies
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !PRIVATE TYPES:
!
  ! Tagged tracers:
  CHARACTER(LEN=128), PARAMETER  :: TaggedConfigFile = 'GEOSChem_TaggedNOx.rc'

  INTEGER, PARAMETER    :: MaxTag = 100
  INTEGER               :: nTagged
  INTEGER, ALLOCATABLE  :: TagID(:)                      ! For Tagged NO
  INTEGER, ALLOCATABLE  :: ParentID_GCC(:)
  INTEGER, ALLOCATABLE  :: ParentID_KPP(:)
  INTEGER, ALLOCATABLE  :: TagIDb(:)                     ! For Tagged NO2
  INTEGER, ALLOCATABLE  :: ParentIDb_GCC(:)
  INTEGER, ALLOCATABLE  :: ParentIDb_KPP(:)
!
! !REVISION HISTORY:
!  19 Mar 2024 - C. Keller / P. Wales - brought into GCv14 from v12 
!  See https://github.com/geoschem/geos-chem for full history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_TaggedSpecies 
!
! !DESCRIPTION: Initialize tagged species chemistry by reading the information
! from GEOSChem_TaggedNOx.rc 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_TaggedSpecies( Input_Opt, State_Chm, State_Diag, RC )
!
! !USE:
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput)                             :: Input_Opt
    TYPE(ChmState)                             :: State_Chm
    TYPE(DgnState)                             :: State_Diag
    INTEGER,             INTENT(INOUT)         :: RC        ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  19 Mar 2024 - C. Keller / P. Wales - brought into GCv14 from v12 
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    TYPE(ESMF_Config)                     :: tagCF     ! ESMF Config obj
    LOGICAL                               :: FileExists
    CHARACTER(LEN=255)                    :: ParentName, ParentNameb
    CHARACTER(LEN=64)                     :: TagName, TagNameb
    INTEGER                               :: KppId, TagInd, N, PartnerTag

    __Iam__('Init_TaggedSpecies')

    !--------------------------------------------------------------------
    ! Initialize tagged tracer chemistry
    !--------------------------------------------------------------------

    ! Initialize 
    nTagged          = 0
                                                   
    TagName = 'NaN'
    ParentName = 'NO'          ! eventually allow this to be read in as a list
    TagNameb = 'NaN'
    ParentNameb = 'NO2'        ! chem that cycles with ParentName
   
    ! Check if file exists 
    INQUIRE( FILE=TRIM(TaggedConfigFile), EXIST=FileExists )

    ! Do the following only if file exits
    IF ( FileExists ) THEN       

       ! Verbose
       IF ( Input_Opt%amIRoot ) THEN
          WRITE(*,*) 'Reading tagged NOx tracer information from '//TRIM(TaggedConfigFile)
       ENDIF 

       tagCF = ESMF_ConfigCreate( __RC__ )
       Call ESMF_ConfigLoadFile(tagCF, TRIM(TaggedConfigFile), __RC__ )

       Call ESMF_ConfigGetAttribute(tagCF, PartnerTag, Label = 'partner_tag:', &
                                    Default = 0, __RC__ )
       Call ESMF_ConfigGetAttribute(tagCF, nTagged, Label  ='nchem_tag:', &
                                    Default = 0, __RC__ )          
       ASSERT_( nTagged <= MaxTag )

       IF ( nTagged > 0) THEN

          ! Now that we know ntagged, can allocate vector arrays 
          ALLOCATE(TagID(nTagged), ParentID_GCC(nTagged), ParentID_KPP(nTagged), &
                   TagIDb(nTagged), ParentIDb_GCC(nTagged), ParentIDb_KPP(nTagged))

          ! Get GEOS-Chem and KPP index of the parent species. This is currently hardcoded
          ! to be NO everywhere, so can do it outside of the nTagged loop below
          ParentID_GCC(:)  = Ind_(TRIM(ParentName))
          ASSERT_(ParentID_GCC(1)>0)

          DO KppID = 1, State_Chm%nKppSpc
             IF ( State_Chm%Map_KppSpc(KppId) == ParentID_GCC(1) ) THEN
                ParentID_KPP(:) = KppID
                EXIT
             ENDIF
          ENDDO
          ! Make sure that KPP parent ID is valid
          ASSERT_(ParentID_KPP(1)>0)

          ! Get GEOS-Chem and KPP index of the parent species for the partner tag. 
          ! This is currently hardcoded to be NO2 everywhere, so can do it outside 
          ! of the nTagged loop below
          IF ( PartnerTag > 0) THEN
             ParentIDb_GCC(:) = Ind_(TRIM(ParentNameb))
             ASSERT_(ParentIDb_GCC(1)>0)

             DO KppID = 1, State_Chm%nKppSpc
                IF ( State_Chm%Map_KppSpc(KppId) == ParentIDb_GCC(1) ) THEN
                   ParentIDb_KPP(:) = KppID
                   EXIT
                ENDIF
             ENDDO
             ASSERT_(ParentIDb_KPP(1)>0)
          ENDIF

          ! Loop over all tagged species and assign species IDs
          DO TagInd = 1, nTagged
     
             ! Species names are currently hardcoded to NO_TAG1, NO_TAG2, etc. 
             WRITE(TagName,"(A5,I0)") "NOTAG", TagInd
             N = Ind_(TRIM(TagName))
             IF ( N > 0 ) THEN
                TagID(TagInd) = N
                IF ( Input_Opt%amIRoot ) THEN
                   WRITE(*,*) 'Tagged species: will apply chem P/L rates of '//TRIM(ParentName)//' to '//TRIM(TagName)
                ENDIF

                ! Check for partner tag, hardcoded to NO2_TAG1, NO2_TAG2, etc.
                IF (PartnerTag > 0) THEN
                   WRITE(TagNameb,"(A6,I0)") "NO2TAG", TagInd
                   TagIDb(TagInd) = Ind_(TRIM(TagNameb))
                   IF ( Input_Opt%amIRoot ) THEN
                      WRITE(*,*) 'Tagged species: will use partner species '//TRIM(TagNameb)//' for '//TRIM(TagName)
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ENDIF ! nTagged>0
    ENDIF ! FileExists

    ! All done
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE Init_TaggedSpecies 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Run_TaggedSpecies 
!
! !DESCRIPTION: Run tagged species chemistry. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Run_TaggedSpecies( I, J, L, C, PRESS, TEMP, State_Chm, RC ) 
!
! !USE:
!
    USE GcKpp_Parameters
    USE ERROR_MOD
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(IN)            :: I,J,L
    REAL(dp),            INTENT(IN)            :: C(NSPEC)
    REAL(dp),            INTENT(IN)            :: PRESS 
    REAL(dp),            INTENT(IN)            :: TEMP 
    TYPE(ChmState)                             :: State_Chm
    INTEGER,             INTENT(INOUT)         :: RC        ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  19 Mar 2024 - C. Keller / P. Wales - brought into GCv14 from v12 
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER  :: N, itagID, SpcID, KppID, iTagIDb, SpcIDb, KppIDb
    REAL(fp) :: Tag0, Tag1, Parent0, Parent1, Prat, TagThres, NOxrat

    __Iam__('Run_TaggedSpecies')

    !--------------------------------------------------------------------
    ! Run tagged tracer chemistry
    !--------------------------------------------------------------------

    !=====================================================================
    ! Check for tagged tracers and adjust those before updating the
    ! Species arrays
    !=====================================================================
    IF ( nTagged > 0 ) THEN

      ! Threshold is 1.E-14 mol/mol, convert to molecular density
      ! density at STP * 298 K / 1013 hPa = 7.24E4
      TagThres = 7.24e4_fp * PRESS / TEMP

      DO N = 1, nTagged
        iTagID  = TagID(N)
        SpcID   = ParentID_GCC(N)
        KppID   = ParentID_KPP(N)
        iTagIDb = TagIDb(N)
        SpcIDb  = ParentIDb_GCC(N)
        KppIDb  = ParentIDb_KPP(N)

        !!! If there is a partner compound for iTagID
        IF ( iTagIDb > 0) THEN
          Tag0    = State_Chm%Species(iTagID)%Conc(I,J,L) + &
                    State_Chm%Species(iTagIDb)%Conc(I,J,L)
          Parent0 = State_Chm%Species(SpcID)%Conc(I,J,L) + &
                    State_Chm%Species(SpcIDb)%Conc(I,J,L)

          IF ( Tag0 > TagThres .and. Parent0 > TagThres) THEN
            Prat    = SAFE_DIV( Tag0, Parent0, 1.0_fp, 1.0_fp, 0.0_fp )
            Prat    = MAX(MIN(Prat,1.0),0.0)

            Parent1 = REAL(MAX(C(KppID),0.0_dp) + MAX(C(KppIDb),0.0_dp),kind=fp)
            NOxrat  = SAFE_DIV(MAX(C(KppID),0.0_dp), Parent1, 1.0_fp, &
                      1.0_fp, 0.0_fp )
            NOxrat  = MAX(MIN(NOxrat,1.0),0.0)

            Tag1 = REAL(Tag0 + ( Parent1 - Parent0 ) * Prat, kind=fp)
            State_Chm%Species(iTagID)%Conc(I,J,L)  = MAX(Tag1 * NOxrat, 0.0_fp)
            State_Chm%Species(iTagIDb)%Conc(I,J,L) = MAX(Tag1 * (1.0_fp - NOxrat), &
                                                     0.0_fp)
          ELSE
            State_Chm%Species(iTagID)%Conc(I,J,L) = 0.0_fp
            State_Chm%Species(iTagIDb)%Conc(I,J,L) = 0.0_fp
          ENDIF

        ! only scale iTagID
        ELSE
          Tag0    = State_Chm%Species(iTagID)%Conc(I,J,L)
          Parent0 = State_Chm%Species(SpcID)%Conc(I,J,L)

          IF ( Tag0 > TagThres .and. Parent0 > TagThres) THEN
            Prat    = SAFE_DIV( Tag0, Parent0, 1.0_fp, 1.0_fp, 0.0_fp )
            Prat    = MAX(MIN(Prat,1.0),0.0)
            Parent1 = REAL(MAX(C(KppID),0.0_dp),kind=fp)

            State_Chm%Species(iTagID)%Conc(I,J,L) = Tag0 + ( Parent1 - Parent0 ) * Prat
            State_Chm%Species(iTagID)%Conc(I,J,L) = MAX(State_Chm%Species(iTagID)%Conc(I,J,L), &
                                                0.0_fp)
          ELSE
            State_Chm%Species(iTagID)%Conc(I,J,L) = 0.0_fp
          ENDIF
        ENDIF
      ENDDO
    ENDIF                                                   

    ! All done
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE Run_TaggedSpecies 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finalize_TaggedSpecies 
!
! !DESCRIPTION: Finalize tagged species arrays. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Finalize_TaggedSpecies( RC ) 
!
! !USE:
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT)         :: RC        ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  19 Mar 2024 - C. Keller / P. Wales - brought into GCv14 from v12 
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    __Iam__('Finalize_TaggedSpecies')

    IF ( ALLOCATED(TagID        ) ) DEALLOCATE(TagID)
    IF ( ALLOCATED(ParentID_GCC ) ) DEALLOCATE(ParentID_GCC)
    IF ( ALLOCATED(ParentID_KPP ) ) DEALLOCATE(ParentID_KPP)
    IF ( ALLOCATED(TagIDb       ) ) DEALLOCATE(TagIDb)
    IF ( ALLOCATED(ParentIDb_GCC) ) DEALLOCATE(ParentIDb_GCC)
    IF ( ALLOCATED(ParentIDb_KPP) ) DEALLOCATE(ParentIDb_KPP)

    ! All done
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE Finalize_TaggedSpecies 
!EOC
END MODULE GEOS_TaggedSpecies
