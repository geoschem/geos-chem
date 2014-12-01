!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: loadHgDeposition
!
! !DESCRIPTION: This code reads in output from geos chem 4x5 degree, 
!  converts it to 1x1 degree, then converts it to (n\_veg,1) for CASA
!\\
!\\
! !INTERFACE:
!
SUBROUTINE loadHgDeposition(LCPLE, DD_Hg0, DD_HgII, WD_HgII)
!
! !USES:
!      
  USE defineConstants
  USE loadCASAinput
  USE defineArrays
  USE CasaRegridModule

  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)
    
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
  LOGICAL, INTENT(IN)           :: LCPLE
  
  REAL(fp), INTENT(IN), OPTIONAL  :: DD_Hg0(72, 46), DD_HgII(72, 46)
  REAL(fp), INTENT(IN), OPTIONAL  :: WD_HgII(72, 46)
!
! !REVISION HISTORY:
! 15 Dec 2009 - C. Carouge  - Add arguments for coupling with GEOS-Chem
  !                         - Change format of emission years file to 
!                             facilitate restart.
! 01 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  REAL(fp) ::                 Hg0dryGEOS(72,46)
  REAL(fp) ::                 HgIIdryGEOS(72,46)
  REAL(fp) ::                 HgIIwetGEOS(72,46)
  REAL(fp) ::                 area_pix(72,46)
  REAL(fp) ::                 Hg0dry1x1(360,180)
  REAL(fp) ::                 HgIIdry1x1(360,180)
  REAL(fp) ::                 HgIIwet1x1(360,180)
  REAL(fp) ::                 LAT_LONG(3,n_veg)
  INTEGER ::                i, j, ios, v,w, k
  CHARACTER(LEN=f_len+14) :: filename1
  CHARACTER(LEN=f_len+6) :: filename2
  CHARACTER(3), DIMENSION(12) :: months
  INTEGER, DIMENSION(2):: em_years
  INTEGER              :: year

  CHARACTER(4)   :: read_year

  LOGICAL, SAVE :: FIRST=.TRUE.  ! Indicate if we read the data for the first
                                 ! time. GTMM offline only.

  !=====================================================================
  ! loadHgDepostion begins here !
  !=====================================================================

  filename2(1:f_len)=filepath
  filename2(f_len+1:f_len+6)='geosm2'
  OPEN(UNIT=3, FILE=filename2, FORM='FORMATTED',IOSTAT=ios)
  READ(3,*) area_pix
  CLOSE(3)
  
  filename2(f_len+1:f_len+6)='months'
  OPEN(UNIT=3, FILE=filename2, FORM='FORMATTED', IOSTAT=ios)
  READ(3,*) months
  CLOSE(3)
  
  filename2(f_len+1:f_len+6)='emyear'
  OPEN(UNIT=3, FILE=filename2, FORM='FORMATTED', IOSTAT=ios)
  READ(3,*) em_years
  CLOSE(3)

  IF ( LCPLE ) THEN
     IF (.NOT.( PRESENT( DD_Hg0  ) )) THEN
        stop 'Dry deposition of Hg0 is missing'

     ELSE
        ! Use deposition from GEOS-Chem
        Hg0dryGEOS = DD_Hg0
        
     ENDIF
     
     IF (.NOT.( PRESENT( DD_HgII  ) )) THEN
        stop 'Dry deposition of HgII is missing'
        
     ELSE
        ! Use deposition from GEOS-Chem
        HgIIdryGEOS = DD_HgII
        
     ENDIF
     
     IF (.NOT.( PRESENT( WD_HgII  ) )) THEN
        stop 'Wet deposition of HgII is missing'
        
     ELSE
        ! Use deposition from GEOS-Chem
        HgIIwetGEOS = WD_HgII
     ENDIF
     
100  FORMAT('Warning: Read deposition from archive for: ',a)
     
     !convert units
     !wet dep is in units of kg/s
!$OMP PARALLEL DO    &
!$OMP PRIVATE(i, j)
     DO i=1,46
     DO j=1,72
        HgIIwetGEOS(j,i)=HgIIwetGEOS(j,i)/area_pix(j,i)
        
        HgIIwetGEOS(j,i)=HgIIwetGEOS(j,i)*1000.0e+0_fp*2629743.83e+0_fp
        !now in units of g/m2/mo
        
        !dry dep is in unitls of molec/cm2/s
        HgIIdryGEOS(j,i)=(HgIIdryGEOS(j,i)*(10000.0e+0_fp*2629743.83e+0_fp* &
                200.59e+0_fp))/(6.022e+23_fp)
        Hg0dryGEOS(j,i)=(Hg0dryGEOS(j,i)*(10000.0e+0_fp*2629743.83e+0_fp* &
                200.59e+0_fp))/(6.022e+23_fp)
        !now in units of g/m2/mo
     END DO
     END DO
!$OMP END PARALLEL DO  

!!!!First step - regrid geos chem 4x5 to 1x1
     CALL CasaRegridInit
     CALL regrid4x5to1x1(1, Hg0dryGEOS, Hg0dry1x1)
     CALL regrid4x5to1x1(1, HgIIdryGEOS, HgIIdry1x1)
     CALL regrid4x5to1x1(1, HgIIwetGEOS, HgIIwet1x1)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Now convert 1x1 to n_veg x1 for use by CASA
     Hg0dry(:,:)=0.0e+0_fp
     HgIIdry(:,:)=0.0e+0_fp
     HgIIwet(:,:)=0.0e+0_fp
     Hg0dry(:,:) = maskfile(Hg0dry1x1, mask2)
     HgIIdry(:,:) = maskfile(HgIIdry1x1, mask2)
     HgIIwet(:,:) = maskfile(HgIIwet1x1, mask2)

  ELSE
     
     ! Find the year for the emissions. (ccc, 12/15/09)
     IF ( yr <= preindYear ) THEN
!        read_year=em_years(1)
        write(read_year,'(i4)') em_years(1)
     ELSE
        year=em_years(1)+yr-preindYear      ! Number of years in 
                                            ! industrialization period
        year=MIN(year,em_years(2))
        write(read_year,'(i4)'), year

        IF ( mo == 1 ) FIRST = .TRUE.
     ENDIF

     IF ( FIRST ) THEN  ! Read the full year at once.

        filename2(1:f_len)=filepath
        filename2(f_len+1:f_len+6)='latlon'
        OPEN(UNIT=3, FILE=filename2, FORM='FORMATTED', IOSTAT=ios)
        READ(3,FMT="(3F9.2)") LAT_LONG
        CLOSE(3)

        DO k=1,12 

           filename1(1:f_len)=filepath
           filename1(f_len+1:f_len+6)='Hg0dry'
           filename1(f_len+7:f_len+9)=months(k)
!--- Previous to (ccc, 12/15/09)
!     filename1(f_len+10:f_len+14)=em_years(yr)
           filename1(f_len+10:f_len+14)=read_year
     
           OPEN(UNIT=3, FILE=filename1, IOSTAT=ios, FORM="FORMATTED")
           READ(3,FMT="(72E12.5)") Hg0dryGEOS
           CLOSE(3)
     
           filename1(f_len+1:f_len+6)='Hg2dry'
           OPEN(UNIT=3, FILE=filename1, IOSTAT=ios, FORM="FORMATTED")
           READ(3,FMT="(72E12.5)") HgIIdryGEOS
           CLOSE(3)
     !------------------------------------------------------------------------------

           filename1(f_len+1:f_len+6)='Hg2wet'
           OPEN(UNIT=3, FILE=filename1, IOSTAT=ios, FORM="FORMATTED")
           READ(3,FMT="(72E12.5)") HgIIwetGEOS
           CLOSE(3)
     
           !convert units
           !wet dep is in units of kg/s
!$OMP PARALLEL DO    &
!$OMP PRIVATE(i, j)
           DO i=1,46
           DO j=1,72
              HgIIwetGEOS(j,i)=HgIIwetGEOS(j,i)/area_pix(j,i)
              
              HgIIwetGEOS(j,i)=HgIIwetGEOS(j,i)*1000.0e+0_fp*2629743.83e+0_fp
              !now in units of g/m2/mo
              
              !dry dep is in unitls of molec/cm2/s
              HgIIdryGEOS(j,i)=(HgIIdryGEOS(j,i)*(10000.0e+0_fp*2629743.83e+0_fp &
                *200.59e+0_fp))/(6.022e+23_fp)
              Hg0dryGEOS(j,i)=(Hg0dryGEOS(j,i)*(10000.0e+0_fp*2629743.83e+0_fp &
                *200.59e+0_fp))/(6.022e+23_fp)
              !now in units of g/m2/mo
           END DO
           END DO
!$OMP END PARALLEL DO  

!!!!First step - regrid geos chem 4x5 to 1x1
           CALL CasaRegridInit
           CALL regrid4x5to1x1(1, Hg0dryGEOS, Hg0dry1x1)
           CALL regrid4x5to1x1(1, HgIIdryGEOS, HgIIdry1x1)
           CALL regrid4x5to1x1(1, HgIIwetGEOS, HgIIwet1x1)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!Now convert 1x1 to n_veg x1 for use by CASA
           Hg0dry_mo(:,:,k)=0.0e+0_fp
           HgIIdry_mo(:,:,k)=0.0e+0_fp
           HgIIwet_mo(:,:,k)=0.0e+0_fp
           Hg0dry_mo(:,:,k) = maskfile(Hg0dry1x1, mask2)
           HgIIdry_mo(:,:,k) = maskfile(HgIIdry1x1, mask2)
           HgIIwet_mo(:,:,k) = maskfile(HgIIwet1x1, mask2)
        ENDDO
        FIRST = .FALSE.
     ENDIF

     ! Store the data for the current month
     Hg0dry = Hg0dry_mo(:,:,mo)
     HgIIdry = HgIIdry_mo(:,:,mo)
     HgIIwet = HgIIwet_mo(:,:,mo)

  ENDIF
  
END SUBROUTINE loadHgDeposition
!EOC
