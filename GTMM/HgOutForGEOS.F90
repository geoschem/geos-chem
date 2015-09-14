!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HgOutForGEOS
!
! !DESCRIPTION: Subroutine HgOutForGEOS converts the (n\_veg,1) data to 1x1, 
!  then to 4x5 grid then writes out the file for use by GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
SUBROUTINE HgOutForGEOS(LCPLE, Hg0reemit)
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE defineArrays
  USE CasaRegridModule
  USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)
  
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
  LOGICAL, INTENT(IN) :: LCPLE
!
! !OUTPUT PARAMETERS:
!
  REAL(fp), INTENT(OUT) :: Hg0reemit(72, 46)
!
! !REVISION HISTORY:
!  09 Jul 2010 - C. Carouge  - Modified to couple with GEOS-Chem
!  25 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  REAL(fp) ::                 Hg0outGEOS(72,46)
  REAL(fp) ::                 Hg0out1x1(360,180)
  REAL(fp) ::                 Hg0CombGEOS(72,46)
  REAL(fp) ::                 Hg0Comb1x1(360,180)
  REAL(fp) ::                 HgPCombGEOS(72,46)
  REAL(fp) ::                 HgPComb1x1(360,180)
  REAL(fp) ::                 Hg0PhotGEOS(72,46)
  REAL(fp) ::                 Hg0Phot1x1(360,180)
  REAL(fp) ::                 Hg0VoltGEOS(72,46)
  REAL(fp) ::                 Hg0Volt1x1(360,180)
  REAL(fp) ::                 Hg0RespGEOS(72,46)
  REAL(fp) ::                 Hg0Resp1x1(360,180) 
  REAL(fp) ::                 LAT_LONG(3,n_veg)
  REAL(fp) ::                 pi=3.14159265e+0_fp
  REAL(fp) ::                 radius=6378140.00e+0_fp
  REAL(fp) ::                 g, a, apixel
  REAL(fp), dimension(360,180) ::                 gridAreaf
  REAL(fp), dimension(72,46)   ::                 gridAreag 
  CHARACTER(LEN=f_len+6) :: filename1
  CHARACTER(LEN=f_len_output+9) :: filename2 
  CHARACTER(3), DIMENSION(12) :: months 
  INTEGER ::                j,i, v, w
  
  filename1(1:f_len)=filepath
  filename1(f_len+1:f_len+6)='months'
  OPEN(UNIT=3, FILE=filename1, FORM='FORMATTED')
  READ(3,*) months
  CLOSE(3) 
  filename1(f_len+1:f_len+6)='geosm2'
  OPEN(UNIT=3, FILE=filename1, FORM='FORMATTED')
  READ(3,*) gridareag
  CLOSE(3)
  
  
!$OMP PARALLEL         &
!$OMP DEFAULT(SHARED)  
!$OMP WORKSHARE
  Hg0out(:,1)=0.0e+0_fp
  Hg0out(:,1)=(wresp_hg(:,1)*frac_tree(:,1))+  &
       (hresp_hg(:,1)*frac_herb(:,1))+        &
       reemitted(:,1)+photoreduced(:,1)
!$OMP END WORKSHARE
!$OMP END PARALLEL      

  filename1(f_len+1:f_len+6)='latlon'
  OPEN(UNIT=5, FILE=filename1, FORM='FORMATTED')
  READ(5,*) LAT_LONG
  CLOSE(5)
  

  ! Initialize arrays to 0. (ccc, 10/21/09)
  Hg0out1x1   = 0e+0_fp
  Hg0Phot1x1  = 0e+0_fp
  Hg0Volt1x1  = 0e+0_fp
  Hg0Resp1x1  = 0e+0_fp
  Hg0outGEOS  = 0e+0_fp
  Hg0PhotGEOS = 0e+0_fp
  Hg0VoltGEOS = 0e+0_fp
  Hg0RespGEOS = 0e+0_fp

!$OMP PARALLEL DO      &
!$OMP DEFAULT(SHARED)  &
!$OMP PRIVATE(i, v, w)  
  DO i=1, n_veg
     v=lat_long(1,i)
     w=lat_long(2,i)
     Hg0out1x1(v,w)=Hg0out(i,1)
     Hg0Phot1x1(v,w)=photoreduced(i,1)
     Hg0Volt1x1(v,w)=reemitted(i,1)
     Hg0Resp1x1(v,w)=(wresp_hg(i,1)*frac_tree(i,1))+(hresp_hg(i,1)*frac_herb(i,1)) 
  END DO
!$OMP END PARALLEL DO    
      
  CALL CasaRegridInit
  CALL regrid1x1to4x5(1, Hg0out1x1, Hg0outGEOS)
  CALL regrid1x1to4x5(1, Hg0Phot1x1, Hg0PhotGEOS)
  CALL regrid1x1to4x5(1, Hg0Volt1x1, Hg0VoltGEOS)
  CALL regrid1x1to4x5(1, Hg0Resp1x1, Hg0RespGEOS)


  !!CONVERT from g/m2/mo --> kg/s

!$OMP PARALLEL DO      &
!$OMP DEFAULT(SHARED)  &
!$OMP PRIVATE(i, j)  
  DO j=1, 46
     DO i=1, 72 
        Hg0outGEOS(i,j)=Hg0outGEOS(i,j)*gridAreag(i,j)*(1e+0_fp/1000e+0_fp)*   &
             (1e+0_fp/2629743.8e+0_fp)
        Hg0PhotGEOS(i,j)=Hg0PhotGEOS(i,j)*gridAreag(i,j)*(1e+0_fp/1000e+0_fp)* &
             (1e+0_fp/2629743.8e+0_fp)
        Hg0VoltGEOS(i,j)=Hg0VoltGEOS(i,j)*gridAreag(i,j)*(1e+0_fp/1000e+0_fp)* &
             (1e+0_fp/2629743.8e+0_fp)
        Hg0RespGEOS(i,j)=Hg0RespGEOS(i,j)*gridAreag(i,j)*(1e+0_fp/1000e+0_fp)* &
             (1e+0_fp/2629743.8e+0_fp)
     END DO
  END DO
!$OMP END PARALLEL DO    

  filename2(1:f_len_output)=outputpath
  filename2(f_len_output+1:f_len_output+6)='Hg0out'
  filename2(f_len_output+7:f_len_output+9)=months(mo)
 
  ! Output the total Hg0 re-emitted to GEOS-Chem. 
  Hg0reemit=Hg0outGEOS
  
  IF ( .NOT. LCPLE ) THEN
     OPEN(UNIT=5, FILE=filename2)
     WRITE(5,*) Hg0outGEOS
     CLOSE(5)
     
     filename2(f_len_output+1:f_len_output+6)='Hg0Pht'
     OPEN(UNIT=5, FILE=filename2)
     WRITE(5,*) Hg0PhotGEOS
     CLOSE(5)
     
     filename2(f_len_output+1:f_len_output+6)='Hg0Vot'
     OPEN(UNIT=5, FILE=filename2)
     WRITE(5,*) Hg0VoltGEOS
     CLOSE(5)
     
     filename2(f_len_output+1:f_len_output+6)='Hg0Res'
     OPEN(UNIT=5, FILE=filename2)
     WRITE(5,*) Hg0RespGEOS
     CLOSE(5)
     
  ENDIF

  !ok, now data has been written out for geoschem
  
END SUBROUTINE HgOutForGEOS
!EOC
