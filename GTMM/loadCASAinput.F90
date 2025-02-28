!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: loadCASAinput
!
! !DESCRIPTION: Loads input files.
!\\
!\\
! !INTERFACE:
!
MODULE loadCASAinput
!
! !USES:
!
  USE defineConstants
  USE CasaRegridModule

  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
!
!CONTINUOUS FIELD MAPS
  REAL(fp), ALLOCATABLE :: perc_tree(:,:)
  REAL(fp), ALLOCATABLE :: perc_herb(:,:) 
  
!RESIZED CONTINUOUS FIELD MAPS
  REAL(fp), ALLOCATABLE :: perc_tree1(:,:)
  REAL(fp), ALLOCATABLE :: perc_herb1(:,:)
  REAL(fp), ALLOCATABLE :: frac_tree(:,:)
  REAL(fp), ALLOCATABLE :: frac_herb(:,:)
  REAL(fp), ALLOCATABLE :: frac_veg(:,:)

!CLIMATE FILES
  REAL(fp), ALLOCATABLE :: airt(:,:,:)   !monthly air temperature
  REAL(fp), ALLOCATABLE :: ppt(:,:,:)    !monthly precipitation
  REAL(fp), ALLOCATABLE :: solrad(:,:,:) !monthly solar radiation, 
                                       !Taken from Bishop and Rossow
                                       !average 83-99 from Jim Collatz

  REAL(fp), ALLOCATABLE :: NDVI(:,:,:)   !monthly fraction PAR
  REAL(fp), ALLOCATABLE :: BF(:,:,:)     !fraction gridcell tha burns
  REAL(fp), ALLOCATABLE :: ppt_mo(:,:)   !sum of all precip in each mo
   
!RESIZED CLIMATE FILES
  REAL(fp), ALLOCATABLE :: airt1(:,:)
  REAL(fp), ALLOCATABLE :: ppt1(:,:)
  REAL(fp), ALLOCATABLE :: solrad1(:,:)
  REAL(fp), ALLOCATABLE :: NDVI1(:,:)
  REAL(fp), ALLOCATABLE :: BF1(:,:)
  REAL(fp), ALLOCATABLE :: maxt(:,:)
  REAL(fp), ALLOCATABLE :: mint(:,:)
 
!OTHER FILES
  REAL(fp), ALLOCATABLE :: soiltext(:,:)  ! soil type
  REAL(fp), ALLOCATABLE :: veg(:,:)       ! vegetation map
  REAL(fp), ALLOCATABLE :: fuelneed(:,:)  ! fuelwood needed 
                                        !   per capita
  REAL(fp), ALLOCATABLE :: popdens(:,:)   ! population density (/m2)
  REAL(fp), ALLOCATABLE :: gridAreaa(:,:)
  REAL(fp), ALLOCATABLE :: gridAreab(:,:)
   
!RESIZED OTHER FILES
  REAL(fp), ALLOCATABLE :: soiltext1(:,:)
  REAL(fp), ALLOCATABLE :: veg1(:,:)
  REAL(fp), ALLOCATABLE :: fuelneed1(:,:)
  REAL(fp), ALLOCATABLE :: popdens1(:,:)

   
  REAL(fp), ALLOCATABLE :: mask2(:,:)
!
! !REMARKS: For global studies, these files should be given 
!  as 180 by 360 matrix for 1x1 degree, 360 by 720 for half by half 
!  etc.  This routine will construct a 'mask' with vegetated gridcells
!  and will reshape them into an X by 1 matrix, where X is the total
!  number of vegetated grid cells.  
!  The files which contain monthly varying parameters (i.e. precip)
!  should be given as (for 1x1) 180x360x12, and will be reshaped to
!  an X by 12 matrix, one column for each month
!
! !REVISION HISTORY:
! 9 July 2010 - C. Carouge  - Modified for coupled simulations with GEOS-Chem
!                             or to restart offline simulations. 
! 01 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: load_data
!
! !DESCRIPTION: Reads input data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE load_data(LCPLE)
!    
! !USES:
!
    USE defineConstants
!
! !INPUT PARAMETERS:
!
   LOGICAL, INTENT(IN)  :: LCPLE
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   integer   ::   i, j, k, l, m, n, b, c, ios
   character(len=f_len+8)                :: filename
   real(fp), dimension(72,46)              ::   geos
   real(fp), dimension(columns, rows)      ::   geos1x1
   character(3), dimension(12)           ::   months     ! Months names
   
   !<<<<<<<< ALLOCATE ARRAYS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   

   IF (.NOT. ALLOCATED(perc_tree) ) &
        ALLOCATE (perc_tree(columns, rows), stat=ios)
   
   IF (.NOT. ALLOCATED(perc_herb) ) &
        ALLOCATE (perc_herb(columns, rows), stat=ios)
   
   IF (.NOT. ALLOCATED(airt) ) &
        ALLOCATE (airt(columns, rows,12), stat=ios)
   
   IF (.NOT. ALLOCATED(ppt) ) &
        ALLOCATE (ppt(columns, rows,12), stat=ios)
   
   IF (.NOT. ALLOCATED(ppt_mo) ) &
        ALLOCATE (ppt_mo(1,12))
   
   IF (.NOT. ALLOCATED(solrad) ) &
        ALLOCATE (solrad(columns, rows,12), stat=ios)
   
   IF (.NOT. ALLOCATED(NDVI) ) &
        ALLOCATE (NDVI(columns, rows,12), stat=ios)
   
   IF (.NOT. ALLOCATED(BF) ) &
        ALLOCATE (BF(columns, rows,12), stat=ios)

   IF (.NOT. ALLOCATED(soiltext) ) &
        ALLOCATE (soiltext(columns, rows), stat=ios)
   
   IF (.NOT. ALLOCATED(veg) ) &
        ALLOCATE (veg(columns, rows), stat=ios)
   
   IF (.NOT. ALLOCATED(fuelneed) ) &
        ALLOCATE (fuelneed(columns, rows), stat=ios)
   
   IF (.NOT. ALLOCATED(popdens) ) &
        ALLOCATE (popdens(columns, rows), stat=ios)
   

   !<<<<<<<<READ IN DATA>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   print*, 'Reading in data arrays'
   filename(1:f_len)=filepath
   filename(f_len+1:f_len+8)='mon_.csv'
   open(unit=2, file=filename, form='FORMATTED', iostat=ios)
   read(2,*) months
   close(2)

   filename(1:f_len)=filepath
   filename(f_len+1:f_len+8)='perc_tre'
   open(unit=2, file=filename, form='FORMATTED', iostat=ios)
   IF (IOS .ne. 0) THEN
      print '(a)', 'There was an error reading the percent treecover file...&
                    &aborting run'
      stop
   ELSE
      read(2, *) perc_tree
      close(2)
   END IF

   filename(f_len+1:f_len+8)='perc_her'
   open(unit=2, file=filename, form='FORMATTED', iostat=ios)
   IF (IOS /= 0) THEN
      print '(a)', 'There was an error reading the percent herbaceous file...&
                    &aborting run'
      stop
   ELSE
      read(2, *) perc_herb
      close(2)
   END IF

   CALL CasaRegridInit
   IF ( .NOT. (LCPLE) ) THEN
      DO i=1,12
         filename(f_len+1:f_len+5)='airtg'
         filename(f_len+6:f_len+8)=months(i)
         open(unit=2, file=filename, form='FORMATTED', iostat=ios)
         IF (IOS /= 0) THEN
            print '(a)', 'There was an error reading the air temp file... &
                 &aborting run'
            stop
         ENDIF
      
         read(2,*) geos
         close(2)
      
         CALL regrid4x5to1x1(1, geos, geos1x1)
         airt(:,:,i)=geos1x1

         geos(:,:)=0.0e+0_fp
         geos1x1(:,:)=0.0e+0_fp

         airt(:,:,i)=airt(:,:,i)-273.15
      
         filename(f_len+1:f_len+5)='pptg_'
         filename(f_len+6:f_len+8)=months(i)
         open(unit=2, file=filename, form='FORMATTED', iostat=ios)
         IF (IOS /= 0) THEN
            print '(a)', 'There was an error reading the precipitation file...&
                 &aborting run'
            stop
         ENDIF
         
         read(2,*) geos
         close(2)
         
         CALL regrid4x5to1x1(1, geos, geos1x1)
         ppt(:,:,i)=geos1x1
         geos(:,:)=0.0e+0_fp
         geos1x1(:,:)=0.0e+0_fp
   
         ppt_mo(1,i)=sum(ppt(:,:,i))

         filename(f_len+1:f_len+5)='radg_'
         filename(f_len+6:f_len+8)=months(i)
         open(unit=2, file=filename, form='FORMATTED', iostat=ios)
         IF (IOS /= 0) THEN
            print '(a)', 'There was an error reading the solar radiation file&
                 &...aborting run'
            stop
         ENDIF

         read(2,*) geos
         close(2)

         CALL regrid4x5to1x1(1, geos, geos1x1)
         solrad(:,:,i)=geos1x1
         geos(:,:)=0.0e+0_fp
         geos1x1(:,:)=0.0e+0_fp
      END DO
   ENDIF

   DO i=1,12
      filename(f_len+1:f_len+5)='NDVI_'
      filename(f_len+6:f_len+8)=months(i)
      open(unit=2, file=filename, form='FORMATTED', iostat=ios)
      IF (IOS /= 0) THEN
         print '(a)', 'There was an error reading the NDVI file...&
                       &aborting run'
         stop
      ELSE
         read(2,*) NDVI(:,:,i)
         close(2)
      END IF
   END DO
 
   DO i=1,12
      filename(f_len+1:f_len+5)='burn_'
      filename(f_len+6:f_len+8)=months(i)
      open(unit=2, file=filename, form='FORMATTED', iostat=ios)
      IF (IOS /= 0) THEN
         print '(a)', 'There was an error reading the fraction burned file&
                        &...aborting run'
         stop
      ELSE
         read(2,*) BF(:,:,i)
         close(2)
      END IF
   END DO

   filename(f_len+1:f_len+8)='soiltext'
   open(unit=2, file=filename, form='FORMATTED', iostat=ios)
   IF (IOS /= 0) THEN
      print '(a)', 'There was an error reading the soil texture file...&
                    &aborting run'
      stop
   ELSE
      read(2,*) soiltext
      close(2)
   END IF
   
   filename(f_len+1:f_len+8)='vegetati'
   open(unit=2, file=filename, form='FORMATTED', iostat=ios)
   IF (IOS /= 0) THEN
      print '(a)', 'There was an error reading the vegetation file...&
                    &aborting run'
      stop
   ELSE
      read(2,*) veg
      close(2)
   END IF

   filename(f_len+1:f_len+8)='fuelneed'
   open(unit=2, file=filename, form='FORMATTED', iostat=ios)
   IF (IOS /= 0) THEN
      print '(a)', 'There was an error reading the needed fuel file...&
                    &aborting run'
      stop
   ELSE
      read(2,*) fuelneed
      close(2)
   END IF
   
   filename(f_len+1:f_len+8)='popdenst'
   open(unit=2, file=filename, form='FORMATTED', iostat=ios)
   IF (IOS /= 0) THEN
      print '(a)', 'There was an error reading the population density file...&
                    &aborting run'
      stop
   ELSE
      read(2,*) popdens
      close(2)
   END IF
   print *, 'Finished reading in data...now resizing arrays'
   print *, ''

 END SUBROUTINE load_data
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: conv_to_1d
!
! !DESCRIPTION: Convert maps into one column with only vegetated grid cells
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CONV_TO_1D
!
! !USES:
!
   USE defineConstants
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   integer :: i
   character(len=f_len+8)                :: filename
   
   real(fp), dimension(columns, rows)      ::   mask1, dummy

   real(fp)  :: radius=6378140.000e+0_fp !m at equator
   real(fp)  :: pi=3.14159265e+0_fp
   real(fp)  :: g, a, apixel
   real(fp), dimension(columns, rows)  :: gridAreac
   real(fp)  :: testa(180)

   IF (.NOT. ALLOCATED(mask2) ) ALLOCATE(mask2(columns, rows))
   mask1=(perc_tree+perc_herb)
   mask2=(mask1*veg)
!<<<<<<<<<<<<REGRID DATA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


   IF (.NOT. ALLOCATED(veg1      ) ) ALLOCATE(veg1(n_veg,1))
   IF (.NOT. ALLOCATED(soiltext1 ) ) ALLOCATE(soiltext1(n_veg,1))
   IF (.NOT. ALLOCATED(fuelneed1 ) ) ALLOCATE(fuelneed1(n_veg,1))
   IF (.NOT. ALLOCATED(popdens1  ) ) ALLOCATE(popdens1(n_veg,1))
   IF (.NOT. ALLOCATED(perc_tree1) ) ALLOCATE(perc_tree1(n_veg,1))
   IF (.NOT. ALLOCATED(perc_herb1) ) ALLOCATE(perc_herb1(n_veg,1))
   IF (.NOT. ALLOCATED(frac_tree ) ) ALLOCATE(frac_tree(n_veg, 1))
   IF (.NOT. ALLOCATED(frac_herb ) ) ALLOCATE(frac_herb(n_veg, 1))
   IF (.NOT. ALLOCATED(frac_veg  ) ) ALLOCATE(frac_veg(n_veg, 1))
   IF (.NOT. ALLOCATED(airt1     ) ) ALLOCATE(airt1(n_veg, 12))
   IF (.NOT. ALLOCATED(ppt1      ) ) ALLOCATE(ppt1(n_veg, 12))
   IF (.NOT. ALLOCATED(solrad1   ) ) ALLOCATE(solrad1(n_veg, 12))
   IF (.NOT. ALLOCATED(NDVI1     ) ) ALLOCATE(NDVI1(n_veg, 12))
   IF (.NOT. ALLOCATED(BF1       ) ) ALLOCATE(BF1(n_veg, 12))
   IF (.NOT. ALLOCATED(gridAreaa ) ) ALLOCATE(gridAreaa(columns, rows))
   IF (.NOT. ALLOCATED(gridAreab ) ) ALLOCATE(gridAreab(n_veg,1))
   IF (.NOT. ALLOCATED(maxt      ) ) ALLOCATE(maxt(n_veg, 1))
   IF (.NOT. ALLOCATED(mint      ) ) ALLOCATE(mint(n_veg, 1))

   veg1            =maskfile(veg, mask2)
   soiltext1       =maskfile(soiltext, mask2)
   fuelneed1       =maskfile(fuelneed, mask2)
   popdens1        =maskfile(popdens, mask2)
   perc_tree1      =maskfile(perc_tree, mask2)
   perc_herb1      =maskfile(perc_herb, mask2)
   airt1           =mask12file(airt, mask2)
   ppt1            =mask12file(ppt, mask2)
   solrad1         =mask12file(solrad, mask2)
   NDVI1           =mask12file(NDVI, mask2)
   BF1             =mask12file(BF, mask2)
   !BF1             =BF1*1.10
   frac_veg(:,1)=perc_tree1(:,1)+perc_herb1(:,1)
   
   DO i=1, n_veg
      IF (frac_veg(i,1) .gt. 0e+0_fp) THEN
         frac_tree(i,1)=perc_tree1(i,1)/frac_veg(i,1)
         frac_herb(i,1)=perc_herb1(i,1)/frac_veg(i,1)
      ELSE
         frac_tree(i,1)=0.0e+0_fp
         frac_herb(i,1)=0.0e+0_fp
      ENDIF
   END DO
   
   DO i=1, n_veg
      maxt(i,1)=maxval(airt1(i,:))
      mint(i,1)=minval(airt1(i,:))
      !  DO j=1, 12
      !  IF (BF1(i,j) .gt. 1.0) THEN
      !      BF1(i,j)=1.0
      !  ENDIF
      !  END DO
   END DO
   !makes a grid area map depending on the resolution 
   
   g=0.0e+0_fp
   a=90.0e+0_fp
   DO i=1,rows
      apixel=2.00e+0_fp*pi*radius
      apixel=apixel/columns
      apixel=apixel*apixel
      g=a*0.0174532925e+0_fp
      testa(i)=cos(g)
      gridAreac(:,i)=apixel*abs(testa(i))
      a=a-1!(180.0/rows)
   END DO
   print *, sum(gridAreac)
   gridAreab       =maskfile(gridAreac, mask2)
   gridAreaa       =gridAreac
   
   
   filename(1:f_len)=filepath
   filename(f_len+1:f_len+8)='grida1x1'
   OPEN(UNIT=4, file=filename, FORM='FORMATTED')
   WRITE(4,FMT="(360E12.5)") gridAreaa
   CLOSE(4)
   
   filename(f_len+1:f_len+8)='gridarea'
   OPEN(UNIT=4, file=filename, FORM='FORMATTED')
   WRITE(4,FMT="(1E12.5)") gridAreab
   CLOSE(4)

   filename(f_len+1:f_len+8)='mask_veg'
   OPEN(UNIT=4, file=filename, FORM='FORMATTED')
   WRITE(4,FMT="(1E12.5)") mask2
   CLOSE(4)
   
 END SUBROUTINE CONV_TO_1D
!EOC 
!-----------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: maskfile
!
! !DESCRIPTION: This subroutine ...
!\\
!\\
! !INTERFACE:
!
 function maskfile(dummy,mask3) result (masked_file)
!
! !USES:
!
   USE defineConstants

   implicit none
!
! !INPUT PARAMETERS:
!
   real(fp), dimension(columns, rows) :: dummy, mask3
!
! !RETURN VALUE:
!
   real(fp), dimension(n_veg,1)  :: masked_file
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   integer                 :: a,b,c,g,i,j
   character(len=f_len_output+3)       :: filename3
   real(fp), dimension(3,n_veg)           :: key
   
   filename3(1:f_len_output)=outputpath
   filename3(f_len_output+1:f_len_output+3)='key'
   
   g=1
   DO i=1,columns
      DO j=1,rows
         IF (mask3(i,j) > 0) THEN
            masked_file(g,1)=dummy(i,j)
            key(1,g)=i
            key(2,g)=j
            key(3,g)=veg(i,j)
            g=g+1
         END IF
      END DO
   END DO
   
   OPEN(UNIT=4, file=filename3, FORM='FORMATTED')
   WRITE(4,FMT="(3F9.2)") key
   CLOSE(4)
 end function maskfile
!EOC 
!----------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: mask12file
!
! !DESCRIPTION: This subroutine ...
!\\
!\\
! !INTERFACE:
!
 function mask12file(dummy12, mask3) result (masked_12file)
!
! !USES:
!
   USE defineConstants

   implicit none
!
! !INPUT PARAMETERS:
!
   real(fp), dimension(columns, rows, 12)  ::  dummy12
   real(fp), dimension(columns, rows)      ::  mask3
!
! !RETURN VALUE:
!
   real(fp), dimension(n_veg, 12)              ::  masked_12file
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   integer                               ::  a, g, i, j, k

   g=1
   DO i=1,columns
      DO j=1, rows
         IF (mask3(i,j) .gt. 0e+0_fp) THEN
            DO k=1,12  ! months
               masked_12file(g,k)=dummy12(i,j,k)
            END DO
            g=g+1
         END IF
      END DO
   END DO
   
 end function mask12file
!EOC
END MODULE loadCASAinput


