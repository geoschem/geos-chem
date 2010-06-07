MODULE loadCASAinput

! Loads input files.  For global studies, these files should be given 
! as 180 by 360 matrix for 1x1 degree, 360 by 720 for half by half 
! etc.  This routine will construct a 'mask' with vegetated gridcells
! and will reshape them into an X by 1 matrix, where X is the total
! number of vegetated grid cells.  
! The files which contain monthly varying parameters (i.e. precip)
! should be given as (for 1x1) 180x360x12, and will be reshaped to
! an X by 12 matrix, one column for each month
!


USE defineConstants
USE CasaRegridModule

implicit none

!CONTINUOUS FIELD MAPS
   REAL*8, ALLOCATABLE :: perc_tree(:,:)
   REAL*8, ALLOCATABLE :: perc_herb(:,:) 

!RESIZED CONTINUOUS FIELD MAPS
   REAL*8, ALLOCATABLE :: perc_tree1(:,:)
   REAL*8, ALLOCATABLE :: perc_herb1(:,:)
   REAL*8, ALLOCATABLE :: frac_tree(:,:)
   REAL*8, ALLOCATABLE :: frac_herb(:,:)
   REAL*8, ALLOCATABLE :: frac_veg(:,:)

!CLIMATE FILES
   REAL*8, ALLOCATABLE :: airt(:,:,:)   !monthly air temperature
   REAL*8, ALLOCATABLE :: ppt(:,:,:)    !monthly precipitation
   REAL*8, ALLOCATABLE :: solrad(:,:,:) !monthly solar radiation, 
                          !Taken from Bishop and Rossow
                          !average 83-99 from Jim Collatz

   REAL*8, ALLOCATABLE :: NDVI(:,:,:)   !monthly fraction PAR
   REAL*8, ALLOCATABLE :: BF(:,:,:)     !fraction gridcell tha burns
   REAL*8, ALLOCATABLE :: ppt_mo(:,:)   !sum of all precip in each mo
   
!RESIZED CLIMATE FILES
   REAL*8, ALLOCATABLE :: airt1(:,:)
   REAL*8, ALLOCATABLE :: ppt1(:,:)
   REAL*8, ALLOCATABLE :: solrad1(:,:)
   REAL*8, ALLOCATABLE :: NDVI1(:,:)
   REAL*8, ALLOCATABLE :: BF1(:,:)
   REAL*8, ALLOCATABLE :: maxt(:,:)
   REAL*8, ALLOCATABLE :: mint(:,:)
 
!OTHER FILES
   REAL*8, ALLOCATABLE :: soiltext(:,:)  ! soil type
   REAL*8, ALLOCATABLE :: veg(:,:)       ! vegetation map
   REAL*8, ALLOCATABLE :: fuelneed(:,:)  ! fuelwood needed 
                                         !   per capita
   REAL*8, ALLOCATABLE :: popdens(:,:)   ! population density (/m2)
   REAL*8, ALLOCATABLE :: gridAreaa(:,:)
   REAL*8, ALLOCATABLE :: gridAreab(:,:)
   
!RESIZED OTHER FILES
   REAL*8, ALLOCATABLE :: soiltext1(:,:)
   REAL*8, ALLOCATABLE :: veg1(:,:)
   REAL*8, ALLOCATABLE :: fuelneed1(:,:)
   REAL*8, ALLOCATABLE :: popdens1(:,:)

   
   REAL*8, ALLOCATABLE :: mask2(:,:)
   
!<<<<<<<<<END DECLARE VARIABLES>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

CONTAINS

!-----------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: read\_CASA\_output
!
! !DESCRIPTION: Subroutine read\_CASA\_output reads data saved previously from
! an equilibrium run of GTMM. (ccc, 10/27/09)
!
! !INTERFACE:
!
  SUBROUTINE read_CASA_output()
!
! !USES:
!
    USE defineConstants
    USE defineArrays
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
   CHARACTER(len=250) :: FILENAME
   
   ! READ LAI
   FILENAME = OUTPUTPATH // 'flai'
   OPEN(UNIT=4, file=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4,FMT="(12F9.2)") LAI
   CLOSE(4)

   ! READ herb_seasonality, grass_herbivory, trees_herbivory
   FILENAME = OUTPUTPATH // 'fhsn'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, *) herb_seasonality
   CLOSE(4)

   FILENAME = OUTPUTPATH // 'fgrh'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, *) grass_herbivory
   CLOSE(4)

   FILENAME = OUTPUTPATH // 'ftrh'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, *) trees_herbivory
   CLOSE(4)

   ! READ LTCON, LTVARSUM, LAI, rootlitscalar, litterscalar, hlitterscalar
   FILENAME = OUTPUTPATH // 'fltc'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, *) LTCON
   CLOSE(4)

   FILENAME = OUTPUTPATH // 'flvr'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, *) LTVARSUM
   CLOSE(4)

   FILENAME = OUTPUTPATH // 'fl_i'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, FMT="(12F9.2)") LAI
   CLOSE(4)
   
   FILENAME = OUTPUTPATH // 'frls'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, FMT="(12F9.2)") rootlitscalar
   CLOSE(4)
   
   FILENAME = OUTPUTPATH // 'f_ls'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, FMT="(12F9.2)") litterscalar
   CLOSE(4)
   
   FILENAME = OUTPUTPATH // 'fhls'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, FMT="(12F9.2)") hlitterscalar
   CLOSE(4)
   
   FILENAME = OUTPUTPATH // 'fali'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, *) AVELAI
   CLOSE(4)
   
   ! READ NPP, abiotic
   FILENAME = OUTPUTPATH // 'NPP'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, *) NPP
   CLOSE(4)
   
   FILENAME = OUTPUTPATH // 'ABIO'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, *) abiotic
   CLOSE(4)
   
   ! READ lais, topt
   FILENAME = OUTPUTPATH // 'flis'
   OPEN(UNIT=4, file=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4,FMT="(13F9.2)") lais
   CLOSE(4)
   
   FILENAME = OUTPUTPATH // 'ftop'
   OPEN(UNIT=4, file=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4,*) topt
   CLOSE(4)
   
   ! READ ccWood, ccLeaf, ccFineLitter, ccCwd, mortality_tree
   FILENAME = OUTPUTPATH // 'fccw'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, FMT="(12F9.2)") ccWood
   CLOSE(4)

   FILENAME = OUTPUTPATH // 'fccl'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, FMT="(12F9.2)") ccLeaf
   CLOSE(4)

   FILENAME = OUTPUTPATH // 'fcfl'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, FMT="(12F9.2)") ccFineLitter
   CLOSE(4)

   FILENAME = OUTPUTPATH // 'fcwd'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, FMT="(12F9.2)") ccCwd
   CLOSE(4)

   FILENAME = OUTPUTPATH // 'fmor'
   OPEN(UNIT=4, FILE=FILENAME, STATUS="OLD", FORM="FORMATTED")
   READ(4, *) mortality_tree
   CLOSE(4)
   

  END SUBROUTINE read_CASA_output
!EOC
!------------------------------------------------------------------------------

  SUBROUTINE load_data(year, month, LCPLE, TS, PREACC, RADSWG)

   USE defineConstants

   ! ARGUMENTS
   INTEGER, INTENT(IN)  :: year, month
   LOGICAL, INTENT(IN)  :: LCPLE

   REAL*8, INTENT(IN), OPTIONAL, DIMENSION(72, 46)  :: TS, &
        PREACC, RADSWG

   ! OTHER VARIABLES
   
   integer   ::   i, j, k, l, m, n, b, c, ios
   character(len=f_len+8)                :: filename
   real*8, dimension(columns, rows)      ::   mask1, dummy
   real*8, dimension(columns, rows, 12)  ::   dummy12
   real*8, dimension(72,46)              ::   geos
   real*8, dimension(columns, rows)      ::   geos1x1
   character(3), dimension(12)           ::   months     ! Months names
   
   real*8  :: radius=6378140.000d0 !m at equator
   real*8  :: pi=3.14159265d0
   real*8  :: g, a, apixel
   real*8, dimension(columns, rows)  :: gridAreac
   real*8  :: testa(180)

   !<<<<<<<<READ IN DATA>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   print*, 'Reading in data arrays'
   filename(1:f_len)=filepath
   filename(f_len+1:f_len+8)='mon_.csv'
   open(unit=2, file=filename, form='FORMATTED', iostat=ios)
   read(2,*) months
   close(2)

   IF (.NOT. ALLOCATED(perc_tree) ) &
        ALLOCATE (perc_tree(columns, rows), stat=ios)
   
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

   IF (.NOT. ALLOCATED(perc_herb) ) &
        ALLOCATE (perc_herb(columns, rows), stat=ios)
   
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

   IF (.NOT. ALLOCATED(airt) ) &
        ALLOCATE (airt(columns, rows,12), stat=ios)
   
   CALL CasaRegridInit
!ccc   DO i=1,12
   IF ( LCPLE ) THEN
      CALL regrid4x5to1x1(1, TS, geos1x1)
   ELSE
      filename(f_len+1:f_len+5)='airtg'
      !ccc   filename(f_len+6:f_len+8)=months(i)
      filename(f_len+6:f_len+8)=months(month)
      open(unit=2, file=filename, form='FORMATTED', iostat=ios)
      IF (IOS /= 0) THEN
         print '(a)', 'There was an error reading the air temp file... &
              &aborting run'
         stop
      ENDIF
      
      read(2,*) geos
      close(2)
      
      CALL regrid4x5to1x1(1, geos, geos1x1)
   ENDIF
   
   !ccc   airt(:,:,i)=geos1x1
   airt(:,:,month)=geos1x1
   geos(:,:)=0.0d0
   geos1x1(:,:)=0.0d0
   !ccc END DO
   !ccc   airt(:,:,:)=airt(:,:,:)-273.15
   airt(:,:,month)=airt(:,:,month)-273.15d0
   
   IF (.NOT. ALLOCATED(ppt) ) &
        ALLOCATE (ppt(columns, rows,12), stat=ios)
   
!ccc   DO i=1,12   
   IF ( LCPLE ) THEN
      CALL regrid4x5to1x1(1, PREACC, geos1x1)
   ELSE
      filename(f_len+1:f_len+5)='pptg_'
      !ccc      filename(f_len+6:f_len+8)=months(i)
      filename(f_len+6:f_len+8)=months(month)
      open(unit=2, file=filename, form='FORMATTED', iostat=ios)
      IF (IOS /= 0) THEN
         print '(a)', 'There was an error reading the precipitation file...&
              &aborting run'
         stop
      ENDIF
   
      read(2,*) geos
      close(2)
      
      CALL regrid4x5to1x1(1, geos, geos1x1)
   ENDIF
!ccc      ppt(:,:,i)=geos1x1
   ppt(:,:,month)=geos1x1
   geos(:,:)=0.0d0
   geos1x1(:,:)=0.0d0
!ccc   END DO
   
   IF (.NOT. ALLOCATED(ppt_mo) ) &
        ALLOCATE (ppt_mo(1,12))
   
   DO i=1,12
      ppt_mo(1,i)=sum(ppt(:,:,i))
   END DO
      
   IF (.NOT. ALLOCATED(solrad) ) &
        ALLOCATE (solrad(columns, rows,12), stat=ios)
   
!ccc   DO i=1,12
   IF ( LCPLE ) THEN
      CALL regrid4x5to1x1(1, RADSWG, geos1x1)
   ELSE
      filename(f_len+1:f_len+5)='radg_'
!ccc      filename(f_len+6:f_len+8)=months(i)
      filename(f_len+6:f_len+8)=months(month)
      open(unit=2, file=filename, form='FORMATTED', iostat=ios)
      IF (IOS /= 0) THEN
         print '(a)', 'There was an error reading the solar radiation file&
                    &...aborting run'
         stop
      ENDIF

      read(2,*) geos
      close(2)

      CALL regrid4x5to1x1(1, geos, geos1x1)
   ENDIF
!ccc      solrad(:,:,i)=geos1x1
      solrad(:,:,month)=geos1x1
      geos(:,:)=0.0d0
      geos1x1(:,:)=0.0d0
!ccc   END DO

   IF (.NOT. ALLOCATED(NDVI) ) &
        ALLOCATE (NDVI(columns, rows,12), stat=ios)
   
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
 
   IF (.NOT. ALLOCATED(BF) ) &
        ALLOCATE (BF(columns, rows,12), stat=ios)

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

   IF (.NOT. ALLOCATED(soiltext) ) &
        ALLOCATE (soiltext(columns, rows), stat=ios)
   
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
   
   IF (.NOT. ALLOCATED(veg) ) &
        ALLOCATE (veg(columns, rows), stat=ios)
   
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

   IF (.NOT. ALLOCATED(fuelneed) ) &
        ALLOCATE (fuelneed(columns, rows), stat=ios)
   
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
   
   IF (.NOT. ALLOCATED(popdens) ) &
        ALLOCATE (popdens(columns, rows), stat=ios)
   
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
!<<<<<<<<<<END READ IN DATA>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


!<<<<<<<<<<CONVERT MAPS INTO ONE COLUMN WITH ONLY VEGETATED GRID
!CELLS>>>>>>
   IF (.NOT. ALLOCATED(mask2) ) &
    ALLOCATE(mask2(columns, rows))
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
   IF (frac_veg(i,1) .gt. 0d0) THEN
        frac_tree(i,1)=perc_tree1(i,1)/frac_veg(i,1)
        frac_herb(i,1)=perc_herb1(i,1)/frac_veg(i,1)
   ELSE
        frac_tree(i,1)=0.0d0
        frac_herb(i,1)=0.0d0
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

g=0.0d0
a=90.0d0
DO i=1,rows
   apixel=2.00d0*pi*radius
   apixel=apixel/columns
   apixel=apixel*apixel
   g=a*0.0174532925d0
   testa(i)=cos(g)
   gridAreac(:,i)=apixel*abs(testa(i))
   a=a-1!(180.0/rows)
END DO
print *, sum(gridAreac)
gridAreab       =maskfile(gridAreac, mask2)
gridAreaa       =gridAreac


filename(f_len+1:f_len+8)='grida1x1'
OPEN(UNIT=4, file=filename, FORM='FORMATTED')
WRITE(4,FMT="(360E12.5)") gridAreaa
CLOSE(4)

filename(f_len+1:f_len+8)='gridarea'
OPEN(UNIT=4, file=filename, FORM='FORMATTED')
WRITE(4,FMT="(1E12.5)") gridAreab
CLOSE(4)
END SUBROUTINE load_data

!-----------------------------------------------------------------------------
function maskfile(dummy,mask3) result (masked_file)

USE defineConstants

implicit none

real*8, dimension(columns, rows) :: dummy, mask3
real*8, dimension(n_veg,1)  :: masked_file

integer                 :: a,b,c,g,i,j
character(len=f_len_output+3)       :: filename3
real*8, dimension(3,n_veg)           :: key

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

!----------------------------------------------------------------------------
function mask12file(dummy12, mask3) result (masked_12file)

USE defineConstants

implicit none

real*8, dimension(columns, rows, 12)  ::  dummy12
real*8, dimension(n_veg, 12)              ::  masked_12file
real*8, dimension(columns, rows)      ::  mask3
integer                               ::  a, g, i, j, k

g=1
DO i=1,columns
   DO j=1, rows
      IF (mask3(i,j) .gt. 0d0) THEN
         DO k=1,12  ! months
            masked_12file(g,k)=dummy12(i,j,k)
         END DO
         g=g+1
      END IF
   END DO
END DO
          
end function mask12file

!----------------------------------------------------------------------------

END MODULE loadCASAinput
