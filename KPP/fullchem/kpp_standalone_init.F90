module kpp_standalone_init
  implicit none
  public
contains

subroutine read_input(filename, R, C, SPC_NAMES, Hstart, Hexit, cosSZA, level, fileTotSteps, OperatorTimestep)
USE gckpp_Parameters

  IMPLICIT NONE

  real(dp), intent(out) :: C(NSPEC)
  real(dp), intent(out) :: R(NREACT)
  real(dp), intent(out) :: Hstart
  real(dp), intent(out) :: Hexit
  real(dp), intent(out) :: cosSZA
  real(dp), intent(out) :: OperatorTimestep
  integer, intent(out)  :: level
  integer, intent(out)  :: fileTotSteps
  integer :: SPC_MAP(NSPEC)



  character(len=*), intent(in) :: SPC_NAMES(NSPEC)
  character(len=*), intent(in) :: filename
  integer :: i, ierr, NHEADER, idx
  character(200) :: line
  logical :: existbool

  ! Declare variables for file I/O
  integer :: file_unit

  ! Open the file for reading
  file_unit = 999
  inquire(file=filename, exist=existbool)
  if (existbool .neqv. .TRUE.) then
     print *, "Error: input file does not exist: ", trim(filename)
     stop
  end if
  open(unit=file_unit, file=filename, iostat=ierr)
  if (ierr /= 0) then
     print *, "Error opening input file"
     stop
  end if

  ! Read the number of header lines
  read(file_unit, *) NHEADER

  ! Read the header lines
  do i = 1, NHEADER
     read(file_unit, '(A)', iostat=ierr) line
     if (ierr /= 0)  then
      print *, "Error reading line", i 
      exit
     end if
     ! Get level
     if (index(line, 'GEOS-Chem Vertical Level:') > 0 ) then
        idx = index(line, ':') + 1
        read(line(idx:), *) level
     endif
     ! Get cosSZA
     if (index(line, 'Cosine of solar zenith angle:') > 0 ) then
      idx = index(line, ':') + 1
      read(line(idx:), *) cosSZA
     end if
     ! get Hstart
     if (index(line, 'Init KPP Timestep (seconds):') > 0 ) then
      idx = index(line, ':') + 1
      read(line(idx:), *) Hstart
     end if

     ! get Hexit
     if (index(line, 'Exit KPP Timestep (seconds):') > 0 ) then
      idx = index(line, ':') + 1
      read(line(idx:), *) Hexit
     end if

     ! get fileTotSteps
     if (index(line, 'Number of internal timesteps:') > 0 ) then
      idx = index(line, ':') + 1
      read(line(idx:), *) fileTotSteps
     end if

     ! Get value of operator splitting timestep
      if (index(line, 'Chemistry operator timestep (seconds):') > 0 ) then
        idx = index(line, ':') + 1
        read(line(idx:), *) OperatorTimestep
      end if
  end do

  ! Read the species and their concentrations
  do i = 1, NSPEC
     read(file_unit, '(A)', iostat=ierr) line
     if (ierr /= 0)  then
      print *, "Error reading line", i+NHEADER 
      exit
     end if
     idx = index(line, ',') + 1
     read(line(idx:), *) C(i)
    !  Check if the species name matches the expected SPC_NAMES(i)
      if (trim(line(1:idx-2)) /= trim(SPC_NAMES(i))) then
          print *, "Error: species name mismatch"
          print *, "Expected: ", SPC_NAMES(i)
          print *, "Found: ", line(1:idx-2)
          stop
      end if
  end do
  ! Read the rate constants
  do i = 1, NREACT
     read(file_unit, '(A)', iostat=ierr) line
     if (ierr /= 0)  then
      print *, "Error reading line", i+NSPEC+NHEADER 
      exit
     end if
     idx = index(line, ',') + 1
     read(line(idx:), *) R(i)
  end do

  ! Close the file
  close(file_unit)
end subroutine read_input

end module kpp_standalone_init

