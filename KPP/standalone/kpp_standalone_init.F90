module kpp_standalone_init

  implicit none
  public

contains

subroutine read_input( filename,     R,                C,      SPC_NAMES,    &
                       Hstart,       Hexit,            cosSZA, level,        &
                       fileTotSteps, OperatorTimestep, ICNTRL, RCNTRL,       &
                       ATOL                                                 )

  USE gckpp_Parameters

  IMPLICIT NONE

  ! Inputs
  character(len=*), intent(in)  :: filename          ! Input file name
  character(len=*), intent(in)  :: SPC_NAMES(NSPEC)  ! Chemical species names

  ! Outputs
  real(dp),         intent(out) :: C(NSPEC)          ! Species conc (molec/cm3)
  real(dp),         intent(out) :: R(NREACT)         ! Rates (molec/cm3/s)
  real(dp),         intent(out) :: Hstart            ! Init KPP timestep (s)
  real(dp),         intent(out) :: Hexit             ! Final KPP timestep (s)
  real(dp),         intent(out) :: cosSZA            ! COS( solar zen angle )
  real(dp),         intent(out) :: OperatorTimestep  ! External timestep
  real(dp),         intent(out) :: RCNTRL(20)        ! Integrator options
  real(dp),         intent(out) :: ATOL(NVAR)        ! Abs. tolerance
  integer,          intent(out) :: level             ! Model level
  integer,          intent(out) :: fileTotSteps      ! Total integration steps
  integer,          intent(out) :: ICNTRL(20)        ! Integrator options

  ! Local variables
  integer                       :: SPC_MAP(NSPEC)
  integer                       :: i, ierr, NHEADER, idx, jdx
  integer                       :: file_unit
  integer                       :: i1, i2
  integer                       :: r1, r2
  logical                       :: existbool
  logical                       :: parse_icntrl
  logical                       :: parse_rcntrl
  character(len=255)            :: line

  ! Initialize outputs for safety's sake
  ATOL             = -1.0_dp   ! Values not specified will get a default value
  C                =  0.0_dp
  cosSZA           =  0.0_dp
  fileTotSteps     =  0
  Hexit            =  0.0_dp
  Hstart           =  0.0_dp
  ICNTRL           =  0
  level            =  0
  OperatorTimestep =  0.0_dp
  R                =  0.0_dp
  RCNTRL           =  0.0_dp

  ! Set ICNTRL options to Rosenbrock inputs by default
  ! These will be overwritten if ICNTRL is found in the input file
  ICNTRL(1)        =  1
  ICNTRL(3)        =  4
  ICNTRL(7)        =  1
  ICNTRL(15)       = -1

  ! For reading ICNTRL and RCNTRL
  parse_icntrl     = .false.
  parse_rcntrl     = .false.
  i1               = 1
  i2               = 10
  r1               = 1
  r2               = 5

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

     ! Get ICNTRL integrator options
     ! Read 10 integer values (width=6) from the next 2 lines
     if (index(line, 'ICNTRL integrator options used:') > 0 ) then
        parse_icntrl = .true.
        cycle
     end if
     if ( parse_icntrl ) then
        read( line, '(10i6)' ) ICNTRL(i1:i2)
        i1 = i1 + 10
        i2 = i2 + 10
        if ( i1 > 20 ) then
           parse_icntrl = .false.
           cycle
        end if
     end if

     ! Get RCNTRL integrator options
     ! Read 5 real values (width=13.6) from the next 2 lines
     if (index(line, 'RCNTRL integrator options used:') > 0 ) then
        parse_rcntrl = .true.
        cycle
     end if
     if ( parse_rcntrl ) then
        read( line, '(5F13.6)' ) RCNTRL(r1:r2)
        r1 = r1 + 5
        r2 = r2 + 5
        if ( r1 > 20 ) then
           parse_rcntrl = .false.
           cycle
        end if
     end if

  end do


  ! Read the species and their concentrations
  do i = 1, NSPEC
     read(file_unit, '(A)', iostat=ierr) line
     if (ierr /= 0)  then
        print *, "Error reading line", i+NHEADER
        exit
     end if

     ! Read species concentration (molec/cm3)
     idx = index( line, ',' ) + 1
     read( line(idx:), * ) C(i)

     ! If present, also read absolute tolerance
     jdx = index( line(idx+1:), ',' ) + idx + 1
     if ( jdx > idx + 1 )  read( line(jdx:), * ) ATOL(i)

     ! Check if the species name matches the expected SPC_NAMES(i)
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
