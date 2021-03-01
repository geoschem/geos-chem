!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: m_do_err_out.F90
!
! !INTERFACE:
!
module m_Do_Err_Out
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public  Do_Err_Out
!
! !DESCRIPTION: Provides a routine to print an error message and exit the code.
!\\
!\\
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  10 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  10 Jul 2014 - R. Yantosca - Cosmetic changes to ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Err_Out
!
! !INTERFACE:
!
  subroutine Do_Err_Out  &
       (err_msg, err_do_stop, err_num_ints, err_int1, err_int2,  &
       err_num_reals, err_real1, err_real2)
!
    implicit none
!
! !INPUT PARAMETERS:
!!     err_msg       : error message to be printed out
!!     err_do_stop   : do stop on error?
!!     err_num_ints  : number of integers to be printed out (0, 1, or 2)
!!     err_int1      : integer 1 to print out
!!     err_int2      : integer 2 to print out
!!     err_num_reals : number of reals to be printed out (0, 1, or 2)
!!     err_real1     : real 1 to print out
!!     err_real2     : real 2 to print out
    character (len=*), intent(in) :: err_msg
    logical          , intent(in) :: err_do_stop
    integer          , intent(in) :: err_num_ints
    integer          , intent(in) :: err_int1
    integer          , intent(in) :: err_int2
    integer          , intent(in) :: err_num_reals
    real*8           , intent(in) :: err_real1
    real*8           , intent(in) :: err_real2
!
! !DESCRIPTION: Outputs error messages, and exits if requested.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REMARKS:
!  NOTE: SHOULD PROPAGATE ERROR CODE TO MAIN PROGRAM LEVEL!
!
! !REVISION HISTORY:
!  Initial code.
!  07 Mar 2017 - R. Yantosca - Now exit with error code 999, to avoid an
!                              inadvertent error code of 0 being returned
!  08 Nov 2017 - R. Yantosca - Now flush the buffer after writing,
!                              to be visilble after stop (esp. w/ gfortran)
!EOP
!-------------------------------------------------------------------------
!BOC

    ! Write separator
    WRITE( 6, '(/,a,/)' ) REPEAT( '!', 79 )

    ! Write error message
    WRITE( 6,'(a)' ) TRIM( err_msg )

    ! Write error codes
    IF ( err_num_ints == 1 ) THEN
       WRITE( 6,'(i10)'   ) err_int1
    ELSE IF ( err_num_ints == 2 ) then
       WRITE( 6, '(2i10)' ) err_int1, err_int2
    ENDIF

    IF ( err_num_reals == 1 ) THEN
       WRITE( 6, '(f13.6  )' ) err_real1
    ELSE IF ( err_num_reals == 2 ) THEN
       WRITE( 6, '(2f13.6 )' ) err_real1, err_real2
    ENDIF

    ! Write separator
    WRITE( 6, '(/,a,/)' ) REPEAT( '!', 79 )

    ! Flush the buffer
    CALL Flush( 6 )

    ! Stop with error (if requested)
    ! NOTE: We should pass back the error code to the main routine
    IF ( err_do_stop ) THEN
        WRITE( 6, '(a,/)' ) 'Code stopped from DO_ERR_OUT '               // &
                            '(in module NcdfUtil/m_do_err_out.F90) '
        WRITE( 6, '(a)'   ) 'This is an error that was encountered '      // &
                            'in one of the netCDF I/O modules,'
        WRITE( 6, '(a)'   ) 'which indicates an error in writing to '     // &
                            'or reading from a netCDF file!'

        ! Write separator
        WRITE( 6, '(/,a,/)' ) REPEAT( '!', 79 )

        ! Flush stdout buffer
        CALL Flush( 6 )

        ! NOTE: Should not exit but pass error code up
        ! work on this for a future version
        CALL Exit( 999 )
    ENDIF

    RETURN

  end subroutine Do_Err_Out
!EOC
end module m_Do_Err_Out
