
!******************************************************************************************************
!>
!  A better `get_command_argument` for modern Fortran.

    module argv_module

    implicit none

    private

    public :: argv

    contains
!******************************************************************************************************

!******************************************************************************************************
!>
!  Returns the specified command line argument in an allocatable string.

    function argv(arg_num) result(arg)

    integer,intent(in) :: arg_num !! argument number (>=0)
    character(len=:),allocatable :: arg !! the command line argument

    integer :: ilen  !! command-line argument string length

    call get_command_argument(arg_num, length=ilen)
    arg = repeat(' ',ilen)
    call get_command_argument(arg_num, value=arg)

    end function argv
!******************************************************************************************************

!******************************************************************************************************
    end module argv_module
!******************************************************************************************************