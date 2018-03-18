module constants

!   use mpi

	implicit none

	integer, parameter :: mp = SELECTED_REAL_KIND(15,307)
!	integer, parameter :: mp = SELECTED_REAL_KIND(33,4931)
	real(mp), parameter :: pi = 4.0_mp*ATAN(1.0_mp)
	complex(mp), parameter :: eye = (0.0_mp,1.0_mp)
    logical :: print_pm_output = .false.

contains

end module
