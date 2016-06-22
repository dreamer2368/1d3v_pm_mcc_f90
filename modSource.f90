module modSource

	use modPM1D
	use random

	implicit none

contains

	subroutine Null_source(pm)
		type(PM1D), intent(inout) :: pm
	end subroutine

end module