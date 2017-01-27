module nullMCC

	use modPM1D
   use modRecord
	use random
	implicit none

contains



!=======================================================
!	MCC global subroutine
!=======================================================
	subroutine mcc_collision(pm,n_coll)
		type(PM1D), intent(inout) :: pm
      integer, intent(out), optional :: n_coll(7)

	end subroutine

end module
