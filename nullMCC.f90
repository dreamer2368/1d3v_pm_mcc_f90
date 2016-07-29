module nullMCC

	use modPM1D
   use modRecord
	use random
	implicit none

contains

	subroutine set_null_discharge(r)
      type(recordData), intent(inout), optional :: r

      if( present(r) ) then
         allocate(r%n_coll(1,r%nt))
			r%n_coll = 0
      end if
	end subroutine

!=======================================================
!	MCC global subroutine
!=======================================================
	subroutine mcc_collision(pm,n_coll)
		type(PM1D), intent(inout) :: pm
      integer, intent(out), optional :: n_coll(7)

	end subroutine

end module
