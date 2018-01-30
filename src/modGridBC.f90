module modGridBC

	use modSpecies
	use modMesh

	implicit none

	abstract interface
		subroutine adjustGrid(ng,g,frac)
			use constants
			integer, intent(in) :: ng
			integer, intent(inout) :: g(:)
			real(mp), intent(inout) :: frac(:)
		end subroutine
	end interface

contains

!======================= AdjustGrid for BC =============================================

	subroutine adjustGrid_periodic(ng,g,frac)
		integer, intent(in) :: ng
		integer, intent(inout) :: g(:)
		real(mp), intent(inout) :: frac(:)

			if( g(1)<1 ) then
				g(1) = g(1) + ng
			elseif( g(1)>ng ) then
				g(1) = g(1) - ng
			end if
			if( g(2)<1 ) then
				g(2) = g(2) + ng
			elseif( g(2)>ng ) then
				g(2) = g(2) - ng
			end if

		if( MINVAL(g)<1 .or. MAXVAL(g)>ng ) then
			print *, MINVAL(g), MAXVAL(g)
			print *, 'Boundary handling is failed. particle is way outside BC. stopped time stepping.'
			stop
		end if
	end subroutine

	subroutine adjustGrid_absorbing(ng,g,frac)
		integer, intent(in) :: ng
		integer, intent(inout) :: g(:)
		real(mp), intent(inout) :: frac(:)

		if( MINVAL(g)<1 .or. MAXVAL(g)>ng ) then
			print *, MINVAL(g), MAXVAL(g)
			print *, 'Boundary handling is failed. particle is way outside BC. stopped time stepping.'
			stop
		end if

		!adjustment for boundary : charge/(dx/2)
		if( g(1).eq.1 ) then
			frac(1) = frac(1)*2.0_mp
		end if
		if( g(2).eq.ng ) then
			frac(2) = frac(2)*2.0_mp
		end if
	end subroutine

end module
