module modQoI

	use modAdj

	implicit none

contains

!=============Mean Kinetic Energy====================

	subroutine MKE(pm,k,J)
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: k
		real(mp), intent(inout) :: J

		if( k .ge. pm%ni ) then
			J = J + 1.0_mp/pm%p(1)%np/(pm%nt-pm%ni)*SUM( pm%p(1)%vp(:,1)**2 )
		end if
	end subroutine

	subroutine dMKE(adj,pm,nk)
		type(adjoint), intent(inout) :: adj
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: nk

		if( nk .ge. pm%ni ) then
			adj%dp(1)%vp(:,1) = adj%dp(1)%vp(:,1) + 2.0_mp/pm%p(1)%np/(pm%nt-pm%ni)*pm%p(1)%vp(:,1)
		end if
	end subroutine

!============Mean Efield Energy======================

	subroutine MPE(pm,k,J)
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: k
		real(mp), intent(inout) :: J

		if( k.ge.1 .and. k.le.pm%ni ) then
			J = J + 1.0_mp/pm%ng/pm%ni*SUM( pm%m%E**2 )
		end if
	end subroutine

	subroutine dMPE(adj,pm,nk)
		type(adjoint), intent(inout) :: adj
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: nk

		if( nk.ge.1 .and. nk.le.pm%ni ) then
			adj%dm%E = adj%dm%E - 2.0_mp/pm%ng/pm%ni*pm%m%E
		end if
	end subroutine

!============Testmodule : Final timestep 1D Total Kinetic Energy

	subroutine TestQoI(pm,k,J)
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: k
		real(mp), intent(inout) :: J

		if( k.eq.pm%nt ) then
			J = J + sum( pm%p(1)%vp(:,1)**2 )
		end if
	end subroutine

	subroutine dTestQoI(adj,pm,nk)
		type(adjoint), intent(inout) :: adj
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: nk

		if( nk.eq.pm%nt ) then
			adj%dp(1)%vp(:,1) = adj%dp(1)%vp(:,1) + 2.0_mp*pm%p(1)%vp(:,1)
		end if
	end subroutine

!==============Debye-length?

	subroutine Debye(pm,k,J)
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: k
		real(mp), intent(inout) :: J

		J = J + 1.0_mp/pm%nt*SUM( pm%p(1)%spwt*(pm%p(1)%xp-0.5_mp*pm%L)**2 )
	end subroutine

end module
