module modQoI

	use modAdj

	implicit none

contains

!=============Mean Kinetic Energy====================

!	subroutine MKE(this,pm,i)
!		type(adjoint), intent(inout) :: this
!		type(PM3D), intent(in) :: pm
!		integer, intent(in), optional :: i
!		integer :: input
!		if( present(i) ) then
!			input = i
!		else
!			input = 0
!		end if
!
!		if( input.eq.0 ) then
!			this%J0 = 1.0_mp/pm%n/(pm%nt-pm%ni)*SUM( pm%r%vpdata(:,:,pm%ni+1:pm%nt)**2 )
!			print *, 'J0 = ', this%J0
!		elseif( input.eq.1 ) then
!			this%J1 = 1.0_mp/pm%n/(pm%nt-pm%ni)*SUM( pm%r%vpdata(:,:,pm%ni+1:pm%nt)**2 )
!			print *, 'J1 = ', this%J1
!		end if
!	end subroutine
!
!	subroutine dMKE(adj,pm,k)
!		type(adjoint), intent(inout) :: adj
!		type(PM3D), intent(in) :: pm
!		integer, intent(in) :: k
!
!		if( k >= pm%ni+1 ) then
!			adj%dvps = adj%dvps + 2.0_mp/pm%n/(pm%nt-pm%ni)*pm%r%vpdata(:,:,k)
!		end if
!	end subroutine

!============Mean Efield Energy======================

!	subroutine MPE(this,pm,i)
!		type(adjoint), intent(inout) :: this
!		type(PM3D), intent(in) :: pm
!		integer, intent(in), optional :: i
!		integer :: input
!		if( present(i) ) then
!			input = i
!		else
!			input = 0
!		end if
!
!		if( input.eq.0 ) then
!			this%J0 = 1.0_mp/PRODUCT(pm%ng)/(pm%nt-pm%ni)*SUM( pm%r%Edata(:,:,:,:,pm%ni+1:pm%nt)**2 )
!			print *, 'J0 = ', this%J0
!		elseif( input.eq.1 ) then
!			this%J1 = 1.0_mp/PRODUCT(pm%ng)/(pm%nt-pm%ni)*SUM( pm%r%Edata(:,:,:,:,pm%ni+1:pm%nt)**2 )
!			print *, 'J1 = ', this%J1
!		end if
!	end subroutine
!
!	subroutine dMPE(adj,pm,k)
!		type(adjoint), intent(inout) :: adj
!		type(PM3D), intent(in) :: pm
!		integer, intent(in) :: k
!
!		if( k >= pm%ni+1 ) then
!			adj%dEs = - 2.0_mp/PRODUCT(pm%ng)/(pm%nt-pm%ni)*pm%r%Edata(:,:,:,:,k)
!		end if
!	end subroutine

!============Testmodule : Final timestep 1D Total Kinetic Energy

	subroutine TestQoI(pm,k,J)
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: k
		real(mp), intent(inout) :: J
		integer :: input

		if( k.eq.pm%nt ) then
			J = J + sum( pm%p(1)%vp(:,1)**2 )
print *, '======QoI at ',k,'-th step======'
print *, 'J=',J
print *, 'vp'
print *, pm%p(1)%vp(:,1)
		end if
	end subroutine

	subroutine dTestQoI(adj,pm,nk)
		type(adjoint), intent(inout) :: adj
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: nk

		if( nk.eq.pm%nt ) then
			adj%dp(1)%vp(:,1) = adj%dp(1)%vp(:,1) + 2.0_mp*pm%p(1)%vp(:,1)
print *, '======DQoI at ',nk,'-th step======'
print *, 'dJ'
print *, adj%dp(1)%vp(:,1)
print *, 'vp'
print *, pm%p(1)%vp(:,1)
		end if
	end subroutine

end module