module modSource

	use modAdj

	implicit none

contains
!==============Increase thermal velocity===================================

	subroutine V_Th(this,k,str)
		type(PM3D), intent(inout) :: this
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		SELECT CASE (str)
			CASE('vp')
				if(k==this%ni) then
					this%p%vp = this%p%vp + this%B0*this%p%vp
				end if
		END SELECT
	end subroutine

	subroutine dV_Th(adj,pm,k,str)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		select case (str)
			case ('vp')
				if( k .eq. pm%ni-1 ) then
					adj%dvps = adj%vps*pm%B0/(-pm%dt)
				end if
		end select
	end subroutine

	subroutine dV_Th_dB(adj,pm,dJdA)
		type(adjoint), intent(in) :: adj
		type(PM3D), intent(in) :: pm
		real(mp), intent(out) :: dJdA

		dJdA = SUM( - pm%r%vpdata(:,:,pm%ni-1)*pm%r%vpsdata(:,:,pm%ni)/pm%dt )
	end subroutine

!==============Wave perturbation on position===============================

	subroutine IC_wave(this,k,str)
		type(PM3D), intent(inout) :: this
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		SELECT CASE (str)
			CASE('xp')
				if(k==this%ni) then
					this%p%xp(:,1) = this%p%xp(:,1) + this%dt*this%B0*this%L(1)/this%n*SIN(4.0_mp*pi*this%p%xp(:,1)/this%L(1)*mode)		!xp_(Ni-1) : k=Ni, xp_Ni : k=(Ni+1)
				end if
		END SELECT
	end subroutine

	subroutine dIC_wave(adj,pm,k,str)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		select case (str)
			case ('xp')
				if( k .eq. pm%ni-1 ) then
					adj%dxps = adj%dxps - adj%xps*(pm%B0*pm%L(1)/pm%n*4.0_mp*pi*mode/pm%L(1))*COS(4.0_mp*pi*mode*pm%r%xpdata(:,:,k)/pm%L(1))
				end if
		end select
	end subroutine

	subroutine dIC_wave_dB(adj,pm,dJdA)
		type(adjoint), intent(in) :: adj
		type(PM3D), intent(in) :: pm
		real(mp), intent(out) :: dJdA

		dJdA = SUM( - pm%L(1)/pm%n*SIN( 4.0_mp*pi*mode*pm%r%xpdata(:,1,pm%ni-1)/pm%L(1) )*pm%r%xpsdata(:,1,pm%ni) )
	end subroutine

!================Two Particle Orbit radius===========================

	subroutine Orbit_radius(pm,k,str)
		type(PM3D), intent(inout) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		select case (str)
			case('xp')
				if( k.eq.2 ) then
					pm%p%xp(1,1) = pm%p%xp(1,1) - pm%B0*0.5_mp*(pm%p%xp(2,1)-pm%p%xp(1,1))
					pm%p%xp(2,1) = pm%p%xp(2,1) + pm%B0*0.5_mp*(pm%p%xp(2,1)-pm%p%xp(1,1))
				end if
		end select
	end subroutine

	subroutine dOrbit(adj,pm,k,str)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

	end subroutine

	subroutine dOrbit_dB(adj,pm,dJdA)
		type(adjoint), intent(in) :: adj
		type(PM3D), intent(in) :: pm
		real(mp), intent(out) :: dJdA
		real(mp) :: r0

		r0 = pm%r%xpdata(2,1,1) - pm%r%xpdata(1,1,1)
		dJdA = - ( -0.5_mp*r0*pm%r%xpsdata(1,1,2) + 0.5_mp*r0*pm%r%xpsdata(2,1,2) )/pm%dt
	end subroutine

end module