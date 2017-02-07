module modTarget

	use modAdj
	use random

	implicit none

contains
!==============Default=====================================================

	subroutine Null_input(this,k,str)
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
	end subroutine

	subroutine Null_Dinput(adj,pm,k,str)
		type(adjoint), intent(inout) :: adj
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
	end subroutine

	subroutine Null_dJdA(adj,pm,k,str,grad)
		type(adjoint), intent(in) :: adj
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
		real(mp), intent(inout) :: grad(:)
	end subroutine

!==============Initial electron temperature================================

	subroutine Te(pm,k,str)
		type(PM1D), intent(inout) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		SELECT CASE (str)
			CASE('xp')
				if(k.eq.1) then
					pm%p(1)%vp = pm%p(1)%vp*( 1.0_mp + pm%A0(2) )
				end if
		END SELECT
	end subroutine

	subroutine dTe(adj,pm,k,str)
		type(adjoint), intent(inout) :: adj
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
	end subroutine

	subroutine dTedA(adj,pm,k,str,grad)
		type(adjoint), intent(in) :: adj
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
		real(mp), intent(inout) :: grad(:)

		if( str.eq.'after' .and. k.eq.0 ) then
			grad(1) = -sum( adj%p(1)%vp(:,1)*pm%p(1)%vp(:,1) )/pm%dt
		end if
	end subroutine

!==============Wave perturbation on position===============================

	subroutine IC_wave(this,k,str)
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		SELECT CASE (str)
			CASE('xp')
				if(k==this%ni) then
					this%p(1)%xp = this%p(1)%xp + this%dt*this%A0(2)*this%L/this%p(1)%np*SIN(4.0_mp*pi*this%p(1)%xp/this%L)		!xp_(Ni-1) : k=Ni, xp_Ni : k=(Ni+1)
				end if
		END SELECT
	end subroutine

	subroutine dIC_wave(adj,pm,k,str)
		type(adjoint), intent(inout) :: adj
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		select case (str)
			case ('xp')
				if( k .eq. pm%ni-1 ) then
					adj%dp(1)%xp = adj%dp(1)%xp - adj%p(1)%xp*(pm%A0(2)*pm%L/pm%p(1)%np*4.0_mp*pi/pm%L)*COS(4.0_mp*pi*pm%p(1)%xp/pm%L)
				end if
		end select
	end subroutine

	subroutine dIC_wave_dB(adj,pm,k,str,grad)
		type(adjoint), intent(in) :: adj
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
		real(mp), intent(inout) :: grad(:)

		if( str.eq.'before' .and. k.eq.pm%ni-1 ) then
			grad(1) = grad(1) + SUM( - pm%L/pm%p(1)%np*SIN( 4.0_mp*pi*pm%p(1)%xp/pm%L )*adj%p(1)%xp )
		end if
	end subroutine

!==============RF discharge current source=================================

   subroutine RF_current(this,k,str)
      type(PM1D), intent(inout) :: this
      integer, intent(in) :: k
      character(len=*), intent(in) :: str
      real(mp) :: dQwall

      !jwall = A0(3)*sin( 2*pi*A0(4)*t )
      dQwall = this%A0(3)*SIN( 2.0_mp*pi*this%A0(4)*k*this%dt )*this%dt

      SELECT CASE (str)
         CASE('rho_back')
            !current flows from x=L(ng) to x=0(1)
     			this%m%rho_back(1) = this%m%rho_back(1) + dQwall
     			this%m%rho_back(this%ng) = this%m%rho_back(this%ng) - dQwall
      END SELECT
   end subroutine

!==================Injection from left wall=================================

	subroutine Injection(this,k,str)
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
		real(mp) :: vT_e, vT_i
		integer, parameter :: NInject=1E3
		real(mp) :: temp_x(NInject), temp_v(NInject,3), temp_w(NInject)
		real(mp), allocatable :: new_x(:), new_v(:,:), new_w(:)
		integer :: np

		SELECT CASE (str)
			CASE('xp')
!				vT_e = SQRT(this%A0(1)*2.0_mp/this%p(1)%ms)
!				vT_i = SQRT(this%A0(1)*2.0_mp/this%p(2)%ms)
!
!				temp_v = ABS(randn(NInject,3))*vT_e
!				call RANDOM_NUMBER(temp_x)
!				temp_x = temp_x*temp_v(:,1)*this%dt
!				temp_w = this%p(1)%spwt(1)
!
!				np = this%p(1)%np
!
!				allocate(new_x(np+NInject))
!				allocate(new_v(np+NInject,3))
!				allocate(new_w(np+NInject))
!
!				new_x(1:np) = this%p(1)%xp
!				new_v(1:np,:) = this%p(1)%vp
!				new_w(1:np) = this%p(1)%spwt
!			
!				new_x(np+1:np+NInject) = temp_x
!				new_v(np+1:np+NInject,:) = temp_v
!				new_w(np+1:np+NInject) = temp_w
!
!				call this%p(1)%setSpecies(np+NInject,new_x,new_v,new_w)
!
!				deallocate(new_x)
!				deallocate(new_v)
!				deallocate(new_w)
!
!				temp_v = ABS(randn(NInject,3))*vT_i
!				call RANDOM_NUMBER(temp_x)
!				temp_x = temp_x*temp_v(:,1)*this%dt
!				temp_w = this%p(2)%spwt(1)
!
!				np = this%p(2)%np
!
!				allocate(new_x(np+NInject))
!				allocate(new_v(np+NInject,3))
!				allocate(new_w(np+NInject))
!
!				new_x(1:np) = this%p(2)%xp
!				new_v(1:np,:) = this%p(2)%vp
!				new_w(1:np) = this%p(2)%spwt
!			
!				new_x(np+1:np+NInject) = temp_x
!				new_v(np+1:np+NInject,:) = temp_v
!				new_w(np+1:np+NInject) = temp_w
!
!				call this%p(2)%setSpecies(np+NInject,new_x,new_v,new_w)
!
!				deallocate(new_x)
!				deallocate(new_v)
!				deallocate(new_w)

			CASE('rho_back')
				if( k>this%ni ) then
					this%m%rho_back(1) = 0.0_mp
					this%m%rho_back(this%m%ng) = -0.1_mp
				else
					this%m%rho_back = 0.0_mp
				end if
		END SELECT	
	end subroutine

end module
