module modTarget

	use modAdj

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

end module
