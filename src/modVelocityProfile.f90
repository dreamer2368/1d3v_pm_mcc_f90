module modVelocityProfile

	use constants

	implicit none

	abstract interface
		function velocityProfile(vp) result(fp)
			use constants
			real(mp), intent(in) :: vp(:)
            real(mp) :: fp(size(vp))
		end function
	end interface

contains

    function GaussianProfile(vp) result(fp)
        real(mp), intent(in) :: vp(:)
        real(mp) :: fp(size(vp))

        fp = EXP( - vp**2/2.0_mp )
    end function

end module
