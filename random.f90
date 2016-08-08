MODULE random

	use constants

	implicit none

	interface randn
		procedure randn_rank1
		procedure randn_rank2
	end interface

	interface randr
		procedure randr_rank1
		procedure randr_rank2
	end interface

contains

   subroutine init_random_seed(input)
      integer :: i, nseed, clock, add
      integer, allocatable :: seed(:)
	  integer, intent(in), optional :: input
		if( present(input) ) then
			add = input
		else
			add = 0
		end if

		call RANDOM_SEED(size=nseed)
		allocate(seed(nseed))
!      call SYSTEM_CLOCK(clock)
!		seed = clock + 3433*(/ ( i, i=1,nseed ) /)
		seed = 3433*(/ ( i, i=1,nseed ) /) + add
		call RANDOM_SEED(put=seed)
		deallocate(seed)
   end subroutine

	function randn_rank1(N) result(x)
		integer, intent(in) :: N
		real(mp) :: r1(N), r2(N), x(N)

		call RANDOM_NUMBER(r1)
		call RANDOM_NUMBER(r2)
		x = sqrt( -2.0_mp*log(r1) )*cos( 2.0_mp*pi*r2 )
	end function

	function randn_rank2(N1,N2) result(x)
		integer, intent(in) :: N1, N2
		real(mp) :: r1(N1,N2), r2(N1,N2), x(N1,N2)

		call RANDOM_NUMBER(r1)
		call RANDOM_NUMBER(r2)
		x = sqrt( -2.0_mp*log(r1) )*cos( 2.0_mp*pi*r2 )
	end function

	function randr_rank1(N) result(x)
		integer, intent(in) :: N
		real(mp) :: r1(N/2), r2(N-N/2), x(N)

		call RANDOM_NUMBER(r1)
		call RANDOM_NUMBER(r2)
		x(1:N/2) = sqrt( -2.0_mp*log(r1) )
		x(N/2+1:N) = -sqrt( -2.0_mp*log(r2) )
	end function

	function randr_rank2(N1,N2) result(x)
		integer, intent(in) :: N1, N2
		real(mp) :: r1(N1/2,N2), r2(N1-N1/2,N2), x(N1,N2)

		call RANDOM_NUMBER(r1)
		call RANDOM_NUMBER(r2)
		x(1:N1/2,:) = sqrt( -2.0_mp*log(r1) )
		x(N1/2+1:N1,:) = -sqrt( -2.0_mp*log(r2) )
	end function

END MODULE random
