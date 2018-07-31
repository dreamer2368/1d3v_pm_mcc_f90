module modSource

	use modPM1D
	use random

	implicit none

	abstract interface
		subroutine source(pm,k)
			use modPM1D
			class(PM1D), intent(inout) :: pm
            integer, intent(in) :: k
		end subroutine
	end interface

contains

	subroutine Null_source(pm,k)
		class(PM1D), intent(inout) :: pm
        integer, intent(in) :: k
	end subroutine

!======= [Spatial]_[Velocity] source distribution ==================

	!Uniform on x in [(0.5-A(1))*L, (0.5+A(1))*L]
	!Rayleigh on v with \sigma = A(2),A(3) (electron,ion)
	!Keep A(4) number of ion
	subroutine PartialUniform_Rayleigh(pm,k)
		class(PM1D), intent(inout) :: pm
        integer, intent(in) :: k
		integer :: Nadd, newN, i
		real(mp), allocatable :: xp_add(:), vp_add(:,:), spwt_add(:)

		Nadd = floor(pm%A0(4)-pm%p(2)%np)
		allocate(xp_add(Nadd))
		allocate(spwt_add(Nadd))
		allocate(vp_add(Nadd,3))

		call RANDOM_NUMBER(xp_add)
		xp_add = xp_add*2.0_mp*pm%A0(1)*pm%L + (0.5_mp - pm%A0(1))*pm%L
		spwt_add = pm%p(1)%spwt(1)
		vp_add = pm%A0(2)*randr(Nadd,3)
		call pm%p(1)%appendSpecies(Nadd,xp_add,vp_add,spwt_add)
	
		spwt_add = pm%p(2)%spwt(1)
		vp_add = pm%A0(3)*randr(Nadd,3)
		call pm%p(2)%appendSpecies(Nadd,xp_add,vp_add,spwt_add)

		deallocate(xp_add)
		deallocate(vp_add)
		deallocate(spwt_add)
	end subroutine

	!Uniform on x in [0, A(3)]
	!Rayleigh on v with \sigma = A(1),A(2) (electron,ion)
	!Keep A(4) number of ion
	subroutine PartialUniform_Rayleigh2(pm,k)
		class(PM1D), intent(inout) :: pm
        integer, intent(in) :: k
		integer :: Nadd, newN, i
		real(mp), allocatable :: xp_add(:), vp_add(:,:), spwt_add(:)

		Nadd = floor(pm%A0(4)-pm%p(2)%np)
		allocate(xp_add(Nadd))
		allocate(spwt_add(Nadd))
		allocate(vp_add(Nadd,3))

		call RANDOM_NUMBER(xp_add)
		xp_add = xp_add*pm%A0(3)*pm%L
		spwt_add = pm%p(1)%spwt(1)
		vp_add = pm%A0(1)*randr(Nadd,3)
		call pm%p(1)%appendSpecies(Nadd,xp_add,vp_add,spwt_add)

		spwt_add = pm%p(2)%spwt(1)
		vp_add = pm%A0(2)*randr(Nadd,3)
		call pm%p(2)%appendSpecies(Nadd,xp_add,vp_add,spwt_add)

		deallocate(xp_add)
		deallocate(vp_add)
		deallocate(spwt_add)
	end subroutine

	!Uniform on x in [(0.5-A(1))*L, (0.5+A(1))*L]
	!Maxwellian on v with \sigma = A(2),A(3) (electron,ion)
	!Keep A(4) number of ion
	subroutine PartialUniform_Maxwellian(pm,k)
		class(PM1D), intent(inout) :: pm
        integer, intent(in) :: k
		integer :: Nadd, newN, i
		real(mp), allocatable :: xp_add(:), vp_add(:,:), spwt_add(:)

		Nadd = floor(pm%A0(4)-pm%p(2)%np)
		allocate(xp_add(Nadd))
		allocate(spwt_add(Nadd))
		allocate(vp_add(Nadd,3))

		call RANDOM_NUMBER(xp_add)
		xp_add = xp_add*2.0_mp*pm%A0(1)*pm%L + (0.5_mp - pm%A0(1))*pm%L
		spwt_add = pm%p(1)%spwt(1)
		vp_add = pm%A0(2)*randn(Nadd,3)
		call pm%p(1)%appendSpecies(Nadd,xp_add,vp_add,spwt_add)

		spwt_add = pm%p(2)%spwt(1)
		vp_add = pm%A0(3)*randn(Nadd,3)
		call pm%p(2)%appendSpecies(Nadd,xp_add,vp_add,spwt_add)

		deallocate(xp_add)
		deallocate(vp_add)
		deallocate(spwt_add)
	end subroutine

	!Uniform on x in [0, A(3)]
	!Maxwellian on v with \sigma = A(1),A(2) (electron,ion)
	!Keep A(4) number of ion
	subroutine PartialUniform_Maxwellian2(pm,k)
		class(PM1D), intent(inout) :: pm
        integer, intent(in) :: k
		integer :: Nadd, newN, i
		real(mp), allocatable :: xp_add(:), vp_add(:,:), spwt_add(:)

		Nadd = floor(pm%A0(4)-pm%p(2)%np)
		allocate(xp_add(Nadd))
		allocate(spwt_add(Nadd))
		allocate(vp_add(Nadd,3))

		call RANDOM_NUMBER(xp_add)
		xp_add = xp_add*pm%A0(3)*pm%L
		spwt_add = pm%p(1)%spwt(1)
		vp_add = pm%A0(1)*randn(Nadd,3)
		call pm%p(1)%appendSpecies(Nadd,xp_add,vp_add,spwt_add)

		spwt_add = pm%p(2)%spwt(1)
		vp_add = pm%A0(2)*randn(Nadd,3)
		call pm%p(2)%appendSpecies(Nadd,xp_add,vp_add,spwt_add)

		deallocate(xp_add)
		deallocate(vp_add)
		deallocate(spwt_add)
	end subroutine

	subroutine Modified_Maxwellian(pm,k)
		class(PM1D), intent(inout) :: pm
        integer, intent(in) :: k
		integer :: Nadd, newN, i
		real(mp), allocatable :: xp_add(:), vp_add(:,:), spwt_add(:)
        real(mp) :: time, convectionVelocity, timeStart, timeEnd, Period

		Nadd = floor(pm%A0(4)-pm%p(2)%np)
		allocate(xp_add(Nadd))
		allocate(spwt_add(Nadd))
		allocate(vp_add(Nadd,3))

		call RANDOM_NUMBER(xp_add)
		xp_add = xp_add*2.0_mp*pm%A0(1)*pm%L + (0.5_mp - pm%A0(1))*pm%L
		spwt_add = pm%p(1)%spwt(1)
		vp_add = pm%A0(2)*randn(Nadd,3)
		call pm%p(1)%appendSpecies(Nadd,xp_add,vp_add,spwt_add)

        time = k*pm%dt
        Period = 0.2_mp*pm%nt*pm%dt
        timeStart = 0.1_mp*pm%nt*pm%dt
        timeEnd = timeStart + 4.0_mp*Period
        convectionVelocity = 0.0_mp
        if( (time.ge.timeStart) .and. (time.le.timeEnd) ) then
            convectionVelocity = pm%A0(5)*SIN( 2.0_mp*pi*(time-timeStart)/Period )
        end if

		spwt_add = pm%p(2)%spwt(1)
		vp_add = pm%A0(3)*randn(Nadd,3) + convectionVelocity
		call pm%p(2)%appendSpecies(Nadd,xp_add,vp_add,spwt_add)

		deallocate(xp_add)
		deallocate(vp_add)
		deallocate(spwt_add)
	end subroutine

	subroutine Modified_Maxwellian2(pm,k)
		class(PM1D), intent(inout) :: pm
        integer, intent(in) :: k
		integer :: Nadd, newN, i
		real(mp), allocatable :: xp_add(:), vp_add(:,:), spwt_add(:)
        real(mp) :: time, convectionVelocity, timeStart, timeEnd, Period

		Nadd = floor(pm%A0(4)-pm%p(2)%np)
		allocate(xp_add(Nadd))
		allocate(spwt_add(Nadd))
		allocate(vp_add(Nadd,3))

		call RANDOM_NUMBER(xp_add)
		xp_add = xp_add*pm%A0(3)*pm%L
		spwt_add = pm%p(1)%spwt(1)
		vp_add = pm%A0(1)*randn(Nadd,3)
		call pm%p(1)%appendSpecies(Nadd,xp_add,vp_add,spwt_add)

        time = k*pm%dt
        Period = 0.2_mp*pm%nt*pm%dt
        timeStart = 0.1_mp*pm%nt*pm%dt
        timeEnd = timeStart + 4.0_mp*Period
        convectionVelocity = 0.0_mp
        if( (time.ge.timeStart) .and. (time.le.timeEnd) ) then
            convectionVelocity = pm%A0(5)*SIN( 2.0_mp*pi*(time-timeStart)/Period )
        end if

		spwt_add = pm%p(2)%spwt(1)
		vp_add = pm%A0(2)*randn(Nadd,3) + convectionVelocity
		call pm%p(2)%appendSpecies(Nadd,xp_add,vp_add,spwt_add)

		deallocate(xp_add)
		deallocate(vp_add)
		deallocate(spwt_add)
	end subroutine

end module
