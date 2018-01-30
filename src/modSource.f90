module modSource

	use modPM1D
	use random

	implicit none

	abstract interface
		subroutine source(pm)
			use modPM1D
			class(PM1D), intent(inout) :: pm
		end subroutine
	end interface

contains

	subroutine Null_source(pm)
		class(PM1D), intent(inout) :: pm
	end subroutine

!======= [Spatial]_[Velocity] source distribution ==================

	!Uniform on x in [(0.5-A(1))*L, (0.5+A(1))*L]
	!Rayleigh on v with \sigma = A(2),A(3) (electron,ion)
	!Keep A(4) number of ion
	subroutine PartialUniform_Rayleigh(pm)
		class(PM1D), intent(inout) :: pm
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
	subroutine PartialUniform_Rayleigh2(pm)
		class(PM1D), intent(inout) :: pm
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
	subroutine PartialUniform_Maxwellian(pm)
		class(PM1D), intent(inout) :: pm
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
	subroutine PartialUniform_Maxwellian2(pm)
		class(PM1D), intent(inout) :: pm
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

end module
