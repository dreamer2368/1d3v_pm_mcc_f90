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

	!Uniform on x in [(0.5-A(3))*L, (0.5+A(3))*L]
	!Rayleigh on v with \sigma = A(1),A(2) (electron,ion)
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
		xp_add = xp_add*2.0_mp*pm%A0(3)*pm%L + (0.5_mp - pm%A0(3))*pm%L
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

	!Uniform on x in [(0.5-A(3))*L, (0.5+A(3))*L]
	!Maxwellian on v with \sigma = A(1),A(2) (electron,ion)
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
		xp_add = xp_add*2.0_mp*pm%A0(3)*pm%L + (0.5_mp - pm%A0(3))*pm%L
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
        Period = 120.0_mp
        timeStart = 30.0_mp
!        Period = 0.4_mp*pm%nt*pm%dt
!        timeStart = 0.1_mp*pm%nt*pm%dt
        timeEnd = timeStart + 2.0_mp*Period
        convectionVelocity = 0.0_mp
        if( (time.ge.timeStart) .and. (time.le.timeEnd) ) then
            convectionVelocity = pm%A0(5)*( 0.5_mp - 0.5_mp*COS( 2.0_mp*pi*(time-timeStart)/Period ) )
        end if

		spwt_add = pm%p(2)%spwt(1)
		vp_add = pm%A0(2)*randn(Nadd,3) + convectionVelocity
		call pm%p(2)%appendSpecies(Nadd,xp_add,vp_add,spwt_add)

		deallocate(xp_add)
		deallocate(vp_add)
		deallocate(spwt_add)
	end subroutine

	subroutine Modified_Maxwellian2_Sensitivity(pm,k,NIonFluxR,ionFluxR,NSensitivityFluxR,sensitivityFluxR)
		class(PM1D), intent(inout) :: pm
        integer, intent(in) :: k, NIonFluxR, NSensitivityFluxR
        real(mp), intent(in) :: ionFluxR, sensitivityFluxR
		integer :: newN, Nadd, i
		real(mp), allocatable :: xp_add(:), vp_add(:,:), spwt_add(:)
        real(mp) :: time, convectionVelocity, timeStart, timeEnd, Period

        time = k*pm%dt
        Period = 0.2_mp*pm%nt*pm%dt
!        Period = 30.0_mp
        timeStart = 0.1_mp*pm%nt*pm%dt
!        timeStart = 10.0_mp
        timeEnd = timeStart + 4.0_mp*Period
        if( (time.ge.timeStart) .and. (time.le.timeEnd) ) then
            Nadd = 2*((NIonFluxR+1)/2)
        else
            Nadd = 0
        end if
        newN = NSensitivityFluxR + Nadd
		allocate(xp_add(newN))
		allocate(spwt_add(newN))
		allocate(vp_add(newN,3))

		call RANDOM_NUMBER(xp_add)
		xp_add = xp_add*pm%A0(3)*pm%L
		spwt_add(1:NSensitivityFluxR) = sensitivityFluxR/NSensitivityFluxR
		vp_add(1:NSensitivityFluxR,1) = pm%A0(1)*randn(NSensitivityFluxR)
		call pm%p(1)%appendSpecies(newN-Nadd,                                   &
                                    xp_add(1:newN-Nadd),                        &
                                    vp_add(1:newN-Nadd,:),                      &
                                    spwt_add(1:newN-Nadd) )

        convectionVelocity = 0.0_mp
        if( (time.ge.timeStart) .and. (time.le.timeEnd) ) then
            convectionVelocity = pm%A0(5)*SIN( 2.0_mp*pi*(time-timeStart)/Period )
        end if

		vp_add(1:NSensitivityFluxR,1) = pm%A0(2)*randn(NSensitivityFluxR) + convectionVelocity
        vp_add(newN-Nadd+1:newN-Nadd/2,1) = pm%A0(2)*randr(Nadd/2) + convectionVelocity
        vp_add(newN-Nadd/2+1:newN,1) = -pm%A0(2)*randr(Nadd/2) + convectionVelocity

        spwt_add(newN-Nadd+1:newN-Nadd/2) = ionFluxR*SQRT(2.0_mp/pi)/pm%A0(2)*SIN( 2.0_mp*pi*(time-timeStart)/Period )/Nadd
        spwt_add(newN-Nadd/2+1:newN) = -ionFluxR*SQRT(2.0_mp/pi)/pm%A0(2)*SIN( 2.0_mp*pi*(time-timeStart)/Period )/Nadd

		call pm%p(2)%appendSpecies(newN,xp_add,vp_add,spwt_add)

		deallocate(xp_add)
		deallocate(vp_add)
		deallocate(spwt_add)
	end subroutine

end module
