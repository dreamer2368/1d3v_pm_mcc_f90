module modSpecies

	use constants

	implicit none

	type species
		integer :: np
		real(mp), allocatable :: xp(:)
		real(mp), allocatable :: vp(:,:)
		real(mp), allocatable :: Ep(:)
		real(mp), allocatable :: spwt(:)

		real(mp) :: ms, qs
	contains
		procedure, pass(this) :: buildSpecies
		procedure, pass(this) :: setSpecies
		procedure, pass(this) :: destroySpecies
		procedure, pass(this) :: moveSpecies
		procedure, pass(this) :: accelSpecies
		procedure, pass(this) :: appendSpecies
	end type

contains

	subroutine buildSpecies(this,qs,ms)
		class(species), intent(out) :: this
		real(mp), intent(in) :: ms, qs

		this%ms = ms
		this%qs = qs

        if( print_pm_output ) print *, 'Species built up: ms=',ms,', qs=',qs
	end subroutine

	subroutine setSpecies(this,np0,xp0,vp0,spwt0)
		class(species), intent(inout) :: this
		integer, intent(in) :: np0
		real(mp), intent(in) :: xp0(np0), spwt0(np0)
		real(mp), intent(in) :: vp0(np0,3)

		if( allocated(this%xp) ) deallocate(this%xp)
		if( allocated(this%vp) ) deallocate(this%vp)
		if( allocated(this%Ep) ) deallocate(this%Ep)
		if( allocated(this%spwt) ) deallocate(this%spwt)

		this%np = np0
		allocate(this%xp(np0))
		allocate(this%spwt(np0))
		allocate(this%vp(np0,3))
		allocate(this%Ep(np0))

		this%xp = xp0
		this%vp = vp0
		this%Ep = 0.0_mp
		this%spwt = spwt0
	end subroutine

	subroutine destroySpecies(this)
		class(species), intent(inout) :: this

		if( allocated(this%xp) ) deallocate(this%xp)
		if( allocated(this%vp) ) deallocate(this%vp)
		if( allocated(this%Ep) ) deallocate(this%Ep)
		if( allocated(this%spwt) ) deallocate(this%spwt)
		this%np = 0
	end subroutine

	subroutine moveSpecies(this,dt)
		class(species), intent(inout) :: this
		real(mp), intent(in) :: dt
        real(mp) :: time1, time2

		call CPU_TIME(time1)
		this%xp = this%xp + dt*this%vp(:,1)
		call CPU_TIME(time2)
        timeProfile(1) = timeProfile(1) + (time2-time1)
        functionCalls(1) = functionCalls(1) + 1
	end subroutine

	subroutine accelSpecies(this,dt)
		class(species), intent(inout) :: this
		real(mp), intent(in) :: dt
        real(mp) :: time1, time2

		call CPU_TIME(time1)

		this%vp(:,1) = this%vp(:,1) + dt*this%qs/this%ms*this%Ep

		call CPU_TIME(time2)
        timeProfile(7) = timeProfile(7) + (time2-time1)
        functionCalls(7) = functionCalls(7) + 1
	end subroutine

	subroutine appendSpecies(this,dN,xp0,vp0,spwt0)
		class(species), intent(inout) :: this
		integer, intent(in) :: dN
		real(mp), intent(in) :: xp0(dN), vp0(dN,3), spwt0(dN)
		integer :: newN
		real(mp), allocatable :: temp_x(:), temp_v(:,:), temp_w(:)

        real(mp) :: time1, time2

		newN = this%np+dN
		allocate(temp_x(this%np))
		allocate(temp_v(this%np,3))
		allocate(temp_w(this%np))
		temp_x = this%xp
		temp_v = this%vp
		temp_w = this%spwt

	    deallocate(this%xp)
		deallocate(this%vp)
		deallocate(this%Ep)
		deallocate(this%spwt)

		allocate(this%xp(newN))
		allocate(this%spwt(newN))
		allocate(this%vp(newN,3))
		allocate(this%Ep(newN))

		this%xp(1:this%np) = temp_x
		this%vp(1:this%np,:) = temp_v
		this%spwt(1:this%np) = temp_w

		deallocate(temp_x)
		deallocate(temp_v)
		deallocate(temp_w)

        call CPU_TIME(time1)

		this%xp(this%np+1:newN) = xp0
		this%vp(this%np+1:newN,:) = vp0
		this%spwt(this%np+1:newN) = spwt0
		this%Ep = 0.0_mp
		this%np = newN

        call CPU_TIME(time2)
        timeProfile(12) = timeProfile(12) + (time2-time1)
        functionCalls(12) = functionCalls(12) + 1
	end subroutine

end module
