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

		this%xp = this%xp + dt*this%vp(:,1)
	end subroutine

	subroutine accelSpecies(this,dt)
		class(species), intent(inout) :: this
		real(mp), intent(in) :: dt

		this%vp(:,1) = this%vp(:,1) + dt*this%qs/this%ms*this%Ep
	end subroutine

	subroutine appendSpecies(this,dN,xp0,vp0,spwt0)
		class(species), intent(inout) :: this
		integer, intent(in) :: dN
		real(mp), intent(in) :: xp0(dN), vp0(dN,3), spwt0(dN)
		integer :: newN
		real(mp), allocatable :: temp_x(:), temp_v(:,:), temp_w(:)
		newN = this%np+dN
		allocate(temp_x(newN))
		allocate(temp_v(newN,3))
		allocate(temp_w(newN))
		temp_x = 0.0_mp
		temp_v = 0.0_mp
		temp_w = 0.0_mp

		temp_x(1:this%np) = this%xp
		temp_x(this%np+1:newN) = xp0

		temp_v(1:this%np,:) = this%vp
		temp_v(this%np+1:newN,:) = vp0

		temp_w(1:this%np) = this%spwt
		temp_w(this%np+1:newN) = spwt0

		call this%setSpecies(newN,temp_x,temp_v,temp_w)

		deallocate(temp_x)
		deallocate(temp_v)
		deallocate(temp_w)
	end subroutine

end module
