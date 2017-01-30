module modAdj

	use modPM1D

	implicit none

	type adjoint
		integer :: nt, ni, n, ng
		real(mp) :: dts

		type(species), allocatable :: p(:), dp(:)
		type(mesh) :: m, dm
	contains
		procedure, pass(this) :: buildAdjoint
		procedure, pass(this) :: destroyAdjoint
		procedure, pass(this) :: reset_Dadj
		procedure, pass(this) :: Adj_accel
		procedure, pass(this) :: Adj_move
	end type

contains

	subroutine buildAdjoint(this,pm)
		class(adjoint), intent(out) :: this
		type(PM1D), intent(in) :: pm
		integer :: i
		real(mp), allocatable :: vec(:), vec2(:,:)

		this%ng = pm%ng
		this%n = pm%n
		this%nt = pm%nt
		this%ni = pm%ni
		this%dts = -pm%dt

		allocate(this%p(pm%n))
		allocate(this%dp(pm%n))
		do i=1,pm%n
			call buildSpecies(this%p(i),pm%p(i)%qs,pm%p(i)%ms)
			call buildSpecies(this%dp(i),0.0_mp,0.0_mp)
			allocate(vec(pm%p(i)%np))
			allocate(vec2(pm%p(i)%np,3))
			vec = 0.0_mp
			vec2 = 0.0_mp
			call setSpecies(this%p(i),pm%p(i)%np,vec,vec2,pm%p(i)%spwt)
			call setSpecies(this%dp(i),pm%p(i)%np,vec,vec2,pm%p(i)%spwt)
			deallocate(vec)
			deallocate(vec2)
		end do
		call buildMesh(this%m,pm%L,pm%ng,pm%mBCindex)
		call buildMesh(this%dm,pm%L,pm%ng,pm%mBCindex)
	end subroutine

	subroutine destroyAdjoint(this)
		class(adjoint), intent(inout) :: this
		integer :: i

		call destroyMesh(this%m)
		call destroyMesh(this%dm)
		do i=1,this%n
			call destroySpecies(this%p(i))
			call destroySpecies(this%dp(i))
		end do
		deallocate(this%p)
		deallocate(this%dp)
	end subroutine

!	subroutine recordAdjoint(this,adj,k)
!		type(recordData), intent(inout) :: this
!		type(adjoint), intent(in) :: adj
!		integer, intent(in) :: k
!
!		this%xpsdata(:,k) = adj%xps
!		this%vpsdata(:,k) = adj%vps
!	end subroutine

	subroutine reset_Dadj(this)
		class(adjoint), intent(inout) :: this
		integer :: i

		do i=1,this%n
			this%dp(i)%xp = 0.0_mp
			this%dp(i)%vp = 0.0_mp
			this%dp(i)%Ep = 0.0_mp
		end do
		this%dm%E = 0.0_mp
		this%dm%rho = 0.0_mp
		this%dm%phi = 0.0_mp
	end subroutine

	subroutine Adj_accel(this)
		class(adjoint), intent(inout) :: this
		integer :: i

		do i=1,this%n
			this%p(i)%vp(:,1) = this%p(i)%vp(:,1) + this%dts*( -this%p(i)%xp + this%dp(i)%vp(:,1) )
		end do
	end subroutine

	subroutine Adj_move(this)
		class(adjoint), intent(inout) :: this
		integer :: i

		do i=1,this%n
			this%p(i)%xp = this%p(i)%xp + this%dts*( this%dp(i)%xp )
		end do
	end subroutine

end module
