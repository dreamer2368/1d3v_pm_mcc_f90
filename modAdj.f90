module modAdj

	use modPM1D

	implicit none

	type adjoint
		integer :: nt, ni, n, ng
		real(mp) :: dts

		type(species), allocatable :: p(:), dp(:)
		type(mesh) :: m, dm
	end type

contains

	subroutine buildAdjoint(this,pm)
		type(adjoint), intent(out) :: this
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
		type(adjoint), intent(inout) :: this
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

	subroutine reset_Dadj(adj)
		type(adjoint), intent(inout) :: adj
		integer :: i

		do i=1,adj%n
			adj%dp(i)%xp = 0.0_mp
			adj%dp(i)%vp = 0.0_mp
			adj%dp(i)%Ep = 0.0_mp
		end do
		adj%dm%E = 0.0_mp
		adj%dm%rho = 0.0_mp
		adj%dm%phi = 0.0_mp
	end subroutine

	subroutine Adj_chargeAssign(this,p,m,rhos,xps)
		type(pmAssign), intent(inout) :: this
		type(species), intent(in) :: p
		type(mesh), intent(in) :: m
		real(mp), intent(in) :: rhos(this%ng)
		real(mp), intent(out) :: xps(this%np)
		real(mp) :: dxps(this%np)
		real(mp) :: dV
		integer :: i

		dV = m%dx
		dxps = 0.0_mp
		dxps = 1.0_mp/m%dx*( rhos( this%g(:,2) ) - rhos( this%g(:,1) ) )
		dxps = - p%qs*p%spwt/dV*dxps
		xps = xps + dxps
	end subroutine

	subroutine Adj_forceAssign_E(this,Eps,Es)
		type(pmAssign), intent(inout) :: this
		real(mp), intent(in) :: Eps(this%np)
		real(mp), intent(inout) :: Es(this%ng)
		integer :: i, g(2)

!		Es = 0.0_mp
		do i=1,this%np
			g = this%g(i,:)
			Es( g ) = Es( g ) + Eps(i)*this%frac(i,:)
		end do
	end subroutine

	subroutine Adj_forceAssign_xp(this,m,E,Eps,xps)
		type(pmAssign), intent(inout) :: this
		type(mesh), intent(in) :: m
		real(mp), intent(in) :: E(this%ng)
		real(mp), intent(in) :: Eps(this%np)
		real(mp), intent(out) :: xps(this%np)
		real(mp) :: dxps(this%np)
		integer :: i

		dxps = 0.0_mp
		!sum : sum in each direction --- this rank will be added by the gradient of assignment
		do i=1,this%np
			dxps(i) = dxps(i) + Eps(i)/m%dx*( E(this%g(i,2)) - E(this%g(i,1)) )			!gradient of fraction = +/- 1/dx
		end do
		dxps = - dxps

		xps = xps + dxps
	end subroutine

	!================   Adjoint Mesh solver   ========================================
	subroutine solveMesh_Adj(m,eps)
		type(mesh), intent(inout) :: m
		real(mp), intent(in) :: eps

		select case(m%BCindex)
			case(0)
				call solveMesh_Adj_periodic(m,eps)
	!		case(1)
	!			call solveMesh_D_D(this,eps)
	!		case(2)
	!			call solveMesh_D_N(this,eps)
		end select
	end subroutine

	subroutine solveMesh_Adj_periodic(m,eps)
		type(mesh), intent(inout) :: m
		real(mp), intent(in) :: eps
		real(mp) :: rhs(m%ng), rho1(m%ng-1)

		rhs = multiplyD(m%E,m%dx,m%BCindex)
		call CG_K(multiplyK,rho1,rhs(1:m%ng-1),m%dx)
		m%rho(1:m%ng-1) = rho1
		m%rho(m%ng) = 0.0_mp
		m%rho = - m%rho/eps
	end subroutine

	subroutine Adj_accel(adj)
		type(adjoint), intent(inout) :: adj
		integer :: i

		do i=1,adj%n
			adj%p(i)%vp(:,1) = adj%p(i)%vp(:,1) + adj%dts*( -adj%p(i)%xp + adj%dp(i)%vp(:,1) )
		end do
	end subroutine

	subroutine Adj_move(adj)
		type(adjoint), intent(inout) :: adj
		integer :: i

		do i=1,adj%n
			adj%p(i)%xp = adj%p(i)%xp + adj%dts*( adj%dp(i)%xp )
		end do
	end subroutine

end module
