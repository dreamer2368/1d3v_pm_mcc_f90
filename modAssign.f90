module modAssign

	use modSpecies
	use modMesh

	implicit none

	type pmAssign
		integer :: np
		integer :: ng
		integer :: order

		integer, allocatable :: g(:,:)
		real(mp), allocatable :: frac(:,:)
		real(mp), allocatable :: h(:)
	contains
		procedure, pass(this) :: buildAssign
		procedure, pass(this) :: destroyAssign
		procedure, pass(this) :: assignMatrix
!		procedure, pass(this) :: chargeAssign				!We handle chargeAssign subroutine globally, since it handles multiple species
		procedure, pass(this) :: forceAssign
		procedure, pass(this) :: Adj_chargeAssign
		procedure, pass(this) :: Adj_forceAssign_E
		procedure, pass(this) :: Adj_forceAssign_xp
	end type

contains

	subroutine buildAssign(this,ng,order)
		class(pmAssign), intent(out) :: this
		integer, intent(in) :: ng, order

		this%ng = ng
		this%order = order

		allocate(this%g(1,1))
		allocate(this%frac(1,1))
		allocate(this%h(1))
		this%g = 0
		this%frac = 0.0_mp
		this%h = 0.0_mp
	end subroutine

	subroutine destroyAssign(this)
		class(pmAssign), intent(inout) :: this

		deallocate(this%g)
		deallocate(this%frac)
		deallocate(this%h)
	end subroutine

	subroutine assignMatrix(this,m,xp)
		class(pmAssign), intent(inout) :: this
		type(mesh), intent(inout) :: m
		real(mp), intent(inout) :: xp(:)

		SELECT CASE(this%order)
			CASE(1)
				call assign_CIC(this,m,xp)
			CASE(2)
				call assign_TSC(this,m,xp)
		END SELECT
	end subroutine

	subroutine assign_CIC(this,m,xp)	!apply BC and create assignment matrix
		type(pmAssign), intent(inout) :: this
		type(mesh), intent(inout) :: m
		real(mp), intent(inout) :: xp(:)
		integer :: i, np
		integer :: g1, gl, gr
		real(mp) :: fracl, fracr		!fraction for left grid point
		real(mp) :: h

		np = size(xp)
		this%np = np
		deallocate(this%g)
		deallocate(this%frac)
		deallocate(this%h)
		allocate(this%g(np,this%order+1))
		allocate(this%frac(np,this%order+1))
		allocate(this%h(np))

		!assignment matrix
		do i=1,this%np
			g1 = FLOOR(xp(i)/m%dx - 0.5_mp)+1
			gl = g1
			gr = gl+1

			h = xp(i)/m%dx - g1 + 0.5_mp
			fracl = 1.0_mp - ABS(h)
			fracr = 1.0_mp - fracl

			this%h(i) = h
			this%g(i,:) = (/ gl, gr /)
			this%frac(i,:) = (/ fracl, fracr /)
		end do
	end subroutine

	subroutine assign_TSC(this,m,xp)	!apply BC and create assignment matrix
		type(pmAssign), intent(inout) :: this
		type(mesh), intent(inout) :: m
		real(mp), intent(inout) :: xp(:)
		integer :: i, np
		integer :: g		!nearest grid point
		integer :: gl		!left grid point
		integer :: gr		!right grid point
		real(mp) :: fracl		!fraction for left grid point
		real(mp) :: frac0		!fraction for nearest grid point
		real(mp) :: fracr		!fraction for right grid point
		real(mp) :: h			!distance from nearest grid point

		np = size(xp)
		this%np = np
		deallocate(this%g)
		deallocate(this%frac)
		deallocate(this%h)
		allocate(this%g(np,this%order+1))
		allocate(this%frac(np,this%order+1))
		allocate(this%h(np))
		!assignment matrix
		do i=1,this%np
			g = FLOOR(xp(i)/m%dx + 0.5_mp)+1
			gl = g-1
			gr = g+1
			h = xp(i)/m%dx - g + 1.0_mp

			frac0 = 0.75_mp - h*h
			fracl = 0.5_mp*(0.5_mp-h)*(0.5_mp-h)
			fracr = 0.5_mp*(0.5_mp+h)*(0.5_mp+h)

			this%g(i,:) = (/gl,g,gr/)
			this%h(i) = h
			this%frac(i,:) = (/ fracl, frac0, fracr /)
		end do
	end subroutine

	subroutine chargeAssign(this,p,m)
		type(pmAssign), intent(inout) :: this(:)
		type(species), intent(inout) :: p(:)
		type(mesh), intent(inout) :: m
		integer :: i,ip

		m%rho = 0.0_mp
		do ip = 1, size(p)
			do i=1,p(ip)%np
				m%rho( this(ip)%g(i,:) ) = m%rho( this(ip)%g(i,:) ) + p(ip)%spwt(i)*p(ip)%qs/m%dx*this(ip)%frac(i,:)
			end do
		end do
	end subroutine

	subroutine forceAssign(this,p,m)
		class(pmAssign), intent(inout) :: this
		type(species), intent(inout) :: p
		type(mesh), intent(in) :: m
		integer :: i

		p%Ep = 0.0_mp
		do i=1,this%np
			p%Ep(i) = sum( this%frac(i,:)*m%E(this%g(i,:)) )
		end do
	end subroutine

!======================= Adjoint assignment ==========================================================

	subroutine Adj_chargeAssign(this,p,m,rhos,xps)
		class(pmAssign), intent(inout) :: this
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
		class(pmAssign), intent(inout) :: this
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
		class(pmAssign), intent(inout) :: this
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

!======================= Forward sensitivity calculation =============================================


end module
