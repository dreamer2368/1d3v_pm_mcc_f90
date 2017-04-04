module modAssign

	use modSpecies
	use modMesh

	implicit none

	type pmAssign
		integer :: np,ng,order,mBCidx

		integer, dimension(:,:), pointer :: g
		real(mp), dimension(:,:), pointer :: frac
		real(mp), allocatable :: h(:)
		procedure(assignMatrix), nopass, pointer :: assignMatrix
		procedure(adjustGrid), nopass, pointer :: adjustGrid
	contains
		procedure, pass(this) :: buildAssign
		procedure, pass(this) :: destroyAssign
		procedure, pass(this) :: chargeAssign				!We handle chargeAssign subroutine globally, since it handles multiple species
		procedure, pass(this) :: forceAssign
		procedure, pass(this) :: Adj_chargeAssign
		procedure, pass(this) :: Adj_forceAssign_E
		procedure, pass(this) :: Adj_forceAssign_xp
	end type

	abstract interface
		subroutine assignMatrix(xp,dx,g,frac)
			use constants
			real(mp), intent(in) :: xp, dx
			integer, intent(out) :: g(:)
			real(mp), intent(out) :: frac(:)
		end subroutine
	end interface

	abstract interface
		subroutine adjustGrid(ng,g,frac)
			use constants
			integer, intent(in) :: ng
			integer, intent(inout) :: g(:)
			real(mp), intent(inout) :: frac(:)
		end subroutine
	end interface

contains

	subroutine buildAssign(this,ng,order,mBCidx)
		class(pmAssign), intent(out) :: this
		integer, intent(in) :: ng, order, mBCidx

		this%ng = ng
		this%order = order
		this%mBCidx = mBCidx

		SELECT CASE(order)
			CASE(1)
				this%assignMatrix=>assign_CIC
			CASE(2)
				this%assignMatrix=>assign_TSC
		END SELECT

		select case(mBCidx)
			case(0)	!periodic
				this%adjustGrid=>adjustGrid_periodic
			case(1)	!Dirichlet-Dirichlet
				this%adjustGrid=>adjustGrid_absorbing
			case(2)	!Dirichlet-Neumann
				this%adjustGrid=>adjustGrid_absorbing
		end select

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

	subroutine assign_CIC(xp,dx,g,frac)	!apply BC and create assignment matrix
		real(mp), intent(in) :: xp, dx
		integer, intent(out) :: g(:)
		real(mp), intent(out) :: frac(:)
		integer :: i, np
		integer :: g1, gl, gr
		real(mp) :: fracl, fracr		!fraction for left grid point
		real(mp) :: h

		!assignment matrix
!		g1 = FLOOR(xp/dx - 0.5_mp)+1				!For periodic BC
		g1 = FLOOR(xp/dx) + 1						!For bounded BC
		gl = g1
		gr = gl+1

!		h = xp/dx - g1 + 0.5_mp						!For periodic BC
		h = xp/dx - g1 + 1.0_mp						!For bounded BC
		fracl = 1.0_mp - ABS(h)
		fracr = 1.0_mp - fracl

		g = (/ gl, gr /)
		frac = (/ fracl, fracr /)
	end subroutine

	subroutine assign_TSC(xp,dx,g,frac)	!apply BC and create assignment matrix
		real(mp), intent(in) :: xp, dx
		integer, intent(out) :: g(:)
		real(mp), intent(out) :: frac(:)
		integer :: i, np
		integer :: g0		!nearest grid point
		integer :: gl		!left grid point
		integer :: gr		!right grid point
		real(mp) :: fracl		!fraction for left grid point
		real(mp) :: frac0		!fraction for nearest grid point
		real(mp) :: fracr		!fraction for right grid point
		real(mp) :: h			!distance from nearest grid point

		!assignment matrix
		g0 = FLOOR(xp/dx + 0.5_mp)+1
		gl = g0-1
		gr = g0+1
		h = xp/dx - g0 + 1.0_mp

		frac0 = 0.75_mp - h*h
		fracl = 0.5_mp*(0.5_mp-h)*(0.5_mp-h)
		fracr = 0.5_mp*(0.5_mp+h)*(0.5_mp+h)

		g = (/gl,g0,gr/)
		frac = (/ fracl, frac0, fracr /)
	end subroutine

	subroutine chargeAssign(this,p,m)
		class(pmAssign), intent(inout) :: this
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		integer :: g(this%order+1)
		real(mp) :: frac(this%order+1)
		integer :: i,ip

		DEALLOCATE(this%g)
		DEALLOCATE(this%frac)
		this%np = p%np
		ALLOCATE(this%g(this%order+1,p%np))
		ALLOCATE(this%frac(this%order+1,p%np))
		do i=1,p%np
			CALL this%assignMatrix(p%xp(i),m%dx,g,frac)
			CALL this%adjustGrid(m%ng,g,frac)
			this%g(:,i) = g
			this%frac(:,i) = frac
			m%rho( g ) = m%rho( g ) + p%spwt(i)*p%qs/m%dx*frac
		end do
	end subroutine

	subroutine forceAssign(this,p,m)
		class(pmAssign), intent(inout) :: this
		type(species), intent(inout) :: p
		type(mesh), intent(in) :: m
		integer :: i

		p%Ep = 0.0_mp
		do i=1,this%np
			p%Ep(i) = sum( this%frac(:,i)*m%E(this%g(:,i)) )
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
		dxps = 1.0_mp/m%dx*( rhos( this%g(2,:) ) - rhos( this%g(1,:) ) )
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
			g = this%g(:,i)
			Es( g ) = Es( g ) + Eps(i)*this%frac(:,i)
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
			dxps(i) = dxps(i) + Eps(i)/m%dx*( E(this%g(2,i)) - E(this%g(1,i)) )			!gradient of fraction = +/- 1/dx
		end do
		dxps = - dxps

		xps = xps + dxps
	end subroutine

!======================= AdjustGrid for BC =============================================

	subroutine adjustGrid_periodic(ng,g,frac)
		integer, intent(in) :: ng
		integer, intent(inout) :: g(:)
		real(mp), intent(inout) :: frac(:)

			if( g(1)<1 ) then
				g(1) = g(1) + ng
			elseif( g(1)>ng ) then
				g(1) = g(1) - ng
			end if
			if( g(2)<1 ) then
				g(2) = g(2) + ng
			elseif( g(2)>ng ) then
				g(2) = g(2) - ng
			end if

		if( MINVAL(g)<1 .or. MAXVAL(g)>ng ) then
			print *, MINVAL(g), MAXVAL(g)
			print *, 'Boundary handling is failed. particle is way outside BC. stopped time stepping.'
			stop
		end if
	end subroutine

	subroutine adjustGrid_absorbing(ng,g,frac)
		integer, intent(in) :: ng
		integer, intent(inout) :: g(:)
		real(mp), intent(inout) :: frac(:)

		if( MINVAL(g)<1 .or. MAXVAL(g)>ng ) then
			print *, MINVAL(g), MAXVAL(g), frac
			print *, 'Boundary handling is failed. particle is way outside BC. stopped time stepping.'
			stop
		end if

		!adjustment for boundary : charge/(dx/2)
		if( g(1).eq.1 ) then
			frac(1) = frac(1)*2.0_mp
		end if
		if( g(2).eq.ng ) then
			frac(2) = frac(2)*2.0_mp
		end if
	end subroutine

end module
