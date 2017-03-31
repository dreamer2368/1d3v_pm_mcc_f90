module modAssign

	use modGridBC

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
		procedure, pass(this) :: AssignPhaseSpace
		procedure, pass(this) :: NumberDensity
		procedure, pass(this) :: PhaseSpaceDensity
	end type

	abstract interface
		subroutine assignMatrix(xp,dx,g,frac)
			use constants
			real(mp), intent(in) :: xp, dx
			integer, intent(out) :: g(:)
			real(mp), intent(out) :: frac(:)
		end subroutine
	end interface

!	abstract interface
!		subroutine assignV(vp,dv,gv,fracv)
!			use constants
!			real(mp), intent(in) :: vp, dv
!			integer, intent(out) :: gv(:)
!			real(mp), intent(out) :: fracv(:)
!		end subroutine
!	end interface

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
		g1 = FLOOR(xp/dx - 0.5_mp)+1
		gl = g1
		gr = gl+1

		h = xp/dx - g1 + 0.5_mp
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
		integer :: i

		this%np = p%np
		do i=1,p%np
			CALL this%assignMatrix(p%xp(i),m%dx,g,frac)
			CALL this%adjustGrid(m%ng,g,frac)
			m%rho( g ) = m%rho( g ) + p%spwt(i)*p%qs/m%dx*frac
		end do
	end subroutine

	subroutine forceAssign(this,p,m)
		class(pmAssign), intent(inout) :: this
		type(species), intent(inout) :: p
		type(mesh), intent(in) :: m
		integer :: g(this%order+1)
		real(mp) :: frac(this%order+1)
		integer :: i

		p%Ep = 0.0_mp
		do i=1,this%np
			CALL this%assignMatrix(p%xp(i),m%dx,g,frac)
			CALL this%adjustGrid(m%ng,g,frac)
			p%Ep(i) = sum( frac*m%E( g ) )
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

!=======================  Phase-Space Assignment  ======================================

	subroutine AssignPhaseSpace(this,xp,vp,m,g,gv,frac)
		class(pmAssign), intent(in) :: this
		real(mp), intent(in) :: xp, vp
		type(mesh), intent(in) :: m
		integer, intent(out) :: g(this%order+1), gv(2)
		real(mp), dimension(this%order+1,2), intent(out) :: frac
		integer :: vgl,vgr
		real(mp) :: f_x(this%order+1), hv

		vgl = FLOOR(vp/m%dv) + m%ngv+1
		vgr = vgl+1
		if( vgl<0 .or. vgr>2*m%ngv+2 ) then
			g = 2
			gv = m%ngv
			frac = 0.0_mp
			RETURN
		end if
		hv = vp/m%dv - FLOOR(vp/m%dv)

		call this%assignMatrix(xp,m%dx,g,f_x)
		call this%adjustGrid(m%ng,g,f_x)

		if( vgl.eq.0 ) then
			gv(1) = m%ngv
			gv(2) = vgr
			frac(:,1) = 0.0_mp
			frac(:,2) = hv*f_x
		elseif( vgr.eq.2*m%ngv+2 ) then
			gv(1) = vgl
			gv(2) = m%ngv
			frac(:,1) = (1.0_mp-hv)*f_x
			frac(:,2) = 0.0_mp
		else
			gv = (/ vgl,vgr /)
			frac(:,1) = (1.0_mp-hv)*f_x
			frac(:,2) = hv*f_x
		end if
	end subroutine

	subroutine NumberDensity(this,p,m)
		class(pmAssign), intent(in) :: this
		type(species), intent(in) :: p
		type(mesh), intent(inout) :: m
		integer :: k, g(this%order+1), gv(2)
		real(mp) :: frac(this%order+1,2)
		real(mp) :: xp, vp
		integer :: vgl,vgr
		real(mp) :: f_x(this%order+1), hv

		!N_A to phase space
		do k = 1, p%np
			xp = p%xp(k)
			vp = p%vp(k,1)
			vgl = FLOOR(vp/m%dv) + m%ngv+1
			vgr = vgl+1
			if( vgl<0 .or. vgr>2*m%ngv+2 ) then
!				g = 2
!				gv = m%ngv
!				frac = 0.0_mp
				CYCLE
			end if
			hv = vp/m%dv - FLOOR(vp/m%dv)
	
			call this%assignMatrix(xp,m%dx,g,f_x)
			call this%adjustGrid(m%ng,g,f_x)
	
			if( vgl.eq.0 ) then
				gv(1) = m%ngv
				gv(2) = vgr
				frac(:,1) = 0.0_mp
				frac(:,2) = hv*f_x
			elseif( vgr.eq.2*m%ngv+2 ) then
				gv(1) = vgl
				gv(2) = m%ngv
				frac(:,1) = (1.0_mp-hv)*f_x
				frac(:,2) = 0.0_mp
			else
				gv = (/ vgl,vgr /)
				frac(:,1) = (1.0_mp-hv)*f_x
				frac(:,2) = hv*f_x
			end if
!			call this%AssignPhaseSpace(p%xp(k),p%vp(k,1),m,g,gv,frac)
			m%f(g,gv) = m%f(g,gv) + frac/m%dx/m%dv
		end do
	end subroutine

	subroutine PhaseSpaceDensity(this,p,m)
		class(pmAssign), intent(in) :: this
		type(species), intent(in) :: p
		type(mesh), intent(inout) :: m
		integer :: k, g(this%order+1),gv(2)
		real(mp) :: frac(this%order+1,2)
		real(mp) :: xp, vp
		integer :: vgl,vgr
		real(mp) :: f_x(this%order+1), hv

		!N_A to phase space
		do k = 1, p%np
			xp = p%xp(k)
			vp = p%vp(k,1)
			vgl = FLOOR(vp/m%dv) + m%ngv+1
			vgr = vgl+1
			if( vgl<0 .or. vgr>2*m%ngv+2 ) then
!				g = 2
!				gv = m%ngv
!				frac = 0.0_mp
				CYCLE
			end if
			hv = vp/m%dv - FLOOR(vp/m%dv)
	
			call this%assignMatrix(xp,m%dx,g,f_x)
			call this%adjustGrid(m%ng,g,f_x)
	
			if( vgl.eq.0 ) then
				gv(1) = m%ngv
				gv(2) = vgr
				frac(:,1) = 0.0_mp
				frac(:,2) = hv*f_x
			elseif( vgr.eq.2*m%ngv+2 ) then
				gv(1) = vgl
				gv(2) = m%ngv
				frac(:,1) = (1.0_mp-hv)*f_x
				frac(:,2) = 0.0_mp
			else
				gv = (/ vgl,vgr /)
				frac(:,1) = (1.0_mp-hv)*f_x
				frac(:,2) = hv*f_x
			end if
!			call this%AssignPhaseSpace(p%xp(k),p%vp(k,1),m,g,gv,frac)
			m%f(g,gv) = m%f(g,gv) + frac*p%spwt(k)/m%dx/m%dv
		end do
	end subroutine

end module
