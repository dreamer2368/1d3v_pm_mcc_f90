module modPhaseSpaceMesh

    use modSpecies
    use modAssign

	implicit none

    type phaseSpaceMesh
        real(mp) :: Lx, Lv, dx, dv
        integer :: ng, ngv

        real(mp), allocatable :: J(:,:), n_A(:,:), Dvf(:,:)
        real(mp), allocatable :: f_A(:,:)

        type(pmAssign) :: a
    contains
        procedure, pass(this) :: buildPhaseSpaceMesh
        procedure, pass(this) :: destroyPhaseSpaceMesh
        procedure, pass(this) :: FSensSourceTerm
        procedure, pass(this) :: createDistribution
        procedure, pass(this) :: FDistribution
        procedure, pass(this) :: FVelocityGradient
        procedure, pass(this) :: numberDensity
        procedure, pass(this) :: updateWeight
    end type

contains

	subroutine buildPhaseSpaceMesh(this,Lx0,Lv0,Ng,Ngv,dx)
		class(phaseSpaceMesh), intent(inout) :: this
		real(mp), intent(in) :: Lx0, Lv0, dx
		integer, intent(in) :: Ng, Ngv

        this%Lx = Lx0
		this%Lv = Lv0
        this%ng = Ng
		this%ngv = Ngv
        this%dx = dx
		this%dv = Lv0/Ngv
		allocate(this%J(Ng,2*Ngv+1))
		allocate(this%f_A(Ng,2*Ngv+1))
		allocate(this%n_A(Ng,2*Ngv+1))
		allocate(this%Dvf(Ng,2*Ngv+1))
		this%J = 0.0_mp
		this%n_A = 0.0_mp
		this%f_A = 0.0_mp
		this%Dvf = 0.0_mp
        call this%a%buildAssign(2*Ngv+1,1,1)
	end subroutine

	subroutine destroyPhaseSpaceMesh(this)
		class(phaseSpaceMesh), intent(inout) :: this
		integer :: i

		deallocate(this%J)
		deallocate(this%n_A)
		deallocate(this%f_A)
		deallocate(this%Dvf)

        call this%a%destroyAssign
	end subroutine

	subroutine FSensSourceTerm(this,qs,ms,E,dt)			    		!Multiply qs/ms*dE*dt
        class(phaseSpaceMesh), intent(inout) :: this
		real(mp), intent(in) :: qs,ms,dt
        real(mp), intent(in) :: E(this%ng)
		integer :: i,j

		!Multiply E_A
	    do i=1,2*this%ngv+1
	    	this%J(:,i) = this%Dvf(:,i)*E
        end do

		!Multiply dt
		this%J = -this%J*qs/ms*dt
	end subroutine

	subroutine createDistribution(this,a,N,xp0,vp0,g,gv,frac)
		class(phaseSpaceMesh), intent(inout) :: this
        type(pmAssign), intent(inout) :: a
		integer, intent(inout) :: N
		real(mp), intent(out), allocatable :: xp0(:), vp0(:,:)
		integer, intent(out), allocatable :: g(:,:), gv(:,:)
		real(mp), intent(out), allocatable :: frac(:,:,:)
		integer :: Nx, newN
		integer :: i,i1,i2,k,Np,nk
		integer :: vgl, vgr
		real(mp) :: w, h, dx,dv
		real(mp), allocatable :: fx(:,:)
		dx = this%dx
		dv = this%dv
		Nx = INT(SQRT(N*1.0_mp))
		newN = Nx*Nx
		N = newN

		allocate(xp0(newN))
		allocate(vp0(newN,3))
		allocate(g(a%order+1,newN))
		allocate(gv(2,newN))
		allocate(frac(a%order+1,2,newN))
		allocate(fx(a%order+1,newN))
	
		!Uniform grid distribution on phase space
		xp0 = 0.0_mp
		vp0 = 0.0_mp
		do i2=1,Nx
			do i1=1,Nx
				xp0(i1+Nx*(i2-1)) = (i1-0.5_mp)*this%Lx/Nx
				vp0(i1+Nx*(i2-1),:) = (i2-0.5_mp)*2.0_mp*this%Lv/Nx - 1.0_mp*this%Lv
			end do
		end do

		do k=1,newN
			!X-direction Interpolation
			CALL a%assignMatrix(xp0(k),dx,g(:,k),fx(:,k))
			CALL a%adjustGrid(this%ng,g(:,k),fx(:,k))

			!V-direction Interpolation
			vgl = FLOOR(vp0(k,1)/this%dv) + this%ngv+1
			vgr = vgl+1

			if( vgl<1 .or. vgr>2*this%ngv+1 ) then
				gv(:,k) = this%ngv
				frac(:,:,k) = 0.0_mp
				cycle
			end if

			gv(:,k) = (/vgl,vgr/)
			h = vp0(k,1)/dv - FLOOR(vp0(k,1)/dv)
			frac(:,1,k) = fx(:,k)*(1.0_mp-h)
			frac(:,2,k) = fx(:,k)*h
		end do
	end subroutine

	subroutine FDistribution(this,p,a)
		class(phaseSpaceMesh), intent(inout) :: this
		type(species), intent(in) :: p
		type(pmAssign), intent(in) :: a
		integer :: i,k, g(a%order+1)
		integer :: vgl, vgr
		real(mp) :: vp, h

		!F to phase space
		this%f_A = 0.0_mp
		do k = 1, p%np
			vp = p%vp(k,1)
			vgl = FLOOR(vp/this%dv) + this%ngv+1
			vgr = vgl+1
			if( vgl<1 .or. vgr>2*this%ngv+1 )	cycle
			g = a%g(:,k)
			h = vp/this%dv - FLOOR(vp/this%dv)
			this%f_A(g,vgl) = this%f_A(g,vgl) + (1.0_mp-h)*p%spwt(k)*a%frac(:,k)/this%dv/this%dx
			this%f_A(g,vgr) = this%f_A(g,vgr) + h*p%spwt(k)*a%frac(:,k)/this%dv/this%dx
		end do
        if( a%mBCidx .ne. 0 ) then
            this%f_A((/1,a%ng/),:) = this%f_A((/1,a%ng/),:)*2.0_mp
        end if
		this%f_A(:,1) = this%f_A(:,1)*2.0_mp
		this%f_A(:,2*this%ngv+1) = this%f_A(:,2*this%ngv+1)*2.0_mp
	end subroutine

	subroutine FDistribution_sync(this,p,a)
		class(phaseSpaceMesh), intent(inout) :: this
		type(species), intent(in) :: p
		type(pmAssign), intent(in) :: a
		real(mp), dimension(this%ng,2*this%ngv+3) :: n_temp, f_temp
		integer :: i,k, g(a%order+1),gv(2)
		real(mp) :: frac(a%order+1,2)
		integer :: vgl, vgr
		real(mp) :: vp, h

		!F to phase space
		n_temp = 0.0_mp
		f_temp = 0.0_mp
		do k = 1, p%np
			vp = p%vp(k,1)
			vgl = FLOOR(vp/this%dv) + this%ngv+2
			vgr = vgl+1
			if( vgl<1 .or. vgr>2*this%ngv+3 )	cycle
			g = a%g(:,k)
			gv = (/vgl,vgr/)
			h = vp/this%dv - FLOOR(vp/this%dv)
			frac(:,1) = (1.0_mp-h)*a%frac(:,k)
			frac(:,2) = h*a%frac(:,k)
			n_temp(g,gv) = n_temp(g,gv) + frac/this%dx/this%dv
			f_temp(g,gv) = f_temp(g,gv) + p%spwt(k)*frac/this%dx/this%dv
		end do
        if( a%mBCidx .ne. 0 ) then
            n_temp((/1,a%ng/),:) = n_temp((/1,a%ng/),:)*2.0_mp
            f_temp((/1,a%ng/),:) = f_temp((/1,a%ng/),:)*2.0_mp
        end if
		this%n_A = n_temp(:,2:2*this%ngv+2)
		this%f_A = f_temp(:,2:2*this%ngv+2)
	end subroutine

	subroutine FVelocityGradient(this,p,a)
		class(phaseSpaceMesh), intent(inout) :: this
		type(species), intent(in) :: p
		type(pmAssign), intent(in) :: a
		integer :: i,k, g(a%order+1), gv(3)
		real(mp) :: vp, fracv(3)

		!DvF to phase space
		this%Dvf = 0.0_mp
		do k = 1, p%np
			vp = p%vp(k,1)
            call assign_TSC_derivative(vp,this%dv,gv,fracv)	!for velocity derivative interpolation
            where( abs(gv)>this%ngv )
                gv = 0
                fracv = 0.0_mp
            elsewhere
                gv = gv + this%ngv + 1
            end where
			g = a%g(:,k)

            do i=1,3
			    this%Dvf(g,gv(i)) = this%Dvf(g,gv(i)) + p%spwt(k)*a%frac(:,k)*fracv(i)/this%dv/this%dx
            end do
		end do
        if( a%mBCidx .ne. 0 ) then
            this%Dvf((/1,a%ng/),:) = this%Dvf((/1,a%ng/),:)*2.0_mp
        end if
		this%Dvf(:,1) = this%Dvf(:,1)*2.0_mp
		this%Dvf(:,2*this%ngv+1) = this%Dvf(:,2*this%ngv+1)*2.0_mp
	end subroutine

	subroutine numberDensity(this,p,a)
		class(phaseSpaceMesh), intent(inout) :: this
		type(species), intent(in) :: p
		type(pmAssign), intent(in) :: a
		real(mp), dimension(this%ng,2*this%ngv+3) :: n_temp
		integer :: i,k, g(a%order+1), g_v(2)
		real(mp) :: frac(2,2)
		integer :: vgl, vgr, k_temp
		real(mp) :: vp, h
		integer, dimension(:,:), pointer :: g_x
		real(mp), dimension(:,:), pointer :: frac_x

		!F_A to phase space
		n_temp = 0.0_mp
		g_x=>a%g
		frac_x=>a%frac
		do k = 1, p%np
			vp = p%vp(k,1)
			vgl = FLOOR(vp/this%dv) + this%ngv+2
			vgr = vgl+1
			if( vgl<1 .or. vgr>2*this%ngv+3 )	then
				cycle
			end if
			g = g_x(:,k)
			g_v = (/ vgl, vgr /)
			h = vp/this%dv - FLOOR(vp/this%dv)
			frac(:,1) = (1.0_mp-h)*frac_x(:,k)
			frac(:,2) = h*frac_x(:,k)
			n_temp(g,g_v) = n_temp(g,g_v) + frac/this%dx/this%dv
		end do
        if( a%mBCidx .ne. 0 ) then
            n_temp((/1,a%ng/),:) = n_temp((/1,a%ng/),:)*2.0_mp
        end if
        this%n_A = n_temp(:,2:2*this%ngv+2)
	end subroutine

	subroutine updateWeight(this,p,a)
		class(phaseSpaceMesh), intent(inout) :: this
		type(species), intent(inout) :: p
		type(pmAssign), intent(in) :: a
		integer :: i,k, g(a%order+1), g_v(2)
		real(mp) :: frac(2,2)
		integer :: vgl, vgr, k_temp
		real(mp) :: vp, h
		integer, dimension(:,:), pointer :: g_x
		real(mp), dimension(:,:), pointer :: frac_x

		!Update weight
		g_x=>a%g
		frac_x=>a%frac
		do k = 1, p%np
			vp = p%vp(k,1)
			vgl = FLOOR(vp/this%dv) + this%ngv+1
			vgr = vgl+1
			if( vgl<1 .or. vgr>2*this%ngv+1 )	then
				cycle
			end if
			g = g_x(:,k)
			g_v = (/ vgl, vgr /)
			h = vp/this%dv - FLOOR(vp/this%dv)
			frac(:,1) = (1.0_mp-h)*frac_x(:,k)
			frac(:,2) = h*frac_x(:,k)
			p%spwt(k) = p%spwt(k) + SUM( this%J(g,g_v)*frac/this%n_A(g,g_v) )
		end do
	end subroutine

end module
