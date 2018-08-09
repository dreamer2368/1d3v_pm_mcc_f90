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
        call this%a%buildAssign(2*Ngv+1,1,3)
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
		integer :: i,k, g(a%order+1), gv(2)
		real(mp) :: vp, fracv(2)

		!DvF to phase space
		this%Dvf = 0.0_mp
		do k = 1, p%np
			vp = p%vp(k,1)
            call assign_CIC_derivative(vp,this%dv,gv,fracv)	!for velocity derivative interpolation
!            call assign_TSC_derivative(vp,this%dv,gv,fracv)	!for velocity derivative interpolation
            if( gv(1)<-this%ngv .or. gv(2)>this%ngv ) cycle
            gv = gv + this%ngv+1
			g = a%g(:,k)

            do i=1,2
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
		real(mp), dimension(this%ng,2*this%ngv+1) :: n_temp
		integer :: i,k, g(a%order+1), g_vp(2)
		real(mp) :: frac_vp(2), frac(2,2)
		integer :: vgl, vgr, k_temp
		real(mp) :: vp, h
		integer, dimension(:,:), pointer :: g_x, g_v
		real(mp), dimension(:,:), pointer :: frac_x, frac_v

		!F_A to phase space
		g_x=>a%g
		frac_x=>a%frac

		DEALLOCATE(this%a%g)
		DEALLOCATE(this%a%frac)
		this%a%np = p%np
		ALLOCATE(this%a%g(this%a%order+1,p%np))
		ALLOCATE(this%a%frac(this%a%order+1,p%np))
		g_v=>this%a%g
		frac_v=>this%a%frac

		n_temp = 0.0_mp
		do k = 1, p%np
			vp = p%vp(k,1)
            call this%a%assignMatrix(vp,this%dv,g_vp,frac_vp)
            g_vp = g_vp + this%ngv
!            call this%a%adjustGrid(this%ngv,g_vp,frac_vp)

            if( g_vp(1)<1 .or. g_vp(2)>2*this%ngv+1 ) then
                g_v(:,k) = this%ngv+1
                frac_v(:,k) = 0.0_mp
                cycle
            else
                g_v(:,k) = g_vp
                frac_v(:,k) = frac_vp
            end if
                
			g = g_x(:,k)
!            g_v(:,k) = g_vp
!            frac_v(:,k) = frac_vp
            do i=1,this%a%order+1
			    n_temp(g,g_vp(i)) = n_temp(g,g_vp(i)) + frac_vp(i)*frac_x(:,k)/this%dx/this%dv
            end do
		end do
        if( a%mBCidx .ne. 0 ) then
            n_temp((/1,a%ng/),:) = n_temp((/1,a%ng/),:)*2.0_mp
        end if
        n_temp(:,(/1,2*this%ngv+1/)) = n_temp(:,(/1,2*this%ngv+1/))*2.0_mp
        
        this%n_A = n_temp
	end subroutine

	subroutine updateWeight(this,p,a)
		class(phaseSpaceMesh), intent(inout) :: this
		type(species), intent(inout) :: p
		type(pmAssign), intent(in) :: a
		integer :: i,j,k, g(a%order+1), g_vp(2)
		real(mp) :: frac(2,2), frac_vp(2)
		integer :: vgl, vgr, k_temp
		real(mp) :: vp, h
		integer, dimension(:,:), pointer :: g_x, g_v
		real(mp), dimension(:,:), pointer :: frac_x, frac_v
        real(mp), dimension(size(this%J,1),size(this%J,2)) :: H_temp
        H_temp = this%J/this%n_A

		!Update weight
		g_x=>a%g
		frac_x=>a%frac
		g_v=>this%a%g
		frac_v=>this%a%frac

		do k = 1, p%np
            if( SUM(frac_v(:,k)).eq.0.0_mp ) cycle
			g = g_x(:,k)
            g_vp = g_v(:,k)
!			vp = p%vp(k,1)
!            call this%a%assignMatrix(vp,this%dv,g_vp,frac_vp)
!            g_vp = g_vp + this%ngv
!            call this%a%adjustGrid(this%ngv,g_vp,frac_vp)

            do i=1,this%a%order+1
!                p%spwt(k) = p%spwt(k) + SUM( H_temp(g,g_vp(i))*frac_vp(i)*frac_x(:,k) )
                p%spwt(k) = p%spwt(k) + SUM( H_temp(g,g_vp(i))*frac_v(i,k)*frac_x(:,k) )
            end do
		end do
	end subroutine

end module
