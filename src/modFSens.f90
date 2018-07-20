module modFSens

	use modBC
	use modRecord
    use modAssign
    use modInputHelper

	implicit none

    integer, parameter :: COLLOCATED = 492, NONCOLLOCATED = 493, INJECTION = 494

	type, extends(PM1D) :: FSens
		real(mp) :: Lv, dv
		integer :: ngv
		real(mp), allocatable :: J(:,:), n_A(:,:), Dvf(:,:)     !For source term
		real(mp), allocatable :: f_A(:,:)                       !For remapping procedure

        integer :: scheme                                       !COLLOCATED, NONCOLLOCATED, INJECTION        

		integer :: NInject										!Number of injecting particle per direction
		real(mp), allocatable :: xp_inject(:), vp_inject(:,:), frac_inject(:,:,:)
		integer, allocatable :: g_inject(:,:), gv_inject(:,:)

		integer :: NLimit											!Number limit of particle for redistribution
		real(mp), allocatable :: xp_remesh(:), vp_remesh(:,:), frac_remesh(:,:,:)
		integer, allocatable :: g_remesh(:,:), gv_remesh(:,:)
	contains
		procedure, pass(dpm) :: buildFSens
		procedure, pass(this) :: destroyFSens
		procedure, pass(this) :: FDistribution
		procedure, pass(this) :: FVelocityGradient
		procedure, pass(dpm) :: FSensSourceTerm
		procedure, pass(this) :: InjectSource
		procedure, pass(this) :: Redistribute
		procedure, pass(this) :: numberDensity
		procedure, pass(this) :: updateWeight
	end type

contains

	subroutine buildFSens(dpm,pm,Lv,Ngv,NInject,NLimit)
		class(FSens), intent(out) :: dpm
		type(PM1D), intent(in) :: pm
		real(mp), intent(in) :: Lv
		integer, intent(in) :: Ngv, NInject, NLimit
		integer :: i
		real(mp), allocatable :: temp(:)

		call dpm%buildPM1D(pm%nt*pm%dt,pm%ni*pm%dt,pm%ng,pm%n,	&
							pm%pBCindex,pm%mBCindex,pm%a(1)%order,	&
							pm%dt,pm%L,pm%A0,pm%eps0)
		do i=1,pm%n
			call dpm%p(i)%buildSpecies(pm%p(i)%qs,pm%p(i)%ms)
		end do
		call dpm%m%setMesh(pm%m%rho_back)

		dpm%NInject = INT(SQRT(1.0_mp*NInject))**2
		dpm%NLimit = INT(SQRT(1.0_mp*NLimit))**2
		dpm%Lv = Lv
		dpm%ngv = Ngv
		dpm%dv = Lv/Ngv
		allocate(dpm%J(pm%ng,2*Ngv+1))
		allocate(dpm%f_A(pm%ng,2*Ngv+1))
		allocate(dpm%n_A(pm%ng,2*Ngv+1))
		allocate(dpm%Dvf(pm%ng,2*Ngv+1))
		dpm%J = 0.0_mp
		dpm%n_A = 0.0_mp
		dpm%f_A = 0.0_mp
		dpm%Dvf = 0.0_mp

        select case( getOption('sensitivity_pdf/discretization','non-collocated') )
            case('collocated')
                dpm%scheme = COLLOCATED
            case('non-collocated')
                dpm%scheme = NONCOLLOCATED
            case('injection')
                dpm%scheme = INJECTION
        end select

        if( dpm%scheme.eq.INJECTION ) then
    		call createDistribution(dpm,dpm%NInject,dpm%xp_inject,dpm%vp_inject,	&
	    									dpm%g_inject,dpm%gv_inject,dpm%frac_inject)
	    	call createDistribution(dpm,dpm%NLimit,dpm%xp_remesh,dpm%vp_remesh,	    &
		    								dpm%g_remesh,dpm%gv_remesh,dpm%frac_remesh)
        end if
	end subroutine

	subroutine destroyFSens(this)
		class(FSens), intent(inout) :: this
		integer :: i

		call destroyPM1D(this)
		deallocate(this%J)
		deallocate(this%n_A)
		deallocate(this%f_A)
		deallocate(this%Dvf)
        if( this%scheme.eq.INJECTION ) then
    		deallocate(this%xp_inject)
    		deallocate(this%vp_inject)
    		deallocate(this%g_inject)
    		deallocate(this%gv_inject)
    		deallocate(this%frac_inject)
    		deallocate(this%xp_remesh)
    		deallocate(this%vp_remesh)
    		deallocate(this%g_remesh)
    		deallocate(this%gv_remesh)
    		deallocate(this%frac_remesh)
        end if
	end subroutine

	subroutine FSensSourceTerm(dpm,qs,ms,nk,fsr)							!Take v-derivative, Multiply dE
		class(FSens), intent(inout) :: dpm
		real(mp), intent(in) :: qs,ms
		integer, intent(in), optional :: nk
		type(recordData), intent(in), optional :: fsr
		integer :: ng,ngv,kr,i
		character(len=100) :: kstr
		real(mp) :: dx,dv
		real(mp) :: E_temp(dpm%m%ng)
		real(mp) :: Dvf(dpm%m%ng,2*dpm%ngv+1)
		ng = dpm%m%ng
		ngv = dpm%ngv
		dx = dpm%m%dx
		dv = dpm%dv

		!Multiply E_A
		E_temp = qs/ms*dpm%m%E
		do i=1,2*ngv+1
			dpm%J(:,i) = dpm%Dvf(:,i)*E_temp 
        end do

		!Multiply dt
		dpm%J = -dpm%J*dpm%dt
	end subroutine

	subroutine InjectSource(this,p,J)
		class(FSens), intent(inout) :: this
		type(species), intent(inout) :: p
		real(mp), intent(in), dimension(this%m%ng,2*this%ngv+1) :: J
		real(mp) :: spwt0(this%NInject)
		logical, dimension(this%NInject) :: IsNonTrivial
		real(mp) :: maxspwt
		integer :: N0
		real(mp), allocatable :: v0(:,:)

		call SourceAssign(J,this%L,this%Lv,this%g_inject,this%gv_inject,this%frac_inject,spwt0)
		maxspwt = MAXVAL(ABS(spwt0))
		IsNonTrivial = (ABS(spwt0).ge.maxspwt*1e-8)
		N0 = COUNT(IsNonTrivial)
		allocate(v0(N0,3))
		v0 = 0.0_mp
		v0(:,1) = PACK(this%vp_inject(:,1),IsNonTrivial)
		call p%appendSpecies(N0,PACK(this%xp_inject,IsNonTrivial),	&
									v0,PACK(spwt0,IsNonTrivial))
		deallocate(v0)
	end subroutine

	subroutine SourceAssign(J,Lx,Lv,g,gv,frac,spwt0)
		real(mp), intent(in) :: J(:,:)
		real(mp), intent(in) :: Lx,Lv
		integer, intent(in) :: g(:,:), gv(:,:)
		real(mp), intent(in) :: frac(:,:,:)
		real(mp), intent(out) :: spwt0(size(g,2))
		integer :: k,N
		N = size(g,2)

		do k=1,N
			spwt0(k) = SUM( J(g(:,k),gv(:,k))*frac(:,:,k) )*Lx*2.0_mp*Lv/N
		end do
	end subroutine

	subroutine createDistribution(this,N,xp0,vp0,g,gv,frac)
		class(FSens), intent(inout) :: this
		integer, intent(inout) :: N
		real(mp), intent(out), allocatable :: xp0(:), vp0(:,:)
		integer, intent(out), allocatable :: g(:,:), gv(:,:)
		real(mp), intent(out), allocatable :: frac(:,:,:)
		integer :: Nx, newN
		integer :: i,i1,i2,k,Np,nk
		integer :: vgl, vgr
		real(mp) :: w, h, dx,dv
		real(mp), allocatable :: fx(:,:)
		dx = this%m%dx
		dv = this%dv
		Nx = INT(SQRT(N*1.0_mp))
		newN = Nx*Nx
		N = newN

		allocate(xp0(newN))
		allocate(vp0(newN,3))
		allocate(g(this%a(1)%order+1,newN))
		allocate(gv(2,newN))
		allocate(frac(this%a(1)%order+1,2,newN))
		allocate(fx(this%a(1)%order+1,newN))
	
		!Uniform grid distribution on phase space
		xp0 = 0.0_mp
		vp0 = 0.0_mp
		do i2=1,Nx
			do i1=1,Nx
				xp0(i1+Nx*(i2-1)) = (i1-0.5_mp)*this%L/Nx
				vp0(i1+Nx*(i2-1),:) = (i2-0.5_mp)*2.0_mp*this%Lv/Nx - 1.0_mp*this%Lv
			end do
		end do

		do k=1,newN
			!X-direction Interpolation
			CALL this%a(1)%assignMatrix(xp0(k),dx,g(:,k),fx(:,k))
			CALL this%a(1)%adjustGrid(this%m%ng,g(:,k),fx(:,k))

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
		class(FSens), intent(inout) :: this
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
			this%f_A(g,vgl) = this%f_A(g,vgl) + (1.0_mp-h)*p%spwt(k)*a%frac(:,k)/this%dv/this%m%dx
			this%f_A(g,vgr) = this%f_A(g,vgr) + h*p%spwt(k)*a%frac(:,k)/this%dv/this%m%dx
		end do
        if( a%mBCidx .ne. 0 ) then
            this%f_A((/1,a%ng/),:) = this%f_A((/1,a%ng/),:)*2.0_mp
        end if
		this%f_A(:,1) = this%f_A(:,1)*2.0_mp
		this%f_A(:,2*this%ngv+1) = this%f_A(:,2*this%ngv+1)*2.0_mp
	end subroutine

	subroutine FVelocityGradient(this,p,a)
		class(FSens), intent(inout) :: this
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
			    this%Dvf(g,gv(i)) = this%Dvf(g,gv(i)) + p%spwt(k)*a%frac(:,k)*fracv(i)/this%dv/this%m%dx
            end do
		end do
        if( a%mBCidx .ne. 0 ) then
            this%Dvf((/1,a%ng/),:) = this%Dvf((/1,a%ng/),:)*2.0_mp
        end if
		this%Dvf(:,1) = this%Dvf(:,1)*2.0_mp
		this%Dvf(:,2*this%ngv+1) = this%Dvf(:,2*this%ngv+1)*2.0_mp
	end subroutine

	subroutine numberDensity(this,p,a)
		class(FSens), intent(inout) :: this
		type(species), intent(in) :: p
		type(pmAssign), intent(in) :: a
		real(mp), dimension(this%m%ng,2*this%ngv+3) :: n_temp
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
			n_temp(g,g_v) = n_temp(g,g_v) + frac/this%m%dx/this%dv
		end do
        if( a%mBCidx .ne. 0 ) then
            n_temp((/1,a%ng/),:) = n_temp((/1,a%ng/),:)*2.0_mp
        end if
        this%n_A = n_temp(:,2:2*this%ngv+2)
	end subroutine

	subroutine updateWeight(this,p,a)
		class(FSens), intent(inout) :: this
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

	subroutine FDistribution_sync(this,p,a)
		class(FSens), intent(inout) :: this
		type(species), intent(in) :: p
		type(pmAssign), intent(in) :: a
		real(mp), dimension(this%m%ng,2*this%ngv+3) :: n_temp, f_temp
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
			n_temp(g,gv) = n_temp(g,gv) + frac/this%m%dx/this%dv
			f_temp(g,gv) = f_temp(g,gv) + p%spwt(k)*frac/this%m%dx/this%dv
		end do
        if( a%mBCidx .ne. 0 ) then
            n_temp((/1,a%ng/),:) = n_temp((/1,a%ng/),:)*2.0_mp
            f_temp((/1,a%ng/),:) = f_temp((/1,a%ng/),:)*2.0_mp
        end if
		this%n_A = n_temp(:,2:2*this%ngv+2)
		this%f_A = f_temp(:,2:2*this%ngv+2)
	end subroutine

	subroutine Redistribute(this,p)  !periodic BC only!!
		class(FSens), intent(inout) :: this
		type(species), intent(inout) :: p
		real(mp), dimension(this%NLimit) :: spwt0
		logical, dimension(this%NLimit) :: IsNonTrivial
		real(mp) :: xp, vp, dx, dv, Lv
		integer :: k, k1, k2, wk, nx, gx(4),gv(4)
		real(mp) :: fx(4), fv(4), hx, hv
		real(mp) :: maxspwt
		integer :: N0
		real(mp), allocatable :: v0(:,:)
		spwt0 = 0.0_mp
		nx = INT(SQRT(1.0_mp*this%NLimit))
		dx = this%m%L/nx
		dv = 2.0_mp*this%Lv/nx
		Lv = this%Lv

		if( p%np>3*this%NLimit ) then
			do k=1,p%np
				xp=p%xp(k)
				vp=p%vp(k,1)

				gx = FLOOR(xp/dx-0.5_mp) + (/0,1,2,3/)
				gv = FLOOR((vp+Lv)/dv-0.5_mp) + (/0,1,2,3/)

				hx = xp/dx-0.5_mp - FLOOR(xp/dx-0.5_mp)
				hv = (vp+Lv)/dv-0.5_mp - FLOOR((vp+Lv)/dv-0.5_mp)
				fx = ABS((/hx+1.0_mp,hx,1.0_mp-hx,2.0_mp-hx/))
				fv = ABS((/hv+1.0_mp,hv,1.0_mp-hv,2.0_mp-hv/))
				fx(2:3) = 1.0_mp - 2.5_mp*fx(2:3)**2 + 1.5_mp*fx(2:3)**3
				fx((/1,4/)) = 0.5_mp*(2.0_mp-fx((/1,4/)))**2*(1.0_mp-fx((/1,4/)))
				fv(2:3) = 1.0_mp - 2.5_mp*fv(2:3)**2 + 1.5_mp*fv(2:3)**3
				fv((/1,4/)) = 0.5_mp*(2.0_mp-fv((/1,4/)))**2*(1.0_mp-fv((/1,4/)))

				where( gx<1 )
					gx = gx+nx
				elsewhere( gx>nx )
					gx = gx-nx
				end where
				where( gv<1 )
					gv = 1
					fv = 0.0_mp
				elsewhere( gv>nx )
					gv = nx
					fv = 0.0_mp
				end where

				do k2=1,4
					do k1=1,4
						wk = gx(k1) + nx*(gv(k2)-1)
						spwt0(wk) = spwt0(wk) + p%spwt(k)*fx(k1)*fv(k2)
					end do
				end do
			end do
			maxspwt = MAXVAL(abs(spwt0))
			IsNonTrivial = (ABS(spwt0).ge.maxspwt*1e-6)
			N0 = COUNT(IsNonTrivial)
!print *, maxspwt, MINVAL(PACK(ABS(spwt0),IsNonTrivial)),MINVAL(ABS(spwt0))
!print *, N0, this%NLimit
			allocate(v0(N0,3))
			v0=0.0_mp
			v0(:,1) = PACK(this%vp_remesh(:,1),IsNonTrivial)

!			call p%setSpecies(this%NLimit,this%xp_remesh,this%vp_remesh,spwt0)
			call p%setSpecies(N0,PACK(this%xp_remesh,IsNonTrivial),	&
										v0,	&
										PACK(spwt0,IsNonTrivial))
			deallocate(v0)
		end if
	end subroutine

end module
