module modFSens

	use modBC

	implicit none

	type FSens
		type(PM1D) :: dpm
		real(mp) :: Lv, dv
		integer :: ngv
		real(mp), allocatable :: j(:,:), f_A(:,:)
		integer :: NInject										!Number of injecting particle per direction
		integer :: NLimit											!Number limit of particle for redistribution
	contains
		procedure, pass(this) :: buildFSens
		procedure, pass(this) :: destroyFSens
		procedure, pass(this) :: FSensSourceTerm
		procedure, pass(this) :: InjectSource
		procedure, pass(this) :: FSensDistribution
		procedure, pass(this) :: Redistribute
	end type

contains

	subroutine buildFSens(this,pm,Lv,Ngv,NInject,NLimit)
		class(FSens), intent(out) :: this
		type(PM1D), intent(in) :: pm
		real(mp), intent(in) :: Lv
		integer, intent(in) :: Ngv, NInject, NLimit
		integer :: i

		call this%dpm%buildPM1D(pm%nt*pm%dt,pm%ni*pm%dt,pm%ng,pm%n,	&
							pm%pBCindex,pm%mBCindex,pm%a(1)%order,	&
							pm%dt,pm%L,pm%A0,pm%eps0)
		do i=1,pm%n
			call this%dpm%p(i)%buildSpecies(pm%p(i)%qs,pm%p(i)%ms)
		end do
		call this%dpm%m%setMesh(pm%m%rho_back)

		this%NInject = NInject
		this%NLimit = NLimit
		this%Lv = Lv
		this%ngv = Ngv
		this%dv = Lv/Ngv
		allocate(this%j(pm%ng,2*Ngv+1))
		allocate(this%f_A(pm%ng,2*Ngv+1))
	end subroutine

	subroutine destroyFSens(this)
		class(FSens), intent(inout) :: this
		integer :: i

		call destroyPM1D(this%dpm)
		deallocate(this%j)
		deallocate(this%f_A)
	end subroutine

	subroutine FSensSourceTerm(this,pm)
		class(FSens), intent(inout) :: this
		type(PM1D), intent(in) :: pm
		integer :: i,k, g(pm%a(1)%order+1)
		integer :: vgl, vgr
		real(mp) :: vp, g0, h

		!F to phase space
		this%j = 0.0_mp
		do i = 1, pm%n
			do k = 1, pm%p(i)%np
				vp = pm%p(i)%vp(k,1)
				g0 = FLOOR(vp/this%dv)
				vgl = g0 + this%ngv+1
				vgr = vgl+1
				if( vgl<1 .or. vgr>2*this%ngv+1 )	cycle
				g = pm%a(i)%g(k,:)
				h = vp/this%dv - g0
				this%j(g,vgl) = this%j(g,vgl) + (1.0_mp-h)*pm%p(i)%spwt(k)/pm%m%dx/this%dv*pm%a(i)%frac(k,:)
				this%j(g,vgr) = this%j(g,vgr) + h*pm%p(i)%spwt(k)/pm%m%dx/this%dv*pm%a(i)%frac(k,:)
			end do
		end do

		!Gradient in v direction
		do i=1,pm%m%ng
			this%j(i,2:2*this%ngv) = ( this%j(i,3:2*this%ngv+1)-this%j(i,1:2*this%ngv-1) )/2.0_mp/this%dv
		end do
		this%j(:,1) = 0.0_mp
		this%j(:,pm%m%ng) = 0.0_mp

		!Multiply E_A
		do i=1,2*this%ngv+1
			this%j(:,i) = this%j(:,i)*this%dpm%m%E
		end do

		!Multiply dt
		this%j = this%j*pm%dt
	end subroutine

	subroutine InjectSource(this,f,N)
		class(FSens), intent(inout) :: this
		real(mp), intent(in), dimension(this%dpm%m%ng,2*this%ngv+1) :: f
		integer, intent(in) :: N
		integer :: Nx, newN
		real(mp), allocatable :: xp0(:), vp0(:,:), spwt0(:)
		real(mp), allocatable :: temp_x(:), temp_v(:,:), temp_w(:)
		integer :: i,i1,i2,k,Np,nk
		integer :: vgl, vgr
		real(mp) :: w, h
		real(mp), allocatable :: frac(:,:)
		Nx = INT(SQRT(N*1.0_mp))
!		newN = Nx*Nx
		newN = N

		do i=1,this%dpm%n
			Np = this%dpm%p(i)%np + newN

			allocate(xp0(newN))
			allocate(vp0(newN,3))
			allocate(spwt0(newN))
			allocate(temp_x(Np))
			allocate(temp_v(Np,3))
			allocate(temp_w(Np))

			!Spatial distribution: Uniformly-random x dimension
			call RANDOM_NUMBER(xp0)
			xp0 = xp0*this%dpm%L
			!Apply periodic BC
			do k=1,newN
				if( xp0(k)<0.0_mp ) then
					xp0(k) = xp0(k) + this%dpm%L
				elseif( xp0(k)>this%dpm%L ) then
					xp0(k) = xp0(k) - this%dpm%L
				end if
			end do

			!Velocity distribution: Gaussian-random v dimension
			vp0 = randn(newN,3)
			w = 0.4_mp*this%Lv
			vp0 = vp0*w
			!Velocity distribution: Uniformly-random x dimension
!			call RANDOM_NUMBER(vp0)
!			vp0 = (2.0_mp*vp0-1.0_mp)*this%Lv

			!Uniform grid distribution on phase space
!			xp0 = 0.0_mp
!			vp0 = 0.0_mp
!			spwt0 = 0.0_mp
!			do i2=1,Nx
!				do i1=1,Nx
!					xp0(i1+Nx*(i2-1)) = (i1-0.5_mp)*this%dpm%L/Nx
!					vp0(i1+Nx*(i2-1),:) = (i2-0.5_mp)*2.0_mp*this%Lv/Nx - 1.0_mp*this%Lv
!				end do
!			end do

			!X-direction Interpolation
			call this%dpm%a(i)%assignMatrix(this%dpm%m,xp0)
			!Adjust grid for periodic BC
			do k=1,newN
				if( this%dpm%a(i)%g(k,1)<1 ) then
					this%dpm%a(i)%g(k,1) = this%dpm%a(i)%g(k,1) + this%dpm%m%ng
				elseif( this%dpm%a(i)%g(k,1)>this%dpm%m%ng ) then
					this%dpm%a(i)%g(k,1) = this%dpm%a(i)%g(k,1) - this%dpm%m%ng
				end if
				if( this%dpm%a(i)%g(k,2)<1 ) then
					this%dpm%a(i)%g(k,2) = this%dpm%a(i)%g(k,2) + this%dpm%m%ng
				elseif( this%dpm%a(i)%g(k,2)>this%dpm%m%ng ) then
					this%dpm%a(i)%g(k,2) = this%dpm%a(i)%g(k,2) - this%dpm%m%ng
				end if
			end do

			!V-direction Interpolation and determine spwt
			allocate(frac(this%dpm%a(i)%order+1,2))
			frac = 0.0_mp
			do k=1,newN
				vgl = FLOOR(vp0(k,1)/this%dv) + this%ngv+1
				vgr = vgl+1

				if( vgl<1 .or. vgr>2*this%ngv+1 ) then
					spwt0(k) = 0.0_mp
					cycle
				end if

				h = vp0(k,1)/this%dv - FLOOR(vp0(k,1)/this%dv)
				frac(:,1) = this%dpm%a(i)%frac(k,:)*(1.0_mp-h)
				frac(:,2) = this%dpm%a(i)%frac(k,:)*h
				spwt0(k) = SUM( f(this%dpm%a(i)%g(k,:),(/vgl,vgr/))*frac )	&
								*this%dpm%L*sqrt(2.0_mp*pi)*w/EXP( -vp0(k,1)**2/2.0_mp/w/w )
!								*this%dpm%L*2.0_mp*this%Lv
			end do
			spwt0 = spwt0/newN

			!Append to sensitivity particles
			temp_x(1:this%dpm%p(i)%np) = this%dpm%p(i)%xp
			temp_v(1:this%dpm%p(i)%np,:) = this%dpm%p(i)%vp
			temp_w(1:this%dpm%p(i)%np) = this%dpm%p(i)%spwt
			temp_x(this%dpm%p(i)%np+1:Np) = xp0
			temp_v(this%dpm%p(i)%np+1:Np,:) = vp0
			temp_w(this%dpm%p(i)%np+1:Np) = spwt0
			call this%dpm%p(i)%setSpecies(Np,temp_x,temp_v,temp_w)

			deallocate(temp_x)
			deallocate(temp_v)
			deallocate(temp_w)
		end do
	end subroutine

	subroutine FSensDistribution(this)
		class(FSens), intent(inout) :: this
		integer :: i,k, g(this%dpm%a(1)%order+1)
		integer :: vgl, vgr
		real(mp) :: vp, h

		!F to phase space
		this%f_A = 0.0_mp
		do i = 1, this%dpm%n
			do k = 1, this%dpm%p(i)%np
				vp = this%dpm%p(i)%vp(k,1)
				vgl = FLOOR(vp/this%dv) + this%ngv+1
				vgr = vgl+1
				if( vgl<1 .or. vgr>2*this%ngv+1 )	cycle
				g = this%dpm%a(i)%g(k,:)
				h = vp/this%dv - FLOOR(vp/this%dv)
				this%f_A(g,vgl) = this%f_A(g,vgl) + (1.0_mp-h)*this%dpm%p(i)%spwt(k)/this%dpm%m%dx/this%dv*this%dpm%a(i)%frac(k,:)
				this%f_A(g,vgr) = this%f_A(g,vgr) + h*this%dpm%p(i)%spwt(k)/this%dpm%m%dx/this%dv*this%dpm%a(i)%frac(k,:)
			end do
		end do
	end subroutine

	subroutine Redistribute(this)
		class(FSens), intent(inout) :: this
print *, 'error'
print *, 'Np: ', this%dpm%p(1)%np
print *, 'LIMIT: ', INT(1.5_mp*this%NLimit)
		if( this%dpm%p(1)%np.ge.INT(1.5_mp*this%NLimit) ) then
			call this%FSensDistribution
			call this%dpm%p(1)%destroySpecies
			call this%InjectSource(this%f_A,INT(0.5_mp*this%NLimit))
			print *, INT(0.5_mp*this%NLimit)
		end if
	end subroutine

end module
