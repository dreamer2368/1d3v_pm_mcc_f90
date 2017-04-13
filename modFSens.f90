module modFSens

	use modBC
	use modRecord

	implicit none

	type, extends(PM1D) :: FSens
		real(mp) :: Lv, dv
		integer :: ngv
		real(mp), allocatable :: j(:,:), f_A(:,:)

		integer :: NInject										!Number of injecting particle per direction
		real(mp), allocatable :: xp_inject(:), vp_inject(:,:), frac_inject(:,:,:)
		integer, allocatable :: g_inject(:,:), gv_inject(:,:)

		integer :: NLimit											!Number limit of particle for redistribution
		real(mp), allocatable :: xp_remesh(:), vp_remesh(:,:), frac_remesh(:,:,:)
		integer, allocatable :: g_remesh(:,:), gv_remesh(:,:)
	contains
		procedure, pass(dpm) :: buildFSens
		procedure, pass(this) :: destroyFSens
		procedure, pass(this) :: FSensDistribution
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
		allocate(dpm%j(pm%ng,2*Ngv+1))
		allocate(dpm%f_A(pm%ng,2*Ngv+1))
		dpm%j = 0.0_mp
		dpm%f_A = 0.0_mp

		call createDistribution(dpm,dpm%NInject,dpm%xp_inject,dpm%vp_inject,	&
										dpm%g_inject,dpm%gv_inject,dpm%frac_inject)
		call createDistribution(dpm,dpm%NLimit,dpm%xp_remesh,dpm%vp_remesh,	&
										dpm%g_remesh,dpm%gv_remesh,dpm%frac_remesh)
	end subroutine

	subroutine destroyFSens(this)
		class(FSens), intent(inout) :: this
		integer :: i

		call destroyPM1D(this)
		deallocate(this%j)
		deallocate(this%f_A)
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
	end subroutine

	subroutine FSensSourceTerm(dpm,qs,ms,f,E,nk,fsr)							!Take v-derivative, Multiply dE
		class(FSens), intent(inout) :: dpm
		real(mp), intent(in) :: qs,ms,f(dpm%m%ng,2*dpm%ngv+1), E(dpm%m%ng)
		integer, intent(in), optional :: nk
		type(recordData), intent(in), optional :: fsr
		integer :: ng,ngv,kr,i
		character(len=100) :: kstr
		real(mp) :: dx,dv
		real(mp) :: E_temp(dpm%m%ng)
		ng = dpm%m%ng
		ngv = dpm%ngv
		dx = dpm%m%dx
		dv = dpm%dv

		!Gradient in v direction
		do i=1,ng
			dpm%j(i,2:2*ngv) = ( f(i,3:2*ngv+1)-f(i,1:2*ngv-1) )/2.0_mp/dv
		end do
		dpm%j(:,1) = 0.0_mp
		dpm%j(:,2*ngv+1) = 0.0_mp
!		if( present(nk) ) then
!			if( (fsr%mod.eq.1) .or. (mod(nk,fsr%mod).eq.0) ) then
!				kr = merge(nk,nk/fsr%mod,fsr%mod.eq.1)
!				write(kstr,*) kr
!				open(unit=305,file='data/'//fsr%dir//'/Dvf_'//trim(adjustl(kstr))//'.bin',	&
!						status='replace',form='unformatted',access='stream')
!				write(305) this%j
!				close(305)
!			end if
!		end if

		!Multiply E_A, E
		!v_T,Q
		E_temp = qs/ms*dpm%m%E
		!qp
!		E_temp = qs/ms*dpm%m%E + 1.0_mp/ms*E
		do i=1,2*ngv+1
			dpm%j(:,i) = dpm%j(:,i)*E_temp
		end do

		!Multiply dt
		dpm%j = -dpm%j*dpm%dt
	end subroutine

	subroutine InjectSource(this,p,f)
		class(FSens), intent(inout) :: this
		type(species), intent(inout) :: p
		real(mp), intent(in), dimension(this%m%ng,2*this%ngv+1) :: f
		real(mp) :: spwt0(this%NInject)

		call SourceAssign(f,this%L,this%Lv,this%g_inject,this%gv_inject,this%frac_inject,spwt0)
		call p%appendSpecies(this%NInject,this%xp_inject,this%vp_inject,spwt0)
	end subroutine

	subroutine SourceAssign(f,Lx,Lv,g,gv,frac,spwt0)
		real(mp), intent(in) :: f(:,:)
		real(mp), intent(in) :: Lx,Lv
		integer, intent(in) :: g(:,:), gv(:,:)
		real(mp), intent(in) :: frac(:,:,:)
		real(mp), intent(out) :: spwt0(size(g,2))
		integer :: k,N
		N = size(g,2)

		do k=1,N
			spwt0(k) = SUM( f(g(:,k),gv(:,k))*frac(:,:,k) )*Lx*2.0_mp*Lv/N
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

	subroutine FSensDistribution(this,p,a)
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
			this%f_A(g,vgl) = this%f_A(g,vgl) + (1.0_mp-h)*p%spwt(k)/this%m%dx/this%dv*a%frac(:,k)
			this%f_A(g,vgr) = this%f_A(g,vgr) + h*p%spwt(k)/this%m%dx/this%dv*a%frac(:,k)
		end do
		this%f_A(:,1) = this%f_A(:,1)*2.0_mp
		this%f_A(:,2*this%ngv+1) = this%f_A(:,2*this%ngv+1)*2.0_mp
	end subroutine

	subroutine Redistribute(this,p,a)
		class(FSens), intent(inout) :: this
		type(species), intent(inout) :: p
		type(pmAssign), intent(in) :: a
		real(mp) :: spwt0(this%NLimit)
		integer :: i
		if( p%np.ge.3*this%NLimit ) then
			call this%FSensDistribution(p,a)
			call SourceAssign(this%f_A,this%L,this%Lv,this%g_remesh,this%gv_remesh,this%frac_remesh,spwt0)
			call p%setSpecies(this%NLimit,this%xp_remesh,this%vp_remesh,spwt0)
		end if
	end subroutine

	subroutine numberDensity(this,p,a,N_A)
		class(FSens), intent(inout) :: this
		type(species), intent(in) :: p
		type(pmAssign), intent(in) :: a
		real(mp), dimension(this%m%ng, 2*this%ngv+1), intent(inout) :: N_A
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
!				this%gv(:,k) = this%ngv
!				this%frac(:,:,k) = 0.0_mp
				cycle
			end if
			g = g_x(:,k)
			g_v = (/ vgl, vgr /)
			h = vp/this%dv - FLOOR(vp/this%dv)
			frac(:,1) = (1.0_mp-h)*frac_x(:,k)
			frac(:,2) = h*frac_x(:,k)
			n_temp(g,g_v) = n_temp(g,g_v) + frac/this%m%dx/this%dv

!			this%gv(:,k) = g_v-1
!			this%frac(:,:,k) = frac
!			if( vgl.eq.1 ) then
!				this%frac(:,1,k) = 0.0_mp
!				this%gv(1,k) = this%ngv+2
!			end if
!			if( vgr.eq.2*this%ngv+3 ) then
!				this%frac(:,2,k) = 0.0_mp
!				this%gv(2,k) = this%ngv+2
!			end if
		end do
		N_A = N_A + n_temp(:,2:2*this%ngv+2)
	end subroutine

	subroutine updateWeight(this,p,a,n_A,j)
		class(FSens), intent(inout) :: this
		type(species), intent(inout) :: p
		type(pmAssign), intent(in) :: a
		real(mp), dimension(this%m%ng,2*this%ngv+1), intent(in) :: n_A, j
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
			p%spwt(k) = p%spwt(k) + SUM( j(g,g_v)*frac/n_A(g,g_v) )
!			g=g_x(:,k)
!			p%spwt(k) = p%spwt(k) + SUM( j(g,this%gv(:,k))*this%frac(:,:,k)/n_A(g,this%gv(:,k)) )
		end do
	end subroutine

end module
