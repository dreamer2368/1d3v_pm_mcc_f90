module modFSens

	use modRecord

	implicit none

!	type FSens
!		type(PM1D) :: dpm
!		real(mp) :: Lv, dv
!		integer :: ngv
!		real(mp), allocatable :: j(:,:), f_A(:,:), frac(:,:,:)
!		integer, allocatable :: gv(:,:)
!		integer :: NInject										!Number of injecting particle per direction
!		integer :: NLimit											!Number limit of particle for redistribution
!	contains
!		procedure, pass(this) :: buildFSens
!		procedure, pass(this) :: destroyFSens
!		procedure, pass(this) :: FSensSourceTerm
!		procedure, pass(this) :: InjectSource
!		procedure, pass(this) :: FSensDistribution
!		procedure, pass(this) :: Redistribute
!		procedure, pass(this) :: updateWeight
!	end type

contains

	subroutine buildSensitivity(dpm,pm)
		class(PM1D), intent(out) :: dpm
		type(PM1D), intent(in) :: pm
		integer :: i

		call dpm%buildPM1D(pm%nt*pm%dt,pm%ni*pm%dt,pm%ng,pm%n,	&
							pm%pBCindex,pm%mBCindex,pm%a(1)%order,	&
							pm%dt,pm%L,pm%A0,pm%eps0,pm%ngv,pm%Lv)
		do i=1,pm%n
			call dpm%p(i)%buildSpecies(pm%p(i)%qs,pm%p(i)%ms)
		end do
		call dpm%m%setMesh(pm%m%rho_back)
	end subroutine

	subroutine FSensSourceTerm(dpm,pm,nk,fsr)
		class(PM1D), intent(inout) :: dpm
		type(PM1D), intent(in) :: pm
		integer, intent(in), optional :: nk
		type(recordData), intent(in), optional :: fsr
		integer :: kr
		character(len=100) :: kstr
		integer :: i,k
		real(mp) :: E_temp(pm%m%ng)

		!Dv_F to phase space
		dpm%m%f = 0.0_mp
!		do i=1,pm%n
			call dpm%a(1)%PhaseSpaceDensity(pm%p(1),dpm%m)
			dpm%m%j = D_V(dpm%m%f,dpm%m%dv)
!		end do

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
		E_temp = pm%p(1)%qs/pm%p(1)%ms*dpm%m%E
		!qp
!		E_temp = pm%p(1)%qs/pm%p(1)%ms*dpm%m%E + 1.0_mp/pm%p(1)%ms*pm%m%E
		do i=1,2*dpm%m%ngv+1
			dpm%m%j(:,i) = dpm%m%j(:,i)*E_temp
		end do

		!Multiply dt
		dpm%m%j = -dpm%m%j*pm%dt
	end subroutine

	subroutine InjectSource(dpm,f,N)
		class(PM1D), intent(inout) :: dpm
		real(mp), intent(in), dimension(dpm%m%ng,2*dpm%m%ngv+1) :: f
		integer, intent(in) :: N
		real(mp), allocatable :: xp0(:), vp0(:,:), spwt0(:)

		call createDistribution(dpm,f,N,xp0,vp0,spwt0)
		call dpm%p(1)%appendSpecies(size(xp0),xp0,vp0,spwt0)
	end subroutine

	subroutine createDistribution(dpm,f,N,xp0,vp0,spwt0)
		class(PM1D), intent(inout) :: dpm
		real(mp), intent(in), dimension(dpm%m%ng,2*dpm%m%ngv+1) :: f
		integer, intent(in) :: N
		real(mp), intent(out), allocatable :: xp0(:), vp0(:,:), spwt0(:)
		integer :: Nx, newN, i1,i2,k
		real(mp) :: dx,dv
		integer :: g(dpm%a(1)%order+1),gv(2)
		real(mp) :: frac(dpm%a(1)%order+1,2)
		dx = dpm%m%dx
		dv = dpm%m%dv
		Nx = INT(SQRT(N*1.0_mp))
		newN = Nx*Nx
!		newN = N

		allocate(xp0(newN))
		allocate(vp0(newN,3))
		allocate(spwt0(newN))

			!Spatial distribution: Uniformly-random x dimension
!			call RANDOM_NUMBER(xp0)
!			xp0 = xp0*this%dpm%L
!			!Apply periodic BC
!			do k=1,newN
!				if( xp0(k)<0.0_mp ) then
!					xp0(k) = xp0(k) + this%dpm%L
!				elseif( xp0(k)>this%dpm%L ) then
!					xp0(k) = xp0(k) - this%dpm%L
!				end if
!			end do

			!Velocity distribution: Gaussian-random v dimension
!			vp0 = randn(newN,3)
!			w = 0.4_mp*this%Lv
!			vp0 = vp0*w
			!Velocity distribution: Uniformly-random x999999 dimension
!			call RANDOM_NUMBER(vp0)
!			vp0 = (2.0_mp*vp0-1.0_mp)*this%Lv

		!Uniform grid distribution on phase space
		xp0 = 0.0_mp
		vp0 = 0.0_mp
		spwt0 = 0.0_mp
		do i2=1,Nx
			do i1=1,Nx
				xp0(i1+Nx*(i2-1)) = (i1-0.5_mp)*dpm%L/Nx
				vp0(i1+Nx*(i2-1),:) = (i2-0.5_mp)*2.0_mp*dpm%Lv/Nx - 1.0_mp*dpm%Lv
			end do
		end do

		!Phase-Space Interpolation and determine spwt
		do k=1,newN
			call dpm%a(1)%AssignPhaseSpace(xp0(k),vp0(k,1),dpm%m,g,gv,frac)
			spwt0(k) = SUM( f( g, gv )*frac )	&
!							*this%dpm%L*sqrt(2.0_mp*pi)*w/EXP( -vp0(k,1)**2/2.0_mp/w/w )
							*dpm%L*2.0_mp*dpm%Lv
		end do
		spwt0 = spwt0/newN
	end subroutine

	subroutine Redistribute(pm,NLimit)
		class(PM1D), intent(inout) :: pm
		integer, intent(in) :: NLimit
		real(mp), allocatable :: xp0(:), vp0(:,:), spwt0(:)
		integer :: i
		if( pm%p(1)%np.ge.INT(1.5_mp*NLimit) ) then
			pm%m%f = 0.0_mp
			call pm%a(1)%PhaseSpaceDensity(pm%p(1),pm%m)
			call createDistribution(pm,pm%m%f,INT(0.5_mp*NLimit),xp0,vp0,spwt0)
			call pm%p(1)%setSpecies(size(xp0),xp0,vp0,spwt0)
		end if
	end subroutine

	subroutine updateWeight(p,a,m,j)
		class(species), intent(inout) :: p
		type(pmAssign), intent(in) :: a
		type(mesh), intent(inout) :: m
		real(mp), intent(in), dimension(m%ng,2*m%ngv+1) :: j
		integer :: k, g(a%order+1),gv(2)
		real(mp) :: frac(a%order+1,2)

		!Get N_A
		m%f = 0.0_mp
		call a%NumberDensity(p,m)

		!Update Weight
		do k = 1, p%np
			call a%AssignPhaseSpace(p%xp(k),p%vp(k,1),m,g,gv,frac)
			p%spwt(k) = p%spwt(k) + SUM( j(g,gv)*frac/m%f(g,gv) )
		end do
	end subroutine

end module
