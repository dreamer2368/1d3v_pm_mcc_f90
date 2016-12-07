module init

	use modPM1D
	use random
	implicit none

contains

	subroutine Debye_initialize(pm,Np,Q)
		type(PM1D), intent(inout) :: pm
		integer, intent(in) :: Np
		real(mp), intent(in) :: Q
		real(mp) :: xp0(Np), vp0(Np,3),rho_back(pm%ng)
		real(mp) :: L
		integer :: i,j,N
		L=pm%L
		N=pm%n

!		call init_random_seed
		vp0 = randn(Np,3)
		call RANDOM_NUMBER(xp0)
		xp0 = xp0*L

		rho_back = 1.0_mp
		rho_back(pm%ng/2) = rho_back(pm%ng/2) + Q/pm%m%dx

		call setSpecies(pm%p(1),Np,xp0,vp0)
		call setMesh(pm%m,rho_back)
	end subroutine

	subroutine Landau_initialize(pm,Np,vT)
		type(PM1D), intent(inout) :: pm
		integer, intent(in) :: Np
		real(mp), intent(in) :: vT
		real(mp) :: xp0(Np), vp0(Np,3),rho_back
		real(mp) :: L,qs,ms
		integer :: i,j,N
		L=pm%L
		N=pm%n

		qs = -pm%wp*pm%wp/(pm%n*Np/L)
		ms = -qs
		rho_back = -qs*pm%n*Np/L

!		call init_random_seed
		vp0 = vT*randn(Np,3)
		call RANDOM_NUMBER(xp0)
		xp0 = xp0*L
		xp0 = xp0 + pm%A0(1)*SIN( 2.0_mp*pi*xp0/L )
		do i=1,N
			call buildSpecies(pm%p(i),qs,ms,1.0_mp)
			call setSpecies(pm%p(i),Np,xp0,vp0)
		end do
		call setMesh(pm%m,rho_back*(/ ( 1.0_mp, i=1,pm%m%ng) /))
	end subroutine

	subroutine twostream_initialize(this,Np,v0,vT,mode)		!generate initial distribution
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: Np, mode
		real(mp), intent(in) :: v0, vT
		real(mp) :: xp0(Np), vp0(Np,3),rho_back
		real(mp) :: L,qs,ms
		integer :: i,j,N,pm(Np)
		L = this%L
		N = this%n

		qs = -this%wp*this%wp/(this%n*Np/L)
		ms = -qs
		rho_back = -qs*this%n*Np/L

		!velocity distribution initialize
		vp0 = vT*randn(Np,3)
		pm = (/ ( i, i=1,Np ) /)
		pm = 1 - 2*MOD(pm,2)
		vp0(:,1) = vp0(:,1) + pm*v0

		do i=1,this%n
			!spatial distribution initialize
!			xp0 = (/ ( j*L/Np, j=0,Np-1 ) /) + 0.5_mp*(i-1)*L/Np
			call RANDOM_NUMBER(xp0)
			xp0 = xp0*L
			xp0 = xp0 + this%A0(1)*L/Np*SIN( 2.0_mp*pi*xp0/L*mode )

			call buildSpecies(this%p(i),qs,ms,1.0_mp)
			call setSpecies(this%p(i),Np,xp0,vp0)
		end do
		call setMesh(this%m,rho_back*(/ ( 1, i=1,this%m%ng) /))
	end subroutine

	subroutine sheath_initialize(this,Ne,Ni,Te,Ti,Kb)
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: Ne, Ni
		real(mp), intent(in) :: Te, Ti, Kb
		real(mp) :: Vth_e, Vth_i
		real(mp) :: xpe(Ne), vpe(Ne,3), xpi(Ni), vpi(Ni,3)
		integer :: i,nseed,clock
		integer, allocatable :: seed(:)

		call RANDOM_SEED(size=nseed)
		allocate(seed(nseed))
		call SYSTEM_CLOCK(COUNT=clock)
		seed = clock + 127*(/ ( i, i=1,nseed ) /)
		call RANDOM_SEED(put=seed)
		call RANDOM_NUMBER(xpe)
		call RANDOM_NUMBER(xpi)
		xpe = xpe*this%L
		xpi = xpi*this%L

		Vth_e = sqrt(2.0_mp*Kb*Te/this%p(1)%ms)
		Vth_i = sqrt(2.0_mp*Kb*Ti/this%p(2)%ms)
		print *, 'Vth_e: ',Vth_e,', Vth_i: ',Vth_i
		vpe = Vth_e/sqrt(2.0_mp)*randn(Ne,3)
		vpi = Vth_i/sqrt(2.0_mp)*randn(Ni,3)

		call setSpecies(this%p(1),Ne,xpe,vpe)
		call setSpecies(this%p(2),Ni,xpi,vpi)

		call setMesh(this%m, (/ (0.0_mp, i=1,this%m%ng) /))
	end subroutine

end module
