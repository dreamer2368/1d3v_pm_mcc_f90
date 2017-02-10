module init

	use modFSens

	implicit none

contains

	subroutine Debye_sensitivity_init(fs,Np,vT)
		type(FSens), intent(inout) :: fs
		integer, intent(in) :: Np
		real(mp), intent(in) :: vT
		real(mp), allocatable :: xp0(:), vp0(:,:), spwt0(:)
		real(mp) :: w, rho_back(fs%dpm%ng)
		integer :: newN,Nx,i1, i2
		Nx = INT(SQRT(Np*1.0_mp))
		newN = Nx*Nx
		allocate(xp0(newN))
		allocate(vp0(newN,3))
		allocate(spwt0(newN))
		xp0 = 0.0_mp
		vp0 = 0.0_mp
		spwt0 = 0.0_mp

!		call RANDOM_NUMBER(xp0)
!		xp0 = fs%dpm%L*xp0

!		vp0 = randn(newN,3)
!		w = vT*1.5_mp
!		vp0 = vp0*w

		do i2=1,Nx
			do i1=1,Nx
				xp0(i1+Nx*(i2-1)) = (i1-0.5_mp)*fs%dpm%L/Nx
				vp0(i1+Nx*(i2-1),:) = (i2-0.5_mp)*2.0_mp*fs%Lv/Nx - 1.0_mp*fs%Lv
			end do
		end do

		spwt0 = ( vp0(:,1)**2/vT/vT - 1.0_mp )/vT/vT*EXP( -vp0(:,1)**2/2.0_mp/vT/vT )	&
!					*fs%dpm%L*w/EXP( -vp0(:,1)**2/2.0_mp/w/w )/newN
					*fs%dpm%L*2.0_mp*fs%Lv/newN

!		call RANDOM_NUMBER(vp0)
!		w = vT*5.0_mp
!		vp0 = (2.0_mp*vp0-1.0_mp)*w
!		spwt0 = 2.0_mp*w/Np*( vp0(:,1)**2/vT/vT - 1.0_mp )/SQRT(2.0_mp*pi)/vT/vT*EXP( -vp0(:,1)**2/2.0_mp/vT/vT )

		call fs%dpm%p(1)%setSpecies(Np,xp0,vp0,spwt0)

		rho_back = 0.0_mp
		call fs%dpm%m%setMesh(rho_back)

		deallocate(xp0)
		deallocate(vp0)
		deallocate(spwt0)
	end subroutine

	subroutine Debye_initialize(pm,Np,Q)
		type(PM1D), intent(inout) :: pm
		integer, intent(in) :: Np
		real(mp), intent(in) :: Q
		real(mp) :: xp0(Np), vp0(Np,3), spwt0(Np), rho_back(pm%ng), xg(pm%ng)
		real(mp) :: L,w
		integer :: i,j,N
		L=pm%L
		N=pm%n
		w=L/10.0_mp
		xg=(/ ((i-0.5_mp)*pm%m%dx,i=1,pm%ng) /)

!		call init_random_seed
		vp0 = randn(Np,3)
		vp0 = pm%A0(1)*vp0
		call RANDOM_NUMBER(xp0)
		xp0 = L*(xp0 - 0.5_mp) + 0.5_mp*L
		spwt0 = L/Np

		rho_back = 1.0_mp - Q/L + Q/sqrt(2.0_mp*pi)/w*exp( -(xg-0.5_mp*L)**2/2.0_mp/w/w )
!		rho_back = 1.0_mp
!		rho_back = rho_back - Q/L
!		pm%p(1)%qs = pm%p(1)%qs - Q/L
!		rho_back(pm%ng/2) = rho_back(pm%ng/2) + Q/pm%m%dx
!		rho_back(pm%ng) = 0.0_mp

		call pm%p(1)%setSpecies(Np,xp0,vp0,spwt0)
		call pm%m%setMesh(rho_back)
	end subroutine

	subroutine Landau_initialize(pm,Np,vT)
		type(PM1D), intent(inout) :: pm
		integer, intent(in) :: Np
		real(mp), intent(in) :: vT
		real(mp) :: xp0(Np), vp0(Np,3), spwt0(Np), rho_back
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
		spwt0 = 1.0_mp
		do i=1,N
			call pm%p(i)%buildSpecies(qs,ms)
			call pm%p(i)%setSpecies(Np,xp0,vp0,spwt0)
		end do
		call pm%m%setMesh(rho_back*(/ ( 1.0_mp, i=1,pm%m%ng) /))
	end subroutine

	subroutine twostream_initialize(this,Np,v0,vT,mode)		!generate initial distribution
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: Np, mode
		real(mp), intent(in) :: v0, vT
		real(mp) :: xp0(Np), vp0(Np,3), spwt0(Np), rho_back
		real(mp) :: L,qs,ms
		integer :: i,j,N,pm(Np)
		L = this%L
		N = this%n

		qs = -this%wp*this%wp/(this%n*Np/L)
		ms = -qs
		rho_back = -qs*this%n*Np/L

		!velocity distribution initialize
		print *, 'vT: ',vT
		vp0 = vT*randn(Np,3)
		pm = (/ ( i, i=1,Np ) /)
		pm = 1 - 2*MOD(pm,2)
		vp0(:,1) = vp0(:,1) + pm*v0

		spwt0 = 1.0_mp

		do i=1,this%n
			!spatial distribution initialize
			xp0 = (/ ( j*L/Np, j=0,Np-1 ) /) + 0.5_mp*(i-1)*L/Np
!			call RANDOM_NUMBER(xp0)
!			xp0 = xp0*L
			xp0 = xp0 + this%A0(1)*L/Np*SIN( 2.0_mp*pi*xp0/L*mode )

			call this%p(i)%buildSpecies(qs,ms)
			call this%p(i)%setSpecies(Np,xp0,vp0,spwt0)
		end do
		call this%m%setMesh(rho_back*(/ ( 1, i=1,this%m%ng) /))
	end subroutine

	subroutine sheath_initialize(this,Ne,Ni,Te,Ti,Kb,n0)
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: Ne, Ni
		real(mp), intent(in) :: Te, Ti, Kb, n0
		real(mp) :: Vth_e, Vth_i
		real(mp) :: xpe(Ne), vpe(Ne,3), spwt_e(Ne), xpi(Ni), vpi(Ni,3), spwt_i(Ni)
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

		spwt_e = n0*this%L/Ne
		spwt_i = n0*this%L/Ni

		call this%p(1)%setSpecies(Ne,xpe,vpe,spwt_e)
		call this%p(2)%setSpecies(Ni,xpi,vpi,spwt_i)

		call this%m%setMesh((/ (0.0_mp, i=1,this%m%ng) /))
	end subroutine

end module
