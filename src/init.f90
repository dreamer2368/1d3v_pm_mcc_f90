module init

	use modFSens

	implicit none

contains

	subroutine Debye_sensitivity_init_sync(dpm,pm,vT,input_str)
		type(FSens), intent(inout) :: dpm
		type(PM1D), intent(in) :: pm
		real(mp), intent(in) :: vT
		character(len=*), intent(in), optional :: input_str
		real(mp) :: xp0(pm%p(1)%np), vp0(pm%p(1)%np,3), spwt0(pm%p(1)%np)
		real(mp) :: w, rho_back(dpm%ng), xg(dpm%ng)
		integer :: N,i1, i2
		N = pm%p(1)%np
		xp0 = pm%p(1)%xp
		vp0 = pm%p(1)%vp

		if( PRESENT(input_str) ) then
			SELECT CASE(input_str)
				CASE('vT')
					spwt0 = ( vp0(:,1)**2/vT/vT - 1.0_mp )/SQRT(2.0_mp*pi)/vT/vT*EXP( -vp0(:,1)**2/2.0_mp/vT/vT )	&
								*dpm%L*2.4_mp*dpm%psM(1)%Lv/N
!								*dpm%L*SQRT(2.0_mp*pi)*vT*vT/EXP( -vp0(:,1)**2/2.0_mp/vT/vT )/N
					rho_back = 0.0_mp
				CASE('Q')
					spwt0 = 0.0_mp
					w = dpm%L/10.0_mp
					xg = (/ ((i1-0.5_mp)*dpm%m%dx, i1=1,dpm%ng) /)
					rho_back = -1.0_mp/dpm%L + 1.0_mp/SQRT(2.0_mp*pi)/w*EXP( -(xg-0.5_mp*dpm%L)**2/2.0_mp/w/w  )
				CASE('qp')
					spwt0 = 0.0_mp
					rho_back = -1.0_mp
			END SELECT
		else		!Default case: vT
			spwt0 = ( vp0(:,1)**2/vT/vT - 1.0_mp )/SQRT(2.0_mp*pi)/vT/vT*EXP( -vp0(:,1)**2/2.0_mp/vT/vT )	&
						*dpm%L*2.4_mp*dpm%psM(1)%Lv/N
			rho_back = 0.0_mp
		end if

		call dpm%p(1)%setSpecies(N,xp0,vp0,spwt0)
		call dpm%m%setMesh(rho_back)
	end subroutine

	subroutine Debye_sensitivity_init(fs,Np,vT,input_str)
		type(FSens), intent(inout) :: fs
		integer, intent(in) :: Np
		real(mp), intent(in) :: vT
		character(len=*), intent(in), optional :: input_str
		real(mp), allocatable :: xp0(:), vp0(:,:), spwt0(:)
		real(mp) :: w, rho_back(fs%ng), xg(fs%ng)
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
				xp0(i1+Nx*(i2-1)) = (i1-0.5_mp)*fs%L/Nx
				vp0(i1+Nx*(i2-1),:) = (i2-0.5_mp)*2.4_mp*fs%psM(1)%Lv/Nx - 1.2_mp*fs%psM(1)%Lv
			end do
		end do

		if( PRESENT(input_str) ) then
			SELECT CASE(input_str)
				CASE('vT')
					spwt0 = ( vp0(:,1)**2/vT/vT - 1.0_mp )/SQRT(2.0_mp*pi)/vT/vT*EXP( -vp0(:,1)**2/2.0_mp/vT/vT )	&
			!					*fs%L*w/EXP( -vp0(:,1)**2/2.0_mp/w/w )/newN
								*fs%L*2.4_mp*fs%psM(1)%Lv/newN
					rho_back = 0.0_mp
				CASE('Q')
					spwt0 = 0.0_mp
					w = fs%L/10.0_mp
					xg = (/ ((i1-0.5_mp)*fs%m%dx, i1=1,fs%ng) /)
					rho_back = -1.0_mp/fs%L + 1.0_mp/SQRT(2.0_mp*pi)/w*EXP( -(xg-0.5_mp*fs%L)**2/2.0_mp/w/w  )
				CASE('qp')
					spwt0 = 0.0_mp
					rho_back = -1.0_mp
			END SELECT
		else		!Default case: vT
			spwt0 = ( vp0(:,1)**2/vT/vT - 1.0_mp )/SQRT(2.0_mp*pi)/vT/vT*EXP( -vp0(:,1)**2/2.0_mp/vT/vT )	&
!					*fs%L*w/EXP( -vp0(:,1)**2/2.0_mp/w/w )/newN
						*fs%L*2.4_mp*fs%psM(1)%Lv/newN
			rho_back = 0.0_mp
		end if

!		call RANDOM_NUMBER(vp0)
!		w = vT*5.0_mp
!		vp0 = (2.0_mp*vp0-1.0_mp)*w
!		spwt0 = 2.0_mp*w/Np*( vp0(:,1)**2/vT/vT - 1.0_mp )/SQRT(2.0_mp*pi)/vT/vT*EXP( -vp0(:,1)**2/2.0_mp/vT/vT )

		call fs%p(1)%setSpecies(newN,xp0,vp0,spwt0)
		call fs%m%setMesh(rho_back)

		deallocate(xp0)
		deallocate(vp0)
		deallocate(spwt0)
	end subroutine

	subroutine Debye_initialize(pm,Np,Q)
		type(PM1D), intent(inout) :: pm
		integer, intent(in) :: Np
		real(mp), intent(in) :: Q
		real(mp) :: rho_back(pm%ng), xg(pm%ng)
		real(mp), allocatable :: xp0(:), vp0(:,:), spwt0(:)
		real(mp) :: L,w
		real(mp) :: vT, Lv
		integer :: i,i1,i2,j,N,Nx,newN
		L=pm%L
		N=pm%n
		w=L/10.0_mp
		xg=(/ ((i-0.5_mp)*pm%m%dx,i=1,pm%ng) /)

		!Background charge + External charge
		rho_back = -pm%p(1)%qs - Q/L + Q/sqrt(2.0_mp*pi)/w*exp( -(xg-0.5_mp*L)**2/2.0_mp/w/w )
		call pm%m%setMesh(rho_back)

        if( getOption('sensitivity_pdf/discretization','non-collocated').eq.'collocated' ) then
    		!Uniform grid distribution
    		vT = pm%A0(1)
    		Lv = 6.0_mp*vT
    		Nx = INT(SQRT(Np*1.0_mp))
    		newN = Nx*Nx
    		allocate(xp0(newN))
    		allocate(vp0(newN,3))
    		allocate(spwt0(newN))
    		xp0 = 0.0_mp
    		vp0 = 0.0_mp
    		spwt0 = 0.0_mp
    
    		do i2=1,Nx
    			do i1=1,Nx
    				xp0(i1+Nx*(i2-1)) = (i1-0.5_mp)*pm%L/Nx
    				vp0(i1+Nx*(i2-1),:) = (i2-0.5_mp)*2.4_mp*Lv/Nx - 1.2_mp*Lv
    			end do
    		end do
    
    		spwt0 = 1/SQRT(2.0_mp*pi)/vT*EXP( -vp0(:,1)**2/2.0_mp/vT/vT )	&
    					*pm%L*2.4_mp*Lv/newN

    		call pm%p(1)%setSpecies(newN,xp0,vp0,spwt0)
        else
    		!Gaussian Random distribution
    		allocate(xp0(Np))
    		allocate(vp0(Np,3))
    		allocate(spwt0(Np))
    !		call init_random_seed
    		vp0 = randn(Np,3)
    		vp0 = pm%A0(1)*vp0
    		call RANDOM_NUMBER(xp0)
    		xp0 = L*xp0
    		spwt0 = L/Np

    		call pm%p(1)%setSpecies(Np,xp0,vp0,spwt0)
        end if

		deallocate(xp0)
		deallocate(vp0)
		deallocate(spwt0)
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

	subroutine sheath_initialize(this,Ne,Ni,tau,mu)
		class(PM1D), intent(inout) :: this
		integer, intent(in) :: Ne, Ni
		real(mp), intent(in) :: tau, mu
		real(mp) :: Vth_e, Vth_i
		real(mp) :: xpe(Ne), vpe(Ne,3), spwt_e(Ne), xpi(Ni), vpi(Ni,3), spwt_i(Ni)
        integer :: i

		call RANDOM_NUMBER(xpe)
		call RANDOM_NUMBER(xpi)
		xpe = xpe*this%L
		xpi = xpi*this%L

		Vth_e = 1.0_mp
		Vth_i = sqrt(1.0_mp/mu/tau)
		print *, 'Vth_e: ',Vth_e,', Vth_i: ',Vth_i
		vpe = Vth_e*randn(Ne,3)
		vpi = Vth_i*randn(Ni,3)

		spwt_e = this%L/Ne
		spwt_i = this%L/Ni

		call this%p(1)%setSpecies(Ne,xpe,vpe,spwt_e)
		call this%p(2)%setSpecies(Ni,xpi,vpi,spwt_i)

		call this%m%setMesh((/ (0.0_mp, i=1,this%m%ng) /))
	end subroutine

	subroutine sheath_DerivativeToTau_initialize(this,Ne,Ni,tau,mu)
		class(FSens), intent(inout) :: this
		integer, intent(in) :: Ne, Ni
		real(mp), intent(in) :: tau, mu
		real(mp) :: Vth_e, Vth_i
        real(mp) :: Z1, Z2
		real(mp) :: xpe(Ne), vpe(Ne,3), spwt_e(Ne), xpi(Ni), vpi(Ni,3), spwt_i(Ni)
        integer :: i

		call RANDOM_NUMBER(xpe)
		call RANDOM_NUMBER(xpi)
		xpe = xpe*this%L
		xpi = xpi*this%L

		Vth_e = 1.0_mp
		Vth_i = sqrt(1.0_mp/mu/tau)
		print *, 'Vth_e: ',Vth_e,', Vth_i: ',Vth_i
		call RANDOM_NUMBER(vpe)
		call RANDOM_NUMBER(vpi)
		vpe = 5.5_mp*Vth_e*(2.0_mp*vpe-1.0_mp)
		vpi = 5.5_mp*Vth_i*(2.0_mp*vpi-1.0_mp)

		spwt_e = 0.0_mp

        Z1 = SUM( EXP( -mu*tau*vpi(:,1)**2/2.0_mp ) )
        Z2 = SUM( vpi(:,1)**2*EXP( -mu*tau*vpi(:,1)**2/2.0_mp ) )
		spwt_i = this%L/2.0_mp/tau*(1.0_mp/Z1 - vpi(:,1)**2/Z2)*EXP( -mu*tau*vpi(:,1)**2/2.0_mp )

		call this%p(1)%setSpecies(Ne,xpe,vpe,spwt_e)
		call this%p(2)%setSpecies(Ni,xpi,vpi,spwt_i)

		call this%m%setMesh((/ (0.0_mp, i=1,this%m%ng) /))
	end subroutine

end module
