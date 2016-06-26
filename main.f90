program main

	use init
	use timeStep

	implicit none

	! print to screen
	print *, 'calling program main'

!	call cross_section
!	call Procassini
!	call test_refluxing_boundary
!	call test_anewvel_Ar
!	call test_mcc_electron
	call test_mcc_Argon
!	call test_ext_voltage_Poisson

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

	subroutine test_ext_voltage_Poisson
		type(mesh) :: m
		integer, parameter :: N=128
		real(mp) :: xg(N), rho_back(N), rho(N), sol(N)
		integer :: i

		xg = (/ (i-1, i=1,N) /)*1.0_mp/(N-1)
		rho = exp(3.0_mp*xg)
		rho_back = 0.0_mp
		rho_back(N) = 1.0_mp

		sol = ( 1.0_mp-exp(3.0_mp*xg) )/9.0_mp + ( 1.0_mp - (1.0_mp-exp(3.0_mp))/9.0_mp )*xg

		call buildMesh(m,1.0_mp,N,1)
		call setMesh(m,rho_back)
		m%rho = rho
		call solveMesh(m,1.0_mp)

		call system('mkdir -p data/test_DD_poisson')
		open(unit=301,file='data/test_DD_poisson/xg.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/test_DD_poisson/phi.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/test_DD_poisson/phi_sol.bin',status='replace',form='unformatted',access='stream')
		write(301) xg
		write(302) m%phi
		write(303) sol
		close(301)
		close(302)
		close(303)

		print *, 'error: ', maxval( abs(sol - m%phi) )
	end subroutine

	subroutine test_mcc_Argon
		type(PM1D) :: pm
		integer, parameter :: np = 100000, Nsample=10000
		real(mp) :: dt = log(100.0_mp/99.0_mp)
		real(mp) :: gden = 1.0_mp/max_sigmav_Ar, TN = 0.026_mp		!neutral temperature TN: scale in eV
		real(mp) :: energy = 6.0_mp, vel										!argon energy
		real(mp) :: xp0(np), vp0(np,3)
		integer :: i, N_coll(Nsample,3)

      call init_random_seed

		call null_collision(gden,dt)
		call buildPM1D(pm,30.0_mp, 15.0_mp,16,2,0,0,1,dt,L=1.0_mp)
		call set_Ar_discharge(pm,(/1.0_mp, 1.0_mp/),(/TN,gden,0.0_mp,0.0_mp/))

		vel = sqrt(2.0_mp/pm%p(2)%ms*q_e*energy)
		xp0 = (/ (i-0.5_mp,i=1,np) /)*(1.0_mp/np)
		vp0(:,1) = 0.0_mp
		vp0(:,2) = vel
		vp0(:,3) = 0.0_mp
		call setSpecies(pm%p(2),np,xp0,vp0)
		vp0 = 0.0_mp
		call setSpecies(pm%p(1),np,xp0,vp0)

		call system('mkdir -p data/test_mcc_argon')

		call system('mkdir -p data/test_mcc_argon/before')
		open(unit=301,file='data/test_mcc_argon/before/np.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/test_mcc_argon/before/xp_e.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/test_mcc_argon/before/vp_e.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/test_mcc_argon/before/xp_Ar.bin',status='replace',form='unformatted',access='stream')
		open(unit=305,file='data/test_mcc_argon/before/vp_Ar.bin',status='replace',form='unformatted',access='stream')
		write(301) np
		write(302) pm%p(1)%xp
		write(303) pm%p(1)%vp
		write(304) pm%p(2)%xp
		write(305) pm%p(2)%vp
		close(301)
		close(302)
		close(303)
		close(304)
		close(305)

		call mcc_argon(pm)

		call system('mkdir -p data/test_mcc_argon/after')
		open(unit=301,file='data/test_mcc_argon/after/np_e.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/test_mcc_argon/after/xp_e.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/test_mcc_argon/after/vp_e.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/test_mcc_argon/after/xp_Ar.bin',status='replace',form='unformatted',access='stream')
		open(unit=305,file='data/test_mcc_argon/after/vp_Ar.bin',status='replace',form='unformatted',access='stream')
		open(unit=306,file='data/test_mcc_argon/after/np_Ar.bin',status='replace',form='unformatted',access='stream')
		write(301) pm%p(1)%np
		write(302) pm%p(1)%xp
		write(303) pm%p(1)%vp
		write(304) pm%p(2)%xp
		write(305) pm%p(2)%vp
		write(306) pm%p(2)%np
		close(301)
		close(302)
		close(303)
		close(304)
		close(305)
		close(306)

		call destroyPM1D(pm)

		call null_collision(gden,dt)
		vel = sqrt(2.0_mp/m_Ar*q_e*energy)
		xp0 = (/ (i-0.5_mp,i=1,np) /)*(1.0_mp/np)
		vp0(:,1) = 0.0_mp
		vp0(:,2) = vel
		vp0(:,3) = 0.0_mp

		open(unit=301,file='data/test_mcc_Argon/prob.bin',status='replace',form='unformatted',access='stream')
		write(301)	col_prob_Ar,	&
                 col_prob_Ar*( asigma4(energy)*vel/max_sigmav_Ar ),   &
                 col_prob_Ar*( asigma5(energy)*vel/max_sigmav_Ar )
		close(301)

		do i=1,Nsample
         print *, ' '
         print *, '================',i,'-th Sample=================='
         print *, ' '
   		call buildPM1D(pm,30.0_mp, 15.0_mp,16,2,0,0,1,dt,L=1.0_mp)
   		call set_Ar_discharge(pm,(/1.0_mp, 1.0_mp/),(/TN,gden,0.0_mp,0.0_mp/))
			call setSpecies(pm%p(1),np,xp0,vp0)
			call setSpecies(pm%p(2),np,xp0,vp0)
         pm%p(1)%vp = 0.0_mp

			call mcc_argon(pm,N_coll(i,:))

			call destroyPM1D(pm)
		end do

		open(unit=301,file='data/test_mcc_Argon/coll_sample.bin',status='replace',form='unformatted',access='stream')
      write(301) N_coll
      close(301)
	end subroutine

	subroutine test_mcc_electron
		type(PM1D) :: pm
		integer, parameter :: np = 100000, Nsample = 10000
		real(mp) :: dt = log(100.0_mp/99.0_mp)
		real(mp) :: gden = 1.0_mp/max_sigmav_e, TN = 0.026_mp		!neutral temperature TN: scale in eV
		real(mp) :: energy = 20.0_mp, vel										!electron energy
		real(mp) :: xp0(np), vp0(np,3)
		integer :: i, N_coll(Nsample,4)
		character(len=100) :: istr

      call init_random_seed

		call null_collision(gden,dt)
		call buildPM1D(pm,30.0_mp, 15.0_mp,16,2,0,0,1,dt,L=1.0_mp)
		call set_Ar_discharge(pm,(/1.0_mp, 1.0_mp/),(/TN,gden,0.0_mp,0.0_mp/))

		vel = sqrt(2.0_mp/pm%p(1)%ms*q_e*energy)
		xp0 = (/ (i-0.5_mp,i=1,np) /)*(1.0_mp/np)
		vp0(:,1) = 0.0_mp
		vp0(:,2) = vel
		vp0(:,3) = 0.0_mp
		call setSpecies(pm%p(1),np,xp0,vp0)
		vp0 = 0.0_mp
		call setSpecies(pm%p(2),np,xp0,vp0)

		call system('mkdir -p data/test_mcc_electron/before')
		open(unit=301,file='data/test_mcc_electron/before/np.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/test_mcc_electron/before/xp_e.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/test_mcc_electron/before/vp_e.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/test_mcc_electron/before/xp_Ar.bin',status='replace',form='unformatted',access='stream')
		open(unit=305,file='data/test_mcc_electron/before/vp_Ar.bin',status='replace',form='unformatted',access='stream')
		write(301) np
		write(302) pm%p(1)%xp
		write(303) pm%p(1)%vp
		write(304) pm%p(2)%xp
		write(305) pm%p(2)%vp
		close(301)
		close(302)
		close(303)
		close(304)
		close(305)

		call mcc_electron(pm)

		call system('mkdir -p data/test_mcc_electron/after')
		open(unit=301,file='data/test_mcc_electron/after/np_e.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/test_mcc_electron/after/xp_e.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/test_mcc_electron/after/vp_e.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/test_mcc_electron/after/xp_Ar.bin',status='replace',form='unformatted',access='stream')
		open(unit=305,file='data/test_mcc_electron/after/vp_Ar.bin',status='replace',form='unformatted',access='stream')
		open(unit=306,file='data/test_mcc_electron/after/np_Ar.bin',status='replace',form='unformatted',access='stream')
		write(301) pm%p(1)%np
		write(302) pm%p(1)%xp
		write(303) pm%p(1)%vp
		write(304) pm%p(2)%xp
		write(305) pm%p(2)%vp
		write(306) pm%p(2)%np
		close(301)
		close(302)
		close(303)
		close(304)
		close(305)
		close(306)
		call destroyPM1D(pm)

		call null_collision(gden,dt)
		vel = sqrt(2.0_mp/m_e*q_e*energy)
		xp0 = (/ (i-0.5_mp,i=1,np) /)*(1.0_mp/np)
		vp0(:,1) = 0.0_mp
		vp0(:,2) = vel
		vp0(:,3) = 0.0_mp

		open(unit=301,file='data/test_mcc_electron/prob.bin',status='replace',form='unformatted',access='stream')
		write(301)	col_prob_e,	&
                 col_prob_e*( asigma1(energy)*vel/max_sigmav_e ),   &
                 col_prob_e*( asigma2(energy)*vel/max_sigmav_e ),   &
                 col_prob_e*( asigma3(energy)*vel/max_sigmav_e )
		close(301)

		do i=1,Nsample
         print *, ' '
         print *, '================',i,'-th Sample=================='
         print *, ' '
   		call buildPM1D(pm,30.0_mp, 15.0_mp,16,2,0,0,1,dt,L=1.0_mp)
   		call set_Ar_discharge(pm,(/1.0_mp, 1.0_mp/),(/TN,gden,0.0_mp,0.0_mp/))
			call setSpecies(pm%p(1),np,xp0,vp0)
			call setSpecies(pm%p(2),np,xp0,vp0)
         pm%p(2)%vp = 0.0_mp

			call mcc_electron(pm,N_coll(i,:))

			call destroyPM1D(pm)
		end do

		open(unit=301,file='data/test_mcc_electron/coll_sample.bin',status='replace',form='unformatted',access='stream')
      write(301) N_coll
      close(301)
	end subroutine

	subroutine test_anewvel_Ar
		integer, parameter :: N = 10000
		real(mp) :: input(3) = (/ 0.0_mp, 1.0_mp, 0.0_mp /)
		real(mp), dimension(N,3) :: output
		integer :: i

		do i=1,N
			output(i,:) = input
			call anewvel_Ar(output(i,:))
		end do

		call system('mkdir -p data/scattering')
		open(unit=301,file='data/scattering/output_Ar.bin',status='replace',form='unformatted',access='stream')
		write(301) output
		close(301)
	end subroutine

	subroutine test_anewvel_e
		integer, parameter :: N = 100000
		real(mp) :: input(3) = (/ 0.0_mp, 1.0_mp, 0.0_mp /)
		real(mp), dimension(N,3) :: output
		real(mp) :: energy = 100.0_mp			!eV
		integer :: i

		do i=1,N
			output(i,:) = input
			call anewvel_e(energy, 1.0_mp, 1.0_mp, output(i,:),.false.)
		end do

		call system('mkdir -p data/scattering')
		open(unit=301,file='data/scattering/output.bin',status='replace',form='unformatted',access='stream')
		write(301) output
		close(301)
	end subroutine

	subroutine cross_section
		integer, parameter :: N=10000
		real(mp), dimension(N) :: energy, sig1, sig2, sig3, sig4, sig5
		integer :: i

		energy = exp( log(10.0_mp)*( (/ (i,i=1,N) /)/(0.2_mp*N) - 2.0_mp ) )
		do i=1,N
			sig1(i) = asigma1(energy(i))
			sig2(i) = asigma2(energy(i))
			sig3(i) = asigma3(energy(i))
			sig4(i) = asigma4(energy(i))
			sig5(i) = asigma5(energy(i))
		end do

		call system('mkdir -p data/cross_section')
		open(unit=301,file='data/cross_section/sig1.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/cross_section/sig2.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/cross_section/sig3.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/cross_section/sig4.bin',status='replace',form='unformatted',access='stream')
		open(unit=305,file='data/cross_section/sig5.bin',status='replace',form='unformatted',access='stream')
		open(unit=306,file='data/cross_section/energy.bin',status='replace',form='unformatted',access='stream')
		write(301) sig1
		write(302) sig2
		write(303) sig3
		write(304) sig4
		write(305) sig5
		write(306) energy
		close(301)
		close(302)
		close(303)
		close(304)
		close(305)
		close(306)
	end subroutine

	subroutine Procassini
		type(PM1D) :: sheath
		type(recordData) :: r
		real(mp), parameter :: Kb = 1.38065E-23, EV_TO_K = 11604.52_mp, eps = 8.85418782E-12
		real(mp), parameter :: Te = 50.0_mp*EV_TO_K, tau = 100.0_mp
		real(mp), parameter :: me = 9.10938215E-31, qe = 1.602176565E-19, mu = 1836
		real(mp), parameter :: n0 = 2.00000000E14
		integer, parameter :: Ne = 10000, Ni = 10000
		real(mp) :: mi, Ti, wp0, lambda0, dt, dx, L
		real(mp) :: ve0, vi0, Time_f
		real(mp) :: A(4)
		integer :: i

		mi = mu*me
		Ti = Te/tau
		wp0 = sqrt(n0*qe*qe/me/eps)
		lambda0 = sqrt(eps*Kb*Te/n0/qe/qe)
		L = 2.0_mp*lambda0

		print *, 'L = ',L,', lambda0 = ',lambda0,' e = lambda/L = ',lambda0/L

		dt = 0.1_mp/wp0
		dx = 0.2_mp*lambda0
		!		dt = 0.5_mp*dx/(lambda0*wp0)

		ve0 = sqrt(Kb*Te/me)
		vi0 = sqrt(Kb*Ti/mi)
		Time_f = 1.0_mp*L/vi0

		A = (/ ve0, vi0, 0.2_mp, 1.0_mp*Ni /)
		call buildPM1D(sheath,Time_f,0.0_mp,ceiling(L/dx),2,pBC=2,mBC=2,order=1,A=A,L=L,dt=dt,eps=eps)
		sheath%wp = wp0
		call buildRecord(r,sheath%nt,2,sheath%L,sheath%ng,'test',20)

		call buildSpecies(sheath%p(1),-qe,me,n0*L/Ne)
		call buildSpecies(sheath%p(2),qe,mi,n0*L/Ni)

		call sheath_initialize(sheath,Ne,Ni,Te,Ti,Kb)
		call forwardsweep(sheath,r,Null_input,Null_source)

		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(sheath)
	end subroutine

	subroutine test_refluxing_boundary
		type(PM1D) :: reflux
		type(recordData) :: r
		integer, parameter :: Ng=64, N=10000, order=1
		real(mp) :: Ti=20, Tf = 40
		real(mp) :: xp0(N), vp0(N,3), rho_back(Ng), qe, me
		integer :: i

		call buildPM1D(reflux,Tf,Ti,Ng,1,pBC=2,mBC=2,order=order,A=(/ 1.0_mp, 1.0_mp /))
		call buildRecord(r,reflux%nt,1,reflux%L,Ng,'test_reflux',1)

		xp0 = -0.5_mp*reflux%L
		vp0 = 0.0_mp
		rho_back = 0.0_mp
		qe = -(0.1_mp)**2/(N/reflux%L)
		me = -qe
		rho_back(Ng) = -qe
		call buildSpecies(reflux%p(1),qe,me,1.0_mp)
		call setSpecies(reflux%p(1),N,xp0,vp0)
		call setMesh(reflux%m,rho_back)

		call applyBC(reflux)
		call recordPlasma(r,reflux,1)
		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(reflux)
	end subroutine

end program
