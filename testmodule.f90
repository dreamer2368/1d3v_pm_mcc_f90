module testmodule

	use init
	use timeStep

	implicit none

contains

	subroutine test_backward_sweep
		type(PM1D) :: pm
		type(adjoint) :: adj
		type(recordData) :: r
		integer, parameter :: Ng=64, Np=2, N=2
		real(mp) :: qs = 3.4_mp, ms = 2.7_mp, spwt = 1.9_mp, xp(Np), vp(Np,3)
		real(mp) :: rho_back(Ng), J0=0.0_mp,J1=0.0_mp
		real(mp) :: fxp = (0.1_mp)**9, dxp
		integer :: j
		character(len=1000)::dir1,dir2

		call buildPM1D(pm,2.0_mp,1.0_mp,Ng,N,0,0,1,eps=2.3_mp)
		dir1='test_backward_sweep/before'
		call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir1),3)
		call set_null_discharge(r)
		call buildSpecies(pm%p(1),qs,ms,spwt)
		call buildSpecies(pm%p(2),-qs,ms,spwt)

		print *, 'Original xp'
		xp = 0.4_mp*pm%L*(/ (j,j=1,Np) /)
		vp = 0.1_mp
		print *, xp
		print *, vp
		call setSpecies(pm%p(1),Np,xp,vp)
		call setSpecies(pm%p(2),Np,pm%L-xp,0.2_mp*vp)
		rho_back = 0.0_mp
		call setMesh(pm%m,rho_back)

		call forwardsweep(pm,r,Null_input,Null_source,TestQoI,J0)
		call printPlasma(r)

		call buildAdjoint(adj,pm)
		call backward_sweep(adj,pm,r,dTestQoI,Null_Dinput,Null_input,Null_source)

		print *, 'dJdxp1:', -adj%p(1)%xp/pm%dt
		print *, 'dJdxp2:', -adj%p(2)%xp/pm%dt
		print *, 'dJdvp1:', -adj%p(1)%vp(:,1)/pm%dt
		print *, 'dJdvp2:', -adj%p(2)%vp(:,1)/pm%dt

		call destroySpecies(pm%p(1))
		call destroySpecies(pm%p(2))
		call setSpecies(pm%p(2),Np,pm%L-xp,0.2_mp*vp)
		dxp = xp(2)*fxp
		xp(2) = xp(2) + dxp
		call setSpecies(pm%p(1),Np,xp,vp)
		print *, 'Perturbed xp'
		do j=1,pm%n
		print *, j,'-th species'
		print *, pm%p(j)%xp
		print *, pm%p(j)%vp(:,1)
		end do

		call destroyRecord(r)
		dir2='test_backward_sweep/after'
		call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir2),3)
		call set_null_discharge(r)
		call forwardsweep(pm,r,Null_input,Null_source,TestQoI,J1)
		call printPlasma(r)

		print *, 'J0=',J0,', J1=',J1
		print *, 'dJdxp1(2)=',(J1-J0)/dxp

		call destroyAdjoint(adj)
		call destroyPM1D(pm)
		call destroyRecord(r)
	end subroutine

	subroutine test_particle_adj(Ng,Np)
		integer, intent(in) :: Ng, Np
		type(PM1D) :: pm
		type(adjoint) :: adj
		integer :: i,j,k
		real(mp) :: qs = 3.4_mp, ms = 2.7_mp, spwt = 1.9_mp, xp(Np), vp(Np,3)
		real(mp) :: rho_back(Ng)
		real(mp) :: weight(Ng), J0, J1
		real(mp) :: rhs(Ng), rho1(Ng-1), dxp1(Np), dxp2(Np)
		real(mp) :: fxp = (0.1_mp)**9, dxp

		call buildPM1D(pm,40.0_mp,20.0_mp,Ng,1,0,0,1,eps=2.3_mp)
		call buildSpecies(pm%p(1),qs,ms,spwt)
		!particle, mesh setup
		print *, 'Original xp'
		xp = 0.4_mp*pm%L*(/ (j,j=1,Np) /)
		vp = 0.1_mp
		print *, xp
		print *, vp
		call setSpecies(pm%p(1),Np,xp,vp)
		rho_back = -qs*Np/pm%L
		call setMesh(pm%m,rho_back)
		!one time-step
		call moveSpecies(pm%p(1),pm%dt)
		call applyBC(pm)
		call assignMatrix(pm%a(1),pm%m,pm%p(1)%xp)
		call adjustGrid(pm)
		call chargeAssign(pm%a,pm%p,pm%m)
		call solveMesh(pm%m,pm%eps0)
		pm%m%E = - multiplyD(pm%m%phi,pm%m%dx,pm%m%BCindex)
		call forceAssign(pm%a(1), pm%p(1), pm%m)
		call accelSpecies(pm%p(1),pm%dt)
		!QoI
		J0 = sum( pm%p(1)%vp(:,1)**2 )
		!weight = 0.0_mp
		!weight(2*Ng/5:3*Ng/5) = 1.0_mp
		!J0 = pm%m%dx*sum( weight*pm%m%E**2 )
		print *, J0

		!Adjoint solver
		call buildAdjoint(adj,pm)
		adj%p(1)%vp(:,1) = -2.0_mp*pm%p(1)%vp(:,1)*pm%dt
		adj%p(1)%Ep = adj%p(1)%qs/pm%p(1)%ms*adj%p(1)%vp(:,1)
		call Adj_forceAssign_E(pm%a(1),adj%p(1)%Ep,adj%m%E)
		!adj%m%E = -2.0_mp*adj%m%dx*weight*pm%m%E
		call solveMesh_Adj(adj%m,pm%eps0)

		dxp1 = 0.0_mp
		dxp2 = 0.0_mp
		call Adj_chargeAssign(pm%a(1),pm%p(1),pm%m,adj%m%rho,dxp1)
		call Adj_forceAssign_xp(pm%a(1),pm%m,pm%m%E,adj%p(1)%Ep,dxp2)
		adj%p(1)%xp = - pm%dt*( dxp1 + dxp2 )

		print *, 'dJdxp'
		print *, -adj%p(1)%xp/pm%dt
		print *, 'dJdvp'
		print *, -adj%p(1)%vp(:,1)/pm%dt - adj%p(1)%xp

		!FD approximation - choose the component that you want to measure
		k = 2
		dxp = vp(k,1)*fxp
		vp(k,1) = vp(k,1) + dxp
		print *, 'Perturbed xp'
		print *, xp
		call destroySpecies(pm%p(1))
		call setSpecies(pm%p(1),Np,xp,vp)

		call moveSpecies(pm%p(1),pm%dt)
		call applyBC(pm)
		call assignMatrix(pm%a(1),pm%m,pm%p(1)%xp)
		call adjustGrid(pm)
		call chargeAssign(pm%a,pm%p,pm%m)
		call solveMesh(pm%m,pm%eps0)
		pm%m%E = - multiplyD(pm%m%phi,pm%m%dx,pm%m%BCindex)
		call forceAssign(pm%a(1), pm%p(1), pm%m)
		call accelSpecies(pm%p(1),pm%dt)
		J1 = sum( pm%p(1)%vp(:,1)**2 )
		!J1 = pm%m%dx*sum( weight*pm%m%E**2 )
		print *, 'J1 = ',J1
		print *, 'dJdxp(',k,')=', (J1-J0)/dxp

		call destroyPM1D(pm)
		call destroyAdjoint(adj)
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

!	subroutine test_anewvel_Ar
!		integer, parameter :: N = 10000
!		real(mp) :: input(3) = (/ 0.0_mp, 1.0_mp, 0.0_mp /)
!		real(mp), dimension(N,3) :: output
!		integer :: i
!
!		do i=1,N
!			output(i,:) = input
!			call anewvel_Ar(output(i,:))
!		end do
!
!		call system('mkdir -p data/scattering')
!		open(unit=301,file='data/scattering/output_Ar.bin',status='replace',form='unformatted',access='stream')
!		write(301) output
!		close(301)
!	end subroutine
!
!	subroutine test_anewvel_e
!		integer, parameter :: N = 100000
!		real(mp) :: input(3) = (/ 0.0_mp, 1.0_mp, 0.0_mp /)
!		real(mp), dimension(N,3) :: output
!		real(mp) :: energy = 100.0_mp			!eV
!		integer :: i
!
!		do i=1,N
!			output(i,:) = input
!			call anewvel_e(energy, 1.0_mp, 1.0_mp, output(i,:),.false.)
!		end do
!
!		call system('mkdir -p data/scattering')
!		open(unit=301,file='data/scattering/output.bin',status='replace',form='unformatted',access='stream')
!		write(301) output
!		close(301)
!	end subroutine
!
!	subroutine test_mcc_Argon
!		type(PM1D) :: pm
!		integer, parameter :: np = 100000, Nsample=10000
!		real(mp) :: dt = log(100.0_mp/99.0_mp)
!		real(mp) :: gden = 1.0_mp/max_sigmav_Ar, TN = 0.026_mp		!neutral temperature TN: scale in eV
!		real(mp) :: energy = 6.0_mp, vel										!argon energy
!		real(mp) :: xp0(np), vp0(np,3)
!		integer :: i, N_coll(Nsample,3)
!
!      call init_random_seed
!
!		call null_collision(gden,dt)
!		call buildPM1D(pm,30.0_mp, 15.0_mp,16,2,0,0,1,dt,L=1.0_mp)
!		call set_Ar_discharge(pm,(/1.0_mp, 1.0_mp/),(/TN,gden,0.0_mp,0.0_mp/))
!
!		vel = sqrt(2.0_mp/pm%p(2)%ms*q_e*energy)
!		xp0 = (/ (i-0.5_mp,i=1,np) /)*(1.0_mp/np)
!		vp0(:,1) = 0.0_mp
!		vp0(:,2) = vel
!		vp0(:,3) = 0.0_mp
!		call setSpecies(pm%p(2),np,xp0,vp0)
!		vp0 = 0.0_mp
!		call setSpecies(pm%p(1),np,xp0,vp0)
!
!		call system('mkdir -p data/test_mcc_argon')
!
!		call system('mkdir -p data/test_mcc_argon/before')
!		open(unit=301,file='data/test_mcc_argon/before/np.bin',status='replace',form='unformatted',access='stream')
!		open(unit=302,file='data/test_mcc_argon/before/xp_e.bin',status='replace',form='unformatted',access='stream')
!		open(unit=303,file='data/test_mcc_argon/before/vp_e.bin',status='replace',form='unformatted',access='stream')
!		open(unit=304,file='data/test_mcc_argon/before/xp_Ar.bin',status='replace',form='unformatted',access='stream')
!		open(unit=305,file='data/test_mcc_argon/before/vp_Ar.bin',status='replace',form='unformatted',access='stream')
!		write(301) np
!		write(302) pm%p(1)%xp
!		write(303) pm%p(1)%vp
!		write(304) pm%p(2)%xp
!		write(305) pm%p(2)%vp
!		close(301)
!		close(302)
!		close(303)
!		close(304)
!		close(305)
!
!		call mcc_argon(pm)
!
!		call system('mkdir -p data/test_mcc_argon/after')
!		open(unit=301,file='data/test_mcc_argon/after/np_e.bin',status='replace',form='unformatted',access='stream')
!		open(unit=302,file='data/test_mcc_argon/after/xp_e.bin',status='replace',form='unformatted',access='stream')
!		open(unit=303,file='data/test_mcc_argon/after/vp_e.bin',status='replace',form='unformatted',access='stream')
!		open(unit=304,file='data/test_mcc_argon/after/xp_Ar.bin',status='replace',form='unformatted',access='stream')
!		open(unit=305,file='data/test_mcc_argon/after/vp_Ar.bin',status='replace',form='unformatted',access='stream')
!		open(unit=306,file='data/test_mcc_argon/after/np_Ar.bin',status='replace',form='unformatted',access='stream')
!		write(301) pm%p(1)%np
!		write(302) pm%p(1)%xp
!		write(303) pm%p(1)%vp
!		write(304) pm%p(2)%xp
!		write(305) pm%p(2)%vp
!		write(306) pm%p(2)%np
!		close(301)
!		close(302)
!		close(303)
!		close(304)
!		close(305)
!		close(306)
!
!		call destroyPM1D(pm)
!
!		call null_collision(gden,dt)
!		vel = sqrt(2.0_mp/m_Ar*q_e*energy)
!		xp0 = (/ (i-0.5_mp,i=1,np) /)*(1.0_mp/np)
!		vp0(:,1) = 0.0_mp
!		vp0(:,2) = vel
!		vp0(:,3) = 0.0_mp
!
!		open(unit=301,file='data/test_mcc_Argon/prob.bin',status='replace',form='unformatted',access='stream')
!		write(301)	col_prob_Ar,	&
!                 col_prob_Ar*( asigma4(energy)*vel/max_sigmav_Ar ),   &
!                 col_prob_Ar*( asigma5(energy)*vel/max_sigmav_Ar )
!		close(301)
!
!		do i=1,Nsample
!         print *, ' '
!         print *, '================',i,'-th Sample=================='
!         print *, ' '
!   		call buildPM1D(pm,30.0_mp, 15.0_mp,16,2,0,0,1,dt,L=1.0_mp)
!   		call set_Ar_discharge(pm,(/1.0_mp, 1.0_mp/),(/TN,gden,0.0_mp,0.0_mp/))
!			call setSpecies(pm%p(1),np,xp0,vp0)
!			call setSpecies(pm%p(2),np,xp0,vp0)
!         pm%p(1)%vp = 0.0_mp
!
!			call mcc_argon(pm,N_coll(i,:))
!
!			call destroyPM1D(pm)
!		end do
!
!		open(unit=301,file='data/test_mcc_Argon/coll_sample.bin',status='replace',form='unformatted',access='stream')
!      write(301) N_coll
!      close(301)
!	end subroutine
!
!	subroutine test_mcc_electron
!		type(PM1D) :: pm
!		integer, parameter :: np = 100000, Nsample = 10000
!		real(mp) :: dt = log(100.0_mp/99.0_mp)
!		real(mp) :: gden = 1.0_mp/max_sigmav_e, TN = 0.026_mp		!neutral temperature TN: scale in eV
!		real(mp) :: energy = 20.0_mp, vel										!electron energy
!		real(mp) :: xp0(np), vp0(np,3)
!      real(mp) :: rnd, pr(3)
!		integer :: i, j, N_coll(Nsample,4), N_exp(Nsample,3)
!
!      call init_random_seed
!
!		call null_collision(gden,dt)
!		call buildPM1D(pm,30.0_mp, 15.0_mp,16,2,0,0,1,dt,L=1.0_mp)
!		call set_Ar_discharge(pm,(/1.0_mp, 1.0_mp/),(/TN,gden,0.0_mp,0.0_mp/))
!
!		vel = sqrt(2.0_mp/pm%p(1)%ms*q_e*energy)
!		xp0 = (/ (i-0.5_mp,i=1,np) /)*(1.0_mp/np)
!		vp0(:,1) = 0.0_mp
!		vp0(:,2) = vel
!		vp0(:,3) = 0.0_mp
!		call setSpecies(pm%p(1),np,xp0,vp0)
!		vp0 = 0.0_mp
!		call setSpecies(pm%p(2),np,xp0,vp0)
!
!		call system('mkdir -p data/test_mcc_electron/before')
!		open(unit=301,file='data/test_mcc_electron/before/np.bin',status='replace',form='unformatted',access='stream')
!		open(unit=302,file='data/test_mcc_electron/before/xp_e.bin',status='replace',form='unformatted',access='stream')
!		open(unit=303,file='data/test_mcc_electron/before/vp_e.bin',status='replace',form='unformatted',access='stream')
!		open(unit=304,file='data/test_mcc_electron/before/xp_Ar.bin',status='replace',form='unformatted',access='stream')
!		open(unit=305,file='data/test_mcc_electron/before/vp_Ar.bin',status='replace',form='unformatted',access='stream')
!		write(301) np
!		write(302) pm%p(1)%xp
!		write(303) pm%p(1)%vp
!		write(304) pm%p(2)%xp
!		write(305) pm%p(2)%vp
!		close(301)
!		close(302)
!		close(303)
!		close(304)
!		close(305)
!
!		call mcc_electron(pm)
!
!		call system('mkdir -p data/test_mcc_electron/after')
!		open(unit=301,file='data/test_mcc_electron/after/np_e.bin',status='replace',form='unformatted',access='stream')
!		open(unit=302,file='data/test_mcc_electron/after/xp_e.bin',status='replace',form='unformatted',access='stream')
!		open(unit=303,file='data/test_mcc_electron/after/vp_e.bin',status='replace',form='unformatted',access='stream')
!		open(unit=304,file='data/test_mcc_electron/after/xp_Ar.bin',status='replace',form='unformatted',access='stream')
!		open(unit=305,file='data/test_mcc_electron/after/vp_Ar.bin',status='replace',form='unformatted',access='stream')
!		open(unit=306,file='data/test_mcc_electron/after/np_Ar.bin',status='replace',form='unformatted',access='stream')
!		write(301) pm%p(1)%np
!		write(302) pm%p(1)%xp
!		write(303) pm%p(1)%vp
!		write(304) pm%p(2)%xp
!		write(305) pm%p(2)%vp
!		write(306) pm%p(2)%np
!		close(301)
!		close(302)
!		close(303)
!		close(304)
!		close(305)
!		close(306)
!		call destroyPM1D(pm)
!
!		call null_collision(gden,dt)
!		vel = sqrt(2.0_mp/m_e*q_e*energy)
!		xp0 = (/ (i-0.5_mp,i=1,np) /)*(1.0_mp/np)
!		vp0(:,1) = 0.0_mp
!		vp0(:,2) = vel
!		vp0(:,3) = 0.0_mp
!
!      pr = (/ asigma1(energy), asigma2(energy), asigma3(energy) /)*vel/max_sigmav_e
!		open(unit=301,file='data/test_mcc_electron/prob.bin',status='replace',form='unformatted',access='stream')
!		write(301)	col_prob_e,	&
!                 col_prob_e*pr(1),   &
!                 col_prob_e*pr(2),   &
!                 col_prob_e*pr(3)
!		close(301)
!
!      N_exp = 0
!		do i=1,Nsample
!         print *, ' '
!         print *, '================',i,'-th Sample=================='
!         print *, ' '
!   		call buildPM1D(pm,30.0_mp, 15.0_mp,16,2,0,0,1,dt,L=1.0_mp)
!   		call set_Ar_discharge(pm,(/1.0_mp, 1.0_mp/),(/TN,gden,0.0_mp,0.0_mp/))
!			call setSpecies(pm%p(1),np,xp0,vp0)
!			call setSpecies(pm%p(2),np,xp0,vp0)
!         pm%p(2)%vp = 0.0_mp
!
!			call mcc_electron(pm,N_coll(i,:))
!
!			call destroyPM1D(pm)
!
!         !Do the simple sampling comparison, if needed
!!         do j=1,N_coll(i,1)
!!            call RANDOM_NUMBER(rnd)
!!            if( rnd.le.pr(1) ) then
!!               N_exp(i,1) = N_exp(i,1)+1
!!            elseif( rnd.le.(pr(1)+pr(2)) ) then
!!               N_exp(i,2) = N_exp(i,2)+1
!!            elseif( rnd.le.(pr(1)+pr(2)+pr(3) ) then
!!               N_exp(i,3) = N_exp(i,3)+1
!!            end if
!!         end do
!		end do
!
!		open(unit=301,file='data/test_mcc_electron/coll_sample.bin',status='replace',form='unformatted',access='stream')
!      write(301) N_coll
!      close(301)
!	end subroutine
!
!	subroutine test_ext_voltage_Poisson
!		type(mesh) :: m
!		integer, parameter :: N=128
!		real(mp) :: xg(N), rho_back(N), rho(N), sol(N)
!		integer :: i
!
!		xg = (/ (i-1, i=1,N) /)*1.0_mp/(N-1)
!		rho = exp(3.0_mp*xg)
!		rho_back = 0.0_mp
!		rho_back(N) = 1.0_mp
!
!		sol = ( 1.0_mp-exp(3.0_mp*xg) )/9.0_mp + ( 1.0_mp - (1.0_mp-exp(3.0_mp))/9.0_mp )*xg
!
!		call buildMesh(m,1.0_mp,N,1)
!		call setMesh(m,rho_back)
!		m%rho = rho
!		call solveMesh(m,1.0_mp)
!
!		call system('mkdir -p data/test_DD_poisson')
!		open(unit=301,file='data/test_DD_poisson/xg.bin',status='replace',form='unformatted',access='stream')
!		open(unit=302,file='data/test_DD_poisson/phi.bin',status='replace',form='unformatted',access='stream')
!		open(unit=303,file='data/test_DD_poisson/phi_sol.bin',status='replace',form='unformatted',access='stream')
!		write(301) xg
!		write(302) m%phi
!		write(303) sol
!		close(301)
!		close(302)
!		close(303)
!
!		print *, 'error: ', maxval( abs(sol - m%phi) )
!	end subroutine

end module