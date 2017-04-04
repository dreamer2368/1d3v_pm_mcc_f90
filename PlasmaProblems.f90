module PlasmaProblems

	use init
	use timeStep
	use modMPI

	implicit none

contains

	subroutine twostream
		type(PM1D) :: pm
		type(adjoint) :: adj
		type(recordData) :: r
		integer, parameter :: Ng=64, Np=10**5, N=1
		real(mp) :: v0 = 0.2_mp, vT = 0.0_mp
		integer :: mode=1
		real(mp) :: J0,J1, grad(1)
		character(len=100)::dir1

		call buildPM1D(pm,70.0_mp,65.0_mp,Ng,N,0,0,1,A=(/1.0_mp,0.0_mp/))
		dir1='twostream_test'
		call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir1),5)
		call set_null_discharge(r)
		call twostream_initialize(pm,Np,v0,vT,mode)
		call forwardsweep(pm,r,Null_input,Null_source,MKE,J0)
		call printPlasma(r)
		print *, 'J0=',J0

		call destroyPM1D(pm)
		call destroyRecord(r)
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

		call buildSpecies(sheath%p(1),-qe,me)
		call buildSpecies(sheath%p(2),qe,mi)

		call sheath_initialize(sheath,Ne,Ni,Te,Ti,Kb,n0)
		call forwardsweep(sheath,r,Null_input,Null_source)

		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(sheath)
	end subroutine

	subroutine debye_shielding
		type(PM1D) :: d
		type(recordData) :: r
		real(mp) :: n0 = 1.0e10, lambda0 = 1.0e-2, vT = 1.0_mp
		integer :: N = 100000, Ng = 128
		real(mp) :: L = 20.0_mp, Wp, Q = 2.0_mp
		real(mp) :: dt = 0.05_mp
		real(mp) :: Time = 150.0_mp
		real(mp) :: A(2),J

		A = (/ vT, lambda0 /)
		call buildPM1D(d,Time,0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=dt)
		call buildRecord(r,d%nt,1,d%L,d%ng,'debye',10)

		call buildSpecies(d%p(1),-1.0_mp,1.0_mp)
		call Debye_initialize(d,N,Q)

		call forwardsweep(d,r,Null_input,Null_source,Debye,J)

		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(d)
	end subroutine

end module
