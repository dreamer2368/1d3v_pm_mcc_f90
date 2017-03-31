program main

	use testmodule

	implicit none

	real(mp) :: output(2) = (/(0.1_mp)**3,0.0_mp/)

	! print to screen
	print *, 'calling program main'

!	call cross_section
!	call Procassini
!	call test_refluxing_boundary
!	call test_anewvel_Ar
!	call test_mcc_electron
!	call test_mcc_Argon
!	call test_ext_voltage_Poisson
!   call Ar_discharge
!	call test_particle_adj(64,2)
!	call test_backward_sweep
!	call twostream_adj(output(1),output(2))
!	call Landau(0.0_mp, 6.0_mp ,'Landau', 1,output )
!	call random_test
!   call Landau_adjoint_sampling
!   call twostream_adjoint_sampling
!	call debye_shielding
!	call debye_characterization
!	call InjectionTest
!	call MPITest
!	call SensitivityInitializeTest
	call Debye_sensitivity
!	call forYeoh
!	call RedistributionTest
!	call updateWeightTest
!	call debye_sensitivity_curve
!	call adj_convergence(twostream_grad)
!	call adjoint_convergence_in_time(debye_adj)
!	call debye_sampling

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

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

	subroutine twostream_adjoint_sampling
		integer :: ierr, my_rank, s
		type(adjoint) :: adj
		type(PM1D) :: pm
		type(recordData) :: r
		integer, parameter :: Nsample = 10000
		integer :: sample_per_core, sendcnt
		integer, allocatable :: recvcnt(:), displc(:)
		real(mp), allocatable :: sendbuf(:,:)
		real(mp) :: recvbuf(Nsample,3)                     !(/J0, J1, dJdA/)
		real(mp) :: Tf = 30.1_mp, Ti = 30.0_mp
		integer, parameter :: Ng=64, Np=10**5, N=1
		real(mp) :: xp0(Np), vp0(Np,3), spwt0(Np)
		real(mp) :: v0 = 0.2_mp, vT = 0.01_mp
		integer :: mode = 1
		real(mp) :: dt=0.1_mp
		character(len=100)::dir,rank_str
		real(mp) :: J0,J1,grad(1)
		integer :: i

		call MPI_INIT(ierr)
		call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
		call MPI_COMM_SIZE(MPI_COMM_WORLD,s,ierr)
		sample_per_core = Nsample/s
		write(rank_str,*) my_rank

		if( my_rank.eq.s-1 ) then
			print *, 'size: ',s
			print *, 'sample/core: ',sample_per_core
			print *, 'remainder: ',MOD(Nsample,s)
			allocate(recvcnt(0:s-1))
			allocate(displc(0:s-1))
			recvcnt(0:MOD(Nsample,s)-1) = sample_per_core+1
			recvcnt(MOD(Nsample,s):s-1) = sample_per_core
			print *, recvcnt
			displc = 0
			do i=0,s-1
				displc(i) = SUM(recvcnt(0:i-1))
			end do
			print *, displc
		end if

		if(my_rank<MOD(Nsample,s) ) then
			sendcnt = sample_per_core+1
			allocate(sendbuf(sendcnt,3))
		else
			sendcnt = sample_per_core
			allocate(sendbuf(sendcnt,3))
		end if
		sendbuf = 0.0_mp

		call init_random_seed(my_rank)

		do i=1,sendcnt
			call buildPM1D(pm,Tf,Ti,Ng,N,0,0,1,dt=dt,A=(/0.1_mp,0.0_mp/))
			dir = 'twostream_sampling2/before'//trim(adjustl(rank_str))
			call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir),10)
			call set_null_discharge(r)
			call twostream_initialize(pm,Np,v0,vT,mode)
			xp0 = pm%p(1)%xp
			vp0 = pm%p(1)%vp
			spwt0 = pm%p(1)%spwt

			call forwardsweep(pm,r,Te,Null_source,MPE,J0)
			!call printPlasma(r)
			print *, 'J0=',J0

			call buildAdjoint(adj,pm)
			call backward_sweep(adj,pm,r,grad,dMPE,dTe,dTedA,Te,Null_source)

			print *, 'dJdA=',grad

			call destroyAdjoint(adj)
			call destroyRecord(r)

			dir = 'twostream_sampling2/after'//trim(adjustl(rank_str))
			pm%A0 = (/0.1_mp,(0.1_mp)**10 /)
			call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir),10)
			call destroySpecies(pm%p(1))
			call setSpecies(pm%p(1),Np,xp0,vp0,spwt0)
			call forwardsweep(pm,r,Te,Null_source,MPE,J1)
			!call printPlasma(r)
			print *, 'J1=',J1

			call destroyRecord(r)
			call destroyPM1D(pm)

			sendbuf(i,:) = (/ J0, J1, grad(1) /)
		end do

		do i=1,3
			call MPI_GATHERV(sendbuf(:,i),sendcnt,MPI_DOUBLE,recvbuf(:,i),recvcnt,displc,MPI_DOUBLE,s-1,MPI_COMM_WORLD,ierr)
		end do

		call MPI_FINALIZE(ierr)
		if( my_rank.eq.s-1 ) then
			open(unit=301,file='data/twostream_sampling2/sample_J0.bin',status='replace',form='unformatted',access='stream')
			open(unit=302,file='data/twostream_sampling2/sample_J1.bin',status='replace',form='unformatted',access='stream')
			open(unit=303,file='data/twostream_sampling2/sample_grad.bin',status='replace',form='unformatted',access='stream')
			write(301) recvbuf(:,1)
			write(302) recvbuf(:,2)
			write(303) recvbuf(:,3)
			close(301)
			close(302)
			close(303)
		end if

		deallocate(sendbuf)
		if( my_rank.eq.s-1) then
			deallocate(recvcnt)
			deallocate(displc)
		end if
	end subroutine

	subroutine Landau_adjoint_sampling
		integer :: ierr, my_rank, s
		type(adjoint) :: adj
		type(PM1D) :: pm
		type(recordData) :: r
		integer, parameter :: Nsample = 10000
		integer :: sample_per_core, sendcnt
		integer, allocatable :: recvcnt(:), displc(:)
		real(mp), allocatable :: sendbuf(:,:)
		real(mp) :: recvbuf(Nsample,3)                     !(/J0, J1, dJdA/)
		real(mp) :: Tf = 20.1_mp, Ti = 20.0_mp
		integer, parameter :: Ng=64, Np=3*10**6, N=1
		real(mp) :: xp0(Np), vp0(Np,3), spwt0(Np)
		real(mp) :: vT = 1.0_mp, L=4.0_mp*pi
		real(mp) :: dt=0.1_mp
		character(len=100)::dir,rank_str
		real(mp) :: J0,J1,grad(1)
		integer :: i

		call MPI_INIT(ierr)
		call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
		call MPI_COMM_SIZE(MPI_COMM_WORLD,s,ierr)
		sample_per_core = Nsample/s
		write(rank_str,*) my_rank

		if( my_rank.eq.s-1 ) then
		   print *, 'size: ',s
		   print *, 'sample/core: ',sample_per_core
		   print *, 'remainder: ',MOD(Nsample,s)
			allocate(recvcnt(0:s-1))
		   allocate(displc(0:s-1))
			recvcnt(0:MOD(Nsample,s)-1) = sample_per_core+1
			recvcnt(MOD(Nsample,s):s-1) = sample_per_core
		   print *, recvcnt
		   displc = 0
		   do i=0,s-1
			  displc(i) = SUM(recvcnt(0:i-1))
		   end do
		   print *, displc
		end if

		if(my_rank<MOD(Nsample,s) ) then
			sendcnt = sample_per_core+1
			allocate(sendbuf(sendcnt,3))
		else
			sendcnt = sample_per_core
			allocate(sendbuf(sendcnt,3))
		end if
		sendbuf = 0.0_mp

		call init_random_seed(my_rank)

		do i=1,sendcnt
			call buildPM1D(pm,Tf,Ti,Ng,N,0,0,1,dt=dt,L=L,A=(/0.1_mp,0.0_mp/))
			dir = 'Landau_sampling/before'//trim(adjustl(rank_str))
			call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir),10)
			call set_null_discharge(r)
			call Landau_initialize(pm,Np,vT)
			xp0 = pm%p(1)%xp
			vp0 = pm%p(1)%vp
			spwt0 = pm%p(1)%spwt

			call forwardsweep(pm,r,Te,Null_source,MPE,J0)
			!call printPlasma(r)
			print *, 'J0=',J0

			call buildAdjoint(adj,pm)
			call backward_sweep(adj,pm,r,grad,dMPE,dTe,dTedA,Te,Null_source)

			print *, 'dJdA=',grad

			call destroyAdjoint(adj)
			call destroyRecord(r)

			dir = 'Landau_sampling/after'//trim(adjustl(rank_str))
			pm%A0 = (/0.1_mp,(0.1_mp)**9 /)
			call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir),10)
		   call destroySpecies(pm%p(1))
			call setSpecies(pm%p(1),Np,xp0,vp0,spwt0)
			call forwardsweep(pm,r,Te,Null_source,MPE,J1)
			!call printPlasma(r)
			print *, 'J1=',J1

			call destroyRecord(r)
			call destroyPM1D(pm)

			sendbuf(i,:) = (/ J0, J1, grad(1) /)
		end do

		do i=1,3
			call MPI_GATHERV(sendbuf(:,i),sendcnt,MPI_DOUBLE,recvbuf(:,i),recvcnt,displc,MPI_DOUBLE,s-1,MPI_COMM_WORLD,ierr)
		end do

		call MPI_FINALIZE(ierr)
		if( my_rank.eq.s-1 ) then
			open(unit=301,file='data/Landau_sampling/sample_J0.bin',status='replace',form='unformatted',access='stream')
			open(unit=302,file='data/Landau_sampling/sample_J1.bin',status='replace',form='unformatted',access='stream')
			open(unit=303,file='data/Landau_sampling/sample_grad.bin',status='replace',form='unformatted',access='stream')
		   write(301) recvbuf(:,1)
		   write(302) recvbuf(:,2)
		   write(303) recvbuf(:,3)
		   close(301)
		   close(302)
		   close(303)
		end if

		deallocate(sendbuf)
		if( my_rank.eq.s-1) then
		   deallocate(recvcnt)
		   deallocate(displc)
		end if
	end subroutine

	subroutine adjoint_convergence_in_time(problem)
		integer, parameter :: N=23, Nt=1
		real(mp) :: fk(N)
		real(mp) :: Tk(Nt)
		real(mp) :: ek(N,Nt)
		character(len=100) :: dir
		integer :: i,j
		real(mp) :: J0,J1,grad(Nt),temp(2)
		interface
			subroutine problem(fk,Ti,str,k,output)
				use modPM1D
				use modAdj
				use modRecord
				real(mp), intent(in) :: fk, Ti
				character(len=*), intent(in) ::str
				integer, intent(in) :: k
				real(mp), intent(out) :: output(:)
				type(adjoint) :: adj
				type(PM1D) :: pm
				type(recordData) :: r
			end subroutine
		end interface

		fk = (/ (EXP(-1.0_mp*(i-1)),i=1,N) /)
		Tk = (/ 3.0_mp /)
		ek = 0.0_mp

		dir = 'debye_adj_test'

		do j=1,Nt
			call problem(fk(i),Tk(j),trim(dir),0,temp)
			J0 = temp(1)
			grad(j) = temp(2)
			do i=1,N
				call problem(fk(i),Tk(j),trim(dir),1,temp)
				J1 = temp(1)
				ek(i,j) = ABS( ((J1-J0)/fk(i) - grad(j))/grad(j) )
			end do
		end do

		open(unit=301,file='data/'//trim(dir)//'/grad_convergence.dat',status='replace')
		do i=1,N
			write(301,*) fk(i), ek(i,:)
		end do
		close(301)
		open(unit=301,file='data/'//trim(dir)//'/grad_in_time.bin',status='replace')
		do i=1,Nt
			write(301,*) Tk(i),grad(i)
		end do
		close(301)
	end subroutine

	subroutine adj_convergence(problem)
		interface
			subroutine problem(fk,Time,str,k,output)
				use modPM1D
				use modAdj
				use modRecord
				real(mp), intent(in) :: fk
				real(mp), intent(in) :: Time
				character(len=*), intent(in) :: str
				integer, intent(in) :: k
				real(mp), intent(out) :: output(:)
			end subroutine
		end interface
		character(len=100) :: dir
		integer, parameter :: N=20
		real(mp) :: temp(2), J0, J1, grad
		real(mp), dimension(N) :: fk,ek
		real(mp) :: Time=0.1_mp
		integer :: i
		dir = 'adj_test'

		fk = (/ (EXP(-1.0_mp*(i-1)),i=1,N) /)
		ek = 0.0_mp

		call problem(fk(i),Time,trim(dir),0,temp)
		J0 = temp(1)
		grad = temp(2)
		do i=1,N
			call problem(fk(i),Time,trim(dir),1,temp)
			J1 = temp(1)
			ek(i) = ABS( ((J1-J0)/fk(i) - grad)/grad )
		end do

		open(unit=301,file='data/'//trim(dir)//'/adj_convergence.dat',status='replace')
		do i=1,N
			write(301,*) fk(i), ek(i)
		end do
		close(301)
	end subroutine

   subroutine Ar_discharge
      type(PM1D) :: pm
      type(recordData) :: r
      real(mp), parameter :: Kb = 1.38065E-23, EV_TO_K = 11604.52_mp, eps = 8.85418782E-12, mTorr_to_Pa = 0.13332237_mp
      real(mp) :: I0 = 25.6_mp, I_f = 13.56E6                !I0(A/m2), If(Hz)
      real(mp) :: TN = 0.026_mp, PN = 50.0_mp, gden        !TN(eV), PN(mTorr), gden(m-3)
	  real(mp), parameter :: n0=10.0_mp**15!, v0_e = 5.9E5, v0_Ar = 2.19E3
	  real(mp) :: T0=1.0_mp, wp0, lambda0, v0_e, v0_Ar
	  real(mp), parameter :: L = 0.02_mp, area = 0.016_mp               !L(m), area(m2)
	  integer, parameter :: nc2p = 10**6, Np=CEILING(n0*L*area/nc2p), Ng = 300
	  real(mp) :: spwt(Np), xp0(Np), vp0(Np,3)
	  real(mp) :: dt
	  integer :: i
	  gden = (PN*mTorr_to_Pa)/(q_e*TN)

	  spwt = n0*L/Np                      !spwt = nc2p/area

	  print *, 'gden(m-3): ',gden,', n0(m-3): ',n0,', spwt(m-2): ',spwt

!      T0 = 0.5_mp*m_e*v0_e**2/q_e*EV_TO_K
	  v0_e = sqrt(2.0_mp*T0/EV_TO_K*q_e/m_e)
	  v0_Ar = sqrt(2.0_mp*T0/EV_TO_K*q_e/m_Ar)
	  wp0 = sqrt(n0*q_e*q_e/m_e/eps)
	  lambda0 = sqrt(eps*T0/n0/q_e)

	  print *, 'L = ',L,', lambda0 = ',lambda0,' e = lambda/L = ',lambda0/L

	  dt = 7.20179E-11
	  print *, 'dt = ',dt,', wp*dt=',wp0*dt
	  call init_random_seed
	  call null_collision(gden,dt)
	  print *, 'P_e = ',col_prob_e,', P_Ar = ', col_prob_Ar

	  call buildPM1D(pm,40000.0_mp*dt,100.0_mp*dt,Ng,2,pBC=1,mBC=2,order=1,dt=dt,L=L,eps=eps)
	  call buildRecord(r,pm%nt,2,pm%L,pm%ng,'rf_Ar4',20)
	  open(unit=301,file='data/rf_Ar4/input',status='replace',form='unformatted',access='stream')
	  write(301) TN, PN, v0_e, v0_Ar, I0, L, area, dt
	  close(301)

	  call set_Ar_discharge(pm,(/TN,gden,I0,I_f/),r)
	  call RANDOM_NUMBER(xp0)
	  vp0 = randn(Np,3)*v0_e
	  call setSpecies(pm%p(1),Np,xp0*L,vp0,spwt)
	  call RANDOM_NUMBER(xp0)
	  vp0 = randn(Np,3)*v0_Ar
	  call setSpecies(pm%p(2),Np,xp0*L,vp0,spwt)

      call forwardsweep(pm,r,RF_current,Null_source)

      call printPlasma(r)

      call destroyPM1D(pm)
      call destroyRecord(r)
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
		real(mp) :: Time = 30.0_mp
		real(mp) :: A(2),J

		A = (/ vT, lambda0 /)
		call buildPM1D(d,Time,0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=dt)
		call buildRecord(r,d%nt,1,d%L,d%ng,'debye',20)

		call buildSpecies(d%p(1),-1.0_mp,1.0_mp)
		call Debye_initialize(d,N,Q)

		call forwardsweep(d,r,Null_input,Null_source,Debye,J)

		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(d)
	end subroutine

	subroutine debye_characterization
		type(PM1D) :: d
		type(recordData) :: r
		type(mpiHandler) :: mpih
		real(mp) :: vT(1001)
		integer :: N = 100000, Ng = 512
		real(mp) :: L = 20.0_mp, Wp, Q = 2.0_mp
		real(mp) :: dx
		real(mp) :: Time = 30.0_mp
		real(mp) :: A(2),J
		integer :: i
		character(len=100)::dir
		vT = (/ (3.0_mp*(i-1)/(1001-1)+0.5_mp,i=1,1001) /)

		call buildMPIHandler(mpih)
		call allocateBuffer(1001,2,mpih)

		call init_random_seed
		do i=1,mpih%sendcnt
			A = (/ vT(mpih%displc(mpih%my_rank)+i), 0.0_mp /)
			call buildPM1D(d,Time,0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=0.01_mp)
			dir = 'Debye_characterization/'//trim(adjustl(mpih%rank_str))
			call buildRecord(r,d%nt,1,d%L,d%ng,trim(dir),20)

			call buildSpecies(d%p(1),-1.0_mp,1.0_mp)
			call Debye_initialize(d,N,Q)

			call forwardsweep(d,r,Null_input,Null_source,Debye,J)

			mpih%sendbuf(i,:) = (/vT(mpih%displc(mpih%my_rank)+i),J/)

			call destroyRecord(r)
			call destroyPM1D(d)
		end do

		call gatherData(mpih)

		if( mpih%my_rank.eq.mpih%size-1 ) then
			open(unit=301,file='data/Debye_characterization/Ak.bin',status='replace',form='unformatted',access='stream')
			open(unit=302,file='data/Debye_characterization/Jk.bin',status='replace',form='unformatted',access='stream')
		   write(301) mpih%recvbuf(:,1)
		   write(302) mpih%recvbuf(:,2)
		   close(301)
		   close(302)
		end if

		call destroyMPIHandler(mpih)
	end subroutine

	subroutine Debye_sensitivity
		type(PM1D) :: pm, dpm
		type(recordData) :: r, fsr
		integer :: N=1E5, Ng=64
		integer :: NInit=5E4, Ngv=32, NInject=5E3, NLimit=3E5
		real(mp) :: L = 20.0_mp, Lv, Q = 2.0_mp
		real(mp) :: dt=0.05_mp, dx
		real(mp) :: Time = 3.0_mp, vT = 1.0_mp
		real(mp) :: A(2), J, grad
		character(len=100)::dir
		A = (/ vT, 0.0_mp /)
		Lv = vT*6.0_mp

		call buildPM1D(pm,Time,0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=dt,Ngv=Ngv,Lv=Lv)
		dir = 'Debye_sensitivity'
		call buildRecord(r,pm%nt,1,pm%L,pm%ng,trim(dir),1)

		call buildSpecies(pm%p(1),-1.0_mp,1.0_mp)
		call Debye_initialize(pm,N,Q)
		
!		NInject = 5*N/pm%nt
		call buildSensitivity(dpm,pm)
		dir = 'Debye_sensitivity/f_A'
		call buildRecord(fsr,dpm%nt,1,dpm%L,dpm%ng,trim(dir),1)
		call Debye_sensitivity_init(dpm,N,vT,'vT')
!		call Debye_sensitivity_init_sync(dpm,pm,vT,'vT')

		call forwardsweep_sensitivity(pm,r,dpm,fsr,NInject,NLimit,Null_input,Null_source,Debye,J,grad)

		print *, 'J: ',J
		print *, 'grad: ',grad

		call printPlasma(r)
		call printPlasma(fsr)

		call destroyPM1D(pm)
		call destroyRecord(r)
		call destroyPM1D(dpm)
		call destroyRecord(fsr)
	end subroutine

	subroutine debye_sensitivity_curve
		type(PM1D) :: d,dpm
		type(recordData) :: r,fsr
		type(mpiHandler) :: mpih
		integer, parameter  :: Nsample=101
		real(mp) :: vT(Nsample)
		integer :: N = 100000, Ng = 64
		integer :: NInit=5E4, Ngv=32, NInject=5E3, NLimit=3E5
		real(mp) :: L = 20.0_mp, Lv, Q = 2.0_mp
		real(mp) :: dt=0.05_mp, dx
		real(mp) :: Time = 750.0_mp
		real(mp) :: A(2),J,grad
		integer :: i
		character(len=100)::dir
		vT = (/ (1.0_mp*(i-1)/(Nsample-1)+1.0_mp,i=1,Nsample) /)

		call buildMPIHandler(mpih)
		call allocateBuffer(Nsample,3,mpih)

		call init_random_seed
		do i=1,mpih%sendcnt
			A = (/ vT(mpih%displc(mpih%my_rank)+i), 0.0_mp /)
			Lv = A(1)*6.0_mp

			print *, 'vT(',mpih%my_rank,')=',A(1)
			call buildPM1D(d,Time,0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=dt,Ngv=Ngv,Lv=Lv)
			dir = 'Debye_sensitivity_curve/'//trim(adjustl(mpih%rank_str))
			call buildRecord(r,d%nt,1,d%L,d%ng,trim(dir),20)

			call buildSpecies(d%p(1),-1.0_mp,1.0_mp)
			call Debye_initialize(d,N,Q)
		
			call buildSensitivity(dpm,d)
			dir = 'Debye_sensitivity_curve/'//trim(adjustl(mpih%rank_str))//'/f_A'
			call buildRecord(fsr,dpm%nt,1,dpm%L,dpm%ng,trim(dir),20)
			call Debye_sensitivity_init(dpm,2*N,A(1))
	
			call forwardsweep_sensitivity(d,r,dpm,fsr,NInject,NLimit,Null_input,Null_source,Debye,J,grad)

			mpih%sendbuf(i,:) = (/vT(mpih%displc(mpih%my_rank)+i),J,grad/)

			call destroyRecord(r)
			call destroyPM1D(d)
			call destroyRecord(fsr)
			call destroyPM1D(dpm)
		end do

		call gatherData(mpih)

		if( mpih%my_rank.eq.mpih%size-1 ) then
			open(unit=301,file='data/Debye_sensitivity_curve/Ak.bin',status='replace',form='unformatted',access='stream')
			open(unit=302,file='data/Debye_sensitivity_curve/Jk.bin',status='replace',form='unformatted',access='stream')
			open(unit=303,file='data/Debye_sensitivity_curve/gradk.bin',status='replace',form='unformatted',access='stream')
		   write(301) mpih%recvbuf(:,1)
		   write(302) mpih%recvbuf(:,2)
			write(303) mpih%recvbuf(:,3)
		   close(301)
		   close(302)
			close(303)
		end if

		call destroyMPIHandler(mpih)
	end subroutine

	subroutine debye_sampling
		type(PM1D) :: d, dpm
		type(adjoint) :: adj
		type(recordData) :: r,fsr
		type(mpiHandler) :: mpih
		integer, parameter  :: Nsample=5
		real(mp), parameter :: vT=1.5_mp
		integer, parameter :: N = 100000, Ng = 64
		integer, parameter :: NInit=5E4, Ngv=32, NInject=5E3, NLimit=3E5
		real(mp) :: L = 20.0_mp, Lv, Q = 2.0_mp
		real(mp) :: dt=0.05_mp, dx
		real(mp) :: Time(4), A(2)
		real(mp) :: J,grad,adj_grad(1)
		integer :: i,k
		character(len=100)::dir,Time_str
		Time = (/ 0.1_mp, 0.2_mp, 0.3_mp, 0.4_mp /)
		A = (/ vT, 0.0_mp /)
		Lv = A(1)*6.0_mp

		call buildMPIHandler(mpih)
		call allocateBuffer(Nsample,3,mpih)

		call init_random_seed(mpih%my_rank)
		do k=1,size(Time)
			J = 0.0_mp
			grad = 0.0_mp
			adj_grad = 0.0_mp
			do i=1,mpih%sendcnt
				call buildPM1D(d,Time(k),0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=dt,Ngv=Ngv,Lv=Lv)
				dir = 'Debye_sampling/'//trim(adjustl(mpih%rank_str))
				call buildRecord(r,d%nt,1,d%L,d%ng,trim(dir),20)

				call buildSpecies(d%p(1),-1.0_mp,1.0_mp)
				call Debye_initialize(d,N,Q)
		
				call buildSensitivity(dpm,d)
				dir = 'Debye_sampling/'//trim(adjustl(mpih%rank_str))//'/f_A'
				call buildRecord(fsr,dpm%nt,1,dpm%L,dpm%ng,trim(dir),20)
				call Debye_sensitivity_init(dpm,2*N,A(1))
	
				call forwardsweep_sensitivity(d,r,dpm,fsr,NInject,NLimit,Null_input,Null_source,Debye,J,grad)

				call buildAdjoint(adj,d)
				call adj%m%setMesh(d%m%rho_back)
				call backward_sweep(adj,d,r,adj_grad,dDebye,Null_dinput,dDebye_dvT,Null_input,Null_source)

				mpih%sendbuf(i,:) = (/J,grad,adj_grad(1)/)

				call destroyRecord(r)
				call destroyPM1D(d)
				call destroyAdjoint(adj)
				call destroyRecord(fsr)
				call destroyPM1D(dpm)
			end do

			call gatherData(mpih)

			if( mpih%my_rank.eq.mpih%size-1 ) then
				write(Time_str,'(F5.2)') Time(k)
				open(unit=301,file='data/Debye_sampling/Jk_'//trim(Time_str)//'.bin',	&
						status='replace',form='unformatted',access='stream')
				open(unit=302,file='data/Debye_sampling/gradk_'//trim(Time_str)//'.bin',	&
						status='replace',form='unformatted',access='stream')
				open(unit=303,file='data/Debye_sampling/adjk_'//trim(Time_str)//'.bin',	&
						status='replace',form='unformatted',access='stream')
			   write(301) mpih%recvbuf(:,1)
			   write(302) mpih%recvbuf(:,2)
				write(303) mpih%recvbuf(:,3)
			   close(301)
			   close(302)
				close(303)
			end if
		end do

		call destroyMPIHandler(mpih)
	end subroutine

end program
