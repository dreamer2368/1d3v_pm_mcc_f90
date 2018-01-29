module testmodule

	use init
	use timeStepAdj
	use timeStepFSens
	use modMPI

	implicit none

contains

	subroutine redistribute_temp_test
		type(PM1D) :: pm
		type(FSens) :: fs
		real(mp), parameter :: L=2.7_mp, Lv=4.6_mp, wv=0.4_mp, dt=0.05_mp
		integer, parameter :: Ng=64, NInject=5E3, NLimit=5E4

		call pm%buildPM1D(1.0_mp,1.0_mp,Ng,N=1,pBC=0,mBC=0,order=1,L=L,dt=dt)
		call pm%p(1)%buildSpecies(1.0_mp,1.0_mp)
		
		allocate(pm%p(1)%xp(NLimit))
		allocate(pm%p(1)%vp(NLimit,3))
		allocate(pm%p(1)%spwt(NLimit))
		call RANDOM_NUMBER(pm%p(1)%xp)
		pm%p(1)%xp = pm%p(1)%xp*L
		pm%p(1)%vp = wv*randn(NLimit,3)
		pm%p(1)%spwt = 1.0_mp/NLimit

		call fs%buildFSens(pm,Lv,Ng/2,NInject,NLimit)
		call fs%p(1)%buildSpecies(1.0_mp,1.0_mp)
		call fs%p(1)%setSpecies(NLimit,pm%p(1)%xp,pm%p(1)%vp,pm%p(1)%spwt)
		call Redistribute_temp(fs,fs%p(1))

		open(unit=301,file='data/redistribution/xp0.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/redistribution/vp0.bin',status='replace',form='unformatted',access='stream')
		write(301) pm%p(1)%xp
		write(302) pm%p(1)%vp(:,1)
		close(301)
		close(302)

		open(unit=301,file='data/redistribution/xp1.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/redistribution/vp1.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/redistribution/spwt1.bin',status='replace',form='unformatted',access='stream')
		write(301) fs%p(1)%xp
		write(302) fs%p(1)%vp(:,1)
		write(303) fs%p(1)%spwt
		close(301)
		close(302)
		close(303)

		call destroyPM1D(pm)
		call destroyFSens(fs)
	end subroutine

	subroutine debye_adj(fk,Time,str,k,output)
		real(mp), intent(in) :: fk, Time
		character(len=*), intent(in) ::str
		integer, intent(in) :: k
		real(mp), intent(out) :: output(:)
		type(PM1D) :: d
		type(adjoint) :: adj
		type(recordData) :: r
		real(mp) :: vT = 1.5_mp
		integer :: N = 1E5, Ng=64
		real(mp) :: L = 20.0_mp, Wp, Q = 2.0_mp
		real(mp) :: dt = 0.05_mp
		real(mp) :: A(2),J0,J1,grad(1)

		SELECT CASE(k)
			CASE(0)
				A = (/ vT, 0.0_mp /)
				call buildPM1D(d,Time,0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=dt)
				call buildRecord(r,d%nt,1,d%L,d%ng,trim(str)//'/before',20)

				call buildSpecies(d%p(1),-1.0_mp,1.0_mp)
				call Debye_initialize(d,N,Q)

				call forwardsweep(d,r,Null_input,Null_source,Debye,J0)

				call printPlasma(r)
				print *, 'J0=',J0

				call buildAdjoint(adj,d)
				call adj%m%setMesh(d%m%rho_back)
				call backward_sweep(adj,d,r,grad,dDebye,Null_dinput,dDebye_dvT,Null_input,Null_source)

				print *, 'dJdA=',grad

				call destroyPM1D(d)
				call destroyRecord(r)
				call destroyAdjoint(adj)

				output = (/ J0, grad(1) /)
			CASE(1)
				A = (/ vT+fk, 0.0_mp /)
				call buildPM1D(d,Time,0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=dt)
				call buildRecord(r,d%nt,1,d%L,d%ng,trim(str)//'/after',20)

				call buildSpecies(d%p(1),-1.0_mp,1.0_mp)
				call Debye_initialize(d,N,Q)

				call forwardsweep(d,r,Null_input,Null_source,Debye,J1)
				print *, 'J1=',J1

				call destroyRecord(r)
				call destroyPM1D(d)
				
				output(1) = J1
		END SELECT
	end subroutine

	subroutine updateWeightTest
		type(PM1D) :: pm
		type(FSens) :: fs
		type(recordData) :: r
		real(mp), parameter :: Tf=1.0_mp, Ti=0.5_mp, dt = 1.0_mp
		real(mp), parameter :: L = 1.0_mp, Lv=0.5_mp, w = 0.1_mp
		integer, parameter :: Ng=64, N=2E6, NInject=1E6
		integer :: i,k
		integer :: g(2)
		real(mp) :: frac(2)
		real(mp) :: xp0(N), vp0(N,3), spwt0(N)
		real(mp), dimension(Ng,Ng+1) :: j

		call pm%buildPM1D(Tf,Ti,Ng,N=1,pBC=0,mBC=0,order=1,L=L,dt=dt)
		call pm%p(1)%buildSpecies(1.0_mp,1.0_mp)
		call fs%buildFSens(pm,Lv,Ng/2,NInject,NInject)
		call fs%p(1)%setSpecies(1,(/0.0_mp/),(/0.0_mp,0.0_mp,0.0_mp/),(/0.0_mp/))
		call r%buildRecord(pm%nt,N,pm%L,Ng,'updateWeightTest',1)

		call RANDOM_NUMBER(xp0)
		xp0 = xp0*L
		call RANDOM_NUMBER(vp0)
		vp0 = (2.0_mp*vp0-1.0_mp)*Lv
		xp0(1:N/2) = randn(N/2)*w + 0.5_mp*L
		vp0(1:N/2,:) = randn(N/2,3)*w
		spwt0 = 1.0_mp/N
		call fs%p(1)%setSpecies(N,xp0,vp0,spwt0)

		call fs%applyBC(fs%p(1),fs%m,fs%dt,fs%A0(1))
		!X-direction Interpolation
		do k=1,N
			CALL fs%a(1)%assignMatrix(xp0(k),pm%m%dx,g,frac)
			CALL fs%a(1)%adjustGrid(fs%m%ng,g,frac)
			fs%a(1)%g(:,k) = g
			fs%a(1)%frac(:,k) = frac
		end do

		call fs%FSensDistribution(fs%p(1),fs%a(1))

		open(unit=300,file='data/updateWeightTest/record.bin',status='replace',form='unformatted',access='stream')
		open(unit=301,file='data/updateWeightTest/xp.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/updateWeightTest/vp.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/updateWeightTest/spwt.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/updateWeightTest/f_A_before.bin',status='replace',form='unformatted',access='stream')
		write(300) Ng, Ng/2, N
		write(301) xp0
		write(302) vp0
		write(303) spwt0
		write(304) fs%f_A
		close(300)
		close(301)
		close(302)
		close(303)
		close(304)

		j = -1.0_mp

		call fs%updateWeight(fs%p(1),fs%a(1),fs%f_A,j)
		open(unit=304,file='data/updateWeightTest/n_A.bin',status='replace',form='unformatted',access='stream')
		write(304) fs%f_A
		close(304)
		call fs%FSensDistribution(fs%p(1),fs%a(1))

		open(unit=304,file='data/updateWeightTest/f_A_after.bin',status='replace',form='unformatted',access='stream')
		write(304) fs%f_A
		close(304)

		call pm%destroyPM1D
		call fs%destroyFSens
		call r%destroyRecord
	end subroutine

	subroutine RedistributionTest
		type(PM1D) :: pm
		type(FSens) :: fs
		type(recordData) :: r
		real(mp), parameter :: Tf=1.0_mp, Ti=0.5_mp, dt = 1.0_mp
		real(mp), parameter :: L = 1.0_mp, Lv=0.5_mp, w = 0.1_mp
		integer, parameter :: Ng=64, N=2E6, NInject=1E6
		integer :: i,k
		integer :: g(2)
		real(mp) :: frac(2)
		real(mp) :: xp0(N), vp0(N,3), spwt0(N)

		call pm%buildPM1D(Tf,Ti,Ng,N=1,pBC=0,mBC=0,order=1,L=L,dt=dt)
		call pm%p(1)%buildSpecies(1.0_mp,1.0_mp)
		call fs%buildFSens(pm,Lv,Ng/2,NInject,NInject)
		call fs%p(1)%setSpecies(1,(/0.0_mp/),(/0.0_mp,0.0_mp,0.0_mp/),(/0.0_mp/))
		call r%buildRecord(pm%nt,N,pm%L,Ng,'RedistributionTest',1)

		call RANDOM_NUMBER(xp0)
		xp0 = xp0*pm%L
		vp0 = randn(N,3)
		vp0 = vp0*w
		spwt0 = pm%L*SIN(2.0_mp*pi*xp0/pm%L)/N
		call fs%p(1)%setSpecies(N,xp0,vp0,spwt0)

		call fs%applyBC(fs%p(1),fs%m,fs%dt,fs%A0(1))
		!X-direction Interpolation
		do k=1,N
			CALL fs%a(1)%assignMatrix(xp0(k),pm%m%dx,g,frac)
			CALL fs%a(1)%adjustGrid(fs%m%ng,g,frac)
			fs%a(1)%g(:,k) = g
			fs%a(1)%frac(:,k) = frac
		end do

		call fs%FSensDistribution(fs%p(1),fs%a(1))

		open(unit=300,file='data/RedistributionTest/record.bin',status='replace',form='unformatted',access='stream')
		open(unit=301,file='data/RedistributionTest/xp.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/RedistributionTest/vp.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/RedistributionTest/spwt.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/RedistributionTest/f_A.bin',status='replace',form='unformatted',access='stream')
		write(300) Ng, Ng/2, N
		write(301) xp0
		write(302) vp0
		write(303) spwt0
		write(304) fs%f_A
		close(300)
		close(301)
		close(302)
		close(303)
		close(304)

		call fs%Redistribute(fs%p(1),fs%a(1))
		call fs%FSensDistribution(fs%p(1),fs%a(1))

		open(unit=301,file='data/RedistributionTest/xp_rdst.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/RedistributionTest/vp_rdst.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/RedistributionTest/spwt_rdst.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/RedistributionTest/f_A_rdst.bin',status='replace',form='unformatted',access='stream')
		write(301) fs%p(1)%xp
		write(302) fs%p(1)%vp
		write(303) fs%p(1)%spwt
		write(304) fs%f_A
		close(301)
		close(302)
		close(303)
		close(304)

		call pm%destroyPM1D
		call fs%destroyFSens
		call r%destroyRecord
	end subroutine

	subroutine SensitivityInitializeTest
		type(PM1D) :: pm
		type(FSens) :: fs
		type(recordData) :: r
		real(mp), parameter :: Tf=1.0_mp, Ti=0.5_mp, vT=1.0_mp, dt = 1.0_mp
		real(mp), parameter :: L = 20.0_mp, Lv=5.0_mp, w = 1.0_mp
		integer, parameter :: Ng=256, N=1000000, NInit=1000000
		integer :: i,k
		integer :: g(2)
		real(mp) :: frac(2)
		real(mp) :: xp0(N), vp0(N,3), spwt0(N)

		call pm%buildPM1D(Tf,Ti,Ng,N=1,pBC=0,mBC=0,order=1,L=L,dt=dt,A=(/vT/))
		call pm%p(1)%buildSpecies(1.0_mp,1.0_mp)
		call fs%buildFSens(pm,Lv,Ng/2,NInit,NInit)
		call r%buildRecord(pm%nt,N,pm%L,Ng,'SensitivityInitTest',1)

		open(unit=300,file='data/SensitivityInitTest/record.bin',status='replace',form='unformatted',access='stream')
		write(300) Ng, Ng/2, N
		close(300)

		call Debye_sensitivity_init(fs,NInit,fs%A0(1))
		call fs%applyBC(fs%p(1),fs%m,fs%dt,fs%A0(1))
		!X-direction Interpolation
		do k=1,N
			CALL fs%a(1)%assignMatrix(fs%p(1)%xp(k),pm%m%dx,g,frac)
			CALL fs%a(1)%adjustGrid(fs%m%ng,g,frac)
			fs%a(1)%g(:,k) = g
			fs%a(1)%frac(:,k) = frac
		end do
		fs%m%E = 1.0_mp
		call fs%FSensDistribution(pm%p(1),pm%a(1))
		call fs%FSensSourceTerm(pm%p(1)%qs,pm%p(1)%ms,fs%f_A,pm%m%E)

		open(unit=301,file='data/SensitivityInitTest/xp.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/SensitivityInitTest/vp.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/SensitivityInitTest/spwt.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/SensitivityInitTest/j.bin',status='replace',form='unformatted',access='stream')
		write(301) fs%p(1)%xp
		write(302) fs%p(1)%vp
		write(303) fs%p(1)%spwt
		write(304) fs%j
		close(301)
		close(302)
		close(303)
		close(304)

		call pm%destroyPM1D
		call fs%destroyFSens
		call r%destroyRecord
	end subroutine

	subroutine MPITest
		type(mpiHandler) :: mpih
		integer, parameter :: Nsample=10000, Ndata=3
      integer :: i

		call mpih%buildMPIHandler
		call mpih%allocateBuffer(Nsample,Ndata)

		do i=1,mpih%sendcnt
			mpih%sendbuf(i,:) = (/i*0.2_mp,i*0.7_mp,i*1.7_mp/)
		end do

		call mpih%gatherData

		call mpih%destroyMPIHandler
	end subroutine

	subroutine InjectionTest
		type(PM1D) :: pm
		type(FSens) :: fs
		type(recordData) :: r
		real(mp), parameter :: Tf=1.0_mp, Ti=0.5_mp, dt = 1.0_mp
		real(mp), parameter :: L = 1.0_mp, Lv=0.5_mp, w = 0.1_mp
		integer, parameter :: Ng=64, N=1000000, NInject=1000000
		integer :: i,k
		integer :: g(2)
		real(mp) :: frac(2)
		real(mp) :: xp0(N), vp0(N,3), spwt0(N)

		call pm%buildPM1D(Tf,Ti,Ng,N=1,pBC=0,mBC=0,order=1,L=L,dt=dt)
		call pm%p(1)%buildSpecies(1.0_mp,1.0_mp)
		call fs%buildFSens(pm,Lv,Ng/2,NInject,NInject)
		call fs%p(1)%setSpecies(1,(/0.0_mp/),(/0.0_mp,0.0_mp,0.0_mp/),(/0.0_mp/))
		call r%buildRecord(pm%nt,N,pm%L,Ng,'InjectionTest',1)

		call RANDOM_NUMBER(xp0)
		xp0 = xp0*pm%L
		vp0 = randn(N,3)
		vp0 = vp0*w
		spwt0 = 1.0_mp/N
!		spwt0 = SQRT(2.0_mp*pi)*w/EXP( -vp0(:,1)**2/2.0_mp/w/w )/N
		call pm%p(1)%setSpecies(N,xp0,vp0,spwt0)
!		fs%m%E = 1.0_mp
		fs%m%E = (/ (SIN( 2.0_mp*pi*i/Ng ),i=1,Ng) /)

		call pm%applyBC(pm%p(1),pm%m,pm%dt,pm%A0(1))
		!X-direction Interpolation
		do k=1,N
			CALL pm%a(1)%assignMatrix(pm%p(1)%xp(k),pm%m%dx,g,frac)
			CALL pm%a(1)%adjustGrid(pm%m%ng,g,frac)
			pm%a(1)%g(:,k) = g
			pm%a(1)%frac(:,k) = frac
		end do

		call fs%FSensDistribution(pm%p(1),pm%a(1))
		call fs%FSensSourceTerm(pm%p(1)%qs,pm%p(1)%ms,fs%f_A,pm%m%E)

		open(unit=300,file='data/InjectionTest/record.bin',status='replace',form='unformatted',access='stream')
		open(unit=301,file='data/InjectionTest/xp.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/InjectionTest/vp.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/InjectionTest/spwt.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/InjectionTest/j.bin',status='replace',form='unformatted',access='stream')
		write(300) Ng, Ng/2, N
		write(301) xp0
		write(302) vp0
		write(303) spwt0
		write(304) fs%j
		close(300)
		close(301)
		close(302)
		close(303)
		close(304)

		do k=1,Ng+1
			do i=1,Ng
				fs%j(i,k) = SIN( 2.0_mp*pi*i/Ng )/SQRT(2.0_mp*pi)/w*EXP( -((k-Ng/2-1)*fs%dv)**2/2.0_mp/w/w )
			end do
		end do
		open(unit=304,file='data/InjectionTest/j_source.bin',status='replace',form='unformatted',access='stream')
		write(304) fs%j
		close(304)
		call fs%InjectSource(fs%p(1),fs%j)

		call fs%applyBC(fs%p(1),fs%m,fs%dt,fs%A0(1))
		!X-direction Interpolation
		do k=1,N
			CALL fs%a(1)%assignMatrix(fs%p(1)%xp(k),pm%m%dx,g,frac)
			CALL fs%a(1)%adjustGrid(fs%m%ng,g,frac)
			fs%a(1)%g(:,k) = g
			fs%a(1)%frac(:,k) = frac
		end do
		call fs%FSensDistribution(fs%p(1),fs%a(1))

		open(unit=301,file='data/InjectionTest/xp_inject.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/InjectionTest/vp_inject.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/InjectionTest/spwt_inject.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/InjectionTest/j_inject.bin',status='replace',form='unformatted',access='stream')
		write(301) fs%p(1)%xp
		write(302) fs%p(1)%vp
		write(303) fs%p(1)%spwt
		write(304) fs%f_A
		close(301)
		close(302)
		close(303)
		close(304)

		call pm%destroyPM1D
		call fs%destroyFSens
		call r%destroyRecord
	end subroutine

	subroutine random_test
		real(mp) :: test(5)
		integer :: ierr, my_rank, s

		call MPI_INIT(ierr)
		call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

		call init_random_seed(my_rank)
		call RANDOM_NUMBER(test)
      print *, 'My rank: ', my_rank
      print *, test

		call MPI_FINALIZE(ierr)
	end subroutine

	subroutine twostream_grad(fk,Ti,str,k,output)
		real(mp), intent(in) :: fk, Ti
		character(len=*), intent(in) ::str
		integer, intent(in) :: k
		real(mp), intent(out) :: output(:)
		type(adjoint) :: adj
		type(PM1D) :: pm
		type(recordData) :: r
		real(mp) :: Tf
		integer, parameter :: Ng=64, Np=100000, N=1
		real(mp) :: v0 = 0.2_mp, vT = 0.01_mp, dt=0.1_mp
		integer :: mode=1
		character(len=100)::dir
		real(mp) :: J0,J1,grad(1)
		output = 0.0_mp
		Tf = 0.1_mp + Ti

		SELECT CASE(k)
			CASE(0)
				call buildPM1D(pm,Tf,Ti,Ng,N,0,0,1,dt=dt,A=(/1.0_mp,0.0_mp/))
				dir = str//'/before'
				call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir),10)
				call set_null_discharge(r)
				call init_random_seed
				call twostream_initialize(pm,Np,v0,vT,mode)
				call forwardsweep(pm,r,Te,Null_source,MPE,J0)
!				call printPlasma(r)
				print *, 'J0=',J0

				call buildAdjoint(adj,pm)
				call backward_sweep(adj,pm,r,grad,dMPE,dTe,dTedA,Te,Null_source)

				print *, 'dJdA=',grad

				output(1:2) = (/J0, grad(1)/)

				call destroyAdjoint(adj)
				call destroyPM1D(pm)
				call destroyRecord(r)
			CASE(1)
				dir = str//'/after'
				call buildPM1D(pm,Tf,Ti,Ng,N,0,0,1,dt=dt,A=(/1.0_mp,fk/))
				call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir),10)
				call init_random_seed
				call twostream_initialize(pm,Np,v0,vT,mode)
				call forwardsweep(pm,r,Te,Null_source,MPE,J1)
!				call printPlasma(r)
				print *, 'J1=',J1
				output(1) = J1

				call destroyRecord(r)
				call destroyPM1D(pm)
		END SELECT
	end subroutine

	subroutine Landau(fk,Ti,str,k,output)
		real(mp), intent(in) :: fk, Ti
		character(len=*), intent(in) ::str
		integer, intent(in) :: k
		real(mp), intent(out) :: output(:)
		type(adjoint) :: adj
		type(PM1D) :: pm
		type(recordData) :: r
		real(mp) :: Tf
		integer, parameter :: Ng=64, Np=3*10**5, N=1
		real(mp) :: vT = 1.0_mp, L=4.0_mp*pi
		real(mp) :: dt=0.1_mp
		character(len=100)::dir
		real(mp) :: J0,J1,grad(1)
		output = 0.0_mp
		Tf = 0.1_mp + Ti

		SELECT CASE(k)
			CASE(0)
				call buildPM1D(pm,Tf,Ti,Ng,N,0,0,1,dt=dt,L=L,A=(/0.1_mp,0.0_mp/))
				dir = str//'/before'
				call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir),10)
				call set_null_discharge(r)
				call init_random_seed
				call Landau_initialize(pm,Np,vT)
				call forwardsweep(pm,r,Te,Null_source,MPE,J0)
!				call printPlasma(r)
				print *, 'J0=',J0

				call buildAdjoint(adj,pm)
				call backward_sweep(adj,pm,r,grad,dMPE,dTe,dTedA,Te,Null_source)

				print *, 'dJdA=',grad

				output(1:2) = (/J0, grad(1)/)

				call destroyAdjoint(adj)
				call destroyPM1D(pm)
				call destroyRecord(r)
			CASE(1)
				dir = str//'/after'
				call buildPM1D(pm,Tf,Ti,Ng,N,0,0,1,dt=dt,L=L,A=(/0.1_mp,fk/))
				call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir),10)
				call init_random_seed
				call Landau_initialize(pm,Np,vT)
				call forwardsweep(pm,r,Te,Null_source,MPE,J1)
!				call printPlasma(r)
				print *, 'J1=',J1
				output(1) = J1

				call destroyRecord(r)
				call destroyPM1D(pm)
		END SELECT
	end subroutine

	subroutine twostream_adj(fk,ek)
		real(mp), intent(out) :: ek
		real(mp), intent(in) :: fk
		type(PM1D) :: pm
		type(adjoint) :: adj
		type(recordData) :: r
		integer, parameter :: Ng=64, Np=10**4, N=1
		real(mp) :: v0 = 0.2_mp, vT = 0.0_mp
		integer :: mode=1
		real(mp) :: J0,J1, grad(1)
		character(len=100)::dir1
		ek = 0.0_mp

		call buildPM1D(pm,45.0_mp,40.0_mp,Ng,N,0,0,1,A=(/1.0_mp,0.0_mp/))
		dir1='twostream/before'
		call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir1),5)
		call set_null_discharge(r)
		call twostream_initialize(pm,Np,v0,vT,mode)
		call forwardsweep(pm,r,IC_wave,Null_source,MKE,J0)
		call printPlasma(r)
		print *, 'J0=',J0

		call buildAdjoint(adj,pm)
		call backward_sweep(adj,pm,r,grad,dMKE,dIC_wave,dIC_wave_dB,IC_wave,Null_source)

		print *, 'dJdA=',grad

		call destroyPM1D(pm)
		call destroyRecord(r)
		dir1 = 'twostream/after'
		call buildPM1D(pm,45.0_mp,40.0_mp,Ng,N,0,0,1,A=(/1.0_mp,fk/))
		call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir1),5)
		call set_null_discharge(r)
		call twostream_initialize(pm,Np,v0,vT,mode)
		call forwardsweep(pm,r,IC_wave,Null_source,MKE,J1)
		print *, 'J1=',J1
		print *, 'dJdA=',(J1-J0)/fk

		ek = ABS( (grad(1) - (J1-J0)/fk)/grad(1) )
		print *, 'ek=',ek

		call destroyAdjoint(adj)
		call destroyRecord(r)
		call destroyPM1D(pm)
	end subroutine

	subroutine test_backward_sweep
		type(PM1D) :: pm
		type(adjoint) :: adj
		type(recordData) :: r
		integer, parameter :: Ng=64, Np=2, N=2
		real(mp) :: qs = 3.4_mp, ms = 2.7_mp, spwt = 1.9_mp, xp(Np), vp(Np,3)
		real(mp) :: rho_back(Ng), J0=0.0_mp,J1=0.0_mp,grad(1)
		real(mp) :: fxp = (0.1_mp)**9, dxp
		integer :: j
		character(len=1000)::dir1,dir2

		call buildPM1D(pm,2.0_mp,1.0_mp,Ng,N,0,0,1,eps=2.3_mp)
		dir1='test_backward_sweep/before'
		call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir1),3)
		call set_null_discharge(r)
		call buildSpecies(pm%p(1),qs,ms)
		call buildSpecies(pm%p(2),-qs,ms)

		print *, 'Original xp'
		xp = 0.4_mp*pm%L*(/ (j,j=1,Np) /)
		vp = 0.1_mp
		print *, xp
		print *, vp
		call setSpecies(pm%p(1),Np,xp,vp,(/spwt,spwt/))
		call setSpecies(pm%p(2),Np,pm%L-xp,0.2_mp*vp,(/spwt,spwt/))
		rho_back = 0.0_mp
		call setMesh(pm%m,rho_back)

		call forwardsweep(pm,r,Null_input,Null_source,TestQoI,J0)
		call printPlasma(r)

		call buildAdjoint(adj,pm)
		call backward_sweep(adj,pm,r,grad,dTestQoI,Null_Dinput,Null_dJdA,Null_input,Null_source)

		print *, 'dJdxp1:', -adj%p(1)%xp/pm%dt
		print *, 'dJdxp2:', -adj%p(2)%xp/pm%dt
		print *, 'dJdvp1:', -adj%p(1)%vp(:,1)/pm%dt
		print *, 'dJdvp2:', -adj%p(2)%vp(:,1)/pm%dt

		call destroySpecies(pm%p(1))
		call destroySpecies(pm%p(2))
		call setSpecies(pm%p(2),Np,pm%L-xp,0.2_mp*vp,(/spwt,spwt/))
		dxp = xp(2)*fxp
		xp(2) = xp(2) + dxp
		call setSpecies(pm%p(1),Np,xp,vp,(/spwt,spwt/))
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
		real(mp) :: qs = 3.4_mp, ms = 2.7_mp, xp(Np), vp(Np,3), spwt(Np)
		real(mp) :: rho_back(Ng)
		real(mp) :: weight(Ng), J0, J1
		real(mp) :: rhs(Ng), rho1(Ng-1), dxp1(Np), dxp2(Np)
		real(mp) :: fxp = (0.1_mp)**9, dxp

		call buildPM1D(pm,40.0_mp,20.0_mp,Ng,1,0,0,1,eps=2.3_mp)
		call buildSpecies(pm%p(1),qs,ms)
		!particle, mesh setup
		print *, 'Original xp'
		xp = 0.4_mp*pm%L*(/ (j,j=1,Np) /)
		vp = 0.1_mp
		spwt = 1.9_mp
		print *, xp
		print *, vp
		call setSpecies(pm%p(1),Np,xp,vp,spwt)
		rho_back = -qs*Np/pm%L
		call setMesh(pm%m,rho_back)
		!one time-step
		call moveSpecies(pm%p(1),pm%dt)
		call pm%applyBC(pm%p(1),pm%m,pm%dt,pm%A0(1))
		pm%m%rho = 0.0_mp
		call pm%a(1)%chargeAssign(pm%p(i),pm%m)
		call solveMesh(pm%m,pm%eps0)
		pm%m%E = - multiplyD(pm%m%phi,pm%m%dx,pm%m%BCindex)
		call pm%a(1)%forceAssign(pm%p(1), pm%m)
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
		print *, 'Perturbed vp'
		print *, vp
		call destroySpecies(pm%p(1))
		call setSpecies(pm%p(1),Np,xp,vp,spwt)

		call moveSpecies(pm%p(1),pm%dt)
		call pm%applyBC(pm%p(1),pm%m,pm%dt,pm%A0(1))
		pm%m%rho = 0.0_mp
		call pm%a(1)%chargeAssign(pm%p(i),pm%m)
		call solveMesh(pm%m,pm%eps0)
		pm%m%E = - multiplyD(pm%m%phi,pm%m%dx,pm%m%BCindex)
		call forceAssign(pm%a(1), pm%p(1), pm%m)
		call accelSpecies(pm%p(1),pm%dt)
		J1 = sum( pm%p(1)%vp(:,1)**2 )
		!J1 = pm%m%dx*sum( weight*pm%m%E**2 )
		print *, 'J1 = ',J1
		print *, 'dJdvp(',k,')=', (J1-J0)/dxp

		call destroyPM1D(pm)
		call destroyAdjoint(adj)
	end subroutine

	subroutine test_refluxing_boundary
		type(PM1D) :: reflux
		type(recordData) :: r
		integer, parameter :: Ng=64, N=10000, order=1
		real(mp) :: Ti=20, Tf = 40
		real(mp) :: xp0(N), vp0(N,3), spwt0(N) = 1.0_mp, rho_back(Ng), qe, me
		integer :: i

		call buildPM1D(reflux,Tf,Ti,Ng,1,pBC=2,mBC=2,order=order,A=(/ 1.0_mp, 1.0_mp /))
		call buildRecord(r,reflux%nt,1,reflux%L,Ng,'test_reflux',1)

		xp0 = -0.5_mp*reflux%L
		vp0 = 0.0_mp
		rho_back = 0.0_mp
		qe = -(0.1_mp)**2/(N/reflux%L)
		me = -qe
		rho_back(Ng) = -qe
		call buildSpecies(reflux%p(1),qe,me)
		call setSpecies(reflux%p(1),N,xp0,vp0,spwt0)
		call setMesh(reflux%m,rho_back)

		call reflux%applyBC(reflux%p(1),reflux%m,reflux%dt,reflux%A0(1))
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

	subroutine forYeoh
		type(PM1D) :: pm
		type(recordData) :: r
		real(mp), parameter :: n0=10.0_mp**15
		real(mp), parameter :: L = 0.02_mp, area = 0.016_mp               !L(m), area(m2)
		integer, parameter :: nc2p = 10**7, Np=CEILING(n0*L*area/nc2p), Ng = 300
		real(mp), parameter :: T=300.0_mp
		real(mp), parameter :: me=9.11E-31, mi=100.0_mp*me
		real(mp), parameter :: qe=-1.60E-19, qi=-qe
		real(mp), parameter :: spwt=n0*L/Np					!spwt=nc2p/area
      real(mp), parameter :: Kb = 1.38065E-23, eps = 8.85418782E-12
		real(mp) :: wp0, lambda0, v0_e, v0_i
		real(mp) :: spwt0(Np), xp0(Np), vp0(Np,3), rho_back(Ng)
		real(mp) :: dt
		integer :: i

		print *, 'n0(m-3): ',n0,', spwt(m-2): ',spwt

		v0_e = sqrt(2.0_mp*Kb*T/me)
		v0_i = sqrt(2.0_mp*Kb*T/mi)
		wp0 = sqrt(n0*q_e*q_e/m_e/eps)
		lambda0 = sqrt(eps*Kb*T/n0/q_e/q_e)

		print *, 'L = ',L,', lambda0 = ',lambda0,' e = lambda/L = ',lambda0/L

		dt = 7.20179E-11
		print *, 'dt = ',dt,', wp*dt=',wp0*dt

		call buildPM1D(pm,400.0_mp*dt,100.0_mp*dt,Ng,2,pBC=0,mBC=0,order=1,dt=dt,L=L,eps=eps,A=(/Kb*T/))
		call buildRecord(r,pm%nt,2,pm%L,pm%ng,'forYeoh',1)

		call RANDOM_NUMBER(xp0)
		xp0 = xp0*L
		vp0 = ABS(randn(Np,3))*v0_e
		spwt0 = spwt
		call pm%p(1)%buildSpecies(qe,me)
		call pm%p(1)%setSpecies(Np,xp0,vp0,spwt0)

		call RANDOM_NUMBER(xp0)
		xp0 = xp0*L
		vp0 = ABS(randn(Np,3))*v0_i
		spwt0 = spwt
		call pm%p(2)%buildSpecies(qi,mi)
		call pm%p(2)%setSpecies(Np,xp0,vp0,spwt0)

		rho_back = 0.0_mp
		call pm%m%setMesh(rho_back)

		call forwardsweep(pm,r,Null_input,Null_source)
		call printPlasma(r)

		call destroyPM1D(pm)
		call destroyRecord(r)
	end subroutine

    subroutine MPI_write_test
		type(mpiHandler) :: mpih
		integer, parameter  :: Nsample=50
		integer :: i,k
		character(len=100) :: prefix, dir_temp
        real(mp), dimension(Nsample) :: J, dJf, dJadj
        prefix = 'MPIwrite'

		call buildMPIHandler(mpih)
        call allocateBuffer(Nsample,3,mpih)

		call init_random_seed(mpih%my_rank)
        do i=1,mpih%sendcnt
            call RANDOM_NUMBER(mpih%writebuf)
            call writeData(mpih,trim(prefix))
        end do

        if( mpih%my_rank.eq.mpih%size-1 ) then
			dir_temp=trim(prefix)//'_J.bin' 
            open(unit=305,file=trim(dir_temp),form='unformatted',access='stream')
			dir_temp=trim(prefix)//'_dJf.bin' 
			open(unit=306,file=trim(dir_temp),form='unformatted',access='stream')
			dir_temp=trim(prefix)//'_dJadj.bin' 
			open(unit=307,file=trim(dir_temp),form='unformatted',access='stream')
			read(305) J
			read(306) dJf
			read(307) dJadj
			close(305)
			close(306)
			close(307)

            print *, 'J'
            print *, J
            print *, 'dJf'
            print *, dJf
            print *, 'dJadj'
            print *, dJadj
        end if

        call destroyMPIHandler(mpih)
    end subroutine

end module
