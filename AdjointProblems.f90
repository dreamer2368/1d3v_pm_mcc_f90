module AdjointProblems

	use init
	use timeStepAdj
	use modMPI

	implicit none

contains

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

end module
