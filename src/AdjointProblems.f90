module AdjointProblems

	use init
	use timeStepAdj
	use modMPI
    use modInputHelper

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

		call init_random_seed(my_rank,addSystemTime=.true.)

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
		type(mpiHandler) :: mpih
		type(adjoint) :: adj
		type(PM1D) :: pm
		type(recordData) :: r
		integer, parameter :: Nsample = 1E4, Nt = 5
        real(mp) :: Time(Nt)
		integer, parameter :: Ng=64, Np=3*10**5, N=1
		real(mp) :: vT = 1.0_mp, L=4.0_mp*pi
		real(mp) :: dt=0.1_mp
		real(mp) :: J0,grad(1)
		integer :: i,j, thefile
		character(len=100) :: rank_str,prefix,dir,Time_str

        Time = (/ (240.0_mp**(1.0_mp*(i-1)/(Nt-1)),i=1,Nt) /)

		call mpih%buildMPIHandler
		call mpih%allocateBuffer(Nsample,2)

		call init_random_seed(mpih%my_rank,addSystemTime=.true.)

!        do j=1,1
         j=5
            write(Time_str,'(I03.3)') INT(Time(j))
            prefix = 'Landau_sampling/T'//trim(Time_str)
            dir = 'data/'//trim(prefix)
            thefile = mpih%MPIWriteSetup(dir)
            do i=1,mpih%sendcnt
			    call buildPM1D(pm,Time(j)+0.1_mp,Time(j),Ng,N,0,0,1,dt=dt,L=L,A=(/0.1_mp,0.0_mp/))
    		    dir = trim(prefix)//'/'//trim(adjustl(mpih%rank_str))
    			call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir),10)
    			call set_null_discharge(r)
    			call Landau_initialize(pm,Np,vT)
    
    			call forwardsweep(pm,r,Te,Null_source,MPE,J0)
    
    			call buildAdjoint(adj,pm)
    			call backward_sweep(adj,pm,r,grad,dMPE,dTe,dTedA,Te,Null_source)

                mpih%writebuf = (/ J0,grad(1) /)
                call MPI_FILE_WRITE(thefile, mpih%writebuf, 2, MPI_DOUBLE, & 
                                    MPI_STATUS_IGNORE, mpih%ierr)
                call MPI_FILE_SYNC(thefile,mpih%ierr)
    
    			call destroyAdjoint(adj)
    			call destroyRecord(r)
                print ('(A,I5,A,F8.3,A,I5,A)'), 'Rank-',mpih%my_rank,    &
                                               'Time: ',Time(j),        &
                                               ', Sample-',i,' is collected.'
            end do
            call MPI_FILE_CLOSE(thefile, mpih%ierr)            
!        end do

        call mpih%destroyMPIHandler
	end subroutine

	subroutine adjoint_convergence_in_time(problem)
		type(mpiHandler) :: mpih
		integer, parameter :: N=70, Nt=9
		real(mp) :: fk, Tk(Nt), ek
		character(len=100) :: dir, Tstr, filename
		integer :: i,j
		real(mp) :: J0,J1,grad,temp(2)
        integer :: thefile
        integer(kind=MPI_OFFSET_KIND) :: disp
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

		call mpih%buildMPIHandler
        if( mpih%my_rank.eq.0 .and.                                          &
            mpih%size.ne.(N+1) ) then
            print *, 'Required same number of processors!'
            print *, 'finite difference iteration: ', N+1
            print *, 'simulation times: ', Nt
            print *, 'required processors s=N*Nt: ',(N+1)

            call mpih%destroyMPIHandler
            stop
        end if

		call mpih%allocateBuffer(N+1,2)

        fk = 1.0_mp*MOD(mpih%my_rank,N+1)
        fk = 10.0_mp*EXP( -0.75_mp*(fk-1) )
        Tk = (/ (60.0_mp*(i-1)/(Nt-1),i=1,Nt) /)
        Tk(1) = 0.03_mp
        Tk = Tk*2.0_mp*pi
		ek = 0.0_mp

		dir = 'Landau_adj_test_dp'

        do j=1,Nt
            if( MOD(mpih%my_rank,N+1).eq.0 ) then
                call problem(fk,Tk(j),trim(dir),0,temp)
                mpih%sendbuf(1,:) = temp
            else
                call problem(fk,Tk(j),trim(dir),1,temp)
                mpih%sendbuf(1,:) = (/ fk, temp(1) /)
            end if

			call mpih%gatherData

			if( mpih%my_rank.eq.mpih%size-1 ) then
                write(Tstr,'(I01)') j
                filename = 'data/'//trim(dir)//'/grad_convergence.'//trim(Tstr)//'.dat'
				open(unit=301, file=filename, status='replace')
                do i=1,N+1
			        write(301,*) mpih%recvbuf(i,:)
                end do
			    close(301)
                print ('(A,F8.3,A)'), 'Simulation time: ',Tk(j),' is complete.'
			end if
        end do

        call mpih%destroyMPIHandler
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
		type(mpiHandler) :: mpih
		character(len=100) :: dir, Tstr, filename
		integer, parameter :: N=70
		real(mp) :: temp(2), J0, grad
		real(mp) :: fk
        real(mp), dimension(N) :: fk_array, ek, J1
		real(mp) :: T(6), Time
		integer :: i
        T = (/ ( 0.1_mp*1500.0_mp**(1.0_mp*(i-1)/(size(T)-1)),i=1,size(T) ) /)
!        Time = T(5)
        Time = getOption('adjoint_convergence/time',T(6))

		call mpih%buildMPIHandler
        if( mpih%my_rank.eq.0 .and.                                          &
            mpih%size.ne.(N+1) ) then
            print *, 'Required same number of processors!'
            print *, 'finite difference iteration: ', N+1
            print *, 'simulation times: ', Time
            print *, 'required processors: ',(N+1)

            call mpih%destroyMPIHandler
            stop
        end if

		call mpih%allocateBuffer(N+1,2)
        fk = 1.0_mp*MOD(mpih%my_rank,N+1)
		fk = 10.0_mp**(1.0_mp-0.25_mp*(fk-1))

		dir = getOption('adjoint_convergence/directory','debye_adj_test/dp')

        if( MOD(mpih%my_rank,N+1).eq.0 ) then
           call problem(0.0_mp,Time,trim(dir),0,temp)
            mpih%sendbuf(1,:) = temp
        else
            call problem(fk,Time,trim(dir),1,temp)
            mpih%sendbuf(1,:) = (/ fk, temp(1) /)
        end if

		call mpih%gatherData

		if( mpih%my_rank.eq.mpih%size-1 ) then
            fk_array = mpih%recvbuf(2:N+1,1)
            J0 = mpih%recvbuf(1,1)
            grad = mpih%recvbuf(1,2)
            J1 = mpih%recvbuf(2:N+1,2)
			ek = ABS( ((J1-J0)/fk_array - grad)/grad )
            mpih%recvbuf(2:N+1,2) = ek

            write(Tstr,'(F4.1)') Time
            filename = 'data/'//trim(dir)//'/grad_convergence.T'//trim(adjustl(Tstr))//'.dat'
			open(unit=301, file=filename, status='replace')
            do i=1,N+1
		        write(301,*) mpih%recvbuf(i,:)
            end do
		    close(301)
            print ('(A,F8.3,A)'), 'Simulation time: ',Time,' is complete.'
		end if

        call mpih%destroyMPIHandler
	end subroutine

end module
