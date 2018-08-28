module FSensProblems

	use init
	use timeStepAdj
	use timeStepFSens
	use modMPI
    use modInputHelper

	implicit none

contains

	subroutine debye_characterization
		type(PM1D) :: d
		type(recordData) :: r
		type(adjoint) :: adj
		real(mp) ::  vT_min, vT_max, vT_target
        real(mp), allocatable :: vT(:)
		integer :: Nsample, N = 100000, Ng = 64
		real(mp) :: L = 20.0_mp, Wp, Q = 2.0_mp
		real(mp) :: dx
		real(mp) :: Time
		real(mp) :: A(2),J,adj_grad(1)
		integer :: i, thefile, idx, input
		character(len=100):: dir, filename
        logical :: sensitivity
        Time = getOption('QoI_curve/time',150.0_mp)
        dir = getOption('QoI_curve/directory','Debye_curve')
        filename = getOption('QoI_curve/filename','J.bin')
        input = getOption('QoI_curve/random_seed',0)
        sensitivity = getOption('QoI_curve/sensitivity', .false.)
        vT_target = getOption('QoI_curve/sensitivity_measurement_point',1.5_mp)
        vT_min = getOption('QoI_curve/min_parameter_value',1.49_mp)
        vT_max = getOption('QoI_curve/max_parameter_value',1.51_mp)
        Nsample = getOption('QoI_curve/number_of_sample',1001)

        allocate(vT(Nsample))
		vT = (/ ((vT_max-vT_min)*(i-1)/(Nsample-1)+vT_min,i=1,Nsample) /)
        idx = MINLOC( ABS(vT-vT_target), DIM=1 )

		call allocateBuffer(1001,2,mpih)
        thefile = MPIWriteSetup(mpih,'data/'//trim(dir),filename)

		do i=1,mpih%sendcnt
			A = (/ vT(mpih%displc(mpih%my_rank)+i), 0.0_mp /)
            if( A(1).ge.0.5 ) then
                Ng = 64*3
                N = 3E5
            end if
			call buildPM1D(d,Time,0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=0.05_mp)
			call buildRecord(r,d%nt,1,d%L,d%ng,trim(dir)//'/'//trim(adjustl(mpih%rank_str)),20)

			call buildSpecies(d%p(1),-1.0_mp,1.0_mp)
		    call init_random_seed(input=input)
			call Debye_initialize(d,N,Q)

			call forwardsweep(d,r,Null_input,Null_source,Debye,J)

			mpih%writebuf = (/vT(mpih%displc(mpih%my_rank)+i),J/)

            call MPI_FILE_WRITE(thefile, mpih%writebuf, 2, MPI_DOUBLE, & 
                                MPI_STATUS_IGNORE, mpih%ierr)
            call MPI_FILE_SYNC(thefile,mpih%ierr)

            if( sensitivity .and.                                           &
                idx.eq.mpih%displc(mpih%my_rank)+i ) then
				call buildAdjoint(adj,d)
				call adj%m%setMesh(d%m%rho_back)

				call backward_sweep(adj,d,r,adj_grad,dDebye,Null_dinput,dDebye_dvT,Null_input,Null_source)

				call destroyAdjoint(adj)

                print *, 'grad(vT=',mpih%writebuf(1),', idx=',idx,')=',adj_grad(1)
            end if

			call destroyRecord(r)
			call destroyPM1D(d)

            print ('(A,I5,A,I5,A,F8.3,A,F8.3)'), 'Rank-',mpih%my_rank,      &
                                                 ' Sample-',i,              &
                                                 ', vT=',mpih%writebuf(1),  &
                                                 ', J=',mpih%writebuf(2)
		end do

        deallocate(vT)

        call MPI_FILE_CLOSE(thefile, mpih%ierr)            
	end subroutine

	subroutine Debye_sensitivity
		type(PM1D) :: pm
		type(FSens) :: fs
		type(recordData) :: r, fsr
		integer :: N, Ng
		integer :: NInit=5E4, Ngv(1), NInject, NLimit
		real(mp) :: L = 20.0_mp, Lv(1), Q = 2.0_mp
		real(mp) :: dt=0.05_mp, dx
		real(mp) :: Time, vT = 1.5_mp
		real(mp) :: A(2), J, grad
		character(len=100)::dir
		A = (/ vT, 0.0_mp /)
        Time = getOption('simulation_time',30.0_mp)
        N = getOption('number_of_particles',100000)
        Ng = getOption('number_of_grids',64)
        Ngv = Ng/2
        NInject = getOption('number_of_injecting_particles',N/20)
        NLimit = getOption('population_limit',N/2)
        dir = getOption('base_directory','Debye_sensitivity')

		call buildPM1D(pm,Time,0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=dt)
		call buildRecord(r,pm%nt,1,pm%L,pm%ng,trim(dir),10)

		call buildSpecies(pm%p(1),-1.0_mp,1.0_mp)
		call Debye_initialize(pm,N,Q)
		
		Lv = vT*6.0_mp
!		NInject = 5*N/pm%nt
		call buildFSens(fs,pm,Lv,Ngv,NInject,NLimit)
		call buildRecord(fsr,fs%nt,1,fs%L,fs%ng,trim(dir)//'/f_A',10)
        if( fs%scheme.eq.COLLOCATED ) then
    		call Debye_sensitivity_init_sync(fs,pm,vT,'vT')
        else
	    	call Debye_sensitivity_init(fs,N,vT,'vT')
        end if

		call forwardsweep_sensitivity(pm,r,fs,fsr,Debye,J,grad)

		print *, 'J: ',J
		print *, 'grad: ',grad

		call printPlasma(r)
		call printPlasma(fsr)

		call destroyPM1D(pm)
		call destroyRecord(r)
		call destroyFSens(fs)
		call destroyRecord(fsr)
	end subroutine

	subroutine debye_sensitivity_curve
		type(PM1D) :: d
		type(FSens) :: fs
		type(recordData) :: r,fsr
		type(mpiHandler) :: mpih
		integer, parameter  :: Nsample=101
		real(mp) :: vT(Nsample)
		integer :: N = 100000, Ng = 64
		integer :: NInit=5E4, Ngv(1)=32, NInject=5E3, NLimit=5E4
		real(mp) :: L = 20.0_mp, Lv(1), Q = 2.0_mp
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
			print *, 'vT(',mpih%my_rank,')=',A(1)
			call buildPM1D(d,Time,0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=dt)
			dir = 'Debye_sensitivity_curve/'//trim(adjustl(mpih%rank_str))
			call buildRecord(r,d%nt,1,d%L,d%ng,trim(dir),20)

			call buildSpecies(d%p(1),-1.0_mp,1.0_mp)
			call Debye_initialize(d,N,Q)
		
			Lv = A(1)*6.0_mp
			call buildFSens(fs,d,Lv,Ngv,NInject,NLimit)
			dir = 'Debye_sensitivity_curve/'//trim(adjustl(mpih%rank_str))//'/f_A'
			call buildRecord(fsr,fs%nt,1,fs%L,fs%ng,trim(dir),20)
			call Debye_sensitivity_init(fs,2*N,A(1))
	
			call forwardsweep_sensitivity(d,r,fs,fsr,Debye,J,grad)

			mpih%sendbuf(i,:) = (/vT(mpih%displc(mpih%my_rank)+i),J,grad/)

			call destroyRecord(r)
			call destroyPM1D(d)
			call destroyRecord(fsr)
			call destroyFSens(fs)
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
		type(PM1D) :: d
		type(adjoint) :: adj
		type(FSens) :: fs
		type(recordData) :: r,fsr
		integer, parameter  :: Nsample=1E4
		real(mp), parameter :: vT=1.5_mp
		integer :: N = 1E5, Ng = 64
		integer :: NInit=5E4, Ngv(1)=32, NInject=5E3, NLimit=3E5
		real(mp) :: L = 20.0_mp, Lv(1), Q = 2.0_mp
		real(mp) :: dt=0.05_mp, dx
		real(mp) :: Time, A(2)
		real(mp) :: J,grad,adj_grad(1)
		integer :: i,thefile
		character(len=100) :: prefix,dir,Time_str
        dir = getOption('sensitivity_sampling/directory','debye_sampling')
        Time = getOption('sensitivity_sampling/time',150.0_mp)
        N = getOption('sensitivity_sampling/number_of_particles',100000)
        Ng = getOption('sensitivity_sampling/number_of_grids',64)
        dt = getOption('sensitivity_sampling/timestep_size',0.05_mp)
!        write(Time_str,'(F8.3)'), Time
!        dir = trim(dir)//'/T'//trim(adjustl(Time_str))
		A = (/ vT, 0.0_mp /)

!		call buildMPIHandler(mpih)
		call mpih%allocateBuffer(Nsample,3)
!        call allocateBuffer(Nsample,2,mpih)

        thefile = MPIWriteSetup(mpih,'data/'//trim(dir))

		call init_random_seed(mpih%my_rank,addSystemTime=.true.)
		J = 0.0_mp
		grad = 0.0_mp
		adj_grad = 0.0_mp
		do i=1,mpih%sendcnt
			call buildPM1D(d,Time,0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=dt)
			call buildRecord(r,d%nt,1,d%L,d%ng,trim(dir)//'/'//trim(adjustl(mpih%rank_str)),20)

			call buildSpecies(d%p(1),-1.0_mp,1.0_mp)
			call Debye_initialize(d,N,Q)
	
			Lv = A(1)*6.0_mp
			call buildFSens(fs,d,Lv,Ngv,NInject,NLimit)
			dir = trim(dir)//'/f_A'
			call buildRecord(fsr,fs%nt,1,fs%L,fs%ng,trim(dir),20)
			call Debye_sensitivity_init(fs,2*N,A(1))

			call forwardsweep_sensitivity(d,r,fs,fsr,Debye,J,grad)

			call buildAdjoint(adj,d)
			call adj%m%setMesh(d%m%rho_back)
			call backward_sweep(adj,d,r,adj_grad,dDebye,Null_dinput,dDebye_dvT,Null_input,Null_source)

!			mpih%sendbuf(i,:) = (/J,grad,adj_grad(1)/)
!             mpih%sendbuf(i,:) = (/J,grad/)
            mpih%writebuf = (/ J,grad,adj_grad(1) /)

            call MPI_FILE_WRITE(thefile, mpih%writebuf, 3, MPI_DOUBLE, & 
                                MPI_STATUS_IGNORE, mpih%ierr)
            call MPI_FILE_SYNC(thefile,mpih%ierr)

			call destroyRecord(r)
			call destroyPM1D(d)
			call destroyAdjoint(adj)
			call destroyRecord(fsr)
			call destroyFSens(fs)

            print ('(A,I5,A,I5,A)'), 'Rank-',mpih%my_rank,  &
                                    ' Sample-',i,' is collected.'
        end do

!		call gatherData(mpih)
!
!		if( mpih%my_rank.eq.mpih%size-1 ) then
!			write(Time_str,'(F5.2)') Time(k)
!			open(unit=301,file='data/Async1/Jk_'//trim(Time_str)//'_N105.bin',	&
!					status='replace',form='unformatted',access='stream')
!			open(unit=302,file='data/Async1/gradk_'//trim(Time_str)//'_N105.bin',	&
!					status='replace',form='unformatted',access='stream')
!			open(unit=303,file='data/Debye_sampling/adjk_'//trim(Time_str)//'.bin',	&
!					status='replace',form='unformatted',access='stream')
!		   write(301) mpih%recvbuf(:,1)
!		   write(302) mpih%recvbuf(:,2)
!			write(303) mpih%recvbuf(:,3)
!		   close(301)
!		   close(302)
!			close(303)
!		end if

        call MPI_FILE_CLOSE(thefile, mpih%ierr)            

!		call destroyMPIHandler(mpih)
	end subroutine

end module
