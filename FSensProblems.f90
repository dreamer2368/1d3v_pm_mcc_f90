module FSensProblems

	use init
	use timeStepAdj
	use timeStepFSens
	use modMPI

	implicit none

contains

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
		type(PM1D) :: pm
		type(FSens) :: fs
		type(recordData) :: r, fsr
		integer :: N=1E5, Ng=64
		integer :: NInit=5E4, Ngv=32, NInject=5E4, NLimit=5E4
		real(mp) :: L = 20.0_mp, Lv, Q = 2.0_mp
		real(mp) :: dt=0.05_mp, dx
		real(mp) :: Time = 30.0_mp, vT = 1.5_mp
		real(mp) :: A(2), J, grad
		character(len=100)::dir
		A = (/ vT, 0.0_mp /)

		call buildPM1D(pm,Time,0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=dt)
		dir = 'Debye_sensitivity'
		call buildRecord(r,pm%nt,1,pm%L,pm%ng,trim(dir),20)

		call buildSpecies(pm%p(1),-1.0_mp,1.0_mp)
		call Debye_initialize(pm,N,Q)
		
		Lv = vT*6.0_mp
!		NInject = 5*N/pm%nt
		call buildFSens(fs,pm,Lv,Ngv,NInject,NLimit)
		dir = 'Debye_sensitivity/f_A'
		call buildRecord(fsr,fs%nt,1,fs%L,fs%ng,trim(dir),20)
		call Debye_sensitivity_init(fs,N,vT,'vT')
!		call Debye_sensitivity_init_sync(fs,pm,vT,'vT')

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
		integer :: NInit=5E4, Ngv=32, NInject=5E3, NLimit=5E4
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
		type(mpiHandler) :: mpih
		integer, parameter  :: Nsample=1
		real(mp), parameter :: vT=1.5_mp
		integer, parameter :: N = 1E5, Ng = 64
		integer, parameter :: NInit=5E4, Ngv=32, NInject=5E3, NLimit=3E5
		real(mp) :: L = 20.0_mp, Lv, Q = 2.0_mp
		real(mp) :: dt=0.05_mp, dx
		real(mp) :: Time(1), A(2)
		real(mp) :: J,grad,adj_grad(1)
		integer :: i,k
		character(len=100)::dir,Time_str
		Time = (/ 30.0_mp /)
		A = (/ vT, 0.0_mp /)

		call buildMPIHandler(mpih)
!		call allocateBuffer(Nsample,3,mpih)
      call allocateBuffer(Nsample,2,mpih)

		call init_random_seed(mpih%my_rank)
		do k=1,size(Time)
			J = 0.0_mp
			grad = 0.0_mp
			adj_grad = 0.0_mp
			do i=1,mpih%sendcnt
				call buildPM1D(d,Time(k),0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=dt)
				dir = 'Async1/'//trim(adjustl(mpih%rank_str))
				call buildRecord(r,d%nt,1,d%L,d%ng,trim(dir),20)

				call buildSpecies(d%p(1),-1.0_mp,1.0_mp)
				call Debye_initialize(d,N,Q)
		
				Lv = A(1)*6.0_mp
				call buildFSens(fs,d,Lv,Ngv,NInject,NLimit)
				dir = trim(dir)//'/f_A'
				call buildRecord(fsr,fs%nt,1,fs%L,fs%ng,trim(dir),20)
				call Debye_sensitivity_init(fs,2*N,A(1))
	
				call forwardsweep_sensitivity(d,r,fs,fsr,Debye,J,grad)

!				call buildAdjoint(adj,d)
!				call adj%m%setMesh(d%m%rho_back)
!				call backward_sweep(adj,d,r,adj_grad,dDebye,Null_dinput,dDebye_dvT,Null_input,Null_source)

!				mpih%sendbuf(i,:) = (/J,grad,adj_grad(1)/)
            mpih%sendbuf(i,:) = (/J,grad/)

				call destroyRecord(r)
				call destroyPM1D(d)
!				call destroyAdjoint(adj)
				call destroyRecord(fsr)
				call destroyFSens(fs)
			end do

			call gatherData(mpih)

			if( mpih%my_rank.eq.mpih%size-1 ) then
				write(Time_str,'(F5.2)') Time(k)
				open(unit=301,file='data/Async1/Jk_'//trim(Time_str)//'_N105.bin',	&
						status='replace',form='unformatted',access='stream')
				open(unit=302,file='data/Async1/gradk_'//trim(Time_str)//'_N105.bin',	&
						status='replace',form='unformatted',access='stream')
!				open(unit=303,file='data/Debye_sampling/adjk_'//trim(Time_str)//'.bin',	&
!						status='replace',form='unformatted',access='stream')
			   write(301) mpih%recvbuf(:,1)
			   write(302) mpih%recvbuf(:,2)
!				write(303) mpih%recvbuf(:,3)
			   close(301)
			   close(302)
!				close(303)
			end if
		end do

		call destroyMPIHandler(mpih)
	end subroutine

end module
