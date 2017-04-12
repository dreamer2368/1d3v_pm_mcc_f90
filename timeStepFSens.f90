module timeStepFSens

	use timeStep

	implicit none

contains
!===============Forward, continuum Sensitivity

	subroutine forwardsweep_sensitivity(this,r,fs,fsr,inputQoI,J,grad,inputControl,inputSource)
		type(PM1D), intent(inout) :: this
		type(FSens), intent(inout) :: fs
		type(recordData), intent(inout) :: r, fsr
		procedure(QoI) :: inputQoI
		real(mp), intent(out) :: J,grad
		procedure(control), optional :: inputControl
		procedure(source), optional :: inputSource
		integer :: i,k,kr
		character(len=100) :: kstr
		real(mp), dimension(this%nt) :: J_hist, grad_hist
		real(mp) :: time1, time2
		procedure(control), pointer :: PtrControl=>Null_input
		procedure(source), pointer :: PtrSource=>Null_source
		if( PRESENT(inputControl) ) PtrControl=>inputControl
		if( PRESENT(inputSource) ) PtrSource=>inputSource
		J = 0.0_mp
		grad = 0.0_mp
		J_hist = 0.0_mp
		grad_hist = 0.0_mp
		k=0

		!Time stepping
		call halfStep(this,PtrControl)
		call inputQoI(this,k,J)
		call r%recordPlasma(this, k)										!record for n=1~Nt

		call halfStep_Sensitivity(fs%dpm,this,PtrControl)
		call fsr%recordPlasma(fs%dpm, k)
		call fs%FSensDistribution

		if( (fsr%mod.eq.1) .or. (mod(k,fsr%mod).eq.0) ) then
			kr = merge(k,k/fsr%mod,fsr%mod.eq.1)
			write(kstr,*) kr
			open(unit=305,file='data/'//fsr%dir//'/'//trim(adjustl(kstr))//'.bin',	&
					status='replace',form='unformatted',access='stream')
			write(305) fs%f_A
			close(305)
		end if
		do k=1,this%nt
			call updatePlasma(this,PtrControl,PtrSource,k,r)
			call inputQoI(this,k,J)
			J_hist(k) = J
			call r%recordPlasma(this, k)									!record for n=1~Nt

			call updateSensitivity(fs%dpm,this,PtrControl,PtrSource,k,fsr)
!			call fs%FSensDistribution
!			call fs%Redistribute

			call CPU_TIME(time1)
			call fs%FSensSourceTerm(this)
			call CPU_TIME(time2)
			fsr%cpt_temp(8) = fsr%cpt_temp(8) + (time2-time1)/fsr%mod

			fs%f_A = 0.0_mp
			do i=1,fs%dpm%n
				call numberDensity(fs,fs%dpm%p(i),fs%dpm%a(i),fs%f_A)
			end do
			do i=1,fs%dpm%n
				call updateWeight_temp(fs,fs%dpm%p(i),fs%dpm%a(i),fs%f_A,fs%j)
			end do
!			call fs%updateWeight(fs%j)
			call CPU_TIME(time1)
			fsr%cpt_temp(9) = fsr%cpt_temp(9) + (time1-time2)/fsr%mod

!			call fs%InjectSource(fs%j,fs%NInject)
			call inputQoI(fs%dpm,k,grad)
			grad_hist(k) = grad
			call fsr%recordPlasma(fs%dpm, k)

			if( (fsr%mod.eq.1) .or. (mod(k,fsr%mod).eq.0) ) then
				kr = merge(k,k/fsr%mod,fsr%mod.eq.1)
				write(kstr,*) kr
				open(unit=305,file='data/'//fsr%dir//'/'//trim(adjustl(kstr))//'.bin',	&
						status='replace',form='unformatted',access='stream')
				call fs%FSensDistribution
				write(305) fs%f_A
				close(305)
				open(unit=305,file='data/'//fsr%dir//'/j_'//trim(adjustl(kstr))//'.bin',	&
						status='replace',form='unformatted',access='stream')
				write(305) fs%j
				close(305)
			end if
		end do
		open(unit=305,file='data/'//fsr%dir//'/grad_hist.bin',	&
					status='replace',form='unformatted',access='stream')
		write(305) grad_hist
		close(305)
		open(unit=305,file='data/'//fsr%dir//'/J_hist.bin',	&
					status='replace',form='unformatted',access='stream')
		write(305) J_hist
		close(305)
	end subroutine

	subroutine halfStep_Sensitivity(dpm,pm,PtrControl)
		type(PM1D), intent(inout) :: dpm
		type(PM1D), intent(in) :: pm
		procedure(control), pointer :: PtrControl
		integer :: i, j
		real(mp) :: rhs(dpm%ng-1)
		real(mp) :: phi1(dpm%ng-1)
		real(mp) :: dt, L
		integer :: N,Ng
		dt = dpm%dt
		L = dpm%L
		N = dpm%N
		Ng = dpm%ng

!		do i=1, dpm%n
!			dpm%p(i)%xp = pm%p(i)%xp
!			dpm%p(i)%vp = pm%p(i)%vp
!			dpm%a(i)%g = pm%a(i)%g
!			dpm%a(i)%frac = pm%a(i)%frac
!		end do

		do i=1,dpm%n
			call dpm%applyBC(dpm%p(i),dpm%m,dt,dpm%A0(i))
		end do

		!charge assignment
		dpm%m%rho = 0.0_mp
		do i=1, dpm%n
			call dpm%a(i)%chargeAssign(dpm%p(i),dpm%m)
		end do

		call PtrControl(dpm,0,'rho_back')

		!for control parameter qp
!		dpm%m%rho = dpm%m%rho + pm%m%rho/pm%p(1)%qs

		call dpm%m%solveMesh(dpm%eps0)

		!Electric field : -D*phi
		dpm%m%E = - multiplyD(dpm%m%phi,dpm%m%dx,dpm%m%BCindex)

		!Force assignment : mat'*E
		do i=1, dpm%n
			call dpm%a(i)%forceAssign(dpm%p(i),pm%m)
		end do

		!Half time step advancement in velocity
		do i=1, dpm%n
			call dpm%p(i)%accelSpecies(dt/2.0_mp)
		end do
	end subroutine

	subroutine updateSensitivity(dpm,pm,PtrControl,PtrSource,k,r)
		type(PM1D), intent(inout) :: dpm
		type(PM1D), intent(in) :: pm
		procedure(control), pointer :: PtrControl
		procedure(source), pointer :: PtrSource
		type(recordData), intent(inout), optional :: r
		integer, intent(in) :: k
		real(mp) :: rhs(dpm%ng-1), phi1(dpm%ng-1)
		real(mp) :: dt, L
		integer :: N, Ng, i
		real(mp) :: time1, time2
		dt = dpm%dt
		L = dpm%L
		N = dpm%n
		Ng = dpm%ng

		call PtrControl(dpm,k,'xp')

		call PtrSource(dpm)

		call CPU_TIME(time1)
		do i=1,dpm%n
			call dpm%p(i)%moveSpecies(dt)
		end do

!		do i=1, dpm%n
!			dpm%p(i)%xp = pm%p(i)%xp
!			dpm%p(i)%vp = pm%p(i)%vp
!			dpm%a(i)%g = pm%a(i)%g
!			dpm%a(i)%frac = pm%a(i)%frac
!		end do

		call CPU_TIME(time2)
		r%cpt_temp(1) = r%cpt_temp(1) + (time2-time1)/r%mod

		do i=1,dpm%n
			call dpm%applyBC(dpm%p(i),dpm%m,dt,dpm%A0(i))
		end do
		call CPU_TIME(time1)
		r%cpt_temp(2) = r%cpt_temp(2) + (time1-time2)/r%mod

		!charge assignment
		dpm%m%rho = 0.0_mp
		do i=1, dpm%n
			call dpm%a(i)%chargeAssign(dpm%p(i),dpm%m)
		end do
		call CPU_TIME(time2)
		r%cpt_temp(3) = r%cpt_temp(3) + (time2-time1)/r%mod

		call PtrControl(dpm,k,'rho_back')

		!for control parameter qp
!		dpm%m%rho = dpm%m%rho + pm%m%rho/pm%p(1)%qs

		call dpm%m%solveMesh(dpm%eps0)
		call CPU_TIME(time1)
		r%cpt_temp(4) = r%cpt_temp(4) + (time1-time2)/r%mod

		!Electric field : E_A = -D*phi_A
		dpm%m%E = - multiplyD(dpm%m%phi,dpm%m%dx,dpm%m%BCindex)
		call CPU_TIME(time2)
		r%cpt_temp(5) = r%cpt_temp(5) + (time2-time1)/r%mod

		!Force assignment : mat'*E  (NOT E_A !!)
		do i=1, dpm%n
			call dpm%a(i)%forceAssign(dpm%p(i), pm%m)
		end do
		call CPU_TIME(time1)
		r%cpt_temp(6) = r%cpt_temp(6) + (time1-time2)/r%mod

		do i=1, dpm%n
			call dpm%p(i)%accelSpecies(dt)
		end do
		call CPU_TIME(time2)
		r%cpt_temp(7) = r%cpt_temp(7) + (time2-time1)/r%mod

!		if( present(r) ) then
!			call mcc_collision(this,r%n_coll(:,k))
!		else
!			call mcc_collision(this)
!		end if
	end subroutine

end module