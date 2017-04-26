module timeStepFSens

	use timeStep

	implicit none

contains
!===============Forward, continuum Sensitivity

	subroutine forwardsweep_sensitivity(this,r,dpm,dr,inputQoI,J,grad,inputControl,inputSource)
		type(PM1D), intent(inout) :: this
		type(FSens), intent(inout) :: dpm
		type(recordData), intent(inout) :: r, dr
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

		call halfStep_Sensitivity(dpm,this,PtrControl)
		call dr%recordPlasma(dpm, k)

!		call dpm%FSensDistribution
!		kr = merge(k,k/fsr%mod,fsr%mod.eq.1)
!		write(kstr,*) kr
!		open(unit=305,file='data/'//fsr%dir//'/'//trim(adjustl(kstr))//'.bin',	&
!				status='replace',form='unformatted',access='stream')
!		write(305) fs%f_A
!		close(305)
		do k=1,this%nt
			call updatePlasma(this,PtrControl,PtrSource,k,r)
			call inputQoI(this,k,J)
			J_hist(k) = J
			call r%recordPlasma(this, k)									!record for n=1~Nt

			call updateSensitivity(dpm,this,PtrControl,PtrSource,k,dr)

			do i=1,dpm%n
				!For Synchronized weight-updating
!				dpm%p(i)%xp=this%p(i)%xp
!				dpm%p(i)%vp=this%p(i)%vp

				call CPU_TIME(time1)
				call dpm%FSensDistribution(this%p(i),this%a(i))
				!Synchronized
!				call dpm%FSensSourceTerm(this%p(i)%qs,this%p(i)%ms,dpm%j,this%m%E)
				!Asynchronous
				call dpm%FSensSourceTerm(this%p(i)%qs,this%p(i)%ms,dpm%f_A,this%m%E)
				call CPU_TIME(time2)
				dr%cpt_temp(8) = dr%cpt_temp(8) + (time2-time1)

				!InjectSource+Remeshing
!				call dpm%Redistribute(dpm%p(i),dpm%a(i))
				call dpm%InjectSource(dpm%p(i),dpm%j)
				call CPU_TIME(time1)
				dr%cpt_temp(9) = dr%cpt_temp(9) + (time1-time2)

				call dpm%Redistribute_temp(dpm%p(i))

				!Weight updating
            !Synchronized
!				call dpm%updateWeight(dpm%p(i),this%a(i),dpm%f_A,dpm%j)
            !Asynchronous
!				dpm%f_A = 0.0_mp
!				call dpm%numberDensity(dpm%p(i),dpm%a(i),dpm%f_A)
!				call dpm%updateWeight(dpm%p(i),dpm%a(i),dpm%f_A,dpm%j)
				call CPU_TIME(time2)
				dr%cpt_temp(10) = dr%cpt_temp(10) + (time2-time1)
			end do

			call inputQoI(dpm,k,grad)
			grad_hist(k) = grad
			call dr%recordPlasma(dpm, k)

			if( (dr%mod.eq.1) .or. (mod(k,dr%mod).eq.0) ) then
				kr = merge(k,k/dr%mod,dr%mod.eq.1)
				write(kstr,*) kr
				open(unit=305,file='data/'//dr%dir//'/'//trim(adjustl(kstr))//'.bin',	&
						status='replace',form='unformatted',access='stream')
!				call dpm%FSensDistribution(dpm%p(1),this%a(1))
				call dpm%a(1)%chargeAssign(dpm%p(1),dpm%m)
				call dpm%FSensDistribution(dpm%p(1),dpm%a(1))
				write(305) dpm%f_A
				close(305)
				open(unit=305,file='data/'//dr%dir//'/j_'//trim(adjustl(kstr))//'.bin',	&
						status='replace',form='unformatted',access='stream')
				write(305) dpm%j
				close(305)
			end if
		end do
		open(unit=305,file='data/'//dr%dir//'/grad_hist.bin',	&
					status='replace',form='unformatted',access='stream')
		write(305) grad_hist
		close(305)
		open(unit=305,file='data/'//dr%dir//'/J_hist.bin',	&
					status='replace',form='unformatted',access='stream')
		write(305) J_hist
		close(305)
	end subroutine

	subroutine halfStep_Sensitivity(dpm,pm,PtrControl)
		type(FSens), intent(inout) :: dpm
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
		type(FSens), intent(inout) :: dpm
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

		call CPU_TIME(time2)
		r%cpt_temp(1) = r%cpt_temp(1) + (time2-time1)

		do i=1,dpm%n
			call dpm%applyBC(dpm%p(i),dpm%m,dt,dpm%A0(i))
		end do
		call CPU_TIME(time1)
		r%cpt_temp(2) = r%cpt_temp(2) + (time1-time2)

		!charge assignment
		dpm%m%rho = 0.0_mp
		do i=1, dpm%n
			call dpm%a(i)%chargeAssign(dpm%p(i),dpm%m)
		end do
		call CPU_TIME(time2)
		r%cpt_temp(3) = r%cpt_temp(3) + (time2-time1)

		call PtrControl(dpm,k,'rho_back')

		!for control parameter qp
!		dpm%m%rho = dpm%m%rho + pm%m%rho/pm%p(1)%qs

		call dpm%m%solveMesh(dpm%eps0)
		call CPU_TIME(time1)
		r%cpt_temp(4) = r%cpt_temp(4) + (time1-time2)

		!Electric field : E_A = -D*phi_A
		dpm%m%E = - multiplyD(dpm%m%phi,dpm%m%dx,dpm%m%BCindex)
		call CPU_TIME(time2)
		r%cpt_temp(5) = r%cpt_temp(5) + (time2-time1)

		!Force assignment : mat'*E  (NOT E_A !!)
		do i=1, dpm%n
			call dpm%a(i)%forceAssign(dpm%p(i), pm%m)
		end do
		call CPU_TIME(time1)
		r%cpt_temp(6) = r%cpt_temp(6) + (time1-time2)

		do i=1, dpm%n
			call dpm%p(i)%accelSpecies(dt)
		end do
		call CPU_TIME(time2)
		r%cpt_temp(7) = r%cpt_temp(7) + (time2-time1)

!		if( present(r) ) then
!			call mcc_collision(this,r%n_coll(:,k))
!		else
!			call mcc_collision(this)
!		end if
	end subroutine

	subroutine updateSensitivity_sync(dpm,pm,PtrControl,PtrSource,k,r)
		type(FSens), intent(inout) :: dpm
		type(PM1D), intent(in) :: pm
		procedure(control), pointer :: PtrControl
		procedure(source), pointer :: PtrSource
		type(recordData), intent(inout), optional :: r
		integer, intent(in) :: k
		real(mp) :: rhs(dpm%ng-1), phi1(dpm%ng-1)
		real(mp) :: dt, L
		integer :: N, Ng, i,j
		real(mp) :: time1, time2
		integer, dimension(:,:), pointer :: g
		real(mp), dimension(:,:), pointer :: frac
		real(mp), allocatable :: spwt(:)
		dt = dpm%dt
		L = dpm%L
		N = dpm%n
		Ng = dpm%ng

		call PtrControl(dpm,k,'xp')

		call PtrSource(dpm)

		call CPU_TIME(time1)

		call CPU_TIME(time2)
		r%cpt_temp(1) = r%cpt_temp(1) + (time2-time1)/r%mod

		call CPU_TIME(time1)
		r%cpt_temp(2) = r%cpt_temp(2) + (time1-time2)/r%mod

		!charge assignment
		dpm%m%rho = 0.0_mp
		do i=1, dpm%n
			g=>pm%a(i)%g
			frac=>pm%a(i)%frac
			allocate(spwt(dpm%p(i)%np))
			spwt=dpm%p(i)%spwt
			do j=1,dpm%p(i)%np
				dpm%m%rho( g(:,j) ) = dpm%m%rho( g(:,j) ) + spwt(j)*dpm%p(i)%qs/dpm%m%dx*frac(:,j)
			end do
			deallocate(spwt)
!			call dpm%a(i)%chargeAssign(dpm%p(i),dpm%m)
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

		call CPU_TIME(time1)
		r%cpt_temp(6) = r%cpt_temp(6) + (time1-time2)/r%mod

		call CPU_TIME(time2)
		r%cpt_temp(7) = r%cpt_temp(7) + (time2-time1)/r%mod
	end subroutine

end module
