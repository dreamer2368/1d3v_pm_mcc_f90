module timeStep

	use modSource
	use modTarget
	use modFSens
	use modRecord
	use ArMCC
	use modAdj
	use modQoI

	implicit none

contains

	subroutine forwardsweep(this,r,inputControl,inputSource,inputQoI,J)
		type(PM1D), intent(inout) :: this
		type(recordData), intent(inout) :: r
		procedure(control), optional :: inputControl
		procedure(source), optional :: inputSource
		procedure(QoI), optional :: inputQoI
		real(mp), intent(out), optional :: J
		procedure(control), pointer :: PtrControl=>Null_input
		procedure(source), pointer :: PtrSource=>Null_source
		procedure(QoI), pointer :: PtrQoI=>Null_QoI
		integer :: i,k
		real(mp) :: Jtemp, J_hist(this%nt)
		if( PRESENT(inputControl) ) PtrControl=>inputControl
		if( PRESENT(inputSource) ) PtrSource=>inputSource
		if( PRESENT(inputQoI) ) PtrQoI=>inputQoI
		if( present(J) ) then
			J = 0.0_mp
		end if
		k=0
		Jtemp = 0.0_mp
		J_hist = 0.0_mp

		!Time stepping
!		call halfStep(this,target_input)
		call PtrQoI(this,k,Jtemp)
		call r%recordPlasma(this, k)									!record for n=1~Nt
		do k=1,this%nt
			call updatePlasma(this,PtrControl,PtrSource,k,r)
			call PtrQoI(this,k,Jtemp)
			J_hist(k) = Jtemp
			call r%recordPlasma(this, k)									!record for n=1~Nt
		end do
		if( PRESENT(J) ) J=Jtemp
		open(unit=305,file='data/'//r%dir//'/J_hist.bin',	&
					status='replace',form='unformatted',access='stream')
		write(305) J_hist
		close(305)
	end subroutine

	subroutine halfStep(this,PtrControl)
		type(PM1D), intent(inout) :: this
		procedure(control), pointer :: PtrControl
		integer :: i, j
		real(mp) :: rhs(this%ng-1)
		real(mp) :: phi1(this%ng-1)
		real(mp) :: dt, L
		integer :: N,Ng
		integer :: g(this%p(1)%np,2)
		real(mp) :: frac(this%p(1)%np,2)
		dt = this%dt
		L = this%L
		N = this%N
		Ng = this%ng

		do i=1,this%n
			call this%applyBC(this%p(i),this%m,this%dt,this%A0(i))
		end do

		!charge assignment
		this%m%rho = 0.0_mp
		do i=1, this%n
			call this%a(i)%chargeAssign(this%p(i),this%m)
		end do

		call PtrControl(this,0,'rho_back')
		call this%m%solveMesh(this%eps0)

		!Electric field : -D*phi
		this%m%E = - multiplyD(this%m%phi,this%m%dx,this%m%BCindex)

		!Force assignment : mat'*E
		do i=1, this%n
			call this%a(i)%forceAssign(this%p(i),this%m)
		end do

		!Half time step advancement in velocity
		do i=1, this%n
			call this%p(i)%accelSpecies(dt/2.0_mp)
		end do
	end subroutine

	subroutine updatePlasma(this,PtrControl,PtrSource,k,r)
		type(PM1D), intent(inout) :: this
		procedure(control), pointer :: PtrControl
		procedure(source), pointer :: PtrSource
		type(recordData), intent(inout), optional :: r
		integer, intent(in) :: k
		real(mp) :: rhs(this%ng-1), phi1(this%ng-1)
		real(mp) :: dt, L
		integer :: N, Ng, i
		real(mp) :: time1, time2, cpt_temp(7)
		dt = this%dt
		L = this%L
		N = this%n
		Ng = this%ng
		cpt_temp = 0.0_mp

		call PtrControl(this,k,'xp')

		call PtrSource(this)

		call CPU_TIME(time1)
		do i=1,this%n
			call this%p(i)%moveSpecies(dt)
		end do
		call CPU_TIME(time2)
		cpt_temp(1) = (time2-time1)

		do i=1,this%n
			call this%applyBC(this%p(i),this%m,this%dt,this%A0(i))
		end do
		call CPU_TIME(time1)
		cpt_temp(2) = (time1-time2)

		!charge assignment
		this%m%rho = 0.0_mp
		do i=1, this%n
			call this%a(i)%chargeAssign(this%p(i),this%m)
		end do
		call CPU_TIME(time2)
		cpt_temp(3) = (time2-time1)

		call PtrControl(this,k,'rho_back')
		call this%m%solveMesh(this%eps0)
		call CPU_TIME(time1)
		cpt_temp(4) = (time1-time2)

		!Electric field : -D*phi
		this%m%E = - multiplyD(this%m%phi,this%m%dx,this%m%BCindex)
		call CPU_TIME(time2)
		cpt_temp(5) = (time2-time1)

		!Force assignment : mat'*E
		do i=1, this%n
			call this%a(i)%forceAssign(this%p(i), this%m)
		end do
		call CPU_TIME(time1)
		cpt_temp(6) = (time1-time2)

		do i=1, this%n
			call this%p(i)%accelSpecies(dt)
		end do
		call CPU_TIME(time2)
		cpt_temp(7) = (time2-time1)

		if( present(r) ) then
			call this%mcc_collision(this%p,this%A0,r%n_coll(:,k))
			r%cpt_temp(1:7) = r%cpt_temp(1:7) + cpt_temp/r%mod
		else
			call this%mcc_collision(this%p,this%A0)
		end if
	end subroutine

!===================Adjoint time stepping==========================

	subroutine backward_sweep(adj,pm,r, grad, inputDJ, inputDControl,inputGrad,inputControl,inputSource)
		type(adjoint), intent(inout) :: adj
		type(PM1D), intent(inout) :: pm
		type(recordData), intent(inout) :: r
		real(mp), intent(out) :: grad(:)
		procedure(Adj_DJ) :: inputDJ
		procedure(Adj_Dcontrol) :: inputDControl
		procedure(Adj_grad) :: inputGrad
		procedure(control) :: inputControl
		procedure(source) :: inputSource
		integer :: k, nk, i
		grad = 0.0_mp

		do k=1,pm%nt
			nk = pm%nt+1-k
			call adj%reset_Dadj

			!=====  Checkpointing  =====
			call checkpoint(pm,r,nk,inputControl,inputSource)

			!===== dJdA : 1st sensitivity calculation ======
			call inputGrad(adj,pm,nk,'before',grad)

			!======= dJ : source term ==========
			call inputDJ(adj,pm,nk)

			!======= dv_p =============
			call inputDControl(adj,pm,nk,'vp')
			call adj%Adj_accel

!			!Check when adjoint reach to the initial step
!			if( k .eq. pm%nt ) then
!				adj%vps = 2.0_mp*adj%vps
!			end if

			!======= dE_p =============
			do i=1,adj%n
				adj%p(i)%Ep = adj%p(i)%qs/pm%p(i)%ms*adj%p(i)%vp(:,1)
			end do

			!======= dE_g =============
			adj%m%E = 0.0_mp
			do i=1,adj%n
				call pm%a(i)%Adj_forceAssign_E(adj%p(i)%Ep,adj%m%E)
			end do
			adj%m%E = adj%m%E + adj%dm%E

			!======= dPhi_g, dRho_g =============
			call solveMesh_Adj(adj%m,pm%eps0)

			!======= dx_p =============
			do i=1,adj%n
				call pm%a(i)%Adj_chargeAssign(pm%p(i),pm%m,adj%m%rho,adj%dp(i)%xp)
				call pm%a(i)%Adj_forceAssign_xp(pm%m,pm%m%E,adj%p(i)%Ep,adj%dp(i)%xp)
			end do
			call inputDControl(adj,pm,nk,'xp')
			call adj%Adj_move

			!===== dJdA : 2nd sensitivity calculation ======
			call inputGrad(adj,pm,nk,'after',grad)
		end do
		nk = 0
		call reset_Dadj(adj)
		!=====  Checkpointing  =====
		call checkpoint(pm,r,nk,inputControl,inputSource)

		!===== dJdA : 1st sensitivity calculation ======
		call inputGrad(adj,pm,nk,'before',grad)

		!======= dJ : source term ==========
		call inputDJ(adj,pm,nk)

		!======= dv_p =============
		call inputDControl(adj,pm,nk,'vp')
		call adj%Adj_accel

		!======= dx_p =============
		call inputDControl(adj,pm,nk,'xp')
		call adj%Adj_move

		!===== dJdA : 2nd sensitivity calculation ======
		call inputGrad(adj,pm,nk,'after',grad)
	end subroutine

	subroutine checkpoint(pm,r,nk,inputControl,inputSource)
		type(PM1D), intent(inout) :: pm
		type(recordData), intent(inout) :: r
		integer, intent(in) :: nk
		procedure(control) :: inputControl
		procedure(source) :: inputSource
		procedure(control), pointer :: PtrControl
		procedure(source), pointer :: PtrSource
		integer :: kr,i
		character(len=1000) :: istr, kstr, dir_temp
		real(mp), allocatable :: xp0(:), vp0(:,:), spwt0(:)
		PtrControl=>inputControl
		PtrSource=>inputSource

		kr = merge(nk,nk/r%mod,r%mod.eq.1)
		write(kstr,*) kr
		do i=1,pm%n
			allocate(xp0(r%np(i,kr+1)))
			allocate(vp0(r%np(i,kr+1),3))
			allocate(spwt0(r%np(i,kr+1)))
			write(istr,*) i
			dir_temp='data/'//r%dir//'/xp/'//trim(adjustl(kstr))//'_'//trim(adjustl(istr))//'.bin'
			open(unit=305,file=trim(dir_temp),form='unformatted',access='stream')
			dir_temp='data/'//r%dir//'/vp/'//trim(adjustl(kstr))//'_'//trim(adjustl(istr))//'.bin'
			open(unit=306,file=trim(dir_temp),form='unformatted',access='stream')
			dir_temp='data/'//r%dir//'/spwt/'//trim(adjustl(kstr))//'_'//trim(adjustl(istr))//'.bin'
			open(unit=307,file=trim(dir_temp),form='unformatted',access='stream')
			read(305) xp0
			read(306) vp0
			read(307) spwt0
			close(305)
			close(306)
			close(307)
			call pm%p(i)%destroySpecies
			call pm%p(i)%setSpecies(r%np(i,kr+1),xp0,vp0,spwt0)
			deallocate(xp0)
			deallocate(vp0)
			deallocate(spwt0)
		end do
		if( nk-kr*r%mod.eq.0 ) then
			!charge assignment
			pm%m%rho = 0.0_mp
			do i=1, pm%n
				call pm%a(i)%chargeAssign(pm%p(i),pm%m)
			end do

			call inputControl(pm,nk,'rho_back')
			call pm%m%solveMesh(pm%eps0)

			!Electric field : -D*phi
			pm%m%E = - multiplyD(pm%m%phi,pm%m%dx,pm%m%BCindex)

			!Force assignment : mat'*E
			do i=1, pm%n
				call pm%a(i)%forceAssign(pm%p(i), pm%m)
			end do
		else
			do i=1,nk-kr*r%mod
				call updatePlasma(pm,PtrControl,PtrSource,i)
			end do
		end if
	end subroutine

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

!			fs%f_A = 0.0_mp
!			do i=1,fs%dpm%n
!				call numberDensity(fs,fs%dpm%p(i),fs%dpm%a(i),fs%f_A)
!			end do
!			do i=1,fs%dpm%n
!				call updateWeight_temp(fs,fs%dpm%p(i),fs%dpm%a(i),fs%f_A,fs%j)
!			end do
			call fs%updateWeight(fs%j)
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
