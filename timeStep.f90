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

	subroutine forwardsweep(this,r,target_input,source,QoI,J)
		type(PM1D), intent(inout) :: this
		type(recordData), intent(inout) :: r
		integer :: i,k
		real(mp), intent(out), optional :: J
		real(mp) :: J_hist(this%nt)
		interface
			subroutine target_input(pm,k,str)
				use modPM1D
				type(PM1D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine source(pm)
				use modPM1D
				type(PM1D), intent(inout) :: pm
			end subroutine
		end interface
		interface
			subroutine QoI(pm,k,J)
				use modQoI
				type(PM1D), intent(in) :: pm
				real(mp), intent(inout) :: J
				integer, intent(in) :: k
			end subroutine
		end interface
		optional :: QoI
		if( present(J) ) then
			J = 0.0_mp
		end if
		k=0
		J_hist = 0.0_mp

		!Time stepping
!		call halfStep(this,target_input)
		if( present(J) ) then
			call QoI(this,k,J)
		end if
		call r%recordPlasma(this, k)									!record for n=1~Nt
		do k=1,this%nt
			call updatePlasma(this,target_input,source,k,r)
			if( present(J) ) then
				call QoI(this,k,J)
			end if
			call r%recordPlasma(this, k)									!record for n=1~Nt
			J_hist(k) = J
		end do
		open(unit=305,file='data/'//r%dir//'/J_hist.bin',	&
					status='replace',form='unformatted',access='stream')
		write(305) J_hist
		close(305)
	end subroutine

	subroutine halfStep(this,target_input)
		type(PM1D), intent(inout) :: this
		integer :: i, j
		real(mp) :: rhs(this%ng-1)
		real(mp) :: phi1(this%ng-1)
		real(mp) :: dt, L
		integer :: N,Ng
		interface
			subroutine target_input(pm,k,str)
				use modPM1D
				type(PM1D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		dt = this%dt
		L = this%L
		N = this%N
		Ng = this%ng

		call applyBC(this)
		do i=1,this%n
			call this%a(i)%assignMatrix(this%m,this%p(i)%xp)
		end do
		call adjustGrid(this)

		!charge assignment
		call chargeAssign(this%a,this%p,this%m)

		call target_input(this,0,'rho_back')
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

	subroutine updatePlasma(this,target_input,source,k,r)
		type(PM1D), intent(inout) :: this
		type(recordData), intent(inout), optional :: r
		integer, intent(in) :: k
		real(mp) :: rhs(this%ng-1), phi1(this%ng-1)
		real(mp) :: dt, L
		integer :: N, Ng, i
		real(mp) :: time1, time2
		interface
			subroutine target_input(pm,k,str)
				use modPM1D
				type(PM1D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine source(pm)
				use modPM1D
				type(PM1D), intent(inout) :: pm
			end subroutine
		end interface
		dt = this%dt
		L = this%L
		N = this%n
		Ng = this%ng

		call target_input(this,k,'xp')

		call source(this)

		call CPU_TIME(time1)
		do i=1,this%n
			call this%p(i)%moveSpecies(dt)
		end do
		call CPU_TIME(time2)
		r%cpt_temp(1) = r%cpt_temp(1) + (time2-time1)/r%mod

		call applyBC(this)
		do i=1, this%n
			call this%a(i)%assignMatrix(this%m,this%p(i)%xp)
		end do
		call adjustGrid(this)
		call CPU_TIME(time1)
		r%cpt_temp(2) = r%cpt_temp(2) + (time1-time2)/r%mod

		!charge assignment
		call chargeAssign(this%a,this%p,this%m)
		call CPU_TIME(time2)
		r%cpt_temp(3) = r%cpt_temp(3) + (time2-time1)/r%mod

		call target_input(this,k,'rho_back')
		call this%m%solveMesh(this%eps0)
		call CPU_TIME(time1)
		r%cpt_temp(4) = r%cpt_temp(4) + (time1-time2)/r%mod

		!Electric field : -D*phi
		this%m%E = - multiplyD(this%m%phi,this%m%dx,this%m%BCindex)
		call CPU_TIME(time2)
		r%cpt_temp(5) = r%cpt_temp(5) + (time2-time1)/r%mod

		!Force assignment : mat'*E
		do i=1, this%n
			call this%a(i)%forceAssign(this%p(i), this%m)
		end do
		call CPU_TIME(time1)
		r%cpt_temp(6) = r%cpt_temp(6) + (time1-time2)/r%mod

		do i=1, this%n
			call this%p(i)%accelSpecies(dt)
		end do
		call CPU_TIME(time2)
		r%cpt_temp(7) = r%cpt_temp(7) + (time2-time1)/r%mod

		if( present(r) ) then
			call mcc_collision(this,r%n_coll(:,k))
		else
			call mcc_collision(this)
		end if
	end subroutine
!
!	subroutine QOI(this,J)
!		type(plasma), intent(in) :: this
!		real(mp), intent(out) :: J
!		J = 1.0_mp/Ng/(Nt-Ni)*SUM(this%Edata(:,Ni+1:Nt)**2)
!	end subroutine
!
!	subroutine QOItwo(this,A,B,J)
!	type(plasma), intent(in) :: this
!	integer, intent(in) :: A, B
!	real(mp), intent(out) :: J
!	J = 1.0_mp/Ng/(B-A)*SUM(this%Edata(:,A+1:B)**2)		!omitted 1/N/T for the sake of machine precision
!	end subroutine
!
!===================Adjoint time stepping==========================

	subroutine backward_sweep(adj,pm,r, grad, dJ, Dtarget_input,dJdA,target_input,source)
		type(adjoint), intent(inout) :: adj
		type(PM1D), intent(inout) :: pm
		type(recordData), intent(inout) :: r
		real(mp), intent(out) :: grad(:)
		integer :: k, nk, i
		interface
			subroutine dJ(adj,pm,k)
				use modPM1D
				use modAdj
				type(adjoint), intent(inout) :: adj
				type(PM1D), intent(in) :: pm
				integer, intent(in) :: k
			end subroutine
		end interface
		interface
			subroutine Dtarget_input(adj,pm,k,str)
				use modPM1D
				use modAdj
				type(adjoint), intent(inout) :: adj
				type(PM1D), intent(in) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine dJdA(adj,pm,k,str,grad)
				use modPM1D
				use modAdj
				type(adjoint), intent(in) :: adj
				type(PM1D), intent(in) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
				real(mp), intent(inout) :: grad(:)
			end subroutine
		end interface
		interface
			subroutine target_input(pm,k,str)
				use modPM1D
				type(PM1D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine source(pm)
				use modPM1D
				type(PM1D), intent(inout) :: pm
			end subroutine
		end interface
		grad = 0.0_mp

		do k=1,pm%nt
			nk = pm%nt+1-k

			call adj%reset_Dadj

			!=====  Checkpointing  =====
			call checkpoint(pm,r,nk,target_input,source)

			!===== dJdA : 1st sensitivity calculation ======
			call dJdA(adj,pm,nk,'before',grad)

			!======= dJ : source term ==========
			call dJ(adj,pm,nk)

			!======= dv_p =============
			call Dtarget_input(adj,pm,nk,'vp')
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
			call Dtarget_input(adj,pm,nk,'xp')
			call adj%Adj_move

			!===== dJdA : 2nd sensitivity calculation ======
			call dJdA(adj,pm,nk,'after',grad)
		end do
		nk = 0
		call reset_Dadj(adj)
		!=====  Checkpointing  =====
		call checkpoint(pm,r,nk,target_input,source)

		!===== dJdA : 1st sensitivity calculation ======
		call dJdA(adj,pm,nk,'before',grad)

		!======= dJ : source term ==========
		call dJ(adj,pm,nk)

		!======= dv_p =============
		call Dtarget_input(adj,pm,nk,'vp')
		call adj%Adj_accel

		!======= dx_p =============
		call Dtarget_input(adj,pm,nk,'xp')
		call adj%Adj_move

		!===== dJdA : 2nd sensitivity calculation ======
		call dJdA(adj,pm,nk,'after',grad)
	end subroutine

	subroutine checkpoint(pm,r,nk,target_input,source)
		type(PM1D), intent(inout) :: pm
		type(recordData), intent(inout) :: r
		integer, intent(in) :: nk
		integer :: kr,i
		character(len=1000) :: istr, kstr, dir_temp
		real(mp), allocatable :: xp0(:), vp0(:,:), spwt0(:)
		interface
			subroutine target_input(pm,k,str)
				use modPM1D
				type(PM1D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine source(pm)
				use modPM1D
				type(PM1D), intent(inout) :: pm
			end subroutine
		end interface

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
			do i=1, pm%n
				call pm%a(i)%assignMatrix(pm%m,pm%p(i)%xp)
			end do
			call adjustGrid(pm)

			!charge assignment
			call chargeAssign(pm%a,pm%p,pm%m)

			call target_input(pm,nk,'rho_back')
			call pm%m%solveMesh(pm%eps0)

			!Electric field : -D*phi
			pm%m%E = - multiplyD(pm%m%phi,pm%m%dx,pm%m%BCindex)

			!Force assignment : mat'*E
			do i=1, pm%n
				call pm%a(i)%forceAssign(pm%p(i), pm%m)
			end do
		else
			do i=1,nk-kr*r%mod
				call updatePlasma(pm,target_input,source,i)
			end do
		end if
	end subroutine

!===============Forward, continuum Sensitivity

	subroutine forwardsweep_sensitivity(this,r,fs,fsr,target_input,source,QoI,J,grad)
		type(PM1D), intent(inout) :: this
		type(FSens), intent(inout) :: fs
		type(recordData), intent(inout) :: r, fsr
		integer :: i,k,kr
		character(len=100) :: kstr
		real(mp), intent(out) :: J,grad
		real(mp), dimension(this%nt) :: J_hist, grad_hist
		real(mp) :: time1, time2
		interface
			subroutine target_input(pm,k,str)
				use modPM1D
				type(PM1D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine source(pm)
				use modPM1D
				type(PM1D), intent(inout) :: pm
			end subroutine
		end interface
		interface
			subroutine QoI(pm,k,J)
				use modQoI
				type(PM1D), intent(in) :: pm
				real(mp), intent(inout) :: J
				integer, intent(in) :: k
			end subroutine
		end interface
		J = 0.0_mp
		grad = 0.0_mp
		J_hist = 0.0_mp
		grad_hist = 0.0_mp
		k=0

		!Time stepping
		call halfStep(this,target_input)
		call QoI(this,k,J)
		call r%recordPlasma(this, k)										!record for n=1~Nt

		call halfStep_Sensitivity(fs%dpm,this,target_input)
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
			call updatePlasma(this,target_input,source,k,r)
			call QoI(this,k,J)
			J_hist(k) = J
			call r%recordPlasma(this, k)									!record for n=1~Nt

			call updateSensitivity(fs%dpm,this,target_input,source,k,fsr)
!			call fs%FSensDistribution
!			call fs%Redistribute

			call CPU_TIME(time1)
			call fs%FSensSourceTerm(this)
			call CPU_TIME(time2)
			fsr%cpt_temp(8) = fsr%cpt_temp(8) + (time2-time1)/fsr%mod

			call fs%updateWeight(fs%j)
			call CPU_TIME(time1)
			fsr%cpt_temp(9) = fsr%cpt_temp(9) + (time1-time2)/fsr%mod

!			call fs%InjectSource(fs%j,fs%NInject)
			call QoI(fs%dpm,k,grad)
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
!				open(unit=305,file='data/'//fsr%dir//'/j_'//trim(adjustl(kstr))//'.bin',	&
!						status='replace',form='unformatted',access='stream')
!				write(305) fs%j
!				close(305)
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

	subroutine halfStep_Sensitivity(dpm,pm,Dtarget_input)
		type(PM1D), intent(inout) :: dpm
		type(PM1D), intent(in) :: pm
		integer :: i, j
		real(mp) :: rhs(dpm%ng-1)
		real(mp) :: phi1(dpm%ng-1)
		real(mp) :: dt, L
		integer :: N,Ng
		interface
			subroutine Dtarget_input(pm,k,str)
				use modPM1D
				type(PM1D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		dt = dpm%dt
		L = dpm%L
		N = dpm%N
		Ng = dpm%ng

		call applyBC(dpm)
		do i=1,dpm%n
			call dpm%a(i)%assignMatrix(dpm%m,dpm%p(i)%xp)
		end do
		call adjustGrid(dpm)

		!charge assignment
		call chargeAssign(dpm%a,dpm%p,dpm%m)

		call Dtarget_input(dpm,0,'rho_back')
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

	subroutine updateSensitivity(dpm,pm,Dtarget_input,Dsource,k,r)
		type(PM1D), intent(inout) :: dpm
		type(PM1D), intent(in) :: pm
		type(recordData), intent(inout), optional :: r
		integer, intent(in) :: k
		real(mp) :: rhs(dpm%ng-1), phi1(dpm%ng-1)
		real(mp) :: dt, L
		integer :: N, Ng, i
		real(mp) :: time1, time2
		interface
			subroutine Dtarget_input(pm,k,str)
				use modPM1D
				type(PM1D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine Dsource(pm)
				use modPM1D
				type(PM1D), intent(inout) :: pm
			end subroutine
		end interface
		dt = dpm%dt
		L = dpm%L
		N = dpm%n
		Ng = dpm%ng

		call Dtarget_input(dpm,k,'xp')

		call Dsource(dpm)

		call CPU_TIME(time1)
		do i=1,dpm%n
			call dpm%p(i)%moveSpecies(dt)
		end do
		call CPU_TIME(time2)
		r%cpt_temp(1) = r%cpt_temp(1) + (time2-time1)/r%mod

		call applyBC(dpm)
		do i=1, dpm%n
			call dpm%a(i)%assignMatrix(dpm%m,dpm%p(i)%xp)
		end do
		call adjustGrid(dpm)
		call CPU_TIME(time1)
		r%cpt_temp(2) = r%cpt_temp(2) + (time1-time2)/r%mod

		!charge assignment: rho_A
		call chargeAssign(dpm%a,dpm%p,dpm%m)
		call CPU_TIME(time2)
		r%cpt_temp(3) = r%cpt_temp(3) + (time2-time1)/r%mod

		call Dtarget_input(dpm,k,'rho_back')
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
