module timeStep

	use modSource
	use modTarget
	use modBC
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

		!Time stepping
!		call halfStep(this,target_input)
		if( present(J) ) then
			call QoI(this,k,J)
		end if
		call recordPlasma(r, this, k)									!record for n=1~Nt
		do k=1,this%nt
			call updatePlasma(this,target_input,source,k,r)
			if( present(J) ) then
				call QoI(this,k,J)
			end if
			call recordPlasma(r, this, k)									!record for n=1~Nt
		end do
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
			call assignMatrix(this%a(i),this%m,this%p(i)%xp)
		end do
		call adjustGrid(this)

		!charge assignment
		call chargeAssign(this%a,this%p,this%m)

		call target_input(this,0,'rho_back')
		call solveMesh(this%m,this%eps0)

		!Electric field : -D*phi
		this%m%E = - multiplyD(this%m%phi,this%m%dx,this%m%BCindex)

		!Force assignment : mat'*E
		do i=1, this%n
			call forceAssign(this%a(i),this%p(i),this%m)
		end do

		!Half time step advancement in velocity
		do i=1, this%n
			call accelSpecies(this%p(i),dt/2.0_mp)
		end do
	end subroutine

	subroutine updatePlasma(this,target_input,source,k,r)
		type(PM1D), intent(inout) :: this
		type(recordData), intent(inout), optional :: r
		integer, intent(in) :: k
		real(mp) :: rhs(this%ng-1), phi1(this%ng-1)
		real(mp) :: dt, L
		integer :: N, Ng, i
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

		do i=1,this%n
			call moveSpecies(this%p(i),dt)
		end do

		call applyBC(this)
		do i=1, this%n
			call assignMatrix(this%a(i),this%m,this%p(i)%xp)
		end do
		call adjustGrid(this)

		!charge assignment
		call chargeAssign(this%a,this%p,this%m)

		call target_input(this,k,'rho_back')
		call solveMesh(this%m,this%eps0)

		!Electric field : -D*phi
		this%m%E = - multiplyD(this%m%phi,this%m%dx,this%m%BCindex)

		!Force assignment : mat'*E
		do i=1, this%n
			call forceAssign(this%a(i), this%p(i), this%m)
		end do

		do i=1, this%n
			call accelSpecies(this%p(i),dt)
		end do

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

			call reset_Dadj(adj)

			!=====  Checkpointing  =====
			call checkpoint(pm,r,nk,target_input,source)

			!===== dJdA : 1st sensitivity calculation ======
			call dJdA(adj,pm,nk,'before',grad)

			!======= dJ : source term ==========
			call dJ(adj,pm,nk)

			!======= dv_p =============
			call Dtarget_input(adj,pm,nk,'vp')
			call Adj_accel(adj)

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
				call Adj_forceAssign_E(pm%a(i),adj%p(i)%Ep,adj%m%E)
			end do
			adj%m%E = adj%m%E + adj%dm%E

			!======= dPhi_g, dRho_g =============
			call solveMesh_Adj(adj%m,pm%eps0)

			!======= dx_p =============
			do i=1,adj%n
				call Adj_chargeAssign(pm%a(i),pm%p(i),pm%m,adj%m%rho,adj%dp(i)%xp)
				call Adj_forceAssign_xp(pm%a(i),pm%m,pm%m%E,adj%p(i)%Ep,adj%dp(i)%xp)
			end do
			call Dtarget_input(adj,pm,nk,'xp')
			call Adj_move(adj)

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
		call Adj_accel(adj)

		!======= dx_p =============
		call Dtarget_input(adj,pm,nk,'xp')
		call Adj_move(adj)

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
			call destroySpecies(pm%p(i))
			call setSpecies(pm%p(i),r%np(i,kr+1),xp0,vp0,spwt0)
			deallocate(xp0)
			deallocate(vp0)
			deallocate(spwt0)
		end do
		if( nk-kr*r%mod.eq.0 ) then
			do i=1, pm%n
				call assignMatrix(pm%a(i),pm%m,pm%p(i)%xp)
			end do
			call adjustGrid(pm)

			!charge assignment
			call chargeAssign(pm%a,pm%p,pm%m)

			call target_input(pm,nk,'rho_back')
			call solveMesh(pm%m,pm%eps0)

			!Electric field : -D*phi
			pm%m%E = - multiplyD(pm%m%phi,pm%m%dx,pm%m%BCindex)

			!Force assignment : mat'*E
			do i=1, pm%n
				call forceAssign(pm%a(i), pm%p(i), pm%m)
			end do
		else
			do i=1,nk-kr*r%mod
				call updatePlasma(pm,target_input,source,i)
			end do
		end if
	end subroutine

end module
