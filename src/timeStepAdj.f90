module timeStepAdj

	use timeStep

	implicit none

contains
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

end module
