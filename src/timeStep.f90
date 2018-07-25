module timeStep

	use modSource, only: source, Null_source
	use modTarget, only: control, Null_input
	use modRecord
	use ArMCC
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
		procedure(control), pointer :: PtrControl=>NULL()
		procedure(source), pointer :: PtrSource=>NULL()
		procedure(QoI), pointer :: PtrQoI=>NULL()
		integer :: i,k
		real(mp) :: Jtemp, J_hist(this%nt)
		if( PRESENT(inputControl) ) then
                        PtrControl=>inputControl
                else
                        PtrControl=>Null_input
                end if
		if( PRESENT(inputSource) ) then
                        PtrSource=>inputSource
                else 
                        PtrSource=>Null_source
                end if       
		if( PRESENT(inputQoI) ) then
                        PtrQoI=>inputQoI
                else
                        PtrQoI=>Null_QoI
                end if
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

		call PtrSource(this,k)

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
			r%cpt_temp(1:7) = r%cpt_temp(1:7) + cpt_temp
		else
			call this%mcc_collision(this%p,this%A0)
		end if
	end subroutine

end module
