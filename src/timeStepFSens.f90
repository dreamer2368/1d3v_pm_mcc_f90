module timeStepFSens

	use timeStep
    use modFSens

	implicit none

	abstract interface
		subroutine updateSensitivity(dpm,pm,PtrControl,PtrSource,k,r)
            use modFSens
			use modPM1D
            use modRecord
            use modSource, only: source, Null_source
            use modTarget, only: control, Null_input
            type(FSens), intent(inout) :: dpm
			type(PM1D), intent(in) :: pm
            type(recordData), intent(inout) :: r
            procedure(control), pointer :: PtrControl=>NULL()
            procedure(source), pointer :: PtrSource=>NULL()
			integer, intent(in) :: k
		end subroutine
	end interface

    abstract interface
        subroutine integrateG(dpm,pm,r)
            use modFSens
            use modPM1D
            use modRecord
            type(FSens), intent(inout) :: dpm
            type(PM1D), intent(in) :: pm
            type(recordData), intent(inout) :: r
        end subroutine
    end interface

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
		integer :: i,k,kr,fileUnit
		character(len=100) :: kstr
		real(mp), dimension(this%nt) :: J_hist, grad_hist
		real(mp) :: time1, time2
		procedure(control), pointer :: PtrControl=>NULL()
		procedure(source), pointer :: PtrSource=>NULL()
        procedure(updateSensitivity), pointer :: PtrUpdate=>NULL()
        procedure(integrateG), pointer :: PtrIntegrateG=>NULL()
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
        select case(dpm%scheme)
            case(COLLOCATED)
                PtrUpdate=>updateSensitivityCollocated
                PtrIntegrateG=>integrateGCollocated
            case(NONCOLLOCATED)
                PtrUpdate=>updateSensitivityNonCollocated
                PtrIntegrateG=>integrateGNonCollocated
            case(INJECTION)
                PtrUpdate=>updateSensitivityNonCollocated
                PtrIntegrateG=>integrateGInjection
        end select 
		J = 0.0_mp
		grad = 0.0_mp
		J_hist = 0.0_mp
		grad_hist = 0.0_mp
		k=0
        fileUnit = mpih%my_rank+305
        timeProfile = 0.0_mp
        functionCalls = 0

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

			call PtrUpdate(dpm,this,PtrControl,PtrSource,k,dr)

!			if( (dr%mod.eq.1) .or. (mod(k,dr%mod).eq.0) ) then
!				kr = merge(k,k/dr%mod,dr%mod.eq.1)
!				write(kstr,*) kr
!				open(unit=fileUnit,file='data/'//dr%dir//'/'//trim(adjustl(kstr))//'.bin',	&
!						status='replace',form='unformatted',access='stream')
!                if( dpm%scheme.eq.COLLOCATED ) then
!				    call dpm%psM(1)%FDistribution(dpm%p(1),this%a(1))
!                else
!				    call dpm%psM(1)%FDistribution(dpm%p(1),dpm%a(1))
!                end if
!				write(fileUnit) dpm%psM(1)%f_A
!				close(fileUnit)
!!				open(unit=305,file='data/'//dr%dir//'/j_'//trim(adjustl(kstr))//'.bin',	&
!!						status='replace',form='unformatted',access='stream')
!!				write(305) dpm%j
!!				close(305)
!			end if

            call PtrIntegrateG(dpm,this,dr)

			call inputQoI(dpm,k,grad)
			grad_hist(k) = grad
			call dr%recordPlasma(dpm, k)
		end do
		open(unit=fileUnit,file='data/'//dr%dir//'/grad_hist.bin',	&
					status='replace',form='unformatted',access='stream')
		write(fileUnit) grad_hist
		close(fileUnit)
		open(unit=fileUnit,file='data/'//dr%dir//'/J_hist.bin',	&
					status='replace',form='unformatted',access='stream')
		write(fileUnit) J_hist
		close(fileUnit)
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

!==========================updateSensitivity procedures===========================

	subroutine updateSensitivityNonCollocated(dpm,pm,PtrControl,PtrSource,k,r)
		type(FSens), intent(inout) :: dpm
		type(PM1D), intent(in) :: pm
		procedure(control), pointer :: PtrControl
		procedure(source), pointer :: PtrSource
		type(recordData), intent(inout) :: r
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

		call PtrSource(dpm,k)

		do i=1,dpm%n
			call dpm%p(i)%moveSpecies(dt)
		end do

		do i=1,dpm%n
			call dpm%applyBC(dpm%p(i),dpm%m,dt,dpm%A0(i))
		end do

		!charge assignment
		dpm%m%rho = 0.0_mp
		do i=1, dpm%n
			call dpm%a(i)%chargeAssign(dpm%p(i),dpm%m)
		end do

		call PtrControl(dpm,k,'rho_back')

		!for control parameter qp
!		dpm%m%rho = dpm%m%rho + pm%m%rho/pm%p(1)%qs

        call CPU_TIME(time2)
		call dpm%m%solveMesh(dpm%eps0)
		call CPU_TIME(time1)
		timeProfile(4) = timeProfile(4) + (time1-time2)
        functionCalls(4) = functionCalls(4) + 1

		!Electric field : E_A = -D*phi_A
		dpm%m%E = - multiplyD(dpm%m%phi,dpm%m%dx,dpm%m%BCindex)
		call CPU_TIME(time2)
		timeProfile(5) = timeProfile(5) + (time2-time1)
        functionCalls(5) = functionCalls(5) + 1

		!Force assignment : mat'*E  (NOT E_A !!)
		do i=1, dpm%n
			call dpm%a(i)%forceAssign(dpm%p(i), pm%m)
		end do

		do i=1, dpm%n
			call dpm%p(i)%accelSpecies(dt)
		end do

!		if( present(r) ) then
!			call mcc_collision(this,r%n_coll(:,k))
!		else
!			call mcc_collision(this)
!		end if
	end subroutine

	subroutine updateSensitivityCollocated(dpm,pm,PtrControl,PtrSource,k,r)
		type(FSens), intent(inout) :: dpm
		type(PM1D), intent(in) :: pm
		procedure(control), pointer :: PtrControl
		procedure(source), pointer :: PtrSource
		type(recordData), intent(inout) :: r
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

		call PtrSource(dpm,k)

        call CPU_TIME(time1)
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
		timeProfile(3) = timeProfile(3) + (time2-time1)
        functionCalls(3) = functionCalls(3) + 1

		call PtrControl(dpm,k,'rho_back')

		!for control parameter qp
!		dpm%m%rho = dpm%m%rho + pm%m%rho/pm%p(1)%qs

		call dpm%m%solveMesh(dpm%eps0)
		call CPU_TIME(time1)
		timeProfile(4) = timeProfile(4) + (time1-time2)
        functionCalls(4) = functionCalls(4) + 1

		!Electric field : E_A = -D*phi_A
		dpm%m%E = - multiplyD(dpm%m%phi,dpm%m%dx,dpm%m%BCindex)
		call CPU_TIME(time2)
		timeProfile(5) = timeProfile(5) + (time2-time1)
        functionCalls(5) = functionCalls(5) + 1

	end subroutine

!=======================integrateG procedure============================

    subroutine integrateGNonCollocated(dpm,pm,r)
        type(FSens), intent(inout) :: dpm
        type(PM1D), intent(in) :: pm
        type(recordData), intent(inout) :: r
        integer :: i
        real(mp) :: time1, time2

		do i=1,dpm%n
			call dpm%psM(i)%FVelocityGradient(pm%p(i),pm%a(i))
			call dpm%psM(i)%FSensSourceTerm(pm%p(i)%qs,pm%p(i)%ms,dpm%m%E,dpm%dt)

			call dpm%psM(i)%numberDensity(dpm%p(i),dpm%a(i))
			call dpm%psM(i)%updateWeight(dpm%p(i),dpm%a(i))
		end do
    end subroutine

    subroutine integrateGInjection(dpm,pm,r)
        type(FSens), intent(inout) :: dpm
        type(PM1D), intent(in) :: pm
        type(recordData), intent(inout) :: r
        integer :: i
        real(mp) :: time1, time2

        do i=1,dpm%n
			call dpm%psM(i)%FVelocityGradient(pm%p(i),pm%a(i))
			call dpm%psM(i)%FSensSourceTerm(pm%p(i)%qs,pm%p(i)%ms,dpm%m%E,dpm%dt)

			!InjectSource+Remeshing
			call dpm%InjectSource(dpm%L,dpm%psM(i)%Lv,dpm%p(i),dpm%psM(i)%J)

			call dpm%Redistribute(dpm%p(i),dpm%psM(i)%Lv)
		end do
    end subroutine

    subroutine integrateGCollocated(dpm,pm,r)
        type(FSens), intent(inout) :: dpm
        type(PM1D), intent(in) :: pm
        type(recordData), intent(inout) :: r
        integer :: i
        real(mp) :: time1, time2

		do i=1,dpm%n
			dpm%p(i)%xp=pm%p(i)%xp
			dpm%p(i)%vp=pm%p(i)%vp

			call dpm%psM(i)%FVelocityGradient(pm%p(i),pm%a(i))
			call dpm%psM(i)%FSensSourceTerm(pm%p(i)%qs,pm%p(i)%ms,dpm%m%E,dpm%dt)

			call dpm%psM(i)%numberDensity(dpm%p(i),pm%a(i))
			call dpm%psM(i)%updateWeight(dpm%p(i),pm%a(i))
		end do
    end subroutine

end module
