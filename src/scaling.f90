program scaling

	use testmodule
	use PlasmaProblems
	use AdjointProblems
	use MCCProblems
	use FSensProblems
	use modMPI
	use modInputHelper

	implicit none

	real(mp) :: output(2) = (/(0.1_mp)**3,0.0_mp/)
    character(len=STRING_LENGTH), parameter :: PROJECT_NAME='PASS'
    character(len=STRING_LENGTH) :: filename

    ! initiate MPI
    call mpih%buildMPIHandler

	! print to screen
    if( mpih%my_rank .eq. 0 )           &
    	print *, 'calling program main'

    ! Parse options from the input file.
    filename = trim(PROJECT_NAME) // ".inp"
    call parseInputFile(filename)
    print_pm_output = getOption('print_simulation_detail',.false.)

		call chargeAssignTest

	! print to screen
    if( mpih%my_rank .eq. 0) then
	    print *, 'program main...done.'
    end if

    call mpih%destroyMPIHandler

contains

	! You can add custom subroutines/functions here later, if you want

	subroutine chargeAssignTest
		type(PM1D) :: pm
		real(mp), parameter :: L = 20.0_mp, Lv(1)=5.0_mp, w = 1.0_mp
		integer, parameter :: Np=1E6, rep=20
		real(mp) :: xp0(Np), vp0(Np,3), spwt0(Np)
		integer :: Ng(4), i, j,k
		real(mp) :: time1, time2

		integer :: g(2), gp
		real(mp) :: frac(2)
		real(mp), allocatable :: rho(:)

		Ng = (/ 10, 100, 1000, 10000 /)*1000

		do k = 1,4
			call pm%buildPM1D(1.0_mp,1.0_mp,Ng(k),N=1,pBC=0,mBC=0,order=1,L=L,dt=1.0_mp)
			call pm%p(1)%buildSpecies(1.0_mp,1.0_mp)

			call RANDOM_NUMBER(xp0)
			xp0 = L*xp0
			vp0 = randn(Np,3)
			spwt0 = L/Np
			call pm%p(1)%setSpecies(Np,xp0,vp0,spwt0)

			if( ALLOCATED(rho) ) DEALLOCATE(rho)
			ALLOCATE(rho(Ng(k)))

			CALL RANDOM_NUMBER(frac)
			g = floor( frac*Ng(k) )+1

			do j = 1,rep
				call CPU_TIME(time1)

        ! pm%m%rho = 0.0_mp
				rho = 0.0_mp

				call CPU_TIME(time2)
        timeProfile(1) = timeProfile(1) + (time2-time1)
        functionCalls(1) = functionCalls(1) + 1

				DEALLOCATE(pm%a(1)%g)
				DEALLOCATE(pm%a(1)%frac)
				pm%a(1)%np = pm%p(1)%np
				ALLOCATE(pm%a(1)%g(pm%a(1)%order+1,pm%p(1)%np))
				ALLOCATE(pm%a(1)%frac(pm%a(1)%order+1,pm%p(1)%np))

				call CPU_TIME(time1)

				do i=1,pm%p(1)%np
					CALL pm%a(1)%assignMatrix(pm%p(1)%xp(i),pm%m%dx,g,frac)
					CALL pm%a(1)%adjustGrid(pm%m%ng,g,frac)
					! CALL RANDOM_NUMBER(frac)
					! g = floor( frac*Ng(k) )+1

					pm%a(1)%g(:,i) = g
					pm%a(1)%frac(:,i) = frac
					rho( g ) = rho( g ) + pm%p(1)%spwt(i)*frac
					! rho( g ) = pm%p(1)%spwt(i)*frac
					! rho( g(i) ) = rho( g(i) ) + 1.0_mp
					! gp = g(k)
					! rho( gp ) = 1.0_mp !This one doesn't scale with mesh.
				end do
				call CPU_TIME(time2)
				timeProfile(2) = timeProfile(2) + (time2-time1)
				functionCalls(2) = functionCalls(2) + 1


				call CPU_TIME(time1)

        rho = rho*pm%p(1)%qs/pm%m%dx
        if( pm%a(1)%mBCidx .ne. 0 ) then
            rho( (/1,pm%m%ng/) ) = 2.0_mp*rho( (/1,pm%m%ng/) )
        end if

        call CPU_TIME(time2)
        timeProfile(3) = timeProfile(3) + (time2-time1)
        functionCalls(3) = functionCalls(3) + 1
			end do
			if( mpih%my_rank .eq. 0 ) print *, 'Ng: ',Ng(k),timeProfile(1:3)/functionCalls(1:3)
		end do
	end subroutine

end program
