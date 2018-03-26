module PlasmaProblems

	use init
	use timeStep
	use modMPI
    use modInputHelper

	implicit none

contains

	subroutine Landau_Jtheta
		type(PM1D) :: pm
		type(recordData) :: r
		type(mpiHandler) :: mpih
		integer, parameter :: N = 1, Np = 3E5, Ng = 64
        integer, parameter :: Nsample = 1001, seed = 1
		real(mp), parameter :: vT = 1.0_mp, L = 4.0_mp, dt = 0.1_mp
        real(mp) :: fk(Nsample)
		real(mp) :: dx
		real(mp) :: Time = 20.0_mp
		real(mp) :: A(2),J
		integer :: i, thefile
		character(len=100) :: prefix,dir,filename
		fk = (/ (0.2_mp*(i-1)/(Nsample-1)-0.1_mp,i=1,Nsample) /)

		call buildMPIHandler(mpih)
		call allocateBuffer(Nsample,2,mpih)

        prefix = 'Landau_Jtheta'
        dir = 'data/'//trim(prefix)
        filename = 'Jtheta.bin'
        thefile = mpih%MPIWriteSetup(dir,filename)

		do i=1,mpih%sendcnt
			A = (/ 0.1_mp, fk(mpih%displc(mpih%my_rank)+i) /)
		    call buildPM1D(pm,Time,0.0_mp,Ng,N,0,0,1,dt=dt,L=L,A=A)
			dir = trim(prefix)//'/'//trim(adjustl(mpih%rank_str))
			call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir),20)
			call set_null_discharge(r)
			call init_random_seed(input=seed)
			call Landau_initialize(pm,Np,vT)

			call forwardsweep(pm,r,Te,Null_source,MPE,J)

			mpih%writebuf = (/ fk(mpih%displc(mpih%my_rank)+i),J /)
            call MPI_FILE_WRITE(thefile, mpih%writebuf, 2, MPI_DOUBLE, & 
                                MPI_STATUS_IGNORE, mpih%ierr)
            call MPI_FILE_SYNC(thefile,mpih%ierr)

			call destroyRecord(r)
			call destroyPM1D(pm)
		end do

        call MPI_FILE_CLOSE(thefile, mpih%ierr)            
		call destroyMPIHandler(mpih)
	end subroutine

	subroutine twostream
		type(PM1D) :: pm
		type(adjoint) :: adj
		type(recordData) :: r
		integer, parameter :: Ng=64, Np=10**5, N=1
		real(mp) :: v0 = 0.2_mp, vT = 0.0_mp
		integer :: mode=1
		real(mp) :: J0,J1, grad(1)
		character(len=100)::dir1

		call buildPM1D(pm,150.0_mp,150.0_mp,Ng,N,0,0,1,A=(/1.0_mp,0.0_mp/))
		dir1='twostream_test'
		call buildRecord(r,pm%nt,N,pm%L,Ng,trim(dir1),1)
		call set_null_discharge(r)
		call twostream_initialize(pm,Np,v0,vT,mode)
		call forwardsweep(pm,r,Null_input,Null_source,MPE,J0)
		call printPlasma(r)
		print *, 'J0=',J0

		call destroyPM1D(pm)
		call destroyRecord(r)
	end subroutine

	subroutine Procassini
		type(PM1D) :: sheath
		type(recordData) :: r
		real(mp), parameter :: Kb = 1.38065E-23, EV_TO_K = 11604.52_mp, eps = 8.85418782E-12
		real(mp), parameter :: Te = 10.0_mp*EV_TO_K, tau = 100.0_mp
		real(mp), parameter :: me = 9.10938215E-31, qe = 1.602176565E-19, mu = 1836
		real(mp), parameter :: n0 = 2.0E14
		integer, parameter :: Ne = 5E4, Ni = 5E4
		real(mp) :: mi, Ti, wp0, lambda0, dt, dx, L
		real(mp) :: ve0, vi0, Time_f
		real(mp) :: A(4)
		integer :: i

		mi = mu*me
		Ti = Te/tau
		wp0 = sqrt(n0*qe*qe/me/eps)
		lambda0 = sqrt(eps*Kb*Te/n0/qe/qe)
		L = 20.0_mp*lambda0

		print *, 'L = ',L,', lambda0 = ',lambda0,' e = lambda/L = ',lambda0/L

		dt = 0.1_mp/wp0
		dx = 0.2_mp*lambda0
		!		dt = 0.5_mp*dx/(lambda0*wp0)

		ve0 = sqrt(Kb*Te/me)
		vi0 = sqrt(Kb*Ti/mi)
		Time_f = 1.0_mp*L/vi0

		A = (/ ve0, vi0, 0.2_mp, 1.0_mp*Ni /)
		call buildPM1D(sheath,Time_f,0.0_mp,ceiling(L/dx),2,pBC=2,mBC=2,order=1,A=A,L=L,dt=dt,eps=eps)
		sheath%wp = wp0
		call buildRecord(r,sheath%nt,2,sheath%L,sheath%ng,'Procassini',20)

		call buildSpecies(sheath%p(1),-qe,me)
		call buildSpecies(sheath%p(2),qe,mi)

		call sheath_initialize(sheath,Ne,Ni,Te,Ti,Kb,n0)
		call forwardsweep(sheath,r,Null_input,PartialUniform_Maxwellian2)

		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(sheath)
	end subroutine

	subroutine debye_shielding
		type(PM1D) :: d
		type(recordData) :: r
		real(mp) :: vT
		integer :: i, N, Ng, Nparam
		real(mp) :: L = 20.0_mp, Wp, Q = 2.0_mp
		real(mp) :: dt = 0.05_mp
		real(mp) :: Time = 150.0_mp
		real(mp) :: J
        real(mp), allocatable :: A(:)
        character(len=STRING_LENGTH) :: dir, option
        N = getOption('number_of_particles',100000)
        Ng = getOption('number_of_grids',64)
        Nparam = getOption('number_of_parameters',2)
        dir = getOption('base_directory','Debye')

        allocate(A(Nparam))
        do i=1,Nparam
            write(option,'(A,I02)') 'parameters_of_interest/',i
            A(i) = getOption(trim(option),0.0_mp)
        end do

        A(1) = A(1) + A(2)

		call buildPM1D(d,Time,0.0_mp,Ng,1,pBC=0,mBC=0,order=1,A=A,L=L,dt=dt)
		call buildRecord(r,d%nt,1,d%L,d%ng,trim(dir),10)

		call buildSpecies(d%p(1),-1.0_mp,1.0_mp)
        call init_random_seed
		call Debye_initialize(d,N,Q)

		call forwardsweep(d,r,Null_input,Null_source,Debye,J)

		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(d)

        deallocate(A)
	end subroutine

end module
