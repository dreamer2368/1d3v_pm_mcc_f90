program main

	use testmodule
	use PlasmaProblems
	use AdjointProblems
	use MCCProblems
	use FSensProblems
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

!	call cross_section
!	call Procassini
!	call test_refluxing_boundary
!	call test_anewvel_Ar
!	call test_mcc_electron
!	call test_mcc_Argon
!	call test_ext_voltage_Poisson
!	call Ar_discharge
!	call test_particle_adj(64,2)
!	call test_backward_sweep
!   call twostream
!	call twostream_adj(output(1),output(2))
!	call Landau(0.0_mp, 6.0_mp ,'Landau', 1,output )
!	call random_test
!   call Landau_adjoint_sampling
!   call twostream_adjoint_sampling
!	call debye_shielding
!	call debye_characterization
!	call InjectionTest
!	call MPITest
!	call SensitivityInitializeTest
!	call Debye_sensitivity
!	call forYeoh
!	call RedistributionTest
!	call updateWeightTest
!	call debye_sensitivity_curve
!	call adj_convergence(debye_adj)
!	call adjoint_convergence_in_time(Landau)
!	call debye_sampling
!	call redistribute_temp_test
!    call MPI_write_test

	! print to screen
    if( mpih%my_rank .eq. 0) then
	    print *, 'program main...done.'
    end if

contains

	! You can add custom subroutines/functions here later, if you want

end program
