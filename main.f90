program main

	use testmodule
	use PlasmaProblems
	use AdjointProblems
	use MCCProblems
	use FSensProblems

	implicit none

	real(mp) :: output(2) = (/(0.1_mp)**3,0.0_mp/)

	! print to screen
	print *, 'calling program main'

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
!	call adj_convergence(twostream_grad)
!	call adjoint_convergence_in_time(debye_adj)
	call debye_sampling

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

end program
