# 1d3v_pm_mcc_f90
At this version, it is possible for newly created particles to experience additional collisions by collision subroutines for subsequent species.

This is not desired, since more than 2 collision per particle per time step is prohibited and truncated, with overestimation of 1 collision probability per particle.

Checkpointing PIC-adjoint is added for periodic BC.

# NOTE
1. Do not pass a string as itself into the subroutine buildRecord. Instead:

	dir1='example'
	call buildRecord(recordData, Nt, N, L, Ng, trim(dir1), mod)
