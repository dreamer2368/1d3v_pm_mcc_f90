module modMPI

   use constants

	implicit none

	include 'mpif.h'

	type mpiHandler
		integer :: ierr, my_rank, size
		character(len=100) :: rank_str
		integer :: sendcnt
		real(mp), allocatable :: sendbuf(:,:), recvbuf(:,:)
		integer, allocatable :: recvcnt(:), displc(:)
	end type

contains

	subroutine buildMPIHandler(this)
		type(mpiHandler), intent(out) :: this

		call MPI_INIT(this%ierr)
		call MPI_COMM_RANK(MPI_COMM_WORLD,this%my_rank,this%ierr)
		call MPI_COMM_SIZE(MPI_COMM_WORLD,this%size,this%ierr)
		write(this%rank_str,*) this%my_rank
		call MPI_FINALIZE(this%ierr)
	end subroutine

	subroutine destroyMPIHandler(this)
		type(mpiHandler), intent(inout) :: this
		if( allocated(this%sendbuf) ) deallocate(this%sendbuf)
		if( allocated(this%sendbuf) ) deallocate(this%recvbuf)
		if( allocated(this%sendbuf) ) deallocate(this%recvcnt)
		if( allocated(this%sendbuf) ) deallocate(this%displc)
	end subroutine

	subroutine allocateBuffer(Nsample,Ndata,this)
		integer, intent(in) :: Nsample, Ndata
		type(mpiHandler), intent(inout) :: this
		integer :: sample_per_core, i

		sample_per_core = Nsample/this%size

		if( this%my_rank.eq.this%size-1 ) then
			print *, 'size: ',this%size
			print *, 'sample/core: ',sample_per_core
			print *, 'remainder: ',MOD(Nsample,this%size)
			allocate(this%recvcnt(0:this%size-1))
			allocate(this%displc(0:this%size-1))
			this%recvcnt(0:MOD(Nsample,this%size)-1) = sample_per_core+1
			this%recvcnt(MOD(Nsample,this%size):this%size-1) = sample_per_core
			print *, this%recvcnt
			this%displc = 0
			do i=0,this%size-1
				this%displc(i) = SUM(this%recvcnt(0:i-1))
			end do
			print *, this%displc
		end if

		if( this%my_rank<MOD(Nsample,this%size) ) then
			this%sendcnt = sample_per_core+1
			allocate(this%sendbuf(this%sendcnt,Ndata))
		else
			this%sendcnt = sample_per_core
			allocate(this%sendbuf(this%sendcnt,Ndata))
		end if
		this%sendbuf = 0.0_mp
	end subroutine

	subroutine gatherData(this)
		type(mpiHandler), intent(inout) :: this
		integer :: i

		call MPI_INIT(this%ierr)
		do i=1,size(this%sendbuf,2)
			call MPI_GATHERV(this%sendbuf(:,i),this%sendcnt,MPI_DOUBLE,	&
							this%recvbuf(:,i),this%recvcnt,this%displc,MPI_DOUBLE,	&
							this%size-1,MPI_COMM_WORLD,this%ierr)
		end do
		call MPI_FINALIZE(this%ierr)
	end subroutine

end module
