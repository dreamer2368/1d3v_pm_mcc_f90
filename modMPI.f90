module modMPI

	use constants

	implicit none

	include 'mpif.h'

	type mpiHandler
		integer :: ierr, my_rank, size
		character(len=100) :: rank_str
		integer :: sendcnt
		real(mp), allocatable :: sendbuf(:,:), recvbuf(:,:), writebuf(:)
		integer, allocatable :: recvcnt(:), displc(:)
	contains
		procedure, pass(this) :: buildMPIHandler
		procedure, pass(this) :: destroyMPIHandler
		procedure, pass(this) :: allocateBuffer
		procedure, pass(this) :: gatherData
	end type

contains

	subroutine buildMPIHandler(this)
		class(mpiHandler), intent(out) :: this

		call MPI_INIT(this%ierr)
		call MPI_COMM_RANK(MPI_COMM_WORLD,this%my_rank,this%ierr)
		call MPI_COMM_SIZE(MPI_COMM_WORLD,this%size,this%ierr)
		write(this%rank_str,*) this%my_rank

	end subroutine

	subroutine destroyMPIHandler(this)
		class(mpiHandler), intent(inout) :: this

		call MPI_FINALIZE(this%ierr)
		if( allocated(this%sendbuf) ) deallocate(this%sendbuf)
		if( allocated(this%recvbuf) ) deallocate(this%recvbuf)
		if( allocated(this%recvcnt) ) deallocate(this%recvcnt)
		if( allocated(this%displc) ) deallocate(this%displc)
	end subroutine

	subroutine allocateBuffer(Nsample,Ndata,this)
		integer, intent(in) :: Nsample, Ndata
		class(mpiHandler), intent(inout) :: this
		integer :: sample_per_core, i

		sample_per_core = Nsample/this%size

		allocate(this%recvcnt(0:this%size-1))
		allocate(this%displc(0:this%size-1))
		this%recvcnt(0:MOD(Nsample,this%size)-1) = sample_per_core+1
		this%recvcnt(MOD(Nsample,this%size):this%size-1) = sample_per_core
		this%displc = 0
		do i=0,this%size-1
			this%displc(i) = SUM(this%recvcnt(0:i-1))
		end do
		if( this%my_rank.eq.this%size-1 ) then
			print *, 'size: ',this%size
			print *, 'sample/core: ',sample_per_core
			print *, 'remainder: ',MOD(Nsample,this%size)
         allocate(this%recvbuf(Nsample,Ndata))
		   print *, this%recvcnt
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

        allocate(this%writebuf(Ndata))
	end subroutine

	subroutine gatherData(this)
		class(mpiHandler), intent(inout) :: this
		integer :: i
		do i=1,size(this%sendbuf,2)
			call MPI_GATHERV(this%sendbuf(:,i),this%sendcnt,MPI_DOUBLE,	&
							this%recvbuf(:,i),this%recvcnt,this%displc,MPI_DOUBLE,	&
							this%size-1,MPI_COMM_WORLD,this%ierr)
		end do
	end subroutine

    function MPIWriteSetup(this, dir) result(thefile)
        class(mpiHandler), intent(inout) :: this
        character(len=*), intent(in) :: dir
        character(len=100) :: filename
        integer :: thefile
        integer(kind=MPI_OFFSET_KIND) :: disp
        integer :: i

		call system('mkdir -p data/'//trim(dir))
        filename = trim(dir)//'/sampling.bin'
            
        call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(filename), & 
                           MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                           MPI_INFO_NULL, thefile, this%ierr)
        disp = this%displc(this%my_rank)*8*size(this%writebuf)
        call MPI_FILE_SET_VIEW(thefile,disp,MPI_DOUBLE,MPI_DOUBLE,   &
                                'native',MPI_INFO_NULL,this%ierr)
    end function

end module
