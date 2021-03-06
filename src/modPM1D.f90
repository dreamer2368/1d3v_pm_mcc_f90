module modPM1D

	use modBC
	use modAssign
	use ArMCC

	implicit none

	type PM1D
		integer :: nt, ni, n, ng, pBCindex, mBCindex
		real(mp) :: eps0, wp
		real(mp) :: dt, L
		real(mp), allocatable :: A0(:)

		type(species), allocatable :: p(:)
		type(mesh) :: m
		type(pmAssign), allocatable :: a(:)

		procedure(applyBC), nopass, pointer :: applyBC
		procedure(mcc_collision), nopass, pointer :: mcc_collision
	contains
		procedure, pass(this) :: buildPM1D
		procedure, pass(this) :: destroyPM1D
	end type

contains

	subroutine buildPM1D(this,Tf,Ti,Ng,N,pBC,mBC,order,dt,L,A,eps)
		class(PM1D), intent(out) :: this
		real(mp), intent(in) :: Tf,Ti
		integer, intent(in) :: Ng, N, pBC, mBC, order
		real(mp), intent(in), optional :: dt, A(:), L, eps
		real(mp) :: L0
		integer :: i
		if( present(dt) ) then
			this%dt = dt
		else
			this%dt = 0.2_mp
		end if
		if( present(A) ) then
			allocate(this%A0(size(A)))
			this%A0 = A
		else
			allocate(this%A0(1))
			this%A0 = 1.0_mp
		end if
		if( present(L) ) then
			this%L = L
		else
			this%L = 2*pi/( sqrt(3.0_mp)/2.0_mp/sqrt(2.0_mp)/0.2_mp )
		end if
		if( present(eps) ) then
			this%eps0 = eps
		else
			this%eps0 = 1.0_mp
		end if
		this%n = N
		this%ng = Ng
		this%pBCindex = pBC
		this%mBCindex = mBC
		this%nt = CEILING(Tf/this%dt)
		this%dt = Tf/this%nt
		this%ni = FLOOR(Ti/this%dt)
		this%wp = 1.0_mp
        if( print_pm_output ) then
	    	print *, 'Plasma is created'
    		print *, 'L = (',this%L,')'
    		print *, 'Ng = (',this%ng,'), N = ',this%n
    		print *, 'Particle BC : ', this%pBCindex
    		print *, 'Mesh BC : ', this%mBCindex
    		print *, 'A = ',this%A0
    		print *, 'Ni = ',this%ni,', Nt = ',this%nt,', dt = ',this%dt
        end if

		allocate(this%p(N))
		call buildMesh(this%m,this%L,Ng,this%mBCindex)
		allocate(this%a(N))
		do i=1,N
			call buildAssign(this%a(i),Ng,order,this%mBCindex)
		end do

		select case(pBC)
			case(0)	!periodic
				this%applyBC=>applyBC_periodic
			case(1) !absorbing-absorbing
				this%applyBC=>applyBC_absorbing
			case(2) !refluxing-absorbing
				this%applyBC=>applyBC_refluxing_absorbing
			case(3) !refluxing-refluxing
				this%applyBC=>applyBC_refluxing_refluxing
		end select

		this%mcc_collision=>no_collision
	end subroutine

	subroutine destroyPM1D(this)
		class(PM1D), intent(inout) :: this
		integer :: i

		deallocate(this%A0)
		do i=1,this%n
			call destroySpecies(this%p(i))
			call destroyAssign(this%a(i))
		end do
		deallocate(this%p)
		call destroyMesh(this%m)
		deallocate(this%a)
	end subroutine


end module
