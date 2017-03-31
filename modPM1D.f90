module modPM1D

	use modParticleBC
	use modAssign

	implicit none

	type PM1D
		integer :: nt, ni, n, ng, ngv, pBCindex, mBCindex
		real(mp) :: eps0, wp
		real(mp) :: dt, L, Lv
		real(mp), allocatable :: A0(:)

		type(species), allocatable :: p(:)
		type(mesh) :: m
		type(pmAssign), allocatable :: a(:)

		procedure(applyBC), nopass, pointer :: applyBC
	contains
		procedure, pass(this) :: buildPM1D
		procedure, pass(this) :: destroyPM1D
	end type

contains

	subroutine buildPM1D(this,Tf,Ti,Ng,N,pBC,mBC,order,dt,L,A,eps,Ngv,Lv)
		class(PM1D), intent(out) :: this
		real(mp), intent(in) :: Tf,Ti
		integer, intent(in) :: Ng, N, pBC, mBC, order
		real(mp), intent(in), optional :: dt, A(:), L, Lv, eps
		integer, intent(in), optional :: Ngv
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
		if( PRESENT(Ngv) ) then
			this%ngv = Ngv
			this%Lv = Lv
		end if
		this%n = N
		this%ng = Ng
		this%pBCindex = pBC
		this%mBCindex = mBC
		this%nt = CEILING(Tf/this%dt)
		this%dt = Tf/this%nt
		this%ni = FLOOR(Ti/this%dt)
		this%wp = 1.0_mp
		print *, 'Plasma is created'
		print *, 'L = (',this%L,')',', Ng_x = (',this%ng,')'
		if( PRESENT(Ngv) ) print *, 'Lv = (',this%Lv,')',', Ng_v = (',this%ngv,')'
		print *, 'Number of species = ',this%n
		print *, 'Particle BC : ', this%pBCindex
		print *, 'Mesh BC : ', this%mBCindex
		print *, 'A = ',this%A0
		print *, 'Ni = ',this%ni,', Nt = ',this%nt,', dt = ',this%dt

		!Allocate array of species
		allocate(this%p(N))

		!Build Mesh
		if( PRESENT(Ngv) ) then
			call buildMesh(this%m,this%L,Ng,this%mBCindex,Lv,Ngv)
		else
			call buildMesh(this%m,this%L,Ng,this%mBCindex)
		end if

		!Build Assignment
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
