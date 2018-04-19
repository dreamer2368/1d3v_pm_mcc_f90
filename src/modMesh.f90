module modMesh

	use MatrixVector

	implicit none

	type mesh
		integer :: ng, BCindex
		real(mp) :: L, dx

		real(mp), allocatable :: E(:)
!		real(mp), dimension(:), pointer :: rho=>NULL()	!For future devolopment of pointer-use
		real(mp), allocatable :: rho(:)
		real(mp), allocatable :: rho_back(:)				!1D sheath: surface charge
		real(mp), allocatable :: phi(:)

        complex(mp), allocatable :: W(:)

		procedure(MeshSolver), pass(this), pointer :: solveMesh=>solveMesh_periodic
		procedure(MeshSolver), pass(this), pointer :: solveMesh_Adj=>solveMesh_Adj_periodic
	contains
		procedure, pass(this) :: buildMesh
		procedure, pass(this) :: setMesh
		procedure, pass(this) :: destroyMesh
	end type

	abstract interface
		subroutine MeshSolver(this,eps)
			use constants
            import mesh
            class(mesh), intent(inout) :: this
            real(mp), intent(in) :: eps
		end subroutine
	end interface

contains

	subroutine buildMesh(this,L,ng,BC)
		class(mesh), intent(out) :: this
		integer, intent(in) :: ng, BC
		real(mp), intent(in) :: L

		this%L = L
		this%ng = ng
		this%BCindex = BC
		select case (this%BCindex)
			case(0)						!periodic
				this%dx = L/ng
                this%solveMesh=>solveMesh_periodic
                this%solveMesh_Adj=>solveMesh_Adj_periodic
			case(1)						!Dirichlet-Dirichlet
				this%dx = L/(ng-1)
                this%solveMesh=>solveMesh_D_D
                this%solveMesh_Adj=>solveMesh_Adj_periodic
			case(2)						!Dirichlet-Neumann
				this%dx = L/(ng-1)
                this%solveMesh=>solveMesh_D_N
                this%solveMesh_Adj=>solveMesh_Adj_periodic
		end select

		allocate(this%E(ng))
		allocate(this%phi(ng))
		allocate(this%rho(ng))
		allocate(this%rho_back(ng))

		this%E = 0.0_mp
		this%phi = 0.0_mp
		this%rho = 0.0_mp
		this%rho_back = 0.0_mp

        allocate(this%W(ng))
		call DSTPoisson_setup(this%ng,this%L,this%W)
	end subroutine

	subroutine setMesh(this,rho_back)
		class(mesh), intent(inout) :: this
		real(mp), intent(in) :: rho_back(this%ng)

		this%rho_back = rho_back
	end subroutine

	subroutine destroyMesh(this)
		class(mesh), intent(inout) :: this

		deallocate(this%E)
		deallocate(this%rho)
		deallocate(this%rho_back)
		deallocate(this%phi)

        deallocate(this%W)
	end subroutine

!===========Mesh Solver===============================================

	subroutine solveMesh_periodic(this,eps)
		class(mesh), intent(inout) :: this
		real(mp), intent(in) :: eps
		real(mp), dimension(this%ng-1) :: rhs, phi1

		rhs = -( this%rho(1:this%ng-1) + this%rho_back(1:this%ng-1) )/eps
		call CG_K(multiplyK,phi1,rhs,this%dx)
		this%phi(1:this%ng-1) = phi1
		this%phi(this%ng) = 0.0_mp
	end subroutine

	subroutine solveMesh_D_D(this,eps)
		class(mesh), intent(inout) :: this
		real(mp), intent(in) :: eps
		real(mp), dimension(this%ng-1) :: rhs, phi1, co1, co2, co3
		co1 = 1.0_mp/this%dx/this%dx
		co2 = -2.0_mp/this%dx/this%dx
		co3 = 1.0_mp/this%dx/this%dx
		co2(this%ng-1) = 1.0_mp
		co1(this%ng-1) = 0.0_mp

		rhs(1:this%ng-2) = -( this%rho(2:this%ng-1) + this%rho_back(2:this%ng-1) )/eps
		!Don't need surface charge: rho_back(ng) will be equal to external voltage
		rhs(this%ng-1) = this%rho_back(this%ng)
		call solve_tridiag(co1,co2,co3,rhs,phi1,this%ng-1)
		this%phi(2:this%ng) = phi1
		this%phi(1) = 0.0_mp
	end subroutine

	subroutine solveMesh_D_N(this,eps)						!D(i=1), N(i=N)
		class(mesh), intent(inout) :: this
		real(mp), intent(in) :: eps
		real(mp), dimension(this%ng-1) :: rhs, phi1, co1, co2, co3
		co1 = 1.0_mp/this%dx/this%dx
		co2 = -2.0_mp/this%dx/this%dx
		co3 = 1.0_mp/this%dx/this%dx
		co2(this%ng-1) = 1.0_mp/this%dx
		co1(this%ng-1) = -1.0_mp/this%dx

		rhs(1:this%ng-2) = -( this%rho(2:this%ng-1) + this%rho_back(2:this%ng-1) )/eps
		rhs(this%ng-1) = ( this%rho_back(this%ng) + 0.5_mp*this%dx*this%rho(this%ng) )/eps
!		call TTA(K_DN,phi1,rhs,this%dx)
		call solve_tridiag(co1,co2,co3,rhs,phi1,this%ng-1)
		this%phi(2:this%ng) = phi1
		this%phi(1) = 0.0_mp
	end subroutine

	!================   Adjoint Mesh solver   ======================================

	subroutine solveMesh_Adj_periodic(m,eps)
		class(mesh), intent(inout) :: m
		real(mp), intent(in) :: eps
		real(mp) :: rhs(m%ng), rho1(m%ng-1)

		rhs = multiplyD(m%E,m%dx,m%BCindex)
		call CG_K(multiplyK,rho1,rhs(1:m%ng-1),m%dx)
		m%rho(1:m%ng-1) = rho1
		m%rho(m%ng) = 0.0_mp
		m%rho = - m%rho/eps
	end subroutine

end module
