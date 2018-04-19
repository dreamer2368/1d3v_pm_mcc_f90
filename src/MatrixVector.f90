module MatrixVector

	use constants

	implicit none

	include 'fftw3.f'

contains

	function multiplyD(x,dx,idx) result(y)						!Derivative with BC
		real(mp), intent(in) :: x(:)
		real(mp), intent(in) :: dx
		integer, intent(in) :: idx								!Boundary index
		real(mp) :: y(size(x))
		integer :: i

		y=0.0_mp
		do i=2,size(x)-1
			y(i) = 0.5_mp/dx*( x(i+1) - x(i-1) )
		end do

		select case(idx)
			case(0)												!periodic BC
				y(1) = 0.5_mp/dx*( x(2) - x(size(x)) )
				y(size(x)) = 0.5_mp/dx*( x(1) - x(size(x)-1) )
			case(1)												!D_D BC
				y(1) = 0.5_mp/dx*( -3.0_mp*x(1) + 4.0_mp*x(2) - x(3) )
				y(size(x)) = 0.5_mp/dx*( x(size(x)-2) - 4.0_mp*x(size(x)-1) + 3.0_mp*x(size(x)) )
			case(2)												!D_N BC
				y(1) = 0.5_mp/dx*( -3.0_mp*x(1) + 4.0_mp*x(2) - x(3) )
				y(size(x)) = 0.5_mp/dx*( x(size(x)-2) - 4.0_mp*x(size(x)-1) + 3.0_mp*x(size(x)) )
		end select
	end function

!===================== Matrices =====================================================

	function multiplyK(x,dx) result(y)		!periodic BC + Dirichlet BC
		real(mp), intent(in) :: x(:)
		real(mp), intent(in) :: dx
		real(mp) :: y(size(x))
		integer :: i

		y=0.0_mp
		do i=2,size(x)-1
			y(i) = 1.0_mp/dx/dx*( x(i+1) - 2.0_mp*x(i) + x(i-1) )
		end do
		y(1) = 1.0_mp/dx/dx*( x(2) - 2.0_mp*x(1) )
		y(size(x)) = 1.0_mp/dx/dx*( - 2.0_mp*x(size(x)) + x(size(x)-1) )
	end function

	function K_DN(i,N,dx) result(co)			!coefficients for 2nd-order K with Dirichlet(i=0) + Neumann(i=N) BC
		integer, intent(in) :: i, N
		real(mp), intent(in) :: dx
		real(mp) :: co(3)						!co = (/ a, b, c /)

		co = 0.0_mp
		if( i.eq.1 ) then
			co(2) = -2.0_mp/dx/dx
			co(3) = 1.0_mp/dx/dx
		elseif( i.eq.N ) then
			co(1) = -1.0_mp/dx
			co(2) = 1.0_mp/dx
		else
			co(1) = 1.0_mp/dx/dx
			co(2) = -2.0_mp/dx/dx
			co(3) = 1.0_mp/dx/dx
		end if
	end function

!===================== Solver =====================================================

	subroutine CG_K(K,x,b,dx)								!Kx = b
		real(mp), intent(in) :: b(:)
		real(mp), intent(in) :: dx
		real(mp), intent(out) :: x(:)
		real(mp) :: r(size(x)), p(size(x)), r1(size(x))
		real(mp) :: alpha, beta
		real(mp) :: tol
		integer :: i = 0
		interface
			function K(x,dx) result(y)
				use constants
				real(mp), intent(in) :: x(:)
				real(mp), intent(in) :: dx
				real(mp) :: y(size(x))
			end function
		end interface

		select case (mp)
			case (SELECTED_REAL_KIND(4))
				tol = 1.0e-16
			case (SELECTED_REAL_KIND(15))
				tol = 1.0e-32
			case (SELECTED_REAL_KIND(33))
				tol = 10.0_mp**(-64.0_mp)
			case default
				tol = 1.0e-32
		end select

		if( size(b)/=size(x) ) then
			print *, '===================================================='
			print *, '====================  FAULT  ======================='
			print *, '=========  x AND b HAVE NOT EQUAL SIZES  ==========='
			print *, '====================  Ax=b  ========================'
			print *, '===================================================='
		end if

		x = 0.0_mp
		r = b - K(x,dx)
		p = r

		i = 0
		do while( DOT_PRODUCT(r,r) > tol )
			alpha = DOT_PRODUCT(r,r)/DOT_PRODUCT(p,K(p,dx))
			x = x + alpha*p
			r1 = r
			r = r - alpha*K(p,dx)

			beta = DOT_PRODUCT(r,r)/DOT_PRODUCT(r1,r1)
			p = r + beta*p
			i = i+1
			if( i > 1e8 ) then
				print *, '====================================='
				print *, '=========  CG METHOD FAILS   ========'
				print *, '====================================='
				exit
			end if
		end do
	end subroutine

	subroutine TTA(K,x,b,dx)												!Thomas Tridiagonal Algorithm
		real(mp), intent(in) :: b(:)
		real(mp), intent(in) :: dx
		real(mp), intent(out) :: x(:)
		real(mp), dimension(size(b)-1) :: w, g
		real(mp) :: co(3)
		integer :: i, N
		interface
			function K(i,N,dx) result(y)
				use constants
				integer, intent(in) :: i,N
				real(mp), intent(in) :: dx
				real(mp) :: y(3)
			end function
		end interface
		N = size(b)

		co = K(1,N,dx)
		w(1) = co(3)/co(2)
		g(1) = b(1)/co(2)
		do i=2,N-1
			co = K(i,N,dx)
			w(i) = co(3)/(co(2) - co(1)*w(i-1))
			g(i) = (b(i) - co(1)*g(i-1))/(co(2) - co(1)*w(i-1))
		end do
		co = K(N,N,dx)
		x(N) = (b(N) - co(1)*g(N-1))/(co(2) - co(1)*w(N-1))

		do i=1,N-1
			x(N-i) = g(N-i) - w(N-i)*x(N-i+1)
		end do
	end subroutine

	subroutine solve_tridiag(a,b,c,d,x,n)
		!	 a - sub-diagonal (means it is the diagonal below the main diagonal)
		!	 b - the main diagonal
		!	 c - sup-diagonal (means it is the diagonal above the main diagonal)
		!	 d - right part
		!	 x - the answer
		!	 n - number of equations
		integer,intent(in) :: n
		real(mp),dimension(n),intent(in) :: a,b,c,d
		real(mp),dimension(n),intent(out) :: x
		real(mp),dimension(n) :: cp,dp
		real(mp) :: m
		integer i

		! initialize c-prime and d-prime
		cp(1) = c(1)/b(1)
		dp(1) = d(1)/b(1)
		! solve for vectors c-prime and d-prime
		do i = 2,n
			m = b(i)-cp(i-1)*a(i)
			cp(i) = c(i)/m
			dp(i) = (d(i)-dp(i-1)*a(i))/m
		enddo
		! initialize x
		x(n) = dp(n)
		! solve for x from the vectors c-prime and d-prime
		do i = n-1, 1, -1
			x(i) = dp(i)-cp(i)*x(i+1)
		end do
	end subroutine

	subroutine DSTPoisson_setup(N,L,W)
		integer, intent(in) :: N
		real(mp), intent(in) :: L
		complex(mp), intent(out) :: W(N)
		integer :: i,wi
		complex(mp) :: wx

		wx = 2.0_mp*pi*eye/L
        do i=0,N-1
            if( i.le.N/2 ) then
                wi = i 
            else
                wi = - ( N-i )
            end if
            W(i+1) = + (wx*wi)**2.0_mp
        end do
        W(1) = 1.0_mp
	end subroutine

	subroutine DSTPoisson(x,rhs,W)
		real(mp), intent(in) :: rhs(:)
		complex(mp), intent(in) :: W(:)
		real(mp), intent(out) :: x(size(rhs))
		complex(mp), dimension(size(rhs)) :: rhsFFT, xFFT, rhsb, xb
		integer(mp) :: plan
		integer :: N
		N = size(rhs)

        rhsb = rhs
		call dfftw_plan_dft_1d(plan,N,rhsb,rhsFFT,FFTW_FORWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,rhsb,rhsFFT)
		call dfftw_destroy_plan(plan)

		xFFT = rhsFFT/W
        xFFT(1) = (0.0_mp,0.0_mp)

		call dfftw_plan_dft_1d(plan,N,xFFT,xb,FFTW_BACKWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,xFFT,xb)
		call dfftw_destroy_plan(plan)

		x = REALPART(xb)*1.0_mp/N
	end subroutine

end module
