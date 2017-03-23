module modBC

	use modSpecies
	use modMesh
	use random

	implicit none

	abstract interface
		subroutine applyBC(p,m,dt,A0)
			use modSpecies
			use modMesh
			type(species), intent(inout) :: p
			type(mesh), intent(inout) :: m
			real(mp), intent(in) :: dt, A0
		end subroutine
	end interface

contains

!==============================particle BC=====================================

	subroutine applyBC_periodic(p,m,dt,A0)
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		real(mp), intent(in) :: dt, A0
		integer :: i

		!apply BC
		do i=1,p%np
			if( p%xp(i)<0 ) then
				p%xp(i) = p%xp(i) + m%L
			elseif( p%xp(i)>=m%L ) then
				p%xp(i) = p%xp(i) - m%L
			end if
		end do
	end subroutine

	subroutine applyBC_absorbing(p,m,dt,A0)
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		real(mp), intent(in) :: dt, A0
		integer :: i, np1
		real(mp), allocatable :: vec(:), vec2(:,:)

		np1 = p%np
		i = 1
		!apply BC
		do while( i .le. np1 )
			if( p%xp(i).le.0.0_mp ) then
				p%xp(i) = p%xp(np1)					!replacement with the last particle
				p%vp(i,:) = p%vp(np1,:)
				p%Ep(i) = p%Ep(np1)
				m%rho_back(1) = m%rho_back(1) + p%spwt(i)*p%qs
				np1 = np1-1
				i = i-1
			elseif( p%xp(i).ge.m%L ) then
				p%xp(i) = p%xp(np1)
				p%vp(i,:) = p%vp(np1,:)
				p%Ep(i) = p%Ep(np1)
				m%rho_back(m%ng) = m%rho_back(m%ng) + p%spwt(i)*p%qs
				np1 = np1-1
				i = i-1
			end if
			i = i+1
		end do
		p%np = np1

		allocate(vec(np1))
		allocate(vec2(np1,3))

		vec = p%xp(1:np1)
		deallocate(p%xp)
		allocate(p%xp(np1))
		p%xp = vec

		vec2 = p%vp(1:np1,:)
		deallocate(p%vp)
		allocate(p%vp(np1,3))
		p%vp = vec2

		vec = p%Ep(1:np1)
		deallocate(p%Ep)
		allocate(p%Ep(np1))
		p%Ep = vec

		deallocate(vec)
		deallocate(vec2)
	end subroutine

	subroutine applyBC_refluxing_absorbing(p,m,dt,vT)			!refluxing at the left plane, absorbing at the right plane
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		real(mp), intent(in) :: dt, vT
		real(mp) :: temp(3)
		integer :: i, np1
		real(mp), allocatable :: vec(:), vec2(:,:)

		!apply refluxing BC
		do i=1,p%np
			if( p%xp(i).le.0.0_mp ) then
				temp = (vT*randn(3))
				temp(1) = abs(temp(1))
				p%vp(i,:) = temp
				call RANDOM_NUMBER(temp)
				p%xp(i) = temp(1)*dt*p%vp(i,1)
			end if
		end do

		np1 = p%np
		i = 1
		!apply absorbing BC
		do while( i .le. np1 )
			if( p%xp(i).ge.m%L ) then
				p%xp(i) = p%xp(np1)
				p%vp(i,:) = p%vp(np1,:)
				p%Ep(i) = p%Ep(np1)
				m%rho_back(m%ng) = m%rho_back(m%ng) + p%spwt(i)*p%qs
				np1 = np1-1
				i = i-1
			end if
			i = i+1
		end do
		p%np = np1

		allocate(vec(np1))
		allocate(vec2(np1,3))

		vec = p%xp(1:np1)
		deallocate(p%xp)
		allocate(p%xp(np1))
		p%xp = vec

		vec2 = p%vp(1:np1,:)
		deallocate(p%vp)
		allocate(p%vp(np1,3))
		p%vp = vec2

		vec = p%Ep(1:np1)
		deallocate(p%Ep)
		allocate(p%Ep(np1))
		p%Ep = vec

		deallocate(vec)
		deallocate(vec2)
	end subroutine

	subroutine applyBC_refluxing_refluxing(p,m,dt,vT)			!refluxing at both planes
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		real(mp), intent(in) :: dt, vT
		real(mp) :: temp(3)
		integer :: i, np1
		real(mp), allocatable :: vec(:), vec2(:,:)

		!apply refluxing BC
		do i=1,p%np
			if( p%xp(i).le.0.0_mp ) then
				temp = (vT*randn(3))
				temp(1) = abs(temp(1))
				p%vp(i,:) = temp
				call RANDOM_NUMBER(temp)
				p%xp(i) = temp(1)*dt*p%vp(i,1)
			end if
			if( p%xp(i).ge.m%L ) then
				temp = (vT*randn(3))
				temp(1) = -abs(temp(1))
				p%vp(i,:) = temp
				call RANDOM_NUMBER(temp)
				p%xp(i) = m%L + temp(1)*dt*p%vp(i,1)
			end if			
		end do
	end subroutine

end module
