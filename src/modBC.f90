module modBC

    use modVelocityProfile
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

	abstract interface
		subroutine applySensBC(p,m,dt,Lv,Nadd,inputVP)
            use modSpecies
            use modMesh
            use modVelocityProfile
			type(species), intent(inout) :: p
			type(mesh), intent(inout) :: m
			real(mp), intent(in) :: dt, Lv
            integer, intent(in) ::  Nadd
            procedure(velocityProfile) :: inputVP
		end subroutine
	end interface

contains

!==============================particle BC=====================================

	subroutine applyBC_periodic(p,m,dt,A0)
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		real(mp), intent(in) :: dt, A0
		integer :: i
        real(mp) :: time1, time2

		call CPU_TIME(time1)
		!apply BC
		do i=1,p%np
			if( p%xp(i)<0 ) then
				p%xp(i) = p%xp(i) + m%L
			elseif( p%xp(i)>=m%L ) then
				p%xp(i) = p%xp(i) - m%L
			end if
		end do
		call CPU_TIME(time2)
		timeProfile(2) = timeProfile(2) + (time2-time1)
        functionCalls(2) = functionCalls(2) + 1
	end subroutine

	subroutine applyBC_absorbing(p,m,dt,A0)
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		real(mp), intent(in) :: dt, A0
		integer :: i, np1
		real(mp), allocatable :: vec(:), vec2(:,:)
        real(mp) :: time1, time2

		call CPU_TIME(time1)
		np1 = p%np
		i = 1
		!apply BC
		do while( i .le. np1 )
			if( p%xp(i).le.0.0_mp ) then
				m%rho_back(1) = m%rho_back(1) + p%spwt(i)*p%qs
				p%xp(i) = p%xp(np1)					!replacement with the last particle
				p%vp(i,:) = p%vp(np1,:)
                p%spwt(i) = p%spwt(np1)
				p%Ep(i) = p%Ep(np1)
				np1 = np1-1
				i = i-1
			elseif( p%xp(i).ge.m%L ) then
				m%rho_back(m%ng) = m%rho_back(m%ng) + p%spwt(i)*p%qs
				p%xp(i) = p%xp(np1)
				p%vp(i,:) = p%vp(np1,:)
                p%spwt(i) = p%spwt(np1)
				p%Ep(i) = p%Ep(np1)
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

		vec = p%spwt(1:np1)
		deallocate(p%spwt)
		allocate(p%spwt(np1))
		p%spwt = vec

		vec = p%Ep(1:np1)
		deallocate(p%Ep)
		allocate(p%Ep(np1))
		p%Ep = vec

		deallocate(vec)
		deallocate(vec2)
		call CPU_TIME(time2)
		timeProfile(2) = timeProfile(2) + (time2-time1)
        functionCalls(2) = functionCalls(2) + 1
	end subroutine

	subroutine applyBC_refluxing_absorbing(p,m,dt,vT)			!refluxing at the left plane, absorbing at the right plane
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		real(mp), intent(in) :: dt, vT
		real(mp) :: temp(3)
		integer :: i, np1
		real(mp), allocatable :: vec(:), vec2(:,:)
        real(mp) :: time1, time2

		call CPU_TIME(time1)
		!apply refluxing BC
		do i=1,p%np
			if( p%xp(i).le.0.0_mp ) then
				temp = (vT*randr(3))
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
				m%rho_back(m%ng) = m%rho_back(m%ng) + p%spwt(i)*p%qs
				p%xp(i) = p%xp(np1)
				p%vp(i,:) = p%vp(np1,:)
                p%spwt(i) = p%spwt(np1)
				p%Ep(i) = p%Ep(np1)
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

		vec = p%spwt(1:np1)
		deallocate(p%spwt)
		allocate(p%spwt(np1))
		p%spwt = vec

		vec = p%Ep(1:np1)
		deallocate(p%Ep)
		allocate(p%Ep(np1))
		p%Ep = vec

		deallocate(vec)
		deallocate(vec2)
		call CPU_TIME(time2)
		timeProfile(2) = timeProfile(2) + (time2-time1)
        functionCalls(2) = functionCalls(2) + 1
	end subroutine

	subroutine applyBC_refluxing_refluxing(p,m,dt,vT)			!refluxing at both planes
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		real(mp), intent(in) :: dt, vT
		real(mp) :: temp(3)
		integer :: i, np1
		real(mp), allocatable :: vec(:), vec2(:,:)
        real(mp) :: time1, time2

		call CPU_TIME(time1)
		!apply refluxing BC
		do i=1,p%np
			if( p%xp(i).le.0.0_mp ) then
				temp = (vT*randr(3))
				temp(1) = abs(temp(1))
				p%vp(i,:) = temp
				call RANDOM_NUMBER(temp)
				p%xp(i) = temp(1)*dt*p%vp(i,1)
			end if
			if( p%xp(i).ge.m%L ) then
				temp = (vT*randr(3))
				temp(1) = -abs(temp(1))
				p%vp(i,:) = temp
				call RANDOM_NUMBER(temp)
				p%xp(i) = m%L + temp(1)*dt*p%vp(i,1)
			end if			
		end do
		call CPU_TIME(time2)
		timeProfile(2) = timeProfile(2) + (time2-time1)
        functionCalls(2) = functionCalls(2) + 1
	end subroutine

!==================== Shock particle BC ====================================

	subroutine applyBC_injecting_reflecting(p,m,dt,vC,vT)			!refluxing at both planes
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		real(mp), intent(in) :: dt, vC, vT
		real(mp) :: temp(3)
		integer :: i, np1
		real(mp), allocatable :: vec(:), vec2(:,:)
        real(mp) :: time1, time2

		call CPU_TIME(time1)
		!apply refluxing BC
		do i=1,p%np
			if( p%xp(i).le.0.0_mp ) then
				temp = vT*randn(3) + vC
				p%vp(i,:) = temp
				call RANDOM_NUMBER(temp)
				p%xp(i) = temp(1)*dt*p%vp(i,1)
			end if
            !reflecting
			if( p%xp(i).ge.m%L ) then
                p%vp(i,:) = -p%vp(i,:)
				p%xp(i) = 2.0_mp*m%L - p%xp(i)
			end if			
		end do
		call CPU_TIME(time2)
		timeProfile(2) = timeProfile(2) + (time2-time1)
        functionCalls(2) = functionCalls(2) + 1
	end subroutine

!=================== sensitivity particle BC ===============================

    subroutine uniformParticleCustomRefluxing(p,m,dt,Lv,Nadd,inputVP)
        type(species), intent(inout) :: p
        type(mesh), intent(inout) :: m
        real(mp), intent(in) :: dt, Lv
        integer, intent(in) :: Nadd
        procedure(velocityProfile) :: inputVP

        integer :: Ninside, NxMax, NaddReal
        real(mp) :: leftFlux, rightFlux
        real(mp), allocatable :: xp(:), vp(:,:), spwt(:)
        logical, dimension(p%np) :: left, right, inside

        integer :: k, ksum, i
        real(mp) :: dxp, dvp

        NxMax = FLOOR(SQRT(2.0_mp*Nadd))
        NaddReal = NxMax*(NxMax+1)/2

        ! Out flux
        left = p%xp.le.0.0_mp
        right = p%xp.ge.m%L
        leftFlux = SUM( PACK(p%spwt, left) )
        rightFlux = SUM( PACK(p%spwt, right) )
        inside = .not.(left.or.right)
        Ninside = COUNT(inside)

        ! allocate new particle array
        allocate(xp(Ninside+2*NaddReal))
        allocate(vp(Ninside+2*NaddReal,3))
        allocate(spwt(Ninside+2*NaddReal))
        vp = 0.0_mp

        ! particle inside domain
        xp(1:Ninside) = PACK( p%xp, inside )
        vp(1:Ninside,1) = PACK( p%vp(:,1), inside )
        spwt(1:Ninside) = PACK( p%spwt, inside )

        dxp = Lv*dt/NxMax
        dvp = Lv/NxMax

        ! influx particles from left boundary
        ksum = Ninside
        do k=1,NxMax
            xp(ksum+1:ksum+k) = (/ ((i-0.5_mp)*dxp,i=1,k) /)
            vp(ksum+1:ksum+k,1) = k*dvp
            ksum = ksum + k
        end do
        spwt(Ninside+1:Ninside+NaddReal) = leftFlux*inputVP(vp(Ninside+1:Ninside+NaddReal,1))         &
                                            /SUM(inputVP(vp(Ninside+1:Ninside+NaddReal,1)))

        ! influx particles from right boundary
        ksum = Ninside+NaddReal
        do k=1,NxMax
            xp(ksum+1:ksum+k) = (/ (m%L-(i-0.5_mp)*dxp,i=1,k) /)
            vp(ksum+1:ksum+k,1) = -k*dvp
            ksum = ksum + k
        end do
        spwt(Ninside+NaddReal+1:Ninside+2*NaddReal) = rightFlux*inputVP(vp(Ninside+NaddReal+1:Ninside+2*NaddReal,1))         &
                                                        /SUM(inputVP(vp(Ninside+NaddReal+1:Ninside+2*NaddReal,1)))

        call p%setSpecies(Ninside+2*NaddReal,xp,vp,spwt)

        ! deallocate particle array
        deallocate(xp)
        deallocate(vp)
        deallocate(spwt)
    end subroutine

    subroutine uniformParticleRefluxingAbsorbing(p,m,dt,vT)
        type(species), intent(inout) :: p
        type(mesh), intent(inout) :: m
        real(mp), intent(in) :: dt, vT

        integer :: Ninside, Nadd, NxMaxL, NxMaxR, NaddReal
        real(mp) :: leftFlux, rightFlux
        real(mp), allocatable :: xp(:), vp(:,:), spwt(:)
        logical, dimension(p%np) :: left, right, inside

        integer :: k, ksum, i
        real(mp) :: dxp, dvp

        ! Out flux
        left = p%xp.le.0.0_mp
        right = p%xp.ge.m%L
        leftFlux = SUM( PACK(p%spwt, left) )
        rightFlux = SUM( PACK(p%spwt, right) )
        inside = .not.(left.or.right)
        Ninside = COUNT(inside)

        Nadd = COUNT(left)
        NxMaxL = FLOOR(SQRT(2.0_mp*Nadd))
        Nadd = COUNT(right)
        NxMaxR = FLOOR(SQRT(2.0_mp*Nadd))
        NaddReal = NxMaxL*(NxMaxL+1)/2 + NxMaxR*(NxMaxR+1)/2

        ! allocate new particle array
        allocate(xp(Ninside+NaddReal))
        allocate(vp(Ninside+NaddReal,3))
        allocate(spwt(Ninside+NaddReal))
        vp = 0.0_mp

        ! particle inside domain
        xp(1:Ninside) = PACK( p%xp, inside )
        vp(1:Ninside,1) = PACK( p%vp(:,1), inside )
        spwt(1:Ninside) = PACK( p%spwt, inside )


        ! influx particles from left boundary
        dxp = 5.0_mp*vT*dt/NxMaxL
        dvp = 5.0_mp*vT/NxMaxL

        ksum = Ninside
        Nadd = NxMaxL*(NxMaxL+1)/2
        do k=1,NxMaxL
            xp(ksum+1:ksum+k) = (/ ((i-0.5_mp)*dxp,i=1,k) /)
            vp(ksum+1:ksum+k,1) = k*dvp
            ksum = ksum + k
        end do
        spwt(Ninside+1:Ninside+Nadd) = leftFlux*EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT )         &
                                            /SUM(EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT ))

        ! influx particles from right boundary
        dxp = 5.0_mp*vT*dt/NxMaxR
        dvp = 5.0_mp*vT/NxMaxR

        ksum = Ninside+Nadd
        do k=1,NxMaxR
            xp(ksum+1:ksum+k) = (/ (m%L-(i-0.5_mp)*dxp,i=1,k) /)
            vp(ksum+1:ksum+k,1) = -k*dvp
            ksum = ksum + k
        end do
        spwt(Ninside+Nadd+1:Ninside+NaddReal) = 0.0_mp

		m%rho_back(m%ng) = m%rho_back(m%ng) + rightFlux*p%qs

        call p%setSpecies(Ninside+NaddReal,xp,vp,spwt)

        ! deallocate particle array
        deallocate(xp)
        deallocate(vp)
        deallocate(spwt)
    end subroutine

    subroutine uniformParticleRefluxingAbsorbing2(p,m,dt,vT)
        type(species), intent(inout) :: p
        type(mesh), intent(inout) :: m
        real(mp), intent(in) :: dt, vT

        integer :: Ninside, Nadd, NxMaxL, NxMaxR, NaddReal
        real(mp) :: leftFlux, rightFlux
        real(mp), allocatable :: xp(:), vp(:,:), spwt(:), xpL(:), xpR(:), vpL(:), vpR(:)
        logical, allocatable :: influxL(:), influxR(:)
        logical, dimension(p%np) :: left, right, inside

        integer :: k, ksum, i
        real(mp) :: dxp, dvp

        ! Out flux
        left = p%xp.le.0.0_mp
        right = p%xp.ge.m%L
        leftFlux = SUM( PACK(p%spwt, left) )
        rightFlux = SUM( PACK(p%spwt, right) )
        inside = .not.(left.or.right)
        Ninside = COUNT(inside)

        Nadd = COUNT(left)
        allocate(xpL(2*Nadd))
        allocate(vpL(2*Nadd))
        allocate(influxL(2*Nadd))
        call RANDOM_NUMBER(xpL)
        call RANDOM_NUMBER(vpL)
        xpL = 5.5_mp*vT*dt*xpL
        vpL = 5.5_mp*vT*vpL
        influxL = (xpL>0.0_mp) .and. (vpL>0.0_mp) .and. ( xpL/vpL.le.dt )

        Nadd = COUNT(right)
        allocate(xpR(2*Nadd))
        allocate(vpR(2*Nadd))
        allocate(influxR(2*Nadd))
        call RANDOM_NUMBER(xpR)
        call RANDOM_NUMBER(vpR)
        xpR = m%L - 5.5_mp*vT*dt*xpR
        vpR = - 5.5_mp*vT*vpR
        influxR = (xpR<m%L) .and. (vpR<0.0_mp) .and. ( (xpR-m%L)/vpR.le.dt )

        NaddReal = COUNT(influxL) + COUNT(influxR)

        ! allocate new particle array
        allocate(xp(Ninside+NaddReal))
        allocate(vp(Ninside+NaddReal,3))
        allocate(spwt(Ninside+NaddReal))
        vp = 0.0_mp

        ! particle inside domain
        xp(1:Ninside) = PACK( p%xp, inside )
        vp(1:Ninside,1) = PACK( p%vp(:,1), inside )
        spwt(1:Ninside) = PACK( p%spwt, inside )


        ! influx particles from left boundary
        Nadd = COUNT(influxL)
        xp(Ninside+1:Ninside+Nadd) = PACK( xpL, influxL )
        vp(Ninside+1:Ninside+Nadd,1) = PACK( vpL, influxL )
        spwt(Ninside+1:Ninside+Nadd) = leftFlux*EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT )         &
                                            /SUM(EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT ))

        ! influx particles from right boundary
        xp(Ninside+Nadd+1:Ninside+NaddReal) = PACK( xpR, influxR )
        vp(Ninside+Nadd+1:Ninside+NaddReal,1) = PACK( vpR, influxR )
        spwt(Ninside+Nadd+1:Ninside+NaddReal) = 0.0_mp

		m%rho_back(m%ng) = m%rho_back(m%ng) + rightFlux*p%qs

        call p%setSpecies(Ninside+NaddReal,xp,vp,spwt)

        ! deallocate particle array
        deallocate(xpL)
        deallocate(vpL)
        deallocate(xpR)
        deallocate(vpR)
        deallocate(influxL)
        deallocate(influxR)
        deallocate(xp)
        deallocate(vp)
        deallocate(spwt)
    end subroutine

    subroutine uniformParticleRefluxingRefluxing(p,m,dt,vT)
        type(species), intent(inout) :: p
        type(mesh), intent(inout) :: m
        real(mp), intent(in) :: dt, vT

        integer :: Ninside, Nadd, NxMaxL, NxMaxR, NaddReal
        real(mp) :: leftFlux, rightFlux
        real(mp), allocatable :: xp(:), vp(:,:), spwt(:), xpL(:), xpR(:), vpL(:), vpR(:)
        logical, allocatable :: influxL(:), influxR(:)
        logical, dimension(p%np) :: left, right, inside

        integer :: k, ksum, i
        real(mp) :: dxp, dvp

        ! Out flux
        left = p%xp.le.0.0_mp
        right = p%xp.ge.m%L
        leftFlux = SUM( PACK(p%spwt, left) )
        rightFlux = SUM( PACK(p%spwt, right) )
        inside = .not.(left.or.right)
        Ninside = COUNT(inside)

        Nadd = COUNT(left)
        allocate(xpL(2*Nadd))
        allocate(vpL(2*Nadd))
        allocate(influxL(2*Nadd))
        call RANDOM_NUMBER(xpL)
        call RANDOM_NUMBER(vpL)
        xpL = 5.5_mp*vT*dt*xpL
        vpL = 5.5_mp*vT*vpL
        influxL = (xpL>0.0_mp) .and. (vpL>0.0_mp) .and. ( xpL/vpL.le.dt )

        Nadd = COUNT(right)
        allocate(xpR(2*Nadd))
        allocate(vpR(2*Nadd))
        allocate(influxR(2*Nadd))
        call RANDOM_NUMBER(xpR)
        call RANDOM_NUMBER(vpR)
        xpR = m%L - 5.5_mp*vT*dt*xpR
        vpR = - 5.5_mp*vT*vpR
        influxR = (xpR<m%L) .and. (vpR<0.0_mp) .and. ( (xpR-m%L)/vpR.le.dt )

        NaddReal = COUNT(influxL) + COUNT(influxR)

        ! allocate new particle array
        allocate(xp(Ninside+NaddReal))
        allocate(vp(Ninside+NaddReal,3))
        allocate(spwt(Ninside+NaddReal))
        vp = 0.0_mp

        ! particle inside domain
        xp(1:Ninside) = PACK( p%xp, inside )
        vp(1:Ninside,1) = PACK( p%vp(:,1), inside )
        spwt(1:Ninside) = PACK( p%spwt, inside )


        ! influx particles from left boundary
        Nadd = COUNT(influxL)
        xp(Ninside+1:Ninside+Nadd) = PACK( xpL, influxL )
        vp(Ninside+1:Ninside+Nadd,1) = PACK( vpL, influxL )
        spwt(Ninside+1:Ninside+Nadd) = leftFlux*EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT )         &
                                            /SUM(EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT ))

        ! influx particles from right boundary
        xp(Ninside+Nadd+1:Ninside+NaddReal) = PACK( xpR, influxR )
        vp(Ninside+Nadd+1:Ninside+NaddReal,1) = PACK( vpR, influxR )
        spwt(Ninside+Nadd+1:Ninside+NaddReal) = 0.0_mp
        spwt(Ninside+Nadd+1:Ninside+NaddReal) = rightFlux*EXP( -vp(Ninside+Nadd+1:Ninside+NaddReal,1)**2     &
                                            /2.0_mp/vT/vT )/SUM(EXP( -vp(Ninside+Nadd+1:Ninside+NaddReal,1)**2        &
                                            /2.0_mp/vT/vT ))

        call p%setSpecies(Ninside+NaddReal,xp,vp,spwt)

        ! deallocate particle array
        deallocate(xpL)
        deallocate(vpL)
        deallocate(xpR)
        deallocate(vpR)
        deallocate(influxL)
        deallocate(influxR)
        deallocate(xp)
        deallocate(vp)
        deallocate(spwt)
    end subroutine

    subroutine ionBoundarySensitivityToTau(p,m,dt,vT,mu,tau,ionFluxL,sensitivityFluxR)
        type(species), intent(inout) :: p
        type(mesh), intent(inout) :: m
        real(mp), intent(in) :: dt, vT, mu, tau, ionFluxL
        real(mp), intent(out) :: sensitivityFluxR
        

        integer :: Ninside, Nadd, NxMaxL, NxMaxR, NaddReal
        real(mp) :: leftFlux, rightFlux
        real(mp) :: Z1, Z2
        real(mp), allocatable :: xp(:), vp(:,:), spwt(:), xpL(:), xpR(:), vpL(:), vpR(:)
        logical, allocatable :: influxL(:), influxR(:)
        logical, dimension(p%np) :: left, right, inside

        integer :: k, ksum, i
        real(mp) :: dxp, dvp

        ! Out flux
        left = p%xp.le.0.0_mp
        right = p%xp.ge.m%L
        leftFlux = SUM( PACK(p%spwt, left) )
        rightFlux = SUM( PACK(p%spwt, right) )
        inside = .not.(left.or.right)
        Ninside = COUNT(inside)

        Nadd = COUNT(left)
        allocate(xpL(2*Nadd))
        allocate(vpL(2*Nadd))
        allocate(influxL(2*Nadd))
        call RANDOM_NUMBER(xpL)
        call RANDOM_NUMBER(vpL)
        xpL = 5.5_mp*vT*dt*xpL
        vpL = 5.5_mp*vT*vpL
        influxL = (xpL>0.0_mp) .and. (vpL>0.0_mp) .and. ( xpL/vpL.le.dt )

        Nadd = COUNT(right)
        allocate(xpR(2*Nadd))
        allocate(vpR(2*Nadd))
        allocate(influxR(2*Nadd))
        call RANDOM_NUMBER(xpR)
        call RANDOM_NUMBER(vpR)
        xpR = m%L - 5.5_mp*vT*dt*xpR
        vpR = - 5.5_mp*vT*vpR
        influxR = (xpR<m%L) .and. (vpR<0.0_mp) .and. ( (xpR-m%L)/vpR.le.dt )

        NaddReal = COUNT(influxL) + COUNT(influxR)

        ! allocate new particle array
        allocate(xp(Ninside+NaddReal))
        allocate(vp(Ninside+NaddReal,3))
        allocate(spwt(Ninside+NaddReal))
        vp = 0.0_mp

        ! particle inside domain
        xp(1:Ninside) = PACK( p%xp, inside )
        vp(1:Ninside,1) = PACK( p%vp(:,1), inside )
        spwt(1:Ninside) = PACK( p%spwt, inside )


        ! influx particles from left boundary
        Nadd = COUNT(influxL)
        xp(Ninside+1:Ninside+Nadd) = PACK( xpL, influxL )
        vp(Ninside+1:Ninside+Nadd,1) = PACK( vpL, influxL )
        Z1 = SUM(EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT ))
        Z2 = SUM(vp(Ninside+1:Ninside+Nadd,1)**2*EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT ))
        spwt(Ninside+1:Ninside+Nadd) = leftFlux*EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT )         &
                                            /SUM(EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT ))       &
                                     + ionFluxL/tau*EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT )      &
                                              *( 1.0_mp/Z1 - vp(Ninside+1:Ninside+Nadd,1)**2/Z2 )

        ! influx particles from right boundary
        xp(Ninside+Nadd+1:Ninside+NaddReal) = PACK( xpR, influxR )
        vp(Ninside+Nadd+1:Ninside+NaddReal,1) = PACK( vpR, influxR )
        spwt(Ninside+Nadd+1:Ninside+NaddReal) = 0.0_mp

		m%rho_back(m%ng) = m%rho_back(m%ng) + rightFlux*p%qs
        sensitivityFluxR = rightFlux

        call p%setSpecies(Ninside+NaddReal,xp,vp,spwt)

        ! deallocate particle array
        deallocate(xpL)
        deallocate(vpL)
        deallocate(xpR)
        deallocate(vpR)
        deallocate(influxL)
        deallocate(influxR)
        deallocate(xp)
        deallocate(vp)
        deallocate(spwt)
    end subroutine

    subroutine testIonBoundarySensitivity(p,m,dt,vT,mu,tau,ionFlux)
        type(species), intent(inout) :: p
        type(mesh), intent(inout) :: m
        real(mp), intent(in) :: dt, vT, mu, tau, ionFlux

        integer :: Ninside, Nadd, NxMaxL, NxMaxR, NaddReal
        real(mp) :: leftFlux, rightFlux
        real(mp) :: Z1, Z2
        real(mp), allocatable :: xp(:), vp(:,:), spwt(:), xpL(:), xpR(:), vpL(:), vpR(:)
        logical, allocatable :: influxL(:), influxR(:)
        logical, dimension(p%np) :: left, right, inside

        integer :: k, ksum, i
        real(mp) :: dxp, dvp

        ! Out flux
        left = p%xp.le.0.0_mp
        right = p%xp.ge.m%L
        leftFlux = SUM( PACK(p%spwt, left) )
        rightFlux = SUM( PACK(p%spwt, right) )
        inside = .not.(left.or.right)
        Ninside = COUNT(inside)

        Nadd = COUNT(left)
        allocate(xpL(2*Nadd))
        allocate(vpL(2*Nadd))
        allocate(influxL(2*Nadd))
        call RANDOM_NUMBER(xpL)
        call RANDOM_NUMBER(vpL)
        xpL = 5.5_mp*vT*dt*xpL
        vpL = 5.5_mp*vT*vpL
        influxL = (xpL>0.0_mp) .and. (vpL>0.0_mp) .and. ( xpL/vpL.le.dt )

        Nadd = COUNT(right)
        allocate(xpR(2*Nadd))
        allocate(vpR(2*Nadd))
        allocate(influxR(2*Nadd))
        call RANDOM_NUMBER(xpR)
        call RANDOM_NUMBER(vpR)
        xpR = m%L - 5.5_mp*vT*dt*xpR
        vpR = - 5.5_mp*vT*vpR
        influxR = (xpR<m%L) .and. (vpR<0.0_mp) .and. ( (xpR-m%L)/vpR.le.dt )

        NaddReal = COUNT(influxL) + COUNT(influxR)

        ! allocate new particle array
        allocate(xp(Ninside+NaddReal))
        allocate(vp(Ninside+NaddReal,3))
        allocate(spwt(Ninside+NaddReal))
        vp = 0.0_mp

        ! particle inside domain
        xp(1:Ninside) = PACK( p%xp, inside )
        vp(1:Ninside,1) = PACK( p%vp(:,1), inside )
        spwt(1:Ninside) = PACK( p%spwt, inside )


        ! influx particles from left boundary
        Nadd = COUNT(influxL)
        xp(Ninside+1:Ninside+Nadd) = PACK( xpL, influxL )
        vp(Ninside+1:Ninside+Nadd,1) = PACK( vpL, influxL )
        Z1 = SUM(EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT ))
        Z2 = SUM(vp(Ninside+1:Ninside+Nadd,1)**2*EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT ))
        spwt(Ninside+1:Ninside+Nadd) = leftFlux*EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT )         &
                                            /SUM(EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT ))       &
                                     + ionFlux/tau*EXP( -vp(Ninside+1:Ninside+Nadd,1)**2/2.0_mp/vT/vT )      &
                                              *( 1.0_mp/Z1 - vp(Ninside+1:Ninside+Nadd,1)**2/Z2 )

        ! influx particles from right boundary
        xp(Ninside+Nadd+1:Ninside+NaddReal) = PACK( xpR, influxR )
        vp(Ninside+Nadd+1:Ninside+NaddReal,1) = PACK( vpR, influxR )
        Z1 = SUM(EXP( -vp(Ninside+Nadd+1:Ninside+NaddReal,1)**2/2.0_mp/vT/vT ))
        Z2 = SUM(vp(Ninside+Nadd+1:Ninside+NaddReal,1)**2*EXP( -vp(Ninside+Nadd+1:Ninside+NaddReal,1)**2/2.0_mp/vT/vT ))
        spwt(Ninside+Nadd+1:Ninside+NaddReal) = rightFlux*EXP( -vp(Ninside+Nadd+1:Ninside+NaddReal,1)**2/2.0_mp/vT/vT )/Z1      &
                                     + ionFlux/tau*EXP( -vp(Ninside+Nadd+1:Ninside+NaddReal,1)**2/2.0_mp/vT/vT )      &
                                              *( 1.0_mp/Z1 - vp(Ninside+Nadd+1:Ninside+NaddReal,1)**2/Z2 )

        call p%setSpecies(Ninside+NaddReal,xp,vp,spwt)

        ! deallocate particle array
        deallocate(xpL)
        deallocate(vpL)
        deallocate(xpR)
        deallocate(vpR)
        deallocate(influxL)
        deallocate(influxR)
        deallocate(xp)
        deallocate(vp)
        deallocate(spwt)
    end subroutine

end module
