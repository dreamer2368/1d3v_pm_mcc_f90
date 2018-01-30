module ArMCC

	use modSpecies
	use random
	implicit none

	real(mp), parameter :: q_e = 1.602E-19, m_e = 9.10938356e-31, m_Ar = 6.68E-26!6.6335209e-26			!kg
	real(mp), parameter :: max_sigmav_e = 6.6038e-13, max_sigmav_Ar = 2.8497e-14				!m3/s
	real(mp), parameter :: extengy0 = 11.55_mp, ionengy0 = 15.76_mp									!eV (excitation, ionization)
	real(mp) :: col_prob_e = 0.0_mp, col_prob_Ar = 0.0_mp

	abstract interface
		subroutine mcc_collision(p,A0,n_coll)
			use modSpecies
			class(species), intent(inout) :: p(:)
			real(mp), intent(in) :: A0(:)
			integer, intent(out), optional :: n_coll(7)
		end subroutine
	end interface

contains

!=======================================================
!	MCC global subroutine interface: mcc_collision
!=======================================================

	!No collision --- Use for collisionless PIC
	subroutine no_collision(p,A0,n_coll)
		class(species), intent(inout) :: p(:)
		real(mp), intent(in) :: A0(:)
		integer, intent(out), optional :: n_coll(7)

	end subroutine

	!Argon-Electron collision
	subroutine Argon_Electron(p,A0,n_coll)
		class(species), intent(inout) :: p(:)
		real(mp), intent(in) :: A0(:)
      integer, intent(out), optional :: n_coll(7)

		if( present(n_coll) ) then
			call mcc_Argon(p(2),A0(1),n_coll(5:7))
			call mcc_electron(p(1),p(2),A0(1),n_coll(1:4))
		else
			call mcc_Argon(p(2),A0(1))
			call mcc_electron(p(1),p(2),A0(1))
		end if
	end subroutine


!=======================================================
!	Compute null collision for each species (e, Ar+)
!=======================================================
	subroutine null_collision(gden, dt)
		real(mp), intent(in) :: gden, dt

		col_prob_e = 1.0_mp - exp( -max_sigmav_e*gden*dt )
		col_prob_Ar = 1.0_mp - exp( -max_sigmav_Ar*gden*dt )
	end subroutine

!=======================================================
!	MCC for each species
!=======================================================
	subroutine mcc_electron(Electron,Ar,T_eV,n_diag)						!electron species
		type(species), intent(inout) :: Electron, Ar
		real(mp), intent(in) :: T_eV									!Temperature in eV
		integer, intent(out), optional :: n_diag(4)
		integer :: n_coll, nnp, idx
		real(mp) :: rnd, rnd_ion, temp_x, temp_v(3)
		real(mp) :: engy, rengy, vel, nu_total_vel, sum_sigma(3)
		integer :: new_e, new_Ar
		real(mp), allocatable :: vec1_e(:), vec2_e(:,:), vec1_Ar(:), vec2_Ar(:,:), temp1(:), temp2(:,:)
		real(mp) :: vT
		integer :: i
		integer :: n_elastic, n_excite, n_ionize
      if( present(n_diag) ) then
   		n_coll = 0
	   	n_elastic=0
		   n_excite=0
		   n_ionize=0
      end if

		!Pick particles for collisions
		n_coll = floor( Electron%np*col_prob_e )
		nnp = Electron%np
		do i=1,n_coll
			call RANDOM_NUMBER(rnd)
			idx = ceiling( nnp*rnd )
			if( idx>nnp ) then
				idx = nnp
			end if

			!sort them into the tail of arrays
			temp_x = Electron%xp(nnp)
			temp_v = Electron%vp(nnp,:)
			Electron%xp(nnp) = Electron%xp(idx)
			Electron%vp(nnp,:) = Electron%vp(idx,:)
			Electron%xp(idx) = temp_x
			Electron%vp(idx,:) = temp_v
			nnp = nnp-1
		end do

		!Pick the type of collision for each particle
		!Note: we don't consider the velocity of the neutral, assuming that it is much smaller than that of electron species
		vT = sqrt( T_eV*q_e/m_Ar )
		new_e = 0
		new_Ar = 0
		allocate(vec1_e(n_coll))
		allocate(vec2_e(n_coll,3))
		allocate(vec1_Ar(n_coll))
		allocate(vec2_Ar(n_coll,3))
		do i=nnp+1,nnp+n_coll
			vel = sqrt( sum( Electron%vp(i,:)**2 ) )
			engy = 0.5_mp*Electron%ms*vel**2/q_e			!scale in eV
			nu_total_vel = max_sigmav_e/vel

			call RANDOM_NUMBER(rnd)
			sum_sigma(1) = asigma1(engy)
			sum_sigma(2) = asigma1(engy) + asigma2(engy)
			sum_sigma(3) = asigma1(engy) + asigma2(engy) + asigma3(engy)
			!Elastic
			if( rnd .le. sum_sigma(1)/nu_total_vel ) then

				call anewvel_e(engy,m_e,m_Ar,Electron%vp(i,:),.true.)
            if( present(n_diag) ) then
   				n_elastic = n_elastic+1
            end if

			!Excitation
			elseif( (engy.ge.extengy0) .and. (rnd.le.sum_sigma(2)/nu_total_vel) ) then

				engy = engy - extengy0
				Electron%vp(i,:) = Electron%vp(i,:)/vel
				vel = sqrt( 2.0_mp/Electron%ms*q_e*engy )
				Electron%vp(i,:) = Electron%vp(i,:)*vel
				call anewvel_e(engy,m_e,m_Ar,Electron%vp(i,:),.false.)
            if( present(n_diag) ) then
   				n_excite = n_excite+1
            end if

			!Ionization
			elseif( (engy.ge.ionengy0) .and. (rnd.le.sum_sigma(3)/nu_total_vel) ) then

				!subtract ionization energy
				!and partition the energy between created and scattered electron
				engy = engy-ionengy0
				call RANDOM_NUMBER(rnd_ion)
				rengy = 10.0_mp*tan( rnd_ion*atan(engy/20.0_mp) )
				engy = engy - rengy

				!scatter the created electron
				new_e = new_e+1
				vec1_e(new_e) = Electron%xp(i)
				vec2_e(new_e,:) = Electron%vp(i,:)/vel*sqrt( 2.0_mp/Electron%ms*q_e*rengy )
				call anewvel_e(rengy,m_e,m_Ar,vec2_e(new_e,:),.false.)

				!assign velocity to the created Ar ion
				new_Ar = new_Ar+1
				vec1_Ar(new_Ar) = Electron%xp(i)
				vec2_Ar(new_Ar,:) = vT*randn(3)

				!scatter the incident electron
				Electron%vp(i,:) = Electron%vp(i,:)/vel*sqrt( 2.0_mp/Electron%ms*q_e*engy )
				call anewvel_e(engy,m_e,m_Ar,Electron%vp(i,:),.false.)
            if( present(n_diag) ) then
   				n_ionize = n_ionize+1
            end if
			end if
		end do
      if( present(n_diag) ) then
   		n_diag = (/ n_coll, n_elastic, n_excite, n_ionize /)
      end if

		!Add newly created particles
      !Electron species
		Electron%np = Electron%np + new_e
      !xp
		allocate(temp1(Electron%np))
		temp1(1:Electron%np-new_e) = Electron%xp
		temp1(Electron%np-new_e+1:Electron%np) = vec1_e(1:new_e)
		deallocate(Electron%xp)
		allocate(Electron%xp(Electron%np))
		Electron%xp = temp1
      !spwt
		temp1(1:Electron%np-new_e) = Electron%spwt
		temp1(Electron%np-new_e+1:Electron%np) = Electron%spwt(1)				!Assume uniform spwt!!
		deallocate(Electron%spwt)
		allocate(Electron%spwt(Electron%np))
		Electron%spwt = temp1
      !Ep
		temp1(1:Electron%np-new_e) = Electron%Ep
		temp1(Electron%np-new_e+1:Electron%np) = 0.0_mp
		deallocate(Electron%Ep)
		allocate(Electron%Ep(Electron%np))
		Electron%Ep = temp1
		deallocate(temp1)
      !vp
		allocate(temp2(Electron%np,3))
		temp2(1:Electron%np-new_e,:) = Electron%vp
		temp2(Electron%np-new_e+1:Electron%np,:) = vec2_e(1:new_e,:)
		deallocate(Electron%vp)
		allocate(Electron%vp(Electron%np,3))
		Electron%vp = temp2
		deallocate(temp2)
      !Argon species
		Ar%np = Ar%np + new_Ar
      !xp
		allocate(temp1(Ar%np))
		temp1(1:Ar%np-new_Ar) = Ar%xp
		temp1(Ar%np-new_Ar+1:Ar%np) = vec1_Ar(1:new_Ar)
		deallocate(Ar%xp)
		allocate(Ar%xp(Ar%np))
		Ar%xp = temp1
      !spwt
		temp1(1:Ar%np-new_Ar) = Ar%spwt
		temp1(Ar%np-new_Ar+1:Ar%np) = Ar%spwt(1)						!Assume uniform spwt!!!
		deallocate(Ar%spwt)
		allocate(Ar%spwt(Ar%np))
		Ar%spwt = temp1
      !Ep
		temp1(1:Ar%np-new_Ar) = Ar%Ep
		temp1(Ar%np-new_Ar+1:Ar%np) = 0.0_mp
		deallocate(Ar%Ep)
		allocate(Ar%Ep(Ar%np))
		Ar%Ep = temp1
		deallocate(temp1)
      !vp
		allocate(temp2(Ar%np,3))
		temp2(1:Ar%np-new_Ar,:) = Ar%vp
		temp2(Ar%np-new_Ar+1:Ar%np,:) = vec2_Ar(1:new_Ar,:)
		deallocate(Ar%vp)
		allocate(Ar%vp(Ar%np,3))
		Ar%vp = temp2
		deallocate(temp2)
      !deallocate temporary arrays
		deallocate(vec1_e)
		deallocate(vec2_e)
		deallocate(vec1_Ar)
		deallocate(vec2_Ar)
	end subroutine

	subroutine mcc_Argon(Ar,T_eV,N_diag)
		type(species), intent(inout) :: Ar
		real(mp), intent(in) :: T_eV								!Temperature in eV
      integer, intent(out), optional :: N_diag(3)
		integer :: n_coll, nnp, idx
		real(mp) :: rnd, temp_x, temp_v(3), vT, vn(3)
		real(mp) :: vel, engy, nu_total_vel, sum_sigma(2)
		integer :: i
		integer :: n_elastic, n_exchange
      if( present(N_diag) ) then
         n_coll = 0
         n_elastic = 0
         n_exchange = 0
      end if

		!Pick particles for collisions
		n_coll = floor( Ar%np*col_prob_Ar )
		nnp = Ar%np
		do i=1,n_coll
			call RANDOM_NUMBER(rnd)
			idx = ceiling( nnp*rnd )
			if( idx>nnp ) then
				idx = nnp
			end if

			temp_x = Ar%xp(nnp)
			temp_v = Ar%vp(nnp,:)
			Ar%xp(nnp) = Ar%xp(idx)
			Ar%vp(nnp,:) = Ar%vp(idx,:)
			Ar%xp(idx) = temp_x
			Ar%vp(idx,:) = temp_v
			nnp = nnp-1
		end do

		!Pick the type of collision for each particle
		!Note: we consider the velocity of the neutral here
		vT = sqrt( T_eV*q_e/m_Ar )
		do i=nnp+1,nnp+n_coll
			vn = vT*randn(3)
			Ar%vp(i,:) = Ar%vp(i,:) - vn
			vel = sqrt( sum( Ar%vp(i,:)**2 ) )
			engy = 0.5_mp*Ar%ms*vel**2/q_e			!scale in eV
			nu_total_vel = max_sigmav_Ar/vel

			call RANDOM_NUMBER(rnd)
			sum_sigma(1) = asigma4(engy)
			sum_sigma(2) = asigma4(engy) + asigma5(engy)
			!Charge Exchange
			if( rnd .le. sum_sigma(1)/nu_total_vel ) then

				Ar%vp(i,:) = 0.0_mp
            if( present(N_diag) ) then
   				n_exchange = n_exchange+1
            end if

			!Elastic scattering
			elseif( rnd .le. sum_sigma(2)/nu_total_vel ) then

				call anewvel_Ar(Ar%vp(i,:))
            if( present(N_diag) ) then
   				n_elastic = n_elastic+1
            end if

			end if
			Ar%vp(i,:) = Ar%vp(i,:) + vn
		end do

      if( present(N_diag) ) then
		   N_diag = (/ n_coll, n_exchange, n_elastic /)
      end if
	end subroutine

!=======================================================
!  Ar+ + Ar Differential cross section(isotropic, hard-sphere): determine scattered velocity
!=======================================================
	subroutine anewvel_Ar(vp)
		real(mp), intent(inout) :: vp(3)
		real(mp) :: random, coschi, sinchi, phi1, cosphi, sinphi, up(3)
		real(mp) :: vel, r1(3), r2(3), r3(3)

		call RANDOM_NUMBER(random)
		coschi = sqrt(random)
		sinchi = sqrt( abs(1.0_mp - coschi**2) )

		call RANDOM_NUMBER(random)
		phi1 = 2.0_mp*pi*random
		cosphi = cos(phi1)
		sinphi = sin(phi1)

		vel = sqrt( sum( vp**2 ) )
		r3 = vp/vel

		up = 0.0_mp
		if( r3(3) .eq. 1.0_mp ) then
			up(2) = 1.0_mp
		else
			up(3) = 1.0_mp
		end if

		r2(1) = r3(2)*up(3) - r3(3)*up(2)
		r2(2) = r3(3)*up(1) - r3(1)*up(3)
		r2(3) = r3(1)*up(2) - r3(2)*up(1)
		r2 = r2/sqrt( sum( r2**2 ) )

		r1(1) = r2(2)*r3(3) - r2(3)*r3(2)
		r1(2) = r2(3)*r3(1) - r2(1)*r3(3)
		r1(3) = r2(1)*r3(2) - r2(2)*r3(1)

		vp = vel*coschi*(r1*sinchi*cosphi + r2*sinchi*sinphi + r3*coschi)
	end subroutine

!=======================================================
!  e + Ar Differential cross section: determine scattered velocity
!=======================================================
	subroutine anewvel_e(energy, m1, m2, vp, elastic_flag)					!m1: electron mass, m2: Argon mass
		real(mp), intent(in) :: energy, m1, m2
		real(mp), intent(inout) :: vp(3)
		logical, intent(in) :: elastic_flag
		real(mp) :: random
		real(mp) :: coschi, sinchi, phi1, cosphi, sinphi, up(3)
		real(mp) :: vel
		real(mp) :: r1(3), r2(3), r3(3)

		if( energy < 1E-30 ) then
			coschi = 1.0_mp
		else
			call RANDOM_NUMBER(random)
			coschi = ( 2.0_mp + energy - 2.0_mp*( 1.0_mp+energy )**random )/energy
		end if
		sinchi = sqrt( abs( 1.0_mp - coschi**2 ) )

		call RANDOM_NUMBER(random)
		phi1 = 2.0_mp*pi*random
		cosphi = cos(phi1)
		sinphi = sin(phi1)

		vel = sqrt( sum( vp**2 ) )
		r3 = vp/vel

		if( elastic_flag ) then
			vel = vel*sqrt( 1.0_mp - 2.0_mp*m1/m2*(1.0_mp-coschi) )
		end if

		up = 0.0_mp
		if( r3(3) .eq. 1.0_mp ) then
			up(2) = 1.0_mp
		else
			up(3) = 1.0_mp
		end if

		r2(1) = r3(2)*up(3) - r3(3)*up(2)
		r2(2) = r3(3)*up(1) - r3(1)*up(3)
		r2(3) = r3(1)*up(2) - r3(2)*up(1)
		r2 = r2/sqrt( sum( r2**2 ) )

		r1(1) = r2(2)*r3(3) - r2(3)*r3(2)
		r1(2) = r2(3)*r3(1) - r2(1)*r3(3)
		r1(3) = r2(1)*r3(2) - r2(2)*r3(1)

		vp = vel*(r1*sinchi*cosphi + r2*sinchi*sinphi + r3*coschi)
	end subroutine

!=======================================================
!  e + Ar -> e + Ar  Elastic      (with Ramseur minimum)
!=======================================================
	function asigma1(energy) result(sig1)
		real(mp), intent(in) :: energy
		real(mp) :: sig1

		if( energy < 0.2_mp ) then
			sig1 = 1.0_mp/( 10.0_mp**(19.0_mp+energy/0.11_mp) )
		else
			sig1 = 9.070000E-19*(energy**1.55_mp)*( (energy+70.0_mp)**1.10_mp )/( (14.0_mp+energy)**3.25_mp )
		end if
	end function

!=======================================================
!  e + Ar -> e + Ar  Excitation
!=======================================================
	function asigma2(energy) result(sig2)
		real(mp), intent(in) :: energy
		real(mp) :: sig2

		if( energy > 12.0_mp ) then
			sig2 = ( 3.85116E-19*log(energy/3.4015_mp) - 4.85227E-19 )/energy
		else
			sig2 = 0.0_mp
		end if
	end function

!=======================================================
!  e + Ar -> e + e + Ar+  Ion.
!=======================================================
	function asigma3(energy) result(sig3)
		real(mp), intent(in) :: energy
		real(mp) :: sig3

		if( energy > 15.76_mp ) then
			sig3 = 1.3596E-18/energy*log( (energy+120.0_mp/energy)/15.76_mp )*	&
					( atan( (energy**2 - 9.76_mp*energy + 2.4_mp)/(20.6_mp*energy + 206.0_mp) ) +	&
						atan( (2.0_mp*energy - 80.0_mp)/(10.3_mp*energy+103.0_mp) ) )
		else
			sig3 = 0.0_mp
		end if
	end function

!=======================================================
!  Ar + Ar+ -> Ar+ + Ar  Charge X
!=======================================================
	function asigma4(energy) result(sig4)
		real(mp), intent(in) :: energy
		real(mp) :: sig4

		if( energy > 4.0_mp ) then
			sig4 = 2.0E-19 + 5.5E-19/sqrt(energy)
		else
			sig4 = -2.95E-19*sqrt(energy) + 10.65E-19
		end if
	end function

!=======================================================
!  Ar + Ar+ -> Ar + Ar+   Scat.
!=======================================================
	function asigma5(energy) result(sig5)
		real(mp), intent(in) :: energy
		real(mp) :: sig5

		if( energy > 4.0_mp ) then
			sig5 = 1.8E-19 + 4.0E-19/sqrt(energy)
		else
			sig5 = -2.0E-19*sqrt(energy) + 7.8E-19
		end if
	end function

end module
