module ArMCC

	use modPM1D
	use random
	implicit none

	real(mp), parameter :: max_sig_e = 1.6022E-19, max_sig_Ar = 1.7955E-18

contains

!=======================================================
!  Ar+ + Ar Differential cross section(isotropic): determine scattered velocity
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
	subroutine anewvel_e(energy, m1, m2, vp, elastic_flag)
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