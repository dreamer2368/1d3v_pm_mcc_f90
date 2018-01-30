module MCCProblems

	use init
	use timeStep
	use modMPI

	implicit none

contains

!===========   MCC collision setting   ==============================

	subroutine set_Ar_discharge(pm, A, r)
		type(PM1D), intent(inout) :: pm
      type(recordData), intent(inout), optional :: r
		real(mp), intent(in) :: A(4)						!A(1): temperature of neutral(eV),	A(2): density of neutral(m-3),
                                                   !A(3): discharge current density(A/m2), A(4): discharge frequency(Hz)
		if( pm%n .ne. 2 ) then
			print *, 'ERROR : the number of species should be two corresponding to electon and Argon+. stopped the simulation.'
			stop
		end if

		!Electron species
		call buildSpecies(pm%p(1),-q_e,m_e)
		!Argon cation species
		call buildSpecies(pm%p(2),q_e,m_Ar)

		deallocate(pm%A0)
		allocate(pm%A0(4))
		pm%A0 = A

      if( present(r) ) then
         allocate(r%n_coll(7,r%nt))
      end if

		pm%mcc_collision=>Argon_Electron
	end subroutine

   subroutine Ar_discharge
      type(PM1D) :: pm
      type(recordData) :: r
      real(mp), parameter :: Kb = 1.38065E-23, EV_TO_K = 11604.52_mp, eps = 8.85418782E-12, mTorr_to_Pa = 0.13332237_mp
      real(mp) :: I0 = 25.6_mp, I_f = 13.56E6                !I0(A/m2), If(Hz)
      real(mp) :: TN = 0.026_mp, PN = 50.0_mp, gden        !TN(eV), PN(mTorr), gden(m-3)
	  real(mp), parameter :: n0=10.0_mp**15!, v0_e = 5.9E5, v0_Ar = 2.19E3
	  real(mp) :: T0=1.0_mp, wp0, lambda0, v0_e, v0_Ar
	  real(mp), parameter :: L = 0.02_mp, area = 0.016_mp               !L(m), area(m2)
	  integer, parameter :: nc2p = 10**6, Np=CEILING(n0*L*area/nc2p), Ng = 300
	  real(mp) :: spwt(Np), xp0(Np), vp0(Np,3)
	  real(mp) :: dt
	  integer :: i
	  gden = (PN*mTorr_to_Pa)/(q_e*TN)

	  spwt = n0*L/Np                      !spwt = nc2p/area

	  print *, 'gden(m-3): ',gden,', n0(m-3): ',n0,', spwt(m-2): ',spwt(1)

!      T0 = 0.5_mp*m_e*v0_e**2/q_e*EV_TO_K
	  v0_e = sqrt(2.0_mp*T0/EV_TO_K*q_e/m_e)
	  v0_Ar = sqrt(2.0_mp*T0/EV_TO_K*q_e/m_Ar)
	  wp0 = sqrt(n0*q_e*q_e/m_e/eps)
	  lambda0 = sqrt(eps*T0/n0/q_e)

	  print *, 'L = ',L,', lambda0 = ',lambda0,' e = lambda/L = ',lambda0/L

	  dt = 7.20179E-11
	  print *, 'dt = ',dt,', wp*dt=',wp0*dt
	  call init_random_seed
	  call null_collision(gden,dt)
	  print *, 'P_e = ',col_prob_e,', P_Ar = ', col_prob_Ar

	  call buildPM1D(pm,40000.0_mp*dt,100.0_mp*dt,Ng,2,pBC=1,mBC=2,order=1,dt=dt,L=L,eps=eps)
	  call buildRecord(r,pm%nt,2,pm%L,pm%ng,'rf_Ar4',20)
	  open(unit=301,file='data/rf_Ar4/input',status='replace',form='unformatted',access='stream')
	  write(301) TN, PN, v0_e, v0_Ar, I0, L, area, dt
	  close(301)

	  call set_Ar_discharge(pm,(/TN,gden,I0,I_f/),r)
	  call RANDOM_NUMBER(xp0)
	  vp0 = randn(Np,3)*v0_e
	  call setSpecies(pm%p(1),Np,xp0*L,vp0,spwt)
	  call RANDOM_NUMBER(xp0)
	  vp0 = randn(Np,3)*v0_Ar
	  call setSpecies(pm%p(2),Np,xp0*L,vp0,spwt)

      call forwardsweep(pm,r,RF_current,Null_source)

      call printPlasma(r)

      call destroyPM1D(pm)
      call destroyRecord(r)
   end subroutine

	subroutine cross_section
		integer, parameter :: N=10000
		real(mp), dimension(N) :: energy, sig1, sig2, sig3, sig4, sig5
		integer :: i

		energy = exp( log(10.0_mp)*( (/ (i,i=1,N) /)/(0.2_mp*N) - 2.0_mp ) )
		do i=1,N
			sig1(i) = asigma1(energy(i))
			sig2(i) = asigma2(energy(i))
			sig3(i) = asigma3(energy(i))
			sig4(i) = asigma4(energy(i))
			sig5(i) = asigma5(energy(i))
		end do

		call system('mkdir -p data/cross_section')
		open(unit=301,file='data/cross_section/sig1.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/cross_section/sig2.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/cross_section/sig3.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/cross_section/sig4.bin',status='replace',form='unformatted',access='stream')
		open(unit=305,file='data/cross_section/sig5.bin',status='replace',form='unformatted',access='stream')
		open(unit=306,file='data/cross_section/energy.bin',status='replace',form='unformatted',access='stream')
		write(301) sig1
		write(302) sig2
		write(303) sig3
		write(304) sig4
		write(305) sig5
		write(306) energy
		close(301)
		close(302)
		close(303)
		close(304)
		close(305)
		close(306)
	end subroutine

end module
