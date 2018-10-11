module modRecord

	use modPM1D
    use modMPI

	implicit none

	type recordData
		integer :: nt, n, ng, mod
		real(mp) :: L
		character(len=:), allocatable :: dir

		integer, allocatable :: np(:,:)
		real(mp), allocatable :: phidata(:,:)
		real(mp), allocatable :: Edata(:,:)
		real(mp), allocatable :: rhodata(:,:)
		real(mp), allocatable :: PE(:), KE(:,:)
		integer, allocatable :: n_coll(:,:)               !(species*collision type, time)

		real(mp), allocatable :: cpt_time(:,:)
	contains
		procedure, pass(this) :: buildRecord
		procedure, pass(this) :: destroyRecord
		procedure, pass(this) :: recordPlasma
		procedure, pass(this) :: printPlasma
	end type

contains

	subroutine buildRecord(this,nt,n,L,ng,input_dir,mod)
		class(recordData), intent(out) :: this
		integer, intent(in) :: nt, n, ng, mod
		real(mp), intent(in) :: L
		character(len=*), intent(in), optional :: input_dir
		integer :: nr
		nr = nt/mod+1

		this%nt = nt
		this%n = n
		this%L = L
		this%ng = ng
		this%mod = mod

		allocate(this%np(n,nr))
		allocate(this%phidata(ng,nr))
		allocate(this%Edata(ng,nr))
		allocate(this%rhodata(ng,nr))
		allocate(this%PE(nr))
		allocate(this%KE(n,nr))

		this%np = 0
		this%phidata = 0.0_mp
		this%Edata = 0.0_mp
		this%rhodata = 0.0_mp
		this%PE = 0.0_mp
		this%KE = 0.0_mp

		if( present(input_dir) ) then
			allocate(character(len=len(input_dir)) :: this%dir)
			this%dir = input_dir
		else
			allocate(character(len=0) :: this%dir)
			this%dir = ''
		end if
        if( print_pm_output ) then
    		print *, 'Directory: ','data/'//trim(this%dir)
        end if

		call system('mkdir -p data/'//this%dir//'/xp')
		call system('mkdir -p data/'//this%dir//'/vp')
		call system('mkdir -p data/'//this%dir//'/spwt')

		call system('rm data/'//this%dir//'/xp/*.*')
		call system('rm data/'//this%dir//'/vp/*.*')
		call system('rm data/'//this%dir//'/spwt/*.*')
	end subroutine

	subroutine destroyRecord(this)
		class(recordData), intent(inout) :: this
		integer :: i,j

		deallocate(this%np)
		deallocate(this%phidata)
		deallocate(this%Edata)
		deallocate(this%rhodata)

		deallocate(this%PE)
		deallocate(this%KE)

		if( allocated(this%n_coll) ) then
			deallocate(this%n_coll)
		end if

		deallocate(this%dir)
	end subroutine

	subroutine recordPlasma(this,pm,k)
		class(recordData), intent(inout) :: this
		class(PM1D), intent(in) :: pm
		integer, intent(in) :: k					!k : time step
		integer :: n,j, kr, fileUnit								!n : species
		character(len=100) :: nstr, kstr
		real(mp) :: qe = 1.602176565E-19
        fileUnit = mpih%my_rank+305

		if( (this%mod.eq.1) .or. (mod(k,this%mod).eq.0) ) then
			kr = merge(k,k/this%mod,this%mod.eq.1)
			do n=1,pm%n
				write(nstr,*) n
				write(kstr,*) kr
				open(unit=fileUnit,file='data/'//this%dir//'/xp/'//trim(adjustl(kstr))//'_'	&
					//trim(adjustl(nstr))//'.bin',status='replace',form='unformatted',access='stream')
				write(fileUnit) pm%p(n)%xp
				close(fileUnit)
				open(unit=fileUnit,file='data/'//this%dir//'/vp/'//trim(adjustl(kstr))//'_'	&
					//trim(adjustl(nstr))//'.bin',status='replace',form='unformatted',access='stream')
				write(fileUnit) pm%p(n)%vp
				close(fileUnit)
				open(unit=fileUnit,file='data/'//this%dir//'/spwt/'//trim(adjustl(kstr))//'_'	&
					//trim(adjustl(nstr))//'.bin',status='replace',form='unformatted',access='stream')
				write(fileUnit) pm%p(n)%spwt
				close(fileUnit)

				!time step: 0~Nt, in array: 1~(Nt+1) (valgrind prefers this way of allocation)
				this%np(n,kr+1) = pm%p(n)%np
				this%KE(n,kr+1) = 0.5_mp*pm%p(n)%ms*SUM( pm%p(n)%spwt*pm%p(n)%vp(:,1)**2 )
			end do
			this%phidata(:,kr+1) = pm%m%phi
			this%Edata(:,kr+1) = pm%m%E
			this%rhodata(:,kr+1) = pm%m%rho
			this%PE(kr+1) = 0.5_mp*SUM(pm%m%E**2)*pm%m%dx
!			this%cpt_time(:,kr+1) = this%cpt_temp

            if( print_pm_output ) then
    			print *, '============= ',k,'-th Time Step ================='
    			do n=1,pm%n
    				print *, 'Species(',n,'): ',pm%p(n)%np, ', KE: ', this%KE(n,kr+1),'J'
    			end do
                print *, 'Voltage = ',pm%m%phi(pm%ng),'V'
            end if
		end if
	end subroutine

	subroutine printPlasma(this)
		class(recordData), intent(in) :: this
		character(len=100) :: s
		integer :: i,j,fileUnit
		real(mp) :: total, mean, pct
        fileUnit = mpih%my_rank+305

		open(unit=fileUnit,file='data/'//this%dir//'/record',status='replace')
        print *, this%n, this%ng, this%nt, this%L, this%mod
		write(fileUnit,*) this%n, this%ng, this%nt, this%L, this%mod
		close(fileUnit)

		open(unit=fileUnit,file='data/'//this%dir//'/Ncoll.bin',status='replace',form='unformatted',access='stream')
        write(fileUnit) this%n_coll
        close(fileUnit)

		open(unit=fileUnit,file='data/'//this%dir//'/cpt_time.bin',status='replace',form='unformatted',access='stream')
		write(fileUnit) this%cpt_time
		close(fileUnit)

		open(unit=fileUnit,file='data/'//this%dir//'/E.bin',status='replace',form='unformatted',access='stream')
		do i = 1,this%nt/this%mod+1
			write(fileUnit) this%Edata(:,i)
		end do
		close(fileUnit)

		open(unit=fileUnit,file='data/'//this%dir//'/rho.bin',status='replace',form='unformatted',access='stream')
		do i = 1,this%nt/this%mod+1
			write(fileUnit) this%rhodata(:,i)
		end do
		close(fileUnit)

		open(unit=fileUnit,file='data/'//this%dir//'/PE.bin',status='replace',form='unformatted',access='stream')
		do i = 1,this%nt/this%mod+1
			write(fileUnit) this%PE(i)
		end do
		close(fileUnit)

		open(unit=fileUnit,file='data/'//this%dir//'/Np.bin',status='replace',form='unformatted',access='stream')
		do i = 1,this%nt/this%mod+1
			write(fileUnit) this%np(:,i)
		end do
		close(fileUnit)

		open(unit=fileUnit,file='data/'//this%dir//'/phi.bin',status='replace',form='unformatted',access='stream')
		do i = 1,this%nt/this%mod+1
			write(fileUnit) this%phidata(:,i)
		end do
		close(fileUnit)

		do i=1,this%n
			write(s,*) i
			open(unit=fileUnit,file='data/'//this%dir//'/KE_'//trim(adjustl(s))//'.bin',                &
                 status='replace',form='unformatted',access='stream')
    		do j = 1,this%nt/this%mod+1
    	        write(fileUnit) this%KE(i,j)
    		end do
		    close(fileUnit)
		end do

701	FORMAT	(A, E10.3,'	', I6,'	', F10.3,'%')
		if( SUM(functionCalls(8:13)).eq.0 ) then
            total = SUM(timeProfile)
			open(unit=fileUnit,file='data/'//this%dir//'/original_cpt_summary.dat',status='replace')
			write(fileUnit,*) 'Step	Total	Mean	Percentage'
			print *, "================ Computation Time Summary ==================================="
			print *, "Original simulation	   	     Total        Calls	 Percentage	"

			pct = timeProfile(1)/total*100.0_mp
			print 701, "Particle Move			", timeProfile(1), functionCalls(1), pct
			write(fileUnit,701) 'Particle-Move	', timeProfile(1), functionCalls(1), pct

			pct = timeProfile(2)/total*100.0_mp
			print 701, "ApplyBC     			", timeProfile(2), functionCalls(2), pct
			write(fileUnit,701) 'ApplyBC    	', timeProfile(2), functionCalls(2), pct

			pct = timeProfile(3)/total*100.0_mp
			print 701, "ChargeAssign			", timeProfile(3), functionCalls(3), pct
			write(fileUnit,701) 'Charge-Assign	', timeProfile(3), functionCalls(3), pct

			pct = timeProfile(4)/total*100.0_mp
			print 701, "Poisson Solver			", timeProfile(4), functionCalls(4), pct
			write(fileUnit,701) 'Poisson-Solver	', timeProfile(4), functionCalls(4), pct

			pct = timeProfile(5)/total*100.0_mp
			print 701, "Efield Gradient			", timeProfile(5), functionCalls(5), pct
			write(fileUnit,701) 'Efield-Gradient	', timeProfile(5), functionCalls(5), pct

			pct = timeProfile(6)/total*100.0_mp
			print 701, "Force Assign			", timeProfile(6), functionCalls(6), pct
			write(fileUnit,701) 'Force-Assign	', timeProfile(6), functionCalls(6), pct

			pct = timeProfile(7)/total*100.0_mp
			print 701, "Particle Accel			", timeProfile(7), functionCalls(7), pct
			write(fileUnit,701) 'Particle-Accel	', timeProfile(7), functionCalls(7), pct
			print *, "============================================================================="
			close(fileUnit)
		else
            total = SUM(timeProfile)

			open(unit=fileUnit,file='data/'//this%dir//'/sensitivity_cpt_summary.dat',status='replace')
			write(fileUnit,*) 'Step	Total	Mean	Percentage'
			print *, "================ Computation Time Summary ==================================="
			print *, "Sensitivity simulation	  	     Total        Calls   	 Percentage	"

			pct = timeProfile(1)/total*100.0_mp
			print 701, "Particle Move			", timeProfile(1), functionCalls(1), pct
			write(fileUnit,701) 'Particle-Move	', timeProfile(1), functionCalls(1), pct

			pct = timeProfile(2)/total*100.0_mp
			print 701, "ApplyBC     			", timeProfile(2), functionCalls(2), pct
			write(fileUnit,701) 'ApplyBC    	', timeProfile(2), functionCalls(2), pct

			pct = timeProfile(3)/total*100.0_mp
			print 701, "ChargeAssign			", timeProfile(3), functionCalls(3), pct
			write(fileUnit,701) 'Charge-Assign	', timeProfile(3), functionCalls(3), pct

			pct = timeProfile(4)/total*100.0_mp
			print 701, "Poisson Solver			", timeProfile(4), functionCalls(4), pct
			write(fileUnit,701) 'Poisson-Solver	', timeProfile(4), functionCalls(4), pct

			pct = timeProfile(5)/total*100.0_mp
			print 701, "Efield Gradient			", timeProfile(5), functionCalls(5), pct
			write(fileUnit,701) 'Efield-Gradient	', timeProfile(5), functionCalls(5), pct

			pct = timeProfile(6)/total*100.0_mp
			print 701, "Force Assign			", timeProfile(6), functionCalls(6), pct
			write(fileUnit,701) 'Force-Assign	', timeProfile(6), functionCalls(6), pct

			pct = timeProfile(7)/total*100.0_mp
			print 701, "Particle Accel			", timeProfile(7), functionCalls(7), pct
			write(fileUnit,701) 'Particle-Accel	', timeProfile(7), functionCalls(7), pct

			pct = timeProfile(8)/total*100.0_mp
			print 701, "f Velocity Gradient		", timeProfile(8), functionCalls(8), pct
			write(fileUnit,701) 'GradVf	', timeProfile(8), functionCalls(8), pct

			pct = timeProfile(9)/total*100.0_mp
			print 701, "Source Evaluation		", timeProfile(9), functionCalls(9), pct
			write(fileUnit,701) 'Source-Evaluation	', timeProfile(9), functionCalls(9), pct

			pct = timeProfile(10)/total*100.0_mp
			print 701, "Number Density			", timeProfile(10), functionCalls(10), pct
			write(fileUnit,701) 'Number-Density	', timeProfile(10), functionCalls(10), pct

			pct = timeProfile(11)/total*100.0_mp
			print 701, "HoverN/SourceAssign	", timeProfile(11), functionCalls(11), pct
			write(fileUnit,701) 'HoverN-SourceAssign	', timeProfile(11), functionCalls(11), pct

			pct = timeProfile(12)/total*100.0_mp
			print 701, "SourceAssign/Addition   		", timeProfile(12), functionCalls(12), pct
			write(fileUnit,701) 'SourceAssign-Addition	', timeProfile(12), functionCalls(12), pct

			pct = timeProfile(13)/total*100.0_mp
			print 701, "Redistribution			", timeProfile(13), functionCalls(13), pct
			write(fileUnit,701) 'Redistribution	', timeProfile(13), functionCalls(13), pct

			print *, "============================================================================="
			close(fileUnit)

            print *, "Factor: ", SUM(timeProfile)/SUM(timeProfile(1:7))
		end if
	end subroutine

!=======================================================
!	For no_collision: don't compute n_coll
!=======================================================

	subroutine set_null_discharge(r)
		type(recordData), intent(inout), optional :: r

		if( present(r) ) then
			allocate(r%n_coll(1,r%nt))
			r%n_coll = 0
		end if
	end subroutine

end module
