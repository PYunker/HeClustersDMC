program test
	use dmc_utils
	implicit none
    
	integer, parameter :: DP = kind(0.0D0)
	integer, parameter :: n_blocks = 10, n_steps = 1000, N0 = 10000, d = 1
	integer :: walker_count, i_blok, i_krok, m, j, step_count
	real(kind=DP), parameter :: x0(d) = 0, dt = 1e-3, dt2 = sqrt(dt), alfa = 0.1d0
	real(kind=DP) :: Er, avg_pe, pe, u(d), rep(d), v
	
	type(walker(d)), pointer :: root, last, current, previous

	write(*,"(1x,'Blok',4x,'N',5x,'<H>',6x,'<V>')")
	
	walker_count = N0

	! walker initialization
	allocate(root)
	root%pos = -1
	root%next => null()
	last => root
	do j = 1,N0
		! spawn all the walkers at zero
		call append(x0, last, d)
	end do
	! obtain initial energy
	call getEnergy(x0,Er)

	do i_blok = 1,n_blocks
		avg_pe = 0
		do i_krok = 1,n_steps
			previous => root
			current => root%next
			step_count = walker_count
			do j = 1,step_count
				! obtain roughly normal random 'u'				
				call random_number(u)
				u = 0.61477929423D0 * log(u/(1-u))
				rep = current%pos + dt2 * u
				call getEnergy(rep, pe)
				! perform the branching
				call random_number(v)
				m = floor(v + exp((Er - pe)*dt))
				if (m == 0) then
					walker_count = walker_count - 1
					call remove(previous, current, d)
					cycle
				else if (m == 2) then
					walker_count = walker_count + 1
					call append(rep, last, d)
				else if (m /= 1) then
					walker_count = walker_count + 2
					call append(rep, last, d)
					call append(rep, last, d)
				end if
				current%pos = rep
				avg_pe = avg_pe + pe
				previous => current
				current => current%next
			end do
			! at the end of every step
			avg_pe = avg_pe / walker_count
			Er = Er + alfa * (1 - real(walker_count,8) / N0) ! bylo k mist N0
		end do
		! at th end of every block
		avg_pe = avg_pe / n_steps
		write(*,"(I3,3x,I7,f9.5,f9.5)") i_blok, walker_count, Er, avg_pe
	end do

	! output into file
	current => root%next
	open(1, file='populace_dmc', status='replace')
	write(1,*) '# name: populace'
	write(1,*) '# type: matrix'
	write(1,*) '# rows: ', d
	write(1,*) '# columns: ', walker_count

	do j = 1,walker_count
		write(1,*) current%pos
		current => current%next
	end do

	print*, 'walker_count : ', walker_count
	print*, 'i : ', j
	
	! manual deallocation to please the almighty leak sanitizer
	! =========================================================
!	previous => root                                         !|
!	current => root%next                                     !|
!	deallocate(last)                                         !|
!	! do while(associated(current))                          !|
!	do j=1,walker_count-1                                    !|
!		call remove(previous, current, d)                    !|
!		current => previous%next                             !|
!	end do                                                   !|
!	deallocate(root)                                         !|
	! =========================================================

contains

	! get energy of a walker
	pure subroutine getEnergy(x,energy)
		real(kind=DP), intent(in) :: x(d)
		real(kind=DP), intent(out) :: energy
		energy = 0.5 * norm2(x)**2
	end subroutine getEnergy

end program test
