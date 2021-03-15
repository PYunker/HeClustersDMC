program akde_mse_test
	use metropolis
	use grid_utils
	use kde
	use distributions
	use mpi_f08
	use omp_lib
	! use iso_fortran_env, only : int64
	implicit none
	integer, parameter :: DP = KIND(0.0D0), DP_SIZE = sizeof(0.0D0), &
	                      INT_SIZE = sizeof(1), d = 1, grid_n = 150
	integer :: world_size, world_rank, run, node_runs, runs, gam, sample_size
	real(kind=DP) :: Ef(grid_n) = 0, Ef2(grid_n) = 0, estimate(grid_n), &
	                 Es(grid_n) = 0, Es2(grid_n) = 0
	real(kind=DP) :: grid(d,grid_n), h, lam
	real(kind=DP), allocatable :: sample(:,:), mu(:,:), w(:), sig(:,:,:), det(:)
	procedure(rn_to_r), pointer :: pdf => comb
	character(len=10) :: arg
	! variables concerning MPI IO
	integer :: offset_runs
	! includes : grid_n, Ef, Ef2, Es, Es2, grid, exact (when available)
	! includes : d, gam, runs, mu, w, sig, det		
	type(MPI_FILE) :: estimator_file
	type(MPI_STATUS) :: estimator_file_stat
	integer :: estimate_file_unit, mu_offset, w_offset, sig_offset, det_offset

	call MPI_INIT()
	
	! loads the value of gam (kernel count), sample_size, lam (scaling constant),
	! and grid_n respectively, they must all be present when invoking the executable
	call load_command_line_args()

	! set up communicators and balance the load between procesess
	call mpi_setup()

	! initialize grid    
	!call grid_init_2d(grid, -5.0D0, 5.0D0, -5.0D0, 5.0D0, 50, 50)
	call grid_init_1d(grid,-5.0D0,5.0D0,grid_n)

	call omp_set_num_threads(4)	
	
	! main cycle - for every run, draw new sample and obtain estimate
	do run=1,node_runs
		call draw_sample(sample, sample_size, d, pdf)
        ! obtain the estimate
		call botev_get_estimate(mu, w, sig, det, d, gam, sample_size, 5.0d0, &
								sample, grid, grid_n, runs, estimate, h)
		Ef = Ef + estimate / runs
		Ef2 = Ef2 + estimate**2 / runs

		call kde_get_estimate(d, sample_size, sample, grid, grid_n, estimate)

		Es = Es + estimate / runs
		Es2 = Es2 + estimate**2 / runs

		! output estimator variables into file
		call append_estimator()		

	end do
    
	! reuse 'run' & 'estimate' to avoid having to declare more one-time
	! variables
	!$OMP PARALLEL DO	
	do run=1,grid_n
		estimate(run) = pdf(grid(:,run))
	end do
	!$OMP END PARALLEL DO
	! normalize (this isn't ideal way to go about it)	
	estimate = estimate / sum(estimate) * sum(Ef)

	if (world_rank == 0) then
		call MPI_REDUCE(MPI_IN_PLACE, Ef, grid_n, MPI_REAL8, &
		                MPI_SUM, 0, MPI_COMM_WORLD)
		call MPI_REDUCE(MPI_IN_PLACE, Ef2, grid_n, MPI_REAL8, &
		                MPI_SUM, 0, MPI_COMM_WORLD)
	else
		call MPI_REDUCE(Ef, Ef, grid_n, MPI_REAL8, &
		                MPI_SUM, 0, MPI_COMM_WORLD)
		call MPI_REDUCE(Ef2, Ef2, grid_n, MPI_REAL8, &
		                MPI_SUM, 0, MPI_COMM_WORLD)
	end if

	if (world_rank == 0) then
		call MPI_REDUCE(MPI_IN_PLACE, Es, grid_n, MPI_REAL8, &
		                MPI_SUM, 0, MPI_COMM_WORLD)
		call MPI_REDUCE(MPI_IN_PLACE, Es2, grid_n, MPI_REAL8, &
		                MPI_SUM, 0, MPI_COMM_WORLD)
	else
		call MPI_REDUCE(Es, Es, grid_n, MPI_REAL8, &
		                MPI_SUM, 0, MPI_COMM_WORLD)
		call MPI_REDUCE(Es2, Es2, grid_n, MPI_REAL8, &
		                MPI_SUM, 0, MPI_COMM_WORLD)
	end if
	
	! output Ef, Ef2, Es, Es2, grid, exact (if available) into file
	if (world_rank == 0) then	
		call save_estimate('estimate', estimate)
	end if

	! disconnect output files, deallocate memory etc...
	call cleanup()

	call MPI_FINALIZE()

contains

!============= utilities converted to subroutines for convenience =============
    
	! sets harvest to num/bins or num/bins+1. if called for 'ranks'
	! from 0 to bins-1, 'harvests' sum to 'num'
	subroutine load_command_line_args()
		
		if (command_argument_count() /= 4) then
			if (world_rank == 0) then
				print*, ' 4 cmd. line arguments need to be passed along &
with the executable : kernel count (int), sample size (int), number independent &
runs(int), scaling constant (real)'
				print*, ' terminating'
			end if
			call MPI_FINALIZE()
			call exit(1)
		end if

		call get_command_argument(1,arg)
		read(arg,*) gam

		call get_command_argument(2,arg)
		read(arg,*) sample_size

		call get_command_argument(3,arg)
		read(arg,*) runs

		call get_command_argument(4,arg)
		read(arg,*) lam

	end subroutine

	! set up MPI environment (ie. communicators and load balancing between
	! processes)
	subroutine mpi_setup()
		! obtaing rank of each process and their total number
		call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size)
		call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank)

		allocate(mu(d, gam))
		allocate(w(gam))
		allocate(sig(d, d, gam))
		allocate(det(gam))

		! allocate the sample
		allocate(sample(d,sample_size))

		! balance the load between processes
		call load_balance_count(world_rank, world_size, runs, node_runs)

		! for IO purposes, runs need to be strictly ordered.
		! 'offset runs' keeps track about position of node's first run in such
		! ordering (indexed from zero). other runs from the same node 
		! follow immediately after
		offset_runs = world_rank * (runs / world_size) + &
					  min(mod(runs, world_size), world_rank)

		mu_offset = 3*INT_SIZE + offset_runs * d * gam * DP_SIZE
		w_offset = mu_offset + (runs * d + offset_runs) * gam * DP_SIZE
		sig_offset = w_offset + (runs + offset_runs * d * d) * gam * DP_SIZE
		det_offset = sig_offset + (runs * d * d + offset_runs) * gam * DP_SIZE

		! connect output files (using MPI IO)
		call MPI_FILE_OPEN(MPI_COMM_WORLD, 'estimator', &
		                   ior(MPI_MODE_CREATE, MPI_MODE_WRONLY), &
		                   MPI_INFO_NULL, estimator_file)

		call MPI_FILE_WRITE(estimator_file, d, 1, MPI_INTEGER, &
		                    estimator_file_stat)

		call MPI_FILE_WRITE(estimator_file, gam, 1, MPI_INTEGER, &
		                    estimator_file_stat)

		call MPI_FILE_WRITE(estimator_file, runs, 1, MPI_INTEGER, &
		                    estimator_file_stat)


	end subroutine

	! subroutine to clean up after 'mpi_setup' above, it's mostly not neccesary
	! (OS would have done it) but even if nothing else at least leak sanitizer
	! is not fed false positives
	subroutine cleanup()
		
		deallocate(mu)
		deallocate(w)
		deallocate(sig)
		deallocate(det)

		deallocate(sample)

		call MPI_FILE_CLOSE(estimator_file)
		
	end subroutine

	! outputs grid_n, Ef, Ef2, Es, Es2, grid, exact (if available) into
	! unformatted binary file (named 'estimate' unless specified otherwise)
	subroutine save_estimate(filename, exact_density)

		real(kind=DP), intent(in), optional :: exact_density(:)
		character(len=*), intent(in) :: filename
		open(newunit=estimate_file_unit, file=filename, access='stream', &
		     form='unformatted', status='replace')

		write(estimate_file_unit) grid_n, Ef, Ef2, Es, Es2, grid

		if (present(exact_density)) write(estimate_file_unit) exact_density

		close(estimate_file_unit)

	end subroutine
	
	! write estimator variables to predefined locations in output file
	subroutine append_estimator()
		
		call MPI_FILE_SEEK(estimator_file, int(mu_offset,MPI_OFFSET_KIND), MPI_SEEK_SET)
		call MPI_FILE_WRITE(estimator_file, mu, gam*d, MPI_REAL8, &
		                    estimator_file_stat)

		call MPI_FILE_SEEK(estimator_file, int(w_offset,MPI_OFFSET_KIND), MPI_SEEK_SET)
		call MPI_FILE_WRITE(estimator_file, w, gam, MPI_REAL8, &
		                    estimator_file_stat)

		call MPI_FILE_SEEK(estimator_file, int(sig_offset,MPI_OFFSET_KIND), MPI_SEEK_SET)
		call MPI_FILE_WRITE(estimator_file, sig, gam*d, MPI_REAL8, &
		                    estimator_file_stat)

		call MPI_FILE_SEEK(estimator_file, int(det_offset,MPI_OFFSET_KIND), MPI_SEEK_SET)
		call MPI_FILE_WRITE(estimator_file, det, gam, MPI_REAL8, &
		                    estimator_file_stat)

		mu_offset = mu_offset + d * gam * DP_SIZE
		w_offset = w_offset + gam * DP_SIZE
		sig_offset = sig_offset + d * d * gam * DP_SIZE
		det_offset = det_offset + gam * DP_SIZE

	end subroutine

	pure subroutine load_balance_count(rank, bins, num, harvest)
		integer, intent(in) :: rank, bins, num
		integer, intent(out) :: harvest
		if (rank < mod(num, bins)) then
			harvest = num / bins + 1
		else
			harvest = num / bins
		end if
	end subroutine

end program akde_mse_test
