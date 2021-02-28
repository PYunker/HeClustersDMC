program akde_mse_test
    use metropolis
    use grid_utils
    use akde_botev
    use mpi_f08
    use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
    implicit none
    integer, parameter :: DP = KIND(0.0D0), d = 2, sample_size = 1e3, &
                          grid_n = 121, runs = 10, gam = 5
    integer :: i, j, world_size, node_size, run, nnodes, node_runs, n
    real(kind=DP) :: Ef(grid_n) = 0, Ef2(grid_n) = 0, estimate(grid_n)!, &
    !                 grid(d, grid_n)
    real(kind=DP), pointer :: grid(:,:)
    real(kind=DP), allocatable :: sample(:,:)

    type(process_type) :: process
    type(MPI_WIN) :: grid_win
    type(c_ptr) :: baseptr

    call MPI_INIT()
    ! set up communicators and balance the load between procesess
    call mpi_setup()
    write(*,"('node rank : ', i0,', / ', i0 , ' ,n_proc : ', i0)") &
        process%node_rank, node_size, process%n_proc

    ! initialize grid
    call win_allocate(d*grid_n, dp_size, process, grid_win, baseptr)
    call c_f_pointer(baseptr, grid, shape=[d, grid_n])
    
    if (process%node_rank == 0) then
        call grid_init_2d(grid, -5.0D0, 5.0D0, -5.0D0, 5.0D0, 11, 11)
    end if

    ! main cycle - for every run, draw new sample and obtain estimate
    do run=1,node_runs
        call draw_sample(sample, process%n_proc, d, pdf)
        ! output (estimate) is only relevant in node root and garbage elsewhere
        call botev_get_estimate(d, gam, sample_size, 9.2d0, process, sample, &
            grid, grid_n, estimate)
        if (process%node_rank == 0) then
            Ef = Ef + estimate / runs
            Ef2 = Ef2 + estimate**2 / runs
        end if
    end do
    
    if (process%world_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE, Ef, grid_n, MPI_REAL8, &
                        MPI_SUM, 0, process%caste_comm)
        call MPI_REDUCE(MPI_IN_PLACE, Ef2, grid_n, MPI_REAL8, &
                        MPI_SUM, 0, process%caste_comm)

        ! save Ef and Ef2 into file for further examination
        call save_ef_ef2('output')
    else if (process%node_rank == 0) then
        call MPI_REDUCE(Ef, Ef, grid_n, MPI_REAL8, &
                        MPI_SUM, 0, process%caste_comm)
        call MPI_REDUCE(Ef2, Ef2, grid_n, MPI_REAL8, &
                        MPI_SUM, 0, process%caste_comm)
    end if

    call MPI_FINALIZE()

contains

    ! the benchmark PDF to be approximated using AKDE
    real(kind=DP) pure function pdf(x)
        real(kind=DP), intent(in) :: x(:)
        pdf = exp(-0.5*norm2(x)**2)
        !out = exp(-1.5*(x-1.5)**2) + exp(-1.5*(x+1.5)**2)
        ! pdf = (0.65d0+0.35d0*cos(6.5d0*x(1)))*exp(-0.5d0*x(1)**2)
        return
    end function

!============= utilities converted to subroutines for convenience =============
    
    ! sets harvest to num/bins or num/bins+1. if called for 'ranks'
    ! from 0 to bins-1, 'harvests' sum to 'num'
    pure subroutine load_balance_count(rank, bins, num, harvest)
        integer, intent(in) :: rank, bins, num
        integer, intent(out) :: harvest
        if (rank < mod(num, bins)) then
            harvest = num / bins + 1
        else
            harvest = num / bins
        end if
    end subroutine

    ! this subroutine divides the range 1,2,...,num into 'bins' consecutive
    ! segments. harvest is set to 'rank'-th such segment (indexing from 0)
    pure subroutine load_balance_vec(rank, bins, num, harvest)
        integer, intent(in) :: rank, bins, num
        integer, allocatable, intent(out) :: harvest(:)
        integer :: remainder, base, low, high
        remainder = mod(num, bins)
        base = num / bins
        low = rank * base + min(remainder, rank) + 1
        high = (rank + 1) * base + min(remainder, rank + 1)
        allocate(harvest(high - low + 1))
        ! use 'base' to iterate through implied do, to avoid
        ! having to declare another integer just for that purpose
        harvest = [(base, base=low,high)]
    end subroutine

    ! set up MPI environment (ie. communicators and load balancing between
    ! processes)
    subroutine mpi_setup()
        ! obtaing rank of each process and their total number
        call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size)
        call MPI_COMM_RANK(MPI_COMM_WORLD, process%world_rank)
        
        ! setting up node communicator
        call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, &
            process%world_rank, MPI_INFO_NULL, process%node_comm)
        call MPI_COMM_SIZE(process%node_comm, node_size)
        call MPI_COMM_RANK(process%node_comm, process%node_rank)
        
        ! it is enough to have Ef and Ef2 allocated in node master procesess
        if (process%node_rank == 0) then
            i = 1
        else
            i = 0
        end if
        ! setting up 'caste' comm (ie. comm of node masters/slaves)
        call MPI_COMM_SPLIT(MPI_COMM_WORLD, i, process%world_rank, &
                            process%caste_comm)

        ! and load balancing runs between nodes    
        if (process%node_rank == 0) then
            call MPI_COMM_SIZE(process%caste_comm, nnodes)
            ! caste ranked is stored in 'j', since its needed only once and 
            ! it would be unnecessary to declare a saparate variable for it
            call MPI_COMM_RANK(process%caste_comm, j)
            call load_balance_count(j, nnodes, runs, node_runs)
        end if
        call MPI_BCAST(node_runs, 1, MPI_INTEGER, 0, process%node_comm)
        
        ! spliting n (sample_size) samples across procesess in one node
        call load_balance_count(process%node_rank, node_size, sample_size, &
                                process%n_proc)

        ! allocating memory for part of full sample assigned to calling process
        allocate(sample(d, process%n_proc))
        
        ! assign a process its portion of do i=1,gam loops
        call load_balance_vec(process%node_rank, node_size, gam, &
                              process%lb_idx)
    end subroutine
    
    ! prints output (mean value Ef and Ef**2) into file with Octave headers
    subroutine save_ef_ef2(jmeno)
        character(len=*), intent(in) :: jmeno
        open(1, file=jmeno)
        write(1,*) '# name: Ef'
        write(1,*) '# type: matrix'
        write(1,*) '# rows: 1'
        write(1,*) '# columns', sample_size
        write(1,*) Ef

        write(1,*) '# name: Ef2'
        write(1,*) '# type: matrix'
        write(1,*) '# rows: 1'
        write(1,*) '# columns', grid_n
        write(1,*) Ef2

        write(1,*) '# name: grid'
        write(1,*) '# type: matrix'
        write(1,*) '# rows: ', d
        write(1,*) '# columns', grid_n
        write(1,*) grid

        close(1)
    end subroutine

end program akde_mse_test