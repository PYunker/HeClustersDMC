program standart
    use mpi_f08
    use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
    implicit none
    
    integer, parameter :: DP = KIND(0.0D0), eps = epsilon(0.0D0)
    integer, parameter :: dp_size = sizeof(0.0D0)

    integer, parameter :: d = 2, grid_n = 50, grid_n2 = grid_n**d

    real(kind=DP), parameter :: konst = 1d0/sqrt(8d0*atan(1d0))**d

    integer :: world_rank, world_size, node_rank, node_size, disp_unit, &
               seed(4), i, j, k, idx(100) = [(i, i=1,100)], nnodes, &
               n, n_proc
    
    real(kind=DP), pointer :: grid(:,:)

    real(kind=DP) :: seed_real(4), sig(d,d), Ef(grid_n2), Ef2(grid_n2), &
                     mesh(d,grid_n2), temp(grid_n2), odhad(grid_n2), det, h

    real(kind=DP), allocatable :: X(:,:)

    type(MPI_COMM) :: node_comm
    type(MPI_WIN) :: win_mu
    type(c_ptr) :: baseptr

    call MPI_INIT()
    call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size)
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank)
    ! -------------------- node komunikator --------------------
    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, &
                             world_rank, MPI_INFO_NULL, node_comm)
    call MPI_COMM_SIZE(node_comm, node_size)
    call MPI_COMM_RANK(node_comm, node_rank)
    ! ---------------- kastovy komunikator ------------------
    ! nacte a rozdeli n
    if (world_rank == 0) then
        write(*,'(A)',advance='no') 'n : '; read(*,*) n
    end if
    call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)
    n_proc = n / world_size
    if (node_rank < mod(n, world_size)) n_proc = n_proc + 1
    write(*,"('world rank : ', i0,', / ', i0 , ' ,n_proc : ', i0)") &
        world_rank, world_size, n_proc
    allocate(X(d,n_proc))
    ! =============== incicializace a sdileni mrizky ===============
    call win_allocate(d*grid_n2, dp_size)
    call c_f_pointer(baseptr, grid, shape=[d,grid_n2])
    
    if (node_rank == 0) then
        ! ------------------------ mrizka --------------------------
        k = 1
        do i=1,grid_n
            do j=1,grid_n
                grid(:,k) = [i,j]
                k = k + 1
            end do
        end do
        ! ax + b
        ! a = (horni - dolni)/(grid_n - 1)
        ! b = (grid_n*dolni-horni)/(grid_n - 1)
        grid(1,:) = (10d0 * grid(1,:)  - 5d0*(grid_n + 1)) / (grid_n - 1)
        grid(2,:) = (10d0 * grid(2,:)  - 5d0*(grid_n + 1)) / (grid_n - 1)
    end if
    ! =================== incializace dat =======================
    call random_number(seed_real)
        seed = floor(4095*seed_real)
        if (mod(seed(4),2) == 0) seed(4) = seed(4) + 1
        ! write(*,"('RANK : ', i2, ' ,seed :', 4i6)") world_rank, seed ! dale musi vypisovat do souboru
    call dlarnv(3, seed, d*n_proc, X)
    ! ================== incializace polosirky ===================
    h = (4/(d+2)/n)**(1/(d+4))
    !h = h/10
    !h = (4.d0/3/n)**(0.2d0)
    !do j=1,d
    !    sig(j,j) = h**2
    !end do
    !call chol(sig)
    !call dpotrf('L', d, sig, d, j)
    !call inverse(sig)
    !call dtrtri('L','N',d, sig, d, j)
    det = konst/n/h**d
    call MPI_BARRIER(MPI_COMM_WORLD)
    odhad = 0
    do i=1,n_proc
        do concurrent (j=1:grid_n2)
            mesh(:,j) = grid(:,j) - X(:,i)
        end do
        !call dtrmm('L','L','N','N',d,grid_n2,1d0,sig,d,mesh,d) ! mesh := L * mesh
        temp = -0.5d0 * (norm2(mesh,dim=1)/h)**2d0
        odhad = odhad + det * exp(temp)
    end do
    
    if (world_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE, odhad, grid_n2, MPI_REAL8, &
                        MPI_SUM, 0, MPI_COMM_WORLD)
        Ef2 = odhad**2
        call uloz_ef('data')
        print*, grid(:,grid_n2)
    else
        call MPI_REDUCE(odhad, odhad, grid_n2, MPI_REAL8, &
                        MPI_SUM, 0, MPI_COMM_WORLD)
    end if

    deallocate(X)
    call MPI_WIN_FREE(win_mu)
    call MPI_COMM_FREE(node_comm)
    call MPI_FINALIZE()

contains
    
    subroutine win_allocate(numel, btsize)
            
        integer, intent(in) :: numel, btsize
        
        integer(kind=MPI_ADDRESS_KIND) :: window_size

        if (node_rank == 0) then 
            window_size = numel*btsize ! 
        else
            window_size = 0
        end if
        call MPI_WIN_ALLOCATE_SHARED(window_size, btsize, MPI_INFO_NULL, &
                                    node_comm, baseptr, win_mu)
        if (node_rank /= 0) then
            call MPI_WIN_SHARED_QUERY(win_mu, 0, window_size, disp_unit, baseptr)
        end if

    end subroutine win_allocate

    pure subroutine inverse(a)
        real(kind=DP), intent(inout) :: a(2,2)
        real(kind=DP) :: detinv, tmp
        detinv = 1/(a(1,1)*a(2,2)-a(1,2)*a(2,1))
        tmp = a(1,1)
        a(1,1) = detinv * a(2,2)
        a(1,2) = -detinv * a(1,2)
        a(2,1) = -detinv * a(2,1)
        a(2,2) = detinv * tmp
    end subroutine inverse

    pure subroutine chol(a)
        real(kind=DP), intent(inout) :: a(2,2)
        a(1,1) = sqrt(a(1,1))
        a(2,1) = a(2,1) / a(1,1)
        a(2,2) = sqrt(a(2,2)-(a(1,2)/a(1,1))**2)
        a(1,2) = 0
    end subroutine chol

    subroutine uloz_ef(jmeno)
        character(len=*), intent(in) :: jmeno
        open(1, file=jmeno)
        write(1,*) '# name: Ef'
        write(1,*) '# type: matrix'
        write(1,*) '# rows: 1'
        write(1,*) '# columns', grid_n2
        write(1,*) odhad

        write(1,*) '# name: Ef2'
        write(1,*) '# type: matrix'
        write(1,*) '# rows: 1'
        write(1,*) '# columns', grid_n2
        write(1,*) Ef2
        close(1)
    end subroutine uloz_ef

end program standart