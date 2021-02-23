#ifndef d
#error "dimension 'd' not specified"
#endif
#ifndef gam       
#error "kernel count 'gam' not specified"
#endif
program neco
    use mpi_f08
    use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
    implicit none

    integer, parameter :: DP = KIND(0.0D0), eps = epsilon(0.0D0)
    integer, parameter :: dp_size = sizeof(0.0D0), l_size = sizeof(.true.)

    integer, parameter :: max_iter = 1000, runs = 100, &
                          grid_n = 50, grid_n2 = grid_n**2
                          ! gam = 5, d = 2
    !include 'kde_hlavicka.f95' ! deklarace d & gam

    real(kind=DP), parameter :: lam = 9.2d0, lam2=sqrt(1/lam)**d, &
                                konst2 = 1d0/sqrt(8d0*atan(1d0))**d

    integer :: world_rank, world_size, node_rank, node_size, disp_unit, ierror, &
               seed(4), i, j, k, run, idx(100) = [(i, i=1,100)], nnodes, &
               node_runs, n, n_proc
    
    real(kind=DP), pointer :: mu(:,:), w(:), sig(:,:,:), h, det(:), stopa(:), &
                              grid(:,:)
    logical, pointer :: ent

    real(kind=DP) :: mu_castecne(d,gam), w_castecne(gam), sig_castecne(d,d,gam)

    !real(kind=DP) :: seed_real(4), X(d,nmax), xRinv(d,nmax), curv, dlantr, &
    !                 p(nmax,gam), psig(nmax,gam), density(nmax), xSig(nmax), &
    !                 mesh(d,grid_n2), temp(grid_n2), odhad(grid_n2), ent_c, konst

    real(kind=DP) :: seed_real(4), curv, dlantr, ent_c, konst, &
                     mesh(d,grid_n2), temp(grid_n2), odhad(grid_n2)

    real(kind=DP), allocatable :: Ef(:), Ef2(:), X(:,:), xRinv(:,:), &
                                  p(:,:), psig(:,:), density(:), xSig(:)

    type(MPI_COMM) :: node_comm
    type(MPI_COMM) :: caste_comm
    type(MPI_WIN) :: win_mu
    type(c_ptr) :: baseptr

    ! =============== inicializace komunikatoru ================
    ! ------------------ globalni komunikator ------------------
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, ierror)
    ! -------------------- node komunikator --------------------
    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, &
                             world_rank, MPI_INFO_NULL, node_comm, ierror)
    call MPI_COMM_SIZE(node_comm, node_size, ierror)
    call MPI_COMM_RANK(node_comm, node_rank, ierror)
    ! --------- kastovy komunikator + load balancing runu ------------------
    if (node_rank == 0) then
        i = 1
    else
        i = 0
    end if
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, i, world_rank, caste_comm)
    if (node_rank == 0) then
        call MPI_COMM_SIZE(caste_comm, nnodes)
        call MPI_COMM_RANK(caste_comm, j) ! caste rank docasne ulozen v j
        node_runs = runs / nnodes
        if (j < mod(runs,nnodes)) node_runs = node_runs + 1
    end if
    call MPI_BCAST(node_runs, 1, MPI_INTEGER, 0, node_comm, ierror)
    ! ----------------------------------------------------------
    if (world_rank == 0) then
        write(*,'(A)',advance='no') 'n : '; read(*,*) n
    end if
    call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)
    konst = 16d0*n*sqrt(atan(1d0))**2d0
    n_proc = n / node_size
    if (node_rank < mod(n, node_size)) n_proc = n_proc + 1
    !if (n_proc > nmax) then
    !    write(*,"('nodostatecne nmax, ', i6, ' vs ', i6)") n_proc, nmax
    !    preteceni_lokalni = .true.
    !end if
    !call MPI_ALLREDUCE(preteceni_lokalni, preteceni_globalni, 1, &
    !                   MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD)
    !if (preteceni_globalni) then ! pokud je nejakemu procesu prideleno vice nez n_max vzorku
    !    call MPI_FINALIZE(ierror)
    !    call exit(0)
    !end if
    write(*,"('node rank : ', i0,', / ', i0 , ' ,n_proc : ', i0)") node_rank, node_size, n_proc
    allocate(X(d,n_proc))
    allocate(xRinv(d,n_proc))
    allocate(p(n_proc,gam))
    allocate(psig(n_proc,gam))
    allocate(density(n_proc))
    allocate(xSig(n_proc))
    ! ================= alokace sdilene pameti =================
    ! -------------------------- mu ----------------------------
    call win_allocate(d*gam, dp_size)
    call c_f_pointer(baseptr, mu, shape=[d, gam])
    ! -------------------------- Sig ---------------------------
    call win_allocate(d*d*gam, dp_size)
    call c_f_pointer(baseptr, sig, shape=[d,d,gam])
    ! -------------------------- w -----------------------------
    call win_allocate(gam, dp_size)
    call c_f_pointer(baseptr, w, shape=[gam])
    ! -------------------------- h -----------------------------
    call win_allocate(1, dp_size)
    call c_f_pointer(baseptr, h)
    ! ------------------------- ent ----------------------------
    call win_allocate(1, l_size)
    call c_f_pointer(baseptr, ent)
    ! ------------------------ stopa ---------------------------
    call win_allocate(gam, dp_size)
    call c_f_pointer(baseptr, stopa, shape=[gam])
    ! ------------------------- det ----------------------------
    call win_allocate(gam, dp_size)
    call c_f_pointer(baseptr, det, shape=[gam])
    ! ------------------------ grid ----------------------------
    call win_allocate(d*grid_n2, dp_size)
    call c_f_pointer(baseptr, grid, shape=[d,grid_n2])
    ! ================ incializace mrizky, ... =================
    if (node_rank == 0) then
        allocate(Ef(grid_n2)) ! protoze mimo rank 0 nejsou zapotrebi
        allocate(Ef2(grid_n2))
        Ef = 0
        Ef2 = 0
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
    ! ##################### hlavni cyklus ######################
    do run=1,node_runs
        ! ================ incializace mu, w, sig, X ===============
        if (node_rank == 0) then
            ! -------------------------- mu ----------------------------
            call random_number(seed_real)
            seed = floor(4095*seed_real)
            if (mod(seed(4),2) == 0) seed(4) = seed(4) + 1
            call dlarnv(3, seed, d*gam, mu)
        ! ---------------------- w, sig & h ------------------------
            sig = 0
            h = 0.1d0*lam*n**(-real(d,8)/(d+4))
            do i = 1,d
                call random_number(w) ! tohle nemusi byt napsane pres w
                w = 0.9d0 * w + 0.1d0
                w = w * h
                sig(i,i,:) = w
            end do
            sig = sig * lam**2
            call random_number(w)
            w = w / sum(w)
            det = 1 ! incializuje cele det na jednicky
            ent = .false. ! jinak vsechno hned skonci
        end if
        ! --------------------------- X ----------------------------
        ! if (world_rank == 0) write(*,"('============= vygenerovani X =============')")
        ! call MPI_BARRIER(MPI_COMM_WORLD) ! jenom, aby byla zprava na zacatku
        call random_number(seed_real)
        seed = floor(4095*seed_real)
        if (mod(seed(4),2) == 0) seed(4) = seed(4) + 1
        ! write(*,"('RANK : ', i2, ' ,seed :', 4i6)") world_rank, seed ! dale musi vypisovat do souboru
        call dlarnv(3, seed, d*n_proc, X)
        density(idx) = 0 ! kvuli zastavovaci podmince
        ! ===================== samotny vypocet ====================
        call MPI_BARRIER(node_comm, ierror) ! nejspis neni nutna
        k = 1 ! NEZAPOMENOUT
        do! k=1,max_iter
            det = 1
            do i=1,gam
                if (mod(i,node_size) == node_rank) then ! prepsat pomoci integer division
                    !call dpotrf('L', d, sig(:,:,i), d, j) ! LDA by melo byt dobre
                    call chol(sig(:,:,i))
                    do j=1,d ! asi muze byt concurrent
                        det(i) = det(i) * sig(j,j,i)
                    end do
                    !call dtrtri('L','N',d,sig(:,:,i),d,j)
                    call inverse(sig(:,:,i))
                    !stopa(i) = dlantr('F','L','N',d,d,sig(:,:,i),d,curv)**2
                    stopa(i) = norma(sig(:,:,i))
                end if
            end do
            call MPI_BARRIER(node_comm) ! master uz ma hotovo
            
            if (ent .or. k == max_iter) exit ! je to tady a ne na konci, protoze proto
    
            do i=1,gam ! mozna muze byt concurrent
                do concurrent (j=1:n_proc) ! asi muze byt do concurrent (j=1:n)
                    xRinv(:,j) = X(:,j) - mu(:,i)
                end do
                call dtrmm('L','L','N','N',d,n_proc,1d0,sig(:,:,i),d,xRinv,d)
                xSig = norm2(xRinv,dim=1)**2d0 ! tohle je asi zbytecne rozepsane
                p(:,i) = w(i) / det(i) * exp(-0.5 * xSig - 0.5 * stopa(i) * h**2d0)
                call dtrmm('L','L','T','N',d,n_proc,lam2,sig(:,:,i),d,xRinv,d)
                xSig = norm2(xRinv,dim=1)**2d0 ! melo by byt +eps
                psig(:,i) = p(:,i) * xSig
            end do
            ! misto idx muze byt proste 1:100
            ent_c = product(density(idx)) ! VNIMANIE VNIANIE
            
            density = sum(p,dim=2) ! + eps ! eps mozna nemusi byt
            do i=1,gam ! asi muze byt concurrent
                p(:,i) = p(:,i) / density
            end do
            ! misto idx muze byt proste 1:100
            ent_c = ent_c / product(density(idx)) ! VYPOCET (CASTECNEHO) ENT, trochu oklikou
    
            w_castecne = sum(p,dim=1)
            do concurrent (i=1:gam)
                mu_castecne(:,i) = matmul(X,p(:,i))
            end do
    
            if (node_rank == 0) then
                call MPI_REDUCE(MPI_IN_PLACE, mu_castecne, d*gam, &
                                MPI_REAL8 , MPI_SUM, 0, node_comm, ierror)
                call MPI_REDUCE(MPI_IN_PLACE, w_castecne, gam, &
                                MPI_REAL8 , MPI_SUM, 0, node_comm, ierror)
                do concurrent (i=1:gam, w_castecne(i) /= 0) ! nezapomenout na podminku
                    mu(:,i) = mu_castecne(:,i) / w_castecne(i)
                end do
            else
                call MPI_REDUCE(mu_castecne, mu_castecne, d*gam, &
                                MPI_REAL8, MPI_SUM, 0, node_comm, ierror)
                call MPI_REDUCE(w_castecne, w_castecne, gam, &
                                MPI_REAL8, MPI_SUM, 0, node_comm, ierror)
            end if
            
            call MPI_BARRIER(node_comm) ! delnici uz maji hotovo
    
            do i=1,gam
                do concurrent (j=1:n_proc) ! asi muze byt do concurrent (j=1:n)
                    xRinv(:,j) = X(:,j) - mu(:,i)
                end do
                do concurrent (j=1:d)
                    xRinv(j,:) = xRinv(j,:) * sqrt(p(:,i))
                end do
                sig_castecne(:,:,i) = matmul(xRinv, transpose(xRinv)) ! tohle nevypada nejlepe
            end do
            
            if (node_rank == 0) then
                call MPI_REDUCE(MPI_IN_PLACE, sig_castecne, d*d*gam, &
                                MPI_REAL8 , MPI_SUM, 0, node_comm, ierror)
                do concurrent (i=1:gam, w_castecne(i) /= 0) ! nezapomenout na podminku
                    sig_castecne(:,:,i) = sig_castecne(:,:,i) / w_castecne(i)
                    do concurrent (j=1:d)
                        sig(j,j,i) = sig_castecne(j,j,i) + h**2
                    end do
                end do
            else
                call MPI_REDUCE(sig_castecne, sig_castecne, d*d*gam, &
                                MPI_REAL8 , MPI_SUM, 0, node_comm, ierror)
            end if
            
            xSig = sum(psig,dim=2)/density        
            curv = sum(xSig)/n
    
            if (node_rank == 0) then
                call MPI_REDUCE(MPI_IN_PLACE, ent_c, 1, &
                                MPI_REAL8 , MPI_PROD, 0, node_comm, ierror)
                call MPI_REDUCE(MPI_IN_PLACE, curv, 1, &
                                MPI_REAL8 , MPI_SUM, 0, node_comm, ierror)
                ent = abs(ent_c - 1) < 1e-6 ! TOLERANCE JE NASTAVENA NA PEVNO
                h = (curv*konst)**(-real(1,kind=DP)/(d+2))
                w = w_castecne / sum(w_castecne)
            else
                call MPI_REDUCE(ent_c, ent_c, 1, &
                                MPI_REAL8 , MPI_PROD, 0, node_comm, ierror)
                call MPI_REDUCE(curv, curv, 1, &
                                MPI_REAL8 , MPI_SUM, 0, node_comm, ierror)
            end if
            k = k + 1
        end do
        ! ============= vypocet mse Ef & Ef2 na mrizce =============
        odhad = 0
        do i=1,gam
            if (mod(i,node_size) == node_rank) then
                do concurrent (j=1:grid_n2)
                    mesh(:,j) = grid(:,j) - mu(:,i)
                end do
                call dtrmm('L','L','N','N',d,grid_n2,1d0,sig(:,:,i),d,mesh,d) ! mesh := L * mesh
                temp = -0.5d0 * norm2(mesh,dim=1)**2d0
                odhad = odhad + w(i) * konst2 / det(i) * exp(temp)
            end if
        end do
        ! =================== DEBUG =====================
        ! if (any(isnan(odhad))) write(*,"('v odhad pojavilsja NaN, proces : ', i2)") node_rank
        if (node_rank == 0) then
            call MPI_REDUCE(MPI_IN_PLACE, odhad, grid_n2, MPI_REAL8, &
                            MPI_SUM, 0, node_comm, ierror)
            Ef = Ef + odhad/runs
            Ef2 = Ef2 + (odhad**2)/runs
        else
            call MPI_REDUCE(odhad, odhad, grid_n2, MPI_REAL8, &
                            MPI_SUM, 0, node_comm, ierror)
        end if
    end do ! ted uz ma kazdy node vypoctene Ef ze svojeho runu (z jednech parammetru)
    ! -------  globalni root posbira data za vsech nodu --------
    if (world_rank == 0) then
        call MPI_REDUCE(MPI_IN_PLACE, Ef, grid_n2, MPI_REAL8, &
                        MPI_SUM, 0, caste_comm, ierror)
        call MPI_REDUCE(MPI_IN_PLACE, Ef2, grid_n2, MPI_REAL8, &
                        MPI_SUM, 0, caste_comm, ierror)
    else if (node_rank == 0) then
        call MPI_REDUCE(Ef, Ef, grid_n2, MPI_REAL8, &
                        MPI_SUM, 0, caste_comm, ierror)
        call MPI_REDUCE(Ef2, Ef2, grid_n2, MPI_REAL8, &
                        MPI_SUM, 0, caste_comm, ierror)
    end if
    ! ===================== zapis vysledku =====================
    if (world_rank == 0) then
        call uloz_ef('vystup')
        call uloz_kde('vystup_kde')
    end if
    ! ======================== cleanup =========================
    deallocate(X)
    deallocate(xRinv)
    deallocate(p)
    deallocate(psig)
    deallocate(density)
    deallocate(xSig)
    if (node_rank == 0) then 
        deallocate(Ef)
        deallocate(Ef2)
    end if
    call MPI_WIN_FREE(win_mu, ierror)
    call MPI_COMM_FREE(caste_comm, ierror)
    call MPI_COMM_FREE(node_comm, ierror)
    call MPI_FINALIZE(ierror)
    
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
            ! call MPI_WIN_SHARED_QUERY(win_mu, 0, var_size, disp_unit, baseptr)
            call MPI_WIN_SHARED_QUERY(win_mu, 0, window_size, disp_unit, baseptr)
        end if
    
    end subroutine win_allocate
    ! -------------------- 2d rutinky ----------------------
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

    pure function norma(a) result(v)
        real(kind=DP), intent(in) :: a(2,2)
        real(kind=DP) :: v
        v = a(1,1)**2 + a(2,1)**2 + a(2,2)**2
    end function norma
    ! ------------------------------------------------------
    subroutine uloz_kde(jmeno)
        character(len=*), intent(in) :: jmeno
        open(1, file=jmeno)
        write(1,*) ! vypis mu
        write(1,*) '# name: mu'
        write(1,*) '# type: matrix'
        write(1,*) '# rows:', d
        write(1,*) '# columns', gam
        do i=1,d
            write(1,*) mu(i,:)
        end do
        
        write(1,*) ! vypis w
        write(1,*) '# name: w'
        write(1,*) '# type: matrix'
        write(1,*) '# rows: 1'
        write(1,*) '# columns', gam
        do i=1,gam
            write(1,*) w(i)
        end do

        write(1,*) ! vypis det
        write(1,*) '# name: det'
        write(1,*) '# type: matrix'
        write(1,*) '# rows: 1'
        write(1,*) '# columns', gam
        write(1,*) det
        
        write(1,*) ! vypis sig
        write(1,*) '# name: sig'
        write(1,*) '# type: matrix'
        write(1,*) '# ndims:', 3
        write(1,*) d, d, gam
        do k=1,gam
            do j=1,d
                do i=1,d
                    write(1,*) sig(i,j,k)
                end do
            end do
        end do
    end subroutine uloz_kde

    subroutine uloz_ef(jmeno)
        character(len=*), intent(in) :: jmeno
        open(1, file=jmeno)
        write(1,*) '# name: Ef'
        write(1,*) '# type: matrix'
        write(1,*) '# rows: 1'
        write(1,*) '# columns', grid_n2
        write(1,*) Ef

        write(1,*) '# name: Ef2'
        write(1,*) '# type: matrix'
        write(1,*) '# rows: 1'
        write(1,*) '# columns', grid_n2
        write(1,*) Ef2
        close(1)
    end subroutine uloz_ef

    !subroutine bezmezer(string)
    !    character(len=*) :: string
    !    integer :: delka 
    !    integer :: konec, kurzor
    !    delka = len(string)
    !    konec = 1
    !    kurzor = 1
    !    do while (kurzor < delka)
    !        if (string(konec:konec) == ' ') then
    !            kurzor = kurzor + 1
    !            string(konec:konec) = string(kurzor:kurzor)
    !            string(kurzor:kurzor) = ' '
    !        else
    !            konec = konec + 1
    !            if (kurzor < konec) kurzor = konec
    !        endif
    !    end do
    !end subroutine

end program neco