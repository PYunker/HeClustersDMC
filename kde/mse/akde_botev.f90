!------------------------------------------------------------------------------
! MODULE: akde_botev
!
!> @author
!> Pavel Junker
!
! DESCRIPTION: 
!>  implementation of Zdravko Botev's adaptive kernel density estimation algor-
! ithm (https://www.mathworks.com/matlabcentral/fileexchange/58312-kernel-densi
! ty-estimator-for-high-dimensions)
!
! REVISION HISTORY:
! 28.02.2021 - Initial Version
!------------------------------------------------------------------------------
module akde_botev
    use mpi_f08
    use metropolis
    use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
    ! use grid
    implicit none
    private

    integer, parameter :: DP = kind(0.D0), dp_size = sizeof(0.0D0), &
                          l_size = sizeof(.true.)

    real(kind=DP), parameter :: eps = epsilon(0.0D0)

    type :: process_type
        integer :: world_rank
        integer :: node_rank
        type(MPI_COMM) :: node_comm
        type(MPI_COMM) :: caste_comm
        integer :: n_proc
        integer :: lb_from
        integer :: lb_to
    end type process_type

    public :: botev_get_estimate, win_allocate, process_type, dp_size

contains
    
    subroutine botev_get_estimate(d, gam, n, lam, process, X, grid, &
                                  grid_n, estimate)
        integer, intent(in) :: d, gam, n, grid_n
        type(process_type), intent(in) :: process
        real(kind=DP), intent(in) :: lam, grid(d,grid_n)
        real(kind=DP), intent(in) :: X(d,process%n_proc)
        real(kind=DP), intent(out) :: estimate(grid_n)
        ! internal variables
        type(MPI_WIN) :: shared_win
        type(c_ptr) :: baseptr
        integer :: i, j, k, idx(100) = [(i, i=1,100)], max_iter = 1000
        real(kind=DP) :: xRinv(d,process%n_proc), p(process%n_proc,gam), &
                         psig(process%n_proc,gam), density(process%n_proc), &
                         xSig(process%n_proc), temp(grid_n), curv, dlantr, &
                         ent_c, mesh(d,grid_n)!, seed_real(4)
        real(kind=DP) :: mu_castecne(d,gam), w_castecne(gam), &
                         sig_castecne(d,d,gam)
        ! these two act as parameters and are never changed
        real(kind=DP) :: norm_const, konst, lam2
        real(kind=DP), pointer :: mu(:,:), w(:), sig(:,:,:), h, det(:), stopa(:)
        logical, pointer :: ent

        ! some neccessary initialization
        density(idx) = 0
        norm_const = 1d0/sqrt(8d0*atan(1d0))**d
        konst = 16d0*n*sqrt(atan(1d0))**2d0
        lam2 = sqrt(1/lam)**d

        ! allocation of memory shared across node
        ! I have no idea how does this work, it just does (at least seems to)
        ! -------------------------- mu ----------------------------
        call win_allocate(d*gam, dp_size, process, shared_win, baseptr)
        call c_f_pointer(baseptr, mu, shape=[d, gam])
        ! -------------------------- Sig ---------------------------
        call win_allocate(d*d*gam, dp_size, process, shared_win, baseptr)
        call c_f_pointer(baseptr, sig, shape=[d,d,gam])
        ! -------------------------- w -----------------------------
        call win_allocate(gam, dp_size, process, shared_win, baseptr)
        call c_f_pointer(baseptr, w, shape=[gam])
        ! -------------------------- h -----------------------------
        call win_allocate(1, dp_size, process, shared_win, baseptr)
        call c_f_pointer(baseptr, h)
        ! ------------------------- ent ----------------------------
        call win_allocate(1, l_size, process, shared_win, baseptr)
        call c_f_pointer(baseptr, ent)
        ! ------------------------ stopa ---------------------------
        call win_allocate(gam, dp_size, process, shared_win, baseptr)
        call c_f_pointer(baseptr, stopa, shape=[gam])
        ! ------------------------- det ----------------------------
        call win_allocate(gam, dp_size, process, shared_win, baseptr)
        call c_f_pointer(baseptr, det, shape=[gam])

        ! inititalization of variables shared across nodes, done by node roots
        if (process%node_rank == 0) then
            ! it's assumed samples are pefectly shuffled
            mu = X(:,1:gam)
            h = 0.1d0*lam*n**(-real(d,8)/(d+4))
            sig = 0
            do i=1,d
                call random_number(sig(i,i,:))
                sig(i,i,:) = (0.1 + 0.9 * sig(i,i,:)) * h * lam**2
            end do
            call random_number(w)
            w = w / sum(w)
            det = 1
            ent = .false.
        end if

        call MPI_BARRIER(process%node_comm)
        ! aliasing process%node_rank to node_rank etc. for better readability
        associate (node_rank  => process%node_rank, &
                   node_comm  => process%node_comm, &
                   n_proc     => process%n_proc, &
                   lb_from    => process%lb_from, &
                   lb_to      => process%lb_to)

        k = 1
        do
            det = 1
            do i=lb_from,lb_to
                call dpotrf('L', d, sig(:,:,i), d, j)
                do j=1,d
                    det(i) = det(i) * sig(j,j,i)
                end do
                call dtrtri('L','N',d,sig(:,:,i),d,j)
                !stopa(i) = dlantr('F','L','N',d,d,sig(:,:,i),d,curv)**2
                stopa(i) = norma(sig(:,:,i))
            end do
            call MPI_BARRIER(node_comm)
            
            if (ent .or. k == max_iter) exit
    
            do i=1,gam
                do concurrent (j=1:n_proc)
                    xRinv(:,j) = X(:,j) - mu(:,i)
                end do
                
                ! =============== DEBUG ================

                if (any(isnan(sig))) then
                    write(*, "('v sig je nan (', i0, '), i = ', i0)") & 
                        count(isnan(sig)), i
                end if
                
                call dtrmm('L','L','N','N',d,n_proc,1d0,sig(:,:,i),d,xRinv,d)
                
                if (any(isnan(xRinv))) then
                    write(*, "('v xrinv je nan (', i0, '), i = ', i0)") & 
                        count(isnan(xRinv)), i
                end if

                xSig = norm2(xRinv,dim=1)**2d0
                p(:,i) = w(i) / det(i) * &
                         exp(-0.5 * xSig - 0.5 * stopa(i) * h**2d0)
                call dtrmm('L','L','T','N',d,n_proc,lam2,sig(:,:,i),d,xRinv,d)
                xSig = norm2(xRinv,dim=1)**2d0 ! + eps
                psig(:,i) = p(:,i) * xSig
            end do

            ent_c = product(density(idx))
            
            density = sum(p,dim=2) ! + eps
            do i=1,gam ! asi muze byt concurrent
                p(:,i) = p(:,i) / density
            end do

            ent_c = ent_c / product(density(idx))
    
            w_castecne = sum(p,dim=1)
            do concurrent (i=1:gam)
                mu_castecne(:,i) = matmul(X,p(:,i))
            end do
    
            if (node_rank == 0) then
                call MPI_REDUCE(MPI_IN_PLACE, mu_castecne, d*gam, &
                                MPI_REAL8 , MPI_SUM, 0, node_comm)
                call MPI_REDUCE(MPI_IN_PLACE, w_castecne, gam, &
                                MPI_REAL8 , MPI_SUM, 0, node_comm)
                do concurrent (i=1:gam, w_castecne(i) /= 0)
                    mu(:,i) = mu_castecne(:,i) / w_castecne(i)
                end do
            else
                call MPI_REDUCE(mu_castecne, mu_castecne, d*gam, &
                                MPI_REAL8, MPI_SUM, 0, node_comm)
                call MPI_REDUCE(w_castecne, w_castecne, gam, &
                                MPI_REAL8, MPI_SUM, 0, node_comm)
            end if
            
            call MPI_BARRIER(node_comm)
    
            do i=1,gam
                do concurrent (j=1:n_proc)
                    xRinv(:,j) = X(:,j) - mu(:,i)
                end do
                do concurrent (j=1:d)
                    xRinv(j,:) = xRinv(j,:) * sqrt(p(:,i))
                end do
                sig_castecne(:,:,i) = matmul(xRinv, transpose(xRinv))
            end do
            
            if (node_rank == 0) then
                call MPI_REDUCE(MPI_IN_PLACE, sig_castecne, d*d*gam, &
                                MPI_REAL8 , MPI_SUM, 0, node_comm)
                do concurrent (i=1:gam, w_castecne(i) /= 0)
                    sig_castecne(:,:,i) = sig_castecne(:,:,i) / w_castecne(i)
                    do concurrent (j=1:d)
                        sig(j,j,i) = sig_castecne(j,j,i) + h**2
                    end do
                end do
            else
                call MPI_REDUCE(sig_castecne, sig_castecne, d*d*gam, &
                                MPI_REAL8 , MPI_SUM, 0, node_comm)
            end if
            
            xSig = sum(psig,dim=2)/density        
            curv = sum(xSig)/n
    
            if (node_rank == 0) then
                call MPI_REDUCE(MPI_IN_PLACE, ent_c, 1, &
                                MPI_REAL8 , MPI_PROD, 0, node_comm)
                call MPI_REDUCE(MPI_IN_PLACE, curv, 1, &
                                MPI_REAL8 , MPI_SUM, 0, node_comm)
                ! tolerance is hardcoded
                ent = abs(ent_c - 1) < 1e-6
                h = (curv*konst)**(-real(1,kind=DP)/(d+2))
                w = w_castecne / sum(w_castecne)
            else
                call MPI_REDUCE(ent_c, ent_c, 1, &
                                MPI_REAL8 , MPI_PROD, 0, node_comm)
                call MPI_REDUCE(curv, curv, 1, &
                                MPI_REAL8 , MPI_SUM, 0, node_comm)
            end if
            k = k + 1
        end do
        ! ============= vypocet mse Ef & Ef2 na mrizce =============
        estimate = 0
        do i=lb_from,lb_to
            do concurrent (j=1:grid_n)
                mesh(:,j) = grid(:,j) - mu(:,i)
            end do
            call dtrmm('L','L','N','N',d,grid_n,1d0,sig(:,:,i),d,mesh,d)
            temp = -0.5d0 * norm2(mesh,dim=1)**2d0
            estimate = estimate + w(i) * norm_const / det(i) * exp(temp)
        end do

        if (node_rank == 0) then
            call MPI_REDUCE(MPI_IN_PLACE, estimate, grid_n, MPI_REAL8, &
                            MPI_SUM, 0, node_comm)
        else
            call MPI_REDUCE(estimate, estimate, grid_n, MPI_REAL8, &
                            MPI_SUM, 0, node_comm)
        end if

        end associate
        
        ! deallocate shared memory to prevent leaks
        deallocate(mu)
        deallocate(w)
        deallocate(sig)
        deallocate(h)
        deallocate(ent)
        deallocate(stopa)
        deallocate(det)
        
    end subroutine botev_get_estimate

    subroutine win_allocate(numel, btsize, process, window, baseptr)
            
        integer, intent(in) :: numel, btsize
        type(process_type) :: process
        integer :: disp_unit
        integer(kind=MPI_ADDRESS_KIND) :: window_size
        ! does it haave to be intent(inout) ? not sure
        type(MPI_WIN), intent(inout) :: window
        type(c_ptr), intent(out) :: baseptr 

        if (process%node_rank == 0) then 
            window_size = numel*btsize ! 
        else
            window_size = 0
        end if
        call MPI_WIN_ALLOCATE_SHARED(window_size, btsize, MPI_INFO_NULL, &
                                     process%node_comm, baseptr, window)
        if (process%node_rank /= 0) then
            call MPI_WIN_SHARED_QUERY(window, 0, window_size, disp_unit, &
                                      baseptr)
        end if

    end subroutine win_allocate

    ! temporary, for easier debugging
    pure function norma(a) result(v)
        real(kind=DP), intent(in) :: a(2,2)
        real(kind=DP) :: v
        v = a(1,1)**2 + a(2,1)**2 + a(2,2)**2
    end function norma

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

end module akde_botev