!------------------------------------------------------------------------------
! MODULE: kde
!
!> @author
!> Pavel Junker
!
! DESCRIPTION: 
!> implementation of vanilla kernel density estimation and Zdravko Botev's ada-
! ptive kernel density estimation algorithm (https://www.mathworks.com/matlabce
! ntral/fileexchange/58312-kernel-density-estimator-for-high-dimensions)
!
! REVISION HISTORY:
! 13.03.2021 - Initial Version
!------------------------------------------------------------------------------
module kde
	use omp_lib
    implicit none
    private

    integer, parameter :: DP = kind(0.D0), DP_SIZE = sizeof(0.0D0)

    real(kind=DP), parameter :: eps = epsilon(0.0D0)
	
    public :: botev_get_estimate, kde_get_estimate

contains
    
    subroutine botev_get_estimate(mu, w, sig, det, d, gam, n, lam, X, grid, &
								  grid_n, runs, estimate, h)
		integer, intent(in) :: d, gam, n, grid_n, runs
        real(kind=DP), intent(in) :: lam, grid(d,grid_n), X(d,n)
        real(kind=DP), intent(out) :: estimate(grid_n), h, mu(d,gam), w(gam), &
									  sig(d,d,gam), det(gam)
		real(kind=DP), external :: dlantr
		! locals
		integer, parameter :: ent_count = 100, max_iter = 1000
        integer :: i, j, k, idx(ent_count) = [(i, i=1,ent_count)]
        real(kind=DP) :: xRinv(d,n), p(n,gam), psig(n,gam), density(n), &
                         xSig(n), temp(grid_n), curv, ent_c, mesh(d,grid_n), &
						 density_old(ent_count), stopa(gam)      
        
        ! these act as parameters and are never changed
        real(kind=DP) :: norm_const, konst, lam2

        ! some neccessary initialization
        density(idx) = 1
        norm_const = 1d0/sqrt(8d0*atan(1d0))**d
        konst = 16d0*n*sqrt(atan(1d0))**2d0
        lam2 = sqrt(1/lam)**d

        ! it's assumed samples are ordered randomly
		mu = X(:,1:gam)
        h = 0.1d0*lam*n**(-real(d,kind=DP)/(d+4))
		sig = 0
        do i=1,d
            call random_number(sig(i,i,:))
			sig(i,i,:) = (0.1 + 0.9 * sig(i,i,:)) * h * lam**2
        end do
        call random_number(w)
        w = w / sum(w)
        det = 1
        ent_c = 10

		!$OMP PARALLEL PRIVATE(i,j,k,xRinv,xSig)
		k = 1
        main_loop : do

			!$OMP DO
	        do i=1,gam
	            call dpotrf('L', d, sig(:,:,i), d, j)
				det(i) = 1
	            do j=1,d
	                det(i) = det(i) * sig(j,j,i)
	            end do
	            call dtrtri('L','N',d,sig(:,:,i),d,j)
	            stopa(i) = dlantr('F','L','N',d,d,sig(:,:,i),d,curv)**2
			    do j=1,n
	                xRinv(:,j) = X(:,j) - mu(:,i)
	            end do
	            call dtrmm('L','L','N','N',d,n,1d0,sig(:,:,i),d,xRinv,d)
	            xSig = norm2(xRinv,dim=1)**2d0
	            p(:,i) = w(i) / det(i) * &
	                     exp(-0.5 * xSig - 0.5 * stopa(i) * h**2d0)
	            call dtrmm('L','L','T','N',d,n,lam2,sig(:,:,i),d,xRinv,d)
	            xSig = norm2(xRinv,dim=1)**2d0 ! + eps
	            psig(:,i) = p(:,i) * xSig
	        end do
			!$OMP END DO

			if (abs(ent_c - 1) < 1e-6 .or. k == max_iter) exit main_loop

			!$OMP BARRIER
			
			!$OMP SINGLE
	        ! ent_c = product(density(idx))
			density_old = density(idx)
	        density = sum(p,dim=2) + eps
			!$OMP END SINGLE
			!$OMP DO            
			do i=1,gam
	            p(:,i) = p(:,i) / density
	        end do
			!$OMP END DO

			!$OMP SINGLE
	        ! ent_c = ent_c / max(product(density(idx)),eps)
			ent_c = 1
			! dividing product(density(idx)) before and after updating
			! density does not cut it, since it tends to underflow for
			! large enough idx			
			do j=1,ent_count
				ent_c = ent_c * density(idx(j)) / density_old(j)
			end do
	        w = sum(p,dim=1)
			!$OMP END SINGLE
			!$OMP DO
	        do i=1,gam
				! adjust for possibility of w(i) being zero
				! openMP does not permit simple use of 'cycle'
				if (w(i) /= 0) then
					mu(:,i) = matmul(X,p(:,i)) / w(i)

			        do concurrent (j=1:n)
			            xRinv(:,j) = X(:,j) - mu(:,i)
			        end do
			        do concurrent (j=1:d)
			            xRinv(j,:) = xRinv(j,:) * sqrt(p(:,i))
			        end do

				    sig(:,:,i) = matmul(xRinv, transpose(xRinv))

				    sig(:,:,i) = sig(:,:,i) / w(i)
				    do j=1,d
				        sig(j,j,i) = sig(j,j,i) + h**2
				    end do
				end if
	        end do
			!$OMP END DO

			!$OMP SINGLE
	        xSig = sum(psig,dim=2)/density        
	        curv = sum(xSig)/n
	        h = (curv*konst)**(-real(1,kind=DP)/(d+2))
	        w = w / sum(w)
	        !$OMP END SINGLE
	        k = k + 1

		end do main_loop
		!$OMP MASTER
		write(*,"('done, k : ', i0)") k
		!$OMP END MASTER
		!$OMP END PARALLEL

        ! actually obtaining the estimate
        estimate = 0
		!$OMP PARALLEL DO PRIVATE(j,mesh,temp) REDUCTION(+ : estimate)
        do i=1,gam
            do concurrent (j=1:grid_n)
                mesh(:,j) = grid(:,j) - mu(:,i)
            end do
            call dtrmm('L','L','N','N',d,grid_n,1d0,sig(:,:,i),d,mesh,d)
            temp = -0.5d0 * norm2(mesh,dim=1)**2d0
            estimate = estimate + w(i) * norm_const / det(i) * exp(temp)
        end do
		!$OMP END PARALLEL DO

		!end associate

    end subroutine botev_get_estimate

    ! vanilla KDE with fixed bandwidth isotropic kernel
    subroutine kde_get_estimate(d, n, X, grid, grid_n, estimate, gam, bandwidth)
        integer, intent(in) :: d, n, grid_n
        integer, optional, intent(in) :: gam
        real(kind=DP), intent(in) :: grid(d,grid_n), X(d,n)
        real(kind=DP), optional, intent(in) :: bandwidth
        real(kind=DP), intent(out) :: estimate(grid_n)
		! locals
        real(kind=DP), parameter :: konst = 1d0/sqrt(8d0*atan(1d0))
        real(kind=DP) :: h, det, temp(grid_n), mesh(d,grid_n)
        integer :: num_kernels, i, j

        if(present(gam)) then
            num_kernels = gam
        else
            num_kernels = n
        end if

        if(present(bandwidth)) then
            h = bandwidth
        else
            h = (4.0D0/(d+2)/n)**(1.0D0/(d+4))
        end if

        det = konst**d/n/h**d

        estimate = 0
		!$OMP PARALLEL DO PRIVATE(j,mesh,temp) REDUCTION(+ : estimate)
		do i=1,num_kernels
            do concurrent(j=1:grid_n)
                mesh(:,j) = grid(:,j) - X(:,i)
            end do
        	temp = -0.5d0 * (norm2(mesh,dim=1)/h) ** 2
        	estimate = estimate + det * exp(temp)
        end do
		!$OMP END PARALLEL DO

    end subroutine kde_get_estimate

end module kde
