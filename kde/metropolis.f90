!------------------------------------------------------------------------------
! MODULE: metropolis
!
!> @author
!> Pavel Junker
!
! DESCRIPTION: 
!>  a simple implementation of Metropolis-Hastings sampling algorithm
!
! REVISION HISTORY:
! 13.03.2021 - Initial Version
!------------------------------------------------------------------------------
module metropolis
    
    implicit none
    private
    integer, parameter :: DP = kind(0.D0), separation = 100

    abstract interface
    real(kind=DP) pure function rn_to_r(x)
        import :: DP
        real(kind=DP), intent(in) :: x(:)
    end function rn_to_r
    end interface
	
	public :: draw_sample, rn_to_r

contains

    subroutine draw_sample(sample, sample_count, sample_dim, pdf)
        integer, intent(in) :: sample_count, sample_dim
        real(kind=DP), intent(out) :: sample(sample_dim, sample_count)
        integer :: i, j
        real(kind=DP) :: u(sample_dim, separation), v(separation), &
                         walker(sample_dim), proposed(sample_dim)
        procedure(rn_to_r) :: pdf
        
        do i=1,sample_count
            call random_number(v)
            call random_number(u)
			! tranformation below generates rougly normal sample from uniform
			! one. it is in fact needlesly precise, so approximating the log-
			! arithm for higher performance should be OK
            u = 0.61477929423D0 * log(u/(1-u))
            do j=1,separation
                proposed = walker + u(:,j)
                if (pdf(proposed)/pdf(walker) > v(j)) walker = proposed
            end do
            sample(:,i) = walker
        end do

    end subroutine draw_sample

end module metropolis
