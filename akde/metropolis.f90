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
! 28.02.2021 - Initial Version
!------------------------------------------------------------------------------
module metropolis
    
    implicit none
    private
    integer, parameter :: DP = kind(0.D0), separation = 100
    public :: draw_sample
    interface
    real(kind=DP) pure function rn2r(x)
        import :: DP
        real(kind=DP), intent(in) :: x(:)
    end function rn2r
    end interface

contains

    subroutine draw_sample(sample, sample_count, sample_dim, pdf)
        integer, intent(in) :: sample_count, sample_dim
        real(kind=DP), intent(out) :: sample(sample_dim, sample_count)
        integer :: i, j
        real(kind=DP) :: u(sample_dim, separation), v(separation), &
                         walker(sample_dim), proposed(sample_dim)
        procedure(rn2r) :: pdf
        
        do i=1,sample_count
            call random_number(v)
            call random_number(u)
            u = 0.61477929423D0 * log(u/(1-u)) ! it might be desirable &
            ! to replace exact logarithm with faster approximation
            do j=1,separation
                proposed = walker + u(:,j)
                if (pdf(proposed)/pdf(walker) > v(j)) walker = proposed
            end do
            sample(:,i) = walker
        end do

    end subroutine draw_sample

end module metropolis