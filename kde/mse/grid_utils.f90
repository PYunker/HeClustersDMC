!------------------------------------------------------------------------------
! MODULE: grid_utils
!
!> @author
!> Pavel Junker
!
! DESCRIPTION: 
!> basically just Fortran analogs of Octave/Matlab 'linspace' and 'meshgrid'
!
! REVISION HISTORY:
! 28.02.2021 - Initial Version
!------------------------------------------------------------------------------
module grid_utils
    implicit none
    private
    integer, parameter :: DP = kind(0.D0)
    public :: grid_init_1d, grid_init_2d
contains
    
    pure subroutine grid_init_1d(grid, low, high, n)
        integer, intent(in) :: n
        real(kind=DP), intent(out) :: grid(1,n)
        real(kind=DP), intent(in) :: low, high
        real(kind=DP) :: a, b
        integer :: i
        
        ! obtain axis by linearly transforming sequence 1, 2, ..., n
        a = (high - low)/(n - 1)
        b = (low*n - high)/(n - 1)
        grid(1,:) = a * [(i, i=1,n)] + b

    end subroutine grid_init_1d

    pure subroutine grid_init_2d(grid, x_low, x_high, y_low, y_high, nx, ny)
        integer, intent(in) :: nx, ny
        real(kind=DP), intent(out) :: grid(2,nx*ny)
        real(kind=DP), intent(in) :: x_low, x_high, y_low, y_high
        real(kind=DP) :: a, b
        integer :: i, j, k

        k = 1
        do i=1,nx
            do j=1,ny
                grid(:,k) = [i,j]
                k = k + 1
            end do
        end do
        ! transform x axis
        a = (x_high - x_low)/(nx - 1)
        b = (x_low*nx - x_high)/(nx - 1)
        grid(1,:) = a * grid(1,:) + b
        ! transform y axis
        a = (y_high - y_low)/(ny - 1)
        b = (y_low*ny - y_high)/(ny - 1)
        grid(2,:) = a * grid(2,:) + b
    
    end subroutine grid_init_2d

end module grid_utils