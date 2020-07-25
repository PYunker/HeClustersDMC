module system
    implicit none
contains

    function potencial(x) result(v) ! zatim jen zkusebne harmonicky oscilator, po kazde zmene je kod zapotrebi znovu zkompilovat

        !real, dimension(:) :: x
        real :: x(1)
        real :: v

        v = 0.5 * x(1) ** 2

        return

    end

end module system