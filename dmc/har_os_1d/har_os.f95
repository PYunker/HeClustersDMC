program har_os_1d

    use dmc
    
    implicit none
    ! ================= metaparametry =====================
    integer :: p_bloku, p_kroku, N0, N_max, N_min
    real :: dt, alfa, Er, x0(1)
    
    print*, "==================================="
    write(*,'(a14)',advance='no') "pocet bloku : "
    read(*,*) p_bloku
    write(*,'(a14)',advance='no') "pocet kroku : "
    read(*,*) p_kroku
    write(*,'(a14)',advance='no') "N_0 : " ! pocatecni pocet replik
    read(*,*) N0
    write(*,'(a14)',advance='no') "N_max : " ! jejich maximalni pocet (je-li prekrocen, program spadne, do budoucna by to chtelo osetrit)
    read(*,*) N_max
    write(*,'(a14)',advance='no') "N_min : " ! jejich minimalni pocet (zatim nic nedala)
    read(*,*) N_min
    write(*,'(a14)',advance='no') "dt : " ! casovy krok
    read(*,*) dt
    write(*,'(a14)',advance='no') "alfa : " ! 'pruznost' referencni enrgie - vetsi hodnota, vetsi opravy
    read(*,*) alfa
    write(*,'(a14)',advance='no') "x0 : " ! pocatecni poloha walkeru
    read(*,*) x0
    print*, "==================================="
    ! =====================================================

    call simulace(p_bloku, p_kroku, N0, N_max, N_min, alfa, dt, x0, Er, potencial) ! samotna simulace

    contains

    function potencial(x) result(v)

        !real, dimension(:) :: x
        real :: x(1)
        real :: v

        v = 2 * x(1) ** 2

        return

    end

    subroutine jacobi(a,l,x,abserr,n)

        integer :: i, j, k
        integer, intent(in) :: n
        real :: b2, bar, beta, coeff, c, s, cs, sc
        real, intent(in) :: a(n,n), abserr
        real, intent(out) :: l(n,n), x(n,n)

        l = a
        
        ! initialize x(i,j)=0, x(i,i)=1
        ! *** the array operation x=0.0 is specific for Fortran 90/95
        x = 0.0
        do i=1,n
          x(i,i) = 1.0
        end do
        
        ! find the sum of all off-diagonal elements (squared)
        b2 = 0.0
        do i=1,n
          do j=1,n
            if (i.ne.j) b2 = b2 + a(i,j)**2
          end do
        end do
        
        if (b2 <= abserr) return
        
        ! average for off-diagonal elements /2
        bar = 0.5*b2/real(n*n)
        
        do while (b2.gt.abserr)
          do i=1,n-1
            do j=i+1,n
              if (l(j,i)**2 <= bar) cycle  ! do not touch small elements
              b2 = b2 - 2.0*l(j,i)**2
              bar = 0.5*b2/real(n*n)
        ! calculate coefficient c and s for Givens matrix
              beta = (l(j,j)-l(i,i))/(2.0*l(j,i))
              coeff = 0.5*beta/sqrt(1.0+beta**2)
              s = sqrt(max(0.5+coeff,0.0))
              c = sqrt(max(0.5-coeff,0.0))
        ! recalculate rows i and j
              do k=1,n
                cs =  c*l(i,k)+s*l(j,k)
                sc = -s*l(i,k)+c*l(j,k)
                l(i,k) = cs
                l(j,k) = sc
              end do
        ! new matrix a_{k+1} from a_{k}, and eigenvectors 
              do k=1,n
                cs =  c*l(k,i)+s*l(k,j)
                sc = -s*l(k,i)+c*l(k,j)
                l(k,i) = cs
                l(k,j) = sc
                cs =  c*x(k,i)+s*x(k,j)
                sc = -s*x(k,i)+c*x(k,j)
                x(k,i) = cs
                x(k,j) = sc
              end do
            end do
          end do
        end do
        return

    end subroutine

end program har_os_1d
