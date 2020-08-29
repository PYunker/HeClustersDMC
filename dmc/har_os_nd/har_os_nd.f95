program har_os_1d

    use dmc
    
    implicit none
    ! ================= metaparametry =====================
    integer :: p_bloku, p_kroku, N0, N_max, N_min, i
    integer, parameter :: n = 2
    real :: dt, alfa, Er, x0(2,1), xi(2,1)
    real :: Kmat(2,2), Q(2,2), L(2,2)

    Kmat = reshape([4,-3,-3,4], shape = [2,2])

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

    call simulace(p_bloku, p_kroku, N0, N_max, N_min, alfa, dt, [x0], Er, potencial) ! samotna simulace
    
    call jacobi(Kmat, L, Q, 1.0e-09, 2)
    Q = transpose(Q)
    open(1, file="p_tran")
    open(2, file="populace_dmc")
    do
        read(2,*, end = 10) xi
        write(1,*) matmul(Q,xi)
    end do
10  continue
    close(2)
    close(1)

    print*, "==================================="
    write(*,*) (L(i,i), i=1,2) ! vypise vlastni hodnoty K

    contains

    function potencial(x) result(v)

        real :: x(n), v
        v = 0.5 * dot_product(x,[matmul(Kmat,reshape(x, shape = [n,1]))])
        !v = 0.5 * x(1) ** 2 + 0.5 * x(2) ** 2
        return

    end

    subroutine jacobi(a,L,E,eps,n)
        ! kod z https://ww2.odu.edu/~agodunov/computing/programs/book2/Ch07/Jacobi.f90
        integer :: i, j, k
        integer, intent(in) :: n ! rozmer A, musi byt bohuzel zadan zvlast
        real :: md, mdp, t, kfc, c, s, cs, sc
        real, intent(in) :: a(n,n), eps
        real, intent(out) :: L(n,n), E(n,n)

        L = a ! 'L' je (na konci algoritmu) diagonalni matice podobna 'A'
        
        E = 0.0
        do i=1,n
          E(i,i) = 1.0
        end do
        
        ! vypocte soucet ctvercu prvku mimo diagonalu
        md = 0.0
        do i=1,n
          do j=1,n
            if (i.ne.j) md = md + a(i,j)**2
          end do
        end do
        
        if (md <= eps) return
        
        ! vypocte prumer 'md'
        mdp = 0.5*md/real(n*n)
        
        do while (md.gt.eps)
          do i=1,n-1
            do j=i+1,n
              if (L(j,i)**2 <= mdp) cycle  ! zacne s vetsimi odchylkami
              md = md - 2.0*L(j,i)**2
              mdp = 0.5*md/real(n*n)
        ! vypocte parametry rotace
              t = (L(j,j)-L(i,i))/(2.0*L(j,i))
              kfc = 0.5*t/sqrt(1.0+t**2)
              s = sqrt(max(0.5+kfc,0.0))
              c = sqrt(max(0.5-kfc,0.0))
        ! L -> L * G^T
              do k=1,n
                cs =  c*L(i,k)+s*L(j,k)
                sc = -s*L(i,k)+c*L(j,k)
                L(i,k) = cs
                L(j,k) = sc
              end do
        ! L -> G * L & E -> G * E
              do k=1,n
                cs =  c*L(k,i)+s*L(k,j)
                sc = -s*L(k,i)+c*L(k,j)
                L(k,i) = cs
                L(k,j) = sc
                cs =  c*E(k,i)+s*E(k,j)
                sc = -s*E(k,i)+c*E(k,j)
                E(k,i) = cs
                E(k,j) = sc
              end do
            end do
          end do
        end do
        return

    end subroutine

end program har_os_1d