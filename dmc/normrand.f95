module normrand
  implicit none
contains ! generator nahodnych cisel stazeny s internetu, ten bych pozdeji zkusil dodelat 'na miru' pro toto pouziti
  FUNCTION grnd() RESULT (ran_norm) ! vraci nahodne cislo z normalnih rozdeleni

      IMPLICIT NONE
      REAL :: ran_norm
      
      !     Local variables
      REAL, PARAMETER :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,   &
                        half = 0.5, r1 = 0.27597, r2 = 0.27846
      REAL            :: u, v, x, y, q
      
      !     Generate P = (u,v) uniform in rectangle enclosing acceptance region
      
      DO
        CALL RANDOM_NUMBER(u)
        CALL RANDOM_NUMBER(v)
        v = 1.7156 * (v - half)
      
      !     Evaluate the quadratic form
        x = u - s
        y = ABS(v) - t
        q = x**2 + y*(a*y - b*x)
      
      !     Accept P if inside inner ellipse
        IF (q < r1) EXIT
      !     Reject P if outside outer ellipse
        IF (q > r2) CYCLE
      !     Reject P if outside acceptance region
        IF (v**2 < -4.0*LOG(u)*u**2) EXIT
      END DO
      
      !     Return ratio of P's coordinates as the normal deviate
      ran_norm = v/u
      RETURN
      
  END FUNCTION grnd

  SUBROUTINE gasdev_s(harvest)
      ! Numerical Recipes routine for generating a single normal random deviate,
      ! adapted to use the compiler's random number generator.
      
      IMPLICIT NONE
      REAL, INTENT(OUT) :: harvest
      
      ! Local variables
      REAL          :: rsq, v1, v2
      REAL, SAVE    :: g
      LOGICAL, SAVE :: gaus_stored = .false.
      
      IF (gaus_stored) THEN
      harvest = g
      gaus_stored = .false.
      ELSE
      DO
          CALL RANDOM_NUMBER(v1)
          CALL RANDOM_NUMBER(v2)
          v1 = 2.0*v1 - 1.0
          v2 = 2.0*v2 - 1.0
          rsq = v1**2 + v2**2
          if (rsq > 0.0 .and. rsq < 1.0) EXIT
      END DO
      rsq = SQRT(-2.0*LOG(rsq)/rsq)
      harvest = v1*rsq
      g = v2*rsq
      gaus_stored = .true.
      END IF
      
      RETURN
  END SUBROUTINE gasdev_s

! vrati vektor (1,n) nahodnych cisel z normaliho rozdeleni, rozhodne potrebuje prepsat
  function normrnd(n) result(x)
    integer :: n, i
    real, dimension(n) :: x
    
    do i = 1,n
    x(i) = grnd()
    end do
  end function normrnd

  function randf(m,n) result(v)
    integer :: m, n
    real, dimension(m,n) :: v
    call random_number(v)
    return
  end function randf

end module normrand