PROGRAM Heat_Eq_TD_1D
    
    IMPLICIT NONE

    INTEGER ,PARAMETER :: Nx = 10, Nt = 400
    REAL(KIND = 8) ,PARAMETER :: deltaX = 1.d0 / (1 + Nx), deltaT = 1.d0 / Nt
    REAL(KIND = 8) ,PARAMETER :: coef = deltaT/(deltaX*deltaX)
    REAL(KIND = 8) :: t
    REAL(KIND = 8), DIMENSION(0:Nx + 1) :: old_solution, solution, x
    REAL(KIND = 8) :: error = 0.d0
    INTEGER ::  i

    DO i = 0, Nx + 1
        x(i) = i * deltaX
    END DO

    old_solution = 0.d0
    t = 0.d0

    DO WHILE ( t < 1.d0)
        t = t + deltaT
        solution(0) = 0.d0
        solution(1:Nx) = old_solution(1:Nx) * (1.d0 - 2.d0 * coef ) + coef*(old_solution(0:Nx - 1) + old_solution(2:Nx + 1)) + deltaT * sin(x(1:Nx))*(sin(t)+cos(t))
        solution(Nx + 1) = sin(1.d0) * sin(t)
        old_solution = solution

    END DO

    error = NORM2(solution - sin(x)*sin(t))/(Nx * Nt)

    WRITE(*,*) error

END






    
