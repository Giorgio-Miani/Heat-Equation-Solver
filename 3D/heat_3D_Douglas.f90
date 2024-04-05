
PROGRAM Heat_Eq_TD_1D
   IMPLICIT NONE

   INTEGER ,PARAMETER :: Nx = 300, Nt = 100
   INTEGER, PARAMETER ::dim = Nx*Nx*Nx
   INTEGER, EXTERNAL :: findex
   REAL(KIND = 8) ,PARAMETER :: deltax = 1.d0/ (1 + Nx), deltat = 1.d0 / Nt, deltax2 = deltaX*deltaX
   REAL(KIND = 8) ,PARAMETER :: diag = 1.d0 +  deltaT/(deltaX*deltaX)
   REAL(KIND = 8), PARAMETER :: extra_diag = -deltaT/(2*deltaX*deltaX)
   REAL(KIND = 8) :: t
   REAL(KIND = 8) :: w
   REAL(KIND = 8), DIMENSION(1:dim) :: rhs, old_solution, solution, b, real_solution
   real(KIND = 8) :: a, sin_1
   REAL(KIND = 8) :: error = 0.d0
   real(kind = 8), DIMENSION(Nx) :: xsin
   REAL :: start, finish, startx, finishx, totalx, starty, finishy, totaly, startz, finishz, totalz
   INTEGER ::  i
   INTEGER :: j
   INTEGER :: k


   DO i = 1, Nx
      xsin(i) = sin(i * deltaX)
   END DO
   sin_1 = sin(1.d0)

   old_solution = 0.d0
   a = extra_diag
   t = 0.d0 + deltaT

   call cpu_time(start)

   DO WHILE ( t <= 1.d0)

      call cpu_time(startx)

      !RHS
      !THOMAS ON X
      b = diag
      !Forcing Term
      do k = 1, Nx
         do j = 1, Nx
            do i = 1, Nx
               rhs(findex(i, j, k, Nx)) = deltaT*((xsin(i)*xsin(j)*xsin(k))*(3*sin(t + deltaT/2) + cos(t + deltaT/2)) - (1/(deltax2)*(6 * old_solution(findex(i, j, k, Nx)))))
            end do
         end do
      end do

      !BCS on X axis
      do k = 1, Nx
         do j = 1, Nx
            rhs(findex(Nx, j, k, Nx)) = rhs(findex(Nx, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(Nx - 1, j, k, Nx))+ deltaT*(1/(deltax2))*(sin_1*xsin(j)*xsin(k)*sin(t))
            rhs(findex(1, j, k, Nx)) = rhs(findex(1, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(2, j, k, Nx))
            do i = 2, Nx - 1
               rhs(findex(i, j, k, Nx)) = rhs(findex(i, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i - 1, j, k, Nx))
               rhs(findex(i, j, k, Nx)) = rhs(findex(i, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i + 1, j, k, Nx))
            end do
         end do
      end do

      !BCS on Y axis
      do k = 1, Nx
         do i = 1, Nx
            rhs(findex(i, Nx, k, Nx)) = rhs(findex(i, Nx, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i, Nx - 1, k, Nx)) + deltaT*(1/(deltax2))*(sin_1*xsin(i)*xsin(k)*sin(t))
            rhs(findex(i, 1, k, Nx)) = rhs(findex(i, 1, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i, 2, k, Nx))
            do j = 2, Nx - 1
               rhs(findex(i, j, k, Nx)) = rhs(findex(i, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i, j - 1, k, Nx))
               rhs(findex(i, j, k, Nx)) = rhs(findex(i, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i, j + 1, k, Nx))
            end do
         end do
      end do

      !BCS on Z axis
      do j = 1, Nx
         do i = 1, Nx
            rhs(findex(i, j, Nx, Nx)) = rhs(findex(i, j, Nx, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i, j, Nx - 1, Nx)) + deltaT*(1/(deltax2))*(sin_1*xsin(j)*xsin(i)*sin(t))
            rhs(findex(i, j, 1, Nx)) = rhs(findex(i, j, 1, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i, j, 2, Nx))
            do k = 2, Nx - 1
               rhs(findex(i, j, k, Nx)) = rhs(findex(i, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i, j, k - 1, Nx))
               rhs(findex(i, j, k, Nx)) = rhs(findex(i, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i, j, k + 1, Nx))
            end do
         end do
      end do

      
      ! do k = 1, Nx
      !    do j = 1, Nx
      !       do i = 1, Nx
      !          write(*, '(F15.9)', advance='no') rhs(findex(i, j, k, Nx))
      !       end do
      !       print *
      !    end do
      !    print *
      ! end do




      do k = 1, Nx
         do j = 1, Nx
            !BCS
            rhs(findex(Nx, j, k, Nx)) = rhs(findex(Nx, j, k, Nx)) - extra_diag*sin_1*xsin(j)*xsin(k)*(sin(t + deltaT) - sin(t))
            do i = 2, Nx
               w = a/b(findex(i - 1, j, k, Nx))
               rhs(findex(i, j, k ,Nx)) = rhs(findex(i, j, k ,Nx)) - w*rhs(findex(i - 1, j, k, Nx))
               b(findex(i, j, k, Nx)) = b(findex(i, j, k, Nx)) - w*a
            end do
         end do
      end do

      do k = 1, Nx
         do j = 1, Nx
            solution(findex(Nx, j, k, Nx)) = rhs(findex(Nx, j, k, Nx))/b(findex(Nx, j, k, Nx))
            do i = Nx - 1,1, -1
               solution(findex(i, j, k, Nx)) = (rhs(findex(i, j, k, Nx)) - a * solution(findex(i + 1, j, k, Nx)))/b(findex(i, j, k, Nx))
            end do
         end do
      end do


      call cpu_time(finishx)

      totalx = totalx + finishx - startx

      call cpu_time(starty)

      !THOMAS ON Y
      rhs = solution
      b = diag
      do k = 1, Nx
         do j = 2 , Nx
            do i = 1, Nx
               if(j == Nx) then
                  !BCS
                  rhs(findex(i, Nx, k, Nx)) = rhs(findex(i, Nx, k, Nx)) - extra_diag*sin_1*xsin(i)*xsin(k)*(sin(t + deltaT) - sin(t))
               end if
               w = a/b(findex(i, j - 1, k, Nx))
               rhs(findex(i, j, k ,Nx)) = rhs(findex(i, j, k ,Nx)) - w*rhs(findex(i, j - 1, k, Nx))
               b(findex(i, j, k, Nx)) = b(findex(i, j, k, Nx)) - w*a
            end do
         end do
      end do

      do k = 1, Nx
         do j = Nx - 1, 1, -1
            do i = 1, Nx
               solution(findex(i, Nx, k, Nx)) = rhs(findex(i, Nx, k, Nx))/b(findex(i, Nx, k, Nx))
               solution(findex(i, j, k, Nx)) = (rhs(findex(i, j, k, Nx)) - a * solution(findex(i, j + 1, k, Nx)))/b(findex(i, j, k, Nx))
            end do
         end do
      end do

      call cpu_time(finishy)

      totaly = totaly + finishy - starty

      call cpu_time(startz)
      !THOMAS ON Z
      rhs = solution
      b = diag
      do k = 2, Nx
         do j = 1, Nx
            do i = 1, Nx
               if ( k == Nx ) then
                  !BCS
                  rhs(findex(i, j,Nx, Nx)) = rhs(findex(i, j, Nx, Nx)) - extra_diag*sin_1*xsin(j)*xsin(i)*(sin(t + deltaT) - sin(t))
               end if
               w = a/b(findex(i, j,k - 1, Nx))
               rhs(findex(i, j, k ,Nx)) = rhs(findex(i, j, k ,Nx)) - w*rhs(findex(i, j, k - 1, Nx))
               b(findex(i, j, k, Nx)) = b(findex(i, j, k, Nx)) - w*a
            end do
         end do
      end do


      do k = Nx -1, 1, -1
         do j = 1, Nx
            do i = 1, Nx
               solution(findex(i, j, Nx, Nx)) = rhs(findex(i, j, Nx, Nx))/b(findex(i, j, Nx, Nx))
               solution(findex(i, j, k, Nx)) = (rhs(findex(i, j, k, Nx)) - a * solution(findex(i, j, k + 1, Nx)))/b(findex(i, j, k, Nx))
            end do
         end do
      end do



      old_solution = solution + old_solution

      t = t + deltaT

      call cpu_time(finishz)

      totalz = totalz + finishz - startz

   END DO

      ! !print rhs
      ! do k = 1, Nx
      !    do j = 1, Nx
      !       do i = 1, Nx
      !          write(*, '(F15.9)', advance='no') old_solution(findex(i, j, k, Nx))
      !       end do
      !       print *
      !    end do
      !    print *
      ! end do


   call cpu_time(finish)

   Write (*,*) "Time x axis: "
   Write (*,*) totalx

   Write (*,*) "Time y axis: "
   Write (*,*) totaly

   Write (*,*) "Time z axis: "
   Write (*,*) totalz

   Write (*,*) "Total time: "
   Write (*,*) finish - start

   do k = 1, Nx
      do j = 1, Nx
         do i = 1, Nx
            real_solution(findex(i, j, k, Nx)) = xsin(i)*xsin(j)*xsin(k)*sin_1
         end do
      end do
   end do

   error = NORM2(old_solution - real_solution)/dim

   WRITE(*,*) "Error of "
   WRITE(*,*) error

END

integer function findex(o,p,q, N)
   integer, intent (in) :: o, p, q , N
   findex = o + (p - 1) * N + (q - 1) * N * N
end function findex
