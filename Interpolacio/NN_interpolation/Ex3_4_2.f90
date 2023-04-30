! -------------------------------------------------------------------
! NN_interpolation
! -------------------------------------------------------------------
program NN_interpolation
implicit none
integer                                 :: i, j, k, l, m, n, l_int
real                                    :: start, finish
real, dimension(:), allocatable         :: x, y, xInt, yInt
real                                    :: xInit, xFin, pas, pas_int

call cpu_time(start)

l = 20; xInit = 0; xFin = 10
allocate(x(l), y(l))

pas = (xFin - xInit)/l

do i = 1, l, 1
        x(i) = pas*i
        y(i) = cos(x(i))
end do


l_int = 101
pas_int = (xFin - xInit)/l_int
allocate(xInt(l_int), yInt(l_int))

do i = 1, l_int, 1
        xInt(i) = pas_int*i
end do


do i = 1, l_int, 1
        n = int(0.5 + xInt(i)/pas)
        write(6,*) n
        yInt(i) = y(n)
end do

open(100,file="interp.log")
open(200,file="known.log")

do i = 1, l, 1
        write(200,*) x(i), y(i)
enddo
do i = 1, l_int, 1
        write(100,*) xInt(i), yInt(i)
enddo

deallocate(x, y, xInt, yInt)
call cpu_time(finish)
print*, "CPU TIME: ", finish - start
end program NN_interpolation
