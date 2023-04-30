! -------------------------------------------------------------------
! linear_interpolation
! -------------------------------------------------------------------
program linear_interpolation
implicit none
integer					:: i, j, k, l, m, n, l_int
real					:: start, finish
real, dimension(:), allocatable		:: x, y, xInt, yInt
real                                    :: xInit, xFin
real                                    :: pas, pas_int
integer                                 :: i1, i2

call cpu_time(start)
xInit = 0; xFin = 10

l = 20
allocate(x(l), y(l))
pas = (xFin - xInit)/l
do i = 1, l, 1
        x(i) = pas*i
        y(i) = cos(x(i))
end do

l_int = 100
allocate(xInt(l_int), yInt(l_int))
pas_int = (xFin - xInit)/l_int

do i = 1, l_int, 1
        xInt(i) = pas_int*i
        i1 = int(xInt(i)/pas)
        i2 = int(xInt(i)/pas) + 1
        yInt(i) = y(i1) + ((xInt(i) - x(i1))/(x(i2) - x(i1)))*(y(i2) - y(i1))
end do

open(1,file="known.log")
open(2,file="interpolated.log")
do i = 1, l_int, 1
        write(2,*) xInt(i), yInt(i)
end do
do i = 1, l, 1
        write(1,*) x(i), y(i)
end do


call cpu_time(finish)
print*, "CPU TIME: ", finish - start
end program linear_interpolation
