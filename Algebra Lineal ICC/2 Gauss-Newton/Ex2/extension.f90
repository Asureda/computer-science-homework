module extension
    implicit none

    private :: gauss_jordan
    public :: invert_matrix
contains
    subroutine invert_matrix(a)
        double precision a(0:, 0:)
        double precision, allocatable :: b(:, :), a0(:, :)
        integer i, len_max, len_min
        len_min = lbound(a,1)
        len_max = ubound(a,1)
        allocate( a0(len_min:len_max,len_min:len_max), b(len_min:len_max, len_min:len_max) )
        a0 = a
        b = 0.0d0
        do i = len_min, len_max
            b(i, i) = 1.0d0
        end do
        do i = len_min, len_max
            a = a0
            call gauss_jordan(a, b(:,i), len_min, len_max)
        end do
        a(:, :) = b(:, :)
        deallocate(b)
    end subroutine

    subroutine gauss_jordan(a, b, len_min, len_max)
        integer, intent(in) :: len_min, len_max
        double precision,intent(inout) :: a(len_min:len_max,len_min:len_max), b(len_min:len_max)
        integer i, j, k
        double precision ar
        do k = len_min, len_max
            if(a(k, k) == 0.0d0) stop 'pivot = 0'
            ar = 1.0d0 / a(k, k)
            a(k, k) = 1.0d0
            do j = k + 1, len_max
                a(k, j) = ar * a(k, j)
            end do
            b(k) = ar * b(k)
            do i = len_min, len_max
                if(i /= k) then
                    do j = k + 1, len_max
                        a(i, j) = a(i, j) - a(i, k) * a(k, j)
                    end do
                    b(i) = b(i) - a(i,k) * b(k)
                    a(i, k) = 0.0d0
                end if
            end do
        end do
    end subroutine gauss_jordan

end module
