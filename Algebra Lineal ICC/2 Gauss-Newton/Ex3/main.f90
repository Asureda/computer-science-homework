program main
    use extension
    implicit none
    integer,parameter :: fi = 10, num = 1
    integer :: i, j, jmax,io,n,k,l
    double precision ::  a , b,a_old,b_old,alpha = 1.0d0, error, r_sum,err_old,sd,x1,x2,t1,t2,t3,t4,r_sq
    double precision, DIMENSION(:), ALLOCATABLE :: x(:),y(:),delta(:),err(:),rate(:),parameter1(:),parameter2(:)
    double precision, DIMENSION(:,:), ALLOCATABLE :: jacob(:,:),jacob_t(:,:),d(:,:),var(:,:),hessian_approx(:,:),Bmat(:,:),Bm(:,:)
    double precision,parameter :: epsilon = 1.0d-6
    integer :: count = 0 ! count=Gauss-Newton
    !parameters
    a=1.d0
    b=1.d0
    jmax = 30


  OPEN (10,FILE='data.txt')
    n = 0
    DO
      READ(10,*,iostat=io)
      IF (io/=0) EXIT
      n = n + 1
    END DO
    print*,n
close(10)

  allocate(x(n))
  allocate(y(n))
  allocate(rate(jmax))
  allocate(err(jmax))
  allocate(parameter1(jmax))
  allocate(parameter2(jmax))
  allocate(delta(num))
  allocate(jacob(n,num))
  allocate(jacob_t(num,n))
  allocate(var(num,num))
  allocate(d(num,num))
  allocate(Bmat(num,num))
  allocate(Bm(num,num))
  allocate(hessian_approx(num,num))

    do
        write (*, '(a)', advance='no') "alpha?(default=1) :"
        read (*, *) alpha
        if (1.0d0 < alpha .or. alpha <= 0.0d0) then
            write (*, *) "alpha value must be between 0< and <=1"
        else
            exit
        end if
    end do

    open(fi, file='data.txt')
    do i = 1, n
        read (fi, *) x(i), y(i)
    end do
    close(fi)

    a_old = a
  !  b_old = b

    do j = 1, jmax
      error = 0.d0
      r_sum = 0.d0


      r_sum =sum(r(x,y,a,b))
      error = error + r_sum*r_sum

      print*,'r square',error
        do i = 1, n
            jacob(i, 1) = df_da(x(i), a, b)
            !jacob(i, 2) = df_db(x(i), a, b)
        end do

        jacob_t = transpose(jacob)
        d = matmul(jacob_t, jacob)
        hessian_approx = d
        call invert_matrix(d)
        delta = -matmul(matmul(d, jacob_t), r(x, y, a, b))

        a = a_old + alpha * delta(1)
        !b = b_old + alpha * delta(2)

        err(j) = a-a_old

        count = count + 1

        a_old = a
        !b_old = b
        parameter1(j) = a
      !  parameter2(j) = b

        if (delta(1) < epsilon ) then
            go to 4
        end if


    end do

    4 continue

    if (count == jmax) then
        write (*, *) "Calculation loop exceeds the limit. Result may not be correct. Try change alpha lower."
    end if

    print*,'----------------------------------------------------------------'
    print*,'Parameters'
    print*,'----------------------------------------------------------------'

    write (*, '(a, i0)') "Iteration count =", count
    write (*, '(2(a, f15.5))') "a=", a, ", b=", b
    !
    ! !
     do i = 1, n
         write(*,*) x(i), y(i), f(x(i),a,b)
     end do

    print*,'----------------------------------------------------------------'
    print*,'Errors'
    print*,'----------------------------------------------------------------'

    do i = 1, n
        write(*,*) y(i) - abs(f(x(i),a,b))
    end do

    do i = 1, n
        t1 = t1 + y(i)*dexp(-b/x(i))
        t2 = t2 + dexp(-2*b*x(i))
        t3 = t3 + y(i)*y(i)
        t4 = t4 + y(i)*dexp(-b/x(i))
    end do

      r_sq = sum(r(x,y,a,b)*r(x,y,a,b))

    x1 = t1/t2
    x2 = t3 - x1*t4


    print*,'x1',x1
    print*,'sum r',x2,r_sq




contains

    double precision function f(x, a, b)
        double precision,intent(in) :: x, a, b
        !f = a*cosh(b*x)
        f = a*exp(-b*x)
    end function

    double precision function df_da(x, a, b)
        double precision,intent(in) :: x, a, b
        !df_da = -cosh(b*x)
        df_da = exp(-b*x)
    end function

    double precision function df_db(x, a, b)
        double precision,intent(in) :: x, a, b
        !df_db = -a*x*sinh(b*x)
        df_db = -a*x*exp(-b*x)
    end function

    double precision function df_dada(x, a, b)
        double precision,intent(in) :: x, a, b
        !df_da = -cosh(b*x)
        df_dada = 0.d0
    end function

    double precision function df_dadb(x, a, b)
        double precision,intent(in) :: x, a, b
        !df_da = -cosh(b*x)
        df_dadb = -x*exp(-b*x)
    end function

    double precision function df_dbdb(x, a, b)
        double precision,intent(in) :: x, a, b
        !df_da = -cosh(b*x)
        df_dbdb = -x*exp(-b*x)
    end function

    double precision function df_dbda(x, a, b)
        double precision,intent(in) :: x, a, b
        !df_da = -cosh(b*x)
        df_dbda = a*x*x*exp(-b*x)
    end function



    function r(x, y, a, b)
        double precision,intent(in) :: x(:), y(:), a, b
        double precision r(size(x)),error
        integer i
        do i = 1, size(x)
            r(i) = y(i) - f(x(i), a, b)
        end do

    end function


    subroutine display(matrix)
        double precision,intent(in) :: matrix(:, :)
        integer i
        integer imax
        imax = ubound(matrix, 1)
        do i = 1, imax
            write (*, '(100(F0.3, x))') matrix(i, :)
        end do
    end subroutine
end program
