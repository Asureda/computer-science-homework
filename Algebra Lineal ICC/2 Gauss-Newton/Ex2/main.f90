program main
    use extension
    implicit none
    integer,parameter :: fi = 10, num = 2
    integer :: i, j, jmax,io,n,k,l
    double precision :: a , b,a_old,b_old,alpha = 1.0d-6, c, error, r_sum,err_old,sd,suma1,suma2,suma3,suma4,alpha1,alpha2
    double precision, DIMENSION(:), ALLOCATABLE :: x(:),y(:),delta(:)
    double precision, DIMENSION(:,:), ALLOCATABLE :: jacob(:,:),jacob_t(:,:),d(:,:),d2(:,:),jacob_old(:,:),d_old(:,:)
    double precision,parameter :: epsilon = 0.0001
    integer :: count = 0 ! count=Gauss-Newton
    !parameters
    a=-0.1d0
    b=2400.d0
    jmax = 100
    alpha1=alpha
    alpha2=alpha

  OPEN (10,FILE='data.txt')
    n = 0
    DO
      READ(10,*,iostat=io)
      IF (io/=0) EXIT
      n = n + 1
    END DO
    print*,n
close(10)

  open(11,file='Resposta.dat')
  allocate(x(n))
  allocate(y(n))
  allocate(delta(num))
  allocate(jacob(n,num))
  allocate(jacob_old(n,num))
  allocate(jacob_t(num,n))
  allocate(d(num,num))
  allocate(d_old(num,num))
  allocate(d2(num,num))





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
        x(i) = x(i)
    end do
    close(fi)



    do j = 1, jmax

        do i = 1, n
            jacob(i, 1) = df_da(x(i), a, b)
            jacob(i, 2) = df_db(x(i), a, b)
        end do
        suma1 = sum(r(x, y, a, b))/(1d0*n)
        jacob_t = transpose(jacob)
        d = matmul(jacob_t, jacob)
        !call invert_matrix(d)
        call inverse(d,d2,num)
        delta = -matmul(matmul(d2, jacob_t), r(x, y, a, b))
        a = a + alpha * delta(1)
        b = b + alpha * delta(2)
        suma2 = sum(r(x, y, a, b))/(1d0*n)
        count = count + 1

        d_old = d
        jacob_old = jacob
        if(suma1>suma2) then
          alpha = alpha/10.d0
          print*,alpha
        else
          go to 50
        end if

        50 continue
        if (delta(1) < epsilon .and. delta(2) < epsilon) then
            exit
        end if


    end do

  !  4 continue

    if (count == jmax) then
        write (*, *) "Calculation loop exceeds the limit. Result may not be correct. Try change alpha lower."
    end if

    print*,'----------------------------------------------------------------'
    print*,'Parameters'
    print*,'----------------------------------------------------------------'

    write (*, '(a, i0)') "Iteration count =", count
    write (*, '(2(a, f15.5))') "a=", a, ", b=", b

    do i = 1, n
        write(*,*) y(i) , f(x(i),a,b)
    end do

    write (11, '(a, i0)') "Iteration count =", count
    write (11, '(2(a, f15.5))') "a=", a, ", b=", b

    write(11, *) 'data' , '                        fit'

    do i = 1, n
        write(11,*) y(i) , f(x(i),a,b)
    end do

contains

    double precision function f(x, a, b)
        double precision,intent(in) :: x, a, b
        double precision :: c
        c = 96.05
        f = (1 - (a*x)/b)**(-1 + 1/(a*c))
        !f = (1-a*x/b)**((1/a*c)-1)
    end function

    double precision function df_da(x, a, b)
        double precision,intent(in) :: x, a, b
        double precision :: term1,term2,c
        c = 96.05
        term1 = b*(1 - (a*x)/b)**(1/(a*c))*(a*(-1 + a*c)*x + (-b + a*x)*log(1 - (a*x)/b))
        term2 = (a**2)*c*(b - a*x)**2
        !df_da = (((1-a*x/b)**((1/a*c)-1)))*(-(-1+c/a)*x/(b*(1-a*x/b))-c*log(1-a*x/b)/a**2)
        df_da = term1/term2
    end function

    double precision function df_db(x, a, b)
        double precision,intent(in) :: x, a, b
        double precision :: c
        c = 96.05
        df_db = (a*(-1 + 1/(a*c))*x*(1 - (a*x)/b)**(-2 + 1/(a*c)))/b**2
        !df_db = (a*(-1+c/a)*x)*((1-a*x/b)**(-2+c/a))/b**2
    end function


    function r(x, y, a, b)
        double precision,intent(in) :: x(:), y(:), a, b
        double precision r(size(x)),error
        integer i
        do i = 1, size(x)
            r(i) = y(i) - abs(f(x(i), a, b))
        end do

    end function

          subroutine inverse(a,c,n)
      !============================================================
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      ! Alex
      !-----------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - inverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed
      ! during the calculation
      !===========================================================
      implicit none
      integer n
      double precision a(n,n), c(n,n)
      double precision L(n,n), U(n,n), b(n), d(n), x(n)
      double precision coeff
      integer i, j, k

      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
      L=0.0
      U=0.0
      b=0.0

      ! step 1: forward elimination
      do k=1, n-1
         do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
               a(i,j) = a(i,j)-coeff*a(k,j)
            end do
         end do
      end do

      ! Step 2: prepare L and U matrices
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
        L(i,i) = 1.0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
        do i=1,j
          U(i,j) = a(i,j)
        end do
      end do

      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
        b(k)=1.0
        d(1) = b(1)
      ! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
          d(i)=b(i)
          do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
          end do
        end do
      ! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
          x(i) = d(i)
          do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
          end do
          x(i) = x(i)/u(i,i)
        end do
      ! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
          c(i,k) = x(i)
        end do
        b(k)=0.0
      end do
      end subroutine inverse


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
