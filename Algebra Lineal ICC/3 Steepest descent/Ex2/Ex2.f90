program steepest_descent
  IMPLICIT NONE
  integer, parameter :: pr = selected_real_kind(15,3)
  integer :: i,num_var
  real(kind=16) :: f,f_prime,term1,term2,alpha,df1,df2,A,aa,rate1,rate2,rate,rate_sum,norm
  real(pr) :: fp, h,quality,rel_tol,er,ea
  real(pr) :: x1_0,x2_0,z0,tmp_z0
  real(pr) :: x1_new,x2_new
  integer :: nb_iter,nb_max_iter
  real(pr) :: eps,cond
  double precision, DIMENSION(:), ALLOCATABLE :: jacob(:),d(:),d_old(:),b(:),Q_epsil(:)
  double precision, DIMENSION(:,:), ALLOCATABLE :: Q(:,:),Q_eiv(:,:)
  real(8) :: dot
  logical :: normalized


  alpha = 0.1 ! learning rate
  nb_max_iter = 100 ! Nb max d'iteration
  eps = 0.00001 ! stop condition
  num_var = 2
  er=1d-6
  ea=1d-6
  normalized = (.FALSE.)
  open(12,file='results.dat')
  allocate(jacob(num_var))
  allocate(d(num_var))
  allocate(d_old(num_var))
  allocate(b(num_var))
  allocate(Q(num_var,num_var))
  allocate(Q_eiv(num_var,num_var))
  allocate(Q_epsil(num_var))



  x1_0 = 1/10.0 ! start point
  x2_0 = 1.0
  z0 = f(x1_0,x2_0)

  cond = eps + 10.0 ! start with cond greater than eps (assumption)
  nb_iter = 0
  tmp_z0 = z0
  do i = 1,num_var
    jacob(i) = f_prime(i,x1_0,x2_0)
  end do
  do i = 1,num_var
    d(i) = -f_prime(i,x1_0,x2_0)
  end do
rate = 0.d0
do while(cond>eps.and.nb_iter<nb_max_iter)
  d_old = d
  if(normalized.eqv.(.TRUE.)) then
    do i = 1,num_var
      d(i) = f_prime(i,x1_0,x2_0)
      norm = norm + d(i)*d(i)
    end do
    d(1:num_var) = -d(1:num_var)/sqrt(norm)
    norm = norm/num_var
    else
      do i = 1,num_var
        d(i) = -f_prime(i,x1_0,x2_0)
      end do
    end if


  term1 = d(1)**2 + d(2)**2
  term2 = 2*(10*d(1)**2+d(2)**2)
  alpha = term1/term2

  x1_new = x1_0 + alpha*d(1)
  x2_new = x2_0 + alpha*d(2)
  dot = dot_product(d_old,d )
  write(12,*)'k',nb_iter ,'  s(k+1) and s(k) orthogonal : ',dot
  rate1 = abs(x1_new/x1_0)
  rate2 = abs(x2_new/x2_0)
  rate = rate+(rate1+rate2)/2.d0
  x1_0 = x1_new
  x2_0 = x2_new

  z0 = f(x1_0,x2_0)
  nb_iter = nb_iter + 1
  cond = abs( tmp_z0 - z0 )
  tmp_z0 = z0

  end do
  rate = rate/nb_iter

  Q(1,1) = 20
  Q(1,2) = 0
  Q(2,1) = 0
  Q(2,2) = 2



  call jacobi(Q,num_var,Q_epsil,Q_eiv,eps,quality)

  A = Q_epsil(1)
  if(A > Q_epsil(2)) then
    A = Q_epsil(1)
    aa = Q_epsil(2)
  end if

  write(12,*)'numerical convergence rate ',rate
  write(12,*)'convergence rate calculated from eigenvalues',(A-aa)/(A+aa)

  end program steepest_descent

  real(kind=16) function f(x1,x2)
    implicit none
    integer, parameter :: pr = selected_real_kind(15,3)
    real(pr), intent(in) :: x1,x2
    !f = 5*x1**2+x2**2+4*x1*x2-14*x1-6*x2+20
    f = 10*x1**2+x2**2
  end function f

  real(kind=16) function df1(x1,x2)
    implicit none
    integer, parameter :: pr = selected_real_kind(15,3)
    real(pr), intent(in) :: x1,x2
    df1 = -10*x1-4*x2+14
  end function df1

  real(kind=16) function df2(x1,x2)
    implicit none
    integer, parameter :: pr = selected_real_kind(15,3)
    real(pr), intent(in) :: x1,x2
    df2 = -2*x2-4*x1+6
  end function df2



  real(kind=16) function f_prime(dpi,x1,x2)
    implicit none
    integer, parameter :: pr = selected_real_kind(15,3)
    integer, intent(in) :: dpi
    real(pr), intent(in) :: x1,x2
    real(pr) :: x1_p,x2_p
    real(kind=16) :: f
    real(pr) :: h
    h = 0.00001
    x1_p = x1
    x2_p = x2
    if( dpi == 1 ) x1_p = x1_p + h
    if( dpi == 2 ) x2_p = x2_p + h
    f_prime = ( f(x1_p,x2_p) - f(x1,x2) ) / h
  end function f_prime



  SUBROUTINE jacobi(H,order,epsil,U,eps,quality)

    integer     :: i,j,k,l,iter,imax,jmax,order
    real(8)     :: eps,temp1,temp2,Smax,quality,theta,pi
    parameter   (pi = dacos(-1.d0))
    real(8),    dimension(order)          ::epsil
    real(8),    dimension(order,order)    :: H,U,H2
    real(8),    dimension(order)          ::col1,col2

    ! Step1 : Copy H to H2
    DO i = 1,order
      DO j = 1,order
        H2(i,j) = H(i,j)
      END DO
    END DO
    ! Step 2 : Initialization of U
    DO i = 1,order
      DO j = 1,order
        U(i,j) = 0.d0
      END DO
      U(i,i) = 1.d0
    END DO

    !open(14,file='jacobi.out')

    ! Step 3 : Find the largest off-diagonal element from H2

    DO iter = 1,500000
      ! Search the off-diagonal element with max module
      Smax = abs(H2(1,2))
      imax = i
      jmax = j

      DO i =1,order-1
        DO j = i+1,order
          if(abs(H(i,j)).ge.Smax) then
            Smax = H(i,j)
            imax=i
            jmax=j
          end if
        END DO
      END DO
      write(14,*) 'diagonalization in progress Smax = ', Smax
      if (Smax.le.eps) then
        go to 4
      end if

      ! Step 4 : Search the rotational angle
      theta = 0.d0
      if ( H2(imax,imax).eq.H2(jmax,jmax)) then
        if (H2(imax,jmax).ge.0.d0) then
          theta = pi/4.d0
        else
          theta = -pi/4.d0
        end if
      else
        theta = 0.5*atan(2.d0*H2(imax,jmax)/(H2(imax,imax)-H2(jmax,jmax)))
      end if

      ! Step 5 : Compute matrix U
      DO i = 1,order
        col1(i) = U(i,imax)
        col2(i) = U(i,jmax)
      END DO

      DO i = 1,order
        U(i,imax) = cos(theta)*col1(i)+sin(theta)*col2(i)
        U(i,jmax) = cos(theta)*col2(i)-sin(theta)*col1(i)
      END DO

      DO i = 1,order
        col1(i) = H2(i,imax)
        col2(i) = H2(i,jmax)
      END DO
      DO i =1,order
        H2(i,imax) = cos(theta)*col1(i)+sin(theta)*col2(i)
        H2(i,jmax) = cos(theta)*col2(i)-sin(theta)*col1(i)
      END DO
      DO i = 1,order
        H2(imax,i) = H2(i,imax)
        H2(jmax,i) = H2(i,jmax)
      END DO
      H2(imax,imax)=cos(theta)**2*col1(imax)+sin(theta)**2*col2(jmax)  &
                  +2.d0*sin(theta)*cos(theta)*col1(jmax)
      H2(jmax,jmax)=cos(theta)**2*col2(jmax)+sin(theta)**2*col1(imax)  &
                  -2.d0*sin(theta)*cos(theta)*col1(jmax)
      H2(imax,jmax)=0.d0
      H2(jmax,imax)=0.d0

    END DO

  4 continue

  do i=1,order
       epsil(i)=H2(i,i)
  enddo


      quality=0.d0
      do i=1,order
       do j=1,order
        temp1=0.d0
        do k=1,order
         temp1=temp1+H(j,k)*U(k,i)
        enddo
        temp1=temp1-epsil(i)*U(j,i)
        quality=quality+temp1**2
       enddo
      enddo
      quality=dble(order)/quality

      return

    end subroutine jacobi

    SUBROUTINE HOUSEHOLDER(A_IN,A_OUT,n)

        IMPLICIT NONE
        INTEGER :: n,i,j,k,z
        REAL*8 :: alpha,r
        REAL*8, DIMENSION(n,n) :: A_IN, A_OUT,A_COPY,P,AUX
        REAL*8, DIMENSION(n) :: v
        !Copying the input matrix
        A_COPY=A_IN
        IF(n.ge.3)THEN
        !--------------------------------------------------------
        ! STARTING HOUSEHOLDER ROTATIONS
        !--------------------------------------------------------
        DO i=1,n-2
            !Computing alpha
            alpha=0d0
            DO j=i+1,n
                alpha=alpha+A_COPY(j,i)*A_COPY(j,i)
            END DO
            alpha=-dsign(1d0,A_COPY(i+1,i))*dsqrt(alpha)
            !Computing r
            r=dsqrt(0.5*(alpha*alpha-A_COPY(i+1,i)*alpha))
            DO k=1,i
                v(k)=0d0
            END DO
            !Computing v vector
            v(i+1)=(A_COPY(i+1,i)-alpha)/(2d0*r)
            DO k=i+2,n
                v(k)=A_COPY(k,i)/(2d0*r)
            END DO
            !Computing rotation matrix
            DO j=1,n
                DO k=1,n
                    IF (j.eq.k) THEN
                        P(j,k)=1-2d0*v(j)*v(k)
                    ELSE
                        P(j,k)=-2d0*v(j)*v(k)
                    END IF
                END DO
            END DO
            !Multiplying 3 maatrices to end the k-essim rotation
            !Multiply last two matrices A P
            AUX=0d0
            DO j=1,n
                DO k=1,n
                    DO z=1,n
                        AUX(j,k)=AUX(j,k)+A_COPY(j,z)*P(z,k)
                    END DO
                END DO
            END DO
            !Multiply P (A P) and giveing the next iteration matrix
            A_COPY=0d0
            DO j=1,n
                DO k=1,n
                    DO z=1,n
                        A_COPY(j,k)=A_COPY(j,k)+P(j,z)*AUX(z,k)
                    END DO
                END DO
            END DO
        END DO
        !End of Householder reduction
        !Changing the final tridiagonal to a clean from
        A_OUT=0d0
        A_OUT(1,1)=A_COPY(1,1)
        A_OUT(1,2)=A_COPY(1,2)
        A_OUT(n,n)=A_COPY(n,n)
        A_OUT(n,n-1)=A_COPY(n,n-1)
        DO j=2,n-1
            A_OUT(j,j-1)=A_COPY(j,j-1)
            A_OUT(j,j)=A_COPY(j,j)
            A_OUT(j,j+1)=A_COPY(j,j+1)
        END DO
        END IF
        RETURN
    END SUBROUTINE HOUSEHOLDER
