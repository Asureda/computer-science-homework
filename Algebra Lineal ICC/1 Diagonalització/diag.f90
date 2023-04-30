MODULE diag
  IMPLICIT NONE
contains
!Diagonalization modules
!Created by: Alexandre Sureda Croguennoc ICC 2020
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

    SUBROUTINE LANCZOS(A_IN,A_OUT,n)
        IMPLICIT NONE
        INTEGER n,i,j,k,h,n_basis
        REAL*8 :: alpha,beta,proj_vu,proj_uu,norm
        REAL*8, DIMENSION(n) :: q,z,q_aux1,q_aux2,q_old,q_old1,q_check
        REAL*8, DIMENSION(n,n) :: A_IN,A_OUT,gram_schmid
        A_OUT=0D0
        !Guessing first q array as the cannonical vector
        q=0D0
        q(n)=1d0
        gram_schmid(1,:)=q
        n_basis=1
        DO h=n,1,-1
            !Computing intermediate z to do one less matrix product
            z=0d0
            DO i=1,n
                DO j=1,n
                    z(i)=z(i)+q(j)*A_IN(j,i)
                END DO
            END DO
            q_old1=z
            !Computing the diagonal element for the tridiagonal matrix
            beta=0d0
            DO i=1,n
                beta=beta+z(i)*q(i)
            END DO
            A_OUT(h,h)=beta
            IF(h.ge.2) THEN
                !Computing the real z array
                DO i=1,n
                    z(i)=z(i)-beta*q(i)
                END DO
                !Computing the non-diagonal element alpha
                alpha=0d0
                DO i=1,n
                    alpha=alpha+z(i)*z(i)
                END DO
                A_OUT(h-1,h)=sqrt(alpha)
                A_OUT(h,h-1)=A_OUT(h-1,h)
                !Computing q for the next iteration
                q=z/A_OUT(h-1,h)
                q_aux1=q
                q_aux2=q_aux1
                !starting gram schmid decomposition to orthonormalize the new q
                DO k=1,n_basis-1
                    proj_vu=0D0
                    proj_uu=0d0
                    DO j=1,n
                        proj_uu=proj_uu+gram_schmid(k,j)*gram_schmid(k,j)
                        proj_vu=proj_vu+gram_schmid(k,j)*q_aux1(j)
                    END DO
                    q_aux2=q_aux2-proj_vu*gram_schmid(k,:)/proj_uu
                END DO
                q_aux2=q_aux1
                !normalizing the vector
                norm=0d0
                DO j=1,n
                    norm=norm+q_aux2(j)*q_aux2(j)
                END DO
                !Defining the normalized array q for next iteration
                q=q_aux2/sqrt(norm)
                n_basis=n_basis+1
                !Saving the q vector for the next gram-schmid
                gram_schmid(n_basis,:)=q
            END IF
        END DO
        RETURN
    END SUBROUTINE LANCZOS


  END MODULE diag
