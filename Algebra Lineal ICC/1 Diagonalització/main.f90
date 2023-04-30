PROGRAM MAIN
    !----------------------------------------------------------------
    !  PROGRAM TO CHECK THE ROUTINES FOR SYMETRIC MATRICES
    !  WITH THE HILBERT MATRIX
    !
    !  CHANGE N FOR THE DIMENSIONS OF THE MATRICES
    !
    use diag
    IMPLICIT NONE
    INTEGER :: order,i,j
    REAL*8 eps,quality
    REAL*8, DIMENSION(:), ALLOCATABLE :: M_epsil(:)
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: M,M_eiv(:,:),M_trid(:,:)
    order=5
    ALLOCATE(M(order,order),M_epsil(order),M_eiv(order,order),M_trid(order,order))
    !Creating Hilbert MATRIX
    DO i=1,order
        DO j=1,order
            IF(i.eq.j) THEN
                M(i,j)=-1d0/(2d0*i-1d0)
            ELSE
                M(i,j)=-1d0/(10.d0*(i*1d0+j*1d0-1d0))
            END IF
        END DO
    END DO

    WRITE(11,*)
    WRITE(11,*)'-----------------------------------------------------'
    WRITE(11,*)'                HILBERT MATRIX                        '
    WRITE(11,*)'-----------------------------------------------------'
    DO i=1,order
        WRITE(11,*)M(i,:)
    END DO
    WRITE(11,*)
    WRITE(11,*)
    WRITE(11,*)'-----------------------------------------------------'
    WRITE(11,*)'                JACOBI METHOD                        '
    WRITE(11,*)'-----------------------------------------------------'
    WRITE(11,*)'EIGENVALUES MATRIX'
    eps=1d-8
    call jacobi(M,order,M_epsil,M_eiv,eps,quality)
    DO i=1,order
      WRITE(11,*) M_eiv(i,:)*M_epsil(i)
    END DO

    write(11,*)'quality ',quality
    !----------------------------------------------------------------

    WRITE(11,*)
    WRITE(11,*)
    WRITE(11,*)'-----------------------------------------------------'
    WRITE(11,*)'                HOUSEHOLDER METHOD                        '
    WRITE(11,*)'-----------------------------------------------------'
    WRITE(11,*)'TRIDIAGONAL MATRIX'
    CALL HOUSEHOLDER(M,M_trid,order)
    DO i=1,order
        WRITE(11,*)M_trid(i,:)
    END DO
    WRITE(11,*)'EIGENVALUES MATRIX (with Jacobi)'

    call jacobi(M,order,M_epsil,M_eiv,eps,quality)
    DO i=1,order
      WRITE(11,*) M_eiv(i,:)*M_epsil(i)
    END DO

    WRITE(11,*)
    WRITE(11,*)
    WRITE(11,*)'-----------------------------------------------------'
    WRITE(11,*)'                LANCZOS-DAVIDSON                        '
    WRITE(11,*)'-----------------------------------------------------'
    WRITE(11,*)'TRIDIAGONAL MATRIX'
    CALL LANCZOS(M,M_trid,order)
    DO i=1,order
        WRITE(11,*)M_trid(i,:)
    END DO
    WRITE(11,*)'EIGENVALUES MATRIX (with Jacobi)'
    call jacobi(M,order,M_epsil,M_eiv,eps,quality)
    DO i=1,order
      WRITE(11,*) M_eiv(i,:)*M_epsil(i)
    END DO
    close(11)



END PROGRAM MAIN
