PROGRAM FICK_1D
    IMPLICIT NONE
    INTEGER n_x,n_t,i,j,wri
    REAL*8 C_oA,C_oB,L,dx,dt,t,D,t_a,t_b,r,k,exp
    REAL*8, DIMENSION(:), ALLOCATABLE :: C_Anew,C_Aold,C_Bnew,C_Bold,a,b,c
    INTEGER, DIMENSION(:), ALLOCATABLE :: cpc
    !Defining the physical parameters of the problem
    L=10.
    D=1.
    k=0.01
    C_oA=7.5
    C_oB=7.0
    n_x=100
    dt=1d-5  !Defining time increments to have a good integration
    t_a=0d0
    t_b=0.001
    dx=L/(1d0*n_x)
    n_t=int((t_b-t_a)/(1d0*dt))
    wri=int(n_t/500)
    r=D*dt/(dx*dx)
    ALLOCATE(C_Anew(1:n_x),C_Aold(1:n_x),C_Bnew(1:n_x),C_Bold(1:n_x),cpc(0:n_x+1),a(1:n_x),b(1:n_x),c(1:n_x))
    OPEN(11,FILE='evolution_Fick_4.dat')
    OPEN(12,FILE='evolution_Fick_42.dat')
    !Definig the initial conditions and boundary values
    C_Anew=0d0
    C_Bnew=0d0
    DO i=1,n_x
        cpc(i)=i
        t=i*dx
        C_Anew(i)=C_oA
        C_Bnew(i)=C_oB
        IF((t.gt.4.95).and.(t.lt.5.05))THEN
            C_Anew(i)=C_oA+1.5
            C_Bnew(i)=C_oB-1.25
        END IF
    END DO
    cpc(0)=n_x
    cpc(n_x+1)=1
    DO j=1,n_x
        WRITE(11,*)j*dx,C_Anew(j),C_Bnew(j)
    END DO
    WRITE(11,*)
    WRITE(11,*)
    !Definig the tridiagonal matrix
    DO i=1,n_x
        B(i)=2d0*(1+r)
        A(i)=-r
    END DO
    C=A
    !A(1)=
    !Starting time loop
    DO i=1,n_t
        t=t_a+i*dt
        !Creating the right concentration array for the matrix equation
        DO j=1,n_x
            C_Aold(j)=r*C_Anew(cpc(j-1))+r*C_Anew(cpc(j+1))+2d0*(1-r-k*dt*C_Bnew(j))*C_Anew(j)
            C_Bold(j)=r*C_Bnew(cpc(j-1))+r*C_Bnew(cpc(j+1))+2d0*(1-r-k*dt*C_Anew(j))*C_Bnew(j)
        END DO
        !print*,C_Aold
        !Solving the matrix equation for the two concentrations
        CALL TRIDIAG(a,b,c,C_Aold,C_Anew,n_x)
        CALL TRIDIAG(a,b,c,C_Bold,C_Bnew,n_x)
        !print*,C_Anew
        ! Writting the concentration values every few iterations
        IF(mod(i,wri).eq.0) THEN
            DO j=1,n_x
                WRITE(11,*)j*dx,C_Anew(j),C_Bnew(j)
            END DO
            WRITE(11,*)
            WRITE(11,*)
            WRITE(12,*)t,sum(C_Anew),sum(C_Bnew)

        END IF
    END DO
END PROGRAM
SUBROUTINE TRIDIAG(A,B,C,R,PSI,IMAX)
! c Solves the problem T psi =R
! c
! c where T is a tridiagonal matrix, A (lower), B (central), C (upper)
! c  A   0 a1, a2, a3, ...., aIMAX
! c  B   b1 b2, b3, b4, ...., bIMAX
! c  C   c1, c2, c3, c4, ....,0

! c       IMPLICIT double precision (A-H,K,O-Z)
! c       IMPLICIT INTEGER (I-J , L-N)
        IMPLICIT NONE
        integer imax,j
        double precision  BET
        double precision GAM(4001)
        double precision A(IMAX),B(IMAX),C(IMAX),R(IMAX),PSI(IMAX)

        IF(B(1).EQ.0.) STOP
        BET=B(1)
        PSI(1)=R(1)/BET
        DO 11 J=2,IMAX
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
        IF(BET.EQ.0) EXIT
        PSI(J)=(R(J)-A(J)*PSI(J-1))/BET
11      CONTINUE

        DO 12 J=IMAX-1,1,-1
        PSI(J)=PSI(J)-GAM(J+1)*PSI(J+1)
12      CONTINUE

       RETURN
END
