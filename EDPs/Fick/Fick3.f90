PROGRAM FICK_1D
    IMPLICIT NONE
    INTEGER n_x,n_t,i,j,wri
    REAL*8 C_o,L,dx,dt,t,D,t_a,t_b
    REAL*8, DIMENSION(:), ALLOCATABLE :: C_new,C_old
    INTEGER, DIMENSION(:), ALLOCATABLE :: cpc
    !Defining the physical parameters of the problem
    L=10.
    D=1.
    C_o=7.5
    n_x=100
    dt=1d-5  !Defining time increments to have a good integration
    t_a=0d0
    t_b=10d0
    dx=L/(1d0*n_x)
    n_t=int((t_b-t_a)/(1d0*dt))
    wri=int(n_t/500)
    ALLOCATE(C_new(0:n_x+1),C_old(0:n_x+1),cpc(0:n_x+1))
    OPEN(11,FILE='evolution_Fick_3.dat')
    !Definig the initial conditions and boundary values
    C_old=0d0
    DO i=1,n_x
        !Array for periodic boundary conditions
        cpc(i)=i
        t=i*dx
        !Inintial concentration as a central
        IF((t.gt.4.8).and.(t.lt.5.2))THEN
            C_old(i)=C_o
        END IF
    END DO
    cpc(0)=n_x
    cpc(n_x+1)=1
    C_new=C_old
    DO j=0,n_x+1
        WRITE(11,*)j*dx,C_old(j)
    END DO
    WRITE(11,*)
    WRITE(11,*)
    !Starting time loop
    DO i=1,n_t
        t=t_a+i*dt
        !Updating the concentrations for the new time
        !print*,i,'---------------------------------------------'
        DO j=1,n_x
            C_new(j)=C_old(j)+D*dt*(C_old(cpc(j+1))+C_old(cpc(j-1))-2d0*C_old(cpc(j)))/(dx*dx)
            !print*,C_old(j-1),C_old(j),C_old(j+1),C_new(j),j
        END DO
        C_old=C_new
        ! Writting the concentration values every few iterations
        IF(mod(i,wri).eq.0) THEN
            DO j=1,n_x
                WRITE(11,*)j*dx,C_old(j)
            END DO
            WRITE(11,*)
            WRITE(11,*)
        END IF
    END DO
END PROGRAM
