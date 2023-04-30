PROGRAM FICK_1D
    IMPLICIT NONE
    INTEGER n_x,n_t,i,j,wri
    REAL*8 C_o,L,dx,dt,t,D,t_a,t_b
    REAL*8, DIMENSION(:), ALLOCATABLE :: C_new,C_old
    !Defining the physical parameters of the problem
    L=10.
    D=1.
    C_o=7.5
    n_x=100
    dt=1d-5  !Defining time increments to have a good integration
    t_a=0d0
    t_b=15d0
    dx=L/(1d0*n_x)
    n_t=int((t_b-t_a)/(1d0*dt))
    wri=int(n_t/500)
    ALLOCATE(C_new(0:n_x+1),C_old(0:n_x+1))
    OPEN(11,FILE='evolution_Fick_1.dat')
    !Definig the initial conditions and boundary values
    C_old=0d0
    print*,C_o,dx,dt
    C_old(0)=C_o
    C_old(n_x+1)=0
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

        DO j=1,n_x-1
            C_new(j)=C_old(j)+D*dt*(C_old(j+1)+C_old(j-1)-2d0*C_old(j))/(dx*dx)
            !print*,C_old(j-1),C_old(j),C_old(j+1),C_new(j),j
        END DO
        !Wall boundary condition for x=L
        C_new(n_x)=C_old(n_x)+D*dt*(C_old(n_x-1)-C_old(n_x))/(dx*dx)
        !Concentrations for the next iterations
        C_old=C_new
        !Checking for the boundary conditions
        C_old(0)=C_o
        C_old(n_x+1)=0
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
