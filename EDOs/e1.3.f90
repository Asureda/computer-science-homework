!Author: Alexandre Sureda Croguennoc
PROGRAM CHEM
	IMPLICIT NONE
	INTEGER n
	real*8 acc,ta,tb,h,a,b,g,d,t,v_o,alpha
	REAL*8, DIMENSION(4) :: y1,y2,y3
	!REAL*8, DIMENSION(2),EXTERNAL :: dC_dt
	COMMON/PARAMS/acc
	COMMON/EULER/a,b,g,d
	ta=0d0
	tb=100d0
	n=1000
	h=1d-3
	v_o=5d0 !m/s
	alpha=30d0*acos(-1.)/180d0 ! degrees to radians
	acc=-9.81
	y1=(/0d0,0d0,v_o*sin(alpha),v_o*cos(alpha)/)
	y2=y1;y3=y1
	OPEN(11,FILE='euler2.dat')
	t=ta
	WRITE(11,*) t,y1(1),y1(2),y2(1),y2(2),y3(1),y3(2)
	DO WHILE ((t.lt.tb).and.(y1(1).ge.0d0))
		!print*,'hola'
		t=t+h
		a=1d0;b=0d0;g=0d0;d=0d0
		CALL EULER_GENERAL(t,y1,h)!Euler simple
		a=1d0;b=0d0;g=0d0;d=0d0
		CALL EULER_GENERAL(t,y2,h)!Euler modificat
		a=0d0;b=1d0;g=5d-1;d=5d-1
		CALL EULER_GENERAL(t,y3,h)!Euler millorat
		WRITE(11,*) t,y1(1),y1(2),y2(1),y2(2),y3(1),y3(2)
	END DO
	CLOSE(11)
END PROGRAM
SUBROUTINE EULER_GENERAL(x,y,h)
	IMPLICIT NONE
	REAL*8 x,x_2,h,a,b,g,d
	REAL*8, DIMENSION(4) :: y, y_1,y_2
	COMMON/EULER/a,b,g,d
	y_1=y
	CALL dF_dT(y_1,x)
	y_2=y+y_1*d*h
	x_2=x+g*h
	CALL dF_dT(y_2,x_2)
	y=y+h*(a*y_1+b*y_2)
RETURN
END SUBROUTINE
SUBROUTINE dF_dt(y,x)
	IMPLICIT NONE
	REAL*8 x,v,acc
	REAL*8, DIMENSION(4) :: y,aux
	COMMON/PARAMS/acc
	aux=y
	y(1)=aux(3)
	y(2)=aux(4)
	y(3)=acc
	y(4)=0.0
RETURN
END SUBROUTINE
