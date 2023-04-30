!Author: Alexandre Sureda Croguennoc
PROGRAM CHEM
	IMPLICIT NONE
	INTEGER n
	real*8 k1,k2,ta,tb,h,a,b,g,d,t
	REAL*8, DIMENSION(2) :: C1,C2,C3
	COMMON/PARAMS/k1,k2
	COMMON/EULER/a,b,g,d
	ta=0d0
	tb=500d0
	n=1000
	h=2d0
	k1=1d-4
	k2=8
	C1=(/1d-3,0d0/)
	C2=C1;C3=C1
	OPEN(11,FILE='euler1.dat')
	t=ta
	DO WHILE ((t.le.tb).and.(C1(1).gt.0d0))
		t=t+h
		a=1d0;b=0d0;g=0d0;d=0d0
		CALL EULER_GENERAL(t,C1,h)!Euler simple
		a=1d0;b=0d0;g=0d0;d=0d0
		CALL EULER_GENERAL(t,C2,h)!Euler modificat
		a=0d0;b=1d0;g=5d-1;d=5d-1
		CALL EULER_GENERAL(t,C3,h)!Euler millorat
		WRITE(11,*) t,C1,C2,C3
	END DO
	CLOSE(11)
END PROGRAM
SUBROUTINE EULER_GENERAL(x,y,h)
	IMPLICIT NONE
	REAL*8 x,x_2,h,a,b,g,d
	REAL*8, DIMENSION(2) :: y, y_1,y_2
	COMMON/EULER/a,b,g,d
	y_1=y
	CALL dC_dT(y_1,x)
	y_2=y+y_1*d*h
	x_2=x+g*h
	CALL dC_dT(y_2,x_2)
	y=y+h*(a*y_1+b*y_2)
RETURN
END SUBROUTINE
SUBROUTINE dC_dt(y,x)
	IMPLICIT NONE
	REAL*8 x,k1,k2,v
	REAL*8, DIMENSION(2) :: y
	COMMON/PARAMS/k1,k2
	v=k1*y(1)+k2*y(2)*y(1)
	y(1)=-2d0*v
	y(2)=+2d0*v
RETURN
END SUBROUTINE
