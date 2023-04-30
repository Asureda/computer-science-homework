!Author: Alexandre Sureda Croguennoc
PROGRAM CHEM
	IMPLICIT NONE
	INTEGER n
	real*8 al,be,ga,de,ta,tb,h,a,b,g,d,t
	REAL*8, DIMENSION(2) :: y1,y2,y3
	!REAL*8, DIMENSION(2),EXTERNAL :: dC_dt
	COMMON/PARAMS/al,be,ga,de
	COMMON/EULER/a,b,g,d
	ta=0d0
	tb=20d0
	n=1000
	h=1d-3
	al=2.;be=1.1;ga=1.;de=0.9;
	y1=(/1.,0.5/)
	y2=y1;y3=y1
	OPEN(11,FILE='euler3.dat')
	t=ta
	DO WHILE (t.lt.tb)
		t=t+h
		a=1d0;b=0d0;g=0d0;d=0d0
		CALL EULER_GENERAL(t,y1,h)!Euler simple
		a=1d0;b=0d0;g=0d0;d=0d0
		CALL EULER_GENERAL(t,y2,h)!Euler modificat
		a=0d0;b=1d0;g=5d-1;d=5d-1
		CALL EULER_GENERAL(t,y3,h)!Euler millorat
		WRITE(11,*) t,y1,y2,y2,-al*ga+de*be*y1(1)*y1(2)
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
	REAL*8 x,al,be,ga,de
	REAL*8, DIMENSION(2) :: y,aux
	COMMON/PARAMS/al,be,ga,de
	aux=y
	y(1)=al*aux(1)-be*aux(1)*aux(2)
	y(2)=-ga*aux(2)+de*aux(1)*aux(2)
RETURN
END SUBROUTINE
