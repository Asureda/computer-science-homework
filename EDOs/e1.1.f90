!Author:Alexandre Sureda Croguennoc
PROGRAM EULER
	IMPLICIT NONE
	INTEGER n
	real*8 t0,T_F,k,ta,tb,h,Temp1,Temp2,Temp3,EULER_GENERAL,t
	REAL*8, EXTERNAL :: dT_dt
	COMMON/PARAMS/T_F,k
	ta=0d0
	tb=100d0
	n=1000
	h=1d0
	Temp1=298d0;Temp2=298d0;Temp3=298d0
	T_F=273d0
	k=0.05
	t=ta
	OPEN(11,FILE='euler.dat')
	DO WHILE(t.le.tb)
		t=t+h
		print*,t
		Temp1=EULER_GENERAL(t,Temp1,h,dT_dt,1d0,0d0,0d0,0d0)!Euler simple
		Temp2=EULER_GENERAL(t,Temp2,h,dT_dt,0d0,1d0,5d-1,5d-1)!Euler modificat
		Temp3=EULER_GENERAL(t,Temp3,h,dT_dt,5d-1,5d-1,1d0,1d0)!Euler millorat
		WRITE(11,*) t,Temp1,Temp2,Temp3
	END DO
	CLOSE(11)
END PROGRAM
FUNCTION EULER_GENERAL(x,y,h,dF_dX,a,b,g,d)
	IMPLICIT NONE
	REAL*8 :: EULER_GENERAL,x,y,h,dF_dX,a,b,g,d
	EULER_GENERAL=y+h*(a*dF_dX(y,x)+b*dF_dX(y+dF_dX(y,x)*d*h,x+g*h))
RETURN
END FUNCTION
FUNCTION dT_dt(y,x)
	IMPLICIT NONE
	REAL*8 :: dT_dt, y,x,T_F,k
	COMMON/PARAMS/T_F,k
	dT_dt=k*(T_F-y)
RETURN
END
