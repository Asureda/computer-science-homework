!Author: Alexandre Sureda Croguennoc
PROGRAM STEEPEST_DESCENT
	IMPLICIT NONE
	INTEGER i,j,nx,ny,file
	REAL*8 t,h,ENERGY,ax,bx,ay,by,hx,hy,k
	REAL*8, DIMENSION(2) :: y,vo
	REAL*8, DIMENSION(2,2) :: hes
	REAL*8, EXTERNAL :: GRADIENT,PATH_HESIAN

	!------------------------------------
	!     PLOTTING THE ENERGY SURFACE
	!------------------------------------
	!TO PLOT THE SURFACE PUT FILE=1
	file=0
	if(file.eq.1)then
	ax=-1.5; bx=0.5
	ay=-0.5; by=1.5

	nx=500
	ny=500
	hx=(bx-ax)/(1d0*nx)
	hy=(by-ay)/(1d0*ny)
	print*,'enter'

	OPEN(11,FILE='surface.dat')
	DO i=0,nx
		DO j=0,ny
			WRITE(11,*)ax+i*hx,ay+j*hy,ENERGY(ax+i*hx,ay+j*hy)
		END DO
		WRITE(11,*)
	END DO
	print*,'out'
	close(11)
	end  if
	!---------------------------------------
	!INITAL DIRECTION
	!vo=(/-0.493,0.493/) !hessien negative eigenvalue
	vo=(/-0.493,-0.493/)
	k=0.1
	!----------------------------------------
	!   STEEPEST DESCENT PATH 1
	!---------------------------------------
	OPEN(12,FILE='steepest_t.dat')
	h=0.00001
	t=0.0
	y=(/-0.822,0.624/)
	WRITE(12,*)y,ENERGY(y(1),y(2))
	CALL HESSIAN(y,hes)
	print*,'xx',hes(1,1),'xy',hes(1,2)
	print*,'yx',hes(2,1),'yy',hes(2,2)
	!COMPUTING THE EIGENVCTOR FOR THE HESSIAN
	!OUT OF THE PROGRAM
	y=y+k*vo
	WRITE(12,*)y,ENERGY(y(1),y(2))
	do i=1,500
		t=t+h
		CALL RK_4(y,t,h,GRADIENT)
		WRITE(12,*)y,ENERGY(y(1),y(2))
	end do
	WRITE(12,*)
	WRITE(12,*)
	!----------------------------------------
	!   STEEPEST DESCENT PATH 2
	!---------------------------------------
	y=(/-0.822,0.624/)
	WRITE(12,*)y,ENERGY(y(1),y(2))
	CALL HESSIAN(y,hes)
	y=y-k*vo
	WRITE(12,*)y,ENERGY(y(1),y(2))
	do i=1,500
		t=t+h
		CALL RK_4(y,t,h,GRADIENT)
		WRITE(12,*)y,ENERGY(y(1),y(2))
	end do
	close(12)
	!----------------------------------------
	!   REACTION PATH 1 (HESSIAN)
	!---------------------------------------
	OPEN(13,FILE='hessianpath_t.dat')
	h=0.00000001
	t=0.0
	y=(/-0.822,0.624/)
	WRITE(13,*)y,ENERGY(y(1),y(2))
	y=y+k*vo
	WRITE(13,*)y,ENERGY(y(1),y(2))
	do i=1,100
		t=t+h
		CALL RK_4(y,t,h,PATH_HESIAN)
		WRITE(13,*)y,ENERGY(y(1),y(2))
		!print*,y
	end do
	WRITE(13,*)
	WRITE(13,*)
	!----------------------------------------
	!   REACTION PATH 2 (HESSIAN)
	!----------------------------------------
	t=0.0
	y=(/-0.822,0.624/)
	WRITE(13,*)y,ENERGY(y(1),y(2))
	y=y-k*vo
	WRITE(13,*)y,ENERGY(y(1),y(2))
	do i=1,100
		t=t+h
		CALL RK_4(y,t,h,PATH_HESIAN)
		WRITE(13,*)y,ENERGY(y(1),y(2))
	end do
	close(13)
END PROGRAM

SUBROUTINE RK_4(y,t,h,dF)
	REAL*8 h,t,eps, t1,t2,t3,t4
	REAL*8,DIMENSION(2) :: y,y1,y2,y3,y4,a1,a2,a3,a4
	y1=y
	t1=t
	CALL dF(y1,t1)!f_o
	y2=y+h*y1*5d-1
	t2=t+h*5d-1
	CALL dF(y2,t2)!f_1
	y3=y+h*y2*5d-1
	t3=t+h*5d-1
	CALL dF(y3,t3)!f_2
	y4=y+y3*h
	t4=t*h
	CALL dF(y4,t4)!f_2
	y=y+h*(+y1+2d0*y2+2d0*y3+Y4)*6d0
	RETURN
END SUBROUTINE
SUBROUTINE GRADIENT(y,t)
	IMPLICIT NONE
	REAL*8 ENERGY,t,h
	REAL*8, DIMENSION(2) :: y,aux
	h=1d-3
	aux=y
	y(1)=-(ENERGY(aux(1)+h,aux(2))-ENERGY(aux(1)-h,aux(2)))/(2d0*h)
	y(2)=-(ENERGY(aux(1),aux(2)+h)-ENERGY(aux(1),aux(2)-h))/(2d0*h)
	RETURN
END SUBROUTINE
SUBROUTINE PATH_HESIAN(y,t)
	IMPLICIT NONE
	REAL*8 ENERGY,t,h
	REAL*8, DIMENSION(2) :: y,aux,grad
	REAL*8, DIMENSION(2,2) :: hess,adjoint
	aux=y
	CALL HESSIAN(aux,hess)
	adjoint(1,1)=hess(2,2);adjoint(1,2)=-hess(2,1)
	adjoint(2,1)=-hess(1,2);adjoint(2,2)=hess(1,1)
	grad=y
	CALL GRADIENT(grad,t)
	aux=y
	y(1)=adjoint(1,1)*grad(1)+adjoint(1,2)*grad(2)
	y(2)=adjoint(2,1)*grad(1)+adjoint(2,2)*grad(2)
	RETURN
END


FUNCTION ENERGY(x,y)
	IMPLICIT NONE
	INTEGER i
	REAL*8 X,Y,ENERGY
	REAL*8,DIMENSION(4)::Ai,a,b,c,xo,yo
	Ai=(/-200.,-100.,-170.,15./)
	a=(/-1.,-1.,-6.5,0.7/)
	b=(/0.,0.,11.,0.6/)
	c=(/-10.,-10.,-6.5,0.7/)
	xo=(/1.,0.,-0.5,-1./)
	yo=(/0.,0.5,1.5,1./)
	ENERGY=0d0
	DO I=1,4
		ENERGY=ENERGY+Ai(i)*dexp(a(i)*(x-xo(i))**2d0+b(i)*(x-xo(i))*(y-yo(i))+c(i)*(y-yo(i))**2d0)
	END DO
	RETURN
END FUNCTION
SUBROUTINE HESSIAN(y,hes)
	REAL*8 y(2),hes(2,2),ENERGY,h
	h=1d-3
	hes(1,1)=(ENERGY(y(1)+2d0*h,y(2))+ENERGY(y(1)-2d0*h,y(2))-2d0*ENERGY(y(1),y(2)))/(4d0*h**2d0)
	hes(1,2)=(ENERGY(y(1)+h,y(2)+h)-ENERGY(y(1)-h,y(2)+h)-ENERGY(y(1)+h,y(2)-h)+ENERGY(y(1)-h,y(2)-h))/(4d0*h**2d0)
	hes(2,1)=(ENERGY(y(1)+h,y(2)+h)-ENERGY(y(1)+h,y(2)-h)-ENERGY(y(1)-h,y(2)+h)+ENERGY(y(1)-h,y(2)-h))/(4d0*h**2d0)
	hes(2,2)=(ENERGY(y(1),y(2)+2d0*h)+ENERGY(y(1),y(2)-2d0*h)-2d0*ENERGY(y(1),y(2)))/(4d0*h**2d0)
	RETURN
END SUBROUTINE
