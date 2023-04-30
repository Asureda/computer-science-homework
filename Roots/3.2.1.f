      PROGRAM ROOT_METHODS
      IMPLICIT NONE
      INTEGER n
      REAL*8 root
      open(11,FILE='method_results.dat')
      open(12,FILE='method_convergence.dat')

      CALL SEQUENTIAL(root,7d-1,1d-8,n)
      WRITE(*,*)n
      WRITE(11,200)'ROOT WITH SEQUENTIAL METHOD: X=',root,
     .'steps= ',real(n)
      WRITE(11,*)
      WRITE(12,*)
      WRITE(12,*)
      CALL BISECTION(root,5d-1,1d0,1d-8,n)
      WRITE(11,100)'ROOT WITH BISECTION METHOD: X=',root, 'steps= ',n
      WRITE(11,*)
      WRITE(12,*)
      WRITE(12,*)
      CALL SECANT(root,25d-2,5d-1,1d-8,n)
      WRITE(11,100)'ROOT WITH SECANT METHOD: X=',root, 'steps= ',n
      WRITE(11,*)
      WRITE(12,*)
      WRITE(12,*)
      CALL REGULA_FALSI(root,25d-2,1d0,1d-8,n)
      WRITE(11,100)'ROOT WITH REGUA FALSI METHOD: X=',root, 'steps= ',n
      WRITE(11,*)
      WRITE(12,*)
      WRITE(12,*)
      CALL 	NEWTON_RAPHSON(root,5d-1,1d-8,n)
      WRITE(11,100)'ROOT WITH NEWTON_RAPHSON METHOD: X=',root, 'steps= ',n
      WRITE(11,*)
      WRITE(12,*)
      WRITE(12,*)
      CALL 	MULLER(root,5d-1,1d0,7d-1,1d-8,n)
      WRITE(11,100)'ROOT WITH MULLER`S METHOD: X=',root, 'steps= ',n

      close(11)
      close(12)
  100 FORMAT('',A,F10.8,A,I5)
  200 FORMAT('',A,F10.8,A,ES11.2)
      end


      REAL*8 FUNCTION f_x(x)
      IMPLICIT NONE
      REAL*8 x
      f_x=cos(x)-x
      return
      end
      REAL*8 FUNCTION df_x(x)
      IMPLICIT NONE
      REAL*8 x
      df_x=-sin(x)-1
      return
      end

      SUBROUTINE SEQUENTIAL(root,x,eps,n)
      IMPLICIT NONE
      INTEGER out,n
      REAL*8 root, x , x_0, x_1, eps,f_x
      out=0

      x_0=x
      n=1
      do while(out.eq.0)
      	x_1=x_0+eps
      	if ((f_x(x_0)*f_x(x_1)).lt.0) then
      	out=1
      	root=x_0+(abs(x_1-x_0)/2d0)
      	WRITE(*,*)n,x_1
      	endif
      	x_0=x_1
      	n=n+1
      	
      enddo
      return
      end

      SUBROUTINE BISECTION(root,x1,x2,eps,n)
      IMPLICIT NONE
      INTEGER n
      REAL*8 root, x_0, x_1, x_2, eps,f_x,diff,x1,x2
c      WRITE(*,*)x1,x1

      x_0=x1
      x_2=x2
      n=1
      diff=1.0
      do while((diff.ge.eps).and.(n.le.100))
      	x_1=x_0+abs(x_2-x_0)/2.0
      	if ((f_x(x_0)*f_x(x_1)).lt.0) then
      		x_2=x_1
      	elseif ((f_x(x_1)*f_x(x_2)).lt.0) then
      		x_0=x_1
      	endif
      	diff=abs(x_2-x_0)
      	n=n+1
      	WRITE(12,*)n,abs(x_2-x_0),abs((x_2-x_0)/x_2)
c      	write(*,*)n,x_0,x_1,x_2
      enddo
      root=x_1
      return
      end

      SUBROUTINE SECANT(root,x1,x2,eps,n)
      IMPLICIT NONE
      INTEGER n
      REAL*8 root, x_0, x_1, x_2, eps,f_x,diff,x1,x2,x_2old
c      WRITE(*,*)x1,x1

      x_0=x1
      x_1=x2
      n=1
      diff=1.0
      do while((diff.ge.eps).and.(n.le.100))
      	x_2= x_1-f_x(x_1)*(x_1-x_0)/(f_x(x_1)-f_x(x_0))
      	n=n+1
      	if (n.ge.3) then
      		diff=abs(x_2old-x_2)
      		WRITE(12,*)n,abs(x_2-x_2old),abs((x_2-x_2old)/x_2)
      	endif
      	x_2old=x_2
      	x_0=x_1
      	x_1=x_2
c      	write(*,*)n,x_0,x_1,x_2
      enddo
      root=x_2
      return
      end
      SUBROUTINE REGULA_FALSI(root,x1,x2,eps,n)
      IMPLICIT NONE
      INTEGER n
      REAL*8 root, x_0, x_1, x_2, eps,f_x,diff,x1,x2,x_2old
c      WRITE(*,*)x1,x1

      x_0=x1
      x_1=x2
      n=1
      diff=1.0
      do while((diff.ge.eps).and.(n.le.100))
      	x_2= x_1-f_x(x_1)*(x_1-x_0)/(f_x(x_1)-f_x(x_0))
      	if ((f_x(x_0)*f_x(x_2)).gt.0) then
      		x_0=x_2
      	elseif ((f_x(x_1)*f_x(x_2)).gt.0) then
      		x_1=x_2
      	endif
      	n=n+1
      	if (n.ge.3) then
      		diff=abs(x_2old-x_2)
      		WRITE(12,*)n,abs(x_2-x_2old),abs((x_2-x_2old)/x_2)
      	endif
      	x_2old=x_2
c      	write(*,*)n,x_0,x_1,x_2
      enddo
      root=x_2
      return
      end
      SUBROUTINE NEWTON_RAPHSON(root,x1,eps,n)
      IMPLICIT NONE
      INTEGER n
      REAL*8 root, x_0, x_1, eps,f_x,df_x,diff,x1,x2
c      WRITE(*,*)x1,x1

      x_0=x1
      n=1
      diff=1.0
      do while((diff.ge.eps).and.(n.le.100))
      	x_1= x_0-f_x(X_0)/df_x(x_0)
      	diff=abs(x_1-x_0)
      	WRITE(12,*)n,abs(x_1-x_0),abs((x_1-x_0)/x_1)
      	n=n+1
      	x_0=x_1
c      	write(*,*)n,x_0,x_1,x_2
      enddo
      root=x_1
      return
      end
      SUBROUTINE MULLER(root,x1,x2,x3,eps,n)
      IMPLICIT NONE
      INTEGER n
      REAL*8 root, x_0, x_1, x_2, eps,f_x,diff,x1,x2,x3,a,b
      REAL*8 disc,x_31,x_32,x_3,x_2old
c      WRITE(*,*)x1,x1

      x_0=x1
      x_1=x2
      x_2=x3
      n=1
      diff=1.0
      do while((diff.ge.eps).and.(n.le.100))
      	a=((x_1-x_2)*(f_x(x_1)-f_x(x_2))-(x_0-x_2)*(f_x(x_1)-f_x(x_2)))/
     .((x_0-x_1)*(x_0-x_2)*(x_1-x_2))
      	b=((x_0-x_2)**2*(f_x(x_1)-f_x(x_2))-(x_1-x_2)**2
     .*(f_x(x_0)-f_x(x_2)))/((x_0-x_1)*(x_0-x_2)*(x_1-x_2))
      	disc=b**2-4*a*f_x(x_2)
c      	WRITE(*,*)a,b,disc
      	if (disc.le.0)then
      		WRITE(*,*) 'MULLER METHOD ERROR: COMPLEX SOLUTION FOUND'
      		exit
      	elseif (disc.gt.0) then
      		x_31=x_2+(-b+sqrt(disc))/(2*a)
      		x_32=x_2+(-b-sqrt(disc))/(2*a)
c      		WRITE(*,*)x_31,x_32
c      		if ((x_31.gt.x_0).and.(x_31.lt.x_1))then
c      			x_3=x_31
c      		elseif ((x_32.gt.x_0).and.(x_32.lt.x_1))then
c      			x_3=x_32
c      		else
c      			WRITE(*,*) 'MULLER METHOD ERROR: COMPUTED POSITION 
c     .OUT OF RANGE'
c      			exit
c      		endif
            x_3=max(x_31,x_32)
c      		WRITE(*,*)x_31,x_32,x_3
      		if ((x_3.lt.x_2).and.((f_x(x_0)*f_x(x_3))).gt.0)then
      			x_0=x_3
      		elseif ((x_3.gt.x_2).and.((f_x(x_0)*f_x(x_3))).gt.0)then
      			x_0=x_2
      			x_2=x_3
      		elseif ((x_3.gt.x_2).and.((f_x(x_1)*f_x(x_3))).lt.0)then
      			x_1=x_3
      		endif
c      		x_0=x_1
c      		x_1=x_2
c      		x_2=x_3
      	endif
      	n=n+1
      	if (n.ge.3) then
      		diff=abs(x_2old-x_2)
      		WRITE(12,*)n,abs(x_2-x_2old),abs((x_2-x_2old)/x_2)
      	endif
      	x_2old=x_2
c      	write(*,*)n,x_0,x_1,x_2
      enddo
      root=x_3
      return
      end

      SUBROUTINE MULLER1(root,x1,x2,x3,eps,n)
      IMPLICIT NONE
      INTEGER n
      REAL*8 root, x_0, x_1, x_2, eps,f_x,diff,x1,x2,x3,a,b
      REAL*8 disc,x_31,x_32,x_3
c      WRITE(*,*)x1,x1

      x_0=x1
      x_1=x2
      x_2=x3
      n=0
      diff=1.0
      do while((diff.ge.eps).and.(n.le.100))
      	a=((x_1-x_2)*(f_x(x_1)-f_x(x_2))-(x_0-x_2)*(f_x(x_1)-f_x(x_2)))/
     .((x_0-x_1)*(x_0-x_2)*(x_1-x_2))
      	b=((x_0-x_2)**2*(f_x(x_1)-f_x(x_2))-(x_1-x_2)**2
     .*(f_x(x_0)-f_x(x_2)))/((x_0-x_1)*(x_0-x_2)*(x_1-x_2))
      	disc=b**2-4*a*f_x(x_2)
      	WRITE(*,*)a,b,disc
      	if (disc.le.0)then
      		WRITE(*,*) 'MULLER METHOD ERROR: COMPLEX SOLUTION FOUND'
      		exit
      	elseif (disc.gt.0) then
      		x_31=x_2+(-b+sqrt(disc))/(2*a)
      		x_32=x_2+(-b-sqrt(disc))/(2*a)
      		x_3=max(x_31,x_32)
      		WRITE(*,*)x_31,x_32,x_3
c      		if ((x_31.gt.x_0).and.(x_31.lt.x_1))then
c      			x_3=x_31
c      		elseif ((x_32.gt.x_0).and.(x_32.lt.x_1))then
c      			x_3=x_32
c      		else
c      			WRITE(*,*) 'MULLER METHOD ERROR: COMPUTED POSITION 
c     .OUT OF RANGE'
c      			exit
c      		endif
c      		if (x_3.lt.x_2)then
c      			x_0=x_3
      		x_0=x_1
      		x_1=x_2
      		x_2=x_3
      	endif
      	diff=abs(x_1-x_0)
      	n=n+1
c      	write(*,*)n,x_0,x_1,x_2
      enddo
      root=x_3
      return
      end
