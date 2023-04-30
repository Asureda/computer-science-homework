      PROGRAM CUADRATIC_SOLUTIONS
      IMPLICIT NONE
      INTEGER deep,k
      REAL*8 b,c,dis,eps,root1,root2,x1,x2

      WRITE(*,*)'PROGRAM TO SOLVE CUADRATIC EQUATIONS x**2+bx+c=0'
      WRITE(*,*)'INPUT A VALUE FOR b'
c      READ(*,*) b
      WRITE(*,*)'INPUT A VALUE FOR c'
c      READ(*,*) c
      b=4.0
      c=-1.0

      dis=(b**2)-4*c
      eps=0.001
      deep=20

      if ((dis.gt.0).and.(dis.lt.eps))then
      	WRITE(*,*)'ONE ROOT FOR THE EQUATION'
      	WRITE(*,*) 'x=',-b/2.d0
      elseif(dis.ge.eps) then
      	WRITE(*,*)'TWO ROOTS FOR THE FINTION'
      	x1=root1(b,c,deep)
      	x2=root2(b,c,deep)
      	WRITE(*,*) 'x1=', x2 , 'exact=',-2+sqrt(5.0)
      	WRITE(*,*) 'x2=', x1 , 'exact=',-2-sqrt(5.0)
      else
      	WRITE(*,*)'NO REAL ROOTS FOR THE EQUATION'
      endif
      end

      REAL*8 RECURSIVE FUNCTION root1(b,c,n) result(res)
      IMPLICIT NONE
      INTEGER n
      REAL*8 b,c,res
c      WRITE(*,*) 'enter function root1',n
      res=1
      if(n>0) then
      	res=-b-c/(1.0*root1(b,c,n-1))
c      	WRITE(*,*) 'root1',n
      else
      	res=-b-c/(-b)
      endif
      return
      end

      REAL*8 RECURSIVE FUNCTION root2(b,c,n) result(res)
      IMPLICIT NONE
      INTEGER n
      REAL*8 b,c,res
c      WRITE(*,*) 'enter function root2',n
      res=1

      if(n>0) then
      	res=-c/(b+root2(b,c,n-1))
c      	WRITE(*,*) 'root2',n
      else
      	res=c/b
      endif
      return
      end