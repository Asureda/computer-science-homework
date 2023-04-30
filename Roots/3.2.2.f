      PROGRAM ROOT_METHODS
      IMPLICIT NONE
      INTEGER n
      REAL*8 root

      CALL SECANT(root,25d-2,5d-1,1d-8,n)
      WRITE(*,*)'ROOT WITH SECANT METHOD: X=',root, 'steps= ',n
      CALL  NEWTON_RAPHSON(root,5d-1,1d-8,n)
      WRITE(*,*)'ROOT WITH NEWTON_RAPHSON METHOD: X=',root, 'steps= ',n
      end

      REAL*8 FUNCTION f_x(x)
      IMPLICIT NONE
      REAL*8 x
      f_x=x*sin(x)-cos(x)
      return
      end
      REAL*8 FUNCTION df_x(x)
      IMPLICIT NONE
      REAL*8 x
      df_x=sin(x)+x*cos(x)+sin(x)
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
c           write(*,*)n,x_0,x_1,x_2
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
c           write(*,*)n,x_0,x_1,x_2
      enddo
      root=x_1
      return
      end