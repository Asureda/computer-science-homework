      PROGRAM TAYLOR_EXP
      IMPLICIT NONE
      INTEGER n
      REAL*8 x,result,TAY_EXP

      x=2d0
      n=15

      result=TAY_EXP(x,n)

      WRITE(*,*)'VALUE OF EXP(X) IN TAYLOR SERIES: exp(x)=',result
      end

      REAL*8 FUNCTION TAY_EXP(x,n)
      IMPLICIT NONE
      INTEGER n,i
      REAL*8 x,sum,ord

      sum=1

      if (n.ge.1) then
      	do i=1,n
      		sum=sum+ord(x,i)
      	enddo
      endif
      TAY_EXP=sum
      return
      end

      REAL*8 FUNCTION ord(x,i)
      IMPLICIT NONE
      INTEGER i,k
      REAL*8 sum,x,FACTORIAL

      sum=1
      do k=1,i
      	sum=sum*x
      enddo
      ord=sum/FACTORIAL(i)
      return
      end

      REAL*8 FUNCTION FACTORIAL(n)
      IMPLICIT NONE
      INTEGER n,k
      REAL*8 ans
      ans=1

      if (n.eq.0)then
      	continue
      elseif(n.ge.1)then
      	do k=1,n
      		ans=ans*k
      	enddo
      endif
      FACTORIAL=ans
      return
      end
