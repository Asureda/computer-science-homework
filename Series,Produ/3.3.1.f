      PROGRAM NUMBER_E

      IMPLICIT NONE
      INTEGER n,n_max
      REAL*8 e,precision,e_old,eps,e_new,FACTORIAL
      PARAMETER(n_max=15)

      WRITE(*,*)'INPUT A PRESICION TO COPUTE THE NUMBER e'
      READ(*,*) precision
      WRITE(*,*)'THE PRESICION IS:', precision

      e_old=0
      n=0
      eps=2d0
      do while((eps.ge.precision).and.(n.le.n_max))

      	e_new=e_old+(1d0/(1d0*FACTORIAL(n)))
      	eps=abs(e_old-e_new)
      	WRITE(*,*)FACTORIAL(n),eps,precision,e_new,e_old,n,n_max

      	e_old=e_new
      	n=n+1
      enddo
      WRITE(*,*) 'THE VALUE FOR IS e=',e_new, eps

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
