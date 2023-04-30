      PROGRAM EXP
      IMPLICIT NONE
      integer n,i,j
      REAL*4 en1,en_1,en10,en10_1,x
      x=1.0

      en_1=1.0/(1.0*dexp(1d0))
      WRITE(*,*) 'VALOR INICIAL',en_1
      n=0
      do while (en_1>0)
      	n=n+1
      	en1=1-n*en_1
      	en_1=en1
      enddo
      WRITE(*,*) 'ITERACIO',n,'VALOR FINAL',en1

      OPEN(11,FILE='en_conv.dat')
      OPEN(12,FILE='e_conv.dat')

      do j=8,15

      en_1=1d0
c      WRITE(*,*) j
      n=j
      do i=2,j
      	en1=(1d0-en_1)/(dble(n))
c      	WRITE(*,*)n,'VALOR',en1
      	n=n-1
      	en_1=en1
            
      enddo
      WRITE(*,*) j,en1,1d0/en1
      enddo
      WRITE(*,*) 'final value for integral: E=',en1
      WRITE(*,*) 'final value for e: e=',1d0/en1

      
      end