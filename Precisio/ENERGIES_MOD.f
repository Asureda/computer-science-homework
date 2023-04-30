      PROGRAM EIGEN_ENERGIES
      IMPLICIT NONE
      INTEGER j,k,n_e,n
      REAL*8 x, etha, pi, a, b, E,result,up,down
      REAL*8 e_min,e_max,h_e,a_res
      REAL*8 x_0,x_1,x_2,f_0,f_1,eps,f_2,f_x,h
      open(11,FILE='integral_conv.dat')
      open(12,FILE='result.dat')

      pi=acos(-1.0)
      write(*,*) pi
      a=0.0
      b=50.0
      

c-----------------------------------------------------------------------
c     Convergence of the integral: finding the best n
C
      do j=10,1000
      E=0.6
      CALL rectan(result,a,b,j,(-E*0.5+3.0/3.4))
      WRITE(11,*)result,j,(-E*0.5+3.0/3.4)
      enddo
C-----------------------------------------------------------------------
C     Root of the function
c     C
      result=0
      E=0.1
c      etha=0.5
c----------------------------------
      e_min=-0.05
      e_max=0.7
      n_e=10000
      h_e=(e_max-e_min)/n_e
      x_0=-20.0
      x_1=21.0     
      do k=0,n_e
      	E=e_min+k*h_e
      CALL rectan(up,a,b,2000,(-E*0.5+3.0/4.0))
      CALL rectan(down,a,b,2000,(-E*0.5+1.0/4.0))
      n=0
      eps=1.0
      do while((eps>0.0001).and.(n<50))
      	x_2=x_0+(x_1-x_0)/2
      	f_0=(1.0/x_0)-(up/down)
      	f_1=(1.0/x_1)-(up/down)
      	f_2=(1.0/x_2)-(up/down)
      	if(f_0*f_2.le.0) then
      		f_1=f_2
      	elseif (f_2*f_1.le.0)then
      		f_0=f_2
      	endif


c      	x_2=x_1-(f_1*(x_1-x_0)/(f_1-f_0))
      	write(*,*)n,x_0,x_2,x_1
      	write(*,*)'--',f_0,f_2,f_1
      	eps=x_2-x_1
c      	x_1=x_2
      	n=n+1
      	
      enddo
      WRITE(12,*)x_2,E
c      WRITE(*,*)'(a,E)',x_2,',',E
      enddo

      end






      SUBROUTINE rectan(result,a,b,n,t)
      IMPLICIT NONE
      INTEGER n,i
      real*8 suma,y_vec,result,h,t,a,b
      real*8 f_x,aux
      suma=0
      result=0
      h=(b-a)/n
c      OPEN(20,FILE='rectan_func.dat')
      do i=1,n
c      	WRITE(*,*)i,h,(a+i*h),(f_x((a+i*h),E,etha)),E,etha
            suma=suma+f_x(a+i*h,t)
      enddo
      result=suma*h
      return
      end

      REAL*8 FUNCTION f_x(x,t)
      IMPLICIT NONE
      REAL*8 x,t
      f_x= (x**(t-1))*exp(-x)
      return
      end