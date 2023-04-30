      PROGRAM EIGEN_ENERGIES
      IMPLICIT NONE
      INTEGER j,k,n_e,n
      REAL*8 x, etha, pi, a, b, E,result
      REAL*8 e_min,e_max,h_e,a_res
      REAL*8 x_0,x_1,x_2,f_0,f_1,eps,f_2,f_x,h
      open(11,FILE='integral_conv.dat')
      open(12,FILE='result.dat')

c-----------------------------------------------------------------------
c     Plot of the function
      a=0
      b=20
      n=500
      h=(b-a)/n

      E=0.1
      etha=0.5
      
      do j=1,n-1
      	WRITE(11,*)(j*h),(f_x(j*h,E,etha))
c      	WRITE(*,*)j,h,j*h,E,etha
      enddo
      WRITE(11,*)
      WRITE(11,*)



      pi=acos(-1.0)
      write(*,*) pi
      a=0.0
      b=50.0
      

c-----------------------------------------------------------------------
c     Convergence of the integral: finding the best n
C
      do j=10,1000
      CALL rectan(result,a,b,j,E*0.5,etha)
      WRITE(11,*)result,j
      enddo
C-----------------------------------------------------------------------
C     Root of the function
c     C
      result=0
      E=0.1
c      etha=0.5
c----------------------------------
c      e_min=-0.05
c      e_max=0.7
c      n_e=10000
c      h_e=(e_max-e_min)/n_e
      etha=0.5
      x_0=0.5
      x_1=10.5      
c      do k=0,n_e
c      	E=e_min+k*h_e
      CALL rectan(result,a,b,5,E*0.5,etha)
      n=0
      do while((eps>0.0001).or.(n<10))
      	x_2=x_0+(x_1-x_0)/2
      	f_0=(-sqrt(2*pi)/x_0)-result
      	f_1=(-sqrt(2*pi)/x_1)-result
      	f_2=(-sqrt(2*pi)/x_2)-result
      	if(f_0*f_2.le.0) then
      		f_1=f_2
      	elseif (f_2*f_1.le.0)then
      		f_0=f_2
      	endif


c      	x_2=x_1-(f_1*(x_1-x_0)/(f_1-f_0))
      	write(*,*)n,x_0,x_2,x_1,result
      	write(*,*)'--',f_0,f_2,f_1,(-(sqrt(2*pi)/result))
      	eps=x_2-x_1
c      	x_1=x_2
      	n=n+1
      	
      enddo
c      WRITE(12,*)x_2,E
c      WRITE(*,*)'(a,E)',x_2,',',E
c      enddo

      WRITE(12,*)
      WRITE(12,*)





C-----------------------------------------------------------------------
C     For a given etha, finding E(a)
c      Best n is n=40 (aprox)
      e_min=-0.05
      e_max=0.7
      n_e=10000
      h_e=(e_max-e_min)/n_e
      etha=0.5
      do k=0,n_e
      	E=e_min+k*h_e
c      	CALL rectan(result,a,b,40,-E*0.5,etha)
      	a_res=-(sqrt(2*pi)/result)
c      	WRITE(*,*)a_res,E,result
      	WRITE(12,*)a_res,E
      enddo
      end






      SUBROUTINE rectan(result,a,b,n,E,etha)
      IMPLICIT NONE
      INTEGER n,i
      real*8 suma,y_vec,result,h,E,etha,a,b
      real*8 f_x,aux
      suma=0
      result=0
      h=(b-a)/n
c      OPEN(20,FILE='rectan_func.dat')
      do i=1,n
c      	WRITE(*,*)i,h,(a+i*h),(f_x((a+i*h),E,etha)),E,etha
            suma=suma+f_x(a+i*h,E,etha)
      enddo
      result=suma*h
      return
      end

      REAL*8 FUNCTION f_x(x,E,etha)
      IMPLICIT NONE
      REAL*8 x,E,etha
      f_x=((etha*exp(-E*x))/(sqrt(1-exp(-x))*(1-exp(-etha*x))))
     #-(1/(x**(3/2)))
      return
      end