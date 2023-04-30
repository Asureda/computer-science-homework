      PROGRAM EPSILON
      implicit none
      REAL*4 eps4
      REAL*8 eps8

      eps4=1.0
      do while ((eps4 +1.0).gt.1.0)
      	eps4=eps4/2.0
      enddo
      WRITE(*,*)'REAL*4', eps4

      eps8=1.0
      do while ((eps8 +1.0).gt.1.0)
      	eps8=eps8/2.0
      enddo
      WRITE(*,*)'REAL*8', eps8
      end