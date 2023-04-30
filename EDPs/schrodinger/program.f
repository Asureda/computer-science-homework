      PROGRAM PRACTICA11




      IMPLICIT NONE
      INTEGER nx,i,TEMPS,nidex,tempsmax,ntemps
      DOUBLE PRECISION Lx,Lx1,h,h1,font,deltat,k1,r1,promig,potencial,
     #gausiana,x_tranmiss,T,Ref,norm
      DOUBLE COMPLEX eye,k,r
      PARAMETER(eye=(0.d0,1.d0))
      PARAMETER (Lx=2*17.d0)
      PARAMETER (nx=200)
      PARAMETER (h=Lx/dble(nx))
      PARAMETER (tempsmax=20)
      PARAMETER (ntemps=9000)
      PARAMETER (deltat=tempsmax/dble(ntemps))
      PARAMETER (k=eye/4.d0)
      PARAMETER (r=k*deltat/h**2.d0)
      DOUBLE COMPLEX Psiinicial(0:nx),a(nx+1),b(nx+1),c(nx+1)
      DOUBLE COMPLEX Psinew(nx+1), Taux(nx+1),Taux1(nx+1)


      COMMON/PARAMETRES/Lx1,h1
      OPEN(10,FILE='program-res1.dat')
      OPEN(11,FILE='program-coeficients.dat')
      Lx1=Lx
      h1=h
C     Construim la funció dona inicial, gaussiana amb moment
      WRITE(*,*)'CALL CONTORN'
      CALL CONTORN(nx,Psiinicial)
C     Reescrivim vector per per coincidir indexos
      norm=0d0
      DO i=1,nx
          norm=norm+h*(abs(Psiinicial(i)))**2
      enddo
      Psiinicial=Psiinicial/sqrt(norm)
      Do i=1,nx+1
            WRITE(10,*)0,(-Lx/2)+h*(i-1),(abs(Psiinicial(i-1)))**2
            Taux1(i)=Psiinicial(i-1)
      enddo
      WRITE(10,*)
      WRITE(10,*)
C     Creem els tre vectors de la esquerra per a la subrutina
      DO i=1,nx+1
            a(i)=-r
            b(i)=1+(eye*deltat*potencial(-(Lx/2.)+(i-1)*h)/2.)+2*r
            c(i)=-r
      enddo

      a(1)=0
      c(nx+1)=0
C     Comensem bucle de temps per evolucio temporal
      DO TEMPS=1,ntemps
            !WRITE(*,*)'pas de temps', TEMPS
            Taux=0
            Psinew=0
C     CALCULEM A CADA PAS EL VECTOR DE LA DRETA PER A LA SUBRUTINA TRIDIAG

            DO i=2,nx
                  Taux(i)=(1-(eye*deltat*potencial(-(Lx/2.)+(i-2)*h)/
     #2.)-2*r)*Taux1(i)+r*(Taux1(i-1)+Taux1(i+1))
            enddo
            Taux(1)=(1-(eye*deltat*potencial(-Lx/2.)/2.)-2*r)*Taux1(1)
     #+r*(Taux1(2))+2*r*Psiinicial(0)
            Taux(nx+1)=(1-(eye*deltat*potencial(Lx/2.)/2.)-2*r)*
     #Taux1(nx+1)+r*(Taux1(nx))+2*r*Psiinicial(nx)
C           Cridem subrutina tridiag
            CALL TRIDIAG(A,B,C,Taux,Psinew,nx+1)
            !Normalitzacio de la funció d'ona
            norm=0d0
            DO i=1,nx
                  norm=norm+h*(abs(Psinew(i)))**2
            enddo
            Psinew=Psinew/sqrt(norm)
C           Escrivim dedes        
            if(mod(temps,100).eq.0) then
            !WRITE(*,*)'WRITE'            
            DO i=1,nx
                  WRITE(10,*) TEMPS*deltat,(-Lx/2.)+h*(i-1),
     #(abs(Psinew(i)))**2
            enddo
            WRITE(10,*)
            WRITE(10,*)
          ENDIF
          T=0d0
            Ref=0d0
            DO i=1,nx
            x_tranmiss=(-Lx/2.)+h*(i-1)
            
            if(x_tranmiss.lt.0d0)then
                T=T+h*(abs(Psinew(i)))**2
            elseif(x_tranmiss.gt.0d0)then
                Ref=Ref+h*(abs(Psinew(i)))**2
            ENDIF
           END DO
            WRITE(11,*)TEMPS*deltat,T,Ref,T+Ref
C     Imposem les condicions de contorn de la barrera de potencial.
            Psinew(1)=Psiinicial(0)
            Psinew(nx+1)=Psiinicial(nx)
C     Fem el canvi de psi vella a la nova
            do i=1,nx+1
            Taux1(i)=Psinew(i)
            enddo

      enddo
      OPEN(15,FILE='potential.dat')
      DO i=1,nx
        x_tranmiss=(-Lx/2.)+h*(i-1)
        WRITE(15,*)x_tranmiss,potencial(x_tranmiss)
      END DO


      end


      SUBROUTINE CONTORN(dim,Tempini)
      INTEGER dim,i
      DOUBLE COMPLEX Tempini(0:dim),gausiana
      DOUBLE PRECISION Lx1,h1
      COMMON/PARAMETRES/Lx1,h1


      do i=0,dim
            Tempini(i)=gausiana((-Lx1/2.d0)+i*h1)
      enddo
      Tempini(0)=0
      Tempini(dim)=0
      return
      end


      DOUBLE COMPLEX FUNCTION gausiana(x)
      DOUBLE COMPLEX gaussiana,eye,aux1,aux2
      DOUBLE PRECISION x,pi
      eye=(0.d0,1.d0)
      pi=acos(-1.d0)
      aux1=(eye*200.d0*x)/17.d0
      aux2=((x-7.d0)**2.d0)/4.d0
      gausiana=(1.d0/(2*pi)**(1./4.))*exp(-aux2)*exp(-aux1)

      return
      end

      DOUBLE PRECISION  FUNCTION potencial(x)
      DOUBLE PRECISION x, potencial, pi
      pi=acos(-1.)
      potencial=0d0
      IF((x.gt.-2.5).and.(x.lt.2.5))then
        potencial=50.0
      ENDIF
      return
      end


      SUBROUTINE TRIDIAG(A,B,C,R,PSI,IMAX)
c Solves the problem T psi =R
c 
c where T is a tridiagonal matrix, A (lower), B (central), C (upper)
c  A   0 a1, a2, a3, ...., aIMAX
c  B   b1 b2, b3, b4, ...., bIMAX
c  C   c1, c2, c3, c4, ....,0 
        
c       IMPLICIT double precision (A-H,K,O-Z)
c       IMPLICIT INTEGER (I-J , L-N)
        double complex  BET
        double complex GAM(4001)
        double complex A(IMAX),B(IMAX),C(IMAX),R(IMAX),
     #PSI(IMAX)

        IF(B(1).EQ.0.) STOP
        BET=B(1)
        PSI(1)=R(1)/BET
        DO 11 J=2,IMAX
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
        IF(BET.EQ.0) EXIT
        PSI(J)=(R(J)-A(J)*PSI(J-1))/BET
11      CONTINUE

        DO 12 J=IMAX-1,1,-1
        PSI(J)=PSI(J)-GAM(J+1)*PSI(J+1)
12      CONTINUE

       RETURN
       END