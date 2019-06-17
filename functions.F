c***********************************************************************
c
c     The subroutines in this file were taken from the T-matrix code 
c     written by Michael Mishchenko for axisymmetric particles 
c     (M. Mishchenko, JQSRT 60, 309-324, 1998). Thanks to 
c     Michal Mishchenko for helpful discussions, and for his kind 
c     permission to include these routines in the distribution of Tsym!    
c
c***********************************************************************



c***********************************************************************
c     Calculation of the Wigner d-functions
c     dv1(n)=dvig(0,m,n,arccos x)
c     and their derivatives
c     dv2(n)=[d/d(arccos x)] dvig(0,m,n,arccos x)
c     1.le.n.le.nmmax
c     -1.le.x.le.1
c***********************************************************************
      subroutine vig (x0, nmmax, m, dv1, dv2)
      implicit real*8 (a-h,o-z)
      integer nmmax,m
      real*8 x0,dv1(nmmax), dv2(nmmax),x,fac

      x=dabs(x0)

      a=1d0
      qs=dsqrt(1d0-x*x)
      qs1=1d0/qs
      do n=1,nmmax
         dv1(n)=0d0
         dv2(n)=0d0
      enddo
      if (m.eq.0) then
         d1=1d0
         d2=x
         do n=1,nmmax
            qn=dble(n)
            qn1=dble(n+1)
            qn2=dble(2*n+1)
            d3=(qn2*x*d2-qn*d1)/qn1
            der=qs1*(qn1*qn/qn2)*(-d1+d3)
            dv1(n)=d2
            dv2(n)=der
            d1=d2
            d2=d3
         enddo
      else
         qmm=dble(m*m)
         do i=1,m
            i2=i*2
            a=a*dsqrt(dble(i2-1)/dble(i2))*qs
         enddo
         d1=0d0
         d2=a
         do n=m,nmmax
            qn=dble(n)
            qn2=dble(2*n+1)
            qn1=dble(n+1)
            qnm=dsqrt(qn*qn-qmm)
            qnm1=dsqrt(qn1*qn1-qmm)
            d3=(qn2*x*d2-qnm*d1)/qnm1
            der=qs1*(-qn1*qnm*d1+qn*qnm1*d3)/qn2
            dv1(n)=d2
            dv2(n)=der
            d1=d2
            d2=d3
         enddo
      endif
      
      if(x0.lt.0d0)then
         do n=max(1,m),nmmax
            fac=(-1d0)**(n-m)
            dv1(n)=fac*dv1(n)
            dv2(n)=-fac*dv2(n)
         enddo
      endif

      end



c***********************************************************************
C     CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE FIRST KIND J        
C     OF REAL-VALUED ARGUMENT X OF ORDERS FROM 1 TO NMMAX BY USING         
C     BACKWARD RECURSION. PARAMETER NNMAX DETERMINES NUMERICAL ACCURACY.   
C     U - FUNCTION (1/X)(D/DX)(X*J(X))                                    
C***********************************************************************
      SUBROUTINE RJB(X,Y,U,nmmax)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer nmmax
      REAL*8 Y(nmmax),U(nmmax),Z(1000)
      NNMAX=12+INT(1.2d0*X)
      L=nmmax+NNMAX
      if(L.gt.1000)then
         write(6,*)'ERROR in RJB: Z-array too small'
         stop
      endif
      XX=1d0/X
      Z(L)=1d0/(dble(2*L+1)*XX)
      L1=L-1
      DO 5 I=1,L1
         I1=L-I
         Z(I1)=1d0/(dble(2*I1+1)*XX-Z(I1+1))
    5 CONTINUE
      Z0=1d0/(XX-Z(1))
      Y0=Z0*dCOS(X)*XX
      Y1=Y0*Z(1)
      U(1)=Y0-Y1*XX
      Y(1)=Y1
      DO 10 I=2,nmmax
         YI1=Y(I-1)
         YI=YI1*Z(I)
         U(I)=YI1-dble(I)*YI*XX
         Y(I)=YI
   10 CONTINUE
      RETURN
      END

C************************************************************************
C     CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE SECOND KIND Y        
C     OF REAL-VALUED ARGUMENT X OF ORDERS FROM 1 TO nmmax BY USING FORWARD  
C     RECURSION. V - FUNCTION (1/X)(D/DX)(X*Y(X))          
C************************************************************************
      SUBROUTINE RYB(X,Y,V,nmmax)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer nmmax
      REAL*8 Y(nmmax),V(nmmax)
      C=dCOS(X)
      S=dSIN(X)
      X1=1d0/X
      X2=X1*X1
      X3=X2*X1
      Y1=-C*X2-S*X1
      Y(1)=Y1
      Y(2)=(-3d0*X3+X1)*C-3d0*X2*S
      NMAX1=nmmax-1
      DO 5 I=2,NMAX1
    5     Y(I+1)=dble(2*I+1)*X1*Y(I)-Y(I-1)
      V(1)=-X1*(C+Y1)
      DO 10 I=2,nmmax
  10       V(I)=Y(I-1)-dble(I)*X1*Y(I)
      RETURN
      END

C***********************************************************************
C     CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE FIRST KIND          
C     Y=J(Z) OF COMPLEX ARGUMENT Z OF ORDERS FROM 1 TO nmmax BY USING      
C     BACKWARD RECURSION.                                                  
C     U=U(Z) - FUNCTION (1/Z)(D/DZ)(Z*J(Z))                                
C***********************************************************************
      SUBROUTINE CJB (CX,CY,CU,nmmax)
      IMPLICIT NONE
      integer nmmax,l,nnmax,l1,i,i1
      COMPLEX*16 CY(nmmax),CZ(1000),CU(nmmax)
      COMPLEX*16 CX,C1,CXX,CZ0,CY0,CY1,CU1,CYI1,CYI,CUI

      NNMAX=12+INT( 1.3d0*ABS(CX) )      
      L=nmmax+NNMAX
      if(L.gt.1000.or.L.lt.0)then
         write(6,*)'ERROR in CJB: CZ-array dimension error, L=',L
         write(6,*)'nnmax,nmmax,cx:',nnmax,nmmax,cx
         stop
      endif
      C1=(0d0,-1d0)
      CXX=1d0/CX
      CZ(L)=1d0/(dble(2*L+1)*CXX)
      L1=L-1
      DO 5 I=1,L1
         I1=L-I
         CZ(I1)=1d0/(dble(2*I1+1)*CXX-CZ(I1+1))
    5 CONTINUE
      CZ0=1d0/(CXX-CZ(1))
      CY0=CZ0*CdCOS(CX)*CXX
      CY1=CY0*CZ(1)
      CU1=CY0-CY1*CXX
      CY(1)=CY1
      CU(1)=CU1
      DO 10 I=2,nmmax
         CYI1=CY(I-1)
         CYI=CYI1*CZ(I)
         CUI=CYI1-dble(I)*CYI*CXX
         CY(I)=CYI
         CU(I)=CUI
   10 CONTINUE
      RETURN
      END


C***********************************************************************
C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE          *
C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON       *
C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED.                *
C    N - NUMBER OF POINTS                                              *
C    Z - DIVISION POINTS                                               *
C    W - WEIGHTS                                                       *
C*********************************************************************** 
      SUBROUTINE GAUSS ( N,IND1,IND2,Z,W )
      IMPLICIT REAL*8 (A-H,P-Z)
      REAL*8 Z(N),W(N)
      DATA A,B,C /1D0,2D0,3D0/
      IND=MOD(N,2)
      K=N/2+IND
      F=DFLOAT(N)
      DO 100 I=1,K
          M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0D0
          NITER=0
          CHECK=1D-16
   10     PB=1D0
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
          CHECK=CHECK*10D0
   15     PC=X
          DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          IF(DABS(PB).GT.check*DABS(X)) GO TO 10
          Z(M)=X
          W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M)
          IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
          Z(I)=-Z(M)
          W(I)=W(M)
  100 CONTINUE
      IF(IND2.NE.1) GO TO 110
      PRINT 1100,N
 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',
     * ' OF ',I4,'-TH ORDER')
      DO 105 I=1,K
          ZZ=-Z(I)
  105     PRINT 1200,I,ZZ,I,W(I)
 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
      GO TO 115
  110 CONTINUE
C     PRINT 1300,N
 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
  115 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO 120 I=1,N
  120     Z(I)=(A+Z(I))/B
  140 CONTINUE
      RETURN
      END
 


