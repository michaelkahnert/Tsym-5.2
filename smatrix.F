
c***********************************************************************  
c The subroutines in this file were taken from the T-matrix code written 
c by Heikki Laitinen and Kari Lumme for general star-shaped particles
c (JQSRT 60, 325-334, 1998). Thanks to Heikki Laitinen
c for some helpful discussions on how to use the code. Thanks also to
c Kari Lumme for his kind permission to include these
c subroutines in the distribution of Tsym!
c***********************************************************************  



c***********************************************************************  
C  Compute the amplitude scattering matrix for given incidence and
C  scattering directions using previously computed T-matrix. 
C
C  Arguments:
C    NAME    I/O  TYPE        DESCRIPTION
C    Tmat     I   complex*16  T-matrix, 
C    nmax     I   integer     maximum of index n used in T
C    lmax     I   integer     dimension of the blocks T11 etc.
C    thi      I   real*8      angle theta of the incindence direction
C    phi      I   real*8      angle phi of the incidence direction
C    ths      I   real*8      angle theta of the scattering direction
C    phs      I   real*8      angle phi of the scattering direction
C    S        O   complex*16  amplitude scattering matrix, dimension (2,2)
C
C  Local variables:
C    n,nn,m,mm,l,ll -- index variables
C    nmin,nnmin     -- lower limits for variables n and nn
c    dvx1(theta) = d_{0,m}^n(theta)/sin(theta)
c    dvx2(z) = (d/d theta) d_{0,m}^n(theta), 
c    where: d_{0,m}^n(theta) : Wigner d-function, and
c           x =i,s denotes:
c              i: theta=theta_i=polar angle of incident radiation
c              s: theta=theta_s=polar angle of scattered radiation
c    
c
c  NOTE: current implementation uses Cn-symmetry of T-matrix only,
c        but not (yet) DIHEDRAL- and SIGMAh-symmetry
c***********************************************************************
      subroutine AmplScat(Tmat,nmax,lmax,thi,phi,ths,phs,
     &                    C_N,DIHEDRAL,S)
      implicit none

      integer n,m,nn,mm,l,ll,nmin,nnmin,nmax,lmax,pidx,ppidx,C_N
      real*8 thi,phi,ths,phs,dvs1(nmax),dvs2(nmax),dvi1(nmax),
     &     dvi2(nmax),fr,dn,sign,mdvs1,mdvi1
      complex*16 Tmat(2,lmax,2,lmax),S(2,2),fc,fc1,fc2,tmp1,tmp2,tmp3,
     &     tmp4,eim,Im,czero
      parameter (czero=(0d0,0d0),Im=(0d0,1d0))
      logical Cn,DIHEDRAL
      
      Cn=.false.
      if (C_N.ge.3) Cn=.true.
      
c-----------------------------------------------------------------------
c     Amplitude scattering matrix:
c-----------------------------------------------------------------------
      S(1,1)=czero
      S(1,2)=czero
      S(2,1)=czero
      S(2,2)=czero
      do 10 m=-nmax,nmax
         nmin=max(abs(m),1)
         call vigampl(ths,nmax,abs(m),dvs1,dvs2)
         if (m.lt.0) then
            do 100 n=nmin,nmax
               dvs1(n)=(-1d0)**(m)*dvs1(n)
               dvs2(n)=(-1d0)**(m)*dvs2(n)
 100        continue
         end if
         do 9 mm=-nmax,nmax
            if (Cn.and.mod(mm-m,C_N).ne.0) goto 991
            nnmin=max(abs(mm),1)
            sign=(-1d0)**(m+mm)
            eim=exp(Im*(m*phs-mm*phi))
            call vigampl(thi,nmax,abs(mm),dvi1,dvi2)
            if (mm.lt.0) then
               do 200 nn=nnmin,nmax
                  dvi1(nn)=(-1d0)**(mm)*dvi1(nn)
                  dvi2(nn)=(-1d0)**(mm)*dvi2(nn)
 200           continue
            end if
            do 8 n=nmin,nmax
               l=n*(n+1)+m
               dn=dble(2*n+1)/dble(n*(n+1))
               mdvs1=m*dvs1(n)
               do 7 nn=nnmin,nmax
                  ll=nn*(nn+1)+mm
                  fr=sign*sqrt( dn*dble(2*nn+1)/dble(nn*(nn+1)) )
                  fc=Im**(nn-n-1)*eim*dcmplx(fr,0d0)
                  mdvi1=mm*dvi1(nn)
                  fc1=mdvi1*fc
                  fc2=dvi2(nn)*fc
                  tmp1=Tmat(1,l,1,ll)*mdvs1+Tmat(2,l,1,ll)*dvs2(n)
                  tmp2=Tmat(2,l,2,ll)*dvs2(n)+Tmat(1,l,2,ll)*mdvs1
                  tmp3=Tmat(1,l,1,ll)*dvs2(n)+Tmat(2,l,1,ll)*mdvs1
                  tmp4=Tmat(2,l,2,ll)*mdvs1+Tmat(1,l,2,ll)*dvs2(n)
                  S(1,1)=S(1,1)+tmp1*fc1+tmp2*fc2
                  S(1,2)=S(1,2)+tmp1*fc2+tmp2*fc1
                  S(2,1)=S(2,1)+tmp3*fc1+tmp4*fc2
                  S(2,2)=S(2,2)+tmp3*fc2+tmp4*fc1
 7             continue
 8          continue
 991        continue
 9       continue
 10   continue

      S(1,1)=S(1,1)
      S(1,2)=-Im*S(1,2)
      S(2,1)=Im*S(2,1)
      S(2,2)=S(2,2)

      
      end


c***********************************************************************
c     Calculation of the functions
c     dv1(n)=dvig(0,m,n,theta)/sin(theta)
c     and
c     dv2(n)=[d/d(theta)] dvig(0,m,n,theta)
c     1.le.n.le.nmax
c     0 <= theta <= 360
c
c     Modifications by M. Kahnert:
c     Make use of symmetry relations of the Wigner d-functions.
c
c***********************************************************************
      subroutine vigampl (theta, nmax, m, dv1, dv2)
      implicit none
      real*8 dv1(nmax), dv2(nmax),x,theta,pi

      integer n,nmax,m,i,i2
      real*8 dsi,a,qs,qs1,d1,d2,d3,der,dn,qn,qn1,qn2,qmm,qnm,qnm1

      pi=dacos(-1d0)
      x=dabs(dcos(theta))
      do 1 n=1,nmax
         dv1(n)=0d0
         dv2(n)=0d0
    1 continue

      if (dabs(1d0-x).le.1d-10) go to 100
      dsi=1d0/dsqrt(1d0-x*x)
      a=1d0
      qs=dsqrt(1d0-x*x)
      qs1=1d0/qs
      if (m.ne.0) go to 20
      d1=1d0
      d2=x
      do 5 n=1,nmax
         qn=dble(n)
         qn1=dble(n+1)
         qn2=dble(2*n+1)
         d3=(qn2*x*d2-qn*d1)/qn1
         der=qs1*(qn1*qn/qn2)*(-d1+d3)
         dv1(n)=d2*dsi
         dv2(n)=der
         d1=d2
         d2=d3
    5 continue
      goto 1000
   20 qmm=dble(m*m)
      do 25 i=1,m
         i2=i*2
         a=a*dsqrt(dble(i2-1)/dble(i2))*qs
   25 continue
      d1=0d0
      d2=a
      do 30 n=m,nmax
         qn=dble(n)
         qn2=dble(2*n+1)
         qn1=dble(n+1)
         qnm=dsqrt(qn*qn-qmm)
         qnm1=dsqrt(qn1*qn1-qmm)
         d3=(qn2*x*d2-qnm*d1)/qnm1
         der=qs1*(-qn1*qnm*d1+qn*qnm1*d3)/qn2
         dv1(n)=d2*dsi
         dv2(n)=der
         d1=d2
         d2=d3
   30 continue
      goto 1000
  100 if (m.ne.1) goto 1000
      do 110 n=1,nmax
         dn=dble(n*(n+1))
         dn=0.5d0*dsqrt(dn)
         if (x.lt.0d0) dn=dn*(-1)**(n+1)
         dv1(n)=dn
         if (x.lt.0d0) dn=-dn
         dv2(n)=dn
  110 continue
c-----------------------------------------------------------------------
c     For scattering angles p/2 < Theta <= pi, use the following symmetry
c     relations of the Wigner d-functions (see Varshalovich, Ch. 4.4):
c     d_{0,m}^n(pi-Theta)=(-1)**(n-m) * d_{0,m}^n(Theta),
c     which implies 
c        dv1(n)(pi-Theta)= (-1)**(n+m) * dv1(n)(Theta); and
c        dv2(n)(pi-Theta)=-(-1)**(n+m) * dv2(n)(Theta)
c-----------------------------------------------------------------------
 1000 if(theta.gt.(pi/2d0).and.theta.le.pi)then
         do n=1,nmax
            dv1(n)= (-1d0)**(n+m) * dv1(n)
            dv2(n)=-(-1d0)**(n+m) * dv2(n)
         enddo
c-----------------------------------------------------------------------
c     For scattering angles pi < Theta <= 3*pi/2, use the following symmetry
c     relations of the Wigner d-functions (see Varshalovich, Ch. 4.4):
c     d_{0,m}^n(Theta-pi)=(-1)**n * d_{0,m}^n(Theta),
c     which implies 
c        dv1(n)(2*pi-Theta)=-(-1)**n * dv1(n)(Theta); and
c        dv2(n)(2*pi-Theta)= (-1)**n * dv2(n)(Theta)
c-----------------------------------------------------------------------
      elseif(theta.gt.pi.and.theta.le.1.5d0*pi)then
         do n=1,nmax
            dv1(n)=-(-1d0)**n * dv1(n)
            dv2(n)= (-1d0)**n * dv2(n)
         enddo
c-----------------------------------------------------------------------
c     For scattering angles 3*pi/2 < Theta <= 2*pi, use the following symmetry
c     relations of the Wigner d-functions (see Varshalovich, Ch. 4.4):
c     d_{0,m}^n(2*pi-Theta)=(-1)**m * d_{0,m}^n(Theta),
c     which implies 
c        dv1(n)(2*pi-Theta)=-(-1)**m * dv1(n)(Theta); and
c        dv2(n)(2*pi-Theta)=-(-1)**m * dv2(n)(Theta)
c-----------------------------------------------------------------------
      elseif(theta.gt.1.5d0*pi.and.theta.le.2d0*pi)then
         do n=1,nmax
            dv1(n)=-(-1d0)**m * dv1(n)
            dv2(n)=-(-1d0)**m * dv2(n)
         enddo
      endif

      return
      end


   
c***********************************************************************
C  This subroutine computes the Mueller matrix from given amplitude 
c  scattering matrix for particles in fixed orientation.  
c     
c  Modification by M. Kahnert: 
C  Added computation of the polarised differential scattering 
c  cross sections .  
C
C  Arguments:
C    NAME         I/O  TYPE        DESCRIPTION
C    S             I   complex*16  amplitude scatt. matrix, dim. (2,2)
C    F             O   real*8      Mueller matrix, dimension (4,4)
C    dCsca         O   real*8      differential scattering cross section
C                                  multiplied by k**2, dimension (4)
C
c***********************************************************************
      subroutine fmatrix(S,F,dCsca)
      implicit none

      real*8 F(4,4),dCsca(4),m1,m2,m3,m4,s23,s41,s24,s31,s21,s34,
     &       d23,d41,d24,d31,d21,d34
      complex*16 S(2,2),tmp

      m1 = S(2,2)*conjg(S(2,2))
      m2 = S(1,1)*conjg(S(1,1))
      m3 = S(1,2)*conjg(S(1,2))
      m4 = S(2,1)*conjg(S(2,1))
      tmp= S(1,1)*conjg(S(1,2))
      s23= dreal(tmp)
      d23=-dimag(tmp)
      tmp= S(2,1)*conjg(S(2,2))
      s41= dreal(tmp)
      d41=-dimag(tmp)
      tmp= S(1,1)*conjg(S(2,1))
      s24= dreal(tmp)
      d24=-dimag(tmp)
      tmp= S(1,2)*conjg(S(2,2))
      s31= dreal(tmp)
      d31=-dimag(tmp)
      tmp= S(1,1)*conjg(S(2,2))
      s21= dreal(tmp)
      d21=-dimag(tmp)
      tmp= S(1,2)*conjg(S(2,1))
      s34= dreal(tmp)
      d34=-dimag(tmp)
      
      F(1,1)= 0.5d0*(m1+m2+m3+m4)
      F(1,2)= 0.5d0*(m2-m3+m4-m1)
      F(1,3)= s23+s41
      F(1,4)=-(d23+d41)
      F(2,1)= 0.5d0*(m2+m3-m4-m1)
      F(2,2)= 0.5d0*(m2-m3-m4+m1)
      F(2,3)= s23-s41
      F(2,4)=-d23+d41
      F(3,1)= s24+s31
      F(3,2)= s24-s31
      F(3,3)= s21+s34
      F(3,4)=-d21+d34
      F(4,1)= d24+d31
      F(4,2)= d24-d31
      F(4,3)= d21+d34
      F(4,4)= s21-s34

      dCsca(1)=m2 ! k**2 * [dCsca/dOmega]_{hh}
      dCsca(2)=m4 ! k**2 * [dCsca/dOmega]_{hv}
      dCsca(3)=m3 ! k**2 * [dCsca/dOmega]_{vh}
      dCsca(4)=m1 ! k**2 * [dCsca/dOmega]_{vv}
      end
