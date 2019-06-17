c***********************************************************************
c     This file contains the code GSPHERE for generating 
c     Gaussian random sphere geometries. The code has been developed
c     by Karri Muinonen and Timo Nousiainen, University of Helsinki.
c     A description of the theory is found in K. Muinonen et al.,
c     JQSRT 55, 577-601, 1996. Thanks to Karri Muinonen and 
c     Timo Nousiainen for helpful discussions about their code, and
c     for their kind permission to include it in the distribution of
c     Tsym.
c     
c     The code has been slightly adopted to Tsym. Most importantly,
c     the partial derivatives of the particle's surface parameterisation
c     have been implemented. (These partial derivatives are needed 
c     for evaluating the surface integrals over the particle surface 
c     in Waterman's T-matrix method).
c
c***********************************************************************      
      subroutine gsphere(gflg,cflg,sig,nuc,gami,lmin,lmax,nthe,
     &     nphi,theta,phis,nss,radius,drdth,drdphi)

c     GSPHERE generates sample Gaussian spheres. Version 2002-12-16.
c     
c     
c     Free software license 
c     
c     G-sphere is free software for the generation of sample Gaussian
c     spheres. It is available under the GNU General Public License that you
c     can find on the World Wide Web (http://www.gnu.org/licenses/gpl.txt)
c     and in the file G-sphere/GPL/gpl.txt.
c     
c     Contact addresses for G-sphere Authors:
c     
c     Karri Muinonen
c     Observatory, University of Helsinki
c     Kopernikuksentie 1, P.O. Box 14
c     FIN-00014 U. Helsinki
c     Finland
c     E-mail: Karri.Muinonen@helsinki.fi
c     
c     Timo Nousiainen
c     Geophysical Research, Finnish Meteorological Institute
c     Vuorikatu 15 A, P.O. Box 503
c     FIN-00101 Helsinki
c     Finland
c     E-mail: Timo.Nousiainen@helsinki.fi
c     
c     G-sphere, Copyright (C) 2002 by the G-sphere Authors Karri Muinonen
c     and Timo Nousiainen. 
c     
c     This program is free software; you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as published by
c     the Free Software Foundation; either version 2 of the License, or
c     (at your option) any later version.
c     
c     This program is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c     GNU General Public License for more details.
c     
c     You should have received a copy of the GNU General Public License
c     along with this program; if not, write to the Free Software
c     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

      implicit none
      integer nss,lmin,lmax,
     -     j0,j1,j2,gflg,cflg
      integer phimax,thmax,legmax
      parameter (legmax=512)
c      parameter (thmax=180,phimax=360,legmax=256)

      integer nthe,nphi,nnod
      double precision THETA(nthe,nphi),
     -     MUS(nthe),PHIS(nphi),RADIUS(nthe,nphi),
     -     DRDTH(nthe,nphi),DRDPHI(nthe,nphi)

      double precision ACF(0:legmax,0:legmax),BCF(0:legmax,0:legmax),
     -     SCFSTD(0:legmax,0:legmax),CSCF(0:legmax),
     -     EU(3),CEU(3),SEU(3),a,sig,beta,gami,elli,
     -     gam,nuc,ell,cs2d,cs4d,the,mu,grid,rmax,pi,rd
      character ca*2

      integer irnd
      double precision RNDU
      common irnd

c     Initializations:

      pi=4.0d0*atan(1.0d0)
      rd=pi/180.0d0

      irnd=1
      irnd=4*irnd+1
      a=RNDU(irnd)

c     Input parameters from option file:

      a=1.0d0

c     Input check:

      if (gflg.ne.1 .and. gflg.ne.2) stop
     -     'Trouble in gsphere: general or axisymmetric spheres.'
      if (cflg.ne.1 .and. cflg.ne.2 .and. cflg.ne.3) stop
     -     'Trouble in gsphere: correlation function unknown.'

      if (sig.le.0.0d0) stop
     -     'Trouble in gsphere: standard deviation .le. 0.'

      if (cflg.eq.2) then
         if (gami.le.0.0d0 .or. gami.gt.180.0d0) stop
     -        'Trouble in gsphere: input angle .le. 0. .or.  .gt. 180'
         if (lmin.gt.0 .or. lmax.lt.int(300.0d0/gami)) then
            write(6,*)
     &           'Warning in GSPHERE: correlation angle will differ '
            write(6,*)'from input value. Set minimum degree to 0 and '
            write(6,*)'maximum degree .gt. (300 deg)/(input value).'
         endif
      endif

      if (lmax.gt.legmax) stop
     -     'Trouble in gsphere: maximum degree .gt. legmax.'
      if (lmin.lt.0) stop
     -     'Trouble in gsphere: minimum degree .lt. 0.'
      if (lmin.gt.lmax) stop
     -     'Trouble in gsphere: minimum degree .lt. maximum degree.'
      if (cflg.eq.1 .and. lmin.lt.2) stop
     -     'Trouble in gsphere: minimum degree .lt.2.'

      if (nss.le.0) stop
     -     'Trouble in gsphere: sphere identification number .lt. 0.'

c     Miscellaneous:

      gami=gami*rd
      elli=2.0d0*sin(0.5d0*gami)
 
c     Initialization of the Gaussian random sphere:

      beta=sqrt(log(sig**2+1.0d0))
      if     (cflg.eq.1) then
         call CS1CF(CSCF,nuc,lmin,lmax,legmax)
      elseif (cflg.eq.2) then
         call CS2CF(CSCF,elli,lmin,lmax,legmax)
      elseif (cflg.eq.3) then
         call CS3CF(CSCF,lmin,lmax,legmax)
      endif

      do 40 j1=lmin,lmax
         if (CSCF(j1).lt.0.0d0) stop
     -        'Trouble in gsphere: negative Legendre coefficient.'
 40   continue

      call SGSCFSTD(SCFSTD,CSCF,beta,lmin,lmax,legmax)

c     Generate a sample Gaussian sphere with identification number
c     nss, then move to discretize and output:

      do 100 j0=1,nss
         if (gflg.eq.1) then
            call SGSCF(ACF,BCF,SCFSTD,lmin,lmax,legmax)
         else
            call SGSAXCF(ACF,CEU,SEU,SCFSTD,lmin,lmax,legmax)
         endif
 100  continue
      
c     Spherical-coordinate representation for general and axisymmetric shapes:

      call SPHDS(THETA,MUS,PHIS,nthe,nphi,nthe,nphi)

      if (gflg.eq.1) then
         call RGSSD(RADIUS,DRDTH,DRDPHI,MUS,PHIS,ACF,BCF,rmax,beta,
     -        nthe,nphi,lmin,lmax,nthe,nphi,legmax)
      else
         call RGSAXSD(RADIUS,DRDTH,DRDPHI,MUS,PHIS,ACF,CEU,SEU,rmax,
     -        beta,nthe,nphi,lmin,lmax,nthe,nphi,legmax)
      endif

      end


c***********************************************************************
c Discretization:
c
c SPHDS: spherical coordinates 
c***********************************************************************




       subroutine SPHDS(THETA,MU,PHI,nthe,nphi,thmax,phimax)

c SPHDS discretizes the spherical surface into a polar angle -azimuth angle
c grid. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer nthe,nphi,j1,j2,thmax,phimax
       double precision theta(thmax,phimax),MU(thmax),PHI(phimax),
     &      dthe,dphi,pi

c$$$       pi=4.0d0*atan(1.0d0)
c$$$       dthe=pi/nthe
c$$$       dphi=2.0d0*pi/nphi
c$$$
c$$$       do 20 j1=0,nthe
c$$$        MU(j1)=cos(j1*dthe)
c$$$        do 10 j2=0,nphi
c$$$         PHI(j2)=(j2+0.5d0)*dphi
c$$$10      continue
c$$$20     continue

       do j1=1,nthe
          mu(j1)=dcos(theta(j1,1))
       enddo

       end



c***********************************************************************
c Gaussian sphere generator
c
c RGSSD:   discrete spherical-coordinate representation
c RGS:     radial distance
c SGS:     logarithm of radial distance
c SGSCF:   spherical harmonics coefficient generation
c***********************************************************************




        subroutine RGSSD(RADIUS,DRDTH,DRDPHI,MU,PHI,ACF,BCF,rmax,beta,
     -                   nthe,nphi,lmin,lmax,thmax,phimax,legmax)

c Discrete spherical-coordinate representation for a sample G-sphere.
c Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer nthe,nphi,lmin,lmax,j1,j2,thmax,phimax,legmax
       double precision MU(thmax),PHI(phimax),
     -  ACF(0:legmax,0:legmax),BCF(0:legmax,0:legmax),rmax,beta,
     -  r,rth,rphi,radius(thmax,phimax),drdth(thmax,phimax),
     -  drdphi(thmax,phimax),nu

       rmax=0.0d0
       do 20 j1=1,nthe
        nu=sqrt(1.0d0-MU(j1)**2)
        do 10 j2=1,nphi
         call RGS(ACF,BCF,MU(j1),PHI(j2),beta,lmin,lmax,legmax,
     &          r,rth,rphi)
         RADIUS(j1,j2)=r
         DRDTH(j1,j2) =rth
         DRDPHI(j1,j2)=rphi
         if (r.gt.rmax) rmax=r
10      continue
20     continue
       end



       subroutine RGS(ACF,BCF,mu,phi,beta,lmin,lmax,
     &     legmax,r,rth,rphi)

c Radial distance in a given direction for a sample G-sphere.
c Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer lmin,lmax,legmax
       double precision ACF(0:legmax,0:legmax),BCF(0:legmax,0:legmax),
     -  mu,phi,beta ,s,sth,sphi,r,rth,rphi

       call SGS(ACF,BCF,mu,phi,lmin,lmax,legmax,s,sth,sphi)

       r=exp(s-0.5d0*beta**2)
       rth=sth*r
       rphi=sphi*r

       end



       subroutine SGS(ACF,BCF,mu,phi,lmin,lmax,legmax,s,sth,sph)

c Logarithmic radial distance s in a given direction for a sample G-sphere.
c Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen
c --
c Added ds/d theta and ds/d phi; Michael Kahnert, 2011

       implicit none
       integer l,m,lmin,lmax,legmax
       double precision ACF(0:legmax,0:legmax),BCF(0:legmax,0:legmax),
     -  LEGP(0:legmax,0:legmax),DLEGP(0:legmax,0:legmax),
     -  CPHI(legmax),SPHI(legmax),mu,phi,s,sth,sph

       if (lmax.eq.0) then
        S=ACF(0,0)
        return
       endif

c Precomputation of sines, cosines, associated Legendre functions,
c and their derivatives

       call LEGA(LEGP,mu,lmax,0,legmax)
       do 10 m=1,lmax
        call LEGA(LEGP,mu,lmax,m,legmax)
        CPHI(m)=cos(m*phi)
        SPHI(m)=sin(m*phi)
10     continue
       LEGP(0,0)=1.0d0

c Precomputed the derivative of the associated Legendre functions wrt the
c polar coordinate:
       call dlega(legp,dlegp,mu,lmax,legmax)

c Compute s:
       S=0.0d0
       do 20 l=lmin,lmax
        S=S+LEGP(l,0)*ACF(l,0)
20     continue
       do 40 m=1,lmax
        do 30 l=max(m,lmin),lmax
         S=S+LEGP(l,m)*(ACF(l,m)*CPHI(m)+BCF(l,m)*SPHI(m))
30      continue
40     continue

c Compute ds/d theta:
       sth=0d0
       do l=lmin,lmax
          sth=sth+DLEGP(l,0)*ACF(l,0)
       enddo
       do m=1,lmax
          do l=max(m,lmin),lmax
             sth=sth+DLEGP(l,m)*(ACF(l,m)*CPHI(m)+BCF(l,m)*SPHI(m))
          enddo
       enddo
       
c Compute ds/d phi:
       sph=0d0
       do m=1,lmax
          do l=max(m,lmin),lmax
             sph=sph+LEGP(l,m)*dble(m)*
     &            (-ACF(l,m)*SPHI(m)+BCF(l,m)*CPHI(m))
          enddo
       enddo
       
       end



       subroutine SGSCF(ACF,BCF,SCFSTD,lmin,lmax,legmax)

c Generates the sample spherical harmonics coefficients for the 
c logarithmic radial distance of the G-sphere. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,lmin,lmax,m,legmax
       double precision ACF(0:legmax,0:legmax),BCF(0:legmax,0:legmax),
     -  SCFSTD(0:legmax,0:legmax),rn

       do 10 l=lmin,lmax
        call RNDG(rn)
        ACF(l,0)=rn*SCFSTD(l,0)
        BCF(l,0)=0.0d0
10     continue
       do 30 m=1,lmax
        do 20 l=max(m,lmin),lmax
         call RNDG(rn)
         ACF(l,m)=rn*SCFSTD(l,m)
         call RNDG(rn)
         BCF(l,m)=rn*SCFSTD(l,m)
20      continue
30     continue
       end


c***********************************************************************
c Gaussian axisymmetric sphere generator
c
c RGSAXSD: discrete spherical-coordinate representation
c RGSAX:   radial distance
c SGSAX:   logarithm of radial distance
c SGSAXCF: spherical harmonics coefficient and random orientation generator
c***********************************************************************

        subroutine RGSAXSD(RADIUS,DRDTH,DRDPHI,MU,PHI,ACF,CEU,SEU,rmax,
     -                     beta,nthe,nphi,lmin,lmax,thmax,phimax,legmax)

c Discrete spherical-coordinate representation for a sample axisymmetric
c G-sphere. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer nthe,nphi,lmin,lmax,j1,j2,thmax,phimax,legmax
       double precision radius(thmax,phimax),drdth(thmax,phimax),
     -  drdphi(thmax,phimax),MU(thmax),PHI(phimax),
     -  ACF(0:legmax,0:legmax),CEU(3),SEU(3),rmax,beta,
     -  r,rth,rphi,nu,phi0

       phi0=0d0
       rmax=0.0d0
       do 20 j1=1,nthe
        nu=sqrt(1.0d0-MU(j1)**2)
        call RGSAX(ACF,CEU,SEU,MU(j1),phi0,beta,lmin,lmax,legmax,
     &       r,rth,rphi)
        do 10 j2=1,nphi
           radius(j1,j2)=r
           drdth(j1,j2)=rth
           drdphi(j1,j2)=rphi
           if (r.gt.rmax) rmax=r
 10     continue
20     continue
       end



       subroutine RGSAX(ACF,CEU,SEU,mu,phi,beta,lmin,lmax,legmax,
     &     r,rth,rphi)

c Radial distance in a given direction for a sample axisymmetric G-sphere.
c Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer lmin,lmax,legmax
       double precision ACF(0:legmax,0:legmax),CEU(3),SEU(3),
     -  mu,phi,beta,r,rth,rphi,s,sth


       call SGSAX(ACF,CEU,SEU,mu,phi,lmin,lmax,legmax,s,sth)
       r=exp(s-0.5d0*beta**2)
       rth=r*sth
       rphi=0d0
       end



       subroutine SGSAX(ACF,CEU,SEU,mu,phi,lmin,lmax,
     &     legmax,s,sth)

c Logarithmic radial distance in a given direction for a sample 
c axisymmetric G-sphere. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,m,lmin,lmax,legmax
       double precision ACF(0:legmax,0:legmax),LEGP(0:legmax,0:legmax),
     -  DLEGP(0:legmax,0:legmax),CEU(3),SEU(3),X(3),mu,nu,phi,cphi,sphi,
     -  r0,mu0,s,sth

       if (lmax.eq.0) then
        s=ACF(0,0)
        sth=0d0
        return
       endif

c Euler rotation:

       nu=sqrt(1.0d0-mu**2)
       cphi=cos(phi)
       sphi=sin(phi)
       X(1)=nu*cphi
       X(2)=nu*sphi
       X(3)=mu
       call VROTEU(X,CEU,SEU)
       
       r0=sqrt(X(1)**2+X(2)**2+X(3)**2)
       if (r0.lt.1.0d-12) 
     -  stop 'Trouble in SGSAX: radial distance too small.'
       if(abs(r0-1d0).gt.1d-12)
     &     write(6,*)'Warning SGSAX: Norm not preserved after rotation'
       mu0=X(3)/r0

c Precomputation of Legendre polynomials and their derivatives:

       call LEGA(LEGP,mu0,lmax,0,legmax)
       LEGP(0,0)=1.0d0
       call dlega(legp,dlegp,mu0,lmax,legmax)


c Sum up:

       s=0.0d0
       sth=0.0d0
       do 10 l=lmin,lmax
          s=s+LEGP(l,0)*ACF(l,0)
          sth=sth+DLEGP(l,0)*ACF(l,0)
10     continue
       end



       subroutine SGSAXCF(ACF,CEU,SEU,SCFSTD,lmin,lmax,legmax)

c Generates the sample spherical harmonics coefficients and Euler rotation
c angle sines and cosines for the logarithmic radial distance of 
c the axisymmetric G-sphere. Note that fixed orientation are currently used.
c To generate random Euler orientations, use the code that is presently 
c commented out.
c Version 2002-12-16, changed by Michael Kahnert 2012-03-07.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer irnd,l,lmin,lmax,legmax
       double precision ACF(0:legmax,0:legmax),CEU(3),SEU(3),
     -  SCFSTD(0:legmax,0:legmax),EU(3),RNDU,rn,pi
       parameter (pi=3.1415926535898d0)
       common irnd

       do 10 l=lmin,lmax
        call RNDG(rn)
        ACF(l,0)=rn*SCFSTD(l,0)*sqrt(dble(2*l+1))
10     continue
c       EU(1)=2.0d0*pi*RNDU(irnd)
c       EU(3)=2.0d0*pi*RNDU(irnd)
c       CEU(1)=cos(EU(1))
c       CEU(2)=1.0d0-2.0d0*RNDU(irnd)
c       CEU(3)=cos(EU(3))
c       SEU(1)=sin(EU(1))
c       SEU(2)=sqrt(1.0d0-CEU(2)**2)
c       SEU(3)=sin(EU(3))

       EU(1)=0.0d0
       EU(3)=0.0d0
       CEU(1)=1.0d0
       CEU(2)=1.0d0
       CEU(3)=1.0d0
       SEU(1)=0.0d0
       SEU(2)=0.0d0
       SEU(3)=0.0d0
       end



c***********************************************************************
c Correlation functions and related coefficients:
c
c SGSCFSTD: standard deviations for the spherical harmonics coeffiecients
c CSELL   : correlation length and angle, derivatives
c CSLEGP  : Legendre expansion of the correlation function
c CS1CF   : power law correlation function
c CS2CF   : modified Gaussian correlation function
c CS3CF   : input correlation function coefficients from file
c***********************************************************************


       subroutine SGSCFSTD(SCFSTD,CSCF,beta,lmin,lmax,legmax)

c Generates the standard deviations for the spherical harmonics 
c coefficients of the logarithmic radial distance. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,lmin,lmax,m,legmax
       double precision SCFSTD(0:legmax,0:legmax),CSCF(0:legmax),FACTI,
     &      beta

       if (lmin.eq.0) then
        SCFSTD(0,0)=beta*sqrt(CSCF(0))
        do 20 l=1,lmax
         SCFSTD(l,0)=beta*sqrt(CSCF(l))
         do 10 m=1,l
          SCFSTD(l,m)=SCFSTD(l,0)*
     -     sqrt(2.0d0*FACTI(l-m)/FACTI(l+m))
10       continue
20      continue
       else
        do 40 l=lmin,lmax
         SCFSTD(l,0)=beta*sqrt(CSCF(l))
         do 30 m=1,l
          SCFSTD(l,m)=SCFSTD(l,0)*
     -     sqrt(2.0d0*FACTI(l-m)/FACTI(l+m))
30       continue
40      continue
       endif
       end



       subroutine CSELL(CSCF,gam,ell,cs2d,cs4d,lmin,lmax,legmax)

c Computes the correlation angle and length, and second and fourth 
c derivatives for a correlation function expressed in Legendre series with
c normalized coefficients. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,lmin,lmax,legmax
       double precision CSCF(0:legmax),gam,ell,cs2d,cs4d

       cs2d=0.0d0
       cs4d=0.0d0
       do 10 l=lmin,lmax
        cs2d=cs2d-CSCF(l)*l*(l+1.0d0)/2.0d0
        cs4d=cs4d+CSCF(l)*l*(l+1.0d0)*(3.0d0*l**2+3.0d0*l-2.0d0)/8.0d0
10     continue
       if (cs2d.eq.0.0d0) then
        ell=2.0d0
       else
        ell=1.0d0/sqrt(-cs2d)
       endif
       gam=2.0d0*asin(0.5d0*ell)
       end



       double precision function CSLEGP(CSCF,xi,lmin,lmax,legmax)

c Computes the Legendre series expansion for the 
c correlation function. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,lmin,lmax,legmax
       double precision LEGP(0:legmax,0:legmax),CSCF(0:legmax),xi

       call LEGA(LEGP,xi,lmax,0,legmax)
       LEGP(0,0)=1.0d0

       CSLEGP=0.0d0
       do 10 l=lmin,lmax
        CSLEGP=CSLEGP+LEGP(l,0)*CSCF(l)
10     continue
       end



       subroutine CS1CF(CSCF,nu,lmin,lmax,legmax)

c Returns the Legendre coefficients for the correlation 
c function with power-law Legendre coefficients. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,lmin,lmax,legmax
       double precision CSCF(0:legmax),nu,norm

       do 10 l=0,lmin-1
        CSCF(l)=0.0d0
10     continue

       norm=0.0d0
       do 20 l=lmin,lmax
        CSCF(l)=1.0d0/l**nu
        norm=norm+CSCF(l)
20     continue

       do 30 l=lmin,lmax
        CSCF(l)=CSCF(l)/norm
30     continue
       end



       subroutine CS2CF(CSCF,ell,lmin,lmax,legmax)

c Returns the Legendre coefficients for the modified Gaussian correlation 
c function. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,lmin,lmax,legmax
       double precision BESISE(0:legmax),CSCF(0:legmax),ell,z,norm

       z=1.0d0/ell**2
       call BESMS(BESISE,z,lmax,legmax)

       do 10 l=0,lmin-1
        CSCF(l)=0.0d0
10     continue

       norm=0.0d0
       do 20 l=lmin,lmax
        CSCF(l)=(2*l+1)*BESISE(l)
        norm=norm+CSCF(l)
20     continue

       do 30 l=lmin,lmax
        CSCF(l)=CSCF(l)/norm
30     continue
       end



       subroutine CS3CF(CSCF,lmin,lmax,legmax)

c Inputs Legendre coefficients for the correlation 
c function from a file (and normalizes). Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer l,ll,lmin,lmax,legmax
       double precision CSCF(0:legmax),norm

       open(unit=1, file='cscf.dat', status='old')

       do 10 l=0,lmin-1
        CSCF(l)=0.0d0
10     continue

       norm=0.0d0
       do 20 l=lmin,lmax
        read (1,*) ll,CSCF(l) 
        if (l.ne.ll) stop
     -   'Trouble in CS3CF: degree inconsistency.'
        if (CSCF(ll).lt.0.0d0) stop
     -   'Trouble in CS3CF: negative Legendre coefficient.'
        norm=norm+CSCF(l)
20     continue
       close(unit=1)

       if (norm.eq.0.0d0) stop
     -  'Trouble in CS3CF: no nonzero Legendre coefficients.'

       do 30 l=lmin,lmax
        CSCF(l)=CSCF(l)/norm
30     continue
       end



c***********************************************************************
c Vector rotations and product:
c
c VROTEU: vector rotation using Euler angles
c VROTX:  vector rotation about the x-axis
c VROTY:  vector rotation about the y-axis
c VROTZ:  vector rotation about the z-axis




       subroutine VROTEU(X,CA,SA)

c Vector rotation using Euler angles. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       double precision X(3),CA(3),SA(3)

       call VROTZ(X,CA(1),SA(1))
       call VROTY(X,CA(2),SA(2))
       call VROTZ(X,CA(3),SA(3))
       end



       subroutine VROTX(X,c,s)

c Vector rotation about the x-axis. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       double precision X(3),c,s,q

       q   = c*X(2)+s*X(3)
       X(3)=-s*X(2)+c*X(3)
       X(2)=q
       end



       subroutine VROTY(X,c,s)

c Vector rotation about the y-axis. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       double precision X(3),c,s,q

       q   = c*X(3)+s*X(1)
       X(1)=-s*X(3)+c*X(1)
       X(3)=q
       end



       subroutine VROTZ(X,c,s)

c Vector rotation about the z-axis. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       double precision X(3),c,s,q

       q   = c*X(1)+s*X(2)
       X(2)=-s*X(1)+c*X(2)
       X(1)=q
       end



***********************************************************************
c Special functions:
c
c BESMS:  modified spherical Bessel functions multiplied by exponential
c LEGA:   associated Legendre functions
c SPHGEN: generalized spherical functions 
c GAMLN:  logarithmic Gamma function 
c FACTI:  factorial function
c***********************************************************************


       subroutine BESMS(BESISE,x,n,legmax)

c Generates modified spherical Bessel functions multiplied by an
c exponential: i_0(x)*exp(-x),...,i_n(x)*exp(-x). Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer j1,j2,n,nini,legmax,ifac
       double precision BESISE(0:legmax),x,b,b0,b1,renorm,fac
      
c     empirical scaling to accomodate angular correlations down to 1.8 degrees:
       if(x.lt.210)then
          ifac=1
       elseif(x.lt.910)then
          ifac=2
       else
          ifac=3
       endif
       fac=dble(ifac)


c Orders n=0 and n=1:

       if (n.le.1) then
          BESISE(0)=exp(-x)*sinh(x)/x
          BESISE(1)=exp(-x)*(-sinh(x)/x**2)+cosh(x)/x
          return
       endif

c     Downward recurrence:

       nini=max(n+4,int(1.5d0*x))
       b1=0.0d0
       b0=exp(-x/fac)*2.0d0*x
       do j1=nini,n,-1
          b=(2*j1+1)*b0/x+b1
          b1=b0
          b0=b
       enddo
       BESISE(n)=b1
       BESISE(n-1)=b0
       do j1=n,2,-1
          BESISE(j1-2)=(2*j1-1)*BESISE(j1-1)/x+BESISE(j1)
       enddo

c     Renormalization:

c     renorm=exp(-x)*(sinh(x)/x)/BESISE(0)
       renorm=(1d0-exp(-2d0*x))/(2d0*x)/BESISE(0)
       do j1=0,n
          BESISE(j1)=renorm*BESISE(j1)
          do j2=1,ifac-1
             BESISE(j1)=BESISE(j1)*exp(-x/fac)
          enddo
       enddo


       end



       subroutine LEGA(LEGP,x,lmax,m,legmax)

c Computes associated Legendre functions from degree l=m
c up to l=lmax. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer lmax,l,m,legmax
       double precision LEGP(0:legmax,0:legmax),FACTI,x
       double complex GSP(0:legmax,0:legmax,-2:2),i

       i=dcmplx(0.0d0,1.0d0)

c Check degree, orders, and argument:

       if (lmax.lt.0)
     -  stop 'Trouble in LEGA: degree negative.'

       if (m.gt.lmax .or. m.lt.0)
     -  stop 'Trouble in LEGA: order out of range.'

       if (abs(x).gt.1.0d0) then
          write(6,*)'Trouble in LEGA: argument out of range, x=',x
          stop
       endif

c Compute associated Legendre functions with the help of
c the generalized spherical functions:
   
       call SPHGEN(GSP,x,lmax,m,0,legmax)

       do 10 l=m,lmax
        LEGP(l,m)=dreal(i**m*
     *   sqrt(FACTI(l+m)/FACTI(l-m))*GSP(l,m,0))
10     continue
       end


      subroutine dlega(legp,dlegp,x,lmax,legmax)
c     
c     Compute derivative dP_n^m/d(theta) of the associated Legendre functions.

c     The following method has been commented out, because it is inaccurate for small theta:
c     Recurrence relations:
c        d P_0^0/d\theta = 0
c        d P_n^m/d\theta = x/\sqrt(1-x^2) n P_n^m + \sqrt(1-x^2) (n+m) P_{n-1}^m
c     (see Abramowitz and Stegun, Eq. 8.5.4); and
c        d P_n^m/d\theta = -(n+m)(n-m+1) P_n^(m-1) - m x/sqrt(1-x^2) P_n^m,
c     where x=\cos\theta
c     The last recurrence relation can be derived from the following relation of 
c     the Wigner D-functions:
c        \partial/\partial\theta D_{0,m}^n(0,\theta,\phi) 
c        = 1/2 \sqrt((n+m)(n-m+1))\exp(-i\phi) D_{m,m-1}^n(0,\theta,\phi)
c        - 1/2 \sqrt((n-m)(n+m+1))\exp(+i\phi) D_{m,m+1}^n(0,\theta,\phi)
c     (see Varshalovich et al., Quantum Theory of Angular Momentum, World Scientific 1988,
c     p 94, Eq. (3)), together with the relation with the Wigner d-functions
c        D_{0,m}^n(0,\theta,\phi) = d_{0,m}^n(\theta) \exp(-im\phi), 
c     and the relation to the associated Legendre functions
c        P_n^m(\theta) = (-1)^m \sqrt((n+m)!/(n-m)!) d_{0,m}^n(\theta).
c
c     Instead, the computation uses the routines for computing the derivatives of
c     the Wigner d-functions, and converts the results to associated Legendre functions:
c       d P_n^m(theta)/d theta = sqrt((n+m)!/(n-m)!) * d d_{0,m}^n(theta)/d theta
c
c     M. Kahnert, 2011-01-24 
c     corrected:  2011-11-12 --- a rainy Saturday :-( 

      implicit none

      integer lmax,l,m,legmax
      double precision legp(0:legmax,0:legmax),
     &     dlegp(0:legmax,0:legmax),x,s,ptmp(legmax),dptmp(legmax),
     &     FACTI

c      s=sqrt(1d0-x*x)
c      if(s.le.0d0)then
c         write(6,*)'FATAL ERROR in dlega: sin(theta) <= 0 ',s,x
c         stop
c      endif
c
c      dlegp(0,0)=0d0
c      do l=1,lmax
c         dlegp(l,0)=x/s*dble(l)*legp(l,0)+s*dble(l)*legp(l-1,0)
c      enddo
c      do m=1,lmax
c         do l=m,lmax
c            dlegp(l,m)=-dble((l+m)*(l-m+1))*legp(l,m-1) - dble(m)*x/s*
c     &           legp(l,m)
c         enddo
c      enddo

      do l=0,lmax
         do m=0,lmax
            dlegp(l,m)=0d0
         enddo
      enddo

      dlegp(0,0)=0d0
      call vig( x, legmax, 0, Ptmp, dPtmp)
      do l=1,lmax
         dlegp(l,0) = dPtmp(l)
      end do
      do  m = 1, lmax
         call vig( x, legmax, m, Ptmp, dPtmp)
         do  l = m, lmax
            dlegp(l,m) = sqrt(FACTI(l+m)/FACTI(l-m))*dPtmp(l)
         enddo
      enddo

      end



      subroutine SPHGEN(GSP,x,lmax,m1,m2,legmax)

c     Computes generalized spherical functions from degree max(abs(m1),abs(m2))
c     up to lmax. Version 2002-12-16.
c     
c     Copyright (C) 2002 Karri Muinonen

      implicit none
      integer lmax,l,m0,m1,m2,m12,p12,legmax
      double precision FACTI,x
      double complex GSP(0:legmax,0:legmax,-2:2),i

      i=dcmplx(0.0d0,1.0d0)

c     Check degree, orders, and argument:

      if (lmax.lt.0)
     -     stop 'Trouble in SPHGEN: degree negative.'

      if (abs(m1).gt.lmax .or. abs(m2).gt.min(2,lmax) .or. m1.lt.0)
     -     stop 'Trouble in SPHGEN: order out of range.'

      if (abs(x).gt.1.0d0)
     -     stop 'Trouble in SPHGEN: argument out of range.'

c     Compute generalized spherical functions:
      
      m0=max(abs(m1),abs(m2))
      m12=abs(m1-m2)
      p12=abs(m1+m2)

      if (m0.gt.0) then

         if (m12.ne.0 .and. p12.ne.0) then
            GSP(m0,m1,m2)=(-i)**m12/2.0d0**m0*
     -           sqrt(FACTI(2*m0)/(FACTI(m12)*FACTI(p12))*
     -           (1.0d0-x)**m12*(1.0d0+x)**p12)
         elseif (m12.eq.0) then
            GSP(m0,m1,m2)=1.0d0/2.0d0**m0*
     -           sqrt(FACTI(2*m0)/FACTI(p12)*(1.0d0+x)**p12)
         else
            GSP(m0,m1,m2)=(-i)**m12/2.0d0**m0*
     -           sqrt(FACTI(2*m0)/FACTI(m12)*(1.0d0-x)**m12)
         endif

         if (m0.eq.lmax) return

         GSP(m0+1,m1,m2)=(2*m0+1)*(m0*(m0+1)*x-m1*m2)*GSP(m0,m1,m2)/
     -      (m0*sqrt(dble((m0+1)**2-m2**2))*sqrt(dble((m0+1)**2-m1**2)))

         if (m0+1.eq.lmax) return

         do 10 l=m0+1,lmax-1
            GSP(l+1,m1,m2)=((2*l+1)*(l*(l+1)*x-m1*m2)*GSP(l,m1,m2)
     -     -(l+1)*sqrt(dble((l**2-m1**2)*(l**2-m2**2)))*GSP(l-1,m1,m2))/
     -           (l*sqrt(dble(((l+1)**2-m1**2)*((l+1)**2-m2**2))))
 10      continue

      else

         GSP(0,0,0)=1.0d0
         if (lmax.eq.0) return
         GSP(1,0,0)=x
         if (lmax.eq.1) return

         do 20 l=m0+1,lmax-1
            GSP(l+1,0,0)=((2*l+1)*x*GSP(l,0,0)-l*GSP(l-1,0,0))/(l+1)
 20      continue

      endif
      end



       function GAMLN(x)

c Function GAMLN computes the natural logarithm of the
c Gamma function with positive real argument. The algorithm
c is based on the approximation by Lanczos (1964: SIAM
c Journal on Numerical Analysis, ser. B, vol 1., pp. 86-96)
c Version 2002-12-16.
c
c Copyright (C) 2002 Timo Nousiainen

       implicit none
       double precision gamln,stp,pi,x,xx,temp,arg1,arg2
       parameter (pi=3.1415926535898d0)

       if (x.le.0.0d0) then
        write (6,*) 'Trouble in GAMLN: x must be positive! Exiting.'
        stop
       endif

       xx=x-1.0d0

       stp=sqrt(2.0d0*pi)
       arg1=xx+0.5d0
       arg2=xx+5.5d0
 
       temp=1.000000000190015d0+76.18009172947146d0/(xx+1.0d0)-
     -      86.50532032941677d0/(xx+2.0d0)+
     -      24.01409824083091d0/(xx+3.0d0)-
     -      1.231739572450155d0/(xx+4.0d0)+
     -      0.1208650973866179d-2/(xx+5.0d0)-
     -      0.5395239384953d-5/(xx+6.0d0)

       temp=dlog(temp*stp)+(arg1)*dlog(arg2)-(arg2)
       GAMLN=temp
       end



       function FACTI(x)

c Function FACTI returns a factorial of an integer argument x
c Values are precomputed for small values for speed. Although
c a factorial of an integer argument is also an integer, it is
c handled as real number with double precision to handle
c large values. Version 2002-12-16.
c
c Copyright (C) 2002 Timo Nousiainen

       implicit none
       integer x 
       double precision FACTI,GAMLN,xx

       if (x.lt.0) then
        write (6,*) 'Trouble in FACTI: x must be non-negative! Exiting.'
        stop
       endif

       if (x.le.14) then
        if (x.eq.0)  FACTI=1.0d0
        if (x.eq.1)  FACTI=1.0d0
        if (x.eq.2)  FACTI=2.0d0
        if (x.eq.3)  FACTI=6.0d0
        if (x.eq.4)  FACTI=24.0d0
        if (x.eq.5)  FACTI=120.0d0
        if (x.eq.6)  FACTI=720.0d0
        if (x.eq.7)  FACTI=5040.0d0
        if (x.eq.8)  FACTI=40320.0d0
        if (x.eq.9)  FACTI=362880.0d0
        if (x.eq.10) FACTI=3628800.0d0
        if (x.eq.11) FACTI=39916800.0d0
        if (x.eq.12) FACTI=479001600.0d0
        if (x.eq.13) FACTI=6227020800.0d0
        if (x.eq.14) FACTI=87178291200.0d0
       else
        xx=x*1.0d0
        FACTI=exp(GAMLN(xx+1.0d0))
       endif
       end



c***********************************************************************
c Random deviates:
c
c RNDU: uniform distribution within [0, 1)
c RNDG: Gaussian distribution with zero mean and unit standard deviation
c***********************************************************************


       double precision function RNDU(irnd)

c Quick and dirty uniform distribution [0,1) based on the TOMS
c library routine SUNIF. Initialize by calling RNDU(irnd) with 
c x = 4*n+1. Calls with negative argument return a random number.
c Version 2002-12-16.
c
c Copyright (C) 2002 Timo Nousiainen

       implicit none
       integer irnd
       double precision r,factor,two28
       save r
       data factor /41475557.0d0/, two28 /268435456.d0/

       if (irnd.ge.0) goto 10

       r=dmod(r*factor,1.0d0)
       RNDU=dble(r)
       return

c Initialization

10     r=dble(float(irnd))/two28
       r=dmod(r*factor,1.0d0)
       RNDU=dble(r)
       irnd=-1
       return
       end



       subroutine RNDG(r1)

c Returns a normally distributed random deviate with zero mean and 
c unit variance. Version 2002-12-16.
c
c Copyright (C) 2002 Karri Muinonen

       implicit none
       integer flg,irnd,xrandom
       double precision RNDU,q1,q2,r1,r2
       save flg,r2
       data flg/0/
       common irnd

       if (flg.eq.1) then
        r1=r2
        flg=0
        return
       endif

       flg=0
10     r1=2.0d0*RNDU(irnd)-1.0d0
       r2=2.0d0*RNDU(irnd)-1.0d0
       q1=r1**2+r2**2
       if (q1.ge.1.0d0 .or. q1.le.0.0d0) goto 10

       q2=sqrt(-2.0d0*log(q1)/q1)
       r1=r1*q2
       r2=r2*q2
       flg=1
       end


