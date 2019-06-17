c***********************************************************************
c The subroutines in this file were taken from the T-matrix code written 
c by Dan Mackowski for clusters of spheres (Mackowski and Mishchenko, 
c JOSA A 13, 2266-2278, 1996). Thanks to Dan Mackowski and to 
c Michal Mishchenko for some helpful discussions on how to use 
c their code, and for their kind permission to include these
c subroutines in the distribution of Tsym!
c***********************************************************************



c***********************************************************************
c Subroutine for calculating the orientationally-averaged Mueller matrix
c (also known as the Stokes scattering matrix), and the extintion and 
c total scattering cross sections from the T-matrix of an 
c arbritrarily-shaped particle.
c
c Some minor modifications were implemented for exploiting
c Cn-symmetry in the orientation-averaging (M. Kahnert)
c
c INPUT PARAMETERS:
c
c Tmat   : T-matrix, where T_{n,m,n',m'}^{i,j} = 
c                         Tmat(i,n*(n+1)+m,j,n'*(n'+1)+m'),
c                         n,n' = 1,2,...,nmax
c                         m    = -n,...,n
c                         m'   = -n',...,n'
c Pgroup : symmetry group
c Nsym   : index of main rotational symmetry
c nmax   : dimension of the T-matrix (see above)
c Forder : max. expansion order of the Stokes scattering matrix
c xv     : size parameter of the particle
c Rv     : volume-equivalent radius
c fout   : character string containing the name of the output file
c          of the Mueller matrix
c cout   : ditto for cross sections, SSA, asymmetry parameter, and
c          backscattering depolarisation ratios
c eout   : ditto for expansion coefficients of Mueller matrix elements
c nt     : number of scattering angles theta at which the scattering
c          matrix will be printed in the output file
c tmin   : min. angle theta
c tmax   : max. angle theta
c
c OUTPUT: printed to file
c      
c***********************************************************************
      subroutine ranprops(Tmat,Pgroup,Nsym,lambda,nmax,lmax2,
     &     Forder,nfd,nbc,nmax2,nbt1,xv,Rv,pi,fout,cout,eout,
     &     nt,tmin,tmax)
      implicit real*8(a-h,o-z)

      integer Nsym,nmax,nfd,lmax2,Forder,nbc,nmax2,nbt1,nt
      real*8 dc(-2:2,0:nfd),pi,lambda,tmin,tmax
      complex*16 aw(0:2,-1:1,0:nmax2),bw(0:2,-1:1,0:nmax2),
     1           cw(0:nmax2),dw(0:nmax2)
      complex*16 ci,Tmat(2,lmax2,2,lmax2)
      real*8 bcof(0:nbc,0:nbc),fnr(0:2*nbc),
     &       a1(nmax2+1), a2(nmax2+1),a3(nmax2+1),a4(nmax2+1),
     &       b1(nmax2+1),b2(nmax2+1)
      real*8 ct,Rv
      real*8 cext,csca,cabs,cbak,w,g,qext,qsca,qabs,qbak,deltaL,deltaC
      character fout*10,cout*10,eout*10
      character*3 Pgroup
      data ci/(0.d0,1.d0)/
      
c-----------------------------------------------------------------------
c calculation of constants
c-----------------------------------------------------------------------
      do n=1,2*nbc
         fnr(n)=dsqrt(dble(n))
      enddo
      bcof(0,0)=1.d0
      do n=0,nbc-1
         do l=n+1,nbc
            bcof(n,l)=fnr(n+l)*bcof(n,l-1)/fnr(l)
            bcof(l,n)=bcof(n,l)
         enddo
         bcof(n+1,n+1)=fnr(n+n+2)*fnr(n+n+1)*bcof(n,n)/fnr(n+1)/fnr(n+1)
      enddo
      
      nblkt=nmax*(nmax+2)
      qet=0.
      qst=0.
      do n=1,nblkt
         do ip=1,2
            qet=qet-Tmat(ip,n,ip,n)
            do l=1,nblkt
               do iq=1,2
                  qst=qst+Tmat(ip,n,iq,l)*conjg(Tmat(ip,n,iq,l))
               enddo
            enddo
         enddo
      enddo
      qet = qet*lambda**2/2.d0/pi
      qst = qst*lambda**2/2.d0/pi
      qat=qet-qst

      if(Forder.eq.-1) then
         Forder=2*nmax
      endif
      call tmsm(Tmat,Pgroup,Nsym,nmax,Forder,lmax2,nmax2,nfd,
     &          nbc,bcof,fnr,aw,bw,cw,dw)
c-----------------------------------------------------------------------
c The following lines have been changed to adapt the definition
c of the scattering matrix in the code to the convention
c used in the codes by Mishchenko, Laitinen et al., Flatau et al.,
c or Kahnert et al.:
c-----------------------------------------------------------------------
      a110=(aw(0,-1,0)+aw(0,1,0))
      do iw=0,Forder
         do k=-1,1
            do i=0,2
               aw(i,k,iw)=aw(i,k,iw)/a110
               bw(i,k,iw)=bw(i,k,iw)/a110
            enddo
         enddo
         cw(iw)=cw(iw)/a110
         dw(iw)=dw(iw)/a110
         a11=aw(0,-1,iw)+aw(0,1,iw)
         snr=a11/dble(iw+iw+1)
         if (iw.eq.1) a111=snr
      enddo
c-----------------------------------------------------------------------
c output file operations
c-----------------------------------------------------------------------
      if(fout.ne.' ') then
         open(1,file=fout,access='sequential')
         open(2,file=cout,access='sequential')
         open(3,file=eout,access='sequential')

         cext=qet
         csca=qst
         cabs=cext-csca
         w=csca/cext
         g=a111
      endif
c-----------------------------------------------------------------------
c scattering matrix calculation
c-----------------------------------------------------------------------
100   do i=1,nt
         if(nt.eq.1) then
            th=tmin
         else
            th=tmin+(tmax-tmin)*dble(i-1)/dble(nt-1)
         endif
         ct=dcos(th*pi/180.d0)
c-----------------------------------------------------------------------
c  dc is the normalized generalized spherical function
c  dc(k,n*(n+1)+m) = ((n-k)!(n+m)!/(n+k)!/(n-m)!)^(1/2) D^k_{mn},
c  where D^k_{mn} is defined in M&M JOSA 96
c-----------------------------------------------------------------------
         call rotcoef(ct,2,Forder,dc,2,nmax2,nbc,nfd,bcof,fnr)
c         write(6,*)'i,ct:',i,ct
c         if (i.eq.3) then
c            do ll=0,Forder*(Forder+1)
c               write(6,*)'ll,dc: ',ll,dc(0,ll)
c            end do
c         end if
         s11t=0.
         s12t=0.
         s13t=0.
         s14t=0.
         s21t=0.
         s22t=0.
         s23t=0.
         s24t=0.
         s31t=0.
         s32t=0.
         s33t=0.
         s34t=0.
         s41t=0.
         s42t=0.
         s43t=0.
         s44t=0.
         do n=0,Forder
            nn0=n*(n+1)
            nnp2=nn0+2
            nnm2=nn0-2
c            write(6,*)'aw-1,aw1,dc:',aw(0,-1,n),aw(0,1,n),dc(0,nn0)
            s11t=s11t+dc(0,nn0)*(aw(0,-1,n)+aw(0,1,n))
            s14t=s14t+dc(0,nn0)*(aw(0,1,n)-aw(0,-1,n))
            s44t=s44t+dc(0,nn0)*(bw(0,1,n)-bw(0,-1,n))
            s41t=s41t+dc(0,nn0)*(bw(0,-1,n)+bw(0,1,n))
            if(n.ge.2) then
               s12t=s12t+dc(2,nn0)*(aw(2,-1,n)+aw(2,1,n))
               s24t=s24t+dc(2,nn0)*(aw(2,1,n)-aw(2,-1,n))
               s34t=s34t+dc(2,nn0)*dimag(bw(2,-1,n)-bw(2,1,n))
               s31t=s31t-dc(2,nn0)*dimag(bw(2,-1,n)+bw(2,1,n))
               s13t=s13t+2.*dc(2,nn0)*dimag(aw(2,0,n))
               s42t=s42t+dc(2,nn0)*2.*bw(2,0,n)
               s22t=s22t+dc(2,nnp2)*cw(n)+dc(2,nnm2)*dw(n)
               s23t=s23t+dimag(dc(2,nnp2)*cw(n)+dc(2,nnm2)*dw(n))
               s33t=s33t+dc(2,nnp2)*cw(n)-dc(2,nnm2)*dw(n)
               s32t=s32t-dimag(dc(2,nnp2)*cw(n)-dc(2,nnm2)*dw(n))
            endif
         enddo
c-----------------------------------------------------------------------
c here are the VV and HH differential cross sections
c
c         gvv=.25*(s11t+s22t-2.*s12t)
c         ghh=.25*(s11t+s22t+2.*s12t)
c
c only the diagonal and upper triangular elements of the
c scattering matrix are written
c-----------------------------------------------------------------------
         if(fout.ne.' ')
     1    write(1,'(f9.3,2x,6e12.4)') th,s11t,s22t/s11t,s33t/s11t,
     1      s44t/s11t,s12t/s11t,-s34t/s11t
c     1    write(1,'(f9.3,2x,10e12.4)') th,s11t,s22t,s33t,s44t,s12t,s32t,
c     1    s34t,s31t,s42t,s41t

         if(dabs(th-180d0).lt.1e-6) then
            cbak=csca*s11t/4d0/pi
            deltaL=(s11t-s22t)/(s11t+s22t)
            deltaC=(s11t+s44t)/(s11t-s44t)
         endif
      enddo
      if(fout.ne.' ') then
c         write(1,'('' scattering matrix expansion coefficients'')')
c         write(1,'(''    w  a11         a22m        a22p        '',
c     1    ''a23m        a23p        a44         a12         '',
c     1    ''a34         a13         a24         a14'')')
         do iw=0,Forder
            a11=dble(aw(0,-1,iw)+aw(0,1,iw))
            a22=0.5d0*(cw(iw)+dw(iw))
            a33=0.5d0*(cw(iw)-dw(iw))
            a22mr=dw(iw)
            a22pr=cw(iw)
            a23mi=dimag(dw(iw))
            a23pi=dimag(cw(iw))
            a44=-dble(bw(0,-1,iw)-bw(0,1,iw))
            a12=dble(aw(2,-1,iw)+aw(2,1,iw))
            a34=dimag(bw(2,-1,iw)-bw(2,1,iw))
            a13=2.*dimag(aw(2,0,iw))
            a24=-dble(aw(2,-1,iw)-aw(2,1,iw))
            a14=-dble(aw(0,-1,iw)-aw(0,1,iw))
            iww=iw+1
            a1(iww)=a11
            a2(iww)=2.d0*a22
            a3(iww)=2.d0*a33
            a4(iww)=a44
            b1(iww)=-a12
            b2(iww)=-a34
            write(3,'(i5,6e24.14)') iw,a1(iww),a2(iww),a3(iww),
     &           a4(iww),b1(iww),b2(iww)
         enddo

         qext=cext/pi/Rv**2
         qsca=csca/pi/Rv**2
         qabs=cabs/pi/Rv**2
         qbak=cbak/pi/Rv**2

         write(2,7001) cext
         write(2,7002) csca
         write(2,7003) cabs
         write(2,7004) cbak
         write(2,7005) qext
         write(2,7006) qsca
         write(2,7007) qabs
         write(2,7008) qbak
         write(2,7009) w
         write(2,7010) g
         write(2,7011) deltaL
         write(2,7012) deltaC

         close(1)
         close(2)
         close(3)
      endif

 7001 format(e17.8,'  Cext')
 7002 format(e17.8,'  Csca')
 7003 format(e17.8,'  Cabs')
 7004 format(e17.8,'  Cbak')
 7005 format(e17.8,'  Qext')
 7006 format(e17.8,'  Qsca')
 7007 format(e17.8,'  Qabs')
 7008 format(e17.8,'  Qbak')
 7009 format(e17.8,'  SSA')
 7010 format(e17.8,'  g')
 7011 format(e17.8,'  delta_L')
 7012 format(e17.8,'  delta_C')
      return
      end



c***********************************************************************
c Calculation of the expansion coefficients of the random-orientation 
c averaged scattering matrix 
c
c Some minor modifications were implemented for exploiting
c Cn-symmetry in the orientation-averaging (M. Kahnert)
c
c***********************************************************************
      subroutine tmsm(Tmat,Pgroup,Nsym,nodr,nodrw,lmax2,
     &           nod2,nbd2,nbc,bcof,fnr,aw,bw,cw,dw)
      implicit real*8(a-h,o-z)
      integer Nsym,nodr,nodrw,nod2,nbd2,nbc,lmax2
c      integer Tmat_sym(2,lmax2,2,lmax2)
      complex*16 aw(0:2,-1:1,0:nod2),bw(0:2,-1:1,0:nod2),
     1           cw(0:nod2),dw(0:nod2)
      complex*16 ci,cin,bm(2,lmax2,2),pp(nodr),
     1           dm(-nodr-1:nodr+1,3,nodr,2,nodr,2),a,am(2,nodr+1,2)
      complex*16 fm(3,nodr,2,nodr,2)
      complex*16 Tmat(2,lmax2,2,lmax2)
      real*8 vc(0:2*nod2+2)
      real*8 bcof(0:nbc,0:nbc),fnr(0:2*nbc)
      character*3 Pgroup
      data ci/(0.d0,1.d0)/
c
      nodr2=nodr+nodr
      nblk=nodr*(nodr+2)
      do n=1,nodr
         do ip=1,2
            do l=1,nodr
               do iq=1,2
                  do k=1,3
                     do iu=-nodr-1,nodr+1
                        dm(iu,k,n,ip,l,iq)=0.
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      do n=1,nodr
         cin=ci**(n+1)
         pp(n) =-.5d0*cin*fnr(n+n+1)
      enddo
      nblk=nodr*(nodr+2)
      do n=1,nblk
         do l=1,nodr
            do k=-l,l
               kl=l*(l+1)+k
               do ip=1,2
                  a=Tmat(ip,n,1,kl)
                  Tmat(ip,n,1,kl)=pp(l)*
     1                 (Tmat(ip,n,1,kl)-Tmat(ip,n,2,kl))
                  Tmat(ip,n,2,kl)=pp(l)*(-a-Tmat(ip,n,2,kl))
               enddo
            enddo
         enddo
      enddo
c
      write(*,'(''      calculating D matrix...'')')
      do iw=0,nodr2
c         write
c     &      (*,'(''         expansion order:'',i3,''/'',i3)') iw,nodr2
         do iv=-iw,iw
c----------------------------------------------------------------------- 
c     For CN-symmetry only those am- and bm-terms are nonzero, for which 
c     |iv| = 0,N,2N,...
c----------------------------------------------------------------------- 
            if(Nsym.gt.1.and.Pgroup.ne.'Sn'.and.mod(iv,Nsym).ne.0) 
     &           goto 900
c----------------------------------------------------------------------- 
            do n=1,nblk
               do ip=1,2
                  do k=1,2
                     bm(k,n,ip)=0.
                  enddo
               enddo
            enddo
            do n=1,nodr
               nn1=n*(n+1)
               do l=max(1,abs(iw-n)),min(nodr,iw+n)
                  am(1,l,1)=0.
                  am(1,l,2)=0.
                  am(2,l,1)=0.
                  am(2,l,2)=0.
               enddo
               do it=-n,n
                  itn=nn1+it
                  lmax=min(nodr,iw+n)
                  call vcfunc(nbc,bcof,fnr,iv,iw,-it,n,lmax,vc)
                  do l=max(1,abs(iv-it),abs(n-iw)),lmax
                     ll1=l*(l+1)
                     itvl=ll1+it-iv
                     do ik=1,2
                        do ip=1,2
                           am(ik,l,ip)=am(ik,l,ip)
     1                          +vc(l)*Tmat(ip,itn,ik,itvl)
                        enddo
                     enddo
                  enddo
               enddo
c
               do m=-n,n
                  mn=nn1+m
                  do ik=1,2
                     k=-3+2*ik
                     iu=m-k
                     if(abs(iu).le.iw) then
                        lmax=min(nodr,iw+n)
                        call vcfunc(nbc,bcof,fnr,-iu,iw,m,n,
     &                       lmax,vc)
                        do l=max(1,abs(iw-n)),lmax
                           fl=-(-1)**l*vc(l)/dble(l+l+1)
                           do ip=1,2
                              bm(ik,mn,ip)=bm(ik,mn,ip)
     1                         +am(ik,l,ip)*fl
                           enddo
                        enddo
                     endif
                  enddo
               enddo
            enddo
            do iu=-min(iw,nodr+1),min(iw,nodr+1)
               do iku=1,3
                  if(iku.eq.1) then
                     k=-1
                     k1=-1
                  elseif(iku.eq.2) then
                     k=1
                     k1=1
                  else
                     k=1
                     k1=-1
                  endif
                  m=iu+k
                  ns=max(1,abs(m))
                  ik=(k+1)/2+1
                  ik1=(k1+1)/2+1
                  m1=iu+k1
                  do n=ns,nodr
                     nu=n*(n+1)+m
                     n1s=max(1,abs(m1),n-nodrw)
                     n1e=min(nodr,n+nodrw)
                     do n1=n1s,n1e
                        cin=ci**(n-n1)
                        nu1=n1*(n1+1)+m1
                        fnn1=-fnr(n+n+1)*fnr(n1+n1+1)*dble(iw+iw+1)
                        do ip=1,2
                           do ip1=1,2
                              a=bm(ik,nu,ip)*cin*fnn1
     1                             *conjg(bm(ik1,nu1,ip1))
                              dm(iu,iku,n,ip,n1,ip1)
     1                            = dm(iu,iku,n,ip,n1,ip1)+a
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
 900        continue
         enddo
      enddo
c
      write(*,'(''      calculating scattering matrix expansion...'')')
      do iw=0,nodrw
         do k=-1,1
            do i=0,2
               aw(i,k,iw)=0.
               bw(i,k,iw)=0.
            enddo
         enddo
         cw(iw)=0.
         dw(iw)=0.
      enddo
      do iw=0,nodrw
         iu2=1
         do n=1,nodr
            n1s=max(1,abs(n-iw))
            n1e=min(nodr,n+iw)
            do n1=n1s,n1e
               do ik=1,3
                  do ip=1,2
                     do ip1=1,2
                        fm(ik,n,ip,n1,ip1)=0.
                     enddo
                  enddo
               enddo
            enddo
         enddo
         do iu=-nodr-1,nodr+1
            do k=-1,1,2
               m=iu+k
               ik=(k+1)/2+1
               ns=max(1,abs(m))
               do n=ns,nodr
                  n1max=min(iw+n,nodr)
                  call vcfunc(nbc,bcof,fnr,m,n,0,iw,n1max,vc)
                  do n1=ns,nodr
                     if((n+n1.lt.iw).or.(abs(n-n1).gt.iw)) goto 60
                     fc1=-(-1)**n*vc(n1)*fnr(iw+iw+1)/fnr(n1+n1+1)
                     do ip=1,2
                        do ip1=1,2
                           fm(ik,n,ip,n1,ip1)=fm(ik,n,ip,n1,ip1)
     1                     +dm(iu,ik,n,ip,n1,ip1)*fc1
                        enddo
                     enddo
60                enddo
               enddo
            enddo
            if(iw.lt.2) goto 75
            m=iu+1
            m1=iu-1
            ns=max(1,abs(m))
            n1s=max(1,abs(m1))
            do n=ns,nodr
               n1max=min(iw+n,nodr)
               call vcfunc(nbc,bcof,fnr,m,n,-2,iw,n1max,vc)
               do n1=n1s,nodr
                  if((n+n1.lt.iw).or.(abs(n-n1).gt.iw)) goto 70
                  fc1=-(-1)**n*vc(n1)*fnr(iw+iw+1)/fnr(n1+n1+1)
                  do ip=1,2
                     do ip1=1,2
                        fm(3,n,ip,n1,ip1)=fm(3,n,ip,n1,ip1)
     1                   +dm(iu,3,n,ip,n1,ip1)*fc1
                     enddo
                  enddo
70             enddo
            enddo
            iu2=1
75       enddo
63       do n=1,nodr
            n1s=max(1,abs(n-iw))
            n1e=min(nodr,n+iw)
            in=(-1)**n
            n1max=min(iw+n,nodr)
            call vcfunc(nbc,bcof,fnr,1,n,0,iw,n1max,vc)
            do n1=n1s,n1e
               fc1=2.d0*in*vc(n1)*fnr(iw+iw+1)/fnr(n1+n1+1)
               i=mod(n+n1-iw,2)+1
               do ip=1,2
                  ip1=(2-i)*ip+(i-1)*(3-ip)
                  do k=-1,1,2
                     ik=(k+1)/2+1
                     aw(0,k,iw)=aw(0,k,iw)+fm(ik,n,ip,n1,ip1)*fc1
                     bw(0,k,iw)=bw(0,k,iw)+fm(ik,n,ip,n1,3-ip1)*fc1
                  enddo
                  bw(2,0,iw)=bw(2,0,iw)+fm(3,n,ip,n1,3-ip1)*fc1
                  aw(2,0,iw)=aw(2,0,iw)+fm(3,n,ip,n1,ip1)*fc1
               enddo
            enddo
            if(iw.lt.2) goto 80
            call vcfunc(nbc,bcof,fnr,1,n,-2,iw,n1max,vc)
            do n1=n1s,n1e
               fc2=2.d0*in*vc(n1)*fnr(iw+iw+1)/fnr(n1+n1+1)
               i=mod(n+n1-iw,2)+1
               do ip=1,2
                  ip1=(2-i)*ip+(i-1)*(3-ip)
                  do k=-1,1,2
                     ik=(k+1)/2+1
                     aw(2,k,iw)=aw(2,k,iw)+fm(ik,n,ip,n1,ip1)*fc2
     1                          *(-1)**ip1
                     bw(2,k,iw)=bw(2,k,iw)+fm(ik,n,ip,n1,3-ip1)*fc2
     1                          *(-1)**(3-ip1)
                  enddo
               enddo
               fc3=2.*(-1)**(n1+iw)*vc(n1)*fnr(iw+iw+1)/fnr(n1+n1+1)
               do ip=1,2
                  do ip1=1,2
                     cw(iw)=cw(iw)+fm(3,n,ip,n1,ip1)*fc2*(-1)**ip1
                     dw(iw)=dw(iw)+fm(3,n,ip,n1,ip1)*fc3*(-1)**ip
                  enddo
               enddo
            enddo
80       enddo
      enddo
      return
      end


c***********************************************************************
c     Calculation of the functions dc, where
c     dc is the normalized generalized spherical function
c     dc(k,n*(n+1)+m) = ((n-k)!(n+m)!/(n+k)!/(n-m)!)^(1/2) D^k_{mn},
c     where D^k_{mn} is defined in M&M JOSA 96.
c     This routine calculates dc(k,n*(n+1)+m) for k=-ndim..ndim and
c     n=0...nmax (where we have nmax=ndim throughout this programme).
c***********************************************************************
      subroutine rotcoef(cbe,nmax,Forder,dc,ndim,nmax2,nbc,nfd,bcof,fnr)
      implicit none 
      integer nmax,Forder,ndim,nmax2,nbc,nfd
      real*8 cbe,dc(-ndim:ndim,0:nfd)
      real*8 dk0(-nmax2:nmax2),dk01(-nmax2:nmax2)
      real*8 bcof(0:nbc,0:nbc),fnr(0:2*nbc)
      integer in,n,kmax,nn1,k,kn,im,m,m1
      real*8 sbe,sben,dkt,fmn,dkm0,dkm1,dkn1
      complex*16 ci,cbe2,sbe2
      data ci/(0.d0,1.d0)/
      
c-----------------------------------------------------------------------
c calculation of rotation functions.
c D^k_{mn}(cos beta) = dc(k,n*(n+1)+m)
c-----------------------------------------------------------------------
      sbe=dsqrt((1.d0+cbe)*(1.d0-cbe))
      cbe2=.5d0*(1.d0+cbe)
      sbe2=.5d0*(1.d0-cbe)
      in=1
      dk0(0)=1.d0
      sben=1.d0
      dc(0,0)=1.d0
      dk01(0)=0.d0

      do n=1,Forder
         kmax=min(n,nmax)
         nn1=n*(n+1)
         in=-in
         sben=sben*sbe/2.d0
         dk0(n)=in*sben*bcof(n,n)
         dk0(-n)=in*dk0(n)
         dk01(n)=0.
         dk01(-n)=0.
         dc(0,nn1+n)=dk0(n)
         dc(0,nn1-n)=dk0(-n)
         do k=-n+1,n-1
            kn=nn1+k
            dkt=dk01(k)
            dk01(k)=dk0(k)
            dk0(k)=(cbe*dble(n+n-1)*dk01(k)-fnr(n-k-1)*fnr(n+k-1)*dkt)
     1             /(fnr(n+k)*fnr(n-k))
            dc(0,kn)=dk0(k)
         enddo
         im=1
         do m=1,kmax
            im=-im
            fmn=1.d0/fnr(n-m+1)/fnr(n+m)
            m1=m-1
            dkm0=0.
            do k=-n,n
               kn=nn1+k
               dkm1=dkm0
               dkm0=dc(m1,kn)
               if(k.eq.n) then
                  dkn1=0.
               else
                  dkn1=dc(m1,kn+1)
               endif
               dc(m,kn)=(fnr(n+k)*fnr(n-k+1)*cbe2*dkm1
     1           -fnr(n-k)*fnr(n+k+1)*sbe2*dkn1
     1              -dble(k)*sbe*dc(m1,kn))*fmn
               dc(-m,nn1-k)=dc(m,kn)*(-1)**(k)*im
            enddo
         enddo
      enddo
      return
      end



c***********************************************************************
c     Vector coupling coefficients vc(iw) = C(m,n|k,l|m+k,iw)
c     uses an upwards recurrence
c***********************************************************************
      subroutine vcfunc(nbc,bcof,fnr,m,n,k,l,wmax,vcn)
      implicit none 
      integer nbc,m,n,k,l
      real*8 vcn(0:*),bcof(0:nbc,0:nbc),fnr(0:2*nbc)
      integer w,wmax,w1,w2,mk,nl,m1,n1,l1,k1
      real*8 vc1,t1,t2,t3
      
      mk=abs(m+k)
      nl=abs(n-l)
      if(nl.ge.mk) then
         w=nl
         if(n.ge.l) then
            m1=m
            n1=n
            l1=l
            k1=k
         else
            m1=k
            n1=l
            k1=m
            l1=n
         endif
         vc1=(-1)**(k1+l1)*bcof(l1+k1,w-m1-k1)
     1     *bcof(l1-k1,w+m1+k1)/bcof(l1+l1,w+w+1)
      else
         w=mk
         if(m+k.ge.0) then
            vc1=(-1)**(n+m)*bcof(n-l+w,l-k)*bcof(l-n+w,n-m)
     1          /bcof(w+w+1,n+l-w)
         else
            vc1=(-1)**(l+k)*bcof(n-l+w,l+k)*bcof(l-n+w,n+m)
     1          /bcof(w+w+1,n+l-w)
         endif
      endif
      w1=w
      vcn(w)=vc1
      w=w1+1
      mk=m+k
      w2=min(wmax,n+l)
      if(w2.gt.w1) then
         t1=2*w*fnr(w+w+1)*fnr(w+w-1)/(fnr(w+mk)*fnr(w-mk)
     1     *fnr(n-l+w)*fnr(l-n+w)*fnr(n+l-w+1)*fnr(n+l+w+1))
         if(w1.eq.0) then
            t2=.5*dble(m-k)
         else
            t2=dble((m-k)*w*(w-1)-mk*n*(n+1)+mk*l*(l+1))
     1       /dble(2*w*(w-1))
         endif
         vcn(w)=t1*t2*vcn(w1)
      endif
      do w=w1+2,w2
         t1=2*w*fnr(w+w+1)*fnr(w+w-1)/(fnr(w+mk)*fnr(w-mk)
     1     *fnr(n-l+w)*fnr(l-n+w)*fnr(n+l-w+1)*fnr(n+l+w+1))
         t2=dble((m-k)*w*(w-1)-mk*n*(n+1)+mk*l*(l+1))
     1    /dble(2*w*(w-1))
         t3=fnr(w-mk-1)*fnr(w+mk-1)*fnr(l-n+w-1)*fnr(n-l+w-1)
     1     *fnr(n+l-w+2)*fnr(n+l+w)/(dble(2*(w-1))*fnr(2*w-3)
     1     *fnr(2*w-1))
         vcn(w)=t1*(t2*vcn(w-1)-t3*vcn(w-2))
      enddo
      return
      end



