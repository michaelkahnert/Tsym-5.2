c=======================================================================
c     Tsym Version 5.2
c     
c     (C) Michael Kahnert
c     2012-09-27
c
c     T-matrix program based on Waterman's null-field method, 
c     exploiting point-group symmetries. 
c
c----------------------------------------------------------------------
c
c     SHORT DESCRIPTION:
c     The null-field routines in this program can, in principle, compute 
c     scattering by homogeneous, star-shaped particles of arbitrary 
c     shape. Currently, the following geometries are implemented:
c     (i)   Regular polygonal prisms with an arbitrary No. of corners
c           (includes circular cylinders as limiting case).
c     (ii)  2D Chebyshev particles of arbitrary order
c     (iii) 3D Chebyshev particles of arbitrary order
c     (iv)  2D Gaussian random spheres
c     (v)   3D Gaussian random spheres 
c     (vi)  Gaussian random sphere/Chebyshev hybrid particles
c     (vii) Spheres
c
c     USER GUIDE:
c     Documentation on how to use this program can be found in the 
c     subdirectory ./DOC/. 
c----------------------------------------------------------------------
c
c     Free software license 
c     
c     The routines in this file (Tsym.F) are free software for 
c     computing light scattering by star-shaped, homogeneous particles. 
c     These routines are available under 
c     the GNU General Public License that you can find on 
c     the World Wide Web (http://www.gnu.org/licenses/gpl.txt).
c     
c     Contact address for the author:
c     
c     Michael Kahnert
c     Swedish Meteorological and Hydrological Institute
c     Research Department
c     Folkborgsvägen 1
c     601 76 Norrköping
c     Sweden
c     E-mail: michael.kahnert@smhi.se
c     
c     The routines in this file are free software; you can redistribute 
c     them and/or modify them under the terms of the GNU General Public 
c     License as published by the Free Software Foundation; 
c     either version 2 of the License, or (at your option) any later 
c     version.
c     
c     This program is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c     GNU General Public License for more details.
c     
c     You should have received a copy of the GNU General Public License
c     along with this program; if not, go to the webpage
c     http://www.gnu.org/licenses/gpl.txt, or write to the Free
c     Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
c     MA  02111-1307  USA
c
c-----------------------------------------------------------------------
c     NOTE: The Tsym code makes use of other program packages, which,
c     for the user's convenience, are included in the distribution of
c     Tsym. The license terms of these packages can be found in the
c     respective files (see acknowledgements below).
c-----------------------------------------------------------------------
c
c     HOW TO CITE THIS PROGRAM:
c     When publishing results obtained with this program, please cite 
c     the following reference:
c
c        Michael Kahnert. Irreducible representations of finite groups
c        in the T-matrix formulation of the electromagnetic scattering
c        problem. J. Opt. Soc. Am. A 22, 1187-1199, 2005.
c
c     Other key papers are:
c
c     Iterative Lippmann-Schwinger approach for particles with small-
c     scale surface roughness:

c        Michael Kahnert and Tom Rother. Modeling optical properties of 
c        particles with small-scale surface roughness: combination of 
c        group theory with a perturbation approach," Opt. Express 19, 
c        11138-11151, 2011. 
c
c     Application of Tsym to Chebyshev particles:
c
c        Michael Kahnert, Timo Nousiainen, and Päivi Mauno. On the 
c        impact of non-sphericity and small-scale surface roughness 
c        on the optical properties of hematite aerosols. JQSRT 112,
c        1815–1824, 2011.
c
c     Application of Tsym to Gaussian random spheres
c
c        Michael Kahnert, Timo Nousiainen, Manu Anna Thomas, and
c        Jani Tyynelä. Light scattering by particles with small-scale
c        surface roughness: Comparison of four classes of model 
c        geometries. JQSRT 2012, 
c        DOI: http://dx.doi.org/10.1016/j.jqsrt.2012.03.017.
c
c-----------------------------------------------------------------------
c
c     Recent changes and additions:
c
c     V 5.2: September 2012
c            Added version handling of params file;
c            Corrected error in rotation of the T-matrix;
c            Replaced lpropack with integer flag Pmatmethod to allow
c            for three options in the construction of the transformation
c            into the irreducible basis
c     V 5.1: May 2012
c            Implementation of PROPACK singular value decomposition
c            for sparse matrices. Introduction of a flag 'lpropack'
c            to choose between Propack and Lapack. Subroutine
c            IrredRep moved to the start of T-matrix computations
c            for more efficient use of computer memory.
c            Old character table format removed (only GAP format used now!)
c     V 5.0: May 2012
c            Improved interface to GAP character tables. (The old version
c            did not work for high-order point groups, since the ordering
c            of the conjugacy classes in GAP does not follow a simple scheme
c            anymore.) Improved dymamic memory allocation in IrredRep. P-matrix
c            inversion replaced by adjoint of P-matrix.
c     V 4.9: April 2012
c            Use of LAPACK singular value decomposition routine
c            to construct the transformation into the irreducible
c            basis from the projectors into the invariant subspaces.
c     V 4.8: April 2012
c            Character tables are read in GAP output format. A script 
c            for generating character tables with GAP has been written.
c     V 4.7: March 2012
c            Parameter file params now read in while running, rather 
c            than included while compiling.
c     V 4.6: March 2012
c            Added GRS/Chebyshev hyrid particles (geom='GRSCHB'). 
c            Commented out random rotation in GRS.F.
c     V 4.5: January 2012
c            Added extinction and scattering cross section for particles
c            in fixed orientations.
c     V 4.4: November 2011
c            Added 2D-Gaussian random spheres (geom='GRSP2D'). Improved
c            computation of derivatives of Legendre functions in GRS.F.
c     V 4.3: March 2011
c            Minor streamlining of code.
c     V 4.2: January 2011
c            Added Gaussian random spheres (adapted the GSPHERE code
c            by Karri Muinonen and Timo Nousiainen). 
c     V 4.1: November 2010
c            Changed Gauss-Legendre scheme from global scheme to 
c            sub-interval scheme. Added case-distinction in quadrature
c            scheme for axisymmetric particles. 
c     V 4.0: November 2010
c            Introduced parameter Geom to discriminate among different
c            classes of geometries. Added 2D- and 3D-Chebyshev particles
c            and spheres (Geom='CHEB2D','CHEB3D','SPHERE').
c            Added iterative computation of T-matrix via Lippmann
c            Schwinger equation. 
c     V 3.1: August 2010
c            Added computation of polarised differential scattering
c            cross sections for particle in fixed orientation; extended
c            computation of Wigner d-functions to \theta\in (180,360].
c            Added surface parameterisation of cube in standard orientation.
c            Added reducible representations of cubic point groups
c            (not fully tested yet!). 
c     V 2.0: 2004 
c            Added irreducible representations
c
c-----------------------------------------------------------------------
c
c     Acknowledgements
c
c     To run Tsym, you need several other routines that are not part
c     of the Tsym code proper. For the user's convenience, these
c     programs have been included in the distribution of Tsym in
c     the following files. (Note that some of these programs have
c     been slightly edited, extended, or modified to fit the needs 
c     of Tsym. Details of these modifications are given in the comments
c     of the respective subroutines.)
c
c     functions.F : Routines for computing Wigner d-functions, Bessel
c                   functions, and Gaussian quadrature points, taken
c                   from the T-matrix code by Michael Mishchenko
c                   (M. Mishchenko, JQSRT 60, 309-324, 1998).
c     gsphere.F   : Program for computing Gaussian random sphere 
c                   geometries, written by Karri Muinonen and Timo
c                   Nousiainen (K. Muinonen et al.,JQSRT 55, 577-601, 1996)
c     lapack.F    : Lapack routines for LU decomposion and SVD 
c                   (www.netlib.org/lapack).
c     OrientAvg.F : Routines for computing ensemble-averaged scattering
c                   matrix and vector coupling coefficients, written by
c                   Dan Mackowski (Mackowski and Mishchenko, JOSA A 13, 
c                   2266-2278, 1996).
c     propack.F   : Propack package for SVD of large sparse matrices,
c                   written by Rasmus Munk Larsen 
c                   (http://sun.stanford.edu/~rmunk/PROPACK).
c     smatrix.F   : Routines for computing the amplitude scattering matrix
c                   and Mueller matrix for particles in fixed orientation,
c                   written by Heikki Laitinen and Kari Lumme
c                   (Laitinen and Lumme, JQSRT 60, 325-334, 1998)
c
c     Special thanks to Kari Lumme, Dan Mackowski, Michael Mishchenko,
c     and Karri Muinonen for their kind permission to include these 
c     subroutines in the distribution of Tsym!
c  
c-----------------------------------------------------------------------
c
c     Short overview over the source code contained in this file:
c
c     Part I  : Tsym main program
c     Part II : Null-field method
c     Part III: Group theory
c     Part IV : Particle geometries
c     Part V  : Miscellaneous routines
c
c=======================================================================



c=======================================================================
c     PART I :  T S Y M   M A I N   P R O G R A M: 
c=======================================================================

c***********************************************************************
      program Tsym
      implicit none
c***********************************************************************
c
c-----------------------------------------------------------------------
c     variables read in from input file 'params' 
c-----------------------------------------------------------------------
      real*8 Tver
      character*6 Geom
      logical lprtb,lreadgeom
      integer cheborder,prtb_order,ngrs
      character*3 Pgroup
      integer nmax,mmax
      integer th_nint,phi_nint
      real*8 lambda
      complex*16 mr
      integer Nsym
      integer igeom
      real*8 gpar1,gpar2,gpar3
      integer norient
      integer neuler
      real*8, dimension(:), allocatable :: alpha, beta, gamma
      logical Cyl_Quad
      logical lirredrep
      character*250 chtabfile
      integer Pmatmethod
      integer ntout
      real*8 tmin,tmax
c-----------------------------------------------------------------------
c     Local variables:
c-----------------------------------------------------------------------
      integer i,ip1,ip2,thptsx,phiptsx,ng1,ng2,Cn,
     & ng1dim,ng2dim,Forder,nmax2,nfd,nbt1,
     & thpts,thpts2,phipts,nbc,lmax2,ieuler,err,iostat,
     & nmmax,lmax,lpmax,dim,unit,nlines
      real*8 k,Rv,F(4,4),dCsca(4),thi,ths,phi,phs,cf
      real*8 Cext,Csca,Qext,Qsca
      real*8 pi, eps, h2,a
      real*8 alpha0,beta0,gamma0
      real*8 dtime
      character fout*10,cout*10,eout*10,dout*10
      logical lrotation,AXISYM
      logical DIHEDRAL,SIGMAh
      integer, dimension(:), allocatable :: nsi,nmx,mmx,mmn
      integer, dimension(:,:), allocatable :: lsi,lpsi
      complex*16, dimension(:,:,:,:), pointer     :: Tmat2
      complex*16, dimension(:,:,:,:), allocatable :: Tmat3
      complex*16 S1(2,2)

      parameter ( pi=3.1415926535897932385d0)
c-----------------------------------------------------------------------
c     The following settings usually do not need to be changed
c     unless certain problems with, e.g., insufficient memory occur:
c-----------------------------------------------------------------------
      integer maxNc,dim1,dim2,ngauss_per_int
      real*8 prec
      parameter (maxNc=999,ngauss_per_int=16)
      parameter (prec=1e-14)
c-----------------------------------------------------------------------
      INTERFACE
         SUBROUTINE NFM(Tmatc,nmaxc,mmaxc,nmmaxc,nmxc,lmaxc,lmax2c,
     &        lpmaxc,dimc,ng1c,ng2c,ng1dimc,ng2dimc,
     &        ngauss_per_intc,th_nintc,phi_nintc,
     &        thptsc,thpts2c,phiptsc,thptsxc,phiptsxc,Nsymc,
     &        Pgroupc,maxNcc,precc,
     &        ac,hc,mrc,kc,epsc,pic,nsic,lsic,lpsic,lirredrepc,
     &        chtabfilec,Pmatmethodc,AXISYMc,DIHEDRALc,SIGMAhc,
     &        Cyl_Quadc,Geomc,cheborderc,prtb_orderc,lprtbc,lreadgeomc,
     &        ngrsc,dtimec)
         complex*16, dimension(:,:,:,:), pointer     :: Tmatc
         integer nmaxc,mmaxc,nmmaxc,lmaxc,lmax2c,lpmaxc,dimc,thptsc,
     &        ngauss_per_intc,th_nintc,phi_nintc,
     &        thpts2c,phiptsc,thptsxc,phiptsxc,Nsymc,ng1c,ng2c,ng1dimc,
     &        ng2dimc,nsic(nmmaxc),lsic(-mmaxc:mmaxc,nmmaxc),
     &        lpsic(0:mmaxc,nmmaxc),nmxc(-mmaxc:mmaxc),maxNcc,
     &        prtb_orderc,cheborderc,ngrsc,Pmatmethodc
         real*8 kc,epsc,pic,ac,hc,dtimec,precc
         complex*16 mrc
         logical DIHEDRALc,SIGMAhc,AXISYMc,Cyl_Quadc,lirredrepc,lprtbc,
     &        lreadgeomc
         character*3 Pgroupc
         character*6 Geomc
         character*(*) chtabfilec
         END SUBROUTINE NFM
      END INTERFACE
c-----------------------------------------------------------------------


c      call logmsg('Start of program                    ')
      open(6,file="logfile")
c-----------------------------------------------------------------------
c     READ INPUT DATA:
c-----------------------------------------------------------------------
      unit=10
      nlines=5
      open(unit,file='params')
      call skip_lines(unit,nlines-1)
      read(unit,*) Tver
      call check_version(Tver,unit,nlines)
c
      rewind(unit)
      call skip_lines(unit,19)
      read(unit,*) Geom
      call skip_lines(unit,11)
      read(unit,*) cheborder
      read(unit,*) ngrs
      read(unit,*) lreadgeom
      call skip_lines(unit,46)
      read(unit,*) igeom
      read(unit,*) gpar1
      read(unit,*) gpar2
      read(unit,*) gpar3
      call skip_lines(unit,8)
      read(unit,*) Pgroup
      read(unit,*) Nsym
      call skip_lines(unit,3)
      read(unit,*) lambda
      read(unit,*) mr
      call skip_lines(unit,10)
      read(unit,*) norient
      call skip_lines(unit,17)
      read(unit,*) neuler
      allocate(alpha(neuler),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for alpha'
         stop
      endif
      allocate(beta(neuler),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for beta'
         stop
      endif
      allocate(gamma(neuler),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for gamma'
         stop
      endif
      read(unit,*) (alpha(i),i=1,neuler)
      read(unit,*) (beta(i),i=1,neuler)
      read(unit,*) (gamma(i),i=1,neuler)
      call skip_lines(unit,10)
      read(unit,*) nmax
      read(unit,*) mmax
      call skip_lines(unit,10)
      read(unit,*) th_nint
      read(unit,*) phi_nint
      call skip_lines(unit,11)
      read(unit,*) lprtb
      read(unit,*) prtb_order
      call skip_lines(unit,6)
      read(unit,*) Cyl_Quad
      call skip_lines(unit,14)
      read(unit,*) lirredrep
      read(unit,*) chtabfile
      read(unit,*) Pmatmethod
      call skip_lines(unit,8)
      read(unit,*) ntout
      read(unit,*) tmin
      read(unit,*) tmax

      close(unit)
c-----------------------------------------------------------------------
      lmax=mmax*(2*nmax-mmax+1)+nmax
      nmmax=nmax
      lpmax=nmax*(mmax+1)-(mmax*(mmax-1))/2
      dim=2*lmax
      nbc=6*nmmax+6
      nmax2=2*nmmax
      nfd=nmax2*(nmax2+2)
      nbt1=(nmmax+1)*(nmmax+3)
      lmax2=nmax*(nmax+2)      
      dtime=0d0      
      k=2*pi/lambda
      Forder = nmax2      
c-----------------------------------------------------------------------
c     allocate arrays:
c-----------------------------------------------------------------------
      allocate(nsi(nmmax),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for nsi'
         stop
      endif
      allocate(mmx(nmmax),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for mmx'
         stop
      endif
      allocate(mmn(nmmax),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for mmn'
         stop
      endif
      allocate(nmx(-mmax:mmax),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for nmx'
         stop
      endif
      allocate(lsi(-mmax:mmax,nmmax),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for lsi'
         stop
      endif
      allocate(lpsi(0:mmax,nmmax),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for lpsi'
         stop
      endif
c-----------------------------------------------------------------------
c     set up some parameters:
c-----------------------------------------------------------------------
      call setup(Geom,cheborder,prtb_order,lprtb,lirredrep,Nsym,
     &     chtabfile,igeom,gpar1,gpar2,gpar3,ip1,ip2,mr,k,
     &     lambda,nmax,mmax,nmmax,thpts,phipts,thpts2,
     &     ngauss_per_int,th_nint,phi_nint,thptsx,phiptsx,
     &     Pgroup,SIGMAh,DIHEDRAL,AXISYM,Cyl_Quad,pi,a,norient,h2,
     &     eps,ng1,ng2,ng1dim,ng2dim,neuler,alpha,beta,gamma,nsi,
     &     lsi,lpsi,nmx,mmx,mmn,ngrs,Rv)       
c-----------------------------------------------------------------------
c     Null-field method:
c-----------------------------------------------------------------------
      call NFM(Tmat2,nmax,mmax,nmmax,nmx,lmax,lmax2,lpmax,
     &     dim,ng1,ng2,ng1dim,ng2dim,
     &     ngauss_per_int,th_nint,phi_nint,
     &     thpts,thpts2,phipts,thptsx,phiptsx,Nsym,
     &     Pgroup,maxNc,prec,
     &     a,h2,mr,k,eps,pi,nsi,lsi,lpsi,lirredrep,chtabfile,Pmatmethod,
     &     AXISYM,DIHEDRAL,SIGMAh,Cyl_Quad,Geom,cheborder,prtb_order,
     &     lprtb,lreadgeom,ngrs,dtime)


c-----------------------------------------------------------------------
c     Calculate Stokes scattering matrix and optical cross sections:
c-----------------------------------------------------------------------
      if(norient.ne.1) then ! fixed orientation(s)
c-----------------------------------------------------------------------
c     Calculate Stokes scattering matrix for particle in fixed
c     orientations:
c-----------------------------------------------------------------------
         allocate(Tmat3(2,lmax2,2,lmax2),stat=err)
         if(err.ne.0)then
            write(6,*)'FATAL ERROR: allocation failed for Tmat3'
            stop
         endif

         do ieuler=1,neuler     ! loop over Euler angles
            write(6,*)'   Orientation No. ',ieuler,':'
            call foutname(fout,cout,eout,dout,ieuler)
c-----------------------------------------------------------------------
c     rotate the T-matrix
c-----------------------------------------------------------------------
            if (dabs(alpha(ieuler)-pi).ge.prec.or.
     &           dabs(beta(ieuler)).ge.prec.or.
     &           dabs(gamma(ieuler)-pi).ge.prec) then
               lrotation=.true.
               call logmsg('   rotating the T-matrix...         ')
               call timeit(' ',dtime)
               alpha0=alpha(ieuler)
               beta0=beta(ieuler)
               gamma0=gamma(ieuler)

               call rot_T
     &              (Tmat2,Tmat3,alpha0,beta0,
     &              gamma0,lmax2,nmmax,nbc,nfd,nmax2,
     &              Nsym,DIHEDRAL,SIGMAh)
               call timeit
     &              ('   ...done with T-matrix rotation       : ',dtime)
            else
               lrotation=.false.
            end if

            open(10,file=fout,status='unknown',
     &           access='sequential')
            open(12,file=cout,status='unknown',
     &           access='sequential')
            open(13,file=dout,status='unknown',
     &           access='sequential')

            call logmsg('   calculating diff`l scatt. CS...: ')
            call timeit(' ',dtime)

            cf=pi/180d0
            do i=1,ntout
               if(ntout.eq.1) then
                  ths=tmin
               else
                  ths=tmin+(tmax-tmin)*dble(i-1)/dble(ntout-1)
               endif
               ths=cf*ths
               Cn=1
               if (dabs(beta(ieuler)-0.d0).lt.1d-6) Cn=Nsym
               
               thi=0.d0*cf
               phi=0.d0*cf
               phs=0.d0*cf
c-----------------------------------------------------------------------
c     Amplitude scattering matrix:
c-----------------------------------------------------------------------
               if(lrotation)then
                  call AmplScat(Tmat3,nmax,lmax2,thi,phi,ths,phs,
     &                 Cn,DIHEDRAL,S1)
               else
                  call AmplScat(Tmat2,nmax,lmax2,thi,phi,ths,phs,
     &                 Cn,DIHEDRAL,S1)
               endif
c-----------------------------------------------------------------------
c     Stokes scattering matrix and polarised differential scattering
c     cross section:
c-----------------------------------------------------------------------
               call fmatrix(S1,F,dCsca)
               write (10,1001) ths/cf,F(1,1),F(2,2)/F(1,1),
     &              F(3,3)/F(1,1),F(4,4)/F(1,1),
     &              F(1,2)/F(1,1),F(3,4)/F(1,1)
               write (13,1002)ths/cf,dCsca(1),dCsca(2),
     &              dCsca(3),dCsca(4)
              
            end do

c-----------------------------------------------------------------------
c     Extinction and scattering cross sections:
c-----------------------------------------------------------------------
            if(lrotation)then
               call crosssection(Tmat3,nmax,lmax2,thi,phi,ths,phs,k,
     &           Cext,Csca)
            else
               call crosssection(Tmat2,nmax,lmax2,thi,phi,ths,phs,k,
     &           Cext,Csca)
            endif
            Qext=Cext/pi/Rv**2
            Qsca=Csca/pi/Rv**2

            write (12,1003) Cext
            write (12,1004) Csca
            write (12,1005) Cext-Csca
            write (12,1006) Qext
            write (12,1007) Qsca
            write (12,1008) Qext-Qsca

            call timeit('   ...done with Mueller matrix          : ',
     &           dtime)
         enddo                  ! loop over Euler orientations
         deallocate(Tmat3)
      end if
c-----------------------------------------------------------------------
      if (norient.ne.2) then ! random orientations
c-----------------------------------------------------------------------
c     Calculate Stokes scattering matrix for ensemble of randomly
c     oriented particles
c-----------------------------------------------------------------------
         call foutname(fout,cout,eout,dout,-1)
         call logmsg('   calculating Stokes matrix...     ')
         call timeit(' ',dtime)
         
         call ranprops(Tmat2,Pgroup,Nsym,lambda,
     *        nmmax,lmax2,Forder,nfd,nbc,nmax2,nbt1,
     *        a,Rv,pi,fout,cout,eout,ntout,tmin,tmax)
         
         call timeit('    ...done with Stokes scattering matrix :',
     &        dtime)
      endif

      close(10)
      close(13)
      
      call logmsg('End of program                      ')
      close(6)
 1001 format(e13.4,2x,6(e15.6,2x))
 1002 format(e13.4,2x,4(e15.6,2x))
 1003 format(e17.8,'  Cext')
 1004 format(e17.8,'  Csca')
 1005 format(e17.8,'  Cabs')
 1006 format(e17.8,'  Qext')
 1007 format(e17.8,'  Qsca')
 1008 format(e17.8,'  Qabs')
      
      end


c=======================================================================
c     PART II :   N U L L - F I E L D   M E T H O D :
c=======================================================================

      module pmat
      implicit none
      complex*16, pointer, dimension(:)       :: P,Pinv,P1
      complex*16, dimension(:,:), allocatable :: P2
      integer, pointer, dimension(:)          :: ip,ipinv,jp,ip1,jp1
      end module pmat


c***********************************************************************
c     Null-field method
c
c     Author: Michael Kahnert
c***********************************************************************
      subroutine NFM(Tmat2,nmax,mmax,nmmax,nmx,lmax,lmax2,lpmax,
     &    dim,ng1,ng2,ng1dim,ng2dim,
     &    ngauss_per_int,th_nint,phi_nint,
     &    thpts,thpts2,phipts,thptsx,phiptsx,Nsym,
     &    Pgroup,maxNc,prec,apoly,hpoly,
     &    mr,k,eps,pi,nsi,lsi,lpsi,lirredrep,chtabfile,Pmatmethod,
     &    AXISYM,DIHEDRAL,SIGMAh,Cyl_Quad,Geom,cheborder,prtb_order,
     &    lprtb,lreadgeom,ngrs,dtime)

      use pmat
      implicit none
c-----------------------------------------------------------------------
c     INPUT:
c-----------------------------------------------------------------------      
      integer nmax,mmax,nmmax,lmax,lmax2,lpmax,dim,thpts,thpts2,phipts,
     &     thptsx,phiptsx,Nsym,ng1,ng2,ng1dim,ng2dim,
     &     ngauss_per_int,th_nint,phi_nint,
     &     nsi(nmmax),lsi(-mmax:mmax,nmmax),lpsi(0:mmax,nmmax),
     &     nmx(-mmax:mmax),maxNc,cheborder,prtb_order,ngrs,
     &     Pmatmethod
      real*8 k,eps,pi,apoly,hpoly,dtime,prec
      complex*16 mr
      logical DIHEDRAL,SIGMAh,AXISYM,Cyl_Quad,lirredrep,lprtb,lreadgeom
      character*3 Pgroup
      character*6 Geom
      character*(*) chtabfile
c-----------------------------------------------------------------------
c     OUTPUT:
c-----------------------------------------------------------------------
      complex*16, dimension(:,:,:,:), pointer :: Tmat2
c-----------------------------------------------------------------------
c     LOCALLY DEFINED:
c-----------------------------------------------------------------------
      integer lwork,info,ipiv(dim),   i,j,
     & ngstart,n,nn,m,mm,l,ll,TQsym(lmax,lmax),lmax3,
     & Nc, No, Nclass(maxNc), Norder(maxNc,maxNc),
     & ne, neinv, ir_size(maxNc),bdim,err

      real*8 xth1(ng1dim),wth1(ng1dim),
     & xth2(ng2dim,phipts),wth2(ng2dim,phipts),
     & xphi(phipts),wphi(phipts),phic1(ng1dim),phic2(ng1dim),
     & r1(ng1dim),rth1(ng1dim),rphi1(ng1dim),
     & r2(ng2dim,phipts),
     & rth2(ng2dim,phipts),rphi2(ng2dim,phipts),
     & r0(ng2dim,phipts),
     & rth0(ng2dim,phipts),rphi0(ng2dim,phipts),
     & zero,alpha
     
      complex*16 ks,work(5*dim),czero,cone,Im
      complex*16, dimension(:,:), allocatable :: Q,Qb
      complex*16, dimension(:,:,:,:), allocatable :: Tmat,RgQ
      complex*16, dimension(:), allocatable :: Q0

      parameter (czero=(0d0,0d0), cone=(1d0,0d0), Im=(0d0,1d0)) 
      parameter ( zero=0d0)
      
c-----------------------------------------------------------------------
      call logmsg('   starting subroutine NFM...       ')
                  
      ks=k*mr

c=======================================================================
c     Transformation matrix into irreducible basis:
c-----------------------------------------------------------------------
      if(lirredrep)then
         call logmsg('     Compute P-matrix...                        ')
         call timeit(' ',dtime)
c-----------------------------------------------------------------------
c     compute the transformation matrix P for transforming the Q-matrix
c     into the irreducible basis
c-----------------------------------------------------------------------
         call IrredRep(Pgroup,Nsym,maxNc,dim,mmax,nmmax,lmax,nmx,
     &        lsi,Nc,No,Nclass,Norder,prec,ne,neinv,ir_size,chtabfile,
     &        Pmatmethod)
         call timeit('     ...done with P-matrix              : ',dtime)
      endif
c=======================================================================
c     Geometry calculations:
c=======================================================================
      call logmsg('     geometry calculations...       ')
      call timeit(' ',dtime)
      alpha = dble(Nsym-2)*pi/dble(2*Nsym)      
c-----------------------------------------------------------------------
c     Gauss-Legendre Quadrature points:
c-----------------------------------------------------------------------
      if (Cyl_Quad) then ! NB: It is made sure in setup that Cyl_Quad=.true.
                         ! for Geom='PRISMS' only.
         call GL_pts_cyl(eps,alpha,pi,ng1,ng2,
     &        ng1dim,ng2dim,SIGMAh,DIHEDRAL,          
     &        ngauss_per_int,th_nint,phi_nint,
     &        phipts,Nsym,xth1,wth1,xth2,wth2,xphi,wphi,phic1,
     &        phic2,ngstart)
      else
         call GL_pts(Geom,AXISYM,eps,alpha,pi,ng2dim,
     &        ngauss_per_int,th_nint,phi_nint,thpts2,
     &        phipts,Nsym,xth2,wth2,xphi,wphi,ng1,ng2)
      end if
c-----------------------------------------------------------------------
c     Parametrization r(theta,phi), as well as d r/d theta,
c     d r/d phi:
c----------------------------------------------------------------------- 
c     Polygonal prisms:
      if(Geom.eq.'PRISMS')then
         if (Cyl_Quad) then       
            call R_poly_cyl(xth1,xth2,xphi,alpha,ng1,ng2,
     &           ng1dim,ng2dim,phipts,Nsym,pi,apoly,hpoly,SIGMAh,
     &           r1,rth1,rphi1,r2,rth2,rphi2)
         elseif(Pgroup.eq.'Oh')then
            call R_cube (xth2,xphi,alpha,ng2dim,thpts2,ng1,ng2,
     &           phipts,Nsym,pi,apoly,hpoly,r2,rth2,rphi2)
         else
            call R_poly (xth2,xphi,alpha,ng2dim,thpts2,ng1,ng2,phipts,
     &           Nsym,pi,apoly,hpoly,r2,rth2,rphi2)
         end if

c     2D Chebyshev particles:
      elseif(Geom.eq.'CHEB2D')then 
         call R_Cheb2D (xth2,xphi,ng2dim,thpts2,
     &        phipts,apoly,hpoly,r2,rth2,rphi2,cheborder)  

c     3D Chebyshev particles:
      elseif(Geom.eq.'CHEB3D')then 
         call R_Cheb3D (xth2,xphi,ng2dim,thpts2,
     &        phipts,apoly,hpoly,r2,rth2,rphi2,cheborder)  

c     Spheres:
      elseif(Geom.eq.'SPHERE')then 
         call R_Sphere (xth2,xphi,ng2dim,thpts2,
     &        phipts,apoly,r2,rth2,rphi2)      
 
c     Gaussian Random Spheres:
      elseif(Geom.eq.'GRSPHR')then 
         call R_GRS (xth2,xphi,ng2dim,
     &        phipts,apoly,hpoly,eps,r2,rth2,rphi2,lreadgeom,ngrs)  
 
c     2D-Gaussian Random Spheres:
      elseif(Geom.eq.'GRSP2D')then 
         call R_GRS2D (xth2,xphi,ng2dim,
     &        phipts,apoly,hpoly,eps,r2,rth2,rphi2,lreadgeom,ngrs)
 
c     Gaussian Random Spheres/Chebyshev hybrids:
      elseif(Geom.eq.'GRSCHB')then 
         call R_GRSCheb (xth2,xphi,ng2dim,phipts,
     &        apoly,hpoly,eps,r2,rth2,rphi2,lreadgeom,ngrs,cheborder)  
      endif

c     also need base-sphere for perturbation approach:
      if(lprtb.and.(Geom.eq.'GRSPHR'.or.Geom.eq.'GRSP2D'
     &     .or.Geom.eq.'CHEB2D'.or.Geom.eq.'CHEB3D'
     &     .or.Geom.eq.'GRSCHB'))then 
         call R_Sphere (xth2,xphi,ng2dim,thpts2,
     &        phipts,apoly,r0,rth0,rphi0) 
      endif
   
c=======================================================================
c     Write geometry file for plotting with Matlab:
c=======================================================================
      if(.not.Cyl_Quad)then
         call write_matlab(xth2,xphi,ng2dim,phipts,r2,Nsym)
      else
         write(6,*)
     & '     NB: No geometry plotting files written for Cyl_Quad=.true.'
      endif
      call timeit('     ...done with geometry              : ',dtime)


c=======================================================================
c     calculation of the Q-matrices:
c=======================================================================  
      allocate(Q(dim,dim),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for Q'
         stop
      endif
      allocate(RgQ(2,lmax,2,lmax),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for RgQ'
         stop
      endif

      call logmsg('     Q-matrix calculations...       ')
      call timeit(' ',dtime)
      if(lprtb)then
         allocate(Q0(dim),stat=err)
         if(err.ne.0)then
            write(6,*)'FATAL ERROR: allocation failed for Q0'
            stop
         endif
         call Qmat_sphere(mmax,nmmax,nmx,lmax,lpmax,dim,thpts,thptsx,
     &          phipts,phiptsx,Nsym,DIHEDRAL,SIGMAh,ng1,ng2,ng1dim,
     &          ng2dim,xth1,xth2,xphi,r1,rth1,rphi1,r0,rth0,rphi0,
     &          wth1,wth2,wphi,k,ks,mr,pi,nsi,lsi,lpsi,zero,cone,czero,
     &          Im,dtime,Q,RgQ,TQsym) 
         do i=1,dim
            Q0(i)=Q(i,i)
         enddo
      endif

      if (Cyl_Quad) then       
         call Qmat_cyl(mmax,nmmax,nmx,lmax,lpmax,dim,phipts,
     &          phiptsx,Nsym,DIHEDRAL,SIGMAh,ng1,ng2,ng1dim,ng2dim,xth1,
     &          xth2,xphi,phic1,phic2,ngstart,r1,rth1,rphi1,r2,rth2,
     &          rphi2,wth1,wth2,wphi,k,ks,mr,pi,nsi,lsi,lpsi,zero,cone,
     &          czero,Im,dtime,Q,RgQ,TQsym)   
      elseif(Geom.eq.'SPHERE')then
         call Qmat_sphere(mmax,nmmax,nmx,lmax,lpmax,dim,thpts,thptsx,
     &          phipts,phiptsx,Nsym,DIHEDRAL,SIGMAh,ng1,ng2,ng1dim,
     &          ng2dim,xth1,xth2,xphi,r1,rth1,rphi1,r2,rth2,rphi2,
     &          wth1,wth2,wphi,k,ks,mr,pi,nsi,lsi,lpsi,zero,cone,czero,
     &          Im,dtime,Q,RgQ,TQsym) 
      elseif(Geom.eq.'CHEB2D'.or.Geom.eq.'GRSP2D')then
         call Qmat_axisym(mmax,nmmax,nmx,lmax,lpmax,dim,thpts,thptsx,
     &          phipts,phiptsx,Nsym,DIHEDRAL,SIGMAh,ng1,ng2,ng1dim,
     &          ng2dim,xth1,xth2,xphi,r1,rth1,rphi1,r2,rth2,rphi2,
     &          wth1,wth2,wphi,k,ks,mr,pi,nsi,lsi,lpsi,zero,cone,czero,
     &          Im,dtime,Q,RgQ,TQsym) 
      else
         call Qmat(mmax,nmmax,nmx,lmax,lpmax,dim,thpts,thptsx,
     &          phipts,phiptsx,Nsym,DIHEDRAL,SIGMAh,ng1,ng2,ng1dim,
     &          ng2dim,xth1,xth2,xphi,r1,rth1,rphi1,r2,rth2,rphi2,
     &          wth1,wth2,wphi,k,ks,mr,pi,nsi,lsi,lpsi,zero,cone,czero,
     &          Im,dtime,Q,RgQ,TQsym)         
      end if
      call timeit('     ...done with Q-matrix              : ',dtime)

c=======================================================================
c     Inversion of the Q-matrix
c
c     This can be done in three different ways:
c     1. Preconditioning the Q-matrix by transforming into the 
c        irreducible basis (you need to provide the Group character table
c        for this!)
c     2. Perturbation theory
c     3. Direct inversion of the Q-matrix
c=======================================================================
c     1. Preconditioning by use of irreducible representations:
c-----------------------------------------------------------------------
      if(lirredrep)then
         call logmsg('     Qirr into irreducible basis... ')
         call timeit('                                    ',dtime)
c-----------------------------------------------------------------------
c     forward transformation into irreducible basis: Qirr = P*Q*P^{-1}:
c-----------------------------------------------------------------------
         call Q2Qirr(dim,ne,neinv,Q)
         call check_Qirr(Q,ir_size,dim,maxNc,Nc,prec,lmax3)
         call timeit('     ...done with Qirr            : ',dtime)
c-----------------------------------------------------------------------
c     inversion of Qirr-matrix by LU decomposition:
c-----------------------------------------------------------------------
         call logmsg('     inversion of Qirr matrix...    ')
         call timeit(' ',dtime)                  
c-----------------------------------------------------------------------
c     loop over all irreducible representations (block matrices) of Qirr
c-----------------------------------------------------------------------
         lwork=5*dim
         allocate(Qb(lmax3,lmax3),stat=err)
         if(err.ne.0)then
            write(6,*)'FATAL ERROR: allocation failed for Q'
            stop
         endif

         do l=1,Nc
            if(ir_size(l).ne.0)then
               call cpy_Q2Qb(Q,Qb,dim,lmax3,maxNc,ir_size,l)
               bdim=ir_size(l)
               lwork=5*bdim
               call zgetrf(bdim,bdim,Qb,lmax3,ipiv,info)
               if (info.lt.0) then
                  call errmsg
     &           ("Illegal argument in LU-decomposition ZGETRF",.false.)
                  write(6,*) 'info = ',info
               else if (info.gt.0) then
                  call errmsg
     &            ("Singular value in LU-decomposition ZGETRF",.false.)
                  write(6,*) 'info = ',info
               end if
               call zgetri(bdim,Qb,lmax3,ipiv,work,lwork,info)
               call cpy_Qb2Q(Q,Qb,dim,lmax3,maxNc,ir_size,l)
            endif
         end do
         deallocate(Qb)
         call timeit('     ...done with Qirr-inversion        : ',dtime)
c-----------------------------------------------------------------------
c     backward transformation into reducible basis: 
c     Q^{-1} = P^{-1}*Q^{-1}*P:
c-----------------------------------------------------------------------
         call logmsg('     back transformation of Qirr-1: ')
         call timeit(' ',dtime)
         call Qirr2Q(dim,ne,neinv,Q)
         call timeit('     ...done with back-transform        : ',dtime)

c-----------------------------------------------------------------------
c     2. Use of perturbation theory:
c     Compute delta Q = Q-Q0, store in Q;
c     compute Q0^{-1}, store in Q0.
c-----------------------------------------------------------------------
      elseif(lprtb) then
         call logmsg('     Q-matrices for perturbation... ')
         call timeit(' ',dtime)
         call Qperturb(dim,Q,Q0) 
         call timeit('     ...done with Q-matrices for perturb: ',dtime)
c-----------------------------------------------------------------------
c     3. Direct inversion of the Q-matrix by LU decomposition:
c-----------------------------------------------------------------------
      else 
c-----------------------------------------------------------------------
c     inversion of Q-matrix:
c-----------------------------------------------------------------------
         call logmsg('     inversion of Q-matrix...       ')
         call timeit(' ',dtime)                  
         lwork=5*dim
         call zgetrf(dim,dim,Q,dim,ipiv,info)
         call zgetri(dim,Q,dim,ipiv,work,lwork,info)
         call timeit('     ...done with Q-inversion           : ',dtime)
        
      endif
c=======================================================================
c     calculation of T-matrix:
c=======================================================================
      call timeit(' ',dtime)
c-----------------------------------------------------------------------
c     Direct calculation:
c     T = -Rg Q * Q^(-1) 
c-----------------------------------------------------------------------
      if(.not.lprtb) then
         call logmsg('     T-matrix calculation (direct)..')
         call Tmatrix(mmax,nmmax,nmx,lmax,dim,lsi,DIHEDRAL,SIGMAh,
     &        czero,Q,TQsym,RgQ)
         deallocate(Q)
         allocate(Tmat2(2,lmax2,2,lmax2),stat=err)
         if(err.ne.0)then
            write(6,*)'FATAL ERROR: allocation failed for Tmat2'
            stop
         endif
         call cpy_Tmat(RgQ,Tmat2,lmax,lmax2,lsi,nmx,nmax,mmax,nmmax)
         deallocate(RgQ)
c-----------------------------------------------------------------------
c    Iterative calculation:
c     T^{0} = -Rg Q * Q_0^(-1)
c     T^{n} = -(Rg Q + T^{n-1}*\Delta Q) * Q_0^(-1)
c-----------------------------------------------------------------------
      else
         call logmsg('     T-matrix calculation (pertrb.).')
         allocate(Tmat(2,lmax,2,lmax),stat=err)
         if(err.ne.0)then
            write(6,*)'FATAL ERROR: allocation failed for RgQ'
            stop
         endif
         call Tperturb(mmax,nmmax,nmx,lmax,dim,lsi,DIHEDRAL,
     &        SIGMAh,czero,Q,Q0,TQsym,RgQ,Tmat,prtb_order)
         deallocate(RgQ)
         deallocate(Q)
         deallocate(Q0)
         allocate(Tmat2(2,lmax2,2,lmax2),stat=err)
         if(err.ne.0)then
            write(6,*)'FATAL ERROR: allocation failed for Tmat2'
            stop
         endif
         call cpy_Tmat(Tmat,Tmat2,lmax,lmax2,lsi,nmx,nmax,mmax,nmmax)
         deallocate(Tmat)
      endif
      call timeit('     ...done with T-matrix              : ',dtime)
      

c=======================================================================
      call logmsg('   ...done with subroutine NFM      ')

      return 
      end
      


C***********************************************************************
c     check if symmetry properties apply for Q or T-matrix elements: 
c
c     Author: Michael Kahnert
C***********************************************************************     
      subroutine sym(Nsym,SIGMAh,m,mm,indxsum,Cn,Sh,TQsym)
      
c-----------------------------------------------------------------------
c     CN-symmetry implies:
c     T_(m,n,m',n')=0 unless mod(m-m',N)=0
c     sigma_h-symmetry implies:
c     T_(m,n,m',n')^(k,k')=0 unless (n+n'+m+m'+k+k') even
c     (analogous relations hold for the Q-matrix)
c
c     Cn = .true. if we have CN-Symmetry .AND. mod(m-m',N).eq.0
c     Sh = 0      if we have no sigma_h symmetry; otherwise Sh is
c        =-1      if n+n'+m+m' odd, 
c        = 1      if n+n'+m+m' even.
c     TQsym= 3    if Cn=.true.
c          = 1    if Cn=.false. .AND. n+n'+m+m' even
c          =-1    if Cn=.false. .AND. n+n'+m+m' odd
c-----------------------------------------------------------------------
      implicit none
      
      integer Nsym,m,mm,indxsum, Sh,TQsym
      logical SIGMAh, Cn
      
      
c-----------------------------------------------------------------------
c     CN-symmetry:
c-----------------------------------------------------------------------           
      if (Nsym.eq.1) then
         Cn = .false.
      else
         if (mod(m-mm,Nsym).eq.0) then
            Cn = .false.
         else
            Cn = .true.
         end if
      end if
      
c-----------------------------------------------------------------------
c     sigma_h-symmetry:
c-----------------------------------------------------------------------  
      if (mod(indxsum,2).eq.0) then
	 Sh =1
      else
	 Sh =-1
      end if
      TQsym=Sh
      if (.not.SIGMAh) Sh=0
      
      if (Cn) TQsym=3
            
      return
      end
      

      
C***********************************************************************
c
c     The following routines for computing the matrices Q, Rg, and T 
c     have stongly evolved over the years and grown in complexity, due
c     to the use of symmetries, special quadrature schemes, perturbation
c     methods, etc. However, I want to mention that the first version
c     of this code has very much benefitted from the work by Heikki
c     Laitinen and Kari Lumme (JQSRT 60, 325-334, 1998). Their T-matrix 
c     code for general, star-shaped particles was very helpful for 
c     validating the first versions of Tsym. Wherever possible, I used 
c     the same notation for variables in the Q-matrix routines as in 
c     the code by Laitinen and Lumme. 
c
C***********************************************************************
c     Calculation of the matrices Q, Rg Q (using the special quadrature 
c     scheme for cylindrical particles):
C***********************************************************************
      subroutine Qmat_cyl(mmax,nmmax,nmx,lmax,lpmax,dim,
     & phipts,phiptsx,Nsym,DIHEDRAL,SIGMAh,ng1,ng2,ng1dim,ng2dim,
     & xth1,xth2,xphi,phic1,phic2,ngstart,r1,rth1,rphi1,r2,rth2,rphi2,
     & wth1,wth2,wphi,k,ks,mr,pi,nsi,lsi,lpsi,zero,cone,czero,Im,dtime,
     & Q,RgQ,TQsym)
      
      implicit none
      
c--------------------------------------------------------------------
c     Subroutine for calculating the Q- and RgQ-matrix of a 
c     polyhedral prism.
c     C A U T I O N :
c     This routine has been written for calculating the T-matrix
c     from the Q- and RgQ-matrix. Thus any constant factors appearing 
c     in front of both the Q- and the RgQ-matrix have been omitted,
c     since they cancel when forming T = -RgQ * Q^{-1}. The output
c     of this routine is therefore (depending on the particle's 
c     symmetries) only identical with the Q- and RgQ-matrices up
c     to a constant factor. This is important when comparing Q- and
c     RgQ-matrices computed with different methods, or when one
c     uses this code for calculating the internal field inside the
c     particle.
c--------------------------------------------------------------------
      
c-----------------------------------------------------------------------
c     input variables:
c-----------------------------------------------------------------------
      integer mmax,nmmax,lmax,lpmax,dim,phipts,phiptsx,Nsym,
     & ng1,ng2,ng1dim,ng2dim,ngstart,
     & nsi(nmmax),lsi(-mmax:mmax,nmmax),lpsi(0:mmax,nmmax),
     & nmx(-mmax:mmax)
      real*8 k,pi,zero,     
     & r1(ng1dim),rth1(ng1dim),rphi1(ng1dim),
     & r2(ng2dim,phipts),rth2(ng2dim,phipts),rphi2(ng2dim,phipts),
     & xth1(ng1dim),xth2(ng2dim,phipts),
     & wth1(ng1dim),wth2(ng2dim,phipts),wphi(phipts),
     & xphi(phipts),phic1(ng1dim),phic2(ng1dim)
      real*8 dtime
      complex*16 cone,czero,Im,ks,mr
      logical DIHEDRAL,SIGMAh

c-----------------------------------------------------------------------
c     output variables:
c-----------------------------------------------------------------------
      complex*16 Q(dim,dim),RgQ(2,lmax,2,lmax)
      integer TQsym(lmax,lmax)

c-----------------------------------------------------------------------
c     locally defined variables:
c-----------------------------------------------------------------------
      integer n,n1,nn,nn1,nindxsum,m,l,lp,lp1,lmx,mindxsum,ln,Mup,mm,
     &  indxsum,ll,llp,llp1,llmx,i,j,lnx,lln,llnx , mfac,mmfac ,Sh
      real*8 sth1(ng1dim),sth2(ng2dim,phiptsx),
     & a,fac,fac1,fac2,factmp,dm,dmm,mma,mmm,phihi,
     & tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,
     & P1(ng1dim,lpmax),dP1(ng1dim,lpmax),
     & P2(ng2dim,phiptsx,lpmax),dP2(ng2dim,phiptsx,lpmax),dfac2(nmmax)
      complex*16 b,cfac,mb,weight,
     & J11,J12,J21,J22,RgJ11,RgJ12,RgJ21,RgJ22,
     & Jt11,Jt12,Jt21,Jt22,RgJt11,RgJt12,RgJt21,RgJt22,
     & Jt12a,Jt21a,Jt22a,RgJt12a,RgJt21a,RgJt22a,
     & Btmp1,Btmp2,Btmp3,Btmp4,Btmp1a,Btmp2a,Btmp3a,Btmp4a,      
     & hn1(ng1dim,nmmax),wn1(ng1dim,nmmax),
     & jnn1(ng1dim,nmmax),
     & unn1(ng1dim,nmmax),
     & hn2(phiptsx,ng2dim,nmmax),wn2(phiptsx,ng2dim,nmmax),
     & jnn2(phiptsx,ng2dim,nmmax),
     & unn2(phiptsx,ng2dim,nmmax),eimphi(phiptsx,-2*mmax:2*mmax)
      logical Cn,lmsign
 
      if (DIHEDRAL) then
         phihi=pi/dble(Nsym)
      else
         phihi=2.d0*pi/dble(Nsym)
      end if
c-----------------------------------------------------------------------
c     calculation of Bessel-,Hankel-, and Wigner functions:
c-----------------------------------------------------------------------
      call logmsg('        Function calculations...    ')       
      call bessel(phipts,phiptsx,mmax,nmmax,
     &          ng1,ng2,ng1dim,ng2dim,ng2dim,r1,r2,k,ks,
     &          xth1,xth2,xphi,cone,Im,hn1,wn1,jnn1,unn1,
     &          hn2,wn2,jnn2,unn2,eimphi,sth1,sth2)
      call wigner(phipts,phiptsx,nmmax,mmax,lpsi,lpmax,
     &          ng1,ng2,ng1dim,ng2dim,ng2dim,xth1,xth2,P1,dP1,P2,dP2)

c-----------------------------------------------------------------------
c     preliminary settings
c-----------------------------------------------------------------------
      do n=1,nmmax
         dfac2(n) = dsqrt(dble(2*n+1)/dble(nsi(n)))
      end do
c-----------------------------------------------------------------------
c     calculation of Q and Rg Q:
c-----------------------------------------------------------------------
      
      if (DIHEDRAL) then
         Mup=0
      else
         Mup=mmax
      end if

      do m = -mmax,mmax
         if (m.lt.0) then
            mfac=(-1)**m
         else
            mfac=1
         end if
         dm = dble(m)
         do mm=-mmax,Mup
            if (mm.lt.0) then
               mmfac=(-1)**mm
            else
               mmfac=1
            end if
            if (mfac*mmfac.lt.0) then
               lmsign=.true.
            else
               lmsign=.false.
            end if
            dmm = dble(mm)
            mmm = dm*dmm
            do n=max(1,abs(m)),nmx(m)
               n1 = nsi(n)
               a = dble(n1)/k
               mma = dmm*a
               fac2 = dfac2(n) ! /(4.d0*pi)
               lp1 = lpsi(0,n)
               lp=lpsi(abs(m),n)
               do nn = max(1,abs(mm)),nmx(mm)
                  nn1=nsi(nn)
                  b = dcmplx(dble(nn1),zero)/ks
                  nindxsum=n+nn
                  factmp=fac2*dfac2(nn)
                  llp1 = lpsi(0,nn)
                  l=lsi(m,n)
                  lmx = l+lmax
                  mb = dcmplx(dm,zero)*b
                  mindxsum=nindxsum+m
                  llp=lpsi(abs(mm),nn)
	          indxsum=mindxsum+mm
	          ll=lsi(mm,nn)
		  llmx = ll+lmax
		  call sym(Nsym,SIGMAh,m,mm,indxsum,
     &                     Cn,Sh,TQsym(l,ll))
		  if (Cn) then
                     Q(l,ll)=czero
                     Q(lmx,ll)=czero
                     Q(l,llmx)=czero
                     Q(lmx,llmx)=czero
                     RgQ(1,l,1,ll)=czero
                     RgQ(2,l,1,ll)=czero
                     RgQ(1,l,2,ll)=czero
                     RgQ(2,l,2,ll)=czero
		     goto 19		     
		  end if
c-----------------------------------------------------------------------
c                 Integration over top (and bottom) facets:
c-----------------------------------------------------------------------                 
                  J11=czero
                  J12=czero
                  J21=czero
                  J22=czero
                  RgJ11=czero
                  RgJ12=czero
                  RgJ21=czero
                  RgJ22=czero
		  do i=1,ng1dim
		     if (m.ne.mm.and.
     &                  (i.lt.ngstart.or.i.gt.(2*ng1-ngstart+1))) 
     &                  goto 10
                     if (i.lt.ngstart.or.i.gt.(2*ng1-ngstart+1)) then  
                        ! inner circle on top or bottom facet:
                        weight = phihi   ! upper phi-interval boundary
                     else ! outside inner circle on top facet
                        if (m.eq.mm) then
                           if (DIHEDRAL) then
                              weight = phic1(i) 
                           else
                              weight = phic1(i) + phic2(i)
                           end if
                        else
                           if (DIHEDRAL) then  
                              weight = 1.d0/dble(mm-m)*
     &                           dsin((mm-m)*phic1(i))
                           else
                              weight = Im/dble(mm-m)*(exp(-Im*(mm-m)*
     &                           phic2(i))-
     &                           exp(Im*(mm-m)*phic1(i)))
                           end if
                       end if
                     end if
                     if (lmsign) then
		        tmp1 =- P1(i,lp)* P1(i,llp)
			tmp2 =- P1(i,lp)*dP1(i,llp)
			tmp3 =-dP1(i,lp)* P1(i,llp)
                     else
		        tmp1 =  P1(i,lp)* P1(i,llp)
			tmp2 =  P1(i,lp)*dP1(i,llp)
			tmp3 = dP1(i,lp)* P1(i,llp)
                     end if
                     if (Sh.ge.0) then
			if (lmsign) then
			   tmp4 =-dP1(i,lp)*dP1(i,llp)
			else
			   tmp4 = dP1(i,lp)*dP1(i,llp)
			end if
			tmp6 = r1(i)*(mmm*tmp1/sth1(i)+
     &                      tmp4*sth1(i))
                     end if
                     if (Sh.le.0) then			
			tmp5 = r1(i)*(dm*tmp2+dmm*tmp3)
			tmp8 = rth1(i)*tmp1
			Btmp4  = wn1(i,n)*unn1(i,nn)
			Btmp4a = dble(wn1(i,n))*unn1(i,nn)
                     end if
                     Btmp1  = hn1(i,n)*jnn1(i,nn)
                     Btmp2  = hn1(i,n)*unn1(i,nn)
                     Btmp3  = wn1(i,n)*jnn1(i,nn)
                     Btmp1a = dble(hn1(i,n))*jnn1(i,nn)
                     Btmp2a = dble(hn1(i,n))*unn1(i,nn)
                     Btmp3a = dble(wn1(i,n))*jnn1(i,nn)
                     weight = wth1(i)*weight
                     if (Sh.ge.0) then
			 J12 = J12 + weight*(Btmp3*tmp6+a*Btmp1*
     &                        rth1(i)*tmp2*sth1(i))
                         J21 = J21 + weight*
     &                         (b*Btmp1*rth1(i)*tmp3*sth1(i)
     &                         +Btmp2*tmp6)
			 RgJ12 = RgJ12 + weight*
     &                         (Btmp3a*tmp6+a*Btmp1a*
     &                        rth1(i)*tmp2*sth1(i))
                         RgJ21 = RgJ21 + weight*(b*Btmp1a*
     &                         rth1(i)*tmp3*sth1(i)
     &                         +Btmp2a*tmp6)
                     end if
                     if (Sh.le.0) then
			 J11 = J11 + weight*Btmp1*tmp5
                         J22 = J22 + weight*(Im*Btmp4*tmp5+b*
     &                         Btmp3*Im*dm*tmp8 +
     &                        a*Btmp2*Im*dmm*tmp8)
			 RgJ11 = RgJ11 + weight*Btmp1a*tmp5
                         RgJ22 = RgJ22 + weight*
     &                         (Im*Btmp4a*tmp5+b*Btmp3a*
     &                         Im*dm*tmp8 +
     &                         a*Btmp2a*Im*dmm*tmp8)
                     end if
 10                  continue
		  end do ! i=1,ng1dim
c-----------------------------------------------------------------------
c                 Integration over side facets:
c-----------------------------------------------------------------------                
		  do j=1,phiptsx
                     Jt11=czero
                     Jt12=czero
                     Jt21=czero
                     Jt22=czero
                     RgJt11=czero
                     RgJt12=czero
                     RgJt21=czero
                     RgJt22=czero
                     if (DIHEDRAL) then
                        Jt12a=czero
                        Jt21a=czero
                        Jt22a=czero
                        RgJt12a=czero
                        RgJt21a=czero
                        RgJt22a=czero                     
                     end if
		     do i=1,ng2dim
			if (lmsign) then
		           tmp1 = -P2(i,j,lp)* P2(i,j,llp)
			   tmp2 = -P2(i,j,lp)*dP2(i,j,llp)
			   tmp3 = -dP2(i,j,lp)* P2(i,j,llp)
			else
		           tmp1 =  P2(i,j,lp)* P2(i,j,llp)
			   tmp2 =  P2(i,j,lp)*dP2(i,j,llp)
			   tmp3 = dP2(i,j,lp)* P2(i,j,llp)
			end if
			if (Sh.ge.0) then
			   if (lmsign) then
			      tmp4 =-dP2(i,j,lp)*dP2(i,j,llp)
			   else
			      tmp4 = dP2(i,j,lp)*dP2(i,j,llp)
			   end if
			   tmp6 = r2(i,j)*(mmm*tmp1/sth2(i,j)+
     &                         tmp4*sth2(i,j))
			   tmp7 = rphi2(i,j)*tmp1/sth2(i,j)
			end if
			if (Sh.le.0) then			
			   tmp5 = r2(i,j)*(dm*tmp2+dmm*tmp3)
			   tmp8 = rth2(i,j)*tmp1
			   Btmp4  =wn2(j,i,n)*unn2(j,i,nn)
			   Btmp4a =dble(wn2(j,i,n))*unn2(j,i,nn)
			end if
			Btmp1  = hn2(j,i,n)*jnn2(j,i,nn)
			Btmp2  = hn2(j,i,n)*unn2(j,i,nn)
			Btmp3  = wn2(j,i,n)*jnn2(j,i,nn)
			Btmp1a = dble(hn2(j,i,n))*jnn2(j,i,nn)
			Btmp2a = dble(hn2(j,i,n))*unn2(j,i,nn)
			Btmp3a = dble(wn2(j,i,n))*jnn2(j,i,nn)
                        weight = wth2(i,j)		
			if (.not.DIHEDRAL) then
                          if (Sh.ge.0) then
			    Jt12 = Jt12 + weight*(Btmp3*tmp6+a*Btmp1*
     &                           (rth2(i,j)*tmp2*sth2(i,j)
     &                           +Im*dmm*tmp7))
                            Jt21 = Jt21 + weight*
     &                            (b*Btmp1*(rth2(i,j)
     &                            *tmp3*sth2(i,j)
     &                            -Im*dm*tmp7)+Btmp2*tmp6)
			    RgJt12 = RgJt12 + weight*
     &                            (Btmp3a*tmp6+a*Btmp1a*
     &                           (rth2(i,j)*tmp2*sth2(i,j)
     &                           +Im*dmm*tmp7))
                            RgJt21 = RgJt21 + weight*(b*Btmp1a*
     &                            (rth2(i,j)*tmp3*sth2(i,j)
     &                            -Im*dm*tmp7)+Btmp2a*tmp6)
                          end if
			  if (Sh.le.0) then
			    Jt11 = Jt11 + weight*Btmp1*tmp5
                            Jt22 = Jt22 + weight*(Im*Btmp4*tmp5+b*
     &                            Btmp3*(Im*dm*tmp8+rphi2(i,j)
     &                            *tmp3) + a*Btmp2
     &                            *(Im*dmm*tmp8-rphi2(i,j)*tmp2))
			    RgJt11 = RgJt11 + weight*Btmp1a*tmp5
                            RgJt22 = RgJt22 + weight*
     &                            (Im*Btmp4a*tmp5+b*Btmp3a*
     &                            (Im*dm*tmp8+rphi2(i,j)*tmp3) +
     &                            a*Btmp2a*(Im*dmm*tmp8-
     &                            rphi2(i,j)*tmp2))
                          end if
                        else  ! if (DIHEDRAL)
                          if (Sh.ge.0) then
			    Jt12 = Jt12 + weight*(Btmp3*tmp6+a*
     &                            Btmp1*rth2(i,j)*tmp2
     &                            *sth2(i,j))
			    Jt12a = Jt12a + weight*a*Btmp1*dmm*tmp7
                            Jt21 = Jt21 + weight*
     &                            (b*Btmp1*rth2(i,j)
     &                            *tmp3*sth2(i,j)
     &                            +Btmp2*tmp6)
                            Jt21a = Jt21a - weight*b*Btmp1*dm*tmp7
			    RgJt12 = RgJt12 + weight*
     &                            (Btmp3a*tmp6+a*Btmp1a*rth2(i,j)*
     &                            tmp2*sth2(i,j))
			    RgJt12a = RgJt12a + weight*a*Btmp1a*
     &                            dmm*tmp7
                            RgJt21 = RgJt21 + weight*(b*Btmp1a*
     &                            rth2(i,j)*tmp3*sth2(i,j)
     &                            +Btmp2a*tmp6)
                            RgJt21a = RgJt21a - weight*b*Btmp1a*
     &                            dm*tmp7
                          end if
			  if (Sh.le.0) then
			    Jt11 = Jt11 + weight*Btmp1*tmp5
                            Jt22 = Jt22 + weight*(Btmp4*tmp5+
     &                            tmp8*(b*Btmp3*dm +
     &                            a*Btmp2*dmm))
                            Jt22a = Jt22a + weight*rphi2(i,j)*
     &                            (b*Btmp3*tmp3 - a*Btmp2*tmp2)
			    RgJt11 = RgJt11 + weight*Btmp1a*tmp5
                            RgJt22 = RgJt22 + weight*
     &                            (Btmp4a*tmp5 + tmp8*(b*Btmp3a*dm +
     &                            a*Btmp2a*dmm))
                            RgJt22a = RgJt22a + weight*rphi2(i,j)*
     &                            (b*Btmp3a*tmp3 - a*Btmp2a*tmp2)
                          end if                        
                        end if
		     end do ! i=1,ng2dim
		     weight = wphi(j)*eimphi(j,mm-m)
		     if (.not.DIHEDRAL) then
		        if (Sh.ge.0) then
		           J12 = J12 + weight*Jt12
		           J21 = J21 + weight*Jt21
		           RgJ12 = RgJ12 + weight*RgJt12
		           RgJ21 = RgJ21 + weight*RgJt21
		        end if
		        if (Sh.le.0) then	     
		           J11 = J11 + weight*Jt11
		           J22 = J22 + weight*Jt22
		           RgJ11 = RgJ11 + weight*RgJt11
		           RgJ22 = RgJ22 + weight*RgJt22
		        end if
		     else  ! if (DIHEDRAL)	     
		        if (Sh.ge.0) then
		           J12 = J12 + dreal(weight)*Jt12
     &                           - dimag(weight)*Jt12a
		           J21 = J21 + dreal(weight)*Jt21
     &                           - dimag(weight)*Jt21a
		           RgJ12 = RgJ12 + dreal(weight)*RgJt12
     &                           - dimag(weight)*RgJt12a
		           RgJ21 = RgJ21 + dreal(weight)*RgJt21
     &                           - dimag(weight)*RgJt21a
		        end if
		        if (Sh.le.0) then	     
		           J11 = J11 + dreal(weight)*Jt11
		           J22 = J22 + Im*(dreal(weight)*Jt22
     &                           +dimag(weight)*Jt22a)
		           RgJ11 = RgJ11 + dreal(weight)*RgJt11
		           RgJ22 = RgJ22 + Im*(dreal(weight)*RgJt22
     &                           +dimag(weight)*RgJt22a)
		        end if
		     end if
		  end do ! j=1,phiptsx
c-----------------------------------------------------------------------
c  factor (-1)^(m+mm) * sqrt((2n+1)(2n'+1)/(n(n+1)n'(n'+1)))/(4*pi)
c  in front of the integrals, where the constant 1/(4*pi) can be 
c  ommited, since Q^{-1} will get a factor 4*pi such that 
c  T = -Rg Q * Q^{-1} will be the same:
c-----------------------------------------------------------------------
                  fac=factmp
                  if (mod(m+mm,2).ne.0) fac = -fac
                  cfac=-Im*fac
c-----------------------------------------------------------------------
		  if (Sh.ge.0) then
                     J12= fac*J12
                     J21=-fac*J21
                     RgJ12= fac*RgJ12
                     RgJ21=-fac*RgJ21
                     Q(l,ll)=J12+mr*J21
                     Q(lmx,llmx)=J21+mr*J12
                     RgQ(1,l,1,ll)=RgJ12+mr*RgJ21
                     RgQ(2,l,2,ll)=RgJ21+mr*RgJ12
		  end if
                  if (Sh.le.0) then
                     J11=cfac*J11
                     J22=-fac*J22
                     RgJ11=cfac*RgJ11
                     RgJ22=-fac*RgJ22
                     Q(lmx,ll)=J11+mr*J22
                     Q(l,llmx)=J22+mr*J11
                     RgQ(2,l,1,ll)=RgJ11+mr*RgJ22
                     RgQ(1,l,2,ll)=RgJ22+mr*RgJ11
		  end if
		  if (Sh.eq.1) then
                     Q(lmx,ll)=czero
                     Q(l,llmx)=czero
                     RgQ(2,l,1,ll)=czero
                     RgQ(1,l,2,ll)=czero
		  else if (Sh.eq.-1) then
                     Q(l,ll)=czero
                     Q(lmx,llmx)=czero
                     RgQ(1,l,1,ll)=czero
                     RgQ(2,l,2,ll)=czero		  
		  end if
 19               continue
               end do           ! nn = max(abs(mm),1),nmx(mm)
	    end do              ! n  = max(abs(m),1),nmx(m)
         end do                 ! mm = -mmax,Mup
      end do                    ! m  = -mmax,mmax
c-----------------------------------------------------------------------
c      The following relations are a consequence of the dihedral
c      symmetry:  
c      Q_(n,-m,n',-m')^(k,k') = 
c              (-1)^(m+m'+k+k') * Q_(n,m,n',m')^(k,k')    
c      CAUTION: The dihedral symmetry relation in this form assumes that
c      the coordinate system is chosen such that one of the dihedral 
c      axis coincides with the x-axis. This is important when 
c      implementing new geometries with dihedral symmetries!
c-----------------------------------------------------------------------  
      if (DIHEDRAL) then
         do m = -mmax,mmax
            do mm = 1,mmax
               do n=max(1,abs(m)),nmx(m)
                  do nn = max(1,mm), nmx(mm)
                     nindxsum=n+nn
                     l=lsi(m,n)
                     lmx = l+lmax
                     ln=lsi(-m,n)
                     lnx= ln+lmax    
                     mindxsum=m+mm
                     indxsum=nindxsum+mindxsum
                     if (mod(mindxsum,2).eq.0) then
                        fac1=1.d0
                     else
                        fac1=-1.d0
                     end if
                     ll=lsi(mm,nn)
                     llmx = ll+lmax
		     call sym(Nsym,SIGMAh,m,mm,indxsum,
     &                     Cn,Sh,TQsym(l,ll))
		     if (Cn) then
                        Q(l,ll)=         czero
                        Q(lmx,llmx)=     czero
                        RgQ(1,l,1,ll)=   czero
                        RgQ(2,l,2,ll)=   czero
                        Q(lmx,ll)=       czero
                        Q(l,llmx)=       czero
                        RgQ(2,l,1,ll)=   czero
                        RgQ(1,l,2,ll)=   czero		     
		     else
                        lln = lsi(-mm,nn)
			llnx= lln+lmax
                        Q(l,ll)=         fac1*Q(ln,lln)
                        Q(lmx,llmx)=     fac1*Q(lnx,llnx)
                        RgQ(1,l,1,ll)=   fac1*RgQ(1,ln,1,lln)
                        RgQ(2,l,2,ll)=   fac1*RgQ(2,ln,2,lln)
                        Q(lmx,ll)=      -fac1*Q(lnx,lln)
                        Q(l,llmx)=      -fac1*Q(ln,llnx)
                        RgQ(2,l,1,ll)=  -fac1*RgQ(2,ln,1,lln)
                        RgQ(1,l,2,ll)=  -fac1*RgQ(1,ln,2,lln)
		     end if
	          end do        ! nn = max(1,abs(mm),nmx(mm)
	       end do           ! n  = max(1,abs(m)),nmx(m)
	    end do              ! mm = 1,mmax
         end do                 ! m  = -mmax,mmax  
      end if	  

      return
      end
      
     
      
C***********************************************************************
c     calculation of the matrices Q, Rg Q (using the double-quadrature 
c     scheme):
C***********************************************************************

      subroutine Qmat(mmax,nmmax,nmx,lmax,lpmax,dim,thpts,thptsx,
     & phipts,phiptsx,Nsym,DIHEDRAL,SIGMAh,ng1,ng2,ng1dim,ng2dim,
     & xth1,xth2,xphi,r1,rth1,rphi1,r2,rth2,rphi2,
     & wth1,wth2,wphi,k,ks,mr,pi,nsi,lsi,lpsi,zero,cone,czero,Im,dtime,
     & Q,RgQ,TQsym)
      
      implicit none
      
c--------------------------------------------------------------------
c     Subroutine for calculating the Q- and RgQ-matrix.
c     C A U T I O N :
c     This routine has been written for calculating the T-matrix
c     from the Q- and RgQ-matrix. Thus any constant factors appearing 
c     in front of both the Q- and the RgQ-matrix have been omitted,
c     since they cancel when forming T = -RgQ * Q^{-1}. The output
c     of this routine is therefore (depending on the particle's 
c     symmetries) only identical with the Q- and RgQ-matrices up
c     to a constant factor. This is important when comparing Q- and
c     RgQ-matrices computed with different methods, or when one
c     uses this code for calculating the internal field inside the
c     particle.
c--------------------------------------------------------------------
      
c--------------------------------------------------------------------
c     input variables:
c--------------------------------------------------------------------
      integer mmax,nmmax,lmax,lpmax,dim,thpts,thptsx,phipts,
     &        phiptsx,Nsym,ng1dim,ng2dim,ng1,ng2,
     & nsi(nmmax),lsi(-mmax:mmax,nmmax),lpsi(0:mmax,nmmax),
     & nmx(-mmax:mmax)
      real*8 k,pi,zero,     
     & r1(ng1dim),rth1(ng1dim),rphi1(ng1dim),
     & r2(ng2dim,phipts),rth2(ng2dim,phipts),rphi2(ng2dim,phipts),
     & xth1(ng1dim),xth2(ng2dim,phipts),
     & wth1(ng1dim),wth2(ng2dim,phipts),
     & wphi(phipts),xphi(phipts),dtime
      complex*16 cone,czero,Im,ks,mr
      logical DIHEDRAL,SIGMAh

c--------------------------------------------------------------------
c     output variables:
c--------------------------------------------------------------------
      complex*16 Q(dim,dim),RgQ(2,lmax,2,lmax)
      integer TQsym(lmax,lmax)

c--------------------------------------------------------------------
c     locally defined variables:
c--------------------------------------------------------------------
      integer n,n1,nn,nn1,nindxsum,m,l,lp,lmx,mindxsum,ln,Mup,mm,
     &  indxsum,ll,llp,llmx,i,j,lnx,lln,llnx , mfac,mmfac,Sh
      real*8 sth1(ng1dim),sth2(thptsx,phiptsx),
     & a,fac,fac1,fac2,factmp,dm,dmm,mma,mmm,
     & tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,
     & P1(ng1dim,lpmax),dP1(ng1dim,lpmax),
     & P2(thptsx,phiptsx,lpmax),dP2(thptsx,phiptsx,lpmax),
     & dfac2(nmmax)  
      complex*16 b,cfac,mb,weight,
     & J11,J12,J21,J22,RgJ11,RgJ12,RgJ21,RgJ22,
     & Jt11,Jt12,Jt21,Jt22,RgJt11,RgJt12,RgJt21,RgJt22,
     & Jt12a,Jt21a,Jt22a,RgJt12a,RgJt21a,RgJt22a,
     & Btmp1,Btmp2,Btmp3,Btmp4,Btmp1a,Btmp2a,Btmp3a,Btmp4a,      
     & hn1(ng1dim,nmmax),wn1(ng1dim,nmmax),
     & jnn1(ng1dim,nmmax),
     & unn1(ng1dim,nmmax),
     & hn2(phiptsx,thptsx,nmmax),wn2(phiptsx,thptsx,nmmax),
     & jnn2(phiptsx,thptsx,nmmax),
     & unn2(phiptsx,thptsx,nmmax),eimphi(phiptsx,-2*mmax:2*mmax)
      logical Cn,lmsign
      
c-----------------------------------------------------------------------
c     preliminary settings
c-----------------------------------------------------------------------
      do n=1,nmmax
         dfac2(n) = dsqrt(dble(2*n+1)/dble(nsi(n)))
      end do
c-----------------------------------------------------------------------
c     calculation of Bessel-,Hankel-, and Wigner functions:
c-----------------------------------------------------------------------
      call logmsg('        Function calculations...    ')       
      call bessel(phipts,phiptsx,mmax,nmmax,
     &          ng1,ng2,ng1dim,ng2dim,thptsx,r1,r2,k,ks,
     &          xth1,xth2,xphi,cone,Im,hn1,wn1,jnn1,unn1,
     &          hn2,wn2,jnn2,unn2,eimphi,sth1,sth2)

      call wigner(phipts,phiptsx,nmmax,mmax,lpsi,lpmax,
     &          ng1,ng2,ng1dim,ng2dim,thptsx,xth1,xth2,P1,dP1,P2,dP2)

c-----------------------------------------------------------------------
c     calculation of Q and Rg Q:
c-----------------------------------------------------------------------
      
      if (DIHEDRAL) then
         Mup=0
      else
         Mup=mmax
      end if

      do m = -mmax,mmax
         if (m.lt.0) then
            mfac=(-1)**m
         else
            mfac=1
         end if
         dm = dble(m)
         do mm=-mmax,Mup
            if (mm.lt.0) then
               mmfac=(-1)**mm
            else
               mmfac=1
            end if
            if (mfac*mmfac.lt.0) then
               lmsign=.true.
            else
               lmsign=.false.
            end if
            dmm = dble(mm)
            mmm = dm*dmm
            do n=max(1,abs(m)),nmx(m)
               n1 = nsi(n)
               a = dble(n1)/k
               mma = dmm*a
               fac2 = dfac2(n) ! /(4.d0*pi)
               lp=lpsi(abs(m),n)
               do nn = max(1,abs(mm)),nmx(mm)
                  nn1=nsi(nn)
                  b = dcmplx(dble(nn1),zero)/ks
                  nindxsum=n+nn
                  factmp=fac2*dfac2(nn)
                  l=lsi(m,n)
                  lmx = l+lmax
                  mb = dcmplx(dm,zero)*b
                  mindxsum=nindxsum+m
                  llp=lpsi(abs(mm),nn)
	          indxsum=mindxsum+mm
	          ll=lsi(mm,nn)
		  llmx = ll+lmax
		  call sym(Nsym,SIGMAh,m,mm,indxsum,
     &                     Cn,Sh,TQsym(l,ll))
		  if (Cn) then
                     Q(l,ll)=czero
                     Q(lmx,ll)=czero
                     Q(l,llmx)=czero
                     Q(lmx,llmx)=czero
                     RgQ(1,l,1,ll)=czero
                     RgQ(2,l,1,ll)=czero
                     RgQ(1,l,2,ll)=czero
                     RgQ(2,l,2,ll)=czero
		     goto 19		     
		  end if
                  J11=czero
                  J12=czero
                  J21=czero
                  J22=czero
                  RgJ11=czero
                  RgJ12=czero
                  RgJ21=czero
                  RgJ22=czero
		  do j=1,phiptsx
                     Jt11=czero
                     Jt12=czero
                     Jt21=czero
                     Jt22=czero
                     RgJt11=czero
                     RgJt12=czero
                     RgJt21=czero
                     RgJt22=czero
                     if (DIHEDRAL) then
                        Jt12a=czero
                        Jt21a=czero
                        Jt22a=czero
                        RgJt12a=czero
                        RgJt21a=czero
                        RgJt22a=czero                     
                     end if
		     do i=1,thptsx
			if (lmsign) then
		           tmp1 = -P2(i,j,lp) * P2(i,j,llp)
			   tmp2 = -P2(i,j,lp) *dP2(i,j,llp)
			   tmp3 = -dP2(i,j,lp)* P2(i,j,llp)
			else
		           tmp1 =  P2(i,j,lp)* P2(i,j,llp)
			   tmp2 =  P2(i,j,lp)*dP2(i,j,llp)
			   tmp3 = dP2(i,j,lp)* P2(i,j,llp)
			end if
			if (Sh.ge.0) then
			   if (lmsign) then
			      tmp4 =-dP2(i,j,lp)*dP2(i,j,llp)
			   else
			      tmp4 = dP2(i,j,lp)*dP2(i,j,llp)
			   end if
			   tmp6 = r2(i,j)*(mmm*tmp1/sth2(i,j)+
     &                         tmp4*sth2(i,j))
			   tmp7 = rphi2(i,j)*tmp1/sth2(i,j)
			end if
			if (Sh.le.0) then			
			   tmp5 = r2(i,j)*(dm*tmp2+dmm*tmp3)
			   tmp8 = rth2(i,j)*tmp1
			   Btmp4  = wn2(j,i,n)*unn2(j,i,nn)
			   Btmp4a = dble(wn2(j,i,n))*unn2(j,i,nn)
			end if
			Btmp1  = hn2(j,i,n)*jnn2(j,i,nn)
			Btmp2  = hn2(j,i,n)*unn2(j,i,nn)
			Btmp3  = wn2(j,i,n)*jnn2(j,i,nn)
			Btmp1a = dble(hn2(j,i,n))*jnn2(j,i,nn)
			Btmp2a = dble(hn2(j,i,n))*unn2(j,i,nn)
			Btmp3a = dble(wn2(j,i,n))*jnn2(j,i,nn)
			if (.not.DIHEDRAL) then
                          if (Sh.ge.0) then
			    Jt12 = Jt12 + wth2(i,j)*(Btmp3*tmp6+a*Btmp1*
     &                           (rth2(i,j)*tmp2*sth2(i,j)+Im*dmm*tmp7))
                            Jt21 = Jt21 + wth2(i,j)*
     &                           (b*Btmp1*(rth2(i,j)*tmp3*sth2(i,j)
     &                           -Im*dm*tmp7)+Btmp2*tmp6)
			    RgJt12 = RgJt12 + wth2(i,j)*
     &                           (Btmp3a*tmp6+a*Btmp1a*
     &                           (rth2(i,j)*tmp2*sth2(i,j)+Im*dmm*tmp7))
                            RgJt21 = RgJt21 + wth2(i,j)*(b*Btmp1a*
     &                           (rth2(i,j)*tmp3*sth2(i,j)
     &                           -Im*dm*tmp7)+Btmp2a*tmp6)
                          end if
			  if (Sh.le.0) then
			    Jt11 = Jt11 + wth2(i,j)*Btmp1*tmp5
                            Jt22 = Jt22 + wth2(i,j)*(Im*Btmp4*tmp5+b*
     &                            Btmp3*(Im*dm*tmp8+rphi2(i,j)*tmp3) +
     &                            a*Btmp2*(Im*dmm*tmp8-rphi2(i,j)*tmp2))
			    RgJt11 = RgJt11 + wth2(i,j)*Btmp1a*tmp5
                            RgJt22 = RgJt22 + wth2(i,j)*
     &                            (Im*Btmp4a*tmp5+b*Btmp3a*
     &                            (Im*dm*tmp8+rphi2(i,j)*tmp3) +
     &                            a*Btmp2a*(Im*dmm*tmp8-
     &                            rphi2(i,j)*tmp2))
                          end if
                        else  ! if (DIHEDRAL)
                          if (Sh.ge.0) then
			    Jt12 = Jt12 + wth2(i,j)*(Btmp3*tmp6+a*Btmp1*
     &                            rth2(i,j)*tmp2*sth2(i,j))
			    Jt12a = Jt12a + wth2(i,j)*a*Btmp1*dmm*tmp7
                            Jt21 = Jt21 + wth2(i,j)*
     &                            (b*Btmp1*rth2(i,j)*tmp3*sth2(i,j)
     &                            +Btmp2*tmp6)
                            Jt21a = Jt21a - wth2(i,j)*b*Btmp1*dm*tmp7
			    RgJt12 = RgJt12 + wth2(i,j)*
     &                            (Btmp3a*tmp6+a*Btmp1a*rth2(i,j)*
     &                            tmp2*sth2(i,j))
			    RgJt12a = RgJt12a + wth2(i,j)*a*Btmp1a*
     &                            dmm*tmp7
                            RgJt21 = RgJt21 + wth2(i,j)*(b*Btmp1a*
     &                            rth2(i,j)*tmp3*sth2(i,j)
     &                            +Btmp2a*tmp6)
                            RgJt21a = RgJt21a - wth2(i,j)*b*Btmp1a*
     &                            dm*tmp7
                          end if
			  if (Sh.le.0) then
			    Jt11 = Jt11 + wth2(i,j)*Btmp1*tmp5
                            Jt22 = Jt22 + wth2(i,j)*(Btmp4*tmp5+
     &                            tmp8*(b*Btmp3*dm +
     &                            a*Btmp2*dmm))
                            Jt22a = Jt22a + wth2(i,j)*rphi2(i,j)*
     &                            (b*Btmp3*tmp3 - a*Btmp2*tmp2)
			    RgJt11 = RgJt11 + wth2(i,j)*Btmp1a*tmp5
                            RgJt22 = RgJt22 + wth2(i,j)*
     &                            (Btmp4a*tmp5 + tmp8*(b*Btmp3a*dm +
     &                            a*Btmp2a*dmm))
                            RgJt22a = RgJt22a + wth2(i,j)*rphi2(i,j)*
     &                            (b*Btmp3a*tmp3 - a*Btmp2a*tmp2)
                          end if                        
                        end if
		     end do ! i=1,thptsx
c-----------------------------------------------------------------------
		     weight = wphi(j)*eimphi(j,mm-m)
		     if (.not.DIHEDRAL) then
		        if (Sh.ge.0) then
		           J12 = J12 + weight*Jt12
		           J21 = J21 + weight*Jt21
		           RgJ12 = RgJ12 + weight*RgJt12
		           RgJ21 = RgJ21 + weight*RgJt21
		        end if
		        if (Sh.le.0) then	     
		           J11 = J11 + weight*Jt11
		           J22 = J22 + weight*Jt22
		           RgJ11 = RgJ11 + weight*RgJt11
		           RgJ22 = RgJ22 + weight*RgJt22
		        end if
		     else  ! if (DIHEDRAL)	     
		        if (Sh.ge.0) then
		           J12 = J12 + dreal(weight)*Jt12
     &                           - dimag(weight)*Jt12a
		           J21 = J21 + dreal(weight)*Jt21
     &                           - dimag(weight)*Jt21a
		           RgJ12 = RgJ12 + dreal(weight)*RgJt12
     &                           - dimag(weight)*RgJt12a
		           RgJ21 = RgJ21 + dreal(weight)*RgJt21
     &                           - dimag(weight)*RgJt21a
		        end if
		        if (Sh.le.0) then	     
		           J11 = J11 + dreal(weight)*Jt11
		           J22 = J22 + Im*(dreal(weight)*Jt22
     &                           +dimag(weight)*Jt22a)
c		           J22 = J22 + weight*Jt22
		           RgJ11 = RgJ11 + dreal(weight)*RgJt11
c		           RgJ22 = RgJ22 + weight*RgJt22
		           RgJ22 = RgJ22 + Im*(dreal(weight)*RgJt22
     &                           +dimag(weight)*RgJt22a)
		        end if
		     end if
		  end do ! j=1,phiptsx
 10               continue
c-----------------------------------------------------------------------
c  factor (-1)^(m+mm) * sqrt((2n+1)(2n'+1)/(n(n+1)n'(n'+1)))/(4*pi)
c  in front of the integrals, where the constant 1/(4*pi) can be 
c  ommited, since Q^{-1} will get a factor 4*pi such that 
c  T = -Rg Q * Q^{-1} will be the same:
c-----------------------------------------------------------------------
                  fac=factmp
                  if (mod(m+mm,2).ne.0) fac = -fac
                  cfac=-Im*fac
c-----------------------------------------------------------------------

		  if (Sh.ge.0) then
                     J12= fac*J12
                     J21=-fac*J21
                     RgJ12= fac*RgJ12
                     RgJ21=-fac*RgJ21
                     Q(l,ll)       = J12+mr*J21
                     Q(lmx,llmx)   = J21+mr*J12
                     RgQ(1,l,1,ll) = RgJ12+mr*RgJ21
                     RgQ(2,l,2,ll) = RgJ21+mr*RgJ12
		  end if
                  if (Sh.le.0) then
                     J11=cfac*J11
                     J22=-fac*J22
                     RgJ11=cfac*RgJ11
                     RgJ22=-fac*RgJ22
                     Q(lmx,ll)     = J11+mr*J22
                     Q(l,llmx)     = J22+mr*J11
                     RgQ(2,l,1,ll) = RgJ11+mr*RgJ22
                     RgQ(1,l,2,ll) = RgJ22+mr*RgJ11
		  end if
		  if (Sh.eq.1) then
                     Q(lmx,ll)     = czero
                     Q(l,llmx)     = czero
                     RgQ(2,l,1,ll) = czero
                     RgQ(1,l,2,ll) = czero
		  else if (Sh.eq.-1) then
                     Q(l,ll)       = czero
                     Q(lmx,llmx)   = czero
                     RgQ(1,l,1,ll) = czero
                     RgQ(2,l,2,ll) = czero		  
		  end if
 19               continue
               end do           ! nn = max(abs(mm),1),nmx(mm)
	    end do              ! n  = max(abs(m),1),nmx(m)
         end do                 ! mm = -mmax,Mup
      end do                    ! m  = -mmax,mmax
c-----------------------------------------------------------------------
c      The following relations are a consequence of the dihedral
c      symmetry:  
c      Q_(n,-m,n',-m')^(k,k') = 
c              (-1)^(m+m'+k+k') * Q_(n,m,n',m')^(k,k')    
c      CAUTION: The dihedral symmetry relation in this form assumes that
c      the coordinate system is chosen such that one of the dihedral 
c      axis coincides with the x-axis. This is important when 
c      implementing new geometries with dihedral symmetries!
c-----------------------------------------------------------------------
      if (DIHEDRAL) then
         do m = -mmax,mmax
            do mm = 1,mmax
               do n=max(1,abs(m)),nmx(m)
                  do nn = max(1,mm),nmx(mm)
                     nindxsum=n+nn
                     l=lsi(m,n)
                     lmx = l+lmax
                     ln=lsi(-m,n)
                     lnx= ln+lmax    
                     mindxsum=m+mm
                     indxsum=nindxsum+mindxsum
                     if (mod(mindxsum,2).eq.0) then
                        fac1=1.d0
                     else
                        fac1=-1.d0
                     end if
                     ll = lsi(mm,nn)
                     llmx = ll+lmax
		     call sym(Nsym,SIGMAh,m,mm,indxsum,
     &                     Cn,Sh,TQsym(l,ll))
		     if (Cn) then
                        Q(l,ll)=         czero
                        Q(lmx,llmx)=     czero
                        RgQ(1,l,1,ll)=   czero
                        RgQ(2,l,2,ll)=   czero
                        Q(lmx,ll)=       czero
                        Q(l,llmx)=       czero
                        RgQ(2,l,1,ll)=   czero
                        RgQ(1,l,2,ll)=   czero		     
		     else
                        lln = lsi(-mm,nn)
			llnx= lln+lmax
                        Q(l,ll)=         fac1*Q(ln,lln)
                        Q(lmx,llmx)=     fac1*Q(lnx,llnx)
                        RgQ(1,l,1,ll)=   fac1*RgQ(1,ln,1,lln)
                        RgQ(2,l,2,ll)=   fac1*RgQ(2,ln,2,lln)
                        Q(lmx,ll)=      -fac1*Q(lnx,lln)
                        Q(l,llmx)=      -fac1*Q(ln,llnx)
                        RgQ(2,l,1,ll)=  -fac1*RgQ(2,ln,1,lln)
                        RgQ(1,l,2,ll)=  -fac1*RgQ(1,ln,2,lln)
		     end if
	          end do        ! nn = abs(mm), nmx(mm)
	       end do           ! n=abs(m),nmx(m)
            end do              ! mm = 1,mmax
         end do                 ! m = -mmax,mmax
      end if

      return
      end
      
C***********************************************************************
c     calculation of the matrices Q, Rg Q for axisymmetric particles
C***********************************************************************

      subroutine Qmat_axisym(mmax,nmmax,nmx,lmax,lpmax,dim,thpts,thptsx,
     & phipts,phiptsx,Nsym,DIHEDRAL,SIGMAh,ng1,ng2,ng1dim,ng2dim,
     & xth1,xth2,xphi,r1,rth1,rphi1,r2,rth2,rphi2,
     & wth1,wth2,wphi,k,ks,mr,pi,nsi,lsi,lpsi,zero,cone,czero,Im,dtime,
     & Q,RgQ,TQsym)
      
      implicit none
      
c--------------------------------------------------------------------
c     Subroutine for calculating the Q- and RgQ-matrix.
c     C A U T I O N :
c     This routine has been written for calculating the T-matrix
c     from the Q- and RgQ-matrix. Thus any constant factors appearing 
c     in front of both the Q- and the RgQ-matrix have been omitted,
c     since they cancel when forming T = -RgQ * Q^{-1}. The output
c     of this routine is therefore (depending on the particle's 
c     symmetries) only identical with the Q- and RgQ-matrices up
c     to a constant factor. This is important when comparing Q- and
c     RgQ-matrices computed with different methods, or when one
c     uses this code for calculating the internal field inside the
c     particle.
c--------------------------------------------------------------------
      
c--------------------------------------------------------------------
c     input variables:
c--------------------------------------------------------------------
      integer mmax,nmmax,lmax,lpmax,dim,thpts,thptsx,phipts,
     &        phiptsx,Nsym,ng1dim,ng2dim,ng1,ng2,
     & nsi(nmmax),lsi(-mmax:mmax,nmmax),lpsi(0:mmax,nmmax),
     & nmx(-mmax:mmax)
      real*8 k,pi,zero,     
     & r1(ng1dim),rth1(ng1dim),rphi1(ng1dim),
     & r2(ng2dim,phipts),rth2(ng2dim,phipts),rphi2(ng2dim,phipts),
     & xth1(ng1dim),xth2(ng2dim,phipts),
     & wth1(ng1dim),wth2(ng2dim,phipts),
     & wphi(phipts),xphi(phipts),dtime
      complex*16 cone,czero,Im,ks,mr
      logical DIHEDRAL,SIGMAh

c--------------------------------------------------------------------
c     output variables:
c--------------------------------------------------------------------
      complex*16 Q(dim,dim),RgQ(2,lmax,2,lmax)
      integer TQsym(lmax,lmax)

c--------------------------------------------------------------------
c     locally defined variables:
c--------------------------------------------------------------------
      integer n,n1,nn,nn1,nindxsum,m,l,lp,lmx,mindxsum,ln,Mup,mm,
     &  indxsum,ll,llp,llmx,i,j,lnx,lln,llnx , mfac,mmfac,Sh
      real*8 sth1(ng1dim),sth2(thptsx,phiptsx),
     & a,fac,fac1,fac2,factmp,dm,dmm,mma,mmm,
     & tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,
     & P1(ng1dim,lpmax),dP1(ng1dim,lpmax),
     & P2(thptsx,phiptsx,lpmax),dP2(thptsx,phiptsx,lpmax),
     & dfac2(nmmax)  
      complex*16 b,cfac,mb,weight,
     & J11,J12,J21,J22,RgJ11,RgJ12,RgJ21,RgJ22,
     & Jt11,Jt12,Jt21,Jt22,RgJt11,RgJt12,RgJt21,RgJt22,
     & Jt12a,Jt21a,Jt22a,RgJt12a,RgJt21a,RgJt22a,
     & Btmp1,Btmp2,Btmp3,Btmp4,Btmp1a,Btmp2a,Btmp3a,Btmp4a,      
     & hn1(ng1dim,nmmax),wn1(ng1dim,nmmax),
     & jnn1(ng1dim,nmmax),
     & unn1(ng1dim,nmmax),
     & hn2(phiptsx,thptsx,nmmax),wn2(phiptsx,thptsx,nmmax),
     & jnn2(phiptsx,thptsx,nmmax),
     & unn2(phiptsx,thptsx,nmmax),eimphi(phiptsx,-2*mmax:2*mmax)
      logical Cn,lmsign
      
c-----------------------------------------------------------------------
c     preliminary settings
c-----------------------------------------------------------------------
      do n=1,nmmax
         dfac2(n) = dsqrt(dble(2*n+1)/dble(nsi(n)))
      end do
c-----------------------------------------------------------------------
c     calculation of Bessel-,Hankel-, and Wigner functions:
c-----------------------------------------------------------------------
      call logmsg('        Function calculations...    ')       
      call bessel(phipts,phiptsx,mmax,nmmax,
     &          ng1,ng2,ng1dim,ng2dim,thptsx,r1,r2,k,ks,
     &          xth1,xth2,xphi,cone,Im,hn1,wn1,jnn1,unn1,
     &          hn2,wn2,jnn2,unn2,eimphi,sth1,sth2)

      call wigner(phipts,phiptsx,nmmax,mmax,lpsi,lpmax,
     &          ng1,ng2,ng1dim,ng2dim,thptsx,xth1,xth2,P1,dP1,P2,dP2)

c-----------------------------------------------------------------------
c     calculation of Q and Rg Q:
c-----------------------------------------------------------------------
      
      lmsign=.false.
      if (DIHEDRAL) then
         Mup=0
      else
         Mup=mmax
      end if

      do m = -mmax,mmax
         dm = dble(m)
         do mm=-mmax,Mup
            dmm = dble(mm)
            mmm = dm*dmm
            do n=max(1,abs(m)),nmx(m)
               n1 = nsi(n)
               a = dble(n1)/k
               mma = dmm*a
               fac2 = dfac2(n) ! /(4.d0*pi)
               lp=lpsi(abs(m),n)
               do nn = max(1,abs(mm)),nmx(mm)
                  nn1=nsi(nn)
                  b = dcmplx(dble(nn1),zero)/ks
                  nindxsum=n+nn
                  factmp=fac2*dfac2(nn)
                  l=lsi(m,n)
                  lmx = l+lmax
                  mb = dcmplx(dm,zero)*b
                  mindxsum=nindxsum+m
                  llp=lpsi(abs(mm),nn)
	          indxsum=mindxsum+mm
	          ll=lsi(mm,nn)
		  llmx = ll+lmax
		  call sym(Nsym,SIGMAh,m,mm,indxsum,
     &                     Cn,Sh,TQsym(l,ll))
                  if(m.ne.mm) Cn=.true. ! symmetry property for axysymmetric particles
		  if (Cn) then
                     Q(l,ll)=czero
                     Q(lmx,ll)=czero
                     Q(l,llmx)=czero
                     Q(lmx,llmx)=czero
                     RgQ(1,l,1,ll)=czero
                     RgQ(2,l,1,ll)=czero
                     RgQ(1,l,2,ll)=czero
                     RgQ(2,l,2,ll)=czero
		     goto 19		     
		  end if
                  J11=czero
                  J12=czero
                  J21=czero
                  J22=czero
                  RgJ11=czero
                  RgJ12=czero
                  RgJ21=czero
                  RgJ22=czero
		  do j=1,phiptsx
                     Jt11=czero
                     Jt12=czero
                     Jt21=czero
                     Jt22=czero
                     RgJt11=czero
                     RgJt12=czero
                     RgJt21=czero
                     RgJt22=czero
                     if (DIHEDRAL) then
                        Jt12a=czero
                        Jt21a=czero
                        Jt22a=czero
                        RgJt12a=czero
                        RgJt21a=czero
                        RgJt22a=czero                     
                     end if
		     do i=1,thptsx
			if (lmsign) then
		           tmp1 = -P2(i,j,lp) * P2(i,j,llp)
			   tmp2 = -P2(i,j,lp) *dP2(i,j,llp)
			   tmp3 = -dP2(i,j,lp)* P2(i,j,llp)
			else
		           tmp1 =  P2(i,j,lp)* P2(i,j,llp)
			   tmp2 =  P2(i,j,lp)*dP2(i,j,llp)
			   tmp3 = dP2(i,j,lp)* P2(i,j,llp)
			end if
			if (Sh.ge.0) then
			   if (lmsign) then
			      tmp4 =-dP2(i,j,lp)*dP2(i,j,llp)
			   else
			      tmp4 = dP2(i,j,lp)*dP2(i,j,llp)
			   end if
			   tmp6 = r2(i,j)*(mmm*tmp1/sth2(i,j)+
     &                         tmp4*sth2(i,j))
			   tmp7 = rphi2(i,j)*tmp1/sth2(i,j)
			end if
			if (Sh.le.0) then			
			   tmp5 = r2(i,j)*(dm*tmp2+dmm*tmp3)
			   tmp8 = rth2(i,j)*tmp1
			   Btmp4  = wn2(j,i,n)*unn2(j,i,nn)
			   Btmp4a = dble(wn2(j,i,n))*unn2(j,i,nn)
			end if
			Btmp1  = hn2(j,i,n)*jnn2(j,i,nn)
			Btmp2  = hn2(j,i,n)*unn2(j,i,nn)
			Btmp3  = wn2(j,i,n)*jnn2(j,i,nn)
			Btmp1a = dble(hn2(j,i,n))*jnn2(j,i,nn)
			Btmp2a = dble(hn2(j,i,n))*unn2(j,i,nn)
			Btmp3a = dble(wn2(j,i,n))*jnn2(j,i,nn)
			if (.not.DIHEDRAL) then
                          if (Sh.ge.0) then
			    Jt12 = Jt12 + wth2(i,j)*(Btmp3*tmp6+a*Btmp1*
     &                           (rth2(i,j)*tmp2*sth2(i,j)+Im*dmm*tmp7))
                            Jt21 = Jt21 + wth2(i,j)*
     &                           (b*Btmp1*(rth2(i,j)*tmp3*sth2(i,j)
     &                           -Im*dm*tmp7)+Btmp2*tmp6)
			    RgJt12 = RgJt12 + wth2(i,j)*
     &                           (Btmp3a*tmp6+a*Btmp1a*
     &                           (rth2(i,j)*tmp2*sth2(i,j)+Im*dmm*tmp7))
                            RgJt21 = RgJt21 + wth2(i,j)*(b*Btmp1a*
     &                           (rth2(i,j)*tmp3*sth2(i,j)
     &                           -Im*dm*tmp7)+Btmp2a*tmp6)
                          end if
			  if (Sh.le.0) then
			    Jt11 = Jt11 + wth2(i,j)*Btmp1*tmp5
                            Jt22 = Jt22 + wth2(i,j)*(Im*Btmp4*tmp5+b*
     &                            Btmp3*(Im*dm*tmp8+rphi2(i,j)*tmp3) +
     &                            a*Btmp2*(Im*dmm*tmp8-rphi2(i,j)*tmp2))
			    RgJt11 = RgJt11 + wth2(i,j)*Btmp1a*tmp5
                            RgJt22 = RgJt22 + wth2(i,j)*
     &                            (Im*Btmp4a*tmp5+b*Btmp3a*
     &                            (Im*dm*tmp8+rphi2(i,j)*tmp3) +
     &                            a*Btmp2a*(Im*dmm*tmp8-
     &                            rphi2(i,j)*tmp2))
                          end if
                        else  ! if (DIHEDRAL)
                          if (Sh.ge.0) then
			    Jt12 = Jt12 + wth2(i,j)*(Btmp3*tmp6+a*Btmp1*
     &                            rth2(i,j)*tmp2*sth2(i,j))
			    Jt12a = Jt12a + wth2(i,j)*a*Btmp1*dmm*tmp7
                            Jt21 = Jt21 + wth2(i,j)*
     &                            (b*Btmp1*rth2(i,j)*tmp3*sth2(i,j)
     &                            +Btmp2*tmp6)
                            Jt21a = Jt21a - wth2(i,j)*b*Btmp1*dm*tmp7
			    RgJt12 = RgJt12 + wth2(i,j)*
     &                            (Btmp3a*tmp6+a*Btmp1a*rth2(i,j)*
     &                            tmp2*sth2(i,j))
			    RgJt12a = RgJt12a + wth2(i,j)*a*Btmp1a*
     &                            dmm*tmp7
                            RgJt21 = RgJt21 + wth2(i,j)*(b*Btmp1a*
     &                            rth2(i,j)*tmp3*sth2(i,j)
     &                            +Btmp2a*tmp6)
                            RgJt21a = RgJt21a - wth2(i,j)*b*Btmp1a*
     &                            dm*tmp7
                          end if
			  if (Sh.le.0) then
			    Jt11 = Jt11 + wth2(i,j)*Btmp1*tmp5
                            Jt22 = Jt22 + wth2(i,j)*(Btmp4*tmp5+
     &                            tmp8*(b*Btmp3*dm +
     &                            a*Btmp2*dmm))
                            Jt22a = Jt22a + wth2(i,j)*rphi2(i,j)*
     &                            (b*Btmp3*tmp3 - a*Btmp2*tmp2)
			    RgJt11 = RgJt11 + wth2(i,j)*Btmp1a*tmp5
                            RgJt22 = RgJt22 + wth2(i,j)*
     &                            (Btmp4a*tmp5 + tmp8*(b*Btmp3a*dm +
     &                            a*Btmp2a*dmm))
                            RgJt22a = RgJt22a + wth2(i,j)*rphi2(i,j)*
     &                            (b*Btmp3a*tmp3 - a*Btmp2a*tmp2)
                          end if                        
                        end if
		     end do ! i=1,thptsx
c-----------------------------------------------------------------------
		     weight = wphi(j)*eimphi(j,mm-m)
		     if (.not.DIHEDRAL) then
		        if (Sh.ge.0) then
		           J12 = J12 + weight*Jt12
		           J21 = J21 + weight*Jt21
		           RgJ12 = RgJ12 + weight*RgJt12
		           RgJ21 = RgJ21 + weight*RgJt21
		        end if
		        if (Sh.le.0) then	     
		           J11 = J11 + weight*Jt11
		           J22 = J22 + weight*Jt22
		           RgJ11 = RgJ11 + weight*RgJt11
		           RgJ22 = RgJ22 + weight*RgJt22
		        end if
		     else  ! if (DIHEDRAL)	     
		        if (Sh.ge.0) then
		           J12 = J12 + dreal(weight)*Jt12
     &                           - dimag(weight)*Jt12a
		           J21 = J21 + dreal(weight)*Jt21
     &                           - dimag(weight)*Jt21a
		           RgJ12 = RgJ12 + dreal(weight)*RgJt12
     &                           - dimag(weight)*RgJt12a
		           RgJ21 = RgJ21 + dreal(weight)*RgJt21
     &                           - dimag(weight)*RgJt21a
		        end if
		        if (Sh.le.0) then	     
		           J11 = J11 + dreal(weight)*Jt11
		           J22 = J22 + Im*(dreal(weight)*Jt22
     &                           +dimag(weight)*Jt22a)
		           RgJ11 = RgJ11 + dreal(weight)*RgJt11
		           RgJ22 = RgJ22 + Im*(dreal(weight)*RgJt22
     &                           +dimag(weight)*RgJt22a)
		        end if
		     end if
		  end do ! j=1,phiptsx
 10               continue
c-----------------------------------------------------------------------
c  factor (-1)^(m+mm) * sqrt((2n+1)(2n'+1)/(n(n+1)n'(n'+1)))/(4*pi)
c  in front of the integrals, where the constant 1/(4*pi) can be 
c  ommited, since Q^{-1} will get a factor 4*pi such that 
c  T = -Rg Q * Q^{-1} will be the same:
c-----------------------------------------------------------------------
                  fac=factmp
                  cfac=-Im*fac
c-----------------------------------------------------------------------

		  if (Sh.ge.0) then
                     J12= fac*J12
                     J21=-fac*J21
                     RgJ12= fac*RgJ12
                     RgJ21=-fac*RgJ21
                     Q(l,ll)       = J12+mr*J21
                     Q(lmx,llmx)   = J21+mr*J12
                     RgQ(1,l,1,ll) = RgJ12+mr*RgJ21
                     RgQ(2,l,2,ll) = RgJ21+mr*RgJ12
		  end if
                  if (Sh.le.0) then
                     J11=cfac*J11
                     J22=-fac*J22
                     RgJ11=cfac*RgJ11
                     RgJ22=-fac*RgJ22
                     Q(lmx,ll)     = J11+mr*J22
                     Q(l,llmx)     = J22+mr*J11
                     RgQ(2,l,1,ll) = RgJ11+mr*RgJ22
                     RgQ(1,l,2,ll) = RgJ22+mr*RgJ11
		  end if
		  if (Sh.eq.1) then
                     Q(lmx,ll)     = czero
                     Q(l,llmx)     = czero
                     RgQ(2,l,1,ll) = czero
                     RgQ(1,l,2,ll) = czero
		  else if (Sh.eq.-1) then
                     Q(l,ll)       = czero
                     Q(lmx,llmx)   = czero
                     RgQ(1,l,1,ll) = czero
                     RgQ(2,l,2,ll) = czero		  
		  end if
 19               continue
               end do           ! nn = max(abs(mm),1),nmx(mm)
	    end do              ! n  = max(abs(m),1),nmx(m)
         end do                 ! mm = -mmax,Mup
      end do                    ! m  = -mmax,mmax
c-----------------------------------------------------------------------
c      The following relations are a consequence of the dihedral
c      symmetry:  
c      Q_(n,-m,n',-m')^(k,k') = 
c              (-1)^(m+m'+k+k') * Q_(n,m,n',m')^(k,k')    
c      CAUTION: The dihedral symmetry relation in this form assumes that
c      the coordinate system is chosen such that one of the dihedral 
c      axis coincides with the x-axis. This is important when 
c      implementing new geometries with dihedral symmetries!
c-----------------------------------------------------------------------
      if (DIHEDRAL) then
         do m = -mmax,mmax
            do mm = 1,mmax
               do n=max(1,abs(m)),nmx(m)
                  do nn = max(1,mm),nmx(mm)
                     nindxsum=n+nn
                     l=lsi(m,n)
                     lmx = l+lmax
                     ln=lsi(-m,n)
                     lnx= ln+lmax    
                     mindxsum=m+mm
                     indxsum=nindxsum+mindxsum
                     if (mod(mindxsum,2).eq.0) then
                        fac1=1.d0
                     else
                        fac1=-1.d0
                     end if
                     ll = lsi(mm,nn)
                     llmx = ll+lmax
		     call sym(Nsym,SIGMAh,m,mm,indxsum,
     &                     Cn,Sh,TQsym(l,ll))
                  if(m.ne.mm) Cn=.true. ! symmetry property for axysymmetric particles  
		     if (Cn) then
                        Q(l,ll)=         czero
                        Q(lmx,llmx)=     czero
                        RgQ(1,l,1,ll)=   czero
                        RgQ(2,l,2,ll)=   czero
                        Q(lmx,ll)=       czero
                        Q(l,llmx)=       czero
                        RgQ(2,l,1,ll)=   czero
                        RgQ(1,l,2,ll)=   czero		     
		     else
                        lln = lsi(-mm,nn)
			llnx= lln+lmax
                        Q(l,ll)=         fac1*Q(ln,lln)
                        Q(lmx,llmx)=     fac1*Q(lnx,llnx)
                        RgQ(1,l,1,ll)=   fac1*RgQ(1,ln,1,lln)
                        RgQ(2,l,2,ll)=   fac1*RgQ(2,ln,2,lln)
                        Q(lmx,ll)=      -fac1*Q(lnx,lln)
                        Q(l,llmx)=      -fac1*Q(ln,llnx)
                        RgQ(2,l,1,ll)=  -fac1*RgQ(2,ln,1,lln)
                        RgQ(1,l,2,ll)=  -fac1*RgQ(1,ln,2,lln)
		     end if
	          end do        ! nn = abs(mm), nmx(mm)
	       end do           ! n=abs(m),nmx(m)
            end do              ! mm = 1,mmax
         end do                 ! m = -mmax,mmax
      end if

      return
      end
      
    
C***********************************************************************
c     calculation of the matrices Q, Rg Q for a sphere:
C***********************************************************************

      subroutine Qmat_sphere(mmax,nmmax,nmx,lmax,lpmax,dim,thpts,thptsx,
     & phipts,phiptsx,Nsym,DIHEDRAL,SIGMAh,ng1,ng2,ng1dim,ng2dim,
     & xth1,xth2,xphi,r1,rth1,rphi1,r2,rth2,rphi2,
     & wth1,wth2,wphi,k,ks,mr,pi,nsi,lsi,lpsi,zero,cone,czero,Im,dtime,
     & Q,RgQ,TQsym)
      
      implicit none
      
c--------------------------------------------------------------------
c     Subroutine for calculating the Q- and RgQ-matrix.
c     C A U T I O N :
c     This routine has been written for calculating the T-matrix
c     from the Q- and RgQ-matrix. Thus any constant factors appearing 
c     in front of both the Q- and the RgQ-matrix have been omitted,
c     since they cancel when forming T = -RgQ * Q^{-1}. The output
c     of this routine is therefore (depending on the particle's 
c     symmetries) only identical with the Q- and RgQ-matrices up
c     to a constant factor. This is important when comparing Q- and
c     RgQ-matrices computed with different methods, or when one
c     uses this code for calculating the internal field inside the
c     particle.
c--------------------------------------------------------------------
      
c--------------------------------------------------------------------
c     input variables:
c--------------------------------------------------------------------
      integer mmax,nmmax,lmax,lpmax,dim,thpts,thptsx,phipts,
     &        phiptsx,Nsym,ng1dim,ng2dim,ng1,ng2,
     & nsi(nmmax),lsi(-mmax:mmax,nmmax),lpsi(0:mmax,nmmax),
     & nmx(-mmax:mmax)
      real*8 k,pi,zero,     
     & r1(ng1dim),rth1(ng1dim),rphi1(ng1dim),
     & r2(ng2dim,phipts),rth2(ng2dim,phipts),rphi2(ng2dim,phipts),
     & xth1(ng1dim),xth2(ng2dim,phipts),
     & wth1(ng1dim),wth2(ng2dim,phipts),
     & wphi(phipts),xphi(phipts),dtime
      complex*16 cone,czero,Im,ks,mr
      logical DIHEDRAL,SIGMAh

c--------------------------------------------------------------------
c     output variables:
c--------------------------------------------------------------------
      complex*16 Q(dim,dim),RgQ(2,lmax,2,lmax)
      integer TQsym(lmax,lmax)

c--------------------------------------------------------------------
c     locally defined variables:
c--------------------------------------------------------------------
      integer n,n1,nn,nn1,nindxsum,m,l,lp,lmx,mindxsum,ln,Mup,mm,
     &  indxsum,ll,llp,llmx,i,j,lnx,lln,llnx , mfac,mmfac,Sh
      real*8 sth1(ng1dim),sth2(thptsx,phiptsx),
     & a,fac,fac1,fac2,factmp,dm,dmm,mma,mmm,
     & tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,
     & P1(ng1dim,lpmax),dP1(ng1dim,lpmax),
     & P2(thptsx,phiptsx,lpmax),dP2(thptsx,phiptsx,lpmax),
     & dfac2(nmmax)  
      complex*16 b,cfac,mb,weight,
     & J12,J21,RgJ12,RgJ21,Jt12,Jt21,RgJt12,RgJt21,
     & Jt12a,Jt21a,RgJt12a,RgJt21a,
     & Btmp1,Btmp2,Btmp3,Btmp4,Btmp1a,Btmp2a,Btmp3a,Btmp4a,      
     & hn1(ng1dim,nmmax),wn1(ng1dim,nmmax),
     & jnn1(ng1dim,nmmax),
     & unn1(ng1dim,nmmax),
     & hn2(phiptsx,thptsx,nmmax),wn2(phiptsx,thptsx,nmmax),
     & jnn2(phiptsx,thptsx,nmmax),
     & unn2(phiptsx,thptsx,nmmax),eimphi(phiptsx,-2*mmax:2*mmax)
      logical Cn,lmsign
      
c-----------------------------------------------------------------------
c     preliminary settings
c-----------------------------------------------------------------------
      do n=1,nmmax
         dfac2(n) = dsqrt(dble(2*n+1)/dble(nsi(n)))
      end do
c-----------------------------------------------------------------------
c     calculation of Bessel-,Hankel-, and Wigner functions:
c-----------------------------------------------------------------------
      call logmsg('        Function calculations...    ')       
      call bessel(phipts,phiptsx,mmax,nmmax,
     &          ng1,ng2,ng1dim,ng2dim,thptsx,r1,r2,k,ks,
     &          xth1,xth2,xphi,cone,Im,hn1,wn1,jnn1,unn1,
     &          hn2,wn2,jnn2,unn2,eimphi,sth1,sth2)

      call wigner(phipts,phiptsx,nmmax,mmax,lpsi,lpmax,
     &          ng1,ng2,ng1dim,ng2dim,thptsx,xth1,xth2,P1,dP1,P2,dP2)

c-----------------------------------------------------------------------
c     calculation of Q and Rg Q:
c-----------------------------------------------------------------------
      
      Mup=0
      lmsign=.false.

      do m = -mmax,mmax
         dm = dble(m)
         do mm=-mmax,Mup
            dmm = dble(mm)
            mmm = dm*dmm
            do n=max(1,abs(m)),nmx(m)
               n1 = nsi(n)
               a = dble(n1)/k
               mma = dmm*a
               fac2 = dfac2(n) ! /(4.d0*pi)
               lp=lpsi(abs(m),n)
               do nn = max(1,abs(mm)),nmx(mm)
                  nn1=nsi(nn)
                  b = dcmplx(dble(nn1),zero)/ks
                  nindxsum=n+nn
                  factmp=fac2*dfac2(nn)
                  l=lsi(m,n)
                  lmx = l+lmax
                  mb = dcmplx(dm,zero)*b
                  mindxsum=nindxsum+m
                  llp=lpsi(abs(mm),nn)
	          indxsum=mindxsum+mm
	          ll=lsi(mm,nn)
		  llmx = ll+lmax
                  Q(lmx,ll)=czero
                  Q(l,llmx)=czero
                  RgQ(2,l,1,ll)=czero
                  RgQ(1,l,2,ll)=czero
		  call sym(Nsym,SIGMAh,m,mm,indxsum,
     &                     Cn,Sh,TQsym(l,ll))
                  if(n.ne.nn.or.m.ne.mm) Cn=.true. ! symmetry for spheres
		  if (Cn) then
                     Q(l,ll)=czero
                     Q(lmx,llmx)=czero
                     RgQ(1,l,1,ll)=czero
                     RgQ(2,l,2,ll)=czero
		     goto 19		     
		  end if
                  J12=czero
                  J21=czero
                  RgJ12=czero
                  RgJ21=czero
		  do j=1,phiptsx
                     Jt12=czero
                     Jt21=czero
                     RgJt12=czero
                     RgJt21=czero
                     Jt12a=czero
                     Jt21a=czero
                     RgJt12a=czero
                     RgJt21a=czero
		     do i=1,thptsx
                        tmp1 =  P2(i,j,lp)* P2(i,j,llp)
                        tmp2 =  P2(i,j,lp)*dP2(i,j,llp)
                        tmp3 = dP2(i,j,lp)* P2(i,j,llp)
			if (Sh.ge.0) then
                           tmp4 = dP2(i,j,lp)*dP2(i,j,llp)
			   tmp6 = r2(i,j)*(mmm*tmp1/sth2(i,j)+
     &                         tmp4*sth2(i,j))
			   tmp7 = rphi2(i,j)*tmp1/sth2(i,j)
			end if
			Btmp1  = hn2(j,i,n)*jnn2(j,i,nn)
			Btmp2  = hn2(j,i,n)*unn2(j,i,nn)
			Btmp3  = wn2(j,i,n)*jnn2(j,i,nn)
			Btmp1a = dble(hn2(j,i,n))*jnn2(j,i,nn)
			Btmp2a = dble(hn2(j,i,n))*unn2(j,i,nn)
			Btmp3a = dble(wn2(j,i,n))*jnn2(j,i,nn)
                        if (Sh.ge.0) then
                           Jt12 = Jt12 + wth2(i,j)*(Btmp3*tmp6+a*Btmp1*
     &                          rth2(i,j)*tmp2*sth2(i,j))
                           Jt12a = Jt12a + wth2(i,j)*a*Btmp1*dmm*tmp7
                           Jt21 = Jt21 + wth2(i,j)*
     &                          (b*Btmp1*rth2(i,j)*tmp3*sth2(i,j)
     &                          +Btmp2*tmp6)
                           Jt21a = Jt21a - wth2(i,j)*b*Btmp1*dm*tmp7
                           RgJt12 = RgJt12 + wth2(i,j)*
     &                          (Btmp3a*tmp6+a*Btmp1a*rth2(i,j)*
     &                          tmp2*sth2(i,j))
                           RgJt12a = RgJt12a + wth2(i,j)*a*Btmp1a*
     &                          dmm*tmp7
                           RgJt21 = RgJt21 + wth2(i,j)*(b*Btmp1a*
     &                          rth2(i,j)*tmp3*sth2(i,j)
     &                          +Btmp2a*tmp6)
                           RgJt21a = RgJt21a - wth2(i,j)*b*Btmp1a*
     &                          dm*tmp7
                        end if
		     end do ! i=1,thptsx
c-----------------------------------------------------------------------
		     weight = wphi(j)*eimphi(j,mm-m)
                     if (Sh.ge.0) then
                        J12 = J12 + dreal(weight)*Jt12
     &                       - dimag(weight)*Jt12a
                        J21 = J21 + dreal(weight)*Jt21
     &                       - dimag(weight)*Jt21a
                        RgJ12 = RgJ12 + dreal(weight)*RgJt12
     &                       - dimag(weight)*RgJt12a
                        RgJ21 = RgJ21 + dreal(weight)*RgJt21
     &                       - dimag(weight)*RgJt21a
                     end if		    
		  end do ! j=1,phiptsx
 10               continue
c-----------------------------------------------------------------------
c  factor (-1)^(m+mm) * sqrt((2n+1)(2n'+1)/(n(n+1)n'(n'+1)))/(4*pi)
c  in front of the integrals, where the constant 1/(4*pi) can be 
c  ommited, since Q^{-1} will get a factor 4*pi such that 
c  T = -Rg Q * Q^{-1} will be the same:
c-----------------------------------------------------------------------
                  fac=factmp
                  if (mod(m+mm,2).ne.0) fac = -fac
                  cfac=-Im*fac
c-----------------------------------------------------------------------

		  if (Sh.ge.0) then
                     J12= fac*J12
                     J21=-fac*J21
                     RgJ12= fac*RgJ12
                     RgJ21=-fac*RgJ21
                     Q(l,ll)       = J12+mr*J21
                     Q(lmx,llmx)   = J21+mr*J12
                     RgQ(1,l,1,ll) = RgJ12+mr*RgJ21
                     RgQ(2,l,2,ll) = RgJ21+mr*RgJ12
		  end if
		  if (Sh.eq.-1) then
                     Q(l,ll)       = czero
                     Q(lmx,llmx)   = czero
                     RgQ(1,l,1,ll) = czero
                     RgQ(2,l,2,ll) = czero		  
		  end if
 19               continue
               end do           ! nn = max(abs(mm),1),nmx(mm)
	    end do              ! n  = max(abs(m),1),nmx(m)
         end do                 ! mm = -mmax,Mup
      end do                    ! m  = -mmax,mmax
c-----------------------------------------------------------------------
c      The following relations are a consequence of the dihedral
c      symmetry:  
c      Q_(n,-m,n',-m')^(k,k') = 
c              (-1)^(m+m'+k+k') * Q_(n,m,n',m')^(k,k')    
c      CAUTION: The dihedral symmetry relation in this form assumes that
c      the coordinate system is chosen such that one of the dihedral 
c      axis coincides with the x-axis. This is important when 
c      implementing new geometries with dihedral symmetries!
c-----------------------------------------------------------------------
      do m = -mmax,mmax
         do mm = 1,mmax
            do n=max(1,abs(m)),nmx(m)
               do nn = max(1,mm),nmx(mm)
                  nindxsum=n+nn
                  l=lsi(m,n)
                  lmx = l+lmax
                  ln=lsi(-m,n)
                  lnx= ln+lmax    
                  mindxsum=m+mm
                  indxsum=nindxsum+mindxsum
                  fac1=1.d0
                  ll = lsi(mm,nn)
                  llmx = ll+lmax
                  Q(lmx,ll)=       czero
                  Q(l,llmx)=       czero
                  RgQ(2,l,1,ll)=   czero
                  RgQ(1,l,2,ll)=   czero		     
                  call sym(Nsym,SIGMAh,m,mm,indxsum,
     &                 Cn,Sh,TQsym(l,ll))
                  if(m.ne.mm.or.n.ne.nn) Cn=.true.
                  if (Cn) then
                     Q(l,ll)=         czero
                     Q(lmx,llmx)=     czero
                     RgQ(1,l,1,ll)=   czero
                     RgQ(2,l,2,ll)=   czero
                  else
                     lln = lsi(-mm,nn)
                     llnx= lln+lmax
                     Q(l,ll)=         fac1*Q(ln,lln)
                     Q(lmx,llmx)=     fac1*Q(lnx,llnx)
                     RgQ(1,l,1,ll)=   fac1*RgQ(1,ln,1,lln)
                     RgQ(2,l,2,ll)=   fac1*RgQ(2,ln,2,lln)
                  end if
               end do           ! nn = abs(mm), nmx(mm)
            end do              ! n=abs(m),nmx(m)
         end do                 ! mm = 1,mmax
      end do                    ! m = -mmax,mmax

      return
      end
      
 
   
C***********************************************************************
c     preparing the Q-matrices for perturbation approach:
C***********************************************************************
      subroutine Qperturb(dim,Q,Q0)
      implicit none

      integer dim,i
      complex*16 Q(dim,dim),Q0(dim),cone

      cone=cmplx(1d0,0d0)

c     compute delta Q:
      do i=1,dim
         Q(i,i)=Q(i,i)-Q0(i)
      enddo
c     invert Q0 (trivial, since Q0 is a diagonal matrix):
      do i=1,dim
         if(zabs(Q0(i)).eq.0d0)then
            write(6,*)'FATAL ERROR: matrix Q0 singular'
            stop
         endif
         Q0(i)=cone/Q0(i)
      enddo

      end


C***********************************************************************
c     direct calculation of the T-matrix:
C***********************************************************************
      subroutine Tmatrix(mmax,nmmax,nmx,lmax,dim,lsi,DIHEDRAL,
     & SIGMAh,czero,Q,TQsym,Tmat)
      
      implicit none
      
c-----------------------------------------------------------------------
c     input variables (Tmat contains on input the matrix RgQ):
c-----------------------------------------------------------------------
      integer mmax,nmmax,lmax,dim,TQsym(lmax,lmax),
     &        lsi(-mmax:mmax,nmmax),nmx(-mmax:mmax)
      complex*16 czero,Q(dim,dim)
      logical DIHEDRAL,SIGMAh

c-----------------------------------------------------------------------
c     output variables:
c-----------------------------------------------------------------------
      complex*16 Tmat(2,lmax,2,lmax)

c-----------------------------------------------------------------------
c     locally defined variables:
c-----------------------------------------------------------------------
      integer n,n1,nn,nn1,m,l,lmx,ln,Mup,mm,
     &  ll,llmx,lln,lx,llx,k1,k2,mindxsum
      real*8 fac      
      complex*16 sum11,sum12,sum21,sum22,Ttmp(2,2,lmax)
      integer Sh,Shx
      
      if (DIHEDRAL) then
         Mup = 0
      else
         Mup = mmax
      end if                        

      do m=-mmax,mmax
         do n=max(1,abs(m)),nmx(m)
            l = lsi(m,n)
            ln= lsi(-m,n)
            lmx=l+lmax
            do mm=-mmax,Mup
               do nn=max(1,abs(mm)),nmx(mm)
                  ll = lsi(mm,nn)
                  llmx=ll+lmax
		  if (TQsym(l,ll).eq.3) then
		     do k1=1,2
		        do k2=1,2
                           Ttmp(k1,k2,ll)=czero
                        end do
                     end do
		     goto 29	
		  else
		     if (SIGMAh) then 
		        Sh=TQsym(l,ll)	     
                     else
		        Sh=0
	             end if
		  end if
                  sum11=czero
                  sum12=czero
                  sum21=czero
                  sum22=czero
		  if (Sh.ge.0) then
 		     do lx=1,lmax
		        llx=lx+lmax
		        if (TQsym(l,lx).ne.3) then
			   if (SIGMAh) then 
			      Shx=TQsym(l,lx)
			   else
			      Shx=0
			   end if
			   if (Shx.ge.0) then
     		              sum11=sum11-Tmat(1,l,1,lx)*q(lx,ll)
			      sum22=sum22-Tmat(2,l,2,lx)*q(llx,llmx)
			   end if
                           if (Shx.le.0) then
			      sum11=sum11-Tmat(1,l,2,lx)*q(llx,ll)
                              sum22=sum22-Tmat(2,l,1,lx)*q(lx,llmx)
			   end if
			end if
		     end do
	          end if
		  if (Sh.le.0) then
		     do lx=1,lmax
		        llx=lx+lmax
		        if (TQsym(l,lx).ne.3) then
			   if (SIGMAh) then 
			      Shx=TQsym(l,lx)
			   else
			      Shx=0
			   end if
			   if (Shx.ge.0) then
                              sum12=sum12-Tmat(1,l,1,lx)*q(lx,llmx)
                              sum21=sum21-Tmat(2,l,2,lx)*q(llx,ll)
			   end if
			   if (Shx.le.0) then
                              sum12=sum12-Tmat(1,l,2,lx)*q(llx,llmx)
                              sum21=sum21-Tmat(2,l,1,lx)*q(lx,ll)
			   end if
			end if
		     end do
		  end if
		  if (Sh.eq.1) then
		     sum12=czero
		     sum21=czero
		  else if (Sh.eq.-1) then
		     sum11=czero
		     sum22=czero
		  end if		  
                  Ttmp(1,1,ll)=sum11
                  Ttmp(1,2,ll)=sum12
                  Ttmp(2,1,ll)=sum21
                  Ttmp(2,2,ll)=sum22
 29               continue
               end do           ! nn=max(1,abs(mm)),nmx(mm)
	    end do              ! mm=-mmax,Mup
            do mm=-mmax,Mup
               do nn=max(1,abs(mm)),nmx(mm)
                  ll = lsi(mm,nn)
	         do k1=1,2
	            do k2=1,2
	               Tmat(k1,l,k2,ll) = Ttmp(k1,k2,ll)
	            end do
	         end do
	      end do            ! nn=max(1,abs(mm)),nmx(mm)
           end do               ! mm=-mmax,Mup
        end do                  ! n=max(1,abs(m)),nmx(m)
      end do                    ! m=-mmax,mmax
c-----------------------------------------------------------------------
c      The following relations are a consequence of the dihedral
c      symmetry:  
c      Q_(n,-m,n',-m')^(k,k') = 
c              (-1)^(m+m'+k+k') * Q_(n,m,n',m')^(k,k')    
c      CAUTION: The dihedral symmetry relation in this form assumes that
c      the coordinate system is chosen such that one of the dihedral 
c      axis coincides with the x-axis. This is important when 
c      implementing new geometries with dihedral symmetries!
c-----------------------------------------------------------------------
      if (DIHEDRAL) then
      do m=-mmax,mmax
         do n=max(1,abs(m)),nmx(m)
            l = lsi(m,n)
            ln= lsi(-m,n)
            do mm=1,mmax
               do nn=max(1,abs(mm)),nmx(mm)
                  ll = lsi(mm,nn)
		     if (TQsym(l,ll).eq.3) then
		        do k1=1,2
		           do k2=1,2
                              Tmat(k1,l,k2,ll)=czero
                           end do
                        end do
		     else
		        mindxsum=m+mm
	                if (mod(mindxsum,2).eq.0) then
	                   fac=1.d0
	                else
	                   fac=-1.d0
	                end if
		        lln=lsi(-mm,nn)
                        Tmat(1,l,1,ll)=  fac*Tmat(1,ln,1,lln)
                        Tmat(2,l,2,ll)=  fac*Tmat(2,ln,2,lln)
                        Tmat(2,l,1,ll)= -fac*Tmat(2,ln,1,lln)
                        Tmat(1,l,2,ll)= -fac*Tmat(1,ln,2,lln)
		     end if
	          end do        ! nn=max(1,abs(mm)),nmx(mm)
	       end do           ! mm=-mmax,Mup
	    end do              ! n=max(1,abs(m)),nmx(m)
         end do                 ! m=-mmax,mmax	  
      end if	  

      return
      end


C***********************************************************************
c     perturbative calculation of the T-matrix:
C***********************************************************************
      subroutine Tperturb(mmax,nmmax,nmx,lmax,dim,lsi,DIHEDRAL,
     & SIGMAh,czero,delQ,Q0inv,TQsym,RgQ,Tmat,prtb_order)
      
      implicit none
      
c-----------------------------------------------------------------------
c     input variables (Tmat contains on input the matrix RgQ):
c-----------------------------------------------------------------------
      integer mmax,nmmax,lmax,dim,prtb_order,TQsym(lmax,lmax),
     &        lsi(-mmax:mmax,nmmax),nmx(-mmax:mmax)
      complex*16 czero,delQ(dim,dim),RgQ(2,lmax,2,lmax),Q0inv(dim)
      logical DIHEDRAL,SIGMAh

c-----------------------------------------------------------------------
c     output variables:
c-----------------------------------------------------------------------
      complex*16 Tmat(2,lmax,2,lmax)

c-----------------------------------------------------------------------
c     locally defined variables:
c-----------------------------------------------------------------------
      integer n,n1,nn,nn1,m,l,lmx,ln,Mup,mm,
     &  ll,llmx,lln,lx,llx,k1,k2,mindxsum,porder
      real*8 fac      
      complex*16 sum11,sum12,sum21,sum22,Ttmp(2,2,lmax)
      integer Sh,Shx
      
c-----------------------------------------------------------------------
c     preliminary setups:
      if (DIHEDRAL) then
         Mup = 0
      else
         Mup = mmax
      end if
      do l=1,lmax
         do lx=1,lmax
            do k1=1,2
               do k2=1,2
                  Tmat(k1,l,k2,lx)=czero
               enddo
            enddo
         enddo
      enddo       
c-----------------------------------------------------------------------
      do porder=0,prtb_order    ! loop over perturbation order
         write(6,*)'        perturbation order ',porder
         if(porder.eq.0) goto 1100
c-----------------------------------------------------------------------
c     Compute T^{porder-1}*delQ
c-----------------------------------------------------------------------
         do m=-mmax,mmax
            do n=max(1,abs(m)),nmx(m)
               l = lsi(m,n)
               ln= lsi(-m,n)
               lmx=l+lmax
               do mm=-mmax,Mup
                  do nn=max(1,abs(mm)),nmx(mm)
                     ll = lsi(mm,nn)
                     llmx=ll+lmax
                     if (TQsym(l,ll).eq.3) then
                        do k1=1,2
                           do k2=1,2
                              Ttmp(k1,k2,ll)=czero
                           end do
                        end do
                        goto 29	
                     else
                        if (SIGMAh) then 
                           Sh=TQsym(l,ll)	     
                        else
                           Sh=0
                        end if
                     end if
                     sum11=czero
                     sum12=czero
                     sum21=czero
                     sum22=czero
                     if (Sh.ge.0) then
                        do lx=1,lmax
                           llx=lx+lmax
                           if (TQsym(l,lx).ne.3) then
                              if (SIGMAh) then 
                                 Shx=TQsym(l,lx)
                              else
                                 Shx=0
                              end if
                              if (Shx.ge.0) then
                                 sum11=sum11-Tmat(1,l,1,lx)*delQ(lx,ll)
                                 sum22=sum22-
     -                                Tmat(2,l,2,lx)*delQ(llx,llmx)
                              end if
                              if (Shx.le.0) then
                                 sum11=sum11-Tmat(1,l,2,lx)*delQ(llx,ll)
                                 sum22=sum22-
     -                                Tmat(2,l,1,lx)*delQ(lx,llmx)
                              end if
                           end if
                        end do
                     end if
                     if (Sh.le.0) then
                        do lx=1,lmax
                           llx=lx+lmax
                           if (TQsym(l,lx).ne.3) then
                              if (SIGMAh) then 
                                 Shx=TQsym(l,lx)
                              else
                                 Shx=0
                              end if
                              if (Shx.ge.0) then
                                 sum12=sum12-
     -                                Tmat(1,l,1,lx)*delQ(lx,llmx)
                                 sum21=sum21-Tmat(2,l,2,lx)*delQ(llx,ll)
                              end if
                              if (Shx.le.0) then
                                 sum12=sum12-
     -                                Tmat(1,l,2,lx)*delQ(llx,llmx)
                                 sum21=sum21-Tmat(2,l,1,lx)*delQ(lx,ll)
                              end if
                           end if
                        end do
                     end if
                     if (Sh.eq.1) then
                        sum12=czero
                        sum21=czero
                     else if (Sh.eq.-1) then
                        sum11=czero
                        sum22=czero
                     end if		  
                     Ttmp(1,1,ll)=sum11
                     Ttmp(1,2,ll)=sum12
                     Ttmp(2,1,ll)=sum21
                     Ttmp(2,2,ll)=sum22
 29                  continue
                  end do        ! nn=max(1,abs(mm)),nmx(mm)
               end do           ! mm=-mmax,Mup
               do mm=-mmax,Mup
                  do nn=max(1,abs(mm)),nmx(mm)
                     ll = lsi(mm,nn)
                     do k1=1,2
                        do k2=1,2
                           Tmat(k1,l,k2,ll) = Ttmp(k1,k2,ll)
                        end do
                     end do
                  end do        ! nn=max(1,abs(mm)),nmx(mm)
               end do           ! mm=-mmax,Mup
            end do              ! n=max(1,abs(m)),nmx(m)
         end do                 ! m=-mmax,mmax
c-----------------------------------------------------------------------
c     Compute -(RgQ+T^{porder-1} * delQ) * Q0^{-1}, 
c     where Q0^{-1} is a diagonal matrix
c-----------------------------------------------------------------------
 1100    do m=-mmax,mmax
            do n=max(1,abs(m)),nmx(m)
               l = lsi(m,n)
               do mm=-mmax,Mup
                  do nn=max(1,abs(mm)),nmx(mm)
                     ll = lsi(mm,nn)
                     do k1=1,2
                        do k2=1,2
                           llmx=ll+(k2-1)*lmax
                           Tmat(k1,l,k2,ll) = (Tmat(k1,l,k2,ll) -
     -                          RgQ(k1,l,k2,ll)) * Q0inv(llmx)
                        end do
                     end do
                  end do        
               end do           
            end do              
         end do                 
c-----------------------------------------------------------------------
c      The following relations are a consequence of the dihedral
c      symmetry:  
c      Q_(n,-m,n',-m')^(k,k') = 
c              (-1)^(m+m'+k+k') * Q_(n,m,n',m')^(k,k')    
c      CAUTION: The dihedral symmetry relation in this form assumes that
c      the coordinate system is chosen such that one of the dihedral 
c      axis coincides with the x-axis. This is important when 
c      implementing new geometries with dihedral symmetries!
c-----------------------------------------------------------------------
         if (DIHEDRAL) then
            do m=-mmax,mmax
               do n=max(1,abs(m)),nmx(m)
                  l = lsi(m,n)
                  ln= lsi(-m,n)
                  do mm=1,mmax
                     do nn=max(1,abs(mm)),nmx(mm)
                        ll = lsi(mm,nn)
                        if (TQsym(l,ll).eq.3) then
                           do k1=1,2
                              do k2=1,2
                                 Tmat(k1,l,k2,ll)=czero
                              end do
                           end do
                        else
                           mindxsum=m+mm
                           if (mod(mindxsum,2).eq.0) then
                              fac=1.d0
                           else
                              fac=-1.d0
                           end if
                           lln=lsi(-mm,nn)
                           Tmat(1,l,1,ll)=  fac*Tmat(1,ln,1,lln)
                           Tmat(2,l,2,ll)=  fac*Tmat(2,ln,2,lln)
                           Tmat(2,l,1,ll)= -fac*Tmat(2,ln,1,lln)
                           Tmat(1,l,2,ll)= -fac*Tmat(1,ln,2,lln)
                        end if
                     end do     ! nn=max(1,abs(mm)),nmx(mm)
                  end do        ! mm=-mmax,Mup
               end do           ! n=max(1,abs(m)),nmx(m)
            end do              ! m=-mmax,mmax	  
         end if	  
      enddo                     ! loop over perturbation order

      return
      end


C***********************************************************************
c     copy the l-th block of matrix Qirr into Qb
C***********************************************************************
      subroutine cpy_Q2Qb(Qirr,Qb,dim,lmax3,maxNc,ir_size,l)
      implicit none
      integer dim,lmax3,maxNc
      complex*16 Qirr(dim,dim),Qb(lmax3,lmax3)
      integer ir_size(maxNc),l,i,j,ii,jj,bstart

      bstart=1
      do i=1,l-1
         bstart=bstart+ir_size(i)
      end do

      do i=1,ir_size(l)
         ii=i+bstart-1
         do j=1,ir_size(l)
            jj=j+bstart-1
            Qb(i,j)=Qirr(ii,jj)
         end do
      end do

      return
      end


C***********************************************************************
c     copy the l-th block of matrix Qbinv into Qirrinv
C***********************************************************************
      subroutine cpy_Qb2Q(Qirrinv,Qbinv,dim,lmax3,maxNc,ir_size,l)
      implicit none
      integer dim,lmax3,maxNc
      complex*16 Qirrinv(dim,dim),Qbinv(lmax3,lmax3)
      integer ir_size(maxNc),l,i,j,ii,jj,bstart

      bstart=1
      do i=1,l-1
         bstart=bstart+ir_size(i)
      end do

      do i=1,ir_size(l)
         ii=i+bstart-1
         do j=1,ir_size(l)
            jj=j+bstart-1
            Qirrinv(ii,jj)=Qbinv(i,j)
         end do
      end do

      return
      end


C***********************************************************************
c     copy the T-matrix into another array with different indexing
c     convention
C***********************************************************************
      subroutine cpy_Tmat(Tmat,Tmat2,lmax,lmax2,lsi,nmx,nmax,mmax,nmmax)
      implicit none

c     input:
      integer lmax,lmax2,nmax,mmax,nmmax,
     &        lsi(-mmax:mmax,nmmax),nmx(-mmax:mmax)
      complex*16 Tmat(2,lmax,2,lmax)
c     output:
      complex*16 Tmat2(2,lmax2,2,lmax2)
c     local variables
      integer n,nn,m,mm,l,ll,l2,ll2,k,kk
      complex*16 czero
      parameter(czero=(0d0,0d0))


      do n=1,nmmax
         do m=-n,n
            l2=n*(n+1)+m
            do nn=1,nmmax
               do mm=-nn,nn
                  ll2=nn*(nn+1)+mm
                  do k=1,2
                     do kk=1,2
                        Tmat2(k,l2,kk,ll2)=czero
                     enddo
                  enddo
                  if(abs(m).le.mmax.and.abs(mm).le.mmax.and.
     &              n.le.nmax+abs(m).and.nn.le.nmax+abs(mm))then
                     l=lsi(m,n)
                     ll=lsi(mm,nn)
                     do k=1,2
                        do kk=1,2
                           Tmat2(k,l2,kk,ll2)=Tmat(k,l,kk,ll)
                        enddo
                     enddo                     
                  endif
               enddo
            enddo
         enddo
      enddo

      end


c***********************************************************************
C  
C    Compute total scattering and extinction cross sections for particle
C    in fixed orientation 
C    (hardcoded for thi=phi=0).
c
c     Author: Michael Kahnert
C  
c***********************************************************************
      subroutine crosssection(Tmat,nmax,lmax,thi,phi,ths,phs,k,
     &           Cext,Csca)
      implicit none

      integer n,m,nn,mm,l,ll,nmin,nnmin,nmax,lmax,pidx,ppidx
      integer nnn,lll,pppidx
      real*8 thi,phi,ths,phs,k,fr,Cext,Csca,pi,tcsca
      complex*16 Tmat(2,lmax,2,lmax),tmp1,Im,czero,cone,ccext
      parameter (czero=(0d0,0d0),cone=(1d0,0d0),Im=(0d0,1d0))

      parameter ( pi=3.1415926535897932385d0)
c-----------------------------------------------------------------------
c     Cext:
c-----------------------------------------------------------------------
      if(thi.ne.0.0.or.phi.ne.0.0)then
         write(6,*)
     &        'WARNING: Cext, Csca in fixed orientation only implemented
     & for theta_i=phi_i=0'
      endif
      cCext=czero
      do n=1,nmax
         l=n*(n+1)+1
         do nn=1,nmax
            ll=nn*(nn+1)+1
            tmp1=czero
            do pidx=1,2
               do ppidx=1,2
                  tmp1=tmp1 + conjg(Tmat(pidx,l,ppidx,ll))
               enddo
            enddo
            cCext=cCext-Im**(n-nn)*sqrt(dble((2*n+1)*(2*nn+1)))*tmp1
         enddo
      enddo
      Cext=2d0*pi/(k**2)*real(cCext)
      

c-----------------------------------------------------------------------
c     Csca:
c-----------------------------------------------------------------------
      Csca=0d0
      do m=-nmax,nmax
         nmin=max(abs(m),1)
         do n=nmin,nmax
            l=n*(n+1)+m
            do pidx=1,2

               tmp1=czero               
               do nn=1,nmax
                  ll=nn*(nn+1)+1
                  do ppidx=1,2
                     tmp1=tmp1+Tmat(pidx,l,ppidx,ll)*Im**(nn)*
     *                    sqrt(dble((2*nn+1)))
                  enddo
               enddo
               Csca=Csca+zabs(tmp1)**2
            enddo
         enddo
      enddo
      Csca=2d0*pi/(k**2) * Csca

      end
    
c***********************************************************************
c   The following is a rotation procedure for the T-matrix. 
c
c     Author: Michael Kahnert
c***********************************************************************
      subroutine rot_T(Tmat,Tmat3,alpha,beta,gamma,lmax,nmax,nbc,nfd,
     &                 nmax2,Nsym,DIHEDRAL,SIGMAh)
       
      implicit none
      
      integer i,j,n,m,nn,mm,m1,m2,nx,nnx,nxm,nnxmm,nxm1,nnxm2,
     &        l,nsum,m1sum,idxsum
      integer lmax,nmax,nbc,nmax2,Nsym,nfd
      real*8 alpha,beta,gamma,cbeta
      complex*16 Tmat(2,lmax,2,lmax),Ttmp(2,-nmax:nmax,2,-nmax:nmax),
     &           expalph(-nmax2:nmax2),expgamm(-nmax2:nmax2),
     &           Im,czero,cone,cfac0,cfac1
      complex*16 Tmat3(2,lmax,2,lmax)
      real*8 dc(-nmax:nmax,0:nfd)
      parameter (czero=(0d0,0d0),Im=(0d0,1d0),cone=(1d0,0d0))
      real*8 bcof(0:nbc,0:nbc),fnr(0:2*nbc)
      logical DIHEDRAL,SIGMAh,seven
      
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
c------------------------------------------------------------------------
c  calculation of Wigner d-functions
c------------------------------------------------------------------------
      cbeta = dcos(beta) 
      if (beta.ne.0.0) then     
         call rotcoef(cbeta,nmax,nmax2,dc,nmax,nmax2,nbc,nfd,bcof,fnr)
      end if

      expalph(0) = cone
      expgamm(0) = cone
      do m = 1, nmax2
         if (alpha.eq.0.0) then
            expalph(m)  = cone
            expalph(-m) = cone
         else
            expalph(m)  = exp( Im*m*alpha )
            expalph(-m) = conjg( expalph(m) )
         end if
         if (gamma.eq.0.0) then
            expgamm(m)  = cone
            expgamm(-m) = cone
         else
            expgamm(m)  = exp( Im*m*gamma )
            expgamm(-m) = conjg( expgamm(m) )
         end if
      end do
      
c------------------------------------------------------------------------
c  rotation of the T-matrix
c------------------------------------------------------------------------
      if (beta.ne.0.0) then 
       do n=1,nmax
         nx = n*(n+1)
         do nn=1,nmax
            nnx = nn*(nn+1)
            nsum=n+nn
            do m=-n,n
               nxm = nx+m
               do mm=-nn,nn
                  nnxmm = nnx+mm
                  do i=1,2
                     do j=1,2
                        Ttmp(i,m,j,mm) = czero
                     end do
                  end do
                  cfac0 = expalph(mm-m) 
                  do m1=-n,n
                     nxm1 = nx+m1
                     m1sum=nsum+m1
                     do m2=-nn,nn
                        if (mod(m1-m2,Nsym).ne.0) goto 1001
                        idxsum=m1sum+m2
                        if (mod(idxsum,2).eq.0) then
                           seven = .true.
                        else
                           seven = .false.
                        end if
                        nnxm2 = nnx+m2
                        cfac1 = cfac0 * expgamm(m2-m1) * dc(m,nxm1)
     &                          *dc(mm,nnxm2) 
                        if (seven.or. .not.SIGMAh) then
                           Ttmp(1,m,1,mm) = Ttmp(1,m,1,mm) +
     &                            Tmat(1,nxm1,1,nnxm2) * cfac1
                           Ttmp(2,m,2,mm) = Ttmp(2,m,2,mm) +
     &                            Tmat(2,nxm1,2,nnxm2) * cfac1
                        end if
                        if (.not.seven .or. .not.SIGMAh) then
                           Ttmp(1,m,2,mm) = Ttmp(1,m,2,mm) +
     &                            Tmat(1,nxm1,2,nnxm2) * cfac1
                           Ttmp(2,m,1,mm) = Ttmp(2,m,1,mm) +
     &                            Tmat(2,nxm1,1,nnxm2) * cfac1
                        end if
 1001                   continue
                     end do !loop over m2
                  end do !loop over m1
               end do !loop over mm
            end do !loop over m
            do m=-n,n
               nxm = nx+m
               do mm=-nn,nn
                  nnxmm = nnx+mm
                  do i=1,2
                     do j=1,2
                        Tmat3(i,nxm,j,nnxmm) = Ttmp(i,m,j,mm)
                     end do
                  end do
               end do !loop over mm
            end do !loop over m
         end do !loop over nn
       end do !loop over n
      else
       do n=1,nmax
         nx = n*(n+1)
         do nn=1,nmax
            nnx = nn*(nn+1)
            do m=-n,n
               nxm = nx+m
               do mm=-nn,nn
                  nnxmm = nnx+mm
                  if (mod(mm-m,Nsym).eq.0) then
                    cfac0 = expalph(mm-m)*expgamm(mm-m)  
                    do i=1,2
                       do j=1,2
                         Tmat3(i,nxm,j,nnxmm) = 
     &                         Tmat(i,nxm,j,nnxmm)*cfac0
                         if(i.eq.1.and.nxm.eq.1.and.j.eq.1.and.
     &                      nnxmm.eq.3) print*,'t:',Tmat3(i,nxm,j,nnxmm)
                        end do
                    end do
                  else
                    do i=1,2
                       do j=1,2
                          Tmat3(i,nxm,j,nnxmm) =czero
                       enddo
                    enddo
                  end if
               end do
            end do
         end do
       end do
      end if     
      
 1002 continue
      return
      end

      
      
c=======================================================================
c     PART III :  G R O U P   T H E O R Y 
c
c     IIIa.  Reducible representation:
c=======================================================================

c*********************************************************************** 
c     Tsym rogramme package for computing the reducible representations
c     of discrete symmetry operations in vector spherical functions
c
c     Author: Michael Kahnert
c     Date:   29.02.2004
c
c     Key paper:
c     
c        Michael Kahnert. Irreducible representations of finite groups
c        in the T-matrix formulation of the electromagnetic scattering
c        problem. J. Opt. Soc. Am. A 22, 1187-1199, 2005.
c
c     TODO: The special symmetry operations for the I and Ih
c           groups are not yet implemented
c***********************************************************************      

c=======================================================================
c     A. Reducible representations for non-cubic finite point groups
c=======================================================================

      subroutine sigma_v(dj,dNsym,pi,m,n,p,mt,fac)
c-----------------------------------------------------------------------
c     transformation by the j-th sigma_v reflection operation: 
c     V_{m,n,p} -> V_{mt,n,pt}=(-1)**(m+p)*exp(-i*4*pi*j*m/Nsym) * V_{-m,n,p}
c-----------------------------------------------------------------------
      implicit none
      integer m,n,p,mt
      complex*16 fac,arg
      real*8 dj,dNsym,pi,rfac,ifac,prec

      prec=1d-10

      mt=-m
      arg=4.d0*pi*dj*dble(m)/dNsym
      rfac=cos(arg)
      ifac=-sin(arg)
      if(dabs(rfac).lt.prec) rfac=0.d0
      if(dabs(ifac).lt.prec) ifac=0.d0
      fac=(-1.0)**(m+p)*cmplx(rfac,ifac)

      end


      subroutine sigma_d(dj,dNsym,pi,m,n,p,mt,fac)
c-----------------------------------------------------------------------
c     transformation by the j-th sigma_d reflection operation: 
c     V_{m,n,p} -> V_{mt,n,pt}=(-1)**(m+p)*exp(-i*4*pi*(j+1/2)*m/Nsym) * V_{-m,n,p}
c-----------------------------------------------------------------------
      implicit none
      integer m,n,p,mt
      complex*16 fac,arg
      real*8 dj,dNsym,pi,rfac,ifac,prec

      prec=1d-10

      mt=-m
      arg=4.d0*pi*(dj+0.5d0)*dble(m)/dNsym
      rfac=cos(arg)
      ifac=-sin(arg)
      if(dabs(rfac).lt.prec) rfac=0.d0
      if(dabs(ifac).lt.prec) ifac=0.d0
      fac=(-1.0)**(m+p)*cmplx(rfac,ifac)
      end


      subroutine Cn(dj,dNsym,pi,m,n,p,mt,fac)
c-----------------------------------------------------------------------
c     transformation by the rotation operation Cn^j: 
c     V_{m,n,p} -> V_{mt,n,pt}=exp(-i*2*pi*j*m/Nsym) * V_{m,n,p}
c-----------------------------------------------------------------------
      implicit none
      integer m,n,p,mt
      complex*16 fac,arg
      real*8 dj,dNsym,pi,rfac,ifac,prec

      prec=1d-10

      mt=m
      arg=2.d0*pi*dj*dble(m)/dNsym
      rfac=cos(arg)
      ifac=-sin(arg)
      if(dabs(rfac).lt.prec) rfac=0.d0
      if(dabs(ifac).lt.prec) ifac=0.d0
      fac=cmplx(rfac,ifac)
      end

      subroutine C2p(dj,dNsym,pi,m,n,p,mt,fac)
c-----------------------------------------------------------------------
c     transformation by the j-th dihedral rotation operation C'2: 
c     V_{m,n,p} -> V_{mt,n,pt}=(-1)**n * exp(-i*4*pi*j*m/Nsym) * V_{-m,n,p}
c-----------------------------------------------------------------------
      implicit none
      integer m,n,p,mt
      complex*16 fac,arg
      real*8 dj,dNsym,pi,rfac,ifac,prec

      prec=1d-10

      mt=-m
      arg=4.d0*pi*dj*dble(m)/dNsym
      rfac=cos(arg)
      ifac=-sin(arg)
      if(dabs(rfac).lt.prec) rfac=0.d0
      if(dabs(ifac).lt.prec) ifac=0.d0
      fac=(-1.0)**n * cmplx(rfac,ifac)
      end

      subroutine C2pp(dj,dNsym,pi,m,n,p,mt,fac)
c-----------------------------------------------------------------------
c     transformation by the j-th dihedral rotation operation C''2: 
c     V_{m,n,p} -> V_{mt,n,pt}=(-1)**n * exp(-i*4*pi*(j+1/2)*m/Nsym) * V_{-m,n,p}
c-----------------------------------------------------------------------
      implicit none
      integer m,n,p,mt
      complex*16 fac,arg
      real*8 dj,dNsym,pi,rfac,ifac,prec

      prec=1d-10

      mt=-m
      arg=4.d0*pi*(dj+0.5d0)*dble(m)/dNsym
      rfac=cos(arg)
      ifac=-sin(arg)
      if(dabs(rfac).lt.prec) rfac=0.d0
      if(dabs(ifac).lt.prec) ifac=0.d0
      fac=(-1.0)**n * cmplx(rfac,ifac)
      end


      subroutine inv(dj,dNsym,pi,m,n,p,mt,fac)
c-----------------------------------------------------------------------
c     transformation by the inversion operation i:
c     V_{m,n,p} -> V_{mt,n,pt}=(-1)**(n+p) * V_{m,n,p}
c-----------------------------------------------------------------------
      implicit none
      integer m,n,p,mt
      complex*16 fac
      real*8 dj,dNsym,pi

      mt=m
      fac=(-1.0)**(n+p)
      end

      subroutine sigma_h(dj,dNsym,pi,m,n,p,mt,fac)
c-----------------------------------------------------------------------
c     transformation by the reflection operation sigma_h:
c     V_{m,n,p} -> V_{mt,n,pt}=(-1)**(n+m+p) * V_{m,n,p}
c-----------------------------------------------------------------------
      implicit none
      integer m,n,p,mt
      complex*16 fac
      real*8 dj,dNsym,pi

      mt=m
      fac=(-1.0)**(n+m+p)
      end

      subroutine Sn(dj,dNsym,pi,m,n,p,mt,fac)
c-----------------------------------------------------------------------
c     transformation by the rotation-reflection operation Sn^j:
c     V_{m,n,p} -> V_{mt,n,pt}=(-1)**(n+m+p) * exp(-i*2*pi*j*m/Nsym) * V_{m,n,p}
c-----------------------------------------------------------------------
      implicit none
      integer m,n,p,mt
      complex*16 fac,arg
      real*8 dj,dNsym,pi,rfac,ifac,prec

      mt=m
      arg=2.d0*pi*dj*dble(m)/dNsym
      rfac=cos(arg)
      ifac=-sin(arg)
      fac=(-1.0)**(n+m+p)*cmplx(rfac,ifac)
      end

c=======================================================================
c     Characters of the reducible representations
c=======================================================================
      subroutine chi_sigma(chi)
c-----------------------------------------------------------------------
c     Characters of the reducible representations of sigma-operations
c-----------------------------------------------------------------------
      implicit none
      complex*16 chi

      chi=cmplx(0d0,0d0)

      end
      
      subroutine chi_Sn(chi)
c-----------------------------------------------------------------------
c     Characters of the reducible representation of Sn^(j)
c-----------------------------------------------------------------------
      implicit none
      complex*16 chi

      chi=cmplx(0d0,0d0)
      end
      

      subroutine chi_inv(chi)
c-----------------------------------------------------------------------
c     Characters of the reducible representation of the inversion operation
c-----------------------------------------------------------------------
      implicit none
      complex*16 chi

      chi=cmplx(0d0,0d0)
      end

      
      subroutine chi_C2(chi,mmax,nmx)
c-----------------------------------------------------------------------
c     Characters of the reducible representation of C'2^(j)
c-----------------------------------------------------------------------
      implicit none
      complex*16 chi
      integer mmax,nmx(-mmax:mmax)

      if(mod(nmx(0),2).eq.0)then
         chi=cmplx(0d0,0d0)
      else
         chi=cmplx(-2d0,0d0)
      endif
      end

      
      subroutine chi_E(chi,dim)
c-----------------------------------------------------------------------
c     Characters of the reducible representation of E
c-----------------------------------------------------------------------
      implicit none
      complex*16 chi
      integer dim

      chi=cmplx(dim,0d0)
      end

      
      subroutine chi_Cn(chi,dNsym,pi,nmx,mmax,j)
c-----------------------------------------------------------------------
c     Characters of the reducible representation Sn^(j)
c-----------------------------------------------------------------------
      implicit none
      integer mmax, nmx(-mmax:mmax),j
      complex*16 fac,chi
      real*8 dj,dNsym,pi,rfac,prec,rchi
      integer tau,m,n,mt

      rchi=0d0
      do m=1,mmax
         rfac=dble(4*(nmx(m)+1-m))
         dj=dble(j)
         rchi=rchi+rfac*cos(2d0*pi*dj*dble(m)/dNsym)
      enddo
      rchi=rchi+2d0*nmx(0)
      chi=cmplx(rchi,0d0)

c$$$      chi=cmplx(0d0,0d0)
c$$$      do tau=1,2
c$$$         do m=-mmax,mmax
c$$$            do n=max(1,abs(m)),nmx(m)
c$$$               dj=dble(j)
c$$$               call Cn(dj,dNsym,pi,m,n,tau,mt,fac)
c$$$               if(m.eq.mt)then
c$$$                  chi=chi+fac
c$$$               endif
c$$$            enddo
c$$$         enddo
c$$$      enddo

      end
      

      


c=======================================================================
c     B. Reducible representations for cubic finite point groups
c
c     (This part is work in progress; it has not been tested yet, 
c     nor is it currently referenced in the upper level routines.)
c
c     (currently implemented groups: T, Td, Th, O, Oh)
c         
c     All reducible representations contain a factor 
c     delta_{n,n'}delta(p,p'). This factor is suppressed in all 
c     comments below. d_{m,m'}^{(n)}(pi/2) is abbreviated by dmm'
c=======================================================================
 
      subroutine cubic_C2p(j,q,d,m,n,p,mm,mt,fac)
c-----------------------------------------------------------------------
c     representation of the rotation operation C'2^{(q)}:
c     q=0: (-1)**n        * i**(-m)   * delta_m,-m'              
c     q=1: (-1)**n        * i**m      * delta_m,-m'
c     q=2: (-1)**m'                   * dmm'         
c     q=3:                * i**(m+m') * dmm'              
c     q=4: (-1)**m                    * dmm'           
c     q=5:                * i**(-m-m')* dmm'              
c-----------------------------------------------------------------------
      implicit none
      integer j,q,m,n,p,mm,mt
      real*8 d
      complex*16 fac
      complex*16 CI
      parameter (CI=(0d0,1d0))

      mt=m
      if(q.eq.0)then
         fac=(-1d0)**n *CI**(-m)
         mt=-m
      elseif(q.eq.1)then
         fac=(-1d0)**n*CI**m
         mt=-m
      elseif(q.eq.2)then
         fac=(-1d0)**mm*d
      elseif(q.eq.3)then
         fac=CI**(m+mm)*d
      elseif(q.eq.4)then
         fac=(-1d0)**m*d
      elseif(q.eq.5)then
         fac=CI**(-m-mm)*d
      endif
      end



      subroutine cubic_C3(j,q,d,m,n,p,mm,mt,fac)
c-----------------------------------------------------------------------
c     representation of the rotation operation C3^{(q)j}:
c     q=0, j=1: (-1)**m'      * i**m   * dmm'
c     q=0, j=2:                 i**m'  * dmm'
c     q=1, j=1:                 i**-m  * dmm'
c     q=1, j=2: (-1)**m       * i**-m' * dmm'
c     q=2, j=1:                 i**m   * dmm'
c     q=2, j=2: (-1)**m       * i**m'  * dmm'
c     q=3, j=1: (-1)**(m+m')  * i**m   * dmm'
c     q=3, j=2:                 i**-m' * dmm'
c-----------------------------------------------------------------------
      implicit none
      integer j,q,m,n,p,mm,mt
      real*8 d
      complex*16 fac
      complex*16 CI
      parameter (CI=(0d0,1d0))

      mt=m
      if(q.eq.0)then
         if(j.eq.1)then
            fac=(-1d0)**mm*CI**m*d
         elseif(j.eq.2)then
            fac=CI**mm*d
         endif
      elseif(q.eq.1)then
         if(j.eq.1)then
            fac=CI**(-m)*d
         elseif(j.eq.2)then
            fac=(-1d0)**m*CI**(-mm)*d
         endif
      elseif(q.eq.2)then
         if(j.eq.1)then
            fac=CI**m*d
         elseif(j.eq.2)then
            fac=(-1d0)**m*CI**mm*d
         endif
      elseif(q.eq.3)then
         if(j.eq.1)then
            fac=(-1d0)**(m+mm)*CI**m*d
         elseif(j.eq.2)then
            fac=CI**(-mm)*d
         endif
      endif
      end



      subroutine cubic_C4(j,q,d,m,n,p,mm,mt,fac)
c-----------------------------------------------------------------------
c     representation of the rotation operation C4^{(q)j}:
c     q=0, j=1:               * i**-m      * delta_m,m'
c     q=0, j=2: (-1)**m                    * delta_m,m'
c     q=0, j=3:                 i**m       * delta_m,m'
c     q=1, j=1:                 i**(m-m')  * dmm'
c     q=1, j=2: (-1)**n                    * delta_m,-m'     
c     q=1, j=3:                 i**(m'-m)  * dmm'
c     q=2, j=1:                              dmm'
c     q=2, j=2: (-1)**(n+m)                  delta_m,-m'     
c     q=2, j=3: (-1)**(m+m')                 dmm'
c-----------------------------------------------------------------------
      implicit none
      integer j,q,m,n,p,mm,mt
      real*8 d
      complex*16 fac
      complex*16 CI
      parameter (CI=(0d0,1d0))

      mt=m
      if(q.eq.0)then
         if(j.eq.1)then
            fac=CI**(-m)
         elseif(j.eq.2)then
            fac=(-1d0)**m
         elseif(j.eq.3)then
            fac=CI**m
         endif
      elseif(q.eq.1)then
         if(j.eq.1)then
            fac=CI**(m-mm)*d
         elseif(j.eq.2)then
            fac=(-1d0)**n
            mt=-m
         elseif(j.eq.3)then
            fac=CI**(mm-m)*d
         endif
      elseif(q.eq.2)then
         if(j.eq.1)then
            fac=d
         elseif(j.eq.2)then
            fac=(-1d0)**(n+m)
            mt=-m
         elseif(j.eq.3)then
            fac=(-1d0)**(m+mm)*d
         endif
      endif
      end



      subroutine cubic_sigmad(j,q,d,m,n,p,mm,mt,fac)
c-----------------------------------------------------------------------
c     representation of the rotation operation sigma_d^{(q)}:
c     q=0: (-1)**(m+p)        * i**m       * delta_m,-m'              
c     q=1: (-1)**p            * i**m       * delta_m,-m'
c     q=2: (-1)**(n+m+p)                   * dmm'         
c     q=3: (-1)**(n+m'+p)                  * dmm'              
c     q=4: (-1)**(n+p)        * i**(-m-m') * dmm'           
c     q=5: (-1)**(n+p)        * i**(m+m')  * dmm'              
c-----------------------------------------------------------------------
      implicit none
      integer j,q,m,n,p,mm,mt
      real*8 d
      complex*16 fac
      complex*16 CI
      parameter (CI=(0d0,1d0))

      mt=m
      if(q.eq.0)then
         fac=(-1d0)**(m+p)*CI**m
         mt=-m
      elseif(q.eq.1)then
         fac=(-1d0)**p*CI**m
         mt=-m
      elseif(q.eq.2)then
         fac=(-1d0)**(n+m+p)*d
      elseif(q.eq.3)then
         fac=(-1d0)**(n+mm+p)*d
      elseif(q.eq.4)then
         fac=(-1d0)**(n+p)*CI**(-m-mm)*d
      elseif(q.eq.5)then
         fac=(-1d0)**(n+p)*CI**(m+mm)*d
      endif
      end



      subroutine cubic_sigmah(j,q,d,m,n,p,mm,mt,fac)
c-----------------------------------------------------------------------
c     representation of the rotation operation sigma_h^{(q)}:
c     q=0: (-1)**(n+m+p)      * i**m       * delta_m,m'              
c     q=1: (-1)**(m+p)                     * delta_m,-m'
c     q=2: (-1)**p                         * delta_m,-m'
c-----------------------------------------------------------------------
      implicit none
      integer j,q,m,n,p,mm,mt
      real*8 d
      complex*16 fac
      complex*16 CI
      parameter (CI=(0d0,1d0))

      mt=m
      if(q.eq.0)then
         fac=(-1d0)**(n+m+p)
      elseif(q.eq.1)then
         fac=(-1d0)**(m+p)
         mt=-m
      elseif(q.eq.2)then
         fac=(-1d0)**p
         mt=-m
      endif
      end



      subroutine cubic_S4(j,q,d,m,n,p,mm,mt,fac)
c-----------------------------------------------------------------------
c     representation of the rotation operation S4^{(q)j}:
c     q=0, j=1: (-1)**(n+p)     * i**m       * delta_m,m'
c     q=0, j=3: (-1)**(n+p)     * i**-m      * delta_m,m'
c     q=1, j=1: (-1)**(n+p)     * i**(-m+m') * dmm'
c     q=1, j=3: (-1)**(n+p)     * i**(m-m')  * dmm'
c     q=2, j=1: (-1)**(n+m+m'+p)             * dmm'
c     q=2, j=3: (-1)**(n+p)                  * dmm'
c-----------------------------------------------------------------------
      implicit none
      integer j,q,m,n,p,mm,mt
      real*8 d
      complex*16 fac
      complex*16 CI
      parameter (CI=(0d0,1d0))

      mt=m
      if(q.eq.0)then
         if(j.eq.1)then
            fac=(-1d0)**(n+p)*CI**m
         elseif(j.eq.3)then
            fac=(-1d0)**(n+p)*CI**(-m)
         endif
      elseif(q.eq.1)then
         if(j.eq.1)then
            fac=(-1d0)**(n+p)*CI**(-m+mm)*d
         elseif(j.eq.3)then
            fac=(-1d0)**(n+p)*CI**(m-mm)*d
         endif
      elseif(q.eq.2)then
         if(j.eq.1)then
            fac=(-1d0)**(n+m+mm+p)*d
         elseif(j.eq.3)then
            fac=(-1d0)**(n+p)*d
         endif
      endif
      end



      subroutine cubic_S6(j,q,d,m,n,p,mm,mt,fac)
c-----------------------------------------------------------------------
c     representation of the rotation operation S6^{(q)j}:
c     q=0, j=1: (-1)**(n+p)     * i**m'   * dmm'
c     q=0, j=5: (-1)**(n+mm+p)  * i**m    * dmm'  
c     q=1, j=1: (-1)**(n+p)     * i**m    * dmm'
c     q=1, j=5: (-1)**(n+m+p)   * i**m'   * dmm'
c     q=2, j=1: (-1)**(n+m+p)   * i**(-m')* dmm'
c     q=2, j=5: (-1)**(n+p)     * i**(-m) * dmm'
c     q=3, j=1: (-1)**(n+m'+p)  * i**(-m) * dmm'
c     q=3, j=5: (-1)**(n+p)     * i**(-m')* dmm'
c-----------------------------------------------------------------------
      implicit none
      integer j,q,m,n,p,mm,mt
      real*8 d
      complex*16 fac
      complex*16 CI
      parameter (CI=(0d0,1d0))

      mt=m
      if(q.eq.0)then
         if(j.eq.1)then
            fac=(-1d0)**(n+p)*CI**mm*d
         elseif(j.eq.5)then
            fac=(-1d0)**(n+mm+p)*CI**m*d
         endif
      elseif(q.eq.1)then
         if(j.eq.1)then
            fac=(-1d0)**(n+p)*CI**m*d
         elseif(j.eq.5)then
            fac=(-1d0)**(n+m+p)*CI**mm*d
         endif
      elseif(q.eq.2)then
         if(j.eq.1)then
            fac=(-1d0)**(n+m+p)*CI**(-mm)*d
         elseif(j.eq.5)then
            fac=(-1d0)**(n+p)*CI**(-m)*d
         endif
      elseif(q.eq.3)then
         if(j.eq.1)then
            fac=(-1d0)**(n+mm+p)*CI**(-m)*d
         elseif(j.eq.5)then
            fac=(-1d0)**(n+p)*CI**(-mm)*d
         endif
      endif
      end



      subroutine cubic_inv(j,q,d,m,n,p,mm,mt,fac)
c-----------------------------------------------------------------------
c     representation of the inversion operation I:
c     (-1)**(n+p)
c-----------------------------------------------------------------------
      implicit none
      integer j,q,m,n,p,mm,mt
      real*8 d
      complex*16 fac
      complex*16 CI
      parameter (CI=(0d0,1d0))

      mt=m
      fac=(-1.0)**(n+p)
      end


c=======================================================================
c     IIIb.  Irreducible representation:
c=======================================================================

c*********************************************************************** 
c     Tsym program package for irreducible representations of finite
c     groups to be used in T-matrix computations
c
c     Author: Michael Kahnert
c     Date:   29.02.2004
c     Update: 16.04.2012
c
c     Key papers:
c     
c        Michael Kahnert. Irreducible representations of finite groups
c        in the T-matrix formulation of the electromagnetic scattering
c        problem. J. Opt. Soc. Am. A 22, 1187-1199, 2005.
c
c        Michael Kahnert. T-matrix computations for particles with 
c        high-order finite symmetries. JQSRT 2012, 
c        http://dx.doi.org/10.1016/j.jqsrt.2012.08.004
c
c     TODO:    The cubic groups T, Th, Td, O, Oh, I, and Ih 
c              are not yet implemented
c***********************************************************************      
      subroutine IrredRep(Pgroup,Nsym,maxNc,dim,mmax,nmmax,lmax,nmx,
     &     lsi,Nc,No,Nclass,Norder,prec,ne,neinv,ir_size,chtabfile,
     &     Pmatmethod)
      use pmat
      implicit none
c-----------------------------------------------------------------------
c     INPUT:
c
c     Pgroup:   character string representing the kind of point group 
c     Nsym:     integer characterising the main symmetry axis (for Cn,
c               Cnv, Cnh, Dnh, Dnd, and Sn groups)
c     chtabfile:Name of file containing character table
c     Pmatmethod: If .eq.1, construct the transformation into the 
c               irreducible basis by finding the linearly independent
c               row vectors of the projection matrices; if .eq.2,
c               use LAPACK SVD; if .eq. 3, use PROPACK SVD.
c
c     OUTPUT:   (most of these variables are defined in module pmat)
c
c     P:        transformation matrix for transforming the Q matrix 
c               from the reducible into the irreducible basis:
c               Qirr=P*Q*P^{-1}. P is highly sparse and therefore
c               stored in compact format: P_{i,j}=P(k), where
c     ip:       contains indexing information of P: ip(k)=i, ip(k+ne)=j
c     ne:       number of non-zero elements of P
c     Pinv:     inverse transformation P^{-1}
c     ipinv:    indexing information for Pinv
c     neinv:    number of non-zero elements of Pinv
c     ir_size:  integer array, where ir_size(l)=dimension of
c               the l-th invariant subspace
c   
c     LOCAL VARIABLES:
c   
c     Nc:       number of conjugacy classes in the point group
c     No:       order of the point group (=number of elements)
c     Nclass:   integer array containing the number of different types 
c               of conjugacy classes (see comment in subroutine
c               get_class_order) 
c     Norder:   integer array containing the number of elements in each 
c               class of the group (see comment in 
c               subroutine get_class_order)
c-----------------------------------------------------------------------
      character*3 Pgroup
      character*(*) chtabfile
      integer Nsym, maxNc, dim,mmax,nmmax,lmax,
     &        nmx(-mmax:mmax),lsi(-mmax:mmax,nmmax), 
     &        ne, neinv,alpha(maxNc),ir_size(maxNc),Pmatmethod
      real*8 prec
c-----------------------
      integer i,j,k, sdim
      integer Nc, No, Nclass(maxNc), Norder(maxNc,maxNc)
      complex*16 chi(maxNc,maxNc) 

c-----------------------------------------------------------------------
c     compute some properties of the group
c-----------------------------------------------------------------------
      call get_class_order(Pgroup,Nsym,Nc,No,Nclass,Norder,maxNc)
c-----------------------------------------------------------------------
c     read in the character table
c-----------------------------------------------------------------------
      call get_characters(Pgroup,Nsym,maxNc,Nc,chi,chtabfile)
c-----------------------------------------------------------------------
c     Number of times an irreducible representation occurs in a 
c     reducible one, as well as dimensions of the invariant subpaces
c-----------------------------------------------------------------------
      call n_irr(Pgroup,dim,Nsym,No,Nc,Nclass,Norder,
     &     chi,maxNc,mmax,nmx,alpha,ir_size)
c-----------------------------------------------------------------------
c     compute the P matrix that transforms from the reducible into the
c     irreducible basis
c-----------------------------------------------------------------------
      sdim=dim*No
      call get_Pmat(Pgroup,dim,Nsym,maxNc,sdim,
     &              Nc,No,Nclass,Norder,mmax,nmmax,lmax,nmx,lsi,
     &              chi,prec,ne,ir_size,Pmatmethod)
c-----------------------------------------------------------------------
c     compute P^{-1} 
c-----------------------------------------------------------------------
      call invP(dim,ne,sdim,neinv)

      return
      end

      

c*********************************************************************** 
c     determine the number of conjugacy classes, the order, and the 
c     symmetry elements of the point group
c
c     Author: Michael Kahnert
c     Date:   17.07.2004
c***********************************************************************      
      subroutine get_class_order(Pgroup,Nsym,Nc,No,Nclass,Norder,
     &                           maxNc)
      implicit none
c-----------------------------------------------------------------------
c     INPUT:
c     Pgroup:   character string representing the kind of point group 
c     Nsym:     integer characterising the main symmetry axis (for Cn,
c               Cnv, Cnh, Dnh, Dnd, and Sn groups)
c   
c     OUTPUT:
c     Nc:       total number of classes in the point group
c     No:       order of the point group
c     Nclass:   integer array, where
c               Nclass(1)=number of E classes=1 (trivial, but convenient
c                         to include)
c               Nclass(2)=number of C_Nsym^j classes
c               Nclass(3)=number of C'_2^j classes
c               Nclass(4)=number of C''_2^j classes
c               Nclass(5)=number of sigma_v^j classes
c               Nclass(6)=number of sigma_d^j classes
c               Nclass(7)=1 if sigma_h symmetry exists, otherwise 0
c               Nclass(8)=1 if inversion centre exists, otherwise 0
c               Nclass(9)=number of S_Nsym^(j) classes, where
c                         S_Nsym^(j)=sigma_h C_Nsym^j
c     Norder:   integer array, where
c               Norder(j,1),...,Norder(j,Nclass(j)): number of elements
c                  in the classes of type j, where 
c                  j=1: E  
c                  j=2: Cn 
c                  j=3: C'2
c                  j=4: C''2
c                  j=5: sigma_v
c                  j=6: sigma_d
c                  j=7: sigma_h
c                  j=8: i
c                  j=9: Sn
c               Example: In D8h we have 4 C_8^j classes, namely
c               {C_8, C_8^7}, {C_4, C_4^3}, {C_8^3, C_8^5}, {C_2},
c               so Nclass(2)=4, and 
c               Norder(2,1)=2, Norder(2,2)=2, Norder(2,3)=2, Norder(2,4)=1
c
c               Larger indices of Nclass and Norder will be used in 
c               coming versions for special symmetries of the groups 
c               Td, O, Oh, and Ih.
c
c     NOTE: The class containing the identity operation, which is contained
c           in any point group, is not explicitly
c           stored in the arrays above, but it is, of course, accounted
c           for in Nc and No.
c     NOTE: if Nsym is even then S_Nsym^(j)=i (inversion element). 
c     This element is counted as i, not as S_Nsym^(j) in Nclass and Norder!
c     NOTE: for the group Dnd, the S-symmetry is S_{2*Nsym}^(j) where j is odd;
c           for the group Sn the C-symmetry is C_{Nsym/2}^j
c-----------------------------------------------------------------------
      character*3 Pgroup
      integer maxNc
      integer Nsym, No, Nc, Nclass(maxNc), Norder(maxNc,maxNc)
      integer j,k

c-----------------------------------------------------------------------
c     preliminary checks and settings:
c-----------------------------------------------------------------------
      if (Pgroup.eq.'Sn'.and.(mod(Nsym,2).ne.0)) then
         Pgroup="Cnh"
         write(6,*) 
     &   "NOTE: Sn group for odd n is isomorphic to Cnh group!"
      end if
      if (Pgroup.eq.'T') then
         write(6,*) 
     &   "WARNING: T group not yet implemented - using subgroup C3 inste
     &ad"
         Pgroup="Cn"
         Nsym=3
      end if
      if (Pgroup.eq.'Td') then
         write(6,*) 
     &   "WARNING: Td group not yet implemented - using subgroup C3 inst
     &ead"
         Pgroup="Cn"
         Nsym=3
      end if
      if (Pgroup.eq.'Th') then
         write(6,*) 
     &   "WARNING: Th group not yet implemented - using subgroup C3 inst
     &ead"
         Pgroup="Cn"
         Nsym=3
      end if
      if (Pgroup.eq.'O') then
         write(6,*)  
     &   "WARNING: O group not yet implemented - using subgroup C4 inste
     &ad"
         Pgroup="Cn"
         Nsym=4
      end if
      if (Pgroup.eq.'Oh') then
         write(6,*)  
     &   "WARNING: Oh group not yet implemented - using subgroup D4h ins
     &tead"
         Pgroup="Dnh"
         Nsym=4
      end if
      if (Pgroup.eq.'Ih') then
         write(6,*)  
     &   "ERROR: Icosahedral group Ih not yet implemented"
         stop
      end if

      if (Pgroup.ne."Cs".and.Pgroup.ne."Ci".and.Nsym.lt.2) then
         write(6,*) "ERROR: Nsym must be >= 2"
         stop
      end if

c-----------------------------------------------------------------------
c     number of classes and symmetry elements:
c-----------------------------------------------------------------------
      do j=1,maxNc
         Nclass(j)=0
         do k=1,maxNc
            Norder(j,k)=0
         end do
      end do
c-----------------------------------------------------------------------
c     class containing the unity element E (the same for all groups):
c-----------------------------------------------------------------------
      Nclass(1)=1
      Norder(1,1)=1
c-----------------------------------------------------------------------
      if (Pgroup.eq."Cs") then
         Nc=2
         No=2
         Nclass(7)=1 
         Norder(7,1)=1                  ! sigma_h
c-----------------------------------------------------------------------
      else if(Pgroup.eq."Ci") then
         Nc=2
         No=2
         Nclass(8)=1
         Norder(8,1)=1                  ! i
c-----------------------------------------------------------------------
      else if(Pgroup.eq."Cn") then
         Nc=Nsym
         No=Nsym     
         Nclass(2)=Nsym-1
         do j=1,Nclass(2)
            Norder(2,j)=1               ! Cn classes
         end do
c-----------------------------------------------------------------------
      else if(Pgroup.eq."Dn") then
         No=2*Nsym         
         if (mod(Nsym,2).ne.0) then
            Nc=(Nsym+3)/2
            Nclass(2)=(Nsym-1)/2
            Nclass(3)=1
            do j=1,Nclass(2)
               Norder(2,j)=2            ! Cn classes
            end do
            do j=1,Nclass(3)
               Norder(3,j)=Nsym         ! C'2 class
            end do
         else
            Nc=Nsym/2+3
            Nclass(2)=Nsym/2
            Nclass(3)=1
            Nclass(4)=1
            do j=1,Nclass(2)-1
               Norder(2,j)=2            ! Cn classes
            end do
            Norder(2,Nclass(2))=1
            do j=1,Nclass(3)
               Norder(3,j)=Nsym/2       ! C'2 class
            end do
            do j=1,Nclass(4)
               Norder(4,j)=Nsym/2       ! C''2 class
            end do
         end if
c-----------------------------------------------------------------------
      else if(Pgroup.eq."Cnv") then
         No=2*Nsym         
         if (mod(Nsym,2).ne.0) then
            Nc=(Nsym+3)/2
            Nclass(2)=(Nsym-1)/2
            Nclass(5)=1
            do j=1,Nclass(2)
               Norder(2,j)=2            ! Cn classes
            end do
            do j=1,Nclass(5)
               Norder(5,j)=Nsym         ! sigma_v class
            end do
         else
            Nc=Nsym/2+3
            Nclass(2)=Nsym/2
            Nclass(5)=1
            Nclass(6)=1
            do j=1,Nclass(2)-1
               Norder(2,j)=2            ! Cn classes
            end do
            Norder(2,Nclass(2))=1
            do j=1,Nclass(5)
               Norder(5,j)=Nsym/2       ! sigma_v class
            end do
            do j=1,Nclass(6)
               Norder(6,j)=Nsym/2       ! sigma_d class
            end do
         end if
c-----------------------------------------------------------------------
      else if(Pgroup.eq."Cnh") then
         Nc=2*Nsym
         No=2*Nsym    
         Nclass(2)=Nsym-1
         Nclass(7)=1
         do j=1,Nclass(2)  
            Norder(2,j)=1               ! Cn classes
         end do
         Norder(7,1)=1                  ! sigma_h
         if (mod(Nsym,2).ne.0) then
            Nclass(9)=Nsym-1
            do j=1,Nclass(9)
               Norder(9,j)=1            ! Sn classes
            end do
         else
            Nclass(8)=1
            Nclass(9)=Nsym-2
            Norder(8,1)=1               ! i
            do j=1,Nclass(9)
               Norder(9,j)=1            ! Sn classes
            end do
         end if
c-----------------------------------------------------------------------
      else if(Pgroup.eq."Dnh") then
         No=4*Nsym         
         if (mod(Nsym,2).ne.0) then
            Nc=Nsym+3
            Nclass(2)=(Nsym-1)/2
            Nclass(3)=1
            Nclass(5)=1
            Nclass(7)=1
            Nclass(9)=(Nsym-1)/2
            do j=1,Nclass(2)  
               Norder(2,j)=2            ! Cn classes
            end do
            do j=1,Nclass(3)
               Norder(3,j)=Nsym         ! C'2 classes
            end do
            do j=1,Nclass(5)
               Norder(5,j)=Nsym         ! sigma_v classes
            end do
            Norder(7,1)=1               ! sigma_h
            do j=1,Nclass(9)
               Norder(9,j)=2            ! Sn classes
            end do            
         else
            Nc=Nsym+6
            Nclass(2)=Nsym/2
            Nclass(3)=1
            Nclass(4)=1
            Nclass(5)=1
            Nclass(6)=1
            Nclass(7)=1
            Nclass(8)=1
            Nclass(9)=Nsym/2-1
            do j=1,Nclass(2)-1
               Norder(2,j)=2            ! Cn classes
            end do
            Norder(2,Nclass(2))=1
            do j=1,Nclass(3)
               Norder(3,j)=Nsym/2       ! C'2 classes
            end do
            do j=1,Nclass(4)
               Norder(4,j)=Nsym/2       ! C''2 classes
            end do
            do j=1,Nclass(5)
               Norder(5,j)=Nsym/2       ! sigma_v classes
            end do
            do j=1,Nclass(6)
               Norder(6,j)=Nsym/2       ! sigma_d classes
            end do
            Norder(7,1)=1               ! sigma_h
            Norder(8,1)=1               ! i
            do j=1,Nclass(9)
               Norder(9,j)=2              ! Sn classes
            end do
         end if
c-----------------------------------------------------------------------
      else if(Pgroup.eq."Dnd") then
         No=4*Nsym         
         Nc=Nsym+3
         if (mod(Nsym,2).ne.0) then
            Nclass(2)=(Nsym-1)/2
            Nclass(3)=1
            Nclass(6)=1
            Nclass(8)=1
            Nclass(9)=(Nsym-1)/2
            do j=1,Nclass(2)  
               Norder(2,j)=2            ! Cn classes
            end do
            do j=1,Nclass(3)
               Norder(3,j)=Nsym         ! C'2 classes
            end do
            do j=1,Nclass(6)
               Norder(6,j)=Nsym         ! sigma_d classes
            end do
            Norder(8,1)=1               ! i
            do j=1,Nclass(9)
               Norder(9,j)=2            ! S_{2n} classes
            end do
         else
            Nclass(2)=Nsym/2
            Nclass(3)=1
            Nclass(6)=1
            Nclass(9)=Nsym/2
            do j=1,Nclass(2)-1
               Norder(2,j)=2            ! Cn classes
            end do
            Norder(2,Nclass(2))=1
            do j=1,Nclass(3)
               Norder(3,j)=Nsym         ! C'2 classes
            end do
            do j=1,Nclass(6)
               Norder(6,j)=Nsym         ! sigma_d classes
            end do
            do j=1,Nclass(9)
               Norder(9,j)=2            ! S_{2n} classes
            end do
         end if
c-----------------------------------------------------------------------
      else if(Pgroup.eq."Sn") then
         No=Nsym         
         Nc=Nsym
         Nclass(2)=Nsym/2-1
         do j=1,Nclass(2)+1
            Norder(2,j)=1               ! C_{n/2} classes
         end do
         if (mod(Nsym/2,2).ne.0) then
            Nclass(8)=1                 ! note: i=Sn^(n/2)
            Nclass(9)=Nsym/2-1
            Norder(8,1)=1               ! i
         else
            Nclass(9)=Nsym/2
         end if
         do j=1,Nclass(9)
            Norder(9,j)=1               ! Sn classes
         end do
c-----------------------------------------------------------------------
      else if(Pgroup.eq."T") then ! still needs to be implemented
         Nc=4
         No=12
c-----------------------------------------------------------------------
      else if(Pgroup.eq."Td") then ! still needs to be implemented
         Nc=5
         No=24
c-----------------------------------------------------------------------
      else if(Pgroup.eq."Th") then ! still needs to be implemented
         Nc=8
         No=24
c-----------------------------------------------------------------------
      else if(Pgroup.eq."O") then ! still needs to be implemented
         Nc=5
         No=24
c-----------------------------------------------------------------------
      else if(Pgroup.eq."Oh") then ! still needs to be implemented
         Nc=10
         No=48
c-----------------------------------------------------------------------
      else if(Pgroup.eq."Ih") then ! still needs to be implemented
         Nc=10
         No=120
         continue
      else
         call errmsg("invalid group in get_class_order",.true.)
      end if

      end


c*********************************************************************** 
c     get character table of the point group
c
c     Author: Michael Kahnert
c     Date:   17.07.2004
c     Update: 16.04.2012
c
c***********************************************************************      
      subroutine get_characters(Pgroup,Nsym,maxNc,Nc,chi,chtabfile)      
      implicit none
c-----------------------------------------------------------------------
c     INPUT:
c     Pgroup:   character string representing the kind of point group 
c     Nsym:     integer characterising the main symmetry axis (for Cn,
c               Cnv, Cnh, Dnh, Dnd, and Sn groups)
c     Nc:       number of classes in the group
c     chtabfile: Name of file containing character table
c   
c     OUTPUT:
c     chi(i,j): character of the i-th irreducible representation (table row)
c               and the j-th conjugacy class (character table column)
c
c     NOTE:     The program assumes that the columns in the character 
c               tables are arranged in the following order:
c               1st column: class containing the unity element
c               next:       Cn^j classes in ascending order of j (if applicable)
c               next:       C'2 class (if applicable)
c               next:       C''2 class (if applicable)
c               next:       sigma_v class (if applicable)
c               next:       sigma_d class (if applicable)
c               next:       sigma_h class (if applicable)
c               next:       i class (if applicable)
c               next:       Sn^{j} classes in ascending order of j
c                           (if applicable)
c               This is the same order as that in the arrays Norder and 
c               Nc (see comments under get_class_order)
c               
c-----------------------------------------------------------------------
      character*3 Pgroup
      character*(*) chtabfile
      integer Nsym, maxNc, Nc, i,j
      complex*16 chi(maxNc,maxNc)
c-----------------------------------------------------------------------
c     Read character table computed with GAP and post-processed with
c     the scripts contained in the sub-directory GAP
c-----------------------------------------------------------------------      
      call get_gap_characters(maxNc,Nc,chi,chtabfile)
c-----------------------------------------------------------------------
c     need complex conjugate of characters for computing P matrix
c-----------------------------------------------------------------------
      do i=1,Nc
         do j=1,Nc
            chi(i,j)=conjg(chi(i,j))
         end do
      end do

      return
      end


c*********************************************************************** 
c     Read in characters in post-processed GAP format    
c     Author: Michael Kahnert
c*********************************************************************** 
      subroutine get_gap_characters(maxNc,Nc,chi,chtabfile)
      implicit none
c-----------------------------------------------------------------------
      integer maxNc, Nc
      complex*16 chi(maxNc,maxNc)
      character*(*) chtabfile

      integer linelen,linelen0
      parameter(linelen0=9999)
      character*(linelen0) line,line2
      integer err,i,j,k,Nc_gap,lnblnk,ii
      logical lread


      open(10,file=chtabfile,status='OLD',err=1000)

      lread=.true.
      do i=1,Nc
         read(10,'(a))',err=100,end=200) line
         linelen=lnblnk(line)
         if(linelen0.lt.linelen)then
            write(6,*)'FATAL ERROR: linelen0 < ',linelen
            stop
         endif
         call char2char(line,linelen,chi,i,Nc,maxNc)
      enddo

 200  close(10)

      return
 100  write(6,*)'ERROR reading file in get_character'
      stop
 1000 write(6,*)'ERROR opening character file'
      stop
      end



c*********************************************************************** 
c     Convert table entries from character strings to complex numbers
c
c     Author: Michael Kahnert
c     Date:   16.04.2012
c*********************************************************************** 
      subroutine char2char(line,linelen,chi,i,Nc,maxNc)
      implicit none

      integer linelen,nc,i,maxNc
      character*(*) line
      complex*16 chi(maxNc,maxNc)

      real*8 sgn,pi,eval,pval,iarr2real
      complex*16 Im,fac,val,val0
      integer j,k,nterms
      logical lexp,lpow,lfirstchar
      integer isnumber,po10,iarr(100)
      
      Im=cmplx(0d0,1d0)
      pi=acos(-1d0)
      fac=2d0*pi*Im

      j=1
      k=0
      sgn=1d0
      lexp=.false.
      lpow=.false.
      po10=-1
      val0=cmplx(0d0,0d0)
      lfirstchar=.true.

      do while (line(j:j).eq.' '.and.j.le.linelen) ! fast forward to 1st entry
         j=j+1
      enddo
      do while(j.le.linelen)
         if(line(j:j).eq.' ')then  ! new character entry
            do while (line(j+1:j+1).eq.' '.and.j+1.le.linelen) ! fast forward to next entry
               j=j+1
            enddo
            k=k+1
            if(lexp)then
               val=exp(fac/eval)
               if(lpow)then
                  pval=iarr2real(iarr,po10)
                  val=val**pval
               endif
               val=sgn*val
            else
               val=iarr2real(iarr,po10)
               val=sgn*val
            endif
            val0=val0+val
            chi(i,k)=val0
            sgn=1d0
            lexp=.false.
            lpow=.false.
            po10=-1
            lfirstchar=.true.
            val0=cmplx(0d0,0d0)
            nterms=1
         elseif(line(j:j).eq.'-'.and.lfirstchar)then
            sgn=-1d0
            lfirstchar=.false.
         elseif(line(j:j).eq.'E')then ! power of e-function
            j=j+1                     ! hop over '('
            lexp=.true.
            lfirstchar=.false.
         elseif(isnumber(line(j:j)).ge.0)then
            po10=po10+1               ! powers of 10
            iarr(po10+1)=isnumber(line(j:j))
            lfirstchar=.false.
         elseif(line(j:j).eq.')')then
            eval=iarr2real(iarr,po10)
            po10=-1
            lfirstchar=.false.
         elseif(line(j:j).eq.'^')then
            lpow=.true.
            lfirstchar=.false.
         elseif((line(j:j).eq.'+'.or.line(j:j).eq.'-')
     &           .and..not.lfirstchar)then ! new summation term
            if(lexp)then
               val=exp(fac/eval)
               if(lpow)then
                  pval=iarr2real(iarr,po10)
                  val=val**pval
               endif
               val=sgn*val
            else
               val=iarr2real(iarr,po10)
               val=sgn*val
            endif
            val0=val0+val
            lexp=.false.
            lpow=.false.
            po10=-1
            lfirstchar=.false.
            if(line(j:j).eq.'+')then
               sgn=1d0
            else
               sgn=-1d0
            endif
            nterms=nterms+1
         endif
         j=j+1
c     last entry:
         if(j.gt.linelen)then
            k=k+1
            if(lexp)then
               val=exp(fac/eval)
               if(lpow)then
                  pval=iarr2real(iarr,po10)
                  val=val**pval
               endif
               val=sgn*val
            else
               val=iarr2real(iarr,po10)
               val=sgn*val
            endif
            val0=val0+val
            chi(i,k)=val0
         endif
      enddo
      
      end




c*********************************************************************** 
      double precision function iarr2real(iarr,po10)
      implicit none
      integer iarr(100)
      integer po10,i,j
      real*8 val

      j=po10
      val=0d0
      do i=1,po10+1
         val=val+iarr(i)*10d0**j
         j=j-1
      enddo
      iarr2real=val

      end



c*********************************************************************** 
      integer function isnumber(c)
      implicit none

      character c
      integer n

      if(c.eq.'0')then
         n=0
      elseif(c.eq.'1')then
         n=1
      elseif(c.eq.'2')then
         n=2
      elseif(c.eq.'3')then
         n=3
      elseif(c.eq.'4')then
         n=4
      elseif(c.eq.'5')then
         n=5
      elseif(c.eq.'6')then
         n=6
      elseif(c.eq.'7')then
         n=7
      elseif(c.eq.'8')then
         n=8
      elseif(c.eq.'9')then
         n=9
      else
         n=-1
      endif

      isnumber=n
      end


c*********************************************************************** 
c     check if using correct character table
c     (under construction)
c
c     Author: Michael Kahnert
c     Date:   17.07.2004
c***********************************************************************      
      subroutine check_TabG(Pgroup,Nc,chtabfile)
      implicit none
c-----------------------------------------------------------------------
      character*3 Pgroup
      character*(*) chtabfile
      integer Nc,lines,lnblnk
      complex*16 dummy
      character cdummy
c-----------------------------------------------------------------------
      open(14,file=chtabfile)
      read(14,*) cdummy
      lines=0
      do while(.true.)
         read(14,*,end=100,err=200) dummy
         lines=lines+1
      enddo
 200  write(6,*)'Error reading character Table'
      close(14)
      return
 100  close(14)
      if(lines.ne.Nc)then
         write(6,*)"FATAL ERROR: wrong character table"
         write(*,1000) Pgroup,Nc
         write(*,2000) chtabfile(1:lnblnk(chtabfile)),lines
         stop
      endif

 1000 format(" Pgroup=",a," should have ",i4," classes")
 2000 format(" Group in table file ",a," has ",i4,
     &     " irreducible representations")
      return
      end


c*********************************************************************** 
c     Compute:
c     - number of times \alpha_\mu the \mu-th irreducible representation 
c       occurs in a reducible representation; 
c     - dimension of the \mu-th invariant subspace
c
c     Author: Michael Kahnert
c     Date:   16.04.2012
c*********************************************************************** 
      subroutine n_irr(Pgroup,dim,Nsym,No,Nc,Nclass,Norder,
     &     chi,maxNc,mmax,nmx,alpha,ir_size)
      implicit none
c-----------------------------------------------------------------------
c     Let D^{(red)} denote a reducible representation, and D^{(\mu)} 
c     denote the \mu-th irreducible representation. D^{(red)} can be 
c     written as a direct sum \oplus_{i=1}^Nc \alpha_\mu D^{(\mu)}.
c     Nc is the number of irreducible representations (which is equal
c     to the number of conjugacy classes in the group). The cardinal 
c     numbers \alpha_\mu can be computed as follows:
c
c     \alpha_\mu = No^{-1} *
c     \sum_{i=1}^{Nc} Norder(i) * chi^{(red)}(g_i) * \chi^{(\mu)\ast}(g_i)
c
c     where 
c     Norder(i):             number of group elements in class i
c     g_i:                   representative group element of class i
c     chi^{(red)}(g_i):      reducible characters
c     \chi^{(\mu)\ast}(g_i): complex conjugate of irreducible characters
c
c     From \alpha_\mu we immediately get the dimension of the 
c     \mu-th invariant subspace T_\mu:
c
c     dim(T_\mu) = \alpha_\mu * \chi^{(\mu)}(E),
c
c     where the character of the identity operation \chi^{(\mu)}(E) is
c     equal to the dimension of the \mu-th irreducible representation.
c
c     dim(T_\mu) is useful when constructing the P-matrix from the
c     projectors into the \mu-th invariant subspace; The number of
c     non-zero singular values of each projector should be equal
c     to dim(T_\mu)!!
c
c     Author: Michael Kahnert
c     Date:   12.04.2012
c-----------------------------------------------------------------------
      integer maxNc,ir_size(maxNc),alpha(maxNc)
      complex*16 chi(maxNc,maxNc),calpha,chired,csum
      character*3 Pgroup
      integer dim,mmax,nmx(-mmax:mmax)
      integer Nsym, No, Nc, Nclass(maxNc), Norder(maxNc,maxNc)
      real*8 dNsym,pi,epsilon
      integer i,j,k,mu,chidx,sum
      parameter ( pi=3.1415926535897932385d0,epsilon=1d-6)      

      ir_size(Nc+1)=0
      sum=0
      csum=cmplx(0d0,0d0)
      dNsym=dble(Nsym)
      do mu=1,Nc
         call chi_E(chired,dim)
         calpha=chired*chi(mu,1)
c-----------------------------------------------------------------------
c     Cs group:
c-----------------------------------------------------------------------
         if (Pgroup.eq."Cs") then
            call chi_sigma(chired)
            calpha=calpha+dble(Norder(7,1))*chired*chi(mu,2)
c-----------------------------------------------------------------------
c     Ci group:
c-----------------------------------------------------------------------
         elseif (Pgroup.eq."Ci") then
            call chi_inv(chired)
            calpha=calpha+dble(Norder(8,1))*chired*chi(mu,2)
c-----------------------------------------------------------------------
c     Cn group:
c-----------------------------------------------------------------------
         elseif (Pgroup.eq."Cn") then
            do j=1,Nclass(2)
               call chi_Cn(chired,dNsym,pi,nmx,mmax,j)
               calpha=calpha+dble(Norder(2,j))*chired*chi(mu,j+1)
            enddo
c-----------------------------------------------------------------------
c     Cnv group:
c-----------------------------------------------------------------------
         elseif (Pgroup.eq."Cnv") then
            do j=1,Nclass(2)
               call chi_Cn(chired,dNsym,pi,nmx,mmax,j)
               calpha=calpha+dble(Norder(2,j))*chired*chi(mu,j+1)
            enddo
c     no need to add contribution from sigma_v and sigma_d, since the
c     reducible characters of those are identically zero
c-----------------------------------------------------------------------
c     Dn group:
c-----------------------------------------------------------------------
         elseif (Pgroup.eq."Dn") then
            do j=1,Nclass(2)
               call chi_Cn(chired,dNsym,pi,nmx,mmax,j)
               calpha=calpha+dble(Norder(2,j))*chired*chi(mu,j+1)
            enddo
            chidx=Nclass(2)+1
            call chi_C2(chired,mmax,nmx)
            calpha=calpha+dble(Norder(3,1))*chired*chi(mu,chidx+1)
            if(mod(Nsym,2).eq.0)then
               calpha=calpha+dble(Norder(4,1))*chired*chi(mu,chidx+2)
            endif
c-----------------------------------------------------------------------
c     Cnh group:
c-----------------------------------------------------------------------
         elseif (Pgroup.eq."Cnh") then
            do j=1,Nclass(2)
               call chi_Cn(chired,dNsym,pi,nmx,mmax,j)
               calpha=calpha+dble(Norder(2,j))*chired*chi(mu,j+1)
            enddo
c     no need to add contribution from sigma_h, i, or Sn, since the
c     reducible characters of those are identically zero
c-----------------------------------------------------------------------
c     Dnh group:
c-----------------------------------------------------------------------
         elseif (Pgroup.eq."Dnh") then
            do j=1,Nclass(2)
               call chi_Cn(chired,dNsym,pi,nmx,mmax,j)
               calpha=calpha+dble(Norder(2,j))*chired*chi(mu,j+1)
            enddo
            chidx=Nclass(2)+1
            call chi_C2(chired,mmax,nmx)
            calpha=calpha+dble(Norder(3,1))*chired*chi(mu,chidx+1)
            if(mod(Nsym,2).eq.0)then
               calpha=calpha+dble(Norder(4,1))*chired*chi(mu,chidx+2)
            endif
c     no need to add contribution from sigma_h, sigma_v, sigma_d, i,
c     or Sn, since the reducible characters of those are identically zero
c-----------------------------------------------------------------------
c     Dnd group:
c-----------------------------------------------------------------------
         elseif (Pgroup.eq."Dnd") then
            do j=1,Nclass(2)
               call chi_Cn(chired,dNsym,pi,nmx,mmax,j)
               calpha=calpha+dble(Norder(2,j))*chired*chi(mu,j+1)
            enddo
            chidx=Nclass(2)+1
            call chi_C2(chired,mmax,nmx)
            calpha=calpha+dble(Norder(3,1))*chired*chi(mu,chidx+1)
            if(mod(Nsym,2).eq.0)
     &           calpha=calpha+dble(Norder(4,1))*chired*chi(mu,chidx+2)
c     no need to add contribution from sigma_d, i, or Sn, 
c     since the reducible characters of those are identically zero
c-----------------------------------------------------------------------
c     Sn group:
c-----------------------------------------------------------------------
         elseif (Pgroup.eq."Sn") then
            do j=1,Nclass(2)
               call chi_Cn(chired,dNsym/2d0,pi,nmx,mmax,j)
               calpha=calpha+dble(Norder(2,j))*chired*chi(mu,j+1)
            enddo
            chidx=Nclass(2)+1
c     no need to add contribution from i or Sn, 
c     since the reducible characters of those are identically zero
         endif
c-----------------------------------------------------------------------
         calpha=calpha/dble(No)
         alpha(mu)=int(calpha+epsilon)
         ir_size(mu)=int(chi(mu,1))*alpha(mu)
         sum=sum+ir_size(mu)
         csum=csum+calpha
         if(ir_size(mu).lt.0)then
            ir_size(Nc+1)=-1
            ir_size(mu)=0
         endif
      enddo
      if(ir_size(Nc+1).eq.-1)then
         call errmsg("dim(T_mu) < 0",.true.)
      endif
      if(sum.ne.dim)then
         ir_size(Nc+1)=-2
         write(6,*)'sum,csum,dim:',sum,csum,dim
         call errmsg("sum of dim(T_mu) is not equal to dim",.true.)
      endif

      end


c*********************************************************************** 
c     compute the transformation matrix P to transform the Q matrix
c     from the reducible into the irreducible basis (Qirr=P*Q*P^{-1})
c
c     Author: Michael Kahnert
c     Date:   31.07.2004
c     Update: 16.04.2012
c*********************************************************************** 
      subroutine get_Pmat(Pgroup,dim,Nsym,maxNc,sdim,
     &                    Nc,No,Nclass,Norder,mmax,nmmax,lmax,nmx,lsi,
     &                    chi,prec,ne,ir_size,Pmatmethod)
      use pmat
      implicit none
c-----------------------------------------------------------------------
c     INPUT:
c
c     Pgroup: character string characterising the point group
c     Nsym:   order of main symmetry axis (if applicable)
c     Nc:     number of classes in the group
c     No:     order of the group
c     Nclass: see comments under get_class_order
c     Norder: see comments under get_class_order
c     chi:    character table of the group
c     prec:   computational precision
c
c     OUTPUT:
c
c     P:      transformation matrix in compact format, where
c     ip:     P_{i,j}=P(k), i=ip(k), j=ip(k+ne), and 
c     ne:     number of non-zero elements of P (P is highly sparse!)
c     ir_size:dimensions of the invariant subspaces
c-----------------------------------------------------------------------
      character*3 Pgroup
      integer dim, Nsym, maxNc, sdim, Nc, No, 
     &        Nclass(maxNc), Norder(maxNc,maxNc), mmax, nmmax, lmax, 
     &        nmx(-mmax:mmax), lsi(-mmax:mmax,nmmax), ne,
     &        ir_size(maxNc),Pmatmethod

      complex*16 chi(maxNc,maxNc)
      real*8 prec
c-------------------------
      integer ne1, lirr, nrow,ncol, NumLinIndep,k,  i,j,err,sumdim,info
      real*8 rn
      double precision, dimension(dim) :: S
      logical lskip

c     Lapack arrays:
      integer lwork,lrwork
      complex*16, dimension(:), allocatable :: Work
      double precision, dimension(:), allocatable :: Rwork
      complex*16  dummy1(1,1),dummy2(1,1)

c     Propack arrays:
      integer lzwork,ldwork,liwork,neig,kmax
      integer, dimension(:), allocatable :: ip2
      double precision, dimension(:), allocatable :: bnd
      complex*16, dimension(:,:), allocatable :: U
      complex*16, dimension(:,:), allocatable :: V
      complex*16, dimension(:), allocatable :: zwork
      double precision, dimension(:), allocatable :: dwork
      double precision, dimension(4) :: doption
      integer, dimension(2) :: ioption
      integer, dimension(:), allocatable :: iwork
      real*8 prec2,eps
      double precision dlamch
      external dlamch, pp_matvecmult

c-----------------------------------------------------------------------
c     allocate arrays:
c-----------------------------------------------------------------------
c     Lapack:
      if(Pmatmethod.eq.2) then
         call logmsg('     SVD of projectors using LAPACK ')
         lwork=70*dim
         lrwork=5*dim
         allocate(work(lwork),stat=err)
         if(err.ne.0) goto 1130
         allocate(rwork(lrwork),stat=err)
         if(err.ne.0) goto 1130
         goto 1140
 1130    write(6,*)'ERROR: Allocation failure for LAPACK arrays'
         stop
 1140    continue
c     Propack:
      elseif(Pmatmethod.eq.3) then
         call logmsg('     SVD of projectors using PROPACK')
         prec2=1d-6
         kmax=dim
         lzwork=34*dim
         ldwork=34*dim+13*kmax+8*kmax**2+8
         liwork=8*kmax
         eps = dlamch('e')
         ioption(1) = 0
         ioption(2) = 1      
         doption(1) = sqrt(eps)
         doption(2) = eps**(3.0/4.0)
         doption(3) = 0d0
         doption(4) = 0.002
         allocate(ip2(sdim),stat=err)
         if(err.ne.0) goto 1110
         allocate(bnd(dim),stat=err)
         if(err.ne.0) goto 1110
         allocate(U(dim,dim+1),stat=err)
         if(err.ne.0) goto 1110
         allocate(V(dim,dim),stat=err)
         if(err.ne.0) goto 1110
         allocate(zwork(lzwork),stat=err)
         if(err.ne.0) goto 1110
         allocate(dwork(ldwork),stat=err)
         if(err.ne.0) goto 1110
         allocate(iwork(liwork),stat=err)
         goto 1120
 1110    write(6,*)'ERROR: Allocation failure for PROPACK arrays'
         stop
 1120    continue
         do i=1,dim
            do j=1,dim
               V(i,j)=cmplx(0d0,0d0)
            enddo
         enddo
      endif

c     P-matrix computation:
      allocate(P(sdim),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for P'
         stop
      endif
      allocate(ip(2*sdim),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for ip'
         stop
      endif
      allocate(jp(sdim),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for jp'
         stop
      endif
      allocate(ip1(sdim),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for ip1'
         stop
      endif
      allocate(jp1(sdim),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for jp1'
         stop
      endif
      allocate(P1(sdim),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for P1'
         stop
      endif
      allocate(P2(dim,dim),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for P2'
         stop
      endif
c-----------------------------------------------------------------------
c     loop over all irreducible representations of the group
c-----------------------------------------------------------------------
      ne=0
      nrow=0
      ncol=0
      sumdim=0
      lskip=.false.
      do lirr=1,Nc
         if(ir_size(lirr).eq.0)then
            lskip=.true. 
         else
            lskip=.false.
         endif

         if(.not.lskip)then
            call zero_par1(ne1,sdim)
c-----------------------------------------------------------------------
c     Projection matrix P1 into the lirr-th invariant subspace:
c-----------------------------------------------------------------------
            if (Pgroup.eq."Cs") then
               call get_Pmat_Cs(lirr,dim,Nsym,maxNc,sdim,
     &              mmax,nmmax,lmax,nmx,lsi,chi,prec, ne1)
            else if (Pgroup.eq."Ci") then
               call get_Pmat_Ci(lirr,dim,Nsym,maxNc,sdim,
     &              mmax,nmmax,lmax,nmx,lsi,chi,prec, ne1)
            else if (Pgroup.eq."Cn") then
               call get_Pmat_Cn(lirr,dim,Nsym,maxNc,sdim,
     &              mmax,nmmax,lmax,nmx,lsi,chi,prec,ne1)
            else if (Pgroup.eq."Cnv") then
               call  get_Pmat_Cnv(lirr,dim,Nsym,maxNc,sdim,
     &              mmax,nmmax,lmax,nmx,lsi,chi,prec, ne1)
            else if (Pgroup.eq."Dn") then
               call  get_Pmat_Dn(lirr,dim,Nsym,maxNc,sdim,
     &              mmax,nmmax,lmax,nmx,lsi,chi,prec, ne1)
            else if (Pgroup.eq."Cnh") then
               call  get_Pmat_Cnh(lirr,dim,Nsym,maxNc,sdim,
     &              mmax,nmmax,lmax,nmx,lsi,chi,prec,ne1)
            else if (Pgroup.eq."Dnh") then
               call  get_Pmat_Dnh(lirr,dim,Nsym,maxNc,sdim,
     &              mmax,nmmax,lmax,nmx,lsi,chi,prec, ne1)
            else if (Pgroup.eq."Sn") then
               call  get_Pmat_Sn(lirr,dim,Nsym,maxNc,sdim,
     &              mmax,nmmax,lmax,nmx,lsi,chi,prec,ne1)
            else
               call errmsg("invalid group in get_Pmat",.true.)
            end if

c-----------------------------------------------------------------------
c
c     Construction of the Transformation matrix P from the projection
c     matrices
c
c-----------------------------------------------------------------------
c     Option 1:
c     Determine the range of the projection matrix P1 by finding
c     the linearly independent row vectors.
c-----------------------------------------------------------------------
            if(Pmatmethod.eq.1)then
               call LinIndep(P1,sdim,ne1,ip1,jp1,prec,NumLinIndep)
c     Check that the number or linearly independent row vectors
c     of P1 is equal to the dimension of the invariant subspace:
               if(ir_size(lirr).ne.NumLinIndep) then
                  write(6,*)'ERROR for invariant subspace T_mu, mu=',
     &                 lirr
                  write(6,*)'  dim Range(P1)=',NumLinIndep
                  write(6,*)'  dim T_mu=',ir_size(lirr)
                  write(6,*)'Try a different choice of Pmatmethod!'
                  stop
               endif
c-----------------------------------------------------------------------
c     Options 2 and 3:
c     Singular value decomposition of projection matrix
c-----------------------------------------------------------------------
c     Option 2:
c     Call of SVD LAPACK routine (uses less work memory, is numerically
c     robust, but uses more CPU time than PROPACK):
c-----------------------------------------------------------------------
            elseif(Pmatmethod.eq.2) then
               do i=1,dim
                  do j=1,dim
                     P2(i,j)=cmplx(0d0,0d0)
                  enddo
               enddo
               do k=1,ne1
                  P2(ip1(k),jp1(k))=P1(k)
               enddo
               call zgesvd('N','O',dim,dim,P2,dim,S,dummy1,1,dummy2,
     &              1,Work,lwork,Rwork,info)
               if(info.ne.0)then
                  write(6,*) "zgesvd info=",info
                  call errmsg("ZGESVD failed",.true.)
               endif
c-----------------------------------------------------------------------
c     Option 3:
c     Call of PROPACK (this is faster than LAPACK for sparse matrices,
c     but it uses more work memory, and it seems to be numerically less
c     robust):
c-----------------------------------------------------------------------
            elseif(Pmatmethod.eq.3) then
c     store integer arrays ip1 and jp1 in ip2:
               do k=1,sdim+1
                  ip2(k)=0
               enddo
               ip2(1)=ne1
               do k=2,ne1+1
                  ip2(k)=ip1(k-1)
                  ip2(k+ne1)=jp1(k-1)
               enddo
c$$$               call zlansvd_irl('L','N','Y',dim,dim,dim,1,
c$$$     &              ir_size(lirr)+1,1000,pp_matvecmult,U,dim,S,bnd,
c$$$     &              V,dim,prec2,dwork,ldwork,zwork,lzwork,iwork,
c$$$     &              liwork,doption,ioption,info,P1,ip2)
               call zlansvd('N','Y',dim,dim,ir_size(lirr)+1,kmax,
     &              pp_matvecmult,U,dim,S,bnd,V,dim,prec2,dwork,ldwork,
     &              zwork,lzwork,iwork,liwork,doption,ioption,info,P1,
     &              ip2)
               if(info.ne.0)
     &              write(6,*)'WARNING: zlansvd returned info=',info
            endif
c-----------------------------------------------------------------------
c     Check the ratio of non-zero singular values 
c     S(i+1)/S(i), where i=ir_size(lirr)=dimension of lirr-th invariant
c     subspace. (The number of non-zero singular values should be equal 
c     to the dimension of the lirr-th invariant subspace. Those 
c     dimensions are stored in the array ir_size).
c-----------------------------------------------------------------------
            if(Pmatmethod.eq.2.or.Pmatmethod.eq.3) then
               i=1
               if (S(1).le.0d0)then
                  write(6,*)'lirr=',lirr,':'
                  call errmsg("No non-zero singular values found",
     &                 .true.)
               endif
               do i=1,ir_size(lirr)
                  if(S(i).le.0d0)then
                     write(6,*)'lirr=',lirr,':'
                     call errmsg
     &                 ("Too few non-zero singular values found",.true.)
                  endif
               enddo
               i=ir_size(lirr)
               if(S(i+1)/S(i).gt.1d0)then
                  write(6,*)'lirr=',lirr,':'
                  call errmsg("Too many singular values found",
     &                 .false.)               
               else
                  write(6,*)"         lirr=",lirr," :S(i+1)/S(i)=",
     &                 S(i+1)/S(i)
               endif
            endif
c-----------------------------------------------------------------------
c     Copy the right-singular vectors V^H belonging to the non-zero
c     singular values into P1 (^H = hermitian adjoint)
c-----------------------------------------------------------------------
            if(Pmatmethod.eq.2)then
               call zero_par1(ne1,sdim)
               k=0
               do i=1,ir_size(lirr)
                  do j=1,dim
                     if(abs(P2(i,j)).gt.prec)then
                        k=k+1
                        P1(k)=P2(i,j)
                        ip1(k)=i
                        jp1(k)=j
                     endif
                  enddo
               enddo
               ne1=k
            elseif(Pmatmethod.eq.3)then
               call zero_par1(ne1,sdim)
               k=0
               do i=1,ir_size(lirr)
                  do j=1,dim
                     if(abs(V(j,i)).gt.prec)then
                        k=k+1
                        P1(k)=conjg(V(j,i))
                        ip1(k)=i
                        jp1(k)=j
                     endif
                  enddo
               enddo
               ne1=k
            endif


 1000       format(18f5.1)
 1001       format(18i3)
c-----------------------------------------------------------------------
c     combine the transformation matrices P1 from different 
c     irreducible representations into one large matrix (which allows 
c     us to transform to the irreducible basis in one step)
c-----------------------------------------------------------------------
            if(ne+ne1.gt.sdim)then
               write(6,*) "ERROR: sdim too small",ne,ne1,sdim
               write(6,*)
     &              "NOTE: PROPACK should only be used for large Nsym!"
               stop
            endif
            call cpy_mat(P1,ne1,ip1,jp1,P,ip,jp,ir_size(lirr),nrow,ne)
         endif
      end do ! loop over irreducible representations/invariant subspaces

      do k=1,ne
         ip(k+ne)=jp(k)
      end do
c-----------------------------------------------------------------------
c     free memory
c-----------------------------------------------------------------------
      deallocate(P1,stat=err)
      deallocate(P2,stat=err)
      deallocate(jp,stat=err)
      deallocate(ip1,stat=err)
      deallocate(jp1,stat=err)
      if(Pmatmethod.eq.2)then
         deallocate(work,stat=err)
         deallocate(rwork,stat=err)
      elseif(Pmatmethod.eq.3)then
         deallocate(ip2,stat=err)
         deallocate(bnd,stat=err)
         deallocate(U,stat=err)
         deallocate(V,stat=err)
         deallocate(zwork,stat=err)
         deallocate(dwork,stat=err)
         deallocate(iwork,stat=err)
      endif
 
      return
      end


c*********************************************************************** 
c     compute the transformation matrix P1 for the 
c     lirr-th irreducible representation of the Cs group
c
c     Author: Michael Kahnert
c     Date:   31.07.2004
c*********************************************************************** 
      subroutine get_Pmat_Cs(lirr,dim,Nsym,maxNc,sdim,
     &        mmax,nmmax,lmax,nmx,lsi,chi,prec, ne1)
      use pmat
      implicit none
c-----------------------------------------------------------------------
      integer lirr, dim, Nsym, sdim, mmax, nmmax, lmax, 
     &        maxNc, nmx(-mmax:mmax), lsi(-mmax:mmax,nmmax), ne1
      complex*16 chi(maxNc,maxNc)
      real*8 prec
c----------------
      complex*16 cfac,cfac1
      real*8 rn,dNsym,pi,dj
      integer m, n, l, k, tau, tau1, mt
c-----------------------------------------------------------------------
      ne1=0
      k=0
      do tau=1,2
         tau1=tau-1
         do m=-mmax,mmax
            do n=max(1,abs(m)),nmx(m)
               l=lsi(m,n)+tau1*lmax
c-----------------------------------------------------------------------
c     E contribution to P1 matrix
c-----------------------------------------------------------------------
               cfac1=chi(lirr,1)
c-----------------------------------------------------------------------
c     sigma_h contribution to P1 matrix
c-----------------------------------------------------------------------
               call sigma_h(dj,dNsym,pi,m,n,tau,mt,cfac)
               cfac1=cfac1+chi(lirr,2)*cfac
c-----------------------------------------------------------------------
c     store non-zero contributions in the P1 matrix
c-----------------------------------------------------------------------
               call cnorm(cfac1,rn)
               if (rn.ge.prec) then
                  k=k+1
                  if(k.gt.sdim)then 
                     stop "ERROR: sdim too small"
                  endif
                  P1(k)=cfac1
                  ip1(k)=l
                  jp1(k)=l
               end if
            end do
         end do
      end do
      ne1=k

      return
      end


c*********************************************************************** 
c     compute the transformation matrix P1 for the 
c     lirr-th irreducible representation of the Ci group
c
c     Author: Michael Kahnert
c     Date:   31.07.2004
c*********************************************************************** 
      subroutine get_Pmat_Ci(lirr,dim,Nsym,maxNc,sdim,
     &        mmax,nmmax,lmax,nmx,lsi,chi,prec, ne1)
      use pmat
      implicit none
c-----------------------------------------------------------------------
      integer lirr, dim, Nsym, sdim, mmax, nmmax, lmax, 
     &        maxNc, nmx(-mmax:mmax), lsi(-mmax:mmax,nmmax),ne1
      complex*16 chi(maxNc,maxNc)
      real*8 prec
c----------------
      complex*16 cfac,cfac1
      real*8 rn,dNsym,pi,dj
      integer m, n, l, k, tau, tau1, mt 
c-----------------------------------------------------------------------
      ne1=0
      k=0
      do tau=1,2
         tau1=tau-1
         do m=-mmax,mmax
            do n=max(1,abs(m)),nmx(m)
               l=lsi(m,n)+tau1*lmax
c-----------------------------------------------------------------------
c     E contribution to P1 matrix
c-----------------------------------------------------------------------
               cfac1=chi(lirr,1)
c-----------------------------------------------------------------------
c     i contribution to P1 matrix
c-----------------------------------------------------------------------
               call inv(dj,dNsym,pi,m,n,tau,mt,cfac)
               cfac1=cfac1+chi(lirr,2)*cfac
c-----------------------------------------------------------------------
c     store non-zro contributions in the P1 matrix
c-----------------------------------------------------------------------
               call cnorm(cfac1,rn)
               if (rn.ge.prec) then
                  k=k+1
                  if(k.gt.sdim)then 
                     stop "ERROR: sdim too small"
                  endif
                  P1(k)=cfac1
                  ip1(k)=l
                  jp1(k)=l
               end if
            end do
         end do
      end do
      ne1=k

      return
      end


c*********************************************************************** 
c     compute the transformation matrix P1 for the 
c     lirr-th irreducible representation of the Cn groups
c
c     Author: Michael Kahnert
c     Date:   31.07.2004
c*********************************************************************** 
      subroutine get_Pmat_Cn(lirr,dim,Nsym,maxNc,sdim,
     &        mmax,nmmax,lmax,nmx,lsi,chi,prec,ne1)
      use pmat
      implicit none
c-----------------------------------------------------------------------
      integer lirr, dim, Nsym, sdim, mmax, nmmax, lmax, nmx(-mmax:mmax),
     &        lsi(-mmax:mmax,nmmax), ne1,maxNc
      complex*16 chi(maxNc,maxNc)
      real*8 prec
c----------------------
      complex*16 cfac,cfac1
      real*8 rn,dNsym,pi,dj
      integer m, n, l, k, tau, tau1, mt, j
c-----------------------------------------------------------------------
      ne1=0
      k=0
      dNsym=dble(Nsym)
      pi=dacos(-1.d0)

      do tau=1,2
         tau1=tau-1
         do m=-mmax,mmax
            do n=max(1,abs(m)),nmx(m)
               l=lsi(m,n)+tau1*lmax
c-----------------------------------------------------------------------
c     E contribution to P1 matrix
c-----------------------------------------------------------------------
               cfac1=chi(lirr,1)
               do j=1,Nsym-1
                  dj=dble(j)
c-----------------------------------------------------------------------
c     Cn contributions to P1 matrix
c-----------------------------------------------------------------------
                  call Cn(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac1=cfac1+chi(lirr,j+1)*cfac
               end do
c-----------------------------------------------------------------------
c     store non-zero contributions in the P1 matrix
c-----------------------------------------------------------------------
               call cnorm(cfac1,rn)
               if (rn.ge.prec) then
                  k=k+1
                  if(k.gt.sdim)then 
                     stop "ERROR: sdim too small"
                  endif
                  P1(k)=cfac1
                  ip1(k)=l
                  jp1(k)=l
               end if
            end do
         end do
      end do
      ne1=k

      return
      end



c*********************************************************************** 
c     compute the transformation matrix P1 for the 
c     lirr-th irreducible representation of the Cnh groups
c
c     Author: Michael Kahnert
c     Date:   31.07.2004
c*********************************************************************** 
      subroutine get_Pmat_Cnh(lirr,dim,Nsym,maxNc,sdim,
     &        mmax,nmmax,lmax,nmx,lsi,chi,prec, ne1)
      use pmat
      implicit none
c-----------------------------------------------------------------------
      integer lirr, dim, Nsym, sdim, mmax, nmmax, lmax, nmx(-mmax:mmax),
     &        lsi(-mmax:mmax,nmmax), ne1,maxNc
      complex*16 chi(maxNc,maxNc)
      real*8 prec
c----------------------
      complex*16 cfac,cfac1
      real*8 rn,dNsym,pi,dj
      integer m, n, l, k, tau, tau1, mt, j, chidx
      logical Nodd
c-----------------------------------------------------------------------
      ne1=0
      k=0
      dNsym=dble(Nsym)
      pi=dacos(-1.d0)
      if (mod(Nsym,2).eq.0) then 
         Nodd=.false.
      else
         Nodd=.true.
      end if

      do tau=1,2
         tau1=tau-1
         do m=-mmax,mmax
            do n=max(1,abs(m)),nmx(m)
               l=lsi(m,n)+tau1*lmax
c-----------------------------------------------------------------------
c     E contribution to P1 matrix
c-----------------------------------------------------------------------
               cfac1=chi(lirr,1)
               chidx=1
c-----------------------------------------------------------------------
c     Cn contributions to P1 matrix
c-----------------------------------------------------------------------
               do j=1,Nsym-1
                  chidx=chidx+1
                  dj=dble(j)
                  call Cn(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
               end do
c-----------------------------------------------------------------------
c     sigma_h contribution to P1 matrix
c-----------------------------------------------------------------------
               chidx=chidx+1
               call sigma_h(dj,dNsym,pi,m,n,tau,mt,cfac)
               cfac1=cfac1+chi(lirr,chidx)*cfac
c-----------------------------------------------------------------------
c     i contribution to P1 matrix (only for Nsym even)
c-----------------------------------------------------------------------
               if (.not. Nodd) then
                  chidx=chidx+1
                  call inv(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
               end if
c-----------------------------------------------------------------------
c     Sn contributions to P1 matrix
c-----------------------------------------------------------------------
               do j=1,Nsym-1
c-----------------------------------------------------------------------
c     for even Nsym, jump over S_Nsym^{Nsym/2}, since it is identical
c     with the inversion i and thus already accounted for!
c-----------------------------------------------------------------------
                  if (.not.Nodd .and. j.eq.(Nsym/2)) goto 2100
                  chidx=chidx+1
                  dj=dble(j)
                  call Sn(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
 2100             continue
               end do
c-----------------------------------------------------------------------
c     store non-zero contributions in the P1 matrix
c-----------------------------------------------------------------------
               call cnorm(cfac1,rn)
               if (rn.ge.prec) then
                  k=k+1
                  if(k.gt.sdim)then 
                     stop "ERROR: sdim too small"
                  endif
                  P1(k)=cfac1
                  ip1(k)=l
                  jp1(k)=l
               end if
            end do
         end do
      end do
      ne1=k

      return
      end



c*********************************************************************** 
c     compute the transformation matrix P1 for the 
c     lirr-th irreducible representation of the Cnv groups
c
c     Author: Michael Kahnert
c     Date:   31.07.2004
c*********************************************************************** 
      subroutine get_Pmat_Cnv(lirr,dim,Nsym,maxNc,sdim,
     &        mmax,nmmax,lmax,nmx,lsi,chi,prec, ne1)
      use pmat
      implicit none
c-----------------------------------------------------------------------
      integer lirr, dim, Nsym, sdim, mmax, nmmax, lmax, nmx(-mmax:mmax),
     &        lsi(-mmax:mmax,nmmax), ne1,maxNc
      complex*16 chi(maxNc,maxNc)
      real*8 prec
c----------------------
      complex*16 cfac,cfac1,cfac2
      real*8 rn1,rn2,dNsym,pi,dj
      integer m,n,l,lx,k,tau,tau1, mt, j,jcmax,jsmax, chidx
      logical Nodd
c-----------------------------------------------------------------------
      ne1=0
      k=0
      dNsym=dble(Nsym)
      pi=dacos(-1.d0)
c-----------------------------------------------------------------------
c     maximum superscripts of Cn^j and sigma_v^{j} symmetries
c-----------------------------------------------------------------------
      if (mod(Nsym,2).eq.0) then
         jcmax=Nsym/2
         jsmax=jcmax-1
         Nodd=.false.
      else
         jcmax=(Nsym-1)/2
         jsmax=Nsym-1
         Nodd=.true.
      end if

      do tau=1,2
         tau1=tau-1
         do m=-mmax,mmax
            do n=max(1,abs(m)),nmx(m)
c-----------------------------------------------------------------------
c     E contribution to P1 matrix
c-----------------------------------------------------------------------
               cfac1=chi(lirr,1)
               cfac2=(0.d0,0.d0)
               chidx=1
c-----------------------------------------------------------------------
c     Cn contributions to P1 matrix
c-----------------------------------------------------------------------
               do j=1,jcmax
                  chidx=chidx+1
                  dj=dble(j)
                  call Cn(dj,dNsym,pi,m,n,tau,mt,cfac)
c-----------------------------------------------------------------------
c     add the contributions from both elements in each Cn^{j} class:
c     they are complex conjugate to each other!
c-----------------------------------------------------------------------
                  if (Nodd.or.(j.lt.jcmax)) 
     &               cfac=cmplx(2d0*dreal(cfac),0d0)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
               end do
c-----------------------------------------------------------------------
c     sigma_v contributions to P1 matrix
c-----------------------------------------------------------------------
               chidx=chidx+1
               do j=0,jsmax
                  dj=dble(j)
                  call sigma_v(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac2=cfac2+chi(lirr,chidx)*cfac
              end do
c-----------------------------------------------------------------------
c     sigma_d contributions to P1 matrix (only for even Nsym)
c-----------------------------------------------------------------------
             if (.not.Nodd) then
                chidx=chidx+1
                do j=0,jsmax
                   dj=dble(j)
                  call sigma_d(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac2=cfac2+chi(lirr,chidx)*cfac
                 end do
              end if
c-----------------------------------------------------------------------
c     store non-zero contributions in the correct order in the P1 matrix
c-----------------------------------------------------------------------
               call cnorm(cfac1,rn1)
               call cnorm(cfac2,rn2)
               l=lsi(m,n)+tau1*lmax
               lx=lsi(-m,n)+tau1*lmax
               if (rn1.ge.prec.and.rn2.ge.prec) then
                  if (l.lt.lx) then 
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac1
                     ip1(k)=l
                     jp1(k)=l
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac2
                     ip1(k)=l
                     jp1(k)=lx
                  else if (l.gt.lx) then
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac2
                     ip1(k)=l
                     jp1(k)=lx
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac1
                     ip1(k)=l
                     jp1(k)=l
                  else 
                     call cnorm(cfac1+cfac2,rn1)
                     if (rn1.ge.prec) then
                        k=k+1
                        if(k.gt.sdim)then 
                           stop "ERROR: sdim too small"
                        endif
                        P1(k)=cfac1+cfac2
                        ip1(k)=l
                        jp1(k)=l
                     end if
                  end if
               else if (rn1.ge.prec) then
                  k=k+1
                  if(k.gt.sdim)then 
                     stop "ERROR: sdim too small"
                  endif
                  P1(k)=cfac1
                  ip1(k)=l
                  jp1(k)=l
               else if (rn2.ge.prec) then
                  k=k+1
                  if(k.gt.sdim)then 
                     stop "ERROR: sdim too small"
                  endif
                  P1(k)=cfac2
                  ip1(k)=l
                  jp1(k)=lx
              end if
            end do
         end do
      end do
      ne1=k

      return
      end




c*********************************************************************** 
c     compute the transformation matrix P1 for the 
c     lirr-th irreducible representation of the Dn groups
c
c     Author: Michael Kahnert
c     Date:   31.07.2004
c*********************************************************************** 
      subroutine get_Pmat_Dn(lirr,dim,Nsym,maxNc,sdim,
     &        mmax,nmmax,lmax,nmx,lsi,chi,prec,ne1)
      use pmat
      implicit none
c-----------------------------------------------------------------------
      integer lirr, dim, Nsym, sdim, mmax, nmmax, lmax, nmx(-mmax:mmax),
     &        lsi(-mmax:mmax,nmmax), ne1,maxNc
      complex*16 chi(maxNc,maxNc)
      real*8 prec
c----------------------
      complex*16 cfac,cfac1,cfac2
      real*8 rn1,rn2,dNsym,pi,dj
      integer m,n,l,lx,k,tau,tau1, mt, j,jcmax,j2max,chidx
      logical Nodd
c-----------------------------------------------------------------------
      ne1=0
      k=0
      dNsym=dble(Nsym)
      pi=dacos(-1.d0)
c-----------------------------------------------------------------------
c     maximum superscripts of Cn^j and C'2^j, C''2^j symmetries
c-----------------------------------------------------------------------
      if (mod(Nsym,2).eq.0) then
         jcmax=Nsym/2
         j2max=jcmax-1
         Nodd=.false.
      else
         jcmax=(Nsym-1)/2
         j2max=Nsym-1
         Nodd=.true.
      end if

      do tau=1,2
         tau1=tau-1
         do m=-mmax,mmax
            do n=max(1,abs(m)),nmx(m)
c-----------------------------------------------------------------------
c     E contribution to P1 matrix
c-----------------------------------------------------------------------
               cfac1=chi(lirr,1)
               cfac2=(0.d0,0.d0)
               chidx=1
c-----------------------------------------------------------------------
c     Cn contributions to P1 matrix
c-----------------------------------------------------------------------
               do j=1,jcmax
                  chidx=chidx+1
                  dj=dble(j)
                  call Cn(dj,dNsym,pi,m,n,tau,mt,cfac)
c-----------------------------------------------------------------------
c     add the contributions from both elements in each Cn^{j} class:
c     they are complex conjugate to each other!
c-----------------------------------------------------------------------
                  if (Nodd.or.(j.lt.jcmax)) 
     &               cfac=cmplx(2d0*dreal(cfac),0d0)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
               end do
c-----------------------------------------------------------------------
c     C'2 contributions to P1 matrix
c-----------------------------------------------------------------------
               chidx=chidx+1
               do j=0,j2max
                  dj=dble(j)
                  call C2p(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac2=cfac2+chi(lirr,chidx)*cfac
              end do
c-----------------------------------------------------------------------
c     C''2 contributions to P1 matrix (only for even Nsym)
c-----------------------------------------------------------------------
             if (.not.Nodd) then
                chidx=chidx+1
                do j=0,j2max
                   dj=dble(j)
                  call C2pp(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac2=cfac2+chi(lirr,chidx)*cfac
                 end do
              end if
c-----------------------------------------------------------------------
c     store non-zero contributions in the correct order in the P1 matrix
c-----------------------------------------------------------------------
               call cnorm(cfac1,rn1)
               call cnorm(cfac2,rn2)
               l=lsi(m,n)+tau1*lmax
               lx=lsi(-m,n)+tau1*lmax
               if (rn1.ge.prec.and.rn2.ge.prec) then
                  if (l.lt.lx) then 
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac1
                     ip1(k)=l
                     jp1(k)=l
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac2
                     ip1(k)=l
                     jp1(k)=lx
                  else if (l.gt.lx) then
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac2
                     ip1(k)=l
                     jp1(k)=lx
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac1
                     ip1(k)=l
                     jp1(k)=l
                  else 
                     call cnorm(cfac1+cfac2,rn1)
                     if (rn1.ge.prec) then
                        k=k+1
                        if(k.gt.sdim)then 
                           stop "ERROR: sdim too small"
                        endif
                        P1(k)=cfac1+cfac2
                        ip1(k)=l
                        jp1(k)=l
                     end if
                  end if
               else if (rn1.ge.prec) then
                  k=k+1
                  if(k.gt.sdim)then 
                     stop "ERROR: sdim too small"
                  endif
                  P1(k)=cfac1
                  ip1(k)=l
                  jp1(k)=l
               else if (rn2.ge.prec) then
                  k=k+1
                  if(k.gt.sdim)then 
                     stop "ERROR: sdim too small"
                  endif
                  P1(k)=cfac2
                  ip1(k)=l
                  jp1(k)=lx
              end if
            end do
         end do
      end do
      ne1=k

      return
      end



c*********************************************************************** 
c     compute the transformation matrix P1 for the 
c     lirr-th irreducible representation of the Dnh groups
c
c     Author: Michael Kahnert
c     Date:   31.07.2004
c*********************************************************************** 
      subroutine get_Pmat_Dnh(lirr,dim,Nsym,maxNc,sdim,
     &        mmax,nmmax,lmax,nmx,lsi,chi,prec,ne1)
      use pmat
      implicit none
c-----------------------------------------------------------------------
      integer lirr, dim, Nsym, sdim, mmax, nmmax, lmax, nmx(-mmax:mmax),
     &        lsi(-mmax:mmax,nmmax), ne1,maxNc
      complex*16 chi(maxNc,maxNc)
      real*8 prec
c----------------------
      complex*16 cfac,cfac1,cfac2
      real*8 rn1,rn2,dNsym,pi,dj
      integer m,n,l,lx,k,tau,tau1, mt, j,jcmax,j2max,jsmax,chidx
      logical Nodd
c-----------------------------------------------------------------------
      ne1=0
      k=0
      dNsym=dble(Nsym)
      pi=dacos(-1.d0)
c-----------------------------------------------------------------------
c     maximum superscript jcmax of Cn^j symmetry j2max of 
c     C'2^j, C''2^j, sigma_v and sigma_d symmetries, and jsmax of 
c     Sn symmetry
c-----------------------------------------------------------------------
      if (mod(Nsym,2).eq.0) then
         jcmax=Nsym/2
         j2max=jcmax-1
         jsmax=j2max
         Nodd=.false.
      else
         jcmax=(Nsym-1)/2
         j2max=Nsym-1
         jsmax=jcmax
         Nodd=.true.
      end if

      do tau=1,2
         tau1=tau-1
         do m=-mmax,mmax
            do n=max(1,abs(m)),nmx(m)
c-----------------------------------------------------------------------
c     E contribution to P1 matrix
c-----------------------------------------------------------------------
               cfac1=chi(lirr,1)
               cfac2=(0.d0,0.d0)
               chidx=1
c-----------------------------------------------------------------------
c     Cn contributions to P1 matrix
c-----------------------------------------------------------------------
               do j=1,jcmax
                  chidx=chidx+1
                  dj=dble(j)
                  call Cn(dj,dNsym,pi,m,n,tau,mt,cfac)
c-----------------------------------------------------------------------
c     add the contributions from both elements in each Cn^{j} class:
c     they are complex conjugate to each other!
c-----------------------------------------------------------------------
                  if (Nodd.or.(j.lt.jcmax)) 
     &               cfac=cmplx(2d0*dreal(cfac),0d0)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
               end do
c-----------------------------------------------------------------------
c     C'2 contributions to P1 matrix
c-----------------------------------------------------------------------
               chidx=chidx+1
               do j=0,j2max
                  dj=dble(j)
                  call C2p(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac2=cfac2+chi(lirr,chidx)*cfac
              end do
c-----------------------------------------------------------------------
c     C''2 contributions to P1 matrix (only for even Nsym)
c-----------------------------------------------------------------------
             if (.not.Nodd) then
                chidx=chidx+1
                do j=0,j2max
                   dj=dble(j)
                  call C2pp(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac2=cfac2+chi(lirr,chidx)*cfac
                 end do
              end if
c-----------------------------------------------------------------------
c     sigma_v contributions to P1 matrix
c-----------------------------------------------------------------------
               chidx=chidx+1
               do j=0,j2max
                  dj=dble(j)
                  call sigma_v(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac2=cfac2+chi(lirr,chidx)*cfac
              end do
c-----------------------------------------------------------------------
c     sigma_d contributions to P1 matrix (only for even Nsym)
c-----------------------------------------------------------------------
             if (.not.Nodd) then
                chidx=chidx+1
                do j=0,j2max
                   dj=dble(j)
                  call sigma_d(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac2=cfac2+chi(lirr,chidx)*cfac
                 end do
              end if
c-----------------------------------------------------------------------
c     sigma_h contribution to P1 matrix
c-----------------------------------------------------------------------
               chidx=chidx+1
               call sigma_h(dj,dNsym,pi,m,n,tau,mt,cfac)
               cfac1=cfac1+chi(lirr,chidx)*cfac
c-----------------------------------------------------------------------
c     i contribution to P1 matrix (only for Nsym even)
c-----------------------------------------------------------------------
               if (.not. Nodd) then
                  chidx=chidx+1
                  call inv(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
               end if
c-----------------------------------------------------------------------
c     Sn contributions to P1 matrix
c-----------------------------------------------------------------------
               do j=1,jsmax
                  chidx=chidx+1
                  dj=dble(j)
                  call Sn(dj,dNsym,pi,m,n,tau,mt,cfac)
c-----------------------------------------------------------------------
c     add the contributions from both elements in each Sn^{j} class:
c     they are complex conjugate to each other!
c-----------------------------------------------------------------------
                  cfac=cmplx(2d0*dreal(cfac),0d0)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
               end do
c-----------------------------------------------------------------------
c     store non-zero contributions in the correct order in the P1 matrix
c-----------------------------------------------------------------------
               call cnorm(cfac1,rn1)
               call cnorm(cfac2,rn2)
               l=lsi(m,n)+tau1*lmax
               lx=lsi(-m,n)+tau1*lmax
               if (rn1.ge.prec.and.rn2.ge.prec) then
                  if (l.lt.lx) then 
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac1
                     ip1(k)=l
                     jp1(k)=l
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac2
                     ip1(k)=l
                     jp1(k)=lx
                  else if (l.gt.lx) then
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac2
                     ip1(k)=l
                     jp1(k)=lx
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac1
                     ip1(k)=l
                     jp1(k)=l
                  else 
                     call cnorm(cfac1+cfac2,rn1)
                     if (rn1.ge.prec) then
                        k=k+1
                        if(k.gt.sdim)then 
                           stop "ERROR: sdim too small"
                        endif
                        P1(k)=cfac1+cfac2
                        ip1(k)=l
                        jp1(k)=l
                     end if
                  end if
               else if (rn1.ge.prec) then
                  k=k+1
                  if(k.gt.sdim)then 
                     stop "ERROR: sdim too small"
                  endif
                  P1(k)=cfac1
                  ip1(k)=l
                  jp1(k)=l
               else if (rn2.ge.prec) then
                  k=k+1
                  if(k.gt.sdim)then 
                     stop "ERROR: sdim too small"
                  endif
                  P1(k)=cfac2
                  ip1(k)=l
                  jp1(k)=lx
              end if
            end do
         end do
      end do
      ne1=k

      return
      end



c*********************************************************************** 
c     compute the transformation matrix P1 for the 
c     lirr-th irreducible representation of the Dnd groups
c     NOTE: This subroutine has not been tested yet, since Tsym does
c     not currently have any geometries with Dnd symmetry!
c
c     Author: Michael Kahnert
c     Date:   12.04.2012
c*********************************************************************** 
      subroutine get_Pmat_Dnd(lirr,dim,Nsym,maxNc,sdim,
     &        mmax,nmmax,lmax,nmx,lsi,chi,prec, ne1)
      use pmat
      implicit none
c-----------------------------------------------------------------------
      integer lirr, dim, Nsym, sdim, mmax, nmmax, lmax, nmx(-mmax:mmax),
     &        lsi(-mmax:mmax,nmmax), ne1,maxNc
      complex*16 chi(maxNc,maxNc)
      real*8 prec
c----------------------
      complex*16 cfac,cfac1,cfac2
      real*8 rn1,rn2,dNsym,pi,dj
      integer m,n,l,lx,k,tau,tau1, mt, j,jcmax,j2max,jsmax,chidx
      logical Nodd
c-----------------------------------------------------------------------
      ne1=0
      k=0
      dNsym=dble(Nsym)
      pi=dacos(-1.d0)
c-----------------------------------------------------------------------
c     maximum superscript jcmax of Cn^j symmetry j2max of 
c     C'2^j and sigma_d symmetries, and jsmax of 
c     Sn symmetry
c-----------------------------------------------------------------------
      j2max=Nsym-1
      if (mod(Nsym,2).eq.0) then
         jcmax=Nsym/2
         jsmax=Nsym-1
         Nodd=.false.
      else
         jcmax=(Nsym-1)/2
         jsmax=Nsym-2
         Nodd=.true.
      end if

      do tau=1,2
         tau1=tau-1
         do m=-mmax,mmax
            do n=max(1,abs(m)),nmx(m)
c-----------------------------------------------------------------------
c     E contribution to P1 matrix
c-----------------------------------------------------------------------
               cfac1=chi(lirr,1)
               cfac2=(0.d0,0.d0)
               chidx=1
c-----------------------------------------------------------------------
c     Cn contributions to P1 matrix
c-----------------------------------------------------------------------
               do j=1,jcmax
                  chidx=chidx+1
                  dj=dble(j)
                  call Cn(dj,dNsym,pi,m,n,tau,mt,cfac)
c-----------------------------------------------------------------------
c     add the contributions from both elements in each Cn^{j} class:
c     they are complex conjugate to each other!
c-----------------------------------------------------------------------
                  if (Nodd.or.(j.lt.jcmax)) 
     &               cfac=cmplx(2d0*dreal(cfac),0d0)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
               end do
c-----------------------------------------------------------------------
c     C'2 contributions to P1 matrix
c-----------------------------------------------------------------------
               chidx=chidx+1
               do j=0,j2max
                  dj=dble(j)
                  call C2p(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac2=cfac2+chi(lirr,chidx)*cfac
              end do
c-----------------------------------------------------------------------
c     sigma_d contributions to P1 matrix (only for even Nsym)
c-----------------------------------------------------------------------
              chidx=chidx+1
              do j=0,j2max
                 dj=dble(j)
                 call sigma_d(dj,dNsym,pi,m,n,tau,mt,cfac)
                 cfac2=cfac2+chi(lirr,chidx)*cfac
              end do
c-----------------------------------------------------------------------
c     sigma_h contribution to P1 matrix
c-----------------------------------------------------------------------
               chidx=chidx+1
               call sigma_h(dj,dNsym,pi,m,n,tau,mt,cfac)
               cfac1=cfac1+chi(lirr,chidx)*cfac
c-----------------------------------------------------------------------
c     i contribution to P1 matrix (only for Nsym even)
c-----------------------------------------------------------------------
               if (Nodd) then
                  chidx=chidx+1
                  call inv(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
               end if
c-----------------------------------------------------------------------
c     Sn contributions to P1 matrix
c-----------------------------------------------------------------------
               do j=1,jsmax,2
                  chidx=chidx+1
                  dj=dble(j)
                  call Sn(dj,2d0*dNsym,pi,m,n,tau,mt,cfac)
c-----------------------------------------------------------------------
c     add the contributions from both elements in each Sn^{j} class:
c     they are complex conjugate to each other!
c-----------------------------------------------------------------------
                  cfac=cmplx(2d0*dreal(cfac),0d0)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
               end do
c-----------------------------------------------------------------------
c     store non-zero contributions in the correct order in the P1 matrix
c-----------------------------------------------------------------------
               call cnorm(cfac1,rn1)
               call cnorm(cfac2,rn2)
               l=lsi(m,n)+tau1*lmax
               lx=lsi(-m,n)+tau1*lmax
               if (rn1.ge.prec.and.rn2.ge.prec) then
                  if (l.lt.lx) then 
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac1
                     ip1(k)=l
                     jp1(k)=l
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac2
                     ip1(k)=l
                     jp1(k)=lx
                  else if (l.gt.lx) then
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac2
                     ip1(k)=l
                     jp1(k)=lx
                     k=k+1
                     if(k.gt.sdim)then 
                        stop "ERROR: sdim too small"
                     endif
                     P1(k)=cfac1
                     ip1(k)=l
                     jp1(k)=l
                  else 
                     call cnorm(cfac1+cfac2,rn1)
                     if (rn1.ge.prec) then
                        k=k+1
                        if(k.gt.sdim)then 
                           stop "ERROR: sdim too small"
                        endif
                        P1(k)=cfac1+cfac2
                        ip1(k)=l
                        jp1(k)=l
                     end if
                  end if
               else if (rn1.ge.prec) then
                  k=k+1
                  if(k.gt.sdim)then 
                     stop "ERROR: sdim too small"
                  endif
                  P1(k)=cfac1
                  ip1(k)=l
                  jp1(k)=l
               else if (rn2.ge.prec) then
                  k=k+1
                  if(k.gt.sdim)then 
                     stop "ERROR: sdim too small"
                  endif
                  P1(k)=cfac2
                  ip1(k)=l
                  jp1(k)=lx
              end if
            end do
         end do
      end do
      ne1=k

      return
      end



c*********************************************************************** 
c     compute the transformation matrix P1 for the 
c     lirr-th irreducible representation of the Sn groups
c
c     Author: Michael Kahnert
c     Date:   31.07.2004
c*********************************************************************** 
      subroutine get_Pmat_Sn(lirr,dim,Nsym,maxNc,sdim,
     &        mmax,nmmax,lmax,nmx,lsi,chi,prec, ne1)
      use pmat
      implicit none
c-----------------------------------------------------------------------
      integer lirr, dim, Nsym, sdim, mmax, nmmax, lmax, nmx(-mmax:mmax),
     &        lsi(-mmax:mmax,nmmax), ne1,maxNc
      complex*16 chi(maxNc,maxNc)
      real*8 prec
c----------------------
      complex*16 cfac,cfac1
      real*8 rn,dNsym,pi,dj,dcidx
      integer m, n, l, k, tau, tau1, mt, j, chidx
      logical Nodd
c-----------------------------------------------------------------------
      ne1=0
      k=0
      dNsym=dble(Nsym)
      dcidx=dNsym/2d0
      pi=dacos(-1.d0)
      if (mod(Nsym/2,2).eq.0) then 
         Nodd=.false.
      else
         Nodd=.true.
      end if

      do tau=1,2
         tau1=tau-1
         do m=-mmax,mmax
            do n=max(1,abs(m)),nmx(m)
               l=lsi(m,n)+tau1*lmax
c-----------------------------------------------------------------------
c     E contribution to P1 matrix
c-----------------------------------------------------------------------
               cfac1=chi(lirr,1)
               chidx=1
c-----------------------------------------------------------------------
c     Cn contributions to P1 matrix
c-----------------------------------------------------------------------
               do j=1,Nsym/2-1
                  chidx=chidx+1
                  dj=dble(j)
                  call Cn(dj,dcidx,pi,m,n,tau,mt,cfac)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
               end do
c-----------------------------------------------------------------------
c     i contribution to P1 matrix (only for Nsym even)
c-----------------------------------------------------------------------
               if (Nodd) then
                  chidx=chidx+1
                  call inv(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
               end if
c-----------------------------------------------------------------------
c     Sn contributions to P1 matrix
c-----------------------------------------------------------------------
               do j=1,Nsym-1,2
c-----------------------------------------------------------------------
c     for odd Nsym/2, jump over S_Nsym^{Nsym/2}, since it is identical
c     with the inversion i and thus already accounted for!
c-----------------------------------------------------------------------
                  if (Nodd .and. j.eq.(Nsym/2)) goto 2100
                  chidx=chidx+1
                  dj=dble(j)
                  call Sn(dj,dNsym,pi,m,n,tau,mt,cfac)
                  cfac1=cfac1+chi(lirr,chidx)*cfac
 2100             continue
               end do
c-----------------------------------------------------------------------
c     store non-zero contributions in the P1 matrix
c-----------------------------------------------------------------------
               call cnorm(cfac1,rn)
               if (rn.ge.prec) then
                  k=k+1
                  if(k.gt.sdim)then 
                     stop "ERROR: sdim too small"
                  endif
                  P1(k)=cfac1
                  ip1(k)=l
                  jp1(k)=l
               end if
            end do
         end do
      end do
      ne1=k

      return
      end



c*********************************************************************** 
c     set some parameters to zero
c
c     Author: Michael Kahnert
c     Date:   31.07.2004
c*********************************************************************** 
      subroutine zero_par1(ne1,sdim)
      use pmat
      implicit none
c-----------------------------------------------------------------------
      integer sdim,ne1,k
      complex*16 c0
      parameter (c0=(0.d0,0.d0))
c-----------------------------------------------------------------------
      ne1=0
      do k=1,sdim
         P1(k)=c0
         ip1(k)=0
         jp1(k)=0
      end do

      return
      end



c*********************************************************************** 
c     set some parameters to zero
c
c     Author: Michael Kahnert
c     Date:   31.07.2004
c*********************************************************************** 
      subroutine zero_par2(sdim)
      use pmat
      implicit none
c-----------------------------------------------------------------------
      integer sdim,k
      complex*16 c0
      parameter (c0=(0.d0,0.d0))
c-----------------------------------------------------------------------
      do k=1,sdim
         Pinv(k)=c0
      end do
      do k=1,2*sdim
         ipinv(k)=c0
      end do

      return
      end



c***********************************************************************
c     copy matrix P1(1..ne1) (in compact vector storage format)
c     into matrix P(nstart,...,nstart+ne1)
c
c     Author: Michael Kahnert
c     24.07.2004
c
c     Input:
c
c     P1:  matrix in compact storage scheme, such that 
c          p1_{i,j}=P1(k), i=ip1(k), j=jp1(k), k=1,...,ne1
c     ne1: size of P1
c     ip1: see P1
c     jp1: see P1
c     ir_size: number of conjugacy classes =
c              number of linearly independent row vectors in P1
c
c     Input/Output:
c
c     P:   matrix in compact storage scheme, such that (on input)
c          p_{i,j}=P(k), i=ip(k), j=jp(k), k=1,...,ne;
c          on output, the matrix P1 is appended to P
c     ne:  on input, dimension of the original matrix P;
c          on output, dimension of the new matrix P, i.e. ne=ne+ne1
c     nrow:number of non-zero row vectors of P
c***********************************************************************
      subroutine cpy_mat(P1,ne1,ip1,jp1,P,ip,jp,
     &                   ir_size,nrow,ne)
      implicit none
      integer ne1,ne,k,l
      complex*16 P1(ne1),P(ne+ne1)
      integer ip1(ne1), jp1(ne1), ip(ne+ne1), jp(ne+ne1),
     &        ir_size,nrow

      do k=1,ne1
         l=ne+k
         P(l)=P1(k)
         ip(l)=ip1(k)+nrow
         jp(l)=jp1(k)
      end do
      ne=ne+ne1
      nrow=nrow+ir_size

      return
      end


c***********************************************************************
c     complex norm:
c
c     Author: Michael Kahnert
c     24.07.2004
c***********************************************************************
      subroutine cnorm(z,r)
      implicit none
      complex*16 z
      real*8 r

      r=dsqrt(dreal(z)**2+dimag(z)**2)

      return
      end


c***********************************************************************
c     Invert the P-matrix:
c     The P-matrix is unitary, so the inverse is the hermetian adjoint.
c***********************************************************************
      subroutine invP(n,ne,sdim,neinv)
      use pmat
      implicit none
c-----------------------------------------------------------------------
c
c     Author: Michael Kahnert
c     Date:   25.07.2004
c     Update: 07.05.2012
c
c     Input: (some of it in module pmat)
c
c     n:  size of matrices P, P^{-1}
c     ne: number of non-zero elements of matrix P
c     P:  transformation matrix P
c     ip: contains indices of P in the same way as ipinv for Pinv
c     
c
c     Output: (some of it in module pmat)
c
c     Pinv:  contains inverse of transformation matrix P^{-1} 
c            in compact storage scheme, such that 
c            P^{-1}_{i,j}=Pinv(k), i=ipinv(k), j=ipinv(ne+k), 
c            k=1,...,neinv
c     ipinv: see Pinv
c     neinv: number of non-zero elements of Pinv
c-----------------------------------------------------------------------

      integer n,ne,neinv,sdim
      integer dim,nelu
      integer i,j,k,err

c-----------------------------------------------------------------------
c     allocate arrays:
c-----------------------------------------------------------------------
      allocate(Pinv(sdim),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for Pinv'
         stop
      endif
      allocate(ipinv(2*sdim),stat=err)
      if(err.ne.0)then
         write(6,*)'FATAL ERROR: allocation failed for ipinv'
         stop
      endif
      call zero_par2(sdim)

      do k=1,ne
         i=ip(k)
         j=ip(k+ne)
         ipinv(k)=j
         ipinv(k+ne)=i
         Pinv(k)=conjg(P(k))
      end do
      neinv=ne
      

      return
      end



c***********************************************************************
c     forward transformation into irreducible basis: Qirr = P*Q*P^{-1}:
c
c     Author: Michael Kahnert
c     25.07.2004
c
c     Input:
c
c     P:     contains transformation matrix P in compact storage scheme, 
c            such that 
c            P_{i,j}=P(k), i=ip(k), j=ip(ne+k), k=1,...,ne
c     ne:    number of non-zero elements of P
c     neinv: number of non-zero elements of Pinv
c     ip:    see P
c     Pinv:  inverse of P
c     ipinv: as ip, but for Pinv
c     n:     size of matrices P, Q, P^{-1}
c
c     Input/Output:
c
c     Q:  on input, Q-matrix in the reducible basis of vector spherical
c         functions. On output, Qirr in the irreducible basis
c***********************************************************************
      subroutine Q2Qirr(n,ne,neinv,Q)
      use pmat
      implicit none

      integer n,ne,neinv
      complex*16 Q(n,n), ctmp(n)
      integer i,j,k

c-----------------------------------------------------------------------
c     compute P*Q and store the result in Q:
c-----------------------------------------------------------------------
      do j=1,n
         do i=1,n
            ctmp(i)=cmplx(0d0,0d0)
         end do
         do k=1,ne
            ctmp(ip(k)) = ctmp(ip(k))+P(k)*Q(ip(k+ne),j)
         end do
         do i=1,n
            Q(i,j)=ctmp(i)
         end do
      end do
c-----------------------------------------------------------------------
c     compute P*Q*P^{-1} 
c-----------------------------------------------------------------------
      do i=1,n
         do j=1,n
            ctmp(j)=cmplx(0d0,0d0)
         end do
         do k=1,neinv
            ctmp(ipinv(k+neinv)) = ctmp(ipinv(k+neinv))+
     &                             Q(i,ipinv(k))*Pinv(k)
         end do
         do j=1,n
            Q(i,j)=ctmp(j)
         end do
      end do

      return
      end



c*********************************************************************** 
c
c     matrix-vector multiplication routine used in propack SVD
c
c     Author: Michael Kahnert
c*********************************************************************** 
      subroutine pp_matvecmult(transa,m,n,x,y,P1,ip2)
      implicit none
      
      character*1 transa
      integer m,n,ip2(*)
      complex*16 x(*),y(*),P1(*)
      integer i,j,k,ne1
      complex*16 aij


      ne1=ip2(1)
      if(transa.eq.'N'.or.transa.eq.'n')then
         do i=1,m
            y(i) = cmplx(0d0,0d0)
         enddo
         do k=1,ne1
            i=ip2(k+1)
            j=ip2(k+1+ne1)
            aij = P1(k)
            y(i) = y(i) + aij*x(j)
         enddo
      else
         do i=1,n
            y(i) = cmplx(0d0,0d0)
         enddo
         do k=1,ne1
            i=ip2(k+1)
            j=ip2(k+1+ne1)
            aij= conjg(P1(k))
            y(j) = y(j) + x(i)*aij
         enddo
      endif

      end



c***********************************************************************
c     back transformation into reducible basis: 
c     Q^{-1} = P^{-1}*Qirr^{-1}*P:
c
c     Author: Michael Kahnert
c     25.07.2004
c
c     Input:
c
c     P:     contains transformation matrix P in compact storage scheme, 
c            such that 
c            P_{i,j}=P(k), i=ip(k), j=ip(ne+k), k=1,...,ne
c     ne:    number of non-zero elements of P
c     neinv: number of non-zero elements of Pinv
c     ip:    see P
c     Pinv:  inverse of P
c     ipinv: as ip, but for Pinv
c     n:     size of matrices P, Q, P^{-1}
c
c     Input/Output:
c
c     Q:  on input, Qirr^{-1}-matrix in the irreducible 
c         basis of vector spherical functions. 
c         On output, Q^{-1} in the reducible basis
c***********************************************************************
      subroutine Qirr2Q(n,ne,neinv,Q)
      use pmat
      implicit none

      integer n,ne,neinv
      complex*16 Q(n,n), ctmp(n)
      integer i,j,k,err

c-----------------------------------------------------------------------
c     compute P^{-1}*Qirr^{-1} and store the result in Q:
c-----------------------------------------------------------------------
      do j=1,n
         do i=1,n
            ctmp(i)=cmplx(0d0,0d0)
         end do
         do k=1,neinv
            ctmp(ipinv(k)) = ctmp(ipinv(k))+Pinv(k)*Q(ipinv(k+neinv),j)
         end do
         do i=1,n
            Q(i,j)=ctmp(i)
         end do
      end do
c-----------------------------------------------------------------------
c     compute P^{-1}*Qirr^{-1}*P
c-----------------------------------------------------------------------
      do i=1,n
         do j=1,n
            ctmp(j)=cmplx(0d0,0d0)
         end do
         do k=1,ne
            ctmp(ip(k+ne)) = ctmp(ip(k+ne))+Q(i,ip(k))*P(k)
         end do
        do j=1,n
           Q(i,j)=ctmp(j)
        end do
      end do

      deallocate(P,stat=err)
      if(err.ne.0)then
         write(6,*)"Deallocation error for array P"
      endif
      deallocate(Pinv,stat=err)
      if(err.ne.0)then
         write(6,*)"Deallocation error for array Pinv"
      endif
      deallocate(ip,stat=err)
      if(err.ne.0)then
         write(6,*)"Deallocation error for array ip"
      endif
      deallocate(ipinv,stat=err)
      if(err.ne.0)then
         write(6,*)"Deallocation error for array ipinv"
      endif
      

      return
      end



c***********************************************************************
c     check that Qirr is block diagonal
c     Author: Michael Kahnert
c***********************************************************************
      subroutine check_Qirr(Q,ir_size,dim,maxNc,Nc,prec,lmax3)
      implicit none
      integer dim,maxNc
      complex*16 Q(dim,dim)
      real*8 prec
      integer ir_size(maxNc),Nc,lmax3
      integer i,j,l,bstart,bend
      real*8 rn,rnmx
      logical bd

      bd=.true.
      rnmx=0d0
      bstart=1
      l=1
      do while (ir_size(l).eq.0)
         l=l+1
      enddo
      bend=bstart+ir_size(l)-1
      l=1
      do i=1,dim
         do j=1,dim
            if (i.ge.bstart.and.i.le.bend.and.j.ge.bstart.and.
     &           j.le.bend) then
               continue
            else
               call cnorm(Q(i,j),rn)
               if (rn.ge.prec) then
                  bd=.false.
                  rnmx=max(rn,rnmx)
               end if
            end if
         end do
         if (i.eq.bend) then
            l=l+1
            do while (ir_size(l).eq.0)
               l=l+1
            enddo            
            bstart=bend+1
            bend=bstart+ir_size(l)-1
         end if
      end do

      if (bd) then
         write(6,*) "         Qirr in block-diagonal form!"
      else
         write(6,*)"WARNING: off-block diagonal elements found in Qirr"
         write(6,*)"with max norm |Q| =",rnmx,', precision = ',prec
      end if

      lmax3=1
      do l=1,Nc
         lmax3=max(lmax3,ir_size(l))
      end do

      return
      end



c***********************************************************************
c     subroutine to determine linearly independent row vectors of A
c
c     Author: Michael Kahnert
c     24.07.2004
c
c     Input:
c
c     sdim:  max size of P1
c     prec:  precision
c
c     Input/Output:
c
c     ne:    number of nonzero matrix elements
c     ip1,jp1: integer arrays: if matrix P1_{i,j}.ne.ZERO then 
c            ip1(k)=i, jp1(k)=j, k=1,...,ne
c     p1:     p1(k) contains matrix element P1_{i,j}, where i=ip1(k), j=jp1(k)
c
c     on output the linearly dependent row vectors are eliminated and
c     the linearly independent row vectors are shifted up
c
c     Output:
c
c     NumLinIndep:
c            number of linearly independent row vectors
c
c***********************************************************************
      subroutine LinIndep(P1,sdim,ne,ip1,jp1,prec,NumLinIndep)
      implicit none
c-----------------------------------------------------------------------
      integer sdim,ne
      complex*16 P1(sdim)
      real*8 prec
      integer ip1(sdim),jp1(sdim),NumLinIndep
c--------------------------------------
      integer i,j,k,iu,nel,nel1,nenew
      logical compare,lindep
      real*8 fac1,fac2
      complex*16 c0, c1
      parameter (c0=(0.d0,0.d0), c1=(1.d0,0.d0))

      NumLinIndep=0
c-----------------------------------------------------------------------
c     loop over 1st row vector
c-----------------------------------------------------------------------
      iu=0
      i=1
      do while (i.le.ne)
c-----------------------------------------------------------------------
c        number of elements in the row vector:
c-----------------------------------------------------------------------
         call get_nel(ip1(i),sdim-i,ne,nel)
c-----------------------------------------------------------------------
c        shift row vector up if there are nullvectors between it and the
c        nearest upper row vector:
c-----------------------------------------------------------------------
         if (ip1(i).gt.iu) then
            iu=iu+1
            call shift_row_vector(sdim,ne,ip1,i,iu,nel)
         end if
c-----------------------------------------------------------------------
c        loop over all row vectors below the current row vector:
c-----------------------------------------------------------------------
         j=i+nel
         do while(j.le.ne)
c-----------------------------------------------------------------------
c           only need to compare those vectors that have an equal  
c           number of elements:
c-----------------------------------------------------------------------
            compare=.false.
            call get_nel(ip1(j),sdim-j,ne,nel1)
            if (nel.eq.nel1) then
               compare=.true.
c-----------------------------------------------------------------------
c           only need to compare those vectors that have elements at 
c           the same position:
c-----------------------------------------------------------------------
               do k=0,nel-1
                  if (jp1(i+k).ne.jp1(j+k)) then
                     compare=.false.
                     goto 1100
                  end if
               end do
            end if
 1100       continue
c-----------------------------------------------------------------------
c           if vectors have elements at identical column indices,
c           check for linear dependence:
c-----------------------------------------------------------------------
            if (compare) then
               call check_lin_dep(P1,sdim,ne,ip1,jp1,i,j,nel,lindep,
     &                            prec)
c-----------------------------------------------------------------------
c              delete linearly dependent row vectors:
c-----------------------------------------------------------------------
               if (lindep) then
                  call delete_row(P1,sdim,ne,ip1,jp1,j,nel1)
               end if
            end if
            j=j+nel1
         end do
         i=i+nel
         NumLinIndep=NumLinIndep+1
      end do

      return
      end


c***********************************************************************
c     shift up row vector from position i to position iu
c
c     Author: Michael Kahnert
c     24.07.2004
c***********************************************************************
      subroutine shift_row_vector(sdim,ne,ip1,i,iu,nel)
      implicit none
      integer ne,sdim
      integer ip1(sdim),i,iu,nel,j

      do j=1,nel
         ip1(i+j-1)=iu
      end do

      return
      end


c***********************************************************************
c     check if two row vectors are linearly dependent
c
c     Author: Michael Kahnert
c     24.07.2004
c***********************************************************************
      subroutine check_lin_dep(P1,sdim,ne,ip1,jp1,i,j,nel,lindep,prec)
      implicit none
      integer ne,sdim
      complex*16 P1(sdim)
      real*8 prec
      integer ip1(sdim),jp1(sdim),i,j,nel
      logical lindep

      complex*16 cfac0,cfac1
      real*8 rn
      integer k
      
      lindep=.true.
      cfac0=P1(i)/P1(j)
      do k=1,nel-1
         cfac1=P1(i+k)/P1(j+k)
         call cnorm(cfac1-cfac0,rn)
         if (rn.ge.prec) then
            lindep=.false.
            return
         end if
      end do

      return
      end


c***********************************************************************
c     delete row vector j from the matrix P1
c
c     Author: Michael Kahnert
c     24.07.2004
c***********************************************************************
      subroutine delete_row(P1,sdim,ne,ip1,jp1,j,nel1)
      implicit none
      integer ne,sdim
      complex*16 P1(sdim)
      integer ip1(sdim),jp1(sdim),j,nel1,i,k
      
      do i=j+nel1,ne
         k=i-nel1
         ip1(k)=ip1(i)-1
         jp1(k)=jp1(i)
         P1(k)=P1(i)
      end do
      ne=ne-nel1

      j=j-nel1

      return
      end


c***********************************************************************
c     number of elements in row vector (complex storage scheme)
c
c     Author: Michael Kahnert
c     24.07.2004
c***********************************************************************
      subroutine get_nel(ip1,sdimi,ne,nel)
      implicit none
      integer sdimi,ne
      integer ip1(sdimi),nel

      nel=0
      do while(ip1(1).eq.ip1(1+nel))
         nel=nel+1
      end do
      return
      end


c=======================================================================
c     PART IV :  P A R T I C L E   G E O M E T R I E S : 
c=======================================================================


c***********************************************************************
c     Gauss-Legendre Quadrature points for a polyhedral cylinder,
c     using a special quadrature scheme for cylindrical particles)
c
c     Author: Michael Kahnert
c***********************************************************************
      subroutine GL_pts_cyl(eps,alpha,pi,ng1,ng2,
     &     ng1dim,ng2dim,SIGMAh,DIHEDRAL,       
     &     ngauss_per_int,th_nint,phi_nint,
     &     phipts,N,th1,wth1,th2,wth2,phi,wphi,phic1,phic2,ngstart)
     
      implicit none
      integer phipts,N,i,j,ng1,ng2,ngj,
     &        ng1dim,ng2dim,ngstart,                                            
     &        ngauss_per_int,th_nint,phi_nint,iint,k
      real*8  th1(ng1dim), wth1(ng1dim),
     &        th2(ng2dim,phipts), wth2(ng2dim,phipts),
     &        phi(phipts),wphi(phipts),phic1(ng1dim),
     &        phic2(ng1dim),x1(ng1),w1(ng1),
     *        x2(ng2),w2(ng2),
     &        x3(phipts),w3(phipts),rho,
     &        eps,alpha,pi,pih,pi2,theta0,dn,
     &        sinphi,cosphi,sinal,phi0,r1,                                               
     &        x0(ngauss_per_int),w0(ngauss_per_int),fac,xold
      logical SIGMAh,DIHEDRAL
     
c=======================================================================
c     Determination of polar and azimuthal quadrature points
c     (the polar quadrature points xth1 lie on the top and
c     bottom facet of the cylinder, the xth2 lie on the 
c     sides of the cylinder):
c=======================================================================
      dn = dble(N)
c-----------------------------------------------------------------------         
c     Double quadrature in phi and x=cos(theta) for polygonial
c     cylinders: 
c-----------------------------------------------------------------------   
      call gauss(ngauss_per_int,1,0,x0,w0)
      xold=0.
      fac=2d0*pi/dn/dble(phi_nint)
      do iint=1,phi_nint
c------------------------------------------------------------
c     The quadrature points 0 < x(i) < 1 and their weights  
c     w(i) are redistributed onto the azimuth-interval      
c     [0,2*pi/N/phi_nint], as follows:
c      phi(i) = 2*pi/N/phi_nint * x(i) + phi1
c      wphi(i) = 2*pi/N/phi_nint * w(i)
c     where phi1 is the offset point, i.e. phi(max) from the
c     preceeding sub-interval
c------------------------------------------------------------
         do i=1,ngauss_per_int
            k=(iint-1)*ngauss_per_int+i
            phi(k)  = fac * x0(i) + xold             
            wphi(k) = fac * w0(i)
         end do
         xold=xold+fac
      enddo 

      ng1 = int(dble(th_nint)*eps/(2.d0*(1.d0+eps)))               
      ng2 = max(th_nint-2*ng1,0)                                   
      if(ng2.eq.0)                                                    
     &     write(6,*)'WARNING: ng2=0 in GL_pts'                       
c------------------------------------------------------------         
c     For PRSIMS, for a given azimuth phi, the interval [0,pi] of     
c     theta is split into the partial intervals                       
c     [0,theta0], [theta0,pi-theta0], [pi-theta0,pi]. The             
c     interval [0,theta0] is split into ng1 quadrature intervals,     
c     where ng1 = arctan(eps * (cos(phi)-cot(alpha)*sin(phi))/        
c     (1-sin(phi)**2/cos(alpha)**2), and 0<=phi<pi/N.                 
c     The quadrature points 0 < x(i) < 1 and their weights            
c     w(i) are redistributed onto the polar-subintervals              
c     as follows:                                                     
c      theta(i)    = theta0/ng1 * x(i)  + theta1                      
c      theta(N1+i) = (pi-2*theta0)/ng2*x(i) + theta1                  
c      wtheta(i)   = theta0/ng1*w(i)                                  
c      wtheta(N1+i)= (pi-2*theta0)/ng2*w(i),                          
c     where theta1 is the offset-point, i.e. theta(max) from          
c     the preceeding subinterval.                                     
c                                                                     
c     The number of quadrature intervals in the                       
c     subintervals [0,theta0], [theta0,pi-theta0] and                 
c     [pi-theta0,pi] is ng1, ng2 and ng1, where these numbers
c     are chosen such that ng1 = eps/2 * ng2, 2*ng1+ng2=thpts/th_nint.
c------------------------------------------------------------        
      pih = pi*0.5d0                                         
      sinal  = dsin(alpha)                                   
      pi2 = 2.d0*pi                                          
      
      xold=0.   
      theta0=datan(eps)
      fac=theta0/dble(ng1)                                
      do iint=1,ng1                                       
         do j=1,ngauss_per_int                            
            k=(iint-1)*ngauss_per_int+j                   
            wth1(k)=fac * w0(j)                       
            th1(k) =fac * x0(j) + xold                
         enddo                                            
         xold=xold+fac                                    
      enddo                                               

      do i=1,phipts        
         phi0 = phi(i)                                       
         sinphi = dsin(phi0)                                 
         cosphi = sqrt(1.d0-sinphi**2)                       
         theta0 = datan(eps * rho(alpha,sinphi,cosphi,sinal)) 
         fac=(pi-2*theta0)/dble(ng2) 
         xold=theta0   
         do iint=1,ng2                                       
            do j=1,ngauss_per_int                            
               k=(iint-1)*ngauss_per_int+j
               if(k.gt.ng2dim) goto 1000
               wth2(k,i)=fac * w0(j)                       
               th2(k,i) =fac * x0(j) + xold  
            end do                                           
            xold=xold+fac    
         enddo  
 1000    continue
      end do        
                                         
      if (.not.SIGMAh) then
         theta0=datan(eps)
         fac=theta0/dble(ng1)  
         xold=pi-theta0
         do iint=1,ng1                                         
            do j=1,ngauss_per_int                              
               k=(ng1+iint-1)*ngauss_per_int+j             
               wth1(k)=fac*w0(j)                           
               th1(k)=fac*x0(j) + xold                     
            enddo                                              
            xold=xold+fac                                      
         enddo    
      endif                                             
c=======================================================================
c       The azimuthal surface integration can be carried out
c       analytically on the top and bottom faces. The phi-
c       integration goes over the full azimuthal range for 
c       those polar angles that lie inside a circle of
c       maximum radius on the top and bottom faces. For those
c       polar points lying outside this circle but still on
c       the top (or bottom) face, the phi-integration only
c       extends up to a critical azimuth phic:
c=======================================================================

c-----------------------------------------------------------------------
c       determine the polar angle theta0 of a circle with
c       maximum radius lying entirely inside the top facet of
c       the particle:
c-----------------------------------------------------------------------
        sinphi = dsin(pi/dble(N))
        cosphi = dcos(pi/dble(N))
        r1  = rho(alpha,sinphi,cosphi,sinal)
        theta0 = datan(eps*r1)
        ng1=ng1*ngauss_per_int
        ng2=ng2*ngauss_per_int
c-----------------------------------------------------------------------
c       determine the index ngstart, starting from which the
c       polar quadrature points xth1 lie outside the circle
c       characterized by the angle th1:
c-----------------------------------------------------------------------
        ngstart=ng1+1
        do i=1,ng1
           if (th1(i).gt.theta0) goto 1155
        end do
 1155   ngstart = i
c-----------------------------------------------------------------------
c       critical azimuthal angles phic for each polar quadrature
c       angle on the top (and bottom) facet of the particle:
c-----------------------------------------------------------------------
        if (ngstart.le.ng1) then
           call phi_crit(eps,dcos(alpha),sinal,th1,pi,ngstart,
     &           ng1,phic1)
        else
           write(6,*)'>>>WARNING: low GL-point density on top facet'
           write(6,*)'th0,th1:',th1(ng1),theta0
        end if
c-----------------------------------------------------------------------
c       The following lines merely allow to test the code for
c       cylinders with polygonal cross sections with the
c       settings DIHEDRAL=.false. and/or SIGMAh=.false.
c       However, for actual particles w/o sigma_h symmetry
c       the critical azimuths phic1, phic2 on the bottom facet
c       would be different from those on the top facet. Analogously,
c       for a particle w/o dihedral symmetry, the critical angles
c       phic1 in the interval [0,pi/N] are different from the phic2
c       in the interval [pi/N,2*pi/N]. Therefore, if the user wants to
c       use the code for particles w/o sigma_h or dihedral symmetry,
c       the following lines need to be changed accordingly.
c-----------------------------------------------------------------------
        if (.not.DIHEDRAL) then
           do i=ngstart,ng1
             phic2(i) = phic1(i)
           end do
        end if
        if (.not.SIGMAh) then
           do i=ngstart,ng1
              j=2*ng1-i+1
              phic1(j)=phic1(i)
              phic2(j)=phic2(i)
           end do
        end if

      return
      end
      
      
c***********************************************************************
c     Gauss-Legendre Quadrature points, using
c     a double-quadrature scheme:
c
c     Author: Michael Kahnert
c***********************************************************************
      subroutine GL_pts(Geom,AXISYM,eps,alpha,pi,thpts,
     &     ngauss_per_int,th_nint,phi_nint,thpts2,
     &     phipts,N,theta,wtheta,phi,wphi,ng1,ng2)
     
      implicit none
      integer thpts,thpts2,phipts,N,i,j,ng1,ng2,ngj,
     &        ngauss_per_int,th_nint,phi_nint,iint,k
      real*8  theta(thpts,phipts), wtheta(thpts,phipts),
     &        phi(phipts),wphi(phipts),
     &        x1(thpts),w1(thpts),
     *        x2(thpts),w2(thpts),
     &        x3(phipts),w3(phipts),rho,
     &        eps,alpha,pi,pih,pi2,theta0,dn,
     &        sinphi,cosphi,sinal,phi0,xold,
     &        x0(ngauss_per_int),w0(ngauss_per_int),fac
      character*6 Geom
      logical AXISYM
     
      dn = dble(N)
c================================================================
c     phi-quadrature:
c================================================================
      if(AXISYM)then
         call gauss(2,1,0,x0,w0)
         fac=pi/dn
         phi(1)= fac*x0(1)
         phi(2)= fac*x1(1)
         wphi(1)  = fac*w0(1)
         wphi(2)  = fac*w0(2)
c     second call to gauss, needed later for theta-points:
         call gauss(ngauss_per_int,1,0,x0,w0)
      else
         call gauss(ngauss_per_int,1,0,x0,w0)
c------------------------------------------------------------
c     Double quadrature in phi and x=cos(theta) for polygonial
c     cylinders: 
c------------------------------------------------------------
         xold=0.
         fac=2d0*pi/dn/dble(phi_nint)
         do iint=1,phi_nint
c------------------------------------------------------------
c     The quadrature points 0 < x(i) < 1 and their weights
c     w(i) are redistributed onto the azimuth-interval
c     [0,2*pi/N/phi_nint], as follows:
c     phi(i) = 2*pi/N/phi_nint * x(i) + phi1
c     wphi(i) = 2*pi/N/phi_nint * w(i)
c
c     where phi1 is the offset point, i.e. phi(max) from the 
c     preceeding sub-interval
c------------------------------------------------------------
            do i=1,ngauss_per_int
               k=(iint-1)*ngauss_per_int+i
               phi(k)  = fac * x0(i) + xold
               wphi(k) = fac * w0(i)
            end do
            xold=xold+fac
         enddo
      endif

c================================================================
c     theta-quadrature:
c================================================================
      if(Geom.eq.'PRISMS')then
         ng1 = int(dble(th_nint)*eps/(2.d0*(1.d0+eps)))
         ng2 = max(th_nint-2*ng1,0)
      else
         ng1=0
         ng2=thpts
      endif
      if(ng2.eq.0)
     &     write(6,*)'WARNING: ng2=0 in GL_pts'
c------------------------------------------------------------
c     For PRSIMS, for a given azimuth phi, the interval [0,pi] of 
c     theta is split into the partial intervals
c     [0,theta0], [theta0,pi-theta0], [pi-theta0,pi]. The
c     interval [0,theta0] is split into ng1 quadrature intervals,
c     where ng1 = arctan(eps * (cos(phi)-cot(alpha)*sin(phi))/
c     (1-sin(phi)**2/cos(alpha)**2), and 0<=phi<pi/N.  
c     The quadrature points 0 < x(i) < 1 and their weights
c     w(i) are redistributed onto the polar-subintervals 
c     as follows:
c      theta(i)    = theta0/ng1 * x(i)  + theta1
c      theta(N1+i) = (pi-2*theta0)/ng2*x(i) + theta1
c      wtheta(i)   = theta0/ng1*w(i)
c      wtheta(N1+i)= (pi-2*theta0)/ng2*w(i),
c     where theta1 is the offset-point, i.e. theta(max) from
c     the preceeding subinterval.
c     
c     The number of quadrature intervals in the
c     subintervals [0,theta0], [theta0,pi-theta0] and
c     [pi-theta0,pi] is ng1, ng2 and ng1, where these numbers
c     are chosen such that ng1 = eps/2 * ng2, 2*ng1+ng2=thpts/th_nint.
c------------------------------------------------------------
      if(Geom.eq.'PRISMS')then
         pih = pi*0.5d0
         sinal  = dsin(alpha)
         pi2 = 2.d0*pi
         
         do i=1,phipts
            phi0 = phi(i)
            sinphi = dsin(phi0)
            cosphi = sqrt(1.d0-sinphi**2)
            theta0 = datan(eps * rho(alpha,sinphi,cosphi,sinal))

            xold=0.
            fac=theta0/dble(ng1)
            do iint=1,ng1
               do j=1,ngauss_per_int
                  k=(iint-1)*ngauss_per_int+j
                  wtheta(k,i)=fac * w0(j)
                  theta(k,i) =fac * x0(j) + xold
               enddo
               xold=xold+fac
            enddo
            fac=(pi-2*theta0)/dble(ng2)
            do iint=1,ng2
               do j=1,ngauss_per_int
                  k=(ng1+iint-1)*ngauss_per_int+j
                  wtheta(k,i)=fac * w0(j)
                  theta(k,i) =fac * x0(j) + xold
               end do
               xold=xold+fac
            enddo
            fac=theta0/dble(ng1)
            do iint=1,ng1
               do j=1,ngauss_per_int
                  k=(ng1+ng2+iint-1)*ngauss_per_int+j
                  wtheta(k,i)=fac*w0(j)
                  theta(k,i)=fac*x0(j) + xold
               enddo
               xold=xold+fac
            enddo
         end do
         ng1=ng1*ngauss_per_int
         ng2=ng2*ngauss_per_int
      else ! non-prismatic geometries
         fac=pi/dble(th_nint)
         do i=1,phipts
            xold=0.
            do iint=1,th_nint
               do j=1,ngauss_per_int
                  k=(iint-1)*ngauss_per_int+j
                  theta(k,i) =fac * x0(j) + xold 
                  wtheta(k,i)=fac * w0(j)
               enddo
               xold=xold+fac
            enddo
         enddo         
      endif

      return
      end



c***********************************************************************
c     This function returns the projection of r(theta,phi) in the
c     xy-plane divided by the 'width' a of the polygonial cylinder
c
c     Author: Michael Kahnert
c***********************************************************************    
      double precision function rho(alpha,sinphi,cosphi,sinal)
      implicit none
      real *8 alpha,sinphi,cosphi,sinal,diff
       
      diff=sinphi-sinal
      if (diff.lt.0.0) diff=-diff
      if (diff.lt.1e-12) then
         rho=0.5d0/cosphi
      else      
         rho = (cosphi-sinphi/dtan(alpha))/
     &    (1.d0-sinphi**2/sinal**2)
      end if
      return
      
      end
      


c***********************************************************************    
c     This function returns the partial derivative  of rho(theta,phi) 
c     wrt phi 
c
c     Author: Michael Kahnert
c***********************************************************************
      double precision function rhophi(alpha,sinphi,cosphi,sinal)
      implicit none
      real *8 alpha,sinphi,cosphi,sinal,diff,rho     
      
      diff=sinphi-sinal
      if (diff.lt.0.0) diff=-diff
      if (diff.lt.1e-12) then
         rhophi = sinphi/cosphi * rho(alpha,sinphi,cosphi,sinal)
      else      
         rhophi = (2.d0*sinphi*cosphi*rho(alpha,sinphi,cosphi,sinal)
     &    - sinal**2*(sinphi+cosphi/dtan(alpha)))/(sinal**2-sinphi**2)
      end if
      return
      
      end
      


c***********************************************************************
c     Surface parameterization r(theta,phi) and its
c     partial derivatives for a polyhedral cylinder, 
c     adapted to the special quadratrue scheme
c     for cylindrical particles:
c
c     Author: Michael Kahnert
C***********************************************************************
      subroutine R_poly_cyl(th1,th2,phi,alpha,ng1,ng2,
     &    ng1dim,ng2dim,phipts,N,pi,a,h,SIGMAh,
     &    r1,rth1,rphi1,r2,rth2,rphi2)
      implicit none
      integer ng1,ng2,ng1dim,ng2dim,phipts,
     &       N,ngj,i,j
      real*8 th1(ng1dim),th2(ng2dim,phipts),
     &       phi(phipts),r1(ng1dim),
     &       rth1(ng1dim),rphi1(ng1dim),
     &       r2(ng2dim,phipts),
     &       rth2(ng2dim,phipts),rphi2(ng2dim,phipts),pi,
     &       h,a,phi0,cosphi,sinphi,
     &       sitheta,cotheta,alpha,sinal,rho,rhophi
      logical SIGMAh

c-----------------------------------------------------------------------
c     Assign surface parameterization r(cos(theta),phi),
c     (d/d theta)r, (d/d phi)r:     
c-----------------------------------------------------------------------
      do j=1,ng1
	sitheta = dsin(th1(j))
	cotheta = dcos(th1(j))
	r1(j)   = h/cotheta
	rth1(j) = r1(j)*sitheta/cotheta
	rphi1(j)   = 0.d0
	if (.not.SIGMAh) then
	   ngj = ng1dim-j+1
	   r1(ngj)      = r1(j)
	   rth1(ngj)    =-rth1(j)
	   rphi1(ngj)   = rphi1(j)
	end if
      end do
c-----------------------------------------------------------------------
      do i=1,phipts
      	phi0  = phi(i)
	cosphi = dcos(phi0)
	sinphi = dsin(phi0)
	sinal = dsin(alpha)	
	do j=1,ng2/2
	   sitheta = dsin(th2(j,i))
	   cotheta = dcos(th2(j,i))
	   r2(j,i) = a/sitheta * rho(alpha,sinphi,cosphi,sinal)
	   rth2(j,i) =-r2(j,i)*cotheta/sitheta
	   rphi2(j,i)   = a/sitheta * rhophi(alpha,sinphi,cosphi,sinal)
	   if (.not.SIGMAh) then
	      ngj = ng2dim-j+1
	      r2(ngj,i)      = r2(j,i)
	      rth2(ngj,i)    =-rth2(j,i)
	      rphi2(ngj,i)   = rphi2(j,i)
	   end if
	end do
      end do
     
      return
      end
      
            
      
c***********************************************************************
c     Surface parameterization r(theta,phi) and its
c     partial derivatives for a polyhedral cylinder, 
c     adapted to the double quadratrue scheme:
c
c     Author: Michael Kahnert
C***********************************************************************
      subroutine R_poly (theta,phi,alpha,thpts,thpts2,ng1,ng2,
     &    phipts,N,pi,a,h,r,rth,rphi)
      implicit none
      integer thpts,thpts2,ng1,ng2,phipts,
     &       N,ngj,i,j
      real*8 theta(thpts,phipts),phi(phipts),r(thpts,phipts),
     &       rth(thpts,phipts),rphi(thpts,phipts),pi,
     &       h,a,phi0,cosphi,sinphi,
     &       sitheta,cotheta,alpha,sinal,rho,rhophi


c------------------------------------------------------------
c     Assign surface parameterization r(cos(theta),phi),
c     (d/d theta)r, (d/d phi)r:     
c------------------------------------------------------------
      do i=1,phipts
      	phi0  = phi(i)
	cosphi = dcos(phi0)
	sinphi = dsin(phi0)
	sinal = dsin(alpha)	
	do j=1,ng1
	   sitheta = dsin(theta(j,i))
	   cotheta = dcos(theta(j,i))
	   r(j,i)   = h/cotheta
	   rth(j,i) = r(j,i)*sitheta/cotheta
	   rphi(j,i)   = 0.d0
	   ngj = thpts-j+1
	   r(ngj,i)      = r(j,i)
	   rth(ngj,i)    =-rth(j,i)
	   rphi(ngj,i)   = rphi(j,i)
	end do
	do j=ng1+1,thpts2
	   sitheta = dsin(theta(j,i))
	   cotheta = dcos(theta(j,i))
	   r(j,i) = a/sitheta * rho(alpha,sinphi,cosphi,sinal)
	   rth(j,i) =-r(j,i)*cotheta/sitheta
	   rphi(j,i)   = a/sitheta * rhophi(alpha,sinphi,cosphi,sinal)
	   ngj = thpts-j+1
	   r(ngj,i)      = r(j,i)
	   rth(ngj,i)    =-rth(j,i)
	   rphi(ngj,i)   = rphi(j,i)
	end do
      end do
      
      return
      end
      
            
      
c***********************************************************************
c     Surface parameterization r(theta,phi) and its
c     partial derivatives for a cube, 
c     adapted to the double quadratrue scheme:
c
c     Author: Michael Kahnert
C***********************************************************************
      subroutine R_cube (theta,phi,alpha,thpts,thpts2,ng1,ng2,
     &    phipts,N,pi,a,h,r,rth,rphi)
      implicit none
      integer thpts,thpts2,ng1,ng2,phipts,phipts2,
     &       N,ngj,i,ii,j
      real*8 theta(thpts,phipts),phi(phipts),r(thpts,phipts),
     &       rth(thpts,phipts),rphi(thpts,phipts),pi,
     &       h,a,phi0,cosphi,sinphi,
     &       sitheta,cotheta,alpha,sinal,rho,rhophi

      if(N.ne.4)then
         write(6,*)'FATAL ERROR: set Nsym=4 for cube!!'
         stop 
      endif 
c------------------------------------------------------------
c     Assign surface parameterization r(cos(theta),phi),
c     (d/d theta)r, (d/d phi)r.
c     Note: for cubic groups the base-cube orientation
c     is rotated about the z-axis by 45 deg compared to 
c     polyhedrons with non-cubic point-group symmetries!
c------------------------------------------------------------
      phipts2=phipts/2
      do i=1,phipts
         if(i.le.phipts2)then
            ii=i+phipts2
         else
            ii=i-phipts2+1
         endif
      	phi0  = phi(i)
	cosphi = dcos(phi0)
	sinphi = dsin(phi0)
	sinal = dsin(alpha)	
	do j=1,ng1
	   sitheta = dsin(theta(j,i))
	   cotheta = dcos(theta(j,i))
	   r(j,ii)   = h/cotheta
	   rth(j,ii) = r(j,ii)*sitheta/cotheta
	   rphi(j,ii)   = 0.d0
	   ngj = thpts-j+1
	   r(ngj,ii)      = r(j,ii)
	   rth(ngj,ii)    =-rth(j,ii)
	   rphi(ngj,ii)   = rphi(j,ii)
	end do
	do j=ng1+1,thpts2
	   sitheta = dsin(theta(j,i))
	   cotheta = dcos(theta(j,i))
	   r(j,ii) = a/sitheta * rho(alpha,sinphi,cosphi,sinal)
	   rth(j,ii) =-r(j,ii)*cotheta/sitheta
	   rphi(j,ii)   = a/sitheta * rhophi(alpha,sinphi,cosphi,sinal)
	   ngj = thpts-j+1
	   r(ngj,ii)      = r(j,ii)
	   rth(ngj,ii)    =-rth(j,ii)
	   rphi(ngj,ii)   = rphi(j,ii)
	end do
      end do
      
      return
      end
      
            
c***********************************************************************
c     Surface parameterization r(theta,phi) and its
c     partial derivatives for a 2D Chebyshev particle, 
c     using the double quadratrue scheme:
c
c     Author: Michael Kahnert
C***********************************************************************
      subroutine R_Cheb2D (theta,phi,thpts,thpts2,
     &    phipts,a,h,r,rth,rphi,cheborder)
      implicit none
      integer thpts,thpts2,phipts,cheborder,
     &       ngj,i,ii,j
      real*8 theta(thpts,phipts),phi(phipts),r(thpts,phipts),
     &       rth(thpts,phipts),rphi(thpts,phipts),
     &       h,a

c------------------------------------------------------------
c     Assign surface parameterization r(cos(theta),phi),
c     (d/d theta)r, (d/d phi)r.
c     In this subroutine, a=radius of base sphere, h=perturbation
c     amplitude
c------------------------------------------------------------
      do i=1,phipts
	do j=1,thpts
	   r(j,i)   =  a*(1d0+h*dcos(cheborder*theta(j,i)))
	   rth(j,i) = -a*cheborder*h*dsin(cheborder*theta(j,i))
	   rphi(j,i)= 0d0
	end do
      end do
      
      return
      end
      
            
c***********************************************************************
c     Surface parameterization r(theta,phi) and its
c     partial derivatives for a 3D Chebyshev particle, 
c     using the double quadratrue scheme:
c
c     Author: Michael Kahnert
C***********************************************************************
      subroutine R_Cheb3D (theta,phi,thpts,thpts2,
     &    phipts,a,h,r,rth,rphi,cheborder)
      implicit none
      integer thpts,thpts2,phipts,cheborder,
     &       ngj,i,ii,j
      real*8 theta(thpts,phipts),phi(phipts),r(thpts,phipts),
     &       rth(thpts,phipts),rphi(thpts,phipts),
     &       h,a

c------------------------------------------------------------
c     Assign surface parameterization r(cos(theta),phi),
c     (d/d theta)r, (d/d phi)r.
c     In this subroutine, a=radius of base sphere, h=perturbation
c     amplitude
c------------------------------------------------------------
      do i=1,phipts
	do j=1,thpts
	   r(j,i)   =  a*(1d0+h*dcos(cheborder*theta(j,i)) * 
     *          dcos(cheborder*phi(i)))
	   rth(j,i) = -a*cheborder*h*dsin(cheborder*theta(j,i)) * 
     *          dcos(cheborder*phi(i))
	   rphi(j,i)= -a*cheborder*h*dcos(cheborder*theta(j,i)) * 
     *          dsin(cheborder*phi(i))
	end do
      end do
      
      return
      end
      
            
c***********************************************************************
c     Surface parameterization r(theta,phi) and its
c     partial derivatives for a homogeneous sphere, 
c     using the double quadratrue scheme:
c
c     Author: Michael Kahnert
C***********************************************************************
      subroutine R_Sphere (theta,phi,thpts,thpts2,
     &    phipts,a,r,rth,rphi)
      implicit none
      integer thpts,thpts2,phipts,phipts2,
     &       ngj,i,ii,j
      real*8 theta(thpts,phipts),phi(phipts),r(thpts,phipts),
     &       rth(thpts,phipts),rphi(thpts,phipts),
     &       h,a

c------------------------------------------------------------
c     Assign surface parameterization r(cos(theta),phi),
c     (d/d theta)r, (d/d phi)r.
c     In this subroutine, a=radius of base sphere
c------------------------------------------------------------
      phipts2=phipts/2
      do i=1,phipts
	do j=1,thpts
	   r(j,i)   =  a
	   rth(j,i) = 0d0
	   rphi(j,i)= 0d0
	end do
      end do
      
      return
      end
      
            
c***********************************************************************
c     Surface parameterization r(theta,phi) and its
c     partial derivatives for a Gaussian random sphere, 
c     using the double quadratrue scheme:
c
c     Author: Michael Kahnert
C***********************************************************************
      subroutine R_GRS (theta,phi,thpts,
     &    phipts,a,amplsig,anglcorr,r,rth,rphi,lreadgeom,ngrs)
      implicit none
      integer thpts,phipts,phipts2,i,j,lmin,lmax,ngrs
      real*8 theta(thpts,phipts),phi(phipts),r(thpts,phipts),
     &       rth(thpts,phipts),rphi(thpts,phipts),
     &       a,amplsig,anglcorr
      logical lreadgeom,lcompute
      integer thptsc,phiptsc


      open(15,file='sfcpar.dat')

      if(lreadgeom)then
         lcompute=.false.
      else
         lcompute=.true.
      endif

      if(lreadgeom)then
         read(15,*,err=1000,end=1000)thptsc,phiptsc 
         if(thptsc.ne.thpts.or.phiptsc.ne.phipts)goto 1000
         read(15,*,end=1000,err=1000)r
         read(15,*,end=1000,err=1000)rth
         read(15,*,end=1000,err=1000)rphi
         goto 2000
 1000    write(6,*)
     &  '      ERROR reading GRS file; computing geometry instead...'
         lcompute=.true.
      endif
 
      if(lcompute)then
         write(6,*)'      Computing GRS geometry...   '
         lmin=0
         lmax=int(300/anglcorr)+1
         
         call gsphere(1,2,amplsig,1d0,anglcorr,lmin,lmax,thpts,phipts,
     &        theta,phi,ngrs,r,rth,rphi)
         do i=1,thpts
            do j=1,phipts
               r(i,j)=a*r(i,j)
               rth(i,j)=a*rth(i,j)
               rphi(i,j)=a*rphi(i,j)
            enddo
         enddo
         
         write(15,*)thpts,phipts
         write(15,*)r
         write(15,*)rth
         write(15,*)rphi
      endif

 2000 close(15)

      end
      


c***********************************************************************
c     Surface parameterization r(theta,phi) and its
c     partial derivatives for a 2D-Gaussian random sphere:
c
c     Author: Michael Kahnert
C***********************************************************************
      subroutine R_GRS2D (theta,phi,thpts,
     &    phipts,a,amplsig,anglcorr,r,rth,rphi,lreadgeom,ngrs)
      implicit none
      integer thpts,phipts,phipts2,i,j,lmin,lmax,ngrs
      real*8 theta(thpts,phipts),phi(phipts),r(thpts,phipts),
     &       rth(thpts,phipts),rphi(thpts,phipts),
     &       a,amplsig,anglcorr
      logical lreadgeom,lcompute
      integer thptsc,phiptsc


      open(15,file='sfcpar.dat')

      if(lreadgeom)then
         lcompute=.false.
      else
         lcompute=.true.
      endif

      if(lreadgeom)then
         read(15,*,err=1000,end=1000)thptsc,phiptsc 
         if(thptsc.ne.thpts.or.phiptsc.ne.phipts)goto 1000
         read(15,*,end=1000,err=1000)r
         read(15,*,end=1000,err=1000)rth
         read(15,*,end=1000,err=1000)rphi
         goto 2000
 1000    write(*,*)
     &  '      ERROR reading GRS file; computing geometry instead...'
         lcompute=.true.
      endif
 
      if(lcompute)then
         write(*,*)'      Computing 2D-GRS geometry...   '
         lmin=0
         lmax=int(300/anglcorr)+1
         
         call gsphere(2,2,amplsig,1d0,anglcorr,lmin,lmax,thpts,phipts,
     &        theta,phi,ngrs,r,rth,rphi)
         do i=1,thpts
            r(i,1)=a*r(i,1)
            rth(i,1)=a*rth(i,1)
            rphi(i,1)=0d0
            do j=2,phipts
               r(i,j)=r(i,1)
               rth(i,j)=rth(i,1)
               rphi(i,j)=0d0
            enddo
         enddo

         write(15,*)thpts,phipts
         write(15,*)r
         write(15,*)rth
         write(15,*)rphi
      endif

 2000 close(15)

      end


c***********************************************************************
c     Surface parameterization r(theta,phi) and its
c     partial derivatives for a Gaussian random sphere/Chebyshev hybrid, 
c     using the double quadratrue scheme.
c     In this hybrid particle, the surface perturbation in the 
c     polar direction is identical to a 2D Gaussian random sphere,
c     while the azimuthal perturbation is a Chebyshev polynomial.
c
c     Author: Michael Kahnert
C***********************************************************************
      subroutine R_GRSCheb (theta,phi,thpts,phipts,a,amplsig,
     &    anglcorr,r,rth,rphi,lreadgeom,ngrs,cheborder)
      implicit none
      integer thpts,phipts,phipts2,i,j,lmin,lmax,ngrs,cheborder
      real*8 theta(thpts,phipts),phi(phipts),r(thpts,phipts),
     &       rth(thpts,phipts),rphi(thpts,phipts),
     &       a,amplsig,anglcorr,h,dl,dr(thpts),rth1(thpts)
      logical lreadgeom,lcompute
      integer thptsc,phiptsc

      open(15,file='sfcpar.dat')

      if(lreadgeom)then
         lcompute=.false.
      else
         lcompute=.true.
      endif

      if(lreadgeom)then
         read(15,*,err=1000,end=1000)thptsc,phiptsc 
         if(thptsc.ne.thpts.or.phiptsc.ne.phipts)goto 1000
         read(15,*,end=1000,err=1000)r
         read(15,*,end=1000,err=1000)rth
         read(15,*,end=1000,err=1000)rphi
         goto 2000
 1000    write(*,*)
     &  '      ERROR reading GRS file; computing geometry instead...'
         lcompute=.true.
      endif

      if(lcompute)then
         write(*,*)'      Computing GRS/Cheb geometry...   '
         lmin=0
         lmax=int(300/anglcorr)+1
         h=amplsig
         
         call gsphere(2,2,amplsig,1d0,anglcorr,lmin,lmax,thpts,phipts,
     &        theta,phi,ngrs,r,rth,rphi)

         dl=dble(cheborder)

         do i=1,thpts
            dr(i)=r(i,1)-1d0
            rth1(i)=a*rth(i,1)
         enddo

         do i=1,thpts
            do j=1,phipts
               r(i,j)   = a* (1d0 + dr(i) *dcos(dl*phi(j)))
               rth(i,j) =    rth1(i)      *dcos(dl*phi(j))
               rphi(i,j)=-a*dl     *dr(i) *dsin(dl*phi(j))
            enddo
         enddo

        


         write(15,*)thpts,phipts
         write(15,*)r
         write(15,*)rth
         write(15,*)rphi
      endif

 2000 close(15)

      

      end


C***********************************************************************
c     This subroutine calculates the critical azimuths phic
c     for each th1 < theta < th2
c
c     Author: Michael Kahnert
C***********************************************************************
      subroutine phi_crit(eps,cosal,sinal,theta,pi,ngstart,ng1,phic)
      implicit none
      integer ngstart,ng1,i,j
      real*8 eps,pi,cosal,sinal,theta(ng1),phic(ng1),
     &       A,B,C,D,tanth,r1,prec,imr,rro,rr
      complex coeff(5),roots(4)
     
      prec = 1.0d-6
      A = -2.d0*eps*cosal*sinal
      B = sinal**2
      C = -A*B
      D = B**2
      do i=ngstart,ng1
         tanth=dtan(theta(i))
         r1=(eps/tanth)**2
         coeff(1)=cmplx(D*(1.d0-r1),0.d0)
         coeff(2)=cmplx(C/tanth,0.d0)
         coeff(3)=cmplx(B*(r1-2.d0),0.d0)
         coeff(4)=cmplx(A/tanth,0.d0)
         coeff(5)=cmplx(1.d0,0.d0)
         call zroots(coeff,4,roots,.true.)
c-----------------------------------------------------------------------
c        take the smallest real root:
c-----------------------------------------------------------------------
         rro=pi
         do j=1,4
            imr=imag(roots(j))
            if (imr.lt.0.0) imr=-imr
            if (imr.lt.prec) then
               rr=real(roots(j))
               if (rr.lt.rro.and.rr.ge.0.0) rro=rr
            end if
         end do
         if (rro.eq.pi) then
            write(6,*)'>>>ERROR: crit. angle phic could not be found'
            stop
         end if
         phic(i) = asin(rro)
 
c         write(6,*)'     i,th,phic:',i,theta(i)*180.0/pi,phic(i)*180.0/pi
      end do
      
      return 
      end
      

C***********************************************************************
c     Write geometry in a format that can be 3D-plotted by Matlab.
c
c     Author: Michael Kahnert
C***********************************************************************
      subroutine write_matlab(theta,phi,thpts,phipts,r,Nsym)
c     Write data file for plotting particle geometry with Matlab
      implicit none
      integer thpts,phipts,i,j,k,kk,Nsym
      real*8 theta(thpts,phipts),phi(phipts),r(thpts,phipts)
      real*8 mu,nu,phip,pi
      real*8, dimension(:,:,:), allocatable :: xs
      integer phidim,err
      
      pi=acos(-1d0)

      phidim=phipts*Nsym
      allocate(xs(thpts,phidim,3),stat=err)
      if(err.ne.0)then
         write(*,*)'FATAL ERROR: allocation failed for xs'
         stop
      endif

      
      do i=1,thpts
         do k=1,Nsym
            do j=1,phipts
               kk=j+(k-1)*phipts
               mu=cos(theta(i,j))
               nu=sqrt(1d0-mu**2)
               phip=(k-1)*2d0*pi/dble(Nsym)+phi(j)
               xs(i,kk,1)=r(i,j)*nu*cos(phip)
               xs(i,kk,2)=r(i,j)*nu*sin(phip)
               xs(i,kk,3)=r(i,j)*mu           
            enddo
         enddo
      enddo

      open(unit=1, file='matlabx.out') ! Matlab
      open(unit=2, file='matlaby.out')
      open(unit=3, file='matlabz.out')
      do i=1,thpts
         write (1,1000) (XS(i,j,1),j=1,phidim)
         write (2,1000) (XS(i,j,2),j=1,phidim)
         write (3,1000) (XS(i,j,3),j=1,phidim)
      enddo
      close(unit=3)
      close(unit=2)
      close(unit=1)

      deallocate(xs)

 1000 format(9999e18.9)

      end



c=======================================================================
c     PART V :   M I S C E L L A N E O U S   R O U T I N E S :
c=======================================================================

c***********************************************************************
c     set some parameters and print input information:
c
c     Author: Michael Kahnert
c***********************************************************************
      subroutine setup(Geom,cheborder,prtb_order,lprtb,lirredrep,Nsym,
     &     chtabfile,igeom,gpar1,gpar2,gpar3,ip1,ip2,mr,k,lambda,
     &     nmax,mmax,nmmax,thpts,phipts,thpts2,ngauss_per_int,
     &     th_nint,phi_nint,thptsx,phiptsx,
     &     Pgroup,SIGMAh,DIHEDRAL,AXISYM,
     &     Cyl_Quad,pi,a,norient,h2,eps,ng1,ng2,ng1dim,ng2dim,neuler,
     &     alpha,beta,gamma,nsi,lsi,lpsi,nmx,mmx,mmn,ngrs,Rv)     
     
      implicit none
      character*6 Geom
      character*(*) chtabfile,Pgroup
      integer cheborder,prtb_order,norient
      integer Nsym,igeom,thpts,phipts,ng1,ng2,ng1dim,ng2dim,
     &     thpts2,thptsx,phiptsx,nmax,mmax,nmmax,ip1,ip2,
     &     neuler,ngauss_per_int,th_nint,phi_nint,ngrs
      integer nsi(nmmax),lsi(-mmax:mmax,nmmax),lpsi(0:mmax,nmmax),
     &     nmx(-mmax:mmax),mmx(nmmax),mmn(nmmax),nn,mm,jj,ieuler
      real*8  gpar1,gpar2,gpar3,h,l,a,h2,eps,pi,rv,ra,dN,cf,
     &     alpha(neuler),beta(neuler),gamma(neuler),k,lambda,xv
      complex*16 mr
      logical SIGMAh,DIHEDRAL,AXISYM,Cyl_Quad,lprtb,lirredrep
      
      dN = dble(Nsym)
c-----------------------------------------------------------------------
c     geometrical parameters of particle
c-----------------------------------------------------------------------
      if(Geom.eq.'PRISMS')then
         if (igeom.eq.1) then
            l   = gpar1
            h   = gpar2
            a   = l/(2.d0*dsin(pi/dN))
            h2  = h/2.d0
            eps = a/h2
            rv  = (3.d0/(16.d0*pi)*dN*l*l*h/dtan(pi/dN))**(1.d0/3.d0)
            ra  = dsqrt(dN*l/(4.d0*pi)*(l/2.d0/dtan(pi/dN)+h))
         else if (igeom.eq.2) then
            a   = gpar1
            h2  = gpar2
            h   = 2.d0*h2
            eps = a/h2
            l   = 2.d0*a*dsin(pi/dN)
            rv  = (3.d0/(16.d0*pi)*dN*l*l*h/dtan(pi/dN))**(1.d0/3.d0)
            ra  = dsqrt(dN*l/(4.d0*pi)*(l/2.d0/dtan(pi/dN)+h))
         else if (igeom.eq.3) then
            rv  = gpar1
            eps = gpar2
            a   = ((2.d0*pi*eps/
     &           (3.d0*dN*dsin(pi/dN)*dcos(pi/dN)))**(1.d0/3.d0))*rv
            h   = 2.d0*a/eps
            h2  = h/2.d0
            l   = 2.d0*a*dsin(pi/dN)
            ra  = dsqrt(dN*l/(4.d0*pi)*(l/2.d0/dtan(pi/dN)+h))      
         else if (igeom.eq.4) then
            ra  = gpar1
            eps = gpar2
            a   = dsqrt(2.d0*pi/(dN*dsin(pi/dN)*
     &           (dcos(pi/dN)+2.d0/eps))) * ra
            h   = 2.d0*a/eps
            h2  = h/2.d0
            l   = 2.d0*a*dsin(pi/dN)
            rv  = (3.d0/(16.d0*pi)*dN*l*l*h/dtan(pi/dN))**(1.d0/3.d0)
         end if
      elseif(Geom.eq.'CHEB2D'.or.Geom.eq.'CHEB3D'.or.Geom.eq.'SPHERE'
     &        .or.Geom.eq.'GRSPHR'.or.Geom.eq.'GRSP2D'
     &        .or.Geom.eq.'GRSCHB')
     &        then
         a=gpar1
         h2=gpar2
         eps=gpar3
         Rv=a
      endif
      xv=k*Rv
      if((Geom.eq.'CHEB3D'.or.Geom.eq.'GRSCHB').and.cheborder.ne.Nsym)
     &then
         write(6,*)'WARNING: Nsym not equal to cheborder'
      endif
      if(Geom.eq.'GRSPHR'.and.Nsym.ne.1)then
         write(6,*)'FATAL ERROR: Nsym must be 1 for GRSPHR'
         stop
      endif
      if((Geom.eq.'CHEB3D'.or.Geom.eq.'CHEB2D'.or.Geom.eq.'GRSPHR'
     &     .or.Geom.eq.'GRSP2D'.or.Geom.eq.'GRSCHB')
     &     .and.(h2.lt.0d0.or.h2.gt.1))then
         write(6,*)'FATAL ERROR: gpar2 must be between 0 and 1'
         stop
      endif
      if((Geom.eq.'GRSPHR'.or.Geom.eq.'GRSP2D')
     &     .and.(gpar3.lt.0d0.or.gpar3.gt.360d0))then
         write(6,*)'FATAL ERROR: gpar3 must be between 0 and 360'
         stop
      endif
      if(Geom.eq.'CHEB2D'.or.Geom.eq.'SPHERE'.or.Geom.eq.'GRSP2D') then
         AXISYM=.true.
      else
         AXISYM=.false.
      endif
c-----------------------------------------------------------------------
c     Symmetries
c-----------------------------------------------------------------------
      if(Pgroup.eq.'Dn'.or.Pgroup.eq.'Dnh'.or.Pgroup.eq.'Dnd'.or.
     &     Pgroup.eq.'T'.or.Pgroup.eq.'Td'.or.Pgroup.eq.'Th'.or.
     &     Pgroup.eq.'O'.or.Pgroup.eq.'Oh'.or.Pgroup.eq.'I'.or.
     &     Pgroup.eq.'Ih') then
         DIHEDRAL=.true.
      else
         DIHEDRAL=.false.
      endif
      if(Pgroup.eq.'Cs'.or.Pgroup.eq.'Cnh'.or.Pgroup.eq.'Dnh'.or.
     &     Pgroup.eq.'Th'.or.Pgroup.eq.'Oh'.or.Pgroup.eq.'Ih') then
         SIGMAh=.true.
      else
         SIGMAh=.false.
      endif
c-----------------------------------------------------------------------
c     Integration scheme
c-----------------------------------------------------------------------
c     If you want to enable Cyl_Quad=.true., then you have to comment
c     out the following lines:
c>>>>>>>>>>>>>>>>>>>
      if(Cyl_Quad)then
         write(6,*)
     &'NOTE: Use of Cyl_Quad=.true. requires the subroutines zroots'
         write(6,*)
     &'and laguer from Numerical Recipes => currently not implemented!'
         write(6,*)'To use Cyl_Quad=.true.,'
         write(6,*)' - Get a License for Numerical Recipes'
         write(6,*)
     &        ' - Include the subroutines zroots and laguer in NumRec.F'
         write(6,*)' - Edit setup.F and disable the automatic setting of
     & Cyl_Quad=.false.'
         write(6,*)
     & '(This is found under the comment heading "Integration scheme")'
         Cyl_Quad=.false.
      endif
c<<<<<<<<<<<<<<<<<<<
      if(Geom.ne.'PRISMS'.and.Cyl_Quad)then
         write(6,*)
     &        'Cyl_Quad must be false for geometries other than PRISMS'
         Cyl_Quad=.false.
      endif
c-----------------------------------------------------------------------
c     Q-matrix inversion method
c-----------------------------------------------------------------------
      if(lirredrep.and.lprtb)then
         write(6,*)
     &        'FATAL ERROR: lirredrep and lprtb must not both be true'
         stop
      endif
      if(lprtb.and.(Geom.ne.'CHEB2D'.and.Geom.ne.'CHEB3D'
     &     .and.Geom.ne.'GRSPHR'.and.Geom.ne.'GRSP2D'
     &     .and.Geom.ne.'GRSCHB'))then
         write(6,*)
     &   'FATAL ERROR: Perturbation approach not implemented for ',Geom
         stop
      endif
c-----------------------------------------------------------------------
c     The surface integrals are evaluated by Gauss-Kronrad quadrature.
c     We use th_nint intervals to split up the polar angular range, 
c     and phi_nint intervals for the azimuthal angular range. 
c     In each interval we compute ngauss_per_int
c     quadrature points and weights (currently ngauss_per_int=16).
c     So in total there are thpts=ngauss_per_int*th_nint discrete 
c     quadrature points and weights for the polar integration, and
c     phipts=ngauss_per_int*phi_nint for the azimuthal integration.
c     For axisymmetric particles, we set phipts=2, phiptsx=1.
c-----------------------------------------------------------------------
      thpts=ngauss_per_int*th_nint
      phipts=ngauss_per_int*phi_nint
      thpts2=thpts/2
      if(AXISYM) phipts=2
c-----------------------------------------------------------------------
c     The number of polar quadrature integrals in the
c     subintervals [0,theta0] is ng1/ngauss_per_int. So the
c     number of quadrature points in that interval is ng1.
c     The number of quadrature points in the intervals [theta0,pi-theta0] and
c     [pi-theta0,pi] are ng2 and ng1, resp., where these numbers
c     are chosen such that ng1/ngauss_per_int = eps/2 * ng2/ngauss_per_int, 
c     2*ng1+ng2=thpts.
c-----------------------------------------------------------------------
      if(Geom.eq.'PRISMS')then
         ng1 = int(dble(th_nint)*eps/(2.d0*(1.d0+eps)))
         ng1=ng1*ngauss_per_int
         ng2 = thpts-2*ng1
      else
         ng1=0
         ng2=thpts
      endif
      if (Cyl_Quad) then
         if (SIGMAh) then
            ng1dim=ng1
            ng2dim=ng2/2
         else
            ng1dim=2*ng1
            ng2dim=ng2
         end if
      else
         ng1dim=1
         ng2dim=thpts
         if (SIGMAh) then
            thptsx=thpts2
         else
            thptsx=thpts
         end if
      end if
c-----------------------------------------------------------------------
c     In a similar way, the range of the phi-integration is
c     reduced for dihedral symmetry:
c-----------------------------------------------------------------------
      if (DIHEDRAL) then
         if (mod(phipts,2).ne.0) then
          write(6,*)'>>>ERROR: phipts must be even when DIHEDRAL=.true.'
          write(6,*)'   Change setting for phipts in file params !'
          stop            
         end if
         phiptsx=phipts/2
      else
         phiptsx=phipts
      end if
c-----------------------------------------------------------------------
c     setting of some super-indices:
c-----------------------------------------------------------------------
c     index convention for expansion of the fields: 
c     \sum_{m=-M}^M \sum_{n=|m|}^{N}...
c--------------------
      do nn=1,nmmax
         nsi(nn)=nn*(nn+1)
c--------------------------------
c     min and max |m| for each n:
c--------------------------------
         mmn(nn)=0
         if (nn.le.mmax) then
            mmx(nn)=nn
         else
            mmx(nn)=mmax
         end if
      end do
      do mm=-mmax,mmax
         do nn=max(1,abs(mm)),nmmax
            lsi(mm,nn)=nn-max(1,abs(mm))+1+(mm+mmax)*(nmax+1)
     &           -(mmax*(mmax+1))/2-(mm*(abs(mm-1)))/2
            if(mm.gt.0) lsi(mm,nn)=lsi(mm,nn)-1
         end do
      end do
      do mm=0,mmax
         do nn=max(1,abs(mm)),nmmax
            lpsi(mm,nn)=mm*(nmax+1)-(mm*(mm-1))/2+nn-max(1,abs(mm))
            if(mm.eq.0) lpsi(mm,nn)=lpsi(mm,nn)+1
         end do
      end do
      do mm=-mmax,mmax
         nmx(mm)=nmax
      end do
c-----------------------------------------------------------------------
c     print input information:
c-----------------------------------------------------------------------
      write(6,*)
      write(6,1037)
      write(6,*)
      if(Geom.eq.'PRISMS')then
         write(6,1001) Nsym
      elseif(Geom.eq.'SPHERE')then
         write(6,1061)
      elseif(Geom.eq.'GRSPHR')then
         write(6,1062)ngrs 
      elseif(Geom.eq.'GRSP2D')then
         write(6,1064)ngrs
      elseif(Geom.eq.'CHEB2D')then
         write(6,1065) cheborder
      elseif(Geom.eq.'CHEB3D')then
         write(6,1067) cheborder
      elseif(Geom.eq.'GRSCHB')then
         write(6,1068) cheborder
      else
         write(6,1069)
         stop
      endif 
      write(6,*)
      write(6,1037)
      write(6,1033) lambda
      write(6,1034) k
      write(6,1037)
      write(6,*) 'Particle Properties:'
      if(Geom.eq.'PRISMS')then
         write(6,1003) l
         write(6,1004) h
         write(6,1041) a
         write(6,1005) rv
         write(6,1050) ra
         write(6,1009) k*rv
         write(6,1006) eps
      elseif(Geom.eq.'SPHERE')then
         write(6,1042) a
      elseif(Geom.eq.'CHEB2D'.or.Geom.eq.'CHEB3D')then
         write(6,1043) a
         write(6,1044) cheborder
         write(6,1045) h2
         write(6,1009) k*a
      elseif(Geom.eq.'GRSPHR'.or.Geom.eq.'GRSP2D')then
         write(6,1043) a
         write(6,1046) h2
         write(6,1047) eps
         write(6,1009) k*a
      elseif(Geom.eq.'GRSCHB')then
         write(6,1043) a
         write(6,1045) h2
         write(6,1046) h2
         write(6,1044) cheborder
         write(6,1047) eps
         write(6,1009) k*a
      endif

      write(6,1051) dreal(mr),imag(mr)
      
      write(6,1010) Pgroup,Nsym

      if (norient.eq.1) then
         write(6,1066)
      elseif(norient.eq.2) then
         write(6,1063)
      elseif(norient.eq.3)then
         write(6,1002)
      endif
      if (norient.ge.2) then
         write(6,1037)
         write(6,*) 'Particle orientations (Euler angles):'
         do ieuler=1,neuler
            write(6,1030) alpha(ieuler),beta(ieuler),gamma(ieuler)
            if (beta(ieuler).lt.0.0.or.beta(ieuler).gt.180.0) then
               write(6,*)'>>>ERROR: beta must be between 0.0 and 180.0'
               write(6,*)'   Change setting for beta in file params !'
               stop  
            end if 
c-----------------------------------------------------------------------
c     CONVERSION OF THE EULER ANGLES:
c
c     This programme uses the convention of PASSIVE rotation:
c     - First the coordinate system {x,y,z} is rotated about an angle 
c       alpha' about the z-axis => new coordinate system {x',y',z'}, 
c       z=z'
c     - Then the coordinate system {x',y',z'} is rotated about and 
c       angle beta' about the new y'-axis =>  new coordinate system 
c       {x'',y'',z''}, y'=y''
c     - finally the coordinate system {x'',y'',z''} is rotated about an 
c       angle gamma' about the z''-axis =>  new coordinate system 
c       {x''',y''',z'''}, z''=z'''
c
c     The active rotation of the particle with the Euler angles
c     (alpha,beta,gamma) as described in the comments in the parameter
c     input file corresponds to these passive Euler angles according to:
c     
c     alpha' = 180 deg - alpha
c     beta'  = beta
c     gamma' = 180 deg - gamma
c-----------------------------------------------------------------------
            alpha(ieuler) = 180.d0-alpha(ieuler)
            beta(ieuler)  = beta(ieuler)
            gamma(ieuler) = 180.d0-gamma(ieuler)
c-----------------------------------------------------------------------
c     conversion into radians
c-----------------------------------------------------------------------
            cf=pi/180d0
            alpha(ieuler)=cf*alpha(ieuler)
            beta(ieuler)=cf*beta(ieuler)
            gamma(ieuler)=cf*gamma(ieuler)
         enddo
      end if
      write(6,1037)
      write(6,*) 'Parameters concerning the calculation:'
      write(6,1007) nmax,mmax,ng1,ng2
      write(6,1008) thpts,phipts
      if (Cyl_Quad) then
         write(6,1036)
      else
         write(6,1035)
      end if
      if (lirredrep)then
         write(6,1039)
      endif
      
      if(lirredrep)then
         open(14,file=chtabfile,err=100)
         goto 110
 100     write(6,*)
     &  'FATAL ERROR: character table ',chtabfile,' not found'
         stop
 110     close(14)
      endif

      if(Geom.eq.'CHEB2D'.or.Geom.eq.'CHEB3D'.or.Geom.eq.'GRSPHR'
     &     .or.Geom.eq.'GRSP2D'.or.Geom.eq.'GRSCHB')then
         if(.not.lprtb)then
            write(6,1070)
         elseif(prtb_order.eq.1)then
            write(6,1071) prtb_order
         elseif(prtb_order.eq.2)then
            write(6,1072) prtb_order
         elseif(prtb_order.eq.3)then
            write(6,1073) prtb_order
         else
            write(6,1074) prtb_order
         endif
      endif
      write(6,1038)
      write(6,*)
      
c-----------------------------------------------------------------------
 1001 format(' SCATTERING BY ',i4,'-HEDRAL PRISMS')     
 1061 format(' SCATTERING BY HOMOGENEOUS SPHERES')     
 1062 format(' SCATTERING BY GAUSSIAN RANDOM SPHERES, SEED=',i4)     
 1064 format(' SCATTERING BY 2D-GAUSSIAN RANDOM SPHERES, SEED=',i4)     
 1065 format(' SCATTERING BY 2D-CHEBYSHEV PARTICLES OF ORDER ',i4)     
 1067 format(' SCATTERING BY 3D-CHEBYSHEV PARTICLES OF ORDER ',i4)     
 1068 format(' SCATTERING BY GRS-CHEBYSHEV HYBRIDS OF ORDER ',i4)     
 1069 format('FATAL ERROR: Invalid geometry')
 1002 format(' Computations for both random and fixed orientations')
 1066 format(' Computations for randomly orientated particles')
 1063 format(' Computations for particles in fixed orientation')
 1003 format(' polygon side length             l  : ',f15.9)
 1004 format(' cylinder height                 h  : ',f15.9) 
 1041 format(' polygon radius                  a  : ',f15.9) 
 1042 format(' sphere radius                   r  : ',f15.9) 
 1043 format(' radius of unperturbed sphere    r  : ',f15.9) 
 1044 format(' Chebyshev order                 n  : ',i5) 
 1045 format(' deformation parameter           e  : ',f15.9) 
 1046 format(' relative radial standard deviation : ',f15.9)
 1047 format(' angular correlation (degrees)      : ',f15.9)
 1005 format(' volume-equivalent sphere radius Rv : ',f15.9)
 1050 format(' area - equivalent sphere radius Ra : ',f15.9)
 1009 format(' size parameter                k*Rv : ',f15.9)
 1006 format(' aspect ratio                    eps: ',f15.9)
 1007 format(' nmax:',i5,'    mmax:',i7,'     ng1:',i7,'     ng2:',i4)
 1008 format(' thpts:',i4,'     phipts:',i4)


 1010 format(' particle symmetry: ',a3,12x,'n:',i6)

 1030 format(' (alpha,beta,gamma) in degrees: ',3(f6.1,1x,','))
 1033 format(' wavelength                   lambda: ',f15.9)
 1034 format(' wavenumber                   k     : ',f15.9)
 1035 format(' integration method: double Gauss-Legendre quadrature')
 1036 format(
     &' integration method: special quadrature scheme for cylinders')
 1039  format
     &  (' Q-matrix inverted by exploiting irreducible representations')
 1070 format (' Scattering problem solved non-perturbatively')
 1071 format (' Scattering problem solved to ',i3,'st order perturbation
     & theory')
 1072 format (' Scattering problem solved to ',i3,'nd order perturbation
     & theory')
 1073 format (' Scattering problem solved to ',i3,'rd order perturbation
     & theory')
 1074 format (' Scattering problem solved to ',i3,'th order perturbation
     & theory')
 1037 format(
     &'=================================================================
     &====')
 1038 format(
     &'-----------------------------------------------------------------
     &----')
 1051 format
     & (' refractive index                mr : ',f15.9,' + i*',f12.9)
      return 
      end

C***********************************************************************
c     calculation of assciated Legendre functions from 
c     Wigner d-functions:
c
C***********************************************************************
      subroutine wigner(phipts,phiptsx,nmmax,mmax,lpsi,lpmax,
     &          ng1,ng2,ng1dim,ng2dim,thptsx,xth1,xth2,P1,dP1,P2,dP2)
     
      implicit none
      
      integer i,j,n,m,l,nmx
      integer phipts,phiptsx,nmmax,mmax,lpmax,
     & ng1,ng2,ng1dim,ng2dim,thptsx,lpsi(0:mmax,nmmax)
      logical Cyl_quad
      real*8 xth1(ng1dim), xth2(ng2dim,phipts),
     & Ptmp(nmmax),dPtmp(nmmax),cth,
     & P1(ng1dim,lpmax),dP1(ng1dim,lpmax),
     & P2(thptsx,phiptsx,lpmax),dP2(thptsx,phiptsx,lpmax)

      if (ng2dim.lt.thptsx) then
         write(6,*)'>>>ERROR: ng2dim.lt.thptsx in subroutine wigner'
         stop
      end if
      do  i = 1, ng1dim
         cth = dcos(xth1(i))
c-----------------------------------------------------------------------
c        m=0:
c-----------------------------------------------------------------------
         call vig( cth, nmmax, 0, Ptmp, dPtmp)
         do n = 1, nmmax
            l = lpsi(0,n)
            P1(i,l) = Ptmp(n)
            dP1(i,l) = dPtmp(n)
         end do
         do  m = 1, mmax
c-----------------------------------------------------------------------
c           m>0 (m<0 is given by symmetry relations and therefore not
c                computed):
c-----------------------------------------------------------------------
            call vig( cth, nmmax, m, Ptmp, dPtmp)
            do  n = m, nmmax
               l = lpsi(m,n)
               P1(i,l) = Ptmp(n)
               dP1(i,l) = dPtmp(n)
	    end do
	 end do
      end do
     
      do  i = 1, thptsx
         do j =1, phiptsx
            cth = dcos(xth2(i,j))
            call vig( cth, nmmax, 0, Ptmp, dPtmp)
            do n = 1, nmmax
               l = lpsi(0,n)
               P2(i,j,l) = Ptmp(n)
               dP2(i,j,l) = dPtmp(n)
            end do
            do  m = 1, mmax
               call vig( cth, nmmax, m, Ptmp, dPtmp)
               do  n = m, nmmax
                  l = lpsi(m,n)
                  P2(i,j,l) = Ptmp(n)
                  dP2(i,j,l) = dPtmp(n)
	       end do
	    end do
         end do
      end do

      
 1000 format(4i4,f15.10)
      return
      end



C***********************************************************************
c     Calculation of Bessel- and Hankel functions:
c
C***********************************************************************
      subroutine bessel(phipts,phiptsx,mmax,nmmax,
     &          ng1,ng2,ng1dim,ng2dim,thptsx,r1,r2,k,ks,
     &          xth1,xth2,xphi,cone,Im,hn1,wn1,jnn1,unn1,
     &          hn2,wn2,jnn2,unn2,eimphi,sth1,sth2)     
      implicit none
      
      integer i,j,n,m,ng1,ng2,ng1dim,ng2dim,thptsx
      integer phipts,phiptsx,mmax,nmmax
      
      real*8 k,r1(ng1dim),r2(ng2dim,phipts),
     & xth1(ng1dim),xth2(ng2dim,phipts),
     & xphi(phipts),
     & cth,sth1(ng1dim),sth2(thptsx,phiptsx),
     & jtmp(nmmax),ytmp(nmmax),utmp(nmmax),vtmp(nmmax)
     
      complex*16 ks,hn1(ng1dim,nmmax),wn1(ng1dim,nmmax),
     & hn2(phiptsx,thptsx,nmmax),wn2(phiptsx,thptsx,nmmax),
     & cjtmp(nmmax),cutmp(nmmax),jnn1(ng1dim,nmmax),
     & unn1(ng1dim,nmmax),jnn2(phiptsx,thptsx,nmmax),
     & unn2(phiptsx,thptsx,nmmax),
     & eimphi(phiptsx,-2*mmax:2*mmax),
     & cone,Im
            
      if (ng2dim.lt.thptsx) then
         write(6,*)'>>>ERROR: ng2dim.lt.thptsx in subroutine bessel'
         stop
      end if

      do i = 1, ng1dim
         cth = dcos(xth1(i))
         sth1(i) = dsqrt(1-cth**2)
         call cjb(ks*r1(i),cjtmp,cutmp,nmmax)
         call rjb(k*r1(i),jtmp,utmp,nmmax)
         call ryb(k*r1(i),ytmp,vtmp,nmmax)
         r1(i)=r1(i)**2
         do n = 1, nmmax
            jnn1(i,n) = cjtmp(n)
            unn1(i,n) = cutmp(n)
            hn1(i,n) = dcmplx( jtmp(n), ytmp(n) )
            wn1(i,n) = dcmplx( utmp(n), vtmp(n) )
	 end do
      end do

      do i = 1, thptsx
         do j = 1, phiptsx
            cth = dcos(xth2(i,j))
            sth2(i,j) = dsqrt(1-cth**2)
            call cjb(ks*r2(i,j),cjtmp,cutmp,nmmax)
            call rjb(k*r2(i,j),jtmp,utmp,nmmax)
            call ryb(k*r2(i,j),ytmp,vtmp,nmmax)
            r2(i,j)=r2(i,j)**2
            do n = 1, nmmax
               jnn2(j,i,n) = cjtmp(n)
               unn2(j,i,n) = cutmp(n)
               hn2(j,i,n) = dcmplx( jtmp(n), ytmp(n) )
               wn2(j,i,n) = dcmplx( utmp(n), vtmp(n) )
	    end do
	 end do
      end do
     
      do m = 1, 2*mmax
         do j = 1, phiptsx
            eimphi(j,m) = exp( Im*m*xphi(j) )
            eimphi(j,-m) = conjg( eimphi(j,m) )
         end do
      end do

      do j = 1, phiptsx
         eimphi(j,0) = cone
      end do

      return
      end


c***********************************************************************      
c     set up names of output files 
c***********************************************************************
      subroutine foutname(fout,cout,eout,dout,nr)
      implicit none

      character*10 fout, cout, eout,dout
      character*6 c
      integer nr

      if (nr.lt.1)then
         fout='F.dat'
         cout='C.dat'
         eout='E.dat'
         dout='D.dat'
      else
         call i2c(nr,c)
         fout='F'//c(1:6)
         cout='C'//c(1:6)
         dout='D'//c(1:6)
      endif

      end


c***********************************************************************
c     convert integer to character string
c***********************************************************************
      subroutine i2c(i,c)
      implicit none

      integer i,i0
      character*6 c
      character nr(10)
      data nr/'0','1','2','3','4','5','6','7','8','9'/
      
      if(i.gt.999999)then
         call errmsg('ERROR in i2c: i too high')
      endif

      i0=i
      c(1:1)=nr(i0/100000+1)
      i0=i0-(i0/100000)*100000
      c(2:2)=nr(i0/10000+1)
      i0=i0-(i0/10000)*10000
      c(3:3)=nr(i0/1000+1)
      i0=i0-(i0/1000)*1000
      c(4:4)=nr(i0/100+1)
      i0=i0-(i0/100)*100
      c(5:5)=nr(i0/10+1)
      i0=i0-(i0/10)*10
      c(6:6)=nr(i0+1)

      end


c***********************************************************************
c     print warnings and error messages
c***********************************************************************
      subroutine errmsg(msg,fatal)
      implicit none
      character*(*) msg
      logical fatal

      if(fatal) then
         write(6,1001) msg
         stop
      else
         write(6,1000) msg
      end if

 1000 format(">>> WARNING: ",a)
 1001 format(">>> ERROR: ",a)
      end


c***********************************************************************
c     Skip comment lines when reading file
c***********************************************************************
      subroutine skip_lines(unit,nlines)
      implicit none
      integer unit,nlines,i

      do i=1,nlines
         read(unit,*)
      enddo

      end



C***********************************************************************
C TIMING RESULTS 
C***********************************************************************

      SUBROUTINE TIMEIT(MSG,DTIME)
      implicit none

      REAL*8 DTIME
      CHARACTER MSG*42
      
! External system calls:

!     EXTERNAL CPU_TIME

! Local variables:

      CHARACTER :: CSTA*3,CMSGNM*70
      REAL :: T1,T2

! Local variables to be saved:

      SAVE CSTA,T1

! Data statements:

      DATA CSTA/'ONE'/

      IF(CSTA.EQ.'ONE')THEN
         CSTA='TWO'
         CALL CPU_TIME(T1)
      ELSEIF(CSTA.EQ.'TWO')THEN
         CSTA='ONE'
         CALL CPU_TIME(T2)
         DTIME=T2-T1
          IF(DTIME.GT.1.E4)THEN
            WRITE(6,1000) msg,DTIME,' = CPU time (sec)'
         ELSE
            WRITE(6,2000) msg,DTIME,' = CPU time (sec)'
         ENDIF
      ENDIF

 1000 format(a,f9.0,a)
 2000 format(a,f9.3,a)

      RETURN
      END
      
      


C***********************************************************************
C Generic TIMING RESULTS 
c (works on all machines, but does not give out timing results)
C***********************************************************************

      SUBROUTINE LOGMSG(MSG)
      
      CHARACTER MSG*36
      
      WRITE(6,*) msg

      RETURN
      END



c***********************************************************************      
c     Quick-check of input parameter file 
c***********************************************************************
      subroutine check_version(Tver,unit,nlines)
      implicit none
      
      real*8 Tver
      integer unit, nlines

      write(6,2000)Tver
      do while(.true.)
         read(unit,*,end=1000)
         nlines=nlines+1
      enddo
      
 1000 if(abs(Tver-5.2).lt.1e-6 .and. nlines.ne.204)then
         write(6,*)
     &        "ERROR: Format of params file does not match Tsym version"
         stop
      endif

 2000 format("Tsym Version ",f5.1)
      end

