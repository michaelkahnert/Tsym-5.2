c**************************************************************************
c$$$      BSD License for PROPACK numerical software
c$$$      http://www.opensource.org/licenses/bsd-license.php
c$$$
c$$$      This notice applies to the software packages
c$$$
c$$$      PROPACK77 (Fortran 77) and PROPACK (Matlab)
c$$$      
c$$$      made available at
c$$$
c$$$      http://sun.stanford.edu/~rmunk/PROPACK
c$$$
c$$$      by 
c$$$
c$$$      Solar Oscillations Investigation,
c$$$      W.W. Hansen Experimental Physics Laboratory (HEPL),
c$$$      452  Lomita Mall,
c$$$      Stanford University,
c$$$      Stanford, CA 94305-4085, USA
c$$$
c$$$
c$$$      Copyright (c) 2005, Rasmus Munk Larsen, Stanford University
c$$$      All rights reserved.
c$$$
c$$$      Redistribution and use in source and binary forms, with or without
c$$$      modification, are permitted provided that the following conditions are
c$$$      met:
c$$$
c$$$      * Redistributions of source code must retain the above copyright
c$$$      notice, this list of conditions and the following disclaimer.
c$$$
c$$$      * Redistributions in binary form must reproduce the above
c$$$      copyright notice, this list of conditions and the following
c$$$      disclaimer in the documentation and/or other materials provided
c$$$      with the distribution.
c$$$
c$$$      * Neither the name of Stanford University nor the names of its
c$$$      contributors may be used to endorse or promote products derived
c$$$      from this software without specific prior written permission.
c$$$
c$$$      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
c$$$      "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
c$$$      LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
c$$$      A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
c$$$      OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
c$$$      SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
c$$$      LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
c$$$      DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
c$$$      THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
c$$$      (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
c$$$      OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c**************************************************************************


      subroutine zlansvd_irl(which,jobu,jobv,m,n,dim,p,neig,maxiter,
     c     aprod,U,ldu,Sigma,bnd,V,ldv,tolin,work,lwork,zwork,lzwrk,
     c     iwork,liwork,doption,ioption,info,zparm,iparm)

c     DLANSVD_IRL: Compute the leading singular triplets of a large
c     and sparse matrix A by implicitly restarted Lanczos bidiagonalization 
c     with partial reorthogonalization.
c
c     Parameters:
c
c     WHICH: CHARACTER*1. Decides which singular triplets to compute. 
c            If WHICH.EQ.'L' then compute triplets corresponding to the K
c            largest singular values. 
c            If WHICH.EQ.'S' then compute triplets corresponding to the K
c            smallest singular values. 
c     JOBU: CHARACTER*1. If JOBU.EQ.'Y' then compute the left singular vectors.
c           Otherwise the array U is not touched.
c     JOBV: CHARACTER*1. If JOBV.EQ.'Y' then compute the right singular 
c           vectors. Otherwise the array V is not touched.
c     M:    INTEGER. Number of rows of A.
c     N:    INTEGER. Number of columns of A.
c     DIM:  INTEGER. Dimension of the Krylov subspace.
c     P:    INTEGER. Number of shift per restart.
c     NEIG: INTEGER. Number of desired singular triplets.  
c           NEIG <= MIN(DIM-P,M,N)
c     MAXITER: INTEGER. Maximum number of restarts.
c     APROD: Subroutine defining the linear operator A. 
c            APROD should be of the form:
c
c           SUBROUTINE DAPROD(TRANSA,M,N,X,Y,DPARM,IPARM)
c           CHARACTER*1 TRANSA
c           INTEGER M,N,IPARM(*)
c           DOUBLE COMPLEX X(*),Y(*),DPARM(*)
c
c           If TRANSA.EQ.'N' then the function should compute the matrix-vector
c           product Y = A * X.
c           If TRANSA.EQ.'T' then the function should compute the matrix-vector
c           product Y = A^T * X.
c           The arrays IPARM and DPARM are a means to pass user supplied
c           data to APROD without the use of common blocks.
c     U(LDU,KMAX+1): DOUBLE COMPLEX array. On return the first K columns of U
c               will contain approximations to the left singular vectors 
c               corresponding to the K largest or smallest (depending on the 
c               value of WHICH)  singular values of A.
c               On entry the first column of U contains the starting vector
c               for the Lanczos bidiagonalization. A random starting vector
c               is used if U is zero.
c     LDU: INTEGER. Leading dimension of the array U. LDV >= M.
c     SIGMA(K): DOUBLE PRECISION array. On return Sigma contains approximation
c               to the K largest or smallest (depending on the 
c               value of WHICH) singular values of A.
c     BND(K)  : DOUBLE PRECISION array. Error estimates on the computed 
c               singular values. The computed SIGMA(I) is within BND(I)
c               of a singular value of A.
c     V(LDV,KMAX): DOUBLE COMPLEX array. On return the first K columns of V
c               will contain approximations to the right singular vectors 
c               corresponding to the K largest or smallest (depending on the 
c               value of WHICH) singular values of A.
c     LDV: INTEGER. Leading dimension of the array V. LDV >= N.
c     TOLIN: DOUBLE PRECISION. Desired relative accuracy of computed singular 
c            values. The error of SIGMA(I) is approximately 
c            MAX( 16*EPS*SIGMA(1), TOLIN*SIGMA(I) )
c     WORK(LWORK): DOUBLE PRECISION array. Workspace of dimension LWORK.
c     LWORK: INTEGER. Dimension of WORK. 
c            If JOBU.EQ.'N' and JOBV.EQ.'N' then  LWORK should be at least
c            M + N + 10*KMAX + 2*KMAX**2 + 5 + MAX(M,N,4*KMAX+4).
c            If JOBU.EQ.'Y' or JOBV.EQ.'Y' then LWORK should be at least
c            M + N + 10*KMAX + 5*KMAX**2 + 4 + 
c            MAX(3*KMAX**2+4*KMAX+4, NB*MAX(M,N)), where NB>0 is a block 
c            size, which determines how large a fraction of the work in
c            setting up the singular vectors is done using fast BLAS-3 
c            operation. 
c     ZWORK: DOUBLE COMPLEX array of dimension ????.
c     IWORK: INTEGER array. Integer workspace of dimension LIWORK.
c     LIWORK: INTEGER. Dimension of IWORK. Should be at least 8*KMAX if
c             JOBU.EQ.'Y' or JOBV.EQ.'Y' and at least 2*KMAX+1 otherwise.
c     DOPTION: DOUBLE PRECISION array. Parameters for LANBPRO.
c        doption(1) = delta. Level of orthogonality to maintain among
c          Lanczos vectors.
c        doption(2) = eta. During reorthogonalization, all vectors with
c          with components larger than eta along the latest Lanczos vector
c          will be purged.
c        doption(3) = anorm. Estimate of || A ||.
c        doption(4) = min relgap. Smallest relgap allowed between any shift
c                     the smallest requested Ritz value.
c
c     IOPTION: INTEGER array. Parameters for LANBPRO.
c        ioption(1) = CGS.  If CGS.EQ.1 then reorthogonalization is done
c          using iterated classical Gram-Schmidt. IF CGS.EQ.0 then 
c          reorthogonalization is done using iterated modified Gram-Schmidt.
c        ioption(2) = ELR. If ELR.EQ.1 then extended local orthogonality is
c          enforced among u_{k}, u_{k+1} and v_{k} and v_{k+1} respectively.
c
c     INFO: INTEGER. 
c         INFO = 0  : The K largest or smallest (depending on the 
c                     value of WHICH) singular triplets were computed 
c                     succesfully.
c         INFO = J>0, J<K: An invariant subspace of dimension J was found.
c         INFO = -1 : K singular triplets did not converge within KMAX
c                     iterations.   
c     ZPARM: DOUBLE COMPLEX array. Array used for passing data to the APROD
c         function.   
c     IPARM: INTEGER array. Array used for passing data to the APROD
c         function.   
c
c
c     (C) Rasmus Munk Larsen, Stanford University, 2000,2004
c


c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      character*1 which,jobu,jobv
      integer m,n,p,neig,maxiter,ldu,ldv,iter,liwork,lzwrk
      integer iwork(liwork),lwork,info,ioption(*)
      complex*16 U(ldu,*),V(ldv,*), zwork(lzwrk),zparm(*)
      double precision Sigma(*),bnd(*),work(lwork)
      double precision tolin,doption(*)
      integer iparm(*)
      external aprod

c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero, FUDGE
      parameter(one = 1.0, zero = 0.0, FUDGE = 1.01)
            
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer k,i,ibnd,iwrk,ierr,ip,iq,nconv,lwrk,kold,dim
      integer ialpha,ibeta,ialpha1,ibeta1,ishift,nshft,lapinfo
      double precision eps,eps34,epsn2,epsn,sfmin,anorm,rnorm,tol
      double precision shift, relgap
      real t0,t1,t2,t3
      integer st,cnk,wst,wcnk, tid, nt

c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external dzero,izero,dcopy,daxpy,dbdsqr
      
c     %--------------------%
c     | External Functions |
c     %--------------------%
      logical lsame
      double precision dlamch,pdznrm2
      external pdznrm2,lsame
      external dlamch
#ifdef _OPENMP
      integer omp_get_thread_num, omp_get_num_threads
      external omp_get_thread_num, omp_get_num_threads
#endif

c-------------------- Here begins executable code ---------------------
      
c     %-------------%
c     | Start timer |
c     %-------------%
      call second(t0)

c     %---------------------------------%
c     | Set machine dependent constants |
c     %---------------------------------%
      eps = dlamch('e')
      eps34 = eps**(3.0/4.0)
      epsn = dble(max(m,n))*eps/2.0
      epsn2 = sqrt(dble(max(m,n)))*eps/2.0
      sfmin = dlamch('s')
      

c     %--------------------------------%
c     | Guard against absurd arguments |
c     %--------------------------------%
      dim = min(dim,n+1,m+1)
      k = dim-p
      tol = min(one,max(16.0*eps,tolin))
      anorm = zero

c     %------------------------------%
c     | Set pointers into work array |
c     %------------------------------%
      ibnd = 1
      ialpha = ibnd + dim+1
      ibeta = ialpha + dim
      ialpha1 = ibeta + dim
      ibeta1 = ialpha1 + dim
      ishift = ibeta1 + dim
      ip = ishift + dim
      iq = ip + (dim+1)**2
      iwrk = iq + dim**2
      lwrk = lwork-iwrk+1
      call dzero(8*dim + 3 + 2*dim**2,work,1)

c     %---------------------------------------------------------------%
c     | Set up random starting vector if none is provided by the user |
c     %---------------------------------------------------------------%
      rnorm = pdznrm2(m,U(1,1),1)
      if (rnorm.eq.zero) then
         call zgetu0('n',m,n,0,1,U,rnorm,U,ldu,aprod,
     c        zparm,iparm, ierr,ioption(1),anorm,zwork)     
      endif

      iter = 0
      info = 0
      nconv = 0
      kold = 0       
c     %------------------------------%
c     | Iterate until convergence... |
c     %------------------------------%
      do while (nconv.lt.neig .and. iter.lt.maxiter)


c     %---------------------------------------------------%
c     | Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j} |
c     %---------------------------------------------------%

         call zlanbpro(m, n, kold, dim, aprod, U, ldu, V, ldv,
     c        work(ialpha),dim,rnorm,doption(1),ioption(1),
     c        work(iwrk), zwork, iwork, zparm, iparm, ierr)         
         kold = k
c     %---------------------------------------------%
c     | Compute and analyze SVD(B) and error bounds |
c     %---------------------------------------------%
         call dcopy(dim, work(ialpha),1,work(ialpha1),1)
         call dcopy(dim, work(ibeta),1,work(ibeta1),1)
         call dzero(dim+1,work(ibnd),1)

         call second(t2)
         call dbdqr((dim.eq.min(m,n)),'N',dim,work(ialpha1),
     c        work(ibeta1),work(ibnd+dim-1),work(ibnd+dim),
     c        work(ip),dim+1)
         
         call dbdsqr('u',dim,0,1,0,work(ialpha1),work(ibeta1),work,1,
     c        work(ibnd),1,work,1,work(iwrk),lapinfo)
         call  second(t3)
         tbsvd = tbsvd + (t3-t2)
         nbsvd = nbsvd + 1

         if (dim.gt.5) then
            anorm = work(ialpha1)
         else
            anorm = max(anorm,work(ialpha1))
         endif
         do i=1,dim
            work(ibnd+i-1) = abs(rnorm*work(ibnd+i-1))
c            write (6,*) 'bnd(',i,') = ',work(ibnd+i-1)
         enddo

c     %---------------------------------------------%
c     | Refine error bounds using the "Gap theorem" |
c     %---------------------------------------------%
         if (lsame(which,'s')) then
            call drefinebounds(min(m,n),dim,work(ialpha1),work(ibnd),
     c           epsn*anorm,eps34) 
         else
            call drefinebounds(min(m,n), min(dim,neig),work(ialpha1),
     c           work(ibnd),epsn*anorm,eps34) 
         endif
c         do i=1,dim
c            write (6,*) 'bnd(',i,') = ',work(ibnd+i-1)
c         enddo
         
c     %----------------------------------------------------%
c     | Determine the number of converged singular values  |
c     %----------------------------------------------------%
c         do i=1,min(dim,neig)
c            write(6,*)  iter,i,work(ialpha1+i-1),work(ibnd+i-1)
c         enddo
         if (lsame(which,'s')) then
            i = dim-neig+1
            nconv = 0
            do while(i.le.dim)
c               write(6,*) 'iter = ',iter,'sigma = ',work(ialpha1+i-1),
c     c              'error = ',work(ibnd+i-1)
c               write(6,*)  iter,work(ialpha1+i-1),work(ibnd+i-1)
               if (work(ibnd+i-1).le.tol*work(ialpha1)) then
                  nconv = nconv + 1
                  sigma(nconv) = work(ialpha1+i-1)
                  bnd(nconv) = work(ibnd+i-1)
               endif
               i = i+1
            enddo
         else
            i = 1
            nconv = 0
            do while(i.le.min(dim,neig))
c               write(6,*)  iter,work(ialpha1+i-1),work(ibnd+i-1)
               if (work(ibnd+i-1).le.tol*work(ialpha1+i-1)) then
                  nconv = nconv + 1
                  sigma(nconv) = work(ialpha1+i-1)
                  bnd(nconv) = work(ibnd+i-1)
                  i = i+1
               else
                  i = k+1
               endif
            enddo
         endif
            
         
c     %-----------------------------------------------%
c     | Test if an invariant subspace have been found |
c     %-----------------------------------------------%
         if (ierr.lt.0) then
            if (dim.lt.k) then
               write(6,*) 'WARNING: Invariant subspace found.',
     c              ' Dimension = ',dim
               info = dim
            endif
            goto 50               
         endif
         if (nconv.lt.neig) then

c     %---------------------------------------------------------------------%
c     | Implicit restart:                                                   |
c     |   Apply  shifts mu_1, mu_2,...,mu_p to "back up" the          |
c     |   bidiagonalization (dim-k) steps to                                |
c     |                                                                     |
c     |      A V_{k}^{+} = U_{k+1}^{+} B_{k}^{+}                            |
c     | corresponding to starting vector                                    |
c     |      u_1^{+} = \prod_{i=1}^{p} (A A^T - mu_i^2) u_1             |
c     |                                                                     |
c     | We use exact shifts mu_i for which the relative gap between mu_i    |
c     | and the lower bound on the k'th Ritzvalue is larger than doption(4) |
c     %---------------------------------------------------------------------%
            call second(t2)
            call dzero(dim-k,work(ishift),1)
            nshft = 0
            if (lsame(which,'s')) then
               do i=1,k
                  relgap = (work(ialpha1+i-1)-work(ibnd+i-1) -
     c                 work(ialpha1+dim-neig-1))
                  if (relgap.gt.doption(4)*work(ialpha1+dim-neig-1))
     c                 then
                     work(ishift + nshft) = work(ialpha1+i-1)
                  else
                     work(ishift + nshft) = work(ialpha1)
                  endif
c                  write (6,*) 'shift = ',work(ishift + nshft)
                  nshft = nshft + 1
               enddo      
            else
               do i=dim,k+1,-1
                  relgap = work(ialpha1+k-1) -
     c                 (work(ialpha1+i-1)+work(ibnd+i-1))
                  if (relgap.gt.doption(4)*work(ialpha1+k-1)) then
                     work(ishift + nshft) = work(ialpha1+i-1)
                  else
                     work(ishift + nshft) = 0d0
                  endif
                  nshft = nshft + 1
               enddo      
            endif

c     %--------------------------------------------------%
c     | Apply shifts and accumulate rotations such that  |
c     |   B_{dim}^{+} = P * B_{dim} * Q^T                |
c     %--------------------------------------------------%
            call dzero((dim+1)*(dim+1),work(ip),1)
            call dzero(dim*dim,work(iq),1)
            do i=1,dim+1
               work(ip+(i-1)*(dim+2)) = one
            enddo
            do i=1,dim
               work(iq+(i-1)*(dim+1)) = one
            enddo

            do i=dim,k+1,-1
               shift = work(ishift+dim-i)
               call dbsvdstep('y','y',dim+1,dim,i,shift,work(ialpha),
     c              work(ibeta),work(ip),dim+1,work(iq), dim)
            enddo

c     %---------------------------------------------------%
c     | Compute first k+1 left and first k right updated  |
c     | Lanczos vectors                                   |
c     |   U_{dim+1}^{+} = U_{dim+1} * P(:,1:k+1)          |
c     |   V_{dim}^{+} = V_{dim} * Q(:,1:k)                |
c     %---------------------------------------------------%
c$OMP PARALLEL private(tid,nt,cnk,st,wcnk,wst) 
#ifdef _OPENMP
            tid = omp_get_thread_num()
            nt = omp_get_num_threads()
#else
            tid = 0
            nt = 1
#endif
            wcnk = lzwrk/nt
            wst = tid*wcnk+1
            cnk = m/nt
            st = tid*cnk+1
            if (tid.eq.nt-1) then
               wcnk = lzwrk-wst+1
               cnk = m-st+1
            endif
            call zdgemm_ovwr_left('n',cnk,k+1,dim+1,U(st,1),ldu,
     c           work(ip),dim+1,zwork(wst),wcnk)            
            cnk = n/nt
            st = tid*cnk+1
            if (tid.eq.nt-1) then
               cnk = n-st+1
            endif
            call zdgemm_ovwr_left('n',cnk,k,dim,V(st,1),ldv,work(iq),
     c           dim,zwork(wst),wcnk)
c$OMP END PARALLEL
            rnorm = work(ibeta+k-1)
            call second(t3)
            trestart = trestart + (t3-t2)
            nrestart = nrestart + 1
         endif         
         iter = iter + 1
      enddo

 50   if ((nconv.ge.neig .or. info.gt.0) .and. 
     c     (lsame(jobu,'y') .or. lsame(jobv,'y'))) then
c     %-----------------------------------------%
c     | Calculate singular vectors if requested %
c     %-----------------------------------------%
         call dcopy(dim, work(ialpha),1,work(ialpha1),1)
         call dcopy(dim, work(ibeta),1,work(ibeta1),1)
         lwrk = lwrk + dim**2 + (dim+1)**2
         call zritzvec(which, jobu,jobv,m,n,nconv,dim,work(ialpha1),
     c        work(ibeta1),work(ialpha1),U,ldu,V,ldv,work(ip),
     c        lwrk,zwork,lzwrk,iwork)
      endif
      neig = nconv
      nlandim = dim
      call second(t1)
      tlansvd = t1-t0
      end



      subroutine zlansvd(jobu,jobv,m,n,k,kmax,aprod,U,ldu,Sigma,bnd,
     c     V,ldv,tolin,work,lwork,zwork,lzwrk,iwork,liwork,doption,
     c     ioption,info,zparm,iparm)

c     DLANSVD: Compute the leading singular triplets of a large and
c     sparse matrix by Lanczos bidiagonalization with partial
c     reorthogonalization.
c
c     Parameters:
c
c     JOBU: CHARACTER*1. If JOBU.EQ.'Y' then compute the left singular vectors.
c           Otherwise the array U is not touched.
c     JOBV: CHARACTER*1. If JOBV.EQ.'Y' then compute the right singular 
c           vectors. Otherwise the array V is not touched.
c     M: INTEGER. Number of rows of A.
c     N: INTEGER. Number of columns of A.
c     K: INTEGER. Number of desired singular triplets. K <= MIN(KMAX,M,N)
c     KMAX: INTEGER. maximal number of iterations / maximal dimension of
c           Krylov subspace.
c     APROD: Subroutine defining the linear operator A. 
c            APROD should be of the form:
c
c           SUBROUTINE DAPROD(TRANSA,M,N,X,Y,ZPARM,IPARM)
c           CHARACTER*1 TRANSA
c           INTEGER M,N,IPARM(*)
c           DOUBLE PRECISION X(*),Y(*),ZPARM(*)
c
c           If TRANSA.EQ.'N' then the function should compute the matrix-vector
c           product Y = A * X.
c           If TRANSA.EQ.'T' then the function should compute the matrix-vector
c           product Y = A^T * X.
c           The arrays IPARM and ZPARM are a means to pass user supplied
c           data to APROD without the use of common blocks.
c     U(LDU,KMAX+1): DOUBLE COMPLEX array. On return the first K columns of U
c               will contain approximations to the left singular vectors 
c               corresponding to the K largest singular values of A.
c               On entry the first column of U contains the starting vector
c               for the Lanczos bidiagonalization. A random starting vector
c               is used if U is zero.
c     LDU: INTEGER. Leading dimension of the array U. LDU >= M.
c     SIGMA(K): DOUBLE PRECISION array. On return Sigma contains approximation
c               to the K largest singular values of A.
c     BND(K)  : DOUBLE PRECISION array. Error estimates on the computed 
c               singular values. The computed SIGMA(I) is within BND(I)
c               of a singular value of A.
c     V(LDV,KMAX): DOUBLE COMPLEX array. On return the first K columns of V
c               will contain approximations to the right singular vectors 
c               corresponding to the K largest singular values of A.
c     LDV: INTEGER. Leading dimension of the array V. LDV >= N.
c     TOLIN: DOUBLE PRECISION. Desired relative accuracy of computed singular 
c            values. The error of SIGMA(I) is approximately 
c            MAX( 16*EPS*SIGMA(1), TOLIN*SIGMA(I) )
c     WORK(LWORK): DOUBLE PRECISION array. Workspace of dimension LWORK.
c     LWORK: INTEGER. Dimension of WORK.
c            If JOBU.EQ.'N' and JOBV.EQ.'N' then  LWORK should be at least
c            M + N + 9*KMAX + 2*KMAX**2 + 4 + MAX(M,N,4*KMAX+4).
c            If JOBU.EQ.'Y' or JOBV.EQ.'Y' then LWORK should be at least
c            M + N + 9*KMAX + 5*KMAX**2 + 4 + 
c            MAX(3*KMAX**2+4*KMAX+4, NB*MAX(M,N)), where NB>0 is a block 
c            size, which determines how large a fraction of the work in
c            setting up the singular vectors is done using fast BLAS-3 
c            operation. 
c     ZWORK: DOUBLE COMPLEX array of dimension ????.
c     IWORK: INTEGER array. Integer workspace of dimension LIWORK.
c     LIWORK: INTEGER. Dimension of IWORK. Should be at least 8*KMAX if
c             JOBU.EQ.'Y' or JOBV.EQ.'Y' and at least 2*KMAX+1 otherwise.
c     DOPTION: DOUBLE PRECISION array. Parameters for LANBPRO.
c        doption(1) = delta. Level of orthogonality to maintain among
c          Lanczos vectors.
c        doption(2) = eta. During reorthogonalization, all vectors with
c          with components larger than eta along the latest Lanczos vector
c          will be purged.
c        doption(3) = anorm. Estimate of || A ||.
c     IOPTION: INTEGER array. Parameters for LANBPRO.
c        ioption(1) = CGS.  If CGS.EQ.1 then reorthogonalization is done
c          using iterated classical GRAM-SCHMIDT. IF CGS.EQ.0 then 
c          reorthogonalization is done using iterated modified Gram-Schmidt.
c        ioption(2) = ELR. If ELR.EQ.1 then extended local orthogonality is
c          enforced among u_{k}, u_{k+1} and v_{k} and v_{k+1} respectively.
c     INFO: INTEGER. 
c         INFO = 0  : The K largest singular triplets were computed succesfully
c         INFO = J>0, J<K: An invariant subspace of dimension J was found.
c         INFO = -1 : K singular triplets did not converge within KMAX
c                     iterations.   
c     ZPARM: DOUBLE COMPLEX array. Array used for passing data to the APROD
c         function.   
c     IPARM: INTEGER array. Array used for passing data to the APROD
c         function.   
c
c     (C) Rasmus Munk Larsen, Stanford, 1999, 2004 
c


c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      character*1 jobu,jobv
      integer info,liwork
      integer m,n,k,kmax,lanmax,ldu,ldv,iwork(liwork),lwork,lzwrk
      complex*16 U(ldu,*),V(ldv,*), zwork(lzwrk)
      double precision Sigma(*),bnd(*),work(lwork)
      complex*16 zparm(*)
      double precision tolin,doption(*)
      integer iparm(*),ioption(*)
      external aprod

c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero, FUDGE
      parameter(one = 1.0, zero = 0.0, FUDGE = 1.01)
            
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j,dj,jold,ibnd,ib,ib1,iwrk,ierr,ip,iq,neig,lwrk,lapinfo
      double precision eps,eps34,epsn2,epsn,sfmin,anorm,rnorm,tol
      real t0,t1,t2,t3

c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external dzero,izero,dcopy,dbdsqr
      external zzero,zcopy
      
c     %--------------------%
c     | External Functions |
c     %--------------------%
      logical lsame
      double precision dlamch,pdznrm2
      external pdznrm2,lsame
      external dlamch

c-------------------- Here begins executable code ---------------------
      
c     %-------------%
c     | Start timer |
c     %-------------%
      call second(t0)

c     %---------------------------------%
c     | Set machine dependent constants |
c     %---------------------------------%
      eps = dlamch('e')
      eps34 = eps**(3.0/4.0)
      epsn = dble(max(m,n))*eps/2.0
      epsn2 = sqrt(dble(max(m,n)))*eps/2.0
      sfmin = dlamch('s')
      

c     %--------------------------------%
c     | Guard against absurd arguments |
c     %--------------------------------%
      lanmax = min(n+1,m+1,kmax)
      tol = min(one,max(16.0*eps,tolin))
      anorm = zero

c     %------------------------------%
c     | Set pointers into work array |
c     %------------------------------%
      ibnd = 1
      ib = ibnd + lanmax+1
      ib1 = ib + 2*lanmax
      ip = ib1 + 2*lanmax
      iq = ip + (lanmax+1)**2
      iwrk = iq + lanmax**2
      lwrk = lwork-iwrk+1

      call dzero(7*lanmax + 2 + 2*lanmax**2,work,1)
      call zzero(lzwrk,zwork,1)

c     %---------------------------------------------------------------%
c     | Set up random starting vector if none is provided by the user |
c     %---------------------------------------------------------------%
c      print *,(U(i,1),i=1,m)
      rnorm = pdznrm2(m,U(1,1),1)
c      print *,'rnorm = ',rnorm
      if (rnorm.eq.zero) then
         call zgetu0('n',m,n,0,1,U,rnorm,U,ldu,aprod,
     c        zparm,iparm, ierr,ioption(1),anorm,zwork)     
      endif
c      print *,'U = '
c      print *,(U(i,1),i=1,m)
      
      nsing = k
      info = 0
      neig = 0
      jold = 0
      j = min(k+max(8,k)+1,lanmax)
       
c     %------------------------------%
c     | Iterate until convergence... |
c     %------------------------------%
      do while (neig.lt.k)

c     %---------------------------------------------------%
c     | Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j} |
c     %---------------------------------------------------%
         call zlanbpro(m, n, jold, j, aprod, U, ldu, V, ldv,
     c        work(ib),lanmax,rnorm,doption(1),ioption(1),
     c        work(iwrk), zwork, iwork, zparm, iparm, ierr)
         jold = j

c     %---------------------------------------------%
c     | Compute and analyze SVD(B) and error bounds |
c     %---------------------------------------------%
         call dcopy(2*lanmax, work(ib),1,work(ib1),1)
         call dzero(j+1,work(ibnd),1)
         
         call second(t2)
         call dbdqr((j.eq.min(m,n)),'N',j,work(ib1),work(ib1+lanmax),
     c        work(ibnd+j-1),work(ibnd+j),work(ip),lanmax+1)
         call dbdsqr('u',j,0,1,0,work(ib1),work(ib1+lanmax),work,1,
     c        work(ibnd),1,work,1,work(iwrk),lapinfo)
         call  second(t3)
         tbsvd = tbsvd + (t3-t2)
         nbsvd = nbsvd + 1

         if (j.gt.5) then
            anorm = work(ib1)
         else
            anorm = max(anorm,work(ib1))
         endif
         do i=1,j
            work(ibnd+i-1) = abs(rnorm*work(ibnd+i-1))
         enddo

c     %---------------------------------------------%
c     | Refine error bounds using the "Gap theorem" |
c     %---------------------------------------------%
c         if (lsame(jobu,'n') .and. lsame(jobv,'n')) then
            call drefinebounds(min(m,n),j,work(ib1),work(ibnd),
     c           epsn*anorm,eps34)
c         endif

c     %----------------------------------------------------%
c     | Determine the number of converged singular values  |
c     %----------------------------------------------------%
         do i=1,min(j,k)
            bnd(i) = work(ibnd+i-1)
c            write(6,*) 'sigma,bnd = ',work(ib1+i),bnd(i)
         enddo
         i = 0
         neig = 0
         do while(i.lt.min(j,k))
            if (work(ibnd+i).le.tol*work(ib1+i)) then
               neig = neig + 1
               sigma(neig) = work(ib1+i)
               i = i+1
            else
               i = k
            endif
         enddo
         
c     %--------------------------------------------------%
c     | Test if an invariant subspace has been found or |
c     | if the workspace has been exhausted.             |
c     %--------------------------------------------------%
         if (ierr.lt.0) then
            if (j.lt.k) then
               write(6,*) 'WARNING: Invariant subspace found.',
     c              ' Dimension = ',j
               info = j
            endif
            goto 50               
         endif
         if (j.ge.lanmax) then
            if (neig.lt.k) then
               write(6,*) 'WARNING: Maximum dimension of Krylov',
     c              ' subspace exceeded prior to convergence.',
     c              ' Try increasing KMAX.'
               write(6,*) 'neig = ',neig
               info = -1
            endif
            goto 50
         endif

c     %----------------------------------------------------%
c     | Increase the dimension of the Krylov subspace.     |
c     | If any Ritz values have converged then try to      | 
c     | estimate the average number of iterations per      |
c     | converged Ritz value.                              |
c     | Else increase the dimension by 50%.                |
c     %----------------------------------------------------%
         if (neig.gt.1) then
            dj = min(j/2,((k-neig)*(j-6))/(2*neig+1))
            dj = min(100,max(2,dj))
         else
            dj = j/2
            dj = min(100,max(10,dj))
        endif
         j = min(j + dj,lanmax)
      enddo

 50   if ((neig.ge.k .or. info.gt.0)  .and. 
     c     (lsame(jobu,'y') .or. lsame(jobv,'y'))) then
c     %-----------------------------------------%
c     | Calculate singular vectors if requested %
c     %-----------------------------------------%
c         print *,'computing vectors: neig = ',neig,' jold = ',jold
         lwrk = lwrk + lanmax**2 + (lanmax+1)**2
         call zritzvec('L',jobu,jobv,m,n,neig,jold,work(ib),
     c        work(ib+lanmax),work(ib1),U,ldu,V,ldv,work(ip),
     c        lwrk,zwork,lzwrk,iwork)
      endif
      k = neig
      nlandim = j
      call second(t1)
      tlansvd = t1-t0
      end
c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

#ifdef _OPENMP
      double precision function pdnrm2(n, x, incx)
      implicit none
      integer n, incx
      double precision x(*)
      
      integer i
      double precision sum

      if ((n.gt.0).and.(incx.ne.0)) then    
         sum = 0d0
         if (incx.eq.1) then
c$OMP PARALLEL DO  reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + x(i)**2
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + x(1+(i-1)*incx)**2
            enddo
         endif
         pdnrm2 = sqrt(sum)
      else
         pdnrm2 = 0d0
      endif   
      return
      end
c
c****************************************************************************
c 
      

      subroutine pdscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      double precision alpha,x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
            do i=1,n
               x(i) = alpha*x(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(alpha, incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = alpha*x(1+(i-1)*incx)
            enddo
         endif
      endif
      return
      end
      
c
c****************************************************************************
c 
         

      subroutine pdcopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision x(*),y(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then         
         if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               y(i) = x(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx, incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = x(1+(i-1)*incx)
            enddo
         endif
      endif
      return
      end

c
c****************************************************************************
c 
      subroutine pdaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision alpha,x(*),y(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then         
         if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha)  schedule(static)
           do i=1,n
               y(i) = alpha*x(i) + y(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(alpha,incx,incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c              y(1+(i-1)*incy)
            enddo
         endif
      endif
      return
      end
      
c
c****************************************************************************
c 
     
         
      double precision function pddot(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision x(*),y(*)
      
      integer i
      double precision sum

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then    
         if (incx.eq.1 .and. incy.eq.1) then
            sum = 0d0
c$OMP PARALLEL DO reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + x(i) * y(i)
            enddo
         else
            sum = 0d0
c$OMP PARALLEL DO firstprivate(incx, incy) reduction(+:sum) 
c$OMP& schedule(static) 
            do i=1,n
               sum = sum + x(1+(i-1)*incx) * y(1+(i-1)*incy)
            enddo
         endif
         pddot = sum
      else
         pddot = 0d0
      endif   
      return
      end


#else

      double precision function pdnrm2(n, x, incx)
      implicit none
      integer n, incx
      double precision x(*), dnrm2
      external dnrm2

      pdnrm2 = dnrm2(n, x, incx)
      end

      subroutine pdscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      double precision alpha,x(*)

      call dscal(n, alpha, x , incx)
      end

      subroutine pdcopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision x(*),y(*)

      call dcopy(n, x , incx, y, incy)
      end

      subroutine pdaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision alpha,x(*),y(*)
     
      call daxpy(n, alpha, x , incx, y, incy)
      end

      
      double precision function pddot(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision x(*),y(*),ddot
      external ddot
      
      pddot = ddot(n, x , incx, y, incy)
      end
#endif

c
c****************************************************************************
c 
         
      subroutine pdzero(n, x , incx)
      implicit none
      integer n, incx
      double precision x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO 
            do i=1,n
               x(i) = 0d0
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = 0d0
            enddo
         endif
      endif
      return
      end


      subroutine pizero(n, x , incx)
      implicit none
      integer n, incx
      integer x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               x(i) = 0
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = 0
            enddo
         endif
      endif
      return
      end


      subroutine pdset(n, alpha, x , incx)
      implicit none
      integer n, incx
      double precision alpha,x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
            do i=1,n
               x(i) = alpha
            enddo
         else
c$OMP PARALLEL DO firstprivate(alpha, incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = alpha
            enddo
         endif
      endif
      return
      end


  
      subroutine pdaxpby(n,alpha,x,incx,beta,y,incy)
c
c     Y = alpha*X + beta*Y
c     

      implicit none
      double precision one,zero
      parameter(one = 1d0,zero = 0d0)
      integer n,incx,incy,i
      double precision alpha,beta,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.zero .and. beta.eq.zero) then
         if (incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               y(i) = zero
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = zero
            enddo
         endif
         
      else if (alpha.eq.zero .and. beta.ne.zero) then
         
         call pdscal(n,beta,y,incy)

      else if (alpha.ne.zero .and. beta.eq.zero) then

         if (alpha.eq.one) then
            call pdcopy(n,x,incx,y,incy)
         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
               do i=1,n
                  y(i) = alpha*x(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(incx, incy, alpha) 
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)
               enddo
            endif
         endif

      else

         if (beta.eq.one) then
c DAXPY
            call pdaxpy(n,alpha,x,incx,y,incy)
         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha,beta) 
c$OMP& schedule(static) 
               do i=1,n
                  y(i) = alpha*x(i) + beta*y(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(alpha,beta,incx,incy)  
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c                 beta*y(1+(i-1)*incy)
               enddo
            endif
         endif
      endif
      return
      end



      subroutine pdaxty(n,alpha,x,incx,y,incy)
c
c     Y = alpha*X*Y
c     

      implicit none
      double precision one,zero
      parameter(one = 1d0,zero = 0d0)
      integer n,incx,incy,i
      double precision alpha,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.zero) then
         if (incy.eq.1) then
c$OMP PARALLEL DO schedule(static)
           do i=1,n
               y(i) = zero
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
           do i=1,n
               y(1+(i-1)*incy) = zero
            enddo
         endif
         
      else if (alpha.ne.zero) then

         if (alpha.eq.one) then
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
               do i=1,n
                  y(i) = x(i)*y(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(incx,incy) schedule(static)
               do i=1,n
                  y(1+(i-1)*incy) = x(1+(i-1)*incx)*y(1+(i-1)*incy)
               enddo
            endif

         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
               do i=1,n
                  y(i) = alpha*x(i)*y(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(alpha,incx,incy) 
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)*
     c                 y(1+(i-1)*incy)
               enddo
            endif
         endif
      endif
      return
      end




      subroutine szero(n, x , incx)
      implicit none
      integer n, incx
      real x(*),zero
      parameter (zero = 0)          
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = zero
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = zero
            enddo
         endif
      endif
      return
      end

                 
      subroutine dzero(n, x , incx)
      implicit none
      integer n, incx
      double precision x(*),zero
      parameter (zero = 0.0)          
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = zero
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = zero
            enddo
         endif
      endif
      return
      end


      subroutine izero(n, x , incx)
      implicit none
      integer n, incx
      integer x(*),zero
      parameter (zero = 0)          
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = zero
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = zero
            enddo
         endif
      endif
      return
      end

         
c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

      subroutine dbsvdstep(jobu,jobv,m,n,k,sigma,D,E,U,ldu,V,ldv)
c
c     Perform one implicit LQ SVD sweep with shift SIGMA.
c

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      character*1 jobu,jobv
      integer m,n,k,ldu,ldv
      double precision D(*),E(*),U(ldu,*),V(ldv,*),sigma

c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero
      parameter(one = 1.0, zero = 0.0)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i
      double precision c,s,x,y,r
      logical dou,dov

      logical lsame
      external lsame,dlartg,drot
      
c-------------------- Here begins executable code ---------------------

      if (k.le.1) return

      dou = lsame(jobu,'y')
      dov = lsame(jobv,'y')

c     Compute the initial rotation based on B*B^T-sigma^2
      x = D(1)*D(1) - sigma*sigma
      y = E(1)*D(1)

c     Chase the "bulge" down the lower bidiagonal with Givens rotations.
c     Below 'y' is the "bulge" and 'x' is the element used to eliminate it.
      do i=1,k-1
        if (i.gt.1) then
           call dlartg(x,y,c,s,E(i-1))
        else
           call dlartg(x,y,c,s,r)
        endif
        x = c*D(i) + s*E(i)
        E(i) = -s*D(i) + c*E(i)
        D(i) = x
        y = s*D(i+1)
        D(i+1) = c*D(i+1)

        if (dou .and. m.gt.0) then
           call drot(m,U(1,i),1,U(1,i+1),1,c,s)
        endif

        call dlartg(x,y,c,s,D(i))
        x = c*E(i) + s*D(i+1)
        D(i+1) = -s*E(i) + c*D(i+1)
        E(i) = x
        y = s*E(i+1)
        E(i+1) = c*E(i+1)

        if (dov .and. n.gt.0) then
           call drot(n,V(1,i),1,V(1,i+1),1,c,s)
        endif
      enddo
      call dlartg(x,y,c,s,E(k-1))
      x = c*D(k) + s*E(k)
      E(k) = -s*D(k) + c*E(k)
      D(k) = x
      if (dou .and. m.gt.0) then
         call drot(m,U(1,k),1,U(1,k+1),1,c,s)
      endif
      return
      end




      subroutine dbdqr(ignorelast, jobq, n, D, E, c1, c2, Qt, ldq)
      implicit none

c Compute QR factorization B = Q*R of (n+1) x n lower bidiagonal matrix 
c with diagonal elements d(1)...d(n) and first subdiagonal elements
c e(1)...e(n). On return [0 ... 0 c1 c2]' = Q'*[0 ... 0 1]'.
c If ignorelast.eq..true. then e(n) is assumed to be zero.
c
c If jobq=='Y' then on return Qt contains Q^T.

c     %------------%
c     | Parameters |
c     %------------%
      character*1 jobq
      logical ignorelast
      integer n,ldq
      double precision D(*),E(*),c1,c2,Qt(ldq,*)
      
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j
      double precision cs,sn,r

c     %------------------------------------%
c     | External Functions and Subroutines |
c     %------------------------------------%
      logical lsame
      external lsame

c-------------------- Here begins executable code ---------------------

      if (n.lt.1) return
      if (lsame(jobq,'Y')) then
         do j=1,n+1
            do i=1,n+1
               Qt(i,j) = 0.0
            enddo
            Qt(j,j) = 1.0
         enddo
      endif
      do i=1,n-1
         call dlartg(d(i),e(i),cs,sn,r)
         d(i) = r
         e(i) = sn*d(i+1)
         d(i+1) = cs*d(i+1)
         if (lsame(jobq,'Y')) then
            do j=1,i
               Qt(i+1,j) = -sn*Qt(i,j)
               Qt(i,j) = cs*Qt(i,j)
            enddo
            Qt(i,i+1) = sn
            Qt(i+1,i+1) = cs
         endif
      enddo
      if (.not.ignorelast) then
         call dlartg(d(n),e(n),cs,sn,r)
         d(n) = r
         e(n) = 0.0
         c1 = sn
         c2 = cs
         if (lsame(jobq,'Y')) then
            do j=1,i
               Qt(i+1,j) = -sn*Qt(i,j)
               Qt(i,j) = cs*Qt(i,j)
            enddo
            Qt(i,i+1) = sn
            Qt(i+1,i+1) = cs
         endif
      endif
      end




      subroutine drefinebounds(n,k,theta,bound,tol,eps34)
c
c     Refine Lanczos error bounds using the gap theorem.
c     
c     Input arguments: 
c              n:     smallest dimension of original matrix
c              k:     number of Ritz values to refine
c              theta: array of Ritz values
c              bound: array of unrefined error bounds
c              tol:   clustering tolerance
c              eps34: machine epsilon to the power 3/4.

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer n,k
      double precision theta(*), bound(*), tol, eps34

c     %-----------------%
c     | Local variables |
c     %-----------------%
      double precision gap
      integer i,l

c     %------------------------------------%
c     | External Functions and Subroutines |
c     %------------------------------------%
      double precision dlapy2
      external dlapy2

c-------------------- Here begins executable code ---------------------
      if (k.le.1) return
      do i=1,k
         do l=-1,1,2
            if ((l.eq.1.and.i.lt.k) .or. (l.eq.-1.and.i.gt.1)) then
               if (abs(theta(i)-theta(i+l)) .lt. eps34*(theta(i))) then
                  if (bound(i).gt.tol .and. bound(i+l).gt.tol) then
                     bound(i+l) = dlapy2(bound(i),bound(i+l))
                     bound(i) = 0.0
                  endif
               endif
            endif
         enddo         
      enddo
      do i=1,k
         if (i.lt.k .or. k.eq.n) then
c
c     We cannot compute a reliable value for the gap of the last
c     Ritz value unless we know it is an approximation to the 
c     smallest singular value (k.eq.n). In this case we can take the 
c     distance to the next bigger one as the gap, which can really 
c     save us from getting stuck on matrices with a single isolated tiny 
c     singular value.
c
            if (i.eq.1) then
               gap = abs(theta(i)-theta(i+1))-max(bound(i),bound(i+1))
            else if (i.eq.n) then
               gap = abs(theta(i-1)-theta(i))-max(bound(i-1),bound(i))
            else
               gap = abs(theta(i)-theta(i+1))-max(bound(i),bound(i+1)) 
               gap = min(gap,abs(theta(i-1) - theta(i)) - 
     c              max(bound(i-1),bound(i)))
            endif
            if (gap.gt.bound(i)) then
               bound(i) = bound(i) * (bound(i)/gap)
            endif
         endif
      enddo
      end

c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

      subroutine dgemm_ovwr(transa,m,n,k,alpha,A,lda,beta,B,ldb,
     c     dwork,ldwork)
c     
c     compute B <- alpha*op(A)*B + beta*B
c
      implicit none
      character*1 transa
      integer m,n,k,lda,ldb,ldwork
      double precision alpha,beta,A(lda,*),B(ldb,*),dwork(ldwork)
      integer i,j,l,blocksize

      if((m.le.0).or.(n.le.0).or.(k.le.0)) return
      if (ldwork.lt.m) stop 'Too little workspace in DGEMM_OVWR'
      if (m.gt.ldb) stop 'm>ldb in DGEMM_OVWR'
      blocksize = int(ldwork/m)
      do i=1,n-blocksize+1,blocksize
         call dgemm(transa,'N',m,blocksize,k,alpha,A,lda,
     c              B(1,i),ldb,0D0,dwork,m)
         if (beta.eq.0D0) then
            do j=0,blocksize-1
               do l=1,m
                  B(l,i+j)  = dwork(j*m+l)
               enddo
            enddo
         else
            do j=0,blocksize-1
               do l=1,m
                  B(l,i+j)  = dwork(j*m+l) + beta*B(l,i+j)
               enddo
            enddo
         endif
      enddo
      call dgemm(transa,'N',m,n-i+1,k,alpha,A,lda,
     c           B(1,i),ldb,0D0,dwork,m)
      if (beta.eq.0D0) then
         do j=0,n-i
            do l=1,m
               B(l,i+j)  = dwork(j*m+l)
            enddo
         enddo
      else
         do j=0,n-i
            do l=1,m
               B(l,i+j)  = dwork(j*m+l) + beta*B(l,i+j)
            enddo            
         enddo
      endif      
      return
      end


      subroutine dgemm_ovwr_left(transb,m,n,k,alpha,A,lda,beta,B,ldb,
     c     dwork,ldwork)
c     
c     compute  A <- alpha*A*op(B)
c
      implicit none
      character*1 transb
      integer m,n,k,lda,ldb,ldwork
      double precision alpha,beta,A(lda,*),B(ldb,*),dwork(ldwork)
      integer i,j,l,blocksize

      if((m.le.0).or.(n.le.0).or.(k.le.0)) return
      if (ldwork.lt.n) stop 'Too little workspace in DGEMM_OVWR_LEFT'
      blocksize = int(ldwork/n)
      do i=1,m-blocksize+1,blocksize
         call dgemm('n',transb,blocksize,n,k,alpha,A(i,1),lda,
     c              B,ldb,0d0,dwork,blocksize)
         do j=0,n-1
            do l=0,blocksize-1
               A(i+l,j+1) = dwork(j*blocksize+1+l)
            enddo
         enddo
      enddo
      call dgemm('n',transb,m-i+1,n,k,alpha,A(i,1),lda,
     c           B,ldb,0d0,dwork,m-i+1)
      do j=0,n-1
         do l=0,m-i
            A(i+l,j+1) = dwork(j*(m-i+1)+1+l)
         enddo
      enddo
      return
      end


c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

c
c****************************************************************************
C     This simple version of MGS is faster on Pentium machines.

      subroutine dmgs(n,k,V,ldv,vnew,index)
      implicit none
      include 'stat.h'
      integer n,k,ldv,index(*)
      double precision V(ldv,*),vnew(*)
      integer i,j,p,q,iblck
      double precision s

c     Check for quick return
      if ((k.le.0).or.(n.le.0)) return
      iblck = 1
      p = index(iblck)
      q = index(iblck+1)
      do while(p.le.k .and.p .gt.0 .and. p.le.q)
         ndot = ndot + (q-p+1)
         do i=p,q
            s = 0d0
CDIR$ LOOP COUNT(10000)
c$OMP PARALLEL DO reduction(+:s)
            do j=1,n
               s = s + V(j,i)*vnew(j)
            enddo
CDIR$ LOOP COUNT(10000)
c$OMP PARALLEL DO firstprivate(s)
            do j=1,n
               vnew(j) = vnew(j) - s*V(j,i)
            enddo
         enddo
         iblck = iblck + 2
         p = index(iblck)
         q = index(iblck+1)
      enddo
      end
c
c     (C) Rasmus Munk Larsen, Stanford University, 2000
c
      subroutine clearstat
      implicit none
      include 'stat.h'
      nopx = 0
      nreorth = 0
      ndot = 0
      nitref = 0
      nbsvd = 0
      nrestart = 0
      tmvopx = 0
      tgetu0 = 0
      tupdmu = 0
      tupdnu = 0
      tintv = 0
      tlanbpro = 0
      treorth = 0      
      treorthu = 0
      treorthv = 0
      telru = 0
      telrv = 0
      tbsvd = 0
      tnorm2 = 0
      tdot = 0
      tlansvd = 0
      nlandim = 0
      nsing = 0
      tritzvec = 0
      trestart = 0
      end

      
      subroutine printstat
      implicit none
      include 'stat.h'

      print *,'+---------------------------------------------------',
     c     '--------+'
      print *,'Dimension of Lanczos basis                  = ',nlandim
      print *,'Number of singular values requested         = ',nsing
      print *,'Number of restarts                          = ',nrestart
      print *,'Number of matrix-vector multiplications     = ',nopx
      print *,'Number of reorthogonalizations              = ',nreorth
      print *,'Number of inner products in reorth.         = ',ndot
c      print *,'Number of iterative refinement steps        = ',nitref
      print *,'Number of bidiagonal SVDs calculated        = ',nbsvd
      print *

      print *
      print *,'  Time spent doing matrix-vector multiply   = ',tmvopx
      print *,'  Time spent generating starting vectors    = ',tgetu0
      print *,'    Time spent reorthogonalizing U_{j+1}    = ',treorthu
      print *,'    Time spent reorthogonalizing V_{j}      = ',treorthv
      print *,'  Time spent reorthogonalizing              = ',treorth
      print *,'Total Time spent in LANBPRO                 = ',tlanbpro

c      print *
c      print *,'Time spent updating mu-recurrence           = ',tupdmu
c      print *,'Time spent updating nu-recurrence           = ',tupdnu
c      print *,'Time spent on local reorth. on U_{j+1}      = ',telru
c      print *,'Time spent on local reorth. on V_{j+1}      = ',telrv
c      print *,'Time spent in PDNORM2                       = ',tnorm2
c      print *,'Time spent in PDDOT                         = ',tdot
      print *
      print *,'  Time spent in LANBPRO                     = ',tlanbpro
      print *,'  Time spent computing bidiagonal SVDs      = ',tbsvd
      print *,'  Time spent doing implicit restarts        = ',trestart
      print *,'  Time spent computing Ritz vectors         = ',tritzvec
      print *
      print *,'Total Time spent in LANSVD                  = ',tlansvd
      print *,'+----------------------------------------------------',
     c     '-------+'
      end

      
      SUBROUTINE SECOND( T )
*
      REAL       T
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     July 26, 1991
*
*  Purpose
*  =======
*
*  SECOND returns the user time for a process in seconds.
*  This version gets the time from the system function ETIME.
*
*     .. Local Scalars ..
      REAL               T1
*     ..
*     .. Local Arrays ..
      REAL               TARRAY( 2 )
*     ..
*     .. External Functions ..

c      REAL               ETIME
c      EXTERNAL           ETIME

*     ..
*     .. Executable Statements ..
*

c      T1 = ETIME( TARRAY )
c      T  = TARRAY( 1 )

      call CPU_TIME(T)

      RETURN
*
*     End of SECOND
*
      END
c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

#ifdef _OPENMP
      double precision function pdznrm2(n, x, incx)
      implicit none
      integer n, incx
      complex*16 x(*)
      
      integer i
      double precision sum

      if ((n.gt.0).and.(incx.ne.0)) then    
         sum = 0d0
         if (incx.eq.1) then
c$OMP PARALLEL DO  reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + dconjg(x(i))*x(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + dconjg(x(1+(i-1)*incx))*x(1+(i-1)*incx)
            enddo
         endif
         pdznrm2 = sqrt(sum)
      else
         pdznrm2 = 0d0
      endif   
      return
      end
c
c****************************************************************************
c 
      subroutine pzscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      complex*16 alpha,x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
            do i=1,n
               x(i) = alpha*x(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(alpha, incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = alpha*x(1+(i-1)*incx)
            enddo
         endif
      endif
      return
      end


      subroutine pzdscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      double precision alpha
      complex*16 x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
            do i=1,n
               x(i) = alpha*x(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(alpha, incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = alpha*x(1+(i-1)*incx)
            enddo
         endif
      endif
      return
      end
      
c
c****************************************************************************
c 
      subroutine pzcopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 x(*),y(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then         
         if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               y(i) = x(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx, incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = x(1+(i-1)*incx)
            enddo
         endif
      endif
      return
      end

c
c****************************************************************************
c 
      subroutine pzaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 alpha,x(*),y(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then         
         if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha)  schedule(static)
           do i=1,n
               y(i) = alpha*x(i) + y(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(alpha,incx,incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c              y(1+(i-1)*incy)
            enddo
         endif
      endif
      return
      end


      subroutine pzdaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision alpha
      complex*16 x(*),y(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then         
         if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha)  schedule(static)
           do i=1,n
               y(i) = alpha*x(i) + y(i)
            enddo
         else
c$OMP PARALLEL DO firstprivate(alpha,incx,incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c              y(1+(i-1)*incy)
            enddo
         endif
      endif
      return
      end

      
c
c****************************************************************************
c 
      complex*16 function pzdotc(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 x(*),y(*)
      
      integer i
      complex*16 sum

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then    
         if (incx.eq.1 .and. incy.eq.1) then
            sum = (0d0,0d0)
c$OMP PARALLEL DO reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + dconjg(x(i)) * y(i)
            enddo
         else
            sum = (0d0,0d0)
c$OMP PARALLEL DO firstprivate(incx, incy) reduction(+:sum) 
c$OMP& schedule(static) 
            do i=1,n
               sum = sum + dconjg(x(1+(i-1)*incx)) * y(1+(i-1)*incy)
            enddo
         endif
         pzdotc = sum
      else
         pzdotc = (0d0,0d0)
      endif   
      return
      end

      complex*16 function pzdotu(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 x(*),y(*)
      
      integer i
      complex*16 sum

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then    
         if (incx.eq.1 .and. incy.eq.1) then
            sum = (0d0,0d0)
c$OMP PARALLEL DO reduction(+:sum) schedule(static)
            do i=1,n
               sum = sum + x(i) * y(i)
            enddo
         else
            sum = (0d0,0d0)
c$OMP PARALLEL DO firstprivate(incx, incy) reduction(+:sum) 
c$OMP& schedule(static) 
            do i=1,n
               sum = sum + x(1+(i-1)*incx) * y(1+(i-1)*incy)
            enddo
         endif
         pzdotu = sum
      else
         pzdotu = (0d0,0d0)
      endif   
      return
      end

c
c****************************************************************************
c 

#else

      double precision function pdznrm2(n, x, incx)
      implicit none
      integer n, incx
      complex*16 x(*)
      double precision dznrm2
      external dznrm2

      pdznrm2 = dznrm2(n, x, incx)
      end
c
c****************************************************************************
c 


      subroutine pzscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      complex*16 alpha,x(*)

      call zscal(n, alpha, x , incx)
      end

c
c****************************************************************************
c 

      subroutine pzdscal(n, alpha, x , incx)
      implicit none
      integer n, incx
      double precision alpha
      complex*16 x(*)

      call zdscal(n, alpha, x , incx)
      end

c
c****************************************************************************
c 

      subroutine pzcopy(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 x(*),y(*)

      call zcopy(n, x , incx, y, incy)
      end
c
c****************************************************************************
c 

      subroutine pzaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 alpha,x(*),y(*)
     
      call zaxpy(n, alpha, x , incx, y, incy)
      end

      subroutine pzdaxpy(n, alpha, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      double precision alpha
      complex*16 x(*),y(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0).and.(incy.ne.0)) then         
         if (incx.eq.1 .and. incy.eq.1) then
           do i=1,n
               y(i) = alpha*x(i) + y(i)
            enddo
         else
            do i=1,n
               y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c              y(1+(i-1)*incy)
            enddo
         endif
      endif
      return
      end

c
c****************************************************************************
c 
      
      complex*16 function pzdotc(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 x(*),y(*)
#ifndef __APPLE__
      complex*16 zdotc
      external zdotc

      pzdotc = zdotc(n, x , incx, y, incy)
#else
      integer i
      pzdotc = (0.0,0.0)
      if (incx.eq.1 .and. incy.eq.1) then
         do i=1,n
            pzdotc = pzdotc + conjg(x(i))*y(i)
         enddo
      else
         do i=1,n
            pzdotc = pzdotc + conjg(x(1+(i-1)*incx))*y(1+(i-1)*incy)
         enddo
      endif
#endif      
      end

      complex*16 function pzdotu(n, x , incx, y, incy)
      implicit none
      integer n, incx, incy
      complex*16 x(*),y(*)
#ifndef __APPLE__
      complex*16 zdotu
      external zdotu
      
      pzdotu = zdotu(n, x , incx, y, incy)
#else
      integer i
      pzdotu = (0.0,0.0)
      if (incx.eq.1 .and. incy.eq.1) then
         do i=1,n
            pzdotu = pzdotu + x(i)*y(i)
         enddo
      else
         do i=1,n
            pzdotu = pzdotu + x(1+(i-1)*incx)*y(1+(i-1)*incy)
         enddo
      endif
#endif
      end
#endif

c
c****************************************************************************
c 
         
      subroutine pzzero(n, x , incx)
      implicit none
      integer n, incx
      complex*16 x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO 
            do i=1,n
               x(i) = (0d0,0d0)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = (0d0,0d0)
            enddo
         endif
      endif
      return
      end

c
c****************************************************************************
c 

      subroutine pzset(n, alpha, x , incx)
      implicit none
      integer n, incx
      complex*16 alpha,x(*)
      
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
            do i=1,n
               x(i) = alpha
            enddo
         else
c$OMP PARALLEL DO firstprivate(alpha, incx) schedule(static)
            do i=1,n
               x(1+(i-1)*incx) = alpha
            enddo
         endif
      endif
      return
      end

c
c****************************************************************************
c 

  
      subroutine pzdaxpby(n,alpha,x,incx,beta,y,incy)
c
c     Y = alpha*X + beta*Y
c     

      implicit none
      double precision one,zero
      parameter(one = 1d0,zero = 0d0)
      integer n,incx,incy,i
      double precision alpha,beta
      complex*16 x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.zero .and. beta.eq.zero) then
         if (incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               y(i) = dcmplx(zero,zero)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = dcmplx(zero,zero)
            enddo
         endif
         
      else if (alpha.eq.zero .and. beta.ne.zero) then
         
         call pzdscal(n,beta,y,incy)

      else if (alpha.ne.zero .and. beta.eq.zero) then

         if (alpha.eq.one) then
            call pzcopy(n,x,incx,y,incy)
         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
               do i=1,n
                  y(i) = alpha*x(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(incx, incy, alpha) 
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)
               enddo
            endif
         endif

      else

         if (beta.eq.one) then
c DAXPY
            call pzdaxpy(n,alpha,x,incx,y,incy)
         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha,beta) 
c$OMP& schedule(static) 
               do i=1,n
                  y(i) = alpha*x(i) + beta*y(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(alpha,beta,incx,incy)  
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c                 beta*y(1+(i-1)*incy)
               enddo
            endif
         endif
      endif
      return
      end


      subroutine pzaxpby(n,alpha,x,incx,beta,y,incy)
c
c     Y = alpha*X + beta*Y
c     

      implicit none
      double precision one,zero
      parameter(one = 1d0,zero = 0d0)
      integer n,incx,incy,i
      complex*16 alpha,beta,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.dcmplx(zero,zero) .and. 
     c     beta.eq.dcmplx(zero,zero)) then
         if (incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
            do i=1,n
               y(i) = dcmplx(zero,zero)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
            do i=1,n
               y(1+(i-1)*incy) = dcmplx(zero,zero)
            enddo
         endif
         
      else if (alpha.eq.dcmplx(zero,zero) .and. 
     c        beta.ne.dcmplx(zero,zero)) then
         
         call pzscal(n,beta,y,incy)

      else if (alpha.ne.dcmplx(zero,zero) .and. 
     c        beta.eq.dcmplx(zero,zero)) then

         if (alpha.eq.one) then
            call pzcopy(n,x,incx,y,incy)
         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
               do i=1,n
                  y(i) = alpha*x(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(incx, incy, alpha) 
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)
               enddo
            endif
         endif

      else

         if (beta.eq.one) then
c DAXPY
            call pzaxpy(n,alpha,x,incx,y,incy)
         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha,beta) 
c$OMP& schedule(static) 
               do i=1,n
                  y(i) = alpha*x(i) + beta*y(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(alpha,beta,incx,incy)  
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
     c                 beta*y(1+(i-1)*incy)
               enddo
            endif
         endif
      endif
      return
      end


c
c****************************************************************************
c 

      subroutine pzaxty(n,alpha,x,incx,y,incy)
c
c     Y = alpha*X*Y
c     

      implicit none
      double precision one,zero
      parameter(one = 1d0,zero = 0d0)
      integer n,incx,incy,i
      complex*16 alpha,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.dcmplx(zero,zero)) then
         if (incy.eq.1) then
c$OMP PARALLEL DO schedule(static)
           do i=1,n
               y(i) = dcmplx(zero,zero)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
           do i=1,n
               y(1+(i-1)*incy) = dcmplx(zero,zero)
            enddo
         endif
         
      else if (alpha.ne.dcmplx(zero,zero)) then

         if (alpha.eq.one) then
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
               do i=1,n
                  y(i) = x(i)*y(i)
               enddo
            else
c$OMP PARALLEL DO shared(y,x) schedule(static)
               do i=1,n
                  y(1+(i-1)*incy) = x(1+(i-1)*incx)*y(1+(i-1)*incy)
               enddo
            endif

         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
               do i=1,n
                  y(i) = alpha*x(i)*y(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(alpha,incx,incy) 
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)*
     c                 y(1+(i-1)*incy)
               enddo
            endif
         endif
      endif
      return
      end


      subroutine pzdaxty(n,alpha,x,incx,y,incy)
c
c     Y = alpha*X*Y
c     

      implicit none
      double precision one,zero
      parameter(one = 1d0,zero = 0d0)
      integer n,incx,incy,i
      double precision alpha,x(n),y(n)

      if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
      if (alpha.eq.zero) then
         if (incy.eq.1) then
c$OMP PARALLEL DO schedule(static)
           do i=1,n
               y(i) = dcmplx(zero,zero)
            enddo
         else
c$OMP PARALLEL DO firstprivate(incy) schedule(static)
           do i=1,n
               y(1+(i-1)*incy) = dcmplx(zero,zero)
            enddo
         endif
         
      else if (alpha.ne.zero) then

         if (alpha.eq.one) then
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO  schedule(static)
               do i=1,n
                  y(i) = x(i)*y(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(incx,incy) schedule(static)
               do i=1,n
                  y(1+(i-1)*incy) = x(1+(i-1)*incx)*y(1+(i-1)*incy)
               enddo
            endif

         else
            if (incx.eq.1 .and. incy.eq.1) then
c$OMP PARALLEL DO firstprivate(alpha) schedule(static)
               do i=1,n
                  y(i) = alpha*x(i)*y(i)
               enddo
            else
c$OMP PARALLEL DO firstprivate(alpha,incx,incy) 
c$OMP& schedule(static) 
               do i=1,n
                  y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)*
     c                 y(1+(i-1)*incy)
               enddo
            endif
         endif
      endif
      return
      end


c
c****************************************************************************
c 
                 
      subroutine zzero(n, x , incx)
      implicit none
      integer n, incx
      complex*16 x(*)
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = dcmplx(0d0,0d0)
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = dcmplx(0d0,0d0)
            enddo
         endif
      endif
      return
      end

c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c


      subroutine zdgemm_ovwr_left(transb,m,n,k,A,lda,B,ldb,zwork,lzwork)
c     
c     compute  A <- A*op(B)
c
      implicit none
      character*1 transb
      integer m,n,k,lda,ldb,lzwork
      complex*16 A(lda,*),zwork(lzwork)
      double precision B(ldb,*)
      integer i,j,l,blocksize

      if((m.le.0).or.(n.le.0).or.(k.le.0)) return
      if (lzwork.lt.n) stop 'Too little workspace in ZDGEMM_OVWR_LEFT'
      blocksize = int(lzwork/n)
      i = 1
      do i=1,m-blocksize+1,blocksize
         call zdgemm(transb,blocksize,n,k, A(i,1),lda,
     c              B,ldb,zwork,blocksize)
         do j=0,n-1
            do l=0,blocksize-1
               A(i+l,j+1) = zwork(j*blocksize+1+l)
            enddo
         enddo
      enddo
      blocksize = m-i+1
      call zdgemm(transb,blocksize,n,k,A(i,1),lda,
     c           B,ldb,zwork,blocksize)
      do j=0,n-1
         do l=0,blocksize-1
            A(i+l,j+1) = zwork(j*(m-i+1)+1+l)
         enddo
      enddo
      return
      end

      subroutine  zdgemm(transb,m,n,k,A,lda,B,ldb,C,ldc)
      implicit none
      character*1 transb
      integer m,n,k,lda,ldb,ldc
      complex*16 A(lda,*), C(ldc,*)
      double precision B(ldb,*),btmp
      integer i,j,l

      do i=1,m
         do j=1,n
            C(i,j) = dcmplx(0d0,0d0)
         enddo
      enddo
      do l=1,k
         do j=1,n
            do i=1,m
               C(i,j) = C(i,j) + A(i,l)*B(j,l)
            enddo
         enddo
      enddo
      end


      subroutine  zdgemmblk(A,lda,B,ldb,C,ldc)
      implicit none
      integer blksz
      parameter (blksz=96)
      integer lda,ldb,ldc
      complex*16 A(lda,blksz), C(ldc,blksz)
      double precision B(ldb,blksz)
      integer i,j,l, i2,j2,l2

      do l=1,blksz
         do j=1,blksz
            do i=1,blksz
               C(i,j) = dcmplx(dreal(A(i,l))*B(j,l)+dreal(C(i,j)),
     c              dimag(A(i,l))*B(j,l)+dimag(C(i,j)))
            enddo
         enddo
      enddo
      end


      subroutine  zdgemm1(transb,m,n,k,A,lda,B,ldb,C,ldc)
c     
c     compute C = A * OP(B)
c     
      implicit none
      character*1 transb
      integer m,n,k,lda,ldb,ldc
      complex*16 A(lda,*), C(ldc,*)
      double precision B(ldb,*)

      integer blksz
      parameter (blksz=96)
      integer i,j,l,iblk,jblk,lblk
      complex*16 CC(blksz,blksz)
      double precision BB(blksz,blksz),btmp
      common/BBcom/BB,CC
      logical lsame
      external lsame

      if (lsame('T',transb)) then
c
c     C = A*B^T
c        
c     Comment: This manually blocked version runs at ~1.15 GFlops 
c     on a 3 GHz Pentium 4. Not very impressive.
c     Even worse on an 1.3 GHz Itanium2: 400 MFlops - awful!
c
         do lblk=1,k-blksz+1,blksz
            do jblk=1,n-blksz+1,blksz
               do l=1,blksz
                  do j=1,blksz
                     BB(j,l) = B(jblk-1+j,lblk-1+l)
                  enddo
               enddo
               do iblk=1,m-blksz+1,blksz
                  if (lblk.eq.1) then
                     do j=jblk,jblk+blksz-1
                        do i=iblk,iblk+blksz-1
                           C(i,j) = dcmplx(0d0,0d0)
                        enddo
                     enddo
                  endif
                  call zdgemmblk(A(iblk,lblk),lda,BB,blksz,
     c                 C(iblk,jblk),ldc)
               enddo
c     
c     clean up loops for i
c     
               if (lblk.eq.1) then
                  do j=jblk,jblk+blksz-1
                     do i=iblk,m
                        C(i,j) = dcmplx(0d0,0d0)
                     enddo
                  enddo
               endif
               do l=lblk,lblk+blksz-1
                  do j=jblk,jblk+blksz-1
                     btmp = B(j,l)
                     do i=iblk,m
                        C(i,j) =  dcmplx(
     c                       dreal(A(i,l))*btmp+dreal(C(i,j)),
     c                       dimag(A(i,l))*btmp+dimag(C(i,j)))
                     enddo
                  enddo
               enddo
            enddo
c     
c     clean up loops for j
c     
            if (lblk.eq.1) then
               do j=jblk,n
                  do i=1,m
                     C(i,j) = dcmplx(0d0,0d0)
                  enddo
               enddo
            endif
            do l=lblk,lblk+blksz-1
               do j=jblk,n
                  btmp = B(j,l)
                  do i=1,m
                     C(i,j) =  dcmplx(
     c                    dreal(A(i,l))*btmp+dreal(C(i,j)),
     c                    dimag(A(i,l))*btmp+dimag(C(i,j)))
                  enddo
               enddo
            enddo 
         enddo
c     
c     clean up loop for l
c     
         do l=lblk,k
            if (l.eq.1) then
               do j=1,n
                  do i=1,m
                     C(i,j) = dcmplx(0d0,0d0)
                  enddo
               enddo      
            endif
            do jblk=1,n-blksz+1,blksz
               do iblk=1,m-blksz+1,blksz
                  do j=jblk,jblk+blksz-1
                     btmp = B(j,l)
                     do i=iblk,iblk+blksz-1
                     C(i,j) =  dcmplx(
     c                    dreal(A(i,l))*btmp+dreal(C(i,j)),
     c                    dimag(A(i,l))*btmp+dimag(C(i,j)))
                     enddo
                  enddo
               enddo
               do j=jblk,jblk+blksz-1
                  btmp = B(j,l)
                  do i=iblk,m
                     C(i,j) =  dcmplx(
     c                    dreal(A(i,l))*btmp+dreal(C(i,j)),
     c                    dimag(A(i,l))*btmp+dimag(C(i,j)))
                  enddo
               enddo               
            enddo
            do j=jblk,n
               btmp = B(j,l)
               do i=1,m
                      C(i,j) =  dcmplx(
     c                    dreal(A(i,l))*btmp+dreal(C(i,j)),
     c                    dimag(A(i,l))*btmp+dimag(C(i,j)))
               enddo
            enddo
         enddo
      else

c
c     C = A*B
c        
         do iblk=1,m-blksz+1,blksz
            do jblk=1,n-blksz+1,blksz
               do j=1,blksz
                  do i=1,blksz
                     CC(i,j) = dcmplx(0d0,0d0)
                  enddo
               enddo
               do j=1,blksz
                  do l=1,k
                     do i=1,blksz
                        CC(i,j) = A(iblk-1+i,l)*B(l,jblk-1+j) + CC(i,j)
                     enddo
                  enddo
               enddo
               do j=1,blksz
                  do i=1,blksz
                     C(iblk-1+i,jblk-1+j) = CC(i,j)
                  enddo
               enddo            
            enddo
            do j=jblk,n
               do i=iblk,iblk+blksz-1
                  C(i,j) = dcmplx(0d0,0d0)
               enddo
            enddo
c     
c     clean up loop for j
c     
            do j=jblk,n
               do l=1,k
                  do i=iblk,iblk+blksz-1
                     C(i,j) = A(i,l)*B(l,j) + C(i,j)
                  enddo
               enddo
            enddo 
         enddo
         do j=1,n
            do i=iblk,m
               C(i,j) = dcmplx(0d0,0d0)
            enddo
         enddo
c     
c     clean up loop for i
c     
         do j=1,n
            do l=1,k
               do i=iblk,m
                  C(i,j) = A(i,l)*B(l,j) + C(i,j)
               enddo
            enddo
         enddo
      endif
      end

      

c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

c 
c     DGETU0: Attempt to generate a pseudo-random vector in SPAN(Op(A)) 
c     orthogonal to span(U(:,1:j)), where Op(A) = A if transa='n' and
c     Op(A) = A^H if transa='c'.
c

      subroutine zgetu0(transa, m, n, j, ntry, u0, u0norm, U, ldu,
     c     aprod, dparm, iparm, ierr, icgs, anormest, work)

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      character*1 transa
      integer m, n, j, ntry, ldu, ierr,icgs
      integer iparm(*)
      double precision u0norm,anormest
      complex*16 u0(*),U(*),work(*),dparm(*)
      external aprod

c     %------------%
c     | Parameters |
c     %------------%
      double precision kappa
      parameter(kappa = 0.717d0)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer itry,idist,iseed(4),rsize,usize,index(3)
      double precision nrm
      real t1,t2,t3

c     %--------------------%
c     | External Functions |
c     %--------------------%
      logical lsame
      double precision pdznrm2
      external pdznrm2,lsame

c-------------------- Here begins executable code ---------------------

      call second(t1)
      iseed(1) = 1
      iseed(2) = 3
      iseed(3) = 5
      iseed(4) = 7

      if (lsame(transa,'n')) then
c     %-------------------------%
c     | u0 is to be an m-vector |
c     %-------------------------%
         rsize = n
         usize = m
      else
c     %-------------------------%
c     | u0 is to be an n-vector |
c     %-------------------------%
         rsize = m
         usize = n
      endif

      idist = 2
      ierr = 0
      do itry=1,ntry
         call zlarnv(idist, iseed, rsize, work)
         nrm = pdznrm2(rsize,work,1)
         call second(t2)
         call aprod(transa,m,n,work,u0,dparm,iparm)
         call second(t3)
         tmvopx = tmvopx + (t3-t2)
         nopx = nopx+1

         u0norm = pdznrm2(usize,u0,1)
         anormest = u0norm/nrm
         
         if (j.ge.1) then
            index(1) = 1
            index(2) = j
            index(3) = j+1            
            call zreorth(usize,j,U,ldu,u0,u0norm,index,kappa,
     c           work,icgs)
         endif

         if (u0norm.gt.0) goto 9999
      enddo
      ierr = -1
 9999 call second(t2)
      tgetu0 = tgetu0 + (t2-t1)
      return
      end
      
      subroutine zlanbpro( m, n, k0, k, APROD, U, ldu, V, ldv, B, ldb,
     c     rnorm, doption, ioption, dwork, zwork, iwork, zparm, iparm, 
     c     ierr)

c     DLANBPRO: Computes K steps of the Lanczos bidiagonalization (LBD) 
c     algorithm with partial reorthogonalization (BPRO) with M-by-1 starting 
c     vector U(:,k0+1), producing a lower bidiagonal K+1-by-K matrix B_k, an 
c     N-by-K matrix V_k, an M-by-K+1 matrix U_{k+1} such that
c           A*V_k = U_{k+1}*B_k
c     Partial reorthogonalization is used to keep the columns of V_K and U_k
c     semiorthogonal to a level prescribed in DOPTION(1), i.e.
c           MAX(DIAG((EYE(K) - V_K'*V_K))) <= DOPTION(1)
c     and
c           MAX(DIAG((EYE(K) - U_K'*U_K))) <= DOPTION(1).
c
c     If K0>0 and K>K0 an existing K0-step LBD of stored in U, V and B is 
c     extended to a K-step LBD.
c
c     Parameters:
c
c     M: INTEGER. Number of rows of A.
c     N: INTEGER. Number of columns of A.
c     K0: INTEGER. The dimension of the previously computed Lanczos
c                  bidiagonalization stored in U, V, and B.
c     K: INTEGER. On entry: The desired dimension of the Lanczos 
c           bidiagonalization. On exit: the actual size of the LBD computed.
c           This can be smaller than the input value if an invariant subspace
c           is computed.
c     APROD: Subroutine defining the linear operator A. 
c            APROD should be of the form:
c
c           SUBROUTINE DAPROD(TRANSA,M,N,X,Y,DPARM,IPARM)
c           CHARACTER*1 TRANSA
c           INTEGER M,N,IPARM(*)
c           DOUBLE PRECISION X(*),Y(*),DPARM(*)
c
c           If TRANSA.EQ.'N' then the function should compute the matrix-vector
c           product Y = A * X.
c           If TRANSA.EQ.'C' then the function should compute the matrix-vector
c           product Y = A^H * X.
c           The arrays IPARM and DPARM are a means to pass user supplied
c           data to APROD without the use of common blocks.
c     U(LDU,K+1): DOUBLE PRECISION array. On return the first K+1 columns of
c               U will contain the left Lanczos vectors.
c               On entry: 
c                  If K0==0 the first column of U contains the starting
c                  vector for the Lanczos bidiagonalization. A random 
c                  starting vector is used if the first column is U is zero.
c                  If K0>0 the first K0+1 columns of U are assumed to 
c                  contain the first K0+1 left Lanczos vectors of an
c                  existing LBD.
c
c     LDU: INTEGER. Leading dimension of the array U. LDU >= M.
c     V(LDV,K): DOUBLE PRECISION array. On return the first K columns of
c               V will contain the right Lanczos vectors.
c               On entry: 
c                  If K0>0 the first K0 columns of V are assumed to 
c                  contain the first K0 right Lanczos vectors of an
c                  existing LBD.
c     LDV: INTEGER. Leading dimension of the array V. LDV >= N.
c     B(K,2): DOUBLE PRECISION array. On return the first columns of
c               B will contain the K diagonal elements of B_k, and
c               the second column of B will contain the K elements
c               of the first sub-diagonal of B_k. 
c     LDB: INTEGER. Leading dimension of the array B. LDB >= K.
c     RNORM: DOUBLE PRECISION. On entry RNORM must contain the norm of
c               the K0+1st column of U.
c               On exit RNORM contains the value of the (K+1,K) element
c               of B_k.
c     DOPTION: DOUBLE PRECISION array. 
c        doption(1) = delta. Level of orthogonality to maintain among
c          Lanczos vectors.
c        doption(2) = eta. During reorthogonalization, all vectors with
c          with components larger than eta along the latest Lanczos vector
c          will be purged.
c        doption(3) = anorm. Estimate of || A ||.
c     IOPTION: INTEGER array. 
c        ioption(1) = CGS.  If CGS.EQ.1 then reorthogonalization is done
c          using iterated classical GRAM-SCHMIDT. IF CGS.EQ.0 then 
c          reorthogonalization is done using iterated modified Gram-Schmidt.
c        ioption(2) = ELR. If ELR.EQ.1 then extended local orthogonality is
c          enforced among u_{k}, u_{k+1} and v_{k} and v_{k+1} respectively.
c     DWORK(LWORK): DOUBLE PRECISION array of dimension >= m+n+2*k+2+max(m,n).
c     ZWORK(max(M,N)): DOUBLE COMPLEX array of dimension max(M,N): Workspace.
c     IWORK(2*K+1): INTEGER ARRAY. Integer workspace.
c     DPARM: DOUBLE PRECISION array. Array used for passing data to the APROD
c         function.   
c     IPARM: INTEGER array. Array used for passing data to the APROD
c         function.   
c     IERR: INTEGER. Error status code.
c         IERR < 0  : An invariant subspace of dimension -J was found.
c         IERR == 0 : The computation succeeded.
c         IERR > 0  : The computation succeeded, but the algorithm 
c                     came close to computing an invariant subspace after 
c                     IERR steps. In a previous version this would have caused 
c                     the algorithm to switch to full reorthogonalization 
c                     after IERR steps, but that is no longer the case. 
c                     It is probably safe to ignore.
c
c     (C) Rasmus Munk Larsen, Stanford University, 1999, 2004
c


c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer m, n, k0, k, ldb, ldu, ldv, ierr
      integer ioption(*), iwork(*), iparm(*)
      double precision rnorm,B(ldb,*), doption(*), dwork(*)
      complex*16 U(ldu,*),V(ldv,*), zwork(*), zparm(*)
      external APROD

c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero, FUDGE, kappa
      parameter(one = 1.0d0, zero = 0.0d0, FUDGE = 1.01d0)
      parameter (kappa = 0.717d0)
      
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j,inu,imu,is,iidx,j0
      double precision eps,eps34,epsn2,epsn,sfmin,delta,eta,anorm,nrm
      complex*16 s
      double precision mumax,numax,alpha,beta,a1,b1,amax,anormest
      logical force_reorth,full_reorth
      real t1,t2,t3

c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external dzero,izero,pdcopy,pzaxpy
      external zgetu0,zreorth,zsafescal,zzero,pzcopy,pzdaxpy
      external dcompute_int,dupdate_nu,dupdate_mu
      
c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dlamch,pdnrm2,pddot,dlapy2,pdznrm2
      complex*16 pzdotc
      external pdnrm2,pddot,pdznrm2,pzdotc
      external dlamch,dlapy2

c-------------------- Here begins executable code ---------------------
      call second(t1)

c     %---------------------------------%
c     | Set machine dependent constants |
c     %---------------------------------%
      eps = dlamch('e')
      eps34 = eps**(3d0/4d0)
      epsn = dble(max(m,n))*eps
      epsn2 = sqrt(dble(max(m,n)))*eps
      sfmin = dlamch('s')
      
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
      if (doption(1).lt.zero) then
         delta = sqrt(eps/k)
      else
         delta = doption(1)
      endif
      if (doption(2).lt.zero) then
         eta = eps34/sqrt(dble(k))
      else
         eta = doption(2)
      endif
      if (delta.le.eta .or. delta.eq.zero) then
         full_reorth = .true.
      else
         full_reorth = .false.
      endif
      if (doption(3).gt.zero) then
         anorm = doption(3)
      else if (k0.gt.0) then
         anorm = dlapy2(B(1,1),B(1,2))
         if (anorm.le.zero) then
            ierr = -1
            goto 9999         
         endif
      else
         anorm = zero
      endif
      ierr = 0

c     %---------------------%
c     | Get starting vector |
c     %---------------------%
      if (rnorm .eq. zero) then
         call zgetu0('n',m, n, k0, 3, U(1,k0+1), rnorm, U,
     c        ldu, aprod, zparm, iparm, ierr, ioption(1), anormest,
     c        zwork)
         anorm = max(anorm,anormest)
      endif

c     %------------------------------%
c     | Set pointers into work array |
c     %------------------------------%
      imu = 1
      inu = imu+k+1
      iidx = 1
      is = 1
      call pdzero(2*k+2,dwork,1)
      call pizero(2*k+1,iwork,1)
      call pzzero(max(m,n),zwork,1)

c     %---------------------------%
c     | Prepare Lanczos iteration |
c     %---------------------------%
      if (k0.eq.0) then
         amax = zero
         alpha = zero
         beta = rnorm
         force_reorth = .false.

c     %---------------------------------------------------%
c     | Compute ||A x|| / ||x|| for a random vector x     |
c     | to make it less likely that ||A|| is grossly      |
c     | underestimated at the beginning of the iteration. |
c     %---------------------------------------------------%
         if (n.gt.m) then            
            call zgetu0('n',m,n,0,1,zwork(is),nrm,U,ldu,aprod,zparm,
     c           iparm,ierr,ioption(1),anormest,zwork(is+m))
         else            
            call zgetu0('c',m,n,0,1,zwork(is),nrm,V,ldv,aprod,zparm,
     c           iparm,ierr,ioption(1),anormest,zwork(is+n))
         endif
         anorm = max(anorm,FUDGE*anormest)
c         print *,'anorm = ',anorm
         j0 = 1
         if (beta.ne.zero) then
            call zsafescal(m,beta,U(1,1))
         endif
c         print *,'U0 = '
c         print *,(U(i,1),i=1,m)
         dwork(imu) = one
         dwork(inu) = one
      else
         force_reorth = .true.
         alpha = B(k0,1)
         beta = rnorm
         if (k0.lt.k .and. beta*delta.lt.anorm*eps) then
            full_reorth = .true.
            ierr = k0
         endif

         iwork(iidx) = 1
         iwork(iidx+1) = k0
         iwork(iidx+2) = k0+1

         call pzdscal(m,rnorm,U(1,k0+1),1)
         call second(t2)         
         call zreorth(m,k0,U,ldu,U(1,k0+1),rnorm,iwork(iidx),kappa,
     c        zwork(is),ioption(1))
         call second(t3)
         treorthu = treorthu+(t3-t2)         
         call zsafescal(m,rnorm,U(1,k0+1))
         call dset_mu(k0,dwork(imu),iwork(iidx),epsn2)         
         call dset_mu(k0,dwork(inu),iwork(iidx),epsn2)         
         beta = rnorm

c     %--------------------------------------%
c     | Estimate ||B||_2^2 as ||B^T * B||_1  |   
c     %--------------------------------------%
         B(k0,2) = beta
         amax = zero
         do j=1,k0
            amax = max(amax,B(j,1),B(j,2))
            if (j.eq.1) then
               anorm = max(anorm,FUDGE*alpha)
            else if (j.eq.2) then
               a1 = B(1,2)/amax
               a1 = FUDGE*amax*sqrt((B(1,1)/amax)**2 + a1**2 +
     c              B(2,1)/amax*a1)
               anorm = max(anorm,a1)
            else
               a1 = B(j-1,1)/amax
               b1 = B(j-1,2)/amax 
               a1 = FUDGE*amax*sqrt( a1**2 + b1**2 +
     c              a1*B(j-2,2)/amax + B(j,1)/amax*b1)
               anorm = max(anorm,a1)
            endif
         enddo
         j0 = k0+1
      endif
      numax = zero
      mumax = zero
      
c     %-------------------------------------------%
c     | Start Lanczos bidiagonalization iteration |
c     %-------------------------------------------%            

      do j=j0,k

c     %---------------------------------------------%
c     | alpha_{j} v_{j} = A'*u_{j} - beta_{j} v_{j} |
c     %---------------------------------------------%
         call second(t2)
         call aprod('c',m,n,U(1,j),V(1,j),zparm,iparm)
         call second(t3)
         tmvopx = tmvopx + (t3-t2)
         nopx = nopx+1

c         print *,'V1 = '
c         print *,(V(i,j),i=1,m)



         if (j.eq.1) then
            alpha = pdznrm2(n,V(1,j),1)
            anorm = max(anorm,FUDGE*alpha)
c            print *,'alpha = ',alpha
         else            
            call pzdaxpy(n,-beta,V(1,j-1),1,V(1,j),1)
            alpha = pdznrm2(n,V(1,j),1)

c            print *,'alpha = ',alpha

c     %------------------------------------%
c     | Extended local reorthogonalization |
c     %------------------------------------%
            call second(t2)
            if (j.gt.1 .and. ioption(2).gt.0 .and.
     c           alpha.lt.kappa*beta) then
               do i=1,ioption(2)
                  s = pzdotc(n,V(1,j-1),1,V(1,j),1)
c                  print *,'s = ',s
                  call pzaxpy(n,-s,V(1,j-1),1,V(1,j),1)
                  nrm = pdznrm2(n,V(1,j),1)
c                  print *, 'nrm = ',nrm
                  if (nrm .ge. kappa*alpha) goto 10                  
                  alpha = nrm
               enddo
 10            dwork(inu+j-2) = eps
               alpha = nrm
            endif
            call second(t3)
            telrv = telrv + (t3-t2)

            B(j,1) = alpha
            amax = max(amax,alpha)

c     %----------------------------%
c     | Update estimate of ||A||_2 |         
c     %----------------------------%
            if (j.eq.2) then
               a1 = B(1,2)/amax
               a1 = FUDGE*amax*sqrt((B(1,1)/amax)**2 + a1**2 +
     c              B(2,1)/amax*a1)
            else
               a1 = B(j-1,1)/amax
               b1 = B(j-1,2)/amax 
               a1 = FUDGE*amax*sqrt( a1**2 + b1**2 +
     c              a1*B(j-2,2)/amax + B(j,1)/amax*b1)
            endif
            anorm = max(anorm,a1)
         endif
c         print *,'j, alpha = ',j,alpha

c     %--------------------------%
c     | Update the nu recurrence |
c     %--------------------------%
         if (.not.full_reorth .and. alpha.ne.zero) then
            call dupdate_nu(numax,dwork(imu),dwork(inu),j,
     c           B(1,1),B(1,2),anorm,epsn2)
         endif

c     %------------------------------%
c     | Reorthogonalize if necessary |
c     %------------------------------%
         if ( (full_reorth .or. numax.gt.delta .or. force_reorth) 
     c        .and. alpha.ne.zero) then
c            print *,'***REORTH***'
            if (full_reorth .or. eta.eq.zero) then
               iwork(iidx) = 1
               iwork(iidx+1) = j-1
               iwork(iidx+2) = j
            else if (.not. force_reorth) then
               call dcompute_int(dwork(inu),j-1,delta,eta,iwork(iidx))
            endif
            call second(t2)
            call zreorth(n,j-1,V,ldv,V(1,j),alpha,iwork(iidx),
     c           kappa,zwork(is),ioption(1))
            call second(t3)
            treorthv = treorthv+(t3-t2)

            call dset_mu(j-1,dwork(inu),iwork(iidx),eps)
            numax = eta
            if (force_reorth) then
               force_reorth = .false.
            else
               force_reorth = .true.
            endif
         endif

c     %-----------------------------------------------%
c     | Check whether an invariant subspace was found |
c     %-----------------------------------------------% 
         if (alpha .lt. anorm*epsn .and. j.lt.k) then
c            print *,'UNVARIANT SUBSPACE'
            rnorm = alpha
            alpha = zero
c     %------------------------------------------------%
c     | Try to build an orthogonal subspace, starting  |
c     | with a random vector.                          |
c     %------------------------------------------------%
            call zgetu0('c', m, n, j-1, 3, V(1,j), alpha, V, ldv,
     c           aprod, zparm, iparm, ierr,ioption(1),anormest,
     c           zwork(is))
            if (alpha .eq. zero) then
c     %------------------------------------------------%
c     | We failed to generate a new random vector      |
c     | in span(A^T) orthogonal to span(V(:,1:j-1)).   |
c     | Most likely span(V(:,1:j-1)) is an invariant   |
c     | subspace.                                      |
c     %------------------------------------------------%
               k = j-1
               ierr = -j
               goto 9999
            else
c     %-------------------------------------------------%
c     | We have managed to generate a random vector     |
c     | in span(A^T) orthogonal to V(:,1:j-1), so we    |
c     | can continue the LBD and "deflate" the subspace |
c     | by setting alpha_{j} = 0.                       |
c     %-------------------------------------------------%
              call zsafescal(n,alpha,V(1,j))
               alpha = zero
               force_reorth = .true.
               if (delta.gt.zero) then
                  full_reorth = .false.
               endif
            endif
         else if (j.gt.1 .and. .not. full_reorth .and. j.lt.k .and.
     c           (delta*alpha .lt. anorm*eps)) then
            ierr = j
         endif            
         B(j,1) = alpha

         if (alpha.ne.zero) then
            call zsafescal(n,alpha,V(1,j))
         endif

c         print *,'V1_scaled = '
c         print *,(V(i,j),i=1,m)



c     %------------------------------------------------%
c     | beta_{j+1} u_{j+1} = A*v_{j} - alpha_{j} u_{j} |
c     %------------------------------------------------%
         call second(t2)
         call aprod('n',m,n,V(1,j),U(1,j+1),zparm,iparm)
         call second(t3)
         tmvopx = tmvopx + (t3-t2)
         nopx = nopx+1

         call pzdaxpy(m,-alpha,U(1,j),1,U(1,j+1),1)
         beta = pdznrm2(m,U(1,j+1),1)
c         print *,'beta = ',beta
         
c     %------------------------------------%
c     | Extended local reorthogonalization |
c     %------------------------------------%
         call second(t2)
         if (ioption(2).gt.0 .and. beta.lt.kappa*alpha) then
            do i=1,ioption(2)
               s = pzdotc(m,U(1,j),1,U(1,j+1),1)
c               print *,'s = ',s
               call pzaxpy(m,-s,U(1,j),1,U(1,j+1),1)
               nrm = pdznrm2(m,U(1,j+1),1)       
c               print *,'nrm = ',nrm
               if (nrm .ge. kappa*beta) goto 20
               beta = nrm
            enddo
 20         dwork(imu+j-1) = eps
            beta = nrm
         endif
         call second(t3)
         telru = telru + (t3-t2)

c         print *,'beta = ',beta

         B(j,2) = beta
         amax = max(amax,beta)
      
c     %----------------------------%
c     | Update estimate of ||A||_2 |         
c     %----------------------------%
         if (j.le.1) then
            a1 = dlapy2(B(1,1), B(1,2))
         else
            a1 = B(j,1)/amax
            a1 = amax*sqrt(a1**2 + (B(j,2)/amax)**2 +
     c           a1*B(j-1,2)/amax)
         endif      
         anorm = max(anorm,a1)
         
c     %--------------------------%
c     | Update the mu recurrence |
c     %--------------------------%
         if (.not.full_reorth .and. beta.ne.zero) then
            call dupdate_mu(mumax,dwork(imu),dwork(inu),j,B(1,1),
     c           B(1,2),anorm,epsn2)
         endif

c     %--------------------------------------%
c     | Reorthogonalize u_{j+1} if necessary |
c     %--------------------------------------%
         if ( (full_reorth .or. mumax.gt.delta .or. force_reorth) 
     c        .and. beta.ne.zero) then
            if (full_reorth .or. eta.eq.zero) then
               iwork(iidx) = 1
               iwork(iidx+1) = j
               iwork(iidx+2) = j+1
            else if (.not. force_reorth) then
               call dcompute_int(dwork(imu),j,delta,eta,iwork(iidx))
            else
               do i=1,2*j+1
                  if (iwork(iidx+i-1).eq.j) then
                     iwork(iidx+i-1) = j+1
                     goto 25
                  endif
               enddo
            endif
 
 25         call second(t2)        
            call zreorth(m,j,U,ldu,U(1,j+1),beta,iwork(iidx),
     c              kappa, zwork(is),ioption(1))
            call second(t3)
            treorthu = treorthu+(t3-t2)

            call dset_mu(j,dwork(imu),iwork(iidx),eps)
            mumax = eta
            if (force_reorth) then
               force_reorth = .false.
            else
               force_reorth = .true.
            endif
         endif

c     %-----------------------------------------------%
c     | Check whether an invariant subspace was found |
c     %-----------------------------------------------%
         if (beta .lt. anorm*epsn .and. j.lt.k) then
            rnorm = beta
            beta = zero
c     %-----------------------------------------------%
c     | Try to build an orthogonal subspace, starting |
c     | with a random vector.                         |
c     %-----------------------------------------------%
            call zgetu0('n', m, n, j, 3, U(1,j+1), beta, U, ldu, aprod,
     c           zparm, iparm, ierr,ioption(1),anormest,zwork(is))
            if (beta .eq. zero) then
c     %-----------------------------------------------%
c     | We failed to generate a new random vector     |
c     | in span(A) orthogonal to span(U(:,1:j)).      |
c     | Most likely span(U(:,1:j)) is an invariant    |
c     | subspace.                                     |
c     %-----------------------------------------------%
               k = j
               ierr = -j
               goto 9999
            else
c     %------------------------------------------------%
c     | We have managed to generate a random vector    |
c     | in span(A) orthogonal to U(:,1:j), so we can   |
c     | continue the LBD and "deflate" the subspace by |
c     | setting beta_{j+1} = 0.                        |
c     %------------------------------------------------%
               call zsafescal(n,beta,U(1,j+1))
               beta = zero
               force_reorth = .true.
               if (delta .gt. zero) then
                  full_reorth = .false.
               endif
            endif
         else if (.not.full_reorth .and. j.lt.k .and. 
     c           (delta*beta .lt. anorm*eps)) then
            ierr = j
         endif            

         B(j,2) = beta
         if (beta.ne.zero .and. beta.ne.one) then
            call zsafescal(m,beta,U(1,j+1))
         endif
         rnorm = beta
         call second(t2)
      enddo
 9999 doption(3) = anorm      
      call second(t2)
      tlanbpro = tlanbpro + (t2-t1)
      return
      end

c
c**********************************************************************
c

      subroutine dset_mu(k,mu,index,val)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer k,index(*)
      double precision mu(*), val

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j,p,q

      i=1
      do while(index(i).le.k .and. index(i).gt.0)
         p = index(i)
         q = index(i+1)
         do j=p,q
            mu(j) = val
         enddo
         i = i+2
      enddo      
      end
c
c**********************************************************************
c
      subroutine dcompute_int(mu,j,delta,eta,index)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer j,index(*)
      double precision mu(*)
      double precision delta,eta

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,k,s,ip
      real t1,t2

      call second(t1)
      if (delta.lt.eta) then
         write (6,*) 'Warning delta<eta in dcompute_int'
         return
      endif

      ip = 0
      index(1) = 0
      i=0
      do while(i.lt.j)
c     find the next mu(k), k>i where abs(mu(k)) > delta
         do k=i+1,j
            if (abs(mu(k)).gt.delta) goto 10
         enddo         
         goto 40
c     find smallest i<k such that for all j=i,..,k, m(j) >= eta
 10      do s=k,max(i,1),-1
            if (abs(mu(s)).lt.eta) goto 20
         enddo
 20      ip= ip+1
         index(ip) = s+1
         do i=s+1,j            
            if (abs(mu(i)).lt.eta) goto 30
         enddo         
 30      ip= ip+1
         index(ip) = i-1
      enddo
 40   ip = ip+1
      index(ip) = j+1
      call second(t2)
      tintv = tintv + (t2-t1)
      end
c
c**********************************************************************
c
      subroutine dupdate_mu(mumax,mu,nu,j,alpha,beta,anorm,eps1)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer j
      double precision mumax,eps1,anorm
      double precision mu(*),nu(*),alpha(*),beta(*)

c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero, FUDGE
      parameter(one = 1.0, zero = 0.0, FUDGE = 1.01)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      double precision d
      integer k
      real t1,t2
      
c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dlamch,dlapy2
      external dlamch,dlapy2

      call second(t1)
      if (j.eq.1) then
         d = eps1*(dlapy2(alpha(j), beta(j)) + alpha(1)) + eps1*anorm
         mu(1) = eps1/beta(1)
         mumax = abs(mu(1))
      else
         mu(1) = alpha(1)*nu(1)-alpha(j)*mu(1)
         d = eps1*(dlapy2(alpha(j), beta(j)) + alpha(1)) + eps1*anorm
         mu(1) = (mu(1) + dsign(d,mu(1))) / beta(j)
         mumax = abs(mu(1))
         do k=2,j-1
            mu(k) = alpha(k)*nu(k) +beta(k-1)*nu(k-1)-alpha(j)*mu(k)
            d = eps1*(dlapy2(alpha(j), beta(j)) + 
     c           dlapy2(alpha(k), beta(k-1))) + eps1*anorm
            mu(k) = (mu(k) + dsign(d,mu(k))) / beta(j)
            mumax = max(mumax,abs(mu(k)))
         enddo
         mu(j) = beta(j-1)*nu(j-1)
         d = eps1*(dlapy2(alpha(j), beta(j)) + 
     c        dlapy2(alpha(j), beta(j-1))) + eps1*anorm
         mu(j) = (mu(j) + sign(d,mu(j))) / beta(j)
         mumax = max(mumax,abs(mu(j)))
      endif
      mu(j+1) = one
      call second(t2)
      tupdmu = tupdmu + (t2-t1)
      end
c
c**********************************************************************
c

      subroutine dupdate_nu(numax,mu,nu,j,alpha,beta,anorm,eps1)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer j
      double precision numax,eps1,anorm
      double precision mu(*),nu(*),alpha(*),beta(*)

c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero, FUDGE
      parameter(one = 1.0, zero = 0.0, FUDGE = 1.01)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      double precision d
      integer k
      real t1,t2

c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dlamch,dlapy2
      external dlamch,dlapy2

      call second(t1)
      if (j.gt.1) then
         numax = zero
         do k=1,j-1
            nu(k) = beta(k)*mu(k+1) + alpha(k)*mu(k) -beta(j-1)*nu(k)
            d = eps1*(dlapy2(alpha(k),beta(k)) +
     c           dlapy2(alpha(j),beta(j-1))) + eps1*anorm
            nu(k) = (nu(k) + dsign(d,nu(k))) / alpha(j)
            numax = max(numax,abs(nu(k)))
         enddo
         nu(j) = one
      endif
      call second(t2)
      tupdnu = tupdnu + (t2-t1)
      end
c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

c
c****************************************************************************
C     This simple version of MGS is faster on Pentium machines.

      subroutine zmgs(n,k,V,ldv,vnew,index)
      implicit none
      include 'stat.h'
      integer n,k,ldv,index(*)
      complex*16 V(ldv,*),vnew(*)
      integer i,j,p,q,iblck
      complex*16 s

c     Check for quick return
      if ((k.le.0).or.(n.le.0)) return
      iblck = 1
      p = index(iblck)
      q = index(iblck+1)
      do while(p.le.k .and.p .gt.0 .and. p.le.q)
         ndot = ndot + (q-p+1)
         do i=p,q
            s = dcmplx(0d0,0d0)
CDIR$ LOOP COUNT(10000)
c$OMP PARALLEL DO reduction(+:s)
            do j=1,n
               s = s + dconjg(V(j,i))*vnew(j)
            enddo
CDIR$ LOOP COUNT(10000)
c$OMP PARALLEL DO firstprivate(s)
            do j=1,n
               vnew(j) = vnew(j) - s*V(j,i)
            enddo
         enddo
         iblck = iblck + 2
         p = index(iblck)
         q = index(iblck+1)
      enddo
      end


      subroutine pzmgs(n,k,V,ldv,vnew,index)
c     
c     Modified Gram-Schmidt orthogonalization:
c     Orthogalizes vnew against the k vectors in V by the
c     iterative process     
c     
c     FOR i= [s_1:e_1 s_2:e_2 ... s_l:e_l] DO          
c       vnew = vnew - DOT( V(:,i), vnew ) * V(:,i) 
c
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer n,k,ldv,index(*)
      complex*16 V(ldv,*),vnew(*)
c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero
      parameter(one = 1.0, zero = 0.0)
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer iblck,i,j,p,q
      complex*16 vn0,newcoef,coef

c     Check for quick return
      if ((k.le.0).or.(n.le.0)) return
      iblck = 1
      p = index(iblck)
      q = index(iblck+1)
      do while(p.le.k .and.p .gt.0 .and. p.le.q)
c     Select the next block of columns from V
         ndot = ndot + (q-p+1)
         coef = dcmplx(zero,zero)
CDIR$ LOOP COUNT(10000)
c$OMP PARALLEL DO reduction(+:coef) firstprivate(p) schedule(static) 
         do j=1,n
            coef = coef + dconjg(V(j,p))*vnew(j)
         enddo         
c   interleaved (software pipelined) loops improve performance
c   of inner loop on machines with fused multiply-add.
         do i=p+1,q
            newcoef = dcmplx(zero,zero)
CDIR$ IVDEP
CDIR$ LOOP COUNT(10000)
c$OMP PARALLEL DO reduction(+:newcoef) private(vn0) firstprivate(i,coef) 
c$OMP& schedule(static) 
            do j=1,n
               vn0 = vnew(j) - coef*V(j,i-1)
               newcoef = newcoef + dconjg(V(j,i))*vn0 
               vnew(j) = vn0
            enddo
            coef = newcoef
         enddo
CDIR$ LOOP COUNT(10000)
c$OMP PARALLEL DO firstprivate(coef,q) schedule(static)
         do j=1,n
            vnew(j) = vnew(j) - coef*V(j,q)
         enddo
         iblck = iblck + 2
         p = index(iblck)
         q = index(iblck+1)
      enddo
      end
c
c     Rasmus Munk Larsen, Stanford University, 1999, 2004.
c

      subroutine zreorth(n,k,V,ldv,vnew,normvnew,index,alpha,work,
     c     iflag)
c
c     Orthogonalize the N-vector VNEW against a subset of the columns of
c     the N-by-K matrix V(1:N,1:K) using iterated classical or modified
c     Gram-Schmidt. LDV is the leading dimension of the array containing
c     V.
c     
c     Which columns to orthogonalize against is decided by the integer
c     array INDEX = [s_1,e_1, s_2,e_2,..., s_k,e_l, s_{l+1}], which
c     selects the columns V(:,[s_1:e_1 s_2:e_2 ... s_l:e_l]). s_{l+1}
c     must be larger than k and marks the end of INDEX.
c
c     The reorthogonalization is repeated until
c
c       ||VNEW'|| > ALPHA * ||VNEW|| , 
c
c     where VNEW' is the vector obtained by orthogonalizing VNEW.  If
c     VNEW' fails to satisfy this after 4 tries, VNEW is deemed to lie
c     numerically in the span of V(:,[s_1:e_1 s_2:e_2 ... s_l:e_l]), and
c     is set to the zero vector.
c
c     On return NORMVNEW contains ||VNEW||.
c
c     WORK is a workspace array of length at least 
c
c       max e_i-s_i+1, i=1...l.
c     
c     WORK is only used if IFLAG==1.
c
c     If IFLAG==0 then iterated modified Gram-Schmidt is used.
c     If IFLAG==1 then iterated classical Gram-Schmidt is used.
c

c     References: 
c       Aake Bjorck, "Numerical Methods for Least Squares Problems",
c       SIAM, Philadelphia, 1996, pp. 68-69.
c     
c       J.~W. Daniel, W.~B. Gragg, L. Kaufman and G.~W. Stewart, 
c       ``Reorthogonalization and Stable Algorithms Updating the
c       Gram-Schmidt QR Factorization'', Math. Comp.,  30 (1976), no.
c       136, pp. 772-795.
c
c       B. N. Parlett, ``The Symmetric Eigenvalue Problem'', 
c       Prentice-Hall, Englewood Cliffs, NJ, 1980. pp. 105-109

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer n,k,ldv,iflag,index(*)
      complex*16 V(ldv,*),vnew(*),work(*)
      double precision alpha,normvnew

c     %------------%
c     | Parameters |
c     %------------%
      integer NTRY
      parameter(NTRY=5)
      
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer itry
      double precision normvnew_0
      real t2,t3
      
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external zmgs,zcgs
      
c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision pdznrm2
      external pdznrm2,pzzero

      if (k.le.0 .or. n.le.0) return

      call second(t2)

      do itry=1,NTRY
         normvnew_0 = normvnew         
         if ( iflag.eq.1 ) then
            call zcgs(n,k,V,ldv,vnew,index,work)
         else
            call zmgs(n,k,V,ldv,vnew,index)
         endif
         ndot = ndot + k
         normvnew = pdznrm2(n,vnew,1)
         if (normvnew.gt.alpha*normvnew_0) goto 9999
      enddo
      normvnew = 0d0
c     vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
      call pzzero(n,vnew,1)
 9999 call second(t3)
      treorth = treorth + (t3-t2)
      nreorth = nreorth + 1
      return
      end
c
c****************************************************************************
c

      subroutine zcgs(n,k,V,ldv,vnew,index,work)

c     Block  Gram-Schmidt orthogonalization:
c     FOR i= 1:l
c         vnew = vnew - V(:,[s_i:e_i])*(V(:,[s_i:e_i])'*vnew)
c      
c     If l=1 and s_1=1 and e_1=k then this becomes classical Gram-Schmidt.

 
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer n,k,ldv,index(*)
      complex*16 V(ldv,*),vnew(*),work(*)
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j,p,q,l, ld
      integer tid, nt, cnk,st
      complex*16 ylocal(n)
c     %--------------------%
c     | External Functions |
c     %--------------------%
      external dgemv
#ifdef _OPENMP
      integer omp_get_thread_num, omp_get_num_threads
      external omp_get_thread_num, omp_get_num_threads
#endif

c The local variable ld was introduced to circumvent an apparent
c bug in Intel's OpenMP implementation.
      ld = ldv

c$OMP PARALLEL private(i,p,q,l,tid,nt,cnk,st,j,ylocal) 
c$OMP& firstprivate(ld) shared(ndot)
#ifdef _OPENMP
      tid = omp_get_thread_num()
      nt = omp_get_num_threads()
#else
      tid = 0
      nt = 1
#endif
      cnk = n/nt
      st = tid*cnk+1

      i=1
      do while(index(i).le.k .and. index(i).gt.0)
c
c     Select the next block of columns from V
c         
         p = index(i)
         q = index(i+1)
         l = q-p+1
         if (tid.eq.0) then
            ndot = ndot + l
         endif
c     Classical Gram-Schmidt: vnew = vnew - V(:,p:q)*(V(:,p:q)'*vnew)
         if (l.gt.0) then
            if (tid.eq.nt-1) then
               cnk = n-st+1
            endif
            call zgemv('C',cnk,l,dcmplx(1d0,0d0),V(st,p),ld,vnew(st),1,
     c           dcmplx(0d0,0d0),ylocal,1)

            if (tid.eq.0) then
               do j=1,l
                  work(j) = ylocal(j)
               enddo
            endif
c$OMP BARRIER
            if (.not.tid.eq.0) then
c$OMP CRITICAL
               do j=1,l
                  work(j) = work(j) + ylocal(j)
               enddo
c$OMP END CRITICAL
            endif
c$OMP BARRIER
            call zgemv('N',cnk,l,dcmplx(-1d0,0d0),V(st,p),ld,work,1,
     c           dcmplx(0d0,0d0),ylocal,1)
            do j=1,cnk
               vnew(st+j-1) = vnew(st+j-1) + ylocal(j)
            enddo
         endif
         i = i+2
      enddo
c$OMP END PARALLEL
      end      

      subroutine zritzvec(which,jobu,jobv,m,n,k,dim,D,E,S,U,ldu,
     c     V,ldv,work,in_lwrk,zwork,lzwrk,iwork)

c
c     DRITZVEC: Compute Ritz-vectors corresponding to the K largest 
c               or smallest (depending on the value of WHICH) 
c               Ritz-values for A from the Lanczos bidiagonalization 
c               A*V_{dim} = U_{dim+1}*B_{dim}.
c
c     Parameters:
c
c     WHICH: CHARACTER*1. Decides which singular triplets to compute. 
c            If WHICH.EQ.'L' then compute triplets corresponding to the K
c            largest singular values. 
c            If WHICH.EQ.'S' then compute triplets corresponding to the K
c            smallest singular values. 
c     JOBU: CHARACTER*1. If JOBU.EQ.'Y' then compute the left singular vectors.
c           Otherwise the array U is not touched.
c     JOBV: CHARACTER*1. If JOBV.EQ.'Y' then compute the right singular 
c           vectors. Otherwise the array V is not touched.
c     M:    INTEGER. Number of rows of A.
c     N:    INTEGER. Number of columns of A.
c     K:    INTEGER. Number of desired singular triplets. K <= MIN(DIM,M,N)
c     DIM:  INTEGER. Dimension of the Krylov subspace.
c     D(DIM): DOUBLE PRECISION array. Contains the diagonal of B.
c     E(DIM): DOUBLE PRECISION array. Contains the first sub-diagonal of B.
c     S(K): DOUBLE PRECISION array. On return S contains approximation
c               to the K largest or smallest (depending on the 
c               value of WHICH) singular values of A.
c     U(LDU,DIM+1): DOUBLE COMPLEX array. On return the first K columns of U
c               will contain approximations to the left singular vectors 
c               corresponding to the K largest or smallest (depending on the 
c               value of WHICH)  singular values of A.
c               On entry the first column of U contains the starting vector
c               for the Lanczos bidiagonalization. A random starting vector
c               is used if U is zero.
c     LDU: INTEGER. Leading dimension of the array U. LDV >= M.
c     V(LDV,DIM): DOUBLE PRECISION array. On return the first K columns of V
c               will contain approximations to the right singular vectors 
c               corresponding to the K largest or smallest (depending on the 
c               value of WHICH) singular values of A.
c     LDV: INTEGER. Leading dimension of the array V. LDV >= N.
c     WORK(LWORK): DOUBLE PRECISION array. Workspace of dimension LWORK.
c     IN_LWORK: INTEGER. Dimension of WORK. 
c            LWORK should be at least 3*DIM**2 + 
c            MAX(3*DIM**2+4*DIM+4, NB*MAX(M,N)), where NB>0 is a block 
c            size, which determines how large a fraction of the work in
c            setting up the singular vectors is done using fast BLAS-3 
c            operations. NB should probably be at least 32 to achieve good
c            performance.
c     IWORK(8*DIM): INTEGER array. Integer workspace of dimension >= 8*DIM. 
c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

      
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      character*1 which, jobu,jobv
      integer m,n,k,dim,ldu,ldv,in_lwrk,iwork(*), lzwrk
      double precision D(*),E(*),S(*),work(*)
      complex*16  U(ldu,*),V(ldv,*), zwork(lzwrk)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer lwrk, mstart
      integer i,j,ip,iqt,imt,iwrk,id(1),info
      double precision c1,c2,dd(1)
      real t0,t1
      integer st,cnk,wst,wcnk, tid, nt

c     %--------------------%
c     | External Functions |
c     %--------------------%
      logical lsame
      external lsame, dbdqr, dbdsdc, dgemm_ovwr_left
#ifdef _OPENMP
      integer omp_get_thread_num, omp_get_num_threads
      external omp_get_thread_num, omp_get_num_threads
#endif
    
c-------------------- Here begins executable code ---------------------

c     %----------------------------------------------------------------------%
c     | The bidiagonal SVD is computed in a two-stage procedure:             
c     |                                                   
c     | 1. Compute a QR-factorization M^T*B = [R; 0] of the (k+1)-by-k lower 
c     |    bidiagonal matrix B.
c     | 2. Compute the SVD of the k-by-k upper bidiagonal matrix 
c     |    R = P*S*Q^T. The SVD of B is then (M*P)*S*Q^T.
c     %----------------------------------------------------------------------%


      call second(t0)

c      print *,'U = '
c      print *,((U(i,j),i=1,m),j=1,dim+1)
c      print *,'V = '
c      print *,((V(i,j),i=1,n),j=1,dim)
c      print *,'D = ',(D(i),i=1,dim)
c      print *,'E = ',(E(i),i=1,dim)

c     %-----------------------------------------%
c     | Set pointers into workspace array
c     %-----------------------------------------%
      iwrk = 1
      lwrk = in_lwrk
      imt = 1
      iqt = imt + (dim+1)**2 
      ip = iqt + dim**2
      iwrk = ip + dim**2
      lwrk = lwrk - iwrk + 1  
c     %-----------------------------------------%
c     | Compute QR-factorization
c     |   B = M * [R; 0]
c     %-----------------------------------------%
c      print *,'info = ',info

      call dbdqr((dim.eq.min(m,n)),jobu,dim,D,E,c1,c2,work(imt),dim+1)

c      print *,'D = ',(D(i),i=1,dim)
c      print *,'E = ',(E(i),i=1,dim)
c      print *,'M^T = '
c      print *,(work(imt+i-1),i=1,(dim+1)*(dim+1))


c     %-----------------------------------------%
c     | Compute SVD of R:
c     |      R = P * S * Q^T, 
c     | using the Divide-and-conquer SVD
c     %-----------------------------------------%
      call dbdsdc('u','I',dim,D,E,work(ip),dim,work(iqt),dim,dd,id,
     c     work(iwrk),iwork,info)

c      print *,'P = '
c      print *,(work(ip-1+i),i=1,dim*dim)

c      print *,'Q^T = '
c      print *,(work(iqt-1+i),i=1,dim*dim)
c      print *,'info = ',info
c      print *,(D(i),i=1,dim)

c     %-----------------------------------------%
c     | Compute left singular vectors for B
c     |    X = P^T * M^T 
c     %-----------------------------------------%
      call dgemm_ovwr('t',dim,dim+1,dim,1d0,work(ip),dim,0d0,
     c     work(imt),dim+1,work(iwrk),lwrk)
      
      
      if (lsame(jobu,'y')) then
c     %-----------------------------------------%
c     | Form left Ritz-vectors                  |
c     |   U = U * X^T                           |
c     %-----------------------------------------%
         if (lsame(which,'s')) then
            mstart = dim-k+1
         else
            mstart = 1
         endif


c         print *,'U = '
c         print *,((U(i,j),i=1,m),j=1,dim+1)
c         print *,'X = '
c         print *,(work(imt+i-1),i=1,(dim+1)*k)

c$OMP PARALLEL private(tid,nt,cnk,st,wcnk,wst) 
#ifdef _OPENMP
         tid = omp_get_thread_num()
         nt = omp_get_num_threads()
#else
         tid = 0
         nt = 1
#endif
         wcnk = lzwrk/nt
         wst = tid*wcnk+1
         cnk = m/nt
         st = tid*cnk+1
         if (tid.eq.nt-1) then
            wcnk = lzwrk-wst+1
            cnk = m-st+1
         endif
         call zdgemm_ovwr_left('t',cnk,k,dim+1,U(st,1),
     c        ldu,work(imt+mstart-1),dim+1,
     c        zwork(wst),wcnk)
c$OMP END PARALLEL 
      endif
c      print *,'U = '
c      print *,((U(i,j),i=1,m),j=1,k)



      if (lsame(jobv,'y')) then
         if (lsame(which,'s')) then
            mstart = dim-k+1
         else
            mstart = 1
         endif
c     %-----------------------------------------%
c     | Form right Ritz-vectors
c     |   V = V * Q
c     %-----------------------------------------%
c$OMP PARALLEL private(tid,nt,cnk,st,wcnk,wst) 
#ifdef _OPENMP
         tid = omp_get_thread_num()
         nt = omp_get_num_threads()
#else
         tid = 0
         nt = 1
#endif
         wcnk = lzwrk/nt
         wst = tid*wcnk+1
         cnk = n/nt
         st = tid*cnk+1
         if (tid.eq.nt-1) then
            wcnk = lzwrk-wst+1
            cnk = n-st+1
         endif
         call zdgemm_ovwr_left('t',cnk,k,dim,V(st,1),
     c        ldv,work(iqt+mstart-1),dim,
     c        zwork(wst),wcnk)
c$OMP END PARALLEL 
      endif    
c      print *,'V = '
c      print *,((V(i,j),i=1,n),j=1,k)
     
      call  second(t1)           
      tritzvec = t1-t0
      end
      
c
c     (C) Rasmus Munk Larsen, Stanford University, 2000,2004
c
      subroutine zsafescal(n,alpha,x)
c
c     Scale the vector x by 1/alpha avoiding unnecessary under- and overflow.
c

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer n
      double precision alpha
      complex*16 x(*)

c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero
      parameter(one = 1.0, zero = 0.0)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,info
      double precision sfmin

c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external pzdscal,zlascl

c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dlamch
      external dlamch

c     %-----------------%
c     | Data statements |
c     %-----------------%
      save
      data sfmin /-1d0/
      
      if (sfmin.eq.-1d0) then         
         sfmin = dlamch('s')
      endif

      if (abs(alpha).ge.sfmin) then
         call pzdscal(n,one/alpha, x, 1)
      else
         call zlascl('General',i,i,alpha,one,n,1,x,n,info)
      endif
      
      end
