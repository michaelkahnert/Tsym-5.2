C-----------------------------------------------------------------------
c     In this file you can add the subroutines 'zroots' and 'laguer'
c     from Numerical Recipes, provided you have a license to use that
c     software. These routines are needed for the special quadrature
c     scheme for prisms. If you do not have a license, you need to 
c     use the less efficient Gauss-Legendre double quadrature scheme.
c
c     If you replace the two dummy routines below by the proper 
c     Numerical Recipes routines, then you can use the special 
c     quadrature scheme for prisms by setting Cyl_Quad=.true. in the
c     file "params". NOTE: You also need to edit the subroutine setup. 
c     Look for the comment heading "Integration scheme", and comment
c     out the lines as indicated there.
C-----------------------------------------------------------------------

C***********************************************************************
C     finding the roots of an m-th order polynomial with Laguerre's
C     method:
C***********************************************************************
      SUBROUTINE zroots(a,m,roots,polish)
      
      write(6,*) 'FATAL ERROR: zroots currently not implemented'
      stop
      END



C***********************************************************************
c     Laguerre's method:
C***********************************************************************
      SUBROUTINE laguer(a,m,x,its)

      write(6,*) 'FATAL ERROR: laguer currently not implemented'
      stop
      END
      
      
