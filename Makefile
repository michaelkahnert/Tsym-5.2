# Makefile for Tsym

# FC specifies the fortran compiler
# FFLAGS is used to compile the modules in "OBJS"
# LDFLAGS are flags for linking

#**********************************************************************

FC		= gfortran 
FFLAGS		= -O3
LDFLAGS	        =

#**********************************************************************


#----------------------------------------------------------------------
SRCS	      = Tsym.F\
		gsphere.F\
		functions.F\
		smatrix.F\
                lapack.F\
                propack.F\
                OrientAvg.F\
                NumRec.F

OBJS	      = Tsym.o\
		gsphere.o\
		functions.o\
		smatrix.o\
                lapack.o\
                propack.o\
                OrientAvg.o\
                NumRec.o



#-----------------------------------------------------------------------------
TSYM=tsym.x

$(TSYM):	$(OBJS)
		$(FC) $(OBJS) $(LDFLAGS) $(FFLAGS) $(LIBNETCDF) \
		$(LINKNETCDF) -o $(TSYM)


clean:;		rm -f *.o 

veryclean:	clean
		rm -f tsym.x *~


