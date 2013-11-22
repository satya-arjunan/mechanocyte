SRC=\
		iolibswL\
		lupack\
		mxlibswD\
		cmlibswD\
		molibswE\
		avlibsw45\
		mechanocyte

GMVSRC=\
			 iolibswL\
			 gmvswU

#FFLAGS = -O3 -march=native -Werror -Wfatal-errors -Wall #-fprofile-use #-pg -fprofile-generate
FFLAGS = -O3 -march=native 
FEFLAGS = -Wextra 
GFF = gfortran
OBJECTS=${SRC:=.o}
GMVOBJS=${GMVSRC:=.o}
MECHANOCYTE = mechanocyte
XGMV = xgmv

all: $(MECHANOCYTE) $(XGMV)

$(MECHANOCYTE): $(OBJECTS)
		$(GFF) $(FFLAGS) $(FEFLAGS) -o $@ $(OBJECTS)

$(XGMV): $(GMVOBJS)
		$(GFF) $(FFLAGS) $(FEFLAGS) -o $@ $(GMVOBJS)

%.o: %.f
		$(GFF) $(FFLAGS) -c -o $@ $<

clean:
		rm -f $(MECHANOCYTE) *.o *.mod
		rm -f $(MECHANOCYTE).dat
	 	rm -f $(MECHANOCYTE).[0-9]*

distclean:
		rm -f $(MECHANOCYTE) *.o *.mod *.init xgmv *.mov test_out
		rm -f $(MECHANOCYTE).dat
	 	rm -f $(MECHANOCYTE).[0-9]*
		rm -f *.jpg
