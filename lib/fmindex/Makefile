
CCP=g++
#CFLAGS=-g -W -Wall -O0 -DDEBUG -fPIC
CFLAGS=-W -Wall -O3 -fPIC

INCCDS=./libcds/includes/
INCDIVSUF=./libdivsufsort/include/

PYTHON_INCLUDE= -I/usr/include/python2.6 -I/usr/lib/python2.6/config
PYTHON_DLNK= -lpython2.6

all: ./libcds/lib/libcds.a libdivsufsort/lib/libdivsufsort.a fmindex_wrap.cxx _fmindex.so fmbuild

./libcds/lib/libcds.a: 
	$(MAKE) -C libcds
	
libdivsufsort/lib/libdivsufsort.a: 
	cd libdivsufsort;cmake .;$(MAKE);cd ..


_fmindex.so: ./libcds/lib/libcds.a libdivsufsort/lib/libdivsufsort.a util.o FM.o fmindex_wrap.o
	$(CCP) -shared -o _fmindex.so fmindex_wrap.o FM.o util.o libdivsufsort/lib/libdivsufsort.a ./libcds/lib/libcds.a

fmindex_wrap.cxx:
	swig -c++ -python fmindex.i

fmbuild: fmbuild.cpp ./libdivsufsort/lib/libdivsufsort.a ./libcds/lib/libcds.a FM.o util.o
	$(CCP) -I $(INCCDS) -I $(INCDIVSUF) $(CFLAGS) -o fmbuild fmbuild.cpp FM.o util.o ./libdivsufsort/lib/libdivsufsort.a ./libcds/lib/libcds.a

# pattern rule for all objects files
%.o: %.c *.h
	$(CCP) -I $(INCCDS) -I $(INCDIVSUF) -c $(CFLAGS) $< -o $@

%.o: %.cpp *.h
	$(CCP) -I $(INCCDS) -I $(INCDIVSUF) -c $(CFLAGS) $< -o $@

%.o: %.cxx
	$(CCP) -I $(INCCDS) -I $(INCDIVSUF) $(PYTHON_INCLUDE) -c $(CFLAGS) $< -o $@

clean:
	rm -f fmbuild *.a *~ *.o *.cxx *.py *.so;
	
cleanall:
	rm -f *~ *.o ; cd libdivsufsort; make -f Makefile clean ; cd .. ; cd libcds ; make -f Makefile clean ; cd ..
