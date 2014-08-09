ROOTCFLAGS  = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS    = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS   = $(shell $(ROOTSYS)/bin/root-config --glibs)

#http://stackoverflow.com/questions/21676759/after-update-to-mac-os-10-8-undefined-symbols-for-architecture-x86-64-clang
#-stdlib=libstdc++

#CXX         = g++ -m64
CXX         = g++
#CXXFLAGS    = -g -fPIC -Wno-deprecated -0 -ansi -D_GNU_SOURCE -g -02
CXXFLAGS    = -Wno-deprecated
CXXFLAGS   +="-stdlib=libstdc++"
CXXFLAGS   += $(ROOTCFLAGS)
LIBS        = $(ROOTLIBS)

NGLIBS      = $(ROOTGLIBS)
NGLIBS     += -lMinuit
GLIBS     = $(filter-out -lNew -lz , $(NGLIBS))

#INCLUDEDIR  = ./include/
#SRCDIR      = ./src/
#CXX        += -I$(INCLUDEDIR) -I.
#OUTLIB      = ./lib/


#.SUFFIXES: .cc,.C,.cxx,.hh,.h
.SUFFIXES: .cc,.cxx,.hh,.h
#.PREFIXES: ./lib/

all: runsmear

runsmear: CSA14.cc \
   mydict.o 
	$(CXX) $(CXXFLAGS) -o runsmear *.o $(GLIBS) $ $<
	touch runsmear

mydict.cxx: myLinkDef.h
	@echo "Generating dictionary ..."
	@rootcint mydict.cxx -c myLinkDef.h

mydict.o: mydict.cxx
	$(CXX) $(CXXFLAGS) -c -o mydict.o $<

# clean
clean: 
	rm -f *.o 
	rm -f runcsa
