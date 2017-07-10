
CXX=g++
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWKHEADERS), $(wildcard *.h))

INCFLAGS = -g -I. -fPIC
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell larcv-config --includes)
INCFLAGS += -I$(LARLITE_BASEDIR)/../
INCFLAGS += -I$(LARLITE_BASEDIR)/UserDev
INCFLAGS += -I$(LARLITE_BASEDIR)/UserDev/BasicTool
INCFLAGS += -I$(LARLITECV_BASEDIR)/app
INCFLAGS += $(shell larlitecv-config --includes)
INCFLAGS += -I$(LARCV_BASEDIR)/app/ann_1.1.2/include
INCFLAGS += `root-config --cflags`
INCFLAGS += -DUSE_OPENCV=1

LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlite-config --libs) -lBasicTool_GeoAlgo
LDFLAGS += $(shell larlitecv-config --libs)
LDFLAGS += `root-config --ldflags --libs`

LIB=libCrossingPtsAna.so
BINARIES = analyze_crossingpts bmt_unittest
BINSRCS = $(addsuffix .cxx,$(BINARIES))
SOURCES = BMTCV.cxx
OBJS = $(SOURCES:.cxx=.o)

all: $(LIB) $(BINARIES)

%.o: %.cxx %.h
	$(CXX) $(INCFLAGS) -c -o $@ $*.cxx

libCrossingPtsAna.so: $(OBJS)
	@echo "binsrcs: $(BINSRCS)"
	@echo "objects: $(OBJS)"
	$(CXX) -shared -o $@ $^ $(LDFLAGS)

analyze_crossingpts: analyze_crossingpts.cxx $(LIB)
	$(CXX) $(INCFLAGS) -o $@ $@.cxx $(LDFLAGS) $(LIB)

bmt_unittest: bmt_unittest.cxx $(LIB)
	$(CXX) $(INCFLAGS) -o $@ $@.cxx $(LDFLAGS) $(LIB)

clean:
	@rm $(BINARIES) *.o $(LIB)
