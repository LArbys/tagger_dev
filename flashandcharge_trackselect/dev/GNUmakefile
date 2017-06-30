
CXX=g++
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWKHEADERS), $(wildcard *.h))

INCFLAGS = -g -I. -fPIC
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell larcv-config --includes) 
INCFLAGS += -I$(LARLITE_BASEDIR)/UserDev
INCFLAGS += -I$(LARLITE_BASEDIR)/UserDev/BasicTool
INCFLAGS += -I$(LARLITECV_BASEDIR)/app
INCFLAGS += $(shell larlitecv-config --includes)
INCFLAGS += -I$(LARCV_BASEDIR)/app/ann_1.1.2/include
INCFLAGS += `root-config --cflags`

LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlite-config --libs) -lBasicTool_GeoAlgo
LDFLAGS += $(shell larlitecv-config --libs)
LDFLAGS += `root-config --ldflags --libs`

LIB=libFQDev.so
BINARIES = compare_trackq compare_flashmatch
BINSRCS = $(addsuffix .cxx,$(BINARIES))
#SOURCES = $(filter-out $(BINSRCS), $(wildcard *.cxx))
SOURCES = crossingPointsAnaMethods.cxx
OBJS = $(SOURCES:.cxx=.o)

all: $(LIB) $(BINARIES)

%.o: %.cxx %.h
	$(CXX) $(INCFLAGS) -c -o $@ $*.cxx

libFQDev.so: $(OBJS)
	@echo "binsrcs: $(BINSRCS)"
	@echo "objects: $(OBJS)"
	$(CXX) -shared -o $@ $^ $(LDFLAGS)

compare_trackq: compare_trackq.cxx $(LIB)
	$(CXX) $(INCFLAGS) -o $@ $@.cxx $(LDFLAGS)

compare_flashmatch: compare_flashmatch.cxx $(LIB)
	$(CXX) $(INCFLAGS) -o $@ $@.cxx $(LDFLAGS)

clean:
	@rm $(BINARIES) *.o $(LIB)
