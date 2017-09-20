
CXX=g++
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWKHEADERS), $(wildcard *.h))

INCFLAGS = -g -I. -fPIC
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell larcv-config --includes)
INCFLAGS += -I$(LARLITE_BASEDIR)/../
INCFLAGS += -I$(LARLITE_BASEDIR)/UserDev
INCFLAGS += -I$(LARLITE_BASEDIR)/UserDev/BasicTool
INCFLAGS += -I$(LARLITE_BASEDIR)/UserDev/SelectionTool
INCFLAGS += -I$(LARLITECV_BASEDIR)/app
INCFLAGS += $(shell larlitecv-config --includes)
INCFLAGS += -I$(LARCV_BASEDIR)/app/ann_1.1.2/include
INCFLAGS += `root-config --cflags`
INCFLAGS += -DUSE_OPENCV=1

LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlite-config --libs) -lBasicTool_GeoAlgo
LDFLAGS += $(shell larlite-config --libs) -lSelectionTool_OpT0FinderAna -lSelectionTool_OpT0FinderApp \
        -lSelectionTool_OpT0PhotonLibrary -lSelectionTool_OpT0FinderAlgorithms -lSelectionTool_OpT0FinderBase
LDFLAGS += $(shell larlitecv-config --libs)
LDFLAGS += `root-config --ldflags --libs`

LIB=libFQDev.so
BINARIES = compare_flashmatch
BINSRCS = $(addsuffix .cxx,$(BINARIES))
#SOURCES = $(filter-out $(BINSRCS), $(wildcard *.cxx))
SOURCES = 
OBJS = $(SOURCES:.cxx=.o)

all: $(BINARIES)

%.o: %.cxx %.h
	$(CXX) $(INCFLAGS) -c -o $@ $*.cxx

libFQDev.so: $(OBJS)
	@echo "binsrcs: $(BINSRCS)"
	@echo "objects: $(OBJS)"
	$(CXX) -shared -o $@ $^ $(LDFLAGS)

compare_trackq: compare_trackq.cxx $(LIB)
	$(CXX) $(INCFLAGS) -o $@ $@.cxx $(LDFLAGS) $(LIB)

compare_flashmatch: compare_flashmatch.cxx
	$(CXX) $(INCFLAGS) -o $@ $@.cxx $(LDFLAGS)

clean:
	@rm $(BINARIES) *.o $(LIB)
