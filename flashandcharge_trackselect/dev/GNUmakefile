
SOURCES = $(wildcard *.cxx)
FMWK_HEADERS = LinkDef.h
HEADERS = $(filter-out $(FMWKHEADERS), $(wildcard *.h))

INCFLAGS = -g -I
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell larcv-config --includes) 
INCFLAGS += -I$(LARLITE_BASEDIR)/UserDev
INCFLAGS += -I$(LARLITE_BASEDIR)/UserDev/BasicTool
INCFLAGS += -I$(LARLITECV_BASEDIR)/app
INCFLAGS += $(shell larlitecv-config --includes)
INCFLAGS += -I$(LARCV_BASEDIR)/app/ann_1.1.2/include

LDFLAGS += $(shell larcv-config --libs)
LDFLAGS += $(shell larlite-config --libs) -lBasicTool_GeoAlgo
LDFLAGS += $(shell larlitecv-config --libs)


all: compare_trackq compare_flashmatch

compare_trackq: compare_trackq.cxx
	$(CXX) $(INCFLAGS) -o $@ $@.cxx $(LDFLAGS)

compare_flashmatch: compare_flashmatch.cxx
	$(CXX) $(INCFLAGS) -o $@ $@.cxx $(LDFLAGS)


