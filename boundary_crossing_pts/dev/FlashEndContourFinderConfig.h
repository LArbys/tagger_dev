#ifndef __FLASHENDCOUNTOURFINDERCONFIG__
#define __FLASHENDCOUNTOURFINDERCONFIG__

#include <string>
#include <vector>

// larcv
#include "Base/PSet.h"

namespace larlitecv {

  class FlashEndContourFinderConfig {

  public:

    FlashEndContourFinderConfig() { setdefaults(); };
    virtual ~FlashEndContourFinderConfig() {};
    
    // values per plane
    std::vector<float>  pixel_value_threshold;
    std::vector<int>    clustering_minpoints;
    std::vector<double> clustering_radius;
    std::vector<int>    endpoint_time_neighborhood;
    int                 verbosity;
    float               trigger_tick;
    float               usec_per_tick;
    float               drift_distance;
    float               drift_velocity; 
    int                 search_row_radius; ///< neighborhood to gather charge for clustering
    float               flash_zrange_extension;
    int                 flash_pixelcluster_minsize; //< min. number of pixels required to make a flash-tagged cluster
    float               max_triarea;
    float               max_triarea_tight;
    float               cathode_drift_tick_correction;
    float               anode_drift_tick_correction;
    bool                make_debug_image;

    void setdefaults();

    static FlashEndContourFinderConfig FromPSet( const larcv::PSet& flashtagger_pset );
    
  };
  


}

#endif
