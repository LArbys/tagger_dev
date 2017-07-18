#ifndef __CONTOUR_CLUSTER_H__
#define __CONTOUR_CLUSTER_H__

#include <vector>

// larcv
#include "DataFormat/Image2D.h"

#include "ContourShapeMeta.h"

namespace larlitecv {

  class ContourClusterAlgo {
  public:

    ContourClusterAlgo();
    virtual ~ContourClusterAlgo();


    void buildCluster( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v,
		       const std::vector<float>& pos3d, const std::vector<ContourShapeMeta>& contours_v );
    

  };

}


#endif
