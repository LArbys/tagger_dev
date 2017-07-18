#ifndef __CONTOUR_CLUSTER_H__
#define __CONTOUR_CLUSTER_H__

#include <vector>

// larcv
#include "DataFormat/Image2D.h"

#include "ContourShapeMeta.h"

namespace larlitecv {

  class ContourCluster : std::vector< std::vector<ContourShapeMeta> > {
    friend class ContourClusterAlgo;
  public:
    ContourCluster() {};
    virtual ~ContourCluster() {};

    void addContours( const std::vector< const ContourShapeMeta*>& plane_contours ) {};
    
    std::vector< cv::Point > earlyEnd;
    std::vector< std::vector<float> > earlyDir;
    std::vector< cv::Point > lateEnd;
    std::vector< std::vector<float> > lateDir;
    
  };
  
  class ContourClusterAlgo {
  public:

    ContourClusterAlgo();
    virtual ~ContourClusterAlgo();


    void buildCluster( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v,
		       const std::vector<float>& pos3d, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v );
    

  };

}


#endif
