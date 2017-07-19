#ifndef __CONTOUR_CLUSTER_H__
#define __CONTOUR_CLUSTER_H__

#include <vector>

// larcv
#include "DataFormat/Image2D.h"

#include "ContourShapeMeta.h"

#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#endif


namespace larlitecv {

  class ContourCluster : std::vector< std::vector<ContourShapeMeta> > {
    friend class ContourClusterAlgo;
  public:
    ContourCluster( const std::vector< const ContourShapeMeta* >& plane_contours );
    virtual ~ContourCluster() {};

    void addEarlyContours( const std::vector< const ContourShapeMeta*>& plane_contours ) {};
    void addLateContours( const std::vector< const ContourShapeMeta*>& plane_contours ) {};
    
    std::vector< cv::Point > earlyEnd;
    std::vector< std::vector<float> > earlyDir;
    std::vector< const ContourShapeMeta* > earlyContours;
    
    std::vector< cv::Point > lateEnd;
    std::vector< std::vector<float> > lateDir;
    std::vector< const ContourShapeMeta* > lateContours;

    std::vector< std::set<int> > indices;
    
  };
  
  class ContourClusterAlgo {
  public:

    ContourClusterAlgo();
    virtual ~ContourClusterAlgo();


    void buildCluster( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v,
		       const std::vector<float>& pos3d, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v );

    bool extendClusterGroup( ContourCluster& cluster, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
			     const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v );
    
    std::vector<float> calculateContourIntersection( const std::vector< cv::Point >& cnt1, const std::vector< cv::Point >& cnt2 );
    
  };

}


#endif
