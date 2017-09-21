#ifndef __CONTOUR_ASTAR_CLUSTER_H__
#define __CONTOUR_ASTAR_CLUSTER_H__

/* ---------------------------------------------------
 * ContourAStarCluster
 * 
 * This algorithm clusters contours using astar to build
 * 3D model. The model is used to extend the track into
 * other contours. Takes BMTCV, the split contours as input.
 *
 * ---------------------------------------------------*/

#include <vector>
#include <set>
#include "DataFormat/Image2D.h"

#include "TaggerContourTools/ContourShapeMeta.h"

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>

namespace larlitecv {

  class ContourAStarCluster {
    friend class ContourAStarClusterAlgo; // builds/manipulates these objects
  public:
    ContourAStarCluster() {
      makeDebugImage(false);
    };
    ContourAStarCluster( const std::vector<larcv::Image2D>& img_v, bool make_debug_img=false ) {
      makeDebugImage(make_debug_img);      
      setImageMeta(img_v);
    };
    virtual ~ContourAStarCluster();

    int numPlanes() { return m_nplanes; };
    
    //protected:
  public: // temporary for debug
    std::vector< std::set<int> > m_bmtcv_indices; //< store indices of contours we've used
    std::vector< std::vector< const ContourShapeMeta*> > m_plane_contours; //< contours we've added to the cluster
    std::vector< larcv::Image2D > m_clusterimg_v; //< contains a binary image of our cluster (should be a cv::Mat)
    std::vector< cv::Mat > m_cvimg_v;  //< stores binary image of pixels that are a part of the cluster
    std::vector< cv::Mat > m_cvpath_v; //< stores binary image of pixels that are a part of the path
    cv::Mat m_cvimg_debug;

    std::vector< std::vector< ContourShapeMeta > > m_current_contours;
    std::vector< std::vector< ContourShapeMeta > > m_path_contours;

    std::vector< std::vector<float> > m_path3d;

    int m_nplanes;
    int m_current_min;
    int m_current_max;
    bool fMakeDebugImage;

    void setImageMeta( const std::vector<larcv::Image2D>& img_v ); // set the size of the containers which have storage for each plane
    void addContour( int plane, const larlitecv::ContourShapeMeta* ctr, int idx );
    void updateCVImage();
    void updateClusterContour();
    void makeDebugImage( bool make=true ) { fMakeDebugImage = true; };
    void resetDebugImage( const std::vector<larcv::Image2D>& img_v );    
    std::vector<int> getOverlappingRowRange();
    bool getCluster3DPointAtTimeTick( const int row, const std::vector<larcv::Image2D>& img_v,
				      const std::vector<larcv::Image2D>& badch_v, bool use_badch,
				      std::vector<int>& imgcoords, std::vector<float>& pos3d );
    
  };
  
  class ContourAStarClusterAlgo {
    
  public:
    ContourAStarClusterAlgo() {};
    virtual ~ContourAStarClusterAlgo() {};

    ContourAStarCluster buildClusterFromSeed( const std::vector<float>& pos3d, const std::vector<larcv::Image2D>& img_v,
					      const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
					      const float min_dist );
    
    ContourAStarCluster makeSeedClustersFrom3DPoint( const std::vector<float>& pos3d, const std::vector<larcv::Image2D>& img_v,
						     const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
						     const float min_dist );

    ContourAStarCluster makeCluster( const std::vector<float>& pos3d, const std::vector<larcv::Image2D>& img_v,
				     const std::vector<larcv::Image2D>& badch_v,
				     const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
				     const float max_dist2cluster, const int maxloopsteps=3 );
    
    void extendSeedCluster( const std::vector<float>& pos3d, const std::vector<larcv::Image2D>& img_v,
			    const std::vector<larcv::Image2D>& badch_v,
			    const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
			    const float max_dist2cluster, const int maxloopsteps, ContourAStarCluster& cluster );

    std::vector< std::set<int> > extendClusterUsingAStarPath( ContourAStarCluster& cluster, std::vector< std::vector<float> >& path3d,
							      const std::vector<larcv::Image2D>& img_v,
							      const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
							      const float distfromend, const float distextended, const float stepsize );
    
    void fillInClusterImage( ContourAStarCluster& cluster, const std::vector< std::vector<float> >& path3d,
			     const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
			     const std::vector< std::set<int> >& cluster_indices, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
			     const float maxstepsize, const float tag_qthreshold, const int neighborhood );

    void makeDebugImage( bool make=true ) { fMakeDebugImage = make; };

  protected:
    bool fMakeDebugImage;
    
    
  };

}

#endif