#ifndef __CACAEndPointFilter_h__
#define __CACAEndPointFilter_h__

/* -------------------------------------------------------------------------------
 * CACAEndPointFilter: ContourAStarClusterAlgo End Point Filter
 * 
 * Uses the ContourAStar algo to cluster a track around a proposed end point
 * We then use the information from the returned track to determine if
 * the end point is good or not
 *
 *
 * Initial author: Taritree Wongjirad (twongj01@tufts.edu)
 * History:
 *   2017/09/05 - initial writing
 * ------------------------------------------------------------------------------*/

// stdlib
#include <vector>

// larlite
#include "DataFormat/opflash.h"

// larcv
#include "DataFormat/Image2D.h"


// larlitecv
#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerContourTools/ContourShapeMeta.h"
#include "MCTruthTools/crossingPointsAnaMethods.h"

// dev
#include "ContourAStarClusterAlgo.h"

namespace larlitecv {

  class CACAEndPtFilter {
  public:
    CACAEndPtFilter() {
      fTruthInfoLoaded = false;
      fMakeDebugImage = false;
      m_verbosity = 0;
    };
    virtual ~CACAEndPtFilter() {};

    void evaluateEndPoints( const std::vector< const std::vector<larlitecv::BoundarySpacePoint>* >& sp_v, const std::vector< larlite::event_opflash* >& flash_v,
			    const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
			    const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
			    const float dist_from_wall, const float chi2_threshold, const float max_dtick,
			    std::vector< std::vector<int> >& passes_filter );

    
    bool isEndPointGood( const larlitecv::BoundarySpacePoint& pt, const larlite::opflash* associated_flash,
			 const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
			 const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,			 
			 const int endpt_type, const float dist_from_wall, const float chi2_threshold, const float max_dtick );

    larlitecv::ContourAStarClusterAlgo& getAlgo() { return m_caca; };

    larlitecv::ContourAStarCluster& getLastCluster() { return m_last_clusters.back(); };
    void clearClusters() { m_last_clusters.clear(); };

    void setTruthInformation( const std::vector<larlitecv::TruthCrossingPointAna_t>& truthinfo, const std::vector<larlitecv::RecoCrossingPointAna_t>& recoinfo );
    void setVerbosity( int v ) { m_verbosity = v; };

    void makeDebugImage( bool makeit=true ) { fMakeDebugImage = makeit; };
    int numDebugImages() { return m_cvimg_rgbdebug.size(); };
    const cv::Mat& getDebugImage( int index=0 ) { return m_cvimg_rgbdebug[index]; };

  protected:
    
    larlitecv::ContourAStarClusterAlgo m_caca;
    std::vector< larlitecv::ContourAStarCluster > m_last_clusters;
    const std::vector<larlitecv::TruthCrossingPointAna_t>* m_truthinfo_ptr_v;
    const std::vector<larlitecv::RecoCrossingPointAna_t>*  m_recoinfo_ptr_v;
    bool fTruthInfoLoaded;
    int m_verbosity;

    std::vector<cv::Mat> m_cvimg_rgbdebug;
    bool fMakeDebugImage;

  };

}


#endif
