#ifndef __FLASHENDCONTOURFINDER_H__
#define __FLASHENDCONTOURFINDER_H__

/* ================================================================================
 * FlashEndContourFinder: Find Anode/Cathode track ends using contours 
 *  and the contour cluster algo
 * --------------------------------------------------------------------------------
 * 
 * Uses contour tools from TaggerContourTools
 *
 *
 * Initial author: Taritree Wongjirad (twongj01@tufts.edu)
 * History:
 *   2017/09/29 - initial writing
 * ================================================================================*/

#include "FlashEndContourFinderConfig.h"

#include "DataFormat/opflash.h"

#include "Base/PSet.h"
#include "DataFormat/Image2D.h"

#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerContourTools/ContourShapeMeta.h"

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>

namespace larlitecv {

  class FlashEndContourFinder {

    FlashEndContourFinder() {};
    
  public:

    typedef enum { kAnode=0, kCathode, kOutOfImage } SearchMode_t;
    typedef enum { kFront=0, kBack, kNotSet } OutOfImage_t;    

    FlashEndContourFinder( SearchMode_t mode, const FlashEndContourFinderConfig& config );
    FlashEndContourFinder( SearchMode_t mode, const larcv::PSet& pset );    
    virtual ~FlashEndContourFinder() {};

    
    bool flashMatchTrackEnds( const std::vector< larlite::event_opflash* >& opflash_vv, 
			      const std::vector<larcv::Image2D>& img_v,
			      const std::vector<larcv::Image2D>& badch_v,
			      const std::vector< std::vector<larlitecv::ContourShapeMeta> >& plane_contours_v,
			      std::vector< BoundarySpacePoint >& trackendpts,
			      std::vector< int > endpoint_flash_idx);

    void makeDebugImage( bool make=true ) {fMakeDebugImage=make; };
    std::vector<cv::Mat>& getDebugImages() { return m_cvimg_debug; };
    
  protected:

    SearchMode_t fSearchMode;
    int fNPMTs;
    FlashEndContourFinderConfig fConfig;

    BoundaryEnd_t SearchModeToEndType( SearchMode_t mode );

    bool findCandidateEndsAtTick( const float tick, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
				  const std::vector< std::vector<larlitecv::ContourShapeMeta> >& plane_contours_v,
				  std::vector< BoundarySpacePoint >& trackendpts );

    
    // struct used to characterize key-points on a contour that
    // intersect a flash time
    struct CtrXing_t {
      int mincol;
      int maxcol;
      int maxqcol;
      float maxq;
      const ContourShapeMeta* pctr;
    };    
    // method generates candidate key points on a plane using contours and row of interest
    void generateKeyPointList( const int row, const std::vector<larcv::Image2D>& img_v,
			       const std::vector< std::vector<larlitecv::ContourShapeMeta> >& plane_contours_v,
			       std::vector< std::vector<FlashEndContourFinder::CtrXing_t> >& plane_xing_v );

    void findCandidateSpacePoints( const std::vector< std::vector<FlashEndContourFinder::CtrXing_t> >& plane_xing_v,
				   const int row, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
				   std::vector< std::vector<float> >& spacepoints_v );
    
    bool find3Dpoint( const int plane1, const FlashEndContourFinder::CtrXing_t& xingp1,
		      const int plane2, const FlashEndContourFinder::CtrXing_t& xingp2,
		      const int row, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
		      std::vector<float>& pos3d );

    bool fMakeDebugImage;

    std::vector< cv::Mat > m_cvimg_debug;
    
  };
  
}

#endif


