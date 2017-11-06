#ifndef CONTOUR_GAP_FILL_H
#define CONTOUR_GAP_FILL_H

#include <vector>

#include "DataFormat/Image2D.h"
#include "DataFormat/ImageMeta.h"
#include "TaggerContourTools/BMTCV.h"

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>


namespace larlitecv {

  class ContourGapFill {
  public:

    ContourGapFill();
    virtual ~ContourGapFill();


    void makeGapFilledImages( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
			      const std::vector< std::vector<larlitecv::ContourShapeMeta> >& plane_contours_v,
			      std::vector<larcv::Image2D>& gapfilled_v );

    //void mergeWireAndGap( const std::vector<larcv::Image2D>& img_v, std::vector<larcv::Image2D>& gapfilled_v );

    void makeDebugImage( bool debug );

    std::vector<cv::Mat>& getDebugImages() { return m_cvimg_debug_v; };

    struct BadChSpan {
      int start;
      int end;
      int width;
      int ngood;
      int planeid;
      std::set<int> leftctridx;  // index of contours that touch left side of span
      std::set<int> rightctridx; // index of contours that touch right side of span
    };
    
  protected:

    std::vector< BadChSpan > findBadChSpans( const larcv::Image2D& badch, int goodchwidth=2 );
    
    bool fMakeDebugImage;

    std::vector<larcv::ImageMeta> m_meta_v;
    std::vector<cv::Mat> m_cvimg_debug_v;
    std::vector< ContourList_t > m_plane_contour_v;

    void createDebugImage( const std::vector<larcv::Image2D>& img_v );
    void associateContoursToSpans( const std::vector<larlitecv::ContourShapeMeta>& contour_v,
				   const larcv::ImageMeta& meta,
				   std::vector<BadChSpan>& span_v,
				   const int colwidth );
    

  };



}

#endif
