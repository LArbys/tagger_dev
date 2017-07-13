#ifndef __BMTCV_H__
#define __BMTCV_H__

/* -------------------------------------------------------------------------
 * BMTCV: BoundaryMuonTagger CV
 * This class does 1-plane segment shape analysis using the contour tools
 * from opencv.
 * -----------------------------------------------------------------------*/

#include <vector>

#include "DataFormat/Image2D.h"

#include "TaggerTypes/BoundarySpacePoint.h"

#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#endif



namespace larlitecv {

  typedef std::vector<cv::Point> Contour_t;
  typedef std::vector< Contour_t > ContourList_t;
  typedef std::vector< int > ContourIndices_t;
  typedef std::vector< cv::Vec4i > Defects_t;

  class BMTCV {
  public:
    BMTCV(){};
    virtual ~BMTCV() {};


    std::vector<larlitecv::BoundarySpacePoint> findBoundarySpacePoints( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v );
    //#ifdef USE_OPENCV
    //std::vector<cv::Mat>& getImages();
    //#endif    

#ifdef USE_OPENCV
    std::vector<cv::Mat> cvimg_stage0_v; // unchanged images
    std::vector<cv::Mat> cvimg_stage1_v; // contour points over time scan
    std::vector<cv::Mat> cvimg_stage2_v; // 3D-matched contour points
    std::vector<cv::Mat> cvimg_stage3_v; // 3D-matched spacepointso    
#endif
  };


}

#endif