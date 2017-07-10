#include "BMTCV.h"

#ifdef USE_OPENCV
#include "CVUtil/CVUtil.h"
#endif

namespace larlitecv {

  std::vector<larlitecv::BoundarySpacePoint> BMTCV::findBoundarySpacePoints( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v ) {
    std::vector<larlitecv::BoundarySpacePoint> sp_v;
    // ------------------------------------------------------------------------
    // NO OPENCV
#ifndef USE_OPENCV
    throw std::runtime_error( "In order to use BMTCV, you must compile with OpenCV" );
    return sp_v;
#else
    // ------------------------------------------------------------------------
    // HAS OPENCV

    // first convert the images into cv and binarize
    cvimg_stage0_v.clear();
    cvimg_stage1_v.clear();    
    for ( auto const& img : img_v ) {
      cv::Mat cvimg = larcv::as_gray_mat( img, 8.0, 256.0, 1.0 );
      cv::Mat cvrgb = larcv::as_mat_greyscale2bgr( img, 10.0, 100.0 );
      cv::Mat thresh( cvimg );
      cv::threshold( cvimg, thresh, 0, 255, cv::THRESH_BINARY );
      cvimg_stage0_v.emplace_back( std::move(cvrgb) );
      cvimg_stage1_v.emplace_back( std::move(thresh) );
    }

    // chop into time slices, find contours, find ends
    int sliceheight = 24;
    int nslices = 1008/sliceheight;
    //for (int islice=0; islice<nslices; islice++) {
    for (int islice=8; islice<9; islice++) {
      std::cout << "analyze slice " << islice << " for contours" << std::endl;
      for (int p=0; p<3; p++) {
	cv::Mat cvslice = cvimg_stage1_v[p]( cv::Rect(0,islice*sliceheight,3456,sliceheight) );
	// for each slice, we find contours
	std::vector< std::vector<cv::Point> > contour_v;
	cv::findContours( cvslice, contour_v, cv::RETR_LIST, cv::CHAIN_APPROX_SIMPLE, cv::Point(0,islice*sliceheight) );
	// mark stage1
	for ( int idx=0; idx<(int)contour_v.size(); idx++ ) {
	  //for ( auto const& pt: contourpt_v ) {
	  //cv::circle( cvimg_stage0_v[p], pt, 2, cv::Scalar(0,0,255,255), 1 );
	  //}
	  cv::drawContours( cvimg_stage0_v[p], contour_v, idx, cv::Scalar(0,0,255,255), 1 );
	}
      }
    }
    
    return sp_v;
#endif
    
  }// end of findBoundarySpacePoints
  
}
