#include "BMTCV.h"

#ifdef USE_OPENCV
#include "CVUtil/CVUtil.h"
#endif

#include "TRandom3.h"

namespace larlitecv {

  std::vector<larlitecv::BoundarySpacePoint> BMTCV::findBoundarySpacePoints( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v ) {
    
    std::vector<larlitecv::BoundarySpacePoint> sp_v;
    TRandom3 rand(1983);
    
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

    /*
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
    */


    std::vector< ContourList_t > plane_contours_v;
    std::vector< std::vector<ContourIndices_t> > plane_hulls_v;
    std::vector< std::vector<Defects_t> > plane_defects_v;
    for (int p=0; p<3; p++) {
      // dilate image first
      cv::Mat& cvimg = cvimg_stage1_v[p];
      cv::dilate( cvimg, cvimg, cv::Mat(), cv::Point(-1,-1), 2, 1, 1 );
      
      // find contours
      ContourList_t contour_v;
      cv::findContours( cvimg_stage1_v[p], contour_v, cv::RETR_LIST, cv::CHAIN_APPROX_SIMPLE );

      std::cout << "Plane " << p << " number of contours: " << contour_v.size() << std::endl;
      
      // for each contour, find convex hull, find defect points
      std::vector< ContourIndices_t > hull_v( contour_v.size() );
      std::vector< Defects_t > defects_v( contour_v.size() );
      for ( int idx=0; idx<(int)contour_v.size(); idx++ ) {

	Contour_t& contour = contour_v[idx];
	if ( contour.size()<10 )
	  continue;
	
	// draw contours
	cv::drawContours( cvimg_stage0_v[p], contour_v, idx, cv::Scalar( rand.Uniform(10,255),rand.Uniform(10,255),rand.Uniform(10,255),255), 1 );	
	
	// convex hull
	cv::convexHull( cv::Mat( contour ), hull_v[idx], false );

	if ( hull_v[idx].size()<=3 ) {
	  // store 
	  continue; // no defects can be found
	}

	// defects
	cv::convexityDefects( contour, hull_v[idx], defects_v[idx] );

	// plot defect point information
	for ( auto& defectpt : defects_v[idx] ) {
	  float depth = defectpt[3]/256;
	  if ( depth>3 ) {
	    int faridx = defectpt[2];
	    cv::Point ptFar( contour[faridx] );
	    cv::circle( cvimg_stage0_v[p], ptFar, 1, cv::Scalar(0,255,0,255), -1 );
	  }
	}

      }
      plane_contours_v.emplace_back( std::move(contour_v) );
      plane_hulls_v.emplace_back( std::move(hull_v) );
      plane_defects_v.emplace_back( std::move(defects_v) );
    }

    // convex hulls
    
    return sp_v;
#endif
    
  }// end of findBoundarySpacePoints

  //void breakContour( 
  
}
