#include "BMTContourFilterAlgo.h"
#include <assert.h>
#include <sstream>
#include <exception>

#include "UBWireTool/UBWireTool.h"

#ifdef USE_OPENCV
#include "CVUtil/CVUtil.h"
#endif


namespace larlitecv {

  // =======================================================================================================
  BMTContourFilterAlgo::BMTContourFilterAlgo() {};

  BMTContourFilterAlgo::~BMTContourFilterAlgo() {};

  bool BMTContourFilterAlgo::buildCluster( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v,
					   const std::vector<float>& pos3d, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
					   const float max_dist2contour ) {
    // we rely on contours produced by BMTCV::analyzeImages + BMTCV::splitcontours
    
    // a blank image to build merged cluster
    if ( clusterpix_v.size()==0 ) {
      for (int p=0; p<(int)img_v.size(); p++) {
	larcv::Image2D img( img_v[p].meta() );
	img.paint(0.0);
	clusterpix_v.emplace_back( std::move(img) );
      }
    }
    
    std::vector<int> imgcoords;
    try {
      imgcoords = larcv::UBWireTool::getProjectedImagePixel( pos3d, img_v.front().meta(), img_v.size() );
    }
    catch (...) {
      std::cout << __FILE__ << ":" << __LINE__ << " Spacepoint could not project into the image." << std::endl;
      return false;
    }
    
    // first check if boundary point is in a contour past?
    // ----------------------------------------------    
    // std::cout << __FILE__ << ":" << __LINE__ << " "
    // 	      << " pos=" << pos3d[0] << "," << pos3d[1] << "," << pos3d[2]
    // 	      << "  imgcoords=(" << imgcoords[0] << "," << imgcoords[1] << "," << imgcoords[2] << "," << imgcoords[3] << ")"
    // 	      << std::endl;

    // int nplanes_in_contour = 0;
    // for (size_t p=0; p<img_v.size(); p++) {
    //   int row = imgcoords[0];
    //   int col = imgcoords[p+1];
    //   if ( clusterpix_v[p].pixel(row,col)>0 )
    // 	nplanes_in_contour++;
    // }


    // if ( nplanes_in_contour==(int)img_v.size() ) {
    //   std::cout << __FILE__ << ":" << __LINE__ << " already in a contour. skipping." << std::endl;
    //   return;
    // }

    // Convert Image2D position into a cv::Point
    std::vector<cv::Point> imgpt;
    for (size_t p=0; p<img_v.size(); p++) {
      cv::Point pt( imgcoords[p+1], imgcoords[0] );
      imgpt.emplace_back( std::move(pt) );
    }

    // ok now we need the seed cluster
    ContourCluster foundcluster;
    bool foundcontour = isPointInContour( imgpt, img_v, plane_contours_v, max_dist2contour, foundcluster ); 

    /*
    std::vector<cv::Mat> cvimg_v;
    for ( auto const& img : img_v ) {
      cv::Mat cvimg = larcv::as_gray_mat( img, 8.0, 256.0, 1.0 );
      cv::Mat cvrgb = larcv::as_mat_greyscale2bgr( img, 10.0, 100.0 );
      cvimg_v.emplace_back( std::move(cvrgb) );
    }

    // perform the expansion loop    
    int maxnsteps = 20;
    int iter = 0;
    bool extended = true;
    while ( extended && iter<maxnsteps ) {
      std::cout << "Expansion " << iter << " ------------------------------------" << std::endl;
      bool extendbyflailing = false;
      extended = ratchetCluster( cluster, plane_contours_v, img_v, badch_v, clusterpix_v );
      if ( !extended  ) {
	extendbyflailing = true;
	extended = extendClusterGroup( cluster, plane_contours_v, img_v, badch_v, clusterpix_v );
      }

      for (int p=0; p<3; p++) {
	cv::Mat& cvrgb = cvimg_v[p];
	if ( cluster.earlyContours[p].size()>=2 ) {
	  int last     = cluster.earlyContours[p].size()-1;
	  int nextlast = cluster.earlyContours[p].size()-2;
	  cv::line( cvrgb, cluster.earlyContours[p][last]->getFitSegmentEnd(), cluster.earlyContours[p][nextlast]->getFitSegmentStart(), cv::Scalar(255,0,0,255), 1 );
	}	

	if ( cluster.earlyContours[p].size()>0 ) {
	  auto const& contour = *(cluster.earlyContours[p].back());
	  std::vector< std::vector<cv::Point> >  contour_v;
	  contour_v.push_back( contour );
	  cv::Scalar contourcolor(0,0,255,255);
	  if ( extendbyflailing )
	    contourcolor = cv::Scalar(255,0,255,255);
	  cv::drawContours( cvrgb, contour_v, 0, contourcolor, -1 );
	  cv::circle( cvrgb, cluster.lateEnd[p], 3, cv::Scalar(0,255,0,255), 1 );	
	}
	auto& earlypt = cluster.earlyEnd[p].back();
	cv::circle( cvrgb, earlypt, 3, cv::Scalar(0,255,0,255), 1 );	  	
      }
      iter++;
    }//end of expansion loop

    // for debug
    for (int p=0; p<3; p++) {
      std::stringstream imgname;
      imgname << "contourclusterdev_plane" << p << ".png";
      cv::imwrite( imgname.str(), cvimg_v[p] );
    }

    // graph attempt
    //buildContourGraph( cluster, plane_contours_v, img_v, badch_v, clusterpix_v );

    
    // fill out cluster, make solid contour
    */

    return foundcontour;
  }

  bool BMTContourFilterAlgo::isPointInContour( const std::vector<cv::Point>& imgpt, const std::vector<larcv::Image2D>& img_v,
					       const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
					       const float max_dist2contour,
					       ContourCluster& outcluster ) {

    // ok now we need the seed cluster
    std::vector< const ContourShapeMeta* > plane_seedcontours;
    int nplanes_found = 0;
    std::vector<int> seed_idx(3,-1);
    for (size_t p=0; p<img_v.size(); p++) {

      bool contains = false;
      int containing_idx = -1;
      for (int idx=0; idx<(int)plane_contours_v[p].size(); idx++) {
	// test imgpt
	bool bboxcontains = plane_contours_v[p][idx].getBBox().contains( imgpt[p] );
	if ( !bboxcontains )
	  continue;

	// more detailed test
	double dist = cv::pointPolygonTest( plane_contours_v[p][idx], imgpt[p], true );
	if (dist>=-max_dist2contour ) {
	  contains = true;
	  containing_idx = idx;
	}

	if ( contains )
	  break;
      }
	
      if ( contains ) {
	plane_seedcontours.push_back( &(plane_contours_v[p][containing_idx]) );
	nplanes_found++;
      }
      else {
	plane_seedcontours.push_back( NULL );
      }
      seed_idx[p] = containing_idx;
    }//end of loop over planes
    
    // we need at least 2 seed clusters
    if ( nplanes_found<2 ) { // 3 for now. for 2, we need an extension. but let's wait.
      std::cout << __FILE__ << ":" << __LINE__ << " Space point matches to contours in " << nplanes_found << " planes only" << std::endl;      
      return false;
    }

    // make the seed cluster
    // ContourCluster cluster( plane_seedcontours );
    // for (int p=0; p<3; p++)
    //   cluster.indices[p].insert( seed_idx[p] );

    // std::swap( outcluster, cluster );
    //outcluster.addEarlyContours( plane_seedcontours );
    
    return true;
  }    
}
