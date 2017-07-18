#include "ContourClusterAlgo.h"

#include "UBWireTool/UBWireTool.h"

namespace larlitecv {

  ContourClusterAlgo::ContourClusterAlgo() {};

  ContourClusterAlgo::~ContourClusterAlgo() {};

  void ContourClusterAlgo::buildCluster( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v,
					 const std::vector<float>& pos3d, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v ) {
    // a blank image to build merged cluster
    if ( clusterpix_v.size()==0 ) {
      for (int p=0; p<(int)img_v.size(); p++) {
	larcv::Image2D img( img_v[p].meta() );
	img.paint(0.0);
	clusterpix_v.emplace_back( std::move(img) );
      }
    }
    
    // first check if this is in a contour
    std::vector<int> imgcoords;
    try {
      imgcoords = larcv::UBWireTool::getProjectedImagePixel( pos3d, img_v.front().meta(), img_v.size() );
    }
    catch (...) {
      std::cout << __FILE__ << ":" << __LINE__ << " Spacepoint could not project into the image." << std::endl;
      return;
    }

    int nplanes_in_contour = 0;
    for (size_t p=0; p<img_v.size(); p++) {
      int row = imgcoords[0];
      int col = imgcoords[p+1];
      if ( clusterpix_v[p].pixel(row,col)>0 )
	nplanes_in_contour++;
    }

    if ( nplanes_in_contour==(int)img_v.size() ) {
      // already in a contour. skipping.
      return;
    }

    std::vector<cv::Point> imgpt;
    for (size_t p=0; p<img_v.size(); p++) {
      cv::Point pt( imgcoords[p+1], imgcoords[0] );
      imgpt.emplace_back( std::move(pt) );
    }

    // ok now we need the seed cluster
    std::vector< const ContourShapeMeta* > plane_seedcontours;
    int nplanes_found = 0;
    for (size_t p=0; p<img_v.size(); p++) {

      bool contains = false;
      int containing_idx = -1;
      for (int idx=0; idx<(int)plane_contours_v[p].size(); idx++) {
	// test imgpt
	bool bboxcontains = plane_contours_v[p][idx].getBBox().contains( imgpt[p] );
	if ( !bboxcontains )
	  continue;

	// more detailed test
	double dist = cv::pointPolygonTest( plane_contours_v[p][idx], imgpt[p], false );
	if (dist>=0 ) {
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
    }

    // we need at least 2 seed clusters
    if ( nplanes_found<2 )
      return;

    // ok we perform the expansion loop
  }

  //void ContourClusterAlgo::extendClusterGroup

}
