#include "ContourClusterAlgo.h"

#include "UBWireTool/UBWireTool.h"

namespace larlitecv {

  ContourClusterAlgo::ContourClusterAlgo() {};

  ContourClusterAlgo::~ContourClusterAlgo() {};

  void ContourClusterAlgo::buildCluster( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v,
					 const std::vector<float>& pos3d, const std::vector<ContourShapeMeta>& contours_v ) {
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

    
  }

}
