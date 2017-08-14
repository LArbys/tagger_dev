#include "ContourAStarClusterAlgo.h"

#include <sstream>
#include <stdexcept>

#include "UBWireTool/UBWireTool.h"

#include "CVUtil/CVUtil.h"

// LArOpenCV
#include "LArOpenCV/ImageCluster/AlgoClass/DefectBreaker.h"
#include "LArOpenCV/ImageCluster/AlgoData/TrackClusterCompound.h"

namespace larlitecv {

  void ContourAStarCluster::setImageMeta( const std::vector<larcv::Image2D>& img_v ) {
    m_nplanes = img_v.size();
    m_bmtcv_indices.resize(m_nplanes);
    m_plane_contours.resize(m_nplanes);

    // we make blank cv images
    for ( auto const& img : img_v ) {
      larcv::Image2D blank( img.meta() );
      blank.paint(0.0);
      cv::Mat cvimg = larcv::as_gray_mat( blank, 0, 255.0, 1.0 );
      m_cvimg_v.emplace_back( std::move(cvimg) );
      m_clusterimg_v.emplace_back( std::move(blank) );
    }
  }

  void ContourAStarCluster::addContour( int plane, const larlitecv::ContourShapeMeta* ctr, int idx ) {
    if (ctr==NULL) {
      std::stringstream ss;
      ss << __FILE__ << ":" << __LINE__ << "Adding a NULL contour shape meta";
      throw std::runtime_error( ss.str() );
    }
    m_plane_contours[plane].push_back( ctr );
    if ( idx>=0 ) {
      m_bmtcv_indices[plane].insert( idx );
    }
  }
  
  void ContourAStarCluster::updateCVImage() {
    // try to use this sparingly. it's probably slow.
    
    // we loop over all contours, fill in the image. Recontour.
    for ( int p=0; p<m_nplanes; p++) {
      std::vector< std::vector<cv::Point> > contour_list;
      int i = 0;
      for ( auto const& pctr : m_plane_contours[p] ) {
	contour_list.push_back( *pctr );
	cv::drawContours( m_cvimg_v[p], contour_list, i, cv::Scalar(255,0,0), -1 );
      }
      // we copy back into the original image (slow-slow)
      for (size_t r=0; r<m_clusterimg_v[p].meta().rows(); r++) {
	for (size_t c=0; c<m_clusterimg_v[p].meta().cols(); c++) {
	  if ( m_cvimg_v[p].at<uchar>(cv::Point(c,r))>0 )
	    m_clusterimg_v[p].set_pixel(r,c,255.0);
	}
      }
    }
  }

  void ContourAStarCluster::updateClusterContour() {
    // we cluster the pixels in cvimg_v and form a cluster
    // we use this to perform the other analyses
    m_current_contours.clear();
    m_current_contours.resize(m_nplanes);
    for ( int p=0; p<m_nplanes; p++) {
      std::vector< std::vector<cv::Point> > contour_v;
      std::vector<cv::Vec4i> hierarchy;
      cv::findContours( m_cvimg_v[p], contour_v, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, cv::Point(0,0) );
      for ( auto& ctr : contour_v ) {
	ContourShapeMeta ctrmeta( ctr, m_clusterimg_v[p].meta() );
	m_current_contours[p].emplace_back( std::move(ctrmeta) );
      }
    }
  }

  std::vector<int> ContourAStarCluster::getOverlappingRowRange() {
    int min_row = -1;
    int max_row = -1;

    for (int p=0; p<m_nplanes; p++) {
      for ( auto& ctr : m_current_contours[p] ) {
	if ( min_row > ctr.getFitSegmentStart().y || min_row<0 )
	  min_row = ctr.getFitSegmentStart().y;
	if ( max_row < ctr.getFitSegmentStart().y || max_row<0 )
	  max_row = ctr.getFitSegmentStart().y;
	if ( min_row > ctr.getFitSegmentEnd().y || min_row<0 )
	  min_row = ctr.getFitSegmentEnd().y;
	if ( max_row < ctr.getFitSegmentEnd().y || max_row<0 )
	  max_row = ctr.getFitSegmentEnd().y;
      }
    }

    std::vector<int> range(2);
    range[0] = min_row;
    range[1] = max_row;
    return range;
  }
  
  ContourAStarCluster ContourAStarClusterAlgo::makeSeedClustersFrom3DPoint( const std::vector<float>& pos3d, const std::vector<larcv::Image2D>& img_v,
									    const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
									    const float max_dist2cluster ) {
    // we establish the seed cluster


    // first check if this is in a contour
    std::vector<int> imgcoords;
    try {
      imgcoords = larcv::UBWireTool::getProjectedImagePixel( pos3d, img_v.front().meta(), img_v.size() );
    }
    catch (...) {
      std::stringstream ss;
      ss << __FILE__ << ":" << __LINE__ << " Spacepoint could not project into the image." << std::endl;
      throw std::runtime_error( ss.str() );
    }

    std::cout << __FILE__ << ":" << __LINE__ << " "
	      << " pos=" << pos3d[0] << "," << pos3d[1] << "," << pos3d[2]
	      << "  imgcoords=(" << imgcoords[0] << "," << imgcoords[1] << "," << imgcoords[2] << "," << imgcoords[3] << ")"
	      << std::endl;
    
    
    std::vector<cv::Point> imgpt;
    for (size_t p=0; p<img_v.size(); p++) {
      cv::Point pt( imgcoords[p+1], imgcoords[0] );
      imgpt.emplace_back( std::move(pt) );
    }

    std::vector< const ContourShapeMeta* > plane_seedcontours;
    int nplanes_found = 0;
    std::vector<int> seed_idx(3,-1);
    size_t nplanes = plane_contours_v.size();
    for (size_t p=0; p<nplanes; p++) {
      
      bool contains = false;
      int containing_idx = -1;
      float closest_dist = -1;
      for (int idx=0; idx<(int)plane_contours_v[p].size(); idx++) {
	// test imgpt
	//bool bboxcontains = plane_contours_v[p][idx].getBBox().contains( imgpt[p] );
	//if ( !bboxcontains )
	//continue;
	
	// more detailed test
	double dist = cv::pointPolygonTest( plane_contours_v[p][idx], imgpt[p], true );
	if ( dist>=max_dist2cluster && (closest_dist<0 || closest_dist>fabs(dist) ) ) {
	  contains       = true;
	  containing_idx = idx;
	  closest_dist   = fabs(dist);
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
    
    // Make a cluster using the seed clusters
    ContourAStarCluster seed( img_v );

    // make the seed cluster
    for (int p=0; p<3; p++) {
      if ( seed_idx[p]>=0 ) {
	seed.addContour(p, plane_seedcontours[p], seed_idx[p] );
      }
    }
    seed.updateCVImage();
    seed.updateClusterContour();
    return seed;
  }

  ContourAStarCluster ContourAStarClusterAlgo::buildClusterFromSeed( const std::vector<float>& pos3d, const std::vector<larcv::Image2D>& img_v,
								     const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
								     const float max_dist2cluster ) {
    ContourAStarCluster cluster = makeSeedClustersFrom3DPoint( pos3d, img_v, plane_contours_v, max_dist2cluster );
    std::vector<int> rowrange = cluster.getOverlappingRowRange();
    std::cout << "row range: [" << rowrange[0] << "," << rowrange[1] << "]" << std::endl;
    return cluster;
  }
  
}
