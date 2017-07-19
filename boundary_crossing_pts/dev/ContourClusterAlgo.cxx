#include "ContourClusterAlgo.h"
#include <assert.h>
#include <sstream>

#include "UBWireTool/UBWireTool.h"

#ifdef USE_OPENCV
#include "CVUtil/CVUtil.h"
#endif


namespace larlitecv {

  ContourCluster::ContourCluster( const std::vector< const ContourShapeMeta* >& plane_contours ) {
    resize( plane_contours.size() );
    indices.resize( plane_contours.size() );
    for (size_t p=0; p<plane_contours.size(); p++) {
      at(p).push_back( *(plane_contours[p]) );

      earlyEnd.push_back( plane_contours[p]->getFitSegmentStart() );
      earlyContours.push_back( plane_contours[p] );
      earlyDir.push_back( plane_contours[p]->getStartDir() );

      lateEnd.push_back( plane_contours[p]->getFitSegmentEnd() );
      lateContours.push_back( plane_contours[p] );
      lateDir.push_back( plane_contours[p]->getEndDir() );
    }
  }
  
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

    std::cout << __FILE__ << ":" << __LINE__ << " "
	      << " pos=" << pos3d[0] << "," << pos3d[1] << "," << pos3d[2]
	      << "  imgcoords=(" << imgcoords[0] << "," << imgcoords[1] << "," << imgcoords[2] << "," << imgcoords[3] << ")"
	      << std::endl;

    int nplanes_in_contour = 0;
    for (size_t p=0; p<img_v.size(); p++) {
      int row = imgcoords[0];
      int col = imgcoords[p+1];
      if ( clusterpix_v[p].pixel(row,col)>0 )
	nplanes_in_contour++;
    }


    if ( nplanes_in_contour==(int)img_v.size() ) {
      std::cout << __FILE__ << ":" << __LINE__ << " already in a contour. skipping." << std::endl;
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
      seed_idx[p] = containing_idx;
    }//end of loop over planes
    
    // we need at least 2 seed clusters
    if ( nplanes_found<3 ) { // 3 for now. for 2, we need an extension. but let's wait.
      std::cout << __FILE__ << ":" << __LINE__ << " Space point matches to contours in " << nplanes_found << " planes only" << std::endl;      
      return;
    }

    // make the seed cluster
    ContourCluster cluster( plane_seedcontours );
    for (int p=0; p<3; p++)
      cluster.indices[p].insert( seed_idx[p] );

    std::vector<cv::Mat> cvimg_v;
    for ( auto const& img : img_v ) {
      cv::Mat cvimg = larcv::as_gray_mat( img, 8.0, 256.0, 1.0 );
      cv::Mat cvrgb = larcv::as_mat_greyscale2bgr( img, 10.0, 100.0 );
      cvimg_v.emplace_back( std::move(cvrgb) );
    }
    
    int maxnsteps = 20;
    int iter = 0;
    bool extended = true;
    while ( extended && iter<maxnsteps ) {
    
      // ok we perform the expansion loop
      extended = extendClusterGroup( cluster, plane_contours_v, img_v, badch_v, clusterpix_v );

      for (int p=0; p<3; p++) {
	cv::Mat& cvrgb = cvimg_v[p];
	for (auto& idx : cluster.indices[p] ) {
	  auto const& contour = plane_contours_v[p].at( idx );
	  std::vector< std::vector<cv::Point> >  contour_v;
	  contour_v.push_back( contour );
	  //cv::drawContours( cvrgb, plane_contours_v[p], idx, cv::Scalar(0,0,255,255), 1 );
	  cv::drawContours( cvrgb, contour_v, 0, cv::Scalar(0,0,255,255), 1 );
	  cv::circle( cvrgb, cluster.earlyEnd[p], 3, cv::Scalar(0,255,0,255), 1 );
	  cv::circle( cvrgb, cluster.lateEnd[p], 3, cv::Scalar(0,255,0,255), 1 );	
	}
      }
      iter++;
    }//end of expansion loop

    for (int p=0; p<3; p++) {
      std::stringstream imgname;
      imgname << "contourclusterdev_plane" << p << ".png";
      cv::imwrite( imgname.str(), cvimg_v[p] );
    }
  }

  bool ContourClusterAlgo::extendClusterGroup( ContourCluster& cluster, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
					       const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v ) {
    // for each of the active ends of the cluster, we look for the best match on each plane
    // appending requires that we can build a 3D consistent line segment
    // then we update the cluster
    
    bool cluster_extended = false;
    
    class ClusterMatchScore {
    public:
      ClusterMatchScore( int iidx, double fdist, double fcosine, double ftriscore )
      {
	dist   = fdist;
	cosine = fcosine;
	triscore = ftriscore;
	idx = iidx;
      };
      virtual ~ClusterMatchScore() {};
      double dist;
      double cosine;
      double triscore;
      int idx;
      bool operator<( const ClusterMatchScore& rhs ) const {
	/*
	if ( triscore < rhs.triscore ) {
	  return true;
	}
	else if ( triscore==rhs.triscore ) {
	  if ( cosine > rhs.cosine )
	    return true;
	  else
	    return false;
	}
	*/
	if ( dist< 20.0 && rhs.dist>20.0 ) {
	  return true;
	}
	else if ( dist>20 && rhs.dist<20.0 )
	  return false;
	else if ( dist<20 && rhs.dist<20.0 ) {
	  if ( cosine < rhs.cosine )
	    return true;
	  else
	    return false;
	}
	else if ( dist>=20 && rhs.dist>=20 ) {
	  if ( dist < rhs.dist )
	    return true;
	  else
	    return false;
	}
	return false;
      };
      
    };//end of clustermatchscore

    std::vector< std::vector<ClusterMatchScore> > cluster_scores;

    // score is based on triangle formed from end-to-end and intersection of end directions
    for ( size_t p=0; p<img_v.size(); p++ ) {

      std::vector< cv::Point > earlyline;
      const cv::Point& earlypt = cluster.earlyEnd[p]; // early end on the plane
      earlyline.push_back( cluster.earlyContours[p]->getFitSegmentStart() ); // segment start
      earlyline.push_back( cluster.earlyContours[p]->getFitSegmentEnd() );   // segment end

      std::vector< cv::Point > lateline;
      const cv::Point& latept = cluster.lateEnd[p]; // late end on the plane
      lateline.push_back( cluster.lateContours[p]->getFitSegmentStart() ); // segment start
      lateline.push_back( cluster.lateContours[p]->getFitSegmentEnd() );   // segment end

      std::vector< ClusterMatchScore > early_extension_candidates_v;
      std::vector< ClusterMatchScore >  late_extension_candidates_v;
      
      for (int idx=0; idx<(int)plane_contours_v[p].size(); idx++) {

	// we already used this cluster
	if ( cluster.indices[p].find(idx)!=cluster.indices[p].end() )
	  continue;

	// get contour
	const ContourShapeMeta& contour = plane_contours_v[p][idx];

	if ( contour.size()<20 )
	  continue;

	// get closet end to early end of current cluster
	float early2start = sqrt( (contour.getFitSegmentStart().x-earlypt.x)*(contour.getFitSegmentStart().x-earlypt.x)
				  +(contour.getFitSegmentStart().y-earlypt.y)*(contour.getFitSegmentStart().y-earlypt.y) );
	float early2end = sqrt( (contour.getFitSegmentEnd().x-earlypt.x)*(contour.getFitSegmentEnd().x-earlypt.x)
				+(contour.getFitSegmentEnd().y-earlypt.y)*(contour.getFitSegmentEnd().y-earlypt.y) );
	float early2startcos = 0;
	float early2endcos   = 0;
	for (int i=0; i<2; i++) {
	  early2startcos += contour.getStartDir()[i]*cluster.earlyDir[p][i];
	  early2endcos   += contour.getEndDir()[i]*cluster.earlyDir[p][i];
	}
		
	std::vector< cv::Point > earlypair;
	if ( early2start<early2end ) {
	  earlypair.push_back( earlypt );
	  earlypair.push_back( contour.getFitSegmentStart() );
	  ClusterMatchScore matchcand( idx, early2start, early2startcos, 0.0 );
	  early_extension_candidates_v.emplace_back( matchcand );
	}
	else {
	  earlypair.push_back( earlypt );
	  earlypair.push_back( contour.getFitSegmentEnd() );
	  ClusterMatchScore matchcand( idx, early2end, early2endcos, 0.0 );
	  early_extension_candidates_v.emplace_back( matchcand );	  
	}

	// get closet end to late end of current cluster
	float late2start = sqrt( (contour.getFitSegmentStart().x-latept.x)*(contour.getFitSegmentStart().x-latept.x)
				  +(contour.getFitSegmentStart().x-latept.y)*(contour.getFitSegmentStart().y-latept.y) );
	float late2end = sqrt( (contour.getFitSegmentEnd().x-latept.x)*(contour.getFitSegmentEnd().x-latept.x)
				+(contour.getFitSegmentEnd().x-latept.y)*(contour.getFitSegmentEnd().y-latept.y) );
	float late2startcos = 0;
	float late2endcos   = 0;
	for (int i=0; i<2; i++) {
	  late2startcos += contour.getStartDir()[i]*cluster.lateDir[p][i];
	  late2endcos   += contour.getEndDir()[i]*cluster.lateDir[p][i];
	}
	
	std::vector< cv::Point > latepair;
	if ( late2start<late2end ) {
	  latepair.push_back( latept );
	  latepair.push_back( contour.getFitSegmentStart() );
	  ClusterMatchScore matchcand( idx, late2start, late2startcos, 0.0 );
	  late_extension_candidates_v.emplace_back( matchcand );	  	  
	}
	else {
	  latepair.push_back( latept );
	  latepair.push_back( contour.getFitSegmentEnd() );
	  ClusterMatchScore matchcand( idx, late2end, late2endcos, 0.0 );
	  late_extension_candidates_v.emplace_back( matchcand );	  	  
	}
	
      }//end of loop over plane conours

      std::sort( early_extension_candidates_v.begin(), early_extension_candidates_v.end() );
      std::sort( late_extension_candidates_v.begin(), late_extension_candidates_v.end() );      

      std::cout << "Plane " << p << " Early Extension matching" << std::endl;
      for (int iex=0; iex<(int)early_extension_candidates_v.size(); iex++) {
	auto& candidate = early_extension_candidates_v[iex];

	if ( candidate.cosine<-0.5 && candidate.dist<100.0 ) {
	  std::cout << " [" << candidate.idx << "] dist=" << candidate.dist << " cos=" << candidate.cosine << std::endl;

	  // append best to early match
	  auto& best_early = early_extension_candidates_v.front();
	  auto const& early_contour = plane_contours_v[p][best_early.idx];
	  cluster.earlyContours[p] = &early_contour;
	  cluster.earlyDir[p] = early_contour.getStartDir();
	  cluster.earlyEnd[p] = early_contour.getFitSegmentStart();
	  cluster.indices[p].insert( best_early.idx );
	  cluster_extended = true;
	  break;
	}
      }

      std::cout << "Plane " << p << " Late Extension matching" << std::endl;
      for (int iex=0; iex<(int)late_extension_candidates_v.size(); iex++) {
	auto& candidate = late_extension_candidates_v[iex];

	if ( candidate.cosine<-0.5 && candidate.dist<100.0 ) {
	  std::cout << " [" << candidate.idx << "] dist=" << candidate.dist << " cos=" << candidate.cosine << std::endl;
	  // append best late match
	  auto& best_late = late_extension_candidates_v.front();
	  auto const& late_contour = plane_contours_v[p][best_late.idx];
	  cluster.lateContours[p] = &late_contour;
	  cluster.lateDir[p] = late_contour.getEndDir();
	  cluster.lateEnd[p] = late_contour.getFitSegmentEnd();
	  cluster.indices[p].insert( best_late.idx );
	  cluster_extended = true;
	  break;
	}
      }
      

      
    }//end of plane

    return cluster_extended;
  }//end of extension method


  std::vector<float> ContourClusterAlgo::calculateContourIntersection( const std::vector< cv::Point >& cnt1, const std::vector< cv::Point >& cnt2 ) {
    std::vector<float> insec;
    //float Y1 = ls1[1][1] - ls1[0][1];
    //float X1 = ls1[0][0] - ls1[1][0];
    //float C1 = Y1*ls1[0][0] + X1*ls1[0][1];

    //float Y2 = ls2[1][1] - ls2[0][1];
    //float X2 = ls2[0][0] - ls2[1][0];
    //float C2 = Y2*ls2[0][0] + X2*ls2[0][1];

    float Y1 = cnt1[1].y - cnt1[0].y;
    float X1 = cnt1[0].x - cnt1[1].x;
    float C1 = Y1*cnt1[0].x + X1*cnt1[0].y;

    float Y2 = cnt2[1].y - cnt2[0].y;
    float X2 = cnt2[0].x - cnt2[1].x;
    float C2 = Y2*cnt2[0].x - X2*cnt2[0].y;

    float det = Y1*X2 - Y2*X1;
    if ( det==0 ) { 
      return insec;
    }

    insec[0] = (X2*C1 - X1*C2)/det;
    insec[1] = (Y1*C2 - Y2*C1)/det;
    
    return insec;
  }
}
