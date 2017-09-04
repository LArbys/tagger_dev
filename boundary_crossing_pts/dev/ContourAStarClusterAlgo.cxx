#include "ContourAStarClusterAlgo.h"

#include <sstream>
#include <stdexcept>

// larlite
#include "LArUtil/LArProperties.h"
#include "LArUtil/Geometry.h"

#include "UBWireTool/UBWireTool.h"

#include "CVUtil/CVUtil.h"
#include "DataFormat/Pixel2D.h"

#include "ThruMu/AStar3DAlgoConfig.h"
#include "ThruMu/AStar3DAlgo.h"

// LArOpenCV
#include "LArOpenCV/ImageCluster/AlgoClass/DefectBreaker.h"
#include "LArOpenCV/ImageCluster/AlgoData/TrackClusterCompound.h"



namespace larlitecv {

  void ContourAStarCluster::setImageMeta( const std::vector<larcv::Image2D>& img_v ) {
    m_nplanes = img_v.size();
    m_bmtcv_indices.resize(m_nplanes);
    m_plane_contours.resize(m_nplanes);

    // we make blank cv images
    m_cvimg_v.clear();
    m_clusterimg_v.clear();
    for ( auto const& img : img_v ) {
      larcv::Image2D blank( img.meta() );
      blank.paint(0.0);
      cv::Mat cvimg = larcv::as_gray_mat( blank, 0, 255.0, 1.0 );
      m_cvimg_v.emplace_back( std::move(cvimg) );
      m_clusterimg_v.emplace_back( std::move(blank) );
    }

    m_cvimg_debug = larcv::as_mat_greyscale2bgr( img_v.front(), 0, 255.0 );
    for (int p=1; p<m_nplanes; p++) {
      for (size_t r=0; r<img_v[p].meta().rows(); r++) {
	for (size_t c=0; c<img_v[p].meta().cols(); c++) {
	  if ( img_v[p].pixel(r,c)>5.0 ) {
	    for (int i=0; i<3; i++)
	      m_cvimg_debug.at<cv::Vec3b>(cv::Point(c,r))[i] = img_v[p].pixel(r,c);
	  }
	}
      }
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

	// debug image
	if ( p==0 )
	  cv::drawContours( m_cvimg_debug, contour_list, i, cv::Scalar(255,0,0), 1 );
	else if (p==1)
	  cv::drawContours( m_cvimg_debug, contour_list, i, cv::Scalar(0,255,0), 1 );
	else if ( p==2)
	  cv::drawContours( m_cvimg_debug, contour_list, i, cv::Scalar(0,0,255), 1 );
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
	
	if ( min_row < ctr.getFitSegmentStart().y || min_row<0 )
	  min_row = ctr.getFitSegmentStart().y;

	if ( max_row > ctr.getFitSegmentEnd().y || max_row<0 )
	  max_row = ctr.getFitSegmentEnd().y;
      }
    }

    std::vector<int> range(2);
    range[0] = min_row;
    range[1] = max_row;
    return range;
  }

  void ContourAStarCluster::getCluster3DPointAtTimeTick( const int row, const std::vector<larcv::Image2D>& img_v,
							 const std::vector<larcv::Image2D>& badch_v, bool use_badch,
							 std::vector<int>& imgcoords, std::vector<float>& pos3d ) {
    // we get 2D points from each plane. Then we infer 3D point
    struct ctr_pt_t {
      int plane;
      int minc;
      int maxc;
      int maxq;
      float maxqc;
    };
    
    std::vector<ctr_pt_t> contour_points;
    
    for (int p=0; p<m_nplanes; p++) {
      const std::vector<ContourShapeMeta>& ctr_v = m_current_contours[p];
      const larcv::Image2D& img                  = img_v[p];
      int nfound = 0;
      for ( auto const& ctr : ctr_v ) {

	// scan across the cols at a certain time and get the min,max and max-q points inside the cluster
	int minc = -1;
	int maxc = -1;
	int maxqc = -1;
	float maxq = -1;
	bool incontour = false;
	for (int c=0; c<(int)img.meta().cols(); c++) {
	  cv::Point testpt( c, row );
	  double dist = cv::pointPolygonTest( ctr, testpt, false );
	  //std::cout << "point (" << c << "," << row << ") dist=" << dist << " incontour=" << incontour << std::endl;
	  if ( dist<0 ) {
	    if ( incontour ) {
	      // close out a contour crossing
	      ctr_pt_t ctr_xing;
	      ctr_xing.plane = p;
	      ctr_xing.minc  = minc;
	      ctr_xing.maxc  = maxc;
	      ctr_xing.maxq  = maxq;
	      ctr_xing.maxqc = maxqc;
	      contour_points.emplace_back( std::move(ctr_xing) );
	      incontour = false;
	      minc = -1;
	      maxc = -1;
	      maxqc = -1;
	      maxq = -1;
	      nfound++;
	    }
	    continue;	    
	  }

	  incontour = true;
	  
	  if ( minc<0 || c<minc ) {
	    minc = c;
	  }
	  if ( maxc<0 || c>maxc )
	    maxc = c;
	  if ( maxq<0 || img.pixel(row,c)>maxq ) {
	    maxq  = img.pixel(row,c);
	    maxqc = c;
	  }
	}//end of col loops
      }//end of loop ctr on the plane
      std::cout << "Number of contour points of plane #" << p << ": " << nfound << std::endl;
    }//end of plane loop


    // We make 3D points based on 2 plane crossings: we check if point inside the other plane
    // (this is for non-horizontal tracks. for that we have to use the edges)
    // we remove close points as well
    int npts = contour_points.size();
    imgcoords.resize(4,0);
    pos3d.resize(3,0);
    for (int a=0; a<npts; a++) {
      ctr_pt_t& pta = contour_points[a];
      for (int b=a+1; b<npts; b++) {
	ctr_pt_t& ptb = contour_points[b];
	if ( pta.plane==ptb.plane )
	  continue; // don't match the same planes, bruh

	int otherp;
	int otherw;
	int crosses;
	std::vector<float> xsec_zy(2,0);
	larcv::UBWireTool::getMissingWireAndPlane( pta.plane, pta.maxqc, ptb.plane, ptb.maxqc, otherp, otherw, xsec_zy, crosses );

	// bad crossing or out of wire range
	if ( crosses==0 || otherw<0 || otherw>=img_v[otherp].meta().max_x() )
	  continue;

	// check for charge or badch
	if ( img_v[otherp].pixel( row, otherw )>10 || badch_v[otherp].pixel(row, otherw)>0 ) {
	  // good point
	  imgcoords[0] = row;
	  imgcoords[pta.plane+1] = pta.maxqc;
	  imgcoords[ptb.plane+1] = ptb.maxqc;
	  imgcoords[otherp+1]    = otherw;

	  pos3d[0] = img_v.front().meta().pos_x( row ); // tick
	  pos3d[1] = xsec_zy[1];
	  pos3d[2] = xsec_zy[0];

	  std::cout << "Produced 3D Point: imgcoords(" << imgcoords[0] << "," << imgcoords[1] << "," << imgcoords[2] << "," << imgcoords[3] << ")"
		    << " pos3d=(" << pos3d[0] << "," << pos3d[1] << "," << pos3d[2] << ")" << std::endl;

	}
	
      }
    }
    
  }
  
  // =================================================================================
  // ALGO METHODS

  ContourAStarCluster ContourAStarClusterAlgo::makeCluster( const std::vector<float>& pos3d, const std::vector<larcv::Image2D>& img_v,
							    const std::vector<larcv::Image2D>& badch_v,
							    const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
							    const float max_dist2cluster ) {

    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    ContourAStarCluster cluster = makeSeedClustersFrom3DPoint( pos3d, img_v, plane_contours_v, max_dist2cluster );

    // now we enter the buiding loop
    int nloopsteps = 1;
    int iloop = 0;
    while( iloop<nloopsteps ) {
      iloop++;
      
      // first get the overlapping time region
      std::vector<int> timerange = cluster.getOverlappingRowRange();
      timerange[0]++;
      timerange[1]--;
      std::cout << "Time range: [" << timerange[0] << "," << timerange[1] << "]" << std::endl;
      
      // draw line
      cv::line( cluster.m_cvimg_debug, cv::Point(0,timerange[0]), cv::Point(3455,timerange[0]), cv::Scalar(255,255,255), 1 );
      cv::line( cluster.m_cvimg_debug, cv::Point(0,timerange[1]), cv::Point(3455,timerange[1]), cv::Scalar(255,255,255), 1 );      

      // we scan at the time range for a good 3D point (or points?)
      std::vector<int> min_imgcoords;
      std::vector<float> min_pos3d;
      std::vector<int> max_imgcoords;
      std::vector<float> max_pos3d;
      cluster.getCluster3DPointAtTimeTick( timerange[0], img_v, badch_v, true, min_imgcoords, min_pos3d );
      cluster.getCluster3DPointAtTimeTick( timerange[1], img_v, badch_v, true, max_imgcoords, max_pos3d );

      cv::circle( cluster.m_cvimg_debug, cv::Point(min_imgcoords[1],timerange[0]), 2, cv::Scalar(0,255,255,255), -1 );
      cv::circle( cluster.m_cvimg_debug, cv::Point(min_imgcoords[2],timerange[0]), 2, cv::Scalar(0,255,255,255), -1 );
      cv::circle( cluster.m_cvimg_debug, cv::Point(min_imgcoords[3],timerange[0]), 2, cv::Scalar(0,255,255,255), -1 );

      cv::circle( cluster.m_cvimg_debug, cv::Point(max_imgcoords[1],timerange[1]), 2, cv::Scalar(0,255,255,255), -1 );
      cv::circle( cluster.m_cvimg_debug, cv::Point(max_imgcoords[2],timerange[1]), 2, cv::Scalar(0,255,255,255), -1 );
      cv::circle( cluster.m_cvimg_debug, cv::Point(max_imgcoords[3],timerange[1]), 2, cv::Scalar(0,255,255,255), -1 );
      
      // use astar between these points!
      larlitecv::AStar3DAlgoConfig astar_cfg;
      astar_cfg.verbosity = 1;      
      astar_cfg.min_nplanes_w_hitpixel = 2;
      astar_cfg.min_nplanes_w_charge = 2;
      astar_cfg.astar_threshold.resize(3,5);
      astar_cfg.astar_neighborhood.resize(3,3);
      astar_cfg.restrict_path = true;
      astar_cfg.path_restriction_radius = 10.0;
      astar_cfg.accept_badch_nodes = true;
      astar_cfg.astar_start_padding = 3;
      astar_cfg.astar_end_padding = 3;
      astar_cfg.lattice_padding = 3;

      const larcv::ImageMeta& meta = img_v.front().meta();
      larlitecv::AStar3DAlgo algo( astar_cfg );
      std::vector<larlitecv::AStar3DNode> path;
      std::vector< int > start_cols;
      std::vector< int > end_cols;
      for (int i=0; i<3; i++) {
	start_cols.push_back( min_imgcoords[i+1] );
	end_cols.push_back(   max_imgcoords[i+1] );
      }
      int goalhit = 0;
      path = algo.findpath( img_v, badch_v, badch_v, timerange[0], timerange[1], start_cols, end_cols, goalhit );
      
      std::cout << "Goal hit: " << goalhit << " pathsize=" << path.size() << std::endl;
      for ( auto &node : path ) {
	std::vector<int> imgcoords = larcv::UBWireTool::getProjectedImagePixel( node.tyz, meta, img_v.size() );
	imgcoords[0] = meta.row( node.tyz[0] );
	for (int i=0; i<3; i++) {
	  cv::circle( cluster.m_cvimg_debug, cv::Point( imgcoords[i+1], imgcoords[0] ), 1, cv::Scalar(255,0,255), -1 );
	}
      }
      
      // evaluate track.
      
      // if good, extend track
      // we do this by getting a 3d fit of the line, extending past the end points, absorbing new clusters on all three planes
      std::vector< std::vector<float> > v3d;
      for (int idx=0; idx<(int)path.size(); idx++) {
	auto &node = path[idx];
	std::vector<float> v3 = node.tyz;
	// convert tick to x
	v3[0] = (v3[0]-3200.0)*cm_per_tick;
	v3d.push_back( v3 );
      }

      extendClusterUsingAStarPath( cluster, v3d, img_v, 10.0, 10.0, 1.0 );
      
      // absorb clusters
      break;
    }
    
    return cluster;
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

    std::cout << __FILE__ << ":" << __LINE__ << " Seed Point: "
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

  void ContourAStarClusterAlgo::extendClusterUsingAStarPath( ContourAStarCluster& cluster, const std::vector< std::vector<float> >& path3d,
							     const std::vector<larcv::Image2D>& img_v,
							     const float distfromend, const float distextended, const float stepsize ) {

    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
    const larcv::ImageMeta& meta = img_v.front().meta();
    
    // make a Point array list for opencv
    // we gather points some distance from the ends
    std::vector< cv::Point3f > point3d_start;
    std::vector< cv::Point3f > point3d_end;    

    const std::vector<float>& start = path3d.front();
    const std::vector<float>& end   = path3d.back();

    std::cout << "extend path w/ size=" << path3d.size() << std::endl;
    
    for (int idx=0; idx<(int)path3d.size(); idx++) {
      const std::vector<float>& v3 = path3d[idx];

      float dist_start = 0;
      float dist_end   = 0;
      for (int i=0; i<3; i++) {
	dist_start += (start[i]-v3[i])*(start[i]-v3[i]);
	dist_end   += (end[i]-v3[i])*(end[i]-v3[i]);
      }
      dist_start = sqrt(dist_start);
      dist_end   = sqrt(dist_end);

      cv::Point3f pt( v3[0], v3[1], v3[2] );
      std::cout << "path [" << idx << "] start_dist=" << dist_start << " end_dist=" << dist_end << std::endl;
      if ( dist_start<distfromend ) {
	point3d_start.push_back( pt );
      }
      if ( dist_end<distfromend ) {
	point3d_end.push_back( pt );
      }
    }//end of path loop
    

    // fit line using opencv function
    std::vector<float> startline; // 6 parameters (vx, vy, vz, x0, y0, z0 )
    cv::fitLine( point3d_start, startline, CV_DIST_L2, 0.0, 0.01, 0.01 );
    std::vector<float> endline; // 6 parameters (vx, vy, vz, x0, y0, z0 )
    cv::fitLine( point3d_end, endline, CV_DIST_L2, 0.0, 0.01, 0.01 );

    // get direction from start to end to help us orient the lines
    std::vector<float> cluster_dir(3); // min -> max
    float cluster_dir_norm = 0.;
    for (int i=0; i<3; i++) {
      cluster_dir[i] = end[i] - start[i];
      cluster_dir_norm += cluster_dir[i]*cluster_dir[i];
    }
    cluster_dir_norm = sqrt(cluster_dir_norm);
    for (int i=0; i<3; i++)
      cluster_dir[i] /= cluster_dir_norm;

    float cosstart = 0;
    float cosend   = 0;
    for (int i=0; i<3; i++) {
      cosstart += cluster_dir[i]*startline[i];
      cosend   += cluster_dir[i]*endline[i];
    }

    // start should be negative, end should be positive
    if ( cosstart>0 ) {
      for (int i=0; i<3; i++)
	startline[i] *= -1.0;
    }
    if ( cosend<0 ) {
      for (int i=0; i<3; i++)
	endline[i] *= -1.0;
    }

    // walk along the line
    //  check which clusters are touched by it
    int numsteps = distextended/stepsize;
    std::cout << "extend line with " << numsteps << " steps" << std::endl;

    for (int isteps=1; isteps<=numsteps; isteps++ ) {

      std::vector<float> stepposmax(3);
      std::vector<float> stepposmin(3);
      for (int i=0; i<3; i++) {
	stepposmax[i] = float(isteps)*stepsize*endline[i]   + end[i];
	stepposmin[i] = float(isteps)*stepsize*startline[i] + start[i];
      }
      
      // x to tick
      stepposmax[0] = stepposmax[0]/cm_per_tick + 3200.0;
      stepposmin[0] = stepposmin[0]/cm_per_tick + 3200.0;
      
      if ( stepposmax[0]>2400 && stepposmax[0]<8448 && stepposmax[1]>-117 && stepposmax[1]<117 && stepposmax[2]>0.3 && stepposmax[2]<1030 ) {
	// project back into the image
	std::vector<int> imgcoordsmax = larcv::UBWireTool::getProjectedImagePixel( stepposmax, meta, img_v.size() );
	imgcoordsmax[0] = meta.row( stepposmax[0] );
	for (int p=0; p<3; p++) {
	  cv::circle( cluster.m_cvimg_debug, cv::Point( imgcoordsmax[p+1], imgcoordsmax[0] ), 1, cv::Scalar(0,255,0), -1 );
	}
      }
      else {
	std::cout << " stepmax out of bounds: (" << stepposmax[0] << "," << stepposmax[1] << "," << stepposmax[2] << ")" << std::endl;
      }
      if ( stepposmin[0]>2400 && stepposmin[0]<8448 && stepposmin[1]>-117 && stepposmin[1]<117 && stepposmin[2]>0.3 && stepposmin[2]<1030 ) {	
	std::vector<int> imgcoordsmin = larcv::UBWireTool::getProjectedImagePixel( stepposmin, meta, img_v.size() );
	imgcoordsmin[0] = meta.row( stepposmin[0] );
	for (int p=0; p<3; p++) {
	  cv::circle( cluster.m_cvimg_debug, cv::Point( imgcoordsmin[p+1], imgcoordsmin[0] ), 1, cv::Scalar(0,255,0), -1 );
	}
      }
      else {
	std::cout << " stepmin out of bounds: (" << stepposmin[0] << "," << stepposmin[1] << "," << stepposmin[2] << ")" << std::endl;
      }
      
    }//end of step loop
      

    return;
  }//end of function

  
}
