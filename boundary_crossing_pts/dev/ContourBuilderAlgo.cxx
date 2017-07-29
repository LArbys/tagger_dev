#include "ContourBuilderAlgo.h"
#include <assert.h>
#include <sstream>
#include <exception>

#include "UBWireTool/UBWireTool.h"

#ifdef USE_OPENCV
#include "CVUtil/CVUtil.h"
#endif


namespace larlitecv {

  // =======================================================================================================
  ContourBuilderAlgo::ContourBuilderAlgo() {};

  ContourBuilderAlgo::~ContourBuilderAlgo() {};

  void ContourBuilderAlgo::buildCluster( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v,
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
    
  }

  bool ContourBuilderAlgo::extendClusterGroup( ContourCluster& cluster, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
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
	// order by distance
	float nearfardivide = 50.0;
	if ( dist< nearfardivide && rhs.dist>nearfardivide ) {
	  return true;
	}
	else if ( dist>nearfardivide && rhs.dist<nearfardivide )
	  return false;
	else if ( dist<nearfardivide && rhs.dist<nearfardivide ) {
	  if ( cosine < rhs.cosine )
	    return true;
	  else
	    return false;
	}
	else if ( dist>=nearfardivide && rhs.dist>=nearfardivide ) {
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

      const cv::Point& earlypt = cluster.earlyEnd[p].back(); // early end on the plane
      const cv::Point& latept = cluster.lateEnd[p]; // late end on the plane

      std::vector< ClusterMatchScore > early_extension_candidates_v;
      std::vector< ClusterMatchScore >  late_extension_candidates_v;
      
      for (int idx=0; idx<(int)plane_contours_v[p].size(); idx++) {

	// we already used this cluster
	if ( cluster.indices[p].find(idx)!=cluster.indices[p].end() )
	  continue;

	// get contour
	const ContourShapeMeta& contour = plane_contours_v[p][idx];

	if ( contour.size()<10 )
	  continue;

	// get closet end to early end of current cluster
	float early2start = sqrt( (contour.getFitSegmentStart().x-earlypt.x)*(contour.getFitSegmentStart().x-earlypt.x)
				  +(contour.getFitSegmentStart().y-earlypt.y)*(contour.getFitSegmentStart().y-earlypt.y) );
	float early2end = sqrt( (contour.getFitSegmentEnd().x-earlypt.x)*(contour.getFitSegmentEnd().x-earlypt.x)
				+(contour.getFitSegmentEnd().y-earlypt.y)*(contour.getFitSegmentEnd().y-earlypt.y) );
	float early2startcos = 0;
	float early2endcos   = 0;
	for (int i=0; i<2; i++) {
	  early2startcos += contour.getStartDir()[i]*cluster.earlyDir[p].back()[i];
	  early2endcos   += contour.getEndDir()[i]*cluster.earlyDir[p].back()[i];
	}
		
	std::vector< cv::Point > earlypair;
	float mindist = 0;
	float mincos  = 0;
	std::vector<float> earlydir(2,0);
	if ( early2start<early2end ) {
	  earlypair.push_back( earlypt );
	  earlypair.push_back( contour.getFitSegmentStart() );
	  mindist = early2start;
	  mincos  = early2startcos;
	  earlydir = contour.getStartDir();
	}
	else {
	  earlypair.push_back( earlypt );
	  earlypair.push_back( contour.getFitSegmentEnd() );
	  mindist = early2end;
	  mincos  = early2endcos;
	  earlydir = contour.getEndDir();	  
	}

	// cos between connecting line and contour direction
	float econnectdir[2];
	econnectdir[0] = earlypair[1].x-earlypair[0].x;
	econnectdir[1] = earlypair[1].y-earlypair[0].y;
	float econnectnorm = 0;
	for (int i=0; i<2; i++)
	  econnectnorm += econnectdir[i]*econnectdir[i];
	econnectnorm = sqrt(econnectnorm);
	float econnectcos = 0;
	for (int i=0; i<2; i++) {
	  econnectdir[i] /= econnectnorm;
	  econnectcos += econnectdir[i]*earlydir[i];
	}

	// save the candidate info
	ClusterMatchScore earlymatchcand( idx, mindist, mincos, econnectcos );
	early_extension_candidates_v.emplace_back( earlymatchcand );


	// get closet end to late end of current cluster
	float late2start = sqrt( (contour.getFitSegmentStart().x-latept.x)*(contour.getFitSegmentStart().x-latept.x)
				  +(contour.getFitSegmentStart().y-latept.y)*(contour.getFitSegmentStart().y-latept.y) );
	float late2end = sqrt( (contour.getFitSegmentEnd().x-latept.x)*(contour.getFitSegmentEnd().x-latept.x)
				+(contour.getFitSegmentEnd().y-latept.y)*(contour.getFitSegmentEnd().y-latept.y) );
	float late2startcos = 0;
	float late2endcos   = 0;
	for (int i=0; i<2; i++) {
	  late2startcos += contour.getStartDir()[i]*cluster.lateDir[p][i];
	  late2endcos   += contour.getEndDir()[i]*cluster.lateDir[p][i];
	}
	
	std::vector< cv::Point > latepair;
	std::vector<float> latedir(2,0);
	float latemindist = 0;
	float latemincos  = 0;
	if ( late2start<late2end ) {
	  latepair.push_back( latept );
	  latepair.push_back( contour.getFitSegmentStart() );
	  latemindist = late2start;
	  latemincos  = late2startcos;
	  latedir = contour.getStartDir();
	  //ClusterMatchScore matchcand( idx, late2start, late2startcos, 0.0 );
	  //late_extension_candidates_v.emplace_back( matchcand );	  	  
	}
	else {
	  latepair.push_back( latept );
	  latepair.push_back( contour.getFitSegmentEnd() );
	  latemindist = late2end;
	  latemincos  = late2endcos;
	  latedir = contour.getEndDir();
	  //ClusterMatchScore matchcand( idx, late2end, late2endcos, 0.0 );
	  //late_extension_candidates_v.emplace_back( matchcand );	  	  
	}

	// cos between connecting line and contour direction
	float lconnectdir[2];
	lconnectdir[0] = latepair[1].x-latepair[0].x;
	lconnectdir[1] = latepair[1].y-latepair[0].y;
	float lconnectnorm = 0;
	for (int i=0; i<2; i++)
	  lconnectnorm += lconnectdir[i]*lconnectdir[i];
	lconnectnorm = sqrt(lconnectnorm);
	float lconnectcos = 0;
	for (int i=0; i<2; i++) {
	  lconnectdir[i] /= lconnectnorm;
	  lconnectcos += lconnectdir[i]*latedir[i];
	}

	// save the candidate info
	ClusterMatchScore latematchcand( idx, latemindist, latemincos, lconnectcos );
	late_extension_candidates_v.emplace_back( latematchcand );
	
      }//end of loop over plane conours

      std::sort( early_extension_candidates_v.begin(), early_extension_candidates_v.end() );
      std::sort( late_extension_candidates_v.begin(), late_extension_candidates_v.end() );      

      std::cout << "Plane " << p << " Early Extension matching" << std::endl;
      for (int iex=0; iex<(int)early_extension_candidates_v.size(); iex++) {
	auto& candidate = early_extension_candidates_v[iex];
	std::cout << " [" << candidate.idx << "] dist=" << candidate.dist << " cos=" << candidate.cosine << " triscore=" << candidate.triscore << std::endl;
	if ( candidate.cosine<-0.7 && candidate.dist<200.0 && (candidate.triscore<-0.5 || candidate.dist<10) ) {
	  std::cout << "  Match: " << candidate.cosine << " " << candidate.triscore << " " << candidate.dist << std::endl;
	  // append best to early match
	  auto const& early_contour = plane_contours_v[p][candidate.idx];
	  cluster.earlyContours[p].push_back( &early_contour );
	  cluster.earlyDir[p].push_back( early_contour.getStartDir() );
	  cluster.earlyEnd[p].push_back( early_contour.getFitSegmentStart() );
	  cluster.indices[p].insert( candidate.idx );
	  cluster_extended = true;
	  break;
	}
      }

      std::cout << "Plane " << p << " Late Extension matching" << std::endl;
      for (int iex=0; iex<(int)late_extension_candidates_v.size(); iex++) {
	auto& candidate = late_extension_candidates_v[iex];

	if ( candidate.cosine<-0.7 && candidate.dist<200.0 && (candidate.triscore<-0.5 || candidate.dist<10) ) {
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


  int ContourBuilderAlgo::getIndexOfContainingContour( const int row, const int col, const std::vector<ContourShapeMeta>& contours_v, int min_cluster_size, float dist_tolerance ) {
    // first test bbox (faster test)
    // then do point poly test
    
    bool contains = false;
    int containing_idx = -1;
    cv::Point imgpt(col,row);
    for (int idx=0; idx<(int)contours_v.size(); idx++) {
      
      if ( contours_v[idx].size()<min_cluster_size )
	continue;
      
      // test imgpt
      bool bboxcontains = contours_v[idx].getBBox().contains( imgpt );
      if ( !bboxcontains )
	continue;
      
      // more detailed test
      double dist = cv::pointPolygonTest( contours_v[idx], imgpt, false );
      if ( dist>=-fabs(dist_tolerance) )
	return idx;
    }
    return -1;
  }
  
  bool ContourBuilderAlgo::ratchetCluster( ContourCluster& cluster, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
					   const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v ) {
    // what we do is we find the tmax and tmin of the clusters across planes
    // we then extend the direction of the current contours to that time until we run into another contour.
    // we build until we reach. we keep going, ratcheting until no more clusters.
    
    bool cluster_extended = false;

    float row_max = 0;
    float row_min = (float)img_v.front().meta().rows();

    // we get the time bounds of the current contour cluster
    for ( size_t p=0; p<img_v.size(); p++ ) {

      // get the current ends of the cluster
      const cv::Point& earlypt = cluster.earlyEnd[p].back(); // early end on the plane
      const cv::Point& latept  = cluster.lateEnd[p];         // late end on the plane

      if ( earlypt.y>row_max ) {
	row_max = earlypt.y;
      }
      if ( earlypt.y<row_min ) {
	row_min = earlypt.y;
      }
      if ( latept.y>row_max ) {
	row_max = latept.y;
      }
      if ( latept.y<row_min ) {
	row_min = latept.y;
      }
    }

    std::cout << __FILE__ << ":Ratchet Cluster bounds [" << row_min << "," << row_max << "]" << std::endl;
    
    // now for each plane, we extrapolate out from the ends, to the end extremes
    // we step until we run in the a new contour
    for ( size_t p=0; p<img_v.size(); p++ ) {

      // early extension
      const cv::Point& earlypt = cluster.earlyEnd[p].back();
      std::vector<float> earlydir = cluster.earlyDir[p].back();

      if ( earlypt.y==row_min )
	continue; // not extending in this plane

      if ( earlydir[1]==0 )
	continue; // outside the use case unfortunately... need a separate horizontal extension check
      
      float drow = row_min-earlypt.y;
      int nsteps = drow/(earlydir[1]/10.0);

      int lastrow = earlypt.y;
      int lastcol = earlypt.x;
      std::set<int> onpath_indices_set; // for dupilicate search
      std::vector<int> onpath_indices_v; // for listing in order of intersection
      for (int istep=0; istep<=nsteps; istep++) {
	cv::Point testpt( earlypt.x+(float(istep)*earlydir[0]/10.0), earlypt.y+(float(istep)*earlydir[1]/10.0) );
	if ( (int)testpt.x==lastcol && (int)testpt.y==lastrow )
	  continue;

	// update pt
	lastrow = testpt.y;
	lastcol = testpt.x;

	int idx = getIndexOfContainingContour( testpt.y, testpt.x, plane_contours_v[p], 10, 5.0 );
	if ( idx<0 || onpath_indices_set.find(idx)!=onpath_indices_set.end() )
	  continue;

	onpath_indices_set.insert(idx);
	onpath_indices_v.push_back(idx);
      }

      std::cout << __FILE__ << ":RatchetCluster number of on-path candidates " << onpath_indices_v.size() << std::endl;
      
      // evaluate quality of candidates

      for ( auto& idx : onpath_indices_v ) {

	auto const& early_contour = plane_contours_v[p][idx];

	float cosine = 0;
	for (int i=0; i<2; i++)
	  cosine += early_contour.getStartDir()[i]*earlydir[i];

	std::cout << __FILE__ << ":   candidate [" << idx << "] cosine=" << cosine << std::endl;
	
	if ( cosine>0.8 ) {
	  cluster.earlyContours[p].push_back( &early_contour );
	  cluster.earlyDir[p].push_back( early_contour.getStartDir() );
	  cluster.earlyEnd[p].push_back( early_contour.getFitSegmentStart() );
	  cluster.indices[p].insert( idx );
	  cluster_extended = true;
	  break;
	}
      }
    }//end of p loop

    return cluster_extended;
  }//end of extension method
  

  std::vector<float> ContourBuilderAlgo::calculateContourIntersection( const std::vector< cv::Point >& cnt1, const std::vector< cv::Point >& cnt2 ) {
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

  bool ContourBuilderAlgo::buildContourGraph( ContourCluster& cluster, const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
					      const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v, std::vector< larcv::Image2D >& clusterpix_v ) {
    bool cluster_extended = false;
    for (int p=0; p<3; p++) {

      // we make a node bank
      std::vector< ContourGraphNode* > nodebank;
      for (int idx=0; idx<(int)plane_contours_v[p].size(); idx++) {
	nodebank.push_back( new ContourGraphNode(idx) );
      }

      // grab the seed contour
      int seedidx = *(cluster.indices[p].begin());
      ContourGraphNode* seednode = nodebank[seedidx];
      RecursiveSearch( seednode, plane_contours_v[p], nodebank );
      
      cv::Mat cvimg = larcv::as_gray_mat( img_v[p], 8.0, 256.0, 1.0 );
      cv::Mat cvrgb = larcv::as_mat_greyscale2bgr( img_v[p], 10.0, 100.0 );

      for (int idx=0; idx<(int)plane_contours_v[p].size(); idx++) {
	if ( nodebank[idx]->visited ) {
	  auto const& contour = plane_contours_v[p].at( idx );
	  std::vector< std::vector<cv::Point> >  contour_v;
	  contour_v.push_back( contour );
	  cv::drawContours( cvrgb, contour_v, 0, cv::Scalar(0,0,255,255), -1 );
	}
	delete nodebank[idx];
	nodebank[idx] = NULL;
      }

      std::stringstream imgname;
      imgname << "contourgraph_plane" << p << ".png";
      cv::imwrite( imgname.str(), cvrgb );

    }
    
    return cluster_extended;
  }//end of extension method

  
}
