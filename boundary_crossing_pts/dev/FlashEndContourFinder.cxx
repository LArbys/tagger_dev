#include "FlashEndContourFinder.h"

#include <vector>

// larlite
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

// larcv
#include "UBWireTool/UBWireTool.h"

namespace larlitecv {

  FlashEndContourFinder::FlashEndContourFinder( SearchMode_t mode, const FlashEndContourFinderConfig& config )
    : fSearchMode(mode),fNPMTs(32),fConfig(config)
  {

  }

  FlashEndContourFinder::FlashEndContourFinder( SearchMode_t mode, const larcv::PSet& pset )
    : fSearchMode(mode),fNPMTs(32),fConfig( FlashEndContourFinderConfig::FromPSet( pset ) )
  {
    
  }

  bool FlashEndContourFinder::flashMatchTrackEnds( const std::vector< larlite::event_opflash* >& opflash_vv, 
						   const std::vector<larcv::Image2D>& img_v,
						   const std::vector<larcv::Image2D>& badch_v,
						   const std::vector< std::vector<larlitecv::ContourShapeMeta> >& plane_contours_v,
						   std::vector< BoundarySpacePoint >& trackendpts,
						   std::vector< int > endpoint_flash_idx) {
    // Primary method to be called by the user
    // -----------------------------------------
    // input
    // -----
    // opflash_vv: vector of event opflash containers
    // img_v: vector of images, each entry is a plane
    // badch_v: vector of images with badch's marked, each entry is a plane
    // plane_contours_v: list of contours for each plane. this object is results of BMTCV::analyzeImages
    
    // output
    // ------
    //  trackendpts: list of vector of end pts. each (inner vector) is point on each plane.
    //  endpoint_flash_idx: (unrolled) index for opflash_vv that each spacepoint was matched to
    
    const int nplanes = img_v.size();
    
    if ( nplanes==0 )
      return false;
    
    if ( fConfig.verbosity>0 ) {
      std::cout << "Begin FlashEndContourFinder::flashMatchTrackEnds [verbsity=" << fConfig.verbosity << "]" << std::endl;
      std::cout << "  drift_v = " << fConfig.drift_velocity << " (fcl) vs." << larutil::LArProperties::GetME()->DriftVelocity() << std::endl;
      std::cout << "  drift distance =" << fConfig.drift_distance << std::endl;
      std::cout << "  number of opflash containers: " << opflash_vv.size() << std::endl;
      for ( int i=0; i<(int)opflash_vv.size(); i++) {
	std::cout << "    #" << i << ": " << opflash_vv[i]->size() << " flashes" << std::endl;
      }
    }
    
    // get a meta
    const larcv::ImageMeta& meta = img_v.at(0).meta();

    // uboone geometry data
    const larutil::Geometry* pgeo = larutil::Geometry::GetME();

    int unrolled_flash_index = -1;
    int ntrackendpts = trackendpts.size();
    
    // loop over all the flash containers
    for ( auto& ptr_event_flash : opflash_vv ) {
      // loop over flashes
      
      int opflash_idx = -1;
      
      for ( auto& opflash : *ptr_event_flash ) {
	opflash_idx++;	
	unrolled_flash_index++;
	
        // determine time of flash
        float tick_target = 0;
        float flash_tick = 0;
        larlitecv::BoundaryEnd_t point_type;
        std::string modename;
        if ( fSearchMode==kAnode ) {
          flash_tick = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick + fConfig.anode_drift_tick_correction;;
          tick_target = flash_tick;
          point_type = larlitecv::kAnode;
          modename = "anode";
        }
        else if ( fSearchMode==kCathode ) {
          flash_tick = fConfig.trigger_tick + opflash.Time()/fConfig.usec_per_tick;
          tick_target = flash_tick + fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick+fConfig.cathode_drift_tick_correction;
          point_type = larlitecv::kCathode;          
          modename = "cathode";
        }
        else if ( fSearchMode==kOutOfImage ) {
          flash_tick  = opflash.Time();
          tick_target = flash_tick; // dummy opflash gives first or last tick of image
          point_type = larlitecv::kImageEnd;          
          modename = "image end";
        }
        else {
          //std::cout << "[ERROR] wrong search mode" << std::endl;
          return false;
        }
        
        // check if the opflash time occurs within the image
        if ( tick_target<=meta.min_y() || tick_target>=meta.max_y() ) {
	  if ( fConfig.verbosity>0 ) {
	    std::cout << "============================================================================================" << std::endl;	    
	    std::cout << " [opflash search] Op Flash out of time. tick_target=" << tick_target << " flash_tick=" << flash_tick
		      << " Image bounds: [" << meta.min_y() << "," << meta.max_y() << "]" << std::endl;
	  }
          continue;
	}

        // first find the weighted mean and total q     
        float qtot = 0;
        float z_weighted = 0.;
        for (int iopdet=0; iopdet<fNPMTs; iopdet++) {
	  std::vector<double> pmtpos(3,0);
	  pgeo->GetOpDetPosition( iopdet, pmtpos );
          z_weighted += opflash.PE( iopdet )*pmtpos[2];
          qtot += opflash.PE( iopdet );
        }
        if ( qtot>0 ) {
          z_weighted /= qtot;
        }

        // to set the range, we find the first hit above threshold from the mean
        float min_dist_z = 1e9;
        float max_dist_z = 0;
	float z_max = 0;
	float q_max = -1;
        for (int iopdet=0; iopdet<fNPMTs; iopdet++) {
	  std::vector<double> pmtpos(3,0);
	  pgeo->GetOpDetPosition( iopdet, pmtpos );	  
          float pe = opflash.PE(iopdet);
          float dist = pmtpos[2]-z_weighted;
          if ( pe>5.0 ) {
            if ( dist<0 && min_dist_z>dist ) min_dist_z = dist;
            else if ( dist>0 && max_dist_z<dist ) max_dist_z = dist;
	    if ( pe>q_max ) {
	      z_max = pmtpos[2];
	      q_max = pe;
	    }
          }
        }

        // extend by some factor (example 10%). This is to ensure acceptance of tracks.
        float zwidth = fabs(max_dist_z)+fabs(min_dist_z);
        float extension = zwidth*fConfig.flash_zrange_extension*0.5;

        std::vector<float> z_range = { z_weighted+min_dist_z, z_weighted+max_dist_z}; // be more intelligent later
        std::vector<float> y_range = { -120.0, 120.0 };
	
	if ( fSearchMode==kCathode ) {
	  // flash position is not super helpful for cathode muons. biases towards position of end point
	  z_range[0] -= extension;
	  z_range[1] += extension;
	  if ( fConfig.verbosity>1 ) {
	    std::cout << "extending cathode window by " << extension << std::endl;
	  }
	}
        if ( fSearchMode==kOutOfImage ) {
          // accept all
          z_range[0] = 0;
          z_range[1] = 1100;
        }
        else {
          // just go to ends in these cases
          if ( z_range[0]<55.0 ) z_range[0] = 0; 
          if ( z_range[1]>980.0 ) z_range[1] = 1100;
        }
	
        int row_target = meta.row( tick_target );
        
        if ( fConfig.verbosity>0 ) {
          std::cout << "============================================================================================" << std::endl;
          std::cout << "[Opflash search] " << modename << " mode" << std::endl;
          std::cout << "  opflash: "
                    << " flash_tick=" << flash_tick
                    << " tick_target=" << tick_target
                    << " row_target=" << row_target
                    << " qtot= " << qtot
		    << " zmax=" << z_max
		    << " qmax=" << q_max
                    << " z_range=[" << z_range[0] << "," << z_range[1] << "] "
                    << " w_range=[" << int(z_range[0]/0.3/meta.pixel_width()) << "," << int(z_range[1]/0.3/meta.pixel_width()) << "] "
                    << " drift_t=" << fConfig.drift_distance/fConfig.drift_velocity/fConfig.usec_per_tick+200.0 << " ticks"
                    << std::endl;
        }

	// =================================================================================
	// Flash search method: we use the charge clusters at a given time

	findCandidateEndsAtTick( tick_target, img_v, badch_v, plane_contours_v, trackendpts );
	
	// =================================================================================

	// add the unrolled flash index to the end point
	for (int iendpt=ntrackendpts; iendpt<(int)trackendpts.size(); iendpt++) {
	  trackendpts.at(iendpt).setFlashIndex( unrolled_flash_index );
	  std::cout << __FILE__ << ":" << __LINE__ << " set flash index=" << unrolled_flash_index << " (ntrackendpts=" << ntrackendpts << ")" << std::endl;
	}
	ntrackendpts = trackendpts.size();
	
      }//end of opflashes loop
    }//end of opflash_vv

    return true;
  }
  
  BoundaryEnd_t FlashEndContourFinder::SearchModeToEndType( SearchMode_t mode ) {
    switch(mode) {
      case FlashEndContourFinder::kAnode:
        return larlitecv::kAnode;
        break;
      case FlashEndContourFinder::kCathode:
        return larlitecv::kCathode;
        break;
      case FlashEndContourFinder::kOutOfImage:
        return larlitecv::kImageEnd;
        break;
      default:
        throw std::runtime_error("FlashEndContourFinder::SearchModeToEndType unknown search mode.");
        break;
    }
    return larlitecv::kUndefined;
  }

  bool FlashEndContourFinder::findCandidateEndsAtTick( const float tick, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
						       const std::vector< std::vector<larlitecv::ContourShapeMeta> >& plane_contours_v,
						       std::vector< BoundarySpacePoint >& trackendpts ) {
    // -----------------------------------------------------------    
    // Given a tick to search, we look for candidate track ends
    // 
    // we use the contour machinery steps
    // - we look for contours on all three planes
    // - each 2 plane conbination gives us a seed point
    // - only keep combinations which have a 3-plane match
    // - these serve as seed points to the contour-astar-cluster algo
    //  
    
    // inputs
    // ------
    // tick: time at which to search for track ends
    // img_v: vector of images, entries correspond to plane view
    // badch_v: vector of images, entries correspond to badch view
    // plane_contours_v: vector of contours for each plane. made by TaggerContourTools/BMTCV::analyzeImages
    //
    // outputs
    // -------
    // trackendpts: container where candidate end points will be stored

    const larcv::ImageMeta& meta = img_v.front().meta();

    int row = (int)meta.rows();

    // get list of contours for each plane
    // and store info on where the tick crosses the contour
    std::vector< std::vector<CtrXing_t> > plane_xing_v(img_v.size());
    
    generateKeyPointList( row, img_v, plane_contours_v, plane_xing_v );

    // now that crossings gathered, we take combinations in order to make candidate space points
    // candidates require that the intersection formed by the wires are
    // (1) within the detector
    // (2) 3-plane match
    // we will use the max-q point most of the time, but for wide-contours, (from horizontal tracks)
    //   we test points at the min and max locations (possible 4-points created)

    // using the key crossing points, we find candidate 3d points
    std::vector< std::vector<float> > candidate_tyz_points_v;
    findCandidateSpacePoints( plane_xing_v, row, img_v, badch_v, candidate_tyz_points_v );

    // we feed this to the contour end point finding machinery
    
    return true;
  }

  void FlashEndContourFinder::generateKeyPointList( const int row, const std::vector<larcv::Image2D>& img_v,
						    const std::vector< std::vector<larlitecv::ContourShapeMeta> >& plane_contours_v,
						    std::vector< std::vector<FlashEndContourFinder::CtrXing_t> >& plane_xing_v ) {
						    
    plane_xing_v.resize(img_v.size());
    
    for (size_t p=0; p<img_v.size(); p++) {
      plane_xing_v[p].clear();
      
      const std::vector<larlitecv::ContourShapeMeta>& plane_ctrs = plane_contours_v[p];

      CtrXing_t xing;
      xing.mincol = -1;
      xing.maxcol = -1;
      xing.maxqcol = -1;
      xing.maxq = -1;
      xing.pctr = NULL;
      
      for ( auto const& ctr : plane_ctrs ) {

	// does the row cross the ctr?
	if ( row<ctr.getMinY() || row>ctr.getMaxY() )
	  continue;

	// so contour lives in the same time-range, where does it cross?, where is the max charge point?
	for (size_t c=ctr.getMinX(); c<=ctr.getMaxX(); c++) {
	  cv::Point testpt( c, row );
	  double dist = cv::pointPolygonTest( ctr, testpt, false );
	  if ( dist<0 )
	    continue; // outside contour
	  
	  // otherwise inside
	  float pixval = img_v[p].pixel(row,c);
	  if ( pixval<fConfig.pixel_value_threshold[p] )
	    continue; // below threshold
	  
	  if ( xing.mincol<0 || xing.mincol>c )
	    xing.mincol = c;
	  if ( xing.maxcol<0 || xing.maxcol<c )
	    xing.maxcol = c;
	  if ( xing.maxq<0 || xing.maxq<pixval ) {
	    xing.maxq = pixval;
	    xing.maxqcol = c;
	  }
	}
	
	// save the crossing for the contour
	if ( xing.mincol>=0 ) {
	  xing.pctr = &ctr;
	  plane_xing_v[p].push_back(xing);
	}
      }//end of contour loop
    }//end of plane loop

    return;
    
  }//end of function
  
  void FlashEndContourFinder::findCandidateSpacePoints( const std::vector< std::vector<FlashEndContourFinder::CtrXing_t> >& plane_xing_v,
							const int row, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
							std::vector< std::vector<float> >& spacepoints_v ) {
    // we make 2 plane combinations
    for (size_t p1=0; p1<plane_xing_v.size(); p1++) {

      // get contours from this plane
      for ( size_t ictr1=0; ictr1<plane_xing_v[p1].size(); ictr1++ ) {

	const CtrXing_t& xingp1 = plane_xing_v[p1][ictr1];

	// match against contour on other plane
	for ( size_t p2=p1+1; p2<plane_xing_v.size(); p2++) {
	  
	  for ( size_t ictr2=0; ictr2<plane_xing_v[p2].size(); ictr2++ ) {

	    const CtrXing_t& xingp2 = plane_xing_v[p2][ictr2];

	    // ok, we try to find a good 3d point now

	    std::vector<float> tyz;
	    bool pointfound = find3Dpoint( p1, xingp1, p2, xingp2, row, img_v, badch_v, tyz );

	    if ( pointfound )
	      spacepoints_v.push_back( tyz );
	    
	  }//end of p2 ctr loop
	}//end of p2 loop
      }//end of p1 ctr loop
    }//end of p1 loop

    return;
  }

  bool FlashEndContourFinder::find3Dpoint( const int plane1, const FlashEndContourFinder::CtrXing_t& xingp1,
					   const int plane2, const FlashEndContourFinder::CtrXing_t& xingp2,
					   const int row, const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
					   std::vector<float>& pos3d ) {

    // inputs
    // -------
    
    // check if planes make sense
    int plane3 = -1;
    if ( plane1==0 && plane1==1 ) plane3==2;
    else if (plane1==0 && plane1==2) plane3==1;
    else if (plane1==1 && plane2==2) plane3==0;
    else {
      throw std::runtime_error("find3dpoint: unrecognized plane combination");
    }

    pos3d.resize(3,-1);

    // find position in plane3
    plane3 = -1;
    int col3   = -1;
    std::vector<float> intersectionzy;
    int crosses = 0;
    larcv::UBWireTool::getMissingWireAndPlane( plane1, xingp1.maxqcol, plane2, xingp2.maxqcol, plane3, col3, intersectionzy, crosses );

    // no valid intersection based on position
    if ( crosses==0 ) {
      return false;
    }

    // get image meta of plane3
    const larcv::ImageMeta& meta = img_v[plane3].meta();

    // load wire vector
    std::vector<int> wids(3,0);
    wids[plane1] = xingp1.maxqcol;
    wids[plane2] = xingp2.maxqcol;
    wids[plane3] = col3;
    

    // check other wire plane (check a neighborhood)
    bool found_anymatch    = false;
    bool found_chargematch = false;
    std::vector<float> finalzy;
    double min_triarea = -1;
    int min_crosses = 0;
    for (int dc=-3; dc<=3; dc++) {
      int col = col3+dc;
      if ( col<0 || col>=(int)meta.cols() )
	continue;
      
      float pixval = img_v[plane3].pixel( row, col );
      float badval = badch_v[plane3].pixel( row, col );
      bool foundmatch = false;
      if ( pixval>fConfig.pixel_value_threshold[plane3] ) {
	found_anymatch = true;
	found_chargematch = true;
	foundmatch = true;
      }
      else if ( badval > 0 ) {
	found_anymatch = true;
	foundmatch = true;
      }

      if ( foundmatch ) {
	// if match found for this location, we examine the intersection point
	wids[plane3] = col;
	std::vector<float> zy;
	int thiscross = 0;
	double thistri = 0.;
	larcv::UBWireTool::wireIntersection( wids, zy, thistri, thiscross );
	if ( thiscross==1 && ( min_triarea<0 || min_triarea>thistri ) ) {
	  // this point must cross in valid position. if it does, is it also the best intersection found so far (using crossing tri-area)
	  min_crosses = 1;
	  min_triarea = thistri;
	  finalzy = zy;
	}
      }
      
      if ( found_chargematch && min_triarea<0.1) {
	// we find a charge spot AND its pretty good?
	// then we break early
	break;
      }
    }//end of col scan

    // check if we didn't find a good point
    if ( !found_anymatch || min_crosses==0 )
      return false;

    // we did.
    // so set position. we return t,y,z
    float tick = meta.pos_y(row);
    pos3d[0] = tick;
    pos3d[1] = finalzy[1];
    pos3d[2] = finalzy[0];

    // return success
    return true;
  }
  
}
