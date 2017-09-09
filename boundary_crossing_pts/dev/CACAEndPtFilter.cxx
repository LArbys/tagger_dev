#include "CACAEndPtFilter.h"

#include <sstream>
#include <assert.h>

// larlite
#include "LArUtil/LArProperties.h"

// larlitecv
#include "TaggerTypes/dwall.h"

namespace larlitecv {

  
  
  bool CACAEndPtFilter::isEndPointGood( const larlitecv::BoundarySpacePoint& pt, const larlite::opflash* associated_flash,
					const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
					const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,					
					const int endpt_type, const float dist_from_wall, const float chi2_threshold, const float max_dtick ) {

    // evaluates end point goodness
    
    clearClusters();
    
    std::cout << __FILE__ << ":" << __LINE__ << " ----------------------------------" << std::endl;
    
    const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;    
    larlitecv::ContourAStarCluster cluster = m_caca.makeCluster( pt.pos(), img_v, badch_v, plane_contours_v, -10, 3 );
    if ( cluster.m_path3d.size()==0 )
      return false;
    
    const std::vector<float>& start = cluster.m_path3d.front();
    const std::vector<float>& end   = cluster.m_path3d.back();
    float tick_start = start[0]/cm_per_tick + 3200.0;
    float tick_end   = end[0]/cm_per_tick   + 3200.0;

    m_last_clusters.emplace_back( std::move(cluster) );
    
    if ( pt.type()==larlitecv::kAnode || pt.type()==larlitecv::kCathode ) {
      // perform end test, needs to be in-time with flash and be exiting or entering the right direction
      std::cout << "ticks start=" << tick_start << " end=" << tick_end << std::endl;
      float flash_tick = 3200 + associated_flash->Time()/0.5;
      if ( pt.type()==larlitecv::kCathode )
	flash_tick += (256.0)/cm_per_tick;
      std::cout << "flash tick=" << flash_tick << std::endl;

      float dir[3] = { 0, 0, 0};
      std::vector<float> flashend;
      float trackend_tick;
      float norm = 0.;      
      if ( fabs(flash_tick-tick_start) < fabs(flash_tick-tick_end)  ) {
	// start end is closest
	flashend = start;
	trackend_tick = tick_start;
	for (int i=0; i<3; i++) {
	  dir[i] = (end[i]-start[i]);
	  norm += dir[i]*dir[i];
	}
      }
      else {
	// end end if closest
	flashend = end;
	trackend_tick = tick_end;
	for (int i=0; i<3; i++) {
	  dir[i] = (start[i]-end[i]);
	  norm += dir[i]*dir[i];
	}
      }
      norm = sqrt(norm);
      if ( norm>0 )
	for (int i=0; i<3; i++) dir[i] /= norm;

      float dtick = fabs(trackend_tick-flash_tick);
      std::cout << "dtick=" << dtick << std::endl;
      std::cout << "type=" << pt.type() << std::endl;
      std::cout << "dir[0]=" << dir[0] << std::endl;
      
      if ( dtick < max_dtick ) {
	if ( pt.type()==larlitecv::kAnode && dir[0]>0 )
	  return true;
	if ( pt.type()==larlitecv::kCathode && dir[0]<0 )
	  return true;
      }
      
      return false;
    }
    else if (pt.type()<=larlitecv::kDownstream ) {
      // check it is crossing the correct boundary and goes the right direction
      int target_wall = pt.type();
      float fdwall_target_start = 1.0e5;
      float fdwall_target_end   = 1.0e5;      
      fdwall_target_start = larlitecv::dspecificwall( start, target_wall );
      fdwall_target_end   = larlitecv::dspecificwall( end,   target_wall );      

      float dir[3] = { 0, 0, 0};
      float norm = 0.;
      float dwall_closest = 0.;
      if ( fdwall_target_start < fdwall_target_end ) {
	for (int i=0; i<3; i++) {
	  dir[i]  = end[i] - start[i];
	  norm    += dir[i]*dir[i];
	}
	dwall_closest = fdwall_target_start;
      }
      else {
	for (int i=0; i<3; i++) {
	  dir[i]  = start[i] - end[i];
	  norm    += dir[i]*dir[i];
	}
	dwall_closest = fdwall_target_end;	
      }
	
      norm = sqrt(norm);
      if ( norm>0 ) {
	for (int i=0; i<3; i++)
	  dir[i] /= norm;
      }

      if ( dwall_closest<17.0 ) {
	if ( target_wall==larlitecv::kTop && dir[1]<0 )
	  return true;
	else if ( target_wall==larlitecv::kBottom && dir[1]>0 )
	  return true;
	else if ( target_wall==larlitecv::kUpstream && dir[2]>0 )
	  return true;
	else if ( target_wall==larlitecv::kDownstream && dir[2]<0 )
	  return true;
      }
    }
    else if ( pt.type()==larlitecv::kImageEnd ) {
      return true;
    }
    else {
      std::stringstream ss;
      ss << __FILE__ << ":" << __LINE__ << " unrecognized boundary type = " << pt.type() << std::endl;
      throw std::runtime_error( ss.str() );
    }
    
    return false;
  }

  void CACAEndPtFilter::evaluateEndPoints( const std::vector< const std::vector<larlitecv::BoundarySpacePoint>* >& sp_v, const std::vector< const larlite::event_opflash* >& flash_v,
					   const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
					   const std::vector< std::vector<ContourShapeMeta> >& plane_contours_v,
					   const float dist_from_wall, const float chi2_threshold, const float max_dtick,
					   std::vector< std::vector<int> >& passes_filter ) {
    // Function to evaluate a vector of space points
    // Inputs
    // ------
    // sp_v: a vector of vector of points
    // flash_v: a vector of event_opflash containers
    // img_v: vector of image, each element belongs to different planes
    // badch_v: vector badch-tagged images, each element belongs to different planes
    // plane_contours_v: vector of contours for each plane. made by BMTCV class.
    // dist_from_wall: max dist from wall that good end point can be
    // chi2_threshold: how well it matches flashes (not used for now)
    // max_dtick: distance of ends from paired flash

    // Implicit parameters
    // --------------------
    // fTruthInfoLoaded: if true, then will evaluate MC info
    // m_truthinfo_ptr_v: analysis of truth crossing points from larlitecv/app/MCTruthTools/crossingPointsAnaMethods.h/cxx
    // m_recoinfo_ptr_v:  analysis of reco crossing points from larlitecv/app/MCTruthTools/crossingPointsAnaMethods.h/cxx

    // Outputs
    // -------
    // passes_filter: vector of output, 0=fails, 1=passes, follows structure of sp_v    

    // tracking variables
    std::vector<int> good_npasses_caca( larlitecv::kNumEndTypes, 0 );
    std::vector<int> good_nfails_caca(  larlitecv::kNumEndTypes, 0 );
    std::vector<int> bad_npasses_caca(  larlitecv::kNumEndTypes, 0 );
    std::vector<int> bad_nfails_caca(   larlitecv::kNumEndTypes, 0 );
    
    int numreco_good_in_contour = 0;
    int numreco_bad_in_contour = 0;    
    int numreco_good_tot = 0;
    int numreco_bad_tot  = 0;    
    int ireco = -1;
    int numreco_good[7][2] = {0};
    int numreco_bad[7][2]  = {0};

    if ( m_verbosity>0 ) {
      std::cout << "==========================================================================" << std::endl;
      std::cout << "CACAENDPTFILTER: evaluate points" << std::endl;
    }

    passes_filter.clear();
    
    for ( auto const& p_sp_v : sp_v ) {
      std::vector<int> passes_v(p_sp_v->size(),0);

      int isp = -1;
      for (auto const& sp : *p_sp_v ) {
	ireco++; // counter for all spacepoint indices
	isp++;   // counter for sp index of this vector

	if ( sp.type()==larlitecv::kAnode || sp.type()==larlitecv::kCathode ) {

	  int tot_flashidx = sp.getFlashIndex();
	  int loc_flashidx = tot_flashidx;
	  if ( m_verbosity>1 ) {
	    std::cout << "--------------------------------------------------------------" << std::endl;
	    std::cout << "[ipt " << ireco << "] Anode/Cathode space point" << std::endl;
	    std::cout << "  flash index: " << tot_flashidx << std::endl;
	  }
	  
	  const larlite::opflash* popflash = NULL;
	  for ( int ievop=0; ievop<(int)flash_v.size(); ievop++ ) {
	    if ( loc_flashidx >= flash_v.at(ievop)->size() )
	      loc_flashidx -= (int)flash_v.at(ievop)->size();
	    else
	      popflash = &(flash_v.at(ievop)->at(loc_flashidx));
	  }

	  if ( m_verbosity>1 ) {
	    std::cout << "  flash local index: " << loc_flashidx << std::endl;		
	    std::cout << "  flash pointer: " << popflash << std::endl;
	  }
	  bool passes = isEndPointGood( sp, popflash, img_v, badch_v, plane_contours_v, sp.type(), dist_from_wall, chi2_threshold, max_dtick );
	  if ( m_verbosity>0 ) {
	    std::cout << "Result: " << passes << std::endl;
	  }

	  const larlitecv::RecoCrossingPointAna_t* recoinfo   = NULL;
	  const larlitecv::TruthCrossingPointAna_t* truthinfo = NULL;
	  
	  if ( fTruthInfoLoaded ) {

	    try {
	      recoinfo = &(m_recoinfo_ptr_v->at(ireco));
	    }
	    catch (...) {
	      continue;
	    }

	    if ( m_verbosity>0 ) 
	      std::cout << "Has Truth Match: " << recoinfo->truthmatch << std::endl;	    
	    
	    if ( recoinfo->truthmatch==1 ) {
	    
	      const larlitecv::TruthCrossingPointAna_t* truthinfo = NULL;
	      try {
		truthinfo = &(m_truthinfo_ptr_v->at(recoinfo->truthmatch_index));
	      }
	      catch (...) {
		continue;
	      }
	      
	      if ( m_verbosity>1 )  {
		std::cout << "Truth crossing position: "
			  << "(" << truthinfo->crossingpt_detsce[0] << "," << truthinfo->crossingpt_detsce[1] << "," << truthinfo->crossingpt_detsce[2] << ")"
			  << std::endl;
	      }
	    
	      // good reco point
	      if ( passes ) {
		passes_v[isp] = 1;	  
		good_npasses_caca[(int)sp.type()-(int)larlitecv::kAnode]++;
	      }
	      else {
		good_nfails_caca[(int)sp.type()-(int)larlitecv::kAnode]++;
	      }
	    }//end if reco is truth matched
	    else {
	      // bad reco point
	      if ( passes ) {
		passes_v[isp] = 1;	  
		bad_npasses_caca[(int)sp.type()-(int)larlitecv::kAnode]++;
	      }
	      else
		bad_nfails_caca[(int)sp.type()-(int)larlitecv::kAnode]++;
	    }//end of if non-matched reco point
	  }// if truth loaded
	}//if anode/cathode
      }//end of space points loop
      
      passes_filter.emplace_back( std::move(passes_v) );
      
    }//end of spacepoints_vv loop
    
    if ( fTruthInfoLoaded && m_verbosity>0 ) {
      std::cout << "CACA Summary" << std::endl;
      for (int i=0; i<larlitecv::kNumEndTypes; i++)
	std::cout << "  Good " << larlitecv::BoundaryEndNames((larlitecv::BoundaryEnd_t)i) << ":   " << good_npasses_caca[i] << "/" << good_npasses_caca[i]+good_nfails_caca[i] << std::endl;
      for (int i=0; i<larlitecv::kNumEndTypes; i++)
	std::cout << "  Bad "  << larlitecv::BoundaryEndNames((larlitecv::BoundaryEnd_t)i) << ":   " << bad_npasses_caca[i] <<  "/" << bad_npasses_caca[i]+bad_nfails_caca[i] << std::endl;
    }
    
  }
  
  void CACAEndPtFilter::setTruthInformation( const std::vector<larlitecv::TruthCrossingPointAna_t>& truthinfo, const std::vector<larlitecv::RecoCrossingPointAna_t>& recoinfo ) {
    m_truthinfo_ptr_v = &truthinfo;
    m_recoinfo_ptr_v  = &recoinfo;
    fTruthInfoLoaded = true;
  }

}
