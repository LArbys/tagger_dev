#include <iostream>

#include <string>

// ROOT
#include "TFile.h"
#include "TTree.h"

// larcv
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
//#include "DataFormat/EventROI.h"
#include "DataFormat/EventChStatus.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
//#include "UBWireTool/UBWireTool.h"

// larlite
// #include "DataFormat/mctruth.h"
// #include "DataFormat/mcpart.h"
// #include "DataFormat/mctrajectory.h"
// #include "DataFormat/mctrack.h"
// #include "DataFormat/mcshower.h"
// #include "DataFormat/simch.h"
#include "DataFormat/trigger.h"
// #include "DataFormat/track.h"
// #include "LArUtil/LArProperties.h"
// #include "LArUtil/Geometry.h"

// larlitecv
#include "Base/DataCoordinator.h"
//#include "SCE/SpaceChargeMicroBooNE.h"
#include "GapChs/EmptyChannelAlgo.h"
#include "TaggerTypes/BoundaryMuonTaggerTypes.h"
#include "TaggerTypes/BoundaryEndPt.h"
#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerTypes/BMTrackCluster3D.h"
#include "TaggerTypes/Path2Pixels.h"
#include "ThruMu/ThruMuTracker.h"
#include "extractTruthMethods.h"
#include "crossingPointsAnaMethods.h"


// OpenCV
#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "CVUtil/CVUtil.h"
#endif


// larlitecv
#include "Base/DataCoordinator.h"
#include "ChargeSegmentAlgos/PathPixelChargeMethods.h"

int main( int nargs, char** argv ) {
  
  std::cout << "COMPARE TRACK CHARGE" << std::endl;

  // What we're doing.

  // 1) Load information using DataCoordinator
  // 2) Load truth end points for cosmics
  // 3) use the end points to reco the tracks (either linearcharge or astar tracker)
  // 4) use the reco track points to calculate the charge and charge profile
  // 5) make plots comparing the relative values on the planes -- is the distribution good enough?

  // 6) make reco end points. repeat exercise.  But we can now compare charge matching for good and bad tracks


  if ( nargs!=2 ) {
    std::cout << "usage: ./compare_trackq [config file]" << std::endl;
    return 0;
  }
  std::string cfg_file = argv[1];

  larcv::PSet cfg_root = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg_root.get<larcv::PSet>("CompareTrackQ");

  enum SourceTypes_t { kSource=0, kCROIfile, kNumSourceTypes };
  std::string source_param[2] = { "InputSourceFilelist", "InputCROIFilelist" };
  larlitecv::DataCoordinator dataco[kNumSourceTypes];
  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    std::cout << "LOADING " << source_param[isrc] << " FILES" << std::endl;
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArCV"),   "larcv" );
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArLite"), "larlite" );
    dataco[isrc].configure( cfg_file, "StorageManager", "IOManager", "CompareTrackQ" );
    dataco[isrc].initialize();
  }

  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    std::cout << "data[" << source_param[isrc] << "] entries=" << dataco[isrc].get_nentries("larcv") << std::endl;
  }


  // =====================================================================
  // configuration parameters

  std::string outfname   = pset.get<std::string>("OutputAnaFile");
  int verbosity          = pset.get<int>("Verbosity",0);
  bool ismc              = pset.get<bool>("IsMC");
  std::string inputimgs  = pset.get<std::string>("InputLArCVImages");
  
  bool badch_in_file     = pset.get<bool>("BadChImageInFile");
  std::string trigname   = pset.get<std::string>("TriggerProducerName");
  std::vector<std::string> flashprod  = pset.get<std::vector<std::string> >("OpFlashProducer");  
  bool printFlashEnds    = pset.get<bool>("PrintFlashEnds");
  bool printImages       = pset.get<bool>("PrintImages");
  float fthreshold       = pset.get<float>("PixelThreshold");
  std::vector<float> thresholds_v( 3, fthreshold );
  std::vector<float> label_thresholds_v( 3, -10 );  
  std::vector<int> label_neighborhood(3,0);
  
  // bool use_reclustered   = pset.get<bool>("UseReclustered");
  // int tag_neighborhood   = pset.get<int>("TagNeighborhood",10);
  // int fvertex_radius     = pset.get<int>("PixelRadius");  


 // =====================================================================

  // // setup input
  // int kThruMu, kStopMu, kUntagged, kCROI, kNumStages;

  // std::vector<std::string> stages_pixel_producers;
  // std::vector<std::string> stages_track_producers;  
  // if ( use_reclustered ) {
  //   kNumStages = 4;
  //   kThruMu = 0;
  //   kStopMu = 1;    
  //   kUntagged = 2;
  //   kCROI = 3;
  //   stages_pixel_producers.resize(kNumStages);
  //   stages_pixel_producers[0] = "mergedthrumupixels";
  //   stages_pixel_producers[1] = "mergedstopmupixels";
  //   stages_pixel_producers[2] = "mergeduntaggedpixels";
  //   stages_pixel_producers[3] = "croipixels";
  //   stages_track_producers.resize(kNumStages);    
  //   stages_track_producers[0] = "mergedthrumu3d";
  //   stages_track_producers[1] = "mergedstopmu3d";
  //   stages_track_producers[2] = "mergeduntagged3d";
  //   stages_track_producers[3] = "croi3d";    
  // }
  
  // setup output
  TFile* rfile = new TFile(outfname.c_str(), "recreate");
  TTree* tree = new TTree("compareq", "Compare Track Charge");

  // Event Index
  int run, subrun, event;

  // Truth Quantities about interaction and lepton
  larlitecv::TruthData_t truthdata;

  // // Truth Pixel Quantities
  // int nthreshold_pixels[4]; // number of pixels above threshold
  // int ncosmic_pixels[4];    // number of non-neutrino pixels
  // int nnu_pixels[4];        // number of neutrino pixels
  // int nvertex_pixels[4];    // number of neutrino pixels within some pixel radius of vertex
  // int nvertex_badch[4];

  // // Tagged Pixel Quantities
  // int ncosmic_tagged[kNumStages][4];       // number of non-neutrino pixels tagged
  // int ncosmic_tagged_once[kNumStages][4];  // number of non-neutrino pixels tagged
  // int ncosmic_tagged_many[kNumStages][4];  // number of non-neutrino pixels tagged
  // int nnu_tagged[kNumStages][4];           // number of neutrino pixels tagged
  // int nvertex_tagged[kNumStages][4];       // number of neutrino pixels within some pixel radius of vertex tagged
  // int nvertex_incroi[4];                   // number of pixels near neutrino vertex that are in an ROI
  // std::stringstream s_arr;
  // s_arr << "[" << (int)kNumStages << "][4]/I";

  // // ROI quantities
  // int num_rois;     // number of identified ROis
  // int nnu_inroi[4]; // number of nu pixels contained in the CROI
  // int vertex_in_croi; // is vertex in an CROI
  // float closest_dist_to_vertex;
  // int closest_dist_stage;

  // Crossing Point data
  larlitecv::CrossingPointAnaData_t xingptdata;

  // Event
  tree->Branch("run",&run,"run/I");
  tree->Branch("subrun",&subrun,"subrun/I");
  tree->Branch("event",&event,"event/I");

  // Truth
  truthdata.bindToTree( tree );

  // tree->Branch("nthreshold_pixels", nthreshold_pixels, "nthreshold_pixels[4]/I");
  // tree->Branch("ncosmic_pixels",    ncosmic_pixels,    "ncosmic_pixels[4]/I");
  // tree->Branch("nnu_pixels",        nnu_pixels,        "nnu_pixels[4]/I");
  // tree->Branch("nvertex_pixels",    nvertex_pixels,    "nvertex_pixels[4]/I");
  // tree->Branch("nvertex_badch",     nvertex_badch,     "nvertex_badch[4]/I");

  // tree->Branch("ncosmic_tagged",      ncosmic_tagged,      std::string("ncosmic_tagged"+s_arr.str()).c_str() );
  // tree->Branch("ncosmic_tagged_once", ncosmic_tagged_once, std::string("ncosmic_tagged_once"+s_arr.str()).c_str() );
  // tree->Branch("ncosmic_tagged_many", ncosmic_tagged_many, std::string("ncosmic_tagged_many"+s_arr.str()).c_str() );
  // tree->Branch("nnu_tagged",          nnu_tagged,          std::string("nnu_tagged"+s_arr.str()).c_str() );
  // tree->Branch("nvertex_tagged",      nvertex_tagged,      std::string("nvertex_tagged"+s_arr.str()).c_str() );
  // tree->Branch("nvertex_incroi",      nvertex_incroi,      "nvertex_incroi[4]/I" );

  // track charge ana variables
  float plane_trackq[3] = {0.0};
  float plane_goodfrac[3] = {0.0};
  int plane_goodpix[3] = {0};
  float uv_trackqdiff = 0.;
  float uy_trackqdiff = 0.;
  float vy_trackqdiff = 0.;
  TTree* ptqtree = new TTree("PointQ", "Pixel Charge per reco track point");
  ptqtree->Branch("plane_goodpix", plane_goodpix, "plane_goodpix[3]/I");
  ptqtree->Branch("plane_goodfrac", plane_goodfrac, "plane_goodfrac[3]/F");
  ptqtree->Branch("plane_trackq", plane_trackq, "plane_trackq[3]/F");
  ptqtree->Branch("uv_trackqdiff", &uv_trackqdiff, "uv_trackqdiff/F");
  ptqtree->Branch("uy_trackqdiff", &uy_trackqdiff, "uy_trackqdiff/F");
  ptqtree->Branch("vy_trackqdiff", &vy_trackqdiff, "vy_trackqdiff/F");

  float plane_tracktotq[3] = {0.0};
  float uv_totqdiff = 0.;
  float uy_totqdiff = 0.;
  float vy_totqdiff = 0.;    
  TTree* trackqtree = new TTree("TrackQ", "Pixel Charger per track");
  trackqtree->Branch( "plane_tracktotq", plane_tracktotq, "plane_tracktotq[3]" ); 
  trackqtree->Branch("uv_totqdiff", &uv_totqdiff, "uv_totqdiff/F");
  trackqtree->Branch("uy_totqdiff", &uy_totqdiff, "uy_totqdiff/F");
  trackqtree->Branch("vy_totqdiff", &vy_totqdiff, "vy_totqdiff/F");
 
  xingptdata.bindToTree( tree );
  
  // // Space Charge Corrections
  // larlitecv::SpaceChargeMicroBooNE sce;
  /// ---------------------------------------------------------------------------------------
  // ALGO SETUP
  
  // Empty Channel Algo
  larlitecv::EmptyChannelAlgo emptyalgo;

  // ThruMu Tracker
  larcv::PSet tracker_pset = pset.get<larcv::PSet>("ThruMuTracker");
  larlitecv::ThruMuTrackerConfig tracker_cfg = larlitecv::ThruMuTrackerConfig::MakeFromPSet( tracker_pset );
  larlitecv::ThruMuTracker thrumualgo( tracker_cfg );

  // -----------------------------------------------------------------------------------------

  int nentries = dataco[kCROIfile].get_nentries("larcv");

  int user_nentries =   pset.get<int>("NumEntries",-1);
  int user_startentry = pset.get<int>("StartEntry",-1);
  int startentry = 0;
  if ( user_startentry>=0 ) {
    startentry = user_startentry;
  }
  int endentry = nentries;
  if ( user_nentries>=0 ) {
    if ( user_nentries+startentry<nentries )
      endentry = user_nentries+startentry;
  }
  if ( startentry>=nentries ) {
    std::cout << "Starting beyond end of file. Nothing to do." << std::endl;
    return 0;
  }
  std::cout << "Start Entry: " << startentry << std::endl;
  std::cout << "End Entry: " << endentry-1 << std::endl;
  std::cout << "Buckle up!" << std::endl;
  

  for (int ientry=startentry; ientry<endentry; ientry++) {

    // load the entry
    dataco[kCROIfile].goto_entry(ientry,"larcv");
    dataco[kCROIfile].get_id(run,subrun,event);

    // sync up the other files
    for ( int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
      if ( isrc!=kCROIfile ) {
        dataco[isrc].goto_event(run,subrun,event, "larcv");
      }
    }

    if ( ientry%10==0 || verbosity>0 ) {
      std::cout << "entry " << ientry << ": (r,s,e)=(" << run << ", " << subrun << ", " << event << ")" << std::endl;
    }

    // initialize the output variables
    truthdata.clear();
    xingptdata.clear();
    // for (int p=0; p<4; p++) {
    //   ncosmic_pixels[p] = 0;
    //   nnu_pixels[0] = 0;
    //   nvertex_pixels[p] = 0;
    //   nthreshold_pixels[p] = 0;
    //   nnu_inroi[p] = 0;
    //   nvertex_badch[p] = 0;
    //   for (int istage=0; istage<kNumStages; istage++) {
    // 	ncosmic_tagged[istage][p] = 0;
    // 	ncosmic_tagged_once[istage][p] = 0;
    // 	ncosmic_tagged_many[istage][p] = 0;
    // 	nnu_tagged[istage][p] = 0;
    // 	nvertex_tagged[istage][p] = 0;
    //   }
    //   nvertex_incroi[p] = 0;
    // }
    // vertex_in_croi = 0;    


    // ok now to do damage

    // get the original, segmentation, and tagged images
    larcv::EventImage2D* ev_imgs   = (larcv::EventImage2D*)dataco[kSource].get_larcv_data(larcv::kProductImage2D,inputimgs);
    larcv::EventImage2D* ev_badch  = NULL;
    larcv::EventImage2D* ev_segs   = NULL;
    if ( ismc ) {
      ev_segs = (larcv::EventImage2D*)dataco[kSource].get_larcv_data(larcv::kProductImage2D,"segment");
    }
    if ( badch_in_file ) {
      // we either get the badch from the file
      ev_badch = (larcv::EventImage2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductImage2D,"gapchs");
    }
    else {
      // or we have to make the badch image from a ChStatus object
      larcv::EventChStatus* ev_status = (larcv::EventChStatus*)dataco[kSource].get_larcv_data( larcv::kProductChStatus, "tpc" );
      std::vector<larcv::Image2D> chstatus_img_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
      ev_badch = new larcv::EventImage2D;
      ev_badch->Emplace( std::move(chstatus_img_v) );
    }
    
    // // get the output of the tagger
    // larcv::EventPixel2D* ev_pix[kNumStages] = {0};
    // larlite::event_track* ev_track[kNumStages] = {0};
    // try {
    //   ev_pix[kThruMu]    = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,stages_pixel_producers[kThruMu]);
    //   ev_track[kThruMu]  = (larlite::event_track*)dataco[kCROIfile].get_larlite_data(larlite::data::kTrack,stages_track_producers[kThruMu]);      
    // }
    // catch (...) {
    //   ev_pix[kThruMu] = NULL;
    //   ev_track[kThruMu] = NULL;
    // }
    // try {
    //   ev_pix[kStopMu]    = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,stages_pixel_producers[kStopMu]);
    //   ev_track[kStopMu]  = (larlite::event_track*)dataco[kCROIfile].get_larlite_data(larlite::data::kTrack,stages_track_producers[kStopMu]);            
    // }
    // catch (...) {
    //   ev_pix[kStopMu] = NULL;
    // }
    // try {
    //   ev_pix[kUntagged]  = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,stages_pixel_producers[kUntagged]);
    //   ev_track[kUntagged]  = (larlite::event_track*)dataco[kCROIfile].get_larlite_data(larlite::data::kTrack,stages_track_producers[kUntagged]);                  
    // }
    // catch (...) {
    //   ev_pix[kUntagged] = NULL;
    // }
    // try {
    //   ev_pix[kCROI]      = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,stages_pixel_producers[kCROI]);
    //   ev_track[kCROI]  = (larlite::event_track*)dataco[kCROIfile].get_larlite_data(larlite::data::kTrack,stages_track_producers[kCROI]);                  
    // }
    // catch (...) {
    //   ev_pix[kCROI] = NULL;
    // }

    const std::vector<larcv::Image2D>& imgs_v   = ev_imgs->Image2DArray();
    const std::vector<larcv::Image2D>& badch_v  = ev_badch->Image2DArray();
    const std::vector<larcv::Image2D>* segs_v   = NULL;
    if ( ismc ) {
      segs_v = &(ev_segs->Image2DArray());
    }

#ifdef USE_OPENCV
    std::vector<cv::Mat> cvimgs_v;    
    if ( printImages ) {
      for ( auto const& img : imgs_v ) {
	cv::Mat cvimg = larcv::as_mat_greyscale2bgr( img, fthreshold, 100.0 );
	cvimgs_v.emplace_back(std::move(cvimg));
      }
    }
#endif
    
    // // get the result of the contained ROI analysis
    // larcv::EventROI* ev_contained_roi = NULL;
    // try {
    //   ev_contained_roi = (larcv::EventROI*)dataco[kCROIfile].get_larcv_data(larcv::kProductROI,"croi");
    // }
    // catch (...) {
    //   ev_contained_roi = NULL;
    // }
    // std::vector<larcv::ROI> containedrois_v;
    // if ( ev_contained_roi ) {
    //   containedrois_v = ev_contained_roi->ROIArray();
    //   num_rois = (int)containedrois_v.size();
    //   std::cout << "====ROIs===========================" << std::endl;
    //   for ( auto& roi : containedrois_v ) {
    // 	std::cout << " roi: " << roi.dump();
    //   }
    //   std::cout << "===================================" << std::endl;
    // }
      
    // get the boundary end point info (only if have MC info to compare against)
    std::vector<larcv::EventPixel2D*> ev_spacepoints(7,0);
    std::string spacepoint_producers[7] = { "topspacepts", "botspacepts", "upspacepts", "downspacepts", "anodepts", "cathodepts", "imgendpts" };
    for ( int i=0; i<7; i++ ) {
      try {
    	ev_spacepoints[i] = (larcv::EventPixel2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductPixel2D,spacepoint_producers[i]);
      }
      catch (...) {
    	ev_spacepoints[i] = NULL;
      }
      if ( ev_spacepoints[i]!=NULL )
    	std::cout << "number of " << spacepoint_producers[i] << ": " << ev_spacepoints[i]->Pixel2DArray(0).size() << std::endl;	
    }

    // get the opflashes
    std::vector< larlite::event_opflash* > opflash_v;
    for ( auto const& prodname : flashprod ) {
      larlite::event_opflash* ev_flash = (larlite::event_opflash*)dataco[kSource].get_larlite_data(larlite::data::kOpFlash, prodname);
      std::cout << "number of flashes in " << prodname << ": " << ev_flash->size() << std::endl;
      opflash_v.push_back( ev_flash );
    }

    // get other information, e.g. truth
    larlite::event_mctruth* ev_mctruth   = NULL;
    larlite::event_mctrack* ev_mctrack   = NULL;
    larlite::event_mcshower* ev_mcshower = NULL; 
    larlite::trigger* ev_trigger         = (larlite::trigger*) dataco[kSource].get_larlite_data(larlite::data::kTrigger,trigname);
    if ( ismc ) {
      ev_mctruth = (larlite::event_mctruth*) dataco[kSource].get_larlite_data(larlite::data::kMCTruth,"generator");
      ev_mctrack = (larlite::event_mctrack*) dataco[kSource].get_larlite_data(larlite::data::kMCTrack,"mcreco");
      ev_mcshower = (larlite::event_mcshower*)dataco[kSource].get_larlite_data(larlite::data::kMCShower,"mcreco");

      // extract the truth quantities of interest
      larlitecv::extractTruth( *ev_mctruth, *ev_mctrack, truthdata );
    }

    // =======================================================================
    // MC INFO ANALYSIS

    // quantities filled if MC present
    std::vector<larcv::Image2D> nupix_imgs_v;
    // std::vector< std::vector<int> > start_pixels;
    // std::vector< std::vector<float> > start_crossingpts;
    // std::vector< std::vector<int> > end_pixels;
    // std::vector< std::vector<float> > end_crossingpts;
    std::vector<int> vertex_col(3,-1);
    std::vector<double> vtx_sce(3,0);    
    int vertex_row = -1;

    if ( ismc ) {

      // loop over MC tracks, get end points of muons
      larlitecv::analyzeCrossingMCTracks( xingptdata, imgs_v.front().meta(),  ev_trigger, ev_mctrack, opflash_v, printFlashEnds );
      int intime_cosmics = xingptdata.true_intime_thrumu + xingptdata.true_intime_stopmu;
      std::cout << "End points from MC Truth" << std::endl;
      for (int i=0; i<larlitecv::kNumEndTypes; i++) {
	std::cout << "  " << spacepoint_producers[i] << ": " << xingptdata.true_crossingpoints[i] << std::endl;
      }

      // loop over MC tracks. Get the truth end points, run thrumu tagger!
      std::cout << ev_mctrack->size() << " = " << xingptdata.mctrack_imgendpoint_indices.size() << std::endl;
      std::vector< larlitecv::BMTrackCluster3D > mc_recotracks;
      for (int itrack=0; itrack<(int)ev_mctrack->size(); itrack++) {
	// did this track have end points?
	int nendpts = xingptdata.mctrack_imgendpoint_indices.at(itrack).size();
	std::cout << "[mc track #" << itrack << "] num end points=" << nendpts << " ";
	if (nendpts==0) {
	  std::cout << std::endl;
	  continue;
	}
	std::vector< const larlitecv::BoundarySpacePoint* > mctrack_endpts;
	int start_index = -1;
	int end_index = -1;
	if ( nendpts>0 )
	  start_index = xingptdata.mctrack_imgendpoint_indices[itrack][0];
	if ( nendpts>1 )
	  end_index = xingptdata.mctrack_imgendpoint_indices[itrack][1];
	larlitecv::BoundarySpacePoint* pstartpt = NULL;
	larlitecv::BoundarySpacePoint* pendpt   = NULL;
	if (  start_index!=-1 ) {
	  // make a boundary space point object for the start;
	  std::vector< larlitecv::BoundaryEndPt > planepts;
	  for (int p=0; p<3; p++) {
	    int row = xingptdata.start_pixels[start_index][0];
	    int col = xingptdata.start_pixels[start_index][p+1];
	    if ( col<0 ) col = 0; // hack
	    larlitecv::BoundaryEndPt planept( row, col, (larlitecv::BoundaryEnd_t)xingptdata.start_type[start_index] );
	    planepts.emplace_back( std::move(planept) );
	  }
	  pstartpt = new larlitecv::BoundarySpacePoint( (larlitecv::BoundaryEnd_t)xingptdata.start_type[start_index], std::move(planepts), imgs_v.front().meta() );
	  mctrack_endpts.push_back( pstartpt );
	}
	if ( end_index!=-1 ) {
	  // make a boundary space point object for the end;
	  std::vector< larlitecv::BoundaryEndPt > planepts;
	  for (int p=0; p<3; p++) {
	    int row=xingptdata.end_pixels[end_index][0];
	    int col=xingptdata.end_pixels[end_index][p+1];
	    if ( col<0 ) col = 0;
	    larlitecv::BoundaryEndPt planept( row, col, (larlitecv::BoundaryEnd_t)xingptdata.end_type[end_index] );
	    planepts.emplace_back( std::move(planept) );
	  }
	  pendpt = new larlitecv::BoundarySpacePoint( (larlitecv::BoundaryEnd_t)xingptdata.end_type[end_index], std::move(planepts), imgs_v.front().meta() );
	  mctrack_endpts.push_back( pendpt );
	}
	std::cout << " startidx=" << start_index << " endidx=" << end_index << " ";

	std::vector< larlitecv::BMTrackCluster3D > trackclusters;
	std::vector< larcv::Image2D > tagged_v;
	std::vector<int> used_endpoints_indices;
	for ( auto const& img : imgs_v ) {
	  larcv::Image2D tag( img.meta() );
	  tag.paint(0);
	  tagged_v.emplace_back( std::move(tag) );
	}
	if ( start_index>=0 && end_index>=0 ) {
	  std::cout << " start=(" << mctrack_endpts[0]->at(0).row << "," << mctrack_endpts[0]->at(0).col << "," << mctrack_endpts[0]->at(1).col << "," << mctrack_endpts[0]->at(2).col << ")";
	  std::cout << " end(" << mctrack_endpts[1]->at(0).row << "," << mctrack_endpts[1]->at(0).col << "," << mctrack_endpts[1]->at(1).col << "," << mctrack_endpts[1]->at(2).col << ")";
	  try {
	    thrumualgo.makeTrackClusters3D( imgs_v, badch_v, mctrack_endpts, trackclusters, tagged_v, used_endpoints_indices );
	    if ( trackclusters.size()>0 ) {
	      mc_recotracks.emplace_back( std::move(trackclusters.at(0)) );
	    }
	  }
	  catch ( std::exception& e ) {
	    std::cout << std::endl;
	    std::cout << "  ThruMu Failed: " << e.what() << std::endl;
	  }
	}
	mctrack_endpts.clear();

	std::cout << "tracker return " << trackclusters.size() << " thrumu tracks" << std::endl;
	
	delete pstartpt;
	delete pendpt;
	std::cout << std::endl;
      }//end of mctrack loop

      // --------------------------------------------------------------------------------
      // Compare Q sums between planes
      for ( auto const& bmtrack : mc_recotracks ) {
	std::vector<larlitecv::PixelQPt> trackq = larlitecv::getPixelQPts( bmtrack.path3d, imgs_v, badch_v, 5, 5.0, 0.3 );
	for ( auto const& qpt : trackq )  {
	  for (int p=0; p<3; p++) {
	    plane_trackq[p] = qpt.getPlaneCharge()[p];
	    float totpix = qpt.getGoodPixelsPerPlane()[p]+qpt.getBadPixelsPerPlane()[p];
	    plane_goodpix[p] = qpt.getGoodPixelsPerPlane()[p];
	    if ( totpix>0 )
	      plane_goodfrac[p] = float(plane_goodpix[p])/totpix;
	    else
	      plane_goodfrac[p] = 0.0;
	  }
	  uv_trackqdiff = plane_trackq[0]-plane_trackq[1];
	  uy_trackqdiff = plane_trackq[0]-plane_trackq[2];
	  vy_trackqdiff = plane_trackq[1]-plane_trackq[2];
	  ptqtree->Fill();
	}
	std::vector<double> total_track_q = larlitecv::getTrackTotalPixelCharge( bmtrack.path3d, imgs_v, badch_v, 5, 5.0, 0.3 );
	for ( int p=0; p<3; p++) {
	  plane_tracktotq[p] = total_track_q[p];
	}
	uv_totqdiff = plane_tracktotq[0]-plane_tracktotq[1];
	uy_totqdiff = plane_tracktotq[0]-plane_tracktotq[2];
	vy_totqdiff = plane_tracktotq[1]-plane_tracktotq[2];
	trackqtree->Fill();
      }
      // --------------------------------------------------------------------------------
      
      // --------------------------------------------------------------------------------
      // OPEN CV VISUALIZATION
#ifdef USE_OPENCV
      // draw start points
      for ( int ipt=0; ipt<(int)xingptdata.start_pixels.size(); ipt++) {
	int row=xingptdata.start_pixels[ipt][0];
	if (row<0) row = 0;
	if (row>=(int)imgs_v.front().meta().rows() ) row = (int)imgs_v.front().meta().rows()-1;
	for (int p=0; p<3; p++) {
	  int col = xingptdata.start_pixels[ipt][p+1];
	  if ( col<0 ) col = 0;
	  if ( col>=(int)imgs_v.front().meta().cols() ) col = imgs_v.front().meta().cols()-1;
	  cv::circle( cvimgs_v[p], cv::Point(col,row), 6, cv::Scalar( 0, 0, 255, 255 ), -1 );
	}
      }
      // draw end points
      for ( int ipt=0; ipt<(int)xingptdata.end_pixels.size(); ipt++) {
	int row=xingptdata.end_pixels[ipt][0];
	if (row<0) row = 0;
	if (row>=(int)imgs_v.front().meta().rows() ) row = (int)imgs_v.front().meta().rows()-1;
	for (int p=0; p<3; p++) {
	  int col = xingptdata.end_pixels[ipt][p+1];
	  if ( col<0 ) col = 0;
	  if ( col>=(int)imgs_v.front().meta().cols() ) col = imgs_v.front().meta().cols()-1;
	  cv::circle( cvimgs_v[p], cv::Point(col,row), 6, cv::Scalar( 255, 0, 0, 255 ), -1 );
	}
      }
      // draw tracks      
      for (auto const& track : mc_recotracks ) {
	std::vector<larcv::Pixel2DCluster> plane_pixels = larlitecv::getTrackPixelsFromImages( track.path3d,  imgs_v, badch_v, label_thresholds_v, label_neighborhood, 0.3 );
	int p=-1;
	for ( auto const& pixs : plane_pixels ) {
	  p++;
	  for ( auto const& pix : pixs ) {
	    cv::circle( cvimgs_v[p], cv::Point( pix.X(), pix.Y() ), 1, cv::Scalar( 0, 200, 0, 200 ), -1 );
	  }
	}
      }
      // --------------------------------------------------------------------------------      
#endif
      
    } // end of MC functions   

#ifdef USE_OPENCV

    for (int p=0; p<3; p++) {
      cv::Mat& cvimg = cvimgs_v[p];
      std::stringstream path;
      path << "trackqimgs/tracks_r" << run << "_s" << subrun << "_e" << event << "_p" << p << ".png";
      cv::imwrite( path.str(), cvimg );
    }
    
#endif
    
    tree->Fill();

  }//end of entry loop

  rfile->Write();
  //
  return 0;
}
