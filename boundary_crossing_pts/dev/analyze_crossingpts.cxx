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
#include "UBWireTool/UBWireTool.h"

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
#include "GapChs/EmptyChannelAlgo.h"

#include "TaggerTypes/BoundaryMuonTaggerTypes.h"
#include "TaggerTypes/BoundaryEndPt.h"
#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerTypes/BMTrackCluster3D.h"
#include "TaggerTypes/Path2Pixels.h"

#include "ThruMu/BoundaryMuonTaggerAlgoConfig.h"
#include "ThruMu/BoundaryMuonTaggerAlgo.h"
#include "ThruMu/FlashMuonTaggerAlgoConfig.h"
#include "ThruMu/FlashMuonTaggerAlgo.h"
#include "ThruMu/EndPointFilter.h"
#include "ThruMu/RadialEndpointFilter.h"
#include "ThruMu/PushBoundarySpacePoint.h"

#include "MCTruthTools/extractTruthMethods.h"
#include "MCTruthTools/crossingPointsAnaMethods.h"


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
  float fMatchRadius     = pset.get<float>("MatchRadius");
  bool printFlashEnds    = pset.get<bool>("PrintFlashEnds");
  bool printImages       = pset.get<bool>("PrintImages");
  float fthreshold       = pset.get<float>("PixelThreshold");
  std::vector<float> thresholds_v( 3, fthreshold );
  std::vector<float> label_thresholds_v( 3, -10 );  
  std::vector<int> label_neighborhood(3,0);
  
 // =====================================================================

  // // setup input
  // int kThruMu, kStopMu, kUntagged, kCROI, kNumStages;

  // setup output
  TFile* rfile = new TFile(outfname.c_str(), "recreate");
  TTree* tree = new TTree("xingpteventana", "Compare Track Charge");
  TTree* mcxingpt_tree = new TTree("mcxingptana", "Info on MC Crossing Point");  

  // Event Indexf
  int run, subrun, event;

  // Truth Quantities about interaction and lepton
  larlitecv::TruthData_t truthdata;
  
  // Crossing Point data
  larlitecv::CrossingPointAnaData_t xingptdata;

  // Event
  tree->Branch("run",&run,"run/I");
  tree->Branch("subrun",&subrun,"subrun/I");
  tree->Branch("event",&event,"event/I");

  // Truth
  truthdata.bindToTree( tree );

  // Crossing Pt Data
  xingptdata.bindToTree( tree );
  
  // truth end point track reco metrics
  int ntracks_2planeq = 0;
  int ntracks_recod_2planeq = 0;
  int ntracks_all = 0;
  int ntracks_recod_all = 0;
  tree->Branch( "ntracks_2planeq", &ntracks_2planeq, "ntracks_2planeq/I" );
  tree->Branch( "ntracks_recod_2planeq", &ntracks_recod_2planeq, "ntracks_recod_2planeq/I" );
  tree->Branch( "ntracks_all", &ntracks_all, "ntracks_all/I" );
  tree->Branch( "ntracks_recod_all", &ntracks_recod_all, "ntracks_recod_all/I" );  
   
  xingptdata.bindToTree( tree );

  int mcxingpt_type;
  int mcxingpt_matched;
  int mcxingpt_matched_type;
  int mcxingpt_nplaneswcharge;
  int mcxingpt_wire[3];
  float mcxingpt_dist;
  float mcxingpt_pos[3];
  mcxingpt_tree->Branch( "truth_type", &mcxingpt_type, "truth_type/I" )
;  mcxingpt_tree->Branch( "matched", &mcxingpt_matched, "matched/I" );
  mcxingpt_tree->Branch( "matched_type", &mcxingpt_matched_type, "matched_type/I" );  
  mcxingpt_tree->Branch( "nplaneswcharge", &mcxingpt_nplaneswcharge, "nplaneswcharge/I" );
  mcxingpt_tree->Branch( "wire", mcxingpt_wire, "wire[3]/I" );
  mcxingpt_tree->Branch( "dist", &mcxingpt_dist, "dist/F" );
  mcxingpt_tree->Branch( "pos", mcxingpt_pos, "pos[3]/F" );
  
  /// ---------------------------------------------------------------------------------------
  // ALGO SETUP
  
  // Empty Channel Algo
  larlitecv::EmptyChannelAlgo emptyalgo;

  // Boundary Taggers
  larcv::PSet bmt_pset = pset.get<larcv::PSet>("BMTSideTagger");
  larlitecv::BoundaryMuonTaggerAlgoConfig bmt_cfg = larlitecv::MakeBoundaryMuonTaggerAlgoConfigFromPSet( bmt_pset );
  larlitecv::BoundaryMuonTaggerAlgo bmt( bmt_cfg );

  larcv::PSet fmt_pset = pset.get<larcv::PSet>("BMTFlashTagger");
  larlitecv::FlashMuonTaggerAlgoConfig fmt_cfg = larlitecv::MakeFlashMuonTaggerAlgoConfigFromPSet( fmt_pset );
  larlitecv::FlashMuonTaggerAlgo fmt_anode( larlitecv::FlashMuonTaggerAlgo::kAnode );
  larlitecv::FlashMuonTaggerAlgo fmt_cathode( larlitecv::FlashMuonTaggerAlgo::kCathode );
  larlitecv::FlashMuonTaggerAlgo fmt_imageend( larlitecv::FlashMuonTaggerAlgo::kOutOfImage );
  fmt_anode.configure( fmt_cfg );
  fmt_cathode.configure( fmt_cfg );
  fmt_imageend.configure( fmt_cfg );

  // Filters
  
  larlitecv::RadialEndpointFilter radialfilter;  // remove end points that cannot form a 3d segment nearby
  larlitecv::PushBoundarySpacePoint endptpusher; // remove endpt
  larlitecv::EndPointFilter endptfilter; // removes duplicates

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
    ntracks_2planeq = 0;
    ntracks_recod_2planeq = 0;
    ntracks_all = 0;
    ntracks_recod_all = 0;
    

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
    std::vector< larlite::event_opflash* > event_opflash_v;
    for ( auto const& prodname : flashprod ) {
      larlite::event_opflash* ev_flash = (larlite::event_opflash*)dataco[kSource].get_larlite_data(larlite::data::kOpFlash, prodname);
      std::cout << "number of flashes in " << prodname << ": " << ev_flash->size() << std::endl;
      event_opflash_v.push_back( ev_flash );
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
    // RUN RECO ALGO

    std::cout << "== Run Crossing point algorithms ===============================" << std::endl;

    // run side tagger
    //timer = std::clock();
    std::vector< larlitecv::BoundarySpacePoint > side_spacepoint_v;
    std::vector< larcv::Image2D> boundarypixel_image_v;
    std::vector< larcv::Image2D> realspacehit_image_v;
    bmt.searchforboundarypixels3D( imgs_v, badch_v, side_spacepoint_v, boundarypixel_image_v, realspacehit_image_v );
    int nsides[4] = {0};
    for ( auto const& sp : side_spacepoint_v ) {
      nsides[ sp.at(0).type ]++;
    }
    std::cout << "Reconstructed Side Tagger End Points: " << side_spacepoint_v.size() << std::endl;
    std::cout << "   Top: "        << nsides[0] << std::endl;
    std::cout << "   Bottom: "     << nsides[1] << std::endl;
    std::cout << "   Upstream: "   << nsides[2] << std::endl;
    std::cout << "   Downstream: " << nsides[3] << std::endl;

    // Declare a dummy variable for the indices of the flashes on the boundary of the image.

    // run flash tagger
    std::vector< larlitecv::BoundarySpacePoint > anode_spacepoint_v;
    std::vector< int > anode_flash_idx_v;
    std::vector< larlitecv::BoundarySpacePoint > cathode_spacepoint_v;
    std::vector< int > cathode_flash_idx_v;
    std::vector< larlitecv::BoundarySpacePoint > imgends_spacepoint_v;
    std::vector< int > imgends_flash_idx_v;
    fmt_anode.flashMatchTrackEnds(   event_opflash_v, imgs_v, badch_v, anode_spacepoint_v, anode_flash_idx_v );
    fmt_cathode.flashMatchTrackEnds( event_opflash_v, imgs_v, badch_v, cathode_spacepoint_v, cathode_flash_idx_v );
    fmt_imageend.findImageTrackEnds( imgs_v, badch_v, imgends_spacepoint_v, imgends_flash_idx_v );
    int totalflashes = (int)anode_spacepoint_v.size() + (int)cathode_spacepoint_v.size() + (int)imgends_spacepoint_v.size();
    std::cout << "Reco Flash Tagger End Points: " << totalflashes << std::endl;
    std::cout << "  Anode: "      << anode_spacepoint_v.size() << std::endl;
    std::cout << "  Cathode: "    << cathode_spacepoint_v.size() << std::endl;
    std::cout << "  Image Ends: " << imgends_spacepoint_v.size() << std::endl;

    // run end point filters
    // ---------------------------
    
    // we collect pointers to all the end points
    std::vector< const larlitecv::BoundarySpacePoint* > all_endpoints;

    // gather endpoints from space points
    for (int isp=0; isp<(int)side_spacepoint_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(side_spacepoint_v.at( isp ));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)anode_spacepoint_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(anode_spacepoint_v.at(isp));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)cathode_spacepoint_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(cathode_spacepoint_v.at(isp));
      all_endpoints.push_back( pts );
    }
    for (int isp=0; isp<(int)imgends_spacepoint_v.size(); isp++) {
      const larlitecv::BoundarySpacePoint* pts = &(imgends_spacepoint_v.at(isp));
      all_endpoints.push_back( pts );
    }

    std::cout << "number of endpoints pre-filters: " << all_endpoints.size() << std::endl;

    // first filter: radialfilter
    std::vector< const larlitecv::BoundarySpacePoint* > pass_radial_filter;
    for ( auto const& psp : all_endpoints ) {
      bool forms_segment = radialfilter.canFormSegment( (*psp).pos(), imgs_v, badch_v, 5.0, thresholds_v, 1, 2, 0.5, false ); // need parameters here
      if ( forms_segment )
        pass_radial_filter.push_back( psp );
    }
    std::cout << "number of endpoints post-radial filter: " << pass_radial_filter.size() << std::endl;

    // then we push the end points out
    std::vector< larlitecv::BoundarySpacePoint > pushed_endpoints;
    for ( auto const& psp : pass_radial_filter ) {
      const larlitecv::BoundarySpacePoint& sp = (*psp);
      if ( psp==NULL || sp.size()!=3 ) {
        std::stringstream ss;
        ss << "pass_radial_filter end point not well-formed. size=" << sp.size() << std::endl;
        throw std::runtime_error(ss.str());
      }
      larlitecv::BoundarySpacePoint pushed = endptpusher.pushPoint( sp, imgs_v, badch_v );
      if ( pushed.size()!=3 || pushed.pos().size()!=3 ) {
        std::stringstream ss;
        ss << "output of pushPoint not well-formed. size=" << pushed.size() << std::endl;
        throw std::runtime_error(ss.str());
      }
      // we keep if the pushed end point is closer to the boundary in question
      bool closer2wall = false;
      if (sp.type()==0 && sp.pos()[1]<pushed.pos()[1])
        closer2wall=true;
      else if (sp.type()==1 && sp.pos()[1]>pushed.pos()[1] )
        closer2wall=true;
      else if (sp.type()==2 && sp.pos()[2]>pushed.pos()[2] )
        closer2wall=true;
      else if (sp.type()==3 && sp.pos()[2]<pushed.pos()[2])
        closer2wall=true;
      else if (sp.type()==4 && sp.pos()[0]>pushed.pos()[0])
        closer2wall=true;
      else if (sp.type()==5 && sp.pos()[0]<pushed.pos()[0])
        closer2wall=true;
      else if (sp.type()==6 && sp.pos()[0]<=0 && sp.pos()[0]>pushed.pos()[0])
        closer2wall=true;
      else if (sp.type()==6 && sp.pos()[0]>0 && sp.pos()[0]<pushed.pos()[0])
        closer2wall=true;

      std::cout << "Pushed type=" << sp.type()
		<< " (" << sp.pos()[0] << "," << sp.pos()[1] << "," << sp.pos()[2] << ") -> "
		<< " (" << pushed.pos()[0] << "," << pushed.pos()[1] << "," << pushed.pos()[2] << ")"
		<< " closer2wall=" << closer2wall << std::endl;
      
      if ( closer2wall  )
        pushed_endpoints.push_back( pushed );
      else
        pushed_endpoints.push_back( sp );
    }
    
    std::vector< const larlitecv::BoundarySpacePoint* > pushed_ptr_endpoints;
    for ( auto const& sp : pushed_endpoints ) {
      pushed_ptr_endpoints.push_back( &sp );
    }
    
    // debug
    int itest=0;
    for ( auto const& psp : pushed_ptr_endpoints ) {
      if ( pushed_endpoints.at(itest).size()!=3 ) {
        std::stringstream ss;
        ss << __FILE__ << ":" << __LINE__
            << " stored push-point not well-formed. size=" << pushed_endpoints.at(itest).size() << " idx=" << itest << std::endl;
        throw std::runtime_error(ss.str());
      }
      if ( psp==NULL || psp->size()!=3 ) {
        std::stringstream ss;
        ss << __FILE__ << ":" << __LINE__
            << " stored push-point ptr not well-formed. psp=" << psp << " size=" << psp->size() << " idx=" << itest << std::endl;
        throw std::runtime_error(ss.str());
      }
      itest++;
    }

    // remove
    std::vector<int> endpoint_passes( pushed_ptr_endpoints.size(), 1 );
    endptfilter.removeBoundaryAndFlashDuplicates( pushed_ptr_endpoints, imgs_v, badch_v, endpoint_passes );
    endptfilter.removeSameBoundaryDuplicates( pushed_ptr_endpoints, imgs_v, badch_v, endpoint_passes );

    // remove the filtered end points
    std::vector< larlitecv::BoundarySpacePoint > side_filtered_v;
    std::vector< larlitecv::BoundarySpacePoint > anode_filtered_v;
    std::vector< larlitecv::BoundarySpacePoint > cathode_filtered_v;
    std::vector< larlitecv::BoundarySpacePoint > imgends_filtered_v;
    
    for ( size_t idx=0; idx<endpoint_passes.size(); idx++ ) {
      larlitecv::BoundarySpacePoint& sp = pushed_endpoints[idx];

      if ( endpoint_passes.at(idx)==1 ) {
        if (sp.type()<=larlitecv::kDownstream ) {
          side_filtered_v.emplace_back( std::move(sp) );
        }
        else if (sp.type()==larlitecv::kAnode) {
          anode_filtered_v.emplace_back( std::move(sp) );
        }
        else if (sp.type()==larlitecv::kCathode) {
          cathode_filtered_v.emplace_back( std::move(sp) );
        }
        else if (sp.type()==larlitecv::kImageEnd) {
          imgends_filtered_v.emplace_back( std::move(sp) );
        }
        else {
          std::stringstream ss;
          ss << __FILE__ << ":" << __LINE__ << " unrecognized boundary type" << std::endl;
          throw std::runtime_error(ss.str());
        }
      }
    }

    // combined vector
    std::vector< const std::vector<larlitecv::BoundarySpacePoint>* > filtered_spacepoints_v;
    filtered_spacepoints_v.push_back( &side_filtered_v );
    filtered_spacepoints_v.push_back( &anode_filtered_v );
    filtered_spacepoints_v.push_back( &cathode_filtered_v );
    filtered_spacepoints_v.push_back( &imgends_filtered_v );

    std::cout << "== End of Boundary Tagger ===============================" << std::endl;
   
    
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
      larlitecv::analyzeCrossingMCTracks( xingptdata, imgs_v.front().meta(),  imgs_v, ev_trigger, ev_mctrack, event_opflash_v, printFlashEnds );
      int intime_cosmics = xingptdata.true_intime_thrumu + xingptdata.true_intime_stopmu;
      std::cout << "End points from MC Truth" << std::endl;
      for (int i=0; i<larlitecv::kNumEndTypes; i++) {
	std::cout << "  " << spacepoint_producers[i] << ": " << xingptdata.true_crossingpoints[i] << std::endl;
      }
      
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
	  int radius = -1;
	  if ( xingptdata.start_crossing_nplanes_w_charge[ipt]<=1 )
	    radius = 1;
	  cv::circle( cvimgs_v[p], cv::Point(col,row), 6, cv::Scalar( 0, 0, 255, 255 ), radius );
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
	  int radius = -1;
	  if ( xingptdata.end_crossing_nplanes_w_charge[ipt]<=1 )
	    radius = 1;
	  cv::circle( cvimgs_v[p], cv::Point(col,row), 6, cv::Scalar( 255, 0, 0, 255 ), radius );
	}
      }
      // --------------------------------------------------------------------------------      
#endif

      // Compare reco and MC crossing point info
      larlitecv::analyzeCrossingMatches( xingptdata, filtered_spacepoints_v, imgs_v.front().meta(), fMatchRadius );

      // store the data into the tree
      for (int istartpt=0; istartpt<(int)xingptdata.start_type.size(); istartpt++) {
	mcxingpt_type           = xingptdata.start_type[istartpt];
	mcxingpt_matched        = xingptdata.matched_startpoint[istartpt];
	if ( xingptdata.matched_startpoint[istartpt] )
	  mcxingpt_matched      = 1;
	else
	  mcxingpt_matched      = 0;	
	mcxingpt_matched_type   = xingptdata.matched_startpoint_type[istartpt];
	mcxingpt_nplaneswcharge = xingptdata.start_crossing_nplanes_w_charge[istartpt];
	for (int p=0; p<3; p++) {
	  mcxingpt_wire[p]      = xingptdata.start_pixels[istartpt][p];
	  mcxingpt_pos[p]       = xingptdata.start_crossingpts[istartpt][p];
	}
	mcxingpt_dist           = xingptdata.start_closest_match_dist[istartpt];
	mcxingpt_tree->Fill();
      }

      for (int iendpt=0; iendpt<(int)xingptdata.end_type.size(); iendpt++) {
	mcxingpt_type           = xingptdata.end_type[iendpt];
	if ( xingptdata.matched_endpoint[iendpt] )
	  mcxingpt_matched      = 1;
	else
	  mcxingpt_matched      = 0;
	mcxingpt_matched_type   = xingptdata.matched_endpoint_type[iendpt];
	mcxingpt_nplaneswcharge = xingptdata.end_crossing_nplanes_w_charge[iendpt];
	for (int p=0; p<3; p++) {
	  mcxingpt_wire[p]      = xingptdata.end_pixels[iendpt][p];
	  mcxingpt_pos[p]       = xingptdata.end_crossingpts[iendpt][p];
	}
	mcxingpt_dist           = xingptdata.end_closest_match_dist[iendpt];
	mcxingpt_tree->Fill();
      }
      
      
    } // end of MC functions   
    
    

    tree->Fill();

#ifdef USE_OPENCV
    if ( printImages ) {
      for (int p=0; p<3; p++) {
	cv::Mat& cvimg = cvimgs_v[p];
	std::stringstream path;
	path << "boudaryptimgs/tracks_r" << run << "_s" << subrun << "_e" << event << "_p" << p << ".png";
	cv::imwrite( path.str(), cvimg );
      }
    }
#endif

  }//end of entry loop

  rfile->Write();
  //
  return 0;
}
