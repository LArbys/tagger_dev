#include <iostream>

#include <string>

// ROOT
#include "TString.h"
#include "TImage.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"

// larcv
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventPixel2D.h"
#include "DataFormat/EventChStatus.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"

// larlite
#include "DataFormat/mctrack.h"
#include "DataFormat/trigger.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"
#include "OpT0Finder/App/MCQCluster.h"
#include "FhiclLite/FhiclLiteUtilFunc.h"
#include "OpT0Finder/Algorithms/QLLMatch.h"

// larlitecv
#include "Base/DataCoordinator.h"
#include "GapChs/EmptyChannelAlgo.h"
#include "TaggerTypes/BoundaryMuonTaggerTypes.h"
#include "TaggerTypes/BoundaryEndPt.h"
#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerTypes/BMTrackCluster3D.h"
#include "TaggerTypes/Path2Pixels.h"
#include "ThruMu/ThruMuTracker.h"
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
#include "GeneralFlashMatchAlgo/GeneralFlashMatchAlgoConfig.h"
#include "GeneralFlashMatchAlgo/GeneralFlashMatchAlgo.h"
#include "ThruMu/FlashMuonTaggerAlgoConfig.h"
#include "ThruMu/FlashMuonTaggerAlgo.h"


int main( int nargs, char** argv ){
  
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

  const larutil::Geometry* geo = ::larutil::Geometry::GetME();
    
  // setup output
  TFile* rfile = new TFile(outfname.c_str(), "recreate");
  TTree* tree = new TTree("compareq", "Compare Track Charge");

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

  // truth end point track reco metrics
  int ntracks_2planeq = 0;
  int ntracks_recod_2planeq = 0;
  int ntracks_all = 0;
  int ntracks_recod_all = 0;
  tree->Branch( "ntracks_2planeq", &ntracks_2planeq, "ntracks_2planeq/I" );
  tree->Branch( "ntracks_recod_2planeq", &ntracks_recod_2planeq, "ntracks_recod_2planeq/I" );
  tree->Branch( "ntracks_all", &ntracks_all, "ntracks_all/I" );
  tree->Branch( "ntracks_recod_all", &ntracks_recod_all, "ntracks_recod_all/I" );  
  
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

  // Declare variables that will be needed in the future.
  float chi2;
  float min_dchi2;

  // A tree that will contain the chi2 value for all of the tracks.
  // The total number of entries in this tree better equal the sum of the entries of the other two trees.
  TTree* totalchi2tree  = new TTree("TotalChi2Tree", "Chi2 Value for All Matched Tracks");
  totalchi2tree->Branch("chi2", &chi2, "chi2/F");

  TTree* flashmatchtree = new TTree("FlashMatchTree", "Flash Match Info w/ MC Tracks");
  flashmatchtree->Branch("chi2", &chi2, "chi2/F");
  flashmatchtree->Branch("min_dchi2", &min_dchi2, "min_dchi2/F");

  float track_num;
  float num_flashes_with_better_chi2;

  // Declare a tree for the tracks that have a lower chi2 when compared to other tracks from the event
  TTree* worsechi2tree  = new TTree("WorseChi2Tree", "Tracks With a Better Flash w/ other Flashes in the Event");
  worsechi2tree->Branch("run",&run,"run/I");
  worsechi2tree->Branch("subrun",&subrun,"subrun/I");
  worsechi2tree->Branch("event",&event,"event/I");
  worsechi2tree->Branch("track_num", &track_num, "track_num/I");
  worsechi2tree->Branch("num_flashes_with_better_chi2", &num_flashes_with_better_chi2, "num_flashes_with_better_chi2/I");

  int event_idx;
  int track_idx;
  int num_of_tracks_with_better_chi2;
  int total_num_tracks;

  // Declare the canvas up here.
  TCanvas *c = new TCanvas("c","c",600.600);
  
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

  // Flash Matching Algo                                                        
  larcv::PSet flash_match_pset = pset.get<larcv::PSet>("GeneralFlashMatchAlgo");
  // Flash Match Algo Config object.                                            
  larlitecv::GeneralFlashMatchAlgoConfig flash_match_algo_config_object;
  // Make the PSet from this object.                                            
  larlitecv::GeneralFlashMatchAlgoConfig flash_match_algo_config_pset = flash_match_algo_config_object.MakeConfigFromPSet(flash_match_pset);
  // GeneralFlashMatchAlgo object.  This will be used below to find a flash hypothesis for tracks.       
  larlitecv::GeneralFlashMatchAlgo flash_match_obj( flash_match_algo_config_pset );

  // Set up a PSet for the MCQCluster functionality as well.
  // Declare the name for the object of type 'MCQCluster'.                                                                                                                                              
  flashana::MCQCluster MCQCluster_object("MCQCluster");
  
  flashana::QLLMatch* QLLMatch_obj = ::flashana::QLLMatch::GetME();

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
  
  // Set 'event_idx'.
  event_idx = 0;

  // Declare a vector for the number of times that you have looked at a good chi2 track.                                                                                                                
  int num_of_good_chi2_tracks = 0;
  
  for (int ientry=startentry; ientry<endentry; ientry++) {

    // Break after the tenth event.
    //if (ientry > 10) break;

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
    // MC INFO ANALYSIS

    // quantities filled if MC present
    std::vector<larcv::Image2D> nupix_imgs_v;
    std::vector<int> vertex_col(3,-1);
    std::vector<double> vtx_sce(3,0);
    int vertex_row = -1;

    // Declare a vector of all the reconstructed endpoints for the tracks.
    std::vector < larlitecv::BoundarySpacePoint > mctrack_endpts_all_tracks;
    mctrack_endpts_all_tracks.clear();

    // Vector for 'mctrack' indices.
    std::vector < int > mctrack_idx_v;
    mctrack_idx_v.clear();

    if ( ismc ) {

      // loop over MC tracks, get end points of muons
      larlitecv::analyzeCrossingMCTracks( xingptdata, imgs_v.front().meta(),  imgs_v, ev_trigger, ev_mctrack, event_opflash_v, printFlashEnds );
      int intime_cosmics = xingptdata.true_intime_thrumu + xingptdata.true_intime_stopmu;
      std::cout << "End points from MC Truth" << std::endl;
      for (int i=0; i<larlitecv::kNumEndTypes; i++) {
	std::cout << "  " << spacepoint_producers[i] << ": " << xingptdata.true_crossingpoints[i] << std::endl;
      }

      // loop over MC tracks. Get the truth end points, run thrumu tagger!
      std::cout << ev_mctrack->size() << " = " << xingptdata.mctrack_imgendpoint_indices.size() << std::endl;
      std::vector< larlitecv::BMTrackCluster3D > mc_recotracks; // mc reco track
      std::vector<int> mc_recotracks_truthindex;
      std::vector< const larlite::opflash* > truth_flash_ptrs; // the flash that goes with each element in mc_recotracks
      std::cout << "MC track imgendpoint indices: " << xingptdata.mctrack_imgendpoint_indices.size() << std::endl;
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
	int start_nplanes_wcharge = 0;
	int end_nplanes_wcharge   = 0;	
	if (  start_index!=-1 ) {
	  // make a boundary space point object for the start;

	  // get truth crossingpt info
	  const larlitecv::TruthCrossingPointAna_t& truthxing = xingptdata.truthcrossingptinfo_v[start_index];
	  start_nplanes_wcharge = truthxing.nplanes_w_charge;
	  
	  std::vector< larlitecv::BoundaryEndPt > planepts;
	  int row = truthxing.imgcoord[0];
	  if (row<0 || row>=imgs_v.front().meta().rows()) {
	    std::cout << "BAD START POINT: (" << truthxing.imgcoord[0] << "," << truthxing.imgcoord[1] << "," << truthxing.imgcoord[2] << "," << truthxing.imgcoord[3] << ")" << std::endl;
	    continue;
	  }
	  bool goodpt = true;
	  for (int p=0; p<3; p++) {
	    int col = truthxing.imgcoord[p+1];
	    if ( col<0 || col>=imgs_v.front().meta().cols() ) {
	      std::cout << "BAD START POINT: (" << truthxing.imgcoord[0] << "," << truthxing.imgcoord[1] << "," << truthxing.imgcoord[2] << "," << truthxing.imgcoord[3] << ")" << std::endl;
	      goodpt = false;
	      break;
	    }
	    larlitecv::BoundaryEndPt planept( row, col, (larlitecv::BoundaryEnd_t)truthxing.type );
	    planepts.emplace_back( std::move(planept) );
	  }
	  if ( !goodpt )
	    continue;
	  // make the boundary spacepoint
	  pstartpt = new larlitecv::BoundarySpacePoint( (larlitecv::BoundaryEnd_t)truthxing.type, std::move(planepts), imgs_v.front().meta() );
	  mctrack_endpts.push_back( pstartpt );
	}
	if ( end_index!=-1 ) {
	  // make a boundary space point object for the end;

	  // get truth crossing information
	  const larlitecv::TruthCrossingPointAna_t& truthxing = xingptdata.truthcrossingptinfo_v[end_index];
	  end_nplanes_wcharge = truthxing.nplanes_w_charge;	  
	  
	  std::vector< larlitecv::BoundaryEndPt > planepts;
	  int row = truthxing.imgcoord[0];
	  if (row<0 || row>=imgs_v.front().meta().rows()) {
	    std::cout << "BAD END POINT: (" << truthxing.imgcoord[0] << "," << truthxing.imgcoord[1] << "," << truthxing.imgcoord[2] << "," << truthxing.imgcoord[3] << ")" << std::endl;
	    // need to cleanup start point
	    if ( pstartpt )
	      delete pstartpt;
	    continue;
	  }
	  bool goodpt = true;
	  for (int p=0; p<3; p++) {
	    int col = truthxing.imgcoord[p+1];
	    if ( col<0 || col>=imgs_v.front().meta().cols() ) {
	      std::cout << "BAD END POINT: (" << truthxing.imgcoord[0] << "," << truthxing.imgcoord[1] << "," << truthxing.imgcoord[2] << "," << truthxing.imgcoord[3] << ")" << std::endl;
	      goodpt = false;
	      break;
	    }
	    larlitecv::BoundaryEndPt planept( row, col, (larlitecv::BoundaryEnd_t)truthxing.type );
	    planepts.emplace_back( std::move(planept) );
	  }
	  if ( !goodpt ) {
	    delete pstartpt;
	    continue;
	  }
	  pendpt = new larlitecv::BoundarySpacePoint( (larlitecv::BoundaryEnd_t)truthxing.type, std::move(planepts), imgs_v.front().meta() );
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
	  std::cout << " nplanes_qcharge[ start=" << start_nplanes_wcharge << " end=" << end_nplanes_wcharge << "]" << std::endl;
	  ntracks_all++;
	  if ( start_nplanes_wcharge>=2 && end_nplanes_wcharge>=2)
	    ntracks_2planeq++;
	  try {
	    thrumualgo.makeTrackClusters3D( imgs_v, badch_v, mctrack_endpts, trackclusters, tagged_v, used_endpoints_indices );
	    if ( trackclusters.size()>0 ) {
	      // successful reco returned!
	      ntracks_recod_all++;
	      if ( start_nplanes_wcharge>=2 && end_nplanes_wcharge>=2)
		ntracks_recod_2planeq++;
	      
	      // Append the endpoints onto the 'mctrack_endpts_all_tracks' by using the information in 'trackclusters.at(0)'.
	      larlitecv::BoundarySpacePoint start_endpt_mctrack = trackclusters.at(0).start_endpts;
	      larlitecv::BoundarySpacePoint end_endpt_mctrack = trackclusters.at(0).end_endpts;
              mctrack_endpts_all_tracks.push_back(start_endpt_mctrack);
              mctrack_endpts_all_tracks.push_back(end_endpt_mctrack);
	      mc_recotracks.emplace_back( std::move(trackclusters.at(0)) );
	      mc_recotracks_truthindex.push_back( itrack );

	      // Append all of the information for the track endpoints at this point in the algorithm.
	      // Append the information for the two endpoints to 'mctrack_endpts_all_tracks' in a loop.
	      //for ( size_t endpt_iter = 0; endpt_iter < mctrack_endpts.size(); endpt_iter++ ) {
	      //mctrack_endpts_all_tracks.push_back( mctrack_endpts[endpt_iter] );
	      //}

	      // get its flash
	      const std::vector<int>& endpoint_indices = xingptdata.mctrack_imgendpoint_indices[itrack];	      
	      int flashindex = -1;
	      if ( endpoint_indices.size()>0 && endpoint_indices[0]!=-1 ) {
		flashindex = xingptdata.truthcrossingptinfo_v[ endpoint_indices[0] ].flashindex;
	      }

	      const larlite::opflash* popflash = NULL;
	      if ( flashindex!=-1 ) {
		int nflash = 0;
		for ( auto const& evflash_v : event_opflash_v ) {
		  if ( flashindex < nflash+(int)evflash_v->size() ) {
		    popflash = &((*evflash_v)[ flashindex-nflash ]);
		    break;
		  }
		  nflash += (int)evflash_v->size();
		}
	      }
	      
	      mctrack_idx_v.push_back( itrack );
	      truth_flash_ptrs.push_back( popflash );
	    }
	  }
	  catch ( std::exception& e ) {
	    std::cout << std::endl;
	    std::cout << "  ThruMu Failed: " << e.what() << std::endl;
	  }
	}

	// we have to delete the end points as we used new
	for (int i=0; i<(int)mctrack_endpts.size(); i++) {
	  delete mctrack_endpts[i];
	  mctrack_endpts[i] = NULL;	  
	}
	mctrack_endpts.clear();

	std::cout << "tracker return " << trackclusters.size() << " thrumu tracks" << std::endl;
	
	std::cout << std::endl;
      } //end of mctrack loop

      // --------------------------------------------------------------------------------
      // Get Chi-2 of flash hypo from reco track and truth flash

      // Make the 'ev_mctrack' pointer into an 'ev_mctrack' object (to be used in the check here and in the normal method below).

      //larlite::event_mctrack ev_mctrack_val       = *ev_mctrack; // note from tmw. this is a copy. you can do the same thing
      larlite::event_mcshower ev_mcs;
      flashana::LightPath lightpath;

      // Use the 'Configure' function with this object to convert the MCTrack objects into qclusters.
      // This function is a void, so it just sets the input data types to the function.
      // tmw note: you can derefence the pointer in the argument, which will pass by reference. no copy involved.
      MCQCluster_object.Construct( *ev_mctrack, ev_mcs, lightpath );

      // Extract the qclusters that are within 5 us of one of the flashes.
      std::vector < flashana::QCluster_t > qcluster_v;
      qcluster_v.clear();

      qcluster_v = MCQCluster_object.QClusters();

      std::cout << "The size of 'qcluster_v' = " << qcluster_v.size() << "." << std::endl;
      std::cout << "The size of 'qcluster_v' = " << MCQCluster_object.QClusters().size() << "." << std::endl;      
      std::cout << "The size of 'truth_flash_pointers' = " << truth_flash_ptrs.size() << "." << std::endl;

      // Declare a new vector of truth_flash_pointers that only includes those that are matched to a QCluster.
      std::vector< const larlite::opflash* > truth_flash_ptrs_parsed;
      truth_flash_ptrs_parsed.clear();

      // Declare a vector of the indices of the 'truth_flash_ptrs_parsed' that have already been appended.
      std::vector < int > truth_flash_ptrs_parsed_idx;
      truth_flash_ptrs_parsed_idx.clear();

      // Declare a new vector of qclusters that only includes those matched to a flash.
      std::vector < flashana::QCluster_t > qcluster_v_parsed;
      qcluster_v_parsed.clear();

      // Declare a matching parameter.
      double time_match = 5.0;

      // Declare a counter for the number of null flashes that there are.
      int num_null_flashes;
      int num_times_in_flash_loop;

      // Start a loop over the elements in 'qcluster_v' to compare them to the flashes contained in 'truth_flash_ptrs'.
      // Check to see if the flash is NULL here instead of in the loop over the qclusters later.
      for ( size_t i = 0; i < qcluster_v.size(); i++ ) {

	// Declare a variable for if the qcluster has been matched yet.
	bool already_matched = false;

	// Declare a new qcluster for the qcluster at this location in 'qcluster_v'.                                                                                                                      
	const flashana::QCluster_t qcluster             = qcluster_v.at(i);

	num_times_in_flash_loop = 0;
	num_null_flashes        = 0;

	for ( size_t j = 0; j < truth_flash_ptrs.size(); j++ ) {

	  num_times_in_flash_loop++;

	  //std::cout << "Starting a loop over the truth flash pointers." << std::endl;

	  // Generate an individual flash pointer and an individual flash out of the information in 'truth_flash_ptrs.size()'.
	  const larlite::opflash*    data_opflash_pointer = truth_flash_ptrs.at(j);
	  // Moving the line to set this opflash object equal to a pointer to the point after it is checked to be null.

	  // Continue if the flash is NULL.
	  if (truth_flash_ptrs.at( j ) == NULL) {
	    num_null_flashes++;
	    continue;                                                                                                                                                                                      
        }

	  // Set the 'data_opflash' object here.
	  const larlite::opflash     data_opflash         = *data_opflash_pointer;

	  // Compare the time of 'qcluster' to the time of 'data_opflash'.  If they are within a matching parameter (say 5 us), then you can append them (in the case of the flash, the pointer) to the respective vectors.
	  if ( fabs( qcluster.time - data_opflash.Time() ) < time_match ) {

	    std::cout << "The time of these two objects match!!" << std::endl;
	    std::cout << "qcluster.time = " << qcluster.time << "." << std::endl;
	    std::cout << "data_opflash.Time() = " << data_opflash.Time() << "." << std::endl;
	    std::cout << "\n" << std::endl;

	    // Check to see if that 'truth_flash_pointer' has already been appended to the 'truth_flash_ptrs_parsed' list.
	    for ( size_t redundancy_check = 0; redundancy_check < truth_flash_ptrs_parsed_idx.size(); redundancy_check++ ) {

	      if ( j == truth_flash_ptrs_parsed_idx.at( redundancy_check ) )
		continue;

	    }
	    
	    // Append them to their respective vectors.
	    qcluster_v_parsed.push_back( qcluster );
	    truth_flash_ptrs_parsed.push_back ( data_opflash_pointer );
	    truth_flash_ptrs_parsed_idx.push_back ( j );
	    break;

	  }

	}

      }

      std::cout << "The number of entries in the original 'qcluster_v' vector = " << qcluster_v.size() << "." << std::endl;
      std::cout << "The number of entries in the original 'truth_flash_ptrs' vector = " << truth_flash_ptrs.size() << "." << std::endl;
      std::cout << "The number of null flashes = " << num_null_flashes << "." << std::endl;
      
      // Compare the size of 'qcluster_v' and the size of the truth flash pointers vector.
      std::cout << "The size of qcluster_v_parsed = " << qcluster_v_parsed.size() << "." << std::endl;
      std::cout << "The number of truth_flash_pointers_parsed = " << truth_flash_ptrs_parsed.size() << "." << std::endl;
      
      // Declare the number of cm/tick outside the loop.
      const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;
      std::cout << "cm_per_tick:  "<< cm_per_tick << std::endl;

      // Loop through these tracks and find the flash hypothesis for each of the qcluster values. \
      // Set the value of 'track_idx'.
      track_idx = 0;

      for ( size_t qcluster_iter = 0; qcluster_iter <  qcluster_v_parsed.size(); qcluster_iter++ ) {

	// Set to 0 the variables that must be plotted to find the new performance of the flash-matcher.
	min_dchi2 = -1.;
	track_num = 0;
	num_flashes_with_better_chi2 = 0;

	// Set an object equal to 'truth_flash_ptr' at position 'qcluster_iter'.
	// I have to use 'mctrack_idx' here instead of 'qcluster_iter'.  This is predicated on the fact that the 'mctrack_idx_v' vector is filled in the same order as 'truth_flash_ptrs'.
	const larlite::opflash* data_opflash_pointer = truth_flash_ptrs_parsed.at(qcluster_iter);

	// If this is true, then find the object that the pointer points to.
	const larlite::opflash  data_opflash         = *data_opflash_pointer;

	// Convert the opflash object into a data flash.
	flashana::Flash_t truth_flash = flash_match_obj.MakeDataFlash( data_opflash );

	// Convert the qcluster into a data flash.
	flashana::Flash_t hypo_flash   = flash_match_obj.GenerateUnfittedFlashHypothesis( qcluster_v_parsed.at(qcluster_iter) );

	chi2 = QLLMatch_obj->QLL( hypo_flash, truth_flash );

	// Print out the chi2 value down here now.
	std::cout << "The chi2 value = " << chi2 << "." << std::endl;

	// Plot the hypothesis flash information and the truth flash information on the same set of axes.
	if ( chi2 < 3.5 && num_of_good_chi2_tracks < 20)  {

	  std::cout << "Making an image of the PEs inside the PMTs." << std::endl;
	  
	  TCanvas *c     = new TCanvas;
          TH1F    *data  = new TH1F("dataflash", Form("chi2 = %f", chi2), 32, 0, 32);
          TH1F    *hypo  = new TH1F("hypoflash", Form("chi2 = %f", chi2), 32, 0, 32);

          float max_data_plot = -1.0;
          float max_hypo_plot = -1.0;

	  unsigned int opdet = 0;
	  
          // Fill the histograms in a loop over the optical detectors.
	  for ( unsigned int opch = 0; opch < 32; opch++ ) {

	    // Convert the opch into an opdet.
	    opdet = geo->OpDetFromOpChannel(opch); 
	    
            data->SetBinContent( opdet, truth_flash.pe_v[opdet] );
            hypo->SetBinContent( opdet, hypo_flash.pe_v[opdet] ); // Do not multiply by the fudge factor with this hypothesis!!

            // Find the max in the two plots.                                                                                                                                                        
            if ( max_data_plot < truth_flash.pe_v[opdet] || max_data_plot < 0. )
              max_data_plot = truth_flash.pe_v[opdet];


	    if ( max_hypo_plot < hypo_flash.pe_v[opdet] || max_hypo_plot < 0. ) // Do not multiply by the fudge factor with this hypothesis!!
	      max_hypo_plot = hypo_flash.pe_v[opdet]; // Do not multiply by the fudge factor with this hypothesis!!

	  }

	  // Draw the two histograms on the same axes.
	  c->cd();
	  data->SetLineColor(kBlue);
	  hypo->SetLineColor(kRed);

	  if ( max_data_plot > max_hypo_plot ) {
	    // Draw both of the histograms on c.
	    data->Draw();
	    hypo->Draw("Sames");
	  }

	  if ( max_data_plot < max_hypo_plot ) {
	    hypo->Draw();
	    data->Draw("Sames");
	  }
	  
	  // Save the canvas as an image.
	  TImage *img = TImage::Create();
	  
	  img->FromPad(c);
	  
	  img->WriteImage(Form("output_png_images_hypo_truth_comparison/optimal_chi2_hypo_truth_comparison_%d.png",num_of_good_chi2_tracks));

	  // Increment 'good_chi2_tracks_plotted'.
	  num_of_good_chi2_tracks++;

	  std::cout << "The number of good chi2 tracks after making the last image = " << num_of_good_chi2_tracks << "." << std::endl;

          // Delete the root infrastructure that you used in this loop.
	  delete c;
          delete data;
          delete hypo;
          delete img;

        }

	// Loop through the other flashes in the event to see if another one has a better chi2 match with 
	float fake_chi2 = 0.0;

	// Set 'track_num' equal to 'track_idx'.
	track_num = track_idx;

	for ( size_t fake_flash_iter = 0; fake_flash_iter < truth_flash_ptrs_parsed.size(); fake_flash_iter++ ) {

	  // Continue if the flash is the same as the true flash - would produce the same chi2 value as before.
	  if ( fake_flash_iter == qcluster_iter ) continue;

	  const larlite::opflash* fake_opflash_pointer = truth_flash_ptrs_parsed.at(fake_flash_iter);
	  const larlite::opflash  fake_opflash         = *fake_opflash_pointer;

	  // Turn the fake opflash into a data flash.
	  flashana::Flash_t fake_flash = flash_match_obj.MakeDataFlash( fake_opflash );
	
	  fake_chi2   = QLLMatch_obj->QLL( hypo_flash, fake_flash );

	  // Print out the fake chi2 value.
	  std::cout << "The fake chi2 value = " << fake_chi2 << "." << std::endl;

	  // Compare 'fake_chi2' and 'chi2'.
	  if ( fake_chi2 < chi2 ) {

	    num_flashes_with_better_chi2++;

	    std::cout << "A better chi2 match is found with flash # " << fake_flash_iter << "." << std::endl;

	  }

	  if ( chi2 < fake_chi2 || chi2 == fake_chi2 ) {

	    if ( fabs( fake_chi2 - chi2 ) < min_dchi2 || min_dchi2 < 0. ) 
	      min_dchi2 = fabs( fake_chi2 - chi2 );

	  }

	}

	// Fill the tree based on the value of 'num_flashes_with_better_chi2'.  If it is == 0, then fill 'flashmatchtree'.  If it is greater than 1, then fill 'worsechi2tree'.
	if ( num_flashes_with_better_chi2 == 0 ) 
	  flashmatchtree->Fill();

	if ( num_flashes_with_better_chi2 > 0 )
	  worsechi2tree->Fill();

	// Fill 'totalchi2tree'.
	totalchi2tree->Fill();
       
	// Increment 'track_idx'.  This variable is set equal to 'track_num'.
	track_idx++;

      } // End of the loop over qclusters.
    

	// CHRIS' CODE GOES HERE
      // --------------------------------------------------------------------------------
      
      // --------------------------------------------------------------------------------
      // OPEN CV VISUALIZATION
#ifdef USE_OPENCV
      // draw crossing points
      for ( int ipt=0; ipt<(int)xingptdata.truthcrossingptinfo_v.size(); ipt++) {
	
	const larlitecv::TruthCrossingPointAna_t& truthxing = xingptdata.truthcrossingptinfo_v[ipt];
	
	// color of marker depends on end point
	cv::Scalar markercolor(0,0,0,0);
	if ( truthxing.start_or_end==0 )
	  markercolor = cv::Scalar(0,0,255,255);
	else
	  markercolor = cv::Scalar(255,0,0,255);
	
	int row=truthxing.imgcoord[0];
	if (row<0) row = 0;
	if (row>=(int)imgs_v.front().meta().rows() ) row = (int)imgs_v.front().meta().rows()-1;
	for (int p=0; p<3; p++) {
	  int col = truthxing.imgcoord[p+1];
	  if ( col<0 ) col = 0;
	  if ( col>=(int)imgs_v.front().meta().cols() ) col = imgs_v.front().meta().cols()-1;
	  int radius = -1;
	  if ( truthxing.nplanes_w_charge<=1 )
	    radius = 1;
	  cv::circle( cvimgs_v[p], cv::Point(col,row), 6, markercolor, radius );
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

    // Increment 'event_idx'.
    event_idx++;

    // Increment 'total_num_tracks' by 'track_idx + 1'.
    //total_num_tracks += track_idx + 1;
    
  } //end of entry loop

  //std::cout << "The total number of tracks in this sample = " << total_num_tracks << "." << std::endl;
  //std::cout << "The number of tracks with a better chi2 than the original = " << num_tracks_better_new_chi2 << "." << std::endl;
  //std::cout << "The number of tracks with a worse chi2 than the original = " << num_tracks_worse_new_chi2 << "." << std::endl;
  //std::cout << "The number of tracks with a better extended chi2 = " << num_tracks_better_extended_chi2 << std::endl;
      
  rfile->Write();
  //
  return 0;
}
