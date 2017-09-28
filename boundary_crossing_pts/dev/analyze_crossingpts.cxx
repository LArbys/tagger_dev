
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
#include "SCE/SpaceChargeMicroBooNE.h"

// larlite
// #include "DataFormat/mctruth.h"
// #include "DataFormat/mcpart.h"
// #include "DataFormat/mctrajectory.h"
// #include "DataFormat/mctrack.h"
// #include "DataFormat/mcshower.h"
// #include "DataFormat/simch.h"
#include "DataFormat/trigger.h"
// #include "DataFormat/track.h"
#include "LArUtil/LArProperties.h"
// #include "LArUtil/Geometry.h"

// larlitecv
#include "Base/DataCoordinator.h"
#include "GapChs/EmptyChannelAlgo.h"

#include "TaggerTypes/BoundaryMuonTaggerTypes.h"
#include "TaggerTypes/BoundaryEndPt.h"
#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerTypes/BMTrackCluster3D.h"
#include "TaggerTypes/Path2Pixels.h"
#include "TaggerTypes/dwall.h"

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
#include "TaggerContourTools/BMTCV.h"
#include "TaggerContourTools/ContourAStarClusterAlgo.h"
#include "TaggerContourTools/CACAEndPtFilter.h"

// dev
#include "BMTContourFilterAlgo.h"
#include "ContourClusterAlgo.h"


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
  float max_dist2contour = pset.get<float>("MaxDist2Contour");
  bool makeCACADebugImage = pset.get<bool>("MakeCACADebugImage");
  std::vector<float> thresholds_v( 3, fthreshold );
  std::vector<float> label_thresholds_v( 3, -10 );  
  std::vector<int> label_neighborhood(3,0);
  const float cm_per_tick = ::larutil::LArProperties::GetME()->DriftVelocity()*0.5;  
  std::string spacepoint_producers[7] = { "topspacepts", "botspacepts", "upspacepts", "downspacepts", "anodepts", "cathodepts", "imgendpts" };
  
 // =====================================================================

  // // setup input
  // int kThruMu, kStopMu, kUntagged, kCROI, kNumStages;

  // setup output
  TFile* rfile = new TFile(outfname.c_str(), "recreate");
  TTree* tree = new TTree("xingpteventana", "Compare Track Charge");
  TTree* mcxingpt_tree = new TTree("mcxingptana", "Info on MC Crossing Point");
  TTree* mcxingpt_prefilter_tree  = new TTree("mcxingptana_prefilter", "Info on MC Crossing Point");      
  TTree* mcxingpt_postfilter_tree = new TTree("mcxingptana_postfilter", "Info on MC Crossing Point");    

  // Event Indexf
  int run, subrun, event;

  // Truth Quantities about interaction and lepton
  larlitecv::TruthData_t truthdata;
  
  // Crossing Point data
  larlitecv::CrossingPointAnaData_t xingptdata_prefilter;
  larlitecv::CrossingPointAnaData_t xingptdata_postfilter;
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

  // This data stores results from point of truth crossings
  int mcxingpt_type;
  int mcxingpt_matched;
  int mcxingpt_flashmatched;
  int mcxingpt_matched_type;
  int mcxingpt_nplaneswcharge;
  int mcxingpt_wire[3];
  float mcxingpt_dist;
  float mcxingpt_dwall;
  float mcxingpt_pos[3];
  TTree* xingpt_trees[3] = { mcxingpt_tree, mcxingpt_prefilter_tree, mcxingpt_postfilter_tree };
  for ( int i=0; i<3; i++) {
    xingpt_trees[i]->Branch( "truth_type", &mcxingpt_type, "truth_type/I" );
    xingpt_trees[i]->Branch( "matched", &mcxingpt_matched, "matched/I" );
    xingpt_trees[i]->Branch( "flashmatched", &mcxingpt_flashmatched, "flashmatched/I" );    
    xingpt_trees[i]->Branch( "matched_type", &mcxingpt_matched_type, "matched_type/I" );  
    xingpt_trees[i]->Branch( "nplaneswcharge", &mcxingpt_nplaneswcharge, "nplaneswcharge/I" );
    xingpt_trees[i]->Branch( "wire", mcxingpt_wire, "wire[3]/I" );
    xingpt_trees[i]->Branch( "dist", &mcxingpt_dist, "dist/F" );
    xingpt_trees[i]->Branch( "dwall", &mcxingpt_dwall, "dwall/F" );    
    xingpt_trees[i]->Branch( "pos", mcxingpt_pos, "pos[3]/F" );
  }

  // results from point of view of reco. true versus false positives.
  float reco_mindist2true = 0;
  int   reco_nplaneswcharge = 0;
  int   reco_type  = 0;
  int   reco_matched = 0;
  int   reco_matchedtype = 0;
  int   reco_passes = 0;
  TTree* recoxingpt_tree = new TTree("recoxingptana", "Reco crossing points");
  recoxingpt_tree->Branch("type",&reco_type,"type/I");
  recoxingpt_tree->Branch("matched",&reco_matched,"matched/I");  
  recoxingpt_tree->Branch("matchedtype",&reco_matchedtype,"matchedtype/I");
  recoxingpt_tree->Branch("passes",&reco_passes,"passes/I");
  recoxingpt_tree->Branch("nplaneswcharge",&reco_nplaneswcharge,"nplaneswcharge/I");
  recoxingpt_tree->Branch("mindist2true",&reco_mindist2true,"mindist2true/F");
  
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

  // BMT CV
  larlitecv::BMTCV bmtcv_algo;

  // Contour Algo
  larlitecv::ContourClusterAlgo cc_algo;
  
  // Filters  
  larlitecv::RadialEndpointFilter radialfilter;  // remove end points that cannot form a 3d segment nearby
  larlitecv::PushBoundarySpacePoint endptpusher; // remove endpt
  larlitecv::EndPointFilter endptfilter; // removes duplicates

  // space charge correctinos
  larlitecv::SpaceChargeMicroBooNE sce;
  
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
    xingptdata_prefilter.clear();
    xingptdata_postfilter.clear();
    ntracks_2planeq = 0;
    ntracks_recod_2planeq = 0;
    ntracks_all = 0;
    ntracks_recod_all = 0;
    

    // ok now to do damage

    // --------------------------------------------------------------------------------------------------------------------
    // LOAD IMAGES 
    // get the original, segmentation, and tagged images
    larcv::EventImage2D* ev_imgs   = (larcv::EventImage2D*)dataco[kSource].get_larcv_data(larcv::kProductImage2D,inputimgs);
    larcv::EventImage2D* ev_badch  = NULL;
    larcv::EventImage2D* ev_gapch  = NULL;    
    larcv::EventImage2D* ev_segs   = NULL;
    if ( ismc ) {
      ev_segs = (larcv::EventImage2D*)dataco[kSource].get_larcv_data(larcv::kProductImage2D,"segment");
    }
    if ( badch_in_file ) {
      // we either get the badch from the file
      ev_gapch = (larcv::EventImage2D*)dataco[kCROIfile].get_larcv_data(larcv::kProductImage2D,"gapchs");
    }
    else {
      // or we have to make the badch image from a ChStatus object
      larcv::EventChStatus* ev_status = (larcv::EventChStatus*)dataco[kSource].get_larcv_data( larcv::kProductChStatus, "wire" );
      std::vector<larcv::Image2D> chstatus_img_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
      std::vector<larcv::Image2D> gapch_v = emptyalgo.findMissingBadChs( ev_imgs->Image2DArray(), chstatus_img_v, 10.0, 5 );
      ev_badch = new larcv::EventImage2D;
      ev_badch->Emplace( std::move(chstatus_img_v) );
      ev_gapch = new larcv::EventImage2D;
      ev_gapch->Emplace( std::move(gapch_v) );
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
      int p=0;
      for ( auto const& img : imgs_v ) {
	cv::Mat cvimg = larcv::as_mat_greyscale2bgr( img, fthreshold, 100.0 );
	const larcv::Image2D& badch = badch_v[p];
	for (int i=0; i<(int)badch.meta().cols(); i++) {
	  if ( badch.pixel(10,i)>0 ) {
	    for (int r=0; r<badch.meta().rows(); r++) {
	      cvimg.at<cv::Vec3b>(cv::Point(i,r))[0] = 10;
	      cvimg.at<cv::Vec3b>(cv::Point(i,r))[1] = 10;
	      cvimg.at<cv::Vec3b>(cv::Point(i,r))[2] = 10;	      
	    }
	  }
	}
	cvimgs_v.emplace_back(std::move(cvimg));
	p++;
      }
    }
#endif

    // --------------------------------------------------------------------------------------------------------------------
    // OPFLASHES
    
    // get the opflashes
    std::vector< larlite::event_opflash* > event_opflash_v;
    for ( auto const& prodname : flashprod ) {
      larlite::event_opflash* ev_flash = (larlite::event_opflash*)dataco[kSource].get_larlite_data(larlite::data::kOpFlash, prodname);
      std::cout << "number of flashes in " << prodname << ": " << ev_flash->size() << std::endl;
      event_opflash_v.push_back( ev_flash );
    }

    // --------------------------------------------------------------------------------------------------------------------
    // MC TRUTH/TRIGGER

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

    if ( ismc ) {
      
      // loop over MC tracks, get end points of muons
      larlitecv::analyzeCrossingMCTracks( xingptdata, imgs_v.front().meta(),  imgs_v, ev_trigger, ev_mctrack, event_opflash_v, printFlashEnds );
      // copy to other structs
      xingptdata_prefilter = xingptdata;
      xingptdata_postfilter = xingptdata;
      
      int intime_cosmics = xingptdata.true_intime_thrumu + xingptdata.true_intime_stopmu;
      std::cout << "End points from MC Truth" << std::endl;
      for (int i=0; i<larlitecv::kNumEndTypes; i++) {
	std::cout << "  " << spacepoint_producers[i] << ": " << xingptdata.true_crossingpoints[i] << std::endl;
      }

      // --------------------------------------------------------------------------------
      // OPEN CV VISUALIZATION
#ifdef USE_OPENCV
      // draw truth points
      for ( int ipt=0; ipt<(int)xingptdata.truthcrossingptinfo_v.size(); ipt++) {
	larlitecv::TruthCrossingPointAna_t& info = xingptdata.truthcrossingptinfo_v[ipt];
	int row=info.imgcoord[0];
	if (row<0) row = 0;
	if (row>=(int)imgs_v.front().meta().rows() ) row = (int)imgs_v.front().meta().rows()-1;
	for (int p=0; p<3; p++) {
	  int col = info.imgcoord[p+1];
	  if ( col<0 ) col = 0;
	  if ( col>=(int)imgs_v.front().meta().cols() ) col = imgs_v.front().meta().cols()-1;

	  // label
	  std::stringstream ptname;	  
	  if ( info.start_or_end==0 )
	    ptname << "#S" << ipt;
	  else
	    ptname << "#E" << ipt;

	  // color
	  cv::Scalar color( 0, 0, 0, 255 );
	  if ( info.start_or_end==0 ) {
	    if ( info.type==4 || info.type==5 )
	      color = cv::Scalar( 0, 255, 255, 255 );
	    else
	      color = cv::Scalar( 0,   0, 255, 255 );
	  }
	  else {
	    if ( info.type==4 || info.type==5 )
	      color = cv::Scalar( 255, 255, 0, 255 );
	    else
	      color = cv::Scalar( 255, 0, 0, 255 );
	  }	      

	  float radius;
	  
	  if ( info.nplanes_w_charge>=2 ) {
	    // findable truth point
	    cv::circle( cvimgs_v[p], cv::Point(col,row), 6, color, radius );
	  }
	  else {
	    cv::rectangle( cvimgs_v[p], cv::Point(col-3,row-3), cv::Point(col+3,row+3), color, 2 );
	  }
	  cv::putText( cvimgs_v[p], cv::String(ptname.str()), cv::Point( col+2, row+2 ), cv::FONT_HERSHEY_SIMPLEX, 0.5, color );
	}
      }
      // --------------------------------------------------------------------------------      
#endif
      
    } // end of MC functions   
    
    // =======================================================================
    // RUN RECO ALGO

    std::cout << "== Run Crossing point algorithms ===============================" << std::endl;

    // preprocess test
    // binarize and blur image
    std::vector<larcv::Image2D> binblur_v;
    for (auto const& img : imgs_v ) {
      larcv::Image2D bbimg( img );
      bbimg.paint(0);
      for (int row=0; row<(int)img.meta().rows(); row++) {
	for (int col=0; col<(int)img.meta().cols(); col++) {
	  if ( img.pixel(row,col)>fthreshold ) {
	    for (int dr=-2; dr<=2; dr++) {
	      int r = row+dr;
	      if ( r<0 || r>=(int)img.meta().rows() )
		continue;
	      for (int dc=-2; dc<=2; dc++) {
		int c = col+dc;
		if ( c<0 || c>=(int)img.meta().cols() )
		  continue;
		bbimg.set_pixel( r, c, 255 );
	      }
	    }
	  }
	}
      }//end of row loop
      binblur_v.emplace_back( std::move(bbimg) );
    }//end of img loop
    
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
    std::cout << "  Anode end points" << std::endl;    
    for (int i=0; i<(int)anode_spacepoint_v.size(); i++) {
      std::cout << "    [" << i << "] flashidx=" << anode_spacepoint_v.at(i).getFlashIndex() << std::endl;
    }
    
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
    std::vector< const std::vector<larlitecv::BoundarySpacePoint>*  > prefilter_spacepoints_v;
    prefilter_spacepoints_v.push_back( &side_spacepoint_v );
    prefilter_spacepoints_v.push_back( &anode_spacepoint_v );
    prefilter_spacepoints_v.push_back( &cathode_spacepoint_v );
    prefilter_spacepoints_v.push_back( &imgends_spacepoint_v );


    // ---------------------------------------------------------------------------------
    // Compare reco and MC crossing point info
    
    if ( ismc ) {
      larlitecv::analyzeCrossingMatches( xingptdata_prefilter, prefilter_spacepoints_v, imgs_v.front().meta(), fMatchRadius );

      // store the pre-filter data into the tree
      int numpts = xingptdata_prefilter.truthcrossingptinfo_v.size();
      for (int ipt=0; ipt<numpts; ipt++) {
      	larlitecv::TruthCrossingPointAna_t& info = xingptdata_prefilter.truthcrossingptinfo_v[ipt];
	
      	mcxingpt_type           = info.type;
      	mcxingpt_matched        = info.matched;
      	if ( info.flashindex>=0 )
      	  mcxingpt_flashmatched = 1;
      	else
      	  mcxingpt_flashmatched = 0;
      	mcxingpt_matched_type   = info.matched_type;
      	mcxingpt_nplaneswcharge = info.nplanes_w_charge;
      	for (int p=0; p<3; p++) {
      	  mcxingpt_wire[p]      = info.imgcoord[p+1];
      	  mcxingpt_pos[p]       = info.crossingpt_det[p];
      	}
      	mcxingpt_dist           = info.matched_dist;
      	xingpt_trees[1]->Fill();
      }
    }
      
    // --------------------------------------------------------------------------------
    // RECO DEV:
    // BMTCV
    bmtcv_algo.analyzeImages( imgs_v, badch_v, 10.0, 2 );

    // RECO DEV: Contour-based Filters
    
    // Pass Truth End points through the filter    
    larlitecv::BMTContourFilterAlgo contour_truthfilter_algo;
    std::vector<larcv::Image2D> truthfilter_clusterpix_v;

    int numtruth_in_contour = 0;
    int numtruth_valid_contour = 0;
    int numtruth = 0;
    for ( int testindex=0; testindex<(int)xingptdata_prefilter.truthcrossingptinfo_v.size(); testindex++) {
      larlitecv::TruthCrossingPointAna_t& info = xingptdata_prefilter.truthcrossingptinfo_v[testindex];
      std::vector<float> testpt(3,0);
      for (int i=1; i<3; i++)
      	testpt[i] =  info.crossingpt_detsce_tyz[i];
      testpt[0] = (info.crossingpt_detsce_tyz[0]-3200.0)*cm_per_tick;
      bool incontour = contour_truthfilter_algo.buildCluster( imgs_v, badch_v, truthfilter_clusterpix_v, testpt, bmtcv_algo.m_plane_atomicmeta_v, max_dist2contour );
      if ( incontour )
	numtruth_in_contour++;
      if ( contour_truthfilter_algo.cuts[1] )
	numtruth_valid_contour++;
      numtruth++;
    }
    
    std::cout << "BMT Contour Filter on truth points (in contour): " << float(numtruth_in_contour)/float(numtruth) << std::endl;
    std::cout << "BMT Contour Filter on truth points (consitent contours): " << float(numtruth_valid_contour)/float(numtruth) << std::endl;    
    //cc_algo.buildCluster( imgs_v, badch_v, clusterpix_v, testpt,  bmtcv_algo.m_plane_atomicmeta_v );

    // PASS RECO POINTS THROUGH CONTOUR FILTER

    larlitecv::BMTContourFilterAlgo contour_recofilter_algo;
    std::vector<larcv::Image2D> recofilter_clusterpix_v;
    std::vector< int > recofilter_passes( xingptdata_prefilter.recocrossingptinfo_v.size(), 0 );
    std::vector< int > recofilter_good( xingptdata_prefilter.recocrossingptinfo_v.size(), 0 );    
    std::vector< larlitecv::BoundarySpacePoint > contourpassing_spacepoints_v;

    int numreco_good_in_contour = 0;
    int numreco_bad_in_contour = 0;    
    int numreco_good_tot = 0;
    int numreco_bad_tot  = 0;    
    int ireco = -1;
    int numreco_good[7][2] = {0};
    int numreco_bad[7][2]  = {0};    
    for ( auto const& p_sp_v : prefilter_spacepoints_v ) {
      for (auto const& sp : *p_sp_v ) {
	ireco++;
	larlitecv::RecoCrossingPointAna_t& recoinfo = xingptdata_prefilter.recocrossingptinfo_v[ireco];

	std::vector<float> testpt(3,0);
	for (int i=0; i<3; i++)
	  testpt[i] =  sp.pos()[i];
	std::cout << "============ [ RECO #" << ireco << " ] =============" << std::endl;
	bool incontour = contour_recofilter_algo.buildCluster( imgs_v, badch_v, recofilter_clusterpix_v, testpt, bmtcv_algo.m_plane_atomicmeta_v, max_dist2contour );
	bool seedclusterok = contour_recofilter_algo.cuts[1];

	// override pass all
	incontour = true;
	seedclusterok = true;
	
	if ( recoinfo.truthmatch==1 ) {
	  numreco_good_tot++;
	  numreco_good[sp.type()][1]++;
	  recofilter_good[ireco] = 1;
	  if ( incontour && seedclusterok) {
	    numreco_good_in_contour++;
	    numreco_good[sp.type()][0]++;
	  }
	}
	else {
	  numreco_bad_tot++;
	  numreco_bad[sp.type()][1]++;
	  if ( incontour && seedclusterok) {
	    numreco_bad_in_contour++;
	    numreco_bad[sp.type()][0]++;	    
	  }
	}
	
	bool passes = incontour & seedclusterok;
	if ( passes ) {
	  recofilter_passes[ireco] = 1;
	  contourpassing_spacepoints_v.push_back( sp );
	}
	
#ifdef USE_OPENCV
	// Visualize
	std::vector<int> testptpix = larcv::UBWireTool::getProjectedImagePixel( testpt, imgs_v.front().meta(), 3 );
	//std::cout << "testptpix: (" << testptpix[0] << "," << testptpix[1] << "," << testptpix[2] << "," << testptpix[3] << ")" << std::endl;

	std::stringstream ptname;
	if ( sp.type()==0 )
	  ptname << "#T" << ireco;
	else if ( sp.type()==1 )
	  ptname << "#B" << ireco;
	else if ( sp.type()==2 )
	  ptname << "#U" << ireco;
	else if ( sp.type()==3 )
	  ptname << "#D" << ireco;
	else if ( sp.type()==4 )
	  ptname << "#A" << ireco;
	else if ( sp.type()==5 )
	  ptname << "#C" << ireco;
	else
	  ptname << "#I" << ireco;

	for (int p=0; p<3; p++) {
	  if ( passes ) {
	    cv::circle( cvimgs_v[p], cv::Point(testptpix[p+1],testptpix[0]), 3, cv::Scalar( 255, 0, 255, 255 ), 1 );
	    cv::putText( cvimgs_v[p], cv::String(ptname.str()), cv::Point(testptpix[p+1]+2,testptpix[0]), cv::FONT_HERSHEY_SIMPLEX, 0.3, cv::Scalar(255,0,255,255) );
	  }
	  else  {
	    cv::circle( cvimgs_v[p], cv::Point(testptpix[p+1],testptpix[0]), 3, cv::Scalar( 51, 151, 255, 255 ), 1 );
	    cv::putText( cvimgs_v[p], cv::String(ptname.str()), cv::Point(testptpix[p+1]+2,testptpix[0]), cv::FONT_HERSHEY_SIMPLEX, 0.3, cv::Scalar(51,151,255,255) );	    	    
	  }
	}
#endif
	
      }//end of sp loop
    }//end of spacepoint v
    std::cout << "BMT Contour Filter on reco points ----------------- " << std::endl;
    std::cout << "  Good Points in contour: " << numreco_good_in_contour << " of " << numreco_good_tot << std::endl;
    std::cout << "    top: " << numreco_good[0][0] << "/" << numreco_good[0][1] << std::endl;
    std::cout << "    bot: " << numreco_good[1][0] << "/" << numreco_good[1][1] << std::endl;
    std::cout << "    ups: " << numreco_good[2][0] << "/" << numreco_good[2][1] << std::endl;
    std::cout << "    dwn: " << numreco_good[3][0] << "/" << numreco_good[3][1] << std::endl;
    std::cout << "    anode: " << numreco_good[4][0] << "/" << numreco_good[4][1] << std::endl;
    std::cout << "    cathode: " << numreco_good[5][0] << "/" << numreco_good[5][1] << std::endl;
    std::cout << "    imgends: " << numreco_good[6][0] << "/" << numreco_good[6][1] << std::endl;            
    std::cout << "  Bad Points in contour:  " << numreco_bad_in_contour << " of " << numreco_bad_tot << std::endl;
    std::cout << "    top: " << numreco_bad[0][0] << "/" << numreco_bad[0][1] << std::endl;
    std::cout << "    bot: " << numreco_bad[1][0] << "/" << numreco_bad[1][1] << std::endl;
    std::cout << "    ups: " << numreco_bad[2][0] << "/" << numreco_bad[2][1] << std::endl;
    std::cout << "    dwn: " << numreco_bad[3][0] << "/" << numreco_bad[3][1] << std::endl;
    std::cout << "    anode: " << numreco_bad[4][0] << "/" << numreco_bad[4][1] << std::endl;
    std::cout << "    cathode: " << numreco_bad[5][0] << "/" << numreco_bad[5][1] << std::endl;
    std::cout << "    imgends: " << numreco_bad[6][0] << "/" << numreco_bad[6][1] << std::endl;                
    std::cout << "--------------------------------------------------- " << std::endl;
    
    // ===============================================================================================================
    //  DEV: ContourAStarCluster

    std::cout << "================================================================" << std::endl;
    std::cout << "================================================================" << std::endl;
    std::cout << " DEV: ContourAStarCluster" << std::endl;

    std::vector< std::vector<int> > caca_results;
    larlitecv::CACAEndPtFilter cacaalgo;
    cacaalgo.setVerbosity(2);
    cacaalgo.setTruthInformation( xingptdata_prefilter.truthcrossingptinfo_v, xingptdata_prefilter.recocrossingptinfo_v );
    if ( makeCACADebugImage )
      cacaalgo.makeDebugImage();
    cacaalgo.evaluateEndPoints( prefilter_spacepoints_v, event_opflash_v, imgs_v, badch_v, bmtcv_algo.m_plane_atomicmeta_v, 150.0, caca_results );
    std::vector<larlitecv::BoundarySpacePoint> cacapassing_spacepoints_v;
    std::vector<larlitecv::BoundarySpacePoint> cacapassing_moved_v = cacaalgo.regenerateFitleredBoundaryPoints( imgs_v );
    
    std::vector<int> cacapasses_unrolled;
    ireco = -1;
    for (size_t iv=0; iv<prefilter_spacepoints_v.size(); iv++) {
      const std::vector<larlitecv::BoundarySpacePoint>* psp_v = prefilter_spacepoints_v[iv];
      for (size_t ipt=0; ipt<psp_v->size(); ipt++) {
	ireco++;
	if ( caca_results[iv][ipt]==1) {
	  cacapassing_spacepoints_v.push_back( psp_v->at(ipt) );
	  cacapasses_unrolled.push_back(1);
	}
	else
	  cacapasses_unrolled.push_back(0);
	
      }
    }

    if ( makeCACADebugImage ) {

      std::stringstream path;
      path << "boundaryptimgs/astarcluster_goodrecopts_r" << run << "_s" << subrun << "_e" << event << "_p" << 4 << ".png";
      cv::imwrite( path.str(), cacaalgo.getDebugImage(0) );
      path.str("");
      path << "boundaryptimgs/astarcluster_badrecopts_r" << run << "_s" << subrun << "_e" << event << "_p" << 4 << ".png";
      cv::imwrite( path.str(), cacaalgo.getDebugImage(1) );
    }

    
    // // visualize cluster images
    // cv::Mat cvimg_clustimg(imgs_v.front().meta().rows(),imgs_v.front().meta().cols(),CV_8UC3);
    // cv::Mat cvimg_pathimg( imgs_v.front().meta().rows(),imgs_v.front().meta().cols(),CV_8UC3);
    // for (int r=0; r<(int)imgs_v.front().meta().rows(); r++) {
    //   for (int c=0; c<(int)imgs_v.front().meta().cols(); c++) {
    // 	for (int p=0; p<3; p++) {
    // 	  cvimg_clustimg.at<cv::Vec3b>( cv::Point(c,r) )[p] = astar_cluster.m_cvimg_v[p].at<uchar>( cv::Point(c,r) );
    // 	  cvimg_pathimg.at<cv::Vec3b>(  cv::Point(c,r) )[p] = astar_cluster.m_cvpath_v[p].at<uchar>( cv::Point(c,r) );
    // 	}
    //   }
    // }
    // path.str("");
    // path << "boundaryptimgs/astarcluster_clust_r" << run << "_s" << subrun << "_e" << event << "_p" << 4 << ".png";
    // cv::imwrite(path.str(),cvimg_clustimg );
    // path.str("");
    // path << "boundaryptimgs/astarcluster_path_r" << run << "_s" << subrun << "_e" << event << "_p" << 4 << ".png";
    // cv::imwrite(path.str(),cvimg_pathimg );

    std::cout << "===================================================" << std::endl;
    std::cout << "===================================================" << std::endl;
    cacaalgo.printStageTimes();
    
    // =================================================================================================================

    
    // dump the images
    // print
    for (int p=0; p<3; p++) {
      std::stringstream path;
      path << "boundaryptimgs/cvbmt_r" << run << "_s" << subrun << "_e" << event << "_p" << p << ".png";
      cv::imwrite( path.str(), bmtcv_algo.cvimg_stage0_v[p] );
    }
    

    // =================================================================================================================
    // POST-FILTER ANALYSIS
    if ( ismc ) {
      // Evaluate post-filter points from the truth perspective
      
      std::vector< const std::vector<larlitecv::BoundarySpacePoint>*  > postfilter_spacepoints_v;
      //postfilter_spacepoints_v.push_back( &contourpassing_spacepoints_v );
      postfilter_spacepoints_v.push_back( &cacapassing_spacepoints_v );
      larlitecv::analyzeCrossingMatches( xingptdata_postfilter, postfilter_spacepoints_v, imgs_v.front().meta(), fMatchRadius );

      // store the post-filter data into the tree
      int numpts = xingptdata_postfilter.truthcrossingptinfo_v.size();
      for (int ipt=0; ipt<numpts; ipt++) {
	larlitecv::TruthCrossingPointAna_t& info = xingptdata_postfilter.truthcrossingptinfo_v[ipt];
	
	mcxingpt_type           = info.type;
	mcxingpt_matched        = info.matched;
	if ( info.flashindex>=0 )
	  mcxingpt_flashmatched = 1;
	else
	  mcxingpt_flashmatched = 0;
	mcxingpt_matched_type   = info.matched_type;
	mcxingpt_nplaneswcharge = info.nplanes_w_charge;
	for (int p=0; p<3; p++) {
	  mcxingpt_wire[p]      = info.imgcoord[p+1];
	  mcxingpt_pos[p]       = info.crossingpt_det[p];
	}
	mcxingpt_dist           = info.matched_dist;
	xingpt_trees[2]->Fill();
      }

      // store the post-filter reco data into the tree
      numpts = (int)xingptdata_prefilter.recocrossingptinfo_v.size();
      for (int ireco=0; ireco<numpts; ireco++) {
	larlitecv::RecoCrossingPointAna_t& info = xingptdata_prefilter.recocrossingptinfo_v[ireco];
	reco_type = info.type;
	reco_matched = info.truthmatch;
	reco_matchedtype = info.truthmatch_type;
	reco_mindist2true = info.truthmatch_dist;
	reco_nplaneswcharge = 3;
	//reco_passes = recofilter_passes[ireco];
	reco_passes = cacapasses_unrolled[ireco];
	recoxingpt_tree->Fill();
      }
    }
        
    // DRAW THE BOUNDARY POINTS IMAGES
    std::cout << "BoundaryPixel Images: " << boundarypixel_image_v.size() << std::endl;
    std::cout << "RealSpaceHit Images: " << realspacehit_image_v.size()  << std::endl;    
    std::vector<cv::Mat> cvboundarypix_v;
    std::vector<cv::Mat> cvrealspace_v;
 
#ifdef USE_OPENCV         
    if ( printImages ) {
      for ( auto const& img : imgs_v ) {
	cv::Mat cvimg = larcv::as_mat_greyscale2bgr( img, fthreshold, 100.0 );
	cvboundarypix_v.emplace_back(std::move(cvimg));
      }

      int type_colors[4][3] = { {204,204,255}, // top
				{255,229,204}, // bot
				{153,204,255}, // upstream
				{205,255,255} }; // downstream
      for (int itype=0; itype<4; itype++) {
	for (int p=0; p<3; p++) {
	  auto const& bp_img = boundarypixel_image_v.at(itype*3+p);
	  auto & cv_img = cvboundarypix_v.at(p);
	  for (size_t row=0; row<bp_img.meta().rows(); row++) {
	    for (size_t col=0; col<bp_img.meta().cols(); col++) {
	      if ( bp_img.pixel(row,col)>0 ) {
		for (int i=0; i<3; i++)
		  cv_img.at<cv::Vec3b>( cv::Point(col,row) )[i] = type_colors[itype][i];
	      }
	    }//end of col loop
	  }//end of row loop
	}//end of plane loop
      }//end of type loop

      // print
      for (int p=0; p<3; p++) {
	std::stringstream path;
	path << "boundaryptimgs/boundarypix_r" << run << "_s" << subrun << "_e" << event << "_p" << p << ".png";
	cv::imwrite( path.str(), cvboundarypix_v[p] );
      }

      for ( auto const& img : realspacehit_image_v ) {
	cv::Mat cvimg = larcv::as_mat_greyscale2bgr( img, 0, 3.0 );
	cvrealspace_v.emplace_back(std::move(cvimg));
      }

      if ( ismc ) {
	// draw truth points
	for (int ipt=0; ipt<(int)xingptdata.truthcrossingptinfo_v.size(); ipt++) {
	  larlitecv::TruthCrossingPointAna_t& info = xingptdata_prefilter.truthcrossingptinfo_v[ipt];
	  int itype = info.type;
	  if ( itype>=4 )
	    continue;
	  std::vector<float>& crossingpt_sce = info.crossingpt_detsce;
	  std::vector<int>&   crossingpx     = info.imgcoord;

	  float x = 0;
	  if ( itype==0 || itype==1 ) {
	    // top and bottom
	    x = crossingpt_sce[2];
	  }
	  else {
	    // upstream downstream
	    x = crossingpt_sce[1]+117.0;
	  }
	  int row=crossingpx[0];

	  std::stringstream ptname;	  
	  if ( info.start_or_end==0 ) {
	    ptname << "#S" << ipt;	    
	    if ( info.matched==1 )
	      cv::circle( cvrealspace_v[itype], cv::Point(x,row), 6, cv::Scalar( 0, 0, 255, 255 ), 2.0 );
	    else
	      cv::circle( cvrealspace_v[itype], cv::Point(x,row), 6, cv::Scalar( 255, 0, 0, 255 ), 2.0 );
	  }
	  else {
	    ptname << "#E" << ipt;	    
	    if ( info.matched==1 )
	      cv::circle( cvrealspace_v[itype], cv::Point(x,row), 6, cv::Scalar( 0, 0, 255, 255 ), 2.0 );
	    else
	      cv::circle( cvrealspace_v[itype], cv::Point(x,row), 6, cv::Scalar( 255, 0, 0, 255 ), 2.0 );
	  }
	  cv::putText( cvrealspace_v[itype], cv::String(ptname.str()), cv::Point( x+2, row+2 ), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255,255,255,255) );	  	  	  
	}
      }

      // draw reco points
      ireco = -1;
      for ( auto const& p_sp_v : prefilter_spacepoints_v ) {
	for (auto const& sp : *p_sp_v ) {
	  ireco++;
	  if ( sp.type()>=4 )
	    continue;
	  std::stringstream ptname;
	  if ( sp.type()==0 )
	    ptname << "#T" << ireco;
	  else if ( sp.type()==1 )
	    ptname << "#B" << ireco;
	  else if ( sp.type()==2 )
	    ptname << "#U" << ireco;
	  else if ( sp.type()==3 )
	    ptname << "#D" << ireco;
	  else if ( sp.type()==4 )
	    ptname << "#A" << ireco;
	  else if ( sp.type()==5 )
	    ptname << "#C" << ireco;
	  else
	    continue; // not relavent here

	  float x = 0;
	  if ( sp.type()==0 || sp.type()==1 ) {
	    // top and bottom
	    x = sp.pos()[2];
	  }
	  else {
	    // upstream downstream
	    x = sp.pos()[1]+117.0;
	  }
	  int row=sp[0].row;
	  cv::Scalar colour(125,125,255,255);
	  if ( recofilter_passes[ireco]==0 )
	    colour = cv::Scalar(255,125,125,255);
	  cv::circle(  cvrealspace_v[sp.type()], cv::Point(x,row), 2, colour, -1.0 );
	  cv::putText( cvrealspace_v[sp.type()], cv::String(ptname.str()), cv::Point( x+2, row+2 ), cv::FONT_HERSHEY_SIMPLEX, 0.3, colour );
	}
      }
      
      // print
      for (int p=0; p<4; p++) {
	std::stringstream path;
	path << "boundaryptimgs/realspace_r" << run << "_s" << subrun << "_e" << event << "_side" << p << ".png";
	cv::imwrite( path.str(), cvrealspace_v[p] );
      }
      
    }//end of if print images
#endif
    // ----------------------------------------------------------------------------------
    // End of Reco Algorithms
    // =======================================================================    

    
    tree->Fill();

#ifdef USE_OPENCV
    if ( printImages ) {
      for (int p=0; p<3; p++) {
	cv::Mat& cvimg = cvimgs_v[p];
	std::stringstream path;
	path << "boundaryptimgs/tracks_r" << run << "_s" << subrun << "_e" << event << "_p" << p << ".png";
	cv::imwrite( path.str(), cvimg );
      }
    }
#endif

  }//end of entry loop

  rfile->Write();
  //
  return 0;
}
