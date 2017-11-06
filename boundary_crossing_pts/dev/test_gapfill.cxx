#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "DataFormat/ImageMeta.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventChStatus.h"
#include "Base/PSet.h"
#include "Base/LArCVBaseUtilFunc.h"
#include "UBWireTool/UBWireTool.h"

// larlitecv
#include "Base/DataCoordinator.h"
#include "TaggerTypes/BoundaryMuonTaggerTypes.h"
#include "TaggerTypes/BoundaryEndPt.h"
#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerTypes/BMTrackCluster3D.h"
#include "TaggerTypes/Path2Pixels.h"
#include "TaggerTypes/dwall.h"
#include "TaggerContourTools/BMTCV.h"
#include "GapChs/EmptyChannelAlgo.h"

#include "ContourGapFill.h"

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>


int main( int nargs, char** argv ) {

  // We go by boundary type, (top,botton,upstream,downstream), making a 5x5 pixel blob around the test location
  std::cout << "Test Gapfill" << std::endl;

  // =====================================================================
  // Get Config file
  if ( nargs!=2 ) {
    std::cout << "usage: ./test_gapfill [config file]" << std::endl;
    return 0;
  }
  std::string cfg_file = argv[1];

  larcv::PSet cfg_root = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg_root.get<larcv::PSet>("TestGapFill");

  // =====================================================================
  // Setup input files
  
  enum SourceTypes_t { kSource=0, kNumSourceTypes };
  std::string source_param[2] = { "InputSourceFilelist" };
  larlitecv::DataCoordinator dataco[kNumSourceTypes];
  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    std::cout << "LOADING " << source_param[isrc] << " FILES" << std::endl;
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArCV"),   "larcv" );
    dataco[isrc].set_filelist( pset.get<std::string>(source_param[isrc]+"LArLite"), "larlite" );
    dataco[isrc].configure( cfg_file, "StorageManager", "IOManager", "TestGapFill" );
    dataco[isrc].initialize();
  }

  for (int isrc=0; isrc<kNumSourceTypes; isrc++ ) {
    std::cout << "data[" << source_param[isrc] << "] entries=" << dataco[isrc].get_nentries("larcv") << std::endl;
  }
  

  // =====================================================================
  // configuration parameters

  std::string outfname                = pset.get<std::string>("OutputAnaFile");
  int         verbosity               = pset.get<int>("Verbosity",0);
  std::string inputimgs               = pset.get<std::string>("Image2DProducer");
  std::string chstatus                = pset.get<std::string>("ChStatusProducer");
  std::string trigname                = pset.get<std::string>("TriggerProducer");  
  std::vector<std::string> flashprod  = pset.get<std::vector<std::string> >("OpFlashProducerList");
  bool        ismc                    = pset.get<bool>("IsMC");
  std::string instanceproducer        = pset.get<std::string>("InstanceProducer");
  std::string larcv_inputlist         = pset.get<std::string>("InputSourceFilelistLArCV");
  std::string larlite_inputlist       = pset.get<std::string>("InputSourceFilelistLArLite");
  
  
  // =====================================================================
  // ALGO SETUP

  // Contour maker
  larlitecv::BMTCV bmtcv_algo;
  
  // Gap Filling Algo
  larlitecv::ContourGapFill gapfillalgo;
  gapfillalgo.makeDebugImage( true );

  // Empty Channel Algo
  larlitecv::EmptyChannelAlgo emptyalgo;
  

  // =====================================================================

  int nentries = dataco[kSource].get_nentries("larcv");

  for (int ientry=0; ientry<nentries; ientry++) {

    // load the entry
    dataco[kSource].goto_entry(ientry,"larcv");
    int run, subrun, event;
    dataco[kSource].get_id(run, subrun, event);

    // --------------------------------------------------------------------------------------------------------------------
    // LOAD IMAGES 
    // get the original, segmentation, and tagged images
    larcv::EventImage2D* ev_imgs   = (larcv::EventImage2D*)dataco[kSource].get_larcv_data(larcv::kProductImage2D,inputimgs);
    larcv::EventImage2D* ev_badch  = NULL;
    larcv::EventImage2D* ev_gapch  = NULL;    
    larcv::EventImage2D* ev_segs   = NULL;
    if ( ismc ) {
      ev_segs  = (larcv::EventImage2D*)dataco[kSource].get_larcv_data(larcv::kProductImage2D,"segment");
    }

    // or we have to make the badch image from a ChStatus object
    larcv::EventChStatus* ev_status = (larcv::EventChStatus*)dataco[kSource].get_larcv_data( larcv::kProductChStatus, "wire" );
    std::vector<larcv::Image2D> chstatus_img_v = emptyalgo.makeBadChImage( 4, 3, 2400, 6048, 3456, 6, 1, *ev_status );
    std::vector<larcv::Image2D> gapch_v = emptyalgo.findMissingBadChs( ev_imgs->Image2DArray(), chstatus_img_v, 10.0, 5 );
    ev_badch = new larcv::EventImage2D;
    ev_badch->Emplace( std::move(chstatus_img_v) );
    ev_gapch = new larcv::EventImage2D;
    ev_gapch->Emplace( std::move(gapch_v) );
    
    const std::vector<larcv::Image2D>& img_v   = ev_imgs->Image2DArray();
    const std::vector<larcv::Image2D>& badch_v = ev_badch->Image2DArray();
        
    bmtcv_algo.analyzeImages( img_v, badch_v, 10.0, 2 );

    std::vector<larcv::Image2D> gapfilled_v;
    gapfillalgo.makeGapFilledImages( img_v, badch_v, bmtcv_algo.m_plane_atomicmeta_v, gapfilled_v );

    std::vector<cv::Mat>& cvimg_v = gapfillalgo.getDebugImages();

    std::stringstream ss;
    ss << "gapfillimgs/gapalgodebug_run" << run << "_subrun" << subrun << "_ev" << event << ".png";
    cv::imwrite( ss.str(), cvimg_v.front() );
    std::cout << "wrote: " << ss.str() << std::endl;
  }
  
    
  return 0;
}
