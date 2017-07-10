#include <iostream>
#include <vector>
#include <string>

#include "DataFormat/ImageMeta.h"
#include "DataFormat/Image2D.h"
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

#include "ThruMu/BoundaryMuonTaggerAlgoConfig.h"
#include "ThruMu/BoundaryMuonTaggerAlgo.h"


// This tests the ability of the side muon tagger algorithm
// to find deposited charge at certain locations

// We create and image where charge is deposited into a specific (y,z) position corresponding to the boundary locations
// We see if we can return the spacepoint back using the BoundaryMuonTaggerAlgo
// We're testing for bugs in the algorithm and/or missing boundary combinations

int main( int nargs, char** argv ) {

  // We go by boundary type, (top,botton,upstream,downstream), making a 5x5 pixel blob around the test location
  std::cout << "BoundaryMuonTaggerAlgo UnitTest" << std::endl;

  // =====================================================================
  // Get Config file
  if ( nargs!=2 ) {
    std::cout << "usage: ./compare_trackq [config file]" << std::endl;
    return 0;
  }
  std::string cfg_file = argv[1];

  larcv::PSet cfg_root = larcv::CreatePSetFromFile( cfg_file );
  larcv::PSet pset = cfg_root.get<larcv::PSet>("CompareTrackQ");

  // =====================================================================
  // Setup input files
  
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
  // ALGO SETUP
  
  // Boundary Taggers
  larcv::PSet bmt_pset = pset.get<larcv::PSet>("BMTSideTagger");
  larlitecv::BoundaryMuonTaggerAlgoConfig bmt_cfg = larlitecv::MakeBoundaryMuonTaggerAlgoConfigFromPSet( bmt_pset );
  larlitecv::BoundaryMuonTaggerAlgo bmt( bmt_cfg );

  // =====================================================================

  // make images
  std::vector<larcv::Image2D> img_v;
  std::vector<larcv::Image2D> badch_v;  
  for (int p=0; p<3; p++) {
    larcv::ImageMeta meta( 3456.0, 6048.0, 1008, 3456, 0, 8448.0, (larcv::PlaneID_t)p );
    larcv::Image2D img( meta );
    img.paint(0);
    img_v.emplace_back( std::move(img) );
    larcv::Image2D badch(meta);
    badch.paint(0);
    badch_v.emplace_back( std::move(badch) );
  }
  
  // TOP
  std::cout << "Testing top region" << std::endl;
  int radius = 1;
  float ypos = 104.0;
  std::vector<float> pos(3,0);
  pos[0] = 128.0;
  pos[1] = ypos;
  pos[2] = 0.5;
  while ( pos[2]<10.0 ) {

    // mark pixel
    std::vector<int> imgcoords = larcv::UBWireTool::getProjectedImagePixel( pos, img_v.front().meta(), 3 );
    for (int p=0; p<3; p++) {
      auto& img = img_v[p];
      img.paint(0);
      for (int dr=-radius; dr<=radius; dr++) {
	int row = imgcoords[0]+dr;
	if ( row<0 || row>=(int)img_v.front().meta().rows() )
	  continue;
	for (int dc=-radius; dc<=radius; dc++) {
	  int col = imgcoords[p+1]+dc;
	  if ( col<0 || col>=(int)img_v.front().meta().cols() )
	    continue;
	  img.set_pixel(row,col,20.0 );
	}
      }
    }//end of plane loop

    // run BMT
    std::vector< larlitecv::BoundarySpacePoint > side_spacepoint_v;
    std::vector< larcv::Image2D> boundarypixel_image_v;
    std::vector< larcv::Image2D> realspacehit_image_v;
    bmt.searchforboundarypixels3D( img_v, badch_v, side_spacepoint_v, boundarypixel_image_v, realspacehit_image_v );
    int nsides[4] = {0};
    for ( auto const& sp : side_spacepoint_v ) {
      nsides[ sp.at(0).type ]++;
    }
    std::cout << "----------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Test point (y,z)=(" << pos[1] << "," << pos[2] << ") imgcoords=(" << imgcoords[0] << "," << imgcoords[1] << "," << imgcoords[2] << "," << imgcoords[3] << ")" << std::endl;
    std::cout << "Reconstructed Side Tagger End Points: " << side_spacepoint_v.size() << std::endl;
    std::cout << "   Top: "        << nsides[0] << std::endl;
    std::cout << "   Bottom: "     << nsides[1] << std::endl;
    std::cout << "   Upstream: "   << nsides[2] << std::endl;
    std::cout << "   Downstream: " << nsides[3] << std::endl;

    pos[2] += 0.3;
  }
  
  
  return 0;
}
