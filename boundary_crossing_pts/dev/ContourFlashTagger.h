#ifndef __CONTOUR_FLASH_TAGGER_H__
#define __CONTOUR_FLASH_TAGGER_H__

#include <vector>

// larlite data
#include "Base/DataFormatConstants.h"
#include "DataFormat/opflash.h"

// larcv data
#include "DataFormat/EventImage2D.h"

// larlitecv header
#include "TaggerTypes/BoundaryMuonTaggerTypes.h"
#include "TaggerTypes/BoundarySpacePoint.h"
#include "TaggerContourTools/BMTCV.h"


namespace larlitecv {

  class ContourFlashTagger {
  public:
    ContourFlashTagger();
    virtual ~ContourFlashTagger();
    
    bool flashMatchTrackEnds( const std::vector< larlite::event_opflash* >& opflashsets, const std::vector<larcv::Image2D>& tpc_imgs,
			      const std::vector<larcv::Image2D>& badch_imgs, const larlitecv::BMTCV& img_contours,
			      std::vector< BoundarySpacePoint >& trackendpts, std::vector< int > endpoint_flash_idx );

  };
  
}


#endif
