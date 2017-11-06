#include "ContourGapFill.h"

// larlite
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

// larcv
#include "UBWireTool/UBWireTool.h"
#include "CVUtil/CVUtil.h"


namespace larlitecv {

  ContourGapFill::ContourGapFill() {
    fMakeDebugImage = false;
    m_meta_v.clear();
  }

  ContourGapFill::~ContourGapFill() {
  }

  void ContourGapFill::makeDebugImage( bool debug ) {
    fMakeDebugImage = true;
  }
  
  void ContourGapFill::makeGapFilledImages( const std::vector<larcv::Image2D>& img_v, const std::vector<larcv::Image2D>& badch_v,
					    const std::vector< std::vector<larlitecv::ContourShapeMeta> >& plane_contours_v,
					    std::vector<larcv::Image2D>& gapfilled_v ) {

    // Steps
    // ------
    // -- define badch spans
    // -- find contours on the sides of the spans
    // -- match contours across a given span which point to one another
    // -- AStar across the span?
    // -- fill the gap

    // store meta information on the image
    if ( m_meta_v.size()!=img_v.size() ) {
      for ( auto const& img : img_v ) {
	m_meta_v.push_back( img.meta() );
      }
    }

    // allocate a debug image if the flag is on
    if ( fMakeDebugImage ) {
      m_cvimg_debug_v.clear();
      m_plane_contour_v.clear();
      
      createDebugImage( img_v );
      // we need to make a list of contours of type contour, not contourshapemeta

      for ( size_t p=0; p<plane_contours_v.size(); p++) {
	std::vector< std::vector<cv::Point> > contour_v;
	for ( auto const& ctr : plane_contours_v[p] )
	  contour_v.push_back( ctr );
	m_plane_contour_v.emplace_back( std::move(contour_v) );
      }
    }

    // define the bad channel spans on each plan
    std::vector< std::vector<BadChSpan> > plane_span_v;
    for ( auto const& badch : badch_v ) {
      std::vector<BadChSpan> span_v = findBadChSpans( badch, 2 );
      //std::cout << "spans on plane=" << badch.meta().plane() << ": " << span_v.size() << std::endl;
      plane_span_v.emplace_back( std::move(span_v) );      
    }

    // find the contours on each plane that touch the spans
    for (size_t p=0; p<m_meta_v.size(); p++) {
      associateContoursToSpans( plane_contours_v[p], m_meta_v[p], plane_span_v[p], 2 );
    }
    
  }
  
  void ContourGapFill::createDebugImage( const std::vector<larcv::Image2D>& img_v ) {
    // we create a cv::Mat image to draw on
    // we make an rgb image
    
    cv::Mat cvimg = larcv::as_mat_greyscale2bgr( img_v.front(), 0, 255.0 );
    for (size_t p=1; p<img_v.size(); p++) {
      for (size_t r=0; r<m_meta_v[p].rows(); r++) {
	for (size_t c=0; c<m_meta_v[p].cols(); c++) {
	  if ( img_v[p].pixel(r,c)>5.0 ) {
	    for (int i=0; i<3; i++)
	      cvimg.at<cv::Vec3b>(cv::Point(c,r))[i] = img_v[p].pixel(r,c);
	  }
	}
      }
    }
    m_cvimg_debug_v.emplace_back(std::move(cvimg));
  }
  
  std::vector< ContourGapFill::BadChSpan > ContourGapFill::findBadChSpans( const larcv::Image2D& badch, int goodchwidth ) {
    // we get a list of bad channel spans for the image.
    const larcv::ImageMeta& meta = badch.meta();
    size_t row = meta.rows()/2;

    // simple on-off region finder, where end of region requires goodchwidth of good channels to be observed
    std::vector< BadChSpan > span_v;
    bool inregion = false;
    for ( size_t c=0; c<meta.cols(); c++) {

      if ( !inregion ) {
	// not yet in a badch region
	if ( badch.pixel(row,c)>0 ) {
	  // now we're on. create region

	  BadChSpan newspan;
	  newspan.start = c;
	  newspan.end   = c;
	  newspan.width = 0;
	  newspan.ngood = 0;
	  newspan.planeid = (int)meta.plane();
	  span_v.emplace_back( std::move(newspan) );
	  inregion = true;
	}
      }
      else {
	// currently in a badch region
	BadChSpan& last = span_v.back();
	if ( badch.pixel(row,c)>0 ) {
	  // still in bad region
	  last.end = c;
	  last.ngood = 0; //< reset good ch counter
	}
	else {
	  // now in good region
	  last.ngood++;
	  if ( last.ngood>=goodchwidth ) {
	    // we end the region
	    last.end = (int)c-goodchwidth;
	    last.width = last.end-last.start+1;
	    inregion = false; // reset the inregion marker
	  }
	  else {
	    // do nothing, wait for another good col
	  }
	}//end of if in good
      }//end of if inregion
    }//end of col loop

    if ( inregion ) {
      // if still in a region
      BadChSpan& last = span_v.back();
      last.end = (int)meta.cols()-1;
      last.width = last.end-last.start+1;
    }
    
    return span_v;
  }

  void ContourGapFill::associateContoursToSpans( const std::vector<larlitecv::ContourShapeMeta>& contour_v,
						 const larcv::ImageMeta& meta,
						 std::vector<ContourGapFill::BadChSpan>& span_v,
						 const int colwidth ) {
    
    // inputs
    // ------
    // contour_v: contour list for a plane
    // span_v: bad channel span list for the same plane
    // colwidth: number of columns away from span edges to look for contours
    //
    // outputs
    // -------
    // span_v: spans are updated with the left/right-ctrindx vectors filled

    int nleftctrs  = 0;
    int nrightctrs = 0;
    int ispan = -1;
    for ( auto& span : span_v ) {
      ispan++;
      // we scan down the start and end col and look to see if it touches a contour

      //std::cout << "span #" << ispan << ": [" << span.start << "," << span.end << "]" << std::endl;

      // start
      for ( int c=span.start-colwidth+1; c<=span.start; c++) {
	if ( c<0 || c>=meta.cols() ) continue;
	for (size_t r=0; r<meta.rows(); r++) {

	  for (int idx=0; idx<(int)contour_v.size(); idx++) {
	    auto const& ctr = contour_v[idx];

	    cv::Point testpt( c, (int)r );
	    double dist = cv::pointPolygonTest( ctr, testpt, false );
	    if ( dist>0 ) {
	      // inside contour
	      span.leftctridx.insert(idx);
	      nleftctrs++;
	    }	    
	  }
	}
      }//end of start loop

      // end
      for ( int c=span.end; c<span.end+colwidth; c++) {
	if ( c<0 || c>=meta.cols() ) continue;
	for (size_t r=0; r<meta.rows(); r++) {
	  
	  for (int idx=0; idx<(int)contour_v.size(); idx++) {
	    auto const& ctr = contour_v[idx];
	    
	    cv::Point testpt( c, (int)r );
	    double dist = cv::pointPolygonTest( ctr, testpt, false );
	    if ( dist>0 ) {
	      // inside contour
	      span.rightctridx.insert(idx);
	      nrightctrs++;	      
	    }	    
	  }
	}
      }//end of end loop
      //std::cout << "ncontours close to this span: left=" << span.leftctridx.size() << " right=" << span.rightctridx.size() << std::endl;
      
      if ( fMakeDebugImage ) {
	//std::cout << "Draw contours on " << m_cvimg_debug_v.size() << " debug images" << std::endl;
	cv::Mat& cvimg = m_cvimg_debug_v.front();
	
      	cv::Scalar contourcolor;
	if ( meta.plane()==0 )
	  contourcolor = cv::Scalar(255,0,0,255);
	else if ( meta.plane()==1 )
	  contourcolor = cv::Scalar(0,255,0,255);
	else if ( meta.plane()==2 )
	  contourcolor = cv::Scalar(0,0,255,255);
	
	for ( auto const& idx : span.leftctridx ) {
	  cv::drawContours( cvimg, m_plane_contour_v[meta.plane()], idx, contourcolor, -1 );
	}
	for ( auto const& idx : span.rightctridx ) {
	  //std::cout << "  draw ctr idx=" << idx << "(of " << contour_v.size() << ")" << std::endl;
	  cv::drawContours( cvimg, m_plane_contour_v[meta.plane()], idx, contourcolor, -1 );
	}
      }
    }

    return;
  }
						 

}
