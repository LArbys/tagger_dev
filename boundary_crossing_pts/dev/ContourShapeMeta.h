#ifndef __ContourShapeMeta__
#define __ContourShapeMeta__

#include <vector>

#include "DataFormat/ImageMeta.h"

#ifdef USE_OPENCV
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#endif


namespace larlitecv {

 class ContourShapeMeta : public std::vector<cv::Point> {
   // Wrapper around OpenCV contour
   // Stores meta data for contours
   
 public:

   ContourShapeMeta( const std::vector<cv::Point>& contour, const larcv::ImageMeta& img );
   virtual ~ContourShapeMeta() {};

   const larcv::ImageMeta& meta() { return m_meta; };    
   const cv::Point getFitSegmentStart() { return m_start; };
   const cv::Point getFitSegmentEnd() { return m_end; };   

 protected:
   
   // ImageMeta
   const larcv::ImageMeta m_meta;
   
   // 2D PCA

   // Line Fit/Projected End
   std::vector<float> m_dir;
   cv::Point m_start;
   cv::Point m_end;
   void _fill_linefit_members();

   // Bounding Box (for collision detection)
   
 };
 

}

#endif
