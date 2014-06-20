// g++ -Wall -g src/word_features.cpp -I/usr/local/include/leptonica -I/usr/local/include/opencv2/ -L/usr/local/lib/ -llept -lopencv_core -lopencv_nonfree -lopencv_highgui -lopencv_features2d -lopencv_flann
#include <stdio.h>
#include <iostream>
#include "allheaders.h"
#include "opencv2/core/core.hpp"
#include "opencv2/nonfree/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/nonfree/nonfree.hpp"

using namespace cv;

static const l_int32  MIN_WORD_WIDTH = 5;
static const l_int32  MIN_WORD_HEIGHT = 10;
static const l_int32  MAX_WORD_WIDTH = 500;
static const l_int32  MAX_WORD_HEIGHT = 80;

void readme();

/** @function main */
int main( int argc, char** argv )
{
  if( argc != 3 )
      { readme(); return -1; }

  char *orig = argv[1];
  char *target = argv[2];
  static char  mainName[] = "Image Patch matching";
  PIX *pix_orig, *pix_target;
  if ((pix_orig = pixRead(orig)) == NULL)
      return ERROR_INT("pix_original image missing", mainName, 1);
  if ((pix_target = pixRead(target)) == NULL)
      return ERROR_INT("target image missing", mainName, 1);
    
  pix_orig = pixConvertRGBToGray(pix_orig, 0.33, 0.34, 0.33);
  pix_orig = pixThresholdToBinary(pix_orig, 128);

  pix_target = pixConvertRGBToGray(pix_target, 0.33, 0.34, 0.33);
  pix_target = pixThresholdToBinary(pix_target, 128);

  BOXA *boxa_current, *boxa_target;
  NUMA *nai_current, *nai_target;
  char patch_name[50];
  l_int32 i, j;

  pixGetWordBoxesInTextlines(pix_orig, 1, MIN_WORD_WIDTH, MIN_WORD_HEIGHT,
                             MAX_WORD_WIDTH, MAX_WORD_HEIGHT,
                             &boxa_current, &nai_current);
    
  pixGetWordBoxesInTextlines(pix_target, 1, MIN_WORD_WIDTH, MIN_WORD_HEIGHT,
                             MAX_WORD_WIDTH, MAX_WORD_HEIGHT,
                             &boxa_target, &nai_target);
  
  for (i = 0; i < boxa_current->n; i++) {
      BOX *box_current = boxaGetBox(boxa_current, i, L_CLONE);
      PIX *patch_current = pixClipRectangle(pix_orig, box_current, NULL);
      sprintf(patch_name, "source-%d.png", i);
      pixWrite(patch_name, patch_current, IFF_PNG);
      Mat img_source = imread( patch_name, CV_LOAD_IMAGE_GRAYSCALE );
      for (j = 0; j < boxa_target->n; j++) {
          BOX *box_target = boxaGetBox(boxa_target, j, L_CLONE);
          PIX *patch_target = pixClipRectangle(pix_target, box_target, NULL);
          sprintf(patch_name, "target-%d.png", j);
          pixWrite(patch_name, patch_target, IFF_PNG);
          Mat img_target = imread( patch_name, CV_LOAD_IMAGE_GRAYSCALE );
          if( !img_source.data || !img_target.data )
              { std::cout<< " --(!) Error reading images " << std::endl; return -1; }
          //-- Step 1: Detect the keypoints using SURF Detector
          int minHessian = 400;

          SurfFeatureDetector detector( minHessian );

          std::vector<KeyPoint> keypoints_object, keypoints_scene;

          detector.detect( img_source, keypoints_object );
          detector.detect( img_target, keypoints_scene );
          //-- Step 2: Calculate descriptors (feature vectors)
          SurfDescriptorExtractor extractor;

          Mat descriptors_object, descriptors_scene;

          extractor.compute( img_source, keypoints_object, descriptors_object );
          extractor.compute( img_target, keypoints_scene, descriptors_scene );

          //-- Step 3: Matching descriptor vectors using FLANN matcher
          FlannBasedMatcher matcher;
          std::vector< DMatch > matches;
          matcher.match( descriptors_object, descriptors_scene, matches );

          double max_dist = 0; double min_dist = 100;

          //-- Quick calculation of max and min distances between keypoints
          for( int i_dist = 0; i_dist < descriptors_object.rows; i_dist++ )
              { double dist = matches[i_dist].distance;
                  if( dist < min_dist ) min_dist = dist;
                  if( dist > max_dist ) max_dist = dist;
              }
          if ( max_dist > 0.0 && min_dist < 100.0) {
              printf("-- For Patch : %d & %d \n", i, j );
              printf("-- Max dist : %f \n", max_dist );
              printf("-- Min dist : %f \n", min_dist );
          }
      }          
  }

  return 0;
  }

/** @function readme */
void readme()
{ std::cout << " Usage: ./SURF_descriptor <img1> <img2>" << std::endl; }
