/* File: features.h */
#ifndef FEATURES_H
#define FEATURES_H
#include "allheaders.h"

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) < (y) ? (y) : (x))
#define min3(x, y, z) (min(x, min(y, z)))
#define argmin3(x, y, z) ((x < y && x < z) ? 0 : (y < z ? 1 : 2))
#define sgn(x) ((x) < 0 ? -1 : 1)
#define abs(x) ((x)*sgn(x))
#define unused(x) ((void) x)

static const l_int32  MIN_WORD_WIDTH = 1;
static const l_int32  MIN_WORD_HEIGHT = 10;
static const l_int32  MAX_WORD_WIDTH = 500;
static const l_int32  MAX_WORD_HEIGHT = 80;

NUMA* Get_Roof(PIX *);
NUMA* Get_Floor(PIX *);
NUMA* Get_Height(PIX *);
NUMA* Left_Extremes(PIX *);
NUMA* Right_Extremes(PIX *);
l_float32 DTW_Distance(NUMA*, NUMA*);
    
#endif
