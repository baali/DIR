/*
 *  compare_patch.c
 *
 *     based on reference image's patch locate the same patch in other
 *     images.
 * 
 *  Compile: gcc -Wall -g src/compare_patch.c
 *  -I/usr/local/include/leptonica -llept
 */

#include "string.h"
#include "allheaders.h"

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) < (y) ? (y) : (x))
#define min3(x, y, z) (min(x, min(y, z)))
#define argmin3(x, y, z) ((x < y && x < z) ? 0 : (y < z ? 1 : 2))
#define sgn(x) ((x) < 0 ? -1 : 1)
#define abs(x) ((x)*sgn(x))
#define unused(x) ((void) x)

// number of neighbor patches to locate on both sides
l_int32 NUM_NEIGHBOR_PATCHS = 1;
static const l_int32  MIN_WORD_WIDTH = 5;
static const l_int32  MIN_WORD_HEIGHT = 10;
static const l_int32  MAX_WORD_WIDTH = 500;
static const l_int32  MAX_WORD_HEIGHT = 80;

l_float32 DTW_distance(NUMA* na_pix1, NUMA* na_pix2)
{
    // Function to calculate DTW of two na sequences
    l_int32 window = max(10, abs(na_pix1->n - na_pix2->n));
    l_float32 *DTW = (l_float32 *)lept_calloc((na_pix1->n)*(na_pix2->n),
                                              sizeof(l_float32));

    l_int32 i, j;
    for (i=0; i<na_pix1->n; i++) {
        for (j=0; j<na_pix2->n; j++) {
            DTW[(i*na_pix2->n)+j] = 100000;
        }
    }
    DTW[0] = 0;
    for (i=1; i<na_pix1->n; i++) {
        for (j=max(1, i-window); j<min(na_pix2->n, i+window); j++) {
        // for (j=1; j<na_pix2->n; j++) {
            l_int32 s, d;
            numaGetIValue(na_pix1, i, &s);
            numaGetIValue(na_pix2, j, &d);
            l_int32 distance = abs(s - d);
            DTW[(i*na_pix2->n)+j] = distance +  \
                min3(DTW[((i-1)*na_pix2->n)+j],
                     DTW[(i*na_pix2->n)+(j-1)],
                     DTW[((i-1)*na_pix2->n)+(j-1)]);
        }
    }
    return DTW[(na_pix1->n)*(na_pix2->n)-1];
}

l_int32 main(int    argc,
             char **argv)
{
    // static char orig[] = "current-patch-1.png";
    // static char shifted[] = "/home/baali/lancing/LouisTully/data/FAX_20140205_1391633868_46_1_1.jpg";
    char *orig = argv[1];
    char *shifted = argv[2];
    static char  mainName[] = "Image Patch matching";
    PIX *pix_orig, *pix_shifted;
    if ((pix_orig = pixRead(orig)) == NULL)
        return ERROR_INT("pix_original image missing", mainName, 1);
    if ((pix_shifted = pixRead(shifted)) == NULL)
        return ERROR_INT("shifted image missing", mainName, 1);
    
    /* pix_orig = pixConvertRGBToGray(pix_orig, 0.33, 0.34, 0.33); */
    /* pix_orig = pixThresholdToBinary(pix_orig, 128); */

    pix_shifted = pixConvertRGBToGray(pix_shifted, 0.33, 0.34, 0.33);
    pix_shifted = pixThresholdToBinary(pix_shifted, 128);

    BOXA *boxa_shifted;
    NUMA *nai_shifted;
    // char patch_name[50];
    int j;
    
    pixGetWordBoxesInTextlines(pix_shifted, 1, MIN_WORD_WIDTH, MIN_WORD_HEIGHT,
                               MAX_WORD_WIDTH, MAX_WORD_HEIGHT,
                               &boxa_shifted, &nai_shifted);

    NUMA *DTW_distances;
    DTW_distances = numaCreate(boxa_shifted->n);
    
    l_int32     cpatch_i, cpatch_j, cpatch_w, cpatch_h, cwpl;
    l_uint32   *cline, *cdata;
    l_float32  *carray;
    // NUMA *cpatch_na = pixCountByColumn(pix_orig, NULL);
    NUMA *cpatch_na;
    pixGetDimensions(pix_orig, &cpatch_w, &cpatch_h, NULL);
    cpatch_na = numaCreate(cpatch_w);
    numaSetCount(cpatch_na, cpatch_w);
    carray = numaGetFArray(cpatch_na, L_NOCOPY);
    // Initializing
    for (cpatch_j = 0; cpatch_j < cpatch_w; cpatch_j++)
        carray[cpatch_j] = -1;

    cdata = pixGetData(pix_orig);
    cwpl = pixGetWpl(pix_orig);
    // We initialize the array with -1, and at first instant of
    // GET_DATA_BIT set it and avoid that width bin henceforth
    // NOTE: This can be done in better way.
    for (cpatch_i = 0; cpatch_i < cpatch_h; cpatch_i++) {
        cline = cdata + cwpl * cpatch_i;
        for (cpatch_j = 0; cpatch_j < cpatch_w; cpatch_j++) {
            if ((carray[cpatch_j] == -1) && GET_DATA_BIT(cline, cpatch_j)) {
                carray[cpatch_j] = cpatch_i/(l_float32)cpatch_h;
            }
        }
    }

    fprintf(stderr, "Found %d word blocks\n", boxa_shifted->n);
    // Now we have to run DTW over these numa series of current
    // patch and target patch
    for (j = 0; j < boxa_shifted->n; j++) {
        BOX *box_shifted = boxaGetBox(boxa_shifted, j, L_CLONE);
        PIX *patch_shifted = pixClipRectangle(pix_shifted, box_shifted, NULL);

        // Calculating roof of target patch
        l_int32     tpatch_i, tpatch_j, tpatch_w, tpatch_h, twpl;
        l_uint32   *tline, *tdata;
        l_float32  *tarray;
        // NUMA *tpatch_na = pixCountByColumn(patch_shifted, NULL);
        NUMA *tpatch_na;
        pixGetDimensions(patch_shifted, &tpatch_w, &tpatch_h, NULL);
        tpatch_na = numaCreate(tpatch_w);
        numaSetCount(tpatch_na, tpatch_w);
        tarray = numaGetFArray(tpatch_na, L_NOCOPY);
        // Initializing
        for (tpatch_j = 0; tpatch_j < tpatch_w; tpatch_j++)
            tarray[tpatch_j] = -1;
        
        tdata = pixGetData(patch_shifted);
        twpl = pixGetWpl(patch_shifted);
        // We initialize the array with -1, and at first instant of
        // GET_DATA_BIT set it and avoid that width bin henceforth
        // NOTE: This can be done in better way.
        for (tpatch_i = 0; tpatch_i < tpatch_h; tpatch_i++) {
            tline = tdata + twpl * tpatch_i;
            for (tpatch_j = 0; tpatch_j < tpatch_w; tpatch_j++) {
                if ((tarray[tpatch_j] == -1) && GET_DATA_BIT(tline, tpatch_j)) {
                    // Normalizing the height/roof value
                    tarray[tpatch_j] = tpatch_i/(l_float32)tpatch_h;
                }
            }
        }
        // DTW between cpatch_na and tpatch_na
        l_float32 distance = DTW_distance(cpatch_na, tpatch_na);
        fprintf(stderr, "Distance %f with patch %d\n", distance, j);
        numaAddNumber(DTW_distances, distance);
        boxDestroy(&box_shifted);
        pixDestroy(&patch_shifted);
        numaDestroy(&tpatch_na);
    }

    // fprintf(stderr, "Done calculating distances\n");
    l_float32 min_distance;
    l_int32 location;
    numaGetMin(DTW_distances, &min_distance, &location);
    fprintf(stderr, "Mininum Distance %f with patch %d\n", min_distance, location);
    BOX *target_box = boxaGetBox(boxa_shifted, location, L_CLONE);
    PIX *target_pix = pixClipRectangle(pix_shifted, target_box, NULL);
    pixWrite("matched-patch.png", target_pix, IFF_PNG);
    /* fprintf(stderr, "Average Shift %d, %d\n", */
    /*         shift_x/neighbor_patches->n, */
    /*         shift_y/neighbor_patches->n); */
    /* pixWrite("orig-patch.png", pix_orig, IFF_PNG); */
    /* pixWrite("shifted-patch.png", pix_shifted, IFF_PNG); */
    
    numaDestroy(&DTW_distances);
    numaDestroy(&cpatch_na);
    boxaDestroy(&boxa_shifted);
    numaDestroy(&nai_shifted);
    pixDestroy(&pix_orig);
    pixDestroy(&pix_shifted);
    return 0;
}
