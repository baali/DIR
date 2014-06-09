/*
 *  compare_patch.c
 *
 *     based on reference image's patch locate the same patch in other
 *     images.
 * 
 *  Compile: gcc src/compare_patch.c -I/usr/local/include/leptonica
 *  -llept
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

// patch location of total amount for reference image
// FAX_20140423_1398285784_751.jpg
l_int32 x = 2330;
l_int32 y = 1020;
l_int32 w = 180;
l_int32 h = 80;

l_int32 DTW_distance(NUMA* na_pix1, NUMA* na_pix2)
{
    // Function to calculate DTW of two na sequences
    l_int32 window = abs(na_pix1->n - na_pix2->n);
    if (window > min(na_pix1->n, na_pix2->n)/2)
        return -1;
    fprintf(stderr, "size: (%d, %d), window: %d\n", na_pix1->n,
            na_pix2->n,
            window);
    l_int32 *DTW = (l_int32 *)lept_calloc((na_pix1->n)*(na_pix2->n),
                                          sizeof(l_int32));
    l_int32 i, j;
    for (i=1; i<na_pix1->n; i++) {
        DTW[(i*na_pix2->n)] = 10000;
    }
    for (j=1; j<na_pix2->n; j++) {
        DTW[j] = 10000;
    }
    DTW[0] = 0;
    fprintf(stderr, "Done initializing\n");
    for (i=1; i<na_pix1->n; i++) {
        // for (j=max(1, i-window); j<min(na_pix2->n, i+window); j++) {
        for (j=1; j<na_pix2->n; j++) {
            l_int32 s, d;
            numaGetIValue(na_pix1, i, &s);
            numaGetIValue(na_pix2, j, &d);
            l_int32 distance = abs(s - d);
            // fprintf(stderr, "(s, d): (%d, %d)\n", s, d);
            DTW[(i*na_pix2->n)+j] = distance +  \
                min3(DTW[((i-1)*na_pix2->n)+j],
                     DTW[(i*na_pix2->n)+(j-1)],
                     DTW[((i-1)*na_pix2->n)+(j-1)]);
        }
    }
    fprintf(stderr, "Done Calculating DTW\n");
    return DTW[(na_pix1->n)*(na_pix2->n)-1];
}

l_int32 main(int    argc,
             char **argv)
{
    static char orig[] = "FAX_20140423_1398285784_751.jpg";
    static char shifted[] = "FAX_20140423_1398285784_751-shift.jpg";
    static char  mainName[] = "Image Patch matching";
    PIX *pix_orig, *pix_shifted;
    if ((pix_orig = pixRead(orig)) == NULL)
        return ERROR_INT("pix_original image missing", mainName, 1);
    if ((pix_shifted = pixRead(shifted)) == NULL)
        return ERROR_INT("shifted image missing", mainName, 1);
    
    PIX *orig_gray = pixConvertRGBToGray(pix_orig, 0.33, 0.34, 0.33);
    orig_gray = pixThresholdToBinary(orig_gray, 128);

    PIX *shifted_gray = pixConvertRGBToGray(pix_shifted, 0.33, 0.34, 0.33);
    shifted_gray = pixThresholdToBinary(shifted_gray, 128);

    BOXA *boxa_orig, *boxa_shifted;
    NUMA *nai_orig, *nai_shifted;
    char patch_name[50];
    int i, j, count;
    
    pixGetWordBoxesInTextlines(orig_gray, 1, MIN_WORD_WIDTH, MIN_WORD_HEIGHT,
                               MAX_WORD_WIDTH, MAX_WORD_HEIGHT,
                               &boxa_orig, &nai_orig);

    // Finding boxa matching our highlighted patch
    BOX *patch_box = boxaGetNearestToPt(boxa_orig, x+(w/2), y+(h/2));
    PIX *patch_pix = pixClipRectangle(pix_orig, patch_box, NULL); 
    sprintf(patch_name, "patch.png");
    pixWrite(patch_name, patch_pix, IFF_PNG);

    // Getting index of the box found
    l_int32 patch_index;
    for (i = 0; i < boxa_orig->n; i++) {
        BOX *box_orig = boxaGetBox(boxa_orig, i, L_CLONE);
        l_int32 local_same;
        boxEqual(box_orig, patch_box, &local_same);
        if (local_same) {
            patch_index = i;
            break;
        }
        boxDestroy(&box_orig);
    }
    fprintf(stderr, "Patch index: %d\n", patch_index);

    BOXA *neighbor_patches = boxaCreate(2*NUM_NEIGHBOR_PATCHS);
    for (i = 1; i <= NUM_NEIGHBOR_PATCHS; i++) {
        BOX *temp_box = boxaGetBox(boxa_orig, patch_index - i, L_CLONE);
        boxaAddBox(neighbor_patches, temp_box, L_CLONE);
        boxDestroy(&temp_box);
        temp_box = boxaGetBox(boxa_orig, patch_index + i, L_CLONE);
        boxaAddBox(neighbor_patches, temp_box, L_CLONE);
        boxDestroy(&temp_box);
    }
    
    pixGetWordBoxesInTextlines(shifted_gray, 1, MIN_WORD_WIDTH, MIN_WORD_HEIGHT,
                               MAX_WORD_WIDTH, MAX_WORD_HEIGHT,
                               &boxa_shifted, &nai_shifted);

    PIXCMAP *cmap = pixcmapCreateRandom(8, 1, 1);
    l_int32 shift_x = 0, shift_y = 0;

    for (i = 0; i < neighbor_patches->n; i++) {
        BOX *box_patch = boxaGetBox(neighbor_patches, i, L_CLONE);
        PIX *pix_patch = pixClipRectangle(orig_gray, box_patch, NULL);

        l_int32     cpatch_i, cpatch_j, cpatch_w, cpatch_h, cwpl;
        l_uint32   *cline, *cdata;
        l_float32  *carray;
        NUMA *cpatch_na;
        pixGetDimensions(pix_patch, &cpatch_w, &cpatch_h, NULL);        
        cpatch_na = numaCreate(cpatch_w);
        numaSetCount(cpatch_na, cpatch_w);
        carray = numaGetFArray(cpatch_na, L_NOCOPY);
        // Initializing
        for (cpatch_j = 0; cpatch_j < cpatch_w; cpatch_j++)
            carray[cpatch_j] = -1;
        
        cdata = pixGetData(pix_patch);
        cwpl = pixGetWpl(pix_patch);
        // We initialize the array with -1, and at first instant of
        // GET_DATA_BIT set it and avoid that width bin henceforth
        // NOTE: This can be done in better way.
        for (cpatch_i = 0; cpatch_i < cpatch_h; cpatch_i++) {
            cline = cdata + cwpl * cpatch_i;
            for (cpatch_j = 0; cpatch_j < cpatch_w; cpatch_j++) {
                if ((carray[cpatch_j] == -1) && GET_DATA_BIT(cline, cpatch_j)) {
                    carray[cpatch_j] = cpatch_i;
                }
            }
        }
        // Now we have to run DTW over these numa series of current
        // patch and target patch
        sprintf(patch_name, "current-patch-%d.png", i);
        pixWrite(patch_name, pix_patch, IFF_PNG);
        for (j = 0; j < boxa_shifted->n; j++) {
            BOX *box_shifted = boxaGetBox(boxa_shifted, j, L_CLONE);
            PIX *patch_shifted = pixClipRectangle(shifted_gray, box_shifted, NULL);

            // Calculating roof of target patch
            l_int32     tpatch_i, tpatch_j, tpatch_w, tpatch_h, twpl;
            l_uint32   *tline, *tdata;
            l_float32  *tarray;
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
                        tarray[tpatch_j] = tpatch_i;
                    }
                }
            }
            // DTW between cpatch_na and tpatch_na
            l_int32 distance = DTW_distance(cpatch_na, tpatch_na);
            if (distance != -1) 
                fprintf(stderr, "Distance (%d, %d) => %d\n", i, j, distance);
            boxDestroy(&box_shifted);
            pixDestroy(&pix_shifted);
        }
    }

    fprintf(stderr, "Average Shift %d, %d\n",
            shift_x/neighbor_patches->n,
            shift_y/neighbor_patches->n);
    /* pixWrite("orig-patch.png", pix_orig, IFF_PNG); */
    /* pixWrite("shifted-patch.png", pix_shifted, IFF_PNG); */
    boxaDestroy(&boxa_orig);
    boxaDestroy(&boxa_shifted);
    numaDestroy(&nai_orig);
    numaDestroy(&nai_shifted);
    pixDestroy(&pix_orig);
    pixDestroy(&pix_shifted);
    return 0;
}
