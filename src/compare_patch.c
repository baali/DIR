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
l_int32 NUM_NEIGHBOR_PATCHES = 1;
static const l_int32  MIN_WORD_WIDTH = 5;
static const l_int32  MIN_WORD_HEIGHT = 10;
static const l_int32  MAX_WORD_WIDTH = 500;
static const l_int32  MAX_WORD_HEIGHT = 80;

NUMA* Get_Roof(PIX *patch) {
    l_int32     patch_i, patch_j, patch_w, patch_h, wpl;
    l_uint32   *line, *data;
    l_float32  *array;
    /* NUMA *patch_na = pixCountByColumn(patch_current, NULL); */
    /* patch_na = numaTransform(patch_na, */
    /*                          0, */
    /*                          1.0/pixGetHeight(patch)); */
    NUMA *patch_na;
    pixGetDimensions(patch, &patch_w, &patch_h, NULL);
    patch_na = numaCreate(patch_w);
    numaSetCount(patch_na, patch_w);
    array = numaGetFArray(patch_na, L_NOCOPY);
    // Initializing
    for (patch_j = 0; patch_j < patch_w; patch_j++)
        array[patch_j] = -1;

    data = pixGetData(patch);
    wpl = pixGetWpl(patch);
    // We initialize the array with -1, and at first instant of
    // GET_DATA_BIT set it and avoid that width bin henceforth
    // NOTE: This can be done in better way.
    for (patch_i = 0; patch_i < patch_h; patch_i++) {
        line = data + wpl * patch_i;
        for (patch_j = 0; patch_j < patch_w; patch_j++) {
            if ((array[patch_j] == -1) && GET_DATA_BIT(line, patch_j)) {
                array[patch_j] = patch_i/(l_float32)patch_h;
            }
        }
    }
    return patch_na;
}

NUMA* Get_Floor(PIX *patch) {
    l_int32     patch_i, patch_j, patch_w, patch_h, wpl;
    l_uint32   *line, *data;
    l_float32  *array;
    /* NUMA *patch_na = pixCountByColumn(patch_current, NULL); */
    /* patch_na = numaTransform(patch_na, */
    /*                          0, */
    /*                          1.0/pixGetHeight(patch)); */
    NUMA *patch_na;
    pixGetDimensions(patch, &patch_w, &patch_h, NULL);
    patch_na = numaCreate(patch_w);
    numaSetCount(patch_na, patch_w);
    array = numaGetFArray(patch_na, L_NOCOPY);
    // Initializing
    for (patch_j = 0; patch_j < patch_w; patch_j++)
        array[patch_j] = -1;

    data = pixGetData(patch);
    wpl = pixGetWpl(patch);
    // We initialize the array with -1, and at first instant of
    // GET_DATA_BIT set it and avoid that width bin henceforth
    // NOTE: This can be done in better way.
    for (patch_i = (patch_h - 1); patch_i > -1 ; patch_i--) {
        line = data + wpl * patch_i;
        for (patch_j = 0; patch_j < patch_w; patch_j++) {
            if ((array[patch_j] == -1) && GET_DATA_BIT(line, patch_j)) {
                array[patch_j] = patch_i/(l_float32)patch_h;
            }
        }
    }
    return patch_na;
}

NUMA* Get_Height(PIX *patch) {
    l_int32 patch_i, patch_j, patch_w, patch_h, wpl;
    l_uint32 *line, *data;
    l_float32 *array_top, *array_bottom;
    /* NUMA *patch_na = pixCountByColumn(patch_current, NULL); */
    /* patch_na = numaTransform(patch_na, */
    /*                          0, */
    /*                          1.0/pixGetHeight(patch)); */
    NUMA *patch_top;
    NUMA *patch_bottom;
    NUMA *patch_height;
    pixGetDimensions(patch, &patch_w, &patch_h, NULL);
    patch_top = numaCreate(patch_w);
    patch_bottom = numaCreate(patch_w);
    patch_height = numaCreate(patch_w);
    numaSetCount(patch_top, patch_w);
    numaSetCount(patch_bottom, patch_w);
    numaSetCount(patch_height, patch_w);
    array_top = numaGetFArray(patch_top, L_NOCOPY);
    // Initializing
    for (patch_j = 0; patch_j < patch_w; patch_j++)
        array_top[patch_j] = -1;

    data = pixGetData(patch);
    wpl = pixGetWpl(patch);
    // We initialize the array with -1, and at first instant of
    // GET_DATA_BIT set it and avoid that width bin henceforth
    // NOTE: This can be done in better way.
    for (patch_i = 0; patch_i < patch_h; patch_i++) {
        line = data + wpl * patch_i;
        for (patch_j = 0; patch_j < patch_w; patch_j++) {
            if ((array_top[patch_j] == -1) && GET_DATA_BIT(line, patch_j)) {
                array_top[patch_j] = patch_i/(l_float32)patch_h;
            }
        }
    }
    array_bottom = numaGetFArray(patch_bottom, L_NOCOPY);
    // Initializing
    for (patch_j = 0; patch_j < patch_w; patch_j++)
        array_bottom[patch_j] = -1;
    
    for (patch_i = (patch_h - 1); patch_i > -1 ; patch_i--) {
        line = data + wpl * patch_i;
        for (patch_j = 0; patch_j < patch_w; patch_j++) {
            if ((array_bottom[patch_j] == -1) && GET_DATA_BIT(line, patch_j)) {
                array_bottom[patch_j] = patch_i/(l_float32)patch_h;
            }
        }
    }
    
    patch_height = numaArithOp(NULL, patch_top, patch_bottom, L_ARITH_SUBTRACT);
    return patch_height;
}

l_float32 DTW_distance(NUMA* na_pix1, NUMA* na_pix2)
{
    // Function to calculate DTW of two na sequences
    l_int32 window = max(50, abs(na_pix1->n - na_pix2->n));
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
        NUMA *DTW_distances;
        DTW_distances = numaCreate(boxa_target->n);
        BOX *box_current = boxaGetBox(boxa_current, i, L_CLONE);
        PIX *patch_current = pixClipRectangle(pix_orig, box_current, NULL);
        NUMA *cpatch_roof = Get_Roof(patch_current);
        NUMA *cpatch_height = Get_Height(patch_current);
        NUMA *cpatch_floor = Get_Floor(patch_current);
        NUMA *cpatch_pix_row = pixCountByRow(patch_current, NULL);
        cpatch_pix_row = numaTransform(cpatch_pix_row, 0, 1/(l_float32)pixGetWidth(patch_current));
        // Now we have to run DTW over these numa series of current
        // patch and target patch
        for (j = 0; j < boxa_target->n; j++) {
            BOX *box_target = boxaGetBox(boxa_target, j, L_CLONE);
            PIX *patch_target = pixClipRectangle(pix_target, box_target, NULL);
            NUMA *tpatch_roof = Get_Roof(patch_target);
            NUMA *tpatch_height = Get_Height(patch_target);
            NUMA *tpatch_floor = Get_Floor(patch_target);
            NUMA *tpatch_pix_row = pixCountByRow(patch_target, NULL);
            tpatch_pix_row = numaTransform(tpatch_pix_row, 0, 1/(l_float32)pixGetWidth(patch_target));
            // DTW between cpatch_na and tpatch_na
            l_float32 distance = 0;
            /* distance += DTW_distance(cpatch_roof, tpatch_roof); */
            /* distance += DTW_distance(cpatch_height, tpatch_height); */
            /* distance += DTW_distance(cpatch_floor, tpatch_floor); */
            distance += DTW_distance(cpatch_pix_row, tpatch_pix_row);
            // fprintf(stderr, "Distance %f with patch %d\n", distance, j);
            numaAddNumber(DTW_distances, distance);
            boxDestroy(&box_target);
            pixDestroy(&patch_target);
            numaDestroy(&tpatch_roof);
            numaDestroy(&tpatch_height);
            numaDestroy(&tpatch_floor);
        }

        // fprintf(stderr, "Done calculating distances\n");
        l_float32 min_distance;
        l_int32 location;
        numaGetMin(DTW_distances, &min_distance, &location);
        fprintf(stderr, "For %d, mininum distance %f with %d\n", i, min_distance, location);
        BOX *result_box = boxaGetBox(boxa_target, location, L_CLONE);
        PIX *result_pix = pixClipRectangle(pix_target, result_box, NULL);
        sprintf(patch_name, "result-%d-%d.png", i, location);
        pixWrite(patch_name, result_pix, IFF_PNG);
        sprintf(patch_name, "patch-%d.png", i);
        pixWrite(patch_name, patch_current, IFF_PNG);
        numaDestroy(&DTW_distances);
        numaDestroy(&cpatch_floor);
        numaDestroy(&cpatch_height);
        numaDestroy(&cpatch_roof);
        pixDestroy(&result_pix);
        boxDestroy(&result_box);
        pixDestroy(&patch_current);
        boxDestroy(&box_current);
    }

    /* fprintf(stderr, "Average Shift %d, %d\n", */
    /*         shift_x/neighbor_patches->n, */
    /*         shift_y/neighbor_patches->n); */
    /* pixWrite("orig-patch.png", pix_orig, IFF_PNG); */
    /* pixWrite("target-patch.png", pix_target, IFF_PNG); */
    boxaDestroy(&boxa_target);
    boxaDestroy(&boxa_current);
    numaDestroy(&nai_target);
    numaDestroy(&nai_current);
    pixDestroy(&pix_orig);
    pixDestroy(&pix_target);
    return 0;
}
