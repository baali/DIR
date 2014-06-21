/*
 *  alphabet_profiling.c
 *
 *     Comparing different alphabets and trying to identify
 *     differentiating shape features
 * 
 *  Compile: gcc -Wall -g -o alpha_profile src/alphabet_profiling.c
 *  src/features.c -I/usr/local/include/leptonica -I. -llept
 */

#include "string.h"
#include "features.h"

l_int32 main(int    argc,
             char **argv)
{
    char *file = argv[1];
    static char  mainName[] = "Alphabet profiling";
    PIX *pix;
    l_int32 i, j;
    char patch_name[50];

    if ((pix = pixRead(file)) == NULL)
        return ERROR_INT("pix_original image missing", mainName, 1);
    
    pix = pixConvertRGBToGray(pix, 0.33, 0.34, 0.33);
    pix = pixThresholdToBinary(pix, 128);
    BOXA *boxa_current;
    NUMA *nai_current;
    pixGetWordBoxesInTextlines(pix, 1, MIN_WORD_WIDTH, MIN_WORD_HEIGHT,
                               MAX_WORD_WIDTH, MAX_WORD_HEIGHT,
                               &boxa_current, &nai_current);
    fprintf(stderr, "Chracacters: %d\n", boxaGetCount(boxa_current));
    for (i = 0; i < boxa_current->n; i++) {
        // NUMA *DTW_distances;
        // DTW_distances = numaCreate(boxa_current->n);
        BOX *box_current = boxaGetBox(boxa_current, i, L_CLONE);
        PIX *patch_current = pixClipRectangle(pix, box_current, NULL);
        sprintf(patch_name, "bin/alpha-%d.png", i);
        pixWrite(patch_name, patch_current, IFF_PNG);
        NUMA *cpatch_roof = Get_Roof(patch_current);
        NUMA *cpatch_height = Get_Height(patch_current);
        NUMA *cpatch_floor = Get_Floor(patch_current);
        NUMA *cpatch_right = Right_Extremes(patch_current);
        NUMA *cpatch_left = Left_Extremes(patch_current);
        NUMA *cpatch_pix_row = pixCountByRow(patch_current, NULL);
        cpatch_pix_row = numaTransform(cpatch_pix_row, 0, 1/(l_float32)pixGetWidth(patch_current));
        for (j = 0; j < boxa_current->n; j++) {
            BOX *box_target = boxaGetBox(boxa_current, j, L_CLONE);
            PIX *patch_target = pixClipRectangle(pix, box_target, NULL);
            NUMA *tpatch_roof = Get_Roof(patch_target);
            NUMA *tpatch_height = Get_Height(patch_target);
            NUMA *tpatch_floor = Get_Floor(patch_target);
            NUMA *tpatch_right = Right_Extremes(patch_target);
            NUMA *tpatch_left = Left_Extremes(patch_target);            
            NUMA *tpatch_pix_row = pixCountByRow(patch_target, NULL);
            tpatch_pix_row = numaTransform(tpatch_pix_row, 0, 1/(l_float32)pixGetWidth(patch_target));
            l_float32 distance = 0;
            distance += DTW_Distance(cpatch_roof, tpatch_roof);
            distance += DTW_Distance(cpatch_height, tpatch_height);
            distance += DTW_Distance(cpatch_floor, tpatch_floor);
            distance += DTW_Distance(cpatch_pix_row, tpatch_pix_row);
            distance += DTW_Distance(cpatch_left, tpatch_left);
            distance += DTW_Distance(cpatch_right, tpatch_right);
            fprintf(stderr, "%.1f, ", distance);
        }
        fprintf(stderr, "\n");
    }
    return 0;
}
