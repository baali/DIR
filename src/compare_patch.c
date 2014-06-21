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
#include "features.h"

// number of neighbor patches to locate on both sides
l_int32 NUM_NEIGHBOR_PATCHES = 1;

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
        NUMA *cpatch_right = Right_Extremes(patch_current);
        NUMA *cpatch_left = Left_Extremes(patch_current);
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
            NUMA *tpatch_right = Right_Extremes(patch_target);
            NUMA *tpatch_left = Left_Extremes(patch_target);            
            NUMA *tpatch_pix_row = pixCountByRow(patch_target, NULL);
            tpatch_pix_row = numaTransform(tpatch_pix_row, 0, 1/(l_float32)pixGetWidth(patch_target));
            // DTW between cpatch_na and tpatch_na
            l_float32 distance = 0;
            distance += DTW_Distance(cpatch_roof, tpatch_roof);
            distance += DTW_Distance(cpatch_height, tpatch_height);
            distance += DTW_Distance(cpatch_floor, tpatch_floor);
            distance += DTW_Distance(cpatch_pix_row, tpatch_pix_row);
            distance += DTW_Distance(cpatch_left, tpatch_left);
            distance += DTW_Distance(cpatch_right, tpatch_right);
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
        sprintf(patch_name, "bin/result-%d-%d.png", i, location);
        pixWrite(patch_name, result_pix, IFF_PNG);
        sprintf(patch_name, "bin/patch-%d.png", i);
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
