#include "features.h"

NUMA* Top_Extremes(PIX *patch) {
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

NUMA* Right_Extremes(PIX *patch) {
    l_int32     patch_i, patch_j, patch_w, patch_h, wpl;
    l_uint32   *line, *data;
    l_float32  *array;
    /* NUMA *patch_na = pixCountByColumn(patch_current, NULL); */
    /* patch_na = numaTransform(patch_na, */
    /*                          0, */
    /*                          1.0/pixGetHeight(patch)); */
    NUMA *patch_na;
    pixGetDimensions(patch, &patch_w, &patch_h, NULL);
    patch_na = numaCreate(patch_h);
    numaSetCount(patch_na, patch_h);
    array = numaGetFArray(patch_na, L_NOCOPY);
    // Initializing
    for (patch_i = 0; patch_i < patch_h; patch_i++)
        array[patch_i] = -1;

    data = pixGetData(patch);
    wpl = pixGetWpl(patch);
    // We initialize the array with -1, and at first instant of
    // GET_DATA_BIT set it and avoid that width bin henceforth
    // NOTE: This can be done in better way.
    for (patch_i = 0; patch_i < patch_h; patch_i++) {
        line = data + wpl * patch_i;
        for (patch_j = (patch_w-1); patch_j > -1; patch_j--) {
            if ((array[patch_i] == -1) && GET_DATA_BIT(line, patch_j)) {
                array[patch_i] = patch_j/(l_float32)patch_w;
                break;
            }
        }
    }
    return patch_na;
}

NUMA* Left_Extremes(PIX *patch) {
    l_int32     patch_i, patch_j, patch_w, patch_h, wpl;
    l_uint32   *line, *data;
    l_float32  *array;
    /* NUMA *patch_na = pixCountByColumn(patch_current, NULL); */
    /* patch_na = numaTransform(patch_na, */
    /*                          0, */
    /*                          1.0/pixGetHeight(patch)); */
    NUMA *patch_na;
    pixGetDimensions(patch, &patch_w, &patch_h, NULL);
    patch_na = numaCreate(patch_h);
    numaSetCount(patch_na, patch_h);
    array = numaGetFArray(patch_na, L_NOCOPY);
    // Initializing
    for (patch_i = 0; patch_i < patch_h; patch_i++)
        array[patch_i] = -1;

    data = pixGetData(patch);
    wpl = pixGetWpl(patch);
    // We initialize the array with -1, and at first instant of
    // GET_DATA_BIT set it and avoid that width bin henceforth
    // NOTE: This can be done in better way.
    for (patch_i = 0; patch_i < patch_h; patch_i++) {
        line = data + wpl * patch_i;
        for (patch_j = 0; patch_j < patch_w; patch_j++) {
            if ((array[patch_i] == -1) && GET_DATA_BIT(line, patch_j)) {
                array[patch_i] = patch_j/(l_float32)patch_w;
                break;
            }
        }
    }
    return patch_na;
}

NUMA* Bottom_Extremes(PIX *patch) {
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

l_float32 DTW_Distance(NUMA* na_pix1, NUMA* na_pix2)
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
