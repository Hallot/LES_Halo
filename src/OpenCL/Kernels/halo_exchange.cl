// From the received out buffers to the array
void exchange_2_halo_in(
    __global float2 *array,
    __global float2 *t_out,
    __global float2 *r_out,
    __global float2 *b_out,
    __global float2 *l_out
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_h
    ) {
    int i, j, k;
    int lcount = 0, rcount = 0;
    //t_out and b_out are contiguous on i and j, but not on k
    for (k = 0; k < km; k++) {
        memcpy(array.s0 + k * im * jm, t_out.s0, h_h * im * sizeof(*t_out));
        memcpy(array.s0 + (k + 1) * im * jm - (h_h * im), b_out.s0, h_h * im * sizeof(*b_out));
        memcpy(array.s1 + k * im * jm, t_out.s1, h_h * im * sizeof(*t_out));
        memcpy(array.s1 + (k + 1) * im * jm - (h_h * im), b_out.s1, h_h * im * sizeof(*b_out));
        
        // l_out and r_out are not contiguous
        for (i = 0; i < im; i++) {
            for (j = 0; j < jm; j++) {
                if (j > h_h - 1 && j < jm - h_h) {
                    // Left side of the array
                    if (i < h_h) {
                        array.s0[k * im * jm + j * im + i] = l_out.s0[lcount];
                        array.s1[k * im * jm + j * im + i] = l_out.s1[lcount];
                        lcount++;
                    }
                    // Right side of the array
                    if (i > im - h_h - 1) {
                        array.s0[k * im * jm + j * im + i] = r_out.s0[rcount];
                        array.s1[k * im * jm + j * im + i] = r_out.s1[rcount];
                        rcount++;
                    }
                }
            }
        }
    }
}

// From the array to the in buffers
void exchange_2_halo_out(
    __global float2 *array,
    __global float2 *t_in,
    __global float2 *r_in,
    __global float2 *b_in,
    __global float2 *l_in,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_h
    ) {
    int i, j, k;
    int lcount = 0, rcount = 0;
    //t_in and b_in are contiguous on i and j, but not on k
    for (k = 0; k < km; k++) {
        memcpy(t_in + k * im * h_h, array + k * im * jm, h_h * im * sizeof(*array));
        memcpy(b_in + k * im * h_h, array + (k + 1) * im * jm - (h_h * im), h_h * im * sizeof(*array));
        
        // l_out and r_out are not contiguous
        for (i = 0; i < im; i++) {
            for (j = 0; j < jm; j++) {
                // Left side of the array
                if (i < h_h) {
                    l_in[lcount] = array[k * im * jm + j * im + i];
                    lcount++;
                }
                // Right side of the array
                if (i > im - h_h - 1) {
                    r_in[rcount] = array[k * im * jm + j * im + i];
                    rcount++;
                }
            }
        }
    }
}