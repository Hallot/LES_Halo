// Write the recieved buffer to the array
void exchange_2_halo_write(
    __global float2 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km
    ) {
    int v_dim = 2;
    int i, j, k;
    
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