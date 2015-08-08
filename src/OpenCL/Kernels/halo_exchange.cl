// Write the recieved buffer to the array
// Halo order in the buffer: top|bottom|left|right where top and bottom contain the corners and thus have two more points along i each
void exchange_2_halo_write(
    __global float2 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km
    ) {
    const unsigned int v_dim = 2;
    const unsigned int buf_sz = v_dim * 4 * (im+1) * km;
    const unsigned int v_limit = buf_sz / v_dim;
    const unsigned int tp_bound = (buf_sz / 4 * v_dim) + 2;
    const unsigned int bl_bound = (buf_sz / 2 * v_dim) + 3;
    const unsigned int lr_bound = (3 * buf_sz / 4 * v_dim) + 2;
    const unsigned int i;
    // Iterate along buffer
    for (i = 0; i < buf_sz; i++) {
        // Which vector component, ie along v_dim
        if (i < v_limit) {
            // top halo
            if (i < tp_bound) {
                // Can't simplfify im because it relies on integer division!
                array.s0[(i%im) + (i/im)*(im*jm)] = buffer[i];
            }
            // bottom halo
            if (i >= tp_bound && i < bl_bound) {
                // Can't simplfify im because it relies on integer division!
                array.s0[((i-tp_bound)%im) + im*(jm-1) + ((i-tp_bound)/im)*(im*jm)] = buffer[i];
            }
            // left halo
            if (i >= bl_bound && i < lr_bound) {
                array.s0[2*im*((i-bl_bound)/(jm-2)) * ((i-bl_bound)+1)*im] = buffer[i];
            }
            // right halo
            if (i >= lr_bound) {
                array.s0[2*im*((i-lr_bound)/(jm-2)) * ((i-lr_bound)+1)*im + (im-1)] = buffer[i];
            }
        }
        else if (i > v_limit) {
            // top halo
            if (i < tp_bound) {
                // Can't simplfify im because it relies on integer division!
                array.s1[(i%im) + (i/im)*(im*jm)] = buffer[i];
            }
            // bottom halo
            if (i >= tp_bound && i < bl_bound) {
                // Can't simplfify im because it relies on integer division!
                array.s1[((i-tp_bound)%im) + im*(jm-1) + ((i-tp_bound)/im)*(im*jm)] = buffer[i];
            }
            // left halo
            if (i >= bl_bound && i < lr_bound) {
                array.s1[2*im*((i-bl_bound)/(jm-2)) * ((i-bl_bound)+1)*im] = buffer[i];
            }
            // right halo
            if (i >= lr_bound) {
                array.s1[2*im*((i-lr_bound)/(jm-2)) * ((i-lr_bound)+1)*im + (im-1)] = buffer[i];
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