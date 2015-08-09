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
    const unsigned int buf_sz = v_dim * 2 * km * (im + jm - 2);
    const unsigned int v_limit = buf_sz / v_dim;
    const unsigned int tp_bound = im * km;
    const unsigned int bl_bound = 2 * im * km;
    const unsigned int lr_bound = 2 * im * km + (jm - 2) * km;
    float *vector[v_dim];
    unsigned int i, i_off, vec_off;
    
    vector[0] = array.s0;
    vector[1] = array.s1;
    
    // Iterate along buffer
    for (i = 0; i < buf_sz; i++) {
        // Which vector component, ie along v_dim
        vec_off = i / v_limit;
        // Offset for each vector
        i_off = i - (vec_off * v_limit);
        // top halo
        if (i_off < tp_bound) {
            // Can't simplify im because it relies on integer division!
            vector[vec_off][(i_off%im) + (i_off/im)*(im*jm)] = buffer[i];
        }
        // bottom halo
        if (i_off >= tp_bound && i_off < bl_bound) {
            // Can't simplify im because it relies on integer division!
            vector[vec_off][((i_off-tp_bound)%im) + im*(jm-1) + ((i_off-tp_bound)/im)*(im*jm)] = buffer[i];
        }
        // left halo
        if (i_off >= bl_bound && i_off < lr_bound) {
            vector[vec_off][2*im*((i_off-bl_bound)/(jm-2)) + ((i_off-bl_bound)+1)*im] = buffer[i];
        }
        // right halo
        if (i_off >= lr_bound) {
            vector[vec_off][2*im*((i_off-lr_bound)/(jm-2)) + ((i_off-lr_bound)+1)*im + (im-1)] = buffer[i];
        }
    }
}

// Write the data from the array into the buffer
// Halo order in the buffer: top|bottom|left|right where top and bottom do not contain the corners and thus have the same size as left and right
void exchange_2_halo_read(
    __global float2 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km
    ) {
    const unsigned int v_dim = 2;
    const unsigned int buf_sz = v_dim * 2 * km * (im + jm - 2);
    const unsigned int v_limit = buf_sz / v_dim;
    const unsigned int tp_bound = im * km;
    const unsigned int bl_bound = 2 * im * km;
    const unsigned int lr_bound = 2 * im * km + (jm - 2) * km;
    float *vector[v_dim];
    unsigned int i, i_off, vec_off;
    unsigned int i_buf = 0;
    
    vector[0] = array.s0;
    vector[1] = array.s1;
    
    // Iterate along buffer
    for (i = 0; i < buf_sz; i++) {
        // Which vector component, ie along v_dim
        vec_off = i / v_limit;
        // Offset for each vector
        i_off = i - (vec_off * v_limit);
        // top halo
        if (i_off < tp_bound) {
            // We don't need the first and last element when reading
            if (i_off%im == 0) {
                continue;
            }
            // Can't simplify im because it relies on integer division!
            buffer[i_buf] = vector[vec_off][(i_off%im) + (i_off/im)*(im*jm)];
            i_buf++;
        }
        // bottom halo
        if (i_off >= tp_bound && i_off < bl_bound) {
            // We don't need the first and last element when reading
            if (i_off%im == 0) {
                continue;
            }
            // Can't simplify im because it relies on integer division!
            buffer[i_buf] = vector[vec_off][((i_off-tp_bound)%im) + im*(jm-1) + ((i_off-tp_bound)/im)*(im*jm)];
            i_buf++;
        }
        // left halo
        if (i_off >= bl_bound && i_off < lr_bound) {
            buffer[i_buf] = vector[vec_off][2*im*((i_off-bl_bound)/(jm-2)) + ((i_off-bl_bound)+1)*im + 1];
            i_buf++;
        }
        // right halo
        if (i_off >= lr_bound) {
            buffer[i_buf] = vector[vec_off][2*im*((i_off-lr_bound)/(jm-2)) + ((i_off-lr_bound)+1)*im + (im-1) - 1];
            i_buf++;
        }
    }
}
