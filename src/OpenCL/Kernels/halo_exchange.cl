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
    const unsigned int k = get_global_id(0);
    const unsigned int i = get_global_id(1);
    const unsigned int v_sz = 2*km*im + 2*km*(jm-2);

    // Which vector component, ie along v_dim
    for (unsigned int v = 0; v < v_dim; v++) {
        if (i < im && k < km){
            // top halo
            ((__global float*)&array[i + k*im*jm])[v] = buffer[v*v_sz + i + k*im];
            // bottom halo
            ((__global float*)&array[i + k*im*jm + im*(jm-1)])[v] = buffer[v*v_sz + km*im + i + k*im];
        }

        if (i < jm-1 && k < km && i > 0) {
            // left halo
            ((__global float*)&array[i*im + k*im*jm])[v] = buffer[v*v_sz + km*im*2 + i-1 + k*(jm-2)];
            // right halo
            ((__global float*)&array[i*im + k*im*jm + (im-1)])[v] = buffer[v*v_sz + km*im*2 + km*(jm-2) + i-1 + k*(jm-2)];
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
    ) {/*
    const unsigned int v_dim = 2;
    unsigned int i, j, k, v, vec_off, i_buf = 0;
    
    
    for (v = 0; v < v_dim; v++) {
        // Which vector component, ie along v_dim
        vec_off = v * im * jm * km;
        // top halo
        for (k = 0; k < km; k++) {
            for (i = 1; i < im-1; i++) {
                buffer[i_buf] = ((__global float*)&array[i + k*im*jm])[v];
                i_buf++;
            }
        }
        // bottom halo
        for (k = 0; k < km; k++) {
            for (i = 1; i < im-1; i++) {
                buffer[i_buf] = ((__global float*)&array[i + k*im*jm + im*(jm-1)])[v];
                i_buf++;
            }
        }
        // left halo
        for (k = 0; k < km; k++) {
            for (j = 1; j < jm-1; j++) {
                buffer[i_buf] = ((__global float*)&array[j*im + k*im*jm + 1])[v];
                i_buf++;
            }
        }
        // right halo
        for (k = 0; k < km; k++) {
            for (j = 1; j < jm-1; j++) {
                buffer[i_buf] = ((__global float*)&array[j*im + k*im*jm + (im-2)])[v];
                i_buf++;
            }
        }
    }*/
}
