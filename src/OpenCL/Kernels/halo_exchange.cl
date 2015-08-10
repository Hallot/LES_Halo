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
    unsigned int i, j, k, v, i_buf = 0;
    
    // Which vector component, ie along v_dim
    for (v = 0; v < v_dim; v++) {
        // top halo
        for (k = 0; k < km; k++) {
            for (i = 0; i < im; i++) {
                ((__global float*)&array[i + k*im*jm])[v] = buffer[i_buf];
                i_buf++;
            }
        }
        // bottom halo
        for (k = 0; k < km; k++) {
            for (i = 0; i < im; i++) {
                ((__global float*)&array[i + k*im*jm + im*(jm-1)])[v] = buffer[i_buf];
                i_buf++;
            }
        }
        // left halo
        for (k = 0; k < km; k++) {
            for (j = 1; j < jm-1; j++) {
                ((__global float*)&array[j*im + k*im*jm])[v] = buffer[i_buf];
                i_buf++;
            }
        }
        // right halo
        for (k = 0; k < km; k++) {
            for (j = 1; j < jm-1; j++) {
                ((__global float*)&array[j*im + k*im*jm + (im-1)])[v] = buffer[i_buf];
                i_buf++;
            }
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
