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
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
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

void exchange_4_halo_write(
    __global float4 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km
) {
    const unsigned int v_dim = 4;
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
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

void exchange_16_halo_write(
    __global float16 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km
) {
    const unsigned int v_dim = 16;
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
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
    ) {const unsigned int v_dim = 2;
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
    const unsigned int v_sz = 2*km*im + 2*km*(jm-2);

    // Which vector component, ie along v_dim
    for (unsigned int v = 0; v < v_dim; v++) {
        if (i < im && k < km){
            // top halo
            buffer[v*v_sz + i + k*im] = ((__global float*)&array[i + k*im*jm])[v];
            // bottom halo
            buffer[v*v_sz + km*im + i + k*im] = ((__global float*)&array[i + k*im*jm + im*(jm-1)])[v];
        }

        if (i < jm-1 && k < km && i > 0) {
            // left halo
            buffer[v*v_sz + km*im*2 + i-1 + k*(jm-2)] = ((__global float*)&array[i*im + k*im*jm])[v];
            // right halo
            buffer[v*v_sz + km*im*2 + km*(jm-2) + i-1 + k*(jm-2)] = ((__global float*)&array[i*im + k*im*jm + (im-1)])[v];
        }
    }
}

void exchange_4_halo_read(
    __global float4 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km
    ) {const unsigned int v_dim = 4;
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
    const unsigned int v_sz = 2*km*im + 2*km*(jm-2);

    // Which vector component, ie along v_dim
    for (unsigned int v = 0; v < v_dim; v++) {
        if (i < im && k < km){
            // top halo
            buffer[v*v_sz + i + k*im] = ((__global float*)&array[i + k*im*jm])[v];
            // bottom halo
            buffer[v*v_sz + km*im + i + k*im] = ((__global float*)&array[i + k*im*jm + im*(jm-1)])[v];
        }

        if (i < jm-1 && k < km && i > 0) {
            // left halo
            buffer[v*v_sz + km*im*2 + i-1 + k*(jm-2)] = ((__global float*)&array[i*im + k*im*jm])[v];
            // right halo
            buffer[v*v_sz + km*im*2 + km*(jm-2) + i-1 + k*(jm-2)] = ((__global float*)&array[i*im + k*im*jm + (im-1)])[v];
        }
    }
}


void exchange_16_halo_read(
    __global float16 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km
    ) {const unsigned int v_dim = 16;
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
    const unsigned int v_sz = 2*km*im + 2*km*(jm-2);

    // Which vector component, ie along v_dim
    for (unsigned int v = 0; v < v_dim; v++) {
        if (i < im && k < km){
            // top halo
            buffer[v*v_sz + i + k*im] = ((__global float*)&array[i + k*im*jm])[v];
            // bottom halo
            buffer[v*v_sz + km*im + i + k*im] = ((__global float*)&array[i + k*im*jm + im*(jm-1)])[v];
        }

        if (i < jm-1 && k < km && i > 0) {
            // left halo
            buffer[v*v_sz + km*im*2 + i-1 + k*(jm-2)] = ((__global float*)&array[i*im + k*im*jm])[v];
            // right halo
            buffer[v*v_sz + km*im*2 + km*(jm-2) + i-1 + k*(jm-2)] = ((__global float*)&array[i*im + k*im*jm + (im-1)])[v];
        }
    }
}