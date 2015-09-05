// Write the recieved halo buffer to the data array
// Halo order in the buffer: north|south|west|east where north and south contain the corners and thus have two more points along i each
void exchange_1_halo_write(
    __global float *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
) {
    const unsigned int v_dim = 1;
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
    const unsigned int v_sz = 2*km*im*h_w + 2*km*(jm-2*h_w)*h_w;

    // Which vector component, ie along v_dim
    for (unsigned int v = 0; v < v_dim; v++) {
        // index for the halo width
        for (unsigned int w = 0; w < h_w; w++) {
            // index along the row
            // im is one row, times the width of the halo
            if (i < im && k < km){
                // north halo
                // i is the index to copy in the row, k*im*jm is the index of the vertical plane
                // v*v_sz is the vector offset for v_dim, i the index of the row, k*im*h_w the offset for the index of the next plane
                ((__global float*)&array[i + im*w + k*im*jm])[v] = buffer[v*v_sz + i + k*im*h_w];
                // south halo
                // im*(jm-h_w) is the offset to go to the last h_w rows of the plane
                // km*im*h_w is the size of the north halo, so the offset to go to the first element of the south halo
                ((__global float*)&array[i + im*w + k*im*jm + im*(jm-h_w)])[v] = buffer[v*v_sz + km*im*h_w + i + k*im*h_w];
            }

            // index along the column
            // every column except fot the h_w first and last used by the north and south halos
            if (i < jm-h_w && k < km && i > h_w-1) {
                // west halo
                // i*im is the index for the column, w is the index along the row, k*im*jm is the index of the vertical plane
                // km*im*h_w*2 is the size of the north+south halos, so the offset to go to the first element of the east halo
                // (i-h_w)*h_w+w is the index in the halo buffer, since we skip the h_w first rows in the data array and they are not present in the halo buffer, we account for this by removing h_w
                // k*(jm-2*h_w)*h_w is the offset for the index of the vertical plane
                ((__global float*)&array[i*im + w + k*im*jm])[v] = buffer[v*v_sz + km*im*h_w*2 + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w];
                // east halo
                // im-h_w are the last w indexes of the row
                // km*im*h_w*2 + km*(jm-2*h_w)*h_w is the size of the north/south plus west halos, so the offset to go to the first element of the east halo
                ((__global float*)&array[i*im + w + k*im*jm + im-h_w])[v] = buffer[v*v_sz + km*im*h_w*2 + km*(jm-2*h_w)*h_w + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w];
            }
        }
    }
}

void exchange_2_halo_write(
    __global float2 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
) {
    const unsigned int v_dim = 2;
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
    const unsigned int v_sz = 2*km*im*h_w + 2*km*(jm-2*h_w)*h_w;

    // Which vector component, ie along v_dim
    for (unsigned int v = 0; v < v_dim; v++) {
        for (unsigned int w = 0; w < h_w; w++) {
            if (i < im && k < km){
                // north halo
                ((__global float*)&array[i + im*w + k*im*jm])[v] = buffer[v*v_sz + i + k*im*h_w];
                // south halo
                ((__global float*)&array[i + im*w + k*im*jm + im*(jm-h_w)])[v] = buffer[v*v_sz + km*im*h_w + i + k*im*h_w];
            }

            if (i < jm-h_w && k < km && i > h_w-1) {
                // west halo
                ((__global float*)&array[i*im + w + k*im*jm])[v] = buffer[v*v_sz + km*im*h_w*2 + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w];
                // east halo
                ((__global float*)&array[i*im + w + k*im*jm + im-h_w])[v] = buffer[v*v_sz + km*im*h_w*2 + km*(jm-2*h_w)*h_w + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w];
            }
        }
    }
}

void exchange_4_halo_write(
    __global float4 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
) {
    const unsigned int v_dim = 4;
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
    const unsigned int v_sz = 2*km*im*h_w + 2*km*(jm-2*h_w)*h_w;

    // Which vector component, ie along v_dim
    for (unsigned int v = 0; v < v_dim; v++) {
         for (unsigned int w = 0; w < h_w; w++) {
            if (i < im && k < km){
                // north halo
                ((__global float*)&array[i + im*w + k*im*jm])[v] = buffer[v*v_sz + i + k*im*h_w];
                // south halo
                ((__global float*)&array[i + im*w + k*im*jm + im*(jm-h_w)])[v] = buffer[v*v_sz + km*im*h_w + i + k*im*h_w];
            }

            if (i < jm-h_w && k < km && i > h_w-1) {
                // west halo
                ((__global float*)&array[i*im + w + k*im*jm])[v] = buffer[v*v_sz + km*im*h_w*2 + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w];
                // east halo
                ((__global float*)&array[i*im + w + k*im*jm + im-h_w])[v] = buffer[v*v_sz + km*im*h_w*2 + km*(jm-2*h_w)*h_w + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w];
            }
        }
    }
}

void exchange_16_halo_write(
    __global float16 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
) {
    const unsigned int v_dim = 16;
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
    const unsigned int v_sz = 2*km*im*h_w + 2*km*(jm-2*h_w)*h_w;

    // Which vector component, ie along v_dim
    for (unsigned int v = 0; v < v_dim; v++) {
         for (unsigned int w = 0; w < h_w; w++) {
            if (i < im && k < km){
                // north halo
                ((__global float*)&array[i + im*w + k*im*jm])[v] = buffer[v*v_sz + i + k*im*h_w];
                // south halo
                ((__global float*)&array[i + im*w + k*im*jm + im*(jm-h_w)])[v] = buffer[v*v_sz + km*im*h_w + i + k*im*h_w];
            }

            if (i < jm-h_w && k < km && i > h_w-1) {
                // west halo
                ((__global float*)&array[i*im + w + k*im*jm])[v] = buffer[v*v_sz + km*im*h_w*2 + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w];
                // east halo
                ((__global float*)&array[i*im + w + k*im*jm + im-h_w])[v] = buffer[v*v_sz + km*im*h_w*2 + km*(jm-2*h_w)*h_w + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w];
            }
        }
    }
}

// Write the data from the data array into the halo buffer
// Halo order in the buffer: north|south|west|east where north and south do not contain the corners and thus have the same size as west and east
void exchange_1_halo_read(
    __global float *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
) {
    const unsigned int v_dim = 1;
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
    const unsigned int v_sz = 2 * km * (im+jm-4*h_w)*h_w;

    // Which vector component, ie along v_dim
    for (unsigned int v = 0; v < v_dim; v++) {
        for (unsigned int w = 0; w < h_w; w++) {
            if (i < im-h_w && k < km && i > h_w-1){
                // north halo
                buffer[v*v_sz + i-h_w + (im-2*h_w)*w + k*(im-2*h_w)*h_w] = ((__global float*)&array[im*h_w + i + im*w + k*im*jm])[v];
                // south halo
                buffer[v*v_sz + km*(im-2*h_w)*h_w + i-h_w + (im-2*h_w)*w + k*(im-2*h_w)*h_w] = ((__global float*)&array[i + k*im*jm + im*(jm-2*h_w) + im*w])[v];
            }
            
            if (i < jm-h_w && k < km && i > h_w-1) {
                // west halo
                buffer[v*v_sz + km*(im-2*h_w)*h_w*2 + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w] = ((__global float*)&array[i*im + w + h_w + k*im*jm])[v];
                // east halo
                buffer[v*v_sz + km*(im-2*h_w)*h_w*2 + km*(jm-2*h_w)*h_w + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w] = ((__global float*)&array[i*im + w + k*im*jm + im-2*h_w])[v];
            }
        }
    }
}

void exchange_2_halo_read(
    __global float2 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
    ) {
    const unsigned int v_dim = 2;
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
    const unsigned int v_sz = 2*km*im*h_w + 2*km*(jm-2*h_w)*h_w - 4*h_w*h_w*km;

    // Which vector component, ie along v_dim
    for (unsigned int v = 0; v < v_dim; v++) {
        for (unsigned int w = 0; w < h_w; w++) {
            if (i < im-h_w && k < km && i > h_w-1){
                // north halo
                buffer[v*v_sz + i-h_w + (im-2*h_w)*w + k*(im-2*h_w)*h_w] = ((__global float*)&array[im*h_w + i + im*w + k*im*jm])[v];
                // south halo
                buffer[v*v_sz + km*(im-2*h_w)*h_w + i-h_w + (im-2*h_w)*w + k*(im-2*h_w)*h_w] = ((__global float*)&array[i + k*im*jm + im*(jm-2*h_w) + im*w])[v];
            }
            
            if (i < jm-h_w && k < km && i > h_w-1) {
                // west halo
                buffer[v*v_sz + km*(im-2*h_w)*h_w*2 + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w] = ((__global float*)&array[i*im + w + h_w + k*im*jm])[v];
                // east halo
                buffer[v*v_sz + km*(im-2*h_w)*h_w*2 + km*(jm-2*h_w)*h_w + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w] = ((__global float*)&array[i*im + w + k*im*jm + im-2*h_w])[v];
            }
        }
    }
}

void exchange_4_halo_read(
    __global float4 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
    ) {
    const unsigned int v_dim = 4;
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
    const unsigned int v_sz = 2*km*im*h_w + 2*km*(jm-2*h_w)*h_w - 4*h_w*h_w*km;

    // Which vector component, ie along v_dim
    for (unsigned int v = 0; v < v_dim; v++) {
        for (unsigned int w = 0; w < h_w; w++) {
            if (i < im-h_w && k < km && i > h_w-1){
                // north halo
                buffer[v*v_sz + i-h_w + (im-2*h_w)*w + k*(im-2*h_w)*h_w] = ((__global float*)&array[im*h_w + i + im*w + k*im*jm])[v];
                // south halo
                buffer[v*v_sz + km*(im-2*h_w)*h_w + i-h_w + (im-2*h_w)*w + k*(im-2*h_w)*h_w] = ((__global float*)&array[i + k*im*jm + im*(jm-2*h_w) + im*w])[v];
            }
            
            if (i < jm-h_w && k < km && i > h_w-1) {
                // west halo
                buffer[v*v_sz + km*(im-2*h_w)*h_w*2 + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w] = ((__global float*)&array[i*im + w + h_w + k*im*jm])[v];
                // east halo
                buffer[v*v_sz + km*(im-2*h_w)*h_w*2 + km*(jm-2*h_w)*h_w + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w] = ((__global float*)&array[i*im + w + k*im*jm + im-2*h_w])[v];
            }
        }
    }
}


void exchange_16_halo_read(
    __global float16 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
    ) {
    const unsigned int v_dim = 16;
    const unsigned int gl_id = get_global_id(0);
    const unsigned int k = gl_id % km;
    const unsigned int i = gl_id / km;
    const unsigned int v_sz = 2*km*im*h_w + 2*km*(jm-2*h_w)*h_w - 4*h_w*h_w*km;

    // Which vector component, ie along v_dim
    for (unsigned int v = 0; v < v_dim; v++) {
        for (unsigned int w = 0; w < h_w; w++) {
            if (i < im-h_w && k < km && i > h_w-1){
                // north halo
                buffer[v*v_sz + i-h_w + (im-2*h_w)*w + k*(im-2*h_w)*h_w] = ((__global float*)&array[im*h_w + i + im*w + k*im*jm])[v];
                // south halo
                buffer[v*v_sz + km*(im-2*h_w)*h_w + i-h_w + (im-2*h_w)*w + k*(im-2*h_w)*h_w] = ((__global float*)&array[i + k*im*jm + im*(jm-2*h_w) + im*w])[v];
            }
            
            if (i < jm-h_w && k < km && i > h_w-1) {
                // west halo
                buffer[v*v_sz + km*(im-2*h_w)*h_w*2 + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w] = ((__global float*)&array[i*im + w + h_w + k*im*jm])[v];
                // east halo
                buffer[v*v_sz + km*(im-2*h_w)*h_w*2 + km*(jm-2*h_w)*h_w + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w] = ((__global float*)&array[i*im + w + k*im*jm + im-2*h_w])[v];
            }
        }   
    }
}