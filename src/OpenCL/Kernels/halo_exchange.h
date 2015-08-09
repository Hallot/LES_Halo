void exchange_2_halo_write(
    __global float2 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km
    );

void exchange_2_halo_read(
    __global float2 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km
    );