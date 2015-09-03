void exchange_1_halo_write(
    __global float *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
    );

void exchange_2_halo_write(
    __global float2 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
    );

void exchange_4_halo_write(
    __global float4 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
    );

void exchange_16_halo_write(
    __global float16 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
    );

void exchange_1_halo_read(
    __global float *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
    );

void exchange_2_halo_read(
    __global float2 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
    );

void exchange_4_halo_read(
    __global float4 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
    );

void exchange_16_halo_read(
    __global float16 *array,
    __global float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
    );