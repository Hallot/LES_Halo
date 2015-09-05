#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX(a,b) (((a)>(b))?(a):(b))

void exchange_2_halo_write(
    float *array,
    float *buffer,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int h_w
) {
    const unsigned int v_dim = 2;
    const unsigned int v_sz = 2*km*im*h_w + 2*km*(jm-2*h_w)*h_w;
    const unsigned int arr_sz = im * jm * km;
    unsigned int k;
    unsigned int i;

    for (unsigned int gl_id = 0; gl_id < km * MAX(im,jm); gl_id++) {
        k = gl_id % km;
        i = gl_id / km;
        // Which vector component, ie along v_dim
        for (unsigned int v = 0; v < v_dim; v++) {
            // index along the row
            for (unsigned int w = 0; w < h_w; w++) {
                if (i < im && k < km){
                    // north halo
                    array[i + im*w + k*im*jm + v*arr_sz] = buffer[v*v_sz + i + im*w + k*im*h_w];
                    // south halo
                    array[i + im*w + k*im*jm + im*(jm-h_w) + v*arr_sz] = buffer[v*v_sz + km*im*h_w + i + im*w + k*im*h_w];
                }

                if (i < jm-h_w && k < km && i > h_w-1) {
                    // west halo
                    array[i*im + w + k*im*jm + v*arr_sz] = buffer[v*v_sz + km*im*h_w*2 + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w];
                    // east halo
                    array[i*im + w + k*im*jm + im-h_w + v*arr_sz] = buffer[v*v_sz + km*im*h_w*2 + km*(jm-2*h_w)*h_w + (i-h_w)*h_w+w + k*(jm-2*h_w)*h_w];
                }
            }
        }
    }
}


int main(void) {
    const rlim_t kStackSize = 64L * 1024L * 1024L;   // min stack size = 64 Mb
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
        }
    }
    
    const unsigned int h_w = 3;
    const unsigned int im = 20;
    const unsigned int jm = 20;
    const unsigned int km = 2;
    const unsigned int v_dim = 2;
    const unsigned int arr_sz = v_dim * im * jm * km;
    const unsigned int buf_sz = v_dim * 2 * (im+jm-2*h_w)*h_w * km;
    float array[arr_sz];
    float buffer[buf_sz];
    int i;
    
    
    for (i = 0; i < arr_sz / 2; i++) {
        array[i] = 1.0;
    }
    for (i = arr_sz / 2; i < arr_sz; i++) {
        array[i] = 2.0;
    }
    for (i = 0; i < buf_sz; i++) {
        buffer[i] = 5.0;
    }
    
    exchange_2_halo_write(array, buffer, im, jm, km, h_w);
    
    for (i = 0; i < arr_sz; i++) {
        if (i % (arr_sz/v_dim) == 0) {
            printf("\n\n==vect==");
        }
        if (i % (im*jm) == 0) {
            printf("\n==k==");
        }
        if (i % im == 0) {
            printf("\n");
        }
        printf("%g ", array[i]);
    }
    printf("\n");
    
    return 0;
}