// considering that im,jm,km are always identical to ip,jp,kp which are constants,
// I could simply define the constants using macros and remove these pesky arguments
__kernel void velnw_kernel (
		__global float * p,
		__global float4 * uvw,
		__global float4 * fgh,
		__global float * dxs,
		__global float * dys,
		__global float * dzs,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km,
		const float dt
);
__kernel void bondv1_kernel (
		__global float4 * uvw,
		__global float * dxs,
		__global float * dzn,
		__global float * z2,
#if ICAL == 0
		__global unsigned int n_ptr,
#endif
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km,
		const float dt
);
__kernel void velFG_kernel(
		__global float4* uvw,
		__global float4* fgh,
		__global float16* diu,
		__global float* dzs,
		__global float* dx1,
		__global float* dy1,
		__global float* dzn,
		const int im, const int jm, const int km
);
__kernel void feedbf_kernel (
		__global float4 *uvw,
		__global float4 *uvwsum,
		__global float4 *fgh,
		__global float4 *mask1,
		const float dt,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);
__kernel void les_calc_sm_kernel (
		__global float4 *fgh,
		__global float *dx1,
		__global float *dy1,
		__global float *dzn,
		__global float16 *diu,
		__global float *sm,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);
void merged_velfg_feedbf_les_calc_sm_kernel(

        __global float4* uvw,
        __global float4 *uvwsum,
        __global float4* fgh,
        __global float4 *mask1,
        __global float16* diu,
        __global float* dzs,
        __global float* dx1,
        __global float* dy1,
        __global float* dzn,
        __global float *sm,
        const float dt,
        const unsigned int im,
        const unsigned int jm,
        const unsigned int km

        );
__kernel void les_bound_sm_kernel (
		__global float4 *fgh,
		__global float *dx1,
		__global float *dy1,
		__global float *dzn,
		__global float16 *diu,
		__global float *sm,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);
__kernel void les_calc_visc_kernel (
		__global float4 *fgh,
		__global float *dx1,
		__global float *dy1,
		__global float *dzn,
		__global float16 *diu,
		__global float *sm,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);
__kernel void adam_kernel (
		__global float4 *fgh,
		__global float4 *fgh_old,
		const unsigned int nmax,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);

__kernel void merged_les_calc_visc_adam_kernel (
		__global float4 *fgh,
		__global float4 *fgh_old,
		__global float *dx1,__global float *dy1,__global float *dzn,
		__global float16 *diu,
		__global float *sm,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);

__kernel void press_rhsav_kernel (
		__global float4* uvw,
		__global float4* fgh,
		__global float *rhs,
		__global float *dx1,__global float *dy1,__global float *dzn,
		__global float* chunks_num,
		__global float* chunks_denom,
		const float dt,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);

__kernel void press_sor_kernel (
		__global float4* uvw,
		__global float* p,
#ifdef P_SCRATCH
		__global float* p_scratch,
#endif
		__global float *rhs,
		__global float *cn1,__global float *cn2l,__global float *cn2s,__global float *cn3l,__global float *cn3s,__global float *cn4l,__global float *cn4s,
		__global float *chunks_num,
		__global float *chunks_denom,
		__global unsigned int *val_ptr,
		__global unsigned int *nrd_ptr,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);
__kernel void press_pav_kernel (
		__global float* p,
		__global float *dx1,__global float *dy1,__global float *dzn,
		__global float *chunks_num,
		__global float *chunks_denom,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);
__kernel void press_adj_kernel (
		__global float* p,
		__global float *pav_ptr,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);
__kernel void press_boundp_kernel (
		__global float* p,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);
