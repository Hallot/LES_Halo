/*
 The arguments that change between kernel calls are:
 n_state (always)
 nrd_ptr (only for press_sor)
 chunks_num
 chunks_denom
 */
combined_kernel (
		__global float* p,
		__global float* p_scratch, // only used if P_SCRATCH is defined
		__global float4* uvw,
		__global float4* uvwsum,
		__global float4* fgh,
		__global float4* fgh_old,
		__global float* rhs,
		__global float4* mask1,
		__global float16* diu,
		__global float* sm,
		__global float* dxs,
		__global float* dys,
		__global float* dzs,
		__global float* dx1,
		__global float* dy1,
		__global float* dzn,
		__global float* z2,
		__global float* cn1,
		__global float* cn2l,
		__global float* cn2s,
		__global float* cn3l,
		__global float* cn3s,
		__global float* cn4l,
		__global float* cn4s,
		__global float* val_ptr, // we use this for any transfer of scalar values
		__global float* chunks_num,
		__global float* chunks_denom,
		__global unsigned int* n_ptr, // also used for nrd
		__global unsigned int* n_state,
		const float dt,
		const unsigned int nmax, // this is actually unused
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);
