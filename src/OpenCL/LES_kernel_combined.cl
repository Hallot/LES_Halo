// This is the combined kernel (super-kernel)
#ifdef __CDT_PARSER__
#include "OpenCLEclipseCompat.hpp"
#endif
// Physical property set
#define RO 1.1763
#define VN 1.583E-5
// IBM parameter set (Feedback force by Goldstein)
#define ALPHA -10.0
#define BETA -1.0
#define TEST

// enum States { ST_INIT,  ST_VELNW, ST_BONDV1, ST_VELFG, ST_FEEDBF, ST_LES, ST_ADAM, ST_PRESS };
#define ST_INIT 0
#ifndef TEST
// considering that im,jm,km are always identical to ip,jp,kp which are constants,
// I could simply define the constants using macros and remove these pesky arguments
void velnw_kernel (
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
void bondv1_kernel (
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
/*
void velFG_kernel(
		__global float4* uvw,
		__global float4* fgh,
		__global float16* diu,
		__global float* dzs,
		__global float* dx1,
		__global float* dy1,
		__global float* dzn,
		const int im, const int jm, const int km
);
void feedbf_kernel (
		__global float4 *uvw,
		__global float4 *uvwsum,
		__global float4 *fgh,
		__global float4 *mask1,
		const float dt,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);
void les_calc_sm_kernel (
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
*/
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
void les_bound_sm_kernel (
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
/*
void les_calc_visc_kernel (
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
void adam_kernel (
		__global float4 *fgh,
		__global float4 *fgh_old,
		const unsigned int nmax,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);
*/
void merged_les_calc_visc_adam_kernel (
		__global float4 *fgh,
		__global float4 *fgh_old,
		__global float *dx1,__global float *dy1,__global float *dzn,
		__global float16 *diu,
		__global float *sm,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);

void press_rhsav_kernel (
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

void press_sor_kernel (
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
void press_pav_kernel (
		__global float* p,
		__global float *dx1,__global float *dy1,__global float *dzn,
		__global float *chunks_num,
		__global float *chunks_denom,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);
void press_adj_kernel (
		__global float* p,
		__global float *pav_ptr,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);
void press_boundp_kernel (
		__global float* p,
		const unsigned int im,
		const unsigned int jm,
		const unsigned int km
);

#if 0
// Function declarations

// velnw takes d{x,y,z}s, u,v,w, p and f,g,h, and modifies u,v,w. This is why we need initial values for all of these
// velnw(km,jm,im,p,ro,dxs,u,dt,f,dys,v,g,dzs,w,h)

void velnw_kernel (
    __global float * p,
    __global float4 * uvw,
    __global float4 * fgh,
    __global float * dxs,
    __global float * dys,
    __global float * dzs,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const float dt,
        )

// bondv1 takes dxs, dzn and z2 and computes u,v,w on the boundaries
//bondv1(jm,u,z2,dzn,v,w,km,ical,n,im,dt,dxs)

void bondv1_kernel (
    __global float4 * uvw,
    __global float * dxs,
    __global float * dzn, 
    __global float * z2,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int n,
    const float dt
        )
//velFG calculates f,g,h from u,v,w and dx1,dy1,dzn,dzs and returns also cov* and diu* from vel2 (the nou* values are unused)
// I will make cov and diu 9-elt arrays for convenience
//  velfg(km,jm,im,dx1,cov1,cov2,cov3,dfu1,diu1,diu2,dy1,diu3,dzn,vn,f,cov4,cov5,cov6,dfv1, &
//      diu4,diu5,diu6,g,cov7,cov8,cov9,dfw1,diu7,diu8,diu9,dzs,h,nou1,u,nou5,v,nou9,w,nou2,nou3, &
//      nou4,nou6,nou7,nou8)

//  velfg(uvw, fgh, cov, diu, nou, dzs, dx1, dy1, dzn, im, jm, km)
//  dfu1,dfv1,dfw1 can be scalar locals
// However, cov, nou and diu must be global, so I think this requires them to be kernel arguments.
void velFG (
    __global float4 * uvw, 
    __global float4 * fgh, 
    __global float16 * cov, 
    __global float16 * diu, 
    __global float16 * nou, 
    __global float * dzs, 
    __global float * dx1, 
    __global float * dy1, 
    __global float * dzn, 
    const unsigned int im, 
    const unsigned int jm, 
    const unsigned int km        
        );
void velFG_kernel(
        __global float * dx1,
        __global float * dy1,
        __global float * dzn,
        __global float * dzs,
        __global float4 * uvw,
        __global float4 * fgh,
        __global float16 * diu,
        const int im,
        const int jm,
        const int km
);

//	vel2 calculates diu{1..9}, cov{1..9} and nou{1..9} from u,v,w and dx1,dy1,dzn,dzs 
//  call vel2(km,jm,im,nou1,u,diu1,dx1,nou5,v,diu5,dy1,nou9,w,diu9,dzn,cov1,cov5,cov9,nou2,diu2, &
//      cov2,nou3,diu3,dzs,cov3,nou4,diu4,cov4,nou6,diu6,cov6,nou7,diu7,cov7,nou8,diu8,cov8)
// vel2(uvw, cov, diu, nou, 

void vel2(
    __global float4 * uvw, 
    __global float16 * diu, 
	float16* cov_ijk,
	float16* cov_ijk_p1,
	float16* diu_ijk,
	float16* diu_ijk_p1,
    __global float * dzs, 
    __global float * dx1, 
    __global float * dy1, 
    __global float * dzn, 
    const unsigned int im, 
    const unsigned int jm, 
    const unsigned int km,
    int i,int j,int k
	);

//feedbf takes u,v,w and {b,c,d}mask1 and modifies f,g,h and {u,v,w}sum so I guess we need initial values for {u,v,w}sum as well
// call feedbf(km,jm,im,usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,fx,fy,fz,f,g,h)
// call feedbf(bmask1,cmask1,dmask1)
// a..d mask1 becomes a 4-elt mask1 array or vector with a the last element, so b/c/d corresponds to x/y/z
// fx,fy,fz can be scalar locals
 // feedbf(uvw, fgh, mask1, im, im, km, dt);
void feedbf (
        __global float4 * uvw, 
        __global float4 * fgh, 
        __global float4 * mask1,
        const unsigned int im, 
        const unsigned int jm, 
        const unsigned int km,             
        const unsigned int dt
        );
void feedbf_kernel (
        __global float *uvw,
        __global float *uvwsum,
        __global float *fgh,
        __global float *fxyz,
        __global float *mask1,
        const float dt,
        const unsigned int im,
        const unsigned int jm,
        const unsigned int km
        )

//les takes diu* (from velFG), dx1,dy1,dzn and modifies f,g,h and calculates delx1 and sm; sm is only used internally for boundsm, similar for delx1, so in the end it only returns f,g,h 
// call les(km,delx1,dx1,dy1,dzn,jm,im,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,sm,f,g,h)
// call les(fgh, diu, sm, dx1,dy1,dzn, im, jm, km)
// delx1 is small, we should be able to make it a local array
// sm needs to be a kernel argument
void les (
   __global float4 * fgh, 
    __global float16 * diu, 
    __global float * sm, 
    __global float * dzs, 
    __global float * dx1, 
    __global float * dy1, 
    __global float * dzn, 
    const unsigned int im, 
    const unsigned int jm, 
    const unsigned int km                
        );

//  call boundsm(km,jm,sm,im)
void boundsm(
           __global float * sm, 
  const unsigned int im, 
    const unsigned int jm, 
    const unsigned int km    
        );


//adam takes f,g,h and {f,g,h}old and fghold and modifies f,g,h and {f,g,h}old + produces output fold, gold, hold, fghold, this means is that fold, gold, hold,
//fghold must be arguments
//call adam(n,nmax,data21,fold,im,jm,km,gold,hold,fghold,f,g,h)
//call adam(fgh, fgh_old,n, im,jm,km, nmax)

void adam (
        __global float4 * fgh, 
        __global float4 * fgh_old, 
        unsigned int n,
        const unsigned int im, 
        const unsigned int jm, 
        const unsigned int km,  
        const unsigned int nmax
        );

//press takes cn{2,3,4}{s,l} and cn1 (which are constant for the run), dx1, dy1, dzn, f,g,h, p, u,v,w, usum,vsum,wsum and modifies p and f,g,h (through bondFG) + produces output u,v,w,p,usum,vsum,wsum
//  call press(km,jm,im,rhs,u,dx1,v,dy1,w,dzn,f,g,h,dt,cn1,cn2l,p,cn2s,cn3l,cn3s,cn4l,cn4s,n, nmax,data20,usum,vsum,wsum)
//  call press(p, uvw, fgh, rhs,dx1,dy1,dzn, cn, n, im, jm, km, dt, nmax)
//  cn2l,cn2s,cn3l,cn3s,cn4l,cn4s are small and constant for the run. cn1 is ip*jp*kp. So let's create cn234ls.
// rhs must be an argument as well
void press(
    __global float * p,         
    __global float4 * uvw,         
    __global float4 * fgh,         
    __global float * rhs, 
    __global float * dx1,
    __global float * dy1,
    __global float * dzn, 
    __global float * cn1,
    __global float * cn234ls,    
    unsigned int n,
    const unsigned int im,
    const unsigned int jm,
    const unsigned int km,
    const unsigned int dt,
    const unsigned int nmax
        );

//   call bondfg(km,jm,f,im,g,h)
void bondFG (
        __global float4 * fgh,  
   const unsigned int im,
    const unsigned int jm,
    const unsigned int km
        
        );
#endif
#endif
// --------------------------------------------------------------------------------
// The actual kernel
// --------------------------------------------------------------------------------

__kernel void LES_combined_kernel (
// The strange order is because of the automatic generations from module_LES_ocl
        __global float* p,
        __global float* cn1,
        __global float* rhs,
        __global float* sm,
        __global float* dxs,
        __global float* dys,
        __global float* dzs,
        __global float* dx1,
        __global float* dy1,
        __global float* dzn,
        __global float* z2,
        __global float* uvw,
        __global float* fgh,
        __global float* fgh_old,
        __global float* cov,
        __global float* diu,
        __global float* nou,
        __global float* mask1,
        __global float* cn234ls,
        __global unsigned int* n_state,
        const float dt, 
        const unsigned int im, 
        const unsigned int jm, 
        const unsigned int km, 
        const unsigned int nmax 
/*
    __global float * p,         
    __global float4 * uvw,         
    __global float4 * fgh,         
    __global float4 * fgh_old,       
    __global float16 * cov, // 3-D array of 9-elt arrays but 16 elts for alignment
    __global float16 * diu, // 3-D array of 9-elt arrays  
    __global float16 * nou, // 3-D array of 9-elt arrays       
    __global float4 * mask1, 
    __global float * rhs, 
    __global float * sm, 
    __global float * dxs,
    __global float * dys,
    __global float * dzs, 
    __global float * dx1,
    __global float * dy1,
    __global float * dzn, 
    __global float * z2,
    __global float * cn,
    __global unsigned int* n_state,
    const float dt,
    const unsigned int ijkm,
//    const unsigned int im,
//    const unsigned int jm,
//    const unsigned int km,
    const unsigned int nmax
*/	
        ) {
#ifndef TEST            
    unsigned int n = n_state[0];
    unsigned int state = n_state[1];

    switch (state) {
        case ST_INIT: 
            {
                // do nothing
                n=0;
                break;
            }
        case ST_VELNW:
            {
                velnw_kernel(p, uvw, fgh, dxs, dys, dzs, im, jm, km,dt);
                break;
            }
        case ST_BONDV1:
            {
                bondv1_kernel(uvw,im,jm,km,dxs,dzn,z2,n,dt);
                break;
            }
        case ST_VELFG:
            {
                velFG_kernel(dx1,dy1,dzn,dzs,uvw,fgh,diu,im,jm,km);
                break;
            }
        case ST_FEEDBF:
            {
//                feedbf(uvw, fgh, mask1, im, im, km, dt);
                feedbf_kernel(uvw, uvwsum, fgh, mask1, ft, im, im, km);
                break;
            }
        case ST_LES_CALC_SM:
            {
                les_calc_sm(fgh, diu, sm, dx1,dy1,dzn, im, jm, km);
                break;
            }
        case ST_LES_BOUND_SM:
            {
                les_bound_sm(fgh, diu, sm, dx1,dy1,dzn, im, jm, km);
                break;
            }
        case ST_LES_CALC_VISC:
            {
                les_calc_visc(fgh, diu, sm, dx1,dy1,dzn, im, jm, km);
                break;
            }

        case ST_ADAM:
            {
                adam(fgh, fgh_old,n, im,jm,km, nmax);
                break;
            }
        case ST_PRESS_RHSAV:
            {
                press_rhsav_kernel(p, uvw, fgh, rhs,dx1,dy1,dzn, cn, n, im, jm, km, dt, nmax)
                break;
            }
        case ST_PRESS_SOR:
            {
                press_rhsav_kernel(p, uvw, fgh, rhs,dx1,dy1,dzn, cn, n, im, jm, km, dt, nmax)
                break;
            }
        case ST_PRESS_PAV:
            {
                press_rhsav_kernel(p, uvw, fgh, rhs,dx1,dy1,dzn, cn, n, im, jm, km, dt, nmax)
                break;
            }
        case ST_PRESS_ADJ:
            {
                press_rhsav_kernel(p, uvw, fgh, rhs,dx1,dy1,dzn, cn, n, im, jm, km, dt, nmax)
                break;
            }
        case ST_PRESS_BOUNDP:
            {
                press_rhsav_kernel(p, uvw, fgh, rhs,dx1,dy1,dzn, cn, n, im, jm, km, dt, nmax)
                break;
            }

        default:    
            n=1;
            // do nothing
    };
#endif            


} // END of LES_kernel

#ifdef TEST
    // Function definitions

#endif    

