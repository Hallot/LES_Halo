! This is a cleaned-up version of the module template, without debug statements.
! For the full version see Revision: 7743
#define GPU_KERNEL 1
#define CPU_KERNEL 2

!#define FIXED_NTH
!#define APPROX_PAR_SOR
!#define TWINNED_DOUBLE_BUFFER

module module_LES_combined_ocl
    use module_LES_conversions
!    use module_LES_tests
#ifdef USE_NETCDF_OUTPUT
    use module_LES_write_netcdf
#endif

contains
    subroutine initialise_LES_kernel ( &
        p,u,v,w,usum,vsum,wsum,f,g,h,fold,gold,hold, &
        diu1, diu2, diu3, diu4, diu5, diu6, diu7, diu8, diu9, &
        amask1, bmask1, cmask1, dmask1, &
        cn1, cn2l, cn2s, cn3l, cn3s, cn4l, cn4s, &
        rhs, sm, dxs, dys, dzs, dx1, dy1, dzn, &
#ifndef EXTERNAL_WIND_PROFILE
        z2, &
#else
        wind_profile, &
#endif
        dt, im, jm, km &
        !, tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, &
        !t_out, r_out, b_out, l_out &
        )
        use params_common_sn

        implicit none

        ! Arguments
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1), intent(InOut)  :: p

        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(InOut) :: u
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: usum
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(InOut) :: v
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: vsum
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(InOut) :: w
        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(InOut) :: wsum

        real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut)  :: f
        real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut)  :: g
        real(kind=4), dimension(0:ip,0:jp,0:kp), intent(InOut)  :: h
        real(kind=4), dimension(ip,jp,kp), intent(InOut)  :: fold
        real(kind=4), dimension(ip,jp,kp), intent(InOut)  :: gold
        real(kind=4), dimension(ip,jp,kp), intent(InOut)  :: hold

        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2), intent(In)  :: diu9

        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1), intent(In)  :: amask1
        real(kind=4), dimension(-1:ip+1,0:jp+1,0:kp+1), intent(In)  :: bmask1
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(In)  :: cmask1
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1), intent(In)  :: dmask1

        real(kind=4), dimension(ip,jp,kp), intent(In)  :: cn1

        real(kind=4), dimension(ip), intent(In)  :: cn2l
        real(kind=4), dimension(ip), intent(In)  :: cn2s
        real(kind=4), dimension(jp), intent(In)  :: cn3l
        real(kind=4), dimension(jp), intent(In)  :: cn3s
        real(kind=4), dimension(kp), intent(In)  :: cn4l
        real(kind=4), dimension(kp), intent(In)  :: cn4s

        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1), intent(In)  :: rhs
        real(kind=4), dimension(-1:ip+1,-1:jp+1,0:kp+1), intent(In)  :: sm
        real(kind=4), dimension(0:ip), intent(In)  :: dxs
        real(kind=4), dimension(0:jp), intent(In)  :: dys
        real(kind=4), dimension(-1:kp+2), intent(In)  :: dzs
        real(kind=4), dimension(-1:ip+1), intent(In)  :: dx1
        real(kind=4), dimension(0:jp+1), intent(In)  :: dy1
        real(kind=4), dimension(-1:kp+2), intent(In)  :: dzn
#ifndef EXTERNAL_WIND_PROFILE
        real(kind=4), dimension(kp+2)  :: z2        
#endif
        real(kind=4), intent(In) :: dt
        integer, intent(In) :: im
        integer, intent(In) :: jm
        integer, intent(In) :: km
!        integer, intent(In) :: nmax

        integer, parameter :: v_dim = 4
        integer, parameter :: h_w = 1
        integer, parameter :: h_h = 1
        integer, parameter :: h_d = kp + 3
        integer, parameter :: c_w = ip + 3
        integer, parameter :: c_h = jp + 3
        
        real(kind=4), dimension(v_dim, c_w, c_h, h_d)  :: tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, r_out, l_out
        real(kind=4), dimension(v_dim, c_w + 2 * h_w, c_h, h_d) :: t_out, b_out

        ! -----------------------------------------------------------------------
        ! Combined arrays for OpenCL kernels

        real(kind=4), dimension(0:3,0:ip+1,-1:jp+1,-1:kp+1)  :: uvw
        real(kind=4), dimension(0:3,0:ip,0:jp,0:kp) :: uvwsum
        real(kind=4), dimension(0:3,0:ip,0:jp,0:kp)  :: fgh
        real(kind=4), dimension(0:3,ip,jp,kp)  :: fgh_old
!        real(kind=4), dimension(1:16,-1:ip+2,0:jp+2,0:kp+2)  :: cov
        real(kind=4), dimension(1:16,-1:ip+2,0:jp+2,0:kp+2)  :: diu
!        real(kind=4), dimension(1:16,-1:ip+2,0:jp+2,0:kp+2)  :: nou
        real(kind=4), dimension(0:3,-1:ip+1,-1:jp+1,0:kp+1)  :: mask1

        ! Extra arguments for kernel
        real(kind=4), dimension(0:1,0:ip+2,0:jp+2,0:kp+1) :: p_scratch
        real(kind=4), dimension(512) :: chunks_num, chunks_denom
        real(kind=4), dimension(256) :: val_ptr
        integer, dimension(256) :: n_ptr
        integer, dimension(256) :: state_ptr
#ifdef EXTERNAL_WIND_PROFILE
        real(kind=4), dimension(jp,kp), intent(In) :: wind_profile
#endif



        ! Convert to new format
        call convert_to_uvw(u,v,w,uvw)
        call convert_to_fgh(usum,vsum,wsum,uvwsum)
        call convert_to_fgh(f,g,h,fgh)
        call convert_to_fgh_old(fold,gold,hold, fgh_old)
        ! The following are all read-only
        call convert_to_9vec(diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,diu)
        call convert_to_mask1(amask1,bmask1,cmask1,dmask1,mask1)
        state_ptr(1) = 0

        p_scratch(0,:,:,:) = p(:,:,:)
        p_scratch(1,:,:,:) = p(:,:,:)

        chunks_num=0.0
        chunks_denom=0.0
        tl_in=1.0

        !$ACC Kernel(LES_combined_kernel_mono), GlobalRange(1), LocalRange(1)
        call LES_combined(p_scratch, uvw, uvwsum, fgh, fgh_old, &
            rhs, mask1, diu, sm,&
            dxs, dys, dzs, dx1, dy1, dzn, &
#ifndef EXTERNAL_WIND_PROFILE
            z2, &
#else
            wind_profile, &
#endif
            cn1, cn2l, cn2s, cn3l, cn3s, cn4l, cn4s,&
            val_ptr, chunks_num, chunks_denom, n_ptr, state_ptr, dt, im, jm, km &
            , tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, r_out, b_out, l_out &
            )
        !$ACC End Kernel

        ! Following buffers are used in the loop, assign to module-level buffer array for convenience
        oclBuffers(1) = p_scratch_buf ! BOUNDP
        oclBuffers(2) = uvw_buf ! BOUNDP
        oclBuffers(3) = uvwsum_buf ! BOUNDP
        oclBuffers(4) = fgh_buf ! DEBUG
        oclBuffers(5) = fgh_old_buf ! ADAM
        oclBuffers(6) = val_ptr_buf ! SOR, ADJ
        oclBuffers(7) = chunks_num_buf ! RHSAV, SOR, PAV
        oclBuffers(8) = chunks_denom_buf ! RHSAV, PAV
        oclBuffers(9) = n_ptr_buf ! BONDV1, SOR
        oclBuffers(10) = state_ptr_buf ! ALL
#ifdef EXTERNAL_WIND_PROFILE
        oclBuffers(11) = wind_profile_buf ! BONDV1
#endif
        oclBuffers(12) = tl_in_buf ! Halo
        oclBuffers(13) = t_in_buf ! Halo
        oclBuffers(14) = tr_in_buf ! Halo
        oclBuffers(15) = r_in_buf ! Halo
        oclBuffers(16) = br_in_buf ! Halo
        oclBuffers(17) = b_in_buf ! Halo
        oclBuffers(18) = bl_in_buf ! Halo
        oclBuffers(19) = l_in_buf ! Halo
        oclBuffers(20) = t_out_buf ! Halo
        oclBuffers(21) = r_out_buf ! Halo
        oclBuffers(22) = b_out_buf ! Halo
        oclBuffers(23) = l_out_buf ! Halo
        oclNunits = initialise_LES_kernel_nunits
        oclNthreadsHint = initialise_LES_kernel_nthreads

#ifdef VERBOSE
        print *,'CL_DEVICE_MAX_COMPUTE_UNITS:',oclNunits
        print *,'CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE:',oclNthreadsHint
#endif

    end subroutine initialise_LES_kernel 
    ! --------------------------------------------------------------------------------
    ! --------------------------------------------------------------------------------

    subroutine run_LES_kernel ( &
#ifndef NO_FILE_IO
        data20, data21, &
#endif
!        im, jm, km, &
#ifdef VERBOSE
        dt, dx1,dy1,dzn, &
#endif
        n, nmax &
#ifdef EXTERNAL_WIND_PROFILE
        , wind_profile &
#endif
        )

        use oclWrapper
        use params_common_sn

        implicit none
#ifndef NO_FILE_IO
        character(len=70), intent(In) :: data20, data21
#endif
!        integer, intent(In) :: im
!        integer, intent(In) :: jm
!        integer, intent(In) :: km
#ifdef VERBOSE
        real(kind=4), intent(In) :: dt
        real(kind=4), dimension(-1:ip+1), intent(In)  :: dx1
        real(kind=4), dimension(0:jp+1), intent(In)  :: dy1
        real(kind=4), dimension(-1:kp+2), intent(In)  :: dzn
#endif
        integer, intent(In) :: n, nmax
        
        ! -----------------------------------------------------------------------
        ! arrays for OpenCL kernels
!        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1), intent(In)  :: p
!        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: u
!        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
!        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1) , intent(In) :: w
!        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In):: f
!        real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: g
!        real(kind=4), dimension(0:ip,0:jp,0:kp), intent(In)  :: h
!        real(kind=4), dimension(ip,jp,kp) , intent(In) :: fold
!        real(kind=4), dimension(ip,jp,kp), intent(In) :: gold
!        real(kind=4), dimension(ip,jp,kp) , intent(In) :: hold

        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: usumo
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: vsumo
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: wsumo

        real(kind=4), dimension(0:3,0:ip+1,-1:jp+1,-1:kp+1)  :: uvw
        real(kind=4), dimension(0:3,0:ip,0:jp,0:kp) :: uvwsum
        real(kind=4), dimension(0:3,0:ip,0:jp,0:kp)  :: fgh
        real(kind=4), dimension(0:3,ip,jp,kp)  :: fgh_old

        real(kind=4), dimension(0:1,0:ip+2,0:jp+2,0:kp+1)  :: po
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: uo
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: vo
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1)  :: wo
!        real(kind=4), dimension(0:ip,0:jp,0:kp) :: fo
!        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: go
!        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: ho
        real(kind=4), dimension(ip,jp,kp)  :: foldo
        real(kind=4), dimension(ip,jp,kp) :: goldo
        real(kind=4), dimension(ip,jp,kp)  :: holdo

        real(kind=4), dimension(256) :: val_ptr
        real(kind=4), dimension(512) :: chunks_num, chunks_denom
        integer, dimension(256) :: n_ptr
        integer, dimension(256) :: state_ptr
        
        integer, parameter :: v_dim = 4
        integer, parameter :: h_w = 1
        integer, parameter :: h_h = 1
        integer, parameter :: h_d = kp + 3
        integer, parameter :: c_w = ip + 3
        integer, parameter :: c_h = jp + 3
        real(kind=4), dimension(v_dim, c_w, c_h, h_d)  :: tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, r_out, l_out
        real(kind=4), dimension(v_dim, c_w + 2 * h_w, c_h, h_d) :: t_out, b_out
        !===== Test =====!
        !real(kind=4), dimension(v_dim, c_w, c_h, h_d) :: test_halo
        !real(kind=4), dimension(v_dim, c_w + 2 * h_w, c_h, h_d) :: big_test_halo
        
        integer(8) :: p_buf
        integer(8) :: uvw_buf
        integer(8) :: uvwsum_buf
        integer(8) :: fgh_old_buf, fgh_buf
        integer(8) :: chunks_num_buf, chunks_denom_buf
        integer(8) :: val_ptr_buf, n_ptr_buf, state_ptr_buf
        integer(8) :: tl_in_buf, t_in_buf, tr_in_buf, r_in_buf, br_in_buf, b_in_buf, bl_in_buf, l_in_buf, t_out_buf, r_out_buf, b_out_buf, l_out_buf
#ifdef EXTERNAL_WIND_PROFILE
        real(kind=4), dimension(jp,kp) , intent(In) :: wind_profile
        integer(8) :: wind_profile_bufshape(l_out)
        integer, dimension(2) :: wind_profile_sz
#endif

        integer, dimension(4) :: p_sz, uvw_sz, uvwsum_sz, fgh_old_sz,fgh_sz
        integer, dimension(1) :: state_ptr_sz, n_ptr_sz, val_ptr_sz, chunks_num_sz, chunks_denom_sz
        integer, dimension(4) :: tl_in_sz, t_in_sz, tr_in_sz, r_in_sz, br_in_sz, b_in_sz, bl_in_sz, l_in_sz, t_out_sz, r_out_sz, b_out_sz, l_out_sz

        integer :: state, nn, ii
        real(kind=4) :: rhsav, area, pav, pco, sor
#if KERNEL == CPU_KERNEL
        real(kind=4) :: aaa, bbb
#endif
        
#ifdef VERBOSE
        real(kind=4) :: cflu, cflv, cflw
        integer :: i, j, k
#endif
#ifndef APPROX_PAR_SOR
        integer, parameter  :: nmaxp = 50 ! 5 ! 50 ! 100
#else
        integer, parameter  :: nmaxp = 1
#endif
        ! --only used mac method
        real, parameter  :: pjuge = 0.0001

        integer :: ngroups, nrd, iter
        real(kind=4) :: avg_iter
#if KERNEL == GPU_KERNEL
        integer :: max_range
#endif

        integer, parameter :: ST_INIT=0, ST_VELNW__BONDV1_INIT_UVW=1, ST_BONDV1_CALC_UOUT=2
        integer, parameter :: ST_BONDV1_CALC_UVW=3, ST_VELFG__FEEDBF__LES_CALC_SM=4, ST_LES_BOUND_SM=5, ST_LES_CALC_VISC__ADAM=6
        integer, parameter :: ST_BONDV1_CALC_UVW__VELFG__FEEDBF__LES_CALC_SM=30, ST_VELFG=31, ST_FEEDBF__LES_CALC_SM=32

        integer, parameter :: ST_PRESS_RHSAV=7, ST_PRESS_SOR=8, ST_PRESS_PAV=9, ST_PRESS_ADJ=10, ST_PRESS_BOUNDP=11, ST_DONE=12
        
        real (kind=4) :: exectime
#ifdef TIMINGS
        real (kind=4) :: sor_exectime, time_start, time_stop !
        integer :: itim
        real (kind=4), dimension(0:12) :: ktimestamp
#endif

        integer, dimension(4), parameter :: lb = (/1,1,1,1/), ub = (/v_dim, ip + 3, jp + 3, kp + 3/)
        



        foldo=0
        goldo=0
        holdo=0
        p_buf = oclBuffers(1) 
        uvw_buf = oclBuffers(2)
        uvwsum_buf = oclBuffers(3)
        fgh_buf = oclBuffers(4)
        fgh_old_buf = oclBuffers(5)
        val_ptr_buf=oclBuffers(6) ! SOR, ADJ
        chunks_num_buf = oclBuffers(7) ! RHSAV, SOR, PAV
        chunks_denom_buf = oclBuffers(8) ! RHSAV, PAV
        n_ptr_buf = oclBuffers(9) ! BONDV1, SOR
        state_ptr_buf = oclBuffers(10) ! ALL
#ifdef EXTERNAL_WIND_PROFILE
        wind_profile_buf = oclBuffers(11)  ! BONDV1
#endif
        tl_in_buf = oclBuffers(12) ! Halo
        t_in_buf = oclBuffers(13) ! Halo
        tr_in_buf = oclBuffers(14) ! Halo
        r_in_buf = oclBuffers(15) ! Halo
        br_in_buf = oclBuffers(16) ! Halo
        b_in_buf = oclBuffers(17) ! Halo
        bl_in_buf = oclBuffers(18) ! Halo
        l_in_buf = oclBuffers(19) ! Halo
        t_out_buf = oclBuffers(20) ! Halo
        r_out_buf = oclBuffers(21) ! Halo
        b_out_buf = oclBuffers(22) ! Halo
        l_out_buf = oclBuffers(23) ! Halo

        p_sz = shape(po)
        uvw_sz = shape(uvw)
        uvwsum_sz = shape(uvwsum)
        fgh_sz = shape(fgh)
        fgh_old_sz = shape(fgh_old)
        val_ptr_sz = shape(val_ptr)
        chunks_num_sz = shape(chunks_num)
        chunks_denom_sz = shape(chunks_denom)
        n_ptr_sz = shape(n_ptr)
        state_ptr_sz = shape(state_ptr)
#ifdef EXTERNAL_WIND_PROFILE
        wind_profile_sz = shape(wind_profile)  ! BONDV1
#endif
        tl_in_sz = shape(tl_in)
        t_in_sz = shape(t_in)
        tr_in_sz = shape(tr_in)
        r_in_sz = shape(r_in)
        br_in_sz = shape(br_in)
        b_in_sz = shape(b_in)
        bl_in_sz = shape(bl_in)
        l_in_sz = shape(l_in)
        t_out_sz = shape(t_out)
        r_out_sz = shape(r_out)
        b_out_sz = shape(b_out)
        l_out_sz = shape(l_out)

        n_ptr(1)=n
#ifdef TIMINGS
         print *, 'run_LES_kernel: time step = ',n
#endif

        !tl_in = 1.0
        !t_in = 2.0
        !tr_in = 3.0
        !r_in = 4.0
        !br_in = 5.0
        !b_in = 6.0
        !bl_in = 7.0
        !l_in = 8.0
        
        
        !call oclWrite4DFloatArrayBuffer(t_in_buf, t_in_sz, t_in)
        !call oclWrite4DFloatArrayBuffer(r_in_buf, r_in_sz, r_in)
        !call oclWrite4DFloatArrayBuffer(b_in_buf, b_in_sz, b_in)
        !call oclWrite4DFloatArrayBuffer(l_in_buf, l_in_sz, l_in)
        !call oclWrite4DFloatArrayBuffer(tl_in_buf, tl_in_sz, tl_in)
        !call oclWrite4DFloatArrayBuffer(tr_in_buf, tr_in_sz, tr_in)
        !call oclWrite4DFloatArrayBuffer(br_in_buf, br_in_sz, br_in)
        !call oclWrite4DFloatArrayBuffer(bl_in_buf, bl_in_sz, bl_in)
        
        !call oclRead4DFloatArrayBuffer(t_in_buf, t_in_sz, t_in)
        !t_in = 9.0
        !call oclWrite4DFloatArrayBuffer(t_in_buf, t_in_sz, t_in)
        !t_in = 2.0
        !call oclRead4DFloatArrayBuffer(t_in_buf, t_in_sz, t_in)
         
        
        !call read_ext_halos_out(tl_in, tr_in, br_in, bl_in, &
         !   v_dim, h_h, h_d, c_w, &
          !  tl_in_buf, tr_in_buf, br_in_buf, bl_in_buf, tl_in_sz, tr_in_sz, br_in_sz, bl_in_sz)
            
        !call read_halos_out(t_in, r_in, b_in, l_in, &
         !   v_dim, h_h, h_d, c_w, &
          !  t_in_buf, r_in_buf, b_in_buf, l_in_buf, t_in_sz, r_in_sz, b_in_sz, l_in_sz)
        
        !call merge_halos_in(tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, b_out, &
         !   v_dim, h_w, h_h, h_d, c_w, &
          !  t_out_buf, r_out_buf, b_out_buf, l_out_buf, t_out_sz, r_out_sz, b_out_sz, l_out_sz)
            
        !b_out = -2.0
            
        !call oclRead4DFloatArrayBuffer(t_out_buf, t_out_sz, b_out)
        
        !call compare_halos(t_out, b_out, lb, ub)


        ! ========================================================================================================================================================
        ! ========================================================================================================================================================
#ifdef TIMINGS
        call cpu_time(time_start)
#endif
#ifndef NO_COMPUTE
        ! 2. Run the time/state nested loops, copying only time and state
        do state = ST_VELNW__BONDV1_INIT_UVW, ST_PRESS_BOUNDP

            state_ptr(1)=state
            call oclWrite1DIntArrayBuffer(state_ptr_buf,state_ptr_sz, state_ptr)
            select case (state)
! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_VELNW__BONDV1_INIT_UVW)
#ifdef TIMINGS
!                    call cpu_time(ktimestamp(ST_VELNW__BONDV1_INIT_UVW))
#endif
                    oclGlobalRange=(ip+1)*jp*kp
                    oclLocalRange=0
#if ICAL == 0
                    call oclWrite1DIntArrayBuffer(n_ptr_buf,n_ptr_sz, n_ptr)
#endif
#ifdef EXTERNAL_WIND_PROFILE
                    call oclWrite2DFloatArrayBuffer(wind_profile_buf,wind_profile_sz,wind_profile)
#endif

                    call runOcl(oclGlobalRange,oclLocalRange,exectime)
                    
                    ! Read the inner halos, those are going to be exchanged ot the other nodes
                    call read_halos_out(t_in, r_in, b_in, l_in, &
                        v_dim, h_h, h_d, c_w, &
                        t_in_buf, r_in_buf, b_in_buf, l_in_buf, t_in_sz, r_in_sz, b_in_sz, l_in_sz)
                        
                    ! Here we receive the halos from the other nodes
                    ! Since there are no other nodes for now, we just read the missing halos ourselves
                    call read_ext_halos_out(tl_in, tr_in, br_in, bl_in, &
                        v_dim, h_h, h_d, c_w, &
                        tl_in_buf, tr_in_buf, br_in_buf, bl_in_buf, tl_in_sz, tr_in_sz, br_in_sz, bl_in_sz)
                        
                    !========= Test ========!
                    ! Wipe the halos between the reading and writing
                    ! Just to make sure thing are working correctly
                    ! To be removed later once tested
                    !test_halo = -1.0
                    !big_test_halo = -1.0
                    !call oclWrite4DFloatArrayBuffer(t_in_buf, t_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_in_buf, r_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_in_buf, b_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_in_buf, l_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tl_in_buf, tl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tr_in_buf, tr_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(br_in_buf, br_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(bl_in_buf, bl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_out_buf, r_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_out_buf, l_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(t_out_buf, t_out_sz, big_test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_out_buf, b_out_sz, big_test_halo)
                    
                        
                    ! Merge the halos that have been received from the other nodes into the outer halos
                    call merge_halos_in(tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, b_out, &
                        v_dim, h_w, h_h, h_d, c_w, &
                        t_out_buf, r_out_buf, b_out_buf, l_out_buf, t_out_sz, r_out_sz, b_out_sz, l_out_sz)
                    
#ifdef TIMINGS
                    ktimestamp(ST_INIT) = exectime
#endif
! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_BONDV1_CALC_UOUT) ! REDUCTION
#ifdef TIMINGS
!                    call cpu_time(ktimestamp(ST_BONDV1_CALC_UOUT))
#endif

#if KERNEL == GPU_KERNEL
                    oclGlobalRange = jp
                    oclLocalRange = jp
                    ngroups = jp
#elif KERNEL == CPU_KERNEL
                    oclGlobalRange = oclNunits*NTH ! OclNunits and then loop over j, at least for Mac
                    oclLocalRange = NTH
                    ngroups = oclNunits
#else
#error "Value for KERNEL not supported!"
#endif

                    call runOcl(oclGlobalRange,oclLocalRange,exectime)
                    
                    ! Read the inner halos, those are going to be exchanged ot the other nodes
                    call read_halos_out(t_in, r_in, b_in, l_in, &
                        v_dim, h_h, h_d, c_w, &
                        t_in_buf, r_in_buf, b_in_buf, l_in_buf, t_in_sz, r_in_sz, b_in_sz, l_in_sz)
                        
                    ! Here we receive the halos from the other nodes
                    ! Since there are no other nodes for now, we just read the missing halos ourselves
                    call read_ext_halos_out(tl_in, tr_in, br_in, bl_in, &
                        v_dim, h_h, h_d, c_w, &
                        tl_in_buf, tr_in_buf, br_in_buf, bl_in_buf, tl_in_sz, tr_in_sz, br_in_sz, bl_in_sz)
                        
                    !========= Test ========!
                    ! Wipe the halos between the reading and writing
                    ! Just to make sure thing are working correctly
                    ! To be removed later once tested
                    !test_halo = -1.0
                    !big_test_halo = -1.0
                    !call oclWrite4DFloatArrayBuffer(t_in_buf, t_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_in_buf, r_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_in_buf, b_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_in_buf, l_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tl_in_buf, tl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tr_in_buf, tr_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(br_in_buf, br_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(bl_in_buf, bl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_out_buf, r_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_out_buf, l_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(t_out_buf, t_out_sz, big_test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_out_buf, b_out_sz, big_test_halo)
                    
                        
                    ! Merge the halos that have been received from the other nodes into the outer halos
                    call merge_halos_in(tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, b_out, &
                        v_dim, h_w, h_h, h_d, c_w, &
                        t_out_buf, r_out_buf, b_out_buf, l_out_buf, t_out_sz, r_out_sz, b_out_sz, l_out_sz)
#ifdef TIMINGS
                    ktimestamp(ST_BONDV1_CALC_UOUT) = exectime
#endif

#if KERNEL == CPU_KERNEL

                    call oclRead1DFloatArrayBuffer(chunks_num_buf,chunks_num_sz,chunks_num)
                    call oclRead1DFloatArrayBuffer(chunks_denom_buf,chunks_denom_sz, chunks_denom)

                    aaa=0.0
                    bbb=0.0
                    do ii = 1,ngroups
                        aaa = max(aaa,chunks_num(ii))
                        bbb = min(bbb,chunks_denom(ii))
                    end do
                    val_ptr(1) = (aaa + bbb)*0.5
#endif
! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_BONDV1_CALC_UVW)
#ifdef TIMINGS
!                    call cpu_time(ktimestamp(ST_BONDV1_CALC_UVW))
#endif
                    oclGlobalRange=(kp*jp)+(kp+2)*(ip+2)+(ip+3)*(jp+3)
                    oclLocalRange=0
#if KERNEL == CPU_KERNEL
                    call oclWrite1DFloatArrayBuffer(val_ptr_buf,val_ptr_sz, val_ptr)
#endif
                    call runOcl(oclGlobalRange,oclLocalRange,exectime)
                    
                    ! Read the inner halos, those are going to be exchanged ot the other nodes
                    call read_halos_out(t_in, r_in, b_in, l_in, &
                        v_dim, h_h, h_d, c_w, &
                        t_in_buf, r_in_buf, b_in_buf, l_in_buf, t_in_sz, r_in_sz, b_in_sz, l_in_sz)
                        
                    ! Here we receive the halos from the other nodes
                    ! Since there are no other nodes for now, we just read the missing halos ourselves
                    call read_ext_halos_out(tl_in, tr_in, br_in, bl_in, &
                        v_dim, h_h, h_d, c_w, &
                        tl_in_buf, tr_in_buf, br_in_buf, bl_in_buf, tl_in_sz, tr_in_sz, br_in_sz, bl_in_sz)
                        
                    !========= Test ========!
                    ! Wipe the halos between the reading and writing
                    ! Just to make sure thing are working correctly
                    ! To be removed later once tested
                    !test_halo = -1.0
                    !big_test_halo = -1.0
                    !call oclWrite4DFloatArrayBuffer(t_in_buf, t_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_in_buf, r_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_in_buf, b_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_in_buf, l_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tl_in_buf, tl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tr_in_buf, tr_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(br_in_buf, br_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(bl_in_buf, bl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_out_buf, r_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_out_buf, l_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(t_out_buf, t_out_sz, big_test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_out_buf, b_out_sz, big_test_halo)
                    
                        
                    ! Merge the halos that have been received from the other nodes into the outer halos
                    call merge_halos_in(tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, b_out, &
                        v_dim, h_w, h_h, h_d, c_w, &
                        t_out_buf, r_out_buf, b_out_buf, l_out_buf, t_out_sz, r_out_sz, b_out_sz, l_out_sz)
#ifdef TIMINGS
                    ktimestamp(ST_BONDV1_CALC_UVW) = exectime
#endif

! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_VELFG__FEEDBF__LES_CALC_SM)
#ifdef TIMINGS
!                    call cpu_time(ktimestamp(ST_VELFG__FEEDBF__LES_CALC_SM))
#endif
!#define NEW_VELFG
#ifdef NEW_VELFG
                    oclGlobalRange=ip*jp*(kp+1)
#else
                    oclGlobalRange=ip*jp*kp
#endif
                    oclLocalRange=0
                    call runOcl(oclGlobalRange,oclLocalRange,exectime)
                    
                    ! Read the inner halos, those are going to be exchanged ot the other nodes
                    call read_halos_out(t_in, r_in, b_in, l_in, &
                        v_dim, h_h, h_d, c_w, &
                        t_in_buf, r_in_buf, b_in_buf, l_in_buf, t_in_sz, r_in_sz, b_in_sz, l_in_sz)
                        
                    ! Here we receive the halos from the other nodes
                    ! Since there are no other nodes for now, we just read the missing halos ourselves
                    call read_ext_halos_out(tl_in, tr_in, br_in, bl_in, &
                        v_dim, h_h, h_d, c_w, &
                        tl_in_buf, tr_in_buf, br_in_buf, bl_in_buf, tl_in_sz, tr_in_sz, br_in_sz, bl_in_sz)
                        
                    !========= Test ========!
                    ! Wipe the halos between the reading and writing
                    ! Just to make sure thing are working correctly
                    ! To be removed later once tested
                    !test_halo = -1.0
                    !big_test_halo = -1.0
                    !call oclWrite4DFloatArrayBuffer(t_in_buf, t_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_in_buf, r_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_in_buf, b_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_in_buf, l_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tl_in_buf, tl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tr_in_buf, tr_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(br_in_buf, br_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(bl_in_buf, bl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_out_buf, r_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_out_buf, l_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(t_out_buf, t_out_sz, big_test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_out_buf, b_out_sz, big_test_halo)
                    
                        
                    ! Merge the halos that have been received from the other nodes into the outer halos
                    call merge_halos_in(tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, b_out, &
                        v_dim, h_w, h_h, h_d, c_w, &
                        t_out_buf, r_out_buf, b_out_buf, l_out_buf, t_out_sz, r_out_sz, b_out_sz, l_out_sz)
#ifdef TIMINGS
                    ktimestamp(ST_VELFG__FEEDBF__LES_CALC_SM) = exectime
#endif

! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_LES_BOUND_SM)
#ifdef TIMINGS
!                    call cpu_time(ktimestamp(ST_LES_BOUND_SM))
#endif
#if KERNEL == GPU_KERNEL
                    max_range = max(ip+3,jp+3,kp+2)
                    oclGlobalRange = max_range*max_range
                    oclLocalRange = max_range
#elif KERNEL == CPU_KERNEL
                    oclGlobalRange = (jp+3)*(kp+2) + (kp+2)*(ip+2) + (jp+3)*(ip+2)
                    oclLocalRange = 0
#else
#error "Value for KERNEL not supported!"
#endif
                    call runOcl(oclGlobalRange,oclLocalRange,exectime)
                    
                    ! Read the inner halos, those are going to be exchanged ot the other nodes
                    call read_halos_out(t_in, r_in, b_in, l_in, &
                        v_dim, h_h, h_d, c_w, &
                        t_in_buf, r_in_buf, b_in_buf, l_in_buf, t_in_sz, r_in_sz, b_in_sz, l_in_sz)
                        
                    ! Here we receive the halos from the other nodes
                    ! Since there are no other nodes for now, we just read the missing halos ourselves
                    call read_ext_halos_out(tl_in, tr_in, br_in, bl_in, &
                        v_dim, h_h, h_d, c_w, &
                        tl_in_buf, tr_in_buf, br_in_buf, bl_in_buf, tl_in_sz, tr_in_sz, br_in_sz, bl_in_sz)
                        
                    !========= Test ========!
                    ! Wipe the halos between the reading and writing
                    ! Just to make sure thing are working correctly
                    ! To be removed later once tested
                    !test_halo = -1.0
                    !big_test_halo = -1.0
                    !call oclWrite4DFloatArrayBuffer(t_in_buf, t_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_in_buf, r_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_in_buf, b_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_in_buf, l_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tl_in_buf, tl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tr_in_buf, tr_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(br_in_buf, br_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(bl_in_buf, bl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_out_buf, r_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_out_buf, l_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(t_out_buf, t_out_sz, big_test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_out_buf, b_out_sz, big_test_halo)
                    
                        
                    ! Merge the halos that have been received from the other nodes into the outer halos
                    call merge_halos_in(tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, b_out, &
                        v_dim, h_w, h_h, h_d, c_w, &
                        t_out_buf, r_out_buf, b_out_buf, l_out_buf, t_out_sz, r_out_sz, b_out_sz, l_out_sz)
#ifdef TIMINGS
                    ktimestamp(ST_LES_BOUND_SM) = exectime
#endif

! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_LES_CALC_VISC__ADAM)
#ifdef TIMINGS
!                    call cpu_time(ktimestamp(ST_LES_CALC_VISC__ADAM))
#endif
                    oclGlobalRange=ip*jp*kp
                    oclLocalRange=0

                    call runOcl(oclGlobalRange,oclLocalRange,exectime)
                    
                    ! Read the inner halos, those are going to be exchanged ot the other nodes
                    call read_halos_out(t_in, r_in, b_in, l_in, &
                        v_dim, h_h, h_d, c_w, &
                        t_in_buf, r_in_buf, b_in_buf, l_in_buf, t_in_sz, r_in_sz, b_in_sz, l_in_sz)
                        
                    ! Here we receive the halos from the other nodes
                    ! Since there are no other nodes for now, we just read the missing halos ourselves
                    call read_ext_halos_out(tl_in, tr_in, br_in, bl_in, &
                        v_dim, h_h, h_d, c_w, &
                        tl_in_buf, tr_in_buf, br_in_buf, bl_in_buf, tl_in_sz, tr_in_sz, br_in_sz, bl_in_sz)
                        
                    !========= Test ========!
                    ! Wipe the halos between the reading and writing
                    ! Just to make sure thing are working correctly
                    ! To be removed later once tested
                    !test_halo = -1.0
                    !big_test_halo = -1.0
                    !call oclWrite4DFloatArrayBuffer(t_in_buf, t_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_in_buf, r_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_in_buf, b_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_in_buf, l_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tl_in_buf, tl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tr_in_buf, tr_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(br_in_buf, br_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(bl_in_buf, bl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_out_buf, r_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_out_buf, l_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(t_out_buf, t_out_sz, big_test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_out_buf, b_out_sz, big_test_halo)
                    
                        
                    ! Merge the halos that have been received from the other nodes into the outer halos
                    call merge_halos_in(tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, b_out, &
                        v_dim, h_w, h_h, h_d, c_w, &
                        t_out_buf, r_out_buf, b_out_buf, l_out_buf, t_out_sz, r_out_sz, b_out_sz, l_out_sz)
#ifdef TIMINGS
                    ktimestamp(ST_LES_CALC_VISC__ADAM) = exectime
#endif
#ifndef NO_FILE_IO
                    if ((mod(n,1000) == 0 .or. n == nmax)) then
                        call oclRead4DFloatArrayBuffer(fgh_old_buf,fgh_old_sz,fgh_old)
                        call convert_from_fgh_old(fgh_old,foldo,goldo,holdo)
                        call write_fgh_old_to_file(foldo,goldo,holdo,ip,jp,kp,data21)
                    end if
#endif
! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_PRESS_RHSAV)
#ifdef TIMINGS
!                    call cpu_time(ktimestamp(ST_PRESS_RHSAV))
#endif
#if KERNEL == GPU_KERNEL
                    oclGlobalRange=ip*kp
                    oclLocalRange=ip
                    ngroups = ip
#elif KERNEL == CPU_KERNEL
                    oclGlobalRange = NTH*oclNunits
                    oclLocalRange = NTH
                    ngroups = oclNunits
#else
#error "Value for KERNEL not supported!"
#endif
                    call runOcl(oclGlobalRange,oclLocalRange,exectime)
                    
                    ! Read the inner halos, those are going to be exchanged ot the other nodes
                    call read_halos_out(t_in, r_in, b_in, l_in, &
                        v_dim, h_h, h_d, c_w, &
                        t_in_buf, r_in_buf, b_in_buf, l_in_buf, t_in_sz, r_in_sz, b_in_sz, l_in_sz)
                        
                    ! Here we receive the halos from the other nodes
                    ! Since there are no other nodes for now, we just read the missing halos ourselves
                    call read_ext_halos_out(tl_in, tr_in, br_in, bl_in, &
                        v_dim, h_h, h_d, c_w, &
                        tl_in_buf, tr_in_buf, br_in_buf, bl_in_buf, tl_in_sz, tr_in_sz, br_in_sz, bl_in_sz)
                        
                    !========= Test ========!
                    ! Wipe the halos between the reading and writing
                    ! Just to make sure thing are working correctly
                    ! To be removed later once tested
                    !test_halo = -1.0
                    !big_test_halo = -1.0
                    !call oclWrite4DFloatArrayBuffer(t_in_buf, t_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_in_buf, r_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_in_buf, b_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_in_buf, l_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tl_in_buf, tl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(tr_in_buf, tr_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(br_in_buf, br_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(bl_in_buf, bl_in_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(r_out_buf, r_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(l_out_buf, l_out_sz, test_halo)
                    !call oclWrite4DFloatArrayBuffer(t_out_buf, t_out_sz, big_test_halo)
                    !call oclWrite4DFloatArrayBuffer(b_out_buf, b_out_sz, big_test_halo)
                    
                        
                    ! Merge the halos that have been received from the other nodes into the outer halos
                    call merge_halos_in(tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, b_out, &
                        v_dim, h_w, h_h, h_d, c_w, &
                        t_out_buf, r_out_buf, b_out_buf, l_out_buf, t_out_sz, r_out_sz, b_out_sz, l_out_sz)
#ifdef TIMINGS
                    ktimestamp(ST_PRESS_RHSAV) = exectime
#endif

                    call oclRead1DFloatArrayBuffer(chunks_num_buf,chunks_num_sz,chunks_num)
                    call oclRead1DFloatArrayBuffer(chunks_denom_buf,chunks_denom_sz, chunks_denom)
                    ! Calc the average over the compute units
                    rhsav = 0.0
                    area = 0.0
                    do ii = 1,ngroups ! number of work groups
                        !                        print *,'RHSAV CHUNK',ii, chunks_num(ii)
                        !                        print *,ii, chunks_denom(ii)
                        rhsav = rhsav + chunks_num(ii)
                        area = area + chunks_denom(ii)
                    end do
                    rhsav = rhsav / area
                    val_ptr(1) = rhsav

! -----------------------------------------------------------------------------------------------------------------------------
                case (ST_PRESS_SOR)
#ifdef TIMINGS
!                    call cpu_time(ktimestamp(ST_PRESS_SOR))
                    sor_exectime=0.0
#endif
                    call oclRead4DFloatArrayBuffer(p_buf,p_sz,po)
                    !print *, '!!!!!====== TEST ======!!!!'
                    
                    call oclWrite1DFloatArrayBuffer(val_ptr_buf,val_ptr_sz, val_ptr)
                    iter = 0
                    sor = pjuge*1.1 ! just to be larger than pjuge
                    do while (sor > pjuge .and. iter < nmaxp)
#ifndef TWINNED_DOUBLE_BUFFER
                        iter = iter + 1
#else
                        iter = iter + 2
#endif
#ifndef APPROX_PAR_SOR
                        !  This is the "regular" approach
#ifdef TWINNED_DOUBLE_BUFFER
                        do nrd = 0,1
#else
                            do nrd = 0,2
#endif
#else
! def APPROX_PAR_SOR
                                nrd = 1
#endif

#ifdef TWINNED_DOUBLE_BUFFER
#ifdef FIXED_NTH
                                oclGlobalRange = NTH*oclNunits ! OclNunits + loop over kp+2
                                oclLocalRange = NTH
                                ngroups = oclNunits
#else
                                oclGlobalRange = (jp+2)*(kp+2)
                                oclLocalRange = jp+2
                                ngroups = kp+2
#endif
#else
! ndef TWINNED_DOUBLE_BUFFER ! This is the twinned double buffer approach
                                if (nrd < 2) then
#ifndef FIXED_NTH
                                    oclGlobalRange = kp*jp
                                    oclLocalRange =  jp !  iter order
                                    ngroups = kp !  iter order
#else
                                    oclGlobalRange = oclNunits*NTH !  OclNunits + loop over kp
                                    oclLocalRange = NTH
                                    ngroups = oclNunits
#endif
                                else
#ifndef FIXED_NTH
                                    oclGlobalRange = (ip+2)*(jp+2);
                                    oclLocalRange = jp+2
                                    ngroups = ip+2
#else
                                    oclGlobalRange = oclNunits*NTH !  OclNunits + loop over (ip+2)
                                    oclLocalRange = NTH
                                    ngroups = oclNunits
#endif
                                end if
#endif
#ifndef APPROX_PAR_SOR
                                n_ptr(1)=nrd
                                call oclWrite1DIntArrayBuffer(n_ptr_buf,n_ptr_sz, n_ptr)
#endif
                                call runOcl(oclGlobalRange,oclLocalRange,exectime)
                                
                                ! Read the inner halos, those are going to be exchanged ot the other nodes
                                call read_halos_out(t_in, r_in, b_in, l_in, &
                                    v_dim, h_h, h_d, c_w, &
                                    t_in_buf, r_in_buf, b_in_buf, l_in_buf, t_in_sz, r_in_sz, b_in_sz, l_in_sz)
                                    
                                ! Here we receive the halos from the other nodes
                                ! Since there are no other nodes for now, we just read the missing halos ourselves
                                call read_ext_halos_out(tl_in, tr_in, br_in, bl_in, &
                                    v_dim, h_h, h_d, c_w, &
                                    tl_in_buf, tr_in_buf, br_in_buf, bl_in_buf, tl_in_sz, tr_in_sz, br_in_sz, bl_in_sz)
                                    
                                !========= Test ========!
                                ! Wipe the halos between the reading and writing
                                ! Just to make sure thing are working correctly
                                ! To be removed later once tested
                                !test_halo = -1.0
                                !big_test_halo = -1.0
                                !call oclWrite4DFloatArrayBuffer(t_in_buf, t_in_sz, test_halo)
                                !call oclWrite4DFloatArrayBuffer(r_in_buf, r_in_sz, test_halo)
                                !call oclWrite4DFloatArrayBuffer(b_in_buf, b_in_sz, test_halo)
                                !call oclWrite4DFloatArrayBuffer(l_in_buf, l_in_sz, test_halo)
                                !call oclWrite4DFloatArrayBuffer(tl_in_buf, tl_in_sz, test_halo)
                                !call oclWrite4DFloatArrayBuffer(tr_in_buf, tr_in_sz, test_halo)
                                !call oclWrite4DFloatArrayBuffer(br_in_buf, br_in_sz, test_halo)
                                !call oclWrite4DFloatArrayBuffer(bl_in_buf, bl_in_sz, test_halo)
                                !call oclWrite4DFloatArrayBuffer(r_out_buf, r_out_sz, test_halo)
                                !call oclWrite4DFloatArrayBuffer(l_out_buf, l_out_sz, test_halo)
                                !call oclWrite4DFloatArrayBuffer(t_out_buf, t_out_sz, big_test_halo)
                                !call oclWrite4DFloatArrayBuffer(b_out_buf, b_out_sz, big_test_halo)
                                    
                                ! Merge the halos that have been received from the other nodes into the outer halos
                                call merge_halos_in(tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, b_out, &
                                    v_dim, h_w, h_h, h_d, c_w, &
                                    t_out_buf, r_out_buf, b_out_buf, l_out_buf, t_out_sz, r_out_sz, b_out_sz, l_out_sz)
#ifdef TIMINGS
                                sor_exectime = sor_exectime + exectime
#endif
                                if (nrd == 1) then
                                    !read back chunks_num. Only once!
                                    call oclRead1DFloatArrayBuffer(chunks_num_buf,chunks_num_sz,chunks_num)
                                    sor = 0.0
                                    do ii = 1,ngroups
                                        sor = sor + chunks_num(ii)
                                    end do
                                    sor = sqrt(sor) ! This is different from the original code but it is correct
                                end if
#ifndef APPROX_PAR_SOR
                            end do ! nrd
#endif

#ifndef APPROX_PAR_SOR
                            avg_iter = iter
#else
                            avg_iter = avg_iter + sum(chunks_denom)/ngroups
#endif
                        end do ! while iter
#ifdef TIMINGS
                    ktimestamp(ST_PRESS_SOR) = sor_exectime
#endif

#ifdef TIMINGS
                        !print *,"Exec time for SOR kernel = ",sor_exectime/1000," s"
!                        print '("Exec time for SOR kernel = ",f6.3," s")',sor_exectime/1000
#endif
! -----------------------------------------------------------------------------------------------------------------------------
                    case (ST_PRESS_PAV)
#ifdef TIMINGS
!                        call cpu_time(ktimestamp(ST_PRESS_PAV))
#endif
#if KERNEL == GPU_KERNEL
                        oclGlobalRange = ip*kp
                        oclLocalRange = ip
                        ngroups = kp
#elif KERNEL == CPU_KERNEL
                        oclGlobalRange = NTH*oclNunits ! so must loop over kp/oclNunits
                        oclLocalRange = NTH
                        ngroups = oclNunits
#else
#error "Value for KERNEL not supported!"
#endif

                        call runOcl(oclGlobalRange,oclLocalRange,exectime)
                        
                        ! Read the inner halos, those are going to be exchanged ot the other nodes
                        call read_halos_out(t_in, r_in, b_in, l_in, &
                            v_dim, h_h, h_d, c_w, &
                            t_in_buf, r_in_buf, b_in_buf, l_in_buf, t_in_sz, r_in_sz, b_in_sz, l_in_sz)
                            
                        ! Here we receive the halos from the other nodes
                        ! Since there are no other nodes for now, we just read the missing halos ourselves
                        call read_ext_halos_out(tl_in, tr_in, br_in, bl_in, &
                            v_dim, h_h, h_d, c_w, &
                            tl_in_buf, tr_in_buf, br_in_buf, bl_in_buf, tl_in_sz, tr_in_sz, br_in_sz, bl_in_sz)
                            
                        !========= Test ========!
                        ! Wipe the halos between the reading and writing
                        ! Just to make sure thing are working correctly
                        ! To be removed later once tested
                        !test_halo = -1.0
                        !big_test_halo = -1.0
                        !call oclWrite4DFloatArrayBuffer(t_in_buf, t_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(r_in_buf, r_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(b_in_buf, b_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(l_in_buf, l_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(tl_in_buf, tl_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(tr_in_buf, tr_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(br_in_buf, br_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(bl_in_buf, bl_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(r_out_buf, r_out_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(l_out_buf, l_out_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(t_out_buf, t_out_sz, big_test_halo)
                        !call oclWrite4DFloatArrayBuffer(b_out_buf, b_out_sz, big_test_halo)
                        
                            
                        ! Merge the halos that have been received from the other nodes into the outer halos
                        call merge_halos_in(tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, b_out, &
                            v_dim, h_w, h_h, h_d, c_w, &
                            t_out_buf, r_out_buf, b_out_buf, l_out_buf, t_out_sz, r_out_sz, b_out_sz, l_out_sz)
#ifdef TIMINGS
                    ktimestamp(ST_PRESS_PAV) = exectime
#endif

                        call oclRead1DFloatArrayBuffer(chunks_num_buf,chunks_num_sz,chunks_num)
                        call oclRead1DFloatArrayBuffer(chunks_denom_buf,chunks_denom_sz, chunks_denom)
                        ! Calc the average over the compute units
                        pav = 0.0
                        pco = 0.0
                        do ii = 1,ngroups
                            pav = pav + chunks_num(ii)
                            pco = pco + chunks_denom(ii)
                        end do
                        pav = pav / pco
                        val_ptr(1) = pav
! -----------------------------------------------------------------------------------------------------------------------------
                    case (ST_PRESS_ADJ)
#ifdef TIMINGS
!                        call cpu_time(ktimestamp(ST_PRESS_ADJ))
#endif
                        oclGlobalRange=ip*jp*kp
                        oclLocalRange=0

                        call oclWrite1DFloatArrayBuffer(val_ptr_buf,val_ptr_sz, val_ptr)
                        call runOcl(oclGlobalRange,oclLocalRange,exectime)
                        
                        ! Read the inner halos, those are going to be exchanged ot the other nodes
                        call read_halos_out(t_in, r_in, b_in, l_in, &
                            v_dim, h_h, h_d, c_w, &
                            t_in_buf, r_in_buf, b_in_buf, l_in_buf, t_in_sz, r_in_sz, b_in_sz, l_in_sz)
                            
                        ! Here we receive the halos from the other nodes
                        ! Since there are no other nodes for now, we just read the missing halos ourselves
                        call read_ext_halos_out(tl_in, tr_in, br_in, bl_in, &
                            v_dim, h_h, h_d, c_w, &
                            tl_in_buf, tr_in_buf, br_in_buf, bl_in_buf, tl_in_sz, tr_in_sz, br_in_sz, bl_in_sz)
                            
                        !========= Test ========!
                        ! Wipe the halos between the reading and writing
                        ! Just to make sure thing are working correctly
                        ! To be removed later once tested
                        !test_halo = -1.0
                        !big_test_halo = -1.0
                        !call oclWrite4DFloatArrayBuffer(t_in_buf, t_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(r_in_buf, r_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(b_in_buf, b_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(l_in_buf, l_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(tl_in_buf, tl_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(tr_in_buf, tr_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(br_in_buf, br_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(bl_in_buf, bl_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(r_out_buf, r_out_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(l_out_buf, l_out_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(t_out_buf, t_out_sz, big_test_halo)
                        !call oclWrite4DFloatArrayBuffer(b_out_buf, b_out_sz, big_test_halo)
                        
                            
                        ! Merge the halos that have been received from the other nodes into the outer halos
                        call merge_halos_in(tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, b_out, &
                            v_dim, h_w, h_h, h_d, c_w, &
                            t_out_buf, r_out_buf, b_out_buf, l_out_buf, t_out_sz, r_out_sz, b_out_sz, l_out_sz)
#ifdef TIMINGS
                    ktimestamp(ST_PRESS_ADJ) = exectime
#endif

! -----------------------------------------------------------------------------------------------------------------------------
                    case (ST_PRESS_BOUNDP)
#ifdef TIMINGS
!                        call cpu_time(ktimestamp(ST_PRESS_BOUNDP))
#endif

#if KERNEL == GPU_KERNEL
                        max_range = max((ip+2),(jp+2),(kp+2))
                        oclGlobalRange = max_range*max_range
                        oclLocalRange = max_range
#elif KERNEL == CPU_KERNEL
                        oclGlobalRange = (jp+2)*(kp+2) + (kp+2)*(ip+2) + (jp+2)*(ip+2)
                        oclLocalRange = 0
#else
#error "Value for KERNEL not supported!"
#endif
                        call runOcl(oclGlobalRange,oclLocalRange,exectime)
                        
                        ! Read the inner halos, those are going to be exchanged ot the other nodes
                        call read_halos_out(t_in, r_in, b_in, l_in, &
                            v_dim, h_h, h_d, c_w, &
                            t_in_buf, r_in_buf, b_in_buf, l_in_buf, t_in_sz, r_in_sz, b_in_sz, l_in_sz)
                            
                        ! Here we receive the halos from the other nodes
                        ! Since there are no other nodes for now, we just read the missing halos ourselves
                        call read_ext_halos_out(tl_in, tr_in, br_in, bl_in, &
                            v_dim, h_h, h_d, c_w, &
                            tl_in_buf, tr_in_buf, br_in_buf, bl_in_buf, tl_in_sz, tr_in_sz, br_in_sz, bl_in_sz)
                            
                        !========= Test ========!
                        ! Wipe the halos between the reading and writing
                        ! Just to make sure thing are working correctly
                        ! To be removed later once tested
                        !test_halo = -1.0
                        !big_test_halo = -1.0
                        !call oclWrite4DFloatArrayBuffer(t_in_buf, t_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(r_in_buf, r_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(b_in_buf, b_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(l_in_buf, l_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(tl_in_buf, tl_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(tr_in_buf, tr_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(br_in_buf, br_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(bl_in_buf, bl_in_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(r_out_buf, r_out_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(l_out_buf, l_out_sz, test_halo)
                        !call oclWrite4DFloatArrayBuffer(t_out_buf, t_out_sz, big_test_halo)
                        !call oclWrite4DFloatArrayBuffer(b_out_buf, b_out_sz, big_test_halo)
                        
                            
                        ! Merge the halos that have been received from the other nodes into the outer halos
                        call merge_halos_in(tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, b_out, &
                            v_dim, h_w, h_h, h_d, c_w, &
                            t_out_buf, r_out_buf, b_out_buf, l_out_buf, t_out_sz, r_out_sz, b_out_sz, l_out_sz)
#ifdef TIMINGS
                    ktimestamp(ST_PRESS_BOUNDP) = exectime
#endif

#ifdef VERBOSE
                        if (mod(n-1,20) == 0) then
                            print *, '=mac= time step, iteration step, conv =',n,iter,sor
                        end if

                        if (mod(n-1,20) == 0) then
                            call oclRead4DFloatArrayBuffer(uvw_buf,uvw_sz,uvw)
                            call oclRead4DFloatArrayBuffer(p_buf,p_sz,po)
                            ! convert to old format
                            call convert_from_uvw(uvw,uo,vo,wo)
                            do k = 1,kp,10
                                write(6,*) 'Inflow=',k,'u=',uo(0,jp/2,k),vo(0,jp/2,k) ,wo(0,jp/2,k),'p=',po(0,0,jp/2,k)
                            end do
                            do k = 1,kp,10
                                write(6,*) 'Urbanflow=',k,'u=',uo(ip/2,jp/2,k),vo(ip/2,jp/2,k) ,wo(ip/2,jp/2,k),'p=',po(0,ip/2,jp/2,k)
                            end do
                            !
                            cflu = 0.
                            cflv = 0.
                            cflw = 0.
                            do k = 1,kp
                                do j = 1,jp
                                    do i = 1,ip
                                        cflu = amax1(cflu,abs(uo(i,j,k)*dt/dx1(i)))
                                        cflv = amax1(cflv,abs(vo(i,j,k)*dt/dy1(j)))
                                        cflw = amax1(cflw,abs(wo(i,j,k)*dt/dzn(k)))
                                    end do
                                end do
                            end do
                        end if

                        if (mod(n-1,20) == 0) then
                            write(6,*) 'Check_CFL,u*dt/dx,v*dt/dy,w*dt/dz=',cflu,cflv,cflw
                        end if
#endif
                        !
                        if ((mod(n,1000) == 0 .or. n == nmax)) then

                            ! read back results and write to file
#ifndef VERBOSE
                            call oclRead4DFloatArrayBuffer(uvw_buf,uvw_sz,uvw)
                            call convert_from_uvw(uvw,uo,vo,wo)
#endif
                            call oclRead4DFloatArrayBuffer(uvwsum_buf,uvwsum_sz,uvwsum)
                            call convert_from_fgh(uvwsum,usumo,vsumo,wsumo)
                            call oclRead4DFloatArrayBuffer(p_buf,p_sz,po)

                            nn = n/1000
                            print *, 'timestep: ',nn,' pressure at centre: ',po(0,ip/2,jp/2,kp/2), &
                                'vel at centre: ', &
                                uo(ip/2,jp/2,kp/2),vo(ip/2,jp/2,kp/2),wo(ip/2,jp/2,kp/2)
#ifdef USE_NETCDF_OUTPUT
                            call write_to_netcdf_file(p,u,v,w,usum,vsum,wsum,nn)
#endif
#ifndef NO_FILE_IO
                            call write_uvw_p_uvwsum_to_file(uo,vo,wo,po,usumo,vsumo,wsumo,ip,jp,kp,data20)
#endif
                        end if
                end select
            end do ! states loop
#endif
! END of NO_COMPUTE
#ifdef TIMINGS
            call cpu_time(time_stop)
#ifndef NO_COMPUTE
            do itim=1, ST_PRESS_BOUNDP
                print '("Time for state ",i2," = ",f6.3," s")',itim,ktimestamp(itim)/1000 ! -ktimestamp(itim-1)
            end do
            print '("Total Kernel time = ",f6.3," s")',time_stop-time_start
#else
            print '("Time for no compute  = ",f6.3," s")',time_stop-time_start
#endif
#endif
        end subroutine run_LES_kernel
        
        
        ! Write the inner halos received in the outer halos
        subroutine merge_halos_in(tl_in, t_in, tr_in, r_in, br_in, b_in, bl_in, l_in, t_out, b_out, &
            v_dim, h_w, h_h, h_d, c_w, &
            t_out_buf, r_out_buf, b_out_buf, l_out_buf, t_out_sz, r_out_sz, b_out_sz, l_out_sz)
            integer, intent (in) :: v_dim, h_w, h_h, h_d, c_w
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(in)  :: tl_in
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(in)  :: t_in
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(in)  :: tr_in
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(in)  :: r_in
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(in)  :: br_in
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(in)  :: b_in
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(in)  :: bl_in
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(in)  :: l_in
            real(kind=4), dimension(v_dim, c_w + 2 * h_w, h_h, h_d), intent(out)  :: t_out
            real(kind=4), dimension(v_dim, c_w + 2 * h_w, h_h, h_d), intent(out)  :: b_out
            integer :: i_axis, j_axis, k_axis
            
            ! OpenCL buffer declarations
            integer(8) :: t_out_buf
            integer(8) :: r_out_buf
            integer(8) :: b_out_buf
            integer(8) :: l_out_buf
            ! OpenCL buffer size declarations
            integer, dimension(4):: t_out_sz
            integer, dimension(4):: r_out_sz
            integer, dimension(4):: b_out_sz
            integer, dimension(4):: l_out_sz
            
            ! Write the correct data inside the buffers
            do j_axis = 1, h_h
                do k_axis = 1, h_d
                    do i_axis = (c_w - h_w), c_w ! Left side
                        t_out(1:v_dim, i_axis - (c_w - h_w) + 1, j_axis, k_axis) = tl_in(1:v_dim, i_axis, j_axis, k_axis)
                        b_out(1:v_dim, i_axis - (c_w - h_w) + 1, j_axis, k_axis) = bl_in(1:v_dim, i_axis, j_axis, k_axis)
                    end do !i
                    do i_axis = 1, c_w ! Center
                        t_out(1:v_dim, i_axis + h_w, j_axis, k_axis) = t_in(1:v_dim, i_axis, j_axis, k_axis)
                        b_out(1:v_dim, i_axis + h_w, j_axis, k_axis) = b_in(1:v_dim, i_axis, j_axis, k_axis)
                    end do !i
                    do i_axis= 1, h_w ! Right side
                        t_out(1:v_dim, c_w + i_axis, j_axis, k_axis) = tr_in(1:v_dim, i_axis, j_axis, k_axis)
                        b_out(1:v_dim, c_w + i_axis, j_axis, k_axis) = br_in(1:v_dim, i_axis, j_axis, k_axis)
                    end do !i
                end do !k
            end do !j
            
            ! Write the buffers to device
            call oclWrite4DFloatArrayBuffer(t_out_buf, t_out_sz, t_out)
            call oclWrite4DFloatArrayBuffer(r_out_buf, r_out_sz, r_in)
            call oclWrite4DFloatArrayBuffer(b_out_buf, b_out_sz, b_out)
            call oclWrite4DFloatArrayBuffer(l_out_buf, l_out_sz, l_in)
                
        end subroutine merge_halos_in
        
        ! Read the inner halos from the device to the host
        subroutine read_halos_out(t_in, r_in, b_in, l_in, &
            v_dim, h_h, h_d, c_w, &
            t_in_buf, r_in_buf, b_in_buf, l_in_buf, t_in_sz, r_in_sz, b_in_sz, l_in_sz)
            integer, intent (in) :: v_dim, h_h, h_d, c_w
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(out)  :: t_in
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(out)  :: r_in
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(out)  :: b_in
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(out)  :: l_in
            
            ! OpenCL buffer declarations
            integer(8) :: t_in_buf
            integer(8) :: r_in_buf
            integer(8) :: b_in_buf
            integer(8) :: l_in_buf
            ! OpenCL buffer size declarations
            integer, dimension(4):: t_in_sz
            integer, dimension(4):: r_in_sz
            integer, dimension(4):: b_in_sz
            integer, dimension(4):: l_in_sz
            
            ! Read halos from device
            call oclRead4DFloatArrayBuffer(t_in_buf, t_in_sz, t_in)
            call oclRead4DFloatArrayBuffer(r_in_buf, r_in_sz, r_in)
            call oclRead4DFloatArrayBuffer(b_in_buf, b_in_sz, b_in)
            call oclRead4DFloatArrayBuffer(l_in_buf, l_in_sz, l_in)
            
        end subroutine read_halos_out
        
        ! Read the exterior inner halos from the device to the host
        ! Test function for a single node
        subroutine read_ext_halos_out(tl_in, tr_in, br_in, bl_in, &
            v_dim, h_h, h_d, c_w, &
            tl_in_buf, tr_in_buf, br_in_buf, bl_in_buf, tl_in_sz, tr_in_sz, br_in_sz, bl_in_sz)
            integer, intent (in) :: v_dim, h_h, h_d, c_w
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(out)  :: tl_in
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(out)  :: tr_in
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(out)  :: br_in
            real(kind=4), dimension(v_dim, c_w, h_h, h_d), intent(out)  :: bl_in
            
            ! OpenCL buffer declarations
            integer(8) :: tl_in_buf
            integer(8) :: tr_in_buf
            integer(8) :: br_in_buf
            integer(8) :: bl_in_buf
            ! OpenCL buffer size declarations
            integer, dimension(4):: tl_in_sz
            integer, dimension(4):: tr_in_sz
            integer, dimension(4):: br_in_sz
            integer, dimension(4):: bl_in_sz
            
            ! Read halos from device
            call oclRead4DFloatArrayBuffer(tl_in_buf, tl_in_sz, tl_in)
            call oclRead4DFloatArrayBuffer(tr_in_buf, tr_in_sz, tr_in)
            call oclRead4DFloatArrayBuffer(br_in_buf, br_in_sz, br_in)
            call oclRead4DFloatArrayBuffer(bl_in_buf, bl_in_sz, bl_in)
            
        end subroutine read_ext_halos_out
        
        
        
            
        ! Test functions
        ! Print the halos if different
        subroutine compare_halos(halo1, halo2, lb, ub)
            use module_LES_tests
            
            integer, dimension(4) :: ub
            integer, dimension(4) :: lb
            real(kind=4), dimension(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)) :: halo1, halo2
            call compare_4D_real_arrays (halo1, halo2, lb, ub)
        end subroutine
     
        ! --------------------------------------------------------------------------------
        ! --------------------------------------------------------------------------------
        ! Auxiliary subroutines for file I/O
        ! --------------------------------------------------------------------------------
! TODO: could refactor into separate module
        subroutine write_uvw_p_to_file(u,v,w,p,im,jm,km,data20)
            use params_common_sn
            real(kind=4), dimension(0:1,0:ip+2,0:jp+2,0:kp+1)  :: p
            real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: u
            real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: v
            real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1)  :: w
            integer :: im, jm, km
            character(len=70) :: data20
            open(unit=20,file=data20,form='unformatted',status='unknown')
            write(20) (((u(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((v(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((w(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((p(0,i,j,k),i=1,im),j=1,jm),k=1,km)
            close(unit=20)
        end subroutine write_uvw_p_to_file

        subroutine write_uvw_p_uvwsum_to_file(u,v,w,p,usum,vsum,wsum,im,jm,km,data20)
            use params_common_sn
            real(kind=4), dimension(0:1,0:ip+2,0:jp+2,0:kp+1), intent(In)  :: p
            real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1), intent(In)  :: u
            real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1) , intent(In) :: v
            real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1), intent(In)  :: w
            real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: usum
            real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: vsum
            real(kind=4), dimension(0:ip,0:jp,0:kp) , intent(In) :: wsum

            integer :: im, jm, km
            character(len=70) :: data20
            open(unit=20,file=data20,form='unformatted',status='unknown')
            write(20) (((u(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((v(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((w(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((p(0,i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((usum(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((vsum(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(20) (((wsum(i,j,k),i=1,im),j=1,jm),k=1,km)
            close(unit=20)
        end subroutine write_uvw_p_uvwsum_to_file

        subroutine write_fgh_old_to_file(fold,gold,hold,im,jm,km,data21)
            use params_common_sn
            real(kind=4), dimension(ip,jp,kp)  :: fold
            real(kind=4), dimension(ip,jp,kp)  :: gold
            real(kind=4), dimension(ip,jp,kp)  :: hold
            integer :: im, jm, km
            character(len=70) :: data21

            open(unit=21,file=data21,form='unformatted',status='unknown')
            write(21) (((fold(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(21) (((gold(i,j,k),i=1,im),j=1,jm),k=1,km)
            write(21) (((hold(i,j,k),i=1,im),j=1,jm),k=1,km)
            close(unit=21)
        end subroutine write_fgh_old_to_file
        
    
    end module module_LES_combined_ocl

