! #define IFBF 0
! #define IANIME 0
      program main
#ifdef USE_NETCDF_OUTPUT
        use module_LES_write_netcdf
#endif
        use module_init 
        use module_grid 
        use module_set 
        use module_timdata 
        use module_aveflow 
        use module_timseris 
        use module_ifdata 
#if IANIME == 1      
        use module_anime 
#endif      
! OpenCL modules
        use module_LES_tests

        use module_LES_combined_ocl
!        use module_velFG_ocl
        use module_velnw__bondv1_init_uvw_ocl
        use module_bondv1_ocl
!#if IFBF == 1
!        use module_feedbf_ocl
!#endif
        use module_les_ocl
!        use module_adam_ocl
        use module_press_ocl

! Reference modules
        use module_velnw
        use module_bondv1
        use module_velFG 
#if IFBF == 1      
        use module_feedbf
#endif      
        use module_les
        use module_adam
        use module_press
        
        use common_sn 
        real(kind=4) :: alpha
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1)  :: amask1
        real(kind=4), dimension(ip,jp,kp)  :: avel
        real(kind=4), dimension(ip,jp,kp)  :: avep
        real(kind=4), dimension(ip,jp,kp)  :: avesm
        real(kind=4), dimension(ip,jp,kp)  :: avesmsm
        real(kind=4), dimension(ip,kp)  :: avesu
        real(kind=4), dimension(ip,kp)  :: avesuu
        real(kind=4), dimension(ip,kp)  :: avesv
        real(kind=4), dimension(ip,kp)  :: avesvv
        real(kind=4), dimension(ip,kp)  :: avesw
        real(kind=4), dimension(ip,kp)  :: avesww
        real(kind=4), dimension(ip,jp,kp)  :: aveu
        real(kind=4), dimension(ip,jp,kp)  :: aveuu
        real(kind=4), dimension(ip,jp,kp)  :: avev
        real(kind=4), dimension(ip,jp,kp)  :: avevv
        real(kind=4), dimension(ip,jp,kp)  :: avew
        real(kind=4), dimension(ip,jp,kp)  :: aveww
        real(kind=4) :: beta
        real(kind=4), dimension(-1:ip+1,0:jp+1,0:kp+1)  :: bmask1
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: cmask1
        real(kind=4), dimension(ip,jp,kp)  :: cn1
        real(kind=4), dimension(ip)  :: cn2l
        real(kind=4), dimension(ip)  :: cn2s
        real(kind=4), dimension(jp)  :: cn3l
        real(kind=4), dimension(jp)  :: cn3s
        real(kind=4), dimension(kp)  :: cn4l
        real(kind=4), dimension(kp)  :: cn4s
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: cov1,cov1o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov2,cov2o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov3,cov3o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov4,cov4o
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: cov5,cov5o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov6,cov6o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov7,cov7o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov8,cov8o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: cov9,cov9o
        character(len=70) :: data10
        character(len=70) :: data11
        character(len=70) :: data20
        character(len=70) :: data21
        character(len=70) :: data22
        character(len=70) :: data23
        character(len=70) :: data24
        character(len=70) :: data25
        character(len=70) :: data26
        character(len=70) :: data27
        character(len=70) :: data30
        character(len=70) :: data31
        real(kind=4), dimension(kp)  :: delx1
        real(kind=4), dimension(0:ip,jp,kp)  :: dfu1
        real(kind=4), dimension(ip,0:jp,kp)  :: dfv1
        real(kind=4), dimension(ip,jp,kp)  :: dfw1
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: diu1,diu1o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu2,diu2o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu3,diu3o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu4,diu4o
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: diu5,diu5o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu6,diu6o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu7,diu7o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu8,diu8o
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: diu9,diu9o
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1)  :: dmask1
        real(kind=4) :: dt
        real(kind=4), dimension(-1:ip+1)  :: dx1
        real(kind=4), dimension(0:ip)  :: dxl
        real(kind=4), dimension(0:ip)  :: dxs
        real(kind=4), dimension(0:jp+1)  :: dy1
        real(kind=4), dimension(0:jp)  :: dyl
        real(kind=4), dimension(0:jp)  :: dys
        real(kind=4), dimension(-1:kp+2)  :: dzn
        real(kind=4), dimension(-1:kp+2)  :: dzs
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: f
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: fo
        real(kind=4), dimension(ip,jp,kp)  :: fghold
        real(kind=4), dimension(ip,jp,kp)  :: fold
        real(kind=4), dimension(ip,jp,kp)  :: foldo
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: fx
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: fxo
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: fy
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: fyo
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: fz
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: fzo
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: g
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: go
        real(kind=4), dimension(ip,jp,kp)  :: gold
        real(kind=4), dimension(ip,jp,kp)  :: goldo
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: h
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: ho
        real(kind=4), dimension(ip,jp,kp)  :: hold
        real(kind=4), dimension(ip,jp,kp)  :: holdo
        integer :: ianime
        integer :: ical
        integer :: ifbf
        integer :: im
        integer :: jm
        integer :: km
        integer :: n, nn
        integer :: n0
        integer :: n1
        integer :: nmax
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: nou1
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou2
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou3
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou4
        real(kind=4), dimension(-1:ip+2,0:jp+2,0:kp+2)  :: nou5
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou6
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou7
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou8
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+2)  :: nou9
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1)  :: po
        real(kind=4), dimension(0:ip+2,0:jp+2,0:kp+1)  :: p
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1)  :: rhs
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1)  :: rhso
        real(kind=4) :: ro
        real(kind=4), dimension(-1:ip+1,-1:jp+1,0:kp+1)  :: sm
        real(kind=4), dimension(-1:ip+1,-1:jp+1,0:kp+1)  :: smo
        real(kind=4) :: time
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: uo
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: u
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: usum
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: usumo
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: vsum
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: vsumo
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: wsum
        real(kind=4), dimension(0:ip,0:jp,0:kp)  :: wsumo
        real(kind=4), dimension(ip,jp,kp)  :: uwfx
        real(kind=4), dimension(ip,kp)  :: uwfxs
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: vo
        real(kind=4), dimension(0:ip+1,-1:jp+1,0:kp+1)  :: v
        real(kind=4) :: vn
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1)  :: wo
        real(kind=4), dimension(0:ip+1,-1:jp+1,-1:kp+1)  :: w
        real(kind=4), dimension(kp+2)  :: z2
        real(kind=4), dimension(-1:ipmax+1,-1:jpmax+1)  :: zbm
        real(kind=4) :: uout
! -----------------------------------------------------------------------
! 
      print * , 'Calling original routines set, grid, timdata, init, ifdata'

      call set(data10,data11,data20,data21,data22,data23,data24,data25,data26,data27,data30,data31, &
      im,jm,km,ifbf,ianime,ical,n0,n1,nmax,dt,ro,vn,alpha,beta)
      call grid(dx1,dxl,dy1,dyl,z2,dzn,dzs,dxs,dys)
      call timdata()
      call init(km,jm,im,u,v,w,p,ifbf,cn2s,dxs,cn2l,cn3s,dys,cn3l,dzs,cn4s,cn4l,cn1,amask1,bmask1, &
      cmask1,dmask1,zbm,z2,dzn)
      call ifdata(ical,data30,n,time,u,im,jm,km,v,w,p,usum,vsum,wsum,data31,fold,gold,hold,fghold, &
      ifbf,delx1,dx1,dy1,dzn,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,sm,f,g,h,z2,dt,dxs,cov1, &
      cov2,cov3,dfu1,vn,cov4,cov5,cov6,dfv1,cov7,cov8,cov9,dfw1,dzs,nou1,nou5,nou9,nou2,nou3,nou4, &
      nou6,nou7,nou8,bmask1,cmask1,dmask1,alpha,beta,fx,fy,fz,amask1,zbm)

      print * , '-----------------------------------------------'  
      print * , '---------------- Running Tests ----------------'  
      print * , '-----------------------------------------------'  

      print * , "Test uvw conversion: "
      call test_uvw_conversion(u,v,w)
      print * , "Test fgh conversion: "
      call test_fgh_conversion(f,g,h)
      print * , "Test fgh_old conversion: "
      call test_fgh_old_conversion(fold,gold,hold)
      print * , "Test cov* conversion: "
      call test_9vec_conversion(cov1,cov2,cov3,cov4,cov5,cov6,cov7,cov8,cov9)
      print * , "Test cn* conversion: "
      call test_cn234ls_conversion(cn2l,cn2s,cn3l,cn3s,cn4l,cn4s)

      print * , ''
      p=1.1
      u=2.2
      v=3.3
      w=4.4
      fold=5.5
      gold=6.6
      hold=7.7
      print * , "-------- Testing  OpenCL init and data transfer --------"
      print * , "1. OpenCL init"
      call initialise_LES_kernel( &
            p,u,v,w,usum,vsum,wsum,f,g,h,fold,gold,hold, &
            diu1, diu2, diu3, diu4, diu5, diu6, diu7, diu8, diu9, &
            amask1, bmask1, cmask1, dmask1, &
            cn1, cn2l, cn2s, cn3l, cn3s, cn4l, cn4s, &
            rhs, sm, dxs, dys, dzs, dx1, dy1, dzn, z2, &
            dt, im, jm, km, nmax &
            )
      print * , "2. Running the kernel and reading back the data"
      call run_LES_kernel ( &
            data20, data21, &
            im, jm, km, &
            dt, dx1,dy1,dzn, &
            n, nmax &
            ,po,uo,vo,wo,foldo,goldo,holdo &
            )
      print * , "3. Comparing with the original data"
      print * , "3.1 Comparing p"
      call test_kernel_transfer_p(p,po)
      print * , "3.2 Comparing u,v,w"
      call test_kernel_transfer_uvw (u,v,w,uo,vo,wo)    
!      print * , "3.3 Comparing f,g,h"
!      call test_kernel_transfer_fgh (f,g,h,fo,go,ho)
      print * , "3.4 Comparing fold,gold,hold"
      call test_kernel_transfer_fghold (fold,gold,hold,foldo,goldo,holdo)
      print * , 'Status: kernel works correctly in pass-through mode'

! The aim now is to validate if my new velFG kernel works correctly.
! To do so, we must call first init, then velFG

      ! WV: added for sanity check
      print *, 'timestep: ',n,' pressure at centre: ',p(ip/2,jp/2,kp/2), &
          'vel at centre: ', &
          u(ip/2,jp/2,kp/2),v(ip/2,jp/2,kp/2),w(ip/2,jp/2,kp/2), &
          f(ip/2,jp/2,kp/2),g(ip/2,jp/2,kp/2),h(ip/2,jp/2,kp/2)

#ifdef USE_NETCDF_OUTPUT
      call init_netcdf_file
#endif
!      nmax = 2000 !DEBUG

      n = nmax
      time = float(n-1)*dt
! -------calculate turbulent flow--------c
        fo = f+0.0
        go = g+0.0
        ho = h+0.0

        uo = u+0.0
        vo = v+0.0
        wo = w+0.0
        print *, ''
!print *, '-------- Comparing velnw with F2C_ACC version --------'
print *, '-------- Comparing velnw with OpenCL version --------'
        call velnw(km,jm,im,p,ro,dxs,u,dt,f,dys,v,g,dzs,w,h)

!        call velnwC(km,jm,im,p,ro,dxs,uo,dt,f,dys,vo,g,dzs,wo,h)
!        call velnw_ocl(km,jm,im,p,ro,dxs,uo,dt,fo,dys,vo,go,dzs,wo,ho)
        call velnw__bondv1_init_uvw_ocl(p,uo,vo,wo,fo,go,ho,dxs,dys,dzs,dzn,z2,n,ro,dt,im,jm,km)
        
! u,v,w are the only outputs for velnw
    print * , ''
!    print * , 'u'
!    call compare_3D_real_arrays(u,uo,lbound(u),ubound(u))
!    print * , 'v'
!    call compare_3D_real_arrays(v,vo,lbound(v),ubound(v))
!    print * , 'w'
!    call compare_3D_real_arrays(w,wo,lbound(w),ubound(w))
print * , 'VELNW_UVWSUM: ',sum(u)+sum(v)+sum(w),sum(uo)+sum(vo)+sum(wo)
!        uo = u+0.0
!        vo = v+0.0
!        wo = w+0.0
print *, '-------- Comparing bondv1 with OpenCL version --------'
    
        call bondv1(jm,u,z2,dzn,v,w,km,ical,n,im,dt,dxs)

!        call bondv1C(jm,uo,z2,dzn,vo,wo,km,ical,n,im,dt,dxs)
        call bondv1_calc_uout_ocl(uo,vo,wo,uout,im,jm,km)
!        call bondv1_calc_uvw_ocl(uo,vo,wo,dxs,uout,dt,im,jm,km)
!        call bondv1_ocl(jm,uo,z2,dzn,vo,wo,km,uout,n,im,dt,dxs)
!        call bondv1_ocl(jm,uo,z2,dzn,vo,wo,km,ical,n,im,dt,dxs)

! u,v,w are the only outputs for velnw
!    print * , ''
!    print * , 'u'
!    call compare_3D_real_arrays(u,uo,lbound(u),ubound(u))
!    print * , 'v'
!    call compare_3D_real_arrays(v,vo,lbound(v),ubound(v))
!    print * , 'w'
!    call compare_3D_real_arrays(w,wo,lbound(w),ubound(w))
print * , 'BONDV1_UVWSUM: ',sum(u)+sum(v)+sum(w),sum(uo)+sum(vo)+sum(wo)

!        uo = u+0.0
!        vo = v+0.0
!        wo = w+0.0
!
        diu1o = diu1+0.0
        diu2o = diu2+0.0
        diu3o = diu3+0.0
        diu4o = diu4+0.0
        diu5o = diu5+0.0
        diu6o = diu6+0.0
        diu7o = diu7+0.0
        diu8o = diu8+0.0
        diu9o = diu9+0.0
        
        cov1o = cov1+0.0
        cov2o = cov2+0.0
        cov3o = cov3+0.0
        cov4o = cov4+0.0
        cov5o = cov5+0.0
        cov6o = cov6+0.0
        cov7o = cov7+0.0
        cov8o = cov8+0.0
        cov9o = cov9+0.0

        print *, ''
print *, '-------- Comparing velFG with OpenCL version --------'
        
        call velfg(km,jm,im,dx1,cov1,cov2,cov3,dfu1,diu1,diu2,dy1,diu3,dzn,vn,f,cov4,cov5,cov6,dfv1, &
      diu4,diu5,diu6,g,cov7,cov8,cov9,dfw1,diu7,diu8,diu9,dzs,h,nou1,u,nou5,v,nou9,w,nou2,nou3, &
      nou4,nou6,nou7,nou8)

!        call velFG_ocl(km,jm,im,dx1,cov1o,cov2o,cov3o,dfu1,diu1o,diu2o,dy1,diu3o,dzn,vn,fo,cov4o,cov5o,cov6o,dfv1, &
!      diu4o,diu5o,diu6o,go,cov7o,cov8o,cov9o,dfw1,diu7o,diu8o,diu9o,dzs,ho,nou1,uo,nou5,vo,nou9,wo,nou2,nou3, &
!      nou4,nou6,nou7,nou8)

!    print * , ''
!    print * , 'diu1'
!    call compare_3D_real_arrays(diu1,diu1o,lbound(diu1),ubound(diu1))
!    print * , 'diu2'
!    call compare_3D_real_arrays(diu2,diu2o,lbound(diu2),ubound(diu2))
!    print * , 'diu3'
!    call compare_3D_real_arrays(diu3,diu3o,lbound(diu3),ubound(diu3))
!    print * , 'diu4'
!    call compare_3D_real_arrays(diu4,diu4o,lbound(diu4),ubound(diu4))
!    print * , 'diu5'
!    call compare_3D_real_arrays(diu5,diu5o,lbound(diu5),ubound(diu5))
!    print * , 'diu6'
!    call compare_3D_real_arrays(diu6,diu6o,lbound(diu6),ubound(diu6))
!    print * , 'diu7'
!    call compare_3D_real_arrays(diu7,diu7o,lbound(diu7),ubound(diu7))
!    print * , 'diu8'
!    call compare_3D_real_arrays(diu8,diu8o,lbound(diu8),ubound(diu8))
!    print * , 'diu9'
!    call compare_3D_real_arrays(diu9,diu9o,lbound(diu9),ubound(diu9))
!
!    print * , ''
!    print * , 'f'
!    call compare_3D_real_arrays(f,fo,lbound(f),ubound(f))
!    print * , 'g'
!    call compare_3D_real_arrays(g,go,lbound(g),ubound(g))
!    print * , 'h'
!    call compare_3D_real_arrays(h,ho,lbound(h),ubound(h))
!    print * , 'VELFG_FGHSUM: ',sum(f)+sum(g)+sum(h),sum(fo)+sum(go)+sum(ho)
!#if IFBF == 1
!    print * , ''
!    print *, '-------- Comparing feedbf with OpenCL version --------'
!
!        fo = f+0.0
!        go = g+0.0
!        ho = h+0.0
!
!        usumo=usum+0.0
!        vsumo=vsum+0.0
!        wsumo=wsum+0.0
        call feedbf(km,jm,im,usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,fx,fy,fz,f,g,h)
!
!!        call feedbf_ocl(km,jm,im,usumo,u,bmask1,vsumo,v,cmask1,wsumo,w,dmask1,alpha,dt,beta,fx,fy,fz,fo,go,ho)
!
!    print * , ''
!    print * , 'usum'
!    call compare_3D_real_arrays(usum,usumo,lbound(usum),ubound(usum))
!    print * , 'vsum'
!    call compare_3D_real_arrays(vsum,vsumo,lbound(vsum),ubound(vsum))
!    print * , 'wsum'
!    call compare_3D_real_arrays(wsum,wsumo,lbound(wsum),ubound(wsum))
!
!    print * , ''
!    print * , 'f'
!    call compare_3D_real_arrays(f,fo,lbound(f),ubound(f))
!    print * , 'g'
!    call compare_3D_real_arrays(g,go,lbound(g),ubound(g))
!    print * , 'h'
!    call compare_3D_real_arrays(h,ho,lbound(h),ubound(h))
!        print * , 'FEEDBF_FGHSUM: ',sum(f)+sum(g)+sum(h),sum(fo)+sum(go)+sum(ho)
!!    print * , ''
!!    print * , 'fx'
!!    call compare_3D_real_arrays(fx,fxo,lbound(fx),ubound(fx))
!!    print * , 'fy'
!!    call compare_3D_real_arrays(fy,fyo,lbound(fy),ubound(fy))
!!    print * , 'fz'
!!    call compare_3D_real_arrays(fz,fzo,lbound(fz),ubound(fz))
!
!
!#endif
    print * , ''        
!    print *, '-------- Comparing les with OpenCL version --------'
!    print * , ''
    print *, '-------- Comparing velFG/feedbf/les_calc_sm with merged OpenCL version --------'

        fo = f+0.0
        go = g+0.0
        ho = h+0.0
        
        usumo=usum+0.0    
        vsumo=vsum+0.0    
        wsumo=wsum+0.0    

!        fo = f+0.0
!        go = g+0.0
!        ho = h+0.0
        smo = sm+0.0
        call les(km,delx1,dx1,dy1,dzn,jm,im,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,sm,f,g,h)
        
        !call lesC(km,delx1,dx1,dy1,dzn,jm,im,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,sm,fo,go,ho)
!        call les_calc_sm_ocl(km,delx1,dx1,dy1,dzn,jm,im,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,smo,fo,go,ho)

        call bondv1_calc_uvw__velfg__feedbf__les_calc_sm_ocl(fo,go,ho,uo,vo,wo,usumo,vsumo,wsumo,smo, &
        diu1o,diu2o,diu3o,diu4o,diu5o,diu6o,diu7o,diu8o,diu9o,dx1,dy1,dzn,dzs,dxs,uout,bmask1,cmask1,dmask1,dt, im, jm, km)
!        call merged_velfg_feedbf_les_calc_sm_ocl(fo,go,ho,uo,vo,wo,usumo,vsumo,wsumo,smo, &
!        diu1o,diu2o,diu3o,diu4o,diu5o,diu6o,diu7o,diu8o,diu9o,dx1,dy1,dzn,dzs, bmask1,cmask1,dmask1,dt, im, jm, km)

        call les_bound_sm_ocl(km,delx1,dx1,dy1,dzn,jm,im,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,smo,fo,go,ho)
!        call les_calc_visc_ocl(km,delx1,dx1,dy1,dzn,jm,im,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,smo,fo,go,ho)
!    print * , ''
!    print * , 'f'
!    call compare_3D_real_arrays(f,fo,lbound(f),ubound(f))
!    print * , 'g'
!    call compare_3D_real_arrays(g,go,lbound(g),ubound(g))
!    print * , 'h'
!    call compare_3D_real_arrays(h,ho,lbound(h),ubound(h))

    print * , ''
    print * , 'sm'
    call compare_3D_real_arrays(sm,smo,lbound(sm),ubound(sm))
print * , 'LES_FGHSUM: ',sum(f)+sum(g)+sum(h),sum(fo)+sum(go)+sum(ho)
    print * , ''   
!    print *, '-------- Comparing adam with OpenCL version --------'
    print *, '-------- Comparing adam/les_calc_visc with merged OpenCL version --------'
!        fo = f+0.0
!        go = g+0.0
!        ho = h+0.0
!        foldo = fold+0.0
!        goldo = gold+0.0
!        holdo = hold+0.0

        call adam(n,nmax,data21,fold,im,jm,km,gold,hold,fghold,f,g,h)

        !call adamC(n,nmax,data21,foldo,im,jm,km,goldo,holdo,fghold,fo,go,ho)
!        call adam_ocl(n,nmax,data21,foldo,im,jm,km,goldo,holdo,fghold,fo,go,ho)
        call les_calc_visc__adam_ocl(km,delx1,dx1,dy1,dzn,jm,im, &
        diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,smo, &
        foldo,goldo,holdo,fghold,fo,go,ho)

    print * , ''
    print * , 'f'
    call compare_3D_real_arrays(f,fo,lbound(f),ubound(f))
    print * , 'g'
    call compare_3D_real_arrays(g,go,lbound(g),ubound(g))
    print * , 'h'
    call compare_3D_real_arrays(h,ho,lbound(h),ubound(h))
    print * , 'ADAM_UVWSUM: ',sum(u)+sum(v)+sum(w),sum(uo)+sum(vo)+sum(wo)
    print * , 'ADAM_FGHSUM: ',sum(f)+sum(g)+sum(h),sum(fo)+sum(go)+sum(ho)
!    print * , ''
!    print * , 'fold'
!    call compare_3D_real_arrays(fold,foldo,lbound(fold),ubound(fold))
!    print * , 'gold'
!    call compare_3D_real_arrays(gold,goldo,lbound(gold),ubound(gold))
!    print * , 'hold'
!    call compare_3D_real_arrays(hold,holdo,lbound(hold),ubound(hold))

    print * , ''    
    print * , "U_SUM: ", sum(u) , sum(uo)
    print * , "V_SUM: ", sum(v), sum(vo)
    print * , "W_SUM: ", sum(w), sum(wo)
        print * , ''
    print * , "F_SUM: ", sum(f), sum(fo)
    print * , "G_SUM: ", sum(g), sum(go)
    print * , "H_SUM: ", sum(h), sum(ho)
    print * , ""
    print * , 'ADAM_PSUM: ',sum(p),sum(po)
    print *, '-------- Comparing press with OpenCL version --------'
        fo = f+0.0
        go = g+0.0
        ho = h+0.0
        po = p+0.0
        rhs = 0.0
        rhso = 0.0
        call press(km,jm,im,rhs,u,dx1,v,dy1,w,dzn,f,g,h,dt,cn1,cn2l,p,cn2s,cn3l,cn3s,cn4l,cn4s,n, &
      nmax,data20,usum,vsum,wsum)

!        call pressC(km,jm,im,rhs,u,dx1,v,dy1,w,dzn,fo,go,ho,dt,cn1,cn2l,po,cn2s,cn3l,cn3s,cn4l,cn4s,n, &
!      nmax,data20,usum,vsum,wsum)

!        call press_ocl(km,jm,im,rhs,u,dx1,v,dy1,w,dzn,fo,go,ho,dt,cn1,cn2l,po,cn2s,cn3l,cn3s,cn4l,cn4s,n, &
!      nmax,data20,usum,vsum,wsum)

        call press_rhsav_ocl(km,jm,im,rhso,u,dx1,v,dy1,w,dzn,fo,go,ho,dt,cn1,cn2l,po,cn2s,cn3l,cn3s,cn4l,cn4s,n,nmax,data20,usum,vsum,wsum)
        print *, "OCL: RHSAV: ", rhso(0,0,0)
        print *, "OCL: RHS_SUM: ", sum(rhso)
        print *, "OK on GPU before SOR"
        call press_sor_ocl(km,jm,im,rhso,u,dx1,v,dy1,w,dzn,fo,go,ho,dt,cn1,cn2l,po,cn2s,cn3l,cn3s,cn4l,cn4s,n,nmax,data20,usum,vsum,wsum)
        print *, "OCL: P_SUM_SOR: ", sum(po)
        print *, "OCL: SOR DONE!"
        print *, ""
        print *, "OK on GPU before PAV"
        call press_pav_ocl(km,jm,im,rhso,u,dx1,v,dy1,w,dzn,fo,go,ho,dt,cn1,cn2l,po,cn2s,cn3l,cn3s,cn4l,cn4s,n,nmax,data20,usum,vsum,wsum)
        print *, "OCL: PAV: ", rhso(0,0,0)
        print *, "OCL: P_SUM_PAV: ", sum(po)
        call press_adj_ocl(km,jm,im,rhso,u,dx1,v,dy1,w,dzn,fo,go,ho,dt,cn1,cn2l,po,cn2s,cn3l,cn3s,cn4l,cn4s,n,nmax,data20,usum,vsum,wsum)
        print *, "OCL: P_SUM_ADJ: ", sum(po)
        call press_boundp_ocl(km,jm,im,rhso,u,dx1,v,dy1,w,dzn,fo,go,ho,dt,cn1,cn2l,po,cn2s,cn3l,cn3s,cn4l,cn4s,n,nmax,data20,usum,vsum,wsum)
        print *, "OCL: P_SUM_BOUND: ", sum(po)

      if(mod(n,1000) == 0.or.n == nmax) then
        nn = n/1000
        print *, 'timestep: ',nn,' pressure at centre: ',p(ip/2,jp/2,kp/2), &
                'vel at centre: ', &
                u(ip/2,jp/2,kp/2),v(ip/2,jp/2,kp/2),w(ip/2,jp/2,kp/2)
#ifdef USE_NETCDF_OUTPUT
        call write_to_netcdf_file(p,u,v,w,usum,vsum,wsum,nn)
#endif
      end if


    print * , ''
    print * , 'f'
    call compare_3D_real_arrays(f,fo,lbound(f),ubound(f))
    print * , 'g'
    call compare_3D_real_arrays(g,go,lbound(g),ubound(g))
    print * , 'h'
    call compare_3D_real_arrays(h,ho,lbound(h),ubound(h))
    print * , ''
    print * , 'p'
    call compare_3D_real_arrays(p,po,lbound(p),ubound(p))
   print * , 'PRESS_UVWSUM: ',sum(u)+sum(v)+sum(w),sum(uo)+sum(vo)+sum(wo)
   print * , 'PRESS_FGHSUM: ',sum(f)+sum(g)+sum(h),sum(fo)+sum(go)+sum(ho)
   print * , 'PRESS_PSUM: ',sum(p),sum(po)

! -------data output ---------------------c
        call timseris(n,dt,u,v,w)

        call aveflow(n,n1,km,jm,im,aveu,avev,avew,avep,avel,aveuu,avevv,aveww,avesm,avesmsm,uwfx, &
      avesu,avesv,avesw,avesuu,avesvv,avesww,u,v,w,p,sm,nmax,uwfxs,data10,time,data11)
#if IANIME == 1
        call anime(n,n0,nmax,km,jm,im,dxl,dx1,dyl,dy1,z2,data22,data23,u,w,v,amask1)

#endif
!
#ifdef USE_NETCDF_OUTPUT
      call close_netcdf_file
#endif

      end program

