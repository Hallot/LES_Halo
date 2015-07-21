      program main
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

      print * , ':-----------------------------------------------:'
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
!      po = p+0.0
!      fo = f+0.0
!      go = g+0.0
!      ho = h+0.0
!      foldo = fold+0.0
!      goldo = gold+0.0
!      holdo = hold+0.0
!
!      uo = u+0.0
!      vo = v+0.0
!      wo = w+0.0
!      p=1.1
!      u=2.2
!      v=3.3
!      w=4.4
!      fold=5.5
!      gold=6.6
!      hold=7.7

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
      n=0
      call run_LES_kernel ( &
            data20, data21, &
            im, jm, km, &
            dt, dx1,dy1,dzn, &
            n, nmax &
            ,p,u,v,w,f,g,h,fold,gold,hold &
!            ,po,uo,vo,wo,foldo,goldo,holdo &
            )
      print * , "3. Comparing with the original data"
!      print * , "3.1 Comparing p"
!      call test_kernel_transfer_p(p,po)
!      print * , "3.2 Comparing u,v,w"
!      call test_kernel_transfer_uvw (u,v,w,uo,vo,wo)
!      print * , "3.3 Comparing fold,gold,hold"
!      call test_kernel_transfer_fghold (fold,gold,hold,foldo,goldo,holdo)
!      print * , "3.4 Comparing f,g,h"
!      call test_kernel_transfer_fgh (f,g,h,fo,go,ho)
!
!      print * , '3_UVWSUM: ',sum(u)+sum(v)+sum(w),sum(uo)+sum(vo)+sum(wo)
!      print * , '3_FGHSUM: ',sum(f)+sum(g)+sum(h),sum(fo)+sum(go)+sum(ho)
      time = float(n-1)*dt
! -------calculate turbulent flow--------c
        print *, ''
! ------------------------------------------------------------------------------------------------
print *, '4. Comparing OpenCL velnw__bondv1_init_uvw with velnw+bondv1'
!    print * , '4_USUM before velnw: ',sum(u),sum(uo)
    n=n0
    call velnw(km,jm,im,p,ro,dxs,u,dt,f,dys,v,g,dzs,w,h)
!    print * , '4_USUM after velnw: ',sum(u),sum(uo)
    call bondv1(jm,u,z2,dzn,v,w,km,ical,n,im,dt,dxs)
!print * , '4_USUM after bondv1: ',sum(u),sum(uo)
    call run_LES_kernel ( &
            data20, data21, &
            im, jm, km, &
            dt, dx1,dy1,dzn, &
            n, nmax &
            ,p,u,v,w,f,g,h,fold,gold,hold &
!            ,po,uo,vo,wo,foldo,goldo,holdo &
            )
!      print * , "4.1 Comparing u,v,w"
!! u,v,w are the only outputs for velnw
!    print * , ''
!    print * , 'u'
!    call compare_3D_real_arrays(u,uo,lbound(u),ubound(u))
!    print * , 'v'
!    call compare_3D_real_arrays(v,vo,lbound(v),ubound(v))
!    print * , 'w'
!    call compare_3D_real_arrays(w,wo,lbound(w),ubound(w))
!    print * , '4_UVWSUM: OCL:',sum(uo)+sum(vo)+sum(wo)
! ------------------------------------------------------------------------------------------------
    print *, '5. Comparing bondv1_calc_uout with bondv1'
    n=2
          call run_LES_kernel ( &
            data20, data21, &
            im, jm, km, &
            dt, dx1,dy1,dzn, &
            n, nmax &
            ,p,u,v,w,f,g,h,fold,gold,hold &
!            ,po,uo,vo,wo,foldo,goldo,holdo &
            )
! ------------------------------------------------------------------------------------------------
print *, '6. Comparing bondv1_calc_uvw'
n=3
 call run_LES_kernel ( &
            data20, data21, &
            im, jm, km, &
            dt, dx1,dy1,dzn, &
            n, nmax &
            ,p,u,v,w,f,g,h,fold,gold,hold &
!            ,po,uo,vo,wo,foldo,goldo,holdo &
            )
! ------------------------------------------------------------------------------------------------
    print *, '7. Comparing velfg__feedbf__les_calc_sm with velfg/feedbf/les'
    n=4
        call velfg(km,jm,im,dx1,cov1,cov2,cov3,dfu1,diu1,diu2,dy1,diu3,dzn,vn,f,cov4,cov5,cov6,dfv1, &
      diu4,diu5,diu6,g,cov7,cov8,cov9,dfw1,diu7,diu8,diu9,dzs,h,nou1,u,nou5,v,nou9,w,nou2,nou3, &
      nou4,nou6,nou7,nou8)
          call run_LES_kernel ( &
            data20, data21, &
            im, jm, km, &
            dt, dx1,dy1,dzn, &
            n, nmax &
            ,p,u,v,w,f,g,h,fold,gold,hold &
!            ,po,uo,vo,wo,foldo,goldo,holdo &
            )
! ------------------------------------------------------------------------------------------------
    print *, '8. Comparing velfg__feedbf__les_calc_sm with velfg/feedbf/les'
    n=5

        call feedbf(km,jm,im,usum,u,bmask1,vsum,v,cmask1,wsum,w,dmask1,alpha,dt,beta,fx,fy,fz,f,g,h)
        print *,'Calling les():'
        call les(km,delx1,dx1,dy1,dzn,jm,im,diu1,diu2,diu3,diu4,diu5,diu6,diu7,diu8,diu9,sm,f,g,h)

          call run_LES_kernel ( &
            data20, data21, &
            im, jm, km, &
            dt, dx1,dy1,dzn, &
            n, nmax &
            ,p,u,v,w,f,g,h,fold,gold,hold &
!            ,po,uo,vo,wo,foldo,goldo,holdo &
            )

! ------------------------------------------------------------------------------------------------
    print *, '9. Comparing up to adam'
    n=6
          call adam(n,nmax,data21,fold,im,jm,km,gold,hold,fghold,f,g,h)

          call run_LES_kernel ( &
            data20, data21, &
            im, jm, km, &
            dt, dx1,dy1,dzn, &
            n, nmax &
            ,p,u,v,w,f,g,h,fold,gold,hold &
!            ,po,uo,vo,wo,foldo,goldo,holdo &
            )

! ------------------------------------------------------------------------------------------------
    print *, '10. Comparing press in stages'
    n=7
        call press(km,jm,im,rhs,u,dx1,v,dy1,w,dzn,f,g,h,dt,cn1,cn2l,p,cn2s,cn3l,cn3s,cn4l,cn4s,n, &
      nmax,data20,usum,vsum,wsum)

          call run_LES_kernel ( &
            data20, data21, &
            im, jm, km, &
            dt, dx1,dy1,dzn, &
            n, nmax &
            ,p,u,v,w,f,g,h,fold,gold,hold &
!            ,po,uo,vo,wo,foldo,goldo,holdo &
            )
      end program

