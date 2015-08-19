module module_LES_combined_kernel
contains
      subroutine LES_combined_kernel( p_scratch, uvw, uvwsum, fgh, fgh_old, &
      rhs, mask1, diu, sm, &
        dxs, dys, dzs, dx1, dy1, dzn, &
        z2, &
        cn1, cn2l, cn2s, cn3l, cn3s, cn4l, cn4s, &
        val_ptr, chunks_num, chunks_denom, n_ptr, state_ptr, dt, im, jm, km &
        , p_halo, uvw_halo, uvwsum_halo, fgh_halo, fgh_old_halo, diu_halo, rhs_halo, sm_halo &
        )
      use common_sn
        real(kind=4), dimension(00:1,:ip+2,0:jp+2,0:kp+1), intent(In) :: p_scratch
        real(kind=4), dimension(0:3,0:ip+1,-1:jp+1,-1:kp+1), intent(InOut) :: uvw
        real(kind=4), dimension(0:3,0:ip,0:jp,0:kp), intent(InOut) :: uvwsum
        real(kind=4), dimension(0:3,0:ip,0:jp,0:kp), intent(InOut) :: fgh
        real(kind=4), dimension(0:3,ip,jp,kp), intent(InOut) :: fgh_old
        real(kind=4), dimension(0:ip+1,0:jp+1,0:kp+1), intent(In) :: rhs
        real(kind=4), dimension(0:3,-1:ip+1,-1:jp+1,0:kp+1) :: mask1
        real(kind=4), dimension(1:16,-1:ip+2,0:jp+2,0:kp+2) :: diu
        real(kind=4), dimension(-1:ip+1,-1:jp+1,0:kp+1), intent(In) :: sm
        real(kind=4), dimension(0:ip), intent(In) :: dxs
        real(kind=4), dimension(0:jp), intent(In) :: dys
        real(kind=4), dimension(-1:kp+2), intent(In) :: dzs
        real(kind=4), dimension(-1:ip+1), intent(In) :: dx1
        real(kind=4), dimension(0:jp+1), intent(In) :: dy1
        real(kind=4), dimension(-1:kp+2), intent(In) :: dzn
        real(kind=4), dimension(kp+2), intent(In) :: z2
        real(kind=4), dimension(ip,jp,kp), intent(In) :: cn1
        real(kind=4), dimension(ip), intent(In) :: cn2l
        real(kind=4), dimension(ip), intent(In) :: cn2s
        real(kind=4), dimension(jp), intent(In) :: cn3l
        real(kind=4), dimension(jp), intent(In) :: cn3s
        real(kind=4), dimension(kp), intent(In) :: cn4l
        real(kind=4), dimension(kp), intent(In) :: cn4s
        real(kind=4), dimension(256), intent(InOut) :: val_ptr
        real(kind=4), dimension(kp), intent(InOut) :: chunks_num, chunks_denom
        integer, dimension(256), intent(InOut) :: n_ptr, state_ptr
        real(kind=4), intent(In) :: dt
        integer, intent(In) :: im
        integer, intent(In) :: jm
        integer, intent(In) :: km
        ! Halos
        ! Putting the size expressions directly in the arrays crash the combined script
        ! Probably the * symbol
        integer, parameter :: s_p = 4 * (ip+jp+4) * (kp+2)
        integer, parameter :: s_uvw = 8 * (ip+jp+3) * (kp+3)
        integer, parameter :: s_uvwsum = 8 * (ip+jp) * (kp+1)
        integer, parameter :: s_fgh = 8 * (ip+jp) * (kp+1)
        integer, parameter :: s_fgh_old = 8 * (ip+jp-2) * kp
        integer, parameter :: s_diu = 32 * (ip+jp+5) * (kp+3)
        integer, parameter :: s_rhs = 2 * (ip+jp+2) * (kp+2)
        integer, parameter :: s_sm = 2 * (ip+jp+4) * (kp+2)
        real(kind=4), dimension(s_p) :: p_halo
        real(kind=4), dimension(s_uvw) :: uvw_halo
        real(kind=4), dimension(s_uvwsum) :: uvwsum_halo
        real(kind=4), dimension(s_fgh) :: fgh_halo
        real(kind=4), dimension(s_fgh_old) :: fgh_old_halo
        real(kind=4), dimension(s_diu) :: diu_halo
        real(kind=4), dimension(s_rhs) :: rhs_halo
        real(kind=4), dimension(s_sm) :: sm_halo
!
end module module_LES_combined_kernel
