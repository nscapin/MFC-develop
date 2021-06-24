module m_weno

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures

    use openacc

    use m_mpi_proxy
    ! ==========================================================================

    ! implicit none

    private; public :: s_initialize_weno_module, &
                       s_finalize_weno_module, s_weno_alt
                    ! s_weno, &

    type(vector_field), allocatable, dimension(:) :: v_rs_wsL

    type(scalar_field), allocatable, dimension(:) :: v_vf
    type(scalar_field), allocatable, dimension(:) :: vL_vf, vR_vf

    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_L
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_R
    !$acc declare create (poly_coef_L,poly_coef_R)

    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_L
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_R
    !$acc declare create (d_L,d_R)

    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: beta_coef
    !$acc declare create (beta_coef)

contains

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_weno_module() ! --------------------------------

        type(bounds_info) :: ix, iy, iz !< Indical bounds in the x-, y- and z-directions
        integer :: i
        ! integer :: ixb_full, ixe_full
        ix%beg = -buff_size + weno_polyn; ix%end = m - ix%beg

        ! Allocating WENO-stencil for the variables to be WENO-reconstructed
        allocate (v_rs_wsL(-weno_polyn:weno_polyn))

        do i = -weno_polyn, weno_polyn
            allocate (v_rs_wsL(i)%vf(1:sys_size) )
        end do

        ! Populate variable buffers at each point (for full stencil)
        do i = -weno_polyn, weno_polyn
            do j = 1, sys_size
                allocate (v_rs_wsL(i)%vf(j)%sf(ix%beg:ix%end, &
                                               iy%beg:iy%end, &
                                              iz%beg:iz%end))
           end do
        end do

        ! Allocating/Computing WENO Coefficients in x-direction ============
        allocate (poly_coef_L(0:weno_polyn, &
                              0:weno_polyn - 1, &
                              ix%beg:ix%end))
        allocate (poly_coef_R(0:weno_polyn, &
                              0:weno_polyn - 1, &
                              ix%beg:ix%end))

        allocate (d_L(0:weno_polyn, ix%beg:ix%end))
        allocate (d_R(0:weno_polyn, ix%beg:ix%end))

        allocate (beta_coef(0:weno_polyn, &
                            0:2*(weno_polyn - 1), &
                            ix%beg:ix%end))

        call s_compute_weno_coefficients(ix)

    end subroutine s_initialize_weno_module ! ------------------------------



    subroutine s_weno_alt(qK_prim_vf_flat, vL_vf_flat, vR_vf_flat, ix, iy, iz)

        !! This -4 is because you can't pass arrays with negative indexing consistently
        !! and it only applys to WENO5 (buffsize = 4), also the other directions start/end at 0
        real(kind(0d0)), dimension(-4:,0:,0:,1:), intent(INOUT) :: qK_prim_vf_flat, &
            vL_vf_flat, vR_vf_flat

        type(bounds_info), intent(IN) :: ix, iy, iz

        real(kind(0d0)), dimension(-2:1) :: dvd
        real(kind(0d0)), dimension(0:2) ::  poly
        real(kind(0d0)), dimension(0:2) :: alpha
        real(kind(0d0)), dimension(0:2) :: omega
        real(kind(0d0)), dimension(0:2) :: beta
        real(kind(0d0)), pointer :: beta_p(:)

        integer :: i, j, k, l, r, s
        integer :: ixb, ixe
        
        integer :: t1, t2, c_rate, c_max

        ! For MP_WENO
        real(kind(0d0)), dimension(-1:1) :: d
        real(kind(0d0)) :: d_MD, d_LC
        real(kind(0d0)) :: vL_UL, vR_UL
        real(kind(0d0)) :: vL_MD, vR_MD
        real(kind(0d0)) :: vL_LC, vR_LC
        real(kind(0d0)) :: vL_min, vR_min
        real(kind(0d0)) :: vL_max, vR_max
        real(kind(0d0)), parameter :: alpha_mp = 2d0
        real(kind(0d0)), parameter :: beta_mp  = 4d0/3d0

        ixb = ix%beg
        ixe = ix%end

        k = 0; l = 0

        !$acc data present(qK_prim_vf_flat, vL_vf_flat, vR_vf_flat, poly_coef_L, poly_coef_R, D_L, D_R, beta_coef, dx)
        !$acc parallel loop collapse (2) gang vector private(dvd, poly, beta, alpha, omega)
        do j = ixb, ixe
            do i = 1, sys_size

                !!! L Reconstruction
                dvd(1)  = qK_prim_vf_flat(j + 2, k, l, i) &
                        - qK_prim_vf_flat(j + 1, k, l, i)
                dvd(0)  = qK_prim_vf_flat(j + 1, k, l, i) &
                        - qK_prim_vf_flat(j, k, l, i)
                dvd(-1) = qK_prim_vf_flat(j, k, l, i) &
                        - qK_prim_vf_flat(j - 1, k, l, i)
                dvd(-2) = qK_prim_vf_flat(j - 1, k, l, i) &
                        - qK_prim_vf_flat(j - 2, k, l, i)

                poly(0) = qK_prim_vf_flat(j, k, l, i) &
                        + poly_coef_L(0, 0, j)*dvd(1) &
                        + poly_coef_L(0, 1, j)*dvd(0)
                poly(1) = qK_prim_vf_flat(j, k, l, i) &
                        + poly_coef_L(1, 0, j)*dvd(0) &
                        + poly_coef_L(1, 1, j)*dvd(-1)
                poly(2) = qK_prim_vf_flat(j, k, l, i) &
                        + poly_coef_L(2, 0, j)*dvd(-1) &
                        + poly_coef_L(2, 1, j)*dvd(-2)

                beta(0) = beta_coef(0, 0, j)*dvd(1)*dvd(1) &
                        + beta_coef(0, 1, j)*dvd(1)*dvd(0) &
                        + beta_coef(0, 2, j)*dvd(0)*dvd(0) &
                        + weno_eps
                beta(1) = beta_coef(1, 0, j)*dvd(0)*dvd(0) &
                        + beta_coef(1, 1, j)*dvd(0)*dvd(-1) &
                        + beta_coef(1, 2, j)*dvd(-1)*dvd(-1) &
                        + weno_eps
                beta(2) = beta_coef(2, 0, j)*dvd(-1)*dvd(-1) &
                        + beta_coef(2, 1, j)*dvd(-1)*dvd(-2) &
                        + beta_coef(2, 2, j)*dvd(-2)*dvd(-2) &
                        + weno_eps

                alpha = d_L(:, j)/(beta*beta)
                omega = alpha/sum(alpha)

                if (mapped_weno) then
                    call s_map_nonlinear_weights(d_L(:, j), alpha, omega)
                end if

                vL_vf_flat(j, k, l, i) = sum(omega*poly)

                !!! R Reconstruction
                dvd(1)  = qK_prim_vf_flat(j + 2, k, l, i) &
                        - qK_prim_vf_flat(j + 1, k, l, i)
                dvd(0)  = qK_prim_vf_flat(j + 1, k, l, i) &
                        - qK_prim_vf_flat(j, k, l, i)
                dvd(-1) = qK_prim_vf_flat(j, k, l, i) &
                        - qK_prim_vf_flat(j - 1, k, l, i)
                dvd(-2) = qK_prim_vf_flat(j - 1, k, l, i) &
                        - qK_prim_vf_flat(j - 2, k, l, i)

                poly(0) = qK_prim_vf_flat(j, k, l, i) &
                        + poly_coef_R(0, 0, j)*dvd(1) &
                        + poly_coef_R(0, 1, j)*dvd(0)
                poly(1) = qK_prim_vf_flat(j, k, l, i) &
                        + poly_coef_R(1, 0, j)*dvd(0) &
                        + poly_coef_R(1, 1, j)*dvd(-1)
                poly(2) = qK_prim_vf_flat(j, k, l, i) &
                        + poly_coef_R(2, 0, j)*dvd(-1) &
                        + poly_coef_R(2, 1, j)*dvd(-2)

                beta(0) = beta_coef(0, 0, j)*dvd(1)*dvd(1) &
                        + beta_coef(0, 1, j)*dvd(1)*dvd(0) &
                        + beta_coef(0, 2, j)*dvd(0)*dvd(0) &
                        + weno_eps
                beta(1) = beta_coef(1, 0, j)*dvd(0)*dvd(0) &
                        + beta_coef(1, 1, j)*dvd(0)*dvd(-1) &
                        + beta_coef(1, 2, j)*dvd(-1)*dvd(-1) &
                        + weno_eps
                beta(2) = beta_coef(2, 0, j)*dvd(-1)*dvd(-1) &
                        + beta_coef(2, 1, j)*dvd(-1)*dvd(-2) &
                        + beta_coef(2, 2, j)*dvd(-2)*dvd(-2) &
                        + weno_eps

                alpha = d_R(:, j)/(beta*beta)
                omega = alpha/sum(alpha)

                if (mapped_weno) then
                    call s_map_nonlinear_weights(d_R(:, j), alpha, omega)
                end if

                vR_vf_flat(j, k, l, i) = sum(omega*poly)
            end do
        end do
        !$acc end parallel loop 

        if (mp_weno) then
        !$acc parallel loop collapse (2) gang vector private(d)
        do j = ixb, ixe
            do i = 1, sys_size

                ! Left Monotonicity Preserving Bound
                d(-1) = qK_prim_vf_flat(j, k, l, i) &
                      + qK_prim_vf_flat(j - 2, k, l, i) &
                      - qK_prim_vf_flat(j - 1, k, l, i) &
                      * 2d0
                d(0)  = qK_prim_vf_flat(j + 1, k, l, i) &
                      + qK_prim_vf_flat(j - 1, k, l, i) &
                      - qK_prim_vf_flat(j, k, l, i) &
                      * 2d0
                d(1)  = qK_prim_vf_flat(j + 2, k, l, i) &
                      + qK_prim_vf_flat(j, k, l, i) &
                      - qK_prim_vf_flat(j + 1, k, l, i) &
                      * 2d0

                d_MD = (sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, 4d0*d(0) - d(-1))) &
                       *abs((sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, d(-1))) &
                            *(sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, d(0)))) &
                       *min(abs(4d0*d(-1) - d(0)), abs(d(-1)), &
                            abs(4d0*d(0) - d(-1)), abs(d(0)))/8d0

                d_LC = (sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, 4d0*d(1) - d(0))) &
                       *abs((sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, d(0))) &
                            *(sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, d(1)))) &
                       *min(abs(4d0*d(0) - d(1)), abs(d(0)), &
                            abs(4d0*d(1) - d(0)), abs(d(1)))/8d0

                vL_UL = qK_prim_vf_flat(j, k, l, i) &
                        - (qK_prim_vf_flat(j + 1, k, l, i) &
                           - qK_prim_vf_flat(j, k, l, i))*alpha_mp

                vL_MD = (qK_prim_vf_flat(j, k, l, i) &
                         + qK_prim_vf_flat(j - 1, k, l, i) &
                         - d_MD)*5d-1

                vL_LC = qK_prim_vf_flat(j, k, l, i) &
                      - (qK_prim_vf_flat(j + 1, k, l, i) &
                      -  qK_prim_vf_flat(j, k, l, i))*5d-1 + beta_mp*d_LC

                vL_min = max(min(qK_prim_vf_flat(j, k, l, i), &
                                 qK_prim_vf_flat(j - 1, k, l, i), &
                                 vL_MD), &
                             min(qK_prim_vf_flat(j, k, l, i), &
                                 vL_UL, &
                                 vL_LC))

                vL_max = min(max(qK_prim_vf_flat(j, k, l, i), &
                                 qK_prim_vf_flat(j - 1, k, l, i), &
                                 vL_MD), &
                             max(qK_prim_vf_flat(j, k, l, i), &
                                 vL_UL, &
                                 vL_LC))

                vL_vf_flat(j, k, l, i) = vL_vf_flat(j, k, l, i) &
                          + (sign(5d-1, vL_min - vL_vf_flat(j, k, l, i)) &
                             + sign(5d-1, vL_max - vL_vf_flat(j, k, l, i))) &
                          *min(abs(vL_min - vL_vf_flat(j, k, l, i)), &
                               abs(vL_max - vL_vf_flat(j, k, l, i)))

                ! Right Monotonicity Preserving Bound
                d(-1) = qK_prim_vf_flat(j, k, l, i) &
                      + qK_prim_vf_flat(j - 2, k, l, i) &
                      - qK_prim_vf_flat(j - 1, k, l, i)*2d0
                d(0)  = qK_prim_vf_flat(j + 1, k, l, i) &
                      + qK_prim_vf_flat(j - 1, k, l, i) &
                      - qK_prim_vf_flat(j, k, l, i)*2d0
                d(1)  = qK_prim_vf_flat(j + 2, k, l, i) &
                      + qK_prim_vf_flat(j, k, l, i) &
                      - qK_prim_vf_flat(j + 1, k, l, i)*2d0

                d_MD = (sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, 4d0*d(1) - d(0))) &
                       *abs((sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, d(0))) &
                            *(sign(1d0, 4d0*d(0) - d(1)) + sign(1d0, d(1)))) &
                       *min(abs(4d0*d(0) - d(1)), abs(d(0)), &
                            abs(4d0*d(1) - d(0)), abs(d(1)))/8d0

                d_LC = (sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, 4d0*d(0) - d(-1))) &
                       *abs((sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, d(-1))) &
                            *(sign(1d0, 4d0*d(-1) - d(0)) + sign(1d0, d(0)))) &
                       *min(abs(4d0*d(-1) - d(0)), abs(d(-1)), &
                            abs(4d0*d(0) - d(-1)), abs(d(0)))/8d0

                vR_UL = qK_prim_vf_flat(j, k, l, i) &
                      + (qK_prim_vf_flat(j, k, l, i) &
                      - qK_prim_vf_flat(j - 1, k, l, i))*alpha_mp

                vR_MD = (qK_prim_vf_flat(j, k, l, i) &
                      + qK_prim_vf_flat(j + 1, k, l, i) &
                      - d_MD)*5d-1

                vR_LC = qK_prim_vf_flat(j, k, l, i) &
                      + (qK_prim_vf_flat(j, k, l, i) &
                      - qK_prim_vf_flat(j - 1, k, l, i))*5d-1 + beta_mp*d_LC

                vR_min = max(min(qK_prim_vf_flat(j, k, l, i), &
                                 qK_prim_vf_flat(j + 1, k, l, i), &
                                 vR_MD), &
                             min(qK_prim_vf_flat(j, k, l, i), &
                                 vR_UL, &
                                 vR_LC))

                vR_max = min(max(qK_prim_vf_flat(j, k, l, i), &
                                 qK_prim_vf_flat(j + 1, k, l, i), &
                                 vR_MD), &
                             max(qK_prim_vf_flat(j, k, l, i), &
                                 vR_UL, &
                                 vR_LC))

                vR_vf_flat(j, k, l, i) = vR_vf_flat(j, k, l, i) &
                          + (sign(5d-1, vR_min - vR_vf_flat(j, k, l, i)) &
                             + sign(5d-1, vR_max - vR_vf_flat(j, k, l, i))) &
                          * min(abs(vR_min - vR_vf_flat(j, k, l, i)), &
                               abs(vR_max - vR_vf_flat(j, k, l, i)))
            end do
        end do
        !$acc end parallel loop
        end if
        !$acc end data


        ! do i = ixb, ixe
        !     print*, 'v,vL,vR', &
        !         qK_prim_vf_flat(i,0,0,1), &
        !         vL_vf_flat(i,0,0,1), &
        !         vR_vf_flat(i,0,0,1)
        ! end do

        ! do j = 1,sys_size
        !     do i = ixb,ixe
        !         vL_vf(j)%sf(i,0,0) = vL_vf_flat(i,0,0,j)
        !         vR_vf(j)%sf(i,0,0) = vR_vf_flat(i,0,0,j)
        !     end do
        ! end do

        ! do i = ixb, ixe
        !     print*, '** vL, vR ', &
        !         vL_vf(1)%sf(i,0,0), &
        !         vR_vf(1)%sf(i,0,0)
        ! end do

        ! call s_mpi_abort()

    end subroutine s_weno_alt

    subroutine s_map_nonlinear_weights(d_K, alpha_K, omega_K) 
    !$acc routine seq 

        real(kind(0d0)), dimension(0:2), intent(IN)    ::     d_K
        real(kind(0d0)), dimension(0:2), intent(INOUT) :: alpha_K
        real(kind(0d0)), dimension(0:2), intent(INOUT) :: omega_K

        ! Mapping the WENO nonlinear weights to the WENOM nonlinear weights
        if (minval(d_K) == 0d0 .or. maxval(d_K) == 1d0) return

        alpha_K = (d_K*(1d0 + d_K - 3d0*omega_K) + omega_K**2d0) &
                  *(omega_K/(d_K**2d0 + omega_K*(1d0 - 2d0*d_K)))

        omega_K = alpha_K/sum(alpha_K)

    end subroutine s_map_nonlinear_weights ! -------------------------------



    subroutine s_compute_weno_coefficients(is) ! -------

        type(bounds_info), intent(IN) :: is
        real(kind(0d0)), pointer, dimension(:) :: s_cb => null()
        type(bounds_info) :: bc_s
        integer :: i

        s_cb => x_cb; bc_s = bc_x

        do i = is%beg - 1, is%end - 1

            poly_coef_R(0, 0, i + 1) = &
                ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i + 2)))/ &
                ((s_cb(i) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i + 1)))
            poly_coef_R(1, 0, i + 1) = &
                ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i)))/ &
                ((s_cb(i - 1) - s_cb(i + 2))*(s_cb(i + 2) - s_cb(i)))
            poly_coef_R(1, 1, i + 1) = &
                ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i + 2)))/ &
                ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i - 1) - s_cb(i + 2)))
            poly_coef_R(2, 1, i + 1) = &
                ((s_cb(i) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))/ &
                ((s_cb(i - 2) - s_cb(i))*(s_cb(i - 2) - s_cb(i + 1)))
            poly_coef_L(0, 0, i + 1) = &
                ((s_cb(i + 1) - s_cb(i))*(s_cb(i) - s_cb(i + 2)))/ &
                ((s_cb(i) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i + 1)))
            poly_coef_L(1, 0, i + 1) = &
                ((s_cb(i) - s_cb(i - 1))*(s_cb(i) - s_cb(i + 1)))/ &
                ((s_cb(i - 1) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 2)))
            poly_coef_L(1, 1, i + 1) = &
                ((s_cb(i + 1) - s_cb(i))*(s_cb(i) - s_cb(i + 2)))/ &
                ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i - 1) - s_cb(i + 2)))
            poly_coef_L(2, 1, i + 1) = &
                ((s_cb(i - 1) - s_cb(i))*(s_cb(i) - s_cb(i + 1)))/ &
                ((s_cb(i - 2) - s_cb(i))*(s_cb(i - 2) - s_cb(i + 1)))

            poly_coef_R(0, 1, i + 1) = &
                ((s_cb(i) - s_cb(i + 2)) + (s_cb(i + 1) - s_cb(i + 3)))/ &
                ((s_cb(i) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 3)))* &
                ((s_cb(i) - s_cb(i + 1)))
            poly_coef_R(2, 0, i + 1) = &
                ((s_cb(i - 2) - s_cb(i + 1)) + (s_cb(i - 1) - s_cb(i + 1)))/ &
                ((s_cb(i - 1) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 2)))* &
                ((s_cb(i + 1) - s_cb(i)))
            poly_coef_L(0, 1, i + 1) = &
                ((s_cb(i) - s_cb(i + 2)) + (s_cb(i) - s_cb(i + 3)))/ &
                ((s_cb(i) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 3)))* &
                ((s_cb(i + 1) - s_cb(i)))
            poly_coef_L(2, 0, i + 1) = &
                ((s_cb(i - 2) - s_cb(i)) + (s_cb(i - 1) - s_cb(i + 1)))/ &
                ((s_cb(i - 2) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))* &
                ((s_cb(i) - s_cb(i + 1)))

            d_R(0, i + 1) = &
                ((s_cb(i - 2) - s_cb(i + 1))*(s_cb(i + 1) - s_cb(i - 1)))/ &
                ((s_cb(i - 2) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i - 1)))
            d_R(2, i + 1) = &
                ((s_cb(i + 1) - s_cb(i + 2))*(s_cb(i + 1) - s_cb(i + 3)))/ &
                ((s_cb(i - 2) - s_cb(i + 2))*(s_cb(i - 2) - s_cb(i + 3)))
            d_L(0, i + 1) = &
                ((s_cb(i - 2) - s_cb(i))*(s_cb(i) - s_cb(i - 1)))/ &
                ((s_cb(i - 2) - s_cb(i + 3))*(s_cb(i + 3) - s_cb(i - 1)))
            d_L(2, i + 1) = &
                ((s_cb(i) - s_cb(i + 2))*(s_cb(i) - s_cb(i + 3)))/ &
                ((s_cb(i - 2) - s_cb(i + 2))*(s_cb(i - 2) - s_cb(i + 3)))

            d_R(1, i + 1) = 1d0 - d_R(0, i + 1) - d_R(2, i + 1)
            d_L(1, i + 1) = 1d0 - d_L(0, i + 1) - d_L(2, i + 1)

            beta_coef(0, 0, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + (s_cb(i + 1) - s_cb(i))*(s_cb(i + 2) - &
                s_cb(i + 1)) + (s_cb(i + 2) - s_cb(i + 1))**2d0)/((s_cb(i) - &
                s_cb(i + 3))**2d0*(s_cb(i + 1) - s_cb(i + 3))**2d0)

            beta_coef(0, 1, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(19d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 - (s_cb(i + 1) - s_cb(i))*(s_cb(i + 3) - &
                s_cb(i + 1)) + 2d0*(s_cb(i + 2) - s_cb(i))*((s_cb(i + 2) - &
                s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1))))/((s_cb(i) - &
                s_cb(i + 2))*(s_cb(i) - s_cb(i + 3))**2d0*(s_cb(i + 3) - &
                s_cb(i + 1)))

            beta_coef(0, 2, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + (s_cb(i + 1) - s_cb(i))*((s_cb(i + 2) - &
                s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1))) + ((s_cb(i + 2) - &
                s_cb(i)) + (s_cb(i + 3) - s_cb(i + 1)))**2d0)/((s_cb(i) - &
                s_cb(i + 2))**2d0*(s_cb(i) - s_cb(i + 3))**2d0)

            beta_coef(1, 0, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + (s_cb(i) - s_cb(i - 1))**2d0 + (s_cb(i) - &
                s_cb(i - 1))*(s_cb(i + 1) - s_cb(i)))/((s_cb(i - 1) - &
                s_cb(i + 2))**2d0*(s_cb(i) - s_cb(i + 2))**2d0)

            beta_coef(1, 1, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*((s_cb(i) - &
                s_cb(i + 1))*((s_cb(i) - s_cb(i - 1)) + 20d0*(s_cb(i + 1) - &
                s_cb(i))) + (2d0*(s_cb(i) - s_cb(i - 1)) + (s_cb(i + 1) - &
                s_cb(i)))*(s_cb(i + 2) - s_cb(i)))/((s_cb(i + 1) - &
                s_cb(i - 1))*(s_cb(i - 1) - s_cb(i + 2))**2d0*(s_cb(i + 2) - &
                s_cb(i)))

            beta_coef(1, 2, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + (s_cb(i + 1) - s_cb(i))*(s_cb(i + 2) - &
                s_cb(i + 1)) + (s_cb(i + 2) - s_cb(i + 1))**2d0)/ &
                ((s_cb(i - 1) - s_cb(i + 1))**2d0*(s_cb(i - 1) - &
                                                   s_cb(i + 2))**2d0)

            beta_coef(2, 0, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(12d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + ((s_cb(i) - s_cb(i - 2)) + (s_cb(i) - &
                s_cb(i - 1)))**2d0 + 3d0*((s_cb(i) - s_cb(i - 2)) + &
                (s_cb(i) - s_cb(i - 1)))*(s_cb(i + 1) - s_cb(i)))/ &
                ((s_cb(i - 2) - s_cb(i + 1))**2d0*(s_cb(i - 1) - &
                                                   s_cb(i + 1))**2d0)

            beta_coef(2, 1, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(19d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + ((s_cb(i) - s_cb(i - 2))*(s_cb(i) - &
                s_cb(i + 1))) + 2d0*(s_cb(i + 1) - s_cb(i - 1))*((s_cb(i) - &
                s_cb(i - 2)) + (s_cb(i + 1) - s_cb(i - 1))))/((s_cb(i - 2) - &
                s_cb(i))*(s_cb(i - 2) - s_cb(i + 1))**2d0*(s_cb(i + 1) - &
                s_cb(i - 1)))

            beta_coef(2, 2, i + 1) = &
                4d0*(s_cb(i) - s_cb(i + 1))**2d0*(10d0*(s_cb(i + 1) - &
                s_cb(i))**2d0 + (s_cb(i) - s_cb(i - 1))**2d0 + (s_cb(i) - &
                s_cb(i - 1))*(s_cb(i + 1) - s_cb(i)))/((s_cb(i - 2) - &
                s_cb(i))**2d0*(s_cb(i - 2) - s_cb(i + 1))**2d0)

        end do

        !$acc update device( beta_coef, D_L, D_R, poly_coef_L, poly_coef_R )

        nullify(s_cb)

    end subroutine s_compute_weno_coefficients ! ---------------------------

    subroutine s_finalize_weno_module() ! ----------------------------------

        deallocate (v_rs_wsL)
        deallocate (poly_coef_L, poly_coef_R)
        deallocate (d_L, d_R)
        deallocate (beta_coef)
        ! deallocate(v_rs_wsL_flat)

    end subroutine s_finalize_weno_module ! --------------------------------

end module m_weno
