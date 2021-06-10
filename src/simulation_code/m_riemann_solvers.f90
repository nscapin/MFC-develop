! typical conditions are avg state -> 2, wave_speeds->1

module m_riemann_solvers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_bubbles              !< To get the bubble wall pressure function
    ! ==========================================================================

    implicit none

    private; public :: &
        s_initialize_riemann_solvers_module, &
        s_hllc_riemann_solver, &
        s_finalize_riemann_solvers_module

    type(scalar_field), allocatable, dimension(:) :: qL_prim_rs_vf
    type(scalar_field), allocatable, dimension(:) :: qR_prim_rs_vf
    type(scalar_field), allocatable, dimension(:) :: q_prim_rs_vf
    type(scalar_field), allocatable, dimension(:) :: flux_rs_vf, flux_src_rs_vf
    type(scalar_field), allocatable, dimension(:) :: flux_gsrc_rs_vf
    type(scalar_field), allocatable, dimension(:) :: vel_src_rs_vf

    real(kind(0d0)), allocatable, dimension(:)   :: alpha_rho_L, alpha_rho_R
    real(kind(0d0))                              ::       rho_L, rho_R
    real(kind(0d0)), allocatable, dimension(:)   ::       vel_L, vel_R
    real(kind(0d0))                              ::      pres_L, pres_R
    real(kind(0d0))                              ::         E_L, E_R
    real(kind(0d0))                              ::         H_L, H_R
    real(kind(0d0)), allocatable, dimension(:)   ::     alpha_L, alpha_R
    real(kind(0d0))                              ::     gamma_L, gamma_R
    real(kind(0d0))                              ::    pi_inf_L, pi_inf_R
    real(kind(0d0))                              ::         c_L, c_R

    real(kind(0d0))                                 :: rho_avg
    real(kind(0d0)), allocatable, dimension(:)      :: vel_avg
    real(kind(0d0))                                 :: H_avg
    real(kind(0d0))                                 :: gamma_avg
    real(kind(0d0))                                 :: c_avg
    real(kind(0d0)) :: s_L, s_R, s_S
    real(kind(0d0)) :: s_M, s_P
    real(kind(0d0)) :: xi_M, xi_P

    type(bounds_info) :: is1, is2, is3

contains


    !!  @param qL_prim_vf The left WENO-reconstructed cell-boundary values of the
    !!      cell-average primitive variables
    !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
    !!      cell-average primitive variables
    !!  @param dqL_prim_dx_vf The left WENO-reconstructed cell-boundary values of the
    !!      first-order x-dir spatial derivatives
    !!  @param dqR_prim_dx_vf The right WENO-reconstructed cell-boundary values of the
    !!      first-order x-dir spatial derivatives
    !!  @param gm_alphaL_vf Left averaged gradient magnitude
    !!  @param gm_alphaR_vf Right averaged gradient magnitude
    !!  @param flux_vf Intra-cell fluxes
    !!  @param flux_src_vf Intra-cell fluxes sources
    !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
    !!  @param norm_dir Dir. splitting direction
    !!  @param q_prim_vf Cell-averaged primitive variables
    subroutine s_hllc_riemann_solver(qL_prim_vf, dqL_prim_dx_vf, & 
                                     gm_alphaL_vf, &
                                     qR_prim_vf, dqR_prim_dx_vf, &
                                     gm_alphaR_vf, &
                                     q_prim_vf, &
                                     flux_vf, flux_src_vf, &
                                     flux_gsrc_vf, &
                                     norm_dir, ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: qL_prim_vf, qR_prim_vf
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             gm_alphaL_vf, gm_alphaR_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(IN) :: norm_dir
        type(bounds_info), intent(IN) :: ix, iy, iz

        real(kind(0d0)) :: xi_L, xi_R

        integer :: i, j, k, l

        ! Reshaping input data based on dimensional splitting direction
        call s_initialize_riemann_solver(qL_prim_vf, &
                                         qR_prim_vf, &
                                         q_prim_vf, &
                                         flux_vf, flux_src_vf, &
                                         flux_gsrc_vf, &
                                         norm_dir, ix, iy, iz)

        ! Computing HLLC flux and source flux for Euler system of equations
        do l = is3%beg, is3%end
            do k = is2%beg, is2%end
                do j = is1%beg, is1%end

                    call s_compute_arithmetic_average_state(j, k, l)
                    call s_compute_direct_wave_speeds(j, k, l)


                    s_M = min(0d0, s_L)
                    s_P = max(0d0, s_R)

                    xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                    xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)

                    xi_M = (5d-1 + sign(5d-1, s_S))
                    xi_P = (5d-1 - sign(5d-1, s_S))

                    do i = 1, cont_idx%end
                        flux_rs_vf(i)%sf(j, k, l) = &
                            xi_M*alpha_rho_L(i) &
                            *(vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                            + xi_P*alpha_rho_R(i) &
                            *(vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                    end do

                    ! Momentum flux.
                    do i = 1, num_dims
                        flux_rs_vf(cont_idx%end + dir_idx(i))%sf(j, k, l) = &
                            xi_M*(rho_L*(vel_L(dir_idx(1))* &
                                  vel_L(dir_idx(i)) + &
                                  s_M*(xi_L*(dir_flg(dir_idx(i))*s_S + &
                                  (1d0 - dir_flg(dir_idx(i)))* &
                                  vel_L(dir_idx(i))) - vel_L(dir_idx(i)))) + &
                                  dir_flg(dir_idx(i))*(pres_L)) &
                          + xi_P*(rho_R*(vel_R(dir_idx(1))* &
                                  vel_R(dir_idx(i)) + &
                                  s_P*(xi_R*(dir_flg(dir_idx(i))*s_S + &
                                  (1d0 - dir_flg(dir_idx(i)))* &
                                  vel_R(dir_idx(i))) - vel_R(dir_idx(i)))) + &
                                  dir_flg(dir_idx(i))*(pres_R))
                    end do

                    ! Energy flux
                    flux_rs_vf(E_idx)%sf(j, k, l) = &
                        xi_M*(vel_L(dir_idx(1))*(E_L + pres_L) + &
                             s_M*(xi_L*(E_L + (s_S - vel_L(dir_idx(1)))* &
                             (rho_L*s_S + pres_L/ &
                             (s_L - vel_L(dir_idx(1))))) - E_L)) &
                      + xi_P*(vel_R(dir_idx(1))*(E_R + pres_R) + &
                             s_P*(xi_R*(E_R + (s_S - vel_R(dir_idx(1)))* &
                             (rho_R*s_S + pres_R/ &
                             (s_R - vel_R(dir_idx(1))))) - E_R))

                    ! Volume fraction flux
                    do i = adv_idx%beg, adv_idx%end
                        flux_rs_vf(i)%sf(j, k, l) = &
                            xi_M*qL_prim_rs_vf(i)%sf(j, k, l) &
                            *(vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                            + xi_P*qR_prim_rs_vf(i)%sf(j + 1, k, l) &
                            *(vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                    end do

                    ! Source for volume fraction advection equation
                    do i = 1, num_dims
                        vel_src_rs_vf(dir_idx(i))%sf(j, k, l) = &
                            xi_M*(vel_L(dir_idx(i)) + &
                                  dir_flg(dir_idx(i))* &
                                  s_M*(xi_L - 1d0)) &
                          + xi_P*(vel_R(dir_idx(i)) + &
                                  dir_flg(dir_idx(i))* &
                                  s_P*(xi_R - 1d0))
                    end do

                    do i = 1, sys_size
                        flux_gsrc_rs_vf(i)%sf(j, k, l) = 0d0
                    end do
                end do
            end do
        end do

        ! Reshaping outputted data based on dimensional splitting direction
        call s_finalize_riemann_solver(flux_vf, flux_src_vf, &
                                       flux_gsrc_vf, &
                                       norm_dir, ix, iy, iz)

    end subroutine s_hllc_riemann_solver 


    subroutine s_compute_mixture_sound_speeds(j, k, l) 

        integer, intent(IN) :: j, k, l
        real(kind(0d0)) :: blkmod1, blkmod2
        integer :: i 

        if ((alt_soundspeed .or. regularization)) then
            do i = 1, num_fluids
                alpha_L(i) = qL_prim_rs_vf(E_idx + i)%sf(j, k, l)
                alpha_R(i) = qR_prim_rs_vf(E_idx + i)%sf(j + 1, k, l)
            end do

            blkmod1 = ((fluid_pp(1)%gamma + 1d0)*pres_L + &
                       fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
            blkmod2 = ((fluid_pp(2)%gamma + 1d0)*pres_L + &
                       fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
            c_L = 1d0/(rho_L*(alpha_L(1)/blkmod1 + alpha_L(2)/blkmod2))

            blkmod1 = ((fluid_pp(1)%gamma + 1d0)*pres_R + &
                       fluid_pp(1)%pi_inf)/fluid_pp(1)%gamma
            blkmod2 = ((fluid_pp(2)%gamma + 1d0)*pres_R + &
                       fluid_pp(2)%pi_inf)/fluid_pp(2)%gamma
            c_R = 1d0/(rho_R*(alpha_R(1)/blkmod1 + alpha_R(2)/blkmod2))
        else
            do i = 1, num_fluids
                alpha_L(i) = qL_prim_rs_vf(E_idx + i)%sf(j, k, l)
                alpha_R(i) = qR_prim_rs_vf(E_idx + i)%sf(j + 1, k, l)
            end do
            c_L = ((H_L - 5d-1*sum(vel_L**2d0))/gamma_L)
            c_R = ((H_R - 5d-1*sum(vel_R**2d0))/gamma_R)
        end if

        if (mixture_err .and. c_L < 0d0) then
            c_L = 100.d0*sgm_eps
        else
            c_L = sqrt(c_L)
        end if

        if (mixture_err .and. c_R < 0d0) then
            c_R = 100.d0*sgm_eps
        else
            c_R = sqrt(c_R)
        end if

    end subroutine s_compute_mixture_sound_speeds 


    subroutine s_compute_arithmetic_average_state(j, k, l) 

        integer, intent(IN) :: j, k, l
        integer :: i, q 

        do i = 1, cont_idx%end
            alpha_rho_L(i) = qL_prim_rs_vf(i)%sf(j, k, l)
            alpha_rho_R(i) = qR_prim_rs_vf(i)%sf(j + 1, k, l)
        end do

        do i = 1, num_dims
            vel_L(i) = qL_prim_rs_vf(cont_idx%end + i)%sf(j, k, l)
            vel_R(i) = qR_prim_rs_vf(cont_idx%end + i)%sf(j + 1, k, l)
        end do

        call s_convert_species_to_mixture_variables(qL_prim_rs_vf, &
                                            rho_L, gamma_L, &
                                            pi_inf_L, &
                                            j, k, l)
        call s_convert_species_to_mixture_variables(qR_prim_rs_vf, &
                                            rho_R, gamma_R, &
                                            pi_inf_R,  &
                                            j + 1, k, l)

        pres_L = qL_prim_rs_vf(E_idx)%sf(j, k, l)
        pres_R = qR_prim_rs_vf(E_idx)%sf(j + 1, k, l)

        E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*sum(vel_L**2d0)
        E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*sum(vel_R**2d0)

        H_L = (E_L + pres_L)/rho_L
        H_R = (E_R + pres_R)/rho_R

        call s_compute_mixture_sound_speeds(j, k, l)

        ! Arithmetic Average Riemann Problem State =========================
        rho_avg = 5d-1*(rho_L + rho_R)

        vel_avg = 5d-1*(vel_L + vel_R)

        H_avg = 5d-1*(H_L + H_R)

        gamma_avg = 5d-1*(gamma_L + gamma_R)

        if (mixture_err) then
            if ((H_avg - 5d-1*sum(vel_avg**2d0)) < 0d0) then
                c_avg = sgm_eps
            else
                c_avg = sqrt((H_avg - 5d-1*sum(vel_avg**2d0))/gamma_avg)
            end if
        else
            c_avg = sqrt((H_avg - 5d-1*sum(vel_avg**2d0))/gamma_avg)
        end if

    end subroutine s_compute_arithmetic_average_state ! --------------------


    subroutine s_compute_direct_wave_speeds(j, k, l) ! -----------------------

        integer, intent(IN) :: j, k, l

        real(kind(0d0)) :: denom

        integer :: i !< Generic loop iterator


        s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
        s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)

        s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
               (s_L - vel_L(dir_idx(1))) - &
               rho_R*vel_R(dir_idx(1))* &
               (s_R - vel_R(dir_idx(1)))) &
              /(rho_L*(s_L - vel_L(dir_idx(1))) - &
                rho_R*(s_R - vel_R(dir_idx(1))))
        denom = rho_L*(s_L - vel_L(dir_idx(1))) - rho_R*(s_R - vel_R(dir_idx(1)))


    end subroutine s_compute_direct_wave_speeds ! --------------------------


    subroutine s_initialize_riemann_solvers_module() ! ---------------------

        allocate (qL_prim_rs_vf(1:sys_size), qR_prim_rs_vf(1:sys_size))
        allocate (flux_rs_vf(1:sys_size), flux_src_rs_vf(1:sys_size))
        allocate (flux_gsrc_rs_vf(1:sys_size))
        allocate (vel_src_rs_vf(1:num_dims))
        allocate (alpha_rho_L(1:cont_idx%end), vel_L(1:num_dims))
        allocate (alpha_rho_R(1:cont_idx%end), vel_R(1:num_dims))
        allocate (vel_avg(1:num_dims))
        allocate (alpha_L(1:num_fluids))
        allocate (alpha_R(1:num_fluids))

    end subroutine s_initialize_riemann_solvers_module ! -------------------



    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures needed to configure the chosen Riemann
        !!      solver algorithm.
        !!  @param qL_prim_vf The  left WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param qR_prim_vf The right WENO-reconstructed cell-boundary values of the
        !!      cell-average primitive variables
        !!  @param flux_vf Intra-cell fluxes
        !!  @param flux_src_vf Intra-cell fluxes sources
        !!  @param flux_gsrc_vf Intra-cell geometric fluxes sources
        !!  @param norm_dir Dir. splitting direction
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
        !!  @param q_prim_vf Cell-averaged primitive variables
    subroutine s_initialize_riemann_solver(qL_prim_vf, &
                                           qR_prim_vf, &
                                           q_prim_vf, &
                                           flux_vf, flux_src_vf, &
                                           flux_gsrc_vf, &
                                           norm_dir, ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: qL_prim_vf, qR_prim_vf
        type(scalar_field), dimension(sys_size), intent(IN) :: q_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(IN) :: norm_dir

        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: i, j, k ! Generic loop iterators

        integer :: xbeg, xend, ybeg, yend, zbeg, zend
        integer :: s1beg, s1end, s2beg, s2end, s3beg, s3end

        is1 = ix; is2 = iy; is3 = iz
        dir_idx = (/1, 2, 3/); dir_flg = (/1d0, 0d0, 0d0/)

        ! Setting up special bounds for cell-average values
        xbeg = -buff_size; ybeg = 0; zbeg = 0
        if (n > 0) ybeg = -buff_size; if (p > 0) zbeg = -buff_size
        xend = m - xbeg; yend = n - ybeg; zend = p - zbeg

        ! Configuring the coordinate direction indexes
        s1beg = xbeg; s1end = xend; s2beg = ybeg; s2end = yend; s3beg = zbeg; s3end = zend

        ! Allocating Left, Right and Average Riemann Problem States ========
        do i = 1, sys_size
            allocate (qL_prim_rs_vf(i)%sf(is1%beg:is1%end, &
                                          is2%beg:is2%end, &
                                          is3%beg:is3%end))
            allocate (qR_prim_rs_vf(i)%sf(is1%beg + 1:is1%end + 1, &
                                          is2%beg:is2%end, &
                                          is3%beg:is3%end))
        end do
        ! ==================================================================

        ! Allocating Intercell Fluxes and Velocity =========================
        do i = 1, sys_size
            flux_rs_vf(i)%sf => flux_vf(i)%sf
            flux_src_rs_vf(i)%sf => flux_src_vf(i)%sf
            flux_gsrc_rs_vf(i)%sf => flux_gsrc_vf(i)%sf
        end do

        vel_src_rs_vf(dir_idx(1))%sf => flux_src_rs_vf(adv_idx%beg)%sf


        do i = 1, sys_size
            qL_prim_rs_vf(i)%sf = qL_prim_vf(i)%sf(ix%beg:ix%end, &
                                                   iy%beg:iy%end, &
                                                   iz%beg:iz%end)
            qR_prim_rs_vf(i)%sf = qR_prim_vf(i)%sf(ix%beg + 1:ix%end + 1, &
                                                   iy%beg:iy%end, &
                                                   iz%beg:iz%end)
        end do

    end subroutine s_initialize_riemann_solver ! ---------------------------


    subroutine s_finalize_riemann_solver(flux_vf, flux_src_vf, & ! --------
                                         flux_gsrc_vf, &
                                         norm_dir, ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(IN) :: norm_dir

        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: i

        do i = 1, sys_size
            deallocate (qL_prim_rs_vf(i)%sf, qR_prim_rs_vf(i)%sf)
        end do

        do i = 1, sys_size
            flux_rs_vf(i)%sf => null()
            flux_src_rs_vf(i)%sf => null()
            flux_gsrc_rs_vf(i)%sf => null()
        end do
        vel_src_rs_vf(dir_idx(1))%sf => null()

    end subroutine s_finalize_riemann_solver 


    subroutine s_finalize_riemann_solvers_module() 

        deallocate (qL_prim_rs_vf, qR_prim_rs_vf)
        deallocate (flux_rs_vf, flux_src_rs_vf, flux_gsrc_rs_vf)
        deallocate (vel_src_rs_vf)
        deallocate (alpha_rho_L, vel_L)
        deallocate (alpha_rho_R, vel_R)
        deallocate (vel_avg)
        deallocate (alpha_L, alpha_R)

    end subroutine s_finalize_riemann_solvers_module 

end module m_riemann_solvers
