!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /
!!    / /  / / __/ / /___
!!   /_/  /_/_/    \____/
!!
!!  This file is part of MFC.
!!
!!  MFC is the legal property of its developers, whose names
!!  are listed in the copyright file included with this source
!!  distribution.
!!
!!  MFC is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published
!!  by the Free Software Foundation, either version 3 of the license
!!  or any later version.
!!
!!  MFC is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with MFC (LICENSE).
!!  If not, see <http://www.gnu.org/licenses/>.

!>
!! @file m_riemann_solvers.f90
!! @brief Contains module m_riemann_solvers
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module features a database of approximate and exact Riemann
!!              problem solvers for the Navier-Stokes system of equations, which
!!              is supplemented by appropriate advection equations that are used
!!              to capture the material interfaces. The closure of the system is
!!              achieved by the stiffened gas equation of state and any required
!!              mixture relations. Surface tension effects are accounted for and
!!              are modeled by means of a volume force acting across the diffuse
!!              material interface region. The implementation details of viscous
!!              and capillary effects, into the Riemann solvers, may be found in
!!              Perigaud and Saurel (2005). Note that both effects are available
!!              only in the volume fraction model. At this time, the approximate
!!              and exact Riemann solvers that are listed below are available:
!!                  1) Harten-Lax-van Leer (HLL)
!!                  2) Harten-Lax-van Leer-Contact (HLLC)
!!                  3) Exact
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

    abstract interface ! =======================================================

        !>  The abstract interface to the subroutines that are used to calculate
        !!  the Roe and arithmetic average states. For more information refer to:
        !!      1) s_compute_roe_average_state
        !!      2) s_compute_arithmetic_average_state
        !!  @param i First coordinate location index
        !!  @param j Second coordinate location index
        !!  @param k Third coordinate location index
        subroutine s_compute_abstract_average_state(i, j, k)

            integer, intent(IN) :: i, j, k

        end subroutine s_compute_abstract_average_state

        !> The abstract interface to the subroutines that are utilized to compute
        !! the wave speeds of the Riemann problem either directly or by the means
        !! of pressure-velocity estimates. For more information please refer to:
        !!      1) s_compute_direct_wave_speeds
        !!      2) s_compute_pressure_velocity_wave_speeds
        !!  @param i First coordinate location index
        !!  @param j Second coordinate location index
        !!  @param k Third coordinate location index
        subroutine s_compute_abstract_wave_speeds(i, j, k)

            integer, intent(IN) :: i, j, k

        end subroutine s_compute_abstract_wave_speeds

    end interface ! ============================================================

    !> The left (L) and right (R) WENO-reconstructed cell-boundary values of the
    !! cell-average primitive variables that define the left and right states of
    !! the Riemann problem. Variables qK_prim_rs_vf, K = L or R, are obtained by
    !! reshaping (RS) qK_prim_vf in a coordinate direction that is normal to the
    !! cell-boundaries along which the fluxes are to be determined.
    !> @{
    type(scalar_field), allocatable, dimension(:) :: qL_prim_rs_vf
    type(scalar_field), allocatable, dimension(:) :: qR_prim_rs_vf
    type(scalar_field), allocatable, dimension(:) :: q_prim_rs_vf
    !> @}


    !> The cell-boundary values of the fluxes (src - source) that are computed
    !! through the chosen Riemann problem solver, and the direct evaluation of
    !! source terms, by using the left and right states given in qK_prim_rs_vf,
    !! dqK_prim_ds_vf and kappaK_rs_vf, where ds = dx, dy or dz.
    !> @{
    type(scalar_field), allocatable, dimension(:) :: flux_rs_vf, flux_src_rs_vf
    !> @}

    type(scalar_field), allocatable, dimension(:) :: flux_gsrc_rs_vf !<
    !! The cell-boundary values of the geometrical source flux that are computed
    !! through the chosen Riemann problem solver by using the left and right
    !! states given in qK_prim_rs_vf. Currently 2D axisymmetric for inviscid only.

    ! The cell-boundary values of the velocity. vel_src_rs_vf is determined as
    ! part of Riemann problem solution and is used to evaluate the source flux.
    type(scalar_field), allocatable, dimension(:) :: vel_src_rs_vf

    !> @name Left and right, WENO-reconstructed, cell-boundary values of cell-average
    !! partial densities, density, velocity, pressure, internal energy, energy, enthalpy, volume
    !! fractions, mass fractions, the specific heat ratio and liquid stiffness functions, speed
    !! of sound, shear and volume Reynolds numbers and the Weber numbers. These
    !! variables are left and right states of the Riemann problem obtained from
    !! qK_prim_rs_vf and kappaK_rs_vf.
    !> @{
    real(kind(0d0)), allocatable, dimension(:)   :: alpha_rho_L, alpha_rho_R
    real(kind(0d0))                              ::       rho_L, rho_R
    real(kind(0d0)), allocatable, dimension(:)   ::       vel_L, vel_R
    real(kind(0d0))                              ::      pres_L, pres_R
    real(kind(0d0))                              ::         E_L, E_R
    real(kind(0d0))                              ::         H_L, H_R
    real(kind(0d0)), allocatable, dimension(:)   ::     alpha_L, alpha_R
    real(kind(0d0))                              ::         Y_L, Y_R
    real(kind(0d0))                              ::     gamma_L, gamma_R
    real(kind(0d0))                              ::    pi_inf_L, pi_inf_R
    real(kind(0d0))                              ::         c_L, c_R
    real(kind(0d0)), allocatable, dimension(:)   ::     tau_e_L, tau_e_R

    !> @}

    !> @name Left and right, WENO-reconstructed, cell-boundary values of cell-average
    !! bubble density, radius, radial velocity, pressure, wall pressure, and modified
    !! pressure. These variables are left and right states of the Riemann problem obtained from
    !! qK_prim_rs_vf and kappaK_rs_vf.
    !> @{
    real(kind(0d0))                              ::       nbub_L, nbub_R
    real(kind(0d0)), allocatable, dimension(:)   ::         R0_L, R0_R
    real(kind(0d0)), allocatable, dimension(:)   ::         V0_L, V0_R
    real(kind(0d0)), allocatable, dimension(:)   ::         P0_L, P0_R
    real(kind(0d0)), allocatable, dimension(:)   ::        pbw_L, pbw_R
    real(kind(0d0)), allocatable, dimension(:, :) ::       moms_L, moms_R
    real(kind(0d0))                              ::     ptilde_L, ptilde_R
    !> @}

    !> @name Star region pressure and velocity
    !> @{
    real(kind(0d0)) :: pres_S
    real(kind(0d0)) :: vel_S
    !> @}

    !> @name Intercell solution used to calculated intercell flux
    !> @{
    real(kind(0d0)), allocatable, dimension(:)   :: alpha_rho_IC
    real(kind(0d0))                              :: rho_IC
    real(kind(0d0)), allocatable, dimension(:)   :: vel_IC
    real(kind(0d0))                              :: pres_IC
    real(kind(0d0))                              :: E_IC
    real(kind(0d0)), allocatable, dimension(:)   :: alpha_IC
    real(kind(0d0)), allocatable, dimension(:)   :: tau_e_IC
    !> @}

    !> @name Surface tension pressure contribution
    !> @{
    real(kind(0d0)) :: dpres_L, dpres_R
    !> @}

    !> @name Roe or arithmetic average density, velocity, enthalpy, volume fractions,
    !! specific heat ratio function, speed of sound, shear and volume Reynolds
    !! numbers, Weber numbers and curvatures, at the cell-boundaries, computed
    !! from the left and the right states of the Riemann problem
    !> @{
    real(kind(0d0))                                 :: rho_avg
    real(kind(0d0)), allocatable, dimension(:)   :: vel_avg
    real(kind(0d0))                                 :: H_avg
    type(scalar_field), allocatable, dimension(:)   :: alpha_avg_rs_vf
    real(kind(0d0))                                 :: gamma_avg
    real(kind(0d0))                                 :: c_avg
    !> @}

    !> @name Left, right and star (S) region wave speeds
    !> @{
    real(kind(0d0)) :: s_L, s_R, s_S
    !> @}

    !> @name Star region variables (HLLC)
    !> @{
    real(kind(0d0)) :: rho_Star, E_Star, p_Star, p_K_Star
    !> @}

    !> Minus (M) and plus (P) wave speeds
    !> @{
    real(kind(0d0)) :: s_M, s_P
    !> @}

    !> Minus and plus wave speeds functions
    !> @{
    real(kind(0d0)) :: xi_M, xi_P
    !> @}

    procedure(s_abstract_riemann_solver), &
        pointer :: s_riemann_solver => null() !<
    !! Pointer to the procedure that is utilized to calculate either the HLL,
    !! HLLC or exact intercell fluxes, based on the choice of Riemann solver

    procedure(s_compute_abstract_average_state), &
        pointer :: s_compute_average_state => null() !<
    !! Pointer to the subroutine utilized to calculate either the Roe or the
    !! arithmetic average state variables, based on the chosen average state

    procedure(s_compute_abstract_wave_speeds), &
        pointer :: s_compute_wave_speeds => null() !<
    !! Pointer to the subroutine that is utilized to compute the wave speeds of
    !! the Riemann problem either directly or by the means of pressure-velocity
    !! estimates, based on the selected method of estimation of the wave speeds

    !> @name Indical bounds in the s1-, s2- and s3-directions
    !> @{
    type(bounds_info) :: is1, is2, is3
    !> @}

contains


    !> This procedure is the implementation of the Harten, Lax,
        !!      van Leer, and contact (HLLC) approximate Riemann solver,
        !!      see Toro (1999) and Johnsen (2007). The viscous and the
        !!      surface tension effects have been included by modifying
        !!      the exact Riemann solver of Perigaud and Saurel (2005).
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
        !!  @param ix Index bounds in the x-dir
        !!  @param iy Index bounds in the y-dir
        !!  @param iz Index bounds in the z-dir
        !!  @param q_prim_vf Cell-averaged primitive variables
    subroutine s_hllc_riemann_solver(qL_prim_vf, dqL_prim_dx_vf, & ! ------
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

        ! Intercell fluxes
        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(IN) :: norm_dir
        type(bounds_info), intent(IN) :: ix, iy, iz

        real(kind(0d0)) :: xi_L, xi_R !< Left and right wave speeds functions

        integer :: i, j, k, l !< Generic loop iterators

        ! Populating the buffers of the left and right Riemann problem
        ! states variables, based on the choice of boundary conditions
        call s_populate_riemann_states_variables_buffers( &
            qL_prim_vf, dqL_prim_dx_vf, &
            gm_alphaL_vf, &
            qR_prim_vf, dqR_prim_dx_vf, &
            gm_alphaR_vf, &
            norm_dir, ix, iy, iz)

        ! Reshaping inputted data based on dimensional splitting direction
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

                    call s_compute_average_state(j, k, l)
                    call s_compute_wave_speeds(j, k, l)

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

    end subroutine s_hllc_riemann_solver ! ---------------------------------



    !> Compute mixture sound speed
        !! @param j  First coordinate index
        !! @param k Second coordinate index
        !! @param l  Third coordinate index
    subroutine s_compute_mixture_sound_speeds(j, k, l) ! ---------------------

        integer, intent(IN) :: j, k, l

        real(kind(0d0)) :: blkmod1, blkmod2 !< Fluid bulk modulus for alternate sound speed

        integer :: i !< Generic loop iterator

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

    end subroutine s_compute_mixture_sound_speeds ! ------------------------


    !>  The procedure computes the Roe average density, velocity,
        !!      enthalpy, volume fractions, specific heat ratio function,
        !!      speed of sound, shear and volume Reynolds numbers, Weber
        !!      numbers and curvatures, at the cell-boundaries, from the
        !!      left and right states of the Riemann problem.
        !! @param j  First coordinate index
        !! @param k Second coordinate index
        !! @param l  Third coordinate index
    subroutine s_compute_roe_average_state(j, k, l) ! ---------------

        integer, intent(IN) :: j, k, l

        integer :: i

        ! Left and Right Riemann Problem States ============================
        do i = 1, cont_idx%end
            alpha_rho_L(i) = qL_prim_rs_vf(i)%sf(j, k, l)
            alpha_rho_R(i) = qR_prim_rs_vf(i)%sf(j + 1, k, l)
        end do

        do i = 1, num_dims
            vel_L(i) = qL_prim_rs_vf(cont_idx%end + i)%sf(j, k, l)
            vel_R(i) = qR_prim_rs_vf(cont_idx%end + i)%sf(j + 1, k, l)
        end do

        pres_L = qL_prim_rs_vf(E_idx)%sf(j, k, l)
        pres_R = qR_prim_rs_vf(E_idx)%sf(j + 1, k, l)

        call s_convert_to_mixture_variables(qL_prim_rs_vf, &
                                            rho_L, gamma_L, &
                                            pi_inf_L, &
                                            j, k, l)
        call s_convert_to_mixture_variables(qR_prim_rs_vf, &
                                            rho_R, gamma_R, &
                                            pi_inf_R, &
                                            j + 1, k, l)

        E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*sum(vel_L**2d0)
        E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*sum(vel_R**2d0)

        H_L = (E_L + pres_L)/rho_L
        H_R = (E_R + pres_R)/rho_R

        call s_compute_mixture_sound_speeds(j, k, l)

        ! ==================================================================

        ! Roe Average Riemann Problem State ================================
        rho_avg = sqrt(rho_L*rho_R)

        vel_avg = (sqrt(rho_L)*vel_L + sqrt(rho_R)*vel_R)/ &
                  (sqrt(rho_L) + sqrt(rho_R))

        H_avg = (sqrt(rho_L)*H_L + sqrt(rho_R)*H_R)/ &
                (sqrt(rho_L) + sqrt(rho_R))

        gamma_avg = (sqrt(rho_L)*gamma_L + sqrt(rho_R)*gamma_R)/ &
                    (sqrt(rho_L) + sqrt(rho_R))

        if (mixture_err) then
            if ((H_avg - 5d-1*sum(vel_avg**2d0)) < 0d0) then
                c_avg = sgm_eps
            else
                c_avg = sqrt((H_avg - 5d-1*sum(vel_avg**2d0))/gamma_avg)
            end if
        else
            c_avg = sqrt((H_avg - 5d-1*sum(vel_avg**2d0))/gamma_avg)
        end if

        ! ==================================================================

    end subroutine s_compute_roe_average_state ! ---------------------------

    !>  This procedure calculates the arithmetic average density,
        !!      velocity, enthalpy, volume fractions, specIFic heat ratio
        !!      function, sound speed, shear and volume Reynolds numbers,
        !!      Weber numbers and the curvatures, at the cell-boundaries,
        !!      from the left and right states of the Riemann problem.
        !!  @param j  First coordinate index
        !!  @param k Second coordinate index
        !!  @param l  Third coordinate index
    subroutine s_compute_arithmetic_average_state(j, k, l) ! --------

        integer, intent(IN) :: j, k, l

        integer :: i, q !< Generic loop iterator

        ! Left and Right Riemann Problem States ============================
        do i = 1, cont_idx%end
            alpha_rho_L(i) = qL_prim_rs_vf(i)%sf(j, k, l)
            alpha_rho_R(i) = qR_prim_rs_vf(i)%sf(j + 1, k, l)
        end do

        do i = 1, num_dims
            vel_L(i) = qL_prim_rs_vf(cont_idx%end + i)%sf(j, k, l)
            vel_R(i) = qR_prim_rs_vf(cont_idx%end + i)%sf(j + 1, k, l)
        end do

        call s_convert_to_mixture_variables(qL_prim_rs_vf, &
                                            rho_L, gamma_L, &
                                            pi_inf_L, &
                                            j, k, l)
        call s_convert_to_mixture_variables(qR_prim_rs_vf, &
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

        ! ==================================================================

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

        ! ==================================================================

    end subroutine s_compute_arithmetic_average_state ! --------------------

    !>  The direct estimation of the left, right and middle wave
        !!      speeds, proposed by Batten et al. (1997) that results in
        !!      the exact resolution of isolated shock and contact waves.
        !!  @param j  First coordinate index
        !!  @param k Second coordinate index
        !!  @param l  Third coordinate index
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

    !>  Estimation of the left, right and star region wave speeds
        !!      by the approximation of the pressures and velocity in the
        !!      star regions, see Toro (1999). The pressures and velocity
        !!      are approximated by using the primitive variables Riemann
        !!      solver (PVRS) and the wave speeds are then estimated from
        !!      those approximations using the exact wave relations.
        !!  @param j  First coordinate index
        !!  @param k Second coordinate index
        !!  @param l  Third coordinate index
    subroutine s_compute_pressure_velocity_wave_speeds(j, k, l) ! ------------

        integer, intent(IN) :: j, k, l

        ! Left and right pressures in the star region
        real(kind(0d0)) :: pres_SL, pres_SR


        ! Left and right shock Mach numbers
        real(kind(0d0)) :: Ms_L, Ms_R

        integer :: i !< Generic loop iterator


        pres_SL = 5d-1*(pres_L + pres_R + rho_avg*c_avg* &
                        (vel_L(dir_idx(1)) - &
                         vel_R(dir_idx(1))))
        pres_SR = pres_SL

        Ms_L = max(1d0, sqrt(1d0 + ((5d-1 + gamma_L)/(1d0 + gamma_L))* &
                             (pres_SL/pres_L - 1d0)*pres_L/ &
                             ((pres_L + pi_inf_L/(1d0 + gamma_L)))))
        Ms_R = max(1d0, sqrt(1d0 + ((5d-1 + gamma_R)/(1d0 + gamma_R))* &
                             (pres_SR/pres_R - 1d0)*pres_R/ &
                             ((pres_R + pi_inf_R/(1d0 + gamma_R)))))

        s_L = vel_L(dir_idx(1)) - c_L*Ms_L
        s_R = vel_R(dir_idx(1)) + c_R*Ms_R

        s_S = 5d-1*((vel_L(dir_idx(1)) + vel_R(dir_idx(1))) + &
                    (pres_L - pres_R)/ &
                    (rho_avg*c_avg))

    end subroutine s_compute_pressure_velocity_wave_speeds ! ---------------

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_riemann_solvers_module() ! ---------------------

        ! Allocating the variables that will be utilized to formulate the
        ! left, right, and average states of the Riemann problem, as well
        ! the Riemann problem solution
        allocate (qL_prim_rs_vf(1:sys_size), qR_prim_rs_vf(1:sys_size))

        allocate (flux_rs_vf(1:sys_size), flux_src_rs_vf(1:sys_size))

        allocate (flux_gsrc_rs_vf(1:sys_size))

        allocate (vel_src_rs_vf(1:num_dims))

        allocate (alpha_rho_L(1:cont_idx%end), vel_L(1:num_dims))
        allocate (alpha_rho_R(1:cont_idx%end), vel_R(1:num_dims))

        allocate (vel_avg(1:num_dims))

        allocate (alpha_L(1:num_fluids))
        allocate (alpha_R(1:num_fluids))


        if (avg_state == 1) then
            s_compute_average_state => s_compute_roe_average_state
        else
            s_compute_average_state => s_compute_arithmetic_average_state
        end if

        if (wave_speeds == 1) then
            s_compute_wave_speeds => s_compute_direct_wave_speeds
        else
            s_compute_wave_speeds => s_compute_pressure_velocity_wave_speeds
        end if


        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables

        if (model_eqns == 1) then        ! Gamma/pi_inf model
            s_convert_to_mixture_variables => &
                s_convert_mixture_to_mixture_variables
        else                            ! Volume fraction model
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables
        end if

    end subroutine s_initialize_riemann_solvers_module ! -------------------


    subroutine s_populate_riemann_states_variables_buffers( & ! ------------
        qL_prim_vf, dqL_prim_dx_vf, &
        gm_alphaL_vf, &
        qR_prim_vf, dqR_prim_dx_vf, &
        gm_alphaR_vf, &
        norm_dir, ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: qL_prim_vf, qR_prim_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(INOUT) :: dqL_prim_dx_vf, dqR_prim_dx_vf, &
                             gm_alphaL_vf, gm_alphaR_vf

        integer, intent(IN) :: norm_dir
        type(bounds_info), intent(IN) :: ix, iy, iz
        integer :: i

        ! Population of Buffers in x-direction =============================
        if (norm_dir == 1) then

            if (bc_x%beg == -4) then    ! Riemann state extrap. BC at beginning

                do i = 1, sys_size
                    qL_prim_vf(i)%sf(-1, iy%beg:iy%end, iz%beg:iz%end) = &
                        qR_prim_vf(i)%sf(0, iy%beg:iy%end, iz%beg:iz%end)
                end do

            end if

            if (bc_x%end == -4) then    ! Riemann state extrap. BC at end

                do i = 1, sys_size
                    qR_prim_vf(i)%sf(m + 1, iy%beg:iy%end, iz%beg:iz%end) = &
                        qL_prim_vf(i)%sf(m, iy%beg:iy%end, iz%beg:iz%end)
                end do

            end if
        end if


    end subroutine s_populate_riemann_states_variables_buffers ! -----------

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

        ! Configuring the coordinate direction indexes and flags
        if (norm_dir == 1) then
            is1 = ix; is2 = iy; is3 = iz
            dir_idx = (/1, 2, 3/); dir_flg = (/1d0, 0d0, 0d0/)
        end if


        ! Setting up special bounds for cell-average values
        xbeg = -buff_size; ybeg = 0; zbeg = 0
        if (n > 0) ybeg = -buff_size; if (p > 0) zbeg = -buff_size
        xend = m - xbeg; yend = n - ybeg; zend = p - zbeg

        ! Configuring the coordinate direction indexes
        if (norm_dir == 1) then
            s1beg = xbeg; s1end = xend; s2beg = ybeg; s2end = yend; s3beg = zbeg; s3end = zend
        end if

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

        ! Reshaping Inputted Data in x-direction ===========================
        if (norm_dir == 1) then

            do i = 1, sys_size
                qL_prim_rs_vf(i)%sf = qL_prim_vf(i)%sf(ix%beg:ix%end, &
                                                       iy%beg:iy%end, &
                                                       iz%beg:iz%end)
                qR_prim_rs_vf(i)%sf = qR_prim_vf(i)%sf(ix%beg + 1:ix%end + 1, &
                                                       iy%beg:iy%end, &
                                                       iz%beg:iz%end)
            end do

        end if
        ! ==================================================================

    end subroutine s_initialize_riemann_solver ! ---------------------------

    subroutine s_finalize_riemann_solver(flux_vf, flux_src_vf, & ! --------
                                         flux_gsrc_vf, &
                                         norm_dir, ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: flux_vf, flux_src_vf, flux_gsrc_vf

        integer, intent(IN) :: norm_dir

        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: i, j, k !< Generic loop iterators

        ! Deallocating Left, Right and Average Riemann Problem States ======
        do i = 1, sys_size
            deallocate (qL_prim_rs_vf(i)%sf, qR_prim_rs_vf(i)%sf)
        end do

        ! ==================================================================

        ! Deallocating Intercell Fluxes and Velocity =======================
        do i = 1, sys_size
            flux_rs_vf(i)%sf => null()
            flux_src_rs_vf(i)%sf => null()
            flux_gsrc_rs_vf(i)%sf => null()
        end do
        vel_src_rs_vf(dir_idx(1))%sf => null()

    end subroutine s_finalize_riemann_solver ! -----------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_riemann_solvers_module() ! -----------------------

        ! Deallocating the variables that were utilized to formulate the
        ! left, right and average states of the Riemann problem, as well
        ! the Riemann problem solution
        deallocate (qL_prim_rs_vf, qR_prim_rs_vf)

        deallocate (flux_rs_vf, flux_src_rs_vf, flux_gsrc_rs_vf)

        deallocate (vel_src_rs_vf)

        deallocate (alpha_rho_L, vel_L)
        deallocate (alpha_rho_R, vel_R)

        deallocate (vel_avg)

        deallocate (alpha_L, alpha_R)

        ! Disassociating procedural pointer to the subroutine which was
        ! utilized to calculate the solution of a given Riemann problem
        s_riemann_solver => null()

        ! Disassociating the procedural pointers to the procedures that were
        ! utilized to compute the average state and estimate the wave speeds
        s_compute_average_state => null(); s_compute_wave_speeds => null()

        ! Disassociating the pointer to the procedure that was utilized to
        ! to convert mixture or species variables to the mixture variables
        s_convert_to_mixture_variables => null()

    end subroutine s_finalize_riemann_solvers_module ! ---------------------

end module m_riemann_solvers
