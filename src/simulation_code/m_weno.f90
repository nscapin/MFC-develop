!! This file only does 1D cases! Careful.
!!
!!
!!
!!
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
!! @file m_weno.f90
!! @brief Contains module m_weno
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief  Weighted essentially non-oscillatory (WENO) reconstruction scheme
!!              that is supplemented with monotonicity preserving bounds (MPWENO)
!!              and a mapping function that boosts the accuracy of the non-linear
!!              weights (WENOM). MPWENO, see Balsara and Shu (2000), prevents the
!!              reconstructed values to lay outside the range set by the stencil,
!!              while WENOM, see Henrick et al. (2005), recovers the formal order
!!              of accuracy of the reconstruction at critical points. Please note
!!              that the basic WENO approach is implemented according to the work
!!              of Jiang and Shu (1996).
module m_weno

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_variables_conversion !< State variables type conversion procedures
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_weno_module, s_weno, s_finalize_weno_module

    type(vector_field), allocatable, dimension(:) :: v_rs_wsL

    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_L
    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: poly_coef_R

    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_L
    real(kind(0d0)), target, allocatable, dimension(:, :) :: d_R

    real(kind(0d0)), target, allocatable, dimension(:, :, :) :: beta_coef

contains

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_weno_module() ! --------------------------------

        type(bounds_info) :: ix, iy, iz !< Indical bounds in the x-, y- and z-directions
        integer :: i

        ! Allocating WENO-stencil for the variables to be WENO-reconstructed
        allocate (v_rs_wsL(-weno_polyn:weno_polyn))

        do i = -weno_polyn, weno_polyn
            allocate (v_rs_wsL(i)%vf(1:sys_size) )
        end do

        ! Allocating/Computing WENO Coefficients in x-direction ============
        ix%beg = -buff_size + weno_polyn; ix%end = m - ix%beg

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

        nullify(s_cb)

    end subroutine s_compute_weno_coefficients ! ---------------------------


    subroutine s_weno(v_vf, vL_vf, vR_vf, weno_dir_dummy, ix, iy, iz)

        type(scalar_field), dimension(:), intent(IN) :: v_vf
        type(scalar_field), dimension(:), intent(INOUT) :: vL_vf, vR_vf
        integer, intent(IN) :: weno_dir_dummy
        type(bounds_info), intent(IN) :: ix, iy, iz

        real(kind(0d0)), dimension(-weno_polyn:weno_polyn-1) :: dvd 
        real(kind(0d0)), dimension(0:weno_polyn) ::  poly_L, poly_R
        real(kind(0d0)), dimension(0:weno_polyn) :: alpha_L, alpha_R
        real(kind(0d0)), dimension(0:weno_polyn) :: omega_L, omega_R
        real(kind(0d0)), dimension(0:weno_polyn) :: beta 

        integer :: i, j, k, l

        call s_initialize_weno(v_vf, ix, iy, iz)

        do i = 1, sys_size
            do l = iz%beg, iz%end
                do k = iy%beg, iy%end
                    do j = ix%beg, ix%end

                        dvd(1) = v_rs_wsL(2)%vf(i)%sf(j, k, l) &
                                 - v_rs_wsL(1)%vf(i)%sf(j, k, l)
                        dvd(0) = v_rs_wsL(1)%vf(i)%sf(j, k, l) &
                                - v_rs_wsL(0)%vf(i)%sf(j, k, l)
                        dvd(-1) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                                - v_rs_wsL(-1)%vf(i)%sf(j, k, l)
                        dvd(-2) = v_rs_wsL(-1)%vf(i)%sf(j, k, l) &
                                - v_rs_wsL(-2)%vf(i)%sf(j, k, l)

                        poly_L(0) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                                  + poly_coef_L(0, 0, j)*dvd(1) &
                                  + poly_coef_L(0, 1, j)*dvd(0)
                        poly_L(1) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                                  + poly_coef_L(1, 0, j)*dvd(0) &
                                  + poly_coef_L(1, 1, j)*dvd(-1)
                        poly_L(2) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
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

                        alpha_L = d_L(:, j)/(beta*beta)
                        omega_L = alpha_L/sum(alpha_L)

                        dvd(1) = v_rs_wsL(2)%vf(i)%sf(j, k, l) &
                               - v_rs_wsL(1)%vf(i)%sf(j, k, l)
                        dvd(0) = v_rs_wsL(1)%vf(i)%sf(j, k, l) &
                               - v_rs_wsL(0)%vf(i)%sf(j, k, l)
                        dvd(-1) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                                - v_rs_wsL(-1)%vf(i)%sf(j, k, l)
                        dvd(-2) = v_rs_wsL(-1)%vf(i)%sf(j, k, l) &
                                - v_rs_wsL(-2)%vf(i)%sf(j, k, l)

                        poly_R(0) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                                  + poly_coef_R(0, 0, j)*dvd(1) &
                                  + poly_coef_R(0, 1, j)*dvd(0)
                        poly_R(1) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
                                  + poly_coef_R(1, 0, j)*dvd(0) &
                                  + poly_coef_R(1, 1, j)*dvd(-1)
                        poly_R(2) = v_rs_wsL(0)%vf(i)%sf(j, k, l) &
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

                        alpha_R = d_R(:, j)/(beta*beta)
                        omega_R = alpha_R/sum(alpha_R)

                        vL_vf(i)%sf(j, k, l) = sum(omega_L*poly_L)
                        vR_vf(i)%sf(j, k, l) = sum(omega_R*poly_R)

                    end do
                end do
            end do
        end do

        call s_finalize_weno()

    end subroutine s_weno 


    subroutine s_initialize_weno(v_vf, ix, iy, iz)

        type(scalar_field), dimension(:), intent(IN) :: v_vf
        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: i, j, k

        ! Allocate space for full stencil variables
        do i = -weno_polyn, weno_polyn
            do j = 1, sys_size
                allocate (v_rs_wsL(i)%vf(j)%sf(ix%beg:ix%end, &
                                               iy%beg:iy%end, &
                                               iz%beg:iz%end))
            end do
        end do

        ! Populate variable buffers at each point (for full stencil)
        do i = -weno_polyn, weno_polyn
            do j = 1, sys_size
                do k = ix%beg, ix%end
                    v_rs_wsL(i)%vf(j)%sf(k, :, :) = &
                        v_vf(j)%sf(i + k, iy%beg:iy%end, iz%beg:iz%end)
                end do
            end do
        end do

    end subroutine s_initialize_weno ! -------------------------------------


    subroutine s_finalize_weno()

        integer :: i, j

        do i = -weno_polyn, weno_polyn
            do j = 1, sys_size
                deallocate (v_rs_wsL(i)%vf(j)%sf)
            end do
        end do

    end subroutine s_finalize_weno ! ---------------------------------------


    subroutine s_finalize_weno_module() ! ----------------------------------

        deallocate (v_rs_wsL)
        deallocate (poly_coef_L, poly_coef_R)
        deallocate (d_L, d_R)
        deallocate (beta_coef)

    end subroutine s_finalize_weno_module ! --------------------------------

end module m_weno
