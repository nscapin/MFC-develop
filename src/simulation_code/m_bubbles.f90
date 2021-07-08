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
!! @file m_bubbles.f90
!! @brief Contains module m_bubbles
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module is used to compute the ensemble-averaged bubble dynamic variables
module m_bubbles

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    ! ==========================================================================

    implicit none

contains

    subroutine s_compute_bubble_source(q_prim_vf, q_cons_vf, &
                                       bub_adv_src, bub_r_src, bub_v_src, nb)

        real(kind(0d0)), dimension(-4:,0:,0:,1:) :: q_prim_vf, q_cons_vf

        real(kind(0d0)), dimension(0:,0:,0:),    intent(INOUT) ::  bub_adv_src
        real(kind(0d0)), dimension(0:,0:,0:,1:), intent(INOUT) :: bub_r_src, bub_v_src
        integer, intent(in) :: nb

        real(kind(0d0)) :: c_liquid, Cpbw, Cpinf, rddot
        real(kind(0d0)) :: n_tait, B_tait
        real(kind(0d0)) :: nbub

        real(kind(0d0)), dimension(nb)  :: Rtmp, Vtmp, qq_temp
        real(kind(0d0)) :: myR, myV, alf, myP, myRho, R2Vav

        integer :: j, k, l, q

        n_tait = gammas(1)
        B_tait = pi_infs(1)
        n_tait = 1.d0/n_tait + 1.d0

        !$acc data present(q_prim_vf, q_cons_vf, bub_idx_rs, bub_idx_vs, bub_adv_src, bub_r_src, bub_v_src) 
        !$acc parallel loop collapse(3) gang vector private(Rtmp, Vtmp, qq_temp)
        do l = 0, p
            do k = 0, n
                do j = 0, m

                    do q = 1, nb
                        Rtmp(q) = q_prim_vf(j, k, l, bub_idx_rs(q))
                        Vtmp(q) = q_prim_vf(j, k, l, bub_idx_vs(q))
                        qq_temp(q) = (Rtmp(q)**2d0)*Vtmp(q)
                    end do

                    call s_comp_n_from_prim(q_prim_vf(j, k, l, alf_idx), Rtmp, nbub)

                    call s_quad_bub(qq_temp, R2Vav)
                    bub_adv_src(j, k, l) = 4.d0*pi*nbub*R2Vav

                    do q = 1, nb
                        bub_r_src(j, k, l, q) = q_cons_vf(j, k, l, bub_idx_vs(q))

                        myRho = q_prim_vf(j, k, l, 1)
                        myP   = q_prim_vf(j, k, l, E_idx)
                        alf   = q_prim_vf(j, k, l, alf_idx)
                        myR   = q_prim_vf(j, k, l, bub_idx_rs(q))
                        myV   = q_prim_vf(j, k, l, bub_idx_vs(q))

                        Cpbw = f_cpbw_KM(R0(q), myR, myV)
                        rddot = f_rddot_RP(myP, myRho, myR, myV, R0(q), Cpbw)
                        bub_v_src(j, k, l, q) = nbub*rddot
                    end do
                end do
            end do
        end do
        !$acc end parallel loop
        !$acc end data


    end subroutine s_compute_bubble_source


    function f_rddot_RP(fCp, fRho, fR, fV, fR0, fCpbw)
        !$acc routine seq

        real(kind(0d0)), intent(IN) :: fCp, fRho, fR, fV, fR0, fCpbw
        real(kind(0d0))             :: f_rddot_RP

        f_rddot_RP = (-1.5d0*(fV**2d0) + (fCpbw - fCp)/fRho)/fR
        ! if (Re_inv /= dflt_real) f_rddot_RP = f_rddot_RP - 4d0*Re_inv*fv/(fr**2d0)/fRho

    end function f_rddot_RP


    function f_cpbw_KM(fR0, fR, fV)
        !$acc routine seq

        real(kind(0d0)), intent(IN) :: fR0, fR, fV
        real(kind(0d0))             :: f_cpbw_KM

        f_cpbw_KM = Ca*((fR0/fR)**(3.d0*gam)) - Ca + 1d0
            ! if (Web /= dflt_real) f_cpbw_KM = f_cpbw_KM + &
            !                                   (2.d0/(Web*fR0))*((fR0/fR)**(3.d0*gam))

        ! if (Web /= dflt_real) f_cpbw_KM = f_cpbw_KM - 2.d0/(fR*Web)
        ! if (Re_inv /= dflt_real) f_cpbw_KM = f_cpbw_KM - 4.d0*Re_inv*fV/fR

    end function f_cpbw_KM


    function f_rddot_KM(fpbdot, fCp, fCpbw, fRho, fR, fV, fR0, fC)

        real(kind(0d0)), intent(IN) :: fpbdot, fCp, fCpbw
        real(kind(0d0)), intent(IN) :: fRho, fR, fV, fR0, fC

        real(kind(0d0))             :: tmp1, tmp2, cdot_star
        real(kind(0d0))             :: f_rddot_KM

        cdot_star = -3d0*gam*Ca*((fR0/fR)**(3d0*gam))*fV/fR
        if (Web /= dflt_real) cdot_star = cdot_star - &
                                          3d0*gam*(2d0/(Web*fR0))*((fR0/fR)**(3d0*gam))*fV/fR

        if (Web /= dflt_real) cdot_star = cdot_star + (2d0/Web)*fV/(fR**2d0)
        if (Re_inv /= dflt_real) cdot_star = cdot_star + 4d0*Re_inv*((fV/fR)**2d0)

        tmp1 = fV/fC
        tmp2 = 1.5d0*(fV**2d0)*(tmp1/2d0 - 1d0) + &
               (1d0 + tmp1)*(fCpbw - fCp)/fRho + &
               cdot_star*fR/(fRho*fC)

        if (Re_inv == dflt_real) then
            f_rddot_KM = tmp2/(fR*(1d0 - tmp1))
        else
            f_rddot_KM = tmp2/(fR*(1d0 - tmp1) + 4d0*Re_inv/(fRho*fC))
        end if

    end function f_rddot_KM

    subroutine s_quad_bub(func, mom)
        !$acc routine seq

        real(kind(0d0)), dimension(nb), intent(in) :: func
        real(kind(0d0)), intent(out) :: mom
        integer :: i

        mom = 0d0
        do i = 1,nb
            mom = mom + weight(i)*func(i)
        end do

    end subroutine s_quad_bub

end module m_bubbles
