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
!! @file m_variables_conversion.f90
!! @brief Contains module m_variables_conversion
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module features a database of subroutines that allow for the
!!              conversion of state variables from one type into another. At this
!!              time, the state variables type conversions below are available:
!!                             1) Mixture        => Mixture
!!                             2) Species        => Mixture
!!                             3) Conservative   => Primitive
!!                             5) Conservative   => Flux
!!                             6) Primitive      => Conservative
!!                             8) Primitive      => Flux
module m_variables_conversion

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use nvtx
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_variables_conversion_module, &
         s_convert_to_mixture_variables, &
         s_convert_mixture_to_mixture_variables, &
         s_convert_species_to_mixture_variables_bubbles, &
         s_convert_species_to_mixture_variables, &
         s_convert_conservative_to_primitive_variables, &
         s_convert_conservative_to_primitive_variables_acc, &
         s_convert_conservative_to_flux_variables, &
         s_convert_primitive_to_conservative_variables, &
         s_convert_primitive_to_flux_variables, &
         s_convert_primitive_to_flux_variables_bubbles, &
         s_finalize_variables_conversion_module

    integer, dimension(:), allocatable :: bub_idx_rs
    !$acc declare create(bub_idx_rs)

    abstract interface ! =======================================================

        !> The abstract interface to the procedures that are utilized to convert
        !! the mixture and the species variables into the mixture variables. For
        !! more information, refer to:
        !!               1) s_convert_mixture_to_mixture_variables
        !!               2) s_convert_species_to_mixture_variables
        !! @param qK_vf Conservative/primitive variables
        !! @param rho_K Mixture density
        !! @param gamma_K Mixture sp. heat ratio
        !! @param pi_inf_K Mixture stiffness function
        !! @param i Cell location first index
        !! @param j Cell location second index
        !! @param k Cell location third index
        subroutine s_convert_abstract_to_mixture_variables(qK_vf, rho_K, &
                                                           gamma_K, pi_inf_K, &
                                                           i, j, k)

            import :: scalar_field, sys_size, num_fluids

            type(scalar_field), dimension(sys_size), intent(IN) :: qK_vf

            real(kind(0d0)), intent(OUT) :: rho_K, gamma_K, pi_inf_K

            integer, intent(IN) :: i, j, k

        end subroutine s_convert_abstract_to_mixture_variables

    end interface ! ============================================================


    procedure(s_convert_abstract_to_mixture_variables), &
        pointer :: s_convert_to_mixture_variables => null() !<

contains

    !>  The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_variables_conversion_module() ! ----------------

        integer :: i 

        ! Associating the procedural pointer to the appropriate subroutine
        ! that will be utilized in the conversion to the mixture variables
        if (model_eqns == 1) then        ! gamma/pi_inf model
            s_convert_to_mixture_variables => &
                s_convert_mixture_to_mixture_variables
        else                            ! Volume fraction model
            s_convert_to_mixture_variables => &
                s_convert_species_to_mixture_variables
        end if

        allocate(bub_idx_rs(1:nb))

        do i = 1,nb
            bub_idx_rs(i) = bub_idx%rs(i)
        end do

    end subroutine s_initialize_variables_conversion_module ! --------------

    !> This procedure is used alongside with the gamma/pi_inf
        !!      model to transfer the density, the specific heat ratio
        !!      function and liquid stiffness function from the vector
        !!      of conservative or primitive variables to their scalar
        !!      counterparts.
        !! @param qK_vf conservative or primitive variables
        !! @param i cell index to transfer mixture variables
        !! @param j cell index to transfer mixture variables
        !! @param k cell index to transfer mixture variables
        !! @param rho_K density
        !! @param gamma_K  specific heat ratio function
        !! @param pi_inf_K liquid stiffness
    subroutine s_convert_mixture_to_mixture_variables(qK_vf, rho_K, &
                                                      gamma_K, pi_inf_K, &
                                                      i, j, k)

        type(scalar_field), dimension(sys_size), intent(IN) :: qK_vf

        real(kind(0d0)), intent(OUT) :: rho_K, gamma_K, pi_inf_K


        integer, intent(IN) :: i, j, k

        ! Performing the transfer of the density, the specific heat ratio
        ! function as well as the liquid stiffness function, respectively
        rho_K = qK_vf(1)%sf(i, j, k)
        gamma_K = qK_vf(gamma_idx)%sf(i, j, k)
        pi_inf_K = qK_vf(pi_inf_idx)%sf(i, j, k)

    end subroutine s_convert_mixture_to_mixture_variables ! ----------------

    !>  This procedure is used alongside with the gamma/pi_inf
        !!      model to transfer the density, the specific heat ratio
        !!      function and liquid stiffness function from the vector
        !!      of conservative or primitive variables to their scalar
        !!      counterparts. Specifially designed for when subgrid bubbles
        !!      must be included.
        !! @param qK_vf primitive variables
        !! @param rho_K density
        !! @param gamma_K specific heat ratio
        !! @param pi_inf_K liquid stiffness
        !! @param i Cell index
        !! @param j Cell index
        !! @param k Cell index
    subroutine s_convert_species_to_mixture_variables_bubbles(qK_vf, rho_K, &
                                                              gamma_K, pi_inf_K, &
                                                              i, j, k)

        type(scalar_field), dimension(sys_size), intent(IN) :: qK_vf

        real(kind(0d0)), intent(OUT) :: rho_K, gamma_K, pi_inf_K

        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_K, alpha_K  !<
            !! Partial densities and volume fractions

        integer, intent(IN) :: i, j, k
        integer :: l

        ! Constraining the partial densities and the volume fractions within
        ! their physical bounds to make sure that any mixture variables that
        ! are derived from them result within the limits that are set by the
        ! fluids physical parameters that make up the mixture
        ! alpha_rho_K(1) = qK_vf(i)%sf(i,j,k)
        ! alpha_K(1)     = qK_vf(E_idx+i)%sf(i,j,k)

        ! Performing the transfer of the density, the specific heat ratio
        ! function as well as the liquid stiffness function, respectively
        if (model_eqns == 4) then
            rho_K = qK_vf(1)%sf(i, j, k)
            gamma_K = fluid_pp(1)%gamma !qK_vf(gamma_idx)%sf(i,j,k)
            pi_inf_K = fluid_pp(1)%pi_inf !qK_vf(pi_inf_idx)%sf(i,j,k)
        else if ((model_eqns == 2) .and. bubbles .and. adv_alphan) then
            rho_k = 0d0; gamma_k = 0d0; pi_inf_k = 0d0

            if (mpp_lim .and. (num_fluids > 2)) then
                do l = 1, num_fluids
                    rho_k = rho_k + qK_vf(l)%sf(i, j, k)
                    gamma_k = gamma_k + qK_vf(l + E_idx)%sf(i, j, k)*fluid_pp(l)%gamma
                    pi_inf_k = pi_inf_k + qK_vf(l + E_idx)%sf(i, j, k)*fluid_pp(l)%pi_inf
                end do
            else if (num_fluids == 2) then
                rho_K = qK_vf(1)%sf(i, j, k)
                gamma_K = fluid_pp(1)%gamma
                pi_inf_K = fluid_pp(1)%pi_inf
            else if (num_fluids > 2) then
                do l = 1, num_fluids - 1 !leave out bubble part of mixture
                    rho_k = rho_k + qK_vf(l)%sf(i, j, k)
                    gamma_k = gamma_k + qK_vf(l + E_idx)%sf(i, j, k)*fluid_pp(l)%gamma
                    pi_inf_k = pi_inf_k + qK_vf(l + E_idx)%sf(i, j, k)*fluid_pp(l)%pi_inf
                end do
                !rho_K    = qK_vf(1)%sf(i,j,k)
                !gamma_K  = fluid_pp(1)%gamma
                !pi_inf_K = fluid_pp(1)%pi_inf
            else
                rho_K = qK_vf(1)%sf(i, j, k)
                gamma_K = fluid_pp(1)%gamma
                pi_inf_K = fluid_pp(1)%pi_inf
            end if
        end if

    end subroutine s_convert_species_to_mixture_variables_bubbles ! ----------------

    !>  This subroutine is designed for the volume fraction model
        !!              and provided a set of either conservative or primitive
        !!              variables, computes the density, the specific heat ratio
        !!              function and the liquid stiffness function from q_vf and
        !!              stores the results into rho, gamma and pi_inf.
        !! @param qK_vf primitive variables
        !! @param rho_K density
        !! @param gamma_K specific heat ratio
        !! @param pi_inf_K liquid stiffness
        !! @param k Cell index
        !! @param l Cell index
        !! @param r Cell index
    subroutine s_convert_species_to_mixture_variables(qK_vf, rho_K, &
                                                      gamma_K, pi_inf_K, &
                                                      k, l, r)

        !! This is the one that gets called

        type(scalar_field), dimension(sys_size), intent(IN) :: qK_vf

        real(kind(0d0)), intent(OUT) :: rho_K, gamma_K, pi_inf_K

        real(kind(0d0)), dimension(num_fluids) :: alpha_rho_K, alpha_K 

        integer, intent(IN) :: k, l, r

        integer :: i, j

        do i = 1, num_fluids
            alpha_rho_K(i) = qK_vf(i)%sf(k, l, r)
            alpha_K(i) = qK_vf(adv_idx%beg + i - 1)%sf(k, l, r)
        end do

        rho_K = 0d0; gamma_K = 0d0; pi_inf_K = 0d0

        do i = 1, num_fluids
            rho_K = rho_K + alpha_rho_K(i)
            gamma_K = gamma_K + alpha_K(i)*fluid_pp(i)%gamma
            pi_inf_K = pi_inf_K + alpha_K(i)*fluid_pp(i)%pi_inf
        end do

    end subroutine s_convert_species_to_mixture_variables ! ----------------


    subroutine s_convert_conservative_to_primitive_variables_acc( &
                                                             qK_cons_vf_flat, &
                                                             qK_prim_vf_flat, &
                                                             ix, iy, iz)

        !! This -4 is because you can't pass arrays with negative indexing consistently
        !! and it only applys to WENO5 (buffsize = 4), also the other directions start/end at 0
        real(kind(0d0)), dimension(-4:,0:,0:,1:), intent(INOUT) :: qK_cons_vf_flat
        real(kind(0d0)), dimension(-4:,0:,0:,1:), intent(INOUT) :: qK_prim_vf_flat

        type(bounds_info), intent(IN) :: ix, iy, iz

        real(kind(0d0)), dimension(10) :: alpha_rho
        real(kind(0d0)), dimension(10) :: alpha
        ! real(kind(0d0)), dimension(10) :: gammas, pi_infs
        real(kind(0d0)) ::      rho_K
        real(kind(0d0)) :: dyn_pres_K
        real(kind(0d0)) ::    gamma_K
        real(kind(0d0)) ::   pi_inf_K

        integer :: ixb, ixe, iyb, iye, izb, ize
        integer :: cont_idx_e
        integer :: adv_idx_b, adv_idx_e
        integer :: mom_idx_b, mom_idx_e
        integer :: i, j, k, l 

        cont_idx_e = cont_idx%end
        adv_idx_b = adv_idx%beg
        adv_idx_e = adv_idx%end
        mom_idx_b = mom_idx%beg
        mom_idx_e = mom_idx%end

        ixb = ix%beg; ixe = ix%end
        iyb = iy%beg; iye = iy%end
        izb = iz%beg; ize = iz%end

        ! print*, '-1 idx', qK_cons_vf_flat(-1,0,0,1)

        ! print*, 'c2p: ixb,ixe', ixb, ixe
        ! print*, 'shape qK_cons_vf_flat:', shape(qK_cons_vf_flat)


        !$acc data present(qK_cons_vf_flat, qK_prim_vf_flat)
        !$acc parallel loop collapse(3) gang vector private(alpha_rho, alpha) 
        do l = izb, ize
            do k = iyb, iye
                do j = ixb, ixe

                    do i = 1, num_fluids
                        alpha_rho(i) = qK_cons_vf_flat(j, k, l, i)
                        alpha(i) = qK_cons_vf_flat(j, k, l, adv_idx_b + i - 1)
                    end do

                    call s_convert_species_to_mixture_variables_acc( &
                                                      alpha_rho, &
                                                      alpha, &
                                                      rho_K, &
                                                      gamma_K, pi_inf_K, &
                                                      num_fluids &
                                                      )
                    ! print*, 'Var conv', rho_K, pi_inf_K, gamma_K

                    dyn_pres_K = 0d0
                    do i = mom_idx_b, mom_idx_e
                        qK_prim_vf_flat(j, k, l, i) = &
                            qK_cons_vf_flat(j, k, l, i)/max(rho_K, sgm_eps)
                        dyn_pres_K = dyn_pres_K + 5d-1*qK_cons_vf_flat(j, k, l, i) & 
                            * qK_prim_vf_flat(j, k, l, i)
                    end do
                    ! print*, '', rho_K, pi_inf_K, gamma_K, alpha(1), alpha_rho(1), dyn_pres_K,


                    qK_prim_vf_flat(j, k, l, E_idx) = ( &
                        qK_cons_vf_flat(j, k, l, E_idx) - dyn_pres_K - pi_inf_K )/gamma_K

                    do i = 1, cont_idx_e
                        qK_prim_vf_flat(j, k, l, i) = qK_cons_vf_flat(j, k, l, i)
                    end do

                    do i = adv_idx_b, adv_idx_e
                        qK_prim_vf_flat(j, k, l, i) = qK_cons_vf_flat(j, k, l, i)
                    end do

                end do
            end do
        end do
        !$acc end parallel loop 
        !$acc end data


    end subroutine s_convert_conservative_to_primitive_variables_acc


    subroutine s_convert_conservative_to_primitive_variables_bubbles_acc( &
                                                             qK_cons_vf_flat, &
                                                             qK_prim_vf_flat, &
                                                             ix, iy, iz)

        !! This -4 is because you can't pass arrays with negative indexing consistently
        !! and it only applys to WENO5 (buffsize = 4), also the other directions start/end at 0
        real(kind(0d0)), dimension(-4:,0:,0:,1:), intent(INOUT) :: qK_cons_vf_flat
        real(kind(0d0)), dimension(-4:,0:,0:,1:), intent(INOUT) :: qK_prim_vf_flat

        type(bounds_info), intent(IN) :: ix, iy, iz

        real(kind(0d0)), dimension(10) :: alpha_rho
        real(kind(0d0)), dimension(10) :: alpha
        real(kind(0d0)) ::      rho_K
        real(kind(0d0)) :: dyn_pres_K
        real(kind(0d0)) ::    gamma_K
        real(kind(0d0)) ::   pi_inf_K
        real(kind(0d0)), dimension(500) :: nRtmp
        real(kind(0d0)) ::   nbub
        

        integer :: ixb, ixe, iyb, iye, izb, ize
        integer :: cont_idx_e
        integer :: adv_idx_b, adv_idx_e
        integer :: mom_idx_b, mom_idx_e
        integer :: bub_idx_b, bub_idx_e
        integer :: i, j, k, l 

        cont_idx_e = cont_idx%end
        adv_idx_b = adv_idx%beg
        adv_idx_e = adv_idx%end
        mom_idx_b = mom_idx%beg
        mom_idx_e = mom_idx%end
        bub_idx_b = bub_idx%beg
        bub_idx_e = bub_idx%end

        ixb = ix%beg; ixe = ix%end
        iyb = iy%beg; iye = iy%end
        izb = iz%beg; ize = iz%end

        ! print*, '-1 idx', qK_cons_vf_flat(-1,0,0,1)

        ! print*, 'c2p: ixb,ixe', ixb, ixe
        ! print*, 'shape qK_cons_vf_flat:', shape(qK_cons_vf_flat)


        !$acc data present(qK_cons_vf_flat, qK_prim_vf_flat, bub_idx_rs)
        !$acc parallel loop collapse(3) gang vector private(alpha_rho, alpha, nRtmp) 
        do l = izb, ize
            do k = iyb, iye
                do j = ixb, ixe

                    do i = 1, num_fluids
                        alpha_rho(i) = qK_cons_vf_flat(j, k, l, i)
                        alpha(i) = qK_cons_vf_flat(j, k, l, adv_idx_b + i - 1)
                    end do

                    call s_convert_species_to_mixture_variables_bubbles_acc( &
                                                      alpha_rho, &
                                                      alpha, &
                                                      rho_K, &
                                                      gamma_K, pi_inf_K, &
                                                      num_fluids &
                                                      )

                    dyn_pres_K = 0d0
                    do i = mom_idx_b, mom_idx_e
                        qK_prim_vf_flat(j, k, l, i) = &
                            qK_cons_vf_flat(j, k, l, i)/max(rho_K, sgm_eps)
                        dyn_pres_K = dyn_pres_K + 5d-1*qK_cons_vf_flat(j, k, l, i) & 
                            * qK_prim_vf_flat(j, k, l, i)
                    end do

                    ! bubble pressure formulation
                    qK_prim_vf_flat(j, k, l, E_idx) = ( &
                        ( qK_cons_vf_flat(j, k, l, E_idx) - dyn_pres_K ) / &
                        ( 1.d0 - qK_cons_vf_flat(j,k,l,alf_idx) ) &
                        - pi_inf_K )/gamma_K

                    do i = 1, cont_idx_e
                        qK_prim_vf_flat(j, k, l, i) = qK_cons_vf_flat(j, k, l, i)
                    end do

                    do i = adv_idx_b, adv_idx_e
                        qK_prim_vf_flat(j, k, l, i) = qK_cons_vf_flat(j, k, l, i)
                    end do

                    ! bubbles
                    do i = 1,nb
                        nRtmp(i) = qK_cons_vf_flat(j,k,l,bub_idx_rs(i))
                    end  do
                    CALL s_comp_n_from_cons( qK_cons_vf_flat(j, k, l, alf_idx), nRtmp, nbub)
                    do i = bub_idx_b, bub_idx_e
                        qK_prim_vf_flat(j,k,l,i) = qK_cons_vf_flat(j,k,l,i)/nbub
                    end do
                end do
            end do
        end do
        !$acc end parallel loop 
        !$acc end data


    end subroutine s_convert_conservative_to_primitive_variables_bubbles_acc


    subroutine s_convert_species_to_mixture_variables_acc( &
                                                      alpha_rho_K, &
                                                      alpha_K, &
                                                      rho_K, &
                                                      gamma_K, pi_inf_K, &
                                                      num_fluids &
                                                      )
        !$acc routine seq

        real(kind(0d0)), dimension(num_fluids), intent(IN) :: alpha_rho_K, alpha_K 
        real(kind(0d0)), intent(OUT) :: rho_K, gamma_K, pi_inf_K
        integer, intent(IN) :: num_fluids
        integer :: i

        rho_K = 0d0; gamma_K = 0d0; pi_inf_K = 0d0

        do i = 1, num_fluids
            rho_K = rho_K + alpha_rho_K(i)
            gamma_K = gamma_K + alpha_K(i)*gammas(i)
            pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
        end do

    end subroutine s_convert_species_to_mixture_variables_acc

    !> The following procedure handles the conversion between
        !!      the conservative variables and the primitive variables.
        !! @param qK_cons_vf Conservative variables
        !! @param qK_prim_vf Primitive variables
        !! @param gm_alphaK_vf Gradient magnitude of the volume fraction
        !! @param ix Index bounds in first coordinate direction
        !! @param iy Index bounds in second coordinate direction
        !! @param iz Index bounds in third coordinate direction
    subroutine s_convert_conservative_to_primitive_variables(qK_cons_vf, &
                                                             qK_prim_vf, &
                                                             gm_alphaK_vf, &
                                                             ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: qK_cons_vf, qK_prim_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(IN) :: gm_alphaK_vf

        type(bounds_info), intent(IN) :: ix, iy, iz

        real(kind(0d0))                                   ::      rho_K
        real(kind(0d0))                                   :: dyn_pres_K
        real(kind(0d0))                                   ::    gamma_K
        real(kind(0d0))                                   ::   pi_inf_K
        real(kind(0d0))                                   ::       nbub

        integer :: i, j, k, l 

        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg, ix%end

                    dyn_pres_K = 0d0

                    call s_convert_species_to_mixture_variables(qK_cons_vf, rho_K, &
                                                        gamma_K, pi_inf_K, &
                                                        j, k, l)

                    do i = mom_idx%beg, mom_idx%end
                        qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l) &
                                                    /max(rho_K, sgm_eps)
                        dyn_pres_K = dyn_pres_K + 5d-1*qK_cons_vf(i)%sf(j, k, l) &
                                     *qK_prim_vf(i)%sf(j, k, l)
                    end do

                    qK_prim_vf(E_idx)%sf(j, k, l) = ( &
                        qK_cons_vf(E_idx)%sf(j, k, l) - dyn_pres_K - pi_inf_K )/gamma_K
                end do
            end do
        end do

    end subroutine s_convert_conservative_to_primitive_variables ! ---------


    !>  The following subroutine handles the conversion between
        !!      the conservative variables and the Euler flux variables.
        !!  @param qK_cons_vf Conservative variables
        !!  @param FK_vf Flux variables
        !!  @param FK_src_vf Flux source variables
        !!  @param ix Index bounds in the first coordinate direction
        !!  @param iy Index bounds in the second coordinate direction
        !!  @param iz Index bounds in the third coordinate direction
    subroutine s_convert_conservative_to_flux_variables(qK_cons_vf, & ! ---
                                                        FK_vf, &
                                                        FK_src_vf, &
                                                        ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: qK_cons_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: FK_vf, FK_src_vf

        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: j, k, l !< Generic loop iterators

        ! Calculating the flux variables from the conservative ones, without
        ! accounting for the contribution of either viscosity or capillarity
        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg, ix%end

                    if (proc_rank == 0) then
                        print '(A)', 'Conversion from conservative to '// &
                            'flux variables not implemented. '// &
                            'Exiting ...'
                        call s_mpi_abort()
                    end if

                end do
            end do
        end do

    end subroutine s_convert_conservative_to_flux_variables ! --------------

    !>  The following procedure handles the conversion between
        !!      the primitive variables and the conservative variables.
        !!  @param qK_prim_vf Primitive variables
        !!  @param qK_cons_vf Conservative variables
        !!  @param gm_alphaK_vf Gradient magnitude of the volume fractions
        !!  @param ix Index bounds in the first coordinate direction
        !!  @param iy Index bounds in the second coordinate direction
        !!  @param iz Index bounds in the third coordinate direction
    subroutine s_convert_primitive_to_conservative_variables(qK_prim_vf, &
                                                             qK_cons_vf, &
                                                             gm_alphaK_vf, &
                                                             ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: qK_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: qK_cons_vf

        type(scalar_field), &
            allocatable, dimension(:), &
            intent(IN) :: gm_alphaK_vf

        type(bounds_info), intent(IN) :: ix, iy, iz

        integer :: j, k, l !< Generic loop iterators

        ! Calculating the momentum and energy from the velocity and pressure
        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg, ix%end

                    if (proc_rank == 0) then
                        print '(A)', 'Conversion from primitive to '// &
                            'conservative variables not '// &
                            'implemented. Exiting ...'
                        call s_mpi_abort()
                    end if

                end do
            end do
        end do

    end subroutine s_convert_primitive_to_conservative_variables ! ---------


    !>  The following subroutine handles the conversion between
        !!      the primitive variables and the Eulerian flux variables.
        !!  @param qK_prim_vf Primitive variables
        !!  @param FK_vf Flux variables
        !!  @param FK_src_vf Flux source variables
        !!  @param ix Index bounds in the first coordinate direction
        !!  @param iy Index bounds in the second coordinate direction
        !!  @param iz Index bounds in the third coordinate direction
    subroutine s_convert_primitive_to_flux_variables(qK_prim_vf, & ! ------
                                                     FK_vf, &
                                                     FK_src_vf, &
                                                     ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: qK_prim_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: FK_vf, FK_src_vf

        type(bounds_info), intent(IN) :: ix, iy, iz

        ! Partial densities, density, velocity, pressure, energy, advection
        ! variables, the specific heat ratio and liquid stiffness functions,
        ! the shear and volume Reynolds numbers and the Weber numbers
        real(kind(0d0)), dimension(cont_idx%end)          :: alpha_rho_K
        real(kind(0d0))                                   ::       rho_K
        real(kind(0d0)), dimension(num_dims)              ::       vel_K
        real(kind(0d0))                                   ::      pres_K
        real(kind(0d0))                                   ::         E_K
        real(kind(0d0)), dimension(adv_idx%end - E_idx)     ::       adv_K
        real(kind(0d0))                                   ::     gamma_K
        real(kind(0d0))                                   ::    pi_inf_K

        integer :: i, j, k, l !< Generic loop iterators

        ! Computing the flux variables from the primitive variables, without
        ! accounting for the contribution of either viscosity or capillarity
        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg, ix%end

                    do i = 1, cont_idx%end
                        alpha_rho_K(i) = qK_prim_vf(i)%sf(j, k, l)
                    end do

                    do i = 1, num_dims
                        vel_K(i) = qK_prim_vf(cont_idx%end + i)%sf(j, k, l)
                    end do

                    pres_K = qK_prim_vf(E_idx)%sf(j, k, l)

                    call s_convert_to_mixture_variables(qK_prim_vf, rho_K, &
                                                        gamma_K, pi_inf_K, &
                                                        j, k, l)

                    ! Computing the energy from the pressure
                    E_K = gamma_K*pres_K + pi_inf_K &
                          + 5d-1*rho_K*sum(vel_K**2d0)

                    do i = 1, adv_idx%end - E_idx
                        adv_K(i) = qK_prim_vf(E_idx + i)%sf(j, k, l)
                    end do

                    ! mass flux, this should be \alpha_i \rho_i u_i
                    do i = 1, cont_idx%end
                        FK_vf(i)%sf(j, k, l) = alpha_rho_K(i)*vel_K(dir_idx(1))
                    end do

                    do i = 1, num_dims
                        FK_vf(cont_idx%end + dir_idx(i))%sf(j, k, l) = &
                            rho_K*vel_K(dir_idx(1)) &
                            *vel_K(dir_idx(i)) &
                            + pres_K*dir_flg(dir_idx(i))
                    end do

                    ! energy flux, u(E+p)
                    FK_vf(E_idx)%sf(j, k, l) = vel_K(dir_idx(1))*(E_K + pres_K)

                    ! have been using == 2
                    if (riemann_solver == 1) then

                        do i = adv_idx%beg, adv_idx%end
                            FK_vf(i)%sf(j, k, l) = 0d0
                            FK_src_vf(i)%sf(j, k, l) = adv_K(i - E_idx)
                        end do

                    else
                        ! Could be bubbles!
                        do i = adv_idx%beg, adv_idx%end
                            FK_vf(i)%sf(j, k, l) = vel_K(dir_idx(1))*adv_K(i - E_idx)
                        end do

                        FK_src_vf(adv_idx%beg)%sf(j, k, l) = vel_K(dir_idx(1))
                    end if
                end do
            end do
        end do

    end subroutine s_convert_primitive_to_flux_variables ! -----------------

    !>  The following subroutine handles the conversion between
        !!      the primitive variables and the Eulerian flux variables
        !!      for cases with ensemble bubble modeling.
        !!  @param qK_prim_vf Primitive variables
        !!  @param qK_cons_vf Primitive variables
        !!  @param FK_vf Flux variables
        !!  @param FK_src_vf Flux source variables
        !!  @param ix Index bounds in the first coordinate direction
        !!  @param iy Index bounds in the second coordinate direction
        !!  @param iz Index bounds in the third coordinate direction
    subroutine s_convert_primitive_to_flux_variables_bubbles(qK_prim_vf, & ! ------
                                                             qk_cons_vf, &
                                                             FK_vf, &
                                                             FK_src_vf, &
                                                             ix, iy, iz)

        type(scalar_field), &
            dimension(sys_size), &
            intent(IN) :: qK_prim_vf, qK_cons_vf

        type(scalar_field), &
            dimension(sys_size), &
            intent(INOUT) :: FK_vf, FK_src_vf

        type(bounds_info), intent(IN) :: ix, iy, iz

        ! Partial densities, density, velocity, pressure, energy, advection
        ! variables, the specific heat ratio and liquid stiffness functions,
        ! the shear and volume Reynolds numbers and the Weber numbers
        real(kind(0d0)), dimension(cont_idx%end)          :: alpha_rho_K
        real(kind(0d0))                                   ::       rho_K
        real(kind(0d0)), dimension(num_dims)              ::       vel_K
        real(kind(0d0))                                   ::      pres_K
        real(kind(0d0))                                   ::         E_K
        real(kind(0d0)), dimension(adv_idx%end - E_idx)     ::       adv_K
        real(kind(0d0))                                   ::     gamma_K
        real(kind(0d0))                                   ::    pi_inf_K

        ! Generic loop iterators
        integer :: i, j, k, l

        dir_idx = (/1, 2, 3/); dir_flg = (/1d0, 0d0, 0d0/)

        ! Computing the flux variables from the primitive variables, without
        ! accounting for the contribution of either viscosity or capillarity
        do l = iz%beg, iz%end
            do k = iy%beg, iy%end
                do j = ix%beg, ix%end

                    do i = 1, cont_idx%end
                        alpha_rho_K(i) = qK_prim_vf(i)%sf(j, k, l)
                    end do

                    do i = 1, num_dims
                        vel_K(i) = qK_prim_vf(cont_idx%end + i)%sf(j, k, l)
                    end do

                    pres_K = qK_prim_vf(E_idx)%sf(j, k, l)

                    call s_convert_to_mixture_variables(qK_prim_vf, rho_K, &
                                                        gamma_K, pi_inf_K, &
                                                        j, k, l)

                    ! mass flux, \rho u
                    do i = 1, cont_idx%end
                        FK_vf(i)%sf(j, k, l) = alpha_rho_K(i)*vel_K(dir_idx(1))
                    end do

                    ! momentum flux, \rho u u + p I
                    do i = 1, num_dims
                        FK_vf(cont_idx%end + dir_idx(i))%sf(j, k, l) = &
                            rho_K*vel_K(dir_idx(1)) &
                            *vel_K(dir_idx(i)) &
                            + pres_K*dir_flg(dir_idx(i))
                    end do

                    ! energy flux, 0
                    FK_vf(E_idx)%sf(j, k, l) = 0.

                    ! vol. frac, nR, and nRdot fluxes, u{\alpha, nR, nRdot}
                    do i = adv_idx%beg, sys_size
                        FK_vf(i)%sf(j, k, l) = vel_K(dir_idx(1))*qK_cons_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !stop

    end subroutine s_convert_primitive_to_flux_variables_bubbles ! -----------------


    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_variables_conversion_module() ! ------------------

        ! Disassociating the pointer to the procedure that was utilized to
        ! to convert mixture or species variables to the mixture variables
        s_convert_to_mixture_variables => null()

    end subroutine s_finalize_variables_conversion_module ! ----------------

end module m_variables_conversion
