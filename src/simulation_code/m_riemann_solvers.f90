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

    type(bounds_info) :: is1, is2, is3
    type(bounds_info) :: ix, iy, iz


    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: qL_prim_rs_vf_flat, qR_prim_rs_vf_flat
    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: flux_vf_flat, flux_src_vf_flat

contains

    subroutine s_initialize_riemann_solvers_module() 

        integer :: i

        allocate (qL_prim_rs_vf(1:sys_size), qR_prim_rs_vf(1:sys_size))

        ! For dir=1
        ix%beg = -1; iy%beg =  0; iz%beg =  0
        ix%end = m; iy%end = n; iz%end = p

        is1 = ix; is2 = iy; is3 = iz

        ! Allocating Left, Right and Average Riemann Problem States
        do i = 1, sys_size
            allocate (qL_prim_rs_vf(i)%sf(is1%beg:is1%end, &
                                          is2%beg:is2%end, &
                                          is3%beg:is3%end))
            allocate (qR_prim_rs_vf(i)%sf(is1%beg + 1:is1%end + 1, &
                                          is2%beg:is2%end, &
                                          is3%beg:is3%end))
        end do

        allocate(qL_prim_rs_vf_flat(is1%beg  :is1%end  ,is2%beg:is2%end,is3%beg:is3%end,1:sys_size))
        allocate(qR_prim_rs_vf_flat(is1%beg+1:is1%end+1,is2%beg:is2%end,is3%beg:is3%end,1:sys_size))

        allocate(flux_vf_flat(is1%beg:is1%end,is2%beg:is2%end,is3%beg:is3%end,1:sys_size))
        allocate(flux_src_vf_flat(is1%beg:is1%end,is2%beg:is2%end,is3%beg:is3%end,1:sys_size))


    end subroutine s_initialize_riemann_solvers_module 


    subroutine s_hllc_riemann_solver(qL_prim_vf,  & 
                                     qR_prim_vf,  &
                                     flux_vf,     &
                                     flux_src_vf, &
                                     norm_dir)

        type(scalar_field), dimension(sys_size), intent(INOUT) :: qL_prim_vf, qR_prim_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: flux_vf, flux_src_vf
        integer, intent(IN) :: norm_dir

        real(kind(0d0)), dimension(10)  :: alpha_rho_L, alpha_rho_R
        real(kind(0d0))                 ::       rho_L, rho_R
        real(kind(0d0)), dimension(10)  ::       vel_L, vel_R
        real(kind(0d0))                 ::      pres_L, pres_R
        real(kind(0d0))                 ::         E_L, E_R
        real(kind(0d0))                 ::         H_L, H_R
        real(kind(0d0)), dimension(10)  ::     alpha_L, alpha_R
        real(kind(0d0))                 ::     gamma_L, gamma_R
        real(kind(0d0))                 ::    pi_inf_L, pi_inf_R
        real(kind(0d0))                 ::         c_L, c_R
        real(kind(0d0))                 :: rho_avg
        real(kind(0d0)), dimension(10)  :: vel_avg
        real(kind(0d0))                 :: H_avg
        real(kind(0d0))                 :: gamma_avg
        real(kind(0d0))                 :: c_avg
        real(kind(0d0))                 :: s_L, s_R, s_S
        real(kind(0d0))                 :: s_M, s_P
        real(kind(0d0))                 :: xi_M, xi_P
        real(kind(0d0))                 :: xi_L, xi_R
        real(kind(0d0)), dimension(10)  :: gammas, pi_infs

        integer :: ixb, ixe, iyb, iye, izb, ize
        integer :: cont_idx_e
        integer :: adv_idx_b, adv_idx_e

        integer :: i, j, k, l

        integer, dimension(3) :: dir_idx
        real(kind(0d0)), dimension(3) :: dir_flg

        cont_idx_e = cont_idx%end
        adv_idx_b = adv_idx%beg
        adv_idx_e = adv_idx%end

        ixb = is1%beg; ixe = is1%end
        iyb = is2%beg; iye = is2%end
        izb = is3%beg; ize = is3%end

        dir_idx(1) = 1
        dir_idx(2) = 2
        dir_idx(3) = 3

        dir_flg(1) = 1d0
        dir_flg(2) = 0
        dir_flg(3) = 0

        do i = 1,num_fluids
            gammas(i) = fluid_pp(i)%gamma
            pi_infs(i) = fluid_pp(i)%pi_inf
        end do

        ! do i = 1, sys_size
        !     qL_prim_rs_vf(i)%sf = qL_prim_vf(i)%sf(ix%beg:ix%end, &
        !                                            iy%beg:iy%end, &
        !                                            iz%beg:iz%end)
        !     qR_prim_rs_vf(i)%sf = qR_prim_vf(i)%sf(ix%beg + 1:ix%end + 1, &
        !                                            iy%beg:iy%end, &
        !                                            iz%beg:iz%end)
        ! end do

        do i = 1, sys_size
            qL_prim_rs_vf_flat(:,:,:,i) = qL_prim_vf(i)%sf(ix%beg:ix%end, &
                                                   iy%beg:iy%end, &
                                                   iz%beg:iz%end)
            qR_prim_rs_vf_flat(:,:,:,i) = qR_prim_vf(i)%sf(ix%beg + 1:ix%end + 1, &
                                                   iy%beg:iy%end, &
                                                   iz%beg:iz%end)
        end do

        !$ acc data copyin(qL_prim_rs_vf_flat,qR_prim_rs_vf_flat,dir_idx,dir_flg,gammas,pi_infs) copyout(flux_vf_flat,flux_src_vf_flat) 
        !$ acc parallel loop collapse(3) gang vector private(alpha_rho_L, alpha_rho_R, vel_L, vel_R, alpha_L, alpha_R)
        do l = izb, ize
            do k = iyb, iye
                do j = ixb, ixe

                    do i = 1, cont_idx_e
                        alpha_rho_L(i) = qL_prim_rs_vf_flat(j, k, l, i)
                        alpha_rho_R(i) = qR_prim_rs_vf_flat(j + 1, k, l, i)
                    end do

                    do i = 1, num_dims
                        vel_L(i) = qL_prim_rs_vf_flat(j, k, l, cont_idx_e + i)
                        vel_R(i) = qR_prim_rs_vf_flat(j + 1, k, l, cont_idx_e + i)
                    end do

                    do i = 1,num_fluids
                        alpha_L(i) = qL_prim_rs_vf_flat(j, k, l, E_idx + i)
                        alpha_R(i) = qR_prim_rs_vf_flat(j + 1, k, l, E_idx + i)
                    end do

                    pres_L = qL_prim_rs_vf_flat(j, k, l, E_idx)
                    pres_R = qR_prim_rs_vf_flat(j + 1, k, l, E_idx)

                    ! Compute average state
                    call s_convert_species_to_mixture_variables_acc_shb( &
                                                        alpha_rho_L, alpha_L, &
                                                        gammas, pi_infs, &
                                                        rho_L, &
                                                        gamma_L, &
                                                        pi_inf_L, num_fluids)


                    call s_convert_species_to_mixture_variables_acc_shb( &
                                                        alpha_rho_R, alpha_R, &
                                                        gammas, pi_infs, &
                                                        rho_R, &
                                                        gamma_R, &
                                                        pi_inf_R, num_fluids)

                    E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*sum(vel_L**2d0)
                    E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*sum(vel_R**2d0)

                    H_L = (E_L + pres_L)/rho_L
                    H_R = (E_R + pres_R)/rho_R
                    
                    ! Compute sound speeds
                    c_L = (H_L - 5d-1*sum(vel_L**2d0))/gamma_L
                    c_R = (H_R - 5d-1*sum(vel_R**2d0))/gamma_R

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

                    ! Arithmetic Average Riemann Problem State 
                    rho_avg = 5d-1*(rho_L + rho_R)
                    vel_avg = 5d-1*(vel_L + vel_R)
                    H_avg  = 5d-1*(H_L + H_R)
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

                    ! Compute wavespeeds
                    s_L = min(vel_L(dir_idx(1)) - c_L, vel_R(dir_idx(1)) - c_R)
                    s_R = max(vel_R(dir_idx(1)) + c_R, vel_L(dir_idx(1)) + c_L)

                    s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx(1))* &
                           (s_L - vel_L(dir_idx(1))) - &
                           rho_R*vel_R(dir_idx(1))* &
                           (s_R - vel_R(dir_idx(1)))) &
                          /(rho_L*(s_L - vel_L(dir_idx(1))) - &
                            rho_R*(s_R - vel_R(dir_idx(1))))

                    s_M = min(0d0, s_L)
                    s_P = max(0d0, s_R)

                    xi_L = (s_L - vel_L(dir_idx(1)))/(s_L - s_S)
                    xi_R = (s_R - vel_R(dir_idx(1)))/(s_R - s_S)

                    xi_M = (5d-1 + sign(5d-1, s_S))
                    xi_P = (5d-1 - sign(5d-1, s_S))

                    do i = 1, cont_idx_e
                        flux_vf_flat(j, k, l, i) = &
                            xi_M*alpha_rho_L(i) &
                            *(vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                            + xi_P*alpha_rho_R(i) &
                            *(vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                    end do

                    ! Momentum flux.
                    do i = 1, num_dims
                        flux_vf_flat(j, k, l, cont_idx_e + dir_idx(i)) = &
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
                    flux_vf_flat(j, k, l, E_idx) = &
                        xi_M*(vel_L(dir_idx(1))*(E_L + pres_L) + &
                             s_M*(xi_L*(E_L + (s_S - vel_L(dir_idx(1)))* &
                             (rho_L*s_S + pres_L/ &
                             (s_L - vel_L(dir_idx(1))))) - E_L)) &
                      + xi_P*(vel_R(dir_idx(1))*(E_R + pres_R) + &
                             s_P*(xi_R*(E_R + (s_S - vel_R(dir_idx(1)))* &
                             (rho_R*s_S + pres_R/ &
                             (s_R - vel_R(dir_idx(1))))) - E_R))

                    ! Volume fraction flux
                    do i = adv_idx_b, adv_idx_e
                        flux_vf_flat(j, k, l, i) = &
                            xi_M*qL_prim_rs_vf_flat(j, k, l, i) &
                            *(vel_L(dir_idx(1)) + s_M*(xi_L - 1d0)) &
                            + xi_P*qR_prim_rs_vf_flat(j + 1, k, l, i) &
                            *(vel_R(dir_idx(1)) + s_P*(xi_R - 1d0))
                    end do

                    ! Source for volume fraction advection equation
                    do i = 1, num_dims
                        ! This only works in 1D
                        flux_src_vf_flat(j, k, l, adv_idx_b) = &
                            xi_M*(vel_L(dir_idx(i)) + &
                                  dir_flg(dir_idx(i))* &
                                  s_M*(xi_L - 1d0)) &
                          + xi_P*(vel_R(dir_idx(i)) + &
                                  dir_flg(dir_idx(i))* &
                                  s_P*(xi_R - 1d0))
                    end do
                end do
            end do
        end do
        !$ acc end parallel loop 
        !$ acc end data

        ! do i = 1,1 !sys_size
        !     do j = ixb, ixe
        !         print*, 'rank,flux:',proc_rank,j,flux_vf_flat(j,0,0,i)
        !     end do
        ! end do

        ! do j = ixb, ixe
        !     print*, 'rs:',j,flux_src_vf_flat(j,0,0,adv_idx_b)
        ! end do

    end subroutine s_hllc_riemann_solver 


    subroutine s_convert_species_to_mixture_variables_acc_shb(alpha_rho_K, &
                                                      alpha_K, &
                                                      gammas, pi_infs, &
                                                      rho_K, &
                                                      gamma_K, pi_inf_K, &
                                                      num_fluids &
                                                      )
        !$ acc routine seq

        real(kind(0d0)), dimension(num_fluids), intent(IN) :: alpha_rho_K, alpha_K, gammas, pi_infs
        real(kind(0d0)), intent(OUT) :: rho_K, gamma_K, pi_inf_K
        integer :: i
        integer, intent(in) :: num_fluids

        rho_K = 0d0; gamma_K = 0d0; pi_inf_K = 0d0

        do i = 1, num_fluids
            rho_K = rho_K + alpha_rho_K(i)
            gamma_K = gamma_K + alpha_K(i)*gammas(i)
            pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
        end do

    end subroutine s_convert_species_to_mixture_variables_acc_shb


    subroutine s_finalize_riemann_solvers_module() 

        integer :: i 

        do i = 1, sys_size
            deallocate (qL_prim_rs_vf(i)%sf, qR_prim_rs_vf(i)%sf)
        end do

        deallocate (qL_prim_rs_vf, qR_prim_rs_vf)

        deallocate (qL_prim_rs_vf_flat)
        deallocate (qR_prim_rs_vf_flat)

        deallocate (flux_vf_flat)
        deallocate (flux_src_vf_flat)

    end subroutine s_finalize_riemann_solvers_module 

end module m_riemann_solvers
