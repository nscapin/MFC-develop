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
        s_hllc_riemann_solver_bubbles, &
        s_finalize_riemann_solvers_module

    type(scalar_field), allocatable, dimension(:) :: qL_prim_rs_vf
    type(scalar_field), allocatable, dimension(:) :: qR_prim_rs_vf

    type(bounds_info) :: is1, is2, is3
    type(bounds_info) :: ix, iy, iz

    ! These are variables that are only on GPU
    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: qL_prim_rs_vf_flat, qR_prim_rs_vf_flat
    !$acc declare create (qL_prim_rs_vf_flat, qR_prim_rs_vf_flat)

    integer, dimension(3) :: dir_idx_acc
    real(kind(0d0)), dimension(3) :: dir_flg_acc
    !$acc declare create (dir_idx_acc, dir_flg_acc)

    integer, dimension(:), allocatable :: bub_idx_rrs, bub_idx_rvs
    !$acc declare create(bub_idx_rrs, bub_idx_rvs)

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

        dir_idx_acc(1) = 1
        dir_idx_acc(2) = 2
        dir_idx_acc(3) = 3

        dir_flg_acc(1) = 1d0
        dir_flg_acc(2) = 0
        dir_flg_acc(3) = 0
        !$acc update device(dir_flg_acc,dir_idx_acc)

        allocate(bub_idx_rrs(1:nb))
        allocate(bub_idx_rvs(1:nb))

        do i = 1,nb
            bub_idx_rrs(i) = bub_idx%rs(i)
            bub_idx_rvs(i) = bub_idx%vs(i)
        end do
        !$acc update device(bub_idx_rrs, bub_idx_rvs)

    end subroutine s_initialize_riemann_solvers_module 


    subroutine s_hllc_riemann_solver(qL_prim_vf,  & 
                                     qR_prim_vf,  &
                                     flux_vf_flat,     &
                                     flux_src_vf_flat  &
                                     )

        !! This -4 is because you can't pass arrays with negative indexing consistently
        !! and it only applys to WENO5 (buffsize = 4), also the other directions start/end at 0
        real(kind(0d0)), dimension(-4:,0:,0:,1:), intent(INOUT) :: qL_prim_vf, qR_prim_vf
        real(kind(0d0)), dimension(-1:,0:,0:,1:), intent(INOUT) :: flux_vf_flat, flux_src_vf_flat

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

        integer :: ixb, ixe, iyb, iye, izb, ize
        integer :: cont_idx_e
        integer :: adv_idx_b, adv_idx_e

        integer :: i, j, k, l

        cont_idx_e = cont_idx%end
        adv_idx_b = adv_idx%beg
        adv_idx_e = adv_idx%end

        ixb = is1%beg; ixe = is1%end
        iyb = is2%beg; iye = is2%end
        izb = is3%beg; ize = is3%end


        !$acc data present(qL_prim_rs_vf_flat, qR_prim_rs_vf_flat, qL_prim_vf, qR_prim_vf, flux_vf_flat, flux_src_vf_flat) copyout(gammas,pi_infs)

        !$acc parallel loop collapse(4) gang vector
        do l = izb, ize
            do k = iyb, iye
                do j = ixb, ixe
                    do i = 1,sys_size
                        qL_prim_rs_vf_flat(j, k, l, i) = qL_prim_vf(j, k, l, i)
                        qR_prim_rs_vf_flat(j + 1, k, l, i) = qR_prim_vf(j + 1, k, l, i)
                    end do
                end do
            end do
        end do
        !$acc end parallel loop

        ! do i = 1,1 !sys_size
        !     do j = ixb, ixe
        !         print*, 'rank,re:',j,qL_prim_rs_vf_flat(j,0,0,i),qR_prim_rs_vf_flat(j+1,0,0,i)
        !     end do
        ! end do

        !$acc parallel loop collapse(3) gang vector private(alpha_rho_L, alpha_rho_R, vel_L, vel_R, alpha_L, alpha_R, vel_avg)
        do l = izb, ize
            do k = iyb, iye
                do j = ixb, ixe

                    dir_idx_acc(1) = 1
                    dir_idx_acc(2) = 2
                    dir_idx_acc(3) = 3

                    dir_flg_acc(1) = 1d0
                    dir_flg_acc(2) = 0
                    dir_flg_acc(3) = 0

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
                                                        rho_L, &
                                                        gamma_L, &
                                                        pi_inf_L, num_fluids)

                    call s_convert_species_to_mixture_variables_acc_shb( &
                                                        alpha_rho_R, alpha_R, &
                                                        rho_R, &
                                                        gamma_R, &
                                                        pi_inf_R, num_fluids)

                    ! Should be sum(vel) in 2/3D but allocaiton issues..
                    E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*vel_L(1)**2d0
                    E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*vel_R(1)**2d0

                    H_L = (E_L + pres_L)/rho_L
                    H_R = (E_R + pres_R)/rho_R
                    
                    ! Compute sound speeds
                    ! Should be sum(vel) in 2/3D but allocaiton issues on GPU..
                    c_L = (H_L - 5d-1*(vel_L(1)**2d0))/gamma_L
                    c_R = (H_R - 5d-1*(vel_R(1)**2d0))/gamma_R

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
                        ! Should be sum(vel) in 2/3D but allocaiton issues..
                        if ((H_avg - 5d-1*(vel_avg(1)**2d0)) < 0d0) then
                            c_avg = sgm_eps
                        else
                            c_avg = sqrt((H_avg - 5d-1*(vel_avg(1)**2d0))/gamma_avg)
                        end if
                    else
                        c_avg = sqrt((H_avg - 5d-1*(vel_avg(1)**2d0))/gamma_avg)
                    end if

                    ! Compute wavespeeds
                    s_L = min(vel_L(dir_idx_acc(1)) - c_L, vel_R(dir_idx_acc(1)) - c_R)
                    s_R = max(vel_R(dir_idx_acc(1)) + c_R, vel_L(dir_idx_acc(1)) + c_L)

                    s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx_acc(1))* &
                           (s_L - vel_L(dir_idx_acc(1))) - &
                           rho_R*vel_R(dir_idx_acc(1))* &
                           (s_R - vel_R(dir_idx_acc(1)))) &
                          /(rho_L*(s_L - vel_L(dir_idx_acc(1))) - &
                            rho_R*(s_R - vel_R(dir_idx_acc(1))))

                    s_M = min(0d0, s_L)
                    s_P = max(0d0, s_R)

                    xi_L = (s_L - vel_L(dir_idx_acc(1)))/(s_L - s_S)
                    xi_R = (s_R - vel_R(dir_idx_acc(1)))/(s_R - s_S)

                    xi_M = (5d-1 + sign(5d-1, s_S))
                    xi_P = (5d-1 - sign(5d-1, s_S))

                    do i = 1, cont_idx_e
                        flux_vf_flat(j, k, l, i) = &
                            xi_M*alpha_rho_L(i) &
                            *(vel_L(dir_idx_acc(1)) + s_M*(xi_L - 1d0)) &
                            + xi_P*alpha_rho_R(i) &
                            *(vel_R(dir_idx_acc(1)) + s_P*(xi_R - 1d0))
                    end do

                    ! Momentum flux.
                    do i = 1, num_dims
                        flux_vf_flat(j, k, l, cont_idx_e + dir_idx_acc(i)) = &
                            xi_M*(rho_L*(vel_L(dir_idx_acc(1))* &
                                  vel_L(dir_idx_acc(i)) + &
                                  s_M*(xi_L*(dir_flg_acc(dir_idx_acc(i))*s_S + &
                                  (1d0 - dir_flg_acc(dir_idx_acc(i)))* &
                                  vel_L(dir_idx_acc(i))) - vel_L(dir_idx_acc(i)))) + &
                                  dir_flg_acc(dir_idx_acc(i))*(pres_L)) &
                          + xi_P*(rho_R*(vel_R(dir_idx_acc(1))* &
                                  vel_R(dir_idx_acc(i)) + &
                                  s_P*(xi_R*(dir_flg_acc(dir_idx_acc(i))*s_S + &
                                  (1d0 - dir_flg_acc(dir_idx_acc(i)))* &
                                  vel_R(dir_idx_acc(i))) - vel_R(dir_idx_acc(i)))) + &
                                  dir_flg_acc(dir_idx_acc(i))*(pres_R))
                    end do


                    flux_vf_flat(j, k, l, E_idx) = &
                        xi_M*(vel_L(dir_idx_acc(1))*(E_L + pres_L) + &
                             s_M*(xi_L*(E_L + (s_S - vel_L(dir_idx_acc(1)))* &
                             (rho_L*s_S + pres_L/ &
                             (s_L - vel_L(dir_idx_acc(1))))) - E_L)) &
                      + xi_P*(vel_R(dir_idx_acc(1))*(E_R + pres_R) + &
                             s_P*(xi_R*(E_R + (s_S - vel_R(dir_idx_acc(1)))* &
                             (rho_R*s_S + pres_R/ &
                             (s_R - vel_R(dir_idx_acc(1))))) - E_R))

                    ! Volume fraction flux
                    do i = adv_idx_b, adv_idx_e
                        flux_vf_flat(j, k, l, i) = &
                            xi_M*qL_prim_rs_vf_flat(j, k, l, i) &
                            *(vel_L(dir_idx_acc(1)) + s_M*(xi_L - 1d0)) &
                            + xi_P*qR_prim_rs_vf_flat(j + 1, k, l, i) &
                            *(vel_R(dir_idx_acc(1)) + s_P*(xi_R - 1d0))
                    end do

                    ! Source for volume fraction advection equation
                    do i = 1, num_dims
                        ! This only works in 1D
                        flux_src_vf_flat(j, k, l, adv_idx_b) = &
                            xi_M*(vel_L(dir_idx_acc(i)) + &
                                  dir_flg_acc(dir_idx_acc(i))* &
                                  s_M*(xi_L - 1d0)) &
                          + xi_P*(vel_R(dir_idx_acc(i)) + &
                                  dir_flg_acc(dir_idx_acc(i))* &
                                  s_P*(xi_R - 1d0))
                    end do
                end do
            end do
        end do
        !$acc end parallel loop 
        !$acc end data

        ! do i = 1,1 !sys_size
        !     do j = ixb, ixe
        !         print*, 'rank,flux:',j,flux_vf_flat(j,0,0,i)
        !     end do
        ! end do

        ! do j = ixb, ixe
        !     print*, 'rs:',j,flux_src_vf_flat(j,0,0,adv_idx_b)
        ! end do

    end subroutine s_hllc_riemann_solver 


    subroutine s_hllc_riemann_solver_bubbles(qL_prim_vf,  & 
                                     qR_prim_vf,  &
                                     flux_vf_flat,     &
                                     flux_src_vf_flat  &
                                     )

        real(kind(0d0)), dimension(-4:,0:,0:,1:), intent(INOUT) :: qL_prim_vf, qR_prim_vf
        real(kind(0d0)), dimension(-1:,0:,0:,1:), intent(INOUT) :: flux_vf_flat, flux_src_vf_flat

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

        real(kind(0d0))                 ::       nbub_L,     nbub_R
        real(kind(0d0))                 ::     ptilde_L,   ptilde_R
        real(kind(0d0)), dimension(11) ::         R0_L,       R0_R
        real(kind(0d0)), dimension(11) ::         V0_L,       V0_R
        real(kind(0d0)), dimension(11) ::        pbw_L,      pbw_R
        real(kind(0d0)), dimension(11) ::       q_temp

        real(kind(0d0)) :: PbwR3Lbar, Pbwr3Rbar
        real(kind(0d0)) :: R3Lbar, R3Rbar
        real(kind(0d0)) :: R3V2Lbar, R3V2Rbar

        integer :: ixb, ixe, iyb, iye, izb, ize
        integer :: cont_idx_e
        integer :: adv_idx_b, adv_idx_e
        integer :: bub_idx_b, bub_idx_e

        integer :: i, j, k, l

        cont_idx_e = cont_idx%end
        adv_idx_b = adv_idx%beg
        adv_idx_e = adv_idx%end
        bub_idx_b = bub_idx%beg
        bub_idx_e = bub_idx%end

        ixb = is1%beg; ixe = is1%end
        iyb = is2%beg; iye = is2%end
        izb = is3%beg; ize = is3%end

        !$acc data present(qL_prim_rs_vf_flat, qR_prim_rs_vf_flat, qL_prim_vf, qR_prim_vf, flux_vf_flat, flux_src_vf_flat, R0, nb)

        !$acc parallel loop collapse(4) gang vector
        do l = izb, ize
            do k = iyb, iye
                do j = ixb, ixe
                    do i = 1,sys_size
                        qL_prim_rs_vf_flat(j, k, l, i) = qL_prim_vf(j, k, l, i)
                        qR_prim_rs_vf_flat(j + 1, k, l, i) = qR_prim_vf(j + 1, k, l, i)
                    end do
                end do
            end do
        end do
        !$acc end parallel loop

        !$acc parallel loop collapse(3) gang vector private(alpha_rho_L, alpha_rho_R, vel_L, vel_R, alpha_L, alpha_R, vel_avg, R0_L, R0_R, V0_L, V0_R, pbw_L, pbw_R, q_temp)
        do l = izb, ize
            do k = iyb, iye
                do j = ixb, ixe

                    dir_idx_acc(1) = 1
                    dir_idx_acc(2) = 2
                    dir_idx_acc(3) = 3

                    dir_flg_acc(1) = 1d0
                    dir_flg_acc(2) = 0
                    dir_flg_acc(3) = 0

                    do i = 1, cont_idx_e
                        alpha_rho_L(i) = qL_prim_rs_vf_flat(j, k, l, i)
                        alpha_rho_R(i) = qR_prim_rs_vf_flat(j + 1, k, l, i)
                    end do

                    do i = 1, num_dims
                        vel_L(i) = qL_prim_rs_vf_flat(j, k, l, cont_idx_e + i)
                        vel_R(i) = qR_prim_rs_vf_flat(j + 1, k, l, cont_idx_e + i)
                    end do

                    ! Compute mixture average state
                    rho_L = alpha_rho_L(1)
                    gamma_L = gammas(1)
                    pi_inf_L = pi_infs(1)

                    rho_R =  alpha_rho_R(1)
                    gamma_R = gammas(1)
                    pi_inf_R = pi_infs(1)

                    pres_L = qL_prim_rs_vf_flat(j, k, l, E_idx)
                    pres_R = qR_prim_rs_vf_flat(j + 1, k, l, E_idx)

                    ! Should be sum(vel) in 2/3D but allocaiton issues..
                    E_L = gamma_L*pres_L + pi_inf_L + 5d-1*rho_L*vel_L(1)**2d0
                    E_R = gamma_R*pres_R + pi_inf_R + 5d-1*rho_R*vel_R(1)**2d0

                    H_L = (E_L + pres_L)/rho_L
                    H_R = (E_R + pres_R)/rho_R

                    ! bubbles
                    do i = 1,num_fluids
                        alpha_L(i) = qL_prim_rs_vf_flat(j, k, l, E_idx + i)
                        alpha_R(i) = qR_prim_rs_vf_flat(j + 1, k, l, E_idx + i)
                    end do

                    do i = 1,nb
                        R0_L(i) = qL_prim_rs_vf_flat(j  ,k,l,bub_idx_rrs(i))
                        R0_R(i) = qR_prim_rs_vf_flat(j+1,k,l,bub_idx_rrs(i))
                        V0_L(i) = qL_prim_rs_vf_flat(j  ,k,l,bub_idx_rvs(i))
                        V0_R(i) = qR_prim_rs_vf_flat(j+1,k,l,bub_idx_rvs(i))
                    end do 
                
                    call s_comp_n_from_prim(alpha_L(num_fluids),R0_L(1:nb),nbub_L)
                    call s_comp_n_from_prim(alpha_R(num_fluids),R0_R(1:nb),nbub_R)
                    
                    DO i = 1,nb
                        pbw_L(i) = f_cpbw_KM(R0(i),R0_L(i),V0_L(i),0d0)
                        pbw_R(i) = f_cpbw_KM(R0(i),R0_R(i),V0_R(i),0d0)
                    END DO

                    q_temp(1:nb) = pbw_L(1:nb)*(R0_L(1:nb)**3d0)
                    CALL s_quad_acc(q_temp(1:nb), PbwR3Lbar, nb)
                    q_temp(1:nb) = pbw_R(1:nb)*(R0_R(1:nb)**3d0)
                    CALL s_quad_acc(q_temp(1:nb), PbwR3Rbar, nb)

                    q_temp(1:nb) = R0_L(1:nb)**3.d0
                    CALL s_quad_acc(q_temp(1:nb), R3Lbar, nb)
                    q_temp(1:nb) = R0_R(1:nb)**3.d0
                    CALL s_quad_acc(q_temp(1:nb), R3Rbar, nb)
                
                    q_temp(1:nb) = (R0_L(1:nb)**3.d0)*(V0_L(1:nb)**2.d0)
                    CALL s_quad_acc(q_temp(1:nb), R3V2Lbar, nb)
                    q_temp(1:nb) = (R0_R(1:nb)**3.d0)*(V0_R(1:nb)**2.d0)
                    CALL s_quad_acc(q_temp(1:nb), R3V2Rbar, nb)

                    ptilde_L = alpha_L(num_fluids)*(pres_L - PbwR3Lbar/R3Lbar - & 
                        rho_L*R3V2Lbar/R3Lbar )
                    ptilde_R = alpha_R(num_fluids)*(pres_R - PbwR3Rbar/R3Rbar - & 
                        rho_R*R3V2Rbar/R3Rbar )
                    
                    ! Compute sound speeds
                    c_L =   & 
                            (1d0/gamma_L + 1d0) *   &
                            (pres_L + pi_inf_L) /   &
                            (rho_L*(1d0-alpha_L(num_fluids))) 
                    c_R =   & 
                            (1d0/gamma_R + 1d0) *   &
                            (pres_R + pi_inf_R) /   &
                            (rho_R*(1d0-alpha_R(num_fluids))) 

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
                        ! Should be sum(vel) in 2/3D but allocaiton issues..
                        if ((H_avg - 5d-1*(vel_avg(1)**2d0)) < 0d0) then
                            c_avg = sgm_eps
                        else
                            c_avg = sqrt((H_avg - 5d-1*(vel_avg(1)**2d0))/gamma_avg)
                        end if
                    else
                        c_avg = sqrt((H_avg - 5d-1*(vel_avg(1)**2d0))/gamma_avg)
                    end if

                    ! Compute wavespeeds
                    s_L = min(vel_L(dir_idx_acc(1)) - c_L, vel_R(dir_idx_acc(1)) - c_R)
                    s_R = max(vel_R(dir_idx_acc(1)) + c_R, vel_L(dir_idx_acc(1)) + c_L)

                    s_S = (pres_R - pres_L + rho_L*vel_L(dir_idx_acc(1))* &
                           (s_L - vel_L(dir_idx_acc(1))) - &
                           rho_R*vel_R(dir_idx_acc(1))* &
                           (s_R - vel_R(dir_idx_acc(1)))) &
                          /(rho_L*(s_L - vel_L(dir_idx_acc(1))) - &
                            rho_R*(s_R - vel_R(dir_idx_acc(1))))

                    s_M = min(0d0, s_L)
                    s_P = max(0d0, s_R)

                    xi_L = (s_L - vel_L(dir_idx_acc(1)))/(s_L - s_S)
                    xi_R = (s_R - vel_R(dir_idx_acc(1)))/(s_R - s_S)

                    xi_M = (5d-1 + sign(5d-1, s_S))
                    xi_P = (5d-1 - sign(5d-1, s_S))

                    do i = 1, cont_idx_e
                        flux_vf_flat(j, k, l, i) = &
                            xi_M*alpha_rho_L(i) &
                            *(vel_L(dir_idx_acc(1)) + s_M*(xi_L - 1d0)) &
                            + xi_P*alpha_rho_R(i) &
                            *(vel_R(dir_idx_acc(1)) + s_P*(xi_R - 1d0))
                    end do

                    ! Momentum flux.
                    do i = 1, num_dims
                        flux_vf_flat(j, k, l, cont_idx_e + dir_idx_acc(i)) = &
                            xi_M*(rho_L*(vel_L(dir_idx_acc(1))* &
                                  vel_L(dir_idx_acc(i)) + &
                                  s_M*(xi_L*(dir_flg_acc(dir_idx_acc(i))*s_S + &
                                  (1d0 - dir_flg_acc(dir_idx_acc(i)))* &
                                  vel_L(dir_idx_acc(i))) - vel_L(dir_idx_acc(i)))) + &
                                  dir_flg_acc(dir_idx_acc(i))*(pres_L-ptilde_L)) &
                          + xi_P*(rho_R*(vel_R(dir_idx_acc(1))* &
                                  vel_R(dir_idx_acc(i)) + &
                                  s_P*(xi_R*(dir_flg_acc(dir_idx_acc(i))*s_S + &
                                  (1d0 - dir_flg_acc(dir_idx_acc(i)))* &
                                  vel_R(dir_idx_acc(i))) - vel_R(dir_idx_acc(i)))) + &
                                  dir_flg_acc(dir_idx_acc(i))*(pres_R-ptilde_R))
                    end do

                    flux_vf_flat(j, k, l, E_idx) = &
                        xi_M*(vel_L(dir_idx_acc(1))*(E_L + pres_L) + &
                             s_M*(xi_L*(E_L + (s_S - vel_L(dir_idx_acc(1)))* &
                             (rho_L*s_S + pres_L/ &
                             (s_L - vel_L(dir_idx_acc(1))))) - E_L)) &
                      + xi_P*(vel_R(dir_idx_acc(1))*(E_R + pres_R) + &
                             s_P*(xi_R*(E_R + (s_S - vel_R(dir_idx_acc(1)))* &
                             (rho_R*s_S + pres_R/ &
                             (s_R - vel_R(dir_idx_acc(1))))) - E_R))

                    ! Volume fraction flux
                    do i = adv_idx_b, adv_idx_e
                        flux_vf_flat(j, k, l, i) = &
                            xi_M*qL_prim_rs_vf_flat(j, k, l, i) &
                            *(vel_L(dir_idx_acc(1)) + s_M*(xi_L - 1d0)) &
                            + xi_P*qR_prim_rs_vf_flat(j + 1, k, l, i) &
                            *(vel_R(dir_idx_acc(1)) + s_P*(xi_R - 1d0))
                    end do

                    ! Bubbles
                    do i = bub_idx_b,bub_idx_e
                        flux_vf_flat(j,k,l,i) =   &
                            xi_M*nbub_L*qL_prim_rs_vf_flat(j,k,l,i) &
                            * (vel_L(dir_idx_acc(1)) + s_M*(xi_L - 1d0)) &
                            + xi_P*nbub_R*qR_prim_rs_vf_flat(j+1,k,l,i)  &
                            * (vel_R(dir_idx_acc(1)) + s_P*(xi_R - 1d0))
                    end do

                    ! Source for volume fraction advection equation
                    do i = 1, num_dims
                        ! This only works in 1D
                        flux_src_vf_flat(j, k, l, adv_idx_b) = &
                            xi_M*(vel_L(dir_idx_acc(i)) + &
                                  dir_flg_acc(dir_idx_acc(i))* &
                                  s_M*(xi_L - 1d0)) &
                          + xi_P*(vel_R(dir_idx_acc(i)) + &
                                  dir_flg_acc(dir_idx_acc(i))* &
                                  s_P*(xi_R - 1d0))
                    end do
                end do
            end do
        end do
        !$acc end parallel loop 
        !$acc end data

        ! print*, 'done with Riemann'
        ! call s_mpi_abort()


    end subroutine s_hllc_riemann_solver_bubbles


    subroutine s_convert_species_to_mixture_variables_acc_shb(alpha_rho_K, &
                                                      alpha_K, &
                                                      rho_K, &
                                                      gamma_K, pi_inf_K, &
                                                      num_fluids &
                                                      )
        !$acc routine seq

        real(kind(0d0)), dimension(num_fluids), intent(IN) :: alpha_rho_K, alpha_K
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


    subroutine s_quad_acc(func, mom, nb)
        !$acc routine seq

        real(kind(0d0)), dimension(nb), intent(in) :: func
        real(kind(0d0)), intent(out) :: mom
        integer, intent(in) :: nb
        integer :: i

        mom = 0d0
        do i = 1,nb
            mom = mom + weight(i)*func(i)
        end do

    end subroutine s_quad_acc


    subroutine s_finalize_riemann_solvers_module() 

        integer :: i 

        do i = 1, sys_size
            deallocate (qL_prim_rs_vf(i)%sf, qR_prim_rs_vf(i)%sf)
        end do

        deallocate (qL_prim_rs_vf, qR_prim_rs_vf)

        deallocate (qL_prim_rs_vf_flat)
        deallocate (qR_prim_rs_vf_flat)

        ! deallocate (flux_vf_flat)
        ! deallocate (flux_src_vf_flat)

    end subroutine s_finalize_riemann_solvers_module 

end module m_riemann_solvers
