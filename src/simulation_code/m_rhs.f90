module m_rhs

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_variables_conversion !< State variables type conversion procedures

    use m_weno                 !< Weighted and essentially non-oscillatory (WENO)
                               !! schemes for spatial reconstruction of variables

    use m_riemann_solvers      !< Exact and approximate Riemann problem solvers

    use m_cbc                  !< Characteristic boundary conditions (CBC)

    use m_bubbles              !< Bubble dynamic routines

    use m_qbmm                 !< Moment inversion

    use nvtx
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_rhs_module, &
         s_alt_rhs, &
         s_populate_variables_buffers, &
         s_finalize_rhs_module

    type(vector_field) :: q_cons_qp
    type(vector_field) :: q_prim_qp

    type(bounds_info) :: iv !< Vector field indical bounds

    !> @name Indical bounds in the x-, y- and z-directions
    !> @{
    type(bounds_info) :: ix, iy, iz
    !> @}


    character(50) :: file_path !< Local file path for saving debug files

    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: qK_cons_vf_flat
    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: qK_prim_vf_flat

    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: rhs_vf_flat

    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: vL_vf_flat
    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: vR_vf_flat

    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: flux_vf_flat
    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: flux_src_vf_flat

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_rhs_module() ! ---------------------------------

        type(bounds_info) :: is1, is2, is3
        integer :: i, j, k !< Generic loop iterators


        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -buff_size; iy%beg = 0; iz%beg = 0

        if (n > 0) iy%beg = -buff_size; if (p > 0) iz%beg = -buff_size

        ix%end = m - ix%beg; iy%end = n - iy%beg; iz%end = p - iz%beg
        ! ==================================================================


        allocate (q_cons_qp%vf(1:sys_size))
        allocate (q_prim_qp%vf(1:sys_size))

        ! ==================================================================

        ! Configuring Coordinate Direction Indexes =========================
        ix%beg = -1; if (n > 0) iy%beg = -1; if (p > 0) iz%beg = -1

        ix%end = m; iy%end = n; iz%end = p
        ! ==================================================================


        is1%beg = -buff_size
        is1%end = m - is1%beg

        is2%beg = 0
        is2%end = n

        is3%beg = 0
        is3%end = p

        ! Var conv input
        allocate( qK_cons_vf_flat(is1%beg:is1%end,is2%beg:is2%end,is3%beg:is3%end,1:sys_size ) )

        ! Var conv output, WENO input
        allocate( qK_prim_vf_flat(is1%beg:is1%end,is2%beg:is2%end,is3%beg:is3%end,1:sys_size ) )

        ! WENO output
        allocate( vL_vf_flat(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size) )
        allocate( vR_vf_flat(is1%beg:is1%end, is2%beg:is2%end, is3%beg:is3%end, 1:sys_size) )

        ! For Riemann solver
        is1%beg = -1; is2%beg =  0; is3%beg =  0
        is1%end = m; is2%end = n; is3%end = p

        ! Riemann solver output
        allocate(flux_vf_flat(is1%beg:is1%end,is2%beg:is2%end,is3%beg:is3%end,1:sys_size))
        allocate(flux_src_vf_flat(is1%beg:is1%end,is2%beg:is2%end,is3%beg:is3%end,1:sys_size))

        ! Compute rhs
        allocate(rhs_vf_flat(is1%beg:is1%end,is2%beg:is2%end,is3%beg:is3%end,1:sys_size)  )

    end subroutine s_initialize_rhs_module ! -------------------------------


    subroutine s_alt_rhs(q_cons_vf, q_prim_vf, rhs_vf, t_step)

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(INOUT) :: rhs_vf
        integer, intent(IN) :: t_step

        integer :: i, j, k, it_t

        type(bounds_info) :: is1_weno, is2_weno, is3_weno

        integer :: adv_idx_b, adv_idx_e
        real(kind(0d0)) :: start_time, end_time

        ix%beg = -buff_size; ix%end = m - ix%beg; 
        iv%beg = 1; iv%end = adv_idx%end

        adv_idx_b = adv_idx%beg
        adv_idx_e = adv_idx%end

        do i = 1, sys_size
            q_cons_qp%vf(i)%sf => q_cons_vf(i)%sf
            q_prim_qp%vf(i)%sf => q_prim_vf(i)%sf
        end do

        ! call nvtxStartRange("RHS-Pop. var. buffers")
        ! call s_populate_conservative_variables_buffers()
        ! call nvtxEndRange

        ! Flatten
        call nvtxStartRange("RHS-Flatten input")
        do i = 1, sys_size
            do j = ix%beg,ix%end
                qK_cons_vf_flat(j,0,0,i) = q_cons_vf(i)%sf(j,0,0)
            end do
        end do
        call nvtxEndRange

        ! if (proc_rank == 1) then
        !     do k = 1,1
        !         do j = 0,m
        !             print*, 'Rank,Cons: ', &
        !                 proc_rank, j, qK_cons_vf_flat(j,0,0,k) 
        !         end do
        !     end do
        ! end if

        ! print*, 'after flatten proc rank', proc_rank

        !$acc data copyin(qK_cons_vf_flat) copyout(rhs_vf_flat,flux_vf_flat, vL_vf_flat, vR_vf_flat) create(qK_prim_vf_flat, flux_src_vf_flat)
        do it_t = 1,t_step_stop
            call nvtxStartRange("Time step")

            start_time = mpi_wtime()

            call nvtxStartRange("RHS-Pop. var. buffers")
            call s_populate_conservative_variables_buffers()
            call nvtxEndRange

            if (proc_rank == num_procs - 1) then
                do k = 1,sys_size
                    do j = ix%beg,ix%end
                        print*, 'Rank,Cons: ', &
                            proc_rank, j, k, qK_cons_vf_flat(j,0,0,k) 
                    end do
                end do
            end if

            ! print*, 'after cons var buff proc rank', proc_rank

            call nvtxStartRange("RHS-Convert to prim")
            call s_convert_conservative_to_primitive_variables_acc( &
                qK_cons_vf_flat, qK_prim_vf_flat, &
                ix, iy, iz)
            call nvtxEndRange

            ! if (proc_rank == num_procs - 1) then
            !     do k = 1,sys_size
            !         do j = 0,m
            !             print*, 'Rank,Prim: ', &
            !                 proc_rank, j, k, qK_prim_vf_flat(j,0,0,k)
            !         end do
            !     end do
            ! end if
            ! call s_mpi_abort()

            call nvtxStartRange("RHS-WENO")
            is1_weno = ix; is2_weno = iy; is3_weno = iz
            is1_weno%beg = is1_weno%beg + weno_polyn
            is1_weno%end = is1_weno%end - weno_polyn

            call s_weno_alt( &
                        qK_prim_vf_flat, &
                        vL_vf_flat, vR_vf_flat, &
                        is1_weno, is2_weno, is3_weno)
            call nvtxEndRange




            ! do k = 1,sys_size-1
            !     print*, 'Variable ', k 
            !     do j = 0,m
            !         print*, 'Prim, L, R: ', &
            !             qK_prim_vf_flat(j,0,0,k), &
            !             vL_vf_flat(j,0,0,k),  &
            !             vR_vf_flat(j,0,0,k)
            !     end do
            ! end do

            call nvtxStartRange("RHS-Riemann")
            call s_hllc_riemann_solver( &
                                  vR_vf_flat, vL_vf_flat, &
                                  flux_vf_flat, &
                                  flux_src_vf_flat )
            call nvtxEndRange

            ! !$acc data copyout(flux_vf_flat)
            ! if (proc_rank == num_procs - 1) then
            !     do k = 1,sys_size
            !         do j = 0,m
            !             print*, 'Rank,flux: ', &
            !                 proc_rank, j, k, flux_vf_flat(j,0,0,k) 
            !         end do
            !     end do
            ! end if
            ! !$acc end data

            ! if (proc_rank == num_procs - 1) then
            !     do k = 1,sys_size
            !         do j = 0,m
            !             print*, 'Rank,flux: ', &
            !                 proc_rank, j, flux_vf_flat(j,0,0,k) 
            !         end do
            !     end do
            ! end if


            call nvtxStartRange("RHS-Diff fluxes")
            !$acc parallel loop collapse (2) gang vector
            do k = 0, m
                do j = 1, sys_size
                    rhs_vf_flat(k,0,0,j) = 1d0/dx(k)* &
                        (flux_vf_flat(k-1,0,0,j) &
                       - flux_vf_flat(k  ,0,0,j))
                end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

            ! if (proc_rank == num_procs - 1) then
            !     do k = 1,sys_size
            !         do j = 0,m
            !             print*, 'Rank,RHS1: ', &
            !                 proc_rank, j, rhs_vf_flat(j,0,0,k) 
            !         end do
            !     end do
            ! end if

            ! Apply source terms to RHS of advection equations
            call nvtxStartRange("RHS-Add srcs")
            !$acc parallel loop collapse(2) gang vector
            do k = 0, m
                do j = adv_idx_b, adv_idx_e
                    rhs_vf_flat(k,0,0,j) = &
                        rhs_vf_flat(k,0,0,j) + 1d0/dx(k) * &
                        qK_cons_vf_flat(k,0,0,j) * &
                         (flux_src_vf_flat(k  ,0,0,j) &
                        - flux_src_vf_flat(k-1,0,0,j))
                end do
            end do
            !$acc end parallel loop
            call nvtxEndRange
        
            call nvtxStartRange("RHS-Add RHS to Cons")
            !$acc parallel loop collapse (2) gang vector
            do k = 0, m
                do j = 1, sys_size
                    qK_cons_vf_flat(k,0,0,j) = &
                        qK_cons_vf_flat(k,0,0,j) + &
                        dt * rhs_vf_flat(k,0,0,j) 
                end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

            end_time = mpi_wtime()
            if (proc_rank == 0) then
                print*, 'RHS Eval Time [s]:', end_time - start_time
            end if
            call nvtxEndRange

            ! print*, 'end TS proc rank', proc_rank
        end do
        !$acc end data

        ! if (proc_rank == num_procs - 1) then
        !     do k = 1,sys_size
        !         do j = 0,m
        !             print*, 'Rank,Recon: ', &
        !                 proc_rank, j, k, vL_vf_flat(j,0,0,k), vR_vf_flat(j,0,0,k)
        !         end do
        !     end do
        ! end if

        !if (proc_rank == num_procs - 1) then
        !    do k = 3,3
        !        !sys_size
        !        do j = 0,m
        !            print*, 'Rank,flux: ', &
        !                proc_rank, j, k, flux_vf_flat(j,0,0,k) 
        !        end do
        !    end do
        !end if
        ! call s_mpi_abort()

        if (proc_rank == num_procs - 1) then
            do k = 1,sys_size
                do j = 0,m
                    print*, 'Rank,RHS: ', &
                        proc_rank, j, k, rhs_vf_flat(j,0,0,k) 
                end do
            end do
        end if

        ! if (proc_rank == 0) print*, 'out of TS loop in RHS'

        do i = 1, sys_size
            nullify (q_cons_qp%vf(i)%sf, q_prim_qp%vf(i)%sf)
        end do

        ! if (proc_rank == 0) print*, 'end rhs sub'
        

    end subroutine s_alt_rhs


    subroutine s_populate_conservative_variables_buffers() 
        ! This is called every RHS evaluation
        ! SHB: Suspect this needs attention for GPUs

        integer :: i, j, k

        print*, 'In pop cons buff: rank, bcx_b/e', proc_rank, bc_xb, bc_xe

        !$acc data present(qK_cons_vf_flat) 
        if (bc_xb == -1) then
            ! Periodic BC at beginning
            do i = 1, sys_size
                do j = 1, buff_size
                    qK_cons_vf_flat(-j, 0, 0, i) = &
                        qK_cons_vf_flat(m - (j - 1), 0, 0, i)
                end do
            end do
        else
            ! Processor BC at beginning
            call s_mpi_sendrecv_conservative_variables_buffers_acc( &
                qK_cons_vf_flat, -1)
        end if

        if (bc_xe == -1) then
            ! Periodic BC at end
            do i = 1, sys_size
                do j = 1, buff_size
                    qK_cons_vf_flat(m + j, 0, 0, i) = &
                        qK_cons_vf_flat(j - 1, 0, 0, i)
                end do
            end do
        else                            
            ! Processor BC at end
            call s_mpi_sendrecv_conservative_variables_buffers_acc( &
                qK_cons_vf_flat, 1)
        end if
        !$acc end data

        ! call s_mpi_abort()

    end subroutine s_populate_conservative_variables_buffers


    subroutine s_populate_variables_buffers(v_vf)
        ! This is only called from p_main before time stepping starts

        type(scalar_field), dimension(sys_size), intent(INOUT) :: v_vf

        integer :: i, j, k !< Generic loop iterators

        if (bc_x%beg <= -3) then         ! Ghost-cell extrap. BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(-j, 0:n, 0:p) = &
                        v_vf(i)%sf(0, 0:n, 0:p)
                end do
            end do

        elseif (bc_x%beg == -2) then     ! Symmetry BC at beginning

            do j = 1, buff_size

                do i = 1, cont_idx%end
                    v_vf(i)%sf(-j, 0:n, 0:p) = &
                        v_vf(i)%sf(j - 1, 0:n, 0:p)
                end do

                v_vf(mom_idx%beg)%sf(-j, 0:n, 0:p) = &
                    -v_vf(mom_idx%beg)%sf(j - 1, 0:n, 0:p)

                do i = mom_idx%beg + 1, sys_size
                    v_vf(i)%sf(-j, 0:n, 0:p) = &
                        v_vf(i)%sf(j - 1, 0:n, 0:p)
                end do

            end do

        elseif (bc_x%beg == -1) then     ! Periodic BC at beginning

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(-j, 0:n, 0:p) = &
                        v_vf(i)%sf(m - (j - 1), 0:n, 0:p)
                end do
            end do

        else                            ! Processor BC at beginning

            call s_mpi_sendrecv_conservative_variables_buffers( &
                v_vf, 1, -1)

        end if

        if (bc_x%end <= -3) then         ! Ghost-cell extrap. BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(m + j, 0:n, 0:p) = &
                        v_vf(i)%sf(m, 0:n, 0:p)
                end do
            end do

        elseif (bc_x%end == -2) then     ! Symmetry BC at end

            do j = 1, buff_size

                do i = 1, cont_idx%end
                    v_vf(i)%sf(m + j, 0:n, 0:p) = &
                        v_vf(i)%sf(m - (j - 1), 0:n, 0:p)
                end do

                v_vf(mom_idx%beg)%sf(m + j, 0:n, 0:p) = &
                    -v_vf(mom_idx%beg)%sf(m - (j - 1), 0:n, 0:p)

                do i = mom_idx%beg + 1, sys_size
                    v_vf(i)%sf(m + j, 0:n, 0:p) = &
                        v_vf(i)%sf(m - (j - 1), 0:n, 0:p)
                end do

            end do

        elseif (bc_x%end == -1) then     ! Periodic BC at end

            do i = 1, sys_size
                do j = 1, buff_size
                    v_vf(i)%sf(m + j, 0:n, 0:p) = &
                        v_vf(i)%sf(j - 1, 0:n, 0:p)
                end do
            end do

        else                            ! Processor BC at end

            call s_mpi_sendrecv_conservative_variables_buffers( &
                v_vf, 1, 1)

        end if

    end subroutine s_populate_variables_buffers


    subroutine s_finalize_rhs_module()

        deallocate (q_cons_qp%vf, q_prim_qp%vf)

        deallocate( vL_vf_flat )
        deallocate( vR_vf_flat )

        deallocate( qK_prim_vf_flat )
        deallocate( qK_cons_vf_flat )

        deallocate(flux_vf_flat)
        deallocate(flux_src_vf_flat)

        deallocate(rhs_vf_flat)

    end subroutine s_finalize_rhs_module 

end module m_rhs
