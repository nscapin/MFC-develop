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

    real(kind(0d0)), private, allocatable, dimension(:) :: q_cons_buff_send
    real(kind(0d0)), private, allocatable, dimension(:) :: q_cons_buff_recv
    !$acc declare create( q_cons_buff_send, q_cons_buff_recv )

    real(kind(0d0)), allocatable, dimension(:,:,:) :: divu
    real(kind(0d0)), allocatable, dimension(:,:,:) :: bub_adv_src
    real(kind(0d0)), allocatable, dimension(:,:,:,:) :: bub_r_src, bub_v_src, bub_p_src, bub_m_src

    integer, private :: err_code, ierr

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

        ! For MPI send/recv, 1D only!
        allocate(q_cons_buff_send(0:-1 + buff_size*sys_size))
        allocate(q_cons_buff_recv(0:ubound(q_cons_buff_send, 1)))

        if (bubbles) then
            allocate( divu(0:m,0:n,0:p) )
            allocate( bub_adv_src(0:m,0:n,0:p) )
            allocate( bub_r_src(1:nb,0:m,0:n,0:p) )
            allocate( bub_v_src(1:nb,0:m,0:n,0:p) )
            allocate( bub_p_src(1:nb,0:m,0:n,0:p) )
            allocate( bub_m_src(1:nb,0:m,0:n,0:p) )
        end if

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

        !$acc data copyin(qK_cons_vf_flat) copyout(rhs_vf_flat,flux_vf_flat, vL_vf_flat, vR_vf_flat) create(qK_prim_vf_flat, flux_src_vf_flat, bub_adv_src, bub_r_src, bub_v_src, bub_p_src, bub_m_src, divu)
        do it_t = 1,t_step_stop
            call nvtxStartRange("Time step")

            start_time = mpi_wtime()

            !if (proc_rank == num_procs - 1) then
            !    do k = 3,3
            !        !sys_size
            !        do j = ix%beg,ix%end
            !            print*, 'Pre:  Rank: ', &
            !                proc_rank, j, qK_cons_vf_flat(j,0,0,k) 
            !        end do
            !    end do
            !end if

            call nvtxStartRange("RHS-Pop. var. buffers")
            call s_populate_conservative_variables_buffers()
            call nvtxEndRange

            !!$acc update self(qK_cons_vf_flat)
            !if (proc_rank == num_procs - 1) then
            !    do k = 3,3
            !        !sys_size
            !        do j = ix%beg,ix%end
            !            print*, 'Post: Rank: ', &
            !                proc_rank, j, qK_cons_vf_flat(j,0,0,k) 
            !        end do
            !    end do
            !end if

            ! print*, 'after cons var buff proc rank', proc_rank

            call nvtxStartRange("RHS-Convert to prim")
            if (bubbles) then
                call s_convert_conservative_to_primitive_variables_bubbles_acc( &
                    qK_cons_vf_flat, qK_prim_vf_flat, &
                    ix, iy, iz)
            else
                call s_convert_conservative_to_primitive_variables_acc( &
                    qK_cons_vf_flat, qK_prim_vf_flat, &
                    ix, iy, iz)
            end if
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
            if (bubbles) then
                call s_hllc_riemann_solver_bubbles( &
                                      vR_vf_flat, vL_vf_flat, &
                                      flux_vf_flat, &
                                      flux_src_vf_flat )
            else
                call s_hllc_riemann_solver( &
                                      vR_vf_flat, vL_vf_flat, &
                                      flux_vf_flat, &
                                      flux_src_vf_flat )
            end if
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
            !! async(1)
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
            !! async(2)
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
            !! !$acc wait


            if (bubbles) then
                call s_get_divergence()
                call s_compute_bubble_source(1, qK_prim_vf_flat, qK_cons_vf_flat, divu, &
                                             bub_adv_src, bub_r_src, bub_v_src)
                do k = 0, m
                    rhs_vf_flat(k,0,0,alf_idx) = rhs_vf_flat(k,0,0,alf_idx) + bub_adv_src(k,0,0)
                    do j = 1,nb
                        rhs_vf_flat(k,0,0,bub_idx_rs(j)) = rhs_vf_flat(k,0,0,bub_idx_rs(j)) + bub_r_src(k,0,0,j)
                        rhs_vf_flat(k,0,0,bub_idx_vs(j)) = rhs_vf_flat(k,0,0,bub_idx_vs(j)) + bub_v_src(k,0,0,j)
                    end do
                end do
            end if
        
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

            ! !$acc update self(qK_cons_vf_flat)

            ! print*, 'end TS proc rank', proc_rank
        end do
        !$acc end data


        ! Some print statements
        !if (proc_rank == num_procs - 1) then
        !    do k = 3,3
        !        !sys_size
        !        do j = ix%beg,ix%end
        !            print*, 'C, Rank: ', &
        !                proc_rank, j, qK_cons_vf_flat(j,0,0,k) 
        !        end do
        !    end do
        !end if

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

        ! if (proc_rank == num_procs - 1) then
        !     do k = 1,sys_size
        !         do j = 0,m
        !             print*, 'Rank,RHS: ', &
        !                 proc_rank, j, k, rhs_vf_flat(j,0,0,k) 
        !         end do
        !     end do
        ! end if

        do i = 1, sys_size
            nullify (q_cons_qp%vf(i)%sf, q_prim_qp%vf(i)%sf)
        end do

    end subroutine s_alt_rhs

    subroutine s_get_divergence()

        integer :: j,k,l
        
        !$acc data present(dx)
        !$acc parallel loop collapse (3) gang vector 
        do l = 0,p
            do k = 0,n
                do j = 0,m
                    divu(j,k,l) = 0.5d0/dx(j)*(qK_prim_vf_flat(j+1, k, l,cont_idx_e+1) - &
                                               qK_prim_vf_flat(j-1, k, l,cont_idx_e+1)) 
                end do
            end do
        end do
        !$acc end parallel loop
        !$acc end data

    end subroutine s_get_divergence


    subroutine s_populate_conservative_variables_buffers() 
        ! This is called every RHS evaluation

        integer :: i, j, k

        !$acc data present(qK_cons_vf_flat) 
        if (bc_xb == -1) then
            ! Periodic BC at beginning
            !$acc kernels
            do i = 1, sys_size
                do j = 1, buff_size
                    qK_cons_vf_flat(-j, 0, 0, i) = &
                        qK_cons_vf_flat(m - (j - 1), 0, 0, i)
                end do
            end do
            !$acc end kernels
        else
            ! Processor BC at beginning
            call s_mpi_sendrecv_conservative_variables_buffers_acc(-1)
        end if

        if (bc_xe == -1) then
            ! Periodic BC at end

            !$acc kernels
            do i = 1, sys_size
                do j = 1, buff_size
                    qK_cons_vf_flat(m + j, 0, 0, i) = &
                        qK_cons_vf_flat(j - 1, 0, 0, i)
                end do
            end do
            !$acc end kernels
        else                            
            ! Processor BC at end
            call s_mpi_sendrecv_conservative_variables_buffers_acc(1)
        end if
        !$acc end data

    end subroutine s_populate_conservative_variables_buffers


    subroutine s_mpi_sendrecv_conservative_variables_buffers_acc(pbc_loc)

        ! real(kind(0d0)), dimension(-4:,0:,0:,1:), intent(INOUT) :: q_cons_vf_flat
        integer, intent(IN) :: pbc_loc

        integer :: i, j, k, l, r

        !! x-dir only
        ! print*, 'proc rank and pbc loc', proc_rank, pbc_loc

        !$acc data present(q_cons_buff_send,q_cons_buff_recv)
        if (pbc_loc == -1) then
        ! PBC at the beginning

            if (bc_xe >= 0) then
            ! PBC at the beginning and end

                ! Packing buffer to be sent to bc_x%end
                !$acc kernels
                do l = 0, p
                    do k = 0, n
                        do j = m - buff_size + 1, m
                            do i = 1, sys_size
                                r = (i - 1) + sys_size* &
                                    ((j - m - 1) + buff_size*((k + 1) + (n + 1)*l))
                                q_cons_buff_send(r) = qK_cons_vf_flat(j, k, l, i)
                            end do
                        end do
                    end do
                end do
                !$acc end kernels


                ! !$acc host_data use_device( q_cons_buff_recv, q_cons_buff_send )
                !$acc update host(q_cons_buff_send)
                call MPI_SENDRECV( &
                    q_cons_buff_send(0), &
                    buff_size*sys_size, &
                    MPI_DOUBLE_PRECISION, bc_xe, 0, &
                    q_cons_buff_recv(0), &
                    buff_size*sys_size, &
                    MPI_DOUBLE_PRECISION, bc_xb, 0, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                !$acc update device(q_cons_buff_recv)

                ! !$acc end host_data
                ! !$acc wait
            else
                ! PBC at the beginning only

                ! Packing buffer to be sent to bc_x%beg
                !$acc kernels
                do l = 0, p
                    do k = 0, n
                        do j = 0, buff_size - 1
                            do i = 1, sys_size
                                r = (i - 1) + sys_size* &
                                    (j + buff_size*(k + (n + 1)*l))
                                q_cons_buff_send(r) = qK_cons_vf_flat(j, k, l, i)
                            end do
                        end do
                    end do
                end do
                !$acc end kernels


                ! Send/receive buffer to/from bc_x%beg/bc_x%beg
                ! !$acc host_data use_device( q_cons_buff_send, q_cons_buff_recv )
                !$acc update host(q_cons_buff_send)
                call MPI_SENDRECV( &
                    q_cons_buff_send(0), &
                    buff_size*sys_size, &
                    MPI_DOUBLE_PRECISION, bc_xb, 1, &
                    q_cons_buff_recv(0), &
                    buff_size*sys_size, &
                    MPI_DOUBLE_PRECISION, bc_xb, 0, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                !$acc update device(q_cons_buff_recv)
                ! !$acc end host_data
                ! !$acc wait
            end if

            ! Unpacking buffer received from bc_x%beg
            !$acc kernels
            do l = 0, p
                do k = 0, n
                    do j = -buff_size, -1
                        do i = 1, sys_size
                            r = (i - 1) + sys_size* &
                                (j + buff_size*((k + 1) + (n + 1)*l))
                            qK_cons_vf_flat(j, k, l, i) = q_cons_buff_recv(r)
                        end do
                    end do
                end do
            end do
            !$acc end kernels

        else
            ! PBC at the end

            if (bc_xb >= 0) then 
                ! PBC at the end and beginning

                ! Packing buffer to be sent to bc_x%beg
                !$acc kernels
                do l = 0, p
                    do k = 0, n
                        do j = 0, buff_size - 1
                            do i = 1, sys_size
                                r = (i - 1) + sys_size* &
                                    (j + buff_size*(k + (n + 1)*l))
                                q_cons_buff_send(r) = qK_cons_vf_flat(j, k, l, i)
                            end do
                        end do
                    end do
                end do
                !$acc end kernels

                ! Send/receive buffer to/from bc_x%beg/bc_x%end
                ! !$acc host_data use_device( q_cons_buff_send, q_cons_buff_recv )

                !$acc update host(q_cons_buff_send)
                call MPI_SENDRECV( &
                    q_cons_buff_send(0), &
                    buff_size*sys_size, &
                    MPI_DOUBLE_PRECISION, bc_xb, 1, &
                    q_cons_buff_recv(0), &
                    buff_size*sys_size, &
                    MPI_DOUBLE_PRECISION, bc_xe, 1, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                ! !$acc update device(q_cons_buff_recv)

                ! !$acc end host_data
                ! !$acc wait

            else
                ! PBC at the end only

                ! Packing buffer to be sent to bc_x%end
                !$acc kernels
                do l = 0, p
                    do k = 0, n
                        do j = m - buff_size + 1, m
                            do i = 1, sys_size
                                r = (i - 1) + sys_size* &
                                    ((j - m - 1) + buff_size*((k + 1) + (n + 1)*l))
                                q_cons_buff_send(r) = qK_cons_vf_flat(j, k, l, i)
                            end do
                        end do
                    end do
                end do
                !$acc end kernels

                ! Send/receive buffer to/from bc_x%end/bc_x%end

                !!$acc host_data use_device( q_cons_buff_send, q_cons_buff_recv )

                !$acc update host(q_cons_buff_send)
                call MPI_SENDRECV( &
                    q_cons_buff_send(0), &
                    buff_size*sys_size*(n + 1)*(p + 1), &
                    MPI_DOUBLE_PRECISION, bc_xe, 0, &
                    q_cons_buff_recv(0), &
                    buff_size*sys_size*(n + 1)*(p + 1), &
                    MPI_DOUBLE_PRECISION, bc_xe, 1, &
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                !$acc update device(q_cons_buff_recv)

                ! !$acc end host_data
                ! !$acc wait


            end if

            ! Unpacking buffer received from bc_x%end
            !$acc kernels
            do l = 0, p
                do k = 0, n
                    do j = m + 1, m + buff_size
                        do i = 1, sys_size
                            r = (i - 1) + sys_size* &
                                ((j - m - 1) + buff_size*(k + (n + 1)*l))
                            qK_cons_vf_flat(j, k, l, i) = q_cons_buff_recv(r)
                        end do
                    end do
                end do
            end do
            !$acc end kernels

        end if
        !$acc end data


    end subroutine s_mpi_sendrecv_conservative_variables_buffers_acc


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
