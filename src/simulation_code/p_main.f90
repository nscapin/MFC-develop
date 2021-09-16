!>
!! @file p_main.f90
!! @brief Contains program p_main
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief  Quasi-conservative, shock- and interface- capturing finite-volume
!!              scheme for the multicomponent Navier-Stokes equations. The system
!!              is augmented with the relevant advection equations to capture the
!!              material interfaces and closed by the stiffened equation of state
!!              as well as any required mixture relations. The effects of surface
!!              tension are included and modeled through a volume force that acts
!!              across the diffuse material interface regions. The implementation
!!              specifics of surface tension may be found in the work by Perigaud
!!              and Saurel (2005). Note that both viscous and capillarity effects
!!              are only available in the volume fraction model.
program p_main

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_start_up             !< Reading and checking procedures for the input
    !< and the initial condition and grid data files

    use m_variables_conversion !< State variables type conversion procedures

    use m_weno                 !< Weighted and essentially non-oscillatory (WENO)
                               !! schemes for spatial reconstruction of variables

    use m_riemann_solvers      !< Exact and approximate Riemann problem solvers

    use m_cbc                  !< Characteristic boundary conditions (CBC)

    use m_rhs                  !< Right-hand-side (RHS) evaluation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_derived_variables    !< Derived variables evaluation procedures

    use m_time_steppers        !< Time-stepping algorithms

    use m_qbmm                 !< Quadrature MOM

    use openacc

    use nvtx

    ! ==========================================================================

    ! implicit none

    integer :: t_step !< Iterator for the time-stepping loop
    integer :: num_devices, num_nodes, ppn, my_device_num

    call system_clock(COUNT=cpu_start, COUNT_RATE=cpu_rate)

    ! Initializing MPI execution environment
    call s_mpi_initialize()

    ! num_nodes = 2
    ! num_devices = acc_get_num_devices(acc_device_default)
    ! ppn =  num_procs / num_nodes
    ! my_device_num = floor(mod(proc_rank,ppn) * num_devices / (1d0*ppn))
    ! call acc_set_device_num(my_device_num,acc_device_default)
    ! print*, 'rank, device', proc_rank, my_device_num, acc_get_device_num(acc_device_default)

    ! The rank 0 processor assigns default values to the user inputs prior to
    ! reading them in from the input file. Next, the user inputs are read and
    ! their consistency is checked. The identification of any inconsistencies
    ! will result in the termination of the simulation.
    if (proc_rank == 0) then
        call s_assign_default_values_to_user_inputs()
        call s_read_input_file()
        call s_check_input_file()
    end if



    ! Broadcasting the user inputs to all of the processors and performing the
    ! parallel computational domain decomposition. Neither procedure has to be
    ! carried out if the simulation is in fact not truly executed in parallel.
    call s_mpi_bcast_user_inputs()
    call s_initialize_parallel_io()
    call s_mpi_decompose_computational_domain()
    if (proc_rank == 0) print*, 'Number of MPI ranks:', num_procs


    ! Computation of parameters, allocation of memory, association of pointers,
    ! and/or the execution of any other tasks needed to properly configure the
    ! modules. The preparations below DO NOT DEPEND on the grid being complete.
    call s_initialize_global_parameters_module()
    call s_initialize_mpi_proxy_module()
    call s_initialize_variables_conversion_module()
    call s_initialize_start_up_module()
    call s_initialize_riemann_solvers_module()
    call s_initialize_rhs_module()
    call s_initialize_data_output_module()
    call s_initialize_derived_variables_module()
    call s_initialize_time_steppers_module()
    if (qbmm) call s_initialize_qbmm_module()


    ! Associate pointers for serial or parallel I/O
    if (parallel_io .neqv. .true.) then
        s_read_data_files => s_read_serial_data_files
        s_write_data_files => s_write_serial_data_files
    else
        s_read_data_files => s_read_parallel_data_files
        s_write_data_files => s_write_parallel_data_files
    end if

    ! Reading in the user provided initial condition and grid data
    call s_read_data_files(q_cons_ts(1)%vf)
    if (model_eqns == 3) call s_initialize_internal_energy_equations(q_cons_ts(1)%vf)

    ! Populating the buffers of the grid variables using the boundary conditions
    ! print*, 'Pop grid var buffs p_main'
    call s_populate_grid_variables_buffers()
    ! print*, 'done!'

    call s_populate_variables_buffers(q_cons_ts(1)%vf)

    ! Computation of parameters, allocation of memory, association of pointers,
    ! and/or execution of any other tasks that are needed to properly configure
    ! the modules. The preparations below DO DEPEND on the grid being complete.
    call s_initialize_weno_module()
    call s_initialize_cbc_module()

    call s_initialize_derived_variables()

    ! Setting the time-step iterator to the first time-step
    t_step = t_step_start
    if (t_step == 0) then
        mytime = 0d0
    else
        mytime = t_step*dt
    end if
    finaltime = t_step_stop*dt
    dt0 = dt

    ! Time-stepping Loop =======================================================
    do
        ! call nvtxStartRange("Main loop")
        if (proc_rank == 0) then
            if (time_stepper == 23) then
                print *, '------------', mytime/finaltime*100d0, 'percent done'
            else
                print *, '------ Time step ', t_step, 'of', t_step_stop, '----'
            end if
        end if
        mytime = mytime + dt

        call s_compute_derived_variables(t_step)
        if (DEBUG) print *, 'Computed derived vars'

        ! Total-variation-diminishing (TVD) Runge-Kutta (RK) time-steppers
        if (time_stepper == 1) then
            call s_1st_order_tvd_rk(t_step)
        elseif (time_stepper == 2) then
            call s_2nd_order_tvd_rk(t_step)
        elseif (time_stepper == 3) then
            call s_3rd_order_tvd_rk(t_step)
        elseif (time_stepper == 4) then
            call s_4th_order_rk(t_step)
        elseif (time_stepper == 23) then
            call s_23_order_tvd_rk(t_step)
        else
            call s_5th_order_rk(t_step)
        end if

        ! Time-stepping loop controls
        ! call nvtxEndRange

        exit
        ! if (t_step == t_step_stop) then
        !     exit
        ! else
        !     t_step = t_step + 1
        ! end if

        call system_clock(cpu_end)

    end do
    ! ==========================================================================

    ! Disassociate pointers for serial and parallel I/O
    s_read_data_files => null()
    s_write_data_files => null()

    ! Deallocation and/or disassociation procedures for the modules
    call s_finalize_time_steppers_module()
    call s_finalize_derived_variables_module()
    call s_finalize_data_output_module()
    call s_finalize_rhs_module()
    call s_finalize_cbc_module()
    call s_finalize_riemann_solvers_module()
    call s_finalize_weno_module()
    call s_finalize_start_up_module()
    call s_finalize_variables_conversion_module()
    call s_finalize_mpi_proxy_module()
    call s_finalize_global_parameters_module()

    ! print*, 'pmain proc rank', proc_rank

    ! Terminating MPI execution environment
    call s_mpi_finalize()

end program p_main
