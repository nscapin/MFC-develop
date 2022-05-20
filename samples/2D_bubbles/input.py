#!/usr/bin/env python3
# Dependencies and Logistics ===================================================
# Command to navigate between directories
from os import chdir

# Command to acquire directory path
from os.path import dirname

# Command to acquire script name and module search path
from sys import argv, path
import math


ps = 1.2*101325
gam = 1.4
rho = 1.

c_l = math.sqrt( 1.4*ps/rho )

vel = 230.

leng = 1.


Ny = 50
Nx = Ny
dx = leng/Nx

time_end = 4*leng/vel
cfl = 0.1

dt = cfl * dx/c_l 
Nt = int(time_end/dt)

print('c_l, Ma', c_l, vel/c_l)
print('Nt', Nt)


# Navigating to script directory
if len(dirname(argv[0])) != 0: chdir(dirname(argv[0]))

# Adding master_scripts directory to module search path
mfc_dir = '../../src'; path[:0] = [mfc_dir + '/master_scripts']

# Command to execute the MFC components
from m_python_proxy import f_execute_mfc_component

# Serial or parallel computational engine
engine = 'parallel'
# ==============================================================================


# Case Analysis Configuration ==================================================

# Selecting MFC component
comp_name = argv[1].strip()

# Configuring case dictionary
case_dict =                                                                    \
    {                                                                          \
                    # Logistics ================================================
                    'case_dir'                     : '\'.\'',                  \
                    'run_time_info'                : 'T',                      \
                    'nodes'                        : 1,                        \
                    'ppn'                          : 24,                        \
                    'queue'                        : 'normal',                 \
                    'walltime'                     : '24:00:00',               \
                    'mail_list'                    : '',                       \
                    # ==========================================================
                                                                               \
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 :  0.0,                   \
                    'x_domain%end'                 :  leng,                   \
                    'y_domain%beg'                 :  0.,                  \
                    'y_domain%end'                 :  leng,           \
                    'm'                            : int(Nx),                      \
                    'n'                            : int(Ny),                        \
                    'p'                            : 0,                        \
                    'cyl_coord'                    : 'F',                      \
                    'dt'                           : dt,                   \
                    't_step_start'                 : 0,                        \
                    't_step_stop'                  : Nt,                     \
                    't_step_save'                  : int(Nt/100.),                     \
		    # ==========================================================
                                                                               \
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 10,                        \
                    'model_eqns'                   : 2,                        \
                    'alt_soundspeed'               : 'T',                      \
                    'num_fluids'                   : 2,                        \
		    'adv_alphan'                   : 'T',                      \
		    'mpp_lim'                      : 'T',                      \
		    'mixture_err'                  : 'T',                      \
		    'time_stepper'                 : 3,                        \
                    'weno_vars'                    : 2,                        \
                    'weno_order'                   : 5,                        \
                    'weno_eps'                     : 1.E-16,                   \
                    'mapped_weno'                  : 'T',                      \
                    'null_weights'                 : 'F',                      \
                    'mp_weno'                      : 'T',                      \
                    'weno_Re_flux'                 : 'F',                      \
		    'riemann_solver'               : 2,                        \
                    'wave_speeds'                  : 1,                        \
                    'avg_state'                    : 2,                        \
                    'bc_x%beg'                     : -9,                       \
                    'bc_x%end'                     : -9,                       \
                    'bc_y%beg'                     : -1,                      \
                    'bc_y%end'                     : -1,                      \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        \
                    'precision'                    : 2,                        \
                    'prim_vars_wrt'                :'T',                       \
		    'parallel_io'                  :'T',                       \
		    # ==========================================================
                                                                                
		    # Patch 1: Background  ============================
                    'patch_icpp(1)%geometry'       : 3,                        \
                    'patch_icpp(1)%x_centroid'     : 0.5,                  \
                    'patch_icpp(1)%y_centroid'     : 0.5,          \
                    'patch_icpp(1)%length_x'       : leng,                   \
                    'patch_icpp(1)%length_y'       : leng,         \
                    'patch_icpp(1)%vel(1)'         : 0.E+00,                   \
                    'patch_icpp(1)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(1)%pres'           : 101325.,                   \
                    'patch_icpp(1)%alpha_rho(1)'   : 1000.0,                \
                    'patch_icpp(1)%alpha_rho(2)'   : 0.E+00,                   \
                    'patch_icpp(1)%alpha(1)'       : 1.E+00,                   \
                    'patch_icpp(1)%alpha(2)'       : 0.E+00,                   \
                    # ==========================================================

			
                    'patch_icpp(2)%geometry'       : 2,                        \
                    'patch_icpp(2)%x_centroid'     : 0.25,                 \
                    'patch_icpp(2)%y_centroid'     : 0.25,                  \
                    'patch_icpp(2)%radius'         : leng/10.,                  \
                    'patch_icpp(2)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(2)%vel(1)'         : 0.,                   \
                    'patch_icpp(2)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(2)%pres'           : 1.2*101325.,                   \
                    'patch_icpp(2)%alpha_rho(1)'   : 0.E+00,                   \
                    'patch_icpp(2)%alpha_rho(2)'   : 0.167,                  \
                    'patch_icpp(2)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(2)%alpha(2)'       : 1.E+00,                   \

                    'patch_icpp(3)%geometry'       : 2,                        \
                    'patch_icpp(3)%x_centroid'     : 0.5,                 \
                    'patch_icpp(3)%y_centroid'     : 0.25,                  \
                    'patch_icpp(3)%radius'         : leng/10.,                  \
                    'patch_icpp(3)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(3)%vel(1)'         : 0.,                   \
                    'patch_icpp(3)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(3)%pres'           : 0.8*101325.,                   \
                    'patch_icpp(3)%alpha_rho(1)'   : 0.E+00,                   \
                    'patch_icpp(3)%alpha_rho(2)'   : 0.167,                  \
                    'patch_icpp(3)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(3)%alpha(2)'       : 1.E+00,                   \

                    'patch_icpp(4)%geometry'       : 2,                        \
                    'patch_icpp(4)%x_centroid'     : 0.75,                 \
                    'patch_icpp(4)%y_centroid'     : 0.25,                  \
                    'patch_icpp(4)%radius'         : leng/10.,                  \
                    'patch_icpp(4)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(4)%vel(1)'         : 0.,                   \
                    'patch_icpp(4)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(4)%pres'           : 1.2*101325.,                   \
                    'patch_icpp(4)%alpha_rho(1)'   : 0.E+00,                   \
                    'patch_icpp(4)%alpha_rho(2)'   : 0.167,                  \
                    'patch_icpp(4)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(4)%alpha(2)'       : 1.E+00,                   \
                    
                    'patch_icpp(5)%geometry'       : 2,                        \
                    'patch_icpp(5)%x_centroid'     : 0.25,                 \
                    'patch_icpp(5)%y_centroid'     : 0.5,                  \
                    'patch_icpp(5)%radius'         : leng/10.,                  \
                    'patch_icpp(5)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(5)%vel(1)'         : 0.,                   \
                    'patch_icpp(5)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(5)%pres'           : 0.8*101325.,                   \
                    'patch_icpp(5)%alpha_rho(1)'   : 0.E+00,                   \
                    'patch_icpp(5)%alpha_rho(2)'   : 0.167,                  \
                    'patch_icpp(5)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(5)%alpha(2)'       : 1.E+00,                   \


                    'patch_icpp(6)%geometry'       : 2,                        \
                    'patch_icpp(6)%x_centroid'     : 0.5,                 \
                    'patch_icpp(6)%y_centroid'     : 0.5,                  \
                    'patch_icpp(6)%radius'         : leng/10.,                  \
                    'patch_icpp(6)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(6)%vel(1)'         : 0.,                   \
                    'patch_icpp(6)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(6)%pres'           : 1.2*101325.,                   \
                    'patch_icpp(6)%alpha_rho(1)'   : 0.E+00,                   \
                    'patch_icpp(6)%alpha_rho(2)'   : 0.167,                  \
                    'patch_icpp(6)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(6)%alpha(2)'       : 1.E+00,                   \

                    'patch_icpp(7)%geometry'       : 2,                        \
                    'patch_icpp(7)%x_centroid'     : 0.75,                 \
                    'patch_icpp(7)%y_centroid'     : 0.5,                  \
                    'patch_icpp(7)%radius'         : leng/10.,                  \
                    'patch_icpp(7)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(7)%vel(1)'         : 0.,                   \
                    'patch_icpp(7)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(7)%pres'           : 0.8*101325.,                   \
                    'patch_icpp(7)%alpha_rho(1)'   : 0.E+00,                   \
                    'patch_icpp(7)%alpha_rho(2)'   : 0.167,                  \
                    'patch_icpp(7)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(7)%alpha(2)'       : 1.E+00,                   \

                    'patch_icpp(8)%geometry'       : 2,                        \
                    'patch_icpp(8)%x_centroid'     : 0.25,                 \
                    'patch_icpp(8)%y_centroid'     : 0.75,                  \
                    'patch_icpp(8)%radius'         : leng/10.,                  \
                    'patch_icpp(8)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(8)%vel(1)'         : 0.,                   \
                    'patch_icpp(8)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(8)%pres'           : 1.2*101325.,                   \
                    'patch_icpp(8)%alpha_rho(1)'   : 0.E+00,                   \
                    'patch_icpp(8)%alpha_rho(2)'   : 0.167,                  \
                    'patch_icpp(8)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(8)%alpha(2)'       : 1.E+00,                   \
                    
                    'patch_icpp(9)%geometry'       : 2,                        \
                    'patch_icpp(9)%x_centroid'     : 0.5,                 \
                    'patch_icpp(9)%y_centroid'     : 0.75,                  \
                    'patch_icpp(9)%radius'         : leng/10.,                  \
                    'patch_icpp(9)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(9)%vel(1)'         : 0.,                   \
                    'patch_icpp(9)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(9)%pres'           : 0.8*101325.,                   \
                    'patch_icpp(9)%alpha_rho(1)'   : 0.E+00,                   \
                    'patch_icpp(9)%alpha_rho(2)'   : 0.167,                  \
                    'patch_icpp(9)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(9)%alpha(2)'       : 1.E+00,                   \


                    'patch_icpp(10)%geometry'       : 2,                        \
                    'patch_icpp(10)%x_centroid'     : 0.75,                 \
                    'patch_icpp(10)%y_centroid'     : 0.75,                  \
                    'patch_icpp(10)%radius'         : leng/10.,                  \
                    'patch_icpp(10)%alter_patch(1)' : 'T',                      \
                    'patch_icpp(10)%vel(1)'         : 0.,                   \
                    'patch_icpp(10)%vel(2)'         : 0.E+00,                  \
                    'patch_icpp(10)%pres'           : 1.2*101325.,                   \
                    'patch_icpp(10)%alpha_rho(1)'   : 0.E+00,                   \
                    'patch_icpp(10)%alpha_rho(2)'   : 0.167,                  \
                    'patch_icpp(10)%alpha(1)'       : 0.E+00,                   \
                    'patch_icpp(10)%alpha(2)'       : 1.E+00,                   \



                    # Fluids Physical Parameters ===============================

                    'fluid_pp(1)%gamma'            : 1.E+00/(4.4E+00-1.E+00),  \
                    'fluid_pp(1)%pi_inf'           : 4.4E+00*6.E+08/(4.4E+00-1.E+00), \
                    'fluid_pp(2)%gamma'            : 1.E+00/(1.4E+00-1.E+00),  \
                    'fluid_pp(2)%pi_inf'           : 0.E+00,                   \

	            # ==========================================================
                    'Monopole'                  : 'T',                  \
                    'num_mono'                  : 2,                  \
                    'Mono(1)%loc(1)'            : 0.55,  \
                    'Mono(1)%loc(2)'            : 0.55,  \
                    'Mono(1)%npulse'            : 1, \
                    'Mono(1)%dir'               : 1., \
                    'Mono(1)%pulse'             : 1, \
                    'Mono(1)%mag'               : 101325, \
                    'Mono(1)%length'            : 0.05, \
                    'Mono(1)%support'           : 1, \

                    'Mono(2)%loc(1)'            : 0.45,  \
                    'Mono(2)%loc(2)'            : 0.45,  \
                    'Mono(2)%npulse'            : 1, \
                    'Mono(2)%dir'               : -1., \
                    'Mono(2)%pulse'             : 1, \
                    'Mono(2)%mag'               : 101325, \
                    'Mono(2)%length'            : 0.05, \
                    'Mono(2)%support'           : 1, \
    }

# Executing MFC component
f_execute_mfc_component(comp_name, case_dict, mfc_dir, engine)

# ==============================================================================
