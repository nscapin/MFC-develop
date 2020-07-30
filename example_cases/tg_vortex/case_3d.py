#!/usr/bin/python

# Dependencies and Logistics ===================================================

# Command to navigate between directories
from os import chdir

# Command to acquire directory path
from os.path import dirname

# Command to acquire script name and module search path
from sys import argv, path

# Navigating to script directory
if len(dirname(argv[0])) != 0: chdir(dirname(argv[0]))

# Adding master_scripts directory to module search path
mfc_dir = '../../'; path[:0] = [mfc_dir + '/master_scripts']

# Command to execute the MFC components
from m_python_proxy import f_execute_mfc_component

# Serial or parallel computational engine
#engine = 'serial'
engine = 'parallel'

x0      = 1.e-03
rho0    = 1.e+03
c0      = 1475.
p0      = rho0*c0*c0

Gam     = 7.1
pi_inf  = (306.e+06)/p0
patm    = 101325./p0
pa      = 0.1*(1.E+06) / p0
# ==============================================================================

print 'dt ', 0.25*3.0E-08*c0/x0

# Case Analysis Configuration ==================================================

# Selecting MFC component
comp_name = argv[1].strip()

np = 24
if (comp_name=='post_process'): np = 1 

# Configuring case dictionary
case_dict =                                                                     \
    {                                                                           \
                    # Logistics ================================================
                    'case_dir'                     : '\'.\'',                   \
                    'run_time_info'                : 'T',                       \
                    'nodes'                        : 1,                         \
                    'ppn'                          : np,                        \
                    'queue'                        : 'normal',                  \
                    'walltime'                     : '24:00:00',                \
                    'mail_list'                    : '',      \
                    # ==========================================================
                                                                                \
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 : -4.0E-03/x0,        \
                    'x_domain%end'                 :  4.0E-03/x0,        \
                    'y_domain%beg'                 : -2.5E-03/x0,        \
                    'y_domain%end'                 :  2.5E-03/x0,        \
                    'z_domain%beg'                 : -2.5E-03/x0,        \
                    'z_domain%end'                 :  2.5E-03/x0,        \
                    'stretch_x'                    : 'T',                       \
                    'a_x'                          : 2.000000000000E-00, \
                    'x_a'                          : -2.5,         \
                    'x_b'                          :  2.5,         \
                    'stretch_y'                    : 'F',                       \
                    'stretch_z'                    : 'F',                       \
                    'm'                            : 500,                       \
                    'n'                            : 50,                       \
                    'p'                            : 50,                       \
                    'dt'                           : 0.25*3.0E-08*c0/x0,   \
                    't_step_start'                 : 0,                         \
                    't_step_stop'                  : 10000,                       \
                    't_step_save'                  : 50,                        \
                    # ==========================================================
                                                                               \
                    # Simulation Algorithm Parameters ==========================
                    'model_eqns'                   : 2,                         \
                    'num_fluids'                   : 1,                         \
                    'num_patches'                  : 1,                         \
                    'adv_alphan'                   : 'T',                       \
                    'mpp_lim'                      : 'F',                       \
                    'time_stepper'                 : 3,                         \
                    'weno_vars'                    : 2,                         \
                    'weno_order'                   : 5,                         \
                    'weno_eps'                     : 1.00000000000000E-16,      \
                    'char_decomp'                  :'T',                        \
                    'mapped_weno'                  :'T',                        \
                    'riemann_solver'               : 2,                         \
                    'wave_speeds'                  : 1,                         \
                    'avg_state'                    : 2,                         \
                    'commute_err'                  :'F',                        \
                    'split_err'                    :'F',                        \
                    'Ac_src'                       :'T',                        \
                    'Ac_src_id'                    : 2,                         \
                    'Ac_src_f'                     : 300000.*x0/c0,  \
                    'Ac_src_a'                     : pa,                   \
                    'Ac_src_c'                     : 1.,                        \
                    'Ac_src_fp'                    : -7.E-03/x0,               \
                    'Ac_src_res'                   : 2.5 * 2.*4.7E-03/x0/172.,  \
                    #'Ac_src_res'                   : 5.0E-04/x0/2.0,           \
                    'bc_x%beg'                     :-6,                         \
                    'bc_x%end'                     :-6,                         \
                    'bc_y%beg'                     :-6,                         \
                    'bc_y%end'                     :-6,                         \
                    'bc_z%beg'                     :-6,                         \
                    'bc_z%end'                     :-6,                         \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        \
                    'precision'                    : 2,                        \
                    'prim_vars_wrt'                :'T',                       \
                    'c_wrt'                        :'T',                       \
                    'parallel_io'                  :'T',                       \
                    'fd_order'                     : 1,                     \
                    #'schlieren_wrt'                :'T',                   \
                    'probe_wrt'                    :'T',                    \
                    'num_probes'                   : 1,                     \
                    'probe(1)%x'                   : 0.,                    \
                    'probe(1)%y'                   : 0.,                    \
                    'probe(1)%z'                   : 0.,                    \
                    # ==========================================================
                                                                               \
                    # Patch 1: Water (left) ====================================
                    'patch_icpp(1)%geometry'       : 9,                         \
                    'patch_icpp(1)%x_centroid'     : 0.,            \
                    'patch_icpp(1)%y_centroid'     : 0.,            \
                    'patch_icpp(1)%z_centroid'     : 0.,            \
                    'patch_icpp(1)%length_x'       : 5.E-02/x0,            \
                    'patch_icpp(1)%length_y'       : 5.E-03/x0,            \
                    'patch_icpp(1)%length_z'       : 5.E-03/x0,            \
                    'patch_icpp(1)%vel(1)'         : 0.,      \
                    'patch_icpp(1)%vel(2)'         : 0.,      \
                    'patch_icpp(1)%vel(3)'         : 0.,      \
                    'patch_icpp(1)%pres'           : patm,     \
                    'patch_icpp(1)%alpha_rho(1)'   : 1.,            \
                    'patch_icpp(1)%alpha(1)'       : 1.,      \
                    # ==========================================================
                    
                    # Fluids Physical Parameters ===============================
                    'fluid_pp(1)%gamma'            : 1.0/(Gam-1.0),                \
                    'fluid_pp(1)%pi_inf'           : Gam*pi_inf/(Gam-1.0),   \
                    # ==========================================================
    }

# Executing MFC component
f_execute_mfc_component(comp_name, case_dict, mfc_dir, engine)

# ==============================================================================
