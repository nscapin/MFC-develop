#!/usr/bin/env python2

# Dependencies and Logistics ===================================================

# Command to navigate between directories
from os import chdir

# Command to acquire directory path
from os.path import dirname

# Command to acquire script name and module search path
from sys import argv, path

import math
# Navigating to script directory
if len(dirname(argv[0])) != 0: chdir(dirname(argv[0]))

# Adding master_scripts directory to module search path
mfc_dir = '../../src'; path[:0] = [mfc_dir + '/master_scripts']

# Command to execute the MFC components
from m_python_proxy import f_execute_mfc_component

# Serial or parallel computational engine
engine = 'serial'
# engine = 'parallel'

#Dimensional Variables
Gam     = 7.1
rho0    = 1000.0
x0      = 1.0
p0      = 101325.0

L0      = 10.0
c0      = 1475.0
u0      = 0.1*c0

#Dimensionless Variables via (rho,p,x)
rho     = 1.0
x       = 1.0
patm    = 1.0

L       = L0/x0
u       = u0/math.sqrt(p0/rho0)
c       = c0/u0
pi_inf  = (306.E+06)/p0

Nx = 100
Ny = Nx
dx = L/float(Nx)

cfl = 0.2
dt = 0.3*cfl*dx/c

# T = 1*L/u
Nt = 30000
# Nt = int(T/dt)
Nfiles = 300

print 'dt ', dt
print 'Nt ', Nt

# Case Analysis Configuration ==================================================

# Selecting MFC component
comp_name = argv[1].strip()

np = 1
# if (comp_name=='post_process'): np = 1 

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
                    'debug' : 'T', \
                    # ==========================================================
                                                                                \
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 : -L/2,        \
                    'x_domain%end'                 :  L/2,        \
                    'y_domain%beg'                 : -L/2,        \
                    'y_domain%end'                 :  L/2,        \
                    # 'z_domain%beg'                 : -2.5E-03/x0,        \
                    # 'z_domain%end'                 :  2.5E-03/x0,        \
                    # 'stretch_x'                    : 'F',                       \
                    # 'a_x'                          : 2.000000000000E-00, \
                    # 'x_a'                          : -2.5,         \
                    # 'x_b'                          :  2.5,         \
                    # 'stretch_y'                    : 'F',                       \
                    # 'stretch_z'                    : 'F',                       \
                    'm'                            : Nx,                       \
                    'n'                            : Ny,                       \
                    'p'                            : 0,                       \
                    'dt'                           : dt,   \
                    # 'dt'                           : 0.25*3.0E-08*c0/x0,   \
                    't_step_start'                 : 0,                         \
                    't_step_stop'                  : Nt,                       \
                    't_step_save'                  : 100, \
                    # 't_step_save'                  : math.floor(float(Nt)/float(Nfiles)),                        \
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
                    # 'weno_Re_flux'                 : 'T',                      \
                    'char_decomp'                  :'F',                        \
                    'mapped_weno'                  :'T',                        \
                    'riemann_solver'               : 2,                         \
                    'wave_speeds'                  : 1,                         \
                    'avg_state'                    : 2,                         \
                    'commute_err'                  :'F',                        \
                    'split_err'                    :'F',                        \
                    'bc_x%beg'                     :-1,                         \
                    'bc_x%end'                     :-1,                         \
                    'bc_y%beg'                     :-1,                         \
                    'bc_y%end'                     :-1,                         \
                    # 'bc_z%beg'                     :-6,                         \
                    # 'bc_z%end'                     :-6,                         \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        \
                    'precision'                    : 2,                        \
                    'prim_vars_wrt'                :'T',                       \
                    # 'c_wrt'                        :'T',                       \
                    # 'omega_wrt(3)'                 :'T',                       \
                    'parallel_io'                  :'T',                       \
                    # 'fd_order'                     : 1,                        \
                    #'schlieren_wrt'                :'T',                     \
                    # 'probe_wrt'                    :'T',                    \
                    # 'num_probes'                   : 1,                     \
                    # 'probe(1)%x'                   : 0.,                    \
                    # 'probe(1)%y'                   : 0.,                    \
                    # 'probe(1)%z'                   : 0.,                    \
                    # ==========================================================
                                                                               \
                    # Patch 1: Water (left) ====================================
                    'patch_icpp(1)%geometry'       : 7,                         \
                    'patch_icpp(1)%x_centroid'     : 0.,            \
                    'patch_icpp(1)%y_centroid'     : 0.,            \
                    'patch_icpp(1)%length_x'       : L,            \
                    'patch_icpp(1)%length_y'       : L,            \
                    'patch_icpp(1)%vel(1)'         : u,      \
                    'patch_icpp(1)%vel(2)'         : u,      \
                    'patch_icpp(1)%pres'           : patm,     \
                    'patch_icpp(1)%alpha_rho(1)'   : rho*1.0,            \
                    'patch_icpp(1)%alpha(1)'       : rho,      \
                    # ==========================================================
                    
                    # Fluids Physical Parameters ===============================
                    'fluid_pp(1)%gamma'            : 1.0/(Gam-1.0),                \
                    'fluid_pp(1)%pi_inf'           : Gam*pi_inf/(Gam-1.0),   \
                    # 'fluid_pp(1)%Re(1)'            : 10000.0,      \
                    # 'fluid_pp(1)%Re(2)'            : 0.01,      \
                    # ==========================================================
    }

# Executing MFC component
f_execute_mfc_component(comp_name, case_dict, mfc_dir, engine)

# ==============================================================================
