import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy_sz import Class
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import math
from decimal import Decimal

font = {'family':'STIXGeneral'}
axislabelfontsize='large'
matplotlib.rc('font', **font)

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

import os 
path_to_class_sz = os.getcwd() + '/../../class_sz/'

websky_Omega_M = 0.31
websky_Omega_B = 0.049
websky_Omega_L = 0.69
websky_h = 0.68
websky_sigma_8 = 0.81
websky_n_s = 0.965

z_min = 0.3
z_max = 0.5
m_min = 5e13
m_max = 5e15
M_min = Decimal(m_min)
M_max = Decimal(m_max)

Mmin_websky_msun = 1.3e12
Mmax_websky_msun = 1e16

cut_path = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_m_"+'{:.2e}'.format(M_max)

websky_cosmo = {
    'omega_b': websky_Omega_B*websky_h**2.,
    'omega_cdm': (websky_Omega_M-websky_Omega_B)*websky_h**2.,
    'h': websky_h,
    'tau_reio': 0.0543,
    'sigma8': websky_sigma_8,
    'n_s': websky_n_s, 
    'N_ncdm': 1,
    'N_ur': 2.0328,
    'm_ncdm': 0,
    'use_websky_m200m_to_m200c_conversion': 1,  
}

unWISE_common = {
    #'pressure profile': 'B12',  # check source/input.c for default parameter values of Battaglia et al profile (B12)
    'concentration parameter': 'D08',  # B13: Bhattacharya et al 2013
    'ell_max' : 25000,
    'ell_min' : 2,
    'dlogell': 0.1,
    'units for tSZ spectrum': 'dimensionless',
    'n_l_pressure_profile' : 100,
    'n_m_pressure_profile' : 100,
    'n_z_pressure_profile' : 100,
    'x_outSZ': 4.,
    'truncate_gas_pressure_wrt_rvir':0,
    'hm_consistency':0,
    'pressure_profile_epsrel':1e-3,
    'pressure_profile_epsabs':1e-40,
    'redshift_epsrel': 1e-4,
    'redshift_epsabs': 1e-100,
    'mass_epsrel':1e-4,
    'mass_epsabs':1e-100,
    'mass function' : 'T08',
}


HOD = {
    'sigma_log10M_HOD': 0.00001,
    'M_min_HOD': m_min*websky_h, #Msun/h
    #‘M_max_HOD’: M_min*websky_h,
    'galaxy_sample': 'custom',
    'full_path_to_dndz_gal': '/moto/hill/users/mr4449/websky_calcs/norm_websky_dndz_' + cut_path + '.txt',
    #'x_out_truncated_nfw_profile': 1.0,
    'z_min': z_min-0.1,
    'z_max': z_max+0.1,
    'M_min': Mmin_websky_msun*websky_h, # all masses in Msun/h
    'M_max': Mmax_websky_msun*websky_h,
}


pp_B12 = {
         'pressure profile': 'B12',
         'x_outSZ': 4,
         'concentration parameter': 'B13',
         'delta for electron pressure':'200c',
         #pressure precision
         'n_z_pressure_profile': 80,
}

pp_fft = {
         'n_z_pressure_profile': 100,
         'n_m_pressure_profile' : 100,
         'n_l_pressure_profile' : 100,
         'l_min_gas_pressure_profile' :  1.e-2,
         'l_max_gas_pressure_profile' :  5.e4,
         'use_fft_for_profiles_transform' : 1,
         'N_samp_fftw' : 1024,
         'x_min_gas_pressure_fftw' : 1e-4,
         'x_max_gas_pressure_fftw' : 1e3,
}

M = Class()
M.set(HOD)
M.set(websky_cosmo)
M.set(unWISE_common)
M.set(pp_B12) 
#M.set({'output':'gal_gal_1h, gal_gal_2h'})
M.set({'output' : 'tSZ_gal_1h, tSZ_2h',
'x_out_truncated_nfw_profile': 2.0,
'hm_consistency': 0,
'include_gk_counterterms_in_gk': 0,
'include_k_counterterms_in_gk': 1,
'include_g_counterterms_in_gk': 0,
'hm_consistency_ngbar' : 0
})
      
M.compute_class_szfast()

cl_yg_b12 = M.cl_yg()

cl_ell_b12=np.array(cl_yg_b12['ell'])
cl_1h_b12=np.array(cl_yg_b12['1h'])
cl_2h_b12=np.array(cl_yg_b12['2h'])

np.save("/moto/hill/users/mr4449/websky_calcs/class_sz_calcs/cl_tsz_1h_"+cut_path, cl_1h_b12)
np.save("/moto/hill/users/mr4449/websky_calcs/class_sz_calcs/cl_tsz_2h_"+cut_path, cl_2h_b12)
np.save("/moto/hill/users/mr4449/websky_calcs/class_sz_calcs/ell_tsz_"+cut_path, cl_ell_b12)

print("TSZ Calculation Complete " + cut_path)