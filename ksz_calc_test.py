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

z_min = 0.5
z_max = 1
m_min = 5e13
m_max = 5e15
M_min = Decimal(m_min)
M_max = Decimal(m_max)

Mmin_websky_msun = 1.3e12
Mmax_websky_msun = 1e16

cut_path = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_m_"+'{:.2e}'.format(M_max)

cosmo_params = {
    'omega_b': websky_Omega_B*websky_h**2.,
    'omega_cdm': (websky_Omega_M-websky_Omega_B)*websky_h**2.,
    'h': websky_h,
    'tau_reio': 0.0543,
    'sigma8': websky_sigma_8,
    'n_s': websky_n_s, 
    'k_pivot': 0.05,
    'N_ncdm': 1,
    'N_ur': 2.0328,
    'm_ncdm': 0, 
    'use_websky_m200m_to_m200c_conversion': 1,  

}

# best-fit from Kusiak et al. https://arxiv.org/pdf/2203.12583.pdf

HOD = {
    'sigma_log10M_HOD': 0.00001,
    'M_min_HOD': m_min*websky_h, #Msun/h
    #‘M_max_HOD’: M_min*websky_h,
    'galaxy_sample': 'custom',
    'full_path_to_dndz_gal': '/moto/hill/users/mr4449/websky_calcs/norm_websky_dndz_' + cut_path + '.txt',
    #'x_out_truncated_nfw_profile': 1.0,
    'z_min': z_min-0.01,
    'z_max': z_max+0.1,
    'M_min': Mmin_websky_msun*websky_h, # all masses in Msun/h
    'M_max': Mmax_websky_msun*websky_h,
}

unWISE_common = {
'galaxy_sample': 'custom',
'M0 equal M_min (HOD)':'no',
'x_out_truncated_nfw_profile': 4.0,
    
'nfw_profile_epsabs' : 1e-33,
'nfw_profile_epsrel' : 0.001,
    
    
# 'use_fft_for_profiles_transform' : 1,
    
    
'x_min_gas_density_fftw' : 1e-5,
'x_max_gas_density_fftw' : 1e4,
    
    
'redshift_epsabs': 1.0e-40,
'redshift_epsrel': 0.001,
'mass_epsabs': 1.0e-40,
'mass_epsrel': 0.001,



'hm_consistency': 1,


'delta for galaxies': "200c",
'delta for matter density': "200c",
'mass function': 'T08M200c',
'concentration parameter': 'B13' ,
    

}

ksz_params = {
#fiducial ksz params

#'k_min_for_pk_in_vrms2' : 15.4e-4,
#'k_max_for_pk_in_vrms2' : 1.67,
#'k_per_decade_for_vrms2' : 40,

'k_min_for_pk_class_sz' : 0.001,
'k_max_for_pk_class_sz' : 50.0,
'k_per_decade_class_sz' : 50,
'P_k_max_h/Mpc' : 50.0,

'nfw_profile_epsabs' : 1e-33,
'nfw_profile_epsrel' : 0.001,


'ndim_masses' : 80,
'ndim_redshifts' : 50,




'n_k_density_profile' : 50,
'n_m_density_profile' : 50,
'n_z_density_profile' : 50,
'k_per_decade_for_pk' : 50,
'z_max_pk' : 4.0,
    
# some settings to try more points to avoid numerical noise in some cases:
# 'ndim_masses' : 100,
# 'ndim_redshifts' : 100,
# 'n_ell_density_profile' : 100,
# 'n_m_density_profile' : 100,
# 'n_z_density_profile' : 100,
    

# slow:
# 'n_z_psi_b1g' : 100,
# 'n_l_psi_b1g' : 400,

# 'n_z_psi_b2g' : 100,
# 'n_l_psi_b2g' : 400,

# 'n_z_psi_b2t' : 100,
# 'n_l_psi_b2t' : 400,

# 'n_z_psi_b1t' : 100,
# 'n_l_psi_b1t' : 100,

# 'n_z_psi_b1gt' : 100,
# 'n_l_psi_b1gt' : 100,
    
    
# fast:
'n_z_psi_b1g' : 50,
'n_l_psi_b1g' : 50,

'n_z_psi_b2g' : 50,
'n_l_psi_b2g' : 50,

'n_z_psi_b2t' : 50,
'n_l_psi_b2t' : 50,

'n_z_psi_b1t' : 50,
'n_l_psi_b1t' : 50,

'n_z_psi_b1gt' : 50,
'n_l_psi_b1gt' : 50,

'N_samp_fftw' : 1024, # fast: 800 ;  slow: 2000
'l_min_samp_fftw' : 1e-9,
'l_max_samp_fftw' : 1e9,
    
}

M = Class()
M.set(cosmo_params)
M.set(HOD)
M.set(unWISE_common)
M.set(ksz_params)
M.set({
'output':'mean_galaxy_bias,kSZ_kSZ_gal fft (1h),kSZ_kSZ_gal fft (2h),kSZ_kSZ_gal fft (3h)',    
    
# for effective approach calculation of kSZ2g, i.e.,kSZ_kSZ_gal_hf also set:
# 'N_kSZ2_gal_multipole_grid' :  70,
# 'N_kSZ2_gal_theta_grid' :  70,
# 'ell_min_kSZ2_gal_multipole_grid' : 2.,
# 'ell_max_kSZ2_gal_multipole_grid' : 2e5,

'ksz_filter_file' :'/moto/home/mr4449/planck_fl_A_170422.txt',           

'dlogell' : 0.1,
'ell_max' : 10000.0,
'ell_min' : 10.0,

'gas profile' : 'B16', # set NFW profile
'gas profile mode' : 'agn',
'normalize_gas_density_profile' : 0,
'use_xout_in_density_profile_from_enclosed_mass' : 1,
    
# 'use_fft_for_profiles_transform' : 1,    

'use_bg_at_z_in_ksz2g_eff' : 1,
'non_linear' : 'halofit',
      })

M.compute()
cl_kSZ_kSZ_g = M.cl_kSZ_kSZ_g()

path = "/moto/hill/users/mr4449/websky_calcs/class_sz_calcs/"
fac = (2.726e6)**2*np.asarray(cl_kSZ_kSZ_g['ell'])*(np.asarray(cl_kSZ_kSZ_g['ell'])+1.)/2./np.pi
np.save(path+"cl_ksz_1h_" + cut_path,fac*np.asarray(cl_kSZ_kSZ_g['1h']))
np.save(path+"cl_ksz_2h_" + cut_path,fac*np.asarray(cl_kSZ_kSZ_g['2h']))
np.save(path+"cl_ksz_3h_" + cut_path,fac*np.asarray(cl_kSZ_kSZ_g['3h']))
np.save(path+"ell_ksz_" + cut_path, cl_kSZ_kSZ_g['ell'])