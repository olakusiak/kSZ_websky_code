import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy_sz import Class
import os
import time
from decimal import Decimal

def bin_avg(data, w):
    ell = np.arange(len(data))
    data1 = data
    ell_avg = []
    data_avg = []
    i = 0
    if w == 1:
        data_avg = data
        ell_avg = np.arange(len(data))
    else:
        while i + w < len(data1):
            if i + w > len(data1):
                data_avg.append(np.sum(data1[i:len(data)])/(w + 1))
                ell_avg.append(i)
            else:                   
                data_avg.append(np.sum(data1[i:i+w+1])/(w + 1))
                ell_avg.append(i)
            i = i + w + 1
    return np.array(data_avg), np.array(ell_avg)


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})


# the parameters needed for cosmology:
# use the last column of Planck 2018 (https://arxiv.org/pdf/1807.06209.pdf) Table 2
# TT,TE,EE+lowE+lensing+BAO

websky_Omega_M = 0.31
websky_Omega_B = 0.049
websky_Omega_L = 0.69
websky_h = 0.68
websky_sigma_8 = 0.81
websky_n_s = 0.965

z_min = 0.5
z_max = 0.7
m_min = 5e13
m_max = 5e15
M_min = Decimal(m_min)
M_max = Decimal(m_max)

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
    'm_ncdm': 0.06   

}

HOD_blue = {
'sigma_log10M_HOD': 0.00000000001, 
'alpha_s_HOD':    1.76,
'M1_prime_HOD': m_max*websky_h, # Msun/h
'M_min_HOD': m_min*websky_h, # Msun/h
'M0_HOD' :0,
'x_out_truncated_nfw_profile_satellite_galaxies': 1.02,
'f_cen_HOD' : 1.,
}   


HOD_common = {
'z_min': z_min - 0.05,
'z_max': z_max + 1,
'M_min': (m_min - 1e10)*websky_h,
'M_max': (m_max + 1e10)*websky_h,
    
'dlogell': 0.2,
'ell_max': 8000.0,
'ell_min': 2.0,
    
# precisions params:
# 'k_min_for_pk_class_sz' :  0.001,
# 'k_max_for_pk_class_sz' :  60.0,
# 'k_per_decade_class_sz' :  50,
# 'P_k_max_h/Mpc' :  50.0,

'redshift_epsabs': 1.0e-40,
'redshift_epsrel': 0.0001,
'mass_epsabs': 1.0e-40,
'mass_epsrel': 0.001,
# 'ndim_masses': 150,
'ndim_redshifts': 50,


'hm_consistency': 1,


'delta for galaxies': "200m",
'delta for matter density': "200m",
'mass function': 'T08',
'concentration parameter': 'B13' ,
    
'M0 equal M_min (HOD)':'no',
'x_out_truncated_nfw_profile': 1.0,
}


M = Class()
M.set(cosmo_params)
M.set(HOD_common)
M.set(HOD_blue)
M.set({
'output' : 'gal_gal_1h,gal_gal_2h',
'galaxy_sample': 'custom',
'full_path_to_dndz_gal': '/moto/hill/users/mr4449/websky_calcs/norm_websky_dndz_'+cut_path+'.txt'
})



M.compute_class_szfast()
cl_gal_lens = M.cl_gg()

ell = np.asarray(cl_gal_lens['ell'])
fac = ell*(ell+1.)/2./np.pi
cl_kg_1h = np.asarray(cl_gal_lens['1h'])
cl_kg_2h = np.asarray(cl_gal_lens['2h'])

np.save("/moto/hill/users/mr4449/websky_calcs/class_sz_calcs/cl_gal_gal_1h_"+cut_path, cl_kg_1h)
np.save("/moto/hill/users/mr4449/websky_calcs/class_sz_calcs/cl_gal_gal_2h_"+cut_path, cl_kg_2h)
np.save("/moto/hill/users/mr4449/websky_calcs/class_sz_calcs/ell_gal_gal_test_"+cut_path, ell)




