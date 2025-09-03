#!/usr/bin/env python
# coding: utf-8

# # Intialize

# In[1]:

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy_sz import Class
import os
import time
from scipy.interpolate import interp1d
import scipy.optimize as opt

def bin_avg(data, w):
    ell_avg = []
    data_avg = []
    std = []
    i = 0
    if w ==1:
        data_avg = data
        ell_avg = np.arange(len(data))
    else:
        while i + w < len(data):
            if i + w > len(data):
                data_avg.append(np.sum(data[i:len(data)])/(w + 1))
                ell_avg.append(i)
                std.append(np.std(np.array(data[i:len(data)])))
            else:                   
                data_avg.append(np.sum(data[i:i+w+1])/(w + 1))
                ell_avg.append(i)
                std.append(np.std(np.array(data[i:i+w+1])))
            i = i + w + 1
    return np.array(data_avg), np.array(ell_avg), np.array(std)

    
# the parameters needed for cosmology:
# use the last column of Planck 2018 (https://arxiv.org/pdf/1807.06209.pdf)
# TT,TE,EE+lowE+lensing+BAO

h = 0.6777

cosmo_params = {
'H0':67.77,
'omega_b':0.048206*h**2,
'omega_cdm':0.307115*h**2 - 0.048206*h**2,
'sigma8': 0.818,
#'ln10^{10}A_s': 3.047,
'n_s': 0.9665,
'tau_reio':0.0925,
# Take fixed value for primordial Helium (instead of automatic BBN adjustment)
# 'YHe':0.246,
 
}


common_params = {

'z_min' : 0.005,
'z_max' : 3.0,
'M_min' : 1.0e10, 
'M_max' : 3.5e15,
    

# 'delta for galaxies' : "200c",
# 'delta for matter density' : "200c",
# 'delta for electron density':"200c",    

'mass function' : 'T08M200c',
'concentration parameter' : 'B13',

'redshift_epsabs': 1.0e-40,
'redshift_epsrel': 0.0005,
'mass_epsabs': 1.0e-40,
'mass_epsrel': 0.0005,


'ell_max': 1600.0,
'ell_min': 2.0,
'dell': 10,

'non_linear' : 'hmcode',

'hm_consistency' : 1,
# 'x_outSZ': 4.,
# 'truncate_wrt_rvir':0,
}


import os 
path_to_class_sz = os.getcwd() + '/../../class_sz/'

# a simple conversion from cl's to dl's
def l_to_dl(lp):
    return lp*(lp+1.)/2./np.pi


# In[2]:


import classy_sz
classy_sz.__file__
def chisq(x):
    HOD = "red"
    bahamas = "76"

    path_agora = "/moto/hill/users/mr4449/class_sz_tau/agora_tau/cross_tau" + bahamas + "_" + HOD +".npy"
    agora_true = np.load(path_agora)
    ell_agora_true = np.arange(len(agora_true))
    agora_avg = bin_avg(agora_true, 20)

    # In[3]:

    ksz_params = {
    'output': 'tau_gal_1h,tau_gal_2h',

    #'pressure profile':'B16',
    #'gas profile mode' :'agn',
    
    "ell_min" : 2,
    "ell_max" : 10000,
    'dell': 0,
    'dlogell': 0.2,
    
    'M_min' : 1.0e10, 
    'M_max' : 5e15,

    'gas profile':'B16',
    'gas profile mode' : 'custom', # important to read values of parameters

    "A_rho0" : x[0],
    "A_alpha" : x[1],
    "A_beta" : x[2],

    "alpha_m_rho0" : x[3],
    "alpha_m_alpha" : x[4],
    "alpha_m_beta" : x[5],

    "alpha_z_rho0" : x[6],
    "alpha_z_alpha" : x[7],
    "alpha_z_beta" : x[8],

    #'use_xout_in_density_profile_from_enclosed_mass' : 1,
    'x_out_truncated_density_profile (electrons)': x[9],
    'n_z_m_to_xout' : 30,
    'n_mass_m_to_xout' : 30,

    

    'n_m_density_profile' :30, # default: 100, decrease for faster
    'n_z_density_profile' :30, # default: 100, decrease for faster


    
    'k_min_samp_fftw' : 1e-3,
    'k_max_samp_fftw' : 1e3,
    'N_samp_fftw' : 1024,
    
    
    'hm_consistency' : 1,
    
    
    'use_fft_for_profiles_transform' : 1,
    
    
    'x_min_gas_density_fftw' : 1e-6,
    'x_max_gas_density_fftw' : 1e5,    
    
    }
    if HOD == "blue":
        ###Blue###
        HOD_blue = {
        'sigma_log10M_HOD': 0.01,
        'alpha_s_HOD':    1.06,
        'M1_prime_HOD': 10**12.61, # Msun/h
        'M_min_HOD': 10**11.69, # Msun/h
        'M0_HOD' :0,
        'x_out_truncated_nfw_profile_satellite_galaxies': 1.80,
        'f_cen_HOD' : 1., 
    
        'galaxy_sample': 'unwise',
        'galaxy_sample_id': "blue",
        'UNWISE_dndz_file': "/moto/home/mr4449/class_sz/class_sz_auxiliary_files/normalised_dndz_cosmos_0.txt",
        }
    if HOD == "green":
        ###Green###
        HOD_blue = {
        'sigma_log10M_HOD': 0.03,
        'alpha_s_HOD': 1.14,
        'M1_prime_HOD': 10**13.17, # Msun/h
        'M_min_HOD': 10**12.23, # Msun/h
        'M0_HOD' :0,
        'x_out_truncated_nfw_profile_satellite_galaxies': 0.94,
        'f_cen_HOD' : 1., 
    
        'galaxy_sample': 'unwise',
        'galaxy_sample_id': "blue",
        'UNWISE_dndz_file': "/moto/home/mr4449/class_sz/class_sz_auxiliary_files/normalised_dndz_cosmos_1.txt",
        }
    if HOD == "red":
        ###Red###
        HOD_blue = {
        'sigma_log10M_HOD': 0.07,
        'alpha_s_HOD':    1.76,
        'M1_prime_HOD': 10**13.57, # Msun/h
        'M_min_HOD': 10**12.56, # Msun/h
        'M0_HOD' :0,
        'x_out_truncated_nfw_profile_satellite_galaxies': 1.02,
        'f_cen_HOD' : 1., 
    
        'galaxy_sample': 'unwise',
        'galaxy_sample_id': "blue",
        'UNWISE_dndz_file': "/moto/home/mr4449/class_sz/class_sz_auxiliary_files/normalised_dndz_cosmos_2.txt",
    
        }
    
    # In[38]:


    M = Class()
    M.set(common_params)
    M.set(cosmo_params)
    M.set(ksz_params)
    M.set(HOD_blue)
    M.set({'use_fft_for_profiles_transform' : 1,
    'ndim_redshifts':30, ## precision parameter -- how many z in redshift grid.
    })
    M.compute_class_szfast() # fast mode works
    # M.compute() # slow mode works")
    
    #  In[39]:
    
    l = np.asarray(M.cl_sz()['ell'])
    cl_eg_1h = np.asarray(M.cl_eg()['1h'])
    cl_eg_2h = np.asarray(M.cl_eg()['2h'])
    cl_tot = cl_eg_1h + cl_eg_1h
    
    l = np.insert(l, 0, 0)
    cl_tot = np.insert(cl_tot, 0, 0)
    x_int = ell_agora_true
    interp_func = interp1d(l, cl_tot)
    cl_int = interp_func(x_int)
    cl_bin = bin_avg(cl_int, 20)

    return np.sum(((2*cl_bin[0] - agora_avg[0])**2)/((agora_avg[2])**2))

def tauxhalos(x, HOD, bahamas, chi_squared):
    path_agora = "/moto/hill/users/mr4449/class_sz_tau/agora_tau/cross_tau" + bahamas + "_" + HOD +".npy"
    agora_true = np.load(path_agora)
    ell_agora_true = np.arange(len(agora_true))
    agora_avg = bin_avg(agora_true, 20)

    ksz_params = {
    'output': 'tau_gal_1h,tau_gal_2h',

    #'pressure profile':'B16',
    #'gas profile mode' :'agn',
    
    "ell_min" : 2,
    "ell_max" : 10000,
    'dell': 0,
    'dlogell': 0.2,
    
    'M_min' : 1.0e10, 
    'M_max' : 5e15,

    'gas profile':'B16',
    'gas profile mode' : 'custom', # important to read values of parameters

    "A_rho0" : x[0],
    "A_alpha" : x[1],
    "A_beta" : x[2],

    "alpha_m_rho0" : x[3],
    "alpha_m_alpha" : x[4],
    "alpha_m_beta" : x[5],

    "alpha_z_rho0" : x[6],
    "alpha_z_alpha" : x[7],
    "alpha_z_beta" : x[8],

    #'use_xout_in_density_profile_from_enclosed_mass' : 1,
    'x_out_truncated_density_profile (electrons)': x[9],
    'n_z_m_to_xout' : 30,
    'n_mass_m_to_xout' : 30,

    

    'n_m_density_profile' :30, # default: 100, decrease for faster
    'n_z_density_profile' :30, # default: 100, decrease for faster


    
    'k_min_samp_fftw' : 1e-3,
    'k_max_samp_fftw' : 1e3,
    'N_samp_fftw' : 1024,
    
    
    'hm_consistency' : 1,
    
    
    'use_fft_for_profiles_transform' : 1,
    
    
    'x_min_gas_density_fftw' : 1e-6,
    'x_max_gas_density_fftw' : 1e5,    
    
    }
    if HOD == "blue":
        ###Blue###
        HOD_blue = {
        'sigma_log10M_HOD': 0.01,
        'alpha_s_HOD':    1.06,
        'M1_prime_HOD': 10**12.61, # Msun/h
        'M_min_HOD': 10**11.69, # Msun/h
        'M0_HOD' :0,
        'x_out_truncated_nfw_profile_satellite_galaxies': 1.80,
        'f_cen_HOD' : 1., 
    
        'galaxy_sample': 'unwise',
        'galaxy_sample_id': "blue",
        'UNWISE_dndz_file': "/moto/home/mr4449/class_sz/class_sz_auxiliary_files/normalised_dndz_cosmos_0.txt",
        }
    if HOD == "green":
        ###Green###
        HOD_blue = {
        'sigma_log10M_HOD': 0.03,
        'alpha_s_HOD': 1.14,
        'M1_prime_HOD': 10**13.17, # Msun/h
        'M_min_HOD': 10**12.23, # Msun/h
        'M0_HOD' :0,
        'x_out_truncated_nfw_profile_satellite_galaxies': 0.94,
        'f_cen_HOD' : 1., 
    
        'galaxy_sample': 'unwise',
        'galaxy_sample_id': "blue",
        'UNWISE_dndz_file': "/moto/home/mr4449/class_sz/class_sz_auxiliary_files/normalised_dndz_cosmos_1.txt",
        }
    if HOD == "red":
        ###Red###
        HOD_blue = {
        'sigma_log10M_HOD': 0.07,
        'alpha_s_HOD':    1.76,
        'M1_prime_HOD': 10**13.57, # Msun/h
        'M_min_HOD': 10**12.56, # Msun/h
        'M0_HOD' :0,
        'x_out_truncated_nfw_profile_satellite_galaxies': 1.02,
        'f_cen_HOD' : 1., 
    
        'galaxy_sample': 'unwise',
        'galaxy_sample_id': "blue",
        'UNWISE_dndz_file': "/moto/home/mr4449/class_sz/class_sz_auxiliary_files/normalised_dndz_cosmos_2.txt",
    
        }
    
    # In[38]:


    M = Class()
    M.set(common_params)
    M.set(cosmo_params)
    M.set(ksz_params)
    M.set(HOD_blue)
    M.set({'use_fft_for_profiles_transform' : 1,
    'ndim_redshifts':30, ## precision parameter -- how many z in redshift grid.
    })
    M.compute_class_szfast() # fast mode works
    # M.compute() # slow mode works")
    
    #  In[39]:
    
    l = np.asarray(M.cl_sz()['ell'])
    cl_eg_1h = np.asarray(M.cl_eg()['1h'])
    cl_eg_2h = np.asarray(M.cl_eg()['2h'])
    
    label_size = 17
    title_size = 18
    legend_size = 13
    handle_length = 1.5
    fig, (ax1) = plt.subplots(1,1,figsize=(7,4))
    ax = ax1

    ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
    ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=label_size)
    plt.setp(ax.get_xticklabels(), fontsize=label_size)
    ax.grid( visible=True, which="both", alpha=0.2, linestyle='--')

    plt.title("Agora(Tau) x Halos " r"$\chi^2$" + " = " + str(chi_squared) + " (" + HOD + " " + bahamas + ")")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.plot(l, 2*cl_eg_1h, label = "1h",c='y',lw=0.7)
    ax.plot(l, 2*cl_eg_2h, label = "2h",c='g',lw=0.7)
    ax.plot(l, 2*(cl_eg_1h + cl_eg_2h), label = "Total",c='r',lw=0.7)
    ax.plot(agora_avg[1], agora_avg[0], label = "Agora", c='b',lw=0.7)
    ax.legend(loc=2,ncol = 1,frameon=False,fontsize=14)
    ax.set_xlabel(r"$\ell$",size=title_size)
    ax.set_ylabel(r"$\ell(\ell+1)C_\ell/2\pi\quad $",size=title_size)
    fig.tight_layout()

    plt.savefig("/moto/home/mr4449/chi_sq_min_tau_results/" + HOD + "_" + bahamas + "_test3")

    np.save("/moto/home/mr4449/chi_sq_min_tau_results/1h" + "_" + HOD + "_" + bahamas + "_test3", cl_eg_1h)
    np.save("/moto/home/mr4449/chi_sq_min_tau_results/2h" + "_" + HOD + "_" + bahamas + "_test3", cl_eg_2h)
    np.save("/moto/home/mr4449/chi_sq_min_tau_results/ell" + "_" + HOD + "_" + bahamas + "_test3", l)



bnds = ((2.e3, 5.e3), (0, 5), (1,10), (0,5), (-2,2), (0,5), (-2,2), (-2,2), (-2,2), (1,5))
result = opt.minimize(chisq, [4.e3, 0.88, 3.83, 0.29, -0.03, 0.04, -0.66, 0.19,-0.025, 1], bounds = bnds)
tauxhalos(x=result.x, HOD="red", bahamas="76", chi_squared = result.fun)

print(result.message)

print(result.fun)

print(result.x)

print(result.nit)
print("x0 = [4.e3, 0.88, 3.83, 0.29, -0.03, 0.04, -0.66, 0.19, -0.025, 1]")







