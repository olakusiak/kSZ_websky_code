import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy_sz import Class
import os
import time
import math
from decimal import Decimal
import pickle

websky_Omega_M = 0.31
websky_Omega_B = 0.049
websky_Omega_L = 0.69
websky_h = 0.68
websky_sigma_8 = 0.81
websky_n_s = 0.965

z_min = 1
z_max = 2.5
m_min = 6e13
m_max = 5e15
M_min = Decimal(m_min)
M_max = Decimal(m_max)

Mmin_websky_msun = 1.3e13
Mmax_websky_msun = 1e16

cut_path = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_m_"+'{:.2e}'.format(M_max)

# the parameters needed for cosmology:
# use the last column of Planck 2018 (https://arxiv.org/pdf/1807.06209.pdf)
# TT,TE,EE+lowE+lensing+BAO
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



unWISE_common = {
        'galaxy_sample': 'custom',
        'M0 equal M_min (HOD)':'no',
        #'pressure profile': 'B12',  # check source/input.c for default parameter values of Battaglia et al profile (B12)
        
        'concentration parameter': 'D08',  # B13: Bhattacharya et al 2013
        'hm_consistency':0, #0 for tSZ? and 1 for kSZ

        'use_websky_m200m_to_m200c_conversion': 1,
        'units for tSZ spectrum': 'dimensionless',

        'delta for galaxies': '200m',
        'delta for matter density': '200c',
        'delta for electron density': '200c',

        'redshift_epsrel': 1e-4,
        'redshift_epsabs': 1e-100,
        'mass_epsrel':1e-4,
        'mass_epsabs':1e-100,
        'mass function' : 'T08M200c',


        'n_l_pressure_profile' : 100,
        'n_m_pressure_profile' : 100,
        'n_z_pressure_profile' : 100,
        'x_outSZ': 4.,
        'truncate_gas_pressure_wrt_rvir':0,
        'pressure_profile_epsrel':1e-3,
        'pressure_profile_epsabs':1e-40,

    }

HOD = {
	    'sigma_log10M_HOD': 0.0000001,
        'M_min_HOD': m_min*websky_h, #Class SZ always has units of Msun/h
        'M_max_HOD': m_max*websky_h,
        'galaxy_sample': 'custom',
        'full_path_to_dndz_gal': '/moto/hill/users/mr4449/websky_calcs/norm_websky_dndz_'+cut_path+'.txt',
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

    # In[33]:

ksz_params = {
    #fiducial ksz params

    'k_min_for_pk_in_vrms2' : 2*np.pi/(15.4e3*websky_h),
    'k_max_for_pk_in_vrms2' :2*np.pi/(6*10**(-1)*websky_h) ,
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
M.set(pp_B12)
M.set({
    'output':'m200m_to_m200c,m200c_to_m200m,m200c_to_m500c',
            
    'dlogell' : 0.1,
    'ell_max' : 10000,
    'ell_min' : 10.0,

    'tabulate_rhob_xout_at_m_and_z': 1,
    'gas profile' : 'B16', # set NFW profile
    'gas profile mode' : 'agn',
    #'normalize_gas_density_profile' : 1, #comment out if need be
    'use_xout_in_density_profile_from_enclosed_mass' : 0,
    'x_out_truncated_density_profile (electrons)' : 3.0, 
    'use_bg_eff_in_ksz2g_eff': 0,

    #'effective_galaxy_bias':  'bias_value',
                    
    'use_fft_for_profiles_transform' : 1,    

    'x_out_truncated_nfw_profile': 1, #set to 1 for kSZ calcs, but probably should not do anything for kSZ since B16 is being used
        })

M.compute()

z = 1
M_halo = 6e13
m_halo = Decimal(M_halo)
# m200m = 3e14
# convert to 200c for b16 profile
m200c = M_halo#M.get_m200m_to_m200c_at_z_and_M(z,m200m)
m200m = M.get_m200c_to_m200m_at_z_and_M(z,m200c)


lambda_min = 0.1
lambda_max = 200
n_lambda = 500
lambda_array = np.geomspace(lambda_min,lambda_max,n_lambda)

# store the radial profiles of the gas
rho_gas_b16 = np.vectorize(M.get_gas_profile_at_x_M_z_b16_200c)

# normalized radial array for b16:
x_200c = lambda_array

# dimensonfull radial array:
r200c = M.get_r_delta_of_m_delta_at_z(200,m200c,z)
r = x_200c*r200c
theta_arcmin = M.get_rad_to_arcmin(r/M.get_dA(z))

c200c = M.get_c200c_at_m_and_z_B13(m200c,z)
# rs_200m = r200m/c200m
# xs_200m =  r/rs_200m
rs_200c = r200c/c200c
xs_200c=  r/rs_200c

def x200c_to_theta_arcmin(X):
    return M.get_rad_to_arcmin(r200c*X/M.get_dA(z))
def theta_arcmin_to_x200c(X):
    theta_rad = M.get_arcmin_to_rad(X)
    r = M.get_dA(z)*theta_rad
    return r/r200c

fig = plt.figure()
ax1 = fig.add_subplot(111)


A_rho0 = 4.e3
A_alpha = 0.88
A_beta = 3.83

alpha_m_rho0 = 0.29
alpha_m_alpha = -0.03
alpha_m_beta = 0.04

alpha_z_rho0 = -0.66
alpha_z_alpha = 0.19
alpha_z_beta = -0.025
gamma = -0.2
xc = 0.5

def gas_param(A_0, A_m, A_z):
    return A_0*((m200c/(1e14*websky_h))**A_m)*(1+z)**A_z

alpha_density = gas_param(A_alpha, alpha_m_alpha, alpha_z_alpha)
beta_density = gas_param(A_beta, alpha_m_beta, alpha_z_beta)
rho_density = gas_param(A_rho0, alpha_m_rho0, alpha_z_rho0)

def rho_fit(x):
    return rho_density*((x/xc)**gamma)*(1+(x/xc)**alpha_density)**(-(beta_density+gamma)/alpha_density)
rho_prof = rho_fit(x_200c)

rho_norm_b16 = rho_gas_b16(x_200c,
                           m200c,
                           z,
                           A_rho0 = A_rho0,
                           A_alpha = A_alpha,
                           A_beta = A_beta,
                           alpha_m_rho0 = alpha_m_rho0,
                           alpha_m_alpha = alpha_m_alpha,
                           alpha_m_beta = alpha_m_beta,
                           alpha_z_rho0 = alpha_z_rho0,
                           alpha_z_alpha = alpha_z_alpha,
                           alpha_z_beta = alpha_z_beta,
                           gamma=gamma,
                           xc = xc)/M.get_rho_crit_at_z(z)/M.get_f_b()



path = "/moto/hill/users/mr4449/shivam_maps/B16_gasprofile_test_data.pkl"

# Read the pickle file
with open(path, 'rb') as file:
    data = pickle.load(file)

# Inspect the data structure
print("Type of data:", type(data))
if isinstance(data, dict):
    print("\nKeys in the dictionary:")
    for key in data.keys():
        print(f"- {key}: {type(data[key])}")
elif isinstance(data, (list, tuple)):
    print("\nLength of sequence:", len(data))
    if len(data) > 0:
        print("Type of first element:", type(data[0]))
else:
    print("\nData shape (if available):", getattr(data, 'shape', 'N/A'))

# Now 'data' contains your pickle file contents
# You can print it to see what's inside

r_shiv = np.array(data["r_array"])
m_shiv = np.array(data["M_array"])
z_shiv = np.array(data["z_array"])
a_shiv = np.array(data["scale_fac_array"])
r200c_shiv = np.array(data["r200c_mat"])
rho_shiv = data["rho_fit_mat"]

z_ind = np.where((z_shiv <= z + 0.1) & (z_shiv >= z - 0.1))[0]
z_ind = z_ind[int(len(z_ind)/2)]

m_ind = np.where((m_shiv <= M_halo*(1 + 1/5)) & (m_shiv >= M_halo*(1 - 1/5)))[0]
m_ind = m_ind[int(len(m_ind)/2)]

rho_shiv = rho_shiv[0:, z_ind, m_ind]
r200c_shiv = r200c_shiv[m_ind, z_ind]
a_shiv = a_shiv[z_ind]
x_shiv = a_shiv*r_shiv/r200c_shiv

ax1.plot(x_200c,rho_norm_b16*x_200c**2, label = r'Class SZ',ls='-',c='r',lw=2.)
ax1.plot(x_200c,rho_prof*x_200c**2, label = r'Manual',ls='-.',c='b',lw=2.)
ax1.plot(x_shiv, rho_shiv*x_shiv**2, label = r'Pasted Profile',ls='--',c='purple',lw=2.)


ax1.set_xscale('log')
ax1.legend(frameon=False)
ax1.set_xlabel(r'$x=r/r_{200c}$')
ax1.set_ylabel(r'$\rho_\mathrm{gas}x^2/f_\mathrm{b}\rho_\mathrm{crit}$')
ax1.set_title(r'z = ' + str(z) + r' & M = '+'{:.2e}'.format(m_halo)+r'M'+'$_\u2609$'+r'/h')

plt.grid(which='both',alpha=0.2)
secax = ax1.secondary_xaxis('top', functions=(x200c_to_theta_arcmin, theta_arcmin_to_x200c))

secax.set_xlabel(r'$\theta = r/d_A$ [arcmin]')
ax1.set_xlim(1e-1,1e1)
#ax1.set_xlim(0,100)

fig.tight_layout()
plt.savefig('shivam_plots/Battaglia_density_profiles_z'+str(z)+'_m'+'{:.2e}'.format(m_halo)+'.png')
# plt.savefig('shivam_plots/r_radius.png')

