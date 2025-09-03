#import matplotlib
#import matplotlib.pyplot as plt
import numpy as np
from classy_sz import Class
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import math
from decimal import Decimal

### Creates Gaussian Filter ####
def g_fil(lstar):
     def fil(l,lstar):
         return math.exp(-l*(l+1)/(2*lstar**2))
     i = 2
     filter = []
     ell_f = []
     while i <= 5*lstar:
         filter.append(fil(i, lstar))
         ell_f.append(i)
         i = i + 1
     f = np.array(filter)
     ell = np.array(ell_f)
     np.savetxt("/moto/home/mr4449/gauss_fils/gauss_" + str(lstar) + "_fl.txt", np.c_[ell,f], newline='\n')


# font = {'family':'STIXGeneral'}
# axislabelfontsize='large'
# matplotlib.rc('font', **font)

# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.sans-serif": ["Helvetica"]})



import os 
path_to_class_sz = os.getcwd() + '/../../class_sz/'


# # Settings

# In[32]:

websky_Omega_M = 0.31
websky_Omega_B = 0.049
websky_Omega_L = 0.69
websky_h = 0.68
websky_sigma_8 = 0.81
websky_n_s = 0.965

z_min = 0
z_max = 4.6
m_min = 5e13
m_max = 5e15
M_min = Decimal(m_min)
M_max = Decimal(m_max)

Mmin_websky_msun = 1.3e13
Mmax_websky_msun = 1e16

cut_path = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_m_"+'{:.2e}'.format(M_max)

def kszxhalos(filt, lstar = 0):

    if filt == "pl":
        filter_file = "/moto/home/mr4449/planck_fl_A_170422.txt"

    if filt == "so":
        filter_file = "/moto/home/mr4449/so_fl_A_170422.txt"
    
    if filt == "s4":
        filter_file = "/moto/home/mr4449/class_sz/class_sz_auxiliary_files/s4_fl_A_170422.txt"
    
    if filt == "gauss":
        g_fil(lstar)
        filter_file = "/moto/home/mr4449/gauss_fils/gauss_" + str(lstar) + "_fl.txt"
        filt = "gauss"+str(lstar)

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
        'ell_max' : 25000,
        'ell_min' : 2,
        'dlogell': 0.1,

        'hm_consistency':1,

        'use_websky_m200m_to_m200c_conversion': 1,

        #'delta for galaxies': '200m',
        #'delta for matter density': '200m',
        #'delta for electron density': '200c',

        'redshift_epsrel': 1e-4,
        'redshift_epsabs': 1e-100,
        'mass_epsrel':1e-4,
        'mass_epsabs':1e-100,
        'mass function' : 'T08',
    }

    HOD = {
        'sigma_log10M_HOD': 0.00001,
        'M_min_HOD': m_min*websky_h, #Msun/h
        #'M_max_HOD': m_max*websky_h,
        'galaxy_sample': 'custom',
        'full_path_to_dndz_gal': '/moto/hill/users/mr4449/websky_calcs/norm_websky_dndz_' + cut_path + '.txt',
        'z_min': z_min+0.01,
        'z_max': z_max+0.1,
        'M_min': Mmin_websky_msun*websky_h, # all masses in Msun/h
        'M_max': Mmax_websky_msun*websky_h,
    }


    # In[33]:

    i = 1.0
    while i <= 5.0:
        #j = 300
        #while j <= 500:
            #g_fil(j)
            # the parameters needed for the ksz calculations:
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
            M.set({
            'output':'mean_galaxy_bias,kSZ_kSZ_gal fft (1h),kSZ_kSZ_gal fft (2h),kSZ_kSZ_gal fft (3h)',

            'ksz_filter_file' : filter_file,
            
            'dlogell' : 0.1,
            'ell_max' : 20000,
            'ell_min' : 10.0,

            'gas profile' : 'B16', # set NFW profile
            'gas profile mode' : 'agn',
            'normalize_gas_density_profile' : 0,
            'use_xout_in_density_profile_from_enclosed_mass' : 1,
            'x_out_truncated_density_profile (electrons)' : i, 
                    
            'use_fft_for_profiles_transform' : 1,    

            'x_out_truncated_nfw_profile': 1,
                })

            M.compute()
            cl_kSZ_kSZ_g = M.cl_kSZ_kSZ_g()

                # # Plot

                # In[49]:

            path = "/moto/hill/users/mr4449/websky_calcs/class_sz_calcs/"
            fac = (2.726e6)**2*np.asarray(cl_kSZ_kSZ_g['ell'])*(np.asarray(cl_kSZ_kSZ_g['ell'])+1.)/2./np.pi
            np.save(path+"cl_ksz_1h_new_"+  filt + "_" +str(float(i)) + "_"+cut_path,fac*np.asarray(cl_kSZ_kSZ_g['1h']))
            np.save(path+"cl_ksz_2h_new_" + filt + "_" + str(float(i)) + "_"+cut_path,fac*np.asarray(cl_kSZ_kSZ_g['2h']))
            np.save(path+"cl_ksz_3h_new_" + filt + "_" + str(float(i)) + "_"+cut_path,fac*np.asarray(cl_kSZ_kSZ_g['3h']))
            np.save(path+"ell_ksz_new_" + filt + "_" + str(float(i)) + "_"+cut_path, cl_kSZ_kSZ_g['ell'])
            i = i + 0.5

websky_Omega_M = 0.31
websky_Omega_B = 0.049
websky_Omega_L = 0.69
websky_h = 0.68
websky_sigma_8 = 0.81
websky_n_s = 0.965

z_min = 0.2
z_max = 0.5
m_min = 5e13
m_max = 5e15
M_min = Decimal(m_min)
M_max = Decimal(m_max)

Mmin_websky_msun = 1.3e13
Mmax_websky_msun = 1e16

cut_path = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_m_"+'{:.2e}'.format(M_max)

kszxhalos(filt = "gauss", lstar = 800)