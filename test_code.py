import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy_sz import Class
import os
import time
from scipy.integrate import simps
from scipy.integrate import quad


# the parameters needed for cosmology:
# use the last column of Planck 2018 (https://arxiv.org/pdf/1807.06209.pdf) Table 2
# TT,TE,EE+lowE+lensing+BAO
websky_Omega_M = 0.31
websky_Omega_B = 0.049
websky_Omega_L = 0.69
websky_h = 0.68
websky_sigma_8 = 0.81
websky_n_s = 0.965

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
}

M = Class()
M.set(cosmo_params)
params = {
        'output': 'dndlnM',
        # mass function
        'mass function' : 'T08',
        #integration precision settings
        'ndim_redshifts' :80,
        #redshift and mass bounds
        'z_min' : 0.,
        'z_max' : 3.,
        'M_min' : 6e13,
        'M_max' : 5e15,
    
        'skip_cmb':1,
        'skip_pknl':1
}
M.set(params)
M.compute_class_szfast()

m_arr = np.geomspace(params['M_min'],params['M_max'],500)
dndlnm = np.vectorize(M.get_dndlnM_at_z_and_M)
sigma = np.vectorize(M.get_sigma_at_z_and_m)
nu = np.vectorize(M.get_nu_at_z_and_m)
b1 = np.vectorize(M.get_first_order_bias_at_z_and_nu)
b2 = np.vectorize(M.get_second_order_bias_at_z_and_nu)

path_02 = "/moto/hill/users/mr4449/websky_calcs/websky_dndm_0.1_z_0.3_&_6.00e+13_mh_5.00e+15.txt"
path_05 = "/moto/hill/users/mr4449/websky_calcs/websky_dndm_0.4_z_0.6_&_6.00e+13_mh_5.00e+15.txt"
path_1 = "/moto/hill/users/mr4449/websky_calcs/websky_dndm_0.9_z_1.1_&_6.00e+13_mh_5.00e+15.txt"
path_2 = "/moto/hill/users/mr4449/websky_calcs/websky_dndm_1.9_z_2.1_&_6.00e+13_mh_5.00e+15.txt"
dndm02 = np.loadtxt(path_02)
dndm05 = np.loadtxt(path_05)
dndm1 = np.loadtxt(path_1)
dndm2 = np.loadtxt(path_2)

i = 0
m_02 = []
n_02 = []
while i < len(dndm02):
    m_02.append(dndm02[i][0])
    n_02.append(dndm02[i][1])
    i = i + 1
dm_02 = (m_02[len(m_02)-1]-m_02[0])/len(m_02)

i = 0
m_05 = []
n_05 = []
while i < len(dndm05):
    m_05.append(dndm05[i][0])
    n_05.append(dndm05[i][1])
    i = i + 1
dm_05 = (m_05[len(m_05)-1]-m_02[0])/len(m_05)

i = 0
m_1 = []
n_1 = []
while i < len(dndm1):
    m_1.append(dndm1[i][0])
    n_1.append(dndm1[i][1])
    i = i + 1
dm_1 = (m_1[len(m_1)-1]-m_1[0])/len(m_1)

i = 0
m_2 = []
n_2 = []
while i < len(dndm2):
    m_2.append(dndm2[i][0])
    n_2.append(dndm2[i][1])
    i = i + 1
dm_2 = (m_2[len(m_2)-1]-m_2[0])/len(m_2)

dvdz02 = M.get_volume_dVdzdOmega_at_z(0.2)
dvdz05 = M.get_volume_dVdzdOmega_at_z(0.5)
dvdz1 = M.get_volume_dVdzdOmega_at_z(1)
dvdz2 = M.get_volume_dVdzdOmega_at_z(2)

label_size = 20
title_size = 25
legend_size = 13
handle_length = 1.5
fig, (ax1) = plt.subplots(1,1,figsize=(20,10),sharex=True)
ax = ax1
ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=label_size)
plt.setp(ax.get_xticklabels(), fontsize=label_size)
ax.grid( visible=True, which="both", alpha=0.1, linestyle='--')

plt.title('Halo Mass Function z = 2', size=title_size)
ax.set_xlabel(r'$m\quad[M_\mathrm{sun}/h]$',size=title_size)
#ax.set_ylabel(r'$\mathrm{d}n/\mathrm{dln} m\quad[\mathrm{Mpc/h}]^3$',size=title_size,labelpad=8)
ax.set_ylabel(r'$\mathrm{d}N/\mathrm{d} m$',size=title_size,labelpad=8)

def int_dndlnm(x,y):
    return M.get_dndlnM_at_z_and_M(x,y)*M.get_volume_dVdzdOmega_at_z(x)
partial_int = lambda y: quad(int_dndlnm, 1.9,2.1, args=(y,))[0]
partial_int = np.vectorize(partial_int)
m_arr_1 = np.linspace(6e13, 5e15, 100)
plt.xscale("log")
plt.yscale("log")

plt.plot(m_arr_1, partial_int(m_arr_1), label = "Integrated Class SZ")

# z = 0.2
# ax.plot(m_arr,dndlnm(z,m_arr),label='Class SZ')
#ax.plot(m_arr,dndlnm(z,m_arr),label=r'$z=0.2$',alpha=1.)
# z = 0.2
# ax.plot(m_arr,dndlnm(z,m_arr),label='Class SZ z = 0.2')
# z = 1
# ax.plot(m_arr,dndlnm(z,m_arr),label=r'$z=1$',alpha=1.)
#z = 2
#ax.plot(m_arr,dndlnm(z,m_arr),label=r'$z=2$',alpha=1.)

# ax.plot(lnm_02,np.array(n_02)*np.array(lnm_02)/dvdz02,label=r'$websky z=0.2$',alpha=1.,c='black',ls='-.')
# ax.plot(lnm_05,np.array(n_05)*np.array(lnm_05)/dvdz05,label=r'$websky z=0.5$',alpha=1.,c='red',ls='-.')
# ax.plot(lnm_1,np.array(n_1)*np.array(lnm_1)/dvdz1,label=r'$websky z=1$',alpha=1.,c='blue',ls='-.')
# ax.plot(lnm_2,np.array(n_2)*np.array(lnm_2)/dvdz2,label=r'$websky z=2$',alpha=1.,c='yellow',ls='-.')

ax.plot(m_2,np.array(n_2),label='Websky')
#ax.plot(m_02,np.array(m_02)*np.array(n_02)/(dvdz02*dm_02),label='websky with dm z = 0.2')
# ax.plot(m_05,np.array(n_05)/dvdz05,label='websky z = 0.5',alpha=1.)
# ax.plot(m_05,np.array(m_05)*np.array(n_05)/(dvdz05*dm_05),label='websky z with dm z = 0.5',alpha=1.)
# ax.plot(m_1,np.array(n_1)/dvdz1,label='websky z = 1',alpha=1.)
# ax.plot(m_1,np.array(m_1)*np.array(n_1)/(dvdz05*dm_1),label='websky with dm z = 1',alpha=1.)
# ax.plot(m_2,np.array(n_2)/dvdz2,label='websky z = 2')
# ax.plot(m_2,np.array(m_2)*np.array(n_2)/(dvdz05*dm_2),label='websky with dm z = 2')

ax.loglog()

#ax.set_ylim(1e-6,2e-1)
ax.set_xlim(params['M_min'],params['M_max'])

ax.legend(loc=1,frameon=False,framealpha=1,fontsize=20)
plt.savefig("dndm_plots/z_2_6e13_m_5e15.png")
