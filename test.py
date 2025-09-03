import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import time
import healpy as hp
from scipy.interpolate import interp1d
from decimal import Decimal
import math

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

start_time = time.time()

def unison_shuffled_copies(a, b, c):
    p = np.random.permutation(len(a))
    return a[p], b[p], c[p]

def l_to_dl(lp):
    return lp*(lp+1.)/2./np.pi



t_cmb_mk = 2.726e6
websky_h = 0.68
filt = "pl"

z_min = 1
z_max = 2.5
M_min = 6e13
M_max = 5e15
z_m_cut = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_mh_"+'{:.2e}'.format(M_max)
#z_m_cut = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_m_"+'{:.2e}'.format(M_max)
z_m_cut_title = str(z_min)+" < z < "+str(z_max)+ " & "+'{:.2e}'.format(M_min)+" < M$_\u2609$/h < "+'{:.2e}'.format(M_max)
#z_m_cut_title = str(z_min)+" < z < "+str(z_max)+ " & "+'{:.2e}'.format(M_min)+" < M$_\u2609$ < "+'{:.2e}'.format(M_max)
#filt = "so"

#path_web = "/moto/hill/users/mr4449/websky_calcs/cl_tau_halo_shiv_8192_new_cut_"+z_m_cut+".npy"
path_shiv = "/moto/hill/users/mr4449/websky_calcs/cl_ksz_halo_shiv_v4_8192_new_cut_"+filt+"_"+z_m_cut+".npy"
#path_shiv_v8 = "/moto/hill/users/mr4449/websky_calcs/cl_tau_tau_shiv_v8_8192_new_cut_"+z_m_cut+".npy"
#path_shiv = "/moto/hill/users/mr4449/websky_calcs/cl_tsz_halo_shiv_8192_"+z_m_cut+".npy"
path_class = "/moto/hill/users/mr4449/websky_calcs/class_sz_calcs/"

cl_1h = np.load(path_class+"cl_ksz_1h_unenc_"+filt+"_3.0_"+z_m_cut+".npy")
cl_2h = np.load(path_class+"cl_ksz_2h_unenc_"+filt+"_3.0_"+z_m_cut+".npy")
cl_3h = np.load(path_class+"cl_ksz_3h_unenc_"+filt+"_3.0_"+z_m_cut+".npy")
ell_class = np.load(path_class+"ell_ksz_unenc_"+filt+"_3.0_"+z_m_cut+".npy")


#print("For " +z_m_cut + " number of halos is "+str(len(np.load("/moto/hill/users/mr4449/websky_calcs/x_halo_"+z_m_cut+".npy"))))

# cl_web = bin_avg(np.load(path_web),50)[0]
# ell_web = bin_avg(np.load(path_web),50)[1]
# cl_shiv_v8 = bin_avg(np.load(path_shiv_v8),50)[0]
# ell_shiv_v8 = bin_avg(np.load(path_shiv_v8),50)[1]
cl_shiv = bin_avg(np.load(path_shiv),100)[0]
ell_shiv = bin_avg(np.load(path_shiv),100)[1]


# NSIDE = 8192
# w = hp.sphtfunc.pixwin(NSIDE)
# interp_func = interp1d(np.arange(len(w)), w)
# cl_pix = interp_func(ell_web)

# NSIDE = 4096
# w = hp.sphtfunc.pixwin(NSIDE)
# interp_func = interp1d(np.arange(len(w)), w)
# cl_pix_shiv = interp_func(ell_shiv)

plt.rc('text', usetex=False)
plt.rc('font', family='serif')

label_size = 17
title_size = 18
legend_size = 13
handle_length = 1.5

fig, (ax1) = plt.subplots(1,1,figsize=(9,4))
ax = ax1

# plt.xscale("log")
#plt.yscale("log")
#plt.ylim(0,3.5)
plt.xlim(0,10000)
ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=label_size)
plt.setp(ax.get_xticklabels(), fontsize=label_size)
ax.grid( visible=True, which="both", alpha=0.2, linestyle='--')
plt.title("kSZ$^2$" +r'$\times$'+ "Halos " + z_m_cut_title + " " + filt)
#plt.title("Tau "+r'$\times$'+" Halos  " + z_m_cut_title)
ax.plot(ell_class, cl_1h + cl_2h + cl_3h, label = "Class SZ")
# ax.plot(ell_web, l_to_dl(ell_web)*cl_web, label = "Websky")
#ax.plot(ell_web, #l_to_dl(ell_web)*cl_web/cl_pix, label = "Websky")
#ax.plot(ell_shiv_v8, l_to_dl(ell_shiv_v8)*cl_shiv_v8, label = "Pasted B16 Maps v8")
ax.plot(ell_shiv, (t_cmb_mk**2)*l_to_dl(ell_shiv)*cl_shiv/(websky_h), label = "Pasted B16 Maps")
ax.set_xlabel(r"$\ell$",size=title_size)
#ax.set_ylabel(r"$\ell(\ell+1)C_\ell/2\pi\quad$ $[\mu K^2]$",size=title_size)
ax.set_ylabel(r"$\ell(\ell+1)C_\ell/2\pi\quad$",size=title_size)
#ax.set_ylabel(r"$C_\ell$",size=title_size)
plt.legend()
fig.tight_layout()
plt.savefig("shivam_plots/cl_kszxhalo_8192_v4_shiv_new_cut_"+filt+"_"+z_m_cut+".png")









