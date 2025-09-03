import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from math import pi as pi
from scipy.interpolate import interp1d
from decimal import Decimal


def bin_avg(data, w):
    ell = np.arange(len(data))
    data1 = data
    ell_avg = []
    data_avg = []
    i = 0
    if w ==1:
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

def plot_z_m_cut(z_min, z_max, m_min, m_max, x):

    M_min = Decimal(m_min)
    M_max = Decimal(m_max)

    z_m_cut = "_"+str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_m_"+'{:.2e}'.format(M_max)
    z_m_cut_title = str(z_min)+"<z<"+str(z_max)+" & "+'{:.2e}'.format(M_min)+"<m<"+'{:.2e}'.format(M_max)
    
    sig_type = ["gal_lens", "gal_gal", "tsz", "ksz"]
    sig_type2 = ["kap_halo", "auto_halo", "tsz_halo", "ksz_halo"]

    path_n = "/moto/hill/users/mr4449/websky_calcs/x_halo" + z_m_cut + ".npy"
    r = np.load(path_n)
    N_halos = len(r)

    path = "/moto/hill/users/mr4449/websky_calcs/class_sz_calcs/"
    cl_1h = np.load(path + "cl_" + sig_type[x] +"_1h" + z_m_cut + ".npy")
    cl_2h = np.load(path + "cl_" + sig_type[x] +"_2h" + z_m_cut + ".npy")
    ell = np.load(path + "ell_" + sig_type[x] + z_m_cut + ".npy")


    path_web = "/moto/hill/users/mr4449/websky_calcs/cl_" + sig_type2[x] + z_m_cut + ".npy"
    cl_web = bin_avg(np.load(path_web),20)[0]
    ell_web = bin_avg(np.load(path_web), 20)[1]

    if x == 0:
        NSIDE = 4096
        w = hp.sphtfunc.pixwin(NSIDE)
        interp_func = interp1d(np.arange(len(w)), w)
        cl_pix = interp_func(ell_web)

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
        plt.title("CMBlens x Halos ("+z_m_cut_title+")")
        #plt.yscale("log")
        #plt.xscale("log")
        plt.xlim(0,10000)
        ax.plot(ell_web, ell_web*(ell_web+1)*cl_web/(2*np.pi*(cl_pix)), label = "Pixwin",c='b',lw=0.7)
        ax.plot(ell_web, ell_web*(ell_web+1)*cl_web/(2*np.pi), label = "Websky",c='r',lw=0.7)
        #ax.plot(ell_web, cl_web/(cl_pix)**2, label = "Pixwin",c='r',lw=0.7)
        #ax.plot(ell_web, cl_web, label = "Websky",c='purple',lw=0.7)
        ax.plot(ell, cl_1h + cl_2h, label = "Class SZ",c='purple',lw=0.7)
        #ax.plot(ell, cl_1h, label = "1 Halo",c='yellow',lw=0.7)
        #ax.plot(ell, cl_2h, label = "2 Halo",c='g',lw=0.7)
        ax.legend(loc=2,ncol = 1,frameon=False,fontsize=14)
        ax.set_xlabel(r"$\ell$",size=title_size)
        ax.set_ylabel(r"$\ell(\ell+1)C_\ell/2\pi\quad$",size=title_size)
        fig.tight_layout()
        plt.savefig("/moto/home/mr4449/class_sz_vs_websky_plots/class_vs_websky(cmbxgal" +z_m_cut + ").png")
    
    if x == 1:
        NSIDE = 8192
        w = hp.sphtfunc.pixwin(NSIDE)
        interp_func = interp1d(np.arange(len(w)), w)
        cl_pix = interp_func(ell_web)

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
        N = 10**6
        plt.title("Halos x Halos ("+z_m_cut_title+")")
        #plt.yscale("log")
        #plt.xscale("log")
        plt.ylim(0,1000)
        ax.plot(ell_web, ell_web*(ell_web+1)*cl_web/(2*np.pi*(cl_pix)), label = "Websky",c='b',lw=0.7)
        ax.plot(ell, cl_1h + cl_2h + ell*(ell+1)*(4*np.pi/N_halos)/(2*np.pi), label = "Total",c='r',lw=0.7)
        ax.plot(ell, cl_1h, label = "1 Halo",c='yellow',lw=0.7)
        ax.plot(ell, cl_2h, label = "2 Halo",c='g',lw=0.7)
        ax.plot(ell, ell*(ell+1)*(4*np.pi/N_halos)/(2*np.pi), label = "Shot Noise", c='purple')
        ax.legend(loc=2,ncol = 1,frameon=False,fontsize=14)
        ax.set_xlabel(r"$\ell$",size=title_size)
        ax.set_ylabel(r"$\ell(\ell+1)C_\ell/2\pi\quad$",size=title_size)
        fig.tight_layout()
        plt.savefig("/moto/home/mr4449/class_sz_vs_websky_plots/class_vs_websky(galxgal" + z_m_cut + ").png")

    if x == 2:
        NSIDE = 8192
        w = hp.sphtfunc.pixwin(NSIDE)
        interp_func = interp1d(np.arange(len(w)), w)
        cl_pix = interp_func(ell_web)

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
        N = 10**6
        #plt.yscale("log")
        #plt.xscale("log")
        plt.xlim(0,8000)
        plt.ylim(0,3)
        plt.title("TSZ x Halos ("+z_m_cut_title+")")
        ax.plot(ell_web, N*ell_web*(ell_web+1)*cl_web/(2*np.pi*(cl_pix)), label = "Websky",c='b',lw=0.7)
        ax.plot(ell, (cl_1h + cl_2h), label = "Total",c='r',lw=0.7)
        ax.plot(ell, cl_1h, label = "1 Halo",c='yellow',lw=0.7)
        ax.plot(ell, cl_2h, label = "2 Halo",c='g',lw=0.7)
        ax.legend(loc=2,ncol = 1,frameon=False,fontsize=14)
        ax.set_xlabel(r"$\ell$",size=title_size)
        ax.set_ylabel(r"$\ell(\ell+1)C_\ell/2\pi\quad$",size=title_size)
        fig.tight_layout()
        plt.savefig("/moto/home/mr4449/class_sz_vs_websky_plots/class_vs_websky(tszxgal" + z_m_cut + ").png")
    
    if x == 3:

        cl_3h = np.load(path + "cl_" + sig_type[x] +"_3h" + z_m_cut + ".npy")

        NSIDE = 4096
        w = hp.sphtfunc.pixwin(NSIDE)
        interp_func = interp1d(np.arange(len(w)), w)
        cl_pix = interp_func(ell_web)

        #interp_func_class_sz = interp1d(ell, cl_1h + cl_2h + cl_3h)
        #print(ell[0])
        #print(ell[len(ell)-1])
        #print(len(ell))

        #cl_class = interp_func_class_sz(ell_web[11:7300])

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
        plt.title("KSZ x Halos ("+z_m_cut_title+")")
        #plt.yscale("log")
        #plt.xscale("log")
        plt.xlim(0,8000)
        #plt.ylim(0,3)
        ax.plot(ell_web, ell_web*(ell_web+1)*cl_web/(2*np.pi*cl_pix), label = "Websky",c='b',lw=0.7)
        #ax.plot(ell_web, ell_web*(ell_web+1)/(2*np.pi*cl_pix**2), label = "Pixwin",c='',lw=0.7)
        #ax.plot(ell_web, ell_web*(ell_web+1)*cl_web/(2*np.pi*cl_pix*cl_class), label = "Websky/Class_sz",c='pink',lw=0.7)
        ax.plot(ell, cl_1h + cl_2h + cl_3h, label = "Total",c='black',lw=0.7)
        ax.plot(ell, cl_1h, label = "1 Halo",c='r',lw=0.7)
        ax.plot(ell, cl_2h, label = "2 Halo",c='g',lw=0.7)
        ax.plot(ell, cl_3h, label = "3 Halo",c='yellow',lw=0.7)
        ax.legend(loc=2,ncol = 1,frameon=False,fontsize=14)
        ax.set_xlabel(r"$\ell$",size=title_size)
        ax.set_ylabel(r"$\ell(\ell+1)C_\ell/2\pi\quad$",size=title_size)
        fig.tight_layout()
        plt.savefig("/moto/home/mr4449/class_sz_vs_websky_plots/class_vs_websky(kszxgal" + z_m_cut + ").png")

#plot_z_m_cut(0.5, 1, 5e13, 5e15, 3)
#plot_z_m_cut(0.3, 0.5, 5e13, 5e15, 3)
#plot_z_m_cut(0.5, 0.7, 5e13, 5e15, 3)
#plot_z_m_cut(0.1, 0.3, 5e13, 5e15, 3)



#plot_z_m_cut(0.5, 1, 5e13, 5e15, 3)

#w = hp.sphtfunc.pixwin(8192)
#b = hp.sphtfunc.pixwin(4096)
#print(np.arange(len(b)))
#plt.plot(np.arange(len(w)), w**2, c = 'blue')
#plt.plot(np.arange(len(b)), b**2, c = 'red')
#plt.xlim(0,8000)
#plt.savefig("/moto/home/mr4449/pixwin.png")



z_min = 0.5
z_max = 1
m_min = 5e13
m_max = 5e15

M_min = Decimal(m_min)
M_max = Decimal(m_max)

z_m_cut = "_"+str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_m_"+'{:.2e}'.format(M_max)
z_m_cut_title = str(z_min)+"<z<"+str(z_max)+" & "+'{:.2e}'.format(M_min)+"<m<"+'{:.2e}'.format(M_max)

path = "/moto/hill/users/mr4449/websky_calcs/class_sz_calcs/"
cl_1h_xout_1 = np.load(path + "cl_ksz_1h_xout_1" + z_m_cut + ".npy")
cl_2h_xout_1 = np.load(path + "cl_ksz_2h_xout_1" + z_m_cut + ".npy" )
cl_3h_xout_1 = np.load(path + "cl_ksz_3h_xout_1" + z_m_cut + ".npy")
ell_x_out_1 = np.load(path + "ell_ksz_xout_1" + z_m_cut + ".npy")

cl_1h_xout_5 = np.load(path + "cl_ksz_1h_xout_5" + z_m_cut + ".npy")
cl_2h_xout_5 = np.load(path + "cl_ksz_2h_xout_5" + z_m_cut + ".npy")
cl_3h_xout_5 = np.load(path + "cl_ksz_3h_xout_5" + z_m_cut + ".npy")
ell_x_out_5 = np.load(path + "ell_ksz_xout_5" + z_m_cut + ".npy")

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
plt.title("kSZ x Halos ("+z_m_cut_title+")")
#plt.yscale("log")
#plt.xscale("log")
plt.xlim(0,8000)
#plt.ylim(0,3)
ax.plot(ell_x_out_1, cl_1h_xout_1 + cl_2h_xout_1 + cl_3h_xout_1, label = "xout 1",c='black',lw=0.7)
ax.plot(ell_x_out_5, cl_1h_xout_5 + cl_2h_xout_5 + cl_3h_xout_5, label = "xout 5",c='red',lw=0.7)
ax.legend(loc=2,ncol = 1,frameon=False,fontsize=14)
ax.set_xlabel(r"$\ell$",size=title_size)
ax.set_ylabel(r"$\ell(\ell+1)C_\ell/2\pi\quad$",size=title_size)
fig.tight_layout()
plt.savefig("ksz_test_xout.png")


