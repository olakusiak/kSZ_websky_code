import healpy as hp
from   cosmology import *
import matplotlib
import matplotlib.pyplot as plt
from scipy import integrate 
import numpy as np
from decimal import Decimal
from classy_sz import Class
import os
from scipy import interpolate as interp
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u
import time

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
    


def z_and_mass_cut_cat(z_min, z_max, m_min, m_max, fil, lstar = 0, mh = False, cat_type = "org"):
    M_min = Decimal(m_min)
    M_max = Decimal(m_max)

    rho = 2.775e11*omegam*h**2 # Msun/Mpc^3
    if mh == True:
        z_m_cut = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_mh_"+'{:.2e}'.format(M_max)
    else:
	    z_m_cut = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_m_"+'{:.2e}'.format(M_max)

    path = "/moto/hill/users/mr4449/websky_calcs/x_halo_"+cat_type+"_"+z_m_cut+".txt"

    if os.path.isfile(path) == False:
        print(path + " does not exist")

        if cat_type == "lambda":
            f=open('/moto/hill/projects/halos_lambda.pksc')
        
        if cat_type == "toronto":
            f=open('/moto/hill/projects/halos_toronto.pksc')
        
        if cat_type == "org":
            f=open('/moto/hill/projects/halos.pksc')
            
        N=np.fromfile(f,count=3,dtype=np.int32)[0]


        # only take first five entries for testing (there are ~8e8 halos total...)
        # comment the following line to read in all halos

        catalog=np.fromfile(f,count=N*10,dtype=np.float32)
        catalog=np.reshape(catalog,(N,10))

        x  = catalog[:,0];  y = catalog[:,1];  z = catalog[:,2] # Mpc (comoving)
        vx = catalog[:,3]; vy = catalog[:,4]; vz = catalog[:,5] # km/sec
        R  = catalog[:,6] # Mpc

        # convert to mass, comoving distance, radial velocity, redshfit, RA and DEc
        #M200m    = 4*np.pi/3.*rho*R**3   
        if mh == True:
            M200m    = 4*np.pi/3.*rho*R**3*h  # this is M200m (mean density 200 times mean) in Msun
        else:
            M200m    = 4*np.pi/3.*rho*R**3  # this is M200m (mean density 200 times mean) in Msun
        chi	 = np.sqrt(x**2+y**2+z**2)    # Mpc
        vrad     = (x*vx + y*vy + z*vz) / chi # km/sec
        redshift = zofchi(chi)     
        print("Largest Halo Mass "  + str(np.max(M200m)))

        x_m = x[np.where(M200m < m_max)]
        y_m = y[np.where(M200m < m_max)]
        z_m = z[np.where(M200m < m_max)]
        redshift_m = redshift[np.where(M200m < m_max)]
        M200m_M = M200m[np.where(M200m < m_max)]

        x_m2 = x_m[np.where(M200m_M > m_min)]
        y_m2 = y_m[np.where(M200m_M > m_min)]
        z_m2 = z_m[np.where(M200m_M > m_min)]
        redshift_m2 = redshift_m[np.where(M200m_M > m_min)]
        M200m_M2 = M200m_M[np.where(M200m_M > m_min)]

        x_red = x_m2[np.where(redshift_m2 < z_max)]
        y_red = y_m2[np.where(redshift_m2 < z_max)]
        z_red = z_m2[np.where(redshift_m2 < z_max)]
        redshift_red = redshift_m2[np.where(redshift_m2 < z_max)]
        M200m_red = M200m_M2[np.where(redshift_m2 < z_max)]

        x_red2 = x_red[np.where(redshift_red > z_min)]
        y_red2 = y_red[np.where(redshift_red > z_min)]
        z_red2 = z_red[np.where(redshift_red > z_min)]
        redshift_red2 = redshift_red[np.where(redshift_red > z_min)]
        M200m_red2 = M200m_red[np.where(redshift_red > z_min)]

        hist_z, bins_z = np.histogram(redshift_red2, bins = 500)
        # hist_m, bins_m = np.histogram(M200m_red2, bins = 500)

        # bins_m = bins_m[0:len(bins_m)-1]
        # print("length of bins " + str(len(bins_m)))
        # print("length of hist " + str(len(hist_m)))
        # bins_m = np.insert(bins_m, 0, m_min)
        # bins_m = np.insert(bins_m, len(bins_m), m_max)
        # hist_m = np.insert(hist_m, 0, 0)
        # hist_m = np.insert(hist_m, len(hist_m), 0)
        #data_m = np.column_stack([bins_m, hist_m])
        
        N = integrate.simpson(hist_z, bins_z[0:len(bins_z)-1])
        hist_new = hist_z/N 
        bins_z = bins_z[0:len(bins_z)-1]
        bins_z = np.insert(bins_z, 0, z_min)
        bins_z = np.insert(bins_z, len(bins_z), z_max)
        hist_new = np.insert(hist_new, 0, 0)
        hist_new = np.insert(hist_new, len(hist_new), 0)
        data_z = np.column_stack([bins_z, hist_new])


        np.savetxt("/moto/hill/users/mr4449/websky_calcs/norm_websky_"+cat_type+"_dndz_"+z_m_cut+ ".txt", data_z)
        #np.savetxt("/moto/hill/users/mr4449/websky_calcs/websky_dndm_"+z_m_cut+ ".txt", data_m)
        np.save("/moto/hill/users/mr4449/websky_calcs/x_halo_"+cat_type+"_"+z_m_cut, x_red2)
        np.save("/moto/hill/users/mr4449/websky_calcs/y_halo_"+cat_type+"_"+z_m_cut, y_red2)
        np.save("/moto/hill/users/mr4449/websky_calcs/z_halo_"+cat_type+"_"+z_m_cut, z_red2)
        np.save("/moto/hill/users/mr4449/websky_calcs/mass_halo_"+cat_type+"_"+z_m_cut, M200m_red2)
        np.save("/moto/hill/users/mr4449/websky_calcs/redshift_halo_"+cat_type+"_"+z_m_cut, redshift_red2)
        x, y, z = x_red2, y_red2, z_red2

    else:
        print(path + " exists")
        x = np.load("/moto/hill/users/mr4449/websky_calcs/x_halo_"+cat_type+"_"+z_m_cut+".npy")
        y = np.load("/moto/hill/users/mr4449/websky_calcs/y_halo_"+cat_type+"_"+z_m_cut+ ".npy")
        z = np.load("/moto/hill/users/mr4449/websky_calcs/z_halo_"+cat_type+"_"+z_m_cut+ ".npy")

    theta, phi  = hp.vec2ang(np.column_stack((x,y,z)))
    path_web = "/moto/hill/users/mr4449/shivam_maps/"
    #path_web = "/moto/hill/projects/websky/"

    #websky_tsz = hp.read_map(path_web + "tsz_8192.fits")
    #websky_kap = hp.read_map(path_web + "kap.fits")
    websky_tsz = hp.read_map(path_web + "tSZ_sim_B12_testv8_nside_.fits")
    #websky_tsz = hp.read_map(path_web + "tau_map_sim_B12_testv6_nside_4096.fits")

    nside_tsz = hp.get_nside(websky_tsz)
    # nside_kap = hp.get_nside(websky_kap)
    # nside_ksz = hp.get_nside(websky_ksz)

    halo_tsz = np.zeros((hp.nside2npix(nside_tsz)))
    # halo_kap = np.zeros((hp.nside2npix(nside_kap)))
    # halo_ksz = np.zeros((hp.nside2npix(nside_ksz)))

    pix_tsz = hp.ang2pix(nside_tsz, theta, phi)
    # pix_kap = hp.ang2pix(nside_kap, theta, phi)
    #pix_ksz = hp.ang2pix(nside_ksz, theta, phi)
    
    weight = 1.

    np.add.at(halo_tsz, pix_tsz, weight)
    # np.add.at(halo_kap, pix_kap, weight)
    # np.add.at(halo_ksz, pix_ksz, weight)

    m_tsz = np.mean(halo_tsz)
    # m_kap = np.mean(halo_kap)
    # m_ksz = np.mean(halo_ksz)

    halo_tsz = (halo_tsz - m_tsz)/m_tsz
    # halo_kap = (halo_kap - m_kap)/m_kap
    # halo_ksz = (halo_ksz - m_ksz)/m_ksz

    cl_tsz = hp.anafast(halo_tsz, websky_tsz, lmax = 10000)
    cl = cl_tsz
    # cl_kap = hp.anafast(halo_kap, websky_kap)
    # cl_halo = hp.anafast(halo_tsz)

    # if fil == "gauss_" + str(lstar):
    #     fl = g_fil(lstar)
        
    # else:
    #     if fil == "pl":
    #         path_fil = "/moto/home/mr4449/planck_fl_A_170422.txt"
    #     if fil == "so":
    #         path_fil = "/moto/home/mr4449/so_fl_A_170422.txt"
    #     if fil == "cmbs4":
    #         path_fil = "/insomnia001/home/mr4449/class_sz/class_sz_auxiliary_files/s4_fl_A_170422.txt"

    #     f = np.loadtxt(path_fil)
    #     i = 0
    #     fl = [0,0]
    #     while i < len(f):
    #         fl.append(f[i][1])
    #         i = i + 1

    # alm_ksz = hp.map2alm(websky_ksz)
    # alm_ksz_fl = hp.almxfl(alm_ksz, fl)
    # websky_map_ksz = (hp.alm2map(alm_ksz_fl, nside = nside_ksz))**2
    #cl_ksz = hp.anafast(halo_ksz, websky_ksz, lmax = 10000)

    #np.save("/moto/hill/users/mr4449/websky_calcs/cl_ksz_halo_shiv_v5_4096_"+fil+"_"+z_m_cut,cl)
    #np.save("/moto/hill/users/mr4449/websky_calcs/cl_halo_halo_8192_"+z_m_cut,cl)
    #np.save("/moto/hill/users/mr4449/websky_calcs/cl_tau_halo_shiv_v5_4096_"+z_m_cut,cl)
    np.save("/moto/hill/users/mr4449/websky_calcs/cl_tsz_halo_shiv_4096_"+cat_type+"_"+z_m_cut,cl)
    # np.save("/insomnia001/depts/hill/users/mr4449/websky_calcs/cl_tsz_halo_"+z_m_cut,cl_tsz)
    # np.save("/insomnia001/depts/hill/users/mr4449/websky_calcs/cl_kap_halo_"+z_m_cut,cl_kap)
    # np.save("/insomnia001/depts/hill/users/mr4449/websky_calcs/cl_auto_halo_"+z_m_cut,cl_halo)

    print("Cuts and Cross Correlations Made for " + cat_type + " & " + z_m_cut + " & " + fil)

def z_and_mass_cut(z_min, z_max, m_min, m_max, fil, lstar = 0, mh = False):
    M_min = Decimal(m_min)
    M_max = Decimal(m_max)

    #rho = 2.775e11*omegam*h**2 # Msun/Mpc^3
    if mh == True:
        z_m_cut = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_mh_"+'{:.2e}'.format(M_max)
    else:
	    z_m_cut = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_m_"+'{:.2e}'.format(M_max)

    path = "/moto/hill/users/mr4449/websky_calcs/x_halo_new_cut_"+z_m_cut+".txt"

    if os.path.isfile(path) == False:
        websky_cosmo = {'Omega_M': 0.31, 'Omega_B': 0.049, 'Omega_L': 0.69, 
                                 'h': 0.68, 'sigma_8': 0.81, 'n_s':0.965}

        astropy_cosmo = FlatLambdaCDM(H0=websky_cosmo['h']*100, Om0=websky_cosmo['Omega_M'])

        h = websky_cosmo['h']
        omegam = websky_cosmo['Omega_M']
        rho = 2.775e11 * omegam * (h**2) # Msun/Mpc^3

        print(path + " does not exist")
        
        f=open('/moto/hill/projects/halos.pksc')
            
        N=np.fromfile(f,count=3,dtype=np.int32)[0]


        # only take first five entries for testing (there are ~8e8 halos total...)
        # comment the following line to read in all halos

        catalog=np.fromfile(f,count=N*10,dtype=np.float32)
        catalog=np.reshape(catalog,(N,10))

        x  = catalog[:,0];  y = catalog[:,1];  z = catalog[:,2] # Mpc (comoving)
        vx = catalog[:,3]; vy = catalog[:,4]; vz = catalog[:,5] # km/sec
        R  = catalog[:,6] # Mpc
        
        
        # convert to mass, comoving distance, radial velocity, redshfit, RA and DEc
        if mh == True:
            M200m    = 4*np.pi/3.*rho*R**3*h  # this is M200m (mean density 200 times mean) in Msun
        else:
            M200m    = 4*np.pi/3.*rho*R**3  # this is M200m (mean density 200 times mean) in Msun
        chi	 = np.sqrt(x**2+y**2+z**2)    # Mpc
        vrad     = (x*vx + y*vy + z*vz) / chi # km/sec
        #redshift = zofchi(chi)     
        z_array = np.linspace(0.0001, 10, 1000)
        zofchi  = interp.interp1d(np.log(astropy_cosmo.comoving_distance(z_array).value), z_array, kind='cubic')

        redshift = zofchi(np.log(chi))
        #print("Largest Halo Mass "  + str(np.max(M200m)))
        indsel = np.where((M200m > m_min) & (M200m < m_max) & (redshift > z_min) & (redshift < z_max))[0]
        M200m = M200m[indsel]
        redshift = redshift[indsel]
        vrad = vrad[indsel]
        x = x[indsel]
        y = y[indsel]
        z = z[indsel]
        print(len(x))

        hist_z, bins_z = np.histogram(redshift, bins = 500)
        hist_m, bins_m = np.histogram(M200m, bins = 500)

        bins_m = bins_m[0:len(bins_m)-1]
        print("length of bins " + str(len(bins_m)))
        print("length of hist " + str(len(hist_m)))
        bins_m = np.insert(bins_m, 0, m_min)
        bins_m = np.insert(bins_m, len(bins_m), m_max)
        hist_m = np.insert(hist_m, 0, 0)
        hist_m = np.insert(hist_m, len(hist_m), 0)
        data_m = np.column_stack([bins_m, hist_m])
        
        N = integrate.simpson(hist_z, bins_z[0:len(bins_z)-1])
        hist_new = hist_z/N 
        bins_z = bins_z[0:len(bins_z)-1]
        bins_z = np.insert(bins_z, 0, z_min)
        bins_z = np.insert(bins_z, len(bins_z), z_max)
        hist_new = np.insert(hist_new, 0, 0)
        hist_new = np.insert(hist_new, len(hist_new), 0)
        data_z = np.column_stack([bins_z, hist_new])


        np.savetxt("/moto/hill/users/mr4449/websky_calcs/norm_websky_dndz_new_cut_"+z_m_cut+ ".txt", data_z)
        np.savetxt("/moto/hill/users/mr4449/websky_calcs/websky_dndm_new_cut_"+z_m_cut+ ".txt", data_m)
        np.save("/moto/hill/users/mr4449/websky_calcs/x_halo_new_cut_"+z_m_cut, x)
        np.save("/moto/hill/users/mr4449/websky_calcs/y_halo_new_cut_"+z_m_cut, y)
        np.save("/moto/hill/users/mr4449/websky_calcs/z_halo_new_cut_"+z_m_cut, z)
        np.save("/moto/hill/users/mr4449/websky_calcs/mass_halo_new_cut_"+z_m_cut, M200m)
        np.save("/moto/hill/users/mr4449/websky_calcs/redshift_halo_new_cut_"+z_m_cut, redshift)

    else:
        print(path + " exists")
        x = np.load("/moto/hill/users/mr4449/websky_calcs/x_halo_new_cut_"+z_m_cut+".npy")
        y = np.load("/moto/hill/users/mr4449/websky_calcs/y_halo_new_cut_"+z_m_cut+ ".npy")
        z = np.load("/moto/hill/users/mr4449/websky_calcs/z_halo_new_cut_"+z_m_cut+ ".npy")
    
    theta, phi  = hp.vec2ang(np.column_stack((x,y,z)))
    path_web = "/moto/hill/users/mr4449/shivam_maps/"
    #path_web = "/moto/hill/projects/websky/"

    #websky_tsz = hp.read_map(path_web + "tsz_8192.fits")
    #websky_tsz = hp.read_map(path_web + "tSZ_sim_B12_testv8_nside_8192.fits")
    websky_tsz = hp.read_map(path_web + "kSZ_sim_B12_testv4_nside_8192.fits")

    nside_tsz = hp.get_nside(websky_tsz)
    # nside_kap = hp.get_nside(websky_kap)
    # nside_ksz = hp.get_nside(websky_ksz)

    halo_tsz = np.zeros((hp.nside2npix(nside_tsz)))
    # halo_kap = np.zeros((hp.nside2npix(nside_kap)))
    # halo_ksz = np.zeros((hp.nside2npix(nside_ksz)))

    pix_tsz = hp.ang2pix(nside_tsz, theta, phi)
    # pix_kap = hp.ang2pix(nside_kap, theta, phi)
    # pix_ksz = hp.ang2pix(nside_ksz, theta, phi)
    
    weight = 1.

    np.add.at(halo_tsz, pix_tsz, weight)
    # np.add.at(halo_kap, pix_kap, weight)
    # np.add.at(halo_ksz, pix_ksz, weight)

    m_tsz = np.mean(halo_tsz)
    # m_kap = np.mean(halo_kap)
    # m_ksz = np.mean(halo_ksz)

    halo_tsz = (halo_tsz - m_tsz)/m_tsz
    # halo_kap = (halo_kap - m_kap)/m_kap
    # halo_ksz = (halo_ksz - m_ksz)/m_ksz

    # cl_tsz = hp.anafast(halo_tsz, websky_tsz, lmax = 10000)
    # cl = cl_tsz
    # print(cl)
    # cl_kap = hp.anafast(halo_kap, websky_kap)
    # cl_halo = hp.anafast(halo_tsz)

    if fil == "gauss_" + str(lstar):
        fl = g_fil(lstar)
        
    else:
        if fil == "pl":
            path_fil = "/moto/home/mr4449/planck_fl_A_170422.txt"
        if fil == "so":
            path_fil = "/moto/home/mr4449/so_fl_A_170422.txt"
        if fil == "cmbs4":
            path_fil = "/insomnia001/home/mr4449/class_sz/class_sz_auxiliary_files/s4_fl_A_170422.txt"

        f = np.loadtxt(path_fil)
        i = 0
        fl = [0,0]
        while i < len(f):
            fl.append(f[i][1])
            i = i + 1

    alm_ksz = hp.map2alm(websky_tsz)
    alm_ksz_fl = hp.almxfl(alm_ksz, fl)
    websky_map_ksz = (hp.alm2map(alm_ksz_fl, nside = nside_tsz))**2
    cl_ksz = hp.anafast(halo_tsz, websky_map_ksz, lmax = 10000)
    cl = cl_ksz

    #np.save("/moto/hill/users/mr4449/websky_calcs/cl_ksz_halo_shiv_v5_4096_"+fil+"_"+z_m_cut,cl)
    #np.save("/moto/hill/users/mr4449/websky_calcs/cl_halo_halo_8192_"+z_m_cut,cl)
    #np.save("/moto/hill/users/mr4449/websky_calcs/cl_tau_halo_shiv_v5_4096_"+z_m_cut,cl)
    np.save("/moto/hill/users/mr4449/websky_calcs/cl_ksz_halo_shiv_v4_8192_new_cut_"+fil+"_"+z_m_cut,cl)
    # np.save("/insomnia001/depts/hill/users/mr4449/websky_calcs/cl_tsz_halo_"+z_m_cut,cl_tsz)
    # np.save("/insomnia001/depts/hill/users/mr4449/websky_calcs/cl_kap_halo_"+z_m_cut,cl_kap)
    # np.save("/insomnia001/depts/hill/users/mr4449/websky_calcs/cl_auto_halo_"+z_m_cut,cl_halo)
    
    print("Cuts and Cross Correlations Made for " + z_m_cut + " & " + fil)

def z_and_mass_cut_map(z_min, z_max, m_min, m_max, num, mh = False):
    M_min = Decimal(m_min)
    M_max = Decimal(m_max)

    rho = 2.775e11*omegam*h**2 # Msun/Mpc^3
    if mh == True:
        z_m_cut = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_mh_"+'{:.2e}'.format(M_max)
    else:
	    z_m_cut = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_m_"+'{:.2e}'.format(M_max)


    rho = 2.775e11*omegam*h**2 # Msun/Mpc^3

    f=open('/moto/hill/projects/halos.pksc')
    N=np.fromfile(f,count=3,dtype=np.int32)[0]


    # only take first five entries for testing (there are ~8e8 halos total...)
    # comment the following line to read in all halos

    catalog=np.fromfile(f,count=N*10,dtype=np.float32)
    catalog=np.reshape(catalog,(N,10))

    x  = catalog[:,0];  y = catalog[:,1];  z = catalog[:,2] # Mpc (comoving)
    vx = catalog[:,3]; vy = catalog[:,4]; vz = catalog[:,5] # km/sec
    R  = catalog[:,6] # Mpc

    # convert to mass, comoving distance, radial velocity, redshfit, RA and DEc
    if mh == True:
        M200m    = h*4*np.pi/3.*rho*R**3  # this is M200m (mean density 200 times mean) in Msun
    else:
        M200m    = 4*np.pi/3.*rho*R**3  # this is M200m (mean density 200 times mean) in Msun     
    chi      = np.sqrt(x**2+y**2+z**2)    # Mpc
    vrad     = (x*vx + y*vy + z*vz) / chi # km/sec
    redshift = zofchi(chi)     

    x_m = x[np.where(M200m < m_max)]
    y_m = y[np.where(M200m < m_max)]
    z_m = z[np.where(M200m < m_max)]
    redshift_m = redshift[np.where(M200m < m_max)]
    M200m_M = M200m[np.where(M200m < m_max)]

    x_m2 = x_m[np.where(M200m_M > m_min)]
    y_m2 = y_m[np.where(M200m_M > m_min)]
    z_m2 = z_m[np.where(M200m_M > m_min)]
    redshift_m2 = redshift_m[np.where(M200m_M > m_min)]
    M200m_M2 = M200m_M[np.where(M200m_M > m_min)]

    x_red = x_m2[np.where(redshift_m2 < z_max)]
    y_red = y_m2[np.where(redshift_m2 < z_max)]
    z_red = z_m2[np.where(redshift_m2 < z_max)]
    redshift_red = redshift_m2[np.where(redshift_m2 < z_max)]
    M200m_red = M200m_M2[np.where(redshift_m2 < z_max)]

    x_red2 = x_red[np.where(redshift_red > z_min)]
    y_red2 = y_red[np.where(redshift_red > z_min)]
    z_red2 = z_red[np.where(redshift_red > z_min)]
    redshift_red2 = redshift_red[np.where(redshift_red > z_min)]
    M200m_red2 = M200m_red[np.where(redshift_red > z_min)]

    x_max = x_red2[np.argsort(M200m_red2)[-1*(num)]]
    y_max = y_red2[np.argsort(M200m_red2)[-1*(num)]]
    Z_max = z_red2[np.argsort(M200m_red2)[-1*(num)]]


    theta, phi  = hp.vec2ang(np.column_stack((x_max,y_max,Z_max)))
    print("theta = " + str(theta))
    print("phi = " + str(phi))

    path_web = "/moto/hill/users/mr4449/shivam_maps/tau_map_sim_B12_testv4_nside_4096.fits"
    tsz_map = hp.read_map(path_web)
    nside = hp.get_nside(tsz_map)
    halos_map = np.zeros((hp.nside2npix(nside))) 
    weight = 1.
    pix_halos = hp.ang2pix(nside, theta, phi)
    np.add.at(halos_map, pix_halos, weight)
    M_min = Decimal(m_min)
    M_max = Decimal(m_max)
    #hp.gnomview(halos_map, rot = [np.degrees(phi), 90 - np.degrees(theta)], reso = 0.1, min = 0, max = 1, xsize = 50, ysize = 50,
     #title = "Most Massive Halo ("+str(z_min)+"<z<"+str(z_max)+")")
    #hp.graticule()
    #plt.savefig("/moto/home/mr4449/" + str(num) + "_most_massive_halo_"+str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_m_"+'{:.2e}'.format(M_max)+".png")
    #plt.close()
    hp.gnomview(tsz_map, rot = [np.degrees(phi), 90 - np.degrees(theta)], reso = 0.1, xsize = 50, ysize = 50,
    title = "Tau ("+str(z_min)+"<z<"+str(z_max)+")")
    hp.graticule()
    plt.savefig("/moto/home/mr4449/shivam_plots/tau_shiv_v4_" + str(num) + "_most_massive_halo_"+z_m_cut+".png")
    plt.close()

start_time = time.time()
#LSTAR = 2000
z_and_mass_cut(0.2, 0.5, 6e13, 5e15, fil = "so", mh = True) 

# path_web = "/moto/hill/users/mr4449/shivam_maps/"
# websky_ksz = hp.read_map(path_web + "tau_map_sim_B12_testv4_nside_8192.fits")
# cl = hp.anafast(websky_ksz)
# np.save("/moto/hill/users/mr4449/websky_calcs/cl_tau_tau_shiv_v4_8192.npy",cl) 
print(str((time.time() - start_time)/60) + " minutes")
print(str(time.time() - start_time) + " seconds")

# MAP = hp.read_map("/moto/hill/users/mr4449/shivam_maps/tau_sim_B16_testv8_nside_8192.fits")
# hp.mollview(MAP, min = -10e-3 , max = 10e-3)
# plt.savefig("/moto/home/mr4449/shivam_plots/tsz_8192v8_.png")

#z_and_mass_cut(1, 2.5, 4e13, 1e15, fil = "pl", mh = True) 
# z_min = 0.1
# z_max = 2.5
# m_min = 6e13
# m_max = 5e15
# M_min = Decimal(m_min)
# M_max = Decimal(m_max)
# z_m_cut = str(z_min)+"_z_"+str(z_max)+"_&_"+'{:.2e}'.format(M_min)+"_mh_"+'{:.2e}'.format(M_max)
# path_1 = "/moto/hill/users/mr4449/shivam_maps/tau_sim_B16_testv8_nside_8192.fits"
# path_2 = "/moto/hill/users/mr4449/shivam_maps/tau_sim_B16_testv10_nside_8192.fits"
# map_v8 = hp.read_map(path_1)
# map_v10 = hp.read_map(path_2)
# hp.mollview(map_v8, min = 0, max = 0.0001)
# plt.savefig("pasted_tau_v8.png")
# cl_1 = np.load(path_1)
# cl_2 = np.load(path_2)
# ell_1 = np.arange(len(cl_1))
# ell_2 = np.arange(len(cl_2))
# #plt.plot(ell_1, ell_1*(ell_1+1)*cl_1/(2*np.pi), label = "v 10")
# plt.plot(ell_2, ell_2*(ell_2+1)*cl_2/(2*np.pi), label = "v 8")
# plt.legend()
# plt.savefig("tau_tau_shiv_v8_test.png")