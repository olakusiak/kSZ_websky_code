import healpy as hp 
import numpy as np 
import matplotlib.pyplot as plt  
from math import pi as pi

tau_path = "/moto/hill/projects/agora_tau_unl/"
tau_maps = [tau_path+"agora_utauNG_bahamas76_bnd_unb_1.0e+12_1.0e+18.fits",
            tau_path+"agora_utauNG_bahamas78_bnd_unb_1.0e+12_1.0e+18.fits",
            tau_path+"agora_utauNG_bahamas80_bnd_unb_1.0e+12_1.0e+18.fits"]

gal_path = "/moto/hill/projects/agora_sims/constr_maps/"
gal_maps = [gal_path+"gal_map_blue_dndz_8192_deCIBHOD.fits",
            gal_path+"gal_map_red_dndz_8192_deCIBHOD.fits",
            gal_path+"gal_map_green_dndz_8192_deCIBHOD.fits"]

tsz_map = "/moto/hill/projects/agora_sims/mdpl2_ltszNG_bahamas80_rot_sum_4_176_bnd_unb_1.0e+12_1.0e+18_v103021_lmax24000_nside8192_interp1.0_method1_1_lensed_map.fits"
tsz = hp.read_map(tsz_map)
col = ["blue", "red", "green"]
num_map = ["76", "78", "80"]


LMAX = 8000
i = 0 
while i < len(tau_maps):
    tau = hp.read_map(tau_maps[i])
    cl = hp.anafast(tau, tsz, lmax=LMAX)
    ell = np.arange(len(cl))
    np.save("/moto/home/mr4449/cross_tau/cross__tsz_tau"+num_map[i], ell*(ell + 1)*cl/(2*pi))
    i = i + 1

"""

### Auto Cross Correlation ###

auto_tau_path = "/moto/home/mr4449/cross_tau/cross_tau_"
LMAX = 8000
i = 0
while i < len(tau_maps):
    tau = hp.read_map(tau_maps[i])
    cl = hp.anafast(tau)
    ell = np.arange(len(cl))
    np.save(auto_tau_path+num_map[i]+"_auto", ell*(ell + 1)*cl/(2*pi))
    plt.plot(ell, ell*(ell + 1)*cl/(2*pi), label = "Bahamas " + num_map[i])
    plt.title("Tau Auto Cross Spectrum ") 
    plt.xlabel("$\ell$")
    plt.ylabel("$\ell(\ell+1)C_{\ell}$/2π")
    plt.yscale("log")
    plt.savefig("/moto/home/mr4449/cross_tau/auto_power_spectrum_" + num_map[i])
    i = i + 1
"""
