import healpy as hp 
import numpy as np 
import matplotlib.pyplot as plt  
from math import pi as pi


path1 = "/moto/hill/projects/agora_sims/constr_maps/"
gal_files = ["gal_map_blue_dndz_8192_deCIBHOD"]

### GAUSS ###

"""
path2 = "/moto/hill/users/mr4449/ksz_sq_gauss_fil_maps/ksz"
ksz_sq_fil_maps = []
ksz_sq_name = []

k = 2300
while k <= 2500:
    homepath = path2 + "80_sq_gauss_" + str(k) + ".fits"
    ksz_sq_fil_maps.append(homepath)
    ksz_sq_name.append("_80_sq_gauss_" + str(k))
    k = k + 50
col = ["blue"]


### S4 ###

"""
path2 = "/moto/home/mr4449/"
ksz_sq_fil_maps = [path2+"ksz76_sq_s4_fil.fits", path2+"ksz78_sq_s4_fil.fits", path2+"ksz80_sq_s4_fil.fits"]
ksz_sq_name = ["_76_sq_s4", "_78_sq_s4", "_80_sq_s4"]
col = ["blue", "green", "red"]
"""

### SO ###
"""
path2 = "/moto/hill/users/mr4449/ksz_sq_fil_maps/ksz"
ksz_sq_fil_maps = [path2+"76_sq_so_fil.fits", path2+"78_sq_so_fil.fits", path2+"80_sq_so_fil.fits"]
ksz_sq_name = ["_76_sq_so", "_78_sq_so", "_80_sq_so"]
col = ["blue", "green", "red"]
"""
LMAX = 8000

i = 0
while i < len(gal_files):
    gal_f = gal_files[i]
    gal = hp.read_map(path1+gal_f+".fits")
    j = 0
    while j < len(ksz_sq_name):
        cl = hp.anafast(hp.read_map(ksz_sq_fil_maps[j]), gal, lmax=LMAX)
        ell = np.arange(len(cl))
        np.save("/moto/home/mr4449/cross_ksz_blue_gauss/cross_ksz_blue_gauss"+col[i]+ksz_sq_name[j], ell*(ell + 1)*cl/(2*pi))
        j = j + 1
    i = i + 1

