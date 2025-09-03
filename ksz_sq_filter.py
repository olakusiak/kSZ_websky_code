import healpy as hp 
import numpy as np 
import math

"""
path1 = "/moto/hill/projects/agora_sims/mdpl2_lkszNG_bahamas"
path2 = "_rot_sum_4_176_bnd_unb_1.0e+12_1.0e+18_v103021_lmax24000_nside8192_interp1.6_method1_1_lensed_map.fits"
pathf = "/moto/home/mr4449/gauss_fils"
ksz_fits = [path1+str(76)+path2, path1+str(78)+path2, path1+str(80)+"_rot_sum_4_176_bnd_unb_1.0e+12_1.0e+18_v103021_lmax24000_nside8192_interp1.0_method1_1_lensed_map.fits"]
filters = [pathf+"gauss_300_fl.txt", pathf+"gauss_400_fl.txt",pathf+"gauss_500_fl.txt", pathf+"planck_fl_A_170422.txt"]
nummap = ["76", "78", "80"]
filname = ["gauss_300", "gauss_400", "gauss_500", "planck"]
"""

path1 = "/moto/hill/projects/agora_sims/mdpl2_lkszNG_bahamas"
path2 = "_rot_sum_4_176_bnd_unb_1.0e+12_1.0e+18_v103021_lmax24000_nside8192_interp1.6_method1_1_lensed_map.fits"
ksz_fits = [path1+str(80)+"_rot_sum_4_176_bnd_unb_1.0e+12_1.0e+18_v103021_lmax24000_nside8192_interp1.0_method1_1_lensed_map.fits"]

filters = []
filname = []
fils_data = []


k = 3050
while k <= 3500:
    homepath = '/moto/home/mr4449/gauss_fils'
    string = homepath + "/gauss_" + str(k) + "_fl.txt"
    filters.append(string)
    filname.append("gauss_" + str(k))
    k = k + 50

i = 0
while i < len(filters):
    fil = np.loadtxt(filters[i])
    filter_data = []
    j = 0
    while j < len(fil):
        filter_data.append(fil[j][1])
        j = j + 1
    fils_data.append(np.array(filter_data))
    i = i + 1
filters = fils_data

LMAX = 8000
i = 0
while i < len(ksz_fits):
    ksz_f = hp.read_map(ksz_fits[i])
    ksz_alm = hp.map2alm(ksz_f, lmax = LMAX)
    j = 0
    while j < len(filters):
        NSIDE = 8192
        fil = filters[j]
        ksz_alm_fil = hp.almxfl(ksz_alm, fil)
        ksz_fil_map = hp.alm2map(ksz_alm_fil, NSIDE)
        ksz_sq_fil = ksz_fil_map**2
        hp.write_map("/moto/hill/users/mr4449/ksz_sq_gauss_fil_maps/ksz80_sq_"+filname[j]+".fits", ksz_sq_fil, overwrite = True)
        j = j + 1
    i = i + 1
