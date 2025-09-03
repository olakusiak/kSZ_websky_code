import sys, os
import healpy as hp
import numpy as np
from astropy.io import fits
from colossus.cosmology import cosmology as cosmology_colossus
from colossus.lss import mass_function
from colossus.halo import concentration
from colossus.halo import mass_defs
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u
import matplotlib.pyplot as pl
from scipy import interpolate as interp
from decimal import Decimal

#Load up the halo catalog and process them to get M200c. Then save it in hdf5 format:

websky_cosmo = {'Omega_M': 0.31, 'Omega_B': 0.049, 'Omega_L': 0.69, 
                                 'h': 0.68, 'sigma_8': 0.81, 'n_s':0.965}

astropy_cosmo = FlatLambdaCDM(H0=websky_cosmo['h']*100, Om0=websky_cosmo['Omega_M'])

h = websky_cosmo['h']
omegam = websky_cosmo['Omega_M']
rho = 2.775e11 * omegam * (h**2) # Msun/Mpc^3

f=open('/moto/hill/projects/halos.pksc')
N=np.fromfile(f,count=3,dtype=np.int32)[0]

catalog=np.fromfile(f,count=N*10,dtype=np.float32)
catalog=np.reshape(catalog,(N,10))

x  = catalog[:,0];  y = catalog[:,1];  z = catalog[:,2] # Mpc (comoving)
vx = catalog[:,3]; vy = catalog[:,4]; vz = catalog[:,5] # km/sec
R  = catalog[:,6] # Mpc

# convert to mass, comoving distance, radial velocity, redshfit, RA and DEc
M200m    = 4*np.pi/3.*rho*R**3        # this is M200m (mean density 200 times mean) in Msun
chi      = np.sqrt(x**2+y**2+z**2)    # Mpc
vrad     = (x*vx + y*vy + z*vz) / chi # km/sec

M200m_wh = M200m * h

z_array = np.linspace(0.0001, 10, 1000)
zofchi  = interp.interp1d(np.log(astropy_cosmo.comoving_distance(z_array).value), z_array, kind='cubic')

redshift = zofchi(np.log(chi))

# select halos in the redshift range 0<z<2.5
indsel = np.where((redshift>0.0) &  (redshift<2.5))[0]
M200m = M200m[indsel]
redshift = redshift[indsel]
vrad = vrad[indsel]
x = x[indsel]
y = y[indsel]
z = z[indsel]
print("Shivam's halo code, after 0 < z < 2.5 cut, number of halos " + str(len(x)))
"""
ra, dec  = hp.vec2ang(np.column_stack((x,y,z)), lonlat=True)



params = {'w0':-1.0 ,'flat': True, 'H0': 69.0, 'Om0': 0.31, 'Ob0': 0.049, 'sigma8':0.81 ,'ns': 0.965}
cosmology_colossus.addCosmology('myCosmo', params)
cosmology_colossus.setCosmology('myCosmo')
cosmo_colossus = cosmology_colossus.setCosmology('myCosmo')

zedges_array = np.linspace(0.0,2.5,101)
zmin_array = zedges_array[:-1]
zmax_array = zedges_array[1:]

M200c_wh = np.zeros(len(M200m_wh))
for ji in range(len(zmin_array)):
    zmin_ji = zmin_array[ji]
    zmax_ji = zmax_array[ji]
    indsel = np.where((redshift>zmin_ji) &  (redshift<zmax_ji))[0]
    mhalo_sim = M200m_wh[indsel]
    zhalo_sim = redshift[indsel]
    
    zmean = np.mean(zhalo_sim)
    zmin, zmax = np.min(zhalo_sim), np.max(zhalo_sim)

    c200m = concentration.concentration(mhalo_sim, '200m', zmean, model = 'bhattacharya13')
    M200c, R200c, c200c = mass_defs.changeMassDefinition(mhalo_sim, c200m, zmean, '200m', '200c')

    M200c_wh[indsel] = M200c


# We don't need all the halos for pasting, so let's select only those with M200m > 4e13 Msun/h
indsel2 = np.where((M200m_wh>4e13))[0]
hcat_new = fits.HDUList()
cols = []
cols.append(fits.Column(name='ra', format='D', array=ra[indsel2]))
cols.append(fits.Column(name='dec', format='D', array=dec[indsel2]))
cols.append(fits.Column(name='z', format='D', array=redshift[indsel2]))
cols.append(fits.Column(name='M200c_wh', format='D', array=M200c_wh[indsel2]))
cols.append(fits.Column(name='M200m_wh', format='D', array=M200m_wh[indsel2]))
cols.append(fits.Column(name='vrad', format='D', array=vrad[indsel2]))
hcat_new.append(fits.BinTableHDU.from_columns(cols))
hcat_new.writeto('/moto/hill/users/mr4449/shivam_maps/halo_input_wM200c_zmax2p5_minM200m_4e13_new.fits', overwrite=True)

#Now load up the saved halo catalog, original websky map and the new painted map. Convert it to some desirable nside

df = fits.open('/moto/hill/users/mr4449/shivam_maps/halo_input_wM200c_zmax2p5_minM200m_4e13_new.fits')
ra_all = df[1].data['ra']
dec_all = df[1].data['dec']
z_all = df[1].data['z']
M200c_all = df[1].data['M200c_wh']
M200m_all = df[1].data['M200m_wh']
vlos_all = df[1].data['vrad']
print(len(vlos_all))
print(len(M200m_all))
print(len(M200c_all))
print(len(z_all))
print(len(dec_all))
print(len(ra_all))

nside=4096
"""
# ymap_orig = hp.read_map('/moto/hill/projects/websky/tsz_8192_hp.fits')
# ymap_orig = hp.ud_grade(ymap_orig, nside_out=nside)

# ymap_paint = hp.read_map('/moto/hill/users/mr4449/shivam_maps/tSZ_sim_B12_testv8_nside_8192_final.pkl')
# ymap_paint = hp.ud_grade(ymap_paint, nside_out=nside)

#Subselect the halo sample for testing the cross correlation and create the density map

# M200m_min = 4e13
# M200m_max = 5e20

# z_min = 0
# z_max = 2.5

z_m_cut_title = str(z_min)+" < z < "+str(z_max)+ " & "+'{:.2e}'.format(M200m_min)+" < M$_\u2609$/h < "+'{:.2e}'.format(M200m_max)
indsel = np.where((M200m_all > M200m_min) & (M200m_all < M200m_max) & (z_all>z_min) & (z_all<z_max))[0]
ra_all_deltah = ra_all[indsel]
dec_all_deltah = dec_all[indsel]
z_all_deltah = z_all[indsel]
M200c_all_deltah = M200c_all[indsel]
vlos_all_deltah = vlos_all[indsel]

# z_array_edges = np.linspace(0.0, 3.0, 200)
# z_cens_hist_halos = 0.5*(z_array_edges[:-1] + z_array_edges[1:])
# nz_hist, _ = np.histogram(z_all_deltah, z_array_edges)
# nz_norm = np.trapz(nz_hist, z_cens_hist_halos) 
# nz_hist_halos = nz_hist/float(nz_norm)


# maph   = np.zeros((hp.nside2npix(nside)))
# pix = hp.ang2pix(nside, ra_all_deltah, dec_all_deltah, lonlat=True)
# weight = 1. #1 for number density, array of size(x) for arbitrary
# np.add.at(maph, pix, weight)
# delta_h = maph/np.mean(maph) - 1

# #Do the correlation

# Cl_hy_paint = hp.anafast(delta_h, np.array(ymap_paint))
# Cl_hy_orig = hp.anafast(delta_h, ymap_orig)

# #Plot it

# ell = np.arange(len(Cl_hy_paint))
# fac = ell*(ell+1)/(2*np.pi)
# pl.figure()
# pl.plot(ell, fac*Cl_hy_paint, label='painted')
# pl.plot(ell, fac*Cl_hy_orig, label='original', alpha=0.5)
# pl.grid()
# pl.legend(fontsize=15)
# pl.xlabel(r'$\ell$', fontsize=17)
# pl.ylabel(r'$\ell(\ell+1)C_{\ell}/2\pi$', fontsize=17)
# pl.title(r'1e14 $< M <$ 5e14, 1.4 $< z <$ 2.0', fontsize=17)
