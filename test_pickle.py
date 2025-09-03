import healpy as hp
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
from classy_sz import Class
import math
from decimal import Decimal


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
Z = 1
M = 5e15

r = np.array(data["r_array"])
m = np.array(data["M_array"])
z = np.array(data["z_array"])
a = np.array(data["scale_fac_array"])
r200c = np.array(data["r200c_mat"])
rho = data["rho_fit_mat"]

z_ind = np.where((z <= Z + 0.1) & (z >= Z - 0.1))[0]
z_ind = z_ind[int(len(z_ind)/2)]

m_ind = np.where((m <= M*(1 + 1/5)) & (m >= M*(1 - 1/5)))[0]
m_ind = m_ind[int(len(m_ind)/2)]

rho = rho[0:len(r), z_ind, m_ind]
r200c = r200c[m_ind, 0:len(r)]
a = a[z_ind]
x = a*r/r200c

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

websky_h = 0.68

def gas_param(A_0, A_m, A_z):
    return A_0*((M/(1e14*websky_h))**A_m)*(1+Z)**A_z

alpha_density = gas_param(A_alpha, alpha_m_alpha, alpha_z_alpha)
beta_density = gas_param(A_beta, alpha_m_beta, alpha_z_beta)
rho_density = gas_param(A_rho0, alpha_m_rho0, alpha_z_rho0)

def rho_fit(x):
    return rho_density*((x/xc)**gamma)*(1+(x/xc)**alpha_density)**(-(beta_density+gamma)/alpha_density)
rho_prof = rho_fit(x)



plt.xlim(1e-1, 1e1)
plt.xscale("log")
plt.plot(x,rho*(x**2), label = "Pasted Profiles")
plt.plot(x,rho_prof*(x**2), label = "Class SZ/ Known Form")
plt.legend()

plt.savefig("density_prof_shiv_test_3.png")
