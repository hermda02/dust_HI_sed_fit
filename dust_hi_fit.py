import numpy as np
import healpy as hp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
from scipy import stats

mpl.rcParams['text.usetex'] = True

# For beautiful maps
planck  = np.loadtxt('/home/daniel/GraduateSchool/MastersThesis/tools/Planck_color.txt')/255.
pl_cmap = colors.LinearSegmentedColormap.from_list('planck',planck)

#----------------------------------
# Read othe maps
#----------------------------------

raw_353   = hp.read_map('../dust_data/npipe6v20_353-5_bmap_QUADCOR_n0064_60arcmin_MJy_no_cmb.fits')#,verbose=False)
raw_545   = hp.read_map('../dust_data/npipe6v20_545-1_bmap_QUADCOR_n0064_60arcmin_MJy_no_cmb.fits')#,verbose=False)
raw_857   = hp.read_map('../dust_data/npipe6v20_857-1_bmap_QUADCOR_n0064_60arcmin_MJy_no_cmb.fits')#,verbose=False)

raw_240   = hp.read_map('../dust_data/DIRBE_240micron_1deg_h064_v2.fits')#,verbose=False)
raw_100   = hp.read_map('../dust_data/DIRBE_100micron_Nside064_60a.fits')#,verbose=False)

HI_temp   = hp.read_map('../dust_data/HI_vel_filter_60arcmin_0064.fits')#,verbose=False)
HI_mask   = hp.read_map('../dust_data/HI_mask.fits')#,verbose=False)

#----------------------------------

npix      = hp.nside2npix(64)

niter     = 500

h         = 6.626e-34     # Planck constant, J s
c         = 3.0e+8        # Speed of light, m/s
k         = 1.38e-23      # Boltzmann constant, J/K

def planck(fre,T):
	return ((2*h*fre**3)/c**2)*(1.0/(np.exp((h*fre)/(k*T))-1))

def sample_T(x,y,z,T,mu):

	thetas    = np.zeros(niter+1,dtype='f8')
	thetas[0] = T
	temp      = T

	for it in range(niter):
		r = np.random.normal(temp,mu)
		test  = np.zeros(len(freq),dtype='f8')
		for j in range(len(freq)):

			test[j] = amps[j]*z*planck(freq[j]*1e9,r)

		p = min(1,x/np.sum(abs(test-y)/y))

		if (np.random.uniform() < p):
			if (r > 10. and r < 30.):
				temp = r
		thetas[it+1] = temp

	# print(temp)
	return temp

def create_T():

	global new_T

	for i in range(npix):
		if (False == np.isnan(HI_mask[i])):

			y = np.empty(5,dtype='f8')
			y[0] = masked_353[i]
			y[1] = masked_545[i]
			y[2] = masked_857[i]
			y[3] = masked_240[i]
			y[4] = masked_100[i]

			z    = masked_HI[i] 

			# Computed a weighted difference between signal and HI model:
			x = (abs(HI_dust_353[i]-masked_353[i])/masked_353[i] + abs(HI_dust_545[i]-masked_545[i])/masked_545[i] + \
				abs(HI_dust_857[i]-masked_857[i])/masked_857[i] + abs(HI_dust_240[i]-masked_240[i])/masked_240[i] + \
				abs(HI_dust_100[i]-masked_100[i])/masked_100[i])

			new_T[i] = sample_T(x,y,z,new_T[i],dust_T_sigma)

def sample_A(T):

	global amps

	dummy_353 = masked_HI*planck(353*1e9,T)
	dummy_545 = masked_HI*planck(545*1e9,T)
	dummy_857 = masked_HI*planck(857*1e9,T)
	dummy_240 = masked_HI*planck(1240*1e9,T)
	dummy_100 = masked_HI*planck(2998*1e9,T)

	sum1      = np.empty(5,dtype='f8')
	sum2      = np.empty(5,dtype='f8')

	for i in range(npix):
		if (False == np.isnan(HI_mask[i])):

			# print(masked_353[i])

			sum1[0] = sum1[0] + (masked_353[i]*dummy_353[i])
			sum1[1] = sum1[1] + (masked_545[i]*dummy_545[i])
			sum1[2] = sum1[2] + (masked_857[i]*dummy_857[i])
			sum1[3] = sum1[3] + (masked_240[i]*dummy_240[i])
			sum1[4] = sum1[4] + (masked_100[i]*dummy_100[i])

			print(masked_353[i]*dummy_353[i])

			sum2[0] = sum2[0] + (dummy_353[i]**2.0)
			sum2[1] = sum2[1] + (dummy_545[i]**2.0)
			sum2[2] = sum2[2] + (dummy_857[i]**2.0)
			sum2[3] = sum2[3] + (dummy_240[i]**2.0)
			sum2[4] = sum2[4] + (dummy_100[i]**2.0)

	print(sum1[0])
	print(sum2[0])

	for i in range(5):
		amps[i] = sum1[i]/sum2[i]
		# if (sum1[i] < 0.0):
		# 	print('SUM1')
		# 	print(i)
		# 	print(sum1[i])
		# 	exit()
		# if (sum2[i] < 0.0):
		# 	print('SUM2')
		# 	print(sum2[i])
		# 	exit()
		print(amps[i])

	chisq = 0.0

	HI_dust_353 = dummy_353*amps[0]
	HI_dust_545 = dummy_545*amps[1]
	HI_dust_857 = dummy_857*amps[2]
	HI_dust_240 = dummy_240*amps[3]
	HI_dust_100 = dummy_100*amps[4]



freq = [353., 545., 857., 1240., 2998.]

#--------------------------------------------------|
# Changing units from MJy/sr to W sr^-1 m^-2 Hz^-1 |
#--------------------------------------------------|

masked_HI    = HI_temp
masked_353   = raw_353#*1e-20
masked_545   = raw_545#*1e-20
masked_857   = raw_857#*1e-20
masked_240   = raw_240#*1e-20
masked_100   = raw_100#*1e-20

HI_dust_353  = np.empty(npix,dtype='f8')
HI_dust_545  = np.empty(npix,dtype='f8')
HI_dust_857  = np.empty(npix,dtype='f8')
HI_dust_240  = np.empty(npix,dtype='f8')
HI_dust_100  = np.empty(npix,dtype='f8')

dummy_353    = np.empty(npix,dtype='f8')
dummy_545    = np.empty(npix,dtype='f8')
dummy_857    = np.empty(npix,dtype='f8')
dummy_240    = np.empty(npix,dtype='f8')
dummy_100    = np.empty(npix,dtype='f8')

new_T        = np.empty(npix,dtype='f8')

amps         = np.empty(5,dtype='f8')

pics         = 0

dust_T_init  = input('Dust temperature estimate: ')
dust_T_sigma = input('Std for dust temperature: ')


# Apply the same mask to all used arrays
for i in range(npix):
	if (True == np.isnan(HI_mask[i])):

		masked_HI[i]   = np.nan
		masked_353[i]  = np.nan
		masked_545[i]  = np.nan
		masked_857[i]  = np.nan
		masked_240[i]  = np.nan
		masked_100[i]  = np.nan
		HI_dust_353[i] = np.nan
		HI_dust_545[i] = np.nan
		HI_dust_857[i] = np.nan
		HI_dust_240[i] = np.nan
		HI_dust_100[i] = np.nan
		dummy_353[i]   = np.nan
		dummy_545[i]   = np.nan
		dummy_857[i]   = np.nan
		dummy_240[i]   = np.nan
		dummy_100[i]   = np.nan
		new_T[i]       = np.nan

	else:
		pics           = pics + 1
		new_T[i]       = dust_T_init

times = input("How many times would you like to sample A and T? ")

for i in range(times):

	print('Iteration:', i+1)

	print('------------------')
	sample_A(new_T)
	create_T()

	iterz = str(i+1)

	np.savetxt('amplitudes_'+iterz+'.dat',amps)
	file = 'dust_Td_'+iterz+'.fits'
	hp.write_map(file,new_T)

	print('------------------')

sample_A(new_T)
np.savetxt('amplitudes_final.dat',amps)

plt.scatter(freq,amps)
plt.ylim(1e-8,1e-5)
plt.yscale('log')
plt.show()

hp.mollview(new_T,min=15,max=25,cmap=pl_cmap)
plt.show()