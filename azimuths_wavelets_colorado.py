### This script performs wavelet transform using the mlpy package on 
### Azimuths of the Colorado river as a function of downstream distance

import numpy as np
import matplotlib.pyplot as plt
import mlpy.wavelet as wave
import math as math

print("##############################################")
print("Performing Wavelet Analysis on Colorado River")
print("##############################################")

      
print("### Initialising ###")

      
# Set the x step
dt=2000. # metres

# set wavelet function
wf = 'dog' # Derivative Of Gaussian [woof]
#wf = 'morlet'
#wf = 'paul'

# set morlet frequency (p = omega_0, e.g. 6), paul order (p = m, e.g. 2) or dog derivative (p = m, e.g. 2)
p = 6
dj = 0.05

# Check parameters
if wf == ('morlet') :
 if p != 6:
  print('ERROR P NOT EQUAL TO 6, ICWT WILL LOOK INCORRECT. CONTINUING ANYWAY...')
 recon_factor = (dj*(dt**0.5))/(0.776*(math.pi**-0.25))  # morlet: equation 11 in Torrance and Compo
 
if wf == ('dog') :
 if p == 2 :
  recon_factor = (dj*(dt**0.5))/(3.541*0.867)  # dog: equation 11 in Torrance and Compo
 if p == 6:
  recon_factor = (dj*(dt**0.5))/(1.966*0.884)
 if (p != 2) and (p !=6) :
  print('ERROR P NOT EQUAL TO 2 or 6, ICWT WILL LOOK INCORRECT. CONTINUING ANYWAY...')
  recon_factor = (dj*(dt**0.5))/(3.541*0.867)
  print('go')

# load in data
river = np.loadtxt('colorado.dat') # long, lat, x, theta
theta = river[:,3]
x = river[:,2] # in km's

print("### Performing wavelet transform ### ")

# Transform data by raising to  complex exponent
cmplx_thet = np.exp(2*math.pi*(theta)*complex(0,1)/360) 
cmplx_mean = np.mean(cmplx_thet)
cmplx_thet = cmplx_thet - np.mean(cmplx_thet) # subtract mean

# padding
[cmplx_thet_pad, cmplx_thet_pad_orig] = wave.pad(x,method='zeros')
cmplx_thet_pad = cmplx_thet # no padding

# calculate wavelet scales
scales = wave.autoscales(N=cmplx_thet_pad.shape[0], dt=dt, dj=dj, wf=wf, p=p)
# convert scales to fourier periods
period = wave.fourier_from_scales(scales,wf=wf,p=p)
# perform continuous wavelet transform
transformed = wave.cwt(x=cmplx_thet_pad, dt=dt, scales=scales, wf=wf, p=p)

# Calculate bearing as function of (k,x)
bearings = ((np.angle(transformed+cmplx_mean))*360/(2*math.pi))%360

# Calculate power as function of (k,x)
power = np.abs(transformed+cmplx_mean)**2
# Set axes for figures 
x_axis = x
y_axis = period # y_axis is wavenumbers units (m-1)

print("### Visualising results ###")
      
# Plot the scalogram of bearings
plt.pcolormesh(x_axis,y_axis,bearings,cmap='twilight')
plt.yscale('log')
plt.title("Scalogram of phase (as bearing) of complex transform for Colorado River")
plt.xlabel("Distance/km")
plt.ylabel("Wavelength/m")
plt.colorbar()
plt.show()

# Plot the scalogram of power
plt.pcolormesh(x_axis,y_axis,power,cmap='viridis')
plt.yscale('log')
plt.title("Scalogram of power of complex transform for Colorado River")
plt.xlabel("Distance/km")
plt.ylabel("Wavelength/m")
plt.colorbar()
plt.show()

print("### Filtering and reconstructing original signal ###")

     
reconstructed = wave.icwt(transformed, dt=dt, scales=scales, wf=wf, p=p)*recon_factor + cmplx_mean
reconstructed_1000 = wave.icwt(transformed[period>1000000,:], dt=dt, scales=scales[period>1000000], wf=wf, p=p)*recon_factor+ cmplx_mean
reconstructed_100 = wave.icwt(transformed[period>100000,:], dt=dt, scales=scales[period>100000], wf=wf, p=p)*recon_factor+ cmplx_mean

reconstructed_bear = (np.angle(reconstructed)*360/(2*math.pi))%360
reconstructed_bear_1000 = (np.angle(reconstructed_1000)*360/(2*math.pi))%360
reconstructed_bear_100 = (np.angle(reconstructed_100)*360/(2*math.pi))%360

plt.plot(x,theta,'0.75')
plt.plot(x,reconstructed_bear,'k--')
plt.legend(['Original Signal','Reconstructed Signal'])
plt.xlabel("Distance/km")
plt.ylabel("Bearing")
plt.title("Reconstructed bearings of Colorado River along length")
plt.show() 

# Visualise the effect of different filters
plt.plot(x,theta,'0.75')
plt.plot(x,reconstructed_bear_1000,'--')
plt.plot(x,reconstructed_bear_100,'--')
plt.legend(['Original Signal','>1000km filter','>100km filter'])
plt.xlabel("Distance/km")
plt.ylabel("Bearing")
plt.title("Filtered Colorado River bearings along length")
plt.show()

#### Easting Northings method ####

print("#######")
print("Wavelet analysis using eastings/northings method")
print("#######")

# Now we perform the wavelet transform and filtering by splitting the original
# theta signal into northings and eastings 

# calculate eastings and northings
eastings = np.sin(theta*math.pi/180)
northings = np.cos(theta*math.pi/180)

eastings_mean = np.mean(eastings)
northings_mean = np.mean(northings)

eastings = eastings-eastings_mean
northings = northings-northings_mean

# Perform wavelet transform with same scales as before
transformed_eastings = wave.cwt(x=eastings, dt=dt, scales=scales, wf=wf, p=p)
transformed_northings = wave.cwt(x=northings, dt=dt, scales=scales, wf=wf, p=p)

# Reconstruct original signal
reconstructed_eastings = np.real(wave.icwt(transformed_eastings, dt=dt, scales=scales, wf=wf, p=p))*recon_factor+ eastings_mean
reconstructed_northings = np.real(wave.icwt(transformed_northings, dt=dt, scales=scales, wf=wf, p=p))*recon_factor+ northings_mean
reconstructed_bear_en = ((180/math.pi)*np.arctan2(reconstructed_eastings,reconstructed_northings))%360

print("#######")
print("Visually comparing results from the complex and eastings/northings methods")
print("#######")

plt.scatter(reconstructed_bear_en,reconstructed_bear)
plt.xlabel("Reconstructed bearing using E/N method")
plt.ylabel("Reconstructed bearing using complex method") 
plt.title("Comparing the reconstruction by different methods") 
plt.show() 

reconstructed_east_1000 = np.real(wave.icwt(transformed_eastings[period>1000000,:], dt=dt, scales=scales[period>1000000], wf=wf, p=p))*recon_factor+ eastings_mean
reconstructed_east_100 = np.real(wave.icwt(transformed_eastings[period>100000,:], dt=dt, scales=scales[period>100000], wf=wf, p=p))*recon_factor+ eastings_mean

reconstructed_north_1000 = np.real(wave.icwt(transformed_northings[period>1000000,:], dt=dt, scales=scales[period>1000000], wf=wf, p=p))*recon_factor+ northings_mean
reconstructed_north_100 = np.real(wave.icwt(transformed_northings[period>100000,:], dt=dt, scales=scales[period>100000], wf=wf, p=p))*recon_factor+ northings_mean

reconstructed_bear_1000_en = ((180/math.pi)*np.arctan2(reconstructed_east_1000,reconstructed_north_1000))%360
reconstructed_bear_100_en = ((180/math.pi)*np.arctan2(reconstructed_east_100,reconstructed_north_100))%360

#plt.scatter(reconstructed_bear_1000_en,reconstructed_bear_1000)
#plt.xlabel("Filtered > 1000km bearing using E/N method")
#plt.ylabel("Filtered > 1000km bearing using complex method") 
#plt.title("Comparing the 1000km filter by different methods") 
#plt.show() 

#plt.scatter(reconstructed_bear_100_en,reconstructed_bear_100)
#plt.xlabel("Filtered > 100km bearing using E/N method")
#plt.ylabel("Filtered > 000km bearing using complex method") 
#plt.title("Comparing the 100km filter by different methods") 
#plt.show() 
