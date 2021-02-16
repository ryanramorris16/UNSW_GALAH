import eleanor
import numpy as np
from astropy import units as u
#import matplotlib		###was supposed to fix some weird failing in code, didn't
#matplotlib.use('TKagg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import lightkurve as lk
from astropy.io import fits
import csv
from astropy.timeseries import BoxLeastSquares
from transitleastsquares import transitleastsquares

#hdul = fits.open('GALAH_DR3_main_200331.fits')
#print(hdul[1].data[0][166:170],hdul[1].data[0][2])

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

tic_ids = 296780789
tic_period = 13.577007
tic_period_error = 0.000232
'''
with open('/home/ryan/Downloads/exofop_tess_tois_faint_no_FP.txt', 'r') as file:
	reader = csv.reader(file,delimiter='\t')
	#next(reader)
	#next(reader)
	#next(reader)
	#next(reader)
	#next(reader)
	#next(reader)
	for row in reader:
		tic_ids.append(int(row[0].replace("'","")))
		tic_period.append(float(row[25].replace("'","")))
		if isfloat(row[26]) == True:
			tic_period_error.append(float(row[26].replace("'","")))
		elif isfloat(row[26]) == False:
			tic_period_error.append(float(0.0))
'''
#print(toi_ids, toi_period, len(toi_ids))
'''
star_ids = []
for i in hdul[1].data:
	star_ids.append(i[2])

star_ids_test = [234284556, 234523599]#, 267263253, 410153553,140068425,261136679,179317684,29857954,425997655,307210830,279741379,441462736,410214986,55652896,200723869] #star_ids[9:10]
real_period = [1.106,3.796]#,4.127,0.463,2.281,6.268,4.231,9.477,17.667,2.253,7.790,14.277,8.138,17.089,18.371,]
'''
best_fit_period = []
#well_fit_target = []
#poor_fit_target = []
failed_tic = []
best_fit_uncert = []

try:
	print(tic_ids)
	star = eleanor.multi_sectors(tic=tic_ids, sectors='all')

	#print('Found TIC {0} (Gaia {1}), with TESS magnitude {2}, RA {3}, and Dec {4}'
	#     .format(star.tic, star.gaia, star.tess_mag, star.coords[0], star.coords[1]))

	#data = eleanor.TargetData(star, height=15, width=15, bkg_size=31, do_psf=True, do_pca=True, regressors='corner')

	#plt.figure(figsize=(15,5))

	data = []
	#plot_fmt = ['k.', 'r.','k.', 'r.']

	for s in star:
	    datum = eleanor.TargetData(s, do_psf=False, do_pca=False, aperture_mode='small')
	    data.append(datum)

	plt.figure(num=150)
	plt.scatter(data[0].time, data[0].quality, c='k', marker='.')
	plt.show()

	#plt.imshow(data[0].aperture)
	#plt.show()

	#print(data[0].tpf[400])
	#print(data[0].cen_x, data[0].cen_y)
	#print(data[0].tpf_flux_bkg)


	time = []
	flux = []
	raw_flux = []
	corr_flux = []
	background = []
	for sector, datum in enumerate(data):
	    q = datum.quality == 0
	    time.append(datum.time[q])
	    raw_flux.append(datum.raw_flux[q])
	    corr_flux.append(datum.all_corr_flux[0,q])
	    flux.append(datum.corr_flux[q])#/np.median(datum.corr_flux[q]))
	    background.append(datum.flux_bkg[q])

	print(eleanor.targetdata.get_flattened_sigma(data[0].all_corr_flux[0,q]))
	print(data[0].bkg_type)
	#print(eleanor.__version__)
	#print(star[0].pointing)
	plt.figure(num=15)
	plt.scatter(time, corr_flux, c='k', marker='.')
	plt.show()

	print("Data downloaded correctly")


	#time = np.concatenate((time))
	#flux = np.concatenate((flux))

	#plt.ylabel('Normalized Flux', fontsize=24)
	#plt.xlabel('Time', fontsize=24)

	#plt.show()

	#print(time, flux)
	time_array = np.concatenate((time[0:(len(time))]))
	flux_array = np.concatenate((flux[0:(len(flux))]))
	background_array = np.concatenate((background[0:len(background)]))

	time_bkg = []
	flux_bkg = []
	for num, back in enumerate(background_array):
		bkg_med = np.median(background_array)
		if back < bkg_med*3:
			time_bkg.append(time_array[num])
			flux_bkg.append(flux_array[num])


	print("Bad background points subtracted")

	print(len(time_bkg))

	###create lightcurve object
	lc = lk.LightCurve(time = time_bkg, flux = flux_bkg).remove_outliers(sigma_lower=100, sigma_upper=3).flatten(51)
	#plt.figure(num=str(i) + "_1")
	#plt.scatter(lc.time,lc.flux, c='k', marker='.')
	#plt.title(i)
	#lt.xlabel('Time')
	#plt.ylabel('Normalized Flux')
	#plt.show()

	print('LightKurve object created')

	###Different methods of periodogram, BLS through lightkurve
	#periodogram = lc.to_periodogram(method='bls', period=np.arange(0.3, 75.0, 0.001))

	###TLS through transitleastsquares package
	periodogram = transitleastsquares(lc.time, lc.flux)
	results = periodogram.power(period_max=75.0, oversampling_factor=1, duration_grid_step=1.2)#, use_threads=1)

	#plt.figure(num=str(i) + "_2")
	#plt.plot(results.periods, results.power, c='k', lw='0.5')
	#plt.title(i)
	#plt.xlabel('Frequency')
	#plt.ylabel('TLS Power')
	#plt.show()


	print("Period generated with TLS")


	###this is specific for dimension-ed quantites (from lightkurve/astropy BLS)
	#dimen_per = periodogram.period_at_max_power
	#dimless_per = float(dimen_per / (1. * u.d))
	#best_fit_period.append(dimless_per)
	best_fit_period = results.period
	best_fit_uncert = results.period_uncertainty

	print('Best fit found')
	print(best_fit_period, best_fit_uncert, tic_period, tic_period_error)

	fold_lc = lc.fold(period=results.period)#, t0=periodogram.transit_time_at_max_power)
	fold_lc_real = lc.fold(period=tic_period).bin(binsize=1)
	plt.figure(num=str(tic_ids) + "_3")
	plt.scatter(fold_lc.time, fold_lc.flux, c='k', marker='.')
	plt.scatter(fold_lc_real.time, fold_lc_real.flux, c='r', marker='.')
	plt.title(tic_ids)
	plt.xlabel('Phase')
	plt.ylabel('Normalized Flux')
	plt.legend(['TLS Folded', 'ExoFOP Folded'])
	plt.show()


	print("Lightcurve successfully folded")



	#print(toi_period[j], best_fit_period[j])

	###can do this separately after all the data is processed
	#if results.period < tic_period[j]*1.5 and results.period > tic_period[j]*0.5:
	#	well_fit_target.append(i)
	#if results.period > tic_period[j]*1.5 or results.period < tic_period[j]*0.5:
	#	poor_fit_target.append(i)
except: 
	print('Target {} failed'.format(tic_ids))
	#failed_tic.append(i)
	best_fit_period = 0.0
	best_fit_uncert = 0.0
finally:
	plt.close('all')



#print(best_fit_period[0:(len(best_fit_period))])
#plt.show()

#print(np.shape(toi_period), np.shape(best_fit_period))
#period_arrays = np.stack((tic_period, tic_period_error, best_fit_period, best_fit_uncert), axis=1)
#print(period_arrays)
'''
with open('/home/ryan/Documents/UNSW/UNSW_GALAH/tls_bkg_periods_small_sigma.txt', 'w') as file:
	for number, values in enumerate(tic_period):
		file.write(str(tic_period[number])+'\t'+str(tic_period_error[number])+'\t'+str(best_fit_period[number])+'\t'+str(best_fit_uncert[number])+'\n')
'''
#print(well_fit_target)
#print(poor_fit_target)
#print(len(tic_ids), len(poor_fit_target), len(well_fit_target))

'''
plt.scatter(tic_period,best_fit_period)
plt.xlabel('ExoFOP Period')
plt.ylabel('TLS Period')
plt.title('Published vs Calculated Period')
plt.show()

'''
