import eleanor
import numpy as np
import numpy.ma as ma
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import lightkurve as lk
from astropy.io import fits
import csv
from astropy.timeseries import BoxLeastSquares
from astroquery.mast import Observations
from transitleastsquares import transitleastsquares
import math


def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

tic_ids = []
tic_period = []
tic_period_error = []
tic_ra = []
tic_dec =[]
with open('/home/rmorris/documents/exofop_tess_tois_faint_no_FP.txt', 'r') as file:
	reader = csv.reader(file,delimiter='\t')
	for row in reader:
		tic_ids.append(int(row[0].replace("'","")))
		tic_period.append(float(row[25].replace("'","")))
		if isfloat(row[26]) == True:
			tic_period_error.append(float(row[26].replace("'","")))
		elif isfloat(row[26]) == False:
			tic_period_error.append(float(0.0))
		tic_ra.append(float(row[17]))
		tic_dec.append(float(row[18]))

best_fit_period = []
failed_tic = []
best_fit_uncert = []
for j,i in enumerate(tic_ids):
	try:
		print(i)
		
		target_name = "TIC " + str(i)
		#print(target_name)
		data_table = Observations.query_criteria(provenance_name="SPOC",objectname=target_name,obs_collection="TESS",dataproduct_type="timeseries",radius=0)
		#print(data_table)
		#print(Observations.get_metadata("products"))
		product_list = Observations.get_product_list(data_table)
		#print(product_list)
		download_data = Observations.download_products(product_list,extension="lc.fits",download_dir="/data/wallaby/rmorris/SPOC/{}".format(tic_ids[j]))
		#print(download_data[1][0])

		print("Data downloaded correctly")

	except:
		print("Connection to MAST Not Established")

	try:
		time = []
		norm_flux = []
		###might not need background because already is PDCSAP flux?
		ticker = 0
		while ticker < len(download_data):
			file = fits.open(str(download_data[ticker][0]))
			ind = 0
			flux = []
			while ind < len(file[1].data):
				time.append(file[1].data[ind][0])	###data[X][0] is time
				flux.append(file[1].data[ind][7])	###[X][7] is PDCSAP flux
				ind += 1
			med_flux = np.nanmedian(flux)
			for flux_val in flux:
				flux_val = flux_val/med_flux
				norm_flux.append(flux_val)	
			ticker += 1

		print('Flux Normalized')

		time_array = []
		flux_array = []
		for ind, value in enumerate(norm_flux):
			if value == value:
				time_array.append(time[ind])
				flux_array.append(norm_flux[ind])


		#time_bkg = []
		#flux_bkg = []
		#for num, back in enumerate(background_array):
		#	bkg_med = np.median(background_array)
		#	if back < bkg_med*3:
		#		time_bkg.append(time_array[num])
		#		flux_bkg.append(flux_array[num])


		#print("Bad background points subtracted")
		#print(len(time_bkg))

		###create lightcurve object
		lc = lk.LightCurve(time=time_array, flux=flux_array).remove_outliers(sigma_lower=7, sigma_upper=3).flatten(window_length=51)
		print('LightKurve object created')

		#print(plt)
		#print(plt.scatter)
		#plt.figure(num=str(i) + "_1")
		#plt.scatter(time_array,flux_array, c='k', marker='.')
		#print('scatter worked')
		#plt.title(i)
		#plt.xlabel('Time')
		#plt.ylabel('Normalized Flux')
		#plt.show()

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

		best_fit_period.append(results.period)
		best_fit_uncert.append(results.period_uncertainty)

		fold_lc = lc.fold(period=results.period)#, t0=periodogram.transit_time_at_max_power)
		fold_lc_real = lc.fold(period=tic_period[j])
		#plt.figure(num=str(i) + "_3")
		#plt.scatter(fold_lc.time, fold_lc.flux, c='k', marker='.')
		#plt.scatter(fold_lc_real.time, fold_lc_real.flux, c='r', marker='.')
		#plt.title(i)
		#plt.xlabel('Phase')
		#plt.ylabel('Normalized Flux')
		#plt.legend(['TLS Folded', 'ExoFOP Folded'])


		print("Lightcurve successfully folded")

	except: 
		print('Target {} failed'.format(i))
		failed_tic.append(i)
		best_fit_period.append(0.0)
		best_fit_uncert.append(0.0)
	finally:
		plt.close('all')

#print(np.shape(toi_period), np.shape(best_fit_period))
#period_arrays = np.stack((tic_period, tic_period_error, best_fit_period, best_fit_uncert), axis=1)
#print(period_arrays)

with open('/home/rmorris/documents/spoc_periods.txt', 'w') as file:
	for number, values in enumerate(tic_period):
		file.write(str(tic_period[number])+'\t'+str(tic_period_error[number])+'\t'+str(best_fit_period[number])+'\t'+str(best_fit_uncert[number])+'\n')

#print(well_fit_target)
#print(poor_fit_target)
#print(len(tic_ids), len(poor_fit_target), len(well_fit_target))


plt.scatter(tic_period,best_fit_period)
plt.xlabel('ExoFOP Period')
plt.ylabel('TLS Period')
plt.title('Published vs Calculated Period')
plt.savefig('/home/rmorris/documents/spoc_per_comparison.png')


