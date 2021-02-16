import csv
import sys
import eleanor
import numpy as np
import lightkurve as lk
from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.timeseries import BoxLeastSquares
from transitleastsquares import transitleastsquares
from tess_stars2px import tess_stars2px_function_entry


def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

galah_dr3 = fits.open('GALAH_DR3_main_allspec_v1.fits')

mass_id = []
source_ids = []
target_ra = np.zeros(len(galah_dr3[1].data))
target_dec = np.zeros(len(galah_dr3[1].data))

for i, j in enumerate(galah_dr3[1].data):
	mass_id.append(j[0])
	source_ids.append(j[1])
	target_ra[i] = j[445]
	target_dec[i] = j[447]

beg_ind = int(sys.argv[1])
end_ind = int(sys.argv[2])

best_fit_period = np.zeros(len(mass_id[beg_ind:end_ind]))
best_fit_uncert = np.zeros(len(mass_id[beg_ind:end_ind]))
#data_source = np.zeros(len(mass_id[beg_ind:end_ind]))
for j,i in enumerate(mass_id[beg_ind:end_ind]):
	print(i)
	coords = SkyCoord(ra=float(target_ra[j]), dec=float(target_dec[j]), unit=(u.deg, u.deg))
	outID,outEclipLong,outEclipLat,outSec,outCam,outCcd,outColPix, outRowPix, scinfo = tess_stars2px_function_entry(j, float(target_ra[j]), float(target_dec[j]))
	sector_source = []			   #array of where the data came from length of the amount of sectors
	#lightcurves = lk.LightCurveCollection()
	for num, sector in enumerate(outSec):
		if sector <= 31:
	###first try to find 2-min data for this sector
			try:
				lc = lk.search_lightcurvefile(coords, mission='TESS', sector=sector).download().PDCSAP_FLUX.remove_nans()
				if sector == outSec[0]:
					lightcurves = lk.LightCurveCollection(lc)
				else:
					lightcurves.append(lc)
				sector_source.append('2-min')
				print('Sector {} downloaded with LightKurve'.format(sector))
			except KeyboardInterrupt:
				raise
			except Exception as e:
				print(e)
				print('Sector {} not available with LightKurve, trying eleanor'.format(sector))
				try:
					###then try to see if the fits file already exists
					time = []
					flux = []
					background = []
					data = fits.open('/data/wallaby/rmorris/GALAH/all_target_lc/{}_s{}_eleanor_data.fits'.format(i,sector))
					for x in data[1].data:
						if x[7] == 0:
							time.append(x[0])
							flux.append(x[5])
							background.append(x[12])

					med_corr_flux = np.median(flux)
					for fluxval in flux:
						fluxval = fluxval/med_corr_flux

					time_bkg = np.zeros(len(background))
					flux_bkg = np.zeros(len(background))
					for num1, back in enumerate(background):
						bkg_med = np.median(background)
						if back < bkg_med*3 and back > bkg_med/3:
							time_bkg[num1] = time[num1]
							flux_bkg[num1] = flux[num1]

					lc = lk.LightCurve(time=time_bkg, flux=flux_bkg)
					#print('test line')
					if sector == outSec[0]:
						lightcurves = lk.LightCurveCollection(lc)
					else:
						lightcurves.append(lc)
					#print('test line 2')
					sector_source.append('30-min')
					print("Sector {} found locally with eleanor".format(sector))
				except Exception as e:
					print(e)
					print('File does not already exist, downloading postcards instead')
					try:
						###otherwise download with eleanor and save
						star = eleanor.Source(coords=coords, sector=int(sector), post_dir='/data/wallaby/postcards/')

						datum = eleanor.TargetData(star, height=15, width=15, bkg_size=31, do_psf=False, do_pca=False, aperture_mode='small')
						datum.save(output_fn='{}_s{}_eleanor_data.fits'.format(i,sector),directory='/data/wallaby/rmorris/GALAH/all_target_lc/')

						time = []
						flux = []
						background = []

						q = datum.quality == 0
						time.append(datum.time[q])
						flux.append(datum.corr_flux[q]/np.median(datum.corr_flux[q]))
						background.append(datum.flux_bkg[q])

						time_bkg = np.zeros(len(background))
						flux_bkg = np.zeros(len(background))
						for num2, back in enumerate(background):
							bkg_med = np.median(background)
							if back < bkg_med*3 and back > bkg_med/3:
								time_bkg[num2] = time[num2]
								flux_bkg[num2] = flux[num2]

						lc = lk.LightCurve(time=time_bkg, flux=flux_bkg)
						if sector == outSec[0]:
							lightcurves = lk.LightCurveCollection(lc)
						else:
							lightcurves.append(lc) 
						sector_source.append('30-min')
						print("Sector {} downloaded with eleanor".format(sector))
					except KeyboardInterrupt:
						raise
					except Exception as e:
						print(e)

	#data_source[j] = sector_source
	#collection = lk.LightCurveCollection(lightcurves)
	lc = lightcurves.stitch()

	try:
		###create lightcurve object
		lc = lc.remove_outliers(sigma_lower=7, sigma_upper=3).flatten(window_length=51)
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
		best_fit_period[j] = results.period
		best_fit_uncert[j] = results.period_uncertainty

		fold_lc = lc.fold(period=results.period)#, t0=periodogram.transit_time_at_max_power)
		#fold_lc_real = lc.fold(period=tic_period[j])
		#plt.figure(num=str(i) + "_3")
		#plt.scatter(fold_lc.time, fold_lc.flux, c='k', marker='.')
		#plt.scatter(fold_lc_real.time, fold_lc_real.flux, c='r', marker='.')
		#plt.title(i)
		#plt.xlabel('Phase')
		#plt.ylabel('Normalized Flux')
		#plt.legend(['TLS Folded', 'ExoFOP Folded'])


		print("Lightcurve successfully folded")


		#print(toi_period[j], best_fit_period[j])

		###can do this separately after all the data is processed
		#if results.period < tic_period[j]*1.5 and results.period > tic_period[j]*0.5:
		#	well_fit_target.append(i)
		#if results.period > tic_period[j]*1.5 or results.period < tic_period[j]*0.5:
		#	poor_fit_target.append(i)
	except: 
		print('Target {} failed'.format(i))
		#failed_tic.append(i)
		best_fit_period[j] = 0.0
		best_fit_uncert[j] = 0.0
	try:
		if end_ind-beg_ind > 10000:
			if j % 10000 == 0:
				with open('/home/rmorris/documents/lc_progress/eleanor_tls_periods_0-{}.txt'.format(j), 'w') as file:
					for number_ind, values_ind in enumerate(mass_id[0:j]):
						file.write(str(mass_id[number_ind])+'\t'+str(best_fit_period[number_ind])+'\t'+str(best_fit_uncert[number_ind])+'\n')
			else:
				print('Not writing yet :)')
		else:
			print('Less than 10k targets, will only print at the end')
	except Exception as e:
		print(e)
		print('\nWriting to file failed \n')

	finally:
		plt.close('all')



#print(best_fit_period[0:(len(best_fit_period))])
#plt.show()

#print(np.shape(toi_period), np.shape(best_fit_period))
#period_arrays = np.stack((tic_period, tic_period_error, best_fit_period, best_fit_uncert), axis=1)
#print(period_arrays)

with open('/home/rmorris/documents/eleanor_tls_periods_{}-{}.txt'.format(beg_ind, end_ind), 'w') as file:
	for number, values in enumerate(mass_id[beg_ind:end_ind]):
		file.write(str(mass_id[number])+'\t'+str(best_fit_period[number])+'\t'+str(best_fit_uncert[number])+'\n')

#print(well_fit_target)
#print(poor_fit_target)
#print(len(tic_ids), len(poor_fit_target), len(well_fit_target))


#plt.scatter(tic_period,best_fit_period)
#plt.xlabel('ExoFOP Period')
#plt.ylabel('TLS Period')
#plt.title('Published vs Calculated Period')
#plt.savefig('/home/rmorris/documents/eleanor_tls_per_plot.png')


