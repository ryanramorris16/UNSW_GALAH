import eleanor
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import lightkurve as lk
from astropy.io import fits

hdul = fits.open('GALAH_DR3_main_200331.fits')
#rint(hdul[1].data[0][166:170],hdul[1].data[0][2])

star_ids = []
for i in hdul[1].data:
	star_ids.append(i[2])

star_ids_test = star_ids[9:10]

for i in star_ids_test:
	star = eleanor.multi_sectors(gaia=i, sectors='all')

	#print('Found TIC {0} (Gaia {1}), with TESS magnitude {2}, RA {3}, and Dec {4}'
	#     .format(star.tic, star.gaia, star.tess_mag, star.coords[0], star.coords[1]))

	#data = eleanor.TargetData(star, height=15, width=15, bkg_size=31, do_psf=True, do_pca=True, regressors='corner')

	plt.figure(figsize=(15,5))

	data = []
	plot_fmt = ['k.', 'r.','k.', 'r.']

	for s in star:
	    datum = eleanor.TargetData(s, height=15, width=15, bkg_size=31, do_psf=False, do_pca=False)
	    data.append(datum)

	time = []
	flux = []
	for sector, datum in enumerate(data):
	    q = datum.quality == 0
	    #plt.plot(datum.time[q], datum.corr_flux[q]/np.median(datum.corr_flux[q]), plot_fmt[sector])
	    time.append(datum.time[q])
	    flux.append(datum.corr_flux[q]/np.median(datum.corr_flux[q]))

	#plt.ylabel('Normalized Flux', fontsize=24)
	#plt.xlabel('Time', fontsize=24)

	#plt.show()

	#print(time, flux)

	time_array = np.concatenate((time[0], time[1], time[2], time[3], time[4], time[5], time[6]))
	flux_array = np.concatenate((flux[0], flux[1], flux[2], flux[3], flux[4], flux[5], flux[6]))

	###create lightcurve object
	lc = lk.LightCurve(time = time_array, flux = flux_array)
	lc.plot()
	plt.show()

	periodogram = lc.to_periodogram(method='bls', period=np.arange(0.3, 50.0, 0.01))
	periodogram.plot()
	plt.show()
	#print(periodogram)

