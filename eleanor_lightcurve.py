import eleanor
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import lightkurve as lk
from astropy.io import fits

hdul = fits.open('GALAH_DR3_main_200331.fits')
print(hdul[1].data[0][166:170],hdul[1].data[0][2])

star = eleanor.multi_sectors(gaia=3231972373528734848, sectors='all')

#print('Found TIC {0} (Gaia {1}), with TESS magnitude {2}, RA {3}, and Dec {4}'
#     .format(star.tic, star.gaia, star.tess_mag, star.coords[0], star.coords[1]))

#data = eleanor.TargetData(star, height=15, width=15, bkg_size=31, do_psf=True, do_pca=True, regressors='corner')

plt.figure(figsize=(15,5))

data = []
plot_fmt = ['k.', 'r.','k.', 'r.']

for s in star:
    datum = eleanor.TargetData(s, height=15, width=15, bkg_size=31, do_psf=False, do_pca=False)
    data.append(datum)

for sector, datum in enumerate(data):
    q = datum.quality == 0
    plt.plot(datum.time[q], datum.corr_flux[q]/np.median(datum.corr_flux[q]), plot_fmt[sector])

plt.ylabel('Normalized Flux', fontsize=24)
plt.xlabel('Time', fontsize=24)

plt.show()