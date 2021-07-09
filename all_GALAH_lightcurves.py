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
#from /home/rmorris/miniconda3/envs/python37/lib/python3.7/site-packages/eleanor/mast.py import tic_from_coords
from eleanor import mast
#from mast import tic_from_coords

galah_dr3 = fits.open('GALAH_DR3_main_allstar_v2.fits')

mass_id = []
source_ids = []
target_ra = np.zeros(len(galah_dr3[1].data))
target_dec = np.zeros(len(galah_dr3[1].data))

for i, j in enumerate(galah_dr3[1].data):
    mass_id.append(j[0])
    source_ids.append(j[1])
    target_ra[i] = j['ra_dr2']
    target_dec[i] = j['dec_dr2']

#beg_ind = int(sys.argv[1])
#end_ind = int(sys.argv[2])

for j,i in enumerate(mass_id[0:2]):#[beg_ind:end_ind]):
    print(i)
    coords = SkyCoord(ra=float(target_ra[j]), dec=float(target_dec[j]), unit=(u.deg, u.deg))
    outID,outEclipLong,outEclipLat,outSec,outCam,outCcd,outColPix, outRowPix, scinfo = tess_stars2px_function_entry(j, float(target_ra[j]), float(target_dec[j]))
    TIC_info = mast.tic_from_coords((float(target_ra[j]),float(target_dec[j])))
    TIC = TIC_info[0]
    print(TIC)
    for num, sector in enumerate(outSec):
        if sector <= 37:
    ###first try to find 2-min data for this sector
            try:
                data = fits.open('/data/wallaby/rmorris/GALAH/lightkurve/{}_s{}.fits'.format(TIC,sector))
                print("Sector {} found locally with LightKurve".format(sector))
            except KeyboardInterrupt:
                raise
            except Exception as e:
                print(e)
                print('Data does not already exist @ 2-min cadence checking with Lightkurve')
                try:
                    lc = lk.search_lightcurvefile(coords, mission='TESS', sector=sector).download().PDCSAP_FLUX.remove_nans()
                    lc.to_fits(path='/data/wallaby/rmorris/GALAH/lightkurve/{}_s{}.fits'.format(TIC,sector), overwrite=True)
                    print('Sector {} downloaded with LightKurve'.format(sector))
                except KeyboardInterrupt:
                    raise
                except Exception as e:
                    print(e)
                    print('Sector {} not available with LightKurve, trying eleanor'.format(sector))
                    try:
                        data = fits.open('/data/wallaby/rmorris/GALAH/eleanor/{}_s{}.fits'.format(TIC,sector)) #have to change for katana storage
                        print("Sector {} found locally with eleanor".format(sector))
                    except Exception as e:
                        print(e)
                        print('File does not already exist, downloading postcards from eleanor instead')
                        try:
                            ###otherwise download with eleanor and save
                            star = eleanor.Source(coords=coords, sector=int(sector), post_dir='/data/wallaby/postcards/')
                            datum = eleanor.TargetData(star, height=15, width=15, bkg_size=31, do_psf=False, do_pca=False, aperture_mode='small')
                            datum.save(output_fn='{}_s{}.fits'.format(TIC,sector),directory='/data/wallaby/rmorris/GALAH/eleanor/')
                            print("Sector {} downloaded with eleanor".format(sector))
                        except KeyboardInterrupt:
                            raise
                        except Exception as e:
                            print(e)


