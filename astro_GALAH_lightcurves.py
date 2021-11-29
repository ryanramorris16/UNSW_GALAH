import os
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
import ticgen
#from mast import tic_from_coords

galah_dr3 = fits.open('GALAH_DR3_main_allstar_v2.fits')

mass_id = []
source_ids = []
target_ra = np.zeros(len(galah_dr3[1].data))
target_dec = np.zeros(len(galah_dr3[1].data))
j_mags = np.zeros(len(galah_dr3[1].data))
ks_mags = np.zeros(len(galah_dr3[1].data))
h_mags = np.zeros(len(galah_dr3[1].data))

done_id = []
done_sec = []
with open('/home/z5318114/targets_ran.txt','r') as file:
    reader = csv.reader(file, delimiter='\t')
    for row in reader:
        try:
            done_id.append(row[0])
            done_sec.append(row[1])
        except:
            print(row)

for i, j in enumerate(galah_dr3[1].data):
    mass_id.append(j[0])
    source_ids.append(j[1])
    target_ra[i] = j['ra_dr2']
    target_dec[i] = j['dec_dr2']
    j_mags[i] = j['j_m']
    ks_mags[i] = j['ks_m']
    h_mags[i] = j['h_m']

arr_ind = int(sys.argv[1])
round = arr_ind // 1200
mod = arr_ind % 1200
beg_ind = (mod-1)*600 + (50*round)
end_ind = beg_ind+50
#print(beg_ind, end_ind)
#end_ind = int(sys.argv[2])

for j,i in enumerate(mass_id[beg_ind:end_ind], start=beg_ind):
    print(j)
    coords = SkyCoord(ra=float(target_ra[j]), dec=float(target_dec[j]), unit=(u.deg, u.deg))
    outID,outEclipLong,outEclipLat,outSec,outCam,outCcd,outColPix, outRowPix, scinfo = tess_stars2px_function_entry(j, float(target_ra[j]), float(target_dec[j]))
    TIC_info = mast.tic_from_coords((float(target_ra[j]),float(target_dec[j])))
    TIC = TIC_info[0]
    Tmag_eleanor = TIC_info[1][0]
    print(TIC, Tmag_eleanor)
    star_ticgen = ticgen.Star(Jmag=j_mags[j], Ksmag=ks_mags[j], Hmag=h_mags[j])
    Tmag_ticgen = star_ticgen.Tmag
    if Tmag_eleanor > 13.0:
        eleanor_ap = 'small'
    else:
        eleanor_ap = 'normal'
    for num, sector in enumerate(outSec):
        skip = 0
        if source_ids[j] in done_id:
            done_indices = [qj for qj, qi in enumerate(done_id) if qi == source_ids[j]]   #done_id.index(source_ids[j])
            for asd in done_indices:
                if sector == done_sec[asd]:
                    skip = 1
                    break
                else:
                    skip = 0
        if skip == 0:     
            if sector < 39.5:
    ###first try to find 2-min data for this sector
                try:
                    print(j, sector)
                    data = fits.open('/srv/scratch/astro/z5318114/GALAH/lightkurve/{}_s{}_120s.fits'.format(TIC,sector))
                    source='lightkurve'
                    ending = '120s'
                    #print("Sector {} found locally with LightKurve".format(sector))
                except KeyboardInterrupt:
                    raise
                except Exception as e:
                    #print(e)
                    #print('Data does not already exist @ 2-min cadence with Lightkurve')
                    try:
                        data = fits.open('/srv/scratch/astro/z5318114/GALAH/lightkurve/{}_s{}_20s.fits'.format(TIC,sector))
                        source = 'lightkurve'
                        ending = '20s'
                    except Exception as e:
                        #print('Data does not already exist @ 20-s cadence with Lightkurve')
                        try:
                            data = fits.open('/srv/scratch/astro/z5318114/GALAH/eleanor/{}_s{}_ffi.fits'.format(TIC,sector)) #have to change for katana storage
                            source='eleanor'
                            ending = 'ffi'
                            #print("Sector {} found locally with eleanor".format(sector))
                        except KeyboardInterrupt:
                            raise
                        except Exception as e:
                            #print(e)
                            #print('Data does not already exist @ 30/10-min cadence with eleanor')
                            try:
                                lc = lk.search_lightcurve(coords, mission='TESS', sector=sector, exptime='short',author='SPOC')#.PDCSAP_FLUX.remove_nans()
                                twomin_file = lc.table[0]['productFilename']
                                lc.download(download_dir='/srv/scratch/astro/z5318114/GALAH/lightkurve/')
                                os.rename('/srv/scratch/astro/z5318114/GALAH/lightkurve/mastDownload/TESS/{}/{}'.format(twomin_file[0:-8],twomin_file),'/srv/scratch/astro/z5318114/GALAH/lightkurve/{}_s{}_120s.fits'.format(TIC,sector))
                                source = 'lightkurve'
                                #print('Sector {} downloaded with LightKurve'.format(sector))
                            except Exception as e:
                                #print(e)
                                #print('Data not available with LightKurve, downloading with eleanor instead')
                                try:
                                    lc = lk.search_lightcurve(coords, mission='TESS', sector=sector, exptime='fast',author='SPOC')
                                    twentys_file = lc.table[0]['productFilename']
                                    lc.download(download_dir='/srv/scratch/astro/z5318114/GALAH/lightkurve/')
                                    os.rename('/srv/scratch/astro/z5318114/GALAH/lightkurve/mastDownload/TESS/{}/{}'.format(twentys_file[0:-8],twentys_file), '/srv/scratch/astro/z5318114/GALAH/lightkurve/{}_s{}_20s.fits'.format(TIC,sector))
                                    source = 'lightkurve'	
                                except Exception as e:
                                    #print(e)
                                    #print('This target has 20s data, very cool man, very cool.')
                                    try:
                                        #print(j, sector)
                                        ###otherwise download with eleanor and save
                                        star = eleanor.Source(coords=coords, sector=int(sector), post_dir='/srv/scratch/astro/z5318114/postcards/')
                                        datum = eleanor.TargetData(star, height=15, width=15, bkg_size=31, do_psf=False, do_pca=True, aperture_mode=eleanor_ap, save_postcard=False)
                                        datum.save(output_fn='{}_s{}_ffi.fits'.format(TIC,sector),directory='/srv/scratch/astro/z5318114/GALAH/eleanor/')
                                        source = 'eleanor'
                                        #print("Sector {} downloaded with eleanor".format(sector))
                                    except KeyboardInterrupt:
                                        raise
                                    except Exception as e:
                                        print("Nah, this one isn't happening:", e)

                try:                
                    if os.path.isfile('/srv/scratch/astro/z5318114/GALAH/{}/{}_s{}_{}.fits'.format(source,TIC,sector,ending)):    
                        #os.system('rclone copy /srv/scratch/astro/z5318114/GALAH/{}/{}_s{}.fits CloudStor:/{}/'.format(source,TIC,sector,source))
                        totalfile = open('/home/z5318114/targets_ran.txt',"a")
                        totalfile.write(str(source_ids[j])+'\t'+str(sector)+'\t'+str(TIC)+'\t'+str(Tmag_eleanor)+'\t'+str(Tmag_ticgen)+'\t'+str(source)+'\n')
                        totalfile.close()
                except:
                    print('Something went wrong, file not run correctly')
    #else:
        ### if the target is already in targets_run.txt, don't run again
        #print('Already in list')    

#print(end_ind)
print('Done')

