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

stream = fits.open('s5_pdr1.fits')

gaia_ids = []
stream_ra = []
stream_dec = []
priority = []
i_mag = []

for i, j in enumerate(stream[1].data):
    gaia_ids.append(j['name'])
    stream_ra.append(j['ra'])
    stream_dec.append(j['dec'])
    priority.append(j['priority'])
    i_mag.append(j['i_skm'])
    
best_fit_per = np.zeros(len(gaia_ids))
best_fit_unc = np.zeros(len(gaia_ids))

for j, i in enumerate(gaia_ids):
    if int(priority[j]) > 2 and i_mag[j] < 16:    # only priority 3-9 (halo and stream stars) and bright targets
        coords = SkyCoord(ra=float(stream_ra[j]), dec=float(stream_dec[j]), unit=(u.deg, u.deg))
        outID,outEclipLong,outEclipLat,outSec,outCam,outCcd,outColPix, outRowPix, scinfo = tess_stars2px_function_entry(j, float(stream_ra[j]), float(stream_dec[j]))
        sector_source = []

        for num, sector in enumerate(outSec):
            if sector <= 35:
                try:
                    lc = lk.search_lightcurve(coords,mission='TESS',sector=sector).download().PDCSAP_FLUX.remove_nans()
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
                        time = []
                        flux = []
                        background = []
                        data = fits.open('/data/wallaby/rmorris/GALAH/all_target_lc/{}_s{}_eleanor_data.fits'.format(i,sector))
                        for x in data[1].data:
                            if x[7] == 0:
                                time.append(x[0])
                                flux.append(x[5])
                                background.append(x[12])
                        
                        #plt.title('raw')
                        #plt.plot(time, flux)
                        #plt.show()
                                
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
                                
                        time_ind = [index for index, value in enumerate(time_bkg) if value != 0]
                        time_bkg = time_bkg[time_ind]
                        flux_bkg = flux_bkg[time_ind]
                        
                        #plt.title('background corrected')
                        #plt.plot(time_bkg, flux_bkg)
                        #plt.show()

                        lc = lk.LightCurve(time=time_bkg, flux=flux_bkg)
                        #plt.title('in the lc object')
                        #plt.plot(lc.time.value, lc.flux.value)
                        #plt.show()
                        if sector == outSec[0]:
                            lightcurves = lc
                            #plt.title('as lightcurves')
                            #plt.plot(lightcurves.time.value, lightcurves.flux.value)
                            #plt.show()
                        else:
                            lightcurves.append(lc)
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
                            
                            #plt.title('downloading - raw')
                            #plt.plot(time, flux)
                            #plt.show()

                            time_bkg = np.zeros(len(background))
                            flux_bkg = np.zeros(len(background))
                            for num2, back in enumerate(background):
                                bkg_med = np.median(background)
                                if back < bkg_med*3 and back > bkg_med/3:
                                    time_bkg[num2] = time[num2]
                                    flux_bkg[num2] = flux[num2]
                                    
                            time_ind = [index for index, value in enumerate(time_bkg) if value != 0]
                            time_bkg = time_bkg[time_ind]
                            flux_bkg = flux_bkg[time_ind]

                            lc = lk.LightCurve(time=time_bkg, flux=flux_bkg)
                            if sector == outSec[0]:
                                lightcurves = lc
                            else:
                                lightcurves.append(lc) 
                            sector_source.append('30-min')
                            print("Sector {} downloaded with eleanor".format(sector))
                        except KeyboardInterrupt:
                            raise
                        except Exception as e:
                            print(e)

        #lc = lightcurves.stitch()      # shouldn't need to stitch anymore

        try:
            ###create lightcurve object
            lc = lightcurves.remove_outliers(sigma_lower=7, sigma_upper=3).flatten(window_length=51)
            #plt.figure(num=str(i) + "_1")
            #plt.scatter(lc.time.value,lc.flux.value, c='k', marker='.')
            #plt.title(i)
            #plt.xlabel('Time')
            #plt.ylabel('Normalized Flux')
            #plt.show()
            print('LightKurve object created')

            ###TLS through transitleastsquares package
            periodogram = transitleastsquares(lc.time.value, lc.flux.value)
            results = periodogram.power(period_max=75.0, oversampling_factor=5, duration_grid_step=1.2)#, use_threads=1)

            #plt.figure(num=str(i) + "_2")
            #plt.plot(results.periods, results.power, c='k', lw='0.5')
            #plt.title(i)
            #plt.xlabel('Frequency')
            #plt.ylabel('TLS Power')
            #plt.show()
            print("Period generated with TLS")

            best_fit_per[j] = results.period
            best_fit_unc[j] = results.period_uncertainty
            #print(best_fit_per[j], best_fit_unc[j])

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
        except: 
            print('Target {} failed'.format(i))
            best_fit_per[j] = 0.0
            best_fit_unc[j] = 0.0
        finally:
            plt.close('all') 
        '''
        try:
            ### only writing out to file every 10000 targets
            if end_ind-beg_ind > 10000:
                if j % 10000 == 0:
                    with open('/home/rmorris/documents/lc_progress/steam_periods_0-{}.txt'.format(j), 'w') as file:
                        for number_ind, values_ind in enumerate(gaia_id[0:j]):
                            file.write(str(gaia_id[number_ind])+'\t'+str(best_fit_per[number_ind])+'\t'+str(best_fit_unc[number_ind])+'\n')
                else:
                    print('Not writing yet :)')
            else:
                print('Less than 10k targets, will only print at the end')
        except Exception as e:
            print(e)
            print('\nWriting to file failed \n')
        '''
        
    else:
        print('Target is neither a halo nor a stream star')
    
with open('/home/rmorris/documents/stream_periods.txt', 'w') as file:
    for number, values in enumerate(gaia_ids):
        if int(priority[number]) > 2 and i_mag[number] < 16:
            file.write(str(gaia_ids[number])+'\t'+str(priority[number])+'\t'+str(best_fit_per[number])+'\t'+str(best_fit_unc[number])+'\n')
            print('writing...')
            
print('done')
                    
                    
                    
                  
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                
                
                
                    