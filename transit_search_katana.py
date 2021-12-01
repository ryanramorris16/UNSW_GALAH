import os
import csv
import glob
import math
import numpy as np
import lightkurve as lk
from eleanor import mast
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from transitleastsquares import transit_mask
from target_data_functions import corrected_flux
from transitleastsquares import transitleastsquares
#from transitleastsquares.stats import calculate_stretch
#from transitleastsquares.stats import all_transit_times
#from transitleastsquares.transit import fractional_transit
#import transitleastsquares.tls_constants as tls_constants
from tess_stars2px import tess_stars2px_function_entry



class transit:
    """
    parameters:
        ra: right acension of the target star
        dec: declination of the target star
        directory: directory where folders of light curves are stored, separated by source (eleanor vs lightkurve) (/data/wallaby/rmorris/GALAH/)
        rerun: if True, run every target given, if False, only run a target if it has not been run before
        
    """
    
    def __init__(self, ra, dec, directory, rerun=True):
        TIC_info = mast.tic_from_coords((float(ra), float(dec)))
        self.ra = float(ra)
        self.dec = float(dec)
        self.TIC = TIC_info[0]
        self.dir = directory
        
        ### setting defaults for writeout in case the self.plot() fails in the try loop
        self.SDE = -1
        self.FAP = -1
        self.distinct_transit_count = -1
        self.empty_transit_count = -1
        self.odd_even_mismatch = -1
        self.real_tot = -999
        self.fake_tot = -999

        ### list of sectors to make sure all data downloaded correctly for this target
        outID,outEclipLong,outEclipLat,outSec,outCam,outCcd,outColPix, outRowPix, scinfo = tess_stars2px_function_entry(self.TIC, self.ra, self.dec)
        self.sectors_list = outSec

        if rerun == False:
            self.check_done()
            if self.run == 1:
                try:
                    self.file_finder()
                    if all(ele == "lightkurve" for ele in self.sources):
                        print('This target is only available in 2-minute data and will not be run')
                        self.lkonly = 1
                    else:
                        self.lkonly = 0
                        self.data_read_in()
                except Exception as e:
                    print('{} failed'.format(self.TIC))
                    print(e)
                self.writeout()
            else:
                print('This target was already run and will not be rerun.')
        else:
            try:
                self.file_finder()
                if all(ele == "lightkurve" for ele in self.sources):
                    print('This target is only available in 2-minute data and will not be run')
                    self.lkonly = 1
                else:
                    self.lkonly = 0
                    self.data_read_in()
            except Exception as e:
                print('{} failed'.format(self.TIC))
                print(e)
            self.writeout()
        
    def file_finder(self):
        os.chdir("{}".format(self.dir))
        self.ffi_dir = "eleanor"            ### could add this name in init for future use but fine now
        self.fast_dir = "lightkurve"
        filepaths = []
        sources = []
        sectors = []
        print('Searching for all downloaded files...')
        for file in glob.glob("{}/{}_*".format(self.ffi_dir,self.TIC)):
            #print(file)
            filepaths.append(self.dir+file)
            sources.append(self.ffi_dir)
            x = file.split('s')
            y = x[1].split('_')
            sectors.append(int(y[0]))
        for file in glob.glob("{}/{}_*".format(self.fast_dir,self.TIC)):
            filepaths.append(self.dir+file)
            sources.append(self.fast_dir)
            x = file.split('s')
            y = x[1].split('_')
            sectors.append(int(y[0]))
        self.filepaths = filepaths
        self.sources = sources
        self.sectors_dwn = sectors

    def data_read_in(self):
        if self.filepaths:
            starts = np.zeros(len(self.filepaths))
            ends = np.zeros(len(self.filepaths))
            self.len_filepaths = len(self.filepaths)
            for index, files in enumerate(self.filepaths):
                print(index, files)
                if self.sources[index] == self.ffi_dir:
                    ### this is 30/10 minute cadence data
                    print('Reading in 30-min data...')
                    data = fits.open('{}'.format(files))
                    cx = data[1].data['X_CENTROID']
                    cy = data[1].data['Y_CENTROID']
                    sector = data[0].header['SECTOR']
                    camera = data[0].header['CAMERA']
                    ccd = data[0].header['CCD']
                    time = data[1].data['TIME']
                    raw_flux = data[1].data['RAW_FLUX']
                    flux_err = data[1].data['FLUX_ERR']
                    quality_flags = data[1].data['QUALITY']
                    bkg_flux = data[1].data['FLUX_BKG']
                    tpfs = data[1].data['TPF']
                    biggest_ap = data[2].data['4_circle_exact']
                    used_ap_name = data[0].header['APERTURE']
                    self.used_ap = data[2].data['{}'.format(used_ap_name)]
                    #print(data[1].header)

                    ### this only applies to eleanor light curves at the moment
                    ### will have to adjust later for 2-min data or move this into eleanor only section
                    q = quality_flags == 0
                    med_image = np.nanmedian(tpfs, axis=0)
                    self.med_image = med_image
                    mask_dict = {}
                    for ind, row in enumerate(med_image):
                        for ind2, col in enumerate(row):
                            if biggest_ap[ind][ind2] == 0:
                                value = med_image[ind][ind2]
                                mask_dict[value] = (ind, ind2)
                    mask_dict = dict(sorted(mask_dict.items(), key=lambda item: item[0]))
                    lowest = list(mask_dict)[0:4]
                    first = np.zeros(4, dtype=int)
                    second = np.zeros(4, dtype=int)
                    for ilow, index3 in enumerate(lowest):
                        first[ilow] = mask_dict[index3][0]
                        second[ilow] = mask_dict[index3][1]
                    regressors = np.array([tpfs[:,0,0], tpfs[:,0,-1], tpfs[:,-1,-1], tpfs[:,-1,0], tpfs[:,first[0],second[0]], tpfs[:,first[1],second[1]], tpfs[:,first[2],second[2]], tpfs[:,first[3],second[3]]]).T
                    corr_flux = corrected_flux(data[1], time=time, flux=raw_flux, quality=quality_flags,centroid_xs=cx, centroid_ys=cy, bkg=bkg_flux, regressors=regressors, sector=sector, camera=camera, chip=ccd)
                    norm_flux = corr_flux[q]/np.nanmedian(corr_flux[q])
                    norm_ferr = flux_err[q]/np.nanmedian(corr_flux[q])

                    rolling_avg = np.zeros(len(bkg_flux[q]))
                    for ind, val in enumerate(bkg_flux[q][11:-11]):
                        rolling_avg[ind] = np.mean(bkg_flux[q][ind-11:ind+11])
                    grad = np.gradient(rolling_avg)
                    vals = np.full(shape=np.shape(grad), fill_value=3)
                    use_bool = np.less(abs(grad),vals)

                    lctime = time[q][use_bool]
                    lcflux = norm_flux[use_bool]
                    lcferr = norm_ferr[use_bool]
                
                elif self.sources[index] == self.fast_dir:
                    ### this is 2 minute/20 second data
                    data = fits.open('{}'.format(files))
                    print('Reading in 2-min data...')
                    time = data[1].data['TIME']
                    raw_flux = data[1].data['PDCSAP_FLUX']
                    flux_err = data[1].data['PDCSAP_FLUX_ERR']
                    quality_flags = data[1].data['QUALITY']

                    q = quality_flags == 0
                    normflux = raw_flux[q]/np.nanmedian(raw_flux[q])
                    normferr = flux_err[q]/np.nanmedian(raw_flux[q])

                    lctime = time[q]
                    lcflux = normflux
                    lcferr = normferr
                    #print(lctime)
                
                

                if index == 0:
                    lc = lk.LightCurve(time=lctime, flux=lcflux, flux_err=lcferr).remove_outliers(sigma_lower=7, sigma_upper=3).flatten(window_length=51,break_tolerance=15)
                    starts[index] =  lctime[0]
                    ends[index] = lctime[-1]
                else:
                    try:
                        new_lc = lk.LightCurve(time=lctime, flux=lcflux, flux_err=lcferr).remove_outliers(sigma_lower=7, sigma_upper=3).flatten(window_length=51,break_tolerance=15)
                        lc = lc.append(new_lc)
                        starts[index] =  lctime[0]
                        ends[index] = lctime[-1]
                    except Exception as e:
                        print("This sector could not be added, check to see if the file is corrupt", e)
            ### after all sectors done, put the whole lc object in self?
            self.lc = lc
            self.starts = starts
            self.ends = ends
            self.period_finder()
            tot_pts = np.sum(self.results.per_transit_count)
            depth = 1e6*(1-self.results.depth)
            if not math.isnan(self.FAP):
                if tot_pts > 9.5:
                    if (24*self.results.duration) > 0.30:
                        if self.len_filepaths < 3.5:
                            if depth > 400:
                                self.plot()
                            else:
                                print('TLS fit a depth of less than 400 ppm, unsuitable for this many sectors of data')
                        if self.len_filepaths > 3.5:
                            if depth > 200:
                                self.plot()
                            else:
                                print('TLS fit a depth of less than 200 ppm, unsuitable for this many sectors of data')
                    else:
                        print('The duration is less than 0.3 hours, will not plot.')
                else:
                    print('TLS did not actually fit any points to this "transit". Will not plot.')
            else:
                print('TLS fit this target with a FAP of NaN. Will not plot.')

    def period_finder(self):
        print('Searching for transit...')
        ### performing TLS fit of data, carry over some to plotting I guess
        periodogram = transitleastsquares(self.lc.time.value, self.lc.flux.value, self.lc.flux_err.value)
        results = periodogram.power(period_max=75.0, n_transits_min=1, oversampling_factor=2.5, show_progress_bar=False)
        self.best_period = results.period
        self.best_p_unc = results.period_uncertainty
        self.results = results

        self.SDE = results.SDE
        self.FAP = results.FAP
        self.distinct_transit_count = results.distinct_transit_count
        self.empty_transit_count = results.empty_transit_count
        self.odd_even_mismatch = results.odd_even_mismatch

        total_transit = results.transit_count
        transit_times = results.transit_times
        real_tot = 0
        fake_tot = 0
        for itt in range(total_transit):
            t_time = transit_times[itt]
            for index_s, start_val in enumerate(self.starts):
                s = self.starts[index_s]
                e = self.ends[index_s]
                if float(t_time) > float(s) and float(t_time) < float(e):
                    real_tot += 1
                else:
                    fake_tot += 1

        self.real_tot = real_tot
        self.fake_tot = fake_tot

    def plot(self):
        print('Beginning plotting...')
        ### plotting something like DV report
        ### below from ben's code
        lkf = lk.LightCurve(time=self.results.folded_phase, flux=self.results.folded_y).bin(0.01)

        in_transit = transit_mask(
            self.lc.time.value,
            self.results.period,
            self.results.duration,
            self.results.T0)

        above1 = 0
        below1 = 0
        for i in self.lc.flux.value[in_transit]:
            if i >= 1.00:
                above1 += 1
            else: 
                below1 += 1

        plt.close('all')
        fig = plt.figure(figsize=(14, 16), constrained_layout=True)
        spec = gridspec.GridSpec(ncols=2, nrows=7, figure=fig)            

        f2 = fig.add_subplot(spec[0:2, 0:1])
        transit_number = np.arange(1,len(self.results.transit_times)+1,1)
        plt.errorbar(
            transit_number,
            self.results.transit_depths,
            yerr=self.results.transit_depths_uncertainties,
            fmt='o', alpha=0.5, color='red')
        plt.plot(
            (transit_number[0], transit_number[-1]),
            (np.mean(self.results.transit_depths), np.mean(self.results.transit_depths)),
            color='black', linestyle='dashed')
        #plt.plot((self.lc.time.value.min(), self.lc.time.value.max()), (1, 1), color='black')
        plt.xlabel('Number of Orbits')
        plt.ylabel('Flux');

        f2 = fig.add_subplot(spec[0:2, 1:])
        plt.imshow(self.med_image)  
        plt.imshow(self.used_ap, cmap='Greys', alpha=0.3)

        f3 = fig.add_subplot(spec[2:3, 0:1])
        ax = plt.gca()
        ax.axvline(self.results.period, alpha=0.4, lw=3)
        plt.xlim(np.min(self.results.periods), np.max(self.results.periods))
        for n in range(2, 10):
            ax.axvline(n*self.results.period, alpha=0.4, lw=1, linestyle="dashed")
            ax.axvline(self.results.period / n, alpha=0.4, lw=1, linestyle="dashed")
        plt.ylabel('Power')
        plt.xlabel('Period (days)')
        plt.plot(self.results.periods, self.results.power, color='black', lw=0.5)
        plt.xlim(0, max(self.results.periods));

        f3 = fig.add_subplot(spec[2:3, 1:2])
        f3.set_xticks(np.array([]))
        f3.set_yticks(np.array([]))
        plt.text(0.05, 0.10, 'Period: {0:.5f} days \n\nTransit depth: {1:.1f} ppm \
            \n\nTransit duration: {2:.2f} hours \n\nFalse Alarm Prob: {3:.2f} \
            \n\nPts Above: {4} Below: {5}'.format(self.results.period,
                1e6*(1-self.results.depth),24*self.results.duration, self.FAP, above1, below1), wrap=True)
        plt.text(0.55, 0.10, 'Odd / Even var: {0:.5f} std devs \
            \n\nMean Odd: {1:.3f} Even: {2:.3f} \n\nTransits: {3:.1f} \
            \n\nPts in tranist: {4:.1f} pts \n\nDisplaying {5} of {6} Sec'.format(
                self.results.odd_even_mismatch,self.results.depth_mean_odd[0], 
                self.results.depth_mean_even[0],self.results.transit_count,
                np.sum(self.results.per_transit_count), self.len_filepaths, 
                len(self.sectors_list)),wrap=True)

        f3 = fig.add_subplot(spec[3:4, 0:1])
        plt.plot(
            self.results.model_folded_phase,
            self.results.model_folded_model,
            color='red',
            lw=5)
        plt.scatter(
            self.results.folded_phase,
            self.results.folded_y,
            color='blue',
            s=10,
            alpha=0.5,
            zorder=2)
        plt.plot(lkf.time.value, lkf.flux.value, 'k')
        duration_phase = self.results.duration / self.results.period
        plt.xlim(0.5-(2*duration_phase), 0.5+(2*duration_phase))
        plt.xlabel('Phase')
        plt.ylabel('Relative flux');

        f3 = fig.add_subplot(spec[3:4, 1:2])
        plt.plot(
            self.results.model_folded_phase,
            self.results.model_folded_model,
            color='red',
            lw=5)
        plt.scatter(
            self.results.folded_phase,
            self.results.folded_y,
            color='blue',
            s=10,
            alpha=0.5,
            zorder=2)
        plt.plot(lkf.time.value, lkf.flux.value, 'k')
        plt.xlim(-0.01, 1.01)
        plt.xlabel('Phase')
        plt.ylabel('Relative flux');

        f4 = fig.add_subplot(spec[4:5,0:])
        plt.scatter(
            self.lc.time.value[in_transit],
            self.lc.flux.value[in_transit],
            color='red',
            s=2,
            zorder=0)
        plt.scatter(
            self.lc.time.value[~in_transit],
            self.lc.flux.value[~in_transit],
            color='blue',
            alpha=0.5,
            s=2,
            zorder=0)
        plt.plot(
            self.results.model_lightcurve_time,
            self.results.model_lightcurve_model, alpha=0.5, color='red', zorder=1)
        plt.xlim(min(self.lc.time.value), max(self.lc.time.value))
        plt.ylim(min(self.lc.flux.value), max(self.lc.flux.value))
        plt.xlabel('Time (days)')
        plt.ylabel('Relative flux');

        f4 = fig.add_subplot(spec[5:6,0:])
        plt.plot(lkf.time.value, lkf.flux.value, 'k')
        #print(self.results.depth_mean_even, self.results.depth_mean_odd)

        sbool_list = np.zeros(len(self.sectors_list), dtype=bool)
        for isec,jsec in enumerate(self.sectors_list):
            if jsec in self.sectors_dwn:
                sbool_list[isec] = True

        ssec = []
        esec=[]
        for ibool, jbool in enumerate(sbool_list):
            if jbool:
                ssec.append(ibool)
                esec.append(ibool+1)
        #print(ssec)
        #print(sbool_list)
        #print(self.sectors_dwn, self.sectors_list)
        f5 = fig.add_subplot(spec[6:7,0:])
        plt.xticks(ticks=np.arange(len(self.sectors_list)),labels=self.sectors_list)
        plt.xlim(-0.5,len(self.sectors_list)-0.5)
        plt.yticks([])
        for ist,jst in enumerate(ssec):
            xs = [ssec[ist]-0.5,esec[ist]-0.5]
            plt.fill_between(xs,1,color='green',alpha=0.7)


        plt.savefig(self.dir + 'plots_2_12/' + str(self.TIC) + '.png')

    def writeout(self):
        print('Writing out results...')
        totalfile = open('/home/z5318114/targets_plotted.txt',"a")
        totalfile.write(str(self.TIC)+'\t'+str(self.ra)+'\t'+str(self.dec)+'\t'+str(self.SDE)+'\t'+str(self.FAP)
            +'\t'+str(self.distinct_transit_count)+'\t'+str(self.empty_transit_count)+'\t'+str(self.odd_even_mismatch)
            +'\t'+str(self.lkonly)+'\t'+str(self.real_tot)+'\t'+str(self.fake_tot)+'\n')
        totalfile.close()    
        '''
        writeout_specific = open('/home/z5318114/394584021_data.txt',"a")
        for i,j in enumerate(self.lc.time.value):
            writeout_specific.write(str(self.lc.time.value[i])+'\t'+str(self.lc.flux.value[i])+'\t'+str(self.lc.flux_err.value[i])+'\t'+str(self.results.model_lightcurve_model[i])+'\n')
        writeout_specific.close()
        '''

    def check_done(self):
        print('Checking if target has been run previously...')
        done_id = []
        sdes = []
        with open('/home/z5318114/targets_plotted.txt','r') as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                done_id.append(row[0])
                sdes.append(row[3])

        ### starting index check at 110 because a bunch of early targets might have been run several times
        if self.TIC in done_id and sdes[done_id.index(self.TIC, 110)] == '-1':
            self.run = 0
        else:
            self.run = 1
