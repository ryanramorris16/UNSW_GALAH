import os
import csv
import glob
import numpy as np
import lightkurve as lk
from eleanor import mast
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from transitleastsquares import transit_mask
from target_data_functions import corrected_flux
from transitleastsquares import transitleastsquares
from transitleastsquares.stats import calculate_stretch
from transitleastsquares.stats import all_transit_times
from transitleastsquares.transit import fractional_transit
import transitleastsquares.tls_constants as tls_constants



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

        if rerun == False:
            self.check_done()
            if self.run == 1:
                #try:
                self.file_finder()
                self.data_read_in()
                #except exception as e:
                #    print('{} failed'.format(self.TIC))
                #    print(e)
                self.writeout()
            else:
                print('This target was already run and will not be rerun.')
        else:
            #try:
            self.file_finder()
            self.data_read_in()
            #except:
            #    print('{} failed'.format(self.TIC))
            self.writeout()
        
    def file_finder(self):
        os.chdir("{}".format(self.dir))
        self.ffi_dir = "eleanor"            ### could add this name in init for future use but fine now
        self.fast_dir = "lightkurve"
        filepaths = []
        sources = []
        print('doing file_finder')
        for file in glob.glob("{}/{}_*".format(self.ffi_dir,self.TIC)):
            print(file)
            filepaths.append(self.dir+file)
            sources.append(self.ffi_dir)
        for file in glob.glob("{}/{}_*".format(self.fast_dir,self.TIC)):
            filepaths.append(self.dir+file)
            sources.append(self.fast_dir)
        self.filepaths = filepaths
        self.sources = sources

    def data_read_in(self):
        if self.filepaths:
            for index, files in enumerate(self.filepaths):
                print(index, files)
                if self.sources[index] == self.ffi_dir:
                    ### this is 30/10 minute cadence data
                    print('doing 30-min data')
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

                    rolling_avg = np.zeros(len(bkg_flux[q]))
                    for ind, val in enumerate(bkg_flux[q][11:-11]):
                        rolling_avg[ind] = np.mean(bkg_flux[q][ind-11:ind+11])
                    grad = np.gradient(rolling_avg)
                    vals = np.full(shape=np.shape(grad), fill_value=3)
                    use_bool = np.less(abs(grad),vals)

                    lctime = time[q][use_bool]
                    lcflux = norm_flux[use_bool]
                    lcferr = flux_err[q][use_bool]
                
                elif self.sources[index] == self.fast_dir:
                    ### this is 2 minute/20 second data
                    ### not sure whether to use SAP or PDC because tpfs?
                    ### could just skip and set corr flux to PDC later
                    data = fits.open('{}'.format(files))
                    print('doing 2-min data')
                    time = data[1].data['TIME']
                    raw_flux = data[1].data['PDCSAP_FLUX']
                    flux_err = data[1].data['PDCSAP_FLUX_ERR']
                    quality_flags = data[1].data['QUALITY']

                    q = quality_flags == 0
                    normflux = raw_flux[q]/np.nanmedian(raw_flux[q])

                    lctime = time[q]
                    lcflux = normflux
                    lcferr = flux_err[q]
                
                

                if index == 0:
                    lc = lk.LightCurve(time=lctime, flux=lcflux, flux_err=lcferr).remove_outliers(sigma_lower=7, sigma_upper=3).flatten(window_length=51,break_tolerance=15)
                else:
                    new_lc = lk.LightCurve(time=lctime, flux=lcflux, flux_err=lcferr).remove_outliers(sigma_lower=7, sigma_upper=3).flatten(window_length=51,break_tolerance=15)
                    lc = lc.append(new_lc)
            ### after all sectors done, put the whole lc object in self?
            self.lc = lc
            self.period_finder()
            self.plot()

    def period_finder(self):
        print('Searching for transit...')
        ### performing TLS fit of data, carry over some to plotting I guess
        periodogram = transitleastsquares(self.lc.time.value, self.lc.flux.value, self.lc.flux_err.value)
        results = periodogram.power(period_max=75.0, n_transits_min=1, oversampling_factor=5)
        self.best_period = results.period
        self.best_p_unc = results.period_uncertainty
        self.results = results

        self.SDE = results.SDE
        self.FAP = results.FAP
        self.distinct_transit_count = results.distinct_transit_count
        self.empty_transit_count = results.empty_transit_count
        self.odd_even_mismatch = results.odd_even_mismatch

        maxwidth_in_samples = int(np.max(periodogram.durations) * np.size(periodogram.y))
        transit_times = all_transit_times(results.T0, periodogram.t, results.period)
        stretch = calculate_stretch(periodogram.t, results.period, transit_times)
        internal_samples = (
                int(len(periodogram.y)) * tls_constants.OVERSAMPLE_MODEL_LIGHT_CURVE)
        ### these are folded models
        self.results.even_model = fractional_transit(duration=(results.duration * maxwidth_in_samples),
                maxwidth=maxwidth_in_samples / stretch,
                depth=1 - results.depth_mean_even[0],
                samples=internal_samples,
                per=periodogram.per,
                rp=periodogram.rp,
                a=periodogram.a,
                inc=periodogram.inc,
                ecc=periodogram.ecc,
                w=periodogram.w,
                u=periodogram.u,
                limb_dark=periodogram.limb_dark,
            )
        self.results.odd_model = fractional_transit(duration=(results.duration * maxwidth_in_samples),
                maxwidth=maxwidth_in_samples / stretch,
                depth=1 - results.depth_mean_odd[0],
                samples=internal_samples,
                per=periodogram.per,
                rp=periodogram.rp,
                a=periodogram.a,
                inc=periodogram.inc,
                ecc=periodogram.ecc,
                w=periodogram.w,
                u=periodogram.u,
                limb_dark=periodogram.limb_dark,
            )
        step = (self.results.model_folded_phase[-1]-self.results.model_folded_phase[0])/len(self.results.even_model)
        self.times = np.arange(self.results.model_folded_phase[0], self.results.model_folded_phase[-1], step)


    def plot(self):
        print('Beginning plotting...')
        ### plotting something like DV report
        ### below from ben's code
        lkf = lk.LightCurve(self.results.folded_phase, self.results.folded_y).bin(0.01)

        plt.close('all')
        fig = plt.figure(figsize=(14, 16), constrained_layout=True)
        spec = gridspec.GridSpec(ncols=2, nrows=6, figure=fig)            

        f2 = fig.add_subplot(spec[0:2, 0:1])
        plt.plot(
            self.times,
            self.results.even_model,
            color='red',
            lw=5)
        plt.plot(
            self.times,
            self.results.odd_model,
            color='yellow',
            lw=5)
        plt.scatter(
            self.results.folded_phase,
            self.results.folded_y,
            color='blue',
            s=10,
            alpha=0.5,
            zorder=2)
        plt.plot(lkf.time.value, lkf.flux.value, 'k')
        plt.xlim(0.46, 0.54)
        plt.xlabel('Phase')
        plt.ylabel('Relative flux');   

        f2 = fig.add_subplot(spec[0:2, 1:])
        plt.imshow(self.med_image)  

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
        plt.text(0.05, 0.75, 'Period: {0:.5f} days'.format(self.results.period))
        plt.text(0.05, 0.50, 'Transit depth: {0:.1f} ppm'.format(1e6*(1-self.results.depth)))
        plt.text(0.05, 0.25, 'Transit duration: {0:.2f} days'.format(24*self.results.duration))
        plt.text(0.55, 0.75, 'Odd / Even var: {0:.5f} std devs'.format(self.results.odd_even_mismatch))
        plt.text(0.55, 0.50, 'Transits: {0:.1f}'.format(self.results.transit_count))
        plt.text(0.55, 0.25, 'Pts in tranist: {0:.1f} pts'.format(np.sum(self.results.per_transit_count)))

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
        plt.xlim(0.46, 0.54)
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
        in_transit = transit_mask(
            self.lc.time.value,
            self.results.period,
            self.results.duration,
            self.results.T0)
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
        print(self.results.depth_mean_even, self.results.depth_mean_odd)

        plt.savefig(self.dir + 'plots/' + str(self.TIC) + '.png')

    def writeout(self):
        print('Writing out results...')
        totalfile = open('/home/z5318114/targets_plotted.txt',"a")
        totalfile.write(str(self.TIC)+'\t'+str(self.ra)+'\t'+str(self.dec)+'\t'+str(self.SDE)+'\t'+str(self.FAP)+'\t'+str(self.distinct_transit_count)+'\t'+str(self.empty_transit_count)+'\t'+str(self.odd_even_mismatch)+'\n')
        totalfile.close()    

    def check_done(self):
        print('Checking if target has been run previously...')
        done_id = []
        with open('/home/z5318114/targets_plotted.txt','r') as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                done_id.append(row[0])

        if self.TIC in done_id:
            self.run = 0
        else:
            self.run = 1
