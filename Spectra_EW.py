import os
import sys
import csv
import numpy as np
from scipy import signal
from astropy.io import fits
from scipy import interpolate
from astropy import units as u
import matplotlib.pyplot as plt
from specutils import Spectrum1D
from astropy.stats import sigma_clip
from specutils.spectra import SpectralRegion
from specutils.analysis import equivalent_width

#fibres = np.arange(1, 401, 1)
#fibres_str = []
#for i in fibres:
#    fibres_str.append(str(i).zfill(3))

class target:
    """
    parameters:
        dir: directory where the exposures are kept, ex. '/home/ryan/Downloads/200708/data/''
        exp: list of exposures to use i.e. [0031, 0032, 0033]
        fiber: fiber of the target (~1-400)
        ccd: which # ccd to use (1-4 for galah)
        date: format DDmmm, i.e. 08jul
        coadd: target/exposure that coadded spectra will be from
        beginning: index of first file in folder to analyze
        end: index of last file in folder to analyz
    """

    def __init__(self, date, oversample_factor, beginning):#, end): #dir, ccd, date, exp, fiber, oversample_factor, coadd=None):    
        ### just date, oversample factor, maybe file w/ exposures
        ### skip end if running all files in folder/list
        #self.fiber = fiber
        self.date = date
        self.oversample = oversample_factor
        #self.exp = exp
        #self.dir = os.path.join(dir, 'ccd_{}/'.format(ccd))
        #self.ccd = ccd

        '''
        ### Opening the field logs for GALAH obs up to 2020-10-07
        ### Use date to match date from inputs
        galah_field = []
        galah_obs_date = []
        galah_obs_ra = []
        galah_obs_dec = []
        galah_exp = []
        with open('/home/rmorris/documents/fieldlogs_export.txt', 'r') as file:
            reader = csv.reader(file,delimiter='\t')
            next(reader)
            for row in reader:
                galah_field.append(row[0])
                galah_obs_date.append(row[1])
                galah_obs_ra.append(row[2])
                galah_obs_dec.append(row[3])
                galah_exp.append(row[4])
                ### these exposures formatted like 'XX-YY', need [00XX, 00XX+1, ..., 00XX+N, 00YY]
                
        galah_exp_lists = []
        for i,j in enumerate(galah_exp):
        ### converting exposure format to [XX, XX+1, XX+2, ..., YY]
            exp_string = galah_exp[i]
            print(exp_string)
            ind_exp = exp_string.split('-')
            print(ind_exp)
            try:
                exp_arr = np.arange(int(ind_exp[0]), int(ind_exp[-1])+1, 1)
            except Exception as e:
                print(e)
                try:
                    ind_end = 
                exp_arr = ind_exp[0]
            print(exp_arr)
            galah_exp_lists.append(exp_arr)
        
        exp_match = []
        for i, j in enumerate(galah_obs_date):
            if self.date == galah_obs_date[i]:
                ### finds the rows with the relevant fields for an input date
                exp_match.append(galah_exp_lists[i])
                ### exp_match has all series of exposures on a given night, ex. [[31,32,33],[36,37,38],[42,43,44]]
        '''           
        
        com_files = os.listdir('/data/wallaby/rmorris/GALAH/Simult_GALAH_Data/reductions/results/{}/spectra/com/'.format(self.date))
        com_dict = {"test gaia": "test file"}
        for file_name in com_files:
            com = fits.open('/data/wallaby/rmorris/GALAH/Simult_GALAH_Data/reductions/results/{}/spectra/com/{}'.format(self.date, file_name))  
            gaia = com[1].header['GAIA_ID']
            com_dict['{}'.format(gaia)] = file_name
            com.close()
        
        fits_files = os.listdir('/data/wallaby/rmorris/GALAH/Simult_GALAH_Data/reductions/results/{}/spectra/all/'.format(self.date))
        ones = [s for s in fits_files if "1.fits" in s]
        threes = [s for s in fits_files if "3.fits" in s]
        ew_files = ones+threes
        flux_arr = [] 
        wave_arr = []   
        error_arr = []    
        log_arr = []
        time_arr = []
        RV_arr = []
        RV_err_arr = []
        v_broad_arr = []
        ccd_arr = []
        file_names = []
        gaia_ids = []
        self.beg_ind = beginning
        #self.end_ind = end
        medians = []
        for file_name in ew_files[self.beg_ind::]:#self.end_ind]:            
            #full_path = os.path('/data/wallaby/rmorris/GALAH/Simult_GALAH_Data/reductions/results/{}/spectra/all/{}'.format(self.date, file_name))
            print(file_name)
            ccd_proxy = file_name[-6]       #should be last digit of file before .fits
            hdu = fits.open('/data/wallaby/rmorris/GALAH/Simult_GALAH_Data/reductions/results/{}/spectra/all/{}'.format(self.date, file_name), memmap=False)
            gaia = hdu[1].header['GAIA_ID']
            com = fits.open('/data/wallaby/rmorris/GALAH/Simult_GALAH_Data/reductions/results/{}/spectra/com/{}'.format(self.date, com_dict['{}'.format(gaia)]), memmap=False)          
            flux = hdu[1].data       # normalized flux straight from GALAH fits files
            error = hdu[2].data      # 'relative error of flux or norm flux'
            pix = hdu[0].header['CRPIX1']
            val = hdu[0].header['CRVAL1']
            delt = hdu[0].header['CDELT1']
            median = np.median(flux)
            if com[0].header['RV_OK'] == 1:
                RV_shift = com[0].header['RV']
                RV_error = com[0].header['E_RV']
            else:
                print('RVs are not well characterized')
                RV_shift = 0
                RV_error = 0
            if com[0].header['PAR_OK'] == 1:
                v_broad = com[0].header['VBROAD_R']
            else:
                print('Parameters are not well characterized')
                v_broad = 0
            time_arr.append(hdu[0].header['UTMJD'])
            wave = np.arange(0, len(flux), dtype='float64')
            for inds in range(len(wave)):
                if inds < float(pix)-1:
                    wave[inds] = float(val)-float(delt)*((float(pix)-1)-inds)
                elif inds > float(pix)-1:
                    wave[inds] = float(val)+float(delt)*(inds-(float(pix)-1))
                else:
                    wave[inds] = float(val)
            hdu.close()
            com.close()
            del hdu[0].data
            del hdu[1].data
            del hdu[2].data
            del com[0].data
            
            logwave = np.log(wave)
            flux_arr.append(flux)
            wave_arr.append(wave)
            error_arr.append(error) 
            log_arr.append(logwave)
            RV_arr.append(RV_shift)
            RV_err_arr.append(RV_error)
            v_broad_arr.append(v_broad)
            ccd_arr.append(ccd_proxy)
            file_names.append(file_name)
            gaia_ids.append(gaia)
            medians.append(median)

        self.flux = flux_arr
        self.wave = wave_arr
        self.error = error_arr
        self.logwave = log_arr
        self.time = time_arr
        self.base_rv = RV_arr
        self.base_rv_err = RV_err_arr
        self.vbroad = v_broad_arr
        self.ccd = ccd_arr
        self.file_name = file_names
        self.gaia_ids = gaia_ids
        self.medians = medians

        #self.index_finder()
        #if coadd is None:
            #print('Developing Coadded Spectrum')
            #self.coadd_spectra()
        #elif coadd is not None:
            #print('Using Alternate Coadded Spectrum')
            #self.logcowave = coadd.logcowave
            #self.nancoindex = coadd.nancoindex
            #self.coadd = coadd.coadd
        #self.polyfitter()
        #self.interpolate_spec()
        #self.cross_correlate()
        #self.subpixmax()
        #self.rvel()
        #self.residual()
        #for index in range(len(self.flux)):
            #if int(self.ccd[index]) == 1 or 3:
        self.equivalent_width()
        self.writeout()
    
    def index_finder(self):
        ### input list of flux and waves to output index of unmasked points and zero centered wave
        out_ind_arr = []
        #out_w_arr = []
        #out_w_med_val_arr = []
        for index in range(len(self.flux)):
            flux = self.flux[index]
            logwave = self.logwave[index]
            fflux = sigma_clip(flux, sigma=2, masked=True)

            num4filt = 5
            ind_list = []
            for ind in range(len(flux)):
                for x in np.arange(-num4filt, num4filt+1, 1):
                    if ind+x > 4095:
                        newind = 4095
                    elif ind-x <0:
                        newind = 0
                    else:
                        newind = ind+x
                    if ind not in ind_list:
                        if np.ma.is_masked(flux[newind]) == True:
                            ind_list.append(ind)

            false = np.zeros(len(fflux), dtype='bool')
            for ind in ind_list:
                false[ind] = True

            masked = np.ma.masked_where(false == True, fflux)
            idx = np.zeros(len(masked), dtype='bool')
            for act_id in range(len(idx)):
                idx[act_id] = (not np.ma.is_masked(masked[act_id]) & np.isfinite(logwave[act_id]))

            #wave_med_val = np.nanmedian(wave)  #don't need anymore with log wave
            #wave -= wave_med_val
            #log_wave = np.log(wave)             # turning wavelength into log(wavelength) space
            out_ind_arr.append(idx)
            #out_w_arr.append(log_wave)
            #out_w_med_val_arr.append(wave_med_val)
        #return(out_ind_arr)#, out_w_arr, out_w_med_val_arr) #don't need w_arr or w_med anymore
        self.nanindex = out_ind_arr

    def coadd_spectra(self):#, flux_list, error_list, wave):
        flux_sum = np.zeros(len(self.flux[0]))
        error_sum = np.zeros(len(self.flux[0]))
        #n = range(len(flux_list))
        #i = range(len(flux_list[0]))
        for i in range(len(self.flux[0])):
            for n in range(len(self.flux)):
                f_by_e = self.flux[n][i] / self.error[n][i]
                flux_sum[i] += f_by_e
                error_sum[i] += (1 / self.error[n][i])
        coadd_flux = flux_sum / error_sum
        self.logcowave = self.logwave[0]#[~np.isnan(coadd_flux)]
        self.nancoindex = self.nanindex[0]#[~np.isnan(coadd_flux)]
        coadd_flux = coadd_flux#[~np.isnan(coadd_flux)]          #chops of the nans at the end
        self.coadd = coadd_flux
        print(np.nanstd(coadd_flux))
        #return(coadd_flux)

    def polyfitter(self, deg=4):
        vals_arr = []
        normalize_arr = []
        norm_co_arr = []
        logwave = self.logcowave#[0]#[0:len(self.coadd)]
        index = self.nancoindex#[0]#[0:len(self.coadd)]
        polyfit = np.polynomial.polynomial.polyfit(logwave[index], self.coadd[index], deg)
        for i in range(len(self.flux)):
            wave = self.logwave[i]
            index = self.nanindex[i]
            flux = self.flux[i]
            vals = np.polynomial.polynomial.polyval(wave, polyfit)
            norm_flux = flux/vals
            norm_flux /= np.nanmedian(norm_flux)
            cos_flux = np.ma.masked_where(norm_flux > np.nanmedian(norm_flux)+(2*np.nanstd(norm_flux)), norm_flux)  #getting rid of cosmic rays
            zero_flux = cos_flux-1
            normalize_arr.append(zero_flux)
            vals_arr.append(vals)

        cowave = self.logcowave
        coindex = self.nancoindex
        coflux = self.coadd
        covals = np.polynomial.polynomial.polyval(cowave, polyfit)
        conorm_flux = coflux/covals
        conorm_flux /= np.nanmedian(conorm_flux)
        cocos_flux = np.ma.masked_where(conorm_flux > np.nanmedian(conorm_flux)+(2*np.nanstd(conorm_flux)), conorm_flux)  #getting rid of cosmic rays
        cozero_flux = cocos_flux-1

        self.normcoadd = cozero_flux
        self.normflux = normalize_arr
        self.polyfit = polyfit      #polyfit coefficients for coadd spectra. only one per coadded spectra
        self.polyval = vals_arr

    ### make sure to interpolate only the non-masked points, should be good now
    def interpolate_spec(self):#flux_list, wave):
        interp_flux_list = []
        interp_flux_reg_list = []
        #interp_co_flux_list = []
        logwave = self.logwave[0]     # probably needs to change
        wave = self.wave[0]
        interp_wave_reg_list = []
        interp_wave_list = []
        for index in range(len(self.normflux)):
            flux = self.flux[index]     #was normflux before using GALAH reductions
            #idx = np.zeros(len(self.normflux[index]), dtype=bool)
            interp_wave = np.arange(logwave[0], logwave[-1], (1/self.oversample), dtype='float64')
            #for act_id in range(len(idx)):
            #    idx[act_id] = (not np.ma.is_masked(flux[act_id]) & np.isfinite(logwave[act_id]))
                
            #idx_reg = np.zeros(len(self.normflux[index]), dtype=bool)
            #interp_wave_reg = np.arange(wave[0], wave[-1], (1/self.oversample), dtype='float64')
            #for act_id in range(len(idx_reg)):
            #    idx_reg[act_id] = (not np.ma.is_masked(flux[act_id]) & np.isfinite(wave[act_id]))
                
            #func_reg = interpolate.interp1d(wave[idx_reg], flux[idx_reg], bounds_error = False, fill_value=0)
            func = interpolate.interp1d(logwave[idx], flux[idx], bounds_error = False, fill_value=0)
            #interp_flux_reg = func_reg(interp_wave_reg)
            interp_flux = func(interp_wave)
            #where_nans_reg = np.isnan(interp_flux_reg)
            #where_nans = np.isnan(interp_flux)
            #interp_flux_reg[where_nans_reg] = 0
            #interp_flux[where_nans] = 0
            #interp_flux_reg_list.append(interp_flux_reg)
            #interp_wave_reg_list.append(interp_wave_reg)
            interp_flux_list.append(interp_flux)
            interp_wave_list.append(interp_wave)

        #coflux = self.normcoadd
        #cowave = self.logcowave
        #print(len(wave), len(cowave))
        #coidx = np.zeros(len(self.normcoadd), dtype=bool)
        #interp_co_wave = np.arange(cowave[0], cowave[-1], (1/self.oversample), dtype='float64')
        #for act_id in range(len(coidx)):
        #    coidx[act_id] = (not np.ma.is_masked(coflux[act_id]) & np.isfinite(cowave[act_id]))
        #cofunc = interpolate.interp1d(cowave[coidx], coflux[coidx], bounds_error = False, fill_value=0)#"extrapolate")
        #cointerp_flux = cofunc(interp_co_wave)
        #where_conans = np.isnan(cointerp_flux)
        #cointerp_flux[where_conans] = 0
        #interp_co_wave = interp_co_wave[~np.isnan(cointerp_flux)]
        #cointerp_flux = cointerp_flux[~np.isnan(cointerp_flux)]

        #self.interpcoflux = cointerp_flux
        #self.interpcowave = interp_co_wave
        #self.interpfluxreg = interp_flux_reg_list
        #self.interpwavereg = interp_wave_reg_list
        self.interpflux = interp_flux_list
        self.interpwave = interp_wave
        self.interpwavelist = interp_wave_list
        #return(interp_flux_list, interp_wave)

    def cross_correlate(self):#flux_list):
        ### correlate an array of fluxes where last flux is the one correlating to
        correlate_list = []
        index_list = []
        new_ind_list= []
        max_val_ind = []
        #auto_max_ind = []
        for index in range(len(self.interpflux)):
            #print(len(self.interpflux[index]), len(self.interpcoflux))
            correlate = signal.correlate(self.interpflux[index], self.interpcoflux)
            #auto_corr = signal.correlate(self.interpflux[index], self.interpflux[index])
            #auto_corr /= np.nanmax(auto_corr)
            correlate /= np.nanmax(correlate)
            ind_where = np.argmax(correlate)
            #auto_ind_where = np.argmax(auto_corr)
            max_val_ind.append(ind_where)
            #auto_max_ind.append(auto_ind_where)
            correlate_list.append(correlate)
            inds = np.arange(len(correlate))
            new_inds = ((len(correlate)/2) + 0.5 - inds) * 1/self.oversample * 299792458
            new_ind_list.append(new_inds)
            index_list.append(inds)

        self.maxindex = max_val_ind     
        self.correlate = correlate_list
        self.corrindex = index_list#new_ind_list#
        #self.autoindex = auto_max_ind
        #return(correlate_list, index_list, max_val_ind)

    def subpixmax(self):#correlate_list, index_list, max_val_ind_list):
        act_list = []
        delta_list = []
        logshift_list = []
        for index in range(len(self.correlate)):
            middle_ind = self.maxindex[index]
            corr = self.correlate[index]
            inds = self.corrindex[index]
            corr_pts = corr[middle_ind-1:middle_ind+2]
            ind_pts = inds[middle_ind-1:middle_ind+2]
            parabola = np.polynomial.polynomial.polyfit(ind_pts, corr_pts, 2)       #polyfitting the 3 points
            derivative = np.polynomial.polynomial.polyder(parabola)                 #derivative of parabola
            zeros = np.polynomial.polynomial.polyroots(derivative)                  #finds zeros of function
            corr_max_ind = zeros[0]                               #actual decimal index of correlate with max val
            #delta = self.autoindex[index] - corr_max_ind
            delta = (len(self.correlate[index])/2) + 0.5 - corr_max_ind
            print(corr_max_ind, delta, (len(self.correlate[index])/2)+0.5, len(self.correlate[index]), ind_pts)
            logshift = delta * 1/self.oversample
            act_list.append(corr_max_ind)
            delta_list.append(delta)
            logshift_list.append(logshift)
        self.pixelshift = act_list
        self.pixeldelta = delta_list
        self.logshift = logshift_list
        #return(act_list, delta_list)

    def rvel(self):#delta_list, wave_list):
        c = 299792458   #m/s   or 299792.458 km/s
        velocity_list = []
        for index in range(len(self.logshift)):
            delta = self.logshift[index]
            velocity = (delta)*c
            velocity_list.append(velocity)
        self.velocity = velocity_list
        
    def residual(self):
        residual_list = []
        for index in range(len(self.interpflux)):
            interp_flux = self.interpflux[index]
            coadd_flux = self.interpcoflux
            residual = np.zeros(len(interp_flux))
            for inds in range(len(interp_flux)):
                residual[inds] = (interp_flux[inds]+1) / (coadd_flux[inds]+1)
            residual_list.append(residual)
        self.residual = residual_list
           
    def equivalent_width(self):
        print('Calculating Equivalent Width')
        ew_list = []
        ew_error_list = []
        halpha = 6562.8 # angstroms
        hbeta = 4861.3 # angstroms
        for index in range(len(self.flux)):
            error_sq = []
            flux_vals = []     # just for plotting relevant sections
            wave_vals = []     # just for plotting relevant sections
            ccd = self.ccd[index]
            print(ccd)
            
            if int(ccd) == 1 or int(ccd) == 3:

                ### this block takes into account the radial velocity shift on all the interpolated wavelengths
                rv = self.base_rv[index]
                center_line_alpha = halpha + (halpha*rv/299792.458)
                center_line_beta = hbeta + (hbeta*rv/299792.458)

                vbroad = self.vbroad[index]
                alpha_broad_r = center_line_alpha + 9*(center_line_alpha*vbroad/299792.458)
                alpha_broad_l = center_line_alpha - 9*(center_line_alpha*vbroad/299792.458)
                beta_broad_r = center_line_beta + 9*(center_line_beta*vbroad/299792.458)
                beta_broad_l = center_line_beta - 9*(center_line_beta*vbroad/299792.458)

                if int(ccd) == 1:
                    center_line = center_line_beta
                    left = beta_broad_l
                    right = beta_broad_r
                elif int(ccd) == 3:
                    center_line = center_line_alpha
                    left = alpha_broad_l
                    right = alpha_broad_r

                for ind in range(len(self.flux[index])):  
                    if self.wave[index][ind] >= left and self.wave[index][ind] <= right:
                        flux_vals.append(self.flux[index][ind])
                        wave_vals.append(self.wave[index][ind])
                        error_sq.append(float(self.error[index][ind])**2)
                ew_est = 0    # comparing EW estimation to EW w/ spec utils
                center_est = []
                print(wave_vals)
                for i, fval in enumerate(flux_vals[0:-1]):
                    center = (fval + flux_vals[i+1])/2
                    center_est.append(center)
                    riemann_r = center*(wave_vals[1]-wave_vals[0])
                    ew_est += riemann_r
                ew_error = np.sqrt(np.sum(error_sq))#*u.AA
                ew_list.append((right-left)-ew_est)
                ew_error_list.append(ew_error)
                '''
                print(wave_vals[0:-1])
                plt.plot(self.wave[index],self.flux[index])
                plt.scatter(wave_vals[0:-1], center_est)
                plt.vlines(right, -1,2, colors='black')
                plt.vlines(left, -1,2, colors='black')
                plt.vlines(center_line, -1,2, colors='red')
                plt.xlim(left-2, right+2)
                plt.xlabel('Angstroms')
                plt.ylabel('Normalized Flux')
                plt.title('{}'.format(self.file_name[index]))
                plt.savefig('/data/wallaby/rmorris/GALAH/EW_Plots/{}.png'.format(self.file_name[index][0:-5]))
                plt.close()
                '''
            else:
                ew = 0
                ew_error = 0
                ew_list.append(ew)
                ew_error_list.append(ew_error)
            
        self.equivwidth = ew_list
        self.equiverror = ew_error_list
        
    def writeout(self):
        file = open('/home/rmorris/documents/spectra_ew_all.txt', "a")
        for number, values in enumerate(self.gaia_ids[self.beg_ind::]):#self.end_ind]):
            file.write(str(self.gaia_ids[number])+'\t'+str(self.file_name[number])+'\t'+str(self.equivwidth[number])+'\t'+str(self.equiverror[number])+'\t'+str(self.medians[number])+'\n')
        file.close()
            
        #with open('/home/rmorris/documents/spectra_ew_{}-{}.txt'.format(self.beg_ind, self.end_ind), 'w') as file:
            #for number, values in enumerate(self.gaia_ids[self.beg_ind:self.end_ind]):
                #file.write(str(self.gaia_ids[number])+'\t'+str(self.file_name[number])+'\t'+str(self.equivwidth[number])+'\t'+str(self.equiverror[number])+'\t'+str(self.medians[number])+'\n')
            

#test = target(201003, 10**4, 0)
#181220, 
#dates = [181221, 181222, 181223, 181224,181225,190204,190210,190212,190223,190224,200708,200712,200714,200724,200728,200802,200825,200826,200901,200906]
#for i in dates:
#    target(i, 10**4, 0)
#test.equivalent_width()
#print(test.equivwidth[0]*u.AA, test.equiverror)






