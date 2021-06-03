###copy of eb_spectra for testing

import os
import sys
import numpy as np
from scipy import signal
from astropy.io import fits
from scipy import interpolate
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip


class target:
    '''
    parameters:
        dir: directory where the exposures are kept, ex. '/home/ryan/Downloads/200708/data/''
        exp: list of exposures to use i.e. [0031, 0032, 0033]
        fiber: fiber of the target (~1-400)
        ccd: which # ccd to use (1-4 for galah)
        date: format DDmmm, i.e. 08jul
    '''

    def __init__(self, dir, ccd, date, exp, fiber, oversample_factor):
        self.fiber = fiber
        self.oversample = oversample_factor
        self.exp = exp
        self.dir = os.path.join(dir, 'ccd_{}/'.format(ccd))

        flux_arr = [] 
        wave_arr = []   
        error_arr = []    
        log_arr = []
        for exposure in exp:
            full_path = os.path.join(self.dir, '{}{}{}red.fits'.format(date,ccd,exposure))
            hdu = fits.open(full_path)
            flux = hdu[0].data[self.fiber]
            error = hdu[1].data[self.fiber]
            pix = hdu[0].header[167]
            val = hdu[0].header[165]
            delt = hdu[0].header[166]
            wave = np.arange(0, len(flux), dtype='float64')
            for inds in range(len(wave)):
                if inds < float(pix)-1:
                    wave[inds] = float(val)-float(delt)*((float(pix)-1)-inds)
                elif inds > float(pix)-1:
                    wave[inds] = float(val)+float(delt)*(inds-(float(pix)-1))
                else:
                    wave[inds] = float(val)
            logwave = np.log(wave)
            flux_arr.append(flux)
            wave_arr.append(wave)
            error_arr.append(error)
            log_arr.append(logwave)
                         #want something like self.exp1, self.exp2, self.exp3 ... self.expN

        self.flux = flux_arr
        self.wave = wave_arr
        self.error = error_arr
        self.logwave = log_arr

        self.index_finder()
        self.coadd_spectra()
        self.polyfitter()
        self.interpolate_spec()
        self.cross_correlate()
        self.subpixmax()
        self.rvel()

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
        self.logcowave = self.logwave[0][~np.isnan(coadd_flux)]
        self.nancoindex = self.nanindex[0][~np.isnan(coadd_flux)]
        coadd_flux = coadd_flux[~np.isnan(coadd_flux)]          #chops of the nans at the end
        self.coadd = coadd_flux
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
        interp_co_flux_list = []
        for index in range(len(self.normflux)):
            flux = self.normflux[index]
            wave = self.logwave[index]
            idx = np.zeros(len(self.normflux[index]), dtype=bool)
            interp_wave = np.arange(wave[0], wave[-1], (1/self.oversample), dtype='float64')
            for act_id in range(len(idx)):
                idx[act_id] = (not np.ma.is_masked(flux[act_id]) & np.isfinite(wave[act_id]))

            func = interpolate.interp1d(wave[idx], flux[idx], fill_value="extrapolate")
            interp_flux = func(interp_wave)
            interp_wave = interp_wave[~np.isnan(interp_flux)]
            interp_flux = interp_flux[~np.isnan(interp_flux)]
            interp_flux_list.append(interp_flux)

        coflux = self.normcoadd
        cowave = self.logcowave
        coidx = np.zeros(len(self.normcoadd), dtype=bool)
        interp_co_wave = np.arange(cowave[0], cowave[-1], (1/self.oversample), dtype='float64')
        for act_id in range(len(coidx)):
            coidx[act_id] = (not np.ma.is_masked(coflux[act_id]) & np.isfinite(cowave[act_id]))
        cofunc = interpolate.interp1d(cowave[coidx], coflux[coidx], fill_value="extrapolate")
        cointerp_flux = cofunc(interp_co_wave)
        interp_co_wave = interp_co_wave[~np.isnan(cointerp_flux)]
        cointerp_flux = cointerp_flux[~np.isnan(cointerp_flux)]

        self.interpcoflux = cointerp_flux
        self.interpcowave = interp_co_wave
        self.interpflux = interp_flux_list
        self.interpwave = interp_wave
        #return(interp_flux_list, interp_wave)

    def cross_correlate(self):#flux_list):
        ### correlate an array of fluxes where last flux is the one correlating to
        correlate_list = []
        index_list = []
        max_val_ind = []
        for index in range(len(self.interpflux)):
            correlate = signal.correlate(self.interpflux[index], self.interpcoflux)
            correlate /= np.nanmax(correlate)
            ind_where = np.argmax(correlate)
            max_val_ind.append(ind_where)
            correlate_list.append(correlate)
            inds = np.arange(len(correlate))
            index_list.append(inds)

        self.maxindex = max_val_ind     
        self.correlate = correlate_list
        self.corrindex = index_list
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
            delta = (len(self.correlate[index])/2) + 0.5 - corr_max_ind
            print(corr_max_ind, delta, len(self.correlate[0]), ind_pts)
            logshift = delta * 1/self.oversample
            act_list.append(corr_max_ind)
            delta_list.append(delta)
            logshift_list.append(logshift)
        self.pixelshift = act_list
        self.pixeldelta = delta_list
        self.logshift = logshift_list
        #return(act_list, delta_list)

    def rvel(self):#delta_list, wave_list):
        c = 299792458   #m/s
        velocity_list = []
        for index in range(len(self.logshift)):
            delta = self.logshift[index]
            velocity = (delta)*c
            velocity_list.append(velocity)
        self.velocity = velocity_list


oversample_factor = 100         #smaller grid used to interpolate finer values
exposures = ['0031', '0032', '0033']
blue = target('/home/ryan/Downloads/200708/data/', 1, '08jul', exposures, 225, 10**7)
green = target('/home/ryan/Downloads/200708/data/', 2, '08jul', exposures, 225, 10**7)
red = target('/home/ryan/Downloads/200708/data/', 3, '08jul', exposures, 225, 10**7)
ir = target('/home/ryan/Downloads/200708/data/', 4, '08jul', exposures, 225, 10**7)
#print(blue.fiber, blue.exp, blue.flux[0], blue.wave, blue.error, 'sick')

exposure_0904 = ['0013','0014']
exposure_0906 = ['0040','0041','0042','0043']
exposure_0907 = ['0051','0052','0053']

blue0904 = target('/home/ryan/Downloads/200904/data', 1, '04sep', exposure_0904, 381, 10**7)
green0904 = target('/home/ryan/Downloads/200904/data', 2, '04sep', exposure_0904, 381, 10**7)
red0904 = target('/home/ryan/Downloads/200904/data', 3, '04sep', exposure_0904, 381, 10**7)
ir0904 = target('/home/ryan/Downloads/200904/data', 4, '04sep', exposure_0904, 381, 10**7)

blue0906 = target('/home/ryan/Downloads/200906/data', 1, '06sep', exposure_0906, 381, 10**7)
green0906 = target('/home/ryan/Downloads/200906/data', 2, '06sep', exposure_0906, 381, 10**7)
red0906 = target('/home/ryan/Downloads/200906/data', 3, '06sep', exposure_0906, 381, 10**7)
ir0906 = target('/home/ryan/Downloads/200906/data', 4, '06sep', exposure_0906, 381, 10**7)

blue0907 = target('/home/ryan/Downloads/200907/data', 1, '07sep', exposure_0907, 381, 10**7)
green0907 = target('/home/ryan/Downloads/200907/data', 2, '07sep', exposure_0907, 381, 10**7)
red0907 = target('/home/ryan/Downloads/200907/data', 3, '07sep', exposure_0907, 381, 10**7)
ir0907 = target('/home/ryan/Downloads/200907/data', 4, '07sep', exposure_0907, 381, 10**7)

bluevelo = blue0904.velocity + blue0906.velocity + blue0907.velocity
greenvelo = green0904.velocity + green0906.velocity + green0907.velocity
redvelo = red0904.velocity + red0906.velocity + red0907.velocity
irvelo = ir0904.velocity + ir0906.velocity + ir0907.velocity

blue_shift = blue0904.logshift + blue0906.logshift + blue0907.logshift
print(blue_shift, bluevelo, greenvelo, redvelo, irvelo)


#print(len(blue.corrindex[0]), len(blue.interpwave), len(blue.correlate[0]), len(blue.normflux[0]))
#print(blue.logshift, blue.velocity, green.velocity, red.velocity, ir.velocity)

x_vals = [1,2,6,7,8,9,10,11,12]
plt.plot(x_vals[0:6],bluevelo[0:6], c='C0')
plt.plot(x_vals,greenvelo, c='C2')
plt.plot(x_vals,redvelo, c='C3')
plt.plot(x_vals,irvelo, c='C7')
plt.ylabel('m / s')
plt.xlabel('Exposure')
plt.show()



'''


### plotting cross-correlation results
fig = plt.figure(figsize=(16,16))
gs = fig.add_gridspec(4,3)

ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])
ax3 = fig.add_subplot(gs[0,2])
ax4 = fig.add_subplot(gs[1,0])
ax5 = fig.add_subplot(gs[1,1])
ax6 = fig.add_subplot(gs[1,2])
ax7 = fig.add_subplot(gs[2,0])
ax8 = fig.add_subplot(gs[2,1])
ax9 = fig.add_subplot(gs[2,2])
ax10 = fig.add_subplot(gs[3,0])
ax11 = fig.add_subplot(gs[3,1])
ax12 = fig.add_subplot(gs[3,2])

im1 = ax1.plot(blue0907.corrindex[0],blue0907.correlate[0], c='b')
im2 = ax2.plot(blue0907.corrindex[1],blue0907.correlate[1], c='b')
im3 = ax3.plot(blue0907.corrindex[2],blue0907.correlate[2], c='b')
im4 = ax4.plot(green0907.corrindex[0],green0907.correlate[0], c='g')
im5 = ax5.plot(green0907.corrindex[1],green0907.correlate[1], c='g')
im6 = ax6.plot(green0907.corrindex[2],green0907.correlate[2], c='g')
im7 = ax7.plot(red0907.corrindex[0],red0907.correlate[0], c='r')
im8 = ax8.plot(red0907.corrindex[1],red0907.correlate[1], c='r')
im9 = ax9.plot(red0907.corrindex[2],red0907.correlate[2], c='r')
im10 = ax10.plot(ir0907.corrindex[0],ir0907.correlate[0], c='k')
im11 = ax11.plot(ir0907.corrindex[1],ir0907.correlate[1], c='k')
im12 = ax12.plot(ir0907.corrindex[2],ir0907.correlate[2], c='k')

ax1.set_xlim(blue0907.maxindex[0]-1.2, blue0907.maxindex[0]+1.2)
ax2.set_xlim(blue0907.maxindex[1]-1.2, blue0907.maxindex[1]+1.2)
ax3.set_xlim(blue0907.maxindex[2]-1.2, blue0907.maxindex[2]+1.2)
ax4.set_xlim(green0907.maxindex[0]-1.2, green0907.maxindex[0]+1.2)
ax5.set_xlim(green0907.maxindex[1]-1.2, green0907.maxindex[1]+1.2)
ax6.set_xlim(green0907.maxindex[2]-1.2, green0907.maxindex[2]+1.2)
ax7.set_xlim(red0907.maxindex[0]-1.2, red0907.maxindex[0]+1.2)
ax8.set_xlim(red0907.maxindex[1]-1.2, red0907.maxindex[1]+1.2)
ax9.set_xlim(red0907.maxindex[2]-1.2, red0907.maxindex[2]+1.2)
ax10.set_xlim(ir0907.maxindex[0]-1.2, ir0907.maxindex[0]+1.2)
ax11.set_xlim(ir0907.maxindex[1]-1.2, ir0907.maxindex[1]+1.2)
ax12.set_xlim(ir0907.maxindex[2]-1.2, ir0907.maxindex[2]+1.2)

ax1.set_ylim(.99999, 1.00002)
ax2.set_ylim(.9999, 1.00002)
ax3.set_ylim(.9999, 1.00002)
ax4.set_ylim(.9999, 1.00002)
ax5.set_ylim(.9999, 1.00002)
ax6.set_ylim(.9999, 1.00002)
ax7.set_ylim(.9999, 1.00002)
ax8.set_ylim(.9999, 1.00002)
ax9.set_ylim(.9999, 1.00002)
ax10.set_ylim(.9999, 1.00002)
ax11.set_ylim(.9999, 1.00002)
ax12.set_ylim(.9999, 1.00002)

ax1.axvline(blue0907.pixelshift[0],0,1)
ax2.axvline(blue0907.pixelshift[1],0,1)
ax3.axvline(blue0907.pixelshift[2],0,1)
ax4.axvline(green0907.pixelshift[0],0,1)
ax5.axvline(green0907.pixelshift[1],0,1)
ax6.axvline(green0907.pixelshift[2],0,1)
ax7.axvline(red0907.pixelshift[0],0,1)
ax8.axvline(red0907.pixelshift[1],0,1)
ax9.axvline(red0907.pixelshift[2],0,1)
ax10.axvline(ir0907.pixelshift[0],0,1)
ax11.axvline(ir0907.pixelshift[1],0,1)
ax12.axvline(ir0907.pixelshift[2],0,1)

plt.show()
'''

### plotting the interpolated fluxes
fig = plt.figure(figsize=(16,16))
gs = fig.add_gridspec(4,4)
fig.suptitle('Interpolated Spectra')

ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])
ax3 = fig.add_subplot(gs[0,2])
ax13 = fig.add_subplot(gs[0,3])
ax4 = fig.add_subplot(gs[1,0])
ax5 = fig.add_subplot(gs[1,1])
ax6 = fig.add_subplot(gs[1,2])
ax14 = fig.add_subplot(gs[1,3])
ax7 = fig.add_subplot(gs[2,0])
ax8 = fig.add_subplot(gs[2,1])
ax9 = fig.add_subplot(gs[2,2])
ax15 = fig.add_subplot(gs[2,3])
ax10 = fig.add_subplot(gs[3,0])
ax11 = fig.add_subplot(gs[3,1])
ax12 = fig.add_subplot(gs[3,2])
ax16 = fig.add_subplot(gs[3,3])

im1 = ax1.plot(blue0907.interpwave,blue0907.interpflux[0], c='b')
im2 = ax2.plot(blue0907.interpwave,blue0907.interpflux[1], c='b')
im3 = ax3.plot(blue0907.interpwave,blue0907.interpflux[2], c='b')
im13 = ax13.plot(blue0907.interpcowave,blue0907.interpcoflux, c='b')
im4 = ax4.plot(green0907.interpwave,green0907.interpflux[0], c='g')
im5 = ax5.plot(green0907.interpwave,green0907.interpflux[1], c='g')
im6 = ax6.plot(green0907.interpwave,green0907.interpflux[2], c='g')
im14 = ax14.plot(green0907.interpcowave,green0907.interpcoflux, c='g')
im7 = ax7.plot(red0907.interpwave,red0907.interpflux[0], c='r')
im8 = ax8.plot(red0907.interpwave,red0907.interpflux[1], c='r')
im9 = ax9.plot(red0907.interpwave,red0907.interpflux[2], c='r')
im15 = ax15.plot(red0907.interpcowave,red0907.interpcoflux, c='r')
im10 = ax10.plot(ir0907.interpwave,ir0907.interpflux[0], c='k')
im11 = ax11.plot(ir0907.interpwave,ir0907.interpflux[1], c='k')
im12 = ax12.plot(ir0907.interpwave,ir0907.interpflux[2], c='k')
im16 = ax16.plot(ir0907.interpcowave,ir0907.interpcoflux, c='k')

ax1.axhline(0, xmin=-150, xmax=150, zorder=4, c='k')
ax2.axhline(0, xmin=-150, xmax=150, zorder=4, c='k')
ax3.axhline(0, xmin=-150, xmax=150, zorder=4, c='k')
ax13.axhline(0, xmin=-150, xmax=150, zorder=4, c='k')
ax4.axhline(0, xmin=-150, xmax=150, zorder=4, c='k')
ax5.axhline(0, xmin=-150, xmax=150, zorder=4, c='k')
ax6.axhline(0, xmin=-150, xmax=150, zorder=4, c='k')
ax14.axhline(0, xmin=-150, xmax=150, zorder=4, c='k')
ax7.axhline(0, xmin=-150, xmax=150, zorder=4, c='k')
ax8.axhline(0, xmin=-150, xmax=150, zorder=4, c='k')
ax9.axhline(0, xmin=-150, xmax=150, zorder=4, c='k')
ax15.axhline(0, xmin=-150, xmax=150, zorder=4, c='k')
ax10.axhline(0, xmin=-150, xmax=150, zorder=4, c='y')
ax11.axhline(0, xmin=-150, xmax=150, zorder=4, c='y')
ax12.axhline(0, xmin=-150, xmax=150, zorder=4, c='y')
ax16.axhline(0, xmin=-150, xmax=150, zorder=4, c='y')
ax1.set_xlim(np.nanmin(blue0907.interpwave)-0.001, np.nanmax(blue0907.interpwave)+0.001)
ax2.set_xlim(np.nanmin(blue0907.interpwave)-0.001, np.nanmax(blue0907.interpwave)+0.001)
ax3.set_xlim(np.nanmin(blue0907.interpwave)-0.001, np.nanmax(blue0907.interpwave)+0.001)
ax13.set_xlim(np.nanmin(blue0907.interpcowave)-0.001, np.nanmax(blue0907.interpcowave)+0.001)
ax4.set_xlim(np.nanmin(green0907.interpwave)-0.001, np.nanmax(green0907.interpwave)+0.001)
ax5.set_xlim(np.nanmin(green0907.interpwave)-0.001, np.nanmax(green0907.interpwave)+0.001)
ax6.set_xlim(np.nanmin(green0907.interpwave)-0.001, np.nanmax(green0907.interpwave)+0.001)
ax14.set_xlim(np.nanmin(green0907.interpcowave)-0.001, np.nanmax(green0907.interpcowave)+0.001)
ax7.set_xlim(np.nanmin(red0907.interpwave)-0.001, np.nanmax(red0907.interpwave)+0.001)
ax8.set_xlim(np.nanmin(red0907.interpwave)-0.001, np.nanmax(red0907.interpwave)+0.001)
ax9.set_xlim(np.nanmin(red0907.interpwave)-0.001, np.nanmax(red0907.interpwave)+0.001)
ax15.set_xlim(np.nanmin(red0907.interpcowave)-0.001, np.nanmax(red0907.interpcowave)+0.001)
ax10.set_xlim(np.nanmin(ir0907.interpwave)-0.001, np.nanmax(ir0907.interpwave)+0.001)
ax11.set_xlim(np.nanmin(ir0907.interpwave)-0.001, np.nanmax(ir0907.interpwave)+0.001)
ax12.set_xlim(np.nanmin(ir0907.interpwave)-0.001, np.nanmax(ir0907.interpwave)+0.001)
ax16.set_xlim(np.nanmin(ir0907.interpcowave)-0.001, np.nanmax(ir0907.interpcowave)+0.001)

ax1.set_title('Exp 01')
ax2.set_title('Exp 02')
ax3.set_title('Exp 03')
ax13.set_title('Coadded')

plt.show()

'''
### plotting the interpolated fluxes
fig = plt.figure(figsize=(16,16))
gs = fig.add_gridspec(4,3)
fig.suptitle('Interpolated Spectra')

ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])
ax3 = fig.add_subplot(gs[0,2])
ax4 = fig.add_subplot(gs[1,0])
ax5 = fig.add_subplot(gs[1,1])
ax6 = fig.add_subplot(gs[1,2])
ax7 = fig.add_subplot(gs[2,0])
ax8 = fig.add_subplot(gs[2,1])
ax9 = fig.add_subplot(gs[2,2])
ax10 = fig.add_subplot(gs[3,0])
ax11 = fig.add_subplot(gs[3,1])
ax12 = fig.add_subplot(gs[3,2])

im1 = ax1.plot(bluew[0],bluef[0], c='b')
im2 = ax2.plot(bluew[1],bluef[1], c='b')
im3 = ax3.plot(bluew[2],bluef[2], c='b')
im4 = ax4.plot(greenw[0],greenf[0], c='g')
im5 = ax5.plot(greenw[1],greenf[1], c='g')
im6 = ax6.plot(greenw[2],greenf[2], c='g')
im7 = ax7.plot(redw[0],redf[0], c='r')
im8 = ax8.plot(redw[1],redf[1], c='r')
im9 = ax9.plot(redw[2],redf[2], c='r')
im10 = ax10.plot(irw[0],irf[0], c='k')
im11 = ax11.plot(irw[1],irf[1], c='k')
im12 = ax12.plot(irw[2],irf[2], c='k')

plt.show()
'''



