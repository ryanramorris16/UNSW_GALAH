import numpy as np
import matplotlib.pyplot as plt
from astropy.nddata import Cutout2D
from photutils import CircularAperture, RectangularAperture, aperture_photometry
from photutils import MMMBackground
from lightkurve import lightcurve
from lightkurve.correctors import SFFCorrector
from scipy.optimize import minimize
from astropy import time, coordinates as coord, units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table, Column
from astropy.wcs import WCS
from astropy.stats import SigmaClip, sigma_clip
from time import strftime
from astropy.io import fits
from scipy.stats import mode
from scipy.signal import savgol_filter
from urllib.request import urlopen
import os, sys, copy
import os.path
import warnings
import pickle

def get_cbvs(sector,camera,chip,time):
        """ Obtains the cotrending basis vectors (CBVs) as convolved down from the short-cadence targets.
        Parameters
        ----------
        """
        #print('/data/wallaby/rmorris/.eleanor/metadata/s{0:04d}/cbv_components_s{0:04d}_{1:04d}_{2:04d}.txt'.format(int(sector),int(sector),int(camera),int(chip)))
        #matrix_file = np.loadtxt('/data/wallaby/rmorris/.eleanor/metadata/s{0:04d}/cbv_components_s{0:04d}_{1:04d}_{2:04d}.txt'.format(int(sector),int(camera),int(chip)))
        try:
            #print('Trying CBVs')
            matrix_file = np.loadtxt('/home/z5318114/.eleanor/metadata/s{0:04d}/cbv_components_s{0:04d}_{1:04d}_{2:04d}.txt'.format(int(sector),int(camera),int(chip)))
            cbvs = np.asarray(matrix_file)
            cbvs = np.reshape(cbvs, (len(time), 16))
            #print('CBVs Acquired')

        except:
            cbvs = np.zeros((len(time), 16))
            #print('No CBVs')
        return cbvs

def norm(l, q):
    l = l[q]
    l /= np.nanmedian(l)
    l -= 1
    return l

def fhat(xhat, data):
    return np.dot(data, xhat)

def xhat(mat, lc):
    return np.linalg.lstsq(mat, lc, rcond=None)[0]

def corrected_flux(self, time=None, flux=None, quality=None, centroid_xs=None, centroid_ys=None, skip=0.25, modes=3, pca=False, bkg=None, regressors=None, sector=None, camera=None, chip=None):
        """
        Corrects for jitter in the light curve by quadratically regressing with centroid position.
        Parameters
        ----------
        skip: float, optional
            The length of time in days at the start of each orbit to skip in determining optimal model weights.
            Default is 0.62 days.
        modes: int, optional
            If doing a PCA-based correction, the number of cotrending basis vectors to regress against.
            Default is 3.
        pca : bool, optional
            Allows users to decide whether or not to use PCA analysis. Default is False.
        regressors : numpy.ndarray or str
            Extra data to regress against in the correction. Should be shape (len(data.time), N) or `'corner'`.
            If `'corner'` will use the four corner pixels of the TPF as extra information in the regression.
        """
        self.modes = modes
        self.time = time
        self.quality = quality
        self.sector=sector
        self.camera=camera
        self.chip=chip

        tdelt = np.median(np.diff(self.time))
        skip = int(np.ceil(skip/tdelt))


        if type(regressors) == str:
            if regressors == 'corner':
                self.regressors = np.array([self.tpf[:,0,0], self.tpf[:,0,-1], self.tpf[:,-1,-1], self.tpf[:,-1,0]]).T
                regressors = np.array([self.tpf[:,0,0], self.tpf[:,0,-1], self.tpf[:,-1,-1], self.tpf[:,-1,0]]).T
        '''
        if type(self.regressors) == str:
            if self.regressors == 'corner':
                self.regressors = np.array([self.tpf[:,0,0], self.tpf[:,0,-1], self.tpf[:,-1,-1], self.tpf[:,-1,0]]).T
                regressors = np.array([self.tpf[:,0,0], self.tpf[:,0,-1], self.tpf[:,-1,-1], self.tpf[:,-1,0]]).T
        

        if regressors is None:
            regressors = self.regressors
        
        if flux is None:
            flux = self.raw_flux

        if bkg is None:
            bkg = self.flux_bkg
        '''
        flux = np.array(flux)
        med = np.nanmedian(flux)
        quality = copy.deepcopy(self.quality)

        cx = centroid_xs
        cy = centroid_ys
        t  = self.time-self.time[0]

        def calc_corr(mask, cx, cy, skip):
            nonlocal quality, flux, bkg, regressors

            badx = np.where(np.abs(cx - np.nanmedian(cx)) > 3*np.std(cx))[0]
            bady = np.where(np.abs(cy - np.nanmedian(cy)) > 3*np.std(cy))[0]

            # temp_lc = lightcurve.LightCurve(t, flux).flatten()
            tmp_flux = np.copy(flux[np.isfinite(flux)], order="C")
            tmp_flux[:] /= savgol_filter(tmp_flux, 101, 2)
            SC = sigma_clip(tmp_flux, sigma_upper=3.5, sigma_lower=3.5)

            quality[badx] = -999
            quality[bady] = -999

            quality[SC.mask == True] = -999

            qm = quality[mask] == 0

            medval = np.nanmedian(flux[mask][qm])
            norm_l = norm(flux[mask], qm)

            #cx, cy = rotate_centroids(cx[mask], cy[mask])
            cx = cx[mask]
            cy = cy[mask]
            cx -= np.median(cx)
            cy -= np.median(cy)


            bkg_use = bkg[mask]
            bkg_use -= np.min(bkg_use)

            cbvs = get_cbvs(self.sector, self.camera, self.chip, self.time)
            vv = cbvs[mask][:,0:modes]

            if pca == False:
                cm = np.column_stack((t[mask][qm][skip:], np.ones_like(t[mask][qm][skip:])))
                cm_full = np.column_stack((t[mask], np.ones_like(t[mask])))

                if np.std(vv) > 1e-10:
                    cm = np.column_stack((cm, vv[qm][skip:]))
                    cm_full = np.column_stack((cm_full, vv))


                if np.std(bkg) > 1e-10:
                    cm = np.column_stack((cm, bkg_use[qm][skip:]))
                    cm_full = np.column_stack((cm_full, bkg_use))

                if np.std(cx) > 1e-10:
                    cm = np.column_stack((cm, cx[qm][skip:], cy[qm][skip:], cx[qm][skip:]**2, cy[qm][skip:]**2))
                    cm_full = np.column_stack((cm_full, cx, cy, cx**2, cy**2))

                if regressors is not None:
                    cm = np.column_stack((cm, regressors[mask][qm][skip:]))
                    cm_full = np.column_stack((cm_full, regressors[mask]))

            else:
                cm = np.column_stack((vv[qm][skip:], np.ones_like(t[mask][qm][skip:])))
                cm_full = np.column_stack((vv, np.ones_like(t[mask])))


            x = xhat(cm, norm_l[skip:])
            fmod = fhat(x, cm_full)
            lc_pred = (fmod+1)
            return lc_pred*medval

        def find_break(time):
            t   = np.diff(time)
            ind = np.where( t == np.max(t))[0][0]
            return ind + 1

        brk = find_break(self.time)
        f   = np.arange(0, brk, 1); s = np.arange(brk, len(self.time), 1)

        lc_pred = calc_corr(f, cx, cy, skip)
        corr_f = flux[f]-lc_pred + med

        lc_pred = calc_corr(s, cx, cy, skip)
        corr_s = flux[s]-lc_pred + med

        if pca==True:
            self.pca_flux = np.append(corr_f, corr_s)
        else:
            return np.append(corr_f, corr_s)









