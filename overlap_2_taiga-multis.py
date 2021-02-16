import csv
import eleanor
import numpy as np
import lightkurve as lk
from astropy.io import fits
from astropy import units as u
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from transitleastsquares import transitleastsquares



###Opening the Simultaneously observed GALAH-TESS targets from before Feb 2019
target_id = []
target_ra = []
target_dec = []
target_sector = []
target_start_mjd = []
target_end_mjd = []

with open('/home/rmorris/documents/simult_obs_targets.txt', 'r') as file:
	reader = csv.reader(file,delimiter='\t')
	for row in reader:
		target_id.append(row[0])
		target_ra.append(float(row[1]))
		target_dec.append(float(row[2]))
		target_sector.append(int(row[3]))
		target_start_mjd.append(float(row[4]))
		target_end_mjd.append(float(row[5]))

data_source = []
pre_std = np.zeros(len(target_id))
post_std = np.zeros(len(target_id))
std_ratio = np.zeros(len(target_id))
std_bins = np.zeros(len(target_id))
for index,ids in enumerate(target_id):
	if target_sector[index] <= 31:
		if target_id.count(ids) != 1:
			print(ids)
			coords = SkyCoord(ra=float(target_ra[index]), dec=float(target_dec[index]), unit=(u.deg, u.deg))
			try:
				###first going to try 2-min data
				lc = lk.search_lightcurvefile(coords, mission='TESS', sector=target_sector[index]).download().PDCSAP_FLUX.remove_nans()
				time = lc.time
				time_mjd = np.zeros(len(time))
				for ind1, times in enumerate(time):
					time_mjd[ind1] = times+57000.5
				lc = lk.LightCurve(time=time_mjd, flux=lc.flux)
				data_source.append('2-min')
			except Exception as e:
				print(e)
				try:
					###then try to see if the fits file already exists
					data = fits.open('/data/wallaby/rmorris/GALAH/all_target_lc/{}_s{}_eleanor_data.fits'.format(ids,target_sector[index]))
					print('Fits file opened')
					time = np.zeros(len(data[1].data))
					flux = np.zeros(len(data[1].data))
					background = np.zeros(len(data[1].data))
					for ind2, x in enumerate(data[1].data):
						if x[7] == 0:
							time[ind2] = x[0]
							flux[ind2] = x[5]
							background[ind2] = x[12]
					print('Data input into arrays')
					print(len(time), len(flux))
					med_corr_flux = np.median(flux)
					for fluxval in flux:
						fluxval = fluxval/med_corr_flux

					time_arr_med = np.ma.masked_where(time == 0, time)
					flux_arr_med = np.ma.masked_where(time == 0, flux)
					back_arr_med = np.ma.masked_where(time == 0, background)
					t_arrc = np.ma.MaskedArray.compressed(time_arr_med)
					f_arrc = np.ma.MaskedArray.compressed(flux_arr_med)
					b_arrc = np.ma.MaskedArray.compressed(back_arr_med)
					print('Bad quality points are masked')

					time_mjd = np.zeros(len(t_arrc))
					for ind3, times in enumerate(t_arrc):
						time_mjd[ind3] = times+57000.5

					print(len(time_mjd), len(f_arrc))
					time_bkg = np.zeros(len(time_mjd))
					flux_bkg = np.zeros(len(time_mjd))
					for num, back in enumerate(b_arrc):
						bkg_med = np.median(b_arrc)
						if back < bkg_med*3 and back > bkg_med/3:
							time_bkg[num] = time_mjd[num]
							flux_bkg[num] = f_arrc[num]
					print('Bad background points masked')
					time_arr_fin = np.ma.masked_where(time_bkg == 0, time_bkg)
					flux_arr_fin = np.ma.masked_where(time_bkg == 0, flux_bkg)
					t_arrcf = np.ma.MaskedArray.compressed(time_arr_fin)
					f_arrcf = np.ma.MaskedArray.compressed(flux_arr_fin)
					print('Final lc created') 
					print(len(t_arrcf), len(f_arrcf)) 
					lc = lk.LightCurve(time=t_arrcf, flux=f_arrcf)
					data_source.append('30-min')
					print("Sector {} found locally with eleanor".format(target_sector[index]))
				except Exception as e:
					print(e)
					print('File does not already exist, downloading postcards instead')
					try:
						###then trying 30-min data
						coords = (float(target_ra[index]), float(target_dec[index]))
						source = eleanor.Source(coords=coords, sector=target_sector[index])
						data = eleanor.TargetData(source, aperture_mode='small')
						data.save(output_fn='{}_s{}_eleanor_data.fits'.format(ids,target_sector[index]),directory='/data/wallaby/rmorris/GALAH/all_target_lc/')
						data_source.append('30-min')

						###only using quality points
						time = []
						flux = []
						background = []


						q = data.quality == 0
						time.append(data.time[q])
						flux.append(data.corr_flux[q]/np.median(data.corr_flux[q]))
						background.append(data.flux_bkg[q])

						#not sure why I do this part, seems sus, might need to delete

						time = time[0]
						flux = flux[0]
						background = background[0]

						time_mjd = np.zeros(len(time))
						for ind4, times in enumerate(time):
							time_mjd[ind4] = times+57000.5

						print("Eleanor data downloaded correctly")

						###getting rid of any points that have a background flux > 3x the median background flux
						time_bkg = np.zeros(len(time_mjd))
						flux_bkg = np.zeros(len(time_mjd))
						for num, back in enumerate(background):
							bkg_med = np.median(background)
							if back < bkg_med*3:
								time_bkg[num] = time_mjd[num]
								flux_bkg[num] = flux[num]

						print("Bad background points subtracted")

						###create lightkurve object
						time_arr_fin = np.ma.masked_where(time_bkg == 0, time_bkg)
						flux_arr_fin = np.ma.masked_where(time_bkg == 0, flux_bkg)
						t_arrcf = np.ma.MaskedArray.compressed(time_arr_fin)
						f_arrcf = np.ma.MaskedArray.compressed(flux_arr_fin)

						lc = lk.LightCurve(time=t_arrcf, flux=f_arrcf)

						print('LightKurve object created')

					except KeyboardInterrupt:
						raise
					except Exception as e: 
						print(e)
						print('Target {} failed during data acquisition'.format(ids))

			try:
				###check for rotation periods by difference in std when flatten vs not flattened
				pre_std[index] = np.std(lc.flux)
				pre = np.std(lc.flux)
				print('test1')
				lc = lc.flatten()
				print('test1.1')
				post_std[index] = np.std(lc.flux)
				print('test1.2')
				post = np.std(lc.flux)
				print('test2')
				std_ratio[index] = pre/post
				bin_ratio = [1,2,5,10]
				if pre/post < 2.0:
					std_bins[index] = bin_ratio[0]
					binner = bin_ratio[0]
				elif pre/post < 5.0:
					std_bins[index] = bin_ratio[1]
					binner = bin_ratio[1]
				elif pre/post < 10.0:
					std_bins[index] = bin_ratio[2]
					binner = bin_ratio[2]
				else:
					std_bins[index] = bin_ratio[3]
					binner = bin_ratio[3]

				print('Ratios established')
				print(pre, post, pre/post)
			except KeyboardInterrupt:
				raise
			except Exception as e: 
				print(e)
				print('Target {} failed during standard deviation calcs'.format(ids))

			try:
				fig, (ax1, ax2) = plt.subplots(1,2)
				fig.set_figheight(16)
				fig.set_figwidth(32)
				fig.suptitle(ids + " Lightcurve")
				ax1.scatter(lc.time, lc.flux, c='k', marker='.')
				ax2.scatter(lc.time, lc.flux, c='k', marker='.')
				print('first part of plot check')
				multis_indices = [mi for mi, mx in enumerate(target_id) if mx == ids]
				sectrs = []
				imy = 0
				while imy < len(multis_indices):
					sectrs.append(target_sector[multis_indices[imy]])
					imy += 1
				if sectrs.count(target_sector[index]) > 1:
					for mult in multis_indices:
						if target_sector[index] == target_sector[mult]:
							ax2.axvspan(target_start_mjd[mult],target_end_mjd[mult], color='b', alpha=0.5)
							ax1.axvspan(target_start_mjd[mult],target_end_mjd[mult], color='b', alpha=0.5)
				###rest of code regardless of which source its from
				###plotting lightcurves with galah markers
				#plt.figure(num=str(ids) + " Lightcurve", figsize=(16,16))
				print('mathy part check')
				ax2.set_xlim(left=target_start_mjd[index]-1,right=target_end_mjd[index]+1)
				ax1.set_xlabel('Time')
				ax1.set_ylabel('Normalized Flux')
				print('second part plot check')
				plt.savefig('/data/wallaby/rmorris/GALAH/Simult_Plots_Ratios/Multis/{}_lc_sector_{}_{}_{}.png'.format(ids, target_sector[index], data_source[index], index))

				print("Lightcurve successfully plotted, and saved in Simult_Plots_Ratios")

				with open('/data/wallaby/rmorris/GALAH/Simult_Plots_Ratios/Bin_{}/{}_sector_{}.txt'.format(binner,ids,target_sector[index]), 'w') as file:
					for number, values in enumerate(lc.time):
						file.write(str(values)+'\t'+str(lc.flux[number])+'\n')

				print("Lightcurve data saved to file")


			###if any part of this process fails
			except KeyboardInterrupt:
				raise
			except Exception as e: 
				print(e)
				print('Target {} failed during plotting'.format(ids))

			finally:
			###not cause any open plot problems
				plt.close('all')
		else:
			print("Target {} is not viewed by GALAH multiple times".format(ids))
			data_source.append('None')
	else:
		print("Target {} is from a sector after 31".format(ids))
		data_source.append('None')
