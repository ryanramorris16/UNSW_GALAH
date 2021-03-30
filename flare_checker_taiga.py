import os
import csv
import eleanor
import imageio
import numpy as np
import lightkurve as lk
from astropy.io import fits
from astropy import units as u
from IPython.display import HTML
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astropy.coordinates import SkyCoord
from transitleastsquares import transitleastsquares
import scipy.ndimage

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

### Set the id of the potential flares to make videos for each one
good_shot_flares = ['gaiadr2_4699289920904010240','gaiadr2_4699457493348110976','gaiadr2_4699459348773982976','gaiadr2_4708499533377604096','gaiadr2_4711447873807672192','gaiadr2_4711619496404895616','gaiadr2_4711642143768625280','gaiadr2_4711648053643629696','gaiadr2_4711650454529161088','gaiadr2_4711729108266385024','gaiadr2_4712630054965882624','gaiadr2_4714036158539045760','gaiadr2_4739406633436466816','gaiadr2_4900721791142464384','gaiadr2_6581509791634606080','gaiadr2_6581817242573714560','gaiadr2_6666186480624454912','gaiadr2_6666252279523495680','gaiadr2_6695525814540719872','gaiadr2_6774017750624069760','gaiadr2_6815700270791439616','gaiadr2_6816077265840584064','gaiadr2_6828181136515742336','galahic_225809','galahic_244867','gaiadr2_4694949014637901696','gaiadr2_4900721791142464384','gaiadr2_4708486923353622272','gaiadr2_4711770992787233408','gaiadr2_4712912835612572416','gaiadr2_4714031932291210112','gaiadr2_4714966173577866240','gaiadr2_4739351314257766656','gaiadr2_4739415085932156800','gaiadr2_4900616852206496896','gaiadr2_4900638434417245824','gaiadr2_4900657465416344192','gaiadr2_4900759689933956480','gaiadr2_4901816732925172864','gaiadr2_6401393466129296000','gaiadr2_6401751666402447360','gaiadr2_6451676847286237696','gaiadr2_6452143658691602432','gaiadr2_6452216363897984128','gaiadr2_6469631494210833152','gaiadr2_6470046937807488640','gaiadr2_6470350334299245056','gaiadr2_6472876011289820928','gaiadr2_6473145155415934976','gaiadr2_6473168382599035776','gaiadr2_6497271704704136064','gaiadr2_6521430449068817536','gaiadr2_6576277628114104960','gaiadr2_6581789136307685632','gaiadr2_6667025201837760640','gaiadr2_6695144524523613312','gaiadr2_6695974312205359232','gaiadr2_6774003727553563904','gaiadr2_6814806982017417728','gaiadr2_6814909137815395456','gaiadr2_6815009360876917120','gaiadr2_6815700025977976960','gaiadr2_6815746553358965632','gaiadr2_6815902988952395392','gaiadr2_6816048609818932096','gaiadr2_6816086744833477760','gaiadr2_6816131348068815744','gaiadr2_6816131378132829952','gaiadr2_6816174641339254272','gaiadr2_6817613592822257536','gaiadr2_6817626713947107200','gaiadr2_6817639259546828800','gaiadr2_6828172890178513152','gaiadr2_6828181136515742336','galahic_235582','galahic_5329599']#['gaiadr2_6816048609818932096','gaiadr2_6816131348068815744','gaiadr2_6816131378132829952','gaiadr2_6817613592822257536','gaiadr2_4708486923353622272','gaiadr2_4711880050596718720','gaiadr2_4712859887255955968','gaiadr2_4712912835612572416','gaiadr2_4900655236329321344','gaiadr2_4900759689933956480','gaiadr2_4901673315377052800','gaiadr2_4972046388882790656','gaiadr2_6401393466129296000','gaiadr2_6401751666402447360','gaiadr2_6451676847286237696','gaiadr2_6452216363897984128','gaiadr2_6470350334299245056','gaiadr2_6473168382599035776','gaiadr2_6521313179280807808','gaiadr2_6576063807462249728','gaiadr2_6576277628114104960','gaiadr2_6578645357685227776','gaiadr2_6578652985547127424','gaiadr2_6581618265329053184','gaiadr2_6581889844700683648','gaiadr2_6666230358010367488','gaiadr2_6666386346925881856','gaiadr2_6666405352157548928','gaiadr2_6695127275935016704','gaiadr2_6695149605468630400','gaiadr2_6695453452931502080','gaiadr2_6695545330872170752','gaiadr2_6695557696081160704','gaiadr2_6695869549362898432','gaiadr2_6814219675305265408','gaiadr2_6814449507595558528','gaiadr2_6814455657988747648','gaiadr2_6815688932077662464','gaiadr2_6815744624918377728','gaiadr2_6815902988952395392','gaiadr2_6816050877561660800','gaiadr2_6816068778985180288','gaiadr2_6816086744833477760','gaiadr2_6817626713947107200','gaiadr2_6817639736287944320','gaiadr2_6817712239631141632','gaiadr2_6828172890178513152','gaiadr2_6828181136515742336','gaiadr2_6828220444056545536']
good_shot_ra = []
good_shot_dec = []
good_shot_sector = []
good_shot_start = []
good_shot_end = []

for flare in good_shot_flares:
	try:
		index = target_id.index(flare)
		good_shot_ra.append(target_ra[index])
		good_shot_dec.append(target_dec[index])
		good_shot_sector.append(target_sector[index])
		good_shot_start.append(target_start_mjd[index])
		good_shot_end.append(target_end_mjd[index])
	except Exception as e:
		print(e)
      
###  
for index, target in enumerate(good_shot_flares):
	try:
		print(target)
		coords = SkyCoord(ra=float(good_shot_ra[index]), dec=float(good_shot_dec[index]), unit=(u.deg,u.deg))
		source = eleanor.Source(coords=coords, sector=good_shot_sector[index])
		data = eleanor.TargetData(source, aperture_mode='small')

		time = []
		flux = []
		background = []
		tpf_qual = []

		q = data.quality == 0
		time.append(data.time[q])
		flux.append(data.corr_flux[q]/np.median(data.corr_flux[q]))
		background.append(data.flux_bkg[q])
		tpf_qual.append(data.tpf[q])

		time = time[0]
		flux = flux[0]
		background = background[0]
		tpf_qual = tpf_qual[0]
		#print(len(background),len(tpf_qual)) 

		time_mjd = []
		for times in time:
			time_mjd.append(times+57000.5)

		print("Eleanor data downloaded correctly")

		###getting rid of any points that have a background flux > 3x the median background flux
		time_bkg = []
		flux_bkg = []
		for num, back in enumerate(background):
			bkg_med = np.median(background)
			if back < bkg_med*3:
				time_bkg.append(time_mjd[num])
				flux_bkg.append(flux[num])

		print("Bad background points subtracted")

		###create lightkurve object
		lc = lk.LightCurve(time = time_bkg, flux = flux_bkg).flatten()

		print('LightKurve object created')

		mid_time_ind = 0
		for nums, timestamps in enumerate(lc.time):
			if mid_time_ind == 0 and timestamps != lc.time[-1]:
				if timestamps-0.01 <= good_shot_start[index] <= lc.time[nums+1]+0.01:
					mid_time_ind = np.where(timestamps == lc.time)
					mid_time_ind = mid_time_ind[0][0]

		images_num = np.arange(mid_time_ind-144, mid_time_ind+144, 1)
		med_images = []
		sub_images = []

		#print('Middle time created')

		for images in images_num:
			low = images-144
			high = images+144
			if low <= 0:
				low = 0
			if high >= len(tpf_qual):
				high = len(tpf_qual)
			#print('test1')   
			med_range = tpf_qual[low:high]		  ###+/- 1 day of ten minute exposures
			rolling_med = np.median(med_range, axis=0)
			#print('test2')
			#print(len(tpf_qual)) 
			pre_sub = tpf_qual[images]
			#print('test3') 
			sub = pre_sub-rolling_med
			med_images.append(rolling_med)
			sub_images.append(sub)

		print('Median images created')

		std_for_plot = np.std(lc.time[mid_time_ind-144:mid_time_ind+144])

		fig = plt.figure()
		gs = fig.add_gridspec(2,2)
		ax1 = fig.add_subplot(gs[0,0])
		ax2 = fig.add_subplot(gs[0,1])
		ax3 = fig.add_subplot(gs[1,:])
		im1 = ax1.imshow(med_images[0])#, cmap='Greys_r')
		im2 = ax2.imshow(sub_images[0], vmin=-(2*std_for_plot), vmax=(2*std_for_plot))
		im3 = ax3.scatter(lc.time, lc.flux, c='k', marker='.')
		im3 = ax3.axvspan(good_shot_start[index], good_shot_end[index],color='b', alpha=0.5)
		im3 = ax3.set_xlim(left=good_shot_start[index]-1, right=good_shot_end[index]+1)
		#im3 = ax3.axvline(x=lc.time[images_num[0]], ymin=-0.1, ymax=1.1, c='r')
		fig.tight_layout()
		fig.subplots_adjust(top=0.9)
		def animate(i):
			im1.set_data(med_images[i])
			ax1.set_title('{}'.format(images_num[i]))
			ax2.set_title('MJD {}'.format(lc.time[images_num[i]]))
			im2.set_data(sub_images[i])
			im3 = ax3.clear()
			im3 = ax3.scatter(lc.time, lc.flux, c='k', marker='.')
			im3 = ax3.axvspan(good_shot_start[index], good_shot_end[index],color='b', alpha=0.5)
			im3 = ax3.set_xlim(left=good_shot_start[index]-1, right=good_shot_end[index]+1)
			im3 = ax3.axvline(x=lc.time[images_num[i]], ymin=-0.1, ymax=1.1, c='r')
			im3 = ax3.scatter(lc.time[i], lc.flux[i], c='orange', marker='*')
			return im1, im2, im3,
		anim = animation.FuncAnimation(fig, animate, frames=np.arange(0, len(images_num), 1), blit=True)
		anim.save('/data/wallaby/rmorris/GALAH/Simult_Gifs/Jan28/{}.mp4'.format(target))
		plt.close('all')
	except Exception as e:
		print(e)
		print('Oh no! Target failed. Anyway...')