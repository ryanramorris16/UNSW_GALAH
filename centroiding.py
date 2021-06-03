from photutils.centroids import centroid_sources, centroid_com
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from astropy import units as u
import lightkurve as lk
import numpy as np
'''
4955894425631686912 160222069
4955894528710881408 160222069
6756793523028737024 425561347
5801881413893258112 311179742
6560686965549688960 140067837
6518496707230129920 139456051
6518497016468468224 139456051
5801881413893258112 311179742
5019047693471033088 32451836
6507063646023467520 161169240
'''
lc = lk.search_lightcurve('TIC 160222069', mission='TESS').download()
#lc[0].plot()
#print(np.nanmedian(lc[0].flux.value), np.nanmedian(lc[2].flux.value), np.nanmedian(lc[4].flux.value))
#plt.scatter(lc[0].time.value,lc[0].flux.value/np.nanmedian(lc[0].flux.value), c='k', marker='.')
#plt.scatter(lc[2].time.value,lc[2].flux.value/np.nanmedian(lc[2].flux.value), c='k', marker='.')
#plt.scatter(lc[4].time.value,lc[4].flux.value/np.nanmedian(lc[4].flux.value), c='k', marker='.')
#plt.show()

#time = np.concatenate((lc[0].time.value,lc[2].time.value,lc[4].time.value))
#flux = np.concatenate((lc[0].flux.value/np.nanmedian(lc[0].flux.value),lc[2].flux.value/np.nanmedian(lc[2].flux.value),lc[4].flux.value/np.nanmedian(lc[4].flux.value)))

#plt.scatter(time, flux)
#plt.show()

new_lc = lk.LightCurve(time=lc.time.value, flux=lc.flux.value)
folded = new_lc.fold(period=5.850943)#, epoch_time = 2459135.870906)
plt.scatter(folded.phase.value, folded.flux.value)
plt.show()

tpfs = lk.search_targetpixelfile('TIC 160222069', mission='TESS').download()

tpf_f_list = []
for i in range(len(tpfs)):
    tpf_f_list.append(tpfs[i].flux)
    print(i)

med_image = np.median(tpf_f_list,axis=0)
med_image = np.squeeze(med_image, axis=0)
med_image = med_image * u.electron/u.s
plt.imshow(med_image)
plt.title('median image')
plt.colorbar()
plt.show()

print(type(med_image))
print(type(tpfs[i].flux))
print(tpfs[i].flux[0])
residuals = []
for i in range(len(tpfs)):
    new = tpfs[i].flux-med_image
    squozzed = np.squeeze(new, axis=0)
    residuals.append(squozzed)
    
plt.imshow(np.squeeze(tpfs[0].flux, axis=0))
plt.colorbar()
plt.title('raw')
plt.show()
    
med = np.median(residuals[0].value)
std = np.std(residuals[0].value)
plt.imshow(residuals[0].value)#, vmin=0, vmax=1)
plt.colorbar()
plt.title('residual')
plt.show()
    

fig = plt.figure()
gs = fig.add_gridspec(2,2)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])
ax3 = fig.add_subplot(gs[1,:])
im1 = ax1.imshow(med_image)#, cmap='Greys_r')
im2 = ax2.imshow(residuals[0], vmin=0, vmax=1)
im3 = ax3.scatter(new_lc.time.value, new_lc.flux.value, c='k', marker='.')

#x_init = (3, 7, 11)
#y_init = (9, 6, 3)
#x, y = centroid_sources(med_images[1420], x_init, y_init, centroid_func=centroid_com)
#x1, y1 = centroid_sources(sub_images[1420], x_init, y_init, centroid_func=centroid_com)

#print(x,y,x1,y1)

'''
multis_indices = [mi for mi, mx in enumerate(target_id) if mx == target]
sectrs = []
imy = 0
while imy < len(multis_indices):
    sectrs.append(target_sector[multis_indices[imy]])
    imy += 1
if sectrs.count(good_shot_sector[index]) > 1:
    for mult in multis_indices:
        if good_shot_sector[index] == target_sector[mult]:
            im3 = ax3.axvspan(target_start_mjd[mult],target_end_mjd[mult], color='b', alpha=0.5)
'''
#im3 = ax3.axvspan(good_shot_start[index], good_shot_end[index],color='b', alpha=0.5)
im3 = ax3.set_xlim(left=new_lc.time[0].value, right=new_lc.time[-1].value)
#im3 = ax3.axvline(x=lc.time[images_num[0]], ymin=-0.1, ymax=1.1, c='r')
fig.tight_layout()
fig.subplots_adjust(top=0.9)
def animate(i):
    im1.set_data(med_image)
    #ax1.set_title('{}'.format(images_num[i]))
    #ax2.set_title('MJD {}'.format(lc.time[images_num[i]]))
    im2.set_data(residuals[i])
    im3 = ax3.clear()
    im3 = ax3.scatter(new_lc.time.value, new_lc.flux.value, c='k', marker='.')
    #im3 = ax3.axvspan(good_shot_start[index], good_shot_end[index],color='b', alpha=0.5)
    #if sectrs.count(good_shot_sector[index]) > 1:
    #    for mult in multis_indices:
    #        if good_shot_sector[index] == target_sector[mult]:
    #            im3 = ax3.axvspan(target_start_mjd[mult],target_end_mjd[mult], color='b', alpha=0.5)
    im3 = ax3.set_xlim(left=new_lc.time[0].value, right=new_lc.time[-1].value)
    #im3 = ax3.axvline(x=new_lc.time[images_num[i]], ymin=-0.1, ymax=1.1, c='r')
    im3 = ax3.scatter(new_lc.time[i].value, new_lc.flux[i].value, c='orange', marker='*')
    return im1, im2, im3,
anim = animation.FuncAnimation(fig, animate, frames=np.arange(0, len(residuals), 1), blit=True)
#anim.save('/data/wallaby/rmorris/GALAH/Simult_Gifs/Feb9/{}_std2.mp4'.format(target))
plt.show()
#plt.close('all')
print('Animation Saved')    
    
    
    
    
    
    
    
    
    

print('done')