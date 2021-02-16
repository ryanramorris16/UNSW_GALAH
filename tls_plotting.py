import csv
import matplotlib.pyplot as plt
import numpy as np
import math
import eleanor
import lightkurve as lk

'''
exo_per = []
tls_per = []
with open('/home/rmorris/documents/tlst_bkg_period_small_sigma.txt', 'r') as file:
	reader = csv.reader(file,delimiter='\t')
	for row in reader:
		tls_per.append(float(row[1]))
		exo_per.append(float(row[0]))
'''
tic_ids = []
with open('/home/rmorris/documents/exofop_tess_tois_faint_no_FP.txt', 'r') as file:
	reader = csv.reader(file,delimiter='\t')
	for row in reader:
		tic_ids.append(int(row[0].replace("'","")))


exo_per = []
exo_error = []
tls_per = []
tls_error = []
with open('/home/rmorris/documents/spoc_tls_periods.txt', 'r') as file:
	reader = csv.reader(file,delimiter='\t')
	for row in reader:
		exo_per.append(float(row[0]))
		exo_error.append(float(row[1]))
		tls_per.append(float(row[2]))
		tls_error.append(float(row[3]))

harmonic_id = []
for i, period in enumerate(tic_ids):
	for n in np.arange(2,16,1):
		if tic_ids[i] not in harmonic_id:
			exo_upper = exo_per[i] + ((3)*exo_error[i])
			exo_lower = exo_per[i] - ((3)*exo_error[i])
			if (tls_per[i]*n)+(3*tls_error[i]) > exo_lower and (tls_per[i]*n)-(3*tls_error[i]) < exo_upper:
				print("{} is a harmonic where {} times {} is the true period".format(tic_ids[i],tls_per[i],n))
				harmonic_id.append(tic_ids[i])
			elif (tls_per[i]/n)+(3*tls_error[i]) > exo_lower and (tls_per[i]/n)-(3*tls_error[i]) < exo_upper:
				print("{} is a harmonic where {} divded by {} is the true period".format(tic_ids[i],tls_per[i],n))
				harmonic_id.append(tic_ids[i])

bad_targets = []
bad_tls = []
bad_tls_err = []
bad_exo = []
bad_exo_err = []
good_targets = []
good_tls = []
good_tls_err = []
good_exo = []
good_exo_err = []
harmonic_counter = 0
wtf_targets = []
for i, period in enumerate(tls_per):
	tls_upper = tls_per[i] + (3*tls_error[i])
	tls_lower = tls_per[i] - (3*tls_error[i])
	exo_upper = exo_per[i] + (3*exo_error[i])
	exo_lower = exo_per[i] - (3*exo_error[i])
	if tls_upper >= exo_lower and tls_lower <= exo_upper:
		#if period < exo_per[i]*0.9 or period > exo_per[i]*1.1:
		good_targets.append(tic_ids[i])
		good_tls.append(period)
		good_exo.append(exo_per[i])
		good_tls_err.append(tls_error[i])
		good_exo_err.append(exo_error[i])
		print("Good target added")
	#if tic_ids[i] not in good_targets and tic_ids[i] not in bad_targets:
	elif tic_ids[i] in harmonic_id:
	#check if this target is a harmonic
		good_targets.append(tic_ids[i])
		good_tls.append(period)
		good_exo.append(exo_per[i])
		good_tls_err.append(tls_error[i])
		good_exo_err.append(exo_error[i])
		harmonic_counter += 1
		print("Harmonic target added")

	elif tic_ids[i] not in harmonic_id:
		if tls_upper < exo_lower or tls_lower > exo_upper:
		#if period < exo_per[i]*0.9 or period > exo_per[i]*1.1:
			bad_targets.append(tic_ids[i])
			bad_tls.append(period)
			bad_exo.append(exo_per[i])
			bad_tls_err.append(tls_error[i])
			bad_exo_err.append(exo_error[i])
			print("Bad target added")
			#print(period, exo_per[i])
		else:
			print("Target got lost")

#print(len(tls_per))
#print(wtf_targets)
print("You have found {} good targets where {} are harmonic targets, and {} bad targets out of {} total targets".format(len(good_targets),harmonic_counter,len(bad_targets), len(tic_ids)))
#print(harmonic_id)

#print(len(tic_ids), len(harmonic_id))

plt.figure(num=1)
plt.scatter(bad_exo,bad_tls, c='r')
plt.xlabel('ExoFOP Period')
plt.ylabel('TLS Period')
#plt.show()

bad_per_res = []
for i, period in enumerate(bad_tls):
	bad_per_res.append(period/bad_exo[i])

#print(bad_per_res)

plt.figure(num=2)
plt.scatter(bad_exo, bad_per_res)
plt.xlabel('ExoFOP Period')
plt.ylabel('TLS / Exo')
plt.show()

for i, j in enumerate(bad_targets):
	if bad_tls[i] != 0:
		print(j)
		star = eleanor.multi_sectors(tic=j, sectors='all')
		data = []
		for s in star:
			    datum = eleanor.TargetData(s, do_psf=False, do_pca=False)
			    data.append(datum)

		time = []
		flux = []
		for sector, datum in enumerate(data):
		    q = datum.quality == 0
		    #plt.plot(datum.time[q], datum.corr_flux[q]/np.median(datum.corr_flux[q]), plot_fmt[sector])
		    time.append(datum.time[q])
		    flux.append(datum.corr_flux[q]/np.median(datum.corr_flux[q]))

		time_array = np.concatenate((time[0:(len(time))]))
		flux_array = np.concatenate((flux[0:(len(flux))]))

		lc = lk.LightCurve(time = time_array, flux = flux_array).remove_outliers(sigma_lower=7, sigma_upper=3).flatten(window_length=51)
		#plt.figure(num=str(j) + "_1")
		#plt.scatter(lc.time,lc.flux, c='k', marker='.')
		#plt.title(j)
		#plt.xlabel('Time')
		#plt.ylabel('Normalized Flux')
		#plt.show()

		fold_lc = lc.fold(period=bad_tls[i])#, t0=periodogram.transit_time_at_max_power)
		fold_lc_real = lc.fold(period=bad_exo[i]).bin(binsize=5)
		plt.figure(num=str(j) + "_3")
		plt.scatter(fold_lc.time, fold_lc.flux, c='k', marker='.')
		plt.plot(fold_lc_real.time, fold_lc_real.flux, c='r', marker='.')
		plt.title(i)
		plt.xlabel('Phase')
		plt.ylabel('Normalized Flux')
		plt.legend(['ExoFOP Folded {}'.format(bad_exo[i]), 'TLS Folded {}'.format(bad_tls[i])])
		plt.savefig('/home/rmorris/documents/Faint_Bad_Sigma_QLP_Plots/{}_folded_lc.png'.format(bad_targets[i]))

