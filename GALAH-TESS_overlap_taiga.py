import csv
import glob
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from tess_stars2px import tess_stars2px_function_entry

### Opening the field logs for GALAH obs up to 2020-10-07
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

print('Fieldlogs Obtained')

### Take only the fields from 2018-07-01 and later
galah_date_18 = []
galah_field_18 = []
galah_ind_18 = []	#index of the galah fields to use in place of the field id
galah_ra_18 = []
galah_dec_18 = []
galah_time_utc = []
galah_exp_18 = []
for i, date in enumerate(galah_obs_date):
	year = str(date)[0:2]
	month = str(date)[2:4]
	day = str(date)[4:6]
	if int(date) >= 180701:
		galah_field_18.append(galah_field[i])
		galah_ind_18.append(i)
		galah_date_18.append(date)
		galah_ra_18.append(galah_obs_ra[i])
		galah_dec_18.append(galah_obs_dec[i])
		galah_exp_18.append(galah_exp[i])
		galah_time_utc.append('20{}-{}-{}T00:00:00'.format(year, month, day))

print('2018+ Fields Only')

### Converting RA and DEC to degrees 
c = SkyCoord(ra=galah_ra_18, dec=galah_dec_18, unit=(u.hourangle, u.deg))
galah_ra_18 = c.ra.degree
galah_dec_18 = c.dec.degree

print('RA and DEC Converted')

outID,outEclipLong,outEclipLat,outSec,outCam,outCcd,outColPix,outRowPix,scinfo = tess_stars2px_function_entry(galah_ind_18, galah_ra_18, galah_dec_18)
for i in range(len(outID)):
	print('{0:d} {1:d} {2:d} {3:d}'.format(outID[i],outSec[i],outCam[i],outCcd[i]))

### Opening TESS orbit times to compare with GALAH times
tess_sector = []
tess_orbit = []
tess_start_tjd = []
tess_end_tjd = []
with open('/home/rmorris/documents/orbit_times_20201013_1338.txt', 'r') as file:
	reader = csv.reader(file,delimiter='\t')
	for row in reader:
		tess_sector.append(row[0])
		tess_orbit.append(row[1])
		tess_start_tjd.append(row[4])
		tess_end_tjd.append(row[5])

print('TESS Orbit Information Acquired')

### Converting times from UTC to MJD for easier comparison later
galah_time_mjd = []
tess_start_mjd =[]
tess_end_mjd = []
for i, time in enumerate(galah_time_utc):
	gtime_utc = Time(galah_time_utc[i], format='isot', scale='utc')
	galah_time_mjd.append(gtime_utc.mjd)

for i, time in enumerate(tess_start_tjd):
	tsmjd = float(tess_start_tjd[i])+57000.5
	temjd = float(tess_end_tjd[i])+57000.5
	tess_start_mjd.append(tsmjd)
	tess_end_mjd.append(temjd)

print('Times Converted to MJD')

### Comparing GALAH times with TESS Point observations on field by field basis
simult_field = []
simult_date = []
simult_exp = []
simult_sector = []
for i, ids in enumerate(outID):
	field_index = galah_ind_18.index(ids)
	sector_num = outSec[i]
	ind_low = (2*sector_num)-2
	ind_high = (2*sector_num)-1
	if sector_num <=35:
		if tess_start_mjd[ind_low] <= galah_time_mjd[field_index] <= tess_end_mjd[ind_low]:
			simult_field.append(galah_field_18[field_index])
			simult_date.append(galah_date_18[field_index])
			simult_exp.append(galah_exp_18[field_index])
			simult_sector.append(sector_num)
			print("This target {} is in the first orbit of sector {}".format(galah_field_18[field_index],sector_num))
		elif tess_start_mjd[ind_high] <= galah_time_mjd[field_index] <= tess_end_mjd[ind_high]:
			simult_field.append(galah_field_18[field_index])
			simult_date.append(galah_date_18[field_index])
			simult_exp.append(galah_exp_18[field_index])
			simult_sector.append(sector_num)
			print("This target {} is in the second orbit of sector {}".format(galah_field_18[field_index],sector_num))
		#else:
			#print("This target was not observed during this sector")


### Finding the start and end exposures for each fits
start_exp = []
end_exp = []
for i, exp in enumerate(simult_exp):
	start_exp.append(str(exp.split('-')[0]).zfill(2))
	end_exp.append(str(exp.split('-')[1]).zfill(2))

### Code to find the exact time of GALAH observations
fits_start_mjd = []
fits_end_mjd = []
sequence_start_mjd = []
sequence_end_mjd = []
for i, fields in enumerate(simult_field):
	date = simult_date[i]
	exp_range = np.arange(int(start_exp[i]),int(end_exp[i])+1)
	### Opening up the specific fits files for a given field
	for num in np.arange(1,5):
		for exp in exp_range:
			new_exp = str(exp).zfill(2)
			fits_file = glob.glob('/data/wallaby/rmorris/GALAH/{}/data/ccd_{}/*{}.fits'.format(date,num,new_exp))
			if not fits_file:
				print('There is no fits file for these specifications')
			### Converting UT Start and End times in fits files to MJD
			elif fits_file:
				print('The fits file location:', fits_file)
				path = Path(fits_file[0])
				galah_fits = fits.open(path)
				utdate = galah_fits[0].header['UTDATE'].replace(":","-")
				utstart = galah_fits[0].header['UTSTART']
				utend = galah_fits[0].header['UTEND']
				utc_full_start = str(utdate)+"T"+str(utstart)
				utc_full_end = str(utdate)+"T"+str(utend)
				utc_object_s = Time(utc_full_start, format='isot', scale='utc')
				utc_object_e = Time(utc_full_end, format='isot', scale='utc')
				fits_start_mjd.append(utc_object_s.mjd)
				fits_end_mjd.append(utc_object_e.mjd)
			if exp == exp_range[0] and num == 1:
				sequence_start_mjd.append(utc_object_s.mjd)
			if exp == exp_range[-1] and num == 1:
				sequence_end_mjd.append(utc_object_e.mjd)

print('GALAH observation times acquired')


### Opening fld files for each field to extract targets
s_targets = []
s_ra = []
s_dec = []
s_field_sector = []
s_index = []
s_mjd_s = []
s_mjd_e = []
s_ref_field = []
for i, fields in enumerate(simult_field):
	date = simult_date[i]
	print(date)
	try:
		try:
			### Actually opening each specific fld file to get target information
			field_fld = fields.replace("_p0","").replace("_p1","").replace("_r1","").replace("_r2","")
			with open('/data/wallaby/rmorris/GALAH/{}/fld/{}.fld'.format(date,field_fld), 'r') as file:
				reader = csv.reader(file,delimiter=' ')
				next(reader)
				next(reader)
				next(reader)
				next(reader)
				next(reader)
				next(reader)
				target_names = []
				target_ra = []
				target_dec = []
				target_desig = []
				for row in reader:
					target_nums = []
					for value in row:
						if value:
							target_nums.append(value)
					target_names.append(target_nums[0])
					ra_hms = target_nums[1]+'h'+target_nums[2]+'m'+target_nums[3]+'s'
					dec_dms = target_nums[4]+'d'+target_nums[5]+'m'+target_nums[6]+'s'
					target_desig.append(target_nums[7])
					c = SkyCoord(ra=ra_hms, dec=dec_dms, unit=(u.hourangle, u.deg))
					target_ra.append(c.ra.degree)
					target_dec.append(c.dec.degree)
		except FileNotFoundError as err:
			print(err)	
		### Sorting out the sky pictures from the real targets
		for j, designation in enumerate(target_desig):
			if str(designation) != 'S' and 's':
				s_targets.append(target_names[j]) 
				s_index.append(j)
				s_ra.append(target_ra[j])
				s_dec.append(target_dec[j])
				s_field_sector.append(simult_sector[i])
				s_ref_field.append(fields)
				s_mjd_s.append(sequence_start_mjd[i])
				s_mjd_e.append(sequence_end_mjd[i])
	except:
		print('Oops, something broke')

print('Individual targets acquired')

### Matching targets specifically to sectors in TESS
s_sector = []
for i, target in enumerate(s_targets):
	print(target)
	outID,outEclipLong,outEclipLat,outSec,outCam,outCcd,outColPix,outRowPix,scinfo = tess_stars2px_function_entry(s_index[i], s_ra[i], s_dec[i])
	s_sector.append(outSec)
	print(outSec)

print('Sectors found for all targets')

### Comparing exact times of GALAH and TESS observations
real_s_targets = []
real_s_ra = []
real_s_dec = []
real_s_sector = []
real_s_mjd_s = []
real_s_mjd_e = []
for i, target in enumerate(s_targets):
	sector_num = s_sector[i]
	for sectors in sector_num:
		print(sectors)
		ind_low = (2*sectors)-2
		ind_high = (2*sectors)-1
		if sectors <=35:
			if tess_start_mjd[ind_low] <= s_mjd_s[i] <= tess_end_mjd[ind_low]:
				real_s_targets.append(target)
				real_s_ra.append(s_ra[i])
				real_s_dec.append(s_dec[i])
				real_s_sector.append(sectors)
				real_s_mjd_s.append(s_mjd_s[i])
				real_s_mjd_e.append(s_mjd_e[i])
				print("This target, {}, was simultaneously observed in the first orbit of sector {}".format(target,sectors))
			elif tess_start_mjd[ind_high] <= s_mjd_e[i] <= tess_end_mjd[ind_high]:
				real_s_targets.append(target)
				real_s_ra.append(s_ra[i])
				real_s_dec.append(s_dec[i])
				real_s_sector.append(sectors)
				real_s_mjd_s.append(s_mjd_s[i])
				real_s_mjd_e.append(s_mjd_e[i])
				print("This target, {}, was simultaneously observed in the second orbit of sector {}".format(target,sectors))
			else:
				print("This target, {}, was not simultaneously observed during this sector".format(target))

print(len(real_s_targets))

### Writes the targets and their information to a txt file for lightcurves later
with open('/home/rmorris/documents/simult_obs_targets.txt', 'w') as file:
	for number, values in enumerate(real_s_targets):
		file.write(str(values)+'\t'+str(real_s_ra[number])+'\t'+str(real_s_dec[number])+'\t'+str(real_s_sector[number])+'\t'+str(real_s_mjd_s[number])+'\t'+str(real_s_mjd_e[number])+'\n')
