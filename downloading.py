import os
import csv
import numpy as np
from Spectra_EW import target



galah_obs_date = []
with open('/home/rmorris/documents/fieldlogs_export.txt', 'r') as file:
    reader = csv.reader(file,delimiter='\t')
    next(reader)
    for row in reader:
        galah_obs_date.append(row[1])
        
### Take only the fields from 2018-07-01 and later
galah_date_18 = []
for i, date in enumerate(galah_obs_date):
    if int(date) >= 180701:
        galah_date_18.append(date)
        
#print(galah_date_18)
#print(len(galah_date_18))

#links = ['https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/201003/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/200906/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/200901/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/200826/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/200825/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/200802/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/200728/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/200724/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/200714/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/200712/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/200708/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/190224/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/190223/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/190212/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/190210/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/190204/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/181225/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/181224/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/181223/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/181222/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/181221/spectra/com/com.tar.gz','https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/181220/spectra/com/com.tar.gz']
#dates = [201003,200906,200901,200826,200825,200802,200728,200724,200714,200712,200708,190224,190223,190212,190210,190204,181225,181224,181223,181222,181221,181220]
'''
for j, i in enumerate(galah_date_18[10:11]):
    os.system("wget --user rmorris1 --password 'Ezc23Pyp%5!%' -P /data/wallaby/rmorris/GALAH/Simult_GALAH_Data/ https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/{}/spectra/all/all.tar.gz".format(i))
    os.system("tar -xvzf /data/wallaby/rmorris/GALAH/Simult_GALAH_Data/all.tar.gz -C /data/wallaby/rmorris/GALAH/Simult_GALAH_Data/") # have to be in same folder as above?
    os.system("rm /data/wallaby/rmorris/GALAH/Simult_GALAH_Data/all.tar.gz")       # same as above^
    os.system("wget --user rmorris1 --password 'Ezc23Pyp%5!%' -P /data/wallaby/rmorris/GALAH/Simult_GALAH_Data/ https://cloud.datacentral.org.au/remote.php/webdav/GALAH/obs/reductions/Iraf_6.0/{}/spectra/com/com.tar.gz".format(i))
    os.system("tar -xvzf /data/wallaby/rmorris/GALAH/Simult_GALAH_Data/com.tar.gz -C /data/wallaby/rmorris/GALAH/Simult_GALAH_Data/") # have to be in same folder as above?
    os.system("rm /data/wallaby/rmorris/GALAH/Simult_GALAH_Data/com.tar.gz")       # same as above^
    print(j)

#link = 'https://cloud.datacentral.org.au/remote.php/webdav/GALAH/public/GALAH_DR3/GALAH_DR3_main_allspec_v2.fits'
#os.system("wget --user rmorris1 --password 'Ezc23Pyp%5!%' /home/rmorris/documents/ {}".format(link))

target(galah_date_18[10:11][0], 10**4, 0)
'''
#print(galah_date_18, len(galah_date_18))

res = []
for i in galah_date_18:
    if i not in res:
        res.append(i)
        
print(res, len(res))

print("Done")
    
    
    
    
    