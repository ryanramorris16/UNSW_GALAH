import csv
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord

stream = fits.open('s5_pdr1_light.fits')
#print(stream[1].header)
ra = stream[1].data['ra']
dec = stream[1].data['dec']
priority = stream[1].data['priority']
gaia = stream[1].data['gaia_source_id']

#print(ra,dec)

ra_toi = []
dec_toi = []
tess_mag = []
tess_disp = []
tfop_disp = []
tic_id = []
with open('/home/rmorris/documents/exofop_tess_tois.csv', 'r') as file:
    reader = csv.reader(file,delimiter=',')
    next(reader)
    next(reader)
    next(reader)
    for row in reader:
        tic_id.append(row[0])
        ra_toi.append(row[18])
        dec_toi.append(row[19])
        tess_mag.append(row[14])
        tess_disp.append(row[12])
        tfop_disp.append(row[13])


match_s = []
match_t = []
for i in range(len(ra)):
    ra_s = ra[i]
    dec_s = dec[i]
    for j in range(len(ra_toi)):
        ra_t = ra_toi[j]
        dec_t = dec_toi[j]
        c1 = SkyCoord(ra_s, dec_s, unit = (u.deg, u.deg))
        c2 = SkyCoord(ra_t, dec_t, unit = (u.deg, u.deg))
        sep = c1.separation(c2)
        pixel_sep = sep.arcsecond/21
        if pixel_sep <= 8:
            match_s.append(i)
            match_t.append(j)

            
#print(match_s, match_t)

file = open('/home/rmorris/documents/stream_toi.txt', "a") 
for x in range(len(match_s)):
    s_ind = match_s[x]
    t_ind = match_t[x]
    file.write(str(tic_id[t_ind])+'\t'+str(gaia[s_ind])+'\t'+str(ra[s_ind])+'\t'+str(dec[s_ind])+'\t'+str(ra_toi[t_ind])+'\t'+str(dec_toi[t_ind])+'\t'+str(tfop_disp[t_ind])+'\t'+str(priority[t_ind])+'\n')
file.close()
   
    
    
    