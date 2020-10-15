import eleanor


### shutil.move function does not work right now, so to fix this, run eleanor.Update on any sector to
### download the target file. Move this file to .eleanor/metadata/s00XX and rename to target_s00XX
### then go to home/ryan/miniconda3/lib/python3.7/site-packages/eleanor and open update.py
### in update.py, comment out get_target (this is where shutil.move is called) lines 210-220, 94-97, 133
### open up code and run eleanor.Update again, this should download a cadences file and quality file
### should be good to use that sector now!
eleanor.Update(sector=1)