from astropy.io import fits
import math
hdul = fits.open('GALAH_DR3_main_200331.fits')

hdul.info()
print(hdul[1].header)
print(hdul[1].data[0][2])

#max_field_id = 7365
#for i in hdul[1].data:
#	if i[3] > max_field_id:
#		max_field_id = i[3]

#print(max_field_id)


# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

###prints the gaia id for all targets that do not have a field_id
#for i in hdul[1].data:
#	if i[3] < 0:
#		print(i[2])

###prints the Al_Fe ratio for all targets that have a valid number
#indices = [42,46,50,54,58,62,66,70,74,78,82,86,90,94,98,102,106,110,114,118,122,126,130,134,138,142,146,150,154,158,162]
#for i in hdul[1].data:
#	if math.isnan(i[62]) == False:
		#print(i[2])
		###print out all elemental abundances for targets that have Al_Fe abundance
#		Al_Fe_targets = []
#		j=0
#		while j < len(indices):
#			Al_Fe_targets.append(i[indices[j]])
#			if j == len(indices)-1:
#				print(Al_Fe_targets)
#			j+=1


###indices of elements in table
#indices = [42,46,50,54,58,62,66,70,74,78,82,86,90,94,98,102,106,110,114,118,122,126,130,134,138,142,146,150,154,158,162]
#gaia_ids = []
#i=0
###check if any targets have abundances for all elements
#printProgressBar(0, len(indices), prefix='Progress:', suffix='Complete',length=50)
#while i < len(hdul[1].data):
#	j=0
#	while j < len(indices):
#		if j == len(indices)-1 and math.isnan(hdul[1].data[i][indices[j]]) == False:
#			gaia_ids.append(hdul[1].data[i][2])
#			printProgressBar(i + 1, len(hdul[1].data), prefix='Progress:', suffix='Complete', length=50)
#			j+=10
#		elif math.isnan(hdul[1].data[i][indices[j]]) == False:
#			j+=1
#		else:
#			printProgressBar(i + 1, len(hdul[1].data), prefix='Progress:', suffix='Complete', length =50)
#			#print(hdul[1].data[i][2])
#			j+=len(indices)
#	i+=1

#print(gaia_ids)
#print("Out of {} targets, {} have abundances for {} elements".format(len(hdul[1].data),len(gaia_ids),len(indices)))

#indices_spec = [42,46,50,54,142,146,150,154,158,162]
#gaia_ids_spec = []
#i=0
###check which targets have abundances for X specific elements
#printProgressBar(0, len(indices_spec), prefix='Progress:', suffix='COmplete',length=50)
#while i < len(hdul[1].data):
#	j=0
#	while j < len(indices_spec):
#		if j == len(indices_spec)-1 and math.isnan(hdul[1].data[i][indices_spec[j]]) == False:
#			gaia_ids_spec.append(hdul[1].data[i][2])
#			printProgressBar(i + 1, len(hdul[1].data), prefix='Progress:', suffix='Complete', length=50)
#			j+=10
#		elif math.isnan(hdul[1].data[i][indices_spec[j]]) == False:
#			j+=1
#		else:
			#print(hdul[1].data[i][2])
#			printProgressBar(i + 1, len(hdul[1].data), prefix='Progress:', suffix='Complete', length=50)
#			j+=len(indices_spec)
#	i+=1

#print(gaia_ids_spec)
#print("Out of {} targets, {} have abundances for {} elements".format(len(hdul[1].data),len(gaia_ids_spec),len(indices_spec)))

###creates an array of just the source_id for each target, which includes date information
#source_ids = []
#for i in hdul[1].data:
#    source_ids.append(i[1])

###creates an array of just the 2018 source_ids
#source_ids_2018 = []
#for i in source_ids:
#    if i > 180000000000000:
#        source_ids_2018.append(i)

#print(source_ids_2018, len(source_ids_2018), max(source_ids_2018))



