from astropy.io import fits

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


galah_dr3 = fits.open('GALAH_DR3_main_200331.fits')

###creates arrays of the source_id, ra, dec, and observation date for each target
source_ids = []
target_ra = []
target_dec = []
target_obs_date = []

printProgressBar(0, len(galah_dr3[1].data), prefix='Progress:', suffix='Complete',length=50)

for i, j in enumerate(galah_dr3[1].data):
    source_ids.append(j[1])
    target_ra.append(j[166])
    target_dec.append(j[168])
    printProgressBar(i + 1, len(galah_dr3[1].data), prefix='Progress:', suffix='Complete', length=50)

###observation date in YYMMDD format
printProgressBar(0, len(galah_dr3[1].data), prefix='Progress:', suffix='Complete',length=50)

for i, j in enumerate(source_ids):
	target_obs_date.append(int(str(j)[:6]))
	printProgressBar(i + 1, len(source_ids), prefix='Progress:', suffix='Complete', length=50)

