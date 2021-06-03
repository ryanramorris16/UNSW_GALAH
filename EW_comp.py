import matplotlib.pyplot as plt
import numpy as np
import csv

date = 201003

gaia_ids = []
file_name = []
equivwidth = []
ew_err = []
with open('/home/rmorris/documents/spectra_ew_YYMMDD/spectra_ew_{}.txt'.format(date), 'r') as file:
    reader = csv.reader(file,delimiter='\t')
    for row in reader:
        gaia_ids.append(row[0])
        file_name.append(row[1])
        equivwidth.append(row[2])
        ew_err.append(row[3])

class galah_target:
    
    def __init__(self, EW, EW_err, files, ccd, gaia_id, date):
        self.files = files
        self.ccd = ccd
        self.fullEW = EW
        self.fullEWerr = EW_err
        self.gaia = gaia_id
        self.date = date
        
        self.separate()
        self.outlier()
        if self.gaia in self.marked:
            self.plotting()
            self.writeout()
        
        
    def separate(self):
        ccd_arr = self.ccd
        ccd_1_ind = [index for index, value in enumerate(ccd_arr) if value == 1]
        ccd_3_ind = [index for index, value in enumerate(ccd_arr) if value == 3]
        self.EW1 = EW[ccd_1_ind]
        self.EW1err = EW_err[ccd_1_ind]
        self.EW3 = EW[ccd_3_ind]
        self.EW3err = EW_err[ccd_3_ind]
    
    def outlier(self):
        marked_ids = []
        for k in [1, 3]:
            if k == 1:
                x = self.EW1
                y = self.EW1err
            elif k == 3:
                x = self.EW3
                y = self.EW3err
            ratio = []
            rat_e = []
            for i in range(len(x)):
                not_i = []
                not_i_e = []
                for j in range(len(x)):
                    if x[j] != x[i]:
                        not_i.append(x[j])
                        not_i_e.append(y[j])
                ij_ratio = []
                ij_e = []
                for j in range(len(not_i)):
                    ij_ratio.append(x[i]-not_i[j])    # not actually ratio, difference
                    ij_e.append((y[i]**2+not_i_e[j]**2)**0.5)
                    if ij_ratio[j]+3*ij_e[j] < 0.0 or ij_ratio[j]-3*ij_e[j] > 0.0:
                        print('Target {} is 3 sigma away from 1'.format(self.gaia))
                        if self.gaia not in marked_ids:
                            marked_ids.append(self.gaia)
                ratio.append(ij_ratio)
                rat_e.append(ij_e)
            
            if k == 1:
                self.ratios1 = ratio          # only gets the k=3 case because of overwriting
                self.ratios1_err = rat_e
            if k == 3:
                self.ratios3 = ratio          # only gets the k=3 case because of overwriting
                self.ratios3_err = rat_e
            self.marked = marked_ids
            print(k, "ratios: {}".format(ratio),"errors: {}".format(rat_e))
                
    def plotting(self):
        length = len(self.ratios1)
        fig = plt.figure()
        gs = fig.add_gridspec(2,1)
        ax1 = fig.add_subplot(gs[0,:])
        ax2 = fig.add_subplot(gs[1,:])
        for i in range(len(self.ratios1)):
            xs = np.zeros(len(self.ratios1[i]))
            for x in range(len(xs)):
                xs[x] = x*0.1
            ax1.errorbar(i+xs, self.ratios1[i], self.ratios1_err[i], marker='*', ms=5, capsize=5)
        for i in range(len(self.ratios1)):
            xs = np.zeros(len(self.ratios1[i]))
            for x in range(len(xs)):
                xs[x] = x*0.1
            ax2.errorbar(i+xs, self.ratios3[i], self.ratios3_err[i], marker='*', ms=5, capsize=5)
        plt.xlabel('Exposure')
        plt.ylabel('Ratio')
        ax1.set_title('{}'.format(self.gaia))
        #plt.show()
        plt.savefig('/data/wallaby/rmorris/GALAH/EW_ratios/{}.png'.format(self.gaia))
        plt.close()    
        
    def writeout(self):
        ### problem with writeout since all have dif EW1/3 lengths
        file = open('/home/rmorris/documents/ew_outliers.txt', "a")
        file.write(str(self.gaia)+'\t'+str(self.date)+'\n')
        file.close()
        
for j, i in enumerate(gaia_ids):
    num_exp = gaia_ids.count(i)         # number of exposures of this gaia id
    ind_exp = [index for index, value in enumerate(gaia_ids) if value == '{}'.format(i)]     # array of indices
    if j == ind_exp[0]:
        EW = np.zeros(len(ind_exp))
        EW_err = np.zeros(len(ind_exp))
        files = []
        ccd = np.zeros(len(ind_exp))
        for k, indices in enumerate(ind_exp):
            EW[k] = equivwidth[indices]
            EW_err[k] = ew_err[indices]
            files.append(file_name[indices])
            ccd[k] = file_name[indices][-6]
        #print(EW)
        #print(ccd)
        test = galah_target(EW, EW_err, files, ccd, i, date)
        #print(test.EW1)
    else:
        print('Target already measured')
    
print("All targets finished")
    
    
    
    
    
    