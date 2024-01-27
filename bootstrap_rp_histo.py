from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from scipy.stats import binned_statistic
from scipy import stats
import numpy as np           # to define our table
import sys
import math                  
#import decimal               
#import numpy.ma as ma
#import random
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec

BG_li=[] # barred galaxies
NBG_li=[] # unbarred galaxies
colorB_li=[]
colorNB_li=[]
rpB_li=[]
rpNB_li=[]
MassB_li=[]
MassNB_li=[]
SFRB_li=[]
SFRNB_li=[]
MtypeB_li=[]
MtypeNB_li=[]

d='disk_revised_weight.dat'
data=np.loadtxt(d,dtype={'names':('RA','Dec','WINGS','Cluster','SB','eps', 'pa','Mtype','Bar', 'weight','z', 'logsSFR','logM', 'Mv','BV','r200','logM200','D','rp'),'formats': ('O','O','O','O','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f')} ,skiprows=1)


SB=data['SB']
Bar=data['Bar']
Mtype=data['Mtype']
z=data['z']
eps=data['eps']
BV=data['BV']
rp=data['rp']
SFR=data['logsSFR']
M_cluster=data['logM200']
Mass=data['logM']


for i in range(len(Bar)):
    if Bar[i]==0:
       NBG=Bar[i]
       colorNB=BV[i]
       rpNB=rp[i]
       MassNB=Mass[i]
       SFRNB=SFR[i]
       MtypeNB=Mtype[i]
       NBG_li.append(NBG)
       colorNB_li.append(colorNB)
       rpNB_li.append(rpNB)
       MassNB_li.append(MassNB)
       SFRNB_li.append(SFRNB)
       MtypeNB_li.append(MtypeNB)
    else:
        BG=Bar[i]
        colorB=BV[i]
        rpB=rp[i]
        MassB=Mass[i]
        SFRB=SFR[i]
        MtypeB=Mtype[i]
        BG_li.append(BG)
        colorB_li.append(colorB)
        rpB_li.append(rpB)
        MassB_li.append(MassB)
        SFRB_li.append(SFRB)
        MtypeB_li.append(MtypeB)
        
print('barred',len(BG_li))
print('unbarred',len(NBG_li))
 





## compute max and min values of the original table to subset the whole range into equal intervals#


x=rp
y=Bar

binx=np.linspace(0,2.2,7)  
print('bins',binx)


bin_mean, bin_edges, binnumber = stats.binned_statistic(x, y, statistic='mean', bins=binx)  ## this will give bin_mean, bin_edges, binnumber

##note that mean values represents the bar fraction where barred galaxies have "1" and unbarred are "0"

bin_width = (bin_edges[1] - bin_edges[0])

bin_centers = bin_edges[1:] - bin_width/2

print('Bar_fraction_dat',bin_mean)


x_data_li=[]
y_data_li=[]

for k in range(len(bin_centers)):
    x_data=bin_centers[k]
    y_data=bin_mean[k]

    x_data_li.append(bin_centers[k])
    y_data_li.append(bin_mean[k])

    print ('x_data', x_data_li)
    print ('y_data', y_data_li)

 



## apply bootstrap resampling

#bootarr = np.array(Bar) ## no need cause Bar is already an array#

#test_statistic = lambda x: (np.sum(x), np.mean(x))

Nsamp=1000
Nbin=binx

with NumpyRNGContext(1):
    
    x_resamp = bootstrap(x,Nsamp)  # 100 is the number of repeted sample ## color should be array
    y_resamp = bootstrap(y,Nsamp)

    #print('xbootresult:',x_resamp)
    
    print('xbootshape:',x_resamp.shape) # to make sure that the shape of the new sample is the same as the original table (100,4359)
    
    #print('ybootresult:',y_resamp)
    
    print('ybootshape:',y_resamp.shape)
    
    
    #sys.exit()




    Bin1_li=[]
    Bin2_li=[]
    Bin3_li=[]
    Bin4_li=[]
    Bin5_li=[]
    Bin6_li=[]
    Bin7_li=[]
    Bin8_li=[]
    #for i, x_resamp in enumerate(x_resamp):
    for i in range(Nsamp):
        bin_mean_z,bin_edges_z,binnumber_z=stats.binned_statistic(x_resamp[i], y_resamp[i], statistic='mean', bins=Nbin)
        print('z_mean',bin_mean_z)
        for j in Nbin:
            Bin1=bin_mean_z[0]
            print('bin1',Bin1)
            Bin1_li.append(Bin1)
            
            Bin2=bin_mean_z[1]
            Bin2_li.append(Bin2)
            
            Bin3=bin_mean_z[2]
            Bin3_li.append(Bin3)
            
            Bin4=bin_mean_z[3]
            Bin4_li.append(Bin4)
            
            Bin5=bin_mean_z[4]
            Bin5_li.append(Bin5)
            
            Bin6=bin_mean_z[5]
            Bin6_li.append(Bin6)
            
            
            
    std_bin1=np.nanstd(Bin1_li)
    print('std_bin1',std_bin1)
    std_bin2=np.nanstd(Bin2_li)
    print('std_bin2',std_bin2)
    std_bin3=np.nanstd(Bin3_li)
    print('std_bin3',std_bin3)
    std_bin4=np.nanstd(Bin4_li)
    print('std_bin4',std_bin4)
    std_bin5=np.nanstd(Bin5_li)
    print('std_bin5',std_bin5)
    std_bin6=np.nanstd(Bin6_li)
    print('std_bin6',std_bin6)
    
    
    error_bar=[std_bin1,std_bin2,std_bin3,std_bin4,std_bin5,std_bin6]
    print('error',error_bar)

            
        
#plt.figure()
#plt.errorbar(x_data_li,y_data_li, yerr= error_bar,fmt='ko-')
#plt.ylim(0,0.5)
#plt.xlim(0.1,1.5)
#plt.xlabel('eps')
#plt.ylabel('Bar fraction')
#plt.savefig('fbar_eps2_weight_revised_data.png')
#plt.close()

rp_li=[]

for i in range (len(rp)):
    if 0 < rp[i] <2.2:
        newrp=rp[i] 
        rp_li.append(newrp)  

fig, axs = plt.subplots(figsize=(6,6)) 

gs = gridspec.GridSpec(2,1, height_ratios= [1,2])

ax1 = plt.subplot(gs[1,:])
plt.errorbar(x_data_li,y_data_li, yerr= error_bar,fmt='ks-')
plt.setp(ax1.get_xticklabels(), fontsize=15)
plt.setp(ax1.get_yticklabels(), fontsize=15)
plt.xlabel('(r$_{p}/$r$_{virial})_{cl}$',fontsize=20) 
plt.ylabel('$f_{bar}$',fontsize=20)
plt.xlim(0,2.2)
plt.ylim(0.15,0.48)
plt.xticks(np.arange(0, 2.2, step=0.5),fontsize=15)

# share x only
ax2 = plt.subplot(gs[0, :], sharex=ax1)
data1= (rpB_li)
data2= (rpNB_li)
num_bins=25
plt.hist(data1, num_bins, histtype='step', color='k', label='Barred')
plt.hist(data2, num_bins, histtype='step', linestyle=('dashed'), color='b', label='unbarred')
#plt.xlabel('r$_p$/r$_{virial}$', fontsize=18)
plt.ylabel('$N$',fontsize=20)
plt.legend(fontsize=16,loc='upper right')
#plt.hist(rp_li,bins=20,histtype='step', color='b')  
# make these tick labels invisible
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), fontsize=15)


plt.subplots_adjust(hspace=.0) ## remove the gap between the two plots
plt.subplots_adjust(left=0.2, bottom=0.15, top=0.95)

plt.savefig('fbar_rp_histo.pdf')
plt.close()        
        
            
        
        
            
        

    
    
    
    

   
