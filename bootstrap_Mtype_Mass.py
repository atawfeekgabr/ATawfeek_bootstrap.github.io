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


x=np.array(MtypeB_li)
y=np.array(MassB_li)

binx=np.linspace(np.min(x),np.max(x),8)  
print('bins',binx)


bin_mean, bin_edges, binnumber = stats.binned_statistic(x, y, statistic='median', bins=binx)  ## this will give bin_mean, bin_edges, binnumber

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
    
    #print('xbootshape:',x_resamp.shape) # to make sure that the shape of the new sample is the same as the original table (100,4359)
    
    #print('ybootresult:',y_resamp)
    
    #print('ybootshape:',y_resamp.shape)
    
    
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
        bin_mean_z,bin_edges_z,binnumber_z=stats.binned_statistic(x_resamp[i], y_resamp[i], statistic='median', bins=Nbin)
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
            
            Bin7=bin_mean_z[6]
            Bin7_li.append(Bin7)
            
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
    std_bin7=np.nanstd(Bin7_li)
    print('std_bin7',std_bin7)
    
    error_bar=[std_bin1,std_bin2,std_bin3,std_bin4,std_bin5,std_bin6,std_bin7 ]
    print('error',error_bar)


x2=np.array(MtypeNB_li)
y2=np.array(MassNB_li)

binx2=np.linspace(np.min(x2),np.max(x2),8)  
print('bins',binx)


bin2_mean, bin2_edges, binnumber = stats.binned_statistic(x2, y2, statistic='median', bins=binx2)  ## this will give bin_mean, bin_edges, binnumber

##note that mean values represents the bar fraction where barred galaxies have "1" and unbarred are "0"

bin2_width = (bin2_edges[1] - bin2_edges[0])

bin2_centers = bin2_edges[1:] - bin2_width/2

#print('Bar_fraction_dat',bin_mean)


x2_data_li=[]
y2_data_li=[]

for k in range(len(bin2_centers)):
    x2_data=bin2_centers[k]
    y2_data=bin2_mean[k]

    x2_data_li.append(bin2_centers[k])
    y2_data_li.append(bin2_mean[k])

    print ('x2_data', x2_data_li)
    print ('y2_data', y2_data_li)

 



## apply bootstrap resampling

#bootarr = np.array(Bar) ## no need cause Bar is already an array#

#test_statistic = lambda x: (np.sum(x), np.mean(x))

Nsamp=1000
Nbin=binx2

with NumpyRNGContext(1):
    
    x2_resamp = bootstrap(x2,Nsamp)  # 100 is the number of repeted sample ## color should be array
    y2_resamp = bootstrap(y2,Nsamp)

    #print('xbootresult:',x_resamp)
    
    #print('xbootshape:',x_resamp.shape) # to make sure that the shape of the new sample is the same as the original table (100,4359)
    
    #print('ybootresult:',y_resamp)
    
    #print('ybootshape:',y_resamp.shape)
    
    
    #sys.exit()




    Bin1_li2=[]
    Bin2_li2=[]
    Bin3_li2=[]
    Bin4_li2=[]
    Bin5_li2=[]
    Bin6_li2=[]
    Bin7_li2=[]
    Bin8_li2=[]
    #for i, x_resamp in enumerate(x_resamp):
    for i in range(Nsamp):
        bin2_mean_z,bin2_edges_z,binnumber_z=stats.binned_statistic(x2_resamp[i], y2_resamp[i], statistic='median', bins=Nbin)
        print('z_mean',bin_mean_z)
        for j in Nbin:
            Bin12=bin2_mean_z[0]
            #print('bin1',Bin1)
            Bin1_li2.append(Bin12)
            
            Bin22=bin2_mean_z[1]
            Bin2_li2.append(Bin22)
            
            Bin32=bin2_mean_z[2]
            Bin3_li2.append(Bin32)
            
            Bin42=bin2_mean_z[3]
            Bin4_li2.append(Bin42)
            
            Bin52=bin2_mean_z[4]
            Bin5_li2.append(Bin52)
            
            Bin62=bin2_mean_z[5]
            Bin6_li2.append(Bin62)
            
            Bin72=bin2_mean_z[6]
            Bin7_li2.append(Bin72)
            
    std_bin12=np.nanstd(Bin1_li2)
    #print('std_bin1',std_bin1)
    std_bin22=np.nanstd(Bin2_li2)
    #print('std_bin2',std_bin2)
    std_bin32=np.nanstd(Bin3_li2)
    #print('std_bin3',std_bin3)
    std_bin42=np.nanstd(Bin4_li2)
    #print('std_bin4',std_bin4)
    std_bin52=np.nanstd(Bin5_li2)
    #print('std_bin5',std_bin5)
    std_bin62=np.nanstd(Bin6_li2)
    #print('std_bin6',std_bin6)
    std_bin72=np.nanstd(Bin7_li2)
    #print('std_bin7',std_bin7)
    
    error_bar2=[std_bin12,std_bin22,std_bin32,std_bin42,std_bin52,std_bin62,std_bin72 ]
    print('error2',error_bar2)


        
plt.figure(1)
plt.errorbar(x_data_li,y_data_li, yerr= error_bar,fmt='ko-',label='Barred')
plt.errorbar(x2_data_li,y2_data_li, yerr= error_bar2,fmt='bs-',label='Unbarred')
#plt.ylim(0,0.5)
#plt.xlim(0.1,1.5)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=18)
plt.xlabel('Mtype',fontsize=20)
plt.ylabel('log($M_{\star}/M_{\odot}$)',fontsize=20)
plt.legend(fontsize=10, loc='upper left')
plt.subplots_adjust(left=0.25, bottom=0.2, top=0.85)
plt.savefig('Mtype_Mass_median.pdf')
plt.close()



    
    
    

   
