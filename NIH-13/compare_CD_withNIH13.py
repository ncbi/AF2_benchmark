#! /usr/local/bin/python3.10/

import sys
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

if __name__ == '__main__':

    f = pd.read_csv(sys.argv[1], sep="\t")
    print(f)

    print(f.iloc[:,0][1:].to_numpy())
    colors = ['#C74B40','#00AAAA']
    #m = ['.','o','s','*','^','v']
    for i in range(10):
        if i < 6:
            c = colors[1]
            #mk = m[i]
            #if i == 2:
            #    c = '#555555'
        else:
            c= colors[0]
            #mk = m[i-6]
            #if i == 8:
                #c = "#CFB997"

        if i == 7:
            plt.plot(f.iloc[:,0][1:].to_numpy(),f.iloc[:,i+1][1:].to_numpy(),marker='.',color=c,markersize=3)
        else:
            plt.plot(f.iloc[:,0][1:].to_numpy(),f.iloc[:,i+1][1:].to_numpy(),marker='.',color=c,markersize=3,alpha=0.15)

    
    plt.plot(f.iloc[:,0][1:].to_numpy(),f.iloc[:,11][1:].to_numpy(),marker='o',color='k',markersize=3)   
            
        #plt.title("C-terminal domain only",fontsize=16)
    plt.xlabel("Wavelength (nm)",fontsize=16)
    plt.title("Full-length proteins",fontsize=16)
        #plt.ylabel('['+r'$\Theta$'+r'] ($deg*cm^2*dmol^{-1})$',fontsize=16)
    plt.ylabel('Mean Residue Ellipticity\n'+r'($deg*cm^2*dmol^{-1})$',fontsize=16)
    plt.xlim([195,250])
    plt.ylim([-18000,18000])
    plt.gcf().subplots_adjust(left=0.20)
    plt.legend(list(np.array(range(10))+1)+['13'],fontsize=10,loc='upper center', \
          ncol=5, fancybox=True, shadow=True)
    plt.savefig("Whole_CD_final_withNIH13.png",dpi=600)
   
   
    #plt.plot(x,ER,'-',color="#3CA1AF")
    #plt.plot(x,KR,'-',color="#FC7802")
    #plt.plot(x,KN,'-',color="#6F286B")
    #plt.plot(x,TR,'-',color="#0D7943")
    #plt.title("Whole protein circular dichroism")
    #plt.xlabel("Wavelength (nm)")
    #plt.ylabel('['+r'$\Theta$'+r'] (deg*$cm^2*dmol^{-1})$')
    #plt.xlim([195,280])
    #plt.ylim([-18000,18000])
    #plt.gcf().subplots_adjust(left=0.15)
    #plt.legend(['ERfaH', 'KNusG', 'KRfaH', 'TRfaH'])
    #plt.legend(['ERfaH', 'KRfaH'])
    #plt.savefig("CD_whole_EK.png",dpi=600)

        
