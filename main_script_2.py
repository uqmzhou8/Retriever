#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:54:17 2023

@author: uqmzhou8
"""

#Main script to perform preprocessing of VCF files to generate reference chimeric genomes

#creating n chimeric reference genomes from samples
# n= input('Please key in how many reference samples do you wish to generate '\
#          + 'from your samples.\n*Please note that the number of reference genomes '\
#              + 'is recommended to be 20-30% of the total number of samples.\n')


if __name__ == '__main__':
    windows= [ 100]
    

    import pandas as pd
    import time
    


    import random
    from func import *
    import pandas as pd
    from cyvcf2 import VCF
    import numpy as np
        

    genos,pos= vcf_read('mango_chr20_subset_1.vcf')
    
    for window in windows:
        print(window)

        start_time = time.perf_counter()
        refpos, refgeno = extract_ref(genos, pos, 50, window)

        
        imputed_gts = main_impute(genos, pos, refgeno, refpos, window, 48)
        # imputed_gts = impute_missing(genos[0:1000], pos[0:1000], refgeno, refpos, window)
        
    
        
        from collections import Counter
        findtruth = np.concatenate(imputed_gts) == genos
        cntmissing= Counter(np.concatenate(genos))
        cnt= Counter(np.concatenate(findtruth))
        
        

        
        finish_time = time.perf_counter()
        
        # end = time.time()
        timetaken = finish_time-start_time
        # print(f"Program finished in {finish_time-start_time} seconds")
        tosto = pd.DataFrame([window, cnt, cntmissing, timetaken])
        tosto.to_csv('read.csv', mode='a')