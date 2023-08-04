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

    

    import pandas as pd
    import time
    


    import random
    from func import *
    import pandas as pd
    from cyvcf2 import VCF
    import numpy as np
        

    genos,pos= vcf_read('mango_chr20_subset_1.vcf')
    
    
    start_time = time.perf_counter()
    refpos, refgeno = extract_ref(genos, pos, 50, 100)

    imputed_gts = main_impute(genos, pos, refgeno, refpos, 100, 48)

    from collections import Counter
    findtruth = imputed_gts == genos
    cntmissing= Counter(np.concatenate(genos))
    cnt= Counter(np.concatenate(findtruth))
    
            
    finish_time = time.perf_counter()
    

    timetaken = finish_time-start_time
    
    print(f"Program finished in {timetaken} seconds")
    print(cnt)
    print(cntmissing)
    # tosto = pd.DataFrame([window, cnt, cntmissing, timetaken])
    # tosto.to_csv('read.csv', mode='a')