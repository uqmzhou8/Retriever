#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 10:16:11 2023

@author: uqmzhou8
"""


import pandas as pd
from cyvcf2 import VCF
import numpy as np
import random
from sklearn.ensemble import RandomForestClassifier
import multiprocessing 
import time


#reading of VCF file and storing the genotypes
def vcf_read(namef):
    pos=[]
    genos=[]
    for v in VCF(namef):
    
        pos.append(v.POS)
        gts = v.genotype.array().astype(int)
    
        if gts.shape[1] != 3:
            raise IOError('VCF must be for diploids.')
        #genotypes are converted in the following format: 0|0=1 0|1=2 1|0=3 1|1=4
        #missing data will be -3
        gts[gts[:, 2] == 1, 0] = 0
        gts[gts[:, 0] == 0, 0] = 2
        gts[gts[:, 0] == 1, 0] = 4
        gts = gts.sum(axis=1)
        gts = gts -1
        genos.append(gts)
    genos=np.stack(genos)
    pos=np.array(pos, dtype= "int")
    return genos, pos


# # extracting out a section of the genome and positions
def extract_ref(genos, pos, ref_size, window_size):

    

    i=0
    
    ref_gts=[]
    extract_gts=[]
    while i < len(pos):
    
        #find the range of sequence to extract
        #range determined by the window_size
        #Larger window_size will theoretically be more accurate as it takes a larger sequence into consideration
        nextpos=np.where(pos >= (pos[i] + window_size))
        
        ##if last position is smaller than the next position, the last position will be used for the window
        if len(nextpos[0]) == 0:
            a=len(pos)
        else:
            a=min(nextpos[0])
            
        while i < a:
            extract_gts = genos[i:a+1]
            missing_pos= np.where(extract_gts == -3)
            extract_gts = np.delete(extract_gts, missing_pos[1] , axis=1)
            # print(i)
            if len(extract_gts[0])<ref_size:
                a = a-1

            else:
                randomsam = random.sample(range(len(extract_gts[0])), ref_size)
                tempgts = extract_gts[: , randomsam]
                

                if i == 0:
                    ref_gts=tempgts
                    ref_pos= [pos[i:a+1]]
                else:
                    ref_pos.append(pos[i:a+1])
                    ref_gts= np.concatenate((ref_gts , tempgts), axis=0)
                i=a+1
    ref_pos=np.concatenate(ref_pos)
    if all(ref_pos != pos):
        raise IOError('Error in alignment')
    return ref_pos, ref_gts


def impute_missing(genos, pos, ref_gts, ref_pos, window_size, end_j, j=0):
    
    if all(pos!= ref_pos):
        raise IOError('Reference genome and imputed genome do not contain the same positions')
    imputed_gts=[]


    while j < len(pos):
        
        if j == 0:
            wind_max = np.where ( pos >= (pos[j] + window_size))
            if len(wind_max[0]) == 0:
                win_max=len(pos)
            else:
                wind_max=min(wind_max[0])
            train_fea = ref_gts[j+1 : wind_max]
            train_lab = ref_gts[j]
            test_fea = genos[j+1 : wind_max]
            test_lab = genos [j]
        elif j ==  len(pos):
            wind_min = np.where ( pos >= (pos[j] - window_size))
            wind_min= min(wind_min[0])
            train_fea = ref_gts[wind_min:j-1]
            train_lab = ref_gts[j]
            test_fea = genos[wind_min, j-1]
            test_lab = genos [j]
        else:
            wind_max = np.where ( pos >= (pos[j] + window_size))
            wind_min = np.where ( pos >= (pos[j] - window_size))
            if len(wind_max[0]) == 0:
                    wind_max=len(pos)
            else:
                    wind_max=min(wind_max[0])
            
            wind_min= min(wind_min[0])
            
            #perform segregrating the reference genome and imputed genome into training and predicting data
            train_fea = ref_gts[wind_min :j]
            train_fea = np.append(train_fea, ref_gts[j+1 : wind_max] , axis=0)
            train_lab = ref_gts[j]
            test_fea = genos[wind_min :j]
            test_fea = np.append(test_fea, genos[j+1 : wind_max] , axis=0)
            test_lab = genos[j]
            
        if test_fea.size == 0:
            if j<=21:
                train_fea = ref_gts[0 :j]
                train_fea = np.append(train_fea, ref_gts[j+1 : j+20] , axis=0)
                test_fea = genos[0 :j]
                test_fea = np.append(test_fea, genos[j+1 : j+20] , axis=0)
            elif j>= (len(pos)-21):
                train_fea = ref_gts[j-20 :j]
                train_fea = np.append(train_fea, ref_gts[j+1 : -1] , axis=0)
                test_fea = genos[j-20 :j]
                test_fea = np.append(test_fea, genos[j+1 : -1] , axis=0)

                
            train_fea = ref_gts[j-20 :j]
            train_fea = np.append(train_fea, ref_gts[j+1 : j+20] , axis=0)
            test_fea = genos[j-20 :j]
            test_fea = np.append(test_fea, genos[j+1 : j+20] , axis=0)
        train_fea=train_fea.T
        test_fea=test_fea.T
        test_lab=list(test_lab)
        
        train_lab = list(train_lab)
        rf = RandomForestClassifier(n_estimators = 100, n_jobs=1)
        rf.fit (train_fea, train_lab)
        test_pred = rf.predict(test_fea)
        
        imput_pos= np.where(np.array(test_lab) != -3)

        for x in imput_pos[0]:
            test_pred[int(x)] = test_lab[int(x)]
        #appending the column of imputed values into imputed gts
        imputed_gts.append(test_pred)
    
        j+=1
        
        if j== end_j+1:
            break
    imputed_gts = np.stack(imputed_gts)

    return imputed_gts

def main_impute(genos, pos, ref_gts, ref_pos, window_size, split_portions=10):

    if split_portions == 1:
        end_j=len(pos)
        imputed_gts = impute_missing(genos, pos, ref_gts, ref_pos, window_size, end_j)
        return imputed_gts
    else:

        no_j=np.array_split(np.array([*range(len(pos))]), split_portions)
        
        pool = multiprocessing.Pool(split_portions)
    
        processes = [pool.apply_async(impute_missing, args = (genos, pos, ref_gts, ref_pos, window_size, no_j[i][-1], no_j[i][0])) for i in range(split_portions)]
        result = [p.get() for p in processes]
        result= np.concatenate(result)
    return result
    
    
