#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:54:17 2023

@author: uqmzhou8
"""

#This is a script that shows how to use Retriever to perform the entire pipeline of assembly of chimeric references to the imputation using RF
#In the current stage, we recommend to use Beagle4 to perform the imputation process following Retriever's output of the chimeric reference panel


if __name__ == '__main__':

    
    from Retriever import *
    vcf2imput("dog_chr1.vcf", 50, out_ref='chimeric_ref_gts_chr1.vcf', out_imputed='imputed_gts_chr1.vcf'):
