#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:54:17 2023

@author: uqmzhou8
"""

#This is a script that shows how to use Retriever to perform the entire pipeline of assembly of chimeric references to the imputation using RF
#In the current stage, we recommend to use Beagle4 to perform the imputation process following Retriever's output of the chimeric reference panel


if __name__ == '__main__':

    #This will test the installation and running of Retriever to generate the chimeric reference panel using a dummy dataset
    #It should take less than a second to run as it is a very small dataset!
    from Retriever import *
    vcf2ref('dummy_dataset.vcf', 50, window_size=1000, outfile='chimeric_ref_gts.vcf')
