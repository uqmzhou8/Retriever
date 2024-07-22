# Retriever
The program uses a VCF input file to create a predefined number of chimeric reference genomes. These reference genomes can be used for imputation with either conventional imputation software or a built-in algorithm that employs machine learning. Although Retriever has been designed to be user-friendly, a basic understanding of coding with Python is assumed, requring users to execute the imputation pipeline through simple command lines with minimal additional programming.

# Installation:
Retriever makes use of external packages for some processes so it is recommended to install these packages in a virtual environment:
> conda create -n Retriever
> 
> conda activate Retriever

Please install the packages required by typing the following line  into the powershell or terminal:
>conda install numpy pandas cyvcf2 scikit-learn

These packages are prerequsites for running Retriever and are readily available in both Anaconda or pip.

# List of functions available
Retriever takes in VCF files to create a predefined number of chimeric reference genomes and store these genomes back into a VCF file. The saved chimeric reference genomic VCF file can serve as a checkpoint, in the event that the operating time exceeds the maximum walltime of a supercomputing cluster, or used as a set of reference genomic data for other imputation softwares. 

Prior to using Retriever, please import the Retriever package into the same script for running the imputation as shown below:
>from Retriever import *


## Running entire imputation pipeline
The vcf2imput function takes an input VCF files and perform the whole pipeline of assembling chimeric references and using it to impute the missing data. A standard use of the vcf2imput function is as follows:
>vcf2imput(file_name, num_ref, window_size=100, num_threads=10, out_ref='chimeric_ref_gts.vcf', out_imputed='imputed_gts.vcf')

### Required arguments
**file_name:**  Name of the input vcf file that contains the genomic data needed for creation of chimeric reference genomes and performing imputation. The file should be in a vcf format.

**num_ref:** Number of chimeric reference genomes to create. A recommended number is approximately 25% of the total sample size.

### Optional arguments
**window_size:** The window length to consider for generating the chimeric reference genomes and the window length used for imputation. This number will be on the scale of base pairs. Theoretically, a higher number will result in higher accuracy but too high of a value can result in overfitting. Default= 100

**num_threads:** This parameter is used for the multithreading process so please choose the value according to the CPU and GPU specification of your computing device. Default= 10

**out_ref:** File name of the output chimeric reference genomic data. Default = 'chimeric_ref_gts.vcf'

**out_imputed:** File name of the imputed reference genomic data.Default = 'imputed_gts.vcf'

## Assembly of chimeric reference panel
>vcf2ref(file_name, num_ref, window_size=100, outfile='chimeric_ref_gts.vcf')

### Required arguments
**file_name:**  Name of the input vcf file that contains the genomic data needed for creation of chimeric reference genomes and performing imputation. The file should be in a vcf format.

**num_ref:** Number of chimeric reference genomes to create. A recommended number is approximately 25% of the total sample size.

### Optional arguments
**window_size:** The window length to consider for generating the chimeric reference genomes and the window length used for imputation. This number will be on the scale of base pairs. Theoretically, a higher number will result in higher accuracy but too high of a value can result in overfitting. Default= 100

**out_ref:** File name of the output chimeric reference genomic data. Default = 'chimeric_ref_gts.vcf'

## Imputation from chimeric reference panel
>chi2imput(file_name, ref_file, window_size=100, num_threads=48, out_imputed='imputed_gts.vcf')

### Required arguments
**file_name:**  Name of the input vcf file that contains the genomic data needed for creation of chimeric reference genomes and performing imputation. The file should be in a vcf format.

### Optional arguments
**window_size:** The window length to consider for generating the chimeric reference genomes and the window length used for imputation. This number will be on the scale of base pairs. Theoretically, a higher number will result in higher accuracy but too high of a value can result in overfitting. Default= 100

**num_threads:** This parameter is used for the multithreading process so please choose the value according to the CPU and GPU specification of your computing device. Default= 10

**out_imputed:** File name of the imputed reference genomic data.Default = 'imputed_gts.vcf'
