# Retriever
The program uses a VCF input file to create a predefined number of chimeric reference genomes. The reference genomes can be used for imputation either from softwares available on the market or a built-in algorithm to perform imputation through machine learning techniques.

# Packages used:
Please install the packages required using conda by typing the line below into the powershell or terminal (for Mac computers)
conda install numpy pandas cyvcf2 scikit-learn

# List of functions available
Prior to using Retriever, please import the Retriever package into the same script for running the imputation as shown below:
>import Retriever
>
Retriever takes in VCF files to create a predefined number of chimeric reference genomes and store these genomes back into a VCF file. The saved chimeric reference genomic VCF file can serve as a checkpoint, in the event that the operating time exceeds the maximum walltime of a supercomputing cluster, or used as a set of reference genomic data for other imputation softwares. 
## Entire imputation pipeline
The vcf2imput function takes an input VCF files and perform the whole pipeline of generating chimeric reference genomic data and using it to impute the missing data. A standard use of the vcf2imput function is as follows:
>vcf2imput(file_name, num_ref, window_size=100, num_threads=10, out_ref='chimeric_ref_gts.vcf', out_imputed='imputed_gts.vcf')
### Required arguments
#### file_name: Name of the input vcf file that contains the genomic data needed for creation of chimeric reference genomes and performing imputation. The file should be in a vcf format.
#### num_ref: Number of chimeric reference genomes to create. A recommended number is approximately 25% of the total sample size.
### Optional arguments
#### window_size: The window length to consider for generating the chimeric reference genomes and the window length used for imputation. This number will be on the scale of base pairs. Theoretically, a higher number will result in higher accuracy but too high of a value can result in overfitting. Default= 100
#### num_threads: This parameter is used for the multithreading process so please choose the value according to the CPU and GPU specification of your computing device. Default= 10
#### out_ref: File name of the output chimeric reference genomic data. Default = 'chimeric_ref_gts.vcf'
#### out_imputed: File name of the imputed reference genomic data.Default = 'imputed_gts.vcf'
