# Retriever
Retriever is a framework that bypasses the need for external reference panels. Retriever constructs a chimeric reference panel directly from the target samples using a sliding-window approach to identify and retrieve genomic partitions with complete data. By exploiting the complementary distribution of missing data across samples, Retriever assembles a panel that preserves local patterns of linkage disequilibrium and captures novel variants. Integrated with Beagle, Retriever achieves high genotype imputation accuracy.
The Retriever program uses a VCF input file to create a predefined number of chimeric reference genomes. These reference genomes can be used for imputation with either conventional imputation software or a built-in algorithm that employs machine learning. Although Retriever has been designed to be user-friendly, a basic understanding of coding with Python is assumed, requring users to execute the imputation pipeline through simple command lines with minimal additional programming.

This software is developed under the Ortiz-Barrientos lab.

# Contents
1. [Installation](#Installation)
2. [Input file checklist](#Assumptions-of-the-input-file)
3. [Functions available](#Running-Retriever)



# Installation:

Retriever makes use of external packages for some processes so it is recommended to install these packages in a virtual environment:
```
conda env create -f Retriever.yml
```

Alternatively, if the environment file above does not work, please create the virtual environment manually as follows:

```
conda create -n Retriever
conda activate Retriever
```
Next, please install the packages required by typing the following line into the powershell or terminal:
```
conda install numpy pandas cyvcf2 scikit-learn
```
These packages are prerequsites for running Retriever and are readily available in both Anaconda or pip.

# Assumptions of the input file
## Format of VCF file
Retriever is designed for use in VCF version 4.0 onwards and all samples are aligned to a single reference genome. While imputation can be performed with approximately 200 samples, larger sample sizes generally lead to higher-quality reference panels and more accurate imputation results. Both compressed (file_name.vcf.gz) and uncompressed (file_name.vcf) VCF files are compatible with Retriever.

## Filtering
1. If the dataset contains samples from multiple populations, the VCF files should be split so that each file contains only a single population, as this reduces population structure bias and improves imputation accuracy.
2. Similarly, each chromosome should be contained in its own VCF file to allow for parallel processing during chimeric reference panel assembly, enabling each chromosome to be processed independently using a single CPU.
3. The input VCF should be filtered to retain only biallelic variants, and we recommend excluding variants with more than 20% missing genotypes.
5. Quality check for VCF file is assumed to be performed prior to using Retriever.
6. Diploid individuals to be used.
7. Duplicated positions are to be removed to avoid overlapping positions when assembling the reference panel.
8. The positions are sorted in ascending order.

Filtering of the SNP file can be performed using tools such as VCFtools and bcftools. Here is an example of using bcftools to extract one chromosome for analysis, removal of multiallelic sites and retaining sites that contain less than 80% missing genotypes. Replace name-of-the-chromosome-in-the-vcf-file, name_of_samples.vcf.gz, and filtered_file.vcf with:
- the chromosome name present in your VCF file (e.g., chr1),
- the name or path of your input sample VCF file, and
- the desired name for your output file, respectively.
```
bcftools view -m2 -M2 -v snps --regions 'name-of-the-chromosome-in-the-vcf-file' name_of_samples.vcf.gz -i 'F_MISSING<0.2'  -o filtered_file.vcf
```
After imputation, the vcf files can be merged back into a single file, if required for downstream analysis. An example of the merger is shown below using bcftools:
```
bcftools merge input_file1.vcf.gz input_file2.vcf.gz -o merged_files.vcf.gz
```

# Running Retriever
Retriever takes in VCF files to create a predefined number of chimeric reference genomes and store these genomes back into a VCF file. The saved chimeric reference genomic VCF file can be used as the reference panel in imputation software such as Beagle.

Prior to using Retriever, please import the Retriever package into the same script for running the imputation as shown below:
```
from Retriever import *
```

## Assembly of chimeric reference panel
```
vcf2ref('file_name.vcf', num_ref, window_size=1000, outfile='chimeric_ref_gts.vcf')
```
### Required arguments
**file_name:**  Name of the input vcf file that contains the genomic data needed for creation of chimeric reference genomes and performing imputation. The file should be in a vcf format.

**num_ref:** Number of chimeric reference genomes to create. A recommended number is approximately 25% of the total sample size or 50 individuals.

### Optional arguments
**window_size:** The window length to consider for assembly of the chimeric reference genomes. This number will be on the scale of base pairs. Theoretically, a higher number will result in larger complete partitions but a larger number may also increase the computational time. Default= 1000

**out_ref:** File name of the output chimeric reference genomic data. Default = 'chimeric_ref_gts.vcf'

## Parallelisation across chromosomes
The multiprocessing tool in python allows several chromosomes to be run with Retriever concurrently. An example is provided below for running multiple chromosomes in parallel:
```
chromosome_numbers= ['chr1', 'chr2', 'chr3', 'chr4', 'chr5'] #input the chromosome names in here and enclose all names in apostrophes
def par_Ret(no):
  vcf2ref('file_name_' + str(no) +'.vcf', num_ref, window_size=1000, outfile='chimeric_ref_gts' + str(no) +'.vcf')
from multiprocessing import Pool
with Pool(48) as p:
    p.map(par_Ret, chromosome_numbers)

```
**file_name** when splitting the files based on the chromosomes, naming the files to contain the chromosome number will result in easier 
**chromosome_numbers** can be in either numbers or strings format

Please visit https://docs.python.org/3/library/multiprocessing.html for more information on multiprocessing.

# Special considerations for running Retriever
1. If the number of individuals within any window is less than the number required to assemble the chimeric panel, the program will terminate with an error. To avoid this, ensure that each position includes at least as many genotypes as the chimeric panel size.
2. Buckets are temporarily stored in local memory. To prevent memory-related issues, allocate at least three times the expected data size in available memory when running Retriever.
4. Ensure only 1 chromosome is in each VCF file input to Retriever.
5. Only 1 CPU core is required for computing so more CPU access will not expediate the program,
6. However, availability to more CPU allows parallelisation of files from multiple chromosomes.
