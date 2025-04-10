# Retriever
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
Retriever is designed for use in VCF version 4.0 onwards and all samples are aligned to a single reference genome.

## Filtering
1. If the dataset contains samples from multiple populations, the VCF files should be split so that each file contains only a single population, as this reduces population structure bias and improves imputation accuracy.
2. Similarly, each chromosome should be contained in its own VCF file to allow for parallel processing during chimeric reference panel assembly, enabling each chromosome to be processed independently using a single CPU.
3. The input VCF should be filtered to retain only biallelic variants, and we recommend excluding variants with more than 20% missing genotypes. Additionally, while imputation can be performed with a approximately 200 samples, larger sample sizes generally lead to higher-quality reference panels and more accurate imputation results.

An example of the filtering criteria using bcftools is shown below. Replace name-of-the-chromosome-in-the-vcf-file, name_of_samples.vcf.gz, and filtered_file.vcf with:
- the chromosome name present in your VCF file (e.g., chr1),
- the name or path of your input sample VCF file, and
- the desired name for your output file, respectively.
```
bcftools view -m2 -M2 -v snps --regions 'name-of-the-chromosome-in-the-vcf-file' name_of_samples.vcf.gz -i 'F_MISSING<0.2' -o filtered_file.vcf
```


# Running Retriever
Retriever takes in VCF files to create a predefined number of chimeric reference genomes and store these genomes back into a VCF file. The saved chimeric reference genomic VCF file can be used as the reference panel in imputation software such as Beagle.

Prior to using Retriever, please import the Retriever package into the same script for running the imputation as shown below:
```
from Retriever import *
```

## Assembly of chimeric reference panel
```
vcf2ref(file_name, num_ref, window_size=1000, outfile='chimeric_ref_gts.vcf')
```
### Required arguments
**file_name:**  Name of the input vcf file that contains the genomic data needed for creation of chimeric reference genomes and performing imputation. The file should be in a vcf format.

**num_ref:** Number of chimeric reference genomes to create. A recommended number is approximately 25% of the total sample size or 50 individuals.

### Optional arguments
**window_size:** The window length to consider for assembly of the chimeric reference genomes. This number will be on the scale of base pairs. Theoretically, a higher number will result in larger complete partitions but a larger number may also increase the computational time. Default= 1000

**out_ref:** File name of the output chimeric reference genomic data. Default = 'chimeric_ref_gts.vcf'

